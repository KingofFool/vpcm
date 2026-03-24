!****************************************************************
!
!     Photochemical routine 
!
!     Authors: Franck Lefevre, Francisco Gonzalez-Galindo
!     -------
!
!     Version: 14/11/2020
!
!     ASIS scheme : for details on the method see
!     Cariolle et al., Geosci. Model Dev., 10, 1467-1485, 2017.
!
!     -------
!
!     2022/09/15: adding ion chemistry by Antoine Martinez
!
!*****************************************************************

subroutine photochemistry_venus(nz, n_lon, zlocal, ptimestep,             &
                                ok_jonline, ok_ionchem, tuneupperatm,     &
                                nb_reaction_3_max, nb_reaction_4_max,     &
                                nb_phot_max, nphotion, ig,                &
                                p, t, t_elect, tr, vmr_dens_euv, mumean,  &
                                sza_input, lon, lat, nesp, nespeuv, iter, &
                                prod_tr, loss_tr, em_no, em_o2)  

use chemparam_mod
use photolysis_mod
use param_v4_h, only: jion
use iono_h, only: phdisrate
      
implicit none

!===================================================================
!     input:
!===================================================================

integer, intent(in) :: nz          ! number of atmospheric layers
integer, intent(in) :: nesp        ! number of tracers in traceur.def
integer, intent(in) :: nespeuv     ! number of tracers for jthermcal_e107.F

integer, intent(in) :: nb_reaction_3_max   ! number of quadratic reactions
integer, intent(in) :: nb_reaction_4_max   ! number of bimolecular reactions
integer, intent(in) :: nb_phot_max         ! total number of photolysis+photoionizations+quenching reactions
integer, intent(in) :: nphotion            ! number of photoionizations

logical, intent(in) :: ok_ionchem  ! switch for ion chemistry
logical, intent(in) :: ok_jonline  ! switch for on-line calculation of photolysis rates
logical, intent(in) :: tuneupperatm! upper atmosphere chemistry tuning

real, dimension(nz) :: p                 ! pressure (hpa)
real, dimension(nz) :: t                 ! temperature (k)
real, dimension(nz) :: t_elect           ! electronic temperature (k)
real, dimension(nz) :: zlocal            ! altitude (km)
real, dimension(nz) :: mumean            ! mean molecular mass (g/mol)
real, dimension(nz,nespeuv) :: vmr_dens_euv ! tracer mixing ratio for jthermcal_e107 routine
real :: ptimestep                        ! physics timestep (s)
real :: sza_input                        ! solar zenith angle (degrees)
real :: lon, lat                         ! longitude and latitude (degrees)

integer :: ig                            ! grid point index

integer :: n_lon                         ! for 1D test

!===================================================================
!     input/output:
!===================================================================

real, dimension(nz,nesp) :: tr      ! tracer mixing ratio
real, dimension(nz,nesp) :: prod_tr ! production (cm-3.s-1)
real, dimension(nz,nesp) :: loss_tr ! loss       (cm-3.s-1)
real, dimension (nz)     :: em_no   ! volume emission rate of no
real, dimension (nz)     :: em_o2   ! volume emission rate of o2(deltag)

!===================================================================
!     output:
!===================================================================

integer :: iter(nz)               ! iteration counter

!===================================================================
!     local: 
!===================================================================

! ok_ jonline: see physiq.def
! true : on-line calculation of photodissociation rates ! false : lookup table

logical, save :: firstcall = .true.

real, dimension(nz)  :: conc      ! total number density (molecule.cm-3)
real, dimension(nz)  :: surfice1d, surfdust1d

integer :: ind_norec
integer :: ind_orec

! photolysis lookup table (case jonline = .false.)
! if prior to jvenus.20211025, set nztable = 201 below

integer, parameter :: nj = 23, nztable = 281, nsza = 27, nso2 = 13
real, dimension(nso2,nsza,nztable,nj), save :: jphot ! nj must be equal to nphot
real, dimension(nztable), save :: table_colair
real, dimension(nso2,nztable), save :: table_colso2
real, dimension(nsza), save :: table_sza
real :: dist_sol

! number densities

real (kind = 8), dimension(nesp)        :: cold ! number densities at previous timestep (molecule.cm-3) 
real (kind = 8), dimension(nz,nesp)     :: c    ! number densities at current timestep (molecule.cm-3) 
real (kind = 8), dimension(nz,nespeuv)  :: c_euv! number densities for jthermcal_e107 at current timestep (molecule.cm-3) 
real (kind = 8), dimension(nesp)        :: cnew ! number densities at next timestep (molecule.cm-3) 

! timesteps

real :: ctimestep           ! standard timestep for the chemistry (s) 
real :: dt_guess            ! first-guess timestep (s) 
real :: dt_corrected        ! corrected timestep (s) 
real :: time                ! internal time (between 0 and ptimestep, in s)
integer :: phychemrat

!Tracer indexes for photionization coeffs

integer,parameter :: induv_co2 = 1
integer,parameter :: induv_o2  = 2
integer,parameter :: induv_o   = 3
integer,parameter :: induv_h2o = 4
integer,parameter :: induv_h2  = 5
integer,parameter :: induv_h2o2= 6
integer,parameter :: induv_o3  = 7
integer,parameter :: induv_n2  = 8
integer,parameter :: induv_n   = 9
integer,parameter :: induv_no  = 10
integer,parameter :: induv_co  = 11
integer,parameter :: induv_h   = 12
integer,parameter :: induv_no2 = 13

! reaction rates

real (kind = 8), dimension(nz,      nb_phot_max) :: v_phot
real (kind = 8), dimension(nz,nb_reaction_3_max) :: v_3
real (kind = 8), dimension(nz,nb_reaction_4_max) :: v_4

logical,parameter :: hetero_ice  = .false.
logical,parameter :: hetero_dust = .false.

! matrix

real (kind = 8), dimension(nesp,nesp) :: mat, mat1
integer, dimension(nesp)              :: indx
integer                               :: code

! production and loss terms (for first-guess solution only)

real (kind = 8), dimension(nesp) :: prod, loss, lossconc

! indexes

integer :: i, iesp, iz

if (firstcall) then
!===================================================================
!     initialisation of  photolysis
!===================================================================

   ! ok jonline
   ! true : on-line calculation of photodissociation rates ! false : lookup table
   if (ok_jonline) then
      print*, 'photochemistry: Read UV absorption cross-sections:'
      call init_photolysis
   else
      print*, 'photochemistry: Read photolysis lookup table:'
      call init_chimie(nphot, nztable, nsza, nso2, jphot, table_colair, table_colso2, table_sza)
   end if

!===================================================================
!     initialisation of the reaction indexes
!===================================================================

   call indice(ok_ionchem, nb_phot_max, nb_reaction_3_max, nb_reaction_4_max)

   firstcall = .false.
end if

! cloud and dust surfaces set to zero for the moment

surfice1d(:)  = 0.
surfdust1d(:) = 0.
      
!===================================================================
!   number densities (molecule.cm-3)
!===================================================================

do iz = 1,nz
   conc(iz)    = p(iz)/(1.38E-19*t(iz))
   c(iz,:)     = tr(iz,:)*conc(iz)
end do
      
!===================================================================
!    photodissociations         
!===================================================================

! dist_sol : sun-venus distance (au)

dist_sol = 0.72333

if (ok_jonline) then
   if (sza_input <= 95.) then ! day at 300 km
      call photolysis_online(nz, nb_phot_max,                                       &
                             zlocal, p, t, mumean,                                  &
                             i_co2, i_co, i_o, i_o1d, i_o2, i_o3, i_h2,             &
                             i_oh, i_ho2, i_h2o2, i_h2o, i_h, i_hcl,                &
                             i_cl, i_clo, i_cl2, i_hocl, i_so2, i_so, i_so3, i_s2,  &
                             i_osso_cis, i_osso_trans, i_s2o2_cyc, i_clso2,         &
                             i_cl2so2, i_ocs, i_cocl2, i_h2so4,                     &
                             i_no2, i_no, i_n2, i_n2d,                              &
                             nesp, tr, sza_input, dist_sol, v_phot)                 

      !Calculation of photoionization rates, if needed
      if (ok_ionchem) then
         do iz = 1,nz
           c_euv(iz,:) = vmr_dens_euv(iz,:)*conc(iz)
         end do
      !! FAIRE ULTRA ATTENTION LE c_euv et vmr_dens_euv NE SONT PAS DANS LE MEME REGIME 
      !! DE i_ESPECE QUE C, CNEW, COLD et TR
         call jthermcalc_e107(ig,nz,2,c_euv,nespeuv,t,zlocal,sza_input)
         do iz=1,nz
            call phdisrate(ig,nz,2,sza_input,iz)
         end do
         !CO2 photoionization
         v_phot(:,nphot+ 1) = jion(induv_co2,:,1)
         v_phot(:,nphot+ 2) = jion(induv_co2,:,2)
         v_phot(:,nphot+ 3) = jion(induv_co2,:,2)
         v_phot(:,nphot+ 4) = jion(induv_co2,:,3)
         v_phot(:,nphot+ 5) = jion(induv_co2,:,3)
         v_phot(:,nphot+ 6) = jion(induv_co2,:,4)
         v_phot(:,nphot+ 7) = jion(induv_co2,:,4)
         !O2 photoionization
         v_phot(:,nphot+ 8) = jion(induv_o2,:,1)
         !O photoionization
         v_phot(:,nphot+ 9) = jion(induv_o,:,1)
         !NO photoionization
         v_phot(:,nphot+10) = jion(induv_no,:,1)
         !CO photoionization
         v_phot(:,nphot+11) = jion(induv_co,:,1)
         v_phot(:,nphot+12) = jion(induv_co,:,2)
         v_phot(:,nphot+13) = jion(induv_co,:,2)
         !N2 photoionization
         v_phot(:,nphot+14) = jion(induv_n2,:,1)
         v_phot(:,nphot+15) = jion(induv_n2,:,2)
         v_phot(:,nphot+16) = jion(induv_n2,:,2)
         !N photoionization
         v_phot(:,nphot+17) = jion(induv_n,:,1)
         !H photoionization
         v_phot(:,nphot+18) = jion(induv_h,:,1)
         !D photoionization
         v_phot(:,nphot+19) = jion(induv_h,:,1)
      end if
   else ! night
      v_phot(:,:) = 0.
   end if
else 
   call phot(nj, nztable, nsza, nso2, sza_input, dist_sol, mumean, tr(:,i_co2), tr(:,i_so2),         &
             jphot, table_colair, table_colso2, table_sza, nz, nb_phot_max, t, p, v_phot)
end if

!===================================================================
!     reaction rates                                     
!===================================================================
                   
call krates(hetero_ice,hetero_dust,ok_ionchem, nphotion, nz, nesp, c, conc, t, t_elect, p, nb_phot_max, nb_reaction_3_max, &
            nb_reaction_4_max, tuneupperatm, v_3, v_4, v_phot, sza_input, ind_norec, ind_orec)

!===================================================================
!     ctimestep : standard chemical timestep (s), defined as 
!                 the fraction phychemrat of the physical timestep                           
!===================================================================

phychemrat = 1

ctimestep  = ptimestep/real(phychemrat)

!===================================================================
!     loop over levels         
!===================================================================

do iz = 1,nz

!  initializations

   time = 0.
   iter(iz) = 0
   dt_guess = ctimestep
   cold(:) = c(iz,:)     

!  internal loop for the chemistry
          
   do while (time < ptimestep)
      
   iter(iz) = iter(iz) + 1

!  first-guess: fill matrix

   call fill_matrix(iz, mat1, prod, loss, lossconc, c, nesp, nz,          &
                    nb_reaction_3_max, nb_reaction_4_max, nb_phot_max,    &
                    v_phot, v_3, v_4)

!  adaptative evaluation of the sub time step

   call define_dt(nesp, dt_corrected, dt_guess, ctimestep, cold(:), c(iz,:), &
                  mat1, prod, loss, conc(iz), lon, lat)

   if (time + dt_corrected > ptimestep) then
      dt_corrected = ptimestep - time
   end if

!  form the matrix identity + mat*dt_corrected

   mat(:,:) = mat1(:,:)*dt_corrected
   do iesp = 1,nesp
      mat(iesp,iesp) = 1. + mat(iesp,iesp)
   end do

!  solve the linear system  M*Cn+1 = Cn (RHS in cnew, then replaced by solution)

   cnew(:) = c(iz,:)

#ifdef LAPACK
   call dgesv(nesp,1,mat,nesp,indx,cnew,nesp,code)
#else
   write(*,*) "photochemistry error, missing LAPACK routine dgesv"
   stop
#endif

!  eliminate small values

   where (cnew(:)/conc(iz) < 1.e-30)
      cnew(:) = 0.
   end where

!  update concentrations

   cold(:) = c(iz,:)
   c(iz,:) = cnew(:)
   cnew(:) = 0.

!  force charge neutrality (mod fgg, july 2019)

   if (ok_ionchem) then
      if(c(iz,i_elec).ne.c(iz,i_co2plus)+c(iz,i_oplus)+c(iz,i_o2plus)+  &
           c(iz,i_noplus)+c(iz,i_coplus)+c(iz,i_cplus)+c(iz,i_n2plus)+  &
           c(iz,i_nplus)+c(iz,i_hplus)+c(iz,i_hco2plus)+                &
           c(iz,i_hcoplus)+c(iz,i_h2oplus)+c(iz,i_h3oplus)+             &
           c(iz,i_ohplus)+c(iz,i_dco2plus)+c(iz,i_dcoplus)+             &
           c(iz,i_hdoplus)+c(iz,i_dplus)+c(iz,i_odplus)+                &
           c(iz,i_h2doplus) ) then
         c(iz,i_elec) = c(iz,i_co2plus)+c(iz,i_oplus)+c(iz,i_o2plus)+   &
              c(iz,i_noplus)+c(iz,i_coplus)+c(iz,i_cplus)+              &
              c(iz,i_n2plus)+c(iz,i_nplus)+c(iz,i_hplus)+               &
              c(iz,i_hco2plus)+c(iz,i_hcoplus)+c(iz,i_h2oplus)+         &
              c(iz,i_h3oplus)+c(iz,i_ohplus)+c(iz,i_dco2plus)+          & 
              c(iz,i_dcoplus)+c(iz,i_hdoplus)+c(iz,i_dplus)+            & 
              c(iz,i_odplus)+c(iz,i_h2doplus)
         !      write(*,*)'photochemistry/359'
         !      write(*,*)'Forcing charge neutrality at ilev,',ilev,' ig=',ig
      end if
   end if

!  increment internal time

   time = time + dt_corrected
   dt_guess = dt_corrected     ! first-guess timestep for next iteration

   end do ! while (time < ptimestep)

!  save mixing ratios

   tr(iz,:)  = max(c(iz,:)/conc(iz), 1.e-30)

!  save production and loss terms (for diagnostic only)
   
   prod_tr(iz,:) = prod(:)
   loss_tr(iz,:) = lossconc(:)
                       
end do  ! end of loop over vertical levels
 
! no and o2(delta) emissions

em_no(:) = c(:,i_o)*c(:,i_n)*v_4(:,ind_norec)   !2.8e-17*(300./temp(:)))**0.5
em_o2(:) = c(:,i_o2dg)/4470.    ! 4470 s : lafferty et al., 1998

end subroutine photochemistry_venus

!======================================================================

 subroutine init_chimie(nj, nztable, nsza, nso2, jphot, table_colair, &
                        table_colso2, table_sza)

!======================================================================

implicit none

! photolysis lookup table

integer, INTENT(IN) :: nj, nztable, nsza, nso2
real, INTENT(OUT), dimension(nso2,nsza,nztable,nj) :: jphot
real, INTENT(OUT), dimension(nztable) :: table_colair
real, INTENT(OUT), dimension(nso2,nztable) :: table_colso2
real, INTENT(OUT), dimension(nsza) :: table_sza

integer           :: iz, isza, iozo, iso2, ij
character(len=44) :: jvenus

! lecture de la table des j

jphot(:,:,:,:) = 0.

jvenus = 'jvenus.dat'
open(30, form = 'formatted', status = 'old', file = jvenus)
print*,'lecture de jvenus = ', jvenus

do iso2 = 1,nso2
   do isza = 1,nsza
      do iz = nztable,1,-1
         read(30,*) table_colair(iz), table_colso2(iso2,iz), table_sza(isza)
         read(30,'(7e11.4)') (jphot(iso2,isza,iz,ij), ij = 1,nj)
         do ij = 1,nj
            if (jphot(iso2,isza,iz,ij) == 1.E-30) then 
               jphot(iso2,isza,iz,ij) = 0.
            end if
         end do
      end do
   end do
end do

close(30)
print*,'lecture de la table des j ok.'

end subroutine init_chimie

!======================================================================

 subroutine indice(ok_ionchem, nb_phot_max, nb_reaction_3_max, nb_reaction_4_max)

!================================================================
! set the "indice" arrays used to fill the jacobian matrix      !
!----------------------------------------------------------------
! reaction               type                array              !
!----------------------------------------------------------------
! A + hv   --> B + C     photolysis          indice_phot        ! 
! A + B    --> C + D     bimolecular         indice_4           !
! A + A    --> B + C     quadratic           indice_3           !
! A + C    --> B + C     quenching           indice_phot        !
! A + ice  --> B + C     heterogeneous       indice_phot        !
!================================================================

use types_asis
use chemparam_mod
 
implicit none

! input

integer, intent(in) :: nb_reaction_3_max   ! number of quadratic reactions
integer, intent(in) :: nb_reaction_4_max   ! number of bimolecular reactions
integer, intent(in) :: nb_phot_max         ! number of processes treated numerically as photodissociations
logical, intent(in) :: ok_ionchem          ! True: Ion reaction

! local

integer :: nb_phot, nb_reaction_3, nb_reaction_4
integer :: i_dummy

allocate (indice_phot(nb_phot_max))
allocate (indice_3(nb_reaction_3_max))
allocate (indice_4(nb_reaction_4_max))

i_dummy = 1

nb_phot       = 0
nb_reaction_3 = 0
nb_reaction_4 = 0

!===========================================================
!      O2 + hv -> O + O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_o2, 2.0, i_o, 0.0, i_dummy)

!===========================================================
!      O2 + hv -> O + O(1D)
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_o2, 1.0, i_o, 1.0, i_o1d)

!===========================================================
!      CO2 + hv -> CO + O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_co2, 1.0, i_co, 1.0, i_o)

!===========================================================
!      CO2 + hv -> CO + O(1D)
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_co2, 1.0, i_co, 1.0, i_o1d)

!===========================================================
!      O3 + hv -> O2(Dg) + O(1D)
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_o3, 1.0, i_o2dg, 1.0, i_o1d)

!===========================================================
!      O3 + hv -> O2 + O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_o3, 1.0, i_o2, 1.0, i_o)

!===========================================================
!      H2 + hv -> H + H
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_h2, 2.0, i_h, 0.0, i_dummy)

!===========================================================
!      H2O + hv -> H + OH
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_h2o, 1.0, i_h, 1.0, i_oh)

!===========================================================
!      HO2 + hv -> OH + O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_ho2, 1.0, i_oh, 1.0, i_o)

!===========================================================
!      H2O2 + hv -> OH + OH
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_h2o2, 2.0, i_oh, 0.0, i_dummy)

!===========================================================
!      HCl + hv -> H + Cl
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_hcl, 1.0, i_h, 1.0, i_cl)

!===========================================================
!      Cl2 + hv -> Cl + Cl
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_cl2, 2.0, i_cl, 0.0, i_dummy)

!===========================================================
!      HOCl + hv -> OH + Cl
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_hocl, 1.0, i_oh, 1.0, i_cl)

!===========================================================
!      ClO + hv -> Cl + O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_clo, 1.0, i_cl, 1.0, i_o)

!===========================================================
!      SO2 + hv -> SO + O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_so2, 1.0, i_so, 1.0, i_o)

!===========================================================
!      SO + hv -> S + O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_so, 1.0, i_s, 1.0, i_o)

!===========================================================
!      SO3 + hv -> SO2 + O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_so3, 1.0, i_so2, 1.0, i_o)

!===========================================================
!       S2 + hv -> S + S
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_s2, 2.0, i_s, 0.0, i_dummy)

!===========================================================
!       OSSO_cis + hv -> SO + SO
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_osso_cis, 2.0, i_so, 0.0, i_dummy)

!===========================================================
!       OSSO_trans + hv -> SO + SO
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_osso_trans, 2.0, i_so, 0.0, i_dummy)

!===========================================================
!       S2O2_cyc + hv -> SO + SO
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_s2o2_cyc, 2.0, i_so, 0.0, i_dummy)

!===========================================================
!      ClSO2 + hv -> Cl + SO2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_clso2, 1.0, i_cl, 1.0, i_so2)

!===========================================================
!      Cl2SO2 + hv -> Cl + ClSO2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_cl2so2, 1.0, i_cl, 1.0, i_clso2)

!===========================================================
!      OCS + hv -> CO + S
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_ocs, 1.0, i_co, 1.0, i_s)

!===========================================================
!      COCl2 + hv -> Cl + Cl + CO
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_cocl2, 2.0, i_cl, 1.0, i_co)

!===========================================================
!      H2SO4 + hv -> SO3 + H2O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_h2so4, 1.0, i_so3, 1.0, i_h2o)

!===========================================================
!      NO2 + hv -> NO + O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_no2, 1.0, i_no, 1.0, i_o)

!===========================================================
!      NO + hv -> N + O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_no, 1.0, i_n, 1.0, i_o)

!===========================================================
!      N2 + hv -> N(2D) + N
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_n2, 1.0, i_n, 1.0, i_n2d)
 
!===========================================================
!      HDO + hv -> OD + H
!===========================================================

   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(1.0, i_hdo, 1.0, i_h, 1.0, i_od) 

!===========================================================
!      HDO + hv -> D + OH
!===========================================================

   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(1.0, i_hdo, 1.0, i_d, 1.0, i_oh)

!===========================================================
!      HDSO4 + hv -> HDO + SO3
!===========================================================

   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(1.0, i_hdso4, 1.0, i_hdo, 1.0, i_so3)

!===========================================================
!      HD + hv -> H + D
!===========================================================

   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(1.0, i_hd, 1.0, i_h, 1.0, i_d)

!===========================================================
!      DO2 + hv -> OD + O
!===========================================================

   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(1.0, i_do2, 1.0, i_od, 1.0, i_o)

!===========================================================
!      HDO2 + hv -> OH + OD
!===========================================================

   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(1.0, i_hdo2, 1.0, i_oh, 1.0, i_od)

!===========================================================
!      DCl + hv -> D + Cl
!===========================================================

   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(1.0, i_dcl, 1.0, i_d, 1.0, i_cl) 

!===========================================================
!      DOCl + hv -> OD + Cl
!===========================================================

   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(1.0, i_docl, 1.0, i_od, 1.0, i_cl)

!Only if ion chemistry included
if (ok_ionchem) then

!===========================================================
!      CO2 + hv -> CO2+ + e-
!===========================================================

   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(1.0, i_co2, 1.0, i_co2plus, 1.0, i_elec)

!===========================================================
!      CO2 + hv -> O+ + CO + e-
!===========================================================
!We divide this reaction in two

!0.5 CO2 + hv -> CO
   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(0.5, i_co2, 1.0, i_co, 0.0, i_dummy)

!0.5 CO2 + hv -> O+ + e-
   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(0.5, i_co2, 1.0, i_oplus, 1.0, i_elec)

!===========================================================
!      CO2 + hv -> CO+ + O + e-
!===========================================================
!We divide this reaction in two

!0.5 CO2 + hv -> O
   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(0.5, i_co2, 1.0, i_o, 0.0, i_dummy)

!0.5 CO2 + hv -> CO+ + e-
   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(0.5, i_co2, 1.0, i_coplus, 1.0, i_elec)

!===========================================================
!      CO2 + hv -> C+ + O2 + e-
!===========================================================
!We divide this reaction in two

!0.5 CO2 + hv -> O2
   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(0.5, i_co2, 1.0, i_o2, 0.0, i_dummy)

!0.5 CO2 + hv -> C+ + e-
   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(0.5, i_co2, 1.0, i_cplus, 1.0, i_elec)

!===========================================================
!      O2 + hv -> O2+ + e-
!===========================================================

   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(1.0, i_o2, 1.0, i_o2plus, 1.0, i_elec)

!===========================================================
!      O + hv -> O+ + e-
!===========================================================

   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(1.0, i_o, 1.0, i_oplus, 1.0, i_elec)

!===========================================================
!      NO + hv -> NO+ + e-
!===========================================================

   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(1.0, i_no, 1.0, i_noplus, 1.0, i_elec)

!===========================================================
!      CO + hv -> CO+ + e-
!===========================================================

   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(1.0, i_co, 1.0, i_coplus, 1.0, i_elec)

!===========================================================
!      CO + hv -> C+ + O + e-
!===========================================================
!We divide this reaction in two

!0.5 CO + hv -> O
   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(0.5, i_co, 1.0, i_o, 0.0, i_dummy)

!0.5 CO + hv -> C+ + e-
   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(0.5, i_co, 1.0, i_cplus, 1.0, i_elec)

!===========================================================
!      N2 + hv -> N2+ + e-
!===========================================================

   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(1.0, i_n2, 1.0, i_n2plus, 1.0, i_elec)

!===========================================================
!      N2 + hv -> N+ + N + e-
!===========================================================
!We divide this reaction in two

!0.5 N2 + hv -> N
   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(0.5, i_n2, 1.0, i_n, 0.0, i_dummy)

!0.5 N2 + hv -> N+ + e-
   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(0.5, i_n2, 1.0, i_nplus, 1.0, i_elec)

!===========================================================
!      N + hv -> N+ + e-
!===========================================================

   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(1.0, i_n, 1.0, i_nplus, 1.0, i_elec)

!===========================================================
!      H + hv -> H+ + e-
!===========================================================

   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(1.0, i_h, 1.0, i_hplus, 1.0, i_elec)

end if   !ok_ionchem
 
if (ok_ionchem) then

!===========================================================
!      D + hv -> D+ + e-
!===========================================================

   nb_phot = nb_phot + 1

   indice_phot(nb_phot) = z3spec(1.0, i_d, 1.0, i_dplus, 1.0, i_elec)

end if   !ok_ionchem
 

!===========================================================
!      a001 : O + O2 + (CO2 or M) -> O3 + (CO2 or M)
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o, 1.0, i_o2, 1.0, i_o3, 0.0, i_dummy)

!===========================================================
!      a002 : O + O + (CO2 or M) -> O2(Dg) + (CO2 or M)
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_o, 1.0, i_o2dg, 0.0, i_dummy)

!===========================================================
!      a003 : O + O3 -> O2 + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o, 1.0, i_o3, 2.0, i_o2, 0.0, i_dummy)

!===========================================================
!      b001 : O(1D) + CO2 -> O + CO2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_o1d, 1.0, i_o, 0.0, i_dummy)

!===========================================================
!      b002 : O(1D) + H2O -> OH + OH
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o1d, 1.0, i_h2o, 2.0, i_oh, 0.0, i_dummy)

!===========================================================
!      b003 : O(1D) + H2 -> OH + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o1d, 1.0, i_h2, 1.0, i_oh, 1.0, i_h)

!===========================================================
!      b004 : O(1D) + O2 -> O + O2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_o1d, 1.0, i_o, 0.0, i_dummy)

!===========================================================
!      b005 : O(1D) + O3 -> O2 + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o1d, 1.0, i_o3, 2.0, i_o2, 0.0, i_dummy)

!===========================================================
!      b006 : O(1D) + O3 -> O2 + O + O
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o1d, 1.0, i_o3, 1.0, i_o2, 2.0, i_o)
 
!===========================================================
!      db001 : O(1D) + HDO -> OD + OH
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_o1d, 1.0, i_hdo, 1.0, i_od, 1.0, i_oh)

!===========================================================
!      db002 : O(1D) + HD -> H + OD
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_o1d, 1.0, i_hd, 1.0, i_h, 1.0, i_od)

!===========================================================
!      db003 : O(1D) + HD -> D + OH
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_o1d, 1.0, i_hd, 1.0, i_d, 1.0, i_oh)

!===========================================================
!      c001 : O + HO2 -> OH + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o, 1.0, i_ho2, 1.0, i_oh, 1.0, i_o2)

!===========================================================
!      c002 : O + OH -> O2 + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o, 1.0, i_oh, 1.0, i_o2, 1.0, i_h)

!===========================================================
!      c003 : H + O3 -> OH + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_h, 1.0, i_o3, 1.0, i_oh, 1.0, i_o2)

!===========================================================
!      c004 : H + HO2 -> OH + OH
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_h, 1.0, i_ho2, 2.0, i_oh, 0.0, i_dummy)

!===========================================================
!      c005 : H + HO2 -> H2 + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_h, 1.0, i_ho2, 1.0, i_h2, 1.0, i_o2)

!===========================================================
!      c006 : H + HO2 -> H2O + O
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_h, 1.0, i_ho2, 1.0, i_h2o, 1.0, i_o)

!===========================================================
!      c007 : OH + HO2 -> H2O + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_oh, 1.0, i_ho2, 1.0, i_h2o, 1.0, i_o2)

!===========================================================
!      c008 : HO2 + HO2 -> H2O2 + O2
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_ho2, 1.0, i_h2o2, 1.0, i_o2)

!===========================================================
!      c009 : OH + H2O2 -> H2O + HO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_oh, 1.0, i_h2o2, 1.0, i_h2o, 1.0, i_ho2)

!===========================================================
!      c010 : OH + H2 -> H2O + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_oh, 1.0, i_h2, 1.0, i_h2o, 1.0, i_h)

!===========================================================
!      c011 : H + O2 + CO2 -> HO2 + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_h, 1.0, i_o2, 1.0, i_ho2, 0.0, i_dummy)

!===========================================================
!      c012 : O + H2O2 -> OH + HO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o, 1.0, i_h2o2, 1.0, i_oh, 1.0, i_ho2)

!===========================================================
!      c013 : OH + OH -> H2O + O
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_oh, 1.0, i_h2o, 1.0, i_o)

!===========================================================
!      c014 : OH + O3 -> HO2 + O2 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_oh, 1.0, i_o3, 1.0, i_ho2, 1.0, i_o2)

!===========================================================
!      c015 : HO2 + O3 -> OH + O2 + O2 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ho2, 1.0, i_o3, 1.0, i_oh, 2.0, i_o2)

!===========================================================
!      c016 : HO2 + HO2 + CO2 -> H2O2 + O2 + CO2 
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1
indice_3(nb_reaction_3) = z3spec(2.0, i_ho2, 1.0, i_h2o2, 1.0, i_o2)

!===========================================================
!      c017 : OH + OH + CO2 -> H2O2 + CO2 
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_oh, 1.0, i_h2o2, 0.0, i_dummy)

!===========================================================
!      c018 : H + H + CO2 -> H2 + CO2 
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_h, 1.0, i_h2, 0.0, i_dummy)

!===========================================================
!      d001 : NO2 + O -> NO + O2 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_no2, 1.0, i_o, 1.0, i_no, 1.0, i_o2) 

!===========================================================
!      d002 : NO + O3 -> NO2 + O2 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_no, 1.0, i_o3, 1.0, i_no2, 1.0, i_o2) 

!===========================================================
!      d003 : NO + HO2 -> NO2 + OH 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_no, 1.0, i_ho2, 1.0, i_no2, 1.0, i_oh) 

!===========================================================
!      d004 : N + NO -> N2 + O 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_n, 1.0, i_no, 1.0, i_n2, 1.0, i_o) 

!===========================================================
!      d005 : N + O2 -> NO + O 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_n, 1.0, i_o2, 1.0, i_no, 1.0, i_o) 

!===========================================================
!      d006 : NO2 + H -> NO + OH 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_no2, 1.0, i_h, 1.0, i_no, 1.0, i_oh) 

!===========================================================
!      d007 : N + O -> NO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_n, 1.0, i_o, 1.0, i_no, 0.0, i_dummy) 

!===========================================================
!      d008 : N + HO2 -> NO + OH
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_n, 1.0, i_ho2, 1.0, i_no, 1.0, i_oh) 

!===========================================================
!      d009 : N + OH -> NO + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_n, 1.0, i_oh, 1.0, i_no, 1.0, i_h) 

!===========================================================
!      d010 : N(2D) + O -> N + O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_n2d, 1.0, i_n, 0.0, i_dummy) 

!===========================================================
!      d011 : N(2D) + N2 -> N + N2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_n2d, 1.0, i_n, 0.0, i_dummy) 

!===========================================================
!      d012 : N(2D) + CO2 -> NO + CO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_n2d, 1.0, i_co2, 1.0, i_no, 1.0, i_co) 

!===========================================================
!      d013 : N + O + CO2 -> NO + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_n, 1.0, i_o, 1.0, i_no, 0.0, i_dummy) 

!===========================================================
!      d014 : N(2D) + CO -> N + CO
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_n2d, 1.0, i_n, 0.0, i_dummy)

!===========================================================
!      e001 : CO + OH -> CO2 + H 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_co, 1.0, i_oh, 1.0, i_co2, 1.0, i_h)

!===========================================================
!      e002 : CO + O + M -> CO2 + M 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_co, 1.0, i_o, 1.0, i_co2, 0.0, i_dummy)

!===========================================================
!      f001 : HCl + O(1D) -> OH + Cl
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_hcl, 1.0, i_o1d, 1.0, i_oh, 1.0, i_cl)

!===========================================================
!      f002 : HCl + O(1D) -> H + ClO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_hcl, 1.0, i_o1d, 1.0, i_h, 1.0, i_clo)

!===========================================================
!      f003 : HCl + O -> OH + Cl
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_hcl, 1.0, i_o, 1.0, i_oh, 1.0, i_cl)

!===========================================================
!      f004 : HCl + OH -> H2O + Cl
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_hcl, 1.0, i_oh, 1.0, i_h2o, 1.0, i_cl)

!===========================================================
!      f005 : ClO + O -> Cl + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clo, 1.0, i_o, 1.0, i_cl, 1.0, i_o2)

!===========================================================
!      f006 : ClO + OH -> Cl + HO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clo, 1.0, i_oh, 1.0, i_cl, 1.0, i_ho2)

!===========================================================
!      f007 : ClO + OH -> HCl + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clo, 1.0, i_oh, 1.0, i_hcl, 1.0, i_o2)

!===========================================================
!      f008 : Cl + H2 -> HCl + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl, 1.0, i_h2, 1.0, i_hcl, 1.0, i_h)

!===========================================================
!      f009 : Cl + O3 -> ClO + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl, 1.0, i_o3, 1.0, i_clo, 1.0, i_o2)

!===========================================================
!      f010 : Cl + HO2 -> ClO + OH 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl, 1.0, i_ho2, 1.0, i_clo, 1.0, i_oh)

!===========================================================
!      f011 : Cl + HO2 -> HCl + O2 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl, 1.0, i_ho2, 1.0, i_hcl, 1.0, i_o2)

!===========================================================
!      f012 : Cl + H2O2 -> HCl + HO2 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl, 1.0, i_h2o2, 1.0, i_hcl, 1.0, i_ho2)

!===========================================================
!      f013 : Cl + CO + CO2 -> ClCO + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl, 1.0, i_co, 1.0, i_clco, 0.0, i_dummy)

!===========================================================
!      f014 : ClCO + CO2 -> Cl + CO + CO2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_clco, 1.0, i_cl, 1.0, i_co)

!===========================================================
!      f015 : ClCO + O2 + CO2 -> ClCO3 + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clco, 1.0, i_o2, 1.0, i_clco3, 0.0, i_dummy)

!===========================================================
!      f016 : 0.5 ClCO3 + 0.5 Cl -> Cl
!             0.5 ClCO3 + 0.5 Cl -> ClO + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(0.5, i_clco3, 0.5, i_cl, 1.0, i_cl, 0.0, i_dummy)

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(0.5, i_clco3, 0.5, i_cl, 1.0, i_clo, 1.0, i_co2)

!===========================================================
!      f017 : 0.5 ClCO3 + 0.5 O -> Cl
!             0.5 ClCO3 + 0.5 O -> O2 + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(0.5, i_clco3, 0.5, i_o, 1.0, i_cl, 0.0, i_dummy)

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(0.5, i_clco3, 0.5, i_o, 1.0, i_o2, 1.0, i_co2)

!===========================================================
!      f018 : ClO + HO2 -> HOCl + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clo, 1.0, i_ho2, 1.0, i_hocl, 1.0, i_o2)

!===========================================================
!      f019 : OH + HOCl -> H2O + ClO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_oh, 1.0, i_hocl, 1.0, i_h2o, 1.0, i_clo)

!===========================================================
!      f020 : O + HOCl -> OH + ClO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o, 1.0, i_hocl, 1.0, i_oh, 1.0, i_clo)

!===========================================================
!      f021 : Cl + Cl + CO2 -> Cl2 + CO2 
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_cl, 1.0, i_cl2, 0.0, i_dummy)

!===========================================================
!      f022 : ClCO + O -> Cl + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clco, 1.0, i_o, 1.0, i_cl, 1.0, i_co2)

!===========================================================
!      f023 : Cl2 + O(1D) -> Cl + ClO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl2, 1.0, i_o1d, 1.0, i_cl, 1.0, i_clo)

!===========================================================
!      f024 : Cl2 + H -> HCl + Cl
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl2, 1.0, i_h, 1.0, i_hcl, 1.0, i_cl)

!===========================================================
!      f025 : Cl + ClCO -> Cl2 + CO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl, 1.0, i_clco, 1.0, i_cl2, 1.0, i_co)

!===========================================================
!      f026 : ClCO + ClCO -> COCl2 + CO
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_clco, 1.0, i_cocl2, 1.0, i_co)

!===========================================================
!      f027 : Cl + SO2 + CO2 -> ClSO2 + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl, 1.0, i_so2, 1.0, i_clso2, 0.0, i_dummy)

!===========================================================
!      f028 : ClSO2 + O -> SO3 + Cl
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clso2, 1.0, i_o, 1.0, i_so3, 1.0, i_cl)

!===========================================================
!      f029 : ClSO2 + H -> SO2 + HCl 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clso2, 1.0, i_h, 1.0, i_so2, 1.0, i_hcl)

!===========================================================
!      f030 : ClSO2 + ClSO2 -> Cl2SO2 + SO2
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_clso2, 1.0, i_cl2so2, 1.0, i_so2)

!===========================================================
!      f031 : Cl + O + CO2 -> ClO + CO2 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl, 1.0, i_o, 1.0, i_clo, 0.0, i_dummy)

!===========================================================
!      f032 : Cl2 + O -> ClO + Cl 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl2, 1.0, i_o, 1.0, i_clo, 1.0, i_cl)

!===========================================================
!      f033 : ClCO + OH -> HOCl + CO 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clco, 1.0, i_oh, 1.0, i_hocl, 1.0, i_co)

!===========================================================
!      f034 : Cl2 + OH -> Cl + HOCl 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl2, 1.0, i_oh, 1.0, i_cl, 1.0, i_hocl)

!===========================================================
!      f035 : ClCO + O -> CO + ClO 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clco, 1.0, i_o, 1.0, i_co, 1.0, i_clo)

!===========================================================
!      f036 : ClCO + Cl2 -> COCl2 + Cl
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clco, 1.0, i_cl2, 1.0, i_cocl2, 1.0, i_cl)

!===========================================================
!      f037 : HCl + H -> H2 + Cl
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_hcl, 1.0, i_h, 1.0, i_h2, 1.0, i_cl)

!===========================================================
!      f038 : ClCO + H -> HCl + CO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clco, 1.0, i_h, 1.0, i_hcl, 1.0, i_co)

!===========================================================
!      f039 : Cl + H + M -> HCl + M
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl, 1.0, i_h, 1.0, i_hcl, 0.0, i_dummy)

!===========================================================
!      f040 : ClSO2 + Cl -> Cl2SO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clso2, 1.0, i_cl, 1.0, i_cl2so2, 0.0, i_dummy)
 
!===========================================================
!      df001 : OD + OH -> HDO + O
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_od, 1.0, i_oh, 1.0, i_hdo, 1.0, i_o) 

!===========================================================
!      df002 : OD + H2 -> HDO + H
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_od, 1.0, i_h2, 1.0, i_hdo, 1.0, i_h)

!===========================================================
!      df003 : OD + HO2 -> HDO + O2
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_od, 1.0, i_ho2, 1.0, i_hdo, 1.0, i_o2)

!===========================================================
!      df004 : OD + H2O2 -> HDO + HO2
!===========================================================
   
   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_od, 1.0, i_h2o2, 1.0, i_hdo, 1.0, i_ho2)

!===========================================================
!      df005 : O + OD -> O2 + D
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_o, 1.0, i_od, 1.0, i_o2, 1.0, i_d)

!===========================================================
!      df006 : OD + H2 -> H2O + D
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_h2, 1.0, i_od, 1.0, i_h2o, 1.0, i_d)

!===========================================================
!      df007 : OD + H -> OH + D
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_h, 1.0, i_od, 1.0, i_oh, 1.0, i_d)

!===========================================================
!      df008 : CO + OD -> CO2 + D
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_co, 1.0, i_od, 1.0, i_co2, 1.0, i_d)

!===========================================================
!      df009 : O3 + D -> O2 + OD
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_o3, 1.0, i_d, 1.0, i_o2, 1.0, i_od)

!===========================================================
!      df010 : HO2 + D -> OH + OD
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_ho2, 1.0, i_d, 1.0, i_oh, 1.0, i_od)

!===========================================================
!      df011 : HO2 + D -> HDO + O
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_ho2, 1.0, i_d, 1.0, i_hdo, 1.0, i_o)

!===========================================================
!      df012 : OH + D -> H + OD
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_oh, 1.0, i_d, 1.0, i_h, 1.0, i_od)

!===========================================================
!      df013 : H + D + CO2 -> HD + CO2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_h, 1.0, i_d, 1.0, i_hd, 0.0, i_dummy)

!===========================================================
!      df014 : D + HO2 -> HD + O2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_d, 1.0, i_ho2, 1.0, i_hd, 1.0, i_o2)


!===========================================================
!      df015 : OH + HD -> HDO + H 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_oh, 1.0, i_hd, 1.0, i_hdo, 1.0, i_h)

!===========================================================
!      df016 : OH + HD -> H2O + D 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_oh, 1.0, i_hd, 1.0, i_h2o, 1.0, i_d)

!===========================================================
!      df017 : D + O2 + CO2 -> DO2 + CO2
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_d, 1.0, i_o2, 1.0, i_do2, 0.0, i_dummy)

!===========================================================
!      df018 : OD + O3 -> DO2 + O2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_od, 1.0, i_o3, 1.0, i_do2, 1.0, i_o2)

!===========================================================
!      df019 : D + HO2 -> DO2 + H 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_d, 1.0, i_ho2, 1.0, i_do2, 1.0, i_h)

!===========================================================
!      df020 : O + DO2 -> OD + O2
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_o, 1.0, i_do2, 1.0, i_od, 1.0, i_o2)

!===========================================================
!      df021 : H + DO2 -> OH + OD
!===========================================================
   
   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_h, 1.0, i_do2, 1.0, i_od, 1.0, i_oh)

!===========================================================
!      df022 : H + DO2 -> HD + O2
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_h, 1.0, i_do2, 1.0, i_hd, 1.0, i_o2)

!===========================================================
!      df023 : H + DO2 -> HDO + O
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_h, 1.0, i_do2, 1.0, i_hdo, 1.0, i_o)

!===========================================================
!      df024 : H + DO2 -> HO2 + D
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_h, 1.0, i_do2, 1.0, i_ho2, 1.0, i_d)

!===========================================================
!      df025 : OH + DO2 -> HDO + O2
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_oh, 1.0, i_do2, 1.0, i_hdo, 1.0, i_o2)

!===========================================================
!      df026 : DO2 + O3 -> OD + O2 + O2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_do2, 1.0, i_o3, 1.0, i_od, 2.0, i_o2)

!===========================================================
!      df027 : OD + OH + CO2 -> HDO2 + CO2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_od, 1.0, i_oh, 1.0, i_hdo2, 0.0, i_dummy)

!===========================================================
!      df028 : DO2 + HO2 -> HDO2 + O2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_do2, 1.0, i_ho2, 1.0, i_hdo2, 1.0, i_o2)

!===========================================================
!      df029 : O + HDO2 -> OD + HO2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_o, 1.0, i_hdo2, 1.0, i_od, 1.0, i_ho2)

!===========================================================
!      df030 : O + HDO2 -> OH + DO2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_o, 1.0, i_hdo2, 1.0, i_oh, 1.0, i_do2)

!===========================================================
!      df031 : OH + HDO2 -> HDO + HO2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_oh, 1.0, i_hdo2, 1.0, i_hdo, 1.0, i_ho2)

!===========================================================
!      df032 : OH + HDO2 -> H2O + DO2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_oh, 1.0, i_hdo2, 1.0, i_h2o, 1.0, i_do2)

!===========================================================
!      df033 : OD + H2O2 -> H2O + DO2
!===========================================================
   
   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_od, 1.0, i_h2o2, 1.0, i_h2o, 1.0, i_do2)

!===========================================================
!      df034 : NO + DO2 -> NO2 + OD 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_no, 1.0, i_do2, 1.0, i_no2, 1.0, i_od) 

!===========================================================
!      df035 : NO2 + D -> NO + OD 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_no2, 1.0, i_d, 1.0, i_no, 1.0, i_od) 

!===========================================================
!      df036 : N + DO2 -> NO + OD
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_n, 1.0, i_do2, 1.0, i_no, 1.0, i_od) 

!===========================================================
!      df037 : N + OD -> NO + D
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_n, 1.0, i_od, 1.0, i_no, 1.0, i_d) 

!===========================================================
!      df038 : DO2 + HO2 + CO2 -> HDO2 + O2 + CO2 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_do2, 1.0, i_ho2, 1.0, i_hdo2, 1.0, i_o2)

!===========================================================
!      df039 : DCl + O(1D) -> OD + Cl
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_dcl, 1.0, i_o1d, 1.0, i_od, 1.0, i_cl)

!===========================================================
!      df040 : DCl + O(1D) -> D + ClO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_dcl, 1.0, i_o1d, 1.0, i_d, 1.0, i_clo)

!===========================================================
!      df041 : DCl + O -> OD + Cl
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_dcl, 1.0, i_o, 1.0, i_od, 1.0, i_cl)

!===========================================================
!      df042 : HCl + OD -> HDO + Cl
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_hcl, 1.0, i_od, 1.0, i_hdo, 1.0, i_cl)

!===========================================================
!      df043 : DCl + OH -> HDO + Cl
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_dcl, 1.0, i_oh, 1.0, i_hdo, 1.0, i_cl)

!===========================================================
!      df044 : ClO + OD -> Cl + DO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clo, 1.0, i_od, 1.0, i_cl, 1.0, i_do2)

!===========================================================
!      df045 : ClO + OD -> DCl + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clo, 1.0, i_od, 1.0, i_dcl, 1.0, i_o2)

!===========================================================
!      df046 : Cl + HD -> DCl + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl, 1.0, i_hd, 1.0, i_dcl, 1.0, i_h)

!===========================================================
!      df047 : Cl + HD -> HCl + D
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl, 1.0, i_hd, 1.0, i_hcl, 1.0, i_d)

!===========================================================
!      df048 : Cl + DO2 -> ClO + OD 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl, 1.0, i_do2, 1.0, i_clo, 1.0, i_od)

!===========================================================
!      df049 : Cl + DO2 -> DCl + O2 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl, 1.0, i_do2, 1.0, i_dcl, 1.0, i_o2)

!===========================================================
!      df050 : Cl + HDO2 -> DCl + HO2 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl, 1.0, i_hdo2, 1.0, i_dcl, 1.0, i_ho2)

!===========================================================
!      df051 : Cl + HDO2 -> HCl + DO2 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl, 1.0, i_hdo2, 1.0, i_hcl, 1.0, i_do2)

!===========================================================
!      df052 : ClO + DO2 -> DOCl + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clo, 1.0, i_do2, 1.0, i_docl, 1.0, i_o2)

!===========================================================
!      df053 : OH + DOCl -> HDO + ClO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_oh, 1.0, i_docl, 1.0, i_hdo, 1.0, i_clo)

!===========================================================
!      df054 : OD + HOCl -> HDO + ClO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_od, 1.0, i_hocl, 1.0, i_hdo, 1.0, i_clo)

!===========================================================
!      df055 : O + DOCl -> OD + ClO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_o, 1.0, i_docl, 1.0, i_od, 1.0, i_clo)

!===========================================================
!      df056 : Cl2 + D -> DCl + Cl
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl2, 1.0, i_d, 1.0, i_dcl, 1.0, i_cl)

!===========================================================
!      df057 : ClSO2 + D -> SO2 + DCl 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clso2, 1.0, i_d, 1.0, i_so2, 1.0, i_dcl)

!===========================================================
!      df058 : ClCO + OD -> DOCl + CO 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clco, 1.0, i_od, 1.0, i_docl, 1.0, i_co)

!===========================================================
!      df059 : Cl2 + OD -> Cl + DOCl 
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl2, 1.0, i_od, 1.0, i_cl, 1.0, i_docl)

!===========================================================
!      df060 : DCl + H -> HD + Cl
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_dcl, 1.0, i_h, 1.0, i_hd, 1.0, i_cl)

!===========================================================
!      df061 : HCl + D -> HD + Cl
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_hcl, 1.0, i_d, 1.0, i_hd, 1.0, i_cl)

!===========================================================
!      df062 : ClCO + D -> DCl + CO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clco, 1.0, i_d, 1.0, i_dcl, 1.0, i_co)

!===========================================================
!      df063 : Cl + D + M -> DCl + M
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_cl, 1.0, i_d, 1.0, i_dcl, 0.0, i_dummy)


!===========================================================
!      g001 : S + O2 -> SO + O
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_s, 1.0, i_o2, 1.0, i_so, 1.0, i_o)

!===========================================================
!      g002 : S + O3 -> SO + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_s, 1.0, i_o3, 1.0, i_so, 1.0, i_o2)

!===========================================================
!      g003 : SO + O2 -> SO2 + O
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_so, 1.0, i_o2, 1.0, i_so2, 1.0, i_o)

!===========================================================
!      g004 : SO + O3 -> SO2 + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_so, 1.0, i_o3, 1.0, i_so2, 1.0, i_o2)

!===========================================================
!      g005 : SO + OH -> SO2 + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_so, 1.0, i_oh, 1.0, i_so2, 1.0, i_h)

!===========================================================
!      g006 : S + OH -> SO + H
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_s, 1.0, i_oh, 1.0, i_so, 1.0, i_h)

!===========================================================
!      g007 : SO + O + CO2 -> SO2 + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_so, 1.0, i_o, 1.0, i_so2, 0.0, i_dummy)

!===========================================================
!      g008 : SO + HO2 -> SO2 + OH
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_so, 1.0, i_ho2, 1.0, i_so2, 1.0, i_oh)

!===========================================================
!      g009 : SO2 + O + CO2 -> SO3 + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_so2, 1.0, i_o, 1.0, i_so3, 0.0, i_dummy)

!===========================================================
!      g010 : S + O + CO2 -> SO + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_s, 1.0, i_o, 1.0, i_so, 0.0, i_dummy)

!===========================================================
!      g011 : SO3 + H2O -> H2SO4
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_so3, 1.0, i_h2o, 1.0, i_h2so4, 0.0, i_dummy)

!===========================================================
!      g012 : SO + ClO -> SO2 + Cl
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_so, 1.0, i_clo, 1.0, i_so2, 1.0, i_cl)

!===========================================================
!      g013 : SO + SO3 -> SO2 + SO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_so, 1.0, i_so3, 2.0, i_so2, 0.0, i_dummy)

!===========================================================
!      g014 : SO3 + O -> SO2 + O2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_so3, 1.0, i_o, 1.0, i_so2, 1.0, i_o2)

!===========================================================
!      g017 : 0.5 ClCO3 + 0.5 SO -> Cl  
!             0.5 ClCO3 + 0.5 SO -> SO2 + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(0.5, i_clco3, 0.5, i_so, 1.0, i_cl, 0.0, i_dummy)

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(0.5, i_clco3, 0.5, i_so, 1.0, i_so2, 1.0, i_co2)

!===========================================================
!      g018 : S + CO + CO2 -> OCS + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_s, 1.0, i_co, 1.0, i_ocs, 0.0, i_dummy)

!===========================================================
!      g019 : ClCO + S -> OCS + Cl
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_clco, 1.0, i_s, 1.0, i_ocs, 1.0, i_cl)

!===========================================================
!      g020 : SO2 + OH + CO2 -> HSO3 + CO2
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_so2, 1.0, i_oh, 1.0, i_hso3, 0.0, i_dummy)

!===========================================================
!      g021 : HSO3 + O2 -> HO2 + SO3
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_hso3, 1.0, i_o2, 1.0, i_ho2, 1.0, i_so3)

!===========================================================
!      g022 : S + S + CO2 -> S2 + CO2
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_s, 1.0, i_s2, 0.0, i_dummy)

!===========================================================
!      g023 : S2 + O -> SO + S
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_s2, 1.0, i_o, 1.0, i_so, 1.0, i_s)

!===========================================================
!      g024 : S + OCS -> S2 + CO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_s, 1.0, i_ocs, 1.0, i_s2, 1.0, i_co)

!===========================================================
!      g025 : OCS + O -> SO + CO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_ocs, 1.0, i_o, 1.0, i_so, 1.0, i_co)

!===========================================================
!      g026 : S + SO3 -> SO2 + SO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_s, 1.0, i_so3, 1.0, i_so2, 1.0, i_so)

!===========================================================
!      g027 : S + HO2 -> SO + OH
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_s, 1.0, i_ho2, 1.0, i_so, 1.0, i_oh)

!===========================================================
!      g028 : S + ClO -> SO + Cl
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_s, 1.0, i_clo, 1.0, i_so, 1.0, i_cl)

!===========================================================
!      g029: h2so4 + h2o -> so3 + h2o + h2o
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_h2so4, 1.0, i_so3, 1.0, i_h2o)

!===========================================================
!      g032: so + so -> so2 + s
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_so, 1.0, i_so2, 1.0, i_s)

!===========================================================
!      g033 : SO + SO + CO2 -> OSSO_cis + CO2
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_so, 1.0, i_osso_cis, 0.0, i_dummy)

!===========================================================
!      g034: OSSO_cis + CO2 -> SO + SO + CO2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_osso_cis, 2.0, i_so, 0.0, i_dummy)

!===========================================================
!      g035 : SO + SO + CO2 -> OSSO_trans + CO2
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_so, 1.0, i_osso_trans, 0.0, i_dummy)

!===========================================================
!      g036: OSSO_trans + CO2 -> SO + SO + CO2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_osso_trans, 2.0, i_so, 0.0, i_dummy)

!===========================================================
!      g037 : SO + SO + CO2 -> S2O2_cyc + CO2
!===========================================================

nb_reaction_3 = nb_reaction_3 + 1

indice_3(nb_reaction_3) = z3spec(2.0, i_so, 1.0, i_s2o2_cyc, 0.0, i_dummy)

!===========================================================
!      g038: S2O2_cyc + CO2 -> SO + SO + CO2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_s2o2_cyc, 2.0, i_so, 0.0, i_dummy)

!===========================================================
!      g040: OSSO_cis + O -> SO2 + SO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_osso_cis, 1.0, i_o, 1.0, i_so2, 1.0, i_so)

!===========================================================
!      g042: OSSO_trans + O -> SO2 + SO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_osso_trans, 1.0, i_o, 1.0, i_so2, 1.0, i_so)

!===========================================================
!      g044: S2O2_cyc + O -> SO2 + SO
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_s2o2_cyc, 1.0, i_o, 1.0, i_so2, 1.0, i_so)
 
!===========================================================
!      dg001: SO3 + HDO -> HDSO4
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_so3, 1.0, i_hdo, 1.0, i_hdso4, 0.0, i_dummy)

!===========================================================
!      dg002: HDSO4 + H2O -> SO3 + HDO + H2O
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_hdso4, 1.0, i_so3, 1.0, i_hdo)

!===========================================================
!      dg003: HDSO4 + HDO -> SO3 + HDO + HDO
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_hdso4, 1.0, i_so3, 1.0, i_hdo)

!===========================================================
!      dg004 : SO + OD -> SO2 + D
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_so, 1.0, i_od, 1.0, i_so2, 1.0, i_d)

!===========================================================
!      dg005 : S + OD -> SO + D
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_s, 1.0, i_od, 1.0, i_so, 1.0, i_d)

!===========================================================
!      dg006 : SO + DO2 -> SO2 + OD
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_so, 1.0, i_do2, 1.0, i_so2, 1.0, i_od)

!===========================================================
!      dg007 : S + DO2 -> SO + OD
!===========================================================

nb_reaction_4 = nb_reaction_4 + 1

indice_4(nb_reaction_4) = z4spec(1.0, i_s, 1.0, i_do2, 1.0, i_so, 1.0, i_od)


!===========================================================
!      h001: HO2 + ice -> products
!            treated as
!            HO2 -> 0.5 H2O + 0.75 O2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_ho2, 0.5, i_h2o, 0.75, i_o2)

!===========================================================
!      h002: OH + ice -> products
!            treated as
!            OH -> 0.5 H2O + 0.25 O2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_oh, 0.5, i_h2o, 0.25, i_o2)

!===========================================================
!      h003: H2O2 + ice -> products
!            treated as
!            H2O2 -> H2O + 0.5 O2
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_h2o2, 1.0, i_h2o, 0.5, i_o2)

!Only if ion chemistry
if (ok_ionchem) then

!===========================================================
!      i001 : CO2+ + O2 -> O2+ + CO2
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_co2plus, 1.0, i_o2, 1.0, i_o2plus, 1.0, i_co2)

!===========================================================
!      i002 : CO2+ + O -> O+ + CO2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_co2plus, 1.0, i_o, 1.0, i_oplus, 1.0, i_co2)

!===========================================================
!      i003 : CO2+ + O -> O2+ + CO 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_co2plus, 1.0, i_o, 1.0, i_o2plus, 1.0, i_co)

!===========================================================
!      i004 : O2+ + e- -> O + O 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_o2plus, 1.0, i_elec, 2.0, i_o, 0.0, i_dummy)

!===========================================================
!      i005 : O+ + CO2 -> O2+ + CO 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_oplus, 1.0, i_co2, 1.0, i_o2plus, 1.0, i_co)

!===========================================================
!      i006 : CO2+ + e -> CO + O 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_co2plus, 1.0, i_elec, 1.0, i_co, 1.0, i_o)

!===========================================================
!      i007 : CO2+ + NO -> NO+ + CO2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_co2plus, 1.0, i_no, 1.0, i_noplus, 1.0, i_co2)

!===========================================================
!      i008 : O2+ + NO -> NO+ + O2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_o2plus, 1.0, i_no, 1.0, i_noplus, 1.0, i_o2)

!===========================================================
!      i009 : O2+ + N2 -> NO+ + NO 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_o2plus, 1.0, i_n2, 1.0, i_noplus, 1.0, i_no)

!===========================================================
!      i010 : O2+ + N -> NO+ + O 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_o2plus, 1.0, i_n, 1.0, i_noplus, 1.0, i_o)

!===========================================================
!      i011 : O+ + N2 -> NO+ + N 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_oplus, 1.0, i_n2, 1.0, i_noplus, 1.0, i_n)

!===========================================================
!      i012 : NO+ + e -> N + O 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_noplus, 1.0, i_elec, 1.0, i_n, 1.0, i_o)

!===========================================================
!      i013 : CO+ + CO2 -> CO2+ + CO 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_coplus, 1.0, i_co2, 1.0, i_co2plus, 1.0, i_co)

!===========================================================
!      i014 : CO+ + O -> O+ + CO 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_coplus, 1.0, i_o, 1.0, i_oplus, 1.0, i_co)

!===========================================================
!      i015 : C+ + CO2 -> CO+ + CO 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_cplus, 1.0, i_co2, 1.0, i_coplus, 1.0, i_co)

!===========================================================

!      i016 : N2+ + CO2 -> CO2+ + N2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_n2plus, 1.0, i_co2, 1.0, i_co2plus, 1.0, i_n2)

!===========================================================
!      i017 : N2+ + O -> NO+ + N 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_n2plus, 1.0, i_o, 1.0, i_noplus, 1.0, i_n)

!===========================================================
!      i018 : N2+ + CO -> CO+ + N2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_n2plus, 1.0, i_co, 1.0, i_coplus, 1.0, i_n2)

!===========================================================
!      i019 : N2+ + e -> N + N 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_n2plus, 1.0, i_elec, 2.0, i_n, 0.0, i_dummy)

!===========================================================
!      i020 : N2+ + O -> O+ + N2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_n2plus, 1.0, i_o, 1.0, i_oplus, 1.0, i_n2)

!===========================================================
!      i021 : N+ + CO2 -> CO2+ + N 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_nplus, 1.0, i_co2, 1.0, i_co2plus, 1.0, i_n)

!===========================================================
!      i022 : CO+ + H -> H+ + CO 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_coplus, 1.0, i_h, 1.0, i_hplus, 1.0, i_co)

!===========================================================
!      i023 : O+ + H -> H+ + O 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_oplus, 1.0, i_h, 1.0, i_hplus, 1.0, i_o)

!===========================================================
!      i024 : H+ + O -> O+ + H 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_hplus, 1.0, i_o, 1.0, i_oplus, 1.0, i_h)

!===========================================================
!      i025 : CO2+ + H2 -> HCO2+ + H 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_co2plus, 1.0, i_h2, 1.0, i_hco2plus, 1.0, i_h)

!===========================================================
!      i026 : spare slot (reaction rate set to zero)
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_hco2plus, 1.0, i_elec, 1.0, i_h, 1.0, i_co2)

!===========================================================
!      i027 : HCO2+ + e -> H + O + CO 
!===========================================================
!We divide this reaction in two

!0.5HCO2+ + 0.5e -> H

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(.5, i_hco2plus, 0.5, i_elec, 1.0, i_h, 0.0, i_dummy)

!0.5 HCO2+ + 0.5 e -> O + CO

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(0.5, i_hco2plus, 0.5, i_elec, 1.0, i_o, 1.0, i_co)

!===========================================================
!      i029 : HCO2+ + e -> OH + CO 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_hco2plus, 1.0, i_elec, 1.0, i_oh, 1.0, i_co)


!===========================================================
!      i030 : HCO2+ + e -> H + CO2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_hco2plus, 1.0, i_elec, 1.0, i_h, 1.0, i_co2)


!===========================================================
!      i031 : HCO2+ + O -> HCO+ + O2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(1.0, i_hco2plus, 1.0, i_o, 1.0, i_hcoplus, 1.0, i_o2)


!===========================================================
!      i032 : HCO2+ + CO -> HCO+ + CO2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_hco2plus, 1.0, i_co, 1.0, i_hcoplus, 1.0, i_co2)


!===========================================================
!      i033 : H+ + CO2 -> HCO+ + O 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_hplus, 1.0, i_co2, 1.0, i_hcoplus, 1.0, i_o)


!===========================================================
!      i034 : CO2+ + H -> HCO+ + O 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_co2plus, 1.0, i_h, 1.0, i_hcoplus, 1.0, i_o)


!===========================================================
!      i035 : CO+ + H2 -> HCO+ + H 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_coplus, 1.0, i_h2, 1.0, i_hcoplus, 1.0, i_h)


!===========================================================
!      i036 : HCO+ + e- -> CO + H 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_hcoplus, 1.0, i_elec, 1.0, i_co, 1.0, i_h)

!===========================================================
!      i037 : CO2+ + H2O -> H2O+ + CO2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_co2plus, 1.0, i_h2o, 1.0, i_h2oplus, 1.0, i_co2)

!===========================================================
!      i038 : CO+ + H2O -> H2O+ + CO 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_coplus, 1.0, i_h2o, 1.0, i_h2oplus, 1.0, i_co)

!===========================================================
!      i039 : O+ + H2O -> H2O+ + O 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_oplus, 1.0, i_h2o, 1.0, i_h2oplus, 1.0, i_o)

!===========================================================
!      i040 : N2+ + H2O -> H2O+ + N2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_n2plus, 1.0, i_h2o, 1.0, i_h2oplus, 1.0, i_n2)

!===========================================================
!      i041 : N+ + H2O -> H2O+ + N 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_nplus, 1.0, i_h2o, 1.0, i_h2oplus, 1.0, i_n)

!===========================================================
!      i042 : H+ + H2O -> H2O+ + H 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_hplus, 1.0, i_h2o, 1.0, i_h2oplus, 1.0, i_h)

!===========================================================
!      i043 : H2O+ + O2 -> O2+ + H2O 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h2oplus, 1.0, i_o2, 1.0, i_o2plus, 1.0, i_h2o)

!===========================================================
!      i044 : H2O+ + CO -> HCO+ + OH 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h2oplus, 1.0, i_co, 1.0, i_hcoplus, 1.0, i_oh) 

!===========================================================
!      i045 : H2O+ + O -> O2+ + H2 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h2oplus, 1.0, i_o, 1.0, i_o2plus, 1.0, i_h2)

!===========================================================
!      i046 : H2O+ + NO -> NO+ + H2O 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h2oplus, 1.0, i_no, 1.0, i_noplus, 1.0, i_h2o)

!===========================================================
!      i047 : H2O+ + e- -> H + H + O 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h2oplus, 1.0, i_elec, 2.0, i_h, 1.0, i_o)

!===========================================================
!      i048 : H2O+ + e- -> H + OH 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h2oplus, 1.0, i_elec, 1.0, i_h, 1.0, i_oh)

!===========================================================
!      i049 : H2O+ + e- -> H2 + O 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h2oplus, 1.0, i_elec, 1.0, i_h2, 1.0, i_o)

!===========================================================
!      i050 : H2O+ + H2O -> H3O+ + OH 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h2oplus, 1.0, i_h2o, 1.0, i_h3oplus, 1.0, i_oh)

!===========================================================
!      i051 : H2O+ + H2 -> H3O+ + H 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h2oplus, 1.0, i_h2, 1.0, i_h3oplus, 1.0, i_h)

!===========================================================
!      i052 : HCO+ + H2O -> H3O+ + CO 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_hcoplus, 1.0, i_h2o, 1.0, i_h3oplus, 1.0, i_co)

!===========================================================
!      i053: H3O+ + e -> OH + H + H
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h3oplus, 1.0, i_elec, 1.0, i_oh, 2.0, i_h)

!===========================================================
!      i054: H3O+ + e -> H2O + H
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h3oplus, 1.0, i_elec, 1.0, i_h2o, 1.0, i_h)

!===========================================================
!      i055: H3O+ + e -> HO + H2
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h3oplus, 1.0, i_elec, 1.0, i_oh, 1.0, i_h2)

!===========================================================
!      i056: H3O+ + e -> O + H2 + H
!===========================================================
!We divide this reaction in two

!0.5H3O+ + 0.5e -> O

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(0.5, i_h3oplus, 0.5, i_elec, 1.0, i_o, 0.0, i_dummy)

!0.5H3O+ + 0.5e -> H2 + H

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(0.5, i_h3oplus, 0.5, i_elec, 1.0, i_h2, 1.0, i_h)

!===========================================================
!      i057: O+ + H2 -> OH+ + H
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_oplus, 1.0, i_h2, 1.0, i_ohplus, 1.0, i_h)

!===========================================================
!      i058: OH+ + O -> O2+ + H
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_ohplus, 1.0, i_o, 1.0, i_o2plus, 1.0, i_h)

!===========================================================
!      i059: OH+ + CO2 -> HCO2+ + O
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_ohplus, 1.0, i_co2, 1.0, i_hco2plus, 1.0, i_o)

!===========================================================
!      i060: OH+ + CO -> HCO+ + O
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_ohplus, 1.0, i_co, 1.0, i_hcoplus, 1.0, i_o)

!===========================================================
!      i061: OH+ + NO -> NO+ + OH 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_ohplus, 1.0, i_no, 1.0, i_noplus, 1.0, i_oh)

!===========================================================
!      i062: OH+ + H2 -> H2O+ + H 
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_ohplus, 1.0, i_h2, 1.0, i_h2oplus, 1.0, i_h)

!===========================================================
!      i063: OH+ + O2 -> O2+ + OH
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_ohplus, 1.0, i_o2, 1.0, i_o2plus, 1.0, i_oh)

end if    !ok_ionchem
 
if(ok_ionchem) then
   
!===========================================================
!      di001: CO2+ + HD -> HCO2+ + D
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_co2plus, 1.0, i_hd, 1.0, i_hco2plus, 1.0, i_d)

!===========================================================
!      di002: CO2+ + HD -> DCO2+ + H
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_co2plus, 1.0, i_hd, 1.0, i_dco2plus, 1.0, i_h)

!===========================================================
!      di003: CO+ + HD -> HCO+ + D
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_coplus, 1.0, i_hd, 1.0, i_hcoplus, 1.0, i_d)

!===========================================================
!      di004: CO+ + HD -> DCO+ + H
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_coplus, 1.0, i_hd, 1.0, i_dcoplus, 1.0, i_h)

!===========================================================
!      di005: H2O+ + HD -> H3O+ + D
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h2oplus, 1.0, i_hd, 1.0, i_h3oplus, 1.0, i_d)

!===========================================================
!      di007: O+ + HD -> OH+ + D
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_oplus, 1.0, i_hd, 1.0, i_ohplus, 1.0, i_d)

!===========================================================
!      di008: O+ + HD -> OD+ + H
!===========================================================
   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_oplus, 1.0, i_hd, 1.0, i_odplus, 1.0, i_h)

!===========================================================
!      di009: OH+ + HD -> H2O+ + D
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_ohplus, 1.0, i_hd, 1.0, i_h2oplus, 1.0, i_d)

!===========================================================
!      di010: OH+ + HD -> HDO+ + H
!===========================================================
   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_ohplus, 1.0, i_hd, 1.0, i_hdoplus, 1.0, i_h)

!===========================================================
!      di011: DCO2+ + e -> D + CO2
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_dco2plus, 1.0, i_elec, 1.0, i_d, 1.0, i_co2)

!===========================================================
!      di012: DCO2+ + e -> D + O + CO
!===========================================================
!We divide this reaction in two

!0.5DCO2+ + 0.5e -> D

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(.5, i_dco2plus, 0.5, i_elec, 1.0, i_d, 0.0, i_dummy)

!0.5 DCO2+ + 0.5 e -> O + CO

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(0.5, i_dco2plus, 0.5, i_elec, 1.0, i_o, 1.0, i_co)

!===========================================================
!      di013: DCO2+ + e -> OD + CO
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_dco2plus, 1.0, i_elec, 1.0, i_od, 1.0, i_co)

!===========================================================
!      di014: DCO+ + e -> D + CO
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_dcoplus, 1.0, i_elec, 1.0, i_d, 1.0, i_co)

!===========================================================
!      di015: DCO2+ + O -> DCO+ + O2
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_dco2plus, 1.0, i_o, 1.0, i_dcoplus, 1.0, i_o2)

!===========================================================
!      di016: DCO2+ + CO -> DCO+ + CO2
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_dco2plus, 1.0, i_co, 1.0, i_dcoplus, 1.0, i_co2)

!===========================================================
!      di017: CO2+ + D -> DCO+ + O
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_co2plus, 1.0, i_d, 1.0, i_dcoplus, 1.0, i_o)

!===========================================================
!      di018: CO2+ + HDO -> HDO+ + CO2
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_co2plus, 1.0, i_hdo, 1.0, i_hdoplus, 1.0, i_co2)

!===========================================================
!      di019: CO+ + HDO -> HDO+ + CO
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_coplus, 1.0, i_hdo, 1.0, i_hdoplus, 1.0, i_co)

!===========================================================
!      di020: O+ + HDO -> HDO+ + O
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_oplus, 1.0, i_hdo, 1.0, i_hdoplus, 1.0, i_o)

!===========================================================
!      di021: N2+ + HDO -> HDO+ + N2
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_n2plus, 1.0, i_hdo, 1.0, i_hdoplus, 1.0, i_n2)

!===========================================================
!      di022: N+ + HDO -> HDO+ + N
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_nplus, 1.0, i_hdo, 1.0, i_hdoplus, 1.0, i_n)

!===========================================================
!      di023: H+ + HDO -> HDO+ + H
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_hplus, 1.0, i_hdo, 1.0, i_hdoplus, 1.0, i_h)

!===========================================================
!      di024: HDO+ + O2 -> HDO + O2+
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_hdoplus, 1.0, i_o2, 1.0, i_hdo, 1.0, i_o2plus)

!===========================================================
!      di025: HDO+ + CO -> DCO+ + OH
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_hdoplus, 1.0, i_co, 1.0, i_dcoplus, 1.0, i_oh)

!===========================================================
!      di026: HDO+ + CO -> HCO+ + OD
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_hdoplus, 1.0, i_co, 1.0, i_hcoplus, 1.0, i_od)

!===========================================================
!      di027: HDO+ + O -> O2+ + HD
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_hdoplus, 1.0, i_o, 1.0, i_o2plus, 1.0, i_hd)

!===========================================================
!      di028: HDO+ + NO -> NO+ + HDO
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_hdoplus, 1.0, i_no, 1.0, i_noplus, 1.0, i_hdo)

!===========================================================
!      di029: HDO+ + e -> H + D + O 
!===========================================================
!We divide this reaction in two

!0.5 HDO+ + 0.5e -> D

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(.5, i_hdoplus, 0.5, i_elec, 1.0, i_d, 0.0, i_dummy)

!0.5 H2O+ + 0.5 e -> O + H

   nb_reaction_4 = nb_reaction_4 + 1

   indice_4(nb_reaction_4) = z4spec(0.5, i_hdoplus, 0.5, i_elec, 1.0, i_o, 1.0, i_h)

!===========================================================
!      di030: HDO+ + e- -> D + OH
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_hdoplus, 1.0, i_elec, 1.0, i_d, 1.0, i_oh)

!===========================================================
!      di031: HDO+ + e- -> H + OD
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_hdoplus, 1.0, i_elec, 1.0, i_h, 1.0, i_od)

!===========================================================
!      di032: HDO+ + e- -> HD + O
!===========================================================

   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_hdoplus, 1.0, i_elec, 1.0, i_hd, 1.0, i_o)

!===========================================================
!      di033: D+ + O -> O+ + D
!===========================================================
   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_dplus, 1.0, i_o, 1.0, i_oplus, 1.0, i_d)

!===========================================================
!      di034: D+ + CO2 -> DCO+ + O
!===========================================================
   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_dplus, 1.0, i_co2, 1.0, i_dcoplus, 1.0, i_o)

!===========================================================
!      di035: D+ + H2O -> HDO+ + H
!===========================================================
   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_dplus, 1.0, i_h2o, 1.0, i_hdoplus, 1.0, i_h)

!===========================================================
!      di036: OD+ + O -> O2+ + D
!===========================================================
   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_odplus, 1.0, i_o, 1.0, i_o2plus, 1.0, i_d)

!===========================================================
!      di037: OD+ + CO2 -> DCO2+ + O
!===========================================================
   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_odplus, 1.0, i_co2, 1.0, i_dco2plus, 1.0, i_o)

!===========================================================
!      di038: OD+ + CO -> DCO+ + O
!===========================================================
   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_odplus, 1.0, i_co, 1.0, i_dcoplus, 1.0, i_o)

!===========================================================
!      di039: OD+ + NO -> NO+ + OD 
!===========================================================
   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_odplus, 1.0, i_no, 1.0, i_noplus, 1.0, i_od)

!===========================================================
!      di040: OD+ + H2 -> HDO+ + H 
!===========================================================
   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_odplus, 1.0, i_h2, 1.0, i_hdoplus, 1.0, i_h)

!===========================================================
!      di041: OD+ + O2 -> O2+ + OD
!===========================================================
   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_odplus, 1.0, i_o2, 1.0, i_o2plus, 1.0, i_od)

!===========================================================
!      di042: CO+ + D -> D+ + CO 
!===========================================================
   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_coplus, 1.0, i_d, 1.0, i_dplus, 1.0, i_co)

!===========================================================
!      di043: O+ + D -> D+ + O 
!===========================================================
   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_oplus, 1.0, i_d, 1.0, i_dplus, 1.0, i_o)

!===========================================================
!      di044: H2O+ + HD -> H2DO+ + H
!===========================================================
   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h2oplus, 1.0, i_hd, 1.0, i_h2doplus, 1.0, i_h)

!===========================================================
!      di045: HDO+ + H2O -> H2DO+ + OH
!===========================================================
   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_hdoplus, 1.0, i_h2o, 1.0, i_h2doplus, 1.0, i_oh)

!===========================================================
!      di046: H2DO+ + e- -> HDO + H
!===========================================================
   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h2doplus, 1.0, i_elec, 1.0, i_hdo, 1.0, i_h)

!===========================================================
!      di047: H2DO+ + e- -> H2O + D
!===========================================================
   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h2doplus, 1.0, i_elec, 1.0, i_h2o, 1.0, i_d)

!===========================================================
!      di048: H2DO+ + e- -> OH + HD
!===========================================================
   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h2doplus, 1.0, i_elec, 1.0, i_oh, 1.0, i_hd)

!===========================================================
!      di049: H2DO+ + e- -> OD + H2
!===========================================================
   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h2doplus, 1.0, i_elec, 1.0, i_od, 1.0, i_h2)

!===========================================================
!      di050: H2DO+ + e- -> OD + H + H
!===========================================================
   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(1.0, i_h2doplus, 1.0, i_elec, 1.0, i_od, 2.0, i_h)

!===========================================================
!      di051: H2DO+ + e- -> OH + D + H
!===========================================================
! 0.5 H2DO+ + 0.5e- -> OH
   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(0.5, i_h2doplus, 0.5, i_elec, 1.0, i_oh, 0.0, i_dummy)

! 0.5 H2DO+ + 0.5e- -> D + H
   nb_reaction_4 = nb_reaction_4 + 1
   indice_4(nb_reaction_4) = z4spec(0.5, i_h2doplus, 0.5, i_elec, 1.0, i_d, 1.0, i_h)


end if   !ionchem.and.deutchem

!===========================================================
!      j001: O2(Dg) + (CO2 or O) -> O2 + (CO2 or O)
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_o2dg, 1.0, i_o2, 0.0, i_dummy)

!===========================================================
!      j002: O2(Dg) -> O2 + hv
!===========================================================

nb_phot = nb_phot + 1

indice_phot(nb_phot) = z3spec(1.0, i_o2dg, 1.0, i_o2, 0.0, i_dummy)

!===========================================================
!  check dimensions 
!===========================================================

print*, 'nb_phot       = ', nb_phot
print*, 'nb_phot_max   = ', nb_phot_max
print*, 'nb_reaction_4 = ', nb_reaction_4
print*, 'nb_reaction_4_max = ', nb_reaction_4_max
print*, 'nb_reaction_3 = ', nb_reaction_3
print*, 'nb_reaction_3_max = ', nb_reaction_3_max

!print*, 'check dimension'
if ((nb_phot /= nb_phot_max)             .or.  &
    (nb_reaction_3 /= nb_reaction_3_max) .or.  &
    (nb_reaction_4 /= nb_reaction_4_max)) then
   print*, 'wrong dimensions in indice' 
   stop
end if  

end subroutine indice

!===========================================================

subroutine phot(nj, nztable, nsza, nso2, sza_input, dist_sol, mumean, rmco2, rmso2,                     &
                jphot, table_colair, table_colso2, table_sza, nz, nb_phot_max, t, p, v_phot)

!===========================================================

use clesphys_mod
implicit none
!#include "clesphys.h"

integer, INTENT(IN) :: nz
integer, INTENT(IN) :: nj, nztable, nsza, nso2

real, INTENT(IN), dimension(nz) :: t, p
real, INTENT(IN), dimension(nz) :: mumean        ! [g/mol]
real, INTENT(IN), dimension(nso2,nsza,nztable,nj) :: jphot
real, INTENT(IN), dimension(nso2,nztable) :: table_colso2
real, INTENT(IN), dimension(nztable) :: table_colair
real, INTENT(IN), dimension(nsza) :: table_sza
real, INTENT(IN), dimension(nz) :: rmco2, rmso2
real, INTENT(IN) :: sza_input, dist_sol

real, dimension(nz,nj) :: j
real, dimension(nz) :: coef, col, colso2
real, dimension(nso2) :: colref
real, dimension(2,2,2) :: poids
real :: cicol, cisza, ciso2
real :: avogadro, gvenus, dp

integer :: indcol, indsza, indso2
integer :: isza, iz, i, iso2, ij
integer :: nb_phot_max
real, dimension(nz,nb_phot_max), INTENT(INOUT) :: v_phot

!mugaz    = 43.44E-3
avogadro = 6.022E+23
gvenus   = 8.87

! day/night test

if (sza_input <= 95.) then      ! day

! interpolation in solar zenith angle

indsza = nsza - 1
do isza = 1,nsza
   if (table_sza(isza) >= sza_input) then
      indsza = min(indsza,isza - 1)
      indsza = max(indsza, 1)
   end if
end do

cisza = (sza_input - table_sza(indsza))                       &
       /(table_sza(indsza + 1) - table_sza(indsza))

!    print*, 'indsza    = ', indsza
!    print*, 'table_sza = ', table_sza(indsza)
!    print*, 'cisza     = ', cisza

! co2 and so2 columns

coef(nz)   = avogadro/(gvenus*mumean(nz)*1.E-3)*1.E-4
col(nz)    = coef(nz)*rmco2(nz)*p(nz)*100.
colso2(nz) = coef(nz)*rmso2(nz)*p(nz)*100.

do iz = nz-1, 1, -1
!   print*,"L2490 new_photochemistry", iz,mumean(iz)
   dp = (p(iz) - p(iz+1))*100.
   coef(iz)   = avogadro/(gvenus*mumean(iz)*1.E-3)*1.E-4
   col(iz)    = col(iz+1) + coef(iz)*(rmco2(iz+1) + rmco2(iz))*0.5*dp
   col(iz)    = min(col(iz), table_colair(1))
   colso2(iz) = colso2(iz+1) + coef(iz)*(rmso2(iz+1) + rmso2(iz))*0.5*dp
   colso2(iz) = min(colso2(iz), table_colso2(nso2,1))
end do

! loop over altitude

do iz = 1,nz

! interpolation in co2 column

   indcol = nztable - 1
   cicol  = 0.

   do i = 1,nztable-1
      if (table_colair(i) < col(iz)) then
         cicol = (log(col(iz)) - log(table_colair(i)))           &
                /(log(table_colair(i-1)) - log(table_colair(i)))
         indcol = i - 1
         exit
      end if
   end do

! interpolation in so2 column

! initialize indso2 and ciso2 in case colref is never larger
! than the gcm so2 column.

   indso2 = nso2 - 1
   ciso2 = 1.

! search for the index indso2 between which interpolate

   do iso2 = 1,nso2 
      colref(iso2) = cicol*table_colso2(iso2,indcol)               &
                   + (1.-cicol)*table_colso2(iso2,indcol+1)
      if (colref(iso2) > colso2(iz)) then
         ciso2 = (colso2(iz) - colref(iso2-1))                     &
                /(colref(iso2) - colref(iso2-1))
         indso2 = iso2 - 1
         exit
      end if
   end do

! 4-dimensional interpolation weights

! poids(so2,sza_input,co2)

   poids(1,1,1) = (1.-ciso2)*(1.-cisza)*    cicol 
   poids(1,1,2) = (1.-ciso2)*(1.-cisza)*(1.-cicol)
   poids(1,2,1) = (1.-ciso2)*    cisza *    cicol 
   poids(1,2,2) = (1.-ciso2)*    cisza *(1.-cicol)
   poids(2,1,1) =     ciso2 *(1.-cisza)*    cicol 
   poids(2,1,2) =     ciso2 *(1.-cisza)*(1.-cicol)
   poids(2,2,1) =     ciso2 *    cisza *    cicol 
   poids(2,2,2) =     ciso2 *    cisza *(1.-cicol)

! 4-dimensional interpolation in the lookup table

   do ij = 1,nj
      j(iz,ij) =                                          &
      poids(1,1,1)*jphot(indso2  ,indsza  ,indcol  ,ij)   &
    + poids(1,1,2)*jphot(indso2  ,indsza  ,indcol+1,ij)   &
    + poids(1,2,1)*jphot(indso2  ,indsza+1,indcol  ,ij)   &
    + poids(1,2,2)*jphot(indso2  ,indsza+1,indcol+1,ij)   &
    + poids(2,1,1)*jphot(indso2+1,indsza  ,indcol  ,ij)   &
    + poids(2,1,2)*jphot(indso2+1,indsza  ,indcol+1,ij)   &
    + poids(2,2,1)*jphot(indso2+1,indsza+1,indcol  ,ij)   &
    + poids(2,2,2)*jphot(indso2+1,indsza+1,indcol+1,ij)
   end do

end do           ! end of loop over altitude

else             ! night
   j(:,:) = 0.
end if

! photodissociation rates numbering in the lookup table

!    1     o2 + hv     -> o + o
!    2     o2 + hv     -> o + o(1d)
!    3     co2 + hv    -> co + o
!    4     co2 + hv    -> co + o(1d)
!    5     o3 + hv     -> o2(Dg) + o(1d)
!    6     o3 + hv     -> o2 + o
!    7     h2o + hv    -> h + oh
!    8     ho2 + hv    -> oh + o
!    9     h2o2 + hv   -> oh + oh
!    10    hcl + hv    -> h + cl
!    11    cl2 + hv    -> cl + cl
!    12    hocl + hv   -> oh + cl
!    13    so2 + hv    -> so + o
!    14    so + hv     -> s + o
!    15    so3 + hv    -> so2 + o
!    16    clo + hv    -> cl + o
!    17    ocs + hv    -> co + s
!    18    cocl2 + hv  -> cl + cl + co
!    19    h2so4 + hv  -> so3 + h2o
!    20    no2 + hv    -> no + o
!    21    no + hv     -> n + o  
!    22    n2 + hv     -> n + n 

! fill v_phot array

do ij = 1,nj
   v_phot(:,ij) = j(:,ij)
end do

!!! TEST: artificial increase of CO2 photodissociation
if (tuneupperatm) then
!-- TuneA
!   v_phot(65:78,4) = v_phot(65:78,4)*10.
!   v_phot(60:64,4) = v_phot(60:64,4)*3.
!   v_phot(55:59,3) = v_phot(55:59,3)*2.
!--
!-- TuneB
!   v_phot(65:78,4) = v_phot(65:78,4)*10.
!   v_phot(55:59,3) = v_phot(55:59,3)*5.
!--
!-- TuneC
! VCD 1.1 tuning
!   v_phot(65:78,4) = v_phot(65:78,4)*10.
!   v_phot(52:59,3) = v_phot(52:59,3)*5.
!--
!-- TuneE
! VCD 2.0 tuning
    v_phot(65:nz,4) = v_phot(65:nz,4)*10. ! CO2 + hv ==> O(1D) + CO
!--
!   v_phot(:,4) = v_phot(:,4)*10.
!do ij=3,4
!   v_phot(:,ij) = v_phot(:,ij)*10.
!end do
endif
!!!!!!!!!!!!!!!

!PRINT*,'sza_input: ',sza_input
!IF (sza_input.le.40.5 .AND. sza_input.gt.39.5) THEN
!open(200, form = 'formatted')
!100    format(e20.6)
!write(200,100)(v_phot(:,19))
!stop
!END IF

end subroutine phot

!======================================================================

 subroutine krates(hetero_ice,hetero_dust,ok_ionchem, nphotion,        &
                   nz, nesp, c, conc, t, t_elect, p,                   &
                   nb_phot_max, nb_reaction_3_max, nb_reaction_4_max,  &
                   tuneupperatm, v_3, v_4, v_phot,sza_input,           &
                   ind_norec, ind_orec)
 
!================================================================
! compute reaction rates                                        !
!----------------------------------------------------------------
! reaction               type                array              !
!----------------------------------------------------------------
! A + B    --> C + D     bimolecular         v_4                !
! A + A    --> B + C     quadratic           v_3                !
! A + C    --> B + C     quenching           v_phot             !
! A + ice  --> B + C     heterogeneous       v_phot             !
!================================================================

 USE chemparam_mod
 USE photolysis_mod, only : nphot
implicit none

!----------------------------------------------------------------------
!     input
!----------------------------------------------------------------------

integer, INTENT(IN) :: nesp, nz
real, INTENT(IN)    :: sza_input  ! [degree]
real, INTENT(IN), dimension(nz)  :: t_elect, t, p, conc
logical, INTENT(IN) :: hetero_ice, hetero_dust, tuneupperatm

real, INTENT(IN), dimension(nz,nesp) :: c

integer, intent(in) :: nb_reaction_3_max ! number of quadratic reactions
integer, intent(in) :: nb_reaction_4_max ! number of bimolecular reactions
integer, intent(in) :: nb_phot_max       ! number of reactions treated numerically as photodissociations
integer, intent(in) :: nphotion          ! number of photoionizations
logical, intent(in) :: ok_ionchem        ! if .true. then ionchem reaction used

!----------------------------------------------------------------------
!     output
!----------------------------------------------------------------------

real (kind = 8), dimension(nz,       nb_phot_max) :: v_phot
real (kind = 8), dimension(nz, nb_reaction_3_max) :: v_3
real (kind = 8), dimension(nz, nb_reaction_4_max) :: v_4
integer :: ind_norec
integer :: ind_orec

!----------------------------------------------------------------------
!     local
!----------------------------------------------------------------------

integer :: iz
real    :: ak0, ak1, xpo, rate, rate1, rate2, pi, gam, epsil
real    :: k1a0, k1b0, k1ainf, k1a, k1b, k0, kinf, kf, kint, kca, fc, fx, x, y
integer :: nb_phot, nb_reaction_3, nb_reaction_4
real, dimension(nz)  :: surfice1d, surfdust1d
real, dimension(nz) :: deq
real, dimension(nz) :: a001, a002, a003,                           &
                       b001, b002, b003, b004, b005, b006, b007,   &
                       b008, b009,                                 &
                       db001, db002, db003,                        &
                       c001, c002, c003, c004, c005, c006, c007,   &
                       c008, c009, c010, c011, c012, c013, c014,   &
                       c015, c016, c017, c018,                     &
                       d001, d002, d003, d004, d005, d006, d007,   &
                       d008, d009, d010, d011, d012, d013, d014,   &
                       e001, e002, e003, e004, e005, e006, e007,   &
                       e008, e009, e010, e011, e012, e013, e014,   &
                       e015, e016, e017, e018, e019, e020, e021,   &
                       e022, e023, e024, e025, e026, e027, e028,   &
                       e029, e030, e031, e032, e033, e034, e035,   &
                       e036, e037, e038, e039, e040, e041, e042,   &
                       e043,                                       &
                       f001, f002, f003, f004, f005, f006, f007,   &
                       f008, f009, f010, f011, f012, f013, f014,   &
                       f015, f016, f017, f018, f019, f020, f021,   &
                       f022, f023, f024, f025, f026, f027, f028,   &
                       f029, f030, f031, f032, f033, f034, f035,   &
                       f036, f037, f038, f039, f040,               &
                       df001, df002, df003, df004, df005, df006,   &
                       df007, df008, df009, df010, df011, df012,   &
                       df013, df014, df015, df016, df017, df018,   &
                       df019, df020, df021, df022, df023, df024,   &
                       df025, df026, df027, df028, df029, df030,   &
                       df031, df032, df033, df034, df035, df036,   &
                       df037, df038, df039, df040, df041, df042,   &
                       df043, df044, df045, df046, df047, df048,   &
                       df049, df050, df051, df052, df053, df054,   &
                       df055, df056, df057, df058, df059, df060,   &
                       df061, df062, df063,                        &
                       g001, g002, g003, g004, g005, g006, g007,   &
                       g008, g009, g010, g011, g012, g013, g014,   &
                       g015, g016, g017, g018, g019, g020, g021,   &
                       g022, g023, g024, g025, g026, g027, g028,   &
                       g029, g030, g031, g032, g033, g034, g035,   &
                       dg001, dg002, dg003, dg004, dg005, dg006,   &
                       dg007,                                      &
                       g036, g037, g038, g039, g040, g041, g042,   &
                       g043, g044, g045,                           &
                       h001, h002, h003,                           &
                       i001, i002, i003, i004, i005, i006,         &
                       i007, i008, i009, i010, i011, i012,         &
                       i013, i014, i015, i016, i017, i018, i019,   &
                       i020, i021, i022, i023, i024, i025, i026,   &
                       i027, i028, i029, i030, i031, i032, i033,   &
                       i034, i035, i036, i037, i038, i039, i040,   &
                       i041, i042, i043, i044, i045, i046, i047,   &
                       i048, i049, i050, i051, i052, i053, i054,   &
                       i055, i056, i057, i058, i059, i060, i061,   &
                       di001, di002, di003, di004, di005, di006,   &
                       di007, di008, di009, di010, di011, di012,   &
                       di013, di014, di015, di016, di017, di018,   &
                       di019, di020, di021, di022, di023, di024,   &
                       di025, di026, di027, di028, di029, di030,   &
                       di031, di032, di033, di034, di035, di036,   &
                       di037, di038, di039, di040, di041, di042,   &
                       di043, di044, di045, di046, di047, di048,   &
                       di049, di050, di051, di052, di053, di054,   &
                       i062, i063,                                 &
                       j001, j002                                  
!----------------------------------------------------------------------
!     initialisation
!----------------------------------------------------------------------

      pi = acos(-1.)
      
      nb_phot       = nphot + nphotion ! initialised to the number of photolysis + number of photoionization rates
      nb_reaction_3 = 0
      nb_reaction_4 = 0

!----------------------------------------------------------------------
!        reactions avec ox
!----------------------------------------------------------------------

!---  a001: o + o2 + m -> o3 + m

!     jpl 2019
 
!     co2/n2 efficiency as a third body = 2.075
!     from sehested et al., j. geophys. res., 100, 1995.

      a001(:) = 6.1e-34*(t(:)/298.)**(-2.4)     &
              *(2.075*c(:,i_co2) + 1.0*(conc(:) - c(:,i_co2))) 

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = a001(:)

!---  a002: o + o + m -> o2(delta_g) + m
 
!     baulch et al., 1976 (confirmed by smith and robertson, 2008)

!     epsil : net effective yield 

      epsil = 0.75 ! (crisp et al., 1996; krasnopolsky, 1991)

      a002(:) = 2.76e-34*exp(720./t(:))  &
              *(2.5*c(:,i_co2) + 1.0*(conc(:) - c(:,i_co2)))*epsil
      
      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = a002(:)

      ind_orec = nb_reaction_3

!---  a003: o + o3 -> o2 + o2

!     jpl 2003

      a003(:) = 8.0E-12*exp(-2060./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = a003(:)

!----------------------------------------------------------------------
!        reactions avec o(1d)
!----------------------------------------------------------------------

!---  b001: o(1d) + co2  -> o + co2

!     jpl 2006

      b001(:) = 7.5E-11*exp(115./t(:))

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = b001(:)*c(:,i_co2)

!---  b002: o(1d) + h2o  -> oh + oh

!     jpl 2006
         
      b002(:) = 1.63E-10*exp(60./t(:))
              
      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = b002(:)

!---  b003: o(1d) + h2  -> oh + h

!     jpl 2011      

      b003(:) = 1.2E-10

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = b003(:)

!---  b004: o(1d) + o2  -> o + o2

!     jpl 2006

      b004(:) = 3.3E-11*exp(55./t(:))

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = b004(:)*c(:,i_o2)
            
!---  b005: o(1d) + o3  -> o2 + o2

!     jpl 2003

      b005(:) = 1.2E-10

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = b005(:)
            
!---  b006: o(1d) + o3  -> o2 + o + o

!     jpl 2003

      b006(:) = 1.2E-10

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = b006(:)
 
!---     db001: o(1d) + hdo -> od + oh

      db001(:) = b002(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = db001(:)

!---     db002: o(1d) + hd -> h + od

      !Laurent et al., 1995

      db002(:) = 1.3e-10

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = db002(:)

!---     db003: o(1d) + hd -> d + oh

      !Laurent et al., 1995

      db003(:) = 1.0e-10
      
      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = db003(:)
            
!----------------------------------------------------------------------
!        reactions des hox    
!----------------------------------------------------------------------

!---  c001: o + ho2 -> oh + o2

!     jpl 2003

      c001(:) = 3.0E-11*exp(200./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c001(:)

!---  c002: o + oh -> o2 + h

!     jpl 2011  

      c002(:) = 1.8E-11*exp(180./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c002(:)

!---  c003: h + o3 -> oh + o2

!     jpl 2003

      c003(:) = 1.4E-10*exp(-470./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c003(:)

!---  c004: h + ho2 -> oh + oh

!     jpl 2006

      c004(:) = 7.2E-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c004(:)

!---  c005: h + ho2 -> h2 + o2

!     jpl 2006

      c005(:) = 6.9E-12

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c005(:)

!---  c006: h + ho2 -> h2o + o

!     jpl 2006

      c006(:) = 1.6E-12

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c006(:)

!---  c007: oh + ho2 -> h2o + o2

!     jpl 2003

      c007(:) = 4.8E-11*exp(250./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c007(:)

!---  c008: ho2 + ho2 -> h2o2 + o2

!     jpl 2015

      c008(:) = 3.0E-13*exp(460./t(:))

!     christensen et al., grl, 13, 2002

!     c008(:) = 1.5E-12*exp(19./t(:))

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = c008(:)

!---  c009: oh + h2o2 -> h2o + ho2

!     jpl 2006

      c009(:) = 1.8E-12

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c009(:)

!---  c010: oh + h2 -> h2o + h

!     jpl 2006

      c010(:) = 2.8E-12*exp(-1800./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c010(:)

!---  c011: h + o2 + co2 -> ho2 + co2

!     jpl 2006

!     do iz = 1,nz
!        ak0 = 2.5*4.4E-32*(t(iz)/300.)**(-1.3)
!        ak1 = 4.7E-11*(t(iz)/300.)**(-0.2)

!        rate = (ak0*conc(iz))/(1. + ak0*conc(iz)/ak1)
!        xpo = 1./(1. + alog10((ak0*conc(iz))/ak1)**2)
!        c011(iz) = rate*0.6**xpo
!     end do

!     jpl 2019

      do iz = 1,nz
         ak0 = 2.4*5.3e-32*(t(iz)/298.)**(-1.8)
         ak1 = 9.5e-11*(t(iz)/298.)**(0.4)

         rate = (ak0*conc(iz))/(1. + ak0*conc(iz)/ak1)
         xpo = 1./(1. + alog10((ak0*conc(iz))/ak1)**2)
         c011(iz) = rate*0.6**xpo
      end do

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c011(:)

!---  c012: o + h2o2 -> oh + ho2

!     jpl 2003

      c012(:) = 1.4E-12*exp(-2000./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c012(:)

!---  c013: oh + oh -> h2o + o

!     jpl 2006

      c013(:) = 1.8E-12

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = c013(:)

!---  c014: oh + o3 -> ho2 + o2

!     jpl 2003

      c014(:) = 1.7E-12*exp(-940./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c014(:)

!---  c015: ho2 + o3 -> oh + o2 + o2

!     jpl 2003

      c015(:) = 1.0E-14*exp(-490./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = c015(:)

!---  c016: ho2 + ho2 + co2 -> h2o2 + o2 + co2

!     jpl 2011

      c016(:) = 2.5*2.1E-33*exp(920./t(:))*conc(:)

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = c016(:)

!---  c017: oh + oh + co2 -> h2o2 + co2

!     jpl 2003

      do iz = 1,nz
         ak0 = 2.5*6.9E-31*(t(iz)/300.)**(-1.0)
         ak1 = 2.6E-11*(t(iz)/300.)**(0.0)

         rate = (ak0*conc(iz))/(1. + ak0*conc(iz)/ak1)
         xpo = 1./(1. + alog10((ak0*conc(iz))/ak1)**2)
         c017(iz) = rate*0.6**xpo
      end do

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = c017(:)

!---  c018: h + h + co2 -> h2 + co2

!     baulch et al., 2005

      c018(:) = 2.5*1.8E-30*(t(:)**(-1.0))*conc(:)

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = c018(:)

!----------------------------------------------------------------------
!        reactions des composes azotes
!----------------------------------------------------------------------

!---  d001: no2 + o -> no + o2

!     jpl 2006

!     d001(:) = 5.1e-12*exp(210./t(:))

!     jpl 2019

!     For the sake of simplicity, it is assumed that the association
!     yield (kf) gives the same product as the chemical activation yield
!     (kca). Thus the only products are no + o2. There is no production
!     of no3.

      do iz = 1,nz

!        association

        k0 = 2.5*3.4e-31*(298./t(iz))**(1.6)
        kinf = 2.3e-11*(298./t(iz))**(0.2)

        kf = (kinf*k0*conc(iz)/(kinf + k0*conc(iz))) &
             *0.6**(1. + (log10(k0*conc(iz)/kinf))**2.)**(-1.0)

!       chemical acitvation

         kint = 5.3e-12*exp(200./t(iz))

         kca = kint*(1. - kf/kinf)

!        total : association + chemical activation

         d001(iz) = kf + kca

      end do

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d001(:)

!---  d002: no + o3 -> no2 + o2

!     jpl 2006

      d002(:) = 3.0e-12*exp(-1500./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d002(:)

!---  d003: no + ho2 -> no2 + oh

!     jpl 2011

!     d003(:) = 3.3e-12*exp(270./t(:))

!     jpl 2019

      d003(:) = 3.44e-12*exp(270./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d003(:)


!---  d004: n + no -> n2 + o

!     jpl 2011

      d004(:) = 2.1e-11*exp(100./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d004(:)

!---  d005: n + o2 -> no + o

!     jpl 2011

!     d005(:) = 1.5e-11*exp(-3600./t(:))

!     jpl 2019

      d005(:) = 3.3e-12*exp(-3150./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d005(:)

!---  d006: no2 + h -> no + oh

!     jpl 2011

!     d006(:) = 4.0e-10*exp(-340./t(:))

!     jpl 2019

      d006(:) = 1.35e-10

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d006(:)

!---  d007: n + o -> no

      d007(:) = 1.9e-17*(300./t(:))**(0.5)*exp(1-(0.57/(t(:)**(0.5))))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d007(:)

      ind_norec = nb_reaction_4

!---  d008: n + ho2 -> no + oh

!     brune et al., j. chem. phys., 87, 1983 

      d008(:) = 2.19e-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d008(:)

!---  d009: n + oh -> no + h

!     atkinson et al., j. phys. chem. ref. data, 18, 881, 1989 

      d009(:) = 3.8e-11*exp(85./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d009(:)

!---  d010: n2d + o -> n + o

!     herron, j. phys. chem. ref. data, 1999

      d010(:) = 3.3e-12*exp(-260./t(:))

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = d010(:)*c(:,i_o)

!---  d011: n2d + n2 -> n + n2

!     herron, j. phys. chem. ref. data, 1999

      d011(:) = 1.7e-14

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = d011(:)*c(:,i_n2)

!---  d012: n2d + co2 -> no + co

!     herron, j. phys. chem. ref. data, 1999

      d012(:) = 3.6e-13

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d012(:)

!---  d013: n + o + co2 -> no + co2

!     Campbell & Trush, 1966      

      d013(:) = 2.5 * conc(:) * 1.83e-32 * (298./t(:))**0.5
      
      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = d013(:)

!--- d014: n2d + co -> n + co

!    herron, j. phys. chem. ref. data, 1999

     d014(:) = 1.9e-12

     nb_phot = nb_phot + 1
     v_phot(:,nb_phot) = d014(:)*c(:,i_co)
!----------------------------------------------------------------------
!        reactions des composes carbones
!----------------------------------------------------------------------

!---  e001: oh + co -> co2 + h

!     jpl 2015

!      do iz = 1,nz

!        branch 1 : oh + co -> h + co2

!         rate1 = 1.5e-13*(t(iz)/300.)**(0.0)

!        branch 2 : oh + co + m -> hoco + m

!         ak0 = 5.9e-33*(t(iz)/300.)**(-1.0)
!         ak1 = 1.1e-12*(t(iz)/300.)**(1.3)
!         rate2 = (ak0*conc(iz))/(1. + ak0*conc(iz)/ak1)
!         xpo = 1./(1. + alog10((ak0*conc(iz))/ak1)**2)

!         e001(iz) = rate1 + rate2*0.6**xpo
!      end do

!     joshi et al., 2006

!     do iz = 1,nz
!        k1a0 = 1.34*2.5*conc(iz)                                &
!              *1/(1/(3.62E-26*t(iz)**(-2.739)*exp(-20./t(iz)))  &
!              + 1/(6.48E-33*t(iz)**(0.14)*exp(-57./t(iz))))     ! typo in paper corrected
!        k1b0 = 1.17E-19*t(iz)**(2.053)*exp(139./t(iz))          &
!             + 9.56E-12*t(iz)**(-0.664)*exp(-167./t(iz))
!        k1ainf = 1.52E-17*t(iz)**(1.858)*exp(28.8/t(iz))        &
!               + 4.78E-8*t(iz)**(-1.851)*exp(-318./t(iz))
!        x = k1a0/(k1ainf - k1b0)
!        y = k1b0/(k1ainf - k1b0)
!        fc = 0.628*exp(-1223./t(iz)) + (1. - 0.628)*exp(-39./t(iz))  &
!           + exp(-t(iz)/255.)
!        fx = fc**(1./(1. + (alog(x))**2))                       ! typo in paper corrected
!        k1a = k1a0*((1. + y)/(1. + x))*fx
!        k1b = k1b0*(1./(1.+x))*fx
!        e001(iz) = k1a + k1b
!     end do

!     jpl 2019

      do iz = 1,nz

!       association

        k0 = 2.5*6.9e-33*(298./t(iz))**(2.1)
        kinf = 1.1e-12*(298./t(iz))**(-1.3)

        kf = (kinf*k0*conc(iz)/(kinf + k0*conc(iz))) &
             *0.6**(1. + (log10(k0*conc(iz)/kinf))**2.)**(-1.0)

!       chemical activation

        kint = 1.85e-13*exp(-65/t(iz))

        kca = kint*(1. - kf/kinf)

!       total : association + chemical activation

        e001(iz) = kf + kca

      end do

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = e001(:)

!---  e002: o + co + m -> co2 + m

!     tsang and hampson, 1986.

      e002(:) = 2.5*6.5E-33*exp(-2184./t(:))*conc(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = e002(:)

!----------------------------------------------------------------------
!        reactions des composes chlores
!----------------------------------------------------------------------

!---  f001: hcl + o(1d) -> oh + cl

!     jpl 2011

      f001(:) = 1.0E-10

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f001(:)

!---  f002: hcl + o(1d) -> h + clo

!     jpl 2011

      f002(:) = 3.6E-11
      
      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f002(:)

!---  f003: hcl + o -> oh + cl

!     jpl 2006

      f003(:) = 1.0E-11*exp(-3300./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f003(:)

!---  f004: hcl + oh -> h2o + cl

!     jpl 2009

      f004(:) = 1.8E-12*exp(-250./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f004(:)

!---  f005: clo + o -> cl + o2

!     jpl 2006

      f005(:) = 2.8E-11*exp(85./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f005(:)

!---  f006: clo + oh -> cl + ho2

!     jpl 2006

      f006(:) = 7.4E-12*exp(270./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f006(:)

!---  f007: clo + oh -> hcl + o2

!     jpl 2006

      f007(:) = 6.0E-13*exp(230./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f007(:)

!---  f008: cl + h2 -> hcl + h

!     jpl 2006

      f008(:) = 3.05E-11*exp(-2270./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f008(:)

!---  f009: cl + o3 -> clo + o2

!     jpl 2006

      f009(:) = 2.3E-11*exp(-200./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f009(:)

!---  f010: cl + ho2 -> clo + oh

!     jpl 2009

      f010(:) = 3.6E-11*exp(-375./t(:)) 
      
      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f010(:)

!---  f011: cl + ho2 -> hcl + o2

!     jpl 2009

      f011(:) = 1.4E-11*exp(270./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f011(:)

!---  f012: cl + h2o2 -> hcl + ho2

!     jpl 2006

      f012(:) = 1.1E-11*exp(-980./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f012(:)

!---  f013: cl + co + co2 -> clco + co2

!     jpl 2011 + nicovich et al., j. phys. chem., 1990

      f013(:) = 3.2*1.3E-33*(t(:)/300.)**(-3.8)*conc(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f013(:)

!---  f014: clco + co2 -> cl + co + co2

!     jpl 2011

!     deq(:) = 3.2*3.5E-25*exp(3730./t(:))

!     mills, 1998

      deq(:) = 1.6E-25*exp(4000./t(:))

      f014(:) = f013(:)/(deq(:)*conc(:))

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = f014(:)*conc(:)

!     do iz = 1, nz
!     print*, z(iz), t(iz), f013(iz), f014(iz), v_phot(iz,nb_phot)
!     end do
!     stop

!---  f015: clco + o2 + m -> clco3 + m

!     yung and demore, icarus, 51, 199-247, 1982.

      f015(:) = 5.7E-15*exp(500./t(:))*conc(:)   &
               /(1.e17 + 0.05*conc(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f015(:)

!---  f016: clco3 + cl -> cl + clo + co2

!     yung and demore, icarus, 51, 199-247, 1982.

!     decomposee en :
!     0.5 clco3 + 0.5 cl -> cl + 0.5 co2
!     0.5 clco3 + 0.5 cl -> clo + 0.5 co2

      f016(:) = 1.0E-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f016(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f016(:)

!---  f017: clco3 + o -> cl + o2 + co2

!     yung and demore, icarus, 51, 199-247, 1982.

!     decomposee en :
!     0.5 clco3 + 0.5 o -> cl
!     0.5 clco3 + 0.5 o -> o2 + co2

      f017(:) = 1.0E-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f017(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f017(:)

!---  f018: clo + ho2  -> hocl + o2

!     jpl 2019

      f018(:) = 2.6E-12*exp(290./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f018(:)

!---  f019: oh + hocl -> h2o + clo

      f019(:) = 3.0E-12*exp(-500./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f019(:)

!---  f020: o + hocl -> oh + clo

      f020(:) = 1.7E-13

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f020(:)

!---  f021: cl + cl + co2 -> cl2 + co2

!     donohoue et al., j. phys. chem. a, 109, 7732-7741, 2005

!     f021(:) = 2.5*8.4E-33*exp(850.*(1./t(:) - 1./298.))*conc(:)

!     valeur utilisee par Zhang et al., 2011:

      f021(:) = 2.6E-33*exp(900./t(:))*conc(:)

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = f021(:)

!---  f022: clco + o -> cl + co2

!     yung et al., icarus, 1982 (estimated)

      f022(:) = 3.0E-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f022(:)

!---  f023: cl2 + o(1d) -> cl + clo

!     jpl 2011

      f023(:) = 2.0E-10

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f023(:)

!---  f024: cl2 + h  -> hcl + cl

!     baulch et al., j. phys. chem. ref. data, 1981

      f024(:) = 1.43E-10*exp(-591./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f024(:)

!---  f025: cl + clco  -> cl2 + co

!     baulch et al., j. phys. chem. ref. data, 1981

      f025(:) = 2.16E-9*exp(-1670./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f025(:)

!---  f026: clco + clco  -> cocl2 + co

!     zhang et al., icarus, 2011 (estimated)

      f026(:) = 5.0E-11

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = f026(:)

!---  f027: cl + so2 + co2  -> clso2 + co2

!     mills, phd, 1998

      f027(:) = 1.3E-34*exp(940./t(:))*conc(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f027(:)

!---  f028: clso2 + o  -> cl + so3

!     mills, phd, 1998 (products clo + so2)

!     f028(:) = 1.0E-11

!     croce and cobos, 2018

      f028(:) = 7.69E-11*(t(:)/250.)**0.093

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f028(:)

!---  f029: clso2 + h  -> so2 + hcl

!     mills, phd, 1998

!     f029(:) = 1.0E-11

!     croce and cobos, 2018

      f029(:) = 2.71E-11*(t(:)/250.)**0.47

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f029(:)

!---  f030: clso2 + clso2  -> cl2so2 + so2

!     moses et al. 2002

!     f030(:) = 5.0E-13

!     croce and cobos, 2018

      do iz = 1,nz
         ak1 = 2.6E-14*(t(iz)/250.)**0.61
         ak0 = 3.27E-28*(t(iz)/250.)**(-6.35)

         rate = (ak0*conc(iz))/(1. + ak0*conc(iz)/ak1)
         xpo = 1./(1. + alog10((ak0*conc(iz))/ak1)**2)
         fc = 0.558*exp(-t(iz)/316.) + 0.442*exp(-t(iz)/7442.)
         f030(iz) = rate*fc**xpo
      end do

      f030(:) = 0.  ! temporary FL 15 May 2025

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = f030(:)

!---  f031: cl + o + co2  -> clo + co2

!     yung and demore, 1999 (estimated)

      f031(:) = 5.0E-32*conc(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f031(:)

!---  f032: cl2 + o -> clo + cl

!     mills, phd, 1998

      f032(:) = 7.4E-12*exp(-1650./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f032(:)

!---  f033: clco + oh -> hocl + co

!     mills, phd, 1998

      f033(:) = 1.5E-10

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f033(:)

!---  f034: cl2 + oh -> cl + hocl

!     jpl 2011

      f034(:) = 2.6E-12*exp(-1100./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f034(:)

!---  f035: clco + o -> co + clo

!     yung and demore, 1982

      f035(:) = 3.0E-12

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f035(:)

!---  f036: clco + cl2 -> cocl2 + cl

!     ohta, bull. chem. soc. jpn., 1983

      f036(:) = 6.45E-2*f015(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f036(:)

!---  f037: hcl + h -> h2 + cl

!     mills, phd, 1998

      f037(:) = 1.5E-11*exp(-1750./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f037(:)

!---  f038: clco + h -> hcl + co

!     yung and demore, 1982

      f038(:) = 1.0E-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f038(:)

!---  f039: cl + h + m -> hcl + m

!     yung and demore, 1982 (estimate)

      f039(:) = 1.0E-32*conc(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f039(:)
 
!---     df001: od + oh -> hdo + o

      ! Bedjanian, Y.; Le Bras, G.; and Poulet, G., 
      !Kinetic Study of OH + OH and OD + OD Reactions, 
      !J. Phys. Chem. A, 103, pp. 7017 - 7025, 1999

      df001(:) = 1.5e-12
      
      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df001(:)

!---     df002: od + h2 -> hdo + h

      !Talukdar, R.K. et al., 
      !Kinetics of hydroxyl radical reactions with isotopically labeled hydrogen, 
      !J. Phys. Chem., 100, pp.   3037 - 3043, 1996

      df002(:) = 7.41e-15

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df002(:)

!---     df003: od + ho2 -> hdo + o2
   
      df003(:) = c007(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df003(:)

!---     df004: od + h2o2 -> hdo + ho2

      !Vaghjiani, G.L. & Ravishankara, A.R., 
      !Reactions of OH and OD with H2O2 and D2O2, 
      !J. Phys. Chem., 93, pp. 7833 - 7837, 1989;

      df004(:) = 1.79e-12

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df004(:)

!---     df005: o + od -> o2 + d

      !Following Yung+1998, rate equal to that of O + OH -> O2 + H (c002)

      df005(:) = c002(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df005(:)

!---     df006: od + h2 -> h2o + d

      !Following Yung+1998, rate equal to that of OH + H2 -> H2O + H (c010)

      df006(:) = c010(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df006(:)

!---     df007: od + h -> oh + d

      !Rate following Yung+1988 and the rate of the inverse reaction (df012) from Atahan et al. J. Chem. Phys. 2005

      df007(:) = 1.61e-10*((298./t(:))**0.32)*exp(-16./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df007(:)

!---     df008: co + od -> co2 + d

      !Following Yung+1988 rate equal to that of reaction CO + OH -> CO2 + H

      df008(:) = e001(:) 

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df008(:)

!---     df009: o3 + d -> o2 + od

      !Rate from NIST, Yu & Varandas, J. Chem. Soc. Faraday Trans., 93, 2651-2656, 1997

      df009(:) = 7.41e-11*exp(-379./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df009(:)

!---     df010: HO2 + D -> OH + OD

      !Following Yung+1988, rate equal to 0.71 times the rate for reaction HO2 + H -> OH + OH (c004)

      df010(:) = 0.71*c004(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df010(:)

!---     df011: HO2 + D -> HDO + O

      !Following Yung+1988, rate equal to 0.71 times the rate for reaction HO2 + H -> H2O + O (c006)

      df011(:) = 0.71*c006(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df011(:)

!---     df012: OH + D -> OD + H

      !Rate from NIST, Atahan et al. J. Chem. Phys. 2005

      df012(:) = 1.16e-10*((298./t(:))**0.32)*exp(-16./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df012(:)

!---     df013: h + d + co2 -> hd + co2

      !According to Yung et al. 1988, rate equal to that of H + H + CO2 
      !(reaction c018). Source: baulch et al., 2005

      df013(:) = c018(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df013(:)

!---     df014: HO2 + D -> HD + O2
   
      !According to Yung et al., rate equal to 0.71 times the rate of
      ! H + HO2 -> H2 + O2 (reaction c005, source JPL 2019)

      df014(:) = 0.71*c005(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df014(:)
   
!---     df015: OH + HD -> HDO + H

      !Talukdar et al., Kinetics of hydroxyl radical reactions with 
      !isotopically labeled hydrogen, J. Phys. Chem. 100, 3037-3043, 1996

      df015(:) = 8.5e-13*exp(-2130./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df015(:)

!---     df016: OH + HD -> H2O + D

      !Talukdar et al., 1996

      df016(:) = 4.15e-12*exp(-2130./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df016(:)

!---     df017: D + O2 + CO2 -> DO2 + CO2

      !Breshears et al., Room-temperature rate constants for the reaction 
      !D + O2 + M -> DO2 + M, M=Ar, D2, CO2 and F2. J. Chem. Soc. Faraday 
      !Trans., 87, 2337-2355 (1991)

      df017(:) = 1.59e-31*conc(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df017(:)

!---     df018: OD + O3 -> DO2 + O2

      !According to Yung+1988, rate equal to that of reaccion 
      !OH + O3 -> HO2 + O2, (reaction c014)

      df018(:) = c014(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df018(:)

!---     df019: D + HO2 -> DO2 + H

      !Yung et al., 1988

      df019(:) = 1.0e-10

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df019(:)

!---     df020: O + DO2 -> OD + O2

      !According to Yung+1988, rate equal to that of O + HO2 -> OH + O2
      ! -> reaction c001

      df020(:) = c001(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df020(:)
   
!---     df021: H + DO2 -> OH + OD

      !According to Yung+1988, rate equal to that of H + HO2 -> OH + OH
      ! -> reaction c004

      df021(:) = c004(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df021(:)

!---     df022: H + DO2 -> HD + O2

      !According to Yung+1988, rate equal to that of H + HO2 -> H2 + O2
      ! -> reaction c005

      df022(:) = c005(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df022(:)

!---     df023: H + DO2 -> HDO + O

      !According to Yung+1988, rate equal to that of H + HO2 -> H2O + O
      ! -> reaction c006

      df023(:) = c006(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df023(:)

!---     df024: H + DO2 -> HO2 + D

      !Yung+1988

      df024(:) = 1.85e-10*exp(-890./t(:))
      
      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df024(:)

!---     df025: OH + DO2 -> HDO + O2

      !According to Yung+1988, rate equal to that of OH + HO2 -> H2O + O2
      ! -> reaction c007

      df025(:) = c007(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df025(:)

!---     df026: DO2 + O3 -> OD + O2 + O2

      !According to Yung+1988, rate equal to that of the reaction
      ! HO2 + O3 -> OH + O2 + O2 -> reaction c015

      df026(:) = c015(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df026(:)

!---     df027: OD + OH + CO2 -> HDO2 + CO2

      !According to Yung+1988, rate equal to that of the reaction
      ! OH + OH + CO2 -> H2O2 + CO2 (reaction c017)

      df027(:) = c017(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df027(:)

!---     df028: DO2 + HO2 -> HDO2 + O2

      !According to Yung+1988, rate equal to that of HO2 + HO2 -> H2O2 + O2
      ! (reaction c008)

      df028(:) = c008(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df028(:)

!---     df029: O + HDO2 -> OD + HO2

      !According to Yung+1988, rate half that of O + H2O2 -> OH + HO2
      ! (reaction c012)

      df029(:) = 0.5*c012(:)
   
      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df029(:)

!---     df030: O + HDO2 -> OH + DO2

      !According to Yung+1988, rate half that of O + H2O2 -> OH + HO2
      ! (reaction c012)

      df030(:) = 0.5*c012(:)
   
      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df030(:)

!---     df031: OH + HDO2 -> HDO + HO2

      !According to Yung+1988, rate half that of OH + H2O2 -> H2O + HO2
      ! (reaction c009)

      df031(:) = 0.5*c009(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df031(:)


!---     df032: OH + HDO2 -> H2O + DO2

      !According to Yung+1988, rate half that of OH + H2O2 -> H2O + HO2
      ! (reaction c009)

      df032(:) = 0.5*c009(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df032(:)

!---     df033: OD + H2O2 -> H2O + DO2

      df033(:) = df004(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df033(:)

!---  df034: NO + DO2 -> NO2 + OD
      df034(:) = d003(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df034(:)

!---  df035: no2 + h -> no + oh
      df035(:) = d006(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df035(:)

!---  df036: N + DO2 -> NO + OD
      df036(:) = d008(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df036(:)

!---  d037: N + OD -> NO + D
      df037(:) = d009(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df037(:)

!---  df038: DO2 + HO2 + CO2 -> HDO2 + O2 + CO2
      df038(:) = c016(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df038(:)

!---  df039: dcl + o(1d) -> od + cl
      df039(:) = 1.0E-10

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df039(:)

!---  df040: dcl + o(1d) -> d + clo
      df040(:) = 3.6E-11
      
      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df040(:)

!---  df041: dcl + o -> od + cl
      df041(:) = 1.0E-11*exp(-3300./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df041(:)

!---  df042: hcl + od -> hdo + cl
      df042(:) = 1.8E-12*exp(-250./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df042(:)

!---  df043: dcl + oh -> hdo + cl
      df043(:) = 0.0*1.8E-12*exp(-250./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = 0.0*df043(:)

!---  df044: clo + od -> cl + do2
      df044(:) = 7.4E-12*exp(270./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df044(:)

!---  df045: clo + od -> dcl + o2
      df045(:) = 6.0E-13*exp(230./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df045(:)

!---  df046: cl + hd -> dcl + h
      df046(:) = 3.05E-11*exp(-2270./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df046(:)

!---  df047: cl + hd -> hcl + d
      df047(:) = 3.05E-11*exp(-2270./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df047(:)

!---  df048: cl + do2 -> clo + od
      df048(:) = 3.6E-11*exp(-375./t(:)) 
      
      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df048(:)

!---  df049: cl + do2 -> dcl + o2
      df049(:) = 1.4E-11*exp(270./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df049(:)

!---  df050: cl + hdo2 -> dcl + ho2
      df050(:) = 1.1E-11*exp(-980./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df050(:)

!---  df051: cl + hdo2 -> hcl + do2
      df051(:) = 1.1E-11*exp(-980./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df051(:)

!---  df052: clo + do2  -> docl + o2
      df052(:) = 2.6E-12*exp(290./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df052(:)

!---  df053: oh + docl -> hdo + clo
      df053(:) = 3.0E-12*exp(-500./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df053(:)

!---  df054: od + hocl -> hdo + clo
      df054(:) = 3.0E-12*exp(-500./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df054(:)

!---  df055: o + docl -> od + clo
      df055(:) = 1.7E-13

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df055(:)

!---  df056: cl2 + d  -> dcl + cl
      df056(:) = 1.43E-10*exp(-591./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df056(:)

!---  df057: clso2 + d  -> so2 + dcl
      df057(:) = 1.0E-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df057(:)

!---  df058: clco + od -> docl + co
      df058(:) = 1.5E-10

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df058(:)

!---  df059: cl2 + od -> cl + docl
      df059(:) = 2.6E-12*exp(-1100./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df059(:)

!---  df060: dcl + h -> hd + cl
      df060(:) = 1.5E-11*exp(-1750./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df060(:)

!---  df061: hcl + d -> hd + cl
      df061(:) = 1.5E-11*exp(-1750./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df061(:)

!---  df062: clco + d -> dcl + co
      df062(:) = 1.0E-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df062(:)

!---  df063: cl + d + m -> dcl + m
      df063(:) = 1.0E-32*conc(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = df063(:)

!---  f040: clso2 + cl -> cl2so2

!     croce and cobos, 2018

      f040(:) = 1.44e-11*(t(:)/250.)**0.47

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = f040(:)

!----------------------------------------------------------------------
!        reactions des composes soufres
!----------------------------------------------------------------------

!---  g001: s + o2 -> so + o

!     jpl 2015

      g001(:) = 1.6E-12*exp(100./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g001(:)

!---  g002: s + o3 -> so + o2

!     jpl 2015

      g002(:) = 1.2E-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g002(:)

!---  g003: so + o2 -> so2 + o

!     jpl 2015

      g003(:) = 1.6E-13*exp(-2280./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g003(:)

!---  g004: so + o3 -> so2 + o2

!     jpl 2015

      g004(:) = 3.4E-12*exp(-1100./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g004(:)

!---  g005: so + oh -> so2 + h

!     jpl 2019

      g005(:) = 2.6E-11*exp(330./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g005(:)

!---  g006: s + oh -> so + h

!     jpl 2015

      g006(:) = 6.6E-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g006(:)

!---  g007: so + o + co2 -> so2 + co2

!     singleton and cvetanovic, j. phys. chem. ref. data, 1988
!     measured with co2 as third body

      do iz = 1,nz
         ak0 = 4.2E-30
         ak1 = 5.3E-11

         rate = (ak0*conc(iz))/(1. + ak0*conc(iz)/ak1)
         xpo = 1./(1. + alog10((ak0*conc(iz))/ak1)**2)
         g007(iz) = rate*0.6**xpo
      end do

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g007(:)

!---  g008: so + ho2 -> so2 + oh 

      g008(:) = 2.8E-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g008(:)

!---  g009: so2 + o + co2 -> so3 + co2

!     Naidoo et al., Proceedings of the Combustion Institute,
!                    30, 1219-1225, 2005
!     Also recommended by jpl 2019 with a simpler expression (fc = 0.6)
!
!     Factor of 5 for third-body efficiency of co2 to be confirmed!

      do iz = 1,nz
         ak0 = 5.*9.5e-23*(t(iz)**(-3.0))*exp(-2400./t(iz))
         ak1 = 6.1e-13*exp(-850./t(iz))
         rate = (ak0*conc(iz))/(1. + ak0*conc(iz)/ak1)
         xpo = 1./(1. + alog10((ak0*conc(iz))/ak1)**2)
         fc = 0.558*exp(-t(iz)/316.) + 0.442*exp(-t(iz)/7442.)
         g009(iz) = rate*fc**xpo
      end do

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g009(:)

!---  g010: s + o + co2 -> so + co2

!     zhang et al., icarus, 2011

      g010(:) = 1.5E-34*exp(900./t(:))*conc(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g010(:)

!---  g011: so3 + h2o + h2o -> h2so4 + h2o

!     lovejoy et al., j.phys.chem., 1996

      do iz = 1,nz
         g011(iz) = 2.26E-23*max(t(iz),100.)*exp(6540./max(t(iz),100.))
         g011(iz) = g011(iz)*1.0E-20*c(iz,i_h2o)
      end do

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g011(:)

!---  g012: so + clo -> so2 + cl 

!     jpl 2011

      g012(:) = 2.8E-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g012(:)

!---  g013: so + so3 -> so2 + so2 

!     chung et al., int. j. chem. kinet., 1975

      g013(:) = 2.0E-15

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g013(:)

!---  g014: so3 + o -> so2 + o2 

!     jacob and winkler, j. chem. soc. faraday trans. 1, 1972

      g014(:) = 2.32E-16*exp(-487./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g014(:)

!---  g015: free slot

!---  g016: free slot

!---  g017: clco3 + so -> cl + so2 + co2

!     mills, phd, 1998

!     decomposee en :
!     0.5 clco3 + 0.5 so -> cl
!     0.5 clco3 + 0.5 so -> so2 + co2

      g017(:) = 1.0E-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g017(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g017(:)

!---  g018: s + co + co2 -> ocs + co2

!     zhang et al., icarus, 2011 (estimate?)

      g018(:) = 2.5*4.0E-33*exp(-1940./t(:))*conc(:)
      
!      g018(:) = 0.0E+0

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g018(:)

!---  g019: clco + s -> ocs + cl 

!     zhang et al., icarus, 2011

      g019(:) = 3.0E-12

!      g019(:) = 0.0E+0

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g019(:)

!---  g020: so2 + oh + co2 -> hso3 + co2

!     jpl 2011

      do iz = 1,nz
         ak0 = 2.5*3.3E-31*(t(iz)/300.)**(-4.3)
         ak1 = 1.6E-12

         rate = (ak0*conc(iz))/(1. + ak0*conc(iz)/ak1)
         xpo = 1./(1. + alog10((ak0*conc(iz))/ak1)**2)
         g020(iz) = rate*0.6**xpo
      end do

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g020(:)

!---  g021: hso3 + o2 -> ho2 + so3 

!     jpl 2011

      g021(:) = 1.3E-12*exp(-330./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g021(:)

!---  g022: s + s + co2 -> s2 + co2

!     nicholas et al., j. chem. soc. faraday trans. 1, 1979

      do iz = 1,nz
         ak0 = 1.19E-29
         ak1 = 1.0E-10

         rate = (ak0*conc(iz))/(1. + ak0*conc(iz)/ak1)
         xpo = 1./(1. + alog10((ak0*conc(iz))/ak1)**2)
         g022(iz) = rate*0.6**xpo
      end do

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = g022(:)

!---  g023: s2 + o -> so + s 

!     zhang et al., icarus, 2011

      g023(:) = 2.2E-11*exp(-84./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g023(:)

!---  g024: s + ocs -> s2 +  co

!     lu et al., j. chem. phys., 2006

      g024(:) = 6.63E-20*(t(:)**2.57)*exp(-1180./t(:))
      
      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g024(:)

!---  g025: ocs + o -> so + co

!     atkinson et al., 2004

      g025(:) = 1.60E-11*exp(-2150./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g025(:)

!---  g026: s + so3 -> so2 +  so

!     moses et al., 2002

      g026(:) = 1.0E-16

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g026(:)

!---  g027: s + ho2 -> so +  oh

!     yung and demore, 1982

      g027(:) = 3.0E-11*exp(200./t(:))

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g027(:)

!---  g028: s + clo -> so +  cl

!     moses et al., 2002

      g028(:) = 4.0E-11

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g028(:)

!---  g029: h2so4 + h2o -> so3 + h2o + h2o 

!     krasnopolsky , 2007 

      g029(:) = 7.0E-14*exp(-5170./t(:))

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = g029(:)*c(:,i_h2o)
      
!---  g030: free slot

!---  g031: free slot

!---  g032: so + so -> so2 + s

!     krasnopolsky, 2012

      g032(:) = 1.0e-12*exp(-1700./t(:))
      
      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = g032(:)
 
!      dg001: so3 + hdo -> hdso4
!     lovejoy et al., j.phys.chem., 1996

      do iz = 1,nz
         dg001(iz) = 2.26E-23*max(t(iz),100.)*exp(6540./max(t(iz),100.))
         dg001(iz) = dg001(iz)*1.0E-20*c(iz,i_h2o)
      end do
   
      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = dg001(:)

!      dg002: hdso4 + h2o -> so3 + hdo + h2o
      dg002(:) = 7.0E-14*exp(-5170./t(:))
      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = dg002(:)*c(:,i_h2o)

!      dg003: hdso4 + hdo -> so3 + hdo + hdo
      dg003(:) = 7.0E-14*exp(-5170./t(:))
      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = dg003(:)*c(:,i_hdo)*0.0

!---  dg004: SO + OD -> SO2 + D
      dg004(:) = g005(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = dg004(:)

!---  dg005: S + OD -> SO + D
      dg005(:) = g006(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = dg005(:)

!---  dg006: SO + DO2 -> SO2 + OD
      dg006(:) = g008(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = dg006(:)

!---  dg007: S + DO2 -> SO +  OD
      dg007(:) = g027(:)

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = dg007(:)


!---  g033: so + so + co2 -> osso_cis + co2
                       
!     Egan et al., 2025  
                       
      do iz = 1,nz   
         ak0 = 2.5*1.54e-31*(t(iz)/298.)**(-3.36)
         ak1 = 1.1e-10*(t(iz)/298.)**0.167
         rate = (ak0*conc(iz))/(1. + ak0*conc(iz)/ak1)
         xpo = 1./(1. + alog10((ak0*conc(iz))/ak1)**2)
         g033(iz) = rate*0.42**xpo
      end do

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = g033(:)
            
!---  g034: osso_cis + co2 -> so + so + co2

!     Egan et al., 2025

      deq(:) = 1.02e-27*exp(17231./t(:))
      g034(:) = g033(:)/(deq(:)*conc(:))

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = g034(:)*conc(:)

!---  g035: so + so + co2 -> osso_trans + co2
                       
!     Egan et al., 2025  
                       
      do iz = 1,nz   
         ak0 = 2.5*1.13e-31*(t(iz)/298.)**(-3.38)
         ak1 = 1.1e-10*(t(iz)/298.)**0.167
         rate = (ak0*conc(iz))/(1. + ak0*conc(iz)/ak1)
         xpo = 1./(1. + alog10((ak0*conc(iz))/ak1)**2)
         g035(iz) = rate*0.42**xpo
      end do

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = g035(:)

!---  g036: osso_trans + co2 -> so + so + co2

!     Egan et al., 2025

      deq(:) = 1.73e-27*exp(15395./t(:))
      g036(:) = g035(:)/(deq(:)*conc(:))

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = g036(:)*conc(:)
            
!---  g037: so + so + co2 -> s2o2_cyc + co2
                       
!     Egan et al., 2025  
                       
      do iz = 1,nz   
         ak0 = 2.5*1.70e-32*(t(iz)/298.)**(-3.38)
         ak1 = 1.1e-10*(t(iz)/298.)**0.167
         rate = (ak0*conc(iz))/(1. + ak0*conc(iz)/ak1)
         xpo = 1./(1. + alog10((ak0*conc(iz))/ak1)**2)
         g037(iz) = rate*0.42**xpo
      end do

      nb_reaction_3 = nb_reaction_3 + 1
      v_3(:,nb_reaction_3) = g037(:)

!---  g038: s2o2_cyc + co2 -> so + so + co2

!     Egan et al., 2025

      deq(:) = 6.75e-28*exp(13392./t(:))
      g038(:) = g037(:)/(deq(:)*conc(:))

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = g038(:)*conc(:)

!---  g040: osso_cis + o -> so2 + so

!     Egan et al., 2025

      g040(:) = 1.10e-10*(t(:)/298.)**0.167

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g040(:)

!---  g042: osso_trans + o -> so2 + so

!     Egan et al., 2025

      g042(:) = 1.10e-10*(t(:)/298.)**0.167

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g042(:)

!---  g044: s2o2_cyc + o -> so2 + so

!     Egan et al., 2025

      g044(:) = 1.10e-10*(t(:)/298.)**0.167

      nb_reaction_4 = nb_reaction_4 + 1
      v_4(:,nb_reaction_4) = g044(:)

!----------------------------------------------------------------------
!     heterogeneous chemistry 
!----------------------------------------------------------------------

      if (hetero_ice) then

!        k = (surface*v*gamma)/4 (s-1)
!        v = 100*sqrt(8rt/(pi*m))  (cm s-1)
 
!---     h001: ho2 + ice -> products
 
!        cooper and abbatt, 1996: gamma = 0.025
      
         gam = 0.025
         h001(:) = surfice1d(:)*1.E-8       &
                   *100.*sqrt(8.*8.31*t(:)/(33.E-3*pi))*gam/4.
 
!        h002: oh + ice -> products
 
!        cooper and abbatt, 1996: gamma = 0.03
 
         gam = 0.03
         h002(:) = surfice1d(:)*1.E-8       &
                   *100.*sqrt(8.*8.31*t(:)/(17.E-3*pi))*gam/4.

!---     h003: h2o2 + ice -> products
 
!        gamma = 0.    test value
 
         gam = 0.
         h003(:) = surfice1d(:)*1.E-8        &
                   *100.*sqrt(8.*8.31*t(:)/(34.E-3*pi))*gam/4.
      else
         h001(:) = 0.
         h002(:) = 0.
         h003(:) = 0.
      end if

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = h001(:)

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = h002(:)

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = h003(:)

!     do iz = 1,nz
!        print*, z(iz), surfice1d(iz), h001(iz), h002(iz)
!     end do
!     stop

!     print*, 'krates : nb_phot       = ', nb_phot
!     print*, 'krates : nb_reaction_4 = ', nb_reaction_4
!     print*, 'krates : nb_reaction_3 = ', nb_reaction_3
!     stop

!----------------------------------------------------------------------
!     ionospheric reactions
!     only if ok_ionchem=true
!----------------------------------------------------------------------

      if (ok_ionchem) then

!---     i001: co2+ + o2 -> o2+ + co2

!        aninich, j. phys. chem. ref. data 1993

         i001(:) = 5.5e-11*(300./t_elect(:))**0.5

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i001(:)

!---     i002: co2+ + o -> o+ + co2

!        UMIST database

         i002(:) = 9.6e-11
      
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i002(:)

!---     i003: co2+ + o -> o2+ + co

!        UMIST database

         i003(:) = 1.64e-10

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i003(:)

!---     i004: o2+ + e- -> o + o

!        Alge et al., J. Phys. B, At. Mol. Phys. 1983
         !i004(:) = 2.0e-7*(300./t_elect(:))**0.7
          
         do iz = 1,nz
           if (t_elect(iz)<1200.) then
             i004(iz) = 2.0e-7*(300./t_elect(iz))**0.7
           else
             i004(iz) = 7.4e-8*(1200./t_elect(iz))**0.56
           end if
         end do

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i004(:)

!---     i005: o+ + co2 -> o2+ + co

!        UMIST database

         i005(:) = 9.4e-10

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i005(:)


!---     i006: co2+ + e- -> co + o

!        UMIST database

         i006(:) = 3.8e-7*(300./t_elect(:))**0.5

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i006(:)


!---     i007: co2+ + no -> no+ + co2

!        UMIST database

         i007(:) = 1.2e-10

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i007(:)

!---     i008: o2+ + no -> no+ + o2

!        UMIST database

         i008(:) = 4.6e-10

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i008(:)

!---     i009: o2+ + n2 -> no+ + no
      
!        Fox & Sung 2001

         i009(:) = 1.0e-15
      
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i009(:)

!---     i010: o2+ + n -> no+ + o

!        Fox & Sung 2001

         i010(:) = 1.0e-10

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i010(:)

!---     i011: o+ + n2 -> no+ + n

!        Fox & Sung 2001

         i011(:) = 1.2e-12 * (300./t_elect(:))**0.45

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i011(:)

!---     i012: no+ + e -> n + o

!        UMIST database

         i012(:) = 4.3e-7*(300./t_elect(:))**0.37

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i012(:)


!---     i013: co+ + co2 -> co2+ + co

!        UMIST database

         i013(:) = 1.0e-9

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i013(:)


!---     i014: co+ + o -> o+ + co

!        UMIST database

         i014(:) = 1.4e-10

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i014(:)

!---     i015: c+ + co2 -> co+ + co

!        UMIST database

         i015(:) = 1.1e-9

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i015(:)


!---     i016: N2+ + co2 -> co2+ + N2

!        Fox & Song 2001

         i016(:) = 9.0e-10*(300./t_elect(:))**0.23

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i016(:)


!---     i017: N2+ + o -> no+ + N

!        Fox & Song 2001

         i017(:) = 1.33e-10*(300./t_elect(:))**0.44

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i017(:)

!---     i018: N2+ + co -> co+ + N2

!        UMIST

         i018(:) = 7.4e-11

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i018(:)

!---     i019: N2+ + e -> N + N

!        UMIST

         i019(:) = 1.7e-7*(300./t_elect(:))**0.3

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i019(:)

!---     i020: N2+ + o -> o+ + N2

!        Fox & Song 2001

         i020(:) = 7.0e-12*(300./t_elect(:))**0.23

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i020(:)

!---     i021: N+ + co2 -> co2+ + N

!        UMIST

         i021(:) = 7.5e-10

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i021(:)

!---     i022: CO+ + H -> H+ + CO

!        Fox & Sung 2001

         i022(:) = 4.0e-10

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i022(:)

!---     i023: O+ + H -> H+ + O

!        UMIST

         i023(:) = 5.66e-10*((t_elect(:)/300.)**0.36)*exp(8.6/t_elect(:))

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i023(:)

!---     i024: H+ + O -> O+ + H

!        UMIST

         i024(:) = 6.86e-10*((t_elect(:)/300.)**0.26)*exp(-224.3/t_elect(:))

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i024(:)

!---     i025: CO+ + H2 -> HCO2+ + H

!        UMIST

         i025(:) = 9.5e-10

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i025(:)

!---     i026: spare slot

         i026(:) = 0.

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i026(:)

!---     i027+i028: HCO2+ + e -> H + O + CO

!        UMIST
         !Reaction splitted in 2: i027: 0.5 (HCO2+ + e-) -> H
         !i028: 0.5 (HCO2+ + e-) -> O + CO

         i027(:) = 8.1e-7*((300./t_elect(:))**0.64)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i027(:)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i027(:)

!---     i029: HCO2+ + e -> OH + CO

!        UMIST

         i029(:) = 3.2e-7*((300./t_elect(:))**0.64)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i029(:)

!---     i030: HCO2+ + e -> H + CO2

!        UMIST

         i030(:) = 6.0e-8*((300./t_elect(:))**0.64)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i030(:)

!---     i031: HCO2+ + O -> HCO+ + O2

!        UMIST

         i031(:) = 1.e-9
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i031(:)

!---     i032: HCO2+ + CO -> HCO+ + CO2

!        UMIST, from Prassad & Huntress 1980

         i032(:) = 7.8e-10
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i032(:)

!---     i033: H+ + CO2 -> HCO+ + O

!        UMIST, from Smith et al., Int. J. Mass Spectrom. Ion Proc., 117, 457-473(1992) 

         i033(:) = 3.5e-9
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i033(:)


!---     i034: CO2+ + H -> HCO+ + O

!        Seen in Fox 2015, from Borodi et al., Int. J. Mass Spectrom. 280, 218-225, 2009

         i034(:) = 4.5e-10
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i034(:)

!---     i035: CO+ + H2 -> HCO+ + H

         !UMIST, from Scott et al., J. Chem. Phys., 106, 3982-3987(1997)

         i035(:) = 7.5e-10
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i035(:)

!---     i036: HCO+ + e- -> CO + H

         !UMIST, from Mitchell, Phys. Rep., 186, 215 (1990)

         i036(:) = 2.4e-7 *((300./t_elect(:))**0.69)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i036(:)

!---     i037: CO2+ + H2O -> H2O+ + CO2

         !UMIST, from Karpas, Z., Anicich, V.G., and Huntress, W.T., Chem. Phys. Lett., 59, 84 (1978)

         i037(:) = 2.04e-9 *((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i037(:)

!---     i038: CO+ + H2O -> H2O+ + CO

         !UMIST, from Huntress, W.T., McEwan, M.J., Karpas, Z., and Anicich, V.G., Astrophys. J. Supp. Series, 44, 481 (1980)

         i038(:) = 1.72e-9*((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i038(:)

!---     i039: O+ + H2O -> H2O+ + O

         !UMIST, from Adams, N.G., Smith, D., and Paulson, J.F., J. Chem. Phys., 72, 288 (1980); Smith, D., Adams, N.G., and Miller, T.M., J. Chem. Phys.., 69, 308 (1978)

         i039(:) = 3.2e-9*((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i039(:)

!---     i040: N2+ + H2O -> H2O+ + N2

         !UMIST, from Adams, N.G., Smith, D., and Paulson, J.F., J. Chem. Phys., 72, 288 (1980); Smith, D., Adams, N.G., and Miller, T.M., J. Chem. Phys.., 69, 308 (1978)

         i040(:) = 2.3e-9*((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i040(:)

!---     i041: N+ + H2O -> H2O+ + N

         !UMIST, from Adams, N.G., Smith, D., and Paulson, J.F., J. Chem. Phys., 72, 288 (1980); Smith, D., Adams, N.G., and Miller, T.M., J. Chem. Phys.., 69, 308 (1978)

         i041(:) = 2.8e-9*((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i041(:)


!---     i042: H+ + H2O -> H2O+ + H

         !UMIST, from D. Smith, P. Spanel and C. A. Mayhew, Int. J. Mass Spectrom. Ion Proc., 117, 457-473(1992)

         i042(:) = 6.9e-9*((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i042(:)

!---     i043: H2O+ + O2 -> O2+ + H2O

         !UMIST, from A. B. Raksit and P. Warneck, J. Chem. Soc. Faraday Trans., 76, 1084-1092(1980)

         i043(:) = 4.6e-10
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i043(:)

!---     i044: H2O+ + CO -> HCO+ + OH

         !UMIST, from Jones, J.D.C., Birkinshaw, K., and Twiddy, N.D., Chem. Phys. Lett., 77, 484 (1981)

         i044(:) = 5.0e-10
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i044(:)

!---     i045: H2O+ + O -> O2+ + H2

         !UMIST, from Viggiano, A.A, Howarka, F., Albritton, D.L., Fehsenfeld, F.C., Adams, N.G., and Smith, D., Astrophys. J., 236, 492 (1980)

         i045(:) = 4.0e-11
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i045(:)

!---     i046: H2O+ + NO -> NO+ + H2O

         !UMIST, from A. B. Raksit and P. Warneck, J. Chem. Soc. Faraday Trans., 76, 1084-1092(1980)

         i046(:) = 2.7e-10
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i046(:)

!---     i047: H2O+ + e- -> H + H + O

         !UMIST, from Rosen, S., Derkatch, A., Semaniak, J., et al., 2000, Far. Disc., 115, 295
         
         i047(:) = 3.05e-7*((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i047(:)

!---     i048: H2O+ + e- -> H + OH
         
         !UMIST, from Rosen, S., Derkatch, A., Semaniak, J., et al., 2000, Far. Disc., 115, 295

         i048(:) = 8.6e-8*((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i048(:)

!---     i049: H2O+ + e- -> O + H2

         !UMIST, from Rosen, S., Derkatch, A., Semaniak, J., et al., 2000, Far. Disc., 115, 295

         i049(:) = 3.9e-8*((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i049(:)

!---     i050: H2O+ + H2O -> H3O+ + OH

         !UMIST, from Huntress, W.T. and Pinizzotto, R.F., J. Chem. Phys., 59, 4742 (1973)

         i050(:) = 2.1e-9*((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i050(:)


!---     i051: H2O+ + H2 -> H3O+ + H

         !UMIST, from A. B. Raksit and P. Warneck, J. Chem. Soc. Faraday Trans., 76, 1084-1092(1980)

         i051(:) = 6.4e-10
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i051(:)

!---     i052: HCO+ + H2O -> H3O+ + CO

         !UMIST, from Adams, N.G., Smith, D., and Grief, D., Int. J. Mass Spectrom. Ion Phys., 26, 405 (1978)

         i052(:) = 2.5e-9*((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i052(:)

!---     i053: H3O+ + e -> OH + H + H

         !UMIST, from Novotny, O. et al., J. Phys. Chem. A 2010, 114, 14, 4870-4874

         i053(:) = 3.05e-7*((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i053(:)

!---     i054: H3O+ + e -> H2O + H

         !UMIST, from Novotny, O. et al., J. Phys. Chem. A 2010, 114, 14, 4870-4874
         
         i054(:) = 7.09e-8*((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i054(:)

!---     i055: H3O+ + e -> OH + H2

         !UMIST, from Novotny, O. et al., J. Phys. Chem. A 2010, 114, 14, 4870-4874

         i055(:) = 5.37e-8*((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i055(:)

!---     i056: H3O+ + e -> O + H2 + H

         !UMIST, from Novotny, O. et al., J. Phys. Chem. A 2010, 114, 14, 4870-4874

         i056(:) = 5.6e-9*((300./t_elect(:))**0.5)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i056(:)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i056(:)

!---     i057: O+ + H2 -> OH+ + H 

         !UMIST, from Adams, N.G., Smith, D., and Paulson, J.F., J. Chem. Phys., 72, 288 (1980); Smith, D., Adams, N.G., and Miller, T.M., J. Chem. Phys.., 69, 308 (1978)

         i057(:) = 1.7e-9
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i057(:)

!---     i058: OH+ + O -> O2+ + H

         !UMIST, from Prasad & Huntress, 1980, ApJS, 43, 1

         i058(:) = 7.1e-10
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i058(:)

!---     i059: OH+ + CO2 -> HCO2+ + O

         !UMIST, from Jones, J.D.C., Birkinshaw, K., and Twiddy, N.D., Chem. Phys. Lett., 77, 484 (1981)

         i059(:) = 1.44e-9
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i059(:)

!---     i060: OH+ + CO -> HCO+ + O

         !UMIST, from Jones, J.D.C., Birkinshaw, K., and Twiddy, N.D., Chem. Phys. Lett., 77, 484 (1981)

         i060(:) = 1.05e-9
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i060(:)

!---     i061: OH+ + NO -> NO+ + OH (tasa de reacciÃ³n UMIST 3.59e-10)

         !UMIST, from Jones, J.D.C., Birkinshaw, K., and Twiddy, N.D., Chem. Phys. Lett., 77, 484 (1981)

         i061(:) = 3.59e-10
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i061(:)

!---     i062: OH+ + H2 -> H2O+ + H (tasa de reacciÃ³n UMIST 1.01e-9,

         !UMIST, from Jones, J.D.C., Birkinshaw, K., and Twiddy, N.D., Chem. Phys. Lett., 77, 484 (1981)

         i062(:) = 1.01e-9
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i062(:)

!---     i063: OH+ + O2 -> O2+ + OH (tasa de reacciÃ³n UMIST 5.9e-10

         !UMIST, from Jones, J.D.C., Birkinshaw, K., and Twiddy, N.D., Chem. Phys. Lett., 77, 484 (1981)

         i063(:) = 5.9e-10
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = i063(:)

      end if   ! ok_ionchem
 
      if(ok_ionchem) then

!---     di001: CO2+ + HD -> HCO2+ + D 

         !Rate for i025: CO2+ + H2 -> HCO2+ + H, multiplied by 0.82
         !(Langevin formula, see Krasnopolsky+2002)
         !and divided by two (reaction splitted in two channels HCO2+ + D, di001,
         !and DCO2+ + H, di002

         di001(:) = i025(:) * 0.82 / 2.
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di001(:)

!---     di002: CO2+ + HD -> HCO2+ + D 

         !Rate for i025: CO2+ + H2 -> HCO2+ + H, multiplied by 0.82
         !(Langevin formula, see Krasnopolsky+2002)
         !and divided by two (reaction splitted in two channels HCO2+ + D, di001,
         !and DCO2+ + H, di002

         di002(:) = i025(:) * 0.82 / 2.
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di002(:)

!---     di003: CO+ + HD -> HCO+ + D

         !Rate for i035: CO+ + H2 -> HCO+ + H, multiplied by 0.82
         !(Langevin formula, see Krasnopolsky+2002)
         !and divided by two (reaction splitted in two channels HCO+ + D, di003,
         !and DCO+ + H, di004

         di003(:) = i035(:) * 0.82 / 2.
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di003(:)

!---     di004: CO+ + HD -> DCO+ + H

         !Rate for i035: CO+ + H2 -> HCO+ + H, multiplied by 0.82
         !(Langevin formula, see Krasnopolsky+2002)
         !and divided by two (reaction splitted in two channels HCO+ + D, di003,
         !and DCO+ + H, di004

         di004(:) = i035(:) * 0.82 / 2.
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di004(:)

!---     di005: H2O+ + HD -> H3O+ + D

         !Rate for i051: H2O+ + H2 -> H3O+ + H, multiplied by 0.82
         !(Langevin formula, see Krasnopolsky+2002)

         di005(:) = i051(:) * 0.82
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di005(:)

!---     di007: O+ + HD -> OH+ + D

         !Rate for i057: O+ + H2 -> OH+ + H, multiplied by 0.82
         !(Langevin formula, see Krasnopolsky+2002)

         di007(:) = i057(:) * 0.82
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di007(:)

!---     di008: O+ + HD -> OD+ + H

         !Rate for i057: O+ + H2 -> OH+ + H, multiplied by 0.18?
         di008(:) = i057(:) * 0.18
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di008(:)

!---     di009: OH+ + HD -> H2O+ + D

         !Rate for i062: OH+ + H2 -> H2O+ + H, multiplied by 0.82
         !(Langevin formula, see Krasnopolsky+2002)

         di009(:) = i062(:) * 0.82
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di009(:)

!---     di010: OH+ + HD -> HDO+ + H

         !Rate for i062: OH+ + H2 -> H2O+ + H, multiplied by 0.18?
         !(Langevin formula, see Krasnopolsky+2002)

         di010(:) = i062(:) * 0.18
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di010(:)

!---     di011: DCO2+ + e -> D + CO2

         !Same rate as for i030: HCO2+ + e -> H + CO2

         di011(:) = i030(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di011(:)

!---     di012: DCO2+ + e -> D + O + CO

         !Reaction splitted in 2: di012: 0.5 (DCO2+ + e-) -> D
         !di012: 0.5 (DCO2+ + e-) -> O + CO
         !Same rate as i027 (HCO2+ + e- -> H + O + CO)

         di012(:) = i027(:)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di012(:)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di012(:)         

!---     di013: DCO2+ + e -> OD + CO

         !Same rate as for i029: HCO2+ + e -> OH + CO

         di013(:) = i029(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di013(:)

!---     di014: DCO+ + e -> D + CO

         !Same rate as for i036: HCO+ + e -> H + CO

         di014(:) = i036(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di014(:)

!---     di015: DCO2+ + O -> DCO+ + O2

         !Same rate as for i031: HCO2+ + O -> HCO+ + O2

         di015(:) = i031(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di015(:)

!---     di016: DCO2+ + CO -> DCO+ + CO2

         !Same rate as for i032: HCO2+ + CO -> HCO+ + CO2

         di016(:) = i032(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di016(:)

!---     di017: CO2+ + D -> DCO+ + O

         !Same rate as for i034: CO2+ + H -> HCO+ + O

         di017(:) = i034(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di017(:)

!---     di018: CO2+ + HDO -> HDO+ + CO2

         !Same rate as for i037: CO2+ + H2O -> H2O+ + CO2

         di018(:) = i037(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di018(:)

!---     di019: CO+ + HDO -> HDO+ + CO

         !Same rate as for i038: CO+ + H2O -> H2O+ + CO

         di019(:) = i038(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di019(:)

!---     di020: O+ + HDO -> HDO+ + O

         !Same rate as for i039: O+ + H2O -> H2O+ + O

         di020(:) = i039(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di020(:)

!---     di021: N2+ + HDO -> HDO+ + N2

         !Same rate as for i040: N2+ + H2O -> H2O+ + N2

         di021(:) = i040(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di021(:)

!---     di022: N+ + HDO -> HDO+ + N

         !Same rate as for i041: N+ + H2O -> H2O+ + N

         di022(:) = i041(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di022(:)

!---     di023: H+ + HDO -> HDO+ + H

         !Same rate as for i042: H+ + H2O -> H2O+ + H

         di023(:) = i042(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di023(:)

!---     di024: HDO+ + O2 -> O2+ + HDO

         !Same rate as for i043: H2O+ + O2 -> O2+ + H2O

         di024(:) = i043(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di024(:)

!---     di025: HDO+ + CO -> DCO+ + OH

         !Same rate as for i044: H2O+ + CO -> HCO+ + OH
         !divided by two (to take into account other channel,
         !HDO+ + CO -> HCO+ + OD, di026)

         di025(:) = i044(:)/2.
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di025(:)

!---     di026: HDO+ + CO -> HCO+ + OD

         !Same rate as for i044: H2O+ + CO -> HCO+ + OH
         !divided by two (to take into account other channel,
         !HDO+ + CO -> DCO+ + OH, di025)

         di026(:) = i044(:)/2.
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di026(:)

!---     di027: HDO+ + O -> O2+ + HD

         !Same rate as for i045: H2O+ + O -> O2+ + H2

         di027(:) = i045(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di027(:)

!---     di028: HDO+ + NO -> NO+ + HDO

         !Same rate as for i046: H2O+ + NO -> NO+ + H2O

         di028(:) = i046(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di028(:)

!---     di029: HDO+ + e -> H + D + O

         !Reaction splitted in 2: 0.5 (HDO+ + e-) -> D
         !0.5 (HDO+ + e-) -> H + O
         !Same rate as i047 (H2O+ + e- -> 2H + O)

         di029(:) = i047(:)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di029(:)

         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di029(:)

!---     di030: HDO+ + e- -> D + OH

         !Same rate as for i048: H2O+ + e- -> H + OH
         !divided by 2 (to take into account other channel,
         !HDO+ + e- -> H + OD, di031)

         di030(:) = i048(:)/2.
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di030(:)

!---     di031: HDO+ + e- -> H + OD

         !Same rate as for i048: H2O+ + e- -> H + OH
         !divided by 2 (to take into account other channel,
         !HDO+ + e- -> D + OH, di030)

         di031(:) = i048(:)/2.
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di031(:)

!---     di032: HDO+ + e- -> HD + O

         !Same rate as for i049: H2O+ + e- -> H2 + O

         di032(:) = i049(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di032(:)

!---     di033: D+ + O -> O+ + D
         !Same rate as for i024: H+ + O -> O+ + H
         di033(:) = i024(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di033(:)

!---     di034: D+ + CO2 -> DCO+ + O
         !Same rate as for i033: H+ + CO2 -> HCO+ + O
         di034(:) = i033(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di034(:)

!---     di035: D+ + H2O -> HDO+ + H
         !Same rate as for i042: H+ + H2O -> H2O+ + H
         di035(:) = i042(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di035(:)

!---     di036: OD+ + O -> O2+ + D
         !Same rate as for i058: OH+ + O -> O2+ + H
         di036(:) = i058(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di036(:)

!---     di037: OD+ + CO2 -> DCO2+ + O
         !Same rate as for i059: OH+ + CO2 -> HCO2+ + O
         di037(:) = i059(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di037(:)

!---     di038: OD+ + CO -> DCO+ + O
         !Same rate as for i060: OH+ + CO -> HCO+ + O
         di038(:) = i060(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di038(:)

!---     di039: OD+ + NO -> NO+ + OD 
         !Same rate as for i061: OH+ + NO -> NO+ + OH 
         di039(:) = i061(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di039(:)

!---     di040: OD+ + H2 -> HDO+ + H 
         !Same rate as for i062: OH+ + H2 -> H2O+ + H 
         di040(:) = i062(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di040(:)

!---     di041: OD+ + O2 -> O2+ + OD
         !Same rate as for i063: OH+ + O2 -> O2+ + OH
         di041(:) = i063(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di041(:)

!---     di042: CO+ + D -> D+ + CO 
         !Same rate as for i022: CO+ + H -> H+ + CO 
         di042(:) = i022(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di042(:)

!---     di043: O+ + D -> D+ + O 
         !Same rate as for i023: O+ + H -> H+ + O 
         di043(:) = i023(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di043(:)

!---     di044: H2O+ + HD -> H2DO+ + H
         !Same rate as for i051: H2O+ + H2 -> H3O+ + H 
         di044(:) = i051(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di044(:)

!---     di045: HDO+ + H2O -> H2DO+ + OH
         !Same rate as for i050: H2O+ + H2O -> H3O+ + OH 
         di045(:) = i050(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di045(:)

!---     di046: H2DO+ + e- -> HDO + H
         !Branch of i054: H3O+ + e- -> H2O + H (2/3 stat weight)
         di046(:) = i054(:) * 0.6667
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di046(:)

!---     di047: H2DO+ + e- -> H2O + D
         !Branch of i054: H3O+ + e- -> H2O + H (1/3 stat weight)
         di047(:) = i054(:) * 0.3333
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di047(:)

!---     di048: H2DO+ + e- -> OH + HD
         !Branch of i055: H3O+ + e- -> OH + H2 (2/3 stat weight)
         di048(:) = i055(:) * 0.6667
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di048(:)

!---     di049: H2DO+ + e- -> OD + H2
         !Branch of i055: H3O+ + e- -> OH + H2 (1/3 stat weight)
         di049(:) = i055(:) * 0.3333
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di049(:)

!---     di050: H2DO+ + e- -> OD + H + H
         !Branch of i053: H3O+ + e- -> OH + H + H (1/3 stat weight)
         di050(:) = i053(:) * 0.3333
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di050(:)

!---     di051: H2DO+ + e- -> OH + D + H
         !Branch of i053: H3O+ + e- -> OH + H + H (2/3 stat weight, splitted in 2)
         di051(:) = i053(:) * 0.6667
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di051(:)
         nb_reaction_4 = nb_reaction_4 + 1
         v_4(:,nb_reaction_4) = di051(:)

      endif   !ionchem.and.deutchem

!----------------------------------------------------------------------
!     reactions avec 02(Dg) 
!----------------------------------------------------------------------

!---     j001: O2(Dg) + (CO2 and O) -> O2 + (CO2 and O) + hv

!        Krasnopolsky (2010a) for CO2 & Clark and Wayne, 1969 for O (JPL) 

      j001(:) = 1.E-20 * c(:,i_co2) !+ 2.E-16 * c(:,i_o)

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = j001(:)

!---     j002: O2(Dg) -> O2 + hv

!        Lafferty et al; (1998)

      j002(:) = 2.2E-4

      nb_phot = nb_phot + 1
      v_phot(:,nb_phot) = j002(:)

      !!! TEST: artificial increase of CO2 photodissociation
      if (tuneupperatm) then
      !-- TuneA
      !   v_phot(65:78,4) = v_phot(65:78,4)*10.
      !   v_phot(60:64,4) = v_phot(60:64,4)*3.
      !   v_phot(55:59,3) = v_phot(55:59,3)*2.
      !--
      !-- TuneB
      !   v_phot(65:78,4) = v_phot(65:78,4)*10.
      !   v_phot(55:59,3) = v_phot(55:59,3)*5.
      !--
      !-- TuneC
      ! VCD 1.1 tuning
      !   v_phot(65:78,4) = v_phot(65:78,4)*10.
      !   v_phot(52:59,3) = v_phot(52:59,3)*5.
      !--
      !-- TuneE
      ! VCD 2.0 tuning
      !  v_phot(65:nz,4) = v_phot(65:nz,4)*10. ! CO2 + hv ==> O(1D) + CO
      !--
      ! VCD 2.1 tuning
      !   v_phot(65:nz,4) = v_phot(65:nz,4)*8. ! CO2 + hv ==> O(1D) + CO
      !--
      ! VCD 2.4 tuning
          v_phot(65:nz,4) = v_phot(65:nz,4)*7. ! CO2 + hv ==> O(1D) + CO
      !--
      ! TEST TUNEH
      !     v_phot(62:74,4) = v_phot(62:74,4)*6.5 ! CO2 + hv ==> O(1D) + CO
      !--
      !   v_phot(:,4) = v_phot(:,4)*10.
      !do ij=3,4
      !   v_phot(:,ij) = v_phot(:,ij)*10.
      !end do
      end if
return
end subroutine krates

!======================================================================

 subroutine fill_matrix(ilev, mat, prod, loss, lossconc, c, nesp, nlayer,  &
                        nb_reaction_3_max, nb_reaction_4_max, nb_phot_max, &
                        v_phot, v_3, v_4)

!======================================================================
! filling of the jacobian matrix
!======================================================================

use types_asis

implicit none

! input

integer             :: ilev    ! level index
integer             :: nesp    ! number of species in the chemistry
integer, intent(in) :: nlayer  ! number of atmospheric layers
integer, intent(in) :: nb_reaction_3_max 
                               ! number of quadratic reactions
integer, intent(in) :: nb_reaction_4_max
                               ! number of bimolecular reactions
integer, intent(in) :: nb_phot_max
                               ! number of processes treated numerically as photodissociations

real (kind = 8), dimension(nlayer,nesp)              :: c    ! number densities
real (kind = 8), dimension(nlayer,      nb_phot_max) :: v_phot
real (kind = 8), dimension(nlayer,nb_reaction_3_max) :: v_3
real (kind = 8), dimension(nlayer,nb_reaction_4_max) :: v_4

! output

real (kind = 8), dimension(nesp,nesp), intent(out) :: mat  ! matrix
real (kind = 8), dimension(nesp), intent(out)      :: prod, loss, lossconc

! local

integer :: iesp
integer :: ind_phot_2,ind_phot_4,ind_phot_6
integer :: ind_3_2,ind_3_4,ind_3_6
integer :: ind_4_2,ind_4_4,ind_4_6,ind_4_8
integer :: iphot,i3,i4

real(kind = 8) :: eps, eps_4  ! implicit/explicit coefficient

! initialisations 

mat(:,:) = 0.
prod(:)  = 0.
loss(:)  = 0.
lossconc(:) = 0.

! photodissociations
! or reactions a + c -> b + c
! or reactions a + ice -> b + c
do iphot = 1,nb_phot_max

  ind_phot_2 = indice_phot(iphot)%z2
  ind_phot_4 = indice_phot(iphot)%z4
  ind_phot_6 = indice_phot(iphot)%z6

  mat(ind_phot_2,ind_phot_2) = mat(ind_phot_2,ind_phot_2) + indice_phot(iphot)%z1*v_phot(ilev,iphot)
  mat(ind_phot_4,ind_phot_2) = mat(ind_phot_4,ind_phot_2) - indice_phot(iphot)%z3*v_phot(ilev,iphot)
  mat(ind_phot_6,ind_phot_2) = mat(ind_phot_6,ind_phot_2) - indice_phot(iphot)%z5*v_phot(ilev,iphot)

  loss(ind_phot_2)     = loss(ind_phot_2) + indice_phot(iphot)%z1*v_phot(ilev,iphot)
  lossconc(ind_phot_2) = lossconc(ind_phot_2) + indice_phot(iphot)%z1*v_phot(ilev,iphot)*c(ilev,ind_phot_2)
  
  prod(ind_phot_4)     = prod(ind_phot_4) + indice_phot(iphot)%z3*v_phot(ilev,iphot)*c(ilev,ind_phot_2)
  prod(ind_phot_6)     = prod(ind_phot_6) + indice_phot(iphot)%z5*v_phot(ilev,iphot)*c(ilev,ind_phot_2)

end do

! reactions a + a -> b + c 

do i3 = 1,nb_reaction_3_max

  ind_3_2 = indice_3(i3)%z2
  ind_3_4 = indice_3(i3)%z4
  ind_3_6 = indice_3(i3)%z6

  mat(ind_3_2,ind_3_2) = mat(ind_3_2,ind_3_2) + indice_3(i3)%z1*v_3(ilev,i3)*c(ilev,ind_3_2)
  mat(ind_3_4,ind_3_2) = mat(ind_3_4,ind_3_2) - indice_3(i3)%z3*v_3(ilev,i3)*c(ilev,ind_3_2)
  mat(ind_3_6,ind_3_2) = mat(ind_3_6,ind_3_2) - indice_3(i3)%z5*v_3(ilev,i3)*c(ilev,ind_3_2)

  loss(ind_3_2)     = loss(ind_3_2) + indice_3(i3)%z1*v_3(ilev,i3)*c(ilev,ind_3_2)
  lossconc(ind_3_2) = lossconc(ind_3_2) + indice_3(i3)%z1*v_3(ilev,i3)*c(ilev,ind_3_2)*c(ilev,ind_3_2)
  
  prod(ind_3_4)     = prod(ind_3_4) + indice_3(i3)%z3*v_3(ilev,i3)*c(ilev,ind_3_2)*c(ilev,ind_3_2)
  prod(ind_3_6)     = prod(ind_3_6) + indice_3(i3)%z5*v_3(ilev,i3)*c(ilev,ind_3_2)*c(ilev,ind_3_2)

end do

! reactions a + b -> c + d 

eps = 1.d-10

do i4 = 1,nb_reaction_4_max

  ind_4_2 = indice_4(i4)%z2
  ind_4_4 = indice_4(i4)%z4
  ind_4_6 = indice_4(i4)%z6
  ind_4_8 = indice_4(i4)%z8

  eps_4 = abs(c(ilev,ind_4_2))/(abs(c(ilev,ind_4_2)) + abs(c(ilev,ind_4_4)) + eps)
  eps_4 = min(eps_4,1.0)

  mat(ind_4_2,ind_4_2) = mat(ind_4_2,ind_4_2) + indice_4(i4)%z1*v_4(ilev,i4)*(1. - eps_4)*c(ilev,ind_4_4) 
  mat(ind_4_2,ind_4_4) = mat(ind_4_2,ind_4_4) + indice_4(i4)%z1*v_4(ilev,i4)*eps_4*c(ilev,ind_4_2)
  mat(ind_4_4,ind_4_2) = mat(ind_4_4,ind_4_2) + indice_4(i4)%z3*v_4(ilev,i4)*(1. - eps_4)*c(ilev,ind_4_4)
  mat(ind_4_4,ind_4_4) = mat(ind_4_4,ind_4_4) + indice_4(i4)%z3*v_4(ilev,i4)*eps_4*c(ilev,ind_4_2)   
  mat(ind_4_6,ind_4_2) = mat(ind_4_6,ind_4_2) - indice_4(i4)%z5*v_4(ilev,i4)*(1. - eps_4)*c(ilev,ind_4_4)
  mat(ind_4_6,ind_4_4) = mat(ind_4_6,ind_4_4) - indice_4(i4)%z5*v_4(ilev,i4)*eps_4*c(ilev,ind_4_2)
  mat(ind_4_8,ind_4_2) = mat(ind_4_8,ind_4_2) - indice_4(i4)%z7*v_4(ilev,i4)*(1. - eps_4)*c(ilev,ind_4_4)
  mat(ind_4_8,ind_4_4) = mat(ind_4_8,ind_4_4) - indice_4(i4)%z7*v_4(ilev,i4)*eps_4*c(ilev,ind_4_2)


  loss(ind_4_2)     = loss(ind_4_2) + indice_4(i4)%z1*v_4(ilev,i4)*c(ilev,ind_4_4)
  lossconc(ind_4_2) = lossconc(ind_4_2) + indice_4(i4)%z1*v_4(ilev,i4)*c(ilev,ind_4_4)*c(ilev,ind_4_2)
  loss(ind_4_4)     = loss(ind_4_4) + indice_4(i4)%z3*v_4(ilev,i4)*c(ilev,ind_4_2)
  lossconc(ind_4_4) = lossconc(ind_4_4) + indice_4(i4)%z3*v_4(ilev,i4)*c(ilev,ind_4_2)*c(ilev,ind_4_4)
  
  prod(ind_4_6)     = prod(ind_4_6) + indice_4(i4)%z5*v_4(ilev,i4)*c(ilev,ind_4_2)*c(ilev,ind_4_4)
  prod(ind_4_8)     = prod(ind_4_8) + indice_4(i4)%z7*v_4(ilev,i4)*c(ilev,ind_4_2)*c(ilev,ind_4_4)

end do

end subroutine fill_matrix

!================================================================

 subroutine define_dt(nesp, dtnew, dtold, ctimestep, cold, ccur, mat1, &
                      prod, loss, dens, lon, lat)

!================================================================
! iterative evaluation of the appropriate time step dtnew
! according to curvature criterion based on
! e = 2 Rtol [r Cn+1 -(1-r) Cn + Cn-1 ]/[(1+r) Cn]
! with r = (tn - tn-1)/(tn+1 - tn)
!================================================================

implicit none

! input

integer :: nesp  ! number of species in the chemistry

real :: dtold, ctimestep
real (kind = 8), dimension(nesp)      :: cold, ccur
real (kind = 8), dimension(nesp,nesp) :: mat1
real (kind = 8), dimension(nesp)      :: prod, loss
real                       :: dens
real                       :: lon, lat

! output

real :: dtnew

! local

real (kind = 8), dimension(nesp)      :: cnew
real (kind = 8), dimension(nesp,nesp) :: mat
real (kind = 8) :: atol, ratio, e, es, coef

integer                  :: code, iesp, iter
integer, dimension(nesp) :: indx
integer :: imax

real :: dttest

! parameters

real (kind = 8), parameter :: dtmin   = 10.      ! minimum time step (s)
real (kind = 8), parameter :: vmrtol  = 1.e-11   ! absolute tolerance on vmr
real (kind = 8), parameter :: rtol    = 0.05     ! rtol recommended value : 0.1-0.02
integer,         parameter :: niter   = 3        ! number of iterations
real (kind = 8), parameter :: coefmax = 2.
real (kind = 8), parameter :: coefmin = 0.1 
logical                    :: fast_guess = .true.


dttest = dtold   ! dttest = dtold = dt_guess

atol = vmrtol*dens ! absolute tolerance in molecule.cm-3

do iter = 1,niter

if (fast_guess) then

! first guess : fast semi-implicit method

   do iesp = 1, nesp
      cnew(iesp) = (ccur(iesp) + prod(iesp)*dttest)/(1. + loss(iesp)*dttest)
   end do

else

! first guess : form the matrix identity + mat*dt_guess

   mat(:,:) = mat1(:,:)*dttest
   do iesp = 1,nesp
      mat(iesp,iesp) = 1. + mat(iesp,iesp)
   end do

! form right-hand side (RHS) of the system

   cnew(:) = ccur(:)

! solve the linear system  M*Cn+1 = Cn (RHS in cnew, then replaced by solution)

#ifdef LAPACK
   call dgesv(nesp,1,mat,nesp,indx,cnew,nesp,code)
#else
   write(*,*) "photochemistry error, missing LAPACK routine dgesv"
   stop
#endif

end if

! ratio old/new subtimestep

ratio = dtold/dttest

! e : local error indicator

e = 0.

do iesp = 1,nesp
   es = 2.*abs((ratio*cnew(iesp) - (1. + ratio)*ccur(iesp) + cold(iesp))   &
         /(1. + ratio)/max(ccur(iesp)*rtol,atol))

   if (es > e) then
      e = es
      imax = iesp
   end if
end do

! timestep correction

coef = max(coefmin, min(coefmax,0.8/sqrt(e)))

dttest = max(dtmin,dttest*coef)
dttest = min(ctimestep,dttest)

end do ! iter

! new timestep

dtnew = dttest

end subroutine define_dt

!======================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CE CODE EST OBSOLETE !! A NE PAS UTILISER !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE  rate_save(            &
                           n_lev,       &
                           pres,        &
                           temperature, &
                           traceur,     &
                           nq_max,      &
                           vphot,       &
                           v3,          &
                           v4)      
!==================
!!!!! MODEL 1D !!!! ==> n_lon = 1 !!!!
!==================
! Ici on a les variables pour le modele 1D, surtout pour la sauvegarde des taux de prod/consom
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!PENSER a changer les conditions de time_tot
!time_tot=nbr_pdt*(nbr_jour-1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      USE chemparam_mod
      IMPLICIT none
            

!INTEGER, PARAMETER :: time_tot=6000*1

INTEGER :: unit_loc, ierr_loc           ! unite de lecture de "rcm1d.def"
      
INTEGER, SAVE :: time_tot,nbr_pdt,nbr_jour
INTEGER, SAVE :: cpt_time, cpt_time_rate
DOUBLE PRECISION, DIMENSION(n_lev,126) :: rate_day
DOUBLE PRECISION, DIMENSION(n_lev,126) :: rate_night
DOUBLE PRECISION :: rate_local
DOUBLE PRECISION :: concentration(n_lev)
DOUBLE PRECISION :: pres(n_lev)
DOUBLE PRECISION :: temperature(n_lev)
DOUBLE PRECISION :: traceur(n_lev,nq_max)
      
INTEGER :: n_lev, nq_max
INTEGER :: i_lev, i_react, i_v

INTEGER :: i
 
LOGICAL, SAVE :: f_call = .true.

integer, parameter :: nb_phot_max = 30
integer, parameter :: nb_reaction_3_max = 12
integer, parameter :: nb_reaction_4_max = 87
      
real, dimension(n_lev,nb_phot_max) :: vphot
real, dimension(n_lev,nb_reaction_3_max) :: v3
real, dimension(n_lev,nb_reaction_4_max) :: v4

!PRINT*,"DEBUT subroutine rate_save" 


      IF (f_call) THEN
! ------------------------------------------------------
!  Lecture des parametres dans "rcm1d.def" 
! ------------------------------------------------------

!   Opening parameters file "rcm1d.def"
!   ---------------------------------------
      unit_loc =98
      OPEN(unit_loc,file='rcm1d.def',status='old',form='formatted'  &
          ,iostat=ierr_loc)

      IF(ierr_loc.ne.0) THEN
        write(*,*) 'Problem to open "rcm1d.def'
        write(*,*) 'Is it there ?'
        stop
      ELSE
        write(*,*) 'open rcm1d.def success '
      END IF

      do i=1, 2
        read (unit_loc, *)
      end do

      PRINT *,'nombre de pas de temps par jour ?'
      READ(unit_loc,*) nbr_pdt
      print*,nbr_pdt

      PRINT *,'nombre de jours simules ?'
      READ(unit_loc,*) nbr_jour
      print*,nbr_jour
      
 
      
      time_tot = nbr_pdt*(nbr_jour-1)
      PRINT *,'nombre de PdT avant calcul des taux production/consommation ?'
      PRINT*,time_tot
      
      PRINT*,'nlev',n_lev
           
         cpt_time = 1
         cpt_time_rate = 1
         f_call = .false.
         PRINT*,"f_call: ",f_call
         rate_night(:,:)=0.
         rate_day(:,:)=0.
      
      END IF           
      
!      PRINT*,"P	T"
!      PRINT*,pres,temperature
 	     
      IF (cpt_time .GE. time_tot) THEN

! 	PRINT*,'cpt_time',cpt_time
 	      
         DO i_lev=1, n_lev
         concentration(i_lev) = pres(i_lev)/(1.3806488E-19 * temperature(i_lev))     
         END DO
         
         IF (((cpt_time_rate .GE. 1).AND.(cpt_time_rate .LE. (nbr_pdt/4))).OR. &
         (cpt_time_rate .GT. (3*(nbr_pdt/4)))) THEN
         
!===============================
!        !!!! NUIT !!!!
!===============================
!	PRINT*,'NUIT'
	
           DO i_lev=1, n_lev
           i_react=1
           i_v=1
!===============================
!    1     o2 + hv     -> o + o
!===============================
		rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o2)*concentration(i_lev)
		rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
		i_react=i_react+1
		i_v=i_v+1
!===============================
!    2     o2 + hv     -> o + o(1d)
!===============================
		rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o2)*concentration(i_lev)
		rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
		i_react=i_react+1
		i_v=i_v+1
!===============================
!    3     co2 + hv    -> co + o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_co2)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    4     co2 + hv    -> co + o(1d)
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_co2)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    5     o3 + hv     -> o2(Dg) + o(1d)
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o3)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    6     o3 + hv     -> o2 + o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o3)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    7     h2o + hv    -> h + oh
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_h2o)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!    8     ho2 + hv    -> oh + o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_ho2)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    9     h2o2 + hv   -> oh + oh
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_h2o2)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    10    hcl + hv    -> h + cl
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_hcl)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    11    cl2 + hv    -> cl + cl
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_cl2)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    12    hocl + hv   -> oh + cl
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_hocl)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    13    so2 + hv    -> so + o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_so2)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    14    so + hv     -> s + o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    15    so3 + hv    -> so2 + o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_so3)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    16    clo + hv    -> cl + o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_clo)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    17    ocs + hv    -> co + s
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_ocs)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    18    cocl2 + hv  -> cl + cl + co
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_cocl2)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!    19    h2so4 + hv  -> so3 + h2o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_h2so4)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!--- 20 b001 o(1d) + co2 -> o + co2
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o1d)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 21 b004 o(1d) + o2 -> o + o2
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o1d)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1
!===============================
!--- 22 f014 clco + co2 -> cl + co + co2
!===============================
            rate_local = vphot(i_lev,22)*traceur(i_lev,i_clco)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 24 g023 s2 + co2 -> 2s + co2
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_s2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 25 h001 ICE
!===============================
            i_react=i_react+1
		i_v=i_v+1
!===============================
!--- 26 h002 ICE
!===============================
            i_react=i_react+1
		i_v=i_v+1  
!===============================
!--- 27 h003 ICE
!===============================
            i_react=i_react+1
		i_v=i_v+1  
!===============================
!---  30 i001 o2(Dg) + CO2 -> O2 + CO2 + hv
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o2dg)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!--- 31 i002 o2(Dg) -> O2 + hv
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o2dg)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1

! DEBUT DES REACTION V3
		i_v = i_v - nb_phot_max
!===============================
!--- 32 a002: o + o + co2 -> o2 + co2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_o)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 33 c008: ho2 + ho2 -> h2o2 + o2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_ho2)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 34 c013: oh + oh -> h2o + o
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_oh)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev)
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 35 c016: ho2 + ho2 + co2 -> h2o2 + o2 + co2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_ho2)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 36 c017: oh + oh + co2 -> h2o2 + co2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_oh)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 37 c018: h + h + co2 -> h2 + co2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_h)*concentration(i_lev) &
            *traceur(i_lev,i_h)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 38 f021: cl + cl + co2 -> cl2 + co2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_cl)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 39 f026: clco + clco  -> cocl2 + co
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_clco)*concentration(i_lev) &
            *traceur(i_lev,i_clco)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 40 f030: clso2 + clso2  -> cl2 + so2 + so2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_clso2)*concentration(i_lev) &
            *traceur(i_lev,i_clso2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 41 g015: so + so + co2 -> s2o2 + co2
!===============================
!           rate_local = v3(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
!           *traceur(i_lev,i_so)*concentration(i_lev) 
!           rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
!           i_react=i_react+1
!	i_v=i_v+1 
!===============================
!--- 42 g022: s + s + co2 -> s2 + co2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_s)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 

! DEBUT DES REACTION V4

		i_v = i_v - nb_reaction_3_max

!===============================
!--- 43 a001: o + o2 + co2 -> o3 + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o2)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 44 a003: o + o3 -> o2 + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o3)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 45 b002: o(1d) + h2o  -> oh + oh
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o1d)*concentration(i_lev) &
            *traceur(i_lev,i_h2o)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 46 b003: o(1d) + h2  -> oh + h
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o1d)*concentration(i_lev) &
            *traceur(i_lev,i_h2)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 47 b005: o(1d) + o3  -> o2 + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o1d)*concentration(i_lev) &
            *traceur(i_lev,i_o3)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 48 b006: o(1d) + o3  -> o2 + o + o
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o1d)*concentration(i_lev) &
            *traceur(i_lev,i_o3)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 49 c001: o + ho2 -> oh + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 50 c002: o + oh -> o2 + h
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 51 c003: h + o3 -> oh + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_h)*concentration(i_lev) &
            *traceur(i_lev,i_o3)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 52 c004: h + ho2 -> oh + oh
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_h)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 53 c005: h + ho2 -> h2 + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_h)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 54 c006: h + ho2 -> h2o + o
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_h)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 55 c007: oh + ho2 -> h2o + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_oh)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 56 c009: oh + h2o2 -> h2o + ho2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_oh)*concentration(i_lev) &
            *traceur(i_lev,i_h2o2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 57 c010: oh + h2 -> h2o + h
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_oh)*concentration(i_lev) &
            *traceur(i_lev,i_h2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 58 c011: h + o2 + co2 -> ho2 + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o2)*concentration(i_lev) &
            *traceur(i_lev,i_h)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 59 c012: o + h2o2 -> oh + ho2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o)*concentration(i_lev) &
            *traceur(i_lev,i_h2o2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 60 c014: oh + o3 -> ho2 + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_oh)*concentration(i_lev) &
            *traceur(i_lev,i_o3)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 61 c015: ho2 + o3 -> oh + o2 + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o3)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 62 e001: oh + co -> co2 + h
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 63 e002: o + co + m -> co2 + m
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o)*concentration(i_lev) &
            *traceur(i_lev,i_co)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 64 f001: hcl + o(1d) -> oh + cl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_hcl)*concentration(i_lev) &
            *traceur(i_lev,i_o1d)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 65 f002: hcl + o(1d) -> h + clo
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_hcl)*concentration(i_lev) &
            *traceur(i_lev,i_o1d)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 66 f003: hcl + o -> oh + cl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_hcl)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 67 f004: hcl + oh -> h2o + cl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_hcl)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 68 f005: clo + o -> cl + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clo)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 69 f006: clo + oh -> cl + ho2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clo)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 70 f007: clo + oh -> hcl + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clo)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 71 f008: cl + h2 -> hcl + h
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_h2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 72 f009: cl + o3 -> clo + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_o3)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 73 f010: cl + ho2 -> clo + oh
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 74 f011: cl + ho2 -> hcl + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 75 f012: cl + h2o2 -> hcl + ho2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_h2o2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 76 f013: cl + co + co2 -> clco + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_co)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 77 f015: clco + o2 + m -> clco3 + m
!===============================
		rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o2)*concentration(i_lev) &
            *traceur(i_lev,i_clco)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 78 & 79 f016: clco3 + cl -> cl + clo + co2
!===============================
!     decomposee en :
!     0.5 clco3 + 0.5 cl -> cl + 0.5 co2
!     0.5 clco3 + 0.5 cl -> clo + 0.5 co2
            rate_local = v4(i_lev,i_v)*0.25*traceur(i_lev,i_clco3)*concentration(i_lev) &
            *traceur(i_lev,i_cl)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
            
            rate_local = v4(i_lev,i_v)*0.25*traceur(i_lev,i_clco3)*concentration(i_lev) &
            *traceur(i_lev,i_cl)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 80 & 81 f017: clco3 + o -> cl + o2 + co2
!===============================
!     decomposee en :
!     0.5 clco3 + 0.5 o -> cl
!     0.5 clco3 + 0.5 o -> o2 + co2
            rate_local = v4(i_lev,i_v)*0.25*traceur(i_lev,i_clco3)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
            
            rate_local = v4(i_lev,i_v)*0.25*traceur(i_lev,i_clco3)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 82 f018: clo + ho2  -> hocl + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clo)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 83 f019: oh + hocl -> h2o + clo
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_oh)*concentration(i_lev) &
            *traceur(i_lev,i_hocl)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 84 f020: o + hocl -> oh + clo
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_hocl)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 85 f022: clco + o -> cl + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clco)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 86 f023: cl2 + o(1d) -> cl + clo
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl2)*concentration(i_lev) &
            *traceur(i_lev,i_o1d)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 87 f024: cl2 + h  -> hcl + cl
!==============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl2)*concentration(i_lev) &
            *traceur(i_lev,i_h)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 88 f025: cl + clco  -> cl2 + co
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_clco)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 89 f027: cl + so2 + co2  -> clso2 + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so2)*concentration(i_lev) &
            *traceur(i_lev,i_cl)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 90 f028: clso2 + o  -> so2 + clo
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clso2)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 91 f029: clso2 + h  -> so2 + hcl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clso2)*concentration(i_lev) &
            *traceur(i_lev,i_h)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 92 f031: cl + o + co2  -> clo + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 93 f032: cl2 + o -> clo + cl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl2)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 94 f033: clco + oh -> hocl + co
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clco)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 95 f034: cl2 + oh -> cl + hocl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl2)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 96 f035: clco + o -> co + clo
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clco)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 97 f036: clco + cl2 -> cocl2 + cl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl2)*concentration(i_lev) &
            *traceur(i_lev,i_clco)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 98 f037: hcl + h -> h2 + cl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_hcl)*concentration(i_lev) &
            *traceur(i_lev,i_h)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 99 f038: clco + h -> hcl + co
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s2)*concentration(i_lev) 
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 100 f039: cl + h + m -> hcl + m
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_h)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 101 g001: s + o2 -> so + o
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_o2)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 102 g002: s + o3 -> so + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_o3)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 103 g003: so + o2 -> so2 + o
!===============================
             rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_o2)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 104 g004: so + o3 -> so2 + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_o3)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 105 g005: so + oh -> so2 + h
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 106 g006: s + oh -> so + h
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 107 g007: so + o + co2 -> so2 + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 108 g008: so + ho2 -> so2 + oh 
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 109 g009: so2 + o + co2 -> so3 + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so2)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 110 g010: s + o + co2 -> so + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 111 g011: so3 + h2o -> h2so4
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so3)*concentration(i_lev) &
            *traceur(i_lev,i_h2o)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 112 g012: so + clo -> so2 + cl 
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_clo)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 113 g013: so + so3 -> so2 + so2 
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_so3)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 114 g014: so3 + o -> so2 + o2 
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so3)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 115 & 116 g017: clco3 + so -> cl + so2 + co2
!===============================
!     decomposee en :
!     0.5 clco3 + 0.5 so -> cl
!     0.5 clco3 + 0.5 so -> so2 + co2
            rate_local = v4(i_lev,i_v)*0.25*traceur(i_lev,i_clco3)*concentration(i_lev) &
            *traceur(i_lev,i_so)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
            
            rate_local = v4(i_lev,i_v)*0.25*traceur(i_lev,i_clco3)*concentration(i_lev) &
            *traceur(i_lev,i_so)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 117 g018: s + co + co2 -> ocs + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_co)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 118 g019: clco + s -> ocs + cl 
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clco)*concentration(i_lev) &
            *traceur(i_lev,i_s)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 119 g020: so2 + oh + co2 -> hso3 + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so2)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 120 g021: hso3 + o2 -> ho2 + so3 
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o2)*concentration(i_lev) &
            *traceur(i_lev,i_hso3)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 121 g024: s2 + o -> so + s 
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s2)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 122 g025: s + ocs -> s2 +  co
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_ocs)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1
!===============================
!--- 123 g026: ocs + o -> so + co
!=============================== 
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o)*concentration(i_lev) &
            *traceur(i_lev,i_ocs)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 124 g027: s + so3 -> so2 +  so
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_so3)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 125 g028: s + ho2 -> so +  oh
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 126 g029: s + clo -> so +  cl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_clo)*concentration(i_lev)  
            rate_night(i_lev,i_react) = rate_night(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1   
           END DO       
         ELSE
!===============================         
!        !!!! JOUR !!!!
!===============================
!	PRINT*,'JOUR'
	
           DO i_lev=1, n_lev
           i_react=1
           i_v=1
!===============================
!    1     o2 + hv     -> o + o
!===============================
		rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o2)*concentration(i_lev)
		rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
		i_react=i_react+1
		i_v=i_v+1
!===============================
!    2     o2 + hv     -> o + o(1d)
!===============================
		rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o2)*concentration(i_lev)
		rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
		i_react=i_react+1
		i_v=i_v+1
!===============================
!    3     co2 + hv    -> co + o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_co2)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    4     co2 + hv    -> co + o(1d)
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_co2)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    5     o3 + hv     -> o2(Dg) + o(1d)
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o3)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    6     o3 + hv     -> o2 + o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o3)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    7     h2o + hv    -> h + oh
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_h2o)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!    8     ho2 + hv    -> oh + o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_ho2)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    9     h2o2 + hv   -> oh + oh
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_h2o2)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    10    hcl + hv    -> h + cl
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_hcl)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    11    cl2 + hv    -> cl + cl
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_cl2)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    12    hocl + hv   -> oh + cl
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_hocl)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    13    so2 + hv    -> so + o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_so2)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    14    so + hv     -> s + o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    15    so3 + hv    -> so2 + o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_so3)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    16    clo + hv    -> cl + o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_clo)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    17    ocs + hv    -> co + s
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_ocs)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    18    cocl2 + hv  -> cl + cl + co
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_cocl2)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!    19    h2so4 + hv  -> so3 + h2o
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_h2so4)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    20     o(1d) + co2 -> o + co2
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o1d)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!    21    o(1d) + o2 -> o + o2
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o1d)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    22    clco + co2 -> cl + co + co2
!===============================
            rate_local = vphot(i_lev,22)*traceur(i_lev,i_clco)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!    23    s2o2 + co2 -> 2so + co2
!===============================
!           rate_local = vphot(i_lev,23)*traceur(i_lev,i_s2o2)*concentration(i_lev) 
!           rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
!           i_react=i_react+1
!           i_v=i_v+1 
!===============================
!    24    s2 + co2 -> 2s + co2
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_s2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!    25    ICE
!===============================
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    26    ICE
!===============================
            i_react=i_react+1
		i_v=i_v+1  
!===============================
!    27    ICE
!===============================
            i_react=i_react+1
		i_v=i_v+1  
!===============================
!    28    ICE
!===============================
            i_react=i_react+1
		i_v=i_v+1  
!===============================
!    29    ICE
!===============================
            i_react=i_react+1
		i_v=i_v+1  
!===============================
!    30    o2(Dg) + CO2 -> O2 + CO2
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o2dg)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
!===============================
!    31    o2(Dg) -> O2 + hv
!===============================
            rate_local = vphot(i_lev,i_v)*traceur(i_lev,i_o2dg)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt
            i_react=i_react+1
		i_v=i_v+1
		
! DEBUT DES REACTION V3
		i_v = i_v - nb_phot_max
!===============================
!--- 32 a002: o + o + co2 -> o2 + co2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_o)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 33 c008: ho2 + ho2 -> h2o2 + o2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_ho2)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 34 c013: oh + oh -> h2o + o
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_oh)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev)
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 35 c016: ho2 + ho2 + co2 -> h2o2 + o2 + co2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_ho2)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 36 c017: oh + oh + co2 -> h2o2 + co2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_oh)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 37 c018: h + h + co2 -> h2 + co2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_h)*concentration(i_lev) &
            *traceur(i_lev,i_h)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 38 f021: cl + cl + co2 -> cl2 + co2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_cl)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 39 f026: clco + clco  -> cocl2 + co
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_clco)*concentration(i_lev) &
            *traceur(i_lev,i_clco)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 40 f030: clso2 + clso2  -> cl2 + so2 + so2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_clso2)*concentration(i_lev) &
            *traceur(i_lev,i_clso2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 41 g015: so + so + co2 -> s2o2 + co2
!===============================
!           rate_local = v3(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
!           *traceur(i_lev,i_so)*concentration(i_lev) 
!           rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
!           i_react=i_react+1
!           i_v=i_v+1 
!===============================
!--- 42 g022: s + s + co2 -> s2 + co2
!===============================
            rate_local = v3(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_s)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 

! DEBUT DES REACTION V4

		i_v = i_v - nb_reaction_3_max

!===============================
!--- 43 a001: o + o2 + co2 -> o3 + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o2)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 44 a003: o + o3 -> o2 + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o3)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 45 b002: o(1d) + h2o  -> oh + oh
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o1d)*concentration(i_lev) &
            *traceur(i_lev,i_h2o)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 46 b003: o(1d) + h2  -> oh + h
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o1d)*concentration(i_lev) &
            *traceur(i_lev,i_h2)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 47 b005: o(1d) + o3  -> o2 + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o1d)*concentration(i_lev) &
            *traceur(i_lev,i_o3)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 48 b006: o(1d) + o3  -> o2 + o + o
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o1d)*concentration(i_lev) &
            *traceur(i_lev,i_o3)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 49 c001: o + ho2 -> oh + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 50 c002: o + oh -> o2 + h
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 51 c003: h + o3 -> oh + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_h)*concentration(i_lev) &
            *traceur(i_lev,i_o3)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 52 c004: h + ho2 -> oh + oh
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_h)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 53 c005: h + ho2 -> h2 + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_h)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 54 c006: h + ho2 -> h2o + o
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_h)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 55 c007: oh + ho2 -> h2o + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_oh)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 56 c009: oh + h2o2 -> h2o + ho2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_oh)*concentration(i_lev) &
            *traceur(i_lev,i_h2o2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 57 c010: oh + h2 -> h2o + h
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_oh)*concentration(i_lev) &
            *traceur(i_lev,i_h2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 58 c011: h + o2 + co2 -> ho2 + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o2)*concentration(i_lev) &
            *traceur(i_lev,i_h)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 59 c012: o + h2o2 -> oh + ho2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o)*concentration(i_lev) &
            *traceur(i_lev,i_h2o2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 60 c014: oh + o3 -> ho2 + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_oh)*concentration(i_lev) &
            *traceur(i_lev,i_o3)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 61 c015: ho2 + o3 -> oh + o2 + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o3)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 62 e001: oh + co -> co2 + h
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 63 e002: o + co + m -> co2 + m
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o)*concentration(i_lev) &
            *traceur(i_lev,i_co)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 64 f001: hcl + o(1d) -> oh + cl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_hcl)*concentration(i_lev) &
            *traceur(i_lev,i_o1d)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 65 f002: hcl + o(1d) -> h + clo
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_hcl)*concentration(i_lev) &
            *traceur(i_lev,i_o1d)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 66 f003: hcl + o -> oh + cl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_hcl)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 67 f004: hcl + oh -> h2o + cl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_hcl)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 68 f005: clo + o -> cl + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clo)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 69 f006: clo + oh -> cl + ho2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clo)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 70 f007: clo + oh -> hcl + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clo)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 71 f008: cl + h2 -> hcl + h
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_h2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 72 f009: cl + o3 -> clo + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_o3)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 73 f010: cl + ho2 -> clo + oh
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 74 f011: cl + ho2 -> hcl + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 75 f012: cl + h2o2 -> hcl + ho2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_h2o2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 76 f013: cl + co + co2 -> clco + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_co)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 77 f015: clco + o2 + m -> clco3 + m
!===============================
		rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o2)*concentration(i_lev) &
            *traceur(i_lev,i_clco)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 78 & 79 f016: clco3 + cl -> cl + clo + co2
!===============================
!     decomposee en :
!     0.5 clco3 + 0.5 cl -> cl + 0.5 co2
!     0.5 clco3 + 0.5 cl -> clo + 0.5 co2
            rate_local = v4(i_lev,i_v)*0.25*traceur(i_lev,i_clco3)*concentration(i_lev) &
            *traceur(i_lev,i_cl)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
            
            rate_local = v4(i_lev,i_v)*0.25*traceur(i_lev,i_clco3)*concentration(i_lev) &
            *traceur(i_lev,i_cl)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 80 & 81 f017: clco3 + o -> cl + o2 + co2
!===============================
!     decomposee en :
!     0.5 clco3 + 0.5 o -> cl
!     0.5 clco3 + 0.5 o -> o2 + co2
            rate_local = v4(i_lev,i_v)*0.25*traceur(i_lev,i_clco3)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
            
            rate_local = v4(i_lev,i_v)*0.25*traceur(i_lev,i_clco3)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 82 f018: clo + ho2  -> hocl + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clo)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 83 f019: oh + hocl -> h2o + clo
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_oh)*concentration(i_lev) &
            *traceur(i_lev,i_hocl)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 84 f020: o + hocl -> oh + clo
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_hocl)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 85 f022: clco + o -> cl + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clco)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 86 f023: cl2 + o(1d) -> cl + clo
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl2)*concentration(i_lev) &
            *traceur(i_lev,i_o1d)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 87 f024: cl2 + h  -> hcl + cl
!==============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl2)*concentration(i_lev) &
            *traceur(i_lev,i_h)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 88 f025: cl + clco  -> cl2 + co
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_clco)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 89 f027: cl + so2 + co2  -> clso2 + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so2)*concentration(i_lev) &
            *traceur(i_lev,i_cl)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 90 f028: clso2 + o  -> so2 + clo
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clso2)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 91 f029: clso2 + h  -> so2 + hcl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clso2)*concentration(i_lev) &
            *traceur(i_lev,i_h)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 92 f031: cl + o + co2  -> clo + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 93 f032: cl2 + o -> clo + cl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl2)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 94 f033: clco + oh -> hocl + co
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clco)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 95 f034: cl2 + oh -> cl + hocl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl2)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 96 f035: clco + o -> co + clo
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clco)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 97 f036: clco + cl2 -> cocl2 + cl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl2)*concentration(i_lev) &
            *traceur(i_lev,i_clco)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 98 f037: hcl + h -> h2 + cl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_hcl)*concentration(i_lev) &
            *traceur(i_lev,i_h)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 99 f038: clco + h -> hcl + co
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s2)*concentration(i_lev) 
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 100 f039: cl + h + m -> hcl + m
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_cl)*concentration(i_lev) &
            *traceur(i_lev,i_h)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 101 g001: s + o2 -> so + o
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_o2)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 102 g002: s + o3 -> so + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_o3)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 103 g003: so + o2 -> so2 + o
!===============================
             rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_o2)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 104 g004: so + o3 -> so2 + o2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_o3)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 105 g005: so + oh -> so2 + h
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 106 g006: s + oh -> so + h
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 107 g007: so + o + co2 -> so2 + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 108 g008: so + ho2 -> so2 + oh 
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 109 g009: so2 + o + co2 -> so3 + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so2)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 110 g010: s + o + co2 -> so + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 111 g011: so3 + h2o -> h2so4
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so3)*concentration(i_lev) &
            *traceur(i_lev,i_h2o)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 112 g012: so + clo -> so2 + cl 
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_clo)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 113 g013: so + so3 -> so2 + so2 
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so)*concentration(i_lev) &
            *traceur(i_lev,i_so3)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 114 g014: so3 + o -> so2 + o2 
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so3)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 115 & 116 g017: clco3 + so -> cl + so2 + co2
!===============================
!     decomposee en :
!     0.5 clco3 + 0.5 so -> cl
!     0.5 clco3 + 0.5 so -> so2 + co2
            rate_local = v4(i_lev,i_v)*0.25*traceur(i_lev,i_clco3)*concentration(i_lev) &
            *traceur(i_lev,i_so)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
            
            rate_local = v4(i_lev,i_v)*0.25*traceur(i_lev,i_clco3)*concentration(i_lev) &
            *traceur(i_lev,i_so)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 117 g018: s + co + co2 -> ocs + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_co)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 118 g019: clco + s -> ocs + cl 
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_clco)*concentration(i_lev) &
            *traceur(i_lev,i_s)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 119 g020: so2 + oh + co2 -> hso3 + co2
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_so2)*concentration(i_lev) &
            *traceur(i_lev,i_oh)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 120 g021: hso3 + o2 -> ho2 + so3 
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o2)*concentration(i_lev) &
            *traceur(i_lev,i_hso3)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 121 g024: s2 + o -> so + s 
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s2)*concentration(i_lev) &
            *traceur(i_lev,i_o)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 122 g025: s + ocs -> s2 +  co
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_ocs)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1
!===============================
!--- 123 g026: ocs + o -> so + co
!=============================== 
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_o)*concentration(i_lev) &
            *traceur(i_lev,i_ocs)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 124 g027: s + so3 -> so2 +  so
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_so3)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 125 g028: s + ho2 -> so +  oh
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_ho2)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1 
!===============================
!--- 126 g029: s + clo -> so +  cl
!===============================
            rate_local = v4(i_lev,i_v)*traceur(i_lev,i_s)*concentration(i_lev) &
            *traceur(i_lev,i_clo)*concentration(i_lev)  
            rate_day(i_lev,i_react) = rate_day(i_lev,i_react) + 2.*rate_local/nbr_pdt 
            i_react=i_react+1
		i_v=i_v+1
           END DO
         END IF
         cpt_time_rate = cpt_time_rate + 1
      
      END IF
            
      IF (cpt_time .EQ. (time_tot+nbr_pdt)) THEN
               OPEN(100,file='profile_rate_day.csv')
               DO i_lev=1,n_lev
               write (100,"(128(e15.8,','))")pres(i_lev), temperature(i_lev), (rate_day(i_lev,i_react),i_react=1,126)
               END DO
               
               OPEN(101,file='profile_rate_night.csv')
               DO i_lev=1,n_lev
               write (101,"(128(e15.8,','))") pres(i_lev), temperature(i_lev), (rate_night(i_lev,i_react),i_react=1,126)
               END DO
               
               OPEN(102,file='profile_rate_fullday.csv')
               rate_day=(rate_day+rate_night)/2.
               DO i_lev=1,n_lev
               write (102,"(128(e15.8,','))") pres(i_lev), temperature(i_lev), (rate_day(i_lev,i_react),i_react=1,126)
               END DO
               
               PRINT*,"pression top",pres(n_lev)
               PRINT*,"temp top",temperature(n_lev)
               
      END IF
      
      cpt_time = cpt_time + 1
      
      END     
