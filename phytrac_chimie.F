      subroutine phytrac_chimie (
     $                    debutphy,
     $                    gmtime,
     $                    nqmax,
     $                    nlon,
     $                    lat,
     $                    lon,
     $                    zzlay,
     $                    nlev, 
     $                    pdtphys,
     $                    temp,
     $                    pplay,
     $                    trac,
     $                    d_tr_chem,
     $                    iter,
     $                    prod_tr,
     $                    loss_tr,
     $                    no_emi,
     $                    o2_emi)

      use chemparam_mod
      use conc, only: mmean,rnew
      use photolysis_mod, only: init_photolysis, nphot
      use iono_h, only: temp_elect
      use clesphys_mod
      use YOMCST_mod
#ifdef CPP_XIOS      
      use xios_output_mod, only: send_xios_field
#endif

      implicit none
      
!#include "clesphys.h" 
!#include "YOMCST.h"

!===================================================================
!     input
!===================================================================

      integer :: nlon, nlev                     ! number of gridpoints and levels
      integer :: nqmax                          ! number of tracers

      real :: gmtime                            ! day fraction
      real :: pdtphys                           ! phytrac_chimie timestep (s)
      real, dimension(nlon,nlev) :: temp        ! temperature (k)
      real, dimension(nlon,nlev) :: pplay       ! pressure (pa)
      real, dimension(nlon,nlev) :: zzlay       ! altitude (m)
      real, dimension(nlon,nlev,nqmax) :: trac  ! tracer mass mixing ratio

      logical :: debutphy                       ! first call flag

!===================================================================
!     output
!===================================================================

      integer, dimension(nlon,nlev) :: iter                 ! chemical iterations
      real, dimension(nlon,nlev,nqmax) :: d_tr_chem         ! chemical tendency for each tracer
      real, dimension(nlon,nlev,nqmax) :: prod_tr, loss_tr  ! production and loss terms (for info)
      real, dimension(nlon,nlev)    :: no_emi               ! no emission
      real, dimension(nlon,nlev)    :: o2_emi               ! o2 emission

!===================================================================
!     local
!===================================================================

      real :: sza_local     ! solar zenith angle (deg)
      real :: lon_sun
      real :: zlocal(nlev)  ! altitude for photochem (km)
      real :: t_elec(nlev)  ! electron temperature [K]
      
      integer, parameter :: t_elec_origin=2
      !Electronic temperature. Index 1 -> Theis et al. 1980 - model data ; Index 2-> Theis et al. 1984 - model data
      
      integer :: i, iq
      integer :: ilon, ilev

      real  lat(nlon), lat_local(nlon)
      real  lon(nlon), lon_local(nlon)

      real, dimension(nlon,nlev) :: mrtwv, mrtsa, mrtdw, mrtdsa  ! total water and total sulfuric acid
      real, dimension(nlon,nlev) :: mrwv, mrsa, mrdw, mrdsa     ! gas-phase water and gas-phase sulfuric acid
      real, dimension(nlon,nlev) :: trac_sum
      real, dimension(nlon,nlev,nqmax) :: ztrac  ! local tracer mixing ratio
      real, dimension(nlev) :: no_em
      real, dimension(nlev) :: o2_em

      integer, save :: nb_reaction_3_max     ! number of quadratic reactions
      integer, save :: nb_reaction_4_max     ! number of bimolecular reactions
      integer, save :: nquench               ! number of quenching + heterogeneous reactions
      integer, save :: nphotion              ! number of photoionizations
      integer, save :: nb_reaction_4_ion     ! quadratic reactions for ionosphere
!      integer, save :: nb_reaction_4_deut    ! quadratic reactions for deuterium chem
      integer, save :: nb_phot_max           ! total number of photolysis+photoionizations+quenching reactions

      ! tracer indexes for the EUV heating:
!!! ATTENTION. These values have to be identical to those in euvheat.F90
!!! If the values are changed there, the same has to be done here  !!!

      integer,parameter :: ix_co2  =  1
      integer,parameter :: ix_co   =  2
      integer,parameter :: ix_o    =  3
      integer,parameter :: ix_o1d  =  4
      integer,parameter :: ix_o2   =  5
      integer,parameter :: ix_o3   =  6
      integer,parameter :: ix_h    =  7
      integer,parameter :: ix_h2   =  8
      integer,parameter :: ix_oh   =  9
      integer,parameter :: ix_ho2  = 10
      integer,parameter :: ix_h2o2 = 11
      integer,parameter :: ix_h2o  = 12
      integer,parameter :: ix_n    = 13
      integer,parameter :: ix_n2d  = 14
      integer,parameter :: ix_no   = 15
      integer,parameter :: ix_no2  = 16
      integer,parameter :: ix_n2   = 17

      ! NEED TO BE THE SAME AS IN EUVHEAT.F90
      integer,parameter :: nespeuv = 17    ! Number of species considered (11, 12 or 17 (with nitrogen))

      real :: vmr_dens_euv(nlon,nlev,nespeuv) ! local species density for EUV heating
          
!===================================================================
!     first call : initialisations
!===================================================================

      if (debutphy) then
      
!--- Adjustment of Helium amount
!       if (i_he/=0) then
!          trac(:,:,i_he)=trac(:,:,i_he)/2.
!       endif
!---

!-------------------------------------------------------------------
!        Determination of chemistry reaction number 
!-------------------------------------------------------------------

         if (ok_chem) then

           ! set number of reactions, depending on ion chemistry or not

           nb_reaction_4_ion  = 106
           !nb_reaction_4_deut = 35 
    
           ! default numbers if no ion and no deuterium chemistry
    
           nb_reaction_4_max = 170
           nb_reaction_3_max = 14     ! set number of quadratic reactions
           nquench           = 17
                                      !               quenching 
                                      !               thermal dissociation
                                      !               heterogeneous
           nphotion          = 0      ! set number of photoionizations
    
           if (ok_ionchem) then
              nb_reaction_4_max = nb_reaction_4_max + nb_reaction_4_ion
              nphotion       = 19     ! set number of photoionizations
           endif

           !if(deutchem) then
           !   nb_reaction_4_max = nb_reaction_4_max + nb_reaction_4_deut  
           !end if
    
           ! nb_phot_max is the total number of processes that are treated
           ! numerically as a photolysis (first-order reaction):
    
           nb_phot_max = nphot + nphotion + nquench
    
!-------------------------------------------------------------------
!        case of tracers re-initialisation with chemistry
!-------------------------------------------------------------------
  
           if (reinit_trac .and. ok_chem) then
  
!             convert mass to volume mixing ratio
              do iq = 1,nqmax - nmicro
                 trac(:,:,iq) = trac(:,:,iq)*mmean(:,:)/m_tr(iq)
              end do

              print*, "SO2 is re-initialised"

              if (i_so2 /= 0) then 
                 trac(:,25:,i_so2)  = 100.e-9
                 trac(:,1:24,i_so2) = 15.e-6
                 trac(:,25:,i_h2o)  = 1.e-6
                 trac(:,1:24,i_h2o) = 30.e-6
                 trac(:,:,i_hd) = trac(:,:,i_h2)*0.032
                 trac(:,:,i_hdo) = trac(:,:,i_h2o)*0.032
                 trac(:,:,i_d) = trac(:,:,i_h)*0.032
                 trac(:,:,i_od) = trac(:,:,i_oh)*0.032
                 trac(:,:,i_do2) = trac(:,:,i_ho2)*0.032
                 trac(:,:,i_hdo2) = trac(:,:,i_h2o2)*0.032
                 trac(:,:,i_hdoliq) = trac(:,:,i_h2oliq)*0.032
                 trac(:,:,i_hdso4) = trac(:,:,i_h2so4)*0.032
                 trac(:,:,i_hdso4liq) = trac(:,:,i_h2so4liq)*0.032
                 trac(:,:,i_dcl) = trac(:,:,i_hcl)*0.032
                 trac(:,:,i_docl) = trac(:,:,i_hocl)*0.032
 
                 if(ok_ionchem) then
                    trac(:,:,i_dco2plus) = trac(:,:,i_hco2plus)*0.032
                    trac(:,:,i_dcoplus) = trac(:,:,i_hcoplus)*0.032
                    trac(:,:,i_hdoplus) = trac(:,:,i_h2oplus)*0.032
                    trac(:,:,i_dplus) = trac(:,:,i_hplus)*0.032
                    trac(:,:,i_odplus) = trac(:,:,i_ohplus)*0.032
                    trac(:,:,i_h2doplus) = trac(:,:,i_h3oplus)*0.032
                 endif
 

                 trac(:,:,i_osso_cis)   = 0.
                 trac(:,:,i_osso_trans) = 0.
                 trac(:,:,i_s2o2_cyc)   = 0.
                 trac(:,:,i_cl2so2)     = 0.
  
! AL TRACERS!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           print*, "Tracers are re-initialised"
!           trac(:,:,:) = 1.0e-30
  
!             if ((i_ocs /= 0) .and. (i_co /= 0) .and. (i_hcl /= 0) 
!     $          .and. (i_so2 /= 0) .and. (i_h2o /= 0) 
!     $          .and. (i_n2 /= 0)  .and. (i_co2 /= 0)) then
  
!                trac(:,1:22,i_ocs) = 3.e-6
!                trac(:,1:22,i_co)  = 25.e-6
!                trac(:,:,i_hcl)    = 0.4e-6
!                trac(:,1:22,i_so2) = 7.e-6
!                trac(:,1:22,i_h2o) = 30.e-6
!                trac(:,:,i_n2)     = 0.35e-1
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
!                ensure that sum of mixing ratios = 1
  
                 trac_sum(:,:) = 0.
  
                 do iq = 1,nqmax - nmicro
                    if (iq /= i_co2) then
                       trac_sum(:,:) = trac_sum(:,:) + trac(:,:,iq) 
                    end if
                 end do
  
!                initialise co2
  
                 trac(:,:,i_co2) = 1. - trac_sum(:,:)
  
              else
                 write(*,*) "at least one tracer is missing : stop"
                 stop
              end if
          
!             update mmean
  
              mmean(:,:) = 0.
              do iq = 1,nqmax - nmicro
                 mmean(:,:) = mmean(:,:)+trac(:,:,iq)*m_tr(iq)
              enddo
              rnew(:,:) = 8.314/mmean(:,:)*1.e3     ! J/kg K 
  
!             convert volume to mass mixing ratio
  
              do iq = 1,nqmax - nmicro
                 trac(:,:,iq) = trac(:,:,iq)*m_tr(iq)/mmean(:,:)
              end do
    
           end if  ! reinit_trac

         end if  ! ok_chem

!-------------------------------------------------------------------
!        case of detailed microphysics without chemistry
!-------------------------------------------------------------------

         if (.not. ok_chem .and. ok_cloud .and. cl_scheme == 2) then

!           convert mass to volume mixing ratio

            do iq = 1,nqmax - nmicro
               ztrac(:,:,iq) = trac(:,:,iq)*mmean(:,:)/m_tr(iq)
            end do
          
!           initialise microphysics
 
            call vapors4muphy_ini(nlon,nlev,ztrac)

!           convert volume to mass mixing ratio

            do iq = 1,nqmax - nmicro
               trac(:,:,iq) = ztrac(:,:,iq)*m_tr(iq)/mmean(:,:)
            end do
    
         end if

      end if  ! debutphy 

!===================================================================
!     convert mass to volume mixing ratio : gas phase
!===================================================================

      do iq = 1,nqmax - nmicro
         ztrac(:,:,iq) = max(trac(:,:,iq)*mmean(:,:)/m_tr(iq), 1.e-30)
      end do

!===================================================================
!     microphysics: simplified scheme (phd aurelien stolzenbach)
!===================================================================

      if (ok_cloud .and. cl_scheme == 1) then

!     convert mass to volume mixing ratio : liquid phase

         ztrac(:,:,i_h2so4liq) = max(trac(:,:,i_h2so4liq)
     $                             *mmean(:,:)/m_tr(i_h2so4liq), 1.e-30)
         ztrac(:,:,i_hdso4liq) = max(trac(:,:,i_hdso4liq)
     $                             *mmean(:,:)/m_tr(i_hdso4liq), 1.e-30)
         ztrac(:,:,i_h2oliq) = max(trac(:,:,i_h2oliq)
     $                             *mmean(:,:)/m_tr(i_h2oliq), 1.e-30)
         ztrac(:,:,i_hdoliq) = max(trac(:,:,i_hdoliq)
     $                             *mmean(:,:)/m_tr(i_hdoliq), 1.e-30)
              
!     total water and sulfuric acid (gas + liquid)

         mrtwv(:,:) = ztrac(:,:,i_h2o) + ztrac(:,:,i_h2oliq)
         mrtdw(:,:) = ztrac(:,:,i_hdo) + ztrac(:,:,i_hdoliq)
         mrtsa(:,:) = ztrac(:,:,i_h2so4) + ztrac(:,:,i_h2so4liq)
         mrtdsa(:,:) = ztrac(:,:,i_hdso4) + ztrac(:,:,i_hdso4liq)

!     all water and sulfuric acid is put in the gas-phase

         mrwv(:,:) = mrtwv(:,:)
         mrdw(:,:) = mrtdw(:,:)
         mrsa(:,:) = mrtsa(:,:)
         mrdsa(:,:) = mrtdsa(:,:)

!     call microphysics

          call new_cloud_venus(nb_mode, nlev, nlon, temp, pplay,
     $      mrtwv, mrtdw, mrtsa, mrtdsa, mrwv, mrdw, mrsa, mrdsa)

!     update water vapour and sulfuric acid

         ztrac(:,:,i_h2o) = mrwv(:,:)
         ztrac(:,:,i_h2oliq) = mrtwv(:,:) - ztrac(:,:,i_h2o)
 
         ztrac(:,:,i_hdo) = mrdw(:,:)
         ztrac(:,:,i_hdoliq) = mrtdw(:,:) - ztrac(:,:,i_hdo)
     
         ztrac(:,:,i_h2so4) = mrsa(:,:)
         ztrac(:,:,i_h2so4liq) = mrtsa(:,:) - ztrac(:,:,i_h2so4)
 
         ztrac(:,:,i_hdso4) = mrdsa(:,:)
         ztrac(:,:,i_hdso4liq) = mrtdsa(:,:) - ztrac(:,:,i_hdso4)

      end if  ! simplified scheme

!===================================================================
!     microphysics: detailed scheme (phd sabrina guilbon)
!     !!! to be confirmed whether mad_muphy expects mmr or vmr for h2o and h2so4
!===================================================================

      if (ok_cloud .and. cl_scheme == 2) then

         do iq = nqmax-nmicro+1,nqmax
            ztrac(:,:,iq) = trac(:,:,iq)
         end do

         do ilon = 1,nlon       
            do ilev = 1, nlev
               if (temp(ilon,ilev) < 500.) then
                  call mad_muphy(pdtphys,                               ! timestep
     $                           temp(ilon,ilev),pplay(ilon,ilev),      ! temperature and pressure
     $                           ztrac(ilon,ilev,i_h2o),
     $                           ztrac(ilon,ilev,i_h2so4),      
     $                           ztrac(ilon,ilev,i_m0_aer),
     $                           ztrac(ilon,ilev,i_m3_aer),   
     $                           ztrac(ilon,ilev,i_m0_mode1drop),
     $                           ztrac(ilon,ilev,i_m0_mode1ccn), 
     $                           ztrac(ilon,ilev,i_m3_mode1sa),
     $                           ztrac(ilon,ilev,i_m3_mode1w),    
     $                           ztrac(ilon,ilev,i_m3_mode1ccn),   
     $                           ztrac(ilon,ilev,i_m0_mode2drop),
     $                           ztrac(ilon,ilev,i_m0_mode2ccn),
     $                           ztrac(ilon,ilev,i_m3_mode2sa),
     $                           ztrac(ilon,ilev,i_m3_mode2w), 
     $                           ztrac(ilon,ilev,i_m3_mode2ccn))
               else
                  ztrac(ilon,ilev,i_m0_aer)       = 0. 
                  ztrac(ilon,ilev,i_m3_aer)       = 0.
                  ztrac(ilon,ilev,i_m0_mode1drop) = 0.
                  ztrac(ilon,ilev,i_m0_mode1ccn)  = 0.
                  ztrac(ilon,ilev,i_m3_mode1sa)   = 0.
                  ztrac(ilon,ilev,i_m3_mode1w)    = 0.
                  ztrac(ilon,ilev,i_m3_mode1ccn)  = 0.
                  ztrac(ilon,ilev,i_m0_mode2drop) = 0.
                  ztrac(ilon,ilev,i_m0_mode2ccn)  = 0.
                  ztrac(ilon,ilev,i_m3_mode2sa)   = 0.
                  ztrac(ilon,ilev,i_m3_mode2w)    = 0.
                  ztrac(ilon,ilev,i_m3_mode2ccn)  = 0.
               end if
            end do
         end do

      end if  ! detailed scheme
            
!===================================================================
!     photochemistry
!===================================================================

      if (ok_chem) then

         lon_sun = (0.5 - gmtime)*2.*rpi
         lon_local(:) = lon(:)*rpi/180.
         lat_local(:) = lat(:)*rpi/180.
         
         if (ok_ionchem) then
        
           vmr_dens_euv(:,:,ix_co2) = ztrac(:,:,i_co2)  ! CO2
           vmr_dens_euv(:,:,ix_co)  = ztrac(:,:,i_co)   ! CO
           vmr_dens_euv(:,:,ix_o)   = ztrac(:,:,i_o)    ! O
           vmr_dens_euv(:,:,ix_o1d) = ztrac(:,:,i_o1d)  ! O(1D)
           vmr_dens_euv(:,:,ix_o2)  = ztrac(:,:,i_o2)   ! O2
           vmr_dens_euv(:,:,ix_o3)  = ztrac(:,:,i_o3)   ! O3
           vmr_dens_euv(:,:,ix_h)   = ztrac(:,:,i_h)    ! H
           vmr_dens_euv(:,:,ix_h2)  = ztrac(:,:,i_h2)   ! H2
           vmr_dens_euv(:,:,ix_oh)  = ztrac(:,:,i_oh)   ! OH
           vmr_dens_euv(:,:,ix_ho2) = ztrac(:,:,i_ho2)  ! HO2
           vmr_dens_euv(:,:,ix_h2o2)= ztrac(:,:,i_h2o2) ! H2O2
           vmr_dens_euv(:,:,ix_h2o) = ztrac(:,:,i_h2o)  ! H2O
           vmr_dens_euv(:,:,ix_n)   = ztrac(:,:,i_n)    ! N
           vmr_dens_euv(:,:,ix_n2d) = ztrac(:,:,i_n2d)  ! N(2D)
           vmr_dens_euv(:,:,ix_no)  = ztrac(:,:,i_no)   ! NO
           vmr_dens_euv(:,:,ix_no2) = ztrac(:,:,i_no2)  ! NO2
           vmr_dens_euv(:,:,ix_n2)  = ztrac(:,:,i_n2)   ! N2
           
         end if
         
         do ilon = 1,nlon
            zlocal(:)=zzlay(ilon,:)/1000.
                        
!           solar zenith angle
!            sza_local = acos(cos(lat_local(ilon))*cos(lon_local(ilon))
!     $                 *cos(lon_sun) + cos(lat_local(ilon))
!     $                 *sin(lon_local(ilon))*sin(lon_sun))*180./rpi

            sza_local = cos(lat_local(ilon))*cos(lon_local(ilon))
     $                 *cos(lon_sun) + cos(lat_local(ilon))
     $                 *sin(lon_local(ilon))*sin(lon_sun)

            ! Security - Handle rare cases where |sza_local| > 1
            sza_local = min(sza_local,1.)
            sza_local = max(-1.,sza_local)
            sza_local = acos(sza_local)*180./rpi

!           electron temperature 

            do ilev = 1, nlev
               t_elec(ilev) = temp_elect(zlocal(ilev),temp(ilon,ilev),
     $                                   sza_local,t_elec_origin)  
            end do

            call photochemistry_venus(nlev, nlon, zlocal, pdtphys,
     $                           ok_jonline,ok_ionchem,tuneupperatm,
     $                           nb_reaction_3_max,nb_reaction_4_max,
     $                           nb_phot_max,nphotion,ilon,          
     $                           pplay(ilon,:)/100.,
     $                           temp(ilon,:), t_elec(:),
     $                           ztrac(ilon,:,:),vmr_dens_euv(ilon,:,:),
     $                           mmean(ilon,:),
     $                           sza_local,
     $                           lon(ilon), lat(ilon),
     $                           nqmax, nespeuv, iter(ilon,:),
     $                           prod_tr(ilon,:,:),
     $                           loss_tr(ilon,:,:),
     $                           no_em, o2_em)

             no_emi(ilon,:) = no_em(:)
             o2_emi(ilon,:) = o2_em(:)

         end do         
      end if  ! ok_chem

!===================================================================
!     compute tendencies
!===================================================================

!     update mmean

      mmean(:,:) = 0.
      do iq = 1,nqmax - nmicro
         mmean(:,:) = mmean(:,:)+ztrac(:,:,iq)*m_tr(iq)
      end do
      rnew(:,:) = 8.314/mmean(:,:)*1.e3     ! J/kg K 

!===================================================================
!     convert volume to mass mixing ratio / then tendencies in mmr
!===================================================================

!     gas phase

      do iq = 1,nqmax - nmicro
         ztrac(:,:,iq) = max(ztrac(:,:,iq)*m_tr(iq)/mmean(:,:), 
     $                       1.e-30)
         d_tr_chem(:,:,iq) = (ztrac(:,:,iq) - trac(:,:,iq))/pdtphys
      end do

!     liquid phase or moments

      if (ok_cloud .and. cl_scheme == 1) then
         ztrac(:,:,i_h2so4liq) = max(ztrac(:,:,i_h2so4liq)
     $                               *m_tr(i_h2so4liq)/mmean(:,:),
     $                               1.e-30)
         ztrac(:,:,i_hdso4liq) = max(ztrac(:,:,i_hdso4liq)
     $                                  *m_tr(i_hdso4liq)/mmean(:,:),
     $                                  1.e-30)
         ztrac(:,:,i_h2oliq)   = max(ztrac(:,:,i_h2oliq)
     $                               *m_tr(i_h2oliq)/mmean(:,:),
     $                               1.e-30)
         ztrac(:,:,i_hdoliq)   = max(ztrac(:,:,i_hdoliq)
     $                                  *m_tr(i_hdoliq)/mmean(:,:),
     $                                  1.e-30)
         d_tr_chem(:,:,i_h2so4liq) = (ztrac(:,:,i_h2so4liq) 
     $                              - trac(:,:,i_h2so4liq))/pdtphys
         d_tr_chem(:,:,i_hdso4liq) = (ztrac(:,:,i_hdso4liq)
     $                              - trac(:,:,i_hdso4liq))/pdtphys
         d_tr_chem(:,:,i_h2oliq) = (ztrac(:,:,i_h2oliq) 
     $                            - trac(:,:,i_h2oliq))/pdtphys
         d_tr_chem(:,:,i_hdoliq) = (ztrac(:,:,i_hdoliq)
     $                            - trac(:,:,i_hdoliq))/pdtphys
      end if


      if (ok_cloud .and. cl_scheme == 2) then
         do iq = nqmax-nmicro+1,nqmax
            d_tr_chem(:,:,iq) = (ztrac(:,:,iq) - trac(:,:,iq))/pdtphys
         end do
      end if

      end subroutine phytrac_chimie
