!============================================================================

module chemparam_mod

!============================================================================

! 1) cloud microphysical parameters for the simplified scheme (cl_scheme = 1)
! 2) indexes and molecular mass of chemical species

implicit none

!----------------------------------------------------------------------------
!     chemical tracers
!----------------------------------------------------------------------------

integer, save :: i_co2, i_co, i_h2, i_h2o, i_o1d,        &
                 i_o, i_o2, i_o2dg, i_o3, i_h,           &
                 i_oh, i_ho2, i_h2o2, i_cl, i_clo,       &
                 i_cl2, i_hcl, i_hocl, i_clco, i_clco3,  &
                 i_cocl2, i_s, i_so, i_so2, i_so3,       &
                 i_osso_cis, i_osso_trans, i_s2o2_cyc,   &
                 i_ocs, i_hso3, i_h2so4, i_s2,           &
                 i_clso2, i_cl2so2, i_oscl, i_n2, i_he,  &
                 i_n, i_no, i_no2, i_n2d,                &
                 i_hd, i_hdo, i_d, i_od, i_do2, i_hdo2,  &
                 i_hdso4, i_dcl, i_docl,                 &
                 i_co2plus, i_coplus, i_oplus, i_o2plus, &
                 i_n2plus, i_hplus, i_h2oplus, i_nplus,  &
                 i_ohplus, i_cplus, i_noplus, i_h3oplus, &
                 i_hcoplus, i_hco2plus, i_elec,          &
                 i_dco2plus, i_dcoplus, i_hdoplus,       &
                 i_dplus, i_odplus, i_h2doplus

integer, save :: i_h2oliq, i_hdoliq, i_h2so4liq, i_hdso4liq

integer, save :: i_m0_aer, i_m3_aer,                       &
                 i_m0_mode1drop, i_m0_mode1ccn,            &
                 i_m3_mode1sa, i_m3_mode1w, i_m3_mode1ccn, &
                 i_m0_mode2drop, i_m0_mode2ccn,            &
                 i_m3_mode2sa, i_m3_mode2w, i_m3_mode2ccn

integer, save :: nmicro  ! number of species in the liquid phase

real, dimension(:), save, allocatable :: m_tr           ! molecular mass of tracers
real, dimension(:), save, allocatable :: type_tr        ! type of tracer

real, dimension(:,:), save, allocatable :: no_emission
!$OMP THREADPRIVATE(no_emission)
real, dimension(:,:), save, allocatable :: o2_emission
!$OMP THREADPRIVATE(o2_emission)

!----------------------------------------------------------------------------
!     cloud parameters
!----------------------------------------------------------------------------

integer, save :: cloudmin, cloudmax

!     qrad : ratio radius shell model of mode 3 (cimino, icarus, 1982)
!     if qrad = 0, fully liquid, if qrad = 1 fully solid

real, save :: qrad

!     median radius and standard deviation in each mode

real, save, dimension(:,:,:), allocatable :: r_median, stddev

!     k_mass : defines how the condensed phase is distributed in each mode.
!              sum of k_mass = 1

real, save, dimension(:,:,:), allocatable :: k_mass

real, save, dimension(:,:,:), allocatable :: nbrtot
real, save, dimension(:,:), allocatable :: wh2so4
real, save, dimension(:,:), allocatable :: rho_droplet

contains

!============================================================================

subroutine cloud_ini(nbr_lon, nbr_lev, nbr_mode)

!============================================================================

!     sets cloud microphysical parameters for each mode:
!     radius, standard deviation, mass distribution

integer :: nbr_lon, nbr_lev, nbr_mode
integer :: i_lev, ilon

allocate(nbrtot(nbr_lon,nbr_lev,nbr_mode))
allocate(r_median(nbr_lon,nbr_lev,nbr_mode))
allocate(k_mass(nbr_lon,nbr_lev,nbr_mode))
allocate(stddev(nbr_lon,nbr_lev,nbr_mode))
allocate(wh2so4(nbr_lon,nbr_lev))
allocate(rho_droplet(nbr_lon,nbr_lev))

!     initialisation

r_median(:,:,:)  = 0.     ! median radius
stddev(:,:,:)    = 0.     ! geometric std deviation
k_mass(:,:,:)    = 0.     ! coeff mass multimodal
nbrtot(:,:,:)    = 0.
wh2so4(:,:)      = 0.
rho_droplet(:,:) = 0.

!     minimum and maximum levels for the clouds

cloudmin = 20   ! 20: 38 km
cloudmax = 50   ! 50: 95 km

print*,'================================'
print*,'start initialisation cloud layer'
print*,'================================'

!	===============================================
!	knollenberg & hunten, 1980 and james et al 1997
!	===============================================
!	initialisation unimodale
!	===============================================

!     lower haze: mode 1
!      do i_lev=cloudmin,20
!      r_median(:,i_lev,1)=0.2e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
!      stddev(:,i_lev,1)=1.56
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
!      k_mass(:,i_lev,1)=1.0
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
!      end do

!     lower cloud: mode 3
!      do i_lev=21,23
!      r_median(:,i_lev,1)=3.65e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
!      stddev(:,i_lev,1)=1.28
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
!      k_mass(:,i_lev,1)=1.0
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
!      end do

!     middle cloud: mode 2 prime
!      do i_lev=24,28
!      r_median(:,i_lev,1)=1.4e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
!      stddev(:,i_lev,1)=1.23
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
!      k_mass(:,i_lev,1)=1.0
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
!      end do

!     upper cloud: mode 2
!      do i_lev=29,35
!      r_median(:,i_lev,1)=1.0e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
!      stddev(:,i_lev,1)=1.29
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
!      k_mass(:,i_lev,1)=1.0
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
!      end do

!     upper haze: mode 1
!      do i_lev=36, cloudmax
!      r_median(:,i_lev,1)=0.2e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
!      stddev(:,i_lev,1)=2.16
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
!      k_mass(:,i_lev,1)=1.0
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
!      end do

!	===============================================
!	initialisation trimodale
!	===============================================

!     lower haze: mode 1
!      do i_lev=cloudmin,20
!      r_median(:,i_lev,1)=0.3e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
!      stddev(:,i_lev,1)=1.56
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
!      k_mass(:,i_lev,1)=1.0
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
!      end do

!     lower haze: mode 2
   !   do i_lev=cloudmin,20
   !   r_median(:,i_lev,2)=1.4e-6
   !   print*,'level',i_lev,'r r_median',r_median(1,i_lev,2)
   !   stddev(:,i_lev,2)=1.23
   !   print*,'level',i_lev,'dev std',stddev(1,i_lev,2)
   !   k_mass(:,i_lev,2)=0.0
   !   print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,2)
   !   end do

!     lower haze: mode 3
!      do i_lev=cloudmin,20
!      r_median(:,i_lev,3)=3.65e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,3)
!      stddev(:,i_lev,3)=1.28
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,3)
!      k_mass(:,i_lev,3)=0.
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,3)
!      end do

!     lower cloud: mode 1
!      do i_lev=21,23
!      r_median(:,i_lev,1)=0.3e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
!      stddev(:,i_lev,1)=1.56
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
!      k_mass(:,i_lev,1)=0.1
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
!      end do

!     lower cloud: mode 2 prime
!      do i_lev=21,23
!      r_median(:,i_lev,2)=1.4e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,2)
!      stddev(:,i_lev,2)=1.23
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,2)
!      k_mass(:,i_lev,2)=0.4
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,2)
!      end do

!     lower cloud: mode 3
!      do i_lev=21,23
!      r_median(:,i_lev,3)=3.65e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,3)
!      stddev(:,i_lev,3)=1.28
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,3)
!      k_mass(:,i_lev,3)=0.5
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,3)
!      end do

!     middle cloud: mode 1
!      do i_lev=24,28
!      r_median(:,i_lev,1)=0.3e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
!      stddev(:,i_lev,1)=1.56
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
!      k_mass(:,i_lev,1)=0.0
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
!      end do

!     middle cloud: mode 2 prime
!      do i_lev=24,28
!      r_median(:,i_lev,2)=1.4e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,2)
!      stddev(:,i_lev,2)=1.23
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,2)
!      k_mass(:,i_lev,2)=0.8
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,2)
!      end do

!     middle cloud: mode 3
!      do i_lev=24,28
!      r_median(:,i_lev,3)=3.65e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,3)
!      stddev(:,i_lev,3)=1.28
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,3)
!      k_mass(:,i_lev,3)=0.2
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,3)
!      end do


!     upper cloud: mode 1
!      do i_lev=29,35
!      r_median(:,i_lev,1)=0.3e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
!      stddev(:,i_lev,1)=1.56
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
!      k_mass(:,i_lev,1)=0.15
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
!      end do

!     upper cloud: mode 2
!      do i_lev=29,35
!      r_median(:,i_lev,2)=1.0e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,2)
!      stddev(:,i_lev,2)=1.29
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,2)
!      k_mass(:,i_lev,2)=0.85
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,2)
!      end do

!     upper cloud: mode 3
!      do i_lev=29,35
!      r_median(:,i_lev,3)=3.65e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,3)
!      stddev(:,i_lev,3)=1.28
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,3)
!      k_mass(:,i_lev,3)=0.0
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,3)
!      end do

!     upper haze: mode 1
!      do i_lev=36, cloudmax
!      r_median(:,i_lev,1)=0.3e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
!      stddev(:,i_lev,1)=1.56
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
!      k_mass(:,i_lev,1)=1.0
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
!      end do

!     upper haze: mode 2
!      do i_lev=36, cloudmax
!      r_median(:,i_lev,2)=1.e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,2)
!      stddev(:,i_lev,2)=1.29
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,2)
!      k_mass(:,i_lev,2)=0.0
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,2)
!      end do

!     upper haze: mode 3
!      do i_lev=36, cloudmax
!      r_median(:,i_lev,3)=3.65e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,3)
!      stddev(:,i_lev,3)=2.16
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,3)
!      k_mass(:,i_lev,3)=0.0
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,3)
!      end do
!=============================================================

!	===============================================
!	initialisation trimodale knollenberg
!	===============================================

!     lower haze: mode 1
      ! do i_lev=cloudmin,22
      ! r_median(:,i_lev,1)=0.1e-6
      ! print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
      ! stddev(:,i_lev,1)=1.57
      ! print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
      ! k_mass(:,i_lev,1)=1.0
      ! print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
      ! end do

!     lower haze: mode 2
      ! do i_lev=cloudmin,22
      ! r_median(:,i_lev,2)=1.4e-6
      ! print*,'level',i_lev,'r r_median',r_median(1,i_lev,2)
      ! stddev(:,i_lev,2)=1.23
      ! print*,'level',i_lev,'dev std',stddev(1,i_lev,2)
      ! k_mass(:,i_lev,2)=0.0
      ! print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,2)
      ! end do

!     lower haze: mode 3
      ! do i_lev=cloudmin,22
      ! r_median(:,i_lev,3)=3.65e-6
      ! print*,'level',i_lev,'r r_median',r_median(1,i_lev,3)
      ! stddev(:,i_lev,3)=1.28
      ! print*,'level',i_lev,'dev std',stddev(1,i_lev,3)
      ! k_mass(:,i_lev,3)=0.0
      ! print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,3)
      ! end do

!     pre cloud: mode 1
      ! do i_lev=23,23
      ! r_median(:,i_lev,1)=0.15e-6
      ! print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
      ! stddev(:,i_lev,1)=1.8
      ! print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
      ! k_mass(:,i_lev,1)=0.04
      ! print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
      ! end do

!     pre cloud: mode 2
      ! do i_lev=23,23
      ! r_median(:,i_lev,2)=1.0e-6
      ! print*,'level',i_lev,'r r_median',r_median(1,i_lev,2)
      ! stddev(:,i_lev,2)=1.29
      ! print*,'level',i_lev,'dev std',stddev(1,i_lev,2)
      ! k_mass(:,i_lev,2)=0.96
      ! print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,2)
      ! end do

!     pre cloud: mode 3
      ! do i_lev=23,23
      ! r_median(:,i_lev,3)=3.65e-6
      ! print*,'level',i_lev,'r r_median',r_median(1,i_lev,3)
      ! stddev(:,i_lev,3)=1.28
      ! print*,'level',i_lev,'dev std',stddev(1,i_lev,3)
      ! k_mass(:,i_lev,3)=0.0
      ! print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,3)
      ! end do

!      lower cloud: mode 1
      ! do i_lev=24,24
      ! r_median(:,i_lev,1)=0.2e-6
      ! print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
      ! stddev(:,i_lev,1)=1.8
      ! print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
      ! k_mass(:,i_lev,1)=0.014
      ! print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
      ! end do

!     lower cloud: mode 2
      ! do i_lev=24,24
      ! r_median(:,i_lev,2)=1.0e-6
      ! print*,'level',i_lev,'r r_median',r_median(1,i_lev,2)
      ! stddev(:,i_lev,2)=1.29
      ! print*,'level',i_lev,'dev std',stddev(1,i_lev,2)
      ! k_mass(:,i_lev,2)=0.02
      ! print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,2)
      ! end do

!     lower cloud: mode 3
      ! do i_lev=24,24
      ! r_median(:,i_lev,3)=3.65e-6
      ! print*,'level',i_lev,'r r_median',r_median(1,i_lev,3)
      ! stddev(:,i_lev,3)=1.28
      ! print*,'level',i_lev,'dev std',stddev(1,i_lev,3)
      ! k_mass(:,i_lev,3)=0.966
      ! print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,3)
      ! end do

!     middle cloud: mode 1
      ! do i_lev=25,28
      ! r_median(:,i_lev,1)=0.15e-6
      ! print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
      ! stddev(:,i_lev,1)=1.9
      ! print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
      ! k_mass(:,i_lev,1)=0.0084
      ! print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
      ! end do

!     middle cloud: mode 2 prime
      ! do i_lev=25,28
      ! r_median(:,i_lev,2)=1.4e-6
      ! print*,'level',i_lev,'r r_median',r_median(1,i_lev,2)
      ! stddev(:,i_lev,2)=1.23
      ! print*,'level',i_lev,'dev std',stddev(1,i_lev,2)
      ! k_mass(:,i_lev,2)=0.21
      ! print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,2)
      ! end do

!     middle cloud: mode 3
      ! do i_lev=25,28
      ! r_median(:,i_lev,3)=3.65e-6
      ! print*,'level',i_lev,'r r_median',r_median(1,i_lev,3)
      ! stddev(:,i_lev,3)=1.28
      ! print*,'level',i_lev,'dev std',stddev(1,i_lev,3)
      ! k_mass(:,i_lev,3)=0.7816
      ! print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,3)
      ! end do

!      option: upper haze remplacee par extension upper cloud
!         => 35 remplace par cloudmax et upper haze commentee
!	===============================================

!     upper cloud: mode 1
      ! do i_lev=29,35 !cloudmax
      ! r_median(:,i_lev,1)=0.2e-6
      ! print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
      ! stddev(:,i_lev,1)=2.16
      ! print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
      ! k_mass(:,i_lev,1)=0.72
      ! print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
      ! end do

!     upper cloud: mode 2
      ! do i_lev=29,35 !cloudmax
      ! r_median(:,i_lev,2)=1.0e-6
      ! print*,'level',i_lev,'r r_median',r_median(1,i_lev,2)
      ! stddev(:,i_lev,2)=1.29
      ! print*,'level',i_lev,'dev std',stddev(1,i_lev,2)
      ! k_mass(:,i_lev,2)=0.28
      ! print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,2)
      ! end do

!     upper cloud: mode 3
      ! do i_lev=29,35 !cloudmax
      ! r_median(:,i_lev,3)=3.65e-6
      ! print*,'level',i_lev,'r r_median',r_median(1,i_lev,3)
      ! stddev(:,i_lev,3)=1.28
      ! print*,'level',i_lev,'dev std',stddev(1,i_lev,3)
      ! k_mass(:,i_lev,3)=0.0
      ! print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,3)
      ! end do

!     upper haze: mode 1
      ! do i_lev=36, cloudmax
      ! r_median(:,i_lev,1)=0.2e-6
      ! print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
      ! stddev(:,i_lev,1)=2.16
      ! print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
      ! k_mass(:,i_lev,1)=1.0
      ! print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
      ! end do

!     upper haze: mode 2
      ! do i_lev=36, cloudmax
      ! r_median(:,i_lev,2)=1.e-6
      ! print*,'level',i_lev,'r r_median',r_median(1,i_lev,2)
      ! stddev(:,i_lev,2)=1.29
      ! print*,'level',i_lev,'dev std',stddev(1,i_lev,2)
      ! k_mass(:,i_lev,2)=0.0
      ! print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,2)
      ! end do

!     upper haze: mode 3
      ! do i_lev=36, cloudmax
      ! r_median(:,i_lev,3)=3.65e-6
      ! print*,'level',i_lev,'r r_median',r_median(1,i_lev,3)
      ! stddev(:,i_lev,3)=2.16
      ! print*,'level',i_lev,'dev std',stddev(1,i_lev,3)
      ! k_mass(:,i_lev,3)=0.0
      ! print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,3)
      ! end do

!=============================================================

!	===============================================================
!	initialisation trimodale "knollenberg" sans mode3, mode2 etendu
!	===============================================================

!     lower haze: mode 1
!      do i_lev=cloudmin,22
!      r_median(:,i_lev,1)=0.1e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
!      stddev(:,i_lev,1)=1.57
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
!      k_mass(:,i_lev,1)=1.0
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
!      end do

!     lower haze: mode 2
!      do i_lev=cloudmin,22
!      r_median(:,i_lev,2)=1.4e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,2)
!      stddev(:,i_lev,2)=1.23
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,2)
!      k_mass(:,i_lev,2)=0.0
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,2)
!      end do

!     lower haze: mode 3
!      do i_lev=cloudmin,22
!      r_median(:,i_lev,3)=3.65e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,3)
!      stddev(:,i_lev,3)=1.28
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,3)
!      k_mass(:,i_lev,3)=0.0
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,3)
!      end do

!     pre cloud: mode 1
!      do i_lev=23,23
!      r_median(:,i_lev,1)=0.15e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
!      stddev(:,i_lev,1)=1.8
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
!      k_mass(:,i_lev,1)=0.04
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
!      end do

!     pre cloud: mode 2
!      do i_lev=23,23
!      r_median(:,i_lev,2)=1.0e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,2)
!      stddev(:,i_lev,2)=1.29
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,2)
!      k_mass(:,i_lev,2)=0.96
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,2)
!      end do

!     pre cloud: mode 3
!      do i_lev=23,23
!      r_median(:,i_lev,3)=3.65e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,3)
!      stddev(:,i_lev,3)=1.28
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,3)
!      k_mass(:,i_lev,3)=0.0
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,3)
!      end do

!      lower cloud: mode 1
!      do i_lev=24,24
!      r_median(:,i_lev,1)=0.2e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
!      stddev(:,i_lev,1)=1.8
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
!      k_mass(:,i_lev,1)=0.014
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
!      end do

!     lower cloud: mode 2
!      do i_lev=24,24
!      r_median(:,i_lev,2)=1.0e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,2)
!      stddev(:,i_lev,2)=1.6
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,2)
!      k_mass(:,i_lev,2)=0.986
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,2)
!      end do

!     lower cloud: mode 3
!      do i_lev=24,24
!      r_median(:,i_lev,3)=3.65e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,3)
!      stddev(:,i_lev,3)=1.28
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,3)
!      k_mass(:,i_lev,3)=0.
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,3)
!      end do

!     middle cloud: mode 1
!      do i_lev=25,28
!      r_median(:,i_lev,1)=0.15e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
!      stddev(:,i_lev,1)=1.9
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
!      k_mass(:,i_lev,1)=0.0084
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
!      end do

!     middle cloud: mode 2 prime
!      do i_lev=25,28
!      r_median(:,i_lev,2)=1.4e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,2)
!      stddev(:,i_lev,2)=1.6
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,2)
!      k_mass(:,i_lev,2)=0.9916
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,2)
!      end do

!     middle cloud: mode 3
!      do i_lev=25,28
!      r_median(:,i_lev,3)=3.65e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,3)
!      stddev(:,i_lev,3)=1.28
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,3)
!      k_mass(:,i_lev,3)=0.0
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,3)
!      end do


!     upper cloud: mode 1
!      do i_lev=29,35
!      r_median(:,i_lev,1)=0.2e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
!      stddev(:,i_lev,1)=2.16
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
!      k_mass(:,i_lev,1)=0.72
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
!      end do

!     upper cloud: mode 2
!      do i_lev=29,35
!      r_median(:,i_lev,2)=1.0e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,2)
!      stddev(:,i_lev,2)=1.29
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,2)
!      k_mass(:,i_lev,2)=0.28
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,2)
!      end do

!     upper cloud: mode 3
!      do i_lev=29,35
!      r_median(:,i_lev,3)=3.65e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,3)
!      stddev(:,i_lev,3)=1.28
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,3)
!      k_mass(:,i_lev,3)=0.0
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,3)
!      end do

!     upper haze: mode 1
!      do i_lev=36, cloudmax
!      r_median(:,i_lev,1)=0.2e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
!      stddev(:,i_lev,1)=2.16
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
!      k_mass(:,i_lev,1)=1.0
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
!      end do

!     upper haze: mode 2
!      do i_lev=36, cloudmax
!      r_median(:,i_lev,2)=1.e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,2)
!      stddev(:,i_lev,2)=1.29
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,2)
!      k_mass(:,i_lev,2)=0.0
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,2)
!      end do

!     upper haze: mode 3
!      do i_lev=36, cloudmax
!      r_median(:,i_lev,3)=3.65e-6
!      print*,'level',i_lev,'r r_median',r_median(1,i_lev,3)
!      stddev(:,i_lev,3)=2.16
!      print*,'level',i_lev,'dev std',stddev(1,i_lev,3)
!      k_mass(:,i_lev,3)=0.0
!      print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,3)
!      end do

   ! ========================================================
   ! initialisation bimodale k&h 1980 with mode 3 fully solid
   ! ========================================================
   !    ! mode 3 fully solid
   !    qrad=1
   !    ! normally nb_mode=2 in physiq.def !!!
   !    do ilon=1,nbr_lon
   ! !     mode 1
   !       do i_lev=cloudmin,20
   !          r_median(ilon,i_lev,1)=0.125e-6
   !          print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
   !          stddev(ilon,i_lev,1)=1.57
   !          print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
   !          k_mass(ilon,i_lev,1)=1.0
   !          print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
   !       end do

   !       r_median(ilon,21,1)=0.125e-6
   !       print*,'level',21,'r r_median',r_median(1,21,1)
   !       stddev(ilon,21,1)=1.57
   !       print*,'level',21,'dev std',stddev(1,21,1)
   !       k_mass(ilon,21,1)=0.02
   !       print*,'level',21,'coeff mass: k_mass',k_mass(1,21,1)

   !       r_median(ilon,22,1)=0.2e-6
   !       print*,'level',22,'r r_median',r_median(1,22,1)
   !       stddev(ilon,22,1)=1.8
   !       print*,'level',22,'dev std',stddev(1,22,1)
   !       k_mass(ilon,22,1)=0.02
   !       print*,'level',22,'coeff mass: k_mass',k_mass(1,22,1)

   !       r_median(ilon,23,1)=0.15e-6
   !       print*,'level',23,'r r_median',r_median(1,23,1)
   !       stddev(ilon,23,1)=1.8
   !       print*,'level',23,'dev std',stddev(1,23,1)
   !       k_mass(ilon,23,1)=0.02
   !       print*,'level',23,'coeff mass: k_mass',k_mass(1,23,1)

   !       do i_lev=24,25
   !          r_median(ilon,i_lev,1)=0.15e-6
   !          print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
   !          stddev(ilon,i_lev,1)=1.9
   !          print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
   !          k_mass(ilon,i_lev,1)=0.02
   !          print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
   !       end do

   !       r_median(ilon,26,1)=0.175e-6
   !       print*,'level',26,'r r_median',r_median(1,26,1)
   !       stddev(ilon,26,1)=2.16
   !       print*,'level',26,'dev std',stddev(1,26,1)
   !       k_mass(ilon,26,1)=0.175
   !       print*,'level',26,'coeff mass: k_mass',k_mass(1,26,1)

   !       do i_lev=27,33
   !          r_median(ilon,i_lev,1)=0.175e-6
   !          print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
   !          stddev(ilon,i_lev,1)=2.16
   !          print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
   !          k_mass(ilon,i_lev,1)=0.25
   !          print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
   !       end do

   !       do i_lev=34,cloudmax
   !          r_median(ilon,i_lev,1)=0.175e-6
   !          print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
   !          stddev(ilon,i_lev,1)=2.16
   !          print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
   !          k_mass(ilon,i_lev,1)=1.0
   !          print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
   !       end do

   !    !     mode 2
   !       do i_lev=cloudmin,20
   !          r_median(ilon,i_lev,2)=1.4e-6
   !          print*,'level',i_lev,'r r_median',r_median(1,i_lev,2)
   !          stddev(ilon,i_lev,2)=1.35
   !          print*,'level',i_lev,'dev std',stddev(1,i_lev,2)
   !          k_mass(ilon,i_lev,2)=0.0
   !          print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,2)
   !       end do

   !       do i_lev=21,22
   !          r_median(ilon,i_lev,2)=1.4e-6
   !          print*,'level',i_lev,'r r_median',r_median(1,i_lev,2)
   !          stddev(ilon,i_lev,2)=1.35
   !          print*,'level',i_lev,'dev std',stddev(1,i_lev,2)
   !          k_mass(ilon,i_lev,2)=0.98
   !          print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,2)
   !       end do

   !       r_median(ilon,23,2)=1.35e-6
   !       print*,'level',23,'r r_median',r_median(1,23,2)
   !       stddev(ilon,23,2)=1.25
   !       print*,'level',23,'dev std',stddev(1,23,2)
   !       k_mass(ilon,23,2)=0.98
   !       print*,'level',23,'coeff mass: k_mass',k_mass(1,23,2)


   !       r_median(ilon,24,2)=1.375e-6
   !       print*,'level',24,'r r_median',r_median(1,24,2)
   !       stddev(ilon,24,2)=1.2
   !       print*,'level',24,'dev std',stddev(1,24,2)
   !       k_mass(ilon,24,2)=0.98
   !       print*,'level',24,'coeff mass: k_mass',k_mass(1,24,2)


   !       r_median(ilon,25,2)=1.4e-6
   !       print*,'level',25,'r r_median',r_median(1,25,2)
   !       stddev(ilon,25,2)=1.16
   !       print*,'level',25,'dev std',stddev(1,25,2)
   !       k_mass(ilon,25,2)=0.98
   !       print*,'level',25,'coeff mass: k_mass',k_mass(1,25,2)


   !       r_median(ilon,26,2)=1.15e-6
   !       print*,'level',26,'r r_median',r_median(1,26,2)
   !       stddev(ilon,26,2)=1.34
   !       print*,'level',26,'dev std',stddev(1,26,2)
   !       k_mass(ilon,26,2)=0.825
   !       print*,'level',26,'coeff mass: k_mass',k_mass(1,26,2)


   !       r_median(ilon,27,2)=1.14e-6
   !       print*,'level',27,'r r_median',r_median(1,27,2)
   !       stddev(ilon,27,2)=1.33
   !       print*,'level',27,'dev std',stddev(1,27,2)
   !       k_mass(ilon,27,2)=0.75
   !       print*,'level',27,'coeff mass: k_mass',k_mass(1,27,2)


   !       r_median(ilon,28,2)=1.35e-6
   !       print*,'level',28,'r r_median',r_median(1,28,2)
   !       stddev(ilon,28,2)=1.32
   !       print*,'level',28,'dev std',stddev(1,28,2)
   !       k_mass(ilon,28,2)=0.75
   !       print*,'level',28,'coeff mass: k_mass',k_mass(1,28,2)


   !       r_median(ilon,29,2)=1.125e-6
   !       print*,'level',29,'r r_median',r_median(1,29,2)
   !       stddev(ilon,29,2)=1.31
   !       print*,'level',29,'dev std',stddev(1,29,2)
   !       k_mass(ilon,29,2)=0.75
   !       print*,'level',29,'coeff mass: k_mass',k_mass(1,29,2)


   !       r_median(ilon,30,2)=1.118e-6
   !       print*,'level',30,'r r_median',r_median(1,30,2)
   !       stddev(ilon,30,2)=1.30
   !       print*,'level',30,'dev std',stddev(1,30,2)
   !       k_mass(ilon,30,2)=0.75
   !       print*,'level',30,'coeff mass: k_mass',k_mass(1,30,2)


   !       r_median(ilon,31,2)=1.11e-6
   !       print*,'level',31,'r r_median',r_median(1,31,2)
   !       stddev(ilon,31,2)=1.29
   !       print*,'level',31,'dev std',stddev(1,31,2)
   !       k_mass(ilon,31,2)=0.75
   !       print*,'level',31,'coeff mass: k_mass',k_mass(1,31,2)

   !       do i_lev=32,33
   !          r_median(ilon,i_lev,2)=1.1e-6
   !          print*,'level',i_lev,'r r_median',r_median(1,i_lev,2)
   !          stddev(ilon,i_lev,2)=1.28
   !          print*,'level',i_lev,'dev std',stddev(1,i_lev,2)
   !          k_mass(ilon,i_lev,2)=0.75
   !          print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,2)
   !       end do

   !       do i_lev=34,cloudmax
   !          r_median(ilon,i_lev,2)=1.1e-6
   !          print*,'level',i_lev,'r r_median',r_median(1,i_lev,2)
   !          stddev(ilon,i_lev,2)=1.28
   !          print*,'level',i_lev,'dev std',stddev(1,i_lev,2)
   !          k_mass(ilon,i_lev,2)=0.0
   !          print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,2)
   !       end do
   !    end do

! ======================================================
! knollenberg & hunten (1980) with mode 3 97% solid
! ======================================================

! mode 3 97% solid => qrad=0.97

qrad = 0.97

do ilon = 1,nbr_lon

   ! mode 1

   do i_lev = cloudmin,21
      r_median(ilon,i_lev,1) = 0.125e-6
      stddev(ilon,i_lev,1)   = 1.57
      k_mass(ilon,i_lev,1)   = 1.0
   end do

   r_median(ilon,22,1)       = 0.125e-6
   stddev(ilon,22,1)         = 1.57
   k_mass(ilon,22,1)         = 1.0

   r_median(ilon,23,1)       = 0.2e-6
   stddev(ilon,23,1)         = 1.8
   k_mass(ilon,23,1)         = 0.01

   r_median(ilon,24,1)       = 0.15e-6
   stddev(ilon,24,1)         = 1.8
   k_mass(ilon,24,1)         = 0.01

   do i_lev = 25,26
      r_median(ilon,i_lev,1) = 0.15e-6
      stddev(ilon,i_lev,1)   = 1.9
      k_mass(ilon,i_lev,1)   = 0.01
   end do

   r_median(ilon,27,1)       = 0.175e-6
   stddev(ilon,27,1)         = 2.16
   k_mass(ilon,27,1)         = 0.175

   do i_lev = 28,34
      r_median(ilon,i_lev,1) = 0.175e-6
      stddev(ilon,i_lev,1)   = 2.16
      k_mass(ilon,i_lev,1)   = 0.25
   end do

   do i_lev = 35,cloudmax
      r_median(ilon,i_lev,1) = 0.175e-6
      stddev(ilon,i_lev,1)   = 2.16
      k_mass(ilon,i_lev,1)   = 0.25
   end do

   !  mode 2

   do i_lev = cloudmin,21
      r_median(ilon,i_lev,2) = 1.4e-6
      stddev(ilon,i_lev,2)   = 1.35
      k_mass(ilon,i_lev,2)   = 0.0
   end do

   r_median(ilon,22,2)       = 1.4e-6
   stddev(ilon,22,2)         = 1.25
   k_mass(ilon,22,2)         = 0.0

   r_median(ilon,23,2)       = 1.4e-6
   stddev(ilon,23,2)         = 1.25
   k_mass(ilon,23,2)         = 0.35

   r_median(ilon,24,2)       = 1.35e-6
   stddev(ilon,24,2)         = 1.25
   k_mass(ilon,24,2)         = 0.35

   r_median(ilon,25,2)       = 1.375e-6
   stddev(ilon,25,2)         = 1.2
   k_mass(ilon,25,2)         = 0.35

   r_median(ilon,26,2)       = 1.4e-6
   stddev(ilon,26,2)         = 1.16
   k_mass(ilon,26,2)         = 0.35

   r_median(ilon,27,2)       = 1.15e-6
   stddev(ilon,27,2)         = 1.34
   k_mass(ilon,27,2)         = 0.825

   r_median(ilon,28,2)       = 1.14e-6
   stddev(ilon,28,2)         = 1.33
   k_mass(ilon,28,2)         = 0.75

   r_median(ilon,29,2)       = 1.35e-6
   stddev(ilon,29,2)         = 1.32
   k_mass(ilon,29,2)         = 0.75

   r_median(ilon,30,2)       = 1.125e-6
   stddev(ilon,30,2)         = 1.31
   k_mass(ilon,30,2)         = 0.75

   r_median(ilon,31,2)       = 1.118e-6
   stddev(ilon,31,2)         = 1.31
   k_mass(ilon,31,2)         = 0.75

   r_median(ilon,32,2)       = 1.11e-6
   stddev(ilon,32,2)         = 1.29
   k_mass(ilon,32,2)         = 0.75

   do i_lev = 33,34
      r_median(ilon,i_lev,2) = 1.1e-6
      stddev(ilon,i_lev,2)   = 1.28
      k_mass(ilon,i_lev,2)   = 0.75
   end do

   do i_lev = 35,cloudmax
      r_median(ilon,i_lev,2) = 1.1e-6
      stddev(ilon,i_lev,2)   = 1.28
      k_mass(ilon,i_lev,2)   = 0.75
   end do

   ! mode 3

   do i_lev = cloudmin,22
      r_median(ilon,i_lev,3) = 3.65e-6
      stddev(ilon,i_lev,3)   = 1.28
      k_mass(ilon,i_lev,3)   = 0.0
   end do

   do i_lev = 23,26
      r_median(ilon,i_lev,3) = 3.65e-6
      stddev(ilon,i_lev,3)   = 1.28
      k_mass(ilon,i_lev,3)   = 0.64
   end do

   do i_lev = 27,cloudmax
      r_median(ilon,i_lev,3) = 3.65e-6
      stddev(ilon,i_lev,3)   = 1.28
      k_mass(ilon,i_lev,3)   = 0.0
   end do
end do

      ! ! ==================================================================
      ! ! initialisation bimodale k&h 1980 with mode 3 0% solid fully liquid
      ! ! ==================================================================
      ! ! mode 3 0% solid, fully liquid
      ! qrad=0.0
      ! ! normally nb_mode=3 in physiq.def !!!
      ! do ilon=1,nbr_lon
      !    ! mode 1
      !    do i_lev=cloudmin,20
      !       r_median(ilon,i_lev,1)=0.125e-6
      !       print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
      !       stddev(ilon,i_lev,1)=1.57
      !       print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
      !       k_mass(ilon,i_lev,1)=1.0
      !       print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
      !    end do

      !    r_median(ilon,21,1)=0.125e-6
      !    print*,'level',21,'r r_median',r_median(1,21,1)
      !    stddev(ilon,21,1)=1.57
      !    print*,'level',21,'dev std',stddev(1,21,1)
      !    k_mass(ilon,21,1)=1.0
      !    print*,'level',21,'coeff mass: k_mass',k_mass(1,21,1)

      !    r_median(ilon,22,1)=0.2e-6
      !    print*,'level',22,'r r_median',r_median(1,22,1)
      !    stddev(ilon,22,1)=1.8
      !    print*,'level',22,'dev std',stddev(1,22,1)
      !    k_mass(ilon,22,1)=0.01
      !    print*,'level',22,'coeff mass: k_mass',k_mass(1,22,1)

      !    r_median(ilon,23,1)=0.15e-6
      !    print*,'level',23,'r r_median',r_median(1,23,1)
      !    stddev(ilon,23,1)=1.8
      !    print*,'level',23,'dev std',stddev(1,23,1)
      !    k_mass(ilon,23,1)=0.01
      !    print*,'level',23,'coeff mass: k_mass',k_mass(1,23,1)

      !    do i_lev=24,25
      !       r_median(ilon,i_lev,1)=0.15e-6
      !       print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
      !       stddev(ilon,i_lev,1)=1.9
      !       print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
      !       k_mass(ilon,i_lev,1)=0.01
      !       print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
      !    end do

      !    r_median(ilon,26,1)=0.175e-6
      !    print*,'level',26,'r r_median',r_median(1,26,1)
      !    stddev(ilon,26,1)=2.16
      !    print*,'level',26,'dev std',stddev(1,26,1)
      !    k_mass(ilon,26,1)=0.175
      !    print*,'level',26,'coeff mass: k_mass',k_mass(1,26,1)

      !    do i_lev=27,33
      !       r_median(ilon,i_lev,1)=0.175e-6
      !       print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
      !       stddev(ilon,i_lev,1)=2.16
      !       print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
      !       k_mass(ilon,i_lev,1)=0.25
      !       print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
      !    end do

      !    do i_lev=34,cloudmax
      !       r_median(ilon,i_lev,1)=0.175e-6
      !       print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
      !       stddev(ilon,i_lev,1)=2.16
      !       print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
      !       k_mass(ilon,i_lev,1)=0.25
      !       print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
      !    end do

      !    !  mode 2
      !    do i_lev=cloudmin,20
      !       r_median(ilon,i_lev,2)=1.4e-6
      !       print*,'level',i_lev,'r r_median',r_median(1,i_lev,2)
      !       stddev(ilon,i_lev,2)=1.35
      !       print*,'level',i_lev,'dev std',stddev(1,i_lev,2)
      !       k_mass(ilon,i_lev,2)=0.0
      !       print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,2)
      !    end do

      !    r_median(ilon,21,2)=1.4e-6
      !    print*,'level',21,'r r_median',r_median(1,21,2)
      !    stddev(ilon,21,2)=1.25
      !    print*,'level',21,'dev std',stddev(1,21,2)
      !    k_mass(ilon,21,2)=0.0
      !    print*,'level',21,'coeff mass: k_mass',k_mass(1,21,2)

      !    r_median(ilon,22,2)=1.4e-6
      !    print*,'level',22,'r r_median',r_median(1,22,2)
      !    stddev(ilon,22,2)=1.25
      !    print*,'level',22,'dev std',stddev(1,22,2)
      !    k_mass(ilon,22,2)=0.04
      !    print*,'level',22,'coeff mass: k_mass',k_mass(1,22,2)

      !    r_median(ilon,23,2)=1.35e-6
      !    print*,'level',23,'r r_median',r_median(1,23,2)
      !    stddev(ilon,23,2)=1.25
      !    print*,'level',23,'dev std',stddev(1,23,2)
      !    k_mass(ilon,23,2)=0.04
      !    print*,'level',23,'coeff mass: k_mass',k_mass(1,23,2)


      !    r_median(ilon,24,2)=1.375e-6
      !    print*,'level',24,'r r_median',r_median(1,24,2)
      !    stddev(ilon,24,2)=1.2
      !    print*,'level',24,'dev std',stddev(1,24,2)
      !    k_mass(ilon,24,2)=0.04
      !    print*,'level',24,'coeff mass: k_mass',k_mass(1,24,2)


      !    r_median(ilon,25,2)=1.4e-6
      !    print*,'level',25,'r r_median',r_median(1,25,2)
      !    stddev(ilon,25,2)=1.16
      !    print*,'level',25,'dev std',stddev(1,25,2)
      !    k_mass(ilon,25,2)=0.04
      !    print*,'level',25,'coeff mass: k_mass',k_mass(1,25,2)


      !    r_median(ilon,26,2)=1.15e-6
      !    print*,'level',26,'r r_median',r_median(1,26,2)
      !    stddev(ilon,26,2)=1.34
      !    print*,'level',26,'dev std',stddev(1,26,2)
      !    k_mass(ilon,26,2)=0.825
      !    print*,'level',26,'coeff mass: k_mass',k_mass(1,26,2)


      !    r_median(ilon,27,2)=1.14e-6
      !    print*,'level',27,'r r_median',r_median(1,27,2)
      !    stddev(ilon,27,2)=1.33
      !    print*,'level',27,'dev std',stddev(1,27,2)
      !    k_mass(ilon,27,2)=0.75
      !    print*,'level',27,'coeff mass: k_mass',k_mass(1,27,2)


      !    r_median(ilon,28,2)=1.35e-6
      !    print*,'level',28,'r r_median',r_median(1,28,2)
      !    stddev(ilon,28,2)=1.32
      !    print*,'level',28,'dev std',stddev(1,28,2)
      !    k_mass(ilon,28,2)=0.75
      !    print*,'level',28,'coeff mass: k_mass',k_mass(1,28,2)


      !    r_median(ilon,29,2)=1.125e-6
      !    print*,'level',29,'r r_median',r_median(1,29,2)
      !    stddev(ilon,29,2)=1.31
      !    print*,'level',29,'dev std',stddev(1,29,2)
      !    k_mass(ilon,29,2)=0.75
      !    print*,'level',29,'coeff mass: k_mass',k_mass(1,29,2)


      !    r_median(ilon,30,2)=1.118e-6
      !    print*,'level',30,'r r_median',r_median(1,30,2)
      !    stddev(ilon,30,2)=1.30
      !    print*,'level',30,'dev std',stddev(1,30,2)
      !    k_mass(ilon,30,2)=0.75
      !    print*,'level',30,'coeff mass: k_mass',k_mass(1,30,2)


      !    r_median(ilon,31,2)=1.11e-6
      !    print*,'level',31,'r r_median',r_median(1,31,2)
      !    stddev(ilon,31,2)=1.29
      !    print*,'level',31,'dev std',stddev(1,31,2)
      !    k_mass(ilon,31,2)=0.75
      !    print*,'level',31,'coeff mass: k_mass',k_mass(1,31,2)

      !    do i_lev=32,33
      !       r_median(ilon,i_lev,2)=1.1e-6
      !       print*,'level',i_lev,'r r_median',r_median(1,i_lev,2)
      !       stddev(ilon,i_lev,2)=1.28
      !       print*,'level',i_lev,'dev std',stddev(1,i_lev,2)
      !       k_mass(ilon,i_lev,2)=0.75
      !       print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,2)
      !    end do

      !    ! if k_mass > 0 it means we have a bimodal upper haze.
      !    do i_lev=34,cloudmax
      !       r_median(ilon,i_lev,2)=1.1e-6
      !       print*,'level',i_lev,'r r_median',r_median(1,i_lev,2)
      !       stddev(ilon,i_lev,2)=1.28
      !       print*,'level',i_lev,'dev std',stddev(1,i_lev,2)
      !       k_mass(ilon,i_lev,2)=0.75
      !       print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,2)
      !    end do

      !    ! mode 3
      !    do i_lev=cloudmin,21
      !       r_median(ilon,i_lev,3)=3.65e-6
      !       print*,'level',i_lev,'r r_median',r_median(1,i_lev,3)
      !       stddev(ilon,i_lev,3)=1.28
      !       print*,'level',i_lev,'dev std',stddev(1,i_lev,3)
      !       k_mass(ilon,i_lev,3)=0.0
      !       print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,3)
      !    end do

      !    do i_lev=22,25
      !       r_median(ilon,i_lev,3)=3.65e-6
      !       print*,'level',i_lev,'r r_median',r_median(1,i_lev,3)
      !       stddev(ilon,i_lev,3)=1.28
      !       print*,'level',i_lev,'dev std',stddev(1,i_lev,3)
      !       k_mass(ilon,i_lev,3)=0.95
      !       print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,3)
      !    end do

      !    do i_lev=26,cloudmax
      !       r_median(ilon,i_lev,3)=3.65e-6
      !       print*,'level',i_lev,'r r_median',r_median(1,i_lev,3)
      !       stddev(ilon,i_lev,3)=1.28
      !       print*,'level',i_lev,'dev std',stddev(1,i_lev,3)
      !       k_mass(ilon,i_lev,3)=0.0
      !       print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,3)
      !    end do
      ! end do

      ! print*,'r_median size',size(r_median)

      ! print*,'r_median(1,28,3)',r_median(1,28,3)
      ! print*,'r_median(1,28,3) ne devrait pas exister puisque les dim sont (285,78,2)'
      ! print*,'r_median(1,28,3) a la valeur de k_mass(1,28,1)', k_mass(1,28,1)
      ! print*,'r_median(1,28,1) just to see',r_median(1,28,1)
      ! print*,'r_median size after adding some shenanigans',size(r_median)

      ! ! ==================================================================
      ! ! initialisation stupid with mode 3 100 um
      ! ! ==================================================================
      ! ! mode 3 0% solid, fully liquid
      ! qrad=0.0
      ! ! normally nb_mode=3 in physiq.def !!!
      ! do ilon=1,nbr_lon
      !    ! mode 1
      !    do i_lev=cloudmin,20
      !       r_median(ilon,i_lev,1)=0.125e-6
      !       print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
      !       stddev(ilon,i_lev,1)=1.57
      !       print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
      !       k_mass(ilon,i_lev,1)=1.0
      !       print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
      !    end do

      !    r_median(ilon,21,1)=0.125e-6
      !    print*,'level',21,'r r_median',r_median(1,21,1)
      !    stddev(ilon,21,1)=1.57
      !    print*,'level',21,'dev std',stddev(1,21,1)
      !    k_mass(ilon,21,1)=1.0
      !    print*,'level',21,'coeff mass: k_mass',k_mass(1,21,1)

      !    r_median(ilon,22,1)=0.2e-6
      !    print*,'level',22,'r r_median',r_median(1,22,1)
      !    stddev(ilon,22,1)=1.8
      !    print*,'level',22,'dev std',stddev(1,22,1)
      !    k_mass(ilon,22,1)=0.0
      !    print*,'level',22,'coeff mass: k_mass',k_mass(1,22,1)

      !    r_median(ilon,23,1)=0.15e-6
      !    print*,'level',23,'r r_median',r_median(1,23,1)
      !    stddev(ilon,23,1)=1.8
      !    print*,'level',23,'dev std',stddev(1,23,1)
      !    k_mass(ilon,23,1)=0.0
      !    print*,'level',23,'coeff mass: k_mass',k_mass(1,23,1)

      !    do i_lev=24,25
      !       r_median(ilon,i_lev,1)=0.15e-6
      !       print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
      !       stddev(ilon,i_lev,1)=1.9
      !       print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
      !       k_mass(ilon,i_lev,1)=0.0
      !       print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
      !    end do

      !    r_median(ilon,26,1)=0.175e-6
      !    print*,'level',26,'r r_median',r_median(1,26,1)
      !    stddev(ilon,26,1)=2.16
      !    print*,'level',26,'dev std',stddev(1,26,1)
      !    k_mass(ilon,26,1)=0.0
      !    print*,'level',26,'coeff mass: k_mass',k_mass(1,26,1)

      !    do i_lev=27,33
      !       r_median(ilon,i_lev,1)=0.175e-6
      !       print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
      !       stddev(ilon,i_lev,1)=2.16
      !       print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
      !       k_mass(ilon,i_lev,1)=0.0
      !       print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
      !    end do

      !    do i_lev=34,cloudmax
      !       r_median(ilon,i_lev,1)=0.175e-6
      !       print*,'level',i_lev,'r r_median',r_median(1,i_lev,1)
      !       stddev(ilon,i_lev,1)=2.16
      !       print*,'level',i_lev,'dev std',stddev(1,i_lev,1)
      !       k_mass(ilon,i_lev,1)=0.0
      !       print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,1)
      !    end do

      !    !  mode 2
      !    do i_lev=cloudmin,20
      !       r_median(ilon,i_lev,2)=1.4e-6
      !       print*,'level',i_lev,'r r_median',r_median(1,i_lev,2)
      !       stddev(ilon,i_lev,2)=1.35
      !       print*,'level',i_lev,'dev std',stddev(1,i_lev,2)
      !       k_mass(ilon,i_lev,2)=0.0
      !       print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,2)
      !    end do

      !    r_median(ilon,21,2)=1.4e-6
      !    print*,'level',21,'r r_median',r_median(1,21,2)
      !    stddev(ilon,21,2)=1.25
      !    print*,'level',21,'dev std',stddev(1,21,2)
      !    k_mass(ilon,21,2)=0.0
      !    print*,'level',21,'coeff mass: k_mass',k_mass(1,21,2)

      !    r_median(ilon,22,2)=1.4e-6
      !    print*,'level',22,'r r_median',r_median(1,22,2)
      !    stddev(ilon,22,2)=1.25
      !    print*,'level',22,'dev std',stddev(1,22,2)
      !    k_mass(ilon,22,2)=0.0
      !    print*,'level',22,'coeff mass: k_mass',k_mass(1,22,2)

      !    r_median(ilon,23,2)=1.35e-6
      !    print*,'level',23,'r r_median',r_median(1,23,2)
      !    stddev(ilon,23,2)=1.25
      !    print*,'level',23,'dev std',stddev(1,23,2)
      !    k_mass(ilon,23,2)=0.0
      !    print*,'level',23,'coeff mass: k_mass',k_mass(1,23,2)


      !    r_median(ilon,24,2)=1.375e-6
      !    print*,'level',24,'r r_median',r_median(1,24,2)
      !    stddev(ilon,24,2)=1.2
      !    print*,'level',24,'dev std',stddev(1,24,2)
      !    k_mass(ilon,24,2)=0.0
      !    print*,'level',24,'coeff mass: k_mass',k_mass(1,24,2)


      !    r_median(ilon,25,2)=1.4e-6
      !    print*,'level',25,'r r_median',r_median(1,25,2)
      !    stddev(ilon,25,2)=1.16
      !    print*,'level',25,'dev std',stddev(1,25,2)
      !    k_mass(ilon,25,2)=0.0
      !    print*,'level',25,'coeff mass: k_mass',k_mass(1,25,2)


      !    r_median(ilon,26,2)=1.15e-6
      !    print*,'level',26,'r r_median',r_median(1,26,2)
      !    stddev(ilon,26,2)=1.34
      !    print*,'level',26,'dev std',stddev(1,26,2)
      !    k_mass(ilon,26,2)=0.0
      !    print*,'level',26,'coeff mass: k_mass',k_mass(1,26,2)


      !    r_median(ilon,27,2)=1.14e-6
      !    print*,'level',27,'r r_median',r_median(1,27,2)
      !    stddev(ilon,27,2)=1.33
      !    print*,'level',27,'dev std',stddev(1,27,2)
      !    k_mass(ilon,27,2)=0.0
      !    print*,'level',27,'coeff mass: k_mass',k_mass(1,27,2)


      !    r_median(ilon,28,2)=1.35e-6
      !    print*,'level',28,'r r_median',r_median(1,28,2)
      !    stddev(ilon,28,2)=1.32
      !    print*,'level',28,'dev std',stddev(1,28,2)
      !    k_mass(ilon,28,2)=0.0
      !    print*,'level',28,'coeff mass: k_mass',k_mass(1,28,2)


      !    r_median(ilon,29,2)=1.125e-6
      !    print*,'level',29,'r r_median',r_median(1,29,2)
      !    stddev(ilon,29,2)=1.31
      !    print*,'level',29,'dev std',stddev(1,29,2)
      !    k_mass(ilon,29,2)=0.0
      !    print*,'level',29,'coeff mass: k_mass',k_mass(1,29,2)


      !    r_median(ilon,30,2)=1.118e-6
      !    print*,'level',30,'r r_median',r_median(1,30,2)
      !    stddev(ilon,30,2)=1.30
      !    print*,'level',30,'dev std',stddev(1,30,2)
      !    k_mass(ilon,30,2)=0.0
      !    print*,'level',30,'coeff mass: k_mass',k_mass(1,30,2)


      !    r_median(ilon,31,2)=1.11e-6
      !    print*,'level',31,'r r_median',r_median(1,31,2)
      !    stddev(ilon,31,2)=1.29
      !    print*,'level',31,'dev std',stddev(1,31,2)
      !    k_mass(ilon,31,2)=0.0
      !    print*,'level',31,'coeff mass: k_mass',k_mass(1,31,2)

      !    do i_lev=32,33
      !       r_median(ilon,i_lev,2)=1.1e-6
      !       print*,'level',i_lev,'r r_median',r_median(1,i_lev,2)
      !       stddev(ilon,i_lev,2)=1.28
      !       print*,'level',i_lev,'dev std',stddev(1,i_lev,2)
      !       k_mass(ilon,i_lev,2)=0.0
      !       print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,2)
      !    end do

      !    ! if k_mass > 0 it means we have a bimodal upper haze.
      !    do i_lev=34,cloudmax
      !       r_median(ilon,i_lev,2)=1.1e-6
      !       print*,'level',i_lev,'r r_median',r_median(1,i_lev,2)
      !       stddev(ilon,i_lev,2)=1.28
      !       print*,'level',i_lev,'dev std',stddev(1,i_lev,2)
      !       k_mass(ilon,i_lev,2)=0.0
      !       print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,2)
      !    end do

      !    ! mode 3
      !    do i_lev=cloudmin,21
      !       r_median(ilon,i_lev,3)=100.0e-6
      !       print*,'level',i_lev,'r r_median',r_median(1,i_lev,3)
      !       stddev(ilon,i_lev,3)=1.28
      !       print*,'level',i_lev,'dev std',stddev(1,i_lev,3)
      !       k_mass(ilon,i_lev,3)=1.0
      !       print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,3)
      !    end do

      !    do i_lev=22,25
      !       r_median(ilon,i_lev,3)=100.0e-6
      !       print*,'level',i_lev,'r r_median',r_median(1,i_lev,3)
      !       stddev(ilon,i_lev,3)=1.28
      !       print*,'level',i_lev,'dev std',stddev(1,i_lev,3)
      !       k_mass(ilon,i_lev,3)=1.0
      !       print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,3)
      !    end do

      !    do i_lev=26,cloudmax
      !       r_median(ilon,i_lev,3)=100.0e-6
      !       print*,'level',i_lev,'r r_median',r_median(1,i_lev,3)
      !       stddev(ilon,i_lev,3)=1.28
      !       print*,'level',i_lev,'dev std',stddev(1,i_lev,3)
      !       k_mass(ilon,i_lev,3)=1.0
      !       print*,'level',i_lev,'coeff mass: k_mass',k_mass(1,i_lev,3)
      !    end do
      ! end do

      ! print*,'r_median size',size(r_median)

      ! print*,'r_median(1,28,3)',r_median(1,28,3)
      ! print*,'r_median(1,28,3) ne devrait pas exister puisque les dim sont (285,78,2)'
      ! print*,'r_median(1,28,3) a la valeur de k_mass(1,28,1)', k_mass(1,28,1)
      ! print*,'r_median(1,28,1) just to see',r_median(1,28,1)
      ! print*,'r_median size after adding some shenanigans',size(r_median)

! check if sum of kmass = 1 at each level

print*, 'chemparam_mod: start checking k_mass=1'

do i_lev = cloudmin,cloudmax
   if (nbr_mode == 3) then
      if ((k_mass(1,i_lev,1) + k_mass(1,i_lev,2) + k_mass(1,i_lev,3)) /= 1.0) then
         print*, 'kmass total is not 1.0'
         print*, k_mass(1,i_lev,1) + k_mass(1,i_lev,2) + k_mass(1,i_lev,3)
         print*, 'at level: ',i_lev
         print*, 'check chemparam_mod cloud structure definition'
         stop
      end if
   else
      if ((k_mass(1,i_lev,1) + k_mass(1,i_lev,2)) /= 1.0) then
         print*, 'kmass total is not 1.0'
         print*, k_mass(1,i_lev,1) + k_mass(1,i_lev,2)
         print*, 'at level: ',i_lev
         print*, 'check chemparam_mod cloud structure definition'
      end if
   end if
end do

print*, 'chemparam_mod: end checking k_mass=1'
print*, 'k_mass is fine'
print*
print*,'==============================='
print*,'end initialisation cloud layer'
print*,'==============================='

end subroutine cloud_ini

!============================================================================

  subroutine chemparam_ini

!============================================================================

      use infotrac_phy, only: nqtot, tname

      implicit none

      integer :: i

      allocate(m_tr(nqtot))    ! molecular mass of tracers
      allocate(type_tr(nqtot)) ! type of chemical tracers 1: neutral, 2: ion, 3: liquid

      ! index for chemical tracers

      ! neutrals

      i_co2      = 0
      i_co       = 0
      i_h2       = 0
      i_h2o      = 0
      i_o1d      = 0
      i_o        = 0
      i_o2       = 0
      i_o2dg     = 0
      i_o3       = 0
      i_h        = 0
      i_oh       = 0
      i_ho2      = 0
      i_h2o2     = 0
      i_cl       = 0
      i_clo      = 0
      i_cl2      = 0
      i_hcl      = 0
      i_hocl     = 0
      i_clco     = 0
      i_clco3    = 0
      i_cocl2    = 0
      i_s        = 0
      i_so       = 0
      i_so2      = 0
      i_so3      = 0
      i_osso_cis = 0
      i_osso_trans = 0
      i_s2o2_cyc = 0
      i_ocs      = 0
      i_hso3     = 0
      i_h2so4    = 0
      i_s2       = 0
      i_clso2    = 0
      i_cl2so2   = 0
      i_oscl     = 0
      i_n2       = 0
      i_he       = 0
      i_n2d      = 0
      i_n        = 0
      i_no       = 0
      i_no2      = 0
      i_hd        = 0
      i_d         = 0
      i_od        = 0
      i_do2       = 0
      i_hdo2      = 0
      i_hdo       = 0
      i_hdso4     = 0
      i_dcl       = 0
      i_docl      = 0

      ! ions

      i_co2plus  = 0
      i_coplus   = 0
      i_oplus    = 0
      i_o2plus   = 0
      i_n2plus   = 0
      i_hplus    = 0
      i_h2oplus  = 0
      i_nplus    = 0
      i_ohplus   = 0
      i_cplus    = 0
      i_noplus   = 0
      i_h3oplus  = 0
      i_hcoplus  = 0
      i_hco2plus = 0
      i_elec     = 0
      i_dco2plus  = 0
      i_dcoplus   = 0
      i_hdoplus   = 0
      i_dplus     = 0
      i_odplus    = 0
      i_h2doplus  = 0

      ! liquid

      i_h2oliq   = 0
      i_h2so4liq = 0
      i_hdoliq   = 0
      i_hdso4liq = 0

      do i = 1,nqtot
         print*,'tname(i)',tname(i)
         select case(tname(i))

            ! neutrals

            case('co2')
               i_co2 = i
               print*,'co2',i_co2
               m_tr(i_co2) = 44.0095
               type_tr(i_co2) = 1
            case('co')
               i_co = i
               print*,'co',i_co
               m_tr(i_co) = 28.0101
               type_tr(i_co) = 1
            case('h2')
               i_h2 = i
               print*,'h2',i_h2
               m_tr(i_h2) = 2.01588
               type_tr(i_h2) = 1
            case('h2o')
               i_h2o = i
               print*,'h2o',i_h2o
               m_tr(i_h2o) = 18.0153
               type_tr(i_h2o) = 1
            case('o1d')
               i_o1d = i
               print*,'o1d',i_o1d
               m_tr(i_o1d) = 15.994
               type_tr(i_o1d) = 1
            case('o')
               i_o = i
               print*,'o',i_o
               m_tr(i_o) = 15.994
               type_tr(i_o) = 1
            case('o2')
               i_o2 = i
               print*,'o2',i_o2
               m_tr(i_o2) = 31.9988
               type_tr(i_o2) = 1
            case('o2dg')
               i_o2dg = i
               print*,'o2dg',i_o2dg
               m_tr(i_o2dg) = 31.9988
               type_tr(i_o2dg) = 1
            case('o3')
               i_o3 = i
               print*,'o3',i_o3
               m_tr(i_o3) = 47.9982
               type_tr(i_o3) = 1
            case('h')
               i_h = i
               print*,'h',i_h
               m_tr(i_h) = 1.00794
               type_tr(i_h) = 1
            case('oh')
               i_oh = i
               print*,'oh',i_oh
               m_tr(i_oh) = 17.0073
               type_tr(i_oh) = 1
            case('ho2')
               i_ho2 = i
               print*,'ho2',i_ho2
               m_tr(i_ho2) = 33.0067
               type_tr(i_ho2) = 1
            case('h2o2')
               i_h2o2 = i
               print*,'h2o2',i_h2o2
               m_tr(i_h2o2) = 34.0147
               type_tr(i_h2o2) = 1
            case('cl')
               i_cl = i
               print*,'cl',i_cl
               m_tr(i_cl) = 35.453
               type_tr(i_cl) = 1
            case('clo')
               i_clo = i
               print*,'clo',i_clo
               m_tr(i_clo) = 51.452
               type_tr(i_clo) = 1
            case('cl2')
               i_cl2 = i
               print*,'cl2',i_cl2
               m_tr(i_cl2) = 70.906
               type_tr(i_cl2) = 1
            case('hcl')
               i_hcl = i
               print*,'hcl',i_hcl
               m_tr(i_hcl) = 36.461
               type_tr(i_hcl) = 1
            case('hocl')
               i_hocl = i
               print*,'hocl',i_hocl
               m_tr(i_hocl) = 52.46
               type_tr(i_hocl) = 1
            case('clco')
               i_clco = i
               print*,'clco',i_clco
               m_tr(i_clco) = 63.463
               type_tr(i_clco) = 1
            case('clco3')
               i_clco3 = i
               print*,'clco3',i_clco3
               m_tr(i_clco3) = 95.462
               type_tr(i_clco3) = 1
            case('cocl2')
               i_cocl2 = i
               print*,'cocl2',i_cocl2
               m_tr(i_cocl2) = 98.916
               type_tr(i_cocl2) = 1
            case('s')
               i_s = i
               print*,'s',i_s
               m_tr(i_s) = 32.065
               type_tr(i_s) = 1
            case('so')
               i_so = i
               print*,'so',i_so
               m_tr(i_so) = 48.0644
               type_tr(i_so) = 1
            case('so2')
               i_so2 = i
               print*,'so2',i_so2
               m_tr(i_so2) = 64.064
               type_tr(i_so2) = 1
            case('so3')
               i_so3 = i
               print*,'so3',i_so3
               m_tr(i_so3) = 80.063
               type_tr(i_so3) = 1
            case('osso_cis')
               i_osso_cis = i
               print*,'osso_cis',i_osso_cis
               m_tr(i_osso_cis)= 96.1288
               type_tr(i_osso_cis) = 1
            case('osso_trans')
               i_osso_trans = i
               print*,'osso_trans',i_osso_trans
               m_tr(i_osso_trans)= 96.1288
               type_tr(i_osso_trans) = 1
            case('s2o2_cyc')
               i_s2o2_cyc = i
               print*,'s2o2_cyc',i_s2o2_cyc
               m_tr(i_s2o2_cyc)= 96.1288
               type_tr(i_s2o2_cyc) = 1
            case('ocs')
               i_ocs = i
               print*,'ocs',i_ocs
               m_tr(i_ocs) = 60.0751
               type_tr(i_ocs) = 1
            case('hso3')
               i_hso3 = i
               print*,'hso3',i_hso3
               m_tr(i_hso3) = 81.071
               type_tr(i_hso3) = 1
            case('h2so4')
               i_h2so4 = i
               print*,'h2so4',i_h2so4
               m_tr(i_h2so4) = 98.078
               type_tr(i_h2so4) = 1
            case('s2')
               i_s2 = i
               print*,'s2',i_s2
               m_tr(i_s2) = 64.13
               type_tr(i_s2) = 1
            case('clso2')
               i_clso2 = i
               print*,'clso2',i_clso2
               m_tr(i_clso2) = 99.517
               type_tr(i_clso2) = 1
            case('cl2so2')
               i_cl2so2 = i
               print*,'cl2so2',i_cl2so2
               m_tr(i_cl2so2) = 134.97
               type_tr(i_cl2so2) = 1
            case('oscl')
               i_oscl = i
               print*,'oscl',i_oscl
               m_tr(i_oscl) = 83.517
               type_tr(i_oscl) = 1
            case('n2')
               i_n2 = i
               print*,'n2',i_n2
               m_tr(i_n2) = 28.013
               type_tr(i_n2) = 1
            case('he')
               i_he = i
               print*,'he',i_he
               m_tr(i_he) = 4.0026
               type_tr(i_he) = 1
            case('n2d')
               i_n2d = i
               print*,'n2d',i_n2d
               m_tr(i_n2d) = 14.0067
               type_tr(i_n2d) = 1
            case('n')
               i_n = i
               print*,'n',i_n
               m_tr(i_n) = 14.0067
               type_tr(i_n) = 1
            case('no')
               i_no = i
               print*,'no',i_no
               m_tr(i_no) = 30.0061
               type_tr(i_no) = 1
            case('no2')
               i_no2 = i
               print*,'no2',i_no2
               m_tr(i_no2) = 46.0055
               type_tr(i_no2) = 1
            case('hd')
               i_hd = i
               print*,'hd',i_hd
               m_tr(i_hd) = 3.022042
               type_tr(i_hd) = 1
            case('d')
               i_d = i
               print*,'d',i_d
               m_tr(i_d) = 2.00855
               type_tr(i_d) = 1
            case('od')
               i_od = i
               print*,'od',i_od
               m_tr(i_od) = 18.0134
               type_tr(i_od) = 1
            case('do2')
               i_do2 = i
               print*,'do2',i_do2
               m_tr(i_do2) = 34.0128
               type_tr(i_do2) = 1
            case('hdo2')
               i_hdo2 = i
               print*,'hdo2',i_hdo2
               m_tr(i_hdo2) = 35.0218
               type_tr(i_hdo2) = 1
            case('hdo')
               i_hdo = i
               print*,'hdo',i_hdo
               m_tr(i_hdo) = 19.0214
               type_tr(i_hdo) = 1
            case('hdso4')
               i_hdso4 = i
               print*,'hdso4',i_hdso4
               m_tr(i_hdso4) = 99.084
               type_tr(i_hdso4) = 1
            case('dcl')
               i_dcl = i
               print*,'dcl',i_dcl
               m_tr(i_dcl) = 37.4671
               type_tr(i_dcl) = 1
            case('docl')
               i_docl = i
               print*,'docl',i_docl
               m_tr(i_docl) = 53.4665
               type_tr(i_docl) = 1

            ! ions

            case('co2plus')
               i_co2plus = i
               print*,'co2plus',i_co2plus
               m_tr(i_co2plus) = 44.0095
               type_tr(i_co2plus) = 2
            case('coplus')
               i_coplus = i
               print*,'coplus',i_coplus
               m_tr(i_coplus) = 28.0101
               type_tr(i_coplus) = 2
            case('oplus')
               i_oplus = i
               print*,'oplus',i_oplus
               m_tr(i_oplus) = 15.994
               type_tr(i_oplus) = 2
            case('o2plus')
               i_o2plus = i
               print*,'o2plus',i_o2plus
               m_tr(i_o2plus) = 31.9988
               type_tr(i_o2plus) = 2
            case('n2plus')
               i_n2plus = i
               print*,'n2plus',i_n2plus
               m_tr(i_n2plus) = 28.013
               type_tr(i_n2plus) = 2
            case('hplus')
               i_hplus = i
               print*,'hplus',i_hplus
               m_tr(i_hplus) = 1.00794
               type_tr(i_hplus) = 2
            case('h2oplus')
               i_h2oplus = i
               print*,'h2oplus',i_h2oplus
               m_tr(i_h2oplus) = 18.0153
               type_tr(i_h2oplus) = 2
            case('nplus')
               i_nplus = i
               print*,'nplus',i_nplus
               m_tr(i_nplus) = 14.0067
               type_tr(i_nplus) = 2
            case('ohplus')
               i_ohplus = i
               print*,'ohplus',i_ohplus
               m_tr(i_ohplus) = 17.0073
               type_tr(i_ohplus) = 2
            case('cplus')
               i_cplus = i
               print*,'cplus',i_cplus
               m_tr(i_cplus) = 12.011
               type_tr(i_cplus) = 2
            case('noplus')
               i_noplus = i
               print*,'noplus',i_noplus
               m_tr(i_noplus) = 30.0061
               type_tr(i_noplus) = 2
            case('h3oplus')
               i_h3oplus = i
               print*,'h3oplus',i_h3oplus
               m_tr(i_h3oplus) = 19.0232
               type_tr(i_h3oplus) = 2
            case('hcoplus')
               i_hcoplus = i
               print*,'hcoplus',i_hcoplus
               m_tr(i_hcoplus) = 29.0180
               type_tr(i_hcoplus) = 2
            case('hco2plus')
               i_hco2plus = i
               print*,'hco2plus',i_hco2plus
               m_tr(i_hco2plus) = 45.
               type_tr(i_hco2plus) = 2
            case('elec')
               i_elec = i
               print*,'elec',i_elec
               m_tr(i_elec) = 1./1822.89
               type_tr(i_elec) = 2
            case('dco2plus')
               i_dco2plus = i
               print*,'dco2plus',i_dco2plus
               m_tr(i_dco2plus) = 46.023
               type_tr(i_dco2plus) = 2
            case('dcoplus')
               i_dcoplus = i
               print*,'dcoplus',i_dcoplus
               m_tr(i_dcoplus) = 30.024
               type_tr(i_dcoplus) = 2
            case('hdoplus')
               i_hdoplus = i
               print*,'hdoplus',i_hdoplus
               m_tr(i_hdoplus) = 19.0214
               type_tr(i_hdoplus) = 2
            case('dplus')
               i_dplus = i
               print*,'dplus',i_dplus
               m_tr(i_dplus) = 2.00855
               type_tr(i_dplus) = 2
            case('odplus')
               i_odplus = i
               print*,'odplus',i_odplus
               m_tr(i_odplus) = 18.0134
               type_tr(i_odplus) = 2
            case('h2doplus')
               i_h2doplus = i
               print*,'h2doplus',i_h2doplus
               m_tr(i_h2doplus) = 20.0214
               type_tr(i_h2doplus) = 2

            ! liquid tracers (cl_scheme = 1)

            case('h2oliq')
               i_h2oliq = i
               print*,'h2oliq',i_h2oliq
               m_tr(i_h2oliq) = 18.0153
               type_tr(i_h2oliq) = 3
            case('h2so4liq')
               i_h2so4liq = i
               print*,'h2so4liq',i_h2so4liq
               m_tr(i_h2so4liq) = 98.078
               type_tr(i_h2so4liq) = 3
            case('hdoliq')
               i_hdoliq = i
               print*,'hdoliq',i_hdoliq
               m_tr(i_hdoliq) = 19.0214
               type_tr(i_hdoliq) = 3
            case('hdso4liq')
               i_hdso4liq = i
               print*,'hdso4liq',i_hdso4liq
               m_tr(i_hdso4liq) = 99.084
               type_tr(i_hdso4liq) = 3

            ! liquid tracers (cl_scheme = 2)

            case('m0_aer')
               i_m0_aer = i
               print*,'m0_aer',i_m0_aer
               type_tr(i_m0_aer) = 10
            case('m3_aer')
               i_m3_aer = i
               print*,'m3_aer',i_m3_aer
               type_tr(i_m3_aer) = 10
            case('m0_m1drop')
               i_m0_mode1drop = i
               print*,'m0_m1drop',i_m0_mode1drop
               type_tr(i_m0_mode1drop) = 10
            case('m0_m1ccn')
               i_m0_mode1ccn = i
               print*,'m0_m1ccn',i_m0_mode1ccn
               type_tr(i_m0_mode1ccn) = 10
            case('m3_m1sa')
               i_m3_mode1sa = i
               print*,'m3_m1sa',i_m3_mode1sa
               type_tr(i_m3_mode1sa) = 10
            case('m3_m1w')
               i_m3_mode1w = i
               print*,'m3_m1w',i_m3_mode1w
               type_tr(i_m3_mode1w) = 10
            case('m3_m1ccn')
               i_m3_mode1ccn = i
               print*,'m3_m1ccn',i_m3_mode1ccn
               type_tr(i_m3_mode1ccn) = 10
            case('m0_m2drop')
               i_m0_mode2drop = i
               print*,'m0_m2drop',i_m0_mode2drop
               type_tr(i_m0_mode2drop) = 10
            case('m0_m2ccn')
               i_m0_mode2ccn = i
               print*,'m0_m2ccn',i_m0_mode2ccn
               type_tr(i_m0_mode2ccn) = 10
            case('m3_m2sa')
               i_m3_mode2sa = i
               print*,'m3_m2sa',i_m3_mode2sa
               type_tr(i_m3_mode2sa) = 10
            case('m3_m2w')
               i_m3_mode2w = i
               print*,'m3_m2w',i_m3_mode2w
               type_tr(i_m3_mode2w) = 10
            case('m3_m2ccn')
               i_m3_mode2ccn = i
               print*,'m3_m2ccn',i_m3_mode2ccn
               type_tr(i_m3_mode2ccn) = 10
         end select
      end do

  end subroutine chemparam_ini

!============================================================================

  subroutine vapors4muphy_ini(nlon,nlev,trac)

!============================================================================

  use infotrac_phy, only: nqtot, tname

  integer :: nlon, nlev
  real    :: trac(nlon,nlev,nqtot) ! traceur ( en vmr)

!  integer :: i
!  real    :: trac1d(nlev,2) ! traceur lu ( en vmr)

! lecture d'un fichier texte contenant les profils de trac1d(:1) = h2o et trac1d(:,2) = h2so4
!  do i=1,nlon
!     trac(i,:,i_h2o) = trac1d(:,1)
!     trac(i,:,i_h2so4) = trac1d(:,2)
!  enddo

!  intitialisation profils altitude h2o et h2so4
!  profil h2o initial vap+liq == que vap
   trac(:,1:24,i_h2o) = 30.e-6 !
   trac(:,25:50,i_h2o) = 1.e-6 !
   trac(:,:,i_hdo) = trac(:,:,i_h2o)*0.032

   trac(:,:,i_h2so4) = 3.e-9 ! limite sup sandor 2012
   trac(:,23:50,i_h2so4) = 2.e-6 ! profil h2so4 initial => vap+liq
   trac(:,:,i_hdso4) = trac(:,:,i_h2so4)*0.032

  end subroutine vapors4muphy_ini

end module chemparam_mod
