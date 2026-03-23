module photolysis_mod

  implicit none

! photolysis

  integer, save :: nphot = 37

!$OMP THREADPRIVATE(nphot)

  integer, parameter :: nabs  = 26        ! number of absorbing gases

! spectral grid

  integer, parameter :: nw = 194          ! number of spectral intervals (low-res)
  integer, save :: mopt                   ! high-res/low-res switch

!$OMP THREADPRIVATE(mopt)

  real, dimension(nw), save :: wl, wc, wu ! lower, center, upper wavelength for each interval

!$OMP THREADPRIVATE(wl, wc, wu)

! solar flux

  real, dimension(nw), save :: f          !  solar flux (w.m-2.nm-1) at 1 au

! cross-sections and quantum yields

  real, dimension(nw), save :: xsco2_195, xsco2_295, xsco2_370        ! co2 absorption cross-section at 195-295-370 k (cm2)
  real, dimension(nw), save :: yieldco2                               ! co2 photodissociation yield
  real, dimension(nw), save :: xso2_150, xso2_200, xso2_250, xso2_300 ! o2 absorption cross-section at 150-200-250-300 k (cm2)
  real, dimension(nw), save :: yieldo2                                ! o2 photodissociation yield
  real, dimension(nw), save :: xso3_218, xso3_298                     ! o3 absorption cross-section at 218-298 k (cm2)
  real, dimension(nw), save :: xsh2o                                  ! h2o absorption cross-section (cm2)
  real, dimension(nw), save :: xsh2o2                                 ! h2o2 absorption cross-section (cm2)
  real, dimension(nw), save :: xsho2                                  ! ho2 absorption cross-section (cm2)
  real, dimension(nw), save :: xshcl                                  ! hcl absorption cross-section (cm2)
  real, dimension(nw), save :: xsdcl
  real, dimension(nw), save :: xscl2                                  ! cl2 absorption cross-section (cm2)
  real, dimension(nw), save :: xshocl                                 ! hocl absorption cross-section (cm2)
  real, dimension(nw), save :: xsclo                                  ! clo absorption cross-section (cm2)
  real, dimension(nw), save :: xsso2_200, xsso2_298, xsso2_360        ! so2 absorption cross-section (cm2)
  real, dimension(nw), save :: xsso_150, xsso_250                     ! so absorption cross-section (cm2)
  real, dimension(nw), save :: xsso3                                  ! so3 absorption cross-section (cm2)
  real, dimension(nw), save :: xss2                                   ! s2 absorption cross-section (cm2)
  real, dimension(nw), save :: xsosso_cis                             ! osso_cis absorption cross-section (cm2)
  real, dimension(nw), save :: xsosso_trans                           ! osso_trans absorption cross-section (cm2)
  real, dimension(nw), save :: xss2o2_cyc                             ! s2o2_cyc absorption cross-section (cm2)
  real, dimension(nw), save :: xsclso2                                ! clso2 absorption cross-section (cm2)
  real, dimension(nw), save :: xscl2so2                               ! cl2so2 absorption cross-section (cm2)
  real, dimension(nw), save :: xsocs                                  ! cos absorption cross-section (cm2)
  real, dimension(nw), save :: xscocl2                                ! cocl2 absorption cross-section (cm2)
  real, dimension(nw), save :: xsh2so4                                ! h2so4 absorption cross-section (cm2)
  real, dimension(nw), save :: xsh2                                   ! h2 absorption cross-section (cm2)
  real, dimension(nw), save :: yieldh2                                ! h2 photodissociation yield
  real, dimension(nw), save :: xsno2, xsno2_220, xsno2_294            ! no2 absorption cross-section at 220-294 k (cm2)
  real, dimension(nw), save :: yldno2_248, yldno2_298                 ! no2 quantum yield at 248-298 k
  real, dimension(nw), save :: xsno                                   ! no absorption cross-section (cm2)
  real, dimension(nw), save :: yieldno                                ! no photodissociation yield
  real, dimension(nw), save :: xsn2                                   ! n2 absorption cross-section (cm2)
  real, dimension(nw), save :: yieldn2                                ! n2 photodissociation yield
  real, dimension(nw), save :: xshdo                                  ! hdo absorption cross-section (cm2)
  real, dimension(nw), save :: albedo                                 ! surface albedo

!$OMP THREADPRIVATE(f, xsco2_195, xsco2_295, xsco2_370, yieldco2, xso2_150, xso2_200, xso2_250, xso2_300, yieldo2, xso3_218, xso3_298, xsclso2, xscl2so2, xsh2o, xsh2o2, xsho2, xsh2, yieldh2, xsno2, xsno2_220, xsno2_294, yldno2_248, yldno2_298)
!$OMP THREADPRIVATE(xsno, yieldno, xsn2, yieldn2, xshdo, albedo)

contains

subroutine init_photolysis

! initialise on-line photolysis

! mopt = 1 high-resolution
! mopt = 2 martian low-resolution (recommended for gcm use)
! mopt = 3 venusian low-resolution (recommended for gcm use)

  mopt = 3

! set wavelength grid

  call gridw(nw,wl,wc,wu,mopt)

! read and grid solar flux data
 
  call rdsolarflux(nw,wl,wc,f)

! read and grid o2 cross-sections

  call rdxso2(nw,wl,xso2_150,xso2_200,xso2_250,xso2_300,yieldo2)

! read and grid co2 cross-sections
 
  call rdxsco2(nw,wl,xsco2_195,xsco2_295,xsco2_370,yieldco2)

! read and grid o3 cross-sections

  call rdxso3(nw,wl,xso3_218,xso3_298)

! read and grid h2o cross-sections

  call rdxsh2o(nw,wl,xsh2o)

! read and grid h2o2 cross-sections

  call rdxsh2o2(nw,wl,xsh2o2)

! read and grid ho2 cross-sections

  call rdxsho2(nw,wl,xsho2)

! read and grid hcl cross-sections

  call rdxshcl(nw,wl,xshcl)
 
! read and grid dcl cross-sections
 
  call rdxsdcl(nw,wl,xsdcl)

! read and grid cl2 cross-sections

  call rdxscl2(nw,wl,xscl2)

! read and grid hocl cross-sections

  call rdxshocl(nw,wl,xshocl)

! read and grid clo cross-sections

  call rdxsclo(nw,wl,xsclo)

! read and grid so2 cross-sections

  call rdxsso2(nw,wl,xsso2_200,xsso2_298,xsso2_360)

! read and grid so cross-sections

  call rdxsso(nw,wl,xsso_150,xsso_250)

! read and grid so3 cross-sections

  call rdxsso3(nw,wl,xsso3)

! read and grid s2 cross-sections

  call rdxss2(nw,wl,xss2)

! read and grid osso cross-sections

  call rdxsosso(nw,wl,wc,xsosso_cis,xsosso_trans,xss2o2_cyc)

! read and grid clso2 cross-sections

  call rdxsclso2(nw,wl,wc,xsclso2)

! read and grid cl2so2 cross-sections

  call rdxscl2so2(nw,wl,wc,xscl2so2)

! read and grid ocs cross-sections

  call rdxsocs(nw,wl,xsocs)

! read and grid cocl2 cross-sections

  call rdxscocl2(nw,wl,xscocl2)

! read and grid h2so4 cross-sections

  call rdxsh2so4(nw,wl,xsh2so4)

! read and grid h2 cross-sections

  call rdxsh2(nw,wl,xsh2,yieldh2)

! read and grid no cross-sections

  call rdxsno(nw,wl,xsno,yieldno)

! read and grid no2 cross-sections

  call rdxsno2(nw,wl,xsno2,xsno2_220,xsno2_294,yldno2_248,yldno2_298)

! read and grid n2 cross-sections

  call rdxsn2(nw,wl,xsn2,yieldn2)

! read and grid hdo cross-sections

  call rdxshdo(nw,wl,xshdo)

! set surface albedo

  call setalb(nw,wl,albedo)

end subroutine init_photolysis

!==============================================================================

      subroutine gridw(nw,wl,wc,wu,mopt)

!     Create the wavelength grid for all interpolations and radiative transfer 
!     calculations.  Grid may be irregularly spaced.  Wavelengths are in nm.  
!     No gaps are allowed within the wavelength grid.                        

      implicit none

!     input
 
      integer :: nw      ! number of wavelength grid points
      integer :: mopt    ! high-res/low-res switch

!     output

      real, dimension(nw) :: wl, wc, wu   ! lower, center, upper wavelength for each interval

!     local

      real :: wincr    ! wavelength increment
      integer :: iw, kw

!     mopt = 1    high-resolution mode (3789 intervals)
!
!                   0-108 nm :  1.0  nm
!                 108-124 nm :  0.1  nm
!                 124-175 nm :  0.5  nm
!                 175-205 nm :  0.01 nm
!                 205-365 nm :  0.5  nm
!                 365-850 nm :  5.0  nm  
!
!     mopt = 2    low-resolution mode (162)
!
!                    0-60 nm :  6.0 nm
!                   60-80 nm :  2.0 nm
!                   80-85 nm :  5.0 nm
!                  85-117 nm :  2.0 nm
!                 117-120 nm :  5.0 nm
!                 120-123 nm :  0.2 nm
!                 123-163 nm :  5.0 nm
!                 163-175 nm :  2.0 nm
!                 175-205 nm :  0.5 nm
!                 205-245 nm :  5.0 nm
!                 245-415 nm : 10.0 nm
!                 415-815 nm : 50.0 nm
!
!     mopt = 3    venusian low-resolution mode (194)
!                 (205-245 nm -> 1nm of resolution)
!
!                    0-60 nm :  6.0 nm
!                   60-80 nm :  2.0 nm
!                   80-85 nm :  5.0 nm
!                  85-117 nm :  2.0 nm
!                 117-120 nm :  5.0 nm
!                 120-123 nm :  0.2 nm
!                 123-163 nm :  5.0 nm
!                 163-175 nm :  2.0 nm
!                 175-205 nm :  0.5 nm
!                 205-245 nm :  1.0 nm
!                 245-415 nm : 10.0 nm
!                 415-815 nm : 50.0 nm

      if (mopt == 1) then   ! high-res

! define wavelength intervals of width 1.0 nm from 0 to 108 nm:

      kw = 0
      wincr = 1.0
      do iw = 0, 107
        kw = kw + 1
        wl(kw) = real(iw)
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      end do

! define wavelength intervals of width 0.1 nm from 108 to 124 nm:

      wincr = 0.1
      do iw = 1080, 1239, 1
        kw = kw + 1
        wl(kw) = real(iw)/10.
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      end do

! define wavelength intervals of width 0.5 nm from 124 to 175 nm:

      wincr = 0.5
      do iw = 1240, 1745, 5
        kw = kw + 1
        wl(kw) = real(iw)/10.
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      end do

! define wavelength intervals of width 0.01 nm from 175 to 205 nm:

      wincr = 0.01
      do iw = 17500, 20499, 1
        kw = kw + 1
        wl(kw) = real(iw)/100.
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      end do

! define wavelength intervals of width 0.5 nm from 205 to 365 nm:

      wincr = 0.5
      do iw = 2050, 3645, 5
        kw = kw + 1
        wl(kw) = real(iw)/10.
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      end do

! define wavelength intervals of width 5.0 nm from 365 to 855 nm:

      wincr = 5.0
      do iw = 365, 850, 5
        kw = kw + 1
        wl(kw) = real(iw)
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      end do
      wl(kw+1) = wu(kw)

!============================================================

      else if (mopt == 2) then   ! low-res

! define wavelength intervals of width 6.0 nm from 0 to 60 nm:

      kw = 0
      wincr = 6.0
      DO iw = 0, 54, 6
        kw = kw + 1
        wl(kw) = real(iw)
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      END DO

! define wavelength intervals of width 2.0 nm from 60 to 80 nm:

      wincr = 2.0
      DO iw = 60, 78, 2
        kw = kw + 1
        wl(kw) = real(iw)
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      END DO

! define wavelength intervals of width 5.0 nm from 80 to 85 nm:

      wincr = 5.0
      DO iw = 80, 80
        kw = kw + 1
        wl(kw) = real(iw)
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      END DO

! define wavelength intervals of width 2.0 nm from 85 to 117 nm:

      wincr = 2.0
      DO iw = 85, 115, 2
        kw = kw + 1
        wl(kw) = real(iw)
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      END DO

! define wavelength intervals of width 3.0 nm from 117 to 120 nm:

      wincr = 3.0
      DO iw = 117, 117
        kw = kw + 1
        wl(kw) = real(iw)
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      END DO

! define wavelength intervals of width 0.2 nm from 120 to 123 nm:

      wincr = 0.2
      DO iw = 1200, 1228, 2
        kw = kw + 1
        wl(kw) = real(iw)/10.
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      ENDDO

! define wavelength intervals of width 5.0 nm from 123 to 163 nm:

      wincr = 5.0
      DO iw = 123, 158, 5
        kw = kw + 1
        wl(kw) = real(iw)
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      ENDDO

! define wavelength intervals of width 2.0 nm from 163 to 175 nm:

      wincr = 2.0
      DO iw = 163, 173, 2
        kw = kw + 1
        wl(kw) = real(iw)
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      ENDDO

! define wavelength intervals of width 0.5 nm from 175 to 205 nm:

      wincr = 0.5
      DO iw = 1750, 2045, 5
        kw = kw + 1
        wl(kw) = real(iw)/10.
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      ENDDO

! define wavelength intervals of width 5.0 nm from 205 to 245 nm:

      wincr = 5.
      DO iw = 205, 240, 5
        kw = kw + 1
        wl(kw) = real(iw)
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      ENDDO

! define wavelength intervals of width 10.0 nm from 245 to 415 nm:

      wincr = 10.0
      DO iw = 245, 405, 10
        kw = kw + 1
        wl(kw) = real(iw)
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      ENDDO

! define wavelength intervals of width 50.0 nm from 415 to 815 nm:

      wincr = 50.0
      DO iw = 415, 815, 50
        kw = kw + 1
        wl(kw) = real(iw)
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      ENDDO

      wl(kw+1) = wu(kw)

!============================================================

      else if (mopt == 3) then   ! low-res

! define wavelength intervals of width 6.0 nm from 0 to 60 nm:

      kw = 0
      wincr = 6.0
      DO iw = 0, 54, 6
        kw = kw + 1
        wl(kw) = real(iw)
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      END DO

! define wavelength intervals of width 2.0 nm from 60 to 80 nm:

      wincr = 2.0
      DO iw = 60, 78, 2
        kw = kw + 1
        wl(kw) = real(iw)
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      END DO

! define wavelength intervals of width 5.0 nm from 80 to 85 nm:

      wincr = 5.0
      DO iw = 80, 80
        kw = kw + 1
        wl(kw) = real(iw)
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      END DO

! define wavelength intervals of width 2.0 nm from 85 to 117 nm:

      wincr = 2.0
      DO iw = 85, 115, 2
        kw = kw + 1
        wl(kw) = real(iw)
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      END DO

! define wavelength intervals of width 3.0 nm from 117 to 120 nm:

      wincr = 3.0
      DO iw = 117, 117
        kw = kw + 1
        wl(kw) = real(iw)
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      END DO

! define wavelength intervals of width 0.2 nm from 120 to 123 nm:

      wincr = 0.2
      DO iw = 1200, 1228, 2
        kw = kw + 1
        wl(kw) = real(iw)/10.
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      ENDDO

! define wavelength intervals of width 5.0 nm from 123 to 163 nm:

      wincr = 5.0
      DO iw = 123, 158, 5
        kw = kw + 1
        wl(kw) = real(iw)
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      ENDDO

! define wavelength intervals of width 2.0 nm from 163 to 175 nm:

      wincr = 2.0
      DO iw = 163, 173, 2
        kw = kw + 1
        wl(kw) = real(iw)
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      ENDDO

! define wavelength intervals of width 0.5 nm from 175 to 205 nm:

      wincr = 0.5
      DO iw = 1750, 2045, 5
        kw = kw + 1
        wl(kw) = real(iw)/10.
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      ENDDO

! define wavelength intervals of width 5.0 nm from 205 to 245 nm:

      wincr = 1
      DO iw = 205, 244, 1
        kw = kw + 1
        wl(kw) = real(iw)
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      ENDDO

! define wavelength intervals of width 10.0 nm from 245 to 415 nm:

      wincr = 10.0
      DO iw = 245, 405, 10
        kw = kw + 1
        wl(kw) = real(iw)
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      ENDDO

! define wavelength intervals of width 50.0 nm from 415 to 815 nm:

      wincr = 50.0
      DO iw = 415, 815, 50
        kw = kw + 1
        wl(kw) = real(iw)
        wu(kw) = wl(kw) + wincr
        wc(kw) = (wl(kw) + wu(kw))/2.
      ENDDO

      wl(kw+1) = wu(kw)

      print*, 'number of spectral intervals : ', kw+1

      endif
      
      end subroutine gridw

!==============================================================================
 
      subroutine rdsolarflux(nw,wl,wc,f)

!     Read and re-grid solar flux data.            

      USE mod_phys_lmdz_para, ONLY: is_master
      USE mod_phys_lmdz_transfert_para, ONLY: bcast
      USE clesphys_mod

      implicit none

!#include "clesphys.h" ! fixed_euv_value

!     input
 
      integer :: nw      ! number of wavelength grid points
      real, dimension(nw) :: wl, wc   ! lower and central wavelength for each interval

!     output

      real, dimension(nw) :: f  ! solar flux (w.m-2.nm-1)

!     local 

      integer, parameter :: kdata = 20000    ! max dimension of input solar flux
      integer :: msun   ! choice of solar flux
      integer :: iw, nhead, ihead, n, i, ierr, kin1, kin2

      real, parameter :: deltax = 1.e-4
      !  atlas1_thuillier_tuv.txt
      real, dimension(kdata) :: x1, y1      ! input solar flux - HIGH SOLAR ACTIVITY
      real :: E107_max = 196
      !  atlas3_thuillier_tuv.txt
      real, dimension(kdata) :: x2, y2      ! input solar flux - "LOW" SOLAR ACTIVITY
      real :: E107_min = 97
      ! be careful, x2 need to be equal to x1
      real, dimension(kdata) :: x0, y0      ! input solar flux - interpolated solar activity
      real, dimension(nw)    :: yg0         ! gridded solar flux
      real            :: factor_interp
            
      character(len=100) :: fil1, fil2

      kin1 = 10    ! input logical unit
      kin2 = 11    ! input logical unit
      
! select desired extra-terrestrial solar irradiance, using msun:

! 18 = atlas1_thuillier_tuv.txt  0-900 nm  March 29, 1992
!      Article: F10.7 =  192 s.f.u | <F10.7>81d =  171 s.f.u | <Rz> = 121 (sunsport number)
!      Thuillier et al., Adv. Space. Res., 34, 256-261, 2004 
!      Model SOLAR2000 version 2015/12: E10.7 = 195.8 s.f.u | <E10.7>81d = 196.9 s.f.u
!      Choix de la limite haute: E10.7 = 196 s.f.u

! 20 = atlas3_thuillier_tuv.txt  0-900 nm  November 11, 1994
!      Article: F10.7 = 77.5 s.f.u | <F10.7>81d = 83.5 s.f.u | <Rz> = 20 (sunsport number)   
!      Thuillier et al., Adv. Space. Res., 34, 256-261, 2004 
!      Model SOLAR2000 version 2015/12: E10.7 =  96.9 s.f.u | <E10.7>81d = 100.0 s.f.u
!      Choix de la limite basse: E10.7 =  97 s.f.u

      if (fixed_euv_value .ge. 196) THEN
        msun = 18
        print*, 'Atlas1 spectrum Thuiller chosen'
      else
        msun = 19
        print*, 'interpolated Spectrum case'
      end if

      IF (is_master) THEN
      
        fil1 = 'solar_fluxes/atlas1_thuillier_tuv.txt'
        print*, 'solar flux : ', fil1
        
        open(kin1, file=fil1, status='old', iostat=ierr)
        
        if (ierr /= 0) THEN 
          write(*,*)'cant find solar flux : ', fil1
          write(*,*)'It should be in : INPUT/solar_fluxes'
          write(*,*)'1) You can change this directory address in '
          write(*,*)'   callphys.def with datadir=/path/to/dir'
          write(*,*)'2) If necessary, /solar fluxes (and other datafiles)'
          write(*,*)'   can be obtained online on:'
          write(*,*)' http://www.lmd.jussieu.fr/~lmdz/planets/mars/datadir'
          stop
        end if
          
        nhead = 9
        n = 19193
        DO ihead = 1, nhead
          READ(kin1,*)
        ENDDO
        DO i = 1, n
          READ(kin1,*) x1(i), y1(i)
          y1(i) = y1(i)*1.e-3    ! mw -> w
        ENDDO
        CLOSE (kin1)
          
        fil2 = 'solar_fluxes/atlas3_thuillier_tuv.txt'
        print*, 'solar flux : ', fil2
        
        open(kin2, file=fil2, status='old', iostat=ierr)
        
        if (ierr /= 0) THEN 
          write(*,*)'cant find solar flux : ', fil2
          write(*,*)'It should be in : INPUT/solar_fluxes'
          write(*,*)'1) You can change this directory address in '
          write(*,*)'   callphys.def with datadir=/path/to/dir'
        write(*,*)'2) If necessary, /solar fluxes (and other datafiles)'
          write(*,*)'   can be obtained online on:'
        write(*,*)' http://www.lmd.jussieu.fr/~lmdz/planets/mars/datadir'
          stop
        end if

        nhead = 9
        n = 19193
        DO ihead = 1, nhead
          READ(kin2,*)
        ENDDO
        DO i = 1, n
          READ(kin2,*) x2(i), y2(i)
          y2(i) = y2(i)*1.e-3    ! mw -> w
        ENDDO
        CLOSE (kin2)

        IF (msun .eq. 18) THEN
          DO i = 1, n
            x0(i) = x1(i)
            y0(i) = y1(i)
          ENDDO
        ELSE
          ! INTERPOLATION BETWEEN E107_min to E107_max and extrapolation below E107_min
          factor_interp=(fixed_euv_value-E107_min)/(E107_max-E107_min)
          DO i = 1, n
            x0(i) = x1(i)
            y0(i) = y2(i) + (y1(i) - y2(i))* factor_interp
          ENDDO
        ENDIF ! msun .ne. 18
        
        CALL addpnt(x0,y0,kdata,n,x0(1)*(1.-deltax),0.)
        CALL addpnt(x0,y0,kdata,n,          0.,0.)
        CALL addpnt(x0,y0,kdata,n,x0(n)*(1.+deltax),0.)
        CALL addpnt(x0,y0,kdata,n,      1.e+38,0.)
        CALL inter2(nw,wl,yg0,n,x0,y0,ierr)

        IF (ierr .NE. 0) THEN
          WRITE(*,*) ierr, fil1, fil2
          STOP
        ENDIF

!     convert to photon.s-1.nm-1.cm-2
!     5.039e11 = 1.e-4*1e-9/(hc = 6.62e-34*2.998e8)

        DO iw = 1, nw-1
          f(iw) = yg0(iw)*wc(iw)*5.039e11
        ENDDO

      endif !is_master

      call bcast(f)

      end subroutine rdsolarflux

!==============================================================================

      subroutine addpnt ( x, y, ld, n, xnew, ynew )

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Add a point <xnew,ynew> to a set of data pairs <x,y>.  x must be in      =*
!=  ascending order                                                          =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  X    - REAL vector of length LD, x-coordinates                       (IO)=*
!=  Y    - REAL vector of length LD, y-values                            (IO)=*
!=  LD   - INTEGER, dimension of X, Y exactly as declared in the calling  (I)=*
!=         program                                                           =*
!=  N    - INTEGER, number of elements in X, Y.  On entry, it must be:   (IO)=*
!=         N < LD.  On exit, N is incremented by 1.                          =*
!=  XNEW - REAL, x-coordinate at which point is to be added               (I)=*
!=  YNEW - REAL, y-value of point to be added                             (I)=*
!-----------------------------------------------------------------------------*

      IMPLICIT NONE

! calling parameters

      INTEGER ld, n
      REAL x(ld), y(ld)
      REAL xnew, ynew
      INTEGER ierr

! local variables

      INTEGER insert
      INTEGER i

!-----------------------------------------------------------------------

! initialize error flag

      ierr = 0

! check n<ld to make sure x will hold another point

      IF (n .GE. ld) THEN
         WRITE(0,*) '>>> ERROR (ADDPNT) <<<  Cannot expand array '
         WRITE(0,*) '                        All elements used.'
         STOP
      ENDIF

      insert = 1
      i = 2

! check, whether x is already sorted.
! also, use this loop to find the point at which xnew needs to be inserted
! into vector x, if x is sorted.

 10   CONTINUE
      IF (i .LT. n) THEN
        IF (x(i) .LT. x(i-1)) THEN
           print*, x(i-1), x(i)
           WRITE(0,*) '>>> ERROR (ADDPNT) <<<  x-data must be in ascending order!'
           STOP
        ELSE
           IF (xnew .GT. x(i)) insert = i + 1
        ENDIF
        i = i+1
        GOTO 10
      ENDIF

! if <xnew,ynew> needs to be appended at the end, just do so,
! otherwise, insert <xnew,ynew> at position INSERT

      IF ( xnew .GT. x(n) ) THEN
 
         x(n+1) = xnew
         y(n+1) = ynew
  
      ELSE

! shift all existing points one index up

         DO i = n, insert, -1
           x(i+1) = x(i)
           y(i+1) = y(i)
         ENDDO

! insert new point

         x(insert) = xnew
         y(insert) = ynew
  
      ENDIF

! increase total number of elements in x, y

      n = n+1

      end subroutine addpnt

!==============================================================================

      subroutine inter2(ng,xg,yg,n,x,y,ierr)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Map input data given on single, discrete points onto a set of target     =*
!=  bins.                                                                    =*
!=  The original input data are given on single, discrete points of an       =*
!=  arbitrary grid and are being linearly interpolated onto a specified set  =*
!=  of target bins.  In general, this is the case for most of the weighting  =*
!=  functions (action spectra, molecular cross section, and quantum yield    =*
!=  data), which have to be matched onto the specified wavelength intervals. =*
!=  The average value in each target bin is found by averaging the trapezoi- =*
!=  dal area underneath the input data curve (constructed by linearly connec-=*
!=  ting the discrete input values).                                         =*
!=  Some caution should be used near the endpoints of the grids.  If the     =*
!=  input data set does not span the range of the target grid, an error      =*
!=  message is printed and the execution is stopped, as extrapolation of the =*
!=  data is not permitted.                                                   =*
!=  If the input data does not encompass the target grid, use ADDPNT to      =*
!=  expand the input array.                                                  =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NG  - INTEGER, number of bins + 1 in the target grid                  (I)=*
!=  XG  - REAL, target grid (e.g., wavelength grid);  bin i is defined    (I)=*
!=        as [XG(i),XG(i+1)] (i = 1..NG-1)                                   =*
!=  YG  - REAL, y-data re-gridded onto XG, YG(i) specifies the value for  (O)=*
!=        bin i (i = 1..NG-1)                                                =*
!=  N   - INTEGER, number of points in input grid                         (I)=*
!=  X   - REAL, grid on which input data are defined                      (I)=*
!=  Y   - REAL, input y-data                                              (I)=*
!-----------------------------------------------------------------------------*

      IMPLICIT NONE

! input:
      INTEGER ng, n
      REAL x(n), y(n), xg(ng)

! output:
      REAL yg(ng)

! local:
      REAL area, xgl, xgu
      REAL darea, slope
      REAL a1, a2, b1, b2
      INTEGER ngintv
      INTEGER i, k, jstart
      INTEGER ierr
!_______________________________________________________________________

      ierr = 0

!  test for correct ordering of data, by increasing value of x

      DO 10, i = 2, n
         IF (x(i) .LE. x(i-1)) THEN
            ierr = 1
            WRITE(*,*)'data not sorted'
            WRITE(*,*) x(i), x(i-1)
            RETURN
         ENDIF
   10 CONTINUE     

      DO i = 2, ng
        IF (xg(i) .LE. xg(i-1)) THEN
           ierr = 2
          WRITE(0,*) '>>> ERROR (inter2) <<<  xg-grid not sorted!'
          RETURN
        ENDIF
      ENDDO

! check for xg-values outside the x-range

      IF ( (x(1) .GT. xg(1)) .OR. (x(n) .LT. xg(ng)) ) THEN
          WRITE(0,*) '>>> ERROR (inter2) <<<  Data do not span grid.  '
          WRITE(0,*) '                        Use ADDPNT to expand data and re-run.'
          STOP
      ENDIF

!  find the integral of each grid interval and use this to 
!  calculate the average y value for the interval      
!  xgl and xgu are the lower and upper limits of the grid interval

      jstart = 1
      ngintv = ng - 1
      DO 50, i = 1,ngintv

! initialize:

            area = 0.0
            xgl = xg(i)
            xgu = xg(i+1)

!  discard data before the first grid interval and after the 
!  last grid interval
!  for internal grid intervals, start calculating area by interpolating
!  between the last point which lies in the previous interval and the
!  first point inside the current interval

            k = jstart
            IF (k .LE. n-1) THEN

!  if both points are before the first grid, go to the next point
   30         CONTINUE
                IF (x(k+1) .LE. xgl) THEN
                   jstart = k - 1
                   k = k+1
                   IF (k .LE. n-1) GO TO 30
                ENDIF


!  if the last point is beyond the end of the grid, complete and go to the next
!  grid
   40         CONTINUE
                 IF ((k .LE. n-1) .AND. (x(k) .LT. xgu)) THEN          

                    jstart = k-1

! compute x-coordinates of increment

                    a1 = MAX(x(k),xgl)
                    a2 = MIN(x(k+1),xgu)

! if points coincide, contribution is zero

                    IF (x(k+1).EQ.x(k)) THEN
                       darea = 0.e0
                    ELSE
                       slope = (y(k+1) - y(k))/(x(k+1) - x(k))
                       b1 = y(k) + slope*(a1 - x(k))
                       b2 = y(k) + slope*(a2 - x(k))
                       darea = (a2 - a1)*(b2 + b1)/2.
                    ENDIF

!  find the area under the trapezoid from a1 to a2

                    area = area + darea

! go to next point
              
                    k = k+1
                    GO TO 40

                ENDIF
            ENDIF

!  calculate the average y after summing the areas in the interval

            yg(i) = area/(xgu - xgl)

   50 CONTINUE

      end subroutine inter2

!==============================================================================

      subroutine rdxsco2(nw,wl,xsco2_195,xsco2_295,xsco2_370,yieldco2)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read and grid CO2 absorption cross-sections and photodissociation yield   =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  XSCO2  - REAL, molecular absoprtion cross section (cm^2) of CO2 at    (O)=*
!=           each specified wavelength                                       =*
!-----------------------------------------------------------------------------*

      USE mod_phys_lmdz_para, ONLY: is_master
      USE mod_phys_lmdz_transfert_para, ONLY: bcast

      implicit none

!     input

      integer :: nw               ! number of wavelength grid points
      real, dimension(nw) :: wl   ! lower and central wavelength for each interval

!     output

      real, dimension(nw) :: xsco2_195, xsco2_295, xsco2_370 ! co2 cross-sections (cm2)
      real, dimension(nw) :: yieldco2                        ! co2 photodissociation yield

!     local

      integer, parameter :: kdata = 42000
      real, parameter :: deltax = 1.e-4
      real, dimension(kdata) :: x1, y1, y2, y3, xion, ion
      real, dimension(nw) :: yg
      real :: xl, xu
      integer :: ierr, i, l, n, n1, n2, n3, n4
      CHARACTER*100 fil

      integer :: kin, kout ! input/ouput logical units

      kin  = 10
      kout = 30

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     CO2 absorption cross-sections
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     195K: huestis and berkowitz (2010) + Starck et al. (2006) 
!           + Yoshino et al. (1996) + Parkinson et al. (2003) + extrapolation
!
!     295K: huestis and berkowitz (2010) + Starck et al. (2006) 
!           + Yoshino et al. (1996) + Parkinson et al. (2003) + extrapolation
!
!     370K: huestis and berkowitz (2010) + Starck et al. (2006) 
!           + Lewis and Carver (1983) + extrapolation
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      n1 = 40769
      n2 = 41586
      n3 = 10110

!     195K:

      fil = 'cross_sections/co2_euv_uv_2018_195k.txt'
      print*, 'section efficace CO2 195K: ', fil

      if(is_master) then

      OPEN(UNIT=kin,FILE=fil,STATUS='old')
      DO i = 1,11 
         read(kin,*)
      END DO

      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
      END DO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,          0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,      1.e+38,0.)
      CALL inter2(nw,wl,yg,n1,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      ENDIF
      
      DO l = 1, nw-1
         xsco2_195(l) = yg(l)
      END DO

      endif !is_master

      call bcast(xsco2_195)

!     295K:

      fil = 'cross_sections/co2_euv_uv_2018_295k.txt'
      print*, 'section efficace CO2 295K: ', fil

      if(is_master) then

      OPEN(UNIT=kin,FILE=fil,STATUS='old')
      DO i = 1,11
         read(kin,*)
      END DO

      DO i = 1, n2
         READ(kin,*) x1(i), y1(i)
      END DO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n2,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n2,          0.,0.)
      CALL addpnt(x1,y1,kdata,n2,x1(n2)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n2,      1.e+38,0.)
      CALL inter2(nw,wl,yg,n2,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      ENDIF

      DO l = 1, nw-1
         xsco2_295(l) = yg(l)
      END DO

      endif !is_master

      call bcast(xsco2_295)

!     370K:

      fil = 'cross_sections/co2_euv_uv_2018_370k.txt'
      print*, 'section efficace CO2 370K: ', fil

      if(is_master) then

      OPEN(UNIT=kin,FILE=fil,STATUS='old')
      DO i = 1,11
         read(kin,*)
      END DO

      DO i = 1, n3
         READ(kin,*) x1(i), y1(i)
      END DO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n3,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n3,          0.,0.)
      CALL addpnt(x1,y1,kdata,n3,x1(n3)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n3,      1.e+38,0.)
      CALL inter2(nw,wl,yg,n3,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      ENDIF

      DO l = 1, nw-1
         xsco2_370(l) = yg(l)
      END DO

      endif !is_master

      call bcast(xsco2_370)

!     photodissociation yield:

      fil = 'cross_sections/efdis_co2-o2_schunkandnagy2000.txt'
      print*, 'photodissociation yield CO2: ', fil

      if(is_master) then

      OPEN(UNIT=kin,FILE=fil,STATUS='old')

      do i = 1,3
         read(kin,*)
      end do

      n4 = 17
      do i = 1, n4
         read(kin,*) xl, xu, ion(i)
         xion(i) = (xl + xu)/2.
         ion(i) = max(ion(i), 0.)
      end do
      close(kin)

      CALL addpnt(xion,ion,kdata,n4,xion(1)*(1.-deltax),0.)
      CALL addpnt(xion,ion,kdata,n4,          0.,0.)
      CALL addpnt(xion,ion,kdata,n4,xion(n4)*(1.+deltax),1.)
      CALL addpnt(xion,ion,kdata,n4,      1.e+38,1.)
      CALL inter2(nw,wl,yieldco2,n4,xion,ion,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      ENDIF

      endif !is_master

      call bcast(yieldco2)

!     DO l = 1, nw-1
!        write(kout,*) wl(l), xsco2_195(l),
!    $                        xsco2_295(l),
!    $                        xsco2_370(l),
!    $                        yieldco2(l)
!     END DO

      end subroutine rdxsco2

!==============================================================================
 
      subroutine rdxso2(nw,wl,xso2_150,xso2_200,xso2_250,xso2_300,yieldo2)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read and grid O2 cross-sections and photodissociation yield              =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW      - INTEGER, number of specified intervals + 1 in working       (I)=*
!=            wavelength grid                                                =*
!=  WL      - REAL, vector of lower limits of wavelength intervals in     (I)=*
!=            working wavelength grid                                        =*
!=  XSO2    - REAL, molecular absorption cross section                       =*
!-----------------------------------------------------------------------------*

      USE mod_phys_lmdz_para, ONLY: is_master
      USE mod_phys_lmdz_transfert_para, ONLY: bcast

      implicit none

!     input

      integer :: nw               ! number of wavelength grid points
      real, dimension(nw) :: wl   ! lower and central wavelength for each interval

!     output

      real, dimension(nw) :: xso2_150, xso2_200, xso2_250, xso2_300 ! o2 cross-sections (cm2)
      real, dimension(nw) :: yieldo2                                ! o2 photodissociation yield

!     local

      integer, parameter :: kdata = 18000
      real, parameter :: deltax = 1.e-4
      real, dimension(kdata) :: x1, y1, x2, y2, x3, y3, x4, y4
      real, dimension(kdata) :: xion, ion
      real    :: factor, xl, xu, dummy
      integer :: i, ierr, n, n1, n2, n3, n4, nhead
      integer :: kin, kout ! input/output logical units
      character*100 fil

      kin  = 10
      kout = 30

!     read o2 cross section data

      nhead = 22
      n     = 17434 

      fil = 'cross_sections/o2_composite_2018_150K.txt'
      print*, 'section efficace O2 150K: ', fil

      if(is_master) then

      open(kin, file=fil, status='old', iostat=ierr)

      if (ierr /= 0) THEN 
         write(*,*)'cant find O2 cross-sections : ', fil
         write(*,*)'It should be in :INPUT/cross_sections'
         write(*,*)'1) You can have to setup the link to the dir '
         write(*,*)'2) If necessary, /cross_sections (and other datafiles)'
         write(*,*)'   can be obtained online on:'
         write(*,*)' http://www.lmd.jussieu.fr/~lmdz/planets/mars/datadir'
         stop
      end if

      DO i = 1,nhead 
         read(kin,*)
      END DO

      n1 = n
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
      END DO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,               0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)
      CALL inter2(nw,wl,xso2_150,n1,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      ENDIF

      endif !is_master

      call bcast(xso2_150)

      fil = 'cross_sections/o2_composite_2018_200K.txt'
      print*, 'section efficace O2 200K: ', fil

      if(is_master) then

      OPEN(UNIT=kin,FILE=fil,STATUS='old')

      DO i = 1,nhead 
         read(kin,*)
      END DO

      n2 = n
      DO i = 1, n2
         READ(kin,*) x2(i), y2(i)
      END DO
      CLOSE (kin)

      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,               0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,           1.e+38,0.)
      CALL inter2(nw,wl,xso2_200,n2,x2,y2,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      ENDIF

      endif !is_master

      call bcast(xso2_200)

      fil = 'cross_sections/o2_composite_2018_250K.txt'
      print*, 'section efficace O2 250K: ', fil

      if(is_master) then

      OPEN(UNIT=kin,FILE=fil,STATUS='old')

      DO i = 1,nhead 
         read(kin,*)
      END DO

      n3 = n
      DO i = 1, n3
         READ(kin,*) x3(i), y3(i)
      END DO
      CLOSE (kin)

      CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,               0.,0.)
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,           1.e+38,0.)
      CALL inter2(nw,wl,xso2_250,n3,x3,y3,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      ENDIF

      endif !is_master

      call bcast(xso2_250)

      fil = 'cross_sections/o2_composite_2018_300K.txt'
      print*, 'section efficace O2 300K: ', fil

      if(is_master) then

      OPEN(UNIT=kin,FILE=fil,STATUS='old')

      DO i = 1,nhead 
         read(kin,*)
      END DO

      n4 = n
      DO i = 1, n4
         READ(kin,*) x4(i), y4(i)
      END DO
      CLOSE (kin)

      CALL addpnt(x4,y4,kdata,n4,x4(1)*(1.-deltax),0.)
      CALL addpnt(x4,y4,kdata,n4,               0.,0.)
      CALL addpnt(x4,y4,kdata,n4,x4(n4)*(1.+deltax),0.)
      CALL addpnt(x4,y4,kdata,n4,           1.e+38,0.)
      CALL inter2(nw,wl,xso2_300,n4,x4,y4,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      ENDIF

      endif !is_master

      call bcast(xso2_300)

!     photodissociation yield

      fil = 'cross_sections/efdis_co2-o2_schunkandnagy2000.txt'

      if(is_master) then

      OPEN(UNIT=kin,FILE=fil,STATUS='old')

      do i = 1,11
         read(kin,*)
      end do

      n = 9
      DO i = 1, n
         READ(kin,*) xl, xu, dummy, ion(i)
         xion(i) = (xl + xu)/2.
         ion(i) = max(ion(i), 0.)
      END DO
      CLOSE (kin)

      CALL addpnt(xion,ion,kdata,n,xion(1)*(1.-deltax),0.)
      CALL addpnt(xion,ion,kdata,n,          0.,0.)
      CALL addpnt(xion,ion,kdata,n,xion(n)*(1.+deltax),1.)
      CALL addpnt(xion,ion,kdata,n,      1.e+38,1.)
      CALL inter2(nw,wl,yieldo2,n,xion,ion,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      ENDIF

      endif !is_master

      call bcast(yieldo2)

      end subroutine rdxso2

!==============================================================================
 
      subroutine rdxso3(nw,wl,xso3_218,xso3_298)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read ozone molecular absorption cross section.  Re-grid data to match    =*
!=  specified wavelength working grid.                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  XSO3_218 REAL, molecular absoprtion cross section (cm^2) of O3 at     (O)=*
!=           each specified wavelength (JPL 2006)  218 K                     =*
!=  XSO3_298 REAL, molecular absoprtion cross section (cm^2) of O3 at     (O)=*
!=           each specified wavelength (JPL 2006)  298 K                     =*
!-----------------------------------------------------------------------------*

      USE mod_phys_lmdz_para, ONLY: is_master
      USE mod_phys_lmdz_transfert_para, ONLY: bcast

      implicit none

!     input

      integer :: nw               ! number of wavelength grid points
      real, dimension(nw) :: wl   ! lower and central wavelength for each interval

!     output

      real, dimension(nw) :: xso3_218, xso3_298 ! o3 cross-sections (cm2)

!     local

      integer, parameter :: kdata = 200
      real, parameter :: deltax = 1.e-4
      real, dimension(kdata) :: x1, x2, y1, y2
      real, dimension(nw) :: yg
      real :: a1, a2

      integer :: i, ierr, iw, n, n1, n2
      integer :: kin, kout ! input/output logical units

      character*100 fil

      kin  = 10

!     JPL 2006 218 K

      fil = 'cross_sections/o3_cross-sections_jpl_2006_218K.txt' 
      print*, 'section efficace O3 218K: ', fil

      if(is_master) then

      OPEN(UNIT=kin,FILE=fil,STATUS='old')
      n1 = 167
      DO i = 1, n1
         READ(kin,*) a1, a2, y1(i)
         x1(i) = (a1+a2)/2.
      END DO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,          0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,      1.e+38,0.)
      CALL inter2(nw,wl,yg,n1,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      ENDIF

      DO iw = 1, nw-1
         xso3_218(iw) = yg(iw)
      END DO

      endif !is_master

      call bcast(xso3_218)

!     JPL 2006 298 K

      fil = 'cross_sections/o3_cross-sections_jpl_2006_298K.txt' 
      print*, 'section efficace O3 298K: ', fil

      if(is_master) then

      OPEN(UNIT=kin,FILE=fil,STATUS='old')
      n2 = 167
      DO i = 1, n2
         READ(kin,*) a1, a2, y2(i)
         x2(i) = (a1+a2)/2.
      END DO
      CLOSE (kin)

      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,          0.,0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,      1.e+38,0.)
      CALL inter2(nw,wl,yg,n2,x2,y2,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      ENDIF

      DO iw = 1, nw-1
         xso3_298(iw) = yg(iw)
      END DO

      endif !is_master

      call bcast(xso3_298)

      end subroutine rdxso3

!==============================================================================

      subroutine rdxss2(nw, wl, yg)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read S2 molecular absorption cross section.  Re-grid data to match      =*
!=  specified wavelength working grid.                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  YG     - REAL, molecular absoprtion cross section (cm^2) of H2O at    (O)=*
!=           each specified wavelength                                       =*
!-----------------------------------------------------------------------------*

      USE mod_phys_lmdz_para, ONLY: is_master
      USE mod_phys_lmdz_transfert_para, ONLY: bcast

      IMPLICIT NONE

!     input

      integer :: nw               ! number of wavelength grid points
      real, dimension(nw) :: wl   ! lower and central wavelength for each interval

!     output

      real, dimension(nw) :: yg   ! h2o cross-sections (cm2)

!     local

      integer, parameter :: kdata = 1500
      real, parameter :: deltax = 1.e-4
      REAL x1(kdata)
      REAL y1(kdata)
      INTEGER ierr, dummy
      INTEGER i, n, iopt
      CHARACTER*100 fil
      integer :: kin, kout ! input/output logical units

      kin = 10

!---------------------------------------------------------------------------------
!     iopt = 1       Frank Mills' best estimate (1998)
!     iopt = 2       Benjamin Frandsen calculation (personal communication, 2024)
!---------------------------------------------------------------------------------

      iopt = 2

      if (iopt == 1) then
         fil = 'cross_sections/s2_millsBestEst_1560_5059.txt'
      else if (iopt == 2) then
         fil = 'cross_sections/s2_cross_sections_frandsen.txt'
      end if

      if(is_master) then

      print*, 'S2 cross-section : ', fil
      OPEN(UNIT=kin,FILE=fil,STATUS='old') 
      
      if (iopt == 1) then
         n = 1203 
         do i = 1, n
            read(kin,*) x1(i), y1(i)
            x1(i) = x1(i)/10.
         end do
      else if (iopt == 2) then
         n = 406 
         do i = 1, 4
            read(kin,*)
         end do
         do i = 1, n
            read(kin,*) x1(i), y1(i)
         end do
      end if

      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      ENDIF

      endif !is_master

      call bcast(yg)

      end subroutine rdxss2

!==============================================================================

      subroutine rdxsh2o(nw, wl, yg)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read H2O molecular absorption cross section.  Re-grid data to match      =*
!=  specified wavelength working grid.                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  YG     - REAL, molecular absoprtion cross section (cm^2) of H2O at    (O)=*
!=           each specified wavelength                                       =*
!-----------------------------------------------------------------------------*

      USE mod_phys_lmdz_para, ONLY: is_master
      USE mod_phys_lmdz_transfert_para, ONLY: bcast

      IMPLICIT NONE

!     input

      integer :: nw               ! number of wavelength grid points
      real, dimension(nw) :: wl   ! lower and central wavelength for each interval

!     output

      real, dimension(nw) :: yg   ! h2o cross-sections (cm2)

!     local

      integer, parameter :: kdata = 500
      real, parameter :: deltax = 1.e-4
      REAL x1(kdata)
      REAL y1(kdata)
      INTEGER ierr
      INTEGER i, n
      CHARACTER*100 fil
      integer :: kin, kout ! input/output logical units

      kin = 10

      fil = 'cross_sections/h2o_composite_250K.txt'
      print*, 'section efficace H2O: ', fil

      if(is_master) then

      OPEN(UNIT=kin,FILE=fil,STATUS='old') 

      DO i = 1,26 
         read(kin,*)
      END DO

      n = 420 
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      END DO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      ENDIF

      endif !is_master

      call bcast(yg)
      
      end subroutine rdxsh2o

!==============================================================================

       subroutine rdxshdo(nw, wl, yg)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read HDO molecular absorption cross section.  Re-grid data to match      =*
!=  specified wavelength working grid.                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  YG     - REAL, molecular absoprtion cross section (cm^2) of HDO at    (O)=*
!=           each specified wavelength                                       =*
!-----------------------------------------------------------------------------*

      USE mod_phys_lmdz_para, ONLY: is_master
      USE mod_phys_lmdz_transfert_para, ONLY: bcast

      IMPLICIT NONE

!     input

      integer :: nw               ! number of wavelength grid points
      real, dimension(nw) :: wl   ! lower and central wavelength for each interval

!     output

      real, dimension(nw) :: yg   ! hdo cross-sections (cm2)

!     local

      integer, parameter :: kdata = 900
      real, parameter :: deltax = 1.e-4
      REAL x1(kdata)
      REAL y1(kdata)
      INTEGER ierr
      INTEGER i, n
      CHARACTER*100 fil
      integer :: kin, kout ! input/output logical units

      kin = 10

      fil = 'cross_sections/hdo_composite_295K.txt'
      print*, 'section efficace HDO: ', fil

      if(is_master) then

      OPEN(UNIT=kin,FILE=fil,STATUS='old') 

      DO i = 1,17 
         read(kin,*)
      END DO

      n = 806 
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      END DO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      ENDIF

      endif !is_master

      call bcast(yg)
      
      end subroutine rdxshdo

!==============================================================================

      subroutine rdxsh2o2(nw, wl, xsh2o2)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read and grid H2O2 cross-sections
!=         H2O2 + hv -> 2 OH                                                 =*
!=  Cross section:  Schuergers and Welge, Z. Naturforsch. 23a (1968) 1508    =*
!=                  from 125 to 185 nm, then JPL97 from 190 to 350 nm.       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!-----------------------------------------------------------------------------*

      USE mod_phys_lmdz_para, ONLY: is_master
      USE mod_phys_lmdz_transfert_para, ONLY: bcast

      implicit none

!     input

      integer :: nw               ! number of wavelength grid points
      real, dimension(nw) :: wl   ! lower wavelength for each interval

!     output

      real, dimension(nw) :: xsh2o2   ! h2o2 cross-sections (cm2)

!     local

      real, parameter :: deltax = 1.e-4
      integer, parameter :: kdata = 100
      real, dimension(kdata) :: x1, y1
      real, dimension(nw)    :: yg
      integer :: i, ierr, iw, n, idum
      integer :: kin, kout ! input/output logical units
      character*100 fil

      kin = 10

!     read cross-sections

      fil = 'cross_sections/h2o2_composite.txt'
      print*, 'section efficace H2O2: ', fil

      if(is_master) then

      OPEN(kin,FILE=fil,STATUS='OLD')
      READ(kin,*) idum,n
      DO i = 1, idum-2
         READ(kin,*)
      ENDDO
      DO i = 1, n
         READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,               0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      ENDIF

      DO iw = 1, nw - 1
         xsh2o2(iw) = yg(iw)
      END DO

      endif !is_master

      call bcast(xsh2o2)

      end subroutine rdxsh2o2

!==============================================================================

      subroutine rdxsho2(nw, wl, yg)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read ho2 cross-sections                                                  =*
!=  JPL 2006 recommendation                                                  =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!-----------------------------------------------------------------------------*

      USE mod_phys_lmdz_para, ONLY: is_master
      USE mod_phys_lmdz_transfert_para, ONLY: bcast

      IMPLICIT NONE

!     input

      integer :: nw               ! number of wavelength grid points
      real, dimension(nw) :: wl   ! lower wavelength for each interval

!     output

      real, dimension(nw) :: yg   ! ho2 cross-sections (cm2)

!     local

      real, parameter :: deltax = 1.e-4
      integer, parameter :: kdata = 100
      real, dimension(kdata) :: x1, y1
      integer :: i, n, ierr
      character*100 fil
      integer :: kin, kout ! input/output logical units

      kin = 10

!*** cross sections from Sander et al. [2003]

      fil = 'cross_sections/ho2_jpl2003.txt'
      print*, 'section efficace HO2: ', fil

      if(is_master) then

      OPEN(kin,FILE=fil,STATUS='OLD')
      READ(kin,*) n
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
  
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, fil
        STOP
      ENDIF

      endif !is_master

      call bcast(yg)
 
      end subroutine rdxsho2

!==============================================================================

      subroutine rdxshcl(nw, wl, yg)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read hcl cross-sections                                                  =*
!=  JPL 2006 recommendation                                                  =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!-----------------------------------------------------------------------------*

      USE mod_phys_lmdz_para, ONLY: is_master
      USE mod_phys_lmdz_transfert_para, ONLY: bcast

      IMPLICIT NONE

!     input

      integer :: nw               ! number of wavelength grid points
      real, dimension(nw) :: wl   ! lower wavelength for each interval

!     output

      real, dimension(nw) :: yg   ! hcl cross-sections (cm2)

!     local

      real, parameter :: deltax = 1.e-4
      integer, parameter :: kdata = 100
      real, dimension(kdata) :: x1, y1
      integer :: i, n, ierr
      character*100 fil
      integer :: kin, kout ! input/output logical units

      kin = 10

!*** cross sections from JPL [2006]

      fil = 'cross_sections/hcl_jpl2006.txt'
      print*, 'section efficace HCl: ', fil

      if(is_master) then
      n = 31
      OPEN(kin,FILE=fil,STATUS='OLD')
      DO i = 1,4
         READ(kin,*)
      ENDDO   
      DO i = 1,n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)

  
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, fil
        STOP
      ENDIF

      endif !is_master

      call bcast(yg)
 
      end subroutine rdxshcl
 
!==============================================================================

      subroutine rdxsdcl(nw, wl, yg)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read dcl cross-sections                                                  =*
!=  JPL 2019 recommendation                                                  =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!-----------------------------------------------------------------------------*

      USE mod_phys_lmdz_para, ONLY: is_master
      USE mod_phys_lmdz_transfert_para, ONLY: bcast

      IMPLICIT NONE

!     input

      integer :: nw               ! number of wavelength grid points
      real, dimension(nw) :: wl   ! lower wavelength for each interval

!     output

      real, dimension(nw) :: yg   ! dcl cross-sections (cm2)

!     local

      real, parameter :: deltax = 1.e-4
      integer, parameter :: kdata = 100
      real, dimension(kdata) :: x1, y1
      integer :: i, n, ierr
      character*100 fil
      integer :: kin, kout ! input/output logical units

      kin = 10

!*** cross sections from JPL [2019]

      fil = 'cross_sections/dcl_jpl2019_298K.txt'
      print*, 'section efficace DCl: ', fil

      if(is_master) then
      n = 29
      OPEN(kin,FILE=fil,STATUS='OLD')
      DO i = 1,4
         READ(kin,*)
      ENDDO   
      DO i = 1,n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)

  
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, fil
        STOP
      ENDIF

      endif !is_master

      call bcast(yg)
 
      end subroutine rdxsdcl

!==============================================================================

      subroutine rdxscl2(nw, wl, yg)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read cl2 cross-sections                                                  =*
!=  JPL 2006 recommendation                                                  =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!-----------------------------------------------------------------------------*

      USE mod_phys_lmdz_para, ONLY: is_master
      USE mod_phys_lmdz_transfert_para, ONLY: bcast

      IMPLICIT NONE

!     input

      integer :: nw               ! number of wavelength grid points
      real, dimension(nw) :: wl   ! lower wavelength for each interval

!     output

      real, dimension(nw) :: yg   ! cl2 cross-sections (cm2)

!     local

      real, parameter :: deltax = 1.e-4
      integer, parameter :: kdata = 100
      real, dimension(kdata) :: x1, y1
      integer :: i, n, ierr
      character*100 fil
      integer :: kin, kout ! input/output logical units

      kin = 10

!*** cross sections from JPL [2006]

      fil = 'cross_sections/cl2_jpl2006.txt'
      print*, 'section efficace Cl2: ', fil

      if(is_master) then

      n = 30 
      OPEN(kin,FILE=fil,STATUS='OLD')
      DO i = 1,4
         READ(kin,*)
      ENDDO   
      DO i = 1,n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
  
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, fil
        STOP
      ENDIF

      endif !is_master

      call bcast(yg)
 
      end subroutine rdxscl2

!==============================================================================

      subroutine rdxshocl(nw, wl, yg)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read hocl cross-sections                                                  =*
!=  JPL 2000 recommendation                                                  =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!-----------------------------------------------------------------------------*

      USE mod_phys_lmdz_para, ONLY: is_master
      USE mod_phys_lmdz_transfert_para, ONLY: bcast

      IMPLICIT NONE

!     input

      integer :: nw               ! number of wavelength grid points
      real, dimension(nw) :: wl   ! lower wavelength for each interval

!     output

      real, dimension(nw) :: yg   ! hocl cross-sections (cm2)

!     local

      real, parameter :: deltax = 1.e-4
      integer, parameter :: kdata = 200
      real, dimension(kdata) :: x1, y1
      integer :: i, n, ierr
      character*100 fil
      integer :: kin, kout ! input/output logical units

      kin = 10

!*** cross sections from JPL [2000]

      fil = 'cross_sections/HOCl_jpl2000.txt'
      print*, 'section efficace HOCl: ', fil

      if(is_master) then

         n = 111
      OPEN(kin,FILE=fil,STATUS='OLD')
      DO i = 1,6
         READ(kin,*)
      ENDDO   
      DO i = 1,111
        READ(kin,*) y1(i)
        x1(i) = 200 + real(i-1)*2.
        y1(i) = y1(i) * 1.E-20
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
  
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, fil
        STOP
      ENDIF

      endif !is_master

      call bcast(yg)
 
      end subroutine rdxshocl

!==============================================================================

      subroutine rdxsso2(nw,wl,xsso2_200,xsso2_298,xsso2_360)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read and grid SO2 absorption cross-sections and photodissociation yield   =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  XSSO2  - REAL, molecular absoprtion cross section (cm^2) of SO2 at    (O)=*
!=           each specified wavelength                                       =*
!-----------------------------------------------------------------------------*

      USE mod_phys_lmdz_para, ONLY: is_master
      USE mod_phys_lmdz_transfert_para, ONLY: bcast

      implicit none

!     input

      integer :: nw               ! number of wavelength grid points
      real, dimension(nw) :: wl   ! lower and central wavelength for each interval

!     output

      real, dimension(nw) :: xsso2_200, xsso2_298, xsso2_360 ! so2 cross-sections (cm2)

!     local

      integer, parameter :: kdata = 55000
      real, parameter :: deltax = 1.e-4
      real, dimension(kdata) :: x1, y1, y2, y3, xion, ion
      real, dimension(nw) :: yg
      real :: xl, xu
      integer :: ierr, i, l, n, n1, n2, n3, n4
      CHARACTER*100 fil

      integer :: kin, kout ! input/ouput logical units

      kin  = 10
      kout = 30

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     SO2 absorption cross-sections
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     iopt = 1
!
!     200K: Wu et al. (2000) + Vandaele et al. (2009) + Hermans et al. (2009)
!           7 lignes d'en-tete et n1 = 50276 
!           fichier: so2_composite_200K.txt 
!
!     298K: Wu et al. (2000) + Vandaele et al. (2009) + Hermans et al. (2009)
!           7 lignes d'en-tete et n2 = 49833 
!           fichier: so2_composite_298K.txt 
!
!     360K: Wu et al. (2000) + Vandaele et al. (2009) + Hermans et al. (2009)
!           7 lignes d'en-tete et n3 = 49261 
!           fichier: so2_composite_360K.txt 
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      n1 = 50276
      n2 = 49833
      n3 = 46261

!     200K:

      fil = 'cross_sections/so2_composite_200K.txt'
      print*, 'section efficace SO2 195K: ', fil

      if(is_master) then

      OPEN(UNIT=kin,FILE=fil,STATUS='old')
      DO i = 1,7 
         read(kin,*)
      END DO

      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
      END DO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,          0.,0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,      1.e+38,0.)
      CALL inter2(nw,wl,yg,n1,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      ENDIF
      
      DO l = 1, nw-1
         xsso2_200(l) = yg(l)
      END DO

      endif !is_master

      call bcast(xsso2_200)

!     298K:

      fil = 'cross_sections/so2_composite_298K.txt'
      print*, 'section efficace SO2 295K: ', fil

      if(is_master) then

      OPEN(UNIT=kin,FILE=fil,STATUS='old')
      DO i = 1,7
         read(kin,*)
      END DO

      DO i = 1, n2
         READ(kin,*) x1(i), y1(i)
      END DO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n2,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n2,          0.,0.)
      CALL addpnt(x1,y1,kdata,n2,x1(n2)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n2,      1.e+38,0.)
      CALL inter2(nw,wl,yg,n2,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      ENDIF

      DO l = 1, nw-1
         xsso2_298(l) = yg(l)
      END DO

      endif !is_master

      call bcast(xsso2_298)

!     360K:

      fil = 'cross_sections/so2_composite_360K.txt'
      print*, 'section efficace SO2 370K: ', fil

      if(is_master) then

      OPEN(UNIT=kin,FILE=fil,STATUS='old')
      DO i = 1,7
         read(kin,*)
      END DO

      DO i = 1, n3
         READ(kin,*) x1(i), y1(i)
      END DO
      CLOSE (kin)

      CALL addpnt(x1,y1,kdata,n3,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n3,          0.,0.)
      CALL addpnt(x1,y1,kdata,n3,x1(n3)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n3,      1.e+38,0.)
      CALL inter2(nw,wl,yg,n3,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      ENDIF

      DO l = 1, nw-1
         xsso2_360(l) = yg(l)
      END DO

      endif !is_master

      call bcast(xsso2_360)

!     DO l = 1, nw-1
!        write(kout,*) wl(l), xsso2_200(l),
!    $                        xsso2_298(l),
!    $                        xsso2_360(l),
!     END DO

      end subroutine rdxsso2

!==============================================================================

      subroutine rdxsso(nw,wl,xsso_150,xsso_250)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read SO cross-sections                                                   =*
!=  Phillips (1981) or Heays et al. (2023)                                   =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!-----------------------------------------------------------------------------*

      USE mod_phys_lmdz_para, ONLY: is_master
      USE mod_phys_lmdz_transfert_para, ONLY: bcast

      implicit none

!     input

      integer :: nw               ! number of wavelength grid points
      real, dimension(nw) :: wl   ! lower wavelength for each interval

!     output

      real, dimension(nw) :: xsso_150, xsso_250 ! so cross-sections (cm2)

!     local

      integer, parameter :: kdata = 1000
      integer :: i, l, n, n1, n2, ierr, iopt
      integer :: kin, kout ! input/output logical units
      real, dimension(nw) :: yg, yg1, yg2
      real, dimension(kdata) :: x1, y1, x2, y2
      real :: dummy
      real, parameter :: deltax = 1.e-4
      character*100 fil

      kin = 10

!     iopt = 1 : Phillips (1981)                        
!     iopt = 2 : Heays et al., Molecular Physics, 2023 

      iopt = 2

      if (iopt == 1) then

!        cross sections from Phillips (1981)

         fil = 'cross_sections/so_marcq.txt'

         if (is_master) then

         print*, 'section efficace SO: ', fil

         n = 851
         OPEN(kin,FILE=fil,STATUS='OLD') 
         DO i = 1,n
            READ(kin,*) x1(i), dummy, dummy, dummy, dummy, y1(i)
            y1(i) = y1(i)*0.88  ! from Mills, PhD, 1998
         ENDDO
         CLOSE(kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,          0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,        1E38,0.)

         do i = 1,n
            print*, x1(i), y1(i)
         end do

         CALL inter2(nw,wl,yg,n,x1,y1,ierr)
  
         IF (ierr .NE. 0) THEN
           WRITE(*,*) ierr, fil
           STOP
         ENDIF

         DO l = 1, nw-1
            xsso_150(l) = yg(l)
            xsso_250(l) = yg(l)
         END DO

         endif !is_master

      else if (iopt == 2) then

!        cross sections from Heays et al. (2023)
!        values given at 150 K and 250 K

         fil = 'cross_sections/SO_heays_2023.txt'

         if (is_master) then

         print*, 'section efficace SO: ', fil
         
         n1 = 525
         n2 = n1
         OPEN(kin,FILE=fil,STATUS='OLD') 

         DO i = 1, 3
            READ(kin,*)
         ENDDO

         DO i = 1,n1
            READ(kin,*) x1(i), dummy, dummy, dummy, dummy, dummy, &
                        dummy, dummy, dummy, y1(i), y2(i)
            x2(i) = x1(i)
         ENDDO
         CLOSE(kin)

         CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,               0.,0.)
         CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n1,           1.e+38,0.)

         CALL inter2(nw,wl,yg1,n1,x1,y1,ierr)

         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, fil
            STOP
         ENDIF

         CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
         CALL addpnt(x2,y2,kdata,n2,               0.,0.)
         CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
         CALL addpnt(x2,y2,kdata,n2,           1.e+38,0.)

         CALL inter2(nw,wl,yg2,n2,x2,y2,ierr)

         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, fil
            STOP
         ENDIF

         DO l = 1, nw-1
            xsso_150(l) = yg1(l)
            xsso_250(l) = yg2(l)
         END DO

         endif !is_master

      end if ! iopt

      call bcast(xsso_150)
      call bcast(xsso_250)
 
      end subroutine rdxsso

!==============================================================================

      subroutine rdxsso3(nw, wl, yg)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read SO3 cross-sections                                                  =*
!=  Composite section 140-294 nm -> Pintze al. [2003]                        =*
!=                    296-330 nm -> Burkholder al. [1997]                    =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!-----------------------------------------------------------------------------*

      USE mod_phys_lmdz_para, ONLY: is_master
      USE mod_phys_lmdz_transfert_para, ONLY: bcast

      IMPLICIT NONE

!     input

      integer :: nw               ! number of wavelength grid points
      real, dimension(nw) :: wl   ! lower wavelength for each interval

!     output

      real, dimension(nw) :: yg   ! SO3 cross-sections (cm2)

!     local

      real, parameter :: deltax = 1.e-4
      integer, parameter :: kdata = 200
      real, dimension(kdata) :: x1, y1
      integer :: i, n, ierr
      character*100 fil
      integer :: kin, kout ! input/output logical units

      kin = 10

!*** cross sections

      fil = 'cross_sections/so3_composite.txt'
      print*, 'section efficace SO3: ', fil

      if(is_master) then

         n = 96
      OPEN(kin,FILE=fil,STATUS='OLD')
      DO i = 1,16
         READ(kin,*)
      ENDDO   
      DO i = 1,n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
  
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, fil
        STOP
      ENDIF

      endif !is_master

      call bcast(yg)
 
      end subroutine rdxsso3

!==============================================================================

      subroutine rdxsosso(nw,wl,wc,yg1,yg2,yg3)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  read SO dimer absorption cross-sections                                  =*
!=                                                                           =*
!=           OSSO_cis + hv -> SO + SO                                        =*
!=           OSSO_trans + hv -> SO + SO                                      =*
!=           S2O2_cyc + hv -> SO + SO                                        =*
!=                                                                           =*
!=  Theoritical Cross section: from Egan et al. GRL, 2025                    =*
!=  Quantum yield: assumed to be unity                                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!-----------------------------------------------------------------------------*
!=  EDIT HISTORY:                                                            =*
!=  2025/02/14  Original, Franck Lefevre                                     =*
!-----------------------------------------------------------------------------*

      USE mod_phys_lmdz_para, ONLY: is_master
      USE mod_phys_lmdz_transfert_para, ONLY: bcast

      implicit none

!     input

      integer :: nw                   ! number of wavelength grid points
      real, dimension(nw) :: wl, wc   ! lower and central wavelength for each interval

!     output

      real, dimension(nw) :: yg1  ! osso_cis cross-sections (cm2)
      real, dimension(nw) :: yg2  ! osso_trans cross-sections (cm2)
      real, dimension(nw) :: yg3  ! s2o2_cyc cross-sections (cm2)

!     local

      integer, parameter :: kdata = 1000
      integer :: i, iw, n, n1, n2, n3, ierr
      integer :: kin
      real, parameter :: deltax = 1.e-4
      real, dimension(kdata) :: x1, y1, x2, y2, x3, y3
      real :: qy, lambda
      character(len = 100) :: fil

!**** cross sections from Egan et al., 2025

      fil = 'cross_sections/OSSO_cross_sections_formatted.txt'
      print*, 'OSSO cross-section: ', fil

      if (is_master) then

      kin = 10
      open(kin, file=fil, status='old', iostat=ierr)

      n = 600
      do i = 1,9
         read(kin,*)
      end do
      do i = 1,n
         read(kin,*) x1(i), y1(i), x2(i), y2(i), x3(i), y3(i)
      end do
      close(kin)

      n1 = n
      call addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
      call addpnt(x1,y1,kdata,n1,          0.,0.)
      call addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      call addpnt(x1,y1,kdata,n1,        1E38,0.)
      call inter2(nw,wl,yg1,n1,x1,y1,ierr) 

      n2 = n
      call addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
      call addpnt(x2,y2,kdata,n2,          0.,0.)
      call addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      call addpnt(x2,y2,kdata,n2,        1E38,0.)
      call inter2(nw,wl,yg2,n2,x2,y2,ierr) 
     
      n3 = n
      call addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),0.)
      call addpnt(x3,y3,kdata,n3,          0.,0.)
      call addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
      call addpnt(x3,y3,kdata,n3,        1E38,0.)
      call inter2(nw,wl,yg3,n3,x3,y3,ierr) 

      if (ierr /= 0) THEN
         write(*,*) ierr, fil
         stop
      end if

!**** quantum yield assumed to be unity

      qy = 1.
     
      do iw = 1,nw - 1
         lambda = wc(iw)
         if (lambda < 842.) then
            yg1(iw) = qy*yg1(iw)
         else
            yg1(iw) = 0.
         end if
         if (lambda < 938.) then
            yg2(iw) = qy*yg2(iw)
         else
            yg2(iw) = 0.
         end if
         yg3(iw) = qy*yg3(iw)
!        write(55,'(f8.3,3e12.4)') lambda, yg1(iw), yg2(iw), yg3(iw)
      end do

      end if ! is_master

      call bcast(yg1)
      call bcast(yg2)
      call bcast(yg3)

      end subroutine rdxsosso

!==============================================================================

      subroutine rdxsclso2(nw, wl, wc, yg)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read clso2 cross-sections                                                =*
!=  From Trabelsi et al., J. Chem. Phys. 161, 109901 (2024)                  =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!-----------------------------------------------------------------------------*

      USE mod_phys_lmdz_para, ONLY: is_master
      USE mod_phys_lmdz_transfert_para, ONLY: bcast

      IMPLICIT NONE

!     input

      integer :: nw                   ! number of wavelength grid points
      real, dimension(nw) :: wl, wc   ! lower and central wavelength for each interval

!     output

      real, dimension(nw) :: yg   ! clso2 cross-sections (cm2)

!     local

      real, parameter :: deltax = 1.e-4
      integer, parameter :: kdata = 1500

      real, dimension(kdata) :: x1, y1
      real :: qy, lambda
      integer :: i, iw, n, ierr
      integer :: kin, kout ! input/output logical units
      character*100 fil

      kin = 10

!*** cross sections from Trabelsi et al., J. Chem. Phys. 161, 109901 (2024)

      fil = 'cross_sections/clso2_trabelsi_2024.txt'
      print*, 'section efficace ClSO2: ', fil

      if (is_master) then

      n = 1001 

      OPEN(kin,FILE=fil,STATUS='OLD')

      do i = 1,5
         read(kin,*)
      end do   
      do i = n,1,-1
        read(kin,*) x1(i), y1(i)
      end do
      close(kin)

      call addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      call addpnt(x1,y1,kdata,n,          0.,0.)
      call addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      call addpnt(x1,y1,kdata,n,        1e38,0.)

      call inter2(nw,wl,yg,n,x1,y1,ierr)
  
      if (ierr /= 0) then
        write(*,*) ierr, fil
        stop
      end if

      do iw = 1,nw - 1
         lambda = wc(iw)
         if (lambda > 316.) then  ! photodissociation threshold (tbc)
            yg(iw) = 0.
         end if
!        write(55,'(f8.3,3e12.4)') lambda, yg(iw)
      end do

      end if !is_master

      call bcast(yg)
 
      end subroutine rdxsclso2

!==============================================================================

      subroutine rdxscl2so2(nw, wl, wc, yg)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read cl2so2 cross-sections                                               =*
!=  From Uthman et al., J. Phys. Chem., 1978                                 =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!-----------------------------------------------------------------------------*

      USE mod_phys_lmdz_para, ONLY: is_master
      USE mod_phys_lmdz_transfert_para, ONLY: bcast

      IMPLICIT NONE

!     input

      integer :: nw                   ! number of wavelength grid points
      real, dimension(nw) :: wl, wc   ! lower and central wavelength for each interval

!     output

      real, dimension(nw) :: yg   ! clso2 cross-sections (cm2)

!     local

      real, parameter :: deltax = 1.e-4
      integer, parameter :: kdata = 1500

      real, dimension(kdata) :: x1, y1
      real :: qy, lambda
      integer :: i, iw, n, ierr
      integer :: kin, kout ! input/output logical units
      character*100 fil

      kin = 10

!*** cross sections from Uthman et al., J. Phys. Chem., 1978 

      fil = 'cross_sections/cl2so2_Uthman1978_295K_190_300nm.txt'
      print*, 'section efficace Cl2SO2: ', fil

      if (is_master) then

      n = 56 

      OPEN(kin,FILE=fil,STATUS='OLD')

      do i = 1,7
         read(kin,*)
      end do   
      do i = 1,n
        read(kin,*) x1(i), y1(i)
      end do
      close(kin)

      call addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      call addpnt(x1,y1,kdata,n,          0.,0.)
      call addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      call addpnt(x1,y1,kdata,n,        1e38,0.)

      call inter2(nw,wl,yg,n,x1,y1,ierr)
  
      if (ierr /= 0) then
        write(*,*) ierr, fil
        stop
      end if

      do iw = 1,nw - 1
         lambda = wc(iw)
!        write(55,'(f8.3,3e12.4)') lambda, yg(iw)
      end do

      end if !is_master

      call bcast(yg)
 
      end subroutine rdxscl2so2

!==============================================================================

      subroutine rdxsclo(nw, wl, yg)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read ClO cross-sections                                                  =*
!=  From Trollier and al. [1990]                                             =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!-----------------------------------------------------------------------------*

      USE mod_phys_lmdz_para, ONLY: is_master
      USE mod_phys_lmdz_transfert_para, ONLY: bcast

      IMPLICIT NONE

!     input

      integer :: nw               ! number of wavelength grid points
      real, dimension(nw) :: wl   ! lower wavelength for each interval

!     output

      real, dimension(nw) :: yg   ! hcl cross-sections (cm2)

!     local

      real, parameter :: deltax = 1.e-4
      integer, parameter :: kdata = 100
      real, dimension(kdata) :: x1, y1
      integer :: i, n, ierr
      character*100 fil
      integer :: kin, kout ! input/output logical units

      kin = 10

!*** cross sections

      fil = 'cross_sections/clo_xs_trolier_1990.txt'
      print*, 'section efficace ClO: ', fil

      if(is_master) then

         n = 81
      OPEN(kin,FILE=fil,STATUS='OLD')
      DO i = 1,6
         READ(kin,*)
      ENDDO   
      DO i = 1,n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
  
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, fil
        STOP
      ENDIF

      endif !is_master

      call bcast(yg)
 
      end subroutine rdxsclo

!==============================================================================

      subroutine rdxsocs(nw, wl, yg)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read OCS cross-sections                                                  =*
!=  JPL 2011 recommendation                                                  =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!-----------------------------------------------------------------------------*

      USE mod_phys_lmdz_para, ONLY: is_master
      USE mod_phys_lmdz_transfert_para, ONLY: bcast

      IMPLICIT NONE

!     input

      integer :: nw               ! number of wavelength grid points
      real, dimension(nw) :: wl   ! lower wavelength for each interval

!     output

      real, dimension(nw) :: yg   ! hcl cross-sections (cm2)

!     local

      real, parameter :: deltax = 1.e-4
      integer, parameter :: kdata = 100
      real, dimension(kdata) :: x1, y1
      integer :: i, n, ierr
      character*100 fil
      integer :: kin, kout ! input/output logical units

      kin = 10

!*** cross sections from JPL [2006]

      fil = 'cross_sections/ocs_cross_sections_jpl2011.txt'
      print*, 'section efficace OCS: ', fil

      if(is_master) then

         n = 40
      OPEN(kin,FILE=fil,STATUS='OLD')
      DO i = 1,4
         READ(kin,*)
      ENDDO   
      DO i = 1,n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
  
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, fil
        STOP
      ENDIF

      endif !is_master

      call bcast(yg)
 
      end subroutine rdxsocs

!==============================================================================

      subroutine rdxscocl2(nw, wl, yg)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read COCl2 cross-sections                                                =*
!=  JPL 2011 recommendation                                                  =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!-----------------------------------------------------------------------------*

      USE mod_phys_lmdz_para, ONLY: is_master
      USE mod_phys_lmdz_transfert_para, ONLY: bcast

      IMPLICIT NONE

!     input

      integer :: nw               ! number of wavelength grid points
      real, dimension(nw) :: wl   ! lower wavelength for each interval

!     output

      real, dimension(nw) :: yg   ! COCl2 cross-sections (cm2)

!     local

      real, parameter :: deltax = 1.e-4
      integer, parameter :: kdata = 100
      real, dimension(kdata) :: x1, y1
      integer :: i, n, ierr
      character*100 fil
      integer :: kin, kout ! input/output logical units

      kin = 10

!*** cross sections from JPL [2011]

      fil = 'cross_sections/cocl2_cross_sections_jpl2011.txt'
      print*, 'section efficace COCl2: ', fil

      if(is_master) then

         n = 53
      OPEN(kin,FILE=fil,STATUS='OLD')
      DO i = 1,4
         READ(kin,*)
      ENDDO   
      DO i = 1,n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
  
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, fil
        STOP
      ENDIF

      endif !is_master

      call bcast(yg)
 
      end subroutine rdxscocl2

!==============================================================================

      subroutine rdxsh2so4(nw, wl, yg)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read H2SO4 cross-sections                                                  =*
!=  JPL 2006 recommendation                                                  =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!-----------------------------------------------------------------------------*

      USE mod_phys_lmdz_para, ONLY: is_master
      USE mod_phys_lmdz_transfert_para, ONLY: bcast

      IMPLICIT NONE

!     input

      integer :: nw               ! number of wavelength grid points
      real, dimension(nw) :: wl   ! lower wavelength for each interval

!     output

      real, dimension(nw) :: yg   ! h2so4 cross-sections (cm2)

!     local

      real, parameter :: deltax = 1.e-4
      integer, parameter :: kdata = 100
      real, dimension(kdata) :: x1, y1
      integer :: i, n, ierr
      character*100 fil
      integer :: kin, kout ! input/output logical units

      kin = 10

!*** cross sections from JPL [2006]

      fil = 'cross_sections/h2so4_cross_sections.txt'
      print*, 'section efficace h2so4: ', fil

      if(is_master) then

         n = 22
      OPEN(kin,FILE=fil,STATUS='OLD')
      DO i = 1,3
         READ(kin,*)
      ENDDO   
      DO i = 1,n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
  
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, fil
        STOP
      ENDIF

      endif !is_master

      call bcast(yg)
 
      end subroutine rdxsh2so4

!==============================================================================

       subroutine rdxsh2(nw, wl, yg, yieldh2)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read h2 cross-sections and photodissociation yield                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!-----------------------------------------------------------------------------*

      USE mod_phys_lmdz_para, ONLY: is_master
      USE mod_phys_lmdz_transfert_para, ONLY: bcast

      implicit none

!     input

      integer :: nw                   ! number of wavelength grid points
      real, dimension(nw) :: wl       ! lower and central wavelength for each interval

!     output

      real, dimension(nw) :: yg        ! h2 cross-sections (cm2)
      real, dimension(nw) :: yieldh2   ! photodissociation yield

!     local

      integer, parameter :: kdata = 1000
      real, parameter :: deltax = 1.e-4
      real, dimension(kdata) :: x1, y1, x2, y2
      real :: xl, xu
      integer :: i, iw, n, ierr
      integer :: kin, kout ! input/output logical units
      character*100 fil

      kin = 10

!     h2 cross sections

      fil = 'cross_sections/h2secef.txt'
      print*, 'section efficace H2: ', fil

      if(is_master) then

      OPEN(kin,FILE=fil,STATUS='OLD')

      n = 792 
      read(kin,*) ! avoid first line with wavelength = 0.
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
      ENDDO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)
      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
  
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, fil
        STOP
      ENDIF

      endif !is_master

      call bcast(yg)
 
!     photodissociation yield

      fil = 'cross_sections/h2_ionef_schunknagy2000.txt'

      if(is_master) then

      OPEN(UNIT=kin,FILE=fil,STATUS='old')

      n = 19
      read(kin,*)
      DO i = 1, n
         READ(kin,*) xl, xu, y2(i)
         x2(i) = (xl + xu)/2.
         y2(i) = max(1. - y2(i),0.)
      END DO
      CLOSE (kin)

      CALL addpnt(x2,y2,kdata,n,x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n,          0.,0.)
      CALL addpnt(x2,y2,kdata,n,x2(n)*(1.+deltax),1.)
      CALL addpnt(x2,y2,kdata,n,      1.e+38,1.)
      CALL inter2(nw,wl,yieldh2,n,x2,y2,ierr)
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, fil
        STOP
      ENDIF

      endif !is_master

      call bcast(yieldh2)

      end subroutine rdxsh2

!==============================================================================

      subroutine rdxsno2(nw,wl,xsno2,xsno2_220,xsno2_294,yldno2_248, yldno2_298)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  read and grid cross section + quantum yield for NO2                      =*
!=  photolysis                                                               =*
!=  Jenouvrier et al., 1996  200-238 nm
!=  Vandaele et al., 1998    238-666 nm 220K and 294K 
!=  quantum yield from jpl 2006
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
!=           photolysis reaction defined, at each defined wavelength and     =*
!=           at each defined altitude level                                  =*
!-----------------------------------------------------------------------------*

      USE mod_phys_lmdz_para, ONLY: is_master
      USE mod_phys_lmdz_transfert_para, ONLY: bcast

      implicit none

!     input

      integer :: nw               ! number of wavelength grid points
      real, dimension(nw) :: wl   ! lower and central wavelength for each interval
      
!     output

      real, dimension(nw) :: xsno2, xsno2_220, xsno2_294 ! no2 cross-sections (cm2)
      real, dimension(nw) :: yldno2_248, yldno2_298      ! quantum yields at 248-298 k

!     local

      integer, parameter :: kdata = 28000
      real, parameter :: deltax = 1.e-4
      real, dimension(kdata) :: x1, x2, x3, x4, x5, y1, y2, y3, y4, y5
      real, dimension(nw) :: yg1, yg2, yg3, yg4, yg5
      real :: dum, qy
      integer :: i, iw, n, n1, n2, n3, n4, n5, ierr
      character*100 fil
      integer :: kin, kout ! input/output logical units

      kin = 10

!*************** NO2 photodissociation

!  Jenouvrier 1996 + Vandaele 1998 (JPL 2006)  

      fil = 'cross_sections/no2_xs_jenouvrier.txt'
      print*, 'section efficace NO2: ', fil

      if (is_master) then

      OPEN(UNIT=kin,FILE=fil,status='old')
      DO i = 1, 3
         READ(kin,*)
      END DO
      n1 = 10001
      DO i = 1, n1
         READ(kin,*) x1(i), y1(i)
      end do

      CALL addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax), 0.)
      CALL addpnt(x1,y1,kdata,n1,               0., 0.)
      CALL addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n1,           1.e+38, 0.)
      CALL inter2(nw,wl,yg1,n1,x1,y1,ierr)

      end if !is_master

      call bcast(yg1)

      fil = 'cross_sections/no2_xs_vandaele_294K.txt'
      print*, 'section efficace NO2: ', fil

      if (is_master) then

      OPEN(UNIT=kin,FILE=fil,status='old')
      DO i = 1, 3
        READ(kin,*)
      END DO
      n2 = 27993
      DO i = 1, n2
         READ(kin,*) x2(i), y2(i)
      end do

      CALL addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax), 0.)
      CALL addpnt(x2,y2,kdata,n2,               0., 0.)
      CALL addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      CALL addpnt(x2,y2,kdata,n2,           1.e+38, 0.)
      CALL inter2(nw,wl,yg2,n2,x2,y2,ierr)

      end if !is_master

      call bcast(yg2)

      fil = 'cross_sections/no2_xs_vandaele_220K.txt'
      print*, 'section efficace NO2: ', fil

      if (is_master) then

      OPEN(UNIT=kin,FILE=fil,status='old')
      DO i = 1, 3
         READ(kin,*)
      END DO
      n3 = 27993
      do i = 1, n3
         READ(kin,*) x3(i), y3(i)
      end do

      CALL addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax), 0.)
      CALL addpnt(x3,y3,kdata,n3,               0., 0.)
      CALL addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
      CALL addpnt(x3,y3,kdata,n3,           1.e+38, 0.)
      CALL inter2(nw,wl,yg3,n3,x3,y3,ierr)

      do iw = 1, nw - 1
         xsno2(iw)     = yg1(iw)
         xsno2_294(iw) = yg2(iw)
         xsno2_220(iw) = yg3(iw)
      end do

      end if !is_master

      call bcast(yg3)
      call bcast(xsno2)
      call bcast(xsno2_294)
      call bcast(xsno2_220)

!     photodissociation efficiency from jpl 2006

      fil = 'cross_sections/no2_yield_jpl2006.txt'
      print*, 'quantum yield NO2: ', fil

      if (is_master) then

      OPEN(UNIT=kin,FILE=fil,STATUS='old')
      DO i = 1, 5
         READ(kin,*) 
      END DO
      n = 25
      n4 = n
      n5 = n
      DO i = 1, n
         READ(kin,*) x4(i), y4(i), y5(i)
         x5(i) = x4(i)
      END DO
      CLOSE(kin)

      CALL addpnt(x4,y4,kdata,n4,x4(1)*(1.-deltax),y4(1))
      CALL addpnt(x4,y4,kdata,n4,               0.,y4(1))
      CALL addpnt(x4,y4,kdata,n4,x4(n4)*(1.+deltax),  0.)
      CALL addpnt(x4,y4,kdata,n4,           1.e+38,   0.)
      CALL inter2(nw,wl,yg4,n4,x4,y4,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      END IF

      end if !is_master

      call bcast(yg4)

      if (is_master) then

      CALL addpnt(x5,y5,kdata,n5,x5(1)*(1.-deltax),y5(1))
      CALL addpnt(x5,y5,kdata,n5,               0.,y5(1))
      CALL addpnt(x5,y5,kdata,n5,x5(n5)*(1.+deltax),  0.)
      CALL addpnt(x5,y5,kdata,n5,           1.e+38,   0.)
      CALL inter2(nw,wl,yg5,n5,x5,y5,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      END IF

      do iw = 1, nw - 1
         yldno2_298(iw) = yg4(iw)
         yldno2_248(iw) = yg5(iw)
      end do

      end if !is_master

      call bcast(yg5)
      call bcast(yldno2_298)
      call bcast(yldno2_248)
      
      end subroutine rdxsno2

!==============================================================================

      subroutine rdxsno(nw, wl, yg, yieldno)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read NO cross-sections  and photodissociation efficiency                 =*
!=  Lida et al 1986 (provided by Francisco Gonzalez-Galindo)                 =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!-----------------------------------------------------------------------------*

      USE mod_phys_lmdz_para, ONLY: is_master
      USE mod_phys_lmdz_transfert_para, ONLY: bcast

      implicit none

!     input

      integer :: nw               ! number of wavelength grid points
      real, dimension(nw) :: wl   ! lower wavelength for each interval

!     output

      real, dimension(nw) :: yg        ! no cross-sections (cm2)
      real, dimension(nw) :: yieldno   ! no photodissociation efficiency

!     local

      integer, parameter :: kdata = 110
      real, parameter :: deltax = 1.e-4
      real, dimension(kdata) :: x1, y1, x2, y2
      integer :: i, iw, n, ierr
      character*100 fil
      integer :: kin, kout ! input/output logical units

      kin = 10

!     no cross-sections

      fil = 'cross_sections/no_xs_francisco.txt'
      print*, 'section efficace NO: ', fil

      if (is_master) then

      OPEN(kin,FILE=fil,STATUS='OLD')

      n = 99
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
      END DO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      END IF

      end if !is_master

      call bcast(yg)
  
!     photodissociation yield

      fil = 'cross_sections/noefdis.txt'

      if (is_master) then

      OPEN(UNIT=kin,FILE=fil,STATUS='old')

      n = 33
      DO i = 1, n
         READ(kin,*) x2(n-i+1), y2(n-i+1)
      END DO
      CLOSE (kin)

      CALL addpnt(x2,y2,kdata,n,x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n,          0.,0.)
      CALL addpnt(x2,y2,kdata,n,x2(n)*(1.+deltax),1.)
      CALL addpnt(x2,y2,kdata,n,      1.e+38,1.)
      CALL inter2(nw,wl,yieldno,n,x2,y2,ierr)
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, fil
        STOP
      END IF

      end if !is_master

      call bcast(yieldno)

      end subroutine rdxsno

!==============================================================================

      subroutine rdxsn2(nw, wl, yg, yieldn2)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Read n2 cross-sections and photodissociation yield                       =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
!=           wavelength grid                                                 =*
!-----------------------------------------------------------------------------*

      USE mod_phys_lmdz_para, ONLY: is_master
      USE mod_phys_lmdz_transfert_para, ONLY: bcast

      implicit none

!     input

      integer :: nw               ! number of wavelength grid points
      real, dimension(nw) :: wl   ! lower wavelength for each interval

!     output

      real, dimension(nw) :: yg        ! n2 cross-sections (cm2)
      real, dimension(nw) :: yieldn2   ! n2 photodissociation yield

!     local

      integer, parameter :: kdata = 1100
      real, parameter :: deltax = 1.e-4
      real, dimension(kdata) :: x1, y1, x2, y2
      real :: xl, xu
      integer :: i, iw, n, ierr
      integer :: kin, kout ! input/output logical units
      character*100 fil

      kin = 10

!     n2 cross sections

      fil = 'cross_sections/n2secef_01nm.txt'
      print*, 'section efficace N2: ', fil

      if (is_master) then

      OPEN(kin,FILE=fil,STATUS='OLD')

      n = 1020
      DO i = 1, n
        READ(kin,*) x1(i), y1(i)
      END DO
      CLOSE(kin)

      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      CALL addpnt(x1,y1,kdata,n,          0.,0.)
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,        1E38,0.)

      CALL inter2(nw,wl,yg,n,x1,y1,ierr)
  
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, fil
        STOP
      END IF

      end if !is_master

      call bcast(yg)
 
!     photodissociation yield

      fil = 'cross_sections/n2_ionef_schunknagy2000.txt'

      if (is_master) then

      OPEN(UNIT=kin,FILE=fil,STATUS='old')

      n = 19
      read(kin,*)
      DO i = 1, n
         READ(kin,*) xl, xu, y2(i)
         x2(i) = (xl + xu)/2.
         y2(i) = 1. - y2(i)
      END DO
      CLOSE (kin)

      CALL addpnt(x2,y2,kdata,n,x2(1)*(1.-deltax),0.)
      CALL addpnt(x2,y2,kdata,n,          0.,0.)
      CALL addpnt(x2,y2,kdata,n,x2(n)*(1.+deltax),1.)
      CALL addpnt(x2,y2,kdata,n,      1.e+38,1.)
      CALL inter2(nw,wl,yieldn2,n,x2,y2,ierr)
      IF (ierr .NE. 0) THEN
        WRITE(*,*) ierr, fil
        STOP
      END IF

      end if !is_master

      call bcast(yieldn2)

      end subroutine rdxsn2

!==============================================================================
 
      subroutine setalb(nw,wl,albedo)
 
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Set the albedo of the surface.  The albedo is assumed to be Lambertian,  =*
!=  i.e., the reflected light is isotropic, and idependt of the direction    =*
!=  of incidence of light.  Albedo can be chosen to be wavelength dependent. =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NW      - INTEGER, number of specified intervals + 1 in working       (I)=*
!=            wavelength grid                                                =*
!=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
!=           working wavelength grid                                         =*
!=  ALBEDO  - REAL, surface albedo at each specified wavelength           (O)=*
!-----------------------------------------------------------------------------*

      implicit none

! input: (wavelength working grid data)

      INTEGER nw
      REAL wl(nw)

! output:

      REAL albedo(nw)

! local:

      INTEGER iw
      REAL alb

!     0.015: mean value from clancy et al., icarus, 49-63, 1999.

      alb = 0.015

      do iw = 1, nw - 1
         albedo(iw) = alb
      end do

      end subroutine setalb
 
end module photolysis_mod
