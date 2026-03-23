      SUBROUTINE concentrations2(pplay,t_seri,tr_seri, nqmx)

      use dimphy
      use conc,  only: mmean, rho, Akknew, rnew, cpnew
      use cpdet_phy_mod, only: cpdet                       
      USE chemparam_mod
      USE infotrac_phy, ONLY: tname
      USE clesphys_mod
      USE YOMCST_mod
      implicit none

!=======================================================================
! CALCULATION OF MEAN MOLECULAR MASS, Cp, Akk and R
!
! mmean(klon,klev)	amu
! cpnew(klon,klev)	J/kg/K
! rnew(klon,klev)	J/kg/K
! akknew(klon,klev)	coefficient of thermal conduction
!
! version: April 2012 - Franck Lefevre
!=======================================================================

!     declarations
 
c#include "YOMCST.h"
c#include "clesphys.h"
c#include "comdiurn.h"
c#include "chimiedata.h"
c#include "tracer.h"
c#include "mmol.h"

!     input/output

      real pplay(klon,klev)
      integer,intent(in) :: nqmx    ! number of tracers
      real t_seri(klon, klev)
      real tr_seri(klon,klev,nqmx)

!     local variables

      integer       :: i, l, ig, iq
      integer, save :: nbq
      integer,allocatable,save :: niq(:)
      real          :: ni(nqmx), ntot
      real          :: zt(klon, klev)
      real          :: zq(klon, klev, nqmx)
      real,allocatable,save    :: aki(:)
      real,allocatable,save    :: cpi(:)
      real, save    :: akin,akin2

      logical, save :: firstcall = .true.

      if (firstcall) then

!        initialize thermal conductivity and specific heat coefficients
!        values are taken from the literature [J/kg K]

         ! allocate local saved arrays:
         allocate(aki(nqmx))
         allocate(cpi(nqmx))
         allocate(niq(nqmx))

!        find index of chemical tracers to use
!        initialize thermal conductivity and specific heat coefficients
!        !? values are estimated

         nbq = 0 ! to count number of tracers used in this subroutine

!        aki values comes from Aeronomy part B by P.M. BANKS / G. KOCKARTS – page 12-20-24)

!        Ions are not here because the sum of all ions abundace is not
!        above 10^-4 until 250 km and we don't have their cpi and aki.

!        Heat capacity for H, He, N, N2, O, O2, CO2, CO: 
!        Circular of the bureau of Standards no. 564: tables of thermal properties of gases comprising [...]
!        Tables of Thermodynamic and Transport Properties of Air, Argon, Carbon Dioxide, Carbon Monoxide, 
!        Hydrogen (molecular and atomic), Nitrogen (molecular and atomic), Oxygen (molecular and atomic), and Steam
!        https://www.govinfo.gov/content/pkg/GOVPUB-C13-89baf9f9b4a43e5f25820bd51b0f3f11/pdf/GOVPUB-C13-89baf9f9b4a43e5f25820bd51b0f3f11.pdf
       

         if (i_co2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = i_co2
            aki(nbq) = 3.072e-4
            cpi(nbq) = 0.834e3
         end if
         if (i_co /= 0) then
            nbq = nbq + 1
            niq(nbq) = i_co
            aki(nbq) = 4.87e-4
            cpi(nbq) = 1.034e3
         end if
         if (i_o /= 0) then
            nbq = nbq + 1
            niq(nbq) = i_o
            aki(nbq) = 7.59e-4
            cpi(nbq) = 1.3e3
         end if
         if (i_o1d /= 0) then
            nbq = nbq + 1
            niq(nbq) = i_o1d
            aki(nbq) = 7.59e-4  !?
            cpi(nbq) = 1.3e3    !?
         end if
         if (i_o2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = i_o2
            aki(nbq) = 5.68e-4
            cpi(nbq) = 0.9194e3
         end if
         if (i_o3 /= 0) then
            nbq = nbq + 1
            niq(nbq) = i_o3
            aki(nbq) = 3.00e-4  !?
            cpi(nbq) = 0.800e3  !?
         end if
         if (i_h /= 0) then
            nbq = nbq + 1
            niq(nbq) = i_h
            !aki(nbq) = 0.0      !! valeur d'origine
            aki(nbq) = 37.9e-4   !??? A verifier
            cpi(nbq) = 20.780e3
         end if
         if (i_d /= 0) then
            nbq = nbq + 1
            niq(nbq) = i_d
            !aki(nbq) = 0.0
            aki(nbq) = 37.9e-4
            cpi(nbq) = 20.780e3
         end if
         if (i_h2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = i_h2
            aki(nbq) = 36.314e-4
            cpi(nbq) = 14.266e3
         end if
         if (i_hd /= 0) then
            nbq = nbq + 1
            niq(nbq) = i_hd
            aki(nbq) = 36.314e-4
            cpi(nbq) = 14.266e3
         end if
         if (i_oh /= 0) then
            nbq = nbq + 1
            niq(nbq) = i_oh
            aki(nbq)  = 7.00e-4 !?
            cpi(nbq)  = 1.045e3
         end if
         if (i_od /= 0) then
            nbq = nbq + 1
            niq(nbq) = i_od
            aki(nbq)  = 7.00e-4
            cpi(nbq)  = 1.045e3
         end if
         if (i_ho2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = i_ho2
            aki(nbq) = 0.0
            cpi(nbq) = 1.065e3  !?
         end if
         if (i_do2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = i_do2
            aki(nbq)  = 7.00e-4
            cpi(nbq)  = 1.045e3
         end if
         if (i_n2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = i_n2
            aki(nbq) = 5.6e-4
            cpi(nbq) = 1.034e3
         end if
c         if (i_ar /= 0) then
c            nbq = nbq + 1
c            niq(nbq) = i_ar
c            aki(nbq) = 0.0      !?
c            cpi(nbq) = 1.000e3  !?
c         end if
         if (i_h2o /= 0) then
            nbq = nbq + 1
            niq(nbq) = i_h2o
            aki(nbq) = 0.0
            cpi(nbq) = 1.870e3
         end if
         if (i_hdo /= 0) then
            nbq = nbq + 1
            niq(nbq) = i_hdo
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         end if
         if (i_n /= 0) then
            nbq = nbq + 1
            niq(nbq) = i_n
            aki(nbq) = 0.0
            cpi(nbq) = 1.4844e3
         endif
         if(i_no /= 0) then
            nbq = nbq + 1
            niq(nbq) = i_no
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(i_no2 /= 0) then
            nbq = nbq + 1
            niq(nbq) = i_no2
            aki(nbq) = 0.0
            cpi(nbq) = 0.0
         endif
         if(i_n2d /= 0) then
            nbq = nbq + 1
            niq(nbq) = i_n2d
            aki(nbq) = 0.0
            cpi(nbq) = 1.4844e3 !?
         endif
         if (i_he /= 0) then
            nbq = nbq + 1
            niq(nbq) = i_he
            aki(nbq) = 29.9e-4
            cpi(nbq) = 5.2e3
         end if

         ! tell the world about it:
         write(*,*) "concentrations: firstcall, nbq=",nbq
!         write(*,*) "test M_tr(nbq)=",M_tr(nbq)
!         write(*,*) "  niq(1:nbq)=",niq(1:nbq)
!         write(*,*) "  aki(1:nbq)=",aki(1:nbq)
!         write(*,*) "  cpi(1:nbq)=",cpi(1:nbq)
         do i = 1,nbq
            write(*,*) "tname(i)=",tname(i)
            write(*,*) "tname(niq(i))=",tname(niq(i))
         end do
         firstcall = .false.
      end if ! if (firstcall)

!     update temperature

      do l = 1,klev
         do ig = 1,klon
            zt(ig,l) = t_seri(ig,l) 
         end do
      end do


!     update mass mixing ratio tracers

      do l = 1,klev
         do ig = 1,klon
            do i = 1,nqmx
!               iq = niq(i) 
               zq(ig,l,i) = max(1.e-30, tr_seri(ig,l,i))
            end do
         end do
      end do

!     mmean : mean molecular mass
!     rho   : mass density [kg/m3]
!     rnew  : specific gas constant
   
      mmean(:,:)  = 0.
      rho(:,:) = 0.     

      do l = 1,klev
         do ig = 1,klon
            do i = 1,nqmx-nmicro
c               iq = niq(i) 
               mmean(ig,l) = mmean(ig,l) + zq(ig,l,i)/M_tr(i)
            end do
            mmean(ig,l) = 1./mmean(ig,l)
            rnew(ig,l) = 8.314/mmean(ig,l)*1.e3     ! J/kg K 
c            write(*,*),'Mmean in concentration2: ',ig, l, mmean(ig,l) 
         end do
      end do

!     cpnew  : specific heat
!     akknew : thermal conductivity cofficient
      
      cpnew(:,:)  = 0.
      akknew(:,:) = 0.

      do l = 1,klev
          do ig = 1,klon

            ntot = pplay(ig,l)/(RKBOL*zt(ig,l))*1.e-6  ! in #/cm3 
            rho(ig,l) = (ntot * mmean(ig,l))/RNAVO*1.e3     ! in kg/m3

c            write(*,*),'Air density: ',ig, l, rho(0,l)            

!!  WARNING -> Cp here below doesn't depend on T (cpdet)

            do i = 1,nbq
               iq = niq(i)
               ni(i) = ntot*zq(ig,l,iq)*mmean(ig,l)/M_tr(iq)
!! On a une super formule pour calculer cp_co2 sur Venus
               if (iq == i_co2) then
                  cpnew(ig,l) = cpnew(ig,l) + ni(i)*cpdet(zt(ig,l))
               else
                  cpnew(ig,l) = cpnew(ig,l) + ni(i)*cpi(i)
               end if
               akknew(ig,l) = akknew(ig,l) + ni(i)*aki(i)
            end do 
 
            cpnew(ig,l) = cpnew(ig,l)/ntot
            akknew(ig,l)= akknew(ig,l)/ntot


          end do 
       end do

      return
      end 
