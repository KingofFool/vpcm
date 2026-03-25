c**********************************************************************

      subroutine hrtherm(ig,euvmod,rm,nespeuv,tx,iz,zenit,jtot)


c     feb 2002        fgg           first version
c     nov 2002        fgg           second version

c**********************************************************************
      use dimphy
      use param_v4_h, only: ninter,nabs,jfotsout,fluxtop,freccen
      use clesphys_mod

      implicit none

c     common variables and constants


c#include "clesphys.h"


c    local parameters and variables

      real       xabsi(nabs,klev) 			!densities (cm^-3)
      real       jergs(ninter,nabs,klev)
      
      integer    i,j,k,indexint          !indexes
      character  dn


c     input and output variables

      integer    ig  ,euvmod 
      integer    nespeuv 
      real       rm(klev,nespeuv)              !density matrix (cm^-3)
      real       jtot(klev)                    !output: heating rate(erg/s cm3)
      real       tx(klev)                      !temperature
      real       zenit
      real       iz(klev)

      ! tracer indexes for the EUV heating:
!!! ATTENTION. These values have to be identical to those in euvheat.F90
!!! If the values are changed there, the same has to be done here  !!!

!      integer,parameter :: i_co2=1
!      integer,parameter :: i_n2=13
!      integer,parameter :: i_n=14
!      integer,parameter :: i_o=3
!      integer,parameter :: i_co=4

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

c*************************PROGRAM STARTS*******************************

      !If nighttime, photoabsorption coefficient set to 0
      if(zenit.gt.90.) then  !140 in the martian routine
         dn='n'
         else
         dn='d'
      end if
      if(dn.eq.'n') then
        do i=1,klev                                    
	      jtot(i)=0.
        enddo       
        return
      endif 

      !initializations
      jergs(:,:,:)=0.
      xabsi(:,:)=0.
      jtot(:)=0.
      !All number densities to a single array, xabsi(species,layer)
      ! WARNING xabs(nabs,nlev), j=1,nabs --> the values of j should 
      !         be the same for xabs than for jfotsout(indexint,j,i) 
      !
      do i=1,klev
         xabsi(1,i)  = rm(i,ix_co2)    ! CO2
         xabsi(2,i)  = rm(i,ix_o2)     ! O2
         xabsi(3,i)  = rm(i,ix_o)      ! O(3P)
         xabsi(4,i)  = rm(i,ix_h2o)    ! H2O
         xabsi(5,i)  = rm(i,ix_h2)     ! H2
         xabsi(6,i)  = rm(i,ix_h2o2)   ! H2O2
         !Only if O3, N or ion chemistry requested
         if(euvmod.ge.1) then
            xabsi(7,i)  = rm(i,ix_o3)  ! O3
         endif
         xabsi(8,i)  = rm(i,ix_n2)     ! N2
         !Only if N or ion chemistry requested
         if(euvmod.ge.2) then
            xabsi(9,i)  = rm(i,ix_n)   ! N
            xabsi(10,i) = rm(i,ix_no)  ! NO
            xabsi(13,i) = rm(i,ix_no2) ! NO2
         endif
         xabsi(11,i) = rm(i,ix_co)     ! CO
         xabsi(12,i) = rm(i,ix_h)      ! H
      end do

      !Calculation of photoabsortion coefficient
      call jthermcalc_e107(ig,klev,euvmod,rm,nespeuv,tx,iz,zenit)

      !Total photoabsorption coefficient    !  erg/(s*cm3) 
      do i=1,klev
         jtot(i)=0.
        do j=1,nabs
          do indexint=1,ninter
            jergs(indexint,j,i) = jfotsout(indexint,j,i) 
     $              * xabsi (j,i) * fluxtop(indexint)  
     $              / (0.5e9 * freccen(indexint))
            jtot(i)=jtot(i)+jergs(indexint,j,i)    
 

          end do
        end do
      end do

      return

      end

