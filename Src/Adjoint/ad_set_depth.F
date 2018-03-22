#include "cppdefs.h"

#ifdef ADJOINT
      SUBROUTINE ad_set_depth(LBi,UBi,LBj,UBj,N,Lm,Mm,
     *          z0_r,z0_w,oh,
     *          ad_zeta,ad_z_r,ad_z_w,ad_Hz,ad_Hz_u,ad_Hz_v)

! Note: ad_z_r,ad_z_w,ad_Hz,ad_Hz_u,ad_Hz_v will be used as inputs and
! zero out after used
!       ad_zeta(:,:) will accumulate sensitivities to above mentioned 
! fields and be output

      implicit none

      integer :: i,j,k,LBi,UBi,LBj,UBj,N,Lm,Mm
      real(8) :: z0_w(LBi:UBi,LBj:UBj,0:N)
      real(8) :: z0_r(LBi:UBi,LBj:UBj,1:N)
      real(8) :: oh(LBi:UBi,LBj:UBj)
      real(8) :: ad_z_r(LBi:UBi,LBj:UBj,N)
      real(8) :: ad_z_w(LBi:UBi,LBj:UBj,0:N)
      real(8) :: ad_zeta(LBi:UBi,LBj:UBj)
      real(8) :: ad_Hz(LBi:UBi,LBj:UBj,N)
      real(8) :: ad_Hz_u(LBi:UBi,LBj:UBj,N)
      real(8) :: ad_Hz_v(LBi:UBi,LBj:UBj,N)

      real(8) :: cff

#if defined NS_PERIODIC
      CALL ad_exchange_3d(ad_Hz_v,LBi,UBi,LBj,UBj,1,N,Lm,Mm)
#endif 

! ADJ to TL loop 4:
      DO k=1,N
       DO j=LBj+1,UBj
        DO i=LBi,UBi
! TL: tl_Hz_v(i,j,k)=0.5d0*(tl_Hz(i,j-1,k)+tl_Hz(i,j,k))
         cff=0.5d0*ad_Hz_v(i,j,k)
         ad_Hz(i,j-1,k)=ad_Hz(i,j-1,k)+cff
         ad_Hz(i,j,k)=ad_Hz(i,j,k)+cff
         ad_Hz_v(i,j,k)=0.d0
        END DO
       END DO
      END DO


! ADJ to TL loop 3:
      DO k=1,N
       DO j=LBj,UBj
        DO i=LBi+1,UBi
! TL: tl_Hz_u(i,j,k)=0.5d0*(tl_Hz(i-1,j,k)+tl_Hz(i,j,k))
         cff=0.5d0*ad_Hz_u(i,j,k)
         ad_Hz(i-1,j,k)=ad_Hz(i-1,j,k)+cff
         ad_Hz(i,j,k)=ad_Hz(i,j,k)+cff
         ad_Hz_u(i,j,k)=0.d0
        END DO
       END DO
      END DO

! ADJ to TL loop 2:
      DO k=1,N
       DO j=LBj,UBj
        DO i=LBi,UBi

! TL line 2: tl_Hz(i,j,k)=tl_z_w(i,j,k)-tl_z_w(i,j,k-1)
         ad_z_w(i,j,k)=ad_z_w(i,j,k)+ad_Hz(i,j,k)
         ad_z_w(i,j,k-1)=ad_z_w(i,j,k-1)-ad_Hz(i,j,k)
         ad_Hz(i,j,k)=0.d0
       
! TL line 1: tl_z_r(i,j,k)=(1.d0+z0_r(i,j,k)*oh(i,j))*tl_zeta(i,j)
         ad_zeta(i,j)=ad_zeta(i,j)+
     *                 (1.d0+z0_r(i,j,k)*oh(i,j))*ad_z_r(i,j,k)
         ad_z_r(i,j,k)=0.d0
        END DO
       END DO
      END DO

! ADJ to TL loop 1:
      DO k=0,N
       DO j=LBj,UBj
        DO i=LBi,UBi
! TL:         tl_z_w(i,j,k)=(1.d0+z0_w(i,j,k)*oh(i,j))*tl_zeta(i,j)
         ad_zeta(i,j)=ad_zeta(i,j)+
     *                  (1.d0+z0_w(i,j,k)*oh(i,j))*ad_z_w(i,j,k)
         ad_z_w(i,j,k)=0.d0
        END DO
       END DO
      END DO

      RETURN
      END

#endif /* ADJOINT */
