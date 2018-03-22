#include "cppdefs.h"

#ifdef TANGENT
      SUBROUTINE tl_set_depth(LBi,UBi,LBj,UBj,N,Lm,Mm,
     *                     z0_r,z0_w,oh,
     *                     tl_zeta,tl_z_r,tl_z_w,tl_Hz,tl_Hz_u,tl_Hz_v)

#ifdef MPI
      USE mod_parallel
#endif
      implicit none
      integer :: i,j,k,LBi,UBi,LBj,UBj,N,Lm,Mm
      real(8) :: z0_w(LBi:UBi,LBj:UBj,0:N)
      real(8) :: z0_r(LBi:UBi,LBj:UBj,1:N)
      real(8) :: oh(LBi:UBi,LBj:UBj)
      real(8) :: tl_z_r(LBi:UBi,LBj:UBj,N)
      real(8) :: tl_z_w(LBi:UBi,LBj:UBj,0:N)
      real(8) :: tl_zeta(LBi:UBi,LBj:UBj)
      real(8) :: tl_Hz(LBi:UBi,LBj:UBj,N)
      real(8) :: tl_Hz_u(LBi:UBi,LBj:UBj,N)
      real(8) :: tl_Hz_v(LBi:UBi,LBj:UBj,N)

! TL loop 1:
      DO k=0,N
       DO j=LBj,UBj
        DO i=LBi,UBi
         tl_z_w(i,j,k)=(1.d0+z0_w(i,j,k)*oh(i,j))*tl_zeta(i,j)
        END DO
       END DO
      END DO

! TL loop 2:
      DO k=1,N
       DO j=LBj,UBj
        DO i=LBi,UBi
         tl_z_r(i,j,k)=(1.d0+z0_r(i,j,k)*oh(i,j))*tl_zeta(i,j)
         tl_Hz(i,j,k)=tl_z_w(i,j,k)-tl_z_w(i,j,k-1)
        END DO
       END DO
      END DO

! TL loop 3:
      DO k=1,N
       DO j=LBj,UBj
        DO i=LBi+1,UBi
         tl_Hz_u(i,j,k)=0.5d0*(tl_Hz(i-1,j,k)+tl_Hz(i,j,k))
        END DO
       END DO
      END DO

! TL loop 4:
      DO k=1,N
       DO j=LBj+1,UBj
        DO i=LBi,UBi
         tl_Hz_v(i,j,k)=0.5d0*(tl_Hz(i,j-1,k)+tl_Hz(i,j,k))
        END DO
       END DO
      END DO

#if defined NS_PERIODIC

# ifdef MPI
      STOP('periodic BC not coded in MPI')
# else /* serial */
      CALL exchange_3d(tl_Hz_v,LBi,UBi,LBj,UBj,1,N,Lm,Mm)
# endif /* MPI or serial */

#endif /* NS_PERIODIC */

      RETURN
      END
#endif /* #ifdef TANGENT */
