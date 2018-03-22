#include "cppdefs.h"

      SUBROUTINE set_depth(LBi,UBi,LBj,UBj,N,Lm,Mm,
     *               z0_r,z0_w,oh,
     *               zeta,z_r,z_w,Hz,Hz_u,Hz_v)
#ifdef MPI
       USE mod_parallel     
#endif
      
      implicit none
      integer :: i,j,k,LBi,UBi,LBj,UBj,N,Lm,Mm
      real(8) :: z0_w(LBi:UBi,LBj:UBj,0:N)
      real(8) :: z0_r(LBi:UBi,LBj:UBj,1:N)
      real(8) :: z_r(LBi:UBi,LBj:UBj,N)
      real(8) :: z_w(LBi:UBi,LBj:UBj,0:N)
      real(8) :: zeta(LBi:UBi,LBj:UBj)
      real(8) :: oh(LBi:UBi,LBj:UBj)
      real(8) :: Hz(LBi:UBi,LBj:UBj,N)
      real(8) :: Hz_u(LBi:UBi,LBj:UBj,N)
      real(8) :: Hz_v(LBi:UBi,LBj:UBj,N)

      DO k=0,N
       DO j=LBj,UBj
        DO i=LBi,UBi
         z_w(i,j,k)=z0_w(i,j,k)+(1.d0+z0_w(i,j,k)*oh(i,j))*zeta(i,j)
        END DO
       END DO
      END DO

      DO k=1,N
       DO j=LBj,UBj
        DO i=LBi,UBi
         z_r(i,j,k)=z0_r(i,j,k)+(1.d0+z0_r(i,j,k)*oh(i,j))*zeta(i,j)
         Hz(i,j,k)=z_w(i,j,k)-z_w(i,j,k-1)
        END DO
       END DO
      END DO

      DO k=1,N
       DO j=LBj,UBj
        DO i=LBi+1,UBi
         Hz_u(i,j,k)=0.5d0*(Hz(i-1,j,k)+Hz(i,j,k))
        END DO
       END DO
      END DO

      DO k=1,N
       DO j=LBj+1,UBj
        DO i=LBi,UBi
         Hz_v(i,j,k)=0.5d0*(Hz(i,j-1,k)+Hz(i,j,k))
        END DO
       END DO
      END DO

#if defined NS_PERIODIC
# ifdef MPI
      STOP('periodic BC not coded in MPI')
# else /* serial */
      CALL exchange_3d(Hz_v,LBi,UBi,LBj,UBj,1,N,Lm,Mm)
# endif /* MPI or serial */
#endif /* NS_PERIODIC */

      RETURN
      END
