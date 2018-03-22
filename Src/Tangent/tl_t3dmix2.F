#include "cppdefs.h"

      SUBROUTINE tl_t3dmix2
# if defined TANGENT && defined TS_DIF2
      USE mod_ocean

      implicit none

      integer :: i,j,k,itrc
      real(8) :: cff,tl_cff

      real(8) :: tl_FX(LBi:UBi,LBj:UBj)
      real(8) :: tl_FE(LBi:UBi,LBj:UBj)

      tl_FX=0.d0
      tl_FE=0.d0

      DO itrc=1,Ntr
       DO k=1,N
!
!  Compute XI- and ETA-components of diffusive tracer flux (T m3/s).
!
! TL CYCLE 1:
        DO j=Jstr,Jend
         DO i=Istr,Iend+1
          cff=0.25d0*(diff2(i,j,itrc)+diff2(i-1,j,itrc))*dyodx_u(i,j)
! NL:          FX(i,j)=cff*(Hz(i,j,k)+Hz(i-1,j,k))*
!     *                (t(i,j,k,nrhs,itrc)-t(i-1,j,k,nrhs,itrc))
          tl_FX(i,j)=cff*(
     *                (tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*
     *                (t(i,j,k,nrhs,itrc)-t(i-1,j,k,nrhs,itrc))
     *               +(Hz(i,j,k)+Hz(i-1,j,k))*
     *                (tl_t(i,j,k,nrhs,itrc)-tl_t(i-1,j,k,nrhs,itrc))
     *                   )
! NL:          FX(i,j)=FX(i,j)*DBLE(mask_u(i,j))
          tl_FX(i,j)=tl_FX(i,j)*DBLE(mask_u(i,j))
         END DO
        END DO

! TL CYCLE 2:
        DO j=Jstr,Jend+1
         DO i=Istr,Iend
          cff=0.25d0*(diff2(i,j,itrc)+diff2(i,j-1,itrc))*dxody_v(i,j)
! NL:          FE(i,j)=cff*(Hz(i,j,k)+Hz(i,j-1,k))*
!     *                (t(i,j,k,nrhs,itrc)-t(i,j-1,k,nrhs,itrc))
          tl_FE(i,j)=cff*(
     *                (tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*
     *                (t(i,j,k,nrhs,itrc)-t(i,j-1,k,nrhs,itrc))
     *               +(Hz(i,j,k)+Hz(i,j-1,k))*
     *                (tl_t(i,j,k,nrhs,itrc)-tl_t(i,j-1,k,nrhs,itrc))
     *                   )


! NL:          FE(i,j)=FE(i,j)*DBLE(mask_v(i,j))
          tl_FE(i,j)=tl_FE(i,j)*DBLE(mask_v(i,j))
         END DO
        END DO
!
! Time-step harmonic, S-surfaces diffusion term (m Tunits).
!
! TL CYCLE 3:
        DO j=Jstr,Jend
         DO i=Istr,Iend
          tl_cff=dt*odx(i,j)*ody(i,j)*
     *           (tl_FX(i+1,j)-tl_FX(i,j)+tl_FE(i,j+1)-tl_FE(i,j))
          tl_t(i,j,k,nnew,itrc)=tl_t(i,j,k,nnew,itrc)+tl_cff
         END DO
        END DO

       END DO ! k
      END DO ! itrc
 
# endif /* if defined TANGENT && defined TS_DIF2 */

      RETURN
      END
