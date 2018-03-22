#include "cppdefs.h"

      SUBROUTINE ad_t3dmix2
# if defined ADJOINT && defined TS_DIF2
      USE mod_ocean

      implicit none

      integer :: i,j,k,itrc
      real(8) :: cff,tmp,tmp1,tmp2,ad_cff

      real(8) :: ad_FX(LBi:UBi,LBj:UBj)
      real(8) :: ad_FE(LBi:UBi,LBj:UBj)

      ad_FX=0.d0
      ad_FE=0.d0

      ad_cff=0.d0

      DO itrc=1,Ntr
       DO k=1,N
!
! ADJ to: Time-step harmonic, S-surfaces diffusion term (m Tunits).
!
! ADJ to TL CYCLE 3:
        DO j=Jstr,Jend
         DO i=Istr,Iend
! tl line 2:          tl_t(i,j,k,nnew,itrc)=tl_t(i,j,k,nnew,itrc)+tl_cff
          ad_cff=ad_cff+ad_t(i,j,k,nnew,itrc)

! tl line 1:          tl_cff=dt*odx(i,j)*ody(i,j)*
!     *           (tl_FX(i+1,j)-tl_FX(i,j)+tl_FE(i,j+1)-tl_FE(i,j))
          tmp=dt*odx(i,j)*ody(i,j)*ad_cff
          ad_FX(i+1,j)=ad_FX(i+1,j)+tmp
          ad_FX(i,j)  =ad_FX(i,j)  -tmp
          ad_FE(i,j+1)=ad_FE(i,j+1)+tmp
          ad_FE(i,j)  =ad_FE(i,j)  -tmp
          ad_cff=0.d0
         END DO
        END DO

!
!  ADJ to Compute XI- and ETA-components of diffusive tracer flux (T m3/s).
!
! ADJ to TL CYCLE 2:
        DO j=Jstr,Jend+1
         DO i=Istr,Iend
          cff=0.25d0*(diff2(i,j,itrc)+diff2(i,j-1,itrc))*dxody_v(i,j)

! tl line 2:          tl_FE(i,j)=tl_FE(i,j)*DBLE(mask_v(i,j))
          ad_FE(i,j)=ad_FE(i,j)*DBLE(mask_v(i,j))

! tl line 1:          tl_FE(i,j)=cff*(
!     *                (tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*
!     *                (t(i,j,k,nrhs,itrc)-t(i,j-1,k,nrhs,itrc))
!     *               +(Hz(i,j,k)+Hz(i,j-1,k))*
!     *                (tl_t(i,j,k,nrhs,itrc)-tl_t(i,j-1,k,nrhs,itrc))
!     *                   )

          tmp1=cff*(t(i,j,k,nrhs,itrc)-t(i,j-1,k,nrhs,itrc))*ad_FE(i,j)
          tmp2=cff*(Hz(i,j,k)+Hz(i,j-1,k))*ad_FE(i,j)

          ad_Hz(i,j,k)  =ad_Hz(i,j,k)  +tmp1
          ad_Hz(i,j-1,k)=ad_Hz(i,j-1,k)+tmp1
          ad_t(i,j,k,nrhs,itrc)  =ad_t(i,j,k,nrhs,itrc)  +tmp2
          ad_t(i,j-1,k,nrhs,itrc)=ad_t(i,j-1,k,nrhs,itrc)-tmp2

          ad_FE(i,j)=0.d0

         END DO
        END DO

! ADJ to TL CYCLE 1:
        DO j=Jstr,Jend
         DO i=Istr,Iend+1
          cff=0.25d0*(diff2(i,j,itrc)+diff2(i-1,j,itrc))*dyodx_u(i,j)

! tl line 2:          tl_FX(i,j)=tl_FX(i,j)*DBLE(mask_u(i,j))
          ad_FX(i,j)=ad_FX(i,j)*DBLE(mask_u(i,j))


! tl line 1:          tl_FX(i,j)=cff*(
!     *                (tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*
!     *                (t(i,j,k,nrhs,itrc)-t(i-1,j,k,nrhs,itrc))
!     *               +(Hz(i,j,k)+Hz(i-1,j,k))*
!     *                (tl_t(i,j,k,nrhs,itrc)-tl_t(i-1,j,k,nrhs,itrc))
!     *                   )
          tmp1=cff*(t(i,j,k,nrhs,itrc)-t(i-1,j,k,nrhs,itrc))*ad_FX(i,j)
          tmp2=cff*(Hz(i,j,k)+Hz(i-1,j,k))*ad_FX(i,j)

          ad_Hz(i,j,k)  =ad_Hz(i,j,k)  +tmp1
          ad_Hz(i-1,j,k)=ad_Hz(i-1,j,k)+tmp1
          ad_t(i,j,k,nrhs,itrc)  =ad_t(i,j,k,nrhs,itrc)  +tmp2
          ad_t(i-1,j,k,nrhs,itrc)=ad_t(i-1,j,k,nrhs,itrc)-tmp2
          ad_FX(i,j)=0.d0

         END DO
        END DO


       END DO ! k
      END DO ! itrc
 
# endif /* if defined ADJOINT && defined TS_DIF2 */

      RETURN
      END
