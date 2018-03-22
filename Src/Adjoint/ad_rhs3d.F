#include "cppdefs.h"

#ifdef ADJOINT 
      SUBROUTINE ad_rhs3d
      USE mod_ocean
      USE mod_forces
      implicit none
      
      integer i,j,k
      real(8) :: cff
      real(8) :: q0

!      write (*,*) 'in ad_rhs3d:...'

!-----------------------------------------------------------------------
!  Adj to : Compute horizontal, harmonic mixing of momentum.
!-----------------------------------------------------------------------
# ifdef UV_VIS2
      CALL ad_uv3dmix2
# endif /* UV_VIS2 */

!-----------------------------------------------------------------------
!  Compute ADJ to TL forcing term for the 2D momentum equations.
!-----------------------------------------------------------------------
!

! ADJ to TL loop 4:
      DO k=1,N
       DO j=JstrV,JendV
        DO i=IstrV,IendV
! TL (recurs.): tl_rvfrc(i,j)=tl_rvfrc(i,j)+tl_rv(i,j,k,nrhs)
         ad_rv(i,j,k,nrhs)=ad_rv(i,j,k,nrhs)+ad_rvfrc(i,j)
        END DO
       END DO
      END DO   

! ADJ to TL loop 3:
      DO j=JstrV,JendV
       DO i=IstrV,IendV
! TL:        tl_rvfrc(i,j)=(tl_svstr(i,j)-tl_bvstr(i,j))*Area_v(i,j)
        cff=ad_rvfrc(i,j)*Area_v(i,j)
        ad_svstr(i,j)=ad_svstr(i,j)+cff
        ad_bvstr(i,j)=ad_bvstr(i,j)-cff
        ad_rvfrc(i,j)=0.d0
       END DO
      END DO

! ADJ to TL loop 2:
      DO k=1,N
       DO j=JstrU,JendU
        DO i=IstrU,IendU
! TL (recurs.): tl_rufrc(i,j)=tl_rufrc(i,j)+tl_ru(i,j,k,nrhs)
         ad_ru(i,j,k,nrhs)=ad_ru(i,j,k,nrhs)+ad_rufrc(i,j)
        END DO
       END DO
      END DO   

! ADJ to TL loop 1:
      DO j=JstrU,JendU
       DO i=IstrU,IendU
! TL:        tl_rufrc(i,j)=(tl_sustr(i,j)-tl_bustr(i,j))*Area_u(i,j)
        cff=ad_rufrc(i,j)*Area_u(i,j)
        ad_sustr(i,j)=ad_sustr(i,j)+cff
        ad_bustr(i,j)=ad_bustr(i,j)-cff
        ad_rufrc(i,j)=0.d0
       END DO
      END DO

# ifdef UV_ADV
      CALL ad_uv_adv
# endif

# if defined UV_COR || defined UV_COR_MASK
      CALL ad_coriolis
# endif

!-----------------------------------------------------------------------
! T, S: horizontal diffusion
!-----------------------------------------------------------------------
# ifdef TS_DIF2
      CALL ad_t3dmix2
# endif

      CALL ad_prsgrd 

      CALL ad_pre_step3d 
        
      RETURN
      END SUBROUTINE ad_rhs3d
#endif /* ADJOINT */
