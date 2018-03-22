#include "cppdefs.h"

#ifdef TANGENT 
      SUBROUTINE tl_rhs3d
      USE mod_ocean
      USE mod_forces
#ifdef MPI
      USE mod_parallel
#endif
  
      implicit none
      integer :: i,j,k

      character(len=80) :: fhead

!-----------------------------------------------------------------------
!  Initialize computations for the new time step of the 3D primitive
!  variables.
!-----------------------------------------------------------------------
      CALL tl_pre_step3d 

!-----------------------------------------------------------------------
!  Compute baroclinic pressure gradient.
      CALL tl_prsgrd 
    
!-----------------------------------------------------------------------
! T, S: horizontal diffusion
!-----------------------------------------------------------------------
# ifdef TS_DIF2
      CALL tl_t3dmix2
# endif


!-----------------------------------------------------------------------
!  Compute right-hand-side terms for the 3D momentum equations.
!-----------------------------------------------------------------------
# if defined UV_COR || defined UV_COR_MASK
       CALL tl_coriolis
# endif

# ifdef UV_ADV
      CALL tl_uv_adv
# endif

!-----------------------------------------------------------------------
!  Compute TL forcing term for the 2D momentum equations.
!-----------------------------------------------------------------------
!
!  Vertically integrate baroclinic right-hand-side terms. 
!  Add in the difference between surface and bottom stresses.
!

! TL loop 1:
      DO j=JstrU,JendU
       DO i=IstrU,IendU
        tl_rufrc(i,j)=(tl_sustr(i,j)-tl_bustr(i,j))*Area_u(i,j)
       END DO
      END DO

! TL loop 2:
      DO k=1,N
       DO j=JstrU,JendU
        DO i=IstrU,IendU
         tl_rufrc(i,j)=tl_rufrc(i,j)+tl_ru(i,j,k,nrhs)
        END DO
       END DO
      END DO   

! TL loop 3:
      DO j=JstrV,JendV
       DO i=IstrV,IendV
        tl_rvfrc(i,j)=(tl_svstr(i,j)-tl_bvstr(i,j))*Area_v(i,j)
       END DO
      END DO

! TL loop 4:
      DO k=1,N
       DO j=JstrV,JendV
        DO i=IstrV,IendV
         tl_rvfrc(i,j)=tl_rvfrc(i,j)+tl_rv(i,j,k,nrhs)
        END DO
       END DO
      END DO   

!-----------------------------------------------------------------------
!  Compute horizontal, harmonic mixing of momentum.
!-----------------------------------------------------------------------
# ifdef UV_VIS2
      CALL tl_uv3dmix2
# endif /* UV_VIS2 */

      RETURN
      END SUBROUTINE tl_rhs3d
#endif /* TANGENT */
