#include "cppdefs.h"

#ifdef TANGENT
      SUBROUTINE tl_set_massflux
      USE mod_ocean
#ifdef MPI
      USE mod_parallel
#endif  

      implicit none
      integer :: i,j,k
!
! see lims_moment_flux.txt
!

# define IU_RANGE min(IstrU-1,IstrV),max(IendU+1,IendV+1)
# define JU_RANGE min(JstrU,JstrV-1),max(JendU,JendV)

# define IV_RANGE min(IstrU-1,IstrV),max(IendU,IendV)
# define JV_RANGE min(JstrU,JstrV-1),max(JendU+1,JendV+1)

      DO k=1,N
       DO j=JU_RANGE
        DO i=IU_RANGE
! NL:         Hudy(i,j,k)=Hz_u(i,j,k)*u(i,j,k,nrhs)*dy_u(i,j)
         tl_Hudy(i,j,k)=(tl_Hz_u(i,j,k)*u(i,j,k,nrhs)+
     *                   Hz_u(i,j,k)*tl_u(i,j,k,nrhs))*dy_u(i,j)
        END DO
       END DO

       DO j=JV_RANGE
        DO i=IV_RANGE
! NL:         Hvdx(i,j,k)=Hz_v(i,j,k)*v(i,j,k,nrhs)*dx_v(i,j)
         tl_Hvdx(i,j,k)=(tl_Hz_v(i,j,k)*v(i,j,k,nrhs)+
     *                   Hz_v(i,j,k)*tl_v(i,j,k,nrhs))*dx_v(i,j)
        END DO
       END DO
      END DO

      RETURN
      END

#endif /* TANGENT */
