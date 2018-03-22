#include "cppdefs.h"

# if defined NONLINEAR
      SUBROUTINE set_massflux
      USE mod_ocean
  
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
         Hudy(i,j,k)=Hz_u(i,j,k)*u(i,j,k,nrhs)*dy_u(i,j)
        END DO
       END DO
       DO j=JV_RANGE
        DO i=IV_RANGE
         Hvdx(i,j,k)=Hz_v(i,j,k)*v(i,j,k,nrhs)*dx_v(i,j)
        END DO
       END DO
      END DO

      RETURN
      END

# elif defined TANGENT || defined ADJOINT
      SUBROUTINE fwd_set_massflux(indx)
      USE mod_ocean
  
      implicit none
      integer :: i,j,k,indx
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
         Hudy(i,j,k,indx)=Hz_u(i,j,k)*u(i,j,k,indx)*dy_u(i,j)
        END DO
       END DO
       DO j=JV_RANGE
        DO i=IV_RANGE
         Hvdx(i,j,k,indx)=Hz_v(i,j,k)*v(i,j,k,indx)*dx_v(i,j)
        END DO
       END DO
      END DO

      RETURN
      END

!************************      
      SUBROUTINE fwd_set_massflux_trc
      USE mod_ocean
  
      implicit none
      integer :: i,j,k

      DO k=1,N
       DO j=LBj,UBj
        DO i=LBi,UBi
         Hudy(i,j,k,3)=0.5d0*(Hudy(i,j,k,nstp)+Hudy(i,j,k,nnew))
         Hvdx(i,j,k,3)=0.5d0*(Hvdx(i,j,k,nstp)+Hvdx(i,j,k,nnew))
        END DO
       END DO
      END DO

      RETURN
      END
      
# endif
      
