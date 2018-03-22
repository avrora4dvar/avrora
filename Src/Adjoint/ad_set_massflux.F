#include "cppdefs.h"

#ifdef ADJOINT
      SUBROUTINE ad_set_massflux
      USE mod_ocean

      implicit none
      integer :: i,j,k
      real(8) :: tmp

# define IU_RANGE min(IstrU-1,IstrV),max(IendU+1,IendV+1)
# define JU_RANGE min(JstrU,JstrV-1),max(JendU,JendV)

# define IV_RANGE min(IstrU-1,IstrV),max(IendU,IendV)
# define JV_RANGE min(JstrU,JstrV-1),max(JendU+1,JendV+1)

      DO k=1,N

       DO j=JV_RANGE
        DO i=IV_RANGE
! NL:         Hvdx(i,j,k)=Hz_v(i,j,k)*v(i,j,k,nrhs)*dx_v(i,j)

! TL:         tl_Hvdx(i,j,k)=(tl_Hz_v(i,j,k)*v(i,j,k,nrhs)+
! TL:     *                   Hz_v(i,j,k)*tl_v(i,j,k,nrhs))*dx_v(i,j)

          tmp=ad_Hvdx(i,j,k)*dx_v(i,j)

          ad_Hz_v(i,j,k)=ad_Hz_v(i,j,k)+tmp*v(i,j,k,nrhs)
          ad_v(i,j,k,nrhs)=ad_v(i,j,k,nrhs)+tmp*Hz_v(i,j,k)
        
          ad_Hvdx(i,j,k)=0.d0

        END DO
       END DO

       DO j=JU_RANGE
        DO i=IU_RANGE
! NL:         Hudy(i,j,k)=Hz_u(i,j,k)*u(i,j,k,nrhs)*dy_u(i,j)
! TL:         tl_Hudy(i,j,k)=(tl_Hz_u(i,j,k)*u(i,j,k,nrhs)+
! TL:     *                   Hz_u(i,j,k)*tl_u(i,j,k,nrhs))*dy_u(i,j)

         tmp=ad_Hudy(i,j,k)*dy_u(i,j)

         ad_Hz_u(i,j,k)=ad_Hz_u(i,j,k)+tmp*u(i,j,k,nrhs)
         ad_u(i,j,k,nrhs)=ad_u(i,j,k,nrhs)+tmp*Hz_u(i,j,k)

         ad_Hudy(i,j,k)=0.d0

        END DO
       END DO

      END DO

      RETURN
      END
#endif /* ADJOINT */
