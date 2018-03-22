#include "cppdefs.h"

#ifdef TANGENT
      SUBROUTINE tl_set_vbc
      USE mod_ocean

      implicit none       
      integer :: i,j

! Linear bottom stress:
      DO j=JstrU,JendU
      DO i=IstrU,IendU
       tl_bustr(i,j)=rdrg*tl_u(i,j,1,nrhs)
      END DO
      END DO

      DO j=JstrV,JendV
      DO i=IstrV,IendV
       tl_bvstr(i,j)=rdrg*tl_v(i,j,1,nrhs)
      END DO
      END DO

      RETURN
      END

#endif /* TANGENT */
