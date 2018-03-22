#include "cppdefs.h"

#ifdef TANGENT
      SUBROUTINE tl_init_zeta
      USE mod_ocean

      implicit none
      integer :: i,j

      IF (iic.gt.1) THEN

       DO j=LBj,UBj
        DO i=LBi,UBi
         tl_zeta(i,j,1)=tl_Zt_avg1(i,j)
        END DO
       END DO

      END IF 

      RETURN
      END
#endif
