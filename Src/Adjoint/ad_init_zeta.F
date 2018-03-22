#include "cppdefs.h"

#ifdef ADJOINT
      SUBROUTINE ad_init_zeta
      USE mod_ocean

      implicit none
      integer :: i,j

      IF (iic.gt.1) THEN
       DO j=LBj,UBj
        DO i=LBi,UBi
! TL:         tl_zeta(i,j,1)=tl_Zt_avg1(i,j)
         ad_Zt_avg1(i,j)=ad_Zt_avg1(i,j)+ad_zeta(i,j,1)
         ad_zeta(i,j,1)=0.d0
        END DO
       END DO
      END IF

      RETURN
      END
#endif
