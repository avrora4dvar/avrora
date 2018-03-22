#include "cppdefs.h"

      SUBROUTINE set_visc2
      USE mod_ocean

      implicit none

      integer i,j,itr

# ifdef UV_VIS2
      DO j=LBj,UBj
      DO i=LBi,UBi
       visc2_r(i,j)=vis2
      END DO
      END DO
   
      DO j=LBj+1,UBj
      DO i=LBi+1,UBi
       visc2_p(i,j)=0.25d0*(visc2_r(i-1,j-1)+visc2_r(i-1,j)+
     *                      visc2_r(i,j-1)+visc2_r(i,j))
      END DO
      END DO
            
# endif /* UV_VIS2 */

# ifdef TS_DIF2

     
      DO itr=1,Ntr
      DO j=LBj,UBj
      DO i=LBi,UBi
       diff2(i,j,itr)=dif2
      END DO
      END DO
      END DO
            
# endif /* TS_DIF2 */

      RETURN
      END
   
