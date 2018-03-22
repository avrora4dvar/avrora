#include "cppdefs.h"

#ifdef TANGENT

      SUBROUTINE tl_set_data    
      USE mod_ocean
      USE mod_forces
 
      implicit none
      integer :: i,j
      real(8) :: cff1,cff2

      cff2=(time-sms_timeG(ifrc1))/(sms_timeG(ifrc2)-sms_timeG(ifrc1))
      cff1=1.d0-cff2

      DO j=LBj,UBj
      DO i=LBi,UBi
       tl_sustr(i,j)=cff1*tl_sustrG(i,j,ifrc1)+cff2*tl_sustrG(i,j,ifrc2)  
       tl_svstr(i,j)=cff1*tl_svstrG(i,j,ifrc1)+cff2*tl_svstrG(i,j,ifrc2)  
      END DO
      END DO

      tl_stflx=0.d0

      RETURN
      END

#endif
