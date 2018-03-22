#include "cppdefs.h"

#ifdef ADJOINT

      SUBROUTINE ad_set_data
      USE mod_ocean
      USE mod_forces

      implicit none
      integer :: i,j
      real(8) :: cff1,cff2

      cff2=(time-sms_timeG(ifrc1))/(sms_timeG(ifrc2)-sms_timeG(ifrc1))
      cff1=1.d0-cff2

! TL:      tl_stflx=0.d0
      ad_stflx=0.d0

      DO j=LBj,UBj
      DO i=LBi,UBi
! TL:       tl_sustr(i,j)=cff1*tl_sustrG(i,j,ifrc1)+cff2*tl_sustrG(i,j,ifrc2)
       ad_sustrG(i,j,ifrc1)=ad_sustrG(i,j,ifrc1)+cff1*ad_sustr(i,j)
       ad_sustrG(i,j,ifrc2)=ad_sustrG(i,j,ifrc2)+cff2*ad_sustr(i,j)
       ad_sustr(i,j)=0.d0
! TL:       tl_svstr(i,j)=cff1*tl_svstrG(i,j,1)+cff2*tl_svstrG(i,j,2)
       ad_svstrG(i,j,ifrc1)=ad_svstrG(i,j,ifrc1)+cff1*ad_svstr(i,j)
       ad_svstrG(i,j,ifrc2)=ad_svstrG(i,j,ifrc2)+cff2*ad_svstr(i,j)
       ad_svstr(i,j)=0.d0
!
      END DO
      END DO

      RETURN
      END

#endif


