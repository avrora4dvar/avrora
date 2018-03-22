#include "cppdefs.h"
      MODULE mod_forces
      implicit none

      integer :: nsms,sms_rec,ifrc1,ifrc2

! Forcing:
      real(8), allocatable :: sms_time(:)
      real(8) :: sms_timeG(2)

#if defined NONLINEAR
      real(8), allocatable :: sustr(:,:)
      real(8), allocatable :: svstr(:,:)
      real(8), allocatable :: sustrG(:,:,:)
      real(8), allocatable :: svstrG(:,:,:)

      real(8), allocatable :: stflx(:,:,:)

!
#elif defined TANGENT
      real(8), allocatable :: tl_sustr(:,:)
      real(8), allocatable :: tl_svstr(:,:)
      real(8), allocatable :: tl_sustrG(:,:,:)
      real(8), allocatable :: tl_svstrG(:,:,:)

      real(8), allocatable :: tl_stflx(:,:,:)
!
#elif defined ADJOINT
      real(8), allocatable :: ad_sustr(:,:)
      real(8), allocatable :: ad_svstr(:,:)
      real(8), allocatable :: ad_sustrG(:,:,:)
      real(8), allocatable :: ad_svstrG(:,:,:)

      real(8), allocatable :: ad_stflx(:,:,:)
!
#endif

      END MODULE mod_forces

