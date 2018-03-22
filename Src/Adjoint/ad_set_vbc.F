#include "cppdefs.h"

#ifdef ADJOINT
      SUBROUTINE ad_set_vbc
      USE mod_ocean

      implicit none       
      integer :: i,j

! Linear bottom stress:
      DO j=JstrU,JendU
      DO i=IstrU,IendU
! TL:       tl_bustr(i,j)=rdrg*tl_u(i,j,1,nrhs)
       ad_u(i,j,1,nrhs)=ad_u(i,j,1,nrhs)+rdrg*ad_bustr(i,j)
       ad_bustr(i,j)=0.d0
      END DO
      END DO

      DO j=JstrV,JendV
      DO i=IstrV,IendV
! TL:       tl_bvstr(i,j)=rdrg*tl_v(i,j,1,nrhs)
       ad_v(i,j,1,nrhs)=ad_v(i,j,1,nrhs)+rdrg*ad_bvstr(i,j)
       ad_bvstr(i,j)=0.d0
      END DO
      END DO

      RETURN
      END

#endif /* ADJOINT */
