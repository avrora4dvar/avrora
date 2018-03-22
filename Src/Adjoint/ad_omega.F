#include "cppdefs.h"

# ifdef ADJOINT
      SUBROUTINE ad_omega
      USE mod_ocean
  
      implicit none

      integer :: i,j,k
      real(8) :: wrk(LBi:UBi),ad_wrk(LBi:UBi)
      real(8) :: tmp(LBi:UBi),ad_tmp,tmp1

      ad_tmp=0.d0
      ad_wrk=0.d0

# if defined NS_PERIODIC
      CALL ad_exchange_3d(ad_W,LBi,UBi,LBj,UBj,0,N,Lm,Mm)
# endif 

      DO j=Jstr,Jend

!-----------------------------------------------------------------------
! ADJ to: Account for moving S-surfaces
!-----------------------------------------------------------------------

! Background wrk and tmp used throughout:
        DO i=Istr,Iend
         tmp(i)=z_w(i,j,N)-z_w(i,j,0)
         wrk(i)=WN(i,j)/tmp(i)
        END DO

        DO i=Istr,Iend
          ad_W(i,j,N)=0.d0
        END DO

        DO k=1,N-1
         DO i=Istr,Iend
! TL (recurs): tl_W(i,j,k)=tl_W(i,j,k)-tl_wrk(i)*(z_w(i,j,k)-z_w(i,j,0))
! TL:     *                -wrk(i)*(tl_z_w(i,j,k)-tl_z_w(i,j,0))

          tmp1=ad_W(i,j,k)*wrk(i)

          ad_wrk(i)=ad_wrk(i)-(z_w(i,j,k)-z_w(i,j,0))*ad_W(i,j,k)
          ad_z_w(i,j,k)=ad_z_w(i,j,k)-tmp1
          ad_z_w(i,j,0)=ad_z_w(i,j,0)+tmp1
          
         END DO
        END DO

!-----------------------------------------------------------------------
!  ADJ to: Vertically integrate horizontal mass flux divergence.
!-----------------------------------------------------------------------

! TL loop 3:
        DO i=Istr,Iend

! tl line 2:         tl_wrk(i)=(tl_W(i,j,N)-wrk(i)*tl_tmp)/tmp {(i)}

         tmp1=ad_wrk(i)/tmp(i)

         ad_W(i,j,N)=ad_W(i,j,N)+tmp1
         ad_tmp=ad_tmp-tmp1*wrk(i)
         ad_wrk(i)=0.d0

! tl line 1:         tl_tmp=tl_z_w(i,j,N)-tl_z_w(i,j,0)
         ad_z_w(i,j,N)=ad_z_w(i,j,N)+ad_tmp
         ad_z_w(i,j,0)=ad_z_w(i,j,0)-ad_tmp
         ad_tmp=0.d0

        END DO

! TL loop 2:
        DO k=N,1,-1
         DO i=Istr,Iend
! TL (recurs):          tl_W(i,j,k)=tl_W(i,j,k-1)- 
! TL:    *             (tl_Hudy(i+1,j,k)-tl_Hudy(i,j,k)+
! TL:    *              tl_Hvdx(i,j+1,k)-tl_Hvdx(i,j,k))

          ad_W(i,j,k-1)=ad_W(i,j,k-1)+ad_W(i,j,k)

          ad_Hudy(i+1,j,k)=ad_Hudy(i+1,j,k)-ad_W(i,j,k)
          ad_Hudy(i,j,k)=ad_Hudy(i,j,k)+ad_W(i,j,k)
          ad_Hvdx(i,j+1,k)=ad_Hvdx(i,j+1,k)-ad_W(i,j,k)
          ad_Hvdx(i,j,k)=ad_Hvdx(i,j,k)+ad_W(i,j,k)

          ad_W(i,j,k)=0.d0

         END DO
        END DO

! TL loop 1:
        DO i=Istr,Iend
          ad_W(i,j,0)=0.d0
        END DO


      END DO ! j=..

      RETURN
      END

#endif /* ADJOINT */
