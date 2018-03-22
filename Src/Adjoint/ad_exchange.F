#include "cppdefs.h"



!==========================================

      SUBROUTINE ad_exchange_2d(a,LBi,UBi,LBj,UBj,Lm,Mm)
#if defined ADJOINT && defined NS_PERIODIC

      implicit none
      integer :: i,j,LBi,UBi,LBj,UBj,Lm,Mm
      real(8) :: a(LBi:UBi,LBj:UBj)

      DO i=LBi,UBi      
! tl line 5:       a(i,0)=a(i,Mm)
       a(i,Mm)=a(i,Mm)+a(i,0)
       a(i,0)=0.d0

! tl line 4:       a(i,-1)=a(i,Mm-1)
       a(i,Mm-1)=a(i,Mm-1)+a(i,-1)
       a(i,-1)=0.d0

! tl line 3:       a(i,-2)=a(i,Mm-2)
       a(i,Mm-2)=a(i,Mm-2)+a(i,-2)
       a(i,-2)=0.d0

! tl line 2:       a(i,Mm+2)=a(i,2)
       a(i,2)=a(i,2)+a(i,Mm+2)
       a(i,Mm+2)=0.d0

! tl line 1:       a(i,Mm+1)=a(i,1)
       a(i,1)=a(i,1)+a(i,Mm+1)
       a(i,Mm+1)=0.d0
            
      END DO
#endif /* defined ADJOINT && defined NS_PERIODIC */

      RETURN
      END SUBROUTINE ad_exchange_2d


!==========================================

      SUBROUTINE ad_exchange_3d(a,LBi,UBi,LBj,UBj,Nb,Ns,Lm,Mm)
#if defined ADJOINT && defined NS_PERIODIC

      implicit none
      integer :: i,j,k,LBi,UBi,LBj,UBj,Nb,Ns,Lm,Mm
      real(8) :: a(LBi:UBi,LBj:UBj,Nb:Ns)

      DO k=Nb,Ns
       DO i=LBi,UBi      

! tl line 5:       a(i,0,k)=a(i,Mm,k)
       a(i,Mm,k)=a(i,Mm,k)+a(i,0,k)
       a(i,0,k)=0.d0

! tl line 4:       a(i,-1,k)=a(i,Mm-1,k)
       a(i,Mm-1,k)=a(i,Mm-1,k)+a(i,-1,k)
       a(i,-1,k)=0.d0

! tl line 3:       a(i,-2,k)=a(i,Mm-2,k)
       a(i,Mm-2,k)=a(i,Mm-2,k)+a(i,-2,k)
       a(i,-2,k)=0.d0

! tl line 2:       a(i,Mm+2,k)=a(i,2,k)
       a(i,2,k)=a(i,2,k)+a(i,Mm+2,k)
       a(i,Mm+2,k)=0.d0

! tl line 1:       a(i,Mm+1,k)=a(i,1,k)
       a(i,1,k)=a(i,1,k)+a(i,Mm+1,k)
       a(i,Mm+1,k)=0.d0
            
       END DO
      END DO
#endif /* defined ADJOINT && defined NS_PERIODIC */

      RETURN
      END SUBROUTINE ad_exchange_3d

