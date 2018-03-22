#include "cppdefs.h"


!      INTERFACE exchange

!==========================================

      SUBROUTINE exchange_2d(a,LBi,UBi,LBj,UBj,Lm,Mm)
#if defined NS_PERIODIC

      implicit none
      integer :: i,j,LBi,UBi,LBj,UBj,Lm,Mm
      real(8) :: a(LBi:UBi,LBj:UBj)

      DO i=LBi,UBi      
       a(i,Mm+1)=a(i,1)
       a(i,Mm+2)=a(i,2)

       a(i,-2)=a(i,Mm-2)
       a(i,-1)=a(i,Mm-1)
       a(i,0)=a(i,Mm)
            

      END DO

#endif /* NS_PERIODIC */
      RETURN
      END SUBROUTINE exchange_2d

!==========================================

      SUBROUTINE exchange_2d_int(a,LBi,UBi,LBj,UBj,Lm,Mm)
#if defined NS_PERIODIC

      implicit none
      integer :: i,j,LBi,UBi,LBj,UBj,Lm,Mm
      integer :: a(LBi:UBi,LBj:UBj)

      DO i=LBi,UBi      

       a(i,Mm+1)=a(i,1)
       a(i,Mm+2)=a(i,2)

       a(i,-2)=a(i,Mm-2)
       a(i,-1)=a(i,Mm-1)
       a(i,0)=a(i,Mm)
            
      END DO

#endif /* NS_PERIODIC */
      RETURN
      END SUBROUTINE exchange_2d_int

!==========================================

      SUBROUTINE exchange_3d(a,LBi,UBi,LBj,UBj,Nb,Ns,Lm,Mm)
#if defined NS_PERIODIC

      implicit none
      integer :: i,j,k,LBi,UBi,LBj,UBj,Nb,Ns,Lm,Mm
      real(8) :: a(LBi:UBi,LBj:UBj,Nb:Ns)

      DO k=Nb,Ns
       DO i=LBi,UBi      

        a(i,Mm+1,k)=a(i,1,k)
        a(i,Mm+2,k)=a(i,2,k)

        a(i,-2,k)=a(i,Mm-2,k)
        a(i,-1,k)=a(i,Mm-1,k)
        a(i,0,k)=a(i,Mm,k)
            
       END DO
      END DO
#endif /* NS_PERIODIC */

      RETURN
      END SUBROUTINE exchange_3d

!      END INTERFACE exchange
