#include "cppdefs.h"

!************************
      SUBROUTINE mean_field_2d(LBi,UBi,LBj,UBj,q,q0)

      implicit none
      integer :: LBi,UBi,LBj,UBj,i,j
      real(8) :: q(LBi:UBi,LBj:UBj)
      real(8) :: q0

      q0=0.d0
      DO j=LBj,UBj
      DO i=LBi,UBi
       q0=q0+q(i,j)
      END DO
      END DO
      q0=q0/((UBi-LBi+1)*(UBj-LBj+1))

      RETURN
      END

!************************
      SUBROUTINE max_field_2d(LBi,UBi,LBj,UBj,q,q0)

      implicit none
      integer :: LBi,UBi,LBj,UBj,i,j
      real(8) :: q(LBi:UBi,LBj:UBj)
      real(8) :: q0

      q0=0.d0
      DO j=LBj,UBj
      DO i=LBi,UBi
       IF (q(i,j)*q(i,j).gt.q0*q0) THEN
        q0=q(i,j)
       ENDIF
      END DO
      END DO

      RETURN
      END
      
