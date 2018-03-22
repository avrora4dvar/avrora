#include "cppdefs.h"


#ifdef TANGENT
      SUBROUTINE tl_zetabc(Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj,  
     *             krhs, kstp, knew,tl_zeta)

      implicit none
      integer, intent(in) :: Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj      
      integer, intent(in) :: krhs, kstp, knew   
      real(8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)

      integer :: i,j

! DEFAULT CONDITIONS: WALL
      DO j=LBj,UBj
       tl_zeta(Istr-1,j,knew)=tl_zeta(Istr,j,knew) !  WEST
       tl_zeta(Iend+1,j,knew)=tl_zeta(Iend,j,knew) !  EAST
      END DO      

# ifndef NS_PERIODIC
! DEFAULT CONDITIONS: WALL
      DO i=LBi,UBi
       tl_zeta(i,Jstr-1,knew)=tl_zeta(i,Jstr,knew) !  SOUTH
       tl_zeta(i,Jend+1,knew)=tl_zeta(i,Jend,knew) !  NORTH
      END DO      
# endif

      RETURN
      END
#endif /* TANGENT */
