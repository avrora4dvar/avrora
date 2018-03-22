#include "cppdefs.h"


#ifdef ADJOINT
      SUBROUTINE ad_zetabc(Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj,  
     *             krhs, kstp, knew,ad_zeta)

      implicit none
      integer, intent(in) :: Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj      
      integer, intent(in) :: krhs, kstp, knew   
      real(8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)

      integer :: i,j

# ifndef NS_PERIODIC
! DEFAULT CONDITIONS: WALL
      DO i=LBi,UBi
! TL line 2:       tl_zeta(i,Jend+1,knew)=tl_zeta(i,Jend,knew) !  NORTH
       ad_zeta(i,Jend,knew)=ad_zeta(i,Jend,knew)+
     *                        ad_zeta(i,Jend+1,knew)
       ad_zeta(i,Jend+1,knew)=0.d0

! TL line 1:       tl_zeta(i,Jstr-1,knew)=tl_zeta(i,Jstr,knew) !  SOUTH
       ad_zeta(i,Jstr,knew)=ad_zeta(i,Jstr,knew)+
     *                        ad_zeta(i,Jstr-1,knew)
       ad_zeta(i,Jstr-1,knew)=0.d0
      END DO      
# endif

! DEFAULT CONDITIONS: WALL
      DO j=LBj,UBj
! TL line 2:       tl_zeta(Iend+1,j,knew)=tl_zeta(Iend,j,knew) !  EAST
       ad_zeta(Iend,j,knew)=ad_zeta(Iend,j,knew)+
     *                        ad_zeta(Iend+1,j,knew) 
       ad_zeta(Iend+1,j,knew)=0.d0


! TL line 1:       tl_zeta(Istr-1,j,knew)=tl_zeta(Istr,j,knew) !  WEST
       ad_zeta(Istr,j,knew)=ad_zeta(Istr,j,knew)+
     *                        ad_zeta(Istr-1,j,knew)
       ad_zeta(Istr-1,j,knew)=0.d0
      END DO      


      RETURN
      END
#endif /* TANGENT */
