#include "cppdefs.h"

#ifdef ADJOINT
!=====================================================================
#ifdef MPI
!=====================================================================
      SUBROUTINE ad_t3dbc(Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj,N,Ntr,  
     *             B_West,B_East,B_South,B_North, ! <- new line, comp ot v2.1
     *             nrhs, nstp, nnew,itrc,ad_t)

      implicit none
      integer, intent(in) :: Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj,N,Ntr
      integer, intent(in) :: nrhs, nstp, nnew, itrc   
      real(8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N,3,Ntr)

      logical, intent(in) :: B_West,B_East,B_South,B_North

      integer :: i,j,k

# ifndef NS_PERIODIC
! DEFAULT CONDITIONS: WALL (NO GRADIENT CONDITION)
      IF (B_North) THEN
       DO k=1,N
        DO i=LBi,UBi
! tl:        tl_t(i,Jend+1,k,nnew,itrc)=tl_t(i,Jend,k,nnew,itrc) !  NORTH
         ad_t(i,Jend,k,nnew,itrc)=ad_t(i,Jend,k,nnew,itrc)+
     *                           ad_t(i,Jend+1,k,nnew,itrc)
         ad_t(i,Jend+1,k,nnew,itrc)=0.d0
        END DO
       END DO
      END IF

      IF (B_South) THEN
       DO k=1,N
        DO i=LBi,UBi
! tl:        tl_t(i,Jstr-1,k,nnew,itrc)=tl_t(i,Jstr,k,nnew,itrc) !  SOUTH
        ad_t(i,Jstr,k,nnew,itrc)=ad_t(i,Jstr,k,nnew,itrc)+
     *                           ad_t(i,Jstr-1,k,nnew,itrc)
        ad_t(i,Jstr-1,k,nnew,itrc)=0.d0

        END DO
       END DO
      END IF
# endif

! DEFAULT CONDITIONS: WALL (NO GRADIENT CONDITION)
      IF (B_East) THEN
       DO k=1,N
        DO j=LBj,UBj

! tl:        tl_t(Iend+1,j,k,nnew,itrc)=tl_t(Iend,j,k,nnew,itrc) ! EAST
         ad_t(Iend,j,k,nnew,itrc)=ad_t(Iend,j,k,nnew,itrc)+
     *                           ad_t(Iend+1,j,k,nnew,itrc)
         ad_t(Iend+1,j,k,nnew,itrc)=0.d0
        END DO      
       END DO
      END IF

      IF (B_West) THEN
       DO k=1,N
        DO j=LBj,UBj
! tl:        tl_t(Istr-1,j,k,nnew,itrc)=tl_t(Istr,j,k,nnew,itrc) ! WEST
         ad_t(Istr,j,k,nnew,itrc)=ad_t(Istr,j,k,nnew,itrc)+
     *                           ad_t(Istr-1,j,k,nnew,itrc)
         ad_t(Istr-1,j,k,nnew,itrc)=0.d0
        END DO      
       END DO
      END IF



      RETURN
      END


!=====================================================================
#else /* serial */
!=====================================================================
      SUBROUTINE ad_t3dbc(Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj,N,Ntr,  
     *             nrhs, nstp, nnew,itrc,ad_t)

      implicit none
      integer, intent(in) :: Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj,N,Ntr
      integer, intent(in) :: nrhs, nstp, nnew, itrc   
      real(8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N,3,Ntr)

      integer :: i,j,k

# ifndef NS_PERIODIC
! DEFAULT CONDITIONS: WALL (NO GRADIENT CONDITION)
      DO k=1,N
       DO i=LBi,UBi

! tl:        tl_t(i,Jend+1,k,nnew,itrc)=tl_t(i,Jend,k,nnew,itrc) !  NORTH
        ad_t(i,Jend,k,nnew,itrc)=ad_t(i,Jend,k,nnew,itrc)+
     *                           ad_t(i,Jend+1,k,nnew,itrc)
        ad_t(i,Jend+1,k,nnew,itrc)=0.d0

! tl:        tl_t(i,Jstr-1,k,nnew,itrc)=tl_t(i,Jstr,k,nnew,itrc) !  SOUTH
        ad_t(i,Jstr,k,nnew,itrc)=ad_t(i,Jstr,k,nnew,itrc)+
     *                           ad_t(i,Jstr-1,k,nnew,itrc)
        ad_t(i,Jstr-1,k,nnew,itrc)=0.d0

       END DO
      END DO
# endif

! DEFAULT CONDITIONS: WALL (NO GRADIENT CONDITION)
      DO k=1,N
       DO j=LBj,UBj

! tl:        tl_t(Iend+1,j,k,nnew,itrc)=tl_t(Iend,j,k,nnew,itrc) ! EAST
        ad_t(Iend,j,k,nnew,itrc)=ad_t(Iend,j,k,nnew,itrc)+
     *                           ad_t(Iend+1,j,k,nnew,itrc)
        ad_t(Iend+1,j,k,nnew,itrc)=0.d0

! tl:        tl_t(Istr-1,j,k,nnew,itrc)=tl_t(Istr,j,k,nnew,itrc) ! WEST
        ad_t(Istr,j,k,nnew,itrc)=ad_t(Istr,j,k,nnew,itrc)+
     *                           ad_t(Istr-1,j,k,nnew,itrc)
        ad_t(Istr-1,j,k,nnew,itrc)=0.d0

       END DO      
      END DO



      RETURN
      END
!=====================================================================
#endif /* MPI or serial */
!=====================================================================

#endif /* ADJOINT */
