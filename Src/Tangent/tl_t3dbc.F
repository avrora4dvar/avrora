#include "cppdefs.h"

#ifdef TANGENT
!=====================================================================
#ifdef MPI
!=====================================================================
      SUBROUTINE tl_t3dbc(Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj,N,Ntr,  
     *             B_West,B_East,B_South,B_North, ! <- new line, comp ot v2.1
     *             nrhs, nstp, nnew,itrc,tl_t)

      implicit none
      integer, intent(in) :: Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj,N,Ntr
      integer, intent(in) :: nrhs, nstp, nnew, itrc   
      real(8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N,3,Ntr)

      logical, intent(in) :: B_West,B_East,B_South,B_North

      integer :: i,j,k

! DEFAULT CONDITIONS: WALL (NO GRADIENT CONDITION)
      IF (B_West) THEN
       DO k=1,N
        DO j=LBj,UBj
         tl_t(Istr-1,j,k,nnew,itrc)=tl_t(Istr,j,k,nnew,itrc) ! WEST
        END DO      
       END DO
      END IF

      IF (B_East) THEN
       DO k=1,N
        DO j=LBj,UBj
         tl_t(Iend+1,j,k,nnew,itrc)=tl_t(Iend,j,k,nnew,itrc) ! EAST
        END DO      
       END DO
      END IF

# ifndef NS_PERIODIC
! DEFAULT CONDITIONS: WALL (NO GRADIENT CONDITION)

      IF (B_South) THEN
       DO k=1,N
        DO i=LBi,UBi
         tl_t(i,Jstr-1,k,nnew,itrc)=tl_t(i,Jstr,k,nnew,itrc) !  SOUTH
        END DO
       END DO
      END IF

      IF (B_North) THEN
       DO k=1,N
        DO i=LBi,UBi
         tl_t(i,Jend+1,k,nnew,itrc)=tl_t(i,Jend,k,nnew,itrc) !  NORTH
        END DO
       END DO
      END IF

# endif

      RETURN
      END
!=====================================================================
#else /* serial */
!=====================================================================
      SUBROUTINE tl_t3dbc(Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj,N,Ntr,  
     *             nrhs, nstp, nnew,itrc,tl_t)

      implicit none
      integer, intent(in) :: Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj,N,Ntr
      integer, intent(in) :: nrhs, nstp, nnew, itrc   
      real(8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N,3,Ntr)

      integer :: i,j,k

! DEFAULT CONDITIONS: WALL (NO GRADIENT CONDITION)
      DO k=1,N
       DO j=LBj,UBj
        tl_t(Istr-1,j,k,nnew,itrc)=tl_t(Istr,j,k,nnew,itrc) ! WEST
        tl_t(Iend+1,j,k,nnew,itrc)=tl_t(Iend,j,k,nnew,itrc) ! EAST
       END DO      
      END DO

# ifndef NS_PERIODIC
! DEFAULT CONDITIONS: WALL (NO GRADIENT CONDITION)
      DO k=1,N
       DO i=LBi,UBi
        tl_t(i,Jstr-1,k,nnew,itrc)=tl_t(i,Jstr,k,nnew,itrc) !  SOUTH
        tl_t(i,Jend+1,k,nnew,itrc)=tl_t(i,Jend,k,nnew,itrc) !  NORTH
       END DO
      END DO
# endif

      RETURN
      END

!=====================================================================
#endif /* MPI or serial */
!=====================================================================

#endif /* TANGENT */
