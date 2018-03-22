#include "cppdefs.h"

#ifdef TANGENT

!=====================================================================
#ifdef MPI
!=====================================================================
      SUBROUTINE tl_u3dbc(IstrU,IendU,JstrU,JendU,LBi,UBi,LBj,UBj,N,  
     *             B_West,B_East,B_South,B_North, ! <- new line, comp ot v2.1
     *             nrhs, nstp, nnew,tl_u)

      implicit none
      integer, intent(in) :: IstrU,IendU,JstrU,JendU,LBi,UBi,LBj,UBj,N      
      integer, intent(in) :: nrhs, nstp, nnew   
      real(8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N,2)

      logical, intent(in) :: B_West,B_East,B_South,B_North

      integer :: i,j,k

! DEFAULT CONDITIONS: WALL
      IF (B_West) THEN
       DO k=1,N
        DO j=LBj,UBj
         tl_u(IstrU-1,j,k,nnew)=0.d0 ! WEST
        END DO      
       END DO
      END IF ! B_West

      IF (B_east) THEN
       DO k=1,N
        DO j=LBj,UBj
         tl_u(IendU+1,j,k,nnew)=0.d0 ! EAST
        END DO      
       END DO
      END IF ! B_East

# ifndef NS_PERIODIC
! DEFAULT CONDITIONS: WALL (this is for tangential velocity)

      IF (B_South) THEN
       DO k=1,N
        DO i=LBi,UBi
         tl_u(i,JstrU-1,k,nnew)=tl_u(i,JstrU,k,nnew) !  SOUTH
        END DO
       END DO
      END IF

      IF (B_North) THEN
       DO k=1,N
        DO i=LBi,UBi
         tl_u(i,JendU+1,k,nnew)=tl_u(i,JendU,k,nnew) !  NORTH
        END DO
       END DO
      END IF

# endif

      RETURN
      END

!****************************************************************
      SUBROUTINE tl_v3dbc(IstrV,IendV,JstrV,JendV,LBi,UBi,LBj,UBj,N,
     *             B_West,B_East,B_South,B_North, ! <- new line, comp ot v2.1
     *             nrhs, nstp, nnew,tl_v)

      implicit none
      integer, intent(in) :: IstrV,IendV,JstrV,JendV,LBi,UBi,LBj,UBj,N   
      integer, intent(in) :: nrhs, nstp, nnew   
      real(8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N,2)

      logical, intent(in) :: B_West,B_East,B_South,B_North

      integer :: i,j,k

! DEFAULT CONDITIONS: WALL (this is for tangential velocity)
      IF (B_West) THEN
       DO k=1,N
        DO j=LBj,UBj
         tl_v(IstrV-1,j,k,nnew)=tl_v(IstrV,j,k,nnew) ! WEST
        END DO      
       END DO
      END IF

      IF (B_East) THEN
       DO k=1,N
        DO j=LBj,UBj
         tl_v(IendV+1,j,k,nnew)=tl_v(IendV,j,k,nnew) ! EAST
        END DO      
       END DO
      END IF

# ifndef NS_PERIODIC
! DEFAULT CONDITIONS: WALL 
      IF (B_South) THEN
       DO k=1,N
        DO i=LBi,UBi
         tl_v(i,JstrV-1,k,nnew)=0.d0 !  SOUTH
        END DO     
       END DO
      END IF

      IF (B_North) THEN
       DO k=1,N
        DO i=LBi,UBi
         tl_v(i,JendV+1,k,nnew)=0.d0 !  NORTH
        END DO     
       END DO
      END IF
# endif

      RETURN
      END

!=====================================================================
#else /* serial */
!=====================================================================
      SUBROUTINE tl_u3dbc(IstrU,IendU,JstrU,JendU,LBi,UBi,LBj,UBj,N,  
     *             nrhs, nstp, nnew,tl_u)

      implicit none
      integer, intent(in) :: IstrU,IendU,JstrU,JendU,LBi,UBi,LBj,UBj,N      
      integer, intent(in) :: nrhs, nstp, nnew   
      real(8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N,2)

      integer :: i,j,k

! DEFAULT CONDITIONS: WALL
      DO k=1,N
       DO j=LBj,UBj
        tl_u(IstrU-1,j,k,nnew)=0.d0 ! WEST
        tl_u(IendU+1,j,k,nnew)=0.d0 ! EAST
       END DO      
      END DO

# ifndef NS_PERIODIC
! DEFAULT CONDITIONS: WALL (this is for tangential velocity)
      DO k=1,N
       DO i=LBi,UBi
        tl_u(i,JstrU-1,k,nnew)=tl_u(i,JstrU,k,nnew) !  SOUTH
        tl_u(i,JendU+1,k,nnew)=tl_u(i,JendU,k,nnew) !  NORTH
       END DO
      END DO
# endif

      RETURN
      END

!****************************************************************
      SUBROUTINE tl_v3dbc(IstrV,IendV,JstrV,JendV,LBi,UBi,LBj,UBj,N,
     *             nrhs, nstp, nnew,tl_v)

      implicit none
      integer, intent(in) :: IstrV,IendV,JstrV,JendV,LBi,UBi,LBj,UBj,N   
      integer, intent(in) :: nrhs, nstp, nnew   
      real(8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N,2)

      integer :: i,j,k

! DEFAULT CONDITIONS: WALL (this is for tangential velocity)
      DO k=1,N
       DO j=LBj,UBj
        tl_v(IstrV-1,j,k,nnew)=tl_v(IstrV,j,k,nnew) ! WEST
        tl_v(IendV+1,j,k,nnew)=tl_v(IendV,j,k,nnew) ! EAST
       END DO      
      END DO

# ifndef NS_PERIODIC
! DEFAULT CONDITIONS: WALL 
      DO k=1,N
       DO i=LBi,UBi
        tl_v(i,JstrV-1,k,nnew)=0.d0 !  SOUTH
        tl_v(i,JendV+1,k,nnew)=0.d0 !  NORTH
       END DO     
      END DO
# endif

      RETURN
      END
!=====================================================================
#endif /* MPI or SERIAL */
!=====================================================================

#endif /* TANGENT   */
