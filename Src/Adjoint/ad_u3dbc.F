#include "cppdefs.h"

#ifdef ADJOINT

!=====================================================================
#ifdef MPI
!=====================================================================
      SUBROUTINE ad_u3dbc(IstrU,IendU,JstrU,JendU,LBi,UBi,LBj,UBj,N,  
     *             B_West,B_East,B_South,B_North, ! <- new line, comp ot v2.1
     *             nrhs, nstp, nnew,ad_u)

      implicit none
      integer, intent(in) :: IstrU,IendU,JstrU,JendU,LBi,UBi,LBj,UBj,N      
      integer, intent(in) :: nrhs, nstp, nnew   
      real(8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N,2)

      logical, intent(in) :: B_West,B_East,B_South,B_North

      integer :: i,j,k

# ifndef NS_PERIODIC

! DEFAULT CONDITIONS: WALL (this is for tangential velocity)

      IF (B_North) THEN
       DO k=1,N
        DO i=LBi,UBi
! TL loop 4:        tl_u(i,JendU+1,k,nnew)=tl_u(i,JendU,k,nnew) !  NORTH
         ad_u(i,JendU,k,nnew)=ad_u(i,JendU,k,nnew)+
     *                         ad_u(i,JendU+1,k,nnew)
         ad_u(i,JendU+1,k,nnew)=0.d0
        END DO
       END DO
      END IF

      IF (B_South) THEN
       DO k=1,N
        DO i=LBi,UBi
! TL loop 3:        tl_u(i,JstrU-1,k,nnew)=tl_u(i,JstrU,k,nnew) !  SOUTH
         ad_u(i,JstrU,k,nnew)=ad_u(i,JstrU,k,nnew)+
     *                         ad_u(i,JstrU-1,k,nnew)
         ad_u(i,JstrU-1,k,nnew)=0.d0
        END DO
       END DO
      END IF

# endif

! DEFAULT CONDITIONS: WALL
      IF (B_east) THEN
       DO k=1,N
        DO j=LBj,UBj
! TL loop 2:        tl_u(IendU+1,j,k,nnew)=0.d0 ! EAST
         ad_u(IendU+1,j,k,nnew)=0.d0 
        END DO      
       END DO
      END IF ! B_East

      IF (B_West) THEN
       DO k=1,N
        DO j=LBj,UBj
! TL loop 1:        tl_u(IstrU-1,j,k,nnew)=0.d0 ! WEST
         ad_u(IstrU-1,j,k,nnew)=0.d0
        END DO      
       END DO
      END IF

      RETURN
      END

!****************************************************************
      SUBROUTINE ad_v3dbc(IstrV,IendV,JstrV,JendV,LBi,UBi,LBj,UBj,N,
     *             B_West,B_East,B_South,B_North, ! <- new line, comp ot v2.1
     *             nrhs, nstp, nnew,ad_v)

      implicit none
      integer, intent(in) :: IstrV,IendV,JstrV,JendV,LBi,UBi,LBj,UBj,N   
      integer, intent(in) :: nrhs, nstp, nnew   
      real(8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N,2)

      logical, intent(in) :: B_West,B_East,B_South,B_North

      integer :: i,j,k

# ifndef NS_PERIODIC
! DEFAULT CONDITIONS: WALL 
      IF (B_North) THEN
       DO k=1,N
        DO i=LBi,UBi
! TL loop 4:        tl_v(i,JendV+1,k,nnew)=0.d0 !  NORTH
         ad_v(i,JendV+1,k,nnew)=0.d0
        END DO     
       END DO
      END IF

      IF (B_South) THEN
       DO k=1,N
        DO i=LBi,UBi
! TL loop 3:        tl_v(i,JstrV-1,k,nnew)=0.d0 !  SOUTH
         ad_v(i,JstrV-1,k,nnew)=0.d0
        END DO     
       END DO
      END IF
# endif

! DEFAULT CONDITIONS: WALL (this is for tangential velocity)
      IF (B_East) THEN
       DO k=1,N
        DO j=LBj,UBj
! TL loop 2:        tl_v(IendV+1,j,k,nnew)=tl_v(IendV,j,k,nnew) ! EAST
         ad_v(IendV,j,k,nnew)=ad_v(IendV,j,k,nnew)+
     *                         ad_v(IendV+1,j,k,nnew)
         ad_v(IendV+1,j,k,nnew)=0.d0
        END DO      
       END DO
      END IF

      IF (B_West) THEN
       DO k=1,N
        DO j=LBj,UBj
! TL loop 1:        tl_v(IstrV-1,j,k,nnew)=tl_v(IstrV,j,k,nnew) ! WEST
         ad_v(IstrV,j,k,nnew)=ad_v(IstrV,j,k,nnew)+
     *                         ad_v(IstrV-1,j,k,nnew)
         ad_v(IstrV-1,j,k,nnew)=0.d0
        END DO      
       END DO
      END IF


      RETURN
      END


!=====================================================================
#else /* serial */
!=====================================================================
      SUBROUTINE ad_u3dbc(IstrU,IendU,JstrU,JendU,LBi,UBi,LBj,UBj,N,  
     *             nrhs, nstp, nnew,ad_u)

      implicit none
      integer, intent(in) :: IstrU,IendU,JstrU,JendU,LBi,UBi,LBj,UBj,N      
      integer, intent(in) :: nrhs, nstp, nnew   
      real(8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N,2)

      integer :: i,j,k

# ifndef NS_PERIODIC
! DEFAULT CONDITIONS: WALL (this is for tangential velocity)
      DO k=1,N
       DO i=LBi,UBi
! TL line 2:        tl_u(i,JendU+1,k,nnew)=tl_u(i,JendU,k,nnew) !  NORTH
        ad_u(i,JendU,k,nnew)=ad_u(i,JendU,k,nnew)+
     *                         ad_u(i,JendU+1,k,nnew)
        ad_u(i,JendU+1,k,nnew)=0.d0
! TL line 1:        tl_u(i,JstrU-1,k,nnew)=tl_u(i,JstrU,k,nnew) !  SOUTH
        ad_u(i,JstrU,k,nnew)=ad_u(i,JstrU,k,nnew)+
     *                         ad_u(i,JstrU-1,k,nnew)
        ad_u(i,JstrU-1,k,nnew)=0.d0
       END DO
      END DO
# endif

! DEFAULT CONDITIONS: WALL
      DO k=1,N
       DO j=LBj,UBj
! TL line 2:        tl_u(IendU+1,j,k,nnew)=0.d0 ! EAST
        ad_u(IendU+1,j,k,nnew)=0.d0 
! TL line 1:        tl_u(IstrU-1,j,k,nnew)=0.d0 ! WEST
        ad_u(IstrU-1,j,k,nnew)=0.d0
       END DO      
      END DO

      RETURN
      END

!****************************************************************
      SUBROUTINE ad_v3dbc(IstrV,IendV,JstrV,JendV,LBi,UBi,LBj,UBj,N,
     *             nrhs, nstp, nnew,ad_v)

      implicit none
      integer, intent(in) :: IstrV,IendV,JstrV,JendV,LBi,UBi,LBj,UBj,N   
      integer, intent(in) :: nrhs, nstp, nnew   
      real(8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N,2)

      integer :: i,j,k

# ifndef NS_PERIODIC
! DEFAULT CONDITIONS: WALL 
      DO k=1,N
       DO i=LBi,UBi
! TL line 2:        tl_v(i,JendV+1,k,nnew)=0.d0 !  NORTH
        ad_v(i,JendV+1,k,nnew)=0.d0
! TL line 1:        tl_v(i,JstrV-1,k,nnew)=0.d0 !  SOUTH
        ad_v(i,JstrV-1,k,nnew)=0.d0
       END DO     
      END DO
# endif

! DEFAULT CONDITIONS: WALL (this is for tangential velocity)
      DO k=1,N
       DO j=LBj,UBj
! TL line 2:        tl_v(IendV+1,j,k,nnew)=tl_v(IendV,j,k,nnew) ! EAST
        ad_v(IendV,j,k,nnew)=ad_v(IendV,j,k,nnew)+
     *                         ad_v(IendV+1,j,k,nnew)
        ad_v(IendV+1,j,k,nnew)=0.d0
! TL line 1:        tl_v(IstrV-1,j,k,nnew)=tl_v(IstrV,j,k,nnew) ! WEST
        ad_v(IstrV,j,k,nnew)=ad_v(IstrV,j,k,nnew)+
     *                         ad_v(IstrV-1,j,k,nnew)
        ad_v(IstrV-1,j,k,nnew)=0.d0
       END DO      
      END DO

      RETURN
      END
!=====================================================================
#endif /* MPI or serial */
!=====================================================================

#endif /* ADJOINT */
