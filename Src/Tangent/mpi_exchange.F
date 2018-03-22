#include "cppdefs.h"

#ifdef ALLOW_MPI_EXCH_1 /* switch to explude 1-ghost mpi_exch from v3.2 */
      subroutine mpi_exchange_2d(a,LBi,UBi,LBj,UBj,
     *                     i1,i2,j1,j2,
     *                     comm2d, stridetype, 
     *                     nbrleft, nbrright, nbrtop, nbrbottom) 

! 1 ghost layer on each side of a tile

! Exchanges elements of a(LBi:UBi,LBj:UBj)
! i1,i2,j1,j2: points of each tile to exchange 
! (interior+domain boundary)

#ifdef MPI
 
      implicit none
      include 'mpif.h' 
      integer LBi,UBi,LBj,UBj,i1,i2,j1,j2 
      double precision a(LBi:UBi,LBj:UBj) 
      integer nbrleft,nbrright,nbrtop,nbrbottom,comm2d,stridetype 
      integer status(MPI_STATUS_SIZE), ierr, nx
c 
      nx = i2 - i1 + 1 
c  These are just like the 1-d versions, except for less data 
      call MPI_SENDRECV( a(i1,j2),  nx, MPI_DOUBLE_PRECISION,  
     &                    nbrtop, 0,  
     &                    a(i1,LBj), nx, MPI_DOUBLE_PRECISION,  
     &                    nbrbottom, 0, comm2d, status, ierr ) 
      call MPI_SENDRECV( a(i1,j1),  nx, MPI_DOUBLE_PRECISION, 
     &                    nbrbottom, 1,  
     &                    a(i1,UBj), nx, MPI_DOUBLE_PRECISION,  
     &                    nbrtop, 1, comm2d, status, ierr ) 
c 
c This uses the vector datatype stridetype 
      call MPI_SENDRECV( a(i2,LBj),  1, stridetype, nbrright, 0,  
     &                     a(LBi,LBj), 1, stridetype, nbrleft, 0, 
     &                     comm2d, status, ierr ) 
      call MPI_SENDRECV( a(i1,LBj),  1, stridetype, nbrleft,   1, 
     &                     a(UBi,LBj), 1, stridetype, nbrright, 1, 
     &                     comm2d, status, ierr ) 

#endif
      return 
      end 

!====================================================================
!====================================================================
      subroutine mpi_exchange_3d(a,LBi,UBi,LBj,UBj,N,
     *                     i1,i2,j1,j2,
     *                     comm2d, stride_bt,stride_lr, 
     *                     nbrleft, nbrright, nbrtop, nbrbottom) 

! Exchanges elements of a(LBi:UBi,LBj:UBj,1:N)
! i1,i2,j1,j2: points of each tile to exchange 
! (interior+domain boundary)

#ifdef MPI
 
      implicit none
      include 'mpif.h' 
      integer :: LBi,UBi,LBj,UBj,i1,i2,j1,j2,N 
      real(8) :: a(LBi:UBi,LBj:UBj,N) 
      integer :: nbrleft,nbrright,nbrtop,nbrbottom,comm2d
      integer :: stride_bt,stride_lr 
      integer :: status(MPI_STATUS_SIZE), ierr

      call MPI_SENDRECV( a(i1,j2,1), 1, stride_bt,nbrtop, 0,  
     &                   a(i1,LBj,1),1, stride_bt,nbrbottom,0, 
     &                   comm2d, status, ierr ) 
      call MPI_SENDRECV( a(i1,j1,1),   1, stride_bt,nbrbottom,1,   
     &                   a(i1,UBj,1), 1, stride_bt,nbrtop, 1, 
     &                   comm2d, status, ierr ) 
c 
c This uses the vector datatype stridetype 
      call MPI_SENDRECV( a(i2,LBj,1),  1, stride_lr, nbrright, 0,  
     &                     a(LBi,LBj,1), 1, stride_lr, nbrleft, 0, 
     &                     comm2d, status, ierr ) 
      call MPI_SENDRECV( a(i1,LBj,1),  1, stride_lr, nbrleft,   1, 
     &                     a(UBi,LBj,1), 1, stride_lr, nbrright, 1, 
     &                     comm2d, status, ierr ) 

#endif
      return 
      end 
#endif /* ALLOW_MPI_EXCH_1 */

!====================================================================
!====================================================================
      subroutine mpi_exchange_2d_2(a,LBi,UBi,LBj,UBj,
     *                     i1,i2,j1,j2,
     *                     comm2d, stridetype, 
     *                     nbrleft, nbrright, nbrtop, nbrbottom) 

! New compared to mpi_exchange_2d: 2 ghost layers on each side of a tile

! Exchanges elements of a(LBi:UBi,LBj:UBj)
! i1,i2,j1,j2: points of each tile to exchange 
! (interior+domain boundary). 

#ifdef MPI
 
      implicit none
      include 'mpif.h' 
      integer LBi,UBi,LBj,UBj,i1,i2,j1,j2 
      double precision a(LBi:UBi,LBj:UBj) 
      integer nbrleft,nbrright,nbrtop,nbrbottom,comm2d,stridetype 
      integer status(MPI_STATUS_SIZE), ierr, nx
c 
      nx = i2 - i1 + 1 

      call MPI_SENDRECV( a(i1,j2),  nx, MPI_DOUBLE_PRECISION,  
     &                    nbrtop, 0,  
     &                    a(i1,LBj+1), nx, MPI_DOUBLE_PRECISION,  
     &                    nbrbottom, 0, comm2d, status, ierr ) 
      call MPI_SENDRECV( a(i1,j1+1),  nx, MPI_DOUBLE_PRECISION, 
     &                    nbrbottom, 1,  
     &                    a(i1,UBj), nx, MPI_DOUBLE_PRECISION,  
     &                    nbrtop, 1, comm2d, status, ierr ) 
      call MPI_SENDRECV( a(i1,j2-1),  nx, MPI_DOUBLE_PRECISION,  
     &                    nbrtop, 0,  
     &                    a(i1,LBj), nx, MPI_DOUBLE_PRECISION,  
     &                    nbrbottom, 0, comm2d, status, ierr ) 
      call MPI_SENDRECV( a(i1,j1),  nx, MPI_DOUBLE_PRECISION, 
     &                    nbrbottom, 1,  
     &                    a(i1,UBj-1), nx, MPI_DOUBLE_PRECISION,  
     &                    nbrtop, 1, comm2d, status, ierr ) 
c 
c This uses the vector datatype stridetype 
      call MPI_SENDRECV( a(i2,LBj),  1, stridetype, nbrright, 0,  
     &                     a(LBi+1,LBj), 1, stridetype, nbrleft, 0, 
     &                     comm2d, status, ierr ) 
      call MPI_SENDRECV( a(i1+1,LBj),  1, stridetype, nbrleft,   1, 
     &                     a(UBi,LBj), 1, stridetype, nbrright, 1, 
     &                     comm2d, status, ierr ) 
      call MPI_SENDRECV( a(i2-1,LBj),  1, stridetype, nbrright, 0,  
     &                     a(LBi,LBj), 1, stridetype, nbrleft, 0, 
     &                     comm2d, status, ierr ) 
      call MPI_SENDRECV( a(i1,LBj),  1, stridetype, nbrleft,   1, 
     &                     a(UBi-1,LBj), 1, stridetype, nbrright, 1, 
     &                     comm2d, status, ierr ) 

#endif
      return 
      end 
!====================================================================
!====================================================================
      subroutine mpi_exchange_3d_2(a,LBi,UBi,LBj,UBj,N,
     *                     i1,i2,j1,j2,
     *                     comm2d, stride_bt,stride_lr, 
     *                     nbrleft, nbrright, nbrtop, nbrbottom) 

! New compared to mpi_exchange_2d: 2 ghost layers on each side of a tile

! Exchanges elements of a(LBi:UBi,LBj:UBj,1:N)
! i1,i2,j1,j2: points of each tile to exchange 
! (interior+domain boundary)

#ifdef MPI
      implicit none
      include 'mpif.h' 
      integer :: LBi,UBi,LBj,UBj,i1,i2,j1,j2,N 
      real(8) :: a(LBi:UBi,LBj:UBj,N) 
      integer :: nbrleft,nbrright,nbrtop,nbrbottom,comm2d
      integer :: stride_bt,stride_lr 
      integer :: status(MPI_STATUS_SIZE), ierr

      call MPI_SENDRECV( a(i1,j2,1), 1, stride_bt,nbrtop, 0,  
     &                   a(i1,LBj+1,1),1, stride_bt,nbrbottom,0, 
     &                   comm2d, status, ierr ) 
      call MPI_SENDRECV( a(i1,j1+1,1),   1, stride_bt,nbrbottom,1,   
     &                   a(i1,UBj,1), 1, stride_bt,nbrtop, 1, 
     &                   comm2d, status, ierr ) 
      call MPI_SENDRECV( a(i1,j2-1,1), 1, stride_bt,nbrtop, 0,  
     &                   a(i1,LBj,1),1, stride_bt,nbrbottom,0, 
     &                   comm2d, status, ierr ) 
      call MPI_SENDRECV( a(i1,j1,1),   1, stride_bt,nbrbottom,1,   
     &                   a(i1,UBj-1,1), 1, stride_bt,nbrtop, 1, 
     &                   comm2d, status, ierr ) 
c 
c This uses the vector datatype stridetype 
      call MPI_SENDRECV( a(i2,LBj,1),  1, stride_lr, nbrright, 0,  
     &                     a(LBi+1,LBj,1), 1, stride_lr, nbrleft, 0, 
     &                     comm2d, status, ierr ) 
      call MPI_SENDRECV( a(i1+1,LBj,1),  1, stride_lr, nbrleft,   1, 
     &                     a(UBi,LBj,1), 1, stride_lr, nbrright, 1, 
     &                     comm2d, status, ierr ) 
      call MPI_SENDRECV( a(i2-1,LBj,1),  1, stride_lr, nbrright, 0,  
     &                     a(LBi,LBj,1), 1, stride_lr, nbrleft, 0, 
     &                     comm2d, status, ierr ) 
      call MPI_SENDRECV( a(i1,LBj,1),  1, stride_lr, nbrleft,   1, 
     &                     a(UBi-1,LBj,1), 1, stride_lr, nbrright, 1, 
     &                     comm2d, status, ierr ) 

#endif
      return 
      end 

!====================================================================
!====================================================================
! in adjoint, outer bounds of tiles are further extended than in the 
! tl

#ifdef ADJOINT
      subroutine mpi_exchange_2d_ext(a,LBi,UBi,LBj,UBj,
     *                     LBi0,UBi0,LBj0,UBj0,
     *                     i1,i2,j1,j2,
     *                     comm2d, stridetype, 
     *                     nbrleft, nbrright, nbrtop, nbrbottom) 

! New compared to mpi_exchange_2d: 2+n_ext_adj ghost layers 
! on each side of a tile. 

! Exchanges elements of a(LBi:UBi,LBj:UBj)
! i1,i2,j1,j2: extended interior = interior+domain boundary (excl. ghost pnts)
! LBi0,UBi0,LBj0,UBj0: limits of the TL tiles (incl. 2 layers of ghost points)  

#ifdef MPI
 
      implicit none
      include 'mpif.h' 
      integer LBi,UBi,LBj,UBj,LBi0,UBi0,LBj0,UBj0,i1,i2,j1,j2 
      double precision a(LBi:UBi,LBj:UBj) 
      integer nbrleft,nbrright,nbrtop,nbrbottom,comm2d,stridetype 
      integer status(MPI_STATUS_SIZE), ierr, nx
c 
      nx = i2 - i1 + 1 

      call MPI_SENDRECV( a(i1,j2),  nx, MPI_DOUBLE_PRECISION,  
     &                    nbrtop, 0,  
     &                    a(i1,j1-1), nx, MPI_DOUBLE_PRECISION,  
     &                    nbrbottom, 0, comm2d, status, ierr ) 
      call MPI_SENDRECV( a(i1,j1+1),  nx, MPI_DOUBLE_PRECISION, 
     &                    nbrbottom, 1,  
     &                    a(i1,j2+2), nx, MPI_DOUBLE_PRECISION,  
     &                    nbrtop, 1, comm2d, status, ierr ) 
      call MPI_SENDRECV( a(i1,j2-1),  nx, MPI_DOUBLE_PRECISION,  
     &                    nbrtop, 0,  
     &                    a(i1,j1-2), nx, MPI_DOUBLE_PRECISION,  
     &                    nbrbottom, 0, comm2d, status, ierr ) 
      call MPI_SENDRECV( a(i1,j1),  nx, MPI_DOUBLE_PRECISION, 
     &                    nbrbottom, 1,  
     &                    a(i1,j2+1), nx, MPI_DOUBLE_PRECISION,  
     &                    nbrtop, 1, comm2d, status, ierr ) 
!      call MPI_SENDRECV( a(i1,j2-2),  nx, MPI_DOUBLE_PRECISION,  
!     &                    nbrtop, 0,  
!     &                    a(i1,j1-3), nx, MPI_DOUBLE_PRECISION,  
!     &                    nbrbottom, 0, comm2d, status, ierr ) 
!      call MPI_SENDRECV( a(i1,j1+2),  nx, MPI_DOUBLE_PRECISION, 
!     &                    nbrbottom, 1,  
!     &                    a(i1,j2+3), nx, MPI_DOUBLE_PRECISION,  
!     &                    nbrtop, 1, comm2d, status, ierr ) 
c 
c This uses the vector datatype stridetype 
      call MPI_SENDRECV( a(i2  ,LBj),  1, stridetype, nbrright, 0,  
     &                   a(i1-1,LBj),  1, stridetype, nbrleft, 0, 
     &                     comm2d, status, ierr ) 
      call MPI_SENDRECV( a(i1+1,LBj),  1, stridetype, nbrleft,   1, 
     &                   a(i2+2,LBj),   1, stridetype, nbrright, 1, 
     &                     comm2d, status, ierr ) 
      call MPI_SENDRECV( a(i2-1,LBj),  1, stridetype, nbrright, 0,  
     &                   a(i1-2,LBj), 1, stridetype, nbrleft, 0, 
     &                     comm2d, status, ierr ) 
      call MPI_SENDRECV( a(i1  ,LBj),  1, stridetype, nbrleft,   1, 
     &                   a(i2+1,LBj), 1, stridetype, nbrright, 1, 
     &                     comm2d, status, ierr ) 
!      call MPI_SENDRECV( a(i2-2,LBj),  1, stridetype, nbrright, 0,  
!     &                   a(i1-3,LBj), 1, stridetype, nbrleft, 0, 
!     &                     comm2d, status, ierr ) 
!      call MPI_SENDRECV( a(i1+2,LBj),  1, stridetype, nbrleft,   1, 
!     &                   a(i2+3,LBj), 1, stridetype, nbrright, 1, 
!     &                     comm2d, status, ierr ) 

#endif
      return 
      end 

!====================================================================
      subroutine mpi_exchange_3d_ext(a,LBi,UBi,LBj,UBj,N,
     *                     LBi0,UBi0,LBj0,UBj0,
     *                     i1,i2,j1,j2,
     *                     comm2d, stride_bt,stride_lr, 
     *                     nbrleft, nbrright, nbrtop, nbrbottom) 

! New compared to mpi_exchange_2d: 2 ghost layers on each side of a tile

! Exchanges elements of a(LBi:UBi,LBj:UBj,1:N)
! i1,i2,j1,j2: extended interior = interior+domain boundary (excl. ghost pnts)
! LBi0,UBi0,LBj0,UBj0: limits of the TL tiles (incl. 2 layers of ghost points)  

#ifdef MPI
      implicit none
      include 'mpif.h' 
      integer :: LBi,UBi,LBj,UBj,LBi0,UBi0,LBj0,UBj0,i1,i2,j1,j2,N 
      real(8) :: a(LBi:UBi,LBj:UBj,N) 
      integer :: nbrleft,nbrright,nbrtop,nbrbottom,comm2d
      integer :: stride_bt,stride_lr 
      integer :: status(MPI_STATUS_SIZE), ierr

      call MPI_SENDRECV( a(i1,j2  ,1), 1, stride_bt,nbrtop, 0,  
     &                   a(i1,j1-1,1), 1, stride_bt,nbrbottom,0, 
     &                   comm2d, status, ierr ) 
      call MPI_SENDRECV( a(i1,j1+1,1), 1, stride_bt,nbrbottom,1,   
     &                   a(i1,j2+2,1), 1, stride_bt,nbrtop, 1, 
     &                   comm2d, status, ierr ) 
      call MPI_SENDRECV( a(i1,j2-1,1), 1, stride_bt,nbrtop, 0,  
     &                   a(i1,j1-2,1), 1, stride_bt,nbrbottom,0, 
     &                   comm2d, status, ierr ) 
      call MPI_SENDRECV( a(i1,j1  ,1), 1, stride_bt,nbrbottom,1,   
     &                   a(i1,j2+1,1), 1, stride_bt,nbrtop, 1, 
     &                   comm2d, status, ierr ) 
!      call MPI_SENDRECV( a(i1,j2-2,1), 1, stride_bt,nbrtop, 0,  
!     &                   a(i1,j1-3,1), 1, stride_bt,nbrbottom,0, 
!     &                   comm2d, status, ierr ) 
!      call MPI_SENDRECV( a(i1,j1+2,1), 1, stride_bt,nbrbottom,1,   
!     &                   a(i1,j2+3,1), 1, stride_bt,nbrtop, 1, 
!     &                   comm2d, status, ierr ) 
c 
c This uses the vector datatype stridetype 
      call MPI_SENDRECV( a(i2  ,LBj,1),  1, stride_lr, nbrright, 0,  
     &                   a(i1-1,LBj,1),  1, stride_lr, nbrleft, 0, 
     &                     comm2d, status, ierr ) 
      call MPI_SENDRECV( a(i1+1,LBj,1),  1, stride_lr, nbrleft,   1, 
     &                   a(i2+2,LBj,1),  1, stride_lr, nbrright, 1, 
     &                     comm2d, status, ierr ) 
      call MPI_SENDRECV( a(i2-1,LBj,1),  1, stride_lr, nbrright, 0,  
     &                   a(i1-2,LBj,1),  1, stride_lr, nbrleft, 0, 
     &                     comm2d, status, ierr ) 
      call MPI_SENDRECV( a(i1  ,LBj,1),  1, stride_lr, nbrleft,   1, 
     &                   a(i2+1,LBj,1),  1, stride_lr, nbrright, 1, 
     &                     comm2d, status, ierr ) 
!      call MPI_SENDRECV( a(i2-2,LBj,1),  1, stride_lr, nbrright, 0,  
!     &                   a(i1-3,LBj,1),  1, stride_lr, nbrleft, 0, 
!     &                     comm2d, status, ierr ) 
!      call MPI_SENDRECV( a(i1+2,LBj,1),  1, stride_lr, nbrleft,   1, 
!     &                   a(i2+3,LBj,1),  1, stride_lr, nbrright, 1, 
!     &                     comm2d, status, ierr ) 



#endif
      return 
      end 

#endif /* ADJOINT */
