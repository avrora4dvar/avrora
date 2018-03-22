#include "cppdefs.h"


      SUBROUTINE set_bounds
      USE mod_ocean
#ifdef MPI
      USE mod_parallel
#endif

      implicit none
#ifdef MPI
      include 'mpif.h'
      include 'netcdf.h'
      integer, allocatable :: LBi_all(:),UBi_all(:),
     *                    LBj_all(:),UBj_all(:),
     *                    Istr_all(:),Iend_all(:), 
     *                    IstrU_all(:),IendU_all(:), 
     *                    IstrV_all(:),IendV_all(:), 
     *                    Jstr_all(:),Jend_all(:), 
     *                    JstrU_all(:),JendU_all(:), 
     *                    JstrV_all(:),JendV_all(:), 
     *                    c1_all(:),c2_all(:) 
#endif
      integer :: status,fid,dimid,ierr,nproc,i


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reading dimensions from netcdf file: only on 1 processor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      myrank=0
#ifdef MPI
      call mpi_comm_rank(mpi_comm_world,myrank,ierr)
#endif
      IF (myrank.eq.0) THEN

! Read dimensions from the grid file (Lp,Mp):
       status=nf_open(TRIM(grdfile),nf_nowrite,fid)
       status=nf_inq_dimid(fid, 'xi_rho', dimid)
       status=nf_inq_dimlen(fid,dimid,Lp)
       status=nf_inq_dimid(fid, 'eta_rho', dimid)
       status=nf_inq_dimlen(fid,dimid,Mp)
       status=nf_inq_dimid(fid, 'N', dimid)
       status=nf_inq_dimlen(fid,dimid,N)
       status=nf_close(fid)

      END IF ! IF (myrank.eq.0) THEN

#ifdef MPI
      call mpi_bcast(Lp,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(Mp,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(N ,1,mpi_integer,0,mpi_comm_world,ierr)
#endif
      L=Lp-1     ! following lines: all ranks in MPI regime
      M=Mp-1
      Lm=L-1
      Mm=M-1
      Np=N+1

! 8/23/2010: I have decided to initialy write MPI version
! only for the open boundary case 

#ifdef NS_PERIODIC /* serial only */
      LBi=0
      UBi=L
      LBj=-2
      UBj=Mp

!! Interior points limits (NS_PERIODIC):
! - rho:
      Istr=1
      Iend=Lm
      Jstr=1
      Jend=Mm

! - u:
      IstrU=2
      IendU=Lm
      JstrU=1
      JendU=Mm

! - v:
      IstrV=1
      IendV=Lm
      JstrV=1    ! <- additional v-row included compared to OBC case
      JendV=Mm

#else   /* ALL OPEN/CLOSED BOUNDARIES, MPI or serial */   
# ifdef MPI

      call mpi_comm_rank(comm2d,myrank,ierr)
      call mpi_cart_coords(comm2d,myrank,2,coords,ierr)

! Bounds w/out ghost points, each direction:
      call get_tile_bounds_noghost(LBi,UBi,Lp,coords(1),ntiles_xi)
      call get_tile_bounds_noghost(LBj,UBj,Mp,coords(2),ntiles_eta)

!! Interior points limits:
! - rho:
      Istr=max(1,LBi)
      Iend=min(Lm,UBi)
      Jstr=max(1,LBj)
      Jend=min(Mm,UBj)

! - u:
      IstrU=max(2,LBi)
      IendU=Iend
      JstrU=Jstr
      JendU=Jend

! - v:
      IstrV=Istr
      IendV=Iend
      JstrV=max(2,LBj)
      JendV=Jend

! Add ghost points (so far, in each direction, add 1 layer max)
      call ghost_2(LBi,UBi,LBj,UBj,coords,ntiles_xi,ntiles_eta)

! Logical variables (each side of the tile: boundary or not)
      B_North=.FALSE.
      B_South=.FALSE.
      B_East =.FALSE.
      B_West =.FALSE.
      if (coords(1).eq.0)           B_West=.TRUE.
      if (coords(1).eq.ntiles_xi-1) B_East=.TRUE.
      if (coords(2).eq.0)           B_South=.TRUE.
      if (coords(2).eq.ntiles_eta-1)B_North=.TRUE. 

! Point limits for exchange (== interior point limits + boundary points)
      Istr_E=Istr
      Iend_E=Iend
      Jstr_E=Jstr
      Jend_E=Jend

      IstrU_E=IstrU
      IendU_E=IendU
      JstrU_E=JstrU
      JendU_E=JendU

      IstrV_E=IstrV
      IendV_E=IendV
      JstrV_E=JstrV
      JendV_E=JendV

      IF (B_West) then
       Istr_E =Istr_E -1
       IstrU_E=IstrU_E-1
       IstrV_E=IstrV_E-1
      ENDIF

      IF (B_East) THEN
       Iend_E =Iend_E +1
       IendU_E=IendU_E+1
       IendV_E=IendV_E+1
      ENDIF

      IF (B_South) THEN
       Jstr_E =Jstr_E -1
       JstrU_E=JstrU_E-1
       JstrV_E=JstrV_E-1
      ENDIF

      IF (B_North) THEN
       Jend_E =Jend_E +1
       JendU_E=JendU_E+1
       JendV_E=JendV_E+1
      ENDIF

#ifdef ADJOINT
! extend each tile by N_EXT_ADJ points to facilitate computation in a larger domain
      LBi0=LBi
      UBi0=UBi
      LBj0=LBj
      UBj0=UBj

      if (.not.B_west) then
       LBi  =LBi  -n_ext_adj
       Istr =Istr -n_ext_adj
       IstrU=IstrU-n_ext_adj
       IstrV=IstrV-n_ext_adj
      endif
      if (.not.B_east) then
       UBi  =UBi  +n_ext_adj
       Iend =Iend +n_ext_adj
       IendU=IendU+n_ext_adj
       IendV=IendV+n_ext_adj
      endif
      if (.not.B_south) then
       LBj  =LBj  -n_ext_adj
       Jstr =Jstr -n_ext_adj
       JstrU=JstrU-n_ext_adj
       JstrV=JstrV-n_ext_adj
      endif
      if (.not.B_north) then
       UBj  =UBj  +n_ext_adj
       Jend =Jend +n_ext_adj
       JendU=JendU+n_ext_adj
       JendV=JendV+n_ext_adj
      endif
     
#endif

! Collect all indices on the master node, output to file
      nproc=ntiles_xi*ntiles_eta
      IF (myrank.eq.0) THEN
       allocate(LBi_all(0:nproc-1))
       allocate(UBi_all(0:nproc-1))
       allocate(LBj_all(0:nproc-1))
       allocate(UBj_all(0:nproc-1))
       allocate(Istr_all(0:nproc-1))
       allocate(Iend_all(0:nproc-1))
       allocate(IstrU_all(0:nproc-1))
       allocate(IendU_all(0:nproc-1))
       allocate(IstrV_all(0:nproc-1))
       allocate(IendV_all(0:nproc-1))
       allocate(Jstr_all(0:nproc-1))
       allocate(Jend_all(0:nproc-1))
       allocate(JstrU_all(0:nproc-1))
       allocate(JendU_all(0:nproc-1))
       allocate(JstrV_all(0:nproc-1))
       allocate(JendV_all(0:nproc-1))
       allocate(c1_all(0:nproc-1))
       allocate(c2_all(0:nproc-1))
      END IF
      call mpi_barrier(comm2d,ierr)

      call mpi_gather(LBi,1,MPI_INTEGER,LBi_all,1,
     *                MPI_INTEGER,0,mpi_comm_world,ierr)
      call mpi_gather(UBi,1,MPI_INTEGER,UBi_all,1,
     *                MPI_INTEGER,0,mpi_comm_world,ierr)
      call mpi_gather(LBj,1,MPI_INTEGER,LBj_all,1,
     *                MPI_INTEGER,0,mpi_comm_world,ierr)
      call mpi_gather(UBj,1,MPI_INTEGER,UBj_all,1,
     *                MPI_INTEGER,0,mpi_comm_world,ierr)

      call mpi_gather(Istr,1,MPI_INTEGER,Istr_all,1,
     *                MPI_INTEGER,0,mpi_comm_world,ierr)
      call mpi_gather(Iend,1,MPI_INTEGER,Iend_all,1,
     *                MPI_INTEGER,0,mpi_comm_world,ierr)
      call mpi_gather(IstrU,1,MPI_INTEGER,IstrU_all,1,
     *                MPI_INTEGER,0,mpi_comm_world,ierr)
      call mpi_gather(IendU,1,MPI_INTEGER,IendU_all,1,
     *                MPI_INTEGER,0,mpi_comm_world,ierr)
      call mpi_gather(IstrV,1,MPI_INTEGER,IstrV_all,1,
     *                MPI_INTEGER,0,mpi_comm_world,ierr)
      call mpi_gather(IendV,1,MPI_INTEGER,IendV_all,1,
     *                MPI_INTEGER,0,mpi_comm_world,ierr)

      call mpi_gather(Jstr,1,MPI_INTEGER,Jstr_all,1,
     *                MPI_INTEGER,0,mpi_comm_world,ierr)
      call mpi_gather(Jend,1,MPI_INTEGER,Jend_all,1,
     *                MPI_INTEGER,0,mpi_comm_world,ierr)
      call mpi_gather(JstrU,1,MPI_INTEGER,JstrU_all,1,
     *                MPI_INTEGER,0,mpi_comm_world,ierr)
      call mpi_gather(JendU,1,MPI_INTEGER,JendU_all,1,
     *                MPI_INTEGER,0,mpi_comm_world,ierr)
      call mpi_gather(JstrV,1,MPI_INTEGER,JstrV_all,1,
     *                MPI_INTEGER,0,mpi_comm_world,ierr)
      call mpi_gather(JendV,1,MPI_INTEGER,JendV_all,1,
     *                MPI_INTEGER,0,mpi_comm_world,ierr)

      call mpi_gather(coords(1),1,MPI_INTEGER,c1_all,1,
     *                MPI_INTEGER,0,mpi_comm_world,ierr)
      call mpi_gather(coords(2),1,MPI_INTEGER,c2_all,1,
     *                MPI_INTEGER,0,mpi_comm_world,ierr)

      IF (myrank.eq.0) THEN
       open(unit=91,file='tiles.dat')
       write(91,'(18a6)') '% c1','c2','LBi','UBi','LBj','UBj',
     *                   'Istr',' Iend','Jstr','Jend ',
     *                   'IstrU','IendU','JstrU','JendU ',
     *                   'IstrV','IendV','JstrV','JendV'
       do i=0,nproc-1
        write (91,'(18i6)') c1_all(i),c2_all(i),
     *               LBi_all(i),UBi_all(i),
     *               LBj_all(i),UBj_all(i),
     *               Istr_all(i),Iend_all(i),
     *               Jstr_all(i),Jend_all(i),
     *               IstrU_all(i),IendU_all(i),
     *               JstrU_all(i),JendU_all(i),
     *               IstrV_all(i),IendV_all(i),
     *               JstrV_all(i),JendV_all(i)
       end do
       close(91)
      END IF

 
# else  /* not mpi */  

!! Array limits (for each process):
      LBi=0
      UBi=L
      LBj=0
      UBj=M

!! Interior points limits:
! - rho:
      Istr=1
      Iend=Lm
      Jstr=1
      Jend=Mm

! - u:
      IstrU=2
      IendU=Lm
      JstrU=1
      JendU=Mm

! - v:
      IstrV=1
      IendV=Lm
      JstrV=2
      JendV=Mm
# endif /* ifdef MPI ... else ... */

#endif /* ifdef NS_PERIODIC ... else ... */
      
      RETURN
      END

#ifdef MPI
      SUBROUTINE get_tile_bounds_noghost(LBi,UBi,Lp,coords,ntiles)

      implicit none
      integer :: LBi,UBi,Lp,coords,ntiles,nmod,n1,ncells

      nmod=mod(Lp,ntiles)
      ncells=Lp/ntiles
      n1=ntiles-nmod
      LBi=ncells*coords
      UBi=ncells*(coords+1)-1
      IF (coords.ge.n1) THEN
       LBi=LBi+coords-n1
       UBi=UBi+coords-n1+1
      END IF

      RETURN
      END

      SUBROUTINE ghost_2(LBi,UBi,LBj,UBj,coords,ntiles_xi,ntiles_eta)
      
      implicit none
      integer :: LBi,UBi,LBj,UBj,coords(2),ntiles_xi,ntiles_eta

      IF (coords(1).ne.0) LBi=LBi-2
      IF (coords(1).ne.ntiles_xi-1) UBi=UBi+2
      IF (coords(2).ne.0) LBj=LBj-2
      IF (coords(2).ne.ntiles_eta-1) UBj=UBj+2

      RETURN
      END

!      SUBROUTINE ghost(LBi,UBi,LBj,UBj,coords,ntiles_xi,ntiles_eta)
!      
!      implicit none
!      integer :: LBi,UBi,LBj,UBj,coords(2),ntiles_xi,ntiles_eta
!
!      IF (coords(1).ne.0) LBi=LBi-1
!      IF (coords(1).ne.ntiles_xi-1) UBi=UBi+1
!      IF (coords(2).ne.0) LBj=LBj-1
!      IF (coords(2).ne.ntiles_eta-1) UBj=UBj+1
!
!      RETURN
!      END

#endif /* MPI */
