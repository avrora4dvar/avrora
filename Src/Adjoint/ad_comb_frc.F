#include "cppdefs.h"

#ifdef ADJOINT
      SUBROUTINE ad_comb_frc
      USE mod_ocean
#ifdef MPI
      USE mod_parallel
#endif

      implicit none
      INCLUDE 'netcdf.h'

      integer :: recnum,status,varid,fid,i,j,k
      real(8) :: time_comb
      real(8), allocatable :: wrk2d(:,:)
      real(8), allocatable :: wrk3d(:,:,:)
#ifdef MPI
      include 'mpif.h'
      integer :: ierr,n1,n2,i1,i2
#endif

! Adds snapshots of the comb (data functional) to the state
      IF (MOD(iic-1,nhis).eq.0) THEN
       recnum=(iic-1)/nhis+1.1  ! note: recnum is integer
       
#ifdef MPI
      IF (myrank.eq.0) THEN
#endif

! read and chk time:
       status=nf_open(TRIM(hisfile),nf_nowrite,fid)
        status=nf_inq_varid(fid,'ocean_time',varid)
        status=nf_get_vara_double(fid,varid,recnum,1,time_comb)
       status=nf_close(fid)

#ifdef MPI
      END IF

      CALL mpi_bcast(time_comb,1,mpi_double_precision,0,
     *               mpi_comm_world,ierr)
#endif /* MPI */
 
       IF (DABS(time_comb-time) .gt. 1.d-03*dt) THEN
        write (*,*) 'STOP in ad_comb_frc:  time_comb ~= time'
        write (*,*) 'time_comb=',time_comb
        write (*,*) 'time=',time
        STOP
       END IF

!
! ADD TO AD_ZETA(:,:,1), time index would corr. to that in tl_output.F:
! 
      
       allocate(wrk2d(LBi:UBi,LBj:UBj))   
       wrk2d(LBi:UBi,LBj:UBj)=0.d0
#ifdef MPI
! - in rho locations:
      n1=UBi-LBi+1
      n2=UBj-LBj+1
      CALL get_2dfld_tile(hisfile,'zeta',recnum,LBi+1,LBj+1,
     *                    n1,n2,wrk2d(LBi:UBi,LBj:UBj))

#else /* serial */
       CALL get_2dfld(hisfile,'zeta',recnum,Lp,Mp,
     *                wrk2d(0:L,0:M))
#endif /* MPI or serial */

       DO j=LBj,UBj
       DO i=LBi,UBi
        ad_zeta(i,j,1)=ad_zeta(i,j,1)+wrk2d(i,j)
       END DO
       END DO
       deallocate(wrk2d)
!
! ADD TO 3D VARS:
!
       allocate(wrk3d(LBi:UBi,LBj:UBj,N))   

! ADD TO AD_U:
       wrk3d(LBi:UBi,LBj:UBj,1:N)=0.d0
#ifdef MPI
! - in u locations:
      i1=max(1,LBi)  ! index of 1st element on roms grid
      n1=UBi-i1+1
      n2=UBj-LBj+1
! note: first element is i1, not i1+1, since grid is staggered 
! (1st element in the grid file corresponds to element 1 on roms grid)
      CALL get_3dfld_tile(hisfile,'u',recnum,
     *                    i1,LBj+1,1,n1,n2,N,     ! start (nc file), count
     *                    wrk3d(i1:UBi,LBj:UBj,1:N))  
 
#else /* serial */
       CALL get_3dfld(hisfile,'u'   ,recnum,L ,Mp,N,
     *                wrk3d(1:L,0:M,1:N))
#endif /* MPI or serial */

       DO k=1,N
       DO j=LBj,UBj
       DO i=LBi,UBi
        ad_u(i,j,k,nstp)=ad_u(i,j,k,nstp)+wrk3d(i,j,k)
       END DO
       END DO
       END DO

! ADD TO AD_V:
       wrk3d(LBi:UBi,LBj:UBj,1:N)=0.d0

#ifdef MPI
! in v-locations:
      i2=max(1,LBj)
      n1=UBi-LBi+1
      n2=UBj-i2+1

      CALL get_3dfld_tile(hisfile,'v',recnum,
     *                    LBi+1,i2,1,n1,n2,N,     ! start (nc file), count
     *                    wrk3d(LBi:UBi,i2:UBj,1:N))
#else /* serial */
       CALL get_3dfld(hisfile,'v'   ,recnum,Lp,M ,N,
     *                wrk3d(0:L,1:M,1:N))
#endif /* MPI or serial */

       DO k=1,N
       DO j=LBj,UBj
       DO i=LBi,UBi
        ad_v(i,j,k,nstp)=ad_v(i,j,k,nstp)+wrk3d(i,j,k)
       END DO
       END DO
       END DO

! ADD TO AD_TEMP:


       wrk3d(LBi:UBi,LBj:UBj,1:N)=0.d0
#ifdef MPI
! - in rho locations:
      n1=UBi-LBi+1
      n2=UBj-LBj+1
      CALL get_3dfld_tile(hisfile,'temp',recnum,
     *                    LBi+1,LBj+1,1,n1,n2,N,     ! start (nc file), count
     *                    wrk3d(LBi:UBi,LBj:UBj,1:N))  
#else /* serial */
       CALL get_3dfld(hisfile,'temp',recnum,Lp,Mp,N,
     *                wrk3d(0:L,0:M,1:N))
#endif /* MPI or serial */

       DO k=1,N
       DO j=LBj,UBj
       DO i=LBi,UBi
        ad_t(i,j,k,nstp,1)=ad_t(i,j,k,nstp,1)+wrk3d(i,j,k)
       END DO
       END DO
       END DO

! ADD TO AD_SALT:
       wrk3d(LBi:UBi,LBj:UBj,1:N)=0.d0
#ifdef MPI
! - in rho locations:
      n1=UBi-LBi+1
      n2=UBj-LBj+1
      CALL get_3dfld_tile(hisfile,'salt',recnum,
     *                    LBi+1,LBj+1,1,n1,n2,N,     ! start (nc file), count
     *                    wrk3d(LBi:UBi,LBj:UBj,1:N))  
#else /* serial */
       CALL get_3dfld(hisfile,'salt',recnum,Lp,Mp,N,
     *                wrk3d(0:L,0:M,1:N))
#endif /* MPI or serial */
       DO k=1,N
       DO j=LBj,UBj
       DO i=LBi,UBi
        ad_t(i,j,k,nstp,2)=ad_t(i,j,k,nstp,2)+wrk3d(i,j,k)
       END DO
       END DO
       END DO
       
       deallocate(wrk3d)

      END IF

      RETURN
      END

#endif /* ADJOINT */
