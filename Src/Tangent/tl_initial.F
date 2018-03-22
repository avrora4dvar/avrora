#include "cppdefs.h"

#ifdef TANGENT


      SUBROUTINE tl_initial     
      USE mod_ocean
#ifdef MPI
      USE mod_parallel
#endif
      implicit none
      include 'netcdf.h'
#ifdef MPI
      include 'mpif.h'
      integer :: ierr,n1,n2,i1,i2
#endif
      integer :: recnum,nini,status,fid,varid,i,j,k,itr
      real(8), allocatable :: wrk3d(:,:,:)
      real(8) :: cff,tl_cff
      character (len=80) :: fhead

      recnum=1      
      nini=1

! read time:

#ifdef MPI
      IF (myrank.eq.0) THEN
#endif

      status=nf_open(TRIM(inifile),nf_nowrite,fid)
       status=nf_inq_varid(fid,'ocean_time',varid)
       status=nf_get_vara_double(fid,varid,recnum,1,time0)
      status=nf_close(fid)

#ifdef MPI
      END IF

      CALL mpi_bcast(time0,1,mpi_double_precision,0,
     *               mpi_comm_world,ierr)
#endif /* MPI */


!!! read ini tl fields:
#ifdef MPI

! - in rho locations:
      n1=UBi-LBi+1
      n2=UBj-LBj+1

      CALL get_2dfld_tile(inifile,'zeta',recnum,LBi+1,LBj+1,
     *                    n1,n2,tl_zeta(LBi:UBi,LBj:UBj,nini))

      CALL get_3dfld_tile(inifile,'temp',recnum,
     *                    LBi+1,LBj+1,1,n1,n2,N,     ! start (nc file), count
     *                    tl_t(LBi:UBi,LBj:UBj,1:N,nini,1))  

      CALL get_3dfld_tile(inifile,'salt',recnum,
     *                    LBi+1,LBj+1,1,n1,n2,N,     ! start (nc file), count
     *                    tl_t(LBi:UBi,LBj:UBj,1:N,nini,2))  

! - in u locations:
      i1=max(1,LBi)  ! index of 1st element on roms grid
      n1=UBi-i1+1
      n2=UBj-LBj+1

! note: first element is i1, not i1+1, since grid is staggered 
! (1st element in the grid file corresponds to element 1 on roms grid)
      CALL get_3dfld_tile(inifile,'u',recnum,
     *                    i1,LBj+1,1,n1,n2,N,     ! start (nc file), count
     *                    tl_u(i1:UBi,LBj:UBj,1:N,nini))  

! in v-locations:
      i2=max(1,LBj)
      n1=UBi-LBi+1
      n2=UBj-i2+1

      CALL get_3dfld_tile(inifile,'v',recnum,
     *                    LBi+1,i2,1,n1,n2,N,     ! start (nc file), count
     *                    tl_v(LBi:UBi,i2:UBj,1:N,nini))

#else /* serial */
      CALL get_2dfld(inifile,'zeta',recnum,Lp,Mp,
     *               tl_zeta(0:L,0:M,nini))
      CALL get_3dfld(inifile,'u'   ,recnum,L ,Mp,N,
     *               tl_u(1:L,0:M,1:N,nini))
      CALL get_3dfld(inifile,'v'   ,recnum,Lp,M ,N,
     *               tl_v(0:L,1:M,1:N,nini))
      CALL get_3dfld(inifile,'temp',recnum,Lp,Mp,N,
     *               tl_t(0:L,0:M,1:N,nini,1))
      CALL get_3dfld(inifile,'salt',recnum,Lp,Mp,N,
     *               tl_t(0:L,0:M,1:N,nini,2))

#endif /* MPI or serial */


! Masking:
      DO j=LBj,UBj
      DO i=LBi,UBi
       tl_zeta(i,j,nini)=tl_zeta(i,j,nini)*DBLE(mask_r(i,j))
       tl_u(i,j,1:N,nini)=tl_u(i,j,1:N,nini)*DBLE(mask_u(i,j))
       tl_v(i,j,1:N,nini)=tl_v(i,j,1:N,nini)*DBLE(mask_v(i,j))
       tl_t(i,j,1:N,nini,1:Ntr)=tl_t(i,j,1:N,nini,1:Ntr)*
     *                                       DBLE(mask_r(i,j)) 
      END DO
      END DO

#ifdef NS_PERIODIC

#ifdef MPI
      STOP('periodic BC not coded in MPI')
#else /* serial */
      CALL exchange_2d(tl_zeta(:,:,nini),LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_3d(tl_u(:,:,:,nini),LBi,UBi,LBj,UBj,1,N,Lm,Mm)
      CALL exchange_3d(tl_v(:,:,:,nini),LBi,UBi,LBj,UBj,1,N,Lm,Mm)
      CALL exchange_3d(tl_t(:,:,:,nini,1),LBi,UBi,LBj,UBj,1,N,Lm,Mm)
      CALL exchange_3d(tl_t(:,:,:,nini,2),LBi,UBi,LBj,UBj,1,N,Lm,Mm)
#endif /* MPI or serial */

#endif /* NS_PERIODIC */

!
      CALL tl_set_depth(LBi,UBi,LBj,UBj,N,Lm,Mm,
     *          z0_r,z0_w,oh,
     *          tl_zeta(:,:,nini),tl_z_r,tl_z_w,tl_Hz,tl_Hz_u,tl_Hz_v)

# ifdef TS_FIXED
      DO itr=1,Ntr
       DO k=1,N
        DO j=LBj,UBj
         DO i=LBi,UBi
          tl_t(i,j,k,2,itr)=tl_t(i,j,k,nini,itr)
          tl_t(i,j,k,3,itr)=tl_t(i,j,k,nini,itr)
         END DO
        END DO
       END DO
      END DO     
# endif

      RETURN
      END

#endif    /* TANGENT */   
