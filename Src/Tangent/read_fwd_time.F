#include "cppdefs.h"

#if defined TANGENT || defined ADJOINT

#ifdef MPI

      SUBROUTINE read_fwd_time

      USE mod_ocean
      USE mod_parallel

      implicit none
      include 'mpif.h'
      include 'netcdf.h'
      integer :: fid,varid,dimid,status,ierr

      character (len=14) :: ftmp
      integer :: i

      IF (myrank.eq.0) THEN
! Wind stress time dim - on processor rank=0
       status=nf_open(TRIM(fwdfile),nf_nowrite,fid)
       status=nf_inq_varid(fid,'ocean_time',varid)
       status=nf_inq_vardimid(fid,varid,dimid)
       status=nf_inq_dimlen(fid,dimid,nfwd)
      END  IF ! myrank=0

      CALL mpi_bcast(nfwd,1,mpi_integer,0,mpi_comm_world,ierr)

      allocate(fwd_time(nfwd))

      IF (myrank.eq.0) THEN
! Wind stress time  - on processor rank=0
       status=nf_get_vara_double(fid,varid,1,nfwd,fwd_time)      
       status=nf_close(fid)
      END  IF 

      CALL mpi_bcast(fwd_time,nfwd,mpi_double_precision,0,
     *               mpi_comm_world,ierr)

      RETURN
      END
 
#else /* serial */ 
      SUBROUTINE read_fwd_time

! A. Kurapov, 5/7/2007: 
! this subroutine 
! checks the length of the background time series, 
! allocates the background time vector, 
! reads background times (in fwd_time)

 
      USE mod_ocean
      USE mod_nf_param

      implicit none
      integer :: fid,varid,dimid,status

      status=nf_open(TRIM(fwdfile),nf_nowrite,fid)
      status=nf_inq_varid(fid,'ocean_time',varid)
      status=nf_inq_vardimid(fid,varid,dimid)
      status=nf_inq_dimlen(fid,dimid,nfwd)

      allocate(fwd_time(nfwd))
      status=nf_get_vara_double(fid,varid,1,nfwd,fwd_time)     

      status=nf_close(fid)
      
      RETURN
      END

#endif /* MPI or serial */
 
#endif /* defined TANGENT || defined ADJOINT */
