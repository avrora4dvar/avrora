#include "cppdefs.h"

! A. Kurapov, 4/14/2007: 
! this subroutine 
! checks the length of forcing time series, 
! allocates forcing time vectors, 
! reads forcing times (such as sms_time)

! MPI version: 9/7/2010

#ifdef MPI

      SUBROUTINE read_frc_time
      USE mod_ocean
      USE mod_forces
      USE mod_parallel

      implicit none
      include 'mpif.h'
      include 'netcdf.h'
      integer :: fid,varid,dimid,status,ierr

      character (len=14) :: ftmp
      integer :: i

      IF (myrank.eq.0) THEN
! Wind stress time  - on processor rank=0
       status=nf_open(TRIM(frcfile),nf_nowrite,fid)
        status=nf_inq_varid(fid,'sms_time',varid)
        status=nf_inq_vardimid(fid,varid,dimid)
        status=nf_inq_dimlen(fid,dimid,nsms)

        allocate(sms_time(nsms))
        status=nf_get_vara_double(fid,varid,1,nsms,sms_time)     
 
       status=nf_close(fid)
      END  IF ! myrank=0

      CALL mpi_bcast(nsms,1,mpi_integer,0,mpi_comm_world,ierr)

      IF (myrank.ne.0) THEN
        allocate(sms_time(nsms))
      END  IF 

      CALL mpi_bcast(sms_time,nsms,mpi_double_precision,0,
     *               mpi_comm_world,ierr)
 

      RETURN
      END 

#else /* serial */

      SUBROUTINE read_frc_time
      USE mod_ocean
      USE mod_forces
      USE mod_nf_param

      implicit none
      integer :: fid,varid,dimid,status

! Wind stress time
      status=nf_open(TRIM(frcfile),nf_nowrite,fid)
      status=nf_inq_varid(fid,'sms_time',varid)
      status=nf_inq_vardimid(fid,varid,dimid)
      status=nf_inq_dimlen(fid,dimid,nsms)

      allocate(sms_time(nsms))
      status=nf_get_vara_double(fid,varid,1,nsms,sms_time)     

      status=nf_close(fid)

      RETURN
      END 
#endif /* MPI or serial */
