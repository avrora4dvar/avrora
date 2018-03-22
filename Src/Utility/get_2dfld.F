#include "cppdefs.h"

#ifdef MPI
      SUBROUTINE get_2dfld_tile(fname,varname,recnum,istr1,istr2,
     *                          nwrk1,nwrk2,wrk)
      USE mod_parallel

      implicit none
      include 'netcdf.h'
      include 'mpif.h'

      character(len=*) :: fname,varname
      integer :: istr1,istr2,nwrk1,nwrk2
      real(8) :: wrk(nwrk1,nwrk2)
      integer :: recnum,fid,varid,status,ndims,i,len,n1,n2
      integer,allocatable :: start(:),count(:)
      integer, allocatable :: dimids(:),dimlens(:)
      integer :: ierr

      DO irank=0,comm_size-1
      CALL mpi_barrier(comm2d,ierr)
      IF(myrank.eq.irank) THEN

      status=nf_open(TRIM(fname),nf_nowrite,fid)

       ! inquire about the dimensions of the 2D variable 
       status=nf_inq_varid(fid,TRIM(varname),varid)
       status=nf_inq_varndims(fid,varid,ndims)

       IF (ndims.gt.3) THEN
        write (*,*) 'STOP in get_2dfld, since ndims>3'
        write (*,*) 'fname: ',fname
        write (*,*) 'varname: ',varname
        STOP
       END IF

       allocate (dimids(ndims))
       allocate (dimlens(ndims))           
       allocate (start(ndims))
       allocate (count(ndims))
       status=nf_inq_vardimid(fid,varid,dimids)

       DO i=1,ndims
        status=nf_inq_dimlen(fid,dimids(i),dimlens(i))
       END DO

       n1=dimlens(1)
       n2=dimlens(2)

       IF (istr1.lt.1 .or. istr2.lt.1 .or. 
     *     istr1+nwrk1-1.gt.n1 .or. istr2+nwrk2-1.gt.n2) THEN
        write (*,*) 'STOP in get_2dfld: inconsistent dims'
        STOP
       END IF

       start(1)=istr1
       start(2)=istr2
       count(1)=nwrk1
       count(2)=nwrk2
       IF (ndims.eq.3) THEN 
        start(3)=recnum
        count(3)=1
       END IF
       status=nf_get_vara_double(fid,varid,start,count,wrk)


      status=nf_close(fid)

      deallocate(dimids,dimlens,count,start)

      END IF
      END DO

      RETURN
      END

!******************
      SUBROUTINE get_2dfld_int_tile(fname,varname,recnum,istr1,istr2,
     *                          nwrk1,nwrk2,wrk)

      USE mod_parallel

      implicit none

      include 'mpif.h'
      include 'netcdf.h'

      character(len=*) :: fname,varname
      integer :: istr1,istr2,nwrk1,nwrk2
      integer :: wrk(nwrk1,nwrk2)
      integer :: recnum,fid,varid,status,ndims,i,len,n1,n2
      integer,allocatable :: start(:),count(:)
      integer, allocatable :: dimids(:),dimlens(:)
      integer :: ierr

      DO irank=0,comm_size-1
      CALL mpi_barrier(comm2d,ierr)
      IF(myrank.eq.irank) THEN

      status=nf_open(TRIM(fname),nf_nowrite,fid)

       ! inquire about the dimensions of the 2D variable 
       status=nf_inq_varid(fid,TRIM(varname),varid)
       status=nf_inq_varndims(fid,varid,ndims)

       IF (ndims.gt.3) THEN
        write (*,*) 'STOP in get_2dfld, since ndims>3'
        STOP
       END IF

       allocate (dimids(ndims))
       allocate (dimlens(ndims))           
       allocate (start(ndims))
       allocate (count(ndims))
       status=nf_inq_vardimid(fid,varid,dimids)

       DO i=1,ndims
        status=nf_inq_dimlen(fid,dimids(i),dimlens(i))
       END DO

       n1=dimlens(1)
       n2=dimlens(2)

       IF (istr1.lt.1 .or. istr2.lt.1 .or. 
     *     istr1+nwrk1-1.gt.n1 .or. istr2+nwrk2-1.gt.n2) THEN
        write (*,*) 'STOP in get_2dfld: inconsistent dims'
        STOP
       END IF

       start(1)=istr1
       start(2)=istr2
       count(1)=nwrk1
       count(2)=nwrk2
       IF (ndims.eq.3) THEN 
        start(3)=recnum
        count(3)=1
       END IF
       status=nf_get_vara_int(fid,varid,start,count,wrk)

      status=nf_close(fid)

      deallocate(dimids,dimlens,count,start)

      END IF
      END DO

      RETURN
      END
#endif /* not MPI */

      SUBROUTINE get_2dfld(fname,varname,recnum,nwrk1,nwrk2,wrk)

      implicit none
      include 'netcdf.h'

      character(len=*) :: fname,varname
      integer :: nwrk1,nwrk2
      real(8) :: wrk(nwrk1,nwrk2)
      integer :: recnum,fid,varid,status,ndims,i,len,n1,n2
      integer,allocatable :: start(:),count(:)
      integer, allocatable :: dimids(:),dimlens(:)

      status=nf_open(TRIM(fname),nf_nowrite,fid)

       ! inquire about the dimensions of the 2D variable 
       status=nf_inq_varid(fid,TRIM(varname),varid)
       status=nf_inq_varndims(fid,varid,ndims)

       IF (ndims.gt.3) THEN
        write (*,*) 'STOP in get_2dfld, since ndims>3'
        write (*,*) 'fname: ',fname
        write (*,*) 'varname: ',varname
        STOP
       END IF

       allocate (dimids(ndims))
       allocate (dimlens(ndims))           
       allocate (start(ndims))
       allocate (count(ndims))
       status=nf_inq_vardimid(fid,varid,dimids)

       DO i=1,ndims
        status=nf_inq_dimlen(fid,dimids(i),dimlens(i))
       END DO

       n1=dimlens(1)
       n2=dimlens(2)

       IF (n1.ne.nwrk1 .or. n2.ne.nwrk2) THEN
        write (*,*) 'STOP in get_2dfld: inconsistent dims'
        STOP
       END IF

       start(1:ndims)=1
       count(1:ndims)=1
       count(1)=n1
       count(2)=n2
       IF (ndims.eq.3) THEN 
        start(3)=recnum
       END IF
       status=nf_get_vara_double(fid,varid,start,count,wrk)

      status=nf_close(fid)

      deallocate(dimids,dimlens,count,start)

      RETURN
      END

!********************************************************************
      SUBROUTINE get_2dfld_int(fname,varname,recnum,nwrk1,nwrk2,wrk)

      implicit none
      include 'netcdf.h'

      character(len=*) :: fname,varname
      integer :: nwrk1,nwrk2
      integer :: wrk(nwrk1,nwrk2)
      integer :: recnum,fid,varid,status,ndims,i,len,n1,n2
      integer,allocatable :: start(:),count(:)
      integer, allocatable :: dimids(:),dimlens(:)

      status=nf_open(TRIM(fname),nf_nowrite,fid)

       ! inquire about the dimensions of the 2D variable 
       status=nf_inq_varid(fid,TRIM(varname),varid)
       status=nf_inq_varndims(fid,varid,ndims)

       IF (ndims.gt.3) THEN
        write (*,*) 'STOP in get_2dfld, since ndims>3'
        STOP
       END IF

       allocate (dimids(ndims))
       allocate (dimlens(ndims))           
       allocate (start(ndims))
       allocate (count(ndims))
       status=nf_inq_vardimid(fid,varid,dimids)

       DO i=1,ndims
        status=nf_inq_dimlen(fid,dimids(i),dimlens(i))
       END DO

       n1=dimlens(1)
       n2=dimlens(2)

       IF (n1.ne.nwrk1 .or. n2.ne.nwrk2) THEN
        write (*,*) 'STOP in get_2dfld: inconsistent dims'
        STOP
       END IF

       start(1:ndims)=1
       count(1:ndims)=1
       count(1)=n1
       count(2)=n2
       IF (ndims.eq.3) THEN 
        start(3)=recnum
       END IF
       status=nf_get_vara_int(fid,varid,start,count,wrk)

      status=nf_close(fid)

      deallocate(dimids,dimlens,count,start)

      RETURN
      END
