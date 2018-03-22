#include "cppdefs.h"


#ifdef MPI
      SUBROUTINE get_3dfld_tile(fname,varname,recnum,
     *                    istr1,istr2,istr3,nwrk1,nwrk2,nwrk3,wrk)
      USE mod_parallel

      implicit none
      include 'netcdf.h'
      include 'mpif.h'

      character(len=*) :: fname,varname
      integer :: istr1,istr2,istr3,nwrk1,nwrk2,nwrk3
      real(8) :: wrk(nwrk1,nwrk2,nwrk3)
      integer :: recnum,fid,varid,status,ndims,i,len,n1,n2,n3
      integer,allocatable :: start(:),count(:)
      integer, allocatable :: dimids(:),dimlens(:)
      integer :: ierr

      DO irank=0,comm_size-1
      CALL mpi_barrier(comm2d,ierr)
      IF(myrank.eq.irank) THEN

      status=nf_open(TRIM(fname),nf_nowrite,fid)

       ! inquire about the dimensions of the 3D variable 
       status=nf_inq_varid(fid,TRIM(varname),varid)
       status=nf_inq_varndims(fid,varid,ndims)

       IF (ndims.gt.4) THEN
        write (*,*) 'STOP in get_3dfld, since ndims>4'
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
       n3=dimlens(3)

       IF (istr1.lt.1 .or. istr1+nwrk1-1.gt.n1 .or.
     *     istr2.lt.1 .or. istr2+nwrk2-1.gt.n2 .or.
     *     istr3.lt.1 .or. istr3+nwrk3-1.gt.n3) THEN
        write (*,*) 'STOP in get_3dfld: inconsistent dims'
        write (*,*) 'file ',TRIM(fname)
        write (*,*) 'variable ',varname

        write (*,*) 'istr1,n1,nwrk1=',istr1,n1,nwrk1
        write (*,*) 'istr2,n2,nwrk2=',istr2,n2,nwrk2
        write (*,*) 'istr3,n3,nwrk3=',istr3,n3,nwrk3

        STOP
       END IF

       start(1)=istr1
       start(2)=istr2
       start(3)=istr3
       count(1)=nwrk1
       count(2)=nwrk2
       count(3)=nwrk3
       IF (ndims.eq.4) THEN 
        start(4)=recnum
        count(4)=1
       END IF
       status=nf_get_vara_double(fid,varid,start,count,wrk)

      status=nf_close(fid)

      deallocate(dimids,dimlens,count,start)

      END IF
      END DO

      RETURN
      END

#endif /* MPI */

      SUBROUTINE get_3dfld(fname,varname,recnum,nwrk1,nwrk2,nwrk3,wrk)

      implicit none
      include 'netcdf.h'
      character(len=*) :: fname,varname
      integer :: nwrk1,nwrk2,nwrk3
      real(8) :: wrk(nwrk1,nwrk2,nwrk3)
      integer :: recnum,fid,varid,status,ndims,i,len,n1,n2,n3
      integer,allocatable :: start(:),count(:)
      integer, allocatable :: dimids(:),dimlens(:)

      status=nf_open(TRIM(fname),nf_nowrite,fid)

       ! inquire about the dimensions of the 3D variable 
       status=nf_inq_varid(fid,TRIM(varname),varid)
       status=nf_inq_varndims(fid,varid,ndims)

       IF (ndims.gt.4) THEN
        write (*,*) 'STOP in get_3dfld, since ndims>4'
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
       n3=dimlens(3)

       IF (n1.ne.nwrk1 .or. n2.ne.nwrk2  .or. n3.ne.nwrk3) THEN
        write (*,*) 'STOP in get_3dfld: inconsistent dims'
        write (*,*) 'file ',TRIM(fname)
        write (*,*) 'variable ',varname

        write (*,*) 'n1,nwrk1=',n1,nwrk1
        write (*,*) 'n2,nwrk2=',n2,nwrk2
        write (*,*) 'n3,nwrk3=',n3,nwrk3

        STOP
       END IF

       start(1:ndims)=1
       count(1:ndims)=1
       count(1)=n1
       count(2)=n2
       count(3)=n3
       IF (ndims.eq.4) THEN 
        start(4)=recnum
       END IF
       status=nf_get_vara_double(fid,varid,start,count,wrk)

      status=nf_close(fid)

      deallocate(dimids,dimlens,count,start)

      RETURN
      END
