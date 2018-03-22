      SUBROUTINE nc_read_int_1d(fname,varname,var,n1,start,count)
    
      implicit none
      INCLUDE 'netcdf.h'

      character(len=*) :: fname
      character(len=*) :: varname
      integer :: status,varid,fid,n1
      integer :: start(*),count(*)

      real(8) :: var(n1)

      status=nf_open(TRIM(fname),nf_nowrite,fid)

      status=nf_inq_varid(fid,TRIM(varname),varid)
      status=nf_get_vara_int(fid,varid,start(1),count(1),var)

      status=nf_close(fid)

      RETURN
      END


      SUBROUTINE nc_read_double_1d(fname,varname,var,n1,start,count)

      implicit none
      INCLUDE 'netcdf.h'

      character(len=*) :: fname
      character(len=*) :: varname
      integer :: status,varid,fid,n1
      integer :: start(*),count(*)

      real(8) :: var(n1)

      status=nf_open(TRIM(fname),nf_nowrite,fid)

      status=nf_inq_varid(fid,TRIM(varname),varid)
      status=nf_get_vara_double(fid,varid,start(1),count(1),var)

      status=nf_close(fid)

      RETURN
      END


      SUBROUTINE nc_read_double_2d(fname,varname,var,n1,n2,start,count)

      implicit none
      INCLUDE 'netcdf.h'

      character(len=*) :: fname
      character(len=*) :: varname
      integer :: status,varid,fid,n1,n2
      integer :: start(*),count(*)

      real(8) :: var(n1,n2)

      status=nf_open(TRIM(fname),nf_nowrite,fid)

      status=nf_inq_varid(fid,TRIM(varname),varid)
      status=nf_get_vara_double(fid,varid,start(1:2),count(1:2),var)

      status=nf_close(fid)

      RETURN
      END

      SUBROUTINE nc_read_double_3d(fname,varname,var,n1,n2,n3,
     *                             start,count)

      implicit none
      INCLUDE 'netcdf.h'
      
      character(len=*) :: fname
      character(len=*) :: varname
      integer :: status,varid,fid,n1,n2,n3
      integer :: start(*),count(*)

      real(8) :: var(n1,n2,n3)

      status=nf_open(TRIM(fname),nf_nowrite,fid)

      status=nf_inq_varid(fid,TRIM(varname),varid)
      status=nf_get_vara_double(fid,varid,start(1:3),count(1:3),var)

      status=nf_close(fid)

      RETURN
      END


      SUBROUTINE nc_read_double_4d(fname,varname,var,n1,n2,n3,n4,
     *                             start,count)

      
      implicit none
      INCLUDe 'netcdf.h'

      character(len=*) :: fname
      character(len=*) :: varname
      integer :: status,varid,fid,n1,n2,n3,n4
      integer :: start(*),count(*)

      real(8) :: var(n1,n2,n3,n4)

      status=nf_open(TRIM(fname),nf_nowrite,fid)

      status=nf_inq_varid(fid,TRIM(varname),varid)
      status=nf_get_vara_double(fid,varid,start(1:4),count(1:4),var)

      status=nf_close(fid)

      RETURN
      END
     
