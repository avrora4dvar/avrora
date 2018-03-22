      integer function nc_dim(fname,varname)

      implicit none
      INCLUDE 'netcdf.h'
      
      character(len=*) :: fname
      character(len=*) :: varname
      integer :: status,fid,dimid

      status=nf_open(TRIM(fname),nf_nowrite,fid)
      status=nf_inq_dimid(fid,TRIM(varname),dimid)
      status=nf_inq_dimlen(fid,dimid,nc_dim) 
      status=nf_close(fid)
      

      return
      end
   
       
