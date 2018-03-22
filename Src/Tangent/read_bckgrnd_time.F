#include "cppdefs.h"

#if defined TANGENT || defined ADJOINT
      SUBROUTINE read_fwd_time

! A. Kurapov, 5/7/2007: 
! this subroutine 
! checks the length of the background time series, 
! allocates the background time vector, 
! reads background times (in fwd_time)

 
      USE mod_ocean
      implicit none
      include 'netcdf.h'

      status=nf_open(TRIM(fwdfile),nf_nowrite,fid)
      status=nf_inq_varid(fid,'ocean_time',varid)
      status=nf_inq_vardimid(fid,varid,dimid)
      status=nf_inq_dimlen(fid,dimid,nfwd)

      allocate(fwd_time(nfwd))
      status=nf_get_vara_double(fid,varid,1,nfwd,fwd_time)     

      status=nf_close(fid)
      
      RETURN
      END
 
#endif /* defined TANGENT || defined ADJOINT */
