      SUBROUTINE get_grddim(grdfile,
     *                      xi_rho,eta_rho,
     *                      xi_u,eta_u,
     *                      xi_v,eta_v,
     *                      xi_psi,eta_psi,
     *                      N)

      implicit none     
      include 'netcdf.h'

      character (len=80), intent (in) :: grdfile
      integer, intent (out) :: xi_rho,eta_rho,
     *                         xi_u,eta_u,
     *                         xi_v,eta_v,
     *                         xi_psi,eta_psi,
     *                         N 

      integer :: status,fid,dimid

      status=nf_open(TRIM(grdfile),nf_nowrite,fid)
! - INQUIRE DIMENSION LENTHS:
      status=nf_inq_dimid(fid,'xi_rho',dimid)
      status=nf_inq_dimlen(fid,dimid,xi_rho)

      status=nf_inq_dimid(fid,'eta_rho',dimid)
      status=nf_inq_dimlen(fid,dimid,eta_rho)

      status=nf_inq_dimid(fid,'xi_u',dimid)
      status=nf_inq_dimlen(fid,dimid,xi_u)

      status=nf_inq_dimid(fid,'xi_psi',dimid)
      status=nf_inq_dimlen(fid,dimid,xi_psi)

      status=nf_inq_dimid(fid,'eta_u',dimid)
      status=nf_inq_dimlen(fid,dimid,eta_u)

      status=nf_inq_dimid(fid,'xi_v',dimid)
      status=nf_inq_dimlen(fid,dimid,xi_v)

      status=nf_inq_dimid(fid,'eta_v',dimid)
      status=nf_inq_dimlen(fid,dimid,eta_v)
      
      status=nf_inq_dimid(fid,'eta_v',dimid)
      status=nf_inq_dimlen(fid,dimid,eta_v)
      
      status=nf_inq_dimid(fid,'eta_psi',dimid)
      status=nf_inq_dimlen(fid,dimid,eta_psi)
 
      status=nf_inq_dimid(fid,'N',dimid)
      status=nf_inq_dimlen(fid,dimid,N)

      status=nf_close(fid)

      RETURN
      END
