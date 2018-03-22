! AK, 11/18/2013: fill-in the newly created history file (see create_his_file)
! with zeros and the given ocean_time:

      SUBROUTINE zero_his_file(fname,xi_rho,eta_rho,N,nt,ocean_time)


      implicit none
      INCLUDE 'netcdf.h'

      character(len=80) :: fname
      integer :: xi_rho,eta_rho,xi_u,eta_u,xi_v,eta_v,N,nt
      integer :: start(4),count(4),fid,varid,status,it
      real(8) :: ocean_time(nt)

      real(8), allocatable :: tmp3d(:,:,:)

      xi_u  =xi_rho-1
      xi_v  =xi_rho

      eta_u  =eta_rho
      eta_v  =eta_rho-1

      allocate(tmp3d(xi_rho,eta_rho,N))
      tmp3d=0.d0
       
      status=nf_open(TRIM(fname),nf_write,fid)

      DO it=1,nt
       status=nf_inq_varid(fid,'ocean_time',varid)
       status=nf_put_vara_double(fid,varid,it,1,ocean_time(it))

       status=nf_inq_varid(fid,'zeta',varid)
       start(1)=1
       start(2)=1
       start(3)=it
       count(1)=xi_rho
       count(2)=eta_rho
       count(3)=1
       status=nf_put_vara_double(fid,varid,start,count,
     *                            tmp3d(1:xi_rho,1:eta_rho,1))

       status=nf_inq_varid(fid,'ubar',varid)
       start(1)=1
       start(2)=1
       start(3)=it
       count(1)=xi_u
       count(2)=eta_u
       count(3)=1
       status=nf_put_vara_double(fid,varid,start,count,
     *                            tmp3d(1:xi_u,1:eta_u,1))

       status=nf_inq_varid(fid,'vbar',varid)
       start(1)=1
       start(2)=1
       start(3)=it
       count(1)=xi_v
       count(2)=eta_v
       count(3)=1
       status=nf_put_vara_double(fid,varid,start,count,
     *                            tmp3d(1:xi_v,1:eta_v,1))

       status=nf_inq_varid(fid,'u',varid)
       start(1)=1
       start(2)=1
       start(3)=1
       start(4)=it
       count(1)=xi_u
       count(2)=eta_u
       count(3)=N
       count(4)=1
       status=nf_put_vara_double(fid,varid,start,count,
     *                            tmp3d(1:xi_u,1:eta_u,1:N))

       status=nf_inq_varid(fid,'v',varid)
       start(1)=1
       start(2)=1
       start(3)=1
       start(4)=it
       count(1)=xi_v
       count(2)=eta_v
       count(3)=N
       count(4)=1
       status=nf_put_vara_double(fid,varid,start,count,
     *                           tmp3d(1:xi_v,1:eta_v,1:N))
                              
       status=nf_inq_varid(fid,'temp',varid)
       start(1)=1
       start(2)=1
       start(3)=1
       start(4)=it
       count(1)=xi_rho
       count(2)=eta_rho
       count(3)=N
       count(4)=1
       status=nf_put_vara_double(fid,varid,start,count,tmp3d)

       status=nf_inq_varid(fid,'salt',varid)
       status=nf_put_vara_double(fid,varid,start,count,tmp3d)

      END DO
      status=nf_close(fid)
      
      deallocate(tmp3d)

      RETURN
      END
       
