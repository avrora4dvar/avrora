      SUBROUTINE ad_output_psi_T(fname,
     *                     L,M,N,
     *                     time,ad_psi,ad_T)

      USE mod_output

      implicit none
      include 'netcdf.h'
      character (len=80) :: fname
      integer :: status,ncid,recnum,varid
      integer, dimension(4) :: dimids,start,count

      integer :: L,M,N,Lp,Mp,Np
      real(8) :: time, fillValue
      real(8) :: ad_psi(0:L,0:M)
      real(8) :: ad_T(0:L,0:M,N)

      real(8), allocatable :: wrk2d(:,:)
      real(8), allocatable :: wrk3d(:,:,:)

      dimids=0

      fillValue=DBLE(1e37)

      Lp=L+1
      Mp=M+1
      Np=N+1

! Define the history file structure:
      status= nf_create(TRIM(fname),nf_clobber,ncid)

      status=nf_def_dim(ncid,'xi_rho' ,Lp, his_dimid%xi_rho )
      status=nf_def_dim(ncid,'eta_rho',Mp, his_dimid%eta_rho)
      status=nf_def_dim(ncid,'xi_u'   ,L, his_dimid%xi_u   )
      status=nf_def_dim(ncid,'eta_u'  ,Mp, his_dimid%eta_u  )
      status=nf_def_dim(ncid,'xi_v'   ,Lp, his_dimid%xi_v   )
      status=nf_def_dim(ncid,'eta_v'  ,M, his_dimid%eta_v  )
      status=nf_def_dim(ncid,'N'  ,N, his_dimid%N      )
      status=nf_def_dim(ncid,'Np'    ,Np, his_dimid%Np     )
      status=nf_def_dim(ncid,'ocean_time', nf_unlimited, his_dimid%time)

      dimids(1)=his_dimid%time

      status=nf_def_var(ncid,'ocean_time',nf_double,1,dimids,varid)
      status=nf_put_att_text(ncid,1,'units',30,
     & 'seconds since start simulation')

      dimids(1)=his_dimid%xi_rho
      dimids(2)=his_dimid%eta_rho
      dimids(3)=his_dimid%time

      status=nf_def_var(ncid,'psi',nf_double,3,dimids,varid)
      status=nf_put_att_text(ncid,2,'time',10,'ocean_time')
      status=nf_put_att_double(ncid,2,'_FillValue',
     & nf_double,1,DBLE(fillValue))
      status=nf_put_att_text(ncid,2,'units',5,'meter2 second-1')

      dimids(1)=his_dimid%xi_rho
      dimids(2)=his_dimid%eta_rho
      dimids(3)=his_dimid%N
      dimids(4)=his_dimid%time

      status=nf_def_var(ncid,'temp',nf_double,4,dimids,varid)
      status=nf_put_att_text(ncid,7,'time',10,'ocean_time')
      status=nf_put_att_double(ncid,7,'_FillValue',
     & nf_double,1,DBLE(fillValue))
      status=nf_put_att_text(ncid,7,'units',7,'Celsius')

      status=nf_close(ncid)

! Output fields:

      write (*,*) 'write to ',TRIM(fname)
      status=nf_open(TRIM(fname),nf_write,ncid)

!--------------
! output ocean_time:
!--------------

      start(1)=1
      count(1)=1
      status=nf_inq_varid(ncid,'ocean_time',varid)
      status=nf_put_vara_double(ncid,varid,start,count,time)

!--------------
! output psi:
!--------------
      start(1)=1
      start(2)=1
      start(3)=1
      count(1)=Lp
      count(2)=Mp
      count(3)=1

      allocate(wrk2d(Lp,Mp))
      wrk2d(1:Lp,1:Mp)=ad_psi(0:L,0:M)

      status=nf_inq_varid(ncid,'psi',varid)
      status=nf_put_vara_double(ncid,varid,start,count,
     *                          wrk2d)

      deallocate(wrk2d)

!--------------
! output T:
!--------------
      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=1
      count(1)=Lp
      count(2)=Mp
      count(3)=N
      count(4)=1

      allocate(wrk3d(Lp,Mp,N))
      wrk3d(1:Lp,1:Mp,1:N)=ad_T(0:L,0:M,1:N)

      status=nf_inq_varid(ncid,'temp',varid)
      status=nf_put_vara_double(ncid,varid,start,count,
     *                          wrk3d)

      deallocate(wrk3d)

      status=nf_close(ncid)

      RETURN
      END
