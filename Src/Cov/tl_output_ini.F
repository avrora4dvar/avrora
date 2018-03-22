      SUBROUTINE tl_output_ini(outfile,
     *                   L,M,N,
     *                   time,tl_zeta,tl_ubar,tl_vbar,
     *                   tl_u,tl_v,tl_T,tl_S)

      USE mod_output

      implicit none
      include 'netcdf.h'

      character (len=80) :: outfile
      integer :: status,ncid,recnum
      integer, dimension(4) :: dimids,start,count

      integer :: L,M,N,Lp,Mp,Np

      real(8) :: time, fillValue
      real(8) :: tl_zeta(0:L,0:M)
      real(8) :: tl_ubar(0:L,0:M)
      real(8) :: tl_vbar(0:L,0:M)
      real(8) :: tl_u(0:L,0:M,N)
      real(8) :: tl_v(0:L,0:M,N)
      real(8) :: tl_T(0:L,0:M,N)
      real(8) :: tl_S(0:L,0:M,N)
      
      real(8), allocatable :: wrk2d(:,:)
      real(8), allocatable :: wrk3d(:,:,:)

      dimids=0

      fillValue=DBLE(REAL(1e37))

      Lp=L+1
      Mp=M+1
      Np=N+1

! Define the history file structure:
      status= nf_create(TRIM(outfile),nf_clobber,ncid)

      status=nf_def_dim(ncid,'xi_rho' ,Lp          , his_dimid%xi_rho )
      status=nf_def_dim(ncid,'eta_rho',Mp          ,his_dimid%eta_rho)
      status=nf_def_dim(ncid,'xi_u'   ,L           ,his_dimid%xi_u   )
      status=nf_def_dim(ncid,'eta_u'  ,Mp          ,his_dimid%eta_u  )
      status=nf_def_dim(ncid,'xi_v'   ,Lp          ,his_dimid%xi_v   )
      status=nf_def_dim(ncid,'eta_v'  ,M           ,his_dimid%eta_v  )
      status=nf_def_dim(ncid,'N'  ,N           ,his_dimid%N      )
      status=nf_def_dim(ncid,'Np'    ,Np          ,his_dimid%Np     )
      status=nf_def_dim(ncid,'ocean_time',nf_unlimited,his_dimid%time)

      dimids(1)=his_dimid%time

      status=nf_def_var(ncid,'ocean_time',nf_double,1,dimids,
     *               his_varid%ocean_time)
      status=nf_put_att_text(ncid,1,'units',30,
     &  'seconds since start simulation')

      dimids(1)=his_dimid%xi_rho
      dimids(2)=his_dimid%eta_rho
      dimids(3)=his_dimid%time

      status=nf_def_var(ncid,'zeta',nf_double,3,dimids,
     *               his_varid%zeta)
      status=nf_put_att_text(ncid,his_varid%zeta,'time',10,'ocean_time')
       status=nf_put_att_double(ncid,2,'_FillValue',
     & nf_double,1,DBLE(fillValue))
       status=nf_put_att_text(ncid,2,'units',5,'meter')

      dimids(1)=his_dimid%xi_u
      dimids(2)=his_dimid%eta_u
      dimids(3)=his_dimid%time

      status=nf_def_var(ncid,'ubar',nf_double,3,dimids,
     *               his_varid%ubar)
      status=nf_put_att_text(ncid,his_varid%ubar,'time',10,'ocean_time')
      status=nf_put_att_double(ncid,3,'_FillValue',
     & nf_double,1,DBLE(fillValue))
       status=nf_put_att_text(ncid,3,'units',14,'meter second-1')

      dimids(1)=his_dimid%xi_v
      dimids(2)=his_dimid%eta_v 
      dimids(3)=his_dimid%time

      status=nf_def_var(ncid,'vbar',nf_double,3,dimids,
     *               his_varid%vbar)
      status=nf_put_att_text(ncid,his_varid%vbar,'time',10,'ocean_time')
      status=nf_put_att_double(ncid,4,'_FillValue',
     & nf_double,1,DBLE(fillValue))
       status=nf_put_att_text(ncid,4,'units',14,'meter second-1')


      dimids(1)=his_dimid%xi_u
      dimids(2)=his_dimid%eta_u
      dimids(3)=his_dimid%N
      dimids(4)=his_dimid%time


      status=nf_def_var(ncid,'u',nf_double,4,dimids,
     *               his_varid%u)
      status=nf_put_att_text(ncid,his_varid%u,'time',10,'ocean_time')

      dimids(1)=his_dimid%xi_v
      dimids(2)=his_dimid%eta_v
      dimids(3)=his_dimid%N
      dimids(4)=his_dimid%time

      status=nf_def_var(ncid,'v',nf_double,4,dimids,
     *               his_varid%v)
      status=nf_put_att_text(ncid,his_varid%v,'time',10,'ocean_time')
      status=nf_put_att_double(ncid,5,'_FillValue',
     & nf_double,1,DBLE(fillValue))
      status=nf_put_att_text(ncid,5,'units',14,'meter second-1')

      dimids(1)=his_dimid%xi_rho
      dimids(2)=his_dimid%eta_rho
      dimids(3)=his_dimid%N
      dimids(4)=his_dimid%time

      status=nf_def_var(ncid,'temp',nf_double,4,dimids,
     *               his_varid%temp)
      status=nf_put_att_text(ncid,his_varid%temp,'time',10,'ocean_time')
      status=nf_put_att_double(ncid,7,'_FillValue',
     & nf_double,1,DBLE(fillValue))
      status=nf_put_att_text(ncid,7,'units',7,'Celsius')

      status=nf_def_var(ncid,'salt',nf_double,4,dimids,
     *               his_varid%salt)
      status=nf_put_att_text(ncid,his_varid%salt,'time',10,'ocean_time')
      status=nf_put_att_double(ncid,8,'_FillValue',
     & nf_double,1,DBLE(fillValue))

      status=nf_close(ncid)

      status=nf_open(TRIM(outfile),nf_write,ncid)

!--------------
! output ocean_time:
!--------------
      start(1)=1
      count(1)=1
      status=nf_put_vara_double(ncid,his_varid%ocean_time,start,count,
     *                                                time)

!--------------
! output zeta:
!--------------
      allocate(wrk2d(Lp,Mp))
      wrk2d(1:Lp,1:Mp)=tl_zeta(0:L,0:M)
      start(1)=1
      start(2)=1
      start(3)=1
      count(1)=Lp
      count(2)=Mp
      count(3)=1
      status=nf_put_vara_double(ncid,his_varid%zeta,start,count,
     *                                          wrk2d)
      deallocate(wrk2d)
      
      
!--------------
! output ubar:
!--------------
      allocate(wrk2d(L,Mp))
      wrk2d(1:L,1:Mp)=tl_ubar(1:L,0:M)
      start(1)=1
      start(2)=1
      start(3)=1
      count(1)=L
      count(2)=Mp
      count(3)=1
      status=nf_put_vara_double(ncid,his_varid%ubar,start,count,wrk2d)
      deallocate(wrk2d)

!--------------
! output vbar:
!--------------
      allocate(wrk2d(Lp,M))
      wrk2d(1:Lp,1:M)=tl_vbar(0:L,1:M)
      start(1)=1
      start(2)=1
      start(3)=1
      count(1)=Lp
      count(2)=M
      count(3)=1
      status=nf_put_vara_double(ncid,his_varid%vbar,start,count,wrk2d)
      deallocate(wrk2d)

!--------------
! output u:
!--------------
      allocate(wrk3d(L,Mp,N))
      wrk3d(1:L,1:Mp,1:N)=tl_u(1:L,0:M,1:N)
      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=1
      count(1)=L
      count(2)=Mp
      count(3)=N
      count(4)=1

      status=nf_put_vara_double(ncid,his_varid%u,start,count,wrk3d)
      deallocate(wrk3d)

!--------------
! output v:
!--------------
      allocate(wrk3d(Lp,M,N))
      wrk3d(1:Lp,1:M,1:N)=tl_v(0:L,1:M,1:N)
      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=1
      count(1)=Lp
      count(2)=M
      count(3)=N
      count(4)=1
      status=nf_put_vara_double(ncid,his_varid%v,start,count,wrk3d)
      deallocate(wrk3d)

!--------------
! output T, S:
!--------------
      allocate(wrk3d(Lp,Mp,N))
      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=1
      count(1)=Lp
      count(2)=Mp
      count(3)=N
      count(4)=1

      wrk3d(1:Lp,1:Mp,1:N)=tl_T(0:L,0:M,1:N)
      status=nf_put_vara_double(ncid,his_varid%temp,start,count,wrk3d)

      wrk3d(1:Lp,1:Mp,1:N)=tl_S(0:L,0:M,1:N)
      status=nf_put_vara_double(ncid,his_varid%salt,start,count,wrk3d)

      deallocate(wrk3d)

      status=nf_close(ncid)

      RETURN
      END
