! AK, 11/18/2013: create a history file (to write comb there)

      SUBROUTINE create_his_file(fname,xi_rho,eta_rho,N)
 
      implicit none
      INCLUDE 'netcdf.h'

      real    :: fillValue
      integer :: xi_rho,eta_rho,
     *           xi_u,  eta_u,
     *           xi_v,  eta_v,
     *           N,Np
      integer :: ncid,di_xi_rho,di_eta_rho,
     *                di_xi_u  ,di_eta_u,
     *                di_xi_v  ,di_eta_v,
     *                di_N,di_Np,di_time,vi
      integer, dimension(4) :: dimids
      integer :: status
      character(len=*) :: fname

      dimids=0

      xi_u  =xi_rho-1
      xi_v  =xi_rho

      eta_u  =eta_rho
      eta_v  =eta_rho-1

      Np=N+1

      fillValue=HUGE(fillValue)
      fillValue=DBLE(1e37)

      status= nf_create(TRIM(fname),nf_clobber,ncid)
 
      status=nf_def_dim(ncid,'xi_rho' ,xi_rho      ,di_xi_rho)
      status=nf_def_dim(ncid,'eta_rho',eta_rho     ,di_eta_rho)
      status=nf_def_dim(ncid,'xi_u'   ,xi_u        ,di_xi_u)
      status=nf_def_dim(ncid,'eta_u'  ,eta_u       ,di_eta_u)
      status=nf_def_dim(ncid,'xi_v'   ,xi_v        ,di_xi_v)
      status=nf_def_dim(ncid,'eta_v'  ,eta_v       ,di_eta_v)
      status=nf_def_dim(ncid,'N'      ,N           ,di_N)
      status=nf_def_dim(ncid,'Np'     ,Np          ,di_Np)
      status=nf_def_dim(ncid,'ocean_time',nf_unlimited,di_time)
      

       dimids(1)=di_time

       status=nf_def_var(ncid,'ocean_time',nf_double,1,dimids,vi)
       status=nf_put_att_text(ncid,1,'units',30,
     &  'seconds since start simulation')

       dimids(1)=di_xi_rho
       dimids(2)=di_eta_rho
       dimids(3)=di_time

       status=nf_def_var(ncid,'zeta',nf_double,3,dimids,vi)
       status=nf_put_att_text(ncid,2,'time',10,'ocean_time')
!!       status=nf_put_att_double(ncid,2,'_FillValue',
!!     & nf_double,1,DBLE(fillValue))
       status=nf_put_att_text(ncid,2,'units',5,'meter')

       dimids(1)=di_xi_u
       dimids(2)=di_eta_u
       dimids(3)=di_time

       status=nf_def_var(ncid,'ubar',nf_double,3,dimids,vi)
       status=nf_put_att_text(ncid,3,'time',10,'ocean_time')
!!       status=nf_put_att_double(ncid,3,'_FillValue',
!!     & nf_double,1,DBLE(fillValue))
       status=nf_put_att_text(ncid,3,'units',14,'meter second-1')
 
       dimids(1)=di_xi_v
       dimids(2)=di_eta_v
       dimids(3)=di_time

       status=nf_def_var(ncid,'vbar',nf_double,3,dimids,vi)
       status=nf_put_att_text(ncid,4,'time',10,'ocean_time')
!!      status=nf_put_att_double(ncid,4,'_FillValue',
!!     & nf_double,1,DBLE(fillValue))
       status=nf_put_att_text(ncid,4,'units',14,'meter second-1')

       dimids(1)=di_xi_u
       dimids(2)=di_eta_u
       dimids(3)=di_N
       dimids(4)=di_time

       status=nf_def_var(ncid,'u',nf_double,4,dimids,vi)
       status=nf_put_att_text(ncid,5,'time',10,'ocean_time')
!!       status=nf_put_att_double(ncid,5,'_FillValue',
!!     & nf_double,1,DBLE(fillValue))
       status=nf_put_att_text(ncid,5,'units',14,'meter second-1')

       dimids(1)=di_xi_v
       dimids(2)=di_eta_v
       dimids(3)=di_N
       dimids(4)=di_time

       status=nf_def_var(ncid,'v',nf_double,4,dimids,vi)
       status=nf_put_att_text(ncid,6,'time',10,'ocean_time')
!!       status=nf_put_att_double(ncid,6,'_FillValue',
!!     &  nf_double,1,DBLE(fillValue))
       status=nf_put_att_text(ncid,6,'units',14,'meter second-1')

       dimids(1)=di_xi_rho
       dimids(2)=di_eta_rho
       dimids(3)=di_N
       dimids(4)=di_time

       status=nf_def_var(ncid,'temp',nf_double,4,dimids,vi)
       status=nf_put_att_text(ncid,7,'time',10,'ocean_time')
!!       status=nf_put_att_double(ncid,7,'_FillValue',
!!     & nf_double,1,DBLE(fillValue))
       status=nf_put_att_text(ncid,7,'units',7,'Celsius')

       status=nf_def_var(ncid,'salt',nf_double,4,dimids,vi)
       status=nf_put_att_text(ncid,8,'time',10,'ocean_time')
!!       status=nf_put_att_double(ncid,8,'_FillValue',
!!     & nf_double,1,DBLE(fillValue))
       status=nf_close(ncid)
 

      RETURN
      END
