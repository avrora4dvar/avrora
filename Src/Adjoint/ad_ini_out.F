#include "cppdefs.h"


#ifdef ADJOINT

#ifdef MPI

      SUBROUTINE ad_ini_out

!!-----------------------------------
! OUTPUT adj sensitivity to ini cond
!!-----------------------------------

      USE mod_ocean
      USE mod_output
      USE mod_parallel  

      implicit none
      include 'mpif.h'
      include 'netcdf.h'

      integer :: nini
      integer :: status,ncid,recnum,i,j,k,itr
      integer, dimension(4) :: dimids,start,count
      real(8), allocatable :: wrk2d(:,:)
      real(8), allocatable :: wrk3d(:,:,:)
      real(8) :: q0
      logical::flag_exist
   
      integer::varid_tmp
      integer::Vtransform,Vstretching
      real(8)::spar(3)
      integer,allocatable::sigma(:) 
     

! MPI related:
      integer :: i1,i2,n1,n2,ierr
      integer :: out_done
      integer,allocatable :: out_done_all(:)

      recnum=1
      nini=1

!      write (*,*) 'in ad_ini_out:...'

# ifdef TS_FIXED
      DO itr=1,Ntr
       DO k=1,N
        DO j=LBj,UBj
         DO i=LBi,UBi
! tl line 2:          tl_t(i,j,k,3,itr)=tl_t(i,j,k,nini,itr)
          ad_t(i,j,k,nini,itr)=ad_t(i,j,k,nini,itr)+ad_t(i,j,k,3,itr)
          ad_t(i,j,k,3,itr)=0.d0
! tl line 1:          tl_t(i,j,k,2,itr)=tl_t(i,j,k,nini,itr)
          ad_t(i,j,k,nini,itr)=ad_t(i,j,k,nini,itr)+ad_t(i,j,k,2,itr)
          ad_t(i,j,k,2,itr)=0.d0
         END DO
        END DO
       END DO
      END DO     
# endif

! ADJ to tl_set_depth:
      CALL ad_set_depth(LBi,UBi,LBj,UBj,N,Lm,Mm,
     *          z0_r,z0_w,oh,
     *          ad_zeta(:,:,nini),ad_z_r,ad_z_w,ad_Hz,ad_Hz_u,ad_Hz_v)
    
! Adjoint to masking:

      DO j=LBj,UBj
      DO i=LBi,UBi
       ad_zeta(i,j,nini)=ad_zeta(i,j,nini)*DBLE(mask_r(i,j))
       ad_u(i,j,1:N,nini)=ad_u(i,j,1:N,nini)*DBLE(mask_u(i,j))
       ad_v(i,j,1:N,nini)=ad_v(i,j,1:N,nini)*DBLE(mask_v(i,j))
       ad_t(i,j,1:N,nini,1:Ntr)=ad_t(i,j,1:N,nini,1:Ntr)*
     *                                       DBLE(mask_r(i,j)) 
      END DO
      END DO

      IF (myrank.eq.0) THEN
! Define the IC correction file structure:
       INQUIRE(file=TRIM(inifile),exist=flag_exist)
       write(*,*) 'file:',flag_exist

       IF(.NOT.flag_exist) THEN
       status=nf_open(TRIM(hisfile),nf_nowrite,ncid)
       status=nf_inq_varid(ncid,'Vtransform',varid_tmp)
       status=nf_get_vara_int(ncid,varid_tmp,1,1,Vtransform)       
       status=nf_inq_varid(ncid,'Vstretching',varid_tmp)
       status=nf_get_vara_int(ncid,varid_tmp,1,1,Vstretching)
       status=nf_inq_varid(ncid,'theta_s',varid_tmp)
       status=nf_get_vara_double(ncid,varid_tmp,1,1,spar(1))
       status=nf_inq_varid(ncid,'theta_b',varid_tmp)
       status=nf_get_vara_double(ncid,varid_tmp,1,1,spar(2))
       status=nf_inq_varid(ncid,'Tcline',varid_tmp)
       status=nf_get_vara_double(ncid,varid_tmp,1,1,spar(3)) 
       status=nf_close(ncid)

       status= nf_create(TRIM(inifile),nf_clobber,ncid)
        
       status=nf_def_dim(ncid,'xi_rho' ,Lp          ,dimid_xi_rho )
       status=nf_def_dim(ncid,'eta_rho',Mp          ,dimid_eta_rho)
       status=nf_def_dim(ncid,'xi_u'   ,L           ,dimid_xi_u   )
       status=nf_def_dim(ncid,'eta_u'  ,Mp          ,dimid_eta_u  )
       status=nf_def_dim(ncid,'xi_v'   ,Lp          ,dimid_xi_v   )
       status=nf_def_dim(ncid,'eta_v'  ,M           ,dimid_eta_v  )
       status=nf_def_dim(ncid,'N'      ,N           ,dimid_N      )
       status=nf_def_dim(ncid,'Np'     ,Np          ,dimid_Np     )
       status=nf_def_dim(ncid,'ocean_time' ,nf_unlimited,dimid_time   )

       status=nf_def_var(ncid,'Vtransform',nf_int,0,0,varid_tmp)
       status=nf_def_var(ncid,'Vstretching',nf_int,0,0,varid_tmp)
       status=nf_def_var(ncid,'theta_s',nf_double,0,0,varid_tmp)
       status=nf_def_var(ncid,'theta_b',nf_double,0,0,varid_tmp)
       status=nf_def_var(ncid,'Tcline',nf_double,0,0,varid_tmp)
       status=nf_def_var(ncid,'s_rho',nf_double,1,dimid_N,varid_tmp)
       status=nf_def_var(ncid,'s_w',nf_double,1,dimid_N,varid_tmp)

       dimids(1)=dimid_time

       status=nf_def_var(ncid,'ocean_time',nf_double,1,dimids,
     *               varid_ocean_time)

       dimids(1)=dimid_xi_rho
       dimids(2)=dimid_eta_rho
       dimids(3)=dimid_time

       status=nf_def_var(ncid,'zeta',nf_double,3,dimids,
     *               varid_zeta)
       status=nf_put_att_text(ncid,varid_zeta,'time',
     * 10,'ocean_time')

       dimids(1)=dimid_xi_u
       dimids(2)=dimid_eta_u
       dimids(3)=dimid_time

       status=nf_def_var(ncid,'ubar',nf_double,3,dimids,
     *               varid_ubar)
       status=nf_put_att_text(ncid,varid_ubar,'time',
     * 10,'ocean_time')

       dimids(1)=dimid_xi_v
       dimids(2)=dimid_eta_v
       dimids(3)=dimid_time

       status=nf_def_var(ncid,'vbar',nf_double,3,dimids,
     *               varid_vbar)
       status=nf_put_att_text(ncid,varid_vbar,'time',
     * 10,'ocean_time')

       dimids(1)=dimid_xi_u
       dimids(2)=dimid_eta_u
       dimids(3)=dimid_N
       dimids(4)=dimid_time

       status=nf_def_var(ncid,'u',nf_double,4,dimids,
     *               varid_u)
       status=nf_put_att_text(ncid,varid_u,'time',
     * 10,'ocean_time')

       dimids(1)=dimid_xi_v
       dimids(2)=dimid_eta_v
       dimids(3)=dimid_N
       dimids(4)=dimid_time

       status=nf_def_var(ncid,'v',nf_double,4,dimids,
     *               varid_v)
       status=nf_put_att_text(ncid,varid_v,'time',
     * 10,'ocean_time')

       dimids(1)=dimid_xi_rho
       dimids(2)=dimid_eta_rho
       dimids(3)=dimid_N
       dimids(4)=dimid_time

       status=nf_def_var(ncid,'temp',nf_double,4,dimids,
     *               varid_temp)
       status=nf_put_att_text(ncid,varid_temp,'time',
     * 10,'ocean_time')

       status=nf_def_var(ncid,'salt',nf_double,4,dimids,
     *               varid_salt)
       status=nf_put_att_text(ncid,varid_salt,'time',
     * 10,'ocean_time')

       status=nf_close(ncid)
      ELSE
       status=nf_open(TRIM(inifile),nf_nowrite,ncid)
       status=nf_inq_varid(ncid,'salt',varid_salt)
       status=nf_inq_varid(ncid,'temp',varid_temp)
       status=nf_inq_varid(ncid,'u',varid_u)
       status=nf_inq_varid(ncid,'v',varid_v)
       status=nf_inq_varid(ncid,'zeta',varid_zeta)
       status=nf_inq_varid(ncid,'ubar',varid_ubar)
       status=nf_inq_varid(ncid,'vbar',varid_vbar)
       status=nf_inq_varid(ncid,'ocean_time',varid_ocean_time)
       status=nf_close(ncid)
      END IF
      END IF ! <- if (myrank.eq.0) 
 
      CALL mpi_bcast(varid_ocean_time,1,mpi_integer,0,mpi_comm_world,
     *                ierr)
      CALL mpi_bcast(varid_zeta,1,mpi_integer,0,mpi_comm_world,ierr)
      CALL mpi_bcast(varid_ubar,1,mpi_integer,0,mpi_comm_world,ierr)
      CALL mpi_bcast(varid_vbar,1,mpi_integer,0,mpi_comm_world,ierr)
      CALL mpi_bcast(varid_u,1,mpi_integer,0,mpi_comm_world,ierr)
      CALL mpi_bcast(varid_v,1,mpi_integer,0,mpi_comm_world,ierr)
      CALL mpi_bcast(varid_temp,1,mpi_integer,0,mpi_comm_world,ierr)
      CALL mpi_bcast(varid_salt,1,mpi_integer,0,mpi_comm_world,ierr)

!------------
! OUTPUT:
! - Output arrays from each processor to the single netcdf file, 
! - one at a time
!------------
      start=0
      count=0


!--------------------------------------------------
! output ocean_time on the master thread (myrank=0)
!--------------------------------------------------
      IF (myrank.eq.0) THEN
       status=nf_open(TRIM(inifile),nf_write,ncid)

!Output grid layout 
       status=nf_inq_varid(ncid,'Vtransform',varid_tmp)
       status=nf_put_vara_int(ncid,varid_tmp,1,1,Vtransform)       
       status=nf_inq_varid(ncid,'Vstretching',varid_tmp)
       status=nf_put_vara_int(ncid,varid_tmp,1,1,Vstretching)
       status=nf_inq_varid(ncid,'theta_s',varid_tmp)
       status=nf_put_vara_double(ncid,varid_tmp,1,1,spar(1))
       status=nf_inq_varid(ncid,'theta_b',varid_tmp)
       status=nf_put_vara_double(ncid,varid_tmp,1,1,spar(2))
       status=nf_inq_varid(ncid,'Tcline',varid_tmp)
       status=nf_put_vara_double(ncid,varid_tmp,1,1,spar(3))

       ALLOCATE(sigma(size(ad_u,3)))
       DO i1=1,size(sigma)
         sigma(i1)=-dble(1)+(dble(i1)-.5)/dble(size(sigma))
       END DO
       status=nf_inq_varid(ncid,'s_rho',varid_tmp)
       status=nf_put_vara_double(ncid,varid_tmp,1,size(sigma),sigma)
       DEALLOCATE(sigma)

       ALLOCATE(sigma(size(ad_u,3)+1))
       DO i1=1,size(sigma)
         sigma(i1)=-dble(1)+(dble(i1)-1.0)/dble(size(sigma)-1)
       END DO
       status=nf_inq_varid(ncid,'s_w',varid_tmp)
       status=nf_put_vara_double(ncid,varid_tmp,1,size(sigma),sigma)
       DEALLOCATE(sigma)

! output ocean_time:

       start(1)=recnum
       count(1)=1
       status=nf_put_vara_double(ncid,varid_ocean_time,start,count,
     *                                                time)
       status=nf_close(ncid)
      END IF  ! (myrank.eq.0)

!--------------------------------------------------------------------
! out_done_all(0:MaxTile) is initialized on each tile
! (note out_done_all will gather out_done values from each processor
!--------------------------------------------------------------------

      allocate(out_done_all(0:MaxTile))
      DO k=0,MaxTile
       out_done_all(k)=0
      END DO

      out_done=0  ! -> to be updated to 1 after output's done for the tile

!------------------------------------------------------------------
! - output 1 tile at a time:
!------------------------------------------------------------------
      DO k=0,MaxTile
 !- gather out_done values for each tile and and broadcast them 
       CALL mpi_gather(out_done,1,mpi_integer,
     *                  out_done_all,1,mpi_integer,
     *                  0,mpi_comm_world,ierr)        
             
       CALL mpi_bcast(out_done_all,MaxTile+1,mpi_integer,
     *                 0,mpi_comm_world,ierr)

       IF (myrank.eq.k .and. out_done_all(k).eq.0) THEN

        status=nf_open(TRIM(inifile),nf_write,ncid)

!--------------
! output zeta:
!--------------
! - output excludes ghost points (but includes boundary points)
        n1=Iend_E-Istr_E+1
        n2=Jend_E-Jstr_E+1

        allocate(wrk2d(n1,n2))
        wrk2d(1:n1,1:n2)=ad_zeta(Istr_E:Iend_E,Jstr_E:Jend_E,nini)
        start(1)=Istr_E+1
        start(2)=Jstr_E+1
        start(3)=recnum
        count(1)=n1
        count(2)=n2
        count(3)=1

        status=nf_put_vara_double(ncid,varid_zeta,start,count,
     *                                          wrk2d)
        deallocate(wrk2d)
       
!--------------
! output ubar (dummy: ==0)
!--------------
! - new: output excludes ghost points (but includes boundary points)
        n1=IendU_E-IstrU_E+1
        n2=JendU_E-JstrU_E+1

        allocate(wrk2d(n1,n2))
        wrk2d(1:n1,1:n2)=0.d0
        start(1)=IstrU_E
        start(2)=JstrU_E+1
        start(3)=recnum
        count(1)=n1
        count(2)=n2
        count(3)=1

        status=nf_put_vara_double(ncid,varid_ubar,start,count,wrk2d)
        deallocate(wrk2d)
       
!--------------
! output vbar (dummy: == 0)
!--------------
! - new: output excludes ghost points (but includes boundary points)
        n1=IendV_E-IstrV_E+1
        n2=JendV_E-JstrV_E+1

        allocate(wrk2d(n1,n2))
        wrk2d(1:n1,1:n2)=0.d0 
        start(1)=IstrV_E+1
        start(2)=JstrV_E
        start(3)=recnum
        count(1)=n1
        count(2)=n2
        count(3)=1

        status=nf_put_vara_double(ncid,varid_vbar,start,count,wrk2d)
        deallocate(wrk2d)
       
!--------------
! output u:
!--------------
        n1=IendU_E-IstrU_E+1
        n2=JendU_E-JstrU_E+1

        allocate(wrk3d(n1,n2,N))
        wrk3d(1:n1,1:n2,1:N)=
     *        ad_u(IstrU_E:IendU_E,JstrU_E:JendU_E,1:N,nini)
        start(1)=IstrU_E
        start(2)=JstrU_E+1
        start(3)=1
        start(4)=recnum
        count(1)=n1
        count(2)=n2
        count(3)=N
        count(4)=1

        status=nf_put_vara_double(ncid,varid_u,start,count,wrk3d)
        deallocate(wrk3d)
       
!--------------
! output v:
!--------------
        n1=IendV_E-IstrV_E+1
        n2=JendV_E-JstrV_E+1

        allocate(wrk3d(n1,n2,N))
        wrk3d(1:n1,1:n2,1:N)=
     *     ad_v(IstrV_E:IendV_E,JstrV_E:JendV_E,1:N,nini)
        start(1)=IstrV_E+1
        start(2)=JstrV_E
        start(3)=1
        start(4)=recnum
        count(1)=n1
        count(2)=n2
        count(3)=N
        count(4)=1

        status=nf_put_vara_double(ncid,varid_v,start,count,wrk3d)
        deallocate(wrk3d)
       
!--------------
! output T, S:
!--------------
        n1=Iend_E-Istr_E+1
        n2=Jend_E-Jstr_E+1

        allocate(wrk3d(n1,n2,N))
        start(1)=Istr_E+1
        start(2)=Jstr_E+1
        start(3)=1
        start(4)=recnum
        count(1)=n1
        count(2)=n2
        count(3)=N
        count(4)=1

        wrk3d(1:n1,1:n2,1:N)=
     *     ad_t(Istr_E:Iend_E,Jstr_E:Jend_E,1:N,nini,1)
        status=nf_put_vara_double(ncid,varid_temp,start,count,wrk3d)

        wrk3d(1:n1,1:n2,1:N)=
     *     ad_t(Istr_E:Iend_E,Jstr_E:Jend_E,1:N,nini,2)
        status=nf_put_vara_double(ncid,varid_salt,start,count,wrk3d)

        deallocate(wrk3d)

        status=nf_close(ncid)

        out_done=1

       END IF ! (myrank.eq.k .and. out_done_all(k).eq.0)
      END DO ! k=0:MaxTileNo     

      deallocate(out_done_all)

      RETURN
      END

! END OF THE MPI PIECE

!==========================================================================
!==========================================================================
# else /* serial */

      SUBROUTINE ad_ini_out

!!-----------------------------------
! OUTPUT adj sensitivity to ini cond
!!-----------------------------------

      USE mod_ocean
      USE mod_nf_param
      USE mod_output

      implicit none
      integer :: nini
      integer :: status,ncid,recnum,i,j,k,itr
      integer, dimension(4) :: dimids,start,count
      real(8), allocatable :: wrk2d(:,:)
      real(8), allocatable :: wrk3d(:,:,:)
      real(8) :: q0

      recnum=1
      nini=1

      write (*,*) 'in ad_ini_out:...'

# ifdef TS_FIXED
      DO itr=1,Ntr
       DO k=1,N
        DO j=LBj,UBj
         DO i=LBi,UBi
! tl line 2:          tl_t(i,j,k,3,itr)=tl_t(i,j,k,nini,itr)
          ad_t(i,j,k,nini,itr)=ad_t(i,j,k,nini,itr)+ad_t(i,j,k,3,itr)
          ad_t(i,j,k,3,itr)=0.d0
! tl line 1:          tl_t(i,j,k,2,itr)=tl_t(i,j,k,nini,itr)
          ad_t(i,j,k,nini,itr)=ad_t(i,j,k,nini,itr)+ad_t(i,j,k,2,itr)
          ad_t(i,j,k,2,itr)=0.d0
         END DO
        END DO
       END DO
      END DO     
# endif

! ADJ to tl_set_depth:
      CALL ad_set_depth(LBi,UBi,LBj,UBj,N,Lm,Mm,
     *          z0_r,z0_w,oh,
     *          ad_zeta(:,:,nini),ad_z_r,ad_z_w,ad_Hz,ad_Hz_u,ad_Hz_v)
      
#ifdef NS_PERIODIC

      CALL ad_exchange_2d(ad_zeta(:,:,nini),LBi,UBi,LBj,UBj,Lm,Mm)
      CALL ad_exchange_3d(ad_u(:,:,:,nini),LBi,UBi,LBj,UBj,1,N,Lm,Mm)
      CALL ad_exchange_3d(ad_v(:,:,:,nini),LBi,UBi,LBj,UBj,1,N,Lm,Mm)
      CALL ad_exchange_3d(ad_t(:,:,:,nini,1),LBi,UBi,LBj,UBj,1,N,Lm,Mm)
      CALL ad_exchange_3d(ad_t(:,:,:,nini,2),LBi,UBi,LBj,UBj,1,N,Lm,Mm)

#endif /* NS_PERIODIC */

! Adjoint to masking:

      DO j=LBj,UBj
      DO i=LBi,UBi
       ad_zeta(i,j,nini)=ad_zeta(i,j,nini)*DBLE(mask_r(i,j))
       ad_u(i,j,1:N,nini)=ad_u(i,j,1:N,nini)*DBLE(mask_u(i,j))
       ad_v(i,j,1:N,nini)=ad_v(i,j,1:N,nini)*DBLE(mask_v(i,j))
       ad_t(i,j,1:N,nini,1:Ntr)=ad_t(i,j,1:N,nini,1:Ntr)*
     *                                       DBLE(mask_r(i,j)) 
      END DO
      END DO

! Define the IC correction file structure:       
       status= nf_create(TRIM(inifile),nf_clobber,ncid)
        
       status=nf_def_dim(ncid,'xi_rho' ,Lp          ,his_dimid%xi_rho )
       status=nf_def_dim(ncid,'eta_rho',Mp          ,his_dimid%eta_rho)
       status=nf_def_dim(ncid,'xi_u'   ,L           ,his_dimid%xi_u   )
       status=nf_def_dim(ncid,'eta_u'  ,Mp          ,his_dimid%eta_u  )
       status=nf_def_dim(ncid,'xi_v'   ,Lp          ,his_dimid%xi_v   )
       status=nf_def_dim(ncid,'eta_v'  ,M           ,his_dimid%eta_v  )
       status=nf_def_dim(ncid,'N'      ,N           ,his_dimid%N      )
       status=nf_def_dim(ncid,'Np'     ,Np          ,his_dimid%Np     )
       status=nf_def_dim(ncid,'ocean_time',nf_unlimited,his_dimid%time)

       dimids(1)=his_dimid%time

       status=nf_def_var(ncid,'ocean_time',nf_double,1,dimids,
     *               his_varid%ocean_time)

       dimids(1)=his_dimid%xi_rho
       dimids(2)=his_dimid%eta_rho
       dimids(3)=his_dimid%time

       status=nf_def_var(ncid,'zeta',nf_double,3,dimids,
     *               his_varid%zeta)
       status=nf_put_att_text(ncid,his_varid%zeta,'time',
     * 10,'ocean_time')


       dimids(1)=his_dimid%xi_u
       dimids(2)=his_dimid%eta_u
       dimids(3)=his_dimid%time

       status=nf_def_var(ncid,'ubar',nf_double,3,dimids,
     *               his_varid%ubar)
      status=nf_put_att_text(ncid,his_varid%ubar,'time',
     * 10,'ocean_time')

       dimids(1)=his_dimid%xi_v
       dimids(2)=his_dimid%eta_v
       dimids(3)=his_dimid%time

       status=nf_def_var(ncid,'vbar',nf_double,3,dimids,
     *               his_varid%vbar)
      status=nf_put_att_text(ncid,his_varid%vbar,'time',
     * 10,'ocean_time')


       dimids(1)=his_dimid%xi_u
       dimids(2)=his_dimid%eta_u
       dimids(3)=his_dimid%N
       dimids(4)=his_dimid%time

       status=nf_def_var(ncid,'u',nf_double,4,dimids,
     *               his_varid%u)
      status=nf_put_att_text(ncid,his_varid%u,'time',
     * 10,'ocean_time')


       dimids(1)=his_dimid%xi_v
       dimids(2)=his_dimid%eta_v
       dimids(3)=his_dimid%N
       dimids(4)=his_dimid%time

       status=nf_def_var(ncid,'v',nf_double,4,dimids,
     *               his_varid%v)
      status=nf_put_att_text(ncid,his_varid%v,'time',
     * 10,'ocean_time')


       dimids(1)=his_dimid%xi_rho
       dimids(2)=his_dimid%eta_rho
       dimids(3)=his_dimid%N
       dimids(4)=his_dimid%time

       status=nf_def_var(ncid,'temp',nf_double,4,dimids,
     *               his_varid%temp)
      status=nf_put_att_text(ncid,his_varid%temp,'time',
     * 10,'ocean_time')


       status=nf_def_var(ncid,'salt',nf_double,4,dimids,
     *               his_varid%salt)
      status=nf_put_att_text(ncid,his_varid%salt,'time',
     * 10,'ocean_time')


       status=nf_close(ncid)

!------------
! OUTPUT:
!------------
       status=nf_open(TRIM(inifile),nf_write,ncid)

! output ocean_time:

       start(1)=recnum
       count(1)=1
       status=nf_put_vara_double(ncid,his_varid%ocean_time,start,count,
     *                                                time)

!--------------
! output zeta:
!--------------
       allocate(wrk2d(Lp,Mp))
       wrk2d(1:Lp,1:Mp)=ad_zeta(0:L,0:M,nini)
       start(1)=1
       start(2)=1
       start(3)=recnum
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
       wrk2d(1:L,1:Mp)=0.d0
       start(1)=1
       start(2)=1
       start(3)=recnum
       count(1)=L
       count(2)=Mp
       count(3)=1
       status=nf_put_vara_double(ncid,his_varid%ubar,start,count,wrk2d)
       deallocate(wrk2d)
       
!--------------
! output vbar:
!--------------
       allocate(wrk2d(Lp,M))
       wrk2d(1:Lp,1:M)=0.d0
       start(1)=1
       start(2)=1
       start(3)=recnum
       count(1)=Lp
       count(2)=M
       count(3)=1
       status=nf_put_vara_double(ncid,his_varid%vbar,start,count,wrk2d)
       deallocate(wrk2d)
       
!--------------
! output u:
!--------------
       allocate(wrk3d(L,Mp,N))
       wrk3d(1:L,1:Mp,1:N)=ad_u(1:L,0:M,1:N,nini)
       start(1)=1
       start(2)=1
       start(3)=1
       start(4)=recnum
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
       wrk3d(1:Lp,1:M,1:N)=ad_v(0:L,1:M,1:N,nini)
       start(1)=1
       start(2)=1
       start(3)=1
       start(4)=recnum
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
       start(4)=recnum
       count(1)=Lp
       count(2)=Mp
       count(3)=N
       count(4)=1

       wrk3d(1:Lp,1:Mp,1:N)=ad_t(0:L,0:M,1:N,nini,1)
       status=nf_put_vara_double(ncid,his_varid%temp,start,count,wrk3d)

       wrk3d(1:Lp,1:Mp,1:N)=ad_t(0:L,0:M,1:N,nini,2)
       status=nf_put_vara_double(ncid,his_varid%salt,start,count,wrk3d)

       deallocate(wrk3d)

       status=nf_close(ncid)

      RETURN
      END

# endif /* MPI or serial */

#endif /* ADJOINT */
