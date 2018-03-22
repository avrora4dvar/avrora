#include "cppdefs.h"

#ifdef TANGENT

#ifdef MPI
      SUBROUTINE tl_output
! Note: output 2D variables (zeta, ubar, vbar) w/ index 1 (such as zeta(:,:,1))
! output 3D variables  w/ index nstp (such as zeta(:,:,nstp))
      USE mod_ocean
      USE mod_output
      USE mod_parallel  

      implicit none
      include 'mpif.h'
      include 'netcdf.h'
      integer :: status,ncid,recnum
      integer, dimension(4) :: dimids,start,count
      real(8), allocatable :: wrk2d(:,:)
      real(8), allocatable :: wrk3d(:,:,:)
      real(8) :: q0

! MPI related:
      integer :: i1,i2,n1,n2,ierr,k
      integer :: out_done
      integer,allocatable :: out_done_all(:)

 
      IF (iic.eq.1) THEN

       dimids=0

       IF (myrank.eq.0) THEN

! Define the history file structure:       
        status= nf_create(TRIM(hisfile),nf_clobber,ncid)
        
        status=nf_def_dim(ncid,'xi_rho' ,Lp          ,dimid_xi_rho )
        status=nf_def_dim(ncid,'eta_rho',Mp          ,dimid_eta_rho)
        status=nf_def_dim(ncid,'xi_u'   ,L           ,dimid_xi_u   )
        status=nf_def_dim(ncid,'eta_u'  ,Mp          ,dimid_eta_u  )
        status=nf_def_dim(ncid,'xi_v'   ,Lp          ,dimid_xi_v   )
        status=nf_def_dim(ncid,'eta_v'  ,M           ,dimid_eta_v  )
        status=nf_def_dim(ncid,'N'      ,N           ,dimid_N      )
        status=nf_def_dim(ncid,'Np'     ,Np          ,dimid_Np     )
        status=nf_def_dim(ncid,'ocean_time'   ,nf_unlimited,dimid_time)

        dimids(1)=dimid_time

        status=nf_def_var(ncid,'ocean_time',nf_double,1,dimids,
     *               varid_ocean_time)

        dimids(1)=dimid_xi_rho
        dimids(2)=dimid_eta_rho
        dimids(3)=dimid_time

        status=nf_def_var(ncid,'zeta',nf_double,3,dimids,
     *               varid_zeta)
        status=nf_put_att_text(ncid,varid_zeta,'time',
     *  10,'ocean_time')

        dimids(1)=dimid_xi_u
        dimids(2)=dimid_eta_u
        dimids(3)=dimid_time

        status=nf_def_var(ncid,'ubar',nf_double,3,dimids,
     *               varid_ubar)
        status=nf_put_att_text(ncid,varid_ubar,'time',
     *  10,'ocean_time')

        dimids(1)=dimid_xi_v
        dimids(2)=dimid_eta_v
        dimids(3)=dimid_time

        status=nf_def_var(ncid,'vbar',nf_double,3,dimids,
     *               varid_vbar)
        status=nf_put_att_text(ncid,varid_vbar,'time',
     *  10,'ocean_time')


        dimids(1)=dimid_xi_u
        dimids(2)=dimid_eta_u
        dimids(3)=dimid_N
        dimids(4)=dimid_time

        status=nf_def_var(ncid,'u',nf_double,4,dimids,
     *               varid_u)
        status=nf_put_att_text(ncid,varid_u,'time',
     *  10,'ocean_time')


        dimids(1)=dimid_xi_v
        dimids(2)=dimid_eta_v
        dimids(3)=dimid_N
        dimids(4)=dimid_time

        status=nf_def_var(ncid,'v',nf_double,4,dimids,
     *               varid_v)
        status=nf_put_att_text(ncid,varid_v,'time',
     *  10,'ocean_time')


        dimids(1)=dimid_xi_rho
        dimids(2)=dimid_eta_rho
        dimids(3)=dimid_N
        dimids(4)=dimid_time

        status=nf_def_var(ncid,'temp',nf_double,4,dimids,
     *               varid_temp)
        status=nf_put_att_text(ncid,varid_temp,'time',
     *  10,'ocean_time')


        status=nf_def_var(ncid,'salt',nf_double,4,dimids,
     *               varid_salt)
        status=nf_put_att_text(ncid,varid_salt,'time',
     *  10,'ocean_time')


        status=nf_close(ncid)
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

      END IF                         ! <- end "if (iic.eq.1)"

      IF (MOD(iic-1,nhis).eq.0) THEN
!-------------------------------------------------------------
! Output arrays from each processor to the single netcdf file, 
! one at a time
!-------------------------------------------------------------
       start=0
       count=0


!--------------------------------------------------
! output ocean_time on the master thread (myrank=0)
!--------------------------------------------------
       IF (myrank.eq.0) THEN

        status=nf_open(TRIM(hisfile),nf_write,ncid)

        status=nf_inq_dimlen(ncid,dimid_time,recnum)
        recnum=recnum+1
         
        start(1)=recnum
        count(1)=1
        status=nf_put_vara_double(ncid,varid_ocean_time,start,count,
     *                                                time)
        status=nf_close(ncid)

        write (*,*) 
        write (*,*) 'SUBROUTINE OUTPUT: output recnum=',recnum

       END IF  ! (myrank.eq.0)

       CALL mpi_bcast(recnum,1,mpi_integer,0,mpi_comm_world,ierr)

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
             
        CALL mpi_bcast(out_done_all,MaxTile,mpi_integer,
     *                 0,mpi_comm_world,ierr)

        IF (myrank.eq.k .and. out_done_all(k).eq.0) THEN
      
!--------------
! output zeta:
!--------------
         status=nf_open(TRIM(hisfile),nf_write,ncid)

! - new: output excludes ghost points (but includes boundary points)
         n1=Iend_E-Istr_E+1
         n2=Jend_E-Jstr_E+1

         allocate(wrk2d(n1,n2))
         wrk2d(1:n1,1:n2)=tl_zeta(Istr_E:Iend_E,Jstr_E:Jend_E,1)
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
! output ubar:
!--------------

! - new: output excludes ghost points (but includes boundary points)
         n1=IendU_E-IstrU_E+1
         n2=JendU_E-JstrU_E+1

         allocate(wrk2d(n1,n2))
         wrk2d(1:n1,1:n2)=tl_ubar(IstrU_E:IendU_E,JstrU_E:JendU_E,kout)
         start(1)=IstrU_E
         start(2)=JstrU_E+1
         start(3)=recnum
         count(1)=n1
         count(2)=n2
         count(3)=1

         status=nf_put_vara_double(ncid,varid_ubar,start,count,wrk2d)
         deallocate(wrk2d)
       
!--------------
! output vbar:
!--------------
       
! - new: output excludes ghost points (but includes boundary points)
         n1=IendV_E-IstrV_E+1
         n2=JendV_E-JstrV_E+1

         allocate(wrk2d(n1,n2))
         wrk2d(1:n1,1:n2)=tl_vbar(IstrV_E:IendV_E,JstrV_E:JendV_E,kout)
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
     *        tl_u(IstrU_E:IendU_E,JstrU_E:JendU_E,1:N,nstp)
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
     *      tl_v(IstrV_E:IendV_E,JstrV_E:JendV_E,1:N,nstp)
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
     *     tl_t(Istr_E:Iend_E,Jstr_E:Jend_E,1:N,nstp,1)
         status=nf_put_vara_double(ncid,varid_temp,start,count,wrk3d)

         wrk3d(1:n1,1:n2,1:N)=
     *     tl_t(Istr_E:Iend_E,Jstr_E:Jend_E,1:N,nstp,2)
         status=nf_put_vara_double(ncid,varid_salt,start,count,wrk3d)
 
         deallocate(wrk3d)

         status=nf_close(ncid)

         out_done=1

        END IF ! (myrank.eq.k .and. out_done_all(k).eq.0)

       END DO ! k=0:MaxTileNo     

       deallocate(out_done_all)

      END IF  ! <- end "if (MOD(iic-1,nhis).eq.0)"

      RETURN
      END

# else /* serial */

      SUBROUTINE tl_output
! Note: output 2D variables (zeta, ubar, vbar) w/ index 1 (such as zeta(:,:,1))
! output 3D variables  w/ index nstp (such as zeta(:,:,nstp))
      USE mod_ocean
      USE mod_nf_param
      USE mod_output

      implicit none
      integer :: status,ncid,recnum
      integer, dimension(4) :: dimids,start,count
      real(8), allocatable :: wrk2d(:,:)
      real(8), allocatable :: wrk3d(:,:,:)

      real(8) :: q0
 
      IF (iic.eq.1) THEN

       dimids=0

! Define the history file structure:       
       status= nf_create(TRIM(hisfile),nf_clobber,ncid)
        
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
     *  10,'ocean_time')


       dimids(1)=his_dimid%xi_u
       dimids(2)=his_dimid%eta_u
       dimids(3)=his_dimid%time

       status=nf_def_var(ncid,'ubar',nf_double,3,dimids,
     *               his_varid%ubar)
        status=nf_put_att_text(ncid,his_varid%ubar,'time',
     *  10,'ocean_time')


       dimids(1)=his_dimid%xi_v
       dimids(2)=his_dimid%eta_v
       dimids(3)=his_dimid%time

       status=nf_def_var(ncid,'vbar',nf_double,3,dimids,
     *               his_varid%vbar)
        status=nf_put_att_text(ncid,his_varid%vbar,'time',
     *  10,'ocean_time')


       dimids(1)=his_dimid%xi_u
       dimids(2)=his_dimid%eta_u
       dimids(3)=his_dimid%N
       dimids(4)=his_dimid%time

       status=nf_def_var(ncid,'u',nf_double,4,dimids,
     *               his_varid%u)
        status=nf_put_att_text(ncid,his_varid%u,'time',
     *  10,'ocean_time')


       dimids(1)=his_dimid%xi_v
       dimids(2)=his_dimid%eta_v
       dimids(3)=his_dimid%N
       dimids(4)=his_dimid%time

       status=nf_def_var(ncid,'v',nf_double,4,dimids,
     *               his_varid%v)
        status=nf_put_att_text(ncid,his_varid%v,'time',
     *  10,'ocean_time')


       dimids(1)=his_dimid%xi_rho
       dimids(2)=his_dimid%eta_rho
       dimids(3)=his_dimid%N
       dimids(4)=his_dimid%time

       status=nf_def_var(ncid,'temp',nf_double,4,dimids,
     *               his_varid%temp)
        status=nf_put_att_text(ncid,his_varid%temp,'time',
     *  10,'ocean_time')


       status=nf_def_var(ncid,'salt',nf_double,4,dimids,
     *               his_varid%salt)
        status=nf_put_att_text(ncid,his_varid%salt,'time',
     *  10,'ocean_time')


       status=nf_close(ncid)
      END IF                         ! <- end "if (iic.eq.1)"

      IF (MOD(iic-1,nhis).eq.0) THEN
       start=0
       count=0

       status=nf_open(TRIM(hisfile),nf_write,ncid)

       status=nf_inq_dimlen(ncid,his_dimid%time,recnum)
       recnum=recnum+1
         
!--------------
! output ocean_time:
!--------------
       start(1)=recnum
       count(1)=1
       status=nf_put_vara_double(ncid,his_varid%ocean_time,start,count,
     *                                                time)

!--------------
! output zeta:
!--------------
       allocate(wrk2d(Lp,Mp))
       wrk2d(1:Lp,1:Mp)=tl_zeta(0:L,0:M,1)
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
       wrk2d(1:L,1:Mp)=tl_ubar(1:L,0:M,kout)
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
       wrk2d(1:Lp,1:M)=tl_vbar(0:L,1:M,kout)
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
       wrk3d(1:L,1:Mp,1:N)=tl_u(1:L,0:M,1:N,nstp)
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
       wrk3d(1:Lp,1:M,1:N)=tl_v(0:L,1:M,1:N,nstp)
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

       wrk3d(1:Lp,1:Mp,1:N)=tl_t(0:L,0:M,1:N,nstp,1)
       status=nf_put_vara_double(ncid,his_varid%temp,start,count,wrk3d)

       wrk3d(1:Lp,1:Mp,1:N)=tl_t(0:L,0:M,1:N,nstp,2)
       status=nf_put_vara_double(ncid,his_varid%salt,start,count,wrk3d)

       deallocate(wrk3d)

       status=nf_close(ncid)

      END IF  ! <- end "if (MOD(iic-1,nhis).eq.0)"

      RETURN
      END

#endif /* MPI or Serial */ 
      
#endif /* TANGENT */

!------------------------------------------------
      SUBROUTINE test_output_tiles(a,LBi,UBi,LBj,UBj,
     *                             myrank,coords,fhead)

      implicit none
      integer :: LBi,UBi,LBj,UBj,coords(2),myrank
      real(8) :: a(LBi:UBi,LBj:UBj)
      character(len=80) :: fhead,fname
      character(len=2)  :: cc
      integer :: i,j,idev,flen

      if (coords(1).eq.0 .and. coords(2).eq.0) then
       cc='00'
      elseif (coords(1).eq.0 .and. coords(2).eq.1) then
       cc='01'
      elseif (coords(1).eq.1 .and. coords(2).eq.0) then
       cc='10'
      elseif (coords(1).eq.1 .and. coords(2).eq.1) then
       cc='11'
      endif
       
      fname=''
      flen=len_trim(fhead)
      fname(1:flen)=fhead(1:flen)
      fname(flen+1:flen+2)=cc
      fname(flen+3:flen+6)='.dat'

      write (*,'(a80)') fname

      idev=80+myrank
      open(idev,file=TRIM(fname))
      do i=LBi,UBi
      do j=LBj,UBj
       write(idev,*) i,j,a(i,j)
      end do
      end do
      close(idev)

      RETURN
      END

!------------------------------------------------
      SUBROUTINE test_output_tiles_3d(a,LBi,UBi,LBj,UBj,LBk,UBk,
     *                             myrank,coords,fhead)

      implicit none
      integer :: LBi,UBi,LBj,UBj,LBk,UBk,coords(2),myrank
      real(8) :: a(LBi:UBi,LBj:UBj,LBk:UBk)
      character(len=80) :: fhead,fname
      character(len=2)  :: cc
      integer :: i,j,k,idev,flen

      if (coords(1).eq.0 .and. coords(2).eq.0) then
       cc='00'
      elseif (coords(1).eq.0 .and. coords(2).eq.1) then
       cc='01'
      elseif (coords(1).eq.1 .and. coords(2).eq.0) then
       cc='10'
      elseif (coords(1).eq.1 .and. coords(2).eq.1) then
       cc='11'
      endif
       
      fname=''
      flen=len_trim(fhead)
      fname(1:flen)=fhead(1:flen)
      fname(flen+1:flen+2)=cc
      fname(flen+3:flen+6)='.dat'

      write (*,'(a80)') fname

      idev=80+myrank
      open(idev,file=TRIM(fname))
      do k=LBk,UBk
      do j=LBj,UBj
      do i=LBi,UBi
       write(idev,*) i,j,k,a(i,j,k)
      end do
      end do
      end do
      close(idev)

      RETURN
      END

