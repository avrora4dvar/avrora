#include "cppdefs.h"

#ifdef TANGENT
      PROGRAM TL_AVRORA
      USE mod_ocean
      USE mod_parallel

      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
      integer :: indx
      logical :: PREDICTOR_2D

      integer :: ierr,i,j,k
      character (len=80) :: fhead

#ifdef MPI
      CALL mpi_init(ierr)      
#endif

      CALL inp_par

#ifdef MPI
      MaxTile=ntiles_xi*ntiles_eta-1 ! declared in mod_param.F
      dims_mpi(1)=ntiles_xi
      dims_mpi(2)=ntiles_eta
      periods_mpi(1)=.false.
      periods_mpi(2)=.false.
      reorder=.true.
      CALL mpi_cart_create(mpi_comm_world,2,dims_mpi,periods_mpi,
     *                     reorder,comm2d,ierr)
      call mpi_comm_rank(comm2d,myrank,ierr)
      call mpi_cart_coords(comm2d,myrank,2,coords,ierr)
      CALL mpi_cart_shift(comm2d,0,1,rankleft,rankright,ierr)
      CALL mpi_cart_shift(comm2d,1,1,rankbottom,ranktop,ierr)
      call mpi_comm_size(comm2d,comm_size,ierr)
#endif

      CALL set_bounds

#ifdef MPI
! 2d:
      call mpi_Type_vector( UBj-LBj+1,1,UBi-LBi+1,  
     *                      MPI_DOUBLE_PRECISION, stride, ierr ) 
      call mpi_Type_commit( stride, ierr ) 

! r3d, bottom-top:
      call mpi_Type_vector( N,Iend_E-Istr_E+1,
     *                      (UBj-LBj+1)*(UBi-LBi+1), 
     *                      MPI_DOUBLE_PRECISION, stride_r3d_bt, ierr) 
      call mpi_Type_commit( stride_r3d_bt, ierr ) 

! u3d, bottom-top:
      call mpi_Type_vector( N,IendU_E-IstrU_E+1,
     *                      (UBj-LBj+1)*(UBi-LBi+1),  
     *                      MPI_DOUBLE_PRECISION, stride_u3d_bt, ierr) 
      call mpi_Type_commit( stride_u3d_bt, ierr ) 

! v3d, bottom-top:
      call mpi_Type_vector( N,IendV_E-IstrV_E+1,
     *                      (UBj-LBj+1)*(UBi-LBi+1),  
     *                      MPI_DOUBLE_PRECISION, stride_v3d_bt, ierr) 
      call mpi_Type_commit( stride_v3d_bt, ierr ) 

! 3d, left to right:
      call mpi_Type_vector( N*(UBj-LBj+1),1,UBi-LBi+1,  
     *                      MPI_DOUBLE_PRECISION, stride_3d_lr, ierr) 
      call mpi_Type_commit( stride_3d_lr, ierr ) 
#endif

      CALL allocate_arrays
      CALL read_grid
      CALL read_frc_time
      CALL read_fwd_time
      CALL set_weights

# if defined UV_VIS2 || defined TS_DIF2
      CALL set_visc2
# endif 

      CALL tl_initial   ! ini tl_zeta, tl_u, tl_v, tl_t, tl_z-coord & tl_Hz
      iic=0                      ! <- index to initialize stateG 

! Note:
! set_state sets: 
! fwd_zeta(:,:,indx)
! fwd_ubar(:,:,indx)
! fwd_vbar(:,:,indx)
! u(:,:,:,indx)
! v(:,:,:,indx) 

! FWD FIELDS AT INI TIME:
      indx=1
      CALL get_state             ! finds window for time+dt, init zetaG etc.
      CALL set_state(time0,indx) ! state(indx) at ini "slow" time

      CALL set_depth(LBi,UBi,LBj,UBj,N,Lm,Mm,
     *               z0_r,z0_w,oh,
     *               fwd_zeta(:,:,indx),z_r,z_w,Hz,Hz_u,Hz_v) 

      CALL fwd_set_massflux(indx)

! tl_ubar(:,:,1), tl_vbar(:,:,1):
      CALL tl_init_ubar 

! Next, the output index for tl_ubar,tl_vbar is initialized
! (will be diff depending on whether nfast is odd or even
! (note: probably by overlook, in Avrora_v2.1, tl_ubar and vbar 
! were output at indx=1. This choice would not affect ADJ symmetry
! since tl_ubar,tl_vbar are not part of the TL state vector.

      kout=1  ! <- output index for tl_ubar,tl_vbar, for the 1st output

      DO iic=1,ntimes+1

!       write (*,*) '***'
!       write (*,*) 'iic=',iic
!       write (*,*) '***'

!** Indices, b/c:
!** iic : nstp  nrhs nnew
!    1      1    1    2
!    2      2    2    1
!    3      1    1    2
!    4      2    2    1
! 
       nstp=1+MOD(iic-1,2)
       nnew=3-nstp
       nrhs=nstp

       time=time0+dt*DBLE(iic-1)    ! sim to ROMS:0, dt, 2dt, ...

       CALL tl_get_data
       CALL tl_set_data
       CALL tl_set_massflux
       CALL tl_rho_eos
       CALL tl_set_vbc
       CALL fwd_omega(nstp)  ! -> WN, W
       CALL tl_omega         ! depends on fwd WN, z_w, tl_Hudy, tl_z_w
       CALL set_intermediate_tracer ! -> t(3)=t(nstp), used in tl_pre_step
       CALL set_n_minus_1_tracer ! -> t(nnew)=t(nstp), used in tl_pre_step
       CALL tl_rhs3d         ! use z_w, W, Hudy, etc. at t=time

       CALL tl_init_zeta     ! tl_zeta(i,j,1)=tl_Zt_avg1(i,j)
       CALL tl_output
 
       IF (iic.eq.ntimes+1) THEN
          EXIT
       END IF

! FWD FIELDS, NNEW:    
       CALL get_state
       CALL set_state(time+dt,nnew) ! state(nnew), at t=time+dt

#ifndef ANA_FWD_MIX
! Akv and Akt at time nnew (to be used only in tl_step3d_uv and tl_step3d_t)
       CALL set_fwd_mix(time+dt)
#endif

       CALL set_depth(LBi,UBi,LBj,UBj,N,Lm,Mm,
     *       z0_r,z0_w,oh,
     *       fwd_zeta(:,:,nnew),z_r,z_w,Hz,Hz_u,Hz_v) ! at t=time+dt ("nnew")
 
       CALL fwd_set_massflux(nnew)

       CALL fwd_fast_vars
       
       DO iif=1,nfast+1

! INDICES: kstp krhs knew
!   P1:     1    1    3
!   C1:     1    3    2
!-
!   P2:     1    2    3
!   C2:     2    3    1
!-
!   P3:     2    1    3
!   C3:     1    3    2
!-
!   P4:     1    2    3
!   C4:     2    3    1
!-
!   P5:     2    1    3
!   C5:     1    3    2

        PREDICTOR_2D=.TRUE.
        knew=3
        IF (iif.eq.1) THEN
         kstp=1
         krhs=1
        ELSE
         IF (mod(iif,2).eq.0) THEN 
          kstp=1
          krhs=2
         ELSE
          kstp=2
          krhs=1
         END IF
        END IF

        CALL tl_step2d(PREDICTOR_2D)

        PREDICTOR_2D=.FALSE.
        IF (mod(iif,2).eq.1) THEN ! iif odd, CORR step
         kstp=1
         krhs=3
         knew=2
        ELSE                      ! iif even, CORR step
         kstp=2
         krhs=3
         knew=1
        END IF

        if (iif.eq.nfast) then
         kout=knew
        end if

        CALL tl_step2d(PREDICTOR_2D)

       END DO

#ifndef TS_FIXED 
       CALL fwd_set_massflux_trc
#endif

! Note: tl_step3d_uv
! obtains tl_Hudy, tl_Hvdx for tracer computation, 
! but does not use Hudy(:::3), Hvdx(:::3) updated in 
! fwd_set_massflux_trc
       CALL tl_step3d_uv    

#ifndef TS_FIXED 
! Here, fwd_omega uses Hudy(:::3), Hvdx(:::3) and Hz at nnew time
       CALL fwd_omega(3)  ! uses Hudy(3), z_w(at t=t_nnew)

       CALL tl_omega      ! depends on fwd WN, z_w("nnew"),tl_Hudy, tl_z_w    

       CALL tl_step3d_t
#endif

!        if (iic.eq.1) then      
!         fhead='tmp_1_'
!         CALL test_output_tiles(tl_Zt_avg1(LBi:UBi,LBj:UBj),
!     *                       LBi,UBi,LBj,UBj,
!     *                       myrank,coords,fhead)
!         fhead='tmp_2_'
!         CALL test_output_tiles(tl_ubar(LBi:UBi,LBj:UBj,1),
!     *                       LBi,UBi,LBj,UBj,
!     *                       myrank,coords,fhead)
!         fhead='tmp_3_'
!         CALL test_output_tiles(tl_vbar(LBi:UBi,LBj:UBj,1),
!     *                       LBi,UBi,LBj,UBj,
!     *                       myrank,coords,fhead)
!         fhead='tmp_4_'
!         CALL test_output_tiles(tl_u(LBi:UBi,LBj:UBj,4,nnew),
!     *                       LBi,UBi,LBj,UBj,
!     *                       myrank,coords,fhead)
!         fhead='tmp_5_'
!         CALL test_output_tiles(tl_v(LBi:UBi,LBj:UBj,4,nnew),
!     *                       LBi,UBi,LBj,UBj,
!     *                       myrank,coords,fhead)
!         fhead='tmp_6_'
!         CALL test_output_tiles(tl_t(LBi:UBi,LBj:UBj,4,nnew,1),
!     *                       LBi,UBi,LBj,UBj,
!     *                       myrank,coords,fhead)
!         fhead='tmp_7_'
!         CALL test_output_tiles(tl_t(LBi:UBi,LBj:UBj,4,nnew,2),
!     *                       LBi,UBi,LBj,UBj,
!     *                       myrank,coords,fhead)
!        endif 

      END DO ! iic=1,ntimes+1

#ifdef MPI
        CALL mpi_barrier(comm2d,ierr)
        IF(myrank.eq.0) WRITE(*,*) 'TL_AVRORA done'
        CALL mpi_finalize(ierr)
#else
        WRITE(*,*) 'TL_AVRORA done'
#endif       

      END PROGRAM tl_avrora

#endif /* TANGENT */
