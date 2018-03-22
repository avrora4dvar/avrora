#include "cppdefs.h"

#ifdef ADJOINT
      PROGRAM AD_AVRORA
      USE mod_ocean
#ifdef MPI
      USE mod_parallel
#endif

      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
      integer :: indx
      logical :: PREDICTOR_2D
      real(8) :: q0

#ifdef MPI
      integer :: ierr,i,j,k
      character (len=80) :: fhead
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
      call mpi_comm_size(comm2d,comm_size,ierr)
      call mpi_cart_coords(comm2d,myrank,2,coords,ierr)
      CALL mpi_cart_shift(comm2d,0,1,rankleft,rankright,ierr)
      CALL mpi_cart_shift(comm2d,1,1,rankbottom,ranktop,ierr)
#endif

      IF(myrank.eq.0) write(*,*) 'Setting bounds:'
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

      IF(myrank.eq.0) WRITE(*,*) 'Allocating arrays:'
      CALL allocate_arrays
      IF(myrank.eq.0) WRITE(*,*) 'Reading grid:'
      CALL read_grid
      IF(myrank.eq.0) WRITE(*,*) 'Reading time:'
      CALL read_frc_time
      CALL read_fwd_time
      IF(myrank.eq.0) WRITE(*,*) 'Setting weights:'
      CALL set_weights

# if defined UV_VIS2 || defined TS_DIF2
      CALL set_visc2
# endif 

      time0=0.d0 ! later: arrange to read from "hisfile"

      
      DO iic=ntimes+1,1,-1
       
! note: iic |  nstp | nrhs | nnew
!       1   !  1    | 1    | 2 
 
       nstp=1+MOD(iic-1,2)
       nnew=3-nstp
       nrhs=nstp

#ifdef MPI
      CALL mpi_exchange_3d_ext(ad_u(LBi:UBi,LBj:UBj,1:N,nnew),
     *                     LBi,UBi,LBj,UBj,N,
     *                     LBi0,UBi0,LBj0,UBj0,
     *                     IstrU_E,IendU_E,JstrU_E,JendU_E,
     *                     comm2d, stride_u3d_bt, stride_3d_lr, 
     *                     rankleft,rankright,ranktop,rankbottom) 
      CALL mpi_exchange_3d_ext(ad_v(LBi:UBi,LBj:UBj,1:N,nnew),
     *                     LBi,UBi,LBj,UBj,N,
     *                     LBi0,UBi0,LBj0,UBj0,
     *                     IstrV_E,IendV_E,JstrV_E,JendV_E,
     *                     comm2d, stride_v3d_bt, stride_3d_lr, 
     *                     rankleft,rankright,ranktop,rankbottom) 
      CALL mpi_exchange_3d_ext(ad_t(LBi:UBi,LBj:UBj,1:N,nnew,1),
     *                     LBi,UBi,LBj,UBj,N,
     *                     LBi0,UBi0,LBj0,UBj0,
     *                     Istr_E,Iend_E,Jstr_E,Jend_E,
     *                     comm2d, stride_r3d_bt, stride_3d_lr, 
     *                     rankleft,rankright,ranktop,rankbottom) 
      CALL mpi_exchange_3d_ext(ad_t(LBi:UBi,LBj:UBj,1:N,nnew,2),
     *                     LBi,UBi,LBj,UBj,N,
     *                     LBi0,UBi0,LBj0,UBj0,
     *                     Istr_E,Iend_E,Jstr_E,Jend_E,
     *                     comm2d, stride_r3d_bt, stride_3d_lr, 
     *                     rankleft,rankright,ranktop,rankbottom) 
#endif
 
       time=time0+dt*DBLE(iic-1) ! 1st step (iic=ntimes+1):time=T
                                 ! in other words, time=t(nstp)
       
#ifdef MPI  
       IF(myrank.EQ.0) THEN
#endif
       IF(MOD(iic,10).EQ.0) WRITE(*,*) 'time=',time
#ifdef MPI
       END IF
#endif

       CALL get_state_r     
       CALL set_state(time,nstp)

       CALL set_depth(LBi,UBi,LBj,UBj,N,Lm,Mm,
     *               z0_r,z0_w,oh,
     *               fwd_zeta(:,:,nstp),z_r,z_w,Hz,Hz_u,Hz_v) ! at t=time ("nstp")

       CALL fwd_set_massflux(nstp)

       IF (iic.lt.ntimes+1) THEN

#ifndef TS_FIXED 
        CALL fwd_set_massflux_trc ! ave of fluxes at time iic and iic+1
        CALL set_intermediate_tracer ! -> t(3)=t(nstp)
#endif

! reset z levels, nnew time, to compute omega at time iic+1/2
        CALL set_depth(LBi,UBi,LBj,UBj,N,Lm,Mm,
     *               z0_r,z0_w,oh,
     *               fwd_zeta(:,:,nnew),z_r,z_w,Hz,Hz_u,Hz_v) ! at t=time ("nstp")

#ifndef TS_FIXED 
        CALL fwd_omega(3)
! ifndef ANA_FWD_MIX, uses background Akt set at the end of the prev. time step
        CALL ad_step3d_t  ! => update ad_t(:::3:)
        CALL ad_omega
#endif

        CALL fwd_fast_vars ! zeta, ubar, vbar, DU_avg1 at time nnew

! uses z levels, DU_avg1 at time nnew 
! ifndef ANA_FWD_MIX, uses background Akv set at the end of the prev. time step
        CALL ad_step3d_uv

! Fast time steps: use zeta, ubar, vbar at time nnew
        DO iif=nfast+1,1,-1
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
         CALL ad_step2d(PREDICTOR_2D)

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
         CALL ad_step2d(PREDICTOR_2D)

        END DO ! "iif=nfast+1,1,-1"

       END IF ! "IF (iic.lt.ntimes+1) THEN "

! "adj" to tl_output, read snapshots from hisfile, add to ad_u(:,nstp)
       CALL ad_comb_frc ! -> note: forces ad_zeta(:,:,1)

! ad_init_zeta: adj to   tl_zeta(i,j,1)=tl_Zt_avg1(i,j)
! It is essential to implement ad_init_2d at step iic=ntimes+1
       CALL ad_init_zeta 
 
! Note: computation of rhs on step ntimes+1 does not affect the 
! output on that step.
! Also, pay attention, sms_timeG is not defined @ iic=ntimes+1 (see below)
       IF (iic.lt.ntimes+1) THEN 

        CALL set_depth(LBi,UBi,LBj,UBj,N,Lm,Mm,
     *               z0_r,z0_w,oh,
     *               fwd_zeta(:,:,nstp),z_r,z_w,Hz,Hz_u,Hz_v) ! at t=time ("nstp")

        CALL fwd_omega(nstp)        ! background W
        CALL rho_eos                ! background rho

        CALL set_n_minus_1_tracer ! t(nnew)=t(nstp)

        CALL ad_rhs3d
        CALL ad_omega 
        CALL ad_set_vbc
        CALL ad_rho_eos
        CALL ad_set_massflux
        CALL ad_set_data
       END IF ! "IF (iic.lt.ntimes+1) THEN"

       CALL ad_get_data  ! sets sms_timeG on the 1st step (ntimes+1)

#ifndef ANA_FWD_MIX
! Akv and Akt at time nnew (to be used in ad_step3d_uv and ad_step3d_t
! on the next step; do it here before the FWD interpolation window is
! shifted back time)
       CALL set_fwd_mix(time)
#endif
       
      END DO

      CALL ad_init_ubar 

! Output IC correction, adjoint to tl_initial.F:

      CALL ad_ini_out

#ifdef MPI
        CALL mpi_finalize(ierr)
#endif

      END
#endif /* ADJOINT */
