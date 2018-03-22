#include "cppdefs.h"
      MODULE mod_ocean
      implicit none

      real(8) :: time0,time,dt,dtfast
      integer :: nnew, nstp, nrhs, knew, kstp, krhs,kout
      integer :: iic,iif,ntimes,ndtfast,nfast,nhis
      integer :: Lp,Mp,L,M,Lm,Mm  ! note: xi_rho=Lp, eta_rho=Mp
      integer :: N,Np
      integer :: LBi,UBi,LBj,UBj
      integer :: Istr,Jstr,IstrU,JstrU,IstrV,JstrV
      integer :: Iend,Jend,IendU,JendU,IendV,JendV
#ifdef MPI
! Extended start/end indices for mpi exchange (may include boundary points)
      integer :: Istr_E,Jstr_E,IstrU_E,JstrU_E,IstrV_E,JstrV_E
      integer :: Iend_E,Jend_E,IendU_E,JendU_E,IendV_E,JendV_E
      logical :: B_North,B_South,B_East,B_West
#endif
      integer :: nfwd,fwd_rec
      character(len=80) :: grdfile,inifile,frcfile,fwdfile,hisfile

      integer, parameter :: stdinp=5
      integer, parameter :: Ntr=2

      real(8) :: Akv0,Akt0(Ntr),rdrg,vis2,dif2

      integer :: ntiles_xi,ntiles_eta

! Grid related:
      real(8), allocatable :: dx(:,:)
      real(8), allocatable :: dy(:,:)
#if defined CURVGRID && defined UV_ADV
      real(8), allocatable :: dmde(:,:)
      real(8), allocatable :: dndx(:,:)
#endif
      real(8), allocatable :: h(:,:)
      real(8), allocatable :: oh(:,:)
      integer, allocatable :: mask_r(:,:)
      integer, allocatable :: mask_u(:,:)
      integer, allocatable :: mask_v(:,:)
      real(8), allocatable :: f(:,:)
      real(8), allocatable :: z0_r(:,:,:)
      real(8), allocatable :: z0_w(:,:,:)     
      real(8), allocatable :: z_r(:,:,:)
      real(8), allocatable :: z_w(:,:,:)     

      real(8), allocatable :: Hz(:,:,:)     
      real(8), allocatable :: Hz_u(:,:,:)     
      real(8), allocatable :: Hz_v(:,:,:)     

      real(8), allocatable :: Akv(:,:,:)
      real(8), allocatable :: Akt(:,:,:,:)

      real(8), allocatable :: dx_u(:,:)
      real(8), allocatable :: dy_u(:,:)
      real(8), allocatable :: dx_v(:,:)
      real(8), allocatable :: dy_v(:,:)

      real(8), allocatable :: Area(:,:)
      real(8), allocatable :: Area_u(:,:)
      real(8), allocatable :: Area_v(:,:)
      real(8), allocatable :: oArea(:,:)
      real(8), allocatable :: oArea_u(:,:)
      real(8), allocatable :: oArea_v(:,:)

      real(8), allocatable :: visc2_r(:,:)
      real(8), allocatable :: visc2_p(:,:)

      real(8), allocatable :: dxody(:,:)
      real(8), allocatable :: dyodx(:,:)
      real(8), allocatable :: dxody_p(:,:)
      real(8), allocatable :: dyodx_p(:,:)
      real(8), allocatable :: odx(:,:)
      real(8), allocatable :: ody(:,:)
      real(8), allocatable :: dx_p(:,:)
      real(8), allocatable :: dy_p(:,:)

#if defined TS_DIF2
      real(8), allocatable :: diff2(:,:,:)
      real(8), allocatable :: dxody_v(:,:)
      real(8), allocatable :: dyodx_u(:,:)
#endif

      real(8), allocatable :: pmask(:,:)  ! used in UV_VIS2

#if defined UV_COR || defined UV_COR_MASK 
      real(8), allocatable :: fdxdy(:,:)
      real(8), allocatable :: omask_uu(:,:)  
      real(8), allocatable :: omask_vv(:,:)  
#endif
 
!  Nonlinear model state.
      real(8), allocatable :: ubar(:,:,:)
      real(8), allocatable :: vbar(:,:,:)
      real(8), allocatable :: zeta(:,:,:)
      real(8), allocatable :: t(:,:,:,:,:)
      real(8), allocatable :: u(:,:,:,:)
      real(8), allocatable :: v(:,:,:,:)
      real(8), allocatable :: w(:,:,:)
      real(8), allocatable :: rho(:,:,:)

      real(8), allocatable :: UFx(:,:)
      real(8), allocatable :: VFx(:,:)
      real(8), allocatable :: UFe(:,:)
      real(8), allocatable :: VFe(:,:)

      real(8), allocatable :: DU_avg1(:,:)
      real(8), allocatable :: DV_avg1(:,:)

      real(8), allocatable :: DU(:,:)
      real(8), allocatable :: DV(:,:)

#if defined NONLINEAR      
      real(8), allocatable :: Zt_avg1(:,:)
      real(8), allocatable :: DU_avg2(:,:)
      real(8), allocatable :: DV_avg2(:,:)

      real(8), allocatable :: bustr(:,:)
      real(8), allocatable :: bvstr(:,:)

      real(8), allocatable :: ru(:,:,:,:)
      real(8), allocatable :: rv(:,:,:,:)

      real(8), allocatable :: rufrc(:,:)
      real(8), allocatable :: rvfrc(:,:)
      real(8), allocatable :: rzeta(:,:,:)
      real(8), allocatable :: rubar(:,:,:)
      real(8), allocatable :: rvbar(:,:,:)

      real(8), allocatable :: rhs_zeta(:,:)
      real(8), allocatable :: rhs_ubar(:,:)
      real(8), allocatable :: rhs_vbar(:,:)

      real(8), allocatable :: Hudy(:,:,:) 
      real(8), allocatable :: Hvdx(:,:,:) 

#elif defined TANGENT      
!
      real(8), allocatable :: Hudy(:,:,:,:) 
      real(8), allocatable :: Hvdx(:,:,:,:) 

      real(8), allocatable :: tl_DU(:,:)
      real(8), allocatable :: tl_DV(:,:)

      real(8), allocatable :: tl_Hudy(:,:,:) 
      real(8), allocatable :: tl_Hvdx(:,:,:) 

      real(8), allocatable :: tl_z_w(:,:,:)     
      real(8), allocatable :: tl_z_r(:,:,:)     
      real(8), allocatable :: tl_Hz(:,:,:)     
      real(8), allocatable :: tl_Hz_u(:,:,:)     
      real(8), allocatable :: tl_Hz_v(:,:,:)     

      real(8), allocatable :: tl_ubar(:,:,:)
      real(8), allocatable :: tl_vbar(:,:,:)
      real(8), allocatable :: tl_zeta(:,:,:)
      real(8), allocatable :: tl_t(:,:,:,:,:)
      real(8), allocatable :: tl_u(:,:,:,:)
      real(8), allocatable :: tl_v(:,:,:,:)
      real(8), allocatable :: tl_w(:,:,:)
      real(8), allocatable :: tl_rho(:,:,:)

      real(8), allocatable :: tl_bustr(:,:)
      real(8), allocatable :: tl_bvstr(:,:)

      real(8), allocatable :: tl_UFx(:,:)
      real(8), allocatable :: tl_VFx(:,:)
      real(8), allocatable :: tl_UFe(:,:)
      real(8), allocatable :: tl_VFe(:,:)

      real(8), allocatable :: tl_ru(:,:,:,:)
      real(8), allocatable :: tl_rv(:,:,:,:)

      real(8), allocatable :: tl_rufrc(:,:)
      real(8), allocatable :: tl_rvfrc(:,:)
      real(8), allocatable :: tl_rzeta(:,:,:)
      real(8), allocatable :: tl_rubar(:,:,:)
      real(8), allocatable :: tl_rvbar(:,:,:)

      real(8), allocatable :: tl_rhs_zeta(:,:)
      real(8), allocatable :: tl_rhs_ubar(:,:)
      real(8), allocatable :: tl_rhs_vbar(:,:)

      real(8), allocatable :: tl_Zt_avg1(:,:)
      real(8), allocatable :: tl_DU_avg1(:,:)
      real(8), allocatable :: tl_DU_avg2(:,:)
      real(8), allocatable :: tl_DV_avg1(:,:)
      real(8), allocatable :: tl_DV_avg2(:,:)

#elif defined ADJOINT
!
      real(8), allocatable :: Hudy(:,:,:,:) 
      real(8), allocatable :: Hvdx(:,:,:,:) 

      real(8), allocatable :: ad_DU(:,:)
      real(8), allocatable :: ad_DV(:,:)

      real(8), allocatable :: ad_Hudy(:,:,:) 
      real(8), allocatable :: ad_Hvdx(:,:,:) 

      real(8), allocatable :: ad_z_w(:,:,:)     
      real(8), allocatable :: ad_z_r(:,:,:)     
      real(8), allocatable :: ad_Hz(:,:,:)     
      real(8), allocatable :: ad_Hz_u(:,:,:)     
      real(8), allocatable :: ad_Hz_v(:,:,:)     

      real(8), allocatable :: ad_ubar(:,:,:)
      real(8), allocatable :: ad_vbar(:,:,:)
      real(8), allocatable :: ad_zeta(:,:,:)
      real(8), allocatable :: ad_t(:,:,:,:,:)
      real(8), allocatable :: ad_u(:,:,:,:)
      real(8), allocatable :: ad_v(:,:,:,:)
      real(8), allocatable :: ad_w(:,:,:)
      real(8), allocatable :: ad_rho(:,:,:)

      real(8), allocatable :: ad_bustr(:,:)
      real(8), allocatable :: ad_bvstr(:,:)

      real(8), allocatable :: ad_UFx(:,:)
      real(8), allocatable :: ad_VFx(:,:)
      real(8), allocatable :: ad_UFe(:,:)
      real(8), allocatable :: ad_VFe(:,:)

      real(8), allocatable :: ad_ru(:,:,:,:)
      real(8), allocatable :: ad_rv(:,:,:,:)

      real(8), allocatable :: ad_rufrc(:,:)
      real(8), allocatable :: ad_rvfrc(:,:)
      real(8), allocatable :: ad_rzeta(:,:,:)
      real(8), allocatable :: ad_rubar(:,:,:)
      real(8), allocatable :: ad_rvbar(:,:,:)

      real(8), allocatable :: ad_rhs_zeta(:,:)
      real(8), allocatable :: ad_rhs_ubar(:,:)
      real(8), allocatable :: ad_rhs_vbar(:,:)

      real(8), allocatable :: ad_Zt_avg1(:,:)
      real(8), allocatable :: ad_DU_avg1(:,:)
      real(8), allocatable :: ad_DU_avg2(:,:)
      real(8), allocatable :: ad_DV_avg1(:,:)
      real(8), allocatable :: ad_DV_avg2(:,:)
#endif

#if defined TANGENT || defined ADJOINT
      real(8), allocatable :: fwd_time(:)
      real(8) :: fwd_timeG(2)
      integer :: ifwd1,ifwd2
 
      real(8), allocatable :: ubarG(:,:,:)
      real(8), allocatable :: vbarG(:,:,:)
      real(8), allocatable :: zetaG(:,:,:)
      real(8), allocatable :: tG(:,:,:,:,:)
      real(8), allocatable :: uG(:,:,:,:)
      real(8), allocatable :: vG(:,:,:,:)

# ifndef ANA_FWD_MIX
      real(8), allocatable :: AkvG(:,:,:,:)
      real(8), allocatable :: AktG(:,:,:,:)
# endif

      real(8), allocatable :: WN(:,:) ! background W(:,:,N) w/out zeta-effect

      real(8), allocatable :: fwd_zeta(:,:,:)
      real(8), allocatable :: fwd_ubar(:,:,:)
      real(8), allocatable :: fwd_vbar(:,:,:)

#endif

! Weights for n/t time-stepping:
      real(8), allocatable :: weight(:,:)

! Scalar parameters:
      real(8), parameter :: rho0=1025.0d0
      real(8), parameter :: g=9.81d0

      real(8), parameter :: Falpha = 2.0d0
      real(8), parameter :: Fbeta  = 4.0d0
      real(8), parameter :: Fgamma = 0.284d0

      integer, parameter :: itemp=1
      integer, parameter :: isalt=2

      END MODULE mod_ocean
