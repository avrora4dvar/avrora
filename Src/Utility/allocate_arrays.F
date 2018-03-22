#include "cppdefs.h"

      SUBROUTINE allocate_arrays
      USE mod_ocean
      USE mod_forces

      implicit none

! Allocate grid related arrays
      allocate(dx(LBi:UBi,LBj:UBj))
      allocate(dy(LBi:UBi,LBj:UBj))
      allocate( h(LBi:UBi,LBj:UBj))
      allocate(oh(LBi:UBi,LBj:UBj))
      allocate( f(LBi:UBi,LBj:UBj))
      allocate(mask_r(LBi:UBi,LBj:UBj))
      allocate(mask_u(LBi:UBi,LBj:UBj))
      allocate(mask_v(LBi:UBi,LBj:UBj))

      dx=0.d0
      dy=0.d0
      h=0.d0
      oh=0.d0
      f=0.d0
      mask_r=0
      mask_u=0
      mask_v=0

#if defined CURVGRID && defined UV_ADV
      allocate(dmde(LBi:UBi,LBj:UBj))
      allocate(dndx(LBi:UBi,LBj:UBj))
      dmde=0.d0
      dndx=0.d0
#endif

      allocate(z0_r(LBi:UBi,LBj:UBj,N))
      allocate(z0_w(LBi:UBi,LBj:UBj,0:N))
      allocate(z_r(LBi:UBi,LBj:UBj,N))
      allocate(z_w(LBi:UBi,LBj:UBj,0:N))
      allocate(Hz(LBi:UBi,LBj:UBj,N))
      allocate(Hz_u(LBi:UBi,LBj:UBj,N))
      allocate(Hz_v(LBi:UBi,LBj:UBj,N))

      z0_r=0.d0
      z0_w=0.d0
      z_r=0.d0
      z_w=0.d0
      Hz=0.d0
      Hz_u=0.d0
      Hz_v=0.d0

      allocate(Akv(LBi:UBi,LBj:UBj,0:N))
      allocate(Akt(LBi:UBi,LBj:UBj,0:N,2))
      Akv=Akv0
      Akt(:,:,:,1)=Akt0(1)
      Akt(:,:,:,2)=Akt0(2)

      allocate(dx_u(LBi:UBi,LBj:UBj))
      allocate(dy_u(LBi:UBi,LBj:UBj))
      allocate(dx_v(LBi:UBi,LBj:UBj))
      allocate(dy_v(LBi:UBi,LBj:UBj))

      dx_u=0.d0
      dy_u=0.d0
      dx_v=0.d0
      dy_v=0.d0

      allocate(Area(LBi:UBi,LBj:UBj))
      allocate(Area_u(LBi:UBi,LBj:UBj))
      allocate(Area_v(LBi:UBi,LBj:UBj))
      allocate(oArea(LBi:UBi,LBj:UBj))
      allocate(oArea_u(LBi:UBi,LBj:UBj))
      allocate(oArea_v(LBi:UBi,LBj:UBj))

      Area=0.d0
      Area_u=0.d0
      Area_v=0.d0
      oArea=0.d0
      oArea_u=0.d0
      oArea_v=0.d0

      allocate(visc2_r(LBi:UBi,LBj:UBj))
      allocate(visc2_p(LBi:UBi,LBj:UBj))

      visc2_r=0.d0
      visc2_p=0.d0

      allocate(dxody(LBi:UBi,LBj:UBj))
      allocate(dyodx(LBi:UBi,LBj:UBj))
      allocate(dxody_p(LBi:UBi,LBj:UBj))
      allocate(dyodx_p(LBi:UBi,LBj:UBj))
      allocate(odx(LBi:UBi,LBj:UBj))
      allocate(ody(LBi:UBi,LBj:UBj))
      allocate(dx_p(LBi:UBi,LBj:UBj))
      allocate(dy_p(LBi:UBi,LBj:UBj))

      dxody=0.d0
      dyodx=0.d0
      dxody_p=0.d0
      dyodx_p=0.d0
      odx=0.d0
      ody=0.d0
      dx_p=0.d0
      dy_p=0.d0

      allocate(pmask(LBi:UBi,LBj:UBj))
      pmask=0.d0

#ifdef TS_DIF2
      allocate(diff2(LBi:UBi,LBj:UBj,Ntr))
      allocate(dxody_v(LBi:UBi,LBj:UBj))
      allocate(dyodx_u(LBi:UBi,LBj:UBj))
 
      diff2=0.d0
      dxody_v=0.d0
      dyodx_u=0.d0
#endif

#if defined UV_COR || defined UV_COR_MASK 
      allocate(fdxdy(LBi:UBi,LBj:UBj))
      allocate(omask_uu(LBi:UBi,LBj:UBj))
      allocate(omask_vv(LBi:UBi,LBj:UBj))
      fdxdy=0.d0
      omask_uu=0.d0
      omask_vv=0.d0
#endif

! Allocate state arrays
! - nonlinear / basic  state: 
      allocate(zeta(LBi:UBi,LBj:UBj,3))
      allocate(ubar(LBi:UBi,LBj:UBj,3))
      allocate(vbar(LBi:UBi,LBj:UBj,3))
      allocate(u(LBi:UBi,LBj:UBj,N,2))
      allocate(v(LBi:UBi,LBj:UBj,N,2))
      allocate(t(LBi:UBi,LBj:UBj,N,3,Ntr))
      allocate(w(LBi:UBi,LBj:UBj,0:N))
      allocate(rho(LBi:UBi,LBj:UBj,N))

      zeta=0.d0
      ubar=0.d0
      vbar=0.d0
      u=0.d0
      v=0.d0
      t=0.d0
      w=0.d0
      rho=rho0

      allocate(UFx(LBi:UBi,LBj:UBj))
      allocate(VFx(LBi:UBi,LBj:UBj))
      allocate(UFe(LBi:UBi,LBj:UBj))
      allocate(VFe(LBi:UBi,LBj:UBj))
      UFx=0.d0
      VFx=0.d0
      UFe=0.d0
      VFe=0.d0

      allocate(DU_avg1(LBi:UBi,LBj:UBj))
      allocate(DV_avg1(LBi:UBi,LBj:UBj))

      DU_avg1=0.d0
      DV_avg1=0.d0
 
      allocate(DU(LBi:UBi,LBj:UBj))
      allocate(DV(LBi:UBi,LBj:UBj))
      DU=0.d0
      DV=0.d0

#if defined NONLINEAR  

      allocate(DU_avg2(LBi:UBi,LBj:UBj))
      allocate(DV_avg2(LBi:UBi,LBj:UBj))
      DU_avg2=0.d0
      DV_avg2=0.d0

      allocate(Zt_avg1(LBi:UBi,LBj:UBj))
      Zt_avg1=0.d0

      allocate(Hudy(LBi:UBi,LBj:UBj,N))
      allocate(Hvdx(LBi:UBi,LBj:UBj,N))
      Hudy=0.d0
      Hvdx=0.d0

      allocate(ru(LBi:UBi,LBj:UBj,0:N,2))
      allocate(rv(LBi:UBi,LBj:UBj,0:N,2))
      ru=0.d0
      rv=0.d0

      allocate(rufrc(LBi:UBi,LBj:UBj))
      allocate(rvfrc(LBi:UBi,LBj:UBj))
      rufrc=0.d0
      rvfrc=0.d0

      allocate(rhs_zeta(LBi:UBi,LBj:UBj))
      allocate(rhs_ubar(LBi:UBi,LBj:UBj))
      allocate(rhs_vbar(LBi:UBi,LBj:UBj))
      rhs_zeta=0.d0
      rhs_ubar=0.d0
      rhs_vbar=0.d0

      allocate(rzeta(LBi:UBi,LBj:UBj,3))
      allocate(rubar(LBi:UBi,LBj:UBj,3))
      allocate(rvbar(LBi:UBi,LBj:UBj,3))
      rzeta=0.d0
      rubar=0.d0
      rvbar=0.d0

      allocate(sustr(LBi:UBi,LBj:UBj))
      allocate(svstr(LBi:UBi,LBj:UBj))
      allocate(sustrG(LBi:UBi,LBj:UBj,2))
      allocate(svstrG(LBi:UBi,LBj:UBj,2))
      sustr=0.d0
      svstr=0.d0
      sustrG=0.d0
      svstrG=0.d0

      allocate(bustr(LBi:UBi,LBj:UBj))
      allocate(bvstr(LBi:UBi,LBj:UBj))
      bustr=0.d0
      bvstr=0.d0

      allocate(stflx(LBi:UBi,LBj:UBj,Ntr))
      stflx=0.d0

#elif defined TANGENT
!
      allocate(Hudy(LBi:UBi,LBj:UBj,N,3))
      allocate(Hvdx(LBi:UBi,LBj:UBj,N,3))
      Hudy=0.d0
      Hvdx=0.d0

      allocate(tl_DU(LBi:UBi,LBj:UBj))
      allocate(tl_DV(LBi:UBi,LBj:UBj))
      tl_DU=0.d0
      tl_DV=0.d0

      allocate(tl_Hudy(LBi:UBi,LBj:UBj,N))
      allocate(tl_Hvdx(LBi:UBi,LBj:UBj,N))
      tl_Hudy=0.d0
      tl_Hvdx=0.d0

      allocate(tl_z_w(LBi:UBi,LBj:UBj,0:N))
      allocate(tl_z_r(LBi:UBi,LBj:UBj,N))
      allocate(tl_Hz(LBi:UBi,LBj:UBj,N))
      allocate(tl_Hz_u(LBi:UBi,LBj:UBj,N))
      allocate(tl_Hz_v(LBi:UBi,LBj:UBj,N))

      tl_z_w=0.d0
      tl_z_r=0.d0
      tl_Hz=0.d0
      tl_Hz_u=0.d0
      tl_Hz_v=0.d0

      allocate(tl_ru(LBi:UBi,LBj:UBj,0:N,2))
      allocate(tl_rv(LBi:UBi,LBj:UBj,0:N,2))
      tl_ru=0.d0
      tl_rv=0.d0

      allocate(tl_rufrc(LBi:UBi,LBj:UBj))
      allocate(tl_rvfrc(LBi:UBi,LBj:UBj))
      tl_rufrc=0.d0
      tl_rvfrc=0.d0

      allocate(tl_rhs_zeta(LBi:UBi,LBj:UBj))
      allocate(tl_rhs_ubar(LBi:UBi,LBj:UBj))
      allocate(tl_rhs_vbar(LBi:UBi,LBj:UBj))
      tl_rhs_zeta=0.d0
      tl_rhs_ubar=0.d0
      tl_rhs_vbar=0.d0

      allocate(tl_rzeta(LBi:UBi,LBj:UBj,3))
      allocate(tl_rubar(LBi:UBi,LBj:UBj,3))
      allocate(tl_rvbar(LBi:UBi,LBj:UBj,3))
      tl_rzeta=0.d0
      tl_rubar=0.d0
      tl_rvbar=0.d0

      allocate(tl_zeta(LBi:UBi,LBj:UBj,3))
      allocate(tl_ubar(LBi:UBi,LBj:UBj,3))
      allocate(tl_vbar(LBi:UBi,LBj:UBj,3))
      allocate(tl_u(LBi:UBi,LBj:UBj,N,2))
      allocate(tl_v(LBi:UBi,LBj:UBj,N,2))
      allocate(tl_t(LBi:UBi,LBj:UBj,N,3,Ntr))
      allocate(tl_w(LBi:UBi,LBj:UBj,0:N))
      allocate(tl_rho(LBi:UBi,LBj:UBj,N))
      tl_zeta=0.d0
      tl_ubar=0.d0
      tl_vbar=0.d0
      tl_u=0.d0
      tl_v=0.d0
      tl_t=0.d0
      tl_w=0.d0
      tl_rho=0.d0

      allocate(tl_sustr(LBi:UBi,LBj:UBj))
      allocate(tl_svstr(LBi:UBi,LBj:UBj))
      allocate(tl_sustrG(LBi:UBi,LBj:UBj,2))
      allocate(tl_svstrG(LBi:UBi,LBj:UBj,2))
      tl_sustr=0.d0
      tl_svstr=0.d0
      tl_sustrG=0.d0
      tl_svstrG=0.d0

      allocate(tl_bustr(LBi:UBi,LBj:UBj))
      allocate(tl_bvstr(LBi:UBi,LBj:UBj))
      tl_bustr=0.d0
      tl_bvstr=0.d0

      allocate(tl_stflx(LBi:UBi,LBj:UBj,Ntr))
      tl_stflx=0.d0

      allocate(tl_UFx(LBi:UBi,LBj:UBj))
      allocate(tl_VFx(LBi:UBi,LBj:UBj))
      allocate(tl_UFe(LBi:UBi,LBj:UBj))
      allocate(tl_VFe(LBi:UBi,LBj:UBj))
      tl_UFx=0.d0
      tl_VFx=0.d0
      tl_UFe=0.d0
      tl_VFe=0.d0

      allocate(tl_Zt_avg1(LBi:UBi,LBj:UBj))
      allocate(tl_DU_avg1(LBi:UBi,LBj:UBj))
      allocate(tl_DU_avg2(LBi:UBi,LBj:UBj))
      allocate(tl_DV_avg1(LBi:UBi,LBj:UBj))
      allocate(tl_DV_avg2(LBi:UBi,LBj:UBj))

      tl_Zt_avg1=0.d0
      tl_DU_avg1=0.d0
      tl_DU_avg2=0.d0
      tl_DV_avg1=0.d0
      tl_DV_avg2=0.d0

#elif defined ADJOINT
!
      allocate(Hudy(LBi:UBi,LBj:UBj,N,3))
      allocate(Hvdx(LBi:UBi,LBj:UBj,N,3))
      Hudy=0.d0
      Hvdx=0.d0

      allocate(ad_DU(LBi:UBi,LBj:UBj))
      allocate(ad_DV(LBi:UBi,LBj:UBj))
      ad_DU=0.d0
      ad_DV=0.d0

      allocate(ad_Hudy(LBi:UBi,LBj:UBj,N))
      allocate(ad_Hvdx(LBi:UBi,LBj:UBj,N))
      ad_Hudy=0.d0
      ad_Hvdx=0.d0

      allocate(ad_z_w(LBi:UBi,LBj:UBj,0:N))
      allocate(ad_z_r(LBi:UBi,LBj:UBj,N))
      allocate(ad_Hz(LBi:UBi,LBj:UBj,N))
      allocate(ad_Hz_u(LBi:UBi,LBj:UBj,N))
      allocate(ad_Hz_v(LBi:UBi,LBj:UBj,N))

      ad_z_w=0.d0
      ad_z_r=0.d0
      ad_Hz=0.d0
      ad_Hz_u=0.d0
      ad_Hz_v=0.d0

      allocate(ad_ru(LBi:UBi,LBj:UBj,0:N,2))
      allocate(ad_rv(LBi:UBi,LBj:UBj,0:N,2))

      ad_ru=0.d0
      ad_rv=0.d0

      allocate(ad_rufrc(LBi:UBi,LBj:UBj))
      allocate(ad_rvfrc(LBi:UBi,LBj:UBj))
      ad_rufrc=0.d0
      ad_rvfrc=0.d0

      allocate(ad_rhs_zeta(LBi:UBi,LBj:UBj))
      allocate(ad_rhs_ubar(LBi:UBi,LBj:UBj))
      allocate(ad_rhs_vbar(LBi:UBi,LBj:UBj))

      ad_rhs_zeta=0.d0
      ad_rhs_ubar=0.d0
      ad_rhs_vbar=0.d0

      allocate(ad_rzeta(LBi:UBi,LBj:UBj,3))
      allocate(ad_rubar(LBi:UBi,LBj:UBj,3))
      allocate(ad_rvbar(LBi:UBi,LBj:UBj,3))

      ad_rzeta=0.d0
      ad_rubar=0.d0
      ad_rvbar=0.d0

      allocate(ad_zeta(LBi:UBi,LBj:UBj,3))
      allocate(ad_ubar(LBi:UBi,LBj:UBj,3))
      allocate(ad_vbar(LBi:UBi,LBj:UBj,3))
      allocate(ad_u(LBi:UBi,LBj:UBj,N,2))
      allocate(ad_v(LBi:UBi,LBj:UBj,N,2))
      allocate(ad_t(LBi:UBi,LBj:UBj,N,3,Ntr))
      allocate(ad_w(LBi:UBi,LBj:UBj,0:N))
      allocate(ad_rho(LBi:UBi,LBj:UBj,N))
      ad_zeta=0.d0
      ad_ubar=0.d0
      ad_vbar=0.d0
      ad_u=0.d0
      ad_v=0.d0
      ad_t=0.d0
      ad_w=0.d0
      ad_rho=0.d0

      allocate(ad_sustr(LBi:UBi,LBj:UBj))
      allocate(ad_svstr(LBi:UBi,LBj:UBj))
      allocate(ad_sustrG(LBi:UBi,LBj:UBj,2))
      allocate(ad_svstrG(LBi:UBi,LBj:UBj,2))

      ad_sustr=0.d0
      ad_svstr=0.d0
      ad_sustrG=0.d0
      ad_svstrG=0.d0

      allocate(ad_bustr(LBi:UBi,LBj:UBj))
      allocate(ad_bvstr(LBi:UBi,LBj:UBj))
      ad_bustr=0.d0
      ad_bvstr=0.d0

      allocate(ad_stflx(LBi:UBi,LBj:UBj,Ntr))
      ad_stflx=0.d0

      allocate(ad_UFx(LBi:UBi,LBj:UBj))
      allocate(ad_VFx(LBi:UBi,LBj:UBj))
      allocate(ad_UFe(LBi:UBi,LBj:UBj))
      allocate(ad_VFe(LBi:UBi,LBj:UBj))
      ad_UFx=0.d0
      ad_VFx=0.d0
      ad_UFe=0.d0
      ad_VFe=0.d0

      allocate(ad_Zt_avg1(LBi:UBi,LBj:UBj))
      allocate(ad_DU_avg1(LBi:UBi,LBj:UBj))
      allocate(ad_DU_avg2(LBi:UBi,LBj:UBj))
      allocate(ad_DV_avg1(LBi:UBi,LBj:UBj))
      allocate(ad_DV_avg2(LBi:UBi,LBj:UBj))

      ad_Zt_avg1=0.d0
      ad_DU_avg1=0.d0
      ad_DU_avg2=0.d0
      ad_DV_avg1=0.d0
      ad_DV_avg2=0.d0
#endif

#if defined TANGENT || defined ADJOINT
! Allocate background (fwd) G - fields (snapshots for interpolation)
      allocate(ubarG(LBi:UBi,LBj:UBj,2))
      allocate(vbarG(LBi:UBi,LBj:UBj,2))
      allocate(zetaG(LBi:UBi,LBj:UBj,2))
      allocate(tG(LBi:UBi,LBj:UBj,N,2,Ntr))
      allocate(uG(LBi:UBi,LBj:UBj,N,2))
      allocate(vG(LBi:UBi,LBj:UBj,N,2))

      ubarG=0.d0
      vbarG=0.d0
      zetaG=0.d0
      tG=0.d0
      uG=0.d0
      vG=0.d0

#ifndef ANA_FWD_MIX
      allocate(AkvG(LBi:UBi,LBj:UBj,0:N,2))
      allocate(AktG(LBi:UBi,LBj:UBj,0:N,2))
      AkvG=Akv0
      AktG=Akt0(1)
#endif

      allocate(WN(LBi:UBi,LBj:UBj))
      WN=0.d0

      allocate(fwd_zeta(LBi:UBi,LBj:UBj,2))
      allocate(fwd_ubar(LBi:UBi,LBj:UBj,2))
      allocate(fwd_vbar(LBi:UBi,LBj:UBj,2))
 
      fwd_zeta=0.d0
      fwd_ubar=0.d0
      fwd_vbar=0.d0

#endif

! Weights for b/t time-stepping:
      allocate(weight(2,2*ndtfast))
      weight=0.d0

      RETURN
      END
