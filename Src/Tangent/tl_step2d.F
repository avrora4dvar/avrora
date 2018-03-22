#include "cppdefs.h"

#ifdef TANGENT
      SUBROUTINE tl_step2d(PREDICTOR)
      USE mod_ocean
#ifdef MPI
      USE mod_parallel
#endif

      implicit none
      integer :: ptsk,i,j
      logical :: PREDICTOR
      real(8) :: cff, cff1, cff2, cff3, cff4, cff5,tl_cff,tl_fac
      real(8) :: tl_zeta_new(LBi:UBi,LBj:UBj)
      real(8) :: tl_gzeta(LBi:UBi,LBj:UBj)
      real(8) :: tl_gzeta2(LBi:UBi,LBj:UBj)
      real(8) :: gzeta(LBi:UBi,LBj:UBj)
!      real(8) :: gzeta2(LBi:UBi,LBj:UBj)
      real(8) :: Drhs(LBi:UBi,LBj:UBj)
      real(8) :: Dstp(LBi:UBi,LBj:UBj)
      real(8) :: Dnew(LBi:UBi,LBj:UBj)
# ifdef UV_VIS2
      real(8) :: Drhs_p(LBi:UBi,LBj:UBj)
      real(8) :: tl_Drhs_p(LBi:UBi,LBj:UBj)
# endif
      real(8) :: tl_Drhs(LBi:UBi,LBj:UBj)
      real(8) :: tl_Dstp(LBi:UBi,LBj:UBj)
      real(8) :: tl_Dnew(LBi:UBi,LBj:UBj)

      real(8) :: tmp,tmp1,tmp2,tmp3,tmp4,tmp5,tl_tmp5

      character(len=80) :: fhead

! New:
# define IU_RANGE min(IstrU-1,IstrV),max(IendU+1,IendV+1)
# define JU_RANGE min(JstrU,JstrV-1),max(JendU,JendV)

!# define IV_RANGE min(IstrU-1,IstrV),max(IendU,IendV)
# define IV_RANGE IstrV-1,IendV+1
# define JV_RANGE min(JstrU,JstrV-1),max(JendU+1,JendV+1)

      ptsk=3-kstp
      
!-----------------------------------------------------------------------
!  Compute total depth (m) and vertically integrated mass fluxes.
!-----------------------------------------------------------------------
!
!      write (*,*) 'in step2d: Drhs ...'

#ifdef V31
# ifdef NS_PERIODIC
      DO j=Jstr,Jend
# else
      DO j=Jstr-1,Jend+1
# endif
       DO i=Istr-1,Iend+1

#else

! v3.2:
      DO j=LBj,UBj
       DO i=LBi,UBi

#endif
        Drhs(i,j)=zeta(i,j,krhs)+h(i,j)
        Dstp(i,j)=zeta(i,j,kstp)+h(i,j)
        Dnew(i,j)=zeta(i,j,knew)+h(i,j)

        tl_Drhs(i,j)=tl_zeta(i,j,krhs)
       END DO
      END DO

#ifdef NS_PERIODIC
      CALL exchange_2d(Drhs,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(Dstp,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(Dnew,LBi,UBi,LBj,UBj,Lm,Mm)
#endif

# ifdef UV_VIS2
!-----------------------------------------------------------------------
!  If horizontal mixing, compute total depth at PSI-points.
!-----------------------------------------------------------------------
      DO j=Jstr,Jend+1
       DO i=Istr,Iend+1
        Drhs_p(i,j)=0.25d0*(Drhs(i,j  )+Drhs(i-1,j  )+
     *                      Drhs(i,j-1)+Drhs(i-1,j-1))
       END DO
      END DO
# endif

#ifdef NS_PERIODIC
      CALL exchange_2d(tl_Drhs,LBi,UBi,LBj,UBj,Lm,Mm)
#endif

!      write (*,*) 'in step2d: vert integrated fluxes ...'

      DO j=JU_RANGE
       DO i=IU_RANGE
        DU(i,j)=0.5d0*(Drhs(i,j)+Drhs(i-1,j))*ubar(i,j,krhs)*dy_u(i,j)
        tl_DU(i,j)=0.5d0*dy_u(i,j)*
     *    ( (tl_zeta(i,j,krhs)+tl_zeta(i-1,j,krhs))*ubar(i,j,krhs) +
     *      (Drhs(i,j)+Drhs(i-1,j))*tl_ubar(i,j,krhs) )
       END DO
      END DO

! OLD: before adv is added
!      DO j=J_RANGE_DV
!       DO i=IstrV-1,IendV+1

      DO j=JV_RANGE
       DO i=IV_RANGE
        DV(i,j)=0.5d0*(Drhs(i,j)+Drhs(i,j-1))*vbar(i,j,krhs)*dx_v(i,j)
        tl_DV(i,j)=0.5d0*dx_v(i,j)*
     *    ( (tl_zeta(i,j,krhs)+tl_zeta(i,j-1,krhs))*vbar(i,j,krhs)+
     *      (Drhs(i,j)+Drhs(i,j-1))*tl_vbar(i,j,krhs) )
       END DO
      END DO

# ifdef NS_PERIODIC
      CALL exchange_2d(DU,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(DV,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(tl_DU,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(tl_DV,LBi,UBi,LBj,UBj,Lm,Mm)
# endif

!-----------------------------------------------------------------------
!  Compute time averaged fields over all short time-steps.
!-----------------------------------------------------------------------
!
!      write (*,*) 'in step2d: time-ave fields ...'

      IF (PREDICTOR) THEN
       IF (iif.eq.1) THEN
!
!  Reset arrays for 2D fields averaged within the short time-steps.
!
        cff2=(-1.d0/12.d0)*weight(2,iif+1)

! TL loop 1:
!        DO j=Jstr-1,Jend+1  ! ? roms ?
!         DO i=Istr-1,Iend+1

! Avrora v2.1:
!        DO j=Jstr,Jend
!         DO i=Istr,Iend

! Note (v3.1): new limits 
        DO j=LBj,UBj
         DO i=LBi,UBi
          tl_Zt_avg1(i,j)=0.d0
         END DO
        END DO

! TL loop 2:
        DO j=JU_RANGE
         DO i=IU_RANGE
          tl_DU_avg1(i,j)=0.d0
          tl_DU_avg2(i,j)=cff2*tl_DU(i,j)
         END DO
        END DO

! TL loop 3:
        DO j=JV_RANGE
         DO i=IV_RANGE
          tl_DV_avg1(i,j)=0.d0
          tl_DV_avg2(i,j)=cff2*tl_DV(i,j)
         END DO
        END DO

       ELSE ! (PREDICTOR, iif>1)
!
!  Accumulate field averages of previous time-step after they are
!  computed in the previous corrector step, updated their boundaries,
!  and synchronized.
!
        cff1=weight(1,iif-1)
        cff2=(8.d0/12.d0)*weight(2,iif)-
     *       (1.d0/12.d0)*weight(2,iif+1)

! TL loop 1:
!        DO j=Jstr-1,Jend+1
!         DO i=Istr-1,Iend+1

!        DO j=Jstr,Jend  ! <- v2.1
!         DO i=Istr,Iend

! Note (v3.1): new limits 
        DO j=LBj,UBj
         DO i=LBi,UBi
          tl_Zt_avg1(i,j)=tl_Zt_avg1(i,j)+cff1*tl_zeta(i,j,krhs)
         END DO
        END DO

 ! TL loop 2:
        DO j=JU_RANGE
         DO i=IU_RANGE
          tl_DU_avg1(i,j)=tl_DU_avg1(i,j)+cff1*tl_DU(i,j)
          tl_DU_avg2(i,j)=tl_DU_avg2(i,j)+cff2*tl_DU(i,j)
         END DO
        END DO

! TL loop 3:
        DO j=JV_RANGE
         DO i=IV_RANGE
          tl_DV_avg1(i,j)=tl_DV_avg1(i,j)+cff1*tl_DV(i,j)
          tl_DV_avg2(i,j)=tl_DV_avg2(i,j)+cff2*tl_DV(i,j)
         END DO
        END DO

       END IF ! (iif.eq.1) or else

      ELSE ! .NOT.PREDICTOR: CORRECTOR STEP

       IF (iif.eq.1) THEN
        cff2=weight(2,iif)
       ELSE
        cff2=(5.d0/12.d0)*weight(2,iif)
       END IF

! TL loop 1:
       DO j=JU_RANGE
        DO i=IU_RANGE
         tl_DU_avg2(i,j)=tl_DU_avg2(i,j)+cff2*tl_DU(i,j)
        END DO
       END DO

! TL loop 2:
       DO j=JV_RANGE
        DO i=IV_RANGE
         tl_DV_avg2(i,j)=tl_DV_avg2(i,j)+cff2*tl_DV(i,j)
        END DO
       END DO

      END IF  ! choices PREDICTOR or ELSE


! Finalize averaging, compute new tl_Hz etc. 
!
      IF (PREDICTOR .AND. (iif.eq.nfast+1)) THEN

# if defined NS_PERIODIC
       CALL exchange_2d(tl_Zt_avg1,LBi,UBi,LBj,UBj,Lm,Mm)
       CALL exchange_2d(tl_DU_avg1,LBi,UBi,LBj,UBj,Lm,Mm)
       CALL exchange_2d(tl_DV_avg1,LBi,UBi,LBj,UBj,Lm,Mm)
       CALL exchange_2d(tl_DU_avg2,LBi,UBi,LBj,UBj,Lm,Mm)
       CALL exchange_2d(tl_DV_avg2,LBi,UBi,LBj,UBj,Lm,Mm)
# endif

       CALL tl_set_depth(LBi,UBi,LBj,UBj,N,Lm,Mm,
     *               z0_r,z0_w,oh,
     *               tl_Zt_avg1,tl_z_r,tl_z_w,tl_Hz,tl_Hz_u,tl_Hz_v)

      END IF ! PREDICTOR, iif=nfast+1
       
      IF (iif.eq.nfast+1) THEN
       RETURN
      ENDIF

!=======================================================================
!  Time step free-surface equation.
!=======================================================================
!
!  During the first time-step, the predictor step is Forward-Euler
!  and the corrector step is Backward-Euler. Otherwise, the predictor
!  step is Leap-frog and the corrector step is Adams-Moulton.
!
      IF (iif.eq.1) THEN

!v3.1:
!       DO j=Jstr,Jend
!        DO i=Istr,Iend

!v3.2
       DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
         tl_rhs_zeta(i,j)=
     *              tl_DU(i,j)-tl_DU(i+1,j)+tl_DV(i,j)-tl_DV(i,j+1)
         tl_zeta_new(i,j)=tl_zeta(i,j,kstp)+
     *                  dtfast*oArea(i,j)*tl_rhs_zeta(i,j)
         tl_Dnew(i,j)=tl_zeta_new(i,j)  

! Note local definition of gzeta(i,j), TL only (consistent with
! notion that background zeta=const during b/t stepping), 
! background  gzeta2 is not needed 
         gzeta(i,j)=zeta(i,j,1)
! NL:         gzeta2(i,j)=gzeta(i,j)*gzeta(i,j)
         tl_gzeta(i,j)=0.5d0*(tl_zeta(i,j,kstp)+tl_zeta_new(i,j))
         tl_gzeta2(i,j)=2.d0*tl_gzeta(i,j)*gzeta(i,j)
        END DO
       END DO

      ELSE IF (PREDICTOR) THEN           ! time-step zeta, RPED, iif>1

       cff1=2.d0*dtfast
       cff4=0.16d0                ! 4.0_r8/25.0_r8
       cff5=1.d0-2.d0*cff4

!v3.1:
!       DO j=Jstr,Jend
!        DO i=Istr,Iend

!v3.2
       DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
         tl_rhs_zeta(i,j)=
     *     tl_DU(i,j)-tl_DU(i+1,j)+tl_DV(i,j)-tl_DV(i,j+1)
         tl_zeta_new(i,j)=tl_zeta(i,j,kstp)+
     *                 cff1*oArea(i,j)*tl_rhs_zeta(i,j)
         tl_Dnew(i,j)=tl_zeta_new(i,j) 

! see comment above about gzeta
         gzeta(i,j)=zeta(i,j,1)
! NL:         gzeta2(i,j)=gzeta(i,j)*gzeta(i,j)
         tl_gzeta(i,j)=cff5*tl_zeta(i,j,krhs)+
     *                 cff4*(tl_zeta(i,j,kstp)+tl_zeta_new(i,j))
         tl_gzeta2(i,j)=2.d0*tl_gzeta(i,j)*gzeta(i,j)
        END DO
       END DO

      ELSE IF (.NOT.PREDICTOR) THEN   ! time-step zeta, CORR, iif>1

       cff1=dtfast*5.d0/12.d0
       cff2=dtfast*8.d0/12.d0
       cff3=dtfast*1.d0/12.d0
       cff4=0.4d0                  ! 2.d0/5.d0
       cff5=1.d0-cff4

!v3.1:
!       DO j=Jstr,Jend
!        DO i=Istr,Iend

!v3.2
       DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
         tl_cff=cff1*(tl_DU(i,j)-tl_DU(i+1,j)+tl_DV(i,j)-tl_DV(i,j+1))
         tl_zeta_new(i,j)=tl_zeta(i,j,kstp)+
     *     oArea(i,j)*( tl_cff
     *                 +cff2*tl_rzeta(i,j,kstp)
     *                 -cff3*tl_rzeta(i,j,ptsk))
         tl_Dnew(i,j)=tl_zeta_new(i,j)

! see comment above about gzeta
         gzeta(i,j)=zeta(i,j,1)
! NL:         gzeta2(i,j)=gzeta(i,j)*gzeta(i,j)
         tl_gzeta(i,j)=cff5*tl_zeta_new(i,j)+cff4*tl_zeta(i,j,krhs)
         tl_gzeta2(i,j)=2.d0*tl_gzeta(i,j)*gzeta(i,j)
        END DO
       END DO

      END IF

!
!  Load new free-surface values into shared array at both predictor
!  and corrector steps.
!
      DO j=Jstr,Jend
       DO i=Istr,Iend
        tl_zeta(i,j,knew)=tl_zeta_new(i,j)
       END DO
      END DO
!
!  If predictor step, load right-side-term into shared array.
!
!      write (*,*) 'in step2d: rzeta ...'
      IF (PREDICTOR) THEN

! v3.1:
!       DO j=Jstr,Jend
!        DO i=Istr,Iend
! v.3.2:
       DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
         tl_rzeta(i,j,krhs)=tl_rhs_zeta(i,j)
        END DO
       END DO
      END IF
!
!  Set free-surface lateral boundary conditions.
!
!      write (*,*) 'in step2d: zetabc ...'
      CALL tl_zetabc (Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj,  
     *             krhs,kstp,knew,tl_zeta)

# if defined NS_PERIODIC
!      write (*,*) 'in step2d: exchange ...'
      CALL exchange_2d (tl_zeta(:,:,knew),LBi,UBi,LBj,UBj,Lm,Mm) 
      CALL exchange_2d (tl_Dnew,LBi,UBi,LBj,UBj,Lm,Mm) 
      CALL exchange_2d (gzeta,LBi,UBi,LBj,UBj,Lm,Mm) 
! NL only:      CALL exchange_2d (gzeta2,LBi,UBi,LBj,UBj,Lm,Mm) 
      CALL exchange_2d (tl_gzeta,LBi,UBi,LBj,UBj,Lm,Mm) 
      CALL exchange_2d (tl_gzeta2,LBi,UBi,LBj,UBj,Lm,Mm) 
# endif

!=======================================================================
!  Compute right-hand-side for the 2D momentum equations.
!=======================================================================
!
!      write (*,*) 'in step2d: rhs of mom eq. ...'
!-----------------------------------------------------------------------
!  Compute pressure gradient terms.
!-----------------------------------------------------------------------
!
      cff1=0.5d0*g

      DO j=JstrU,JendU
       DO i=IstrU,IendU
! NL:        rhs_ubar(i,j)=cff1*dy_u(i,j)*(
! NL:     *        (h(i-1,j)+h(i,j))*(gzeta(i-1,j)-gzeta(i,j))+
! NL:     *        gzeta2(i-1,j)-gzeta2(i,j)  )
        tl_rhs_ubar(i,j)=cff1*dy_u(i,j)*(
     *     (h(i-1,j)+h(i,j))*(tl_gzeta(i-1,j)-tl_gzeta(i,j))+
     *     tl_gzeta2(i-1,j)-tl_gzeta2(i,j) )

       END DO
      END DO

      DO j=JstrV,JendV
       DO i=IstrV,IendV
! NL:        rhs_vbar(i,j)=cff1*dx_v(i,j)*(
! NL:     *        (h(i,j-1)+h(i,j))*(gzeta(i,j-1)-gzeta(i,j))+
! NL:     *        gzeta2(i,j-1)-gzeta2(i,j)  )
        tl_rhs_vbar(i,j)=cff1*dx_v(i,j)*(
     *        (h(i,j-1)+h(i,j))*(tl_gzeta(i,j-1)-tl_gzeta(i,j))+
     *        tl_gzeta2(i,j-1)-tl_gzeta2(i,j)  )
       END DO
      END DO

!-----------------------------------------------------------------------
!  Add in Coriolis term
!-----------------------------------------------------------------------
!
# if defined UV_COR || (defined UV_ADV && defined CURVGRID)

! TL loop 1:

#ifdef V31
#  if defined NS_PERIODIC
      DO j=Jstr-1,Jend
#  else
      DO j=Jstr,Jend
#  endif
       DO i=Istr,Iend
#else

      DO j=Jstr-1,Jend
       DO i=Istr-1,Iend

#endif /* V31 */

        IF (mask_r(i,j).ne.0) THEN
 
         tmp5=(
#  ifdef UV_COR
     *        fdxdy(i,j)
#  endif
#  if (defined CURVGRID && defined UV_ADV)
     *       +0.5d0*((vbar(i,j  ,krhs)+vbar(i,j+1,krhs))*dndx(i,j)-
     *                (ubar(i  ,j,krhs)+ubar(i+1,j,krhs))*dmde(i,j))
#  endif
     *       )

#  if (defined CURVGRID && defined UV_ADV)
         tl_tmp5=0.5d0*
     *           ((tl_vbar(i,j  ,krhs)+tl_vbar(i,j+1,krhs))*dndx(i,j)-
     *            (tl_ubar(i  ,j,krhs)+tl_ubar(i+1,j,krhs))*dmde(i,j))         
#  endif

         cff=0.5d0*Drhs(i,j)*tmp5

         tl_cff=0.5d0*(
     *                  tl_Drhs(i,j)*tmp5
#  if (defined CURVGRID && defined UV_ADV)
     *                 +Drhs(i,j)*tl_tmp5
#  endif
     *                )


! NL:         UFx(i,j)=cff*(vbar(i,j,krhs)+vbar(i,j+1,krhs))
         tl_UFx(i,j)=tl_cff*(vbar(i,j,krhs)+vbar(i,j+1,krhs))+
     *                 cff*(tl_vbar(i,j,krhs)+tl_vbar(i,j+1,krhs))               

! NL:         VFe(i,j)=cff*(ubar(i,j,krhs)+ubar(i+1,j,krhs))
         tl_VFe(i,j)=tl_cff*(ubar(i,j,krhs)+ubar(i+1,j,krhs))+
     *               cff*(tl_ubar(i,j,krhs)+tl_ubar(i+1,j,krhs))
        END IF
       END DO
      END DO

! TL loop 2:
      DO j=JstrU,JendU
       DO i=IstrU,IendU
        tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+
     *                     0.5d0*(tl_UFx(i,j)+tl_UFx(i-1,j))
       END DO
      END DO
      
! TL loop 3:
      DO j=JstrV,JendV
       DO i=IstrV,IendV
        tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)-
     *                     0.5d0*(tl_VFe(i,j)+tl_VFe(i,j-1))
       END DO
      END DO

# endif /* defined UV_COR || (defined UV_ADV && defined CURVGRID) */

# ifdef UV_VIS2
!-----------------------------------------------------------------------
!  If horizontal mixing, compute total depth at PSI-points.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend+1
       DO i=Istr,Iend+1
        tl_Drhs_p(i,j)=0.25d0*(tl_Drhs(i,j  )+tl_Drhs(i-1,j  )+
     *                         tl_Drhs(i,j-1)+tl_Drhs(i-1,j-1))
       END DO
      END DO

!-----------------------------------------------------------------------
!  Add in horizontal harmonic viscosity.
!-----------------------------------------------------------------------

! ATTN: NEW COMPARED TO v2.1 
! zero out flux variables before they are used 
!(important to provide BC for the diss term)
!(later, investigation will be necessary to optimize this -- only boundary
!values need to be zeroed out).
      DO j=LBj,UBj
      DO i=LBi,UBi
       tl_UFx(i,j)=0.d0
       tl_UFe(i,j)=0.d0
       tl_VFx(i,j)=0.d0
       tl_VFe(i,j)=0.d0
      END DO
      END DO

!
!  Compute flux-components of the horizontal divergence of the stress
!  tensor (m5/s2) in XI- and ETA-directions.
!
      DO j=JstrV-1,JendV
       DO i=max(1,IstrV-1),IendV    ! new in v3.2, v3.1 had IstrV,IendV

!          cff=visc2_r(i,j)*Drhs(i,j)*0.5d0*
!     *           (dyodx(i,j)*
!     *            ((ody(i  ,j)+ody(i+1,j))*ubar(i+1,j,krhs)-                 
!     *             (ody(i-1,j)+ody(i  ,j))*ubar(i  ,j,krhs))-                
!     *            dxody(i,j)*                                            
!     *            ((odx(i,j  )+odx(i,j+1))*vbar(i,j+1,krhs)-                 
!     *             (odx(i,j-1)+odx(i,j  ))*vbar(i,j  ,krhs)))

        tmp=visc2_r(i,j)*0.5d0*dyodx(i,j)
        tmp1=tmp*(ody(i  ,j)+ody(i+1,j))        
        tmp2=tmp*(ody(i-1,j)+ody(i  ,j))

        tmp=visc2_r(i,j)*0.5d0*dxody(i,j)
        tmp3=tmp*(odx(i,j  )+odx(i,j+1))
        tmp4=tmp*(odx(i,j-1)+odx(i,j  ))

! Rewrite NL using defs above:
!         cff=Drhs(i,j)*( tmp1*ubar(i+1,j,krhs)-tmp2*ubar(i  ,j,krhs)
!     *                  -tmp3*vbar(i,j+1,krhs)+tmp4*vbar(i,j  ,krhs) )

          tl_cff=tl_Drhs(i,j)*
     *            ( tmp1*ubar(i+1,j,krhs)-tmp2*ubar(i  ,j,krhs)
     *             -tmp3*vbar(i,j+1,krhs)+tmp4*vbar(i,j  ,krhs) ) +
     *           Drhs(i,j)*
     *            ( tmp1*tl_ubar(i+1,j,krhs)-tmp2*tl_ubar(i  ,j,krhs)
     *             -tmp3*tl_vbar(i,j+1,krhs)+tmp4*tl_vbar(i,j  ,krhs) )

        tl_UFx(i,j)=dy(i,j)*dy(i,j)*tl_cff
        tl_VFe(i,j)=dx(i,j)*dx(i,j)*tl_cff
       END DO
      END DO

      DO j=JstrV,JendV+1
       DO i=IstrV,IendV+1

! NL:            cff=visc2_p(i,j)*Drhs_p(i,j)*       
! NL:     *          (dyodx_p(i,j)*                                     
! NL:     *           ((ody(i  ,j-1)+ody(i  ,j))*vbar(i  ,j,krhs)-  
! NL:     *            (ody(i-1,j-1)+ody(i-1,j))*vbar(i-1,j,krhs))+ 
! NL:     *           dxody_p(i,j)*                  
! NL:     *           ((odx(i-1,j  )+odx(i,j  ))*ubar(i,j  ,krhs)-
! NL:     *            (odx(i-1,j-1)+odx(i,j-1))*ubar(i,j-1,krhs)))

        tmp=visc2_p(i,j)*dyodx_p(i,j)
        tmp1=tmp*(ody(i  ,j-1)+ody(i  ,j))
        tmp2=tmp*(ody(i-1,j-1)+ody(i-1,j))

        tmp=visc2_p(i,j)*dxody_p(i,j)
        tmp3=tmp*(odx(i-1,j  )+odx(i,j  ))
        tmp4=tmp*(odx(i-1,j-1)+odx(i,j-1))

! Rewrite NL using defs above:
!           cff=Drhs_p(i,j)*( tmp1*vbar(i  ,j,krhs)-tmp2*vbar(i-1,j,krhs)
!     *                      +tmp3*ubar(i,j  ,krhs)-tmp4*ubar(i,j-1,krhs) )

        tl_cff= tl_Drhs_p(i,j)*
     *           ( tmp1*vbar(i  ,j,krhs)-tmp2*vbar(i-1,j,krhs)
     *            +tmp3*ubar(i,j  ,krhs)-tmp4*ubar(i,j-1,krhs) ) 
     *         +Drhs_p(i,j)*
     *           ( tmp1*tl_vbar(i  ,j,krhs)-tmp2*tl_vbar(i-1,j,krhs)
     *            +tmp3*tl_ubar(i,j  ,krhs)-tmp4*tl_ubar(i,j-1,krhs) ) 

        tl_cff=tl_cff*pmask(i,j)

        tl_UFe(i,j)=dx_p(i,j)*dx_p(i,j)*tl_cff
        tl_VFx(i,j)=dy_p(i,j)*dy_p(i,j)*tl_cff
       END DO
      END DO

!  Add in harmonic viscosity.
!
      DO j=JstrU,JendU
       DO i=IstrU,IendU
        tl_fac=0.5d0*((ody(i-1,j)+ody(i,j))* 
     *             (tl_UFx(i,j  )-tl_UFx(i-1,j))+
     *             (odx(i-1,j)+odx(i,j))*
     *             (tl_UFe(i,j+1)-tl_UFe(i  ,j)))
        tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+tl_fac
       END DO
      END DO

      DO j=JstrV,JendV
       DO i=IstrV,IendV
        tl_fac=0.5d0*((ody(i,j-1)+ody(i,j))*     
     *             (tl_VFx(i+1,j)-tl_VFx(i,j  ))-
     *             (odx(i,j-1)+odx(i,j))*   
     *             (tl_VFe(i  ,j)-tl_VFe(i,j-1)))
        tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+tl_fac
       END DO
      END DO

# endif /* UV_VIS2 */

# ifdef UV_ADV
!-------------------------------------------------------------
! Add in advection:
!-------------------------------------------------------------

!
! Second-order, centered differences:
!

! U-flux of u momentum (at rho locations)
      DO j=JstrU,JendU
       DO i=IstrU-1,IendU
! NL:        UFx(i,j)=0.25d0*(ubar(i,j,krhs)+ubar(i+1,j,krhs))*
! NL:     *                  (DU(i,j)+DU(i+1,j))
        tl_UFx(i,j)=0.25d0*(
     *    (tl_ubar(i,j,krhs)+tl_ubar(i+1,j,krhs))*(DU(i,j)+DU(i+1,j))+
     *    (ubar(i,j,krhs)+ubar(i+1,j,krhs))*(tl_DU(i,j)+tl_DU(i+1,j))
     *                     )
       END DO
      END DO

! V-flux of u momentum (at psi locations):

      DO j=JstrU,JendU+1
       DO i=IstrU,IendU
! NL:        UFe(i,j)=0.25d0*(ubar(i,j-1,krhs)+ubar(i,j,krhs))* 
! NL:     *                  (DV(i-1,j)+DV(i,j))
        tl_UFe(i,j)=0.25d0*(
     *    (tl_ubar(i,j-1,krhs)+tl_ubar(i,j,krhs))*(DV(i-1,j)+DV(i,j))+
     *    (ubar(i,j-1,krhs)+ubar(i,j,krhs))*(tl_DV(i-1,j)+tl_DV(i,j))
     *                     )
       END DO
      END DO

! U-flux of v momentum (at psi locations):
      DO j=JstrV,JendV
       DO i=IstrV,IendV+1
! NL:        VFx(i,j)=0.25d0*(vbar(i-1,j,krhs)+vbar(i,j,krhs))*
! NL:     *                  (DU(i,j-1)+DU(i,j))
        tl_VFx(i,j)=0.25d0*(
     *    (tl_vbar(i-1,j,krhs)+tl_vbar(i,j,krhs))*(DU(i,j-1)+DU(i,j))+
     *    (vbar(i-1,j,krhs)+vbar(i,j,krhs))*(tl_DU(i,j-1)+tl_DU(i,j))
     *                     )
       END DO
      END DO

! V-flux of v momentum (at rho locations):
      DO j=JstrV-1,JendV
       DO i=IstrV,IendV
! NL:        VFe(i,j)=0.25d0*(vbar(i,j,krhs)+vbar(i,j+1,krhs))*
! NL:     *                  (DV(i,j)+DV(i,j+1))
        tl_VFe(i,j)=0.25d0*(
     *    (tl_vbar(i,j,krhs)+tl_vbar(i,j+1,krhs))*(DV(i,j)+DV(i,j+1))+
     *    (vbar(i,j,krhs)+vbar(i,j+1,krhs))*(tl_DV(i,j)+tl_DV(i,j+1))
     *                  )
       END DO
      END DO


!
! Add advection to rhs_ubar, rhs_vbar
!
      DO j=JstrU,JendU
       DO i=IstrU,IendU
! NL:        fac=UFx(i,j)-UFx(i-1,j)+UFe(i,j+1)-UFe(i,j)
        tl_fac=tl_UFx(i,j)-tl_UFx(i-1,j)+tl_UFe(i,j+1)-tl_UFe(i,j)
        tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)-tl_fac
       END DO
      END DO

      DO j=JstrV,JendV
       DO i=IstrV,IendV
! NL:        fac=VFx(i+1,j)-VFx(i,j)+VFe(i,j)-VFe(i,j-1)
        tl_fac=tl_VFx(i+1,j)-tl_VFx(i,j)+tl_VFe(i,j)-tl_VFe(i,j-1)
        tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)-tl_fac
       END DO
      END DO

# endif /* UV_ADV */

!      write (*,*) 'in step2d: rhs of mom eq.2d/3d coupling ...'
!-----------------------------------------------------------------------
!  Coupling between 2D and 3D equations.
!-----------------------------------------------------------------------
! 
!  NOTE FROM ROMS:
! 
!  Before the predictor step of the first barotropic time-step,
!  arrays "rufrc" and "rvfrc" contain the vertical integrals of
!  the 3D right-hand-side terms for momentum equations (including
!  surface and bottom stresses, if so prescribed).
!
!  Convert them into forcing terms by subtracting the fast time
!  "rhs_ubar" and "rhs_vbar" from them; Also, immediately apply
!  these forcing terms "rhs_ubar" and "rhs_vbar".
!
!  From now on, these newly computed forcing terms will remain
!  constant during the fast time stepping and will added to
!  "rhs_ubar" and "rhs_vbar" during all subsequent time steps.
!
      IF (iif.eq.1 .and. PREDICTOR) THEN

       IF (iic.eq.1) THEN

        DO j=JstrU,JendU
         DO i=IstrU,IendU
          tl_rufrc(i,j)=tl_rufrc(i,j)-tl_rhs_ubar(i,j)
          tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+tl_rufrc(i,j)
          tl_ru(i,j,0,nstp)=tl_rufrc(i,j)
         END DO
        END DO

        DO j=JstrV,JendV
         DO i=IstrV,IendV
          tl_rvfrc(i,j)=tl_rvfrc(i,j)-tl_rhs_vbar(i,j)
          tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+tl_rvfrc(i,j)
          tl_rv(i,j,0,nstp)=tl_rvfrc(i,j)
         END DO
        END DO

       ELSE IF (iic.eq.2) THEN

        DO j=JstrU,JendU
         DO i=IstrU,IendU
          tl_rufrc(i,j)=tl_rufrc(i,j)-tl_rhs_ubar(i,j)
          tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+
     *                   1.5d0*tl_rufrc(i,j)-0.5d0*tl_ru(i,j,0,nnew)
          tl_ru(i,j,0,nstp)=tl_rufrc(i,j)
         END DO
        END DO

        DO j=JstrV,JendV
         DO i=IstrV,IendV
          tl_rvfrc(i,j)=tl_rvfrc(i,j)-tl_rhs_vbar(i,j)
          tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+
     *                   1.5d0*tl_rvfrc(i,j)-0.5d0*tl_rv(i,j,0,nnew)
          tl_rv(i,j,0,nstp)=tl_rvfrc(i,j)
         END DO
        END DO

       ELSE   ! iic>2

        cff1=23.d0/12.d0
        cff2=16.d0/12.d0
        cff3= 5.d0/12.d0

        DO j=JstrU,JendU
         DO i=IstrU,IendU
          tl_rufrc(i,j)=tl_rufrc(i,j)-tl_rhs_ubar(i,j)
          tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+
     *                     cff1*tl_rufrc(i,j)-
     *                     cff2*tl_ru(i,j,0,nnew)+
     *                     cff3*tl_ru(i,j,0,nstp)
          tl_ru(i,j,0,nstp)=tl_rufrc(i,j)
         END DO
        END DO
        DO j=JstrV,JendV
         DO i=IstrV,IendV
          tl_rvfrc(i,j)=tl_rvfrc(i,j)-tl_rhs_vbar(i,j)
          tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+
     *                     cff1*tl_rvfrc(i,j)-
     *                     cff2*tl_rv(i,j,0,nnew)+
     *                     cff3*tl_rv(i,j,0,nstp)
          tl_rv(i,j,0,nstp)=tl_rvfrc(i,j)
         END DO
        END DO

       END IF  ! IF (iic==1)... ELSEIF (iic==2)... ELSE...

      ELSE  ! iif==1, CORR step, or (iif>1, both PRED and CORR)
 
       DO j=JstrU,JendU
        DO i=IstrU,IendU
         tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+tl_rufrc(i,j)
        END DO
       END DO
       DO j=JstrV,JendV
        DO i=IstrV,IendV
         tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+tl_rvfrc(i,j)
        END DO
       END DO

      END IF  ! (iif.eq.1 .and. PREDICTOR) ... ELSE...

!      write (*,*) 'in step2d: time-step mom eq. ...'
!=======================================================================
!  Time step 2D momentum equations.
!=======================================================================
!  Compute total water column depth.
!

#ifdef MPI
      DO j=LBj,UBj
       DO i=LBi,UBi
#else
      DO j=Jstr,Jend
       DO i=Istr,Iend
#endif /* MPI or serial */

        tl_Dstp(i,j)=tl_zeta(i,j,kstp)

       END DO
      END DO

#ifdef NS_PERIODIC
      CALL exchange_2d(tl_Dstp,LBi,UBi,LBj,UBj,Lm,Mm)
#endif

!
!  During the first time-step, the predictor step is Forward-Euler
!  and the corrector step is Backward-Euler. Otherwise, the predictor
!  step is Leap-frog and the corrector step is Adams-Moulton.
!

      IF (iif.eq.1) THEN

       cff1=2.d0*dtfast

       DO j=JstrU,JendU
        DO i=IstrU,IendU
         IF (mask_u(i,j).ne.0) THEN 
! NL:          ubar(i,j,knew)=( ubar(i,j,kstp)*(Dstp(i,j)+Dstp(i-1,j))
! NL:     *                    +oArea_u(i,j)*cff1*rhs_ubar(i,j) )
! NL:     *                   /(Dnew(i,j)+Dnew(i-1,j))
          cff=1.d0/(Dnew(i,j)+Dnew(i-1,j))

          tl_ubar(i,j,knew)=cff*( 
     *      tl_ubar(i,j,kstp)*(Dstp(i,j)+Dstp(i-1,j))+
     *      ubar(i,j,kstp)*(tl_Dstp(i,j)+tl_Dstp(i-1,j))+
     *      oArea_u(i,j)*cff1*tl_rhs_ubar(i,j)-
     *      (tl_Dnew(i,j)+tl_Dnew(i-1,j))*ubar(i,j,knew) )

         END IF 
        END DO
       END DO

       DO j=JstrV,JendV
        DO i=IstrV,IendV
         IF (mask_v(i,j).ne.0) THEN 
! NL:          vbar(i,j,knew)=( vbar(i,j,kstp)*(Dstp(i,j)+Dstp(i,j-1))+
! NL:     *                     oArea_v(i,j)*cff1*rhs_vbar(i,j) )
! NL:     *                   /(Dnew(i,j)+Dnew(i,j-1))
          cff=1.d0/(Dnew(i,j)+Dnew(i,j-1))
          tl_vbar(i,j,knew)=cff*(
     *      tl_vbar(i,j,kstp)*(Dstp(i,j)+Dstp(i,j-1))+
     *      vbar(i,j,kstp)*(tl_Dstp(i,j)+tl_Dstp(i,j-1))+
     *      oArea_v(i,j)*cff1*tl_rhs_vbar(i,j)-
     *      (tl_Dnew(i,j)+tl_Dnew(i,j-1))*vbar(i,j,knew) )
         END IF 
        END DO
       END DO

      ELSE IF (PREDICTOR) THEN
       cff1=4.d0*dtfast  ! step for 2*dtfast, also factor 2 since 2D used

       DO j=JstrU,JendU
        DO i=IstrU,IendU
         IF (mask_u(i,j).ne.0) THEN 
! NL:          ubar(i,j,knew)=( ubar(i,j,kstp)*(Dstp(i,j)+Dstp(i-1,j))
! NL:     *                    +oArea_u(i,j)*cff1*rhs_ubar(i,j) )
! NL:     *                   /(Dnew(i,j)+Dnew(i-1,j))
          cff=1.d0/(Dnew(i,j)+Dnew(i-1,j))
          tl_ubar(i,j,knew)=cff*( 
     *      tl_ubar(i,j,kstp)*(Dstp(i,j)+Dstp(i-1,j))+
     *      ubar(i,j,kstp)*(tl_Dstp(i,j)+tl_Dstp(i-1,j))+
     *      oArea_u(i,j)*cff1*tl_rhs_ubar(i,j)-
     *      (tl_Dnew(i,j)+tl_Dnew(i-1,j))*ubar(i,j,knew) )    
         END IF 
        END DO
       END DO

       DO j=JstrV,JendV
        DO i=IstrV,IendV
         IF (mask_v(i,j).ne.0) THEN 
! NL:          vbar(i,j,knew)=( vbar(i,j,kstp)*(Dstp(i,j)+Dstp(i,j-1))+
! NL:     *                     oArea_v(i,j)*cff1*rhs_vbar(i,j) )
! NL:     *                   /(Dnew(i,j)+Dnew(i,j-1))
          cff=1.d0/(Dnew(i,j)+Dnew(i,j-1))
          tl_vbar(i,j,knew)=cff*(
     *      tl_vbar(i,j,kstp)*(Dstp(i,j)+Dstp(i,j-1))+
     *      vbar(i,j,kstp)*(tl_Dstp(i,j)+tl_Dstp(i,j-1))+
     *      oArea_v(i,j)*cff1*tl_rhs_vbar(i,j)-
     *      (tl_Dnew(i,j)+tl_Dnew(i,j-1))*vbar(i,j,knew) )
         END IF  MPI
        END DO
       END DO

      ELSE IF (.NOT.PREDICTOR) THEN
       cff1=2.d0*dtfast*5.d0/12.d0 ! Use factor 2 since 
       cff2=2.d0*dtfast*8.d0/12.d0 ! twice D is used below, both
       cff3=2.d0*dtfast*1.d0/12.d0 ! in numerator and denominator

       DO j=JstrU,JendU
        DO i=IstrU,IendU
         IF (mask_u(i,j).ne.0) THEN 
! NL:          ubar(i,j,knew)=( ubar(i,j,kstp)*(Dstp(i,j)+Dstp(i-1,j)) 
! NL:     *                       +oArea_u(i,j)*(cff1*rhs_ubar(i,j)+
! NL:     *                                     cff2*rubar(i,j,kstp)-
! NL:     *                                     cff3*rubar(i,j,ptsk)) )
! NL:     *                     /(Dnew(i,j)+Dnew(i-1,j))
          cff=1.d0/(Dnew(i,j)+Dnew(i-1,j))
          tl_ubar(i,j,knew)=cff*( 
     *      tl_ubar(i,j,kstp)*(Dstp(i,j)+Dstp(i-1,j))+
     *      ubar(i,j,kstp)*(tl_Dstp(i,j)+tl_Dstp(i-1,j))+
     *      oArea_u(i,j)*( cff1*tl_rhs_ubar(i,j)+
     *                     cff2*tl_rubar(i,j,kstp)-
     *                     cff3*tl_rubar(i,j,ptsk) )-
     *      (tl_Dnew(i,j)+tl_Dnew(i-1,j))*ubar(i,j,knew) )
         END IF  MPI
        END DO
       END DO
        
       DO j=JstrV,JendV
        DO i=IstrV,IendV
         IF (mask_v(i,j).ne.0) THEN 
! NL:          vbar(i,j,knew)=( vbar(i,j,kstp)*(Dstp(i,j)+Dstp(i,j-1)) 
! NL:     *                       +oArea_v(i,j)*(cff1*rhs_vbar(i,j)+
! NL:     *                                     cff2*rvbar(i,j,kstp)-
! NL:     *                                     cff3*rvbar(i,j,ptsk)) )
! NL:     *                     /(Dnew(i,j)+Dnew(i,j-1))
          cff=1.d0/(Dnew(i,j)+Dnew(i,j-1))
          tl_vbar(i,j,knew)=cff*(
     *      tl_vbar(i,j,kstp)*(Dstp(i,j)+Dstp(i,j-1))   +
     *      vbar(i,j,kstp)*(tl_Dstp(i,j)+tl_Dstp(i,j-1))+
     *      oArea_v(i,j)*( cff1*tl_rhs_vbar(i,j)  +
     *                     cff2*tl_rvbar(i,j,kstp)-
     *                     cff3*tl_rvbar(i,j,ptsk) )    - 
     *      (tl_Dnew(i,j)+tl_Dnew(i,j-1))*vbar(i,j,knew) )
         END IF 
        END DO
       END DO
      END IF   ! (iif==1)... ELSEIF (PRED)... ELSEIF (.NOT.PRED)...
!
!  If predictor step, load right-side-term into shared arrays for
!  future use during the subsequent corrector step.
!

      IF (PREDICTOR) THEN
       DO j=JstrU,JendU
        DO i=IstrU,IendU
         tl_rubar(i,j,krhs)=tl_rhs_ubar(i,j)
        END DO
       END DO
       DO j=JstrV,JendV
        DO i=IstrV,IendV
         tl_rvbar(i,j,krhs)=tl_rhs_vbar(i,j)
        END DO
       END DO
      END IF

!-----------------------------------------------------------------------
!  Apply lateral boundary conditions.

!      write (*,*) 'in step2d: lateral BC ...'      
      CALL tl_u2dbc (Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj,  
     *               krhs, kstp, knew,
     *               iif,PREDICTOR,g,h,mask_u,
#ifdef MPI
     *               B_West,B_East,B_South,B_North,
#endif
     *               tl_ubar,tl_zeta)

      CALL tl_v2dbc (Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj,  
     *               krhs, kstp, knew,
     *               iif,PREDICTOR,g,h,mask_v,
#ifdef MPI
     *               B_West,B_East,B_South,B_North,
#endif
     *               tl_vbar,tl_zeta)

#ifdef MPI

      CALL mpi_exchange_2d_2(tl_zeta(LBi:UBi,LBj:UBj,knew),
     *                     LBi,UBi,LBj,UBj,
     *                     Istr_E,Iend_E,Jstr_E,Jend_E,
     *                     comm2d, stride, 
     *                     rankleft,rankright,ranktop,rankbottom) 

      CALL mpi_exchange_2d_2(tl_ubar(LBi:UBi,LBj:UBj,knew),
     *                     LBi,UBi,LBj,UBj,
     *                     IstrU_E,IendU_E,JstrU_E,JendU_E,
     *                     comm2d, stride, 
     *                     rankleft,rankright,ranktop,rankbottom) 

      CALL mpi_exchange_2d_2(tl_vbar(LBi:UBi,LBj:UBj,knew),
     *                     LBi,UBi,LBj,UBj,
     *                     IstrV_E,IendV_E,JstrV_E,JendV_E,
     *                     comm2d, stride, 
     *                     rankleft,rankright,ranktop,rankbottom) 

#endif

# ifdef NS_PERIODIC
!      write (*,*) 'in step2d: periodic  BC ...'      
      CALL exchange_2d(tl_ubar(:,:,knew),LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(tl_vbar(:,:,knew),LBi,UBi,LBj,UBj,Lm,Mm)
# endif

      RETURN
      END
#endif /* TANGENT */
