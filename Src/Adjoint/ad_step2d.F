#include "cppdefs.h"

#ifdef ADJOINT
      SUBROUTINE ad_step2d(PREDICTOR)
      USE mod_ocean
#ifdef MPI
      USE mod_parallel
#endif

      implicit none
      integer :: ptsk,i,j
      logical :: PREDICTOR
      real(8) :: cff, cff1, cff2, cff3, cff4, cff5,ad_cff,ad_fac
      real(8) :: tmp,tmp1,tmp2,tmp3,tmp4,tmp5,ad_tmp5
      real(8) :: ad_zeta_new(LBi:UBi,LBj:UBj)
      real(8) :: ad_gzeta(LBi:UBi,LBj:UBj)
      real(8) :: ad_gzeta2(LBi:UBi,LBj:UBj)
      real(8) :: gzeta(LBi:UBi,LBj:UBj)
      real(8) :: Drhs(LBi:UBi,LBj:UBj)
      real(8) :: Dstp(LBi:UBi,LBj:UBj)
      real(8) :: Dnew(LBi:UBi,LBj:UBj)
# ifdef UV_VIS2
      real(8) :: Drhs_p(LBi:UBi,LBj:UBj)
      real(8) :: ad_Drhs_p(LBi:UBi,LBj:UBj)
# endif
      real(8) :: ad_Drhs(LBi:UBi,LBj:UBj)
      real(8) :: ad_Dstp(LBi:UBi,LBj:UBj)
      real(8) :: ad_Dnew(LBi:UBi,LBj:UBj)

#ifdef MPI
      CALL mpi_exchange_2d_ext(ad_zeta(LBi:UBi,LBj:UBj,knew),
     *                     LBi,UBi,LBj,UBj,
     *                     LBi0,UBi0,LBj0,UBj0,
     *                     Istr_E,Iend_E,Jstr_E,Jend_E,
     *                     comm2d, stride, 
     *                     rankleft,rankright,ranktop,rankbottom) 

      CALL mpi_exchange_2d_ext(ad_ubar(LBi:UBi,LBj:UBj,knew),
     *                     LBi,UBi,LBj,UBj,
     *                     LBi0,UBi0,LBj0,UBj0,
     *                     IstrU_E,IendU_E,JstrU_E,JendU_E,
     *                     comm2d, stride, 
     *                     rankleft,rankright,ranktop,rankbottom) 

      CALL mpi_exchange_2d_ext(ad_vbar(LBi:UBi,LBj:UBj,knew),
     *                     LBi,UBi,LBj,UBj,
     *                     LBi0,UBi0,LBj0,UBj0,
     *                     IstrV_E,IendV_E,JstrV_E,JendV_E,
     *                     comm2d, stride, 
     *                     rankleft,rankright,ranktop,rankbottom)

#endif

      ad_zeta_new=0.d0
      ad_gzeta=0.d0
      ad_gzeta2=0.d0
      ad_Drhs=0.d0
      ad_Dstp=0.d0
      ad_Dnew=0.d0
# ifdef UV_VIS2
      ad_Drhs_p=0.d0
# endif
      ad_cff=0.d0
      ad_fac=0.d0

      ad_tmp5=0.d0

! New:
# define IU_RANGE min(IstrU-1,IstrV),max(IendU+1,IendV+1)
# define JU_RANGE min(JstrU,JstrV-1),max(JendU,JendV)

!# define IV_RANGE min(IstrU-1,IstrV),max(IendU,IendV)
# define IV_RANGE IstrV-1,IendV+1
# define JV_RANGE min(JstrU,JstrV-1),max(JendU+1,JendV+1)

      ptsk=3-kstp
!-----------------------------------------------------------------------
!  Compute total depth (m): Background, all time levels
!-----------------------------------------------------------------------
#ifdef V31

# ifdef NS_PERIODIC
      DO j=Jstr,Jend
# else
      DO j=Jstr-1,Jend+1
# endif
       DO i=Istr-1,Iend+1

#else /* not V31 */

! v3.2:
      DO j=LBj,UBj
       DO i=LBi,UBi

#endif
        Drhs(i,j)=zeta(i,j,krhs)+h(i,j)
        Dstp(i,j)=zeta(i,j,kstp)+h(i,j)
        Dnew(i,j)=zeta(i,j,knew)+h(i,j)
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
# endif /* UV_VIS2 */

      DO j=JU_RANGE
       DO i=IU_RANGE
        DU(i,j)=0.5d0*(Drhs(i,j)+Drhs(i-1,j))*ubar(i,j,krhs)*dy_u(i,j)
       END DO
      END DO

      DO j=JV_RANGE
       DO i=IV_RANGE
        DV(i,j)=0.5d0*(Drhs(i,j)+Drhs(i,j-1))*vbar(i,j,krhs)*dx_v(i,j)
       END DO
      END DO

# ifdef NS_PERIODIC
      CALL exchange_2d(DU,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(DV,LBi,UBi,LBj,UBj,Lm,Mm)
# endif

!==================================================================
! After defining fwd vars, begin w/ ADJ, from the bottom of TL code
!==================================================================

!**** The most part of the code is performed only for iif< nfast+1
      IF (iif.ne.nfast+1) THEN    


# ifdef NS_PERIODIC
      CALL ad_exchange_2d(ad_vbar(:,:,knew),LBi,UBi,LBj,UBj,Lm,Mm)
      CALL ad_exchange_2d(ad_ubar(:,:,knew),LBi,UBi,LBj,UBj,Lm,Mm)
# endif

!-----------------------------------------------------------------------
!  ADJ to lateral boundary conditions.
!-----------------------------------------------------------------------
      CALL ad_v2dbc (Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj,  
     *               krhs, kstp, knew,
     *               iif,PREDICTOR,g,h,mask_v,
#ifdef MPI
     *               B_West,B_East,B_South,B_North, ! <- new line, comp ot v2.1
#endif
     *               ad_vbar,ad_zeta)

      CALL ad_u2dbc (Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj,  
     *               krhs, kstp, knew,
     *               iif,PREDICTOR,g,h,mask_u,
#ifdef MPI
     *               B_West,B_East,B_South,B_North, ! <- new line, comp ot v2.1
#endif
     *               ad_ubar,ad_zeta)

!-----------------------------------------------------------------------
!  ADJ TO:
!  If predictor step, load right-side-term into shared arrays for
!  future use during the subsequent corrector step.
!-----------------------------------------------------------------------
!
      IF (PREDICTOR) THEN

       DO j=JstrV,JendV
        DO i=IstrV,IendV
! TL:         tl_rvbar(i,j,krhs)=tl_rhs_vbar(i,j)
         ad_rhs_vbar(i,j)=ad_rhs_vbar(i,j)+ad_rvbar(i,j,krhs)
         ad_rvbar(i,j,krhs)=0.d0
        END DO
       END DO

       DO j=JstrU,JendU
        DO i=IstrU,IendU
! TL:         tl_rubar(i,j,krhs)=tl_rhs_ubar(i,j)
         ad_rhs_ubar(i,j)=ad_rhs_ubar(i,j)+ad_rubar(i,j,krhs)
         ad_rubar(i,j,krhs)=0.d0
        END DO
       END DO

      END IF

!=======================================================================
!  ADJ to time stepping 2D momentum equations.
!=======================================================================

      IF (iif.eq.1) THEN

       cff1=2.d0*dtfast

       DO j=JstrV,JendV
        DO i=IstrV,IendV

         IF (mask_v(i,j).ne.0) THEN 

! NL:          vbar(i,j,knew)=( vbar(i,j,kstp)*(Dstp(i,j)+Dstp(i,j-1))+
! NL:     *                     oArea_v(i,j)*cff1*rhs_vbar(i,j) )
! NL:     *                   /(Dnew(i,j)+Dnew(i,j-1))

          cff=1.d0/(Dnew(i,j)+Dnew(i,j-1))

! TL:          tl_vbar(i,j,knew)=cff*(
! TL:     *      tl_vbar(i,j,kstp)*(Dstp(i,j)+Dstp(i,j-1))+
! TL:     *      vbar(i,j,kstp)*(tl_Dstp(i,j)+tl_Dstp(i,j-1))+
! TL:     *      oArea_v(i,j)*cff1*tl_rhs_vbar(i,j)-
! TL:     *      (tl_Dnew(i,j)+tl_Dnew(i,j-1))*vbar(i,j,knew) )

          tmp=cff*ad_vbar(i,j,knew)
          tmp1=tmp*vbar(i,j,kstp)
          tmp2=tmp*vbar(i,j,knew)

          ad_vbar(i,j,kstp)=ad_vbar(i,j,kstp)+
     *                       tmp*(Dstp(i,j)+Dstp(i,j-1))
          ad_Dstp(i,j)=ad_Dstp(i,j)+tmp1
          ad_Dstp(i,j-1)=ad_Dstp(i,j-1)+tmp1
          ad_rhs_vbar(i,j)=ad_rhs_vbar(i,j)+tmp*oArea_v(i,j)*cff1
          ad_Dnew(i,j)=ad_Dnew(i,j)-tmp2
          ad_Dnew(i,j-1)=ad_Dnew(i,j-1)-tmp2

          ad_vbar(i,j,knew)=0.d0

         END IF 

        END DO
       END DO

       DO j=JstrU,JendU
        DO i=IstrU,IendU
         IF (mask_u(i,j).ne.0) THEN 
! NL:          ubar(i,j,knew)=( ubar(i,j,kstp)*(Dstp(i,j)+Dstp(i-1,j))
! NL:     *                    +oArea_u(i,j)*cff1*rhs_ubar(i,j) )
! NL:     *                   /(Dnew(i,j)+Dnew(i-1,j))

          cff=1.d0/(Dnew(i,j)+Dnew(i-1,j))

! TL:          tl_ubar(i,j,knew)=cff*( 
! TL:     *      tl_ubar(i,j,kstp)*(Dstp(i,j)+Dstp(i-1,j))+
! TL:     *      ubar(i,j,kstp)*(tl_Dstp(i,j)+tl_Dstp(i-1,j))+
! TL:     *      oArea_u(i,j)*cff1*tl_rhs_ubar(i,j)-
! TL:     *      (tl_Dnew(i,j)+tl_Dnew(i-1,j))*ubar(i,j,knew) )

          tmp=ad_ubar(i,j,knew)*cff
          tmp1=tmp*ubar(i,j,kstp)
          tmp2=tmp*ubar(i,j,knew)

          ad_ubar(i,j,kstp)=ad_ubar(i,j,kstp)
     *                      +tmp*(Dstp(i,j)+Dstp(i-1,j))
          ad_Dstp(i,j)=ad_Dstp(i,j)+tmp1
          ad_Dstp(i-1,j)=ad_Dstp(i-1,j)+tmp1
          ad_rhs_ubar(i,j)=ad_rhs_ubar(i,j)+tmp*oArea_u(i,j)*cff1
          ad_Dnew(i,j)=ad_Dnew(i,j)-tmp2
          ad_Dnew(i-1,j)=ad_Dnew(i-1,j)-tmp2

          ad_ubar(i,j,knew)=0.d0

         END IF 
        END DO
       END DO

      ELSE IF (PREDICTOR) THEN
       cff1=4.d0*dtfast  ! step for 2*dtfast, also factor 2 since 2D used

       DO j=JstrV,JendV
        DO i=IstrV,IendV
         IF (mask_v(i,j).ne.0) THEN 
! NL:          vbar(i,j,knew)=( vbar(i,j,kstp)*(Dstp(i,j)+Dstp(i,j-1))+
! NL:     *                     oArea_v(i,j)*cff1*rhs_vbar(i,j) )
! NL:     *                   /(Dnew(i,j)+Dnew(i,j-1))

          cff=1.d0/(Dnew(i,j)+Dnew(i,j-1))

! TL:          tl_vbar(i,j,knew)=cff*(
! TL:     *      tl_vbar(i,j,kstp)*(Dstp(i,j)+Dstp(i,j-1))+
! TL:     *      vbar(i,j,kstp)*(tl_Dstp(i,j)+tl_Dstp(i,j-1))+
! TL:     *      oArea_v(i,j)*cff1*tl_rhs_vbar(i,j)-
! TL:     *      (tl_Dnew(i,j)+tl_Dnew(i,j-1))*vbar(i,j,knew) )

          tmp=ad_vbar(i,j,knew)*cff
          tmp1=tmp*vbar(i,j,kstp)
          tmp2=tmp*vbar(i,j,knew)

          ad_vbar(i,j,kstp)=ad_vbar(i,j,kstp)
     *                          +tmp*(Dstp(i,j)+Dstp(i,j-1))
          ad_Dstp(i,j)=ad_Dstp(i,j)+tmp1
          ad_Dstp(i,j-1)=ad_Dstp(i,j-1)+tmp1
          ad_rhs_vbar(i,j)=ad_rhs_vbar(i,j)+tmp*oArea_v(i,j)*cff1
          ad_Dnew(i,j)=ad_Dnew(i,j)-tmp2
          ad_Dnew(i,j-1)=ad_Dnew(i,j-1)-tmp2

          ad_vbar(i,j,knew)=0.d0

         END IF 
        END DO
       END DO

       DO j=JstrU,JendU
        DO i=IstrU,IendU
         IF (mask_u(i,j).ne.0) THEN 
! NL:          ubar(i,j,knew)=( ubar(i,j,kstp)*(Dstp(i,j)+Dstp(i-1,j))
! NL:     *                    +oArea_u(i,j)*cff1*rhs_ubar(i,j) )
! NL:     *                   /(Dnew(i,j)+Dnew(i-1,j))

          cff=1.d0/(Dnew(i,j)+Dnew(i-1,j))

! TL:          tl_ubar(i,j,knew)=cff*( 
! TL:     *      tl_ubar(i,j,kstp)*(Dstp(i,j)+Dstp(i-1,j))+
! TL:     *      ubar(i,j,kstp)*(tl_Dstp(i,j)+tl_Dstp(i-1,j))+
! TL:     *      oArea_u(i,j)*cff1*tl_rhs_ubar(i,j)-
! TL:     *      (tl_Dnew(i,j)+tl_Dnew(i-1,j))*ubar(i,j,knew) )    


          tmp=ad_ubar(i,j,knew)*cff
          tmp1=tmp*ubar(i,j,kstp)
          tmp2=tmp*ubar(i,j,knew)

          ad_ubar(i,j,kstp)=ad_ubar(i,j,kstp)
     *                      +tmp*(Dstp(i,j)+Dstp(i-1,j))
          ad_Dstp(i,j)=ad_Dstp(i,j)+tmp1
          ad_Dstp(i-1,j)=ad_Dstp(i-1,j)+tmp1
          ad_rhs_ubar(i,j)=ad_rhs_ubar(i,j)+tmp*oArea_u(i,j)*cff1
          ad_Dnew(i,j)=ad_Dnew(i,j)-tmp2
          ad_Dnew(i-1,j)=ad_Dnew(i-1,j)-tmp2

          ad_ubar(i,j,knew)=0.d0

         END IF 
        END DO
       END DO

      ELSE IF (.NOT.PREDICTOR) THEN
       cff1=2.d0*dtfast*5.d0/12.d0 ! Use factor 2 since 
       cff2=2.d0*dtfast*8.d0/12.d0 ! twice D is used below, both
       cff3=2.d0*dtfast*1.d0/12.d0 ! in numerator and denominator

       DO j=JstrV,JendV
        DO i=IstrV,IendV
         IF (mask_v(i,j).ne.0) THEN 
! NL:          vbar(i,j,knew)=( vbar(i,j,kstp)*(Dstp(i,j)+Dstp(i,j-1)) 
! NL:     *                       +oArea_v(i,j)*(cff1*rhs_vbar(i,j)+
! NL:     *                                     cff2*rvbar(i,j,kstp)-
! NL:     *                                     cff3*rvbar(i,j,ptsk)) )
! NL:     *                     /(Dnew(i,j)+Dnew(i,j-1))

          cff=1.d0/(Dnew(i,j)+Dnew(i,j-1))

! TL:          tl_vbar(i,j,knew)=cff*(
! TL:     *      tl_vbar(i,j,kstp)*(Dstp(i,j)+Dstp(i,j-1))   +
! TL:     *      vbar(i,j,kstp)*(tl_Dstp(i,j)+tl_Dstp(i,j-1))+
! TL:     *      oArea_v(i,j)*( cff1*tl_rhs_vbar(i,j)  +
! TL:     *                     cff2*tl_rvbar(i,j,kstp)-
! TL:     *                     cff3*tl_rvbar(i,j,ptsk) )    - 
! TL:     *      (tl_Dnew(i,j)+tl_Dnew(i,j-1))*vbar(i,j,knew) )
          
          tmp=cff*ad_vbar(i,j,knew)
          tmp1=tmp*vbar(i,j,kstp)
          tmp2=tmp*vbar(i,j,knew)
          tmp3=tmp*oArea_v(i,j)
          
          ad_vbar(i,j,kstp)=ad_vbar(i,j,kstp)
     *                       +tmp*(Dstp(i,j)+Dstp(i,j-1))
          ad_Dstp(i,j)=ad_Dstp(i,j)+tmp1
          ad_Dstp(i,j-1)=ad_Dstp(i,j-1)+tmp1
          ad_rhs_vbar(i,j)=ad_rhs_vbar(i,j)+cff1*tmp3
          ad_rvbar(i,j,kstp)=ad_rvbar(i,j,kstp)+cff2*tmp3
          ad_rvbar(i,j,ptsk)=ad_rvbar(i,j,ptsk)-cff3*tmp3
          ad_Dnew(i,j)=ad_Dnew(i,j)-tmp2
          ad_Dnew(i,j-1)=ad_Dnew(i,j-1)-tmp2

          ad_vbar(i,j,knew)=0.d0

         END IF 
        END DO
       END DO

       DO j=JstrU,JendU
        DO i=IstrU,IendU
         IF (mask_u(i,j).ne.0) THEN 
! NL:          ubar(i,j,knew)=( ubar(i,j,kstp)*(Dstp(i,j)+Dstp(i-1,j)) 
! NL:     *                       +oArea_u(i,j)*(cff1*rhs_ubar(i,j)+
! NL:     *                                     cff2*rubar(i,j,kstp)-
! NL:     *                                     cff3*rubar(i,j,ptsk)) )
! NL:     *                     /(Dnew(i,j)+Dnew(i-1,j))

          cff=1.d0/(Dnew(i,j)+Dnew(i-1,j))

! TL:          tl_ubar(i,j,knew)=cff*( 
! TL:     *      tl_ubar(i,j,kstp)*(Dstp(i,j)+Dstp(i-1,j))+
! TL:     *      ubar(i,j,kstp)*(tl_Dstp(i,j)+tl_Dstp(i-1,j))+
! TL:     *      oArea_u(i,j)*( cff1*tl_rhs_ubar(i,j)+
! TL:     *                     cff2*tl_rubar(i,j,kstp)-
! TL:     *                     cff3*tl_rubar(i,j,ptsk) )-
! TL:     *      (tl_Dnew(i,j)+tl_Dnew(i-1,j))*ubar(i,j,knew) )

          tmp=cff*ad_ubar(i,j,knew)
          tmp1=tmp*ubar(i,j,kstp)
          tmp2=tmp*ubar(i,j,knew)
          tmp3=tmp*oArea_u(i,j)

          ad_ubar(i,j,kstp)=ad_ubar(i,j,kstp)
     *                         +tmp*(Dstp(i,j)+Dstp(i-1,j))
          ad_Dstp(i,j)=ad_Dstp(i,j)+tmp1
          ad_Dstp(i-1,j)=ad_Dstp(i-1,j)+tmp1
          ad_rhs_ubar(i,j)=ad_rhs_ubar(i,j)+tmp3*cff1
          ad_rubar(i,j,kstp)=ad_rubar(i,j,kstp)+tmp3*cff2
          ad_rubar(i,j,ptsk)=ad_rubar(i,j,ptsk)-tmp3*cff3
          ad_Dnew(i,j)=ad_Dnew(i,j)-tmp2
          ad_Dnew(i-1,j)=ad_Dnew(i-1,j)-tmp2

          ad_ubar(i,j,knew)=0.d0

         END IF 
        END DO
       END DO
        
      END IF   ! (iif==1)... ELSEIF (PRED)... ELSEIF (.NOT.PRED)...


#ifdef NS_PERIODIC
      CALL ad_exchange_2d(ad_Dstp,LBi,UBi,LBj,UBj,Lm,Mm)
#endif

!  ADj to: Compute total water column depth (tl_Dstp)

      DO j=Jstr,Jend
       DO i=Istr,Iend
! NL (computed above):        Dstp(i,j)=zeta(i,j,kstp)+h(i,j)
! TL:        tl_Dstp(i,j)=tl_zeta(i,j,kstp)

        ad_zeta(i,j,kstp)=ad_zeta(i,j,kstp)+ad_Dstp(i,j)
        ad_Dstp(i,j)=0.d0

       END DO
      END DO

!-----------------------------------------------------------------------
!  ADJ to: Coupling between 2D and 3D equations.
!-----------------------------------------------------------------------
      
      IF (iif.eq.1 .and. PREDICTOR) THEN

       IF (iic.eq.1) THEN

        DO j=JstrV,JendV
         DO i=IstrV,IendV

! TL line 3:          tl_rv(i,j,0,nstp)=tl_rvfrc(i,j)
          ad_rvfrc(i,j)=ad_rvfrc(i,j)+ad_rv(i,j,0,nstp)
          ad_rv(i,j,0,nstp)=0.d0

! TL line 2 (recurs.): tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+tl_rvfrc(i,j)
          ad_rvfrc(i,j)=ad_rvfrc(i,j)+ad_rhs_vbar(i,j)

! TL line 1 (recurs.): tl_rvfrc(i,j)=tl_rvfrc(i,j)-tl_rhs_vbar(i,j)
          ad_rhs_vbar(i,j)=ad_rhs_vbar(i,j)-ad_rvfrc(i,j)

         END DO
        END DO

        DO j=JstrU,JendU
         DO i=IstrU,IendU
! TL line 3:           tl_ru(i,j,0,nstp)=tl_rufrc(i,j)
          ad_rufrc(i,j)=ad_rufrc(i,j)+ad_ru(i,j,0,nstp)
          ad_ru(i,j,0,nstp)=0.d0

! TL line 2 (recurs):           tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+tl_rufrc(i,j)
          ad_rufrc(i,j)=ad_rufrc(i,j)+ad_rhs_ubar(i,j)

! TL line 1 (recurs):          tl_rufrc(i,j)=tl_rufrc(i,j)-tl_rhs_ubar(i,j)
          ad_rhs_ubar(i,j)=ad_rhs_ubar(i,j)-ad_rufrc(i,j)

         END DO
        END DO

       ELSE IF (iic.eq.2) THEN

        DO j=JstrV,JendV
         DO i=IstrV,IendV
! TL line 3:          tl_rv(i,j,0,nstp)=tl_rvfrc(i,j)
          ad_rvfrc(i,j)=ad_rvfrc(i,j)+ad_rv(i,j,0,nstp)
          ad_rv(i,j,0,nstp)=0.d0

! TL line 2 (recurs):         tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+
! TL line 2:     *                   1.5d0*tl_rvfrc(i,j)-0.5d0*tl_rv(i,j,0,nnew)
          ad_rvfrc(i,j)=ad_rvfrc(i,j)+1.5d0*ad_rhs_vbar(i,j)
          ad_rv(i,j,0,nnew)=ad_rv(i,j,0,nnew)-0.5d0*ad_rhs_vbar(i,j)
 
! TL line 1 (recurs):          tl_rvfrc(i,j)=tl_rvfrc(i,j)-tl_rhs_vbar(i,j)
          ad_rhs_vbar(i,j)=ad_rhs_vbar(i,j)-ad_rvfrc(i,j)

         END DO
        END DO

        DO j=JstrU,JendU
         DO i=IstrU,IendU
! TL line 3:          tl_ru(i,j,0,nstp)=tl_rufrc(i,j)
          ad_rufrc(i,j)=ad_rufrc(i,j)+ad_ru(i,j,0,nstp)
          ad_ru(i,j,0,nstp)=0.d0

! TL line 2 (recurs):  tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+
! TL line 2:     *       1.5d0*tl_rufrc(i,j)-0.5d0*tl_ru(i,j,0,nnew)
          ad_rufrc(i,j)=ad_rufrc(i,j)+1.5d0*ad_rhs_ubar(i,j)
          ad_ru(i,j,0,nnew)=ad_ru(i,j,0,nnew)-0.5d0*ad_rhs_ubar(i,j)
 
! TL line 1 (recurs):  tl_rufrc(i,j)=tl_rufrc(i,j)-tl_rhs_ubar(i,j)
          ad_rhs_ubar(i,j)=ad_rhs_ubar(i,j)-ad_rufrc(i,j)

         END DO
        END DO

       ELSE   ! iic>2

        cff1=23.d0/12.d0
        cff2=16.d0/12.d0
        cff3= 5.d0/12.d0

        DO j=JstrV,JendV
         DO i=IstrV,IendV
! TL line 3:          tl_rv(i,j,0,nstp)=tl_rvfrc(i,j)
          ad_rvfrc(i,j)=ad_rvfrc(i,j)+ad_rv(i,j,0,nstp)
          ad_rv(i,j,0,nstp)=0.d0

! TL line 2 (recurs):          tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+
! TL line 2     *                     cff1*tl_rvfrc(i,j)-
! TL line 2     *                     cff2*tl_rv(i,j,0,nnew)+
! TL line 2     *                     cff3*tl_rv(i,j,0,nstp)
          ad_rvfrc(i,j)=ad_rvfrc(i,j)+ad_rhs_vbar(i,j)*cff1
          ad_rv(i,j,0,nnew)=ad_rv(i,j,0,nnew)-ad_rhs_vbar(i,j)*cff2
          ad_rv(i,j,0,nstp)=ad_rv(i,j,0,nstp)+ad_rhs_vbar(i,j)*cff3

! TL line 1 (recurs):          tl_rvfrc(i,j)=tl_rvfrc(i,j)-tl_rhs_vbar(i,j)
          ad_rhs_vbar(i,j)=ad_rhs_vbar(i,j)-ad_rvfrc(i,j)

         END DO
        END DO

        DO j=JstrU,JendU
         DO i=IstrU,IendU
! TL line 3:          tl_ru(i,j,0,nstp)=tl_rufrc(i,j)
          ad_rufrc(i,j)=ad_rufrc(i,j)+ad_ru(i,j,0,nstp)
          ad_ru(i,j,0,nstp)=0.d0
          

! TL line 2 (recurs):          tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+
! TL line 2:     *                     cff1*tl_rufrc(i,j)-
! TL line 2:     *                     cff2*tl_ru(i,j,0,nnew)+
! TL line 2:     *                     cff3*tl_ru(i,j,0,nstp)
          ad_rufrc(i,j)=ad_rufrc(i,j)+cff1*ad_rhs_ubar(i,j)
          ad_ru(i,j,0,nnew)=ad_ru(i,j,0,nnew)-cff2*ad_rhs_ubar(i,j)
          ad_ru(i,j,0,nstp)=ad_ru(i,j,0,nstp)+cff3*ad_rhs_ubar(i,j)

! TL line 1 (recurs):          tl_rufrc(i,j)=tl_rufrc(i,j)-tl_rhs_ubar(i,j)
          ad_rhs_ubar(i,j)=ad_rhs_ubar(i,j)-ad_rufrc(i,j)

         END DO
        END DO

       END IF  ! IF (iic==1)... ELSEIF (iic==2)... ELSE...

      ELSE  ! iif==1, CORR step, or (iif>1, both PRED and CORR)

       DO j=JstrV,JendV
        DO i=IstrV,IendV
! TL (recurs):         tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+tl_rvfrc(i,j)
         ad_rvfrc(i,j)=ad_rvfrc(i,j)+ad_rhs_vbar(i,j)
        END DO
       END DO

       DO j=JstrU,JendU
        DO i=IstrU,IendU
! TL (recurs):         tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+tl_rufrc(i,j)
         ad_rufrc(i,j)=ad_rufrc(i,j)+ad_rhs_ubar(i,j)
        END DO
       END DO

      END IF  ! (iif.eq.1 .and. PREDICTOR) ... ELSE...

!=======================================================================
!  ADJ to: Compute right-hand-side for the 2D momentum equations.
!=======================================================================

# ifdef UV_ADV
!-------------------------------------------------------------
! Add in advection:
!-------------------------------------------------------------

!
! ADj to: Add advection to rhs_ubar, rhs_vbar
!

      DO j=JstrV,JendV
       DO i=IstrV,IendV

! tl line 2:     tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)-tl_fac
        ad_fac=ad_fac-ad_rhs_vbar(i,j)

! tl line 1:     tl_fac=tl_VFx(i+1,j)-tl_VFx(i,j)+tl_VFe(i,j)-tl_VFe(i,j-1)
        ad_VFx(i+1,j)=ad_VFx(i+1,j)+ad_fac
        ad_VFx(i,j)=ad_VFx(i,j)-ad_fac
        ad_VFe(i,j)=ad_VFe(i,j)+ad_fac
        ad_VFe(i,j-1)=ad_VFe(i,j-1)-ad_fac

        ad_fac=0.d0

       END DO
      END DO

      DO j=JstrU,JendU
       DO i=IstrU,IendU

! tl line 2:     tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)-tl_fac
        ad_fac=ad_fac-ad_rhs_ubar(i,j)

! tl line 1:     tl_fac=tl_UFx(i,j)-tl_UFx(i-1,j)+tl_UFe(i,j+1)-tl_UFe(i,j)
        ad_UFx(i,j)=ad_UFx(i,j)+ad_fac
        ad_UFx(i-1,j)=ad_UFx(i-1,j)-ad_fac
        ad_UFe(i,j+1)=ad_UFe(i,j+1)+ad_fac
        ad_UFe(i,j)=ad_UFe(i,j)-ad_fac

        ad_fac=0.d0

       END DO
      END DO

!
! Second-order, centered differences:
!

! V-flux of v momentum (at rho locations):
      DO j=JstrV-1,JendV
       DO i=IstrV,IendV
! NL:        VFe(i,j)=0.25d0*(vbar(i,j,krhs)+vbar(i,j+1,krhs))*
! NL:     *                  (DV(i,j)+DV(i,j+1))

! TL:        tl_VFe(i,j)=0.25d0*(
! TL:     *   (tl_vbar(i,j,krhs)+tl_vbar(i,j+1,krhs))*(DV(i,j)+DV(i,j+1))+
! TL:     *   (vbar(i,j,krhs)+vbar(i,j+1,krhs))*(tl_DV(i,j)+tl_DV(i,j+1))
! TL:     *                  )
        tmp=0.25d0*ad_VFe(i,j)
        tmp1=tmp*(DV(i,j)+DV(i,j+1))
        tmp2=tmp*(vbar(i,j,krhs)+vbar(i,j+1,krhs))

        ad_vbar(i,j,krhs)=ad_vbar(i,j,krhs)+tmp1
        ad_vbar(i,j+1,krhs)=ad_vbar(i,j+1,krhs)+tmp1
        ad_DV(i,j)=ad_DV(i,j)+tmp2
        ad_DV(i,j+1)=ad_DV(i,j+1)+tmp2

        ad_VFe(i,j)=0.d0

       END DO
      END DO

! U-flux of v momentum (at psi locations):
      DO j=JstrV,JendV
       DO i=IstrV,IendV+1
! NL:        VFx(i,j)=0.25d0*(vbar(i-1,j,krhs)+vbar(i,j,krhs))*
! NL:     *                  (DU(i,j-1)+DU(i,j))

! TL:        tl_VFx(i,j)=0.25d0*(
! TL:     *    (tl_vbar(i-1,j,krhs)+tl_vbar(i,j,krhs))*(DU(i,j-1)+DU(i,j))+
! TL:     *    (vbar(i-1,j,krhs)+vbar(i,j,krhs))*(tl_DU(i,j-1)+tl_DU(i,j))
! TL:     *                     )
        tmp=0.25d0*ad_VFx(i,j)
        tmp1=tmp*(DU(i,j-1)+DU(i,j))
        tmp2=tmp*(vbar(i-1,j,krhs)+vbar(i,j,krhs))

        ad_vbar(i-1,j,krhs)=ad_vbar(i-1,j,krhs)+tmp1
        ad_vbar(i,j,krhs)=ad_vbar(i,j,krhs)+tmp1
        ad_DU(i,j-1)=ad_DU(i,j-1)+tmp2
        ad_DU(i,j)=ad_DU(i,j)+tmp2

        ad_VFx(i,j)=0.d0

       END DO
      END DO

! V-flux of u momentum (at psi locations):

      DO j=JstrU,JendU+1
       DO i=IstrU,IendU
! NL:        UFe(i,j)=0.25d0*(ubar(i,j-1,krhs)+ubar(i,j,krhs))* 
! NL:     *                  (DV(i-1,j)+DV(i,j))

! TL:        tl_UFe(i,j)=0.25d0*(
! TL:     *    (tl_ubar(i,j-1,krhs)+tl_ubar(i,j,krhs))*(DV(i-1,j)+DV(i,j))+
! TL:     *    (ubar(i,j-1,krhs)+ubar(i,j,krhs))*(tl_DV(i-1,j)+tl_DV(i,j))
! TL:     *                     )

        tmp=0.25d0*ad_UFe(i,j)
        tmp1=tmp*(DV(i-1,j)+DV(i,j))
        tmp2=tmp*(ubar(i,j-1,krhs)+ubar(i,j,krhs))

        ad_ubar(i,j-1,krhs)=ad_ubar(i,j-1,krhs)+tmp1
        ad_ubar(i,j,krhs)=ad_ubar(i,j,krhs)+tmp1
        ad_DV(i-1,j)=ad_DV(i-1,j)+tmp2
        ad_DV(i,j)=ad_DV(i,j)+tmp2
 
        ad_UFe(i,j)=0.d0

       END DO
      END DO

! U-flux of u momentum (at rho locations)
      DO j=JstrU,JendU
       DO i=IstrU-1,IendU
! NL:        UFx(i,j)=0.25d0*(ubar(i,j,krhs)+ubar(i+1,j,krhs))*
! NL:     *                  (DU(i,j)+DU(i+1,j))

! TL:        tl_UFx(i,j)=0.25d0*(
! TL:     *    (tl_ubar(i,j,krhs)+tl_ubar(i+1,j,krhs))*(DU(i,j)+DU(i+1,j))+
! TL:     *    (ubar(i,j,krhs)+ubar(i+1,j,krhs))*(tl_DU(i,j)+tl_DU(i+1,j))
! TL:     *                     )
        tmp=0.25d0*ad_UFx(i,j)
        tmp1=tmp*(DU(i,j)+DU(i+1,j))
        tmp2=tmp*(ubar(i,j,krhs)+ubar(i+1,j,krhs))

        ad_ubar(i,j,krhs)=ad_ubar(i,j,krhs)+tmp1
        ad_ubar(i+1,j,krhs)=ad_ubar(i+1,j,krhs)+tmp1
        ad_DU(i,j)=ad_DU(i,j)+tmp2
        ad_DU(i+1,j)=ad_DU(i+1,j)+tmp2

        ad_UFx(i,j)=0.d0

       END DO
      END DO

# endif /* UV_ADV */

# ifdef UV_VIS2
!
!  ADj to: add in harmonic viscosity.
!
      DO j=JstrV,JendV
       DO i=IstrV,IendV
! tl line 2:        tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)+tl_fac
        ad_fac=ad_fac+ad_rhs_vbar(i,j)

! tl line 1:        tl_fac=0.5d0*((ody(i,j-1)+ody(i,j))*     
! tl line 1:     *             (tl_VFx(i+1,j)-tl_VFx(i,j  ))-
! tl line 1:     *             (odx(i,j-1)+odx(i,j))*   
! tl line 1:     *             (tl_VFe(i  ,j)-tl_VFe(i,j-1)))

        tmp1=0.5d0*(ody(i,j-1)+ody(i,j))*ad_fac
        tmp2=0.5d0*(odx(i,j-1)+odx(i,j))*ad_fac

        ad_VFx(i+1,j)=ad_VFx(i+1,j)+tmp1
        ad_VFx(i,j  )=ad_VFx(i,j  )-tmp1
        ad_VFe(i  ,j)=ad_VFe(i  ,j)-tmp2
        ad_VFe(i,j-1)=ad_VFe(i,j-1)+tmp2

        ad_fac=0.d0
        
       END DO
      END DO

      DO j=JstrU,JendU
       DO i=IstrU,IendU
! tl line 2:        tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+tl_fac
        ad_fac=ad_fac+ad_rhs_ubar(i,j)

! tl line 1:        tl_fac=0.5d0*((ody(i-1,j)+ody(i,j))* 
! tl line 1:     *             (tl_UFx(i,j  )-tl_UFx(i-1,j))+
! tl line 1:     *             (odx(i-1,j)+odx(i,j))*
! tl line 1:     *             (tl_UFe(i,j+1)-tl_UFe(i  ,j)))

        tmp1=0.5d0*(ody(i-1,j)+ody(i,j))*ad_fac
        tmp2=0.5d0*(odx(i-1,j)+odx(i,j))*ad_fac

        ad_UFx(i,j  )=ad_UFx(i,j  )+tmp1
        ad_UFx(i-1,j)=ad_UFx(i-1,j)-tmp1
        ad_UFe(i,j+1)=ad_UFe(i,j+1)+tmp2
        ad_UFe(i  ,j)=ad_UFe(i  ,j)-tmp2

        ad_fac=0.d0

       END DO
      END DO

!
!  Compute flux-components of the horizontal divergence of the stress
!  tensor (m5/s2) in XI- and ETA-directions.
!
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

! tl line 4:            tl_VFx(i,j)=dy_p(i,j)*dy_p(i,j)*tl_cff
        ad_cff=ad_cff+dy_p(i,j)*dy_p(i,j)*ad_VFx(i,j)
        ad_VFx(i,j)=0.d0

! tl line 3:            tl_UFe(i,j)=dx_p(i,j)*dx_p(i,j)*tl_cff
        ad_cff=ad_cff+dx_p(i,j)*dx_p(i,j)*ad_UFe(i,j)
        ad_UFe(i,j)=0.d0

! tl line 2:            tl_cff=tl_cff*pmask(i,j)
        ad_cff=ad_cff*pmask(i,j)

! tl line 1:            tl_cff= tl_Drhs_p(i,j)*
! tl line 1:     *             ( tmp1*vbar(i  ,j,krhs)-tmp2*vbar(i-1,j,krhs)
! tl line 1:     *              +tmp3*ubar(i,j  ,krhs)-tmp4*ubar(i,j-1,krhs) ) 
! tl line 1:     *           +Drhs_p(i,j)*
! tl line 1:     *             ( tmp1*tl_vbar(i  ,j,krhs)-tmp2*tl_vbar(i-1,j,krhs)
! tl line 1:     *              +tmp3*tl_ubar(i,j  ,krhs)-tmp4*tl_ubar(i,j-1,krhs) ) 

        ad_Drhs_p(i,j)=ad_Drhs_p(i,j)+
     *           ( tmp1*vbar(i  ,j,krhs)-tmp2*vbar(i-1,j,krhs)
     *            +tmp3*ubar(i,j  ,krhs)-tmp4*ubar(i,j-1,krhs) )*ad_cff 

        tmp=ad_cff*Drhs_p(i,j)

        ad_vbar(i  ,j,krhs)=ad_vbar(i  ,j,krhs)+tmp1*tmp
        ad_vbar(i-1,j,krhs)=ad_vbar(i-1,j,krhs)-tmp2*tmp  
        ad_ubar(i,j  ,krhs)=ad_ubar(i,j  ,krhs)+tmp3*tmp
        ad_ubar(i,j-1,krhs)=ad_ubar(i,j-1,krhs)-tmp4*tmp  

        ad_cff=0.d0

       END DO
      END DO

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

! tl line 3:        tl_VFe(i,j)=dx(i,j)*dx(i,j)*tl_cff
        ad_cff=ad_cff+dx(i,j)*dx(i,j)*ad_VFe(i,j)
        ad_VFe(i,j)=0.d0

! tl line 2:        tl_UFx(i,j)=dy(i,j)*dy(i,j)*tl_cff
        ad_cff=ad_cff+dy(i,j)*dy(i,j)*ad_UFx(i,j)
        ad_UFx(i,j)=0.d0

! tl line 1:          tl_cff=tl_Drhs(i,j)*
! tl line 1:     *            ( tmp1*ubar(i+1,j,krhs)-tmp2*ubar(i  ,j,krhs)
! tl line 1:     *             -tmp3*vbar(i,j+1,krhs)+tmp4*vbar(i,j  ,krhs) ) +
! tl line 1:     *           Drhs(i,j)*
! tl line 1:     *            ( tmp1*tl_ubar(i+1,j,krhs)-tmp2*tl_ubar(i  ,j,krhs)
! tl line 1:     *             -tmp3*tl_vbar(i,j+1,krhs)+tmp4*tl_vbar(i,j  ,krhs) )
        ad_Drhs(i,j)=ad_Drhs(i,j)+
     *          ad_cff*( tmp1*ubar(i+1,j,krhs)-tmp2*ubar(i  ,j,krhs)
     *                  -tmp3*vbar(i,j+1,krhs)+tmp4*vbar(i,j  ,krhs) )

        tmp=ad_cff*Drhs(i,j)

        ad_ubar(i+1,j,krhs)=ad_ubar(i+1,j,krhs)+tmp1*tmp
        ad_ubar(i  ,j,krhs)=ad_ubar(i  ,j,krhs)-tmp2*tmp
        ad_vbar(i,j+1,krhs)=ad_vbar(i,j+1,krhs)-tmp3*tmp
        ad_vbar(i,j  ,krhs)=ad_vbar(i,j  ,krhs)+tmp4*tmp
 
        ad_cff=0.d0
      
       END DO
      END DO

! ATTN: NEW COMPARED TO v2.1 
! ADJ TO: zero out flux variables before they are used 
!(important to provide BC for the diss term)
!(later, investigation will be necessary to optimize this -- only boundary
!values need to be zeroed out).

      DO j=LBj,UBj
      DO i=LBi,UBi
       ad_UFx(i,j)=0.d0
       ad_UFe(i,j)=0.d0
       ad_VFx(i,j)=0.d0
       ad_VFe(i,j)=0.d0
      END DO
      END DO
 

!-----------------------------------------------------------------------
!  If horizontal mixing, compute total depth at PSI-points.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend+1
       DO i=Istr,Iend+1
! TL:        tl_Drhs_p(i,j)=0.25d0*(tl_Drhs(i,j  )+tl_Drhs(i-1,j  )+
! TL:     *                         tl_Drhs(i,j-1)+tl_Drhs(i-1,j-1))

        tmp=0.25d0*ad_Drhs_p(i,j)

        ad_Drhs(i,j  )=ad_Drhs(i,j  )+tmp
        ad_Drhs(i-1,j  )=ad_Drhs(i-1,j  )+tmp
        ad_Drhs(i,j-1)=ad_Drhs(i,j-1)+tmp
        ad_Drhs(i-1,j-1)=ad_Drhs(i-1,j-1)+tmp

        ad_Drhs_p(i,j)=0.d0

       END DO
      END DO

# endif /* UV_VIS2 */


!-----------------------------------------------------------------------
!  ADJ to: Add in Coriolis term
!-----------------------------------------------------------------------

# if defined UV_COR || (defined UV_ADV && defined CURVGRID)

! TL loop 3:
      DO j=JstrV,JendV
       DO i=IstrV,IendV
! TL (recurs):        tl_rhs_vbar(i,j)=tl_rhs_vbar(i,j)-
! TL:     *                     0.5d0*(tl_VFe(i,j)+tl_VFe(i,j-1))
        tmp=0.5d0*ad_rhs_vbar(i,j)
        ad_VFe(i,j)=ad_VFe(i,j)-tmp
        ad_VFe(i,j-1)=ad_VFe(i,j-1)-tmp
       END DO
      END DO

! TL loop 2:
      DO j=JstrU,JendU
       DO i=IstrU,IendU
! TL (recurs):        tl_rhs_ubar(i,j)=tl_rhs_ubar(i,j)+
! TL:     *                     0.5d0*(tl_UFx(i,j)+tl_UFx(i-1,j))
        tmp=0.5d0*ad_rhs_ubar(i,j)
        ad_UFx(i,j)=ad_UFx(i,j)+tmp
        ad_UFx(i-1,j)=ad_UFx(i-1,j)+tmp
       END DO
      END DO
      
! TL loop 1:

#ifdef V31
#  if defined NS_PERIODIC
      DO j=Jstr-1,Jend
#  else
      DO j=Jstr,Jend
#  endif
       DO i=Istr,Iend
#else /* not V31 */

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

         cff=0.5d0*Drhs(i,j)*tmp5





! NL:         VFe(i,j)=cff*(ubar(i,j,krhs)+ubar(i+1,j,krhs))
! TL line 3:         tl_VFe(i,j)=tl_cff*(ubar(i,j,krhs)+ubar(i+1,j,krhs))+
! TL line 3     *                cff*(tl_ubar(i,j,krhs)+tl_ubar(i+1,j,krhs)) 

         tmp=ad_VFe(i,j)
         tmp1=cff*tmp

         ad_cff=ad_cff+(ubar(i,j,krhs)+ubar(i+1,j,krhs))*tmp
         ad_ubar(i,j,krhs)=ad_ubar(i,j,krhs)+tmp1
         ad_ubar(i+1,j,krhs)=ad_ubar(i+1,j,krhs)+tmp1
         
         ad_VFe(i,j)=0.d0
         
! NL:         UFx(i,j)=cff*(vbar(i,j,krhs)+vbar(i,j+1,krhs))*omask_vv(i,j)
! TL line 2:         tl_UFx(i,j)=tl_cff*(vbar(i,j,krhs)+vbar(i,j+1,krhs))+
! TL line 2:     *               cff*(tl_vbar(i,j,krhs)+tl_vbar(i,j+1,krhs)) 

         tmp=ad_UFx(i,j)
         tmp1=tmp*cff

         ad_cff=ad_cff+tmp*(vbar(i,j,krhs)+vbar(i,j+1,krhs))
         ad_vbar(i,j,krhs)=ad_vbar(i,j,krhs)+tmp1
         ad_vbar(i,j+1,krhs)=ad_vbar(i,j+1,krhs)+tmp1

         ad_UFx(i,j)=0.d0

! TL line 1:
!         tl_cff=0.5d0*(
!     *                  tl_Drhs(i,j)*tmp5
!#  if (defined CURVGRID && defined UV_ADV)
!     *                 +Drhs(i,j)*tl_tmp5
!#  endif
!     *                )

          ad_Drhs(i,j)=ad_Drhs(i,j)+0.5d0*tmp5*ad_cff
#  if (defined CURVGRID && defined UV_ADV)
          ad_tmp5=ad_tmp5+0.5d0*Drhs(i,j)*ad_cff
#  endif
         ad_cff=0.d0
 
! TL line, extra def, tl_tmp5:
#  if (defined CURVGRID && defined UV_ADV)
! TL:         tl_tmp5=0.5d0*
!     *           ((tl_vbar(i,j  ,krhs)+tl_vbar(i,j+1,krhs))*dndx(i,j)-
!     *            (tl_ubar(i  ,j,krhs)+tl_ubar(i+1,j,krhs))*dmde(i,j))      

         tmp3=0.5d0*ad_tmp5*dndx(i,j)
         tmp4=0.5d0*ad_tmp5*dmde(i,j)
 
         ad_vbar(i,j  ,krhs)=ad_vbar(i,j  ,krhs)+tmp3
         ad_vbar(i,j+1,krhs)=ad_vbar(i,j+1,krhs)+tmp3
         ad_ubar(i  ,j,krhs)=ad_ubar(i  ,j,krhs)-tmp4
         ad_ubar(i+1,j,krhs)=ad_ubar(i+1,j,krhs)-tmp4

         ad_tmp5=0.d0
#  endif

        END IF
       END DO

      END DO

# endif /* defined UV_COR || defined UV_COR_MASK */


!-----------------------------------------------------------------------
!  ADJ to: Compute pressure gradient terms.
!-----------------------------------------------------------------------
!
      cff1=0.5d0*g

      DO j=JstrV,JendV
       DO i=IstrV,IendV
! NL:        rhs_vbar(i,j)=cff1*dx_v(i,j)*(
! NL:     *        (h(i,j-1)+h(i,j))*(gzeta(i,j-1)-gzeta(i,j))+
! NL:     *        gzeta2(i,j-1)-gzeta2(i,j)  )

! TL:        tl_rhs_vbar(i,j)=cff1*dx_v(i,j)*(
! TL:     *        (h(i,j-1)+h(i,j))*(tl_gzeta(i,j-1)-tl_gzeta(i,j))+
! TL:     *        tl_gzeta2(i,j-1)-tl_gzeta2(i,j)  )

        tmp=ad_rhs_vbar(i,j)*cff1*dx_v(i,j)
        tmp1=(h(i,j-1)+h(i,j))*tmp

        ad_gzeta(i,j-1)=ad_gzeta(i,j-1)+tmp1
        ad_gzeta(i,j)=ad_gzeta(i,j)-tmp1
        ad_gzeta2(i,j-1)=ad_gzeta2(i,j-1)+tmp
        ad_gzeta2(i,j)=ad_gzeta2(i,j)-tmp

        ad_rhs_vbar(i,j)=0.d0

       END DO
      END DO

      DO j=JstrU,JendU
       DO i=IstrU,IendU
! NL:        rhs_ubar(i,j)=cff1*dy_u(i,j)*(
! NL:     *        (h(i-1,j)+h(i,j))*(gzeta(i-1,j)-gzeta(i,j))+
! NL:     *        gzeta2(i-1,j)-gzeta2(i,j)  )

! TL:        tl_rhs_ubar(i,j)=cff1*dy_u(i,j)*(
! TL:     *     (h(i-1,j)+h(i,j))*(tl_gzeta(i-1,j)-tl_gzeta(i,j))+
! TL:     *     tl_gzeta2(i-1,j)-tl_gzeta2(i,j) )

        tmp=cff1*dy_u(i,j)*ad_rhs_ubar(i,j)
        tmp1=tmp*(h(i-1,j)+h(i,j))

        ad_gzeta(i-1,j)=ad_gzeta(i-1,j)+tmp1
        ad_gzeta(i,j)=ad_gzeta(i,j)-tmp1
        ad_gzeta2(i-1,j)=ad_gzeta2(i-1,j)+tmp
        ad_gzeta2(i,j)=ad_gzeta2(i,j)-tmp

        ad_rhs_ubar(i,j)=0.d0

       END DO
      END DO

!=======================================================================
!  Time step free-surface equation.
!=======================================================================

# if defined NS_PERIODIC
      CALL ad_exchange_2d (ad_gzeta2,LBi,UBi,LBj,UBj,Lm,Mm) 
      CALL ad_exchange_2d (ad_gzeta,LBi,UBi,LBj,UBj,Lm,Mm) 
      CALL ad_exchange_2d (ad_Dnew,LBi,UBi,LBj,UBj,Lm,Mm) 
      CALL ad_exchange_2d (ad_zeta(:,:,knew),LBi,UBi,LBj,UBj,Lm,Mm) 
# endif

!  ADj to free-surface lateral boundary conditions.
!
      CALL ad_zetabc (Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj,  
     *             krhs,kstp,knew,ad_zeta)

!  ADj to: if predictor step, load right-side-term into shared array.
!
      IF (PREDICTOR) THEN
! v3.1:
!       DO j=Jstr,Jend
!        DO i=Istr,Iend
! v.3.2:
       DO j=Jstr-1,Jend+1
        DO i=Istr-1,Iend+1
! TL:         tl_rzeta(i,j,krhs)=tl_rhs_zeta(i,j)
         ad_rhs_zeta(i,j)=ad_rhs_zeta(i,j)+ad_rzeta(i,j,krhs)
         ad_rzeta(i,j,krhs)=0.d0
        END DO
       END DO
      END IF

!
!  ADJ to: Load new free-surface values into shared array at both predictor
!  and corrector steps.
!
      DO j=Jstr,Jend
       DO i=Istr,Iend
! TL:        tl_zeta(i,j,knew)=tl_zeta_new(i,j)
        ad_zeta_new(i,j)=ad_zeta_new(i,j)+ad_zeta(i,j,knew)
        ad_zeta(i,j,knew)=0.d0
       END DO
      END DO

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

         gzeta(i,j)=zeta(i,j,1)

! TL line 5:         tl_gzeta2(i,j)=2.d0*tl_gzeta(i,j)*gzeta(i,j)
         ad_gzeta(i,j)=ad_gzeta(i,j)+2.d0*gzeta(i,j)*ad_gzeta2(i,j)
         ad_gzeta2(i,j)=0.d0

! TL line 4:         tl_gzeta(i,j)=0.5d0*(tl_zeta(i,j,kstp)+tl_zeta_new(i,j))
         tmp=ad_gzeta(i,j)*0.5d0
         ad_zeta(i,j,kstp)=ad_zeta(i,j,kstp)+tmp
         ad_zeta_new(i,j)=ad_zeta_new(i,j)+tmp
         ad_gzeta(i,j)=0.d0

! TL line 3:         tl_Dnew(i,j)=tl_zeta_new(i,j)
         ad_zeta_new(i,j)=ad_zeta_new(i,j)+ad_Dnew(i,j)
         ad_Dnew(i,j)=0.d0

! TL line 2:         tl_zeta_new(i,j)=tl_zeta(i,j,kstp)+
! TL line 2:     *                  dtfast*oArea(i,j)*tl_rhs_zeta(i,j)
         ad_zeta(i,j,kstp)=ad_zeta(i,j,kstp)+ad_zeta_new(i,j)
         ad_rhs_zeta(i,j)=ad_rhs_zeta(i,j)
     *                     +dtfast*oArea(i,j)*ad_zeta_new(i,j)
         ad_zeta_new(i,j)=0.d0

! TL line 1:         tl_rhs_zeta(i,j)=
! TL line 1:     *              tl_DU(i,j)-tl_DU(i+1,j)+tl_DV(i,j)-tl_DV(i,j+1)
         ad_DU(i,j)=ad_DU(i,j)+ad_rhs_zeta(i,j)
         ad_DU(i+1,j)=ad_DU(i+1,j)-ad_rhs_zeta(i,j)
         ad_DV(i,j)=ad_DV(i,j)+ad_rhs_zeta(i,j)
         ad_DV(i,j+1)=ad_DV(i,j+1)-ad_rhs_zeta(i,j)
         ad_rhs_zeta(i,j)=0.d0

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

         gzeta(i,j)=zeta(i,j,1)

! TL line 5:         tl_gzeta2(i,j)=2.d0*tl_gzeta(i,j)*gzeta(i,j)
         ad_gzeta(i,j)=ad_gzeta(i,j)+2.d0*gzeta(i,j)*ad_gzeta2(i,j)
         ad_gzeta2(i,j)=0.d0

! TL line 4:         tl_gzeta(i,j)=cff5*tl_zeta(i,j,krhs)+
! TL line 4:     *                 cff4*(tl_zeta(i,j,kstp)+tl_zeta_new(i,j))

         tmp=cff4*ad_gzeta(i,j)

         ad_zeta(i,j,krhs)=ad_zeta(i,j,krhs)+cff5*ad_gzeta(i,j)
         ad_zeta(i,j,kstp)=ad_zeta(i,j,kstp)+tmp
         ad_zeta_new(i,j)=ad_zeta_new(i,j)+tmp

         ad_gzeta(i,j)=0.d0

! TL line 3:         tl_Dnew(i,j)=tl_zeta_new(i,j)
         ad_zeta_new(i,j)=ad_zeta_new(i,j)+ad_Dnew(i,j)
         ad_Dnew(i,j)=0.d0

! TL line 2:         tl_zeta_new(i,j)=tl_zeta(i,j,kstp)+
! TL line 2:     *                 cff1*oArea(i,j)*tl_rhs_zeta(i,j)
         ad_zeta(i,j,kstp)=ad_zeta(i,j,kstp)+ad_zeta_new(i,j)
         ad_rhs_zeta(i,j)=ad_rhs_zeta(i,j)
     *                      +cff1*oArea(i,j)*ad_zeta_new(i,j)
         ad_zeta_new(i,j)=0.d0

! TL line 1:         tl_rhs_zeta(i,j)=
! TL line 1:     *     tl_DU(i,j)-tl_DU(i+1,j)+tl_DV(i,j)-tl_DV(i,j+1)

         ad_DU(i,j)=ad_DU(i,j)+ad_rhs_zeta(i,j)
         ad_DU(i+1,j)=ad_DU(i+1,j)-ad_rhs_zeta(i,j)
         ad_DV(i,j)=ad_DV(i,j)+ad_rhs_zeta(i,j)
         ad_DV(i,j+1)=ad_DV(i,j+1)-ad_rhs_zeta(i,j)

         ad_rhs_zeta(i,j)=0.d0

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
         gzeta(i,j)=zeta(i,j,1)

! TL line 5:         tl_gzeta2(i,j)=2.d0*tl_gzeta(i,j)*gzeta(i,j)
         ad_gzeta(i,j)=ad_gzeta(i,j)+2.d0*gzeta(i,j)*ad_gzeta2(i,j)
         ad_gzeta2(i,j)=0.d0
         
! TL line 4:         tl_gzeta(i,j)=cff5*tl_zeta_new(i,j)+cff4*tl_zeta(i,j,krhs)
         ad_zeta_new(i,j)=ad_zeta_new(i,j)+cff5*ad_gzeta(i,j)
         ad_zeta(i,j,krhs)=ad_zeta(i,j,krhs)+cff4*ad_gzeta(i,j)
         ad_gzeta(i,j)=0.d0

! TL line 3:         tl_Dnew(i,j)=tl_zeta_new(i,j)
         ad_zeta_new(i,j)=ad_zeta_new(i,j)+ad_Dnew(i,j)
         ad_Dnew(i,j)=0.d0

! TL line 2:         tl_zeta_new(i,j)=tl_zeta(i,j,kstp)+
! TL line 2:     *     oArea(i,j)*( tl_cff
! TL line 2:     *                 +cff2*tl_rzeta(i,j,kstp)
! TL line 2:     *                 -cff3*tl_rzeta(i,j,ptsk))

         tmp=oArea(i,j)*ad_zeta_new(i,j)

         ad_zeta(i,j,kstp)=ad_zeta(i,j,kstp)+ad_zeta_new(i,j)
         ad_cff=ad_cff+tmp
         ad_rzeta(i,j,kstp)=ad_rzeta(i,j,kstp)+cff2*tmp
         ad_rzeta(i,j,ptsk)=ad_rzeta(i,j,ptsk)-cff3*tmp

         ad_zeta_new(i,j)=0.d0

! TL line 1:         tl_cff=cff1*(tl_DU(i,j)-tl_DU(i+1,j)+tl_DV(i,j)-tl_DV(i,j+1))
 
         tmp=cff1*ad_cff

         ad_DU(i,j)=ad_DU(i,j)+tmp
         ad_DU(i+1,j)=ad_DU(i+1,j)-tmp
         ad_DV(i,j)=ad_DV(i,j)+tmp
         ad_DV(i,j+1)=ad_DV(i,j+1)-tmp     

         ad_cff=0.d0

        END DO
       END DO

      END IF

      END IF !! <- "IF (iif.ne.nfast+1)"

! ADJ to: Finalize averaging, compute new tl_Hz etc. 
!
      IF (PREDICTOR .AND. (iif.eq.nfast+1)) THEN

       CALL ad_set_depth(LBi,UBi,LBj,UBj,N,Lm,Mm,
     *               z0_r,z0_w,oh,
     *               ad_Zt_avg1,ad_z_r,ad_z_w,ad_Hz,ad_Hz_u,ad_Hz_v)

# if defined NS_PERIODIC
       CALL ad_exchange_2d(ad_Zt_avg1,LBi,UBi,LBj,UBj,Lm,Mm)
       CALL ad_exchange_2d(ad_DU_avg1,LBi,UBi,LBj,UBj,Lm,Mm)
       CALL ad_exchange_2d(ad_DV_avg1,LBi,UBi,LBj,UBj,Lm,Mm)
       CALL ad_exchange_2d(ad_DU_avg2,LBi,UBi,LBj,UBj,Lm,Mm)
       CALL ad_exchange_2d(ad_DV_avg2,LBi,UBi,LBj,UBj,Lm,Mm)
# endif

      END IF ! PREDICTOR, iif=nfast+1

!-----------------------------------------------------------------------
!  ADJ to: compute time averaged fields over all short time-steps.
!-----------------------------------------------------------------------
      IF (PREDICTOR) THEN

       IF (iif.eq.1) THEN
        cff2=(-1.d0/12.d0)*weight(2,iif+1)
!
! ADJ: Reset arrays for 2D fields averaged within the short time-steps.
!

! TL loop 3:
        DO j=JV_RANGE
         DO i=IV_RANGE
! tl line 2:          tl_DV_avg2(i,j)=cff2*tl_DV(i,j)
          ad_DV(i,j)=ad_DV(i,j)+cff2*ad_DV_avg2(i,j) 
          ad_DV_avg2(i,j)=0.d0

! tl line 1:          tl_DV_avg1(i,j)=0.d0
          ad_DV_avg1(i,j)=0.d0

         END DO
        END DO

! TL loop 2:
        DO j=JU_RANGE
         DO i=IU_RANGE
! tl line 2:          tl_DU_avg2(i,j)=cff2*tl_DU(i,j)
          ad_DU(i,j)=ad_DU(i,j)+cff2*ad_DU_avg2(i,j)
          ad_DU_avg2(i,j)=0.d0

! tl line 1:          tl_DU_avg1(i,j)=0.d0
          ad_DU_avg1(i,j)=0.d0
 
         END DO
        END DO

! TL loop 1:

! old limits (v2.1)
!        DO j=Jstr,Jend
!         DO i=Istr,Iend

! Note: new limits (consistent with v3.1, 3.2)
        DO j=LBj,UBj
         DO i=LBi,UBi
! TL:          tl_Zt_avg1(i,j)=0.d0
          ad_Zt_avg1(i,j)=0.d0
         END DO
        END DO

       ELSE ! (PREDICTOR, iif>1)
!
!  ADJ to: Accumulate field averages of previous time-step after they are
!  computed in the previous corrector step, updated their boundaries,
!  and synchronized.
!

        cff1=weight(1,iif-1)
        cff2=(8.d0/12.d0)*weight(2,iif)-
     *       (1.d0/12.d0)*weight(2,iif+1)
        
! TL loop 3:
        DO j=JV_RANGE
         DO i=IV_RANGE
! tl line 2 (recurs):          tl_DV_avg2(i,j)=tl_DV_avg2(i,j)+cff2*tl_DV(i,j)
          ad_DV(i,j)=ad_DV(i,j)+cff2*ad_DV_avg2(i,j)  
          
! tl line 1 (recurs):          tl_DV_avg1(i,j)=tl_DV_avg1(i,j)+cff1*tl_DV(i,j)
          ad_DV(i,j)=ad_DV(i,j)+cff1*ad_DV_avg1(i,j)

         END DO
        END DO

! TL loop 2:
        DO j=JU_RANGE
         DO i=IU_RANGE

! tl line 2 (recurs):          tl_DU_avg2(i,j)=tl_DU_avg2(i,j)+cff2*tl_DU(i,j)
          ad_DU(i,j)=ad_DU(i,j)+cff2*ad_DU_avg2(i,j)

! tl line 1 (recurs):          tl_DU_avg1(i,j)=tl_DU_avg1(i,j)+cff1*tl_DU(i,j)
          ad_DU(i,j)=ad_DU(i,j)+cff1*ad_DU_avg1(i,j)

         END DO
        END DO

! TL loop 1:

! old limits (v2.1)
!        DO j=Jstr,Jend
!         DO i=Istr,Iend

! Note (v2.1_test, v3.1): new limits 
        DO j=LBj,UBj
         DO i=LBi,UBi
! tl (recurs):          tl_Zt_avg1(i,j)=tl_Zt_avg1(i,j)+cff1*tl_zeta(i,j,krhs)
          ad_zeta(i,j,krhs)=ad_zeta(i,j,krhs)+cff1*ad_Zt_avg1(i,j)
         END DO
        END DO

       END IF ! (iif.eq.1) or else

      ELSE ! .NOT.PREDICTOR: CORRECTOR STEP

       IF (iif.eq.1) THEN
        cff2=weight(2,iif)
       ELSE
        cff2=(5.d0/12.d0)*weight(2,iif)
       END IF
       
! TL loop 2:
       DO j=JV_RANGE
        DO i=IV_RANGE
! tl (recurs):         tl_DV_avg2(i,j)=tl_DV_avg2(i,j)+cff2*tl_DV(i,j)
         ad_DV(i,j)=ad_DV(i,j)+cff2*ad_DV_avg2(i,j)
        END DO
       END DO

! TL loop 1:
       DO j=JU_RANGE
        DO i=IU_RANGE
! tl (recurs):         tl_DU_avg2(i,j)=tl_DU_avg2(i,j)+cff2*tl_DU(i,j)
         ad_DU(i,j)=ad_DU(i,j)+cff2*ad_DU_avg2(i,j)
        END DO
       END DO

      END IF  ! choices PREDICTOR or ELSE

!-----------------------------------------------------------------------
!  ADJ to: compute total depth (m) and vertically integrated mass fluxes.
!-----------------------------------------------------------------------
      
# ifdef NS_PERIODIC
      CALL ad_exchange_2d(ad_DU,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL ad_exchange_2d(ad_DV,LBi,UBi,LBj,UBj,Lm,Mm)
# endif

      DO j=JV_RANGE
       DO i=IV_RANGE

! NL:        DV(i,j)=0.5d0*(Drhs(i,j)+Drhs(i,j-1))*vbar(i,j,krhs)*dx_v(i,j)

! TL:        tl_DV(i,j)=0.5d0*dx_v(i,j)*
! TL:     *    ( (tl_zeta(i,j,krhs)+tl_zeta(i,j-1,krhs))*vbar(i,j,krhs)+
! TL:     *      (Drhs(i,j)+Drhs(i,j-1))*tl_vbar(i,j,krhs) )
        
        tmp=0.5d0*dx_v(i,j)*ad_DV(i,j)
        tmp1=tmp*vbar(i,j,krhs)
 
        ad_zeta(i,j,krhs)=ad_zeta(i,j,krhs)+tmp1
        ad_zeta(i,j-1,krhs)=ad_zeta(i,j-1,krhs)+tmp1
        ad_vbar(i,j,krhs)=ad_vbar(i,j,krhs)+tmp*(Drhs(i,j)+Drhs(i,j-1))

        ad_DV(i,j)=0.d0

       END DO
      END DO

      DO j=JU_RANGE
       DO i=IU_RANGE

! NL:        DU(i,j)=0.5d0*(Drhs(i,j)+Drhs(i-1,j))*ubar(i,j,krhs)*dy_u(i,j)

! TL:        tl_DU(i,j)=0.5d0*dy_u(i,j)*
! TL:     *    ( (tl_zeta(i,j,krhs)+tl_zeta(i-1,j,krhs))*ubar(i,j,krhs) +
! TL:     *      (Drhs(i,j)+Drhs(i-1,j))*tl_ubar(i,j,krhs) )

        tmp=0.5d0*dy_u(i,j)*ad_DU(i,j)
        tmp1=tmp*ubar(i,j,krhs)

        ad_zeta(i,j,krhs)=ad_zeta(i,j,krhs)+tmp1
        ad_zeta(i-1,j,krhs)=ad_zeta(i-1,j,krhs)+tmp1
        ad_ubar(i,j,krhs)=ad_ubar(i,j,krhs)+tmp*(Drhs(i,j)+Drhs(i-1,j))

        ad_DU(i,j)=0.d0

       END DO
      END DO

#ifdef NS_PERIODIC
      CALL ad_exchange_2d(ad_Drhs,LBi,UBi,LBj,UBj,Lm,Mm)
#endif

# ifdef NS_PERIODIC
      DO j=Jstr,Jend
# else
      DO j=Jstr-1,Jend+1
# endif
       DO i=Istr-1,Iend+1
! TL:        tl_Drhs(i,j)=tl_zeta(i,j,krhs)
        ad_zeta(i,j,krhs)=ad_zeta(i,j,krhs)+ad_Drhs(i,j)
        ad_Drhs(i,j)=0.d0
       END DO
      END DO



      RETURN
      END

#endif /* ADJOINT */
