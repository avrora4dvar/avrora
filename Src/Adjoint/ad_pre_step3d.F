#include "cppdefs.h"

#ifdef ADJOINT 

      SUBROUTINE ad_pre_step3d
      USE mod_ocean
      USE mod_forces

      implicit none
      integer :: indx,i,j,k,itrc
      real(8) :: cff1,cff2,cff,ad_cff1,ad_cff2,tmp,tmp1,tmp2,ad_tmp

      real(8) :: DC(LBi:UBi,N)
      real(8) :: FC(LBi:UBi,0:N)

      real(8) :: ad_FX(LBi:UBi,LBj:UBj)
      real(8) :: ad_FE(LBi:UBi,LBj:UBj)
      real(8) :: ad_DC(LBi:UBi,N)
      real(8) :: ad_FC(LBi:UBi,0:N)

      real(8), parameter :: Gamma = 1.0d0/6.0d0

      ad_cff1=0.d0
      ad_cff2=0.d0
      ad_FX=0.d0
      ad_FE=0.d0
      ad_DC=0.d0
      ad_FC=0.d0 

      ad_tmp=0.d0

      indx=3-nrhs

!=======================================================================
! ADJ TO: adding svstr and bvstr (ATTN: recursive statements in the TL)
!=======================================================================
      DO j=JstrV,JendV
       DO i=IstrV,IendV

        cff=dt*DBLE(mask_v(i,j))

! TL line 2:        tl_v(i,j,N,nnew)=tl_v(i,j,N,nnew)+cff*tl_svstr(i,j)
        ad_svstr(i,j)=ad_svstr(i,j)+cff*ad_v(i,j,N,nnew)
 
! TL line 1:        tl_v(i,j,1,nnew)=tl_v(i,j,1,nnew)-cff*tl_bvstr(i,j)
        ad_bvstr(i,j)=ad_bvstr(i,j)-cff*ad_v(i,j,1,nnew)

       END DO
      END DO

!=======================================================================
!  ADJ TO: 3D momentum equation in the ETA-direction.
!=======================================================================
!

      IF (iic.eq.1) THEN

       DO k=1,N
        DO j=JstrV,JendV
         DO i=IstrV,IendV

! NL:           v(i,j,k,nnew)=v(i,j,k,nstp)*Hz_v(i,j,k)*DBLE(mask_v(i,j))
! TL:           tl_v(i,j,k,nnew)=(tl_v(i,j,k,nstp)*Hz_v(i,j,k)
! TL:     *                +v(i,j,k,nstp)*tl_Hz_v(i,j,k) )*DBLE(mask_v(i,j))

          ad_v(i,j,k,nstp)=ad_v(i,j,k,nstp)+
     *                 DBLE(mask_v(i,j))*Hz_v(i,j,k)*ad_v(i,j,k,nnew)
          ad_Hz_v(i,j,k)=ad_Hz_v(i,j,k)+
     *                 DBLE(mask_v(i,j))*v(i,j,k,nstp)*ad_v(i,j,k,nnew)
          ad_v(i,j,k,nnew)=0.d0

         END DO
        END DO
       END DO


      ELSEIF (iic.eq.2) THEN       

       ! Adam-Bashforth, iic==2: scaling coef. for rhs(n-1)
       cff1=dt*0.5d0

       DO k=1,N
        DO j=JstrV,JendV
         DO i=IstrV,IendV
! NL:          v(i,j,k,nnew)=v(i,j,k,nstp)*Hz_v(i,j,k)-
! NL:     *                   cff1*oArea_v(i,j)*rv(i,j,k,indx)

! TL:           tl_v(i,j,k,nnew)=tl_v(i,j,k,nstp)*Hz_v(i,j,k)
! TL:     *                     +v(i,j,k,nstp)*tl_Hz_v(i,j,k)
! TL:     *                     -cff1*oArea_v(i,j)*tl_rv(i,j,k,indx)


          ad_v(i,j,k,nstp)=ad_v(i,j,k,nstp)
     *                                   +Hz_v(i,j,k)*ad_v(i,j,k,nnew)
          ad_Hz_v(i,j,k)=ad_Hz_v(i,j,k)+v(i,j,k,nstp)*ad_v(i,j,k,nnew)
          ad_rv(i,j,k,indx)=ad_rv(i,j,k,indx)-
     *                        cff1*oArea_v(i,j)*ad_v(i,j,k,nnew)
          ad_v(i,j,k,nnew)=0.d0
 
         END DO
        END DO
       END DO

      ELSE
       ! Adam-Bashforth, iic>2: scaling coef. for rhs(n-2), rhs(n-1)
       cff1=dt* 5.0d0/12.0d0
       cff2=dt*16.0d0/12.0d0

       DO k=1,N
        DO j=JstrV,JendV
         DO i=IstrV,IendV
! NL:          v(i,j,k,nnew)=v(i,j,k,nstp)*Hz_v(i,j,k)+
! NL:     *          oArea_v(i,j)*(cff1*rv(i,j,k,nrhs)-cff2*rv(i,j,k,indx))

! TL:           tl_v(i,j,k,nnew)=tl_v(i,j,k,nstp)*Hz_v(i,j,k)
! TL:     *       +v(i,j,k,nstp)*tl_Hz_v(i,j,k)
! TL:     *       +oArea_v(i,j)*(cff1*tl_rv(i,j,k,nrhs)-cff2*tl_rv(i,j,k,indx))

          ad_v(i,j,k,nstp)=ad_v(i,j,k,nstp)+
     *                                     Hz_v(i,j,k)*ad_v(i,j,k,nnew)
          ad_Hz_v(i,j,k)=ad_Hz_v(i,j,k)+v(i,j,k,nstp)*ad_v(i,j,k,nnew)
          ad_rv(i,j,k,nrhs)=ad_rv(i,j,k,nrhs)+
     *                                cff1*oArea_v(i,j)*ad_v(i,j,k,nnew)
          ad_rv(i,j,k,indx)=ad_rv(i,j,k,indx)-
     *                                cff2*oArea_v(i,j)*ad_v(i,j,k,nnew)
          ad_v(i,j,k,nnew)=0.d0
         END DO
        END DO
       END DO

      END IF

!=======================================================================
! ADJ TO: adding sustr and bustr (ATTN: recursive statements in the TL)
!=======================================================================
      DO j=JstrU,JendU
       DO i=IstrU,IendU

        cff=dt*DBLE(mask_u(i,j))

! TL line 2:        tl_u(i,j,N,nnew)=tl_u(i,j,N,nnew)+cff*tl_sustr(i,j)
        ad_sustr(i,j)=ad_sustr(i,j)+cff*ad_u(i,j,N,nnew)
 
! TL line 1:        tl_u(i,j,1,nnew)=tl_u(i,j,1,nnew)-cff*tl_bustr(i,j)
        ad_bustr(i,j)=ad_bustr(i,j)-cff*ad_u(i,j,1,nnew)

       END DO
      END DO

!=======================================================================
!  ADJ TO: 3D momentum equation in the XI-direction.
!=======================================================================
!

      IF (iic.eq.1) THEN

       DO k=1,N
        DO j=JstrU,JendU
         DO i=IstrU,IendU
! NL:           u(i,j,k,nnew)=u(i,j,k,nstp)*Hz_u(i,j,k)*DBLE(mask_u(i,j))

! TL:          tl_u(i,j,k,nnew)=(tl_u(i,j,k,nstp)*Hz_u(i,j,k)+
! TL:     *                      u(i,j,k,nstp)*tl_Hz_u(i,j,k))
! TL:     *                     *DBLE(mask_u(i,j))

          ad_u(i,j,k,nstp)=ad_u(i,j,k,nstp)+
     *                 DBLE(mask_u(i,j))*Hz_u(i,j,k)*ad_u(i,j,k,nnew)
          ad_Hz_u(i,j,k)=ad_Hz_u(i,j,k)+
     *                 DBLE(mask_u(i,j))*u(i,j,k,nstp)*ad_u(i,j,k,nnew)
          ad_u(i,j,k,nnew)=0.d0

         END DO
        END DO
       END DO


      ELSEIF (iic.eq.2) THEN       

       ! Adam-Bashforth, iic==2: scaling coef. for rhs(n-1)
       cff1=dt*0.5d0

       DO k=1,N
        DO j=JstrU,JendU
         DO i=IstrU,IendU

! NL:           u(i,j,k,nnew)=u(i,j,k,nstp)*Hz_u(i,j,k)-
! NL:     *                    cff1*oArea_u(i,j)*ru(i,j,k,indx)

! TL:           tl_u(i,j,k,nnew)=tl_u(i,j,k,nstp)*Hz_u(i,j,k)+
! TL:     *                      u(i,j,k,nstp)*tl_Hz_u(i,j,k)-
! TL:     *                       cff1*oArea_u(i,j)*tl_ru(i,j,k,indx)

          ad_u(i,j,k,nstp)=ad_u(i,j,k,nstp)
     *                                    +Hz_u(i,j,k)*ad_u(i,j,k,nnew)
          ad_Hz_u(i,j,k)=ad_Hz_u(i,j,k)+u(i,j,k,nstp)*ad_u(i,j,k,nnew)
          ad_ru(i,j,k,indx)=ad_ru(i,j,k,indx)-
     *                           cff1*oArea_u(i,j)*ad_u(i,j,k,nnew)
          ad_u(i,j,k,nnew)=0.d0
 
          END DO
        END DO
       END DO

      ELSE
       ! Adam-Bashforth, iic>2: scaling coef. for rhs(n-2), rhs(n-1)
       cff1=dt* 5.0d0/12.0d0
       cff2=dt*16.0d0/12.0d0

       DO k=1,N
        DO j=JstrU,JendU
         DO i=IstrU,IendU

! NL:           u(i,j,k,nnew)=u(i,j,k,nstp)*Hz_u(i,j,k)-
! NL:     *                    cff1*oArea_u(i,j)*ru(i,j,k,indx)

! TL:           tl_u(i,j,k,nnew)=tl_u(i,j,k,nstp)*Hz_u(i,j,k)+
! TL:     *                      u(i,j,k,nstp)*tl_Hz_u(i,j,k)-
! TL:     *                       cff1*oArea_u(i,j)*tl_ru(i,j,k,indx)


          ad_u(i,j,k,nstp)=ad_u(i,j,k,nstp)+
     *                                    Hz_u(i,j,k)*ad_u(i,j,k,nnew)
          ad_Hz_u(i,j,k)=ad_Hz_u(i,j,k)+u(i,j,k,nstp)*ad_u(i,j,k,nnew)
          ad_ru(i,j,k,nrhs)=ad_ru(i,j,k,nrhs)+
     *                       cff1*oArea_u(i,j)*ad_u(i,j,k,nnew)
          ad_ru(i,j,k,indx)=ad_ru(i,j,k,indx)-
     *                       cff2*oArea_u(i,j)*ad_u(i,j,k,nnew)
          ad_u(i,j,k,nnew)=0.d0
         END DO
        END DO
       END DO

      END IF  ! (xi-eqns:choices for iic, 1, 2, or >2)

#ifndef TS_FIXED 
!=======================================================================
!  ADJ to: Tracer equation(s).
!=======================================================================
!

!
!-----------------------------------------------------------------------
!  ADJ to: Start computation of tracers at n+1 time-step, t(i,j,k,nnew,itrc).
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
       DO itrc=1,Ntr

        DO i=Istr,Iend
! TL:         tl_t(i,j,N,nnew,itrc)=tl_t(i,j,N,nnew,itrc)+
!     *                          dt*tl_stflx(i,j,itrc)
! (recurs)    

         ad_stflx(i,j,itrc)=ad_stflx(i,j,itrc)+dt*ad_t(i,j,N,nnew,itrc)

        END DO

        DO k=1,N
         DO i=Istr,Iend
! TL:          tl_t(i,j,k,nnew,itrc)=tl_Hz(i,j,k)*t(i,j,k,nstp,itrc)+
!     *                          Hz(i,j,k)*tl_t(i,j,k,nstp,itrc)

          ad_Hz(i,j,k)=ad_Hz(i,j,k)+
     *                  t(i,j,k,nstp,itrc)*ad_t(i,j,k,nnew,itrc)
          ad_t(i,j,k,nstp,itrc)=ad_t(i,j,k,nstp,itrc)+
     *                  Hz(i,j,k)*ad_t(i,j,k,nnew,itrc)
          ad_t(i,j,k,nnew,itrc)=0.d0

         END DO
        END DO


       END DO   ! DO itrc=1,Ntr
      END DO   ! DO j=Jstr,Jend

!-----------------------------------------------------------------------
! ADJ to: Boundary conditions for the intermediate tracer:
!-----------------------------------------------------------------------

      DO itrc=1,Ntr
#ifdef NS_PERIODIC
       CALL ad_exchange_3d(ad_t(:,:,:,3,itrc),LBi,UBi,LBj,UBj,1,N,Lm,Mm)
#endif

       CALL ad_t3dbc(Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj,N,Ntr,
#ifdef MPI
     *             B_West,B_East,B_South,B_North, ! <- new line, comp ot v2.1
#endif
     *            nrhs, nstp, 3,itrc,ad_t)

      END DO

!
!-----------------------------------------------------------------------
! ADJ to: Compute intermediate tracer at n+1/2 time-step, t(i,j,k,3,itrc).
!-----------------------------------------------------------------------
!

!
! Pipelined J-loop
!

      J_LOOP1 : DO j=Jstr,Jend

! FWD fields, artificial continuity equation:

       IF (iic.eq.1) THEN
        cff=0.5d0*dt
       ELSE
        cff=(1.0d0-Gamma)*dt
       END IF

       DO k=1,N
        DO i=Istr,Iend
         cff1=cff*oArea(i,j)
         tmp=Hz(i,j,k)-
     *       cff1*
     *        (Hudy(i+1,j,k,nstp)-Hudy(i,j,k,nstp)+
     *         Hvdx(i,j+1,k,nstp)-Hvdx(i,j,k,nstp)+
     *         W(i,j,k)-W(i,j,k-1) )

         DC(i,k)=1.d0/tmp

        END DO
       END DO

!
!  Compute time rate of change of intermediate tracer due to vertical
!  advection.  Impose artificial continuity equation.
!
       T_LOOP2: DO itrc=1,Ntr

! FWD field, vertical advective flux:

        DO k=1,N-1
         DO i=Istr,Iend
          FC(i,k)=W(i,j,k)*   
     *                0.5d0*(t(i,j,k  ,nstp,itrc)+ 
     *                       t(i,j,k+1,nstp,itrc))
         END DO
        END DO
        
        DO i=Istr,Iend
         FC(i,0)=0.0d0
         FC(i,N)=0.0d0
        END DO

!
! ADJ to: Time-step vertical advection of tracers (Tunits).
!
        DO k=1,N
         DO i=Istr,Iend

          cff1=cff*oArea(i,j)

! TL:          tl_t(i,j,k,3,itrc)=
!     *     tl_DC(i,k)*
!     *         ( t(i,j,k,3,itrc)-cff1*(FC(i,k)-FC(i,k-1)) ) +
!     *     DC(i,k)*
!     *         ( tl_t(i,j,k,3,itrc)-cff1*(tl_FC(i,k)-tl_FC(i,k-1)) )
! (recurs., scaled): 

          ad_DC(i,k)=ad_DC(i,k)+
     *     ( t(i,j,k,3,itrc)-cff1*(FC(i,k)-FC(i,k-1)) )*
     *     ad_t(i,j,k,3,itrc)

          tmp1=DC(i,k)*cff1*ad_t(i,j,k,3,itrc)
          ad_FC(i,k)=ad_FC(i,k)-tmp1
          ad_FC(i,k-1)=ad_FC(i,k-1)+tmp1

          ad_t(i,j,k,3,itrc)=ad_t(i,j,k,3,itrc)*DC(i,k)
        
         END DO
        END DO

!
!  ADJ to: Second-order, central differences vertical advective flux.
!
        DO i=Istr,Iend
! TL:         tl_FC(i,N)=0.0d0
         ad_FC(i,N)=0.0d0

! TL:         tl_FC(i,0)=0.0d0
         ad_FC(i,0)=0.0d0
        END DO

        DO k=1,N-1
         DO i=Istr,Iend
! TL:          tl_FC(i,k)= 0.5d0*
!     *    ( tl_W(i,j,k)*(   t(i,j,k  ,nstp,itrc)+   t(i,j,k+1,nstp,itrc))+
!     *         W(i,j,k)*(tl_t(i,j,k  ,nstp,itrc)+tl_t(i,j,k+1,nstp,itrc)) )

          tmp1=0.5d0*ad_FC(i,k)           
          tmp2=tmp1*W(i,j,k)

          ad_W(i,j,k)=ad_W(i,j,k)+
     *                tmp1*(t(i,j,k  ,nstp,itrc)+t(i,j,k+1,nstp,itrc))
          ad_t(i,j,k  ,nstp,itrc)=ad_t(i,j,k  ,nstp,itrc)+tmp2
          ad_t(i,j,k+1,nstp,itrc)=ad_t(i,j,k+1,nstp,itrc)+tmp2
          ad_FC(i,k)=0.d0

         END DO
        END DO

       END DO T_LOOP2

!
!  ADJ to: Compute artificial continuity equation (same for all tracers) and
!  load it into private array DC (1/m).
!
! (NOTE: fwd DC)
       DO k=1,N
        DO i=Istr,Iend

         cff1=cff*oArea(i,j)

! tl:         tl_DC(i,k)=-tl_tmp*DC(i,k)*DC(i,k)
         ad_tmp=ad_tmp-ad_DC(i,k)*DC(i,k)*DC(i,k)
         ad_DC(i,k)=0.d0

! TL:         tl_tmp=tl_Hz(i,j,k)-
!     *          cff1*
!     *           (tl_Hudy(i+1,j,k)-tl_Hudy(i,j,k)+
!     *            tl_Hvdx(i,j+1,k)-tl_Hvdx(i,j,k)+
!     *            tl_W(i,j,k)-tl_W(i,j,k-1) )

         tmp1=cff1*ad_tmp

         ad_Hz(i,j,k)=ad_Hz(i,j,k)+ad_tmp
         ad_Hudy(i+1,j,k)=ad_Hudy(i+1,j,k)-tmp1
         ad_Hudy(i,j,k)=ad_Hudy(i,j,k)    +tmp1
         ad_Hvdx(i,j+1,k)=ad_Hvdx(i,j+1,k)-tmp1
         ad_Hvdx(i,j,k)=ad_Hvdx(i,j,k)    +tmp1
         ad_W(i,j,k)=ad_W(i,j,k)          -tmp1
         ad_W(i,j,k-1)=ad_W(i,j,k-1)      +tmp1
         ad_tmp=0.d0

        END DO
       END DO

      END DO J_LOOP1

!
!  Compute time rate of change of intermediate tracer due to
!  horizontal advection.
!

      T_LOOP1 :DO itrc=1,Ntr
       K_LOOP: DO k=1,N

!
!  Time-step horizontal advection (m Tunits).
!
        IF (iic.eq.1) THEN
         cff=0.5d0*dt
         cff1=1.0d0
         cff2=0.0d0
        ELSE
         cff=(1.0d0-Gamma)*dt
         cff1=0.5d0+Gamma
         cff2=0.5d0-Gamma
        END IF

        DO j=Jstr,Jend
         DO i=Istr,Iend

! TL:          tl_t(i,j,k,3,itrc)=
!     *      tl_Hz(i,j,k)*
!     *       (cff1*t(i,j,k,nstp,itrc)+cff2*t(i,j,k,nnew,itrc))+
!     *      Hz(i,j,k)*
!     *       (cff1*tl_t(i,j,k,nstp,itrc)+cff2*tl_t(i,j,k,nnew,itrc))-
!     *      cff*oArea(i,j)*                
!     *       (tl_FX(i+1,j)-tl_FX(i,j)+tl_FE(i,j+1)-tl_FE(i,j))

          tmp1=Hz(i,j,k)*ad_t(i,j,k,3,itrc)
          tmp2=cff*oArea(i,j)*ad_t(i,j,k,3,itrc)

          ad_Hz(i,j,k)=ad_Hz(i,j,k)+
     *     (cff1*t(i,j,k,nstp,itrc)+cff2*t(i,j,k,nnew,itrc))*
     *     ad_t(i,j,k,3,itrc)
          ad_t(i,j,k,nstp,itrc)=ad_t(i,j,k,nstp,itrc)+cff1*tmp1
          ad_t(i,j,k,nnew,itrc)=ad_t(i,j,k,nnew,itrc)+cff2*tmp1
          ad_FX(i+1,j)=ad_FX(i+1,j)-tmp2
          ad_FX(i,j)=ad_FX(i,j)    +tmp2
          ad_FE(i,j+1)=ad_FE(i,j+1)-tmp2
          ad_FE(i,j)=ad_FE(i,j)    +tmp2
          ad_t(i,j,k,3,itrc)=0.d0  

         END DO
        END DO
!
!  ADJ to: Second-order, centered differences horizontal advective fluxes.
!

        DO j=Jstr,Jend+1
         DO i=Istr,Iend
!          tl_FE(i,j)=0.5d0*
!     *               (tl_Hvdx(i,j,k)*                       
!     *                (   t(i,j-1,k,nstp,itrc)+   t(i,j  ,k,nstp,itrc))+
!     *                Hvdx(i,j,k,nstp)*                       
!     *                (tl_t(i,j-1,k,nstp,itrc)+tl_t(i,j  ,k,nstp,itrc)))

          tmp1=0.5d0*ad_FE(i,j)
          tmp2=tmp1*Hvdx(i,j,k,nstp)

          ad_Hvdx(i,j,k)=ad_Hvdx(i,j,k)+
     *                 tmp1*(t(i,j-1,k,nstp,itrc)+t(i,j  ,k,nstp,itrc))         
          ad_t(i,j-1,k,nstp,itrc)=ad_t(i,j-1,k,nstp,itrc)+tmp2         
          ad_t(i,j  ,k,nstp,itrc)=ad_t(i,j  ,k,nstp,itrc)+tmp2    
          ad_FE(i,j)=0.d0

         END DO
        END DO

        DO j=Jstr,Jend
         DO i=Istr,Iend+1
          
! TL:          tl_FX(i,j)=0.5d0*
!     *              (tl_Hudy(i,j,k)*                                       
!     *               (   t(i-1,j,k,nstp,itrc)+   t(i  ,j,k,nstp,itrc))+
!     *               Hudy(i,j,k,nstp)*                                       
!     *               (tl_t(i-1,j,k,nstp,itrc)+tl_t(i  ,j,k,nstp,itrc))
!     *              )

          tmp1=0.5d0*ad_FX(i,j)
          tmp2=tmp1*Hudy(i,j,k,nstp)

          ad_Hudy(i,j,k)=ad_Hudy(i,j,k)+
     *                 tmp1*(t(i-1,j,k,nstp,itrc)+t(i  ,j,k,nstp,itrc))
          ad_t(i-1,j,k,nstp,itrc)=ad_t(i-1,j,k,nstp,itrc)+tmp2
          ad_t(i  ,j,k,nstp,itrc)=ad_t(i  ,j,k,nstp,itrc)+tmp2
          ad_FX(i,j)=0.d0

         END DO
        END DO

       END DO K_LOOP
      END DO T_LOOP1

#endif /* ifndef TS_FIXED */

      RETURN
      END

#endif /* ADJOINT */
      
