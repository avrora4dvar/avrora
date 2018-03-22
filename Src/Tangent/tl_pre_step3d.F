#include "cppdefs.h"

#ifdef TANGENT 

      SUBROUTINE tl_pre_step3d
      USE mod_ocean
      USE mod_forces
#ifdef MPI
      USE mod_parallel
#endif
  
      implicit none
      integer :: indx,i,j,k,itrc
      real(8) :: cff,cff1,cff2,tl_cff1,tl_cff2,tmp,tl_tmp

      real(8) :: DC(LBi:UBi,N)
      real(8) :: FC(LBi:UBi,0:N)

      real(8) :: tl_FX(LBi:UBi,LBj:UBj)
      real(8) :: tl_FE(LBi:UBi,LBj:UBj)
      real(8) :: tl_DC(LBi:UBi,N)
      real(8) :: tl_FC(LBi:UBi,0:N)

      real(8), parameter :: Gamma = 1.0d0/6.0d0

      indx=3-nrhs

#ifndef TS_FIXED 
!=======================================================================
!  Tracer equation(s).
!=======================================================================
!

!
!-----------------------------------------------------------------------
!  Compute intermediate tracer at n+1/2 time-step, t(i,j,k,3,itrc).
!-----------------------------------------------------------------------
!
!  Compute time rate of change of intermediate tracer due to
!  horizontal advection.
!
      T_LOOP1 :DO itrc=1,Ntr
       K_LOOP: DO k=1,N
!
!  Second-order, centered differences horizontal advective fluxes.
!
        DO j=Jstr,Jend
         DO i=Istr,Iend+1
! NL:         FX(i,j)=Hudy(i,j,k)*                                       
!     *                0.5d0*(t(i-1,j,k,nstp,itrc)+                      
!     *                       t(i  ,j,k,nstp,itrc))
          
          tl_FX(i,j)=0.5d0*
     *              (tl_Hudy(i,j,k)*                                       
     *               (   t(i-1,j,k,nstp,itrc)+   t(i  ,j,k,nstp,itrc))+
     *               Hudy(i,j,k,nstp)*                                       
     *               (tl_t(i-1,j,k,nstp,itrc)+tl_t(i  ,j,k,nstp,itrc))
     *              )

         END DO
        END DO
          
        DO j=Jstr,Jend+1
         DO i=Istr,Iend
! NL:          FE(i,j)=Hvdx(i,j,k)*                       
!     *                0.5d0*(t(i,j-1,k,nstp,itrc)+
!     *                       t(i,j  ,k,nstp,itrc))

          tl_FE(i,j)=0.5d0*
     *               (tl_Hvdx(i,j,k)*                       
     *                (   t(i,j-1,k,nstp,itrc)+   t(i,j  ,k,nstp,itrc))+
     *                Hvdx(i,j,k,nstp)*                       
     *                (tl_t(i,j-1,k,nstp,itrc)+tl_t(i,j  ,k,nstp,itrc)))
         END DO
        END DO
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
! NL:         t(i,j,k,3,itrc)=Hz(i,j,k)*(cff1*t(i,j,k,nstp,itrc)+
!     *                                   cff2*t(i,j,k,nnew,itrc))-
!     *                        cff*oArea(i,j)*                
!     *                        (FX(i+1,j)-FX(i,j)+                 
!     *                         FE(i,j+1)-FE(i,j))

          tl_t(i,j,k,3,itrc)=
     *      tl_Hz(i,j,k)*
     *       (cff1*t(i,j,k,nstp,itrc)+cff2*t(i,j,k,nnew,itrc))+
     *      Hz(i,j,k)*
     *       (cff1*tl_t(i,j,k,nstp,itrc)+cff2*tl_t(i,j,k,nnew,itrc))-
     *      cff*oArea(i,j)*                
     *       (tl_FX(i+1,j)-tl_FX(i,j)+tl_FE(i,j+1)-tl_FE(i,j))

         END DO
        END DO

       END DO K_LOOP
      END DO T_LOOP1

!
! Pipelined J-loop
!

      J_LOOP1 : DO j=Jstr,Jend

!
!  Compute artificial continuity equation (same for all tracers) and
!  load it into private array DC (1/m).
!

       IF (iic.eq.1) THEN
        cff=0.5d0*dt
       ELSE
        cff=(1.0d0-Gamma)*dt
       END IF

       DO k=1,N
        DO i=Istr,Iend

! NL:         DC(i,k)=1.0d0/(Hz(i,j,k)-cff*oArea(i,j)*
!     *                      (Hudy(i+1,j,k)-Hudy(i,j,k)+
!     *                       Hvdx(i,j+1,k)-Hvdx(i,j,k)+
!     *                      (W(i,j,k)-W(i,j,k-1))))

! Rewrite NL:
         cff1=cff*oArea(i,j)
         tmp=Hz(i,j,k)-
     *       cff1*
     *        (Hudy(i+1,j,k,nstp)-Hudy(i,j,k,nstp)+
     *         Hvdx(i,j+1,k,nstp)-Hvdx(i,j,k,nstp)+
     *         W(i,j,k)-W(i,j,k-1) )

         DC(i,k)=1.d0/tmp

         tl_tmp=tl_Hz(i,j,k)-
     *          cff1*
     *           (tl_Hudy(i+1,j,k)-tl_Hudy(i,j,k)+
     *            tl_Hvdx(i,j+1,k)-tl_Hvdx(i,j,k)+
     *            tl_W(i,j,k)-tl_W(i,j,k-1) )

         tl_DC(i,k)=-tl_tmp*DC(i,k)*DC(i,k)
      
        END DO
       END DO
!
!  Compute time rate of change of intermediate tracer due to vertical
!  advection.  Impose artificial continuity equation.
!
       T_LOOP2: DO itrc=1,Ntr
!
!  Second-order, central differences vertical advective flux.
!
        DO k=1,N-1
         DO i=Istr,Iend
! NL:
          FC(i,k)=W(i,j,k)*   
     *                0.5d0*(t(i,j,k  ,nstp,itrc)+ 
     *                       t(i,j,k+1,nstp,itrc))
    
          tl_FC(i,k)= 0.5d0*
     *  ( tl_W(i,j,k)*(   t(i,j,k  ,nstp,itrc)+   t(i,j,k+1,nstp,itrc))+
     *      W(i,j,k)*(tl_t(i,j,k  ,nstp,itrc)+tl_t(i,j,k+1,nstp,itrc)) )
           

         END DO
        END DO
        DO i=Istr,Iend

         FC(i,0)=0.0d0
         FC(i,N)=0.0d0

         tl_FC(i,0)=0.0d0
         tl_FC(i,N)=0.0d0

        END DO
!
! Time-step vertical advection of tracers (Tunits).
!

        DO k=1,N
         DO i=Istr,Iend

          cff1=cff*oArea(i,j)

! NL:          t(i,j,k,3,itrc)=DC(i,k)*         
!     *                     (t(i,j,k,3,itrc)-
!     *                      cff1*(FC(i,k)-FC(i,k-1)))

          tl_t(i,j,k,3,itrc)=
     *     tl_DC(i,k)*
     *         ( t(i,j,k,3,itrc)-cff1*(FC(i,k)-FC(i,k-1)) ) +
     *     DC(i,k)*
     *         ( tl_t(i,j,k,3,itrc)-cff1*(tl_FC(i,k)-tl_FC(i,k-1)) )
        
         END DO
        END DO

       END DO T_LOOP2
      END DO J_LOOP1

!-----------------------------------------------------------------------
! Boundary conditions for the intermediate tracer:
!-----------------------------------------------------------------------
      DO itrc=1,Ntr
       CALL tl_t3dbc(Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj,N,Ntr,
#ifdef MPI
     *             B_West,B_East,B_South,B_North, ! <- new line, comp ot v2.1
#endif
     *            nrhs, nstp, 3,itrc,tl_t)


#ifdef MPI
      CALL mpi_exchange_3d_2(tl_t(LBi:UBi,LBj:UBj,1:N,3,itrc),
     *                     LBi,UBi,LBj,UBj,N,
     *                     Istr_E,Iend_E,Jstr_E,Jend_E,
     *                     comm2d, stride_r3d_bt, stride_3d_lr, 
     *                     rankleft,rankright,ranktop,rankbottom) 
#endif

#ifdef NS_PERIODIC
       CALL exchange_3d(tl_t(:,:,:,3,itrc),LBi,UBi,LBj,UBj,1,N,Lm,Mm)
#endif

      END DO

!
!-----------------------------------------------------------------------
!  Start computation of tracers at n+1 time-step, t(i,j,k,nnew,itrc).
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
       DO itrc=1,Ntr
        DO k=1,N
         DO i=Istr,Iend
          tl_t(i,j,k,nnew,itrc)=tl_Hz(i,j,k)*t(i,j,k,nstp,itrc)+
     *                          Hz(i,j,k)*tl_t(i,j,k,nstp,itrc)
         END DO
        END DO

        DO i=Istr,Iend
         tl_t(i,j,N,nnew,itrc)=tl_t(i,j,N,nnew,itrc)+
     *                          dt*tl_stflx(i,j,itrc)
        END DO

       END DO   ! DO itrc=1,Ntr
      END DO   ! DO j=Jstr,Jend

#endif /* ifndef TS_FIXED */

!=======================================================================
!  3D momentum equation in the XI-direction.
!=======================================================================
! note: oArea_u and oArea_v are scaled by mask_u and mask_v, corresp.

      IF (iic.eq.1) THEN

       DO k=1,N
        DO j=JstrU,JendU
         DO i=IstrU,IendU
! NL:           u(i,j,k,nnew)=u(i,j,k,nstp)*Hz_u(i,j,k)*DBLE(mask_u(i,j))
          tl_u(i,j,k,nnew)=(tl_u(i,j,k,nstp)*Hz_u(i,j,k)+
     *                      u(i,j,k,nstp)*tl_Hz_u(i,j,k))
     *                     *DBLE(mask_u(i,j))
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
          tl_u(i,j,k,nnew)=tl_u(i,j,k,nstp)*Hz_u(i,j,k)+
     *                      u(i,j,k,nstp)*tl_Hz_u(i,j,k)-
     *                       cff1*oArea_u(i,j)*tl_ru(i,j,k,indx)
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
! NL:           u(i,j,k,nnew)=u(i,j,k,nstp)*Hz_u(i,j,k)+
! NL:     *           oArea_u(i,j)*(cff1*ru(i,j,k,nrhs)-cff2*ru(i,j,k,indx))
          tl_u(i,j,k,nnew)=tl_u(i,j,k,nstp)*Hz_u(i,j,k)
     *                       +u(i,j,k,nstp)*tl_Hz_u(i,j,k)
     *     +oArea_u(i,j)*(cff1*tl_ru(i,j,k,nrhs)-cff2*tl_ru(i,j,k,indx))
         END DO
        END DO
       END DO
 
      END IF

! Add sustr and bustr to be used in the implicit vert. diss. term comput.
      DO j=JstrU,JendU
       DO i=IstrU,IendU
        cff=dt*DBLE(mask_u(i,j))
        tl_u(i,j,1,nnew)=tl_u(i,j,1,nnew)-cff*tl_bustr(i,j)
        tl_u(i,j,N,nnew)=tl_u(i,j,N,nnew)+cff*tl_sustr(i,j)
       END DO
      END DO

!=======================================================================
!  3D momentum equation in the ETA-direction.
!=======================================================================
!
      IF (iic.eq.1) THEN

       DO k=1,N
        DO j=JstrV,JendV
         DO i=IstrV,IendV
! NL:           v(i,j,k,nnew)=v(i,j,k,nstp)*Hz_v(i,j,k)*DBLE(mask_v(i,j))
           tl_v(i,j,k,nnew)=(tl_v(i,j,k,nstp)*Hz_v(i,j,k)
     *                +v(i,j,k,nstp)*tl_Hz_v(i,j,k) )*DBLE(mask_v(i,j))
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
           tl_v(i,j,k,nnew)=tl_v(i,j,k,nstp)*Hz_v(i,j,k)
     *                     +v(i,j,k,nstp)*tl_Hz_v(i,j,k)
     *                     -cff1*oArea_v(i,j)*tl_rv(i,j,k,indx)
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
           tl_v(i,j,k,nnew)=tl_v(i,j,k,nstp)*Hz_v(i,j,k)
     *       +v(i,j,k,nstp)*tl_Hz_v(i,j,k)
     *     +oArea_v(i,j)*(cff1*tl_rv(i,j,k,nrhs)-cff2*tl_rv(i,j,k,indx))
         END DO
        END DO
       END DO
 
      END IF  ! ETA-eqn., choices of iic==1, iic==2, or iic>2

! Add tl_svstr and tl_bvstr to be used in implicit vert. diss. term comput.
      DO j=JstrV,JendV
       DO i=IstrV,IendV
        cff=dt*DBLE(mask_v(i,j))
        tl_v(i,j,1,nnew)=tl_v(i,j,1,nnew)-cff*tl_bvstr(i,j)
        tl_v(i,j,N,nnew)=tl_v(i,j,N,nnew)+cff*tl_svstr(i,j)
       END DO
      END DO

      RETURN
      END

#endif /* TANGENT */
