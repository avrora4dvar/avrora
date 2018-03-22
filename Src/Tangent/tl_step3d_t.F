#include "cppdefs.h"

#ifdef TANGENT 
      SUBROUTINE tl_step3d_t

      USE mod_ocean
#ifdef MPI
      USE mod_parallel
#endif
      implicit none

      integer :: i, j, k, itrc,ltrc

      real(8) :: tl_FX(LBi:UBi,LBj:UBj)
      real(8) :: tl_FE(LBi:UBi,LBj:UBj)
      real(8) :: tl_DC(LBi:UBi,N)
      real(8) :: BC(LBi:UBi,0:N)
      real(8) :: FC(LBi:UBi,0:N)
      real(8) :: tl_FC(LBi:UBi,0:N)

      real(8) :: cff,cff1,tl_cff1

      character(80) :: fhead

!
!-----------------------------------------------------------------------
!  Time-step horizontal advection term.
!-----------------------------------------------------------------------
!

!
!  Compute horizontal tracer advection fluxes.
!
      T_LOOP : DO itrc=1,Ntr
       K_LOOP : DO k=1,N
!
!  Second-order, centered differences horizontal advective fluxes.
!
        DO j=JstrU,JendU
         DO i=IstrU-1,IendU+1
! NL:          FX(i,j)=Hudy(i,j,k)*
!     *                0.5d0*(t(i-1,j,k,3,itrc)+
!     *                       t(i  ,j,k,3,itrc))

          tl_FX(i,j)=0.5d0*
     *    ( tl_Hudy(i,j,k)*(t(i-1,j,k,3,itrc)+t(i  ,j,k,3,itrc))+
     *      Hudy(i,j,k,3)*(tl_t(i-1,j,k,3,itrc)+tl_t(i  ,j,k,3,itrc)) )

         END DO
        END DO

        DO j=JstrV-1,JendV+1
         DO i=IstrV,IendV
! NL:          FE(i,j)=Hvdx(i,j,k)*  
!     *                0.5d0*(t(i,j-1,k,3,itrc)+
!     *                       t(i,j  ,k,3,itrc))

          tl_FE(i,j)=0.5d0*
     *    ( tl_Hvdx(i,j,k)*(t(i,j-1,k,3,itrc)+t(i,j  ,k,3,itrc))+
     *      Hvdx(i,j,k,3)*(tl_t(i,j-1,k,3,itrc)+tl_t(i,j  ,k,3,itrc)) )
     

         END DO
        END DO
!
!  Time-step horizontal advection term.
!
        DO j=Jstr,Jend
         DO i=Istr,Iend
          cff=dt*oArea(i,j)
          tl_t(i,j,k,nnew,itrc)=tl_t(i,j,k,nnew,itrc)-
     *          cff*(tl_FX(i+1,j)-tl_FX(i,j)+tl_FE(i,j+1)-tl_FE(i,j))
         END DO
        END DO

       END DO K_LOOP
      END DO T_LOOP

!
! Pipeline j (vert. adv., -tl_Hz(nnew)*t(nnew), vert. diff.
!

      DO j=Jstr,Jend
!
!-----------------------------------------------------------------------
!  Time-step vertical advection term.
!-----------------------------------------------------------------------
!
       DO itrc=1,Ntr
!
!  Second-order, central differences vertical advective flux.
!
        DO k=1,N-1
         DO i=Istr,Iend
          tl_FC(i,k)=0.5d0*
     *    ( tl_W(i,j,k)*(t(i,j,k  ,3,itrc)+t(i,j,k+1,3,itrc))+
     *      W(i,j,k)*(tl_t(i,j,k  ,3,itrc)+tl_t(i,j,k+1,3,itrc)) )
         END DO
        END DO
          
        DO i=Istr,Iend
         tl_FC(i,0)=0.0d0
         tl_FC(i,N)=0.0d0
        END DO
!
!  Time-step vertical advection term.
!
        DO k=1,N
         DO i=Istr,Iend
          tl_cff1=dt*oArea(i,j)*(tl_FC(i,k)-tl_FC(i,k-1))
          tl_t(i,j,k,nnew,itrc)=tl_t(i,j,k,nnew,itrc)-tl_cff1
         END DO
        END DO

!---------------------------
! subtract tl_Hz(n+1) t(n+1) 
! (the part of tl tendency term that is known by now)
!---------------------------
        DO k=1,N
         DO i=Istr,Iend
          tl_t(i,j,k,nnew,itrc)=tl_t(i,j,k,nnew,itrc)-
     *     tl_Hz(i,j,k)*t(i,j,k,nnew,itrc)
         END DO
        END DO

       END DO  ! DO itrc=1,Ntr
!
!-----------------------------------------------------------------------
!  Time-step vertical diffusion term.
!-----------------------------------------------------------------------
!
       DO itrc=1,Ntr
        ltrc=MIN(2,itrc)
!
!  Compute off-diagonal coefficients FC [lambdadt*Akt/Hz] for the
!  implicit vertical diffusion terms at future time step, located
!  at horizontal RHO-points and vertical W-points.
!  Also set FC at the top and bottom levels.
!
!  Fully implicit scheme: lambda=1

! Note: at this time, no tl to Akt/Hz, hence no tl_FC, no tl_BC

        cff=-dt
        DO k=1,N-1
         DO i=Istr,Iend
          cff1=1.0d0/(z_r(i,j,k+1)-z_r(i,j,k))
          FC(i,k)=cff*cff1*Akt(i,j,k,ltrc)
         END DO
        END DO

        DO i=Istr,Iend
         FC(i,0)=0.0d0
         FC(i,N)=0.0d0
        END DO
!
!  Compute diagonal matrix coefficients BC and load right-hand-side
!  terms for the tracer equation into DC.
!
        DO k=1,N
         DO i=Istr,Iend
          BC(i,k)=Hz(i,j,k)-FC(i,k)-FC(i,k-1)
          tl_DC(i,k)=tl_t(i,j,k,nnew,itrc)
         END DO
        END DO
!
!  Solve the tridiagonal system 
! a(i,k)tl_t(i,k-1) + BC(i,k)tl_t(i,k) + c(k)tl_t(i,k+1)=tl_DC(i,k)
!  where a(k)=FC(k-1), c(k)=FC(k)
!

! - Forward elimination:

        DO k=2,N
         DO i=Istr,Iend
          cff1=FC(i,k-1)/BC(i,k-1)
          BC(i,k)=BC(i,k)-cff1*FC(i,k-1)
          tl_DC(i,k)=tl_DC(i,k)-cff1*tl_DC(i,k-1)
         END DO        
        END DO
       
! - Backward substitution:       

        DO i=Istr,Iend
         tl_t(i,j,N,nnew,itrc)=tl_DC(i,N)/BC(i,N)
        END DO

        DO k=N-1,1,-1
         DO i=Istr,Iend
          tl_t(i,j,k,nnew,itrc)=(tl_DC(i,k)-
     *                       FC(i,k)*tl_t(i,j,k+1,nnew,itrc))/BC(i,k)
         END DO
        END DO

       END DO ! itrc
      END DO  ! j

!
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions, mask
!-----------------------------------------------------------------------
!
      DO itrc=1,Ntr
!
!  Set lateral boundary conditions.
!
       CALL tl_t3dbc(Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj,N,Ntr,
#ifdef MPI
     *             B_West,B_East,B_South,B_North, ! <- new line, comp ot v2.1
#endif
     *            nrhs, nstp, nnew,itrc,tl_t)

!
!  Apply Land/Sea mask.
!
       DO k=1,N
        DO j=LBj,UBj
         DO i=LBi,UBi
          tl_t(i,j,k,nnew,itrc)=tl_t(i,j,k,nnew,itrc)*
     *                                              REAL(mask_r(i,j),8)
         END DO
        END DO
       END DO

#ifdef MPI
       CALL mpi_exchange_3d_2(tl_t(LBi:UBi,LBj:UBj,1:N,nnew,itrc),
     *                     LBi,UBi,LBj,UBj,N,
     *                     Istr_E,Iend_E,Jstr_E,Jend_E,
     *                     comm2d, stride_r3d_bt, stride_3d_lr, 
     *                     rankleft,rankright,ranktop,rankbottom) 
#endif


!
!  Apply periodic boundary conditions.
!

#ifdef NS_PERIODIC
       CALL exchange_3d(tl_t(:,:,:,nnew,itrc),LBi,UBi,LBj,UBj,1,N,Lm,Mm)
#endif

      END DO ! itrc

      RETURN
      END
#endif /* TANGENT */
   
