#include "cppdefs.h"

#ifdef ADJOINT 
      SUBROUTINE ad_step3d_t

      USE mod_ocean

      implicit none

      integer :: i, j, k, itrc,ltrc

      real(8) :: ad_FX(LBi:UBi,LBj:UBj)
      real(8) :: ad_FE(LBi:UBi,LBj:UBj)
      real(8) :: ad_DC(LBi:UBi,N)
      real(8) :: ad_tmp(LBi:UBi,N)
      real(8) :: BC(LBi:UBi,0:N)
      real(8) :: FC(LBi:UBi,0:N)
      real(8) :: ad_FC(LBi:UBi,0:N)

      real(8) :: cff,cff1,ad_cff1,tmp1,tmp2

      ad_cff1=0.d0
      ad_FX=0.d0
      ad_FE=0.d0
      ad_DC=0.d0
      ad_FC=0.d0
      ad_tmp=0.d0

!
!-----------------------------------------------------------------------
!  Apply lateral boundary conditions, mask
!-----------------------------------------------------------------------
!

      DO itrc=1,Ntr

!
!  Apply periodic boundary conditions.
!

#ifdef NS_PERIODIC
       CALL ad_exchange_3d(ad_t(:,:,:,nnew,itrc),
     *                      LBi,UBi,LBj,UBj,1,N,Lm,Mm)
#endif

!
!  Apply Land/Sea mask.
!
       DO k=1,N
        DO j=LBj,UBj
         DO i=LBi,UBi
          ad_t(i,j,k,nnew,itrc)=ad_t(i,j,k,nnew,itrc)*
     *                                              REAL(mask_r(i,j),8)
         END DO
        END DO
       END DO


!
!  Set lateral boundary conditions.
!
       CALL ad_t3dbc(Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj,N,Ntr,
#ifdef MPI
     *             B_West,B_East,B_South,B_North, ! <- new line, comp ot v2.1
#endif
     *            nrhs, nstp, nnew,itrc,ad_t)
             
      END DO ! itrc


!
! Pipeline j (vert. adv., -tl_Hz(nnew)*t(nnew), vert. diff.)
!
      DO j=Jstr,Jend

!
!-----------------------------------------------------------------------
!  ADJ: Time-step vertical diffusion term 
!-----------------------------------------------------------------------
! Note: No TL to Akt/Hz, so self-adj operator

       DO itrc=1,Ntr
        ltrc=MIN(2,itrc)

!
!  Compute off-diagonal (fwd) coefficients FC 
!

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
!  Compute (fwd) diagonal matrix coefficients BC:
!
!
! Assign ad_tmp(i,k)=ad_t(i,j,k,nnew,itrc)
!

        DO k=1,N
         DO i=Istr,Iend
          BC(i,k)=Hz(i,j,k)-FC(i,k)-FC(i,k-1)
          ad_tmp(i,k)=ad_t(i,j,k,nnew,itrc)
         END DO
        END DO

! After using ad_t(i,j,k,nnew,itrc) -> zero it out, 
! or not use additive stat. below

!
!  Solve the tridiagonal system 
! a(i,k)ad_DC(i,k-1) + BC(i,k)ad_DC(i,k) + c(k)ad_DC(i,k+1)= ad_tmp(i,k)
!  where a(k)=FC(k-1), c(k)=FC(k)
!

! - Forward elimination:

        DO k=2,N
         DO i=Istr,Iend
          cff1=FC(i,k-1)/BC(i,k-1)
          BC(i,k)=BC(i,k)-cff1*FC(i,k-1)
! tl (self-adj):          tl_DC(i,k)=tl_DC(i,k)-cff1*tl_DC(i,k-1)
          ad_tmp(i,k)=ad_tmp(i,k)-cff1*ad_tmp(i,k-1)
         END DO        
        END DO
       
! - Backward substitution:       

        DO i=Istr,Iend
! tl (self-adj):         tl_t(i,j,N,nnew,itrc)=tl_DC(i,N)/BC(i,N)
         ad_DC(i,N)=ad_tmp(i,N)/BC(i,N)
        END DO

        DO k=N-1,1,-1
         DO i=Istr,Iend
! tl (self-adj):          tl_t(i,j,k,nnew,itrc)=(tl_DC(i,k)-
!     *                       FC(i,k)*tl_t(i,j,k+1,nnew,itrc))/BC(i,k)

          ad_DC(i,k)=(ad_tmp(i,k)-FC(i,k)*ad_DC(i,k+1))/BC(i,k)

         END DO
        END DO

! ADJ TO assignment of the rhs in the vert. integration:
! (ad_t was not zeroed out above, do not use recurs. addition here)
        DO k=1,N
         DO i=Istr,Iend
! tl:          tl_DC(i,k)=tl_t(i,j,k,nnew,itrc)
          ad_t(i,j,k,nnew,itrc)=ad_DC(i,k)
          ad_DC(i,k)=0.d0
         END DO
        END DO

!---------------------------
! subtract tl_Hz(n+1) t(n+1) 
! (the part of tl tendency term that is known by now)
!---------------------------

        DO k=1,N
         DO i=Istr,Iend
! tl:          tl_t(i,j,k,nnew,itrc)=tl_t(i,j,k,nnew,itrc)-
!     *     tl_Hz(i,j,k)*t(i,j,k,nnew,itrc)

          ad_Hz(i,j,k)=ad_Hz(i,j,k)-
     *     t(i,j,k,nnew,itrc)*ad_t(i,j,k,nnew,itrc)

         END DO
        END DO

!
!-----------------------------------------------------------------------
!  Time-step vertical advection term.
!-----------------------------------------------------------------------
!

!
!  Time-step vertical advection term.
!
        DO k=1,N
         DO i=Istr,Iend
! tl line 2:          tl_t(i,j,k,nnew,itrc)=tl_t(i,j,k,nnew,itrc)-tl_cff1
          ad_cff1=ad_cff1-ad_t(i,j,k,nnew,itrc)       

! tl line 1:          tl_cff1=dt*oArea(i,j)*(tl_FC(i,k)-tl_FC(i,k-1))

          tmp1=ad_cff1*dt*oArea(i,j)
          ad_FC(i,k)=ad_FC(i,k)+tmp1
          ad_FC(i,k-1)=ad_FC(i,k-1)-tmp1
          ad_cff1=0.d0
         END DO
        END DO

!
!  Second-order, central differences vertical advective flux.
!
        DO i=Istr,Iend
         ad_FC(i,0)=0.0d0
         ad_FC(i,N)=0.0d0
        END DO

        DO k=1,N-1
         DO i=Istr,Iend
! tl:          tl_FC(i,k)=0.5d0*
!     *    ( tl_W(i,j,k)*(t(i,j,k  ,3,itrc)+t(i,j,k+1,3,itrc))+
!     *      W(i,j,k)*(tl_t(i,j,k  ,3,itrc)+tl_t(i,j,k+1,3,itrc)) )

          tmp1=ad_FC(i,k)*0.5d0
          tmp2=tmp1*W(i,j,k)

          ad_W(i,j,k)=ad_W(i,j,k)+
     *                tmp1*(t(i,j,k  ,3,itrc)+t(i,j,k+1,3,itrc))
          ad_t(i,j,k  ,3,itrc)=ad_t(i,j,k  ,3,itrc)+tmp2
          ad_t(i,j,k+1,3,itrc)=ad_t(i,j,k+1,3,itrc)+tmp2
          ad_FC(i,k)=0.d0

         END DO
        END DO
       
       END DO ! itrc=1,Ntr
      END DO  ! j

!
!-----------------------------------------------------------------------
!  ADJ to: Time-step horizontal advection term.
!-----------------------------------------------------------------------
!

      T_LOOP : DO itrc=1,Ntr
       K_LOOP : DO k=1,N
!
!  Time-step horizontal advection term.
!
        DO j=Jstr,Jend
         DO i=Istr,Iend

          cff=dt*oArea(i,j)

! TL:          tl_t(i,j,k,nnew,itrc)=tl_t(i,j,k,nnew,itrc)-
!     *          cff*(tl_FX(i+1,j)-tl_FX(i,j)+tl_FE(i,j+1)-tl_FE(i,j))

          tmp1=cff*ad_t(i,j,k,nnew,itrc)

          ad_FX(i+1,j)=ad_FX(i+1,j)-tmp1
          ad_FX(i,j)=ad_FX(i,j)    +tmp1
          ad_FE(i,j+1)=ad_FE(i,j+1)-tmp1
          ad_FE(i,j)=ad_FE(i,j)    +tmp1

         END DO
        END DO

!
!  Second-order, centered differences horizontal advective fluxes.
!

        DO j=JstrV-1,JendV+1
         DO i=IstrV,IendV

! TL:          tl_FE(i,j)=0.5d0*
!     *    ( tl_Hvdx(i,j,k)*(t(i,j-1,k,3,itrc)+t(i,j  ,k,3,itrc))+
!     *      Hvdx(i,j,k,3)*(tl_t(i,j-1,k,3,itrc)+tl_t(i,j  ,k,3,itrc)) )
 
          tmp1=0.5d0*ad_FE(i,j)
          tmp2=tmp1*Hvdx(i,j,k,3)

          ad_Hvdx(i,j,k)=ad_Hvdx(i,j,k)+
     *       tmp1*(t(i,j-1,k,3,itrc)+t(i,j  ,k,3,itrc))
          ad_t(i,j-1,k,3,itrc)=ad_t(i,j-1,k,3,itrc)+tmp2
          ad_t(i,j  ,k,3,itrc)=ad_t(i,j  ,k,3,itrc)+tmp2
          ad_FE(i,j)=0.d0
    
         END DO
        END DO

        DO j=JstrU,JendU
         DO i=IstrU-1,IendU+1
! TL:          tl_FX(i,j)=0.5d0*
!     *    ( tl_Hudy(i,j,k)*(t(i-1,j,k,3,itrc)+t(i  ,j,k,3,itrc))+
!     *      Hudy(i,j,k,3)*(tl_t(i-1,j,k,3,itrc)+tl_t(i  ,j,k,3,itrc)) )

          tmp1=0.5d0*ad_FX(i,j)
          tmp2=tmp1*Hudy(i,j,k,3)

          ad_Hudy(i,j,k)=ad_Hudy(i,j,k)+
     *       tmp1*(t(i-1,j,k,3,itrc)+t(i  ,j,k,3,itrc))
          ad_t(i-1,j,k,3,itrc)=ad_t(i-1,j,k,3,itrc)+tmp2
          ad_t(i  ,j,k,3,itrc)=ad_t(i  ,j,k,3,itrc)+tmp2
          ad_FX(i,j)=0.d0

         END DO
        END DO


       END DO K_LOOP
      END DO T_LOOP

      RETURN
      END
#endif /* ADJOINT */
