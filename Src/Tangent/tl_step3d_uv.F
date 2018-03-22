#include "cppdefs.h"

#ifdef TANGENT 
      SUBROUTINE tl_step3d_uv

      USE mod_ocean
#ifdef MPI
      USE mod_parallel
#endif

      implicit none

      integer :: i, j, k

      real(8) :: cff, cff1, cff2, tmp, tl_tmp, AA, tl_AA,tl_cff1

      real(8) :: AK(LBi:UBi,0:N)
      real(8) :: BC(LBi:UBi,0:N)
      real(8) :: CF(LBi:UBi,0:N)
      real(8) :: DC(LBi:UBi,0:N)
      real(8) :: FC(LBi:UBi,0:N)
      real(8) :: Hzk(LBi:UBi,N)

      real(8) :: tl_CF(LBi:UBi,0:N)
      real(8) :: tl_DC(LBi:UBi,0:N)
      real(8) :: tl_FC(LBi:UBi,0:N)
      real(8) :: tl_Hzk(LBi:UBi,N)

!-----------------------------------------------------------------------
!  Time step momentum equations: 
!-----------------------------------------------------------------------
!
      DO j=JstrU,JendU

!-----------------------------------------------------------------------
!  Time step momentum equation in the XI-direction.
!-----------------------------------------------------------------------
!
       DO i=IstrU,IendU
        AK(i,0)=0.5d0*(Akv(i-1,j,0)+Akv(i  ,j,0))
        DO k=1,N
         AK(i,k)=0.5d0*(Akv(i-1,j,k)+Akv(i  ,j,k))
         Hzk(i,k)=Hz_u(i,j,k)  !0.5d0*(Hz(i-1,j,k)+Hz(i  ,j,k))

         tl_Hzk(i,k)=tl_Hz_u(i,j,k)

        END DO
       END DO
!
!  Time step right-hand-side terms.
!
       IF (iic.eq.1) THEN
        cff=dt
       ELSE IF (iic.eq.2) THEN
        cff=dt*3.0d0/2.0d0
       ELSE
        cff=dt*23.0d0/12.0d0
       END IF

       DO i=IstrU,IendU
        DC(i,0)=cff*oArea_u(i,j)
       END DO

! tl_u initialization in step3d_uv:
!
! Notes: (1) tl_Hz_u(i,j,k)*u(i,j,k,nnew) is subtracted here
! (as a result of differentiating of the tendency lhs term):
! (2) DC(i,0) is a constant 

       DO k=1,N
        DO i=IstrU,IendU
! NL:         u(i,j,k,nnew)=u(i,j,k,nnew)+DC(i,0)*ru(i,j,k,nrhs)
         tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)+DC(i,0)*tl_ru(i,j,k,nrhs)
     *        -tl_Hzk(i,k)*u(i,j,k,nnew)
        END DO
       END DO

!
!  Compute off-diagonal coefficients for the
!  implicit vertical viscosity term at horizontal U-points and
!  vertical W-points.
!
       cff=-2.d0*dt
       DO k=1,N-1
        DO i=IstrU,Iend
         FC(i,k)=cff*AK(i,k)/(Hzk(i,k+1)+Hzk(i,k))
        END DO
       END DO

       DO i=IstrU,Iend
        FC(i,0)=0.0d0
        FC(i,N)=0.0d0
       END DO
!
!  Solve the tridiagonal system 
! a(i,k)tl_u(i,k-1)+ BC(i,k)tl_u(i,k)+c(k)tl_u(i,k+1)=tl_DC(i,k)
!  where a(k)=FC(k-1), c(k)=FC(k)
! (no adj to AK)

! - Coefficients:
       DO k=1,N
        DO i=IstrU,IendU
         tl_DC(i,k)=tl_u(i,j,k,nnew)
         BC(i,k)=Hzk(i,k)-FC(i,k)-FC(i,k-1)
        END DO
       END DO

! - Forward elimination:

       DO k=2,N
        DO i=IstrU,IendU
         cff1=FC(i,k-1)/BC(i,k-1)
         BC(i,k)=BC(i,k)-cff1*FC(i,k-1)
         tl_DC(i,k)=tl_DC(i,k)-cff1*tl_DC(i,k-1)
        END DO        
       END DO
       
! - Backward substitution:       

       DO i=IstrU,IendU
        tl_u(i,j,N,nnew)=tl_DC(i,N)/BC(i,N)
       END DO

       DO k=N-1,1,-1
        DO i=IstrU,IendU
         tl_u(i,j,k,nnew)=(tl_DC(i,k)-FC(i,k)*tl_u(i,j,k+1,nnew))
     *                    /BC(i,k)
        END DO
       END DO

!
!  Replace INTERIOR POINTS incorrect vertical mean with more accurate
!  barotropic component, ubar=DU_avg1/(D*on_u). Recall that, D=CF(:,0).
!
!  Note (AK, 10/10/2007): Our choice of u, ubar, and Zt_avg1 will imply
!  that DC(i,0)=D*ubar=DU_avg1/dy, such that fwd tmp below is =0. 
!  Hence, DC(i,0) is not needed. We still compute tl perturbation
!  w/ respect to it.
! 
! Note (AK, 10/10/2007): CF(i,0) should be = D. However, leave it 
! like this, for consistency with the NL ROMS, NL Avrora.

! tl loop 1:  
       DO i=IstrU,Iend
        CF(i,0)=Hzk(i,1)
        tl_CF(i,0)=tl_Hzk(i,1)

! NL:        DC(i,0)=u(i,j,1,nnew)*Hzk(i,1)
        tl_DC(i,0)=tl_u(i,j,1,nnew)*Hzk(i,1)+u(i,j,1,nnew)*tl_Hzk(i,1)

       END DO

! tl loop 2:  
       DO k=2,N
        DO i=IstrU,Iend
         CF(i,0)=CF(i,0)+Hzk(i,k)
         tl_CF(i,0)=tl_CF(i,0)+tl_Hzk(i,k)

! NL:         DC(i,0)=DC(i,0)+u(i,j,k,nnew)*Hzk(i,k)
         tl_DC(i,0)=tl_DC(i,0)+tl_u(i,j,k,nnew)*Hzk(i,k)+
     *                         u(i,j,k,nnew)*tl_Hzk(i,k)

        END DO
       END DO

! tl loop 3:  
       DO i=IstrU,Iend
        cff1=1.0d0/(CF(i,0)*dy_u(i,j))
! Not used:        tl_cff1=-tl_CF(i,0)/(CF(i,0)*CF(i,0)*dy_u(i,j))

! NL:        tmp=DC(i,0)*dy_u(i,j)-DU_avg1(i,j)   ! = 0 (see note above)
        tl_tmp=tl_DC(i,0)*dy_u(i,j)-tl_DU_avg1(i,j)

! NL:        DC(i,0)=tmp*cff1 ! = 0 (see note above)
        tl_DC(i,0)=tl_tmp*cff1 ! +tmp*tl_cff1, second term =0

       END DO
!
!  Couple and update new solution.
!
       DO k=1,N
        DO i=IstrU,Iend
         tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)-tl_DC(i,0)
         tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)*REAL(mask_u(i,j),8)
        END DO
       END DO

!
!-----------------------------------------------------------------------
!  Time step momentum equation in the ETA-direction.
!-----------------------------------------------------------------------
!
       IF (j.ge.JstrV) THEN
        DO i=Istr,Iend
         AK(i,0)=0.5d0*(Akv(i,j-1,0)+Akv(i,j  ,0))
         DO k=1,N
          AK(i,k)=0.5d0*(Akv(i,j-1,k)+Akv(i,j  ,k))
          Hzk(i,k)=Hz_v(i,j,k) 

          tl_Hzk(i,k)=tl_Hz_v(i,j,k) 

         END DO
        END DO
!
!  Time step right-hand-side terms.
!
        IF (iic.eq.1) THEN
         cff=dt
        ELSE IF (iic.eq.2) THEN
         cff=dt*3.0d0/2.0d0
        ELSE
         cff=dt*23.0d0/12.0d0
        END IF

        DO i=Istr,Iend
         DC(i,0)=cff*oArea_v(i,j)
        END DO

! tl_v initialization:
! 
! Notes: (1) tl_Hz_v(i,j,k)*v(i,j,k,nnew) is subtracted here
! (as a result of differentiating of the tendency lhs term):
! (2) DC(i,0) is a constant 
        DO k=1,N
         DO i=Istr,Iend
          tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)+DC(i,0)*tl_rv(i,j,k,nrhs)
     *           -tl_Hzk(i,k)*v(i,j,k,nnew)
         END DO
        END DO

!
!  Compute off-diagonal coefficients [dt*Akv/Hz] for the
!  implicit vertical viscosity term at horizontal V-points and
!  vertical W-points.
!
        cff=-2.0d0*dt
        DO k=1,N-1
         DO i=Istr,Iend
          FC(i,k)=cff*AK(i,k)/(Hzk(i,k+1)+Hzk(i,k))
         END DO
        END DO

        DO i=Istr,Iend
         FC(i,0)=0.0d0
         FC(i,N)=0.0d0
        END DO
!
!  Solve the tridiagonal system:
!  a(i,k)tl_v(i,k-1) + BC(i,k)tl_v(i,k) + c(k)tl_v(i,k+1)=tl_DC(i,k)
!  where a(k)=FC(k-1), c(k)=FC(k)
!

! - Coefficients:
        DO k=1,N
         DO i=IstrV,IendV
          tl_DC(i,k)=tl_v(i,j,k,nnew)
          BC(i,k)=Hzk(i,k)-FC(i,k)-FC(i,k-1)
         END DO
        END DO

! - Forward elimination:

        DO k=2,N
         DO i=IstrV,IendV
          cff1=FC(i,k-1)/BC(i,k-1)
          BC(i,k)=BC(i,k)-cff1*FC(i,k-1)
          tl_DC(i,k)=tl_DC(i,k)-cff1*tl_DC(i,k-1)
         END DO        
        END DO
       
! - Backward substitution:       

        DO i=IstrV,IendV
         tl_v(i,j,N,nnew)=tl_DC(i,N)/BC(i,N)
        END DO

        DO k=N-1,1,-1
         DO i=IstrV,IendV
          tl_v(i,j,k,nnew)=(tl_DC(i,k)-FC(i,k)*tl_v(i,j,k+1,nnew))
     *                     /BC(i,k)
         END DO
        END DO

!
!  Replace INTERIOR POINTS incorrect vertical mean with more accurate
!  barotropic component, vbar=DV_avg1/(D*om_v). Recall that, D=CF(:,0).
! (see notes above, where INTERIOR POINTS vert. mean corrected for u)

! tl loop 1:
        DO i=IstrV,IendV
         CF(i,0)=Hzk(i,1)
         tl_CF(i,0)=tl_Hzk(i,1)

! NL:         DC(i,0)=v(i,j,1,nnew)*Hzk(i,1)
         tl_DC(i,0)=tl_v(i,j,1,nnew)*Hzk(i,1)+v(i,j,1,nnew)*tl_Hzk(i,1)
        END DO

! tl loop 2:
        DO k=2,N
         DO i=Istr,Iend
          CF(i,0)=CF(i,0)+Hzk(i,k)
          tl_CF(i,0)=tl_CF(i,0)+tl_Hzk(i,k)
 
! NL:          DC(i,0)=DC(i,0)+v(i,j,k,nnew)*Hzk(i,k)
          tl_DC(i,0)=tl_DC(i,0)+tl_v(i,j,k,nnew)*   Hzk(i,k)+
     *                             v(i,j,k,nnew)*tl_Hzk(i,k)
         END DO
        END DO

! tl loop 3:
        DO i=Istr,Iend
         cff1=1.0d0/(CF(i,0)*dx_v(i,j))
! not used:         tl_cff1=-tl_CF(i,0)/(CF(i,0)*CF(i,0)*dx_v(i,j))

! NL:         tmp=DC(i,0)*dx_v(i,j)-DV_avg1(i,j)      ! =0
         tl_tmp=tl_DC(i,0)*dx_v(i,j)-tl_DV_avg1(i,j)

! NL:         DC(i,0)=tmp*cff1                        ! =0
         tl_DC(i,0)=tl_tmp*cff1 ! +tmp*tl_cff1, 2nd term = 0

        END DO
!
!  Couple and update new solution.
!
        DO k=1,N
         DO i=Istr,Iend
          tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)-tl_DC(i,0)
          tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)*REAL(mask_v(i,j),8)
         END DO
        END DO

       END IF ! (j.ge.JstrV)
      END DO ! j=JstrU,JendU (pipelined j)

!
!-----------------------------------------------------------------------
! Set lateral boundary conditions.
!-----------------------------------------------------------------------
!
      CALL tl_u3dbc(IstrU,IendU,JstrU,JendU,LBi,UBi,LBj,UBj,N,
#ifdef MPI
     *             B_West,B_East,B_South,B_North, ! <- new line, comp ot v2.1
#endif
     *             nrhs, nstp, nnew,tl_u)
      CALL tl_v3dbc(IstrV,IendV,JstrV,JendV,LBi,UBi,LBj,UBj,N,
#ifdef MPI
     *             B_West,B_East,B_South,B_North, ! <- new line, comp ot v2.1
#endif
     *             nrhs, nstp, nnew,tl_v)


!
!-----------------------------------------------------------------------
!  Couple 2D and 3D momentum equations.
!-----------------------------------------------------------------------
!

# define IU_RANGE IstrU-1,IendU+1
# define IV_RANGE IstrV-1,IendV+1

# ifdef NS_PERIODIC
#  define JU_RANGE JstrU,JendU
#  define JV_RANGE JstrV,JendV
# else
#  define JU_RANGE JstrU-1,JendU+1
#  define JV_RANGE JstrV-1,JendV+1
# endif

      DO j=JU_RANGE      ! pipelined j

!---
!  Couple velocity component in the XI-direction.
!---
       DO i=IU_RANGE
        DC(i,0)=0.0d0
        CF(i,0)=0.0d0
        FC(i,0)=0.0d0

        tl_DC(i,0)=0.0d0
        tl_CF(i,0)=0.0d0
        tl_FC(i,0)=0.0d0
       END DO
!
!  Compute thicknesses of U-boxes DC(i,1:N), total depth of the water
!  column DC(i,0), and incorrect vertical mean CF(i,0).  Notice that
!  barotropic component is replaced with its fast-time averaged
!  values.
!

! tl loop 1:
       DO k=1,N
        DO i=IU_RANGE
         DC(i,k)=dy_u(i,j)*Hz_u(i,j,k)
         tl_DC(i,k)=dy_u(i,j)*tl_Hz_u(i,j,k)

         DC(i,0)=DC(i,0)+DC(i,k)
         tl_DC(i,0)=tl_DC(i,0)+tl_DC(i,k)

! NL:         CF(i,0)=CF(i,0)+DC(i,k)*u(i,j,k,nnew)
         tl_CF(i,0)=tl_CF(i,0)+tl_DC(i,k)*u(i,j,k,nnew)
     *                        +DC(i,k)*tl_u(i,j,k,nnew)
        END DO
       END DO

! tl loop 2:
       DO i=IU_RANGE
! Orig NL:
!        DC(i,0)=1.0d0/DC(i,0)                       
!        CF(i,0)=DC(i,0)*(CF(i,0)-DU_avg1(i,j))
!        ubar(i,j,1)=DC(i,0)*DU_avg1(i,j)

        DC(i,0)=1.0d0/DC(i,0)
        tl_DC(i,0)=-tl_DC(i,0)*DC(i,0)*DC(i,0) ! used below for tl_ubar

! Note, AK, 10/8: if fwd u, zeta, ubar are consistent,
! must obtain cff1=0

! NL:        cff1=CF(i,0)-DU_avg1(i,j)        ! =0
        tl_cff1=tl_CF(i,0)-tl_DU_avg1(i,j)

! NL:        CF(i,0)=DC(i,0)*cff1             ! =0
        tl_CF(i,0)=DC(i,0)*tl_cff1  !+tl_DC(i,0)*cff1, 2nd term=0
             
        tl_ubar(i,j,1)=tl_DC(i,0)*DU_avg1(i,j)+DC(i,0)*tl_DU_avg1(i,j)
       END DO

!
!  Replace BOUNDARY POINTS incorrect vertical mean with more
!  accurate barotropic component, ubar=DU_avg1/(D*on_u). 
!

       IF (B_West) THEN
! tl loop 1:
        DO k=1,N
         tl_u(IstrU-1,j,k,nnew)=tl_u(IstrU-1,j,k,nnew)-tl_CF(IstrU-1,0)
         tl_u(IstrU-1,j,k,nnew)=tl_u(IstrU-1,j,k,nnew)
     *                       *REAL(mask_u(IstrU-1,j),8)
        END DO
       END IF

       IF (B_East) THEN
! tl loop 2:
        DO k=1,N
         tl_u(IendU+1,j,k,nnew)=tl_u(IendU+1,j,k,nnew)-tl_CF(IendU+1,0)
         tl_u(IendU+1,j,k,nnew)=tl_u(IendU+1,j,k,nnew)
     *                       *REAL(mask_u(IendU+1,j),8)
        END DO
       END IF 

# ifndef NS_PERIODIC
! tl loop 3:
       IF (j.eq.0) THEN
        DO k=1,N
         DO i=IstrU,IendU
          tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)-tl_CF(i,0)
          tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)*REAL(mask_u(i,j),8)
         END DO
        END DO
       END IF

! tl loop 4:
       IF (j.eq.Mm+1) THEN
        DO k=1,N
         DO i=IstrU,IendU
          tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)-tl_CF(i,0)
          tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)*REAL(mask_u(i,j),8)
         END DO
        END DO
       END IF
# endif

!
!  Compute mass flux for tracer computation, Hz*u*dy; 
!  note DC(i,k)==dy_u(i,j)*Hz_u(i,j,k)
!

! tl loop 1:
       DO k=N,1,-1
        DO i=IU_RANGE
! NL:       Hudy(i,j,k)=0.5d0*(Hudy(i,j,k)+u(i,j,k,nnew)*DC(i,k))  
        tl_Hudy(i,j,k)=0.5d0*(tl_Hudy(i,j,k)+tl_u(i,j,k,nnew)*DC(i,k)
     *                                      +u(i,j,k,nnew)*tl_DC(i,k))

         
! NL:        FC(i,0)=FC(i,0)+Hudy(i,j,k)
         tl_FC(i,0)=tl_FC(i,0)+tl_Hudy(i,j,k)
         
        END DO
       END DO

! tl loop 2:
       DO i=IU_RANGE
!         tmp=FC(i,0)-DU_avg2(i,j)   ! by assumption, =0
        tl_tmp=tl_FC(i,0)-tl_DU_avg2(i,j)

!        FC(i,0)=DC(i,0)*tmp         !                =0
        tl_FC(i,0)=DC(i,0)*tl_tmp    !+tl_DC(i,0)*tmp, 2nd term =0
       END DO

! tl loop 3:
       DO k=1,N
        DO i=IU_RANGE
! NL:         Hudy(i,j,k)=Hudy(i,j,k)-DC(i,k)*FC(i,0)
         tl_Hudy(i,j,k)=tl_Hudy(i,j,k)-DC(i,k)*tl_FC(i,0) !-tl_DC(i,k)*FC(i,0)
        END DO
       END DO

!---
!  Couple velocity component in the ETA-direction.
!---
       IF (j.ge.Jstr) THEN

!
! Initialize some vars:
!
        DO i=IV_RANGE
         DC(i,0)=0.0d0
         CF(i,0)=0.0d0
         FC(i,0)=0.0d0

         tl_DC(i,0)=0.0d0
         tl_CF(i,0)=0.0d0
         tl_FC(i,0)=0.0d0
        END DO
!
!  Compute thicknesses of V-boxes DC(i,1:N), total depth of the water
!  column DC(i,0), and incorrect vertical mean CF(i,0).  Notice that
!  barotropic component is replaced with its fast-time averaged
!  values.
!

! tl loop 1:
        DO k=1,N
         DO i=IV_RANGE
          DC(i,k)=dx_v(i,j)*Hz_v(i,j,k)
          tl_DC(i,k)=dx_v(i,j)*tl_Hz_v(i,j,k) 
 
          DC(i,0)=DC(i,0)+DC(i,k)
          tl_DC(i,0)=tl_DC(i,0)+tl_DC(i,k)

! NL:          CF(i,0)=CF(i,0)+DC(i,k)*v(i,j,k,nnew)
          tl_CF(i,0)=tl_CF(i,0)+tl_DC(i,k)*v(i,j,k,nnew)
     *                         +DC(i,k)*tl_v(i,j,k,nnew)

         END DO
        END DO

! tl loop 2:
        DO i=IV_RANGE
         DC(i,0)=1.0d0/DC(i,0)
         tl_DC(i,0)=-tl_DC(i,0)*DC(i,0)*DC(i,0)       

! NL:         cff1=CF(i,0)-DV_avg1(i,j)    ! =0 w/ our choice of fwd
          tl_cff1=tl_CF(i,0)-tl_DV_avg1(i,j)

! NL:         CF(i,0)=DC(i,0)*cff1         ! =0
         tl_CF(i,0)=DC(i,0)*tl_cff1        !+tl_DC(i,0)*cff1, 2nd term=0 

         tl_vbar(i,j,1)=tl_DC(i,0)*DV_avg1(i,j)+DC(i,0)*tl_DV_avg1(i,j)

        END DO

!
!  TL: Replace BOUNDARY POINTS incorrect vertical mean with more
!  accurate barotropic component, vbar=DV_avg1/(D*om_v). 
!

        IF (B_West) THEN
! tl loop 1:
         DO k=1,N
          tl_v(IstrV-1,j,k,nnew)=tl_v(IstrV-1,j,k,nnew)-tl_CF(IstrV-1,0)
          tl_v(IstrV-1,j,k,nnew)=tl_v(IstrV-1,j,k,nnew)
     *                         *REAL(mask_v(IstrV-1,j),8)
         END DO
        END IF

        IF (B_East) THEN
! tl loop 2:
         DO k=1,N
          tl_v(IendV+1,j,k,nnew)=tl_v(IendV+1,j,k,nnew)-tl_CF(IendV+1,0)
          tl_v(IendV+1,j,k,nnew)=tl_v(IendV+1,j,k,nnew)
     *                         *REAL(mask_v(IendV+1,j),8)
         END DO
        END IF

# ifndef NS_PERIODIC
! tl loop 3:
        IF (j.eq.1) THEN
         DO k=1,N
          DO i=IstrV,IendV
           tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)-tl_CF(i,0)
           tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)*REAL(mask_v(i,j),8)
          END DO
         END DO
        END IF

! tl loop 4:
        IF (j.eq.Mm+1) THEN
         DO k=1,N
          DO i=IstrV,IendV
           tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)-tl_CF(i,0)
           tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)*REAL(mask_v(i,j),8)
          END DO
         END DO
        END IF
# endif
!
!  Compute mass flux for tracer computation, Hz*v*dx
!

! tl loop 1:
        DO k=N,1,-1
         DO i=IV_RANGE
! NL:          Hvdx(i,j,k)=0.5d0*(Hvdx(i,j,k)+v(i,j,k,nnew)*DC(i,k))
          tl_Hvdx(i,j,k)=0.5d0*(tl_Hvdx(i,j,k)+tl_v(i,j,k,nnew)*DC(i,k)
     *                                        +v(i,j,k,nnew)*tl_DC(i,k))

! NL:          FC(i,0)=FC(i,0)+Hvdx(i,j,k)
          tl_FC(i,0)=tl_FC(i,0)+tl_Hvdx(i,j,k)

         END DO
        END DO

! tl loop 2:
        DO i=IV_RANGE
! NL:         tmp=FC(i,0)-DV_avg2(i,j)    ! by assumption, =0
         tl_tmp=tl_FC(i,0)-tl_DV_avg2(i,j)

! NL:         FC(i,0)=DC(i,0)*tmp         !=0
         tl_FC(i,0)=DC(i,0)*tl_tmp        !+tl_DC(i,0)*tmp, 2nd term=0 
        END DO

! tl loop 3:
        DO k=1,N
         DO i=IV_RANGE
! NL:          Hvdx(i,j,k)=Hvdx(i,j,k)-DC(i,k)*FC(i,0)
          tl_Hvdx(i,j,k)=tl_Hvdx(i,j,k)-DC(i,k)*tl_FC(i,0) !-tl_DC(i,k)*FC(i,0)
         END DO
        END DO

       END IF ! j.ge.Jstr

      END DO  ! j=JU_RANGE

# undef IU_RANGE
# undef JU_RANGE
# undef IV_RANGE
# undef JV_RANGE

#ifdef MPI
! note mpi exchange 2d after tl_ubar, tl_vbar corrected (2d-3d coupling)
      CALL mpi_exchange_2d_2(tl_ubar(LBi:UBi,LBj:UBj,1),
     *                     LBi,UBi,LBj,UBj,
     *                     IstrU_E,IendU_E,JstrU_E,JendU_E,
     *                     comm2d, stride, 
     *                     rankleft,rankright,ranktop,rankbottom) 
      CALL mpi_exchange_2d_2(tl_vbar(LBi:UBi,LBj:UBj,1),
     *                     LBi,UBi,LBj,UBj,
     *                     IstrV_E,IendV_E,JstrV_E,JendV_E,
     *                     comm2d, stride, 
     *                     rankleft,rankright,ranktop,rankbottom) 
      CALL mpi_exchange_3d_2(tl_u(LBi:UBi,LBj:UBj,1:N,nnew),
     *                     LBi,UBi,LBj,UBj,N,
     *                     IstrU_E,IendU_E,JstrU_E,JendU_E,
     *                     comm2d, stride_u3d_bt, stride_3d_lr, 
     *                     rankleft,rankright,ranktop,rankbottom) 
      CALL mpi_exchange_3d_2(tl_v(LBi:UBi,LBj:UBj,1:N,nnew),
     *                     LBi,UBi,LBj,UBj,N,
     *                     IstrV_E,IendV_E,JstrV_E,JendV_E,
     *                     comm2d, stride_v3d_bt, stride_3d_lr, 
     *                     rankleft,rankright,ranktop,rankbottom) 
      CALL mpi_exchange_3d_2(tl_Hudy(LBi:UBi,LBj:UBj,1:N),
     *                     LBi,UBi,LBj,UBj,N,
     *                     IstrU_E,IendU_E,JstrU_E,JendU_E,
     *                     comm2d, stride_u3d_bt, stride_3d_lr, 
     *                     rankleft,rankright,ranktop,rankbottom) 
      CALL mpi_exchange_3d_2(tl_Hvdx(LBi:UBi,LBj:UBj,1:N),
     *                     LBi,UBi,LBj,UBj,N,
     *                     IstrV_E,IendV_E,JstrV_E,JendV_E,
     *                     comm2d, stride_v3d_bt, stride_3d_lr, 
     *                     rankleft,rankright,ranktop,rankbottom) 
#endif

# if defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      CALL exchange_2d(tl_ubar(:,:,1),LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(tl_vbar(:,:,1),LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_3d(tl_u(:,:,:,nnew),LBi,UBi,LBj,UBj,1,N,Lm,Mm)
      CALL exchange_3d(tl_v(:,:,:,nnew),LBi,UBi,LBj,UBj,1,N,Lm,Mm)
      CALL exchange_3d(tl_Hudy,LBi,UBi,LBj,UBj,1,N,Lm,Mm)
      CALL exchange_3d(tl_Hvdx,LBi,UBi,LBj,UBj,1,N,Lm,Mm)
#endif

      RETURN
      END

#endif /* TANGENT */
