#include "cppdefs.h"

#ifdef ADJOINT 
      SUBROUTINE ad_step3d_uv

      USE mod_ocean
      USE mod_forces

      implicit none
      integer :: i,j,k
      real(8) :: cff, cff1, cff2,tmp,ad_tmp,ad_cff1

      real(8) :: AK(LBi:UBi,0:N)
      real(8) :: BC(LBi:UBi,0:N)
      real(8) :: CF(LBi:UBi,0:N)
      real(8) :: DC(LBi:UBi,0:N)
      real(8) :: FC(LBi:UBi,0:N)
      real(8) :: Hzk(LBi:UBi,N)

      real(8) :: ad_CF(LBi:UBi,0:N)
      real(8) :: ad_DC(LBi:UBi,0:N)
      real(8) :: ad_FC(LBi:UBi,0:N)
      real(8) :: ad_Hzk(LBi:UBi,N)

      ad_cff1=0.d0
      ad_tmp=0.d0
      ad_CF=0.d0
      ad_DC=0.d0
      ad_FC=0.d0
      ad_Hzk=0.d0

# if defined NS_PERIODIC
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!
      CALL ad_exchange_3d(ad_Hvdx,LBi,UBi,LBj,UBj,1,N,Lm,Mm)
      CALL ad_exchange_3d(ad_Hudy,LBi,UBi,LBj,UBj,1,N,Lm,Mm)
      CALL ad_exchange_3d(ad_v(:,:,:,nnew),LBi,UBi,LBj,UBj,1,N,Lm,Mm)
      CALL ad_exchange_3d(ad_u(:,:,:,nnew),LBi,UBi,LBj,UBj,1,N,Lm,Mm)
      CALL ad_exchange_2d(ad_vbar(:,:,1),LBi,UBi,LBj,UBj,Lm,Mm)
      CALL ad_exchange_2d(ad_ubar(:,:,1),LBi,UBi,LBj,UBj,Lm,Mm)
#endif

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
!  Couple velocity component in the ETA-direction.
!---
       IF (j.ge.Jstr) THEN

!
! Initialize FWD variables:
!

        DO i=IV_RANGE
         DC(i,0)=0.0d0
         CF(i,0)=0.0d0
         FC(i,0)=0.0d0
        END DO

!
!  FWD: Compute thicknesses of V-boxes DC(i,1:N) and total depth of the water
!  column DC(i,0)
!
        DO k=1,N
         DO i=IV_RANGE
          DC(i,k)=dx_v(i,j)*Hz_v(i,j,k)
          DC(i,0)=DC(i,0)+DC(i,k)
         END DO
        END DO

        DO i=IV_RANGE
         DC(i,0)=1.0d0/DC(i,0)
        END DO

!
!  ADJ TO: Compute mass flux for tracer computation, Hz*v*dx
!
        
! tl loop 3:
        DO k=1,N
         DO i=IV_RANGE

! NL:          Hvdx(i,j,k)=Hvdx(i,j,k)-DC(i,k)*FC(i,0)
! TL:          tl_Hvdx(i,j,k)=tl_Hvdx(i,j,k)-DC(i,k)*tl_FC(i,0) !-tl_DC(i,k)*FC(i,0)

          ad_FC(i,0)=ad_FC(i,0)-DC(i,k)*ad_Hvdx(i,j,k)

         END DO
        END DO

! tl loop 2:
        DO i=IV_RANGE

! NL:         FC(i,0)=DC(i,0)*tmp         !=0
! line 2:         tl_FC(i,0)=DC(i,0)*tl_tmp        !+tl_DC(i,0)*tmp, 2nd term=0 

         ad_tmp=ad_tmp+DC(i,0)*ad_FC(i,0)
         ad_FC(i,0)=0.d0
 
! NL:         tmp=FC(i,0)-DV_avg2(i,j)    ! by assumption, =0
! line 1:         tl_tmp=tl_FC(i,0)-tl_DV_avg2(i,j)

         ad_FC(i,0)=ad_FC(i,0)+ad_tmp
         ad_DV_avg2(i,j)=ad_DV_avg2(i,j)-ad_tmp
         ad_tmp=0.d0

        END DO

! tl loop 1:
        DO k=N,1,-1
         DO i=IV_RANGE

! NL:          FC(i,0)=FC(i,0)+Hvdx(i,j,k)
! line 2:          tl_FC(i,0)=tl_FC(i,0)+tl_Hvdx(i,j,k)

          ad_Hvdx(i,j,k)=ad_Hvdx(i,j,k)+ad_FC(i,0)

! NL:          Hvdx(i,j,k)=0.5d0*(Hvdx(i,j,k)+v(i,j,k,nnew)*DC(i,k))
! line 1:          tl_Hvdx(i,j,k)=0.5d0*(tl_Hvdx(i,j,k)+tl_v(i,j,k,nnew)*DC(i,k)
!      *                                        +v(i,j,k,nnew)*tl_DC(i,k))

          tmp=0.5d0*ad_Hvdx(i,j,k)
          ad_Hvdx(i,j,k)=tmp
          ad_v(i,j,k,nnew)=ad_v(i,j,k,nnew)+tmp*DC(i,k)
          ad_DC(i,k)=ad_DC(i,k)+tmp*v(i,j,k,nnew)

         END DO
        END DO

!
!  TL: Replace BOUNDARY POINTS incorrect vertical mean with more
!  accurate barotropic component, vbar=DV_avg1/(D*om_v). 
!

# ifndef NS_PERIODIC
! tl loop 4:
        IF (j.eq.Mm+1) THEN
         DO k=1,N
          DO i=IstrV,IendV
! tl line 2:         tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)*REAL(mask_v(i,j),8)
           ad_v(i,j,k,nnew)=ad_v(i,j,k,nnew)*REAL(mask_v(i,j),8)

! tl line 1 (recurs):  tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)-tl_CF(i,0)
           ad_CF(i,0)=ad_CF(i,0)-ad_v(i,j,k,nnew)

          END DO
         END DO
        END IF

! tl loop 3:
        IF (j.eq.1) THEN
         DO k=1,N
          DO i=IstrV,IendV
! tl line 2:           tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)*REAL(mask_v(i,j),8)
           ad_v(i,j,k,nnew)=ad_v(i,j,k,nnew)*REAL(mask_v(i,j),8)

! tl line 1 (recurs.): tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)-tl_CF(i,0)
           ad_CF(i,0)=ad_CF(i,0)-ad_v(i,j,k,nnew)

          END DO
         END DO
        END IF
# endif

        IF (B_East) THEN
! tl loop 2:
         DO k=1,N

! tl line 2:         tl_v(IendV+1,j,k,nnew)=tl_v(IendV+1,j,k,nnew)
!     *                         *REAL(mask_v(IendV+1,j),8)
          ad_v(IendV+1,j,k,nnew)=ad_v(IendV+1,j,k,nnew)
     *                          *REAL(mask_v(IendV+1,j),8)

! tl line 1 (recurs.): tl_v(IendV+1,j,k,nnew)=tl_v(IendV+1,j,k,nnew)-tl_CF(IendV+1,0)
          ad_CF(IendV+1,0)=ad_CF(IendV+1,0)-ad_v(IendV+1,j,k,nnew)

         END DO

        END IF

        IF (B_West) THEN
! tl loop 1:
         DO k=1,N

! tl line 2:         tl_v(IstrV-1,j,k,nnew)=tl_v(IstrV-1,j,k,nnew)
!     *                         *REAL(mask_v(IstrV-1,j),8)
          ad_v(IstrV-1,j,k,nnew)=ad_v(IstrV-1,j,k,nnew)
     *                         *REAL(mask_v(IstrV-1,j),8)

! tl line 1 (recurs): tl_v(IstrV-1,j,k,nnew)=tl_v(IstrV-1,j,k,nnew)-tl_CF(IstrV-1,0)
          ad_CF(IstrV-1,0)=ad_CF(IstrV-1,0)-ad_v(IstrV-1,j,k,nnew)

         END DO

        END IF


!
!  ADJ TO: Compute thicknesses of V-boxes DC(i,1:N), total depth of the water
!  column DC(i,0), and incorrect vertical mean CF(i,0).  
!

! tl loop 2:
        DO i=IV_RANGE

! tl line 4:  tl_vbar(i,j,1)=tl_DC(i,0)*DV_avg1(i,j)+DC(i,0)*tl_DV_avg1(i,j)
         ad_DC(i,0)=ad_DC(i,0)+ad_vbar(i,j,1)*DV_avg1(i,j)
         ad_DV_avg1(i,j)=ad_DV_avg1(i,j)+ad_vbar(i,j,1)*DC(i,0)
         ad_vbar(i,j,1)=0.d0         

! tl line 3:  tl_CF(i,0)=DC(i,0)*tl_cff1
         ad_cff1=ad_cff1+ad_CF(i,0)*DC(i,0)
         ad_CF(i,0)=0.d0         

! tl line 2:         tl_cff1=tl_CF(i,0)-tl_DV_avg1(i,j)
         ad_CF(i,0)=ad_CF(i,0)+ad_cff1
         ad_DV_avg1(i,j)=ad_DV_avg1(i,j)-ad_cff1
         ad_cff1=0.d0

! tl line 1 (recurs, scaling):  tl_DC(i,0)=-tl_DC(i,0)*DC(i,0)*DC(i,0)       
         ad_DC(i,0)=-ad_DC(i,0)*DC(i,0)*DC(i,0)

        END DO

! tl loop 1:
        DO k=1,N
         DO i=IV_RANGE
! tl 3:          tl_CF(i,0)=tl_CF(i,0)+tl_DC(i,k)*v(i,j,k,nnew)
! tl 3:     *                         +DC(i,k)*tl_v(i,j,k,nnew)
          ad_DC(i,k)=ad_DC(i,k)+ad_CF(i,0)*v(i,j,k,nnew)
          ad_v(i,j,k,nnew)=ad_v(i,j,k,nnew)+ad_CF(i,0)*DC(i,k)

! tl 2:          tl_DC(i,0)=tl_DC(i,0)+tl_DC(i,k)
          ad_DC(i,k)=ad_DC(i,k)+ad_DC(i,0)

! tl 1:          tl_DC(i,k)=dx_v(i,j)*tl_Hz_v(i,j,k) 
          ad_Hz_v(i,j,k)=ad_Hz_v(i,j,k)+dx_v(i,j)*ad_DC(i,k)
          ad_DC(i,k)=0.d0
 
         END DO
        END DO


!
! ADJ TO: Initialize some vars
!
        DO i=IV_RANGE
         ad_DC(i,0)=0.0d0
         ad_CF(i,0)=0.0d0
         ad_FC(i,0)=0.0d0
        END DO

       END IF ! j.ge.Jstr

!---
!  Couple velocity component in the XI-direction.
!---

       DO i=IU_RANGE
        DC(i,0)=0.0d0
        CF(i,0)=0.0d0
        FC(i,0)=0.0d0
       END DO

       DO k=1,N
        DO i=IU_RANGE
         DC(i,k)=dy_u(i,j)*Hz_u(i,j,k)
         DC(i,0)=DC(i,0)+DC(i,k)
        END DO
       END DO

       DO i=IU_RANGE
        DC(i,0)=1.0d0/DC(i,0)
       END DO

!
!  ADJ TO: Compute mass flux for tracer computation, Hz*u*dy; 
!  note DC(i,k)==dy_u(i,j)*Hz_u(i,j,k)
!

! tl loop 3:
       DO k=1,N
        DO i=IU_RANGE
! NL:         Hudy(i,j,k)=Hudy(i,j,k)-DC(i,k)*FC(i,0)
! TL:         tl_Hudy(i,j,k)=tl_Hudy(i,j,k)-DC(i,k)*tl_FC(i,0) !-tl_DC(i,k)*FC(i,0)

         ad_FC(i,0)=ad_FC(i,0)-DC(i,k)*ad_Hudy(i,j,k)      

        END DO
       END DO

! tl loop 2:
       DO i=IU_RANGE

!        FC(i,0)=DC(i,0)*tmp         !                =0
! tl line 2:        tl_FC(i,0)=DC(i,0)*tl_tmp    !+tl_DC(i,0)*tmp, 2nd term =0
        ad_tmp=ad_tmp+DC(i,0)*ad_FC(i,0)
        ad_FC(i,0)=0.d0

!         tmp=FC(i,0)-DU_avg2(i,j)   ! by assumption, =0
! tl line 1:        tl_tmp=tl_FC(i,0)-tl_DU_avg2(i,j)
        ad_FC(i,0)=ad_FC(i,0)+ad_tmp
        ad_DU_avg2(i,j)=ad_DU_avg2(i,j)-ad_tmp
        ad_tmp=0.d0

       END DO

! tl loop 1:
       DO k=N,1,-1
        DO i=IU_RANGE

! NL:        FC(i,0)=FC(i,0)+Hudy(i,j,k)
! tl line 2:         tl_FC(i,0)=tl_FC(i,0)+tl_Hudy(i,j,k)
         ad_Hudy(i,j,k)=ad_Hudy(i,j,k)+ad_FC(i,0)
         
! NL:       Hudy(i,j,k)=0.5d0*(Hudy(i,j,k)+u(i,j,k,nnew)*DC(i,k))  
! tl line 1:  tl_Hudy(i,j,k)=0.5d0*(tl_Hudy(i,j,k)+tl_u(i,j,k,nnew)*DC(i,k)
!     *                                      +u(i,j,k,nnew)*tl_DC(i,k))
         tmp=0.5d0*ad_Hudy(i,j,k)
         ad_Hudy(i,j,k)=tmp
         ad_u(i,j,k,nnew)=ad_u(i,j,k,nnew)+tmp*DC(i,k)
         ad_DC(i,k)=ad_DC(i,k)+tmp*u(i,j,k,nnew)
         
        END DO
       END DO

!
!  Replace BOUNDARY POINTS incorrect vertical mean with more
!  accurate barotropic component, ubar=DU_avg1/(D*on_u). 
!
# ifndef NS_PERIODIC

! tl loop 4:
       IF (j.eq.Mm+1) THEN
        DO k=1,N
         DO i=IstrU,IendU
! tl line 2:          tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)*REAL(mask_u(i,j),8)
          ad_u(i,j,k,nnew)=ad_u(i,j,k,nnew)*REAL(mask_u(i,j),8)

! tl line 1:          tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)-tl_CF(i,0)
          ad_CF(i,0)=ad_CF(i,0)-ad_u(i,j,k,nnew)

         END DO
        END DO
       END IF

! tl loop 3:
       IF (j.eq.0) THEN
        DO k=1,N
         DO i=IstrU,IendU

!          tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)*REAL(mask_u(i,j),8)
          ad_u(i,j,k,nnew)=ad_u(i,j,k,nnew)*REAL(mask_u(i,j),8)     

!          tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)-tl_CF(i,0)
          ad_CF(i,0)=ad_CF(i,0)-ad_u(i,j,k,nnew)

         END DO
        END DO
       END IF
# endif

       IF (B_East) THEN
! tl loop 2:
        DO k=1,N

!        tl_u(IendU+1,j,k,nnew)=tl_u(IendU+1,j,k,nnew)
!     *                       *REAL(mask_u(IendU+1,j),8)
         ad_u(IendU+1,j,k,nnew)=ad_u(IendU+1,j,k,nnew)
     *                       *REAL(mask_u(IendU+1,j),8)

!        tl_u(IendU+1,j,k,nnew)=tl_u(IendU+1,j,k,nnew)-tl_CF(IendU+1,0)
         ad_CF(IendU+1,0)=ad_CF(IendU+1,0)-ad_u(IendU+1,j,k,nnew)
        END DO
       END IF

       IF (B_West) THEN
! tl loop 1:
        DO k=1,N

!        tl_u(IstrU-1,j,k,nnew)=tl_u(IstrU-1,j,k,nnew)
!     *                       *REAL(mask_u(IstrU-1,j),8)
         ad_u(IstrU-1,j,k,nnew)=ad_u(IstrU-1,j,k,nnew)
     *                       *REAL(mask_u(IstrU-1,j),8)
        

!        tl_u(IstrU-1,j,k,nnew)=tl_u(IstrU-1,j,k,nnew)-tl_CF(IstrU-1,0)
         ad_CF(IstrU-1,0)=ad_CF(IstrU-1,0)-ad_u(IstrU-1,j,k,nnew)

        END DO
       END IF

!
!  Compute thicknesses of U-boxes DC(i,1:N), total depth of the water
!  column DC(i,0), and incorrect vertical mean CF(i,0).  Notice that
!  barotropic component is replaced with its fast-time averaged
!  values.
!

! tl loop 2:
       DO i=IU_RANGE

! tl line 4: tl_ubar(i,j,1)=tl_DC(i,0)*DU_avg1(i,j)+DC(i,0)*tl_DU_avg1(i,j)
        ad_DC(i,0)=ad_DC(i,0)+ad_ubar(i,j,1)*DU_avg1(i,j)
        ad_DU_avg1(i,j)=ad_DU_avg1(i,j)+ad_ubar(i,j,1)*DC(i,0)
        ad_ubar(i,j,1)=0.d0

! tl line 3: tl_CF(i,0)=DC(i,0)*tl_cff1  !+tl_DC(i,0)*cff1, 2nd term=0
        ad_cff1=ad_cff1+ad_CF(i,0)*DC(i,0)
        ad_CF(i,0)=0.d0
             
! tl line 2: tl_cff1=tl_CF(i,0)-tl_DU_avg1(i,j)
        ad_CF(i,0)=ad_CF(i,0)+ad_cff1
        ad_DU_avg1(i,j)=ad_DU_avg1(i,j)-ad_cff1
        ad_cff1=0.d0
 

! tl line 1: tl_DC(i,0)=-tl_DC(i,0)*DC(i,0)*DC(i,0) ! used below for tl_ubar
        ad_DC(i,0)=-ad_DC(i,0)*DC(i,0)*DC(i,0)

       END DO

! tl loop 1:
       DO k=1,N
        DO i=IU_RANGE
! 3:        tl_CF(i,0)=tl_CF(i,0)+tl_DC(i,k)*u(i,j,k,nnew)
!     *                        +DC(i,k)*tl_u(i,j,k,nnew)
         ad_DC(i,k)=ad_DC(i,k)+ad_CF(i,0)*u(i,j,k,nnew)
         ad_u(i,j,k,nnew)=ad_u(i,j,k,nnew)+ad_CF(i,0)*DC(i,k)

! 2:        tl_DC(i,0)=tl_DC(i,0)+tl_DC(i,k)
         ad_DC(i,k)=ad_DC(i,k)+ad_DC(i,0)

! 1:        tl_DC(i,k)=dy_u(i,j)*tl_Hz_u(i,j,k)
         ad_Hz_u(i,j,k)=ad_Hz_u(i,j,k)+dy_u(i,j)*ad_DC(i,k)
         ad_DC(i,k)=0.d0

        END DO
       END DO

       DO i=IU_RANGE
        ad_DC(i,0)=0.0d0
        ad_CF(i,0)=0.0d0
        ad_FC(i,0)=0.0d0
       END DO

      END DO  ! j=JU_RANGE

# undef IU_RANGE
# undef JU_RANGE
# undef IV_RANGE
# undef JV_RANGE

!
!-----------------------------------------------------------------------
! ADJ to Set lateral boundary conditions.
!-----------------------------------------------------------------------
!
      CALL ad_v3dbc(IstrV,IendV,JstrV,JendV,LBi,UBi,LBj,UBj,N,
#ifdef MPI
     *             B_West,B_East,B_South,B_North, ! <- new line, comp ot v2.1
#endif
     *             nrhs, nstp, nnew,ad_v)
      CALL ad_u3dbc(IstrU,IendU,JstrU,JendU,LBi,UBi,LBj,UBj,N,
#ifdef MPI
     *             B_West,B_East,B_South,B_North, ! <- new line, comp ot v2.1
#endif
     *             nrhs, nstp, nnew,ad_u)

!-----------------------------------------------------------------------
!  ADJ TO: Time step momentum equations 
!-----------------------------------------------------------------------
!
      DO j=JstrU,JendU     ! (pipelined j)
!
!-----------------------------------------------------------------------
!  Time step momentum equation in the ETA-direction.
!-----------------------------------------------------------------------
!
       IF (j.ge.JstrV) THEN

!
! FWD Variables: AK, Hzk
!

        DO i=Istr,Iend
         AK(i,0)=0.5d0*(Akv(i,j-1,0)+Akv(i,j  ,0))
         DO k=1,N
          AK(i,k)=0.5d0*(Akv(i,j-1,k)+Akv(i,j  ,k))
          Hzk(i,k)=Hz_v(i,j,k) 
         END DO
        END DO

! FWD variables, CF
        DO i=IstrV,IendV
         CF(i,0)=Hzk(i,1)
        END DO

        DO k=2,N
         DO i=Istr,Iend
          CF(i,0)=CF(i,0)+Hzk(i,k)
         END DO
        END DO

!
!  ADJ TO: Couple and update new solution.
!
        DO k=1,N
         DO i=Istr,Iend
! 2:          tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)*REAL(mask_v(i,j),8)
          ad_v(i,j,k,nnew)=ad_v(i,j,k,nnew)*REAL(mask_v(i,j),8)

! 1:          tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)-tl_DC(i,0)
          ad_DC(i,0)=ad_DC(i,0)-ad_v(i,j,k,nnew)
         END DO
        END DO

!
!  ADJ TO: Replace INTERIOR POINTS incorrect vertical mean with more accurate
!  barotropic component, vbar=DV_avg1/(D*om_v). Recall that, D=CF(:,0).
! (see notes above, where INTERIOR POINTS vert. mean corrected for u)
!

! tl loop 3:
        DO i=Istr,Iend
         cff1=1.0d0/(CF(i,0)*dx_v(i,j))

! tl 2:         tl_DC(i,0)=tl_tmp*cff1
         ad_tmp=ad_tmp+cff1*ad_DC(i,0)
         ad_DC(i,0)=0.d0

! tl 1:         tl_tmp=tl_DC(i,0)*dx_v(i,j)-tl_DV_avg1(i,j)
         ad_DC(i,0)=ad_DC(i,0)+ad_tmp*dx_v(i,j)
         ad_DV_avg1(i,j)=ad_DV_avg1(i,j)-ad_tmp
         ad_tmp=0.d0
  
        END DO

! tl loop 2:
        DO k=2,N
         DO i=Istr,Iend

! tl 2:          tl_DC(i,0)=tl_DC(i,0)+tl_v(i,j,k,nnew)*   Hzk(i,k)+
! tl 2:     *                             v(i,j,k,nnew)*tl_Hzk(i,k)
          ad_v(i,j,k,nnew)=ad_v(i,j,k,nnew)+ad_DC(i,0)*Hzk(i,k)
          ad_Hzk(i,k)=ad_Hzk(i,k)+ad_DC(i,0)*v(i,j,k,nnew)

! tl 1:          tl_CF(i,0)=tl_CF(i,0)+tl_Hzk(i,k)
          ad_Hzk(i,k)=ad_Hzk(i,k)+ad_CF(i,0)
 
         END DO
        END DO

! tl loop 1:
        DO i=IstrV,IendV

! tl 2:         tl_DC(i,0)=tl_v(i,j,1,nnew)*Hzk(i,1)+v(i,j,1,nnew)*tl_Hzk(i,1)
         ad_v(i,j,1,nnew)=ad_v(i,j,1,nnew)+ad_DC(i,0)*Hzk(i,1)
         ad_Hzk(i,1)=ad_Hzk(i,1)+ad_DC(i,0)*v(i,j,1,nnew)
         ad_DC(i,0)=0.d0

! tl 1:         tl_CF(i,0)=tl_Hzk(i,1)
         ad_Hzk(i,1)=ad_Hzk(i,1)+ad_CF(i,0)
         ad_CF(i,0)=0.d0

        END DO

!
!  ADJ TO: Time step right-hand-side terms.
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
!  Solve the tridiagonal system (SELF-ADJOINT):
!  a(i,k)ad_DC(i,k-1) + BC(i,k)ad_DC(i,k) + c(k)ad_DC(i,k+1)=ad_v(i,j,k,nnew)
!  where a(k)=FC(k-1), c(k)=FC(k)
!

! Note: ad_DC(i,k) must be 0 at this point

! - FWD Coefficients:
        DO k=1,N
         DO i=IstrV,IendV
          BC(i,k)=Hzk(i,k)-FC(i,k)-FC(i,k-1)
         END DO
        END DO

! - Forward elimination:

        DO k=2,N
         DO i=IstrV,IendV
          cff1=FC(i,k-1)/BC(i,k-1)
          BC(i,k)=BC(i,k)-cff1*FC(i,k-1)
          ad_v(i,j,k,nnew)=ad_v(i,j,k,nnew)-cff1*ad_v(i,j,k-1,nnew)
         END DO        
        END DO
       
! - Backward substitution:       

        DO i=IstrV,IendV
         ad_DC(i,N)=ad_v(i,j,N,nnew)/BC(i,N)
        END DO

        DO k=N-1,1,-1
         DO i=IstrV,IendV
          ad_DC(i,k)=(ad_v(i,j,k,nnew)-FC(i,k)*ad_DC(i,k+1))
     *                     /BC(i,k)
         END DO
        END DO

! After using ad_v(i,j,k,nnew) -> zero it out, or not use additive stat.
! below

! ADJ TO assignment of the rhs in the vert. integration:
        DO k=1,N
         DO i=IstrV,IendV

! TL:          tl_DC(i,k)=tl_v(i,j,k,nnew)
          ad_v(i,j,k,nnew)=ad_DC(i,k)
          ad_DC(i,k)=0.d0

         END DO
        END DO

! ADJ TO tl_v initialization:
! 
        DO k=1,N
         DO i=Istr,Iend
! TL:          tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)+DC(i,0)*tl_rv(i,j,k,nrhs)
!     *           -tl_Hzk(i,k)*v(i,j,k,nnew)

          ad_rv(i,j,k,nrhs)=ad_rv(i,j,k,nrhs)+DC(i,0)*ad_v(i,j,k,nnew)
          ad_Hzk(i,k)=ad_Hzk(i,k)-v(i,j,k,nnew)*ad_v(i,j,k,nnew)   

         END DO
        END DO

        DO i=Istr,Iend
         DO k=1,N
! TL:          tl_Hzk(i,k)=tl_Hz_v(i,j,k) 
          ad_Hz_v(i,j,k)=ad_Hz_v(i,j,k)+ad_Hzk(i,k)
          ad_Hzk(i,k)=0.d0
         END DO
        END DO


       END IF ! (j.ge.JstrV)

!-----------------------------------------------------------------------
!  Time step momentum equation in the XI-direction.
!-----------------------------------------------------------------------
!

!
! FWD Variables: AK, Hzk
!
       DO i=IstrU,IendU
        AK(i,0)=0.5d0*(Akv(i-1,j,0)+Akv(i  ,j,0))
        DO k=1,N
         AK(i,k)=0.5d0*(Akv(i-1,j,k)+Akv(i  ,j,k))
         Hzk(i,k)=Hz_u(i,j,k)  !0.5d0*(Hz(i-1,j,k)+Hz(i  ,j,k))
        END DO
       END DO

! FWD variables, CF
       DO i=IstrU,Iend
        CF(i,0)=Hzk(i,1)
       END DO

       DO k=2,N
        DO i=IstrU,Iend
         CF(i,0)=CF(i,0)+Hzk(i,k)
        END DO
       END DO

!
!  ADJ TO: Couple and update new solution.
!
       DO k=1,N
        DO i=IstrU,Iend

! tl 2:         tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)*REAL(mask_u(i,j),8)
         ad_u(i,j,k,nnew)=ad_u(i,j,k,nnew)*REAL(mask_u(i,j),8)
         
! tl 1:         tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)-tl_DC(i,0)
         ad_DC(i,0)=ad_DC(i,0)-ad_u(i,j,k,nnew) 

        END DO
       END DO

!
!  ADJ TO: Replace INTERIOR POINTS incorrect vertical mean with more accurate
!  barotropic component, ubar=DU_avg1/(D*on_u). Recall that, D=CF(:,0).
!

! tl loop 3:  
       DO i=IstrU,Iend
        cff1=1.0d0/(CF(i,0)*dy_u(i,j))

! tl 2:        tl_DC(i,0)=tl_tmp*cff1 ! +tmp*tl_cff1, second term =0
        ad_tmp=ad_tmp+cff1*ad_DC(i,0)
        ad_DC(i,0)=0.d0

! tl 1:        tl_tmp=tl_DC(i,0)*dy_u(i,j)-tl_DU_avg1(i,j)
        ad_DC(i,0)=ad_DC(i,0)+ad_tmp*dy_u(i,j)
        ad_DU_avg1(i,j)=ad_DU_avg1(i,j)-ad_tmp
        ad_tmp=0.d0
 
       END DO

! tl loop 2:  
       DO k=2,N
        DO i=IstrU,Iend
! tl 2:         tl_DC(i,0)=tl_DC(i,0)+tl_u(i,j,k,nnew)*Hzk(i,k)+
! tl 2:     *                         u(i,j,k,nnew)*tl_Hzk(i,k)
         ad_u(i,j,k,nnew)=ad_u(i,j,k,nnew)+ad_DC(i,0)*Hzk(i,k)
         ad_Hzk(i,k)=ad_Hzk(i,k)+u(i,j,k,nnew)*ad_DC(i,0)

! tl 1:         tl_CF(i,0)=tl_CF(i,0)+tl_Hzk(i,k)
         ad_Hzk(i,k)=ad_Hzk(i,k)+ad_CF(i,0)

        END DO
       END DO

! tl loop 1:  
       DO i=IstrU,Iend
! tl 2:        tl_DC(i,0)=tl_u(i,j,1,nnew)*Hzk(i,1)+u(i,j,1,nnew)*tl_Hzk(i,1)
        ad_u(i,j,1,nnew)=ad_u(i,j,1,nnew)+ad_DC(i,0)*Hzk(i,1)        
        ad_Hzk(i,1)=ad_Hzk(i,1)          +ad_DC(i,0)*u(i,j,1,nnew)
        ad_DC(i,0)=0.d0

! tl 1:        tl_CF(i,0)=tl_Hzk(i,1)
        ad_Hzk(i,1)=ad_Hzk(i,1)+ad_CF(i,0)
        ad_CF(i,0)=0.d0
       END DO

!
!  ADJ TO: Time step right-hand-side terms.
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
!  Solve the tridiagonal system (SELF-ADJOINT): 
!  a(i,k)ad_DC(i,k-1) + BC(i,k)ad_DC(i,k) + c(k)ad_DC(i,k+1)=ad_u(i,j,k,nnew)
!  where a(k)=FC(k-1), c(k)=FC(k)
!

! - FWD Coefficients:
       DO k=1,N
        DO i=IstrU,IendU
         BC(i,k)=Hzk(i,k)-FC(i,k)-FC(i,k-1)
        END DO
       END DO

! - Forward elimination:

        DO k=2,N
        DO i=IstrU,IendU
          cff1=FC(i,k-1)/BC(i,k-1)
          BC(i,k)=BC(i,k)-cff1*FC(i,k-1)
          ad_u(i,j,k,nnew)=ad_u(i,j,k,nnew)-cff1*ad_u(i,j,k-1,nnew)
         END DO        
        END DO
       
! - Backward substitution:       

        DO i=IstrU,IendU
         ad_DC(i,N)=ad_u(i,j,N,nnew)/BC(i,N)
        END DO

        DO k=N-1,1,-1
         DO i=IstrU,IendU
          ad_DC(i,k)=(ad_u(i,j,k,nnew)-FC(i,k)*ad_DC(i,k+1))
     *                     /BC(i,k)
         END DO
        END DO

! After using ad_u(i,j,k,nnew) -> zero it out, or not use additive statement
! below

! ADJ TO assignment of the rhs in the vert. integration:
       DO k=1,N
        DO i=IstrU,IendU
! tl:         tl_DC(i,k)=tl_u(i,j,k,nnew)
         ad_u(i,j,k,nnew)=ad_DC(i,k)
         ad_DC(i,k)=0.d0
        END DO
       END DO

! ADJ TO tl_u initialization:
! 
       DO k=1,N
        DO i=IstrU,IendU
! TL:         tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)+DC(i,0)*tl_ru(i,j,k,nrhs)
! TL:     *        -tl_Hzk(i,k)*u(i,j,k,nnew)
         ad_ru(i,j,k,nrhs)=ad_ru(i,j,k,nrhs)+DC(i,0)*ad_u(i,j,k,nnew)
         ad_Hzk(i,k)=ad_Hzk(i,k)-u(i,j,k,nnew)*ad_u(i,j,k,nnew)
        END DO
       END DO

       DO i=IstrU,IendU
        DO k=1,N
! TL:         tl_Hzk(i,k)=tl_Hz_u(i,j,k)
         ad_Hz_u(i,j,k)=ad_Hz_u(i,j,k)+ad_Hzk(i,k)
         ad_Hzk(i,k)=0.d0
        END DO
       END DO

      END DO ! j=JstrU,JendU (pipelined j)
      
      RETURN
      END

#endif /* ADJOINT */

