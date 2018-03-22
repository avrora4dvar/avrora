#include "cppdefs.h"

#if defined TANGENT || defined ADJOINT
      SUBROUTINE set_state(t0,indx)

! A. Kurapov, 5/7/2007: obtains the background state by linear 
! interpolation b/ween snapshots zetaG, ubarG etc., at t=t0
! 3D vars: u(:,:,:,indx), 2D vars: ubar(:,:,1) 

      USE mod_ocean
 
      implicit none
      integer :: i,j,k,itrc,indx
      real(8) :: cff1,cff2,t0

      cff2=(t0-fwd_timeG(ifwd1))/(fwd_timeG(ifwd2)-fwd_timeG(ifwd1))
      cff1=1.d0-cff2

      DO j=LBj,UBj
       DO i=LBi,UBi
! NEW:
        fwd_zeta(i,j,indx)=cff1*zetaG(i,j,ifwd1)+cff2*zetaG(i,j,ifwd2)
        fwd_ubar(i,j,indx)=cff1*ubarG(i,j,ifwd1)+cff2*ubarG(i,j,ifwd2)
        fwd_vbar(i,j,indx)=cff1*vbarG(i,j,ifwd1)+cff2*vbarG(i,j,ifwd2)

       END DO
      END DO

      DO k=1,N
       DO j=LBj,UBj
        DO i=LBi,UBi
         u(i,j,k,indx)=cff1*uG(i,j,k,ifwd1)+cff2*uG(i,j,k,ifwd2)
         v(i,j,k,indx)=cff1*vG(i,j,k,ifwd1)+cff2*vG(i,j,k,ifwd2)
        END DO
       END DO
      END DO

      DO itrc=1,Ntr
       DO k=1,N
        DO j=LBj,UBj
         DO i=LBi,UBi
          t(i,j,k,indx,itrc)=cff1*tG(i,j,k,ifwd1,itrc)+
     *                       cff2*tG(i,j,k,ifwd2,itrc)
         END DO
        END DO
       END DO
      END DO

      RETURN
      END

!====================================
      SUBROUTINE set_intermediate_tracer
 
      USE mod_ocean
 
      implicit none
      integer :: i,j,k,itrc

      DO itrc=1,Ntr
       DO k=1,N
        DO j=LBj,UBj
         DO i=LBi,UBi
          t(i,j,k,3,itrc)=t(i,j,k,nstp,itrc)
         END DO
        END DO
       END DO
      END DO

      RETURN
      END


!====================================
      SUBROUTINE set_n_minus_1_tracer


! Sets background tracer at n-1 time level, to be used 
! in tl_pre_step3d. The n-1 level field is written in nnew. 
      USE mod_ocean
 
      implicit none
      integer :: i,j,k,itrc

      DO itrc=1,Ntr
       DO k=1,N
        DO j=LBj,UBj
         DO i=LBi,UBi
          t(i,j,k,nnew,itrc)=t(i,j,k,nstp,itrc)
         END DO
        END DO
       END DO
      END DO

      RETURN
      END


# ifndef ANA_FWD_MIX
!====================================
      SUBROUTINE set_fwd_mix(t0)

! A. Kurapov, 12/1/2007: obtains the background Akv, Akt by linear 
! interpolation b/ween snapshots at t=t0

! Note: w/ fully implicit vertical dissipation, both tl_avrora and ad_avrora 
! will use fwd Akv and Akt at time nnew. 

      USE mod_ocean
 
      implicit none
      integer :: i,j,k,itrc
      real(8) :: cff1,cff2,t0

      cff2=(t0-fwd_timeG(ifwd1))/(fwd_timeG(ifwd2)-fwd_timeG(ifwd1))
      cff1=1.d0-cff2

      DO k=0,N
       DO j=LBj,UBj
        DO i=LBi,UBi
         Akv(i,j,k)=cff1*AkvG(i,j,k,ifwd1)+cff2*AkvG(i,j,k,ifwd2)
         Akt(i,j,k,1)=cff1*AktG(i,j,k,ifwd1)+cff2*AktG(i,j,k,ifwd2)
         Akt(i,j,k,2)=Akt(i,j,k,1)
        END DO
       END DO
      END DO

      RETURN
      END

# endif /* ifndef ANA_FWD_MIX */

#endif
