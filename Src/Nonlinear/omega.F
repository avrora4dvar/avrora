#include "cppdefs.h"

# if defined NONLINEAR
      SUBROUTINE omega
      USE mod_ocean
  
      implicit none
      integer :: i,j,k
      real(8) :: wrk(LBi:UBi)

! W=Area*Hz*Omega, where Omega is the s-coord velocity

      DO j=Jstr,Jend

!-----------------------------------------------------------------------
!  Vertically integrate horizontal mass flux divergence.
!-----------------------------------------------------------------------

        DO i=Istr,Iend
          W(i,j,0)=0.d0
        END DO

        DO k=1,N
          DO i=Istr,Iend
            W(i,j,k)=W(i,j,k-1)- 
     *               (Hudy(i+1,j,k)-Hudy(i,j,k)+
     *                Hvdx(i,j+1,k)-Hvdx(i,j,k))
          END DO
        END DO

        DO i=Istr,Iend
          wrk(i)=W(i,j,N)/(z_w(i,j,N)-z_w(i,j,0))
        END DO

!-----------------------------------------------------------------------
! Account for moving S-surfaces
!-----------------------------------------------------------------------
      
        DO k=N-1,1,-1
          DO i=Istr,Iend
            W(i,j,k)=W(i,j,k)-wrk(i)*(z_w(i,j,k)-z_w(i,j,0))
          END DO
        END DO

        DO i=Istr,Iend
          W(i,j,N)=0.d0
        END DO

      END DO

# if defined NS_PERIODIC
      CALL exchange_3d(W,LBi,UBi,LBj,UBj,0,N,Lm,Mm)
# endif 

      RETURN
      END

# elif defined TANGENT || defined ADJOINT
 
      SUBROUTINE fwd_omega(indx)
      USE mod_ocean
#ifdef MPI
      USE mod_parallel
#endif
  
      implicit none
      integer :: i,j,k,indx
      real(8) :: wrk(LBi:UBi)

! W=Area*Hz*Omega, where Omega is the s-coord velocity

      DO j=Jstr,Jend

!-----------------------------------------------------------------------
!  Vertically integrage horizontal mass flux divergence.
!-----------------------------------------------------------------------

        DO i=Istr,Iend
          W(i,j,0)=0.d0
        END DO

        DO k=1,N
          DO i=Istr,Iend
            W(i,j,k)=W(i,j,k-1)- 
     *               (Hudy(i+1,j,k,indx)-Hudy(i,j,k,indx)+
     *                Hvdx(i,j+1,k,indx)-Hvdx(i,j,k,indx))
          END DO
        END DO

! Note, save non-zero fwd WN(i,j) before correction for moving S-surf. is applied
        DO i=Istr,Iend
! NL:          wrk(i)=W(i,j,N)/(z_w(i,j,N)-z_w(i,j,0)), rewrite:
         WN(i,j)=W(i,j,N)                        
         wrk(i)=WN(i,j)/(z_w(i,j,N)-z_w(i,j,0))
        END DO

!-----------------------------------------------------------------------
! Account for moving S-surfaces
!-----------------------------------------------------------------------
      
        DO k=N-1,1,-1
          DO i=Istr,Iend
            W(i,j,k)=W(i,j,k)-wrk(i)*(z_w(i,j,k)-z_w(i,j,0))
          END DO
        END DO

        DO i=Istr,Iend
          W(i,j,N)=0.d0
        END DO

      END DO

#ifdef MPI
! note: since tl_W(:,:,0)==0, no need to exchange. 
! use the rho-stride for tl_W(:,:,1:N)

!      CALL mpi_exchange_2d(WN(LBi:UBi,LBj:UBj),
!     *                     LBi,UBi,LBj,UBj,
!     *                     Istr_E,Iend_E,Jstr_E,Jend_E,
!     *                     comm2d, stride, 
!     *                     rankleft,rankright,ranktop,rankbottom) 

      CALL mpi_exchange_3d_2(W(LBi:UBi,LBj:UBj,1:N),
     *                     LBi,UBi,LBj,UBj,N,
     *                     Istr_E,Iend_E,Jstr_E,Jend_E,
     *                     comm2d, stride_r3d_bt, stride_3d_lr, 
     *                     rankleft,rankright,ranktop,rankbottom) 
#endif

#ifndef MPI
# if defined NS_PERIODIC
      CALL exchange_3d(W,LBi,UBi,LBj,UBj,0,N,Lm,Mm)
# endif 
#endif

      RETURN
      END

# endif
