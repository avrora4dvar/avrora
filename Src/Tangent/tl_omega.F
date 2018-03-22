#include "cppdefs.h"

# ifdef TANGENT
      SUBROUTINE tl_omega
      USE mod_ocean
#ifdef MPI
      USE mod_parallel
#endif
  
      implicit none
      integer :: i,j,k
      real(8) :: wrk(LBi:UBi),tl_wrk(LBi:UBi),tmp,tl_tmp

      integer :: io
      character (len=80) :: fhead

! W=Area*Hz*Omega, where Omega is the s-coord velocity

      DO j=Jstr,Jend

!-----------------------------------------------------------------------
!  Vertically integrate horizontal mass flux divergence.
!-----------------------------------------------------------------------

! TL loop 1:
        DO i=Istr,Iend
! NL:          W(i,j,0)=0.d0
          tl_W(i,j,0)=0.d0
        END DO

! TL loop 2:
        DO k=1,N
         DO i=Istr,Iend
! NL:          W(i,j,k)=W(i,j,k-1)- 
! NL:     *             (Hudy(i+1,j,k)-Hudy(i,j,k)+
! NL:     *              Hvdx(i,j+1,k)-Hvdx(i,j,k))
          tl_W(i,j,k)=tl_W(i,j,k-1)- 
     *             (tl_Hudy(i+1,j,k)-tl_Hudy(i,j,k)+
     *              tl_Hvdx(i,j+1,k)-tl_Hvdx(i,j,k))
         END DO
        END DO

! TL loop 3:
        DO i=Istr,Iend
         tmp=z_w(i,j,N)-z_w(i,j,0)
 
         wrk(i)=WN(i,j)/tmp
         tl_tmp=tl_z_w(i,j,N)-tl_z_w(i,j,0)
         tl_wrk(i)=(tl_W(i,j,N)-wrk(i)*tl_tmp)/tmp   ! diff in "flux" form
        END DO

!-----------------------------------------------------------------------
! Account for moving S-surfaces
!-----------------------------------------------------------------------
        DO k=N-1,1,-1
          DO i=Istr,Iend
! NL:            W(i,j,k)=W(i,j,k)-wrk(i)*(z_w(i,j,k)-z_w(i,j,0))
            tl_W(i,j,k)=tl_W(i,j,k)-tl_wrk(i)*(z_w(i,j,k)-z_w(i,j,0))
     *                             -wrk(i)*(tl_z_w(i,j,k)-tl_z_w(i,j,0))
          END DO
        END DO

        DO i=Istr,Iend
! NL:          W(i,j,N)=0.d0
          tl_W(i,j,N)=0.d0
        END DO

      END DO ! "j=..."


#ifdef MPI
! note: since tl_W(:,:,0)==0, no need to exchange. 
! use the rho-stride for tl_W(:,:,1:N)
      CALL mpi_exchange_3d_2(tl_W(LBi:UBi,LBj:UBj,1:N),
     *                     LBi,UBi,LBj,UBj,N,
     *                     Istr_E,Iend_E,Jstr_E,Jend_E,
     *                     comm2d, stride_r3d_bt, stride_3d_lr, 
     *                     rankleft,rankright,ranktop,rankbottom) 
#endif



#ifndef MPI
# if defined NS_PERIODIC
! NL:      CALL exchange_3d(W,LBi,UBi,LBj,UBj,0,N,Lm,Mm)
      CALL exchange_3d(tl_W,LBi,UBi,LBj,UBj,0,N,Lm,Mm)
# endif 
#endif

      RETURN
      END

#endif /* TANGENT */
