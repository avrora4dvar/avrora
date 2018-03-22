#include "cppdefs.h"

#ifdef TANGENT
      SUBROUTINE tl_uv_adv

      USE mod_ocean
#ifdef MPI
      USE mod_parallel
#endif

      implicit none
      integer :: i,j,k
      real(8) :: tl_cff
      real(8) :: tl_FC(LBi:UBi,0:N)

      character(len=80) :: fhead

!-----------------------------------------------------------------------
!  Horizontal advection of momentum:
!-----------------------------------------------------------------------

      K_LOOP : DO k=1,N

!
! Second-order, centered differences:
!

! U-flux of u momentum (at rho locations)
       DO j=JstrU,JendU
        DO i=IstrU-1,IendU
! NL:         UFx(i,j)=0.25d0*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))*
! NL:     *                   (Hudy(i,j,k,nrhs)+Hudy(i+1,j,k,nrhs))
         tl_UFx(i,j)=0.25d0*(
     *              (tl_u(i,j,k,nrhs)+tl_u(i+1,j,k,nrhs))*
     *              (Hudy(i,j,k,nrhs)+Hudy(i+1,j,k,nrhs)) + 
     *                    (u(i,j,k,nrhs)+u(i+1,j,k,nrhs))*
     *                    (tl_Hudy(i,j,k)+tl_Hudy(i+1,j,k)) )
        END DO
       END DO

! V-flux of u momentum (at psi locations):

       DO j=JstrU,JendU+1
        DO i=IstrU,IendU
! NL:         UFe(i,j)=0.25d0*(u(i,j-1,k,nrhs)+u(i,j  ,k,nrhs))* 
! NL:     *                   (Hvdx(i-1,j,k)+Hvdx(i  ,j,k))
         tl_UFe(i,j)=0.25d0*(
     *              (tl_u(i,j-1,k,nrhs)+tl_u(i,j  ,k,nrhs))*
     *              (Hvdx(i-1,j,k,nrhs)+Hvdx(i  ,j,k,nrhs))+
     *                       (u(i,j-1,k,nrhs)+u(i,j  ,k,nrhs))*
     *                       (tl_Hvdx(i-1,j,k)+tl_Hvdx(i  ,j,k)) )
        END DO
       END DO

! U-flux of v momentum (at psi locations):
       DO j=JstrV,JendV
        DO i=IstrV,IendV+1
! NL:         VFx(i,j)=0.25d0*(v(i-1,j,k,nrhs)+v(i  ,j,k,nrhs))*
! NL:     *                   (Hudy(i,j-1,k)+Hudy(i,j  ,k))
         tl_VFx(i,j)=0.25d0*( 
     *              (tl_v(i-1,j,k,nrhs)+tl_v(i  ,j,k,nrhs))*
     *              (Hudy(i,j-1,k,nrhs)+Hudy(i,j  ,k,nrhs))+
     *                       (v(i-1,j,k,nrhs)+v(i  ,j,k,nrhs))*
     *                       (tl_Hudy(i,j-1,k)+tl_Hudy(i,j  ,k)) )
        END DO
       END DO

! V-flux of v momentum (at rho locations):
       DO j=JstrV-1,JendV
        DO i=IstrV,IendV
! NL:         VFe(i,j)=0.25d0*(v(i,j  ,k,nrhs)+v(i,j+1,k,nrhs))*
! NL:     *                   (Hvdx(i,j  ,k)+Hvdx(i,j+1,k))
         tl_VFe(i,j)=0.25d0*( 
     *              (tl_v(i,j  ,k,nrhs)+tl_v(i,j+1,k,nrhs))*
     *              (Hvdx(i,j  ,k,nrhs)+Hvdx(i,j+1,k,nrhs))+
     *              (v(i,j  ,k,nrhs)+v(i,j+1,k,nrhs))*
     *              (tl_Hvdx(i,j  ,k)+tl_Hvdx(i,j+1,k)) )
        END DO
       END DO
       
!
! Add in horizontal advection:
!
       DO j=JstrU,JendU
        DO i=IstrU,IendU
         tl_cff=tl_UFx(i,j)-tl_UFx(i-1,j)+tl_UFe(i,j+1)-tl_UFe(i,j)
         tl_ru(i,j,k,nrhs)=tl_ru(i,j,k,nrhs)-tl_cff
        END DO
       END DO

       DO j=JstrV,JendV
        DO i=IstrV,IendV
         tl_cff=tl_VFx(i+1,j)-tl_VFx(i,j)+tl_VFe(i,j)-tl_VFe(i,j-1)
         tl_rv(i,j,k,nrhs)=tl_rv(i,j,k,nrhs)-tl_cff
        END DO
       END DO

      END DO K_LOOP

!-----------------------------------------------------------------------
!  Vertical advection of momentum:
!-----------------------------------------------------------------------

! TL vert mom, loop 1 ( U ):
      DO i=LBi,UBi
       tl_FC(i,0)=0.d0
       tl_FC(i,N)=0.d0
      END DO

! TL vert mom, loop 2 ( U ):
      DO j=JstrU,JendU

! TL vert mom, loop 2.1:
       DO k=1,N-1
        DO i=IstrU,IendU
! NL:         FC(i,k)=0.25d0*(u(i,j,k,nrhs)+u(i,j,k+1,nrhs))*
! NL:     *                   (W(i,j,k)+W(i-1,j,k))
         tl_FC(i,k)=0.25d0*(
     *     (tl_u(i,j,k,nrhs)+tl_u(i,j,k+1,nrhs))*(W(i,j,k)+W(i-1,j,k))+
     *     (u(i,j,k,nrhs)+u(i,j,k+1,nrhs))*(tl_W(i,j,k)+tl_W(i-1,j,k)) 
     *                     )
        END DO
       END DO

! TL vert mom, loop 2.2:
       DO k=1,N
        DO i=IstrU,IendU
         tl_cff=tl_FC(i,k)-tl_FC(i,k-1)
         tl_ru(i,j,k,nrhs)=tl_ru(i,j,k,nrhs)-tl_cff
        END DO
       END DO
       
      END DO ! j=JstrU,JendU, vert adv. (u)

! TL vert mom, loop 3 (V):
      DO j=JstrV,JendV

! TL vert mom, loop 3.1:
       DO k=1,N-1
        DO i=IstrV,IendV
! NL:         FC(i,k)=0.25d0*(v(i,j,k,nrhs)+v(i,j,k+1,nrhs))*
! NL:     *                  (W(i,j,k)+W(i,j-1,k))
         tl_FC(i,k)=0.25d0*(
     *     (tl_v(i,j,k,nrhs)+tl_v(i,j,k+1,nrhs))*(W(i,j,k)+W(i,j-1,k))+
     *     (v(i,j,k,nrhs)+v(i,j,k+1,nrhs))*(tl_W(i,j,k)+tl_W(i,j-1,k))
     *                     )  
        END DO
       END DO

! TL vert mom, loop 3.2:
       DO k=1,N
        DO i=IstrV,IendV
         tl_cff=tl_FC(i,k)-tl_FC(i,k-1)
         tl_rv(i,j,k,nrhs)=tl_rv(i,j,k,nrhs)-tl_cff
        END DO
       END DO

      END DO  ! "j=JstrV,JendV", vert adv. (v)

      RETURN
      END

#endif /* TANGENT */
