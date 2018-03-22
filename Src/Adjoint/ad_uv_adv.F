#include "cppdefs.h"

#ifdef ADJOINT
      SUBROUTINE ad_uv_adv

      USE mod_ocean

      implicit none
      integer :: i,j,k
      real(8) :: tmp,tmp1,tmp2
      real(8) :: ad_cff
      real(8) :: ad_FC(LBi:UBi,0:N)

      ad_cff=0.d0
      ad_FC=0.d0

!-----------------------------------------------------------------------
!  ADJ to vertical advection of momentum:
!-----------------------------------------------------------------------

! ADJ to TL vert mom, loop 3 (V):
      DO j=JstrV,JendV

! ADJ to TL vert mom, loop 3.2:
       DO k=1,N
        DO i=IstrV,IendV
! tl line 2 (recurs):         tl_rv(i,j,k,nrhs)=tl_rv(i,j,k,nrhs)-tl_cff
         ad_cff=ad_cff-ad_rv(i,j,k,nrhs)

! tl line 1:         tl_cff=tl_FC(i,k)-tl_FC(i,k-1)
         ad_FC(i,k)=ad_FC(i,k)+ad_cff
         ad_FC(i,k-1)=ad_FC(i,k-1)-ad_cff
         ad_cff=0.d0
        END DO
       END DO

! ADJ to TL vert mom, loop 3.1:
       DO k=1,N-1
        DO i=IstrV,IendV
! NL:         FC(i,k)=0.25d0*(v(i,j,k,nrhs)+v(i,j,k+1,nrhs))*
! NL:     *                  (W(i,j,k)+W(i,j-1,k))

! TL:         tl_FC(i,k)=0.25d0*(
! TL:     *     (tl_v(i,j,k,nrhs)+tl_v(i,j,k+1,nrhs))*(W(i,j,k)+W(i,j-1,k))+
! TL:     *     (v(i,j,k,nrhs)+v(i,j,k+1,nrhs))*(tl_W(i,j,k)+tl_W(i,j-1,k))
! TL:     *                     )  
 
         tmp=0.25d0*ad_FC(i,k)
         tmp1=tmp*(W(i,j,k)+W(i,j-1,k))
         tmp2=tmp*(v(i,j,k,nrhs)+v(i,j,k+1,nrhs))
 
         ad_v(i,j,k,nrhs)=ad_v(i,j,k,nrhs)+tmp1   
         ad_v(i,j,k+1,nrhs)=ad_v(i,j,k+1,nrhs)+tmp1       
         ad_W(i,j,k)=ad_W(i,j,k)+tmp2 
         ad_W(i,j-1,k)=ad_W(i,j-1,k)+tmp2       

         ad_FC(i,k)=0.d0

        END DO
       END DO

      END DO  ! "j=JstrV,JendV", vert adv. (v)

! ADJ to TL vert mom, loop 2 ( U ):
      DO j=JstrU,JendU
 
! ADJ to TL vert mom, loop 2.2:
       DO k=1,N
        DO i=IstrU,IendU
! tl line 2 (recur): tl_ru(i,j,k,nrhs)=tl_ru(i,j,k,nrhs)-tl_cff
         ad_cff=ad_cff-ad_ru(i,j,k,nrhs)

! tl line 1:         tl_cff=tl_FC(i,k)-tl_FC(i,k-1)
         ad_FC(i,k)=ad_FC(i,k)+ad_cff
         ad_FC(i,k-1)=ad_FC(i,k-1)-ad_cff
         ad_cff=0.d0
         
        END DO
       END DO

! ADJ to TL vert mom, loop 2.1:
       DO k=1,N-1
        DO i=IstrU,IendU
! NL:         FC(i,k)=0.25d0*(u(i,j,k,nrhs)+u(i,j,k+1,nrhs))*
! NL:     *                   (W(i,j,k)+W(i-1,j,k))

! TL:         tl_FC(i,k)=0.25d0*(
! TL:     *     (tl_u(i,j,k,nrhs)+tl_u(i,j,k+1,nrhs))*(W(i,j,k)+W(i-1,j,k))+
! TL:     *     (u(i,j,k,nrhs)+u(i,j,k+1,nrhs))*(tl_W(i,j,k)+tl_W(i-1,j,k)) 
! TL:     *                     )

         tmp=0.25d0*ad_FC(i,k)
         tmp1=tmp*(W(i,j,k)+W(i-1,j,k))
         tmp2=tmp*(u(i,j,k,nrhs)+u(i,j,k+1,nrhs))
         
         ad_u(i,j,k,nrhs)=ad_u(i,j,k,nrhs)+tmp1
         ad_u(i,j,k+1,nrhs)=ad_u(i,j,k+1,nrhs)+tmp1
         ad_W(i,j,k)=ad_W(i,j,k)+tmp2
         ad_W(i-1,j,k)=ad_W(i-1,j,k)+tmp2

         ad_FC(i,k)=0.d0

        END DO
       END DO

      END DO ! j=JstrU,JendU, vert adv. (u)

! TL vert mom, loop 1 ( U ):
      DO i=LBi,UBi
       ad_FC(i,0)=0.d0
       ad_FC(i,N)=0.d0
      END DO

!-----------------------------------------------------------------------
!  ADJ to horizontal advection of momentum:
!-----------------------------------------------------------------------

      K_LOOP : DO k=N,1,-1
!
! Add in horizontal advection:
!
       DO j=JstrV,JendV
        DO i=IstrV,IendV
! tl line 2 (recurs):   tl_rv(i,j,k,nrhs)=tl_rv(i,j,k,nrhs)-tl_cff
         ad_cff=ad_cff-ad_rv(i,j,k,nrhs)

! tl line 1:   tl_cff=tl_VFx(i+1,j)-tl_VFx(i,j)+tl_VFe(i,j)-tl_VFe(i,j-1)
         ad_VFx(i+1,j)=ad_VFx(i+1,j)+ad_cff     
         ad_VFx(i,j)=ad_VFx(i,j)    -ad_cff  
         ad_VFe(i,j)=ad_VFe(i,j)    +ad_cff  
         ad_VFe(i,j-1)=ad_VFe(i,j-1)-ad_cff      
        
         ad_cff=0.d0

        END DO
       END DO

       DO j=JstrU,JendU
        DO i=IstrU,IendU
! tl line 2 (recurs):  tl_ru(i,j,k,nrhs)=tl_ru(i,j,k,nrhs)-tl_cff
         ad_cff=ad_cff-ad_ru(i,j,k,nrhs)

! tl line 1:    tl_cff=tl_UFx(i,j)-tl_UFx(i-1,j)+tl_UFe(i,j+1)-tl_UFe(i,j)
         ad_UFx(i,j)=ad_UFx(i,j)+ad_cff
         ad_UFx(i-1,j)=ad_UFx(i-1,j)-ad_cff
         ad_UFe(i,j+1)=ad_UFe(i,j+1)+ad_cff
         ad_UFe(i,j)=ad_UFe(i,j)-ad_cff

         ad_cff=0.d0

        END DO
       END DO

!
! Second-order, centered differences:
!

! V-flux of v momentum (at rho locations):
       DO j=JstrV-1,JendV
        DO i=IstrV,IendV
! NL:         VFe(i,j)=0.25d0*(v(i,j  ,k,nrhs)+v(i,j+1,k,nrhs))*
! NL:     *                   (Hvdx(i,j  ,k)+Hvdx(i,j+1,k))

! TL:         tl_VFe(i,j)=0.25d0*( 
! TL:     *              (tl_v(i,j  ,k,nrhs)+tl_v(i,j+1,k,nrhs))*
! TL:     *              (Hvdx(i,j  ,k)+Hvdx(i,j+1,k))+
! TL:     *              (v(i,j  ,k,nrhs)+v(i,j+1,k,nrhs))*
! TL:     *              (tl_Hvdx(i,j  ,k)+tl_Hvdx(i,j+1,k)) )

         tmp=0.25d0*ad_VFe(i,j)
         tmp1=tmp*(Hvdx(i,j  ,k,nrhs)+Hvdx(i,j+1,k,nrhs))
         tmp2=tmp*(v(i,j  ,k,nrhs)+v(i,j+1,k,nrhs))

         ad_v(i,j  ,k,nrhs)=ad_v(i,j  ,k,nrhs)+tmp1
         ad_v(i,j+1,k,nrhs)=ad_v(i,j+1,k,nrhs)+tmp1
         ad_Hvdx(i,j  ,k)=ad_Hvdx(i,j  ,k)+tmp2
         ad_Hvdx(i,j+1,k)=ad_Hvdx(i,j+1,k)+tmp2

         ad_VFe(i,j)=0.d0

        END DO
       END DO

! U-flux of v momentum (at psi locations):
       DO j=JstrV,JendV
        DO i=IstrV,IendV+1
! NL:         VFx(i,j)=0.25d0*(v(i-1,j,k,nrhs)+v(i  ,j,k,nrhs))*
! NL:     *                   (Hudy(i,j-1,k)+Hudy(i,j  ,k))

! TL:         tl_VFx(i,j)=0.25d0*( 
! TL:     *              (tl_v(i-1,j,k,nrhs)+tl_v(i  ,j,k,nrhs))*
! TL:     *              (Hudy(i,j-1,k)+Hudy(i,j  ,k))+
! TL:     *                       (v(i-1,j,k,nrhs)+v(i  ,j,k,nrhs))*
! TL:     *                       (tl_Hudy(i,j-1,k)+tl_Hudy(i,j  ,k)) )

         tmp=0.25d0*ad_VFx(i,j)
         tmp1=tmp*(Hudy(i,j-1,k,nrhs)+Hudy(i,j  ,k,nrhs))
         tmp2=tmp*(v(i-1,j,k,nrhs)+v(i  ,j,k,nrhs))

         ad_v(i-1,j,k,nrhs)=ad_v(i-1,j,k,nrhs)+tmp1
         ad_v(i  ,j,k,nrhs)=ad_v(i  ,j,k,nrhs)+tmp1
         ad_Hudy(i,j-1,k)=ad_Hudy(i,j-1,k)+tmp2
         ad_Hudy(i,j  ,k)=ad_Hudy(i,j  ,k)+tmp2

         ad_VFx(i,j)=0.d0
        END DO
       END DO

! V-flux of u momentum (at psi locations):

       DO j=JstrU,JendU+1
        DO i=IstrU,IendU
! NL:         UFe(i,j)=0.25d0*(u(i,j-1,k,nrhs)+u(i,j  ,k,nrhs))* 
! NL:     *                   (Hvdx(i-1,j,k)+Hvdx(i  ,j,k))

! TL:         tl_UFe(i,j)=0.25d0*(
! TL:     *              (tl_u(i,j-1,k,nrhs)+tl_u(i,j  ,k,nrhs))*
! TL:     *              (Hvdx(i-1,j,k)+Hvdx(i  ,j,k))+
! TL:     *                       (u(i,j-1,k,nrhs)+u(i,j  ,k,nrhs))*
! TL:     *                       (tl_Hvdx(i-1,j,k)+tl_Hvdx(i  ,j,k)) )

         tmp=0.25d0*ad_UFe(i,j)
         tmp1=tmp*(Hvdx(i-1,j,k,nrhs)+Hvdx(i  ,j,k,nrhs))
         tmp2=tmp*(u(i,j-1,k,nrhs)+u(i,j  ,k,nrhs))

         ad_u(i,j-1,k,nrhs)=ad_u(i,j-1,k,nrhs)+tmp1
         ad_u(i,j  ,k,nrhs)=ad_u(i,j  ,k,nrhs)+tmp1
         ad_Hvdx(i-1,j,k)=ad_Hvdx(i-1,j,k)+tmp2
         ad_Hvdx(i  ,j,k)=ad_Hvdx(i  ,j,k)+tmp2

         ad_UFe(i,j)=0.d0

        END DO
       END DO

! U-flux of u momentum (at rho locations)
       DO j=JstrU,JendU
        DO i=IstrU-1,IendU
! NL:         UFx(i,j)=0.25d0*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))*
! NL:     *                   (Hudy(i,j,k)+Hudy(i+1,j,k))

! TL:         tl_UFx(i,j)=0.25d0*(
! TL:     *              (tl_u(i,j,k,nrhs)+tl_u(i+1,j,k,nrhs))*
! TL:     *              (Hudy(i,j,k)+Hudy(i+1,j,k)) + 
! TL:     *                       (u(i,j,k,nrhs)+u(i+1,j,k,nrhs))*
! TL:     *                       (tl_Hudy(i,j,k)+tl_Hudy(i+1,j,k)) )

         tmp=0.25d0*ad_UFx(i,j)
         tmp1=tmp*(Hudy(i,j,k,nrhs)+Hudy(i+1,j,k,nrhs))
         tmp2=tmp*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))

         ad_u(i,j,k,nrhs)=ad_u(i,j,k,nrhs)+tmp1
         ad_u(i+1,j,k,nrhs)=ad_u(i+1,j,k,nrhs)+tmp1
         ad_Hudy(i,j,k)=ad_Hudy(i,j,k)+tmp2
         ad_Hudy(i+1,j,k)=ad_Hudy(i+1,j,k)+tmp2

         ad_UFx(i,j)=0.d0

        END DO
       END DO

      END DO K_LOOP

      RETURN
      END

#endif /* ADJOINT */
