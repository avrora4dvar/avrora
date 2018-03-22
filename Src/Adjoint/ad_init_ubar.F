#include "cppdefs.h"

#ifdef ADJOINT
      SUBROUTINE ad_init_ubar
      USE mod_ocean

      implicit none
      integer :: i,j,k,nini
      real(8) :: Dini,ad_flux,tmp
 
      nini=1
      ad_flux=0.d0

#ifdef NS_PERIODIC

      CALL ad_exchange_2d(ad_vbar(:,:,nini),LBi,UBi,LBj,UBj,Lm,Mm)
      CALL ad_exchange_2d(ad_ubar(:,:,nini),LBi,UBi,LBj,UBj,Lm,Mm)

#endif /* NS_PERIODIC */


      DO j=JstrV,JendV
       DO i=IstrV,IendV

        IF (mask_v(i,j).ne.0) THEN
         Dini=0.5d0*(fwd_zeta(i,j,nini)+fwd_zeta(i,j-1,nini)
     *               +h(i,j)+h(i,j-1))

! ADJ to: Divide flux by the depth:
!
! TL:         tl_vbar(i,j,nini)=tl_flux/Dini
         ad_flux=ad_flux+ad_vbar(i,j,nini)/Dini
         ad_vbar(i,j,nini)=0.d0

! ADJ to:  Flux (v):
!
         DO k=1,N
! TL (recurs):          tl_flux=tl_flux+tl_Hz_v(i,j,k)*v(i,j,k,nini)
!      *                   +Hz_v(i,j,k)*tl_v(i,j,k,nini)
          ad_Hz_v(i,j,k)=ad_Hz_v(i,j,k)+ad_flux*v(i,j,k,nini)
          ad_v(i,j,k,nini)=ad_v(i,j,k,nini)+ad_flux*Hz_v(i,j,k)
          
         END DO

! TL:         tl_flux=0.d0
         ad_flux=0.d0

        END IF

       END DO
      END DO
       

! TL:      tl_ubar(IendU+1,:,nini)=0.d0
      ad_ubar(IendU+1,:,nini)=0.d0

! TL:      tl_ubar(IstrU-1,:,nini)=0.d0
      ad_ubar(IstrU-1,:,nini)=0.d0

      DO j=JstrU,JendU
       DO i=IstrU,IendU
        IF (mask_u(i,j).ne.0) THEN
         Dini=0.5d0*(fwd_zeta(i,j,nini)+fwd_zeta(i-1,j,nini)
     *               +h(i,j)+h(i-1,j))

! TL 4:         tl_ubar(i,j,nini)=tl_flux/Dini
         ad_flux=ad_flux+ad_ubar(i,j,nini)/Dini
         ad_ubar(i,j,nini)=0.d0

         DO k=1,N
! TL 2 (recurs):          tl_flux=tl_flux+tl_Hz_u(i,j,k)*u(i,j,k,nini)
! TL 2:     *                   +Hz_u(i,j,k)*tl_u(i,j,k,nini)
          ad_Hz_u(i,j,k)=ad_Hz_u(i,j,k)+ad_flux*u(i,j,k,nini)
          ad_u(i,j,k,nini)=ad_u(i,j,k,nini)+ad_flux*Hz_u(i,j,k)
         END DO

! Flux (u):
! TL 1:         tl_flux=0.d0
         ad_flux=0.d0

        END IF

       END DO
      END DO

      RETURN
      END

#endif /* ADJOINT */

