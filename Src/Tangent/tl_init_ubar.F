#include "cppdefs.h"

      SUBROUTINE tl_init_ubar
#ifdef TANGENT
      USE mod_ocean
# ifdef MPI
      USE mod_parallel
# endif
      implicit none
      integer :: i,j,k,nini
      real(8) :: Dini,tl_flux

#ifdef MPI
      if (myrank.eq.0) then
       write (*,*) 'TL_INIT_UBAR: note simplified tl_ubar=tl_flux/Dini';
      end if 
#endif

      nini=1
 
      DO j=JstrU,JendU
       DO i=IstrU,IendU
        IF (mask_u(i,j).ne.0) THEN

! Flux (u):
         tl_flux=0.d0
         DO k=1,N
          tl_flux=tl_flux+tl_Hz_u(i,j,k)*u(i,j,k,nini)
     *                   +Hz_u(i,j,k)*tl_u(i,j,k,nini)
         END DO

! Divide flux by the depth:
         Dini=0.5d0*(fwd_zeta(i,j,nini)+fwd_zeta(i-1,j,nini)
     *               +h(i,j)+h(i-1,j))
         tl_ubar(i,j,nini)=tl_flux/Dini

        END IF

       END DO
      END DO

! next two lines commented in the MPI version (3.1) (must be 0 anyway)
!      tl_ubar(IstrU-1,:,nini)=0.d0
!      tl_ubar(IendU+1,:,nini)=0.d0

      DO j=JstrV,JendV
       DO i=IstrV,IendV

        IF (mask_v(i,j).ne.0) THEN

! Flux (v):
         tl_flux=0.d0
         DO k=1,N
          tl_flux=tl_flux+tl_Hz_v(i,j,k)*v(i,j,k,nini)
     *                   +Hz_v(i,j,k)*tl_v(i,j,k,nini)
         END DO

! Divide flux by the depth:
         Dini=0.5d0*(fwd_zeta(i,j,nini)+fwd_zeta(i,j-1,nini)
     *               +h(i,j)+h(i,j-1))
         tl_vbar(i,j,nini)=tl_flux/Dini

        END IF

       END DO
      END DO

#ifdef MPI
      CALL mpi_exchange_2d_2(tl_ubar(LBi:UBi,LBj:UBj,nini),
     *                     LBi,UBi,LBj,UBj,
     *                     IstrU_E,IendU_E,JstrU_E,JendU_E,
     *                     comm2d, stride, 
     *                     rankleft,rankright,ranktop,rankbottom) 

      CALL mpi_exchange_2d_2(tl_vbar(LBi:UBi,LBj:UBj,nini),
     *                     LBi,UBi,LBj,UBj,
     *                     IstrV_E,IendV_E,JstrV_E,JendV_E,
     *                     comm2d, stride, 
     *                     rankleft,rankright,ranktop,rankbottom) 
#endif

#ifdef NS_PERIODIC
# ifdef MPI
      STOP('periodic BC not coded in MPI')
# else /* serial */

      CALL exchange_2d(tl_ubar(:,:,nini),LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(tl_vbar(:,:,nini),LBi,UBi,LBj,UBj,Lm,Mm)
# endif /* MPI or serial */

#endif /* NS_PERIODIC */
            
#endif /* TANGENT */
      RETURN
      END



