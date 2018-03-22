#include "cppdefs.h"

#ifdef ADJOINT
# if defined NONLIN_EOS

      SUBROUTINE ad_rho_eos
      USE mod_ocean
      USE mod_eoscoef 

      implicit none

      integer :: i,j,k      

      WRITE (*,*) 'STOP in ad_tho_eos: ADJ for NONLIN_EOS not ready yet'
      STOP
      
      RETURN
      END

# else /* undef NONLIN_EOS, USE LINEAR EQUATION OF STATE */
!
!=======================================================================
!  Linear equation of state.
!=======================================================================

      SUBROUTINE ad_rho_eos

      USE mod_ocean
      USE mod_eoscoef 

      implicit none
      integer :: i,j,k

# ifdef NS_PERIODIC
      CALL ad_exchange_3d(ad_rho,LBi,UBi,LBj,UBj,1,N,Lm,Mm)
# endif
      
      DO j=Jstr,Jend
       DO k=1,N
        DO i=Istr,Iend
         
! tl line 2:         tl_rho(i,j,k)=tl_rho(i,j,k)*REAL(mask_r(i,j),8)
         ad_rho(i,j,k)=ad_rho(i,j,k)*REAL(mask_r(i,j),8)

! tl line 1:         tl_rho(i,j,k)=R0*(-Tcoef*tl_t(i,j,k,nrhs,itemp)
! tl line 1:     *                     +Scoef*tl_t(i,j,k,nrhs,isalt))

         ad_t(i,j,k,nrhs,itemp)=ad_t(i,j,k,nrhs,itemp)
     *                          -R0*Tcoef*ad_rho(i,j,k)
         ad_t(i,j,k,nrhs,isalt)=ad_t(i,j,k,nrhs,isalt)
     *                          +R0*Scoef*ad_rho(i,j,k)
         ad_rho(i,j,k)=0.d0

        END DO
       END DO
      END DO

      RETURN
      END

# endif /* NONLIN_EOS ... else... */
#endif /* ADJOINT */
