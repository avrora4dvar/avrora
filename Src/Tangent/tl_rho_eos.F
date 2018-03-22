#include "cppdefs.h"

#ifdef TANGENT
# if defined NONLIN_EOS

      SUBROUTINE tl_rho_eos
      USE mod_ocean
      USE mod_eoscoef 

      implicit none

      integer :: i,j,k      

      WRITE (*,*) 'STOP in tl_tho_eos: TL for NONLIN_EOS not ready yet'
      STOP
      
      RETURN
      END

# else /* undef NONLIN_EOS, USE LINEAR EQUATION OF STATE */

!
!=======================================================================
!  Linear equation of state.
!=======================================================================

      SUBROUTINE tl_rho_eos

      USE mod_ocean
      USE mod_eoscoef 

      implicit none
      integer :: i,j,k

!
!-----------------------------------------------------------------------
!  Compute "in situ" density anomaly (kg/m3 - 1000) using the linear
!  equation of state.
!-----------------------------------------------------------------------
!

#ifdef MPI
      DO j=LBj,UBj
       DO k=1,N
        DO i=LBi,UBi
#else /* serial */
      DO j=Jstr,Jend
       DO k=1,N
        DO i=Istr,Iend
#endif /* MPI or serial */

! NL (keep since fwd rho will be needed):
         rho(i,j,k)=R0*(1.d0-Tcoef*(t(i,j,k,nrhs,itemp)-T0)
     *                      +Scoef*(t(i,j,k,nrhs,isalt)-S0))
         rho(i,j,k)=rho(i,j,k)-1000.d0
         rho(i,j,k)=rho(i,j,k)*REAL(mask_r(i,j),8)

! TL:         
         tl_rho(i,j,k)=R0*(-Tcoef*tl_t(i,j,k,nrhs,itemp)
     *                     +Scoef*tl_t(i,j,k,nrhs,isalt))
         tl_rho(i,j,k)=tl_rho(i,j,k)*REAL(mask_r(i,j),8)
        END DO
       END DO
      END DO

#ifndef MPI
# ifdef NS_PERIODIC
      CALL exchange_3d(rho,LBi,UBi,LBj,UBj,1,N,Lm,Mm)
      CALL exchange_3d(tl_rho,LBi,UBi,LBj,UBj,1,N,Lm,Mm)
# endif
#endif

      RETURN
      END

# endif /* NONLIN_EOS ... else ... */

#endif /* TANGENT */
