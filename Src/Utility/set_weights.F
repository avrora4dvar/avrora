#include "cppdefs.h"

      SUBROUTINE set_weights

      USE mod_ocean
#ifdef MPI
      USE mod_parallel
#endif

      implicit none

! NOTE (from ROMS):
! ndtfast: Number of b/t timesteps between each baroclinic timestep.
! nfast:   Number of b/t timesteps needed to compute time-averaged 
!          barotropic variables centered at time level n+1.

      integer :: i, j, iter
      real(8) :: gamma, scale
      real(8) :: cff1, cff2
      real(8) :: wsum, shift, cff   ! in roms: r16 (really needed?)

!=======================================================================
!  Compute time-averaging filter for barotropic fields.
!=======================================================================
!
!  Initialize both sets of weights to zero.
!
      nfast=0
      DO i=1,2*ndtfast
        weight(1,i)=0.0d0
        weight(2,i)=0.0d0
      END DO
!
!-----------------------------------------------------------------------
!  Power-law shape filters.
!-----------------------------------------------------------------------
!
!  The power-law shape filters are given by:
!
!     F(xi)=xi^Falpha*(1-xi^Fbeta)-Fgamma*xi 
!
!  where xi=scale*i/ndtfast; and scale, Falpha, Fbeta, Fgamma, and
!  normalization are chosen to yield the correct zeroth-order
!  (normalization), first-order (consistency), and second-order moments,
!  resulting in overall second-order temporal accuracy for time-averaged
!  barotropic motions resolved by baroclinic time step. There parameters
!  are set in "mod_scalars".
!
      scale=(Falpha+1.0d0)*(Falpha+Fbeta+1.0d0)/                      
     *      ((Falpha+2.0d0)*(Falpha+Fbeta+2.0d0)*REAL(ndtfast,8))
!
!  Find center of gravity of the primary weighting shape function and
!  iteratively adjust "scale" to place the  centroid exactly at
!  "ndtfast".
!
      gamma=Fgamma*MAX(0.d0, 1.d0-10.d0/REAL(ndtfast,8))
      DO iter=1,16
        nfast=0
        DO i=1,2*ndtfast
          cff=scale*REAL(i,8)
          weight(1,i)=cff**Falpha-cff**(Falpha+Fbeta)-gamma*cff
          IF (weight(1,i) .gt. 0.d0) nfast=i
          IF ((nfast .gt. 0) .and. (weight(1,i) .lt. 0.d0)) THEN
            weight(1,i)=0.d0
          END IF
        END DO
        wsum=0.d0  ! r16 in roms
        shift=0.d0 ! r16 in roms  
        DO i=1,nfast
          wsum=wsum+weight(1,i)
          shift=shift+weight(1,i)*REAL(i,8)
        END DO
        scale=scale*shift/(wsum*REAL(ndtfast,8))
      END DO
!
!-----------------------------------------------------------------------
!  Post-processing of primary weights.
!-----------------------------------------------------------------------
!
!  Although it is assumed that the initial settings of the primary
!  weights has its center of gravity "reasonably close" to NDTFAST,
!  it may be not so according to the discrete rules of integration.
!  The following procedure is designed to put the center of gravity
!  exactly to NDTFAST by computing mismatch (NDTFAST-shift) and
!  applying basically an upstream advection of weights to eliminate
!  the mismatch iteratively. Once this procedure is complete primary
!  weights are normalized.
!
!  Find center of gravity of the primary weights and subsequently
!  calculate the mismatch to be compensated.
!
      DO iter=1,ndtfast
        wsum=0.0d0   ! r16 in roms
        shift=0.d0   ! r16 in roms
        DO i=1,nfast 
          wsum=wsum+weight(1,i)
          shift=shift+REAL(i,8)*weight(1,i)
        END DO
        shift=shift/wsum
        cff=REAL(ndtfast,8)-shift
!
!  Apply advection step using either whole, or fractional shifts.
!  Notice that none of the four loops here is reversible.
!
        IF (cff .gt. 1.d0) THEN
          nfast=nfast+1
          DO i=nfast,2,-1
            weight(1,i)=weight(1,i-1)
          END DO
          weight(1,1)=0.d0
        ELSE IF (cff.gt.0.d0) THEN     ! r16 in roms
          wsum=1.d0-cff                ! r16 in roms
          DO i=nfast,2,-1
            weight(1,i)=wsum*weight(1,i)+cff*weight(1,i-1)
          END DO
          weight(1,1)=wsum*weight(1,1)
        ELSE IF (cff.lt.-1.d0) THEN
          nfast=nfast-1
          DO i=1,nfast,+1
            weight(1,i)=weight(1,i+1)
          END DO
          weight(1,nfast+1)=0.d0
        ELSE IF (cff.lt.0.d0) THEN
          wsum=1.d0+cff
          DO i=1,nfast-1,+1
            weight(1,i)=wsum*weight(1,i)-cff*weight(1,i+1)
          END DO
          weight(1,nfast)=wsum*weight(1,nfast)
        END IF
      END DO
!
!  Set SECONDARY weights assuming that backward Euler time step is used
!  for free surface.  Notice that array weight(2,i) is assumed to
!  have all-zero status at entry in this segment of code.
!
      DO j=1,nfast
        cff=weight(1,j)
        DO i=1,j
          weight(2,i)=weight(2,i)+cff
        END DO
      END DO
!
!  Normalize both set of weights.
!
      wsum=0.d0
      cff=0.d0
      DO i=1,nfast
        wsum=wsum+weight(1,i)
        cff=cff+weight(2,i)
      END DO
      wsum=1.d0/wsum
      cff=1.d0/cff
      DO i=1,nfast
        weight(1,i)=wsum*weight(1,i)
        weight(2,i)=cff*weight(2,i)
      END DO

#ifdef MPI
      IF (myrank.eq.0) THEN
#endif

!
!  Report weights.
!
      WRITE (*,10) ndtfast, nfast
      cff=0.d0
      cff1=0.d0
      cff2=0.d0
      wsum=0.d0    ! r16 in roms
      shift=0.d0   ! r16 in roms
      DO i=1,nfast
       cff=cff+weight(1,i)
       cff1=cff1+weight(1,i)*REAL(i,8)
       cff2=cff2+weight(1,i)*REAL(i*i,8)
       wsum=wsum+weight(2,i)
       shift=shift+weight(2,i)*(REAL(i,8)-0.5d0)
       WRITE (*,20) i, weight(1,i), weight(2,i), cff, wsum
      END DO
      cff1=cff1/REAL(ndtfast,8)
      cff2=cff2/(REAL(ndtfast,8)*REAL(ndtfast,8))
      shift=shift/REAL(ndtfast,8)
      WRITE (*,30) ndtfast, nfast,
     *                  REAL(nfast,8)/REAL(ndtfast,8)
      WRITE (*,40) cff1, cff2, shift, cff, wsum, Fgamma, gamma
        IF (cff2.lt.1.0001d0) WRITE (*,50)

#ifdef MPI
      END IF
#endif

!
   10 FORMAT (/,' Time Splitting Weights: ndtfast = ',i3,4x,'nfast = ', 
     *        i3,/,/,4x,'Primary',12x,'Secondary',12x,                  
     *        'Accumulated to Current Step',/)
   20 FORMAT (i3,4f19.16)
   30 FORMAT (/,1x,'ndtfast, nfast = ',2i4,3x,'nfast/ndtfast = ',f7.5)
   40 FORMAT (/,1x,'Centers of gravity and integrals ',                  
     *        '(values must be 1, 1, approx 1/2, 1, 1):',/,              
     *        /,3x,5F15.12,/,/,1x,'Power filter parameters, ',           
     *        'Fgamma, gamma = ', f8.5,2x,f8.5)
   50 FORMAT (/,' WARNING: unstable weights, reduce parameter',          
     *          ' Fgamma in mod_scalars.F',/)
      RETURN
      END SUBROUTINE set_weights
