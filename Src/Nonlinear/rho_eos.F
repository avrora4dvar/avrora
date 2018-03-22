#include "cppdefs.h"

# if defined NONLIN_EOS

      SUBROUTINE rho_eos
      
!=======================================================================
!  ROMS COMMENT:                                                       !
!  This routine computes  "in situ" density and other associated       !
!  quantitites as a function of potential temperature,  salinity,      !
!  and pressure from a polynomial expression (Jackett & McDougall,     !
!  1992). The polynomial expression was found from fitting to 248      !
!  values  in the  oceanographic  ranges of  salinity,  potential      !
!  temperature,  and pressure.  It  assumes no pressure variation      !
!  along geopotential surfaces, that is, depth (meters; negative)      !
!  and pressure (dbar; assumed negative here) are interchangeable.     !
!                                                                      !
!  Check Values: (T=3 C, S=35.5 PSU, Z=-5000 m)                        !
!                                                                      !
!     alpha = 2.1014611551470d-04 (1/Celsius)                          !
!     beta  = 7.2575037309946d-04 (1/PSU)                              !
!     gamma = 3.9684764511766d-06 (1/Pa)                               !
!     den   = 1050.3639165364     (kg/m3)                              !
!     den1  = 1028.2845117925     (kg/m3)                              !
!     sound = 1548.8815240223     (m/s)                                !
!     bulk  = 23786.056026320     (Pa)                                 !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!  Jackett, D. R. and T. J. McDougall, 1995, Minimal Adjustment of     !
!    Hydrostatic Profiles to Achieve Static Stability, J. of Atmos.    !
!    and Oceanic Techn., vol. 12, pp. 381-389.                         !
!                                                                      !
!=======================================================================
!
      USE mod_ocean
      USE mod_eoscoef 

      implicit none

      integer :: i,j,k      

      real(8) :: Tt, Ts, sqrtTs, Tp, Tpr10,cff
      real(8) :: C(0:9)
      real(8) :: bulk (LBi:UBi,N)
      real(8) :: bulk0(LBi:UBi,N)
      real(8) :: bulk1(LBi:UBi,N)
      real(8) :: bulk2(LBi:UBi,N)
      real(8) :: den(LBi:UBi,N)
      real(8) :: den1(LBi:UBi,N)
 

!
!=======================================================================
!  Nonlinear equation of state.  Notice that this equation of state
!  is only valid for potential temperature range of -2C to 40C and
!  a salinity range of 0 PSU to 42 PSU.
!=======================================================================
!
      DO j=Jstr,Jend
       DO k=1,N
        DO i=Istr,Iend
!
!  Check temperature and salinity lower values. Assign depth to the
!  pressure.
!  
            Tt=MAX(-2.0d0,t(i,j,k,nrhs,itemp))
            Ts=MAX(0.0d0,t(i,j,k,nrhs,isalt))
            sqrtTs=SQRT(Ts)
            Tp=z_r(i,j,k)
            Tpr10=0.1d0*Tp
!
!-----------------------------------------------------------------------
!  Compute density (kg/m3) at standard one atmosphere pressure.
!-----------------------------------------------------------------------
!
            C(0)=Q00+Tt*(Q01+Tt*(Q02+Tt*(Q03+Tt*(Q04+Tt*Q05))))
            C(1)=U00+Tt*(U01+Tt*(U02+Tt*(U03+Tt*U04)))
            C(2)=V00+Tt*(V01+Tt*V02)
!
            den1(i,k)=C(0)+Ts*(C(1)+sqrtTs*C(2)+Ts*W00)
!
!-----------------------------------------------------------------------
!  Compute secant bulk modulus.
!-----------------------------------------------------------------------
!
            C(3)=A00+Tt*(A01+Tt*(A02+Tt*(A03+Tt*A04)))
            C(4)=B00+Tt*(B01+Tt*(B02+Tt*B03))
            C(5)=D00+Tt*(D01+Tt*D02)
            C(6)=E00+Tt*(E01+Tt*(E02+Tt*E03))
            C(7)=F00+Tt*(F01+Tt*F02)
            C(8)=G01+Tt*(G02+Tt*G03)
            C(9)=H00+Tt*(H01+Tt*H02)
!
            bulk0(i,k)=C(3)+Ts*(C(4)+sqrtTs*C(5))
            bulk1(i,k)=C(6)+Ts*(C(7)+sqrtTs*G00)
            bulk2(i,k)=C(8)+Ts*C(9)
            bulk (i,k)=bulk0(i,k)-Tp*(bulk1(i,k)-Tp*bulk2(i,k))
!
!-----------------------------------------------------------------------
!  Compute local "in situ" density anomaly (kg/m3 - 1000).
!-----------------------------------------------------------------------
!
            cff=1.0d0/(bulk(i,k)+Tpr10)
            den(i,k)=den1(i,k)*bulk(i,k)*cff
            den(i,k)=den(i,k)-1000.0d0
            den(i,k)=den(i,k)*REAL(mask_r(i,j),8)
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Load "in situ" density anomaly (kg/m3 - 1000) into the global
!  array
!-----------------------------------------------------------------------
!
       DO k=1,N
        DO i=Istr,Iend
         rho(i,j,k)=den(i,k)
        END DO
       END DO
      END DO
!
!-----------------------------------------------------------------------
!  Exchange boundary data.
!-----------------------------------------------------------------------
!

# ifdef NS_PERIODIC
      CALL exchange_3d(rho,LBi,UBi,LBj,UBj,1,N,Lm,Mm)
# endif

      RETURN
      END

# else   /* not defined NONLIN_EOS */

!
!=======================================================================
!  Linear equation of state.
!=======================================================================

      SUBROUTINE rho_eos

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
      DO j=Jstr,Jend
       DO k=1,N
        DO i=Istr,Iend
         rho(i,j,k)=R0*(1.d0-Tcoef*(t(i,j,k,nrhs,itemp)-T0)
     *                      +Scoef*(t(i,j,k,nrhs,isalt)-S0))
         rho(i,j,k)=rho(i,j,k)-1000.d0
         rho(i,j,k)=rho(i,j,k)*REAL(mask_r(i,j),8)
        END DO
       END DO
      END DO

# ifdef NS_PERIODIC
      CALL exchange_3d(rho,LBi,UBi,LBj,UBj,1,N,Lm,Mm)
# endif

      RETURN
      END

# endif /* defined ... or not ... NONLIN_EOS */

