#include "cppdefs.h"

      MODULE mod_eoscoef

        implicit none

#ifdef NONLIN_EOS

!
! Modified version of Hernan's code:
!

!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Polynomial  expansion  coefficients for the computation of          !
!  "in situ" density  and other associated quantities via the          !
!  nonlinear equation of state for seawater  as a function of          !
!  potential temperature, salinity, and pressure (Jackett and          !
!  McDougall, 1992).                                                   !
!                                                                      !
!=======================================================================
!

        real(8), parameter :: A00 = +1.909256d+04
        real(8), parameter :: A01 = +2.098925d+02
        real(8), parameter :: A02 = -3.041638d+00
        real(8), parameter :: A03 = -1.852732d-03
        real(8), parameter :: A04 = -1.361629d-05
        real(8), parameter :: B00 = +1.044077d+02
        real(8), parameter :: B01 = -6.500517d+00
        real(8), parameter :: B02 = +1.553190d-01
        real(8), parameter :: B03 = +2.326469d-04
        real(8), parameter :: D00 = -5.587545d+00
        real(8), parameter :: D01 = +7.390729d-01
        real(8), parameter :: D02 = -1.909078d-02
        real(8), parameter :: E00 = +4.721788d-01
        real(8), parameter :: E01 = +1.028859d-02
        real(8), parameter :: E02 = -2.512549d-04
        real(8), parameter :: E03 = -5.939910d-07
        real(8), parameter :: F00 = -1.571896d-02
        real(8), parameter :: F01 = -2.598241d-04
        real(8), parameter :: F02 = +7.267926d-06
        real(8), parameter :: G00 = +2.042967d-03
        real(8), parameter :: G01 = +1.045941d-05
        real(8), parameter :: G02 = -5.782165d-10
        real(8), parameter :: G03 = +1.296821d-07
        real(8), parameter :: H00 = -2.595994d-07
        real(8), parameter :: H01 = -1.248266d-09
        real(8), parameter :: H02 = -3.508914d-09
        real(8), parameter :: Q00 = +9.99842594d+02
        real(8), parameter :: Q01 = +6.793952d-02
        real(8), parameter :: Q02 = -9.095290d-03
        real(8), parameter :: Q03 = +1.001685d-04
        real(8), parameter :: Q04 = -1.120083d-06
        real(8), parameter :: Q05 = +6.536332d-09
        real(8), parameter :: U00 = +8.24493d-01
        real(8), parameter :: U01 = -4.08990d-03
        real(8), parameter :: U02 = +7.64380d-05
        real(8), parameter :: U03 = -8.24670d-07
        real(8), parameter :: U04 = +5.38750d-09
        real(8), parameter :: V00 = -5.72466d-03
        real(8), parameter :: V01 = +1.02270d-04
        real(8), parameter :: V02 = -1.65460d-06
        real(8), parameter :: W00 = +4.8314d-04

#else /* not defined NONLIN_EOS */
       
        real(8) :: R0,Tcoef,Scoef,T0,S0      

#endif /* NONLIN_EOS or not */

      END MODULE mod_eoscoef
