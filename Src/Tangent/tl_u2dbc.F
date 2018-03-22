#include "cppdefs.h"

#ifdef TANGENT
!=====================================================================
#ifdef MPI
!=====================================================================
      SUBROUTINE tl_u2dbc(Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj,  
     *               krhs, kstp, knew,
     *               iif,PREDICTOR,g,h,mask_u,
     *               B_West,B_East,B_South,B_North, ! <- new line, comp ot v2.1
     *               tl_ubar,tl_zeta)

      implicit none
      integer, intent(in) :: Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj      
      integer, intent(in) :: krhs, kstp,knew,iif
      logical, intent(in) :: PREDICTOR
      integer, intent(in) :: mask_u(LBi:UBi,LBj:UBj)
      real(8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(8), intent(in) :: tl_zeta(LBi:UBi,LBj:UBj,3)
      real(8), intent(in) :: g

      logical, intent(in) :: B_West,B_East,B_South,B_North

      real(8) :: cff

      integer :: i,j,know

      IF (iif.eq.1) THEN
        know=krhs
      ELSE IF (PREDICTOR) THEN
        know=krhs
      ELSE
        know=kstp
      END IF

!--------------------------------------
! U2: WEST
!--------------------------------------

      IF (B_West) THEN

#if defined WESTERN_WALL
      DO j=Jstr,Jend
       tl_ubar(Istr,j,knew)=0.d0
      END DO      

#else
# ifdef WEST_M2FLATHER
      DO j=Jstr,Jend
       cff=0.5d0*(h(Istr-1,j)+h(Istr,j))
       tl_ubar(Istr,j,knew)=
     *          -DSQRT(g/cff)
     *          *0.5d0*(tl_zeta(Istr-1,j,know)+tl_zeta(Istr,j,know))
     *          *DBLE(mask_u(Istr,j))
      END DO      
      
# else
      M2, WESTERN BOUNDARY, UNDEFINED OPTION

# endif
#endif /* defined or or not defined  WESTERN_WALL */

      END IF ! B_West

!--------------------------------------
! U2: EAST
!--------------------------------------
      IF (B_East) THEN

#if defined EASTERN_WALL
      DO j=Jstr,Jend
       tl_ubar(Iend+1,j,knew)=0.d0
      END DO      

#else
# ifdef EAST_M2FLATHER
      DO j=Jstr,Jend
       cff=0.5d0*(h(Iend,j)+h(Iend+1,j))
       tl_ubar(Iend+1,j,knew)= 
     *    DSQRT(g/cff)
     *   *0.5d0*(tl_zeta(Iend,j,know)+tl_zeta(Iend+1,j,know))
     *   *DBLE(mask_u(Iend+1,j))
      END DO      
# else
      M2, EASTERN BOUNDARY, UNDEFINED OPTION

# endif
#endif /* defined or or not defined  EASTERN_WALL */

      END IF ! B_East

# ifndef NS_PERIODIC
!--------------------------------------
! U2: NORTH (tangential, no gradient)
!--------------------------------------
      IF (B_North) THEN

      DO i=Istr,Iend     ! v3.1: compute at i=Istr, then overwrite in corners
!      DO i=Istr+1,Iend  ! <- v2.1

       tl_ubar(i,Jend+1,knew)=tl_ubar(i,Jend,knew)
     *                        *DBLE(mask_u(i,Jend+1))
      END DO      

      END IF ! B_North

!--------------------------------------
! U2: SOUTH (tangential, no gradient)
!--------------------------------------
      IF (B_South) THEN

      DO i=Istr,Iend     ! v3.1: compute at i=Istr, then overwrite in corners
!      DO i=Istr+1,Iend  ! <- v2.1

       tl_ubar(i,Jstr-1,knew)=tl_ubar(i,Jstr,knew)
     *                        *DBLE(mask_u(i,Jstr-1))
      END DO      

      END IF ! B_South

!--------------------------------------
! Corners:
!--------------------------------------
      IF (B_North .and. B_West) THEN
       tl_ubar(Istr,Jend+1,knew)=0.5d0*
     *                                (tl_ubar(Istr  ,Jend  ,knew)+
     *                                 tl_ubar(Istr+1,Jend+1,knew))
     *                               *DBLE(mask_u(Istr,Jend+1))
      END IF


      IF (B_North .and. B_East) THEN
       tl_ubar(Iend+1,Jend+1,knew)=0.5d0*
     *                                (tl_ubar(Iend  ,Jend+1,knew)+
     *                                 tl_ubar(Iend+1,Jend  ,knew))
     *                               *DBLE(mask_u(Iend+1,Jend+1))
      END IF


      IF (B_South .and. B_West) THEN
       tl_ubar(Istr,Jstr-1,knew)=0.5d0*
     *                                (tl_ubar(Istr  ,Jstr  ,knew)+
     *                                 tl_ubar(Istr+1,Jstr-1,knew))
     *                               *DBLE(mask_u(Istr,Jstr-1))       
      END IF

      IF (B_South .and. B_East) THEN
       tl_ubar(Iend+1,Jstr-1,knew)=0.5d0*
     *                                (tl_ubar(Iend+1,Jstr  ,knew)+
     *                                 tl_ubar(Iend  ,Jstr-1,knew))
     *                               *DBLE(mask_u(Iend+1,Jstr-1))       
      END IF

# endif /* not defined NS_PERIODIC */

      RETURN
      END

!****************************************************************
! 1/24/11, AK: Note this is the version for MPI

      SUBROUTINE tl_v2dbc(Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj,  
     *               krhs, kstp, knew,
     *               iif,PREDICTOR,g,h,mask_v,
     *               B_West,B_East,B_South,B_North, ! <- new line, comp ot v2.1
     *               tl_vbar,tl_zeta)

      implicit none
      integer, intent(in) :: Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj      
      integer, intent(in) :: krhs, kstp, knew, iif   
      logical, intent(in) :: PREDICTOR
      integer, intent(in) :: mask_v(LBi:UBi,LBj:UBj)
      real(8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
      real(8), intent(in) :: tl_zeta(LBi:UBi,LBj:UBj,3)
      real(8), intent(in) :: g
      logical, intent(in) :: B_West,B_East,B_South,B_North

      real(8) :: cff

      integer :: i,j,know

      IF (iif.eq.1) THEN
        know=krhs
      ELSE IF (PREDICTOR) THEN
        know=krhs
      ELSE
        know=kstp
      END IF

!--------------------------------------
! V2: WEST (tangential, no gradient)
!--------------------------------------
      IF (B_West) THEN

      DO j=Jstr,Jend    ! v3.1: compute at j=Jstr, then overwrite at corners
!      DO j=Jstr+1,Jend ! <- v2.1
       tl_vbar(Istr-1,j,knew)=tl_vbar(Istr,j,knew)
     *                        *DBLE(mask_v(Istr-1,j))
      END DO      
      
      END IF ! B_West

!--------------------------------------
! V2: EAST (tangential, no gradient)
!--------------------------------------
      IF (B_East) THEN

      DO j=Jstr,Jend    ! v3.1: compute at j=Jstr, then overwrite at corners
!      DO j=Jstr+1,Jend ! <- v2.1
       tl_vbar(Iend+1,j,knew)=tl_vbar(Iend,j,knew)
     *                        *DBLE(mask_v(Iend+1,j))
      END DO      

      END IF ! B_East

#ifndef NS_PERIODIC

!--------------------------------------
! V2: NORTH
!--------------------------------------
      IF (B_North) THEN

# if defined NORTHERN_WALL
      DO i=Istr,Iend
       tl_vbar(i,Jend+1,knew)=0.d0
      END DO

# else
#  if defined NORTH_M2FLATHER
      DO i=Istr,Iend
       cff=0.5d0*(h(i,Jend)+h(i,Jend+1))
       tl_vbar(i,Jend+1,knew)=
     *    DSQRT(g/cff)
     *   *0.5d0*(tl_zeta(i,Jend,know)+tl_zeta(i,Jend+1,know))
     *   *DBLE(mask_v(i,Jend+1))
      END DO

#  else
      M2, NORTHERN BOUNDARY, UNDEFINED OPTION
#  endif /* NORTH_M2FLATHER */
# endif /* NORTHERN_WALL */

      END IF ! B_North

!--------------------------------------
! V2: SOUTH
!--------------------------------------
      IF (B_South) THEN

# if defined SOUTHERN_WALL
      DO i=Istr,Iend
       tl_vbar(i,Jstr,knew)=0.d0
      END DO

# else
#  if defined SOUTH_M2FLATHER
      DO i=Istr,Iend
       cff=0.5d0*(h(i,Jstr)+h(i,Jstr-1))
       tl_vbar(i,Jstr,knew)=
     *    -DSQRT(g/cff)
     *   *0.5d0*(tl_zeta(i,Jstr,know)+tl_zeta(i,Jstr-1,know))
     *   *DBLE(mask_v(i,Jstr))
      END DO
#  else
      M2, SOUTHERN BOUNDARY, UNDEFINED OPTION
#  endif /* SOUTH_M2FLATHER */
# endif /* SOUTHERN_WALL */

      END IF ! B_South

!------------
! Corners
!------------
      IF (B_North .and. B_West) THEN
       tl_vbar(Istr-1,Jend+1,knew)=
     *   0.5d0
     *  *(tl_vbar(Istr-1,Jend,knew)+tl_vbar(Istr,Jend+1,knew))
     *  *DBLE(mask_v(Istr-1,Jend+1))
      END IF

      IF (B_North .and. B_East) THEN
       tl_vbar(Iend+1,Jend+1,knew)=
     *   0.5d0
     *  *(tl_vbar(Iend,Jend+1,knew)+tl_vbar(Iend+1,Jend,knew))
     *  *DBLE(mask_v(Iend+1,Jend+1))
      END IF

      IF (B_South .and. B_West) THEN
       tl_vbar(Istr-1,Jstr,knew)=
     *   0.5d0
     *  *(tl_vbar(Istr-1,Jstr+1,knew)+tl_vbar(Istr,Jstr,knew))
     *  *DBLE(mask_v(Istr-1,Jstr))
      END IF

      IF (B_South .and. B_East) THEN
       tl_vbar(Iend+1,Jstr,knew)=
     *   0.5d0
     *  *(tl_vbar(Iend+1,Jstr+1,knew)+tl_vbar(Iend,Jstr,knew))
     *  *DBLE(mask_v(Iend+1,Jstr))
      END IF

#endif /* ifndef NS_PERIODIC */ 

      RETURN
      END

!=====================================================================
#else /* SERIAL, v2.1 untouched */
!=====================================================================
      SUBROUTINE tl_u2dbc(Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj,  
     *               krhs, kstp, knew,
     *               iif,PREDICTOR,g,h,mask_u,
     *               tl_ubar,tl_zeta)

      implicit none
      integer, intent(in) :: Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj      
      integer, intent(in) :: krhs, kstp,knew,iif
      logical, intent(in) :: PREDICTOR
      integer, intent(in) :: mask_u(LBi:UBi,LBj:UBj)
      real(8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(8), intent(in) :: tl_zeta(LBi:UBi,LBj:UBj,3)
      real(8), intent(in) :: g

      real(8) :: cff

      integer :: i,j,know

      IF (iif.eq.1) THEN
        know=krhs
      ELSE IF (PREDICTOR) THEN
        know=krhs
      ELSE
        know=kstp
      END IF

!--------------------------------------
! U2: WEST
!--------------------------------------
#if defined WESTERN_WALL
      DO j=Jstr,Jend
       tl_ubar(Istr,j,knew)=0.d0
      END DO      

#else
# ifdef WEST_M2FLATHER
      DO j=Jstr,Jend
       cff=0.5d0*(h(Istr-1,j)+h(Istr,j))
       tl_ubar(Istr,j,knew)=
     *          -DSQRT(g/cff)
     *          *0.5d0*(tl_zeta(Istr-1,j,know)+tl_zeta(Istr,j,know))
     *          *DBLE(mask_u(Istr,j))
      END DO      
      
# else
      M2, WESTERN BOUNDARY, UNDEFINED OPTION

# endif
#endif /* defined or or not defined  WESTERN_WALL */

!--------------------------------------
! U2: EAST
!--------------------------------------
#if defined EASTERN_WALL
      DO j=Jstr,Jend
       tl_ubar(Iend+1,j,knew)=0.d0
      END DO      

#else
# ifdef EAST_M2FLATHER
      DO j=Jstr,Jend
       cff=0.5d0*(h(Iend,j)+h(Iend+1,j))
       tl_ubar(Iend+1,j,knew)= 
     *    DSQRT(g/cff)
     *   *0.5d0*(tl_zeta(Iend,j,know)+tl_zeta(Iend+1,j,know))
     *   *DBLE(mask_u(Iend+1,j))
      END DO      
# else
      M2, EASTERN BOUNDARY, UNDEFINED OPTION

# endif
#endif /* defined or or not defined  EASTERN_WALL */

# ifndef NS_PERIODIC
!--------------------------------------
! U2: NORTH (tangential, no gradient)
!--------------------------------------
      DO i=Istr+1,Iend
       tl_ubar(i,Jend+1,knew)=tl_ubar(i,Jend,knew)
     *                        *DBLE(mask_u(i,Jend+1))
      END DO      

!--------------------------------------
! U2: SOUTH (tangential, no gradient)
!--------------------------------------
      DO i=Istr+1,Iend
       tl_ubar(i,Jstr-1,knew)=tl_ubar(i,Jstr,knew)
     *                        *DBLE(mask_u(i,Jstr-1))
      END DO      

!--------------------------------------
! Corners:
!--------------------------------------
      tl_ubar(Istr,Jend+1,knew)=0.5d0*
     *                                (tl_ubar(Istr  ,Jend  ,knew)+
     *                                 tl_ubar(Istr+1,Jend+1,knew))
     *                               *DBLE(mask_u(Istr,Jend+1))
      tl_ubar(Iend+1,Jend+1,knew)=0.5d0*
     *                                (tl_ubar(Iend  ,Jend+1,knew)+
     *                                 tl_ubar(Iend+1,Jend  ,knew))
     *                               *DBLE(mask_u(Iend+1,Jend+1))
      tl_ubar(Istr,Jstr-1,knew)=0.5d0*
     *                                (tl_ubar(Istr  ,Jstr  ,knew)+
     *                                 tl_ubar(Istr+1,Jstr-1,knew))
     *                               *DBLE(mask_u(Istr,Jstr-1))       
      tl_ubar(Iend+1,Jstr-1,knew)=0.5d0*
     *                                (tl_ubar(Iend+1,Jstr  ,knew)+
     *                                 tl_ubar(Iend  ,Jstr-1,knew))
     *                               *DBLE(mask_u(Iend+1,Jstr-1))       

# endif /* not defined NS_PERIODIC */

      RETURN
      END

!****************************************************************
      SUBROUTINE tl_v2dbc(Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj,  
     *               krhs, kstp, knew,
     *               iif,PREDICTOR,g,h,mask_v,
     *               tl_vbar,tl_zeta)

      implicit none
      integer, intent(in) :: Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj      
      integer, intent(in) :: krhs, kstp, knew, iif   
      logical, intent(in) :: PREDICTOR
      integer, intent(in) :: mask_v(LBi:UBi,LBj:UBj)
      real(8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
      real(8), intent(in) :: tl_zeta(LBi:UBi,LBj:UBj,3)
      real(8), intent(in) :: g

      real(8) :: cff

      integer :: i,j,know

      IF (iif.eq.1) THEN
        know=krhs
      ELSE IF (PREDICTOR) THEN
        know=krhs
      ELSE
        know=kstp
      END IF

!--------------------------------------
! V2: WEST (tangential, no gradient)
!--------------------------------------
      DO j=Jstr+1,Jend
       tl_vbar(Istr-1,j,knew)=tl_vbar(Istr,j,knew)
     *                        *DBLE(mask_v(Istr-1,j))
      END DO      
      
!--------------------------------------
! V2: EAST (tangential, no gradient)
!--------------------------------------
      DO j=Jstr+1,Jend
       tl_vbar(Iend+1,j,knew)=tl_vbar(Iend,j,knew)
     *                        *DBLE(mask_v(Iend+1,j))
      END DO      


#ifndef NS_PERIODIC

!--------------------------------------
! V2: NORTH
!--------------------------------------
# if defined NORTHERN_WALL
      DO i=Istr,Iend
       tl_vbar(i,Jend+1,knew)=0.d0
      END DO

# else
#  if defined NORTH_M2FLATHER
      DO i=Istr,Iend
       cff=0.5d0*(h(i,Jend)+h(i,Jend+1))
       tl_vbar(i,Jend+1,knew)=
     *    DSQRT(g/cff)
     *   *0.5d0*(tl_zeta(i,Jend,know)+tl_zeta(i,Jend+1,know))
     *   *DBLE(mask_v(i,Jend+1))
      END DO

#  else
      M2, NORTHERN BOUNDARY, UNDEFINED OPTION
#  endif /* NORTH_M2FLATHER */
# endif /* NORTHERN_WALL */

!--------------------------------------
! V2: SOUTH
!--------------------------------------
# if defined SOUTHERN_WALL
      DO i=Istr,Iend
       tl_vbar(i,Jstr,knew)=0.d0
      END DO

# else
#  if defined SOUTH_M2FLATHER
      DO i=Istr,Iend
       cff=0.5d0*(h(i,Jstr)+h(i,Jstr-1))
       tl_vbar(i,Jstr,knew)=
     *    -DSQRT(g/cff)
     *   *0.5d0*(tl_zeta(i,Jstr,know)+tl_zeta(i,Jstr-1,know))
     *   *DBLE(mask_v(i,Jstr))
      END DO
#  else
      M2, SOUTHERN BOUNDARY, UNDEFINED OPTION
#  endif /* SOUTH_M2FLATHER */
# endif /* SOUTHERN_WALL */

!------------
! Corners
!------------
      tl_vbar(Istr-1,Jend+1,knew)=
     *  0.5d0
     * *(tl_vbar(Istr-1,Jend,knew)+tl_vbar(Istr,Jend+1,knew))
     * *DBLE(mask_v(Istr-1,Jend+1))

      tl_vbar(Iend+1,Jend+1,knew)=
     *  0.5d0
     * *(tl_vbar(Iend,Jend+1,knew)+tl_vbar(Iend+1,Jend,knew))
     * *DBLE(mask_v(Iend+1,Jend+1))
      
      tl_vbar(Istr-1,Jstr,knew)=
     *  0.5d0
     * *(tl_vbar(Istr-1,Jstr+1,knew)+tl_vbar(Istr,Jstr,knew))
     * *DBLE(mask_v(Istr-1,Jstr))

      tl_vbar(Iend+1,Jstr,knew)=
     *  0.5d0
     * *(tl_vbar(Iend+1,Jstr+1,knew)+tl_vbar(Iend,Jstr,knew))
     * *DBLE(mask_v(Iend+1,Jstr))


#endif /* ifndef NS_PERIODIC */ 

      RETURN
      END
!=====================================================================
#endif /* MPI or SERIAL */
!=====================================================================

#endif /* TANGENT */
