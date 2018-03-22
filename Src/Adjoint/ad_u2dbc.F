#include "cppdefs.h"

#ifdef ADJOINT
!=====================================================================
#ifdef MPI
!=====================================================================
      SUBROUTINE ad_u2dbc(Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj,  
     *               krhs, kstp, knew,
     *               iif,PREDICTOR,g,h,mask_u,
     *               B_West,B_East,B_South,B_North, ! <- new line, comp ot v2.1
     *               ad_ubar,ad_zeta)

      implicit none

      integer, intent(in) :: Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj      
      integer, intent(in) :: krhs, kstp,knew,iif
      logical, intent(in) :: PREDICTOR
      integer, intent(in) :: mask_u(LBi:UBi,LBj:UBj)
      real(8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
      real(8), intent(in) :: g

      logical, intent(in) :: B_West,B_East,B_South,B_North

      real(8) :: cff,tmp

      integer :: i,j,know

      IF (iif.eq.1) THEN
        know=krhs
      ELSE IF (PREDICTOR) THEN
        know=krhs
      ELSE
        know=kstp
      END IF

# ifndef NS_PERIODIC
!--------------------------------------
! ADJ to Corners:
!--------------------------------------
      IF (B_South .and. B_East) THEN
! tl line 4:      tl_ubar(Iend+1,Jstr-1,knew)=0.5d0*
!     *                                (tl_ubar(Iend+1,Jstr  ,knew)+
!     *                                 tl_ubar(Iend  ,Jstr-1,knew))
!     *                               *DBLE(mask_u(Iend+1,Jstr-1))       

       tmp=0.5d0*DBLE(mask_u(Iend+1,Jstr-1))*ad_ubar(Iend+1,Jstr-1,knew)
       ad_ubar(Iend+1,Jstr  ,knew)=ad_ubar(Iend+1,Jstr  ,knew)+tmp
       ad_ubar(Iend  ,Jstr-1,knew)=ad_ubar(Iend  ,Jstr-1,knew)+tmp
       ad_ubar(Iend+1,Jstr-1,knew)=0.d0
      END IF
      
      IF (B_South .and. B_West) THEN
! tl line 3:      tl_ubar(Istr,Jstr-1,knew)=0.5d0*
!     *                                (tl_ubar(Istr  ,Jstr  ,knew)+
!     *                                 tl_ubar(Istr+1,Jstr-1,knew))
!     *                               *DBLE(mask_u(Istr,Jstr-1))       
      
       tmp=0.5d0*DBLE(mask_u(Istr,Jstr-1))*ad_ubar(Istr,Jstr-1,knew)
       ad_ubar(Istr  ,Jstr  ,knew)=ad_ubar(Istr  ,Jstr  ,knew)+tmp
       ad_ubar(Istr+1,Jstr-1,knew)=ad_ubar(Istr+1,Jstr-1,knew)+tmp
       ad_ubar(Istr,Jstr-1,knew)=0.d0
      END IF


      IF (B_North .and. B_East) THEN
! tl line 2:      tl_ubar(Iend+1,Jend+1,knew)=0.5d0*
!     *                                (tl_ubar(Iend  ,Jend+1,knew)+
!     *                                 tl_ubar(Iend+1,Jend  ,knew))
!     *                               *DBLE(mask_u(Iend+1,Jend+1))

       tmp=0.5d0*DBLE(mask_u(Iend+1,Jend+1))*ad_ubar(Iend+1,Jend+1,knew)
       ad_ubar(Iend  ,Jend+1,knew)=ad_ubar(Iend  ,Jend+1,knew)+tmp
       ad_ubar(Iend+1,Jend  ,knew)=ad_ubar(Iend+1,Jend  ,knew)+tmp
       ad_ubar(Iend+1,Jend+1,knew)=0.d0
      END IF

      IF (B_North .and. B_West) THEN
! tl line 1:      tl_ubar(Istr,Jend+1,knew)=0.5d0*
!     *                                (tl_ubar(Istr  ,Jend  ,knew)+
!     *                                 tl_ubar(Istr+1,Jend+1,knew))
!     *                               *DBLE(mask_u(Istr,Jend+1))

       tmp=0.5d0*DBLE(mask_u(Istr,Jend+1))*ad_ubar(Istr,Jend+1,knew)
       ad_ubar(Istr  ,Jend  ,knew)=ad_ubar(Istr  ,Jend  ,knew)+tmp
       ad_ubar(Istr+1,Jend+1,knew)=ad_ubar(Istr+1,Jend+1,knew)+tmp
       ad_ubar(Istr,Jend+1,knew)=0.d0
      END IF

!--------------------------------------
! ADJ to U2: SOUTH (tangential, no gradient)
!--------------------------------------
      IF (B_South) THEN

       DO i=Istr,Iend     ! v3.1: compute at i=Istr
!       DO i=Istr+1,Iend ! <- v2.1

! TL:       tl_ubar(i,Jstr-1,knew)=tl_ubar(i,Jstr,knew)
!     *                        *DBLE(mask_u(i,Jstr-1))
        ad_ubar(i,Jstr,knew)=ad_ubar(i,Jstr,knew)+
     *    ad_ubar(i,Jstr-1,knew)*DBLE(mask_u(i,Jstr-1))
        ad_ubar(i,Jstr-1,knew)=0.d0

       END DO      
      END IF ! B_South

!--------------------------------------
! ADJ to U2: NORTH (tangential, no gradient)
!--------------------------------------
      IF (B_North) THEN
       DO i=Istr,Iend     ! v3.1: compute at i=Istr, then overwrite in corners
!       DO i=Istr+1,Iend ! v2.1 (see tl for additional comments)

! TL:       tl_ubar(i,Jend+1,knew)=tl_ubar(i,Jend,knew)
!     *                        *DBLE(mask_u(i,Jend+1))
        ad_ubar(i,Jend,knew)=ad_ubar(i,Jend,knew)+
     *    ad_ubar(i,Jend+1,knew)*DBLE(mask_u(i,Jend+1))
        ad_ubar(i,Jend+1,knew)=0.d0

       END DO      
      END IF ! B_North

# endif /* not defined NS_PERIODIC */

!--------------------------------------
! ADJ to U2: EAST
!--------------------------------------
      IF (B_East) THEN

#if defined EASTERN_WALL
      DO j=Jstr,Jend
! TL:       tl_ubar(Iend+1,j,knew)=0.d0
       ad_ubar(Iend+1,j,knew)=0.d0
      END DO      

#else
# ifdef EAST_M2FLATHER
      DO j=Jstr,Jend
       cff=0.5d0*(h(Iend,j)+h(Iend+1,j))

! TL:       tl_ubar(Iend+1,j,knew)= 
!     *    DSQRT(g/cff)
!     *   *0.5d0*(tl_zeta(Iend,j,know)+tl_zeta(Iend+1,j,know))
!     *   *DBLE(mask_u(Iend+1,j))

       tmp=DSQRT(g/cff)*0.5d0*DBLE(mask_u(Iend+1,j))
     *                       *ad_ubar(Iend+1,j,knew)

       ad_zeta(Iend,j,know)=ad_zeta(Iend,j,know)    +tmp
       ad_zeta(Iend+1,j,know)=ad_zeta(Iend+1,j,know)+tmp
       ad_ubar(Iend+1,j,knew)=0.d0

      END DO      

# else
      M2, EASTERN BOUNDARY, UNDEFINED OPTION

# endif
#endif /* defined or or not defined  EASTERN_WALL */
      END IF ! B_East

!--------------------------------------
! ADJ to U2: WEST
!--------------------------------------
      IF (B_West) THEN

#if defined WESTERN_WALL
      DO j=Jstr,Jend
! TL:       tl_ubar(Istr,j,knew)=0.d0
       ad_ubar(Istr,j,knew)=0.d0
      END DO      

#else
# ifdef WEST_M2FLATHER
      DO j=Jstr,Jend
       cff=0.5d0*(h(Istr-1,j)+h(Istr,j))

! TL       tl_ubar(Istr,j,knew)=
!     *          -DSQRT(g/cff)
!     *          *0.5d0*(tl_zeta(Istr-1,j,know)+tl_zeta(Istr,j,know))
!     *          *DBLE(mask_u(Istr,j))
  
       tmp=-DSQRT(g/cff)*0.5d0*DBLE(mask_u(Istr,j))*ad_ubar(Istr,j,knew)

       ad_zeta(Istr-1,j,know)=ad_zeta(Istr-1,j,know)+tmp
       ad_zeta(Istr  ,j,know)=ad_zeta(Istr  ,j,know)+tmp
       ad_ubar(Istr,j,knew)=0.d0

      END DO      
      
# else
      M2, WESTERN BOUNDARY, UNDEFINED OPTION

# endif
#endif /* defined or or not defined  WESTERN_WALL */

      END IF ! B_West


      RETURN
      END

!****************************************************************
      SUBROUTINE ad_v2dbc(Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj,  
     *               krhs, kstp, knew,
     *               iif,PREDICTOR,g,h,mask_v,
     *               B_West,B_East,B_South,B_North, ! <- new line, comp ot v2.1
     *               ad_vbar,ad_zeta)

      implicit none
      integer, intent(in) :: Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj      
      integer, intent(in) :: krhs, kstp, knew,iif   
      logical, intent(in) :: PREDICTOR
      integer, intent(in) :: mask_v(LBi:UBi,LBj:UBj)
      real(8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
      real(8), intent(in) :: g
      logical, intent(in) :: B_West,B_East,B_South,B_North

      real(8) :: cff,tmp

      integer :: i,j,know

      IF (iif.eq.1) THEN
        know=krhs
      ELSE IF (PREDICTOR) THEN
        know=krhs
      ELSE
        know=kstp
      END IF

#ifndef NS_PERIODIC

!------------
! ADJ to Corners
!------------

      IF (B_South .and. B_East) THEN
! tl line 4:      tl_vbar(Iend+1,Jstr,knew)=
!     *  0.5d0
!     * *(tl_vbar(Iend+1,Jstr+1,knew)+tl_vbar(Iend,Jstr,knew))
!     * *DBLE(mask_v(Iend+1,Jstr))

       tmp=0.5d0*DBLE(mask_v(Iend+1,Jstr))*ad_vbar(Iend+1,Jstr,knew)
       ad_vbar(Iend+1,Jstr+1,knew)=ad_vbar(Iend+1,Jstr+1,knew)+tmp
       ad_vbar(Iend  ,Jstr  ,knew)=ad_vbar(Iend  ,Jstr  ,knew)+tmp
       ad_vbar(Iend+1,Jstr,knew)=0.d0
      END IF
 
      IF (B_South .and. B_West) THEN
! tl line 3:      tl_vbar(Istr-1,Jstr,knew)=
!     *  0.5d0
!     * *(tl_vbar(Istr-1,Jstr+1,knew)+tl_vbar(Istr,Jstr,knew))
!     * *DBLE(mask_v(Istr-1,Jstr))

       tmp=0.5d0*DBLE(mask_v(Istr-1,Jstr))*ad_vbar(Istr-1,Jstr,knew)
       ad_vbar(Istr-1,Jstr+1,knew)=ad_vbar(Istr-1,Jstr+1,knew)+tmp
       ad_vbar(Istr  ,Jstr  ,knew)=ad_vbar(Istr  ,Jstr  ,knew)+tmp
       ad_vbar(Istr-1,Jstr,knew)=0.d0
      END IF
 
      IF (B_North .and. B_East) THEN
! tl line 2:      tl_vbar(Iend+1,Jend+1,knew)=
!     *  0.5d0
!     * *(tl_vbar(Iend,Jend+1,knew)+tl_vbar(Iend+1,Jend,knew))
!     * *DBLE(mask_v(Iend+1,Jend+1))
      
       tmp=0.5d0*DBLE(mask_v(Iend+1,Jend+1))*ad_vbar(Iend+1,Jend+1,knew)
       ad_vbar(Iend,Jend+1,knew)=ad_vbar(Iend,Jend+1,knew)+tmp
       ad_vbar(Iend+1,Jend,knew)=ad_vbar(Iend+1,Jend,knew)+tmp
       ad_vbar(Iend+1,Jend+1,knew)=0.d0
      END IF
 
      IF (B_North .and. B_West) THEN
! tl line 1:      tl_vbar(Istr-1,Jend+1,knew)=
!     *  0.5d0
!     * *(tl_vbar(Istr-1,Jend,knew)+tl_vbar(Istr,Jend+1,knew))
!     * *DBLE(mask_v(Istr-1,Jend+1))

       tmp=0.5d0*DBLE(mask_v(Istr-1,Jend+1))*ad_vbar(Istr-1,Jend+1,knew)
       ad_vbar(Istr-1,Jend,knew)=ad_vbar(Istr-1,Jend,knew)+tmp
       ad_vbar(Istr,Jend+1,knew)=ad_vbar(Istr,Jend+1,knew)+tmp
       ad_vbar(Istr-1,Jend+1,knew)=0.d0
      END IF
 
!--------------------------------------
! ADJ to V2: SOUTH
!--------------------------------------
      IF (B_South) THEN

# if defined SOUTHERN_WALL
      DO i=Istr,Iend
! TL:       tl_vbar(i,Jstr,knew)=0.d0
       ad_vbar(i,Jstr,knew)=0.d0
      END DO

# else
#  if defined SOUTH_M2FLATHER
      DO i=Istr,Iend
       cff=0.5d0*(h(i,Jstr)+h(i,Jstr-1))

! TL:       tl_vbar(i,Jstr,knew)=
!     *    -DSQRT(g/cff)
!     *   *0.5d0*(tl_zeta(i,Jstr,know)+tl_zeta(i,Jstr-1,know))
!     *   *DBLE(mask_v(i,Jstr))

       tmp=-DSQRT(g/cff)*0.5d0*DBLE(mask_v(i,Jstr))*ad_vbar(i,Jstr,knew)
       ad_zeta(i,Jstr,know)=ad_zeta(i,Jstr,know)+tmp
       ad_zeta(i,Jstr-1,know)=ad_zeta(i,Jstr-1,know)+tmp
       ad_vbar(i,Jstr,knew)=0.d0

      END DO
#  else
      M2, SOUTHERN BOUNDARY, UNDEFINED OPTION
#  endif /* SOUTH_M2FLATHER */
# endif /* SOUTHERN_WALL */

      END IF ! B_South

!--------------------------------------
! ADJ to V2: NORTH
!--------------------------------------
      IF (B_North) THEN

# if defined NORTHERN_WALL
      DO i=Istr,Iend
! TL:       tl_vbar(i,Jend+1,knew)=0.d0
       ad_vbar(i,Jend+1,knew)=0.d0
      END DO

# else
#  if defined NORTH_M2FLATHER
      DO i=Istr,Iend

       cff=0.5d0*(h(i,Jend)+h(i,Jend+1))

! TL:       tl_vbar(i,Jend+1,knew)=
!     *    DSQRT(g/cff)
!     *   *0.5d0*(tl_zeta(i,Jend,know)+tl_zeta(i,Jend+1,know))
!     *   *DBLE(mask_v(i,Jend+1))

       tmp=DSQRT(g/cff)*0.5d0*DBLE(mask_v(i,Jend+1))
     *                       *ad_vbar(i,Jend+1,knew)
       ad_zeta(i,Jend,know)  =ad_zeta(i,Jend,know)  +tmp
       ad_zeta(i,Jend+1,know)=ad_zeta(i,Jend+1,know)+tmp
       ad_vbar(i,Jend+1,knew)=0.d0

      END DO

#  else
      M2, NORTHERN BOUNDARY, UNDEFINED OPTION
#  endif /* NORTH_M2FLATHER */
# endif /* NORTHERN_WALL */


#endif /* ifndef NS_PERIODIC */ 

      END IF ! B_North


!--------------------------------------
! ADJ to V2: EAST (tangential, no gradient)
!--------------------------------------
      IF (B_East) THEN

       DO j=Jstr,Jend    ! v3.1: compute at j=Jstr, then overwrite at corners
!      DO j=Jstr+1,Jend
! TL:       tl_vbar(Iend+1,j,knew)=tl_vbar(Iend,j,knew)
!     *                        *DBLE(mask_v(Iend+1,j))
       
        ad_vbar(Iend,j,knew)=ad_vbar(Iend,j,knew)+
     *   DBLE(mask_v(Iend+1,j))*ad_vbar(Iend+1,j,knew)
        ad_vbar(Iend+1,j,knew)=0.d0

       END DO      

      END IF ! B_East


!--------------------------------------
! ADJ to V2: WEST (tangential, no gradient)
!--------------------------------------
      IF (B_West) THEN

       DO j=Jstr,Jend    ! v3.1: compute at j=Jstr
!       DO j=Jstr+1,Jend ! v2.1
! TL:       tl_vbar(Istr-1,j,knew)=tl_vbar(Istr,j,knew)
!     *                        *DBLE(mask_v(Istr-1,j))
       
        ad_vbar(Istr,j,knew)=ad_vbar(Istr,j,knew)+
     *   DBLE(mask_v(Istr-1,j))*ad_vbar(Istr-1,j,knew)
        ad_vbar(Istr-1,j,knew)=0.d0

       END DO      

      END IF ! B_East

      RETURN
      END

!=====================================================================
#else /* SERIAL, v2.1 untouched */
!=====================================================================
      SUBROUTINE ad_u2dbc(Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj,  
     *               krhs, kstp, knew,
     *               iif,PREDICTOR,g,h,mask_u,
     *               ad_ubar,ad_zeta)

      implicit none

      integer, intent(in) :: Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj      
      integer, intent(in) :: krhs, kstp,knew,iif
      logical, intent(in) :: PREDICTOR
      integer, intent(in) :: mask_u(LBi:UBi,LBj:UBj)
      real(8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
      real(8), intent(in) :: g

      real(8) :: cff,tmp

      integer :: i,j,know

      IF (iif.eq.1) THEN
        know=krhs
      ELSE IF (PREDICTOR) THEN
        know=krhs
      ELSE
        know=kstp
      END IF

# ifndef NS_PERIODIC
!--------------------------------------
! ADJ to Corners:
!--------------------------------------
! tl line 4:      tl_ubar(Iend+1,Jstr-1,knew)=0.5d0*
!     *                                (tl_ubar(Iend+1,Jstr  ,knew)+
!     *                                 tl_ubar(Iend  ,Jstr-1,knew))
!     *                               *DBLE(mask_u(Iend+1,Jstr-1))       

      tmp=0.5d0*DBLE(mask_u(Iend+1,Jstr-1))*ad_ubar(Iend+1,Jstr-1,knew)
      ad_ubar(Iend+1,Jstr  ,knew)=ad_ubar(Iend+1,Jstr  ,knew)+tmp
      ad_ubar(Iend  ,Jstr-1,knew)=ad_ubar(Iend  ,Jstr-1,knew)+tmp
      ad_ubar(Iend+1,Jstr-1,knew)=0.d0
      
! tl line 3:      tl_ubar(Istr,Jstr-1,knew)=0.5d0*
!     *                                (tl_ubar(Istr  ,Jstr  ,knew)+
!     *                                 tl_ubar(Istr+1,Jstr-1,knew))
!     *                               *DBLE(mask_u(Istr,Jstr-1))       
      
      tmp=0.5d0*DBLE(mask_u(Istr,Jstr-1))*ad_ubar(Istr,Jstr-1,knew)
      ad_ubar(Istr  ,Jstr  ,knew)=ad_ubar(Istr  ,Jstr  ,knew)+tmp
      ad_ubar(Istr+1,Jstr-1,knew)=ad_ubar(Istr+1,Jstr-1,knew)+tmp
      ad_ubar(Istr,Jstr-1,knew)=0.d0


! tl line 2:      tl_ubar(Iend+1,Jend+1,knew)=0.5d0*
!     *                                (tl_ubar(Iend  ,Jend+1,knew)+
!     *                                 tl_ubar(Iend+1,Jend  ,knew))
!     *                               *DBLE(mask_u(Iend+1,Jend+1))

      tmp=0.5d0*DBLE(mask_u(Iend+1,Jend+1))*ad_ubar(Iend+1,Jend+1,knew)
      ad_ubar(Iend  ,Jend+1,knew)=ad_ubar(Iend  ,Jend+1,knew)+tmp
      ad_ubar(Iend+1,Jend  ,knew)=ad_ubar(Iend+1,Jend  ,knew)+tmp
      ad_ubar(Iend+1,Jend+1,knew)=0.d0

! tl line 1:      tl_ubar(Istr,Jend+1,knew)=0.5d0*
!     *                                (tl_ubar(Istr  ,Jend  ,knew)+
!     *                                 tl_ubar(Istr+1,Jend+1,knew))
!     *                               *DBLE(mask_u(Istr,Jend+1))

      tmp=0.5d0*DBLE(mask_u(Istr,Jend+1))*ad_ubar(Istr,Jend+1,knew)
      ad_ubar(Istr  ,Jend  ,knew)=ad_ubar(Istr  ,Jend  ,knew)+tmp
      ad_ubar(Istr+1,Jend+1,knew)=ad_ubar(Istr+1,Jend+1,knew)+tmp
      ad_ubar(Istr,Jend+1,knew)=0.d0

!--------------------------------------
! ADJ to U2: SOUTH (tangential, no gradient)
!--------------------------------------
      DO i=Istr+1,Iend
! TL:       tl_ubar(i,Jstr-1,knew)=tl_ubar(i,Jstr,knew)
!     *                        *DBLE(mask_u(i,Jstr-1))
       ad_ubar(i,Jstr,knew)=ad_ubar(i,Jstr,knew)+
     *    ad_ubar(i,Jstr-1,knew)*DBLE(mask_u(i,Jstr-1))
       ad_ubar(i,Jstr-1,knew)=0.d0
      END DO      

!--------------------------------------
! ADJ to U2: NORTH (tangential, no gradient)
!--------------------------------------
      DO i=Istr+1,Iend
! TL:       tl_ubar(i,Jend+1,knew)=tl_ubar(i,Jend,knew)
!     *                        *DBLE(mask_u(i,Jend+1))
       ad_ubar(i,Jend,knew)=ad_ubar(i,Jend,knew)+
     *    ad_ubar(i,Jend+1,knew)*DBLE(mask_u(i,Jend+1))
       ad_ubar(i,Jend+1,knew)=0.d0
      END DO      

# endif /* not defined NS_PERIODIC */

!--------------------------------------
! ADJ to U2: EAST
!--------------------------------------

#if defined EASTERN_WALL
      DO j=Jstr,Jend
! TL:       tl_ubar(Iend+1,j,knew)=0.d0
       ad_ubar(Iend+1,j,knew)=0.d0
      END DO      

#else
# ifdef EAST_M2FLATHER
      DO j=Jstr,Jend
       cff=0.5d0*(h(Iend,j)+h(Iend+1,j))

! TL:       tl_ubar(Iend+1,j,knew)= 
!     *    DSQRT(g/cff)
!     *   *0.5d0*(tl_zeta(Iend,j,know)+tl_zeta(Iend+1,j,know))
!     *   *DBLE(mask_u(Iend+1,j))

       tmp=DSQRT(g/cff)*0.5d0*DBLE(mask_u(Iend+1,j))
     *                       *ad_ubar(Iend+1,j,knew)

       ad_zeta(Iend,j,know)=ad_zeta(Iend,j,know)    +tmp
       ad_zeta(Iend+1,j,know)=ad_zeta(Iend+1,j,know)+tmp
       ad_ubar(Iend+1,j,knew)=0.d0

      END DO      

# else
      M2, EASTERN BOUNDARY, UNDEFINED OPTION

# endif
#endif /* defined or or not defined  EASTERN_WALL */

!--------------------------------------
! ADJ to U2: WEST
!--------------------------------------
#if defined WESTERN_WALL
      DO j=Jstr,Jend
! TL:       tl_ubar(Istr,j,knew)=0.d0
       ad_ubar(Istr,j,knew)=0.d0
      END DO      

#else
# ifdef WEST_M2FLATHER
      DO j=Jstr,Jend
       cff=0.5d0*(h(Istr-1,j)+h(Istr,j))

! TL       tl_ubar(Istr,j,knew)=
!     *          -DSQRT(g/cff)
!     *          *0.5d0*(tl_zeta(Istr-1,j,know)+tl_zeta(Istr,j,know))
!     *          *DBLE(mask_u(Istr,j))
  
       tmp=-DSQRT(g/cff)*0.5d0*DBLE(mask_u(Istr,j))*ad_ubar(Istr,j,knew)

       ad_zeta(Istr-1,j,know)=ad_zeta(Istr-1,j,know)+tmp
       ad_zeta(Istr  ,j,know)=ad_zeta(Istr  ,j,know)+tmp
       ad_ubar(Istr,j,knew)=0.d0

      END DO      
      
# else
      M2, WESTERN BOUNDARY, UNDEFINED OPTION

# endif
#endif /* defined or or not defined  WESTERN_WALL */


      RETURN
      END

!****************************************************************
      SUBROUTINE ad_v2dbc(Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj,  
     *               krhs, kstp, knew,
     *               iif,PREDICTOR,g,h,mask_v,
     *               ad_vbar,ad_zeta)

      implicit none
      integer, intent(in) :: Istr,Iend,Jstr,Jend,LBi,UBi,LBj,UBj      
      integer, intent(in) :: krhs, kstp, knew,iif   
      logical, intent(in) :: PREDICTOR
      integer, intent(in) :: mask_v(LBi:UBi,LBj:UBj)
      real(8), intent(in) :: h(LBi:UBi,LBj:UBj)
      real(8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
      real(8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
      real(8), intent(in) :: g

      real(8) :: cff,tmp

      integer :: i,j,know

      IF (iif.eq.1) THEN
        know=krhs
      ELSE IF (PREDICTOR) THEN
        know=krhs
      ELSE
        know=kstp
      END IF

#ifndef NS_PERIODIC

!------------
! ADJ to Corners
!------------

! tl line 4:      tl_vbar(Iend+1,Jstr,knew)=
!     *  0.5d0
!     * *(tl_vbar(Iend+1,Jstr+1,knew)+tl_vbar(Iend,Jstr,knew))
!     * *DBLE(mask_v(Iend+1,Jstr))

      tmp=0.5d0*DBLE(mask_v(Iend+1,Jstr))*ad_vbar(Iend+1,Jstr,knew)
      ad_vbar(Iend+1,Jstr+1,knew)=ad_vbar(Iend+1,Jstr+1,knew)+tmp
      ad_vbar(Iend  ,Jstr  ,knew)=ad_vbar(Iend  ,Jstr  ,knew)+tmp
      ad_vbar(Iend+1,Jstr,knew)=0.d0
 
! tl line 3:      tl_vbar(Istr-1,Jstr,knew)=
!     *  0.5d0
!     * *(tl_vbar(Istr-1,Jstr+1,knew)+tl_vbar(Istr,Jstr,knew))
!     * *DBLE(mask_v(Istr-1,Jstr))

      tmp=0.5d0*DBLE(mask_v(Istr-1,Jstr))*ad_vbar(Istr-1,Jstr,knew)
      ad_vbar(Istr-1,Jstr+1,knew)=ad_vbar(Istr-1,Jstr+1,knew)+tmp
      ad_vbar(Istr  ,Jstr  ,knew)=ad_vbar(Istr  ,Jstr  ,knew)+tmp
      ad_vbar(Istr-1,Jstr,knew)=0.d0
 
! tl line 2:      tl_vbar(Iend+1,Jend+1,knew)=
!     *  0.5d0
!     * *(tl_vbar(Iend,Jend+1,knew)+tl_vbar(Iend+1,Jend,knew))
!     * *DBLE(mask_v(Iend+1,Jend+1))
      
      tmp=0.5d0*DBLE(mask_v(Iend+1,Jend+1))*ad_vbar(Iend+1,Jend+1,knew)
      ad_vbar(Iend,Jend+1,knew)=ad_vbar(Iend,Jend+1,knew)+tmp
      ad_vbar(Iend+1,Jend,knew)=ad_vbar(Iend+1,Jend,knew)+tmp
      ad_vbar(Iend+1,Jend+1,knew)=0.d0
 
! tl line 1:      tl_vbar(Istr-1,Jend+1,knew)=
!     *  0.5d0
!     * *(tl_vbar(Istr-1,Jend,knew)+tl_vbar(Istr,Jend+1,knew))
!     * *DBLE(mask_v(Istr-1,Jend+1))

      tmp=0.5d0*DBLE(mask_v(Istr-1,Jend+1))*ad_vbar(Istr-1,Jend+1,knew)
      ad_vbar(Istr-1,Jend,knew)=ad_vbar(Istr-1,Jend,knew)+tmp
      ad_vbar(Istr,Jend+1,knew)=ad_vbar(Istr,Jend+1,knew)+tmp
      ad_vbar(Istr-1,Jend+1,knew)=0.d0
 
!--------------------------------------
! ADJ to V2: SOUTH
!--------------------------------------
# if defined SOUTHERN_WALL
      DO i=Istr,Iend
! TL:       tl_vbar(i,Jstr,knew)=0.d0
       ad_vbar(i,Jstr,knew)=0.d0
      END DO

# else
#  if defined SOUTH_M2FLATHER
      DO i=Istr,Iend
       cff=0.5d0*(h(i,Jstr)+h(i,Jstr-1))

! TL:       tl_vbar(i,Jstr,knew)=
!     *    -DSQRT(g/cff)
!     *   *0.5d0*(tl_zeta(i,Jstr,know)+tl_zeta(i,Jstr-1,know))
!     *   *DBLE(mask_v(i,Jstr))

       tmp=-DSQRT(g/cff)*0.5d0*DBLE(mask_v(i,Jstr))*ad_vbar(i,Jstr,knew)
       ad_zeta(i,Jstr,know)=ad_zeta(i,Jstr,know)+tmp
       ad_zeta(i,Jstr-1,know)=ad_zeta(i,Jstr-1,know)+tmp
       ad_vbar(i,Jstr,knew)=0.d0

      END DO
#  else
      M2, SOUTHERN BOUNDARY, UNDEFINED OPTION
#  endif /* SOUTH_M2FLATHER */
# endif /* SOUTHERN_WALL */

!--------------------------------------
! ADJ to V2: NORTH
!--------------------------------------
# if defined NORTHERN_WALL
      DO i=Istr,Iend
! TL:       tl_vbar(i,Jend+1,knew)=0.d0
       ad_vbar(i,Jend+1,knew)=0.d0
      END DO

# else
#  if defined NORTH_M2FLATHER
      DO i=Istr,Iend

       cff=0.5d0*(h(i,Jend)+h(i,Jend+1))

! TL:       tl_vbar(i,Jend+1,knew)=
!     *    DSQRT(g/cff)
!     *   *0.5d0*(tl_zeta(i,Jend,know)+tl_zeta(i,Jend+1,know))
!     *   *DBLE(mask_v(i,Jend+1))

       tmp=DSQRT(g/cff)*0.5d0*DBLE(mask_v(i,Jend+1))
     *                       *ad_vbar(i,Jend+1,knew)
       ad_zeta(i,Jend,know)  =ad_zeta(i,Jend,know)  +tmp
       ad_zeta(i,Jend+1,know)=ad_zeta(i,Jend+1,know)+tmp
       ad_vbar(i,Jend+1,knew)=0.d0

      END DO

#  else
      M2, NORTHERN BOUNDARY, UNDEFINED OPTION
#  endif /* NORTH_M2FLATHER */
# endif /* NORTHERN_WALL */


#endif /* ifndef NS_PERIODIC */ 

!--------------------------------------
! ADJ to V2: EAST (tangential, no gradient)
!--------------------------------------
      DO j=Jstr+1,Jend
! TL:       tl_vbar(Iend+1,j,knew)=tl_vbar(Iend,j,knew)
!     *                        *DBLE(mask_v(Iend+1,j))
       
       ad_vbar(Iend,j,knew)=ad_vbar(Iend,j,knew)+
     *  DBLE(mask_v(Iend+1,j))*ad_vbar(Iend+1,j,knew)
       ad_vbar(Iend+1,j,knew)=0.d0

      END DO      

!--------------------------------------
! ADJ to V2: WEST (tangential, no gradient)
!--------------------------------------
      DO j=Jstr+1,Jend
! TL:       tl_vbar(Istr-1,j,knew)=tl_vbar(Istr,j,knew)
!     *                        *DBLE(mask_v(Istr-1,j))
       
       ad_vbar(Istr,j,knew)=ad_vbar(Istr,j,knew)+
     *  DBLE(mask_v(Istr-1,j))*ad_vbar(Istr-1,j,knew)
       ad_vbar(Istr-1,j,knew)=0.d0

      END DO      

      RETURN
      END
!=====================================================================
#endif /* MPI or SERIAL */
!=====================================================================

#endif /* ADJOINT */


