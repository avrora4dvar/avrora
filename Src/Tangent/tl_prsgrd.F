#include "cppdefs.h"

#ifdef TANGENT

#ifndef PRSGRD_BT_ONLY

# ifdef DJ_GRADPS   /* NOTE: NOT FINISHED !!! */
      SUBROUTINE tl_prsgrd

      USE mod_ocean
      implicit none

      integer :: i,j,k

      real(8) :: GRho,Grho0,HalfGRho
      real(8) :: dR(LBi:UBi,0:N)
      real(8) :: dZ(LBi:UBi,0:N)
      real(8) :: aux(LBi:UBi,LBj:UBj)
      real(8) :: FC(LBi:UBi,LBj:UBj)
      real(8) :: dRx(LBi:UBi,LBj:UBj)
      real(8) :: dZx(LBi:UBi,LBj:UBj)
      real(8) :: P(LBi:UBi,LBj:UBj,N)

      real(8) :: cff,cff1,cff2

      real(8), parameter :: OneFifth = 0.2d0
      real(8), parameter :: OneTwelfth = 1.0d0/12.d0
      real(8), parameter :: eps = 1.0d-10

!-----------------------------------------------------------------------
!  Preliminary step (same for XI- and ETA-components:
!-----------------------------------------------------------------------
!
      GRho=g/rho0
      GRho0=1000.0d0*GRho
      HalfGRho=0.5d0*GRho

      DO j=JstrV-1,Jend

       DO k=1,N-1
        DO i=IstrU-1,Iend
! NL:
         dR(i,k)=rho(i,j,k+1)-rho(i,j,k)
         dZ(i,k)=z_r(i,j,k+1)-z_r(i,j,k)

! TL:
         tl_dR(i,k)=tl_rho(i,j,k+1)-tl_rho(i,j,k)
         tl_dZ(i,k)=tl_z_r(i,j,k+1)-tl_z_r(i,j,k)

        END DO
       END DO

       DO i=IstrU-1,Iend
! NL:
        dR(i,N)=dR(i,N-1)
        dZ(i,N)=dZ(i,N-1)
        dR(i,0)=dR(i,1)
        dZ(i,0)=dZ(i,1)

! TL:
        tl_dR(i,N)=tl_dR(i,N-1)
        tl_dZ(i,N)=tl_dZ(i,N-1)
        tl_dR(i,0)=tl_dR(i,1)
        tl_dZ(i,0)=tl_dZ(i,1)
       END DO

       DO k=N,1,-1
        DO i=IstrU-1,Iend
! NL:
         cff=2.0d0*dR(i,k)*dR(i,k-1)
! TL:
         tl_cff=2.0d0*(tl_dR(i,k)*dR(i,k-1)+dR(i,k)*tl_dR(i,k-1))

! ORIGINAL NL (note dR, including its units, is changes recursively):
!         IF (cff.gt.eps) THEN
!          dR(i,k)=cff/(dR(i,k)+dR(i,k-1))
!         ELSE
!          dR(i,k)=0.0d0
!         END IF
!         dZ(i,k)=2.0d0*dZ(i,k)*dZ(i,k-1)/(dZ(i,k)+dZ(i,k-1))

! REWRITE NL TO AVOID RECURS STATEMENT:
!         IF (cff.gt.eps) THEN
!          dR1=cff/(dR(i,k)+dR(i,k-1))
!         ELSE
!          dR1=0.0d0
!         END IF
!         dR(i,k)=dR1
!         dZ1=2.0d0*dZ(i,k)*dZ(i,k-1)/(dZ(i,k)+dZ(i,k-1))
!         dZ(i,k)=dZ1

! TL TO THE MODIFIED FORM OF NL, AS ABOVE:
!
         IF (cff.gt.eps) THEN
          otmp=1.d0/(dR(i,k)+dR(i,k-1))
          dR1=cff*otmp
! note: utilize TL for (dR(i,k)+dR(i,k-1))*dR1=cff:
          tl_dR1=(tl_cff-(tl_dR(i,k)+tl_dR(i,k-1))*dR1)*otmp
         ELSE
          dR1=0.0d0
          tl_dR1=0.0d0
         END IF

         dR(i,k)=dR1
         tl_dR(i,k)=tl_dR1

         otmp=1.d0/(dZ(i,k)+dZ(i,k-1))
         dZ1=2.0d0*dZ(i,k)*dZ(i,k-1)*otmp
         tl_dZ1=( 2.0d0*(tl_dZ(i,k)*dZ(i,k-1)+dZ(i,k)*tl_dZ(i,k-1))
                 -dZ1*(tl_dZ(i,k)+tl_dZ(i,k-1)) )*otmp  

         dZ(i,k)=dZ1
         tl_dZ(i,k)=tl_dZ1

        END DO
       END DO

       DO i=IstrU-1,Iend
        cff1=1.0d0/(z_r(i,j,N)-z_r(i,j,N-1))
        tl_cff1=-(tl_z_r(i,j,N)-tl_z_r(i,j,N-1))*cff1*cff1

! ORIGINAL NL cff2:    
        cff2=0.5d0*(rho(i,j,N)-rho(i,j,N-1))*                
     *         (z_w(i,j,N)-z_r(i,j,N))*cff1
! MODIFIED:
        tmp1=rho(i,j,N)-rho(i,j,N-1)
        tmp2=z_w(i,j,N)-z_r(i,j,N)
        cff2=0.5d0*tmp1*tmp2*cff1

        tl_tmp1=tl_rho(i,j,N)-tl_rho(i,j,N-1)
        tl_tmp2=tl_z_w(i,j,N)-tl_z_r(i,j,N)
        tl_cff2=0.5d0*(tl_tmp1*   tmp2*   cff1+
     *                    tmp1*tl_tmp2*   cff1+
     *                    tmp1*   tmp2*tl_cff1)

        P(i,j,N)= GRho0*z_w(i,j,N)                            
     *           +GRho*(rho(i,j,N)+cff2)*(z_w(i,j,N)-z_r(i,j,N))

        tl_P(i,j,N)= GRho0*tl_z_w(i,j,N)                            
     *              +GRho*(
     *                (tl_rho(i,j,N)+tl_cff2)*(z_w(i,j,N)-z_r(i,j,N))
     *               +(rho(i,j,N)+cff2)*(tl_z_w(i,j,N)-tl_z_r(i,j,N))
     *                    )
       END DO

       DO k=N-1,1,-1
        DO i=IstrU-1,Iend
         P(i,j,k)=P(i,j,k+1)                                        
     *             +HalfGRho*((rho(i,j,k+1)+rho(i,j,k))*               
     *                         (z_r(i,j,k+1)-z_r(i,j,k))-               
     *                         OneFifth*                                
     *                         ((dR(i,k+1)-dR(i,k))*                    
     *                          (z_r(i,j,k+1)-z_r(i,j,k)-               
     *                           OneTwelfth*                            
     *                           (dZ(i,k+1)+dZ(i,k)))-                  
     *                          (dZ(i,k+1)-dZ(i,k))*                    
     *                          (rho(i,j,k+1)-rho(i,j,k)-               
     *                           OneTwelfth*                            
     *                           (dR(i,k+1)+dR(i,k)))))

         tl_P(i,j,k)=tl_P(i,j,k+1)                                        
     *    +HalfGRho*(
     *      (tl_rho(i,j,k+1)+tl_rho(i,j,k))*(z_r(i,j,k+1)-z_r(i,j,k))      
     *     +(rho(i,j,k+1)+rho(i,j,k))*(tl_z_r(i,j,k+1)-tl_z_r(i,j,k))      
     *     -OneFifth*(                                
     *        (tl_dR(i,k+1)-tl_dR(i,k))*
     *        (z_r(i,j,k+1)-z_r(i,j,k)-OneTwelfth*(dZ(i,k+1)+dZ(i,k)))
     *       +(dR(i,k+1)-dR(i,k))*
     *        (tl_z_r(i,j,k+1)-tl_z_r(i,j,k)
     *                          -OneTwelfth*(tl_dZ(i,k+1)+tl_dZ(i,k)))
     *       -(tl_dZ(i,k+1)-tl_dZ(i,k))*                    
     *        (rho(i,j,k+1)-rho(i,j,k)-OneTwelfth*(dR(i,k+1)+dR(i,k)))
     *       -(dZ(i,k+1)-dZ(i,k))*                    
     *        (tl_rho(i,j,k+1)-tl_rho(i,j,k)
     *                          -OneTwelfth*(tl_dR(i,k+1)+tl_dR(i,k)))
     *                )
     *              )
        END DO
       END DO

      END DO  ! <- j=JstrV-1,Jend

      WRITE (*,*) 'STOP in tl_prsgrd: option DJ_GRADPS not finished'
      STOP
!
!-----------------------------------------------------------------------
!  Compute XI-component pressure gradient term.
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!  ETA-component pressure gradient term.
!-----------------------------------------------------------------------
!

      RETURN
      END

# endif /* DJ_GRADPS, NOT FINISHED */

!=============================================================
# ifdef DJ_STD2 /* STANDARD 2nd ORDER */

      SUBROUTINE tl_prsgrd

      USE mod_ocean
      implicit none
      integer :: i,j,k

      real(8) :: fac1,fac2,fac3,cff1,cff2,cff3,cff4
      real(8) :: tl_cff1,tl_cff2,tl_cff3,tl_cff4
      real(8) :: phix(LBi:UBi)
      real(8) :: phie(LBi:UBi)
      real(8) :: tl_phix(LBi:UBi)
      real(8) :: tl_phie(LBi:UBi)

!
!-----------------------------------------------------------------------
!  Calculate pressure gradient in the XI-direction (m4/s2).
!-----------------------------------------------------------------------
!
!  Compute surface baroclinic pressure gradient. (XI)
!
      fac1=0.5d0*g/rho0
      fac2=1000.0d0*g/rho0
      fac3=0.25d0*g/rho0

      DO j=Jstr,Jend
       DO i=IstrU,Iend
        cff1=z_w(i  ,j,N)-z_r(i  ,j,N)+ 
     *       z_w(i-1,j,N)-z_r(i-1,j,N)
        tl_cff1=tl_z_w(i  ,j,N)-tl_z_r(i  ,j,N)+ 
     *          tl_z_w(i-1,j,N)-tl_z_r(i-1,j,N)

        phix(i)=fac1*(rho(i,j,N)-rho(i-1,j,N))*cff1
        tl_phix(i)=fac1*( (tl_rho(i,j,N)-tl_rho(i-1,j,N))*cff1+
     *                    (rho(i,j,N)-rho(i-1,j,N))*tl_cff1 )

        phix(i)=phix(i)+
     *            (fac2+fac1*(rho(i,j,N)+rho(i-1,j,N)))*
     *            (z_w(i,j,N)-z_w(i-1,j,N))
        tl_phix(i)=tl_phix(i)+
     *            fac1*(tl_rho(i,j,N)+tl_rho(i-1,j,N))*
     *            (z_w(i,j,N)-z_w(i-1,j,N))
     *           +(fac2+fac1*(rho(i,j,N)+rho(i-1,j,N)))*
     *            (tl_z_w(i,j,N)-tl_z_w(i-1,j,N))

        tl_ru(i,j,N,nrhs)=-0.5d0*dy_u(i,j)*
     *   ( (tl_Hz(i,j,N)+tl_Hz(i-1,j,N))*phix(i)+
     *     (Hz(i,j,N)+Hz(i-1,j,N))*tl_phix(i) )
       END DO

!
!  Compute interior baroclinic pressure gradient.  Differentiate and
!  then vertically integrate. (XI)
!
       DO k=N-1,1,-1
        DO i=IstrU,Iend
         cff1=rho(i,j,k+1)-rho(i-1,j,k+1)+
     *        rho(i,j,k  )-rho(i-1,j,k  )
         cff2=rho(i,j,k+1)+rho(i-1,j,k+1)-
     *        rho(i,j,k  )-rho(i-1,j,k  )
         cff3=z_r(i,j,k+1)+z_r(i-1,j,k+1)-
     *        z_r(i,j,k  )-z_r(i-1,j,k  )
         cff4=z_r(i,j,k+1)-z_r(i-1,j,k+1)+
     *        z_r(i,j,k  )-z_r(i-1,j,k  )

         tl_cff1=tl_rho(i,j,k+1)-tl_rho(i-1,j,k+1)+
     *           tl_rho(i,j,k  )-tl_rho(i-1,j,k  )
         tl_cff2=tl_rho(i,j,k+1)+tl_rho(i-1,j,k+1)-
     *           tl_rho(i,j,k  )-tl_rho(i-1,j,k  )
         tl_cff3=tl_z_r(i,j,k+1)+tl_z_r(i-1,j,k+1)-
     *           tl_z_r(i,j,k  )-tl_z_r(i-1,j,k  )
         tl_cff4=tl_z_r(i,j,k+1)-tl_z_r(i-1,j,k+1)+
     *           tl_z_r(i,j,k  )-tl_z_r(i-1,j,k  )

         phix(i)=phix(i)+fac3*(cff1*cff3-cff2*cff4)
         tl_phix(i)=tl_phix(i)+fac3*( tl_cff1*cff3+cff1*tl_cff3
     *                               -tl_cff2*cff4-cff2*tl_cff4)

         tl_ru(i,j,k,nrhs)=-0.5d0*dy_u(i,j)*
     *     ( (tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*phix(i)+
     *       (Hz(i,j,k)+Hz(i-1,j,k))*tl_phix(i) )
        END DO
       END DO
!
!-----------------------------------------------------------------------
!  Calculate pressure gradient in the ETA-direction (m4/s2).
!-----------------------------------------------------------------------
!
!  Compute surface baroclinic pressure gradient. (ETA)
!
       IF (j.ge.JstrV) THEN
        DO i=Istr,Iend

         cff1=z_w(i,j  ,N)-z_r(i,j  ,N)+
     *        z_w(i,j-1,N)-z_r(i,j-1,N)
         tl_cff1=tl_z_w(i,j  ,N)-tl_z_r(i,j  ,N)+
     *           tl_z_w(i,j-1,N)-tl_z_r(i,j-1,N)

         phie(i)=fac1*(rho(i,j,N)-rho(i,j-1,N))*cff1
         tl_phie(i)=fac1*
     *     ( (tl_rho(i,j,N)-tl_rho(i,j-1,N))*cff1+
     *       (rho(i,j,N)-rho(i,j-1,N))*tl_cff1 )

         phie(i)=phie(i)+
     *              (fac2+fac1*(rho(i,j,N)+rho(i,j-1,N)))*
     *              (z_w(i,j,N)-z_w(i,j-1,N))
         tl_phie(i)= tl_phie(i)
     *              +fac1*(tl_rho(i,j,N)+tl_rho(i,j-1,N))*
     *                (z_w(i,j,N)-z_w(i,j-1,N))
     *              +(fac2+fac1*(rho(i,j,N)+rho(i,j-1,N)))*
     *                (tl_z_w(i,j,N)-tl_z_w(i,j-1,N))

         tl_rv(i,j,N,nrhs)=-0.5d0*dx_v(i,j)*
     *     ( (tl_Hz(i,j,N)+tl_Hz(i,j-1,N))*phie(i)
     *      +(Hz(i,j,N)+Hz(i,j-1,N))*tl_phie(i) )
        END DO
!
!  Compute interior baroclinic pressure gradient.  Differentiate and
!  then vertically integrate. (ETA)
!
        DO k=N-1,1,-1
         DO i=Istr,Iend
          cff1=rho(i,j,k+1)-rho(i,j-1,k+1)+
     *         rho(i,j,k  )-rho(i,j-1,k  )
          cff2=rho(i,j,k+1)+rho(i,j-1,k+1)-
     *         rho(i,j,k  )-rho(i,j-1,k  )
          cff3=z_r(i,j,k+1)+z_r(i,j-1,k+1)-
     *         z_r(i,j,k  )-z_r(i,j-1,k  )
          cff4=z_r(i,j,k+1)-z_r(i,j-1,k+1)+
     *         z_r(i,j,k  )-z_r(i,j-1,k  )

          tl_cff1=tl_rho(i,j,k+1)-tl_rho(i,j-1,k+1)+
     *         tl_rho(i,j,k  )-tl_rho(i,j-1,k  )
          tl_cff2=tl_rho(i,j,k+1)+tl_rho(i,j-1,k+1)-
     *         tl_rho(i,j,k  )-tl_rho(i,j-1,k  )
          tl_cff3=tl_z_r(i,j,k+1)+tl_z_r(i,j-1,k+1)-
     *         tl_z_r(i,j,k  )-tl_z_r(i,j-1,k  )
          tl_cff4=tl_z_r(i,j,k+1)-tl_z_r(i,j-1,k+1)+
     *         tl_z_r(i,j,k  )-tl_z_r(i,j-1,k  )

          phie(i)=phie(i)+fac3*(cff1*cff3-cff2*cff4)
          tl_phie(i)=tl_phie(i)+fac3*( tl_cff1*cff3+cff1*tl_cff3
     *                                -tl_cff2*cff4-cff2*tl_cff4)

          tl_rv(i,j,k,nrhs)=-0.5d0*dx_v(i,j)* 
     *      ( (tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*phie(i)
     *       +(Hz(i,j,k)+Hz(i,j-1,k))*tl_phie(i) )
         END DO
        END DO
       END IF  ! (j.ge.JstrV)
 
      END DO

      RETURN
      END

# endif /* DJ_STD2 */


#endif /* not defined PRSGRD_BT_ONLY */
!=============================================================

#ifdef PRSGRD_BT_ONLY
      SUBROUTINE tl_prsgrd

      USE mod_ocean
      implicit none
      integer :: i,j,k

      real(8) :: P(LBi:UBi,LBj:UBj,N)
      real(8) :: tl_P(LBi:UBi,LBj:UBj,N)
      real(8) :: GRho,Grho0,HalfGRho

      GRho=g/rho0
      GRho0=1000.0d0*GRho

      DO j=JstrV-1,Jend

       DO i=IstrU-1,Iend
        P(i,j,N)=GRho0*z_w(i,j,N)                            
        tl_P(i,j,N)=GRho0*tl_z_w(i,j,N)                            
       END DO

       DO k=N-1,1,-1
        DO i=IstrU-1,Iend
         P(i,j,k)=P(i,j,k+1)     
         tl_P(i,j,k)=tl_P(i,j,k+1)                                        
        END DO
       END DO

      END DO  ! <- j=JstrV-1,Jend

!-----------------------------------------------------------------------
!  Compute XI-component pressure gradient term.
!-----------------------------------------------------------------------
!
      DO k=N,1,-1
!
       DO j=Jstr,Jend
        DO i=IstrU,Iend
! NL:         ru(i,j,k,nrhs)=dy_u(i,j)*Hz_u(i,j,k)*
! NL:     *                     (P(i-1,j,k)-P(i,j,k))  
         tl_ru(i,j,k,nrhs)=dy_u(i,j)*(
     *             tl_Hz_u(i,j,k)*(P(i-1,j,k)-P(i,j,k)) +  
     *             Hz_u(i,j,k)*( tl_P(i-1,j,k)-tl_P(i,j,k) )
     *                               )
        END DO
       END DO

      END DO ! <- "k=N,1,-1"

!
!-----------------------------------------------------------------------
!  ETA-component pressure gradient term.
!-----------------------------------------------------------------------
!
      DO k=N,1,-1

       DO j=JstrV,Jend
        DO i=Istr,Iend
! NL:         rv(i,j,k,nrhs)=dx_v(i,j)*Hz_v(i,j,k)*
! NL:     *                       ( P(i,j-1,k)-P(i,j,k) ) 
         tl_rv(i,j,k,nrhs)=dx_v(i,j)*(
     *             tl_Hz_v(i,j,k)*( P(i,j-1,k)-P(i,j,k) ) +
     *             Hz_v(i,j,k)*( tl_P(i,j-1,k)-tl_P(i,j,k) ) 
     *                               )
        END DO
       END DO

      END DO ! k=N,1,-1, ETA comp

      RETURN
      END 
#endif /* PRSGRD_BT_ONLY */


#endif /* TANGENT */
