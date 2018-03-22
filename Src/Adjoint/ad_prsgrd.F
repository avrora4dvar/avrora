#include "cppdefs.h"

#ifdef ADJOINT

# ifndef PRSGRD_BT_ONLY

#  ifdef DJ_GRADPS   /* NOTE: NOT YET !!! */

#  endif

!=============================================================
#  ifdef DJ_STD2 /* STANDARD 2nd ORDER */
      SUBROUTINE ad_prsgrd

      USE mod_ocean
      implicit none
      integer :: i,j,k

      real(8) :: fac1,fac2,fac3,tmp,tmp1,tmp2
      real(8) :: ad_cff1,ad_cff2,ad_cff3,ad_cff4
      real(8) :: cff1(LBi:UBi,N)
      real(8) :: cff2(LBi:UBi,N)
      real(8) :: cff3(LBi:UBi,N)
      real(8) :: cff4(LBi:UBi,N)
      real(8) :: phix(LBi:UBi,N)
      real(8) :: phie(LBi:UBi,N)
      real(8) :: ad_phix(LBi:UBi)
      real(8) :: ad_phie(LBi:UBi)

      ad_cff1=0.d0
      ad_cff2=0.d0
      ad_cff3=0.d0
      ad_cff4=0.d0

      ad_phix=0.d0
      ad_phie=0.d0

      fac1=0.5d0*g/rho0
      fac2=1000.0d0*g/rho0
      fac3=0.25d0*g/rho0

      DO j=Jstr,Jend
!
!-----------------------------------------------------------------------
!  ADJ TO: Calculate pressure gradient in the ETA-direction (m4/s2).
!-----------------------------------------------------------------------

       IF (j.ge.JstrV) THEN

!
! Background (fwd) cff1,cff2,cff3,cff3 and phie:
! 

! - surface:
        DO i=Istr,Iend
         cff1(i,N)=z_w(i,j  ,N)-z_r(i,j  ,N)+
     *             z_w(i,j-1,N)-z_r(i,j-1,N)
         phie(i,N)=fac1*(rho(i,j,N)-rho(i,j-1,N))*cff1(i,N)
         phie(i,N)=phie(i,N)+
     *              (fac2+fac1*(rho(i,j,N)+rho(i,j-1,N)))*
     *              (z_w(i,j,N)-z_w(i,j-1,N))
        END DO

! - interior:
        DO k=N-1,1,-1
         DO i=Istr,Iend
          cff1(i,k)=rho(i,j,k+1)-rho(i,j-1,k+1)+
     *              rho(i,j,k  )-rho(i,j-1,k  )
          cff2(i,k)=rho(i,j,k+1)+rho(i,j-1,k+1)-
     *              rho(i,j,k  )-rho(i,j-1,k  )
          cff3(i,k)=z_r(i,j,k+1)+z_r(i,j-1,k+1)-
     *              z_r(i,j,k  )-z_r(i,j-1,k  )
          cff4(i,k)=z_r(i,j,k+1)-z_r(i,j-1,k+1)+
     *              z_r(i,j,k  )-z_r(i,j-1,k  )

          phie(i,k)=phie(i,k+1)+
     *     fac3*(cff1(i,k)*cff3(i,k)-cff2(i,k)*cff4(i,k))
         END DO
        END DO

!
!  ADJ TO: Compute interior baroclinic pressure gradient (ETA)  
!  Differentiate and then vertically integrate.
!  
! NOTE (AK): cff1..4 and phie are computed above and saved as 2D arrays

        DO k=1,N-1        ! direction switched 
         DO i=Istr,Iend

! tl line 6:          tl_rv(i,j,k,nrhs)=-0.5d0*dx_v(i,j)* 
! tl line 6:     *      ( (tl_Hz(i,j,k)+tl_Hz(i,j-1,k))*phie(i)
! tl line 6:     *       +(Hz(i,j,k)+Hz(i,j-1,k))*tl_phie(i) )

          tmp=0.5d0*dx_v(i,j)*ad_rv(i,j,k,nrhs)

          ad_Hz(i,j,k)=ad_Hz(i,j,k)-tmp*phie(i,k)          
          ad_Hz(i,j-1,k)=ad_Hz(i,j-1,k)-tmp*phie(i,k)  
          ad_phie(i)=ad_phie(i)-tmp*(Hz(i,j,k)+Hz(i,j-1,k))         

          ad_rv(i,j,k,nrhs)=0.d0

! tl line 5:          tl_phie(i)=tl_phie(i)+fac3*( tl_cff1*cff3+cff1*tl_cff3
! tl line 5:     *                                -tl_cff2*cff4-cff2*tl_cff4)
          tmp=ad_phie(i)*fac3
          ad_cff1=ad_cff1+tmp*cff3(i,k)
          ad_cff3=ad_cff3+tmp*cff1(i,k)
          ad_cff2=ad_cff2-tmp*cff4(i,k)
          ad_cff4=ad_cff4-tmp*cff2(i,k)

! tl line 4:          tl_cff4=tl_z_r(i,j,k+1)-tl_z_r(i,j-1,k+1)+
! tl line 4:     *         tl_z_r(i,j,k  )-tl_z_r(i,j-1,k  )
          ad_z_r(i,j,k+1)=ad_z_r(i,j,k+1)    +ad_cff4
          ad_z_r(i,j-1,k+1)=ad_z_r(i,j-1,k+1)-ad_cff4
          ad_z_r(i,j,k  )=ad_z_r(i,j,k  )    +ad_cff4
          ad_z_r(i,j-1,k  )=ad_z_r(i,j-1,k  )-ad_cff4
          ad_cff4=0.d0

! tl line 3:          tl_cff3=tl_z_r(i,j,k+1)+tl_z_r(i,j-1,k+1)-
! tl line 3:     *         tl_z_r(i,j,k  )-tl_z_r(i,j-1,k  )
          ad_z_r(i,j,k+1)=ad_z_r(i,j,k+1)    +ad_cff3
          ad_z_r(i,j-1,k+1)=ad_z_r(i,j-1,k+1)+ad_cff3
          ad_z_r(i,j,k  )=ad_z_r(i,j,k  )    -ad_cff3
          ad_z_r(i,j-1,k  )=ad_z_r(i,j-1,k  )-ad_cff3
          ad_cff3=0.d0

! tl line 2:          tl_cff2=tl_rho(i,j,k+1)+tl_rho(i,j-1,k+1)-
! tl line 2:     *         tl_rho(i,j,k  )-tl_rho(i,j-1,k  )
          ad_rho(i,j,k+1)=ad_rho(i,j,k+1)    +ad_cff2
          ad_rho(i,j-1,k+1)=ad_rho(i,j-1,k+1)+ad_cff2
          ad_rho(i,j,k  )=ad_rho(i,j,k  )    -ad_cff2
          ad_rho(i,j-1,k  )=ad_rho(i,j-1,k  )-ad_cff2
          ad_cff2=0.d0
          
! tl line 1:          tl_cff1=tl_rho(i,j,k+1)-tl_rho(i,j-1,k+1)+
! tl line 1:     *         tl_rho(i,j,k  )-tl_rho(i,j-1,k  )
          ad_rho(i,j,k+1)=ad_rho(i,j,k+1)    +ad_cff1
          ad_rho(i,j-1,k+1)=ad_rho(i,j-1,k+1)-ad_cff1
          ad_rho(i,j,k  )=ad_rho(i,j,k  )    +ad_cff1
          ad_rho(i,j-1,k  )=ad_rho(i,j-1,k  )-ad_cff1
          ad_cff1=0.d0

         END DO
        END DO

!
!  ADJ TO: Compute surface baroclinic pressure gradient (ETA)
!  (background variables are computed above and stored in 2D arrays)
!
        DO i=Istr,Iend

! tl line 4:         tl_rv(i,j,N,nrhs)=-0.5d0*dx_v(i,j)*
! tl line 4:     *     ( (tl_Hz(i,j,N)+tl_Hz(i,j-1,N))*phie(i)
! tl line 4:     *      +(Hz(i,j,N)+Hz(i,j-1,N))*tl_phie(i) )
         tmp=0.5d0*dx_v(i,j)*ad_rv(i,j,N,nrhs)
         ad_Hz(i,j,N)=ad_Hz(i,j,N)-phie(i,N)*tmp
         ad_Hz(i,j-1,N)=ad_Hz(i,j-1,N)-phie(i,N)*tmp
         ad_phie(i)=ad_phie(i)-(Hz(i,j,N)+Hz(i,j-1,N))*tmp
         ad_rv(i,j,N,nrhs)=0.d0

! tl line 3:         tl_phie(i)= tl_phie(i)
! tl line 3:     *              +fac1*(tl_rho(i,j,N)+tl_rho(i,j-1,N))*
! tl line 3:     *                (z_w(i,j,N)-z_w(i,j-1,N))
! tl line 3:     *              +(fac2+fac1*(rho(i,j,N)+rho(i,j-1,N)))*
! tl line 3:     *                (tl_z_w(i,j,N)-tl_z_w(i,j-1,N))

         tmp1=ad_phie(i)*fac1*(z_w(i,j,N)-z_w(i,j-1,N))        
         tmp2=ad_phie(i)*(fac2+fac1*(rho(i,j,N)+rho(i,j-1,N)))        

         ad_rho(i,j,N)=ad_rho(i,j,N)+tmp1
         ad_rho(i,j-1,N)=ad_rho(i,j-1,N)+tmp1
         ad_z_w(i,j,N)=ad_z_w(i,j,N)+tmp2
         ad_z_w(i,j-1,N)=ad_z_w(i,j-1,N)-tmp2

! tl line 2:         tl_phie(i)=fac1*
! tl line 2:     *     ( (tl_rho(i,j,N)-tl_rho(i,j-1,N))*cff1+
! tl line 2:     *       (rho(i,j,N)-rho(i,j-1,N))*tl_cff1 )

         tmp=ad_phie(i)*fac1*cff1(i,N)

         ad_rho(i,j,N)=ad_rho(i,j,N)+tmp
         ad_rho(i,j-1,N)=ad_rho(i,j-1,N)-tmp
         ad_cff1=ad_cff1+ad_phie(i)*fac1*(rho(i,j,N)-rho(i,j-1,N))

         ad_phie(i)=0.d0

! tl line 1:         tl_cff1=tl_z_w(i,j  ,N)-tl_z_r(i,j  ,N)+
! tl line 1:     *           tl_z_w(i,j-1,N)-tl_z_r(i,j-1,N)
 
         ad_z_w(i,j  ,N)=ad_z_w(i,j  ,N)+ad_cff1
         ad_z_r(i,j  ,N)=ad_z_r(i,j  ,N)-ad_cff1
         ad_z_w(i,j-1,N)=ad_z_w(i,j-1,N)+ad_cff1
         ad_z_r(i,j-1,N)=ad_z_r(i,j-1,N)-ad_cff1

         ad_cff1=0.d0

        END DO

       END IF   ! (j.ge.JstrV) 
!
!-----------------------------------------------------------------------
!  ADJ TO: Calculate pressure gradient in the XI-direction (m4/s2).
!-----------------------------------------------------------------------

!
! Background (fwd) cff1,cff2,cff3,cff3 and phix:
! 

! - surface:
       DO i=IstrU,Iend
        cff1(i,N)=z_w(i  ,j,N)-z_r(i  ,j,N)+ 
     *            z_w(i-1,j,N)-z_r(i-1,j,N)
        phix(i,N)=fac1*(rho(i,j,N)-rho(i-1,j,N))*cff1(i,N)
        phix(i,N)=phix(i,N)+
     *            (fac2+fac1*(rho(i,j,N)+rho(i-1,j,N)))*
     *            (z_w(i,j,N)-z_w(i-1,j,N))
      
       END DO

! - interior:
       DO k=N-1,1,-1
        DO i=IstrU,Iend
         cff1(i,k)=rho(i,j,k+1)-rho(i-1,j,k+1)+
     *             rho(i,j,k  )-rho(i-1,j,k  )
         cff2(i,k)=rho(i,j,k+1)+rho(i-1,j,k+1)-
     *             rho(i,j,k  )-rho(i-1,j,k  )
         cff3(i,k)=z_r(i,j,k+1)+z_r(i-1,j,k+1)-
     *             z_r(i,j,k  )-z_r(i-1,j,k  )
         cff4(i,k)=z_r(i,j,k+1)-z_r(i-1,j,k+1)+
     *             z_r(i,j,k  )-z_r(i-1,j,k  )
       
         phix(i,k)=phix(i,k+1)+
     *            fac3*(cff1(i,k)*cff3(i,k)-cff2(i,k)*cff4(i,k))
        END DO
       END DO

!
!  ADJ TO: Compute interior baroclinic pressure gradient (XI)  
!  Differentiate and then vertically integrate.
!  NOTE (AK): background cff1,cff2,cff3,cff4,phix computed above
!             and saved as 2d arrays
! 
       DO k=1,N-1         ! order switched
        DO i=IstrU,Iend

! tl line 6:         tl_ru(i,j,k,nrhs)=-0.5d0*dy_u(i,j)*
! tl line 6:     *     ( (tl_Hz(i,j,k)+tl_Hz(i-1,j,k))*phix(i)+
! tl line 6:     *       (Hz(i,j,k)+Hz(i-1,j,k))*tl_phix(i) )

         tmp=0.5d0*dy_u(i,j)*ad_ru(i,j,k,nrhs)
         tmp1=tmp*phix(i,k)

         ad_Hz(i,j,k)=ad_Hz(i,j,k)-tmp1
         ad_Hz(i-1,j,k)=ad_Hz(i-1,j,k)-tmp1
         ad_phix(i)=ad_phix(i)-tmp*(Hz(i,j,k)+Hz(i-1,j,k))
         ad_ru(i,j,k,nrhs)=0.d0

! tl line 5:         tl_phix(i)=tl_phix(i)+fac3*( tl_cff1*cff3+cff1*tl_cff3
! tl line 5:     *                               -tl_cff2*cff4-cff2*tl_cff4)
         tmp=fac3*ad_phix(i)
         ad_cff1=ad_cff1+tmp*cff3(i,k)
         ad_cff3=ad_cff3+tmp*cff1(i,k)
         ad_cff2=ad_cff2-tmp*cff4(i,k)
         ad_cff4=ad_cff4-tmp*cff2(i,k)

! tl line 4:         tl_cff4=tl_z_r(i,j,k+1)-tl_z_r(i-1,j,k+1)+
! tl line 4:     *           tl_z_r(i,j,k  )-tl_z_r(i-1,j,k  )
         ad_z_r(i,j,k+1)=ad_z_r(i,j,k+1)    +ad_cff4
         ad_z_r(i-1,j,k+1)=ad_z_r(i-1,j,k+1)-ad_cff4
         ad_z_r(i,j,k  )=ad_z_r(i,j,k  )    +ad_cff4
         ad_z_r(i-1,j,k  )=ad_z_r(i-1,j,k  )-ad_cff4
         ad_cff4=0.d0

! tl line 3:         tl_cff3=tl_z_r(i,j,k+1)+tl_z_r(i-1,j,k+1)-
! tl line 3:     *           tl_z_r(i,j,k  )-tl_z_r(i-1,j,k  )
         ad_z_r(i,j,k+1)=ad_z_r(i,j,k+1)    +ad_cff3
         ad_z_r(i-1,j,k+1)=ad_z_r(i-1,j,k+1)+ad_cff3
         ad_z_r(i,j,k  )=ad_z_r(i,j,k  )    -ad_cff3
         ad_z_r(i-1,j,k  )=ad_z_r(i-1,j,k  )-ad_cff3
         ad_cff3=0.d0

! tl line 2:         tl_cff2=tl_rho(i,j,k+1)+tl_rho(i-1,j,k+1)-
! tl line 2:     *           tl_rho(i,j,k  )-tl_rho(i-1,j,k  )
         ad_rho(i,j,k+1)=ad_rho(i,j,k+1)    +ad_cff2
         ad_rho(i-1,j,k+1)=ad_rho(i-1,j,k+1)+ad_cff2
         ad_rho(i,j,k  )=ad_rho(i,j,k  )    -ad_cff2
         ad_rho(i-1,j,k  )=ad_rho(i-1,j,k  )-ad_cff2
         ad_cff2=0.d0

! tl line 1:         tl_cff1=tl_rho(i,j,k+1)-tl_rho(i-1,j,k+1)+
! tl line 1:     *           tl_rho(i,j,k  )-tl_rho(i-1,j,k  )
         ad_rho(i,j,k+1)=ad_rho(i,j,k+1)    +ad_cff1
         ad_rho(i-1,j,k+1)=ad_rho(i-1,j,k+1)-ad_cff1
         ad_rho(i,j,k  )=ad_rho(i,j,k  )    +ad_cff1
         ad_rho(i-1,j,k  )=ad_rho(i-1,j,k  )-ad_cff1
         ad_cff1=0.d0

        END DO
       END DO


!
!  ADJ TO: Compute surface baroclinic pressure gradient (XI)
!
       DO i=IstrU,Iend
! tl line 4:        tl_ru(i,j,N,nrhs)=-0.5d0*dy_u(i,j)*
! tl line 4:     *   ( (tl_Hz(i,j,N)+tl_Hz(i-1,j,N))*phix(i)+
! tl line 4:     *     (Hz(i,j,N)+Hz(i-1,j,N))*tl_phix(i) )
 
        tmp=0.5d0*dy_u(i,j)*ad_ru(i,j,N,nrhs)
        tmp1=tmp*phix(i,N)

        ad_Hz(i,j,N)=ad_Hz(i,j,N)-tmp1
        ad_Hz(i-1,j,N)=ad_Hz(i-1,j,N)-tmp1
        ad_phix(i)=ad_phix(i)-tmp*(Hz(i,j,N)+Hz(i-1,j,N))

        ad_ru(i,j,N,nrhs)=0.d0

! tl line 3:        tl_phix(i)=tl_phix(i)+
! tl line 3:     *            fac1*(tl_rho(i,j,N)+tl_rho(i-1,j,N))*
! tl line 3:     *            (z_w(i,j,N)-z_w(i-1,j,N))
! tl line 3:     *           +(fac2+fac1*(rho(i,j,N)+rho(i-1,j,N)))*
! tl line 3:     *            (tl_z_w(i,j,N)-tl_z_w(i-1,j,N))

        tmp1=fac1*(z_w(i,j,N)-z_w(i-1,j,N))*ad_phix(i)
        tmp2=(fac2+fac1*(rho(i,j,N)+rho(i-1,j,N)))*ad_phix(i)

        ad_rho(i,j,N)=ad_rho(i,j,N)    +tmp1
        ad_rho(i-1,j,N)=ad_rho(i-1,j,N)+tmp1
        ad_z_w(i,j,N)=ad_z_w(i,j,N)    +tmp2
        ad_z_w(i-1,j,N)=ad_z_w(i-1,j,N)-tmp2

! tl line 2:        tl_phix(i)=fac1*( (tl_rho(i,j,N)-tl_rho(i-1,j,N))*cff1+
! tl line 2:     *                    (rho(i,j,N)-rho(i-1,j,N))*tl_cff1 )
        tmp=ad_phix(i)*fac1
        tmp1=tmp*cff1(i,N)

        ad_rho(i,j,N)=ad_rho(i,j,N)+tmp1
        ad_rho(i-1,j,N)=ad_rho(i-1,j,N)-tmp1
        ad_cff1=ad_cff1+(rho(i,j,N)-rho(i-1,j,N))*tmp
        ad_phix(i)=0.d0

! tl line 1:        tl_cff1=tl_z_w(i  ,j,N)-tl_z_r(i  ,j,N)+ 
! tl line 1:     *          tl_z_w(i-1,j,N)-tl_z_r(i-1,j,N)
        ad_z_w(i  ,j,N)=ad_z_w(i  ,j,N)+ad_cff1
        ad_z_r(i  ,j,N)=ad_z_r(i  ,j,N)-ad_cff1
        ad_z_w(i-1,j,N)=ad_z_w(i-1,j,N)+ad_cff1
        ad_z_r(i-1,j,N)=ad_z_r(i-1,j,N)-ad_cff1
        ad_cff1=0.d0
       END DO

      END DO

      RETURN
      END

#  endif /* DJ_STD2 */

# endif /* not PRSGRD_BT_ONLY */

# ifdef PRSGRD_BT_ONLY
      SUBROUTINE ad_prsgrd

      USE mod_ocean
      implicit none
      integer :: i,j,k

      real(8) :: P(LBi:UBi,LBj:UBj,N)
      real(8) :: ad_P(LBi:UBi,LBj:UBj,N)
      real(8) :: GRho,Grho0,HalfGRho
      real(8) :: cff,cff1

      GRho=g/rho0
      GRho0=1000.0d0*GRho

      ad_P=0.d0

!-----------------------------------------------------------------------
! Background P:
!-----------------------------------------------------------------------
      DO j=JstrV-1,Jend

       DO i=IstrU-1,Iend
        P(i,j,N)=GRho0*z_w(i,j,N)                            
       END DO

       DO k=N-1,1,-1
        DO i=IstrU-1,Iend
         P(i,j,k)=P(i,j,k+1)     
        END DO
       END DO

      END DO  ! <- j=JstrV-1,Jend

!
!-----------------------------------------------------------------------
!  ADJ to ETA-component pressure gradient term.
!-----------------------------------------------------------------------
!
      DO k=1,N  ! Note: order flipped (not essential yet for b/t case)

       DO j=JstrV,Jend
        DO i=Istr,Iend
! NL:         rv(i,j,k,nrhs)=dx_v(i,j)*Hz_v(i,j,k)*
! NL:     *                       ( P(i,j-1,k)-P(i,j,k) ) 

! TL:         tl_rv(i,j,k,nrhs)=dx_v(i,j)*(
! TL:     *             tl_Hz_v(i,j,k)*( P(i,j-1,k)-P(i,j,k) ) +
! TL:     *             Hz_v(i,j,k)*( tl_P(i,j-1,k)-tl_P(i,j,k) ) 
! TL:     *                               )
         cff=ad_rv(i,j,k,nrhs)*dx_v(i,j)
         cff1=cff*Hz_v(i,j,k)
         ad_Hz_v(i,j,k)=ad_Hz_v(i,j,k)+cff*( P(i,j-1,k)-P(i,j,k) )
         ad_P(i,j-1,k)=ad_P(i,j-1,k)+cff1
         ad_P(i,j,k)=ad_P(i,j,k)-cff1
         ad_rv(i,j,k,nrhs)=0.d0
        END DO
       END DO

      END DO ! k=1,N ETA comp

!-----------------------------------------------------------------------
!  ADJ to XI-component pressure gradient term.
!-----------------------------------------------------------------------
!
      DO k=1,N ! Note: order flipped (not essential yet for b/t case)
!
       DO j=Jstr,Jend
        DO i=IstrU,Iend
! NL:         ru(i,j,k,nrhs)=dy_u(i,j)*Hz_u(i,j,k)*
! NL:     *                     (P(i-1,j,k)-P(i,j,k))  

! TL:         tl_ru(i,j,k,nrhs)=dy_u(i,j)*(
! TL:     *             tl_Hz_u(i,j,k)*(P(i-1,j,k)-P(i,j,k)) +  
! TL:     *             Hz_u(i,j,k)*( tl_P(i-1,j,k)-tl_P(i,j,k) )
! TL:     *                               )
         cff=ad_ru(i,j,k,nrhs)*dy_u(i,j)
         cff1=cff*Hz_u(i,j,k)
         ad_Hz_u(i,j,k)=ad_Hz_u(i,j,k)+cff*(P(i-1,j,k)-P(i,j,k))
         ad_P(i-1,j,k)=ad_P(i-1,j,k)+cff1
         ad_P(i,j,k)=ad_P(i,j,k)-cff1
         ad_ru(i,j,k,nrhs)=0.d0
        END DO
       END DO

      END DO ! <- "k=1,N"

!-----------------------------------------------------------------------
! Adj to tl_P computation:
!-----------------------------------------------------------------------
      DO j=JstrV-1,Jend

       DO k=1,N-1  ! Note: order flipped
        DO i=IstrU-1,Iend
! TL:         tl_P(i,j,k)=tl_P(i,j,k+1)
         ad_P(i,j,k+1)=ad_P(i,j,k+1)+ad_P(i,j,k)
         ad_P(i,j,k)=0.d0
        END DO
       END DO

       DO i=IstrU-1,Iend
! TL:        tl_P(i,j,N)=GRho0*tl_z_w(i,j,N)
        ad_z_w(i,j,N)=ad_z_w(i,j,N)+GRho0*ad_P(i,j,N)
        ad_P(i,j,N)=0.d0
       END DO

      END DO  ! <- j=JstrV-1,Jend
       
      RETURN
      END 

# endif /* PRSGRD_BT_ONLY */

#endif /* ADJOINT */
