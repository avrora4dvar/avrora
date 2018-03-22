      SUBROUTINE ad_prsgrd_SDJ(L,M,N,
     *                ad_rho,z_w,z_r,dx,dy,
     *                mask_u,mask_v,
     *                ad_ru,ad_rv)

      implicit none

! A. Kurapov (01/22/2008): compute ADJ to 
! TL baroclinic pressure gradient term
!
! -g/rho0*int_z^0 grad(tl_rho) dz
!
! (units m/s2) for covariance implementation.
! Assume that z_w, z_r do not vary
! Use algorithm ROMS/Nonlinear/prsgrd31.h
! (1) remove bt pressuge grad
! (2) divide tl_ru, tl_rv by the volume

! Open boundary configuration

      real(8), intent(in) :: z_w(0:L,0:M,0:N)
      real(8), intent(in) :: z_r(0:L,0:M,N)
      real(8), intent(in) :: dx(0:L,0:M)
      real(8), intent(in) :: dy(0:L,0:M)
      real(8), intent(in) :: mask_u(0:L,0:M)
      real(8), intent(in) :: mask_v(0:L,0:M)
       
      real(8) :: ad_rho(0:L,0:M,N)
      real(8) :: ad_ru(0:L,0:M,N)
      real(8) :: ad_rv(0:L,0:M,N)

      integer :: L,M,N,i,j,k
      integer :: Istr,Iend,IstrU,IendU,IstrV,IendV
      integer :: Jstr,Jend,JstrU,JendU,JstrV,JendV
      real(8) :: fac1,fac3,ad_cff1,ad_cff2,cff1,cff3,cff4,g,rho0

      real(8) :: odxu(0:L),odyv(0:L),ad_phix(0:L),ad_phie(0:L),tmp

      g=9.81d0
      rho0=1025.d0

      fac1=0.5d0*g/rho0
      fac3=0.25d0*g/rho0

      ad_cff1=0.d0
      ad_cff2=0.d0
      ad_phix=0.d0
      ad_phie=0.d0

!====================================
! Array bounds (interior points):
!====================================

! - rho:
      Istr=1
      Iend=L-1
      Jstr=1
      Jend=M-1

! - u:
      IstrU=2
      IendU=L-1
      JstrU=1
      JendU=M-1

! - v:
      IstrV=1
      IendV=L-1
      JstrV=2
      JendV=M-1

!
!-----------------------------------------------------------------------
!  ADJ: Calculate pressure gradient in the ETA-direction (m/s2).
!-----------------------------------------------------------------------

      DO j=JstrV,JendV

       DO i=IstrV,IendV
        odyv(i)=2.d0/(dy(i,j)+dy(i,j-1))*mask_v(i,j)
       END DO

!
!  Interior layers:
!
       DO k=1,N-1
        DO i=IstrV,IendV

         cff3=z_r(i,j,k+1)+z_r(i,j-1,k+1)- 
     *        z_r(i,j,k  )-z_r(i,j-1,k  )
         cff4=z_r(i,j,k+1)-z_r(i,j-1,k+1)+
     *        z_r(i,j,k  )-z_r(i,j-1,k  )

! tl 4:         tl_rv(i,j,k)=-tl_phie(i)*odyv(i)
         ad_phie(i)=ad_phie(i)-ad_rv(i,j,k)*odyv(i)
         ad_rv(i,j,k)=0.d0

! tl 3:         tl_phie(i)=tl_phie(i)+fac3*(tl_cff1*cff3-tl_cff2*cff4)
! (recurs)
         ad_cff1=ad_cff1+fac3*cff3*ad_phie(i)
         ad_cff2=ad_cff2-fac3*cff4*ad_phie(i)
  

! tl 2:         tl_cff2=tl_rho(i,j,k+1)+tl_rho(i,j-1,k+1)-
!     *        tl_rho(i,j,k  )-tl_rho(i,j-1,k  )
         ad_rho(i,j  ,k+1)=ad_rho(i,j  ,k+1)+ad_cff2
         ad_rho(i,j-1,k+1)=ad_rho(i,j-1,k+1)+ad_cff2
         ad_rho(i,j  ,k  )=ad_rho(i,j  ,k  )-ad_cff2
         ad_rho(i,j-1,k  )=ad_rho(i,j-1,k  )-ad_cff2
         ad_cff2=0.d0

! tl 1:         tl_cff1=tl_rho(i,j,k+1)-tl_rho(i,j-1,k+1)+
!     *        tl_rho(i,j,k  )-tl_rho(i,j-1,k  )
         ad_rho(i,j  ,k+1)=ad_rho(i,j  ,k+1)+ad_cff1
         ad_rho(i,j-1,k+1)=ad_rho(i,j-1,k+1)-ad_cff1
         ad_rho(i,j  ,k  )=ad_rho(i,j  ,k  )+ad_cff1
         ad_rho(i,j-1,k  )=ad_rho(i,j-1,k  )-ad_cff1
         ad_cff1=0.d0

        END DO
       END DO

! Surface layer:
       DO i=IstrV,IendV
        cff1=z_w(i,j  ,N)-z_r(i,j  ,N)+
     *       z_w(i,j-1,N)-z_r(i,j-1,N)

! tl 2:        tl_rv(i,j,N)=-tl_phie(i)*odyv(i)
        ad_phie(i)=ad_phie(i)-ad_rv(i,j,N)*odyv(i)
        ad_rv(i,j,N)=0.d0

! tl 1:        tl_phie(i)=fac1*(tl_rho(i,j,N)-tl_rho(i,j-1,N))*cff1
        tmp=ad_phie(i)*fac1*cff1
        ad_rho(i,j  ,N)=ad_rho(i,j  ,N)+tmp
        ad_rho(i,j-1,N)=ad_rho(i,j-1,N)-tmp
        ad_phie(i)=0.d0
       END DO

      END DO ! "DO j=JstrV,JendV"

!-----------------------------------------------------------------------
!  Calculate pressure gradient in the XI-direction (m/s2)
!-----------------------------------------------------------------------
!

      DO j=JstrU,JendU

       DO i=IstrU,IendU
         odxu(i)=2.d0/(dx(i-1,j)+dx(i,j))*mask_u(i,j)
        END DO

!
!  Interior layers.
!
       DO k=1,N-1
        DO i=IstrU,Iend

         cff3=z_r(i,j,k+1)+z_r(i-1,j,k+1)-
     *        z_r(i,j,k  )-z_r(i-1,j,k  )
         cff4=z_r(i,j,k+1)-z_r(i-1,j,k+1)+
     *        z_r(i,j,k  )-z_r(i-1,j,k  )

! tl 4:          tl_ru(i,j,k)=-tl_phix(i)*odxu(i)
         ad_phix(i)=ad_phix(i)-ad_ru(i,j,k)*odxu(i)
         ad_ru(i,j,k)=0.d0

! tl 3:          tl_phix(i)=tl_phix(i)+fac3*(tl_cff1*cff3-tl_cff2*cff4)
! (recurs)
         ad_cff1=ad_cff1+fac3*cff3*ad_phix(i)
         ad_cff2=ad_cff2-fac3*cff4*ad_phix(i)

! tl 2:          tl_cff2=tl_rho(i,j,k+1)+tl_rho(i-1,j,k+1)-
!     *         tl_rho(i,j,k  )-tl_rho(i-1,j,k  )
         ad_rho(i  ,j,k+1)=ad_rho(i  ,j,k+1)+ad_cff2
         ad_rho(i-1,j,k+1)=ad_rho(i-1,j,k+1)+ad_cff2
         ad_rho(i  ,j,k  )=ad_rho(i  ,j,k  )-ad_cff2
         ad_rho(i-1,j,k  )=ad_rho(i-1,j,k  )-ad_cff2
         ad_cff2=0.d0

! tl 1:          tl_cff1=tl_rho(i,j,k+1)-tl_rho(i-1,j,k+1)+
!     *         tl_rho(i,j,k  )-tl_rho(i-1,j,k  )
         ad_rho(i  ,j,k+1)=ad_rho(i  ,j,k+1)+ad_cff1
         ad_rho(i-1,j,k+1)=ad_rho(i-1,j,k+1)-ad_cff1
         ad_rho(i  ,j,k  )=ad_rho(i  ,j,k  )+ad_cff1
         ad_rho(i-1,j,k  )=ad_rho(i-1,j,k  )-ad_cff1
         ad_cff1=0.d0

        END DO
       END DO

! Surface layer:
       DO i=IstrU,IendU
         cff1=z_w(i  ,j,N)-z_r(i  ,j,N)+
     *        z_w(i-1,j,N)-z_r(i-1,j,N)

! tl 2:         tl_ru(i,j,N)=-tl_phix(i)*odxu(i)
         ad_phix(i)=ad_phix(i)-odxu(i)*ad_ru(i,j,N)
         ad_ru(i,j,N)=0.d0

! tl 1:        tl_phix(i)=fac1*(tl_rho(i,j,N)-tl_rho(i-1,j,N))*cff1
         tmp=fac1*cff1*ad_phix(i)
         ad_rho(i  ,j,N)=ad_rho(i  ,j,N)+tmp
         ad_rho(i-1,j,N)=ad_rho(i-1,j,N)-tmp
         ad_phix(i)=0.d0

       END DO

      END DO ! j=JstrU,JendU


      RETURN
      END
