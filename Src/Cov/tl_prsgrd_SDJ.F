      SUBROUTINE tl_prsgrd_SDJ(L,M,N,
     *                       tl_rho,z_w,z_r,dx,dy,
     *                       mask_u,mask_v,
     *                       tl_ru,tl_rv)

! A. Kurapov (01/05/2008): compute TL baroclinic pressure gradient term
!
! -g/rho0*int_z^0 grad(tl_rho) dz
!
! (units m/s2) for covariance implementation.
! Assume that z_w, z_r do not vary
! Use algorithm ROMS/Nonlinear/prsgrd31.h
! (1) remove bt pressuge grad
! (2) divide tl_ru, tl_rv by the volume

! Open boundary configuration

      implicit none

      real(8), intent(in) :: tl_rho(0:L,0:M,N)
      real(8), intent(in) :: z_w(0:L,0:M,0:N)
      real(8), intent(in) :: z_r(0:L,0:M,N)
      real(8), intent(in) :: dx(0:L,0:M)
      real(8), intent(in) :: dy(0:L,0:M)
      real(8), intent(in) :: mask_u(0:L,0:M)
      real(8), intent(in) :: mask_v(0:L,0:M)
       
      real(8), intent(out) :: tl_ru(0:L,0:M,N)
      real(8), intent(out) :: tl_rv(0:L,0:M,N)

      integer :: L,M,N,i,j,k
      integer :: Istr,Iend,IstrU,IendU,IstrV,IendV
      integer :: Jstr,Jend,JstrU,JendU,JstrV,JendV
      real(8) :: fac1,fac3,tl_cff1,tl_cff2,cff1,cff3,cff4,g,rho0

      real(8) :: odxu(0:L),odyv(0:L),tl_phix(0:L),tl_phie(0:L)

      g=9.81d0
      rho0=1025.d0

      fac1=0.5d0*g/rho0
      fac3=0.25d0*g/rho0

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

!-----------------------------------------------------------------------
!  Calculate pressure gradient in the XI-direction (m/s2)
!-----------------------------------------------------------------------
!
      DO j=JstrU,JendU

! Surface layer:
       DO i=IstrU,IendU
        cff1=z_w(i  ,j,N)-z_r(i  ,j,N)+
     *       z_w(i-1,j,N)-z_r(i-1,j,N)
        tl_phix(i)=fac1*(tl_rho(i,j,N)-tl_rho(i-1,j,N))*cff1

        odxu(i)=2.d0/(dx(i-1,j)+dx(i,j))*mask_u(i,j)

        tl_ru(i,j,N)=-tl_phix(i)*odxu(i)

       END DO
!
!  Interior layers.
!
       DO k=N-1,1,-1
        DO i=IstrU,Iend
         tl_cff1=tl_rho(i,j,k+1)-tl_rho(i-1,j,k+1)+
     *        tl_rho(i,j,k  )-tl_rho(i-1,j,k  )
         tl_cff2=tl_rho(i,j,k+1)+tl_rho(i-1,j,k+1)-
     *        tl_rho(i,j,k  )-tl_rho(i-1,j,k  )
         cff3=z_r(i,j,k+1)+z_r(i-1,j,k+1)-
     *        z_r(i,j,k  )-z_r(i-1,j,k  )
         cff4=z_r(i,j,k+1)-z_r(i-1,j,k+1)+
     *        z_r(i,j,k  )-z_r(i-1,j,k  )

         tl_phix(i)=tl_phix(i)+fac3*(tl_cff1*cff3-tl_cff2*cff4)

         tl_ru(i,j,k)=-tl_phix(i)*odxu(i)
        END DO
       END DO

      END DO ! j=JstrU,JendU

!
!-----------------------------------------------------------------------
!  Calculate pressure gradient in the ETA-direction (m/s2).
!-----------------------------------------------------------------------

      DO j=JstrV,JendV

! Surface layer:
       DO i=IstrV,IendV
        cff1=z_w(i,j  ,N)-z_r(i,j  ,N)+
     *       z_w(i,j-1,N)-z_r(i,j-1,N)
        tl_phie(i)=fac1*(tl_rho(i,j,N)-tl_rho(i,j-1,N))*cff1
        odyv(i)=2.d0/(dy(i,j)+dy(i,j-1))*mask_v(i,j)
        tl_rv(i,j,N)=-tl_phie(i)*odyv(i)
       END DO

!
!  Interior layers:
!
       DO k=N-1,1,-1
        DO i=IstrV,IendV
         tl_cff1=tl_rho(i,j,k+1)-tl_rho(i,j-1,k+1)+
     *        tl_rho(i,j,k  )-tl_rho(i,j-1,k  )
         tl_cff2=tl_rho(i,j,k+1)+tl_rho(i,j-1,k+1)-
     *        tl_rho(i,j,k  )-tl_rho(i,j-1,k  )
         cff3=z_r(i,j,k+1)+z_r(i,j-1,k+1)- 
     *        z_r(i,j,k  )-z_r(i,j-1,k  )
         cff4=z_r(i,j,k+1)-z_r(i,j-1,k+1)+
     *        z_r(i,j,k  )-z_r(i,j-1,k  )

         tl_phie(i)=tl_phie(i)+fac3*(tl_cff1*cff3-tl_cff2*cff4)

         tl_rv(i,j,k)=-tl_phie(i)*odyv(i)
        END DO
       END DO
      END DO

      RETURN
      END
