#include "cpp_balance_3D.h"
! A. Kurapov, 11/25/2008, Balance operator, 
! from
! {psi, T}
! to
! {zeta, ubar, vbar, u, v, T, S}

      PROGRAM tl_balance_3D

      implicit none
      include 'netcdf.h'
      character(len=80) :: univarfile,statefile,fwdfile
      character(len=80) :: grdfile

      integer :: xi_rho,eta_rho,xi_u,eta_u,xi_v,eta_v,xi_psi,eta_psi
      integer :: L,M,i,j,k,N,rec,index1d,ii,jj
      integer :: status,fid,dimid,varid,dimids(4),ncid
      integer :: start(4),count(4)
      integer :: flag_salt

      real(8) :: g,cff,cff1,tl_cff,f1,R0,Tcoef,Scoef,hij

      real(8) :: time

      real(8), allocatable :: f(:,:),h(:,:),dx(:,:),dy(:,:)
      real(8), allocatable :: z_r(:,:,:),z_w(:,:,:),Hz(:,:,:)
      real(8), allocatable :: mask_r(:,:),mask_u(:,:),mask_v(:,:)

      real(8), allocatable :: tl_psi(:,:)
      real(8), allocatable :: tl_zeta(:,:)
      real(8), allocatable :: tl_ubar(:,:)
      real(8), allocatable :: tl_vbar(:,:)
      real(8), allocatable :: tl_u(:,:,:)
      real(8), allocatable :: tl_v(:,:,:)
      real(8), allocatable :: tl_T(:,:,:)
      real(8), allocatable :: tl_S(:,:,:)
      real(8), allocatable :: tl_rho(:,:,:)
      real(8), allocatable :: T(:,:,:)
      real(8), allocatable :: S(:,:,:)
      real(8), allocatable :: tl_Fpsi(:,:)

! b/c prs grad terms:
      real(8), allocatable :: tl_ru(:,:,:)
      real(8), allocatable :: tl_rv(:,:,:)

! fluxes in the rhs of the Poisson eq.
      real(8), allocatable :: tl_Fx(:,:)
      real(8), allocatable :: tl_Fe(:,:)
      real(8), allocatable :: tl_F(:,:)     ! rhs of the Poisson eq. (flux div)
      real(8), allocatable :: cffx(:,:)
      real(8), allocatable :: cffe(:,:)

! definitions for the gauss solver
      integer :: p,q,Nb
      real(8), allocatable :: A(:,:),tl_b(:) 

#ifdef PG_BOTTOM_TREATMENT
      real(8) :: dlz,dlz2,dist,dist0,z_bot,z_u,z_v
#endif
      g=9.81

      start=0
      count=0

!================================================
! READ IN-FILE
!================================================
      READ (*,'(a1)') 
      READ (*,'(a80)') grdfile
      READ (*,'(a1)') 
      READ (*,'(a80)') univarfile
      READ (*,'(a1)') 
      READ (*,'(a80)') statefile
      READ (*,'(a1)') 
      READ (*,'(a80)') fwdfile
      READ (*,'(a1)') 
      READ (*,*) R0,Tcoef,Scoef

      write (*,*) 'Grid file: ',grdfile
      write (*,*) 'Univar file: ',univarfile
      write (*,*) 'Output file: ',statefile
      write (*,*) 'FWD file: ',fwdfile
      write (*,*) 'R0,Tcoef,Scoef: ',R0,Tcoef,Scoef


!================================================
! READ DIMENSIONS
!================================================
      CALL get_grddim(grdfile,
     *                      xi_rho,eta_rho,
     *                      xi_u,eta_u,
     *                      xi_v,eta_v,
     *                      xi_psi,eta_psi,
     *                      N)


      write (*,*) 'xi_rho,eta_rho=',xi_rho,eta_rho
      write (*,*) 'xi_u  ,eta_u  =',xi_u,eta_u
      write (*,*) 'xi_v  ,eta_v  =',xi_v,eta_v
      write (*,*) 'xi_psi,eta_psi=',xi_psi,eta_psi
      write (*,*) 'N=',N

      L=xi_rho-1
      M=eta_rho-1

# define RHO_RANGE 0:L,0:M
# define U_RANGE 1:L,0:M
# define V_RANGE 0:L,1:M
# define I_INTERIOR_RHO 1,L-1
# define J_INTERIOR_RHO 1,M-1
# define I_INTERIOR_U 2,L-1
# define J_INTERIOR_U 1,M-1
# define I_INTERIOR_V 1,L-1
# define J_INTERIOR_V 2,M-1

!================================================
! ALLOCATE VARIABLES
!================================================

! - grid
      allocate(dx(RHO_RANGE))
      allocate(dy(RHO_RANGE))
      allocate(f(RHO_RANGE))
      allocate(h(RHO_RANGE))
      allocate(z_r(RHO_RANGE,N))
      allocate(z_w(RHO_RANGE,0:N))
      allocate(Hz(RHO_RANGE,N))
      allocate(mask_r(RHO_RANGE))
      allocate(mask_u(RHO_RANGE))
      allocate(mask_v(RHO_RANGE))

      dx=0.d0
      dy=0.d0
      f=0.d0
      h=0.d0
      z_r=0.d0
      z_w=0.d0
      Hz=0.d0
      mask_r=0
      mask_u=0
      mask_v=0


! - tl and background variables:    
      allocate(tl_psi(RHO_RANGE))
      allocate(tl_T(RHO_RANGE,N))
      allocate(tl_S(RHO_RANGE,N))
      allocate(tl_zeta(RHO_RANGE))
      allocate(tl_ubar(RHO_RANGE))
      allocate(tl_vbar(RHO_RANGE))
      allocate(tl_u(RHO_RANGE,N))
      allocate(tl_v(RHO_RANGE,N))
      allocate(tl_rho(RHO_RANGE,N))

      allocate(tl_ru(RHO_RANGE,N))
      allocate(tl_rv(RHO_RANGE,N))

      allocate(T(RHO_RANGE,N))
      allocate(S(RHO_RANGE,N))
      T=0.d0
      S=0.d0

      tl_psi=0.d0
      tl_T=0.d0
      tl_S=0.d0
      tl_zeta=0.d0
      tl_ubar=0.d0
      tl_vbar=0.d0
      tl_u=0.d0
      tl_v=0.d0

      tl_ru=0.d0
      tl_rv=0.d0

! - auxilary arrays:
      allocate(cffx(RHO_RANGE))
      allocate(cffe(RHO_RANGE))
      allocate(tl_Fx(RHO_RANGE))
      allocate(tl_Fe(RHO_RANGE))
      allocate(tl_F(RHO_RANGE))
      cffx=0.d0
      cffe=0.d0
      tl_Fx=0.d0
      tl_Fe=0.d0
      tl_F=0.d0

!================================================
! READ GRID
!================================================
      rec=1
      CALL get_2dfld(grdfile,'h' ,rec,xi_rho,eta_rho,h(RHO_RANGE)) 
      CALL get_2dfld(grdfile,'f' ,rec,xi_rho,eta_rho,f(RHO_RANGE))
      CALL get_2dfld(grdfile,'dx',rec,xi_rho,eta_rho,dx(RHO_RANGE))
      CALL get_2dfld(grdfile,'dy',rec,xi_rho,eta_rho,dy(RHO_RANGE))
      write (*,*) 'read mask_r:'
      CALL get_2dfld
     * (grdfile,'mask_rho',rec,xi_rho,eta_rho,mask_r(RHO_RANGE))
      write (*,*) 'read mask_u:'
       CALL get_2dfld
     * (grdfile,'mask_u',rec,xi_u,eta_u,mask_u(U_RANGE))
      write (*,*) 'read mask_v:'
      CALL get_2dfld
     * (grdfile,'mask_v',rec,xi_v,eta_v,mask_v(V_RANGE))

      write (*,*) 'read z_r:'
      CALL get_3dfld
     * (grdfile,'z0_r',rec,xi_rho,eta_rho,N,z_r(RHO_RANGE,1:N))  
      write (*,*) 'read z_w:'
      CALL get_3dfld
     * (grdfile,'z0_w',rec,xi_rho,eta_rho,N+1,z_w(RHO_RANGE,0:N))
  
      write (*,*) 'read grid, ... done'
    
!================================================
! SET LAYER DEPTH:
!================================================
      DO j=0,M
       DO i=0,L
        DO k=1,N         
         Hz(i,j,k)=z_w(i,j,k)-z_w(i,j,k-1)
        END DO
       END DO
      END DO
      
#ifndef CONST_DSDT
!================================================
! READ BACKGROUND T-S (first record in fwdfile)
!================================================
      rec=1
      write (*,*) 'read fwd temp:'
      CALL get_3dfld
     * (fwdfile,'temp',rec,xi_rho,eta_rho,N,T(RHO_RANGE,1:N))  
      write (*,*) 'read fwd salt:'
      CALL get_3dfld
     * (fwdfile,'salt',rec,xi_rho,eta_rho,N,S(RHO_RANGE,1:N))  
#endif


!================================================
! READ TL_PSI, TL_T
!================================================
      rec=1
      write (*,*) 'read tl psi:'
      CALL get_2dfld
     * (univarfile,'psi',rec,xi_rho,eta_rho,tl_psi(RHO_RANGE))
      write (*,*) 'read tl temp:'
      CALL get_3dfld
     * (univarfile,'temp',rec,xi_rho,eta_rho,N,tl_T(RHO_RANGE,1:N))
      write (*,*) 'read tl temp .. done'

      status=nf_open(trim(univarfile),0,ncid)
      flag_salt=nf_inq_varid(ncid,'salt',varid)
      status=nf_close(ncid)
      IF(flag_salt.EQ.0) THEN
       write (*,*) 'read tl salt:'
       CALL get_3dfld
     * (univarfile,'salt',rec,xi_rho,eta_rho,N,tl_S(RHO_RANGE,1:N))
       write(*,*) 'read tl salt .. done'
      END IF


!================================================
! TL_PSI=0 in exterior and masked locations
!================================================

      write (*,*) 'tl_psi in exterior and masked locations ...'

      DO i=0,L
       tl_psi(i,0)=0.d0
       tl_psi(i,M)=0.d0
      END DO

      DO j=0,M
       tl_psi(0,j)=0.d0
       tl_psi(L,j)=0.d0
      END DO

      DO j=J_INTERIOR_RHO
       DO i=I_INTERIOR_RHO
         tl_psi(i,j)=tl_psi(i,j)*mask_r(i,j)
       END DO
      END DO

!================================================
! tl_S 
!================================================

      IF(flag_salt.NE.0) THEN

#ifdef CONST_DSDT 
      cff=-0.16d0
#endif
      DO k=1,N
       DO j=J_INTERIOR_RHO
        DO i=I_INTERIOR_RHO
         IF (abs(mask_r(i,j)).gt.1.d-10) THEN

          ! cff=dS/dT
                   
#ifndef CONST_DSDT
          IF (k.eq.1) THEN
           cff=( S(i,j,k+1)-S(i,j,k) )/
     *         ( T(i,j,k+1)-T(i,j,k) )    
          ELSEIF (k.eq.N) THEN
           cff=( S(i,j,k)-S(i,j,k-1) )/
     *         ( T(i,j,k)-T(i,j,k-1) )    
          ELSE
           cff=( S(i,j,k+1)-S(i,j,k-1) )/
     *         ( T(i,j,k+1)-T(i,j,k-1) )    
          END IF                      
#endif

          tl_S(i,j,k)=cff*tl_T(i,j,k)

         END IF ! mask_r(i,j).ne.0
        END DO
       END DO
      END DO

      END IF

!================================================
! tl_RHO 
!================================================
      write (*,*) 'tl_RHO'

      DO k=1,N
       DO j=J_INTERIOR_RHO
        DO i=I_INTERIOR_RHO

         tl_rho(i,j,k)=R0*(-Tcoef*tl_T(i,j,k)+Scoef*tl_S(i,j,k))

        END DO
       END DO
      END DO

!================================================
! b/c pressure grad term: tl_ru, tl_rv = -(g/rho0) \int_z^0 grad rho |_z dz'
!================================================

      CALL tl_prsgrd_SDJ(L,M,N,
     *                tl_rho,z_w,z_r,dx,dy,
     *                mask_u,mask_v,
     *                tl_ru,tl_rv)


#ifdef PG_BOTTOM_TREATMENT
!================================================
! reduce b/c pressure grad term near bottom
!================================================
      WRITE(*,*) 'Applying bottom treatment'
 
      dist0=20.d0
      dlz=5.d0
      dlz2=dlz*dlz

      DO j=0,M
       DO i=1,L

        z_bot=-0.5d0*(h(i,j)+h(i-1,j))

        DO k=1,N

         z_u=0.5d0*(z_r(i,j,k)+z_r(i-1,j,k))
         dist=z_u-z_bot-dist0

         IF (dist.lt.0.d0) THEN   
          cff=DEXP(-dist*dist/(2.d0*dlz2))
          tl_ru(i,j,k)=tl_ru(i,j,k)*cff
         END IF

        END DO

       END DO
      END DO

      DO j=1,M
       DO i=0,L

        z_bot=-0.5d0*(h(i,j)+h(i,j-1))

        DO k=1,N

         z_v=0.5d0*(z_r(i,j,k)+z_r(i,j-1,k))
         dist=z_v-z_bot-dist0

         IF (dist.lt.0.d0) THEN   
          cff=DEXP(-dist*dist/(2.d0*dlz2))
          tl_rv(i,j,k)=tl_rv(i,j,k)*cff
         END IF

        END DO

       END DO
      END DO

#endif


!================================================
! tl_F=div {f grad(tl_psi)+int_{-H}^0 (tl_ru,tl_rv)}
! intergate div term over the rho grid area, compute boundary fluxes
! shown above in {}. 
! BC: tl_psi=0 , b/c Fx, Fy = 0
!================================================

      write (*,*) 'RHS of Poisson eqn.'

! == Fluxes in xi-dir:
      DO j=0,M
       DO i=1,L

! tl 1:
! b/c (tl_Fx has initially been set to 0)
        DO k=1,N
         cff=0.5d0*(Hz(i-1,j,k)+Hz(i,j,k))
         tl_Fx(i,j)=tl_Fx(i,j)+cff*tl_ru(i,j,k)
        END DO

! tl 2:
! BC (b/c):
        tl_Fx(i,j)=tl_Fx(i,j)*mask_u(i,j)

! tl 3:
! add b/t (psi) part:
        cff=0.5d0*(f(i-1,j)+f(i,j))*(dy(i-1,j)+dy(i,j))
     *             /(dx(i-1,j)+dx(i,j))
        tl_Fx(i,j)=tl_Fx(i,j)+cff*(tl_psi(i,j)-tl_psi(i-1,j))                     
       END DO
      END DO

! == Fluxes in eta-dir:
      DO j=1,M
       DO i=0,L

! tl 1:
! b/c (tl_Fe has initially been set to 0)
        DO k=1,N
         cff=0.5d0*(Hz(i,j-1,k)+Hz(i,j,k))
         tl_Fe(i,j)=tl_Fe(i,j)+cff*tl_rv(i,j,k)
        END DO

! tl 2:
! BC (b/c):
        tl_Fe(i,j)=tl_Fe(i,j)*mask_v(i,j)

! tl 3:
! add b/t (psi) part:
        cff=0.5d0*(f(i,j-1)+f(i,j))*(dx(i,j-1)+dx(i,j))
     *             /(dy(i,j-1)+dy(i,j))
        tl_Fe(i,j)=tl_Fe(i,j)+cff*(tl_psi(i,j)-tl_psi(i,j-1))
       END DO
      END DO

! == rhs of the Poisson eq.

      DO j=J_INTERIOR_RHO
       DO i=I_INTERIOR_RHO
        tl_F(i,j)=tl_Fx(i+1,j)-tl_Fx(i,j)+tl_Fe(i,j+1)-tl_Fe(i,j)
        tl_F(i,j)=tl_F(i,j)*mask_r(i,j)
       END DO
      END DO


!================================================
! tl_zeta: Solve Poisson equation
!================================================
      write (*,*) 'Solve Poisson eqn.'

      q=L+1
      p=L+1
      Nb=(L+1)*(M+1)
      allocate(A(p+q+1,Nb),tl_b(Nb))

      tl_b=0.d0
      A=0.d0

! == H/dx at u-locations, H/dy at v-locations
      write (*,*) 'H/dx at u-locations, H/dy at v-locations'
      DO j=0,M
       DO i=1,L
        cffx(i,j)=(h(i-1,j)+h(i,j))/(dx(i-1,j)+dx(i,j))
       END DO
      END DO

      DO j=1,M
       DO i=0,L
        cffe(i,j)=(h(i,j-1)+h(i,j))/(dy(i,j-1)+dy(i,j))
       END DO
      END DO

! == Matrix coefficients:
      write (*,*) 'Matrix coeff...'
      DO j=0,M
       DO i=0,L

        ! ii = full matrix row number       
        ! jj = full matrix column number
        ii=index1d(i,j,L)

        tl_b(ii)=tl_F(i,j)

        IF (abs(mask_r(i,j))<1.e-8 .or. 
     *          i.eq.0 .or. 
     *          i.eq.L .or.
     *          j.eq.0 .or.
     *          j.eq.M ) THEN
         ! ( masked location, or edge: identity eqn. zeta=0 )
         A(q+1,ii)=1.d0         
        ELSE

         ! coefficient of zeta(i,j):
         A(q+1,ii)=-g*(cffx(i+1,j)+cffx(i,j)+cffe(i,j+1)+cffe(i,j))

         ! coefficient of zeta(i+1,j):
         jj=index1d(i+1,j,L)
         A(ii-jj+q+1,jj)=g*cffx(i+1,j)

         ! coefficient of zeta(i-1,j):
         jj=index1d(i-1,j,L)
         A(ii-jj+q+1,jj)=g*cffx(i,j)

         ! coefficient of zeta(i,j+1):
         jj=index1d(i,j+1,L)
         A(ii-jj+q+1,jj)=g*cffe(i,j+1)

         ! coefficient of zeta(i,j-1):
         jj=index1d(i,j-1,L)
         A(ii-jj+q+1,jj)=g*cffe(i,j)

        END IF

       END DO
      END DO

! == Gauss solver:
      CALL gauss_band(A,tl_b,Nb,p,q)

      DO j=0,M
       DO i=0,L
        ii=index1d(i,j,L)
        tl_zeta(i,j)=tl_b(ii)
       END DO
      END DO


!================================================
! tl_u, tl_v
!================================================

      write (*,*) 'tl_u, tl_v ...'


! == add BT pres grad term to bc pres grad, xi dir:
      DO j=0,M
       DO i=1,L
        cff1=-g*2.d0/(dx(i-1,j)+dx(i,j))
        tl_cff=cff1*(tl_zeta(i,j)-tl_zeta(i-1,j))

        DO k=1,N
         tl_ru(i,j,k)=tl_ru(i,j,k)+tl_cff
        END DO

       END DO
      END DO

! == add BT pres grad term, eta dir:
      DO j=1,M
       DO i=0,L
        cff1=-g*2.d0/(dy(i,j-1)+dy(i,j))
        tl_cff=cff1*(tl_zeta(i,j)-tl_zeta(i,j-1))

        DO k=1,N
         tl_rv(i,j,k)=tl_rv(i,j,k)+tl_cff
        END DO

       END DO
      END DO

! == tl_u

      DO j=J_INTERIOR_U
       DO i=I_INTERIOR_U
               
        cff1=2.d0*mask_u(i,j)/(f(i-1,j)+f(i,j))

        DO k=1,N
         tl_u(i,j,k)=0.25d0*cff1*(tl_rv(i-1,j+1,k)+tl_rv(i,j+1,k)+
     *                            tl_rv(i-1,j  ,k)+tl_rv(i,j  ,k))
        END DO

       END DO
      END DO


! == tl_v
      DO j=J_INTERIOR_V
       DO i=I_INTERIOR_V
               
        cff1=-2.d0*mask_v(i,j)/(f(i,j-1)+f(i,j))

        DO k=1,N
         tl_v(i,j,k)=0.25d0*cff1*(tl_ru(i  ,j  ,k)+tl_ru(i+1,j  ,k)+
     *                            tl_ru(i  ,j-1,k)+tl_ru(i+1,j-1,k))
        END DO

       END DO
      END DO

!================================================
! tl_ubar, tl_vbar (diagnostic)
!================================================

      DO j=J_INTERIOR_U
       DO i=I_INTERIOR_U

        DO k=1,N
         tl_ubar(i,j)=tl_ubar(i,j)+
     *                (Hz(i-1,j,k)+Hz(i,j,k))*tl_u(i,j,k)
        END DO

        tl_ubar(i,j)=tl_ubar(i,j)/(h(i-1,j)+h(i,j))
 
       END DO
      END DO
      
      DO j=J_INTERIOR_V
       DO i=I_INTERIOR_V

        DO k=1,N
         tl_vbar(i,j)=tl_vbar(i,j)+
     *                (Hz(i,j-1,k)+Hz(i,j,k))*tl_v(i,j,k)
        END DO

        tl_vbar(i,j)=tl_vbar(i,j)/(h(i,j-1)+h(i,j))
 
       END DO
      END DO

      
!================================================
! Output:
!================================================
      write (*,*) 'output tl fields...'


      time=0.d0

      CALL tl_output_ini(statefile,
     *                   L,M,N,
     *                   time,tl_zeta,tl_ubar,tl_vbar,
     *                   tl_u,tl_v,tl_T,tl_S)

!      CALL tl_output_ini('test.nc',
!     *                   L,M,N,
!     *                   time,tl_F,0.d0*tl_ubar,0.d0*tl_vbar,
!     *                   tl_ru,tl_rv,tl_rho,0.d0*tl_S)
      
      END

