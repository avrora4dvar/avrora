#include "cpp_balance_3D.h"
! A. Kurapov, 1/16/2008, Adj to balance operator, 
! from
! {psi, T}
! to
! {zeta, ubar, vbar, u, v, T, S}

! Input: statefile
! Output: univarfile

      PROGRAM ad_balance_3D
  
      implicit none
      include 'netcdf.h'
      character(len=80) :: univarfile,statefile,fwdfile
      character(len=80) :: grdfile
      
      integer :: xi_rho,eta_rho,xi_u,eta_u,xi_v,eta_v,xi_psi,eta_psi
      integer :: L,M,i,j,k,N,rec,index1d,ii,jj
      integer :: status,fid,dimid,varid,dimids(4)
      integer :: start(4),count(4)

      real(8) :: g,cff,cff1,ad_cff,f1,R0,Tcoef,Scoef,hij
      real(8) :: tmp,tmp1,tmp2,tmp3,tmp4


      real(8) :: time

      real(8), allocatable :: f(:,:),h(:,:),dx(:,:),dy(:,:)
      real(8), allocatable :: z_r(:,:,:),z_w(:,:,:),Hz(:,:,:)
      real(8), allocatable :: mask_r(:,:),mask_u(:,:),mask_v(:,:)

      real(8), allocatable :: ad_psi(:,:)
      real(8), allocatable :: ad_zeta(:,:)
      real(8), allocatable :: ad_u(:,:,:)
      real(8), allocatable :: ad_v(:,:,:)
      real(8), allocatable :: ad_T(:,:,:)
      real(8), allocatable :: ad_S(:,:,:)
      real(8), allocatable :: ad_rho(:,:,:)
      real(8), allocatable :: T(:,:,:)
      real(8), allocatable :: S(:,:,:)
      real(8), allocatable :: ad_Fpsi(:,:)

! b/c prs grad terms:
      real(8), allocatable :: ad_ru(:,:,:)
      real(8), allocatable :: ad_rv(:,:,:)

! fluxes in the rhs of the Poisson eq.
      real(8), allocatable :: ad_Fx(:,:)
      real(8), allocatable :: ad_Fe(:,:)
      real(8), allocatable :: ad_F(:,:)     ! rhs of the Poisson eq. (flux div)
      real(8), allocatable :: cffx(:,:)
      real(8), allocatable :: cffe(:,:)

! definitions for the gauss solver
      integer :: p,q,Nb
      real(8), allocatable :: A(:,:),ad_b(:) 

#ifdef PG_BOTTOM_TREATMENT
      real(8) :: dlz,dlz2,dist,dist0,z_bot,z_u,z_v
#endif

      ad_cff=0.d0

      g=9.81
      time=0.d0

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
      write (*,*) 'Multivar file: ',statefile
      write (*,*) 'FWD file: ',fwdfile
      write (*,*) 'R0,Tcoef,Scoef: ',R0,Tcoef,Scoef


!================================================
! READ DIMENSIONS
!================================================
      WRITE(*,*) 'Read grid dimensions'
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
      allocate(ad_psi(RHO_RANGE))
      allocate(ad_T(RHO_RANGE,N))
      allocate(ad_S(RHO_RANGE,N))
      allocate(ad_zeta(RHO_RANGE))
      allocate(ad_u(RHO_RANGE,N))
      allocate(ad_v(RHO_RANGE,N))
      allocate(ad_rho(RHO_RANGE,N))

      allocate(ad_ru(RHO_RANGE,N))
      allocate(ad_rv(RHO_RANGE,N))

      allocate(T(RHO_RANGE,N))
      allocate(S(RHO_RANGE,N))
      T=0.d0
      S=0.d0

      ad_psi=0.d0
      ad_T=0.d0
      ad_S=0.d0
      ad_zeta=0.d0
      ad_u=0.d0
      ad_v=0.d0

      ad_ru=0.d0
      ad_rv=0.d0

! - auxilary arrays:
      allocate(cffx(RHO_RANGE))
      allocate(cffe(RHO_RANGE))
      allocate(ad_Fx(RHO_RANGE))
      allocate(ad_Fe(RHO_RANGE))
      allocate(ad_F(RHO_RANGE))
      cffx=0.d0
      cffe=0.d0
      ad_Fx=0.d0
      ad_Fe=0.d0
      ad_F=0.d0

!================================================
! READ GRID
!================================================
      WRITE(*,*) 'Reading grid'
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
  
      write(*,*) 'read grid, ... done'
    
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
! ADJ to Output, input of adjoint variables
!================================================

      write(*,*) 'Reading statefile'
      CALL ad_input_ini(statefile,
     *                   L,M,N,
     *                   time,ad_zeta,
     *                   ad_u,ad_v,ad_T,ad_S)

!================================================
! ADJ to: tl_u, tl_v
!================================================
      write(*,*) 'calculating'
! == adj to: tl_v

      DO j=J_INTERIOR_V
       DO i=I_INTERIOR_V
               
        cff1=-2.d0*mask_v(i,j)/(f(i,j-1)+f(i,j))

        DO k=1,N
! tl:         tl_v(i,j,k)=0.25d0*cff1*(tl_ru(i  ,j  ,k)+tl_ru(i+1,j  ,k)+
!     *                            tl_ru(i  ,j-1,k)+tl_ru(i+1,j-1,k))

         tmp=0.25d0*cff1*ad_v(i,j,k)

         ad_ru(i  ,j  ,k)=ad_ru(i  ,j  ,k)+tmp
         ad_ru(i+1,j  ,k)=ad_ru(i+1,j  ,k)+tmp
         ad_ru(i  ,j-1,k)=ad_ru(i  ,j-1,k)+tmp
         ad_ru(i+1,j-1,k)=ad_ru(i+1,j-1,k)+tmp
         ad_v(i,j,k)=0.d0

        END DO

       END DO
      END DO

! == adj to: tl_u
      DO j=J_INTERIOR_U
       DO i=I_INTERIOR_U
               
        cff1=2.d0*mask_u(i,j)/(f(i-1,j)+f(i,j))

        DO k=1,N
! tl:         tl_u(i,j,k)=0.25d0*cff1*(tl_rv(i-1,j+1,k)+tl_rv(i,j+1,k)+
!     *                                tl_rv(i-1,j  ,k)+tl_rv(i,j  ,k))

         tmp=0.25d0*cff1*ad_u(i,j,k)
         ad_rv(i-1,j+1,k)=ad_rv(i-1,j+1,k)+tmp
         ad_rv(i  ,j+1,k)=ad_rv(i  ,j+1,k)+tmp
         ad_rv(i-1,j  ,k)=ad_rv(i-1,j  ,k)+tmp
         ad_rv(i  ,j  ,k)=ad_rv(i  ,j  ,k)+tmp
         ad_u(i,j,k)=0.d0

        END DO

       END DO
      END DO

! == adj to: add BT pres grad term, eta dir:
      DO j=1,M
       DO i=0,L

        cff1=-g*2.d0/(dy(i,j-1)+dy(i,j))

        DO k=1,N
! tl (recurs):         tl_rv(i,j,k)=tl_rv(i,j,k)+tl_cff
         ad_cff=ad_cff+ad_rv(i,j,k)
        END DO

! tl:        tl_cff=cff1*(tl_zeta(i,j)-tl_zeta(i,j-1))
        tmp=cff1*ad_cff
        ad_zeta(i,j)=ad_zeta(i,j)+tmp
        ad_zeta(i,j-1)=ad_zeta(i,j-1)-tmp
        ad_cff=0.d0

       END DO
      END DO

! == adj to: add BT pres grad term to bc pres grad, xi dir:
      DO j=0,M
       DO i=1,L
        cff1=-g*2.d0/(dx(i-1,j)+dx(i,j))

        DO k=1,N
! tl:         tl_ru(i,j,k)=tl_ru(i,j,k)+tl_cff
         ad_cff=ad_cff+ad_ru(i,j,k)
        END DO
        

! tl:        tl_cff=cff1*(tl_zeta(i,j)-tl_zeta(i-1,j))
        tmp=cff1*ad_cff
        ad_zeta(i,j)=ad_zeta(i,j)+tmp
        ad_zeta(i-1,j)=ad_zeta(i-1,j)-tmp
        ad_cff=0.d0     
   
       END DO
      END DO

!================================================
! ADJ to tl_zeta: Solve Poisson equation
!================================================
      q=L+1
      Nb=(L+1)*(M+1)
      allocate(A(q+q+1,Nb),ad_b(Nb))
     
      ad_b=0.d0
      A=0.d0


      DO j=0,M
       DO i=0,L
        ii=index1d(i,j,L)
! tl:        tl_zeta(i,j)=tl_b(ii)
        ad_b(ii)=ad_b(ii)+ad_zeta(i,j)
        ad_zeta(i,j)=0.d0
       END DO
      END DO
      
! == H/dx at u-locations, H/dy at v-locations
      write(*,*) 'H/dx at u-locations, H/dy at v-locations'
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

! == Matrix coefficients (A transposed):
      DO j=0,M
       DO i=0,L

        ! ii = full matrix row number A, col number A*       
        ! jj = full matrix column number A, row num A
        ii=index1d(i,j,L)

        IF (abs(mask_r(i,j))<1.e-8 .or. 
     *          i.eq.0 .or. 
     *          i.eq.L .or.
     *          j.eq.0 .or.
     *          j.eq.M ) THEN
         ! ( masked location, or edge: identity eqn. zeta=0 )
         A(q+1,ii)=1.d0         
        ELSE

         ! coefficient of zeta(i,j):
         jj=index1d(i,j,L)
         A(q+1,ii)=-g*(cffx(i+1,j)+cffx(i,j)+cffe(i,j+1)+cffe(i,j))

         ! coefficient of zeta(i+1,j):
         jj=index1d(i+1,j,L)
         A(jj-ii+q+1,ii)=g*cffx(i+1,j)

         ! coefficient of zeta(i-1,j):
         jj=index1d(i-1,j,L)
         A(jj-ii+q+1,ii)=g*cffx(i,j)

         ! coefficient of zeta(i,j+1):
         jj=index1d(i,j+1,L)
         A(jj-ii+q+1,ii)=g*cffe(i,j+1)

         ! coefficient of zeta(i,j-1):
         jj=index1d(i,j-1,L)
         A(jj-ii+q+1,ii)=g*cffe(i,j)

        END IF

       END DO
      END DO

! == Gauss solver:
      CALL gauss_band(A,ad_b,Nb,q,q)

      DO j=0,M
       DO i=0,L

        ii=index1d(i,j,L)

! tl:        tl_b(ii)=tl_F(i,j)
        ad_F(i,j)=ad_F(i,j)+ad_b(ii)
        ad_b(ii)=0.d0

       END DO
      END DO

!================================================
! ADJ to:
! tl_F=div {f grad(tl_psi)+int_{-H}^0 (tl_ru,tl_rv)}
! intergate div term over the rho grid area, compute boundary fluxes
! shown above in {}. 
! BC: tl_psi=0 , b/c Fx, Fy = 0
!================================================

! == adj to: rhs of the Poisson eq.

      DO j=J_INTERIOR_RHO
       DO i=I_INTERIOR_RHO

! tl 2:        tl_F(i,j)=tl_F(i,j)*mask_r(i,j)
! tl 1:        tl_F(i,j)=tl_Fx(i+1,j)-tl_Fx(i,j)+tl_Fe(i,j+1)-tl_Fe(i,j)
        tmp=ad_F(i,j)*mask_r(i,j)
        ad_Fx(i+1,j)=ad_Fx(i+1,j)+tmp
        ad_Fx(i,j)  =ad_Fx(i,j)  -tmp
        ad_Fe(i,j+1)=ad_Fe(i,j+1)+tmp
        ad_Fe(i,j)  =ad_Fe(i,j)  -tmp
        ad_F(i,j)   =0.d0

       END DO
      END DO

! == Adj to: Fluxes in eta-dir:
      DO j=1,M
       DO i=0,L

! tl 3:
! add b/t (psi) part:
        cff=0.5d0*(f(i,j-1)+f(i,j))*(dx(i,j-1)+dx(i,j))
     *             /(dy(i,j-1)+dy(i,j))

! tl(recur):        tl_Fe(i,j)=tl_Fe(i,j)+cff*(tl_psi(i,j)-tl_psi(i,j-1))
        tmp=ad_Fe(i,j)*cff
        ad_psi(i,j  )=ad_psi(i,j  )+tmp
        ad_psi(i,j-1)=ad_psi(i,j-1)-tmp

! tl 2:
! BC (b/c):
!tl:        tl_Fe(i,j)=tl_Fe(i,j)*mask_v(i,j)
        ad_Fe(i,j)=ad_Fe(i,j)*mask_v(i,j)


! tl 1:
! b/c (tl_Fe has initially been set to 0)
        DO k=N,1,-1
         cff=0.5d0*(Hz(i,j-1,k)+Hz(i,j,k))
! tl:         tl_Fe(i,j)=tl_Fe(i,j)+cff*tl_rv(i,j,k)
         ad_rv(i,j,k)=ad_rv(i,j,k)+cff*ad_Fe(i,j)
        END DO

       END DO
      END DO
 
! == Fluxes in xi-dir:
      DO j=0,M
       DO i=1,L

! tl 3:
! add b/t (psi) part:
        cff=0.5d0*(f(i-1,j)+f(i,j))*(dy(i-1,j)+dy(i,j))
     *             /(dx(i-1,j)+dx(i,j))
! tl:        tl_Fx(i,j)=tl_Fx(i,j)+cff*(tl_psi(i,j)-tl_psi(i-1,j))
        tmp=cff*ad_Fx(i,j)
        ad_psi(i,j  )=ad_psi(i,j  )+tmp
        ad_psi(i-1,j)=ad_psi(i-1,j)-tmp
 
! tl 2:
! BC (b/c):
! tl:        tl_Fx(i,j)=tl_Fx(i,j)*mask_u(i,j)
        ad_Fx(i,j)=ad_Fx(i,j)*mask_u(i,j)

! tl 1:
! b/c (tl_Fx has initially been set to 0)
        DO k=N,1,-1
         cff=0.5d0*(Hz(i-1,j,k)+Hz(i,j,k))
! tl:         tl_Fx(i,j)=tl_Fx(i,j)+cff*tl_ru(i,j,k)
         ad_ru(i,j,k)=ad_ru(i,j,k)+cff*ad_Fx(i,j)
        END DO

       END DO
      END DO

!================================================
! ADJ to b/c pressure grad term:
!================================================


#ifdef PG_BOTTOM_TREATMENT
!================================================
! adj to: reduce b/c pressure grad term near bottom
!================================================
      dist0=20.d0
      dlz=5.d0
      dlz2=dlz*dlz

      DO j=1,M
       DO i=0,L

        z_bot=-0.5d0*(h(i,j)+h(i,j-1))

        DO k=1,N

         z_v=0.5d0*(z_r(i,j,k)+z_r(i,j-1,k))
         dist=z_v-z_bot-dist0

         IF (dist.lt.0.d0) THEN   
          cff=DEXP(-dist*dist/(2.d0*dlz2))

! TL, recursive:          tl_rv(i,j,k)=tl_rv(i,j,k)*cff
          ad_rv(i,j,k)=ad_rv(i,j,k)*cff
         END IF

        END DO

       END DO
      END DO

      DO j=0,M
       DO i=1,L

        z_bot=-0.5d0*(h(i,j)+h(i-1,j))

        DO k=1,N

         z_u=0.5d0*(z_r(i,j,k)+z_r(i-1,j,k))
         dist=z_u-z_bot-dist0

         IF (dist.lt.0.d0) THEN   
          cff=DEXP(-dist*dist/(2.d0*dlz2))
! TL, recursive:          tl_ru(i,j,k)=tl_ru(i,j,k)*cff
          ad_ru(i,j,k)=ad_ru(i,j,k)*cff
         END IF

        END DO

       END DO
      END DO
#endif
     
      CALL ad_prsgrd_SDJ(L,M,N,
     *                ad_rho,z_w,z_r,dx,dy,
     *                mask_u,mask_v,
     *                ad_ru,ad_rv)


!================================================
! ADJ to: tl_RHO 
!================================================

      DO k=1,N
       DO j=J_INTERIOR_RHO
        DO i=I_INTERIOR_RHO

! tl:         tl_rho(i,j,k)=R0*(-Tcoef*tl_T(i,j,k)+Scoef*tl_S(i,j,k))
         ad_T(i,j,k)=ad_T(i,j,k)-Tcoef*R0*ad_rho(i,j,k)
         ad_S(i,j,k)=ad_S(i,j,k)+Scoef*R0*ad_rho(i,j,k)
         ad_rho(i,j,k)=0.d0
  
        END DO
       END DO
      END DO

     
!================================================
! ADJ to: tl_S 
!================================================

#ifdef CONST_DSDT 
      cff=-0.16d0
#endif

      DO k=1,N
       DO j=J_INTERIOR_RHO
        DO i=I_INTERIOR_RHO
         IF (abs(mask_r(i,j)).gt.1.d-10) THEN
 
#ifndef CONST_DSDT
          ! cff=dS/dT
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

! tl:          tl_S(i,j,k)=cff*tl_T(i,j,k)
          ad_T(i,j,k)=ad_T(i,j,k)+cff*ad_S(i,j,k)
          ad_S(i,j,k)=0.d0

         END IF ! mask_r(i,j).ne.0
        END DO
       END DO
      END DO

!================================================
! ADJ to: TL_PSI=0 in exterior and masked locations
!================================================

      write (*,*) 'tl_psi in exterior and masked locations ...'

! tl loop 3:
      DO j=J_INTERIOR_RHO
       DO i=I_INTERIOR_RHO
! tl:         tl_psi(i,j)=tl_psi(i,j)*mask_r(i,j)
        ad_psi(i,j)=ad_psi(i,j)*mask_r(i,j) 
       END DO
      END DO

! tl loop 2:
      DO j=0,M
! tl:       tl_psi(L,j)=0.d0
       ad_psi(L,j)=0.d0
! tl:       tl_psi(0,j)=0.d0
       ad_psi(0,j)=0.d0
      END DO

! tl loop 1:
      DO i=0,L
! tl:       tl_psi(i,M)=0.d0
       ad_psi(i,M)=0.d0
! tl:       tl_psi(i,0)=0.d0
       ad_psi(i,0)=0.d0
      END DO


!================================================
! Adj to READ TL_PSI, TL_T: write ad_psi, ad_T
!================================================
      CALL ad_output_psi_T(univarfile,
     *                     L,M,N,
     *                     time,ad_psi,ad_T)

      END
