#include "cppdefs_cov.h"
      PROGRAM cov_ini_temp_uni

! Overwrites the entries of the inifile: temp is smoothed, psi is zeroed

      implicit none
      include 'netcdf.h'

      character(len=1024) :: inifile,sig_file
      character(len=1024) :: grdfile
      character(len=1)  :: tmpstr

      integer :: status,fid,varid,dimid,fid_sig
      integer :: xi_rho,eta_rho
      integer :: i,i1,j,j1,k,k1,N
      integer :: start(4),count(4)

      real(8) :: dlx,dly,dlz,dllon,dllat,sig,pi,cff,dx,z0,dz,sig_psi

      real(8), allocatable :: lon_r(:)
      real(8), allocatable :: lat_r(:)
      real(8), allocatable :: mask_r(:,:)
      real(8), allocatable :: z_r(:,:,:)

      real(8), allocatable :: psi(:,:)
      real(8), allocatable :: T(:,:,:),sig_T(:,:,:),sig_S(:,:,:)

      real(8), allocatable :: tmp(:,:,:),tmp2(:,:)
      real(8), allocatable :: C(:,:)

      READ (*,'(a1)') tmpstr
      READ (*,'(a80)') grdfile
      READ (*,'(a1)') tmpstr
      READ (*,'(a80)') inifile
      READ (*,'(a1)') tmpstr
      READ (*,*) dlx
      READ (*,'(a1)') tmpstr
      READ (*,*) dly
      READ (*,'(a1)') tmpstr
      READ (*,*) dlz
      READ (*,'(a1)') tmpstr
      READ (*,'(a80)') sig_file
     
      WRITE (*,*) 'GRID FILE:'
      WRITE (*,*)  TRIM(grdfile)
      WRITE (*,*) 'INI FILE:'
      WRITE (*,*)  TRIM(inifile)
      WRITE (*,*) 'SIG FILE:'
      WRITE (*,*)  TRIM(sig_file)
      WRITE (*,*) 'dlx=',dlx
      WRITE (*,*) 'dly=',dly
      WRITE (*,*) 'dlz=',dlz

      start=1
      count=1
      sig_psi=1.0d3

      pi=DACOS(-1.d0)

!================================================
! READ GRD FILE (DIMS, MASK_R, LON_R(:), LAT_R(:))
!================================================
      status=nf_open(TRIM(grdfile),nf_nowrite,fid)
! - INQUIRE DIMENSION LENTHS:
      status=nf_inq_dimid(fid,'xi_rho',dimid)
      status=nf_inq_dimlen(fid,dimid,xi_rho)

      status=nf_inq_dimid(fid,'eta_rho',dimid)
      status=nf_inq_dimlen(fid,dimid,eta_rho)
      
      status=nf_inq_dimid(fid,'N',dimid)
      status=nf_inq_dimlen(fid,dimid,N)


      write (*,*) 'xi_rho, eta_rho, N=',xi_rho,eta_rho,N

      allocate(lon_r(xi_rho))               ! rect. grid
      allocate(lat_r(eta_rho))              ! rect. grid
      allocate(mask_r(xi_rho,eta_rho))
      allocate(z_r(xi_rho,eta_rho,N))

      status=nf_inq_varid(fid,'lon_rho',varid)
      count(1)=xi_rho
      count(2)=1
      status=nf_get_vara_double(fid,varid,start,count,lon_r)

      status=nf_inq_varid(fid,'lat_rho',varid)
      count(1)=1
      count(2)=eta_rho
      status=nf_get_vara_double(fid,varid,start,count,lat_r)

      status=nf_inq_varid(fid,'mask_rho',varid)
      count(1)=xi_rho
      count(2)=eta_rho
      status=nf_get_vara_double(fid,varid,start,count,mask_r)

      status=nf_inq_varid(fid,'z0_r',varid)
      count(1)=xi_rho
      count(2)=eta_rho
      count(3)=N
      status=nf_get_vara_double(fid,varid,start,count,z_r)

      status=nf_close(fid)
! (^^GRD FILE CLOSED^^)

!================================================
! OPEN INI FILE, INITIALIZE VARS AND READ temp:
!================================================
      status=nf_open(TRIM(inifile),nf_nowrite,fid)

! ALLOCATE VARS, INITIALIZE w/ ZEROS:
      allocate(psi(xi_rho,eta_rho))
      allocate(T(xi_rho,eta_rho,N))
      allocate(sig_T(xi_rho,eta_rho,N))
      allocate(sig_S(xi_rho,eta_rho,N))

      psi=0.d0
      T=0.d0

! READ sig_T
      status=nf_open(TRIM(sig_file),nf_nowrite,fid_sig)
      write(*,*) 'status',status
      status=nf_inq_varid(fid_sig,'temp',varid)
      count(1:3)=[xi_rho,eta_rho,N]
      status=nf_get_vara_double(fid_sig,varid,start,count,sig_T)
      status=nf_close(fid_sig)
      WRITE(*,*) 'Min/max sig_temp:',minval(sig_T),maxval(sig_T)
 
! READ PSI:
      status=nf_inq_varid(fid,'psi',varid)
      count(1)=xi_rho
      count(2)=eta_rho
      count(3)=1
      status=nf_get_vara_double(fid,varid,start,count,psi)
      WRITE(*,*) 'Min/max input psi:',minval(psi),maxval(psi)

! READ T:
      status=nf_inq_varid(fid,'temp',varid)
      count(1)=xi_rho
      count(2)=eta_rho
      count(3)=N
      status=nf_get_vara_double(fid,varid,start,count,T)

      status=nf_close(fid)
      WRITE(*,*) 'Min/max input temp:',minval(T),maxval(T)
! (^^INI FILE CLOSED^^)

!================================================
! SMOOTH PSI (add later, now it is zero):
!================================================

      psi=dble(0) !zero variance streamfield

!================================================
! SMOOTH T 
!================================================

! Decorrelation scales in lat and lon:
      dllon=180.d0*dlx/(pi*6.4d+06*DCOS(pi*lat_r(1)/180.d0))
      dllat=180.d0*dly/(pi*6.4d+06)

! Zero variance at masked locations, scale by sig(z) everywhere else:
      DO k=1,N
       DO j=1,eta_rho
        DO i=1,xi_rho
         T(i,j,k)=sig_T(i,j,k)*T(i,j,k)*mask_r(i,j)
        END DO
       END DO
      END DO

! Bell-shape in lon:
      allocate(tmp(xi_rho,eta_rho,N))
      tmp=0.d0

      allocate(C(xi_rho,xi_rho))

      cff=-0.5d0/(dllon*dllon)
      DO i=1,xi_rho  
       DO i1=1,xi_rho
        dx=lon_r(i)-lon_r(i1)        
        C(i,i1)=DEXP( cff*dx*dx )
       END DO
      END DO

      DO k=1,N
       DO j=1,eta_rho
        DO i=1,xi_rho
         DO i1=1,xi_rho
          tmp(i,j,k)=tmp(i,j,k)+C(i,i1)*T(i1,j,k)
         END DO
        END DO
       END DO 
      END DO

      deallocate(C)

! Bell-shape in lat:
      T=0.d0
      allocate(C(eta_rho,eta_rho))

      cff=-0.5d0/(dllat*dllat)
      DO i=1,eta_rho  
       DO i1=1,eta_rho
        dx=lat_r(i)-lat_r(i1)        
        C(i,i1)=DEXP( cff*dx*dx )
       END DO
      END DO

      DO k=1,N
       DO j=1,eta_rho
        DO i=1,xi_rho
         DO j1=1,eta_rho
          T(i,j,k)=T(i,j,k)+C(j,j1)*tmp(i,j1,k)
         END DO
        END DO
       END DO
      END DO

      deallocate(C,tmp)

! Bell-shape in z:
      allocate(tmp(xi_rho,eta_rho,N))
      tmp=0.d0
   
      cff=-0.5d0/(dlz*dlz)

      DO j=1,eta_rho
       DO i=1,xi_rho
       
        DO k=1,N
         DO k1=1,N 
          dx=z_r(203,132,k)-z_r(203,132,k1) !250m depth point 
          tmp(i,j,k)=tmp(i,j,k)+DEXP( cff*dx*dx )*T(i,j,k1)
         END DO
        END DO

       END DO
      END DO

! Repeat zero variance at masked locations, and scale by sig:
      DO k=1,N
       DO j=1,eta_rho
        DO i=1,xi_rho
         T(i,j,k)=sig_T(i,j,k)*tmp(i,j,k)*mask_r(i,j)
        END DO
       END DO
      END DO

!================================================
! WRITE OUTPUTS (OVERWRITE):
!================================================
      status=nf_open(TRIM(inifile),nf_write,fid)

      status=nf_inq_varid(fid,'psi',varid)
      count(1)=xi_rho
      count(2)=eta_rho
      count(3)=1
      status=nf_put_vara_double(fid,varid,start,count,psi)
      WRITE(*,*) 'Min/max output psi:',minval(psi),maxval(psi)

      status=nf_inq_varid(fid,'temp',varid)
      count(1)=xi_rho
      count(2)=eta_rho 
      count(3)=N
      count(4)=1
      status=nf_put_vara_double(fid,varid,start,count,T)
      WRITE(*,*) 'Min/max output temp:',minval(T),maxval(T)

      deallocate(tmp,psi,T)

      status=nf_close(fid)

      WRITE(*,*) 'cov_ini_temp_uni DONE'
      END
