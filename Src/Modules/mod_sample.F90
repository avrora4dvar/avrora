MODULE mod_sample
USE mod_netcdf
USE mod_interp
IMPLICIT NONE

TYPE(sigma_param)::sample_spar
real(8),allocatable::sample_s(:)

real(8),allocatable,dimension(:,:)::sample_lon_r,sample_lat_r,&
&sample_lon_u,sample_lat_u,sample_lon_v,sample_lat_v,&
&sample_h
logical,allocatable,dimension(:,:)::sample_mask_r,sample_mask_u,&
&sample_mask_v

CONTAINS


!-------------------------------------------------------------------------
 SUBROUTINE sample_read_grid(grdfile)
  IMPLICIT NONE
 
  character(len=*)::grdfile
  integer::ncid,status,dimlen(8),gs(2)

  status=nf_open(TRIM(grdfile),nf_nowrite,ncid)
  IF(status.NE.nf_noerr) THEN
      WRITE(*,*) 'mod_sample: cannot open ',TRIM(grdfile)
   STOP
  END IF

  !Size output
  CALL ncsize(ncid,'lon_rho',dimlen)
  gs=dimlen(1:2)


  !Read rho grid
  ALLOCATE(sample_lon_r(gs(1),gs(2)))
  ALLOCATE(sample_lat_r(gs(1),gs(2)))
  ALLOCATE(sample_h(gs(1),gs(2)))
  ALLOCATE(sample_mask_r(gs(1),gs(2)))
  sample_lon_r=ncread2d(ncid,'lon_rho',[1,1],gs)
  sample_lat_r=ncread2d(ncid,'lat_rho',[1,1],gs)
  sample_h=ncread2d(ncid,'h',[1,1],gs)
  sample_mask_r=(ncread2d(ncid,'mask_rho',[1,1],gs).NE.0)

  !Read u grid
  ALLOCATE(sample_lon_u(gs(1)-1,gs(2)))
  ALLOCATE(sample_lat_u(gs(1)-1,gs(2)))
  ALLOCATE(sample_mask_u(gs(1)-1,gs(2)))
  sample_lon_u=ncread2d(ncid,'lon_u',[1,1],[gs(1)-1,gs(2)])
  sample_lat_u=ncread2d(ncid,'lat_u',[1,1],[gs(1)-1,gs(2)])
  sample_mask_u=(ncread2d(ncid,'mask_u',[1,1],[gs(1)-1,gs(2)]).NE.0)

  !Read v grid
  ALLOCATE(sample_lon_r(gs(1),gs(2)-1))
  ALLOCATE(sample_lat_r(gs(1),gs(2)-1))
  ALLOCATE(sample_mask_r(gs(1),gs(2)-1))
  sample_lon_v=ncread2d(ncid,'lon_v',[1,1],[gs(1),gs(2)-1])
  sample_lat_v=ncread2d(ncid,'lat_v',[1,1],[gs(1),gs(2)-1])
  sample_mask_v=(ncread2d(ncid,'mask_v',[1,1],[gs(1),gs(2)-1]).NE.0)

 END SUBROUTINE sample_read_grid
!--------------------------------------------------------------------------
 SUBROUTINE sample_index1d(x,x1,ix,w)
 !Get indices of x-coordinates lower and higher than x1 and the weights for
 !linear interpolation

  IMPLICIT NONE
  real(8),intent(in)::x(:),x1
  integer::ix(2)
  real(8)::w(2)

  IF(x1.LT.MINVAL(x,1)) THEN
   ix=MINLOC(x,1); w=dble(.5)
  ELSEIF(x1.GE.MAXVAL(x,1)) THEN
   ix=MAXLOC(x,1); w=dble(.5)
  ELSE
   ix(1)=MAXLOC(x,1,x.LE.x1)
   ix(2)=MINLOC(x,1,x.GT.x1)
   w(2)=(x1-x(ix(1)))/(x(ix(2))-x(ix(1)))
   w(1)=dble(1)-w(2)
  END IF
 END SUBROUTINE sample_index1d

!---------------------------------------------------------------------------
 SUBROUTINE sample_index2d(x,y,mask,x1,y1,ix,iy,w)
  !Get indices of the four grid points surrounding x1,y1 with their weights
  !for bilinear interpolation
  IMPLICIT NONE
  real(8),intent(in)::x(:),y(:),x1,y1
  logical,intent(in)::mask(:,:)
  integer::ix(4),iy(4),ix1(2),iy1(2),i0,i1,i2
  real(8)::w(4),wx1(2),wy1(2)

  CALL sample_index1d(x,x1,ix1,wx1)
  CALL sample_index1d(y,y1,iy1,wy1)
  ix=[ix1(1),ix1(2),ix1(1),ix1(2)]
  iy=[iy1(1),iy1(1),iy1(2),iy1(2)]
  
  DO i0=1,4
   i1=MOD(i0-1,2)+1; i2=(i0-1)/2+1
   IF(mask(i1,i2)) THEN
    w(i0)=wx1(i1)*wy1(i2)
   ELSE
    w(i0)=dble(0)
   END IF
  END DO

  IF(SUM(w).GT.dble(0)) THEN
   w=w/SUM(w)
  END If
 END SUBROUTINE sample_index2d

!---------------------------------------------------------------------------
 SUBROUTINE sample_avg2d(x,y,mask,x1,y1,dl,ix,iy,w)
  !Generates weigths for 2d averaging

  IMPLICIT NONE
  real(8),intent(in)::x(:),y(:),x1,y1
  real(8),intent(in)::dl(2)
  logical,intent(in)::mask(:,:)
  integer::i0,i1,i2,ni(2,3)
  real(8),allocatable::w(:,:)
  integer,allocatable::ix(:,:),iy(:,:)
  real(8)::pi

  pi=ACOS(-1.)

  !Find cells in area
  ni(1,1)=MINLOC(x,1,x.GE.x1-dl(1))
  ni(1,2)=MAXLOC(x,1,x.LE.x1+dl(1))
  ni(2,1)=MINLOC(y,1,y.GE.y1-dl(2))
  ni(2,2)=MAXLOC(y,1,y.LE.y1+dl(2))
    
  !Add points outside boundary if possible
  ni(1,1)=MAXVAL([1,ni(1,1)-1])
  ni(1,2)=MINVAL([size(x),ni(1,2)+1])
  ni(2,1)=MAXVAL([1,ni(2,1)-1])
  ni(2,2)=MINVAL([size(y),ni(2,2)+1])
  ni(1,3)=ni(1,2)-ni(1,1)+1
  ni(2,3)=ni(2,2)-ni(2,1)+1

  !Calculate area grid cells
  ALLOCATE(w(ni(1,3),ni(2,3))); w=dble(0)
  DO i2=ni(2,1)+1,ni(2,2)-1
   w(2:ni(1,3)-1,i2-ni(2,1)+1)=.5*(y(i2+1)-y(i2-1))
  END DO
  DO i1=ni(1,1)+1,ni(1,2)-1
   w(i1-ni(1,1)+1,:)=w(i1-ni(1,1)+1,:)*.5*(x(i1+1)-x(i1-1))*&
   &cos(pi*y(ni(2,1):ni(2,2))/180.)
  END DO
 

  !Nomalize weights
  w=ABS(w)
  WHERE(.NOT.mask(ni(1,1):ni(1,2),ni(2,1):ni(2,2)))
   w=dble(0)
  END WHERE
  IF(SUM(w).GT.dble(0)) w=w/SUM(w)

  !Write output
  ALLOCATE(ix(ni(1,3),ni(2,3))); ALLOCATE(iy(ni(1,3),ni(2,3)))
  ix=0; iy=0
  DO i2=1,size(w,2)
  DO i1=1,size(w,1)
   IF(w(i1,i2).EQ.dble(0)) CYCLE
   ix(i1,i2)=(i1-1)+ni(1,1); iy(i1,i2)=(i2-1)+ni(2,1)
  END DO 
  END DO
  
 END SUBROUTINE sample_avg2d

!---------------------------------------------------------------------------
 SUBROUTINE sample_depthAvg(spar,h,zeta,z1,dl,iz,w)
  IMPLICIT NONE
  TYPE(sigma_param),intent(in)::spar
  real(8),intent(in)::h,zeta,z1,dl
  real(8),allocatable::w(:,:),sw(:),sr(:),zw(:),zr(:)
  integer,allocatable::iz(:,:)
  integer::ni(3),i0,i1
  real(8)::bnd(2)

  !sigma coordinates
  ALLOCATE(sr(spar%n_sigma))
  DO i1=1,size(sr)
   sr(i1)=dble(-1)+(dble(i1)-0.5)/dble(spar%n_sigma)
  END DO

  !z-coordinates
  ALLOCATE(zr(size(sr)))
  zr=sigma2z(spar,h,zeta,sr)
  
  !Convert to depths
  zr=zeta-zr

  !Integrated linear interpolant
  ALLOCATE(w(1,spar%n_sigma)); w=dble(0)
  DO i1=1,size(zr)-1
   bnd(2)=MINVAL([z1+dl,zr(i1)]) !largest depth
   bnd(1)=MAXVAL([z1-dl,zr(i1+1)]) !smallest depth
   IF(bnd(2).LE.bnd(1)) CYCLE
   w(1,i1+1)=w(1,i1+1)+(bnd(2)-bnd(1))/(zr(i1)-zr(i1+1))*&
   &(zr(i1)-.5*(bnd(1)+bnd(2)))
   w(1,i1)=w(1,i1)-(bnd(2)-bnd(1))/(zr(i1)-zr(i1+1))*&
   &(zr(i1+1)-.5*(bnd(1)+bnd(2)))
  END DO
  
  !Add bottom and surface
  w(1,1)=w(1,1)+MAXVAL([dble(0),MINVAL([z1+dl,h])-MAXVAL([z1-dl,zr(1)])])
  w(1,size(w))=w(1,size(w))+MAXVAL([dble(0),&
  &MINVAL([z1+dl,zr(size(zr))])-MAXVAL([z1-dl,dble(0)])])
  
  !Normalize
  IF(SUM(w).GT.dble(0)) w=w/SUM(w)
   
  !Output
  ALLOCATE(iz(1,spar%n_sigma)); iz=0
  DO i1=1,size(iz)
   IF(w(1,i1).GT.dble(0)) iz(1,i1)=i1
  END DO
   
 END SUBROUTINE sample_depthAvg
!---------------------------------------------------------------------------
 SUBROUTINE sample_index_delta(t,t1,it,w)
 !Get time indices it and their weights for sampling a series given at times
 !t at time t1
  IMPLICIT NONE
  real(8),intent(in)::t(:),t1
  integer::it(:),itmin,itmax
  real(8)::w(:)  

  it=0; w=dble(0)

  IF(t1.LT.MINVAL(t)) THEN
   IF(size(it).EQ.2) THEN
    it=MINLOC(t); w=dble(.5)
   ELSE
    it(MINLOC(t))=MINLOC(t); w(MINLOC(t))=dble(0)
   END IF
  ELSEIF(t1.GE.MAXVAL(t)) THEN
   IF(size(it).EQ.2) THEN
    it=MAXLOC(t); w=dble(.5)
   ELSE
    it(MAXLOC(t))=MAXLOC(t); w(MAXLOC(t))=dble(1)
   END IF 
  ELSE
   itMin=MAXLOC(t,1,t.LE.t1)
   itMax=MINLOC(t,1,t.GT.t1)
   IF(size(it).EQ.2) THEN
    it(1)=itMin; it(2)=itMax; 
    w(2)=(t1-t(itMin))/(t(itMax)-t(itMin))
    w(1)=dble(1)-w(2)
   ELSE
    it(itMin)=itMin; it(itMax)=itMax
    w(itMax)=(t1-t(itMin))/(t(itMax)-t(itMin))
    w(itMin)=dble(1)-w(itMax)
   END IF
  END IF
 END SUBROUTINE sample_index_delta

!---------------------------------------------------------------------------
 SUBROUTINE sample_index_mean(t,t1,dt,it,w)
 !Get time indices and their weights for calculating the time average over
 !the period [t1-.5*dt,t1+.5*dt]
  IMPLICIT NONE
  real(8),intent(in)::t(:),t1,dt
  integer::it(:)
  real(8)::w(:)
  integer::i0,itmin,itmax,nt
  real(8)::tmin,tmax

  IF(t1-.5*dt.LT.MINVAL(t)) THEN
   itmin=MINLOC(t,1)
  ELSE
   itmin=MAXLOC(t,1,t.LE.t1-.5*dt)
  END IF
  IF(t1+.5*dt.GT.MAXVAL(t)) THEN
   itmax=MAXLOC(t,1)
  ELSE
   itmax=MINLOC(t,1,t.GE.t1+.5*dt)
  END IF

  it=0; w=dble(0)

  DO i0=itmin,itmax-1
   tmin=MAX(t1-.5*dt,t(i0))
   tmax=MIN(t1+.5*dt,t(i0+1))
   IF(tmin.GE.tmax) CYCLE
   it(i0)=i0; it(i0+1)=i0+1
   w(i0)=w(i0)+dble(1)/(t(i0+1)-t(i0))*&
   &( t(i0+1)*(tmax-tmin)-.5*(tmax**2-tmin**2) )
   w(i0+1)=w(i0+1)+dble(1)/(t(i0+1)-t(i0))*&
   &( .5*(tmax**2-tmin**2)-(tmax-tmin)*t(i0) )
  END DO
  IF(SUM(w).GT.dble(0)) w=w/SUM(w)
  
 END SUBROUTINE sample_index_mean

!------------------------------------------------------------------------------
 SUBROUTINE sample_index_lanczos(t,t1,dt,fcut,it,w)
 !Finds time indices and their weighting to determine the value of a time
 !series after low-pass filtering with a cos-Lanczos filter at time t1. Here
 !the filter is zero outside the range [t1-.5*dt,t1+.5*dt] and with cut-off
 !frequency fcut
  IMPLICIT NONE
  real(8),intent(in)::t(:),t1,dt,fcut
  real(8)::w(:)
  integer::it(:)
  integer::i0,itmin,itmax,nt
  real(8)::tmin,tmax,tt1
  real(8)::pi

  pi=ACOS(-1.)

  IF(t1-.5*dt.LT.MINVAL(t)) THEN
   itmin=MINLOC(t,1)
  ELSE
   itmin=MAXLOC(t,1,t.LE.t1-.5*dt)
  END IF
  IF(t1+.5*dt.GT.MAXVAL(t)) THEN
   itmax=MAXLOC(t,1)
  ELSE
   itmax=MINLOC(t,1,t.GE.t1+.5*dt)
  END IF

  it=0; w=dble(0)

  DO i0=itmin,itmax-1
   tmin=MAX(t1-.5*dt,t(i0))
   tmax=MIN(t1+.5*dt,t(i0+1))
   IF(tmin.GE.tmax) CYCLE
   it(i0)=i0; it(i0+1)=i0+1
   
   tt1=ABS(t(i0)-t1)
   IF(tt1.GE..5*dt) THEN
    w(i0)=w(i0)+dble(0)
   ELSEIF(tt1.EQ.dble(0)) THEN
    w(i0)=w(i0)+dble(.5)*(tmax-tmin)
   ELSE
    w(i0)=w(i0)+dble(.5)*(tmax-tmin)*&
    .5*(1.+COS(2*pi*tt1/dt))*SIN(pi*tt1*fcut)/(pi*tt1*fcut)
   END IF

   tt1=ABS(t(i0+1)-t1)
   IF(tt1.GE..5*dt) THEN
    w(i0+1)=w(i0+1)+dble(0)
   ELSEIF(tt1.EQ.dble(0)) THEN
    w(i0+1)=w(i0+1)+dble(.5)*(tmax-tmin)
   ELSE
    w(i0+1)=w(i0+1)+dble(.5)*(tmax-tmin)*&
    .5*(1.+COS(2*pi*tt1/dt))*SIN(pi*tt1*fcut)/(pi*tt1*fcut)
   END IF

  END DO
  IF(SUM(w).GT.dble(0)) w=w/SUM(w)
 
 END SUBROUTINE sample_index_lanczos
!-------------------------------------------------------------------------
 SUBROUTInE sample_index_harmonic(t,t1,f,it,w)
  IMPLICIT NONE
  real(8),intent(in)::t(:),t1,f(:)
  integer::it(:),lwork,bnd(2)
  real(8)::w(:),pi,wmean
  real(8),allocatable::VM(:,:),WM(:,:),tau(:),work(:)
  integer::i0,i1,i2,i3,status

  pi=ACOS(-1.0)

  !Construct Vandermonde matrix
  ALLOCATE(VM(size(t),2*size(f)+3))
  VM(:,1)=dble(1); VM(:,2)=t-t(1); VM(:,3)=(t-t(1))**2
  DO i2=1,size(f)
   VM(:,2*i2+2)=COS(2*pi*t*f(i2))
   VM(:,2*i2+3)=SIN(2*pi*t*f(i2))
  END DO

  !Construct linear interpolation weight matrix
  ALLOCATE(WM(size(t),1)); WM=dble(0)
  DO i2=1,size(WM,2)
   IF(t1.LE.MINVAL(t)) THEN
     WM(MINLOC(t),i2)=dble(1)
   ELSEIF(t1.GE.MAXVAL(t)) THEN
     WM(MAXLOC(t),i2)=dble(1)
   ELSE
     bnd(1)=MAXLOC(t,1,t.LE.t1)
     bnd(2)=MINLOC(t,1,t.GT.t1)
     WM(bnd(2),i2)= (t1-t(bnd(1)))/(t(bnd(2))-t(bnd(1)))
     WM(bnd(1),i2)= (t(bnd(2))-t1)/(t(bnd(2))-t(bnd(1)))
   END IF
  END DO
  w=WM(:,1)

  !Calculate QR-factorization VM
  ALLOCATE(tau(size(VM,2)))
  lwork=-1; ALLOCATE(work(1))
  CALL DGEQRF(size(VM,1),size(VM,2),VM,size(VM,1),tau,&
  &work,lwork,status)
  lwork=work(1); DEALLOCATE(work); ALLOCATE(work(lwork))
  CALL DGEQRF(size(VM,1),size(VM,2),VM,size(VM,1),tau,&
  &work,lwork,status)
  DEALLOCATE(work)
 
  !Calculate Q'w
  lwork=-1; ALLOCATE(work(1))
  CALL DORMQR('L','T',size(WM,1),size(WM,2),size(VM,2),&
  &VM,size(VM,1),tau,WM,size(WM,1),work,lwork,status)
  lwork=work(1); DEALLOCATE(work); ALLOCATE(work(lwork))
  CALL DORMQR('L','T',size(WM,1),size(WM,2),size(VM,2),&
  &VM,size(VM,1),tau,WM,size(WM,1),work,lwork,status)
  DEALLOCATE(work)
  WM(1:3,:)=dble(0)
  WM(4+2*size(f):size(WM,1),:)=dble(0)

 !Calculate QQ'w
  ALLOCATE(work(1)); lwork=-1
  CALL DORMQR('L','N',size(WM,1),size(WM,2),size(Vm,2),&
  &VM,size(VM,1),tau,WM,size(WM,1),work,lwork,status)
  lwork=work(1); DEALLOCATE(work); ALLOCATE(work(lwork))
  CALL DORMQR('L','N',size(WM,1),size(WM,2),size(Vm,2),&
  &VM,size(VM,1),tau,WM,size(WM,1),work,lwork,status)
  DEALLOCATE(work)

  !Calculate (I-QQ')w
  w=w-WM(:,1)
 
  DO i1=1,size(it)
   it(i1)=i1
  END DO
 END SUBROUTINE sample_index_harmonic

!-------------------------------------------------------------------------
 SUBROUTINE sample_index_depth(spar,h,zeta,depth,iz,wz)
 !sample_index_depth calculates given the s-grid parameters spar,
 !depth h and sea-surface height zeta the vertical indices below and 
 !above the depth depth and the weights for linear interpolation. Here
 !depth is defined and the distance below the sea-surface

   IMPLICIT NONE
   TYPE(sigma_param),intent(in)::spar
   real(8),intent(in)::h,zeta,depth
   integer::iz(2),i3
   real(8)::wz(2)
   real(8),allocatable::s(:),z(:)

   !sigma coordinates
   ALLOCATE(s(spar%n_sigma))
   DO i3=1,size(s)
    s(i3)=dble(-spar%n_sigma+i3-.5)/dble(spar%n_sigma)
   END DO

   !z-coordinates
   ALLOCATE(z(spar%n_sigma))
   z=sigma2z(spar,h,zeta,s)
   z=zeta-z

   IF(depth.LE.MINVAL(z)) THEN
    iz=MINLOC(z,1); wz=dble(.5)
   ELSEIF(depth.GE.MAXVAL(z)) THEN
    iz=MAXLOC(z,1); wz=dble(.5)
   ELSE
    iz(1)=MINLOC(z,1,z.GE.depth)
    iz(2)=MAXLOC(z,1,z.LT.depth)
    wz(2)=(depth-z(iz(1)))/(z(iz(2))-z(iz(1)))
    wz(1)=dble(1)-wz(2)
   END IF

 END SUBROUTINE
!-------------------------------------------------------------------------
 SUBROUTINE sample_index_z(spar,h,zeta,z1,iz,wz)
 !sample_index_depth calculates given the s-grid parameters spar,
 !depth h and sea-surface height zeta the vertical indices below and 
 !above the vertical position z1 and the weights for linear interpolation. 
 !Here z1 is defined as the position with respect to the reference surface 
 !surface with z1 becoming negative in direction of the bottom.

   IMPLICIT NONE
   TYPE(sigma_param),intent(in)::spar
   real(8),intent(in)::h,zeta,z1
   integer::iz(2),i3
   real(8)::wz(2)
   real(8),allocatable::s(:),z(:)

   !sigma coordinates
   ALLOCATE(s(spar%n_sigma))
   DO i3=1,size(s)
    s(i3)=dble(-spar%n_sigma+i3-.5)/dble(spar%n_sigma)
   END DO

   !z-coordinates
   ALLOCATE(z(spar%n_sigma))
   z=sigma2z(spar,h,zeta,s)

   IF(z1.LE.MINVAL(z)) THEN
    iz=MINLOC(z,1); wz=dble(.5)
   ELSEIF(z1.GE.MAXVAL(z)) THEN
    iz=MAXLOC(z,1); wz=dble(.5)
   ELSE
    iz(1)=MAXLOC(z,1,z.LE.z1)
    iz(2)=MINLOC(z,1,z.GT.z1)
    wz(2)=(z1-z(iz(1)))/(z(iz(2))-z(iz(1)))
    wz(1)=dble(1)-wz(2)
   END IF

 END SUBROUTINE

END MODULE


