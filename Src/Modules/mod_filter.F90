MODULE mod_filter

 USE mod_interp
 IMPLICIT NONE

 INTERFACE filter_gauss_conv
  module procedure filter_gauss_conv1D,filter_gauss_conv2D,&
  &filter_gauss_conv3D
 END INTERFACE

 INTERFACE filter_gauss_deconv
  module procedure filter_gauss_deconv1D,filter_gauss_deconv2D,&
  &filter_gauss_deconv3D
 END INTERFACE

 REAL(8),PARAMETER::filter_pi=3.141592653589793

 INTEGER::nThreads

 CONTAINS

!--------------------------------------------------------------------------

SUBROUTINE filter_coslanczos(time_out,time,w)
 !Calculate filter weights for a cosine-Lanczos filter (Mooers et al. 1968).
 !Time has to be given in days

 IMPLICIT NONE
 REAL(8),intent(in)::time(:),time_out
 REAL(8),intent(out)::w(lbound(time,1):ubound(time,1))
 REAL::fw(2),t(lbound(time,1):ubound(time,1)),&
 &dt(lbound(time,1):ubound(time,1))

 !filter-characteristics
 fw(1)=dble(61)/dble(24)
 fw(2)=dble(40)/dble(24) 

 t=time-time_out

 !Time steps for trapezium rule integration
 dt=dble(0)
 dt(:ubound(dt,1)-1)=dt(:ubound(dt,1)-1)+&
 &.5*(t(lbound(dt,1)+1:ubound(dt,1))-t(lbound(dt,1):ubound(dt,1)-1))
 dt(lbound(dt,1)+1:)=dt(lbound(dt,1)+1:)+&
 &.5*(t(lbound(dt,1)+1:ubound(dt,1))-t(lbound(dt,1):ubound(dt,1)-1))

 !Calculate weight
 w=dble(0.5)*(dble(1)+COS(filter_pi*t/fw(1)))*&
 &SIN(2*filter_pi*t/fw(2))/&
 &(2*filter_pi*t/fw(2))
 WHERE(ABS(t).LT.1e-5)
  w=dble(1)
 END WHERE

 !Normalize
 w=w*dt
 w=w/SUM(w)
 WHERE(ABS(t).GT.fw(1))
  w=dble(0)
 END WHERE
    
END SUBROUTINE filter_coslanczos

!---------------------------------------------------------------------------

 FUNCTION filter_gauss_conv1D(x,y,ndim,dl,mask_in) RESULT(yi)
 !FILTER_GAUSS_CONV1D Smoothens 1D signal with a Gaussian distribution.
 !Does not correct for boundary effects as to preserve symmetry
 !NB1: for symmetric matrices the spacing of x must be uniform
 !NB2: dl>= grid spacing in x
 !
 !Syntax:
 ! yi=filter_gauss_conv1d(x,y,dl,mask)
 !Input:
 ! x: n double array with coordinates
 ! y: n double array array with values of signal at coordinates x
 ! ndim: dummy argument
 ! dl: double with the standard deviation of Gauss funtion
 ! mask (optional): n logical array with points that should be included
 !                  in the computation
 !Output:
 ! yi: n double array with convoluted signal at coordinates x 

 REAL(8),INTENT(in)::x(:),y(:),dl
 REAL(8)::yi(lbound(y,1):ubound(y,1))
 LOGICAL,INTENT(inout),OPTIONAL::mask_in(lbound(x,1):ubound(x,1))
 LOGICAL::mask(lbound(x,1):ubound(x,1))
 REAL(8)::kernel(lbound(x,1):ubound(x,1))
 REAL(8)::dx(lbound(x,1):ubound(x,1))
 INTEGER::i1,j1,ndim
  
 IF(.NOT.PRESENT(mask_in)) THEN
  mask=.TRUE.
 ELSE
  mask=mask_in
 END IF

 !Initialize
 yi=DBLE(0); kernel=DBLE(0); 
 
 DO i1=lbound(x,1)+1,ubound(x,1)-1
  dx(i1)=0.5*(x(i1+1)-x(i1-1))
 END DO
 dx(lbound(x,1))=dx(lbound(x,1)+1)
 dx(ubound(x,1))=dx(ubound(x,1)-1)
 dx=dx/SQRT(2*filter_pi*dl**2)

 !Convolute
 DO i1=lbound(yi,1),ubound(yi,1)
  IF(.NOT.mask(i1)) CYCLE
  kernel=EXP(-0.5*(x(i1)-x)**2/dl**2) !*dx
  yi(i1)=sum(y*kernel,mask)
 END DO

 WHERE(.NOT.mask)
  yi=y
 END WHERE
 !write(*,*) yi(20), yi(22) 

 END FUNCTION filter_gauss_conv1D

!----------------------------------------------------------------------------

 FUNCTION filter_gauss_deconv1D(x,y,ndim,dl,order,mask_in) RESULT(yi)
 !FILTER_GAUSS_DECONV1D Deconvolutes a signal convoluted with 
 !filter_gauss_conv1D
 !NB1: for symmetric matrices the spacing of x must be uniform
 !NB2: dl>= grid spacing in x
 !
 !Syntax:
 ! yi=filter_gauss_conv1d(x,y,dl,mask)
 !Input:
 ! x: n double array with coordinates
 ! y: n double array array with values of signal at coordinates x
 ! ndim: dummy argument
 ! dl: double with the standard deviation of Gauss funtion
 ! mask (optional): n logical array with points that should be included
 !                  in the computation
 ! order: integer with the order of the approximation of the inverse
 !        convolution operator normally 4-6 will do. 
 !Output:
 ! yi: n double array with convoluted signal at coordinates x 

 INTEGER,INTENT(in)::order
 REAL(8),INTENT(in)::x(:),y(:),dl
 REAL(8)::yi(lbound(y,1):ubound(y,1))
 LOGICAL,INTENT(inout),OPTIONAL::mask_in(lbound(x,1):ubound(x,1))
 LOGICAL::mask(lbound(x,1):ubound(x,1))
 REAL(8)::kernel(lbound(x,1):ubound(x,1))
 REAL(8)::dx(lbound(x,1):ubound(x,1))
 INTEGER::i1,j1,ndim,n,n_fac
  
 IF(.NOT.PRESENT(mask_in)) THEN
  mask=.TRUE.
 ELSE
  mask=mask_in
 END IF

 !Initialize
 yi=DBLE(0); kernel=DBLE(0); 
 
 DO i1=lbound(x,1)+1,ubound(x,1)-1
  dx(i1)=0.5*(x(i1+1)-x(i1-1))
 END DO
 dx(lbound(x,1))=dx(lbound(x,1)+1)
 dx(ubound(x,1))=dx(ubound(x,1)-1)
 dx=dx/SQRT(2*filter_pi*dl**2)

 !Convolute
 DO i1=lbound(yi,1),ubound(yi,1)
  IF(.NOT.mask(i1)) CYCLE

  !Calculate deconvolution operator at these differences
  kernel=DBLE(0); n_fac=1
  DO n=0,order
   IF(n.GT.0) n_fac=n_fac*n
   kernel=kernel+filter_diff_gauss(x(i1)-x,-0.5/dl**2,2*n)*(-dl*dl)**n&
   &/DBLE(n_fac)
  END DO
  kernel=kernel*dx

  yi(i1)=sum(y*kernel,mask)
 END DO

 WHERE(.NOT.mask)
  yi=y
 END WHERE

 END FUNCTION filter_gauss_deconv1D

!---------------------------------------------------------------------------

 FUNCTION filter_gauss_conv2D(x,y,ndim,dl,mask_in) RESULT(yi)
 !FILTER_GAUSS_CONV2D As filter_gauss_conv1d but now for a 2D seperable grid
 !
 !Syntax:
 ! yi=filter_gauss_conv2d(x,y,ndim,dl,mask)
 !Input:
 ! x: n/m double array with coordinates along the 1st dimension
 ! y:  nxm double array with values at coordinates in x1,x2
 ! ndim: integer with dimension along which convolution takes place
 ! dl: double with standard deviations along 1st,2nd direction
 ! mask: nxm logical array with points that should be included in calculation
 !Output
 ! yi: nxm double array with convoluted values at coordinates in x
 !NB1: for symmetric matrices grid spacing must be uniform
 !NB2: dl>= grid spacing in respectively x

 REAL(8),INTENT(in)::x(:),y(:,:),dl
 INTEGER,INTENT(in)::ndim
 LOGICAL,INTENT(in),OPTIONAL::mask_in(lbound(y,1):ubound(y,1),&
 &lbound(y,2):ubound(y,2))
 LOGICAL::mask(lbound(y,1):ubound(y,1),&
 &lbound(y,2):ubound(y,2))
 REAL(8)::yi(lbound(y,1):ubound(y,1),lbound(y,2):ubound(y,2))
 REAL(8)::kernel(lbound(x,1):ubound(x,1)),&
 &dx(lbound(x,1):ubound(x,1)),&
 &y_tmp(lbound(y,1):ubound(y,1),lbound(y,2):ubound(y,2))
 INTEGER::i1,i2,j1,j2

 IF(PRESENT(mask_in)) THEN
  mask=mask_in
 ELSE
  mask=.TRUE.
 END IF

 !Calculate grid spacing
 DO i1=lbound(x,1)+1,ubound(x,1)-1
  dx(i1)=0.5*(x(i1+1)-x(i1-1))
 END DO
 dx(lbound(x,1))=dx(lbound(x,1)+1)
 dx(ubound(x,1))=dx(ubound(x,1)-1)
 dx=dx/SQRT(2*filter_pi*dl**2)
 
 yi=DBLE(0); y_tmp=DBLE(0)
  !write(*,*) 'conv'
  !write(*,*) lbound(yi,1),ubound(yi,1),lbound(yi,2),ubound(yi,2)
 IF(ndim.EQ.1) THEN
  !Convolute along 1st dimension

  IF(SIZE(x,1).NE.SIZE(y,1)) THEN
   WRITE(*,*) 'mod_filter.filter_gauss_conv2D: incongruent dimensions&
   & x,y'
   STOP
  END IF
 
  DO i1=lbound(yi,1),ubound(yi,1)
   !Calculate kernel
   kernel=EXP(-.5*(x(i1)-x)**2/dl**2) !*dx
 
   !Perform convolution
   DO j1=lbound(y,1),ubound(y,1)
    y_tmp(j1,:)=y(j1,:)*kernel(j1)
   END DO
   yi(i1,:)=SUM(y_tmp,1,mask)

  END DO !1st dimension
 ELSEIF(ndim.EQ.2) THEN
  !Convolute along 2nd dimension

  IF(SIZE(x,1).NE.SIZE(y,2)) THEN
   WRITE(*,*) 'mod_filter.filter_gauss_conv2D: incongruent dimensions&
   & x and y'
   STOP
  END IF

  DO i2=lbound(yi,2),ubound(yi,2)
   !Calculate kernel
   kernel=EXP(-.5*(x(i2)-x)**2/dl**2) !*dx

   !Perform convolution
   DO j2=lbound(y,2),ubound(y,2)
    y_tmp(:,j2)=y(:,j2)*kernel(j2)
   END DO
   yi(:,i2)=SUM(y_tmp,2,mask)

  END DO !2nd dimension
 ELSE
  WRITE(*,*) 'mod_filter.gauss_conv2D: ndim must be 1 or 2'
  STOP
 END IF

 !Normalize
 WHERE(.NOT.mask)
  yi=y
 END WHERE

 END FUNCTION filter_gauss_conv2D

!----------------------------------------------------------------------------

FUNCTION filter_gauss_deconv2D(x,y,ndim,dl,order,mask_in) RESULT(yi)
 !FILTER_GAUSS_DECONV2D deconvolutes a signal convoluted with filter_gauss
 !conv2D
 !Syntax:
 ! yi=filter_gauss_conv2d(x,y,ndim,dl,mask)
 !Input:
 ! x: n/m double array with coordinates along the 1st dimension
 ! y:  nxm double array with values at coordinates in x
 ! ndim: integer with dimension along which deconvolution takes place
 ! dl: double with standard deviations
 ! mask: nxm logical array with points that should be included in calculation
 ! order: order of the approximation of the deconvolution operator. Usually
 !        order 4-6 will do. 
 !Output
 ! yi: nxm double array with convoluted values at coordinates in x
 !NB1: for symmetric matrices grid spacing must be uniform
 !NB2: dl>= grid spacing in respectively x

 INTEGER,INTENT(in)::order
 REAL(8),INTENT(in)::x(:),y(:,:),dl
 INTEGER,INTENT(in)::ndim
 LOGICAL,INTENT(in),OPTIONAL::mask_in(lbound(y,1):ubound(y,1),&
 &lbound(y,2):ubound(y,2))
 LOGICAL::mask(lbound(y,1):ubound(y,1),&
 &lbound(y,2):ubound(y,2))
 REAL(8)::yi(lbound(y,1):ubound(y,1),lbound(y,2):ubound(y,2))
 REAL(8)::kernel(lbound(x,1):ubound(x,1)),&
 &dx(lbound(x,1):ubound(x,1)),&
 &y_tmp(lbound(y,1):ubound(y,1),lbound(y,2):ubound(y,2))
 INTEGER::i1,i2,j1,j2,n_fac,n

 IF(PRESENT(mask_in)) THEN
  mask=mask_in
 ELSE
  mask=.TRUE.
 END IF

 !Calculate grid spacing
 DO i1=lbound(x,1)+1,ubound(x,1)-1
  dx(i1)=0.5*(x(i1+1)-x(i1-1))
 END DO
 dx(lbound(x,1))=dx(lbound(x,1)+1)
 dx(ubound(x,1))=dx(ubound(x,1)-1)
 dx=dx/SQRT(2*filter_pi*dl**2)
 
 yi=DBLE(0); y_tmp=DBLE(0)
  !write(*,*) 'conv'
  !write(*,*) lbound(yi,1),ubound(yi,1),lbound(yi,2),ubound(yi,2)
 IF(ndim.EQ.1) THEN
  !Convolute along 1st dimension

  IF(SIZE(x,1).NE.SIZE(y,1)) THEN
   WRITE(*,*) 'mod_filter.filter_gauss_conv2D: incongruent dimensions&
   & x,y'
   STOP
  END IF
 
  DO i1=lbound(yi,1),ubound(yi,1)

   !Calculate deconvolution operator at these differences
   kernel=DBLE(0); n_fac=1
   DO n=0,order
    IF(n.GT.0) n_fac=n_fac*n
    kernel=kernel+filter_diff_gauss(x(i1)-x,-0.5/dl**2,2*n)*(-dl*dl)**n&
    &/DBLE(n_fac)
   END DO
   kernel=kernel*dx
    
   !Perform convolution
   DO j1=lbound(y,1),ubound(y,1)
    y_tmp(j1,:)=y(j1,:)*kernel(j1)
   END DO
   yi(i1,:)=SUM(y_tmp,1,mask)

  END DO !1st dimension
 ELSEIF(ndim.EQ.2) THEN
  !Convolute along 2nd dimension

  IF(SIZE(x,1).NE.SIZE(y,2)) THEN
   WRITE(*,*) 'mod_filter.filter_gauss_conv2D: incongruent dimensions&
   & x and y'
   STOP
  END IF

  DO i2=lbound(yi,2),ubound(yi,2)
   !Calculate deconvolution operator at these differences
   kernel=DBLE(0); n_fac=1
   DO n=0,order
    IF(n.GT.0) n_fac=n_fac*n
    kernel=kernel+filter_diff_gauss(x(i2)-x,-0.5/dl**2,2*n)*(-dl*dl)**n&
    &/DBLE(n_fac)
   END DO
   kernel=kernel*dx

   !Perform convolution
   DO j2=lbound(y,2),ubound(y,2)
    y_tmp(:,j2)=y(:,j2)*kernel(j2)
   END DO
   yi(:,i2)=SUM(y_tmp,2,mask)

  END DO !2nd dimension
 ELSE
  WRITE(*,*) 'mod_filter.gauss_conv2D: ndim must be 1 or 2'
  STOP
 END IF

 !Normalize
 WHERE(.NOT.mask)
  yi=y
 END WHERE

 END FUNCTION filter_gauss_deconv2D

!------------------------------------------------------------------------------

 FUNCTION filter_gauss_conv3D(x,y,ndim,dl,mask_in) RESULT(yi)
 !FILTER_GAUSS_CONV3D As filter_gauss_conv1d but now for a 3D seperable grid
 !
 !Syntax:
 ! yi=filter_gauss_conv3d(x,y,ndim,dl,mask)
 !Input:
 ! x: p/q/r double array with coordinates along the 1st dimension
 ! y:  pxqxr double array with values at coordinates in x1,x2
 ! ndim: integer with dimension along which convolution takes place
 ! dl: double with standard deviation
 ! mask: pxqxr logical array with points that should be included in calculation
 !Output
 ! yi: pxqxr double array with convoluted values at coordinates in x
 !NB1: for symmetric matrices grid spacing must be uniform
 !NB2: dl>= grid spacing in x

 REAL(8),INTENT(in)::x(:),y(:,:,:),dl
 INTEGER,INTENT(in)::ndim
 LOGICAL,INTENT(in),OPTIONAL::mask_in(lbound(y,1):ubound(y,1),&
 &lbound(y,2):ubound(y,2),lbound(y,3):ubound(y,3))
 LOGICAL::mask(lbound(y,1):ubound(y,1),&
 &lbound(y,2):ubound(y,2),lbound(y,3):ubound(y,3))
 REAL(8)::yi(lbound(y,1):ubound(y,1),lbound(y,2):ubound(y,2),&
 &lbound(y,3):ubound(y,3))
 REAL(8)::kernel(lbound(x,1):ubound(x,1)),&
 &dx(lbound(x,1):ubound(x,1)),&
 &y_tmp(lbound(y,1):ubound(y,1),lbound(y,2):ubound(y,2),&
 &lbound(y,3):ubound(y,3))
 INTEGER::i1,i2,i3,j1,j2,j3

 IF(PRESENT(mask_in)) THEN
  mask=mask_in
 ELSE
  mask=.TRUE.
 END IF

 !Calculate grid spacing
 DO i1=lbound(x,1)+1,ubound(x,1)-1
  dx(i1)=0.5*(x(i1+1)-x(i1-1))
 END DO
 dx(lbound(x,1))=dx(lbound(x,1)+1)
 dx(ubound(x,1))=dx(ubound(x,1)-1)
 
 yi=DBLE(0); y_tmp=DBLE(0)
  !write(*,*) 'conv'
  !write(*,*) lbound(yi,1),ubound(yi,1),lbound(yi,2),ubound(yi,2)
 IF(ndim.EQ.1) THEN
  !Convolute along 1st dimension

  IF(SIZE(x,1).NE.SIZE(y,1)) THEN
   WRITE(*,*) 'mod_filter.filter_gauss_conv2D: incongruent dimensions&
   & x,y'
   STOP
  END IF
 
  DO i1=lbound(yi,1),ubound(yi,1)
   !Calculate kernel
   kernel=EXP(-.5*(x(i1)-x)**2/dl**2)*dx
 
   !Perform convolution
   DO j1=lbound(y,1),ubound(y,1)
    y_tmp(j1,:,:)=y(j1,:,:)*kernel(j1)
   END DO
   yi(i1,:,:)=SUM(y_tmp,1,mask)

  END DO !1st dimension
 ELSEIF(ndim.EQ.2) THEN
  !Convolute along 2nd dimension

  IF(SIZE(x,1).NE.SIZE(y,2)) THEN
   WRITE(*,*) 'mod_filter.filter_gauss_conv2D: incongruent dimensions&
   & x and y'
   STOP
  END IF

  DO i2=lbound(yi,2),ubound(yi,2)
   !Calculate kernel
   kernel=EXP(-.5*(x(i2)-x)**2/dl**2)*dx

   !Perform convolution
   DO j2=lbound(y,2),ubound(y,2)
    y_tmp(:,j2,:)=y(:,j2,:)*kernel(j2)
   END DO
   yi(:,i2,:)=SUM(y_tmp,2,mask)

  END DO !2nd dimension
 ELSEIF(ndim.EQ.3) THEN
  !Convolute along 3rd dimension

  IF(SIZE(x,1).NE.SIZE(y,3)) THEN
   WRITE(*,*) 'mod_filter.filter_gauss_conv2D: incongruent dimensions&
   & x and y'
   STOP
  END IF

  DO i3=lbound(yi,3),ubound(yi,3)
   !Calculate kernel
   kernel=EXP(-.5*(x(i3)-x)**2/dl**2)*dx

   !Perform convolution
   DO j3=lbound(y,3),ubound(y,3)
    y_tmp(:,:,j3)=y(:,:,j3)*kernel(j3)
   END DO
   yi(:,:,i3)=SUM(y_tmp,3,mask)

  END DO !3rd dimension
 ELSE
  WRITE(*,*) 'mod_filter.gauss_conv2D: ndim must be 1 or 2'
  STOP
 END IF

 !Normalize
 yi=yi/SQRT(2*filter_pi*dl**2)
 WHERE(.NOT.mask)
  yi=y
 END WHERE

 END FUNCTION filter_gauss_conv3D

!----------------------------------------------------------------------------

 FUNCTION filter_gauss_deconv3D(x,y,ndim,dl,order,mask_in) RESULT(yi)
 !FILTER_GAUSS_DECONV3D As filter_gauss_deconv1d but now for a 3D grid
 !
 !Syntax:
 ! yi=filter_gauss_conv3d(x,y,ndim,dl,mask)
 !Input:
 ! x: p/q/r double array with coordinates along the 1st dimension
 ! y:  pxqxr double array with values at coordinates in x1,x2
 ! ndim: integer with dimension along which deconvolution takes place
 ! dl: double with standard deviations 
 ! mask: pxqxr logical array with points that should be included in calculation
 ! order: integer with the order of the deconvolution operator
 !Output
 ! yi: pxqxr double array with convoluted values at coordinates in x
 !NB1: for symmetric matrices grid spacing must be uniform
 !NB2: dl>= grid spacing in x

 INTEGER,INTENT(in)::order
 REAL(8),INTENT(in)::x(:),y(:,:,:),dl
 INTEGER,INTENT(in)::ndim
 LOGICAL,INTENT(in),OPTIONAL::mask_in(lbound(y,1):ubound(y,1),&
 &lbound(y,2):ubound(y,2),lbound(y,3):ubound(y,3))
 LOGICAL::mask(lbound(y,1):ubound(y,1),&
 &lbound(y,2):ubound(y,2),lbound(y,3):ubound(y,3))
 REAL(8)::yi(lbound(y,1):ubound(y,1),lbound(y,2):ubound(y,2),&
 &lbound(y,3):ubound(y,3) )
 REAL(8)::kernel(lbound(x,1):ubound(x,1)),&
 &dx(lbound(x,1):ubound(x,1)),&
 &y_tmp(lbound(y,1):ubound(y,1),lbound(y,2):ubound(y,2),&
 &lbound(y,3):ubound(y,3))
 INTEGER::i1,i2,i3,j1,j2,j3,n,n_fac

 IF(PRESENT(mask_in)) THEN
  mask=mask_in
 ELSE
  mask=.TRUE.
 END IF

 !Calculate grid spacing
 DO i1=lbound(x,1)+1,ubound(x,1)-1
  dx(i1)=0.5*(x(i1+1)-x(i1-1))
 END DO
 dx(lbound(x,1))=dx(lbound(x,1)+1)
 dx(ubound(x,1))=dx(ubound(x,1)-1)
 dx=dx/SQRT(2*filter_pi*dl**2)
 
 yi=DBLE(0); y_tmp=DBLE(0)
  !write(*,*) 'conv'
  !write(*,*) lbound(yi,1),ubound(yi,1),lbound(yi,2),ubound(yi,2)
 IF(ndim.EQ.1) THEN
  !Convolute along 1st dimension

  IF(SIZE(x,1).NE.SIZE(y,1)) THEN
   WRITE(*,*) 'mod_filter.filter_gauss_conv2D: incongruent &
   &dimensions x,y'
   STOP
  END IF
 
  DO i1=lbound(yi,1),ubound(yi,1)

   !Calculate deconvolution operator at these differences
   kernel=DBLE(0); n_fac=1
   DO n=0,order
    IF(n.GT.0) n_fac=n_fac*n
    kernel=kernel+filter_diff_gauss(x(i1)-x,-0.5/dl**2,2*n)*(-dl*dl)**n&
    &/DBLE(n_fac)
   END DO
   kernel=kernel*dx
 
   !Perform convolution
   DO j1=lbound(y,1),ubound(y,1)
    y_tmp(j1,:,:)=y(j1,:,:)*kernel(j1)
   END DO
   yi(i1,:,:)=SUM(y_tmp,1,mask)

  END DO !1st dimension
 ELSEIF(ndim.EQ.2) THEN
  !Convolute along 2nd dimension

  IF(SIZE(x,1).NE.SIZE(y,2)) THEN
   WRITE(*,*) 'mod_filter.filter_gauss_conv2D: incongruent dimensions&
   & x and y'
   STOP
  END IF

  DO i2=lbound(yi,2),ubound(yi,2)
   !Calculate deconvolution operator at these differences
   kernel=DBLE(0); n_fac=1
   DO n=0,order
    IF(n.GT.0) n_fac=n_fac*n
    kernel=kernel+filter_diff_gauss(x(i2)-x,-0.5/dl**2,2*n)*(-dl*dl)**n&
    &/DBLE(n_fac)
   END DO
   kernel=kernel*dx

   !Perform convolution
   DO j2=lbound(y,2),ubound(y,2)
    y_tmp(:,j2,:)=y(:,j2,:)*kernel(j2)
   END DO
   yi(:,i2,:)=SUM(y_tmp,2,mask)

  END DO !2nd dimension
 ELSEIF(ndim.EQ.3) THEN
  !Convolute along 3rd dimension

  IF(SIZE(x,1).NE.SIZE(y,3)) THEN
   WRITE(*,*) 'mod_filter.filter_gauss_conv2D: incongruent dimensions&
    x and y'
   STOP
  END IF

  DO i3=lbound(yi,3),ubound(yi,3)
   !Calculate deconvolution operator at these differences
   kernel=DBLE(0); n_fac=1
   DO n=0,order
    IF(n.GT.0) n_fac=n_fac*n
    kernel=kernel+filter_diff_gauss(x(i3)-x,-0.5/dl**2,2*n)*(-dl*dl)**n&
    &/DBLE(n_fac)
   END DO
   kernel=kernel*dx

   !Perform convolution
   DO j3=lbound(y,3),ubound(y,3)
    y_tmp(:,:,j3)=y(:,:,j3)*kernel(j3)
   END DO
   yi(:,:,i3)=SUM(y_tmp,3,mask)

  END DO !3rd dimension
 ELSE
  WRITE(*,*) 'mod_filter.gauss_conv2D: ndim must be 1 or 2'
  STOP
 END IF

 !Normalize
 WHERE(.NOT.mask)
  yi=y
 END WHERE

 END FUNCTION filter_gauss_deconv3D

!----------------------------------------------------------------------------

 FUNCTION filter_diff_gauss(x,alpha,order) RESULT(y)
 !FILTER_DIFF_GAUSS Calculates the nth order derivative of 
 !x->exp(alpha*x^2)
 !
 !Syntax:
 ! y=filter_diff_gauss(x,sigma,order)
 !Input:
 ! x: [1D double array] x-coordinates
 ! sigma: [double] sigma as described above
 ! order: the order of the derivative
 !Output:
 ! y: [1D double array] values of derivative evaluated at x

 REAL(8),INTENT(in)::x(:),alpha
 INTEGER,INTENT(in)::order
 REAL(8)::y(lbound(x,1):ubound(x,1))
 REAL(8)::z_tmp(lbound(x,1):ubound(x,1)),z(lbound(x,1):ubound(x,1))
 INTEGER::n_fac,i1_fac,i2_fac
 INTEGER::i1,i2

 y=DBLE(0)

 !Calculate n!
 n_fac=1
 DO i1=1,order
  n_fac=n_fac*i1
 END DO

 i1_fac=1
 !Calculate the sum of x monomials
 DO i1=0,FLOOR(0.5*DBLE(order))
  !Calculate i1_fac
  IF(i1.GT.0) i1_fac=i1_fac*i1

  !Calculate (-2i1+n)!
  i2_fac=1
  DO i2=1,(-2*i1+order)
   i2_fac=i2_fac*i2
  END DO
 
  y=y+(DBLE(2**(order-2*i1))*(alpha)**(order-i1)*&
  &x**(order-2*i1))/DBLE(i1_fac)/DBLE(i2_fac)

 END DO

 !Multiply with exp(-alpha*x^2)*n!
 y=y*EXP(alpha*x**2)*DBLE(n_fac)

 !Finite difference differentiation (not active)
 IF(.FALSE.) THEN
 z_tmp=EXP(alpha*x**2); z=z_tmp
 DO i1=1,order
  DO i2=lbound(z,1)+1,ubound(z,1)-1
   z(i2)=0.5*(z_tmp(i2)-z_tmp(i2-1))/(x(i2)-x(i2-1))+&
   &0.5*(z_tmp(i2+1)-z_tmp(i2))/(x(i2+1)-x(i2))
   END DO
  z(1)=z(2); z(ubound(z,1))=z(ubound(z,1)-1)
  z_tmp=z
 END DO
 END IF!finite difference

 END FUNCTION filter_diff_gauss

!-----------------------------------------------------------------------------
 FUNCTION filter_trapz(x,y,mask) RESULT(int)
 !FILTER_TRAPZ Numerical 1D-integration using trapezium rule

 REAL(8),INTENT(in)::x(:),y(:)
 LOGICAL,INTENT(in),OPTIONAL::mask(:)
 LOGICAL::mask_loc(lbound(x,1):ubound(x,1))
 REAL(8)::y_tmp(lbound(y,1):ubound(y,1))
 REAL(8)::int

 y_tmp=y

 IF(PRESENT(mask)) THEN 
  int=SUM( (0.5*y_tmp(lbound(y,1):ubound(y,1)-1)+&
  &0.5*y_tmp(lbound(y,1)+1:ubound(y,1)))*&
  &(x(lbound(x,1)+1:ubound(x,1))-x(lbound(x,1):ubound(x,1)-1)),&
  &mask(lbound(x,1):ubound(x,1)-1).AND.mask(lbound(x,1)+1:ubound(x,1)) )
 ELSE
  int=SUM( (0.5*y_tmp(lbound(y,1):ubound(y,1)-1)+&
  &0.5*y_tmp(lbound(y,1)+1:ubound(y,1)))*&
  &(x(lbound(x,1)+1:ubound(x,1))-x(lbound(x,1):ubound(x,1)-1)) )
 END IF

 END FUNCTION filter_trapz

!--------------------------------------------------------------------------
END MODULE
