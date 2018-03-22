MODULE mod_interp
!
!Contains:
! get_sigma_param(grd_ncid,ini_ncid,param,h)
! interp1d(x_in,val_in,x_out)     
! mesh_interp2d(x_in,y_in,mask_in,val_in,x_out,y_out,mask_out)   
! mesh_interp3d(x_in,y_in,mask_in,z_in,val_in,x_out,y_out,mask_out,z_out)  
! mesh_extrap2d(x_in,y_in,mask_in,val_in,opt)    
! mesh_extrap2d(x_in,y_in,mask_in,val_in,opt) 
! z2sigma(param,h,zeta,z)   
! sigma2z(param,h,zeta,sigma)                                                  ! zweight(param,h,zeta,var)                                              
 
 USE mod_netcdf
 USE omp_lib
 IMPLICIT NONE

 !INCLUDE 'netcdf.inc'

 TYPE sigma_param
  REAL(8)::theta_s=0.65,theta_b=0.58,t_cline=1e16,alpha=0.0,beta=0.0
  INTEGER::v_transform=1,v_stretching=2
  INTEGER::n_sigma
  REAL(8)::gamma=3.0
  CHARACTER(len=1024)::file
 END TYPE

 INTERFACE interp
  MODULE PROCEDURE interp1d, interp2d, interp3d
 END INTERFACE

 INTERFACE sigma2z
  MODULE PROCEDURE sigma2z_scalar, sigma2z_2d
 END INTERFACE

 INTERFACE z2sigma
  MODULE PROCEDURE z2sigma_scalar, z2sigma_2d
 END INTERFACE

 INTERFACE mesh2mesh_interp
  MODULE PROCEDURE mesh2mesh_interp2d,mesh2mesh_interp3d,&
  &mesh2mesh_interp4d
 END INTERFACE

 INTERFACE mesh2mesh_extrap
  MODULE PROCEDURE mesh2mesh_extrap2d,mesh2mesh_extrap3d,&
  &mesh2mesh_extrap4d
 END INTERFACE

 INTERFACE mesh2mesh_tl_interp
  MODULE PROCEDURE mesh2mesh_tl_interp2d,mesh2mesh_tl_interp3d,&
  &mesh2mesh_tl_interp4d
 END INTERFACE


 INTERFACE mesh2mesh_ad_interp
  MODULE PROCEDURE mesh2mesh_ad_interp2d,mesh2mesh_ad_interp3d,&
  &mesh2mesh_ad_interp4d
 END INTERFACE


 REAL:: interp_fill_single
 !REAL(8)::interp_fill_value=DBLE(HUGE(interp_fill_single))
 REAL(8)::interp_fill_value=DBLE(1e37)
!-------------------------------------------------------------------------------

CONTAINS

 SUBROUTINE get_sigma_param(grd_ncid,ini_ncid,param,h)
 !GET_SIGMA_PARAM Loads bottom and sigma parameters from netcdf files
 !
 !Syntax:
 ! CALL get_sigma_param(grd_file,ini_file,param,h)
 !Input:
 ! grd_file: netcdf file with ROMS grid file
 ! ini_file: netcdf file with model output
 !Output:
 ! param: TYPE(sigma_param) struct with parameters for sigma transformation
 ! h: 2D double array with depth bottom .w.r.t reference level

 USE mod_netcdf
 IMPLICIT NONE

 INTEGER,INTENT(in)::grd_ncid
 INTEGER,INTENT(in)::ini_ncid
 TYPE(sigma_param),INTENT(out)::param
 REAL(8),ALLOCATABLE,INTENT(out),optional::h(:,:)
 INTEGER::status,ncid,varid,dimids(2),diml(2),i1,diml1

 !Read bottom file
 IF(PRESENT(h)) THEN
  status=nf_inq_varid(grd_ncid,'h',varid)
  IF(status.NE.nf_noerr) THEN
   WRITE(*,*) 'mod_interp.get_sigma_param: cannot read bottom.'
   STOP
  END IF
  status=nf_inq_vardimid(grd_ncid,varid,dimids)
  status=nf_inq_dimlen(grd_ncid,dimids(1),diml(1))
  status=nf_inq_dimlen(grd_ncid,dimids(2),diml(2))
  ALLOCATE(h(diml(1),diml(2))); h=0.0
  status=nf_get_vara_double(grd_ncid,varid,[1,1],&
  &diml,h)
 END IF

 !Read sigma parameters
 status=nf_inq_varid(ini_ncid,'theta_s',varid)
 IF(status.NE.nf_noerr) THEN
  WRITE(*,*) 'mod_interp.get_sigma_param: cannot find theta_s.'
  STOP
 END IF
 status=nf_get_vara_double(ini_ncid,varid,1,1,param%theta_s)

 status=nf_inq_varid(ini_ncid,'theta_b',varid)
 IF(status.NE.nf_noerr) THEN
  WRITE(*,*) 'mod_interp.get_sigma_param: cannot find theta_b.'
  STOP
 END IF
 status=nf_get_vara_double(ini_ncid,varid,1,1,param%theta_b)

 status=nf_inq_varid(ini_ncid,'Tcline',varid)
 IF(status.NE.nf_noerr) THEN
  WRITE(*,*) 'mod_interp.get_sigma_param: cannot find Tcline.'
  STOP
 END IF
 status=nf_get_vara_double(ini_ncid,varid,1,1,param%t_cline)

 status=nf_inq_varid(ini_ncid,'Vtransform',varid)
 IF(status.NE.nf_noerr) THEN
  WRITE(*,*) 'mod_interp.get_sigma_param: cannot find Vtransform.'
  STOP
 END IF
 status=nf_get_vara_int(ini_ncid,varid,(/1/),(/1/),param%v_transform)

 status=nf_inq_varid(ini_ncid,'Vstretching',varid)
 IF(status.NE.nf_noerr) THEN
  WRITE(*,*) 'mod_interp.get_sigma_param: cannot find Vstretching.'
  STOP
 END IF
 status=nf_get_vara_int(ini_ncid,varid,1,1,param%v_stretching) 

 !Create the sigma values of the layers
 status=nf_inq_dimid(ini_ncid,'s_rho',varid)
 IF(status.NE.nf_noerr) THEN
  status=nf_inq_dimid(ini_ncid,'N',varid)
 END IF
 IF(status.EQ.nf_noerr) THEN
  status=nf_inq_dimlen(ini_ncid,varid,param%n_sigma)
 ELSE
  param%n_sigma=-999
 END IF

 END SUBROUTINE get_sigma_param

!----------------------------------------------------------------------------
SUBROUTINE extrapz(val)
 IMPLICIT NONE
 REAL(8),intent(inout)::val(:)
 INTEGER::i1
 
 DO i1=lbound(val,1),ubound(val,1)
   IF(val(i1).NE.interp_fill_value) THEN
    val(lbound(val,1):i1)=val(i1)
    EXIT
   END IF
  END DO

  DO i1=ubound(val,1),lbound(val,1),-1
   IF(val(i1).NE.interp_fill_value) THEN
    val(i1:ubound(val,1))=val(i1)
    EXIT
   END IF
  END DO

END SUBROUTINE


!-----------------------------------------------------------------------------

FUNCTION extrap1d(x_in,val_in,mask_in) RESULT(val_out)
!EXTRAP1D Extrapolates to fill NaN values
!
!Syntax:
! val_out=extrap(x_in,val_in,mask_in,opt)
!Input:
! x_in: [1D-double] coordinates grid
! val_in: [1D-double] value at coordinates
! mask_in: [1D-logical] mask
! opt: 'maskfill' extrapolate to all NaN values
!      'gapfill' extrapolate to only NaN values in mask
!Output:
! val_out: [1D-double] values without NaNs

 REAL(8),INTENT(in)::x_in(:),val_in(:)
 LOGICAL,INTENT(in),OPTIONAL::mask_in(:)
 REAL(8)::val_out(lbound(val_in,1):ubound(val_in,1))
 LOGICAL::mask_out(lbound(val_in,1):ubound(val_in,1))
 LOGICAL::mask_val(lbound(val_in,1):ubound(val_in,1))
 INTEGER::i1,iClose


 !Indicate where extrapolation is necessary
 mask_out=val_in.EQ.interp_fill_value
 mask_val=val_in.NE.interp_fill_value
 IF(PRESENT(mask_in)) THEN
  mask_val=mask_val.AND.mask_in
 END IF

 IF(ANY(mask_val)) THEN
  !Extrapolate
  val_out=val_in
  DO i1=lbound(x_in,1),ubound(x_in,1)
   IF(.NOT.mask_out(i1)) CYCLE
   iClose=MINLOC(ABS(x_in(i1)-x_in),1,mask_val)
   val_out(i1)=val_in(iClose) 
  END DO
 ELSE
  val_out=val_in
 END IF

END FUNCTION extrap1d

!-----------------------------------------------------------------------------
FUNCTION extrap1d_lin(x_in,val_in,dx,mask_in) RESULT(val_out)

IMPLICIT NONE 
 
 REAL(8),INTENT(in)::dx
 REAL(8),INTENT(in)::x_in(:),val_in(:)
 LOGICAL,INTENT(in),OPTIONAL::mask_in(:)
 REAL(8)::val_out(lbound(val_in,1):ubound(val_in,1))
 LOGICAL::mask_out(lbound(val_in,1):ubound(val_in,1))
 LOGICAL::mask_val(lbound(val_in,1):ubound(val_in,1))
 INTEGER::i1,i2,j1,j2,j_tmp(2),nj
 REAL(8)::xbnd,coef

 !Indicate where extrapolation is necessary
 mask_out=val_in.EQ.interp_fill_value
 mask_val=val_in.NE.interp_fill_value
 IF(PRESENT(mask_in)) THEN
  mask_val=mask_val.AND.mask_in
 END IF
 
 IF(COUNT(mask_val).LT.2) THEN
  val_out=interp_fill_value
 ELSE
  val_out=val_in

  !Extrapolate 
  i1=MINLOC(x_in,1,mask_val)
  i2=MAXLOC(x_in,1,mask_val)  

  j1=i1
  xbnd=MAXVAL([x_in(j1)+dx,MINVAL(x_in,x_in.GT.x_in(j1))])
  j2=MAXLOC(x_in,1,x_in.LE.xbnd)
  IF(j2.LT.j1) THEN 
   j_tmp=[j1,j2]; j2=j_tmp(1); j1=j_tmp(2)
  END IF
  nj=dble(j2-j1+1)
 
  !regression 
  coef=(nj*SUM(x_in(j1:j2)*val_in(j1:j2))-&
  &sum(x_in(j1:j2))*sum(val_in(j1:j2)))/&
  &(nj*sum(x_in(j1:j2)**2)-&
  &sum(x_in(j1:j2))**2)

  !extrapolation
  IF(i1.LE.i2) THEN  
   val_out(lbound(val_out,1):j1)=val_in(j1)+&
   &(x_in(lbound(x_in,1):j1)-x_in(j1))*coef
  ELSE
   val_out(j2:ubound(val_out,1))=val_in(j2)+&
   &(x_in(j2:ubound(x_in,1))-x_in(j2))*coef
  END IF

  j2=i2
  xbnd=MINVAL([x_in(j2)-dx,MAXVAL(x_in,x_in.LT.x_in(j2))])
  j1=MINLOC(x_in,1,x_in.GE.xbnd)
  IF(j2.LT.j1) THEN
   j_tmp=[j1,j2]; j2=j_tmp(1); j1=j_tmp(2)
  END IF
  nj=dble(j2-j1+1)

  !regression 
  coef=(nj*SUM(x_in(j1:j2)*val_in(j1:j2))-&
  &sum(x_in(j1:j2))*sum(val_in(j1:j2)))/&
  &(nj*sum(x_in(j1:j2)**2)-&
  &sum(x_in(j1:j2))**2)

  !extrapolation
  IF(i1.LE.i2) THEN
   val_out(j2:size(x_in))=val_in(j2)+&
   &(x_in(j2:ubound(val_out,1))-x_in(j2))*coef
  ELSE
   val_out(lbound(val_out,1):j1)=val_in(j1)+&
   &(x_in(lbound(x_in,1):j1)-x_in(j1))*coef
  END IF
 END IF

END FUNCTION extrap1d_lin


!-----------------------------------------------------------------------------

FUNCTION extrap2d(x_in,y_in,val_in,mask_in) RESULT(val_out)
!EXTRAP1D Extrapolates to fill NaN values
!
!Syntax:
! val_out=extrap(x_in,val_in,mask_in,opt)
!Input:
! x_in: [2D-double] coordinates grid 1st dimension
! y_in: [2D-double] coordinates grid 2nd dimension
! val_in: [2D-double] value at coordinates
! mask_in: [2D-logical] mask
!Output:
! val_out: [2D-double] values without NaNs

REAL(8),INTENT(in)::x_in(:,:),y_in(:,:),val_in(:,:)
REAL(8)::val_out(lbound(val_in,1):ubound(val_in,1),&
&lbound(val_in,2):ubound(val_in,2))
LOGICAL,INTENT(in),OPTIONAL::mask_in(:,:)
LOGICAL::mask_out(lbound(val_in,1):ubound(val_in,1),&
&lbound(val_in,2):ubound(val_in,2))
LOGICAL::mask_val(lbound(val_in,1):ubound(val_in,1),&
&lbound(val_in,2):ubound(val_in,2))
REAL(8)::w(lbound(val_in,1):ubound(val_in,1),&
&lbound(val_in,2):ubound(val_in,2))
REAL(8)::w1
INTEGER::i0,i1,i2,iList(size(val_in),2)

!Check dimensions input
IF(ANY(shape(val_in).NE.shape(x_in))) THEN
 WRITE(*,*) 'mod_interp.extrap2d: dimensions val_in and x_in&
 & incongruent'
 STOP
END IF
IF(ANY(shape(val_in).NE.shape(y_in))) THEN
 WRITE(*,*) 'mod_interp.extrap2d: dimensions val_in and y_in&
 & incongruent'
 STOP
END IF

 !Indicate where extrapolation is necessary
 mask_out=val_in.EQ.interp_fill_value
 mask_val=val_in.NE.interp_fill_value
 IF(PRESENT(mask_in)) THEN
  mask_val=mask_in.AND.mask_val
 END IF

 IF(.NOT.ANY(mask_val)) THEN
  val_out=val_in
 ELSE
  !Initialize output
  val_out=val_in
  WHERE(mask_out)
   val_out=dble(0)
  END WHERE
  w=dble(0)

  !Make list with indices om mask_out
  i0=1
  DO i2=lbound(y_in,2),ubound(y_in,2)
  DO i1=lbound(x_in,1),ubound(x_in,1)
   iList(i0,1)=i1; iList(i0,2)=i2;
   i0=i0+1
  END DO
  END DO

  
  DO WHILE(ANY(mask_out))

   DO i0=1,size(iList,1)
    i1=iList(i0,1); i2=iList(i0,2)
    IF(.NOT.mask_out(i1,i2)) CYCLE
   
    !left
    IF(i1-1.GE.lbound(val_in,1)) THEN
    IF(mask_val(i1-1,i2)) THEN
     w1=dble(1)/(x_in(i1-1,i2)-x_in(i1,i2))**2
     w(i1,i2)=w(i1,i2)+w1
     val_out(i1,i2)=val_out(i1,i2)+val_out(i1-1,i2)*w1
    END IF
    END IF  

    !right
    IF(i1+1.LE.ubound(val_in,1)) THEN
    IF(mask_val(i1+1,i2)) THEN
     w1=dble(1)/(x_in(i1+1,i2)-x_in(i1,i2))**2
     w(i1,i2)=w(i1,i2)+w1
     val_out(i1,i2)=val_out(i1,i2)+val_out(i1+1,i2)*w1
    END IF 
    END IF

    !top
    IF(i2+1.LE.ubound(val_in,2)) THEN
    IF(mask_val(i1,i2+1)) THEN
     w1=dble(1)/(y_in(i1,i2+1)-y_in(i1,i2))**2
     w(i1,i2)=w(i1,i2)+w1
     val_out(i1,i2)=val_out(i1,i2)+val_out(i1,i2+1)*w1
    END IF 
    END IF

    !bottom
    IF(i2-1.GE.lbound(val_in,2)) THEN
    IF(mask_val(i1,i2-1)) THEN
     w1=dble(1)/(y_in(i1,i2-1)-y_in(i1,i2))**2
     w(i1,i2)=w(i1,i2)+w1
     val_out(i1,i2)=val_out(i1,i2)+val_out(i1,i2-1)*w1
    END IF 
    END IF
 
   END DO

   !update masks
   WHERE(w.NE.dble(0))
    val_out=val_out/w
    mask_val=.TRUE.
    mask_out=.FALSE.
   END WHERE
   w=dble(0)

  END DO
 END IF


END FUNCTION extrap2d

!------------------------------------------------------------------------------

FUNCTION extrap3d(x_in,y_in,z_in,val_in,mask_in) RESULT(val_out)
!EXTRAP3D Extrapolates to fill NaN values
!
!Syntax:
! val_out=extrap(x_in,val_in,mask_in,opt)
!Input:
! x_in: [3D-double] coordinates grid 1st dimension
! y_in: [3D-double] coordinates grid 2nd dimension
! z_in: [3D-double] coordinates grid 3rd dimension
! val_in: [3D-double] value at coordinates
! mask_in: [3D-logical] mask
!Output:
! val_out: [3D-double] values without NaNs

REAL(8),INTENT(in)::x_in(:,:,:),y_in(:,:,:),z_in(:,:,:),val_in(:,:,:)
REAL(8)::val_out(lbound(val_in,1):ubound(val_in,1),&
&lbound(val_in,2):ubound(val_in,2),lbound(val_in,3):ubound(val_in,3))
LOGICAL,INTENT(in),OPTIONAL::mask_in(:,:,:)
LOGICAL::mask_out(lbound(mask_in,1):ubound(mask_in,1),&
&lbound(mask_in,2):ubound(mask_in,2),&
&lbound(mask_in,3):ubound(mask_in,3))
LOGICAL::mask_val(lbound(mask_in,1):ubound(mask_in,1),&
&lbound(mask_in,2):ubound(mask_in,2),&
&lbound(mask_in,3):ubound(mask_in,3))
REAL(8)::w(lbound(mask_in,1):ubound(mask_in,1),&
&lbound(mask_in,2):ubound(mask_in,2),&
&lbound(mask_in,3):ubound(mask_in,3))
REAL(8)::w1
INTEGER::i0,i1,i2,i3,iList(size(mask_in),3)

!Check dimensions input
IF(ANY(shape(val_in).NE.shape(x_in))) THEN
 WRITE(*,*) 'mod_interp.extrap3d: dimensions val_in and x_in&
 & incongruent'
 STOP
END IF
IF(ANY(shape(val_in).NE.shape(y_in))) THEN
 WRITE(*,*) 'mod_interp.extrap3d: dimensions val_in and y_in&
 & incongruent'
 STOP
END IF
IF(ANY(shape(val_in).NE.shape(z_in))) THEN
 WRITE(*,*) 'mod_interp.extrap3d: dimensions val_in and z_in&
 & incongruent'
 STOP
END IF

 !Indicate where extrapolation is necessary
 mask_out=val_in.EQ.interp_fill_value
 mask_val=val_in.NE.interp_fill_value
 IF(PRESENT(mask_in)) THEN
  mask_val=mask_val.AND.mask_in
 END IF

 !Initialize output
 val_out=val_in
 WHERE(mask_out)
  val_out=dble(0)
 END WHERE
 w=dble(0)

 !Make list with indices om mask_out
 i0=1
 DO i3=lbound(z_in,3),ubound(z_in,3)
 DO i2=lbound(y_in,2),ubound(y_in,2)
 DO i1=lbound(x_in,1),ubound(x_in,1)
  iList(i0,1)=i1; iList(i0,2)=i2; iList(i0,3)=i3
  i0=i0+1
 END DO
 END DO
 END DO

 IF(.NOT.ANY(mask_val)) THEN
  val_out=val_in
 ELSE

  DO WHILE(ANY(mask_out))

   DO i0=1,size(iList,1)
    i1=iList(i0,1); i2=iList(i0,2); i3=iList(i0,3)
    IF(.NOT.mask_out(i1,i2,i3)) CYCLE
   
    !left
    IF(i1-1.GE.lbound(val_in,1)) THEN
    IF(mask_val(i1-1,i2,i3)) THEN
     w1=dble(1)/(x_in(i1-1,i2,i3)-x_in(i1,i2,i3))**2
     w(i1,i2,i3)=w(i1,i2,i3)+w1
     val_out(i1,i2,i3)=val_out(i1,i2,i3)+val_out(i1-1,i2,i3)*w1
    END IF
    END IF  

    !right
    IF(i1+1.LE.ubound(val_in,1)) THEN
    IF(mask_val(i1+1,i2,i3)) THEN
     w1=dble(1)/(x_in(i1+1,i2,i3)-x_in(i1,i2,i3))**2
     w(i1,i2,i3)=w(i1,i2,i3)+w1
     val_out(i1,i2,i3)=val_out(i1,i2,i3)+val_out(i1+1,i2,i3)*w1
    END IF 
    END IF

    !top
    IF(i2+1.LE.ubound(val_in,2)) THEN
    IF(mask_val(i1,i2+1,i3)) THEN
     w1=dble(1)/(y_in(i1,i2+1,i3)-y_in(i1,i2,i3))**2
     w(i1,i2,i3)=w(i1,i2,i3)+w1
     val_out(i1,i2,i3)=val_out(i1,i2,i3)+val_out(i1,i2+1,i3)*w1
    END IF 
    END IF

    !bottom
    IF(i2-1.GE.lbound(val_in,2)) THEN
    IF(mask_val(i1,i2-1,i3)) THEN
     w1=dble(1)/(y_in(i1,i2-1,i3)-y_in(i1,i2,i3))**2
     w(i1,i2,i3)=w(i1,i2,i3)+w1
     val_out(i1,i2,i3)=val_out(i1,i2,i3)+val_out(i1,i2-1,i3)*w1
    END IF 
    END IF

    !below
    IF(i3-1.GE.lbound(val_in,3)) THEN
    IF(mask_val(i1,i2,i3-1)) THEN
     w1=dble(1)/(z_in(i1,i2,i3-1)-z_in(i1,i2,i3))**2
     w(i1,i2,i3)=w(i1,i2,i3)+w1
     val_out(i1,i2,i3)=val_out(i1,i2,i3)+val_out(i1,i2,i3-1)*w1
    END IF 
    END IF

    !above
    IF(i3+1.GE.ubound(val_in,3)) THEN
    IF(mask_val(i1,i2,i3+1)) THEN
     w1=dble(1)/(z_in(i1,i2,i3+1)-z_in(i1,i2,i3))**2
     w(i1,i2,i3)=w(i1,i2,i3)+w1
     val_out(i1,i2,i3)=val_out(i1,i2,i3)+val_out(i1,i2,i3+1)*w1
    END IF 
    END IF
 
   END DO

   !update masks
   WHERE(w.NE.dble(0))
    val_out=val_out/w
    mask_val=.TRUE.
    mask_out=.FALSE.
   END WHERE
   w=dble(0)

  END DO
 END IF

END FUNCTION extrap3d

!------------------------------------------------------------------------------

 FUNCTION interp1d(x_in,val_in,x_out)  RESULT(val_out) 
 !INTERP1D Performs 1D interpolation. For x_out outside domain x_in boundary values
 !are used
 !
 !Syntax: 
 ! val_out=interp(x_in,val_in,x_out)
 !Input:
 ! x_in: 1D double array with x-coordinates where values are given. Values must be
 !       monotomic increasing and grid spacing may vary. 
 ! x_out: 1D double array with x-coordinates where values have to be interpolated.
 ! val_in: 1D double array with the values given
 !Output:
 ! val_out: 1D double array with interpolated values.

 IMPLICIT NONE

 REAL(8),DIMENSION(:),INTENT(in)::x_in,val_in,x_out
 REAL(8),DIMENSION(:)::&
 &val_out(lbound(x_out,1):ubound(x_out,1))
 INTEGER::i1,ibnd(2)
 REAL(8)::bnd(2),w

  !boundaries interval
  bnd=[MINVAL(x_in),MAXVAL(x_in)]
 
  DO i1=lbound(val_out,1),ubound(val_out,1)
   IF(x_out(i1).LT.bnd(1).OR.x_out(i1).GT.bnd(2)) THEN
    val_out(i1)=interp_fill_value
   ELSE
    ibnd(1)=MAXLOC(x_in,1,x_in.LE.x_out(i1))
    ibnd(2)=MINLOC(x_in,1,x_in.GE.x_out(i1))
    IF(ibnd(1).EQ.ibnd(2)) THEN
     w=dble(0.5)
    ELSE
     w=(x_out(i1)-x_in(ibnd(1)))/(x_in(ibnd(2))-x_in(ibnd(1)))
    END IF    
    val_out(i1)=w*val_in(ibnd(2))+(dble(1)-w)*val_in(ibnd(1))
   END IF
  END DO

  WHERE(val_out.GE.0.99*interp_fill_value)
   val_out=interp_fill_value
  END WHERE

 END FUNCTION interp1d

!-----------------------------------------------------------------------------

FUNCTION interp2d(x_in,val_in,x_out,dim) RESULT(val_out)
!INTERP3D Interpolate one 2D grid to another with the 2 grid differing only
!along dimension dim
!
!Syntax:
! val_out=interp(x_in,val_in,x_out,dim)
!Input:
! x_in: [2D-double] input grid coordinate corresponding to dimension dim
! val_in: [2D-double] values on input grid
! x_out: [2D-double] output grid coordinate corresponding to dimension dim
! dim: [integer] dimension along which interpolating will take place
!Output:
! val_out: values on output grid
  
 IMPLICIT NONE
 
 REAL(8),INTENT(in),DIMENSION(:,:)::x_in,val_in,x_out
 INTEGER,INTENT(in)::dim
 REAL(8)::val_out(lbound(x_out,1):ubound(x_out,1),&
 &lbound(x_out,2):ubound(x_out,2))
 INTEGER::i0,i1,i2,shape_in(2),shape_out(2)

 !Test input
 shape_in=SHAPE(x_in)
 shape_out=SHAPE(x_out)
 IF( (dim.EQ.1.AND.shape_in(2).NE.shape_out(2)).OR.&
 & (dim.EQ.2.AND.shape_in(1).NE.shape_out(1)) ) THEN
  WRITE(*,*) 'mod_interp.interp3d: at least 2 dimensions of x_in and&
  &x_out must be equal.'
  STOP
 END IF 

 !Select case based on dimension 
 SELECT CASE(dim)
  CASE(1)
   !$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i0,i1,i2)
   DO i0=0,SIZE(x_out)/SIZE(x_out,dim)-1
    i2=lbound(x_out,2)+i0
    val_out(:,i2)=&
    &interp1d(x_in(:,i2),val_in(:,i2),&
    &x_out(:,i2))
   END DO
   !$OMP END PARALLEL DO
  CASE(2)
   !$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i0,i1,i2)
   DO i0=0,SIZE(x_out)/SIZE(x_out,dim)-1
    i1=lbound(x_out,1)+i0
    val_out(i1,:)=&
    &interp1d(x_in(i1,:),val_in(i1,:),&
    &x_out(i1,:))
   END DO
   !$OMP END PARALLEL DO
  CASE DEFAULT
   WRITE(*,*) "mod_interp.interp2d: dim must be 1 or 2"
   STOP
 END SELECT
  
 END FUNCTION interp2d

!-----------------------------------------------------------------------------

 FUNCTION interp3d(x_in,val_in,x_out,dim) RESULT(val_out)
 !INTERP3D Interpolate one 3D grid to another with the 2 grid differing only
 !along dimension dim
 !
 !Syntax:
 ! val_out=interp(x_in,val_in,x_out,dim)
 !Input:
 ! x_in: [3D-double] input grid coordinate corresponding to dimension dim
 ! val_in: [3D-doublevalues on input grid
 ! x_out: [3D-double] output grid coordinate corresponding to dimension dim
 ! dim: [integer] dimension along which interpolating will take place
 !Output:
 ! val_out: [3D-double] values on output grid
  
 IMPLICIT NONE
 
 REAL(8),INTENT(in),DIMENSION(:,:,:)::x_in,val_in,x_out
 INTEGER,INTENT(in)::dim
 REAL(8)::val_out(lbound(x_out,1):ubound(x_out,1),&
 &lbound(x_out,2):ubound(x_out,2),lbound(x_out,3):ubound(x_out,3))
 INTEGER::i0,i1,i2,i3,shape_in(3),shape_out(3)

 !Test input
 shape_in=SHAPE(x_in)
 shape_out=SHAPE(x_out)
 IF( (dim.EQ.1.AND.ANY(shape_in(2:3).NE.shape_out(2:3))).OR.&
 &(dim.EQ.2.AND.ANY(shape_in([1,3]).NE.shape_out([1,3]))).OR.&
 &(dim.EQ.3.AND.ANY(shape_in(1:2).NE.shape_out(1:2))) ) THEN
  WRITE(*,*) 'mod_interp.interp3d: at least 2 dimensions of x_in and&
  &x_out must be equal.'
  STOP
 END IF

 !Select case based on dimension 
 SELECT CASE(dim)
  CASE(1)
   !$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i0,i1,i2,i3)
   DO i0=0,SIZE(x_out)/SIZE(x_out,dim)-1
    i2=lbound(x_out,2)+MOD(i0,SIZE(x_out,2))
    i3=lbound(x_out,3)+FLOOR(dble(i0)/dble(SIZE(x_out,3)))
    val_out(:,i2,i3)=&
    &interp1d(x_in(:,i2,i3),val_in(:,i2,i3),&
    &x_out(:,i2,i3))
   END DO
   !$OMP END PARALLEL DO
  CASE(2)
   !$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i0,i1,i2,i3)
   DO i0=0,SIZE(x_out)/SIZE(x_out,dim)-1
    i1=lbound(x_out,1)+MOD(i0,SIZE(x_out,1))
    i3=lbound(x_out,3)+FLOOR(dble(i0)/dble(SIZE(x_out,3)))
    val_out(i1,:,i3)=&
    &interp1d(x_in(i1,:,i3),val_in(i1,:,i3),&
    &x_out(i1,:,i3))
   END DO
   !$OMP END PARALLEL DO
  CASE(3)
   !$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i0,i1,i2,i3)
   DO i0=0,SIZE(x_out)/SIZE(x_out,dim)-1
    i1=lbound(x_out,1)+MOD(i0,SIZE(x_out,1))
    i2=lbound(x_out,2)+FLOOR(dble(i0)/dble(SIZE(x_out,2)))
    val_out(i1,i2,:)=&
    &interp1d(x_in(i1,i2,:),val_in(i1,i2,:),&
    &x_out(i1,i2,:))
   END DO
   !$OMP END PARALLEL DO
  CASE DEFAULT
   WRITE(*,*)  "mod_interp.interp3d: dim must be 1,2 or 3"
   STOP
 END SELECT
  
 END FUNCTION interp3d

!----------------------------------------------------------------------------

 FUNCTION mesh2mesh_tl_interp4d(x_in,y_in,mask_in,val_in,&
 &x_out,y_out,mask_out) RESULT(val_out)
 ! val_out=mesh2mesh_tl_interp4d(x_in,y_in,mask_in,val_in,&
 ! &x_out,y_out)
 ! 
 ! Interpolates from fine to coarse grid
 ! x_in [1D double]: coordinates dimension 1 fine input grid
 ! y_in [1D double]: coordinates dimension 2 fine input grid
 ! mask_in [2D logical]: mask input grid
 ! val_in [4D double]: values on input grid
 ! x_out [1D double]: coordinates dimension 1 coarse output grid
 ! y_out [1D double]: coordinates dimension 2 coarse output grid
 ! mask_out [2D logical]: mask output grid

  IMPLICIT NONE
  REAL(8):: x_in(:),y_in(:),x_out(:),y_out(:)
  REAL(8):: val_in(:,:,:,:)
  REAL(8):: val_out(lbound(x_out,1):ubound(x_out,1),&
  &lbound(y_out,1):ubound(y_out,1),lbound(val_in,3):ubound(val_in,3),&
  &lbound(val_in,4):ubound(val_in,4))
  LOGICAL::mask_in(:,:),mask_out(:,:)

  REAL(8)::val_tmp(lbound(x_out,1):ubound(x_out,1),&
  &lbound(y_in,1):ubound(y_in,1),lbound(val_in,3):ubound(val_in,3),&
  &lbound(val_in,4):ubound(val_in,4))
  LOGICAL::mask_tmp(lbound(x_out,1):ubound(x_out,1),&
  &lbound(y_in,1):ubound(y_in,1))
  REAL(8)::w1(lbound(x_in,1):ubound(x_in,1))
  REAL(8)::w2(lbound(y_in,1):ubound(y_in,1))
  REAL(8),allocatable::wtmp(:)
  INTEGER::i1,i2,j1,j2,iL,iU

  !Mask_tmp
  DO i1=lbound(val_tmp,1),ubound(val_tmp,1)
  DO i2=lbound(val_tmp,2),ubound(val_tmp,2)
   IF(i1.EQ.lbound(val_tmp,1)) THEN
     mask_tmp(i1,i2)=ANY(&
     &x_in.GE.x_out(i1).AND.x_in.LE.x_out(i1+1).AND.mask_in(:,i2))
   ELSEIF(i1.EQ.ubound(val_tmp,1)) THEN
     mask_tmp(i1,i2)=ANY(&
     &x_in.GE.x_out(i1-1).AND.x_in.LE.x_out(i1).AND.mask_in(:,i2))
   ELSE
     mask_tmp(i1,i2)=ANY(&
     &x_in.GE.x_out(i1-1).AND.x_in.LE.x_out(i1+1).AND.mask_in(:,i2))
   END IF
  END DO
  END Do

  !Direction 1
  val_tmp=dble(0)
  DO i1=lbound(val_tmp,1),ubound(val_tmp,1)
   !Calculate weights
   IF(i1.EQ.lbound(val_tmp,1)) THEN
    w1=interp1d(x_out(i1:i1+1),dble([1,0]),x_in)
   ELSEIF(i1.EQ.ubound(val_tmp,1)) THEN
    w1=interp1d(x_out(i1-1:i1),dble([0,1]),x_in)
   ELSE
    w1=interp1d(x_out(i1-1:i1+1),dble([0,1,0]),x_in)
   END IF
   WHERE(w1.GE.interp_fill_value); w1=dble(0); END WHERE   

   IF(ALL(w1.EQ.dble(0))) CYCLE 
   iL=MINLOC(x_in,1,w1.NE.dble(0))
   iU=MAXLOC(x_in,1,w1.NE.dble(0))
 
   !Interp
   ALLOCATE(wtmp(iL:iU))
   DO i2=lbound(val_tmp,2),ubound(val_tmp,2)
    wtmp(iL:iU)=w1(iL:iU)
    WHERE(.NOT.mask_in(iL:iU,i2)); wtmp=dble(0); END WHERE
    IF(SUM(wtmp).LE.dble(0)) CYCLE
    wtmp=wtmp/SUM(wtmp)

    DO j1=iL,iU
     val_tmp(i1,i2,:,:)=val_tmp(i1,i2,:,:)+wtmp(j1)*val_in(j1,i2,:,:)
    END DO
   END DO !i2
   DEALLOCATE(wtmp)
  END DO !i1

  !Direction 2
  val_out=dble(0)
  DO i2=lbound(val_out,2),ubound(val_out,2)
   IF(i2.EQ.lbound(val_out,2)) THEN
    w2=interp1d(y_out(i2:i2+1),dble([1,0]),y_in)
   ELSEIF(i2.EQ.ubound(val_out,2)) THEN
    w2=interp1d(y_out(i2-1:i2),dble([0,1]),y_in)
   ELSE
    w2=interp1d(y_out(i2-1:i2+1),dble([0,1,0]),y_in)
   END IF
   WHERE(w2.EQ.interp_fill_value); w2=dble(0); END WHERE
   
   IF(ALL(w2.EQ.dble(0))) CYCLE
   iL=MINLOC(y_in,1,w2.NE.dble(0))
   iU=MAXLOC(y_in,1,w2.NE.dble(0))
   
   ALLOCATE(wtmp(iL:iU))
   DO i1=lbound(val_out,1),ubound(val_out,1)
    wtmp=w2(iL:iU)
    WHERE(.NOT.mask_tmp(i1,iL:iU)); wtmp=dble(0); END WHERE
    IF(SUM(wtmp).LE.dble(0)) CYCLE
    wtmp=wtmp/SUM(wtmp)
   
    DO j2=iL,iU
     val_out(i1,i2,:,:)=val_out(i1,i2,:,:)+wtmp(j2)*val_tmp(i1,j2,:,:)
    END DO
   END DO !i1
   DEALLOCATE(wtmp)
  END DO !i2
  
  !Mask out
  WHERE(SPREAD(SPREAD(.NOT.mask_out,3,size(val_out,3)),&
  &4,size(val_out,4)))
   val_out=dble(0)
  END WHERE

 END FUNCTION mesh2mesh_tl_interp4d

 FUNCTION mesh2mesh_tl_interp3d(x_in,y_in,mask_in,val_in,&
 &x_out,y_out,mask_out) RESULT(val_out)
  IMPLICIT NONE
  REAL(8):: x_in(:),y_in(:),x_out(:),y_out(:)
  REAL(8):: val_in(:,:,:)
  REAL(8):: val_out(lbound(x_out,1):ubound(x_out,1),&
  &lbound(y_out,1):ubound(y_out,1),lbound(val_in,3):ubound(val_in,3))
  LOGICAL::mask_in(:,:),mask_out(:,:)
  INTEGER::s_in(4),s_out(3)
  
  s_in=[size(val_in,1),size(val_in,2),size(val_in,3),1]
  s_out=[size(x_out),size(y_out),size(val_in,3)]
  val_out=RESHAPE(mesh2mesh_tl_interp4d(&
  &x_in,y_in,mask_in,RESHAPE(val_in,s_in),&
  &x_out,y_out,mask_out),s_out)
 END FUNCTION mesh2mesh_tl_interp3d

 FUNCTION mesh2mesh_tl_interp2d(x_in,y_in,mask_in,val_in,&
 &x_out,y_out,mask_out) RESULT(val_out)
  IMPLICIT NONE
  REAL(8):: x_in(:),y_in(:),x_out(:),y_out(:)
  REAL(8):: val_in(:,:)
  REAL(8):: val_out(lbound(x_out,1):ubound(x_out,1),&
  &lbound(y_out,1):ubound(y_out,1))
  LOGICAL::mask_in(:,:),mask_out(:,:)
  INTEGER::s_in(4),s_out(2)
  
  s_in=[size(val_in,1),size(val_in,2),1,1]
  s_out=[size(x_out),size(y_out)]
  val_out=RESHAPE(mesh2mesh_tl_interp4d(&
  &x_in,y_in,mask_in,RESHAPE(val_in,s_in),&
  &x_out,y_out,mask_out),s_out)
 END FUNCTION mesh2mesh_tl_interp2d

!----------------------------------------------------------------------------

 FUNCTION mesh2mesh_ad_interp4d(x_out,y_out,mask_out,val_out,&
 &x_in,y_in,mask_in) RESULT(val_in)
  !val_in=mesh2mesh_ad_interp4d(x_out,y_out,mask_in,val_out,&
  !&x_in,y_in)
  !
  !Adjoint of mesh2mesh_tl_interp4d

  IMPLICIT NONE
  REAL(8):: x_in(:),y_in(:),x_out(:),y_out(:)
  REAL(8):: val_out(:,:,:,:)
  REAL(8):: val_in(lbound(x_in,1):ubound(x_in,1),&
  &lbound(y_in,1):ubound(y_in,1),lbound(val_out,3):ubound(val_out,3),&
  &lbound(val_out,4):ubound(val_out,4))
  LOGICAL::mask_in(:,:),mask_out(:,:)

  REAL(8)::val_tmp(lbound(x_out,1):ubound(x_out,1),&
  &lbound(y_in,1):ubound(y_in,1),lbound(val_out,3):ubound(val_out,3),&
  &lbound(val_out,4):ubound(val_out,4))
  LOGICAL::mask_tmp(lbound(x_out,1):ubound(x_out,1),&
  &lbound(y_in,1):ubound(y_in,1))
  REAL(8)::w1(lbound(x_in,1):ubound(x_in,1))
  REAL(8)::w2(lbound(y_in,1):ubound(y_in,1))
  REAL(8),allocatable::wtmp(:)
  INTEGER::i1,i2,j1,j2,iL,iU

  !Mask_tmp
  DO i1=lbound(val_tmp,1),ubound(val_tmp,1)
  DO i2=lbound(val_tmp,2),ubound(val_tmp,2)
   IF(i1.EQ.lbound(val_tmp,1)) THEN
     mask_tmp(i1,i2)=ANY(&
     &x_in.GE.x_out(i1).AND.x_in.LE.x_out(i1+1).AND.mask_in(:,i2))
   ELSEIF(i1.EQ.ubound(val_tmp,1)) THEN
     mask_tmp(i1,i2)=ANY(&
     &x_in.GE.x_out(i1-1).AND.x_in.LE.x_out(i1).AND.mask_in(:,i2))
   ELSE
     mask_tmp(i1,i2)=ANY(&
     &x_in.GE.x_out(i1-1).AND.x_in.LE.x_out(i1+1).AND.mask_in(:,i2))
   END IF
  END DO
  END Do

  !Mask out
  WHERE(SPREAD(SPREAD(.NOT.mask_out,3,size(val_out,3)),&
  &4,size(val_out,4)))
   val_out=dble(0)
  END WHERE

  !Direction 2
  val_tmp=dble(0)
  DO i2=lbound(val_out,2),ubound(val_out,2)
   IF(i2.EQ.lbound(val_out,2)) THEN
    w2=interp1d(y_out(i2:i2+1),dble([1,0]),y_in)
   ELSEIF(i2.EQ.ubound(val_out,2)) THEN
    w2=interp1d(y_out(i2-1:i2),dble([0,1]),y_in)
   ELSE
    w2=interp1d(y_out(i2-1:i2+1),dble([0,1,0]),y_in)
   END IF
   WHERE(w2.EQ.interp_fill_value); w2=dble(0); END WHERE
   
   IF(ALL(w2.EQ.dble(0))) CYCLE
   iL=MINLOC(y_in,1,w2.NE.dble(0))
   iU=MAXLOC(y_in,1,w2.NE.dble(0))
   
   ALLOCATE(wtmp(iL:iU))
   DO i1=lbound(val_out,1),ubound(val_out,1)
    wtmp=w2(iL:iU)
    WHERE(.NOT.mask_tmp(i1,iL:iU)); wtmp=dble(0); END WHERE
    IF(SUM(wtmp).LE.dble(0)) CYCLE
    wtmp=wtmp/SUM(wtmp)
   
    DO j2=iL,iU
     val_tmp(i1,j2,:,:)=val_tmp(i1,j2,:,:)+wtmp(j2)*val_out(i1,i2,:,:)
    END DO
   END DO !i1
   DEALLOCATE(wtmp)
  END DO !i2

  !Direction 1
  val_in=dble(0)
  DO i1=lbound(val_tmp,1),ubound(val_tmp,1)
   !Calculate weights
   IF(i1.EQ.lbound(val_tmp,1)) THEN
    w1=interp1d(x_out(i1:i1+1),dble([1,0]),x_in)
   ELSEIF(i1.EQ.ubound(val_tmp,1)) THEN
    w1=interp1d(x_out(i1-1:i1),dble([0,1]),x_in)
   ELSE
    w1=interp1d(x_out(i1-1:i1+1),dble([0,1,0]),x_in)
   END IF
   WHERE(w1.GE.interp_fill_value); w1=dble(0); END WHERE   

   IF(ALL(w1.EQ.dble(0))) CYCLE 
   iL=MINLOC(x_in,1,w1.NE.dble(0))
   iU=MAXLOC(x_in,1,w1.NE.dble(0))
 
   !Interp
   ALLOCATE(wtmp(iL:iU))
   DO i2=lbound(val_tmp,2),ubound(val_tmp,2)
    wtmp(iL:iU)=w1(iL:iU)
    WHERE(.NOT.mask_in(iL:iU,i2)); wtmp=dble(0); END WHERE
    IF(SUM(wtmp).LE.dble(0)) CYCLE
    wtmp=wtmp/SUM(wtmp)

    DO j1=iL,iU
     val_in(j1,i2,:,:)=val_in(j1,i2,:,:)+wtmp(j1)*val_tmp(i1,i2,:,:)
    END DO
   END DO !i2
   DEALLOCATE(wtmp)
  END DO !i1
  
 END FUNCTION mesh2mesh_ad_interp4d

 FUNCTION mesh2mesh_ad_interp3d(x_out,y_out,mask_out,val_out,&
 &x_in,y_in,mask_in) RESULT(val_in)
  IMPLICIT NONE
  REAL(8):: x_in(:),y_in(:),x_out(:),y_out(:)
  REAL(8):: val_out(:,:,:)
  REAL(8):: val_in(lbound(x_in,1):ubound(x_in,1),&
  &lbound(y_in,1):ubound(y_in,1),lbound(val_out,3):ubound(val_out,3))
  LOGICAL::mask_in(:,:),mask_out(:,:)
  INTEGER::s_in(3),s_out(4)
  
  s_in=[size(x_in,1),size(y_in,1),size(val_out,3)]
  s_out=[size(x_out,1),size(y_out,1),size(val_out,3),1]
  val_in=RESHAPE(mesh2mesh_ad_interp4d(&
  &x_out,y_out,mask_out,RESHAPE(val_out,s_out),&
  &x_in,y_in,mask_in),s_in)
 END FUNCTION mesh2mesh_ad_interp3d

 FUNCTION mesh2mesh_ad_interp2d(x_out,y_out,mask_out,val_out,&
 &x_in,y_in,mask_in) RESULT(val_in)
  IMPLICIT NONE
  REAL(8):: x_in(:),y_in(:),x_out(:),y_out(:)
  REAL(8):: val_out(:,:)
  REAL(8):: val_in(lbound(x_in,1):ubound(x_in,1),&
  &lbound(y_in,1):ubound(y_in,1))
  LOGICAL::mask_in(:,:),mask_out(:,:)
  INTEGER::s_in(2),s_out(4)
  
  s_in=[size(val_in,1),size(val_in,2)]
  s_out=[size(x_out),size(y_out),1,1]
  val_in=RESHAPE(mesh2mesh_ad_interp4d(&
  &x_out,y_out,mask_out,RESHAPE(val_out,s_out),&
  &x_in,y_in,mask_in),s_in)
 END FUNCTION mesh2mesh_ad_interp2d
 
!------------------------------------------------------------------------------

 FUNCTION mesh2mesh_interp4d(x_in,y_in,mask_in,val_in,&
 &x_out,y_out) RESULT(val_out)
 !MESH2MESH_interp4d: Interpolation from one 2D orthogonal grid to another
 !
 !Syntax:
 ! val_out=mesh2mesh_interp2d(x_in,y_in,mask_in,val_in)
 !Input:
 ! x_in: 1D-double array with x-coordinates input grid
 ! y_in: 1D-double array with y-coordinates input grid
 ! mask-in: 2D logical array with mask input 
 ! x_out: 1D-double array with x-coordinates output grid
 ! y_out: 1D-double array with y-coordinates output grid
 ! val_in: 4D-double array with values on input grid
 !Output:
 ! val_out: 4D-double array with values on output grid

  IMPLICIT NONE
 
  REAL(8),INTENT(in),DIMENSION(:)::x_in,y_in,x_out,y_out
  REAL(8),INTENT(in),DIMENSION(:,:,:,:)::val_in
  REAL(8),INTENT(out)::val_out(size(x_out),size(y_out),&
  &size(val_in,3),size(val_in,4))
  LOGICAL,INTENT(in)::mask_in(:,:)

  LOGICAL::mask_tmp(lbound(val_in,1):ubound(val_in,1),&
  &lbound(y_out,1):ubound(y_out,1))
  LOGICAL::mask_out(lbound(x_out,1):ubound(x_out,1),&
  &lbound(y_out,1):ubound(y_out,1))
  REAL(8)::val_tmp(lbound(val_in,1):ubound(val_in,1),&
  &lbound(y_out,1):ubound(y_out,1),&
  &lbound(val_in,3):ubound(val_in,3),&
  &lbound(val_in,4):ubound(val_in,4))
  INTEGER::i1,i2,iL,iU
  REAL(8)::w(2)
  
  val_out=dble(0)
 
  !Interpolate to new y val
  val_tmp=interp_fill_value
  mask_tmp=.TRUE.
  DO i2=lbound(y_out,1),ubound(y_out,1)
   IF(y_out(i2).LT.MINVAL(y_in)) THEN
    val_tmp(:,i2,:,:)=interp_fill_value
    mask_tmp(:,i2)=.FALSE.
   ELSEIF(y_out(i2).GT.MAXVAL(y_in)) THEN
    val_tmp(:,i2,:,:)=interp_fill_value
    mask_tmp(:,i2)=.FALSE.
   ELSE
    iL=MAXLOC(y_in,1,y_in.LE.y_out(i2))
    iU=MINLOC(y_in,1,y_in.GE.y_out(i2))
    IF(iL.EQ.iU) THEN
     w=dble(1)
    ELSE
     w(1)=ABS(y_out(i2)-y_in(iU))
     w(2)=ABS(y_out(i2)-y_in(iL))
    END IF
    w=w/SUM(w)
   
    val_tmp(:,i2,:,:)=&
    &w(1)*val_in(:,iL,:,:)+w(2)*val_in(:,iU,:,:)
    DO i1=lbound(val_tmp,1),ubound(val_tmp,1)
     IF((.NOT.mask_in(i1,iL)).AND.(.NOT.mask_in(i1,iU))) THEN
      val_tmp(i1,i2,:,:)=interp_fill_value
      mask_tmp(i1,i2)=.FALSE.
     ELSEIF(.NOT.mask_in(i1,iL)) THEN
      val_tmp(i1,i2,:,:)=val_in(i1,iU,:,:)
     ELSEIF(.NOT.mask_in(i1,iU)) THEN
      val_tmp(i1,i2,:,:)=val_in(i1,iL,:,:)
     END IF
    END DO !i1
   END IF
  END DO !i2

  !Interpolate to new x val
  val_out=interp_fill_value
  mask_out=.TRUE.
  DO i1=lbound(x_out,1),ubound(x_out,1)
   IF(x_out(i1).LT.MINVAL(x_in)) THEN
    val_out(i1,:,:,:)=interp_fill_value
    mask_out(i1,:)=.FALSE.
   ELSEIF(x_out(i1).GT.MAXVAL(y_in)) THEN
     val_out(i1,:,:,:)=interp_fill_value
     mask_out(i1,:)=.FALSE.
   ELSE
    iL=MAXLOC(x_in,1,x_in.LE.x_out(i1))
    iU=MINLOC(x_in,1,x_in.GE.x_out(i1))
    IF(iL.EQ.iU) THEN
     w=dble(1)
    ELSE
     w(1)=ABS(x_out(i1)-x_in(iU))
     w(2)=ABS(x_out(i1)-x_in(iL))
    END IF
    w=w/SUM(w)
   
    val_out(i1,:,:,:)=&
    &w(1)*val_tmp(iL,:,:,:)+w(2)*val_tmp(iU,:,:,:)
    DO i2=lbound(val_out,2),ubound(val_out,2)
     IF((.NOT.mask_tmp(iL,i2)).AND.(.NOT.mask_tmp(iU,i2))) THEN
      val_out(i1,i2,:,:)=interp_fill_value
      mask_out(i1,i2)=.FALSE.
     ELSEIF(.NOT.mask_tmp(iL,i2)) THEN
      val_out(i1,i2,:,:)=val_tmp(iU,i2,:,:)
     ELSEIF(.NOT.mask_tmp(iU,i2)) THEN
      val_out(i1,i2,:,:)=val_tmp(iL,i2,:,:)
     END IF
    END DO !i2
    
   END IF
  END DO !i1

  WHERE(val_out.GE..99*interp_fill_value)
   val_out=interp_fill_value
  END WHERE

END FUNCTION mesh2mesh_interp4d

FUNCTION mesh2mesh_interp3d(x_in,y_in,mask_in,val_in,&
&x_out,y_out) RESULT(val_out)

  IMPLICIT NONE
 
  REAL(8),INTENT(in),DIMENSION(:)::x_in,y_in,x_out,y_out
  REAL(8),INTENT(in),DIMENSION(:,:,:)::val_in
  REAL(8),INTENT(out)::val_out(size(x_out),size(y_out),size(val_in,3))
  LOGICAL::mask_in(:,:)
  INTEGER::s_in(4),s_out(3) 

  s_in=1; s_in(1:3)=SHAPE(val_in)
  s_out=1; s_out(1:3)=SHAPE(val_out)

  val_out=RESHAPE(&
  &mesh2mesh_interp4d(x_in,y_in,mask_in,&
  &RESHAPE(val_in,s_in),x_out,y_out),s_out)

END FUNCTION mesh2mesh_interp3d

FUNCTION mesh2mesh_interp2d(x_in,y_in,mask_in,val_in,&
&x_out,y_out) RESULT(val_out)

  IMPLICIT NONE
 
  REAL(8),INTENT(in),DIMENSION(:)::x_in,y_in,x_out,y_out
  REAL(8),INTENT(in),DIMENSION(:,:)::val_in
  REAL(8),INTENT(out)::val_out(size(x_out),size(y_out))
  LOGICAL,INTENT(in)::mask_in(:,:)
  INTEGER::s_in(4),s_out(2) 

  s_in=1; s_in(1)=size(val_in,1); s_in(2)=size(val_in,2)
  s_out=1; s_out(1)=size(val_out,1); s_out(2)=size(val_out,2)
  val_out=dble(0)

  val_out=RESHAPE(&
  &mesh2mesh_interp4d(x_in,y_in,mask_in,&
  &RESHAPE(val_in,s_in),x_out,y_out),s_out)

END FUNCTION mesh2mesh_interp2d
!-----------------------------------------------------------------------------
SUBROUTINE mesh2mesh_extrap4d(x_in,y_in,mask_in,val,mask_out,fill_value)
!MESH2MESH_EXPTRAP4D(x_in,y_in,mask_in,val_in,mask_out) extrapolates values 
!inside mask_in to mask_out

 IMPLICIT NONE
 REAL(8),intent(in)::x_in(:),y_in(:)
 REAL(8),intent(inout)::val(:,:,:,:)
 REAL(8),optional::fill_value
 LOGICAL,intent(in)::mask_in(:,:),mask_out(:,:)

 INTEGER::i0,i1,i2,i3,i4,j0,j1,j2
 INTEGER::ntodo(2),step
 REAL(8)::w,w1
 LOGICAL::&
 &mask1(lbound(val,1):ubound(val,1),lbound(val,2):ubound(val,2)),&
 &mask2(lbound(val,1):ubound(val,1),lbound(val,2):ubound(val,2))

 WHERE(SPREAD(SPREAD(mask_out.AND..NOT.mask_in,3,size(val,3)),&
 &4,size(val,4)))
  val=dble(0)
 END WHERE

 mask1=mask_in
 mask2=mask_in
 
 step=1
 ntodo=COUNT(mask_out.AND..NOT.mask1)
 DO WHILE(ntodo(1).GT.0)

  DO i2=lbound(val,2),ubound(val,2)
  DO i1=lbound(val,1),ubound(val,1)
   IF(.NOT.mask_out(i1,i2).OR.mask1(i1,i2)) CYCLE

   w=0
   DO j2=MAX(i2-step,lbound(val,2)),MIN(i2+step,ubound(val,2))
   DO j1=MAX(i1-step,lbound(val,1)),MAX(i1+step,lbound(val,1))
    IF(.NOT.mask1(j1,j2)) CYCLE
    w1=(x_in(j1)-x_in(i1))**2+(y_in(j2)-y_in(i2))**2
    w=w+dble(1)/w1
    val(i1,i2,:,:)=val(i1,i2,:,:)+val(j1,j2,:,:)/w1    
   END DO
   END DO
   
   IF(w.NE.0) THEN
    val(i1,i2,:,:)=val(i1,i2,:,:)/w;
    mask2(i1,i2)=.TRUE.
   END IF 

  END DO
  END DO
  
  mask1=mask2; ntodo(2)=ntodo(1)
  ntodo(1)=COUNT(mask_out.AND..NOT.mask1)
  IF(ntodo(1).EQ.ntodo(2)) step=step+1
 END DO

 IF(PRESENT(fill_value)) THEN
  WHERE(SPREAD(SPREAD(.NOT.mask_out,3,size(val,3)),4,size(val,4)))
   val=fill_value
  END WHERE
 END IF

END SUBROUTINE mesh2mesh_extrap4d

SUBROUTINE mesh2mesh_extrap3d(x_in,y_in,mask_in,val,mask_out,fill_value)
  IMPLICIT NONE
  REAL(8),intent(in)::x_in(:),y_in(:)
  LOGICAL,intent(in)::mask_in(:,:),mask_out(:,:)
  REAL(8),intent(inout)::val(:,:,:)
  REAL(8),allocatable::val_tmp(:,:,:,:)
  REAL(8),optional::fill_value
  INTEGER::s_in(4),s_out(3)

  ALLOCATE(val_tmp(lbound(val,1):ubound(val,1),&
  &lbound(val,2):ubound(val,2),lbound(val,3):ubound(val,3),1))
  s_in=SHAPE(val_tmp); s_out=SHAPE(val)
  val_tmp=RESHAPE(val,s_in)
  IF(PRESENT(fill_value)) THEN
   CALL mesh2mesh_extrap4d(x_in,y_in,mask_in,val_tmp,mask_out,fill_value)
  ELSE
   CALL mesh2mesh_extrap4d(x_in,y_in,mask_in,val_tmp,mask_out)
  END IF
  val=RESHAPE(val_tmp,s_out)
END SUBROUTINE mesh2mesh_extrap3d

SUBROUTINE mesh2mesh_extrap2d(x_in,y_in,mask_in,val,mask_out,fill_value)
  IMPLICIT NONE
  REAL(8),intent(in)::x_in(:),y_in(:)
  LOGICAL,intent(in)::mask_in(:,:),mask_out(:,:)
  REAL(8),intent(inout)::val(:,:)
  REAL(8),allocatable::val_tmp(:,:,:,:)
  REAL(8),optional::fill_value
  INTEGER::s_in(4),s_out(2)

  ALLOCATE(val_tmp(lbound(val,1):ubound(val,1),&
  &lbound(val,2):ubound(val,2),1,1))
  s_in=SHAPE(val_tmp); s_out=SHAPE(val)
  val_tmp=RESHAPE(val,s_in)
  IF(PRESENT(fill_value)) THEN
   CALL mesh2mesh_extrap4d(x_in,y_in,mask_in,val_tmp,mask_out,fill_value)
  ELSE
   CALL mesh2mesh_extrap4d(x_in,y_in,mask_in,val_tmp,mask_out)
  END IF
  val=RESHAPE(val_tmp,s_out)
END SUBROUTINE mesh2mesh_extrap2d

!------------------------------------------------------------------------------

 FUNCTION mesh_interp2d(x_in,y_in,mask_in,val_in,x_out,y_out,mask_out)&
  & RESULT(val_out)
 !MESH_INTERP2D Interpolation on mesh with constant spacing in x and y direction !
 !Syntax
 ! val_out= mesh_interp2d(x_in,y_in,mask_in,val_in,x_out,y_out,mask_out)
 !Input:
 ! x_in: 2D DOUBLE array with x-coordinates input
 ! y_in: 2D DOUBLE array with y-coordinates input
 ! mask_in: 2D INTEGER array where masked out input points have value 0
 ! x_out: 2D DOUBLE array with x-coordinates output
 ! y_out: 2D DOUBLE array with y-coordinates output
 ! mask_out: 2D INTEGER array where masked out output points have value 0
 ! val_in: 2D DOUBLE  array containing the values on input points
 !Output
 ! val_out: 2D DOUBLE array containing the interpolated values.

  IMPLICIT NONE 

  REAL(8),INTENT(in),DIMENSION(:,:)::x_in,y_in,x_out,y_out,val_in
  INTEGER,INTENT(in):: mask_in(:,:), mask_out(:,:)
  REAL(8)::val_out(lbound(x_out,1):ubound(x_out,1),&
  &lbound(x_out,2):ubound(x_out,2))
  INTEGER::b(2,2),dimX,dimY
  REAL(8)::ind_d1,ind_d2,r1,r2,dx,dy,mask(4)
  INTEGER::ind_d1f,ind_d2f,i1,i2

  !Lower and upper bounds of the indices of grd_in
  b(1,1)=LBOUND(x_in,1)
  b(1,2)=UBOUND(x_in,1)
  b(2,1)=lBOUND(x_in,2)
  b(2,2)=UBOUND(x_in,2)

  !Create/test output and test input
  IF (ANY(SHAPE(val_out).NE.SHAPE(x_out)) ) THEN
   WRITE(*,*) "mod_interp.mesh_interp2d: dimensions val_out incongruent with&
    & dimensions output grid"
   STOP
  END IF 
  IF ( ANY(SHAPE(val_in).NE.SHAPE(x_in)) ) THEN
   WRITE(*,*) "mod_interp.mesh_interp2d: dimensions val_in incongruent with&
    & dimensions input grid"
    STOP
  END IF

  IF (MAXVAL(x_out).GT.MAXVAL(x_in).OR.&
   &MINVAL(x_out).LT.MINVAL(x_in).OR.&
   &MAXVAL(y_out).GT.MAXVAL(y_in).OR.&
   &MINVAL(y_out).LT.MINVAL(y_in)) THEN
   !WRITE(*,*) "mod_interp.mesh_interp2d: output grid is larger than input&
   ! & grid."
  END IF

  !Determine in which dimension x and y change
  IF ( ABS(x_in(b(1,1),b(2,2))-x_in(b(1,1),b(2,1))).GT.&
   & ABS(x_in(b(1,2),b(2,1))-x_in(b(1,1),b(2,1))) ) THEN
   dimX=2
   dimY=1
   dx=(MAXVAL(x_in(:,b(2,2)))-MINVAL(x_in(:,b(2,1))))&
    &/REAL(SIZE(x_in,2)-1)
   dy=(MAXVAL(y_in(b(1,2),:))-MINVAL(y_in(b(1,1),:)))&
    &/REAL(SIZE(y_in,1)-1)
  ELSE
   dimX=1
   dimY=2
   dy=(MAXVAL(y_in(:,b(2,2)))-MINVAL(y_in(:,b(2,1))))&
    &/REAL(SIZE(y_in,2)-1)
   dx=(MAXVAL(x_in(b(1,2),:))-MINVAL(x_in(b(1,1),:)))&
    &/REAL(SIZE(x_in,1)-1)
  END IF

  !Calculate output
  val_out=interp_fill_value
  DO i1=lbound(val_out,1),ubound(val_out,1)
  DO i2=lbound(val_out,2),ubound(val_out,2)

   IF(mask_out(i1,i2).EQ.0) CYCLE 

   !Grid coordinates of the 4 points surrounding the interpolation point
   IF (dimX.EQ.1) THEN
    ind_d1=(x_out(i1,i2)-x_in(b(1,1),b(2,1)))/dx+lbound(val_out,1)
    ind_d2=(y_out(i1,i2)-y_in(b(1,1),b(2,1)))/dy+lbound(val_out,2)
   ELSE
    ind_d2=(x_out(i1,i2)-x_in(b(1,1),b(2,1)))/dx+lbound(val_out,2)
    ind_d1=(y_out(i1,i2)-y_in(b(1,1),b(2,1)))/dy+lbound(val_out,1)
   END IF

   ind_d1f=FLOOR(ind_d1)
   ind_d2f=FLOOR(ind_d2)
   ind_d1f=MINVAL([ind_d1f,SIZE(x_in,1)-1])
   ind_d1f=MAXVAL([ind_d1f,1])
   ind_d2f=MINVAL([ind_d2f,SIZE(x_in,2)-1])
   ind_d2f=MAXVAL([ind_d2f,1])

   !Check if interpolation point lies within boundaries grid
   r1=ind_d1-REAL(ind_d1f)
   r2=ind_d2-REAL(ind_d2f)
   IF (r1.GT.1.0.OR.r2.GT.1.OR.r1.LT.0.0.OR.r2.LT.0.0) THEN
    CYCLE
    !STOP 'mesh_interp2d:interpolation points are located outside grid'
   END IF
    
   !Perform bilinear interpolation
   IF(ANY(mask_in(ind_d1f:ind_d1f+1,ind_d2f:ind_d2f+1).EQ.0)) CYCLE
   val_out(i1,i2)=val_in(ind_d1f,ind_d2f)*(1-r1)*(1-r2)+&
    &val_in(ind_d1f+1,ind_d2f)*r1*(1-r2)+&
    &val_in(ind_d1f+1,ind_d2f+1)*r1*r2+&
    &val_in(ind_d1f,ind_d2f+1)*(1-r1)*r2  
  END DO
  END DO

 END FUNCTION mesh_interp2d


!------------------------------------------------------------------------------

 FUNCTION mesh_interp3d(x_in,y_in,mask_in,z_in,val_in,&
  &x_out,y_out,mask_out,z_out) RESULT(val_out)
 !MESH_INTERP3D interpolates a 3D field from 1 grid to another.
 !
 !Syntax:
 ! val_out=mesh_interp3d(x_in,y_in,mask_in,z_in,val_in,x_out,y_out,
 !         mask_out,z_out)
 !Input:
 ! x_in: 2D DOUBLE array with x-coordinates input
 ! y_in: 2D DOUBLE array with y-coordinate input
 ! mask_in: 2D INTEGER array with 0s at positions where input is masked out
 ! z_in: 3D DOUBLE array with z-coordinates input
 ! val_in: 3D DOUBLE with values at input coordinates
 ! x_out: 2D DOUBLE array with x-coordinates output
 ! y_out: 2D DOUBLE array with y-coordinates output
 ! mask_out: 2D INTEGER array with 0s at positions where outpu is masked out
 ! z_out: 3D DOUBLE aray with z-coordinates output
 !Output:
 ! val_out: 3D DOUBLE array with interpolated values. 

  IMPLICIT NONE

  REAL(8),DIMENSION(:,:),INTENT(in)::x_in,y_in,x_out,y_out
  REAL(8),DIMENSION(:,:,:),INTENT(in)::z_in,z_out,val_in
  INTEGER,DIMENSION(:,:),INTENT(in)::mask_in,mask_out
  REAL(8)::val_out(lbound(x_out,1):ubound(x_out,1),&
  &lbound(x_out,2):ubound(x_out,2),lbound(val_in,3):ubound(val_in,3))
  REAL(8)::val_int(lbound(z_out,1):ubound(z_out,1),&
   &lbound(z_out,1):ubound(z_out,2),lbound(z_in,3):ubound(z_in,3))
  REAL(8)::z_int(lbound(z_out,1):ubound(z_out,1),&
   &lbound(z_out,1):ubound(z_out,2),lbound(z_in,3):ubound(z_in,3))
  REAL(8)::val1(SIZE(z_out,3))
  INTEGER::i1,i2,i3

  !Check dimensions
  IF (size(x_in,1).NE.SIZE(z_in,1).OR.&
   &size(x_in,2).NE.SIZE(z_in,2)) THEN
   WRITE(*,*) "mod_interp.mesh_interp3d: dimensions input grid and &
    & z_in incongruent"
   STOP
  END IF
  IF (size(x_out,1).NE.SIZE(z_out,1).OR.&
   &size(x_out,2).NE.SIZE(z_out,2)) THEN
   WRITE(*,*) "mod_interp.mesh_interp3d: dimensions grd_out and z_out &
    &incongruent" 
   STOP
  END IF
  IF ( ANY(SHAPE(z_in).NE.SHAPE(val_in))  ) THEN
   WRITE(*,*) "mod_interp.mesh_interp3d: dimensions z_in and val_in &
   &incongruent"
   STOP
  END IF
  IF ( ANY(SHAPE(z_out).NE.SHAPE(val_out)) ) THEN
   WRITE(*,*) "mod_interp.mesh_interp3d: dimensions z_out and &
    &val_out incongruent"
   STOP
  END IF

  !Interpolate in the horizontal
  DO i3=LBOUND(z_in,3),UBOUND(z_in,3)
    val_int(:,:,i3)=mesh_interp2d(x_in,y_in,mask_in,val_in(:,:,i3),&
     &x_out,y_out,mask_out)
    z_int(:,:,i3)=mesh_interp2d(x_in,y_in,mask_in,z_in(:,:,i3),&
     &x_out,y_out,mask_out)
  END DO 
  
  !Interpolate in the vertical
  DO i1=LBOUND(z_out,1),UBOUND(z_out,1)
  DO i2=LBOUND(z_out,2),UBOUND(z_out,2)
    val1=interp1d(z_int(i1,i2,:),val_int(i1,i2,:),&
     &z_out(i1,i2,:))
    val_out(i1,i2,:)=val1
  END DO
  END DO

 END FUNCTION mesh_interp3d
!------------------------------------------------------------------------------

 FUNCTION mesh_extrap2d(x,y,mask,val_in,opt) RESULT(val_out)
 !MESH_EXTRAP2D replace the value of unmasked points which have NaN value
 !with the average of the value of the closest non-nan points
 !
 !Syntax:
 ! val_out=mesh_extrap2d(x_in,y_in,mask_in,val_in,opt)
 !Input:
 ! x_in: 2D DOUBLE array with x-coordinates
 ! y_in 2D DOUBLE array with y-coordinates
 ! mask_in: 2D INTEGER array where masked out points have value 0
 ! val_in: 2D DOUBLE array with input values
 ! opt: OPT CHARACTER use 'maskgap' (default) to extrapolate only to unmasked 
 !      points. Use 'maskfill' to extrapolate to all points. 
 !Output:
 ! val_out: extrapolated values are added to the existing values in val_in. 
  
 IMPLICIT NONE
 
 REAL(8),DIMENSION(:,:),INTENT(in)::x,y,val_in
 INTEGER,INTENT(in)::mask(:,:)
 REAL(8)::val_out(lbound(val_in,1):ubound(val_in,1),&
 &lbound(val_in,2):ubound(val_in,2))
 REAL(8)::dist_field(size(x,1),size(x,2))
 CHARACTER(len=*),INTENT(in),OPTIONAL::opt
 CHARACTER(len=32)::opt_def
 REAL(8)::val_point,w_point,dist
 INTEGER::b(2,2),i1,i2,j1,j2,min_index(2),min_val
 CHARACTER(len=1024)::err_str
 
 !Set default for opt
 opt_def='maskgap'
 IF(PRESENT(opt)) THEN
  opt_def=TRIM(opt)
 END IF
  
 !Minimum and maximum indices grd
 b(1,1)=LBOUND(val_in,1)
 b(1,2)=UBOUND(val_in,1)
 b(2,1)=LBOUND(val_in,2)
 b(2,2)=UBOUND(val_in,2)

 !Check sizes
 IF(ANY(SHAPE(x).NE.SHAPE(val_in))) THEN
  WRITE(*,*) 'mod_interp.mesh_extrap2d: dimension x and val&
  & incongruent.'
  STOP
 ELSEIF(ANY(SHAPE(y).NE.SHAPE(val_in))) THEN
  WRITE(*,*) 'mod_interp.mesh_extrap2d: dimension y and val&
  & incongruent.'
  STOP
 ELSEIF(ANY(SHAPE(mask).NE.SHAPE(val_in))) THEN
  WRITE(*,*) 'mod_interp.mesh_extrap2d: dimension mask and val_in&
   & incongruent.'
  STOP
 END IF

 !Copy input to output
 val_out=val_in

 !Perform extrapolation
 DO i2=b(2,1),b(2,2)
 DO i1=b(1,1),b(1,2) 
  
  !If point is masked out continue to nex point
  IF(strcmp(opt_def,'maskgap')) THEN
   IF(mask(i1,i2).EQ.0) CYCLE
  END IF

  !If point has non-nan value continue to next point
  IF(val_in(i1,i2).NE.interp_fill_value) CYCLE

  w_point=DBLE(0); val_point=DBLE(0); 
  !Search for non-nan point to the left
  DO j1=i1,b(1,1),-1
   IF(val_in(j1,i2).NE.interp_fill_value) THEN
    dist=(x(j1,i2)-x(i1,i2))**2+(y(j1,i2)-y(i1,i2))**2
    val_point=val_point+val_in(j1,i2)/dist
    w_point=w_point+1.0/dist
    EXIT    
   END IF
  END DO
  !Search for non-nan point to the right
  DO j1=i1,b(1,2)
   IF(val_in(j1,i2).NE.interp_fill_value) THEN
    dist=(x(j1,i2)-x(i1,i2))**2+(y(j1,i2)-y(i1,i2))**2
    val_point=val_point+val_in(j1,i2)/dist
    w_point=w_point+1.0/dist
    EXIT    
   END IF
  END DO
  !Search for non-nan point above
  DO j2=i2,b(2,1),-1
   IF(val_in(i1,j2).NE.interp_fill_value) THEN
    dist=(x(i1,j2)-x(i1,i2))**2+(y(i1,j2)-y(i1,i2))**2
    val_point=val_point+val_in(i1,j2)/dist
    w_point=w_point+1.0/dist
    EXIT    
   END IF
  END DO
 !Search for non-nan point below
  DO j2=i2,b(2,2)
   IF(val_in(i1,j2).NE.interp_fill_value) THEN
    dist=(x(i1,j2)-x(i1,i2))**2+(y(i1,j2)-y(i1,i2))**2
    val_point=val_point+val_in(i1,j2)/dist
    w_point=w_point+1.0/dist
    EXIT    
   END IF
  END DO

  !IF the previous fails used find the closest point
  
  IF(w_point.EQ.0) THEN
   dist_field= (x-x(i1,i2))**2+(y-y(i1,i2))**2
   min_val=HUGE(min_val) 
   DO j2=b(2,1),b(2,2)
   DO j1=b(1,1),b(1,2)
    IF(dist_field(j1,j2).LT.min_val.AND.&
     &val_in(j1,j2).NE.interp_fill_value) THEN
     min_val=dist_field(j1,j2)
     min_index=[j1,j2]
    END IF 
   END DO
   END DO
   
   IF(min_val.EQ.1e36) THEN
    WRITE(err_str,*) 'mod_interp.mesh_extrap2d: cannot&
    & extrapolate point ',i1,' ',i2
    WRITE(*,*) TRIM(err_str)
    STOP
   ELSE
    val_out(i1,i2)=val_in(min_index(1),min_index(2))
   END IF
  ELSE
   val_out(i1,i2)=val_point/w_point
  END IF

 END DO
 END DO

 END FUNCTION mesh_extrap2d

!-----------------------------------------------------------------------------

 FUNCTION mesh_extrap3d(x,y,mask,val_in,opt) RESULT(val_out)
 !MESH_EXTRAP3D replace the value of unmasked points which have NaN value
 !with the average of the value of the closest non-nan points
 !
 !Syntax:
 ! val_out=mesh_extrap2d(x_in,y_in,mask_in,val_in,opt)
 !Input:
 ! x_in: 3D DOUBLE array with x-coordinates
 ! y_in 3D DOUBLE array with y-coordinates
 ! mask_in: 3D INTEGER array where masked out points have value 0
 ! val_in: 3D DOUBLE array with input values
 ! opt: OPT CHARACTER use 'maskgap' (default) to extrapolate only to unmasked 
 !      points. Use 'maskfill' to extrapolate to all points
 !Output:
 ! val_out: extrapolated values are added to the existing values in val_in
  
 IMPLICIT NONE

 REAL(8),INTENT(in)::x(:,:),y(:,:),val_in(:,:,:)
 CHARACTER(len=*),INTENT(in),OPTIONAL::opt
 REAL(8)::val_out(lbound(val_in,1):ubound(val_in,1),&
 &lbound(val_in,2):ubound(val_in,2),lbound(val_in,3):ubound(val_in,3))
 CHARACTER(len=32)::opt_def
 INTEGER,INTENT(in)::mask(:,:)
 REAL(8)::val_tmp(SIZE(val_in,1),SIZE(val_in,2))
 INTEGER::i3

  !Set default for opt
  opt_def='maskgap'
  IF(PRESENT(opt)) THEN
   opt_def=TRIM(opt)
  END IF

 !Check dimension
 IF(ANY(SHAPE(val_in).NE.SHAPE(val_out))) THEN
  WRITE(*,*) 'mod_interp.mesh_extrap3d: dimensions val_in and&
   &val_out incongruent.'
  STOP
 ELSEIF(SIZE(val_in,1).NE.SIZE(x,1).OR.SIZE(val_in,2).NE.SIZE(x,2)) THEN
  WRITE(*,*) 'mod_interp.mesh_extrap3d: dimensions val_in and&
   &x incongruent.'
  STOP
 ELSEIF(SIZE(val_in,1).NE.SIZE(y,1).OR.SIZE(val_in,2).NE.SIZE(y,2)) THEN
  WRITE(*,*) 'mod_interp.mesh_extrap3d: dimensions val_in and&
   &y incongruent.'
  STOP
 ELSEIF(SIZE(val_in,1).NE.SIZE(mask,1).OR.&
  &SIZE(val_in,2).NE.SIZE(mask,2)) THEN
  WRITE(*,*) 'mod_interp.mesh_extrap3d: dimensions val_in and&
   &mask incongruent.'
  STOP
 END IF

 !Extrapolate layer by layer
 DO i3=LBOUND(val_in,3),UBOUND(val_in,3)
  val_tmp=mesh_extrap2d(x,y,mask,val_in(:,:,i3),opt_def)
  val_out(:,:,i3)=val_tmp
 END DO


 END FUNCTION mesh_extrap3d

!-------------------------------------------------------------------------------
 FUNCTION  z2sigma_2d(param,h,zeta,z) RESULT(sigma)
 !Z2SIGMA Convert z-coordinates to sigma coordinates
 !
 !Syntax:
 ! sigma=z2sigma_2d(param,h,zeta,z)
 !Input:
 ! param: TYPE(sigma_param) containing the parameters for the transformation
 ! z: 3D double array with z-coordinates
 ! h: 2D double array with bottom depth w.r.t. reference level
 ! zeta: 2D double array with water level w.r.t. reference level
 !Output
 ! sigma: 3D double array with sigma-coordinates corresponding to z-coordinates

 IMPLICIT NONE

  TYPE(sigma_param),INTENT(in)::param
  REAL(8),INTENT(in)::z(:,:,:),h(:,:),zeta(:,:)
  REAL(8)::sigma(lbound(z,1):ubound(z,1),&
   &lbound(z,2):ubound(z,2),lbound(z,3):ubound(z,3))
  REAL(8),ALLOCATABLE::sigma_int(:),z_int(:,:,:),sigma1(:)
  INTEGER::dim(3),i1,i2,i3
  INTEGER,PARAMETER::n_int=100 !Number of levels in vertical used for interpolation

  !Dimensions input
  dim(1)=SIZE(z,1)
  dim(2)=SIZE(z,2)
  dim(3)=SIZE(z,3)
  IF (SIZE(h,1).NE.dim(1).OR.SIZE(h,2).NE.dim(2)) THEN
   WRITE(*,*) 'mod_interp.z2sigma: dimensions h and z incongruent'
   STOP
  END IF
  IF (SIZE(zeta,1).NE.dim(1).OR.SIZE(zeta,2).NE.dim(2)) THEN
   WRITE(*,*) 'mod_interp.z2sigma: dimensions zeta and z incongruent'
   STOP
  END IF
  IF (SIZE(sigma,1).NE.dim(1).OR.SIZE(sigma,2).NE.dim(2)&
   &.OR.SIZE(sigma,3).NE.dim(3)) THEN
   WRITE(*,*) 'mod_interp.z2sigma: dimensions sigma and z incongruent' 
   STOP
  END IF

  !Sigma levels for interpolation
  ALLOCATE(sigma_int(n_int))
  DO i3=1,n_int
   sigma_int(i3)=-dble(n_int-i3)/dble(n_int-1)
  END DO

  !z-coordinates at aforementioned sigma-coordinates
  ALLOCATE(z_int(lbound(z,1):ubound(z,1),lbound(z,2):ubound(z,2),n_int))
  z_int=sigma2z(param,h,zeta,sigma_int)

  !For each horizontal point interpolate in the vertical to find sigma values
  sigma=interp_fill_value
  ALLOCATE(sigma1(LBOUND(sigma,3):UBOUND(sigma,3)))
  DO i1=lbound(z,1),ubound(z,1)
  DO i2=lbound(z,2),ubound(z,2)
   IF(zeta(i1,i2).EQ.interp_fill_value) CYCLE
   IF(h(i1,i2).EQ.interp_fill_value) CYCLE
   sigma1=interp1d(z_int(i1,i2,:),sigma_int,z(i1,i2,:))
   sigma(i1,i2,:)=sigma1
  END DO
  END DO

 END FUNCTION z2sigma_2d

!-------------------------------------------------------------------------------
 FUNCTION  z2sigma_scalar(param,h,zeta,z) RESULT(sigma)
 !Z2SIGMA Convert z-coordinates to sigma coordinates
 !
 !Syntax:
 ! sigma=z2sigma_scalar(param,h,zeta,z)
 !Input:
 ! param: TYPE(sigma_param) containing the parameters for the transformation
 ! z: 1D double array with z-coordinates
 ! h: double  with bottom depth w.r.t. reference level
 ! zeta: double with water level w.r.t. reference level
 !Output
 ! sigma: 1D double array with sigma-coordinates corresponding to z-coordinates

 IMPLICIT NONE

  TYPE(sigma_param),INTENT(in)::param
  REAL(8),INTENT(in)::z(:)
  REAL(8),INTENT(in)::h,zeta
  REAL(8)::sigma(lbound(z,1):ubound(z,1))
  REAL(8),ALLOCATABLE::sigma_int(:),z_int(:),sigma1(:)
  INTEGER::dim,i1,i2,i3
  INTEGER,PARAMETER::n_int=100 !Number of levels in vertical used for interpolation

  !Dimensions input
  IF (SIZE(sigma,1).NE.dim) THEN
   WRITE(*,*) 'mod_interp.z2sigma: dimensions sigma and z incongruent' 
   STOP
  END IF

  !Sigma levels for interpolation
  ALLOCATE(sigma_int(n_int))
  DO i3=1,n_int
   sigma_int(i3)=-dble(n_int-i3)/dble(n_int-1)
  END DO

  !z-coordinates at aforementioned sigma-coordinates
  ALLOCATE(z_int(n_int))
  z_int=sigma2z(param,h,zeta,sigma_int)

  !For each horizontal point interpolate in the vertical to find sigma values
  sigma=interp_fill_value
  ALLOCATE( sigma1(LBOUND(sigma,1):UBOUND(sigma,1)) )
  IF(h.NE.interp_fill_value.OR.zeta.NE.interp_fill_value) THEN
   sigma1=interp1d(z_int,sigma_int,z)
   sigma=sigma1
  END IF

 END FUNCTION z2sigma_scalar

!-------------------------------------------------------------------------------  
 FUNCTION sigma2z_2d(param,h,zeta,sigma) RESULT(z) 
 !SIGMA2Z Convert a sigma-coordinates to a grid of z-coordinates
 !
 !Syntax:
 ! z=sigma2z_2d(param,h,zeta,sigma)
 !Input:
 ! param: TYPE(sigma_param) containing the parameters for the transformation
 ! sigma: 1D double array with sigma coordinates for grid
 ! h: 2D double array with distance reference level bottom
 ! zeta: 2D double array with water level w.r.t. reference level
 !Output:
 ! z: 3D double array with z-coordinates

 IMPLICIT NONE

  REAL(8),INTENT(in)::sigma(:)
  TYPE(sigma_param),INTENT(in)::param
  REAL(8),DIMENSION(:,:),INTENT(in)::h,zeta
  REAL(8),ALLOCATABLE::S(:,:),C(:),mu(:)
  REAL(8)::z(lbound(h,1):ubound(h,1),&
   &lbound(h,2):ubound(h,2),lbound(sigma,1):ubound(sigma,1))
  REAL(8)::hc
  INTEGER::diml(3),i1,i2,i3

  !Dimension grid
  diml=(/SIZE(z,1),SIZE(z,2),SIZE(z,3)/)

  !Check dimensions
  IF (diml(1).NE.SIZE(h,1).OR.diml(2).NE.SIZE(h,2)) THEN
   WRITE(*,*) 'mod_interp.sigma2z: dimensions h and z incongruent'
   STOP 
 END IF
  IF (diml(1).NE.SIZE(zeta,1).OR.diml(2).NE.SIZE(zeta,2)) THEN
   WRITE(*,*) 'mod_interp.sigma2z: dimensions zeta and z incongruent'
   STOP
  END IF
  IF (diml(3).NE.SIZE(sigma,1)) THEN
   WRITE(*,*) 'mod_interp.sigma2z: dimensions sigma and z incongruent'
   STOP
  END iF
 
  !Calculate vertical stretching function
  ALLOCATE(C(lbound(sigma,1):ubound(sigma,1)))
  SELECT CASE(param%v_stretching)
  CASE (1) !Vstretching=1
   C=(1-param%theta_b)*SINH(param%theta_s*sigma)/SINH(param%theta_s)+&
    &0.5*param%theta_b*( TANH(param%theta_s*(sigma+0.5))/&
    &TANH(0.5*param%theta_s)-1.0)
  CASE (2) !Vstretching=2
   ALLOCATE(mu(lbound(sigma,1):ubound(sigma,1)))
   mu=(sigma+1)**param%alpha*(1+param%alpha/param%beta*&
    &(1-(sigma+1)**param%beta))
   C=mu*(1-COSH(param%theta_s*sigma))/(COSH(param%theta_s)-1)+&
    &(1-mu)*(SINH(param%theta_b*(sigma+1))/SINH(param%theta_b)-1)
  CASE (3) !Vstretching=3
   ALLOCATE(mu(lbound(sigma,1):ubound(sigma,1)))
   mu=0.5*(1-TANH(param%gamma*(sigma+0.5)))
   C=-(1-mu)*LOG(COSH(param%gamma*ABS(sigma)**param%theta_s))&
    &/LOG(COSH(param%gamma))+mu*&
    &(LOG(COSH(param%gamma*(sigma+1)**param%theta_b))/&
    &LOG(COSH(param%gamma))-1)
   CASE (4) !Vstretching=4
    IF (param%theta_s.GT.0.0) THEN
     C=(1-COSH(param%theta_s*sigma))/(COSH(param%theta_s)-1)
    ELSE
     C=-sigma**2
    END IF
    C=(EXP(param%theta_b*C)-1)/(1-EXP(-param%theta_b))
   CASE DEFAULT
    WRITE(*,*) 'mod_interp.sigma2z: not a valid option for Vstretching.'
    STOP
   END SELECT

  !Calculate z coordinates
  ALLOCATE(S(lbound(z,1):ubound(z,1),lbound(z,2):ubound(z,2)))
  SELECT CASE(param%v_transform)
  CASE (1) !Vtransform=1
   hc=MINVAL((/MINVAL(h),param%t_cline/))
   DO i3=lbound(sigma,1),ubound(sigma,1) 
    S(:,:)=hc*sigma(i3)+(h-hc)*C(i3)
    z(:,:,i3)=S+zeta*(1+S/h)
   END DO
  CASE (2) !Vtransform=2
   hc=param%t_cline
   DO i3=lbound(sigma,1),ubound(sigma,1)
    S=(hc*sigma(i3)+h*C(i3))/(hc+h)
    z(:,:,i3)=zeta+(zeta+h)*S
   END DO
  CASE DEFAULT
   WRITE(*,*) 'mod_interp.sigma2z: not a valid option for Vtransform.'
   STOP
  END SELECT  
 
  DO i1=lbound(z,1),ubound(z,1)
  DO i2=lbound(z,2),ubound(z,2)
    IF(zeta(i1,i2).EQ.interp_fill_value.OR.&
     &h(i1,i2).EQ.interp_fill_value) THEN
     z(i1,i2,:)=interp_fill_value
    END IF
  END DO
  END DO
 
 END FUNCTION sigma2z_2d

!------------------------------------------------------------------------

FUNCTION sigma2z_scalar(param,h,zeta,sigma) RESULT(z) 
 !SIGMA2Z Convert a sigma-coordinates to a grid of z-coordinates
 !
 !Syntax:
 ! z=sigma2z(param,h,zeta,sigma)
 !Input:
 ! param: TYPE(sigma_param) containing the parameters for the transformation
 ! sigma: 1D double array with sigma coordinates for grid
 ! h: double with distance reference level bottom
 ! zeta: double with water level w.r.t. reference level
 !Output:
 ! z: 1D double array with z-coordinates
 !ATTENTION: Tcline must be smaller or equal to min(h)

 IMPLICIT NONE

  REAL(8),INTENT(in)::sigma(:)
  TYPE(sigma_param),INTENT(in)::param
  REAL(8),INTENT(in)::h,zeta
  REAL(8),ALLOCATABLE::C(:),mu(:)
  REAL(8)::z(lbound(sigma,1):ubound(sigma,1))
  REAL(8)::hc,S
  INTEGER::diml,i1,i2,i3

  !Dimension grid
  diml=SIZE(z,1)

  !Check dimensions
  IF (diml.NE.SIZE(sigma,1)) THEN
   WRITE(*,*) 'mod_interp.sigma2z: dimensions sigma and z incongruent'
   STOP
  END IF
 
  !Calculate vertical stretching function
  ALLOCATE(C(lbound(sigma,1):ubound(sigma,1)))
  SELECT CASE(param%v_stretching)
  CASE (1) !Vstretching=1
   C=(1-param%theta_b)*SINH(param%theta_s*sigma)/SINH(param%theta_s)+&
    &0.5*param%theta_b*( TANH(param%theta_s*(sigma+0.5))/&
    &TANH(0.5*param%theta_s)-1.0)
  CASE (2) !Vstretching=2
   ALLOCATE(mu(lbound(sigma,1):ubound(sigma,1)))
   mu=(sigma+1)**param%alpha*(1+param%alpha/param%beta*&
    &(1-(sigma+1)**param%beta))
   C=mu*(1-COSH(param%theta_s*sigma))/(COSH(param%theta_s)-1)+&
    &(1-mu)*(SINH(param%theta_b*(sigma+1))/SINH(param%theta_b)-1)
  CASE (3) !Vstretching=3
   ALLOCATE(mu(lbound(sigma,1):ubound(sigma,1)))
   mu=0.5*(1-TANH(param%gamma*(sigma+0.5)))
   C=-(1-mu)*LOG(COSH(param%gamma*ABS(sigma)**param%theta_s))&
    &/LOG(COSH(param%gamma))+mu*&
    &(LOG(COSH(param%gamma*(sigma+1)**param%theta_b))/&
    &LOG(COSH(param%gamma))-1)
   CASE (4) !Vstretching=4
    IF (param%theta_s.GT.0.0) THEN
     C=(1-COSH(param%theta_s*sigma))/(COSH(param%theta_s)-1)
    ELSE
     C=-sigma**2
    END IF
    C=(EXP(param%theta_b*C)-1)/(1-EXP(-param%theta_b))
   CASE DEFAULT
    WRITE(*,*) 'mod_interp.sigma2z: not a valid option for Vstretching.'
    STOP
   END SELECT

  !Calculate z coordinates
  SELECT CASE(param%v_transform)
  CASE (1) !Vtransform=1
   hc=param%t_cline
   DO i3=lbound(sigma,1),ubound(sigma,1) 
    S=hc*sigma(i3)+(h-hc)*C(i3)
    z(i3)=S+zeta*(1+S/h)
   END DO
  CASE (2) !Vtransform=2
   hc=param%t_cline
   DO i3=lbound(sigma,1),ubound(sigma,1)
    S=(hc*sigma(i3)+h*C(i3))/(hc+h)
    z(i3)=zeta+(zeta+h)*S
   END DO
  CASE DEFAULT
   WRITE(*,*) 'mod_interp.sigma2z: not a valid option for Vtransform.'
   STOP
  END SELECT  
 

  IF(zeta.EQ.interp_fill_value.OR.&
   &h.EQ.interp_fill_value) THEN
   z(:)=interp_fill_value
  END IF
 
 END FUNCTION sigma2z_scalar

!------------------------------------------------------------------------

 FUNCTION zweight(param,h,zeta,dimvar) RESULT(zw)
 !ZWEIGHT calculates w such that SUM(w*val,3) gives the depth-averaged value
 !of val. Use dimvar=1 for u-points, dimvar=2 for v-points and dimvar=0 for
 !rho-points

 TYPE(sigma_param),intent(in)::param
 REAL(8),intent(in)::h(:,:),zeta(:,:)
 INTEGER,intent(in)::dimvar
 REAL(8),allocatable::zw(:,:,:)
 REAL(8),allocatable::hvar(:,:),zetavar(:,:),z(:,:,:),sigma(:)
 INTEGER::nsize(2),i0
 
  !Size grid
  nsize=shape(h)


  IF(dimvar.EQ.1) THEN
    !weights for u-points
    ALLOCATE(hvar(nsize(1)-1,nsize(2)))
    ALLOCATE(zetavar(nsize(1)-1,nsize(2)))
    ALLOCATE(zw(nsize(1)-1,nsize(2),param%n_sigma))
    hvar=.5*h(1:nsize(1)-1,:)+.5*h(2:nsize(1),:)
    zetavar=.5*zeta(1:nsize(1)-1,:)+.5*zeta(2:nsize(1),:)
    WHERE( zeta(1:nsize(1)-1,:).EQ.interp_fill_value .OR. &
    & zeta(2:nsize(1),:).EQ.interp_fill_value )
      zetavar=dble(0)
    END WHERE
  ELSEIF(dimvar.EQ.2) THEN
    !weights for v-points
    ALLOCATE(hvar(nsize(1),nsize(2)-1))
    ALLOCATE(zetavar(nsize(1),nsize(2)-1))
    ALLOCATE(zw(nsize(1),nsize(2)-1,param%n_sigma))
    hvar=.5*h(:,1:nsize(2)-1)+.5*h(:,2:nsize(2))
    zetavar=.5*zeta(:,1:nsize(2)-1)+.5*zeta(:,2:nsize(2))
    WHERE( zeta(:,1:nsize(2)-1).EQ.interp_fill_value .OR. &
    & zeta(:,2:nsize(2)).EQ.interp_fill_value )
      zetavar=dble(0)
    END WHERE
  ELSE
    !weights for rho-points
    ALLOCATE(hvar(nsize(1),nsize(2)))
    ALLOCATE(zetavar(nsize(1),nsize(2)))
    ALLOCATE(zw(nsize(1),nsize(2),param%n_sigma))
    hvar=h
    zetavar=zeta
    WHERE( zeta(:,:).EQ.interp_fill_value )
      zetavar=dble(0)
    END WHERE
  END IF

  !calculate z-coordinate w-points
  ALLOCATE(sigma(param%n_sigma+1))
  DO i0=1,size(sigma,1)
   sigma(i0)=-dble(1)+dble(i0-1)/dble(param%n_sigma)
  END DO
  ALLOCATE(z(size(hvar,1),size(hvar,2),param%n_sigma+1))
  z=sigma2z(param,hvar,zetavar,sigma)

  !Calculate weights
  zw=z(:,:,2:size(z,3))-z(:,:,1:size(z,3)-1)
  zw=zw/SPREAD(SUM(zw,3),3,size(zw,3))
 
 END FUNCTION zweight

!--------------------------------------------------------------------------
 SUBROUTINE loc2index_1d(time,mask,timeo,it,w)
 !Finds for the point with coordinates timeo the indices of the two points
 !surrounding it with the weights for linear interpolation

 IMPLICIT NONE
 real(8),dimension(:),intent(in)::time,timeo
 logical,dimension(:),intent(in)::mask
 integer,dimension(:,:)::it
 real(8),dimension(:,:)::w
 integer::i1,i2
 real(8)::mint,maxt,sumw,wt(2)
  
 !min/max
 mint=MINVAL(time)
 maxt=MAXVAL(time)

 w=dble(0)

 DO i2=1,size(timeo)
  IF(timeo(i2).LT.mint.OR.timeo(i2).GT.maxt) THEN
   it(:,i2)=[1,1]
   w(:,i2)=[0.0,0.0]
  ELSE
   it(1,i2)=MAXLOC(time,1,time.LE.timeo(i2))
   it(2,i2)=MINLOC(time,1,time.GE.timeo(i2))
  
   IF(it(1,i2).EQ.it(2,i2)) THEN
    wt=[.5,.5]
   ELSE
    wt(1)=ABS( time(it(2,i2))-timeo(i2) )
    wt(2)=ABS( time(it(1,i2))-timeo(i2) )
   END IF
   
   DO i1=1,2
    IF(.NOT.mask(it(i1,i2))) wt(i1)=dble(0)
   END DO
   
   sumw=SUM(wt)
   IF(sumw.GT.dble(0)) THEN
    w(:,i2)=wt/sumw
   ELSE
    w(:,i2)=wt
   END IF

  END IF
 END DO
 END SUBROUTINE loc2index_1d


!--------------------------------------------------------------------------

 SUBROUTINE loc2index_2d(lon,lat,mask,lono,lato,ix,iy,w)
 !Finds for the points with x-coordinates lono and y-coordinates lato
 !the grid indices of the four points surrounding it and their weights for
 !bilinear interpolation. Here lon are x-coordinates of a rectangular grid,
 !lat the y-coordinates and mask the mask of the grid. 

  IMPLICIT NONE
  real(8),dimension(:),intent(in)::lon,lat
  logical,dimension(:,:),intent(in)::mask
  real(8),dimension(:),intent(in)::lono,lato
  integer,dimension(:,:)::ix,iy
  real(8),dimension(:,:)::w
  real(8)::lonmin,lonmax,latmin,latmax,sumw(2,2),wx1(2,2),wy1(2,2)  
  integer::i1,i2,j1,j2,ix1(2,2),iy1(2,2)

  lonmin=MINVAL(lon); lonmax=MAXVAL(lon)
  latmin=MINVAL(lat); latmax=MAXVAL(lat)
  
  w=dble(0) 
  DO i2=1,size(lono)
    IF(lono(i2).LE.lonmin.OR.lono(i2).GE.lonmax&
    .OR.lato(i2).LE.latmin.OR.lato(i2).GE.latmax) THEN
     ix(:,i2)=[1,1,1,1]
     iy(:,i2)=[1,1,1,1]
     w(:,i2)=dble(0)
    ELSE
     ix1(1,:)=MAXLOC(lon,1,lon.LE.lono(i2))
     ix1(2,:)=MINLOC(lon,1,lon.GE.lono(i2))
     iy1(:,1)=MAXLOC(lat,1,lat.LE.lato(i2))
     iy1(:,2)=MINLOC(lat,1,lat.GE.lato(i2))

     wx1=dble(0.5); wy1=dble(0.5)
     DO j2=1,2
     DO j1=1,2
      IF(ix1(1,1).NE.ix1(2,1)) THEN
       wx1(j1,j2)=ABS( lon(ix1(mod(j1,2)+1,j2))-lono(i2) )
      END IF
      IF(iy(1,1).NE.iy1(1,2)) THEN
       wy1(j1,j2)=ABS( lat(iy1(j1,mod(j2,2)+1))-lato(i2) )
      END IF
      IF(.NOT.mask(ix1(j1,j2),iy1(j1,j2))) THEN
       wx1(j1,j2)=dble(0); wy1(j1,j2)=dble(0) 
      END IF
     END DO
     END DO

     sumw=SPREAD(SUM(wx1,1),1,2)
     WHERE(sumw.NE.dble(0))
      wx1=wx1/sumw
     END WHERE
     sumw=SPREAD(SUM(wy1,2),2,2)
     WHERE(sumw.NE.dble(0))
      wy1=wy1/sumw
     END WHERE

     w(:,i2)=RESHAPE(wx1*wy1,[4])

    END IF
   END DO
  END SUBROUTINE loc2index_2d

  

END MODULE
