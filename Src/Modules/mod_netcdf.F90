MODULE mod_netcdf

IMPLICIT NONE

TYPE nf_att
 INTEGER::len
 CHARACTER(len=64)::name
 INTEGER::type
END TYPE nf_att

TYPE nf_dim
 INTEGER::len,id
 CHARACTER(len=64)::name
END TYPE nf_dim

TYPE nf_var
 INTEGER:: id
 CHARACTER(len=64)::name
 INTEGER::type,timedim,ndim
 REAL(8)::nan
END TYPE nf_var

INTERFACE get_nf_info
 MODULE PROCEDURE get_nf_info_global, get_nf_info_var
END INTERFACE
 
INTERFACE read_nf_field
 MODULE PROCEDURE read_nf_scalar,read_nf_1dfield,read_nf_2dfield,&
  &read_nf_3dfield
END INTERFACE

INTERFACE write_nf_field
 MODULE PROCEDURE write_nf_scalar,write_nf_1dfield,write_nf_2dfield,&
  &write_nf_3dfield
END INTERFACE

INTERFACE read_nf_points
 MODULE PROCEDURE read_nf_2dpoints,read_nf_3dpoints
END INTERFACE 

INTERFACE read_nf_att
 MODULE PROCEDURE read_nf_att_num, read_nf_att_text
END INTERFACE


INCLUDE 'netcdf.inc'

!----------------------------------------------------------------------------

CONTAINS

 SUBROUTINE ncsize(ncid,var_name,val)
  IMPLICIT NONE
  integer,intent(in)::ncid
  character(len=*),intent(in)::var_name
  integer::val(8)
  integer::status,varid,ndim,i1
  integer,allocatable::dimid(:)

  !Get dimids
  status=nf_inq_varid(ncid,TRIM(var_name),varid)
  status=nf_inq_varndims(ncid,varid,ndim)
  ALLOCATE(dimid(ndim))
  status=nf_inq_vardimid(ncid,varid,dimid)
  
  !Get size for each dimension
  val=-99
  ndim=MIN(ndim,size(val))
  DO i1=1,ndim
   status=nf_inq_dimlen(ncid,dimid(i1),val(i1))
  END DO
 END SUBROUTINE

 SUBROUTINE ncwrite0d(ncid,var_name,val)
  IMPLICIT NONE
  INTEGER,intent(in)::ncid
  CHARACTER(len=*),intent(in)::var_name
  REAL(8),intent(in)::val
  INTEGER::vartype,varid,status

  !Get varid
  status=nf_inq_varid(ncid,TRIM(var_name),varid)
  IF(status.NE.nf_noerr) THEN
   WRITE(*,*) 'ncwrite:cannot find variable ',TRIM(var_name)
   STOP
  END IF

  !Get data type
  status=nf_inq_vartype(ncid,varid,vartype)
  IF(status.NE.nf_noerr) THEN
   WRITE(*,*) 'ncwrite:cannot find type of variable ',TRIM(var_name)
   STOP
  END IF

  !Write data
  SELECT CASE(vartype)
   CASE (nf_int)
    status=nf_put_vara_int(ncid,varid,[1],[1],int(val))
   CASE (nf_float)
    status=nf_put_vara_real(ncid,varid,[1],[1],real(val))
   CASE (nf_double)
    status=nf_put_vara_double(ncid,varid,[1],[1],dble(val))
   CASE DEFAULT
    WRITE(*,*) 'ncwrite:writing type ',vartype,' not supported'
    STOP
  END SELECT
  IF(status.NE.nf_noerr) THEN
   WRITE(*,*) 'ncwrite:error writing variable ',TRIM(var_name)
   STOP
  END IF

 END SUBROUTINE ncwrite0d

!--------------------------------------------------------------------

 SUBROUTINE ncwrite1d(ncid,var_name,ncStart,ncCount,val)
  IMPLICIT NONE
  INTEGER,intent(in)::ncid,ncStart(1),ncCount(1)
  CHARACTER(len=*),intent(in)::var_name
  REAL(8),intent(in)::val(:)
  INTEGER::vartype,varid,status

  !Get varid
  status=nf_inq_varid(ncid,TRIM(var_name),varid)
  IF(status.NE.nf_noerr) THEN
   WRITE(*,*) 'ncwrite:cannot find variable ',TRIM(var_name)
   STOP
  END IF

  !Get data type
  status=nf_inq_vartype(ncid,varid,vartype)
  IF(status.NE.nf_noerr) THEN
   WRITE(*,*) 'ncwrite:cannot find type of variable ',TRIM(var_name)
   STOP
  END IF

  !Write data
  SELECT CASE(vartype)
   CASE (nf_int)
    status=nf_put_vara_int(ncid,varid,ncStart,ncCount,int(val))
   CASE (nf_float)
    status=nf_put_vara_real(ncid,varid,ncStart,ncCount,real(val))
   CASE (nf_double)
    status=nf_put_vara_double(ncid,varid,ncStart,ncCount,dble(val))
   CASE DEFAULT
    WRITE(*,*) 'ncwrite:writing type ',vartype,' not supported'
    STOP
  END SELECT
  IF(status.NE.nf_noerr) THEN
   WRITE(*,*) 'ncwrite:error writing variable ',TRIM(var_name)
   STOP
  END IF

 END SUBROUTINE ncwrite1d

!---------------------------------------------------------------------

 SUBROUTINE ncwrite2d(ncid,var_name,ncStart,ncCount,val)
  IMPLICIT NONE
  INTEGER,intent(in)::ncid,ncStart(2),ncCount(2)
  CHARACTER(len=*),intent(in)::var_name
  REAL(8),intent(in)::val(:,:)
  INTEGER::vartype,varid,status

  !Get varid
  status=nf_inq_varid(ncid,TRIM(var_name),varid)
  IF(status.NE.nf_noerr) THEN
   WRITE(*,*) 'ncwrite:cannot find variable ',TRIM(var_name)
   STOP
  END IF

  !Get data type
  status=nf_inq_vartype(ncid,varid,vartype)
  IF(status.NE.nf_noerr) THEN
   WRITE(*,*) 'ncwrite:cannot find type of variable ',TRIM(var_name)
   STOP
  END IF

  !Write data
  SELECT CASE(vartype)
   CASE (nf_int)
    status=nf_put_vara_int(ncid,varid,ncStart,ncCount,int(val))
   CASE (nf_float)
    status=nf_put_vara_real(ncid,varid,ncStart,ncCount,real(val))
   CASE (nf_double)
    status=nf_put_vara_double(ncid,varid,ncStart,ncCount,dble(val))
   CASE DEFAULT
    WRITE(*,*) 'ncwrite:writing type ',vartype,' not supported'
    STOP
  END SELECT
  IF(status.NE.nf_noerr) THEN
   WRITE(*,*) 'ncwrite:error writing variable ',TRIM(var_name)
   STOP
  END IF

 END SUBROUTINE ncwrite2d

!---------------------------------------------------------------------

 SUBROUTINE ncwrite3d(ncid,var_name,ncStart,ncCount,val)
  IMPLICIT NONE
  INTEGER,intent(in)::ncid,ncStart(3),ncCount(3)
  CHARACTER(len=*),intent(in)::var_name
  REAL(8),intent(in)::val(:,:,:)
  INTEGER::vartype,varid,status

  !Get varid
  status=nf_inq_varid(ncid,TRIM(var_name),varid)
  IF(status.NE.nf_noerr) THEN
   WRITE(*,*) 'ncwrite:cannot find variable ',TRIM(var_name)
   STOP
  END IF

  !Get data type
  status=nf_inq_vartype(ncid,varid,vartype)
  IF(status.NE.nf_noerr) THEN
   WRITE(*,*) 'ncwrite:cannot find type of variable ',TRIM(var_name)
   STOP
  END IF

  !Write data
  SELECT CASE(vartype)
   CASE (nf_int)
    status=nf_put_vara_int(ncid,varid,ncStart,ncCount,int(val))
   CASE (nf_float)
    status=nf_put_vara_real(ncid,varid,ncStart,ncCount,real(val))
   CASE (nf_double)
    status=nf_put_vara_double(ncid,varid,ncStart,ncCount,dble(val))
   CASE DEFAULT
    WRITE(*,*) 'ncwrite:writing type ',vartype,' not supported'
    STOP
  END SELECT
  IF(status.NE.nf_noerr) THEN
   WRITE(*,*) 'ncwrite:error writing variable ',TRIM(var_name)
   STOP
  END IF

 END SUBROUTINE ncwrite3d

!----------------------------------------------------------------------------

 SUBROUTINE ncwrite4d(ncid,var_name,ncStart,ncCount,val)
  IMPLICIT NONE
  INTEGER,intent(in)::ncid,ncStart(4),ncCount(4)
  CHARACTER(len=*),intent(in)::var_name
  REAL(8),intent(in)::val(:,:,:,:)
  INTEGER::vartype,varid,status

  !Get varid
  status=nf_inq_varid(ncid,TRIM(var_name),varid)
  IF(status.NE.nf_noerr) THEN
   WRITE(*,*) 'ncwrite:cannot find variable ',TRIM(var_name)
   STOP
  END IF

  !Get data type
  status=nf_inq_vartype(ncid,varid,vartype)
  IF(status.NE.nf_noerr) THEN
   WRITE(*,*) 'ncwrite:cannot find type of variable ',TRIM(var_name)
   STOP
  END IF

  !Write data
  SELECT CASE(vartype)
   CASE (nf_int)
    status=nf_put_vara_int(ncid,varid,ncStart,ncCount,int(val))
   CASE (nf_float)
    status=nf_put_vara_real(ncid,varid,ncStart,ncCount,real(val))
   CASE (nf_double)
    status=nf_put_vara_double(ncid,varid,ncStart,ncCount,dble(val))
   CASE DEFAULT
    WRITE(*,*) 'ncwrite:writing type ',vartype,' not supported'
    STOP
  END SELECT
  IF(status.NE.nf_noerr) THEN
   WRITE(*,*) 'ncwrite:error writing variable ',TRIM(var_name)
   STOP
  END IF

 END SUBROUTINE ncwrite4d

!--------------------------------------------------------------------------

 FUNCTION ncread0d(ncid,var_name) RESULT(val)
  IMPLICIT NONE
  INTEGER,intent(in)::ncid
  CHARACTER(len=*),intent(in)::var_name
  REAL(8)::val
  INTEGER::varid,status

  !Get varid
  status=nf_inq_varid(ncid,TRIM(var_name),varid)
  IF(status.NE.nf_noerr) THEN
   WRITE(*,*) 'ncread:cannot find ',TRIM(var_name)
   STOP
  END IF

  status=nf_get_vara_double(ncid,varid,[1],&
  &[1],val)
  IF(status.NE.nf_noerr) THEN
   WRITE(*,*) 'ncread: reading ',TRIM(var_name),' failed'
   STOP
  END IF
  
 END FUNCTION ncread0d

!-------------------------------------------------------------------------

 FUNCTION ncread1d(ncid,var_name,ncStart,ncCount) RESULT(val)
  IMPLICIT NONE
  INTEGER,intent(in)::ncid,ncStart(1),ncCount(1)
  CHARACTER(len=*),intent(in)::var_name
  REAL(8)::val(ncCount(1))
  INTEGER::varid,status

  !Get varid
  status=nf_inq_varid(ncid,TRIM(var_name),varid)
  IF(status.NE.nf_noerr) THEN
   WRITE(*,*) 'ncread:cannot find ',TRIM(var_name)
   STOP
  END IF

  status=nf_get_vara_double(ncid,varid,ncStart(:1),&
  &ncCount(:1),val)
  IF(status.NE.nf_noerr) THEN
   WRITE(*,*) 'ncread: reading ',TRIM(var_name),' failed'
   STOP
  END IF
  
 END FUNCTION ncread1d

!-------------------------------------------------------------------------

 FUNCTION ncread2d(ncid,var_name,ncStart,ncCount) RESULT(val)
  IMPLICIT NONE
  INTEGER,intent(in)::ncid,ncStart(2),ncCount(2)
  CHARACTER(len=*),intent(in)::var_name
  REAL(8)::val(ncCount(1),ncCount(2))
  INTEGER::varid,status

  !Get varid
  status=nf_inq_varid(ncid,TRIM(var_name),varid)
  IF(status.NE.nf_noerr) THEN
   WRITE(*,*) 'ncread:cannot find ',TRIM(var_name)
   STOP
  END IF

  status=nf_get_vara_double(ncid,varid,ncStart(:2),&
  &ncCount(:2),val)
  IF(status.NE.nf_noerr) THEN
   WRITE(*,*) 'ncread: reading ',TRIM(var_name),' failed'
   STOP
  END IF
  
 END FUNCTION ncread2d

!---------------------------------------------------------------------------

 FUNCTION ncread3d(ncid,var_name,ncStart,ncCount) RESULT(val)
  IMPLICIT NONE
  INTEGER,intent(in)::ncid,ncStart(3),ncCount(3)
  CHARACTER(len=*),intent(in)::var_name
  REAL(8)::val(ncCount(1),ncCount(2),ncCount(3))
  INTEGER::varid,status

  !Get varid
  status=nf_inq_varid(ncid,TRIM(var_name),varid)
  IF(status.NE.nf_noerr) THEN
   WRITE(*,*) 'ncread:cannot find ',TRIM(var_name)
   STOP
  END IF

  status=nf_get_vara_double(ncid,varid,ncStart,&
  &ncCount,val)
  IF(status.NE.nf_noerr) THEN
   WRITE(*,*) 'ncread: reading ',TRIM(var_name),' failed'
   STOP
  END IF
  
 END FUNCTION ncread3d

!------------------------------------------------------------------------------

 FUNCTION ncread4d(ncid,var_name,ncStart,ncCount) RESULT(val)
  IMPLICIT NONE
  INTEGER,intent(in)::ncid,ncStart(4),ncCount(4)
  CHARACTER(len=*),intent(in)::var_name
  REAL(8)::val(ncCount(1),ncCount(2),ncCount(3),ncCount(4))
  INTEGER::varid,status

  !Get varid
  status=nf_inq_varid(ncid,TRIM(var_name),varid)
  IF(status.NE.nf_noerr) THEN
   WRITE(*,*) 'ncread:cannot find ',TRIM(var_name)
   STOP
  END IF

  status=nf_get_vara_double(ncid,varid,ncStart(:4),&
  &ncCount(:4),val)
  IF(status.NE.nf_noerr) THEN
   WRITE(*,*) 'ncread: reading ',TRIM(var_name),' failed'
   STOP
  END IF
  
 END FUNCTION ncread4d

!--------------------------------------------------------------------------

 SUBROUTINE col_nfopen(dir,name,tBnd,ncid,time)
 !CALL dir_collect(dir,fname,ncid,time)
 !Open stream to files fname in directory dir and read their times
 
  IMPLICIT NONE
  CHARACTER(len=*),intent(in)::dir,name
  INTEGER::ncid(:)
  REAL(8),intent(in)::tBnd(2)
  REAL(8),allocatable,dimension(:)::time,time1
  CHARACTER(len=1024)::filename
  INTEGER::i0,i1,s(8),ibnd(2),fbnd(2),nt,status
  LOGICAL::flag_exist  

  !Open streams
  ncid=0
  DO i0=1,size(ncid)
   WRITE(filename,'(A,A,A,A,I0.4,A)') TRIM(dir),'/',TRIM(name),&
   &'_',i0,'.nc'
   INQUIRE(file=TRIM(filename),exist=flag_exist)   
   IF(flag_exist) status=nf_open(TRIM(filename),nf_nowrite,ncid(i0))
  END DO

  !Count number of times
  nt=0
  DO i0=1,size(ncid)
   IF(ncid(i0).EQ.0) CYCLE
   CALL ncsize(ncid(i0),'ocean_time',s)
   ALLOCATE(time1(s(1)))
   time1=ncread1d(ncid(i0),'ocean_time',[1],[s(1)])
   IF(.NOT.ANY(time1.GE.tBnd(1).AND.time1.LE.tBnd(2))) THEN
    status=nf_close(ncid(i0))
    ncid(i0)=0
   ELSE 
    nt=nt+COUNT(time1.GE.tBnd(1).AND.time1.LE.tBnd(2))
   END IF
   DEALLOCATE(time1)
  END DO

  !Read times
  ALLOCATE(time(nt)); ibnd(1)=1
  DO i0=1,size(ncid)
   IF(ncid(i0).EQ.0) CYCLE

   CALL ncsize(ncid(i0),'ocean_time',s)
   ALLOCATE(time1(s(1)))
   time1=ncread1d(ncid(i0),'ocean_time',[1],[s(1)])

   fbnd(1)=MINLOC(time1,1,time1.GE.tBnd(1).AND.time1.LE.tBnd(2))
   fbnd(2)=MAXLOC(time1,1,time1.GE.tBnd(1).AND.time1.LE.tBnd(2))
   ibnd(2)=ibnd(1)+fbnd(2)-fbnd(1)
  
   time(ibnd(1):ibnd(2))=time1(fbnd(1):fbnd(2))

   ibnd(1)=ibnd(2)+1
   DEALLOCATE(time1)
  END DO  
 
 END SUBROUTINE col_nfopen

!--------------------------------------------------------------------------

 FUNCTION col_ncread4d(ncid,varname,ncstart,nccount,tBnd) RESULT(val)
  IMPLICIT NONE
  INTEGER,intent(in)::ncid(:),ncstart(4),nccount(4)
  CHARACTER(len=*),intent(in)::varname
  REAL(8),intent(in)::tBnd(2)
  REAL(8)::val(ncstart(1):ncstart(1)-1+nccount(1),&
  &ncstart(2):ncstart(2)-1+nccount(2),&
  &ncstart(3):ncstart(3)-1+nccount(3),&
  &ncstart(4):ncstart(4)-1+nccount(4))
  INTEGER::i0,i1,ncend(4),s(8),fbnd1(2),fbnd2(3),ibnd1(2),ibnd2(2)
  INTEGER::status 
  REAL(8),allocatable::time1(:)
  INTEGER,allocatable::it1(:)
 
  !End indices
  ncend=ncstart+nccount-[1,1,1,1]
  
  ibnd1=1
  DO i0=1,size(ncid)
   IF(ncid(i0).EQ.0) CYCLE
   
   !Read times in file
   CALL ncsize(ncid(i0),'ocean_time',s)
   ALLOCATE(time1(s(1))); ALLOCATE(it1(s(1))); it1=0
   time1=ncread1d(ncid(i0),'ocean_time',[1],[s(1)])
   fbnd1(1)=MINLOC(time1,1,time1.GE.tbnd(1).AND.time1.LE.tbnd(2))
   fbnd1(2)=MAXLOC(time1,1,time1.GE.tbnd(1).AND.time1.LE.tbnd(2))
   ibnd1(2)=ibnd1(1)+fbnd1(2)-fbnd1(1)
   DO i1=fbnd1(1),fbnd1(2)
    it1(i1)=ibnd1(1)+(i1-fbnd1(1))
   END DO

   !Read values
   IF(ANY(it1.GE.ncstart(4).AND.it1.LE.ncend(4))) THEN
    fbnd2(1)=MINLOC(time1,1,it1.GE.ncstart(4).AND.it1.LE.ncend(4))
    fbnd2(2)=MAXLOC(time1,1,it1.GE.ncstart(4).AND.it1.LE.ncend(4))
    fbnd2(3)=fbnd2(2)-fbnd2(1)+1
    ibnd2(1)=it1(fbnd2(1)); ibnd2(2)=it1(fbnd2(2))

    val(:,:,:,ibnd2(1):ibnd2(2))=ncread4d(ncid(i0),varname,&
    &[ncstart(1),ncstart(2),ncstart(3),fbnd2(1)],&
    &[nccount(1),nccount(2),nccount(3),fbnd2(3)])
   END IF   

   ibnd1(1)=ibnd1(2)+1
   DEALLOCATE(it1,time1)
  END DO
 END FUNCTION col_ncread4d

 FUNCTION col_ncread3d(ncid,varname,ncstart,nccount,tBnd) RESULT(val)
  IMPLICIT NONE
  INTEGER,intent(in)::ncid(:),ncstart(3),nccount(3)
  CHARACTER(len=*),intent(in)::varname
  REAL(8),intent(in)::tBnd(2)
  REAL(8)::val(ncstart(1):ncstart(1)-1+nccount(1),&
  &ncstart(2):ncstart(2)-1+nccount(2),&
  &ncstart(3):ncstart(3)-1+nccount(3))
  INTEGER::i0,i1,ncend(3),s(8),fbnd1(2),fbnd2(3),ibnd1(2),ibnd2(2)
  INTEGER::status 
  REAL(8),allocatable::time1(:)
  INTEGER,allocatable::it1(:)

  !End indices
  ncend=ncstart+nccount-[1,1,1]
  
  ibnd1=1
  DO i0=1,size(ncid)
   IF(ncid(i0).EQ.0) CYCLE
   
   !Read times in file
   CALL ncsize(ncid(i0),'ocean_time',s)
   ALLOCATE(time1(s(1))); ALLOCATE(it1(s(1))); it1=0
   time1=ncread1d(ncid(i0),'ocean_time',[1],[s(1)])
   fbnd1(1)=MINLOC(time1,1,time1.GE.tbnd(1).AND.time1.LE.tbnd(2))
   fbnd1(2)=MAXLOC(time1,1,time1.GE.tbnd(1).AND.time1.LE.tbnd(2))
   ibnd1(2)=ibnd1(1)+fbnd1(2)-fbnd1(1)
   DO i1=fbnd1(1),fbnd1(2)
    it1(i1)=ibnd1(1)+(i1-fbnd1(1))
   END DO

   !Read values
   IF(ANY(it1.GE.ncstart(3).AND.it1.LE.ncend(3))) THEN
    fbnd2(1)=MINLOC(time1,1,it1.GE.ncstart(3).AND.it1.LE.ncend(3))
    fbnd2(2)=MAXLOC(time1,1,it1.GE.ncstart(3).AND.it1.LE.ncend(3))
    fbnd2(3)=fbnd2(2)-fbnd2(1)+1
    ibnd2(1)=it1(fbnd2(1)); ibnd2(2)=it1(fbnd2(2))

    val(:,:,ibnd2(1):ibnd2(2))=ncread3d(ncid(i0),varname,&
    &[ncstart(1),ncstart(2),fbnd2(1)],&
    &[nccount(1),nccount(2),fbnd2(3)])
   END IF   

   ibnd1(1)=ibnd1(2)+1
   DEALLOCATE(it1,time1)
  END DO
 END FUNCTION col_ncread3d


!---------------------------------------------------------------------------

 FUNCTION ncvarsize(ncid,var_name) RESULT(var_size)
  
  IMPLICIT NONE
  INTEGER,intent(in)::ncid
  CHARACTER(len=*),intent(in)::var_name
  INTEGER::var_size(4)
  INTEGER::i0,ndims,status,varid
  INTEGER,allocatable::dimid(:)

  var_size=0

  !Get varid
  status=nf_inq_varid(ncid,TRIM(var_name),varid)
  IF(status.NE.nf_noerr) THEN
   WRITE(*,*) 'ncvarsize:cannot find ',TRIM(var_name)
   STOP
  END IF

  !Get dimensions
  status=nf_inq_varndims(ncid,varid,ndims)
  IF(ndims.GT.0) THEN
   ALLOCATE(dimid(ndims))
   status=nf_inq_vardimid(ncid,varid,dimid)

   !Get size of dimensions
   DO i0=1,MIN(ndims,4)
    status=nf_inq_dimlen(ncid,dimid(i0),var_size(i0))
   END DO    

  END IF

 END FUNCTION ncvarsize


 SUBROUTINE read_grd(ncid,var_names,x,y,mask)
 !READ_GRD Read grid from netcdf file
 !
 !Syntax:
 ! CALL read_grd(ncid,'x_name y_name mask_name',x,y,mask)
 !Input:
 ! ncid: INTEGER created by nf_open
 ! x_name: fieldname in netcdf file of x-coordinates
 ! y_name: fielname in netcdf file of y-coordinates
 ! mask_name: fieldname in netcdf file of mask
 !Output:
 ! x: ALL DOUBLE 2D array with x-coordinates
 ! y: ALL DOUBLE 2D array with y-coordinates
 ! mask: ALL INTEGER 2D array with mask

  IMPLICIT NONE 

  INTEGER,INTENT(in)::ncid
  CHARACTER(len=*),INTENT(in)::var_names
  CHARACTER(len=64)::var_name(3)
  REAL(8),ALLOCATABLE,INTENT(inout)::x(:,:),y(:,:)
  INTEGER,ALLOCATABLE,INTENT(inout)::mask(:,:)
  INTEGER::var_id(3),dim_id(2,3),dim_len(2)
  INTEGER::status,i1
  CHARACTER(len=1024)::err_str


  read(var_names,*) var_name(1),var_name(2),var_name(3)
  !Find variable id of variables
  DO i1=1,3
  status=nf_inq_varid(ncid,TRIM(var_name(i1)),var_id(i1))
  IF(status.NE.nf_noerr) THEN
   WRITE(err_str,*) 'mod_netcdf.read_grd: variable ',&
    &TRIM(var_name(i1)),' not present in file.'
   WRITE(*,*) TRIM(err_str)
   STOP
  END IF
  status=nf_inq_vardimid(ncid,var_id(i1),dim_id(:,i1))
  END DO
  
  !Read size output
  IF((ANY(dim_id(:,1).NE.dim_id(:,2))).OR.&
   &(ANY(dim_id(:,1).NE.dim_id(:,3)))) THEN 
   WRITE(*,*) 'mod_interp.read_grd: dimensions variables are &
   &incongruent.'
   STOP
  END IF
  status=nf_inq_dimlen(ncid,dim_id(1,1),dim_len(1))
  status=nf_inq_dimlen(ncid,dim_id(2,1),dim_len(2))
  IF(ALLOCATED(x)) THEN
   IF(ANY(SHAPE(x).NE.dim_len)) DEALLOCATE(x)
  END IF
  IF(.NOT.ALLOCATED(x)) ALLOCATE(x(dim_len(1),dim_len(2)))
  IF(ALLOCATED(y)) THEN
   IF(ANY(SHAPE(Y).NE.dim_len)) DEALLOCATE(y)
  END IF
  IF(.NOT.ALLOCATED(y)) ALLOCATE(y(dim_len(1),dim_len(2)))
  IF(ALLOCATED(mask)) THEN
   IF(ANY(SHAPE(mask).NE.dim_len)) DEALLOCATE(mask)
  END IF
  IF(.NOT.ALLOCATED(mask)) ALLOCATE(mask(dim_len(1),dim_len(2)))

  !Read coordinates
  status=nf_get_vara_double(ncid,var_id(1),[1,1],&
   dim_len,x)
  status=nf_get_vara_double(ncid,var_id(2),[1,1],&
   dim_len,y)
  status=nf_get_vara_int(ncid,var_id(3),[1,1],&
   dim_len,mask)

 END SUBROUTINE read_grd

!----------------------------------------------------------------------------
 
 SUBROUTINE get_nf_info_global(ncid,var,dim,att)
 !GET_NF_INFO_GLOBAL Read info about all variables, dimensions and global 
 !attributes in a netcdf file
 !
 !Syntax:
 ! CALL get_nf_info(ncid,var,dim,att)
 !Input:
 ! ncid: INTEGER created with nf_open
 !Output: 
 ! var:  ALL. TYPE(nf_var) array with metadata on variables 
 ! dim: ALL. TYPE(nf_dim) array with metadata on dimensions
 ! att: ALL. TYPE(nf_att) array with metadata on global attributes
  
 TYPE(nf_var),INTENT(inout),ALLOCATABLE::var(:)
 TYPE(nf_dim),INTENT(inout),ALLOCATABLE::dim(:)
 TYPE(nf_att),INTENT(inout),ALLOCATABLE::att(:)
 INTEGER::n_dim,n_var,n_att,n_unlim,var_id,att_type,att_len
 INTEGER::status,i1,i2
 INTEGER,ALLOCATABLE::dim_id(:)
 CHARACTER*64::att_name
 INTEGER,intent(in)::ncid


 !Read info netcdf file
 status=nf_inq(ncid,n_dim,n_var,n_att,n_unlim)
 IF(status.NE.nf_noerr) THEN
  WRITE(*,*) 'mod_nf_netcdf.get_nf_info_global: cannot read netcdf&
  & file.'
  STOP
 END IF

 !Read dimensions
 IF (n_dim.GT.0) THEN
  IF(ALLOCATED(dim)) THEN
   IF (SIZE(dim).NE.n_dim) DEALLOCATE(dim)
  END IF
  IF(.NOT.ALLOCATED(dim))ALLOCATE(dim(n_dim))
 
  DO i1=1,n_dim
   dim(i1)%id=i1
   status=nf_inq_dimname(ncid,i1,dim(i1)%name)
   status=nf_inq_dimlen(ncid,i1,dim(i1)%len)
  END DO

 END IF !dim

 !Read variable information
 IF (n_var.GT.0) THEN
  IF (ALLOCATED(var)) THEN
   IF (SIZE(var).NE.n_var) DEALLOCATE(var)
  END IF
  IF (.NOT.ALLOCATED(var)) ALLOCATE(var(n_var))
  
  DO i1=1,n_var
   !info
   var(i1)%id=i1
   status=nf_inq_varname(ncid,i1,var(i1)%name)
   status=nf_inq_vartype(ncid,i1,var(i1)%type)
   status=nf_inq_varndims(ncid,i1,var(i1)%ndim)

   !time dimension
   status=nf_get_att_text(ncid,i1,'time',att_name)
   status=nf_inq_attlen(ncid,i1,'time',att_len)
   IF (status.NE.nf_noerr.OR.n_dim.EQ.0) THEN
    var(i1)%timedim=0
   ELSE
    ALLOCATE(dim_id(n_dim))
    status=nf_inq_vardimid(ncid,i1,dim_id)
    DO i2=1,n_dim
     IF (strcmp(att_name(1:att_len),dim(i2)%name)) THEN
      var(i1)%timedim=find(dim_id.EQ.i2,'first')
     END IF
    END DO
    DEALLOCATE(dim_id)
   END IF

   !nan value
   !!status=nf_inq_atttype(ncid,i1,'_FillValue',att_type)
   !!IF (status.NE.nf_noerr) THEN
   !! var(i1)%nan=DBLE(1e37)
   !!ELSE
   !! status=nf_get_att_double(ncid,i1,'_FillValue',var(i1)%nan)
   !!END IF
  END DO

  END IF !var

  !Get global attributes
  status=nf_inq_varnatts(ncid,nf_global,n_att)
  IF (n_att.GT.0) THEN

  IF (ALLOCATED(att)) THEN
   IF (SIZE(att).NE.n_att) DEALLOCATE(att)
  END IF
  IF (.NOT.ALLOCATED(att)) ALLOCATE(att(n_att))

  DO i1=1,n_att
   status=nf_inq_attname(ncid,nf_global,i1,att_name)
   att(i1)%name=TRIM(att_name)
   status=nf_inq_att(ncid,nf_global,att_name,att(i1)%type,att(i1)%len)   
  END DO
 END IF !att

 END SUBROUTINE get_nf_info_global

!---------------------------------------------------------------------------

 SUBROUTINE get_nf_info_var(ncid,var_name,var,dim,att)
 !GET_NF_INFO_VAR Read info on 1 variable and its dimensions and attributes
 !
 !Syntax:
 ! CALL get_nf_info(ncid,var,dim,att)
 !Input:
 ! ncid: INTEGER created with nf_open
 ! var_name: CHARACTER containing the name of the variable
 !Output: 
 ! var: TYPE(nf_var) with metadata on variables 
 ! dim: ALL. TYPE(nf_dim) array with metadata on dimensions
 ! att: ALL. TYPE(nf_att) array with metadata on global attributes
 
  INTEGER,INTENT(in)::ncid
  CHARACTER(len=*),INTENT(in)::var_name
  CHARACTER(len=1024)::err_str
  CHARACTER*64::att_name
  TYPE(nf_var),INTENT(out)::var
  TYPE(nf_dim),ALLOCATABLE,INTENT(out),OPTIONAL::dim(:)
  TYPE(nf_att),ALLOCATABLE,INTENT(out),OPTIONAL::att(:)
  TYPE(nf_dim),ALLOCATABLE::dim_loc(:)
  TYPE(nf_att),ALLOCATABLE::att_loc(:)
  INTEGER,ALLOCATABLE::dim_id(:)
  INTEGER::var_id,n_att,att_type,att_len
  INTEGER::status,i1,i2

  !Find var_id
  status=nf_inq_varid(ncid,TRIM(var_name),var_id)
  IF(status.NE.nf_noerr) THEN
   WRITE(err_str,*) 'mod_netcdf.get_nf_info_var: cannot find &
    &variable ',TRIM(var_name),' in netcdf file.'
   WRITE(*,*) TRIM(err_str)
   STOP
  END IF

  !Read var info
  var%id=var_id
  status=nf_inq_varname(ncid,var_id,var%name)
  status=nf_inq_vartype(ncid,var_id,var%type)
  status=nf_inq_varndims(ncid,var_id,var%ndim)

  !Get dimensions
  IF(var%ndim.GT.0) THEN
   ALLOCATE(dim_id(var%ndim))
   ALLOCATE(dim_loc(var%ndim))

   status=nf_inq_vardimid(ncid,var_id,dim_id)
   DO i1=1,var%ndim
    dim_loc(i1)%id=dim_id(i1)
    status=nf_inq_dimname(ncid,dim_id(i1),dim_loc(i1)%name)
    status=nf_inq_dimlen(ncid,dim_id(i1),dim_loc(i1)%len)
   END DO
  END IF !dim

  !Get attributes
  status=nf_inq_varnatts(ncid,var_id,n_att)
  IF (n_att.GT.0) THEN
   ALLOCATE(att_loc(n_att))

   DO i1=1,n_att
    status=nf_inq_attname(ncid,var_id,i1,att_name)
    att_loc(i1)%name=TRIM(att_name)
    status=nf_inq_att(ncid,var_id,att_name,&
     &att_loc(i1)%type,att_loc(i1)%len)   
   END DO
  END IF !att

  !time dimension
  status=nf_get_att_text(ncid,var_id,'time',att_name)
  status=nf_inq_attlen(ncid,var_id,'time',att_len)
  IF (status.NE.nf_noerr.OR.var%ndim.EQ.0) THEN
   var%timedim=0
  ELSE
   var%timedim=-1
   DO i2=1,var%ndim
    IF (strcmp(att_name(1:att_len),dim_loc(i2)%name)) THEN
     var%timedim=i2
    END IF
   END DO
  END IF
  
  
  !nan value
  !!status=nf_inq_atttype(ncid,var_id,'_FillValue',att_type)
  !!IF (status.NE.nf_noerr) THEN
  !! var%nan=DBLE(1e37)
  !!ELSE
  !! status=nf_get_att_double(ncid,var_id,'_FillValue',var%nan)
  !!END IF

  !Copy to optional arguments
  IF(ALLOCATED(dim_id)) DEALLOCATE(dim_id) 
  IF(PRESENT(dim)) THEN
   IF(ALLOCATED(dim)) DEALLOCATE(dim)
   IF(.NOT.ALLOCATED(dim)) ALLOCATE(dim(var%ndim))
   dim=dim_loc
  END IF
  IF(ALLOCATED(dim_loc)) DEALLOCATE(dim_loc)
  IF(PRESENT(att)) THEN
   IF(ALLOCATED(att)) DEALLOCATE(att)
   IF(.NOT.ALLOCATED(att)) ALLOCATE(att(n_att))
   att=att_loc
  END IF
  IF(ALLOCATED(att_loc)) DEALLOCATE(att_loc)
  
 END SUBROUTINE get_nf_info_var

!---------------------------------------------------------------------------

SUBROUTINE read_nf_scalar(ncid,var_name,time,val)
 !READ_NF_SCALAR Read scalar value from netcdf file
 !
 !Syntax:
 ! CALL read_nf_scalar(ncid,var_name,time,val)
 !Input:
 ! ncid: INTEGER generated with nf_open
 ! var_name: CHARACTER with name of scalar in netcdf file
 ! time: INTEGER irrelevant for scalar
 !Output:
 ! val: double with value of scalar

 IMPLICIT NONE

 INTEGER,INTENT(in)::ncid
 INTEGER,INTENT(in)::time
 CHARACTER(len=*),INTENT(in)::var_name
 CHARACTER(len=1024)::err_str
 REAL(8),INTENT(out)::val
 INTEGER::status
 TYPE(nf_var)::var
  
 !Get info
 CALL get_nf_info(ncid,TRIM(var_name),var)

 !Read values
 status=nf_get_vara_double(ncid,var%id,1,1,val)
 IF(status.NE.nf_noerr) THEN
  WRITE(err_str,*) 'mod_netcdf.read_nf_scalar: reading variable '&
   &,trim(var_name),' unsuccessfull.'
  WRITE(*,*) TRIM(err_str)
  STOP
 END IF

 END SUBROUTINE read_nf_scalar 

!----------------------------------------------------------------------------

 SUBROUTINE write_nf_scalar(ncid,var_name,time,val)
 !WRITE_NF_SCALAR Write scalar value to netcdf file
 !
 !Syntax:
 ! CALL write_nf_scalar(ncid,var_name,time,val)
 !Input:
 ! ncid: INTEGER generated with nf_open
 ! var_name: CHARACTER with name of scalar in netcdf file
 ! time: INTEGER irrelevant for scalar
 !Output:
 ! val: double with value of scalar

 IMPLICIT NONE

 INTEGER,INTENT(in)::ncid
 INTEGER,INTENT(in)::time
 CHARACTER(len=*),INTENT(in)::var_name
 TYPE(nf_var)::var
 REAL(8),INTENT(in)::val
 INTEGER::status
 CHARACTER(len=1024)::err_str

 !Get variable info
 CALL get_nf_info(ncid,TRIM(var_name),var)

 !Write values
 IF (var%type.EQ.nf_int) THEN
  status=nf_put_vara_int(ncid,var%id,[1],[1],NINT(val)) 
 ELSEIF (var%type.EQ.nf_double) THEN
  status=nf_put_vara_double(ncid,var%id,1,1,DBLE(val)) 
 ELSEIF (var%type.EQ.nf_float) THEN
  status=nf_put_vara_real(ncid,var%id,1,1,REAL(val)) 
 ELSE
  WRITE(*,*) "mod_netcdf.write_nf_scalar: data type not supported."
  STOP
 END IF
 IF(status.NE.nf_noerr) THEN
  WRITE(err_str,*) 'mod_netcdf.write_nf_scalar: writing variable '&
   &,trim(var_name),' unsuccessfull.'
  WRITE(*,*) TRIM(err_str)
  STOP
 END IF
 
 END SUBROUTINE write_nf_scalar

!----------------------------------------------------------------------------

SUBROUTINE read_nf_1dfield(ncid,var_name,time,val,timedim)
 !READ_NF_1DFIELD Read scalar value from netcdf file
 !
 !Syntax:
 ! CALL read_nf_1dfield(ncid,var_name,time,val,timedim)
 !Input:
 ! ncid: INTEGER generated with nf_open
 ! var_name: CHARACTER with name of scalar in netcdf file
 ! time: INTEGER time index of the field to be read
 ! timedim: OPT INTEGER giving the dimension in which time is stored. The 
 !          routine tries to find this index itself by looking for the time
 !          attribute 
 !Output:
 ! val: ALL 1D DOUBLE array where read values can be stored. 

 IMPLICIT NONE

 INTEGER,INTENT(in)::ncid
 INTEGER,INTENT(in)::time
 INTEGER,INTENT(in),OPTIONAL::timedim
 CHARACTER(len=*),INTENT(in)::var_name
 CHARACTER(len=1024)::err_str
 REAL(8),ALLOCATABLE,INTENT(inout)::val(:)
 INTEGER::status,timedim_loc,i1,i2
 INTEGER,ALLOCATABLE::nf_start(:),nf_count(:),nf_shape(:)
 TYPE(nf_var)::var
 TYPE(nf_dim),ALLOCATABLE::dim(:)

 !Get info
 CALL get_nf_info(ncid,TRIM(var_name),var,dim)


 SELECT CASE(var%ndim)
  CASE(1)
   !No time dependence
   timedim_loc=2

  !Range
   ALLOCATE(nf_start(1)); nf_start=[1]
   ALLOCATE(nf_count(1)); nf_count=dim(1)%len
  CASE(2)
   !Time dependence
   IF(PRESENT(timedim)) THEN
    timedim_loc=timedim
   ELSE
    timedim_loc=var%timedim
   END IF

   !Allocate
   ALLOCATE(nf_start(2)); nf_start=[1,1]
   ALLOCATE(nf_count(2)); nf_count=[dim(1)%len,dim(2)%len]
   
   !Set time dimension
   IF(timedim_loc.LT.1.OR.timedim_loc.GT.2) THEN
     WRITE(err_str,*) 'mod_netcdf.read_nf_1dfield: cannot find time&
     & dimension for variable ',TRIM(var_name)
    WRITE(*,*) TRIM(err_str)
     STOP
   END IF
   nf_start(timedim_loc)=time
   nf_count(timedim_loc)=1
    
  CASE DEFAULT
   WRITE(err_str,*) 'mod_netcdf.write_nf_1dfield: ',TRIM(var_name)&
    &,' is not a 1D field.'
   WRITE(*,*) TRIM(err_str)
   STOP
 END SELECT

 !Allocate val
 ALLOCATE(nf_shape(1))
 i2=1;
 DO i1=1,SIZE(nf_count)
  IF(i1.NE.timedim_loc) THEN
   nf_shape(i2)=nf_count(i1); i2=i2+1
  END IF
 END DO
 IF(ALLOCATED(val)) THEN
  IF(ANY(SHAPE(val).NE.nf_shape)) DEALLOCATE(val)
 END IF
 IF(.NOT.ALLOCATED(val)) ALLOCATE(val(nf_shape(1)))  

 !Read values
 status=nf_get_vara_double(ncid,var%id,nf_start,nf_count,val)
 IF(status.NE.nf_noerr) THEN
  WRITE(err_str,*) 'mod_netcdf.read_nf_1dfield: reading variable '&
   &,trim(var_name),' unsuccessfull.'
  WRITE(*,*) TRIM(err_str)
  STOP
 END IF

 END SUBROUTINE read_nf_1dfield

!------------------------------------------------------------------------------

SUBROUTINE write_nf_1dfield(ncid,var_name,time,val,timedim)
 !WRITE__NF_1DFIELD Write scalar value to netcdf file
 !
 !Syntax:
 ! CALL read_nf_1dfield(ncid,var_name,time,val,timedim)
 !Input:
 ! ncid: INTEGER generated with nf_open
 ! var_name: CHARACTER with name of scalar in netcdf file
 ! time: INTEGER time index of the field to be read
 ! timedim: OPT INTEGER giving the dimension in which time is stored. The 
 !          routine tries to find this index itself by looking for the time
 !          attribute 
 !Output:
 ! val: ALL 1D DOUBLE array with values to be written


 IMPLICIT NONE

 INTEGER,INTENT(in)::ncid
 INTEGER,INTENT(in)::time
 INTEGER,INTENT(in),OPTIONAL::timedim
 CHARACTER(len=*),INTENT(in)::var_name
 CHARACTER(len=1024)::err_str
 REAL(8),ALLOCATABLE,INTENT(in)::val(:)
 INTEGER::status,timedim_loc,i1,i2
 INTEGER,ALLOCATABLE::nf_start(:),nf_count(:),nf_shape(:)
 TYPE(nf_var)::var
 TYPE(nf_dim),ALLOCATABLE::dim(:)

 !Get info
 CALL get_nf_info(ncid,TRIM(var_name),var,dim)


 SELECT CASE(var%ndim)
  CASE(1)
   !No time dependence
   timedim_loc=2

  !Range
   ALLOCATE(nf_start(1)); nf_start=[1]
   ALLOCATE(nf_count(1)); nf_count=dim(1)%len
  CASE(2)
   !Time dependence
   IF(PRESENT(timedim)) THEN
    timedim_loc=timedim
   ELSE
    timedim_loc=var%timedim
   END IF

   !Allocate
   ALLOCATE(nf_start(2)); nf_start=[1,1]
   ALLOCATE(nf_count(2)); nf_count=[dim(1)%len,dim(2)%len]
   
   !Set time dimension
   IF(timedim_loc.LT.1.OR.timedim_loc.GT.2) THEN
    WRITE(err_str,*) 'mod_netcdf.write_nf_1dfield: cannot find time&
     & dimension for variable ',TRIM(var_name)
    WRITE(*,*) TRIM(err_str)
    STOP
   END IF
   nf_start(timedim_loc)=time
   nf_count(timedim_loc)=1
    
  CASE DEFAULT
   WRITE(err_str,*) 'mod_netcdf.write_nf_1dfield: ',TRIM(var_name)&
    &,' is not a 1D field.'
   WRITE(*,*) TRIM(err_str)
   STOP
 END SELECT

 !Allocate nf_shape
 ALLOCATE(nf_shape(1))
 i2=1;
 DO i1=1,SIZE(nf_count)
  IF(i1.NE.timedim_loc) THEN
   nf_shape(i2)=nf_count(i1); i2=i2+1
  END IF
 END DO
 IF(.NOT.ALLOCATED(val)) THEN
 WRITE(err_str,*) 'mod_netcdf.write_nf_1dfield: no input provided.'
  WRITE(*,*) TRIM(err_str)
  STOP
 END IF  

 !Write values
 IF (var%type.EQ.nf_int) THEN
  status=nf_put_vara_int(ncid,var%id,nf_start,nf_count,NINT(val)) 
 ELSEIF (var%type.EQ.nf_double) THEN
  status=nf_put_vara_double(ncid,var%id,nf_start,nf_count,DBLE(val)) 
 ELSEIF (var%type.EQ.nf_float) THEN
  status=nf_put_vara_real(ncid,var%id,nf_start,nf_count,REAL(val)) 
 ELSE
  WRITE(*,*) "mod_netcdf.write_nf_1dfield: data type not supported."
  STOP
 END IF
 IF(status.NE.nf_noerr) THEN
  WRITE(err_str,*) 'mod_netcdf.write_nf_1dfield: writing&
  & variable ',trim(var_name),' unsuccessfull.'
  WRITE(*,*) TRIM(err_str)
  STOP
 END IF

 END SUBROUTINE write_nf_1dfield

!-------------------------------------------------------------------------

SUBROUTINE read_nf_2dfield(ncid,var_name,time,val,timedim)
 !READ_NF_2DFIELD Read scalar value from netcdf file
 !
 !Syntax:
 ! CALL read_nf_2dfield(ncid,var_name,time,val,timedim)
 !Input:
 ! ncid: INTEGER generated with nf_open
 ! var_name: CHARACTER with name of scalar in netcdf file
 ! time: INTEGER time index of the field to be read
 ! timedim: OPT INTEGER giving the dimension in which time is stored. The 
 !          routine tries to find this index itself by looking for the time
 !          attribute 
 !Output:
 ! val: 2D ALL DOUBLE array for storing read values

 IMPLICIT NONE

 INTEGER,INTENT(in)::ncid
 INTEGER,INTENT(in)::time
 INTEGER,INTENT(in),OPTIONAL::timedim
 CHARACTER(len=*),INTENT(in)::var_name
 CHARACTER(len=1024)::err_str
 REAL(8),ALLOCATABLE,INTENT(inout)::val(:,:)
 INTEGER::status,timedim_loc,i1,i2
 INTEGER,ALLOCATABLE::nf_start(:),nf_count(:),nf_shape(:)
 TYPE(nf_var)::var
 TYPE(nf_dim),ALLOCATABLE::dim(:)

 !Get info
 CALL get_nf_info(ncid,TRIM(var_name),var,dim)


 SELECT CASE(var%ndim)
  CASE(2)
   !No time dependence
   timedim_loc=3

  !Range
   ALLOCATE(nf_start(2)); nf_start=[1,1]
   ALLOCATE(nf_count(2)); nf_count=[dim(1)%len,dim(2)%len]
  CASE(3)
   !Time dependence
   IF(PRESENT(timedim)) THEN
    timedim_loc=timedim
   ELSE
    timedim_loc=var%timedim
   END IF

   !Allocate
   ALLOCATE(nf_start(3)); nf_start=[1,1,1]
   ALLOCATE(nf_count(3)); nf_count=[dim(1)%len,dim(2)%len,dim(3)%len]
   
   !Set time dimension
   IF(timedim_loc.LT.1.OR.timedim_loc.GT.3) THEN
    WRITE(err_str,*) 'mod_netcdf.read_nf_2dfield: cannot find time&
     & dimension for variable ',TRIM(var_name)
    WRITE(*,*) TRIM(err_str)
    STOP
   END IF
   nf_start(timedim_loc)=time
   nf_count(timedim_loc)=1
    
  CASE DEFAULT
   WRITE(err_str,*) 'mod_netcdf.read_nf_2dfield: ',TRIM(var_name)&
    &,' is not a 2D field.'
   WRITE(*,*) TRIM(err_str)
   STOP
 END SELECT

 !Allocate val
 ALLOCATE(nf_shape(2))
 i2=1;
 DO i1=1,SIZE(nf_count)
  IF(i1.NE.timedim_loc) THEN
   nf_shape(i2)=nf_count(i1); i2=i2+1
  END IF
 END DO
 IF(ALLOCATED(val)) THEN
  IF(ANY(SHAPE(val).NE.nf_shape)) DEALLOCATE(val)
 END IF
 IF(.NOT.ALLOCATED(val)) ALLOCATE(val(nf_shape(1),nf_shape(2)))  

 !Read values

 status=nf_get_vara_double(ncid,var%id,nf_start,nf_count,val)
 IF(status.NE.nf_noerr) THEN
  WRITE(err_str,*) 'mod_netcdf.read_nf_2dfield: reading variable '&
   &,trim(var_name),' unsuccessfull.'
  WRITE(*,*)TRIM(err_str)
  STOP
 END IF

 END SUBROUTINE read_nf_2dfield

!----------------------------------------------------------------------------

SUBROUTINE write_nf_2dfield(ncid,var_name,time,val,timedim)
 !WRITE__NF_2DFIELD Write scalar value to netcdf file
 !
 !Syntax:
 ! CALL read_nf_1dfield(ncid,var_name,time,val,timedim)
 !Input:
 ! ncid: INTEGER generated with nf_open
 ! var_name: CHARACTER with name of scalar in netcdf file
 ! time: INTEGER time index of the field to be read
 ! timedim: OPT INTEGER giving the dimension in which time is stored. The 
 !          routine tries to find this index itself by looking for the time
 !          attribute 
 !Output:
 ! val: ALL 2D DOUBLE array with values to be written


 IMPLICIT NONE

 INTEGER,INTENT(in)::ncid
 INTEGER,INTENT(in)::time
 INTEGER,INTENT(in),OPTIONAL::timedim
 CHARACTER(len=*),INTENT(in)::var_name
 CHARACTER(len=1024)::err_str
 REAL(8),INTENT(in)::val(:,:)
 INTEGER::status,timedim_loc,i1,i2
 INTEGER,ALLOCATABLE::nf_start(:),nf_count(:),nf_shape(:)
 TYPE(nf_var)::var
 TYPE(nf_dim),ALLOCATABLE::dim(:)

 !Get info
 CALL get_nf_info(ncid,TRIM(var_name),var,dim)

 SELECT CASE(var%ndim)
  CASE(2)
   !No time dependence
   timedim_loc=3

  !Range
   ALLOCATE(nf_start(2)); nf_start=[1,1]
   ALLOCATE(nf_count(2)); nf_count=[dim(1)%len,dim(2)%len]
  CASE(3)
   !Time dependence
   IF(PRESENT(timedim)) THEN
    timedim_loc=timedim
   ELSE
    timedim_loc=var%timedim
   END IF

   !Allocate
   ALLOCATE(nf_start(3)); nf_start=[1,1,1]
   ALLOCATE(nf_count(3)); nf_count=[dim(1)%len,dim(2)%len,dim(3)%len]
   
   !Set time dimension
   IF(timedim_loc.LT.1.OR.timedim_loc.GT.3) THEN
    WRITE(err_str,*) 'mod_netcdf.write_nf_2dfield: cannot find time&
     & dimension for variable ',TRIM(var_name)
    WRITE(*,*) TRIM(err_str)
    STOP
   END IF
   nf_start(timedim_loc)=time
   nf_count(timedim_loc)=1
    
  CASE DEFAULT
   WRITE(err_str,*) 'mod_netcdf.write_nf_2dfield: ',TRIM(var_name)&
    &,' is not a 2D field.'
   WRITE(*,*) TRIM(err_str)
   STOP
 END SELECT

 !Allocate nf_shape
 ALLOCATE(nf_shape(2))
 i2=1;
 DO i1=1,SIZE(nf_count)
  IF(i1.NE.timedim_loc) THEN
   nf_shape(i2)=nf_count(i1); i2=i2+1
  END IF
 END DO 

 !Write values
 IF (var%type.EQ.nf_int) THEN
  status=nf_put_vara_int(ncid,var%id,nf_start,nf_count,NINT(val)) 
 ELSEIF (var%type.EQ.nf_double) THEN
  status=nf_put_vara_double(ncid,var%id,nf_start,nf_count,DBLE(val)) 
 ELSEIF (var%type.EQ.nf_float) THEN
  status=nf_put_vara_real(ncid,var%id,nf_start,nf_count,REAL(val)) 
 ELSE
  WRITE(*,*) "mod_netcdf.write_nf_2dfield: data type not supported."
  STOP
 END IF
 IF(status.NE.nf_noerr) THEN
  WRITE(err_str,*) 'mod_netcdf.write_nf_2dfield: writing&
  & variable ',trim(var_name),' unsuccessfull.'
  WRITE(*,*) TRIM(err_str)
  STOP
 END IF

 END SUBROUTINE write_nf_2dfield

!-----------------------------------------------------------------------------

SUBROUTINE read_nf_3dfield(ncid,var_name,time,val,timedim)
 !READ_NF_3DFIELD Read scalar value from netcdf file
 !
 !Syntax:
 ! CALL read_nf_3dfield(ncid,var_name,time,val,timedim)
 !Input:
 ! ncid: INTEGER generated with nf_open
 ! var_name: CHARACTER with name of scalar in netcdf file
 ! time: INTEGER time index of the field to be read
 ! timedim: OPT INTEGER giving the dimension in which time is stored. The 
 !          routine tries to find this index itself by looking for the time
 !          attribute 
 !Output:
 ! val: 3D ALL DOUBLE array for storing read values

 IMPLICIT NONE

 INTEGER,INTENT(in)::ncid
 INTEGER,INTENT(in)::time
 INTEGER,INTENT(in),OPTIONAL::timedim
 CHARACTER(len=*),INTENT(in)::var_name
 CHARACTER(len=1024)::err_str
 REAL(8),ALLOCATABLE,INTENT(inout)::val(:,:,:)
 INTEGER::status,timedim_loc,i1,i2
 INTEGER,ALLOCATABLE::nf_start(:),nf_count(:),nf_shape(:)
 TYPE(nf_var)::var
 TYPE(nf_dim),ALLOCATABLE::dim(:)

 !Get info
 CALL get_nf_info(ncid,TRIM(var_name),var,dim)

 SELECT CASE(var%ndim)
  CASE(3)
   !No time dependence
   timedim_loc=4

  !Range
   ALLOCATE(nf_start(3)); nf_start=[1,1,1]
   ALLOCATE(nf_count(3)); nf_count=[dim(1)%len,dim(2)%len,dim(3)%len]
  CASE(4)
   !Time dependence
   IF(PRESENT(timedim)) THEN
    timedim_loc=timedim
   ELSE
    timedim_loc=var%timedim
   END IF

   !Allocate
   ALLOCATE(nf_start(4)); nf_start=[1,1,1,1]
   ALLOCATE(nf_count(4))
   nf_count=[dim(1)%len,dim(2)%len,dim(3)%len,dim(4)%len]
   
   !Set time dimension
   IF(timedim_loc.LT.1.OR.timedim_loc.GT.4) THEN
    WRITE(err_str,*) 'mod_netcdf.read_nf_3dfield: cannot find time&
     & dimension for variable ',TRIM(var_name)
    WRITE(*,*) TRIM(err_str)
    STOP
   END IF
   nf_start(timedim_loc)=time
   nf_count(timedim_loc)=1
    
  CASE DEFAULT
   WRITE(err_str,*) 'mod_netcdf.read_nf_3dfield: ',TRIM(var_name)&
    &,' is not a 3D field.'
   WRITE(*,*) TRIM(err_str)
   STOP
 END SELECT

 !Allocate val
 ALLOCATE(nf_shape(3))
 i2=1;
 DO i1=1,SIZE(nf_count)
  IF(i1.NE.timedim_loc) THEN
   nf_shape(i2)=nf_count(i1); i2=i2+1
  END IF
 END DO
 IF(ALLOCATED(val)) THEN
  IF(ANY(SHAPE(val).NE.nf_shape)) DEALLOCATE(val)
 END IF
 IF(.NOT.ALLOCATED(val)) &
  &ALLOCATE(val(nf_shape(1),nf_shape(2),nf_shape(3)))  

 !Read values
 status=nf_get_vara_double(ncid,var%id,nf_start,nf_count,val)
 IF(status.NE.nf_noerr) THEN
  WRITE(err_str,*) 'mod_netcdf.read_nf_3dfield: reading variable '&
   &,trim(var_name),' unsuccessfull.'
  WRITE(*,*) TRIM(err_str)
  STOP
 END IF

 END SUBROUTINE read_nf_3dfield

!----------------------------------------------------------------------------

SUBROUTINE write_nf_3dfield(ncid,var_name,time,val,timedim)
 !WRITE__NF_3DFIELD Write scalar value to netcdf file
 !
 !Syntax:
 ! CALL read_nf_3dfield(ncid,var_name,time,val,timedim)
 !Input:
 ! ncid: INTEGER generated with nf_open
 ! var_name: CHARACTER with name of scalar in netcdf file
 ! time: INTEGER time index of the field to be read
 ! timedim: OPT INTEGER giving the dimension in which time is stored. The 
 !          routine tries to find this index itself by looking for the time
 !          attribute 
 !Output:
 ! val: ALL 3D DOUBLE array with values to be written


 IMPLICIT NONE

 INTEGER,INTENT(in)::ncid
 INTEGER,INTENT(in)::time
 INTEGER,INTENT(in),OPTIONAL::timedim
 CHARACTER(len=*),INTENT(in)::var_name
 CHARACTER(len=1024)::err_str
 REAL(8),INTENT(in)::val(:,:,:)
 INTEGER::status,timedim_loc,i1,i2
 INTEGER,ALLOCATABLE::nf_start(:),nf_count(:),nf_shape(:)
 TYPE(nf_var)::var
 TYPE(nf_dim),ALLOCATABLE::dim(:)

 !Get info
 CALL get_nf_info(ncid,TRIM(var_name),var,dim)
 
 SELECT CASE(var%ndim)
  CASE(3)
   !No time dependence
   timedim_loc=4

  !Range
   ALLOCATE(nf_start(3)); nf_start=[1,1,1]
   ALLOCATE(nf_count(3))
   nf_count=[dim(1)%len,dim(2)%len,dim(3)%len]
  CASE(4)
   !Time dependence
   IF(PRESENT(timedim)) THEN
    timedim_loc=timedim
   ELSE
    timedim_loc=var%timedim
   END IF

   !Allocate
   ALLOCATE(nf_start(4)); nf_start=[1,1,1,1]
   ALLOCATE(nf_count(4))
   nf_count=[dim(1)%len,dim(2)%len,dim(3)%len,dim(4)%len]
   
   !Set time dimension
   IF(timedim_loc.LT.1.OR.timedim_loc.GT.4) THEN
    WRITE(err_str,*) 'mod_netcdf.write_nf_3dfield: cannot find time&
     & dimension for variable ',TRIM(var_name)
    WRITE(*,*) TRIM(err_str)
    STOP
   END IF
   nf_start(timedim_loc)=time
   nf_count(timedim_loc)=1
    
  CASE DEFAULT
   WRITE(err_str,*) 'mod_netcdf.write_nf_3dfield: ',TRIM(var_name)&
    &,' is not a 3D field.'
   WRITE(*,*) TRIM(err_str)
   STOP
 END SELECT

 !Allocate nf_shape
 ALLOCATE(nf_shape(4))
 i2=1;
 DO i1=1,SIZE(nf_count)
  IF(i1.NE.timedim_loc) THEN
   nf_shape(i2)=nf_count(i1); i2=i2+1
  END IF
 END DO

 !Write values
 IF (var%type.EQ.nf_int) THEN
  status=nf_put_vara_int(ncid,var%id,nf_start,nf_count,NINT(val)) 
 ELSEIF (var%type.EQ.nf_double) THEN   
  status=nf_put_vara_double(ncid,var%id,nf_start,nf_count,DBLE(val)) 
 ELSEIF (var%type.EQ.nf_float) THEN
  status=nf_put_vara_real(ncid,var%id,nf_start,nf_count,REAL(val)) 
 ELSE
  WRITE(*,*) "mod_netcdf.write_nf_3dfield: data type not supported."
  STOP
 END IF
 IF(status.NE.nf_noerr) THEN
  WRITE(err_str,*) 'mod_netcdf.write_nf_3dfield: writing&
  & variable ',trim(var_name),' unsuccessfull.'
  WRITE(*,*) TRIM(err_str)
  STOP
 END IF

 END SUBROUTINE write_nf_3dfield

!----------------------------------------------------------------------------
 
 SUBROUTINE read_nf_2dpoints(ncid,var_name,time,i1,i2,val,timedim)
 !READ_NF_2DPOINTS read values at point (i1,i2) from 2D-field
 !
 !Syntax:
 ! CALL read_nf_point(ncid,var_name,time,i1,i2,val,timedim)
 !Input:
 ! ncid: INTEGER generated with nf_open
 ! var_name: CHARACTER with name of scalar in netcdf file
 ! i1: 1D INTEGER array with indices first horizontal coordinate
 ! i2: 1D INTEGER array with indices second horizontal coordinate
 ! time: 1D INTEGER time index of the field to be read
 ! timedim: OPT INTEGER giving the dimension in which time is stored. The 
 !          routine tries to find this index itself by looking for the time
 !          attribute 
 !Output:
 ! val: 1D DOUBLE array with values at the given horizontal indices


  IMPLICIT NONE
  
  !Declare variables
  INTEGER,INTENT(in)::ncid
  INTEGER,INTENT(in),OPTIONAL::timedim
  INTEGER,INTENT(in)::i1(:),i2(:),time(:)
  REAL(8),INTENT(inout),ALLOCATABLE::val(:)
  CHARACTER(len=*),INTENT(in)::var_name
  
  INTEGER::loc(3),i0,status
  INTEGER,ALLOCATABLE::nf_start(:),nf_count(:)
  LOGICAL::flag_filled(lbound(i1,1):ubound(i1,1))
  CHARACTER(len=1024)::err_str
  REAL(8)::val_tmp
  TYPE(nf_var)::var
  TYPE(nf_dim),ALLOCATABLE::dim(:)

  IF(SIZE(i1,1).NE.SIZE(i2,1)) THEN
   WRITE(err_str,*) 'mod_netcdf.write_nf_2dpoints: length i1 and& 
   &must be equal.'
   STOP
  END IF
 
  !Get info
  CALL get_nf_info(ncid,TRIM(var_name),var,dim)

  SELECT CASE(var%ndim)
   CASE(2)
   !No time dependence
   ALLOCATE(nf_start(2))
   ALLOCATE(nf_count(2))
   nf_count=[1,1]
   loc=[1,2,3]   
  CASE(3)
   !Time dimension 
   IF(PRESENT(timedim)) THEN
    loc(3)=timedim
   ELSE
    loc(3)=var%timedim
   END IF
   DO i0=1,3
    IF(i0.LT.loc(3)) loc(i0)=i0
    IF(i0.GT.loc(3)) loc(i0)=i0-1
   END DO

   ALLOCATE(nf_start(3))
   ALLOCATE(nf_count(3))
   nf_count=[1,1,1]    
  CASE DEFAULT
   WRITE(err_str,*) 'mod_netcdf.write_nf_2dpoints: ',TRIM(var_name)&
    &,' is not a 2D field.'
   WRITE(*,*) TRIM(err_str)
   STOP
  END SELECT

  !Read netcdf file
  flag_filled=.FALSE.
  ALLOCATE(val(lbound(i1,1):ubound(i1,1)))
  DO i0=lbound(i1,1),ubound(i1,1)
   IF(.NOT.flag_filled(i0)) THEN
    nf_start(loc(1))=i1(i0)
    nf_start(loc(2))=i2(i0)
    nf_start(loc(3))=time(i0)
    status=nf_get_vara_double(ncid,var%id,nf_start,nf_count,val_tmp)
    !Fill all entries with the same indices at the same time 
    WHERE( i1.EQ.i1(i0).AND.i2.EQ.i2(i0).AND.time.EQ.time(i0) )
     val=val_tmp
     flag_filled=.TRUE.
    END WHERE
   END IF
  END DO

 END SUBROUTINE read_nf_2dpoints

!----------------------------------------------------------------------------

SUBROUTINE read_nf_3dpoints(ncid,var_name,time,i1,i2,i3,val,timedim)
!READ_NF_3DPOINT read values at point (i1,i2) from 3D-field
!
!Syntax:
! CALL read_nf_points(ncid,var_name,time,i1,i2,i3,val,timedim)
!Input:
! ncid: INTEGER generated with nf_open
! var_name: CHARACTER with name of scalar in netcdf file
! i1: 1D INTEGER array with indices first horizontal coordinate
! i2: 1D INTEGER array with indices second horizontal coordinate
! i3: 1D INTEGER array with indices s-coordinate
! time: 1D INTEGER time index of the field to be read
! timedim: OPT INTEGER giving the dimension in which time is stored. The 
!          routine tries to find this index itself by looking for the time
!          attribute 
!Output:
! val: 2D DOUBLE array with values at the given horizontal indices along the 1st dimension
!      and in the vertical along the 2nd dimension 

  IMPLICIT NONE
  
  !Declare variables
  INTEGER,INTENT(in)::ncid
  INTEGER,INTENT(in),OPTIONAL::timedim
  INTEGER,INTENT(in)::i1(:),i2(:),i3(:),time(:)
  REAL(8),INTENT(inout),ALLOCATABLE::val(:)
  CHARACTER(len=*),INTENT(in)::var_name
  
  INTEGER::loc(4),i0,i00,status
  INTEGER,ALLOCATABLE::nf_start(:),nf_count(:)
  LOGICAL::flag_filled(lbound(i1,1):ubound(i1,1))
  CHARACTER(len=1024)::err_str
  REAL(8)::val_tmp
  TYPE(nf_var)::var
  TYPE(nf_dim),ALLOCATABLE::dim(:)

  IF(SIZE(i1,1).NE.SIZE(i2,1)) THEN
   WRITE(err_str,*) 'mod_netcdf.read_nf_3dpoints: length i1 and& 
   &must be equal.'
   STOP
  END IF
 
  !Get info
  CALL get_nf_info(ncid,TRIM(var_name),var,dim)


  SELECT CASE(var%ndim)
   CASE(3)
   !No time dependence
   ALLOCATE(nf_start(3)); 
   ALLOCATE(nf_count(3))
   nf_count=[1,1,1]
   loc=[1,2,3,4]   
  CASE(4)
   !Time dimension 
   IF(PRESENT(timedim)) THEN
    loc(4)=timedim
   ELSE
    loc(4)=var%timedim
   END IF
   DO i0=1,4
    IF(i0.LT.loc(4)) loc(i0)=i0
    IF(i0.GT.loc(4)) loc(i0)=i0-1
   END DO

   ALLOCATE(nf_start(4))
   ALLOCATE(nf_count(4))
   nf_count=[1,1,1,1]    
  CASE DEFAULT
   WRITE(err_str,*) 'mod_netcdf.read_nf_3dpoints: ',TRIM(var_name)&
    &,' is not a 3D field.'
   WRITE(*,*) TRIM(err_str)
   STOP
  END SELECT
 
  !Read netcdf file
  flag_filled=.FALSE.
  DO i0=lbound(i1,1),ubound(i1,1)
   IF(.NOT.flag_filled(i0)) THEN
    nf_start(loc(1))=i1(i0)
    nf_start(loc(2))=i2(i0)
    nf_start(loc(3))=i3(i0)
    nf_start(loc(4))=time(i0)
    status=nf_get_vara_double(ncid,var%id,nf_start,nf_count,val_tmp)
    !Fill all entries with the same indices at the same time 
    DO i00=lbound(i1,1),ubound(i1,1)
     IF(flag_filled(i00)) CYCLE 
     IF( i1(i00).EQ.i1(i0).AND.i2(i00).EQ.i2(i0).AND.&
     &i3(i0).EQ.i3(i00).AND.time(i0).EQ.time(i00) ) THEN
      val(i00)=val_tmp
      flag_filled(i00)=.TRUE.
     END IF
    END DO
   END IF
  END DO

 END SUBROUTINE read_nf_3dpoints

!----------------------------------------------------------------------------

SUBROUTINE add_dim(ncid,dim)
 !ADD_DIM Add a dimension to a netcdf file
 !
 !Syntax:
 ! CALL add_dim(dim,file_out)
 !Input:
 ! dim: TYPE(nf_dim) struct with information about the dimension
 ! ncid: INTEGER generated with nf_open. Must be set to nf_write
  
  IMPLICIT NONE
  
  TYPE(nf_dim),INTENT(in)::dim
  INTEGER,INTENT(in)::ncid
  INTEGER::status,dim_id
  CHARACTER(len=1024)::err_str

  !Enable define mode
  status=nf_redef(ncid)
  IF(status.NE.nf_noerr) THEN
   WRITE(err_str,*) 'mod_netcdf.add_dim: cannot redefine netcdf&
    & file.'
   WRITE(*,*) TRIM(err_str)
   STOP
  END IF

  !Check if dimension already exists and add if not
  status=nf_inq_dimid(ncid,TRIM(dim%name),dim_id)
  IF(status.NE.nf_noerr) THEN
   !Add dimension
   status=nf_def_dim(ncid,TRIM(dim%name),dim%len,dim_id)
   IF(status.NE.nf_noerr) THEN
    WRITE(err_str,*) 'mod_netcdf.add_dim: unable to add &
    &dimension ',TRIM(dim%name)
    WRITE(*,*) TRIM(err_str)
    STOP
   ELSE
   END IF
  ELSE
   WRITE(*,*) 'mod_netcdf.add_dim: dimension ',TRIM(dim%name),&
    &' does already exist.'
  END IF
   
  !End redefine mode
  status=nf_enddef(ncid) 

 END SUBROUTINE add_dim

!-----------------------------------------------------------------------------

SUBROUTINE add_var(ncid,var,dim)
 !ADD_VAR add variable to netcdf file
 !
 !Syntax:
 ! CALL add_var(var,file_out)
 !Input:
 ! var: TYPE(nf_var) struct with information about the variable
 ! ncid: INTEGER generated with nf_open. Must be in nf_write mode. 

  IMPLICIT NONE
 
  TYPE(nf_var),INTENT(in)::var
  TYPE(nf_dim),ALLOCATABLE,INTENT(in)::dim(:)
  INTEGER,INTENT(in)::ncid
  INTEGER::status,i1,i2,dim_id,var_id
  INTEGER::dim_ids(SIZE(dim))
  CHARACTER(len=1024)::err_str
 
 
  !Check if all the necessary dimensions are present
  DO i1=1,SIZE(dim)
   status=nf_inq_dimid(ncid,TRIM(dim(i1)%name),dim_id)
   IF(status.NE.nf_noerr) THEN
    CALL add_dim(ncid,dim(i1))
    status=nf_inq_dimid(ncid,TRIM(dim(i1)%name),dim_id)
   END IF 
   dim_ids(i1)=dim_id
  END DO 

  !Check if variable is not already defined
  status=nf_inq_varid(ncid,TRIM(var%name),var_id)
  IF (status.EQ.nf_noerr) THEN 
   WRITE(*,*) 'mod_netcdf.add_var: variable ',&
     &TRIM(var%name),' has already been defined.'
   var_id=-1
  END IF

  !Enable define mode
  status=nf_redef(ncid)
  IF(status.NE.nf_noerr) THEN
   WRITE(err_str,*) 'mod_netcdf.add_var: cannot redefine netcdf&
    & file.'
   WRITE(*,*) TRIM(err_str)
   STOP
  END IF    

  !Define variable
  IF(var_id.NE.-1) THEN
   status=nf_def_var(ncid,TRIM(var%name),var%type,SIZE(dim_ids),&
    &dim_ids,var_id)
   IF(status.NE.nf_noerr) THEN
    WRITE(err_str,*) 'mod_netcdf.add_var: unable to add &
    &variable ',TRIM(var%name)
    WRITE(*,*) TRIM(err_str)
    STOP
   END IF
  END IF

  !Close netcdf file
  status=nf_enddef(ncid)

 END SUBROUTINE add_var

!----------------------------------------------------------------------------

 SUBROUTINE copy_att(ncid_in,ncid_out,var_in,var_out)
 !COPY_ATT Copy attributes from 1 netcdf file to the other or from 
 !different variables in 1 netcdf file
 ! 
 !Syntax:
 ! CALL copy_att(ncid_in,ncid_out,var_in,var_out)
 !Input:
 ! ncid_in: INTEGER generated with nf_open representing the stream to source
 !          file.
 ! ncid_out: INTEGER generated with nf_open respresenting the stream to the 
 !           target file.
 ! var_in: OPT CHARACTER name of the variable from which the attributes must 
 !         be copied. If var_name='global' only global attributes are copied. 
 !         If neither var_in nor var_out are specified the attributes of all
 !         variables are copied. 
 ! var_out: OPT CHARACTER name of the variable to which the attributes must
 !          be copied. 

  
  INTEGER,INTENT(in)::ncid_in,ncid_out
  CHARACTER(len=*),INTENT(in),OPTIONAL::var_in,var_out
  INTEGER::status,i1,i2,n_att,n_var,var_id1,var_id2
  CHARACTER(len=1024)::err_str
  CHARACTER(len=64)::att_name,var_name
  
 
  !Put output in define mode
  status=nf_redef(ncid_out)
  IF(status.NE.nf_noerr) THEN
   WRITE(err_str,*) 'mod_netcdf.copy_att: cannot redefine netcdf&
    & file.'
   WRITE(*,*) TRIM(err_str)
   STOP
  END IF    

  IF(.NOT.PRESENT(var_in)) THEN
   !Copy attributes of all variables

    !Information on input netcdf file

    !Global attributes
    status=nf_inq_varnatts(ncid_in,nf_global,n_att)
    IF(status.NE.nf_noerr) n_att=0
    DO i1=1,n_att
     status=nf_inq_attname(ncid_in,nf_global,i1,att_name)
     status=nf_copy_att(ncid_in,nf_global,TRIM(att_name),&
      &ncid_out,nf_global)
    END DO

    !Attributes per variable
    status=nf_inq_nvars(ncid_in,n_var)
    DO i1=1,n_var
     status=nf_inq_varname(ncid_in,i1,var_name)
     status=nf_inq_varid(ncid_out,TRIM(var_name),var_id2)
     IF(status.EQ.nf_noerr) THEN
      status=nf_inq_varnatts(ncid_in,i1,n_att)
      DO i2=1,n_att
       status=nf_inq_attname(ncid_in,i1,i2,att_name)
       status=nf_copy_att(ncid_in,i1,TRIM(att_name),&
        &ncid_out,var_id2)
      END DO
     END IF      
    END DO 
    
  ELSE
   !Copy attributes from 1 variable

   status=nf_inq_varid(ncid_in,TRIM(var_in),var_id1)
   IF(status.NE.nf_noerr) THEN
    WRITE(err_str,*) 'mod_netcdf.copy_att: variable ',TRIM(var_in),&
     &' not present.'
    WRITE(*,*) err_str
    STOP
   END IF
   status=nf_inq_varnatts(ncid_in,var_id1,n_att)
   status=nf_inq_varid(ncid_out,TRIM(var_out),var_id2)
   DO i2=1,n_att
    status=nf_inq_attname(ncid_in,var_id1,i2,att_name)
    status=nf_copy_att(ncid_in,var_id1,TRIM(att_name),&
       &ncid_out,var_id2)
    IF(status.NE.nf_noerr) THEN
     WRITE(err_str,*) 'mod_netcdf.copy_att: unable to copy attribute',&
      &TRIM(att_name)
     WRITE(*,*) TRIM(err_str)
     STOP
    END IF
   END DO  

  END IF

  !Close output netcdf file
  status=nf_enddef(ncid_out)
  
 END SUBROUTINE

 !--------------------------------------------------------------------------

 SUBROUTINE copy_dim(ncid_in,ncid_out)
 !COPY_DIM Copy dimensions from 1 netcdf file to the other
 !
 !Syntax:
 ! CALL copy_dim(ncid_in,ncid_out
 !Input:
 ! ncid_in: INTEGER generated with nf_open representing the stream to the
 !          source file
 ! ncid_out: INTEGER generated with nf_open representing the stream to the
 !           target file

 INTEGER,INTENT(in)::ncid_in,ncid_out
 CHARACTER(len=64)::dim_name
 CHARACTER(len=1024)::err_str
 INTEGER::status,i1,dim_len,dim_id,n_dim
 
 status=nf_inq_ndims(ncid_in,n_dim)
 IF(n_dim.GT.0) THEN
  
  !Open define mode
  status=nf_redef(ncid_out)
  IF(status.NE.nf_noerr) THEN
   WRITE(*,*) 'mod_netcdf.copy_dim: cannot open redefine mode.'
   STOP
  END IF

  !Get metadata on dimension
  DO i1=1,n_dim
   status=nf_inq_dim(ncid_in,i1,dim_name,dim_len)
   status=nf_def_dim(ncid_out,TRIM(dim_name),dim_len,dim_id)
   IF(status.NE.nf_noerr) THEN
    WRITE(err_str,*) 'mod_netcdf.copy_dim: cannot copy dimension ',&
     &TRIM(dim_name)
     WRITE(*,*) TRIM(err_str)
     STOP
   END IF
  END DO
  
  status=nf_enddef(ncid_out)

 END IF

 END SUBROUTINE copy_dim

!-----------------------------------------------------------------------------

 SUBROUTINE get_var_shape(ncid,var_name,var_shape,opt)
 !GET_VAR_SHAPE Get the shape of the array of a variable
 !
 !Syntax:
 ! CALL get_var_shape(ncid,var_name,var_shape,opt)
 !Input:
 ! ncid: INTEGER generated with nf_open representing stream to netcdf file
 ! var_name: CHARACTER with name of field in netcdf file
 ! opt: OPT CHARACTER. If opt='f' only the size of the non-time dimensions will
 !      be given. If opt='t' only the size of the time dimension will be given.
 !Output:
 ! var_shape: INTEGER array with the size of the different dimensions. 
 
  INTEGER,INTENT(in)::ncid
  CHARACTER(len=*),INTENT(in)::var_name
  CHARACTER(len=1),INTENT(in),OPTIONAL::opt
  INTEGER,ALLOCATABLE,INTENT(out)::var_shape(:)
  INTEGER::status,n_dim,var_id,time_id,i1,i2,time_name_len
  INTEGER,ALLOCATABLE::dim_id(:)
  CHARACTER(len=1024)::err_str
  CHARACTER(len=64)::time_name,dim_name
  
  !Get var id
  status=nf_inq_varid(ncid,TRIM(var_name),var_id)
  IF(status.NE.nf_noerr) THEN
   WRITE(err_str,*) 'mod_netcdf.get_var_shape: cannot find variable '&
    &,TRIM(var_name)
   WRITE(*,*) TRIM(err_str)
   STOP
  END IF

  !Get dimensions
  status=nf_inq_varndims(ncid,var_id,n_dim)
  ALLOCATE(dim_id(n_dim))
  status=nf_inq_vardimid(ncid,var_id,dim_id)

  !Find time dimension
  time_id=-1
  status=nf_get_att_text(ncid,var_id,'time',time_name)
  status=nf_inq_attlen(ncid,var_id,'time',time_name_len)
  time_name=time_name(1:time_name_len)
  !write(*,*) 'time_name:',time_name,att_len
  IF(status.EQ.nf_noerr) THEN
   DO i1=1,n_dim
    status=nf_inq_dimname(ncid,dim_id(i1),dim_name)
    IF(strcmp(dim_name,time_name)) time_id=dim_id(i1)
   END DO
  END IF
   
  !Get size of each dimension
  IF(.NOT.PRESENT(opt)) THEN
   ALLOCATE(var_shape(n_dim))
   DO i1=1,SIZE(dim_id)
    status=nf_inq_dimlen(ncid,dim_id(i1),var_shape(i1))
   END DO    
  ELSEIF(opt.EQ.'t') THEN
   ALLOCATE(var_shape(1))
   status=nf_inq_dimlen(ncid,time_id,var_shape(1))
  ELSEIF(opt.EQ.'f') THEN
   IF(time_id.NE.-1) THEN
    ALLOCATE(var_shape(n_dim-1))
   ELSE
    ALLOCATE(var_shape(n_dim))
   END IF
   i2=1
   DO i1=1,n_dim
    IF(dim_id(i1).NE.time_id)  THEN
     status=nf_inq_dimlen(ncid,dim_id(i1),var_shape(i2))
     i2=i2+1
    END IF
   END DO
  END IF

 END SUBROUTINE get_var_shape

!-----------------------------------------------------------------------------

 SUBROUTINE read_nf_att_num(ncid,var_name,att_name,att_val)
 !READ_NF_ATT

 IMPLICIT NONE

 INTEGER,INTENT(in)::ncid
 CHARACTER(len=*),INTENT(in)::var_name,att_name
 REAL(8),ALLOCATABLE,DIMENSION(:),INTENT(out)::att_val
 CHARACTER(len=1024)::err_str
 INTEGER::status,var_id,att_len

 !Find varid
 IF(strcmp(var_name,'global')) THEN
  var_id=nf_global
 ELSE
  status=nf_inq_varid(ncid,TRIM(var_name),var_id)
  IF(status.NE.0) THEN
   WRITE(err_str,*) 'mod_netcdf.read_nf_att_num: cannot find &
    &variable ',TRIM(var_name)
   WRITE(*,*) TRIM(err_str)
   STOP
  END IF
 END IF

 !Find length attribute
 status=nf_inq_attlen(ncid,var_id,TRIM(att_name),att_len)
 IF(att_len.GE.1) THEN
  ALLOCATE(att_val(att_len))
 ELSE
  WRITE(err_str,*) 'mod_netcdf.read_nf_att_num: cannot establish &
    &attribute length.'
  WRITE(*,*) TRIM(err_str)
  STOP
 END IF
 
 !Get values attribute 
 status=nf_get_att_double(ncid,var_id,TRIM(att_name),att_val)
 IF(status.NE.0) THEN
  WRITE(err_str,*) 'mod_interp.read_nf_att: error reading attribute ',&
   &TRIM(att_name),' of variable ',TRIM(var_name)
  WRITE(*,*) TRIM(err_str)
  STOP
 END IF

 END SUBROUTINE read_nf_att_num

!--------------------------------------------------------------------------

SUBROUTINE read_nf_att_text(ncid,var_name,att_name,att_val)
 !READ_NF_ATT_TEXT

 IMPLICIT NONE

 INTEGER,INTENT(in)::ncid
 CHARACTER(len=*),INTENT(in)::var_name,att_name
 CHARACTER(len=*),INTENT(out)::att_val
 CHARACTER(len=1024)::err_str
 INTEGER::status,var_id,att_len

 !Find varid
 IF(strcmp(var_name,'global')) THEN
  var_id=nf_global
 ELSE
  status=nf_inq_varid(ncid,TRIM(var_name),var_id)
  IF(status.NE.0) THEN
   WRITE(err_str,*) 'mod_netcdf.read_nf_att_num: cannot find &
    &variable ',TRIM(var_name)
   WRITE(*,*) TRIM(err_str)
   STOP
  END IF
 END IF

 !Find length attribute
 status=nf_inq_attlen(ncid,var_id,TRIM(att_name),att_len)
 IF(att_len.GE.1) THEN
  IF(LEN(att_val).LT.att_len) THEN
   WRITE(err_str,*) 'mod_netcdf.read_nf_att_text: length string att_val &
    &is insufficient.'
   WRITE(*,*) TRIM(err_str)
   STOP
  END IF
 ELSE
  WRITE(err_str,*) 'mod_netcdf.read_nf_att_text: cannot establish &
    &attribute length.'
  WRITE(*,*) TRIM(err_str)
  STOP
 END IF
 
 !Get values attribute 
 status=nf_get_att_text(ncid,var_id,TRIM(att_name),att_val)
 IF(status.NE.0) THEN
  WRITE(err_str,*) 'mod_interp.read_nf_att_text: error reading &
   &attribute ',TRIM(att_name),' of variable ',TRIM(var_name)
  WRITE(*,*) TRIM(err_str)
  STOP
 END IF
 att_val=att_val(1:att_len)

 END SUBROUTINE read_nf_att_text

!-----------------------------------------------------------------------------

 ELEMENTAL FUNCTION strcmp(str1,str2)
 !STRCMP Compare to strings, a string and an array of strings or two equally
 !shaped arrays of strings with each other.
 !
 !Syntax:
 ! result=strcmp(str1,str2)
 !Input:
 ! str1,str2: CHARACTER arrays to be compared on a entry by entry basis
 !Output:
 ! results: LOGICAL array of the same shape as str1/str2 with T if strings
 !          are equal 
 
 CHARACTER(len=*),INTENT(in)::str1,str2
 LOGICAL::strcmp

 !Compare strings
 strcmp=TRIM(str1).EQ.TRIM(str2)
 
 END FUNCTION strcmp

!-----------------------------------------------------------------------------

 FUNCTION find(bool,dir) RESULT(index)
 !FIND Function with finds the first/last index in a 1D array which is true.
 !
 !Syntax:
 ! index=find(bool,dir)
 !Input:
 ! bool: 1D LOGICAL array
 ! dir:indicates whether search start from beginning ('first') or end ('last') 
 !Output:
 ! index: INTEGER given the index where bool is true for the first/last time

 CHARACTER(len=*),INTENT(in)::dir
 LOGICAL,INTENT(in)::bool(:)
 INTEGER::index,i1
  
 index=0
 IF (strcmp(dir,'first')) THEN
  DO i1=LBOUND(bool,1),UBOUND(bool,1)
   IF(bool(i1)) THEN
    index=i1; EXIT
   END IF
  END DO
 ELSEIF (strcmp(dir,'last')) THEN
  DO i1=UBOUND(bool,1),LBOUND(bool,1),-1
   IF(bool(i1)) THEN
    index=i1; EXIT
   END IF
  END DO
 ELSE
  WRITE(*,*) 'mod_netcdf.find: direction must be first or last'
  STOP
 END IF 
   
 END FUNCTION find


END MODULE
