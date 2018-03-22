MODULE mod_roms

USE mod_netcdf

TYPE roms_type
 INTEGER::grid_size(3),Vtransform,Vstretching,nTimes
 REAL(8)::theta_s,theta_b,Tcline      
 CHARACTER(len=1024)::filename
END TYPE roms_type

TYPE state
 REAL(8),allocatable,dimension(:,:,:,:)::temp,salt,u,v
 REAL(8),allocatable,dimension(:,:,:)::zeta
END TYPE state


INTERFACE collect_his
 MODULE PROCEDURE collect_his_2d, collect_his_3d
END INTERFACE collect_his

CONTAINS
!---------------------------------------------------------------------
SUBROUTINE roms_get_type(filename,fileDesign)
!Obtain basic layout from history file and write it to fileDesign

 IMPLICIT NONE

 !Declare
 CHARACTER(len=*),intent(in)::filename
 TYPE(roms_type),intent(out)::fileDesign
 INTEGER::status,ncid,dim_id(4),var_id(5)
 REAL(8),allocatable::sigma(:)

 !Filename
 fileDesign%filename=TRIM(filename)

 !Open netcdf stream
 status=nf_open(TRIM(filename),nf_nowrite,ncid)

 !Dimensions
 status=nf_inq_dimid(ncid,'xi_rho',dim_id(1))
 status=nf_inq_dimid(ncid,'eta_rho',dim_id(2))
 status=nf_inq_dimid(ncid,'s_rho',dim_id(3))
 IF(status.NE.nf_noerr) THEN
  status=nf_inq_dimid(ncid,'N',dim_id(3))
 END IF
 status=nf_inq_dimid(ncid,'ocean_time',dim_id(4)) 

 !Size
 status=nf_inq_dimlen(ncid,dim_id(1),fileDesign%grid_size(1))
 status=nf_inq_dimlen(ncid,dim_id(2),fileDesign%grid_size(2))
 status=nf_inq_dimlen(ncid,dim_id(3),fileDesign%grid_size(3))
 status=nf_inq_dimlen(ncid,dim_id(4),fileDesign%ntimes)

 !s-coordinate parameters
 status=nf_inq_varid(ncid,'Vtransform',var_id(1))
 status=nf_inq_varid(ncid,'Vstretching',var_id(2))
 status=nf_inq_varid(ncid,'theta_s',var_id(3))
 status=nf_inq_varid(ncid,'theta_b',var_id(4))
 status=nf_inq_varid(ncid,'Tcline',var_id(5))

 status=nf_get_var_int1(ncid,var_id(1),fileDesign%Vtransform)
 status=nf_get_var_int1(ncid,var_id(2),fileDesign%Vstretching)
 status=nf_get_var_double(ncid,var_id(3),fileDesign%theta_s)
 status=nf_get_var_double(ncid,var_id(4),fileDesign%theta_b)
 status=nf_get_var_double(ncid,var_id(5),fileDesign%Tcline)
 
 !Close stream
 status=nf_close(ncid)

END SUBROUTINE roms_get_type
 
!----------------------------------------------------------------------
SUBROUTINE roms_create_his_file(filename,fileDesign)
!Create empty history file 

 IMPLICIT NONE 

 !Declare
 CHARACTER(len=*),intent(in)::filename
 TYPE(roms_type),intent(in)::fileDesign
 INTEGER::status,ncid,var_id(20),dim_id(14),data_type
 INTEGER::i1,i2,i3,i0
 REAL(8),allocatable::sigma(:)

 data_type=nf_float
 status=nf_create(TRIM(filename),nf_classic_model,ncid)


 !Create dimensions
 status=nf_def_dim(ncid,'xi_rho',fileDesign%grid_size(1),dim_id(1))
 status=nf_def_dim(ncid,'xi_u',fileDesign%grid_size(1)-1,dim_id(2))
 status=nf_def_dim(ncid,'xi_v',fileDesign%grid_size(1),dim_id(3))
 status=nf_def_dim(ncid,'xi_psi',fileDesign%grid_size(1)-1,dim_id(4))
 status=nf_def_dim(ncid,'eta_rho',fileDesign%grid_size(2),dim_id(5))
 status=nf_def_dim(ncid,'eta_u',fileDesign%grid_size(2),dim_id(6))
 status=nf_def_dim(ncid,'eta_v',fileDesign%grid_size(2)-1,dim_id(7))
 status=nf_def_dim(ncid,'eta_psi',fileDesign%grid_size(2)-1,dim_id(8))
 status=nf_def_dim(ncid,'N',fileDesign%grid_size(3),dim_id(9))
 status=nf_def_dim(ncid,'s_rho',fileDesign%grid_size(3),dim_id(10))
 status=nf_def_dim(ncid,'s_w',fileDesign%grid_size(3)+1,dim_id(11))
 status=nf_def_dim(ncid,'ocean_time',nf_unlimited,dim_id(12))
 status=nf_def_dim(ncid,'tracer',2,dim_id(13))
 status=nf_def_dim(ncid,'boundary',4,dim_id(14))

 !Create variables
 status=nf_def_var(ncid,'Vtransform',nf_int,0,0,var_id(1))
 status=nf_def_var(ncid,'Vstretching',nf_int,0,0,var_id(2))
 status=nf_def_var(ncid,'theta_s',nf_double,0,0,var_id(3))
 status=nf_def_var(ncid,'theta_b',nf_double,0,0,var_id(4))
 status=nf_def_var(ncid,'Tcline',nf_double,0,0,var_id(5))
 status=nf_def_var(ncid,'s_rho',nf_double,1,dim_id(10),var_id(6))
 status=nf_def_var(ncid,'s_w',nf_double,1,dim_id(11),var_id(7))
 status=nf_def_var(ncid,'Cs_r',nf_double,1,dim_id(10),var_id(8))
 status=nf_def_var(ncid,'Cs_w',nf_double,1,dim_id(11),var_id(9))
 status=nf_def_var(ncid,'h',nf_double,2,[dim_id(1),dim_id(5)],var_id(10))
 status=nf_def_var(ncid,'ocean_time',nf_double,1,dim_id(12),var_id(11))
 status=nf_def_var(ncid,'zeta',data_type,3,&
 &[dim_id(1),dim_id(5),dim_id(12)],var_id(12))
 status=nf_def_var(ncid,'ubar',data_type,3,&
 &[dim_id(2),dim_id(6),dim_id(12)],var_id(13))
 status=nf_def_var(ncid,'vbar',data_type,3,&
 &[dim_id(3),dim_id(7),dim_id(12)],var_id(14))
 status=nf_def_var(ncid,'u',data_type,4,&
 &[dim_id(2),dim_id(6),dim_id(10),dim_id(12)],var_id(15))
 status=nf_def_var(ncid,'v',data_type,4,&
 &[dim_id(3),dim_id(7),dim_id(10),dim_id(12)],var_id(16))
 status=nf_def_var(ncid,'temp',data_type,4,&
 &[dim_id(1),dim_id(5),dim_id(10),dim_id(12)],var_id(17))
 status=nf_def_var(ncid,'salt',data_type,4,&
 &[dim_id(1),dim_id(5),dim_id(10),dim_id(12)],var_id(18))
 status=nf_def_var(ncid,'AKv',data_type,4,&
 &[dim_id(1),dim_id(5),dim_id(11),dim_id(12)],var_id(19))
 status=nf_def_var(ncid,'AKt',data_type,4,&
 &[dim_id(1),dim_id(5),dim_id(11),dim_id(12)],var_id(20))

 !Attributes 
 status=nf_put_att_text(ncid,var_id(12),'units',5,'meter')
 status=nf_put_att_text(ncid,var_id(12),'time',10,'ocean_time')
 !!status=nf_put_att_real(ncid,var_id(12),'_FillValue',nf_float,&
 !!&1,real(1e37))
 
 status=nf_put_att_text(ncid,var_id(13),'units',14,'meter second-1')
 status=nf_put_att_text(ncid,var_id(13),'time',10,'ocean_time')
 !!status=nf_put_att_real(ncid,var_id(13),'_FillValue',nf_float,&
 !!&1,real(1e37))

 status=nf_put_att_text(ncid,var_id(14),'units',14,'meter second-1')
 status=nf_put_att_text(ncid,var_id(14),'time',10,'ocean_time')
 !!status=nf_put_att_real(ncid,var_id(14),'_FillValue',nf_float,&
 !!&1,real(1e37))

 status=nf_put_att_text(ncid,var_id(15),'units',14,'meter second-1')
 status=nf_put_att_text(ncid,var_id(15),'time',10,'ocean_time')
 !!status=nf_put_att_real(ncid,var_id(15),'_FillValue',nf_float,&
 !!&1,real(1e37))

 status=nf_put_att_text(ncid,var_id(19),'units',15,'meter2 second-1')
 status=nf_put_att_text(ncid,var_id(19),'time',10,'ocean_time')
 !!status=nf_put_att_real(ncid,var_id(19),'_FillValue',nf_float,1,&
 !!&real(1e37))

 status=nf_put_att_text(ncid,var_id(20),'units',15,'meter2 second-1')
 status=nf_put_att_text(ncid,var_id(20),'time',10,'ocean_time')
 !!status=nf_put_att_real(ncid,var_id(20),'_FillValue',nf_float,&
 !!&1,real(1e37))

 status=nf_put_att_text(ncid,var_id(16),'units',14,'meter second-1')
 status=nf_put_att_text(ncid,var_id(16),'time',10,'ocean_time')
 !!status=nf_put_att_real(ncid,var_id(16),'_FillValue',nf_float,1,&
 !!&real(1e37))

 status=nf_put_att_text(ncid,var_id(17),'units',7,'Celsius')
 status=nf_put_att_text(ncid,var_id(17),'time',10,'ocean_time')
 !!status=nf_put_att_real(ncid,var_id(17),'_FillValue',nf_float,&
 !!&1,real(1e37))

 status=nf_put_att_text(ncid,var_id(12),'units',3,'ppt')
 status=nf_put_att_text(ncid,var_id(12),'time',10,'ocean_time')
 !!status=nf_put_att_real(ncid,var_id(12),'_FillValue',nf_float,&
 !!&1,real(1e37))

 status=nf_put_att_text(ncid,var_id(11),'units',7,'seconds')


 !Write s-grid parameters
 status=nf_enddef(ncid)

 status=nf_put_var_int(ncid,var_id(1),fileDesign%Vtransform)
 status=nf_put_var_int(ncid,var_id(2),fileDesign%Vstretching)
 status=nf_put_var_double(ncid,var_id(3),fileDesign%theta_s)
 status=nf_put_var_double(ncid,var_id(4),fileDesign%theta_b)
 status=nf_put_var_double(ncid,var_id(5),fileDesign%Tcline)
 
 ALLOCATE(sigma(fileDesign%grid_size(3)))
 DO i1=1,size(sigma)
  sigma(i1)=dble(-size(sigma)+i1-.5)/&
  &dble(fileDesign%grid_size(3))
 END DO
 status=nf_put_var_double(ncid,var_id(6),sigma)

 DEALLOCATE(sigma); ALLOCATE(sigma(fileDesign%grid_size(3)+1))
 DO i1=1,size(sigma)
  sigma(i1)=dble(-size(sigma)+1)/dble(fileDesign%grid_size(3))
 END DO
 status=nf_put_var_double(ncid,var_id(7),sigma)

 
 !Close stream
 status=nf_close(ncid)
 IF(status.NE.0) THEN
  WRITE(*,*) 'Error in creating file ',TRIM(filename)
  STOP
 END IF


END SUBROUTINE roms_create_his_file
 
!---------------------------------------------------------------------
SUBROUTINE collect_his_3d(filedir,filebase,var_name,times,val)
 !Collect 'var_name' field at times 'times' from different history files
 !in the directory 'filedir' and store them in 'val'. The subroutine only
 !searches in files of which the file name starts with 'filebase'
 
 IMPLICIT NONE

 !Declare
 CHARACTER(len=1024),intent(in)::filedir,filebase
 CHARACTER(len=64),intent(in)::var_name
 REAL(8),intent(in)::times(:)
 REAL(8),allocatable,intent(out)::val(:,:,:,:)

 CHARACTER(len=1024)::fileContents,command_line,filename
 INTEGER::n_files,n_base
 INTEGER::status,ncid

 REAL(8),allocatable::filetimes1(:),filetimes2(:)
 CHARACTER(len=1024),allocatable::filenames(:)

 INTEGER::itime(2),i0
 REAL(8),allocatable::times1(:),times2(:)
 REAL(8),allocatable,dimension(:,:,:)::val1,val2
 REAL(8)::w 

 !Read all files in directory
 WRITE(*,*) 'Reading files in directory'
 WRITE(fileContents,*) filedir,'/fileContents.txt'
 WRITE(command_line,*) 'ls --width=1 ',TRIM(filedir),&
 &' > ',TRIM(fileContents)
 CALL system(TRIM(command_line))
 
 !Count the number of files in the directory
 n_files=0
 n_base=LEN(TRIM(filebase))
 OPEN(101,file=TRIM(fileContents),action='read')
 DO WHILE(.TRUE.)
  READ(101,*,iostat=status)  filename
  IF(status.LT.0) EXIT
  IF(filename(1:n_base).NE.filebase(1:n_base)) CYCLE
  n_files=n_files+1
 END DO
 CLOSE(101)
 WRITE(*,*) n_files,' files in directory'

 !Read for each suitable file the filename, the 1st time and 
 !the last time
 ALLOCATE(filenames(n_files))
 ALLOCATE(filetimes1(n_files)); ALLOCATE(filetimes2(n_files))
 OPEN(101,file=TRIM(fileContents),action='read',dispose='delete')
 DO i0=1,n_files
  READ(101,*,iostat=status)  filename
  IF(status.LT.0) EXIT
  IF(filename(1:n_base).NE.filebase(1:n_base)) CYCLE
   
  WRITE(filenames(i0),*) TRIM(filedir),'/',TRIM(filename)
  status=nf_open(TRIM(filenames(i0)),nf_nowrite,ncid)
  IF(ALLOCATED(times1)) DEALLOCATE(times1)
  CALL read_nf_field(ncid,'ocean_time',0,times1)
  filetimes1(i0)=MINVAL(times1); filetimes2(i0)=MAXVAL(times1)
  status=nf_close(ncid)

  WRITE(*,*) 'file ',TRIM(filenames(i0)),' with times ',&
  filetimes1(i0),'-',filetimes2(i0)
 END DO
 CLOSE(101)

 !For each time find the model output at time before and after output
 !time
 DO i0=1,size(times,1)
  itime=[0,0]
  itime(1)=MAXLOC(filetimes1,1,filetimes1.LE.times(i0))
  itime(2)=MINLOC(filetimes2,1,filetimes2.GE.times(i0))
 
  IF(ANY(itime.EQ.0) ) CYCLE 

  !Retrieve field at time before output time
  status=nf_open(TRIM(filenames(itime(1))),nf_nowrite,ncid)
  IF(ALLOCATED(times1)) DEALLOCATE(times1)
  CALL read_nf_field(ncid,'ocean_time',0,times1)
  itime(1)=MAXLOC(times1,1,times1.LE.times(i0))
  CALL read_nf_field(ncid,TRIM(var_name),&
  &itime(1),val1)
  status=nf_close(ncid)

  !Retrieve field at time after output time
  status=nf_open(TRIM(filenames(itime(2))),nf_nowrite,ncid)
  IF(ALLOCATED(times2)) DEALLOCATE(times2)
  CALL read_nf_field(ncid,'ocean_time',0,times2)
  itime(2)=MINLOC(times2,1,times2.GE.times(i0))
  CALL read_nf_field(ncid,TRIM(var_name),&
  &itime(2),val2)
  status=nf_close(ncid)
    
  WRITE(*,*) 'Interpolating for time ',times(i0),' using times ',&
  &times1(itime(1)),' and ',times2(itime(2))

  !Perform linear interpolation
  IF(times2(itime(2)).EQ.times1(itime(1))) THEN 
   w=0
  ELSE
   w=(times(i0)-times1(itime(1)))/(times2(itime(2))-times1(itime(1)))
  END IF
  IF(.NOT.ALLOCATED(val)) ALLOCATE(val(size(val1,1),size(val1,2),&
  &size(val1,3),size(times,1)))
  val(:,:,:,i0)=(1-w)*val1+w*val2
 END DO

END SUBROUTINE collect_his_3d

!-------------------------------------------------------------------------
SUBROUTINE collect_his_2d(filedir,filebase,var_name,times,val)
 !Collect 'var_name' field at times 'times' from different history files
 !in the directory 'filedir' and store them in 'val'. The subroutine only
 !searches in files of which the file name starts with 'filebase'
 
 IMPLICIT NONE

 !Declare
 CHARACTER(len=1024),intent(in)::filedir,filebase
 CHARACTER(len=64),intent(in)::var_name
 REAL(8),intent(in)::times(:)
 REAL(8),allocatable,intent(out)::val(:,:,:)

 CHARACTER(len=1024)::fileContents,command_line,filename
 INTEGER::n_files,n_base
 INTEGER::status,ncid

 REAL(8),allocatable::filetimes1(:),filetimes2(:)
 CHARACTER(len=1024),allocatable::filenames(:)

 INTEGER::itime(2),i0
 REAL(8),allocatable::times1(:),times2(:)
 REAL(8),allocatable,dimension(:,:)::val1,val2
 REAL(8)::w 

 !Read all files in directory
 WRITE(*,*) 'Reading files in directory'
 WRITE(fileContents,*) filedir,'/fileContents.txt'
 WRITE(command_line,*) 'ls --width=1 ',TRIM(filedir),&
 &' > ',TRIM(fileContents)
 CALL system(TRIM(command_line))
 
 !Count the number of files in the directory
 n_files=0
 n_base=LEN(TRIM(filebase))
 OPEN(101,file=TRIM(fileContents),action='read')
 DO WHILE(.TRUE.)
  READ(101,*,iostat=status)  filename
  IF(status.LT.0) EXIT
  IF(filename(1:n_base).NE.filebase(1:n_base)) CYCLE
  n_files=n_files+1
 END DO
 CLOSE(101)
 WRITE(*,*) n_files,' files in directory'

 !Read for each suitable file the filename, the 1st time and 
 !the last time
 ALLOCATE(filenames(n_files))
 ALLOCATE(filetimes1(n_files)); ALLOCATE(filetimes2(n_files))
 OPEN(101,file=TRIM(fileContents),action='read',dispose='delete')
 DO i0=1,n_files
  READ(101,*,iostat=status)  filename
  IF(status.LT.0) EXIT
  IF(filename(1:n_base).NE.filebase(1:n_base)) CYCLE
   
  WRITE(filenames(i0),*) TRIM(filedir),'/',TRIM(filename)
  status=nf_open(TRIM(filenames(i0)),nf_nowrite,ncid)
  IF(ALLOCATED(times1)) DEALLOCATE(times1)
  CALL read_nf_field(ncid,'ocean_time',0,times1)
  filetimes1(i0)=MINVAL(times1); filetimes2(i0)=MAXVAL(times1)
  status=nf_close(ncid)

  WRITE(*,*) 'file ',TRIM(filenames(i0)),' with times ',&
  filetimes1(i0),'-',filetimes2(i0)
 END DO
 CLOSE(101)

 !For each time find the model output at time before and after output
 !time
 DO i0=1,size(times,1)
  itime=[0,0]
  itime(1)=MAXLOC(filetimes1,1,filetimes1.LE.times(i0))
  itime(2)=MINLOC(filetimes2,1,filetimes2.GE.times(i0))
 
  IF(ANY(itime.EQ.0)) CYCLE 

  !Retrieve field at time before output time
  status=nf_open(TRIM(filenames(itime(1))),nf_nowrite,ncid)
  IF(ALLOCATED(times1)) DEALLOCATE(times1)
  CALL read_nf_field(ncid,'ocean_time',0,times1)
  itime(1)=MAXLOC(times1,1,times1.LE.times(i0))
  CALL read_nf_field(ncid,TRIM(var_name),&
  &itime(1),val1)
  status=nf_close(ncid)

  !Retrieve field at time after output time
  status=nf_open(TRIM(filenames(itime(2))),nf_nowrite,ncid)
  IF(ALLOCATED(times2)) DEALLOCATE(times2)
  CALL read_nf_field(ncid,'ocean_time',0,times2)
  itime(2)=MINLOC(times2,1,times2.GE.times(i0))
  CALL read_nf_field(ncid,TRIM(var_name),&
  &itime(2),val2)
  status=nf_close(ncid)
    
  WRITE(*,*) 'Interpolating for time ',times(i0),' using times ',&
  &times1(itime(1)),' and ',times2(itime(2))

  !Perform linear interpolation
  IF(times2(itime(2)).EQ.times1(itime(1))) THEN 
   w=0
  ELSE
   w=(times(i0)-times1(itime(1)))/(times2(itime(2))-times1(itime(1)))
  END IF
  IF(.NOT.ALLOCATED(val)) ALLOCATE(val(size(val1,1),size(val1,2),&
  &size(times,1)))
  val(:,:,i0)=(1-w)*val1+w*val2
 END DO

END SUBROUTINE collect_his_2d

!-------------------------------------------------------------------------
END MODULE mod_roms
