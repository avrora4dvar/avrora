PROGRAM create_background
 USE mod_interp
 USE mod_netcdf
 USE mod_roms
 USE mod_time
 IMPLICIT NONE

 TYPE grd_type
  INTEGER::s(8)
  LOGICAL,allocatable,dimension(:,:)::maskr,masku,maskv
  REAL(8),allocatable,dimension(:)::lonr,lonu,lonv
  REAL(8),allocatable,dimension(:)::latr,latu,latv
  REAL(8),allocatable,dimension(:,:)::h
  REAL(8),allocatable,dimension(:,:,:)::zw,zr,zdr,zdu,zdv
 END TYPE grd_type
 
 INTEGER,allocatable,dimension(:,:)::mask1

 !Input
 CHARACTER(len=1024)::grid_file(2),his_dir,outfile
 REAL(8)::T0,S0,Akt0,Akv0,tStart,tEnd,refShift,tStep
 INTEGER::iHisBnd(2),tiles(2),tref(6)
 CHARACTER::flag_filter_in
 LOGICAL::flag_filter

 !Grid
 TYPE(grd_type)::grd(2)
 INTEGER::ilat(3)
 REAL(8),allocatable::dist(:,:),sigma(:)

 !Time reading
 INTEGER::ntimes,iHis,itStart,itEnd
 INTEGER,allocatable::ntFile(:)
 CHARACTER(len=1024)::filename
 LOGICAL::flag_exist
 REAL(8),allocatable::tOut(:)
 REAL(8)::w
 INTEGER::ibnd(2)

 !Netcdf/counters
 INTEGER::i0,i1,i2,i3,i4,j0,j1,j2,j3,j4,s(8)
 INTEGER::status,gid(2),oid
 INTEGER,allocatable::nid(:)
 INTEGER::nccount(4),ncstart(4)

 !Roms file
 TYPE(roms_type)::fdesign
 REAL(8),allocatable::time(:)

 !Data copying
 REAL(8),allocatable,dimension(:,:,:)::val3d1,val3d2,zeta
 REAL(8),allocatable,dimension(:,:,:,:)::val4d1,val4d2

 !MPI
 INTEGER::np,myrank
 INTEGER,allocatable::itp(:,:)
#ifdef MPI
 INCLUDE 'mpif.h'
 INTEGER::iostatus(mpi_status_size)
 INTEGER::mpi_double=mpi_double_precision
#endif

 !Default MPI settings
 myrank=0; np=1
!---------------------------------------------------------------------------
!Read input

IF(myrank.EQ.0) THEN

 READ(*,*) !Grid file
 READ(*,'(A)') grid_file(1)
 READ(*,'(A)') grid_file(2)
 READ(*,*) !Directory with history files
 READ(*,'(A)') his_dir
 READ(*,*) !Number of the first and last ocean_his_####.nc to be read
 READ(*,*) iHisBnd(1),iHisBnd(2)
 READ(*,*) !Output file
 READ(*,'(A)') outfile
 READ(*,*) !First,last time (days since date_ref) and time step (s)
 READ(*,*) tStart, tEnd, tStep
 READ(*,*) !Reference date (date_ref)
 READ(*,*) tref(1),tref(2),tref(3),tref(4),tref(5),tref(6)
 READ(*,*) !Shift reference date output (days)
 READ(*,*) refShift
 READ(*,*) !Background diffusivity
 READ(*,*) T0,S0,AKt0,AKv0
 READ(*,*) !Use of filtering (T/F)
 READ(*,'(A)') flag_filter_in
 READ(*,*) !'Tiles'
 READ(*,*) tiles(1),tiles(2)

 IF(flag_filter_in.EQ.'T') THEN
  flag_filter=.true.
 ELSE
  flag_filter=.false.
 END IF

 !Times in output
 ALLOCATE(tOut( CEILING((tEnd-tStart)*86400/tStep)+1 ))
 DO i1=1,size(tOut)
  tOut(i1)=tStart*86400+dble(i1-1)*tStep
 END DO

 WRITE(*,*) 'Grid file'
 WRITE(*,*) TRIM(grid_file(1))
 WRITE(*,*) TRIM(grid_file(2))
 WRITE(*,*) 'Starting create_background'
 WRITE(*,*) '#Directory with history files'
 WRITE(*,*) TRIM(his_dir)
 WRITE(*,*) '#Number of the first and last ocean_his_####.nc to be read'
 WRITE(*,*) iHisBnd
 WRITE(*,*) '#Output file'
 WRITE(*,*) TRIM(outfile)
 WRITE(*,*) '#First time to be written to output file'
 WRITE(*,*) unix2cal(cal2unix(tRef)+INT(tOut(1)))
 WRITE(*,*) '#Last time to be written to output file'
 WRITE(*,*) unix2cal(cal2unix(tRef)+INT(tOut(size(tOut))))
 WRITE(*,*) '#Reference time input'
 WRITE(*,*) tRef
 WRITE(*,*) '#Reference time output'
 WRITE(*,*) unix2cal(cal2unix(tRef)+INT(refShift*86400))
 WRITE(*,*) '#Use of filtering (T/F)'
 WRITE(*,*) flag_filter
 WRITE(*,*) '#Tiles'
 WRITE(*,*) tiles

END IF !myrank

!---------------------------------------------------------------------------
!Read grid files

IF(myrank.EQ.0) THEN

 DO i0=1,2
  !Open stream
  status=nf_open(TRIM(grid_file(i0)),nf_nowrite,gid(i0))
   
  !Allocate
  CALL ncsize(gid(i0),'z0_r',grd(i0)%s)
  ALLOCATE(grd(i0)%lonr(grd(i0)%s(1)))
  ALLOCATE(grd(i0)%lonu(grd(i0)%s(1)-1))
  ALLOCATE(grd(i0)%lonv(grd(i0)%s(1)))
  ALLOCATE(grd(i0)%latr(grd(i0)%s(2)))
  ALLOCATE(grd(i0)%latu(grd(i0)%s(2)))
  ALLOCATE(grd(i0)%latv(grd(i0)%s(2)-1))
  ALLOCATE(grd(i0)%maskr(grd(i0)%s(1),grd(i0)%s(2)))
  ALLOCATE(grd(i0)%masku(grd(i0)%s(1)-1,grd(i0)%s(2)))
  ALLOCATE(grd(i0)%maskv(grd(i0)%s(1),grd(i0)%s(2)-1))

  !Read mask
  grd(i0)%maskr=ncread2d(gid(i0),'mask_rho',[1,1],&
  &[grd(i0)%s(1),grd(i0)%s(2)]).NE.0
  grd(i0)%masku=ncread2d(gid(i0),'mask_u',[1,1],&
  &[grd(i0)%s(1)-1,grd(i0)%s(2)]).NE.0
  grd(i0)%maskv=ncread2d(gid(i0),'mask_v',[1,1],&
  &[grd(i0)%s(1),grd(i0)%s(2)-1]).NE.0

  !Read lat
  grd(i0)%latr=RESHAPE(ncread2d(gid(i0),'lat_rho',[1,1],&
  &[1,grd(i0)%s(2)]),[grd(i0)%s(2)])
  grd(i0)%latu=RESHAPE(ncread2d(gid(i0),'lat_u',[1,1],&
  &[1,grd(i0)%s(2)]),[grd(i0)%s(2)])
  grd(i0)%latv=RESHAPE(ncread2d(gid(i0),'lat_v',[1,1],&
  &[1,grd(i0)%s(2)-1]),[grd(i0)%s(2)-1])
 
 END DO	

 !Find latitude for which the difference with lon between the grids is minimal
 ALLOCATE(dist(size(grd(2)%latr),3)); dist=dble(0)
 DO i2=1,size(grd(2)%latr)
 dist(i2,1)=MINVAL(abs(grd(1)%latr-grd(2)%latr(i2)),1)
 dist(i2,2)=MINLOC(abs(grd(1)%latr-grd(2)%latr(i2)),1)
 dist(i2,3)=dist(i2,1)**2
 END DO
 ilat(2)=MINLOC(dist(:,3),1)
 ilat(1)=INT(dist(ilat(2),2));
 DEALLOCATE(dist)

 DO i0=1,2
  !Read lon
  grd(i0)%lonr=RESHAPE(ncread2d(gid(i0),'lon_rho',[1,ilat(i0)],&
  &[grd(i0)%s(1),1]),[grd(i0)%s(1)])
  grd(i0)%lonu=RESHAPE(ncread2d(gid(i0),'lon_u',[1,ilat(i0)],&
  &[grd(i0)%s(1)-1,1]),[grd(i0)%s(1)-1])
  grd(i0)%lonv=RESHAPE(ncread2d(gid(i0),'lon_v',[1,ilat(i0)],&
  &[grd(i0)%s(1),1]),[grd(i0)%s(1)])
 END DO

 !Close streams
 DO i0=1,2
  status=nf_close(gid(i0))
 END DO

END IF !myrank

!---------------------------------------------------------------------------
!Count the number of times

IF(myrank.EQ.0) THEN
 WRITE(*,*) 'Collecting time in history files'

 ntimes=0
 ALLOCATE(nid(iHisBnd(1):iHisBnd(2))); nid=0
 ALLOCATE(ntfile(iHisBnd(1):iHisBnd(2))); ntFile=0
 DO ihis=iHisBnd(1),iHisBnd(2)
  !Open file
  WRITE(filename,'(A,A,I4.4,A)') TRIM(his_dir),'/ocean_his_',&
  &iHis,'.nc'
  INQUIRE(file=TRIM(filename),exist=flag_exist)
  IF(.NOT.flag_exist) THEN
   WRITE(*,*) TRIM(filename),' not present.&
   & Ending search.'
   EXIT
  ELSE
   WRITE(*,*) 'Opening ',TRIM(filename)
  END IF

  !Get number of times in input
  status=nf_open(TRIM(filename),nf_nowrite,nid(iHis))
  CALL ncsize(nid(iHis),'ocean_time',s)

  !Update file counter
  ntFile(iHis)=s(1)

  !Get roms layout
  IF(iHis.EQ.iHisBnd(1)) THEN
   CALL roms_get_type(TRIM(filename),fdesign)

  END IF
 END DO
 ntimes=SUM(ntFile)

 IF(ntimes.EQ.0) THEN
  WRITE(*,*) 'No entries within time frame. Ending program.'
  STOP
 END IF

 !Distribution of time over different processes
 ALLOCATE(itp(0:np-1,3)); itp=1
 DO i0=1,np
  IF(i0.GT.1) itp(i0-1,1)=itp(i0-2,2)+1
  itp(i0-1,2)=NINT(dble(i0*ntimes)/dble(np))
  itp(i0-1,3)=itp(i0-1,2)-itp(i0-1,1)+1
 END DO

END IF !myrank

!--------------------------------------------------------------------------
!Output file

IF(myrank.EQ.0) THEN

 WRITE(*,*) 'Creating output file'
 ALLOCATE(time(ntimes)); time=dble(0)
 i0=1; 
 DO iHis=iHisBnd(1),iHisBnd(2)
  time(i0:i0+ntFile(iHis)-1)=ncread1d(nid(iHis),'ocean_time',[1],&
  &[ntFile(iHis)])
  i0=i0+ntFile(iHis)
 END DO
 
 !Count times in output file
 fdesign%ntimes=size(tOut)
 IF(MINVAL(time).GT.MINVAL(tOut).OR.&
 &MAXVAL(time).LT.MAXVAL(tOut)) THEN
   WRITE(*,*) 'Output times outside range input files'
   STOP
 END IF

 !Create output file
 fdesign%grid_size(1)=SIZE(grd(2)%lonr)
 fdesign%grid_size(2)=SIZE(grd(2)%latr)
 CALL roms_create_his_file(TRIM(outfile),fdesign)

 !Open stream to output
 status=nf_open(TRIM(outfile),nf_write,oid)

END IF !myrank

!-------------------------------------------------------------------------
!MPI setup

 !Not supported currently

!--------------------------------------------------------------------------
!Read zeta

 WRITE(*,*) 'Interpolating zeta'
 ALLOCATE(val3d1(grd(1)%s(1),grd(1)%s(2),ntimes))
 ALLOCATE(val3d2(grd(2)%s(1),grd(2)%s(2),ntimes))

 !Read
 i4=1
 DO iHis=iHisBnd(1),iHisBnd(2)
  ncstart=[1,1,1,0]
  nccount=[grd(1)%s(1),grd(1)%s(2),ntFile(iHis),0]
  val3d1(:,:,i4:i4-1+ntFile(iHis))=&
  &ncread3d(nid(iHis),'zeta',ncstart(1:3),nccount(1:3))
  i4=i4+ntFile(iHis)
 END DO

 !Interpolate
 val3d2(:,:,itp(myrank,1):itp(myrank,2))=&
 &mesh2mesh_interp(grd(1)%lonr,grd(1)%latr,grd(1)%maskr,&
 &val3d1(:,:,itp(myrank,1):itp(myrank,2)),grd(2)%lonr,grd(2)%latr)

 !Extrapolate
 CALL mesh2mesh_extrap(grd(2)%lonr,grd(2)%latr,&
 &.NOT.ANY(val3d2.EQ.interp_fill_value,3),val3d2,grd(2)%maskr,dble(0))

 !Filter
 IF(flag_filter) THEN
 END IF

 !Write to output
 DO i1=1,size(tOut)
  ibnd(1)=MAXLOC(time,1,time.LE.tOut(i1))
  ibnd(2)=MINLOC(time,1,time.GE.tOut(i1))
  IF(ibnd(1).EQ.ibnd(2)) THEN
   w=.5
  ELSE
   w=(tOut(i1)-time(ibnd(1)))/(time(ibnd(2))-time(ibnd(1)))
  END IF
  ncstart=[1,1,i1,0]
  nccount=[size(val3d2,1),size(val3d2,2),1,0]
  CALL ncwrite3d(oid,'zeta',ncstart(1:3),nccount(1:3),&
  &(1.0-w)*val3d2(:,:,ibnd(1):ibnd(1))&
  &+w*val3d2(:,:,ibnd(2):ibnd(2)))
 END DO

 !Deallocate
 DEALLOCATE(val3d1,val3d2)

!----------------------------------------------------------------------------
!Read ubar
  
 WRITE(*,*) 'Interpolating ubar'
 ALLOCATE(val3d1(grd(1)%s(1)-1,grd(1)%s(2),ntimes))
 ALLOCATE(val3d2(grd(2)%s(1)-1,grd(2)%s(2),ntimes))
 
 i4=1
 DO iHis=iHisBnd(1),iHisBnd(2)
  ncstart=[1,1,1,0]
  nccount=[grd(1)%s(1)-1,grd(1)%s(2),ntFile(iHis),0]
  val3d1(:,:,i4:i4-1+ntFile(iHis))=&
  &ncread3d(nid(iHis),'ubar',ncstart(1:3),nccount(1:3))
  i4=i4+ntFile(iHis)
 END DO

 !Interpolate
 val3d2=&
 &mesh2mesh_interp(grd(1)%lonu,grd(1)%latu,grd(1)%masku,&
 &val3d1,grd(2)%lonu,grd(2)%latu)

 !Extrapolate
 CALL mesh2mesh_extrap(grd(2)%lonu,grd(2)%latu,&
 &.NOT.ANY(val3d2.EQ.interp_fill_value,3),val3d2,grd(2)%masku,dble(0))
 
 !Filter
 IF(flag_filter) THEN
 END IF

 !Write to output
 DO i1=1,size(tOut)
  ibnd(1)=MAXLOC(time,1,time.LE.tOut(i1))
  ibnd(2)=MINLOC(time,1,time.GE.tOut(i1))
  IF(ibnd(1).EQ.ibnd(2)) THEN
   w=.5
  ELSE
   w=(tOut(i1)-time(ibnd(1)))/(time(ibnd(2))-time(ibnd(1)))
  END IF 
  ncstart=[1,1,i1,0]
  nccount=[size(val3d2,1),size(val3d2,2),1,0]
  CALL ncwrite3d(oid,'ubar',ncstart(1:3),nccount(1:),&
  &(1.0-w)*val3d2(:,:,ibnd(1):ibnd(1))+&
  &w*val3d2(:,:,ibnd(2):ibnd(2)))
 END DO

 !Deallocate
 DEALLOCATE(val3d1,val3d2)

!---------------------------------------------------------------------------
!Read vbar
  
 WRITE(*,*) 'Interpolating vbar'
 ALLOCATE(val3d1(grd(1)%s(1),grd(1)%s(2)-1,ntimes))
 ALLOCATE(val3d2(grd(2)%s(1),grd(2)%s(2)-1,ntimes))
 
 i4=1
 DO iHis=iHisBnd(1),iHisBnd(2)
  ncstart=[1,1,1,0]
  nccount=[grd(1)%s(1),grd(1)%s(2)-1,ntFile(iHis),0]
  val3d1(:,:,i4:i4-1+ntFile(iHis))=&
  &ncread3d(nid(iHis),'vbar',ncstart(1:3),nccount(1:3))
  i4=i4+ntFile(iHis)
 END DO

 !Interpolate
 val3d2(:,:,:)=&
 &mesh2mesh_interp(grd(1)%lonv,grd(1)%latv,grd(1)%maskv,&
 &val3d1(:,:,:),grd(2)%lonv,grd(2)%latv)

 !Extrapolate
 CALL mesh2mesh_extrap(grd(2)%lonv,grd(2)%latv,&
 &.NOT.ANY(val3d2.EQ.interp_fill_value,3),val3d2,grd(2)%maskv,dble(0))

 !Filter
 IF(flag_filter) THEN
 END IF

 !Write to output
 DO i1=1,size(tOut)
  ibnd(1)=MAXLOC(time,1,time.LE.tOut(i1))
  ibnd(2)=MINLOC(time,1,time.GE.tOut(i1))
  IF(ibnd(1).EQ.ibnd(2)) THEN
   w=.5
  ELSE
   w=(tOut(i1)-time(ibnd(1)))/(time(ibnd(2))-time(ibnd(1)))
  END IF
 
  ncstart=[1,1,i1,0]
  nccount=[size(val3d2,1),size(val3d2,2),1,0]
  CALL ncwrite3d(oid,'vbar',ncstart(1:3),nccount(1:),&
  &(1.0-w)*val3d2(:,:,iBnd(1):iBnd(1))+&
  &w*val3d2(:,:,iBnd(2):iBnd(2)))
 END DO
 
 !Deallocate
 DEALLOCATE(val3d1,val3d2)

!-----------------------------------------------------------------------------
!Read temp
  
 WRITE(*,*) 'Interpolating temp'
 ALLOCATE(val4d1(grd(1)%s(1),grd(1)%s(2),grd(1)%s(3),ntimes))
 ALLOCATE(val4d2(grd(2)%s(1),grd(2)%s(2),grd(2)%s(3),ntimes))
 
 i4=1
 DO iHis=iHisBnd(1),iHisBnd(2)
  ncstart=[1,1,1,1]
  nccount=[grd(1)%s(1),grd(1)%s(2),grd(1)%s(3),ntFile(iHis)]
  val4d1(:,:,:,i4:i4-1+ntFile(iHis))=&
  &ncread4d(nid(iHis),'temp',ncstart(1:4),nccount(1:4))
  i4=i4+ntFile(iHis)
 END DO

 !Interpolate
 DO i0=0,np-1
  val4d2(:,:,:,itp(i0,1):itp(i0,2))=&
  &mesh2mesh_interp(grd(1)%lonr,grd(1)%latr,grd(1)%maskr,&
  &val4d1(:,:,:,itp(i0,1):itp(i0,2)),grd(2)%lonr,grd(2)%latr)
 END DO

 !Extrapolate
 CALL mesh2mesh_extrap(grd(2)%lonr,grd(2)%latr,&
 &.NOT.ANY(ANY(val4d2.EQ.interp_fill_value,4),3),&
 &val4d2,grd(2)%maskr,T0)

 !Filter
 IF(flag_filter) THEN
 END IF

 !Write to output
 DO i1=1,size(tOut)
  ibnd(1)=MAXLOC(time,1,time.LE.tOut(i1))
  ibnd(2)=MINLOC(time,1,time.GE.tOut(i1))
  IF(ibnd(1).EQ.ibnd(2)) THEN
   w=.5
  ELSE
   w=(tOut(i1)-time(ibnd(1)))/(time(ibnd(2))-time(ibnd(1)))
  END IF
 
  ncstart=[1,1,1,i1]
  nccount=[size(val4d2,1),size(val4d2,2),size(val4d2,3),&
  &1]
  CALL ncwrite4d(oid,'temp',ncstart(1:4),nccount(1:4),&
  &(1.0-w)*val4d2(:,:,:,iBnd(1):iBnd(1))+&
  &w*val4d2(:,:,:,iBnd(2):iBnd(2)))
 END DO


!-----------------------------------------------------------------------------
!Read salt
 
 WRITE(*,*) 'Interpolating salt'
 val4d1=dble(0); val4d2=dble(0)

 i4=1
 DO iHis=iHisBnd(1),iHisBnd(2)
  ncstart=[1,1,1,1]
  nccount=[grd(1)%s(1),grd(1)%s(2),grd(1)%s(3),ntFile(iHis)]
  val4d1(:,:,:,i4:i4-1+ntFile(iHis))=&
  &ncread4d(nid(iHis),'salt',ncstart(1:4),nccount(1:4))
  i4=i4+ntFile(iHis)
 END DO

 DO i0=0,np-1
  !Interpolate
  val4d2(:,:,:,itp(i0,1):itp(i0,2))=&
  &mesh2mesh_interp(grd(1)%lonr,grd(1)%latr,grd(1)%maskr,&
  &val4d1(:,:,:,itp(i0,1):itp(i0,2)),grd(2)%lonr,grd(2)%latr)
 END DO

 !Extrapolate
 CALL mesh2mesh_extrap(grd(2)%lonr,grd(2)%latr,&
 &.NOT.ANY(ANY(val4d2.EQ.interp_fill_value,4),3),&
 &val4d2,grd(2)%maskr,S0)

 !Filter
 IF(flag_filter) THEN
 END IF

 !Write to output
 DO i1=1,size(tOut)
  ibnd(1)=MAXLOC(time,1,time.LE.tOut(i1))
  ibnd(2)=MINLOC(time,1,time.GE.tOut(i1))
  IF(ibnd(1).EQ.ibnd(2)) THEN
   w=.5
  ELSE
   w=(tOut(i1)-time(ibnd(1)))/(time(ibnd(2))-time(ibnd(1)))
  END IF
 
  ncstart=[1,1,1,i1]
  nccount=[size(val4d2,1),size(val4d2,2),size(val4d2,3),&
  &1]
  CALL ncwrite4d(oid,'salt',ncstart(1:4),nccount(1:4),&
  &(1.0-w)*val4d2(:,:,:,iBnd(1):iBnd(1))&
  &+w*val4d2(:,:,:,iBnd(2):iBnd(2)))
 END DO

 DEALLOCATE(val4d1,val4d2) 

!-----------------------------------------------------------------------------
!Read AKt
  
 WRITE(*,*) 'Interpolating AKt'
 ALLOCATE(val4d1(grd(1)%s(1),grd(1)%s(2),grd(1)%s(3)+1,ntimes))
 ALLOCATE(val4d2(grd(2)%s(1),grd(2)%s(2),grd(2)%s(3)+1,ntimes))
 
 i4=1
 DO iHis=iHisBnd(1),iHisBnd(2)
  ncstart=[1,1,1,1]
  nccount=[grd(1)%s(1),grd(1)%s(2),grd(1)%s(3)+1,ntFile(iHis)]
  val4d1(:,:,:,i4:i4-1+ntFile(iHis))=&
  &ncread4d(nid(iHis),'AKt',ncstart(1:4),nccount(1:4))
  i4=i4+ntFile(iHis)
 END DO

 DO i0=0,np-1
  !Interpolate
  val4d2(:,:,:,itp(i0,1):itp(i0,2))=&
  &mesh2mesh_interp(grd(1)%lonr,grd(1)%latr,grd(1)%maskr,&
  &val4d1(:,:,:,itp(i0,1):itp(i0,2)),grd(2)%lonr,grd(2)%latr)
 END DO

 !Extrapolate
 CALL mesh2mesh_extrap(grd(2)%lonr,grd(2)%latr,&
 &.NOT.ANY(ANY(val4d2.EQ.interp_fill_value,4),3),&
 &val4d2,grd(2)%maskr,AKt0)
 val4d2=MAX(val4d2,Akt0)

 !Filter
 IF(flag_filter) THEN
 END IF

 !Write to output
 DO i1=1,size(tOut)
  ibnd(1)=MAXLOC(time,1,time.LE.tOut(i1))
  ibnd(2)=MINLOC(time,1,time.GE.tOut(i1))
  IF(ibnd(1).EQ.ibnd(2)) THEN
   w=.5
  ELSE
   w=(tOut(i1)-time(ibnd(1)))/(time(ibnd(2))-time(ibnd(1)))
  END IF
 
  ncstart=[1,1,1,i1]
  nccount=[size(val4d2,1),size(val4d2,2),size(val4d2,3),&
  &1]
  CALL ncwrite4d(oid,'AKt',ncstart(1:4),nccount(1:4),&
  &(1.0-w)*val4d2(:,:,:,iBnd(1):iBnd(1))+&
  &w*val4d2(:,:,:,iBnd(2):iBnd(2)))
 END DO
!-----------------------------------------------------------------------------
!Read AKv
  
 WRITE(*,*) 'Interpolating Akv'
 val4d1=dble(0); val4d2=dble(0)
 
 i4=1
 DO iHis=iHisBnd(1),iHisBnd(2)
  ncstart=[1,1,1,1]
  nccount=[grd(1)%s(1),grd(1)%s(2),grd(1)%s(3)+1,ntFile(iHis)]
  val4d1(:,:,:,i4:i4-1+ntFile(iHis))=&
  &ncread4d(nid(iHis),'AKv',ncstart(1:4),nccount(1:4))
  i4=i4+ntFile(iHis)
 END DO

 DO i0=0,np-1
  !Interpolate
  val4d2(:,:,:,itp(i0,1):itp(i0,2))=&
  &mesh2mesh_interp(grd(1)%lonr,grd(1)%latr,grd(1)%maskr,&
  &val4d1(:,:,:,itp(i0,1):itp(i0,2)),grd(2)%lonr,grd(2)%latr)
END DO

 !Extrapolate
 CALL mesh2mesh_extrap(grd(2)%lonr,grd(2)%latr,&
 &.NOT.ANY(ANY(val4d2.EQ.interp_fill_value,4),3),&
 &val4d2,grd(2)%maskr,Akv0)
 val4d2=MAX(val4d2,Akv0)
  
 !Filter
 IF(flag_filter) THEN
 END IF

 !Write to output
 DO i1=1,size(tOut)
  ibnd(1)=MAXLOC(time,1,time.LE.tOut(i1))
  ibnd(2)=MINLOC(time,1,time.GE.tOut(i1))
  IF(ibnd(1).EQ.ibnd(2)) THEN
   w=.5
  ELSE
   w=(tOut(i1)-time(ibnd(1)))/(time(ibnd(2))-time(ibnd(1)))
  END IF

  ncstart=[1,1,1,i1]
  nccount=[size(val4d2,1),size(val4d2,2),size(val4d2,3),&
  &1]
  CALL ncwrite4d(oid,'AKv',ncstart(1:4),nccount(1:4),&
  &(1.0-w)*val4d2(:,:,:,iBnd(1):iBnd(1))+&
  &w*val4d2(:,:,:,iBnd(2):iBnd(2)))
 END DO

 !deallocate
 DEALLOCATE(val4d1,val4d2)

!-----------------------------------------------------------------------------
!Read u
 
 WRITE(*,*) 'Interpolating u'
 ALLOCATE(val4d1(grd(1)%s(1)-1,grd(1)%s(2),grd(1)%s(3),ntimes))
 ALLOCATE(val4d2(grd(2)%s(1)-1,grd(2)%s(2),grd(2)%s(3),ntimes))

 i4=1
 DO iHis=iHisBnd(1),iHisBnd(2)
  ncstart=[1,1,1,1]
  nccount=[grd(1)%s(1)-1,grd(1)%s(2),grd(1)%s(3),ntFile(iHis)]
  val4d1(:,:,:,i4:i4-1+ntFile(iHis))=&
  &ncread4d(nid(iHis),'u',ncstart(1:4),nccount(1:4))
  i4=i4+ntFile(iHis)
 END DO

 DO i0=0,np-1
  !Interpolate
  val4d2(:,:,:,itp(i0,1):itp(i0,2))=&
  &mesh2mesh_interp(grd(1)%lonu,grd(1)%latu,grd(1)%masku,&
  &val4d1(:,:,:,itp(i0,1):itp(i0,2)),grd(2)%lonu,grd(2)%latu)
 END DO

 !Extrapolate
 CALL mesh2mesh_extrap(grd(2)%lonu,grd(2)%latu,&
 &.NOT.ANY(ANY(val4d2.EQ.interp_fill_value,4),3),&
 &val4d2,grd(2)%masku,dble(0))

 !Filter
 IF(flag_filter) THEN
 END IF

 !Write to output
 DO i1=1,size(tOut)
  ibnd(1)=MAXLOC(time,1,time.LE.tOut(i1))
  ibnd(2)=MINLOC(time,1,time.GE.tOut(i1))
  IF(ibnd(1).EQ.ibnd(2)) THEN
   w=.5
  ELSE
   w=(tOut(i1)-time(ibnd(1)))/(time(ibnd(2))-time(ibnd(1)))
  END IF
 
  ncstart=[1,1,1,i1]
  nccount=[size(val4d2,1),size(val4d2,2),size(val4d2,3),&
  &1]
  CALL ncwrite4d(oid,'u',ncstart(1:4),nccount(1:4),&
  &(1.0-w)*val4d2(:,:,:,iBnd(1):iBnd(1))+&
  &w*val4d2(:,:,:,iBnd(2):iBnd(2)))
 END DO

 DEALLOCATE(val4d1,val4d2) 


!---------------------------------------------------------------------------
!Read v

 WRITE(*,*) 'Interpolating v'
 ALLOCATE(val4d1(grd(1)%s(1),grd(1)%s(2)-1,grd(1)%s(3),ntimes))
 ALLOCATE(val4d2(grd(2)%s(1),grd(2)%s(2)-1,grd(2)%s(3),ntimes))

 i4=1
 DO iHis=iHisBnd(1),iHisBnd(2)
  ncstart=[1,1,1,1]
  nccount=[grd(1)%s(1),grd(1)%s(2)-1,grd(1)%s(3),ntFile(iHis)]
  val4d1(:,:,:,i4:i4-1+ntFile(iHis))=&
  &ncread4d(nid(iHis),'v',ncstart(1:4),nccount(1:4))
  i4=i4+ntFile(iHis)
 END DO

 DO i0=0,np-1
  !Interpolate
  val4d2(:,:,:,itp(i0,1):itp(i0,2))=&
  &mesh2mesh_interp(grd(1)%lonv,grd(1)%latv,grd(1)%maskv,&
  &val4d1(:,:,:,itp(i0,1):itp(i0,2)),grd(2)%lonv,grd(2)%latv)
 ENd DO

 !Extrapolate
 CALL mesh2mesh_extrap(grd(2)%lonv,grd(2)%latv,&
 &.NOT.ANY(ANY(val4d2.EQ.interp_fill_value,4),3),&
 &val4d2,grd(2)%maskv,dble(0))

 !Filter
 IF(flag_filter) THEN
 END IF

 !Write to output
 DO i1=1,size(tOut)
  ibnd(1)=MAXLOC(time,1,time.LE.tOut(i1))
  ibnd(2)=MINLOC(time,1,time.GE.tOut(i1))
  IF(ibnd(1).EQ.ibnd(2)) THEN
   w=.5
  ELSE
   w=(tOut(i1)-time(ibnd(1)))/(time(ibnd(2))-time(ibnd(1)))
  END IF
 
  ncstart=[1,1,1,i1]
  nccount=[size(val4d2,1),size(val4d2,2),size(val4d2,3),&
  &1]
  CALL ncwrite4d(oid,'v',ncstart(1:4),nccount(1:4),&
  &(1.0-w)*val4d2(:,:,:,iBnd(1):iBnd(1))+&
  &w*val4d2(:,:,:,iBnd(2):iBnd(2)))
 END DO

 DEALLOCATE(val4d1,val4d2) 

!--------------------------------------------------------------------------
!Write time

 WRITE(*,*) 'Writing time'

 tOut=tOut-refShift*dble(86400)
 CALL ncwrite1d(oid,'ocean_time',[1],[size(tOut)],&
 &tOut)
 
!---------------------------------------------------------------------------
!Close

 !Close streams
 status=nf_close(oid)
 DO iHis=iHisBnd(1),iHisBnd(2)
  status=nf_close(nid(iHis))
 END DO

 WRITE(*,*) 'create_background DONE'
END PROGRAM create_background

