#undef demean_all
PROGRAM tl_sample
 USE mod_netcdf
 USE mod_interp
 USE mod_sample
 
 IMPLICIT NONE

 INTERFACE 
  FUNCTION demean(val_in,mask) RESULT(val_out)
  IMPLICIT NONE
  REAL(8),intent(in)::val_in(:)
  REAL(8)::val_out(lbound(val_in,1):ubound(val_in,1))
  LOGICAL,intent(in)::mask(:)
  REAL(8)::val_mean
  END FUNCTION demean
 END INTERFACE

 !-------------------------------------------------------------------------
 !DECLARE VARIABLES
 
 !Input
 character(len=1024)::commandfile,grdfile,hisfile,obsfile,outfile,&
 &infile,hisDir,inDir
 character(len=64)::varname
 character::char_lpfilter
 real(8)::tbnd(3)

 !Grid
 integer::gs(4),dimlen(8)
 real(8),allocatable,dimension(:,:)::lonr,latr,lonu,latu,lonv,latv,h
 logical,allocatable,dimension(:,:)::maskr,masku,maskv
 TYPE(sigma_param)::spar

 !Observations
 integer::os
 real(8),allocatable,dimension(:,:)::obs

 !Output
 real(8),allocatable::valm(:),timem(:),timein(:),timetmp(:)
 logical,allocatable::flag_todo(:),selection(:)
 integer::flag_file

 !Reading
 logical::flag_exist
 character(len=1024)::filename
 real(8),allocatable::wxy(:,:),wt(:,:),wz(:,:)
 integer,allocatable::ix(:,:),iy(:,:),iobs(:),it(:,:),iz(:,:)
 real(8),allocatable::his4d(:,:,:,:),his3d(:,:,:),tmp1(:) 
 integer::nSelect,nccount(4),ncstart(4),ncend(4),typeo(2),track
 integer::locstart(4),locend(4),loccount(4)
 real(8)::timeo,dt,fcut,tmp,zeta1,h1,track_mean
 real(8)::tfile(2),dl(3)
 real(8),allocatable::tmp4(:,:,:,:)

 !Netcdf
 integer::dimid(20),varid(20),status
 integer::ncidGrd,ncidCom,ncidObs,ncidOut,ncidIn(9999),ncidHis(9999)
 integer::i0,i1,i2,i3,i4,j0,j1,j2
 
 !Tidal frequency
 real(8)::f(8)
 f=[11.9672,12.000,12.4206,12.6583,23.9345,24.0659,25.8193,26.8684]*&
 &dble(3600) 
 f=dble(1)/f
  

 !----------------------------------------------------------------------
 ! Read input

 READ(*,*) !obslist file with location/time/type observations
 READ(*,'(A)') obsfile
 READ(*,*) !background file
 READ(*,'(A)') hisDir
 READ(*,'(A)') hisfile
 READ(*,*) !grid file
 READ(*,'(A)') grdfile
 READ(*,*) !input file from which sampling takes place
 READ(*,'(A)') inDir
 READ(*,'(A)') infile
 READ(*,*) !output file with sampled values from background file
 READ(*,'(A)') outfile
 READ(*,*) !Name variable in output file
 READ(*,'(A)') varname
 READ(*,*) !start time, end time, time shift (days)
 READ(*,*) tBnd(1),tBnd(2),tBnd(3)

 WRITE(*,*) 'observation file:',TRIM(obsfile)
 WRITE(*,*) 'grid file:',TRIM(grdfile)
 WRITE(*,*) 'background file:',TRIM(hisfile)
 WRITE(*,*) 'input file:',TRIM(infile)
 WRITE(*,*) 'output file:',TRIM(outfile)
 WRITE(*,*) 'output variable:',TRIM(varname)
 WRITE(*,*) 'reference time:',tBnd(1),tBnd(2),tBnd(3)

 tBnd=dble(tBnd*86400) !days->seconds

!-------------------------------------------------------------------------
!Read input and history files

 !Open stream to background files
 ncidHis=0
 IF(LEN_TRIM(hisDir).EQ.0) THEN
  status=nf_open(TRIM(hisfile),nf_nowrite,ncidHis(1))
  CALL ncsize(ncidHis(1),'ocean_time',dimlen)
  ALLOCATE(timetmp(dimlen(1)))
  timetmp=ncread1d(ncidHis(1),'ocean_time',[1],[dimlen(1)])

  ncstart(1)=MINLOC(timetmp,1,timetmp.GE.tBnd(1).AND.timetmp.LE.tBnd(2))
  ncend(1)=MAXLOC(timetmp,1,timetmp.GE.tBnd(1).AND.timetmp.LE.tBnd(2))
  nccount(1)=ncend(1)-ncstart(1)+1
  IF(nccount(1).GE.0) THEN
   ALLOCATE(timem(nccount(1))); timem=timetmp(ncstart(1):ncend(1))
  ELSE
   WRITE(*,*) 'tl_sample: no times in background file'; STOP
  END IF
  DEALLOCATE(timetmp)
 ELSE
  CALL col_nfopen(TRIM(hisDir),TRIM(hisfile),tBnd(1:2),ncidHis,timem)
 END IF
 IF(ALL(ncidHis.EQ.0)) THEN
  WRITE(*,*) 'tl_sample: requested time period not in files'; STOP
 END IF
 timem=timem+tBnd(3) 
 write(*,*) 'min/max time:',minval(timem),maxval(timem)

 !Open stream to input files
 ncidIn=0
 IF(LEN_TRIM(inDir).EQ.0) THEN
  status=nf_open(TRIM(infile),nf_nowrite,ncidIn(1))
  CALL ncsize(ncidIn(1),'ocean_time',dimlen)
  ALLOCATE(timetmp(dimlen(1)))
  timetmp=ncread1d(ncidIn(1),'ocean_time',[1],[dimlen(1)])

  ncstart(1)=MINLOC(timetmp,1,timetmp.GE.tBnd(1).AND.timetmp.LE.tBnd(2))
  ncend(1)=MAXLOC(timetmp,1,timetmp.GE.tBnd(1).AND.timetmp.LE.tBnd(2))
  nccount(1)=ncend(1)-ncstart(1)+1
  IF(nccount(1).GE.0) THEN
   ALLOCATE(timein(nccount(1))); timein=timetmp(ncstart(1):ncend(1))
  ELSE
   WRITE(*,*) 'tl_sample: no times in input file'; STOP
  END IF
  DEALLOCATE(timetmp)
 ELSE
  CALL col_nfopen(TRIM(inDir),TRIM(inFile),tBnd(1:2),ncidIn,timein)
 END IF 
 timein=timein+tBnd(3)


 IF(ANY(timem.NE.timein)) THEN
  WRITE(*,*) 'tl_sample: times in input file must be equal to times in&
  &background file'
  STOP
 END IF

 !-----------------------------------------------------------------------
 !READ GRID

 WRITE(*,*) 'reading grid file ',TRIM(grdfile)

 !Open stream
 status=nf_open(TRIM(grdfile),nf_nowrite,ncidGrd)
 
 !Get grid dimensions
 CALL ncsize(ncidGrd,'z0_r',dimlen)
 gs(1:3)=dimlen(1:3)
 gs(4)=size(timem)

 !Read horizontal coordinates
 ALLOCATE(lonr(gs(1),gs(2))); ALLOCATE(latr(gs(1),gs(2)))
 ALLOCATE(lonu(gs(1)-1,gs(2))); ALLOCATE(latu(gs(1)-1,gs(2)))
 ALLOCATE(lonv(gs(1),gs(2)-1)); ALLOCATE(latv(gs(1),gs(2)-1))
 ALLOCATE(maskr(gs(1),gs(2))); ALLOCATE(masku(gs(1)-1,gs(2)))
 ALLOCATE(maskv(gs(1),gs(2)-1))
 lonr=ncread2d(ncidGrd,'lon_rho',[1,1],gs(1:2))
 lonu=ncread2d(ncidGrd,'lon_u',[1,1],[gs(1)-1,gs(2)])
 lonv=ncread2d(ncidGrd,'lon_v',[1,1],[gs(1),gs(2)-1])
 latr=ncread2d(ncidGrd,'lat_rho',[1,1],gs(1:2))
 latu=ncread2d(ncidGrd,'lat_u',[1,1],[gs(1)-1,gs(2)])
 latv=ncread2d(ncidGrd,'lat_v',[1,1],[gs(1),gs(2)-1])
 maskr=(ncread2d(ncidGrd,'mask_rho',[1,1],[gs(1),gs(2)]).NE.0)
 masku=(ncread2d(ncidGrd,'mask_u',[1,1],[gs(1)-1,gs(2)]).NE.0)
 maskv=(ncread2d(ncidGrd,'mask_v',[1,1],[gs(1),gs(2)-1]).NE.0)

 !Read bottom
 ALLOCATE(h(gs(1),gs(2)))
 h=ncread2d(ncidGrd,'h',[1,1],gs(1:2))

 !Read s-grid parameters
 DO i0=1,size(ncidHis)
  IF(ncidHis(i0).EQ.0) CYCLE
  CALL get_sigma_param(ncidGrd,ncidHis(i0),spar,h)
  EXIT
 END DO

 !Close stream
 status=nf_close(ncidGrd)

 !---------------------------------------------------------------------------
 !Read obsfile

 write(*,*) 'reading obsfile',TRIM(obsfile)

 !Get size observation list
 status=nf_open(TRIM(obsfile),nf_nowrite,ncidObs)
 CALL ncsize(ncidObs,'time',dimlen); os=INT(dimlen(1)) 
 ALLOCATE(obs(os,7))
 WRITE(*,*) os,' entries in obsfile'

 !Read observations
 obs(:,1)=ncread1d(ncidObs,'lon',[1],[os])
 obs(:,2)=ncread1d(ncidObs,'lat',[1],[os])
 obs(:,3)=ncread1d(ncidObs,'z',[1],[os])
 obs(:,4)=ncread1d(ncidObs,'time',[1],[os])
 obs(:,5)=dble(ncread1d(ncidObs,'type',[1],[os]))
 obs(:,6)=ncread1d(ncidObs,'obs',[1],[os])
 CALL ncsize(ncidObs,'dir',dimlen)
 IF (dimlen(1).GT.0) THEN
  obs(:,7)=ncread1d(ncidObs,'dir',[1],[os])
  obs(:,7)=2*ACOS(-1.0)/360.*obs(:,7) !nautical deg->nautical rad
 END IF

 !Close stream
 status=nf_close(ncidObs) 

 !---------------------------------------------------------------------------
 !Create output

 write(*,*) 'creating output'

 !Array to store output
 ALLOCATE(valm(os)); valm=dble(0)
 ALLOCATE(flag_todo(os)); flag_todo=.FALSE. 

 OPEN(unit=102,iostat=flag_file,file=TRIM(outfile),status='old')
 IF(flag_file.EQ.0) CLOSE(102,status='delete')

 !Create output file
 status=nf_create(TRIM(outfile),OR(nf_noclobber,nf_classic_model),&
 &ncidOut)
 status=nf_def_dim(ncidOut,'K',nf_unlimited,dimid(1))
 status=nf_def_var(ncidOut,TRIM(varname),nf_double,1,dimid(1),varid(1))
 status=nf_close(ncidOut)
 
 !---------------------------------------------------------------------------
 ! Process type 2: daily-averaged surface zonal current

 typeo=[2,2] !observation type
 dt=dble(24*3600) !averaging period
 WRITE(*,*) 'type 2:'

 !Select entries to process
 flag_todo=INT(obs(:,5)).Ge.typeo(1).AND.INT(obs(:,5)).LE.typeo(2)
 IF(ANY(flag_todo)) THEN

 !Allocate arrays with weights and indices for this type
 nSelect=COUNT(flag_todo)
 WRITE(*,*) nSelect,' samples'
 ALLOCATE(iobs(nSelect))
 ALLOCATE(ix(4,nSelect)); ix=0
 ALLOCATE(iy(4,nSelect)); iy=0
 ALLOCATE(wxy(4,nSelect)); wxy=dble(0)
 ALLOCATE(it(gs(4),1)); it=0
 ALLOCATE(wt(gs(4),1)); wt=dble(0)
 
 !Get horizontal indices
 j0=1
 DO i0=1,os
  IF(flag_todo(i0)) THEN
   iobs(j0)=i0
   CALL sample_index2d(lonu(:,1),latu(1,:),masku,&
   &obs(i0,1),obs(i0,2),ix(:,j0),iy(:,j0),wxy(:,j0))
   j0=j0+1
  END IF
 END DO
 WRITE(*,*) 'min/max ix:',minval(ix),maxval(ix)
 WRITE(*,*) 'min/max iy:',minval(iy),maxval(iy)

 !Read history file
 ncstart=[MINVAL(ix),MINVAL(iy),gs(3),1]
 ncend=[MAXVAL(ix),MAXVAL(iy),gs(3),gs(4)]
 nccount=ncend-ncstart+[1,1,1,1]
 ALLOCATE(his4d(ncstart(1):ncend(1),ncstart(2):ncend(2),&
 &ncstart(3):ncend(3),ncstart(4):ncend(4)) )
 his4d=col_ncread4d(ncidIn,'u',ncstart,nccount,tBnd(1:2)) 
 WHERE(abs(his4d).GE.interp_fill_value); his4d=0.0; END WHERE

 !Get time indices
 DO WHILE(ANY(flag_todo))
  !time
  timeo=MINVAL(obs(:,4),1,flag_todo)

  !time indices and weight
  CALL sample_index_mean(timem,timeo,dt,it(:,1),wt(:,1))
      

  !Sample 
  DO j0=1,nSelect
   IF(obs(iobs(j0),4).EQ.timeo) THEN
    DO i1=MINVAL(it,it.GT.0),MAXVAL(it,it.GT.0)
    DO j1=1,size(ix,1)
     valm(iObs(j0))=valm(iObs(j0))+&
     &wt(i1,1)*wxy(j1,j0)*his4d(ix(j1,j0),iy(j1,j0),gs(3),&
     &it(i1,1))
    END DO !j1
    END DO !i1
    flag_todo(iObs(j0))=.FALSE.
   END IF
  END DO !j0
 
 END DO !do while

 DEALLOCATE(iobs,ix,iy,it,wxy,wt,his4d)
 flag_todo=INT(obs(:,5)).Ge.typeo(1).AND.INT(obs(:,5)).LE.typeo(2)
 WRITE(*,*) 'min/max samples:',MINVAL(valm,1,flag_todo),MAXVAL(valm,1,&
 &flag_todo)

 END IF !type 2

 !---------------------------------------------------------------------------
 ! Process type 3: daily-averaged surface meridional current

 typeo=[3,3] !observation type
 dt=dble(24*3600) !averaging period
 WRITE(*,*) 'type 3:'

 !Select entries to process
 flag_todo=INT(obs(:,5)).Ge.typeo(1).AND.INT(obs(:,5)).LE.typeo(2)
 IF(ANY(flag_todo)) THEN

 !Allocate arrays with weights and indices for this type
 nSelect=COUNT(flag_todo)
 WRITE(*,*) nSelect,' samples'
 ALLOCATE(iobs(nSelect))
 ALLOCATE(ix(4,nSelect)); ix=0
 ALLOCATE(iy(4,nSelect)); iy=0
 ALLOCATE(wxy(4,nSelect)); wxy=dble(0)
 ALLOCATE(it(gs(4),1)); it=0
 ALLOCATE(wt(gs(4),1)); wt=dble(0)
 
 !Get horizontal indices
 j0=1
 DO i0=1,os
  IF(flag_todo(i0)) THEN
   iobs(j0)=i0
   CALL sample_index2d(lonv(:,1),latv(1,:),maskv,&
   &obs(i0,1),obs(i0,2),ix(:,j0),iy(:,j0),wxy(:,j0))
   j0=j0+1
  END IF
 END DO
 WRITE(*,*) 'min/max ix:',minval(ix),maxval(ix)
 WRITE(*,*) 'min/max iy:',minval(iy),maxval(iy)

 !Read history file
 ncstart=[MINVAL(ix),MINVAL(iy),gs(3),1]
 ncend=[MAXVAL(ix),MAXVAL(iy),gs(3),gs(4)]
 nccount=ncend-ncstart+[1,1,1,1]
 ALLOCATE(his4d(ncstart(1):ncend(1),ncstart(2):ncend(2),&
 &ncstart(3):ncend(3),ncstart(4):ncend(4)) )
 his4d=col_ncread4d(ncidIn,'v',ncstart,nccount,tBnd(1:2))
 WHERE(ABS(his4d).GE.interp_fill_value); his4d=0.0; END WHERE

 !Get time indices
 DO WHILE(ANY(flag_todo))
  !time
  timeo=MINVAL(obs(:,4),1,flag_todo)

  !time indices and weight
  CALL sample_index_mean(timem,timeo,dt,it(:,1),wt(:,1))
    
  !Sample 
  DO j0=1,nSelect
   IF(obs(iobs(j0),4).EQ.timeo) THEN
    DO i1=MINVAL(it,it.GT.0),MAXVAL(it,it.GT.0)
    DO j1=1,size(ix,1)
     valm(iObs(j0))=valm(iObs(j0))+&
     &wt(i1,1)*wxy(j1,j0)*his4d(ix(j1,j0),iy(j1,j0),gs(3),&
     &it(i1,1))
    END DO !j1
    END DO !i1
    flag_todo(iObs(j0))=.FALSE.
   END IF
  END DO !j0
 
 END DO !do while

 DEALLOCATE(iobs,ix,iy,it,wxy,wt,his4d)
 flag_todo=INT(obs(:,5)).Ge.typeo(1).AND.INT(obs(:,5)).LE.typeo(2)
 WRITE(*,*) 'min/max samples:',MINVAL(valm,1,flag_todo),MAXVAL(valm,1,&
 &flag_todo)

 END IF !type 3
!---------------------------------------------------------------------------
! Process type 4: instant sea-surface temperature

 typeo=[4,4] !observation type
 WRITE(*,*) 'type 4:'

 !Select entries to process
 flag_todo=INT(obs(:,5)).Ge.typeo(1).AND.INT(obs(:,5)).LE.typeo(2)
 IF(ANY(flag_todo)) THEN

 !Allocate arrays with weights and indices for this type
 nSelect=COUNT(flag_todo)
 WRITE(*,*) nSelect,' samples'
 ALLOCATE(iobs(nSelect))
 ALLOCATE(ix(4,nSelect)); ix=0
 ALLOCATE(iy(4,nSelect)); iy=0
 ALLOCATE(wxy(4,nSelect)); wxy=dble(0)
 ALLOCATE(it(gs(4),1)); it=0
 ALLOCATE(wt(gs(4),1)); wt=dble(0)
 
 !Get horizontal indices
 j0=1
 DO i0=1,os
  IF(flag_todo(i0)) THEN
   iobs(j0)=i0
   CALL sample_index2d(lonr(:,1),latr(1,:),maskr,&
   &obs(i0,1),obs(i0,2),ix(:,j0),iy(:,j0),wxy(:,j0))
   j0=j0+1
  END IF
 END DO
 WRITE(*,*) 'min/max ix:',minval(ix),maxval(ix)
 WRITE(*,*) 'min/max iy:',minval(iy),maxval(iy)

 !Read history file
 ncstart=[MINVAL(ix),MINVAL(iy),gs(3),1]
 ncend=[MAXVAL(ix),MAXVAL(iy),gs(3),gs(4)]
 nccount=ncend-ncstart+[1,1,1,1]
 ALLOCATE(his4d(ncstart(1):ncend(1),ncstart(2):ncend(2),&
 &ncstart(3):ncend(3),ncstart(4):ncend(4)) )
 his4d=col_ncread4d(ncidIn,'temp',ncstart,nccount,tBnd(1:2))
 WHERE(ABS(his4d).GE.interp_fill_value); his4d=dble(0); END WHERE

 !Get time indices
 DO WHILE(ANY(flag_todo))
  !time
  timeo=MINVAL(obs(:,4),1,flag_todo)

  !time indices and weight
  CALL sample_index_delta(timem,timeo,it(:,1),wt(:,1))
    
  !Sample 
  DO j0=1,nSelect
   IF(obs(iobs(j0),4).EQ.timeo) THEN
    DO i1=MINVAL(it,it.GT.0),MAXVAL(it,it.GT.0)
    DO j1=1,size(ix,1)
     valm(iObs(j0))=valm(iObs(j0))+&
     &wt(i1,1)*wxy(j1,j0)*his4d(ix(j1,j0),iy(j1,j0),gs(3),&
     &it(i1,1))
    END DO !j1
    END DO !i1
    flag_todo(iObs(j0))=.FALSE.
   END IF
  END DO !j0
 
 END DO !do while

 DEALLOCATE(iobs,ix,iy,it,wxy,wt,his4d)
 flag_todo=INT(obs(:,5)).GE.typeo(1).AND.INT(obs(:,5)).LE.typeo(2)
 WRITE(*,*) 'min/max samples:',MINVAL(valm,1,flag_todo),MAXVAL(valm,1,&
 &flag_todo)

 END IF !type 4

!---------------------------------------------------------------------------
 ! Process type 5: daily-averaged sea-surface temperature

 typeo=[5,5] !observation type
 dt=dble(24*3600) !Averagin period
 dl=[.125,.125,0.0] !Half-width in zonal/meridional direction (degrees)

 !Select entries to process
 flag_todo=INT(obs(:,5)).Ge.typeo(1).AND.INT(obs(:,5)).LE.typeo(2)
 nSelect=COUNT(flag_todo)
 WRITE(*,*) 'type 5:',nSelect,' samples'
 IF(ANY(flag_todo)) THEN

 !Read history file
 ncstart=[1,1,gs(3),1]
 ncend=[gs(1),gs(2),gs(3),gs(4)]
 nccount=ncend-ncstart+[1,1,1,1]
 ALLOCATE(his4d(ncstart(1):ncend(1),ncstart(2):ncend(2),&
 &ncstart(3):ncend(3),ncstart(4):ncend(4)) )
 his4d=col_ncread4d(ncidIn,'temp',ncstart,nccount,tBnd(1:2))
 WHERE(ABS(his4d).GE.interp_fill_value); his4d=dble(0); END WHERE

 !Allocate time indices
 ALLOCATE(it(gs(4),1)); it=0
 ALLOCATE(wt(gs(4),1)); wt=dble(0)

 !Sample
 DO WHILE(ANY(flag_todo))

  !Get time indices
  timeo=MINVAL(obs(:,4),1,flag_todo)
  CALL sample_index_mean(timem,timeo,dt,it(:,1),wt(:,1))
    
  !Sample in horizontal
  DO j0=1,size(valm)
   IF(.NOT.flag_todo(j0).OR.obs(j0,4).NE.timeo) CYCLE

   !Horizontal weights
   CALL sample_avg2d(lonr(:,1),latr(1,:),maskr,obs(j0,1),obs(j0,2),&
   &dl(1:2),ix,iy,wxy)
  
   !Interpolate
   DO i4=MINVAL(it,it.GT.0),MAXVAL(it,it.GT.0)
   DO i2=1,size(wxy,2)
   DO i1=1,size(wxy,1)
    IF(ix(i1,i2).EQ.0.OR.iy(i1,i2).EQ.0) CYCLE
    valm(j0)=valm(j0)+&
    &wxy(i1,i2)*wt(i4,1)*&
    &his4d(ix(i1,i2),iy(i1,i2),gs(3),it(i4,1))
   END DO !i1
   END DO !i2
   END DO !i4
   
   !Deallocate locationn
   DEALLOCATE(wxy,ix,iy)
  
   !Mark as complete
   flag_todo(j0)=.FALSE.

  END DO !j0  
 END DO !do while

 DEALLOCATE(his4d,wt,it)
 flag_todo=INT(obs(:,5)).GE.typeo(1).AND.INT(obs(:,5)).LE.typeo(2)
 WRITE(*,*) 'min/max samples:',MINVAL(valm,1,flag_todo),MAXVAL(valm,1,&
 &flag_todo)

 END IF !type 5

 !---------------------------------------------------------------------------
 ! Process type 6: instant temperature

 typeo=[6,6] !observation type
 WRITE(*,*) 'type 6:'

 !Select entries to process
 flag_todo=INT(obs(:,5)).Ge.typeo(1).AND.INT(obs(:,5)).LE.typeo(2) 
 IF(ANY(flag_todo)) THEN

 !Allocate memory for indices
 nSelect=COUNT(flag_todo)
 ALLOCATE(selection(nSelect))
 ALLOCATE(iobs(nSelect)); iobs=0
 ALLOCATE(ix(4,nSelect)); ix=0
 ALLOCATE(iy(4,nSelect)); iy=0
 ALLOCATE(iz(2,nSelect)); iz=0
 ALLOCATE(it(2,nSelect)); it=0
 ALLOCATE(wxy(4,nSelect)); wxy=dble(0)
 ALLOCATE(wz(2,nSelect)); wz=dble(0)
 ALLOCATE(wt(2,nSelect)); wt=dble(0)
 
 !Horizontal indices, time indices
 j0=0
 DO i0=1,os
  IF(flag_todo(i0)) THEN
   j0=j0+1
   iobs(j0)=i0
   CALL sample_index2d(lonr(:,1),latr(1,:),maskr,obs(i0,1),obs(i0,2),&
   &ix(:,j0),iy(:,j0),wxy(:,j0))
   CALL sample_index_delta(timem,obs(i0,4),it(:,j0),wt(:,j0))
  END IF
 END DO
 WRITE(*,*) 'min/max ix:',minval(ix),maxval(ix)
 WRITE(*,*) 'min/max iy:',minval(iy),maxval(iy)
 WRITE(*,*) 'min/max it:',minval(it),maxval(it)

 !Get zetas and vertical indices/weights
 ncstart=[MINVAL(ix),MINVAL(iy),MINVAL(it),0]
 ncend=[MAXVAL(ix),MAXVAL(iy),MAXVAL(it),0]
 nccount=ncend-ncstart+[1,1,1,0]
 ALLOCATE( his3d(ncstart(1):ncend(1),ncstart(2):ncend(2),&
 &ncstart(3):ncend(3)) )
 his3d=col_ncread3d(ncidHis,'zeta',ncstart(1:3),nccount(1:3),tBnd(1:2))

 WHERE(ABS(his3d).GE.interp_fill_value); his3d=dble(0); END WHERE

 DO j0=1,size(iobs)
  zeta1=dble(0); h1=dble(0);
  DO i1=1,size(ix,1)
  DO i4=1,size(it,1)
   zeta1=zeta1+wt(i4,j0)*wxy(i1,j0)*&
   &his3d(ix(i1,j0),iy(i1,j0),it(i4,j0))
  END DO
  h1=h1+wxy(i1,j0)*h(ix(i1,j0),iy(i1,j0))
  END DO
  
  CALL sample_index_depth(spar,h1,zeta1,obs(iobs(j0),3),&
  &iz(:,j0),wz(:,j0))
 END DO
 DEALLOCATE(his3d)
 WRITE(*,*) 'min/max iz:',minval(iz),maxval(iz) 
 
 !Sample temp
 DO i1=MINVAL(ix),MAXVAL(ix)
 DO i2=MINVAL(iy),MAXVAL(iy)
  selection=ANY(ix.EQ.i1.AND.iy.EQ.i2,1)
  IF(.NOT.ANY(selection)) CYCLE  

  !Read temperature at (ix,iy)=(i1,i2)
  ncstart(1)=i1; ncstart(2)=i2
  ncstart(3)=MINVAL(iz(1,:),selection)
  ncstart(4)=MINVAL(it(1,:),selection)  

  ncend(1)=i1; ncend(2)=i2
  ncend(3)=MAXVAL(iz(2,:),selection)
  ncend(4)=MAXVAL(it(2,:),selection)
  nccount=ncend-ncstart+[1,1,1,1]

  ALLOCATE(his4d(ncstart(1):ncend(1),ncstart(2):ncend(2),&
  &ncstart(3):ncend(3),ncstart(4):ncend(4))); his4d=dble(0)
  his4d=col_ncread4d(ncidIn,'temp',ncstart,nccount,tBnd(1:2))
  WHERE(ABS(his4d).GE.interp_fill_value); his4d=dble(0); END WHERE

  !sample
  DO j0=1,size(iobs)
   IF(.NOT.selection(j0)) CYCLE
  DO i0=1,size(wxy,1)
  DO i3=1,size(wz,1)
  DO i4=1,size(wt,1)
   IF(ix(i0,j0).EQ.i1.AND.iy(i0,j0).EQ.i2) THEN
    valm(iobs(j0))=valm(iobs(j0))+&
    wxy(i0,j0)*wz(i3,j0)*wt(i4,j0)*&
    &his4d(i1,i2,iz(i3,j0),it(i4,j0))
   END IF
  END DO !i4
  END DO !i3
  END DO !i0
  END DO !j0

  DEALLOCATE(his4d)

 END DO !i2
 END DO !i1

 DEALLOCATE(iobs,ix,iy,iz,it,wxy,wt,wz,selection)
 ALLOCATE(selection(os))
 selection=(INT(obs(:,5)).GE.typeo(1).AND.INT(obs(:,5)).LE.typeo(2))
 WRITE(*,*) 'min/max samples:',MINVAL(valm,1,selection),MAXVAL(valm,1,&
 &selection)
 DEALLOCATE(selection)
 
 END IF !type 6
 !---------------------------------------------------------------------------
 ! Process type 7: instant salinity

 typeo=[7,7] !observation type
 WRITE(*,*) 'type 7:'

 !Select entries to process
 flag_todo=INT(obs(:,5)).Ge.typeo(1).AND.INT(obs(:,5)).LE.typeo(2) 
 IF(ANY(flag_todo)) THEN

 !Allocate memory for indices
 nSelect=COUNT(flag_todo)
 ALLOCATE(selection(nSelect))
 ALLOCATE(iobs(nSelect)); iobs=0
 ALLOCATE(ix(4,nSelect)); ix=0
 ALLOCATE(iy(4,nSelect)); iy=0
 ALLOCATE(iz(2,nSelect)); iz=0
 ALLOCATE(it(2,nSelect)); it=0
 ALLOCATE(wxy(4,nSelect)); wxy=dble(0)
 ALLOCATE(wz(2,nSelect)); wz=dble(0)
 ALLOCATE(wt(2,nSelect)); wt=dble(0)
 
 !Horizontal indices, time indices
 j0=0
 DO i0=1,os
  IF(flag_todo(i0)) THEN
   j0=j0+1
   iobs(j0)=i0
   CALL sample_index2d(lonr(:,1),latr(1,:),maskr,obs(i0,1),obs(i0,2),&
   &ix(:,j0),iy(:,j0),wxy(:,j0))
   CALL sample_index_delta(timem,obs(i0,4),it(:,j0),wt(:,j0))
  END IF
 END DO
 WRITE(*,*) 'min/max ix:',minval(ix),maxval(ix)
 WRITE(*,*) 'min/max iy:',minval(iy),maxval(iy)
 WRITE(*,*) 'min/max it:',minval(it),maxval(it)

 !Get zetas and vertical indices/weights
 ncstart=[MINVAL(ix),MINVAL(iy),MINVAL(it),0]
 ncend=[MAXVAL(ix),MAXVAL(iy),MAXVAL(it),0]
 nccount=ncend-ncstart+[1,1,1,0]
 ALLOCATE( his3d(ncstart(1):ncend(1),ncstart(2):ncend(2),&
 &ncstart(3):ncend(3)) )
 his3d=col_ncread3d(ncidHis,'zeta',ncstart(1:3),nccount(1:3),&
 &tBnd(1:2))
 WHERE(ABS(his3d).GE.interp_fill_value); his3d=dble(0); END WHERE

 DO j0=1,size(iobs)
  zeta1=dble(0); h1=dble(0);
  DO i1=1,size(ix,1)
  DO i4=1,size(it,1)
   zeta1=zeta1+wt(i4,j0)*wxy(i1,j0)*&
   &his3d(ix(i1,j0),iy(i1,j0),it(i4,j0))
  END DO
  h1=h1+wxy(i1,j0)*h(ix(i1,j0),iy(i1,j0))
  END DO
  
  CALL sample_index_depth(spar,h1,zeta1,obs(iobs(j0),3),&
  &iz(:,j0),wz(:,j0))
 END DO
 DEALLOCATE(his3d)
 WRITE(*,*) 'min/max iz:',minval(iz),maxval(iz) 
 
 !Sample temp
 DO i1=MINVAL(ix),MAXVAL(ix)
 DO i2=MINVAL(iy),MAXVAL(iy)
  selection=ANY(ix.EQ.i1.AND.iy.EQ.i2,1)
  IF(.NOT.ANY(selection)) CYCLE 
 
  !Read temperature at (ix,iy)=(i1,i2)
  ncstart(1)=i1; ncstart(2)=i2
  ncstart(3)=MINVAL(iz(1,:),selection)
  ncstart(4)=MINVAL(it(1,:),selection)  

  ncend(1)=i1; ncend(2)=i2
  ncend(3)=MAXVAL(iz(2,:),selection)
  ncend(4)=MAXVAL(it(2,:),selection)
  nccount=ncend-ncstart+[1,1,1,1]

  ALLOCATE(his4d(ncstart(1):ncend(1),ncstart(2):ncend(2),&
  &ncstart(3):ncend(3),ncstart(4):ncend(4)))
  his4d=col_ncread4d(ncidIn,'salt',ncstart,nccount,tBnd(1:2))
  WHERE(ABS(his4d).GE.interp_fill_value); his4d=dble(0); END WHERE

  !sample
  DO j0=1,size(iobs)
   IF(.NOT.selection(j0)) CYCLE
  DO i0=1,size(wxy,1)
  DO i3=1,size(wz,1)
  DO i4=1,size(wt,1)
   IF(ix(i0,j0).EQ.i1.AND.iy(i0,j0).EQ.i2) THEN
    valm(iobs(j0))=valm(iobs(j0))+&
    wxy(i0,j0)*wz(i3,j0)*wt(i4,j0)*&
    &his4d(i1,i2,iz(i3,j0),it(i4,j0))
   END IF
  END DO !i4
  END DO !i3
  END DO !i0
  END DO !j0

  DEALLOCATE(his4d)

 END DO !i2
 END DO !i1

 DEALLOCATE(iobs,ix,iy,iz,it,wxy,wt,wz,selection)
 ALLOCATE(selection(os))
 selection=(INT(obs(:,5)).GE.typeo(1).AND.INT(obs(:,5)).LE.typeo(2))
 WRITE(*,*) 'min/max samples:',MINVAL(valm,1,selection),MAXVAL(valm,1,&
 &selection)
 DEALLOCATE(selection)
 
 END IF !type 7

!--------------------------------------------------------------------------
! Process type 8: daily-averaged radial velocity

 typeo=[8,8] !observation type          
 dt=dble(24*3600)                                    
 WRITE(*,*) 'type 8:'

 !Select entries to process                                                  
 flag_todo=INT(obs(:,5)).Ge.typeo(1).AND.INT(obs(:,5)).LE.typeo(2)
 IF(ANY(flag_todo)) THEN

 !Allocate arrays with weights and indices for this type                       
 nSelect=COUNT(flag_todo)
 WRITE(*,*) nSelect,' samples'
 ALLOCATE(iobs(nSelect))
 ALLOCATE(ix(4,nSelect)); ix=0
 ALLOCATE(iy(4,nSelect)); iy=0
 ALLOCATE(wxy(4,nSelect)); wxy=dble(0)
 ALLOCATE(it(gs(4),1)); it=0
 ALLOCATE(wt(gs(4),1)); wt=dble(0)

 !Get horizontal indices U                                                     
 j0=1
 DO i0=1,os
  IF(flag_todo(i0)) THEN
   iobs(j0)=i0
   CALL sample_index2d(lonu(:,1),latu(1,:),masku,&
   &obs(i0,1),obs(i0,2),ix(:,j0),iy(:,j0),wxy(:,j0))
   j0=j0+1
  END IF
 END DO
 WRITE(*,*) 'min/max ix:',minval(ix),maxval(ix)
 WRITE(*,*) 'min/max iy:',minval(iy),maxval(iy)

 !Read history file U                                                          
 ncstart=[MINVAL(ix),MINVAL(iy),gs(3),1]
 ncend=[MAXVAL(ix),MAXVAL(iy),gs(3),gs(4)]
 nccount=ncend-ncstart+[1,1,1,1]
 ALLOCATE(his4d(ncstart(1):ncend(1),ncstart(2):ncend(2),&
 &ncstart(3):ncend(3),ncstart(4):ncend(4)) )
 his4d=col_ncread4d(ncidIn,'u',ncstart,nccount,tBnd(1:2))
 WHERE(abs(his4d).GE.interp_fill_value); his4d=0.0; END WHERE

 !Get time indices U                                                           
 DO WHILE(ANY(flag_todo))
  !time                                                                        
  timeo=MINVAL(obs(:,4),1,flag_todo)

  !time indices and weight U                
  CALL sample_index_mean(timem,timeo,dt,it(:,1),wt(:,1))

  !Sample U                                                                   
  DO j0=1,nSelect
   IF(obs(iobs(j0),4).EQ.timeo) THEN
    DO i1=MINVAL(it,it.GT.0),MAXVAL(it,it.GT.0)
    DO j1=1,size(ix,1)
     valm(iObs(j0))=valm(iObs(j0))+&
     &SIN(obs(iObs(j0),7))*wt(i1,1)*wxy(j1,j0)*&
     &his4d(ix(j1,j0),iy(j1,j0),gs(3),it(i1,1))
    END DO !j1                                                                 
    END DO !i1                                                                 
    flag_todo(iObs(j0))=.FALSE.
   END IF
  END DO !j0                                                                   
 END DO !do while

 !Allocate for V
 DEALLOCATE(iobs,ix,iy,it,wxy,wt,his4d)
 flag_todo=INT(obs(:,5)).Ge.typeo(1).AND.INT(obs(:,5)).LE.typeo(2)
 nSelect=COUNT(flag_todo)

 ALLOCATE(iobs(nSelect))
 ALLOCATE(ix(4,nSelect)); ix=0
 ALLOCATE(iy(4,nSelect)); iy=0
 ALLOCATE(wxy(4,nSelect)); wxy=dble(0)
 ALLOCATE(it(gs(4),1)); it=0
 ALLOCATE(wt(gs(4),1)); wt=dble(0)

 !Get horizontal indices V                                                    
 j0=1
 DO i0=1,os
  IF(flag_todo(i0)) THEN
   iobs(j0)=i0
   CALL sample_index2d(lonv(:,1),latv(1,:),maskv,&
   &obs(i0,1),obs(i0,2),ix(:,j0),iy(:,j0),wxy(:,j0))
   j0=j0+1
  END IF
 END DO
 WRITE(*,*) 'min/max ix:',minval(ix),maxval(ix)
 WRITE(*,*) 'min/max iy:',minval(iy),maxval(iy)

 !Read history file V                                                          
 ncstart=[MINVAL(ix),MINVAL(iy),gs(3),1]
 ncend=[MAXVAL(ix),MAXVAL(iy),gs(3),gs(4)]
 nccount=ncend-ncstart+[1,1,1,1]
 ALLOCATE(his4d(ncstart(1):ncend(1),ncstart(2):ncend(2),&
 &ncstart(3):ncend(3),ncstart(4):ncend(4)) )
 his4d=col_ncread4d(ncidIn,'v',ncstart,nccount,tBnd(1:2))
 WHERE(abs(his4d).GE.interp_fill_value); his4d=0.0; END WHERE

 !Get time indices V                                                          
 DO WHILE(ANY(flag_todo))
  !time                                                                        
  timeo=MINVAL(obs(:,4),1,flag_todo)

  !time indices and weight V                                           
  CALL sample_index_mean(timem,timeo,dt,it(:,1),wt(:,1))

  !Sample V                                                                   
  DO j0=1,nSelect
   IF(obs(iobs(j0),4).EQ.timeo) THEN
    DO i1=MINVAL(it,it.GT.0),MAXVAL(it,it.GT.0)
    DO j1=1,size(ix,1)
     valm(iObs(j0))=valm(iObs(j0))+&
     &COS(obs(iObs(j0),7))*wt(i1,1)*wxy(j1,j0)*&
     &his4d(ix(j1,j0),iy(j1,j0),gs(3),it(i1,1))
    END DO !j1                                                                 
    END DO !i1                                                                 
    flag_todo(iObs(j0))=.FALSE.
   END IF
  END DO !j0   
 END DO !do while      

 DEALLOCATE(iobs,ix,iy,it,wxy,wt,his4d)
 flag_todo=INT(obs(:,5)).Ge.typeo(1).AND.INT(obs(:,5)).LE.typeo(2)
 WRITE(*,*) 'min/max samples:',MINVAL(valm,1,flag_todo),MAXVAL(valm,1,&
 &flag_todo)

 END IF !type 8         
!---------------------------------------------------------------------------
! Process type 100-199: average salinity

 typeo=[100,199]
 WRITE(*,*) 'type 100-199:'
 !dl=[.6,1.5,25.]
 !-124.28 46.60 2.55

 !Select entries to process
 flag_todo=INT(obs(:,5)).Ge.typeo(1).AND.INT(obs(:,5)).LE.typeo(2) 
 IF(ANY(flag_todo)) THEN

 !Allocate arrays with weights and indices for this type                        
 nSelect=COUNT(flag_todo)
 WRITE(*,*) nSelect,' samples'
 ALLOCATE(it(gs(4),1)); it=0
 ALLOCATE(wt(gs(4),1)); wt=dble(0)
                                                              
 DO WHILE(ANY(flag_todo))
  !time                                                                         
  timeo=MINVAL(obs(:,4),1,flag_todo)

  !time indices and weight                                                      
  CALL sample_index_delta(timem,timeo,it(:,1),wt(:,1))
  WRITE(*,*) 'Min/max it:',MINVAL(it,it.GT.0),MAXVAL(it,it.GT.0)

  !Read zeta                                                                   
  ncstart=[1,1,MINVAL(it,it.GT.0),0]
  ncend=[gs(1),gs(2),MAXVAL(it,it.GT.0),0]
  nccount=ncend-ncstart+[1,1,1,0]
  ALLOCATE( his3d(ncstart(1):ncend(1),ncstart(2):ncend(2),&
  &ncstart(3):ncend(3)) )
  his3d=col_ncread3d(ncidHis,'zeta',ncstart(1:3),nccount(1:3),tBnd(1:2))

  !Read salinity field
  ncstart=[1,1,1,MINVAL(it,it.GT.0)]
  ncend=[gs(1),gs(2),gs(3),MAXVAL(it,it.GT.0)]
  nccount=ncend-ncstart+[1,1,1,1]
  ALLOCATE(his4d(ncstart(1):ncend(1),ncstart(2):ncend(2),&
  &ncstart(3):ncend(3),ncstart(4):ncend(4)))
  his4d=col_ncread4d(ncidIn,'salt',ncstart,nccount,tBnd(1:2))

  DO j0=1,os
   IF(.NOT.flag_todo(j0)) CYCLE
   IF(obs(j0,4).NE.timeo) CYCLE

   !Size based on sample on 4km grid
   dl(1)=.5*6.0610/2.0**(obs(j0,5)-100)
   dl(2)=.5*9.3106/2.0**(obs(j0,5)-100)
   dl(3)=1.0

   !Horizontal weights                                                          
   CALL sample_avg2d(lonr(:,1),latr(1,:),maskr,obs(j0,1),obs(j0,2),&
   &dl(1:2),ix,iy,wxy)

   !Sample                                                                       
   DO i4=MINVAL(it,it.GT.0),MAXVAL(it,it.GT.0)
   DO i2=1,size(wxy,2)
   DO i1=1,size(wxy,1)

    IF(ix(i1,i2).EQ.0.OR.iy(i1,i2).EQ.0) THEN
     flag_todo(j0)=.FALSE.
     CYCLE
    END IF

    !Vertical weights                                                           
    CALL sample_depthAvg(spar,h(ix(i1,i2),iy(i1,i2)),&
    &his3d(ix(i1,i2),iy(i1,i2),it(i4,1)),obs(j0,3),dl(3),iz,wz)

    DO i3=MINVAL(iz,iz.GT.0),MAXVAL(iz,iz.GT.0)
     !Multiply salinity with weights and add                                    
     valm(j0)=valm(j0)+wz(1,i3)*wt(i4,1)*wxy(i1,i2)&
     &*his4d(ix(i1,i2),iy(i1,i2),iz(1,i3),it(i4,1))
    END DO
    DEALLOCATE(iz,wz)

   END DO !i1                                                                   
   END DO !i2                                                                   
   END DO !i4                                                                   
   DEALLOCATE(ix,iy,wxy)

   !Mark as complete
   flag_todo(j0)=.FALSE.

  END DO !j0
  DEALLOCATE(his3d,his4d)
 
 END DO !while
 
 DEALLOCATE(it,wt)
 flag_todo=INT(obs(:,5)).Ge.typeo(1).AND.INT(obs(:,5)).LE.typeo(2)
 WRITE(*,*) 'min/max samples:',MINVAL(valm,1,flag_todo),MAXVAL(valm,1,&
 &flag_todo)

 END IF !type 100-199

 !---------------------------------------------------------------------------
 ! Process type 1000-9999: detided sea-surface height

 typeo=[1000,9999] !observation type
 dt=dble(24*3600) !time length filter 
 WRITE(*,*) 'type 1000-9999:'

 DO track=typeo(1),typeo(2)

 !select observations to process
 flag_todo=INT(obs(:,5)).EQ.track
 IF(ANY(flag_todo)) THEN
 WRITE(*,*) 'track:',track
 
 !Allocate arrays with weights and indices for this type
 nSelect=COUNT(flag_todo)
 WRITE(*,*) nSelect,' samples'
 ALLOCATE(iobs(nSelect))
 ALLOCATE(ix(4,nSelect)); ix=0
 ALLOCATE(iy(4,nSelect)); iy=0
 ALLOCATE(wxy(4,nSelect)); wxy=dble(0)
 ALLOCATE(it(gs(4),1)); it=0
 ALLOCATE(wt(gs(4),1)); wt=dble(0)

 !Get horizontal indices
 j0=1
 DO i0=1,os
  IF(flag_todo(i0)) THEN
   iobs(j0)=i0
   CALL sample_index2d(lonr(:,1),latr(1,:),maskr,&
   &obs(i0,1),obs(i0,2),ix(:,j0),iy(:,j0),wxy(:,j0))
   j0=j0+1
  END IF
 END DO
 WRITE(*,*) 'min/max ix:',minval(ix),maxval(ix)
 WRITE(*,*) 'min/max iy:',minval(iy),maxval(iy)

 !Read history file
 ncstart=[MINVAL(ix),MINVAL(iy),1,0]
 ncend=[MAXVAL(ix),MAXVAL(iy),gs(4),0]
 nccount=ncend-ncstart+[1,1,1,0]

 ALLOCATE(his3d(ncstart(1):ncend(1),ncstart(2):ncend(2),&
 &ncstart(3):ncend(3)) )
 his3d=col_ncread3d(ncidIn,'zeta',ncstart(1:3),nccount(1:3),tBnd(1:2)) 
 WHERE(ABS(his3d).GE.interp_fill_value); his3d=dble(0); END WHERE

 !Get time indices
 DO WHILE(ANY(flag_todo))
  !time
  timeo=MINVAL(obs(:,4),1,flag_todo)

  !time indices and weight
  CALL sample_index_mean(timem,timeo,dt,it(:,1),wt(:,1)) 

  !Sample 
  DO j0=1,nSelect
   IF(obs(iobs(j0),4).EQ.timeo) THEN
    DO i1=MINVAL(it,it.GT.0),MAXVAL(it,it.GT.0)
    DO j1=1,size(ix,1)
     valm(iObs(j0))=valm(iObs(j0))+&
     &wt(i1,1)*wxy(j1,j0)*his3d(ix(j1,j0),iy(j1,j0),&
     &it(i1,1))
    END DO !j1
    END DO !i1
    flag_todo(iObs(j0))=.FALSE.
   END IF
  END DO !j0
 
 END DO !do while

 !demean
 valm=demean(valm,INT(obs(:,5)).EQ.track)

 DEALLOCATE(iobs,ix,iy,it,wxy,wt,his3d)
 flag_todo=INT(obs(:,5)).EQ.track
 WRITE(*,*) 'min/max samples:',MINVAL(valm,1,flag_todo),MAXVAL(valm,1,&
 &flag_todo)

 END IF !type.EQ.track
 END DO !track
!---------------------------------------------------------------------------
! Process type 10000-19999: instant sea-surface temperature

 typeo=[10000,19999] !observation type
 WRITE(*,*) 'type 10000-19999:'

 DO track=typeo(1),typeo(2)
  
  !Select entries to process
  flag_todo=INT(obs(:,5)).EQ.track
  IF(.NOT.ANY(flag_todo)) CYCLE
  WRITE(*,*) 'track:',track
 
  !Allocate arrays with weights and indices for this type
  nSelect=COUNT(flag_todo)
  WRITE(*,*) nSelect,' samples'
  ALLOCATE(iobs(nSelect))
  ALLOCATE(ix(4,nSelect)); ix=0
  ALLOCATE(iy(4,nSelect)); iy=0
  ALLOCATE(wxy(4,nSelect)); wxy=dble(0)
  ALLOCATE(it(gs(4),1)); it=0
  ALLOCATE(wt(gs(4),1)); wt=dble(0)
  ALLOCATE(tmp1(nSelect)); tmp1=dble(0)  

  !Get horizontal indices
  j0=1
  DO i0=1,os
   IF(flag_todo(i0)) THEN
    iobs(j0)=i0
    CALL sample_index2d(lonr(:,1),latr(1,:),maskr,&
    &obs(i0,1),obs(i0,2),ix(:,j0),iy(:,j0),wxy(:,j0))
    
    j0=j0+1
   END IF
  END DO
  WRITE(*,*) 'min/max ix:',minval(ix),maxval(ix)
  WRITE(*,*) 'min/max iy:',minval(iy),maxval(iy)

  !Read history file
  ncstart=[MINVAL(ix),MINVAL(iy),gs(3),1]
  ncend=[MAXVAL(ix),MAXVAL(iy),gs(3),gs(4)]
  nccount=ncend-ncstart+[1,1,1,1]
  ALLOCATE(his4d(ncstart(1):ncend(1),ncstart(2):ncend(2),&
  &ncstart(3):ncend(3),ncstart(4):ncend(4)) )
  his4d=col_ncread4d(ncidIn,'temp',ncstart,nccount,tBnd(1:2))
  WHERE(ABS(his4d).GE.interp_fill_value); his4d=dble(0); END WHERE

  !Get time indices
  DO WHILE(ANY(flag_todo))
   !time
   timeo=MINVAL(obs(:,4),1,flag_todo)
   CALL sample_index_delta(timem,obs(i0,4),it(:,1),wt(:,1))    

   !Sample 
   DO j0=1,nSelect
    IF(obs(iobs(j0),4).NE.timeo) CYCLE
    
    DO i1=MINVAL(it,it.GT.0),MAXVAL(it,it.GT.0)
    DO j1=1,size(ix,1)
     tmp1(j0)=tmp1(j0)+&
     &wt(i1,1)*wxy(j1,j0)*his4d(ix(j1,j0),iy(j1,j0),gs(3),&
     &it(i1,1))
    END DO !j1
    END DO !i1
    
    flag_todo(iObs(j0))=.FALSE.
   END DO !j0
  END DO !do while

  !Seperate in mean and deviations from mean
  DO j0=1,nSelect
   valm(iObs(1))=valm(iObs(1))+tmp1(j0)/dble(nSelect)
  END DO
  DO i0=2,nSelect
   valm(iObs(i0))=valm(iObs(i0))+tmp1(i0)
   DO j0=1,nSelect
    valm(iObs(i0))=valm(iObs(i0))-tmp1(j0)/dble(nSelect)
   END DO
  END DO 
  WRITE(*,*) 'mean track:',valm(iObs(1))

  DEALLOCATE(iobs,ix,iy,it,wxy,wt,his4d,tmp1)
 
 END DO !type 10000-19999

!--------------------------------------------------------------------------
!Demean all types of observations

#ifdef demean_all
 DO i0=MINVAL(obs(:,5)),MAXVAL(obs(:,5)) 
  IF(.NOT.ANY(obs(:,5).EQ.i0)) CYCLE
  valm=demean(valm,obs(:,5).EQ.i0)
 END DO
#endif

!---------------------------------------------------------------------------

 !Close streams
 DO i0=1,size(ncidHis)
  IF(ncidHis(i0).EQ.0) CYCLE
  status=nf_close(ncidHis(i0))
 END DO
 DO i0=1,size(ncidIn)
  IF(ncidIn(i0).EQ.0) CYCLE
  status=nf_close(ncidIn(i0))
 END DO
 !--------------------------------------------------------------------------
 
 !Open stream
 status=nf_open(TRIM(outfile),nf_write,ncidOut)

 !Write output to output file
 WRITE(*,*) 'writing samples to output'
 ncstart=[1,0,0,0]
 nccount=[os,0,0,0]
 CALL ncwrite1d(ncidOut,TRIM(varname),ncstart(1),nccount(1),valm)

 !Close stream
 status=nf_close(ncidOut)

 !-------------------------------------------------------------------------
 WRITE(*,*) 'DONE'
END PROGRAM tl_sample

!-------------------------------------------------------------------------

FUNCTION demean(val_in,mask) RESULT(val_out)
!val_out=demean(val_in,mask) removes mean from subset mask of val_in
!
!val_in [1D double]: input values
!val_out [1D double]: output values with mean removed
!mask [1D logical]: subset from which mean should be removed

 IMPLICIT NONE
 REAL(8),intent(in)::val_in(:)
 REAL(8)::val_out(lbound(val_in,1):ubound(val_in,1))
 LOGICAL,intent(in)::mask(:)
 REAL(8)::val_mean

 val_out=val_in
 IF(ANY(mask)) THEN
  val_mean=SUM(val_in,1,mask)/dble(count(mask))
  WHERE(mask)
   val_out=val_out-val_mean
  END WHERE
 END IF

END FUNCTION demean
