#undef demean_all
PROGRAM ad_sample
 USE mod_netcdf
 USE mod_interp
 USE mod_roms
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
 &impulsfile
 character(len=64)::varname
 character::char_lpfilter

 !Grid
 integer::gs(4),dimlen(8)
 real(8),allocatable,dimension(:,:)::lonr,latr,lonu,latu,lonv,latv,h
 logical,allocatable,dimension(:,:)::maskr,masku,maskv
 TYPE(sigma_param)::spar

 !Observations
 integer::os
 real(8),allocatable,dimension(:,:)::obs

 !Output
 TYPE(roms_type)::fdesign
 real(8),allocatable::valm(:),timem(:)
 logical,allocatable::flag_todo(:),selection(:)
 integer::flag_file

 !Reading
 real(8),allocatable::wxy(:,:),wt(:,:),wz(:,:)
 integer,allocatable::ix(:,:),iy(:,:),iobs(:),it(:,:),iz(:,:)
 real(8),allocatable::his4d(:,:,:,:),his3d(:,:,:)
 real(8),allocatable,dimension(:,:,:,:)::temp,salt,u,v
 real(8),allocatable,dimension(:,:,:)::zeta,ubar,vbar 
 integer::nSelect,nccount(4),ncstart(4),ncend(4),typeo(2),track
 real(8)::timeo,dt,fcut,tmp,zeta1,h1,track_mean,dl(3)
 real(8),allocatable::tmp1(:)

 !Netcdf
 integer::dimid(20),varid(20),status
 integer::ncidGrd,ncidCom,ncidHis,ncidObs,ncidOut,ncidImpuls
 integer::i0,i1,i2,i3,i4,j0,j1,j2

 !----------------------------------------------------------------------
 ! Read input

 READ(*,*) !obslist file with location/time/type observations
 READ(*,'(A)') obsfile
 READ(*,*) !background file
 READ(*,'(A)') hisfile
 READ(*,*) !grid file
 READ(*,'(A)') grdfile
 READ(*,*) !Impulse file
 READ(*,'(A)') impulsfile
 READ(*,*) !Name variable in output file
 READ(*,'(A)') varname
 READ(*,*) !output file with sampled values from background file
 READ(*,'(A)') outfile

 WRITE(*,*) 'observation file:',TRIM(obsfile)
 WRITE(*,*) 'grid file:',TRIM(grdfile)
 WRITE(*,*) 'background file:',TRIM(hisfile)
 WRITE(*,*) 'output file:',TRIM(outfile)
 WRITE(*,*) 'input file:',TRIM(impulsfile)
 WRITE(*,*) 'input variable:',TRIM(varname)


 !-----------------------------------------------------------------------
 !READ GRID

 WRITE(*,*) 'reading grid file ',TRIM(grdfile)

 !Open stream
 status=nf_open(TRIM(grdfile),nf_nowrite,ncidGrd)
 
 !Get grid dimensions
 CALL ncsize(ncidGrd,'z0_r',dimlen)
 gs(1:3)=dimlen(1:3)

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
 status=nf_open(TRIM(hisfile),nf_nowrite,ncidHis)
 CALL get_sigma_param(ncidGrd,ncidHis,spar,h)

 !Close stream
 status=nf_close(ncidGrd)
 status=nf_close(ncidHis)

 !---------------------------------------------------------------------------
 !Read obsfile

 write(*,*) 'reading obsfile ',TRIM(obsfile)

 !Get size observation list
 status=nf_open(TRIM(obsfile),nf_nowrite,ncidObs)
 CALL ncsize(ncidObs,'time',dimlen); os=INT(dimlen(1)) 
 ALLOCATE(obs(os,7))
 WRITE(*,*) os,' entries in obsfile'

 !Read observations
 obs(:,4)=ncread1d(ncidObs,'time',[1],[os])
 obs(:,1)=ncread1d(ncidObs,'lon',[1],[os])
 obs(:,2)=ncread1d(ncidObs,'lat',[1],[os])
 obs(:,3)=ncread1d(ncidObs,'z',[1],[os])
 obs(:,5)=dble(ncread1d(ncidObs,'type',[1],[os]))
 obs(:,6)=ncread1d(ncidObs,'obs',[1],[os])
 CALL ncsize(ncidObs,'dir',dimlen)
 IF (dimlen(1).GT.0) THEN
  obs(:,7)=ncread1d(ncidObs,'dir',[1],[os])
  obs(:,7)=2*ACOS(-1.0)/360.*obs(:,7) !nautical deg->nautical rad              
 END IF

 !todo mask
 ALLOCATE(flag_todo(os))

 !Close stream
 status=nf_close(ncidObs) 

 !------------------------------------------------------------------------
 !Read input

 !Open stream
 status=nf_open(TRIM(impulsfile),nf_nowrite,ncidImpuls)

 ALLOCATE(valm(os))
 valm=ncread1d(ncidImpuls,TRIM(varname),[1],[os])

 !close stream
 status=nf_close(ncidImpuls)
 !---------------------------------------------------------------------------
 !Create output

 write(*,*) 'creating output'

 !Open stream to background file
 status=nf_open(TRIM(hisfile),nf_nowrite,ncidHis)

 !Create file
 CALL roms_get_type(TRIM(hisfile),fdesign)
 CALL roms_create_his_file(TRIM(outfile),fdesign)

 !Read time from history file
 CALL ncsize(ncidHis,'ocean_time',dimlen); gs(4)=dimlen(1)
 ALLOCATE(timem(gs(4)))
 timem=ncread1d(ncidHis,'ocean_time',[1],[gs(4)])
 write(*,*) gs(4),' times in ',TRIM(hisfile)

 !Open stream
 status=nf_open(TRIM(outfile),nf_write,ncidOut)

 !Initiate output
  write(*,*) 'write AKt,AKv'
 ALLOCATE(his4d(gs(1),gs(2),gs(3)+1,gs(4))); his4d=dble(0)
 CALL ncwrite4d(ncidOut,'AKt',[1,1,1,1],[gs(1),gs(2),gs(3)+1,gs(4)],&
 &his4d)
 CALL ncwrite4d(ncidOut,'AKv',[1,1,1,1],[gs(1),gs(2),gs(3)+1,gs(4)],&
 &his4d)
 DEALLOCATE(his4d)

 write(*,*) 'write ubar'
 ALLOCATE(his3d(gs(1)-1,gs(2),gs(4))); his3d=dble(0)
 CALL ncwrite3d(ncidOut,'ubar',[1,1,1],[gs(1)-1,gs(2),gs(4)],his3d)
 DEALLOCATE(his3d)

 write(*,*) 'write vbar'
 ALLOCATE(his3d(gs(1),gs(2)-1,gs(4))); his3d=dble(0)
 CALL ncwrite3d(ncidOut,'vbar',[1,1,1],[gs(1),gs(2)-1,gs(4)],his3d)
 DEALLOCATE(his3d)

 !Fields
 write(*,*) 'allocate other fields'
 ALLOCATE(temp(gs(1),gs(2),gs(3),gs(4))); temp=dble(0)
 ALLOCATE(salt(gs(1),gs(2),gs(3),gs(4))); salt=dble(0)
 ALLOCATE(u(gs(1)-1,gs(2),gs(3),gs(4))); u=dble(0)
 ALLOCATE(v(gs(1),gs(2)-1,gs(3),gs(4))); v=dble(0)
 ALLOCATE(zeta(gs(1),gs(2),gs(4))); zeta=dble(0)

 !---------------------------------------------------------------------------
 !Demean all
 
#ifdef demean_all
 DO i0=MINVAL(obs(:,5)),MAXVAL(obs(:,5))
  IF(.NOT.ANY(obs(:,5).EQ.i0)) CYCLE
  valm=demean(valm,obs(:,5).EQ.i0)
 END DO
#endif

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
     u(ix(j1,j0),iy(j1,j0),gs(3),it(i1,1))=&
     &u(ix(j1,j0),iy(j1,j0),gs(3),it(i1,1))+&
     &wt(i1,1)*wxy(j1,j0)*valm(iObs(j0))
    END DO !j1                                                                 
    END DO !i1                                                                 
    flag_todo(iObs(j0))=.FALSE.
   END IF
  END DO !j0                                                                   

 END DO !do while   

 !clean
 DEALLOCATE(iobs,ix,iy,it,wxy,wt)
 WRITE(*,*) 'min/max u_surface:',MINVAL(u(:,:,gs(3),:)),&
 &MAXVAL(u(:,:,gs(3),:))

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
     v(ix(j1,j0),iy(j1,j0),gs(3),it(i1,1))=&
     &v(ix(j1,j0),iy(j1,j0),gs(3),it(i1,1))+&
     &wt(i1,1)*wxy(j1,j0)*valm(iObs(j0))
    END DO !j1                                                                 
    END DO !i1                                                                 
    flag_todo(iObs(j0))=.FALSE.
   END IF
  END DO !j0                                                                   

 END DO !do while   

 !clean
 DEALLOCATE(iobs,ix,iy,it,wxy,wt)
 WRITE(*,*) 'min/max v_surface:',MINVAL(v(:,:,gs(3),:)),&
 &MAXVAL(v(:,:,gs(3),:))

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
     temp(ix(j1,j0),iy(j1,j0),gs(3),it(i1,1))=&
     &temp(ix(j1,j0),iy(j1,j0),gs(3),it(i1,1))+&
     &wt(i1,1)*wxy(j1,j0)*valm(iObs(j0))
    END DO !j1                                                                 
    END DO !i1                                                                 
    flag_todo(iObs(j0))=.FALSE.
   END IF
  END DO !j0                                                                   

 END DO !do while   

 !clean
 DEALLOCATE(iobs,ix,iy,it,wxy,wt)
 WRITE(*,*) 'min/max temp_surface:',MINVAL(temp(:,:,gs(3),:)),&
 &MAXVAL(temp(:,:,gs(3),:))

 END IF !type 4
 !---------------------------------------------------------------------------
 ! Process type 5: daily-averaged sea-surface temperature

 typeo=[5,5] !observation type
 dt=dble(24*3600) !Averaging period
 dl=[.125,.125,0.0] !Half-width in zonal and meridional direction (degrees)

 !Select entries to process                                                    
 flag_todo=INT(obs(:,5)).Ge.typeo(1).AND.INT(obs(:,5)).LE.typeo(2)
 nSelect=COUNT(flag_todo)
 WRITE(*,*) 'type 5:',nSelect,' samples'
 IF(ANY(flag_todo)) THEN

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
    temp(ix(i1,i2),iy(i1,i2),gs(3),it(i4,1))=&
    &temp(ix(i1,i2),iy(i1,i2),gs(3),it(i4,1))+&
    &wxy(i1,i2)*wt(i4,1)*valm(j0) 
   END DO !i1                                                                  
   END DO !i2                                                                  
   END DO !i4     

   !Deallocate location                                                       
   DEALLOCATE(wxy,ix,iy)

   !Mark as complete                                                           
   flag_todo(j0)=.FALSE.

  END DO !j0                                                                   
 END DO !do while 
            
 !clean
 DEALLOCATE(it,wt)
 WRITE(*,*) 'min/max temp_surface:',MINVAL(temp(:,:,gs(3),:)),&
 &MAXVAL(temp(:,:,gs(3),:))

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
 his3d=ncread3d(ncidHis,'zeta',ncstart(1:3),nccount(1:3))

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

  !sample
  DO j0=1,size(iobs)
   IF(.NOT.selection(j0)) CYCLE
  DO i0=1,size(wxy,1)
  DO i3=1,size(wz,1)
  DO i4=1,size(wt,1)
   IF(ix(i0,j0).EQ.i1.AND.iy(i0,j0).EQ.i2) THEN
    temp(i1,i2,iz(i3,j0),it(i4,j0))=&
    &temp(i1,i2,iz(i3,j0),it(i4,j0))+&
    wxy(i0,j0)*wz(i3,j0)*wt(i4,j0)*valm(iObs(j0)) 
   END IF
  END DO !i4
  END DO !i3
  END DO !i0
  END DO !j0

 END DO !i2
 END DO !i1

 DEALLOCATE(iobs,ix,iy,iz,it,wxy,wt,wz,selection)
 WRITE(*,*) 'min/max subsurface temp:',&
 &MINVAL(temp(:,:,1:gs(3)-1,:)),MAXVAL(temp(:,:,1:gs(3)-1,:)) 

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
 his3d=ncread3d(ncidHis,'zeta',ncstart(1:3),nccount(1:3))

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
 
 !Sample salt
 DO i1=MINVAL(ix),MAXVAL(ix)
 DO i2=MINVAL(iy),MAXVAL(iy)
  selection=ANY(ix.EQ.i1.AND.iy.EQ.i2,1)
  IF(.NOT.ANY(selection)) CYCLE 

  !sample
  DO j0=1,size(iobs)
   IF(.NOT.selection(j0)) CYCLE
  DO i0=1,size(wxy,1)
  DO i3=1,size(wz,1)
  DO i4=1,size(wt,1)
   IF(ix(i0,j0).EQ.i1.AND.iy(i0,j0).EQ.i2) THEN
    salt(i1,i2,iz(i3,j0),it(i4,j0))=&
    &salt(i1,i2,iz(i3,j0),it(i4,j0))+&
    valm(iobs(j0))*wxy(i0,j0)*wz(i3,j0)*wt(i4,j0)
   END IF
  END DO !i4
  END DO !i3
  END DO !i0
  END DO !j0

 END DO !i2
 END DO !i1

 !clean
 DEALLOCATE(iobs,ix,iy,iz,it,wxy,wt,wz,selection)
 WRITE(*,*) 'min/max subsurface salt:',&
 &MINVAL(salt(:,:,1:gs(3)-1,:)),MAXVAL(salt(:,:,1:gs(3)-1,:)) 

 END IF !type 7

!--------------------------------------------------------------------------    
! Process type 8: instant radial velocity                              

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
 ALLOCATE(wt(gs(4), 1)); wt=dble(0)

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
     u(ix(j1,j0),iy(j1,j0),gs(3),it(i1,1))=&
     &u(ix(j1,j0),iy(j1,j0),gs(3),it(i1,1))+&
     &SIN(obs(iObs(j0),7))*wt(i1,1)*wxy(j1,j0)*&
     &valm(iObs(j0))
    END DO !j1                                                                 
    END DO !i1                                                                 
    flag_todo(iObs(j0))=.FALSE.
   END IF
  END DO !j0                                                                   
 END DO !do while          

 !Allocate for V                                                            
 DEALLOCATE(iobs,ix,iy,it,wxy,wt)
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
     v(ix(j1,j0),iy(j1,j0),gs(3),it(i1,1))=&
     &v(ix(j1,j0),iy(j1,j0),gs(3),it(i1,1))+&
     &COS(obs(iObs(j0),7))*wt(i1,1)*wxy(j1,j0)*&
     &valm(iObs(j0))
    END DO !j1                                                                 
    END DO !i1                                                                 
    flag_todo(iObs(j0))=.FALSE.
   END IF
  END DO !j0                                                                   
 END DO !do while   

 DEALLOCATE(iobs,ix,iy,it,wxy,wt)
 WRITE(*,*) 'min/max u:',MINVAL(u(:,:,size(u,3),:)),&
 &MAXVAL(u(:,:,size(u,3),:))
 WRITE(*,*) 'min/max v:',MINVAL(v(:,:,size(v,3),:)),&
 &MAXVAL(v(:,:,size(v,3),:))

 END IF !type 8    
 !--------------------------------------------------------------------------
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
  his3d=ncread3d(ncidHis,'zeta',ncstart(1:3),nccount(1:3))

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
     salt(ix(i1,i2),iy(i1,i2),iz(1,i3),it(i4,1))=&
     &salt(ix(i1,i2),iy(i1,i2),iz(1,i3),it(i4,1))+&
     &valm(j0)*wz(1,i3)*wt(i4,1)*wxy(i1,i2)
    END DO
    DEALLOCATE(iz,wz)

   END DO !i1                                                               
   END DO !i2                                                                       
   END DO !i4                                                                                  
   DEALLOCATE(ix,iy,wxy)

   !Mark as complete                                                                       
   flag_todo(j0)=.FALSE.

  END DO !j0
  DEALLOCATE(his3d)

 END DO !while                                                                                              
 DEALLOCATE(it,wt)
 flag_todo=INT(obs(:,5)).Ge.typeo(1).AND.INT(obs(:,5)).LE.typeo(2)
 WRITE(*,*) 'min/max salt:',MINVAL(salt(:,:,gs(3),1)),&
 &MAXVAL(salt(:,:,gs(3),1))

 END IF !type 100-199        

 !---------------------------------------------------------------------------
 ! Process type 1000-9999: detided sea-surface height

 typeo=[1000,9999] !observation type
 dt=dble(24*3600) !time length filter 
 WRITE(*,*) 'type 1000-9999:'

 DO track=typeo(1),typeo(2)
 !Select entries to process                                                    
 flag_todo=INT(obs(:,5)).EQ.track
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

 !Demean
 valm=demean(valm,obs(:,5).EQ.track)

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
     zeta(ix(j1,j0),iy(j1,j0),it(i1,1))=&
     &zeta(ix(j1,j0),iy(j1,j0),it(i1,1))+&
     &wt(i1,1)*wxy(j1,j0)*valm(iObs(j0))
    END DO !j1                                                                 
    END DO !i1                                                                 
    flag_todo(iObs(j0))=.FALSE.
   END IF
  END DO !j0                                                                   

 END DO !do while   

 !clean
 DEALLOCATE(iobs,ix,iy,it,wxy,wt)
 END IF !type.EQ.track
 END DO !track
 
 WRITE(*,*) 'min/max zeta:',MINVAL(zeta),&
 &MAXVAL(zeta)

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

  !Seperate in mean and deviations from mean
  WRITE(*,*) 'mean track:',valm(iObs(1))
  DO i0=2,nSelect
   DO j0=1,nSelect
    tmp1(j0)=tmp1(j0)-valm(iObs(i0))/dble(nSelect)
   END DO
   tmp1(i0)=tmp1(i0)+valm(iObs(i0)) 
  END DO
  DO j0=1,nSelect
   tmp1(j0)=tmp1(j0)+valm(iObs(1))/dble(nSelect)
  END DO

  DO WHILE(ANY(flag_todo))
  !time                                                                       
  timeo=MINVAL(obs(:,4),1,flag_todo)
  CALL sample_index_delta(timem,obs(i0,4),it(:,1),wt(:,1))

  !Sample                                                                      
  DO j0=1,nSelect
   IF(obs(iobs(j0),4).NE.timeo) CYCLE
   DO i1=MINVAL(it,it.GT.0),MAXVAL(it,it.GT.0)
   DO j1=1,size(ix,1)
    temp(ix(j1,j0),iy(j1,j0),gs(3),it(i1,1))=&
    &temp(ix(j1,j0),iy(j1,j0),gs(3),it(i1,1))+&
    tmp1(j0)*wt(i1,1)*wxy(j1,j0)
   END DO !j1                                                                 
   END DO !i1          

   flag_todo(iObs(j0))=.FALSE.
  END DO !j0                                                                  
  END DO !do while                                                              
  WRITE(*,*) 'min/max SST:',MINVAL(temp(:,:,gs(3),:)),&
  MAXVAL(temp(:,:,gs(3),:))
  DEALLOCATE(iobs,ix,iy,it,wxy,wt,his4d,tmp1)

 END DO !type 10000-19999                                                      
 !--------------------------------------------------------------------------
 !Write fields

 write(*,*) 'writing fields to netcdf file with size',gs
 CALL ncwrite4d(ncidOut,'temp',[1,1,1,1],[gs(1),gs(2),gs(3),gs(4)],&
 &temp)
 CALL ncwrite4d(ncidOut,'salt',[1,1,1,1],[gs(1),gs(2),gs(3),gs(4)],&
 &salt)
 CALL ncwrite4d(ncidOut,'u',[1,1,1,1],[gs(1)-1,gs(2),gs(3),gs(4)],&
 &u)
 CALL ncwrite4d(ncidOut,'v',[1,1,1,1],[gs(1),gs(2)-1,gs(3),gs(4)],&
 &v)
 CALL ncwrite3d(ncidOut,'zeta',[1,1,1],[gs(1),gs(2),gs(4)],&
 &zeta)
 CALL ncwrite1d(ncidOut,'ocean_time',[1],[gs(4)],timem)

 !---------------------------------------------------------------------------

 !Close streams
 status=nf_close(ncidOut)
 status=nf_close(ncidHis)

 !-------------------------------------------------------------------------
 WRITE(*,*) 'DONE'
END PROGRAM ad_sample

!-----------------------------------------------------------------------------

FUNCTION demean(val_in,mask) RESULT(val_out)
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
