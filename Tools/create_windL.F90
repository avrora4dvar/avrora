PROGRAM	create_wind
 USE mod_netcdf
 USE mod_bhm
 IMPLICIT NONE
!--------------------------------------------------------------------------
!Declare
 
 !Input
 CHARACTER(len=1024)::file_in(2),file_out(3),file_svd
 INTEGER::nEOF,dateref(6)
 REAL(8):: time_bnd(2)
 INTEGER::ncid(5),status,var_id(20),dim_id(3)
 LOGICAL::flag_exist
 REAL::eps

 !reading wind fields
 REAL,allocatable::uIn(:,:,:),vIn(:,:,:)
 INTEGER::fsize(3),itime_bnd(2),dimlen(8)
 REAL(8),allocatable::time_in(:)

 !EOF fields
 REAL,allocatable::uMean(:,:),vMean(:,:)
 REAL,allocatable::uEOF(:,:,:),vEOF(:,:,:),lam(:)
 COMPLEX,allocatable::uvEOF(:,:,:),uvMean(:,:)
 INTEGER::i0,i1,i2,i3 
 COMPLEX,allocatable::true_cL(:,:)
 REAL::bias(2,2)

 !small-scale fields
 REAL,allocatable::corr_cS(:),sig_cS(:)

 !Generating output
 INTEGER::step,stepStart

 !Output winds
 REAL,allocatable::vOut(:,:,:),uOut(:,:,:)
 
!---------------------------------------------------------------------------
!Read input
 
 READ(*,*) !Input file with input field
 READ(*,'(A)') file_in(1)
 READ(*,*) !File with EOFs
 READ(*,'(A)') file_svd
 READ(*,*) !Number of EOFs to be used
 READ(*,*) nEOF
 READ(*,*) !start time, end time (days since dateref)
 READ(*,*) time_bnd(1), time_bnd(2)
 READ(*,*) !date ref (yyyy mm dd HH MM SS)
 READ(*,*) dateref
 READ(*,*) !output file (file name will be appended with _u.nc, _v.nc)
 READ(*,'(A)') file_out(1)

 WRITE(*,*) 'CREATE_WIND input:'
 WRITE(*,*) 'Input file with u-wind field'
 WRITE(*,'(A)') TRIM(file_in(1))
 WRITE(*,*) 'Input file with v-wind field'
 WRITE(*,'(A)') TRIM(file_in(1))
 WRITE(*,*) 'File with EOFs'
 WRITE(*,'(A)') TRIM(file_svd)
 WRITE(*,*) 'Number of EOFs to be used'
 WRITE(*,*) nEOF
 WRITE(*,*) 'start time, end time (days since dateref)'
 WRITE(*,*) time_bnd(1), time_bnd(2)
 WRITE(*,*) 'date ref (yyyy mm dd HH MM SS)'
 WRITE(*,*) dateref
 WRITE(*,*) 'output file (file name will be appended with _u.nc, _v.nc)'
 WRITE(*,*) TRIM(file_out(1))

 eps=1e-6

 INQUIRE(file=TRIM(file_svd),exist=flag_exist)
 
!--------------------------------------------------------------------------
!Load u and v date from input files

 WRITE(*,*) 'Reading u,v from input files'

 !Open netcdf streams
 status=nf_open(TRIM(file_in(1)),nf_nowrite,ncid(1))

 !Read size input fields
 CALL ncsize(ncid(1),'Uwind',dimlen)
 fsize=dimlen(1:3)

 !Read time
 ALLOCATE(time_in(fsize(3)))
 time_in=ncread1d(ncid(1),'wind_time',[1],[fsize(3)])
 itime_bnd(1)=MAXLOC(time_in,1,time_in.LE.time_bnd(1))
 itime_bnd(2)=MINLOC(time_in,1,time_in.GE.time_bnd(2))
 fsize(3)=itime_bnd(2)-itime_bnd(1)+1
 WRITE(*,*) 'Reading time steps ',itime_bnd(1),'-',itime_bnd(2)

 !Read wind fields
 WRITE(*,*) 'Read wind fields' 
 ALLOCATE(uIn(fsize(1),fsize(2),fsize(3)))
 ALLOCATE(vIn(fsize(1),fsize(2),fsize(3)))
 uIn=ncread3d(ncid(1),'Uwind',[1,1,itime_bnd(1)],fsize)
 vIn=ncread3d(ncid(1),'Vwind',[1,1,itime_bnd(1)],fsize)

 !Close netcdf stream
 status=nf_close(ncid(1))

 write(*,*) 'min/max u:',minval(uIn),' ',maxval(uIn)
 write(*,*) 'min/max v:',minval(vIn),' ',maxval(vIn)

!--------------------------------------------------------------------------
!Read EOFs

 IF(flag_exist) THEN
 WRITE(*,*) 'Reading EOFs from ',TRIM(file_svd)

 ALLOCATE(uEOF(fsize(1),fsize(2),nEOF))
 ALLOCATE(vEOF(fsize(1),fsize(2),nEOF))
 ALLOCATE(uvEOF(fsize(1),fsize(2),nEOF))
 ALLOCATE(uMean(fsize(1),fsize(2)))
 ALLOCATE(vMean(fsize(1),fsize(2)))
 ALLOCATE(uvMean(fsize(1),fsize(2)))
 ALLOCATE(lam(nEOF))

 status=nf_open(TRIM(file_svd),nf_nowrite,ncid(5))
 uEOF=ncread3d(ncid(5),'Uwind',[1,1,1],[fsize(1),fsize(2),nEOF])
 vEOF=ncread3d(ncid(5),'Vwind',[1,1,1],[fsize(1),fsize(2),nEOF])
 uMean=ncread2d(ncid(5),'Umean',[1,1],[fsize(1),fsize(2)])
 vMean=ncread2d(ncid(5),'Vmean',[1,1],[fsize(1),fsize(2)])
 lam=ncread1d(ncid(5),'lam',[1],[nEOF])
 status=nf_close(ncid(5))

 uvEOF=CMPLX(uEOF,vEOF)
 uvMean=CMPLX(uMean,vMean)
 uvMean=(0.0,0.0)*uvMean !TMP
 DEALLOCATE(uEOF,vEOF)

 WRITE(*,*) 'lam: ',lam
 END IF

!--------------------------------------------------------------------------
!Create output


 IF(flag_exist) THEN
  WRITE(*,*) 'Creating large-scale wind field'
  
  !Allocate BHM fields
  ALLOCATE(bhm_mask(fsize(1),fsize(2))); bhm_mask=.FALSE.
  ALLOCATE(bhm_uvTrue(fsize(1),fsize(2),fsize(3))); bhm_uvTrue=(0.0,0.0)
  ALLOCATE(bhm_uvL(fsize(1),fsize(2),fsize(3))); bhm_uvL=(0.0,0.0)
  ALLOCATE(bhm_uvS(fsize(1),fsize(2),fsize(3))); bhm_uvS=(0.0,0.0)

  !Subsample the input wind field to create 'observations'  
  step=115; stepStart=CEILING(step/2.)
  bhm_mask(stepStart:fsize(1):step,stepStart:fsize(2):step)=.TRUE.
  WHERE( ALL(ABS(uvEOF).LT.eps,3) )
   bhm_mask=.FALSE.
  END WHERE 
  WHERE( SPREAD(bhm_mask,3,fsize(3)) )
    bhm_uvTrue=CMPLX(uIn,vIn)
  END WHERE
  WRITE(*,*) 'number of samples:',COUNT(bhm_mask)

  !Calculate projection on EOF
  ALLOCATE(true_cL(nEOF,fsize(3)))
  DO i2=1,fsize(3)
  DO i1=1,nEOF
   true_cL(i1,i2)=SUM( CONJG(uvEOF(:,:,i1))*&
   &( cmplx(uIn(:,:,i2),vIn(:,:,i2))-uvMean ))
  END DO
  END DO
  WRITE(*,*) 'true_cL:',true_cL(:,1)

  !Add bias
  CALL bhm_draw_bias((0.0,0.0),.4,[0.0,0.0])
  DO i3=1,fsize(3)
   bhm_uvTrue(stepStart:fsize(1):step,stepStart:fsize(2):step,&
   &i3)=bhm_uvTrue(stepStart:fsize(1):step,stepStart:fsize(2):step,i3)-&
   &bhm_bias(i3)
  END DO
  WRITE(*,*) 'bias:',bhm_bias(1)

  !Draw coefficients
  ALLOCATE(bhm_var_cL(nEOF))
  bhm_var_cL=REAL(lam**2)
  bhm_var_uv=[3.4,3.4]
  CALL bhm_draw_cEOF(uvMean,uvEOF)
  WRITE(*,*) 'bhm_cL:',bhm_cL(:,1)

  !Convert back to normal mask
  bhm_mask=.TRUE.
  WHERE( ALL(ABS(uvEOF).LT.eps,3) )
   bhm_mask=.FALSE.
  END WHERE 

  !Construct large scale wind field
  CALL bhm_compose_uvL(uvMean,uvEOF)
    
  WRITE(*,*) 'min/max uL:', minval(REAL(bhm_uvL)),maxval(REAL(bhm_uvL))
 END IF !switch if

!--------------------------------------------------------------------
!Add small-scale wind field

 WRITE(*,*) 'construct small-scale wind field'

 write(*,*) 'min/max uvS:',minval(real(bhm_uvS)),maxval(real(bhm_uvS)) 

!---------------------------------------------------------------------

 WRITE(*,*) 'Creating new field'

 ALLOCATE(uOut(fsize(1),fsize(2),fsize(3))); uOut=uIn
 ALLOCATE(vOut(fsize(1),fsize(2),fsize(3))); vOut=vIn; 

 WHERE( SPREAD(bhm_mask,3,fsize(3)) )
  uOut=REAL(bhm_uvS)+REAL(bhm_uvL)
  vOut=IMAG(bhm_uvS)+IMAG(bhm_uvL)
 END WHERE

  WRITE(*,*) 'min/max uOut:',minval(uOut),maxval(uOut)
  WRITE(*,*) 'min/max vOut:',minval(vOut),maxval(vOut)
!---------------------------------------------------------------------

 IF(.NOT.flag_exist) THEN
  WRITE(*,*) 'Cannot find ',TRIM(file_svd)
  STOP
 END IF

!---------------------------------------------------------------------------- 
! Write output files

 !Remove large discrepencies input and output
 WHERE( SQRT( (uOut-uIn)**2+(vOut-vIn)**2 ).GT.10.0  )
   uOut=uIn+(uOut-uIn)*10.0/sqrt( (uOut-uIn)**2+(vOut-vIn)**2 )
   vOut=vIn+(vOut-vIn)*10.0/sqrt( (uOut-uIn)**2+(vOut-vIn)**2 )
 END WHERE

 WRITE(*,*) 'Writing output'
 !Name files
 WRITE(*,*) 'Create ',TRIM(file_out(1))

 !Open netcdf streams
 status=nf_create(TRIM(file_out(1)),nf_classic_model,ncid(3))

 !Define dimensions
 status=nf_def_dim(ncid(3),'xi_rho',fsize(1),dim_id(1))
 status=nf_def_dim(ncid(3),'eta_rho',fsize(2),dim_id(2))
 status=nf_def_dim(ncid(3),'time',nf_unlimited,dim_id(3))

 !Define variables
 status=nf_open(TRIM(file_in(1)),nf_nowrite,ncid(1))
 status=nf_inq_varid(ncid(1),'wind_time',var_id(3))
 status=nf_inq_varid(ncid(1),'Uwind',var_id(1))
 status=nf_inq_varid(ncid(1),'Vwind',var_id(2))
 status=nf_def_var(ncid(3),'wind_time',nf_double,1,dim_id(3),var_id(5))
 status=nf_def_var(ncid(3),'Uwind',nf_float,3,dim_id,var_id(6))
 status=nf_def_var(ncid(3),'Vwind',nf_float,3,dim_id,var_id(8))

 !Define attributes
 status=nf_copy_att(ncid(1),var_id(3),ncid(3),'units',&
 &ncid(3),var_id(5))
 status=nf_copy_att(ncid(1),var_id(1),'units',&
 &ncid(3),var_id(6))
 status=nf_copy_att(ncid(1),var_id(1),'time',&
 &ncid(3),var_id(6))
 status=nf_copy_att(ncid(1),var_id(2),'units',&
 &ncid(3),var_id(8))
 status=nf_copy_att(ncid(1),var_id(2),'time',&
 &ncid(3),var_id(8))

 !End define mode
 status=nf_close(ncid(3))
 status=nf_close(ncid(1))
 status=nf_open(TRIM(file_out(1)),nf_write,ncid(3))

 WRITE(*,*) 'Write u,v'
 !Write output to file_out
 CALL ncwrite3d(ncid(3),'Uwind',[1,1,1],fsize,dble(uOut))
 CALL ncwrite3d(ncid(3),'Vwind',[1,1,1],fsize,dble(vOut))
 CALL ncwrite1d(ncid(3),'wind_time',[1],fsize(3),&
 &time_in(itime_bnd(1):itime_bnd(2)) )

 !close stream
 status=nf_close(ncid(3))

 WRITE(*,*) 'create_wind done'

END PROGRAM create_wind

!---------------------------------------------------------------------------


