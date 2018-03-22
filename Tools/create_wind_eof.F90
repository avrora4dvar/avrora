PROGRAM	create_wind
 USE mod_netcdf
 USE mod_BHM
 IMPLICIT NONE
!--------------------------------------------------------------------------
!Declare
 
 !Input
 CHARACTER::mode
 CHARACTER(len=1024)::file_in(2),file_out(3),file_svd
 INTEGER::nEOF,dateref(6),n_sample
 REAL(8):: time_bnd(2),K
 INTEGER::ncid(5),status,var_id(20),dim_id(3)
 LOGICAL::flag_exist
 REAL::eps

 !reading wind fields
 REAL,allocatable::uIn(:,:,:),vIn(:,:,:)
 COMPLEX,allocatable::uvIn(:,:,:),uvTmp(:,:,:)
 INTEGER::fsize(3),itime_bnd(2),dimlen(8)
 REAL(8),allocatable::time_in1(:),time_in2(:),time_out(:)
 REAL(8),allocatable::flag_perturb(:)

 !EOF fields
 REAL,allocatable::uMean(:,:),vMean(:,:)
 REAL,allocatable::uEOF(:,:,:),vEOF(:,:,:),sig_EOF(:)
 COMPLEX,allocatable::uvEOF(:,:,:),uvMean(:,:)

 !large-scale fields
 REAL,allocatable::rand(:),rand2(:,:),pdf(:,:),corr_cL(:)
 INTEGER,allocatable::seq(:)
 INTEGER::i0,i1,i2,i3,j0,j1,j2,j3
 COMPLEX,allocatable::true_cL(:,:)
 REAL::bias(2,2),cdf(2)

 !small-scale fields
 REAL,allocatable::corr_cS(:),sig_cS(:)
 REAL::uvS_shift(3),var_cS

 !Generating output
 INTEGER::step,stepStart

 !Output winds
 COMPLEX,allocatable::uvOut(:,:,:)

 REAL(8),allocatable::w(:) 
!---------------------------------------------------------------------------
!Read input
 
 READ(*,*) !Input files with wind fields
 READ(*,'(A)') file_in(1)
 READ(*,'(A)') file_in(2)
 READ(*,*) !File with EOFs
 READ(*,'(A)') file_svd
 READ(*,*) !Number of samples for generating large scale field
 READ(*,*) n_sample
 READ(*,*) !start time, end time (days since dateref)
 READ(*,*) time_bnd(1), time_bnd(2)
 READ(*,*) !date ref (yyyy mm dd HH MM SS)
 READ(*,*) dateref
 READ(*,*) !output file (file name will be appended with _u.nc, _v.nc)
 READ(*,'(A)') file_out(1)
 READ(*,*) !output file for coefficients
 READ(*,'(A)') file_out(2)
 READ(*,*) !Seed number (integer>0)
 READ(*,*) bhm_seed
 mode='R'

 WRITE(*,*) 'CREATE_WIND input:'
 WRITE(*,*) 'Input file with winds that needs to be perturbed'
 WRITE(*,'(A)') TRIM(file_in(1))
 WRITE(*,*) 'Input file with winds that are not perturbed '
 WRITE(*,'(A)') TRIM(file_in(2))
 WRITE(*,*) 'File with EOFs'
 WRITE(*,'(A)') TRIM(file_svd)
 WRITE(*,*) 'start time, end time (days since dateref)'
 WRITE(*,*) time_bnd(1), time_bnd(2)
 WRITE(*,*) 'date ref (yyyy mm dd HH MM SS)'
 WRITE(*,*) dateref
 WRITE(*,*) 'output file (file name will be appended with _u.nc, _v.nc)'
 WRITE(*,*) TRIM(file_out(1))
 WRITE(*,*) 'seed number'
 WRITE(*,*) bhm_seed
 WRITE(*,*) 'mode'
 WRITE(*,*) mode

 eps=1e-6

 INQUIRE(file=TRIM(file_svd),exist=flag_exist)
 
!--------------------------------------------------------------------------
!Load u and v date from first wind input file

 WRITE(*,*) 'Reading u,v from 1st input file'

 !Open netcdf streams
 status=nf_open(TRIM(file_in(1)),nf_nowrite,ncid(1))

 !Read size input fields
 CALL ncsize(ncid(1),'Uwind',dimlen)
 fsize=dimlen(1:3)

 !Read time
 ALLOCATE(time_in1(fsize(3)))
 time_in1=ncread1d(ncid(1),'wind_time',[1],[fsize(3)])
 itime_bnd(1)=MAXLOC(time_in1,1,time_in1.LE.time_bnd(1))
 itime_bnd(2)=MINLOC(time_in1,1,time_in1.GE.time_bnd(2))
 fsize(3)=itime_bnd(2)-itime_bnd(1)+1
 WRITE(*,*) 'Reading time steps ',itime_bnd(1),'-',itime_bnd(2)

 !Output times
 ALLOCATE(time_out(itime_bnd(2)-itime_bnd(1)+1))
 time_out=time_in1(itime_bnd(1):itime_bnd(2))
 ALLOCATE(flag_perturb(size(time_out))); flag_perturb=dble(1)


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

 !Copy to complex array
 ALLOCATE(uvIn(fsize(1),fsize(2),fsize(3)))
 uvIn=cmplx(uIn,vIn)
 DEALLOCATE(uIn,vIn,time_in1)

!--------------------------------------------------------------------------
!Replace wind with files from second input

IF(LEN(TRIM(file_in(2))).GT.0) THEN
 
 WRITE(*,*) 'Reading u,v from 2nd input file'

 !Open netcdf streams
 status=nf_open(TRIM(file_in(2)),nf_nowrite,ncid(1))

 !Read size input fields
 CALL ncsize(ncid(1),'Uwind',dimlen)
 fsize=dimlen(1:3)

 !Read time
 ALLOCATE(time_in2(fsize(3)))
 time_in2=ncread1d(ncid(1),'wind_time',[1],[fsize(3)])
 itime_bnd(1)=MINLOC(time_in2,1,time_in2.GE.MINVAL(time_out))
 itime_bnd(2)=MAXLOC(time_in2,1,time_in2.LE.MAXVAL(time_out))
 fsize(3)=itime_bnd(2)-itime_bnd(1)+1
 WRITE(*,*) 'Reading time steps ',itime_bnd(1),'-',itime_bnd(2)

 !Read wind fields
 WRITE(*,*) 'Read wind fields from 2nd file'
 ALLOCATE(uIn(fsize(1),fsize(2),fsize(3)))
 ALLOCATE(vIn(fsize(1),fsize(2),fsize(3)))
 uIn=ncread3d(ncid(1),'Uwind',[1,1,itime_bnd(1)],fsize)
 vIn=ncread3d(ncid(1),'Vwind',[1,1,itime_bnd(1)],fsize)

 !Close netcdf stream
 status=nf_close(ncid(1))

 write(*,*) 'min/max u:',minval(uIn),' ',maxval(uIn)
 write(*,*) 'min/max v:',minval(vIn),' ',maxval(vIn)

 !Replace wind fields
 DO i3=1,size(time_out)
 DO i0=1,fsize(3)
  IF(time_out(i3).NE.time_in2(i0-1+itime_bnd(1))) CYCLE
  WRITE(*,*) 'Replacing time ',time_out(i3)
  flag_perturb(i3)=dble(0)
  uvIn(:,:,i3)=cmplx(uIn(:,:,i0),vIn(:,:,i0))
 END DO
 END DO
 DEALLOCATE(uIn,vIn,time_in2)
 fsize(3)=size(time_out)

END IF
!--------------------------------------------------------------------------
!Read EOFs

 WRITE(*,*) 'Reading EOFs from ',TRIM(file_svd)

 !Read size fields
 status=nf_open(TRIM(file_svd),nf_nowrite,ncid(5))
 CALL ncsize(ncid(5),'Uwind',dimlen)
 nEOF=dimlen(3)

 !Allocate EOF arrays
 ALLOCATE(uEOF(fsize(1),fsize(2),nEOF))
 ALLOCATE(vEOF(fsize(1),fsize(2),nEOF))
 ALLOCATE(uvEOF(fsize(1),fsize(2),nEOF))
 ALLOCATE(uMean(fsize(1),fsize(2)))
 ALLOCATE(vMean(fsize(1),fsize(2)))
 ALLOCATE(uvMean(fsize(1),fsize(2)))
 ALLOCATE(sig_EOF(nEOF))

 !Read EOFs and mean
 uEOF=ncread3d(ncid(5),'Uwind',[1,1,1],[fsize(1),fsize(2),nEOF])
 vEOF=ncread3d(ncid(5),'Vwind',[1,1,1],[fsize(1),fsize(2),nEOF])
 uMean=ncread2d(ncid(5),'Umean',[1,1],[fsize(1),fsize(2)])
 vMean=ncread2d(ncid(5),'Vmean',[1,1],[fsize(1),fsize(2)])
 sig_EOF=ncread1d(ncid(5),'sig',[1],[nEOF])
 status=nf_close(ncid(5))

 !Combine into complex arrays
 uvEOF=CMPLX(uEOF,vEOF)
 uvMean=CMPLX(uMean,vMean)
 DEALLOCATE(uEOF,vEOF)

 WRITE(*,*) 'sig_EOF: ',sig_EOF

!------------------------------------------------------------------------
!Set BHM 

  WRITE(*,*) 'Allocating  fields'
 
  !Allocate BHM fields
  ALLOCATE(bhm_mask(fsize(1),fsize(2))); bhm_mask=.FALSE.
  ALLOCATE(bhm_uvTrue(fsize(1),fsize(2),fsize(3))); bhm_uvTrue=(0.0,0.0)
  ALLOCATE(bhm_uvL(fsize(1),fsize(2),fsize(3))); bhm_uvL=(0.0,0.0)
  ALLOCATE(bhm_uvS(fsize(1),fsize(2),fsize(3))); bhm_uvS=(0.0,0.0)
  
!--------------------------------------------------------------------
!Draw bias

  !Add bias (sig_bias=1.2)
  CALL bhm_draw_bias((0.0,0.0),.4,[0.0,0.0])
  WRITE(*,*) 'bias:',bhm_bias(1)

!--------------------------------------------------------------------
!Small-scale wind field

 IF(.true.) THEN
  WRITE(*,*) 'construct small-scale wind field'

  !mask
  bhm_mask=.TRUE. 

  !Order in DB2 hierarchy
  CALL bhm_cal_DB_order()

  CALL bhm_cal_DB_order();
  !Standard deviation different wavelet coefficients (upto norm. factor)
  ALLOCATE(sig_cS(9)); sig_cS=0.0
  !DO i1=1,6; sig_cS(i1)=SQRT(1.25**dble(i1)*2.75**dble(2*i1)); END DO
  !Spectrum P( k0/(2**i1) )=P0*(sig_cS(i1))**2
  DO i1=1,9
   sig_cS(i1)=exp(1.3*dble(i1-3))
   sig_cS(i1)=sig_cS(i1)*(.5-.5*&
   &(EXP(.4*(i1-4.0))-EXP(.4*(4.0-i1)))/&
   &(EXP(.4*(i1-4.0))+EXP(.4*(4.0-i1))) )
  END DO 
     
  !Correlation wavelet coefficients at consequetive times
  ALLOCATE(corr_cS(9))
  corr_cS=0.4 

  !Draw small scale wind field
  CALL bhm_draw_cS(corr_cS,sig_cS)

  !Create small-scale wind field
  CALL bhm_compose_uvS()
 
  !Normalize such that real(bhm_uvS) and imag(bhm_uvS) have a std of 1
  var_cS=.3 !2.2
  DO i3=1,size(bhm_uvS,3)
   bhm_uvS(:,:,i3)=bhm_uvS(:,:,i3)*&
   &SQRT( 2.0*dble(COUNT(bhm_mask))/SUM(ABS(bhm_uvS(:,:,i3))**2) )
  END DO
  WHERE(ABS(bhm_uvS).GT.3.0)
   bhm_uvS=bhm_uvS/ABS(bhm_uvS)*3.0
  END WHERE
  IF(mode.EQ.'C') THEN
   bhm_uvS=SQRT(.5*var_cS)*bhm_uvS; 
  ELSE
   bhm_uvS=SQRT(.5*var_cS)*bhm_uvS; 
  END IF
 
  !Shift field other you end up with stripes
  CALL bhm_random_uniform(1,uvS_shift)
  uvS_shift=uvS_shift*[size(bhm_uvS,1),size(bhm_uvS,2),0]
  write(*,*) 'uvS_shift:',INT(uvS_shift)
  bhm_uvS=CSHIFT(bhm_uvS,INT(uvS_shift(1)),1)
  bhm_uvS=CSHIFT(bhm_uvS,INT(uvS_shift(2)),2)

  write(*,*) 'min/max uvS:',minval(real(bhm_uvS)),maxval(real(bhm_uvS)) 
 END IF

!--------------------------------------------------------------------------
!Create large-scale wind field

 IF(.TRUE.) THEN
  WRITE(*,*) 'Creating large-scale wind field'

  bhm_mask=.TRUE.

  !Calculate projection on EOF
  ALLOCATE(true_cL(nEOF,fsize(3)))
  DO i2=1,fsize(3)
  DO i1=1,nEOF
   IF(mode.EQ.'C') THEN
    true_cL(i1,i2)=SUM( CONJG(uvEOF(:,:,i1))*&
    &(uvIn(:,:,i2)-uvMean),bhm_mask )
   ELSE
    true_cL(i1,i2)=SUM( REAL(uvEOF(:,:,i1))*REAL(uvIn(:,:,i2)-uvMean),&
    &bhm_mask )&
    &+SUM( IMAG(uvEOF(:,:,i1))*IMAG(uvIn(:,:,i2)-uvMean),&
    &bhm_mask )*(1.0,0.0)
   END IF
  END DO
  END DO
  WRITE(*,*) 'true_cL:',true_cL(:,1)

  !Draw EOF error coefficients
  IF(mode.EQ.'C') THEN
   ALLOCATE(bhm_var_cL(10))
   ALLOCATE(corr_cL(10))
   !Exp28f_2299-Exp28f_2311
   bhm_var_cL=[4.6e5,2.4e5,2.2e5,2.1e5,4.5e5,&
   &2.2e5,7.5e4,1.1e5,5.7e4,4.1e4]
   bhm_var_cL=.5*bhm_var_cL !Half variance in real, half in imag part
   corr_cL=0.0
  ELSE
   ALLOCATE(bhm_var_cL(10))
   ALLOCATE(corr_cL(10))
   !Exp28f_2314 and later
   bhm_var_cL=[2.1,1.2,.58,.44,.33,.39,.25,.25,.23,.22]*1e5
   corr_cL=.4
  END IF

  ALLOCATE(bhm_cL(nEOF,fsize(3))); bhm_cL=(0.0,0.0)
  ALLOCATE(rand(size(bhm_cL))); rand=0.0
  ALLOCATE(rand2(size(bhm_cL,1),size(bhm_cL,2)))

  !Real part 
  CALL bhm_random_gauss(1,rand)
  rand2=RESHAPE(rand,[size(rand2,1),size(rand2,2)])
  DO i1=1,nEOF
  DO i2=2,size(rand2,2)
   rand2(i1,i2)=corr_cL(i1)*rand2(i1,i2-1)&
   +SQRT(REAL(1.0-corr_cL(i1)**2))*rand2(i1,i2)
  END DO
   bhm_cL(i1,:)=bhm_cL(i1,:)+SQRT(bhm_var_cL(i1))*rand2(i1,:)*(1.0,0.0)
  END DO

  !Imaginary part
  IF(mode.EQ.'C') THEN
  CALL bhm_random_gauss(2,rand)
  rand2=RESHAPE(rand,[size(rand2,1),size(rand2,2)])
  DO i1=1,nEOF
  DO i2=2,size(rand2,2)
   rand2(i1,i2)=corr_cL(i1)*rand2(i1,i2-1)&
   &+SQRT(REAL(1.0-corr_cL(i1)**2))*rand2(i1,i2)
  END DO
   bhm_cL(i1,:)=bhm_cL(i1,:)+SQRT(bhm_var_cL(i1))*rand2(i1,:)*(0.0,1.0)
  END DO
  END IF

  DEALLOCATE(rand,rand2)
  WRITE(*,*) 'bhm_cL:',bhm_cL(:,1)
 
  !Construct large scale wind field
  CALL bhm_compose_uvL((0.0,0.0)*uvMean,uvEOF)

  WRITE(*,*) 'min/max uL:', minval(REAL(bhm_uvL)),maxval(REAL(bhm_uvL))
 END IF !switch if

!---------------------------------------------------------------------

 WRITE(*,*) 'Creating new field'

 ALLOCATE(uvOut(fsize(1),fsize(2),fsize(3)))
 DO i3=1,size(uvOut,3)
  uvOut(:,:,i3)=uvIn(:,:,i3)+&
  &flag_perturb(i3)*(bhm_uvL(:,:,i3)+bhm_uvS(:,:,i3)+&
  &bhm_bias(i3))
 END DO

 !Remove outliers
 WHERE( ABS(uvOut-uvIn).GT.10.0 )
  uvOut=uvIn+10.0*(uvOut-uvIn)/ABS(uvOut-uvIn)
 END WHERE

 WRITE(*,*) 'min/max uOut:',minval(REAL(uvOut)),maxval(REAL(uvOut))
 WRITE(*,*) 'min/max vOut:',minval(IMAG(uvOut)),maxval(IMAG(uvOut))

!--------------------------------------------------------------------------
! Write coefficients

 WRITE(*,*) 'Write coefficients to output'

 OPEN(unit=11,file=file_out(2))
 WRITE(11,'(A)') 'file_in_1'
 WRITE(11,'(A)') TRIM(file_in(1))
 WRITE(11,'(A)') 'file_in_2'
 WRITE(11,'(A)') TRIM(file_in(2))
 WRITE(11,'(A)') 'file_svd'
 WRITE(11,'(A)') TRIM(file_svd)
 WRITE(11,'(A)') 'time_out'
 WRITE(11,*) SHAPE(time_out)
 WRITE(11,*) time_out
 WRITE(11,'(A)') 'flag_perturb'
 WRITE(11,*) SHAPE(flag_perturb)
 WRITE(11,*) flag_perturb
 WRITE(11,'(A)') 'True cL'
 WRITE(11,*) SHAPE(true_cL)
 WRITE(11,*) true_cL
 WRITE(11,'(A)') 'bhm_cL'
 WRITE(11,*) SHAPE(bhm_cL)
 WRITE(11,*) bhm_cL
 WRITE(11,'(A)') 'uvS_shift'
 WRITE(11,*) uvS_shift
 !WRITE(11,'(A)') 'bhm_cS'
 !WRITE(11,*) SHAPE(bhm_cS)
 !WRITE(11,*) bhm_cS
 CLOSE(11)

!---------------------------------------------------------------------------- 
! Write output files

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
 CALL ncwrite3d(ncid(3),'Uwind',[1,1,1],fsize,dble(REAL(uvOut)))
 CALL ncwrite3d(ncid(3),'Vwind',[1,1,1],fsize,dble(IMAG(uvOut)))
 CALL ncwrite1d(ncid(3),'wind_time',[1],fsize(3),time_out)

 !close stream
 status=nf_close(ncid(3))

 WRITE(*,*) 'create_wind done'

!---------------------------------------------------------------------------

CONTAINS
FUNCTION wind_div(uv) 
 IMPLICIT NONE
 COMPLEX::uv(:,:,:)
 COMPLEX::wind_div(size(uv,1),size(uv,2),size(uv,3))
 REAL,dimension(size(uv,1),size(uv,2),size(uv,3))::phi,div,u,v,r
 INTEGER::i0,i1,i2,j1,j2,k1,k2
 REAL(8)::norm0
 REAL::omega

 IF(ANY(SHAPE(uv).NE.SHAPE(wind_div))) THEN
  WRITE(*,*) 'wind_div: shape input and output do not match'
 END IF

 !Interpolate u to u-points
 u=.5*CSHIFT(REAL(uv),-1,1)+.5*REAL(uv)
 !Interpolate v to v-points
 v=.5*CSHIFT(IMAG(uv),-1,2)+.5*IMAG(uv)
 v=IMAG(uv)
 !Calculate divergence in zeta-points
 div=(CSHIFT(u,1,1)-u)+(CSHIFT(v,1,2)-v)
 DO i0=1,size(div,3)
  div(:,:,i0)=div(:,:,i0)-SUM(div(:,:,i0))/dble(size(div,1)*size(div,2))
 END DO
 norm0=SUM(div*div)

 !Solve nabla(phi)=div(u,v) using Jacobi solver
 phi=0.0
 DO i0=1,1e4
  omega=.9

  phi=(1.0-omega)*phi-omega*.25*(div&
  &-CSHIFT(phi,-1,1)-CSHIFT(phi,1,1)&
  &-CSHIFT(phi,-1,2)-CSHIFT(phi,1,2))

  IF(MOD(i0,100).EQ.0) THEN
   !Check solution 
   r=div-(CSHIFT(phi,1,1)+CSHIFT(phi,-1,1)-2*phi)& 
   &-(CSHIFT(phi,1,2)+CSHIFT(phi,-1,2)-2*phi)
   WRITE(*,*) 'Rel. error at ',i0,SQRT(SUM(r*r)/norm0)
   IF(SQRT(SUM(r*r)/norm0)<.01) EXIT
  END IF
 END DO

 !Calculate divergence part u on u-points
 u=(phi-CSHIFT(phi,-1,1))
 !Calculate divergence part v on v-points
 v=(phi-CSHIFT(phi,-1,2))
 !Interpolate to zeta-points
 wind_div=CMPLX(.5*u+.5*CSHIFT(u,1,1),&
 &.5*v+.5*CSHIFT(v,1,2))
 
 !Test divergence
 !Interpolate u to u-points
 u=.5*CSHIFT(REAL(uv),-1,1)+.5*REAL(uv)
 !Interpolate v to v-points
 v=.5*CSHIFT(IMAG(uv),-1,2)+.5*IMAG(uv)
 !Calculate divergence in zeta-points
 div=(CSHIFT(u,1,1)-u)+(CSHIFT(v,1,2)-v)
 DO i0=1,size(div,3)
  div(:,:,i0)=div(:,:,i0)-SUM(div(:,:,i0))/dble(size(div,1)*size(div,2))
 END DO
 WRITE(*,*) 'Divergence prior/post:',SQRT(norm0),SQRT(SUM(div*div))
 
END FUNCTION wind_div


END PROGRAM create_wind
