PROGRAM cov_scale 
 USE mod_netcdf
 USE mod_BHM
 IMPLICIT NONE

!----------------------------------------------------------------------------
!Variables

 !Input
 CHARACTER(len=1024)::rbcg_file,out_file,member_name,ens_dir
 CHARACTER(len=64)::var_name

 !RBCG
 REAL(8),allocatable,dimension(:,:)::r,Br
 REAL(8),allocatable,dimension(:)::alpha,beta,rNorm,obs,sig
 REAL(8),allocatable::matT(:,:)
 REAL(8),allocatable,dimension(:)::x,x1,x2,Bx1
 INTEGER,allocatable::mask(:)
 REAL(8),allocatable::tmp1(:)

 !Scale
 REAL,allocatable::scaleFac(:),tmp2(:)
 REAL(8)::chi2(2),dof
 REAL(8),allocatable,dimension(:,:,:,:)::temp,salt,u,v
 REAL(8),allocatable,dimension(:,:,:)::zeta
 CHARACTER(len=1024)::member_file
 LOGICAL::flag_exist

 !Netcdf/counters
 INTEGER::status,ncid,s(8)
 INTEGER::dimid(10),varid(10)
 INTEGER::i0,i1,i2,i3,i4
 REAL(8)::nan


!----------------------------------------------------------------------------
!Input

 READ(*,*) !Input RBCG file
 READ(*,'(A)') rbcg_file
 READ(*,*) !Ensemble directory
 READ(*,'(A)') ens_dir
 READ(*,*) !Name of file with ensemble perturbation
 READ(*,'(A)') member_name
 READ(*,*) !Name of output file
 READ(*,'(A)') out_file
 READ(*,*) !Name variable that is input for ad_sample
 READ(*,'(A)') var_name

 nan=dble(1e36)
!----------------------------------------------------------------------------
! RBCG file
 
  !Open stream to file
  WRITE(*,*) 'Reading RBCG file: ',TRIM(rbcg_file)
  status=nf_open(TRIM(rbcg_file),nf_nowrite,ncid)

  !Get size
  CALL ncsize(ncid,'r',s)
  ALLOCATE(alpha(s(2)-1))
  ALLOCATE(beta(s(2)-2))
  ALLOCATE(r(s(1),s(2)-1))
  ALLOCATE(Br(s(1),s(2)-1))
  ALLOCATE(rNorm(s(2)-1))
  ALLOCATE(sig(s(1)))
  ALLOCATE(obs(s(1)))
  ALLOCATE(mask(s(1)))
  ALLOCATE(Bx1(s(1)))
  ALLOCATE(x1(s(1)))
  ALLOCATE(x2(s(1)))
  ALLOCATE(x(s(1)))

  !Read
  mask=INT( ncread1d(ncid,'mask',[1],[s(1)]) )
  obs=ncread1d(ncid,'obs',[1],[s(1)])
  sig=ncread1d(ncid,'sig',[1],[s(1)])

  alpha=ncread1d(ncid,'alpha',[1],[s(2)-1])
  beta=ncread1d(ncid,'beta',[1],[s(2)-2])
  Bx1=RESHAPE(ncread2d(ncid,'Bx',[1,s(2)],[s(1),1]),[s(1)])
  x1=RESHAPE(ncread2d(ncid,'x',[1,s(2)],[s(1),1]),[s(1)])
  r=ncread2d(ncid,'r',[1,1],[s(1),s(2)-1])
  Br=ncread2d(ncid,'Br',[1,1],[s(1),s(2)-1])
  rNorm=SQRT(ncread1d(ncid,'rNorm2',[1],[s(2)-1]))

  !Close stream
  status=nf_close(ncid)

  !Normalize
  DO i0=1,size(rNorm)
   r(:,i0)=r(:,i0)/rNorm(i0)
   Br(:,i0)=Br(:,i0)/rNorm(i0)
  END DO
!-------------------------------------------------------------------------
!Construct Lanczos approximation T^{-1} with
! (B+I) = r*T*Br'
! B=Co^{-1/2}*Ts*Tl*Cb*Ad*As*Co^{-1/2}

 WRITE(*,*) 'Calculating Lanczos approximation'

 !T(:,2) is main diagonal T, T(:,1) first diagonal below main,T(:,3)
 !first above
 ALLOCATE(matT(size(rNorm),3)); matT=dble(0)
 DO i1=1,size(matT,1)
  matT(i1,2)=dble(1)/alpha(i1)
  IF(i1.GT.1) matT(i1,2)=matT(i1,2)+beta(i1-1)/alpha(i1-1)
  matT(i1,1)=-SQRT(beta(i1))/alpha(i1)
  matT(i1,3)=matT(i1,1)
 END DO
 
!------------------------------------------------------------------------
!Calculate the 

  !(B+I)^{-1}*(Co^{-1/2}*d)
  x1=x1

  !(B+I)^{-2}*(Co^{-1/2}*d) 
  ALLOCATE(tmp1(size(Br,2))); tmp1=dble(0)
  CALL dgemv('T',size(Br,1),size(Br,2),dble(1),Br,size(Br,1),&
  &x1,1,dble(0),tmp1,1)
  CALL dptsv(SIZE(matT,1),1,matT(:,2),matT(1:size(matT,1)-1,1),&
  &tmp1,size(tmp1),status)
  CALL dgemv('N',size(r,1),size(r,2),dble(1),r,size(r,1),&
  tmp1,1,dble(0),x2,1)
  DEALLOCATE(tmp1)

  WRITE(*,*) '2-norm x1,x2:',SQRT(SUM(x1*x1)),SQRT(SUM(x2*x2)),&
  SUM(x1*x2)/SQRT(SUM(x1*x1)*SUM(x2*x2))

  !d*Co^{-1/2}(B+I)^{-1}*Co^{-1/2}*d
  chi2(1)=rNorm(1)*SUM(Br(:,1)*x1,mask.NE.0)
  
  !d*Co^{-1/2}(B+I)^{-2}*Co^{-1/2}*d
  chi2(2)=rNorm(1)*SUM(Br(:,1)*x2,mask.NE.0)

  WRITE(*,*) 'chi2:',COUNT(mask.NE.0),chi2
!--------------------------------------------------------------------------
!Draw scale factor

 ALLOCATE(scaleFac(1000))
 ALLOCATE(tmp1(1000))
 
 !Draw possible scale factors from inverse gamma distribution with mode 1
 CALL bhm_random_IG(1,50.0,49.8,scaleFac)
 tmp1=chi2(1)/scaleFac-(1.0-scaleFac)/scaleFac**2*chi2(2)
 tmp1=MAX(tmp1,dble(1)/nan)

 !Weight each draw according to chi2 distribution
 dof=dble(COUNT(mask.NE.0))
 WRITE(*,*) 'draws below/above dof:',COUNT(tmp1.LT.dof),&
 &COUNT(tmp1.GT.dof)
 tmp1=(.5*dof-1.0)*LOG(tmp1/chi2(1))-.5*(tmp1-chi2(1))
 tmp1=MAX(MIN(tmp1,dble(700)),dble(-700))
 tmp1=exp(tmp1)
 WRITE(*,*) 'Min/max chi2 weights:',MINVAL(tmp1),MAXVAL(tmp1)
 
 !Draw one scale factor
 DO i1=2,size(tmp1)
   tmp1(i1)=tmp1(i1)+tmp1(i1-1)
 END DO
 tmp1=tmp1/MAXVAL(tmp1)
 ALLOCATE(tmp2(1))
 CALL bhm_random_uniform(1,tmp2)
 scaleFac(1)=scaleFac(MINLOC(tmp1,1,tmp1.GT.tmp2(1))) 

 DEALLOCATE(tmp1,tmp2)

 WRITE(*,*) 'scale factor:',scaleFac(1)
!-------------------------------------------------------------------------
!Rescale ensemble

 DO i1=1,999
  WRITE(member_file,'(A,A,I0.3,A,A)') TRIM(ens_dir),'/Member_',&
  &i1,'/',TRIM(member_name)
  INQUIRE(file=TRIM(member_file),exist=flag_exist)
  IF(.NOT.flag_exist) CYCLE

  !Open stream
  WRITE(*,*) 'Scaling ',TRIM(member_file)
  status=nf_open(TRIM(member_file),nf_write,ncid)
  
  IF(i1.EQ.1) THEN
   CALL ncsize(ncid,'temp',s)
   ALLOCATE(temp(s(1),s(2),s(3),1))
   ALLOCATE(salt(s(1),s(2),s(3),1))
   ALLOCATE(zeta(s(1),s(2),1))
   ALLOCATE(u(s(1)-1,s(2),s(3),1))
   ALLOCATE(v(s(1),s(2)-1,s(3),1))
  END IF
  
  !Read grid
  temp=ncread4d(ncid,'temp',[1,1,1,1],SHAPE(temp))
  salt=ncread4d(ncid,'salt',[1,1,1,1],SHAPE(salt))
  zeta=ncread3d(ncid,'zeta',[1,1,1],SHAPE(zeta))
  u=ncread4d(ncid,'u',[1,1,1,1],SHAPE(u))
  v=ncread4d(ncid,'v',[1,1,1,1],SHAPE(v))
  
  !Scale grid
  temp=SQRT(scaleFac(1))*temp
  salt=SQRT(scaleFac(1))*salt
  zeta=SQRT(scaleFac(1))*zeta
  u=SQRT(scaleFac(1))*u
  v=SQRT(scaleFac(1))*v

  !Write scaled fields
  CALL ncwrite4d(ncid,'temp',[1,1,1,1],SHAPE(temp),temp)
  CALL ncwrite4d(ncid,'salt',[1,1,1,1],SHAPE(salt),salt)
  CALL ncwrite4d(ncid,'u',[1,1,1,1],SHAPE(u),u)
  CALL ncwrite4d(ncid,'v',[1,1,1,1],SHAPE(v),v)
  CALL ncwrite3d(ncid,'zeta',[1,1,1],SHAPE(zeta),zeta)
  
  !Close stream
  status=nf_close(ncid)
 END DO

 IF(ALLOCATED(temp)) DEALLOCATE(temp,salt,u,v,zeta)
!--------------------------------------------------------------------------
!Create output

 WRITE(*,*) 'Creating ',TRIM(out_file)

 !Create file
 status=nf_create(TRIM(out_file),nf_clobber,ncid)
 status=nf_def_dim(ncid,'K',size(x1),dimid(1))
 status=nf_def_var(ncid,TRIM(var_name),nf_double,1,dimid(1:1),varid(1))
 status=nf_def_var(ncid,'x1',nf_double,1,dimid(1:1),varid(2))
 status=nf_def_var(ncid,'x2',nf_double,1,dimid(1:1),varid(3))
 status=nf_def_var(ncid,'chi1',nf_double,0,0,varid(4))
 status=nf_def_var(ncid,'chi2',nf_double,0,0,varid(5))
 status=nf_def_var(ncid,'scaleFac',nf_double,0,0,varid(6))
 status=nf_close(ncid)

 !Calculate output
 WHERE(mask.EQ.0)
  sig=nan
 END WHERE
 x1=x1/sig
 x2=x2/sig
 
 ALLOCATE(x(size(x1)))
 x=x1/scaleFac(1)-x2*(1.0-scaleFac(1))/scaleFac(1)**2

 !Write to output
 status=nf_open(TRIM(out_file),nf_write,ncid)
 CALL ncwrite1d(ncid,TRIM(var_name),[1],[size(x)],x)
 CALL ncwrite1d(ncid,'x1',[1],[size(x1)],x1)
 CALL ncwrite1d(ncid,'x2',[1],[size(x2)],x2)
 CALL ncwrite0d(ncid,'scaleFac',dble(scaleFac(1)))
 CALL ncwrite0d(ncid,'chi1',chi2(1))
 CALL ncwrite0d(ncid,'chi2',chi2(2))

 !Close stream
 status=nf_close(ncid)
 
!----------------------------------------------------------------------------
! Close
  
  WRITE(*,*) 'cov_scale DONE'

END PROGRAM cov_scale
