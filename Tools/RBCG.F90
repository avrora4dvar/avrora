PROGRAM RBCG

 USE mod_netcdf
 IMPLICIT NONE

!-----------------------------------------------------------------------
!DECLARE

 !Input
 CHARACTER(len=1024)::r_file,rbcg_file,obs_file,x_file
 CHARACTER(len=64)::var_name
 
 !Netcdf variables
 INTEGER::ncid,status,s(8)
 REAL(8)::nan

 !Read B*r
 REAL(8),allocatable::Br(:)

 !Read observations
 REAL(8),allocatable::obs(:),sig(:),type(:)
 INTEGER,allocatable::mask(:)
 REAL(8)::sig_max,c
 INTEGER::K

 !Create netcdf file/update
 INTEGER::dim_id(3),var_id(20),n,i0,i1,i2
 REAL(8),allocatable,dimension(:)::rhoOld,rhoNew,&
 &JOld,JNew,aOld,aNew,bOld,bNew,ioVector
 REAL(8),allocatable,dimension(:,:)::xOld,xNew,BxOld,BxNew,&
 rOld,rNew,BrOld,BrNew,p,Bp,rInit,BrAll,rAll
 LOGICAL::flag_exist

 !Output x
 INTEGER::rmin
 REAL(8),allocatable::xmin(:,:)
!---------------------------------------------------------------------
!READ INPUT

 WRITE(*,*) 'Reading input'
 READ(*,*) !Observation file
 READ(*,'(A)') obs_file
 READ(*,*) !RBCG file
 READ(*,'(A)') rbcg_file
 READ(*,*) !Br file
 READ(*,'(A)') r_file
 READ(*,*) !Output file with x
 READ(*,'(A)') x_file
 READ(*,*) !variable name in Br file
 READ(*,'(A)') var_name
 READ(*,*) !Ceiling for which obs cost-function reaches maximum
 READ(*,*) sig_max 

 nan=dble(1e36)
!----------------------------------------------------------------------
!READ observation vector
 
 WRITE(*,*) 'Reading observation vector'
 status=nf_open(TRIM(obs_file),nf_nowrite,ncid)
 CALL ncsize(ncid,'obs',s)
 K=s(1)

 ALLOCATE(obs(K))
 ALLOCATE(type(K))
 ALLOCATE(sig(K))

 obs=ncread1d(ncid,'obs',[1],[K])
 type=ncread1d(ncid,'type',[1],[K])
 sig=ncread1d(ncid,'sig_d',[1],[K])

 status=nf_close(ncid)

!----------------------------------------------------------------------
!READ (Co^(-1/2)*Ts*Tl*Cb*Ad*As*Co^(-1/2))*r with r residual

 WRITE(*,*) 'Reading residual'
 status=nf_open(TRIM(r_file),nf_nowrite,ncid)
 ALLOCATE(ioVector(K))
 ioVector=ncread1d(ncid,TRIM(var_name),[1],[K])
 status=nf_close(ncid)

!----------------------------------------------------------------------
!CREATE NETCDF CONTAINING DATA
 
 INQUIRE(file=TRIM(rbcg_file),exist=flag_exist)
 IF(.NOT.flag_exist) THEN
  WRITE(*,*) 'Create new RBCG-file'

  !Set average salt patches to zero
  WHERE(type.GE.100.AND.type.LT.200)
   obs=ioVector  
  END WHERE

  !Create file
  status=nf_create(TRIM(rbcg_file),nf_noclobber,ncid)

  !Define dimensions with K length vector in observation space
  !and n the number of iterations
  status=nf_def_dim(ncid,'K',K,dim_id(1))
  status=nf_def_dim(ncid,'n',nf_unlimited,dim_id(2))

  !Define variables
  status=nf_def_var(ncid,'r',nf_double,2,[dim_id(1),dim_id(2)],&
  &var_id(1))
  status=nf_def_var(ncid,'Br',nf_double,2,[dim_id(1),dim_id(2)],&
  &var_id(2))
  status=nf_def_var(ncid,'x',nf_double,2,[dim_id(1),dim_id(2)],&
  &var_id(3))
  status=nf_def_var(ncid,'Bx',nf_double,2,[dim_id(1),dim_id(2)],&
  &var_id(4))
  status=nf_def_var(ncid,'p',nf_double,2,[dim_id(1),dim_id(2)],&
  &var_id(5))
  status=nf_def_var(ncid,'Bp',nf_double,2,[dim_id(1),dim_id(2)],&
  &var_id(6))
  status=nf_def_var(ncid,'J',nf_double,1,dim_id(2),var_id(7))
  status=nf_def_var(ncid,'rNorm2',nf_double,1,dim_id(2),var_id(8))
  status=nf_def_var(ncid,'obs',nf_double,1,dim_id(1),var_id(9))
  status=nf_def_var(ncid,'type',nf_double,1,dim_id(1),var_id(10))
  status=nf_def_var(ncid,'sig',nf_double,1,dim_id(1),var_id(11))
  status=nf_def_var(ncid,'mask',nf_int,1,dim_id(1),var_id(12))
  status=nf_def_var(ncid,'alpha',nf_double,1,dim_id(2),var_id(13))
  status=nf_def_var(ncid,'beta',nf_double,1,dim_id(2),var_id(14))
  status=nf_def_var(ncid,'n',nf_int,0,0,var_id(15))

  !Close file
  status=nf_close(ncid)

  !Create 1st input
  ALLOCATE(mask(K))
  ALLOCATE(rNew(K,1)); ALLOCATE(BrNew(K,1))
  ALLOCATE(xNew(K,1)); ALLOCATE(BxNew(K,1))
  ALLOCATE(p(K,1)); ALLOCATE(Bp(K,1))
  ALLOCATE(JNew(1)); ALLOCATE(rhoNew(1))
  ALLOCATE(aNew(1)); ALLOCATE(bNew(1))

  !Adjust observational uncertainty to prevent overfitting
  mask=1
  WHERE( ABS(obs-ioVector)/sig.GT.sig_max )
   sig=ABS(obs-ioVector)/sig_max
   mask=0
  END WHERE

  !Write first input to file
  status=nf_open(TRIM(rbcg_file),nf_write,ncid)
  WRITE(*,*) 'Writing residual 0'
  CALL ncwrite0d(ncid,'n',dble(0))
  CALL ncwrite1d(ncid,'obs',[1],[K],obs)
  CALL ncwrite1d(ncid,'type',[1],[K],type)
  CALL ncwrite1d(ncid,'sig',[1],[K],sig)

  !Values for n=1
  rNew=RESHAPE((obs-ioVector)/sig,[K,1]); BrNew=nan
  xNew=DBLE(0); BxNew=DBLE(0)
  p=nan; Bp=nan
  rhoNew=nan
  aNew=nan; bNew=dble(0)  
  Jnew=0.5*SUM(xNew*BxNew)+&
  &0.5*SUM( (rNew-BxNew)**2 )

  !Write to output
  CALL ncwrite2d(ncid,'r',[1,1],[K,1],rNew)
  CALL ncwrite2d(ncid,'Br',[1,1],[K,1],BrNew)
  CALL ncwrite2d(ncid,'x',[1,1],[K,1],xNew)
  CALL ncwrite2d(ncid,'Bx',[1,1],[K,1],BxNew)
  CALL ncwrite2d(ncid,'p',[1,1],[K,1],p)
  CALL ncwrite2d(ncid,'Bp',[1,1],[K,1],Bp)
  CALL ncwrite1d(ncid,'J',[1],[1],JNew)
  CALL ncwrite1d(ncid,'rNorm2',[1],[1],rhoNew)
  CALL ncwrite1d(ncid,'mask',[1],[K],dble(mask))
  CALL ncwrite1d(ncid,'alpha',[1],[1],aNew)
  CALL ncwrite1d(ncid,'beta',[1],[1],bNew)  

  status=nf_close(ncid)
 END IF

!------------------------------------------------------------------
!UPDATE VALUES

 IF(flag_exist) THEN
  WRITE(*,*) 'Updating previous RBCG-file' 
  !open stream
  status=nf_open(TRIM(rbcg_file),nf_write,ncid)

  !get number of iterations stored
  CALL ncsize(ncid,'r',s)
  n=INT(ncread0d(ncid,'n'))+1
  WRITE(*,*) 'Calculating residual ',n

  !read last entries
  ALLOCATE(rAll(K,n))
  ALLOCATE(BrAll(K,n))

  ALLOCATE(rOld(K,1))
  rOld=ncread2d(ncid,'r',[1,n],[K,1])
  ALLOCATE(BrOld(K,1))
  ALLOCATE(xOld(K,1))
  xOld=ncread2d(ncid,'x',[1,n],[K,1])
  ALLOCATE(BxOld(K,1))
  BxOld=ncread2d(ncid,'Bx',[1,n],[K,1])
  ALLOCATE(p(K,1))
  ALLOCATE(Bp(K,1))
  ALLOCATE(rInit(K,1))
  rInit=ncread2d(ncid,'r',[1,1],[K,1])
  ALLOCATE(JOld(n))
  JOld=ncread1d(ncid,'J',[1],[n])
  ALLOCATE(rhoOld(n))
  rhoOld=ncread1d(ncid,'rNorm2',[1],[n])
  ALLOCATE(mask(K))
  mask=INT(ncread1d(ncid,'mask',[1],[K]))
  ALLOCATE(bOld(n))
  bOld=ncread1d(ncid,'beta',[1],[n])
  ALLOCATE(aOld(n))
  aOld=ncread1d(ncid,'alpha',[1],[n])
  sig=ncread1d(ncid,'sig',[1],[K])

  !Calculate entries for new iteration
  ALLOCATE(rNew(K,1)); rNew=nan
  ALLOCATE(BrNew(K,1)); BrNew=nan
  ALLOCATE(xNew(K,1)); xNew=nan
  ALLOCATE(BxNew(K,1)); BxNew=nan
  ALLOCATE(JNew(n+1)); JNew=nan
  ALLOCATE(rhoNew(n+1)); rhoNew=nan
  ALLOCATE(aNew(n+1)); aNew=nan
  ALLOCATE(bNew(n+1)); bNew=dble(0)

  !Complete calculation Co^(-1/2)*Ts*Tl*Cb*Ad*As*Co^(-1/2)
  BrOld=RESHAPE(ioVector/sig,[K,1])
  CALL ncwrite2d(ncid,'Br',[1,n],[K,1],BrOld)

  !Orthogonalize
  rAll=ncread2d(ncid,'r',[1,1],[K,n])
  BrAll=ncread2d(ncid,'Br',[1,1],[K,n])
  DO i2=n-1,1,-1
   c=SUM(rAll(:,n)*BrAll(:,i2))/SUM(rAll(:,i2)*BrAll(:,i2))
   write(*,*) 'cr:',n,i2,c
   rAll(:,n)=rAll(:,n)-c*rAll(:,i2)
   c=SUM(BrAll(:,n)*rAll(:,i2))/SUM(rAll(:,i2)*BrAll(:,i2))
   write(*,*) 'cBr:',n,i2,c
   BrAll(:,n)=BrAll(:,n)-c*BrAll(:,i2)   
  END DO
  rOld=RESHAPE(rAll(:,n),SHAPE(rOld))
  BrOld=RESHAPE(BrAll(:,n),SHAPE(BrOld))
  DEALLOCATE(rAll,BrAll)
  CALL ncwrite2d(ncid,'r',[1,n],[K,1],rOld)
  CALL ncwrite2d(ncid,'Br',[1,n],[K,1],BrOld)

  !Complete calculation ||r_n||^2
  rhoOld(n)=SUM(rOld*BrOld)
  rhoNew(1:n)=rhoOld
  CALL ncwrite1d(ncid,'rNorm2',[1],[n+1],rhoNew)

  !Complete calculate beta_{n-1}
  IF(n.GT.1) THEN
   bOld(n-1)=rhoNew(n)/rhoNew(n-1)
  END IF
  bNew(1:n)=bOld
  CALL ncwrite1d(ncid,'beta',[1],[n+1],bNew)
 
  !Complete calculation p_n, Bp_n
  IF(n.EQ.1) THEN
   p=rOld
   Bp=BrOld
  ELSE
   p=ncread2d(ncid,'p',[1,n-1],[K,1])
   Bp=ncread2d(ncid,'Bp',[1,n-1],[K,1])

   p=rOld+bOld(n-1)*p
   Bp=BrOld+bOld(n-1)*Bp
  END IF
  CALL ncwrite2d(ncid,'p',[1,n],[K,1],p) 
  CALL ncwrite2d(ncid,'Bp',[1,n],[K,1],Bp)

  !Complete calculation alpha_n
  aOld(n)=rhoOld(n)/SUM( (Bp+p)*Bp )
  aNew(1:n)=aOld
  CALL ncwrite1d(ncid,'alpha',[1],[n+1],aNew)

  !Calculate new solution x_{n+1}
  xNew=xOld+aNew(n)*p
  CALL ncwrite2d(ncid,'x',[1,n+1],[K,1],xNew)

  !Calculate Bx_{n+1}
  BxNew=BxOld+aNew(n)*Bp
  CALL ncwrite2d(ncid,'Bx',[1,n+1],[K,1],BxNew)

  !Calculate new residual
  rNew=rOld-aNew(n)*(Bp+p)
  CALL ncwrite2d(ncid,'r',[1,n+1],[K,1],rNew)

  !Calculate J_{n+1}
  JNew(n+1)=0.5*SUM(xNew*BxNew)+&
  &0.5*SUM( (rInit-BxNew)**2 )
  JNew(1:n)=JOld
  CALL ncwrite1d(ncid,'J',[1],[n+1],Jnew)

  !New iteration number 
  CALL ncwrite0d(ncid,'n',dble(n))

  WRITE(*,*) 'Residual calculated'

  status=nf_close(ncid)
 END IF

!------------------------------------------------------------------
!OUTPUT NEW RESIDUAL Co^(-1/2)*r_{n+1}

 WRITE(*,*) 'Write new residual'
 ioVector=RESHAPE(rNew,[K])/sig
 status=nf_open(TRIM(r_file),nf_write,ncid)
 CALL ncwrite1d(ncid,TRIM(var_name),[1],[K],ioVector)
 status=nf_close(ncid)

!---------------------------------------------------------------------------
!Output new x

 WRITE(*,*) 'Write new x'

 INQUIRE(file=TRIM(x_file),exist=flag_exist)
 IF(.NOT.flag_exist) THEN
   status=nf_create(TRIM(x_file),nf_classic_model,ncid)
   status=nf_def_dim(ncid,'K',K,dim_id(1))
   status=nf_def_var(ncid,TRIM(var_name),nf_double,1,dim_id(1),&
   &var_id(1))
   status=nf_close(ncid)
 END IF   

 ioVector=RESHAPE(xNew,[K])/sig

 status=nf_open(TRIM(x_file),nf_write,ncid)
 CALL ncwrite1d(ncid,TRIM(var_name),[1],[K],ioVector)
 status=nf_close(ncid)
 
 WRITE(*,*) 'RBCG done'

END PROGRAM RBCG
