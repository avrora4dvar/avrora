Program PCG
 USE mod_netcdf
 USE mod_BHM
 IMPLICIT NONE

 INTERFACE
  FUNCTION median(valIn) RESULT(valOut)
   IMPLICIT NONE
   REAL(8)::valIn(:),valOut,bnd(2)
   INTEGER::i0,n0,nMedian
  END FUNCTION median
  SUBROUTINE read_ens(iter_dir,ens_dir,fName,varName,val)  
   IMPLICIT NONE
   CHARACTER(len=*)::iter_dir,ens_dir,fName,varName
   CHARACTER(len=1024)::member_file
   INTEGER::ncid,n_members,s(8),flag_iter,status,i0,i1,i2
   REAL(8),allocatable::val(:,:)
   LOGICAL::flag_member
  END SUBROUTINE read_ens
  SUBROUTINE cal_eigen(M,svdU,svdS)
   IMPLICIT NONE
   REAL(8)::svdU(:,:),M(:),svdS(:)
   REAL(8),allocatable::work(:),svdVT(:,:)
   INTEGER::status,i0,i1,i2 
  END SUBROUTINE cal_eigen
  SUBROUTINE Proj(r,Br,valIn,valOut)
   IMPLICIT NONE
   REAL(8)::r(:,:),Br(:,:),valIn(:,:),valOut(:,:)
   REAL(8),allocatable::work(:),rBr(:,:)
   INTEGER,allocatable::IPIV(:)
  END SUBROUTINE Proj
  SUBROUTINE InvProj(r,valProj,valOut)
   IMPLICIT NONE
   REAL(8)::r(:,:),valProj(:,:),valOut(:,:)
  END SUBROUTINE InvProj
  SUBROUTINE kmeans(nclass,valIn,w,class_out)
   IMPLICIT NONE
   INTEGER::class_out(:),nclass,iterNo,i0,j0
   INTEGER,allocatable::class(:)
   REAL(8)::valIn(:),w(:),maxD(:),varTot
   REAL(8),allocatable::classMean(:),classW(:),classVar(:)
  END SUBROUTINE kmeans
  SUBROUTINE write_ens(iter_dir,ens_dir,fName,varName,val)
   IMPLICIT NONE
   CHARACTER(len=*)::iter_dir,ens_dir,fName,varName
   CHARACTER(len=1024)::member_file
   INTEGER::ncid,s(8),n_members,flag_iter
   INTEGER::status,dimid,varid,i0,i1,i2
   REAL(8)::val(:,:)
   LOGICAL::flag_member
  END SUBROUTINE write_ens
 END INTERFACE

!-----------------------------------------------------------------------
!Declare

 !Counters
 INTEGER::i0,i1,i2,i3,i4,j0,j1,j2,j3,j4
 INTEGER::s(8),ncid,dimid(4),varid,status

 !Lapack 
 REAL(8),allocatable::work(:)
 INTEGER::work_size(1)
 INTEGER(8),allocatable::IPIV(:)

 !Read observations
 INTEGER,allocatable::type(:)
 REAL(8),allocatable::obs(:),sig(:)
 REAL(8)::sig_max,scale_for(3),scale_obs(2),scale_ana(2)

 !Read settings
 INTEGER::n_par,n_iter,n_vec
 REAL(8)::nan
 LOGICAL::flag_out
 CHARACTER(len=1024)::out_file,x_file,x_name,r_file,r_name,iter_dir
 CHARACTER(len=1024)::obs_file,ens_dir,pcg_dir

 !Read input
 INTEGER::n_members,n_elim
 CHARACTER(len=1024)::member_file,iter_file
 REAL(8),allocatable::r0(:,:),input(:,:)
 LOGICAL::flag_member
 CHARACTER::mode
 REAL(8)::rms,norm

 !Calculate innvation vector
 REAL,allocatable::draw1(:),draw2(:,:)

 !Calculate eigenvectors
 REAL(8),allocatable::preU(:,:),preS(:),preP(:,:)

 !CG
 REAL(8)::T(1000,1000)
 REAL(8),allocatable,dimension(:,:)::x,Bx,xobs
 REAL(8),allocatable,dimension(:,:)::prer0,r,Br,robs,Ar
 REAL(8),allocatable,dimension(:,:)::PAr,Pr0,Probs,Px,Pxobs
 REAL(8),allocatable::M0(:),M1(:,:),M2(:,:)
 REAL(8),allocatable::J(:,:)
  
 !r selection
 INTEGER,allocatable::Pclass(:)
 REAL(8),allocatable::Mr0(:,:),weight(:),c0(:),c(:)
 REAL(8),allocatable::svdU(:,:),svdVT(:,:),svdS(:),PAS(:),PAU(:,:)
 REAL(8)::mean_weight

 !Write x,r
 LOGICAL::flag_file

 REAL(8),allocatable::preDia(:)
 scale_for=dble(0)
 scale_obs=dble(0)
 scale_ana=dble(0)
 scale_for=[-1.0,0.0,0.0]
!---------------------------------------------------------------------
!Read settings

 REAd(*,*) !Observation file
 REAd(*,'(A)') obs_file
 READ(*,*) !PCG file
 READ(*,'(A)') out_file
 READ(*,*) !Read input iter dir and ensemble dir 
 READ(*,'(A)') iter_dir
 READ(*,'(A)') ens_dir
 READ(*,'(A)') pcg_dir
 READ(*,*) !Read number of parallel threads
 READ(*,*) n_par
 READ(*,*) !Read r file name and r variable
 READ(*,'(A)') r_file
 READ(*,'(A)') r_name
 REAd(*,*) !Read x file name and x variable
 READ(*,'(A)') x_file
 READ(*,'(A)') x_name 
 READ(*,*) !Read maximum obs
 READ(*,*) sig_max
 READ(*,*) !Read mode
 READ(*,*) mode

 WRITE(*,*) 'Observation file:',TRIM(obs_file)
 WRITE(*,*) 'Input/Output file:',TRIM(out_file)
 WRITE(*,*) 'Iter directory:',TRIM(iter_dir)
 WRITE(*,*) 'Ens directory:',TRIM(ens_dir)
 WRITE(*,*) 'Number of threads:',n_par
 WRITE(*,*) 'r file:',TRIM(r_file)
 WRITE(*,*) 'x file:',TRIM(x_file)
 WRITE(*,*) 'mode:',mode

 !Check if output file already exist
 INQUIRE(file=TRIM(out_file),exist=flag_out)
 IF(flag_out) THEN
  WRITE(*,*) 'Output file found'
 ELSE
  WRITE(*,*) 'Starting new pcg'
  n_iter=0
 END IF

 n_iter=0
 nan=dble(1e36)
 T=dble(0)

!-----------------------------------------------------------------------
!Read previous results

IF(flag_out) THEN

 status=nf_open(TRIM(out_file),nf_nowrite,ncid)
 n_iter=INT(ncread0d(ncid,'n_iter'))+1
 n_members=INT(ncread0d(ncid,'n_members'))
 CALL ncsize(ncid,'r',s)
 n_vec=n_par*n_iter
 WRITE(*,*) 'Reading results from iteration ',n_iter-1

 !Read r
 ALLOCATE(r0(s(1),n_members+1))
 r0=ncread2d(ncid,'r0',[1,1],[s(1),n_members+1])
 ALLOCATE(r(s(1),n_vec+n_par)); r=dble(0)
 r(:,1:n_vec)=ncread2d(ncid,'r',[1,1],[s(1),n_vec])
 ALLOCATE(Br(s(1),n_vec)); Br=dble(0)
 IF(size(Br,2).GT.n_par) THEN
  Br(:,1:n_vec-n_par)=ncread2d(ncid,'Br',[1,1],[s(1),n_vec-n_par])
 END IF
 
 !Read observations
 ALLOCATE(sig(s(1)))
 sig=ncread1d(ncid,'sig',[1],[s(1)])
 ALLOCATE(obs(s(1)))
 obs=ncread1d(ncid,'obs',[1],[s(1)])
 ALLOCATE(type(s(1)))
 type=INT(ncread1d(ncid,'type',[1],[s(1)]))

 !Read preconditioner
 ALLOCATE(preU(s(1),n_members+1))
 ALLOCATE(preS(n_members+1))
 ALLOCATE(preDia(s(1)))
 preU=ncread2d(ncid,'preU',[1,1],[size(preU,1),size(preU,2)])
 preS=ncread1d(ncid,'preS',[1],[size(preS,1)])
 preDia=ncread1d(ncid,'preDia',[1],[size(preDia,1)])
 
 !Close stream
 status=nf_close(ncid); ncid=0

END IF

!----------------------------------------------------------------------
!Read input

 WRITE(*,*) 'Reading input'

 IF(ALLOCATED(input)) DEALLOCATE(input)

 IF(.NOT.flag_out) THEN
  CALL read_ens(iter_dir,ens_dir,r_file,r_name,input)
 
  !Define r
  ALLOCATE(r(size(input,1),n_par)); r=dble(0)
  n_vec=0
  n_iter=0
  n_members=size(input,2)-1
 ELSE 
  CALL read_ens('',pcg_dir,r_file,r_name,input)
 END IF 

!--------------------------------------------------------------------
!Read observations

IF(.NOT.flag_out) THEN

 WRITE(*,*) 'Reading observation file:',TRIM(obs_file)
 status=nf_open(TRIM(obs_file),nf_nowrite,ncid)
 CALL ncsize(ncid,'obs',s)

 ALLOCATE(obs(s(1)))
 obs=ncread1d(ncid,'obs',[1],[s(1)])
 ALLOCATE(sig(s(1)))
 sig=ncread1d(ncid,'sig_d',[1],[s(1)])
 ALLOCATE(type(s(1)))
 type=INT(ncread1d(ncid,'type',[1],[s(1)]))

 !Calculate obs. inflation factor
 IF(scale_for(1).LT.0) THEN
  scale_ana=[1.0,0.0]
 ELSEIF(.NOT.ANY(type.EQ.4)) THEN
  scale_ana=scale_for(1:2)
 ELSE
  !Inflation from innovation vector
  scale_obs(1)=SUM( (obs-input(:,1))**2/sig**2,type.EQ.4 )&
  &-dble(COUNT(type.EQ.4))
  scale_obs(1)=scale_obs(1)/&
  &SUM(scale_for(3)**2/sig**2,type.EQ.4)
  scale_obs(2)=2/dble(COUNT(type.EQ.4))*&
  &(scale_for(1)*SUM((obs-input(:,1))**2/sig**2,type.EQ.4 )&
  &+dble(COUNT(type.EQ.4)))**2/&
  &SUM(scale_for(3)**2/sig**2,type.EQ.4)**2
  !Posteriori estimate inflation vector
  scale_ana(1)=(scale_for(1)*scale_obs(2)+&
  &scale_obs(1)*scale_for(2))/&
  &(scale_obs(2)+scale_for(2))
  scale_ana(2)=scale_for(2)*scale_obs(2)/&
  &(scale_for(2)+scale_obs(2))
  !Reset
  scale_ana(1)=MAX(scale_ana(1),dble(0))
  !Write to output
  WRITE(*,*) 'scale for:',scale_for
  WRITE(*,*) 'scale obs:',scale_obs
  WRITE(*,*) 'scale ana:',scale_ana
 END IF

 !Rescale sigma to eliminate outliers
 WHERE( ABS(obs-input(:,1)).GT.sig_max*sig.AND..NOT.&
  &(type.GE.100.AND.type.LT.200) )
  sig=ABS(obs-input(:,1))/dble(sig_max)
 END WHERE

 
 !Set observations for average salinity
 DO i1=1,size(obs,1)
  IF(type(i1).LT.100.OR.type(i1).GE.200) CYCLE
  IF(mode.EQ.'B') THEN
   obs(i1)=input(i1,1)
  ELSE 
   obs(i1)=SUM(input(i1,:))/dble(size(input,2))
  END IF
 END DO

 !Create preconditioner M with M diagonal
 ALLOCATE(preDia(s(1))); preDia=dble(0)
 WHERE(type.GE.2.AND.type.LE.3)
  preDia=dble(COUNT(type.GE.2.AND.type.LE.3))
 END WHERE
 WHERE(type.GE.4.AND.type.LE.4)
  preDia=dble(COUNT(type.GE.4.AND.type.LE.4))
 END WHERE
 WHERE(type.GE.6.AND.type.LE.7)
  preDia=dble(COUNT(type.GE.6.AND.type.LE.7))
 END WHERE
 WHERE(type.EQ.8)
  preDia=dble(COUNT(type.EQ.8))
 END WHERE
 WHERE(type.GE.100.AND.type.LT.200)
  preDia=dble(COUNT(type.GE.100.AND.type.LT.200))
 END WHERE
 WHERE(type.GE.1000.AND.type.LE.9999)
  preDia=dble(COUNT(type.GE.1000.AND.type.LE.9999))
 END WHERE
 preDia=dble(1)/SQRT(preDia)
 preDia=dble(1)

 WRITE(*,*) 'Min/max observation:',MINVAL(obs),MAXVAL(obs)
 WRITE(*,*) 'Min/max preDia:',minval(preDia),maxval(preDia)

END IF !.not.flag_out

!-----------------------------------------------------------------------
!Calculate Br

IF(flag_out) THEN
 WRITE(*,*) 'Calculating Br'
 DO i2=1,size(input,2)
  !Apply Co^{-1/2}
  input(:,i2)=input(:,i2)/sig
  !Apply inflation
  input(:,i2)=input(:,i2)*scale_ana(1)
  !Apply preDia
  input(:,i2)=input(:,i2)*preDia
  !Save
  Br(:,n_vec-n_par+i2)=input(:,i2)
 END DO
 DEALLOCATE(input)

 !Orthogonalize
 DO i2=1,0 !size(Br,2)
  IF(i2.GT.1) THEN
   ALLOCATE(M1(i2-1,1)); M1=dble(0)
   CALL Proj(Br(:,1:i2-1),Br(:,1:i2-1)+r(:,1:i2-1),Br(:,i2:i2),M1)
   ALLOCATE(M2(size(r,1),1)); M2=dble(0)
   CALL InvProj( r(:,1:i2-1),M1,M2);  r(:,i2:i2)= r(:,i2:i2)-M2
   M2=dble(0)
   CALL InvProj(Br(:,1:i2-1),M1,M2); Br(:,i2:i2)=Br(:,i2:i2)-M2
   DEALLOCATE(M1,M2)
  END IF
  norm=SQRT( SUM((Br(:,i2)+r(:,i2))*Br(:,i2)) )
   r(:,i2)= r(:,i2)/norm
  Br(:,i2)=Br(:,i2)/norm
 END DO

END IF

!-------------------------------------------------------------------------
!Calculate normalized innovation vectors

IF(.NOT.flag_out) THEN
 WRITE(*,*) 'Calculating innovation vectors'

 !Copy input
 ALLOCATE(r0(size(input,1),size(input,2))); r0=input
 DEALLOCATE(input)

 !Create ensemble for balance operator
 IF(mode.EQ.'B') THEN
  !Gaussian draws with zero mean and standard deviation of 1
  ALLOCATE(draw2(n_members,n_members))
  ALLOCATE(draw1(SIZE(draw2)))
  CALL bhm_random_gauss(1,draw1)
  draw2=RESHAPE(draw1,[size(draw2,1),size(draw2,2)])
  !Remove mean
  draw2=draw2-SPREAD(SUM(draw2,2)/REAL(size(draw2,2)),2,size(draw2,2))
  !Normalize standard deviation
  draw2=draw2/SPREAD(SQRT(SUM(draw2**2,2)/REAL(size(draw2,2))),&
  &2,size(draw2,2))

  !Read input to Co^{-1/2}*AD*Cb*TL*Co^{-1/2}
  CALL read_ens('',ens_dir,x_file,x_name,input)
  !Calculate AD*Cb*TL*Co^{-1/2}+I for balance operator
  DO i2=1,size(input,2)
   r0(:,i2+1)=(r0(:,i2+1)/sig*scale_ana(1)+input(:,i2)*sig) !(B+I)
  END DO

  !Remove outlier due to numerical instability
  ALLOCATE(M2(size(r0,1),size(r0,2)-1)); M2=r0(:,2:size(r0,2))
  ALLOCATE(M1(size(M2,1),size(M2,2))); M1=dble(0)  
  DO i0=1,0
  j0=0
  DO WHILE(.TRUE.) 
   j0=MINVAL(type,1,type.GT.j0)
   M1=dble(0)
   DO i1=1,size(type,1)
    IF(type(i1).EQ.j0) M1(i1,:)=dble(1)
   END DO
   
   !Find median for type
   ALLOCATE(M0(COUNT(INT(M1).EQ.1)))
   M0=ABS(PACK(M2,INT(M1).EQ.1))
   rms=median(M0)
   WRITE(*,*) 'Type median/min/max',j0,rms,MINVAL(M0),MAXVAL(M0)
  
   !Set entries to far from zero to zero
   WHERE(ABS(M2).GT.sig_max*rms.AND.INT(M1).EQ.1)
    M2=dble(0)
    M1=dble(2)
   END WHERE
   WRITE(*,*) 'Type',j0,'eliminated',&
   &COUNT(INT(M1).EQ.2),COUNT(INT(M1).EQ.1)

   DEALLOCATE(M0)
   IF(j0.EQ.MAXVAL(type)) EXIT
  END DO !while
  END DO !i0
  r0(:,2:size(r0,2))=M2
  DEALLOCATE(M1,M2) 

  !Find SVD
  ALLOCATE(svdU(size(r0,1),n_members)); svdU=r0(:,2:size(r0,2))
  ALLOCATE(svdS(n_members)); svdS=dble(0)
  ALLOCATE(M0(size(r0,1))); M0=dble(1)
  CALL cal_eigen(M0,svdU,svdS)
  WRITE(*,*) 'svdS balance',SQRT(svdS)*dble(n_members)**.25
  DEALLOCATE(M0)  

  !Generate synthetic innovation vectors
  DO i2=1,n_members
    r0(:,i2+1)=svdU(:,i2)*svdS(i2)**.5*dble(n_members)**.25
    r0(:,i2+1)=obs+r0(:,i2+1)*sig
  END DO

  DEALLOCATE(input,svdU,svdS,draw1,draw2)
 END IF !mode.EQ.B

 !Normalized innovation vector
 DO i2=1,size(r0,2)
  r0(:,i2)=(obs-r0(:,i2))/sig
 END DO

 !Rescale 
 WRITE(*,*) 'Min/max normalized innovation:',MINVAL(r0(:,1)),&
 &MAXVAL(r0(:,1))

 !Calculate eigenvectors and eigenvalues of matix M*A*M
 ALLOCATE(preU(size(r0,1),size(r0,2))); preU=r0/SQRT(dble(size(preU,2)))
 ALLOCATE(preS(size(r0,2))); preS=dble(0)
 CALL cal_eigen(preDia,preU,preS)

 WRITE(*,*) 'Min/max eigenvalues M*A*M:',MINVAL(preS),MAXVAL(preS)

 !Salinity penality
 DO i1=1,size(r0,1)
  IF(type(i1).GE.100.AND.type(i1).LT.200) r0(i1,:)=dble(0)
 END DO


END IF !.NOT.flag_out

!-----------------------------------------------------------------------
!Calculate perturbations observations

ALLOCATE(robs(s(1),n_members+1)); robs=dble(0)
IF(n_members.GT.0) THEN
 !Gaussian draws with zero mean and standard deviation of 1
 ALLOCATE(draw2(size(robs,1),n_members))
 ALLOCATE(draw1(SIZE(draw2)))
 CALL bhm_random_gauss(1,draw1)
 draw2=RESHAPE(draw1,[size(draw2,1),size(draw2,2)])
 !Remove mean
 draw2=draw2-SPREAD(SUM(draw2,2)/REAL(size(draw2,2)),2,size(draw2,2))
 !Normalize standard deviation
 draw2=draw2/SPREAD(SQRT(SUM(draw2**2,2)/REAL(size(draw2,2))),&
 &2,size(draw2,2))
 robs(:,2:size(robs,2))=dble(draw2)
 DEALLOCATE(draw1,draw2)
END IF !n_members.GT.0

 !Set avg. salinity constraint
 DO i1=1,size(type)
  IF(type(i1).GE.100.AND.type(i1).LT.200) THEN
   robs(i1,:)=dble(0)
  END IF
 END DO

!-----------------------------------------------------------------------
!Determine x

 !Calculate J
 ALLOCATE(J(size(r0,2),1))
 DO i2=1,size(r0,2)
  J(i2,1)=.5*SUM( r0(:,i2)**2 )
  WRITE(*,*) 'Prior J member',i2,':',J(i2,1)
 END DO

 !Calculate normalized RMS
 DO i2=1,size(r0,2)
  WRITE(*,*) 'Prior RMS r member',i2,':',&
  SQRT( SUM(r0(:,i2)**2)/dble(size(r0,1)) )
 END DO

 !Apply preconditioner
 ALLOCATE(Mr0(size(r0,1),size(r0,2)))
 Mr0=r0*SPREAD(preDia,2,size(r0,2))
 robs=robs*SPREAD(preDia,2,size(robs,2))

 ALLOCATE(x(size(r0,1),size(r0,2)))
 ALLOCATE(xobs(size(r0,1),size(r0,2)))
 ALLOCATE(Bx(size(r0,1),size(r0,2)))
IF(.NOT.flag_out) THEN
 x=dble(0); Bx=dble(0)
ELSE
 !Calculate Ar
 ALLOCATE(Ar(size(Br,1),size(Br,2)))
 DO i2=1,size(Ar,2)
  Ar(:,i2)=Br(:,i2)+preDia*r(:,i2)*preDia
 END DO

 !Projection of r0 on r
 ALLOCATE(Pr0(size(Br,2),size(Mr0,2)))
 CALL Proj(r,Br,Mr0,Pr0)
 WRITE(*,*) 'Pr0:',Pr0(:,1)
 !Projection of robs on r
 ALLOCATE(Probs(size(Br,2),size(robs,2)))
 CALL Proj(r,Br,robs,Probs)
 WRITE(*,*) 'Min/max Probs:',MINVAL(Probs),MAXVAL(Probs)
 !Projection of Ar on r
 ALLOCATE(PAr(size(Br,2),size(Ar,2)))
 CALL Proj(r,Br,Ar,PAr)
 WRITE(*,*) 'PAr:',PAr(:,1)

 !Calculate Ritz values
 ALLOCATE(M0(size(PAr,1))); M0=dble(1)
 ALLOCATE(svdU(size(PAr,1),size(PAr,2))); svdU=PAr
 ALLOCATE(svdS(size(svdU,2))); svdS=dble(0)
 CALL cal_eigen(M0,svdU,svdS)
 WRITE(*,*) 'Ritz values:',svdS
 IF(size(svdS).GT.2.AND.MINVAL(svdS).LT..1*MAXVAL(svdS)) THEN
  WRITE(member_file,'(A,A)') TRIM(out_file),'.done'
  OPEN(unit=2,file=TRIM(member_file))
  WRITE(2,*) svdS
  CLOSE(2)
 END IF
 DEALLOCATE(svdS,svdU,M0)

 !Calculate Px
 ALLOCATE(IPIV(size(PAr,1)))
 ALLOCATE(M1(size(PAr,1),size(PAr,2))); M1=PAr
 ALLOCATE(Px(size(Pr0,1),size(Pr0,2))); Px=Pr0
 CALL dgesv(size(M1,1),size(Px,2),M1,size(M1,1),IPIV,&
 &Px,size(Px,1),status)
 DEALLOCATE(IPIV,M1)  
 WRITE(*,*) 'Min/max Px0(:,1):',MINVAL(Px(:,1)),MAXVAL(Px(:,1))

 !Calculate Pxobs
 ALLOCATE(IPIV(size(PAr,1)))
 ALLOCATE(M1(size(PAr,1),size(PAr,2))); M1=PAr
 ALLOCATE(Pxobs(size(Probs,1),size(Probs,2))); Pxobs=Probs
 CALL dgesv(size(M1,1),size(Pxobs,2),M1,size(M1,1),IPIV,&
 &Pxobs,size(Pxobs,1),status)
 DEALLOCATE(IPIV,M1)  
 WRITE(*,*) 'Min/max Pxobs:',MINVAL(Pxobs),MAXVAL(Pxobs)

 !Calculate x,xobs,Bx
 CALL InvProj(r,Px,x)
 x=x*SPREAD(preDia,2,size(x,2))
 WRITE(*,*) 'Min/max x:',MINVAL(x(:,1)),MAXVAL(x(:,1))
 CALL InvProj(Br,Px,Bx)
 Bx=Bx*SPREAD(1/preDia,2,size(x,2))
 WRITE(*,*) 'Min/max Bx:',MINVAL(Bx(:,1)),MAXVAL(Bx(:,1))
 CALL InvProj(r,Pxobs,xobs)
 xobs=xobs*SPREAD(preDia,2,size(xobs,2))
 WRITE(*,*) 'Min/max xobs:',MINVAL(xobs),MAXVAL(xobs)
 DEALLOCATE(Px,Pxobs,Pr0,PAr)

 !Calculating new residual
 Mr0=r0-(Bx+x)
 Mr0=Mr0*SPREAD(preDia,2,size(Mr0,2))
 WRITE(*,*) 'Min/max r0:',&
 &MINVAL(Mr0(:,1)/preDia),MAXVAL(Mr0(:,1)/preDia)

 !Calculate J
 DO i2=1,size(r0,2)
  J(i2,1)=.5*SUM(x(:,i2)*Bx(:,i2))&
  &+.5*SUM( (r0(:,i2)-Bx(:,i2))**2 )
  WRITE(*,*) 'Post J member',i2,':',J(i2,1)
 END DO

 !Calculate normalized RMS
 DO i2=1,size(r0,2)
  WRITE(*,*) 'Post RMS r member',i2,':',&
  SQRT(SUM((Mr0(:,i2)/preDia)**2)/dble(size(r0,1)))
 END DO

END IF !.NOT.flag_out


!--------------------------------------------------------------------
!Create new search directions

IF(n_par.EQ.1) THEN
 r(:,n_vec+1)=Mr0(:,1)
ELSE
 !Projection initial residuals on preU
 ALLOCATE(preP(size(preU,2),size(r0,2))); preP=dble(0)
 CALL Proj(preU,preU,r0,preP)
 
 !Create new projection space
 ALLOCATE(svdU(size(Mr0,1),size(preP,2))); svdU=dble(0)
 CALL DGEMM('N','T',size(Mr0,1),size(preP,1),size(Mr0,2),&
 &dble(1),Mr0,size(Mr0,1),preP,size(Mr0,2),dble(0),&
 &svdU,size(Mr0,1),status)
 IF(status.NE.0) CALL xerbla('DGEMM',status)
 DO i2=1,size(svdU,2)
  svdU(:,i2)=svdU(:,i2)/preS(i2)**2/dble(size(preU,2))
 END DO
 WRITE(*,*) 'Singular values Ctot:',preS**2
 
 !Coefficients
 ALLOCATE(c0(size(svdU,2)))
 ALLOCATE(c(size(svdU,2)))
 ALLOCATE(weight(size(svdU,2)))
 DO i2=1,size(svdU,2)
  c0(i2)=SUM( preU(:,i2)*(Mr0(:,1)/preDia) )
  c(i2)=preP(i2,1)
  weight(i2)=(1.0-1.0/preS(i2)**2)*preS(i2)**4*c(i2)**2
 END DO  
 weight=weight/SUM(weight)*dble(size(weight))
 WRITE(*,*) 'Coefficient ratio:',c0/c
 WRITE(*,*) 'Min/max weight:',MINVAL(weight),MAXVAL(weight)
 weight=MAX(weight,MINVAL(weight,1,weight.GT.0.0))

 !Cluster
 ALLOCATE(Pclass(size(preS)))
 CALL kmeans(n_par,preS**(-2)*c0/c,&
 &weight,Pclass)
 DO i1=1,n_par
  WRITE(*,*) 'Class',i1,'count',COUNT(Pclass.EQ.i1),&
  &'mean',1.0/SQRT(SUM(weight/preS**2,Pclass.EQ.i1)/&
  &dble(SUM(weight,Pclass.EQ.i1))),&
  &'SVD',MINVAL(preS,1,Pclass.EQ.i1),MAXVAL(preS,1,Pclass.EQ.i1)
  WRITE(*,*) SUM(1/preS,1,Pclass.EQ.i1)/dble(count(Pclass.EQ.i1))
 END DO
 
 !Compose new residuals
 ALLOCATE(M1(size(svdU,1),n_par)); M1=dble(0)
 ALLOCATE(M0(n_par)); M0=dble(0)
 DO i1=1,size(Pclass)
  M1(:,Pclass(i1))=M1(:,Pclass(i1))+&
  &svdU(:,i1)*preP(i1,1)
 END DO
 WRITE(*,*) 'Rel. difference 1st residual:',&
 &SQRT( SUM( (Mr0(:,1)-SUM(M1,2))**2 )/SUM( Mr0(:,1)**2 ) )

 !Force sum to be equal to residual
 M1(:,1)=M1(:,1)+(Mr0(:,1)-SUM(M1,2))

 !Approximately normalize search directions
 DO i2=1,size(M1,2)
  M0(i2)=dble(0)
  DO i0=1,size(preU,2)
   M0(i2)=M0(i2)+SUM( preU(:,i0)*M1(:,i2) )**2*(preS(i0)**2-1.0)
  END DO
  M1(:,i2)=M1(:,i2)/SQRT(M0(i2))
 END DO
 WRITE(*,*) 'Bnorm r:',SQRT(M0)

 r(:,n_vec+1:n_vec+n_par)=M1
 DEALLOCATE(M0,M1,preP,c0,c)

END IF

!-------------------------------------------------------------------
!Create output file

IF(.NOT.flag_out) THEN

 WRITE(*,*) 'Creating output file:',TRIM(out_file)
 status=nf_create(TRIM(out_file),nf_classic_model,ncid)

 !Define dimensions
 status=nf_def_dim(ncid,TRIM(r_name),size(r,1),dimid(1))
 status=nf_def_dim(ncid,'members',size(r0,2),dimid(2))
 status=nf_def_dim(ncid,'scale',2,dimid(4))
 status=nf_def_dim(ncid,'n',nf_unlimited,dimid(3))

 !Define variables  
 status=nf_def_var(ncid,'n_iter',nf_int,0,0,varid)
 status=nf_def_var(ncid,'n_members',nf_int,0,0,varid)
 status=nf_def_var(ncid,'obs',nf_double,1,[dimid(1)],varid)
 status=nf_def_var(ncid,'sig',nf_double,1,[dimid(1)],varid)
 status=nf_def_var(ncid,'type',nf_int,1,[dimid(1)],varid)
 status=nf_def_var(ncid,'preDia',nf_double,1,[dimid(1)],varid)
 status=nf_def_var(ncid,'r',nf_double,2,[dimid(1),dimid(3)],varid) 
 status=nf_def_var(ncid,'Br',nf_double,2,[dimid(1),dimid(3)],varid)
 status=nf_def_var(ncid,'r0',nf_double,2,[dimid(1),dimid(2)],varid) 
 status=nf_def_var(ncid,'J',nf_double,2,[dimid(2),dimid(3)],varid)
 status=nf_def_var(ncid,'preU',nf_double,2,[dimid(1),dimid(2)],varid)
 status=nf_def_var(ncid,'preS',nf_double,1,[dimid(2)],varid)
 status=nf_def_var(ncid,'scale',nf_double,1,[dimid(4)],varid)

 !Close stream
 status=nf_close(ncid)

END IF

!---------------------------------------------------------------------
!Write to output file

 WRITE(*,*) 'Writing to output'

 !Open write stream
 status=nf_open(TRIM(out_file),nf_write,ncid)

 CALL ncwrite0d(ncid,'n_iter',dble(n_iter))
 CALL ncwrite2d(ncid,'J',[1,n_iter+1],[size(J,1),size(J,2)],J)
 CALL ncwrite2d(ncid,'r',[1,1],[size(r,1),size(r,2)],r)

 IF(.NOT.flag_out) THEN
  CALL ncwrite0d(ncid,'n_members',dble(n_members))

  CALL ncwrite1d(ncid,'obs',[1],[size(obs,1)],obs)
  CALL ncwrite1d(ncid,'sig',[1],[size(sig,1)],sig)
  CALL ncwrite1d(ncid,'type',[1],[size(type,1)],dble(type))
  CALL ncwrite1d(ncid,'preDia',[1],[size(preDia,1)],preDia)

  CALL ncwrite2d(ncid,'r0',[1,1],[size(r0,1),size(r0,2)],r0)

  CALL ncwrite2d(ncid,'preU',[1,1],[size(preU,1),size(preU,2)],&
  &preU)
  CALL ncwrite1d(ncid,'preS',[1],[size(preS,1)],preS)
 ELSE 
  CALL ncwrite2d(ncid,'Br',[1,1],[size(Br,1),size(Br,2)],Br)
 END IF !.NOT.flag_out

 !Close stream
 status=nf_close(ncid)

!-------------------------------------------------------------------
!Write output vector C_{obs}^{-1/2}*chi such that the 4DVAR correction
!is C_{b}*H'*C_{obs}^{-1/2}

 WRITE(*,*) 'Writing x'
 
 !Add perturbations
 x=x+xobs
 
 !Apply Co^{-1/2}
 DO i2=1,size(x,2)
  x(:,i2)=x(:,i2)/sig
 END DO

 !Write output
 CALL write_ens(iter_dir,ens_dir,x_file,x_name,x)
 
!--------------------------------------------------------------------
!Compute Bp with B=C_{obs}^{-1/2}*H*C_{b}*H'*C_{obs}^{-1/2}
!Application of H,C_{b} and H' needs to be done outside this program

 WRITE(*,*) 'Writing r'
 !Apply C_{obs}^{-1/2} and write to member directories

 ALLOCATE(M1(size(r,1),n_par))
 M1=r(:,n_vec+1:n_vec+n_par)

 DO i2=1,size(M1,2)
  !Apply preconditioner
  M1(:,i2)=M1(:,i2)*preDia
  !Apply Co^{-1/2}
  M1(:,i2)=M1(:,i2)/sig
 END DO

 CALL write_ens('',pcg_dir,r_file,r_name,M1)
 DEALLOCATE(M1)

!-------------------------------------------------------------------------

 IF(scale_ana(1).EQ.dble(0)) THEN
  WRITE(*,*) 'NO DA'
 ELSE
  WRITE(*,*) 'pcg DONE'
 END IF

END PROGRAM PCG

!--------------------------------------------------------------------------

SUBROUTINE read_ens(iter_dir,ens_dir,fName,varName,val)
 USE mod_netcdf
 IMPLICIT NONE
 CHARACTER(len=*)::iter_dir,ens_dir,fName,varName
 CHARACTER(len=1024)::member_file
 INTEGER::ncid,n_members,s(8),flag_iter,status,i0,i1,i2
 REAL(8),allocatable::val(:,:)
 LOGICAL::flag_member

 !Initialize
 ncid=0; s=0; n_members=0; flag_iter=0
 IF(ALLOCATED(val)) DEALLOCATE(val)

 !Count number of members
 IF(len_trim(iter_dir).GT.0) flag_iter=1
 DO i0=1,999
  WRITE(member_file,'(A,A,I0.3,A,A)') TRIM(ens_dir),'/Member_',i0,'/',&
  &TRIM(fName)
  INQUIRE(file=TRIM(member_file),exist=flag_member)
  IF(flag_member) n_members=n_members+1
 END DO
  
 !Read members
 DO i0=1,n_members+flag_iter
  !Full filename
  IF(flag_iter.EQ.1.AND.i0.EQ.1) THEN
   WRITE(member_file,'(A,A,A)') TRIM(iter_dir),'/',TRIM(fName)
  ELSE
   WRITE(member_file,'(A,A,I0.3,A,A)') TRIM(ens_dir),'/Member_',&
   &i0-flag_iter,'/',TRIM(fName)
  END IF
  
  !Open stream
  write(*,*) trim(member_file)
  status=nf_open(TRIM(member_file),nf_nowrite,ncid)
  IF(status.NE.0) WRITE(*,*) 'PCG.F90:error opening ',TRIM(member_file)  

  !Find size
  IF(i0.EQ.1) THEN
   CALL ncsize(ncid,TRIM(varName),s)
   ALLOCATE(val(s(1),n_members+flag_iter))
  END IF  

  !Read input
  val(:,i0)=RESHAPE(ncread1d(ncid,TRIM(varName),[1],[s(1)]),[s(1)])
  WRITE(*,*) 'Min/max input ',TRIM(member_file),':',&
  &MINVAL(val(:,i0)),MAXVAL(val(:,i0))

  !Close
  status=nf_close(ncid)
  IF(status.NE.0) WRITE(*,*) 'PCG.F90:error closing ',TRIM(member_file)
 END DO
 
END SUBROUTINE read_ens
!------------------------------------------------------------------------

SUBROUTINE cal_eigen(M,svdU,svdS)
 IMPLICIT NONE
 REAL(8)::svdU(:,:),M(:),svdS(:),valIn(:,:)
 REAL(8),allocatable::work(:),svdVT(:,:) 
 INTEGER::status,i0,i1,i2

 ALLOCATE(svdVT(size(svdU,2),size(svdU,2))); svdVT=dble(0)

 !Eigenvectors and values of A
 ALLOCATE(work(5*size(svdU,1)))
 CALL DGESVD('O','N',size(svdU,1),size(svdU,2),svdU,size(svdU,1),&
 svdS,0,size(svdU,1),svdVT,size(svdVT,1),work,size(work),status)
 DEALLOCATE(work)

 !Apply M
 DO i2=1,size(svdU,2)
  svdU(:,i2)=svdU(:,i2)*M*SQRT(svdS(i2))
 END DO

 !Eigenvectors and values of M*A*M
 ALLOCATE(work(5*size(svdU,1)))
 svdS=dble(0); svdVT=dble(0)
 CALL DGESVD('O','N',size(svdU,1),size(svdU,2),svdU,size(svdU,1),&
 svdS,0,size(svdU,1),svdVT,size(svdVT,1),work,size(work),status)
 DEALLOCATE(work)
 svdS=svdS*svdS
END SUBROUTINE cal_eigen
!----------------------------------------------------------------------

SUBROUTINE Proj(r,Br,valIn,valOut)
 IMPLICIT NONE
 REAL(8)::r(:,:),Br(:,:),valIn(:,:),valOut(:,:)
 REAL(8),allocatable::work(:),rBr(:,:)
 INTEGER,allocatable::IPIV(:)
 INTEGER::status

 !r*B*r
 ALLOCATE(rBr(size(Br,2),size(Br,2)))
 CALL dgemm('T','N',size(Br,2),size(Br,2),size(Br,1),dble(1),&
 &Br,size(Br,1),r(:,1:size(Br,2)),&
 &size(r,1),dble(0),rBr,size(rBr,1),status) 
 IF(status.NE.0) WRITE(*,*) 'ERROR calculating rBr'  

 !Br'*valIn
 CALL dgemm('T','N',size(Br,2),size(valIn,2),size(Br,1),dble(1),&
 &Br,size(Br,1),valIn,&
 &size(valIn,1),dble(0),valOut,size(valOut,1),status)
 IF(status.NE.0) WRITE(*,*) 'ERROR calculating r*'
 ALLOCATE(IPIV(size(rBr,1))); 
 CALL dgesv(size(rBr,1),size(valOut,2),rBr,size(rBr,1),IPIV,&
 &valOut,size(valOut,1),status)
 IF(status.NE.0) WRITE(*,*) 'ERROR calculating inv(rBr)'
 DEALLOCATE(IPIV) 
END SUBROUTINE Proj

!------------------------------------------------------------------------

SUBROUTINE InvProj(r,valProj,valOut)
 IMPLICIT NONE
 REAL(8),intent(in)::r(:,:),valProj(:,:)
 REAL(8),intent(inout)::valOut(:,:)
 INTEGER::status
 REAL(8),allocatable::rIn(:,:)

 ALLOCATE(rIn(size(r,1),size(valProj,1)))
 rIn=r(:,1:size(valProj,1))

 CALL dgemm('N','N',size(rIn,1),size(valProj,2),size(valProj,1),&
 &dble(1),rIn,size(rIn,1),valProj,&
 &size(valProj,1),dble(0),valOut,size(rIn,1),status)   
 IF(status.NE.0) CALL xerbla('DGEMM',status)
END SUBROUTINE InvProj

!-------------------------------------------------------------------------
SUBROUTINE kmeans(nclass,valIn,w,class_out)
 IMPLICIT NONE
  INTEGER::class_out(:),nclass,iterNo,iTry,i0,j0
  INTEGER,allocatable::class(:)
  REAL(8)::valIn(:),w(:),varTot
  REAL(8),allocatable::classMean(:),classVar(:),classW(:),maxD(:)

  !Allocate
  ALLOCATE(classVar(nclass)); classVar=dble(0)
  ALLOCATE(classMean(nclass)); classMean=dble(0)
  ALLOCATE(classW(nclass)); classW=dble(0)
  ALLOCATE(maxD(size(valIn))); maxD=dble(0)
  ALLOCATE(class(size(valIn))); class=dble(0)
  varTot=dble(1e9)

  IF(nclass.EQ.1) THEN
   class=1
  ELSE
   DO iTry=1,2e2
   DO iterNo=1,5e3
    IF(iterNo.EQ.1) THEN
     DO i0=1,nclass
      classMean(i0)=valIn( MOD(iTry*i0,size(valIn))+1 )
     END DO
    END IF

    !Assign class
    DO i0=1,size(valIn)
     class(i0)=MINLOC( ABS(valIn(i0)-classMean) ,1)
    END DO

    !Create new cluster if necessary
    maxD=dble(0)
    DO i0=1,nclass
     IF(COUNT(class.EQ.i0).EQ.0) CYCLE
     maxD=MAX(maxD,(valIn-classMean(i0))**2)
    END DO
    DO i0=1,nclass
     IF(COUNT(class.EQ.i0).GT.0) CYCLE
     class(MAXLOC(maxD))=i0
     maxD(MAXLOC(maxD))=0.0
    END DO

    !Recalculate class mean
    DO i0=1,size(classMean)
     classW(i0)=SUM(w,1,class.EQ.i0)
     classMean(i0)=SUM(w*valIn,1,class.EQ.i0)/classW(i0)
     classVar(i0)=SUM(w*(valIn-classMean(i0))**2,1,class.EQ.i0)
    END DO

  END DO !iterNo
   IF(SUM(classVar).LT.varTot) THEN
    class_out=class
    varTot=SUM(classVar)
    WRITE(*,*) 'varTot',iTry,varTot
   END IF
  END DO !iTry
 
 END IF
END SUBROUTINE kmeans
!------------------------------------------------------------------------

FUNCTION median(valIn) RESULT(valOut)
 IMPLICIT NONE
 REAL(8)::valIn(:),valOut,bnd(2)
 INTEGER::n0,nMedian

 IF(mod(size(valIn,1),2).EQ.0) THEN
  nMedian=size(valIn,1)/2
 ELSE
  nMedian=(size(valIn,1)-1)/2
 END IF
 

 !Find median using bisection
 bnd(1)=MINVAL(valIn); bnd(2)=MAXVAL(valIn)
 DO WHILE(.TRUE.)
  valOut=.5*SUM(bnd)
  n0=COUNT(valIn.LT.valOut)
  IF(n0.EQ.nMedian) THEN
   EXIT
  ELSEIF(n0.LT.nMedian) THEN
   bnd(1)=valOut
  ELSEIF(n0.GT.nMedian) THEN
   bnd(2)=valOut
  END IF
 END DO

END FUNCTION median

!--------------------------------------------------------------------------
SUBROUTINE write_ens(iter_dir,ens_dir,fName,varName,val)
 USE mod_netcdf
 IMPLICIT NONE
 CHARACTER(len=*)::iter_dir,ens_dir,fName,varName
 CHARACTER(len=1024)::member_file
 INTEGER::ncid,n_members,s(8),flag_iter,i0
 INTEGER::status,dimid,varid
 REAL(8)::val(:,:)
 LOGICAL::flag_member

 !Initialize
 ncid=0; s=0; n_members=0; flag_iter=0

 !Check iter dir
 IF(len_trim(iter_dir).GT.0) flag_iter=1
  
 !Write members
 DO i0=1,size(val,2)
  !Full filename
  IF(flag_iter.EQ.1.AND.i0.EQ.1) THEN
   WRITE(member_file,'(A,A,A)') TRIM(iter_dir),'/',TRIM(fName)
  ELSE
   WRITE(member_file,'(A,A,I0.3,A,A)') TRIM(ens_dir),'/Member_',&
   &i0-flag_iter,'/',TRIM(fName)
  END IF
  
  INQUIRE(file=TRIM(member_file),exist=flag_member)
  IF(.NOT.flag_member) THEN
   status=nf_create(TRIM(member_file),nf_classic_model,&
   &ncid)
   status=nf_def_dim(ncid,TRIM(varName),size(val,1),dimid)
   status=nf_def_var(ncid,TRIM(varName),nf_double,1,dimid,varid)
   status=nf_close(ncid)
  END IF

  !Open stream
  status=nf_open(TRIM(member_file),nf_write,ncid)
  IF(status.NE.0) WRITE(*,*) 'PCG.F90:error opening ',TRIM(member_file)  

  !Write output
  WRITE(*,*) 'Min/max output ',TRIM(member_file),':',MINVAL(val(:,i0)),&
  &MAXVAL(val(:,i0))
  CALL ncwrite1d(ncid,TRIM(varName),[1],[size(val,1)],&
  &RESHAPE(val(:,i0),[size(val,1)]))

  !Close
  status=nf_close(ncid)
  IF(status.NE.0) WRITE(*,*) 'PCG.F90:error closing ',TRIM(member_file)
 END DO
 
END SUBROUTINE write_ens
