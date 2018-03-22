PROGRAM ini_scale
 USE mod_netcdf
 USE mod_BHM
 IMPLICIT NONE

 !Input
 CHARACTER(len=1024)::rbcg_file,scale_file,sample_file,out_file
 CHARACTER(len=64)::var_name

 !RBCG and sample
 REAL(8),allocatable,dimension(:,:)::r,Br
 REAL(8),allocatable,dimension(:)::obs,sig,sample,alpha,beta,rNorm,Bx
 REAL(8),allocatable,dimension(:)::oa,ob,ab,cor,corVec
 REAl(8)::cor1
 LOGICAL,allocatable::maskCor(:)
 INTEGER,allocatable::mask(:),type(:)

 !Matrix inversions
 REAL(8),allocatable,dimension(:,:)::matT
 REAL(8),allocatable,dimension(:)::tmp1,Psample,x1,x2
 REAL,allocatable,dimension(:)::draw

 !Scaling
 REAL(8)::scaleFac
 LOGICAL::flag_exist

 !Netcdf
 INTEGER::ncid,status,s(8),dimid,varid
 INTEGER::i0,i1,i2,i3,i4
 REAL(8)::nan

!------------------------------------------------------------------------
!Input 

 READ(*,*) !RBCG file
 READ(*,'(A)') rbcg_file
 READ(*,*) !File with scaling factor
 READ(*,'(A)') scale_file
 READ(*,*) !Output from Ts*Tl*perturbation field
 READ(*,'(A)') sample_file
 READ(*,*) !Output file
 READ(*,'(A)') out_file
 READ(*,*) !Variable name in output file
 READ(*,'(A)') var_name

 nan=dble(1e36)
!-----------------------------------------------------------------------
!Read RBCG file

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
 ALLOCATE(Bx(s(1)))
 ALLOCATE(type(s(1)))
 ALLOCATE(oa(s(1)))
 ALLOCATE(ob(s(1)))
 ALLOCATE(ab(s(1)))
 
 !Read
 mask=ncread1d(ncid,'mask',[1],[s(1)])
 obs=ncread1d(ncid,'obs',[1],[s(1)])
 sig=ncread1d(ncid,'sig',[1],[s(1)])
 type=INT(ncread1d(ncid,'type',[1],[s(1)]))
 Bx=RESHAPE( ncread2d(ncid,'Bx',[1,s(2)],[s(1),1]) ,[s(1)])

 alpha=ncread1d(ncid,'alpha',[1],[s(2)-1])
 beta=ncread1d(ncid,'beta',[1],[s(2)-2])
 r=ncread2d(ncid,'r',[1,1],[s(1),s(2)-1])
 Br=ncread2d(ncid,'Br',[1,1],[s(1),s(2)-1])
 rNorm=SQRT(ncread1d(ncid,'rNorm2',[1],[s(2)-1]))

 !Close stream
 status=nf_close(ncid)

 !Difference observations and prior control run
 ob=sig*r(:,1)
 !Difference post-prior control run
 !(control run uses results form 1 extra CG iteration)
 ab=sig*Bx 
 !Difference observations and post control run
 oa=ob-ab

 !Normalize
 DO i0=1,size(rNorm)
  r(:,i0)=r(:,i0)/rNorm(i0)
  Br(:,i0)=Br(:,i0)/rNorm(i0)
 END DO

!----------------------------------------------------------------------
!Read Ts*Tl*perturbation

 WRITE(*,*) 'Read samples form ',TRIM(sample_file)

 !Read difference ensemble member-background observations 
 status=nf_open(TRIM(sample_file),nf_nowrite,ncid)
 CALL ncsize(ncid,TRIM(var_name),s)
 ALLOCATE(sample(s(1)))
 sample=ncread1d(ncid,TRIM(var_name),[1],[s(1)])
 status=nf_close(ncid)

 !Create member innovation vector
 sample=ob-sample
 mask=1

 !Apply Co^{-1/2} to samples from perturbation file
 sample=sample/sig
 WRITE(*,*) 'Min/max normalized input:',minval(sample),maxval(sample)

 !Salinity conservation
 WHERE(type.GE.100.AND.type.LT.200)
  sample=0.0
 END WHERE

 !Add perturbation observations
 ALLOCATE(draw(size(sample)))
 CALL bhm_random_gauss(1,draw)
 WRITE(*,*) 'Min/max normalized observations:',minval(draw),maxval(draw)
 sample=sample+draw
 DEALLOCATE(draw)

 WRITE(*,*) 'Min/max normalized member innovation vector:',&
 &minval(sample),maxval(sample)
 WRITE(*,*) 'Norm normalized member innovation vector:',&
 &SQRT(SUM(sample*sample,mask.NE.0))


!--------------------------------------------------------------------------
!Calculate the perturbation in the innovation factor
! (B+I) = r*T*Br'
! B=Co^{-1/2}*Ts*Tl*Cb*Ad*As*Co^{-1/2}

 WRITE(*,*) 'Calculating member correction'

 !Calculate projection of innovation vector on vectors in r
 ALLOCATE(Psample(size(r,2)))
 CALL dgemv('T',size(Br,1),size(Br,2),dble(1),Br,size(Br,1),&
 &sample,1,dble(0),Psample,1)

 !T(:,2) is main diagonal T, T(c:,1) first diagonal below main,T(:,3)
 !first above
 ALLOCATE(matT(size(rNorm),3)); matT=dble(0)
 DO i1=1,size(matT,1)
  matT(i1,2)=dble(1)/alpha(i1)
  IF(i1.GT.1) matT(i1,2)=matT(i1,2)+beta(i1-1)/alpha(i1-1)
  matT(i1,1)=-SQRT(beta(i1))/alpha(i1)
  matT(i1,3)=matT(i1,1)
 END DO
 
 !Appply (B+I)^{-1}
 CALL dptsv(SIZE(matT,1),1,matT(:,2),matT(1:size(matT,1)-1,1),&
 &Psample,size(Psample),status)

 !Convert (B+I)^{-1}*(Co^{-1/2}*d) back to full observation space
 ALLOCATE(x1(size(sample)))
 CALL dgemv('N',size(r,1),size(r,2),dble(1),r,size(r,1),&
 Psample,1,dble(0),x1,1)
 WRITE(*,*) 'Norm (B+I)^{-1}*d:',SQRT(SUM(x1*x1,mask.NE.0))

 !T(:,2) is main diagonal T, T(:,1) first diagonal below main,T(:,3)
 !first above
 matT=dble(0)
 DO i1=1,size(matT,1)
  matT(i1,2)=dble(1)/alpha(i1)
  IF(i1.GT.1) matT(i1,2)=matT(i1,2)+beta(i1-1)/alpha(i1-1)
  matT(i1,1)=-SQRT(beta(i1))/alpha(i1)
  matT(i1,3)=matT(i1,1)
 END DO

 !Appply (B+I)^{-2}
 CALL dptsv(SIZE(matT,1),1,matT(:,2),matT(1:size(matT,1)-1,1),&
 &Psample,size(Psample),status)

 !Convert  (B+I)^{-2}*(Co^{-1/2}*d) back to full observation space
 ALLOCATE(x2(size(sample)))
 CALL dgemv('N',size(r,1),size(r,2),dble(1),r,size(r,1),&
 Psample,1,dble(0),x2,1)
 WRITE(*,*) 'Norm (B+I)^{-2}*d:',SQRT(SUM(x2*x2,mask.NE.0))

!--------------------------------------------------------------------
!Calculate 
! Tl*C_{b}*Ad*Co^{-1/2}(B+A)^{-2}Co^{-1/2}ob

 WRITE(*,*) 'Calculating correction factor'

 !Calculate projection of innovation vector on vectors in r
 Psample=dble(0.0)
 CALL dgemv('T',size(Br,1),size(Br,2),dble(1),Br,size(Br,1),&
 &ob/sig,1,dble(0),Psample,1)

 !T(:,2) is main diagonal T, T(c:,1) first diagonal below main,T(:,3)
 !first above
 matT=dble(0)
 DO i1=1,size(matT,1)
  matT(i1,2)=dble(1)/alpha(i1)
  IF(i1.GT.1) matT(i1,2)=matT(i1,2)+beta(i1-1)/alpha(i1-1)
  matT(i1,1)=-SQRT(beta(i1))/alpha(i1)
  matT(i1,3)=matT(i1,1)
 END DO
 
 !Appply (B+I)^{-1}
 CALL dptsv(SIZE(matT,1),1,matT(:,2),matT(1:size(matT,1)-1,1),&
 &Psample,size(Psample),status)

 !T(:,2) is main diagonal T, T(:,1) first diagonal below main,T(:,3)
 !first above
 matT=dble(0)
 DO i1=1,size(matT,1)
  matT(i1,2)=dble(1)/alpha(i1)
  IF(i1.GT.1) matT(i1,2)=matT(i1,2)+beta(i1-1)/alpha(i1-1)
  matT(i1,1)=-SQRT(beta(i1))/alpha(i1)
  matT(i1,3)=matT(i1,1)
 END DO

 !Apply (B+I)^{-2}
 CALL dptsv(SIZE(matT,1),1,matT(:,2),matT(1:size(matT,1)-1,1),&
 &Psample,size(Psample),status)

 !Convert back to full observation space
 ALLOCATE(corVec(size(sample))); corVec=0.0
 CALL dgemv('N',size(r,1),size(r,2),dble(1),Br,size(r,1),&
 Psample,1,dble(0),corVec,1)

 WRITE(*,*) 'Norm correction vector:',SQRT(SUM(corVec*corVec,mask.NE.0))

 !Unnormalize
 corVec=corVec*sig 

!---------------------------------------------------------------------
!Calculate correction

 ALLOCATE(cor(s(1))); cor=0.0
 ALLOCATE(maskCor(s(1)))

 !Type 2:
 maskCor=type.EQ.2
 cor1=(SUM(sig**2,maskCor)-SUM(oa*ob,maskCor))/SUM(corVec*ob,maskCor)
 WHERE(maskCor);cor=cor1;END WHERE
 WRITE(*,*) 'Correction type 2:',cor1
 !Type 3:
 maskCor=type.EQ.3
 cor1=(SUM(sig**2,maskCor)-SUM(oa*ob,maskCor))/SUM(corVec*ob,maskCor)
 WHERE(maskCor);cor=cor1;END WHERE
 WRITE(*,*) 'Correction type 3:',cor1
 !Type 4:
 maskCor=type.EQ.4
 cor1=(SUM(sig**2,maskCor)-SUM(oa*ob,maskCor))/SUM(corVec*ob,maskCor)
 WHERE(maskCor);cor=cor1;END WHERE
 WRITE(*,*) 'Correction type 4:',cor1
 !Type 1000-9999:
 maskCor=type.GE.1000.AND.type.LE.9999
 cor1=(SUM(sig**2,maskCor)-SUM(oa*ob,maskCor))/SUM(corVec*ob,maskCor)
 WHERE(maskCor);cor=cor1;END WHERE
 WRITE(*,*) 'Correction type 1000-9999:',cor1

 DEALLOCATE(maskCor)
!------------------------------------------------------------------------
!Read scaling factor
 
 IF(LEN_TRIM(scale_file).GT.0) THEN
  status=nf_open(TRIM(scale_file),nf_nowrite,ncid)
  scaleFac=ncread0d(ncid,'scaleFac')
  status=nf_close(ncid)
 ELSE
  scaleFac=dble(1)
 END IF

 !Scale (no scaling currently)
 x1=x1-0.0*cor*x2

 !Apply Co^{-1/2} to corrected (B+I)^{-1}*d
 x1=x1/sig

 !Create output Co^{-1/2}*(scaleFac*B+I)^{-1}*(Co^{-1/2}*d)
 status=nf_create(TRIM(out_file),nf_clobber,ncid)
 status=nf_def_dim(ncid,'K',nf_unlimited,dimid)
 status=nf_def_var(ncid,TRIM(var_name),nf_double,1,dimid,varid)
 status=nf_close(ncid)

 !Write to output
 status=nf_open(TRIM(out_file),nf_write,ncid)
 CALL ncwrite1d(ncid,TRIM(var_name),[1],[size(x1)],x1)
 status=nf_close(ncid)

 WRITE(*,*) 'Min/max output:',MINVAL(x1),MAXVAL(x1)

 !Now run As,Ad and the rescaled Cb

END PROGRAM ini_scale
