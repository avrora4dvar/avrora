MODULE mod_BHM

 COMPLEX,allocatable::bhm_uvTrue(:,:,:) !true complex velocity field
 COMPLEX,allocatable::bhm_uvL(:,:,:) !large scale wind field
 COMPLEX,allocatable::bhm_uvS(:,:,:) !small scale wind field 
 COMPLEX,allocatable::bhm_cL(:,:) !coefficient of large scale velocity field
 COMPLEX,allocatable::bhm_cS(:,:,:) !coefficient small scale velocity field
 COMPLEX,allocatable::bhm_bias(:) !area wide error
 INTEGER,allocatable::bhm_cS_order(:,:) !Debauchie level
 LOGICAL,allocatable::bhm_mask(:,:) !mask with active points 
 REAL,allocatable::bhm_var_cL(:) !variance in coef. large scale vel. field
 REAL::bhm_var_uv(2) !variance representing uncertainty in velocity field
 INTEGER,allocatable::bhm_DB_order(:,:)
 INTEGER::bhm_n_DB_orders
 INTEGER::bhm_seed=1
CONTAINS
!------------------------------------------------------------------------
SUBROUTINE bhm_compose_uvL(uv_mean,uv_eof)
 !Compose large-scale wind fields from EOF and cL-coefficients 
 IMPLICIT NONE
 COMPLEX,intent(in)::uv_eof(:,:,:),uv_mean(:,:)
 INTEGER::fsize(3),it,i3

 fsize=SHAPE(bhm_uvTrue)
 IF(.NOT.ALLOCATED(bhm_uvL)) THEN
  ALLOCATE(bhm_uvL(fsize(1),fsize(2),fsize(3)))
 END IF
 bhm_uvL=(0.0,0.0)

 DO it=1,size(bhm_cL,2)
  DO i3=1,size(bhm_cL,1)
    bhm_uvL(:,:,it)=bhm_uvL(:,:,it)+bhm_cL(i3,it)*uv_eof(:,:,i3)
  END DO
  bhm_uvL(:,:,it)=bhm_uvL(:,:,it)+uv_mean
 END DO
 
END SUBROUTINE 
!--------------------------------------------------------------------------
SUBROUTINE bhm_decompose_uvL(uv_mean,uv_eof)
 !Project wind field on the EOFs
 IMPLICIT NONE
 COMPLEX,intent(in)::uv_eof(:,:,:),uv_mean(:,:)
 COMPLEX,allocatable,dimension(:,:)::val2d1,val2d2,mat,imat
 COMPLEX,allocatable,dimension(:)::v 
 INTEGER::it,i1,i2,i3,status
 
 ALLOCATE(val2d1(size(uv_mean,1),size(uv_mean,2)))
 ALLOCATE(val2d2(size(uv_mean,1),size(uv_mean,2)))
 ALLOCATE(mat(size(uv_eof,3),size(uv_eof,3))); mat=(0.0,0.0)
 ALLOCATE(imat(size(uv_eof,3),size(uv_eof,3))); imat=(0.0,0.0)
 ALLOCATE(v(size(uv_eof,3)))

 !Matrix cross-relations EOFs (should be identity for proper EOFs)
 DO i2=1,size(mat,2)
 DO i1=1,size(mat,1)
  val2d1=uv_eof(:,:,i1)
  val2d2=uv_eof(:,:,i2)
  mat(i1,i2)=SUM(CONJG(val2d1)*val2d2,bhm_mask)
  IF(i1.EQ.i2) imat(i1,i2)=(1.0,0.0)
 END DO 
 END DO

 !Invert mat
 CALL CPOSV('U',size(mat,1),size(imat,2),mat,size(mat,1),&
 &imat,size(imat,1),status)

 DO it=1,size(bhm_cL,2)
  val2d2=bhm_uvL(:,:,it)-uv_mean 
  DO i3=1,size(bhm_cL,1)
   val2d1=uv_eof(:,:,i3)
   v(i3)=SUM(CONJG(val2d1)*val2d2,bhm_mask)
  END DO

  v=MATMUL(imat,v)
  bhm_cL(:,it)=v
 END DO

END SUBROUTINE bhm_decompose_uvL
!-----------------------------------------------------------------------
SUBROUTINE bhm_compose_uvS()
 !Compose small-scale wind field using DB2 wavelets and 
 !cS-coefficients
 IMPLICIT NONE
 INTEGER::fsize(3),step,i3,max_order
 REAL,allocatable,dimension(:,:)::valU,valV
 
 !reset
 fsize=SHAPE(bhm_uvTrue)
 IF(.NOT.ALLOCATED(bhm_uvS)) THEN
   ALLOCATE(bhm_uvS(fsize(1),fsize(2),fsize(3)))
 END IF
 bhm_uvS=(0.0,0.0)

 !maximum order Debauchie decomposition
 CALL bhm_cal_DB_order()
 ALLOCATE(valU(fsize(1),fsize(2)))
 ALLOCATE(valV(fsize(1),fsize(2)))
 max_order=MAXVAL(bhm_DB_order)

 !Convert Debauchie coefficients to wind field
 DO i3=1,fsize(3)
  valU=REAL(bhm_cS(:,:,i3))
  valV=IMAG(bhm_cS(:,:,i3))
  step=2**max_order !step size between scaling coefficients
  DO WHILE(step.GT.1)
   step=step/2 !step size between scaling and transform coefficients
   CALL bhm_IDB2(valU(1:fsize(1):step,1:fsize(2):step))
   CALL bhm_IDB2(valV(1:fsize(1):step,1:fsize(2):step))
  END DO
  bhm_uvS(:,:,i3)=(1,0)*valU+(0,1)*valV
 END DO

 DEALLOCATE(valU,valV)
 
END SUBROUTINE bhm_compose_uvS
!-----------------------------------------------------------------------
SUBROUTINE bhm_decompose_uvS()
 !Find decomposition small-scale wind field in DB2 wavelets
 IMPLICIT NONE
 INTEGER::fsize(3),step,i3
 REAL,allocatable,dimension(:,:)::valU,valV

 !reset
 fsize=SHAPE(bhm_uvTrue)
 IF(.NOT.ALLOCATED(bhm_cS)) THEN
   ALLOCATE(bhm_cS(fsize(1),fsize(2),fsize(3)))
 END IF
 bhm_cS=(0.0,0.0)

 !maximum order Debauchie decomposition
 ALLOCATE(valU(fsize(1),fsize(2)))
 ALLOCATE(valV(fsize(1),fsize(2)))

 !Convert Debauchie coefficients to wind field
 step=1 !step between scaling and transform coefficients
 DO i3=1,fsize(3)
  valU=REAL(bhm_uvS(:,:,i3))
  valV=IMAG(bhm_uvS(:,:,i3))
  DO WHILE(ALL(INT(fsize(1:2)/step).GE.2))
   CALL bhm_DB2(valU(1:fsize(1):step,1:fsize(2):step))
   CALL bhm_DB2(valV(1:fsize(1):step,1:fsize(2):step))
   step=step*2
  END DO
  bhm_cS(:,:,i3)=(1.,0.)*valU+(0.,1.)*valV
 END DO

 DEALLOCATE(valU,valV)
 
END SUBROUTINE bhm_decompose_uvS
!------------------------------------------------------------------------
SUBROUTINE bhm_draw_uv(sample_uv,sig2)
 !Draw a true wind field based on observations in sample_uv, large-scale
 !wind field estimate bhm_uvL and small-scale wind field estimate bhm_uvS
 IMPLICIT NONE
 COMPLEX,intent(in)::sample_uv(:,:,:)
 REAL,intent(in)::sig2(:,:,:)

 INTEGER::fsize(4),i0,i1,i2,i3,it,nSample
 REAL,allocatable::Q2(:,:),Q2mu(:,:),Q2draw(:,:),rand(:)

 fsize(1:3)=SHAPE(bhm_uvTrue)
 fsize(4)=fsize(1)*fsize(2)

 !Allocate inverse variance matrix
 ALLOCATE(Q2(fsize(1),fsize(2)))

 !Generate estimate true wind field for every time step
 ALLOCATE(Q2draw(fsize(1),fsize(2)))
 ALLOCATE(Q2mu(fsize(1),fsize(2)))
 ALLOCATE(rand(2*fsize(4)))
 bhm_uvTrue=(0.0,0.0)
 DO it=1,fsize(3)
       
  CALL bhm_random_gauss(it,rand)

  !draw u_true
  Q2=sig2(:,:,it)
  WHERE(Q2.NE.dble(0)); Q2=dble(1)/Q2; END WHERE
  Q2=Q2+1/bhm_var_uv(1)

  Q2mu=sig2(:,:,it)
  WHERE(Q2mu.NE.dble(0)); Q2mu=dble(1)/Q2mu; END WHERE
  Q2mu=Q2mu*REAL(sample_uv(:,:,it)-bhm_bias(it))+&
  &REAL(bhm_uvL(:,:,it)+bhm_uvS(:,:,it))/bhm_var_uv(1)
 
  Q2draw=SQRT(sig2(:,:,it))
  WHERE(Q2draw.NE.dble(0)); Q2draw=1/Q2draw; END WHERE
  Q2draw=Q2draw*&
  &RESHAPE(rand(fsize(4)+1:2*fsize(4)),[fsize(1),fsize(2)])+&
  &RESHAPE(rand(1:fsize(4)),[fsize(1),fsize(2)])*&
  &1/SQRT(bhm_var_uv(1))

  bhm_uvTrue(:,:,it)=bhm_uvTrue(:,:,it)+CMPLX((Q2mu+Q2draw)/Q2,0.0)
 END DO

 !draw v_true
 DO it=1,fsize(3)

  CALL bhm_random_gauss(it,rand)

  Q2=sig2(:,:,it)
  WHERE(Q2.NE.dble(0)); Q2=dble(1)/Q2; END WHERE
  Q2=Q2+1/bhm_var_uv(2)

  Q2mu=sig2(:,:,it)
  WHERE(Q2mu.NE.dble(0)); Q2mu=dble(1)/Q2mu; END WHERE
  Q2mu=Q2mu*IMAG(sample_uv(:,:,it)-bhm_bias(it))+&
  &IMAG(bhm_uvL(:,:,it)+bhm_uvS(:,:,it))/bhm_var_uv(2)

  Q2draw=SQRT(sig2(:,:,it))
  WHERE(Q2draw.NE.dble(0)); Q2draw=1/Q2draw; END WHERE
  Q2draw=RESHAPE(rand(1:fsize(4)),[fsize(1),fsize(2)])*&
  &1/SQRT(bhm_var_uv(2))+&
  &RESHAPE(rand(fsize(4)+1:2*fsize(4)),[fsize(1),fsize(2)])*&
  &Q2draw

  bhm_uvTrue(:,:,it)=bhm_uvTrue(:,:,it)+CMPLX(0.0,(Q2mu+Q2draw)/Q2)
 END DO 

 !Set masked_out area to 0
 DO i2=1,size(bhm_uvTrue,2)
 DO i1=1,size(bhm_uvTrue,1)
  IF(.NOT.bhm_mask(i1,i2)) bhm_uvTrue(i1,i2,:)=(0.0,0.0)
 END DO
 END DO
 
END SUBROUTINE bhm_draw_uv
!------------------------------------------------------------------------
SUBROUTINE bhm_draw_cEOF(uv_mean,uv_eof)
 !BHM_DRAW_cEOF update the large scale coefficients bhm_cL in the Gibbs
 !sampler by approximating the large scale true wind field with EOFs
 !Input:
 ! uv_eof: COMPLEX array with along dim=3 the different EOFs of the complex
 !         velocity fields
 ! sample_mask: LOGICAL array with the locations (dim=1,2) and time (dim=3)
 !              where the true velocity field is sampled. 
 
 IMPLICIT NONE
 COMPLEX,intent(in)::uv_eof(:,:,:),uv_mean(:,:)
 INTEGER::fsize(3),it,i0,i1,i2,i3,nSample,status
 COMPLEX,allocatable::vec_eof(:,:),vec_true(:),vec_s(:),vec_mean(:)
 COMPLEX,allocatable::Q2(:,:),Q2draw_mu(:)
 INTEGER,allocatable::ipiv(:)
 COMPLEX,allocatable::rand1(:),rand2(:)
 REAL,allocatable::rand0(:)

 !size
 fsize=SHAPE(uv_eof)
 
 !Sample wind field
 nSample=COUNT(bhm_mask)
 ALLOCATE(vec_eof(nSample,fsize(3)))
 ALLOCATE(vec_true(nSample))
 ALLOCATE(vec_s(nSample))
 ALLOCATE(vec_mean(nSample))
 vec_mean=PACK( uv_mean,bhm_mask )

 IF(.NOT.ALLOCATED(bhm_cL))ALLOCATE(bhm_cL(fsize(3),size(bhm_uvTrue,3)))   

 !Generate coefficients for different times
 DO it=1,size(bhm_uvTrue,3)

  vec_true=PACK( bhm_uvtrue(:,:,it),bhm_mask )
  vec_s=PACK(bhm_uvS(:,:,it),bhm_mask )
  DO i3=1,fsize(3)
   vec_eof(:,i3)=PACK(uv_eof(:,:,i3),bhm_mask )
  END DO

  !Generate random numbers
  IF(.NOT.ALLOCATED(rand1))ALLOCATE(rand1(size(vec_eof,1)))
  IF(.NOT.ALLOCATED(rand2))ALLOCATE(rand2(fsize(3)))
  IF(.NOT.ALLOCATED(rand0))ALLOCATE(rand0(2*size(vec_eof,1)+2*fsize(3)))
  CALL bhm_random_gauss(it,rand0)
  rand1=CMPLX(rand0(1:size(vec_eof,1)),&
  &rand0(size(vec_eof,1)+1:2*size(vec_eof,1)))
  rand2=CMPLX(rand0(2*size(vec_eof,1)+1:2*size(vec_eof,1)+fsize(3)),&
  &rand0(2*size(vec_eof,1)+fsize(3)+1:2*size(vec_eof,1)+2*fsize(3)))

  !Q2
  IF(.NOT.ALLOCATED(Q2)) ALLOCATE(Q2(fsize(3),fsize(3)))
  Q2=(0.0,0.0)
  DO i3=1,fsize(3)
   Q2(i3,i3)=1/bhm_var_cL(i3) 
  END DO
  Q2=Q2+MATMUL(TRANSPOSE(CONJG(vec_eof)),vec_eof)/SUM(bhm_var_uv)

  !Q2draw, Q2mu
  IF(.NOT.ALLOCATED(Q2draw_mu)) ALLOCATE(Q2draw_mu(fsize(3)))
  Q2draw_mu=MATMUL(TRANSPOSE(CONJG(vec_eof))/SQRT(SUM(bhm_var_uv)),&
  &rand1)+rand2/SQRT(bhm_var_cL)& 
  &+MATMUL(TRANSPOSE(CONJG(vec_eof))/SUM(bhm_var_uv),vec_true-&
  &vec_mean-vec_s) 

  !Apply Q2^(-1)
  IF(.NOT.ALLOCATED(ipiv)) ALLOCATE(ipiv(fsize(3)))
  CALL cgesv(size(Q2,1),1,Q2,size(Q2,1),ipiv,Q2draw_mu,size(Q2,2),&
  &status)
  IF(status.EQ.0) THEN
   bhm_cL(:,it)=Q2draw_mu
  ELSE
   WRITE(*,*) 'Failed to generate bhm_cL for time step ',it
   STOP
  END IF

 END DO !i

 CALL bhm_compose_uvL(uv_mean,uv_eof)
 
END SUBROUTINE bhm_draw_cEOF
!-------------------------------------------------------------------------
SUBROUTINE bhm_draw_cS(corr1,sig)
 !Generate small-scale DB2 wavelet coefficients such that the small-scale
 !wind field energy spectrum follows a wavelength^2 law
 IMPLICIT NONE
 REAL,intent(in)::sig(:),corr1(:)
 REAL,allocatable::rand(:),sig_2d(:,:),corr1_2d(:,:)
 INTEGER::fsize(3),i0,order,max_order,i3

 !Size
 fsize=SHAPE(bhm_uvTrue)
 IF(.NOT.ALLOCATED(bhm_cS)) ALLOCATE(bhm_cS(fsize(1),fsize(2),fsize(3)))

 !Create AR1 series
 ALLOCATE(rand(2*PRODUCT(fsize)))
 CALL bhm_random_gauss(1,rand)
 bhm_cS=cmplx(RESHAPE(rand(:PRODUCT(fsize)),fsize),&
 &RESHAPE(rand(PRODUCT(fsize)+1:),fsize))

 !maximum Debauchie order
 CALL bhm_cal_DB_order
 max_order=MAXVAL(bhm_DB_order)

 !Assign standard deviation and correlation to each point
 ALLOCATE(corr1_2d(fsize(1),fsize(2))); corr1_2d=0.0
 ALLOCATE(sig_2d(fsize(1),fsize(2))); sig_2d=0.0
 DO i0=1,SIZE(corr1)
  WHERE(bhm_DB_order.EQ.(i0-1))
   corr1_2d=corr1(i0)
   sig_2d=SQRT(&
   &dble(COUNT(bhm_DB_order.EQ.(i0-1)))*sig(i0)**2) !0.5**dble(i0-1)
  END WHERE
 END DO
 WHERE(bhm_DB_order.GT.(SIZE(corr1)-1))
  sig_2d=0.0
  corr1_2d=0.0
 END WHERE  

 !Generate AR1
 bhm_cS(:,:,1)=sig_2d*bhm_cS(:,:,1)
 DO i3=2,fsize(3)
  bhm_cS(:,:,i3)=corr1_2d*bhm_cS(:,:,i3-1)+&
  &SQRT(1.0-corr1_2d**2)*sig_2d*bhm_cS(:,:,i3)
 END DO
 
 write(*,*) 'min/max bhm_cS:',minval(real(bhm_cS)),maxval(real(bhm_cS))&
 ,minval(imag(bhm_cS)),maxval(imag(bhm_cS))

END SUBROUTINE bhm_draw_cS
!--------------------------------------------------------------------------
SUBROUTINE bhm_draw_var_cL(mu,sig2)
 !Draw large-scale wind field EOF coefficients 
 IMPLICIT NONE
 REAL,intent(in)::mu(:),sig2(:)
 REAL::par(2),alpha,beta,rand(1)
 INTEGER::i0

 DO i0=1,size(bhm_cL,1)
  par=bhm_IG_stat(mu(i0),sig2(i0))
  alpha=par(1)+0.5*size(bhm_cL,2)
  beta=par(2)+0.5*SUM( CONJG(bhm_cL(i0,:))*bhm_cL(i0,:) )
  CALL bhm_random_IG(i0,alpha,beta,rand)
  bhm_var_cL(i0)=rand(1)
 END DO

END SUBROUTINE
!--------------------------------------------------------------------------
SUBROUTINE bhm_draw_bias(mu,AR1,sig2)
 !Draw bias value from Gaussian distribution
 IMPLICIT NONE
 COMPLEX,intent(in)::mu
 REAL,intent(in)::sig2(2),AR1
 REAL,allocatable::r(:) 
 INTEGER::s,i1

 IF(.NOT.ALLOCATED(bhm_bias)) THEN
  ALLOCATE(bhm_bias(size(bhm_uvTrue,3)))
 END IF
 s=SIZE(bhm_bias)

 !Draw random 
 ALLOCATE(r(2*s))
 CALL bhm_random_gauss(1,r)
 
 !Set bias
 r(1:s)=SQRT(sig2(1))*r(1:s)
 r(s+1:2*s)=SQRT(sig2(2))*r(s+1:2*s)
 bhm_bias(1)=CMPLX(r(1),r(s+1))
 DO i1=2,s
  bhm_bias(i1)=AR1*bhm_bias(i1-1)+&
  &SQRT(1.0-AR1**2)*CMPLX(r(i1),r(s+i1))
 END DO
 
 !Add mean bias
 bhm_bias=bhm_bias+mu

END SUBROUTINE bhm_draw_bias

!---------------------------------------------------------------------------
SUBROUTINE bhm_draw_var_uv(mu,sig2)
 !Draw single-cell variance wind field
 IMPLICIT NONE
 REAL,intent(in)::mu,sig2
 REAL::rand(1),alpha,beta,par(2)
 INTEGER::nPoints,nT

 !Parameters IG distribution
 par=bhm_IG_stat(mu,sig2)

 nPoints=COUNT(bhm_mask)
 nT=size(bhm_uvTrue,3)

 !Variance u
 alpha=par(1)+0.5*nPoints*nT
 beta=par(2)+0.5*SUM( REAL(bhm_uvTrue-bhm_uvL-bhm_uvS)*&
 &REAL(bhm_uvTrue-bhm_uvL-bhm_uvS) )
 CALL bhm_random_ig(1,alpha,beta,rand)
 bhm_var_uv(1)=rand(1)

 !Variance v
 alpha=par(1)+.5*nPoints*nT
 beta=par(2)+.5*SUM( IMAG(bhm_uvTrue-bhm_uvL-bhm_uvS)*&
 IMAG(bhm_uvTrue-bhm_uvL-bhm_uvS) )
 CALL bhm_random_ig(2,alpha,beta,rand)
 bhm_var_uv(2)=rand(1)

END SUBROUTINE bhm_draw_var_uv
!---------------------------------------------------------------------------
SUBROUTINE bhm_random_uniform(seedFac,rand)
 !BHM_RANDOM_UNIFORM fills the input REAL-array rand with values drawn
 !from a uniform distribution on [0,1). If bhm_random_uniform is called
 !multiple times within a short time period use different strictly positive
 !integers for seedFac. Otherwise use seedFac=1
 
 IMPLICIT NONE
 INTEGER,intent(in)::seedFac
 REAL,intent(inout)::rand(:) 
 INTEGER::nSeed,now(8),i0,i1
 INTEGER,allocatable::seed(:)
 INTEGER*8::seed1
 
 !Generate seed
 CALL RANDOM_SEED(size=nSeed); ALLOCATE(seed(nSeed))            
 CALL RANDOM_SEED()
 CALL RANDOM_SEED(get=seed)
 
 DO i0=1,nSeed
  CALL DATE_AND_TIME(values=now)
  seed1=INT8(seed(i0))
  seed1= MOD(seed1*INT8(MAX(ABS(seedFac),1)),9999999)+1
  seed1= MOD(seed1*INT8(bhm_seed),9999999)+1
  DO i1=1,8
   now(i1)=MAX(ABS(now(i1)),1)
   seed1= MOD( seed1*INT8(now(i1)),9999999)+1
  END DO
  seed(i0)=INT4(seed1)
  END DO
  seed=cshift(seed,now(8))
  CALL RANDOM_SEED(put=seed)
  WRITE(*,*) 'seed:',seed
  
  !Draw random numbers
  CALL RANDOM_NUMBER(rand)

END SUBROUTINE bhm_random_uniform

!---------------------------------------------------------------------------
SUBROUTINE bhm_random_gauss(seedFac,rand) 
!BHM_RANDOM_GAUSS fills the input REAL-array rand with values drawn
!from a standard Normal distribution. If bhm_random_gauss is called
!multiple times within a short time period use different strictly positive
!integers for seedFac. Otherwise use seedFac=1

 IMPLICIT NONE
 INTEGER,intent(in)::seedFac
 REAL,intent(inout)::rand(:)
 REAL,allocatable::rand_uni(:)
 INTEGER::fsize
 REAL::pi

 pi=ACOS(-1.0)
 fsize=size(rand,1)

 !Draw from uniform distribution
 ALLOCATE(rand_uni(2*fsize))
 CALL bhm_random_uniform(seedFac,rand_uni)
 WHERE(rand_uni.EQ.0.0)
  rand_uni=1E-18
 END WHERE
  

 !Convert to Gaussian using Box-Mueller
 rand=SQRT(-2*LOG(rand_uni(1:fsize)))*&
 &COS(2*pi*rand_uni(fsize+1:2*fsize))
 
END SUBROUTINE bhm_random_gauss
!-----------------------------------------------------------------------
SUBROUTINE bhm_random_IG(seedFac,alpha,beta,rand) 
 !BHM_RANDOM_IG draws values distributed as a Inverse Gamma(alpha,beta)
 !distribution. Use INTEGER seedFac when called in rapid progression. 
 !For drawing according to Gamma distribution the accept/rejection
 !algorithm by (Best,1978) is used. 

 IMPLICIT NONE
 REAL,intent(inout)::rand(:)
 INTEGER,intent(in)::seedFac
 REAL,intent(in)::alpha,beta
 INTEGER::i0,i1,fsize
 REAL::b,c,w,x,y,z
 REAL,allocatable::rand_uni(:)

 fsize=size(rand,1)
 b=alpha-1.0
 c=3*alpha-0.75 

 !Draw uniform values
 ALLOCATE(rand_uni(4*fsize))
 CALL bhm_random_uniform(seedFac,rand_uni)
 WHERE(rand_uni.EQ.0.0)
  rand_uni=1.0E-8
 ENd WHERE
 
 i0=0
 DO i1=1,fsize
  DO WHILE(.TRUE.)
   !Selector for random numbers
   IF(i0.EQ.2*fsize+1) THEN
    CALL bhm_random_uniform(seedFac,rand_uni)
    WHERE(rand_uni.EQ.0.0)
     rand_uni=1.0E-8
    END WHERE
    i0=1
   ELSE
    i0=i0+1
   END IF

   !temporary variables
   w=rand_uni(i0)*(1.0-rand_uni(i0))
   y=SQRT(c/w)*(rand_uni(i0)-0.5)
   x=b+y
   
   !Rejection
   IF(x.LE.0.0) CYCLE
   z=64*w**3*rand_uni(2*fsize+i0)**2
   IF(z.LE.(1.0-2*y**2/x)) EXIT
   IF(LOG(z).LE.2*(b*LOG(x/b)-y)) EXIT
  END DO
  rand(i1)=x
 END DO
 
 !Now transform gamma(alpha,1)-distribution to gamma(alpha,beta)
 rand=rand/beta
 !Transform gamma(alpha,beta) to IG(alpha,beta)
 rand=1/rand
 
END SUBROUTINE bhm_random_IG
!------------------------------------------------------------------------
FUNCTION bhm_IG_stat(mu,sig2) RESULT(par)
 IMPLICIT NONE
 REAL,intent(in)::mu,sig2
 REAL::par(2)

 par(1)=mu**2/sig2+2.0
 par(2)=mu*(par(1)-1.0)

END FUNCTION bhm_IG_stat
!------------------------------------------------------------------------
SUBROUTINE bhm_DB2(val2d)
 !Decompose the scaling function 1 level
 IMPLICIT NONE
 REAL,intent(inout)::val2d(:,:)
 REAL,allocatable::tmp1(:,:),tmp2(:,:)
 INTEGER::s0(2),s(2),i1,i2,j1,k1,k2
 REAL,dimension(4)::father1,mother1
 REAL,dimension(4,4)::father2,hor,dia,ver

 !1-dimension farther and mother wavelet
 mother1=[1.-sqrt(3.),-(3.-sqrt(3.)),&
 &3.+sqrt(3.),-(1.+sqrt(3.))]/(4.*sqrt(2.))
 father1=[1.+sqrt(3.),3.+sqrt(3.),3.-sqrt(3.),1.-sqrt(3.)]/(4.*sqrt(2.))
 father1=SQRT(0.5)*father1
 mother1=SQRT(0.5)*mother1 

 !2-dimensional farther wavelet and mother wavelets (superfluous)
 father2=SPREAD(father1,2,4)*SPREAD(father1,1,4)
 hor=SPREAD(father1,2,4)*SPREAD(mother1,1,4)
 ver=SPREAD(mother1,2,4)*SPREAD(father1,1,4)
 dia=SPREAD(mother1,2,4)*SPREAD(mother1,1,4)

 !Size grids
 s0=SHAPE(val2d)
 s=s0-MOD(s0,2)+2
 
 !Allocate temporary grid
 ALLOCATE(tmp1(s(1),s(2))); tmp1=dble(0)
 ALLOCATE(tmp2(s(1),s(2))); tmp2=dble(0)
 tmp1(1:s0(1),1:s0(2))=val2d

 !Decompose along dimension 1
 DO i1=1,s(1)-2,2
  k1=i1
  DO j1=1,4
   IF(k1.GT.s(1)-2) k1=1
   tmp2(i1,:)=tmp2(i1,:)+father1(j1)*tmp1(k1,:) !father2,hor
   tmp2(i1+1,:)=tmp2(i1+1,:)+mother1(j1)*tmp1(k1,:) !ver,dia
   k1=k1+1
  END DO !j1
 END DO !i1 

 !Odd row (father only)
 tmp2(s(1)-1,:)=tmp1(s(1)-1,:)

 !Decompose along dimension 2
 tmp1=dble(0)
 DO i2=1,s(2)-2,2
  k2=i2
  DO j1=1,4
   IF(k2.GT.s(2)-2) k2=1
   tmp1(:,i2)=tmp1(:,i2)+father1(j1)*tmp2(:,k2) !i1 odd: father2, even:ver
   tmp1(:,i2+1)=tmp1(:,i2+1)+mother1(j1)*tmp2(:,k2) !i1 odd: hor, even:dia
   k2=k2+1
  END DO !j1
 END DO !i2

 !Odd column (father only)
 tmp1(:,s(2)-1)=tmp2(:,s(2)-1)
 
 !Output
 val2d=tmp1(1:s0(1),1:s0(2))

END SUBROUTINE bhm_DB2
!------------------------------------------------------------------------
SUBROUTINE bhm_IDB2(val2d)
 !Compose a scaling function at level orderIn
 IMPLICIT NONE
 REAL,intent(inout)::val2d(:,:)
 REAL,allocatable::tmp1(:,:),tmp2(:,:)
 INTEGER::s0(2),s(2),i1,i2,j1,k1,k2
 REAL,dimension(4)::father1,mother1
 REAL,dimension(4,4)::father2,hor,dia,ver

 !1-dimension farther and mother wavelet
 mother1=[1.-sqrt(3.),-(3.-sqrt(3.)),&
 &3.+sqrt(3.),-(1.+sqrt(3.))]/(4.*sqrt(2.))
 father1=[1.+sqrt(3.),3.+sqrt(3.),3.-sqrt(3.),1.-sqrt(3.)]/(4.*sqrt(2.))
 father1=SQRT(2.0)*father1
 mother1=SQRT(2.0)*mother1
 
 !2-dimensional farther wavelet and mother wavelets (superfluous)
 father2=SPREAD(father1,2,4)*SPREAD(father1,1,4)
 hor=SPREAD(father1,2,4)*SPREAD(mother1,1,4)
 ver=SPREAD(mother1,2,4)*SPREAD(father1,1,4)
 dia=SPREAD(mother1,2,4)*SPREAD(mother1,1,4)

 !Size grids
 s0=SHAPE(val2d)
 s=s0-MOD(s0,2)+2
 
 !Allocate temporary grid
 ALLOCATE(tmp1(s(1),s(2))); tmp1=dble(0)
 ALLOCATE(tmp2(s(1),s(2))); tmp2=dble(0)
 tmp1(1:s0(1),1:s0(2))=val2d


 !Odd column (father only)
 tmp2(:,s(2)-1)=tmp1(:,s(2)-1)

 !Compose along dimension 2
 DO i2=1,s(2)-2,2
  k2=i2
  DO j1=1,4
   IF(k2.GT.s(2)-2) k2=1
   tmp2(:,k2)=tmp2(:,k2)+father1(j1)*tmp1(:,i2)&
   &                    +mother1(j1)*tmp1(:,i2+1)
   k2=k2+1
  END DO !j1
 END DO !i2
 tmp1=dble(0)

  !Odd row (father only)
 tmp1(s(1)-1,:)=tmp2(s(1)-1,:)

 !Compose along dimension 1
 DO i1=1,s(1)-2,2
  k1=i1
  DO j1=1,4
   IF(k1.GT.s(1)-2) k1=1
   tmp1(k1,:)=tmp1(k1,:)+father1(j1)*tmp2(i1,:)&
                        +mother1(j1)*tmp2(i1+1,:)
   k1=k1+1
  END DO !j1
 END DO !i1 

 !Out
 val2d=tmp1(1:s0(1),1:s0(2))

END SUBROUTINE bhm_IDB2
!--------------------------------------------------------------------------
SUBROUTINE bhm_cal_DB_order()
 !Calculate for each grid point the highest order where the grid point
 !is still a scaling fuction. 
 IMPLICIT NONE
 INTEGER::order,fsize(3),step
 INTEGER::s(2),s0(2)

 s0=[size(bhm_uvTrue,1),size(bhm_uvTrue,2)]

 IF(.NOT.ALLOCATED(bhm_DB_order)) THEN
  ALLOCATE(bhm_DB_order(s0(1),s0(2)))
 bhm_DB_order=0

 order=1; step=2
 DO WHILE(MINVAL(s0)/step.GE.1)
  s=FLOOR(s0/dble(step))*step
  bhm_DB_order(1:s0(1):step,1:s0(2):step)=order
  order=order+1
  step=step*2
 END DO
 WRITE(*,*) 'max DB order:',MAXVAL(bhm_DB_order)

 END IF

END SUBROUTINE bhm_cal_DB_order
!----------------------------------------------------------------------------
END MODULE mod_BHM



