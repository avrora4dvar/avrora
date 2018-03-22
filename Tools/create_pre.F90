PROGRAM create_pre

 USE mod_netcdf 
 USE mod_BHM
 IMPLICIT NONE

 !Counters
 INTEGER::i0,i1,i2,i3,ncid(999),s(8),status,varid,dimid(4)
 LOGICAL::flag_file

 !Observation
 INTEGER,allocatable::type(:)
 REAL(8),allocatable::obs(:),sig(:)

 !Input
 INTEGER::n_members
 CHARACTER(len=1024)::ens_dir,obs_file,member_file,iter_dir  
 CHARACTER(len=1024)::r_file,r_name,x_file,x_name
 REAL(8)::sig_max

 !Draw
 REAL,allocatable::draw1(:)
 REAL(8),allocatable::draw2(:,:)
 
!--------------------------------------------------------------------

 READ(*,*) !Observation file
 READ(*,'(A)') obs_file
 READ(*,*) !Ensemble dir
 READ(*,'(A)') iter_dir
 READ(*,'(A)') ens_dir
 READ(*,*) !Read number members
 READ(*,*) n_members
 READ(*,*) !Name output file and variable name
 READ(*,'(A)') r_file
 READ(*,'(A)') r_name
 READ(*,*) !Name output file and variable name
 READ(*,'(A)') x_file
 READ(*,'(A)') x_name 
 READ(*,*) !Read maximum obs
 READ(*,*) sig_max

 WRITE(*,*) 'Observation file:',TRIM(obs_file)
 WRITE(*,*) 'Iter directory:',TRIM(iter_dir)
 WRITE(*,*) 'Ens directory:',TRIM(ens_dir)
 WRITE(*,*) 'Number of threads:',n_members
 WRITE(*,*) 'r file:',TRIM(r_file)
 WRITE(*,*) 'x file:',TRIM(x_file)

 ncid=0
!--------------------------------------------------------------------
!Read observations


 WRITE(*,*) 'Reading observation file:',TRIM(obs_file)
 status=nf_open(TRIM(obs_file),nf_nowrite,ncid(1))
 CALL ncsize(ncid(1),'obs',s)

 ALLOCATE(obs(s(1)))
 obs=ncread1d(ncid(1),'obs',[1],[s(1)])
 ALLOCATE(sig(s(1)))
 sig=ncread1d(ncid(1),'sig_d',[1],[s(1)])
 ALLOCATE(type(s(1)))
 type=INT(ncread1d(ncid(1),'type',[1],[s(1)]))

 status=nf_close(ncid(1))

!--------------------------------------------------------------------
!Read sample

 WRITE(*,*) 'Reading sample ',TRIM(r_file)

 WRITE(member_file,'(A,A,A)') TRIM(iter_dir),'/',&
 &TRIM(r_file)
 status=nf_open(TRIM(member_file),nf_nowrite,ncid(1))
 obs=obs-ncread1d(ncid(1),TRIM(r_name),[1],[size(obs,1)])
 status=nf_close(ncid(1)); ncid(1)=0

 WRITE(*,*) 'Min/max innoviation:',MINVAL(obs),MAXVAL(obs)
 WHERE(ABS(obs).GT.dble(sig_max)*sig.AND..NOT.&
 &(type.GE.100.AND.type.LT.200) )
  sig=ABS(obs)/dble(sig_max)
 END WHERE
 WRITE(*,*) 'Min/max sig:',MINVAL(sig),MAXVAL(sig)


!--------------------------------------------------------------------
!Draw

 ALLOCATE(draw1(s(1)*n_members))
 ALLOCATE(draw2(s(1),n_members))
 CALL bhm_random_gauss(1,draw1)
 draw2=RESHAPE(dble(draw1),[size(draw2,1),size(draw2,2)])

 DO i0=1,size(draw2,2)
  WRITE(*,*) 'Std ',i0,':',SQRT(SUM(draw2(:,i0)**2)/dble(size(draw2,1)))
 END DO
 !draw2=draw2/SQRT(dble(n_members))
 DEALLOCATE(draw1)

 !Set salinity to zero
 DO i1=1,size(draw2,1)
  IF(type(i1).GE.100.AND.type(i1).LT.200) THEN
   draw2(i1,:)=0.0
  END IF
 END DO

 !Apply sig
 DO i0=1,size(draw2,2)
  draw2(:,i0)=draw2(:,i0)/sig
 END DO

!----------------------------------------------------------------------

 !Open stream to members
 DO i0=1,n_members
  WRITE(member_file,'(A,A,I0.3,A,A)') TRIM(ens_dir),'/Member_',i0,'/',&
  &TRIM(x_file)
  INQUIRE(file=TRIM(member_file),exist=flag_file)
  WRITE(*,*) 'Creating ',TRIM(member_file),flag_file

  IF(.NOT.flag_file) THEN
   status=nf_create(TRIM(member_file),nf_classic_model,&
   &ncid(i0))
   status=nf_def_dim(ncid(i0),TRIM(x_name),s(1),dimid(1))
   status=nf_def_var(ncid(i0),TRIM(x_name),nf_double,1,dimid(1),varid)
   status=nf_close(ncid(i0)); ncid(i0)=0
  END IF

  !Write output
  status=nf_open(TRIM(member_file),nf_write,ncid(i0))
  WRITE(*,*) 'Write min/max:',minval(draw2(:,i0)),maxval(draw2(:,i0))
  CALL ncwrite1d( ncid(i0),TRIM(x_name),[1],[s(1)],&
  &draw2(:,i0) )
  status=nf_close(ncid(i0)); ncid(i0)=0

 END DO

 WRITE(*,*) 'create_pre DONE'
 
END PROGRAM create_pre

