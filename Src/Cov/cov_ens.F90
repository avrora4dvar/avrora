#include "cppdefs_cov.h"
PROGRAM cov_ens

 USE omp_lib
 USE mod_roms
 USE mod_netcdf
 USE mod_interp
 USE mod_constants
 IMPLICIT NONE
#ifdef MPI
 INCLUDE 'mpif.h'
#endif

 INTERFACE 
  SUBROUTINE cov_nf_read(ncid,it,val_out)
   INTEGER,intent(in)::ncid,it
   REAL(8),intent(inout)::val_out(:,:,:,:)
  END SUBROUTINE cov_nf_read
  SUBROUTINE cov_nf_write(ncid,it,val_out)
   INTEGER,intent(in)::ncid,it
   REAL(8),intent(in)::val_out(:,:,:,:)
  END SUBROUTINE cov_nf_write
  SUBROUTINE create_mask_file(mask,nsize,mask_file)
   INTEGER,intent(in)::mask(:,:),nsize(3)
   CHARACTER(len=*),intent(in)::mask_file
  END SUBROUTINE create_mask_file
  SUBROUTINE mask_bnd(xi_cen,eta_cen,d,ncCount,ibnd)
   INTEGER,intent(in)::xi_cen,eta_cen,d(2),ncCount(2)
   INTEGER,intent(out)::ibnd(4,3)
  END SUBROUTINE
 END INTERFACE

 !Variables for parallel code
 INTEGER::myrank,nThreads,irank
#ifdef MPI
 INTEGER::iostatus(mpi_status_size)
#endif

 !Variables for NETCDF and IO
 INTEGER::i0,i1,i2,i3
 INTEGER::status,ncid,fid
 INTEGER::varid(20),dimid(20)
 INTEGER,allocatable::ncids(:)
 CHARACTER(len=1024)::numberStr,memberStr

 !Input parameters                                                     
 CHARACTER(len=1024) :: inFile, grdFile, outFile,pointFile
 CHARACTER(len=1024) :: memberDir,memberFile
 CHARACTER(len=1024) :: command
 INTEGER :: nMembers, tStep(3), nMask, d(2)
 REAL(8)::dl(2)

 !Grid parameter
 REAL(8),allocatable,dimension(:,:)::x,y,h
 INTEGER,allocatable::mask(:,:)
 INTEGER::ncCount(4),ncStart(4),deg2meter(2)

 !Input file
 REAL(8),allocatable,dimension(:,:,:,:)::input

 !Mask
 LOGICAL,allocatable::mMask(:,:,:,:)

 !Mask file
 LOGICAL::flag_exist
 INTEGER::nSeed,now(8),max_seed,ibnd(4,3)
 INTEGER,allocatable::seed(:)
 REAL,allocatable::rNumber(:),pMask(:)
 INTEGER,allocatable::xi_cen(:),eta_cen(:),nCov(:,:,:,:)
 INTEGER,allocatable::nTmp(:,:,:),nTmpMPI(:,:,:)

 !Average
 REAL(8),allocatable,dimension(:,:,:,:)::member,avg

 !Covariance
 REAL(8),allocatable,dimension(:,:,:,:)::cov,tmp,tmpMPI
 REAL(8),allocatable,dimension(:,:,:)::member_input
 REAL(8)::inprod
 
 !Output
 TYPE(roms_type)::par
 TYPE(sigma_param)::spar
 REAL(8),allocatable::zw(:,:,:)
!------------------------------------------------------------------------
! INPUT FILE

#ifdef MPI
 CALL mpi_init(status)
 CALL mpi_comm_size(mpi_comm_world,nThreads,status)
 CALL mpi_comm_rank(mpi_comm_world,myrank,status)
 WRITE(*,*) 'Rank ',myrank,' active'
 IF(myrank.EQ.0) THEN
 WRITE(*,*) 'mpi nThreads:',nThreads
#endif

  !Process input:                                                         
  READ(*,*) !Avrora grid file  
  READ(*,'(A)') grdFile
  READ(*,*) !Input vector (in .nc format) with fields given in section 3       
  READ(*,'(A)') inFile
  READ(*,*) !File to which output will be written                              
  READ(*,'(A)') outFile
  READ(*,*) !Directory containing output of different ensembles runs as
  !subdirectories. Naming must have format Member_###              
  READ(*,'(A)')memberDir
  READ(*,*) !Name of the .nc-file containing model output fields               
  READ(*,'(A)') memberFile
  READ(*,*) !Total number of ensemble members                                 
  READ(*,'(i)') nMembers
  READ(*,*) !Time step to be read from each memberFile                 
  READ(*,*) tStep(1),tStep(2),tStep(3)
  READ(*,*)!Number of localization masks applied to each member              
  READ(*,'(i)') nMask
  READ(*,*) !File containing the centers of the localization masks      
  READ(*,'(A)') pointFile
  READ(*,*) !Half-width localization mask in x,y-direction (in meter)   
  READ(*,*) dl(1),dl(2)
#ifndef MPI
  READ(*,*) !Number of threads                                          
  READ(*,*) nThreads
#endif

  !Display process input:                                                      
  WRITE(*,*) 'Avrora grid file'
  WRITE(*,'(A)') TRIM(grdFile)
  WRITE(*,*) 'Input vector (in .nc format)'
  WRITE(*,'(A)') TRIM(inFile)
  WRITE(*,*) 'File to which output will be written'
  WRITE(*,'(A)') TRIM(outFile)
  WRITE(*,*) 'Directory with members as  subdirectories'
  WRITE(*,'(A)') TRIM(memberDir)
  WRITE(*,*) 'Name of the .nc-file containing model output fields'
  WRITE(*,'(A)') TRIM(memberFile)
  WRITE(*,*) 'Total number of ensemble members'
  WRITE(*,*) nMembers
  WRITE(*,*) 'Time step to be read from each memberFile'
  WRITE(*,*) tStep(1),tStep(2),tStep(3)
  WRITE(*,*) 'Number of localization masks applied to each member'
  WRITE(*,*) nMask
  WRITE(*,*) 'File containing the localization mask'
  WRITE(*,*) TRIM(pointFile)
  WRITE(*,*) 'Half-width localization mask in x-direction (meter)'
  WRITE(*,*) dl(1)
  WRITE(*,*) 'Half-width localization mask in y-direction (meter)'
  WRITE(*,*) dl(2)
  WRITE(*,*)

#ifdef MPI
 END IF !IF(myrank.eq.0)
#endif

!-----------------------------------------------------------------------
! READ GRID

#ifdef MPI
 IF(myrank.EQ.0) THEN
#else
 CALL omp_set_num_threads(nThreads) 
#endif

  WRITE(*,*) 'Reading grid'

  !Open stream
  status=nf_open(TRIM(grdFile),nf_nowrite,ncid)
  varid=0; dimid=0

  !Get size grid
  status=nf_inq_dimid(ncid,'xi_rho',dimid(1))
  status=nf_inq_dimlen(ncid,dimid(1),ncCount(1))
  status=nf_inq_dimid(ncid,'eta_rho',dimid(2))
  status=nf_inq_dimlen(ncid,dimid(2),ncCount(2))
  ncCount(4)=4 !Number of 3D fields in netcdf files
  ncStart=1

  !Create storage grid
  ALLOCATE(x(ncCount(1),ncCount(2)))
  ALLOCATE(y(ncCount(1),ncCount(2)))
  ALLOCATE(mask(ncCount(1),ncCount(2)))
  ALLOCATE(h(ncCount(1),ncCount(2)))

  !Read grid
  status=nf_inq_varid(ncid,'lon_rho',varid(1))
  status=nf_get_vara_double(ncid,varid(1),ncStart(1:2),&
  &ncCount(1:2),x)
  status=nf_inq_varid(ncid,'lat_rho',varid(1))
  status=nf_get_vara_double(ncid,varid(1),ncStart(1:2),&
  &ncCount(1:2),y)
  status=nf_inq_varid(ncid,'mask_rho',varid(1))
  status=nf_get_vara_int(ncid,varid(1),ncStart(1:2),&
  &ncCount(1:2),mask)
  status=nf_inq_varid(ncid,'h',varid(1))
  status=nf_get_vara_double(ncid,varid(1),ncStart(1:2),&
  &ncCount(1:2),h)

  !Close stream
  status=nf_close(ncid)

  !Convert dl to grid points d
  deg2meter(2)=pi*R_earth/dble(180)
  deg2meter(1)=deg2meter(2)*&
  &COS((0.5*MAXVAL(y)+0.5*MINVAL(y))/dble(180)*pi)

  d(1)=INT( dl(1)*dble((size(x,1)-1)*size(x,2))/&
  &SUM(x(2:size(x,1),:)-x(1:size(x,1)-1,:))/deg2meter(1) )
  d(2)=INT( dl(2)*dble(size(y,1)*(size(y,2)-1))/&
  &SUM(y(:,2:size(y,2))-y(:,1:size(y,2)-1))/deg2meter(2) )

  !Write to log
  WRITE(*,*) 'Grid size:',ncCount(1:2)
  WRITE(*,*) 'Half-widths masks are:',d
  WRITE(*,*)

#ifdef MPI
 END IF !IF(myrank.eq.0)
#endif

!------------------------------------------------------------------------
! READ

#ifdef MPI
 IF(myrank.EQ.0) THEN
#endif 

  !Open stream
  status=nf_open(TRIM(inFile),nf_nowrite,ncid)
  varid=0; dimid=0

  !Get size vertical dimension
  status=nf_inq_dimid(ncid,'s_rho',dimid(1))
  IF(status.NE.nf_noerr) THEN
   status=nf_inq_dimid(ncid,'N',dimid(1))
  END IF
  status=nf_inq_dimlen(ncid,dimid(1),ncCount(3))

  !Read input file
  ALLOCATE(input(ncCount(1),ncCount(2),ncCount(3)+1,ncCount(4)))
  CALL cov_nf_read(ncid,tStep(1),input)

  !Close stream
  status=nf_close(ncid)

  WRITE(*,*) 'Input file ',TRIM(inFile),' read'
  WRITE(*,*) 'min/max in file:',minval(input),' ',maxval(input)
  WRITE(*,*)

#ifdef MPI
 END IF !IF(myrank.eq.0)
#endif

!---------------------------------------------------------------------------
! CREATE MASK 

#ifdef MPI
 IF(myrank.EQ.0) THEN
#endif

  !Create mask
  ALLOCATE(mMask(ncCount(1),ncCount(2),ncCount(3)+1,ncCount(4))) 
  mMask=.FALSE.

  !mask zeta
  mMask(:,:,1,1)=mask.EQ.1
  !mask temp
  mMask(:,:,2:,1)=SPREAD(mMask(:,:,1,1),3,ncCount(3))
  !mask salt
  mMask(:,:,2:,2)=mMask(:,:,2:,1)
  !mask u
  mMask(2:,:,2:,3)=mMask(1:ncCount(1)-1,:,2:,1).AND.mMask(2:,:,2:,1)
  !mask v
  mMask(:,2:,2:,4)=mMask(:,1:ncCount(2)-1,2:,1).AND.mMask(:,2:,2:,1)

#ifdef MPI
 END IF !if(myrank.eq.0)
#endif
!--------------------------------------------------------------------------
! CREATE MASK FILE

#ifdef MPI
 IF(myrank.EQ.0) THEN
#endif 

  !Create mask NETCDF file if necessary
  INQUIRE(file=TRIM(pointFile),exist=flag_exist)
  IF(.NOT.flag_exist) THEN
   WRITE(*,*) 'Creating mask file ',TRIM(pointFile)
   CALL create_mask_file(mask,ncCount(1:3),TRIM(pointFile))
  END IF

  !Allocate arrays for mask centers
  ALLOCATE(xi_cen(nMask*nMembers))
  ALLOCATE(eta_cen(nMask*nMembers))
  
  IF(.NOT.flag_exist) THEN
  !Draw center points for mask file
   !Create a seed
   CALL RANDOM_SEED(size=nSeed); ALLOCATE(seed(nSeed))
   CALL RANDOM_SEED(get=seed)
   CALL DATE_AND_TIME(values=now)
   max_seed=MAXVAL(seed)
   DO i0=1,nSeed
    seed(i0)=INT( SQRT( dble(seed(i0))/dble(max_seed)*&
    &dble(PRODUCT(now(7:8)))/dble(60*999) )*dble(max_seed) )+1 
   END DO
   CALL RANDOM_SEED(put=seed)

   !Draw random numbers
   ALLOCATE(rNumber(nMask*nMembers))
   CALL RANDOM_NUMBER(rNumber)

   !Create CDF for field
   ALLOCATE(pMask(ncCount(1)*ncCount(2)))
   pMask=real(RESHAPE(mask,[size(mask)]))
   DO i0=2,size(pMask)
    pMask(i0)=pMask(i0-1)+pMask(i0)
   END DO
   pMask=pMask/MAXVAL(pMask)

   !Find point indices
   DO i0=1,size(rNumber,1)
    i3=MINLOC(pMask,1,pMask.GT.rNumber(i0))
    xi_cen(i0)=MOD(i3-1,ncCount(1))+1
    eta_cen(i0)=INT(dble(i3-1)/dble(ncCount(1)))+1
   END DO

  write(*,*) 'min/max xi_cen',minval(xi_cen),maxval(xi_cen)

   !Write point indices
   status=nf_open(TRIM(pointFile),nf_write,ncid)
   status=nf_inq_varid(ncid,'xi_cen',varid(1))
   status=nf_put_vara_int(ncid,varid(1),[1],[size(xi_cen,1)],xi_cen)
   status=nf_inq_varid(ncid,'eta_cen',varid(2))
   status=nf_put_vara_int(ncid,varid(2),[1],[size(eta_cen,1)],eta_cen)
   status=nf_close(ncid)

   WRITE(*,*) 'xi_cen/eta_cen written to ',TRIM(pointFile)

  ELSE 
   !Read center points from file

   !Open stream
   status=nf_open(TRIM(pointFile),nf_nowrite,ncid)
   varid=0; dimid=0

   !Read
   status=nf_inq_varid(ncid,'xi_cen',varid(1))
   status=nf_get_vara_int(ncid,varid(1),[1],[nMask*nMembers],xi_cen)
   status=nf_inq_varid(ncid,'eta_cen',varid(2))
   status=nf_get_vara_int(ncid,varid(2),[1],[nMask*nMembers],eta_cen)

   !Close stream
   status=nf_close(ncid)
  END IF !IF(.NOT.flag_exist)
  WRITE(*,*) 'min/max xi_cen:',minval(xi_cen),maxval(xi_cen)
  WRITE(*,*) 'min/max eta_cen:',minval(eta_cen),maxval(eta_cen)
 
  !Determine for each grid point in how many masks it is
  ALLOCATE(nCov(ncCount(1),ncCount(2),ncCount(3)+1,ncCount(4))); nCov=0
  IF(flag_exist) THEN
   !Read if from the file
   status=nf_open(TRIM(pointFile),nf_nowrite,ncid)
   status=nf_inq_varid(ncid,'zeta',varid(1))
   status=nf_get_vara_int(ncid,varid(1),ncStart(1:2),ncCount(1:2),&
   &nCov(:,:,1,1))
   status=nf_inq_varid(ncid,'temp',varid(2))
   status=nf_get_vara_int(ncid,varid(2),ncStart(1:3),ncCount(1:3),&
   &nCov(:,:,2:,1))
   status=nf_inq_varid(ncid,'salt',varid(3))
   status=nf_get_vara_int(ncid,varid(3),ncStart(1:3),ncCount(1:3),&
   &nCov(:,:,2:,2))
   status=nf_inq_varid(ncid,'u',varid(4))
   status=nf_get_vara_int(ncid,varid(4),ncStart(1:3),&
   &ncCount(1:3)-[1,0,0],nCov(2:,:,2:,3))
   status=nf_inq_varid(ncid,'v',varid(5))
   status=nf_get_vara_int(ncid,varid(5),ncStart(1:3),&
   &ncCount(1:3)-[0,1,0],nCov(:,2:,2:,4))
   status=nf_close(ncid)
  END IF !IF(flag_exist)

#ifdef MPI
 END IF !IF(myrank.EQ.0)
 CALL mpi_bcast(flag_exist,1,mpi_logical,0,mpi_comm_world,status)
 CALL mpi_bcast(nMask,1,mpi_integer,0,mpi_comm_world,status)
 CALL mpi_bcast(nMembers,1,mpi_integer,0,mpi_comm_world,status)
 CALL mpi_bcast(ncStart,4,mpi_integer,0,mpi_comm_world,status)
 CALL mpi_bcast(ncCount,4,mpi_integer,0,mpi_comm_world,status)
 CALL mpi_bcast(d,2,mpi_integer,0,mpi_comm_world,status) 

 IF(.NOT.ALLOCATED(xi_cen)) ALLOCATE(xi_cen(nMask*nMembers))
 IF(.NOT.ALLOCATED(eta_cen)) ALLOCATE(eta_cen(nMask*nMembers))
 CALL mpi_bcast(xi_cen,size(xi_cen),mpi_integer,0,&
 &mpi_comm_world,status)
 CALL mpi_bcast(eta_cen,size(eta_cen),mpi_integer,0,&
 &mpi_comm_world,status)
#endif

 IF(.NOT.flag_exist) THEN
#ifdef MPI
 ALLOCATE(nTmp(ncCount(1),ncCount(2),3)); nTmp=0
 ALLOCATE(nTmpMPI(ncCount(1),ncCount(2),3)); nTmpMPI=0
#else
!$OMP PARALLEL default(shared) private(i0,nTmp,ibnd)
 IF(omp_get_thread_num().EQ.0) THEN
  WRITE(*,*) 'Entering parallel region with ',omp_get_num_threads(),&
  &' threads'
 END IF
 ALLOCATE(nTmp(ncCount(1),ncCount(2),3)); nTmp=0
#endif

!$OMP DO
 DO i0=1,nMask*nMembers
#ifdef MPI
  IF(MOD(i0,nThreads).NE.myrank) CYCLE
#endif
  CALL mask_bnd(xi_cen(i0),eta_cen(i0),d,ncCount(1:2),ibnd)

  !rho
  nTmp(ibnd(1,1):ibnd(2,1),ibnd(3,1):ibnd(4,1),1)=&
  nTmp(ibnd(1,1):ibnd(2,1),ibnd(3,1):ibnd(4,1),1)+1

  !u
  IF(ibnd(1,2).LE.ibnd(2,2)) THEN
   nTmp(ibnd(1,2):ibnd(2,2),ibnd(3,2):ibnd(4,2),2)=&
   nTmp(ibnd(1,2):ibnd(2,2),ibnd(3,2):ibnd(4,2),2)+1
  END IF

  !v
  IF(ibnd(3,3).LE.ibnd(4,3)) THEN
   nTmp(ibnd(1,3):ibnd(2,3),ibnd(3,3):ibnd(4,3),3)=&
   nTmp(ibnd(1,3):ibnd(2,3),ibnd(3,3):ibnd(4,3),3)+1
  END IF
  
 END DO
!$OMP END DO

#ifdef MPI
 !Collect counts from different threads
 IF(myrank.LT.nThreads-1) CALL mpi_recv(nTmpMPI,size(nTmpMPI),&
 &mpi_integer,myrank+1,0,mpi_comm_world,iostatus,status)
 nTmpMPI=nTmpMPI+nTmp
 IF(myrank.GT.0) CALL mpi_send(nTmpMPI,size(nTmpMPI),&
 &mpi_integer,myrank-1,0,mpi_comm_world,status)

 IF(myrank.EQ.0) THEN
  nCov(:,:,1,1)=nCov(:,:,1,1)+nTmpMPI(:,:,1)
  nCov(:,:,2,3)=nCov(:,:,2,3)+nTmpMPI(:,:,2)
  nCov(:,:,2,4)=nCov(:,:,2,4)+nTmpMPI(:,:,3)

  !Copy masks to all depth layers
  nCov(:,:,2:,1)=SPREAD(nCov(:,:,1,1),3,ncCount(3))
  nCov(:,:,2:,2)=SPREAD(nCov(:,:,1,1),3,ncCount(3))
  nCov(2:,:,2:,3)=SPREAD(nCov(2:,:,2,3),3,ncCount(3))
  nCov(:,2:,2:,4)=SPREAD(nCov(:,2:,2,4),3,ncCount(3))
 END IF

 DEALLOCATE(nTmp,nTmpMPI)
#else
!$OMP CRITICAL
 nCov(:,:,1,1)=nCov(:,:,1,1)+nTmp(:,:,1)
 nCov(:,:,2,3)=nCov(:,:,2,3)+nTmp(:,:,2)
 nCov(:,:,2,4)=nCov(:,:,2,4)+nTmp(:,:,3)
!$OMP END CRITICAL

 DEALLOCATE(nTmp)
!$OMP END PARALLEL
 
 !Copy masks to all depth layers
 nCov(:,:,2:,1)=SPREAD(nCov(:,:,1,1),3,ncCount(3))
 nCov(:,:,2:,2)=SPREAD(nCov(:,:,1,1),3,ncCount(3))
 nCov(2:,:,2:,3)=SPREAD(nCov(2:,:,2,3),3,ncCount(3))
 nCov(:,2:,2:,4)=SPREAD(nCov(:,2:,2,4),3,ncCount(3))

#endif 

#ifdef MPI
  IF(myrank.EQ.0) THEN
#endif
   !Write to file
   status=nf_open(TRIM(pointFile),nf_write,ncid)
   status=nf_inq_varid(ncid,'zeta',varid(1))
   status=nf_put_vara_int(ncid,varid(1),ncStart(1:2),ncCount(1:2),&
   &nCov(:,:,1,1))
   status=nf_inq_varid(ncid,'temp',varid(2))
   status=nf_put_vara_int(ncid,varid(2),ncStart(1:3),ncCount(1:3),&
   &nCov(:,:,2:,1))
   status=nf_inq_varid(ncid,'salt',varid(3))
   status=nf_put_vara_int(ncid,varid(3),ncStart(1:3),ncCount(1:3),&
   &nCov(:,:,2:,2))
   status=nf_inq_varid(ncid,'u',varid(4))
   status=nf_put_vara_int(ncid,varid(4),ncStart(1:3),&
   &ncCount(1:3)-[1,0,0],nCov(2:,:,2:,3))
   status=nf_inq_varid(ncid,'v',varid(5))
   status=nf_put_vara_int(ncid,varid(5),ncStart(1:3),&
   &ncCount(1:3)-[0,1,0],nCov(:,2:,2:,4))
   status=nf_close(ncid)
#ifdef MPI
  END IF !IF(myrank.eq.0)
#endif
 END IF !IF(.NOT.flag_exist)

#ifdef MPI
 IF(myrank.EQ.0) THEN
#endif

  WRITE(*,*) 'min/max coverage:',minval(nCov),maxval(nCov)
  WRITE(*,*)

#ifdef MPI
 END IF !if(myrank.eq.0)
#endif

!--------------------------------------------------------------------------
! OPEN STREAMS TO MEMBERS

  !Allocate array for fields member
  ALLOCATE(member(ncCount(1),ncCount(2),ncCount(3)+1,ncCount(4)));
  member=dble(0)

#ifdef MPI
 IF(myrank.EQ.0) THEN
#endif

  ALLOCATE(ncids(nMembers))

  DO i0=1,nMembers
   WRITE(numberStr,'(I0.3)') i0
   memberStr=Trim( memberDir)//'/Member_'//Trim(numberStr)//'/'//&
   &Trim(memberFile)
   status=nf_open(TRIM(memberStr),nf_nowrite,ncids(i0))
  END DO

#ifdef MPI
 END IF !if(myrank.eq.0)
#endif

!-------------------------------------------------------------------------
! CALCULATE AVERAGE

#ifdef MPI
 IF(myrank.EQ.0) THEN
#endif

  ALLOCATE(avg(ncCount(1),ncCount(2),ncCount(3)+1,ncCount(4)));
  avg=dble(0)

  !Calculate sum members
  DO i0=1,nMembers
   CALL cov_nf_read(ncids(i0),tStep(2),member)
   avg=avg+member
  END DO 

  !Calculate average
  avg=avg/dble(nMembers)
  WHERE(.NOT.mMask)
   avg=dble(0)
  END WHERE

  WRITE(*,*) 'min/max average:',minval(avg),maxval(avg)
  WRITE(*,*)

#ifdef MPI
 END IF !if(myrank.eq.0)
#endif
!----------------------------------------------------------------------------
! CALCULATE COVARIANCE

#ifdef MPI
 IF(myrank.EQ.0) THEN
  !Allocate cov
  ALLOCATE(cov(ncCount(1),ncCount(2),ncCount(3)+1,ncCount(4))); 
  cov=dble(0)
 END IF
 ALLOCATE(member_input(ncCount(1),ncCount(2),ncCount(4)))
 member_input=dble(0)
 ALLOCATE(tmp(ncCount(1),ncCount(2),ncCount(3)+1,ncCount(4)))
 tmp=dble(0)
 ALLOCATE(tmpMPI(ncCount(1),ncCount(2),ncCount(3)+1,ncCount(4)))
 tmpMPI=dble(0)
#else
 ALLOCATE(cov(ncCount(1),ncCount(2),ncCount(3)+1,ncCount(4)))
 cov=dble(0)
 ALLOCATE(member_input(ncCount(1),ncCount(2),ncCount(4)))
 member_input=dble(0) 
!$OMP PARALLEL default(shared) & 
!$OMP& private(i0,i1,i2,ibnd,tmp,inprod,myrank)
 IF(omp_get_thread_num().EQ.0) THEN
  WRITE(*,*) 'Entering parallel region with ',omp_get_num_threads(),&
  &' threads'
 END IF
 ALLOCATE(tmp(ncCount(1),ncCount(2),ncCount(3)+1,ncCount(4)))
 tmp=dble(0)
 myrank=omp_get_thread_num()
#endif

 DO i0=1,nMembers
  IF(myrank.EQ.0) THEN      
   !Read perturbation field of member i0
   CALL cov_nf_read(ncids(i0),tStep(2),member)
   member=member-avg
   WHERE(.NOT.mMask)
    member=dble(0)
   END WHERE
   WRITE(*,*) 'min/max member ',i0,':',minval(member),maxval(member)
   WHERE(nCov.NE.0)
    member=member/SQRT(dble(nCov))
   END WHERE
   member_input=SUM(member*input,3)
  END IF 

#ifdef MPI
  CALL mpi_bcast(member,size(member),mpi_double_precision,&
  &0,mpi_comm_world,status)
  CALL mpi_bcast(member_input,size(member_input),mpi_double_precision,&
  &0,mpi_comm_world,status) 
#else
!$OMP BARRIER 
#endif

!$OMP DO
  DO i1=1,nMask
#ifdef MPI
   IF(MOD(i1,nThreads).NE.myrank) CYCLE
#endif
   i2=nMask*(i0-1)+i1
   CALL mask_bnd(xi_cen(i2),eta_cen(i2),d,ncCount(1:2),iBnd)

   inprod=&
   &SUM(member_input(iBnd(1,1):iBnd(2,1),iBnd(3,1):iBnd(4,1),1:2))+&
   &SUM(member_input(iBnd(1,2):iBnd(2,2),iBnd(3,2):iBnd(4,2),3))+&
   &SUM(member_input(iBnd(1,3):iBnd(2,3),iBnd(3,3):iBnd(4,3),4))

   !zeta
   tmp(iBnd(1,1):iBnd(2,1),iBnd(3,1):iBnd(4,1),1,1)=&
   &tmp(iBnd(1,1):iBnd(2,1),iBnd(3,1):iBnd(4,1),1,1)+&
   &member(iBnd(1,1):iBnd(2,1),iBnd(3,1):iBnd(4,1),1,1)*inprod
   !temp
   tmp(iBnd(1,1):iBnd(2,1),iBnd(3,1):iBnd(4,1),2:,1)=&
   &tmp(iBnd(1,1):iBnd(2,1),iBnd(3,1):iBnd(4,1),2:,1)+&
   &member(iBnd(1,1):iBnd(2,1),iBnd(3,1):iBnd(4,1),2:,1)*inprod
   !salt
   tmp(iBnd(1,1):iBnd(2,1),iBnd(3,1):iBnd(4,1),2:,2)=&
   &tmp(iBnd(1,1):iBnd(2,1),iBnd(3,1):iBnd(4,1),2:,2)+&
   &member(iBnd(1,1):iBnd(2,1),iBnd(3,1):iBnd(4,1),2:,2)*inprod
   !u
   tmp(iBnd(1,2):iBnd(2,2),iBnd(3,2):iBnd(4,2),2:,3)=&
   &tmp(iBnd(1,2):iBnd(2,2),iBnd(3,2):iBnd(4,2),2:,3)+&
   &member(iBnd(1,2):iBnd(2,2),iBnd(3,2):iBnd(4,2),2:,3)*inprod
   !v
   tmp(iBnd(1,3):iBnd(2,3),iBnd(3,3):iBnd(4,3),2:,4)=&
   &tmp(iBnd(1,3):iBnd(2,3),iBnd(3,3):iBnd(4,3),2:,4)+&
   &member(iBnd(1,3):iBnd(2,3),iBnd(3,3):iBnd(4,3),2:,4)*inprod
  END DO !mask
!$OMP END DO

 END DO !members

!Combine data from different threads into a covariance
#ifdef MPI
 IF(myrank.LT.nThreads-1) CALL mpi_recv(tmpMPI,size(tmpMPI),&
 &mpi_double_precision,myrank+1,0,mpi_comm_world,iostatus,status)
 tmpMPI=tmpMPI+tmp
 IF(myrank.GT.0) CALL mpi_send(tmpMPI,size(tmpMPI),&
 &mpi_double_precision,myrank-1,0,mpi_comm_world,status)

 IF(myrank.EQ.0) THEN
  cov=tmpMPI
  cov=dble(nMembers)/dble(nMembers-1)*cov
 END IF

 DEALLOCATE(tmp,tmpMPI,member_input)
#else
!$OMP CRITICAL
 cov=cov+tmp
!$OMP END CRITICAL
 DEALLOCATE(tmp)
!$OMP END PARALLEL
 DEALLOCATE(member_input)
 cov=dble(nMembers)/dble(nMembers-1)*cov
#endif 

#ifdef MPI
 IF(myrank.EQ.0) THEN
#endif
  !Write to log
  WRITE(*,*) 'min/max covariance:',minval(cov),maxval(cov)
  WRITE(*,*)
#ifdef MPI
 END IF
#endif

!---------------------------------------------------------------------------
! CREATE OUTPUT FILE

#ifdef MPI
 IF(myrank.EQ.0) THEN
#endif

  !Get parameters for s-coordinate transformation 
  status=nf_inq_varid(ncids(1),'Vtransform',varid(1))
  status=nf_get_vara_int(ncids(1),varid(1),1,1,par%Vtransform)
  status=nf_inq_varid(ncids(1),'Vstretching',varid(2))
  status=nf_get_vara_int(ncids(1),varid(2),1,1,par%Vstretching)
  status=nf_inq_varid(ncids(1),'theta_s',varid(3))
  status=nf_get_vara_double(ncids(1),varid(3),1,1,par%theta_s)
  status=nf_inq_varid(ncids(1),'Vstretching',varid(4))
  status=nf_get_vara_double(ncids(1),varid(4),1,1,par%theta_b)
  status=nf_inq_varid(ncids(1),'Tcline',varid(5))
  status=nf_get_vara_int(ncids(1),varid(5),1,1,par%Tcline)
  par%grid_size=ncCount(1:3)
  par%nTimes=1
  par%filename=TRIM(outFile)

  !Create output file
  INQUIRE(file=TRIM(outFile),exist=flag_exist)
  IF(flag_exist) THEN
   OPEN(unit=101,file=TRIM(outFile),status='old',iostat=status)
   IF(status.EQ.0) CLOSE(101,status='delete')
  END IF
  CALL roms_create_his_file(TRIM(outFile),par)


  !Write covariance to output file
  status=nf_open(TRIM(outFile),nf_write,ncid)
  status=nf_inq_varid(ncid,'ocean_time',varid(1))
  status=nf_put_vara_double(ncid,varid(1),[1],[1],dble(0))
  CALL cov_nf_write(ncid,1,cov)

  spar%v_transform=par%Vtransform
  spar%v_stretching=par%Vstretching
  spar%theta_s=par%theta_s
  spar%theta_b=par%theta_b
  spar%t_cline=par%Tcline
  spar%n_sigma=ncCount(3)

  ALLOCATE(zw(ncCount(1)-1,ncCount(2),ncCount(3)))
  zw=zweight(spar,h,cov(:,:,1,1),1)
  write(*,*) 'min/max zw_u:',minval(zw),maxval(zw)
  status=nf_inq_varid(ncid,'ubar',varid(2))
  status=nf_put_vara_double(ncid,varid(2),&
  &[ncStart(1),ncStart(2),tStep(3)],&
  &[ncCount(1)-1,ncCount(2),1],SUM(zw*cov(2:,:,2:,3),3) )
  DEALLOCATE(zw)

  ALLOCATE(zw(ncCount(1),ncCount(2)-1,ncCount(3)))
  zw=zweight(spar,h,cov(:,:,1,1),2)
  write(*,*) 'min/max zw_v:',minval(zw),maxval(zw)
  status=nf_inq_varid(ncid,'vbar',varid(3))
  status=nf_put_vara_double(ncid,varid(3),&
  &[ncStart(1),ncStart(2),tStep(3)],&
  &[ncCount(1),ncCount(2)-1,1],SUM(zw*cov(:,2:,2:,4),3) )
  DEALLOCATE(zw)

  status=nf_close(ncid)

  !Write to log
  WRITE(*,*) 'Output file ',TRIM(outFile),' created'
  WRITE(*,*) 
#ifdef MPI
 END IF !if(myrank.eq.0)
#endif 

!---------------------------------------------------------------------------
! FINALIZE

#ifdef MPI
 IF(myrank.EQ.0) THEN
#endif

 DO i0=1,nMembers
  status=nf_close(ncids(i0))
 END DO

 WRITE(*,*) 'cov_ens done'

#ifdef MPI
 END IF !IF(myrank.EQ.0)
 CALL mpi_finalize(status)
#endif

END PROGRAM cov_ens
!-------------------------------------------------------------------------
SUBROUTINE cov_nf_read(ncid,it,val_out)
 IMPLICIT NONE
 INCLUDE 'netcdf.h'
 INTEGER,intent(in)::ncid,it
 REAL(8),intent(inout)::val_out(:,:,:,:)
 INTEGER::varid,ncStart(4),ncCount(4),status

 !Initiate
 val_out=dble(0)

 !count
 ncStart=[1,1,1,it]
 ncCount=[size(val_out,1),size(val_out,2),size(val_out,3)-1,1]

 !Zeta
 status=nf_inq_varid(ncid,'zeta',varid)
 IF(status.NE.nf_noerr) WRITE(*,*) 'Failed to find zeta'
 status=nf_get_vara_double(ncid,varid,&
 &[ncStart(1),ncStart(2),ncStart(4)],&
 &[ncCount(1),ncCount(2),ncCount(4)],&
 val_out(:,:,1,1))
 IF(status.NE.nf_noerr) WRITE(*,*) 'Failed to read zeta'

 !temp
 status=nf_inq_varid(ncid,'temp',varid)
 IF(status.NE.nf_noerr) WRITE(*,*) 'Failed to find temp'
 status=nf_get_vara_double(ncid,varid,ncStart,ncCount,&
 &val_out(:,:,2:,1))
 IF(status.NE.nf_noerr) WRITE(*,*) 'Failed to read temp'

 !salt
 status=nf_inq_varid(ncid,'salt',varid)
 IF(status.NE.nf_noerr) WRITE(*,*) 'Failed to find salt'
 status=nf_get_vara_double(ncid,varid,ncStart,ncCount,&
 &val_out(:,:,2:,2))
 IF(status.NE.nf_noerr) WRITE(*,*) 'Failed to read salt'
 
 !u
 status=nf_inq_varid(ncid,'u',varid)
 IF(status.NE.nf_noerr) WRITE(*,*) 'Failed to find u'
 status=nf_get_vara_double(ncid,varid,ncStart,ncCount-[1,0,0,0],&
 &val_out(2:,:,2:,3))
 IF(status.NE.nf_noerr) WRITE(*,*) 'Failed to read u'
 
 !v
 status=nf_inq_varid(ncid,'v',varid)
 IF(status.NE.nf_noerr) WRITE(*,*) 'Failed to find v'
 status=nf_get_vara_double(ncid,varid,ncStart,ncCount-[0,1,0,0],&
 &val_out(:,2:,2:,4))
 IF(status.NE.nf_noerr) WRITE(*,*) 'Failed to read v'

END SUBROUTINE cov_nf_read
!--------------------------------------------------------------------------
SUBROUTINE cov_nf_write(ncid,it,val_out)
 IMPLICIT NONE
 INCLUDE 'netcdf.h'
 INTEGER,intent(in)::ncid,it
 REAL(8),intent(in)::val_out(:,:,:,:)
 INTEGER::varid,ncStart(4),ncCount(4),status

 !count
 ncStart=[1,1,1,it]
 ncCount=[size(val_out,1),size(val_out,2),size(val_out,3)-1,1]

 !Zeta
 status=nf_inq_varid(ncid,'zeta',varid)
 IF(status.NE.nf_noerr) WRITE(*,*) 'Failed to find zeta'
 status=nf_put_vara_double(ncid,varid,&
 &[ncStart(1),ncStart(2),ncStart(4)],&
 &[ncCount(1),ncCount(2),ncCount(4)],&
 val_out(:,:,1,1))
 IF(status.NE.nf_noerr) WRITE(*,*) 'Failed to write zeta'

 !temp
 status=nf_inq_varid(ncid,'temp',varid)
 IF(status.NE.nf_noerr) WRITE(*,*) 'Failed to find temp'
 status=nf_put_vara_double(ncid,varid,ncStart,ncCount,&
 &val_out(:,:,2:,1))
 IF(status.NE.nf_noerr) WRITE(*,*) 'Failed to write temp'

 !salt
 status=nf_inq_varid(ncid,'salt',varid)
 IF(status.NE.nf_noerr) WRITE(*,*) 'Failed to find salt'
 status=nf_put_vara_double(ncid,varid,ncStart,ncCount,&
 &val_out(:,:,2:,2))
 IF(status.NE.nf_noerr) WRITE(*,*) 'Failed to write salt'
 
 !u
 status=nf_inq_varid(ncid,'u',varid)
 IF(status.NE.nf_noerr) WRITE(*,*) 'Failed to find u'
 status=nf_put_vara_double(ncid,varid,ncStart,ncCount-[1,0,0,0],&
 &val_out(2:,:,2:,3))
 IF(status.NE.nf_noerr) WRITE(*,*) 'Failed to write u'
 
 !v
 status=nf_inq_varid(ncid,'v',varid)
 IF(status.NE.nf_noerr) WRITE(*,*) 'Failed to find v'
 status=nf_put_vara_double(ncid,varid,ncStart,ncCount-[0,1,0,0],&
 &val_out(:,2:,2:,4))
 IF(status.NE.nf_noerr) WRITE(*,*) 'Failed to write v'

END SUBROUTINE cov_nf_write
!---------------------------------------------------------------------------
SUBROUTINE create_mask_file(mask,nsize,mask_file)
!CREATE_MASK_FILE generate the netcdf structure in the file mask      
 IMPLICIT NONE
 INCLUDE 'netcdf.h'
 INTEGER,intent(in)::mask(:,:),nsize(3)
 CHARACTER(len=*),intent(in)::mask_file
 INTEGER::ncid,status,dim_id(12),var_id(9)

 !Create file                                                         
 status=nf_create(TRIM(mask_file),nf_classic_model,ncid)

 !Create dimensions                                                   
 status=nf_def_dim(ncid,'xi_rho',nsize(1),dim_id(1))
 status=nf_def_dim(ncid,'xi_u',nsize(1)-1,dim_id(2))
 status=nf_def_dim(ncid,'xi_v',nsize(1),dim_id(3))
 status=nf_def_dim(ncid,'xi_psi',nsize(1)-2,dim_id(4))
 status=nf_def_dim(ncid,'eta_rho',nsize(2),dim_id(5))
 status=nf_def_dim(ncid,'eta_u',nsize(2),dim_id(6))
 status=nf_def_dim(ncid,'eta_v',nsize(2)-1,dim_id(7))
 status=nf_def_dim(ncid,'eta_psi',nsize(2)-2,dim_id(8))
 status=nf_def_dim(ncid,'N',nsize(3),dim_id(9))
 status=nf_def_dim(ncid,'s_rho',nsize(3),dim_id(10))
 status=nf_def_dim(ncid,'s_w',nsize(3)+1,dim_id(11))
 status=nf_def_dim(ncid,'p',nf_unlimited,dim_id(12))

 !Create variables                                                    
 status=nf_def_var(ncid,'zeta',nf_int,2,&
 &[dim_id(1),dim_id(5),dim_id(12)],var_id(1))
 status=nf_def_var(ncid,'u',nf_int,3,&
&[dim_id(2),dim_id(6),dim_id(10)],var_id(2))
 status=nf_def_var(ncid,'v',nf_int,3,&
 &[dim_id(3),dim_id(7),dim_id(10)],var_id(3))
 status=nf_def_var(ncid,'temp',nf_double,3,&
 &[dim_id(1),dim_id(5),dim_id(10)],var_id(4))
 status=nf_def_var(ncid,'salt',nf_int,3,&
 &[dim_id(1),dim_id(5),dim_id(10)],var_id(5))
 status=nf_def_var(ncid,'xi_cen',nf_int,1,&
 &[dim_id(12)],var_id(6))
 status=nf_def_var(ncid,'eta_cen',nf_int,1,&
 &[dim_id(12)],var_id(7))
 status=nf_def_var(ncid,'w',nf_double,1,&
 &[dim_id(12)],var_id(9))

 !Close netcdf                                                        
 status=nf_close(ncid)
END SUBROUTINE create_mask_file
!---------------------------------------------------------------------------
SUBROUTINE mask_bnd(xi_cen,eta_cen,d,ncCount,ibnd)

 IMPLICIT NONE
 INTEGER,intent(in)::xi_cen,eta_cen,ncCount(2),d(2)
 INTEGER,intent(out)::ibnd(4,3)

  !rho
  ibnd(1,1)=MAXVAL([1,xi_cen-d(1)])
  ibnd(2,1)=MINVAL([ncCount(1),xi_cen+d(1)])
  ibnd(3,1)=MAXVAL([1,eta_cen-d(2)])
  ibnd(4,1)=MINVAL([ncCount(2),eta_cen+d(2)])

  !u
  ibnd(1,2)=MAXVAL([1,xi_cen-d(1)])+1
  ibnd(2,2)=MINVAL([ncCount(1),xi_cen+d(1)])
  ibnd(3,2)=MAXVAL([1,eta_cen-d(2)])
  ibnd(4,2)=MINVAL([ncCount(2),eta_cen+d(2)])

  !v
  ibnd(1,3)=MAXVAL([1,xi_cen-d(1)])
  ibnd(2,3)=MINVAL([ncCount(1),xi_cen+d(1)])
  ibnd(3,3)=MAXVAL([1,eta_cen-d(2)])+1
  ibnd(4,3)=MINVAL([ncCount(2),eta_cen+d(2)]) 

END SUBROUTINE mask_bnd
