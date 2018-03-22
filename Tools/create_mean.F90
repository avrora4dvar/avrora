PROGRAM create_mean

 USE mod_netcdf
 USE mod_roms
 IMPLICIT NONE

!-------------------------------------------------------------------------
!Variables

 !Input
 CHARACTER(len=1024)::grd_file,ens_dir,member_name,per_name,mean_file
 INTEGER::time_in,nMembers
 REAL(8)::time_out(1),default_val(5)

 !Grid
 INTEGER::s(8)
 LOGICAL,allocatable,dimension(:,:)::mask_r,mask_v,mask_u

 !Member files
 CHARACTER(len=1024)::member_file
 TYPE(roms_type)::fpar
 REAL(8),allocatable,dimension(:,:,:,:)::temp,salt,u,v
 REAL(8),allocatable,dimension(:,:,:)::zeta,ubar,vbar

 !Mean
 REAL(8),allocatable,dimension(:,:,:,:)::tempM,saltM,uM,vM
 REAL(8),allocatable,dimension(:,:,:)::zetaM,ubarM,vbarM

 !netcdf
 INTEGER::status,ncid

 !Counter
 INTEGER::iMember 


!------------------------------------------------------------------------
!Input

 READ(*,*) !Grid file
 READ(*,'(A)') grd_file
 READ(*,*) !Ensemble directory
 READ(*,'(A)') ens_dir
 READ(*,*) !Number of ensemble members
 READ(*,*) nMembers
 READ(*,*) !Name of file in Member_XXX directory
 READ(*,'(A)') member_name
 READ(*,*) !Time field to be read from member_name file
 READ(*,*) time_in
 READ(*,*) !output file with ensemble mean
 READ(*,'(A)') mean_file
 READ(*,*) !Default values for temp,salt,zeta,u/ubar,v/vbar outside mask
 READ(*,*) default_val(1),default_val(2),default_val(3),&
 &default_val(4),default_val(5)
 READ(*,*) !time in output file in seconds since dateref
 READ(*,*) time_out
 READ(*,*) !output file with perturbation from mean (opt)
 READ(*,'(A)') per_name

!-----------------------------------------------------------------------
!Read grid file

 status=nf_open(TRIM(grd_file),nf_nowrite,ncid)
 CALL ncsize(ncid,'mask_rho',s)

 !Read masks
 ALLOCATE(mask_r(s(1),s(2)))
 mask_r=ncread2d(ncid,'mask_rho',[1,1],[s(1),s(2)]).EQ.dble(0)
 ALLOCATE(mask_u(s(1)-1,s(2)))
 mask_u=ncread2d(ncid,'mask_u',[1,1],[s(1)-1,s(2)]).EQ.dble(0)
 ALLOCATE(mask_v(s(1),s(2)-1))
 mask_v=ncread2d(ncid,'mask_v',[1,1],[s(1),s(2)-1]).EQ.dble(0)
 
 status=nf_close(ncid)

!------------------------------------------------------------------------
!Read files

 DO iMember=1,nMembers
  WRITE(member_file,'(A,A,I0.3,A,A)') TRIM(ens_dir),'/Member_',&
  &iMember,'/',TRIM(member_name)
  WRITE(*,*) 'Reading from ',TRIM(member_file)

  IF(iMember.EQ.1) THEN
   !Read layout file
   CALL roms_get_type(TRIM(member_file),fpar)
   s(3)=fpar%grid_size(3)   

   !Allocate field storage
   ALLOCATE(temp(s(1),s(2),s(3),1))
   ALLOCATE(salt(s(1),s(2),s(3),1))
   ALLOCATE(u(s(1)-1,s(2),s(3),1))
   ALLOCATE(v(s(1),s(2)-1,s(3),1))
   ALLOCATE(zeta(s(1),s(2),1))
   ALLOCATE(ubar(s(1)-1,s(2),1))
   ALLOCATE(vbar(s(1),s(2)-1,1))
   
   !Allocate storage for means
   ALLOCATE(tempM(s(1),s(2),s(3),1)); tempM=dble(0)
   ALLOCATE(saltM(s(1),s(2),s(3),1)); saltM=dble(0)
   ALLOCATE(uM(s(1)-1,s(2),s(3),1)); uM=dble(0)
   ALLOCATE(vM(s(1),s(2)-1,s(3),1)); vM=dble(0)
   ALLOCATE(zetaM(s(1),s(2),1)); zetaM=dble(0)
   ALLOCATE(ubarM(s(1)-1,s(2),1)); ubarM=dble(0)
   ALLOCATE(vbarM(s(1),s(2)-1,1)); vbarM=dble(0)
  END IF !iMember.EQ.1

  !Open netcdf stream 
  status=nf_open(TRIM(member_file),nf_nowrite,ncid)

  !Read member
  temp=ncread4d(ncid,'temp',[1,1,1,time_in],[s(1),s(2),s(3),1])
  salt=ncread4d(ncid,'salt',[1,1,1,time_in],[s(1),s(2),s(3),1])
  u=ncread4d(ncid,'u',[1,1,1,time_in],[s(1)-1,s(2),s(3),1])
  v=ncread4d(ncid,'v',[1,1,1,time_in],[s(1),s(2)-1,s(3),1]) 
  zeta=ncread3d(ncid,'zeta',[1,1,time_in],[s(1),s(2),1])
  ubar=ncread3d(ncid,'ubar',[1,1,time_in],[s(1)-1,s(2),1])
  vbar=ncread3d(ncid,'vbar',[1,1,time_in],[s(1),s(2)-1,1])

  !Add to mean
  tempM=tempM+temp
  saltM=saltM+salt
  uM=uM+u
  vM=vM+v
  zetaM=zetaM+zeta
  ubarM=ubarM+ubar
  vbarM=vbarM+vbar

  !Close netcdf stream
  status=nf_close(ncid)
 END DO !iMember

!-----------------------------------------------------------------------------
 !Calculate mean
 WRITE(*,*)
 WRITE(*,*) 'Calculating mean'

 tempM=tempM/dble(nMembers)
 saltM=saltM/dble(nMembers)
 WHERE(RESHAPE(SPREAD(mask_r,3,s(3)),[s(1),s(2),s(3),1]))
  tempM=default_val(1)
  saltM=default_val(2)
 END WHERE
  
 zetaM=zetaM/dble(nMembers)
 WHERE(RESHAPE(mask_r,[s(1),s(2),1]))
  zetaM=default_val(3)
 END WHERE
  
 uM=uM/dble(nMembers)
 WHERE(RESHAPE(SPREAD(mask_u,3,s(3)),[s(1)-1,s(2),s(3),1]))
  uM=default_val(4)
 END WHERE
 
 ubarM=ubarM/dble(nMembers)
 WHERE(RESHAPE(mask_u,[s(1)-1,s(2),1]))
  ubarM=default_val(4)
 END WHERE
 
 vM=vM/dble(nMembers)
 WHERE(RESHAPE(SPREAD(mask_v,3,s(3)),[s(1),s(2)-1,s(3),1]))
  vM=default_val(5)
 END WHERE

 vbarM=vbarM/dble(nMembers)
 WHERE(RESHAPE(mask_v,[s(1),s(2)-1,1]))
  vbarM=default_val(5)
 END WHERE

!---------------------------------------------------------------------------
!Create output

 !Remove if exist
 OPEN(unit=999,file=TRIM(mean_file),iostat=status,status='old')
 IF(status.EQ.0) CLOSE(999,status='delete')

 !Create output file
 fpar%nTimes=1
 CALL roms_create_his_file(TRIM(mean_file),fpar)
 WRITE(*,*) 'Writing mean to ',TRIM(mean_file)
 WRITE(*,*) 

 !Open netcdf stream
 status=nf_open(TRIM(mean_file),nf_write,ncid)

 !Write mean to output
 CALL ncwrite4d(ncid,'temp',[1,1,1,1],[s(1),s(2),s(3),1],tempM)
 CALL ncwrite4d(ncid,'salt',[1,1,1,1],[s(1),s(2),s(3),1],saltM)
 CALL ncwrite4d(ncid,'u',[1,1,1,1],[s(1)-1,s(2),s(3),1],uM)
 CALL ncwrite4d(ncid,'v',[1,1,1,1],[s(1),s(2)-1,s(3),1],vM)
 CALL ncwrite3d(ncid,'zeta',[1,1,1],[s(1),s(2),1],zetaM)
 CALL ncwrite3d(ncid,'ubar',[1,1,1],[s(1)-1,s(2),1],ubarM)
 CALL ncwrite3d(ncid,'vbar',[1,1,1],[s(1),s(2)-1,1],vbarM)

 !Write time to output
 CALL ncwrite1d(ncid,'ocean_time',[1],[1],time_out)

 !Close netcdf stream
 status=nf_close(ncid)

!---------------------------------------------------------------------------
!Write perturbation files
 
 IF(LEN_TRIM(per_name).NE.0) THEN

  DO iMember=1,nMembers
   WRITE(member_file,'(A,A,I0.3,A,A)') TRIM(ens_dir),'/Member_',&
   &iMember,'/',TRIM(member_name)
   WRITE(*,*) 'Reading ',TRIM(member_file)

   !Open netcdf stream 
   status=nf_open(TRIM(member_file),nf_nowrite,ncid)

   !Read member
   temp=ncread4d(ncid,'temp',[1,1,1,time_in],[s(1),s(2),s(3),1])
   salt=ncread4d(ncid,'salt',[1,1,1,time_in],[s(1),s(2),s(3),1])
   u=ncread4d(ncid,'u',[1,1,1,time_in],[s(1)-1,s(2),s(3),1])
   v=ncread4d(ncid,'v',[1,1,1,time_in],[s(1),s(2)-1,s(3),1]) 
   zeta=ncread3d(ncid,'zeta',[1,1,time_in],[s(1),s(2),1])
   ubar=ncread3d(ncid,'ubar',[1,1,time_in],[s(1)-1,s(2),1])
   vbar=ncread3d(ncid,'vbar',[1,1,time_in],[s(1),s(2)-1,1])

   !Substract mean
   temp=temp-tempM
   salt=salt-saltM
   u=u-uM
   v=v-vM
   zeta=zeta-zetaM
   ubar=ubar-ubarM
   vbar=vbar-vbarM

   !mask out
   WHERE(RESHAPE(SPREAD(mask_r,3,s(3)),[s(1),s(2),s(3),1]))
    temp=dble(0)
    salt=dble(0)
   END WHERE
   WHERE(RESHAPE(SPREAD(mask_u,3,s(3)),[s(1)-1,s(2),s(3),1]))
    u=dble(0)
   END WHERE
   WHERE(RESHAPE(SPREAD(mask_v,3,s(3)),[s(1),s(2)-1,s(3),1]))
    v=dble(0)
   END WHERE
   WHERE(RESHAPE(mask_r,[s(1),s(2),1]))
    zeta=dble(0)
   END WHERE
   WHERE(RESHAPE(mask_u,[s(1)-1,s(2),1]))
    ubar=dble(0)
   END WHERE
   WHERE(RESHAPE(mask_v,[s(1),s(2)-1,1]))
    vbar=dble(0)
   END WHERE

   !Close netcdf stream to input
   status=nf_close(ncid)

   WRITE(member_file,'(A,A,I0.3,A,A)') TRIM(ens_dir),'/Member_',&
   &iMember,'/',TRIM(per_name)
   WRITE(*,*) 'Creating ',TRIM(member_file)

   !Remove if exist
   OPEN(unit=999,file=TRIM(member_file),iostat=status,status='old')
   IF(status.EQ.0) CLOSE(999,status='delete')

   !Create output file
   fpar%nTimes=1
   CALL roms_create_his_file(TRIM(member_file),fpar)

   !Open stream to output
   status=nf_open(TRIM(member_file),nf_write,ncid)

   !Write mean to output
   CALL ncwrite4d(ncid,'temp',[1,1,1,1],[s(1),s(2),s(3),1],temp)
   CALL ncwrite4d(ncid,'salt',[1,1,1,1],[s(1),s(2),s(3),1],salt)
   CALL ncwrite4d(ncid,'u',[1,1,1,1],[s(1)-1,s(2),s(3),1],u)
   CALL ncwrite4d(ncid,'v',[1,1,1,1],[s(1),s(2)-1,s(3),1],v)
   CALL ncwrite3d(ncid,'zeta',[1,1,1],[s(1),s(2),1],zeta)
   CALL ncwrite3d(ncid,'ubar',[1,1,1],[s(1)-1,s(2),1],ubar)
   CALL ncwrite3d(ncid,'vbar',[1,1,1],[s(1),s(2)-1,1],vbar)

   !Write time to output
   CALL ncwrite1d(ncid,'ocean_time',[1],[1],time_out)

   !Close netcdf stream to output
   status=nf_close(ncid)
   
  END DO !iMember 

 END IF !per_name(1).NE.' '

 
END PROGRAM create_mean
