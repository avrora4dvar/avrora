PROGRAM ad_identity
 USE mod_netcdf
 USE mod_roms

 !input
 CHARACTER(len=1024)::in_file,outi_file,outf_file,grd_file
 INTEGER::ntimes,nhis
 REAL(8)::dt

 !file
 TYPE(roms_type)::fpar
 INTEGER::gs(3),nt

 !content input file
 REAL(8),allocatable,dimension(:)::t
 REAL(8),allocatable,dimension(:,:,:)::zeta
 REAL(8),allocatable,dimension(:,:,:,:)::temp,salt,u,v

 !netcdf
 INTEGER::status,ncid,i0,i1,i2,i3,i4

!---------------------------------------------------------------------------
!Read input

 READ(*,*) !number times 
 READ(*,*) ntimes
 READ(*,*) !B/c time step DT:
 READ(*,*) dt
 READ(*,*) !Nfast
 READ(*,*) 
 READ(*,*) !Output to history file each nhis times
 READ(*,*) nhis
 READ(*,*) !background eddy viscosity
 READ(*,*)
 READ(*,*) !background eddy diffusivity
 READ(*,*) 
 READ(*,*) !bottom friction coefficient
 READ(*,*)
 READ(*,*) !horizontal viscosity
 READ(*,*)
 READ(*,*) !horizontal diffusion
 READ(*,*) 
 READ(*,*) !Parameters of linear EOS
 READ(*,*)
 READ(*,*) !grid file
 READ(*,'(A)') grd_file
 READ(*,*) !Initial condition file (output)
 READ(*,'(A)') outi_file
 READ(*,*) !Initial forcing file (output)
 READ(*,*) outf_file
 READ(*,*) !background file
 READ(*,*) 
 READ(*,*) !input file (input)
 READ(*,'(A)') in_file
 READ(*,*) !tiles
 READ(*,*) 

!---------------------------------------------------------------------------
!Read input

 WRITE(*,*) 'reading input file: ',TRIM(in_file)

 !Read vertical structure
 CALL roms_get_type(TRIM(in_file),fpar)
 gs=fpar%grid_size
 WRITE(*,*) 'grid size: ',gs
 
 !open file 
 status=nf_open(TRIM(in_file),nf_nowrite,ncid)

 !zeta
 ALLOCATE(zeta(gs(1),gs(2),fpar%nTimes)); zeta=dble(0)
 zeta=ncread3d(ncid,'zeta',[1,1,1],[gs(1),gs(2),fpar%nTimes])
  
 !temp
 ALLOCATE(temp(gs(1),gs(2),gs(3),fpar%nTimes)); temp=dble(0)
 temp=ncread4d(ncid,'temp',[1,1,1,1],[gs(1),gs(2),gs(3),fpar%nTimes])
 
 !salt
 ALLOCATE(salt(gs(1),gs(2),gs(3),fpar%nTimes)); salt=dble(0)
 salt=ncread4d(ncid,'salt',[1,1,1,1],[gs(1),gs(2),gs(3),fpar%nTimes])
 
 !u
 ALLOCATE(u(gs(1)-1,gs(2),gs(3),fpar%nTimes)); u=dble(0)
 u=ncread4d(ncid,'u',[1,1,1,1],[gs(1)-1,gs(2),gs(3),fpar%nTimes])

 !v
 ALLOCATE(v(gs(1),gs(2)-1,gs(3),fpar%nTimes)); v=dble(0)
 v=ncread4d(ncid,'v',[1,1,1,1],[gs(1),gs(2)-1,gs(3),fpar%nTimes])
 
 !Time
 ALLOCATE(t(1)); t=dble(0)
 t=ncread1d(ncid,'ocean_time',[1],[1])

 !close stream
 status=nf_close(ncid)

!-------------------------------------------------------------------------
!Create output
 
 WRITE(*,*) 'creating output: ',TRIM(outi_file)

 !Times
 nt=fpar%nTimes
 
 !Create file
 fpar%nTimes=1
 CALL roms_create_his_file(TRIM(outi_file),fpar)

!---------------------------------------------------------------------------
!Write to output

 WRITE(*,*) 'writing to output'

 !Open stream
 status=nf_open(TRIM(outi_file),nf_write,ncid)

 !Time
 CALL ncwrite1d(ncid,'ocean_time',[1],[1],t)
 
 DO it=nt-1,1,-1
  zeta(:,:,it)=zeta(:,:,it)+zeta(:,:,it+1)
  temp(:,:,:,it)=temp(:,:,:,it)+temp(:,:,:,it+1)
  salt(:,:,:,it)=salt(:,:,:,it)+salt(:,:,:,it+1)
  u(:,:,:,it)=u(:,:,:,it)+u(:,:,:,it+1)
  v(:,:,:,it)=v(:,:,:,it)+v(:,:,:,it+1)
 END DO

 CALL ncwrite3d(ncid,'zeta',[1,1,1],[gs(1),gs(2),1],&
 &RESHAPE(zeta(:,:,1),[gs(1),gs(2),1]))
 CALL ncwrite4d(ncid,'temp',[1,1,1,1],[gs(1),gs(2),gs(3),1],&
 &RESHAPE(temp(:,:,:,1),[gs(1),gs(2),gs(3),1]))
 CALL ncwrite4d(ncid,'salt',[1,1,1,1],[gs(1),gs(2),gs(3),1],&
 &RESHAPE(salt(:,:,:,1),[gs(1),gs(2),gs(3),1]))
 CAlL ncwrite4d(ncid,'u',[1,1,1,1],[gs(1)-1,gs(2),gs(3),1],&
 &RESHAPE(u(:,:,:,1),[gs(1)-1,gs(2),gs(3),1]))
 CALL ncwrite4d(ncid,'v',[1,1,1,1],[gs(1),gs(2)-1,gs(3),1],&
 &RESHAPE(v(:,:,:,1),[gs(1),gs(2)-1,gs(3),1]))  

 !Close stream
 status=nf_close(ncid) 

 WRITE(*,*) 'DONE'
END PROGRAM ad_identity
