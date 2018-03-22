PROGRAM tl_identity
 USE mod_netcdf
 USE mod_roms

 !input
 CHARACTER(len=1024)::ini_file,out_file,grd_file
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
 READ(*,*) !background viscosity
 READ(*,*) 
 READ(*,*) !background diffusion
 READ(*,*)
 READ(*,*) !bottom friction coefficient
 READ(*,*)
 READ(*,*) !horizontal viscosity
 READ(*,*)
 READ(*,*) !horizontal diff
 READ(*,*) 
 READ(*,*) !Parameters of linear EOS
 READ(*,*)
 READ(*,*) !grid file
 READ(*,'(A)') grd_file
 READ(*,*) !Initial condition file (input)
 READ(*,'(A)') ini_file
 READ(*,*) !Initial forcing file (input)
 READ(*,*) 
 READ(*,*) !background file
 READ(*,*) 
 READ(*,*) !output file (output)
 READ(*,'(A)') out_file
 READ(*,*) !tiles
 READ(*,*) 

!---------------------------------------------------------------------------
!Read input

 WRITE(*,*) 'reading input file: ',TRIM(ini_file)

 !Read vertical structure
 CALL roms_get_type(TRIM(ini_file),fpar)
 gs=fpar%grid_size
 
 !open file 
 status=nf_open(TRIM(ini_file),nf_nowrite,ncid)

 !zeta
 ALLOCATE(zeta(gs(1),gs(2),1)); zeta=dble(0)
 zeta=ncread3d(ncid,'zeta',[1,1,1],[gs(1),gs(2),1])
  
 !temp
 ALLOCATE(temp(gs(1),gs(2),gs(3),1)); temp=dble(0)
 temp=ncread4d(ncid,'temp',[1,1,1,1],[gs(1),gs(2),gs(3),1])
 
 !salt
 ALLOCATE(salt(gs(1),gs(2),gs(3),1)); salt=dble(0)
 salt=ncread4d(ncid,'salt',[1,1,1,1],[gs(1),gs(2),gs(3),1])
 
 !u
 ALLOCATE(u(gs(1)-1,gs(2),gs(3),1)); u=dble(0)
 u=ncread4d(ncid,'u',[1,1,1,1],[gs(1)-1,gs(2),gs(3),1])

 !v
 ALLOCATE(v(gs(1),gs(2)-1,gs(3),1)); v=dble(0)
 v=ncread4d(ncid,'v',[1,1,1,1],[gs(1),gs(2)-1,gs(3),1])
 
 !close stream
 status=nf_close(ncid)

!-------------------------------------------------------------------------
!Create output

 WRITE(*,*) 'creating output file: ',TRIM(out_file)

 !Times
 nt=INT(ntimes/nhis)+1
 ALLOCATE(t(nt))
 DO i1=1,nt
  t(i1)=nhis*dt*dble(i1-1)
 END DO
 
 !Create file
 CALL roms_create_his_file(TRIM(out_file),fpar)

!---------------------------------------------------------------------------
!Write to output

 WRITE(*,*) 'writing to output'

 !Open stream
 status=nf_open(TRIM(out_file),nf_write,ncid)

 !Time
 CALL ncwrite1d(ncid,'ocean_time',[1],[nt],t)
 
 DO it=1,nt
  CALL ncwrite3d(ncid,'zeta',[1,1,it],[gs(1),gs(2),1],zeta)
  CALL ncwrite4d(ncid,'temp',[1,1,1,it],[gs(1),gs(2),gs(3),1],temp)
  CALL ncwrite4d(ncid,'salt',[1,1,1,it],[gs(1),gs(2),gs(3),1],salt)
  CALL ncwrite4d(ncid,'u',[1,1,1,it],[gs(1)-1,gs(2),gs(3),1],u)
  CALL ncwrite4d(ncid,'v',[1,1,1,it],[gs(1),gs(2)-1,gs(3),1],v)
 END DO
  
 !Close stream
 status=nf_close(ncid) 

 WRITE(*,*) 'DONE'

END PROGRAM tl_identity
