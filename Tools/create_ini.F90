PROGRAM create_ini
 USE mod_interp
 USE mod_roms
 USE mod_netcdf
 IMPLICIT NONE

 TYPE grd_type
  INTEGER::s(8)
  LOGICAL,allocatable,dimension(:,:)::maskr,masku,maskv
  REAL(8),allocatable,dimension(:)::lonr,lonu,lonv
  REAL(8),allocatable,dimension(:)::latr,latu,latv
  REAL(8),allocatable,dimension(:,:)::h
  REAL(8),allocatable,dimension(:,:,:)::zw,zr,zdr,zdu,zdv
 END TYPE grd_type
!-----------------------------------------------------------------------
!Declare


 real(8),allocatable::lon1(:),lat1(:),lon2(:),lat2(:)
 logical,allocatable::mask1(:,:),mask2(:,:)

 real(8),allocatable::val1d(:)
 !Input
 CHARACTER(len=1024)::nc_file(3),grid_file(3)
 INTEGER::itime(2)
 REAL(8)::time
 REAL::fac(2)
 LOGICAL::flag_exist(3),flag_runable
 CHARACTER::str_runable
 REAL(8)::T0,S0,Akv0,Akt0

 !Grid
 TYPE(grd_type)::grd(3)
 INTEGER::ilat(3),idgrid(3)
 REAL(8),allocatable::dist(:,:)
 REAL(8),allocatable::sigma(:)
 REAL(8),allocatable,dimension(:,:)::h1,h2
 REAL(8),allocatable,dimension(:,:,:)::zr1,zr2,zw1,zw2
 REAL(8),allocatable,dimension(:)::sr1,sr2,sr3,sw1,sw2,sw3

 !Netcdf, counters
 INTEGER::status,gid(3),nid(3),i0,i1,i2,i3,i4,nccount(4),ncstart(4)
 REAL::interp_fill
 REAL(8),allocatable,dimension(:,:,:)::val3d1,val3d2,val3d3,val3dt
 REAL(8),allocatable,dimension(:,:)::val2d1,val2d2,val2d3,val2dt
 REAL(8),allocatable,dimension(:,:)::h

 LOGICAL,allocatable::mask(:,:)

 TYPE(roms_type)::fdesign
 TYPE(sigma_param)::spar(3)
!------------------------------------------------------------------------
!Read input

 READ(*,*) !Grid files
 READ(*,'(A)') grid_file(1)
 READ(*,'(A)') grid_file(2)
 READ(*,'(A)') grid_file(3)
 READ(*,*) !Netcdf files
 READ(*,'(A)') nc_file(1)
 READ(*,'(A)') nc_file(2)
 READ(*,'(A)') nc_file(3)
 READ(*,*) !Time steps
 READ(*,*) itime(1),itime(2)
 READ(*,*) !Weight factors
 READ(*,*) fac(1),fac(2)
 READ(*,*) !Output time
 READ(*,*) time
 READ(*,*) !Default values T,S,Akt,Akv
 READ(*,*) T0,S0,Akt0,Akv0
 READ(*,*) !ROMS runable
 READ(*,'(A)') str_runable

 WRITE(*,*) 'grid 1:',TRIM(grid_file(1))
 WRITE(*,*) 'netcdf 1:',TRIM(nc_file(1))
 WRITE(*,*) 'time index 1:',itime(1)
 WRITE(*,*) 'multiplication factor 1:',fac(1)
 WRITE(*,*) 'grid 2:',TRIM(grid_file(2))
 WRITE(*,*) 'netcdf 2:',TRIM(nc_file(2))
 WRITE(*,*) 'time index 2:',itime(2)
 WRITE(*,*) 'multiplication factor 2:',fac(2)
 WRITE(*,*) 'grid out:',TRIM(grid_file(3))
 WRITE(*,*) 'netcdf out:',TRIM(nc_file(3))
 WRITE(*,*) 'time out:',time

 !Check if files exist. 
 INQUIRE(file=TRIM(nc_file(1)),exist=flag_exist(1))
 INQUIRE(file=TRIM(nc_file(2)),exist=flag_exist(2))
 INQUIRE(file=TRIM(nc_file(3)),exist=flag_exist(3))
 WRITE(*,*) 'flag_exist',flag_exist

 IF(str_runable.EQ.'T') THEN
  flag_runable=.TRUE.
 ELSE
  flag_runable=.FALSE.
 END IF

 IF( .NOT.flag_exist(1) ) THEN
  WRITE(*,*) 'Cannot find any input netcdf files.'
  STOP
 END IF

!----------------------------------------------------------------------------
! Read grid files

 WRITE(*,*) 'Reading grid files'

 DO i0=1,3
  
  !Open stream grid file
  IF(LEN_TRIM(grid_file(i0)).EQ.0) THEN
    status=nf_open(TRIM(grid_file(3)),nf_nowrite,gid(i0))
  ELSE
    status=nf_open(TRIM(grid_file(i0)),nf_nowrite,gid(i0))
  END IF
  
  !Allocate grid storage
  CALL ncsize(gid(i0),'z0_r',grd(i0)%s)
  ALLOCATE(grd(i0)%lonr(grd(i0)%s(1)))
  ALLOCATE(grd(i0)%lonu(grd(i0)%s(1)-1))
  ALLOCATE(grd(i0)%lonv(grd(i0)%s(1)))
  ALLOCATE(grd(i0)%latr(grd(i0)%s(2)))
  ALLOCATE(grd(i0)%latu(grd(i0)%s(2)))
  ALLOCATE(grd(i0)%latv(grd(i0)%s(2)-1))
  ALLOCATE(grd(i0)%zr(grd(i0)%s(1),grd(i0)%s(2),grd(i0)%s(3)))
  ALLOCATE(grd(i0)%zdr(grd(i0)%s(1),grd(i0)%s(2),grd(i0)%s(3)))
  ALLOCATE(grd(i0)%zdu(grd(i0)%s(1)-1,grd(i0)%s(2),grd(i0)%s(3)))
  ALLOCATE(grd(i0)%zdv(grd(i0)%s(1),grd(i0)%s(2)-1,grd(i0)%s(3)))
  ALLOCATE(grd(i0)%zw(grd(i0)%s(1),grd(i0)%s(2),grd(i0)%s(3)+1))
  ALLOCATE(grd(i0)%maskr(grd(i0)%s(1),grd(i0)%s(2)))
  ALLOCATE(grd(i0)%masku(grd(i0)%s(1)-1,grd(i0)%s(2)))
  ALLOCATE(grd(i0)%maskv(grd(i0)%s(1),grd(i0)%s(2)-1))

  IF(i0.EQ.3) THEN
  ALLOCATE(zr1(grd(3)%s(1),grd(3)%s(2),grd(3)%s(3)))
  ALLOCATE(zr2(grd(3)%s(1),grd(3)%s(2),grd(3)%s(3)))
  ALLOCATE(zw1(grd(3)%s(1),grd(3)%s(2),grd(3)%s(3)+1))
  ALLOCATE(zw2(grd(3)%s(1),grd(3)%s(2),grd(3)%s(3)+1))
  ALLOCATE(h1(grd(3)%s(1),grd(3)%s(2)))
  ALLOCATE(h2(grd(3)%s(1),grd(3)%s(2)))

  ALLOCATE(sr1(grd(1)%s(3))); sr1=dble(0)
  ALLOCATE(sr2(grd(2)%s(3))); sr2=dble(0)
  ALLOCATE(sr3(grd(3)%s(3))); sr3=dble(0)
  ALLOCATE(sw1(grd(1)%s(3)+1)); sw1=dble(0)
  ALLOCATE(sw2(grd(2)%s(3)+1)); sw2=dble(0)
  ALLOCATE(sw3(grd(3)%s(3)+1)); sw3=dble(0)
  END IF 

 !Read mask
  grd(i0)%maskr=ncread2d(gid(i0),'mask_rho',[1,1],&
  &[grd(i0)%s(1),grd(i0)%s(2)]).NE.0
  grd(i0)%masku=ncread2d(gid(i0),'mask_u',[1,1],&
  &[grd(i0)%s(1)-1,grd(i0)%s(2)]).NE.0
  grd(i0)%maskv=ncread2d(gid(i0),'mask_v',[1,1],&
  &[grd(i0)%s(1),grd(i0)%s(2)-1]).NE.0

  !Read lat
  grd(i0)%latr=RESHAPE(ncread2d(gid(i0),'lat_rho',[1,1],&
  &[1,grd(i0)%s(2)]),[grd(i0)%s(2)])
  grd(i0)%latu=RESHAPE(ncread2d(gid(i0),'lat_u',[1,1],&
  &[1,grd(i0)%s(2)]),[grd(i0)%s(2)])
  grd(i0)%latv=RESHAPE(ncread2d(gid(i0),'lat_v',[1,1],&
  &[1,grd(i0)%s(2)-1]),[grd(i0)%s(2)-1])
  
 END DO

 !Find latitude for which the difference with lon between the grids is minimal
 ALLOCATE(dist(size(grd(3)%latr),5)); dist=dble(0)
 DO i2=1,size(grd(3)%latr)
 dist(i2,1)=MINVAL(abs(grd(1)%latr-grd(3)%latr(i2)),1) 
 dist(i2,2)=MINLOC(abs(grd(1)%latr-grd(3)%latr(i2)),1)
 dist(i2,3)=MINVAL(abs(grd(2)%latr-grd(3)%latr(i2)),1)
 dist(i2,4)=MINLOC(abs(grd(2)%latr-grd(3)%latr(i2)),1)
 dist(i2,5)=dist(i2,1)**2+dist(i2,3)**2
 END DO
 ilat(3)=MINLOC(dist(:,5),1)
 ilat(1)=INT(dist(ilat(3),2)); ilat(2)=INT(dist(ilat(3),4))
 DEALLOCATE(dist)

 DO i0=1,3
  !Read lon
  grd(i0)%lonr=RESHAPE(ncread2d(gid(i0),'lon_rho',[1,ilat(i0)],&
  &[grd(i0)%s(1),1]),[grd(i0)%s(1)])
  grd(i0)%lonu=RESHAPE(ncread2d(gid(i0),'lon_u',[1,ilat(i0)],&
  &[grd(i0)%s(1)-1,1]),[grd(i0)%s(1)-1])
  grd(i0)%lonv=RESHAPE(ncread2d(gid(i0),'lon_v',[1,ilat(i0)],&
  &[grd(i0)%s(1),1]),[grd(i0)%s(1)])
 END DO


!-----------------------------------------------------------------------
!Create output file

 !Get grid parameters
 CALL roms_get_type(TRIM(nc_file(1)),fdesign)
 
 !Adapt to new grid
 WRITE(*,*) 'Creating output ',TRIM(nc_file(3))
 fdesign%ntimes=1
 fdesign%grid_size(1)=grd(3)%s(1)
 fdesign%grid_size(2)=grd(3)%s(2)
 fdesign%grid_size(3)=grd(3)%s(3)
 CALL roms_create_his_file(TRIM(nc_file(3)),fdesign)
 INQUIRE(file=TRIM(nc_file(3)),exist=flag_exist(3))

!-------------------------------------------------------------------------
!Interpolate and write

 !Open stream
 IF(flag_exist(1)) status=nf_open(TRIM(nc_file(1)),nf_nowrite,nid(1))
 IF(flag_exist(2)) status=nf_open(TRIM(nc_file(2)),nf_nowrite,nid(2))
 status=nf_open(TRIM(nc_file(3)),nf_write,nid(3))
 write(*,*) 'Status outfile:',status


 WRITE(*,*) 'Interpolating zeta' 
 !Zeta
 ALLOCATE(val2d1(grd(1)%s(1),grd(1)%s(2))); val2d1=dble(0)
 ALLOCATE(val2d2(grd(2)%s(1),grd(2)%s(2))); val2d2=dble(0)
 ALLOCATE(val2d3(grd(3)%s(1),grd(3)%s(2))); val2d3=interp_fill_value
 ALLOCATE(val2dt(grd(3)%s(1),grd(3)%s(2))); val2dt=dble(0)

 IF(flag_exist(1)) THEN
  ncstart=[1,1,itime(1),0]
  nccount=[grd(1)%s(1),grd(1)%s(2),1,0]
  val2d1=RESHAPE(ncread3d(nid(1),'zeta',ncstart(1:3),nccount(1:3)),&
  &[nccount(1:2)])
  
  ALLOCATE(lon1(size(grd(1)%lonr))); lon1=grd(1)%lonr
  ALLOCATE(lat1(size(grd(1)%latr))); lat1=grd(1)%latr
  ALLOCATE(mask1(size(grd(1)%maskr,1),size(grd(1)%maskr,2))); 
  mask1=grd(1)%maskr
  ALLOCATE(lon2(size(grd(3)%lonr))); lon2=grd(3)%lonr
  ALLOCATE(lat2(size(grd(3)%latr))); lat2=grd(3)%latr

  val2dt=mesh2mesh_interp2d(&
  &lon1,lat1,mask1,val2d1,lon2,lat2)
  WHERE(val2dt.NE.interp_fill_value)
   val2d3=fac(1)*val2dt
  END WHERE
 END IF

 IF(flag_exist(2)) THEN
  ncstart=[1,1,itime(2),0]
  nccount=[grd(2)%s(1),grd(2)%s(2),1,0]
  val2d2=RESHAPE(ncread3d(nid(2),'zeta',ncstart(1:3),nccount(1:3)),&
  &nccount(1:2))
  val2dt=mesh2mesh_interp(&
  &grd(2)%lonr,grd(2)%latr,grd(2)%maskr,val2d2,grd(3)%lonr,grd(3)%latr)
  WHERE(val2dt.NE.interp_fill_value.AND.val2d3.NE.interp_fill_value)
   val2d3=val2d3+fac(2)*val2dt
  ELSEWHERE(val2dt.NE.interp_fill_value.AND.val2d3.EQ.interp_fill_value)
   val2d3=fac(2)*val2dt
  END WHERE
 END IF


 !Extrapolate
 IF(flag_runable) THEN
  CALL mesh2mesh_extrap(grd(3)%lonr,grd(3)%latr,&
  &val2d3.NE.interp_fill_value,val2d3,grd(3)%maskr,dble(0)) 
  WHERE(val2d3.GT.dble(5)); val2d3=dble(5); END WHERE
  WHERE(val2d3.LT.dble(-5)); val2d3=dble(-5); END WHERE
 ELSE
  WHERE(val2d3.EQ.interp_fill_value.OR..NOT.grd(3)%maskr)
   val2d3=dble(0)
  END WHERE
 END IF

 !Write zeta
 write(*,*) 'min/max zeta:',minval(val2d3),maxval(val2d3)
 ncstart=[1,1,1,0]
 nccount=[grd(3)%s(1),grd(3)%s(2),1,0] 

 CALL ncwrite3d(nid(3),'zeta',[ncstart(1),ncstart(2),ncstart(3)],&
 &[nccount(1),nccount(2),nccount(3)],&
 &RESHAPE(val2d3,nccount(1:3)))

 !Sigma coordinates
 DO i1=1,size(sr1)
  sr1(i1)=dble(-1)+(dble(i1)-.5)/dble(size(sr1))
  sw1(i1)=dble(-1)+(dble(i1)-1.0)/dble(size(sr1))
 END DO
 DO i1=1,size(sr2)
  sr2(i1)=dble(-1)+(dble(i1)-.5)/dble(size(sr2))
  sw2(i1)=dble(-1)+(dble(i1)-1.0)/dble(size(sr2))
 END DO
 DO i1=1,size(sr3)
  sr3(i1)=dble(-1)+(dble(i1)-.5)/dble(size(sr3))
  sw3(i1)=dble(-1)+(dble(i1)-1.0)/dble(size(sr3))
 END DO

 !Calculate z for each grid
 DO i0=1,3
 
  IF(.NOT.flag_exist(i0)) CYCLE
  WRITE(*,*) 'Calculating z-grid ',i0

  !Get grid parameters
  status=nf_inq_varid(nid(i0),'Vtransform',i4)
  IF(status.EQ.nf_noerr) THEN
   CALL get_sigma_param(gid(i0),nid(i0),spar(i0),grd(i0)%h)
  ELSE
   spar(i0)=spar(1)
   ALLOCATE(grd(i0)%h(grd(i0)%s(1),grd(i0)%s(2)))
   grd(i0)%h=ncread2d(gid(i0),'h',[1,1],[grd(i0)%s(1),grd(i0)%s(2)])
  END IF   
  status=nf_close(gid(i0))
 
  !h at points grid 3
  h1=mesh2mesh_interp(grd(1)%lonr,grd(1)%latr,&
  &grd(1)%maskr.OR..NOT.grd(1)%maskr,&
  &grd(1)%h,grd(3)%lonr,grd(3)%latr)
  h2=mesh2mesh_interp(grd(2)%lonr,grd(2)%latr,&
  &grd(2)%maskr.OR..NOT.grd(2)%maskr,&
  &grd(1)%h,grd(3)%lonr,grd(3)%latr)
  WHERE(h1.EQ.interp_fill_value); h1=MINVAL(h1); END WHERE
  WHERE(h2.EQ.interp_fill_value); h2=MINVAL(h2); END WHERE

  !z at rho points
  IF(i0.EQ.1) grd(i0)%zr=sigma2z_2d(spar(i0),grd(i0)%h,val2d1,sr1)
  IF(i0.EQ.2) grd(i0)%zr=sigma2z_2d(spar(i0),grd(i0)%h,val2d2,sr2)
  IF(i0.EQ.3) grd(i0)%zr=sigma2z_2d(spar(i0),grd(i0)%h,val2d3,sr3)
  zr1=sigma2z_2d(spar(1),h1,val2d3,sr1)
  zr2=sigma2z_2d(spar(2),h2,val2d3,sr2)

  !z at w points
  IF(.NOT.flag_runable) THEN
    val2d1=dble(0); val2d2=dble(0); val2d3=dble(0)
  END IF
  IF(i0.EQ.1) grd(i0)%zw=sigma2z_2d(spar(i0),grd(i0)%h,val2d1,sw1)
  IF(i0.EQ.2) grd(i0)%zw=sigma2z_2d(spar(i0),grd(i0)%h,val2d2,sw2)
  IF(i0.EQ.3) grd(i0)%zw=sigma2z_2d(spar(i0),grd(i0)%h,val2d3,sw3) 
  zw1=sigma2z_2d(spar(1),h1,val2d3,sw1)
  zw2=sigma2z_2d(spar(2),h2,val2d3,sw2)

  !Vertical dimension grid cells
  grd(i0)%zdr=grd(i0)%zw(:,:,2:grd(i0)%s(3)+1)-&
  &grd(i0)%zw(:,:,1:grd(i0)%s(3))
  grd(i0)%zdu=.5*grd(i0)%zdr(1:grd(i0)%s(1)-1,:,:)+&
  &.5*grd(i0)%zdr(2:grd(i0)%s(1),:,:)
  grd(i0)%zdv=.5*grd(i0)%zdr(:,1:grd(i0)%s(2)-1,:)+&
  &.5*grd(i0)%zdr(:,2:grd(i0)%s(2),:)

 END DO
 
 DEALLOCATE(val2d1,val2d2,val2d3,val2dt)

!------------------------------------------------------------------------
!Interpolate u-velocities

 WRITE(*,*) 'Interpolating u'

 ALLOCATE(val3d1(grd(1)%s(1)-1,grd(1)%s(2),grd(1)%s(3))); val3d1=dble(0)
 ALLOCATE(val3d2(grd(2)%s(1)-1,grd(2)%s(2),grd(2)%s(3))); val3d2=dble(0)
 ALLOCATE(val3d3(grd(3)%s(1)-1,grd(3)%s(2),grd(3)%s(3)));
 val3d3=interp_fill_value
 ALLOCATE(val3dt(grd(3)%s(1)-1,grd(3)%s(2),grd(3)%s(3))); val3dt=dble(0)
 ALLOCATE(val2d1(grd(1)%s(1)-1,grd(1)%s(2))); val2d1=dble(0)
 ALLOCATE(val2d2(grd(2)%s(1)-1,grd(2)%s(2))); val2d2=dble(0)
 ALLOCATE(val2d3(grd(3)%s(1)-1,grd(3)%s(2))); val2d3=interp_fill_value
 ALLOCATE(val2dt(grd(3)%s(1)-1,grd(3)%s(2))); val2dt=dble(0)

 !Read u-velocities
 IF(flag_exist(1)) THEN
  ncstart=[1,1,1,itime(1)]
  nccount=[grd(1)%s(1)-1,grd(1)%s(2),grd(1)%s(3),1] 
  val3d1=RESHAPE(ncread4d(nid(1),'u',ncstart,nccount),&
  &nccount(1:3))
  val3dt=mesh2mesh_interp(&
  &grd(1)%lonu,grd(1)%latu,grd(1)%masku,val3d1,grd(3)%lonu,grd(3)%latu)

   !Write to output
   WHERE(val3dt.NE.interp_fill_value)
    val3d3=fac(1)*val3dt
   END WHERE
 END IF

 IF(flag_exist(2)) THEN
  ncstart=[1,1,1,itime(2)]
  nccount=[grd(2)%s(1)-1,grd(2)%s(2),grd(2)%s(3),1] 
  val3d2=RESHAPE(ncread4d(nid(2),'u',ncstart,nccount),&
  &nccount(1:3))
  val3dt=mesh2mesh_interp(&
  &grd(2)%lonu,grd(2)%latu,grd(2)%masku,val3d2,grd(3)%lonu,grd(3)%latu)

 !Write to output
  WHERE(val3dt.NE.interp_fill_value.AND.val3d3.NE.interp_fill_value)
   val3d3=val3d3+fac(2)*val3dt
  ELSEWHERE(val3dt.NE.interp_fill_value.AND.val3d3.EQ.interp_fill_value)
   val3d3=fac(2)*val3dt
  END WHERE
 END IF

 !Extrapolate u
 IF(flag_runable) THEN
  CALL mesh2mesh_extrap(grd(3)%lonu,grd(3)%latr,&
  &.NOT.ANY(val3d3.EQ.interp_fill_value,3),val3d3,grd(3)%masku,dble(0))
 ELSE
  WHERE(val3d3.EQ.interp_fill_value.OR.&
  &SPREAD(.NOT.grd(3)%masku,3,size(val3d3,3)))
   val3d3=dble(0)
  END WHERE
 END IF

 !Calculate depth averaged currents
 WRITE(*,*) 'Interpolating ubar'
 val2d1=SUM(val3d1*grd(1)%zdu,3)/SUM(grd(1)%zdu,3)
 WHERE(SUM(grd(1)%zdu,3).EQ.dble(0)); val2d1=dble(0); END WHERE
 val2d2=SUM(val3d2*grd(2)%zdu,3)/SUM(grd(2)%zdu,3)
 WHERE(SUM(grd(2)%zdu,3).EQ.dble(0)); val2d2=dble(0); END WHERE
 val2dt=mesh2mesh_interp(grd(1)%lonu,grd(1)%latu,grd(1)%masku,&
 &val2d1,grd(3)%lonu,grd(3)%latu)
 WHERE(val2dt.NE.interp_fill_value)
  val2d3=fac(1)*val2dt
 END WHERE
 val2dt=mesh2mesh_interp(grd(2)%lonu,grd(2)%latu,grd(2)%masku,&
 &val2d2,grd(3)%lonu,grd(3)%latu)
 WHERE(val2dt.NE.interp_fill_value.AND.val2d3.NE.interp_fill_value)
  val2d3=val2d3+fac(2)*val2dt
 ELSEWHERE(val2dt.NE.interp_fill_value.AND.val2d3.EQ.interp_fill_value)
  val2d3=fac(2)*val2dt
 END WHERE

 !Extrapolate ubar
 IF(flag_runable) THEN
  CALL mesh2mesh_extrap(grd(3)%lonu,grd(3)%latu,&
  &val2d3.NE.interp_fill_value,val2d3,grd(3)%masku,dble(0))
  WHERE(val2d3.GT.dble(2.5)); val2d3=dble(2.5); END WHERE
  WHERE(val2d3.LT.dble(-2.5)); val2d3=dble(-2.5); END WHERE
 ELSE
  WHERE(val2d3.EQ.interp_fill_value.OR..NOT.grd(3)%masku)
   val2d3=dble(0)
  END WHERE
 END IF

 write(*,*) 'min/max ubar:',minval(val2d3),maxval(val2d3)
 ncstart=[1,1,1,0]
 nccount=[grd(3)%s(1)-1,grd(3)%s(2),1,0]
 CALL ncwrite3d(nid(3),'ubar',ncstart(1:3),nccount(1:3),&
 &RESHAPE(val2d3,nccount(1:3)))

 !Rescale u-velocities such that depth-averaged currents match
 WRITE(*,*) 'Rescaling u'
 val3dt=SPREAD(SUM(val3d3*grd(3)%zdu,3),&
 &3,SIZE(val3d3,3))
 WHERE(val3dt.NE.dble(0))
  !val3d3=val3d3/val3dt*&
  !&SPREAD(SUM(grd(3)%zdu,3)*val2d3,3,size(val3d3,3))
 END WHERE

 IF(flag_runable) THEN
  WHERE(val3d3.GT.dble(2.5)); val3d3=dble(2.5); END WHERE
  WHERE(val3d3.LT.dble(-2.5)); val3d3=dble(-2.5); END WHERE
 END IF

 write(*,*) 'min/max u:',minval(val3d3),maxval(val3d3)
 ncstart=[1,1,1,1]
 nccount=[grd(3)%s(1)-1,grd(3)%s(2),grd(3)%s(3),1]
 CALL ncwrite4d(nid(3),'u',ncstart,nccount,&
 &RESHAPE(val3d3,nccount))

 !deallocate
 DEALLOCATE(val3d1,val3d2,val3d3,val3dt,val2d1,val2d2,val2d3,val2dt)


!--------------------------------------------------------------------
!Interpolate v-velocities

 WRITE(*,*) 'Interpolating v'

 ALLOCATE(val3d1(grd(1)%s(1),grd(1)%s(2)-1,grd(1)%s(3))); val3d1=dble(0)
 ALLOCATE(val3d2(grd(2)%s(1),grd(2)%s(2)-1,grd(2)%s(3))); val3d2=dble(0)
 ALLOCATE(val3d3(grd(3)%s(1),grd(3)%s(2)-1,grd(3)%s(3))); 
 val3d3=interp_fill_value
 ALLOCATE(val3dt(grd(3)%s(1),grd(3)%s(2)-1,grd(3)%s(3))); val3dt=dble(0)
 ALLOCATE(val2d1(grd(1)%s(1),grd(1)%s(2)-1)); val2d1=dble(0)
 ALLOCATE(val2d2(grd(2)%s(1),grd(2)%s(2)-1)); val2d2=dble(0)
 ALLOCATE(val2d3(grd(3)%s(1),grd(3)%s(2)-1)); val2d3=interp_fill_value
 ALLOCATE(val2dt(grd(3)%s(1),grd(3)%s(2)-1)); val2dt=dble(0) 

 !Read v-velocities
 IF(flag_exist(1)) THEN
  !Read
  ncstart=[1,1,1,itime(1)]
  nccount=[grd(1)%s(1),grd(1)%s(2)-1,grd(1)%s(3),1] 
  val3d1=RESHAPE(ncread4d(nid(1),'v',ncstart,nccount),&
  &nccount(1:3))
  val3dt=mesh2mesh_interp(&
  &grd(1)%lonv,grd(1)%latv,grd(1)%maskv,val3d1,grd(3)%lonv,grd(3)%latv)

  !Write to output
  WHERE(val3dt.NE.interp_fill_value)
   val3d3=fac(1)*val3dt
  END WHERE
 END IF

 IF(flag_exist(2)) THEN
  !Read
  ncstart=[1,1,1,itime(2)]
  nccount=[grd(2)%s(1),grd(2)%s(2)-1,grd(2)%s(3),1] 
  val3d2=RESHAPE(ncread4d(nid(2),'v',ncstart,nccount),&
  &nccount(1:3))
  val3dt=mesh2mesh_interp(&
  &grd(2)%lonv,grd(2)%latv,grd(2)%maskv,val3d2,grd(3)%lonv,grd(3)%latv)
 
  !Write to output
  WHERE(val3dt.NE.interp_fill_value.AND.val3d3.NE.interp_fill_value)
   val3d3=val3d3+fac(2)*val3dt
  ELSEWHERE(val3dt.NE.interp_fill_value.AND.val3d3.EQ.interp_fill_value)
   val3d3=fac(2)*val3dt
  END WHERE
 END IF

 !Extrapolate v
 IF(flag_runable) THEN
  CALL mesh2mesh_extrap(grd(3)%lonv,grd(3)%latv,&
  &.NOT.ANY(val3d3.EQ.interp_fill_value,3),val3d3,grd(3)%maskv,dble(0))
 ELSE
  WHERE(val3d3.EQ.interp_fill_value.OR.SPREAD(.NOT.grd(3)%maskv,&
  &3,size(val3d3,3)))
   val3d3=dble(0)
  END WHERE
 END IF 

 !Calculate depth averaged currents
 WRITE(*,*) 'Interpolating vbar'
 val2d1=SUM(val3d1*grd(1)%zdv,3)/SUM(grd(1)%zdv,3)
 WHERE(SUM(grd(1)%zdv,3).EQ.dble(0)); val2d1=dble(0); END WHERE
 val2d2=SUM(val3d2*grd(2)%zdv,3)/SUM(grd(2)%zdv,3)
 WHERE(SUM(grd(2)%zdv,3).EQ.dble(0)); val2d2=dble(0); END WHERE
 val2dt=mesh2mesh_interp(grd(1)%lonv,grd(1)%latv,grd(1)%maskv,&
 &val2d1,grd(3)%lonv,grd(3)%latv)
 WHERE(val2dt.NE.interp_fill_value)
  val2d3=fac(1)*val2dt
 END WHERE 
 val2dt=mesh2mesh_interp(grd(2)%lonv,grd(2)%latv,grd(2)%maskv,&
 &val2d2,grd(3)%lonv,grd(3)%latv)
 WHERE(val2dt.NE.interp_fill_value.AND.val2d3.NE.interp_fill_value)
  val2d3=val2d3+fac(2)*val2dt
 ELSEWHERE(val2dt.NE.interp_fill_value.AND.val2d3.EQ.interp_fill_value)
  val2d3=fac(2)*val2dt
 END WHERE 

 !Extrapolate vbar
 IF(flag_runable) THEN
  CALL mesh2mesh_extrap(grd(3)%lonv,grd(3)%latv,&
  &val2d3.NE.interp_fill_value,val2d3,grd(3)%maskv,dble(0))
  WHERE(val2d3.GT.dble(2.5)); val2d3=dble(2.5); END WHERE
  WHERE(val2d3.LT.dble(-2.5)); val2d3=dble(-2.5); END WHERE
 ELSE
  WHERE(val2d3.EQ.interp_fill_value.OR..NOT.grd(3)%maskv)
   val2d3=dble(0)
  END WHERE
 END IF

 WRITE(*,*) 'min/max vbar:',minval(val2d3),maxval(val2d3)
 ncstart=[1,1,1,0]
 nccount=[grd(3)%s(1),grd(3)%s(2)-1,1,0]
 CALL ncwrite3d(nid(3),'vbar',ncstart(1:3),nccount(1:3),&
 &RESHAPE(val2d3,nccount(1:3)))

 !Rescale v-velocities such that depth-averaged currents match
 WRITE(*,*) 'Rescaling v'
 val3dt=SPREAD(SUM(val3d3*grd(3)%zdv,3),&
 &3,SIZE(val3d3,3))
 WHERE(val3dt.NE.dble(0))
  !val3d3=val3d3/val3dt*&
  !&SPREAD(val2d3*SUM(grd(3)%zdv,3),3,size(val3d3,3))
 END WHERE

 IF(flag_runable) THEN
  WHERE(val3d3.GT.dble(2.5)); val3d3=dble(2.5); END WHERE
  WHERE(val3d3.LT.dble(-2.5)); val3d3=dble(-2.5); END WHERE
 END IF

 WRITE(*,*) 'min/max v:',minval(val3d3),maxval(val3d3)
 ncstart=[1,1,1,1]
 nccount=[grd(3)%s(1),grd(3)%s(2)-1,grd(3)%s(3),1]
 CALL ncwrite4d(nid(3),'v',ncstart,nccount,&
 &RESHAPE(val3d3,nccount))

 DEALLOCATE(val2d1,val2d2,val2d3,val2dt,val3d1,val3d2,val3d3,val3dt)

!-------------------------------------------------------------------------
!Interpolating temp 
 WRITE(*,*) 'Interpolating temp'

 ALLOCATE(val3d1(grd(1)%s(1),grd(1)%s(2),grd(1)%s(3))); val3d1=dble(0)
 ALLOCATE(val3d2(grd(2)%s(1),grd(2)%s(2),grd(2)%s(3))); val3d2=dble(0)
 ALLOCATE(val3d3(grd(3)%s(1),grd(3)%s(2),grd(3)%s(3)));
 val3d3=interp_fill_value
 ALLOCATE(val3dt(grd(3)%s(1),grd(3)%s(2),grd(3)%s(3))); val3dt=dble(0)
 ALLOCATE(val1d(size(val3d3,3)))

 !Read temp
 IF(flag_exist(1)) THEN
  ncstart=[1,1,1,itime(1)]
  nccount=[grd(1)%s(1),grd(1)%s(2),grd(1)%s(3),1] 
  val3d1=RESHAPE(ncread4d(nid(1),'temp',ncstart,nccount),&
  &nccount(1:3))
  val3dt=mesh2mesh_interp(&
  &grd(1)%lonr,grd(1)%latr,grd(1)%maskr,val3d1,grd(3)%lonr,grd(3)%latr)

  !Write to output
  WHERE(val3dt.NE.interp_fill_value)
   val3d3=fac(1)*val3dt
  END WHERE

 END IF

 IF(flag_exist(2)) THEN
  ncstart=[1,1,1,itime(2)]
  nccount=[grd(2)%s(1),grd(2)%s(2),grd(2)%s(3),1] 
  val3d2=RESHAPE(ncread4d(nid(2),'temp',ncstart,nccount),&
  &nccount(1:3))
  val3dt=mesh2mesh_interp(&
  &grd(2)%lonr,grd(2)%latr,grd(2)%maskr,val3d2,grd(3)%lonr,grd(3)%latr)

  !Write to output
  WHERE(val3d3.NE.interp_fill_value.AND.val3dt.NE.interp_fill_value)
   val3d3=val3d3+fac(2)*val3dt
  ELSEWHERE(val3d3.EQ.interp_fill_value.AND.val3dt.NE.interp_fill_value)
   val3d3=fac(2)*val3dt
  END WHERE
 END IF

 !Extrapolation
 IF(flag_runable) THEN
  CALL mesh2mesh_extrap(grd(3)%lonr,grd(3)%latr,&
  &.NOT.ANY(val3d3.EQ.interp_fill_value,3),val3d3,grd(3)%maskr,T0)
  WHERE(val3d3.GT.dble(40)); val3d3=dble(40); END WHERE
  WHERE(val3d3.LT.dble(0)); val3d3=dble(0); END WHERE
 ELSE
  WHERE(SPREAD(.NOT.grd(3)%maskr,3,size(val3d3,3)).OR.&
  &val3d3.EQ.interp_fill_value)
   val3d3=dble(0)
  END WHERE
 END IF

 !Write temperature
 write(*,*) 'min/max temp:', minval(val3d3),maxval(val3d3)
 ncstart=[1,1,1,1]
 nccount=[grd(3)%s(1),grd(3)%s(2),grd(3)%s(3),1]
 CALL ncwrite4d(nid(3),'temp',ncstart,nccount,&
 &RESHAPE(val3d3,nccount))

!----------------------------------------------------------------------------
!Interpolating salinity

 WRITE(*,*) 'Interpolating salt'
 val3d1=dble(0); val3d2=dble(0); val3d3=interp_fill_value; 
 val3dt=dble(0)

 !Read salinity
 IF(flag_exist(1)) THEN
  ncstart=[1,1,1,itime(1)]
  nccount=[grd(1)%s(1),grd(1)%s(2),grd(1)%s(3),1] 
  val3d1=RESHAPE(ncread4d(nid(1),'salt',ncstart,nccount),&
  &nccount(1:3))
  val3dt=mesh2mesh_interp(&
  &grd(1)%lonr,grd(1)%latr,grd(1)%maskr,val3d1,grd(3)%lonr,grd(3)%latr)

  !Write to output
  WHERE(val3dt.NE.interp_fill_value)
   val3d3=fac(1)*val3dt
  END WHERE
 END IF

 IF(flag_exist(2)) THEN
  ncstart=[1,1,1,itime(2)]
  nccount=[grd(2)%s(1),grd(2)%s(2),grd(2)%s(3),1] 
  val3d2=RESHAPE(ncread4d(nid(2),'salt',ncstart,nccount),&
  &nccount(1:3))
  val3dt=mesh2mesh_interp(&
  &grd(2)%lonr,grd(2)%latr,grd(2)%maskr,val3d2,grd(3)%lonr,grd(3)%latr)

  !Write to output
  WHERE(val3dt.NE.interp_fill_value.AND.val3d3.NE.interp_fill_value)
   val3d3=val3d3+fac(2)*val3dt
  ELSEWHERE(val3dt.NE.interp_fill_value.AND.val3d3.NE.interp_fill_value)
   val3d3=fac(2)*val3dt
  END WHERE

 END IF

 !Extrapolate
 IF(flag_runable) THEN
  CALL mesh2mesh_extrap(grd(3)%lonr,grd(3)%latr,&
  &.NOT.ANY(val3d3.EQ.interp_fill_value,3),val3d3,grd(3)%maskr,S0)
  WHERE(val3d3.GT.dble(40)); val3d3=dble(40); END WHERE
  WHERE(val3d3.LT.dble(0)); val3d3=dble(0); END WHERE
 ELSE
  WHERE(SPREAD(.NOT.grd(3)%maskr,3,size(val3d3,3)).OR.&
  &val3d3.EQ.interp_fill_value)
   val3d3=dble(0)
  END WHERE
 END IF 

 !Write salinity
 WRITE(*,*) 'min/max salt:',minval(val3d3),maxval(val3d3)
 ncstart=[1,1,1,1]
 nccount=[grd(3)%s(1),grd(3)%s(2),grd(3)%s(3),1]
 CALL ncwrite4d(nid(3),'salt',ncstart,nccount,&
 &RESHAPE(val3d3,nccount))

 DEALLOCATE(val3d1,val3d2,val3d3,val3dt)
 
!-----------------------------------------------------------------------
!Interpolating AKt

 WRITE(*,*) 'Interpolating AKt'
 
 ALLOCATE(val3d1(grd(1)%s(1),grd(1)%s(2),grd(1)%s(3)+1)); val3d1=dble(0)
 ALLOCATE(val3d2(grd(2)%s(1),grd(2)%s(2),grd(2)%s(3)+1)); val3d2=dble(0)
 ALLOCATE(val3d3(grd(3)%s(1),grd(3)%s(2),grd(3)%s(3)+1))
 val3d3=interp_fill_value
 ALLOCATE(val3dt(grd(3)%s(1),grd(3)%s(2),grd(3)%s(3)+1)); val3dt=dble(0)

 !Read AKt
 IF(flag_exist(1).AND.nf_inq_varid(nid(1),'AKt',i0).EQ.nf_noerr) THEN
  ncstart=[1,1,1,itime(1)]
  nccount=[grd(1)%s(1),grd(1)%s(2),grd(1)%s(3)+1,1] 
  val3d1=RESHAPE(ncread4d(nid(1),'AKt',ncstart,nccount),&
  &nccount(1:3))
  val3dt=mesh2mesh_interp(&
  &grd(1)%lonr,grd(1)%latr,grd(1)%maskr,val3d1,grd(3)%lonr,grd(3)%latr)

  !Write to output
  WHERE(val3dt.NE.interp_fill_value)
   val3d3=fac(1)*val3d1
  END WHERE

 END IF
 IF(flag_exist(2).AND.nf_inq_varid(nid(2),'AKt',i0).EQ.nf_noerr) THEN
  ncstart=[1,1,1,itime(2)]
  nccount=[grd(2)%s(1),grd(2)%s(2),grd(2)%s(3)+1,1] 
  val3d2=RESHAPE(ncread4d(nid(2),'AKt',ncstart,nccount),&
  &nccount(1:3))
  val3dt=mesh2mesh_interp(&
  &grd(2)%lonr,grd(2)%latr,grd(2)%maskr,val3d2,grd(3)%lonr,grd(3)%latr)

  WHERE(val3dt.NE.interp_fill_value.AND.val3d3.NE.interp_fill_value)
   val3d3=val3d3+fac(2)*val3dt
  ELSEWHERE(val3dt.NE.interp_fill_value.AND.val3d3.NE.interp_fill_value)
   val3d3=fac(2)*val3dt
  END WHERE
 END IF

 !Extrapolate
 IF(flag_runable) THEN
  WHERE(SPREAD(.NOT.grd(3)%maskr,3,size(val3d3,3)).OR.&
  &val3d3.EQ.interp_fill_value)
   val3d3=Akt0
  END WHERE
  val3d3=MAX(val3d3,Akt0)
 ELSE
  WHERE(SPREAD(.NOT.grd(3)%maskr,3,size(val3d3,3)).OR.&
  &val3d3.EQ.interp_fill_value)
   val3d3=dble(0)
  END WHERE
 END IF

 !Write Akt
 WRITE(*,*) 'min/max AKt:',minval(val3d3),maxval(val3d3)
 ncstart=[1,1,1,1]
 nccount=[grd(3)%s(1),grd(3)%s(2),grd(3)%s(3)+1,1]
 CALL ncwrite4d(nid(3),'AKt',ncstart,nccount,&
 &RESHAPE(val3d3,nccount))

!-----------------------------------------------------------------------
!Interpolating AKv

 WRITE(*,*) 'Interpolatin AKv'
 val3d1=dble(0); val3d2=dble(0); val3dt=dble(0)
 val3d3=interp_fill_value 

 !Read AKv
 IF(flag_exist(1).AND.nf_inq_varid(nid(1),'AKv',i0).EQ.nf_noerr) THEN
  ncstart=[1,1,1,itime(1)]
  nccount=[grd(1)%s(1),grd(1)%s(2),grd(1)%s(3)+1,1] 
  val3d1=RESHAPE(ncread4d(nid(1),'AKv',ncstart,nccount),&
  &nccount(1:3))
  val3dt=mesh2mesh_interp(&
  &grd(1)%lonr,grd(1)%latr,grd(1)%maskr,val3d1,grd(3)%lonr,grd(3)%latr)

  !Write to output
  WHERE(val3dt.NE.interp_fill_value)
   val3d3=fac(1)*val3dt
  END WHERE
 END IF

 IF(flag_exist(2).AND.nf_inq_varid(nid(2),'AKv',i0).EQ.nf_noerr) THEN
  ncstart=[1,1,1,itime(2)]
  nccount=[grd(2)%s(1),grd(2)%s(2),grd(2)%s(3)+1,1] 
  val3d2=RESHAPE(ncread4d(nid(2),'AKv',ncstart,nccount),&
  &nccount(1:3))
  val3dt=mesh2mesh_interp(&
  &grd(2)%lonr,grd(2)%latr,grd(2)%maskr,val3d2,grd(3)%lonr,grd(3)%latr)

  !Write to output
  WHERE(val3dt.NE.interp_fill_value.AND.val3d3.NE.interp_fill_value)
   val3d3=val3d3+fac(2)*val3dt
  ELSEWHERE(val3dt.NE.interp_fill_value.AND.val3d3.NE.interp_fill_value)
   val3d3=fac(2)*val3dt
  END WHERE

 END IF

 IF(flag_runable) THEN
  WHERE(SPREAD(.NOT.grd(3)%maskr,3,size(val3d3,3)).OR.&
  &val3d3.EQ.interp_fill_value)
   val3d3=Akv0
  END WHERE
  val3d3=MAX(val3d3,Akv0)
 ELSE
  WHERE(SPREAD(.NOT.grd(3)%maskr,3,size(val3d3,3)).OR.&
  &val3d3.EQ.interp_fill_value)
  val3d3=dble(0)
  END WHERE
 END IF

 !Write Akv
 WRITE(*,*) 'min/max AKv:',minval(val3d3),maxval(val3d3)
 ncstart=[1,1,1,1]
 nccount=[grd(3)%s(1),grd(3)%s(2),grd(3)%s(3)+1,1]
 CALL ncwrite4d(nid(3),'AKv',ncstart,nccount,&
 &RESHAPE(val3d3,nccount))

 DEALLOCATE(val3d1,val3d2,val3d3,val3dt)

!--------------------------------------------------------------------------
!Write time

 CALL ncwrite1d(nid(3),'ocean_time',[1],[1],[time])

!---------------------------------------------------------------------------
!Close open streams

 WRITE(*,*) 'Close streams'

 DO i0=1,3
  status=nf_close(nid(i0))
 END DO
 
 WRITE(*,*) 'create_ini DONE'


END PROGRAM create_ini
