PROGRAM tl_interp
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

 !Input
 CHARACTER(len=1024)::nc_file(3),grid_file(2)
 INTEGER::itime(2)
 REAL(8)::time
 REAL::fac(2)
 LOGICAL::flag_exist(2)
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
 INTEGER::status,gid(2),nid(2),i0,i1,i2,i3,i4,nccount(4),ncstart(4)
 REAL(8),allocatable,dimension(:,:,:)::val3d1,val3d2,val3d3,val3dt
 REAL(8),allocatable,dimension(:,:)::val2d1,val2d2,val2d3,val2dt
 REAL(8),allocatable,dimension(:,:)::h

 LOGICAL,allocatable::mask(:,:)

 TYPE(roms_type)::fdesign
 TYPE(sigma_param)::spar(2)
!------------------------------------------------------------------------
!Read input

 READ(*,*) !Grid files
 READ(*,'(A)') grid_file(1)
 READ(*,'(A)') grid_file(2)
 READ(*,*) !Netcdf files
 READ(*,'(A)') nc_file(1)
 READ(*,'(A)') nc_file(2)
 READ(*,*) !Time steps
 READ(*,*) itime(1),itime(2)
 READ(*,*) !Output time
 READ(*,*) time
 READ(*,*) !Background file
 READ(*,'(A)') nc_file(3)

 WRITE(*,*) 'grid 1:',TRIM(grid_file(1))
 WRITE(*,*) 'netcdf 1:',TRIM(nc_file(1))
 WRITE(*,*) 'time index 1:',itime(1)
 WRITE(*,*) 'grid out:',TRIM(grid_file(2))
 WRITE(*,*) 'netcdf out:',TRIM(nc_file(2))
 WRITE(*,*) 'time out:',time
 WRITE(*,*) 'background file:',TRIM(nc_file(3))

!----------------------------------------------------------------------------
! Read grid files

 WRITE(*,*) 'Reading grid files'

 DO i0=1,2
  
  status=nf_open(TRIM(grid_file(i0)),nf_nowrite,gid(i0))

  !Allocate grid storage
  CALL ncsize(gid(i0),'z0_r',grd(i0)%s)
  ALLOCATE(grd(i0)%lonr(grd(i0)%s(1)))
  ALLOCATE(grd(i0)%lonu(grd(i0)%s(1)-1))
  ALLOCATE(grd(i0)%lonv(grd(i0)%s(1)))
  ALLOCATE(grd(i0)%latr(grd(i0)%s(2)))
  ALLOCATE(grd(i0)%latu(grd(i0)%s(2)))
  ALLOCATE(grd(i0)%latv(grd(i0)%s(2)-1))
  ALLOCATE(grd(i0)%maskr(grd(i0)%s(1),grd(i0)%s(2)))
  ALLOCATE(grd(i0)%masku(grd(i0)%s(1)-1,grd(i0)%s(2)))
  ALLOCATE(grd(i0)%maskv(grd(i0)%s(1),grd(i0)%s(2)-1))
  ALLOCATE(grd(i0)%zr(grd(i0)%s(1),grd(i0)%s(2),grd(i0)%s(3)))
  ALLOCATE(grd(i0)%zw(grd(i0)%s(1),grd(i0)%s(2),grd(i0)%s(3)+1))

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
  
  !Read z
  ncstart=[1,1,1,0]; nccount=[grd(i0)%s(1),grd(i0)%s(2),grd(i0)%s(3),0]
  grd(i0)%zr=ncread3d(gid(i0),'z0_r',ncstart(1:3),nccount(1:3))
  ncstart=[1,1,1,0]; nccount=[grd(i0)%s(1),grd(i0)%s(2),grd(i0)%s(3)+1,0]
  grd(i0)%zw=ncread3d(gid(i0),'z0_w',ncstart(1:3),nccount(1:3))

 END DO

 !Find latitude for which the difference with lon between the grids is minimal
 ALLOCATE(dist(size(grd(2)%latr),3)); dist=dble(0)
 DO i2=1,size(grd(2)%latr)
 dist(i2,1)=MINVAL(abs(grd(1)%latr-grd(2)%latr(i2)),1) 
 dist(i2,2)=MINLOC(abs(grd(1)%latr-grd(2)%latr(i2)),1)
 dist(i2,3)=dist(i2,1)**2
 END DO
 ilat(2)=MINLOC(dist(:,3),1)
 ilat(1)=INT(dist(ilat(2),2)); 
 DEALLOCATE(dist)

 DO i0=1,2
  !Read lon
  grd(i0)%lonr=RESHAPE(ncread2d(gid(i0),'lon_rho',[1,ilat(i0)],&
  &[grd(i0)%s(1),1]),[grd(i0)%s(1)])
  grd(i0)%lonu=RESHAPE(ncread2d(gid(i0),'lon_u',[1,ilat(i0)],&
  &[grd(i0)%s(1)-1,1]),[grd(i0)%s(1)-1])
  grd(i0)%lonv=RESHAPE(ncread2d(gid(i0),'lon_v',[1,ilat(i0)],&
  &[grd(i0)%s(1),1]),[grd(i0)%s(1)])

  status=nf_close(gid(i0))
 END DO

!-----------------------------------------------------------------------
!Create output file
 
  WRITE(*,*) 'Create output'

 !Get grid parameters
 CALL roms_get_type(TRIM(nc_file(3)),fdesign)
 
 !Adapt to new grid
 WRITE(*,*) 'Creating output ',TRIM(nc_file(2))
 fdesign%ntimes=1
 fdesign%grid_size(1:3)=grd(2)%s(1:3)
 CALL roms_create_his_file(TRIM(nc_file(2)),fdesign)

!-------------------------------------------------------------------------
!Interpolate and write


 !Open stream
 status=nf_open(TRIM(nc_file(1)),nf_nowrite,nid(1))
 write(*,*) trim(nc_file(2))
 status=nf_open(TRIM(nc_file(2)),nf_write,nid(2)) 

 WRITE(*,*) 'Interpolating zeta' 
 !Zeta
 ALLOCATE(val2d1(grd(1)%s(1),grd(1)%s(2))); val2d1=dble(0)
 ALLOCATE(val2d2(grd(2)%s(1),grd(2)%s(2))); val2d2=dble(0)

 ncstart=[1,1,itime(1),0]
 nccount=[grd(1)%s(1),grd(1)%s(2),1,0]
 val2d1=RESHAPE(ncread3d(nid(1),'zeta',ncstart(1:3),nccount(1:3)),&
 &[nccount(1:2)])
 WHERE(val2d1.GE.interp_fill_value); val2d1=dble(0); END WHERE
 val2d2=mesh2mesh_tl_interp(&
 &grd(1)%lonr,grd(1)%latr,grd(1)%maskr,val2d1,&
 &grd(2)%lonr,grd(2)%latr,grd(2)%maskr)

 !Write zeta
 write(*,*) 'min/max zeta:',minval(val2d2),maxval(val2d2)
 ncstart=[1,1,1,1]
 nccount=[grd(2)%s(1),grd(2)%s(2),1,1]
 write(*,*) ncstart,nccount
 CALL ncwrite3d(nid(2),'zeta',[ncstart(1),ncstart(2),ncstart(3)],&
 &[nccount(1),nccount(2),nccount(3)],RESHAPE(val2d2,&
 &[nccount(1),nccount(2),nccount(3)]) )
 
 DEALLOCATE(val2d1,val2d2)

!------------------------------------------------------------------------
!Interpolate u-velocities

 WRITE(*,*) 'Interpolating u'

 ALLOCATE(val3d1(grd(1)%s(1)-1,grd(1)%s(2),grd(1)%s(3))); val3d1=dble(0)
 ALLOCATE(val3d2(grd(2)%s(1)-1,grd(2)%s(2),grd(2)%s(3))); val3d2=dble(0)
 ALLOCATE(val2d1(grd(1)%s(1)-1,grd(1)%s(2))); val2d1=dble(0)
 ALLOCATE(val2d2(grd(2)%s(1)-1,grd(2)%s(2))); val2d2=dble(0)

 !Read u-velocities
 ncstart=[1,1,itime(1),1]
 nccount=[grd(1)%s(1)-1,grd(1)%s(2),grd(1)%s(3),1] 
 val3d1=RESHAPE(ncread4d(nid(1),'u',ncstart,nccount),&
 &nccount(1:3))
 WHERE(val3d1.GE.interp_fill_value); val3d1=dble(0); END WHERE
 val3d2=mesh2mesh_tl_interp3d(&
 &grd(1)%lonu,grd(1)%latu,grd(1)%masku,val3d1,&
 &grd(2)%lonu,grd(2)%latu,grd(2)%masku)
  

 !Calculate depth averaged currents
 WRITE(*,*) 'Interpolating ubar'
 val2d2=SUM((&
 &.5*grd(2)%zw(1:grd(2)%s(1)-1,:,2:grd(2)%s(3)+1)+&
 &.5*grd(2)%zw(2:grd(2)%s(1),:,2:grd(2)%s(3)+1)-&
 &.5*grd(2)%zw(1:grd(2)%s(1)-1,:,1:grd(2)%s(3))-&
 &.5*grd(2)%zw(2:grd(2)%s(1),:,1:grd(2)%s(3)))*&
 &val3d2,&
 &3)
 val2d2=val2d2/SUM((&
 &.5*grd(2)%zw(1:grd(2)%s(1)-1,:,2:grd(2)%s(3)+1)+&
 &.5*grd(2)%zw(2:grd(2)%s(1),:,2:grd(2)%s(3)+1)-&
 &.5*grd(2)%zw(1:grd(2)%s(1)-1,:,1:grd(2)%s(3))-&
 &.5*grd(2)%zw(2:grd(2)%s(1),:,1:grd(2)%s(3))),&
 &3)

 write(*,*) 'min/max ubar:',minval(val2d2),maxval(val2d2)
 ncstart=[1,1,1,0]
 nccount=[grd(2)%s(1)-1,grd(2)%s(2),1,0]
 CALL ncwrite3d(nid(2),'ubar',[ncstart(1),ncstart(2),ncstart(3)],&
 &[nccount(1),nccount(2),nccount(3)],&
 &RESHAPE(val2d2,[nccount(1),nccount(2),nccount(3)]))
 
 write(*,*) 'min/max u:',minval(val3d2),maxval(val3d2)
 ncstart=[1,1,1,1]
 nccount=[grd(2)%s(1)-1,grd(2)%s(2),grd(2)%s(3),1]
 CALL ncwrite4d(nid(2),'u',ncstart,nccount,&
 &RESHAPE(val3d2,nccount))

 !deallocate
 DEALLOCATE(val3d1,val3d2,val2d1,val2d2)

!--------------------------------------------------------------------
!Interpolate v-velocities

 WRITE(*,*) 'Interpolating v'

 ALLOCATE(val3d1(grd(1)%s(1),grd(1)%s(2)-1,grd(1)%s(3))); val3d1=dble(0)
 ALLOCATE(val3d2(grd(2)%s(1),grd(2)%s(2)-1,grd(2)%s(3))); val3d2=dble(0)
 ALLOCATE(val2d1(grd(1)%s(1),grd(1)%s(2)-1)); val2d1=dble(0)
 ALLOCATE(val2d2(grd(2)%s(1),grd(2)%s(2)-1)); val2d2=dble(0)

 !Read v-velocities
 ncstart=[1,1,itime(1),1]
 nccount=[grd(1)%s(1),grd(1)%s(2)-1,grd(1)%s(3),1] 
 val3d1=RESHAPE(ncread4d(nid(1),'v',ncstart,nccount),&
 &nccount(1:3))
 WHERE(val3d1.GE.interp_fill_value); val3d1=dble(0); END WHERE
 val3d2=mesh2mesh_tl_interp(&
 &grd(1)%lonv,grd(1)%latv,grd(1)%maskv,val3d1,&
 &grd(2)%lonv,grd(2)%latv,grd(2)%maskv)

 !Calculate depth averaged currents
 WRITE(*,*) 'Interpolating ubar'
 val2d2=SUM((&
 &.5*grd(2)%zw(:,1:grd(2)%s(2)-1,2:grd(2)%s(3)+1)+&
 &.5*grd(2)%zw(:,2:grd(2)%s(2),2:grd(2)%s(3)+1)-&
 &.5*grd(2)%zw(:,1:grd(2)%s(2)-1,1:grd(2)%s(3))-&
 &.5*grd(2)%zw(:,2:grd(2)%s(2),1:grd(2)%s(3)))*&
 &val3d2,&
 &3)
 val2d2=val2d2/SUM((&
 &.5*grd(2)%zw(:,1:grd(2)%s(2)-1,2:grd(2)%s(3)+1)+&
 &.5*grd(2)%zw(:,2:grd(2)%s(2),2:grd(2)%s(3)+1)-&
 &.5*grd(2)%zw(:,1:grd(2)%s(2)-1,1:grd(2)%s(3))-&
 &.5*grd(2)%zw(:,2:grd(2)%s(2),1:grd(2)%s(3))),&
 &3)

 WRITE(*,*) 'min/max vbar:',minval(val2d2),maxval(val2d2)
 ncstart=[1,1,1,0]
 nccount=[grd(2)%s(1),grd(2)%s(2)-1,1,0]
 CALL ncwrite3d(nid(2),'vbar',[ncstart(1),ncstart(2),ncstart(3)],&
 &[nccount(1),nccount(2),nccount(3)],&
 &RESHAPE(val2d2,[nccount(1),nccount(2),nccount(3)]))
 
 WRITE(*,*) 'min/max v:',minval(val3d2),maxval(val3d2)
 ncstart=[1,1,1,1]
 nccount=[grd(2)%s(1),grd(2)%s(2)-1,grd(2)%s(3),1]
 CALL ncwrite4d(nid(2),'v',ncstart,nccount,&
 &RESHAPE(val3d2,nccount))

 DEALLOCATE(val2d1,val2d2,val3d1,val3d2)

!-------------------------------------------------------------------------
!Interpolating temp 
 WRITE(*,*) 'Interpolating temp'

 ALLOCATE(val3d1(grd(1)%s(1),grd(1)%s(2),grd(1)%s(3))); val3d1=dble(0)
 ALLOCATE(val3d2(grd(2)%s(1),grd(2)%s(2),grd(2)%s(3))); val3d2=dble(0)

 !Read temp
 ncstart=[1,1,itime(1),1]
 nccount=[grd(1)%s(1),grd(1)%s(2),grd(1)%s(3),1] 
 val3d1=RESHAPE(ncread4d(nid(1),'temp',ncstart,nccount),&
 &nccount(1:3))
 WHERE(val3d1.GE.interp_fill_value); val3d1=dble(0); END WHERE
 val3d2=mesh2mesh_tl_interp(&
 &grd(1)%lonr,grd(1)%latr,grd(1)%maskr,val3d1,&
 &grd(2)%lonr,grd(2)%latr,grd(2)%maskr)
 !Write temperature
 write(*,*) 'min/max temp:', minval(val3d2),maxval(val3d2)
 ncstart=[1,1,1,1]
 nccount=[grd(2)%s(1),grd(2)%s(2),grd(2)%s(3),1]
 CALL ncwrite4d(nid(2),'temp',ncstart,nccount,&
 &RESHAPE(val3d2,nccount))

!----------------------------------------------------------------------------
!Interpolating salinity

 WRITE(*,*) 'Interpolating salt'
 val3d1=dble(0); val3d2=dble(0); 

 !Read salinity
 ncstart=[1,1,itime(1),1]
 nccount=[grd(1)%s(1),grd(1)%s(2),grd(1)%s(3),1] 
 val3d1=RESHAPE(ncread4d(nid(1),'salt',ncstart,nccount),&
 &nccount(1:3))
 WHERE(val3d1.GE.interp_fill_value); val3d1=dble(0); END WHERE
 val3d2=mesh2mesh_tl_interp(&
 &grd(1)%lonr,grd(1)%latr,grd(1)%maskr,val3d1,&
 &grd(2)%lonr,grd(2)%latr,grd(2)%maskr)

 !Write salinity
 WRITE(*,*) 'min/max salt:',minval(val3d2),maxval(val3d2)
 ncstart=[1,1,1,1]
 nccount=[grd(2)%s(1),grd(2)%s(2),grd(2)%s(3),1]
 CALL ncwrite4d(nid(2),'salt',ncstart,nccount,&
 &RESHAPE(val3d2,nccount))

 DEALLOCATE(val3d1,val3d2)
 
!-----------------------------------------------------------------------
!Interpolating AKt

IF(nf_inq_varid(nid(1),'AKt',i0).EQ.nf_noerr) THEN
 WRITE(*,*) 'Interpolating AKt'
 
 ALLOCATE(val3d1(grd(1)%s(1),grd(1)%s(2),grd(1)%s(3)+1)); val3d1=dble(0)
 ALLOCATE(val3d2(grd(2)%s(1),grd(2)%s(2),grd(2)%s(3)+1)); val3d2=dble(0)

 !Read AKt
 ncstart=[1,1,itime(1),1]
 nccount=[grd(1)%s(1),grd(1)%s(2),grd(1)%s(3)+1,1] 
 val3d1=RESHAPE(ncread4d(nid(1),'AKt',ncstart,nccount),&
 &nccount(1:3))
 WHERE(val3d1.GE.interp_fill_value); val3d1=dble(0); END WHERE
 val3d2=mesh2mesh_tl_interp(&
 &grd(1)%lonr,grd(1)%latr,grd(1)%maskr,val3d1,&
 &grd(2)%lonr,grd(2)%latr,grd(2)%maskr)

 !Write Akt
 WRITE(*,*) 'min/max AKt:',minval(val3d2),maxval(val3d2)
 ncstart=[1,1,1,1]
 nccount=[grd(2)%s(1),grd(2)%s(2),grd(2)%s(3)+1,1]
 CALL ncwrite4d(nid(2),'AKt',ncstart,nccount,&
 &RESHAPE(val3d2,nccount))

 DEALLOCATE(val3d1,val3d2)
END IF

!----------------------------------------------------------------------------

IF(nf_inq_varid(nid(1),'AKt',i0).EQ.nf_noerr) THEN
 WRITE(*,*) 'Interpolating AKv'
 
 ALLOCATE(val3d1(grd(1)%s(1),grd(1)%s(2),grd(1)%s(3)+1)); val3d1=dble(0)
 ALLOCATE(val3d2(grd(2)%s(1),grd(2)%s(2),grd(2)%s(3)+1)); val3d2=dble(0)

 !Read AKv
 ncstart=[1,1,itime(1),1]
 nccount=[grd(1)%s(1),grd(1)%s(2),grd(1)%s(3)+1,1] 
 val3d1=RESHAPE(ncread4d(nid(1),'AKt',ncstart,nccount),&
 &nccount(1:3))
 WHERE(val3d1.GE.interp_fill_value); val3d1=dble(0); END WHERE
 val3d2=mesh2mesh_tl_interp(&
 &grd(1)%lonr,grd(1)%latr,grd(1)%maskr,val3d1,&
 &grd(2)%lonr,grd(2)%latr,grd(2)%maskr)

 !Write Akv
 WRITE(*,*) 'min/max AKv:',minval(val3d2),maxval(val3d2)
 ncstart=[1,1,1,1]
 nccount=[grd(2)%s(1),grd(2)%s(2),grd(2)%s(3)+1,1]
 CALL ncwrite4d(nid(2),'AKv',ncstart,nccount,&
 &RESHAPE(val3d2,nccount))
 
 DEALLOCATE(val3d1,val3d2)
END IF

!--------------------------------------------------------------------------
!Write time

 CALL ncwrite1d(nid(2),'ocean_time',[1],[1],[time])

!---------------------------------------------------------------------------
!Close open streams

 WRITE(*,*) 'Close streams'

 DO i0=1,2
  status=nf_close(nid(i0))
 END DO
 
 WRITE(*,*) 'tl_interp DONE'
END PROGRAM tl_interp
