      SUBROUTINE ad_input_ini(fname,
     *                   L,M,N,
     *                   time,ad_zeta,
     *                   ad_u,ad_v,ad_T,ad_S)

      implicit none
      include 'netcdf.h'
      character (len=80) :: fname
      integer :: status,ncid,recnum,varid
      integer, dimension(4) :: dimids,start,count

      integer :: L,M,N,Lp,Mp,Np
      
      real(8) :: time
      real(8) :: ad_zeta(0:L,0:M)
      real(8) :: ad_u(0:L,0:M,N)
      real(8) :: ad_v(0:L,0:M,N)
      real(8) :: ad_T(0:L,0:M,N)
      real(8) :: ad_S(0:L,0:M,N)

      real(8), allocatable :: wrk2d(:,:)
      real(8), allocatable :: wrk3d(:,:,:)

      dimids=0

      Lp=L+1
      Mp=M+1

      status=nf_open(TRIM(fname),nf_nowrite,ncid)
!--------------
! read ocean_time:
!--------------
      start(1)=1
      count(1)=1

      status=nf_inq_varid(ncid,'ocean_time',varid)
      status=nf_get_vara_double(ncid,varid,start,count,time)

!--------------
! read zeta:
!--------------
      allocate(wrk2d(Lp,Mp))
      start(1)=1
      start(2)=1
      start(3)=1
      count(1)=Lp
      count(2)=Mp
      count(3)=1
      status=nf_inq_varid(ncid,'zeta',varid)
      status=nf_get_vara_double(ncid,varid,start,count,wrk2d)

! tl:      wrk2d(1:Lp,1:Mp)=tl_zeta(0:L,0:M)
      ad_zeta(0:L,0:M)=wrk2d(1:Lp,1:Mp)    

      deallocate(wrk2d)

!--------------
! read ad_u:
!--------------
      allocate(wrk3d(L,Mp,N))
      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=1
      count(1)=L
      count(2)=Mp
      count(3)=N
      count(4)=1

      status=nf_inq_varid(ncid,'u',varid)
      status=nf_get_vara_double(ncid,varid,start,count,wrk3d)

! tl:      wrk3d(1:L,1:Mp,1:N)=tl_u(1:L,0:M,1:N)
      ad_u(1:L,0:M,1:N)=wrk3d(1:L,1:Mp,1:N)

      deallocate(wrk3d)

!--------------
! read ad_v:
!--------------
      allocate(wrk3d(Lp,M,N))
      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=1
      count(1)=Lp
      count(2)=M
      count(3)=N
      count(4)=1

      status=nf_inq_varid(ncid,'v',varid)
      status=nf_get_vara_double(ncid,varid,start,count,wrk3d)

! tl:      wrk3d(1:Lp,1:M,1:N)=tl_v(0:L,1:M,1:N)
      ad_v(0:L,1:M,1:N)=wrk3d(1:Lp,1:M,1:N)

      deallocate(wrk3d)

!--------------
! output T, S:
!--------------
      allocate(wrk3d(Lp,Mp,N))
      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=1
      count(1)=Lp
      count(2)=Mp
      count(3)=N
      count(4)=1

      status=nf_inq_varid(ncid,'temp',varid)
      status=nf_get_vara_double(ncid,varid,start,count,wrk3d)
! tl:      wrk3d(1:Lp,1:Mp,1:N)=tl_T(0:L,0:M,1:N)
      ad_T(0:L,0:M,1:N)=wrk3d(1:Lp,1:Mp,1:N)

      status=nf_inq_varid(ncid,'salt',varid)
      status=nf_get_vara_double(ncid,varid,start,count,wrk3d)
! tl:      wrk3d(1:Lp,1:Mp,1:N)=tl_S(0:L,0:M,1:N)
      ad_S(0:L,0:M,1:N)=wrk3d(1:Lp,1:Mp,1:N)

      deallocate(wrk3d)

      status=nf_close(ncid)

 
      RETURN
      END
