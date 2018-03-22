#include "cppdefs.h"

#ifdef ADJOINT
#ifdef MPI
      SUBROUTINE ad_get_data
!
! 5/10/2007: ADJ to reading sustrG and svstrG from files, 
! writing the ad_sustrG (overwrite the file)
! 
!------------
! 2/24/2011, MPI OUTPUT:
! - Output arrays from each processor to the single netcdf file, 
! - one at a time
!------------
      USE mod_forces
      USE mod_ocean
      USE mod_parallel  

      implicit none
      include 'mpif.h'
      include 'netcdf.h'

      integer :: start(3),count(3)
      integer :: status,varid,ncid,i,j
      real(8) :: orho0
      real(8), allocatable :: wrk2d(:,:)
      logical::flag_io

! MPI related:
      integer :: k,i1,i2,n1,n2,ierr
      integer :: out_done
      integer,allocatable :: out_done_all(:)

      flag_io=.FALSE.

      IF (iic.eq.ntimes+1) THEN

!------------------------------------------------------------
! iic=ntimes+1: Zero out the fields in the output file of adjoint forcing:
!------------------------------------------------------------

       allocate(out_done_all(0:MaxTile))
       DO k=0,MaxTile
        out_done_all(k)=0
       END DO
       out_done=0  ! -> to be updated to 1 after output's done for the tile

!------------------------------------------------------------------
! - output 1 tile at a time:
!------------------------------------------------------------------
       DO k=0,MaxTile
 !- gather out_done values for each tile and and broadcast them 
        CALL mpi_gather(out_done,1,mpi_integer,
     *                  out_done_all,1,mpi_integer,
     *                  0,mpi_comm_world,ierr)        

        CALL mpi_bcast(out_done_all,MaxTile+1,mpi_integer,
     *                 0,mpi_comm_world,ierr)

        IF (myrank.eq.k .and. out_done_all(k).eq.0) THEN

         status=nf_open(TRIM(frcfile),nf_write,ncid)
 
! - sustr=0: output excludes ghost points (but includes boundary points)
         n1=IendU_E-IstrU_E+1
         n2=JendU_E-JstrU_E+1
 
         allocate(wrk2d(n1,n2)) ! sustr dimensions
         wrk2d(1:n1,1:n2)=0.d0
         start(1)=IstrU_E
         start(2)=JstrU_E+1
         count(1)=n1
         count(2)=n2
         count(3)=1

         DO i=1,nsms
          start(3)=i
          status=nf_inq_varid(ncid,'sustr',varid)
          status=nf_put_vara_double(ncid,varid,start,count,wrk2d)
         END DO

         deallocate(wrk2d)
         
! - svstr=0: output excludes ghost points (but includes boundary points)
         n1=IendV_E-IstrV_E+1
         n2=JendV_E-JstrV_E+1

         allocate(wrk2d(n1,n2)) ! svstr dimensions
         wrk2d(1:n1,1:n2)=0.d0 
         start(1)=IstrV_E+1
         start(2)=JstrV_E
         count(1)=n1
         count(2)=n2
         count(3)=1

         DO i=1,nsms
          start(3)=i
          status=nf_inq_varid(ncid,'svstr',varid)
          status=nf_put_vara_double(ncid,varid,start,count,wrk2d)
         END DO
 
         deallocate(wrk2d)

         status=nf_close(ncid)
         out_done=1

        END IF ! (myrank.eq.k .and. out_done_all(k).eq.0)
       END DO ! k=0:MaxTileNo     

       deallocate(out_done_all)

!-----------------------------------------------------
! iic=ntimes+1: find appropriate time window: define sms_rec, sms_timeG(1:2)
!-----------------------------------------------------

       IF (time.lt.sms_time(1)+0.1d0*dt) THEN
        write (*,*) 'STOP in ad_get_data: time < sms_time(1)+0.1dt'
        STOP
       END IF

       IF (time.gt.sms_time(nsms)+0.1d0*dt) THEN
        write (*,*) 'STOP in ad_get_data: time > sms_time(1)+0.1dt'
        write (*,*) 'time=',time,', sms_time(nsms)=',sms_time(nsms)
        STOP
       END IF
  
       sms_rec=nsms-1
       DO WHILE (time .lt. sms_time(sms_rec))                  
        sms_rec=sms_rec-1 
       END DO
       
       ifrc1=1
       ifrc2=2
       sms_timeG(ifrc1)=sms_time(sms_rec)
       sms_timeG(ifrc2)=sms_time(sms_rec+1)


       IF (flag_io.and.(myrank.eq.0)) then
        write (*,*) ' '
        write (*,*) 'ad_get_data: set sms_timeG'
        write (*,*) '  iic=',iic
        write (*,*) '  sms_timeG(ifrc1)=',sms_timeG(ifrc1)
        write (*,*) '  sms_timeG(ifrc2)=',sms_timeG(ifrc2)
        write (*,*) '  sms_rec=',sms_rec
        write (*,*) ' '
       END IF
       
      END IF ! "iic==ntimes+1"

      IF (iic.gt.1 .and. time-dt.lt.sms_timeG(ifrc1)) THEN

       orho0=1/rho0

!------------------------------------------------------------------
! - output 1 tile at a time:
!------------------------------------------------------------------
       allocate(out_done_all(0:MaxTile))
       DO k=0,MaxTile
        out_done_all(k)=0
       END DO
       out_done=0  ! -> to be updated to 1 after output's done for the tile

       DO k=0,MaxTile
!- gather out_done values for each tile and and broadcast them 
        CALL mpi_gather(out_done,1,mpi_integer,
     *                  out_done_all,1,mpi_integer,
     *                  0,mpi_comm_world,ierr)        
             
        CALL mpi_bcast(out_done_all,MaxTile+1,mpi_integer,
     *                 0,mpi_comm_world,ierr)

        IF (myrank.eq.k .and. out_done_all(k).eq.0) THEN
      
         status=nf_open(TRIM(frcfile),nf_write,ncid)

!-----------------
! output ad_sustr:
!-----------------

#if defined NS_PERIODIC
         CALL ad_exchange_2d(ad_sustrG(LBi:UBi,LBj:UBj,ifrc2),
     *                            LBi,UBi,LBj,UBj,Lm,Mm)
         CALL ad_exchange_2d(ad_svstrG(LBi:UBi,LBj:UBj,ifrc2),
     *                            LBi,UBi,LBj,UBj,Lm,Mm)
#endif 

         n1=IendU_E-IstrU_E+1
         n2=JendU_E-JstrU_E+1
 
         allocate(wrk2d(n1,n2)) ! sustr dimensions
         wrk2d(1:n1,1:n2)=
     *    ad_sustrG(IstrU_E:IendU_E,JstrU_E:JendU_E,ifrc2)*orho0
         start(1)=IstrU_E
         start(2)=JstrU_E+1
         start(3)=sms_rec+1
         count(1)=n1
         count(2)=n2
         count(3)=1
         status=nf_inq_varid(ncid,'sustr',varid)
         status=nf_put_vara_double(ncid,varid,start,count,wrk2d)
         deallocate(wrk2d)
       
!-----------------
! output ad_svstr:
!-----------------
         n1=IendV_E-IstrV_E+1
         n2=JendV_E-JstrV_E+1

         allocate(wrk2d(n1,n2)) ! svstr dimensions
         wrk2d(1:n1,1:n2)=
     *    ad_svstrG(IstrV_E:IendV_E,JstrV_E:JendV_E,ifrc2)*orho0
         start(1)=IstrV_E+1
         start(2)=JstrV_E
         start(3)=sms_rec+1
         count(1)=n1
         count(2)=n2
         count(3)=1
         status=nf_inq_varid(ncid,'svstr',varid)
         status=nf_put_vara_double(ncid,varid,start,count,wrk2d)
         deallocate(wrk2d)
     
         status=nf_close(ncid)
         out_done=1

        END IF ! (myrank.eq.k .and. out_done_all(k).eq.0)
       END DO ! k=0:MaxTileNo     

       deallocate(out_done_all)

       ! clear ifrc2 fields:
       DO j=LBj,UBj
       DO i=LBi,UBi
        ad_sustrG(i,j,ifrc2)=0.d0        
        ad_svstrG(i,j,ifrc2)=0.d0        
       END DO
       END DO

       ! Move the window (rec ifrc1 becomes ifrc2):
       sms_rec=sms_rec-1
       ifrc2=ifrc1
       ifrc1=3-ifrc2
       sms_timeG(ifrc1)=sms_time(sms_rec)

       if(flag_io) then
       write (*,*) ' '
       write (*,*) ' ad_get_data: new sms_timeG'
       write (*,*) '  iic=',iic
       write (*,*) ' sms_timeG(ifrc1)=',sms_timeG(ifrc1)
       write (*,*) ' sms_timeG(ifrc2)=',sms_timeG(ifrc2)
       write (*,*) ' sms_rec=',sms_rec
       write (*,*) ' '
       end if
 
      END IF ! (iic.gt.1 .and. time-dt.lt.sms_timeG(ifrc1))

      IF (iic.eq.1) THEN

       IF(flag_io) THEN
       write (*,*) ' ad_get_data: iic=1 output sms_timeG'
       write (*,*) ' sms_timeG(ifrc1)=',sms_timeG(ifrc1)
       write (*,*) ' sms_timeG(ifrc2)=',sms_timeG(ifrc2)
       write (*,*) ' sms_rec=',sms_rec
       END IF

       orho0=1/rho0

!------------------------------------------------------------------
! - output 1 tile at a time:
!------------------------------------------------------------------
       allocate(out_done_all(0:MaxTile))
       DO k=0,MaxTile
        out_done_all(k)=0
       END DO
       out_done=0  ! -> to be updated to 1 after output's done for the tile

       DO k=0,MaxTile
!- gather out_done values for each tile and and broadcast them 
        CALL mpi_gather(out_done,1,mpi_integer,
     *                  out_done_all,1,mpi_integer,
     *                  0,mpi_comm_world,ierr)        
             
        CALL mpi_bcast(out_done_all,MaxTile+1,mpi_integer,
     *                 0,mpi_comm_world,ierr)

        IF (myrank.eq.k .and. out_done_all(k).eq.0) THEN

         status=nf_open(TRIM(frcfile),nf_write,ncid)

#if defined NS_PERIODIC
         CALL ad_exchange_2d(ad_sustrG(LBi:UBi,LBj:UBj,ifrc1),
     *                            LBi,UBi,LBj,UBj,Lm,Mm)
         CALL ad_exchange_2d(ad_svstrG(LBi:UBi,LBj:UBj,ifrc1),
     *                            LBi,UBi,LBj,UBj,Lm,Mm)
         CALL ad_exchange_2d(ad_sustrG(LBi:UBi,LBj:UBj,ifrc2),
     *                            LBi,UBi,LBj,UBj,Lm,Mm)
         CALL ad_exchange_2d(ad_svstrG(LBi:UBi,LBj:UBj,ifrc2),
     *                            LBi,UBi,LBj,UBj,Lm,Mm)
#endif 

!-----------------
! output ad_sustr (iic=1, rec=sms_rec+1)
!-----------------

         n1=IendU_E-IstrU_E+1
         n2=JendU_E-JstrU_E+1
 
         allocate(wrk2d(n1,n2)) ! sustr dimensions
         wrk2d(1:n1,1:n2)=
     *    ad_sustrG(IstrU_E:IendU_E,JstrU_E:JendU_E,ifrc2)*orho0
         start(1)=IstrU_E
         start(2)=JstrU_E+1
         start(3)=sms_rec+1
         count(1)=n1
         count(2)=n2
         count(3)=1
         status=nf_inq_varid(ncid,'sustr',varid)
         status=nf_put_vara_double(ncid,varid,start,count,wrk2d)
         deallocate(wrk2d)

!-----------------
! output ad_svstr (iic=1, rec=sms_rec+1)
!-----------------
         n1=IendV_E-IstrV_E+1
         n2=JendV_E-JstrV_E+1

         allocate(wrk2d(n1,n2)) ! svstr dimensions
         wrk2d(1:n1,1:n2)=
     *    ad_svstrG(IstrV_E:IendV_E,JstrV_E:JendV_E,ifrc2)*orho0
         start(1)=IstrV_E+1
         start(2)=JstrV_E
         start(3)=sms_rec+1
         count(1)=n1
         count(2)=n2
         count(3)=1
         status=nf_inq_varid(ncid,'svstr',varid)
         status=nf_put_vara_double(ncid,varid,start,count,wrk2d)
         deallocate(wrk2d)

!-----------------
! output ad_sustr (iic=1, rec=sms_rec)
!-----------------
         n1=IendU_E-IstrU_E+1
         n2=JendU_E-JstrU_E+1
 
         allocate(wrk2d(n1,n2)) ! sustr dimensions
         wrk2d(1:n1,1:n2)=
     *    ad_sustrG(IstrU_E:IendU_E,JstrU_E:JendU_E,ifrc1)*orho0
         start(1)=IstrU_E
         start(2)=JstrU_E+1
         start(3)=sms_rec
         count(1)=n1
         count(2)=n2
         count(3)=1
         status=nf_inq_varid(ncid,'sustr',varid)
         status=nf_put_vara_double(ncid,varid,start,count,wrk2d)
         deallocate(wrk2d)
       
!-----------------
! output ad_svstr:
!-----------------
         n1=IendV_E-IstrV_E+1
         n2=JendV_E-JstrV_E+1

         allocate(wrk2d(n1,n2)) ! svstr dimensions
         wrk2d(1:n1,1:n2)=
     *    ad_svstrG(IstrV_E:IendV_E,JstrV_E:JendV_E,ifrc1)*orho0
         start(1)=IstrV_E+1
         start(2)=JstrV_E
         start(3)=sms_rec
         count(1)=n1
         count(2)=n2
         count(3)=1
         status=nf_inq_varid(ncid,'svstr',varid)
         status=nf_put_vara_double(ncid,varid,start,count,wrk2d)
         deallocate(wrk2d)
     
         status=nf_close(ncid)
         out_done=1

        END IF ! (myrank.eq.k .and. out_done_all(k).eq.0)
       END DO ! k=0:MaxTileNo     

       deallocate(out_done_all)

      END IF ! (iic.eq.1)

      RETURN
      END

! END OF THE MPI PIECE

!==========================================================================
!==========================================================================
#else /* serial */

      SUBROUTINE ad_get_data
      USE mod_forces
      USE mod_ocean
      USE mod_nf_param

      implicit none
      integer :: start(3),count(3)
      integer :: status,varid,ncid,i,j
      real(8) :: orho0
      real(8), allocatable :: wrk2d(:,:)
      logical :: flag_io

      flag_io=.false.
      
!
! 5/10/2007: ADJ to reading sustrG and svstrG from files, 
! writing the ad_sustrG (overwrite the file)
! 
 
      IF (iic.eq.ntimes+1) THEN

!------------------------------------------------------------
! iic=ntimes+1: Zero out the fields in the output file of adjoint forcing:
!------------------------------------------------------------
       status=nf_open(TRIM(frcfile),nf_write,ncid)

       start(1)=1
       start(2)=1
       count(3)=1

! - sustr=0:       
       allocate(wrk2d(L,Mp)) ! sustr dimensions
       wrk2d=0.d0

       count(1)=L
       count(2)=Mp

       DO i=1,nsms
        start(3)=i
        status=nf_inq_varid(ncid,'sustr',varid)
        status=nf_put_vara_double(ncid,varid,start,count,wrk2d)
       END DO

       deallocate(wrk2d)
         
! - svstr=0:
       allocate(wrk2d(Lp,M)) ! svstr dimensions
       wrk2d=0.d0

       count(1)=Lp
       count(2)=M

       DO i=1,nsms
        start(3)=i
        status=nf_inq_varid(ncid,'svstr',varid)
        status=nf_put_vara_double(ncid,varid,start,count,wrk2d)
       END DO

      status=nf_close(ncid)

      deallocate(wrk2d)

!-----------------------------------------------------
! iic=ntimes+1: find appropriate time window: define sms_rec, sms_timeG(1:2)
!-----------------------------------------------------

       IF (time.lt.sms_time(1)+0.1d0*dt) THEN
        write (*,*) 'STOP in ad_get_data: time < sms_time(1)+0.1dt'
        STOP
       END IF

       IF (time.gt.sms_time(nsms)+0.1d0*dt) THEN
        write (*,*) 'STOP in ad_get_data: time > sms_time(1)+0.1dt'
        write (*,*) 'time=',time,', sms_time(nsms)=',sms_time(nsms)
        STOP
       END IF
  
       sms_rec=nsms-1
       DO WHILE (time .lt. sms_time(sms_rec))                  
        sms_rec=sms_rec-1 
       END DO
       
       ifrc1=1
       ifrc2=2
       sms_timeG(ifrc1)=sms_time(sms_rec)
       sms_timeG(ifrc2)=sms_time(sms_rec+1)

       IF(flag_io) THEN
       write (*,*) ' '
       write (*,*) 'ad_get_data: set sms_timeG'
       write (*,*) '  iic=',iic
       write (*,*) '  sms_timeG(ifrc1)=',sms_timeG(ifrc1)
       write (*,*) '  sms_timeG(ifrc2)=',sms_timeG(ifrc2)
       write (*,*) '  sms_rec=',sms_rec
       write (*,*) ' '
       END IF
       
      END IF ! "iic==ntimes+1"

      IF (iic.gt.1 .and. time-dt.lt.sms_timeG(ifrc1)) THEN

       orho0=1/rho0
      
       status=nf_open(TRIM(frcfile),nf_write,ncid)

!-----------------
! output ad_sustr:
!-----------------

#if defined NS_PERIODIC
       CALL ad_exchange_2d(ad_sustrG(LBi:UBi,LBj:UBj,ifrc2),
     *                            LBi,UBi,LBj,UBj,Lm,Mm)
       CALL ad_exchange_2d(ad_svstrG(LBi:UBi,LBj:UBj,ifrc2),
     *                            LBi,UBi,LBj,UBj,Lm,Mm)
#endif 


       allocate(wrk2d(L,Mp))
       wrk2d(1:L,1:Mp)=ad_sustrG(1:L,0:M,ifrc2)*orho0
       start(1)=1
       start(2)=1
       start(3)=sms_rec+1
       count(1)=L
       count(2)=Mp
       count(3)=1
       status=nf_inq_varid(ncid,'sustr',varid)
       status=nf_put_vara_double(ncid,varid,start,count,wrk2d)
       deallocate(wrk2d)
       
!-----------------
! output ad_svstr:
!-----------------
       allocate(wrk2d(Lp,M))
       wrk2d(1:Lp,1:M)=ad_svstrG(0:L,1:M,ifrc2)*orho0
       start(1)=1
       start(2)=1
       start(3)=sms_rec+1
       count(1)=Lp
       count(2)=M
       count(3)=1
       status=nf_inq_varid(ncid,'svstr',varid)
       status=nf_put_vara_double(ncid,varid,start,count,wrk2d)
       deallocate(wrk2d)
     
       status=nf_close(ncid)

       ! clear ifrc2 fields:
       DO j=LBj,UBj
       DO i=LBi,UBi
        ad_sustrG(i,j,ifrc2)=0.d0        
        ad_svstrG(i,j,ifrc2)=0.d0        
       END DO
       END DO

       ! Move the window (rec ifrc1 becomes ifrc2):
       sms_rec=sms_rec-1
       ifrc2=ifrc1
       ifrc1=3-ifrc2
       sms_timeG(ifrc1)=sms_time(sms_rec)

       IF(flag_io) THEN
       write (*,*) ' '
       write (*,*) ' ad_get_data: new sms_timeG'
       write (*,*) '  iic=',iic
       write (*,*) ' sms_timeG(ifrc1)=',sms_timeG(ifrc1)
       write (*,*) ' sms_timeG(ifrc2)=',sms_timeG(ifrc2)
       write (*,*) ' sms_rec=',sms_rec
       write (*,*) ' '
       END IF
 
      END IF

      IF (iic.eq.1) THEN

       IF(flag_io) THEN
       write (*,*) ' ad_get_data: iic=1 output sms_timeG'
       write (*,*) ' sms_timeG(ifrc1)=',sms_timeG(ifrc1)
       write (*,*) ' sms_timeG(ifrc2)=',sms_timeG(ifrc2)
       write (*,*) ' sms_rec=',sms_rec
       END IF
       
       orho0=1/rho0

      !output both fields
       status=nf_open(TRIM(frcfile),nf_write,ncid)

#if defined NS_PERIODIC
       CALL ad_exchange_2d(ad_sustrG(LBi:UBi,LBj:UBj,ifrc1),
     *                            LBi,UBi,LBj,UBj,Lm,Mm)
       CALL ad_exchange_2d(ad_svstrG(LBi:UBi,LBj:UBj,ifrc1),
     *                            LBi,UBi,LBj,UBj,Lm,Mm)
       CALL ad_exchange_2d(ad_sustrG(LBi:UBi,LBj:UBj,ifrc2),
     *                            LBi,UBi,LBj,UBj,Lm,Mm)
       CALL ad_exchange_2d(ad_svstrG(LBi:UBi,LBj:UBj,ifrc2),
     *                            LBi,UBi,LBj,UBj,Lm,Mm)
#endif 

!-----------------
! output ad_sustr (iic=1, rec=sms_rec+1)
!-----------------


       allocate(wrk2d(L,Mp))
       wrk2d(1:L,1:Mp)=ad_sustrG(1:L,0:M,ifrc2)*orho0
       start(1)=1
       start(2)=1
       start(3)=sms_rec+1
       count(1)=L
       count(2)=Mp
       count(3)=1
       status=nf_inq_varid(ncid,'sustr',varid)
       status=nf_put_vara_double(ncid,varid,start,count,wrk2d)
       deallocate(wrk2d)
       
!-----------------
! output ad_svstr (iic=1, rec=sms_rec+1)
!-----------------
       allocate(wrk2d(Lp,M))
       wrk2d(1:Lp,1:M)=ad_svstrG(0:L,1:M,ifrc2)*orho0
       start(1)=1
       start(2)=1
       start(3)=sms_rec+1
       count(1)=Lp
       count(2)=M
       count(3)=1
       status=nf_inq_varid(ncid,'svstr',varid)
       status=nf_put_vara_double(ncid,varid,start,count,wrk2d)
       deallocate(wrk2d)

!-----------------
! output ad_sustr (iic=1, rec=sms_rec)
!-----------------
       allocate(wrk2d(L,Mp))
       wrk2d(1:L,1:Mp)=ad_sustrG(1:L,0:M,ifrc1)*orho0
       start(1)=1
       start(2)=1
       start(3)=sms_rec
       count(1)=L
       count(2)=Mp
       count(3)=1
       status=nf_inq_varid(ncid,'sustr',varid)
       status=nf_put_vara_double(ncid,varid,start,count,wrk2d)
       deallocate(wrk2d)
       
!-----------------
! output ad_svstr:
!-----------------
       allocate(wrk2d(Lp,M))
       wrk2d(1:Lp,1:M)=ad_svstrG(0:L,1:M,ifrc1)*orho0
       start(1)=1
       start(2)=1
       start(3)=sms_rec
       count(1)=Lp
       count(2)=M
       count(3)=1
       status=nf_inq_varid(ncid,'svstr',varid)
       status=nf_put_vara_double(ncid,varid,start,count,wrk2d)
       deallocate(wrk2d)
     
       status=nf_close(ncid)

      END IF        

      RETURN
      END


# endif /* MPI or serial */

#endif /* ADJOINT */
