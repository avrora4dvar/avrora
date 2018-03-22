#include "cppdefs.h"

#ifdef TANGENT
      SUBROUTINE tl_get_data
      USE mod_forces
      USE mod_ocean
#ifdef MPI
      USE mod_parallel
#endif
      
      implicit none
      include 'netcdf.h'
      integer :: status,fid,varid,dimid,i,j
      real(8) :: orho0
#ifdef MPI
      include 'mpif.h'
      integer :: n1,n2,i1,i2
#endif

!
! Finding the time window, wind stress:
!
      IF (iic.eq.1) THEN

       IF (time.lt.sms_time(1)-0.1d0*dt) THEN
        write (*,*) 
     *   'STOP in tl_get_data, iic=1: time < sms_time(1)-0.1dt'
        STOP
       END IF

       IF (time.gt.sms_time(nsms)-0.1d0) THEN
        write (*,*) 
     *   'STOP in tl_get_data, iic=1: time > sms_time(nsms)-0.1dt'
        STOP
       END IF

! Now, since time>=sms_time(1)-0.1dt
       IF (time.lt.sms_time(1)) THEN 
        ! think time==sms_time(1) with machine accur.
        sms_rec=1
       ELSE
        sms_rec=1
        DO WHILE (time.ge.sms_time(sms_rec))                  
         sms_rec=sms_rec+1
        END DO
        sms_rec=sms_rec-1        
       END IF
   
       ifrc1=1
       ifrc2=2
       sms_timeG(ifrc1)=sms_time(sms_rec)
       sms_timeG(ifrc2)=sms_time(sms_rec+1)


#ifdef MPI
! - in u locations:
      i1=max(1,LBi)  ! index of 1st element on roms grid
      n1=UBi-i1+1
      n2=UBj-LBj+1

      CALL get_2dfld_tile(frcfile,'sustr',sms_rec,i1,LBj+1,
     *                    n1,n2,tl_sustrG(i1:UBi,LBj:UBj,ifrc1))
      CALL get_2dfld_tile(frcfile,'sustr',sms_rec+1,i1,LBj+1,
     *                    n1,n2,tl_sustrG(i1:UBi,LBj:UBj,ifrc2))

! in v-locations:
      i2=max(1,LBj)
      n1=UBi-LBi+1
      n2=UBj-i2+1

      CALL get_2dfld_tile(frcfile,'svstr',sms_rec,LBi+1,i2,n1,n2,
     *                    tl_svstrG(LBi:UBi,i2:UBj,ifrc1))
      CALL get_2dfld_tile(frcfile,'svstr',sms_rec+1,LBi+1,i2,n1,n2,
     *                    tl_svstrG(LBi:UBi,i2:UBj,ifrc2))

#else /* not MPI, serial */
       CALL get_2dfld(frcfile,'sustr',sms_rec,  L,Mp,
     *                tl_sustrG(1:L,0:M,ifrc1))
       CALL get_2dfld(frcfile,'sustr',sms_rec+1,L,Mp,
     *                tl_sustrG(1:L,0:M,ifrc2))

       CALL get_2dfld(frcfile,'svstr',sms_rec,  Lp,M,
     *                tl_svstrG(0:L,1:M,ifrc1))
       CALL get_2dfld(frcfile,'svstr',sms_rec+1,Lp,M,
     *                tl_svstrG(0:L,1:M,ifrc2))

#endif /* MPI or else */

       orho0=1/rho0

       DO j=LBj,UBj
       DO i=LBi,UBi       
        tl_sustrG(i,j,ifrc1)=tl_sustrG(i,j,ifrc1)*orho0
        tl_sustrG(i,j,ifrc2)=tl_sustrG(i,j,ifrc2)*orho0
        tl_svstrG(i,j,ifrc1)=tl_svstrG(i,j,ifrc1)*orho0
        tl_svstrG(i,j,ifrc2)=tl_svstrG(i,j,ifrc2)*orho0
       END DO
       END DO

#ifndef MPI
#if defined NS_PERIODIC
       CALL exchange_2d(tl_sustrG(LBi:UBi,LBj:UBj,ifrc1),
     *                            LBi,UBi,LBj,UBj,Lm,Mm)
       CALL exchange_2d(tl_sustrG(LBi:UBi,LBj:UBj,ifrc2),
     *                            LBi,UBi,LBj,UBj,Lm,Mm)
       CALL exchange_2d(tl_svstrG(LBi:UBi,LBj:UBj,ifrc1),
     *                            LBi,UBi,LBj,UBj,Lm,Mm)
       CALL exchange_2d(tl_svstrG(LBi:UBi,LBj:UBj,ifrc2),
     *                            LBi,UBi,LBj,UBj,Lm,Mm)
#endif 
#endif /* MPI not defined */

      END IF

!-----------
! Updating:
!-----------
      IF (time.ge.sms_timeG(ifrc2)) THEN

       IF (sms_rec.lt.nsms) THEN

        sms_rec=sms_rec+1

        ifrc1=ifrc2
        ifrc2=3-ifrc1
       
        sms_timeG(ifrc2)=sms_time(sms_rec+1)


#ifdef MPI
! - in u locations:
      i1=max(1,LBi)  ! index of 1st element on roms grid
      n1=UBi-i1+1
      n2=UBj-LBj+1
      CALL get_2dfld_tile(frcfile,'sustr',sms_rec+1,i1,LBj+1,
     *                    n1,n2,tl_sustrG(i1:UBi,LBj:UBj,ifrc2))

! in v-locations:
      i2=max(1,LBj)
      n1=UBi-LBi+1
      n2=UBj-i2+1
      CALL get_2dfld_tile(frcfile,'svstr',sms_rec+1,LBi+1,i2,n1,n2,
     *                    tl_svstrG(LBi:UBi,i2:UBj,ifrc2))

#else
        CALL get_2dfld(frcfile,'sustr',sms_rec+1,L,Mp,
     *                 tl_sustrG(1:L,0:M,ifrc2))
        CALL get_2dfld(frcfile,'svstr',sms_rec+1,Lp,M,
     *                 tl_svstrG(0:L,1:M,ifrc2))
#endif /* MPI or else */

        orho0=1/rho0
        DO j=LBj,UBj
        DO i=LBi,UBi       
         tl_sustrG(i,j,ifrc2)=tl_sustrG(i,j,ifrc2)*orho0
         tl_svstrG(i,j,ifrc2)=tl_svstrG(i,j,ifrc2)*orho0
        END DO
        END DO

#ifndef MPI
#if defined NS_PERIODIC
       CALL exchange_2d(tl_sustrG(LBi:UBi,LBj:UBj,ifrc2),
     *                            LBi,UBi,LBj,UBj,Lm,Mm)
       CALL exchange_2d(tl_svstrG(LBi:UBi,LBj:UBj,ifrc2),
     *                            LBi,UBi,LBj,UBj,Lm,Mm)
#endif 
#endif /* MPI not defined */

       ELSE 
        ! end of record, sms_rec==nsms, allow for round-up error
        IF (time.gt.sms_time(nsms)+0.1d0*dt) THEN
         write (*,*) 'STOP in get_data, since time>sms_time(last)'
         STOP
        END IF                
       END IF
 
      END IF

      RETURN
      END
#endif /* TANGENT */
