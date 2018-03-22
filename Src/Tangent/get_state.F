#include "cppdefs.h"

#if defined TANGENT
      SUBROUTINE get_state

! AK, 5/7/2007: finds two appropriate snapshots of background (fwd)
! state, to be used for linear interpolatation between them
! Fields are initialized at iic=0

      USE mod_ocean     ! in part., keeps iic

      implicit none
      include 'netcdf.h'
      real(8) :: t0
#ifdef MPI
      integer :: i1,i2,n1,n2
#endif

      t0=time+dt

! Finding the time window (during first call -- ini fwd state):

      IF (iic.eq.0) THEN

       IF (t0.lt.fwd_time(1)-0.1d0*dt) THEN
        write (*,*)'STOP in get_state, iic=0: t0 < fwd_time(1)-0.1dt'
        STOP
       END IF

       IF (t0.gt.fwd_time(nfwd)-0.1d0) THEN
        write (*,*) 'STOP in get_state, iic=0: t0 > fwd_time(end)-0.1dt'
        STOP
       END IF

! Now, since t0>=fwd_time(1)-0.1dt
       IF (t0.lt.fwd_time(1)) THEN 
        ! think t0==fwd_time(1) with machine accur.
        fwd_rec=1
       ELSE
        fwd_rec=1
        DO WHILE (t0.ge.fwd_time(fwd_rec))                  
         fwd_rec=fwd_rec+1
        END DO
        fwd_rec=fwd_rec-1        
       END IF

       ifwd1=1
       ifwd2=2

       fwd_timeG(ifwd1)=fwd_time(fwd_rec)
       fwd_timeG(ifwd2)=fwd_time(fwd_rec+1)


#ifdef MPI
!!! read copy 1

! - in rho locations 
      n1=UBi-LBi+1
      n2=UBj-LBj+1

      CALL get_2dfld_tile(fwdfile,'zeta',fwd_rec,LBi+1,LBj+1,
     *                    n1,n2,zetaG(LBi:UBi,LBj:UBj,1))
      CALL get_3dfld_tile(fwdfile,'temp',fwd_rec,
     *                    LBi+1,LBj+1,1,n1,n2,N,     ! start (nc file), count
     *                    tG(LBi:UBi,LBj:UBj,1:N,1,1))  
      CALL get_3dfld_tile(fwdfile,'salt',fwd_rec,
     *                    LBi+1,LBj+1,1,n1,n2,N,     ! start (nc file), count
     *                    tG(LBi:UBi,LBj:UBj,1:N,1,2))  

#ifndef ANA_FWD_MIX
      CALL get_3dfld_tile(fwdfile,'AKv',fwd_rec,
     *                    LBi+1,LBj+1,1,n1,n2,Np,     ! start (nc file), count
     *                    AkvG(LBi:UBi,LBj:UBj,0:N,1))  
      CALL get_3dfld_tile(fwdfile,'AKt',fwd_rec,
     *                    LBi+1,LBj+1,1,n1,n2,Np,     ! start (nc file), count
     *                    AktG(LBi:UBi,LBj:UBj,0:N,1))  
#endif

! - in u locations:
      i1=max(1,LBi)  ! index of 1st element on roms grid
      n1=UBi-i1+1
      n2=UBj-LBj+1
      CALL get_2dfld_tile(fwdfile,'ubar',fwd_rec,i1,LBj+1,
     *                    n1,n2,ubarG(i1:UBi,LBj:UBj,1))
      CALL get_3dfld_tile(fwdfile,'u',fwd_rec,
     *                    i1,LBj+1,1,n1,n2,N,     ! start (nc file), count
     *                    uG(i1:UBi,LBj:UBj,1:N,1))  

! in v-locations:
      i2=max(1,LBj)
      n1=UBi-LBi+1
      n2=UBj-i2+1
      CALL get_2dfld_tile(fwdfile,'vbar',fwd_rec,LBi+1,i2,
     *                    n1,n2,vbarG(LBi:UBi,i2:UBj,1))
      CALL get_3dfld_tile(fwdfile,'v',fwd_rec,
     *                    LBi+1,i2,1,n1,n2,N,     ! start (nc file), count
     *                    vG(LBi:UBi,i2:UBj,1:N,1))

!!! read copy 2

! - in rho locations 
      n1=UBi-LBi+1
      n2=UBj-LBj+1

      CALL get_2dfld_tile(fwdfile,'zeta',fwd_rec+1,LBi+1,LBj+1,
     *                    n1,n2,zetaG(LBi:UBi,LBj:UBj,2))
      CALL get_3dfld_tile(fwdfile,'temp',fwd_rec+1,
     *                    LBi+1,LBj+1,1,n1,n2,N,     ! start (nc file), count
     *                    tG(LBi:UBi,LBj:UBj,1:N,2,1))  
      CALL get_3dfld_tile(fwdfile,'salt',fwd_rec+1,
     *                    LBi+1,LBj+1,1,n1,n2,N,     ! start (nc file), count
     *                    tG(LBi:UBi,LBj:UBj,1:N,2,2))  

#ifndef ANA_FWD_MIX
      CALL get_3dfld_tile(fwdfile,'AKv',fwd_rec+1,
     *                    LBi+1,LBj+1,1,n1,n2,Np,     ! start (nc file), count
     *                    AkvG(LBi:UBi,LBj:UBj,0:N,2))  
      CALL get_3dfld_tile(fwdfile,'AKt',fwd_rec+1,
     *                    LBi+1,LBj+1,1,n1,n2,Np,     ! start (nc file), count
     *                    AktG(LBi:UBi,LBj:UBj,0:N,2))  
#endif

! - in u locations:
      i1=max(1,LBi)  ! index of 1st element on roms grid
      n1=UBi-i1+1
      n2=UBj-LBj+1
      CALL get_2dfld_tile(fwdfile,'ubar',fwd_rec+1,i1,LBj+1,
     *                    n1,n2,ubarG(i1:UBi,LBj:UBj,2))
      CALL get_3dfld_tile(fwdfile,'u',fwd_rec+1,
     *                    i1,LBj+1,1,n1,n2,N,     ! start (nc file), count
     *                    uG(i1:UBi,LBj:UBj,1:N,2))  

! in v-locations:
      i2=max(1,LBj)
      n1=UBi-LBi+1
      n2=UBj-i2+1
      CALL get_2dfld_tile(fwdfile,'vbar',fwd_rec+1,LBi+1,i2,
     *                    n1,n2,vbarG(LBi:UBi,i2:UBj,2))
      CALL get_3dfld_tile(fwdfile,'v',fwd_rec+1,
     *                    LBi+1,i2,1,n1,n2,N,     ! start (nc file), count
     *                    vG(LBi:UBi,i2:UBj,1:N,2))

#else /* serial */

!!! read copy 1
       CALL get_2dfld(fwdfile,'zeta',fwd_rec,Lp,Mp,zetaG(0:L,0:M,1))
       CALL get_2dfld(fwdfile,'ubar',fwd_rec,L ,Mp,ubarG(1:L,0:M,1))
       CALL get_2dfld(fwdfile,'vbar',fwd_rec,Lp,M ,vbarG(0:L,1:M,1))
       CALL get_3dfld(fwdfile,'u'   ,fwd_rec,L ,Mp,N,uG(1:L,0:M,1:N,1))
       CALL get_3dfld(fwdfile,'v'   ,fwd_rec,Lp,M ,N,vG(0:L,1:M,1:N,1))
       CALL get_3dfld(fwdfile,'temp',fwd_rec,Lp,Mp,N,
     *                                             tG(0:L,0:M,1:N,1,1))
       CALL get_3dfld(fwdfile,'salt',fwd_rec,Lp,Mp,N,
     *                                             tG(0:L,0:M,1:N,1,2))

#ifndef ANA_FWD_MIX
       CALL get_3dfld(fwdfile,'AKv'  ,fwd_rec,Lp,Mp,Np,
     *                                            AkvG(0:L,0:M,0:N,1))
       CALL get_3dfld(fwdfile,'AKt'  ,fwd_rec,Lp,Mp,Np,
     *                                            AktG(0:L,0:M,0:N,1))
#endif
       
!!! read copy 2

       CALL get_2dfld(fwdfile,'zeta',fwd_rec+1,Lp,Mp,zetaG(0:L,0:M,2))
       CALL get_2dfld(fwdfile,'ubar',fwd_rec+1,L ,Mp,ubarG(1:L,0:M,2))
       CALL get_2dfld(fwdfile,'vbar',fwd_rec+1,Lp,M ,vbarG(0:L,1:M,2))
       CALL get_3dfld(fwdfile,'u'   ,fwd_rec+1,L ,Mp,N,
     *                                             uG(1:L,0:M,1:N,2))
       CALL get_3dfld(fwdfile,'v'   ,fwd_rec+1,Lp,M ,N,
     *                                             vG(0:L,1:M,1:N,2))
       CALL get_3dfld(fwdfile,'temp',fwd_rec+1,Lp,Mp,N,
     *                                             tG(0:L,0:M,1:N,2,1))
       CALL get_3dfld(fwdfile,'salt',fwd_rec+1,Lp,Mp,N,
     *                                             tG(0:L,0:M,1:N,2,2))

       CALL get_3dfld(fwdfile,'v'   ,fwd_rec+1,Lp,M ,N,
     *                                             vG(0:L,1:M,1:N,2))

#ifndef ANA_FWD_MIX
       CALL get_3dfld(fwdfile,'AKv'  ,fwd_rec+1,Lp,Mp,Np,
     *                                            AkvG(0:L,0:M,0:N,2))
       CALL get_3dfld(fwdfile,'AKt'  ,fwd_rec+1,Lp,Mp,Np,
     *                                            AktG(0:L,0:M,0:N,2))
#endif

#endif /* MPI or serial */



#if defined NS_PERIODIC
#ifdef MPI
      STOP('periodic BC not coded in MPI')
#else /* serial */
      CALL exchange_2d(zetaG(:,:,ifwd1),LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(ubarG(:,:,ifwd1),LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(vbarG(:,:,ifwd1),LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_3d(uG(:,:,:,ifwd1),LBi,UBi,LBj,UBj,1,N,Lm,Mm)
      CALL exchange_3d(vG(:,:,:,ifwd1),LBi,UBi,LBj,UBj,1,N,Lm,Mm)
      CALL exchange_3d(tG(:,:,:,ifwd1,1),LBi,UBi,LBj,UBj,1,N,Lm,Mm)
      CALL exchange_3d(tG(:,:,:,ifwd1,2),LBi,UBi,LBj,UBj,1,N,Lm,Mm)

      CALL exchange_2d(zetaG(:,:,ifwd2),LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(ubarG(:,:,ifwd2),LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(vbarG(:,:,ifwd2),LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_3d(uG(:,:,:,ifwd2),LBi,UBi,LBj,UBj,1,N,Lm,Mm)
      CALL exchange_3d(vG(:,:,:,ifwd2),LBi,UBi,LBj,UBj,1,N,Lm,Mm)
      CALL exchange_3d(tG(:,:,:,ifwd2,1),LBi,UBi,LBj,UBj,1,N,Lm,Mm)
      CALL exchange_3d(tG(:,:,:,ifwd2,2),LBi,UBi,LBj,UBj,1,N,Lm,Mm)

# ifndef ANA_FWD_MIX
      CALL exchange_3d(AkvG(:,:,:,ifwd1),LBi,UBi,LBj,UBj,0,N,Lm,Mm)
      CALL exchange_3d(AktG(:,:,:,ifwd1),LBi,UBi,LBj,UBj,0,N,Lm,Mm)
      CALL exchange_3d(AkvG(:,:,:,ifwd2),LBi,UBi,LBj,UBj,0,N,Lm,Mm)
      CALL exchange_3d(AktG(:,:,:,ifwd2),LBi,UBi,LBj,UBj,0,N,Lm,Mm)
# endif

#endif /* MPI or serial */

#endif /* defined NS_PERIODIC */

      END IF ! (iic.eq.0)
      
!-----------
! Updating:
!-----------
      IF (t0.gt.fwd_timeG(ifwd2)) THEN

       IF (fwd_rec.lt.nfwd) THEN

        fwd_rec=fwd_rec+1

! switch indices (instead of copying arrays)
        ifwd1=3-ifwd1
        ifwd2=3-ifwd2

        fwd_timeG(ifwd2)=fwd_time(fwd_rec+1)

#ifdef MPI

! - in rho locations 
        n1=UBi-LBi+1
        n2=UBj-LBj+1

        CALL get_2dfld_tile(fwdfile,'zeta',fwd_rec+1,LBi+1,LBj+1,
     *                    n1,n2,zetaG(LBi:UBi,LBj:UBj,ifwd2))
        CALL get_3dfld_tile(fwdfile,'temp',fwd_rec+1,
     *                    LBi+1,LBj+1,1,n1,n2,N,     ! start (nc file), count
     *                    tG(LBi:UBi,LBj:UBj,1:N,ifwd2,1))  
        CALL get_3dfld_tile(fwdfile,'salt',fwd_rec+1,
     *                    LBi+1,LBj+1,1,n1,n2,N,     ! start (nc file), count
     *                    tG(LBi:UBi,LBj:UBj,1:N,ifwd2,2))  

#ifndef ANA_FWD_MIX
        CALL get_3dfld_tile(fwdfile,'AKv',fwd_rec+1,
     *                    LBi+1,LBj+1,1,n1,n2,Np,     ! start (nc file), count
     *                    AkvG(LBi:UBi,LBj:UBj,0:N,ifwd2))  
        CALL get_3dfld_tile(fwdfile,'AKt',fwd_rec+1,
     *                    LBi+1,LBj+1,1,n1,n2,Np,     ! start (nc file), count
     *                    AktG(LBi:UBi,LBj:UBj,0:N,ifwd2))  
#endif

! - in u locations:
        i1=max(1,LBi)  ! index of 1st element on roms grid
        n1=UBi-i1+1
        n2=UBj-LBj+1
        CALL get_2dfld_tile(fwdfile,'ubar',fwd_rec+1,i1,LBj+1,
     *                    n1,n2,ubarG(i1:UBi,LBj:UBj,ifwd2))
        CALL get_3dfld_tile(fwdfile,'u',fwd_rec+1,
     *                    i1,LBj+1,1,n1,n2,N,     ! start (nc file), count
     *                    uG(i1:UBi,LBj:UBj,1:N,ifwd2))  

! in v-locations:
        i2=max(1,LBj)
        n1=UBi-LBi+1
        n2=UBj-i2+1
        CALL get_2dfld_tile(fwdfile,'vbar',fwd_rec+1,LBi+1,i2,
     *                    n1,n2,vbarG(LBi:UBi,i2:UBj,ifwd2))
        CALL get_3dfld_tile(fwdfile,'v',fwd_rec+1,
     *                    LBi+1,i2,1,n1,n2,N,     ! start (nc file), count
     *                    vG(LBi:UBi,i2:UBj,1:N,ifwd2))

#else /* serial */

        CALL get_2dfld(fwdfile,'zeta',fwd_rec+1,Lp,Mp,
     *                                          zetaG(0:L,0:M,ifwd2))
        CALL get_2dfld(fwdfile,'ubar',fwd_rec+1,L ,Mp,
     *                                          ubarG(1:L,0:M,ifwd2))
        CALL get_2dfld(fwdfile,'vbar',fwd_rec+1,Lp,M ,
     *                                          vbarG(0:L,1:M,ifwd2))
        CALL get_3dfld(fwdfile,'u'   ,fwd_rec+1,L ,Mp,N,
     *                                          uG(1:L,0:M,1:N,ifwd2))
        CALL get_3dfld(fwdfile,'v'   ,fwd_rec+1,Lp,M ,N,
     *                                          vG(0:L,1:M,1:N,ifwd2))
        CALL get_3dfld(fwdfile,'temp',fwd_rec+1,Lp,Mp,N,
     *                                          tG(0:L,0:M,1:N,ifwd2,1))
        CALL get_3dfld(fwdfile,'salt',fwd_rec+1,Lp,Mp,N,
     *                                          tG(0:L,0:M,1:N,ifwd2,2))

#ifndef ANA_FWD_MIX
        CALL get_3dfld(fwdfile,'AKv'  ,fwd_rec+1,Lp,Mp,Np,
     *                                          AkvG(0:L,0:M,0:N,ifwd2))
        CALL get_3dfld(fwdfile,'AKt'  ,fwd_rec+1,Lp,Mp,Np,
     *                                          AktG(0:L,0:M,0:N,ifwd2))
#endif

#endif /* MPI or serial */

#ifdef NS_PERIODIC

#ifdef MPI
        STOP('periodic BC not coded in MPI')
#else /* serial */
        CALL exchange_2d(zetaG(:,:,ifwd2),LBi,UBi,LBj,UBj,Lm,Mm)
        CALL exchange_2d(ubarG(:,:,ifwd2),LBi,UBi,LBj,UBj,Lm,Mm)
        CALL exchange_2d(vbarG(:,:,ifwd2),LBi,UBi,LBj,UBj,Lm,Mm)
        CALL exchange_3d(uG(:,:,:,ifwd2),LBi,UBi,LBj,UBj,1,N,Lm,Mm)
        CALL exchange_3d(vG(:,:,:,ifwd2),LBi,UBi,LBj,UBj,1,N,Lm,Mm)
        CALL exchange_3d(tG(:,:,:,ifwd2,1),LBi,UBi,LBj,UBj,1,N,Lm,Mm)
        CALL exchange_3d(tG(:,:,:,ifwd2,2),LBi,UBi,LBj,UBj,1,N,Lm,Mm)

# ifndef ANA_FWD_MIX
       CALL exchange_3d(AkvG(:,:,:,ifwd2),LBi,UBi,LBj,UBj,0,N,Lm,Mm)
       CALL exchange_3d(AktG(:,:,:,ifwd2),LBi,UBi,LBj,UBj,0,N,Lm,Mm)
# endif
#endif /* NS_PERIODIC */

#endif /* MPI or serial */

       ELSE 
        ! end of record, fwd_rec==nfwd, allow for round-up error
        IF (t0.gt.fwd_time(nfwd)+0.1d0*dt) THEN
         write (*,*) 'STOP in get_state, since t0>fwd_time(last)'
         STOP
        END IF                
       END IF  ! (fwd_rec.lt.nfwd)
      END IF  ! (t0.gt.fwd_timeG(ifwd2))

      RETURN
      END

#endif /* defined TANGENT */
