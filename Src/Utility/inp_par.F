#include "cppdefs.h"

      SUBROUTINE inp_par
      USE mod_ocean
      
#ifndef NONLIN_EOS
      USE mod_eoscoef
#endif

      implicit none
#ifdef MPI
      include 'mpif.h'
#endif
      integer :: myid,ierr


! read from file: on the master processor

      myid=0
#ifdef MPI
      call mpi_comm_rank(mpi_comm_world,myid,ierr)
#endif
      IF (myid.eq.0) THEN

       READ (stdinp,'(a)')
       READ (stdinp,*) ntimes

       READ (stdinp,'(a)')
       READ (stdinp,*) dt

       READ (stdinp,'(a)') 
       READ (stdinp,*) ndtfast

       READ (stdinp,'(a)')
       READ (stdinp,*) nhis 

       READ (stdinp,'(a)')
       READ (stdinp,*) Akv0
       WRITE (*,*) 'Background Eddy Visc. Akv0 = ',Akv0
 
       READ (stdinp,'(a)')
       READ (stdinp,*) Akt0(1:2) 
 
       READ (stdinp,'(a)')
       READ (stdinp,*) rdrg
       WRITE (*,*) 'Friction Coefficient rdrg = ',rdrg
 
       READ (stdinp,'(a)')
       READ (stdinp,*) vis2
       WRITE (*,*) 'Horizontal viscosity = ',vis2
 
       READ (stdinp,'(a)')
       READ (stdinp,*) dif2
       WRITE (*,*) 'Horizontal diffusion = ',dif2
 
       READ (stdinp,'(a)')
#ifdef NONLIN_EOS
       READ (stdinp,'(a)')
#else 
       READ (stdinp,*) R0,Tcoef,Scoef,T0,S0
       WRITE (*,*) 'Linear EOS parameters: '
       WRITE (*,*) 'R0   =',R0
       WRITE (*,*) 'Tcoef=',Tcoef
       WRITE (*,*) 'Scoef=',Scoef
       WRITE (*,*) 'T0=',T0
       WRITE (*,*) 'S0=',S0
#endif

       READ (stdinp,'(a)')
       READ (stdinp,'(a80)') grdfile
 
       READ (stdinp,'(a)')
       READ (stdinp,'(a80)') inifile
 
       READ (stdinp,'(a)')
       READ (stdinp,'(a80)') frcfile
 
       READ (stdinp,'(a)')
       READ (stdinp,'(a80)') fwdfile
 
       READ (stdinp,'(a)')
       READ (stdinp,'(a80)') hisfile

#ifdef MPI
       READ (stdinp,'(a)')
       READ (stdinp,*) ntiles_xi,ntiles_eta
#endif /* MPI */


      ENDIF ! myid.eq.0
      call mpi_barrier(mpi_comm_world,ierr)

#ifdef MPI
      call mpi_bcast(ntimes,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(dt,1,mpi_double_precision,0,mpi_comm_world,ierr)
      call mpi_bcast(ndtfast,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(nhis,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(Akv0,1,mpi_double_precision,0,mpi_comm_world,ierr)
      call mpi_bcast(Akt0(1),2,mpi_double_precision,0,mpi_comm_world,
     *               ierr)
      call mpi_bcast(rdrg, 1,mpi_double_precision,0,mpi_comm_world,ierr)
      call mpi_bcast(vis2, 1,mpi_double_precision,0,mpi_comm_world,ierr)
      call mpi_bcast(dif2, 1,mpi_double_precision,0,mpi_comm_world,ierr)
      call mpi_bcast(R0,   1,mpi_double_precision,0,mpi_comm_world,ierr)
      call mpi_bcast(Tcoef,1,mpi_double_precision,0,mpi_comm_world,ierr)
      call mpi_bcast(Scoef,1,mpi_double_precision,0,mpi_comm_world,ierr)
      call mpi_bcast(T0,1,mpi_double_precision,0,mpi_comm_world,ierr)
      call mpi_bcast(S0,1,mpi_double_precision,0,mpi_comm_world,ierr)
      call mpi_bcast(grdfile,80,mpi_character,0,mpi_comm_world,ierr)
      call mpi_bcast(inifile,80,mpi_character,0,mpi_comm_world,ierr)
      call mpi_bcast(frcfile,80,mpi_character,0,mpi_comm_world,ierr)
      call mpi_bcast(fwdfile,80,mpi_character,0,mpi_comm_world,ierr)
      call mpi_bcast(hisfile,80,mpi_character,0,mpi_comm_world,ierr)

      call mpi_bcast(ntiles_xi,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(ntiles_eta,1,mpi_integer,0,mpi_comm_world,ierr)
#endif /* MPI */

      dtfast=dt/REAL(ndtfast,8)

      RETURN
      END
