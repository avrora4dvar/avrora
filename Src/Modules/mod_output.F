#include "cppdefs.h"
      MODULE mod_output


#ifdef MPI
! do not use structures

! dimid:
      integer :: dimid_xi_rho
      integer :: dimid_eta_rho
      integer :: dimid_xi_u
      integer :: dimid_eta_u
      integer :: dimid_xi_v
      integer :: dimid_eta_v
      integer :: dimid_N
      integer :: dimid_Np
      integer :: dimid_time

!varid
      integer :: varid_ocean_time
      integer :: varid_zeta
      integer :: varid_ubar
      integer :: varid_vbar
      integer :: varid_u
      integer :: varid_v  
      integer :: varid_temp
      integer :: varid_salt

#else /* serial */
! use structure
      TYPE state_dimid
       integer :: xi_rho
       integer :: eta_rho
       integer :: xi_u
       integer :: eta_u
       integer :: xi_v
       integer :: eta_v
       integer :: N
       integer :: Np
       integer :: time
      END TYPE state_dimid

      TYPE state_varid
       integer :: ocean_time
       integer :: zeta
       integer :: ubar
       integer :: vbar
       integer :: u
       integer :: v  
       integer :: temp
       integer :: salt
      END TYPE state_varid

      TYPE (state_dimid) :: his_dimid
      TYPE (state_varid) :: his_varid


#endif /* MPI or serial */

      END MODULE mod_output
