#include "cppdefs.h"
      MODULE mod_parallel

      integer :: comm2d,myrank,MaxTile,comm_size,irank
      integer :: rankleft,rankright,rankbottom,ranktop,stride
      integer :: stride_r3d_bt,stride_u3d_bt,stride_v3d_bt,
     *           stride_3d_lr
      integer :: dims_mpi(2)
      integer :: coords(2)
      logical :: periods_mpi(2)
      logical :: reorder

#ifdef ADJOINT
      integer, parameter :: n_ext_adj=3
      integer :: LBi0,UBi0,LBj0,UBj0
#endif

      END MODULE mod_parallel
