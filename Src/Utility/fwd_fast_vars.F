#include "cppdefs.h"

#if defined TANGENT || defined ADJOINT
      SUBROUTINE fwd_fast_vars

      USE mod_ocean
 
      implicit none
      integer :: i,j

      DO j=LBj,UBj
       DO i=LBi,UBi
        zeta(i,j,1)=fwd_zeta(i,j,nnew)
        zeta(i,j,2)=fwd_zeta(i,j,nnew)
        zeta(i,j,3)=fwd_zeta(i,j,nnew)

        ubar(i,j,1)=fwd_ubar(i,j,nnew)
        ubar(i,j,2)=fwd_ubar(i,j,nnew)
        ubar(i,j,3)=fwd_ubar(i,j,nnew)

        vbar(i,j,1)=fwd_vbar(i,j,nnew)
        vbar(i,j,2)=fwd_vbar(i,j,nnew)
        vbar(i,j,3)=fwd_vbar(i,j,nnew)

       END DO
      END DO

      DO j=LBj,UBj
       DO i=LBi+1,UBi
        DU_avg1(i,j)=0.5d0*dy_u(i,j)*ubar(i,j,1)*
     *               (h(i,j)+h(i-1,j)+zeta(i,j,1)+zeta(i-1,j,1))
       END DO
      END DO

      DO j=LBj+1,UBj
       DO i=LBi,UBi
        DV_avg1(i,j)=0.5d0*dx_v(i,j)*vbar(i,j,1)*
     *               (h(i,j)+h(i,j-1)+zeta(i,j,1)+zeta(i,j-1,1))
       END DO
      END DO

      RETURN
      END

#endif /* TANGENT || ADJOINT */
