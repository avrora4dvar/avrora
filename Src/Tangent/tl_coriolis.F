#include "cppdefs.h"

      SUBROUTINE tl_coriolis
#if defined TANGENT && (defined UV_COR || (defined CURVGRID && defined UV_ADV))
      USE mod_ocean
      implicit none
      integer i,j,k
      real(8) :: cff, tl_cff,tmp5,tl_tmp5

!-----------------------------------------------------------------------
!  Add in Coriolis
!-----------------------------------------------------------------------
      DO k=1,N

       ! compute UFx(i,j) = Hz*v*f*Area at rho(ij)-location
       ! compute VFe(i,j) = Hz*u*f*Area at rho(ij)-location


!=== TL double cycle 1:
#ifdef MPI
       DO j=max(1,Jstr-1),Jend
        DO i=max(1,Istr-1),Iend

#else /* SERIAL */

# if defined NS_PERIODIC
       DO j=Jstr-1,Jend
# else
       DO j=Jstr,Jend
# endif
        DO i=Istr,Iend
#endif /* mpi or serial */

         IF (mask_r(i,j).ne.0) THEN

            tmp5=
#  ifdef UV_COR
     *          fdxdy(i,j)     
#  endif
#  if (defined CURVGRID && defined UV_ADV)
     *         +0.5d0*((v(i,j  ,k,nrhs)+v(i,j+1,k,nrhs))*dndx(i,j)-   
     *                 (u(i  ,j,k,nrhs)+u(i+1,j,k,nrhs))*dmde(i,j))
#  endif

#  if (defined CURVGRID && defined UV_ADV)
            tl_tmp5=0.5d0*
     *               ((tl_v(i,j  ,k,nrhs)+tl_v(i,j+1,k,nrhs))*dndx(i,j)-   
     *                (tl_u(i  ,j,k,nrhs)+tl_u(i+1,j,k,nrhs))*dmde(i,j))
#  endif

            cff=0.5d0*Hz(i,j,k)*tmp5

            tl_cff=0.5d0*(tl_Hz(i,j,k)*tmp5
#  if (defined CURVGRID && defined UV_ADV)
     *                    +Hz(i,j,k)*tl_tmp5  
#  endif
     *                   )

! NL: UFx(i,j)=cff*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))
          tl_UFx(i,j)=cff*(tl_v(i,j,k,nrhs)+tl_v(i,j+1,k,nrhs))+
     *                tl_cff*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))
! NL:          VFe(i,j)=cff*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))
          tl_VFe(i,j)=cff*(tl_u(i,j,k,nrhs)+tl_u(i+1,j,k,nrhs))+
     *                tl_cff*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs)) 
         END IF


        END DO
       END DO

!=== TL double cycle 2:
       DO j=JstrU,JendU
        DO i=IstrU,IendU
         tl_ru(i,j,k,nrhs)=tl_ru(i,j,k,nrhs)
     *                         +0.5d0*(tl_UFx(i,j)+tl_UFx(i-1,j))
        END DO
       END DO

!=== TL double cycle 3:
       DO j=JstrV,JendV
        DO i=IstrV,IendV
         tl_rv(i,j,k,nrhs)=tl_rv(i,j,k,nrhs)
     *                         -0.5d0*(tl_VFe(i,j)+tl_VFe(i,j-1))
        END DO
       END DO

      END DO ! END "DO k=1,N, Coriolis computation"

#endif /* TANGENT && UV_COR */

      RETURN
      END
