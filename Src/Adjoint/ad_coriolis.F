#include "cppdefs.h"

      SUBROUTINE ad_coriolis
#if defined ADJOINT && (defined UV_COR || (defined CURVGRID && defined UV_ADV))
      USE mod_ocean
#ifdef MPI
      USE mod_parallel
#endif
      implicit none
      integer i,j,k
      real(8) :: cff, ad_cff,tmp,tmp1,tmp2,tmp3,tmp4,tmp5,ad_tmp5

      ad_cff=0.d0
      ad_tmp5=0.d0

!#ifdef MPI
!      if (myrank.eq.0) then      
!#endif
!       write (*,*) 'ad_ru(LBi,LBj,1,nrhs)=',ad_ru(LBi,LBj,1,nrhs)
!#ifdef MPI
!      endif
!#endif

!-----------------------------------------------------------------------
!  ADJ TO: add in Coriolis
!-----------------------------------------------------------------------
      DO k=1,N

!=== ADJ to TL double cycle 3:
       DO j=JstrV,JendV
        DO i=IstrV,IendV
! TL (recurs.):         tl_rv(i,j,k,nrhs)=tl_rv(i,j,k,nrhs)
! TL:     *                         -0.5d0*(tl_VFe(i,j)+tl_VFe(i,j-1))

         ad_VFe(i,j)=ad_VFe(i,j)-0.5d0*ad_rv(i,j,k,nrhs)
         ad_VFe(i,j-1)=ad_VFe(i,j-1)-0.5d0*ad_rv(i,j,k,nrhs)

        END DO
       END DO


!=== ADJ to TL bouble cycle 2:
       DO j=JstrU,JendU
        DO i=IstrU,IendU

! TL (recurs.): tl_ru(i,j,k,nrhs)=tl_ru(i,j,k,nrhs)
! TL:     *                  +0.5d0*(tl_UFx(i,j)+tl_UFx(i-1,j))

         ad_UFx(i,j)=ad_UFx(i,j)+0.5d0*ad_ru(i,j,k,nrhs)
         ad_UFx(i-1,j)=ad_UFx(i-1,j)+0.5d0*ad_ru(i,j,k,nrhs)

        END DO
       END DO

!=== ADJ to TL bouble cycle 1:
#if defined NS_PERIODIC
       DO j=Jstr-1,Jend
#else
       DO j=Jstr,Jend
#endif
        DO i=Istr,Iend

         IF (mask_r(i,j).ne.0) THEN


            tmp5=
#  ifdef UV_COR
     *          fdxdy(i,j)     
#  endif
#  if (defined CURVGRID && defined UV_ADV)
     *         +0.5d0*((v(i,j  ,k,nrhs)+v(i,j+1,k,nrhs))*dndx(i,j)-   
     *                 (u(i  ,j,k,nrhs)+u(i+1,j,k,nrhs))*dmde(i,j))
#  endif

            cff=0.5d0*Hz(i,j,k)*tmp5


! NL:          VFe(i,j)=cff*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))

! TL:          tl_VFe(i,j)= cff*(tl_u(i,j,k,nrhs)+tl_u(i+1,j,k,nrhs))+
! TL:     *                  tl_cff*(u(i,j,k,nrhs)+u(i+1,j,k,nrhs)) 

          tmp1=ad_VFe(i,j)
          ad_u(i,j,k,nrhs)=ad_u(i,j,k,nrhs)+cff*tmp1
          ad_u(i+1,j,k,nrhs)=ad_u(i+1,j,k,nrhs)+cff*tmp1
          ad_cff=ad_cff+(u(i,j,k,nrhs)+u(i+1,j,k,nrhs))*tmp1
          ad_VFe(i,j)=0.d0

!==
! NL: UFx(i,j)=cff*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))
!==
! TL:          tl_UFx(i,j)=cff*(tl_v(i,j,k,nrhs)+tl_v(i,j+1,k,nrhs))+
! TL:     *                  tl_cff*(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))
!==
          tmp1=ad_UFx(i,j)
          ad_v(i,j,k,nrhs)=ad_v(i,j,k,nrhs)+cff*tmp1
          ad_v(i,j+1,k,nrhs)=ad_v(i,j+1,k,nrhs)+cff*tmp1
          ad_cff=ad_cff+(v(i,j,k,nrhs)+v(i,j+1,k,nrhs))*tmp1
          ad_UFx(i,j)=0.d0

!-- (NL cff is defined above, background cff is needed)
! TL:            tl_cff=0.5d0*(tl_Hz(i,j,k)*tmp5
!#  if (defined CURVGRID && defined UV_ADV)
!    *                         +Hz(i,j,k)*tl_tmp5  
!#  endif
!     *                   )

          ad_Hz(i,j,k)=ad_Hz(i,j,k)+0.5d0*tmp5*ad_cff
#  if (defined CURVGRID && defined UV_ADV)
          ad_tmp5=ad_tmp5+Hz(i,j,k)*0.5d0*ad_cff         
#  endif
          ad_cff=0.d0

#  if (defined CURVGRID && defined UV_ADV)
! An extra line, tl_tmp5:
! TL:           tl_tmp5=0.5d0*
!     *                ((tl_v(i,j  ,k,nrhs)+tl_v(i,j+1,k,nrhs))*dndx(i,j)-   
!     *                 (tl_u(i  ,j,k,nrhs)+tl_u(i+1,j,k,nrhs))*dmde(i,j))

          tmp=0.5d0*ad_tmp5
          tmp3=tmp*dndx(i,j) 
          tmp4=tmp*dmde(i,j) 

          ad_v(i,j  ,k,nrhs)=ad_v(i,j  ,k,nrhs)+tmp3
          ad_v(i,j+1,k,nrhs)=ad_v(i,j+1,k,nrhs)+tmp3
          ad_u(i  ,j,k,nrhs)=ad_u(i  ,j,k,nrhs)-tmp4
          ad_u(i+1,j,k,nrhs)=ad_u(i+1,j,k,nrhs)-tmp4

          ad_tmp5=0.d0
#  endif

         END IF

        END DO
       END DO


      END DO ! END "DO k=1,N, Coriolis computation"

#endif /* ADJOINT && UV_COR || UV_COR_MASK */

      RETURN
      END
