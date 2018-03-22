#include "cppdefs.h"

#ifdef MPI

      SUBROUTINE read_grid
      USE mod_ocean
      USE mod_parallel

      implicit none
      include 'mpif.h'
      integer :: recnum,i,j,n1,n2,i1,i2

      recnum=1

! in rho-locations (2D):
      n1=UBi-LBi+1
      n2=UBj-LBj+1

! get_2dflrd_tile arguments: grdfile,var,recnum,start index 1,
!                            start index 2, count 1, count 2, buffer
      CALL get_2dfld_tile(grdfile,'h',recnum,LBi+1,LBj+1,
     *                    n1,n2,h(LBi:UBi,LBj:UBj))
      CALL get_2dfld_tile(grdfile,'f',recnum,LBi+1,LBj+1,
     *                    n1,n2,f(LBi:UBi,LBj:UBj))
      CALL get_2dfld_tile(grdfile,'dx',recnum,LBi+1,LBj+1,
     *                    n1,n2,dx(LBi:UBi,LBj:UBj))
      CALL get_2dfld_tile(grdfile,'dy',recnum,LBi+1,LBj+1,
     *                    n1,n2,dy(LBi:UBi,LBj:UBj))

#if defined CURVGRID && defined UV_ADV
      CALL get_2dfld_tile(grdfile,'dmde',recnum,LBi+1,LBj+1,
     *                    n1,n2,dmde(LBi:UBi,LBj:UBj))
      CALL get_2dfld_tile(grdfile,'dndx',recnum,LBi+1,LBj+1,
     *                    n1,n2,dndx(LBi:UBi,LBj:UBj))
#endif

      CALL get_2dfld_int_tile(grdfile,'mask_rho',recnum,LBi+1,LBj+1,
     *                    n1,n2,mask_r(LBi:UBi,LBj:UBj))

! in u-locations (2D):  
      i1=max(1,LBi)  ! index of 1st element on roms grid
      n1=UBi-i1+1
      n2=UBj-LBj+1

! note: first element is i1, not i1+1, since grid is staggered 
! (1st element in the grid file corresponds to element 1 on roms grid)
      CALL get_2dfld_int_tile(grdfile,'mask_u',recnum,
     *                    i1,LBj+1,
     *                    n1,n2,mask_u(i1:UBi,LBj:UBj))

! in v-locations (2D):
      i2=max(1,LBj)
      n1=UBi-LBi+1
      n2=UBj-i2+1
      CALL get_2dfld_int_tile(grdfile,'mask_v',recnum,
     *                    LBi+1,i2,
     *                    n1,n2,mask_v(LBi:UBi,i2:UBj))

! in rho-locations (3D):
      n1=UBi-LBi+1
      n2=UBj-LBj+1
      CALL get_3dfld_tile(grdfile,'z0_r',recnum,
     *                    LBi+1,LBj+1,1,n1,n2,N,     ! start (nc file), count
     *                    z0_r(LBi:UBi,LBj:UBj,1:N))  
! in w-locations (3D):
      CALL get_3dfld_tile(grdfile,'z0_w',recnum,
     *                    LBi+1,LBj+1,1,n1,n2,Np,
     *                    z0_w(LBi:UBi,LBj:UBj,0:N))  


#if defined NS_PERIODIC

      STOP('STOP since NO MPI WRITTEN FOR PERIDIC CHANNEL')

      CALL exchange_2d(h,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(f,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(dx,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(dy,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d_int(mask_r,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d_int(mask_u,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d_int(mask_v,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_3d(z0_r,LBi,UBi,LBj,UBj,1,N,Lm,Mm)
      CALL exchange_3d(z0_w,LBi,UBi,LBj,UBj,1,N,Lm,Mm)

#endif /* NS_PERIODIC */

!--------------
! Metrics: 
!--------------
! (dx and dy are sizes of the rho-cell, measured via middle)
! (dx_v and dy_u are lengths of the edges of the grid)
      
      
      DO j=LBj,UBj   ! note, in serial we had : Jstr,Jend
      DO i=LBi,UBi   ! note, in serial we had : Istr,Iend
       Area(i,j)=dx(i,j)*dy(i,j)
       oArea(i,j)=DBLE(mask_r(i,j))/Area(i,j)
       oh(i,j)=1.d0/h(i,j)

#if defined UV_COR || defined UV_COR_MASK
       fdxdy(i,j)=f(i,j)*Area(i,j)
#endif

#if defined UV_COR
!
       IF (mask_r(i,j).ne.0) THEN
         omask_vv(i,j)=0.5d0
         omask_uu(i,j)=0.5d0
       END IF

#elif defined UV_COR_MASK

       IF (mask_r(i,j).ne.0) THEN

# ifdef NS_PERIODIC
        omask_vv(i,j)=1.d0/DBLE(mask_v(i,j)+mask_v(i,j+1))
# else
        IF (j.eq.1 .or. j.eq.M-1) THEN
         omask_vv(i,j)=1.d0
        ELSE
         omask_vv(i,j)=1.d0/DBLE(mask_v(i,j)+mask_v(i,j+1))
        END IF
# endif 

        IF (i.eq.1 .or. i.eq.L-1) THEN
         omask_uu(i,j)=1.d0
        ELSE
         omask_uu(i,j)=1.d0/DBLE(mask_u(i,j)+mask_u(i+1,j)) 
        END IF
       END IF

#endif /* UV_COR, elif UV_COR_MASK */
        
      END DO
      END DO

! in u-locations each tile (excl. left (ghost) edge of each tile) 
      DO j=LBj,UBj     ! note, in serial, it was: JstrU-1,JendU+1
      DO i=LBi+1,UBi   ! notem in serial, it was: IstrU,IendU+1
       dx_u(i,j)=0.5d0*(dx(i-1,j)+dx(i,j))
       dy_u(i,j)=0.5d0*(dy(i-1,j)+dy(i,j))
       Area_u(i,j)=dx_u(i,j)*dy_u(i,j)
       oArea_u(i,j)=DBLE(mask_u(i,j))/Area_u(i,j)
      END DO
      END DO

! in v-locations each tile (excl. left (ghost) edge of each tile) 
      DO j=LBj+1,UBj   ! note in serial it was: JstrV-1,JendV+1
      DO i=LBi,UBi     ! note in serial it was: IstrV-1,IendV+1
       dx_v(i,j)=0.5d0*(dx(i,j-1)+dx(i,j))
       dy_v(i,j)=0.5d0*(dy(i,j-1)+dy(i,j))
       Area_v(i,j)=dx_v(i,j)*dy_v(i,j)
       oArea_v(i,j)=DBLE(mask_v(i,j))/Area_v(i,j)
      END DO
      END DO
      
      CALL mpi_exchange_2d_2(dx_u(LBi:UBi,LBj:UBj),
     *                     LBi,UBi,LBj,UBj,
     *                     IstrU_E,IendU_E,JstrU_E,JendU_E,
     *                     comm2d, stride, 
     *                     rankleft,rankright,ranktop,rankbottom) 
      CALL mpi_exchange_2d_2(dy_u(LBi:UBi,LBj:UBj),
     *                     LBi,UBi,LBj,UBj,
     *                     IstrU_E,IendU_E,JstrU_E,JendU_E,
     *                     comm2d, stride, 
     *                     rankleft,rankright,ranktop,rankbottom) 

      CALL mpi_exchange_2d_2(dx_v(LBi:UBi,LBj:UBj),
     *                     LBi,UBi,LBj,UBj,
     *                     IstrV_E,IendV_E,JstrV_E,JendV_E,
     *                     comm2d, stride, 
     *                     rankleft,rankright,ranktop,rankbottom) 
      CALL mpi_exchange_2d_2(dy_v(LBi:UBi,LBj:UBj),
     *                     LBi,UBi,LBj,UBj,
     *                     IstrV_E,IendV_E,JstrV_E,JendV_E,
     *                     comm2d, stride, 
     *                     rankleft,rankright,ranktop,rankbottom) 

#if defined NS_PERIODIC
      CALL exchange_2d(dx_u,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(dy_u,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(dx_v,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(dy_v,LBi,UBi,LBj,UBj,Lm,Mm)

      CALL exchange_2d(Area,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(oArea,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(oh,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(Area_u,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(oArea_u,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(Area_v,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(oArea_v,LBi,UBi,LBj,UBj,Lm,Mm)

# if defined UV_COR || defined UV_COR_MASK
      CALL exchange_2d(fdxdy,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(omask_vv,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(omask_uu,LBi,UBi,LBj,UBj,Lm,Mm)
# endif /* defined UV_COR || defined UV_COR_MASK */

#endif /* NS_PERIODIC */

! These were first needed when horizontal dissipation was included:
      DO j=LBj,UBj
      DO i=LBi,UBi

! dyodx=pmon
       dyodx(i,j)=dy(i,j)/dx(i,j)
       dxody(i,j)=dx(i,j)/dy(i,j)
       ody(i,j)=1.d0/dy(i,j)
       odx(i,j)=1.d0/dx(i,j)

      END DO
      END DO

! At Psi locations (w/ 1 ghost layer, 
! this should include the edge of the tile, in other words, the edge
! of interior rho-cells)
      DO j=LBj+1,UBj
      DO i=LBi+1,UBi
       dx_p(i,j)=0.25d0*(dx(i-1,j-1)+dx(i-1,j)+dx(i,j-1)+dx(i,j))
       dy_p(i,j)=0.25d0*(dy(i-1,j-1)+dy(i-1,j)+dy(i,j-1)+dy(i,j))
       dyodx_p(i,j)=dy_p(i,j)/dx_p(i,j)
       dxody_p(i,j)=1.d0/dyodx_p(i,j)

       IF (mask_r(i-1,j-1)*mask_r(i-1,j)*mask_r(i,j-1)*mask_r(i,j)
     *     .eq.0) THEN
        pmask(i,j)=0.d0
       ELSE 
        pmask(i,j)=1.d0
       END IF

      END DO
      END DO

#if defined TS_DIF2
! These are needed for tracer horizontal diffusion:
      DO j=LBj,UBj
      DO i=LBi,UBi

       IF (i.gt.LBi) THEN
        dyodx_u(i,j)=(dy(i-1,j)+dy(i,j))/
     *               (dx(i-1,j)+dx(i,j)) 
       END IF

       IF (j.gt.LBj) THEN
        dxody_v(i,j)=(dx(i,j-1)+dx(i,j))/
     *               (dy(i,j-1)+dy(i,j)) 
       END IF
 
      END DO
      END DO
#endif


#if defined NS_PERIODIC
      CALL exchange_2d(dx_p,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(dy_p,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(dyodx_p,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(dxody_p,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(dyodx_u,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(dxody_v,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(pmask,LBi,UBi,LBj,UBj,Lm,Mm)
#endif /* NS_PERIODIC */

      RETURN
      END

#else /* serial */
!!!!!!!! v2.1 SERIAL CODE 

      SUBROUTINE read_grid
      USE mod_ocean

      implicit none
      integer :: recnum,i,j

      recnum=1

      CALL get_2dfld(grdfile,'h',recnum,Lp,Mp,h(0:L,0:M))

      CALL get_2dfld(grdfile,'f',recnum,Lp,Mp,f(0:L,0:M))
      CALL get_2dfld(grdfile,'dx',recnum,Lp,Mp,dx(0:L,0:M))
      CALL get_2dfld(grdfile,'dy',recnum,Lp,Mp,dy(0:L,0:M))
#if defined CURVGRID && defined UV_ADV 
      CALL get_2dfld(grdfile,'dmde',recnum,Lp,Mp,dmde(0:L,0:M))
      CALL get_2dfld(grdfile,'dndx',recnum,Lp,Mp,dndx(0:L,0:M))
#endif

      CALL get_2dfld_int(grdfile,'mask_rho',recnum,Lp,Mp,
     *                                                mask_r(0:L,0:M))
      CALL get_2dfld_int(grdfile,'mask_u',recnum,L,Mp,mask_u(1:L,0:M))
      CALL get_2dfld_int(grdfile,'mask_v',recnum,Lp,M,mask_v(0:L,1:M))

      CALL get_3dfld(grdfile,'z0_r',recnum,Lp,Mp,N,z0_r(0:L,0:M,1:N))  
      CALL get_3dfld(grdfile,'z0_w',recnum,Lp,Mp,Np,z0_w(0:L,0:M,0:N))

#if defined NS_PERIODIC
      CALL exchange_2d(h,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(f,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(dx,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(dy,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d_int(mask_r,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d_int(mask_u,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d_int(mask_v,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_3d(z0_r,LBi,UBi,LBj,UBj,1,N,Lm,Mm)
      CALL exchange_3d(z0_w,LBi,UBi,LBj,UBj,1,N,Lm,Mm)

#endif /* NS_PERIODIC */

!--------------
! Metrics: 
!--------------
! (dx and dy are sizes of the rho-cell, measured via middle)
! (dx_v and dy_u are lengths of the edges of the grid)

      DO j=Jstr,Jend
      DO i=Istr,Iend
       Area(i,j)=dx(i,j)*dy(i,j)
       oArea(i,j)=DBLE(mask_r(i,j))/Area(i,j)
       oh(i,j)=1.d0/h(i,j)

#if defined UV_COR || defined UV_COR_MASK
       fdxdy(i,j)=f(i,j)*Area(i,j)
#endif

#if defined UV_COR
!
       IF (mask_r(i,j).ne.0) THEN
         omask_vv(i,j)=0.5d0
         omask_uu(i,j)=0.5d0
       END IF

#elif defined UV_COR_MASK

       IF (mask_r(i,j).ne.0) THEN

# ifdef NS_PERIODIC
        omask_vv(i,j)=1.d0/DBLE(mask_v(i,j)+mask_v(i,j+1))
# else
        IF (j.eq.Jstr .or. j.eq.Jend) THEN
         omask_vv(i,j)=1.d0
        ELSE
         omask_vv(i,j)=1.d0/DBLE(mask_v(i,j)+mask_v(i,j+1))
        END IF
# endif 

        IF (i.eq.Istr .or. i.eq.Iend) THEN
         omask_uu(i,j)=1.d0
        ELSE
         omask_uu(i,j)=1.d0/DBLE(mask_u(i,j)+mask_u(i+1,j)) 
        END IF
       END IF

#endif /* UV_COR, elif UV_COR_MASK */
        
      END DO
      END DO

      DO j=JstrU-1,JendU+1
      DO i=IstrU-1,IendU+1
       dx_u(i,j)=0.5d0*(dx(i-1,j)+dx(i,j))
       dy_u(i,j)=0.5d0*(dy(i-1,j)+dy(i,j))
       Area_u(i,j)=dx_u(i,j)*dy_u(i,j)
       oArea_u(i,j)=DBLE(mask_u(i,j))/Area_u(i,j)
      END DO
      END DO

      DO j=JstrV-1,JendV+1
      DO i=IstrV-1,IendV+1
       dx_v(i,j)=0.5d0*(dx(i,j-1)+dx(i,j))
       dy_v(i,j)=0.5d0*(dy(i,j-1)+dy(i,j))
       Area_v(i,j)=dx_v(i,j)*dy_v(i,j)
       oArea_v(i,j)=DBLE(mask_v(i,j))/Area_v(i,j)
      END DO
      END DO


#if defined NS_PERIODIC
      CALL exchange_2d(dx_u,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(dy_u,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(dx_v,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(dy_v,LBi,UBi,LBj,UBj,Lm,Mm)

      CALL exchange_2d(Area,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(oArea,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(oh,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(Area_u,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(oArea_u,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(Area_v,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(oArea_v,LBi,UBi,LBj,UBj,Lm,Mm)

# if defined UV_COR || defined UV_COR_MASK
      CALL exchange_2d(fdxdy,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(omask_vv,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(omask_uu,LBi,UBi,LBj,UBj,Lm,Mm)
# endif /* defined UV_COR || defined UV_COR_MASK */

#endif /* NS_PERIODIC */

! These were first needed when horizontal dissipation was included:
      DO j=LBj,UBj
      DO i=LBi,UBi

! dyodx=pmon
       dyodx(i,j)=dy(i,j)/dx(i,j)
       dxody(i,j)=dx(i,j)/dy(i,j)
       ody(i,j)=1.d0/dy(i,j)
       odx(i,j)=1.d0/dx(i,j)


! At Psi locations:
       IF (i.gt.LBi .and. j.gt.LBj) THEN
        dx_p(i,j)=0.25d0*(dx(i-1,j-1)+dx(i-1,j)+dx(i,j-1)+dx(i,j))
        dy_p(i,j)=0.25d0*(dy(i-1,j-1)+dy(i-1,j)+dy(i,j-1)+dy(i,j))
        dyodx_p(i,j)=dy_p(i,j)/dx_p(i,j)
        dxody_p(i,j)=1.d0/dyodx_p(i,j)

        IF (mask_r(i-1,j-1)*mask_r(i-1,j)*mask_r(i,j-1)*mask_r(i,j)
     *      .eq.0) THEN
         pmask(i,j)=0.d0
        ELSE 
         pmask(i,j)=1.d0
        END IF
 
       END IF

      END DO
      END DO

#if defined TS_DIF2
! These are needed for tracer horizontal diffusion:
      DO j=LBj,UBj
      DO i=LBi,UBi

       IF (i.gt.LBi) THEN
        dyodx_u(i,j)=(dy(i-1,j)+dy(i,j))/
     *               (dx(i-1,j)+dx(i,j)) 
       END IF

       IF (j.gt.LBj) THEN
        dxody_v(i,j)=(dx(i,j-1)+dx(i,j))/
     *               (dy(i,j-1)+dy(i,j)) 
       END IF
 
      END DO
      END DO
#endif


#if defined NS_PERIODIC
      CALL exchange_2d(dx_p,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(dy_p,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(dyodx_p,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(dxody_p,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(dyodx_u,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(dxody_v,LBi,UBi,LBj,UBj,Lm,Mm)
      CALL exchange_2d(pmask,LBi,UBi,LBj,UBj,Lm,Mm)
#endif /* NS_PERIODIC */

      RETURN
      END
#endif /* MPI or serial */
