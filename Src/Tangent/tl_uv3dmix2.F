#include "cppdefs.h"

#ifdef TANGENT
      SUBROUTINE tl_uv3dmix2
      USE mod_ocean
#ifdef MPI
      USE mod_parallel    
#endif

      implicit none
      integer :: i,j,k
      real(8) :: tl_cff,tl_cff1,tl_cff2
      real(8) :: tmp,tmp1,tmp2,tmp3,tmp4

# ifdef UV_VIS2
      K_LOOP : DO k=1,N
!
!  Compute flux-components of the horizontal divergence of the stress
!  tensor (m5/s2) in XI- and ETA-directions.
!
        DO j=JstrV-1,JendV
!         DO i=min(1,IstrV-1),IendV ! new in v3.2, v3.1 had i=IstrV,IendV
          DO i=max(1,IstrV-1),IendV ! new in v3.2, v3.1 had i=IstrV,IendV
! NL:          cff=visc2_r(i,j)*Hz(i,j,k)*0.5d0*
! NL:    *           (dyodx(i,j)*
! NL:    *            ((ody(i  ,j)+ody(i+1,j))*u(i+1,j,k,nrhs)-       
! NL:    *             (ody(i-1,j)+ody(i  ,j))*u(i  ,j,k,nrhs))-        
! NL:    *            dxody(i,j)*                                            
! NL:    *            ((odx(i,j  )+odx(i,j+1))*v(i,j+1,k,nrhs)-      
! NL:    *             (odx(i,j-1)+odx(i,j  ))*v(i,j  ,k,nrhs)))

          tmp=visc2_r(i,j)*0.5d0*dyodx(i,j)
          tmp1=tmp*(ody(i  ,j)+ody(i+1,j))        
          tmp2=tmp*(ody(i-1,j)+ody(i  ,j))

          tmp=visc2_r(i,j)*0.5d0*dxody(i,j)
          tmp3=tmp*(odx(i,j  )+odx(i,j+1))
          tmp4=tmp*(odx(i,j-1)+odx(i,j  ))

! rewrite NL using defs above:
!         cff=Hz(i,j,k)*( tmp1*u(i+1,j,k,nrhs)-tmp2*u(i  ,j,k,nrhs)
!                        -tmp3*v(i,j+1,k,nrhs)+tmp4*v(i,j ,k,nrhs) )

          tl_cff= tl_Hz(i,j,k)*
     *             ( tmp1*u(i+1,j,k,nrhs)-tmp2*u(i  ,j,k,nrhs)
     *              -tmp3*v(i,j+1,k,nrhs)+tmp4*v(i,j+1,k,nrhs) )
     *           +Hz(i,j,k)*
     *             ( tmp1*tl_u(i+1,j,k,nrhs)-tmp2*tl_u(i  ,j,k,nrhs)
     *              -tmp3*tl_v(i,j+1,k,nrhs)+tmp4*tl_v(i,j ,k,nrhs) )
 
          tl_UFx(i,j)=dy(i,j)*dy(i,j)*tl_cff
          tl_VFe(i,j)=dx(i,j)*dx(i,j)*tl_cff
         END DO
        END DO

!#ifdef MPI
!        CALL mpi_exchange_2d(tl_UFx(LBi:UBi,LBj:UBj),
!     *                     LBi,UBi,LBj,UBj,
!     *                     Istr_E,Iend_E,Jstr_E,Jend_E,
!     *                     comm2d, stride, 
!     *                     rankleft,rankright,ranktop,rankbottom) 
!        CALL mpi_exchange_2d(tl_VFe(LBi:UBi,LBj:UBj),
!     *                     LBi,UBi,LBj,UBj,
!     *                     Istr_E,Iend_E,Jstr_E,Jend_E,
!     *                     comm2d, stride, 
!     *                     rankleft,rankright,ranktop,rankbottom) 
!#endif /* MPI */

        DO j=JstrV,JendV+1
         DO i=IstrV,IendV+1

! NL:            cff=visc2_p(i,j)*0.125d0*(Hz(i-1,j  ,k)+Hz(i,j  ,k)+        
! NL:     *                                Hz(i-1,j-1,k)+Hz(i,j-1,k))*       
! NL:     *          (dyodx_p(i,j)*               
! NL:     *           ((ody(i  ,j-1)+ody(i  ,j))*v(i  ,j,k,nrhs)- 
! NL:     *            (ody(i-1,j-1)+ody(i-1,j))*v(i-1,j,k,nrhs))+ 
! NL:     *           dxody_p(i,j)*                   
! NL:     *           ((odx(i-1,j  )+odx(i,j  ))*u(i,j  ,k,nrhs)-    
! NL:     *            (odx(i-1,j-1)+odx(i,j-1))*u(i,j-1,k,nrhs)))

            tmp=visc2_p(i,j)*0.125d0*dyodx_p(i,j)
            tmp1=tmp*(ody(i  ,j-1)+ody(i  ,j))
            tmp2=tmp*(ody(i-1,j-1)+ody(i-1,j))

            tmp=visc2_p(i,j)*0.125d0*dxody_p(i,j)
            tmp3=tmp*(odx(i-1,j  )+odx(i,j  ))
            tmp4=tmp*(odx(i-1,j-1)+odx(i,j-1))

! rewrite NL using defs above:
!           cff=(Hz(i-1,j  ,k)+Hz(i,j  ,k)+Hz(i-1,j-1,k)+Hz(i,j-1,k))*
!    *          (tmp1*v(i  ,j,k,nrhs)-tmp2*v(i-1,j,k,nrhs)+
!    *           tmp3*u(i,j  ,k,nrhs)-tmp4*u(i,j-1,k,nrhs))

            tl_cff=(tl_Hz(i-1,j  ,k)+tl_Hz(i,j  ,k)+
     *              tl_Hz(i-1,j-1,k)+tl_Hz(i,j-1,k))*
     *             (tmp1*v(i  ,j,k,nrhs)-tmp2*v(i-1,j,k,nrhs)+
     *              tmp3*u(i,j  ,k,nrhs)-tmp4*u(i,j-1,k,nrhs)) +
     *            (Hz(i-1,j  ,k)+Hz(i,j  ,k)+Hz(i-1,j-1,k)+Hz(i,j-1,k))*
     *             (tmp1*tl_v(i  ,j,k,nrhs)-tmp2*tl_v(i-1,j,k,nrhs)+
     *              tmp3*tl_u(i,j  ,k,nrhs)-tmp4*tl_u(i,j-1,k,nrhs)) 

            tl_cff=tl_cff*pmask(i,j)

            tl_UFe(i,j)=dx_p(i,j)*dx_p(i,j)*tl_cff
            tl_VFx(i,j)=dy_p(i,j)*dy_p(i,j)*tl_cff
          END DO
        END DO
!
! Time-step harmonic, S-surfaces viscosity term. Notice that momentum
! at this stage is HzU and HzV and has m2/s units. Add contribution for
! barotropic forcing terms.
!
        DO j=JstrU,JendU
          DO i=IstrU,IendU
! roms, oArea_u(i,j): cff=0.25_r8*(pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))
            tl_cff1=0.5d0*((ody(i-1,j)+ody(i,j))*           
     *                  (tl_UFx(i,j  )-tl_UFx(i-1,j))+
     *                  (odx(i-1,j)+odx(i,j))*   
     *                  (tl_UFe(i,j+1)-tl_UFe(i  ,j)))
            tl_cff2=dt*oArea_u(i,j)*tl_cff1
            tl_rufrc(i,j)=tl_rufrc(i,j)+tl_cff1
            tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)+tl_cff2
          END DO
        END DO

        DO j=JstrV,JendV
          DO i=IstrV,IendV
! roms, oArea_v=cff=0.25d0*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
            tl_cff1=0.5d0*((ody(i,j-1)+ody(i,j))*                
     *                  (tl_VFx(i+1,j)-tl_VFx(i,j  ))-              
     *                  (odx(i,j-1)+odx(i,j))*                
     *                  (tl_VFe(i  ,j)-tl_VFe(i,j-1)))
            tl_cff2=dt*oArea_v(i,j)*tl_cff1
            tl_rvfrc(i,j)=tl_rvfrc(i,j)+tl_cff1
            tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)+tl_cff2
          END DO
        END DO
      END DO K_LOOP

# endif /* UV_VIS2 */

      RETURN
      END

#endif /* TANGENT */
