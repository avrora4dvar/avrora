#include "cppdefs.h"

#ifdef ADJOINT
      SUBROUTINE ad_uv3dmix2
      USE mod_ocean

      integer :: i,j,k
      real(8) :: ad_cff,ad_cff1,ad_cff2
      real(8) :: tmp,tmp1,tmp2,tmp3,tmp4

      ad_cff=0.d0
      ad_cff1=0.d0
      ad_cff2=0.d0

# ifdef UV_VIS2
      K_LOOP : DO k=1,N
! ADJ to: Time-step harmonic, S-surfaces viscosity term. Notice that momentum
! at this stage is HzU and HzV and has m2/s units. Add contribution for
! barotropic forcing terms.
!
        DO j=JstrV,JendV
          DO i=IstrV,IendV
! roms, oArea_v=cff=0.25d0*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))

! tl line 4 (recur):   tl_v(i,j,k,nnew)=tl_v(i,j,k,nnew)+tl_cff2
          ad_cff2=ad_cff2+ad_v(i,j,k,nnew)

! tl line 3 (recur):            tl_rvfrc(i,j)=tl_rvfrc(i,j)+tl_cff1
          ad_cff1=ad_cff1+ad_rvfrc(i,j)

! tl line 2:            tl_cff2=dt*oArea_v(i,j)*tl_cff1
          ad_cff1=ad_cff1+dt*oArea_v(i,j)*ad_cff2
          ad_cff2=0.d0

! tl line 1:            tl_cff1=0.5d0*((ody(i,j-1)+ody(i,j))*                
! tl line 1:     *                  (tl_VFx(i+1,j)-tl_VFx(i,j  ))-              
! tl line 1:     *                  (odx(i,j-1)+odx(i,j))*                
! tl line 1:     *                  (tl_VFe(i  ,j)-tl_VFe(i,j-1)))

           tmp1=0.5d0*(ody(i,j-1)+ody(i,j))*ad_cff1
           tmp2=0.5d0*(odx(i,j-1)+odx(i,j))*ad_cff1

           ad_VFx(i+1,j)=ad_VFx(i+1,j)+tmp1
           ad_VFx(i,j  )=ad_VFx(i,j  )-tmp1
           ad_VFe(i  ,j)=ad_VFe(i  ,j)-tmp2
           ad_VFe(i,j-1)=ad_VFe(i,j-1)+tmp2

           ad_cff1=0.d0

          END DO
        END DO

        DO j=JstrU,JendU
          DO i=IstrU,IendU
! roms, oArea_u(i,j): cff=0.25_r8*(pm(i-1,j)+pm(i,j))*(pn(i-1,j)+pn(i,j))
! tl line 4:            tl_u(i,j,k,nnew)=tl_u(i,j,k,nnew)+tl_cff2
           ad_cff2=ad_cff2+ad_u(i,j,k,nnew)

! tl line 3:            tl_rufrc(i,j)=tl_rufrc(i,j)+tl_cff1
           ad_cff1=ad_cff1+ad_rufrc(i,j)

! tl line 2:            tl_cff2=dt*oArea_u(i,j)*tl_cff1
           ad_cff1=ad_cff1+dt*oArea_u(i,j)*ad_cff2
           ad_cff2=0.d0

! tl line 1:            tl_cff1=0.5d0*((ody(i-1,j)+ody(i,j))*           
! tl line 1:     *                  (tl_UFx(i,j  )-tl_UFx(i-1,j))+
! tl line 1:     *                  (odx(i-1,j)+odx(i,j))*   
! tl line 1:     *                  (tl_UFe(i,j+1)-tl_UFe(i  ,j)))

           tmp1=0.5d0*(ody(i-1,j)+ody(i,j))*ad_cff1
           tmp2=0.5d0*(odx(i-1,j)+odx(i,j))*ad_cff1
 
           ad_UFx(i,j  )=ad_UFx(i,j  )+tmp1
           ad_UFx(i-1,j)=ad_UFx(i-1,j)-tmp1
           ad_UFe(i,j+1)=ad_UFe(i,j+1)+tmp2
           ad_UFe(i,j  )=ad_UFe(i,j  )-tmp2

           ad_cff1=0.d0

          END DO
        END DO

!
!  ADJ to: Compute flux-components of the horizontal divergence of the stress
!  tensor (m5/s2) in XI- and ETA-directions.
!
!
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

! tl line 4:            tl_VFx(i,j)=dy_p(i,j)*dy_p(i,j)*tl_cff
            ad_cff=ad_cff+dy_p(i,j)*dy_p(i,j)*ad_VFx(i,j)
            ad_VFx(i,j)=0.d0

! tl line 3:            tl_UFe(i,j)=dx_p(i,j)*dx_p(i,j)*tl_cff
            ad_cff=ad_cff+dx_p(i,j)*dx_p(i,j)*ad_UFe(i,j)
            ad_UFe(i,j)=0.d0

! tl line 2:            tl_cff=tl_cff*pmask(i,j)
            ad_cff=ad_cff*pmask(i,j)

! tl line 1:        tl_cff=(tl_Hz(i-1,j  ,k)+tl_Hz(i,j  ,k)+
! tl line 1:     *          tl_Hz(i-1,j-1,k)+tl_Hz(i,j-1,k))*
! tl line 1:     *         (tmp1*v(i  ,j,k,nrhs)-tmp2*v(i-1,j,k,nrhs)+
! tl line 1:     *          tmp3*u(i,j  ,k,nrhs)-tmp4*u(i,j-1,k,nrhs)) +
! tl line 1:     *         (Hz(i-1,j  ,k)+Hz(i,j  ,k)+Hz(i-1,j-1,k)+Hz(i,j-1,k))*
! tl line 1:     *         (tmp1*tl_v(i  ,j,k,nrhs)-tmp2*tl_v(i-1,j,k,nrhs)+
! tl line 1:     *          tmp3*tl_u(i,j  ,k,nrhs)-tmp4*tl_u(i,j-1,k,nrhs)) 

            tmp=( tmp1*v(i  ,j,k,nrhs)-tmp2*v(i-1,j,k,nrhs)
     *           +tmp3*u(i,j  ,k,nrhs)-tmp4*u(i,j-1,k,nrhs) )*ad_cff

            ad_Hz(i-1,j  ,k)=ad_Hz(i-1,j  ,k)+tmp
            ad_Hz(i,j  ,k)=ad_Hz(i,j  ,k)+tmp
            ad_Hz(i-1,j-1,k)=ad_Hz(i-1,j-1,k)+tmp
            ad_Hz(i,j-1,k)=ad_Hz(i,j-1,k)+tmp

            tmp=(Hz(i-1,j,k)+Hz(i,j,k)+Hz(i-1,j-1,k)+Hz(i,j-1,k))
     *          *ad_cff

            ad_v(i  ,j,k,nrhs)=ad_v(i  ,j,k,nrhs)+tmp*tmp1
            ad_v(i-1,j,k,nrhs)=ad_v(i-1,j,k,nrhs)-tmp*tmp2
            ad_u(i,j  ,k,nrhs)=ad_u(i,j  ,k,nrhs)+tmp*tmp3
            ad_u(i,j-1,k,nrhs)=ad_u(i,j-1,k,nrhs)-tmp*tmp4

            ad_cff=0.d0

          END DO
        END DO

        DO j=JstrV-1,JendV
!        DO i=min(1,IstrV-1),IendV ! new in v3.2, v3.1 had i=IstrV,IendV
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
!                        -tmp3*v(i,j+1,k,nrhs)+tmp4*v(i,j+1,k,nrhs) )

! tl line 3:          tl_VFe(i,j)=dx(i,j)*dx(i,j)*tl_cff
          ad_cff=ad_cff+dx(i,j)*dx(i,j)*ad_VFe(i,j)
          ad_VFe(i,j)=0.d0
 
! tl line 2:          tl_UFx(i,j)=dy(i,j)*dy(i,j)*tl_cff
          ad_cff=ad_cff+dy(i,j)*dy(i,j)*ad_UFx(i,j)
          ad_UFx(i,j)=0.d0

! tl line 1:          tl_cff= tl_Hz(i,j,k)*
! tl line 1:     *             ( tmp1*u(i+1,j,k,nrhs)-tmp2*u(i  ,j,k,nrhs)
! tl line 1:     *              -tmp3*v(i,j+1,k,nrhs)+tmp4*v(i,j+1,k,nrhs) )
! tl line 1:     *           +Hz(i,j,k)*
! tl line 1:     *             ( tmp1*tl_u(i+1,j,k,nrhs)-tmp2*tl_u(i  ,j,k,nrhs)
! tl line 1:     *              -tmp3*tl_v(i,j+1,k,nrhs)+tmp4*tl_v(i,j ,k,nrhs) )
          ad_Hz(i,j,k)=ad_Hz(i,j,k)+
     *      ( tmp1*u(i+1,j,k,nrhs)-tmp2*u(i  ,j,k,nrhs)
     *       -tmp3*v(i,j+1,k,nrhs)+tmp4*v(i,j+1,k,nrhs) )*ad_cff

          tmp=Hz(i,j,k)*ad_cff

          ad_u(i+1,j,k,nrhs)=ad_u(i+1,j,k,nrhs)+tmp1*tmp
          ad_u(i  ,j,k,nrhs)=ad_u(i  ,j,k,nrhs)-tmp2*tmp
          ad_v(i,j+1,k,nrhs)=ad_v(i,j+1,k,nrhs)-tmp3*tmp
          ad_v(i,j  ,k,nrhs)=ad_v(i,j  ,k,nrhs)+tmp4*tmp

          ad_cff=0.d0

         END DO
        END DO


      END DO K_LOOP

# endif /* UV_VIS2 */

      RETURN
      END

#endif /* ADJOINT */
