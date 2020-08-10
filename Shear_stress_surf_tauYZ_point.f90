! A subroutine for calculating the stresses from a pont shear traction stress in hte xz direction
    subroutine Shear_stress_Surf_tauYZ_point(i,j,jl)
    implicit none
        integer      NN, NNX
        integer      i, j, k, l
        integer         ik, jl, jljl, ikik, iki, jlj, ikir, jljr
        real*8        min, pi, r2
        real*8        sigxx, sigyy, rho5, rho, rho3
        real*8      psmall1, psmall2, psmall3
        integer        step, stepy, NNz
        real        Psmall, Psmall4
        save        /Stress_surf/
        include     'inc_Grid.h'
        include     'inc_Nodes_no_pressure.h'
        include     'inc_Stress_surf.h'
        include     'inc_CurrentP.h'
        include     'inc_CurrentH.h'
        include     'inc_Minimum_values.h'
        include     'inc_StressParam.h'
        include     'inc_load_surf.h'
        
    pi=3.14159265
    NN=(NY+1)/2
    NNX=(NX+1)/2
    NNZ=(NZ_cvot-1)/2   
    
     DO K=NXstart,NXstop
        l=j-jl
  
        IF( l .GE. NY)    l=2    ! If outside the area
        IF( l .LT. 2)     l=2
                
        ik=i-k
                    
        IF (abs(ik) .LT. Nr_xy .AND. abs(jl) .LT. Nr_xy) then ! If we're close to the line where the pressure is applied         
            DO jljl=-NNZ,NNZ
                DO ikik=-NNZ,NNZ
                    if (ikik .GT. 0)then 
                        step=1 
                    elseif( ikik .EQ. 0) then 
                        step=0
                    else 
                        step=-1
                    endif
                    if (jljl .GT. 0)then 
                        stepy=1 
                    elseif( jljl .EQ. 0) then 
                        stepy=0
                    else 
                        stepy=-1
                    endif
                        
                    ! No line load is taken into account for the shear in yz. 
                    Psmall1 = (tau_yzl(k,l)       *(NNZ-abs(ikik))     + tau_yzl(k+step,l)       *abs(ikik))   /NNZ   ! Point loas in x
                    Psmall2 = (tau_yzl(k,l+stepy) *(NNZ-abs(ikik))     + tau_yzl(k+step,l+stepy) *abs(ikik))   /NNZ   ! Point load in x one y step away
                    Psmall3 = (Psmall1            *(NNZ-abs(jljl))     + Psmall2                 *abs(jljl))   /NNZ   ! Point load in y given x-values
                    Psmall4 = Psmall3                                                               ! The pointloas-line load 
                    Psmall4 = (Psmall4)/(2.0*pi)
                    Psmall4 = -Psmall4
                        
                        iki = -ikik  +NZ_cvot*ik
                        jlj = -jljl  +NZ_cvot*jl
                
                        r2=iki**2+jlj**2   

                        rho=(r2)**0.5          ! The total distance [m]  ! z=0.0
                        rho3=rho**3
                        rho5=rho**5
                        
                        !Rotating coordinate system
                        ikir = jlj
                        jljr = -iki
                        
                        if ( iki .EQ. 0) THen !If we're at the line where the pressure is applied. Here the line solution is not applicable
                            If ( jlj .EQ. 0) THEN ! If we're at the nod where the shear traction  is applied---------------------------
                                ! Nothing is added
                                tau_yzs(i,j)=tau_yzs(i,j)+tau_yzl(i,j)
                            else
                                !sig_xs=0 since ik=0  
                                !sig_ys=0 since ik=0
                                ! sig_zs=0 since z=0
                                
                                tau_xys(i,j)=tau_xys(i,j)  -Psmall4  * (1.0-2.0*0.3) * (-jljr/rho3)     !Minus since tau_yz instead of tau_xz
                                !tau_xzs since z=0
                                !tau_yzs since z=0
                            endif
                
                        ! The stresses from inside the simulation area  ----------------------------------------------
                        else IF( abs(Psmall4) .GT. 0.00005) then
                    
                            sig_ys(i,j)=sig_ys(i,j)    +Psmall4   *(-3.0*ikir**3/rho5        +(1.0-2.0*0.3) * (ikir/rho3  - 3.0*ikir/rho**3  + ikir**3/rho5      + 2.0*ikir**3/rho5))   !Subtract the lineload since a line load not yields stresses in the perpendiculat direction on the surface. 
                            sig_xs(i,j)=sig_xs(i,j)    +Psmall4   *(-3.0*ikir*jljr**2/rho5   +(1.0-2.0*0.3) * (ikir/rho3  -     ikir/rho**3  + ikir*jljr**2/rho5 + 2.0*ikir*jljr**2/rho5))
                            ! sig_zs=0 since z=0
                 
                            tau_xys(i,j)=tau_xys(i,j)  -Psmall4   *(-3.0*ikir**2*jljr/rho5   +(1.0-2.0*0.3) * (-jljr/rho3                    + ikir**2*jljr/rho5 + 2.0*ikir**2*jljr/rho5))    !Minus since tau_yz instead of tau_xz
                            !tau_xzs only from lubrication shear
                            !tau_yzs only from lubrication shear
                        endif
                    
                enddo
            enddo
        else
                        Psmall4 = (tau_yzl(k,l)-0)/(2*pi)
                        Psmall4 = - Psmall4
                        IF( abs(Psmall4) .GT. 0.0005)THEN   
                            iki=ik
                            jlj=jl
                            !Rotating coordinate system
                            ikir = jlj
                            jljr = -iki
                            
                            r2=iki**2+jlj**2   
                            rho=(r2)**0.5         
                            rho3=rho**3
                            rho5=rho**5
                            
                            sig_ys(i,j)=sig_ys(i,j)    +Psmall4   *(-3.0*ikir**3/rho5        +(1.0-2.0*0.3) * (ikir/rho3  - 3.0*ikir/rho**3  + ikir**3/rho5      + 2.0*ikir**3/rho5))   !Subtract the lineload since a line load not yields stresses in the perpendiculat direction on the surface. 
                            sig_xs(i,j)=sig_xs(i,j)    +Psmall4   *(-3.0*ikir*jljr**2/rho5   +(1.0-2.0*0.3) * (ikir/rho3  -     ikir/rho**3  + ikir*jljr**2/rho5 + 2.0*ikir*jljr**2/rho5))
                            ! sig_zs=0 since z=0
                 
                            tau_xys(i,j)=tau_xys(i,j)  -Psmall4   *(-3.0*ikir**2*jljr/rho5   +(1.0-2.0*0.3) * (-jljr/rho3                    + ikir**2*jljr/rho5 + 2.0*ikir**2*jljr/rho5))    !Minus since tau_yz instead of tau_xz
                            !tau_xzs only from lubrication shear
                            !tau_yzs only from lubrication shear
                        endif
                        
        endif
        
        
        
     enddo
     return
    end
    