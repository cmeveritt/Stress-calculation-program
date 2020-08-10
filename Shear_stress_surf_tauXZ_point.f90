! A subroutine for calculating the stresses from a pont shear traction stress in hte xz direction
    subroutine Shear_stress_Surf_tauXZ_point(i,j,jl)
    implicit none
        integer     NN, NNX
        integer      i, j, k, l
        integer         ik, jl, jljl, ikik, iki, jlj
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
                        
                    Psmall  = (tau_xzl_min(k)     *(NNZ-abs(ikik))     + tau_xzl_min(k+step)     *abs(ikik))   /NNZ   ! Line load
                    Psmall1 = (tau_xzl(k,l)       *(NNZ-abs(ikik))     + tau_xzl(k+step,l)       *abs(ikik))   /NNZ   ! Point loas in x
                    Psmall2 = (tau_xzl(k,l+stepy) *(NNZ-abs(ikik))     + tau_xzl(k+step,l+stepy) *abs(ikik))   /NNZ   ! Point load in x one y step away
                    Psmall3 = (Psmall1            *(NNZ-abs(jljl))     + Psmall2                 *abs(jljl))   /NNZ   ! Point load in y given x-values
                    Psmall4 = Psmall3-Psmall                                                                ! The pointloas-line load 
                    Psmall4 = (Psmall4)/(2.0*pi)
                    Psmall4 =-Psmall4
                           
                    iki = -ikik  +NZ_cvot*ik
                    jlj = -jljl  +NZ_cvot*jl
                
                    r2=iki**2+jlj**2   

                    rho=(r2)**0.5          ! The total distance [m]  ! z=0.0
                    rho3=rho**3
                    rho5=rho**5
               
                    if ( iki .EQ. 0) THen !If we're at the line where the pressure is applied. Here the line solution is not applicable
                        If ( jlj .EQ. 0) THEN ! If we're at the nod where the shear traction  is applied---------------------------
                                ! Nothing is added
                            tau_xzs(i,j)=tau_xzs(i,j)+tau_xzl(i,j)
                        else
                            !sig_xs=0 since ik=0  
                            !sig_ys=0 since ik=0
                            ! sig_zs=0 since z=0
                 
                            tau_xys(i,j)=tau_xys(i,j)  +Psmall4 * (1.0-2.0*0.3) * (-jlj/rho3)
                            !tau_xzs since z=0
                            !tau_yzs since z=0
                        endif
                
                        ! The stresses from inside the simulation area  ----------------------------------------------
                    elseIF( abs(Psmall4) .GT. 0.00005)THEN   
                    
                        sig_xs(i,j)=sig_xs(i,j)    +Psmall4     *(-3.0*iki**3/rho5       +(1.0-2.0*0.3) * (iki/rho3  - 3.0*iki/rho3  + iki**3/rho5     + 2.0*iki**3/rho5))    !Subtract the lineload since a line load not yields stresses in the perpendiculat direction on the surface. 
                        sig_ys(i,j)=sig_ys(i,j)    +Psmall4     *(-3.0*iki*jlj**2/rho5   +(1.0-2.0*0.3) * (iki/rho3  -     iki/rho3  + iki*jlj**2/rho5 + 2.0*iki*jlj**2/rho5))
                        ! sig_zs=0 since z=0
                 
                        tau_xys(i,j)=tau_xys(i,j)  +Psmall4     *(-3.0*iki**2*jlj/rho5   +(1.0-2.0*0.3) * (-jlj/rho3                   + iki**2*jlj/rho5 + 2.0*iki**2*jlj/rho5))
                        !tau_xzs since z=0
                        !tau_yzs since z=0
                     
                    endif
                enddo
            enddo
        else
                        Psmall4 = (tau_xzl(k,l)-tau_xzl_min(k))/(2*pi)
                        Psmall4 = -Psmall4
                        IF( abs(Psmall4) .GT. 0.0005)THEN   
                            iki=ik
                            jlj=jl
                            r2=iki**2+jlj**2   
                            rho=(r2)**0.5         
                            rho3=rho**3
                            rho5=rho**5
                            
                            sig_xs(i,j)=sig_xs(i,j)    +Psmall4     *(-3.0*ik**3/rho5      +(1.0-2.0*0.3) * (ik/rho3  - 3.0*ik/rho3  + ik**3/rho5    + 2.0*ik**3/rho5))    !Subtract the lineload since a line load not yields stresses in the perpendiculat direction on the surface. 
                            sig_ys(i,j)=sig_ys(i,j)    +Psmall4     *(-3.0*ik*jl**2/rho5   +(1.0-2.0*0.3) * (ik/rho3  -     ik/rho3  + ik*jl**2/rho5 + 2.0*ik*jl**2/rho5))
                            ! sig_zs=0 since z=0
                 
                            tau_xys(i,j)=tau_xys(i,j)  +Psmall4     *(-3.0*ik**2*jl/rho5   +(1.0-2.0*0.3) * (-jl/rho3                  + ik**2*jl/rho5 + 2.0*ik**2*jl/rho5))
                            !tau_xzs only from lubrication shear
                            !tau_yzs only from lubrication shear
                        endif
                        
        endif
        
        
        
     enddo
     return
    end
    