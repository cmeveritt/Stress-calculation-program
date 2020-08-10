!A subroutine for calculating the stresses dowin in the material based on the point shear loads on the surface
    Subroutine Shear_Stress_under_XZ_point_load(hh,i,j,jl)
        implicit none
        integer      NN, NNX
        integer      i, j, k, l
        integer      ik, jl, js
        integer      ikik, jljl, step, stepy, iki, jlj
        real*8         Psmall1, Psmall2, Psmall3, Psmall4, Psmall, z20
        real*8         sigxx, sigzz, ik3, ikjl2
        real*8         rho5, rho3, base
        real*8         rhoz2, rho3z2, rho2z3, pi
        real*8         rho, r2, z
        integer      hh, NNZ
        save        /Stress_under/
        include     'inc_Grid.h'
        include     'inc_Nodes_no_pressure.h'
        include     'inc_Stress_surf.h'
        include     'inc_CurrentP.h'
        include     'inc_Minimum_values.h'
        include     'inc_StressParam.h'
        include     'inc_Stress_under.h'
        include     'inc_load_surf.h'
    
        pi=3.14159265
        NN=(NY+1)/2
        NNX=(NX+1)/2
        NNZ=(NZ_cvot-1)/2      
    
        if( hh .LT. Nr_z*NZ_cvot+2 )then
            z=(hh-1.0)/NZ_cvot 
        else
            z=(hh-Nr_z*(NZ_cvot-1)-1)                    ! To have same scale on z as ik and jl 
        endif
        
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
                        Psmall4 = -Psmall4
                        
                        IF( abs(Psmall4) .GT. 0.0005)THEN            ! To gain some speed
                            z20 =         Nz_cvot*z
                            iki = -ikik  +NZ_cvot*ik
                            jlj = -jljl  +NZ_cvot*jl
                
               
                            r2=iki**2+jlj**2                  ! The square of the  horizontal distance from the point [m**2]
                                IF(r2 .EQ. 0.0) r2=0.5**2       ! If straight underneath the pressure, assume we're half a distance away. 
                                rho=(iki**2+jlj**2+z20**2)**0.5     ! The total distance [m]
                        
                                rho5=1.0/rho**5
                                rho3=1.0/rho**3
                                rhoz2=1.0/(rho*(rho+z20))**2
                                rho3z2=1.0/(rho**3*(rho+z20)**2)
                                rho2z3=1.0/(rho**2*(rho+z20)**3)
                                ik3=iki**3
                                ikjl2 = iki*jlj**2
                                ! 1-2*0.3=0.4
                                sig_x(i,hh)=sig_x(i,hh)   +Psmall4* (-3.0*ik3   *rho5 + 0.4 * (iki*rho3 - 3.0*iki*rhoz2 + ik3   *rho3z2 + 2.0*ik3   *rho2z3))
                                sig_y(i,hh)=sig_y(i,hh)   +Psmall4* (-3.0*ikjl2 *rho5 + 0.4 * (iki*rho3 -     iki*rhoz2 + ikjl2 *rho3z2 + 2.0*ikjl2 *rho2z3))
                        
                                sig_z(i,hh)=sig_z(i,hh)   -Psmall4*   3.0*iki*z20**2*rho5
                        
                                ! No shear xy stresses due to symmetry
                                !tau_xy(i,hh)=tau_xy(i,hh) +(tau_xzs(k,l)-tau_xzl_min(k))/(2*pi)* ( 3*ik**2*jl/rho5 + (1-2*0.3) * (        -   jl/rhoz2 + ik**2*jl/rho3z2 + 2*ik**2*jl/rho2z3))
                                !tau_yz(i,hh)=tau_yz(i,hh) -Psmall4*   3.0*iki**2 *z20*rho5
                                 tau_xz(i,hh)=tau_xz(i,hh) -Psmall4*   3.0*iki*jlj*z20*rho5
                                if (sig_x(i,hh) .GT. 0.5) then ! Theses is verry high values which are obtained in the exit region
                                    sig_x(1,1)=sig_x(i,hh)
                                endif
                                
                        ENDIF
                    ENDDO
                    ENDDO
                      
                ELSE
                    
                    Psmall4 = (tau_xzl(k,l))/(2.0*pi)
                    Psmall4 = - Psmall4
                    
                    IF( abs(Psmall4) .GT. 0.0005)THEN            ! To gain some speed
                    r2=ik**2+jl**2                  ! The square of the  horizontal distance from the point [m**2]
                    IF(r2 .EQ. 0.0) r2=0.5**2       ! If straight underneath the pressure, assume we're half a distance away. 
                        rho=(ik**2+jl**2+z**2)**0.5     ! The total distance [m]
                        
                        rho5=1.0/rho**5
                        rho3=1.0/rho**3
                        rhoz2=1.0/((rho*(rho+z))**2)
                        rho3z2=1.0/(rho**3*(rho+z)**2)
                        rho2z3=1.0/(rho**2*(rho+z)**3)
                        Psmall4 = (tau_xzs(k,l)-tau_xzl_min(k))/(2.0*pi)
                        ik3=ik**3
                        ikjl2=ik*jl**2
                        
                        sig_x(i,hh)=sig_x(i,hh)   +Psmall4* (-3.0*ik3*rho5   + 0.4 * (ik*rho3 - 3.0*ik*rhoz2 + ik3   *rho3z2 + 2.0*ik3   *rho2z3))
                        sig_y(i,hh)=sig_y(i,hh)   +Psmall4* (-3.0*ikjl2*rho5 + 0.4 * (ik*rho3 -     ik*rhoz2 + ikjl2 *rho3z2 + 2.0*ikjl2 *rho2z3))
                        
                        sig_z(i,hh)=sig_z(i,hh)   -Psmall4*  3.0*ik*z**2*rho5
                        
                        ! No shear xy stresses due to symmetry
                        !tau_xy(i,hh)=tau_xy(i,hh) +(tau_xzs(k,l)-tau_xzl_min(k))/(2*pi)* ( 3*ik**2*jl/rho5 + (1-2*0.3) * (        -   jl/rhoz2 + ik**2*jl/rho3z2 + 2*ik**2*jl/rho2z3))
                        !tau_yz(i,hh)=tau_yz(i,hh) -(tau_xzs(k,l)-tau_xzl_min(k))/(2.0*pi)*  3.0*ik**2*z/rho5
                        tau_xz(i,hh)=tau_xz(i,hh) -Psmall4 *  3.0*ik*jl*z*rho5
                        if (sig_x(i,hh) .GT. 0.5) then ! Theses is verry high values which are obtained in the exit region
                            sig_x(1,1)=sig_x(i,hh)
                        endif
                        
                    ENDIF
                    
                ENDIF
            ENDDO 
            
            return
    end
    