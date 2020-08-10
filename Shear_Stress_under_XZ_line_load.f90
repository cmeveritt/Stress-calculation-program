!A subroutine for calculating the stresses dowin in the material based on the line shear load tau_xzmin
    Subroutine Shear_Stress_under_XZ_line_load(hh,i,j)
        implicit none
        integer      NN, NNX
        integer      i, j, k, l
        integer      ik, jl, js, jll
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
        include     'inc_CurrentP.h'
        include     'inc_Minimum_values.h'
        include     'inc_StressParam.h'
        include     'inc_Stress_under.h'
    
        pi=3.14159265
        NN=(NY+1)/2
        NNX=(NX+1)/2
        NNZ=(NZ_cvot-1)/2      
    
        if( hh .LT. Nr_z*NZ_cvot+2 )then
            z=(hh-1.0)/NZ_cvot 
        else
            z=(hh-Nr_z*(NZ_cvot-1)-1)                    ! To have same scale on z as ik and jl 
        endif
        
    ! At the surface a line load, Pmin, will not couse any stresses exept than where its applied. Point loads will on the other hand.  
            Do k=NXstart,NXstop
                ik=i-k
                
                if (abs(ik) .LT. Nr_xy ) then ! If we're close to the line where the pressure is applied
                        
                    DO ikik=-NNZ,NNZ
                        if (ikik .GT. 0)then 
                            step=1 
                        elseif( ikik .EQ. 0) then 
                            step=0
                        else 
                            step=-1
                        endif
                    
                        Psmall=(tau_xzl_min(k)*(NNZ-abs(ikik))+tau_xzl_min(k+step)*abs(ikik))/NNZ
                        Psmall=-Psmall
                        IF( abs(Psmall) .GT. 0.0001)THEN            ! To gain some speed
                            
                            z20  =         NZ_cvot*z
                            iki  = -ikik + NZ_cvot*ik
                            base = (iki**2+(z20)**2)**2
                            Psmall      = -2.0*Psmall/(pi*base)
                            
                            sig_x(i,hh) = sig_x(i,hh)    +Psmall*iki**3
                            sig_y(i,hh) = sig_y(i,hh)    +0.3*Psmall*(iki**3-iki*z20**2)  
                            sig_z(i,hh) = sig_z(i,hh)    +Psmall*iki*z20**2
                
                            tau_xz(i,hh)=tau_xz(i,hh)    +Psmall*iki**2*z20
                            !tau_xy(i,hh)=tau_xy(i,hh) 
                            !tau_yz(i,hh)=tau_yz(i,hh)   
                        endif

                    enddo
     
                else
                    IF( abs(tau_xzl_min(k)) .GT. 0.0001)THEN 
                        base        = (ik**2+z**2)**2    
                        Psmall      = -2.0*tau_xzl_min(k)/(pi*base)
                        Psmall=-Psmall
                    
                        sig_x(i,hh) = sig_x(i,hh)   + Psmall*ik**3
                        sig_y(i,hh) = sig_y(i,hh)   + 0.3*Psmall*(ik**3-ik*z**2)   
                        sig_z(i,hh) = sig_z(i,hh)   + Psmall*ik*z**2
                
                        tau_xz(i,hh)=tau_xz(i,hh)   + Psmall*ik**2*z
                        !tau_xy(i,hh)=tau_xy(i,hh)  
                        !tau_yz(i,hh)=tau_yz(i,hh)  
                    endif

                endif
            enddo
            
            return
    end
    