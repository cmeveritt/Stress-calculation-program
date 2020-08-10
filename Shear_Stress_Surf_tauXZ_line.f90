! A subroutine for calculating the stresses from a line shear traction stress in hte xz direction
    subroutine Shear_Stress_Surf_tauXZ_line(i)
    implicit none
        integer      NN, NNX
        integer      i, j, k, l
        Integer         ik, jl
        real*8        min, pi, r2
        real*8        sigxx, sigyy, rho5, rho, rho3
        integer        step, NNz, ikik, iki
        real        Psmall, Psmall4
        save        /Stress_surf/
        include     'inc_Grid.h'
        include     'inc_Nodes_no_pressure.h'
        include     'inc_Stress_surf.h'
        include     'inc_CurrentP.h'
        include     'inc_CurrentH.h'
        include     'inc_Minimum_values.h'
        include     'inc_StressParam.h'
        
        
    pi=3.14159265
    NN=(NY+1)/2
    NNX=(NX+1)/2
    NNZ=(NZ_cvot-1)/2   
    
        !Sum up the line load at line i 
        sigxx=0
        sigyy=0
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
                        
                        iki  = -ikik+NZ_cvot*ik
                        if( iki .NE. 0) then
                            Psmall=(tau_xzl_min(k)*(NNZ-abs(ikik))+tau_xzl_min(k+step)*abs(ikik))/NNZ
                            Psmall=-Psmall
                        
                            IF( abs(Psmall) .GT. 0.0001)THEN            ! To gain some speed
                            
                                sigxx    = sigxx -2.0*Psmall/(pi*iki)                    !* (ik**3)       /base
                                sigyy    = sigyy + 0.3*(-2.0*Psmall/(pi*iki)+0)
                            endif
                        endif
                    enddo
                else
                    iki  = ik
                        
                    Psmall=(tau_xzl_min(k))
                    Psmall=-Psmall
                        
                    IF( abs(Psmall) .GT. 0.0001)THEN            ! To gain some speed
                            
                        sigxx    = sigxx -2.0*Psmall/(pi*iki)                    !* (ik**3)       /base
                        sigyy    = sigyy + 0.3*(-2.0*Psmall/(pi*iki)+0)
                    endif
                        
                endif
                
                    
        enddo
        
        !Distribute the values to al nodes at i
        DO j=1,NN            ! Interested in the whole area
            sig_xs(i,j)  = sig_xs(i,j)  + sigxx  
            sig_ys(i,j)  = sig_ys(i,j)  + sigyy
            !tau_xzs(i,j) = tau_xzs(i,j) + tau_xzl_min(i) Added directly on each node in the point calculations instead
        enddo
        
        return
    end
    