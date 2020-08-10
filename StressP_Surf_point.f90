! Subroutine to fint the stresses at the surface due to a normal load
    SUBROUTINE StressP_Surf_point(i,j,jl)
        implicit none
        integer      NN, NNX, NNz
        integer     jljl, ikik, jlj, iki 
        integer      i, j, k, l
        integer      ik, jl
        integer     stepy, step
        real*8         Psmall4, Psmall, Psmall3, Psmall2, Psmall1
        real*8      min, pi, r2, r4
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
        NNX=(NX+1)/2            !Its enough distance to get converged stress values. Moste of the stress levels are due to the line loads
        NNZ=(NZ_cvot-1)/2      ! The extra number of nodes due to refinement due to the cvote of refinment NZ_cvot
        
        DO K=NXstart,NXstop

            l=j-jl
  
            IF( l .GE. NY)    l=2    ! If outside the area
            IF( l .LT. 2)     l=2
                
            ik=i-k
                
            IF (abs(ik) .LT. Nr_xy .AND. abs(jl) .LT. Nr_xy) then ! If we're close to the point where the stresses are evaluated
                        
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
                        
                        Psmall  = (Pmin(k)     *(NNZ-abs(ikik))     + Pmin(k+step)      *abs(ikik))   /NNZ   ! Line load
                        Psmall1 = (P(k,l)      *(NNZ-abs(ikik))     + P(k+step,l)       *abs(ikik))   /NNZ   ! Point loas in x
                        Psmall2 = (P(k,l+stepy)*(NNZ-abs(ikik))     + P(k+step,l+stepy) *abs(ikik))   /NNZ   ! Point load in x one y step away
                        Psmall3 = (Psmall1     *(NNZ-abs(jljl))     + Psmall2           *abs(jljl))   /NNZ   ! Point load in y given x-values
                        Psmall4 = Psmall3-Psmall                                                             ! The pointloas-line load 
                        Psmall4 = Psmall4/(2.0*pi)
                        
                        iki = -ikik + NZ_cvot * ik
                        jlj = -jljl + NZ_cvot * jl
                            
                        if( iki .EQ. 0) then
                            if(jlj .EQ. 0) then          ! If we're at the nod where the pressure is applied-------------------------------------------------
                                sig_xs(i,j) = sig_xs(i,j)   -P(k,l)*0.5*(1.0+2.0*0.3) !=-P*(2*ny+(2/pi)(1-2*ny)*atan(a/b))=-P*(1+2ny)/2 since a=b =DX=DY
                                sig_ys(i,j) = sig_ys(i,j)   -P(k,l)*0.5*(1.0+2.0*0.3)
                                sig_zs(i,j) = sig_zs(i,j)   -P(k,l) 
                            else                        !If we're at the line where the pressure is applied but not at the node
                                r2          = jlj**2    
                                Psmall4     = Psmall3/(2.0*pi)*(1.0-2.0*0.3)/r2
                    
                                sig_xs(i,j) = sig_xs(i,j)   - Psmall4 !* (-jl**2)/r2  ! IK=0  
                                sig_ys(i,j) = sig_ys(i,j)   + Psmall4 !* (jl**2)/r2   ! IK=0
                            endif
                            
                        elseif(Psmall4 .GT. 0.0005) THen
                                r2=iki**2+jlj**2
                                r4=r2*r2
                                Psmall4 = Psmall4*(1.0-2.0*0.3)/r4
                    
                                sig_xs(i,j)=sig_xs(i,j)   + Psmall4 * (iki**2-jlj**2)   !Subtract the lineload since a line load not yields stresses in the perpendiculat direction on the surface. 
                                sig_ys(i,j)=sig_ys(i,j)   + Psmall4 * (jlj**2-iki**2)
                                tau_xys(i,j)=tau_xys(i,j) + Psmall4 * (2.0*iki*jlj)
                                
                        endif
                        
                        if( abs(sig_xs(i,j)) .GT. 10) then
                            sig_xs(1,1)=sig_xs(i,j)
                        endif  
                    enddo
                enddo
                
            else
                iki=ik
                jlj=jl
                Psmall4 = (P(k,l)-Pmin(k))/(2.0*pi)
                
                                r2=iki**2+jlj**2
                                r4=r2*r2
                                Psmall4 = Psmall4*(1.0-2.0*0.3)/r4
                    
                                sig_xs(i,j)=sig_xs(i,j)   + Psmall4 * (iki**2-jlj**2)   !Subtract the lineload since a line load not yields stresses in the perpendiculat direction on the surface. 
                                sig_ys(i,j)=sig_ys(i,j)   + Psmall4 * (jlj**2-iki**2)
                                tau_xys(i,j)=tau_xys(i,j) + Psmall4 * (2.0*iki*jlj)
 
            endif
            if( abs(sig_xs(i,j)) .GT. 10) then
                sig_xs(1,1)=sig_xs(i,j)
            endif  
            
            if( isnan(sig_xs(i,j))) then
                sig_xs(1,1)=sig_xs(i,j)
            endif
        enddo
 
    
    return
    
end
    
        
