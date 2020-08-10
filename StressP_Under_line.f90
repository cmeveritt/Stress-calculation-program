! Asubroutine for evaluating the stresses down in the material due to a pure line load of normal pressure
    subroutine StressP_Under_line(hh,i,j)
        implicit none
        integer      NN, NNX
        integer      i, j, k, l
        integer      ik, jl, js, jll 
        integer      ikik, jljl, step, stepy, iki, jlj
        real*8         Psmall1, Psmall2, Psmall3, Psmall4
        real*8         sigxx, sigzz
        real*8         rho3, rho5, base, Psmall, z20
        real*8         min, pi
        real*8         rho, r2, z
        integer      hh, NNZ
        save        /Stress_under/
        include     'inc_Grid.h'
        include     'inc_Nodes_no_pressure.h'
        include     'inc_Stress_surf.h'
        include     'inc_Stress_under.h'
        include     'inc_CurrentP.h'
        include     'inc_Minimum_values.h'
        include     'inc_StressParam.h'
        
        pi=3.14159265
        NN=(NY+1)/2
        NNX=(NX+1)/2
        NNZ=(NZ_cvot-1)/2
        
        
        if( hh .LT. Nr_z*NZ_cvot+2 )then
            z=(hh-1.0)/NZ_cvot 
        else
            z=(hh-Nr_z*(NZ_cvot-1)-1)                    ! To have same scale on z as ik and jl 
        endif
        
            ! Line load -----------------------------------------------------------------------------------------------------
            ! At the surface a line load, Pmin, will not couse any stresses exept than where its applied. Point loads will on the other hand.  
            Do k=NXstart,NXstop
                ik=i-k
                
                if (abs(ik) .LT. Nr_xy) then ! If we're close to the line where the pressure is applied
                        
                    DO ikik=-NNZ,NNZ
                        if (ikik .GT. 0)then 
                            step=1 
                        elseif( ikik .EQ. 0) then 
                            step=0
                        else 
                            step=-1
                        endif
                    
                        Psmall=(Pmin(k)*(NNZ-abs(ikik))+Pmin(k+step)*abs(ikik))/NNZ
                        z20=      Nz_cvot*z
                        iki=-ikik+NZ_cvot*ik
                        IF( iki .EQ. 0 .AND. z20 .EQ. 0) THEN
                            ! These values are added later
                        else
                            
                        base=(iki**2+(z20)**2)**2
                        
                        sigxx=-2.0*Psmall/(pi)* (iki**2*z20)  /base
                        sigzz=-2.0*Psmall/(pi)* (z20**3)     /base
                
                        sig_x(i,hh)=sig_x(i,hh)    +sigxx
                        sig_z(i,hh)=sig_z(i,hh)    +sigzz
                        sig_y(i,hh)=sig_y(i,hh)    +0.3*(sigzz+sigxx)  ! FS page 22 for plane deformation
                
                        tau_xz(i,hh)=tau_xz(i,hh)  -2.0*Psmall/(pi)* (iki*z20**2)  /base
                        !tau_xy(i,hh)=tau_xy(i,hh) +0
                        !tau_yz(i,hh)=tau_yz(i,hh) +0
                        ENDIF
                        
                    enddo
                        
                else
                    
                base=(ik**2+z**2)**2    
                sigxx=-2.0*Pmin(k)/(pi)* (ik**2*z)  /base
                sigzz=-2.0*Pmin(k)/(pi)* (z**3)     /base
                
                sig_x(i,hh)=sig_x(i,hh)    +sigxx
                sig_z(i,hh)=sig_z(i,hh)    +sigzz
                sig_y(i,hh)=sig_y(i,hh)    +0.3*(sigzz+sigxx)  ! FS page 22 for plane deformation
                
                tau_xz(i,hh)=tau_xz(i,hh)  -2.0*Pmin(k)/(pi)* (ik*z**2)  /base
                !tau_xy(i,hh)=tau_xy(i,hh) +0
                !tau_yz(i,hh)=tau_yz(i,hh) +0
                endif
                
            enddo
            if(isnan(sig_x(i,hh)) .or. isnan(sig_z(i,hh)) ) then
                Psmall=P(i,j)
            endif
            

    return
     end    