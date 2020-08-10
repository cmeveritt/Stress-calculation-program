! Subrutine to calculate the stresses  in the x-z plane at the centerline of the contact
    SUBROUTINE StressP_near_Under_point(hh,i,j,jl)
        implicit none
        integer      NN,  NNX
        integer      i, j, k, l
        integer      ik, jl, js
        integer      ikik, jljl, step, stepy, iki, jlj
        real*8         Psmall1, Psmall2, Psmall3, Psmall4
        real*8         sigxx, sigzz
        real*8         rho3, rho5, base, Psmall, z20
        real*8         min, pi
        real*8         rho, r2, z, r4
        real*8        term1,term2, term3, term4, term5, term6, term7
        real*8        t1, t2, t3, t4, t5, t6, t7, t8
        real*8        al, be, ga, de
        real*8        Psi1_xx, Psi_xx, Psi_z, Psi1_x
        real*8        tA, tB, tC, tD, Te1, Te2, Te3, Te4
        real*8        s1, s2, s3, s4, s5
        integer      hh, NNZ
        save        /Stress_under/
        include     'inc_Grid.h'
        include     'inc_Nodes_no_pressure.h'
        include     'inc_Stress_surf.h'
        include     'inc_CurrentP.h'
        include     'inc_Minimum_values.h'
        include     'inc_StressParam.h'
        include     'inc_Stress_under.h'
        
    pi=3.14159265
    NN=(NY+1)/2
    NNX=(NX+1)/2
    NNZ=(NZ_cvot-1)/2       
    IF( hh .LT. Nr_z*NZ_cvot+2)Then
        z=(hh-1)/NZ_cvot                   ! To have same scale on z as ik and jl  
    else
        z=(hh-Nr_z*(NZ_cvot-1)-1) 
    endif
    
            DO K=NXstart,NXstop

                l=j-jl
  
                IF( l .GE. NY)    l=2    ! If outside the area
                IF( l .LT. 2)     l=2
                
                ik=i-k
                    
                IF (abs(ik) .LT. Nr_xy .AND. abs(jl) .LT. Nr_xy ) then ! If we're close to the line where the pressure is applied
                    
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
                        Psmall4 = Psmall3-Psmall                                                                ! The pointloas-line load 
                        Psmall4 = Psmall4/(2.0*pi)


                        z20=z*NZ_cvot
                        iki=-ikik + NZ_cvot*ik
                        jlj=-jljl + NZ_cvot*jl
                            
                        IF( iki .EQ. 0 .AND. z20 .EQ. 0) THEN
                            If(jlj .EQ. 0) THEN 
                                Psmall4 = Psmall3       
                            else
                                Psmall4 = Psmall3/(2.0*pi)
                            endif
                        endif
                        
                        IF( Psmall4 .GT. 0.0005)THEN            ! To gain some speed
                            Call StressP_near_Under_point_stress(i,hh, iki, jlj, z20, Psmall4)
               
                            if(isnan(sig_x(i,hh)) .or. isnan(sig_z(i,hh)) ) then
                                Psmall=P(i,j)
                            endif  
                        ENDIF
                        
                    ENDDO
                    ENDDO
                      
                ELSE
                    Psmall4 = (P(k,l)-Pmin(k))/(2.0*pi)
                    IF( Psmall4 .GT. 0.0005) Call StressP_near_Under_point_stress(i,hh, ik, jl, z, Psmall4)
                    
                    if(isnan(sig_x(i,hh)) .or. isnan(sig_z(i,hh)) ) then
                        Psmall=P(i,j)
                    endif  
                ENDIF
            ENDDO 
           
    RETURN
    
end
    