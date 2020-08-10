! Subrutine to calculate the stresses  in the x-z plane at the centerline of the contact
    SUBROUTINE StressP_near_Under_point_stress(i,hh, iki, jlj, z20, Psmall4)
        implicit none
        integer      NN,  NNX
        integer      i, j, k, l
        integer      ik, jl, js, jll 
        integer      ikik, jljl, step, stepy, iki, jlj
        real*8         Psmall1, Psmall2, Psmall3, Psmall4
        real*8         sigxx, sigzz
        real*8         rho3, rho5, base, Psmall, z20
        real*8         min, pi
        real*8         rho, r2, r4, z
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
        include     'inc_load_surf.h'
        
    pi=3.14159265
    NN=(NY+1)/2
    NNX=(NX+1)/2
    NNZ=(NZ_cvot-1)/2 
    
    IF( iki .EQ. 0 .AND. z20 .EQ. 0) THEN
                    If(jlj .EQ. 0) THEN          ! If we're at the nod where the pressure is applied-------------------------------------------------
                        sig_x(i,hh) = sig_x(i,hh)   -Psmall4*0.5*(1.0+2.0*0.3) !=-P*(2*ny+(2/pi)(1-2*ny)*atan(a/b))=-P*(1+2ny)/2 since a=b =DX=DY
                        sig_y(i,hh) = sig_y(i,hh)   -Psmall4*0.5*(1.0+2.0*0.3)
                        sig_z(i,hh) = sig_z(i,hh)   -Psmall4     ! or , P(i,NN), it's the same            
                    else ! We're at the line where the pressure is applied. The line load here is not accounted for.
                        r2=iki**2+jlj**2
                        sig_x(i,hh) = sig_x(i,hh)   - Psmall4*(1.0-2.0*0.3)/r2 !* (-jlj**2)/r2  ! IK=0  
                        sig_y(i,hh) = sig_y(i,hh)   + Psmall4*(1.0-2.0*0.3)/r2 !* ( jlj**2)/r2   ! IK=0
                    endif
    else
    
                                r2=iki**2+jlj**2                  ! The square of the  horizontal distance from the point [m**2]
                                IF(r2 .EQ. 0.0) r2=0.5**2       ! If straight underneath the pressure, assume we're half a distance away. 
                                rho=(iki**2+jlj**2+z20**2)**0.5     ! The total distance [m]
                        
                                rho5=rho**5
                                rho3=rho**3
                                
                                term1 = (1.0-2.0*0.3)/r2
                                term2 = (1.0-z20/rho)
                                term3 = (iki**2-jlj**2)/r2 
                                term4 = z20*jlj**2/rho3
                                term5 = z20*iki**2/rho3

                                
                                sig_x(i,hh)=sig_x(i,hh)   + Psmall4* (  term1 * ( term2 *   term3       + term4 )  -3.0*z20*iki**2/rho5)
                                if (sig_x(i,hh) .GT. 0.5) then ! Theses is verry high values which are obtained in the exit region
                                    sig_x(1,1)=sig_x(i,hh)
                                endif
                                
                                sig_y(i,hh)=sig_y(i,hh)   + Psmall4* (  term1 * ( term2 *   (-term3)    + term5 )  -3.0*z20*jlj**2/rho5)
                    
                                sig_z(i,hh)=sig_z(i,hh)   -3.0*(Psmall4)*z20**3/rho5
                    
                                tau_xz(i,hh)=tau_xz(i,hh) -3.0*(Psmall4)*iki*z20**2/rho5
    endif                     
    
        if(isnan(sig_x(i,hh)) .or. isnan(sig_z(i,hh)) ) then
            Psmall=P(i,1)
        endif  
  
    

    RETURN
    
    end
    
    
     ! Hills Contact mechanics eq (11.X) and (12.X)
                                !al = sqrt( (iki-0.5)**2 + (jlj -0.5)**2 + z20**2 )
                                !be = sqrt( (iki+0.5)**2 + (jlj -0.5)**2 + z20**2 )
                                !ga = sqrt( (iki+0.5)**2 + (jlj +0.5)**2 + z20**2 )
                                !de = sqrt( (iki-0.5)**2 + (jlj +0.5)**2 + z20**2 )
                                !
                                !tA = atan(( iki -0.5 + jlj - 0.5 + al)  /  (z20))
                                !tB = atan(( iki +0.5 + jlj - 0.5 + be)  /  (z20))
                                !tC = atan(( iki +0.5 + jlj + 0.5 + ga)  /  (z20))
                                !tD = atan(( iki -0.5 + jlj + 0.5 + de)  /  (z20))
                                !
                                !Psi_z = 2.0* ( tA + tC- tB- tD)
                                !
                                !t1=atan(        (0.5-jlj)   / (0.5-iki) )
                                !t2=atan(        (0.5+jlj)   / (0.5-iki) )
                                !t3=atan(z20/al* (0.5-jlj)   / (0.5-iki) )
                                !
                                !t4=atan(z20/de* (0.5+jlj)   / (0.5-iki) )
                                !t5=atan(        (0.5-jlj)   / (0.5+iki) )
                                !t6=atan(        (0.5+jlj)   / (0.5+iki) )
                                !
                                !t7=atan(z20/be* (0.5-jlj)   / (0.5+iki) )
                                !t8=atan(z20/ga* (0.5+jlj)   / (0.5+iki) )
                                !
                                !Psi1_xx= ( t1 + t2 - t3 - t4 + t5 +t6 -t7 - t8)
                                !
                                !t1     = ( (0.5 - iki) / ( ( 0.5 - iki)**2+ z20**2 ) )
                                !t2     = ( (0.5 - jlj) / al + ( 0.5 + jlj) / de ) 
                                !t3     = ( (0.5 + iki) / ( ( 0.5 + iki)**2+ z20**2 ) )
                                !t4     = ( (0.5 - jlj) / be + ( 0.5 + jlj) / ga ) 
                                !
                                !Psi_xx = ( t1 * t2 + t3*t4)
                                !
                                !Te1= atan( (z20 + jlj - 0.5 + al) / ( iki - 0.5) )
                                !Te2= atan( (z20 + jlj - 0.5 + be) / ( iki + 0.5) )
                                !Te3= atan( (z20 + jlj + 0.5 + ga) / ( iki + 0.5) )
                                !Te4= atan( (z20 + jlj + 0.5 + de) / ( iki - 0.5) )
                                !
                                !s1 = (jlj + 0.5) * log( (z20+de) / (z20+ga) ) 
                                !s2 = (jlj - 0.5) * log( (z20+be) / (z20+al) )
                                !s3 = z20*log( (jlj + 0.5 + de)*( jlj - 0.5 +be) / ( (jlj - 0.5 + al)*( jlj + 0.5 +ga) ) )
                                !s4 = 2.0 * ( iki - 0.5 ) * ( Te4 - Te1 )
                                !s5 = 2.0 * ( iki + 0.5 ) * ( Te3 - Te2 )
                                !
                                !Psi1_x = s1 + s2 + s3 +s4 - s5
                                !
                                !term6 = ( 2.0*0.3* Psi_z - z20*Psi_xx - (1.0-2.0*0.3) * Psi1_xx )
                                !term7 = (Psmall4)* (  term1 *    ( term2 *   term3       + term4 )  -3.0*z20*iki**2/rho5)