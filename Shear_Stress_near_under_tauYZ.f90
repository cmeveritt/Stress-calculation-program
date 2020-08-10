! Subrutine to calculate the stresses  in the x-z plane due to shear loads at the centerline of the contact
    ! Not fully debuged since not nessesery at symmetry line
    SUBROUTINE Shear_Stress_near_under_tauYZ
        implicit none
 !renaming the surface stresses so they do no collide with the undersurface stresses
        integer       NN, NNX
        integer      i, j, k, l
        integer      ik, jl, js, jll
        integer      ikik, jljl, step, stepy, iki, jlj
        real*8         Psmall1, Psmall2, Psmall3, Psmall4, Psmall, z20
        real*8         sigxx, sigzz
        real*8         rho5, rho3, base
        real*8         rhoz2, rho3z2, rho2z3,ikz
        real*8         min, pi
        real*8         rho, r2, z
        integer      hh, NNZ
        save        /Stress_under/
        include     'inc_Grid.h'
        include     'inc_CurrentP.h'
        include     'inc_StressParam.h'
        include     'inc_Stress_surf.h'
        include     'inc_Stress_under.h'
        include     'inc_Minimum_values.h'
        include     'inc_Nodes_no_pressure.h'
        
    pi=3.14159265
    NN=(NY+1)/2
    NNX=(NX+1)/2
    NNZ=(NZ_cvot-1)/2       
    
    ! Calculats the stresses from the pressure at the area         
         ! Line load -----------------------------------------------------------------------------------------------------
    DO hh=2,Nr_z*NZ_cvot+1                         ! For each height below the surface
        z=(hh-1.0)/NZ_cvot                     ! To have same scale on z as ik and jl 
        DO i=NXstart,NXstop
            j=NN                        ! Only interested in the centerline
      
            ! At the surface a line load, Pmin, will not couse any stresses exept than where its applied. Point loads will on the other hand.  
            Do k=NXstart,NXstop
                ik=(k-i)
                
                if (abs(ik) .LT. 5) then ! If we're close to the line where the pressure is applied
                        
                    DO ikik=-NNZ,NNZ
                        if (ikik .GT. 0)then 
                            step=1 
                        elseif( ikik .EQ. 0) then 
                            step=0
                        else 
                            step=-1
                        endif
                    
                        Psmall=(tau_yzl_min(k)*(NNZ-abs(ikik))+tau_yzl_min(k+step)*abs(ikik))/NNZ
                        z20=(hh-1.0)
                        iki=ikik+NZ_cvot*ik
                        
                        !sig_x(i,hh)=sig_x(i,hh)    -2*tau_xzs(k,l)*ik**3/(pi*base)
                        !sig_z(i,hh)=sig_z(i,hh)    +0.3/(1-2*0.3)*(-2*tau_xzs(k,l)*ik**3/(pi*base)-2*tau_xzs(k,l)*ik*z**2/(pi*base))  
                        !sig_y(i,j)=sig_y(i,j)      -2*tau_xzs(k,l)*ik*z**2/(pi*base)
                
                        !tau_xz(i,hh)=tau_xz(i,hh)  -2*tau_xzs(k,l)ik**2*z/(pi*base)
                        ikz =   (iki**2+z20**2)
                 
                        ! No shear xy stresses due to symmetry
                        !tau_xy(i,hh)=tau_xy(i,hh)   -tau_yzl_min(k)*ik/ (pi*ikz)
                        tau_yz(i,hh)=tau_yz(i,hh)   -Psmall*z20 / (pi*ikz) 

                    enddo
     
                else
                    
                !sig_x(i,hh)=sig_x(i,hh)    -2*tau_xzs(k,l)*ik**3/(pi*base)
                !sig_z(i,hh)=sig_z(i,hh)    +0.3/(1-2*0.3)*(-2*tau_xzs(k,l)*ik**3/(pi*base)-2*tau_xzs(k,l)*ik*z**2/(pi*base))  
                !sig_y(i,j)=sig_y(i,j)      -2*tau_xzs(k,l)*ik*z**2/(pi*base)
                
                !tau_xz(i,hh)=tau_xz(i,hh)  -2*tau_xzs(k,l)ik**2*z/(pi*base)
                 ikz =   (ik**2+z**2)
                 
                ! No shear xy stresses due to symmetry
                !tau_xy(i,hh)=tau_xy(i,hh)   -tau_yzl_min(k)*ik/ (pi*ikz)
                tau_yz(i,hh)=tau_yz(i,hh)   -tau_yzl_min(k)*z / (pi*ikz)

                endif
            enddo

        enddo
    enddo
                
    
 ! Combined with the stresses from pointloads----------------------------------------------------------------------     
            
    DO hh=2,Nr_z*NZ_cvot+1                        ! For each height below the surface
        z=(hh-1)/NZ_cvot                   ! To have same scale on z as ik and jl 

        DO i=NXstart,NXstop
            j=NN     
            DO jll=-NNX,NNX
            DO K=NXstart,NXstop

                l=jll+j
  
                IF( l .GT. NY)    l=2    ! If outside the area
                IF( l .LT. 2)     l=2
                
                ik=(k-i)
                jl=jll 
                    
                IF (abs(ik) .LT. 5 .AND. abs(jl) .LT. 5) then ! If we're close to the line where the pressure is applied
                    
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
                        
                        Psmall  = (tau_yzl_min(k)     *(NNZ-abs(ikik))     + tau_yzl_min(k+step)     *abs(ikik))   /NNZ   ! Line load
                        Psmall1 = (tau_yzs(k,l)       *(NNZ-abs(ikik))     + tau_yzs(k+step,l)       *abs(ikik))   /NNZ   ! Point loas in x
                        Psmall2 = (tau_yzs(k,l+stepy) *(NNZ-abs(ikik))     + tau_yzs(k+step,l+stepy) *abs(ikik))   /NNZ   ! Point load in x one y step away
                        Psmall3 = (Psmall1            *(NNZ-abs(jljl))     + Psmall2                 *abs(jljl))   /NNZ   ! Point load in y given x-values
                        Psmall4 = Psmall3-Psmall                                                                ! The pointloas-line load 
                        
                        IF( Psmall4 .GT. 0.0001)THEN            ! To gain some speed
                            z20=(hh-1.0)
                            iki=ikik+NZ_cvot*ik
                            jlj=jljl+NZ_cvot*jl
                
               
                            r2=iki**2+jlj**2                  ! The square of the  horizontal distance from the point [m**2]
                            if(r2 .LT. 7000) THEN           ! Else the contribution will be really small. 
                                IF(r2 .EQ. 0.0) r2=0.5**2       ! If straight underneath the pressure, assume we're half a distance away. 
                                rho=(iki**2+jlj**2+z20**2)**0.5     ! The total distance [m]
                        
                                rho5=rho**5
                                rho3=rho**3
                                rhoz2=(rho*(rho+z20))**2
                                rho3z2=rho**3*(rho+z20)**2
                                rho2z3=rho**2*(rho+z20)**3
                                

                                ! The coordinate rotation from that the shear traction is applied in y instead of x direction gives:
                                ! sig_x=sig_y, sig_y=sig_x, sig_z=sig_z,  tau_xy=-tau_xy, tau_xz=tau_yz, tau_yz=tau_xz
                        
                                sig_y(i,hh)=sig_y(i,hh)   +(Psmall4)/(2*pi)* ( 3*iki**3/rho5    + (1-2*0.3) * (iki/rho3 - 3*iki/rhoz2 + iki**3   /rho3z2 + 2*iki**3   /rho2z3))
                                sig_x(i,hh)=sig_x(i,hh)   +(Psmall4)/(2*pi)* ( 3*iki*jlj**2/rho5 + (1-2*0.3) * (iki/rho3 -   iki/rhoz2 + iki*jlj**2/rho3z2 + 2*iki*jlj**2/rho2z3))
                        
                                sig_z(i,hh)=sig_z(i,hh)   -(Psmall4)/(2*pi)*  3*iki*z20**2/rho5
                        
                                ! No shear xy stresses due to symmetry
                                !tau_xy(i,hh)=tau_xy(i,hh) -(tau_yzs(k,l)-tau_yzl_min(k))/(2*pi)* ( 3*ik**2*jl/rho5 + (1-2*0.3) * (        -   jl/rhoz2 + ik**2*jl/rho3z2 + 2*ik**2*jl/rho2z3))
                                tau_yz(i,hh)=tau_yz(i,hh) -(Psmall4)/(2*pi)*  3*iki*jlj*z20/rho5
                                tau_xz(i,hh)=tau_xz(i,hh) -(Psmall4)/(2*pi)*  3*iki**2*z20/rho5
                

                        

                            ENDIF
                        ENDIF
                    ENDDO
                    ENDDO
                      
                ELSE
                    r2=ik**2+jl**2                  ! The square of the  horizontal distance from the point [m**2]
                    if(r2 .LT. 7000) THEN           ! Else the contribution will be really small. 
                    IF(r2 .EQ. 0.0) r2=0.5**2       ! If straight underneath the pressure, assume we're half a distance away. 
                        rho=(ik**2+jl**2+z**2)**0.5     ! The total distance [m]
                        
                        rho5=rho**5
                        rho3=rho**3
                        rhoz2=(rho*(rho+z))**2
                        rho3z2=rho**3*(rho+z)**2
                        rho2z3=rho**2*(rho+z)**3
                        ! The coordinate rotation from that the shear traction is applied in y instead of x direction gives:
                        ! sig_x=sig_y, sig_y=sig_x, sig_z=sig_z,  tau_xy=-tau_xy, tau_xz=tau_yz, tau_yz=tau_xz
                        
                        sig_y(i,hh)=sig_y(i,hh)   +(tau_yzs(k,l)-tau_yzl_min(k))/(2*pi)* ( 3*ik**3/rho5    + (1-2*0.3) * (ik/rho3 - 3*ik/rhoz2 + ik**3   /rho3z2 + 2*ik**3   /rho2z3))
                        sig_x(i,hh)=sig_x(i,hh)   +(tau_yzs(k,l)-tau_yzl_min(k))/(2*pi)* ( 3*ik*jl**2/rho5 + (1-2*0.3) * (ik/rho3 -   ik/rhoz2 + ik*jl**2/rho3z2 + 2*ik*jl**2/rho2z3))
                        
                        sig_z(i,hh)=sig_z(i,hh)   -(tau_yzs(k,l)-tau_yzl_min(k))/(2*pi)*  3*ik*z**2/rho5
                        
                        ! No shear xy stresses due to symmetry
                        !tau_xy(i,hh)=tau_xy(i,hh) -(tau_yzs(k,l)-tau_yzl_min(k))/(2*pi)* ( 3*ik**2*jl/rho5 + (1-2*0.3) * (        -   jl/rhoz2 + ik**2*jl/rho3z2 + 2*ik**2*jl/rho2z3))
                        tau_yz(i,hh)=tau_yz(i,hh) -(tau_yzs(k,l)-tau_yzl_min(k))/(2*pi)*  3*ik*jl*z/rho5
                        tau_xz(i,hh)=tau_xz(i,hh) -(tau_yzs(k,l)-tau_yzl_min(k))/(2*pi)*  3*ik**2*z/rho5
                

                        

                    ENDIF
                ENDIF
            ENDDO 
            ENDDO
            
 
        ENDDO
    ENDDO
                    
     
    
    RETURN
    
end
    