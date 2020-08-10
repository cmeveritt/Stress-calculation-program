!Subroutine to fint the stresses at the surface due to a shear load
    SUBROUTINE vonMises_surf
        implicit none
        integer      NN, Nz
        integer      i, j, k, l
        real*8         ik, jl, ikk, jll
        real*8      min, pi, r2
        save        /Stress_surf/
        include     'inc_Grid.h'
        include     'inc_Stress_surf.h'
            
    pi=3.14159265
    NN=(NY+1)/2
    
 
    DO i=2,NX-1
        DO j=2,NN            ! Interested in the whole area
            
            
        ! Effective stress
        sig_vms(i,j)=sqrt(sig_xs(i,j)**2  +sig_ys(i,j)**2  +sig_zs(i,j)**2  -sig_xs(i,j)*sig_ys(i,j)  -sig_xs(i,j)*sig_zs(i,j)  -sig_ys(i,j)*sig_zs(i,j) +3*tau_xys(i,j)**2 +3*tau_xzs(i,j)**2 +3*tau_yzs(i,j)**2) 
 
        if( isnan(sig_vms(i,j))) then
                sig_xs(1,1)=sig_xs(i,j)
                sig_ys(1,1)=sig_ys(i,j)
                sig_zs(1,1)=sig_zs(i,j)
                tau_xys(1,1)=tau_xys(i,j)
                tau_xzs(1,1)=tau_xzs(i,j)
                tau_yzs(1,1)=tau_yzs(i,j)
                WRITE(4,*)'Error in Von Mises surface stress calc at i,j=',i,j
                WRITE(*,*)'Error in Von Mises surface stress calc at i,j=',i,j
        endif
            
            
        enddo
    enddo
    
        
    return
    
end
    
        
