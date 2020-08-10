Subroutine Surface_at_top_depth
        integer     i, NN
        save    /Stress_under/
        include     'inc_Grid.h'
        include     'inc_Stress_surf.h'
        include     'inc_Stress_under.h'

        
        
    NN=(NY+1)/2
    
    !Near the surface the stresses in the z-direction obtains  not exactly great walues
    DO i=1,NX
       sig_x(i,1)  = sig_xs(i,NN)
       sig_y(i,1)  = sig_ys(i,NN)
       sig_z(i,1)  = sig_zs(i,NN)
       tau_xz(i,1) = tau_xzs(i,NN)
       
       !tau_yz = tau_xy=0
    enddo
    
    return
    end
    