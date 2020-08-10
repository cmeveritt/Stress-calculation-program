!A subroutine for smothening the depth data near the surface
    subroutine smothening
        implicit none
        integer    i
        save    /Stress_under/
        include     'inc_Grid.h'
        include     'inc_Stress_under.h'

    
    !Near the surface the stresses in the z-direction obtains  not exactly great walues
    DO i=2,NX-1
       sig_x(i,2)=0.5*(sig_x(i,1)+sig_x(i,3))
       sig_y(i,2)=0.5*(sig_y(i,1)+sig_y(i,3))
       sig_z(i,2)=0.5*(sig_z(i,1)+sig_z(i,3))
       tau_xz(i,2)=0.5*(tau_xz(i,1)+tau_xz(i,3))
       
       !tau_yz = tau_xy=0
    enddo
    
    return
    end
    