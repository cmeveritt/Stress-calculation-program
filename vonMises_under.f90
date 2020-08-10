! Subrutine to calculate the stresses  in the x-z plane due to shear loads at the centerline of the contact
    SUBROUTINE vonMises_under
        implicit none
        integer      NN, i, j
        integer      hh
        save        /Stress_under/
        include     'inc_Grid.h'
        include     'inc_StressParam.h'
        include     'inc_Stress_under.h'
        
    NN=(NY+1)/2
    j=NN  
    
    ! Total depth investigated is the half contact halfwidth b=1 in dimless scale
    DO hh=2,Nz       ! For each height below the surface
        
    
        ! Calculats the stresses from the pressure at the area
        DO i=2,NX-1
                       ! Only interested in the centerline

            ! Effective stress
           sig_vm(i,hh)=sqrt(sig_x(i,hh)**2  +sig_y(i,hh)**2  +sig_z(i,hh)**2  -sig_x(i,hh)*sig_y(i,hh)  -sig_x(i,hh)*sig_z(i,hh)  -sig_y(i,hh)*sig_z(i,hh)  +3*tau_xz(i,hh)**2  +3*tau_yz(i,hh)**2) !tau_xy has to be 0
            
            if( isnan(sig_vm(i,j))) then
                sig_x(1,1)=sig_x(i,j)
                sig_y(1,1)=sig_y(i,j)
                sig_z(1,1)=sig_z(i,j)
                tau_xy(1,1)=tau_xy(i,j)
                tau_yz(1,1)=tau_yz(i,j)
                WRITE(4,*)'Error in Von Mises stress calc at i,j=',i,j
                WRITE(*,*)'Error in Von Mises stress calc at i,j=',i,j
            endif
        enddo
    enddo

    
    RETURN
    
end
    