Subroutine Findley_print
    !Soubroutine for printing the maximum values when calculated
    implicit none
       integer      NN, i
       real         t, j
    include     'inc_grid.h'
    include     'inc_Stress_surf.h'
    include     'inc_Stress_under.h'
    include     'inc_StressParam.h'
    include     'inc_Coord.h'
        
        NN=(NY+1)/2


        ! In EHL for Findley 1-3 is whatever, then fist coordinate, second coordinate used as time, sxx, syy, szz, txy,tyz,txz
        
        t=1 !Setting t as real to se if prints better
        ! Contact viev---------------------------------------------------
        Do j=1,NN
            DO i=1,Nx
                WRITE(101,110)0,0,0,X(i),j,sig_xs(i,j),sig_ys(i,j),sig_zs(i,j), tau_xys(i,j),tau_yzs(i,j),tau_xzs(i,j)     !Surface
            ENDDO 
        enddo
        
        t=1 !Setting t as real to se if prints better
        ! Contact viev---------------------------------------------------
        Do j=1,Nz
            DO i=1,Nx
                WRITE(102,110)0,0,0,X(i),j,sig_x(i,j),sig_y(i,j),sig_z(i,j), 0.0, tau_yz(i,j),tau_xz(i,j)     !Surface Tau_xy is zero at symmetry line
            ENDDO 
        enddo

        
    110 FORMAT(2001(E12.6,1X))
        RETURN
    END