! testprogram for calculating the eigenvalues of a symmetric real 3x3 matrix.
Subroutine Eigenvalunder
    implicit none 
    integer     I(3,3), ii ,jj, ix, jx 
    integer      NN
    real*8        A(3,3), B(3,3)
    real*8        aa(3)
    real*8        p1, eig1, eig2, eig3, q, p2, phi, r, pi, P
    real*8      tau_m
    real*8        acos, trace
    save        /Prin_Stress_under/
    save        /tau_max_values/
    include     'inc_grid.h'
    include     'inc_Stress_under.h'
    include     'inc_Prin_Stress_surf.h'
    include     'inc_Prin_Stress_under.h'
    include     'inc_StressParam.h'
    include     'inc_tau_max_values.h'
    
    NN=(NY+1)/2
    pi=3.14159265
    ! Identity matrix
    I(1,1)= 1
    I(1,2)= 0
    I(1,3)= 0   
    I(2,1)= 0
    I(2,2)= 1
    I(2,3)= 0
    I(3,1)= 0
    I(3,2)= 0
    I(3,3)= 1 
    
    
    do jx=1,Nz
        do ix=1,NX
        

        
        ! Stress matrix
        A(1,1)= sig_x(ix,jx)
        A(1,2)= 0               ! Since at the symmetry line
        A(1,3)= tau_xz(ix,jx) 
        A(2,1)= 0
        A(2,2)= sig_y(ix,jx)
        A(2,3)= tau_yz(ix,jx)  
        A(3,1)= tau_xz(ix,jx) 
        A(3,2)= tau_yz(ix,jx) 
        A(3,3)= sig_z(ix,jx)

        !! from https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3.C3.973_matrices
        p1 = A(1,2)**2 + A(1,3)**2 + A(2,3)**2
        if( p1 .eq. 0) then
            ! A is diagonal.
            eig1 = A(1,1)
            eig2 = A(2,2)
            eig3 = A(3,3)
            
            sig_1(ix,jx)=max(eig1, eig2, eig3)
            sig_3(ix,jx)=min(eig1, eig2, eig3)
            
            sig_2(ix,jx)=eig1
            If( sig_2(ix,jx) .EQ. sig_1(ix,jx) .or. sig_2(ix,jx) .EQ. sig_3(ix,jx)) sig_2(ix,jx)=eig2
            If( sig_2(ix,jx) .EQ. sig_1(ix,jx) .or. sig_2(ix,jx) .EQ. sig_3(ix,jx)) sig_2(ix,jx)=eig3
        else
            !q = trace(A)/3
            q=0
            do ii=1,3
                jj=ii    
                q=q+A(ii,jj)
            enddo
            q=q/3
        
            p2 = (A(1,1) - q)**2 + (A(2,2) - q)**2 + (A(3,3) - q)**2 + 2 * p1
            p = sqrt(p2 / 6)
            B = (1 / p) * (A - q * I)       ! I is the identity matrix
            !r = det(B) / 2
            r=  ( B(1,1)*(B(2,2)*B(3,3) - B(3,2)*B(2,3)) + B(1,2)*(B(3,1)*B(2,3) - B(2,1)*B(3,3))  + B(1,3)*(B(2,1)*B(3,2) - B(3,1)*B(2,2)))/2

            ! In exact arithmetic for a symmetric matrix  -1 <= r <= 1
            ! but computation error can leave it slightly outside this range.
            if (r .LE. -1) then
                phi = pi / 3
            else if (r .gt. 1) then
                phi = 0
            else
                phi = acos(r) / 3
            endif

            ! the eigenvalues satisfy eig3 <= eig2 <= eig1
            eig1 = q + 2 * p * cos(phi)
            eig3 = q + 2 * p * cos(phi + (2*pi/3))
            eig2 = 3 * q - eig1 - eig3    ! since trace(A) = eig1 + eig2 + eig3
            
            sig_1(ix,jx)   = eig1
            sig_2(ix,jx)   = eig2
            sig_3(ix,jx)   = eig3
            tau_m          = eig1-eig3
            tau(ix,jx)     = tau_m ! For the asperity point of view calculations
            tau_max(ix,jx) = max(tau_max(ix,jx),tau_m)  !Tau max is absolute value and since at minimum 0
        endif
        enddo
    enddo
    
    
5    RETURN
END  