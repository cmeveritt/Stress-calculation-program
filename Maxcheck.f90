Subroutine Maxcheck(t,ntime)
    !Subroutine to calculate the maximum stresses the boddy experiences
        implicit none
        integer     NN, ntime, i ,j, aa, asp_pos, t, floor, b, aai
        real*8      s1, sx, sy, sz, txy, txz, tyz,  n_x, n_y, n_z
        save        /Max_val_asp/
        save        /Max_val_cont/
        save        /tau_max_values/
        save        /directions/
   include     'inc_Grid.h'   
   include     'inc_StressParam.h'
   include     'inc_Stress_surf.h'
   include     'inc_Stress_under.h'
   include     'inc_tau_max_values.h'
   include     'inc_Prin_Stress_surf.h'
   include     'inc_Prin_Stress_under.h'
   include     'inc_Max_val_asp.h'
   include     'inc_Max_val_cont.h'
   include     'inc_CurrentP.h'
   include     'inc_directions.h'
   
        
        NN=(NY+1)/2
        
        asp_pos=NINT(1.0+Nx-(1.0*Nx)/Ntime*t)      ! The contact moves from the +side of the asperity to the - side. 
        if (asp_pos .lt.  1) asp_pos=1
        aa=asp_pos
    
        b=asp_pos+Nx-1;
        if (b .gt. 2*Nx) b=2*Nx

        !Contact view
        if (t .EQ. 0) then
                sig_1m   = sig_1                ! Under neeth surface principal stress 1
                sig_1ms  = sig_1s               ! Surface principal stress 1
                sig_vmm  = sig_vm               ! Under neeth surface von Mises stress
                sig_vmms = sig_vms              ! Surface von Mises stress
                
                sig_xm  = sig_x
                sig_xms = sig_xs
                sig_ym  = sig_y
                sig_yms = sig_ys
                sig_zm  = sig_z
                sig_zms = sig_zs
                
        else
            do i=1,NX
                do j=1,NN
                    
                    ! storing the maximum values at each node on the surface
                    sig_1ms(i,j)  =max(sig_1ms(i,j),  sig_1s(i,j))
                    sig_vmms(i,j) =max(sig_vmms(i,j), sig_vms(i,j))

                    sig_xms(i,j)  =max(sig_xms(i,j),  sig_xs(i,j))
                    sig_yms(i,j)  =max(sig_yms(i,j),  sig_ys(i,j))
                    sig_zms(i,j)  =max(sig_zms(i,j),  sig_zs(i,j))

                    
                enddo
                do j=1,NZ
                    ! storing the maximum values at each node underneeth
                    sig_1m(i,j)  = max(sig_1m(i,j),  sig_1(i,j))
                    sig_vmm(i,j) = max(sig_vmm(i,j), sig_vm(i,j))
                    
                    sig_xm(i,j)  = max(sig_xm(i,j),  sig_x(i,j))
                    sig_ym(i,j)  = max(sig_ym(i,j),  sig_y(i,j))
                    sig_zm(i,j)  = max(sig_zm(i,j),  sig_z(i,j))
                enddo
            enddo
        endif
            
        !asperity wiev
        if (t .EQ. 0) then
            
            !initiating P and sig_vm as zero becouse they can not have negative values
            P_asp=0
            sig_vmma=0
            sig_vmmsa=0 
            tau_max_a=0
            tau_max_as=0
            
            ! Sigma_1 can have negative values
            do i=1,Nx
                do j=1,NN
                    sig_1msa(aa+i-1,j) = sig_1s(i,j)
                        s1              = sig_1s(i,j)
                        sx              = sig_xs(i,j)
                        sy              = sig_ys(i,j)
                        sz              = sig_zs(i,j)
                        txy             = tau_xys(i,j)
                        txz             = tau_xzs(i,j)
                        tyz             = tau_yzs(i,j)
                        
                        Call s1_direction(s1, sx, sy, sz, txy, txz, tyz,  n_x, n_y, n_z)

                        nx_s(aa+i-1,j)      =N_x
                        ny_s(aa+i-1,j)      =N_y
                        nz_s(aa+i-1,j)      =N_z
                    sig_xmsa(aa+i-1,j) = sig_xs(i,j)
                    sig_ymsa(aa+i-1,j) = sig_ys(i,j)
                    sig_zmsa(aa+i-1,j) = sig_zs(i,j)
                    
                enddo
                do j=1, NZ
                        

                        s1              = sig_1(i,j)
                        sx              = sig_x(i,j)
                        sy              = sig_y(i,j)
                        sz              = sig_z(i,j)
                        txy             = tau_xy(i,j)
                        txz             = tau_xz(i,j)
                        tyz             = tau_yz(i,j)
                        
                        Call s1_direction(s1, sx, sy, sz, txy, txz, tyz,  n_x, n_y, n_z)

                        nx_d(aa+i-1,j)     =N_x
                        ny_d(aa+i-1,j)     =N_y
                        nz_d(aa+i-1,j)     =N_z
                        
                    sig_1ma(aa+i-1,j)  = sig_1(i,j)
                    
                    sig_xma(aa+i-1,j)  = sig_x(i,j)
                    sig_yma(aa+i-1,j)  = sig_y(i,j)
                    sig_zma(aa+i-1,j)  = sig_z(i,j)
                enddo
            enddo
                    
        else
            do i=1,Nx
                aai=aa+i-1
                if (aai .gt. 2*Nx) aai=2*Nx             !Safety check
                
                ! surface
                do j=1,NN
                    P_asp(aai,j)     = max(P_asp      (aai,j), P(i,j))
                    if ( sig_1s(i,j) .GT. sig_1msa(aai,j))  then
                        sig_1msa(aai,j) = sig_1s(i,j)
                        
                        
                        s1              = sig_1s(i,j)
                        sx              = sig_xs(i,j)
                        sy              = sig_ys(i,j)
                        sz              = sig_zs(i,j)
                        txy             = tau_xys(i,j)
                        txz             = tau_xzs(i,j)
                        tyz             = tau_yzs(i,j)
                        
                        Call s1_direction(s1, sx, sy, sz, txy, txz, tyz,  n_x, n_y, n_z)

                        nx_s(aai,j)     =N_x
                        ny_s(aai,j)     =N_y
                        nz_s(aai,j)     =N_z
                    endif
                
                        
                    sig_1msa(aai,j)  = max(sig_1msa   (aai,j), sig_1s(i,j))
                    sig_vmmsa(aai,j) = max(sig_vmmsa  (aai,j), sig_vms(i,j))
                    tau_max_as(aai,j)= max(tau_max_as (aai,j), taus(i,j))
                    
                    sig_xmsa(aai,j) = max(sig_xmsa   (aai,j), sig_xs(i,j))
                    sig_ymsa(aai,j) = max(sig_ymsa   (aai,j), sig_ys(i,j))
                    sig_zmsa(aai,j) = max(sig_zmsa   (aai,j), sig_zs(i,j))
   
                enddo
                !depth
                do j=1, NZ
                    
                    if ( sig_1(i,j) .GT. sig_1ma(aai,j))  then
                        sig_1ma(aai,j) = sig_1(i,j)

                        s1              = sig_1(i,j)
                        sx              = sig_x(i,j)
                        sy              = sig_y(i,j)
                        sz              = sig_z(i,j)
                        txy             = tau_xy(i,j)
                        txz             = tau_xz(i,j)
                        tyz             = tau_yz(i,j)
                        
                        Call s1_direction(s1, sx, sy, sz, txy, txz, tyz,  n_x, n_y, n_z)

                        nx_d(aai,j)     =N_x
                        ny_d(aai,j)     =N_y
                        nz_d(aai,j)     =N_z
                    endif
                    
                    sig_1ma(aai,j)  = max(sig_1ma    (aai,j), sig_1(i,j))
                    sig_vmma(aai,j) = max(sig_vmma   (aai,j), sig_vm(i,j))
                    tau_max_a(aai,j)= max(tau_max_a  (aai,j), tau(i,j))
                    
                    sig_xma(aai,j)  = max(sig_xma    (aai,j), sig_x(i,j))
                    sig_yma(aai,j)  = max(sig_yma    (aai,j), sig_y(i,j))
                    sig_zma(aai,j)  = max(sig_zma    (aai,j), sig_z(i,j))
                enddo
                
            enddo
        endif
            
        RETURN
    end
    

    ! A subroutine to calculate the direction of s1 in the x,y,z coordinate system. 
    Subroutine s1_direction(s1, sx, sy, sz, txy, txz, tyz,  nx, ny, nz)
    implicit none
    real*8          s1, sx, sy, sz, txy, txz, tyz,  nx, ny, nz
    real*8          nxn, nyn, nzn
    real*8          Cx1, Cx2, Cy1, Cy2, Cz1, Cz2

        if( s1 .EQ. sx) then
            nx=1
            ny=0
            nz=0
        elseif(s1 .EQ. sy)then
            nx=0
            ny=1
            nz=0
        elseif(s1 .EQ. sz) then
            nx=0
            ny=0
            nz=1
        else! From page 5 in Formelsamlingen. If one direction component is close to zero, the calculations gets more messy. Becouse of this we have the if statements

             
            Cx1              =  (txz*tyz/((s1-sz)*(s1-sy)) + txy/(s1-sy)) /( 1.0 - tyz*tyz/(s1-sz)*(s1-sy))
            Cy1              =  (txy*txz/((s1-sx)*(s1-sz)) + tyz/(s1-sz)) /( 1.0 - txz*txz/(s1-sx)*(s1-sz))
            Cz1              =  (tyz*txy/((s1-sy)*(s1-sx)) + txz/(s1-sx)) /( 1.0 - txy*txy/(s1-sy)*(s1-sx))

            Cx2              =  (txz+Cx1*tyz)/(s1-sz)
            Cy2              =  (txy+Cy1*txz)/(s1-sx)
            Cz2              =  (tyz+Cz1*txy)/(s1-sy)
            
            nx               = 1.0/sqrt(1.0**2 +Cx1**2 +Cx2**2)
            ny               = 1.0/sqrt(1.0**2 +Cy1**2 +Cy2**2)
            nz               = 1.0/sqrt(1.0**2 +Cz1**2 +Cz2**2)
            
            ! If rotating around on axis, that n becomes 1 from the above equations, but that should not be the case. 
            ! Not needed when I take the middle value
            !If( Cx1 .eq. 0 .and. Cx2 .EQ. 0) nx=0
            !If( Cy1 .eq. 0 .and. Cy2 .EQ. 0) ny=0
            !If( Cz1 .eq. 0 .and. Cz2 .EQ. 0) nz=0
            
            ! Try with taking the middle value instead since it did not work with the greatest
            If((abs(nx) .GT. abs(ny) .AND. abs(NX) .LE. abs(nz)) .or. (abs(nx) .LT. abs(ny) .AND. abs(NX) .GT. abs(nz)) )then 
                ny          = Cx1*nx
                nz          = Cx2*nx
            elseif( (abs(ny) .GT. abs(nx) .AND. abs(Ny) .LE. abs(nz)) .or. (abs(ny) .LT. abs(nx) .AND. abs(Ny) .GT. abs(nz)) ) then
                nz          = Cy1*ny
                nx          = Cy2*ny
            else
                nx          = Cz1*nz
                ny          = Cz2*nz
            endif
            
        endif
                        
        return
    end
    