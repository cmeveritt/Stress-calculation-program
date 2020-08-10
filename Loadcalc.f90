! A Subroutine for evaluating the minimum pressure and stresses at the surface. 
    subroutine Loadcalc
            implicit none
        integer     I0, I1, I2, J0, J1, J2
        real*8       Pmax(601)
        integer      NN, nstop
        integer      i, j, k, l, jj
        real*8      ik, jl, ikk, jll
        real*8      pi, r2, t2(97), t3(97), t4(97)
        real*8      h_div_b, dpdx, dpdy
        real        u_slip, T40,Tcvot, h_real, tau_xzl_max
        integer     SS, NYs, t
        real*8      dxss_Phb
        real        tau_lim_c, tau_gamma_r
        save        /Minimum_values/
        save        /Load_surf/
        save        /Stress_under/
        save        /Stress_surf/
        save        /Nodes_no_pressure/
        save        /CurrentP/
        include     'inc_load_surf.h'
        include     'inc_CurrentP.h'
        include     'inc_CurrentH.h'
        include     'inc_CurrentT.h'
        include     'inc_grid.h'
        include     'inc_Stress_under.h'
        include     'inc_Stress_surf.h'
        include     'inc_outp.h'
        include     'inc_Current.h'
        include     'inc_Method.h'
        include     'inc_asp.h'
        include     'inc_Minimum_values.h'
        include     'inc_Nodes_no_pressure.h'
        include     'inc_Visc.h'
        include     'inc_Temp_reduction.h'
        
    pi=3.14159265
    NN=(NY+1)/2
    ub=2*um-Ua
    u_slip=ua-ub         !ua has to be the top one
    SS=1
    NYs=Ny
    T40=40
    Tcvot=1
    t=1  
    
    dxss_Phb=1.0/(2.0*DX*SS)*Ph/b
    
    !Reset stress arrays
    sig_xs   =0.0
    sig_ys   =0.0
    sig_zs   =0.0
    tau_xzs  =0.0
    tau_yzs  =0.0
    tau_xys  =0.0
    
    ! Resetting the stress arrays
    sig_x  =0.0
    sig_y  =0.0
    sig_z  =0.0
    tau_xz =0.0
    tau_yz =0.0
    tau_xy =0.0
    
      
    CALL Newtoinan(SS, NYs, T40, Tcvot, NN, t, 0, 0)
    
    !dPdy= minval(H(1:200,1:49)) ! For debugging
    ! Calculate the shear stresses tau_xz and tau_xy
    DO j=1,NN 
        J0=J-SS
        IF( J0 .LE. 0)  J0=J+SS
        J2=J-2*SS
        IF( J2 .LE. 0)  J2=J0
        J1=J+1*SS
        IF( J1 .GT. NN) J1=NN-SS 
        
        DO i=1,NX
        ! Interested in the whole area
        IF (H(I,J) .LE. Hminimum*1.1 ) THEN     !Include a tiny margine since the EHL program might add a bit to high pressure
            IF( U_slip .NE. 0) then
                tau_xzl(i,j) = sign(0.3*P(i,j),U_slip)
            else
                tau_xzl(i,j) = 0 !If pure rolling no friction in metal contact
            endif
            
            !h_real=P(i,j)          ! Only for debugging
            !dPdX = tau_xzl(i,j)     ! Only for debugging
            tau_yzl(i,j) = 0.0
                
        else
            
            I0=I-1*SS
            IF( I0 .LE. 0)  I0=1
            I2=I-2*SS
            IF( I2 .LE. 0)  I2=I0
            I1=I+1*SS
            IF( I1 .GE. NX) I1=NX


            dPdx=0.5*(P(I1,J)-P(I0,J))*dxss_Phb!  (term4*P(I1,J)+term3*P(I,J)+term2*P(I0,J)+term1*P(I2,J))/(2*DX*SS)*Ph/b
            dPdy=0.5*(P(I,J1)-P(I,J0))*dxss_Phb!  (term4*P(I,J1)+term3*P(I,J)+term2*P(I,J0)+term1*P(I,J2))/(2*DX*SS)*Ph/b
                
            h_real=H(I,J)*b**2/Rx
            tau_xzl(i,j) = ( h_real/2*dPdx+EDAX(I,j)*u_slip/h_real)/Ph      ! The stresses ar calculated according to Article 1 where the asperity is on the top surface and z is directed into the material. 
            tau_yzl(i,j) = ( h_real/2*dPdy)/Ph
            
            
            ! Check so shear stress not above shear limit
            if (temp_red .GT. 0) then                                       ! Take calculate temperature reductio of the shear limit
                tau_gamma_r = tau_gamma - temp_fac * temp(I,J)              ! tau_gamma and temp_fac anready includes PH
                if( tau_gamma_r .LT. temp_red) then                         ! Just to make sure that not negative
                    tau_gamma_r = temp_red
                endif
                        
            else
                tau_gamma_r = tau_gamma
            endif
                    
            if (temp_param==2)then
                tau_lim_c = tau_lim + tau_gamma_r * P(i,j)
            elseif ( temp_param == 6) then
                tau_lim_c = tau_lim + tau_gamma_r * P(i,j) + tau_xzl(i,j)*p_red
            endif
                    
            if (abs(tau_xzl(i,j)) .gt. tau_lim_c .and. temp_param .NE. 0) then
                tau_xzl(i,j)=sign( tau_lim_c , h_real/2*dPdx)               ! Take the sign of the pressure derivative
            endif
            
                        
            
            if( Temp_param == 0)then    !The viscosity gets to high close to the contact
                if(abs(tau_xzl(i,j)) .GT. 0.3*P(i,j)) then
                    tau_xzl(i,j) = sign( 0.3*P(i,j) , tau_xzl(i,j) )
                endif
                
                if(abs(tau_yzl(i,j)) .GT. 0.3*P(i,j)) then
                    tau_yzl(i,j) = sign( 0.3*P(i,j) , tau_yzl(i,j) )
                endif
            endif
            
                    
        ENDIF
            
        IF (H(I,J) .GT. Hminimum*1.1 .and. H(I,J) .LT. Hminimum*1.5 ) THEN     !Include a transistion zone between contact and no contact
            IF( U_slip .NE. 0) then
                tau_xzl(i,j) = sign(0.3*P(i,j),U_slip)*( 1- (H(i,j)-1.1*Hminimum ) / (0.4*hminimum) )  +  tau_xzl(i,j)*( (H(i,j)-1.1*Hminimum) / (0.4*hminimum) )
            else
                !tau_xzl(i,j) = 0 !If pure rolling do not adjust the shear stress
            endif
        ENDIF
        
            
        !IF ( abs(tau_xzl(i,j)) .GT. 0.3*max(P(i0,j), P(i,j), P(i1,j))) THEN
        !    if( P(i,j) .GT. 0.01) then !The sear stress is also limited in the outlet region due to cavitaions
        !        WRITE(4,*)'To high tau_xz stresses at i,j = ', i, j
        !        WRITE(4,*)' Old stress was = ', tau_xzl(i,j) 
        !    endif
        !    
        !    tau_xzl(i,j) =     sign(0.3*max(P(i0,j), P(i,j), P(i1,j)) , U_slip)
        !    WRITE(4,*)' Limit stress is = ', 0.3*max(P(i0,j), P(i,j), P(i1,j))
        !endif
        !IF ( abs(tau_yzl(i,j)) .GT. 0.3*max(P(i,j0), P(i,j), P(i,j1))) THEN
        !    if (P(i,j) .gt. 0.01) then
        !        WRITE(4,*)'To high tau_yz stresses at i,j = ', i, j
        !        WRITE(4,*)' Old stress was = ', tau_yzl(i,j) 
        !    endif
        !    tau_yzl(i,j) =     sign(0.3*max(P(i,j0), P(i,j), P(i,j1)) , tau_yzl(i,j))
        !    WRITE(4,*)' Limit stress is = ', 0.3*max(P(i,j0), P(i,j), P(i,j1))
        !endif
             
        ENDDO
    ENDDO
    
  ! P=0*P
  !tau_yzl=0*tau_yzl
  !tau_xzl=0*tau_xzl
  !!
  !tau_xzl(100,20)=1;
  !tau_xzl(100,21)=1;
  !tau_xzl(100,22)=1;
  !tau_xzl(101,20)=1;
  !tau_xzl(101,21)=2;
  !tau_xzl(101,22)=1;
  !tau_xzl(102,20)=1;
  !tau_xzl(102,21)=1;
  !tau_xzl(102,22)=1;
  !
  !tau_xzl(200,20)=1;
  !tau_xzl(200,21)=1;
  !tau_xzl(200,22)=1;
  !tau_xzl(201,20)=1;
  !tau_xzl(201,21)=1;
  !tau_xzl(201,22)=1;
  !tau_xzl(202,20)=1;
  !tau_xzl(202,21)=1;
  !tau_xzl(202,22)=1;
  do j=1,NY
      jj=NY+1-j
      do i=1,NX
          P(i,jj)=P(i,j)
          tau_xzl(i,jj)=tau_xzl(i,j)
          tau_yzl(i,jj)=-tau_yzl(i,j)
      enddo
      
  enddo
  
    nstop=0
    NXstart=2
    NXstop=NX-1
    DO i=2,NX
        ! Finds the minimum sear traction for each line
        tau_xzl_min(i)      = sign(minval(abs(tau_xzl(i,2:NY-1))),U_slip) ! The sign is changed during load calculations
        tau_xzl_max         = maxval(abs(tau_xzl(i,2:NY-1)))
        !t2=tau_xzl(i,2:NY)
        !t3=P(i,2:NY)!
        !t4=H(i,2:NY)!
        !tau_yzs_min(i)      =sign(minval(abs(tau_yzs(i,2:NY))),minval(tau_yzs(i,2:NY))) !Shold be zero due to geometry of the problem
        
        ! Finds the minimum lineload for each line
        Pmin(i)      =minval(P(i,2:NN))
        Pmax(i)      =maxval(P(i,2:NN))
        
        ! To limit the investigated region to where the pressure is above a surtain level
        IF (Pmax(i) .LT. 0.001 .and. abs(tau_xzl_max) .LT. 0.001 .and. nstop .Eq. 0) then
            NXstart=i
            nstop=0
        else IF (Pmax(i) .LT. 0.001 .and. abs(tau_xzl_max) .LT. 0.001 .and. nstop .Eq. 1) then
            NXstop=i
            nstop=2
        else IF (Pmax(i) .GT. 0.001 .or. abs(tau_xzl_max) .GT. 0.001 ) then
            nstop=1
            NXstop=i
        endif
        ! Easy way to disable the lineloads
        !Pmin(i)      =0*Pmin(i)
        !Pmax(i)      =0*Pmax(i)
        !tau_xzs_min(i) =0* tau_xzs_min(i) 
    enddo
    ! Add a bit of margine . 
    NXstop=NXstop+2
    NXstart=NXstart-2
    ! Ensure that within resonable range
    iF( NXstop .GT. NX-1) NXstop=NX-1
    if( NXstart .LT. 2) NXstart=2
    !Pmin=0
    return
    end
    