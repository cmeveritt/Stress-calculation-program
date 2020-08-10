! Subrutine for updating the film thicknes, the viscocity and the density based on the given pressure
	SUBROUTINE Newtoinan(SS, NYs, T40, Tcvot, NN, t, k, M_conv)
        implicit none 
        include     'inc_Current.h'
	    include     'inc_CurrentH.h'
	    include     'inc_CurrentP.h'
	    include     'inc_CurrentT.h'
        include     'inc_Contact_mat.h'
        
        include     'inc_EDA_contact.h'
        include     'inc_Grid.h'
        include     'inc_G0DT.h'
        include     'inc_Holmes.h'
        
	    include     'inc_Method.h'
        
        include     'inc_Outp.h'
        include     'inc_Ref.h'
        include     'inc_Rho.h'
        include     'inc_RLarsson.h'
        
        include     'inc_Shear_lim.h'
	    include     'inc_Temp_param.h'
        include     'inc_Temp_reduction.h'
        include     'inc_Visc.h'
        include     'inc_Yasutomi.h'
	    include     'inc_Y_liu.h'
 
        ! Input
        integer     NN, NYs, SS, t, k, M_conv
        real        T40, Tcvot

        ! Calculations
        integer     I,J, JJ, temp_iter, err, printer
        real        EDA_1, EDA_2, EDA_3
        real        EDA22, EDA33, EDA1
        real        pg, Tg, YF 
        real        temp_0, u_slip, h_real, tau
        real        tau_x_f, tau_y_f, tau_x_c, tau_y_c, tau_lim_c
        real        dPdx, dPdy, dpds
        integer     J0, J1, J2, I0, I1, I2
        real        eda01,eda02,eda03,eda04,delta_s
        real        tau_gamma_r, dxss_Phb
        real        epst

        ! Output
        SAVE      /Current/                    ! Current timestep
        SAVE      /CurrentT/
        
          
        
    ! Calculate the viscocity, the density and the dimentionles parameter EPS -----------------------------------------------------------------------------------------
    IF( Lub_param .EQ. 1) THEN
        ! Newtonian acc X.Tans paper ref 30 31, Cylinder
        DO J=1,NN,SS
            DO I=1,NX,SS

                ! Roelands equation acc X.Tan and Venner
                EDAx(I,J)=EXP(alpha*Pref/Z*(-1+(1+P(I,J)*PH/Pref)**Z)) 
                if (EDAx(I,J) .GT. 1E32) EDAx(I,J)=1E32                 ! To limit the viscosity to som reasonable? values. Intended to increase convergense. 
                EDAy(I,J)=EDAx(I,J)
                
            
                ! D-H Formulation acc X.Tan
                 !RO(I,J)=(5.9E8+1.34*P(I,J)*PH)/(5.9E8+P(I,J)*PH) 
                 RO(I,J)=1+RA1*PH*P(I,J)/(1+RA2*PH*P(I,J)) 
                 xi(I,J)=0.0                                        ! Should be able to remove xi becouse we're never changing it
                
                
                IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                    WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'P = ', P(I,J),  'EDAO = ',EDA0
                    EDAx(I,J)=0.1
                    Call Stop_to_large_out(t)
                ENDIF

                IF ( RO(I,J) .LT. 1.0 .OR. isnan(RO(I,J))) THEN
                    WRITE(4,*)'BAD RO', RO(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT
                    RO(I,J)=1.0
                    Call Stop_to_large_out(t)
                ENDIF
        
                IF(contact(I,J) .EQ. 2 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                EPSy(I,J)=EPSx(I,J)!RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))
                
                
            ENDDO
        ENDDO
      ELSE IF( Lub_param .EQ. 11) THEN
        ! Barus equation with extra Z
        DO J=1,NN,SS
            DO I=1,NX,SS

                ! Roelands equation acc X.Tan and Venner
                EDAx(I,J)=EXP(alpha*(P(I,J)*PH)**Z)
                if (EDAx(I,J) .GT. 1E32) EDAx(I,J)=1E32                 ! To limit the viscosity to som reasonable? values. Intended to increase convergense. 
                EDAy(I,J)=EDAx(I,J)
                
            
                ! D-H Formulation acc X.Tan
                 !RO(I,J)=(5.9E8+1.34*P(I,J)*PH)/(5.9E8+P(I,J)*PH) 
                 RO(I,J)=1+RA1*PH*P(I,J)/(1+RA2*PH*P(I,J)) 
                 xi(I,J)=0.0
                
                
                
                IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                    WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0
                    EDAx(I,J)=0.1
                    Call Stop_to_large_out(t)
                ENDIF

                IF ( RO(I,J) .LT. 1.0 .OR. isnan(RO(I,J))) THEN
                    WRITE(4,*)'BAD RO', RO(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT
                    RO(I,J)=1.0
                    Call Stop_to_large_out(t)
                ENDIF
        
                IF(contact(I,J) .EQ. 2 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                EPSy(I,J)=EPSx(I,J)!RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))
                
                
            ENDDO
        ENDDO  
    ELSEIF( Lub_param .EQ. 2) THEN
        ! Newtonian acc Gohars book Elastohydrodynamics
        
        DO J=1,NN,SS
            DO I=1,NX,SS

                ! Roelands equation acc Gohar
                EDA_1=5.1*PH*P(I,J)/10**(9)
                EDA_2=(1.0+EDA_1)**Z
                EDA_3=(log(EDA0)+9.67)
                EDAx(I,J)=EXP(EDA_3*(EDA_2-1.0))
                
                if (EDAx(I,J) .GT. 1E32) EDAx(I,J)=1E32                 ! To limit the viscosity to som reasonable? values. Intended to increase convergense. 
                EDAy(I,J)=EDAx(I,J)
                
            
                ! D-H Formulation acc X.Tan
                 RO(I,J)=1+RA1*PH*P(I,J)/(1+RA2*PH*P(I,J)) 
                 xi(I,J)=0.0
                
                
                IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                    WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0
                    EDAx(I,J)=0.1
                    Call Stop_to_large_out(t)
                ENDIF

                IF ( RO(I,J) .LT. 1.0 .OR. isnan(RO(I,J))) THEN
                    WRITE(4,*)'BAD RO', RO(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT
                    RO(I,J)=1.0
                    Call Stop_to_large_out(t)
                ENDIF
        
                IF(contact(I,J) .EQ. 2 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                EPSy(I,J)=EPSx(I,J)!RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))
                
                
        ENDDO
        ENDDO
        
    ELSE IF( Lub_param .EQ. 3)THEN
        ! Newtonian input parameters ac Holmes et al Transient EHL point contact analysis 2003, Ball 
        DO J=1,NN,SS
            DO I=1,NX,SS

                ! Roelands equation acc Holmes
                 EDAx(I,J)=EXP(LOG(EDA0/kH)*((1+xH*PH*P(I,J))**Z-1)) 
                 EDAy(I,J)=EDAx(I,J) !EXP(LOG(EDA0/kH)*((1+xH*PH*P(I,J))**Z-1)) 
                
                ! D-H Formulation acc Homes et al
                 !RO(I,J)=(1 + gH*PH*P(I,J))/(1+lH*P(I,J)*PH) 
                 RO(I,J)=1+RA1*PH*P(I,J)/(1+RA2*PH*P(I,J)) 
                 xi(I,J)=0.0
                
                
                IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                    WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0
                    EDAx(I,J)=0.1
                    Call Stop_to_large_out(t)
                ENDIF

                IF ( RO(I,J) .LT. 1.0 .OR. isnan(RO(I,J))) THEN
                    WRITE(4,*)'BAD RO', RO(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT
                    RO(I,J)=1.0
                    Call Stop_to_large_out(t)
                ENDIF
        
                IF(contact(I,J) .EQ. 2 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                EPSy(I,J)=EPSx(I,J)!RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))
        ENDDO
    ENDDO
        
    ELSE IF( Lub_param .EQ. 4 .and. lub_temp .GE. 0)THEN      ! No temperature adjustments
        ! R.Larssons 2000 formulation
        ! This part is copied to EDA_calc
        if ( temp_param == 2 ) then
            ub= 2*um-ua
            u_slip=ua-ub
        endif
        printer=0
        dxss_Phb=1.0/(2*DX*SS)*Ph/b 
        
        DO J=1,NN,SS
            DO I=1,NX,SS
899             err=0
                If ( Temp_param == 3) then
                    RL_Ta=temp(1,1) +( temp(I,J) - temp(1,1))*temp_fac
                    IF( RL_Ta .LT. temp(1,1)) then
                        RL_ta=temp(1,1)
                    endif
                    
                else
                    RL_Ta=temp(I,J)
                endif
                
                ! Roelands equation acc R.Larsson
                EDA0 = 10**(-4.2+RL_G0*(1.0+RL_Ta/135.0)**S0)             ! Eq (3)
                Z=Dz+Cz*log10(1.0+RL_Ta/135)
                
                EDA1=EDA0*EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
                
                
                IF( U_slip == 0) then
                    EDAx(i,j) = EDA1
                    
                elseif ( Temp_param == 1) then
                    ! Reducing eta
                    EDAx(i,j) = EDA1 *  1/ (1 + temp_fac * max( 0.0, (P(i,j)-p_red)/p_red) * ( max(0.0 , (RL_Ta-temp_red)/temp_red ) ) )    
                    
                elseif( temp_param == 2 .or. temp_param==6) then
                    ! applying a shear limmit  based on the pressure according to data from Höglund 1989 and 1999 with the temperature reduction of the shear limit
                    h_real=H(I,J)*b**2/Rx
                    
                    ! Define nodenumbers for past and next nodes
                    J0=J-SS
                    IF( J0 .LE. 0)  J0=J+SS
                    J2=J-2*SS
                    IF( J2 .LE. 0)  J2=J0
                    J1=J+1*SS
                    IF( J1 .GE. NYs) J1=NYs-SS
                    JJ=NYs+1-J
                    IF( JJ .LE. 1)  JJ=1+SS
                
                    I0=I-1*SS
                    IF( I0 .LE. 0)  I0=1
                    I2=I-2*SS
                    IF( I2 .LE. 0)  I2=I0
                    I1=I+1*SS
                    IF( I1 .GE. NX) I1=NX
                    
                    dPdx=(term4*P(I1,J)+term3*P(I,J)+term2*P(I0,J)+term1*P(I2,J))* dxss_Phb ! /(2*DX*SS)*Ph/b
                    tau = abs(H_real/2*dpdx) + abs(EDA1*U_slip/H_real) 
                    
                    dPdy=(term4*P(I,J1)+term3*P(I,J)+term2*P(I,J0)+term1*P(I,J2))* dxss_Phb 
                    tau= sqrt( tau**2 + (H_real/2*dpdy)**2)
                    
                    if (temp_red .GT. 0) then                                       ! Take calculate temperature reductio of the shear limit
                        tau_gamma_r = tau_gamma - temp_fac * RL_TA     ! tau_gamma and temp_fac anready includes PH
                        if( tau_gamma_r .LT. temp_red) then                             ! Just to make sure that not negative
                            tau_gamma_r = temp_red
                        endif
                        
                    else
                        tau_gamma_r = tau_gamma
                    endif
                    
                    if (temp_param==2)then
                        tau_lim_c = tau_lim + tau_gamma_r * P(i,j)
                    elseif ( temp_param == 6) then
                        tau_lim_c = tau_lim + tau_gamma_r * P(i,j) + tau*p_red
                    endif
                    
                    if (tau .gt. tau_lim_c) then
                        dpds= sqrt( dPdx**2 + dPdy**2)
                        
                        if( tau_lim_c - abs(H_real/2*dpds) .LT. lim_factor*tau_lim ) then
                            
                            ! To increase numerical stability the maximum pressure of the neighbouring nodes is used
                            if (temp_param==2)then
                                tau_lim_c = tau_lim + tau_gamma_r *max( P(i0,j),P(i,j),P(i1,j), P(i,j0), P(i,j1))
                            elseif ( temp_param == 6) then
                                tau_lim_c = tau_lim + tau_gamma_r *max( P(i0,j),P(i,j),P(i1,j), P(i,j0), P(i,j1)) + tau*p_red
                            endif
                            
                            if( tau_lim_c - abs(H_real/2*dpds) .LT. lim_factor*tau_lim ) then
                                EDA1 = abs( (lim_factor*tau_lim)*H_real/U_slip)    ! Adding a factor of 2 here in hope of increased numerical stability
                                printer=printer+1
                            else
                                EDA1 = abs( (tau_lim_c - abs(H_real/2*dpds))*H_real/U_slip)
                            endif
                            
                        else
                            EDA1 = abs( (tau_lim_c - abs(H_real/2*dpds))*H_real/U_slip)
                        endif
                    endif
                    
                    EDAx(i,j) = EDA1
                elseif( temp_param== 4 .or. temp_param== 5) then ! affects the empt_calc_ij subroutine
                    EDAx(i,j) = EDA1
                else
                    EDAx(i,j) = EDA1
                endif
                
                ! D-H Formulation acc R.Larsson
                EpsT       = EpsT0*exp(-RL_c*P(I,J)*PH)
                RO(I,J)    = (1.0+RA1*PH*P(I,J)/(1.0+RA2*PH*P(I,J)))*(1.0-EpsT*(RL_Ta-RL_T0)) 
                !xi(I,J)    = 0.0
                
                
                IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) .or. (EDAx(i,j) .gt. 1e20) ) THEN
                    err=err+1
                    if (err==1) then 
                        WRITE(4,*)'BAD EDAx =', EDAx(I,J), '. For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'temp = ',temp(i,j)
                        I0=I-1*SS
                        IF( I0 .LE. 0)  I0=1
                        temp(i,j)=temp(i0,j)
                        go to 899
                    else
                        WRITE(4,*)'BAD EDAx =', EDAx(I,J), ', after temperature adjustment'
                        Write(4,*)'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'temp = ',temp(i,j)
                        err=1
                        EDAx(I,J)=0.1
                    endif
                    Call Stop_to_large_out(t)
                ENDIF

                IF ( RO(I,J) .LT. 0.0 .OR. isnan(RO(I,J))) THEN
                    err=err+1
                    if (err==1) then
                        WRITE(4,*)'BAD RO =', RO(I,J), '. For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT, 'temp = ',temp(i,j)
                        Call Stop_to_large_out(t)
                        I0=I-1*SS
                        IF( I0 .LE. 0)  I0=1
                        temp(i,j)=temp(i0,j)
                        go to 899
                    else
                        WRITE(4,*)'BAD RO =', RO(I,J), ', after temperature adjustment'
                        WRITE(4,*)'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT, 'temp = ',temp(i,j)
                        RO(I,J)=1.0
                    endif
                    Call Stop_to_large_out(t)
                ENDIF
                
                IF(contact(I,J) .EQ. 2 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                EPSy(I,J)=EPSx(I,J)!RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))
                
            ENDDO
        ENDDO
     if (printer .gt. 100 .and. k == 1) then
            WRITE(4,*)'Bad equation for shear limmit ', printer ,' number of times'
            Call Stop_to_large_out(t)
            if ( printer .gt. 800 .and. SS == 1 .and. M_conv .gt. 10) then  ! This is a clear sign that the code is diverging. Howevere, at the first itterations of a timestep it might not be true. 
                WRITE(4,*) 'Breaking due to many bad shear limit equations'
                WRITE(4,*) 'The time is t=',t
                OPEN(20,FILE='STOP_Too_many_bad_Shear.DAT',STATUS='UNKNOWN')
                stop
            endif
            
            
    endif
    
    ELSE IF(lub_param .EQ. 4 .and. lub_temp == -1) then ! simple Themperature rise
        
        ub= 2*um-ua
        u_slip=ua-ub
        
        DO J=1,NN,SS
            DO I=1,NX,SS
                temp_iter=0             ! Resetting counter
                temp_0=temp(1,1)                        ! The global temperature
                
                ! Ensure that the lubrication can not cool off downstream. 
                !IF( L_stab .EQ. 2 .AND. t .LE. -2 .and. I .GT. 1 ) Then
                !    temp(I,j)=max(temp(i,j),temp(I-SS,J))       ! Do not update if higher temp 
                !    temp_0=temp(I-SS,J)       
                !ENDIF
            
777             RL_Ta=temp(I,J)         ! Extract the temperature
                
                ! Roelands equation acc R.Larsson
                EDA0 = 10**(-4.2+RL_G0*(1+RL_Ta/135.0)**S0)             ! Eq (3)
                Z=Dz+Cz*log10(1.0+RL_Ta/135.0)
                
                EDAx(i,j)=EDA0*EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
                EDAy(I,J)=EDAx(I,J)!EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
                
                ! D-H Formulation acc R.Larsson
                EpsT=EpsT0*exp(-RL_c*P(I,J)*PH)
                RO(I,J)=(1.0+RA1*PH*P(I,J)/(1.0+RA2*PH*P(I,J)))*(1.0-EpsT*(RL_Ta-RL_T0)) 
                xi(I,J)=0.0
                
                
                IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                    WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0
                    EDAx(I,J)=0.1
                    Call Stop_to_large_out(t)
                ENDIF

                IF ( RO(I,J) .LT. 0.0 .OR. isnan(RO(I,J))) THEN
                    WRITE(4,*)'BAD RO', RO(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT
                    RO(I,J)=1.0
                    Call Stop_to_large_out(t)
                ENDIF
        
                IF(contact(I,J) .EQ. 2 .and. EDAx(I,J) .LT. EDA_cont) then
                    EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                    
                ELSE IF( temp_iter .LE. 10 .AND. SS .LE. 2) THEN     ! Might be good to introduce this after the coursest solution is obtained. 
                    
                    ! Define nodenumbers for past and next nodes
                    J0=J-SS
                    IF( J0 .LE. 0)  J0=J+SS
                    J2=J-2*SS
                    IF( J2 .LE. 0)  J2=J0
                    J1=J+1*SS
                    IF( J1 .GE. NYs) J1=NYs-SS
                    JJ=NYs+1-J
                    IF( JJ .LE. 1)  JJ=1+SS
                
                    I0=I-1*SS
                    IF( I0 .LE. 0)  I0=1
                    I2=I-2*SS
                    IF( I2 .LE. 0)  I2=I0
                    I1=I+1*SS
                    IF( I1 .GE. NX) I1=NX
                
                    dPdx=(term4*P(I1,J)+term3*P(I,J)+term2*P(I0,J)+term1*P(I2,J))/(2*DX*SS)*Ph/b
                    dPdy=(term4*P(I,J1)+term3*P(I,J)+term2*P(I,J0)+term1*P(I,J2))/(2*DX*SS)*Ph/b
                
                    h_real=H(I,J)*b**2/Rx
                    tau_x_f = abs(-h_real/2*dPdx+EDAX(I,j)*u_slip/h_real)
                    tau_y_f = abs(-h_real/2*dPdy)
                    tau_x_c = abs( h_real/2*dPdx+EDAX(I,j)*u_slip/h_real)
                    tau_y_c = abs( h_real/2*dPdy)
                    
                    
                    if( (tau_x_f .GT. shear_max .or. tau_y_f .gt. shear_max .or. tau_x_c .gt. shear_max .OR. tau_y_c .GT. shear_max) .AND. temp(i,j) .LT. temp_max) then
                        
                        temp(I,J)=temp(I,J)+5
                        RL_Ta=temp(I,J)
                        ! Roelands equation acc R.Larsson
                        EDA0 = 10**(-4.2+RL_G0*(1+RL_Ta/135.0)**S0)             ! Eq (3)

                        EDAx(I,J)=EDA0*EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
                        EDAy(I,J)=EDAx(I,J)!EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
                
                        ! D-H Formulation acc R.Larsson
                        EpsT=EpsT0*exp(-RL_c*P(I,J)*PH)
                        RO(I,J)=(1.0+RA1*PH*P(I,J)/(1.0+RA2*PH*P(I,J)))*(1.0-EpsT*(RL_Ta-RL_T0)) 
                
                
                        IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                            WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0
                            EDAx(I,J)=0.1
                            Call Stop_to_large_out(t)
                        ENDIF

                        IF ( RO(I,J) .LT. 0.0 .OR. isnan(RO(I,J))) THEN
                            WRITE(4,*)'BAD RO', RO(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT
                            RO(I,J)=1.0
                            Call Stop_to_large_out(t)
                        ENDIF
        
                        IF(contact(I,J) .EQ. 2 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                        EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                        EPSy(I,J)=EPSx(I,J)!RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))
                        
                        
                        
                        temp_iter=temp_iter+1
                        GO TO 777
                    elseif( (tau_x_f .LT. shear_min .AND. tau_y_f .LT. shear_min .AND. tau_x_c .LT. shear_min .AND. tau_y_c .LT. shear_min) .AND. temp(i,j) .GT. temp_0) then
                        
                        temp(I,J)=temp(I,J)-2
                        RL_Ta=temp(I,J)
                        ! Roelands equation acc R.Larsson
                        EDA0 = 10**(-4.2+RL_G0*(1.0+RL_Ta/135.0)**S0)             ! Eq (3)
                        ENDA=12  *Um * RX**2 / ( B**3*PH) ! *EDA0 but EDA0 is included in EPSx

                        EDAx(I,J)=EDA0*EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
                        EDAy(I,J)=EDAx(I,J)!EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
                
                        ! D-H Formulation acc R.Larsson
                        EpsT=EpsT0*exp(-RL_c*P(I,J)*PH)
                        RO(I,J)=(1.0+RA1*PH*P(I,J)/(1.0+RA2*PH*P(I,J)))*(1.0-EpsT*(RL_Ta-RL_T0)) 
                
                
                        IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                            WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0
                            EDAx(I,J)=0.1
                            Call Stop_to_large_out(t)
                        ENDIF

                        IF ( RO(I,J) .LT. 0.0 .OR. isnan(RO(I,J))) THEN
                            WRITE(4,*)'BAD RO', RO(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT
                            RO(I,J)=1.0
                            Call Stop_to_large_out(t)
                        ENDIF
        
                        IF(contact(I,J) .EQ. 2 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                        EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                        EPSy(I,J)=EPSx(I,J)!RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))

                        
                        temp_iter=temp_iter+1
                        GO TO 777
                    endif
                ENDIF 
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                EPSy(I,J)=EPSx(I,J)!RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))


            ENDDO
        ENDDO
        
        !Mirroring the temperature
        DO J=1,NN,SS
            JJ=NYs-J+1
            DO I=1,NX,SS   
                temp(I,JJ)=temp(I,J)
            ENDDO
        ENDDO   
        
        
    ELSE IF( Lub_param .EQ. 5) THEN
        ! Yasutomi lubrication
        DO J=1,NN,SS
            DO I=1,NX,SS

                ! Yasutomi equations
                Tg=Tg0+YA1*LOG(1.0+YA2*PH*P(I,J))
                pg=1.0/YA2*(EXP((1.0/YA1)*((T40*(1.0-Tcvot)+Tcvot*Temp(I,J))-Tg0))-1.0)
                YF=1.0-YB1*LOG(1.0+YB2*PH*P(I,J))
                
                ! Since for some specific combinations of p and YF EDA33 gest biger than EDA22 eaven do p < pg this formulation should be better. 
                IF (P(I,J)*PH .GT. pg) THEN
                    EDA1=Yedag*EXP(Yalfag*(P(I,J)*PH-pg))
                ELSE if(Temp(i,j) .LT. TG) then
                    EDA1=Yedag*EXP(Yalfag*(P(I,J)*PH-pg))
                    
                    WRITE(4,*)'Selfcomponed case for temp < Tg. Temp and Tg are', Temp, Tg
                    WRITE(4,*)'at i,j, =', i,j
                    Call Stop_to_large_out(t)
                    
                ELSE
                    EDA33=Yedag*10**(-(YC1*((T40*(1-Tcvot)+Tcvot*Temp(I,J))-Tg)*YF) / (YC2+((T40*(1-Tcvot)+Tcvot*Temp(I,J))-Tg)*YF)) 
                    
                    EDA22=Yedag*EXP(Yalfag*1)
                    if (EDA22 .lt. EDA33) then
                    EDA1 = EDA22
                    else
                    EDA1 = EDA33
                    endif
                ENDIF
                
                EDAx(I,J)=EDA1                                              ! No Non-newtonian reduction. 
                
                 IF (EDAx(I,J) .LE. 0.001*EDA0 ) THEN ! With themp increase outside of pressure zone EDAx can decrease
                     IF( J .GT. 5) then
                        WRITE(4,*)'BAD EDAx. EDAx lt 0.1*EDA0. EDAx=', EDAx(I,J), 'For I J = ', I ,J ,'EDA1 = ', EDA1
                        WRITE(4,*)'EDA33 = ', EDA33, ' EDA22 = ', EDA22 ,'Tcvot = ', Tcvot, ' Temp(i,j) = ', temp(i,j)
                     ENDIF
                     
                    EDAx(I,J)=0.001*EDA0
                    Call Stop_to_large_out(t)
                ENDIF
                
                IF( isnan(EDAX(i,j))) then
                    EDAX(i,j)=EDA1
                endif
                IF( (EDAX(i,j)) .gt. 1e15) then
                    !temp(1,1)=temp(i,j)
                    !P(1,1)=P(i,j)
                    EDAX(i,j)=1e15
                endif
            
                ! D-H Formulation acc P.Ehret D.Dowsin and C.M. Taylor
                RO(I,J)=1+RA1*PH*P(I,J)/(1+RA2*PH*P(I,J)) 
                
                IF(contact(I,J) .EQ. 2 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                EPSy(I,J)=EPSx(I,J)!RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                
               
                
            ENDDO
        ENDDO
    ELSE IF( Lub_param .EQ. 6) THEN
        ! Newtonian acc N Deolalikers contact paper from 2008
        ! Does not exist in Lubrication_def
        DO J=1,NN,SS
            DO I=1,NX,SS

                ! Roelands equation acc X.Tan and Venner
                EDAx(I,J)=EXP(log(EDA0+9.67)*((1+P(I,J)*PH/Pref)**Z-1)) ! Should include EDA0 due to fomulation in main
                if (EDAx(I,J) .GT. 1E32) EDAx(I,J)=1E32                 ! To limit the viscosity to som reasonable? values. Intended to increase convergense. 
                EDAy(I,J)=EDAx(I,J)
                
                ! D-H Formulation acc X.Tan and N Deolaliker
                ! RO(I,J)=(5.9E8+1.34*P(I,J)*PH)/(5.9E8+P(I,J)*PH) 
                 RO(I,J)=1+RA1*PH*P(I,J)/(1+RA2*PH*P(I,J)) 
                 xi(I,J)=0.0
                
                
                IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                    WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0
                    EDAx(I,J)=0.1
                    Call Stop_to_large_out(t)
                ENDIF

                IF ( RO(I,J) .LT. 1.0 .OR. isnan(RO(I,J))) THEN
                    WRITE(4,*)'BAD RO', RO(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT
                    RO(I,J)=1.0
                    Call Stop_to_large_out(t)
                ENDIF
        
                IF(contact(I,J) .EQ. 2 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                EPSy(I,J)=EPSx(I,J)!RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))
                
                
            ENDDO
        ENDDO
        
    ELSE IF( Lub_param .EQ. 7 .or. Lub_param .EQ. 8)THEN      
        ! Wang 2015 formulation  with temperature calculations or Liu 2005 formulation  with temperature calculations
        
         DO J=1,NN,SS
            J0=J-SS
            IF( J0 .LE. 0)  J0=J+SS
            J2=J-2*SS
            IF( J2 .LE. 0)  J2=J0
            J1=J+1*SS
            IF( J1 .GE. NYs) J1=NYs-SS
            JJ=NYs+1-J
            IF( JJ .LE. 1)  JJ=1+SS
                    
            DO I=1,NX,SS
                err=0
888             RL_Ta=temp(I,J)

                eda01=log(EDA0)+9.67
                eda02=(1.0+5.1E-9*P(I,J)*PH)**Z
                eda03=((RL_Ta-138.0)/(RL_T0-138.0))**(-S0)
                
                EDAx(I,J)=EDA0*exp(eda01 * (-1.0 + eda02 * eda03))
                
                !EDAy(I,J)=EDAx(I,J)!EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
                
                ! D-H Formulation acc R.Larsson
                RO(I,J)=1.0 + RA1*PH*P(I,J)/(1.0+RA2*PH*P(I,J)) - EpsT0*(RL_Ta-RL_T0) 
                
                
                IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                    err=err+1
                    if (err==1) then 
                        WRITE(4,*)'BAD EDAx =', EDAx(I,J), '. For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'temp = ',temp(i,j)
                        I0=I-1*SS
                        IF( I0 .LE. 0)  I0=1
                        temp(i,j)=temp(i0,j)
                        go to 888
                    else
                        WRITE(4,*)'BAD EDAx =', EDAx(I,J), ', after temperature adjustment'
                        Write(4,*)'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'temp = ',temp(i,j)
                        err=1
                        EDAx(I,J)=0.1
                    endif
                    Call Stop_to_large_out(t)
                ENDIF

                IF ( RO(I,J) .LT. 0.0 .OR. isnan(RO(I,J))) THEN
                    err=err+1
                    if (err==1) then
                        WRITE(4,*)'BAD RO =', RO(I,J), '. For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT, 'temp = ',temp(i,j)
                        Call Stop_to_large_out(t)
                        I0=I-1*SS
                        IF( I0 .LE. 0)  I0=1
                        temp(i,j)=temp(i0,j)
                        go to 888
                    else
                        WRITE(4,*)'BAD RO =', RO(I,J), ', after temperature adjustment'
                        WRITE(4,*)'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT, 'temp = ',temp(i,j)
                        RO(I,J)=1.0
                    endif
                    Call Stop_to_large_out(t)
                ENDIF
        
                IF(contact(I,J) .EQ. 2 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                EPSy(I,J)=EPSx(I,J)
            enddo
         enddo
         
    ELSE IF( Lub_param .EQ. 9)THEN      
        ! Bruyerer 2012 formulation
        ! Also used for lub_param=10, Hartinger 2008
         DO J=1,NN,SS
                    J0=J-SS
                    IF( J0 .LE. 0)  J0=J+SS
                    J2=J-2*SS
                    IF( J2 .LE. 0)  J2=J0
                    J1=J+1*SS
                    IF( J1 .GE. NYs) J1=NYs-SS
                    JJ=NYs+1-J
                    IF( JJ .LE. 1)  JJ=1+SS
                    
            DO I=1,NX,SS
                err=0
855             RL_Ta=temp(I,J)

                eda01=log(EDA0)+9.67
                eda02=(1.0+P(I,J)*PH/(1.98*1E8))**Z
                eda03=((RL_Ta-138.0)/(RL_T0-138.0))
                
                delta_s=eda01*eda02*S0/(RL_T0-138.0)
                
                eda04=delta_s*(RL_Ta-RL_T0)
                
                EDAx(I,J)=EDA0*exp(eda01 * ((-1.0 + eda02 )* eda03)-eda04)
                
                !EDAy(I,J)=EDAx(I,J)!EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
                
                ! D-H Formulation acc R.Larsson
                RO(I,J)=(5.9*1E8 + RA1*PH*P(I,J))/(5.9*1E8 + RA2*PH*P(I,J)) - EpsT0*(RL_Ta-RL_T0) 
                
                
                IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                    err=err+1
                    if (err==1) then 
                        WRITE(4,*)'BAD EDAx =', EDAx(I,J), '. For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'temp = ',temp(i,j)
                        I0=I-1*SS
                        IF( I0 .LE. 0)  I0=1
                        temp(i,j)=temp(i0,j)
                        go to 855
                    else
                        WRITE(4,*)'BAD EDAx =', EDAx(I,J), ', after temperature adjustment'
                        Write(4,*)'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'temp = ',temp(i,j)
                        err=1
                        EDAx(I,J)=0.1
                    endif
                    Call Stop_to_large_out(t)
                ENDIF

                IF ( RO(I,J) .LT. 0.0 .OR. isnan(RO(I,J))) THEN
                    err=err+1
                    if (err==1) then
                        WRITE(4,*)'BAD RO =', RO(I,J), '. For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT, 'temp = ',temp(i,j)
                        Call Stop_to_large_out(t)
                        I0=I-1*SS
                        IF( I0 .LE. 0)  I0=1
                        temp(i,j)=temp(i0,j)
                        go to 855
                    else
                        WRITE(4,*)'BAD RO =', RO(I,J), ', after temperature adjustment'
                        WRITE(4,*)'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT, 'temp = ',temp(i,j)
                        RO(I,J)=1.0
                    endif
                    Call Stop_to_large_out(t)
                ENDIF
        
                IF(contact(I,J) .EQ. 2 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                EPSy(I,J)=EPSx(I,J)
            enddo
         enddo
         
    ELSE IF( Lub_param .EQ. 12 .or. Lub_param .EQ. 13)THEN      
        ! Wang 2015 formulation  with temperature calculations or Liu 2005 formulation  with temperature calculations
        
         DO J=1,NN,SS
            J0=J-SS
            IF( J0 .LE. 0)  J0=J+SS
            J2=J-2*SS
            IF( J2 .LE. 0)  J2=J0
            J1=J+1*SS
            IF( J1 .GE. NYs) J1=NYs-SS
            JJ=NYs+1-J
            IF( JJ .LE. 1)  JJ=1+SS
                    
            DO I=1,NX,SS
                err=0
889             RL_Ta=temp(I,J)

                eda01 = log(EDA0)+9.67
                eda02 = (1.0+5.1E-9*P(I,J)*PH)**Z
                eda03 = S0*(RL_Ta-RL_T0)
                
                EDAx(I,J) = EDA0*exp(eda01 * (-1.0 + eda02) - eda03)
                
                !EDAy(I,J)=EDAx(I,J)!EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
                
                ! D-H Formulation acc R.Larsson
                RO(I,J)=(1.0 + RA1*PH*P(I,J)/(1.0+RA2*PH*P(I,J)))*( 1- EpsT0*(RL_Ta-RL_T0)) 
                
                
                IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                    err=err+1
                    if (err==1) then 
                        WRITE(4,*)'BAD EDAx =', EDAx(I,J), '. For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'temp = ',temp(i,j)
                        I0=I-1*SS
                        IF( I0 .LE. 0)  I0=1
                        temp(i,j)=temp(i0,j)
                        go to 889
                    else
                        WRITE(4,*)'BAD EDAx =', EDAx(I,J), ', after temperature adjustment'
                        Write(4,*)'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'temp = ',temp(i,j)
                        err=1
                        EDAx(I,J)=0.1
                    endif
                    Call Stop_to_large_out(t)
                ENDIF

                IF ( RO(I,J) .LT. 0.0 .OR. isnan(RO(I,J))) THEN
                    err=err+1
                    if (err==1) then
                        WRITE(4,*)'BAD RO =', RO(I,J), '. For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT, 'temp = ',temp(i,j)
                        Call Stop_to_large_out(t)
                        I0=I-1*SS
                        IF( I0 .LE. 0)  I0=1
                        temp(i,j)=temp(i0,j)
                        go to 889
                    else
                        WRITE(4,*)'BAD RO =', RO(I,J), ', after temperature adjustment'
                        WRITE(4,*)'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT, 'temp = ',temp(i,j)
                        RO(I,J)=1.0
                    endif
                    Call Stop_to_large_out(t)
                ENDIF
        
                IF(contact(I,J) .EQ. 2 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                EPSy(I,J)=EPSx(I,J)
            enddo
         enddo
    else
        WRITE(4,*)'Bad lubrication number in subroutine Newtonian. lub_param ='
        WRITE(4,*) lub_param
        stop 'Bad Lubrication'
    ENDIF
    
    RETURN
    END
