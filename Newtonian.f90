! Subrutine for updating the film thicknes, the viscocity and the density based on the given pressure
    ! Almost the same as in the EHL calculations but does not change the temperature, nor call the stop to large out
	SUBROUTINE Newtonian(SS, NYs, T40, Tcvot, NN, t)
        implicit none 
        COMMON  /Current/       RO,EPSx,EPSy,EDAx,EDAy,xi,W,Wside                      ! Current timestep
        COMMON  /CurrentP/      P
        COMMON  /CurrentH/      H
        COMMON  /Contact_mat/   contact
        COMMON  /Visc/          ENDA,A1,A2,A3,Z,HM0r, PH, Pref, alpha, EDA0         ! Lubrication parameters
        COMMON  /Rho/           RA1,RA2                                             ! Density parameters
        COMMON  /Grid/          NX,NY,X0,XE,DX                                      ! Grid parameters
        COMMON  /Holmes/        kH, xH, gH, lH                                   ! Viscosity and denisty param acc Holems et al
        COMMON  /Ref/           PAIAK, meth, tmeth, Geom, Lub_param, asp_shape, contact_alg      ! Coise of equations based on referense
        COMMON  /CurrentT/      Temp
        COMMON  /RLarsson/      EpsT0,RL_Ta,S0, RL_G0, Dz, Cz,RL_T0, RL_c                              ! Parameters for ref 10, R. Larssons formulation
        COMMON  /outp/          W0,EE,RX,Um, Ua, B, U2_over_Us, BY, Ry 
        COMMON /Yasutomi/       Ta,Tg0,YA1,YA2,YB1,YB2,YC1,YC2,Yedag,Yalfag         ! Lubrication parameters acc Yasutomi
        COMMON  /Y_Liu/         L_n, L_G, L_h_limit, L_iter, L_stab                 ! Lubrication and convergence parameters acc Y.Liu 2007 
        COMMON  /Method/        term1,term2,term3,term4                             ! Controling the numerical method
        COMMON /shear_lim/      shear_max, shear_min, temp_max
        ! Input
        real        P(1:601,1:601),H(1:601,1:601), Temp(1:601,1:601)
        integer     contact(1:601,1:601)
        real        RX, ENDA, Um, Ua
        integer     NX, NN, SS, Lub_param, t
        Real        Z,  EDA0, Pref, alpha
        real        B, PH
        real        kH, xH, gH, lH
        real        EpsT0,RL_Ta,S0, RL_G0, Dz, Cz,RL_T0, RL_c, EpsT
        real        T40, Tcvot, L_stab
        real        Ta,Tg0,YA1,YA2,YB1,YB2,YC1,YC2,Yedag,Yalfag  
        integer     term1,term2,term3,term4
        real        shear_max, shear_min, temp_max
        ! Calculations
        integer     I,J, JJ, temp_iter
        real        EDA_1, EDA_2, EDA_3
        real        EDA_cont, EDA22, EDA33, EDA1
        real        Emat, pg, Tg, YF 
        real        temp_0, ub, u_slip, h_real
        real        tau_x_f, tau_y_f, tau_x_c, tau_y_c
        real        dPdx, dPdy
        integer     J0, J1, J2, I0, I1, I2
        ! Other
        real        W(1:601,1:601),Wside(1:601,1:601)                      ! Current timestep
        Integer     NY, tmeth, NYs, Geom, meth, contact_alg
        Integer     asp_shape 
        Real        X0, XE, PAIAK
        Real        W0, RA1, RA2, BY, Ry, EE, U2_over_Us
        Real        DX
        real        HM0r
        real        A1, A2, A3
        real        L_n, L_G, L_h_limit, L_iter 
        ! Output
        real        RO(1:601,1:601),EPSx(1:601,1:601),EPSy(1:601,1:601),EDAx(1:601,1:601),EDAy(1:601,1:601),xi(1:601,1:601)
        SAVE      /Current/                    ! Current timestep

            
        !Minimum viscosity if contact
        If( lub_param .EQ. 1)THEN    
            EDA_cont=EXP(alpha*Pref/Z*(-1+(1+3.0*PH/Pref)**Z))
        ELSE If( lub_param .EQ. 11)THEN    
            EDA_cont=EXP(alpha*(3.0*PH)**Z)
        else If( lub_param .EQ. 2)THEN     
                ! Roelands equation acc Gohar
                EDA_1=5.1*PH*3.0/10**(9)
                EDA_2=(1.0+EDA_1)**Z
                EDA_3=(log(EDA0)+9.67)
                EDA_cont=EXP(EDA_3*(EDA_2-1.0))
        else If( lub_param .EQ. 3)THEN        
            EDA_cont=EXP(LOG(EDA0/kH)*((1+xH*PH*3.0)**Z-1)) 
        else If( lub_param .EQ. 4 .OR. lub_param .EQ. 60)THEN
            EDA_cont=EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*3.0*PH)**Z))  !Roland Larsson pressure visc relation
        elseif( Lub_param .eq.5) then
            EDA_cont=Yedag
        else If( lub_param .EQ. 6)THEN
            EDA_cont=EXP(log(EDA0+9.67)*((1+3.0*PH/Pref)**Z-1)) 
        !else
            !EDA_cont=Yedag
        ENDIF
        
        
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
                
                
                !IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                !    WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'P = ', P(I,J),  'EDAO = ',EDA0
                !    EDAx(I,J)=0.1
                !    Call Stop_to_large_out
                !ENDIF
                !
                !IF ( RO(I,J) .LT. 1.0 .OR. isnan(RO(I,J))) THEN
                !    WRITE(4,*)'BAD RO', RO(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT
                !    RO(I,J)=1.0
                !    Call Stop_to_large_out
                !ENDIF
        
                IF(contact(I,J) .EQ. 1 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
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
                
                
                !
                !IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                !    WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0
                !    EDAx(I,J)=0.1
                !    Call Stop_to_large_out
                !ENDIF
                !
                !IF ( RO(I,J) .LT. 1.0 .OR. isnan(RO(I,J))) THEN
                !    WRITE(4,*)'BAD RO', RO(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT
                !    RO(I,J)=1.0
                !    Call Stop_to_large_out
                !ENDIF
                !
                IF(contact(I,J) .EQ. 1 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
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
                
                !
                !IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                !    WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0
                !    EDAx(I,J)=0.1
                !    Call Stop_to_large_out
                !ENDIF
                !
                !IF ( RO(I,J) .LT. 1.0 .OR. isnan(RO(I,J))) THEN
                !    WRITE(4,*)'BAD RO', RO(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT
                !    RO(I,J)=1.0
                !    Call Stop_to_large_out
                !ENDIF
        
                IF(contact(I,J) .EQ. 1 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
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
                
                !
                !IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                !    WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0
                !    EDAx(I,J)=0.1
                !    Call Stop_to_large_out
                !ENDIF
                !
                !IF ( RO(I,J) .LT. 1.0 .OR. isnan(RO(I,J))) THEN
                !    WRITE(4,*)'BAD RO', RO(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT
                !    RO(I,J)=1.0
                !    Call Stop_to_large_out
                !ENDIF
                !
                IF(contact(I,J) .EQ. 1 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                EPSy(I,J)=EPSx(I,J)!RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))
        ENDDO
    ENDDO
        
    ELSE IF( Lub_param .EQ. 4 )THEN      ! No temperature adjustments !This one is changed from the EHL vesion
        ! R.Larssons 2000 formulation
        DO J=1,NN,SS
            DO I=1,NX,SS
                RL_Ta=temp(I,J)
                ! Roelands equation acc R.Larsson
                EDA0 = 10**(-4.2+RL_G0*(1.0+RL_Ta/135.0)**S0)             ! Eq (3)
                Z=Dz+Cz*log10(1.0+RL_Ta/135)
                ENDA=12  *Um * RX**2 / ( B**3*PH) ! *EDA0 but EDA0 is included in EPSx

                EDAx(I,J)=EDA0*EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
                EDAy(I,J)=EDAx(I,J)!EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
                
                ! D-H Formulation acc R.Larsson
                EpsT       = EpsT0*exp(-RL_c*P(I,J)*PH)
                RO(I,J)    = (1.0+RA1*PH*P(I,J)/(1.0+RA2*PH*P(I,J)))*(1.0-EpsT*(RL_Ta-RL_T0)) 
                xi(I,J)=0.0
                
                !
                !IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                !    WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0
                !    EDAx(I,J)=0.1
                !    Call Stop_to_large_out
                !ENDIF
                !
                !IF ( RO(I,J) .LT. 0.0 .OR. isnan(RO(I,J))) THEN
                !    WRITE(4,*)'BAD RO', RO(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT
                !    RO(I,J)=1.0
                !    Call Stop_to_large_out
                !ENDIF
        
                IF(contact(I,J) .EQ. 1 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                EPSy(I,J)=EPSx(I,J)!RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))
            ENDDO
        ENDDO
    
!    ELSE IF(lub_param .EQ. 4 .and. l_stab .NE. 0) then ! Themperature rise
!        
!        ub= 2*um-ua
!        u_slip=ua-ub
!        
!        DO J=1,NN,SS
!            DO I=1,NX,SS
!                temp_iter=0             ! Resetting counter
!                temp_0=temp(1,1)                        ! The global temperature
!                
!                ! Ensure that the lubrication can not cool off downstream. 
!777             IF( L_stab .EQ. 2 .AND. t .LE. -2 .and. I .GT. 1 ) Then
!                    temp(I,j)=max(temp(i,j),temp(I-SS,J))       ! Do not update if higher temp 
!                    temp_0=temp(I-SS,J)       
!                ENDIF
!            
!                RL_Ta=temp(I,J)         ! Extract the temperature
!                
!                ! Roelands equation acc R.Larsson
!                EDA0 = 10**(-4.2+RL_G0*(1+RL_Ta/135)**S0)             ! Eq (3)
!                ENDA=12  *Um * RX**2 / ( B**3*PH) ! *EDA0 but EDA0 is included in EPSx
!
!                EDAx(I,J)=EDA0*EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
!                EDAy(I,J)=EDAx(I,J)!EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
!                
!                ! D-H Formulation acc R.Larsson
!                EpsT=EpsT0*exp(-RL_c*P(I,J)*PH)
!                RO(I,J)=(1+RA1*PH*P(I,J)/(1+RA2*PH*P(I,J)))*(1-EpsT*(RL_Ta-RL_T0)) 
!                xi(I,J)=0.0
!                
!                
!                IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
!                    WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0
!                    EDAx(I,J)=0.1
!                    Call Stop_to_large_out
!                ENDIF
!
!                IF ( RO(I,J) .LT. 0.0 .OR. isnan(RO(I,J))) THEN
!                    WRITE(4,*)'BAD RO', RO(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT
!                    RO(I,J)=1.0
!                    Call Stop_to_large_out
!                ENDIF
!        
!                IF(contact(I,J) .EQ. 1 .and. EDAx(I,J) .LT. EDA_cont) then
!                    EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
!                    
!                ELSE IF( temp_iter .LE. 10 .AND. SS .LE. 2) THEN     ! Might be good to introduce this after the coursest solution is obtained. 
!                    
!                    ! Define nodenumbers for past and next nodes
!                    J0=J-SS
!                    IF( J0 .LE. 0)  J0=J+SS
!                    J2=J-2*SS
!                    IF( J2 .LE. 0)  J2=J0
!                    J1=J+1*SS
!                    IF( J1 .GE. NYs) J1=NYs-SS
!                    JJ=NYs+1-J
!                    IF( JJ .LE. 1)  JJ=1+SS
!                
!                    I0=I-1*SS
!                    IF( I0 .LE. 0)  I0=1
!                    I2=I-2*SS
!                    IF( I2 .LE. 0)  I2=I0
!                    I1=I+1*SS
!                    IF( I1 .GE. NX) I1=NX
!                
!                    dPdx=(term4*P(I1,J)+term3*P(I,J)+term2*P(I0,J)+term1*P(I2,J))/(2*DX*SS)*Ph/b
!                    dPdy=(term4*P(I,J1)+term3*P(I,J)+term2*P(I,J0)+term1*P(I,J2))/(2*DX*SS)*Ph/b
!                
!                    h_real=H(I,J)*b**2/Rx
!                    tau_x_f = abs(-h_real/2*dPdx+EDAX(I,j)*u_slip/h_real)
!                    tau_y_f = abs(-h_real/2*dPdy)
!                    tau_x_c = abs( h_real/2*dPdx+EDAX(I,j)*u_slip/h_real)
!                    tau_y_c = abs( h_real/2*dPdy)
!                    
!                    
!                    if( (tau_x_f .GT. shear_max .or. tau_y_f .gt. shear_max .or. tau_x_c .gt. shear_max .OR. tau_y_c .GT. shear_max) .AND. temp(i,j) .LT. temp_max) then
!                        
!                        temp(I,J)=temp(I,J)+5
!                        RL_Ta=temp(I,J)
!                        ! Roelands equation acc R.Larsson
!                        EDA0 = 10**(-4.2+RL_G0*(1+RL_Ta/135)**S0)             ! Eq (3)
!                        ENDA=12  *Um * RX**2 / ( B**3*PH) ! *EDA0 but EDA0 is included in EPSx
!
!                        EDAx(I,J)=EDA0*EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
!                        EDAy(I,J)=EDAx(I,J)!EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
!                
!                        ! D-H Formulation acc R.Larsson
!                        EpsT=EpsT0*exp(-RL_c*P(I,J)*PH)
!                        RO(I,J)=(1+RA1*PH*P(I,J)/(1+RA2*PH*P(I,J)))*(1-EpsT*(RL_Ta-RL_T0)) 
!                
!                
!                        IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
!                            WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0
!                            EDAx(I,J)=0.1
!                            Call Stop_to_large_out
!                        ENDIF
!
!                        IF ( RO(I,J) .LT. 0.0 .OR. isnan(RO(I,J))) THEN
!                            WRITE(4,*)'BAD RO', RO(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT
!                            RO(I,J)=1.0
!                            Call Stop_to_large_out
!                        ENDIF
!        
!                        IF(contact(I,J) .EQ. 1 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
!                        EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
!                        EPSy(I,J)=EPSx(I,J)!RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))
!                        
!                        
!                        
!                        temp_iter=temp_iter+1
!                        GO TO 777
!                    elseif( (tau_x_f .LT. shear_min .AND. tau_y_f .LT. shear_min .AND. tau_x_c .LT. shear_min .AND. tau_y_c .LT. shear_min) .AND. temp(i,j) .GT. temp_0) then
!                        
!                        temp(I,J)=temp(I,J)-2
!                        RL_Ta=temp(I,J)
!                        ! Roelands equation acc R.Larsson
!                        EDA0 = 10**(-4.2+RL_G0*(1+RL_Ta/135)**S0)             ! Eq (3)
!                        ENDA=12  *Um * RX**2 / ( B**3*PH) ! *EDA0 but EDA0 is included in EPSx
!
!                        EDAx(I,J)=EDA0*EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
!                        EDAy(I,J)=EDAx(I,J)!EXP((LOG(EDA0)+9.67)*(-1.0+(1.0+5.1E-9*P(I,J)*PH)**Z)) 
!                
!                        ! D-H Formulation acc R.Larsson
!                        EpsT=EpsT0*exp(-RL_c*P(I,J)*PH)
!                        RO(I,J)=(1+RA1*PH*P(I,J)/(1+RA2*PH*P(I,J)))*(1-EpsT*(RL_Ta-RL_T0)) 
!                
!                
!                        IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
!                            WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0
!                            EDAx(I,J)=0.1
!                            Call Stop_to_large_out
!                        ENDIF
!
!                        IF ( RO(I,J) .LT. 0.0 .OR. isnan(RO(I,J))) THEN
!                            WRITE(4,*)'BAD RO', RO(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT
!                            RO(I,J)=1.0
!                            Call Stop_to_large_out
!                        ENDIF
!        
!                        IF(contact(I,J) .EQ. 1 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
!                        EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
!                        EPSy(I,J)=EPSx(I,J)!RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))
!
!                        
!                        temp_iter=temp_iter+1
!                        GO TO 777
!                    endif
!                ENDIF 
!                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
!                EPSy(I,J)=EPSx(I,J)!RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))
!
!
!            ENDDO
!        ENDDO
!        
!        !Mirroring the temperature
!        DO J=1,NN,SS
!            JJ=NYs-J+1
!            DO I=1,NX,SS   
!                temp(I,JJ)=temp(I,J)
!            ENDDO
        !ENDDO         
                    
    ELSE IF( Lub_param .EQ. 5) THEN
        ! Yasutomi lubrication
        DO J=1,NN,SS
            DO I=1,NX,SS

                ! Yasutomi equations
                Tg=Tg0+YA1*LOG(1+YA2*PH*P(I,J))
                pg=1/YA2*(EXP((1/YA1)*((T40*(1-Tcvot)+Tcvot*Temp(I,J))-Tg0))-1)
                YF=1-YB1*LOG(1+YB2*PH*P(I,J))
                
                ! Since for some specific combinations of p and YF EDA33 gest biger than EDA22 eaven do p < pg this formulation should be better. 
                IF (P(I,J)*PH .GT. pg) THEN
                    EDA22=Yedag*EXP(Yalfag*(P(I,J)*PH-pg))
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
                
                ! IF (EDAx(I,J) .LE. 0.0) THEN
                !    WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J ,'EDA1 = ', EDA1
                !    EDAx(I,J)=0.1
                !    Call Stop_to_large_out
                !ENDIF
                !
                
                ! D-H Formulation acc P.Ehret D.Dowsin and C.M. Taylor
                RO(I,J)=1+RA1*PH*P(I,J)/(1+RA2*PH*P(I,J)) 
                
                IF(contact(I,J) .EQ. 1 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                EPSy(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                
               
                
            ENDDO
        ENDDO
    ELSE IF( Lub_param .EQ. 6) THEN
        ! Newtonian acc N Deolalikers contact paper from 2008
        DO J=1,NN,SS
            DO I=1,NX,SS

                ! Roelands equation acc X.Tan and Venner
                EDAx(I,J)=EXP(log(EDA0+9.67)*((1+P(I,J)*PH/Pref)**Z-1)) 
                if (EDAx(I,J) .GT. 1E32) EDAx(I,J)=1E32                 ! To limit the viscosity to som reasonable? values. Intended to increase convergense. 
                EDAy(I,J)=EDAx(I,J)
                
                ! D-H Formulation acc X.Tan and N Deolaliker
                ! RO(I,J)=(5.9E8+1.34*P(I,J)*PH)/(5.9E8+P(I,J)*PH) 
                 RO(I,J)=1+RA1*PH*P(I,J)/(1+RA2*PH*P(I,J)) 
                 xi(I,J)=0.0
                
                
                !IF (EDAx(I,J) .LE. 0.0 .OR. isnan(EDAx(I,J)) ) THEN
                !    WRITE(4,*)'BAD EDAx', EDAx(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0
                !    EDAx(I,J)=0.1
                !    Call Stop_to_large_out
                !ENDIF
                !
                !IF ( RO(I,J) .LT. 1.0 .OR. isnan(RO(I,J))) THEN
                !    WRITE(4,*)'BAD RO', RO(I,J), 'For I J = ', I ,J , 'P = ', P(I,J), 'H = ', H(I,J), 'EDAO = ',EDA0, 'EpsT =', EpsT
                !    RO(I,J)=1.0
                !    Call Stop_to_large_out
                !ENDIF
        
                IF(contact(I,J) .EQ. 1 .and. EDAx(I,J) .LT. EDA_cont) EDAx(I,J)=EDA_cont    !If contact ensure high viscosity
                EPSx(I,J)=RO(I,J)*H(I,J)**3/(ENDA*EDAx(I,J))
                EPSy(I,J)=EPSx(I,J)!RO(I,J)*H(I,J)**3/(ENDA*EDAy(I,J))
                
                
            ENDDO
        ENDDO
    ENDIF
     
    
    !        !Mirroring the viscosity
        DO J=1,NN,SS
            JJ=NYs-J+1
            DO I=1,NX,SS   
                EDAx(I,JJ)=EDAx(I,J)
                EDAy(I,JJ)=EDAy(I,J)
            ENDDO
        ENDDO   
    
    RETURN
    END
