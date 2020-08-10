Subroutine PrintMax
    !Soubroutine for printing the maximum values when calculated
    implicit none
    integer     NN,I,J, JJ
    include     'inc_Grid.h' 
    include     'inc_StressParam.h'
    include     'inc_Max_val_asp.h'
    include     'inc_Max_val_cont.h'
    include     'inc_tau_max_values.h'
    include     'inc_directions.h'


    
    NN=(NY+1)/2
        
        
    ! Mirroring the surface stresses. 
        DO J=1,NN
            JJ=NY-J+1
            DO I=1,NX
                sig_1ms(I,JJ)  =sig_1ms(i,j)
                sig_vmms(I,JJ) =sig_vmms(i,J)
                tau_maxs(I,JJ) =tau_maxs(i,j)
                
                sig_xms(I,JJ)  =sig_xms(i,j)
                sig_yms(I,JJ)  =sig_yms(i,J)
                sig_zms(I,JJ)  =sig_zms(i,j)
            ENDDO
                        
            DO I=1,NX*2
                sig_1msa(I,JJ)  =sig_1msa(i,j)
                sig_vmmsa(I,JJ) =sig_vmmsa(i,J)
                P_asp(I,JJ)     =P_asp(i,j)
                tau_max_as(I,JJ)=tau_max_as(i,j)
                
                sig_xmsa(I,JJ) =sig_xmsa(i,j)
                sig_ymsa(I,JJ) =sig_ymsa(i,J)
                sig_zmsa(I,JJ) =sig_zmsa(i,j)
            ENDDO
        ENDDO
         
        ! Principal stresses and pressure-----------------------------------
        OPEN(21,FILE='Sig_1surf_max.DAT',STATUS='UNKNOWN')      ! To Write
        OPEN(22,FILE='Sig_1depth_max.DAT',STATUS='UNKNOWN')      ! To Write
        
        OPEN(23,FILE='Sig_vm_surf_max.DAT',STATUS='UNKNOWN')
        OPEN(24,FILE='Sig_vm_depth_max.DAT',STATUS='UNKNOWN')
        
        OPEN(25,FILE='Sig_1surf_max_asp.DAT',STATUS='UNKNOWN')      ! To Write
        OPEN(26,FILE='Sig_1depth_max_asp.DAT',STATUS='UNKNOWN')      ! To Write
        
        OPEN(27,FILE='Sig_vm_surf_max_asp.DAT',STATUS='UNKNOWN')
        OPEN(28,FILE='Sig_vm_depth_max_asp.DAT',STATUS='UNKNOWN')
        
        OPEN(29,FILE='PRESSURE_max_asp.DAT',STATUS='UNKNOWN')
        
        OPEN(45,FILE='Tau_maxs.DAT',STATUS='UNKNOWN')           ! Surface
        OPEN(46,FILE='Tau_max_asps.DAT',STATUS='UNKNOWN') 
        
        OPEN(47,FILE='Tau_max.DAT',STATUS='UNKNOWN')            ! Down under
        OPEN(48,FILE='Tau_max_asp.DAT',STATUS='UNKNOWN')
        
        OPEN(50,FILE='Direction_x_s.DAT',STATUS='UNKNOWN')
        OPEN(52,FILE='Direction_y_s.DAT',STATUS='UNKNOWN')
        OPEN(54,FILE='Direction_z_s.DAT',STATUS='UNKNOWN')
        
        OPEN(51,FILE='Direction_x_d.DAT',STATUS='UNKNOWN')
        OPEN(53,FILE='Direction_y_d.DAT',STATUS='UNKNOWN')
        OPEN(55,FILE='Direction_z_d.DAT',STATUS='UNKNOWN')
        
        ! In xyz- coordinates--------------------------------------------------
        OPEN(30,FILE='Sig_xsurf_max.DAT',STATUS='UNKNOWN')      ! To Write
        OPEN(31,FILE='Sig_ysurf_max.DAT',STATUS='UNKNOWN')      ! To Write
        OPEN(32,FILE='Sig_zsurf_max.DAT',STATUS='UNKNOWN')      ! To Write
        
        OPEN(33,FILE='Sig_x_depth_max.DAT',STATUS='UNKNOWN')
        OPEN(34,FILE='Sig_y_depth_max.DAT',STATUS='UNKNOWN')
        OPEN(35,FILE='Sig_z_depth_max.DAT',STATUS='UNKNOWN')
        
        OPEN(36,FILE='Sig_xsurf_max_asp.DAT',STATUS='UNKNOWN')      ! To Write
        OPEN(37,FILE='Sig_ysurf_max_asp.DAT',STATUS='UNKNOWN')      ! To Write
        OPEN(38,FILE='Sig_zsurf_max_asp.DAT',STATUS='UNKNOWN')      ! To Write
        
        OPEN(39,FILE='Sig_x_depth_max_asp.DAT',STATUS='UNKNOWN')
        OPEN(40,FILE='Sig_y_depth_max_asp.DAT',STATUS='UNKNOWN')
        OPEN(41,FILE='Sig_z_depth_max_asp.DAT',STATUS='UNKNOWN')
        
        ! Contact viev---------------------------------------------------
        DO i=1,Nx
            WRITE(21,110)(sig_1ms(I,J),j=1,Ny)     !Surface
        ENDDO     
        DO J=1,Nz
            WRITE(22,110)(sig_1m(I,NZ-J+1),I=1,Nx)      !Depth
        ENDDO
        
        DO i=1,Nx
            WRITE(23,110)(sig_vmms(I,J),j=1,Ny)     !surface
        ENDDO
        DO J=1,Nz
            WRITE(24,110)(sig_vmm(I,NZ-J+1),I=1,Nx)      !Depth
        ENDDO
        
        DO i=1,Nx
            WRITE(45,110)(tau_maxs(I,J),j=1,Ny)     !surface
        ENDDO
        DO J=1,Nz
            WRITE(47,110)(tau_max(I,NZ-J+1),I=1,Nx)      !Depth
        ENDDO
        
        ! Asperity viev -------------------------------------------------
        DO i=1,Nx*2
            WRITE(25,110)(sig_1msa(I,J),j=1,Ny)             !Surface
        ENDDO     
        DO J=1,Nz
            WRITE(26,110)(sig_1ma(I,NZ-J+1),I=1,Nx*2)      !Depth
        ENDDO
        
        DO i=1,Nx*2
            WRITE(27,110)(sig_vmmsa(I,J),j=1,Ny)     !surface
        ENDDO
        DO J=1,Nz
            WRITE(28,110)(sig_vmma(I,NZ-J+1),I=1,Nx*2)      !Depth
        ENDDO
        
        DO i=1,Nx*2
            WRITE(29,110)(P_asp(I,J),j=1,Ny)      !surface
        ENDDO
        
        DO i=1,Nx*2
            WRITE(46,110)(tau_max_as(I,J),j=1,Ny)     !surface
        ENDDO
                
        
        DO J=1,Nz
            WRITE(48,110)(tau_max_a(I,NZ-J+1),I=1,Nx*2)      !Depth
        ENDDO
        
        !! Directions ---------------------------------------------------------
        DO i=1,Nx*2
            WRITE(50,110)(Nx_s(I,J),j=1,NN)             !surface
        ENDDO
                
        DO J=1,Nz
            WRITE(51,110)(Nx_d(I,NZ-J+1),I=1,Nx*2)      !Depth
        ENDDO
        
        DO i=1,Nx*2
            WRITE(52,110)(Ny_s(I,J),j=1,NN)             !surface
        ENDDO
                
        DO J=1,Nz
            WRITE(53,110)(Ny_d(I,NZ-J+1),I=1,Nx*2)      !Depth
        ENDDO
        
        DO i=1,Nx*2
            WRITE(54,110)(Nz_s(I,J),j=1,NN)             !surface
        ENDDO
                
        DO J=1,Nz
            WRITE(55,110)(Nz_d(I,NZ-J+1),I=1,Nx*2)      !Depth
        ENDDO
        
        
       !xyz coordinate system ----------------------------------------
        !Contact view ------------------------------------------------
        DO i=1,Nx
            WRITE(30,110)(sig_xms(I,J),j=1,Ny)             !Surface
        ENDDO     
        DO i=1,Nx
            WRITE(31,110)(sig_yms(I,J),j=1,Ny)             !Surface
        ENDDO   
        DO i=1,Nx
            WRITE(32,110)(sig_zms(I,J),j=1,Ny)             !Surface
        ENDDO     
                                      
        DO J=1,Nz
            WRITE(33,110)(sig_xm(I,NZ-J+1),I=1,Nx)      !Depth
        ENDDO
        DO J=1,Nz
            WRITE(34,110)(sig_ym(I,NZ-J+1),I=1,Nx)      !Depth
        ENDDO
        DO J=1,Nz
            WRITE(35,110)(sig_zm(I,NZ-J+1),I=1,Nx)      !Depth
        ENDDO
        
        !Asperity view ------------------------------------------------
        DO i=1,Nx*2
            WRITE(36,110)(sig_xmsa(I,J),j=1,Ny)             !Surface
        ENDDO     
        DO i=1,Nx*2
            WRITE(37,110)(sig_ymsa(I,J),j=1,Ny)             !Surface
        ENDDO     
        DO i=1,Nx*2
            WRITE(38,110)(sig_zmsa(I,J),j=1,Ny)             !Surface
        ENDDO     
                                
                                
        DO J=1,Nz
            WRITE(39,110)(sig_xma(I,NZ-J+1),I=1,Nx*2)      !Depth
        ENDDO
        DO J=1,Nz
            WRITE(40,110)(sig_yma(I,NZ-J+1),I=1,Nx*2)      !Depth
        ENDDO
        DO J=1,Nz
            WRITE(41,110)(sig_zma(I,NZ-J+1),I=1,Nx*2)      !Depth
        ENDDO
                
        
        
        
        
        
    110 FORMAT(2001(E12.6,1X))
        RETURN
    END