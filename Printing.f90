subroutine Printing
        integer      NN, I, J
        include     'inc_Grid.h'
        include     'inc_StressParam.h'
        include     'inc_Stress_surf.h'
        include     'inc_Stress_under.h'
        include     'inc_Prin_Stress_surf.h'
        include     'inc_Prin_Stress_under.h'

        
        NN=(NY+1)/2
        
        ! Mirroring the surface stresses. 
        DO J=1,NN
            JJ=NY-J+1
            DO I=1,NX
                !sig_xs(I,JJ)=sig_xs(i,j)
                !sig_ys(I,JJ)=sig_ys(i,J)
                !sig_zs(I,JJ)=sig_zs(i,j)
                !tau_xys(I,JJ)=tau_xys(i,j)
                !tau_xzs(I,JJ)=tau_xzs(i,j)
                !tau_yzs(I,JJ)=tau_yzs(i,j)
                
                !sig_1s(I,JJ)=sig_1s(i,j)
                !sig_2s(I,JJ)=sig_2s(i,j)
                !sig_3s(I,JJ)=sig_3s(i,j)
                
                !sig_vms(I,JJ)=sig_vms(i,j)
            ENDDO
        ENDDO
            
        !DO J=1,Nz
        !    WRITE(13,110)(sig_1(I,NZ-J+1),I=1,Nx)      !Depth
        !ENDDO
        
        DO i=1,Nx
            sig_1s(i,1)=sig_1s(i,2)
            WRITE(14,110)(sig_1s(I,J),j=1,NN)     !Surface
        ENDDO
        
        DO i=1,Nx
            sig_vms(i,1)=sig_vms(i,2)
            WRITE(15,110)(sig_vms(I,J),j=1,NN)     !surface
        ENDDO
        
        !DO J=1,Nz
        !    WRITE(16,110)(sig_vm(I,NZ-J+1),I=1,Nx)      !Depth
        !ENDDO
        
        DO i=1,Nx
            WRITE(8,110)(sig_xs(I,J),j=1,NN)     !Surface
        ENDDO
        !
        !DO i=1,Nx
        !    WRITE(9,110)(sig_ys(I,J),j=1,Ny)     !Surface
        !ENDDO
        !        
        !DO i=1,Nx
        !    WRITE(10,110)(sig_zs(I,J),j=1,Ny)     !Surface
        !ENDDO
        !
        !DO i=1,Nx
        !    WRITE(30,110)(tau_xys(I,J),j=1,Ny)     !Surface
        !ENDDO
        !
        DO i=1,Nx
            tau_xzs(i,1)=tau_xzs(i,2)
            WRITE(31,110)(tau_xzs(I,J),j=1,NN)     !Surface
        ENDDO
        !        
        DO i=1,Nx
            tau_yzs(i,1)=tau_yzs(i,2)
            WRITE(32,110)(tau_yzs(I,J),j=1,NN)     !Surface
        ENDDO
        !
        !!!----------Depth stress state--------------------
        !
        !DO J=1,Nz
        !    WRITE(41,110)(sig_x(I,NZ-J+1),I=1,Nx)      !Depth
        !ENDDO
        !
        !DO J=1,Nz
        !    WRITE(42,110)(sig_y(I,NZ-J+1),I=1,Nx)      !Depth
        !ENDDO
        !        
        !DO J=1,Nz
        !    WRITE(43,110)(sig_z(I,NZ-J+1),I=1,Nx)      !Depth
        !ENDDO
        !                
        !                        
        !DO J=1,Nz
         !   WRITE(44,110)(tau_xz(I,NZ-J+1),I=1,Nx)      !Depth
        !ENDDO
        !                                
        !DO J=1,Nz
        !    WRITE(45,110)(tau_yz(I,NZ-J+1),I=1,Nx)      !Depth
        !ENDDO
        !        
        
    110 FORMAT(2001(E12.6,1X))
        RETURN
    END