!Subroutine to fint the stresses at the surface due to a shear load
    SUBROUTINE Shear_Stress_Surf_tauYZ
        implicit none
        integer       NN, Nz, NNX
        integer      i, j, k, l
        integer         ik, jl, ikk, jll
        real*8      min, pi, r2
        real*8        tau_xxyy, rho5, rho, rho3
        include     'inc_Grid.h'

        
    pi=3.14159265
    NN=(NY+1)/2
    NNX=(NX+1)/2        
    
    ! Calculats the stresses from the pressure at the area
    
    ! The stresses from a pure line load------------------------------------------
    ! Takes no lineload into account due to symmetry 
    
    ! The stresses from point loads ----------------------------------------------------   -
    
    DO i=2,NX-1
        DO j=2,NN            ! Interested in the whole area

            DO jl=-NNX,-1
                CALL Shear_Stress_Surf_tauYZ_point(i,j,jl)
                CALL Shear_Stress_Surf_tauYZ_point(i,j,-jl)
                
            enddo
            
            CALL Shear_Stress_Surf_tauYZ_point(i,j,0)
        enddo

    enddo
        
    return
    
    end
        
