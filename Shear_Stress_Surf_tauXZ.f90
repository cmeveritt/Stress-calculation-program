!Subroutine to fint the stresses at the surface due to a shear load
    SUBROUTINE Shear_Stress_Surf_tauXZ
        implicit none
        integer      NN, NNX
        integer      i, j, k, l
        integer     ik, jl
        real*8        min, pi, r2 
        real*8        sigxx, sigyy, rho5, rho, rho3
        integer        step, ikik, iki
        real        Psmall, Psmall4
        include     'inc_Grid.h'
        
    pi=3.14159265
    NN=(NY+1)/2
    NNX=(NX+1)/2 
    
    ! Calculats the stresses from the pressure at the area
    ! The stresses from a pure line load------------------------------------------
    DO i=2,NX-1
        Call Shear_Stress_Surf_tauXZ_line(i)
        
    enddo

    ! The stresses from point loads ----------------------------------------------------   - 
    DO i=2,NX-1
        DO j=2,NN            ! Interested in the whole area

            DO jl=-NNX,-1
                CALL Shear_Stress_Surf_tauXZ_point(i,j,jl)
                CALL Shear_Stress_Surf_tauXZ_point(i,j,-jl)
                
            enddo
            
            CALL Shear_Stress_Surf_tauXZ_point(i,j,0)
        enddo

    enddo
        
    return
    
end
    
        
