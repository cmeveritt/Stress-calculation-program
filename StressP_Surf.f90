! Subroutine to fint the stresses at the surface due to a normal load
    SUBROUTINE StressP_Surf
        implicit none
        integer      NN, NNX
        integer      i, j, jl
        real        pi
        include     'inc_Grid.h'


        
    pi=3.14159265
    NN=(NY+1)/2   
    NNX=(NX+1)/2            !Its enough distance to get converged stress values. Moste of the stress levels are due to the line loads
    ! Calculats the stresses from the pressure at the area
    DO i=2,NX-1
        DO j=2,NN               ! Interested in the whole area   
                                !  from the pressure at the area
            
        ! The stresses from a pure line load------------------------------------------
        ! At the surface a line load, Pmin, will not couse any stresses exept than where its applied. Point loads will on the other hand.  
            ! To mimimize nummerical summation error
            DO jl=-NNX,-1
                Call StressP_Surf_point(i,j,jl)
                Call StressP_Surf_point(i,j,-jl)
            ENDDO
                Call StressP_Surf_point(i,j,0)       
            
        enddo
    enddo
    
    
    
        
    return
    
end
    
        
