! Subrutine to calculate the stresses  in the x-z plane at the centerline of the contact
    SUBROUTINE StressP_under
        implicit none
        integer      NN, NNX
        integer      i, j, k, l
        integer      ik, jl, js
        real*8         sigxx, sigzz
        real*8         rho3, rho5, base
        real*8         min, pi
        real*8         rho, r2, z
        integer      hh
        include     'inc_Grid.h'
        include     'inc_Nodes_no_pressure.h'
        include     'inc_StressParam.h'
        
    pi=3.14159265
    NN=(NY+1)/2
    NNX=(NX+1)/2
            
    ! Calculats the stresses from the pressure at the area         
   !Near under
    
        DO hh=1,Nr_z*NZ_cvot+1                         ! For each height below the surface
        z=(hh-1.0)/NZ_cvot                     ! To have same scale on z as ik and jl 
        DO i=2,NX-1
            j=NN                        ! Only interested in the centerline
            Call StressP_Under_line(hh,i,j)
            
        enddo
    enddo
                
    ! Combined with the stresses from pointloads----------------------------------------------------------------------     
            
    DO hh=1,Nr_z*NZ_cvot+1                        ! For each height below the surface

        DO i=2,NX-1
            j=NN     
            DO jl=-NNX,-1
                    CALL StressP_near_Under_point(hh,i,j,jl)
                    Call StressP_near_Under_point(hh,i,j,-jl)
            
            enddo
            CALL StressP_near_Under_point(hh,i,j,0)
            
        ENDDO
    ENDDO
    
    
    !Far under
    
    DO hh=Nr_z*NZ_cvot+2,Nz                          ! For each height below the surface
        z=(hh-Nr_z*(NZ_cvot-1)-1)                    ! To have same scale on z as ik and jl 
        DO i=2,NX-1
            j=NN                        ! Only interested in the centerline
            ! Line load -----------------------------------------------------------------------------------------------------
            Call StressP_Under_line(hh,i,j)
        enddo
    enddo
                
    ! Combined with the stresses from pointloads----------------------------------------------------------------------     
            
    DO hh=Nr_z*NZ_cvot+2,Nz                          ! For each height below the surface
        z=(hh-Nr_z*(NZ_cvot-1)-1)                    ! To have same scale on z as ik and jl 

        DO i=2,NX-1
            j=NN 
            DO jl=-NNX,-1                          ! To sum symmetric to mimimize nummerical summation errors
                    CALL StressP_near_Under_point(hh,i,j,jl)
                    Call StressP_near_Under_point(hh,i,j,-jl)
            enddo
            CALL StressP_near_Under_point(hh,i,j,0)
                
        enddo
    enddo

    
    RETURN
    
end
    