! Subrutine to calculate the stresses  in the x-z plane due to shear loads at the centerline of the contact
    SUBROUTINE Shear_Stress_under_tauXZ
        implicit none
        integer      NN, NNX
        integer      i, j, k, l
        integer      ik, jl, js
        real*8         sigxx, sigzz, ik3, Psmall, ikjl2
        real*8         rho5, rho3, base
        real*8         rhoz2, rho3z2, rho2z3
        real*8        min, pi
        real*8         rho, r2, z
        integer      hh
        include     'inc_Grid.h'
        include     'inc_Nodes_no_pressure.h'
        include     'inc_Minimum_values.h'
        include     'inc_StressParam.h'

        
    pi=3.14159265
    NN=(NY+1)/2
    NNX=(NX+1)/2
    
    !Near under
    ! Calculats the stresses from the pressure at the area         
         ! Line load -----------------------------------------------------------------------------------------------------
    DO hh=2,Nr_z*NZ_cvot+1                         ! For each height below the surface
        z=(hh-1.0)/NZ_cvot                     ! To have same scale on z as ik and jl 
        DO i=2,NX-1
            j=NN                        ! Only interested in the centerline
      
            CALL Shear_Stress_under_XZ_line_load(hh,i,j)
            

        enddo
    enddo
                
    ! Combined with the stresses from pointloads----------------------------------------------------------------------     
            
    DO hh=2,Nr_z*NZ_cvot+1                        ! For each height below the surface
        z=(hh-1)/NZ_cvot                   ! To have same scale on z as ik and jl 

        DO i=2,NX-1
            j=NN     
            DO jl=-NNX,-1            ! This loop s not changet since I did not get any improvements by changing the loops for the pressure
                Call Shear_Stress_under_XZ_point_load(hh,i,j,jl)
                Call Shear_Stress_under_XZ_point_load(hh,i,j,-jl)
            ENDDO
                Call Shear_Stress_under_XZ_point_load(hh,i,j,0)
            
 
        ENDDO
    ENDDO
    
    ! Deep under
    ! Line load -----------------------------------------------------------------------------------------------------        
    DO hh=Nr_z*NZ_cvot+2,Nz                          ! For each height below the surface
        z=(hh-Nr_z*(NZ_cvot-1)-1)                    ! To have same scale on z as ik and jl 
        DO i=2,NX-1
            j=NN                        ! Only interested in the centerline
            
            CALL Shear_Stress_under_XZ_line_load(hh,i,j)
                    
               
        enddo
    enddo
    
 ! Combined with the stresses from pointloads----------------------------------------------------------------------     
            
    DO hh=Nr_z*NZ_cvot+2,Nz                          ! For each height below the surface
        z=(hh-Nr_z*(NZ_cvot-1)-1)                    ! To have same scale on z as ik and jl 

        DO i=2,NX-1
            j=NN     
            DO jl=-NNX,-1            ! This loop s not changet since I did not get any improvements by changing the loops for the pressure
                Call Shear_Stress_under_XZ_point_load(hh,i,j,jl)
                Call Shear_Stress_under_XZ_point_load(hh,i,j,-jl)
            ENDDO
                Call Shear_Stress_under_XZ_point_load(hh,i,j,0)
            
            
        enddo
    enddo


            
    RETURN
    
end
    