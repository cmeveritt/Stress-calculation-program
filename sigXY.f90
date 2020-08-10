Subroutine sigXY
! Subrutine to calculate the stresses  in the x-z plane at the centerline of the contact
    
    CALL Loadcalc
    
    ! Calculate the stresses due to the normal pressure
    CALL StressP_Surf
    
    !CALL StressP_Under
    
    ! Calculate the stresses due to shear loads. They read the stresses from the normal pressure and att the stresses from the shear loads
    CALL Shear_Stress_Surf_tauXZ
    CALL Shear_Stress_Surf_tauYZ
    
    
    !CALL Shear_Stress_under_tauXZ
    
    ! Due to symmetry this two will not contribute
    ! CALL Shear_Stress_under_tauYZ
    
    !CALL Surface_at_top_depth
    !CALL Smothening
    
    CALL vonMises_surf
    
    !CALL vonMises_under
    
    
    
    return
    
end
    