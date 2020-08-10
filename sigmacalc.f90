Subroutine sigmacalc
    !Subroutine for calculating the stresses
    
    Call sigXY
    
    Call Findley_print
    
    Call Eigenvalsurf
    !Call Eigenvalunder
    
    Call Printing
    
    RETURN
END