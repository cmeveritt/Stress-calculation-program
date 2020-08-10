! Subroutine to calculate the influence of point and line loads on stresses at a point ik, jl distance away. 
    !SUBROUTINE term_calc
    !    implicit none
    !    COMMON
    !    COMMON
    !    real*8 term_xx_l_s(-601:601), term_yy_l_s(-601:601), term_zz_l_s(-601:601), term_xy_l_s(-601:601), term_xz_l_s(-601:601), term_zy_l_s(-601:601)
    !    real*8 term_xx_p_s(-601:601,-601:601), term_yy_p_s(-601:601,-601:601), term_zz_p_s(-601:601,-601:601), term_xy_p_s(-601:601,-601:601), term_xz_p_s(-601:601,-601:601), term_zy_p_s(-601:601,-601:601)
    !    real*8 term_xx_l_d(-601:601,0:301), term_yy_l_d(-601:601,0:301), term_zz_l_d(-601:601,0:301), term_xy_l_d(-601:601,0:301), term_xz_l_d(-601:601,0:301), term_zy_l_d(-601:601,0:301)
    !    real*8 term_xx_p_d(-301:301,-301:301,0:301), term_yy_p_d, term_zz_p_d, term_xy_p_d, term_xz_p_d, term_zy_p_d