Subroutine File_definitions
        
!------ Output data files-------------------------------------------------------------------------------
        OPEN(14,FILE='Sig_1surf.DAT',STATUS='UNKNOWN')      ! To Write
        OPEN(13,FILE='Sig_1depth.DAT',STATUS='UNKNOWN')     ! To Write
        OPEN(15,FILE='Sig_vm_surf.DAT',STATUS='UNKNOWN')
        OPEN(16,FILE='Sig_vm_depth.DAT',STATUS='UNKNOWN')
        OPEN(4,FILE='Stress_OUT.DAT',STATUS='UNKNOWN')
        
        OPEN(8,FILE='Sig_xsurf.DAT',STATUS='UNKNOWN')       ! To Write
        OPEN(9,FILE='Sig_ysurf.DAT',STATUS='UNKNOWN')       ! To Write
        OPEN(10,FILE='Sig_zsurf.DAT',STATUS='UNKNOWN')
        
        ! file 20-50 are used by the maxprint
        OPEN(30,FILE='Tau_xy_surf.DAT',STATUS='UNKNOWN')
        OPEN(31,FILE='Tau_xz_surf.DAT',STATUS='UNKNOWN')
        OPEN(32,FILE='Tau_yz_surf.DAT',STATUS='UNKNOWN')
        
       ! OPEN(41,FILE='Sig_xdepth.DAT',STATUS='UNKNOWN')       ! To Write
        !OPEN(42,FILE='Sig_ydepth.DAT',STATUS='UNKNOWN')       ! To Write
        !OPEN(43,FILE='Sig_zdepth.DAT',STATUS='UNKNOWN')
        
        !!OPEN(30,FILE='Tau_xy_surf.DAT',STATUS='UNKNOWN')
        OPEN(44,FILE='Tau_xz_depth.DAT',STATUS='UNKNOWN')
        !OPEN(45,FILE='Tau_yz_depth.DAT',STATUS='UNKNOWN')
        
        
        OPEN(101,FILE='EHL_for_Findley.DAT',STATUS='UNKNOWN')      ! To Write
        OPEN(102,FILE='EHL_for_Findley_depth.DAT',STATUS='UNKNOWN')      ! To Write
                
        return 
    end
