	PROGRAM Principal_stresses
    ! A program to calculate the stresses based on the pressure height and temperature from TEHL simulation
    ! This program is bassed upone the line and point load solutions on infinite halfplanes. 
    ! This program has been developed with the goal of achiving a working algorithm, not a fast of beuatifull code.
    ! The structure was based upone the TEHL program so that as many variables would have the same names as possible
    ! This in order to minimize the risks of bugs. 
    
    ! The program caluclates the loads on each line and point in the load calc subroutine. 
    ! A good way to debug the program is to change the loads in this subroutine to a single line or point load, in order to validate the soultion agains the literature. 
    implicit none                                                   ! För att kräva initiering av samtliga variabler http://www.obliquity.com/computer/fortran/common.html
    Integer L1NX,   L1NY, I_NUMBER, nr_N   
    integer lub_temp, solid_material        ! Used in the TEHL code
    Real    PAI,      Y0
    Real    Elast1, Elast2
    Real    asph_real,      aspw_real
    real    asph_real2,     aspw_real2
    Real    width, HM0f,  alpha_bar
    real    Tauc_real, H00,    U,          SRR       
    real    W,      ALFA,   G,      AHM,        HM0,        UTL 
    real    sig_d,  B_ref           
    real    M,      L 
   include     'inc_Grid.h'
   include     'inc_itt.h'
   include     'inc_asp.h'
   include     'inc_outp.h'
   include     'inc_ref.h'
   include     'inc_Visc.h'    
   include     'inc_Geom5.h'
   include     'inc_Method.h'
   include     'inc_Higginson.h'
   include     'inc_Com_H00.h'
   include     'inc_G0DT.h'
   include     'inc_Pres_ave_param.h'
   include     'inc_C_method.h'
   include     'inc_DZ_com.h'

    
    
    ! Input data -----------------------------------------------------------------------------------
    ! Nummerical method and equatios -----------------------
    OPEN(4,FILE='OUT.DAT',STATUS='UNKNOWN')
    open (unit = 7, file = "Input8_1.csv")
    
    read (7,*)
    read (7,*) meth, tmeth, Geom, asp_shape, contact_alg, p_ave_param, shift_y, Multi_grid_param
    ! p_ave_param is only used in EHL calculations to obtain a more stable solution
    PAI=3.14159265                      ! Only first 7 digits that counts since single pressition
    
    IF( meth .EQ. 1) THEN               ! 1st dx
        term1   =  0
        term2   = -2
        term3   =  2
        term4   =  0
    ELSE IF( meth .EQ. 2)THEN           ! Backspace
        term1   =  1
        term2   = -4
        term3   =  3
        term4   =  0
    ELSE                                ! Central space Not working correctly
        term1   =  0
        term2   = -1
        term3   =  0
        term4   =  1
    ENDIF
    
    ! Model geometry --------------------------------------
    !Example
    !NX=150
    !NY=70                              ! NY <= NX 
    !X0=-2.0
    !XE=1.5
    !F=2.0 !0.0625                      ! Kvot between time and space steps
    !Ntime=450                          ! Number of steps in time direction. The code is based on that the asperity travels from X0 to XE
    
    read (7,*)
    read (7,*)
    read (7,*) NX, NY, X0, XE, F, DZ_method ! DZ_method only accounts for the location of the nodes in the metal and therefore not nessesery for the stresses. 
    sig_d=0.2 ! Old parameter which is no longer used
    
    L1NX = NX/8                         ! Number of nodes at coursest level of the simulation. 
    L1NY = NY/8
    
    L1NY = L1NY/2                       ! Ensuring odd number of nodes in Y-direction to ensure a single middel line. 
    L1NY = L1NY*2+1
    
    NX = (L1NX-1)*8+1             
    NY = (L1NY-1)*8+1
    
    DX = (XE-X0)/(NX-1.)                ! Space increment
    Y0 = -0.5*DX*NY+0.5*DX              ! Width in Y-direction
    width = DX*NY                       ! Width of model
    

    ! EHL parameters -------------------------------------------------------------
    !RX=0.0136/2 ! 0.0136/2             ! Radius  acc Dave Hannes                           [m]
    !W0=10.0E5  ![N]                    ! Load parameter  to get ish Pherts=2GPa            [N/m]
    !Ua=23.0                            ! =U2 The speed of the surface with the asperity    [m/s]
    !Ub=17.0                            ! Velocity of lower surface                         [m/s]       ! Mean velocity is US=(Ua+Ub)/2
    read (7,*)
    read (7,*)
    IF(geom .LE. 5)THEN
        read (7,*) RX, W0, Ua, Ub, lub_param, shear_thin, lub_temp, solid_material  
        !  lub_temp  controles int temp simulations. Here the temp is read as input so not nessesery 
        ! solid_material controles the setup of the metal. 
        read (7,*)
        read (7,*)
        IF( geom .eq. 3) then
            read (7,*) b, by, Ph, RY    ! For geom 3 which is elliptical
            ph          = ph*1E9
            
        ELSE IF( Geom .EQ. 4)THEN
            read (7,*) B_ref,  by, PH
            Ph          =Ph*1e9
            B_ref       =B_ref*1e-4
            
        ELSE IF( Geom .EQ. 5)THEN
            read (7,*) PH_new, by, PH
            PH_new      =PH_new*1e9
            Ph          =Ph*1e9
            
        else
            read(7,*)                   ! No extra geometry input
        endif

    ELSE
        WRITE(4,*)'Bad Geometric entry. Geom ='
	    WRITE(4,*) geom
        stop 'Bad Geometric entry'
    endif
    
    Um=(Ua+Ub)/2                        ! Mean velocity of fluid                            [m/s]

    ! Asperity Data -------------------------------------------------------------
    !asph_real= 1.5E-6                  ! Real asperity height                              [m]
    !aspw_real= 50E-6                   ! Real asperity radius                              [m]
    read (7,*)
    read (7,*)
    read (7,*)  asph_real, aspw_real,asph_real2, aspw_real2, asp_ratio, Hminimum
    ! Converting
    asph_real  = asph_real*1E-6
    aspw_real  = aspw_real*1E-6
    asph_real2 = asph_real2*1E-6
    aspw_real2 = aspw_real2*1E-6
    Hminimum   = Hminimum *1E-6
    
    
    ! Elastic modulus -------------------------------------------------------------
    ! Elast1=206E9                      ! Elastic modulus of first body                     [Pa]
    ! Elast2=206E9                      ! Elastic modulus of second body                    [Pa]
    ! EE=2.26E11                        ! Equvivalent E =E'                                 [N/m^2]
    read (7,*)
    read (7,*)
    read (7,*) Elast1, Elast2, EE

    Elast1=Elast1*1E9
    Elast2=Elast2*1E9
    
    IF(EE .EQ. 0) THEN
        EE= 2/((1-0.3**2)/Elast1+(1-0.3**2)/Elast2)             ! If asperity case, ref2, Equvivalent Elastic modulus [Pa] Se Contact mech, KL Johnsson page 92. The book writes without a 2
    ElSE
        EE = EE*1E9                                             ! Equvivalent Elastic modulus [Pa]
    ENDIF
    
    ! Parameter calculations----------------------------------------------------------
    IF(Geom .EQ. 1)THEN         ! Cylinder
        W0    = W0*1E5
        nr_N  = 2*NX*NY/(XE-X0)                                 ! Number of nodes in the contact
        PAIAK = 0.15915489024
        G0    = width*PAI/2                                     ! G0 is a scale factors for the load, G0=2/3*pi = ball G0=pi/2=1D cylinder - see Tan X
        B     = SQRT(8*W0*RX/(PAI*EE))                          ! Contact halfwidth                     [m]
	    PH    = 2*W0/(PAI*B)                                    ! Hertzian pressure                     [N/m^2]
        
    ELSE IF(Geom .EQ. 2) THEN   ! Ball
        PAIAK = 0.2026423
        nr_N  = 2*NX*NY/(XE-X0)*1.0/(-Y0)                       ! Number of nodes in the contact
        
        G0    = 2.0/3.0*PAI                                     ! G0 is a scale factors for the load, G0=2/3*pi = ball G0=pi/2=1D cylinder - see Tan or Contact mech by KL Johnsson page 92.
        PH    = (3*W0*EE**2/(2*PAI**3*RX**2))**0.3333333        ! Hertzian pressure                     [N/m^2] se Contact mech by KL Johnsson page 93.
        B     = PAI*PH*RX/EE                                    ! Contact halfwidth                     [m] se Contact mech by KL Johnsson page 92.
        
    ELSE IF(Geom .EQ. 3) THEN   ! Elliptical Ball.              ! So far just based on EQs for a ball but with different radious
        PAIAK = 0.2026423
        nr_N  = 2*NX*NY/(XE-X0)*1.0/(-Y0)                       ! Number of nodes in the contact
        
        !G0=2.0/3.0*PAI                                         ! For elliptical contacts, this is defined later in Initi
        !Read as input                                          ! Hertzian pressure                     [N/m^2]
        !Read as input                                          ! Contact halfwidth                     [m]
        W0    = 2.0/3.0*Ph*pai*b*by                             ! Applied load, see equation (4.27) page 92 in KL Johnsson Contact mech
    Else IF(Geom .EQ. 4) THEN         ! Cylinder with given Ph
          
          ! PH    = 2*W0/(PAI*B)                                  ! Hertzian pressure   !Given in input.  
          PAIAK = 0.15915489024
          W0    = 2* PH**2 * PAI * RX / EE                        ! Applied load acc Johansson page 101 given that EE = 2*E'
          B     = SQRT(8*W0*RX/(PAI*EE))                          ! Contact halfwidth                     [m]
           
          XE    = XE * B_ref/B                                    ! Rescaling so that the real dx has the same length, not the dimensionless DX
          X0    = X0 * B_ref/B
          DX    = (XE-X0)/(NX-1.)                                 ! Space increment
          Y0    = -0.5*DX*NY+0.5*DX                               ! Width in Y-direction
          width = DX*NY                                           ! Width of model
          G0    = width*PAI/2                                     ! G0 is a scale factors for the load, G0=2/3*pi = ball G0=pi/2=1D cylinder - see Tan X
          nr_N  = 2*NX*NY/(XE-X0)                                 ! Number of nodes in the contact
          Geom  = 1                                               ! To continiue the calculations with a cylinder    
    Else IF(Geom .EQ. 5) THEN           ! Cylinder with given Ph - only load change
            
            W0    = 2* PH**2 * PAI * RX / EE
            nr_N  = 2*NX*NY/(XE-X0)                                 ! Number of nodes in the contact
            PAIAK = 0.15915489024
            G0    = width*PAI/2*(PH_new/PH)**2                                     ! G0 is a scale factors for the load, G0=2/3*pi = ball G0=pi/2=1D cylinder - see Tan X
            B     = SQRT(8*W0*RX/(PAI*EE))                          ! Contact halfwidth of Ph=2.3893356                    [m]
	        !PH    = 2*W0/(PAI*B)                                    ! Hertzian pressure                     [N/m^2]

    ELSE
        WRITE(4,*)'Bad Geometric entry. Geom ='
	    WRITE(4,*) geom
        stop 'Bad Geometric entry'
    ENDIF
    
    ! Accurace paramters ------------------------------------------------------------
    ! MK_stat = 200      ! Maximum numbers of itterations for static solution
    ! MK_time = 200      ! Maximum numbers of itterations for timedep solution
    ! ER_stat = 1E-4     ! Maximum numbers of itterations for static solution
    ! ER_time = 2E-5     ! Maximum numbers of itterations for timedep solution
    ! KK=15              ! Number of internal itterations
    
    read (7,*)
    read (7,*)
    read (7,*)     MK_stat, MK_time, ER_stat, ER_time, KK_stat, KK_time
    ER_stat = ER_stat*1E-8*nr_N               ! The error limit is set proportional to the number of nodes in the contact area
    ER_time = ER_time*1E-8*nr_N
    
    read (7,*)
    read (7,*)
    read (7,*)     C_meth, C_loc, C_glob, C_min, C_max, umax_p
    
    ! Load ballance ------------------------------------------------------------------
    ! H00=-0.1                  ! Initial offset for lubrication height 
    ! HM0r                      ! Reduction factor for uppd of H00
    read (7,*)
    read (7,*)
    read (7,*)  H00, HM0f!, H00_method, temp_r !Since these two parameters are not needed
    
    If ( geom .EQ. 6) Then ! Rescale h00 fore easier input
        H00=H00*Rx/b**2
    endif
        
    ! Choice of lubrication equation ------------------------------------------------
    Call Lubrication_def(HM0r,  Z, EDA0, alpha)
    
    ! End Input data -------------------------------------------------------------------------------------------------------------
    
    ! Pre calculations to generate dimentionless paramters ------------------------------------------------------------------------
    
    ! Geometry setup
    DT          = (Um*(XE-X0))/(Ntime*Ua)                           ! The time increment. T=t*Us/a. t_end=(XE-X0)*a/Ua => T_end=(XE-X0)Us/Ua
    SRR         = (Ua-Ub)/Um                                        ! Slide to Roll Ratio
    U2_over_Us  = Ua/(2.0*Um)                                       ! Alternative SRR number used by Venner
             
	U           = EDA0*Um/(EE*RX)                                   ! Dimensionless speed parameter. Venner writed U=EDA0*US/(2E*Rx) and US=2*Um.  
    
	W           = 2.*PAI*PH/(3.*EE)*(B/RX)**2
	
	AHM         = 1.0-EXP(-0.68*1.03)
	HM0         = 1.0!                                              ! Parameter for uppdating H0, the fim thicknes ofset value
    ! Alternative, 
    ! ALFA=Z*5.1E-9*A1
	! G=ALFA*EE
    ! HM0=3.63*(RX/B)**2*G**0.49*U**0.68*W**(-0.073)*AHM      
    HM0r        = HM0*HM0f                                          ! Reduced lubrication uppdation. 
    ENDA        = 12 *EDA0 *Um * RX**2 / ( B**3*PH)                 ! Epsilon in Reynolds equation. = 12.*(2*U)*(EE/PH)*(RX/B)**3 
	UTL         = EDA0*Um*RX/(B*B*2.E7)
    H00past     = H00
    
    ! Moes and Bosma non-dimentional paramters
    ! L = alpha * EE * (2*EDA0*Um/(EE*RX))**0.25
    ! M = (W0/(EE*RX**2))*(EE*RX/(2*EDA0*Um))**0.75
    
    ! Higginson and Dowsson non-dim param from Spikes review from 2006 for line load
    Uhd         = Um*EDA0/(EE*RX)                                   ! Dimensionless speed parameter. Venner writed U=EDA0*US/(2E*Rx) and US=2*Um.
    Ghd         = alpha*EE                                          ! Dimensionless material parameter
    alpha_bar   = alpha*PH
    L           = Ghd*(2.0*Uhd)**0.25                               ! same as L above
    
   IF(Geom .EQ. 1 .OR. Geom .EQ. 5)THEN
        !Acc C. H. Venner 1994 Transient Analysis of Surface line contact          
       Whd      = W0/(EE*RX) 
       M        = Whd/(2.0*Uhd)**(0.5) 
   ELSE
        !Acc C. H. Venner 1994 Nummerical simulations of a transverse ridge in a circular .......        
       Whd      = W0/(EE*RX**2) 
       M        = Whd/(2.0*Uhd)**(0.75)
   ENDIF
   
	
    ! Asperity  rescaling data
    asph        = asph_real*RX/B**2                                 ! Normalized height of asperity            [-]   
    aspw        = aspw_real/B                                       ! Normalized radius or wavelength of asperity            [-]
    asph2       = asph_real2*RX/B**2   
    aspw2       = aspw_real2/B        
    
    
    ! Only normal main file down to here, exept first 2 line
    
    ! Setup output files
    CALL file_definitions
    
    ! Start simulation setup
    CALL Setup
    
!---------When done -----------------------------------------------------------------------
    OPEN(20,FILE='1_Finished_stress_Wihoooo.DAT',STATUS='UNKNOWN')
        
    ! To Output file
    WRITE(20,*)'This is awesome'
    
        STOP
    END 