Subroutine Setup
        implicit none
        ! Input
        ! Calculations
        integer     NX1, NY1, temp_time, NN
        integer     Sstart
        real        tcvot, T40
        real*8      APlb(258,258), AHlb(258,258), ATlb(258,258), Yb(257)
        real*8      APl(258,98), AHl(258,98), ATl(258,98), X(257), Y(97)  !A needs to contain the full P matrix
        real*8      APly(258,50), AHly(258,50), ATly(258,50)
        real*8      AP(202,98), AH(202,98), AT(202,98)   !A needs to contain the full P matrix
        real*8      APm(154,50), AHm(154,50), ATm(154,50) 
        real*8      APs(106,50), AHs(106,50), ATs(106,50) !As(106,50), AAs(106,50)
        integer     ix, ii
        integer     Ntime_tot, nrows, extra_time
        ! Output 
        save        /CurrentP/
        save        /CurrentH/
        save        /CurrentT/
        save        /StressParam/ 
        include     'inc_Grid.h'
        include     'inc_itt.h'
        include     'inc_ref.h'    
        include     'inc_Yasutomi.h'
        include     'inc_CurrentP.h'
        include     'inc_CurrentH.h'
        include     'inc_CurrentT.h'
        include     'inc_StressParam.h'
   
        
        
        
        if(NX .NE. 105 .AND. NX .NE. 201 .AND. NX .NE. 257)then
        WRITE(*,*) ' Warning for not standard NX'
        WRITE(4,*) ' Warning for not standard NX'
     endif
     
     IF (NX .GT. 105) then
        Nz=21
        NX1=NX+1
        NY1=Ny+1
        Nr_xy=4                                 ! Number of nodes in each spatial direction for which we increase the stress evaluation
        Nr_z=Nr_xy                                 ! Number of DX depth that we refine the mesh. Has not to be equal to Nr_xy since rescaling the z direction if close to the surface
        Nz_cvot=5                              ! Has to be an odd number
     else
        Nz=26
        NX1=NX+1
        NY1=Ny+1
        Nr_xy=2                                 ! Number of nodes in each spatial direction for which we increase the stress evaluation
        Nr_z=Nr_xy                                 ! Number of DX depth that we refine the mesh
        Nz_cvot=3!13                             ! Has to be an odd number
        
     endif
     NN=(NY+1)/2
     dh=DX
     sig_d=dh*(Nz-(Nz_cvot-1)*NR_z)         !For each Nr_z it takes Nz_cvot steps to go 1 step in DX
       
        !Compatible with version 5917 and later for lub_param eq. 5
        IF(NX .GT. 200 .AND. NY .GT. 160)THEN
            Sstart=1                                ! Start level of simulation
        ELSE
            Sstart=2                                ! The coursest solution should not be too course. 
        ENDIF
        
        T40=40
        tcvot=0
        If(lub_param .EQ. 5 .and. Ta .GT. T40) then
                temp_time=0
            DO WHILE ((Ta-T40)*tcvot .LT. Ta-T40)
                IF( (Ta-T40)*tcvot .LT. 20) THEN
                    tcvot   = tcvot+4/(Ta-T40)      ! Increasing 4 degreas at a time
                    temp_time=temp_time+1
                elseif ( (Ta-T40)*tcvot .LT. 40) THEN
                    tcvot   = tcvot+2/(Ta-T40)      ! Increasing 2 degreas at a time
                    temp_time=temp_time+1
                else
                    tcvot   = tcvot+1/(Ta-T40)      ! Increasing 1 degrea at a time
                    temp_time=temp_time+1
                ENDIF
            enddo
            
            
            extra_time=5-Sstart+temp_time+10+1      ! 5-sstart itterations for refining the results, 20 itterations for increasing the temperature, 10 itterations for adding time dependencd +1 for that time simulations start at t=0 and ii starts at 1
            WRITE(4,*)'Non Newtioian simulation with temperature increase. extra time = ',extra_time
        else
            extra_time=5-Sstart+10+1
            WRITE(4,*)'No temperature increase. extra time = ', extra_time
        endif
        
        !If asp_shape = 0, calculate the stresses for only the last time step
        if( asp_shape .EQ. 0) then
            Ntime=extra_time-(5-Sstart)
            extra_time=5-Sstart-1
        endif

        ntime_tot=Ntime+extra_time
        if( asp_shape .ge. 200 .and. asp_shape .le. 209 ) ntime_tot = ntime_tot + (135-15)*Ntime/Nx  ! To account for the larger model if rough surfaces
        
        CALL Short_printout(extra_time, Ntime_tot)
        
!-------Input data --------------------------------------------------------------------------------
        open (unit = 2, file = "PRESSURE.DAT")              ! To read
        open (unit = 3, file = "FILM.DAT")                  ! To read
        open (unit = 5, file = "Temp.DAT")                  ! To read
    
!-------Start the stresscalculations---------------------------------------------------------------
        nrows=NX1
    
        if( NY .GT. 201 .and. NY .eq. NX)THEN
            WRITE(4,*)'Ball calculation with assumed NX=NY=257'
            WRITE(4,*)'We have Nx = ', NX, ' and NY = ', NY
        do ii = 1,ntime_tot               
            
            do ix = 1,nrows
                read(2,*) APlb(ix,:)
                read(3,*) AHlb(ix,:)
                read(5,*) ATlb(ix,:)
            end do
            
            if (ii .eq. extra_time) then
                Yb=APlb(1,2:NY1)
                X=APlb(2:NX1,1)
            endif
            
            if(ii .GE. extra_time)then
              !if(ii .GE. 131+extra_time)then  
                P(1:NX,1:NN)=APlb(2:NX1,2:NY1/2)
                H(1:NX,1:NN)=AHlb(2:NX1,2:NY1/2)
                Temp(1:NX,1:NN)=ATlb(2:NX1,2:NY1/2)
                CALL sigmacalc
            
                Call Maxcheck(ii-extra_time,ntime)
            endif
            
            WRITE(4,*)'Now finished with time = ', ii-extra_time
            WRITE(*,*)'Now finished with time = ', ii-extra_time
            
        enddo
        else if( NX .GT. 201)THEN
            if( NY .GT. 50) then
        do ii = 1,ntime_tot               
            
            do ix = 1,nrows
                read(2,*) APl(ix,:)
                read(3,*) AHl(ix,:)
                read(5,*) ATl(ix,:)
            end do
            
            if (ii .eq. extra_time) then
                Y=APl(1,2:NY1)
                X=APl(2:NX1,1)
            endif
            
            if(ii .GE. extra_time)then
              !if(ii .GE. 131+extra_time)then  
                P(1:NX,1:NN)=APl(2:NX1,2:NY1/2)
                H(1:NX,1:NN)=AHl(2:NX1,2:NY1/2)
                Temp(1:NX,1:NN)=ATl(2:NX1,2:NY1/2)
                CALL sigmacalc
            
                Call Maxcheck(ii-extra_time,ntime)
            endif
            
            WRITE(4,*)'Now finished with time = ', ii-extra_time
            WRITE(*,*)'Now finished with time = ', ii-extra_time
            
            enddo
        else 
            do ii = 1,ntime_tot               
            
            do ix = 1,nrows
                read(2,*) APly(ix,:)
                read(3,*) AHly(ix,:)
                read(5,*) ATly(ix,:)
            end do
            
            if (ii .eq. extra_time) then
               ! Y=APly(1,2:NY1)
                X=APly(2:NX1,1)
            endif
            
            if(ii .GE. extra_time)then
              !if(ii .GE. 131+extra_time)then  
                P(1:NX,1:NN)=APly(2:NX1,2:NY1/2)
                H(1:NX,1:NN)=AHly(2:NX1,2:NY1/2)
                Temp(1:NX,1:NN)=ATly(2:NX1,2:NY1/2)
                CALL sigmacalc
            
                Call Maxcheck(ii-extra_time,ntime)
            endif
            
            WRITE(4,*)'Now finished with time = ', ii-extra_time
            WRITE(*,*)'Now finished with time = ', ii-extra_time
            
        enddo
        endif
        
        else if( NX .GT. 153)THEN
        do ii = 1,ntime_tot               
            
            do ix = 1,nrows
                read(2,*) AP(ix,:)
                read(3,*) AH(ix,:)
                read(5,*) AT(ix,:)
            end do
            
            if (ii .eq. extra_time) then
                Y=AP(1,2:NY1)
                X=AP(2:NX1,1)
            endif
            
            if(ii .GE. extra_time)then
              !if(ii .GE. 131+extra_time)then  
                P(1:NX,1:NN)=AP(2:NX1,2:NY1/2)
                H(1:NX,1:NN)=AH(2:NX1,2:NY1/2)
                Temp(1:NX,1:NN)=AT(2:NX1,2:NY1/2)
                CALL sigmacalc
            
                Call Maxcheck(ii-extra_time,ntime)
            endif
            
            WRITE(4,*)'Now finished with time = ', ii-extra_time
            WRITE(*,*)'Now finished with time = ', ii-extra_time
            
        enddo
        
        else if( NX .GT. 105)THEN
        do ii = 1,ntime_tot               
            
            do ix = 1,nrows
                read(2,*) APm(ix,:)
                read(3,*) AHm(ix,:)
                read(5,*) ATm(ix,:)
            end do
            
            if (ii .eq. extra_time) then
                Y=APm(1,2:NY1)
                X=APm(2:NX1,1)
            endif
            
            if(ii .GE. extra_time)then
              !if(ii .GE. 131+extra_time)then  
                P(1:NX,1:NN)=APm(2:NX1,2:NY1/2)
                H(1:NX,1:NN)=AHm(2:NX1,2:NY1/2)
                Temp(1:NX,1:NN)=ATm(2:NX1,2:NY1/2)
                CALL sigmacalc
            
                Call Maxcheck(ii-extra_time,ntime)
            endif
            
            WRITE(4,*)'Now finished with time = ', ii-extra_time
            WRITE(*,*)'Now finished with time = ', ii-extra_time
            
        enddo
        
        
        else        ! For smaller mesh size
            do ii = 1,ntime_tot               
            
            do ix = 1,nrows
                read(2,*) APs(ix,:)
                read(3,*) AHs(ix,:)
                read(5,*) ATs(ix,:)
            end do
            
            if (ii .eq. extra_time) then
                Y(1:NY)=APs(1,2:NY1)
                X(1:NX)=APs(2:NX1,1)
            endif
            
            if(ii .GE. extra_time)then
                P(1:NX,1:NN)=APs(2:NX1,2:NY1/2)
                H(1:NX,1:NN)=AHs(2:NX1,2:NY1/2)
                Temp(1:NX,1:NN)=ATs(2:NX1,2:NY1/2)
                
                ! For debugging
                !H=0*H+1
                !P=0*P
                !P(50,25)=1
                !P(40:70,:)=1
                !P(10,25)=0.01
                !P(90,25)=0.01
                !do i=1,NX
                !    P(i,:)=max(1.0-X(i)*X(i),0.0)
                !enddo
                !P(10,1)=1
                
                CALL sigmacalc
            
                Call Maxcheck(ii-extra_time,ntime)
            endif
            
            WRITE(4,*)'Now finished with time = ', ii-extra_time
            WRITE(*,*)'Now finished with time = ', ii-extra_time
            
            enddo
        endif
    
            
        Call Printmax
            
        return
    end
    
        