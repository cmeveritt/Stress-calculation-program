subroutine Short_printout(extra_time, ntime_tot)
    implicit none
        integer  extra_time, ntime_tot
        include     'inc_Grid.h'
        include     'inc_StressParam.h'
        include     'inc_Itt.h'
        
    !------- Write in outputs -------------------------------------------------------------------        
        WRITE(4,*)'Grid infromation'
        WRITE(4,*)'Nz, sig_d, NX, NY, Ntime'
	    WRITE(4,*) Nz, sig_d, NX, NY, Ntime
        WRITE(4,*)'Nr_z, NZ_cvot'
	    WRITE(4,*) Nr_z, NZ_cvot
        WRITE(4,*)'extra_time, total time'
	    WRITE(4,*) extra_time, ntime_tot

        WRITE(*,*)'Grid infromation'
        WRITE(*,*)'Nz, sig_d, NX, NY, Ntime'
	    WRITE(*,*) Nz, sig_d, NX, NY, Ntime
        WRITE(*,*)'Nr_z, NZ_cvot'
	    WRITE(*,*) Nr_z, NZ_cvot
        WRITE(*,*) ' Now we calculate'
        WRITE(*,*)'extra_time'
	    WRITE(*,*) extra_time
        return
    end
    