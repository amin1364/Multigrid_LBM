	!This is the source code for the article "Application of Multiple Time-Step 
	!Grid Lattice Boltzmann Method to Simulate Natural Convection in a Wide." 
	!This program is written by Seyed Amin Nabavizadeh


	!
	!	 ========================================================
    PROGRAM LBLENS
    USE D2Q9Const
    USE simParam
    USE m_veloc
    USE m_temp
    IMPLICIT NONE
	DOUBLE PRECISION::nu_ave_left,nu_ave_right
    INTEGER:: tStep,inner_loop,y
    INTEGER*4 now(3)
	OPEN(2,FILE='time_of_Simulation.txt')
	OPEN(3,FILE='average_nu_right.txt')
	OPEN(4,FILE='average_nu_left.txt')
	OPEN(5,FILE='Nu_left_last_iteration.txt')
	OPEN(6,FILE='Nu_right_last_iteration.txt')
	write(4,*)'FO,ave_Nu_left'
	write(3,*)'FO,ave_Nu_right'
	!This subroutine calculates lattice constants cs (speed of sound) for D2Q9 lattice.
    CALL calculatecons() 
	
	!This subroutine located in simParam.f90 load the variable and calculate parameters 
	!such as relaxation time, etc.
    CALL ComputeParam ()
	
	!This subroutine initializes the fluid flow variable and allocates memory to them.
    CALL veloc_init ()
	
	!This subroutine defines the wall boundary condition.
    CALL soliddeff(solid)
	
	!This subroutine initializes the heat transfer variable and allocates memory to them.
    CALL temp_init
	
	! "tpghost" is temperature field in fluid flow solver. This variable is required for calculating force term.
    IF (xDim.eq.xDim_thermal) THEN
		tpghost=tp
    ELSEIF (xDim_thermal<xDim) THEN
		CALL interpolate2d_usedin_coarse_to_fine(tp,tpghost)
	ELSEIF (xDim_thermal>xDim) THEN
		CALL interpolate2d_usedin_fine_to_coarse(tp,tpghost)
    END IF


	! WRITE the initial condition in the output file.
    CALL writeoutput(tpghost,uTot,0)
	CALL itime(now)
	WRITE(2,1000)now
	

	! The main part of the code
    DO tStep = 1, tMax
		DO inner_loop=1,n_ratio_u
			call veloc_main
		END DO
		DO inner_loop=1,n_ratio
			CALL temp_main()
		END DO
		!
		IF (xDim.eq.xDim_thermal)THEN
			tpghost=tp
		ELSEIF(xDim<xDim_thermal)THEN
			CALL interpolate2d_usedin_fine_to_coarse(tp,tpghost)
		ELSEIF(xDim>xDim_thermal)THEN
			CALL interpolate2d_usedin_coarse_to_fine(tp,tpghost)
		END IF



		IF(MOD(tStep,kPlot).eq.0)THEN
			CALL writeoutput(tpghost,uTot,tStep)
			CALL calculate_nu(nu_left,nu_right,nu_ave_left,nu_ave_right,tp)
			CALL outputcoarse(tp,tStep)
			CALL itime(now)
			WRITE(2,3000)tstep,now
			if(n_ratio>1) then
				write(3,"(f9.5,',',f9.6)")(tau_t-.5d0)/3.0d0*dble(tstep*n_ratio)/dble(xDim_thermal)/dble(xDim_thermal),nu_ave_right
				write(4,"(f9.5,',',f9.6)")(tau_t-.5d0)/3.0d0*dble(tstep*n_ratio)/dble(xDim_thermal)/dble(xDim_thermal),nu_ave_left
			else
				write(3,"(f9.5,',',f9.6)")(tau_t-.5d0)/3.0d0*dble(tstep)/dble(xDim_thermal)/dble(xDim_thermal),nu_ave_right
				write(4,"(f9.5,',',f9.6)")(tau_t-.5d0)/3.0d0*dble(tstep)/dble(xDim_thermal)/dble(xDim_thermal),nu_ave_left
			endif
		ENDIF
		
		IF (tstep.eq.tMax) THEN
			write(5,*)'Y,Nu_left'
			write(6,*)'Y,Nu_right'
			DO y=1,yDim_thermal
				write(6,"(i4,',',f9.6)")y,nu_right(y)
				write(5,"(i4,',',f9.6)")y,nu_left(y)
			END DO
		ENDIF
		
		WRITE(*,2000)tStep
	3000 FORMAT('time at Iteration Step =',i9,':     ',(i2.2, ':', i2.2, ':', i2.2))
	1000 FORMAT('time at Iteration Step =0000000:     ',(i2.2, ':', i2.2, ':', i2.2))
	2000 FORMAT('Iteration Step =',I8)
    END DO

	close(2)
	close(3)
	close(4)
	close(5)
	close(6)

    END PROGRAM LBLENS
	
	SUBROUTINE calculate_nu(nu_left,nu_right,nu_ave_left,nu_ave_right,tp)
	USE simParam, ONLY: xDim_thermal, yDim_thermal
	IMPLICIT NONE
	
	DOUBLE PRECISION, INTENT(OUT):: nu_left(yDim_thermal), nu_right(yDim_thermal)
    DOUBLE PRECISION, INTENT(OUT):: nu_ave_left,nu_ave_right
	DOUBLE PRECISION, INTENT(IN):: tp(yDim_thermal,xDim_thermal)
	integer:: y,x
	DO y=1,yDim_thermal
		nu_left(y)=abs(tp(y,2)-tp(y,1))*dble(xDim_thermal);
		nu_right(y)=abs(tp(y,xDim_thermal-1)-tp(y,xDim_thermal))*dble(xDim_thermal);
	END DO
	nu_ave_left=sum(nu_left)/size(nu_left);
	nu_ave_right=sum(nu_right)/size(nu_right);
	
	
	END SUBROUTINE calculate_nu



