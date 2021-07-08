	! All subroutines related to energy equation and interpolation scheme

    ! initialize temperature variables and parametrs for calculating interpolation coefficents in each node
    SUBROUTINE temp_init
    !---------------------------------------------------------------------
    
    USE m_temp
    USE m_veloc, ONLY: uTot
    USE simParam, ONLY: xDim_thermal, yDim_thermal, U_ratio,xDim, yDim
    REAL*8::nn
    INTEGER::x,y
    DOUBLE PRECISION::utemp(yDim_thermal,xDim_thermal,2)
    ALLOCATE(ft(yDim_thermal,xDim_thermal,9))
    ALLOCATE(fEqt(yDim_thermal,xDim_thermal,9))
    ALLOCATE(tp(yDim_thermal,xDim_thermal))
    ALLOCATE(tpghost(yDim,xDim))
	ALLOCATE(uSqr_t(yDim_thermal,xDim_thermal))
    ALLOCATE(location_thermal_x(xDim_thermal))
    ALLOCATE(location_thermal_y(yDim_thermal))
    ALLOCATE(location_x(xDim))
    ALLOCATE(location_y(yDim))
	ALLOCATE(nu_left(yDim_thermal))
    ALLOCATE(nu_right(yDim_thermal))
    ALLOCATE(weight_t_fine(yDim*xDim,6))
    ALLOCATE(weight_t_coarse(yDim_thermal*xDim_thermal,6))
    weight_t_fine=0;
	
	
    tp=0.0d0;
	
	! find the constants for interpolation at each node
    nn=1.0/real(yDim-1)
    location_y(1)=0
    DO y=2,yDim
        location_y(y)=location_y(y-1)+nn
    END DO

    nn=1.0/real(xDim-1)
    location_x(1)=0
    DO x=2,xDim
        location_x(x)=location_x(x-1)+nn
    END DO



    nn=1.0/real(yDim_thermal-1)
    location_thermal_y(1)=0
    DO y=2,yDim_thermal
        location_thermal_y(y)=location_thermal_y(y-1)+nn
    END DO

    nn=1.0/real(xDim_thermal-1)
    location_thermal_x(1)=0
    DO x=2,xDim_thermal
        location_thermal_x(x)=location_thermal_x(x-1)+nn
    END DO

    IF (xDim.ne.xDim_thermal)THEN
        CALL findweight_usedin_coarse_to_fine
        CALL findweight_usedin_fine_to_coarse
    ENDIF



	! Temperature Initialization
    DO x=1,xDim_thermal
        tp(:,x)=DBLE(x-1)/DBLE(xDim_thermal-1)
    END DO
	
	

    IF (xDim.eq.xDim_thermal)THEN
        utemp=uTot
    ELSEIF (xDim>xDim_thermal) THEN
		! calculate the velocity field for coarser grid "utemp"
        CALL interpolate2d_usedin_fine_to_coarse_vector(uTot,utemp)
		
	ELSEIF (xDim<xDim_thermal) THEN
		CALL interpolate2d_usedin_coarse_to_fine_vector(uTot,utemp)
    ENDIF
	uSqr_t(:,:)=utemp(:,:,1)**2+utemp(:,:,2)**2;
    CALL computeFeqT(fEqt,tp,utemp,uSqr_t*U_ratio**2)

    ft = fEqt
    END SUBROUTINE temp_init

    !	 ========================================================
    !	 The main scheme of implementing temprature in LB
    !	 ========================================================

    SUBROUTINE temp_main

    USE m_temp
    USE m_veloc, ONLY: uTot
    USE simParam, ONLY:  tau_t,U_ratio,yDim_thermal,xDim_thermal,xDim

    IMPLICIT NONE
    DOUBLE PRECISION::utemp(yDim_thermal,xDim_thermal,2)
    IF (xDim.eq.xDim_thermal)THEN
        utemp=uTot
    ELSEIF (xDim>xDim_thermal) THEN
		! calculate the velocity field for coarser grid "utemp"
        CALL interpolate2d_usedin_fine_to_coarse_vector(uTot,utemp)
		
	ELSEIF (xDim<xDim_thermal) THEN
		CALL interpolate2d_usedin_coarse_to_fine_vector(uTot,utemp)
    ENDIF
	uSqr_t(:,:)=utemp(:,:,1)**2+utemp(:,:,2)**2;
	CALL computeFeqT(fEqt,tp,utemp*U_ratio,uSqr_t*U_ratio**2)
	CALL collideT(ft,fEqt,tau_t)
	CALL streamT(ft)
	CALL boundariesT(ft,tp)
	
	CALL computeMacroT(ft,tp)
	
	
    
    
    
    
    

    END SUBROUTINE temp_main





    !	 ========================================================
    !	 Compute equilibrium distribution
    !	 ========================================================
    SUBROUTINE computeFeqT(fEqt,tp,uTot,uSqr_t)
    USE D2Q9COnst, ONLY: t, v
    USE simParam, ONLY: xDim_thermal, yDim_thermal, tp_u
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: uSqr_t(yDim_thermal,xDim_thermal), uTot(yDim_thermal,xDim_thermal,0:1), tp(yDim_thermal,xDim_thermal)
    DOUBLE PRECISION, INTENT(INOUT):: fEqt(yDim_thermal,xDim_thermal,0:8)
    INTEGER:: i, x, y
    DOUBLE PRECISION:: uxy

    DO i = 0, 8
        DO x = 1, xDim_thermal
            DO y = 1, yDim_thermal
                uxy = uTot(y,x,0) * v(i,1) + uTot(y,x,1) * v(i,2)
                IF(tp_u.eq.1)then
                    fEqt(y,x,i) = t(i) * tp(y,x) * (1.0d0 + 3.0d0 * uxy + 4.5d0 * uxy * uxy - 1.5d0 * uSqr_t(y,x))
                ELSE
                    fEqt(y,x,i) = t(i) * tp(y,x) * (1.0d0)
                endif
            END DO
        END DO
    END DO
    END SUBROUTINE computeFeqT


    !	 ========================================================
    !	 Compute Temprature from distribution functions
    !	 ========================================================
    SUBROUTINE computeMacroT(ft,tp)
    USE simParam, ONLY: xDim_thermal, yDim_thermal
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: ft(yDim_thermal,xDim_thermal,0:8)
    DOUBLE PRECISION, INTENT(INOUT):: tp(yDim_thermal,xDim_thermal)
    INTEGER:: x,y

    DO x = 1, xDim_thermal
        DO y = 1, yDim_thermal
            tp(y,x)  = ft(y,x,0) + ft(y,x,1) + ft(y,x,2) + ft(y,x,3) + ft(y,x,4) + ft(y,x,5) + ft(y,x,6) + ft(y,x,7) + ft(y,x,8)
        END DO
    END DO
    END SUBROUTINE computeMacroT


    !	 ========================================================
    !	 Compute Boundary Condition
    !	 ========================================================
    SUBROUTINE boundariesT(ft,tp)
    USE simParam

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(INOUT):: ft(yDim_thermal,xDim_thermal,0:8), tp(yDim_thermal,xDim_thermal)


    INTEGER:: x, y

    DO x = 1, xDim_thermal
        DO y = 1, yDim_thermal
            IF (x == 1) then
                tp(y,x) = 0.0;
                CALL leftwall(ft(y,x,:),tp(y,x))
            ELSE IF (x == xDim_thermal) then
                tp(y,x) = 1.0;
                CALL right_wall(ft(y,x,:),tp(y,x))
            ELSEIF (y == 1) then
                tp(y,x) = tp(y+1,x);
                CALL bottomwall(ft(y,x,:),tp(y,x))
            ELSE IF (y == yDim_thermal) then
                tp(y,x) = tp(y-1,x);
                CALL topwall(ft(y,x,:),tp(y,x))
            END IF
        END DO
    END DO



    CONTAINS


    SUBROUTINE leftwall(ft,tp)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(INOUT):: ft(0:8)
    DOUBLE PRECISION, INTENT(IN):: tp(1)
    DOUBLE PRECISION:: t_bar

    t_bar = 6.0d0*(tp(1)-(ft(0)+ft(2)+ft(4)+ft(3)+ft(6)+ft(7)))
    ft(1) = (1.0d0/9.0d0)*t_bar
    ft(5) = (1.0d0/36.0d0)*t_bar
    ft(8) = (1.0d0/36.0d0)*t_bar
    END SUBROUTINE leftwall


    SUBROUTINE right_wall(ft,tp)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(INOUT):: ft(0:8)
    DOUBLE PRECISION, INTENT(IN):: tp(1)
    DOUBLE PRECISION:: t_bar

    t_bar = 6.0d0*(tp(1)-(ft(0)+ft(2)+ft(4)+ft(1)+ft(5)+ft(8)))
    ft(3) = (1.0d0/9.0d0)*t_bar
    ft(6) = (1.0d0/36.0d0)*t_bar
    ft(7) = (1.0d0/36.0d0)*t_bar
    END SUBROUTINE right_wall




    SUBROUTINE topwall(ft,tp)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(INOUT):: ft(0:8)
    DOUBLE PRECISION, INTENT(IN):: tp(1)
    DOUBLE PRECISION:: t_bar

    t_bar = 6.0d0*(tp(1)-(ft(0)+ft(1)+ft(2)+ft(3)+ft(5)+ft(6)))
    ft(4) = (1.0d0/9.0d0)*t_bar
    ft(7) = (1.0d0/36.0d0)*t_bar
    ft(8) = (1.0d0/36.0d0)*t_bar
    END SUBROUTINE topwall



    SUBROUTINE bottomwall(ft,tp)
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(INOUT):: ft(0:8)
    DOUBLE PRECISION, INTENT(IN):: tp(1)
    DOUBLE PRECISION:: t_bar

    t_bar = 6.0d0*(tp(1)-(ft(0)+ft(1)+ft(3)+ft(4)+ft(7)+ft(8)))
    ft(2) = (1.0d0/9.0d0)*t_bar
    ft(5) = (1.0d0/36.0d0)*t_bar
    ft(6) = (1.0d0/36.0d0)*t_bar
    END SUBROUTINE bottomwall
    END SUBROUTINE boundariesT


    !	 ========================================================
    !	 Streaming step: 
    !	 ========================================================
SUBROUTINE streamT(ft)
    USE simParam
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(INOUT):: ft(yDim_thermal,xDim_thermal,0:8)
	DOUBLE PRECISION::ft_t(yDim_thermal,xDim_thermal,0:8)
    INTEGER :: x,y,y_n,x_e,y_s,x_w

    ft_t=ft;
    DO x = 1, xDim_thermal
       DO y = 1, yDim_thermal
       
          y_n = MOD(y+1,yDim_thermal);  IF(y_n.eq.0)y_n = yDim_thermal
          x_e = MOD(x+1,xDim_thermal);  IF(x_e.eq.0)x_e = xDim_thermal
          y_s = MOD(y-1,yDim_thermal);  IF(y_s.eq.0)y_s = yDim_thermal
          x_w = MOD(x-1,xDim_thermal);  IF(x_w.eq.0)x_w = xDim_thermal
                 
		  ft(y,x,0) = ft_t(y,x,0)       !.........zero: just copy
		  ft(y,x_e,1) = ft_t(y,x,1)     !.........east
		  ft(y_n,x,2) = ft_t(y,x,2)     !.........north
		  ft(y,x_w,3) = ft_t(y,x,3)     !.........west
		  ft(y_s,x,4) = ft_t(y,x,4)     !.........south
		  ft(y_n,x_e,5) = ft_t(y,x,5)   !.........north-east
		  ft(y_n,x_w,6) = ft_t(y,x,6)   !.........north-west
		  ft(y_s,x_w,7) = ft_t(y,x,7)   !.........south-west
		  ft(y_s,x_e,8) = ft_t(y,x,8)   !.........south-east   
		   
       END DO
    END DO
    
END SUBROUTINE streamT


    !	 ========================================================
    !	 LBGK collision step
    !	 ========================================================
    SUBROUTINE collideT(ft,fEqt,tau_t)
    USE simParam, ONLY: xDim_thermal, yDim_thermal
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: fEqt(yDim_thermal,xDim_thermal,0:8)
    DOUBLE PRECISION, INTENT(INOUT):: ft(yDim_thermal,xDim_thermal,0:8)
    DOUBLE PRECISION :: tau_t
    INTEGER:: x,y,i

    DO i = 0, 8
        DO x = 1, xDim_thermal
            DO y = 1, yDim_thermal
                ft(y,x,i) = (1.0d0 - 1.0d0/tau_t) * ft(y,x,i) + 1.0d0/tau_t * feqt(y,x,i)
            END DO
        END DO
    END DO
    END SUBROUTINE collideT








    !	 ========================================================
    !	 Interpolation Subroutines. The method of implementation was discussed in the article.
    !	 ========================================================
    SUBROUTINE findweight_usedin_coarse_to_fine
    USE simParam, ONLY:xDim_thermal, yDim_thermal,xDim, yDim
    USE m_temp, ONLY: weight_t_fine,location_x,location_y,location_thermal_x,location_thermal_y
    IMPLICIT NONE
    INTEGER::y,x,ii,jj
    REAL*8::d1,d2,d3,d4
    REAL*8::npx,dummy,npy
    INTEGER::counter
    weight_t_fine=0.0
    counter=0
    npy=real(yDim-1)/real(yDim_thermal-1)
    npx=real(xDim-1)/real(xDim_thermal-1)
    DO ii=1,xDim !fine mesh
        DO jj=1,yDim !fine mesh
            IF((jj.lt.yDim).and.(ii.lt.xDim))THEN
                counter=counter+1
                y=FLOOR(real(jj-1)/real(npy))+1 !coarse mesh
                x=FLOOR(real(ii-1)/real(npx))+1 !coarse mesh
                d1=((location_y(jj)-location_thermal_y(y))**2+(location_x(ii)-location_thermal_x(x))**2)**0.5
                d2=((location_y(jj)-location_thermal_y(y+1))**2+(location_x(ii)-location_thermal_x(x))**2)**0.5
                d3=((location_y(jj)-location_thermal_y(y))**2+(location_x(ii)-location_thermal_x(x+1))**2)**0.5
                d4=((location_y(jj)-location_thermal_y(y+1))**2+(location_x(ii)-location_thermal_x(x+1))**2)**0.5
                weight_t_fine(counter,1)=x
                weight_t_fine(counter,2)=y
                dummy=1.0/(1+d1/d2+d1/d3+d1/d4)
                weight_t_fine(counter,3)=dummy
                weight_t_fine(counter,4)=dummy*d1/d2
                weight_t_fine(counter,5)=dummy*d1/d3
                weight_t_fine(counter,6)=dummy*d1/d4
            ELSEIF((jj.eq.yDim).and.(ii.lt.xDim))THEN
                counter=counter+1
                y=FLOOR(real(jj-1)/real(npy))+1 !coarse mesh
                x=FLOOR(real(ii-1)/real(npx))+1 !coarse mesh
                d1=((location_y(jj)-location_thermal_y(y))**2+(location_x(ii)-location_thermal_x(x))**2)**0.5
                d2=((location_y(jj)-location_thermal_y(y))**2+(location_x(ii)-location_thermal_x(x+1))**2)**0.5
                weight_t_fine(counter,1)=x
                weight_t_fine(counter,2)=y
                dummy=1.0/(1+d1/d2)
                weight_t_fine(counter,3)=dummy
                weight_t_fine(counter,4)=dummy*d1/d2
                weight_t_fine(counter,5)=0
                weight_t_fine(counter,6)=0
            ELSEIF((jj.lt.yDim).and.(ii.eq.xDim))THEN
                counter=counter+1
                y=FLOOR(real(jj-1)/real(npy))+1 !coarse mesh
                x=FLOOR(real(ii-1)/real(npx))+1 !coarse mesh
                d1=((location_y(jj)-location_thermal_y(y))**2+(location_x(ii)-location_thermal_x(x))**2)**0.5
                d2=((location_y(jj)-location_thermal_y(y+1))**2+(location_x(ii)-location_thermal_x(x))**2)**0.5
                weight_t_fine(counter,1)=x
                weight_t_fine(counter,2)=y
                dummy=1.0/(1+d1/d2)
                weight_t_fine(counter,3)=dummy
                weight_t_fine(counter,4)=dummy*d1/d2
                weight_t_fine(counter,5)=0
                weight_t_fine(counter,6)=0
            ELSEIF((jj.eq.yDim).and.(ii.eq.xDim))THEN
                counter=counter+1
                y=FLOOR(real(jj-1)/real(npy))+1 !coarse mesh
                x=FLOOR(real(ii-1)/real(npx))+1 !coarse mesh
                weight_t_fine(counter,1)=x
                weight_t_fine(counter,2)=y
                weight_t_fine(counter,3)=1
                weight_t_fine(counter,4)=0
                weight_t_fine(counter,5)=0
                weight_t_fine(counter,6)=0
            ENDIF
        END DO
    END DO

    END SUBROUTINE findweight_usedin_coarse_to_fine





    SUBROUTINE interpolate2d_usedin_coarse_to_fine_vector(utot_coarse,utot_fine)
    USE simParam, ONLY:xDim_thermal, yDim_thermal,xDim, yDim
    USE m_temp, ONLY: weight_t_coarse
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN):: utot_coarse(yDim,xDim,2)
    DOUBLE PRECISION, INTENT(OUT):: utot_fine(yDim_thermal,xDim_thermal,2)
    INTEGER::y,x,counter,ii,jj
    DOUBLE PRECISION::w1,w2,w3,w4
    counter=0
    DO ii=1,xDim_thermal
        DO jj=1,yDim_thermal
            IF((jj.lt.yDim_thermal).and.(ii.lt.xDim_thermal))THEN
                counter=counter+1
                
                y=nint(weight_t_coarse(counter,2))
                x=nint(weight_t_coarse(counter,1))
                w1=weight_t_coarse(counter,3)
                w2=weight_t_coarse(counter,4)
                w3=weight_t_coarse(counter,5)
                w4=weight_t_coarse(counter,6)
                utot_fine(jj,ii,:)=w1*utot_coarse(y,x,:)+w2*utot_coarse(y+1,x,:)+w3*utot_coarse(y,x+1,:)+w4*utot_coarse(y+1,x+1,:)
            ELSEIF((jj.eq.yDim_thermal).and.(ii.lt.xDim_thermal))THEN
                counter=counter+1
                y=nint(weight_t_coarse(counter,2))
                x=nint(weight_t_coarse(counter,1))
                w1=weight_t_coarse(counter,3)
                w2=weight_t_coarse(counter,4)
                utot_fine(jj,ii,:)=w1*utot_coarse(y,x,:)+w2*utot_coarse(y,x+1,:)
            ELSEIF((jj.lt.yDim_thermal).and.(ii.eq.xDim_thermal))THEN
                counter=counter+1
                y=nint(weight_t_coarse(counter,2))
                x=nint(weight_t_coarse(counter,1))
                w1=weight_t_coarse(counter,3)
                w2=weight_t_coarse(counter,4)
                utot_fine(jj,ii,:)=w1*utot_coarse(y,x,:)+w2*utot_coarse(y+1,x,:)
            ELSEIF((jj.eq.yDim_thermal).and.(ii.eq.xDim_thermal))THEN
                counter=counter+1
                y=nint(weight_t_coarse(counter,2))
                x=nint(weight_t_coarse(counter,1))
                utot_fine(jj,ii,:)=utot_coarse(y,x,:)
            ENDIF
        END DO
    END DO
    END SUBROUTINE interpolate2d_usedin_coarse_to_fine_vector






    SUBROUTINE interpolate2d_usedin_fine_to_coarse(tp,tpghost)
    USE simParam, ONLY:xDim_thermal, yDim_thermal,xDim, yDim
    USE m_temp, ONLY: weight_t_fine
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(OUT):: tpghost(yDim,xDim)
    DOUBLE PRECISION, INTENT(IN):: tp(yDim_thermal,xDim_thermal)
    INTEGER::y,x,counter,ii,jj
    DOUBLE PRECISION::w1,w2,w3,w4
    counter=0
    DO ii=1,xDim
        DO jj=1,yDim
            IF((jj.lt.yDim).and.(ii.lt.xDim))THEN
                counter=counter+1
                y=nint(weight_t_fine(counter,2))
                x=nint(weight_t_fine(counter,1))
                w1=weight_t_fine(counter,3)
                w2=weight_t_fine(counter,4)
                w3=weight_t_fine(counter,5)
                w4=weight_t_fine(counter,6)
                tpghost(jj,ii)=w1*tp(y,x)+w2*tp(y+1,x)+w3*tp(y,x+1)+w4*tp(y+1,x+1)
            ELSEIF((jj.eq.yDim).and.(ii.lt.xDim))THEN
                counter=counter+1
                y=nint(weight_t_fine(counter,2))
                x=nint(weight_t_fine(counter,1))
                w1=weight_t_fine(counter,3)
                w2=weight_t_fine(counter,4)
                tpghost(jj,ii)=w1*tp(y,x)+w2*tp(y,x+1)
            ELSEIF((jj.lt.yDim).and.(ii.eq.xDim))THEN
                counter=counter+1
                y=nint(weight_t_fine(counter,2))
                x=nint(weight_t_fine(counter,1))
                w1=weight_t_fine(counter,3)
                w2=weight_t_fine(counter,4)
                tpghost(jj,ii)=w1*tp(y,x)+w2*tp(y+1,x)
            ELSEIF((jj.eq.yDim).and.(ii.eq.xDim))THEN
                counter=counter+1
                y=nint(weight_t_fine(counter,2))
                x=nint(weight_t_fine(counter,1))
                tpghost(jj,ii)=tp(y,x)
            ENDIF
        END DO
    END DO
    END SUBROUTINE interpolate2d_usedin_fine_to_coarse





    SUBROUTINE interpolate2d_usedin_coarse_to_fine(tp,tpghost)
    USE simParam, ONLY:xDim_thermal, yDim_thermal,xDim, yDim
    USE m_temp, ONLY: weight_t_fine
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(OUT):: tpghost(yDim,xDim)
    DOUBLE PRECISION, INTENT(IN):: tp(yDim_thermal,xDim_thermal)
    INTEGER::y,x,counter,ii,jj
    DOUBLE PRECISION::w1,w2,w3,w4
    counter=0
    DO ii=1,xDim
        DO jj=1,yDim
            IF((jj.lt.yDim).and.(ii.lt.xDim))THEN
                counter=counter+1
                y=nint(weight_t_fine(counter,2))
                x=nint(weight_t_fine(counter,1))
                w1=weight_t_fine(counter,3)
                w2=weight_t_fine(counter,4)
                w3=weight_t_fine(counter,5)
                w4=weight_t_fine(counter,6)
                tpghost(jj,ii)=w1*tp(y,x)+w2*tp(y+1,x)+w3*tp(y,x+1)+w4*tp(y+1,x+1)
            ELSEIF((jj.eq.yDim).and.(ii.lt.xDim))THEN
                counter=counter+1
                y=nint(weight_t_fine(counter,2))
                x=nint(weight_t_fine(counter,1))
                w1=weight_t_fine(counter,3)
                w2=weight_t_fine(counter,4)
                tpghost(jj,ii)=w1*tp(y,x)+w2*tp(y,x+1)
            ELSEIF((jj.lt.yDim).and.(ii.eq.xDim))THEN
                counter=counter+1
                y=nint(weight_t_fine(counter,2))
                x=nint(weight_t_fine(counter,1))
                w1=weight_t_fine(counter,3)
                w2=weight_t_fine(counter,4)
                tpghost(jj,ii)=w1*tp(y,x)+w2*tp(y+1,x)
            ELSEIF((jj.eq.yDim).and.(ii.eq.xDim))THEN
                counter=counter+1
                y=nint(weight_t_fine(counter,2))
                x=nint(weight_t_fine(counter,1))
                tpghost(jj,ii)=tp(y,x)
            ENDIF
        END DO
    END DO
    END SUBROUTINE interpolate2d_usedin_coarse_to_fine

    SUBROUTINE findweight_usedin_fine_to_coarse
    USE simParam, ONLY:xDim_thermal, yDim_thermal,xDim, yDim
    USE m_temp, ONLY: weight_t_coarse,location_x,location_y,location_thermal_x,location_thermal_y
    IMPLICIT NONE
    INTEGER::y,x,ii,jj
    REAL*8::d1,d2,d3,d4
    REAL*8::npx,dummy,npy
    INTEGER::counter
    weight_t_coarse=0.0
    counter=0
    npy=real(yDim-1)/real(yDim_thermal-1)
    npx=real(xDim-1)/real(xDim_thermal-1)
    DO ii=1,xDim_thermal !coarse mesh
        DO jj=1,yDim_thermal !coarse mesh
            IF((jj.lt.yDim_thermal).and.(ii.lt.xDim_thermal))THEN
                counter=counter+1
                y=FLOOR(real(jj-1)*real(npy))+1 !fine mesh
                x=FLOOR(real(ii-1)*real(npx))+1 !fine mesh
                d1=((location_y(y)-location_thermal_y(jj))**2+(location_x(x)-location_thermal_x(ii))**2)**0.5
                d2=((location_y(y+1)-location_thermal_y(jj))**2+(location_x(x)-location_thermal_x(ii))**2)**0.5
                d3=((location_y(y)-location_thermal_y(jj))**2+(location_x(x+1)-location_thermal_x(ii))**2)**0.5
                d4=((location_y(y+1)-location_thermal_y(jj))**2+(location_x(x+1)-location_thermal_x(ii))**2)**0.5
                weight_t_coarse(counter,1)=x
                weight_t_coarse(counter,2)=y
                dummy=1.0/(1+d1/d2+d1/d3+d1/d4)
                weight_t_coarse(counter,3)=dummy
                weight_t_coarse(counter,4)=dummy*d1/d2
                weight_t_coarse(counter,5)=dummy*d1/d3
                weight_t_coarse(counter,6)=dummy*d1/d4
            ELSEIF((jj.eq.yDim_thermal).and.(ii.lt.xDim_thermal))THEN
                counter=counter+1
                y=FLOOR(real(jj-1)*real(npy))+1 !fine mesh
                x=FLOOR(real(ii-1)*real(npx))+1 !fine mesh
                d1=((location_y(y)-location_thermal_y(jj))**2+(location_x(x)-location_thermal_x(ii))**2)**0.5
                d2=((location_y(y)-location_thermal_y(jj))**2+(location_x(x+1)-location_thermal_x(ii))**2)**0.5
                weight_t_coarse(counter,1)=x
                weight_t_coarse(counter,2)=y
                dummy=1.0/(1+d1/d2)
                weight_t_coarse(counter,3)=dummy
                weight_t_coarse(counter,4)=dummy*d1/d2
                weight_t_coarse(counter,5)=0
                weight_t_coarse(counter,6)=0
            ELSEIF((jj.lt.yDim_thermal).and.(ii.eq.xDim_thermal))THEN
                counter=counter+1
                y=FLOOR(real(jj-1)*real(npy))+1 !fine mesh
                x=FLOOR(real(ii-1)*real(npx))+1 !fine mesh
                d1=((location_y(y)-location_thermal_y(jj))**2+(location_x(x)-location_thermal_x(ii))**2)**0.5
                d2=((location_y(y+1)-location_thermal_y(jj))**2+(location_x(x)-location_thermal_x(ii))**2)**0.5
                weight_t_coarse(counter,1)=x
                weight_t_coarse(counter,2)=y
                dummy=1.0/(1+d1/d2)
                weight_t_coarse(counter,3)=dummy
                weight_t_coarse(counter,4)=dummy*d1/d2
                weight_t_coarse(counter,5)=0
                weight_t_coarse(counter,6)=0
            ELSEIF((jj.eq.yDim_thermal).and.(ii.eq.xDim_thermal))THEN
                counter=counter+1
                y=FLOOR(real(jj-1)*real(npy))+1 !fine mesh
                x=FLOOR(real(ii-1)*real(npx))+1 !fine mesh
                weight_t_coarse(counter,1)=x
                weight_t_coarse(counter,2)=y
                weight_t_coarse(counter,3)=1
                weight_t_coarse(counter,4)=0
                weight_t_coarse(counter,5)=0
                weight_t_coarse(counter,6)=0
            ENDIF
        END DO
    END DO

    END SUBROUTINE findweight_usedin_fine_to_coarse




    SUBROUTINE interpolate2d_usedin_fine_to_coarse_vector(utot_fine,utot_coarse)
    USE simParam, ONLY:xDim_thermal, yDim_thermal,xDim, yDim
    USE m_temp, ONLY: weight_t_coarse
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN):: utot_fine(yDim,xDim,2)
    DOUBLE PRECISION, INTENT(OUT):: utot_coarse(yDim_thermal,xDim_thermal,2)
    INTEGER::y,x,counter,ii,jj
    DOUBLE PRECISION::w1,w2,w3,w4
    counter=0
    DO ii=1,xDim_thermal
        DO jj=1,yDim_thermal
            IF((jj.lt.yDim_thermal).and.(ii.lt.xDim_thermal))THEN
                counter=counter+1
                y=nint(weight_t_coarse(counter,2))
                x=nint(weight_t_coarse(counter,1))
                w1=weight_t_coarse(counter,3)
                w2=weight_t_coarse(counter,4)
                w3=weight_t_coarse(counter,5)
                w4=weight_t_coarse(counter,6)
                utot_coarse(jj,ii,:)=w1*utot_fine(y,x,:)+w2*utot_fine(y+1,x,:)+w3*utot_fine(y,x+1,:)+w4*utot_fine(y+1,x+1,:)
            ELSEIF((jj.eq.yDim_thermal).and.(ii.lt.xDim_thermal))THEN
                counter=counter+1
                y=nint(weight_t_coarse(counter,2))
                x=nint(weight_t_coarse(counter,1))
                w1=weight_t_coarse(counter,3)
                w2=weight_t_coarse(counter,4)
                utot_coarse(jj,ii,:)=w1*utot_fine(y,x,:)+w2*utot_fine(y,x+1,:)
            ELSEIF((jj.lt.yDim_thermal).and.(ii.eq.xDim_thermal))THEN
                counter=counter+1
                y=nint(weight_t_coarse(counter,2))
                x=nint(weight_t_coarse(counter,1))
                w1=weight_t_coarse(counter,3)
                w2=weight_t_coarse(counter,4)
                utot_coarse(jj,ii,:)=w1*utot_fine(y,x,:)+w2*utot_fine(y+1,x,:)
            ELSEIF((jj.eq.yDim_thermal).and.(ii.eq.xDim_thermal))THEN
                counter=counter+1
                y=nint(weight_t_coarse(counter,2))
                x=nint(weight_t_coarse(counter,1))
                utot_coarse(jj,ii,:)=utot_fine(y,x,:)
            ENDIF
        END DO
    END DO
    END SUBROUTINE interpolate2d_usedin_fine_to_coarse_vector