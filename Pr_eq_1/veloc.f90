	! All subroutines related to fluid flow equation

    !---------------------------------------------------------------------
    SUBROUTINE veloc_init
    !---------------------------------------------------------------------
    ! initialize flow variables
    ! last modification: 03-Mar-2010
    USE m_veloc
    USE simParam, ONLY: xDim, yDim,rho0
    IMPLICIT NONE
    INTEGER::y

    ALLOCATE(fu(yDim,xDim,9))
    ALLOCATE(fEqu(yDim,xDim,9))
    ALLOCATE(uTot(yDim,xDim,2))
    ALLOCATE(uSqr(yDim,xDim))
    ALLOCATE(rho(yDim,xDim))
	ALLOCATE(solid(yDim,xDim))
	ALLOCATE(fyexternal(yDim,xDim))

    DO y = 1, yDim
        uTot(y,:,1) = 0.0d0
        uTot(y,:,2) = 0.0d0
    END DO
    rho  = rho0
    uSqr = uTot(:,:,1) * uTot(:,:,1) + uTot(:,:,2) * uTot(:,:,2)

    CALL computeFeq(fEqu,rho,uTot,uSqr)

    fu = fEqu

    END SUBROUTINE veloc_init

    !	 ========================================================
    !	 The main scheme of implementing Velocity field in LB
    !	 ========================================================

    SUBROUTINE veloc_main

    USE m_veloc
    USE simParam, ONLY: xDim, yDim, tau,gforce
    USE m_temp, ONLY:tpghost

    IMPLICIT NONE
    INTEGER::x,y
    CALL stream(fu)
	CALL computeMacros(fu,rho,uTot,uSqr)
	DO x = 1, xDim
		DO y = 1, yDim
			fyexternal(y,x)=-gforce*rho(y,x)*(tpghost(y,x)-0.5d0)
		END DO
	END DO
    CALL computeFeq(fEqu,rho,uTot,uSqr)

    CALL collide(fu,fEqu,tau)
    CALL bounceback(fu)
     

    !CALL bounceback(fu,ff,image) 


    END SUBROUTINE veloc_main


	! Streaming subroutine
    SUBROUTINE stream(fu)
    USE simParam
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(INOUT):: fu(yDim,xDim,0:8)
	DOUBLE PRECISION f_hlp(yDim,xDim,0:8)
    INTEGER :: x,y,y_n,x_e,y_s,x_w

	f_hlp=fu;
    DO x = 1, xDim
        DO y = 1, yDim

            y_n = MOD(y,yDim) + 1
            x_e = MOD(x,xDim) + 1
            y_s = yDim-MOD(yDim + 1-y, yDim)
            x_w = xDim-MOD(xDim + 1-x, xDim)

            fu(y,x,0) = f_hlp(y,x,0)       !.........zero: just copy
            fu(y,x_e,1) = f_hlp(y,x,1)     !.........east
            fu(y_n,x,2) = f_hlp(y,x,2)     !.........north
            fu(y,x_w,3) = f_hlp(y,x,3)     !.........west
            fu(y_s,x,4) = f_hlp(y,x,4)     !.........south
            fu(y_n,x_e,5) = f_hlp(y,x,5)   !.........north-east
            fu(y_n,x_w,6) = f_hlp(y,x,6)   !.........north-west
            fu(y_s,x_w,7) = f_hlp(y,x,7)   !.........south-west
            fu(y_s,x_e,8) = f_hlp(y,x,8)   !.........south-east   

        enddo
    enddo

    END SUBROUTINE stream

    !	 ========================================================
    !	 Implement Bounce-back on solid==1 
    !	 ========================================================
    SUBROUTINE bounceback(fu)
    USE D2Q9Const, ONLY: opposite
    USE simParam, ONLY: xDim, yDim
    USE m_veloc, ONLY: solid

    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(INOUT):: fu(yDim,xDim,0:8)
    DOUBLE PRECISION::f_help(yDim,xDim,0:8)

    INTEGER:: i,x,y
    f_help=fu;
    DO x = 1 , xDim 
        DO y = 1,yDim 

            ! finish judging

            IF (solid(y,x).eq.1) then   !bounce back at solid-liquid interface
                DO i = 0, 8
                    fu(y,x,i) = f_help(y,x,opposite(i))
                END DO
            end IF

        END DO
    END DO

    END SUBROUTINE bounceback


    !	 ========================================================
    !	 Compute macroscopic variables density and velocity from distribution functions
    !	 ========================================================
    SUBROUTINE computeMacros(ff,rho,uTot,uSqr)

    USE simParam, ONLY: xDIm, yDim, rho0
    USE m_veloc, ONLY: solid

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(INOUT):: ff(yDim,xDim,0:8)
    DOUBLE PRECISION, INTENT(INOUT):: uTot(yDim,xDim,2), rho(yDim, xDim), uSqr(yDim, xDim)
    INTEGER:: x,y

    !---------------------------- this is for liquid when fs<1.0. If fs>=1.0, u(),v() are zero and rho=rho0?
    DO y = 1, yDim
        DO x = 1, xDim
            IF(solid(y,x).lt.1.0)then 
                rho(y,x)  = ff(y,x,0) + ff(y,x,1) + ff(y,x,2) + ff(y,x,3) + ff(y,x,4) + ff(y,x,5) + ff(y,x,6) + ff(y,x,7) + ff(y,x,8)
                if (rho(y,x) .ne. rho(y,x) ) STOP 'the problem diverged'
                uTot(y,x,1)  = (ff(y,x,1) - ff(y,x,3) + ff(y,x,5) - ff(y,x,6) - ff(y,x,7) + ff(y,x,8)) / rho(y,x)
                uTot(y,x,2)  = (ff(y,x,2) - ff(y,x,4) + ff(y,x,5) + ff(y,x,6) - ff(y,x,7) - ff(y,x,8)) / rho(y,x)
                uSqr(y,x) = uTot(y,x,1) * uTot(y,x,1) + uTot(y,x,2) * uTot(y,x,2)
            ELSE
                rho(y,x) = rho0
                uTot(y,x,1) = 0.0d0
                uTot(y,x,2) = 0.0d0
                uSqr(y,x) = 0.0d0
            endif
        END DO

        uTot(y,1,1) = 0.0  !----------------------------------------?????
        uTot(y,xDim,2) = 0.0
    END DO
    END SUBROUTINE computeMacros



    !	 ========================================================
    !	 Compute equilibrium distribution
    !	 ========================================================
    SUBROUTINE computeFeq(fEqu,rho,uTot,uSqr)
    USE D2Q9COnst, ONLY: t, v
    USE simParam, ONLY: xDim, yDim
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN):: rho(yDim,xDim), uSqr(yDim,xDim), uTot(yDim,xDim,2)
    DOUBLE PRECISION, INTENT(INOUT):: fEqu(yDim,xDim,0:8)
    INTEGER:: i, x, y
    DOUBLE PRECISION:: uxy

    DO i = 0, 8
        DO x = 1, xDim
            DO y = 1, yDim
                uxy = uTot(y,x,1) * v(i,1) + uTot(y,x,2) * v(i,2)
                fEqu(y,x,i) = t(i) * rho(y,x) * (1.0d0 + 3.0d0 * uxy + 4.5d0 * uxy * uxy - 1.5d0 * uSqr(y,x))
            END DO
        END DO
    END DO
    END SUBROUTINE computeFeq




    !	 ========================================================
    !	 LBGK collision step
    !	 ========================================================
    SUBROUTINE collide(fu,fEqu,tau)

    USE simParam, ONLY: xDim, yDim
    USE m_veloc, ONLY: solid,fyexternal
	USE D2Q9Const

    IMPLICIT NONE


    DOUBLE PRECISION, INTENT(IN):: fEqu(yDim,xDim,0:8)
    DOUBLE PRECISION, INTENT(INOUT):: fu(yDim,xDim,0:8)
    DOUBLE PRECISION :: tau
	DOUBLE PRECISION::fy_local(yDim,xDim,0:8)
    INTEGER:: x,y,i

	fy_local=0.0d0;
    DO x = 1, xDim
        DO y = 1, yDim
            DO i = 0, 8
				IF (i.ne.0) fy_local(y,x,i)=3.d0*t(i)*fyexternal(y,x)*v(i,2)
                IF(solid(y,x).eq.0) fu(y,x,i) = (1.0d0 - 1.0d0/tau) * fu(y,x,i) + 1.0d0/tau * fequ(y,x,i)+fy_local(y,x,i)
            END DO
        END DO
    END DO

    END SUBROUTINE collide





	! Assign solid nodes

    SUBROUTINE soliddeff(solid)
    USE simParam, ONLY: xDim, yDim
    IMPLICIT NONE
    INTEGER:: solid(yDim,xDim) ! a matrix that define fluid/solid interface flid=0 interface=1 solid=2
    INTEGER::x,y
    solid=0
    DO x = 1,xDim
        DO y = 1 , yDim

            IF((x==1).or.(x==xDim).or.(y==1).or.(y==yDim))THEN
                solid(y,x)=1
            ENDIF
        ENDDO
    ENDDO

    END SUBROUTINE soliddeff

