
    !	 Write the result file in a .tec format that can be viewed by Paraview or Tecplot.
    SUBROUTINE writeoutput(tp,u,tStep)

    USE simParam, ONLY: xDim, yDim,uchar_lb
    IMPLICIT NONE

    INTEGER, INTENT(IN):: tStep
    DOUBLE PRECISION, INTENT(IN):: u(yDim,xDim,1:2)
    DOUBLE PRECISION, INTENT(IN)::  tp(yDim,xDim)
    CHARACTER(len=24) :: FN
    INTEGER:: x,y, tDim, i
    DOUBLE PRECISION, ALLOCATABLE:: pu(:)
    ALLOCATE(pu(yDim*xDim))

     WRITE (FN, "(A4,I8,A8)") "re",tstep," sec.tec"
    OPEN(250,FILE=FN)
    WRITE(250,1021)tstep
    WRITE(250,1022)xDim,yDim

1021 format( 'TITLE = "Time = ',I10,'"',/, &
        'VARIABLES = "X", "Y" ,"U Vel", "V Vel" , "T"')

1022 format( 'ZONE I=',i6,', J=',i6,', F=BLOCK')
    tDim=xDim*yDim
    DO y=1, yDim
        DO x=1, xDim
            i=x+(y-1)*xDim
            pu(i)=x
        END DO
    END DO
    WRITE(250,310)(int(pu(i)),i=1,tDim)
    DO y=1, yDim
        DO x=1, xDim
            i=x+(y-1)*xDim
            pu(i)=y
        END DO
    END DO
    WRITE(250,310)(int(pu(i)),i=1,tDim)


    DO y=1, yDim
        DO x=1, xDim
            i=x+(y-1)*xDim
            pu(i)=u(y,x,1)/uchar_lb

        END DO
    END DO
    WRITE(250,320)(pu(i),i=1,tDim)

    DO y=1, yDim
        DO x=1, xDim
            i=x+(y-1)*xDim
            pu(i)=u(y,x,2)/uchar_lb
        END DO
    END DO
    WRITE(250,320)(pu(i),i=1,tDim)








    DO y=1, yDim
        DO x=1, xDim
            i=x+(y-1)*xDim
            pu(i)=tp(y,x)
        END DO
    END DO
    WRITE(250,320)(pu(i),i=1,tDim)





310 format(1000I10)
320 format(1000f20.10)
    CLOSE(250)


    RETURN

    END SUBROUTINE writeoutput
