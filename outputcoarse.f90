    !	 ========================================================
    !	 Write the temperature component for a coarse mesh in the output file.
    !	 ========================================================
    SUBROUTINE outputcoarse(tp,tStep)

    USE simParam, ONLY: yDim_thermal, xDim_thermal
    IMPLICIT NONE

    INTEGER, INTENT(IN):: tStep

    DOUBLE PRECISION, INTENT(IN):: tp(yDim_thermal,xDim_thermal)
    CHARACTER(len=24) :: FN
    INTEGER:: x,y, tDim, i
    DOUBLE PRECISION, ALLOCATABLE:: pu(:)

    ALLOCATE(pu(yDim_thermal*xDim_thermal))





    WRITE (FN, "(A4,I8,A8)") "temp",tstep,"sec.tec"
    OPEN(250,FILE=FN)
    WRITE(250,1021)tstep
    WRITE(250,1022)xDim_thermal,yDim_thermal

1021 format( 'TITLE = "Time = ',I10,'"',/, &
    'VARIABLES = "X", "Y","T"')
1022 format( 'ZONE I=',i6,', J=',i6,', F=BLOCK')
    tDim=xDim_thermal*yDim_thermal
    DO y=1, yDim_thermal
        DO x=1, xDim_thermal
            i=x+(y-1)*xDim_thermal
            pu(i)=x
        END DO
    END DO
    WRITE(250,310)(int(pu(i)),i=1,tDim)
    DO y=1, yDim_thermal
        DO x=1, xDim_thermal
            i=x+(y-1)*xDim_thermal
            pu(i)=y
        END DO
    END DO
    WRITE(250,310)(int(pu(i)),i=1,tDim)




    DO y=1, yDim_thermal
        DO x=1, xDim_thermal
            i=x+(y-1)*xDim_thermal
            pu(i)=tp(y,x)
        END DO
    END DO
    WRITE(250,320)(pu(i),i=1,tDim)






310 format(5I10)
320 format(5f20.10)
    CLOSE(250)


    RETURN

    END SUBROUTINE outputcoarse