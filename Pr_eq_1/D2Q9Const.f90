    !	 ========================================================
    !	 Lattice constants for the D2Q9 lattice
    !	 ========================================================
    MODULE D2Q9Const
    IMPLICIT NONE
    !	 D2Q9 Weights
    DOUBLE PRECISION,PARAMETER:: t(0:8) = (/4.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0,1.0d0/9.0d0&
        &,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0,1.0d0/36.0d0/)
    !	D2Q9 Directions
    INTEGER:: v(0:8,2) = (/(/0,1,0,-1,0,1,-1,-1,1/),(/0,0,1,0,-1,1,1,-1,-1/)/)

    INTEGER, PARAMETER:: opposite(0:8) = (/0,3,4,1,2,7,8,5,6/)
    DOUBLE PRECISION, PARAMETER::c_squ=1.0d0/3.0d0
    DOUBLE PRECISION::con2
    CONTAINS
    SUBROUTINE calculatecons ()
    IMPLICIT NONE
    con2=(2.d0 * c_squ *c_squ)


    END SUBROUTINE calculatecons
    END MODULE D2Q9Const
