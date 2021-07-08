    !	 ========================================================
    !	 This part code implementing a flow past a cylinder
    !	 ========================================================
    MODULE m_veloc

    IMPLICIT NONE

    DOUBLE PRECISION, ALLOCATABLE:: fu(:,:,:), fEqu(:,:,:), uTot(:,:,:)
    DOUBLE PRECISION, ALLOCATABLE:: rho(:,:), uSqr(:,:)
    DOUBLE PRECISION, ALLOCATABLE::fyexternal(:,:)
    INTEGER, ALLOCATABLE::solid(:,:)

    END MODULE m_veloc