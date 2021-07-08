    !	 ========================================================
    !	 This part code implementing a flow past a cylinder
    !	 ========================================================
    MODULE m_temp

    IMPLICIT NONE

    DOUBLE PRECISION, ALLOCATABLE:: ft(:,:,:), fEqt(:,:,:), tp(:,:),tpghost(:,:)
    DOUBLE PRECISION, ALLOCATABLE:: location_x(:),location_y(:),location_thermal_x(:),location_thermal_y(:),weight_t_fine(:,:)
    DOUBLE PRECISION, ALLOCATABLE::weight_t_coarse(:,:),uSqr_t(:,:)
	DOUBLE PRECISION, ALLOCATABLE::nu_left(:),nu_right(:)

    END MODULE m_temp