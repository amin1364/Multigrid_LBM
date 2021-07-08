    !	 ========================================================
    !	 Basic Simulation Parameter/Input Parameters
    !	 ========================================================
    MODULE simParam
    IMPLICIT NONE
    INTEGER:: xDim, yDim,yDim_thermal, xDim_thermal, tMax, kPlot,n_ratio,n_ratio_u
    INTEGER:: tp_u, cl_u
    DOUBLE PRECISION:: uchar_lb_fine,uchar_lb_coarse,uchar_lb
    DOUBLE PRECISION:: U_ratio
    DOUBLE PRECISION::rho0
    DOUBLE PRECISION::tau,Ra_no,tau_t,tau_coarse,tau_fine,tau_t_coarse,tau_t_fine
    DOUBLE PRECISION::Length_ref,Length_ref_temp,gforce
    DOUBLE PRECISION::visc_lb_fine,visc_lb_coarse
    DOUBLE PRECISION::thermal_diff_fine,thermal_diff_coarse
    DOUBLE PRECISION::Pr_no
    CONTAINS

    SUBROUTINE ComputeParam
    USE D2Q9Const,ONLY:c_squ
    IMPLICIT NONE

	U_ratio=1;
	n_ratio_u=1;
	n_ratio=1;
	
    tp_u=1;
    tMax=100000*20;
    kplot=100000
    Pr_no=0.01d0;
    rho0=1.0d0;

    Ra_no=1.0e6;
	OPEN(1021,FILE='check_input.txt');

	xDim=100;
	yDim=100;
	xDim_thermal=100;
	yDim_thermal=100;
	uchar_lb_fine=0.05d0;
	Length_ref=DBLE(xDim);
	Length_ref_temp=DBLE(xDim_thermal);
	visc_lb_fine=uchar_lb_fine*Length_ref*SQRT(Pr_no/Ra_no);
	thermal_diff_fine=visc_lb_fine/Pr_no;
	gforce = -uchar_lb_fine*uchar_lb_fine/Length_ref;	
	tau=3.0d0*visc_lb_fine+0.5;
	visc_lb_coarse=visc_lb_fine;	
	thermal_diff_coarse=thermal_diff_fine;
	tau_t=3.0d0*thermal_diff_coarse+0.5d0;
	U_ratio=1.0
	uchar_lb=uchar_lb_fine;

	tau_coarse=3.0d0*visc_lb_coarse+0.5d0;
	tau_fine=3.0d0*visc_lb_fine+0.5d0;
	tau_t_coarse=3.0d0*thermal_diff_coarse+0.5d0;
	tau_t_fine=3.0d0*thermal_diff_fine+0.5d0;
	uchar_lb_coarse=uchar_lb_fine;
	WRITE(1021,*)'Pr_no',Pr_no
	WRITE(1021,*)'Ra_no',Ra_no
	WRITE(1021,*)'xDim',xDim
	WRITE(1021,*)'yDim',yDim
	WRITE(1021,*)'xDim_thermal',xDim_thermal
	WRITE(1021,*)'xDim_thermal',xDim_thermal
	WRITE(1021,*)'tauf',tau
	WRITE(1021,*)'taug',tau_t
	WRITE(1021,*)'uchar_lb',uchar_lb






	
	CLOSE(1021)



    END SUBROUTINE ComputeParam
    END MODULE simParam
