/*========================================*/
/*                 CBASE                  */
/*========================================*/
#ifndef ___CBASE___
#define ___CBASE___

#include "cgrid.h"

typedef double Param;
enum Num_Value { //...result values identification;
										ERR_VALUE = -1, //...results are absent;
										SPL_VALUE,
									  HEAT_VALUE,
									  FLUX_VALUE,
							 FLUX_COMPOS_VALUE,
							 HEAT_ANALYT_VALUE,
							 FLUX_ANALYT_VALUE,
									  LINK_VALUE,
									 DISPL_VALUE,
									 DILAT_VALUE,
							  DILAT_GRAD_VALUE,
						  DILAT_CLASSIC_VALUE,
									 CROSS_VALUE,
									 FRACT_VALUE,
									 TESTI_VALUE,
									ENERGY_VALUE,
									MOMENT_VALUE,
								 PRESSURE_VALUE,
								 VELOCITY_VALUE,
								 ROTATION_VALUE,
								POTENTIAL_VALUE,
								PROCESSOR_VALUE,
								PUREFRACT_VALUE,
									FLUX_R_VALUE,
									FLUX_X_VALUE,
									FLUX_Y_VALUE,
									FLUX_Z_VALUE,
									THIN_R_VALUE,
									THIN_X_VALUE,
									THIN_Y_VALUE,
									THIN_Z_VALUE,
								 MOMENT_R_VALUE,
								 MOMENT_X_VALUE,
								 MOMENT_Y_VALUE,
								 MOMENT_Z_VALUE,
								 STRESS_R_VALUE,
								 STRESS_X_VALUE,
								 STRESS_Y_VALUE,
								 STRESS_Z_VALUE,
								 NORMAL_R_VALUE,
								 NORMAL_X_VALUE,
								 NORMAL_Y_VALUE,
								 NORMAL_Z_VALUE,
						  NORMAL_R_GRAD_VALUE,
						  NORMAL_X_GRAD_VALUE,
						  NORMAL_Y_GRAD_VALUE,
						  NORMAL_Z_GRAD_VALUE,
					  NORMAL_R_CLASSIC_VALUE,
					  NORMAL_X_CLASSIC_VALUE,
					  NORMAL_Y_CLASSIC_VALUE,
					  NORMAL_Z_CLASSIC_VALUE,
			  NORMAL_R_GRAD_COHESION_VALUE,
			  NORMAL_X_GRAD_COHESION_VALUE,
			  NORMAL_Y_GRAD_COHESION_VALUE,
			  NORMAL_Z_GRAD_COHESION_VALUE,
								 STRAIN_R_VALUE,
								 STRAIN_X_VALUE,
								 STRAIN_Y_VALUE,
								 STRAIN_Z_VALUE,
								HEAT_ESHE_VALUE,
								FLUX_ESHE_VALUE,
							  DISPL_JUMP_VALUE,
							  DISPL_UNIX_VALUE,
							  DISPL_INDI_VALUE,
							  DISPL_ESHE_VALUE,
							  DISPL_GRAD_VALUE,
						  DISPL_CLASSIC_VALUE,
						 DISPL_COHESION_VALUE,
				  DISPL_GRAD_COHESION_VALUE,
						  STRESS_R_JUMP_VALUE,
						  STRESS_X_JUMP_VALUE,
						  STRESS_Y_JUMP_VALUE,
						  STRESS_Z_JUMP_VALUE,
						  STRESS_R_UNIX_VALUE,
						  STRESS_X_UNIX_VALUE,
						  STRESS_Y_UNIX_VALUE,
						  STRESS_Z_UNIX_VALUE,
						  STRESS_X_ESHE_VALUE,
						  STRESS_Y_ESHE_VALUE,
						  STRESS_Z_ESHE_VALUE,
						  STRESS_X_GRAD_VALUE,
						  STRESS_Y_GRAD_VALUE,
						  STRESS_Z_GRAD_VALUE,
						  STRESS_R_INDI_VALUE,
						  NORMAL_R_JUMP_VALUE,
						  NORMAL_X_JUMP_VALUE,
						  NORMAL_Y_JUMP_VALUE,
						  NORMAL_Z_JUMP_VALUE,
					  STRESS_R_CLASSIC_VALUE,
					  STRESS_X_CLASSIC_VALUE,
					  STRESS_Y_CLASSIC_VALUE,
					  STRESS_Z_CLASSIC_VALUE,
					 STRESS_R_COHESION_VALUE,
					 STRESS_X_COHESION_VALUE,
					 STRESS_Y_COHESION_VALUE,
					 STRESS_Z_COHESION_VALUE,
									ANALYT_VALUE = 1000,
									ANALYT2VALUE,
									ANALYT3VALUE,
									ANALYT4VALUE,
									ANALYT5VALUE,
									ANALYT6VALUE,
									ANALYT7VALUE,
									ANALYT8VALUE,
									ANALYT9VALUE,
							TESTI_ANALYT_VALUE,
					  SPL_ANALYT_RIGID_VALUE,
					PRESS_ANALYT_RIGID_VALUE,
				  FLUX_X_ANALYT_RIGID_VALUE,
				  FLUX_Z_ANALYT_RIGID_VALUE,
					 SPL_ANALYT_ABSORB_VALUE,
				  PRESS_ANALYT_ABSORB_VALUE,
				 FLUX_Z_ANALYT_ABSORB_VALUE,
	 HEAT_HOMOG_ALONG_LAYER_ANALYT_VALUE,
	 HEAT_HOMOG_ACROSSLAYER_ANALYT_VALUE,
	 FLUX_HOMOG_ALONG_LAYER_ANALYT_VALUE,
	 FLUX_HOMOG_ACROSSLAYER_ANALYT_VALUE,
};

/////////////////////////////////////////////////////////
//...bases class for serialization of multipoles objects;
class CBase {
protected:
		Param * param; //...data base parameters;
		static unsigned long COD; //...data base current code key;
public:
		virtual int  size_of_param() { return(0);}
		virtual void init() {
			delete_struct(param);
			param = new_struct<Param>(size_of_param());
		}
		void set_param(int k, Param value) { 
			if (0 <= k && k < size_of_param()) param[k] = value;
		}
		void get_param(int k, Param & value) { 
			if (0 <= k && k < size_of_param()) value = param[k];
		}
		Param get_param(int k) { 
			Param value(0.); get_param(k, value); 
			return(value);
		}
		virtual void release() {	
			delete_struct(param);
		}
public:
static int NUM_MPLS, NUM_QUAD, NUM_TIME, NUM_LOCAL, NUM_ADHES, NUM_VIBRO, NUM_GEOMT, NUM_PHASE, MAX_PHASE, BOX_LINK_PERIOD;
////////////////////////////////////////
//...filling block structure parameters;
		virtual void set_fasa_hmg(double CC){}
		virtual void set_fasa_hmg(double CC, double CC2){}
		virtual void set_fasa_hmg(double Hz, double Ro1, double C1){}
		virtual void set_fasa_hmg(double nju1, double nju2, double G1, double G2){}
		virtual void set_fasa_hmg(double Hz, double Ro1, double Ro2, double C1, double C2){}
		virtual void set_fasa_hmg(double nju1, double nju2, double G1, double G2, double C1, double C2){}
		virtual void set_fasa_hmg(double GM_re, double GM_im, double KM_re, double KM_im, double GI_re, double GI_im, double KI_re, double KI_im){}
		virtual void set_fasa_hmg(double R1, double R2, double nju1, double nju2, double nju3, double G1, double G2, double G3, double alpha){}
		virtual void set_fasa_hmg(double nju1, double nju2, double G1, double G2, double l11, double l12, double l21, double l22, double A, double B){}
		virtual void set_fasa_hmg(double K1, double K2, double G1, double G2, double l1, double l2, double d1, double d2, double A, double B, double d0){}
		virtual void set_fasa_hmg(double GM_re, double GM_im, double KM_re, double KM_im, double GI_re, double GI_im, double KI_re, double KI_im, double GL_re, double GL_im, double KL_re, double KL_im){}
		virtual void set_time	 (double d_time = 0.1, double end_time = 1.) {}
		//	set_param(NUM_TIME+1, 0.); set_param(NUM_TIME+2, d_time); set_param(NUM_TIME+3, end_time);
		//}
		virtual void set_adhesion(double AA, double BB){}
		virtual void set_geometry(double rad, double layer = 0.) { set_param(NUM_GEOMT, rad); set_param(NUM_GEOMT+1, rad+layer); }
		virtual void set_mpls    (double mpls) {set_param(NUM_MPLS, mpls);}
		virtual void set_quad	 (double quad) {set_param(NUM_QUAD, quad);}
		virtual void set_vibro   (double Hz)   {set_param(NUM_VIBRO,  Hz);}
		virtual void set_lattice (double lattice)  {set_param(NUM_TIME, lattice);}
		virtual void set_normaliz(double normaliz) {set_param(NUM_MPLS+1,		normaliz);}
		virtual void set_lagrange(double lagrange) {set_param(size_of_param()-1, lagrange);}
		virtual void set_local	 (double local)	 {set_param(NUM_LOCAL, local);}
////////////////////////////
//...extracting sample data;
		virtual double get_Young(int num_phase){ return(0.);}
		virtual double get_Shear(int num_phase){ return(0.);}
		virtual double get_Poisn(int num_phase){ return(0.);}
		virtual double get_geometry(int m = 0) {
				if (size_of_param() <= NUM_GEOMT+m) return(0.);
				else return(param[NUM_GEOMT+m]); 
		}
////////////////////////////
//...sample parametrization;
		virtual Param * ID_EPS()        { return(param);}
		virtual Param * ID_CUSTOMER()   { return(param);}
		virtual Param * ID_MAIN()       { return(param);}
		virtual Param * ID_YOUNG_MODUL(){ return(param);}
		virtual Param * ID_INERTIA_MOM(){ return(param);}
		virtual Param * ID_TORSION_MOM(){ return(param);}
		virtual Param * ID_ENERGY_MOM() { return(param);}
public:
//... output format for Golden Surfer;
		virtual void GetSurferFormat(char * SURF_FILE, CGrid * nd, Num_Value _FMF = ERR_VALUE, int id_variant = 0, int id_axis = AXIS_Z, int iparam = 0) {};
		virtual void GetSurferFormat(const char * SURF_FILE, CGrid * nd, Num_Value _FMF = ERR_VALUE, int id_variant = 0, int id_axis = AXIS_Z, int iparam = 0) {};
		virtual void GetDataFormat(char * DATA_FILE, CGrid * nd, Num_Value _FMF = ERR_VALUE, int id_variant = 0, int id_axis = AXIS_Z, int iparam = 0) {};
		virtual void GetDataFormat(const char * DATA_FILE, CGrid * nd, Num_Value _FMF = ERR_VALUE, int id_variant = 0, int id_axis = AXIS_Z, int iparam = 0) {};
//... output format on the block structure for CSV tables;
		virtual void GetCsvFormat(char * CSV_FILE, CGrid * nd, int id_variant = 0, int id_centroid = 0, CGrid * bnd = NULL) {};
		virtual void GetCsvFormat(const char * CSV_FILE, CGrid * nd, int id_variant = 0, int id_centroid = 0, CGrid * bnd = NULL) {};
//...analytical models;
		virtual void rgradf_matrix(double ** T, int n, int shift_m = 0, int shift_n = 0, double f = 1.){};
		virtual void rdivrf_matrix(double ** T, int n, int shift_m = 0, int shift_n = 0, double f = 1.){};
		virtual void unit_f_matrix(double ** T, int n, int shift_m = 0, int shift_n = 0, double f = 1.){};
		virtual void toreal_matrix(double ** T, int n, int shift_m = 0, int shift_n = 0){};
		virtual void grdivf_forwrd(double ** T, int n, int shift_m = 0, int shift_n = 0, double f = 1.){};
		virtual void toreal_forwrd(double ** T, int n, int shift_m = 0, int shift_n = 0){};
		virtual void grdivf_transf(double ** T, int n, int shift_m = 0, int shift_n = 0, double f = 1.){};
		virtual void toreal_transf(double ** T, int n, int shift_m = 0, int shift_n = 0){};
		virtual double  TakeLayer				(double ff)						  { return 0.;}
		virtual double  TakeLayer_E1			(double ff)						  { return 0.;}
		virtual double  TakeLayer_E2			(double ff)						  { return 0.;}
		virtual double  TakeLayer_G1			(double ff)						  { return 0.;}
		virtual double  TakeLayer_G2			(double ff)						  { return 0.;}
		virtual double  TakeCylinder			(double ff, double eps = EE) { return 0.;}
		virtual double  TakeEshelby				(double ff, double ff_l)  { return 0.;}
		virtual double  TakeEshelby_two		(double ff)						  { return 0.;}
		virtual double  TakeEshelby_grad		(double ff)						  { return 0.;}
		virtual double  TakeEshelby_volm		(double ff, double ff_l)	  { return 0.;}
		virtual double  TakeEshelby_volm_two (double ff)					  { return 0.;}
		virtual double  TakeEshelby_volm_sym (double ff)					  { return 0.;}
		virtual double  TakeEshelby_shear		(double ff, double ff_l, double eps = EE, int max_iter = 100) { return 0.;}
		virtual double  TakeEshelby_shear_two(double ff, double eps = EE, int max_iter = 100)				  { return 0.;}
		virtual double  TakeEshelby_shear_sym(double ff, double eps = EE, int max_iter = 100)				  { return 0.;}
		virtual double  TakeEshelby_shear_det(double ff, double alpha = -1)										  { return 0.;}
		virtual double  TakeLayer_kk(int N, double * ff, double * kk)												  { return 0.;}
		virtual double  TakeLayer_kk(int N, double * ff, double * kk, double * mu)								  { return 0.;}
		virtual double  TakeLayer_kk(int N, double * ff, double * kk, double * mu, double * nj)			  { return 0.;}
		virtual double  TakeLayer_GG(int N, double * ff, double * mu, double * nj)								  { return 0.;}
		virtual double  TakeLayer_GG(int N, double * ff, double * kk, double * mu, double eps = EE, int max_iter = 100) { return 0.;}
		virtual double  TakeLayer_GG(int N, double * ff, double * kp, double * mu, double * nj, double eps = EE, int max_iter = 100) { return 0.;}
		virtual double  TakeCrack_two(double X, double Y, double kappa, double * AA, double * BB, double * AD, double * BD, int max_iter = 100) { return 0.; }
		virtual complex TakeEshelby_visco_volm_two (double ff)													{ return comp(0.);}
		virtual complex TakeEshelby_visco_shear_two(double ff, double eps = EE, int max_iter = 100)  { return comp(0.);}
//...analytical models with blocks installing;
		virtual void TakeEshelbyModel	 (double ff, double ff_l, double fK = 0.5, double fG = 0.5){};
		virtual void TakeEshelbyModel_two(double ff, double fK = 0.5, double fG = 0.5){};
		virtual void TakeEshelbyGradModel	 (double ff, double ff_l, double fK = 0.5, double fG = 0.5){};
		virtual void TakeEshelbyGradModel_two(double ff, double fK = 0.5, double fG = 0.5){};
		virtual void TakeEshelbyGradIncluModel(double ff, double ff_l, double fK = 0.5, double fG = 0.5){};
//...recomputing of the temperature pole accordance to stabilization scheme;
		virtual void TakeStabStep(double * Temp, int NN, double alpha){};
		virtual void TakeStabStep_layer(double * Temp, int N_SC, int N_CU, int N_cells, double alpha){};
public:
//...constructor/destructor;
		CBase() {
			param = NULL;
		}
      virtual ~CBase(void)	{
			release();
		}
};
#endif
