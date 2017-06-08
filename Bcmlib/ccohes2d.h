/*===========================================*/
/*                  CCOHES2D                 */
/*===========================================*/
#ifndef ___CCohes2D___
#define ___CCohes2D___

#include "ccomput2d.h"

////////////////////////////////////////////////////////
//...class of blocks partition for double plane problem;
class CCohes2D : public CComput2D<double> {
public:
		int size_of_param() { return(26);}
public:
		void init() {
			delete_struct(param);
			param = new_struct<Param>(size_of_param());
		}
public:
//...constructor;
		CCohes2D (int num_phase = 8) {
			NUM_PHASE = num_phase;
			MAX_PHASE = 4;
			init();
		};
protected:
double KK[11], GG[22];
static int NUM_SHEAR, NUM_SHIFT, NUM_ADHES, MAX_PHASE, NUM_HESS, regul;
		int  block_shape_init(Block<double> & B, Num_State id_free);
//...auxilliary operations with block matrix;
		void jump1_common_x	(double * P, int i, int m);
		void jump1_compress_x(double * P, int i, int m);
		void jump1_classic_x (double * P, int i, int m);
		void jump1_cohesion_x(double * P, int i, int m);
		void jump1_common_y	(double * P, int i, int m);
		void jump1_compress_y(double * P, int i, int m);
		void jump1_classic_y (double * P, int i, int m);
		void jump1_cohesion_y(double * P, int i, int m);
		void jump2_common_x	(double * P, int i, int m);
		void jump2_compress_x(double * P, int i, int m);
		void jump2_classic_x	(double * P, int i, int m);
		void jump2_cohesion_x(double * P, int i, int m);
		void jump2_common_y	(double * P, int i, int m);
		void jump2_compress_y(double * P, int i, int m);
		void jump2_classic_y	(double * P, int i, int m);
		void jump2_cohesion_y(double * P, int i, int m);
		void jump4_common_x	(double * P, int i, int m);
		void jump4_compress_x(double * P, int i, int m);
		void jump4_classic_x (double * P, int i, int m);
		void jump4_common_y	(double * P, int i, int m);
		void jump4_compress_y(double * P, int i, int m);
		void jump4_classic_y (double * P, int i, int m);
		void hessian_deriv_N (int k, double * P, int i, int id_dop);
		void hessian_deriv_N (int k, double * P, int i) {
			  hessian_deriv_N (k,   P, i, 0);
			  hessian_deriv_N (k+1, P, i, 1);
		}
		void jump_admittance (int l, int i, int m, double adm_re, int k = -1, double adm_im = 0.);
		void jump_make_local (int i, int m);
		void jump_make_common(int i, int m);
//...forming block matrix elements;
		Num_State gram1    (CGrid * nd, int i, int id_local){ return gram4(nd, i, id_local);};
		Num_State gram2    (CGrid * nd, int i, int id_local){ return gram3(nd, i, id_local);};
		Num_State gram3    (CGrid * nd, int i, int id_local);
		Num_State gram4    (CGrid * nd, int i, int id_local);
		Num_State transfer1(CGrid * nd, int i, int k, int id_local){ return transfer4(nd, i, k, id_local);};
		Num_State transfer2(CGrid * nd, int i, int k, int id_local){ return transfer3(nd, i, k, id_local);};
		Num_State transfer3(CGrid * nd, int i, int k, int id_local);
		Num_State transfer4(CGrid * nd, int i, int k, int id_local);
		Num_State rigidy1  (CGrid * nd, int i, double * K);
		Num_State computing_header(Num_Comput Num);
public:
//...параметры задачи;
		void set_fasa_hmg(double nju1, double nju2, double G1, double G2, double C1, double C2);
		void set_fasa_hmg(double nju1, double nju2, double nju3, double G1, double G2, double G3, double C1, double C2, double C3);
//...извлечение параметров задачи;
		double get_Young(int num_phase) {
			if (size_of_param() <= NUM_SHEAR+3+NUM_SHIFT*num_phase) return(0.);
			else return(param[NUM_SHEAR+NUM_SHIFT*num_phase]*2.*(1.+param[NUM_SHEAR+1+NUM_SHIFT*num_phase])); 
		}
		double get_Shear(int num_phase) {
			if (size_of_param() <= NUM_SHEAR+3+NUM_SHIFT*num_phase) return(0.);
			else return(param[NUM_SHEAR+NUM_SHIFT*num_phase]); 
		}
		double get_Poisn(int num_phase) {
			if (size_of_param() <= NUM_SHEAR+3+NUM_SHIFT*num_phase) return(0.);
			else return(param[NUM_SHEAR+1+NUM_SHIFT*num_phase]); 
		}
//...результаты решения задачи;
		void GetFuncAllValues(double X, double Y, double Z, double * F, int id_block, Num_Value id_F, int id_variant = 0, int iparam = 0);
		void GetEnergyValue	(int k, double * energy);
//...аналитические модели с межфазным слоем и адгезией;
		void	 TakeLayerModel(double L, double H, double l, double nju1, double nju2, double G1, double G2, double l1, double l2, double A = 0.);
		double TakeLayer_E1	(double ff);
		double TakeLayer_G1	(double ff);
//...прямые аналитические методы;
		double TakeEshelby_volm		 (double ff, double ff_l);
		double TakeEshelby_volm_two (double ff);
		double TakeEshelby_shear	 (double ff, double ff_l,  double eps = EE, int max_iter = 100);
		double TakeEshelby_shear_two(double ff, double & det, double eps = EE, int max_iter = 100);
//...четырехфазная модель с градиентным включением и слоем;
		double TakeEshelbyGradInclu_volm (double ff, double ff_l);
		double TakeEshelbyGradInclu_shear(double ff, double ff_l,  double eps = EE, int max_iter = 100);
//...трехфазная модель с градиентной матрицей и классическим включением;
		double TakeGradMatrix_k0(double ff);
		double TakeGradMatrix_k1(double ff);
		double TakeGradMatrix_G1(double ff, double & det, double eps = EE, int max_iter = 100);
		double TakeGradMatrix_G1_classic(double ff, double & det, double eps = EE, int max_iter = 100);
		//double TakeGradMatrix_k2(double ff);
		//double TakeGradMatrix_G2(double ff);
		//double TakeGradMatrix_G2_simple(double ff);
//...четырехфазная модель с градиентным слоем;
		double TakeGradLayer_k0(double ff, double ff_l);
		double TakeGradLayer_k1(double ff, double ff_l);
		double TakeGradLayer_G1(double ff, double ff_l, double eps = EE, int max_iter = 100);
		double TakeGradLayer_k2(double ff, double ff_l);
		double TakeGradLayer_G2(double ff, double ff_l);
		double TakeGradLayer_G2_simple(double ff, double ff_l);
//...аналитические модели с установкой блоков;
		void TakeEshelbyModel		   (double ff, double ff_l, double fK = 0.5, double fG = 0.5);
		void TakeEshelbyModel_two	   (double ff, double fK = 0.5, double fG = 0.5);
		void TakeEshelbyGradModel	   (double ff, double ff_l, double fK = 0.5, double fG = 0.5);
		void TakeEshelbyGradModel_two (double ff, double fK = 0.5, double fG = 0.5);
		void TakeEshelbyGradIncluModel(double ff, double ff_l, double fK = 0.5, double fG = 0.5);
//...конечная трещина;
		double TakeCrack_two(double X, double Y, double kappa, double * AA, double * BB, double * AD, double * BD, int max_iter = 100);
};

/////////////////////////////////////////////////
//...parametrization of the sample (preliminary);
#define DRAFT_N                     2    //...degree of multipoles;
#undef  DRAFT_N                          //...param(0);
#define DRAFT_Q                     0.92 //...normalization coefficient;
#undef  DRAFT_Q                          //...param(1);
#define DRAFT_N_quad                0.   //...parameters of quadrature;
#undef  DRAFT_N_quad                     //...param(2);
#define DRAFT_A							0.   //...normal surface adhesion parameter;
#undef  DRAFT_A                          //...param(3);
#define DRAFT_B                     0.   //...shear surface adhesion parameter;
#undef  DRAFT_B                          //...param(4);
#define DRAFT_C1                    0.   //...cohesion parameter in matrix;
#undef  DRAFT_C1                         //...param(5);
#define DRAFT_G1                    1.   //...shear modulus in matrix;
#undef  DRAFT_G1                         //...param(6);
#define DRAFT_nju1                  0.3  //...Poisson coefficient in matrix;
#undef  DRAFT_nju1                       //...param(7);
#define DRAFT_kappa_kk					0.   //...compression scale factor in matrix;
#undef  DRAFT_kappa_kk						  //...param(8);
#define DRAFT_kappa_mu              0.   //...shear scale factor in matrix;
#undef  DRAFT_kappa_mu                   //...param(9);
#define DRAFT_C2                    0.   //...cohesion parameter in inclusion;
#undef  DRAFT_C2                         //...param(10);
#define DRAFT_G2                    1.   //...shear modulus in inclusion;
#undef  DRAFT_G2                         //...param(11);
#define DRAFT_nju2                  0.3  //...Poisson coefficient in inclusion;
#undef  DRAFT_nju2                       //...param(12);
#define DRAFT_kappa_kk2					0.   //...compression scale factor in inclusion;
#undef  DRAFT_kappa_kk2						  //...param(13);
#define DRAFT_kappa_mu2             0.   //...shear scale factor in inclusion;
#undef  DRAFT_kappa_mu2                  //...param(14);
#define DRAFT_C3                    0.   //...cohesion parameter in layer;
#undef  DRAFT_C3                         //...param(15);
#define DRAFT_G3                    1.   //...shear modulus in layer;
#undef  DRAFT_G3                         //...param(16);
#define DRAFT_nju3                  0.3  //...Poisson coefficient in layer;
#undef  DRAFT_nju3                       //...param(17);
#define DRAFT_kappa_kk3					0.   //...compression scale factor in layer;
#undef  DRAFT_kappa_kk3						  //...param(18);
#define DRAFT_kappa_mu3             0.   //...shear scale factor in layer;
#undef  DRAFT_kappa_mu3                  //...param(19);
#define DRAFT_C4                    0.   //...cohesion parameter in layer;
#undef  DRAFT_C4                         //...param(20);
#define DRAFT_G4                    1.   //...shear modulus in layer;
#undef  DRAFT_G4                         //...param(21);
#define DRAFT_nju4                  0.3  //...Poisson coefficient in layer;
#undef  DRAFT_nju4                       //...param(22);
#define DRAFT_kappa_kk4					0.   //...compression scale factor in layer;
#undef  DRAFT_kappa_kk4						  //...param(23);
#define DRAFT_kappa_mu4             0.   //...shear scale factor in layer;
#undef  DRAFT_kappa_mu4                  //...param(24);
#define DRAFT_lagrange              1.   //...Lagrange coefficient for energy in block functional;
#undef  DRAFT_lagrange                   //...param(25);
#endif
