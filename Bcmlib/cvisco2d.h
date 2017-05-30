/*==========================================*/
/*                  CVISCO2D                 */
/*==========================================*/
#ifndef ___CVisco2D___
#define ___CVisco2D___

#include "ccomput2d.h"

/////////////////////////////////////////////////
//...class of blocks partition for plane problem;
class CVisco2D : public CComput2D<complex> {
public:
		int size_of_param() { return(19);}
public:
		void init() {
			delete_struct(param);
			param = new_struct<Param>(size_of_param());
		}
public:
//...constructor;
		CVisco2D (int num_phase = 8) {
			NUM_PHASE = num_phase;
			init();
		};
protected:
static int NUM_BASIC, NUM_SHIFT, NUM_SHEAR;//...временная величина;
		int  block_shape_init(Block<complex> & B, Num_State id_free);
//...auxilliary operations with block matrix;
		void jump1_classic_x (double * P, int i, int m);
		void jump1_classic_y (double * P, int i, int m);
		void jump2_classic_x (double * P, int i, int m);
		void jump2_classic_y (double * P, int i, int m);
		void jump3_classic_x (double * P, int i, int m);
		void jump3_classic_y (double * P, int i, int m);
		void jump4_classic_x (double * P, int i, int m);
		void jump4_classic_y (double * P, int i, int m);
		void jump5_classic_x (double * P, int i, int m);
		void jump5_classic_y (double * P, int i, int m);
		void jump_make_local (int i, int m);
		void jump_make_common(int i, int m);
//...forming block matrix elements;
		Num_State  gram1    (CGrid * nd, int i, int id_local);
		Num_State  gram2    (CGrid * nd, int i, int id_local);
		Num_State  gram3    (CGrid * nd, int i, int id_local);
		Num_State  gram4    (CGrid * nd, int i, int id_local);
		Num_State  transfer1(CGrid * nd, int i, int k, int id_local);
		Num_State  transfer2(CGrid * nd, int i, int k, int id_local);
		Num_State  transfer3(CGrid * nd, int i, int k, int id_local);
		Num_State  transfer4(CGrid * nd, int i, int k, int id_local){ return transfer3(nd, i, k, id_local);};
		Num_State  rigidy1  (CGrid * nd, int i, complex * K);
		Num_State  counting_header (int Num);
public:
//...operations with block matrix;
		void set_fasa_hmg(double GM_re, double GM_im, double KM_re, double KM_im, double GI_re, double GI_im, double KI_re, double KI_im) {
			  set_fasa_hmg(GM_re, GM_im, KM_re, KM_im, GI_re, GI_im, KI_re, KI_im, 0., 0., 1., 0.);
		}
		void set_fasa_hmg(double GM_re, double GM_im, double KM_re, double KM_im, double GI_re, double GI_im, double KI_re, double KI_im, double GL_re, double GL_im, double KL_re, double KL_im);
		void GetFuncAllValues(double X, double Y, double Z, complex * F, int id_block, Num_Value id_F, int id_variant = 0, int iparam = 0);
		void GetRigidy(complex * K, Num_Comput Num = BASIC_COMPUT);
		void GetEnergy(complex * energy, Num_Value _FMF = ENERGY_VALUE);
		void GetEnergyValue(int k, double * energy);
//...аналитические модели с межфазным слоем и адгезией;
		complex TakeSphereEshelby(int nn, double * RR, complex * KK, complex * MU);
		complex TakeSphereEshelby_grad(double R, double c0, double * par);
//...аналитические модели (самосогласованные);
		complex TakeEshelby_visco_volm_two (double ff);
		complex TakeEshelby_visco_shear_two(double ff, double eps = EE, int max_iter = 100);
};

/////////////////////////////////////////////////
//...parametrization of the sample (preliminary);
#define DRAFT_N                     2    //...degree of multipoles;
#undef  DRAFT_N                          //...param(0);
#define DRAFT_Q                     0.92 //...normalization coefficient;
#undef  DRAFT_Q                          //...param(1);
#define DRAFT_N_quad                0.   //...parameters of quadrature;
#undef  DRAFT_N_quad                     //...param(2);
#define DRAFT_requl						1.   //...normalization coefficient in reqularization;
#undef  DRAFT_requl							  //...param(3);
#define DRAFT_R1							1.	  //...geometry of inclusion;
#undef  DRAFT_R1								  //...param(4);
#define DRAFT_R2							1.	  //...geometry of intermediate layer;
#undef  DRAFT_R2								  //...param(5);
#define DRAFT_GM_re                 1.   //...shear modulus in matrix;
#undef  DRAFT_GM_re                      //...param(6);
#define DRAFT_GM_im                 1.   //...shear modulus in matrix;
#undef  DRAFT_GM_im                      //...param(7);
#define DRAFT_KM_re                 1.   //...volume modulus in matrix;
#undef  DRAFT_KM_re                      //...param(8);
#define DRAFT_KM_im                 1.   //...volume modulus in matrix;
#undef  DRAFT_KM_im                      //...param(9);
#define DRAFT_GI_re                 1.   //...shear modulus in inclusion;
#undef  DRAFT_GI_re                      //...param(10);
#define DRAFT_GI_im                 1.   //...shear modulus in inclusion;
#undef  DRAFT_GI_im                      //...param(11);
#define DRAFT_KI_re                 1.   //...volume modulus in inclusion;
#undef  DRAFT_KI_re                      //...param(12);
#define DRAFT_KI_im                 1.   //...volume modulus in inclusion;
#undef  DRAFT_KI_im                      //...param(13);
#define DRAFT_GL_re                 1.   //...shear modulus in interphase layer;
#undef  DRAFT_GL_re                      //...param(14);
#define DRAFT_GL_im                 1.   //...shear modulus in interphase layer;
#undef  DRAFT_GL_im                      //...param(15);
#define DRAFT_KL_re                 1.   //...volume modulus in interphase layer;
#undef  DRAFT_KL_re                      //...param(16);
#define DRAFT_KL_im                 1.   //...volume modulus in interphase layer;
#undef  DRAFT_KL_im                      //...param(17);
#define DRAFT_lagrange              1.   //...Lagrange coefficient for energy in block functional;
#undef  DRAFT_lagrange                   //...param(18);
#endif
