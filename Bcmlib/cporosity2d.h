/*==============================================*/
/*                  CPOROSITY2D                 */
/*==============================================*/
#ifndef ___CPorocity2D___
#define ___CPorocity2D___

#include "ccomput2d.h"

///////////////////////////////////////////////////////////////////
//...class of blocks partition for plane problem (one phase model);
class CPorosity2D : public CComput2D<double> {
public:
		int size_of_param() { return(13);}
public:
		void init() {
			delete_struct(param);
			param = new_struct<Param>(size_of_param());
		}
public:
//...constructor;
		CPorosity2D (int num_phase = 7) {
			NUM_PHASE = num_phase;
			MAX_PHASE = 1;
			init();
		};
protected:
static int NUM_SHEAR;
		int  block_shape_init(Block<double> & B, Num_State id_free);
//...auxilliary operations with block matrix;
		void jump1_common_x(double * P, int i, int m);
		void jump1_common_y(double * P, int i, int m);
		void jump1_common_h(double * P, int i, int m);
		void jump2_common_x(double * P, int i, int m);
		void jump2_common_y(double * P, int i, int m);
		void jump2_common_h(double * P, int i, int m);
		void jump4_common_x(double * P, int i, int m);
		void jump4_common_y(double * P, int i, int m);
		void jump4_common_h(double * P, int i, int m);
		void jump_make_local (int i, int m);
		void jump_make_common(int i, int m);
//...forming block matrix elements;
		Num_State gram1    (CGrid * nd, int i, int id_local);
		Num_State gram4    (CGrid * nd, int i, int id_local);
		Num_State transfer1(CGrid * nd, int i, int k, int id_local);
		Num_State transfer4(CGrid * nd, int i, int k, int id_local);
		Num_State computing_header(Num_Comput Num);
public:
//...problem parameters;
		void set_fasa_hmg(double nju, double G1, double C0, double alpha, double gamma);
//...problem solving results;
		void GetFuncAllValues(double X, double Y, double Z, double * F, int id_block, Num_Value id_F, int id_variant = 0, int iparam = 0);
};

/////////////////////////////////////////////////
//...parametrization of the sample (preliminary);
#define DRAFT_N                     2    //...degree of multipoles;
#undef  DRAFT_N                          //...param(0);
#define DRAFT_Q                     0.92 //...normalization coefficient;
#undef  DRAFT_Q                          //...param(1);
#define DRAFT_N_quad                0.   //...parameters of quadrature;
#undef  DRAFT_N_quad                     //...param(2);
#define DRAFT_local						0.   //...using gradient porosity;
#undef  DRAFT_local                      //...param(3);
#define DRAFT_alpha                 0.   //...volume porosity parameter;
#undef  DRAFT_alpha                      //...param(4);
#define DRAFT_gamma                 0.   //...gradient porosity parameter;
#undef  DRAFT_gamma                      //...param(5);
#define DRAFT_C0							1.   //...parameter of gradient defectness;
#undef  DRAFT_C0								  //...param(6);
#define DRAFT_G1                    1.   //...shear modulus;
#undef  DRAFT_G1                         //...param(7);
#define DRAFT_nju                   0.3  //...Poisson coefficient;
#undef  DRAFT_nju                        //...param(8);
#define DRAFT_kappa                 0    //...parameter of porosity remission;
#undef  DRAFT_kappa                      //...param(9);
#define DRAFT_kappa_beta            0    //...additional parameter of porosity remission;
#undef  DRAFT_kappa_beta                 //...param(10);
#define DRAFT_kappa_PN              1    //...parameter of porosity remission in representation;
#undef  DRAFT_kappa_PN                   //...param(11);
#define DRAFT_lagrange              1.   //...Lagrange coefficient for energy in block functional;
#undef  DRAFT_lagrange                   //...param(12);
#endif
