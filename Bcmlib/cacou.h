/*=========================================*/
/*                 CACOU                   */
/*=========================================*/
#ifndef ___CACOU___
#define ___CACOU___

#include "cshapes.h"

//////////////////////////////////////////////
//														  //
//  BASIC SYSTEM of ELLIPSOIDAL MULTIPOLES  //
//														  //
//////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//...class of acoustic multipoles, based on plane wave in Z-direction;
class CAcou3DEll : public CShape<double> {
public:
		int freedom (int m) { return (m+1)*(m+2);}
		int size_of_param() { return(5);};
public:
//...initialization and calculation of multipoles;
		void parametrization     (double * P = NULL, int m_dop = 0);
		void parametrization_grad(double * P = NULL, int m_dop = 0);
		void parametrization_hess(double * P = NULL, int m_dop = 0);
//...differentiation;
		void deriv_X(double * deriv, double f = 1.);
		void deriv_Y(double * deriv, double f = 1.);
		void deriv_Z(double * deriv, double f = 1.);
public:
		void init() {
			delete_struct(param);
			param = new_struct<Param>(size_of_param());
		}
public:
//...constructor;
		CAcou3DEll() {
			init();
		};
};

///////////////////////////////////////
//...parametrization of the multipoles;
#define SHAPE_norm_kappa                0. //...normalized wave number;
#undef  SHAPE_norm_kappa                   //...param(0);
#define SHAPE_norm_kappa_inv            1. //...inverse of normalized wave number;
#undef  SHAPE_norm_kappa_inv               //...param(1);
#define SHAPE_norm_kappa_sqr            0. //...square of normalized wave number;
#undef  SHAPE_norm_kappa_sqr               //...param(2);
#define SHAPE_kappa                     0. //...non-normalized wave number;
#undef  SHAPE_kappa                        //...param(3);
#define SHAPE_kappa_shift               0. //...shift-normalized wave number;
#undef  SHAPE_kappa_shift                  //...param(4);


///////////////////////////////////////////////////////////
//																			//
//  ACOUSTIC MULTIPOLES with POLYNOMIAL CHARACTERISTIC   //
//																			//
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
//...class of acoustic multipoles with polynomial behaviour along the line;
class CAcou3DPoly : public CShape<complex> {
protected:
		double * au, * E1;
public:
		int freedom (int m) { return sqr(m+1);}
		int size_of_param() { return(4);}
public:
//...calculation of multipoles;
		void parametrization(double * P = NULL, int m_dop = 0);
//...differentiation;
		void deriv_X(complex * deriv, double f = 1.);
		void deriv_Y(complex * deriv, double f = 1.);
		void deriv_Z(complex * deriv, double f = 1.);
public:
		void init() {
			delete_struct(param);
			param = new_struct<Param>(size_of_param());
		}
public:
//...constructor and destructor;
		CAcou3DPoly() {
			au = E1 = NULL;
			init();
		};
		virtual ~CAcou3DPoly(void) { 
			delete_struct(au);
			delete_struct(E1);
		}
};

///////////////////////////////////////
//...parametrization of the multipoles;
#define SHAPE_norm_kappa                0. //...normalized wave number;
#undef  SHAPE_norm_kappa                   //...param(0);
#define SHAPE_norm_kappa_inv            1. //...inverse of normalized wave number;
#undef  SHAPE_norm_kappa_inv               //...param(1);
#define SHAPE_norm_kappa_sqr            0. //...square of normalized wave number;
#undef  SHAPE_norm_kappa_sqr               //...param(2);
#define SHAPE_kappa                     0. //...non-normalized wave number;
#undef  SHAPE_kappa                        //...param(3);

////////////////////////////////////////////////////////////
//																			 //
//  ACOUSTIC MULTIPOLES for SPHERE INTERIOUR OR EXTERIOR  //
//																			 //
////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////
//...class of classic acoustic multipoles for interior or exterior;
class CAcou3DZoom : public CShape<double> {
protected:
		double * au;
public:
		int freedom (int m) { return sqr(m+1);}
		int size_of_param() { return(4);}
public:
//...initialization and calculation of multipoles;
		void parametrization     (double * P = NULL, int m_dop = 0);
		void parametrization_grad(double * P = NULL, int m_dop = 0);
		void parametrization_hess(double * P = NULL, int m_dop = 0);
//...differentiation;
		void deriv_X(double * deriv, double f = 1.);
		void deriv_Y(double * deriv, double f = 1.);
		void deriv_Z(double * deriv, double f = 1.);
public:
		void init() {
			delete_struct(param);
			param = new_struct<Param>(size_of_param());
			//id_cmpl = 1; //...we can install id_cmpl by void change()!!!
		}
public:
//...constructor and destructor;
		CAcou3DZoom(int inverse = 0) {
			au = NULL;
			id_inverse = inverse;
			init();
		};
		virtual ~CAcou3DZoom(void) { 
			delete_struct(au);
		}
};

///////////////////////////////////////
//...parametrization of the multipoles;
#define SHAPE_norm_kappa                0. //...normalized wave number;
#undef  SHAPE_norm_kappa                   //...param(0);
#define SHAPE_norm_kappa_inv            1. //...inverse of normalized wave number;
#undef  SHAPE_norm_kappa_inv               //...param(1);
#define SHAPE_norm_kappa_sqr            0. //...square of normalized wave number;
#undef  SHAPE_norm_kappa_sqr               //...param(2);
#define SHAPE_kappa                     0. //...non-normalized wave number;
#undef  SHAPE_kappa                        //...param(3);


//////////////////////////////////////////
//													 //
//     SUPERPOSITION of PLANE WAVES     //
//													 //
//////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
//...class of exponential plane waves (separated as complex potential with real and imaginary part);
class CAcou3DWave : public CShape<double> {
protected:
       double * cz, * sz, * cy, * sy;
public:
		 int freedom (int m) { return sqr(m+1);}
       int size_of_param() { return(4);}
public:
//...initialization and calculation of multipoles;
       void degree_init1(int N, int dim) { degree_init2(2, N, dim);};
       void degree_init2(int id_flag, int N, int dim);
       void parametrization(double * P, int m_dop = 1);
       int  add_expo_power (double X, double Y, double Z, int flag = 0);
//...differentiation;
       void deriv_X(double * deriv, double f = 1.);
       void deriv_Y(double * deriv, double f = 1.);
       void deriv_Z(double * deriv, double f = 1.);
       void deriv_N();
public:
		void init() {
			delete_struct(param);
			param = new_struct<Param>(size_of_param());
			id_cmpl = 1;
		}
public:
//...constructor and destructor;
		CAcou3DWave() {
			cz = sz = cy = sy = NULL;
			init();
		};
		virtual ~CAcou3DWave(void) { 
			 delete_struct(cz);
			 delete_struct(sz);
			 delete_struct(cy);
			 delete_struct(sy);
		}
};

///////////////////////////////////////
//...parametrization of the multipoles;
#define SHAPE_norm_kappa                0. //...normalized wave number;
#undef  SHAPE_norm_kappa                   //...param(0);
#define SHAPE_norm_kappa_inv            1. //...inverse of normalized wave number;
#undef  SHAPE_norm_kappa_inv               //...param(1);
#define SHAPE_norm_kappa_sqr            0. //...square of normalized wave number;
#undef  SHAPE_norm_kappa_sqr               //...param(2);
#define SHAPE_kappa                     0. //...non-normalized wave number;
#undef  SHAPE_kappa                        //...param(3);


////////////////////////////////////////////////////////////
//																			 //
//  BASIC SYSTEM of ACOUSTIC MULTIPOLES (acoustic beam)   //
//																			 //
////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
//...class of acoustic multipoles of beams type with two oscillation characteristics;
class CAcou3DBeam : public CShape<double> {
protected:
		double  * au, * az, * pim, * pxim, * pyim, * pzim;
		complex * E1;
public:
		Num_Shape	 type() { return AU3D_BEAM_SHAPE;}
		int freedom (int m) { return sqr(m+1);}
		int size_of_param() { return(6);}
public:
//...initialization of multipoles;
		void set_shape(double R0, double kk = 0., double kk_dop = 0., double L1 = 0., double L2 = 0.);
public:
//...calculation of multipoles;
		void parametrization     (double * P = NULL, int m_dop = 0);
		void parametrization_grad(double * P = NULL, int m_dop = 0);
		void parametrization_hess(double * P = NULL, int m_dop = 0);
//...differentiation;
		void deriv_X(double * deriv, double f = 1.);
		void deriv_Y(double * deriv, double f = 1.);
		void deriv_Z(double * deriv, double f = 1.);
public:
		void init() {
			delete_struct(param);
			param = new_struct<Param>(size_of_param());
		}
public:
//...constructor and destructor;
		CAcou3DBeam() {
			au = az = pim = pxim = pyim = pzim = NULL;
			E1 = NULL;
			init();
		}
		virtual ~CAcou3DBeam(void) { 
			delete_struct(pim);
			delete_struct(pxim);	delete_struct(pyim);	delete_struct(pzim);
			delete_struct(au);	delete_struct(az);
			delete_struct(E1);
		}
};

///////////////////////////////////////
//...parametrization of the multipoles;
#define SHAPE_norm_kappa                0. //...normalized cylindrical wave number;
#undef  SHAPE_norm_kappa                   //...param(0);
#define SHAPE_norm_kappa_inv            1. //...inverse of normalized wave number;
#undef  SHAPE_norm_kappa_inv               //...param(1);
#define SHAPE_norm_kappa_sqr            0. //...square of normalized wave number;
#undef  SHAPE_norm_kappa_sqr               //...param(2);
#define SHAPE_kappa                     0. //...non-normalized wave number;
#undef  SHAPE_kappa                        //...param(3);
#define SHAPE_norm_kappa_dop            0. //...normalized directional wave number;
#undef  SHAPE_norm_kappa_dop               //...param(4);
#define SHAPE_kappa_dop                 0. //...non-normalized wave number;
#undef  SHAPE_kappa_dop                    //...param(5);

/////////////////////////////////////////////////////////////////////////////////////////
//...class of acoustic multipoles of beams type with oscillation characteristics along Z;
class CAcou3DBeamZ : public CShape<double> {
protected:
		double  * az, * pim, * pxim, * pyim, * pzim;
		complex * E1;
public:
		int freedom (int m) { return sqr(m+1);}
		int size_of_param() { return(4);};
public:
//...initialization and calculation of multipoles;
		void parametrization     (double * P = NULL, int m_dop = 0);
		void parametrization_grad(double * P = NULL, int m_dop = 0);
		void parametrization_hess(double * P = NULL, int m_dop = 0);
//...differentiation;
		void deriv_X(double * deriv, double f = 1.);
		void deriv_Y(double * deriv, double f = 1.);
		void deriv_Z(double * deriv, double f = 1.);
public:
		void init() {
			delete_struct(param);
			param = new_struct<Param>(size_of_param());
		}
public:
//...constructor and destructor;
		CAcou3DBeamZ() {
			az = pim = pxim = pyim = pzim = NULL;
			E1 = NULL;
			init();
		};
		virtual ~CAcou3DBeamZ(void) { 
			delete_struct(pim);
			delete_struct(pxim);	delete_struct(pyim);	delete_struct(pzim);
			delete_struct(az);
			delete_struct(E1);
		}
};

///////////////////////////////////////
//...parametrization of the multipoles;
#define SHAPE_norm_kappa                0. //...normalized wave number (directional);
#undef  SHAPE_norm_kappa                   //...param(0);
#define SHAPE_norm_kappa_inv            1. //...inverse of normalized wave number;
#undef  SHAPE_norm_kappa_inv               //...param(1);
#define SHAPE_norm_kappa_sqr            0. //...square of normalized wave number;
#undef  SHAPE_norm_kappa_sqr               //...param(2);
#define SHAPE_kappa                     0. //...non-normalized wave number;
#undef  SHAPE_kappa                        //...param(3);
#endif
