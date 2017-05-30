/*========================================*/
/*                 CMAPI                  */
/*========================================*/
#ifndef ___CMAPI___
#define ___CMAPI___

#include "cshapes.h"

//////////////////////////////////////////////
//														  //
//      SYSTEMS of SPHEROID MULTIPOLES      //
//														  //
//////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//...class of acoustic multipoles, based on three plane waves;
class CMapi3DSpheroid : public CShape<double>  {
public:
		int freedom (int m) { return 4;} //...вычисляем только четыре функции; 
		int size_of_param() { return(2);}
public:
//...initialization of multipoles;
		void set_shape(double R0, double eps = 1., double p1 = 1., double p2 = 0., double p3 = 0., double p4 = 0.);
public:
//...calculation of multipoles;
		void parametrization		 (double * P = NULL, int m = 0);
		void parametrization_grad(double * P = NULL, int m = 0) { parametrization(P, m);} 
		void parametrization_hess(double * P = NULL, int m = 0) { parametrization(P, m);} 
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
		CMapi3DSpheroid() {
			init();
		}
};

class CMapi3DSpheroidFull : public CShape<double>  {
public:
		int freedom(int m)  { return sqr(m+1);}
		int size_of_param() { return(2);}
public:
//...initialization of multipoles;
		void set_shape(double R0, double eps = 1., double p1 = 1., double p2 = 0., double p3 = 0., double p4 = 0.);
public:
//...calculation of multipoles;
		void parametrization		 (double * P = NULL, int m_dop = 0);
		void parametrization_grad(double * P = NULL, int m_dop = 0) { parametrization(P, m_dop);} 
		void parametrization_hess(double * P = NULL, int m_dop = 0) { parametrization(P, m_dop);} 
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
		CMapi3DSpheroidFull() {
			init();
	}
};

///////////////////////////////////////
//...parametrization of the multipoles;
#define SHAPE_spheroid_eps					0. //...параметр сфероида -- соотношение полуосей;
#undef  SHAPE_spheroid_eps						//...param(0);
#define SHAPE_spheroid_A					1. //...большая полуось сфероида;
#undef  SHAPE_spheroid_A						//...param(1);
#endif
