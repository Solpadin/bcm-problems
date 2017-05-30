/*=========================================*/
/*                 CRADII                  */
/*=========================================*/
#ifndef ___CRADII___
#define ___CRADII___

#include "cshapes.h"

/////////////////////////////////////////////////////////
//          													    //
//   RADII MULTIPOLES WITH POLYNOMIAL CHARACTERISTIC   //
//          													    //
/////////////////////////////////////////////////////////
///////////////////////////////////////////////////////
//...class of skin multipoles with polynomial behaviour;
template <typename T>
class CSkin2DRadiiPoly : public CShape<T> {
protected:
      T * hh;
public:
		int freedom (int m) { return m*2+1;}
		int size_of_param() { return(4);}
public:
//...initialization and calculation of multipoles;
		void parametrization     (double * P = NULL, int m_dop = 0);
		void parametrization_grad(double * P = NULL, int m_dop = 0);
		void parametrization_hess(double * P = NULL, int m_dop = 0);
//...differentiation;
		void deriv_X(T * deriv, double f = 1.);
		void deriv_Y(T * deriv, double f = 1.);
public:
		void init() {
			delete_struct(this->param);
			this->param = new_struct<Param>(size_of_param());
		}
public:
//...constructor and destructor;
		CSkin2DRadiiPoly() {
			hh = NULL;
			init();
		}
		virtual ~CSkin2DRadiiPoly(void) { 
			delete_struct(hh);
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

///////////////////////////////////////
//...parametrization of the multipoles;
template <typename T>
void CSkin2DRadiiPoly<T>::parametrization(double * P, int m_dop)
{
	if (! this->p) {
		this->p = new_struct<T>(this->NN_dop = freedom(this->N+m_dop));
///////////////////////
//...technical support;
		delete_struct(hh); hh = new_struct<T>(max(this->N+1+m_dop, 2));
	}
	if (P && this->p) {
		complex w = comp(P[0], P[1])*this->R0_inv;
		T		 zm = 1., zm_i = 0., zz;
		double rr = abs(w);
		this->cs[0] = P[3];
		this->cs[1] = P[4];

///////////////////////////////////
//...calculation of the multipoles;
		if (hh) {
			T dd = sqr(.5*this->param[0]*rr), f = 1.;
			for (int k, m = 0; m <= this->N+m_dop; m++) { 
				k = 0; hh[m] = (f = 1.);
				do {
					k += 1; hh[m] += (f *= dd/(k*(k+m)));
				} 
				while (fabs(f) > EE);
			}
		}
		int m;
		for ( this->p[0] = hh[0], m = 1; m <= this->N+m_dop; m++) {
				this->p[m*2-1] = hh[m]*(zm_i = zm*imag(w)+(zz = zm_i)*real(w));
				this->p[m*2]   = hh[m]*(zm = zm*real(w)-zz*imag(w));
		}
	}
}

////////////////////////////////////////////
//...parametrization grad of the multipoles;
template <typename T>
void CSkin2DRadiiPoly<T>::parametrization_grad(double * P, int m_dop)
{
	parametrization(P, m_dop+1);
	if (! P) {
		delete_struct(this->px);
		delete_struct(this->py);
	}
	if (P && ! this->px && ! this->py) {
		this->px = new_struct<T>(this->NN_grad = freedom(this->N+m_dop));
		this->py = new_struct<T>(this->NN_grad);
	}
	if (P && this->px && this->py) {
		memset (this->px, 0, this->NN_grad*sizeof(T));
		memset (this->py, 0, this->NN_grad*sizeof(T));
		this->N += m_dop;
		deriv_X(this->px);
		deriv_Y(this->py);
		this->N -= m_dop;
	}
}

//////////////////////////////////////////////
//...parametrization hessian of the multipoles;
template <typename T>
void CSkin2DRadiiPoly<T>::parametrization_hess(double * P, int m_dop)
{
	parametrization_grad(P, m_dop+1);
	if (! P) {
		delete_struct(this->pxx);
		delete_struct(this->pxy);
		delete_struct(this->pyy);
	}
	if (P && ! this->pxx && ! this->pxy && ! this->pyy) {
		this->pxx = new_struct<T>(this->NN_hess = freedom(this->N+m_dop));
		this->pxy = new_struct<T>(this->NN_hess);
		this->pyy = new_struct<T>(this->NN_hess);
	}
	if (P && this->pxx && this->pxy && this->pyy) {
		memset (this->pxx, 0, this->NN_hess*sizeof(T));
		memset (this->pxy, 0, this->NN_hess*sizeof(T));
		memset (this->pyy, 0, this->NN_hess*sizeof(T));
		this->N += m_dop;
		swap(this->p, this->px);
			deriv_X(this->pxx);
		swap(this->p, this->px);
		swap(this->p, this->py);
			deriv_X(this->pxy);
			deriv_Y(this->pyy);
		swap(this->p, this->py);
		this->N -= m_dop;
	}
}

///////////////////////////////////////
//...differentiation of the multipoles;
template <typename T>
void CSkin2DRadiiPoly<T>::deriv_X(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;
	int m;
	for (m = this->N; m > 1; m--) {
			deriv[m*2-1] += (m*this->p[m*2-3]+.25*this->param[2]*this->p[m*2+1]/(m+1.))*f;
			deriv[m*2]   += (m*this->p[m*2-2]+.25*this->param[2]*this->p[m*2+2]/(m+1.))*f;
	}
	if (m > 0) {
			deriv[m*2-1] += .25*this->param[2]*this->p[m*2+1]/(m+1.)*f;
			deriv[m*2]   += (m*this->p[m*2-2]+.25*this->param[2]*this->p[m*2+2]/(m+1.))*f;
	}
	deriv[0] += .5*this->param[2]*this->p[2]*f;
}

///////////////////////////////////////
//...differentiation of the multipoles;
template <typename T>
void CSkin2DRadiiPoly<T>::deriv_Y(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;
	int m;
	for (m = this->N; m > 1; m--) {
			deriv[m*2-1] += (m*this->p[m*2-2]-.25*this->param[2]*this->p[m*2+2]/(m+1.))*f;
			deriv[m*2]   -= (m*this->p[m*2-3]-.25*this->param[2]*this->p[m*2+1]/(m+1.))*f;
	}
	if (m > 0) {
			deriv[m*2-1] += (m*this->p[m*2-2]-.25*this->param[2]*this->p[m*2+2]/(m+1.))*f;
			deriv[m*2]   += .25*this->param[2]*this->p[m*2+1]/(m+1.)*f;
	}
	deriv[0] += .5*this->param[2]*this->p[1]*f;
}

/////////////////////////////////////////////////////////
//...class of skin multipoles with logarothmic behaviour;
template <typename T>
class CSkin2DRadiiZoom : public CShape<T> {
protected:
      T * hh;
public:
		int size_of_param() { return(4);}
		int freedom (int m) { return m*2+1;}
public:
//...initialization and calculation of multipoles;
		void parametrization     (double * P = NULL, int m_dop = 0);
		void parametrization_grad(double * P = NULL, int m_dop = 0);
		void parametrization_hess(double * P = NULL, int m_dop = 0);
//...differentiation;
		void deriv_X(T * deriv, double f = 1.);
		void deriv_Y(T * deriv, double f = 1.);
public:
		void init() {
			delete_struct(this->param);
			this->param = new_struct<Param>(size_of_param());
		}
public:
//...constructor and destructor;
		CSkin2DRadiiZoom() {
			hh = NULL;
			init();
		}
		virtual ~CSkin2DRadiiZoom(void) { 
			delete_struct(hh);
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

///////////////////////////////////////
//...parametrization of the multipoles;
template <typename T>
void CSkin2DRadiiZoom<T>::parametrization(double * P, int m_dop)
{
	if (! this->p) {
		this->p = new_struct<T>(this->NN_dop = freedom(this->N+m_dop));
///////////////////////
//...technical support;
		delete_struct(hh); hh = new_struct<T>(this->N+1+m_dop);
	}
	if (P && this->p) {
		complex w = comp(P[0], P[1])*this->R0_inv;
		T		 zm = 1., zm_i = 0., zz;
		double rr = abs(w);
		this->cs[0] = P[3];
		this->cs[1] = P[4];

///////////////////////////////////
//...calculation of the multipoles;
		if (hh) {
			T dd = sqr(.5*this->param[0]*rr), g_euler = T(0.577215664901532860606512), dd_inv = 1./dd, h_m = 0., h_k, f;
			for ( int k, m = 0; m <= this->N+m_dop; m++, h_m +=  T(1.)/m) { 
				for (hh[m]  = (f = -.5*dd_inv*m), k = 1; k < m; k++) 
					  hh[m] += (f *= -k*(m-k)*dd_inv);

				k = 0;  hh[m] += (h_k = -g_euler+(h_m-log(to_double(dd)))*.5)*(f = 1.);
				do {
					k += 1; hh[m] += (h_k += (k+m*.5)/(k*(k+m)))*(f *= dd/(k*(k+m)));
				} 
				while (to_double(f*(1.+fabs(h_k))) > EE);
			}
		}
		int m;
		for ( this->p[0] = hh[0], m = 1; m <= this->N+m_dop; m++) {
				this->p[m*2-1] = hh[m]*(zm_i = zm*imag(w)+(zz = zm_i)*real(w));
				this->p[m*2]   = hh[m]*(zm = zm*real(w)-zz*imag(w));
		}
	}
}

////////////////////////////////////////////
//...parametrization grad of the multipoles;
template <typename T>
void CSkin2DRadiiZoom<T>::parametrization_grad(double * P, int m_dop)
{
	parametrization(P, m_dop+1);
	if (! P) {
		delete_struct(this->px);
		delete_struct(this->py);
	}
	if (P && ! this->px && ! this->py) {
		this->px = new_struct<T>(this->NN_grad = freedom(this->N+m_dop));
		this->py = new_struct<T>(this->NN_grad);
	}
	if (P && this->px && this->py) {
		memset (this->px, 0, this->NN_grad*sizeof(T));
		memset (this->py, 0, this->NN_grad*sizeof(T));
		this->N += m_dop;
		deriv_X(this->px);
		deriv_Y(this->py);
		this->N -= m_dop;
	}
}

//////////////////////////////////////////////
//...parametrization hessian of the multipoles;
template <typename T>
void CSkin2DRadiiZoom<T>::parametrization_hess(double * P, int m_dop)
{
	parametrization_grad(P, m_dop+1);
	if (! P) {
		delete_struct(this->pxx);
		delete_struct(this->pxy);
		delete_struct(this->pyy);
	}
	if (P && ! this->pxx && ! this->pxy && ! this->pyy) {
		this->pxx = new_struct<T>(this->NN_hess = freedom(this->N+m_dop));
		this->pxy = new_struct<T>(this->NN_hess);
		this->pyy = new_struct<T>(this->NN_hess);
	}
	if (P && this->pxx && this->pxy && this->pyy) {
		memset (this->pxx, 0, this->NN_hess*sizeof(T));
		memset (this->pxy, 0, this->NN_hess*sizeof(T));
		memset (this->pyy, 0, this->NN_hess*sizeof(T));
		this->N += m_dop;
		swap(this->p, this->px);
			deriv_X(this->pxx);
		swap(this->p, this->px);
		swap(this->p, this->py);
			deriv_X(this->pxy);
			deriv_Y(this->pyy);
		swap(this->p, this->py);
		this->N -= m_dop;
	}
}

///////////////////////////////////////
//...differentiation of the multipoles;
template <typename T>
void CSkin2DRadiiZoom<T>::deriv_X(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;
	int m;
	for (m = this->N; m > 1; m--) {
			deriv[m*2-1] += (m*this->p[m*2-3]+.25*this->param[2]*this->p[m*2+1]/(m+1.))*f;
			deriv[m*2]   += (m*this->p[m*2-2]+.25*this->param[2]*this->p[m*2+2]/(m+1.))*f;
	}
	if (m > 0) {
			deriv[m*2-1] += .25*this->param[2]*this->p[m*2+1]/(m+1.)*f;
			deriv[m*2]   += (m*this->p[m*2-2]+.25*this->param[2]*this->p[m*2+2]/(m+1.))*f;
	}
	deriv[0] += .5*this->param[2]*this->p[2]*f;
}

///////////////////////////////////////
//...differentiation of the multipoles;
template <typename T>
void CSkin2DRadiiZoom<T>::deriv_Y(T * deriv, double f)
{
	if (! deriv || ! this->p) return;
	f *= this->R0_inv;
	int m;
	for (m = this->N; m > 1; m--) {
			deriv[m*2-1] += (m*this->p[m*2-2]-.25*this->param[2]*this->p[m*2+2]/(m+1.))*f;
			deriv[m*2]   -= (m*this->p[m*2-3]-.25*this->param[2]*this->p[m*2+1]/(m+1.))*f;
	}
	if (m > 0) {
			deriv[m*2-1] += (m*this->p[m*2-2]-.25*this->param[2]*this->p[m*2+2]/(m+1.))*f;
			deriv[m*2]   += .25*this->param[2]*this->p[m*2+1]/(m+1.)*f;
	}
	deriv[0] += .5*this->param[2]*this->p[1]*f;
}
#endif
