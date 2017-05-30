#include "stdafx.h"

#include "cshapes.h"
#include "ccohes2d.h"

#define  Message(Msg)   { printf("%s", Msg);  printf("\n");}

int CCohes2D::NUM_ADHES = 3;
int CCohes2D::NUM_SHEAR = 6;
int CCohes2D::NUM_SHIFT = 5;
int CCohes2D::MAX_PHASE = 3;
int CCohes2D::NUM_HESS  = 8;
int CCohes2D::regul = 1;
#define __NORMAL_DERIV_FIRST0__
#define n__NORMAL_DERIV_SECOND__

//////////////////////////////////
//...initialization of the blocks;
int CCohes2D::block_shape_init(Block<double> & B, Num_State id_free)
{
	int k, m;
   if (  B.shape && id_free == INITIAL_STATE) delete_shapes(B.shape);
   if (! B.shape && B.mp) {
		B.shape = new CShapeMixer<double>;
		if ((B.type & ERR_CODE) == POLY_BLOCK) {
			B.shape->add_shape(CreateShape<double>(MP2D_POLY_SHAPE));

			B.shape->degree_init1(UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));
			B.shape->set_shape(0, get_param(NUM_MPLS+1)*fabs(B.mp[7]));
		}
		else
		if ((B.type & ERR_CODE) == GRAD_POLY_BLOCK) {
			B.shape->add_shape(CreateShape<double>(MP2D_POLY_SHAPE));
			B.shape->add_shape(CreateShape<double>(SK2D_RADII_POLY_SHAPE));
			B.shape->add_shape(CreateShape<double>(SK2D_RADII_POLY_SHAPE));

			B.shape->degree_init1(UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));
			B.shape->set_shape(0, get_param(NUM_MPLS+1)*fabs(B.mp[7]));
			B.shape->set_shape(1, get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_SHEAR+3+(-B.link[NUM_PHASE]-1)*NUM_SHIFT));
			B.shape->set_shape(2, get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_SHEAR+2+(-B.link[NUM_PHASE]-1)*NUM_SHIFT));
		}
		else
		if ((B.type & ERR_CODE) == ZOOM_BLOCK && B.mp[0] == ID_MAP(1, SPHEROID_GENUS)) {
			B.shape->add_shape(CreateShape<double>(MP2D_POLY_SHAPE));
			B.shape->add_shape(CreateShape<double>(MP2D_ZOOM_SHAPE));

			B.shape->degree_init1(UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));
			B.shape->set_shape(0, get_param(NUM_MPLS+1)*fabs(B.mp[7]));
			B.shape->set_shape(1, fabs(B.mp[8]));
		}
		else
		if ((B.type & ERR_CODE) == GRAD_ZOOM_BLOCK && B.mp[0] == ID_MAP(1, SPHEROID_GENUS)) {
			B.shape->add_shape(CreateShape<double>(MP2D_POLY_SHAPE));
			B.shape->add_shape(CreateShape<double>(SK2D_RADII_POLY_SHAPE));
			B.shape->add_shape(CreateShape<double>(SK2D_RADII_POLY_SHAPE));
			B.shape->add_shape(CreateShape<double>(MP2D_ZOOM_SHAPE));
			B.shape->add_shape(CreateShape<double>(SK2D_RADII_ZOOM_SHAPE));
			B.shape->add_shape(CreateShape<double>(SK2D_RADII_ZOOM_SHAPE));

			B.shape->degree_init1(UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));
			B.shape->set_shape(0, get_param(NUM_MPLS+1)*fabs(B.mp[7]));
			B.shape->set_shape(1, get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_SHEAR+3+(-B.link[NUM_PHASE]-1)*NUM_SHIFT));
			B.shape->set_shape(2, get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_SHEAR+2+(-B.link[NUM_PHASE]-1)*NUM_SHIFT));
			B.shape->set_shape(3, fabs(B.mp[8]));
			B.shape->set_shape(4, fabs(B.mp[8]), get_param(NUM_SHEAR+3+(-B.link[NUM_PHASE]-1)*NUM_SHIFT));
			B.shape->set_shape(5, fabs(B.mp[8]), get_param(NUM_SHEAR+2+(-B.link[NUM_PHASE]-1)*NUM_SHIFT));
		}
		else {
			B.shape->add_shape(CreateShape<double>(MP2D_POLY_SHAPE));
			extern int gradient_model;
			if (gradient_model)	{//...using gradient displacements;
				if ((B.type & ERR_CODE) == ELLI_BLOCK) B.shape->add_shape(CreateShape<double>(SK2D_ELLI_SHAPE));
				else											   B.shape->add_shape(CreateShape<double>(SK2D_POLY_SHAPE/*SK2D_BEAMZ_SHAPE*/));
			}
////////////////////////
//...setting parameters;
			B.shape->set_shape(0, get_param(NUM_MPLS+1)*fabs(B.mp[7]));
			if ((B.type & ERR_CODE) == ELLI_BLOCK) B.shape->set_shape(1, get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_MPLS+1)*fabs(B.mp[7]));
			else												B.shape->set_shape(1, get_param(NUM_MPLS+1)*fabs(B.mp[7])); 
      
			B.shape->degree_init1(UnPackInts(get_param(NUM_MPLS)), solver.id_norm, 2);
			if (B.link[NUM_PHASE] == -2) //...another degree of multipoles for nclusion!!!
			B.shape->degree_init1(UnPackInts(get_param(NUM_MPLS), 1), solver.id_norm, 2);
		}

//////////////////////////////////////////////////////////////////////////////
//...setting acselerator, local system of coordinate and init parametrization;
         B.shape->set_local(B.mp+1);
         B.shape->release  ();
   }

///////////////////////////////////////////////
//...setting cohesion parameter and potentials;
   if (B.shape && id_free != INITIAL_STATE) {
		if (id_free == SPECIAL_STATE) { //...переустановка радиуса и центра мультиполей;
			B.shape->set_local(B.mp+1);
			//if (B.link[NUM_PHASE] == -1) B.mp[4] = M_PI/3.;
			if ((B.type & ERR_CODE) == POLY_BLOCK) B.shape->set_shape(get_param(NUM_MPLS+1)*fabs(B.mp[7]));	else
			if ((B.type & ERR_CODE) == GRAD_POLY_BLOCK) {
				B.shape->set_shape(0, get_param(NUM_MPLS+1)*fabs(B.mp[7]));
				B.shape->set_shape(1, get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_SHEAR+3+(-B.link[NUM_PHASE]-1)*NUM_SHIFT));
				B.shape->set_shape(2, get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_SHEAR+2+(-B.link[NUM_PHASE]-1)*NUM_SHIFT));
			}
			else
			if ((B.type & ERR_CODE) == ZOOM_BLOCK && B.mp[0] == ID_MAP(1, SPHEROID_GENUS)) {
				B.shape->set_shape(0, get_param(NUM_MPLS+1)*fabs(B.mp[7]));
				B.shape->set_shape(1, fabs(B.mp[8]));
			}
			else
			if ((B.type & ERR_CODE) == GRAD_ZOOM_BLOCK && B.mp[0] == ID_MAP(1, SPHEROID_GENUS)) {
				B.shape->set_shape(0, get_param(NUM_MPLS+1)*fabs(B.mp[7]));
				B.shape->set_shape(1, get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_SHEAR+3+(-B.link[NUM_PHASE]-1)*NUM_SHIFT));
				B.shape->set_shape(2, get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_SHEAR+2+(-B.link[NUM_PHASE]-1)*NUM_SHIFT));
				B.shape->set_shape(3, fabs(B.mp[8]));
				B.shape->set_shape(4, fabs(B.mp[8]), get_param(NUM_SHEAR+3+(-B.link[NUM_PHASE]-1)*NUM_SHIFT));
				B.shape->set_shape(5, fabs(B.mp[8]), get_param(NUM_SHEAR+2+(-B.link[NUM_PHASE]-1)*NUM_SHIFT));
			}
			else {
				B.shape->set_shape(NUM_MPLS, get_param(1)*fabs( B.mp[7]));
				if ((B.type & ERR_CODE) == ELLI_BLOCK) B.shape->set_shape(1, get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_MPLS+1)*fabs(B.mp[7]));
				else												B.shape->set_shape(1, get_param(NUM_MPLS+1)*fabs(B.mp[7])); 
			}
		}
		else
		if (id_free == OK_STATE)
			for (m = 0; m < solver.id_norm; m++)
				B.shape->set_potential(solver.hh[k = (int)(&B-B.B)][0][m], m);
		else
		if (id_free == NO_STATE)
			for (m = 0; m < solver.id_norm; m++)
				B.shape->get_potential(solver.hh[k = (int)(&B-B.B)][0][solver.id_norm+m], m);
		else
		if (id_free == NULL_STATE) //...переустановка потенциалов (в случае перемены степени, например);
				B.shape->init_potential();
   }                    
   return(B.shape != NULL);
}

//////////////////////////////////////////////
//...realization of common displacements (Rx);
void CCohes2D::jump1_common_x(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, k, l = 0;
	double alpha = .25/(get_param(NUM_SHEAR+1+shift)-1.), G0 = 1./get_param(NUM_SHEAR+shift),
			 C0 = 1./get_param(NUM_SHEAR-1+shift), * ptr; m += solver.id_norm;

	for (k = 0; k < 2; k++) {
		B[i].shape->cpy_x     (l, B[i].shape->deriv);
		B[i].shape->cpy       (l, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->p_cpy, B[i].shape->deriv, (alpha+1.)*G0, P[0]*alpha*G0);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[1]*alpha*G0, 0.);

		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l, 1));

		if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) {
			jump_admittance	(l, i, m-solver.id_norm, 0.);
			B[i].shape->adm_xy(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), C0);
			B[i].shape->adm_xx(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l), C0);
			B[i].shape->adm	(l, ptr, -G0);
		
			if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) {
				jump_admittance	(l, i, m-solver.id_norm, 0.);
				B[i].shape->adm_xy(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), -C0);
				B[i].shape->adm_xx(l, B[i].shape->FULL(solver.hh[i][0][m], l), -C0);
				l++;
			}
		}
	}
}

////////////////////////////////////////////////////////////////////
//...additional inclusion of cohesion compression displacement (ux);
void CCohes2D::jump1_compress_x(double * P, int i, int m)
{
	int   shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
	double	C0 = 1./get_param(NUM_SHEAR-1+shift); m += solver.id_norm;
	B[i].shape->adm_xy(1, B[i].shape->FULL(solver.hh[i][0][m], 1, 1), -C0);
	B[i].shape->adm_xx(1, B[i].shape->FULL(solver.hh[i][0][m], 1), -C0);
}

///////////////////////////////////////////////
//...realization of classic displacements (Rx);
void CCohes2D::jump1_classic_x(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, k, l = 0;
	double alpha = .25/(get_param(NUM_SHEAR+1+shift)-1.), G0 = 1./get_param(NUM_SHEAR+shift),
			 C0 = 1./get_param(NUM_SHEAR-1+shift); m += solver.id_norm;

	for (k = 0; k < 2; k++) {
		B[i].shape->cpy_x     (l, B[i].shape->deriv);
		B[i].shape->cpy       (l, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->p_cpy, B[i].shape->deriv, (alpha+1.)*G0, P[0]*alpha*G0);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[1]*alpha*G0, 0.);

		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l, 1));

		if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) 
		if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) l++;
	}
}

////////////////////////////////////////////////
//...realization of cohesion displacements (ux);
void CCohes2D::jump1_cohesion_x(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, k, l = 0;
	double alpha = .25/(get_param(NUM_SHEAR+1+shift)-1.), G0 = 1./get_param(NUM_SHEAR+shift),
			 C0 = 1./get_param(NUM_SHEAR-1+shift), * ptr; m += solver.id_norm;

	for (k = 0; k < 2; k++) {
		jump_admittance	(l, i, m-solver.id_norm, 0.);
		if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) {
			jump_admittance	(l, i, m-solver.id_norm, 0.);
			B[i].shape->adm_xy(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), -C0);
			B[i].shape->adm_xx(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l), -C0);
			B[i].shape->adm	(l, ptr, G0);
		
			if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) {
				jump_admittance	(l, i, m-solver.id_norm, 0.);
				B[i].shape->adm_xy(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), C0);
				B[i].shape->adm_xx(l, B[i].shape->FULL(solver.hh[i][0][m], l), C0);
				l++;
			}
		}
	}
}

//////////////////////////////////////////////
//...realization of common displacements (Ry);
void CCohes2D::jump1_common_y(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, k, l = 0;
	double alpha = .25/(get_param(NUM_SHEAR+1+shift)-1.), G0 = 1./get_param(NUM_SHEAR+shift),
			 C0 = 1./get_param(NUM_SHEAR-1+shift), * ptr; m += solver.id_norm;
	for (k = 0; k < 2; k++) {
		B[i].shape->cpy_y     (l, B[i].shape->deriv);
		B[i].shape->cpy       (l, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->p_cpy, B[i].shape->deriv, (alpha+1.)*G0, P[1]*alpha*G0);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[0]*alpha*G0, 0.);

		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], l, 1));

		if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) {
			jump_admittance	(l, i, m-solver.id_norm, 0.);
			B[i].shape->adm_xy(l, B[i].shape->FULL(solver.hh[i][0][m], l), C0);
			B[i].shape->adm_yy(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1), C0);
			B[i].shape->adm	(l, ptr, -G0);

			if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) {
				jump_admittance	(l, i, m-solver.id_norm, 0.);
				B[i].shape->adm_xy(l, B[i].shape->FULL(solver.hh[i][0][m], l), -C0);
				B[i].shape->adm_yy(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), -C0);
				l++;
			}
		}
	}
}

////////////////////////////////////////////////////////////////////
//...additional inclusion of cohesion compression displacement (uy);
void CCohes2D::jump1_compress_y(double * P, int i, int m)
{
	int   shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
	double	C0 = 1./get_param(NUM_SHEAR-1+shift); m += solver.id_norm;
	B[i].shape->adm_xy(1, B[i].shape->FULL(solver.hh[i][0][m], 1), -C0);
	B[i].shape->adm_yy(1, B[i].shape->FULL(solver.hh[i][0][m], 1, 1), -C0);
}

///////////////////////////////////////////////
//...realization of classic displacements (Ry);
void CCohes2D::jump1_classic_y(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, k, l = 0;
	double alpha = .25/(get_param(NUM_SHEAR+1+shift)-1.), G0 = 1./get_param(NUM_SHEAR+shift),
			 C0 = 1./get_param(NUM_SHEAR-1+shift); m += solver.id_norm;
	for (k = 0; k < 2; k++) {
		B[i].shape->cpy_y     (l, B[i].shape->deriv);
		B[i].shape->cpy       (l, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->p_cpy, B[i].shape->deriv, (alpha+1.)*G0, P[1]*alpha*G0);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[0]*alpha*G0, 0.);

		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], l, 1));

		if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) 
		if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) l++;
	}
}

////////////////////////////////////////////////
//...realization of cohesion displacements (uy);
void CCohes2D::jump1_cohesion_y(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, k, l = 0;
	double alpha = .25/(get_param(NUM_SHEAR+1+shift)-1.), G0 = 1./get_param(NUM_SHEAR+shift),
			 C0 = 1./get_param(NUM_SHEAR-1+shift), * ptr; m += solver.id_norm;
	for (k = 0; k < 2; k++) {
		jump_admittance	(l, i, m-solver.id_norm, 0.);
		if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) {
			jump_admittance	(l, i, m-solver.id_norm, 0.);
			B[i].shape->adm_xy(l, B[i].shape->FULL(solver.hh[i][0][m], l), -C0);
			B[i].shape->adm_yy(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1), -C0);
			B[i].shape->adm	(l, ptr, G0);

			if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) {
				jump_admittance	(l, i, m-solver.id_norm, 0.);
				B[i].shape->adm_xy(l, B[i].shape->FULL(solver.hh[i][0][m], l), C0);
				B[i].shape->adm_yy(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), C0);
				l++;
			}
		}
	}
}

////////////////////////////////////////////////////////////////////
//...realization of normal derivative for common displacements (Rx);
void CCohes2D::jump2_common_x(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm, k, l = 0;
	double alpha = .25/(get_param(NUM_SHEAR+1+shift)-1.), G0 = 1./get_param(NUM_SHEAR+shift),
			 C0 = 1./get_param(NUM_SHEAR-1+shift), * ptr; m += solver.id_norm;
	for (k = 0; k < 2; k++) {
		B[i].shape->cpy_xx	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[3], 0.);
		B[i].shape->adm_xy    (l, B[i].shape->deriv, P[4]);

		B[i].shape->cpy       (l, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[0]*alpha*G0, 0.);
		B[i].shape->adm_x     (l, B[i].shape->deriv, P[3]*(2.*alpha+1.)*G0);
		B[i].shape->adm_y     (l, B[i].shape->deriv, P[4]*(1.+alpha)*G0);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, P[1]*alpha*G0, 0.);
		B[i].shape->adm_x     (l, B[i].shape->p_cpy, P[4]*alpha*G0);

		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], l, 1));

		if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) {
			B[i].shape->admittance(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l), B[i].shape->FULL(solver.hh[i][0][num_hess], l), 0., C0);
			B[i].shape->adm_x     (l, ptr, -P[3]*G0);
			B[i].shape->adm_y     (l, ptr, -P[4]*G0);
			B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->FULL(solver.hh[i][0][num_hess], l, 1), 0., C0);

			if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) {
				B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l), B[i].shape->FULL(solver.hh[i][0][num_hess+2], l), 0., -C0);
				B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->FULL(solver.hh[i][0][num_hess+2], l, 1), 0., -C0);
				l++;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//...additional inclusion of normal derivatives of cohesion compression displacement (ux);
void CCohes2D::jump2_compress_x(double * P, int i, int m)
{
	int   shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm;
	double	C0 = 1./get_param(NUM_SHEAR-1+shift); m += solver.id_norm;

	B[i].shape->admittance(1, B[i].shape->FULL(solver.hh[i][0][m], 1), B[i].shape->FULL(solver.hh[i][0][num_hess], 1), 1., -C0);
	B[i].shape->admittance(1, B[i].shape->FULL(solver.hh[i][0][m], 1, 1), B[i].shape->FULL(solver.hh[i][0][num_hess], 1, 1), 1., -C0);
}

/////////////////////////////////////////////////////////////////////
//...realization of normal derivative for classic displacements (Rx);
void CCohes2D::jump2_classic_x(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm, k, l = 0;
	double alpha = .25/(get_param(NUM_SHEAR+1+shift)-1.), G0 = 1./get_param(NUM_SHEAR+shift),
			 C0 = 1./get_param(NUM_SHEAR-1+shift); m += solver.id_norm;
	for (k = 0; k < 2; k++) {
		B[i].shape->cpy_xx	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[3], 0.);
		B[i].shape->adm_xy    (l, B[i].shape->deriv, P[4]);

		B[i].shape->cpy       (l, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[0]*alpha*G0, 0.);
		B[i].shape->adm_x     (l, B[i].shape->deriv, P[3]*(2.*alpha+1.)*G0);
		B[i].shape->adm_y     (l, B[i].shape->deriv, P[4]*(1.+alpha)*G0);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, P[1]*alpha*G0, 0.);
		B[i].shape->adm_x     (l, B[i].shape->p_cpy, P[4]*alpha*G0);

		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], l, 1));

		if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) 
		if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) l++;
	}
}

//////////////////////////////////////////////////////////////////////
//...realization of normal derivative for cohesion displacements (Rx);
void CCohes2D::jump2_cohesion_x(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm, k, l = 0;
	double alpha = .25/(get_param(NUM_SHEAR+1+shift)-1.), G0 = 1./get_param(NUM_SHEAR+shift),
			 C0 = 1./get_param(NUM_SHEAR-1+shift), * ptr; m += solver.id_norm;
	for (k = 0; k < 2; k++) {
		if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) {
			B[i].shape->admittance(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l), B[i].shape->FULL(solver.hh[i][0][num_hess], l), 0., -C0);
			B[i].shape->adm_x     (l, ptr, P[3]*G0);
			B[i].shape->adm_y     (l, ptr, P[4]*G0);
			B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->FULL(solver.hh[i][0][num_hess], l, 1), 0., -C0);

			if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) {
				B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l), B[i].shape->FULL(solver.hh[i][0][num_hess+2], l), 0., C0);
				B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->FULL(solver.hh[i][0][num_hess+2], l, 1), 0., C0);
				l++;
			}
		}
	}
}

////////////////////////////////////////////////////////////////////
//...realization of normal derivative for common displacements (Ry);
void CCohes2D::jump2_common_y(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm, k, l = 0;
	double alpha = .25/(get_param(NUM_SHEAR+1+shift)-1.), G0 = 1./get_param(NUM_SHEAR+shift),
			 C0 = 1./get_param(NUM_SHEAR-1+shift), * ptr; m += solver.id_norm;
	for (k = 0; k < 2; k++) {
		B[i].shape->cpy_yy	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[4], 0.);
		B[i].shape->adm_xy    (l, B[i].shape->deriv, P[3]);

		B[i].shape->cpy       (l, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[1]*alpha*G0, 0.);
		B[i].shape->adm_x     (l, B[i].shape->deriv, P[3]*(1.+alpha)*G0);
		B[i].shape->adm_y     (l, B[i].shape->deriv, P[4]*(2.*alpha+1.)*G0);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, P[0]*alpha*G0, 0.);
		B[i].shape->adm_y     (l, B[i].shape->p_cpy, P[3]*alpha*G0);

		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l, 1));

		if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) {
			B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l), B[i].shape->FULL(solver.hh[i][0][num_hess], l, 1), 0., C0);
			B[i].shape->admittance(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->FULL(solver.hh[i][0][num_hess+1], l), 0., C0);
			B[i].shape->adm_x     (l, ptr, -P[3]*G0);
			B[i].shape->adm_y     (l, ptr, -P[4]*G0);

			if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) {
				B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l), B[i].shape->FULL(solver.hh[i][0][num_hess+2], l, 1), 0., -C0);
				B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->FULL(solver.hh[i][0][num_hess+3], l), 0., -C0);
				l++;
			}
		}		
	}
}

//////////////////////////////////////////////////////////////////////////////////////////
//...additional inclusion of normal derivatives of cohesion compression displacement (uy);
void CCohes2D::jump2_compress_y(double * P, int i, int m)
{
	int   shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm;
	double	C0 = 1./get_param(NUM_SHEAR-1+shift); m += solver.id_norm;

	B[i].shape->admittance(1, B[i].shape->FULL(solver.hh[i][0][m], 1), B[i].shape->FULL(solver.hh[i][0][num_hess], 1, 1), 1., -C0);
	B[i].shape->admittance(1, B[i].shape->FULL(solver.hh[i][0][m], 1, 1), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 1), 1., -C0);
}

/////////////////////////////////////////////////////////////////////
//...realization of normal derivative for classic displacements (Ry);
void CCohes2D::jump2_classic_y(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm, k, l = 0;
	double alpha = .25/(get_param(NUM_SHEAR+1+shift)-1.), G0 = 1./get_param(NUM_SHEAR+shift),
			 C0 = 1./get_param(NUM_SHEAR-1+shift); m += solver.id_norm;
	for (k = 0; k < 2; k++) {
		B[i].shape->cpy_yy	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[4], 0.);
		B[i].shape->adm_xy    (l, B[i].shape->deriv, P[3]);

		B[i].shape->cpy       (l, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[1]*alpha*G0, 0.);
		B[i].shape->adm_x     (l, B[i].shape->deriv, P[3]*(1.+alpha)*G0);
		B[i].shape->adm_y     (l, B[i].shape->deriv, P[4]*(2.*alpha+1.)*G0);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, P[0]*alpha*G0, 0.);
		B[i].shape->adm_y     (l, B[i].shape->p_cpy, P[3]*alpha*G0);

		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l, 1));

		if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) 
		if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) l++;
	}
}

//////////////////////////////////////////////////////////////////////
//...realization of normal derivative for cohesion displacements (Ry);
void CCohes2D::jump2_cohesion_y(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm, k, l = 0;
	double alpha = .25/(get_param(NUM_SHEAR+1+shift)-1.), G0 = 1./get_param(NUM_SHEAR+shift),
			 C0 = 1./get_param(NUM_SHEAR-1+shift), * ptr; m += solver.id_norm;
	for (k = 0; k < 2; k++) {
		if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) {
			B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l), B[i].shape->FULL(solver.hh[i][0][num_hess], l, 1), 0., -C0);
			B[i].shape->admittance(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->FULL(solver.hh[i][0][num_hess+1], l), 0., -C0);
			B[i].shape->adm_x     (l, ptr, P[3]*G0);
			B[i].shape->adm_y     (l, ptr, P[4]*G0);

			if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) {
				B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l), B[i].shape->FULL(solver.hh[i][0][num_hess+2], l, 1), 0., C0);
				B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->FULL(solver.hh[i][0][num_hess+3], l), 0., C0);
				l++;
			}
		}		
	}
}

/////////////////////////////////////////////////////////////////
//...realization of surface forces for common displacements (Px);
void CCohes2D::jump4_common_x(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm, k, l = 0;
	double alpha = .5/(get_param(NUM_SHEAR+1+shift)-1.), G0 = get_param(NUM_SHEAR+shift),
			 C0 = G0/get_param(NUM_SHEAR-1+shift)*2., * ptr; m += solver.id_norm;
	for (k = 0; k < 2; k++) {
		B[i].shape->cpy_xx	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[3], 0.);
		B[i].shape->adm_xy    (l, B[i].shape->deriv, P[4]);

		B[i].shape->admittance(l, B[i].shape->deriv, NULL, alpha, 0.);
		B[i].shape->cpy       (l, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[0], 0.);
		B[i].shape->adm_x     (l, B[i].shape->deriv, P[3]);
		B[i].shape->adm_y     (l, B[i].shape->deriv, P[4]*(1.+alpha));

		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, P[1], 0.);
		B[i].shape->adm_y     (l, B[i].shape->p_cpy, P[3]*(-2.*alpha-1.));
		B[i].shape->adm_x     (l, B[i].shape->p_cpy, P[4]*( 1.+alpha));

		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], l, 1));

		if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) {
			B[i].shape->admittance(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l), B[i].shape->FULL(solver.hh[i][0][num_hess], l), 0., C0);
			B[i].shape->adm_x     (l, ptr, -P[3]*2.);
			B[i].shape->adm_y     (l, ptr, -P[4]);
			B[i].shape->admittance(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->FULL(solver.hh[i][0][num_hess], l, 1), 0., C0);
			B[i].shape->adm_x     (l, ptr, -P[4]);

			if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) {
				double alpha2 = get_param(NUM_SHEAR+1+shift)/(get_param(NUM_SHEAR+1+shift)-1.);
				B[i].shape->admittance(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l), B[i].shape->FULL(solver.hh[i][0][num_hess+2], l), 0., -C0);
				B[i].shape->adm_x     (l, ptr, P[3]*alpha2);
				B[i].shape->admittance(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->FULL(solver.hh[i][0][num_hess+2], l, 1), 0., -C0);
				B[i].shape->adm_y     (l, ptr, P[3]*alpha2);
				l++;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
//...additional inclusion of surface forces for compression displacement (ux);
void CCohes2D::jump4_compress_x(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm;
	double alpha = get_param(NUM_SHEAR+1+shift)/(get_param(NUM_SHEAR+1+shift)-1.), G0 = get_param(NUM_SHEAR+shift),
			 C0 = G0/get_param(NUM_SHEAR-1+shift)*2., * ptr; m += solver.id_norm;

	B[i].shape->admittance(1, ptr = B[i].shape->FULL(solver.hh[i][0][m], 1), B[i].shape->FULL(solver.hh[i][0][num_hess], 1), 1., -C0);
	B[i].shape->adm_x     (1, ptr, P[3]*alpha);
	B[i].shape->admittance(1, ptr = B[i].shape->FULL(solver.hh[i][0][m], 1, 1), B[i].shape->FULL(solver.hh[i][0][num_hess], 1, 1), 1., -C0);
	B[i].shape->adm_y     (1, ptr, P[3]*alpha);
}

//////////////////////////////////////////////////////////////////
//...realization of surface forces for classic displacements (Px);
void CCohes2D::jump4_classic_x(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm, k, l = 0;
	double alpha = .5/(get_param(NUM_SHEAR+1+shift)-1.), G0 = get_param(NUM_SHEAR+shift),
			 C0 = G0/get_param(NUM_SHEAR-1+shift)*2.; m += solver.id_norm;
	for (k = 0; k < 2; k++) {
		B[i].shape->cpy_xx	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[3], 0.);
		B[i].shape->adm_xy    (l, B[i].shape->deriv, P[4]);

		B[i].shape->admittance(l, B[i].shape->deriv, NULL, alpha, 0.);
		B[i].shape->cpy       (l, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[0], 0.);
		B[i].shape->adm_x     (l, B[i].shape->deriv, P[3]);
		B[i].shape->adm_y     (l, B[i].shape->deriv, P[4]*(1.+alpha));

		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, P[1], 0.);
		B[i].shape->adm_y     (l, B[i].shape->p_cpy, P[3]*(-2.*alpha-1.));
		B[i].shape->adm_x     (l, B[i].shape->p_cpy, P[4]*( 1.+alpha));

		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], l, 1));

		if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) 
		if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) l++;
	}
}

/////////////////////////////////////////////////////////////////
//...realization of surface forces for common displacements (Py);
void CCohes2D::jump4_common_y(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm, k, l = 0;
	double alpha = .5/(get_param(NUM_SHEAR+1+shift)-1.), G0 = get_param(NUM_SHEAR+shift),
			 C0 = G0/get_param(NUM_SHEAR-1+shift)*2., * ptr; m += solver.id_norm;
	for (k = 0; k < 2; k++) {
		B[i].shape->cpy_yy	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[4], 0.);
		B[i].shape->adm_xy    (l, B[i].shape->deriv, P[3]);

		B[i].shape->admittance(l, B[i].shape->deriv, NULL, alpha, 0.);
		B[i].shape->cpy       (l, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[1], 0.);
		B[i].shape->adm_x     (l, B[i].shape->deriv, P[3]*(1.+alpha));
		B[i].shape->adm_y     (l, B[i].shape->deriv, P[4]);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, P[0], 0.);
		B[i].shape->adm_x     (l, B[i].shape->p_cpy, P[4]*(-2.*alpha-1.));
		B[i].shape->adm_y     (l, B[i].shape->p_cpy, P[3]*( 1.+alpha));

		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l, 1));

		if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) {
			B[i].shape->admittance(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l), B[i].shape->FULL(solver.hh[i][0][num_hess], l, 1), 0., C0);
			B[i].shape->adm_y     (l, ptr, -P[3]);
			B[i].shape->admittance(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->FULL(solver.hh[i][0][num_hess+1], l), 0., C0);
			B[i].shape->adm_x     (l, ptr, -P[3]);
			B[i].shape->adm_y     (l, ptr, -P[4]*2.);

			if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) {
				double alpha2 = get_param(NUM_SHEAR+1+shift)/(get_param(NUM_SHEAR+1+shift)-1.);
				B[i].shape->admittance(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l), B[i].shape->FULL(solver.hh[i][0][num_hess+2], l, 1), 0., -C0);
				B[i].shape->adm_x     (l, ptr, P[4]*alpha2);
				B[i].shape->admittance(l, ptr = B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->FULL(solver.hh[i][0][num_hess+3], l), 0., -C0);
				B[i].shape->adm_y     (l, ptr, P[4]*alpha2);
				l++;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
//...additional inclusion of surface forces for compression displacement (uy);
void CCohes2D::jump4_compress_y(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm;
	double alpha = get_param(NUM_SHEAR+1+shift)/(get_param(NUM_SHEAR+1+shift)-1.), G0 = get_param(NUM_SHEAR+shift),
			 C0 = G0/get_param(NUM_SHEAR-1+shift)*2., * ptr; m += solver.id_norm;

	B[i].shape->admittance(1, ptr = B[i].shape->FULL(solver.hh[i][0][m], 1), B[i].shape->FULL(solver.hh[i][0][num_hess], 1, 1), 1., -C0);
	B[i].shape->adm_x     (1, ptr, P[4]*alpha);
	B[i].shape->admittance(1, ptr = B[i].shape->FULL(solver.hh[i][0][m], 1, 1), B[i].shape->FULL(solver.hh[i][0][num_hess+1], 1), 1., -C0);
	B[i].shape->adm_y     (1, ptr, P[4]*alpha);
}

//////////////////////////////////////////////////////////////////
//...realization of surface forces for classic displacements (Py);
void CCohes2D::jump4_classic_y(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, num_hess = NUM_HESS+solver.id_norm, k, l = 0;
	double alpha = .5/(get_param(NUM_SHEAR+1+shift)-1.), G0 = get_param(NUM_SHEAR+shift),
			 C0 = G0/get_param(NUM_SHEAR-1+shift)*2.; m += solver.id_norm;
	for (k = 0; k < 2; k++) {
		B[i].shape->cpy_yy	 (l, B[i].shape->deriv);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[4], 0.);
		B[i].shape->adm_xy    (l, B[i].shape->deriv, P[3]);

		B[i].shape->admittance(l, B[i].shape->deriv, NULL, alpha, 0.);
		B[i].shape->cpy       (l, B[i].shape->deriv, B[i].shape->p_cpy);
		B[i].shape->admittance(l, B[i].shape->deriv, NULL, P[1], 0.);
		B[i].shape->adm_x     (l, B[i].shape->deriv, P[3]*(1.+alpha));
		B[i].shape->adm_y     (l, B[i].shape->deriv, P[4]);
		B[i].shape->admittance(l, B[i].shape->p_cpy, NULL, P[0], 0.);
		B[i].shape->adm_x     (l, B[i].shape->p_cpy, P[4]*(-2.*alpha-1.));
		B[i].shape->adm_y     (l, B[i].shape->p_cpy, P[3]*( 1.+alpha));

		B[i].shape->cpy(l, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], l));
		B[i].shape->cpy(l, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], l, 1));

		if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) 
		if (B[i].shape->type(++l) == SK2D_RADII_POLY_SHAPE || B[i].shape->type(l) == SK2D_RADII_ZOOM_SHAPE || B[i].shape->type(l) == SK2D_POLY_SHAPE) l++;
	}
}

///////////////////////////////////////////////
//...realization normal derivatives of hessian;
void CCohes2D::hessian_deriv_N(int k, double * P, int i, int id_dop)
{
	int num_hess = NUM_HESS+solver.id_norm+id_dop*2;
	B[i].shape->set_norm_cs(P);

	B[i].shape->cpy_xx(k);
	B[i].shape->deriv_N(k);
	B[i].shape->cpy_xx(k);
	B[i].shape->cpy(k, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess], k));

	B[i].shape->cpy_xy(k);
	B[i].shape->deriv_N(k);
	B[i].shape->cpy_xy(k);
	B[i].shape->cpy(k, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess], k, 1));

	B[i].shape->cpy_yy(k);
	B[i].shape->deriv_N(k);
	B[i].shape->cpy_yy(k);
	B[i].shape->cpy(k, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][num_hess+1], k));
}

///////////////////////////////////////
//...composition of collocation vector;
void CCohes2D::jump_admittance(int l, int i, int m, double adm_re, int k, double adm_im)
{
	m += solver.id_norm;
	if (m >= 0 && adm_re != 1.) {
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l),    NULL, adm_re, 0.);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), NULL, adm_re, 0.);
	}
	if (k >= 0) {
		k += solver.id_norm;
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l),	 B[i].shape->FULL(solver.hh[i][0][k], l),	   1., adm_im);
		B[i].shape->admittance(l, B[i].shape->FULL(solver.hh[i][0][m], l, 1), B[i].shape->FULL(solver.hh[i][0][k], l, 1), 1., adm_im);
	}
}

//////////////////////////////////////////////
//...transformation of the collocation vector;
void CCohes2D::jump_make_local(int i, int m)
{
	if (! solver.dim || ! solver.hh[i][0].GetMatrix()) return;
	m  += solver.id_norm;
	for (int j = 0; j < solver.dim[i]; j++) {
		double P[] = { solver.hh[i][0][m  ][j], 
							solver.hh[i][0][m+1][j], 0. };
		B[i].shape->norm_local(P);
		solver.hh[i][0][m  ][j] = P[0];
		solver.hh[i][0][m+1][j] = P[1];
	}
}

void CCohes2D::jump_make_common(int i, int m)
{
	if (! solver.dim || ! solver.hh[i][0].GetMatrix()) return;
	m  += solver.id_norm;
	for (int j = 0; j < solver.dim[i]; j++) {
		double P[] = { solver.hh[i][0][m  ][j], 
							solver.hh[i][0][m+1][j], 0.};
		B[i].shape->norm_common(P);
		solver.hh[i][0][m  ][j] = P[0];
		solver.hh[i][0][m+1][j] = P[1];
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////
//...формирование матрицы Грама с учетом функционала энергии (условие периодического скачка);
Num_State CCohes2D::gram3(CGrid * nd, int i, int id_local)
{
	if (nd && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double G1 = get_param(NUM_SHEAR), AX, AY, f, P[6], TX, TY, hx;
      int	 m  = solver.id_norm, id_dir, k, j, shift;

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			AX = nd->get_param(0, l); 
			AY = nd->get_param(1, l);
			f  = nd->get_param(3, l);

			for (k  = (int)nd->get_param(4, l), j = 0; j < solver.JR[i][0]; j++) 
			if ( k == solver.JR[i][j+solver.JR_SHIFT]) {
				P[0] = nd->X[l]; P[3] = nd->nX[l];
				P[1] = nd->Y[l]; P[4] = nd->nY[l];
				P[2] = 0.;       P[5] = 0.;
				B[i].shape->make_local(P);
				B[i].shape->norm_local(P+3);

////////////////////////////////////////////////////////////////////////////////////////
//...вычисляем граничные условия периодического скачка и сдвигаем соответственные блоки;
				TX = TY = hx = 0.;
				switch (abs(id_dir = (int)nd->get_param(2, l))) {
					case 1: TX =  AX; hx = -G1*AX; break;
					case 2: TX = -AX; hx =  G1*AX; break;
					case 3: TY =  AY; break;
					case 4: TY = -AY; break;
				}
				B[k].mp[1] -= TX;
				B[k].mp[2] -= TY; B[k].shape->set_local_P0(B[k].mp+1);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m; num < solver.n; num++) {
					 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}

////////////////////////////////////////////////////////////
//...jump of all neaded displacement moments for this block;
				shift = ( -B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, i);
				jump1_common_x (P, i, 0);
				jump1_common_y (P, i, 1);

				jump4_common_x (P, i, 2);
				jump4_common_y (P, i, 3);

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump1_compress_x(P, i, 0); solver.admittance(i, 0, G1);
				jump1_compress_y(P, i, 1); solver.admittance(i, 1, G1);

				jump4_compress_x(P, i, 2); 
				jump4_compress_y(P, i, 3); 

				jump_make_common(i, 0);
				jump_make_common(i, 2);

////////////////////////////////////////////////////////////////////
//...jump of all neaded displacement moments for neighbouring block;
				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);
				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				shift = ( -B[k].link[NUM_PHASE]-1)*NUM_SHIFT;
				B[k].shape->set_shape(1, B[k].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[k].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, k);
				jump1_common_x (P, k, 0);
				jump1_common_y (P, k, 1);

				jump4_common_x (P, k, 2);
				jump4_common_y (P, k, 3);

				B[k].shape->set_shape(1, B[k].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[k].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, k);
				jump1_compress_x(P, k, 0); solver.admittance(k, 0, G1);
				jump1_compress_y(P, k, 1); solver.admittance(k, 1, G1);

				jump4_compress_x(P, k, 2); 
				jump4_compress_y(P, k, 3); 

				jump_make_common(k, 0);
				jump_make_common(k, 2);

////////////////////////////////////////////////////////////////////////////////////
//...условие скачка для классической составляющей поля методом наименьших квадратов;
				solver.admittance(i, 4, 0., 0, 1.); jump_admittance(1, i, 4, 0.);
				solver.admittance(k, 4, 0., 0, 1.); jump_admittance(1, k, 4, 0.); 
				solver.to_transferTR(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);

				solver.admittance(i, 5, 0., 1, 1.); jump_admittance(1, i, 5, 0.);
				solver.admittance(k, 5, 0., 1, 1.); jump_admittance(1, k, 5, 0.); 
				solver.to_transferTR(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);

				if (fabs(hx) > EE) {
				  solver.to_equationHH(i, 0, solver.hh[i][0][m+4],  hx*f);
				  solver.to_equationHH(i, 1, solver.hh[i][0][m+5],  hx*f);

				  solver.to_equationHL(k, 0, solver.hh[k][0][m+4], -hx*f);
				  solver.to_equationHL(k, 1, solver.hh[k][0][m+5], -hx*f);
				}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//...граничное условие для классической составляющей поля перемещений (регуляризация в случае симметричного включения);
				if (regul) {
					if (id_dir == 1 || id_dir == 2) {
						solver.to_equationDD(i, solver.hh[i][0][m+4], solver.hh[i][0][m+4], f);
						solver.to_equationHH(i, 0, solver.hh[i][0][m+4], hx*f*.5);

						solver.to_equationDD(k, solver.hh[k][0][m+4], solver.hh[k][0][m+4], f);
						solver.to_equationHH(k, 0, solver.hh[k][0][m+4], -hx*f*.5);
					}
					if (id_dir == 3 || id_dir == 4) {
						solver.to_equationDD(i, solver.hh[i][0][m+5], solver.hh[i][0][m+5], f);
						solver.to_equationDD(k, solver.hh[k][0][m+5], solver.hh[k][0][m+5], f);
					}
					//аналогичная регуляризация второй задачи???
				}

//////////////////////////////////////////////////////////
//...сшивка когезионого поля методом наименьших квадратов;
				solver.admittance(i, 4, 0., 0, 1.); jump_admittance(0, i, 4, 0.);
				solver.admittance(k, 4, 0., 0, 1.); jump_admittance(0, k, 4, 0.); 
				solver.to_transferTR(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);

				solver.admittance(i, 5, 0., 1, 1.); jump_admittance(0, i, 5, 0.);
				solver.admittance(k, 5, 0., 1, 1.); jump_admittance(0, k, 5, 0.); 
				solver.to_transferTR(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);

/////////////////////////////////////////////////
//...энергетические слагаемые (нормированные G1);
				solver.to_equationER(i, solver.hh[i][0][m+4], solver.hh[i][0][m+2], -f); jump_admittance(1, i, 2, 0.);
				solver.to_equationER(i, solver.hh[i][0][m+5], solver.hh[i][0][m+3], -f); jump_admittance(1, i, 3, 0.); 

				solver.to_equationER(i, solver.hh[i][0][m],	 solver.hh[i][0][m+2], -f);
				solver.to_equationER(i, solver.hh[i][0][m+1], solver.hh[i][0][m+3], -f);

				solver.to_equationEL(k, solver.hh[k][0][m+4], solver.hh[k][0][m+2],  f); jump_admittance(1, k, 2, 0.);
				solver.to_equationEL(k, solver.hh[k][0][m+5], solver.hh[k][0][m+3],  f); jump_admittance(1, k, 3, 0.); 

				solver.to_equationEL(k, solver.hh[k][0][m],	 solver.hh[k][0][m+2],  f);
				solver.to_equationEL(k, solver.hh[k][0][m+1], solver.hh[k][0][m+3],  f);
	
				B[k].mp[1] += TX;
				B[k].mp[2] += TY; B[k].shape->set_local_P0(B[k].mp+1);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

//////////////////////////////////////////////////////////////////////////////
//...формирование матриц перехода с учетом функционала энергии на границе фаз;
Num_State CCohes2D::transfer3(CGrid * nd, int i, int k, int id_local)
{
	if (nd && 0 <= i && i < solver.N) {
      double G1 = get_param(NUM_SHEAR), f, P[6], Ai, Ak, Bi, Bk;
      int m = solver.id_norm, shift;

////////////////////////////////////////
//...коэффициенты поверхностной адгезии;
		Ai = Ak = get_param(NUM_ADHES)*.5-(Bi = Bk = get_param(NUM_ADHES+1)*.5);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int j = 0; j < solver.JR[i][0]; j++) if (k == solver.JR[i][j+solver.JR_SHIFT]) {
			for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
				P[0] = nd->X[l];  P[3] = nd->nX[l];
				P[1] = nd->Y[l];  P[4] = nd->nY[l];
				P[2] = 0.;        P[5] = 0.;
				f = nd->get_param(0, l);
				B[i].shape->make_local(P);
				B[i].shape->norm_local(P+3);

				for (int num = m; num < solver.n; num++) {
					 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}

////////////////////////////////////////////////////////////
//...jump of all neaded displacement moments for this block;
				shift = ( -B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, i);
				jump1_common_x (P, i, 0);
				jump1_common_y (P, i, 1);

				jump4_common_x (P, i, 2);
				jump4_common_y (P, i, 3);

				extern int gradient_model;
				if (gradient_model)	{//...using gradient displacements;
					jump2_common_x (P, i, 4);
					jump2_common_y (P, i, 5);
				}
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump1_compress_x(P, i, 0); solver.admittance(i, 0, G1);
				jump1_compress_y(P, i, 1); solver.admittance(i, 1, G1);

				jump4_compress_x(P, i, 2); 
				jump4_compress_y(P, i, 3); 

				extern int gradient_model;
				if (gradient_model)	{//...using gradient displacements;
					jump2_compress_x(P, i, 4); solver.admittance(i, 4, G1);
					jump2_compress_y(P, i, 5); solver.admittance(i, 5, G1); 
				}
				jump_make_common(i, 0);
				jump_make_common(i, 2);
				jump_make_common(i, 4);

////////////////////////////////////////////////////////////////////
//...jump of all neaded displacement moments for neighbouring block;
				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);
				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				shift = ( -B[k].link[NUM_PHASE]-1)*NUM_SHIFT;
				B[k].shape->set_shape(1, B[k].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[k].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, k);
				jump1_common_x (P, k, 0);
				jump1_common_y (P, k, 1);

				jump4_common_x (P, k, 2);
				jump4_common_y (P, k, 3);

				extern int gradient_model;
				if (gradient_model)	{//...using gradient displacements;
					jump2_common_x (P, k, 4);
					jump2_common_y (P, k, 5);
				}
				B[k].shape->set_shape(1, B[k].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[k].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, k);
				jump1_compress_x(P, k, 0); solver.admittance(k, 0, G1);
				jump1_compress_y(P, k, 1); solver.admittance(k, 1, G1);

				jump4_compress_x(P, k, 2); 
				jump4_compress_y(P, k, 3); 

				extern int gradient_model;
				if (gradient_model)	{//...using gradient displacements;
					jump2_compress_x(P, k, 4); solver.admittance(k, 4, G1); 
					jump2_compress_y(P, k, 5); solver.admittance(k, 5, G1); 
				}
				jump_make_common(k, 0);
				jump_make_common(k, 2);
				jump_make_common(k, 4);

//////////////////////////////////////////////////////////////////////////
//...сшивка функций и нормальных производных методом наименьших квадратов;
				solver.to_transferTR(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);

////////////////////////////////////////////////////////////////
//...добавляем поверхностную энергию адгезии (нормированную G1);
				extern int gradient_model;
				if (gradient_model)	{//...using gradient displacements;
					solver.to_equationER(i, solver.hh[i][0][m+4], solver.hh[i][0][m+4], -f*(Ai*sqr(P[3])+Bi*sqr(P[4])));
					solver.to_equationEL(k, solver.hh[k][0][m+4], solver.hh[k][0][m+4], -f*(Ak*sqr(P[3])+Bk*sqr(P[4])));

					solver.to_equationER(i, solver.hh[i][0][m+4], solver.hh[i][0][m+5], -f*Ai*Bi*P[3]*P[4]*2.);
					solver.to_equationEL(k, solver.hh[k][0][m+4], solver.hh[k][0][m+5], -f*Ak*Bk*P[3]*P[4]*2.);

					solver.to_equationER(i, solver.hh[i][0][m+5], solver.hh[i][0][m+5], -f*(Ai*sqr(P[4])+Bi*sqr(P[3])));
					solver.to_equationEL(k, solver.hh[k][0][m+5], solver.hh[k][0][m+5], -f*(Ak*sqr(P[4])+Bk*sqr(P[3])));
				}

/////////////////////////////////////////////////
//...энергетические слагаемые (нормированные G1);
				if (B[i].link[NUM_PHASE] != B[k].link[NUM_PHASE] && B[i].link[NUM_PHASE] < -1) f = -f;
				solver.admittance(i, 4, 0., 0, 1.); jump_admittance(0, i, 4, 0.);
				solver.admittance(i, 5, 0., 1, 1.); jump_admittance(0, i, 5, 0.);
				solver.to_equationER(i, solver.hh[i][0][m+4], solver.hh[i][0][m+2], -f); jump_admittance(1, i, 2, 0.);
				solver.to_equationER(i, solver.hh[i][0][m+5], solver.hh[i][0][m+3], -f); jump_admittance(1, i, 3, 0.); 

				solver.to_equationER(i, solver.hh[i][0][m],	 solver.hh[i][0][m+2], -f);
				solver.to_equationER(i, solver.hh[i][0][m+1], solver.hh[i][0][m+3], -f);

				solver.admittance(k, 4, 0., 0, 1.); jump_admittance(0, k, 4, 0.); 
				solver.admittance(k, 5, 0., 1, 1.); jump_admittance(0, k, 5, 0.); 
				solver.to_equationEL(k, solver.hh[k][0][m+4], solver.hh[k][0][m+2],  f); jump_admittance(1, k, 2, 0.);
				solver.to_equationEL(k, solver.hh[k][0][m+5], solver.hh[k][0][m+3],  f); jump_admittance(1, k, 3, 0.); 

				solver.to_equationEL(k, solver.hh[k][0][m],	 solver.hh[k][0][m+2],  f);
				solver.to_equationEL(k, solver.hh[k][0][m+1], solver.hh[k][0][m+3],  f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

///////////////////////////////////////////////////////////////////////////
//...inclusion of the boundary condition data to the solver for all blocks;
Num_State CCohes2D::gram4(CGrid * nd, int i, int id_local)
{
	if (nd && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double G1 = get_param(NUM_SHEAR), hx, hy, p3, f, P[6];
		int 	 m  = solver.id_norm, shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;

////////////////////////////////////////////////////////////////////
//...realization of pattern cell boundary condition for Gram matrix;
//		nd->TestGrid("nodes.bln", 0.0005, 0., 0., 0., AXIS_Z, 1);
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = nd->nX[l];
			P[1] = nd->Y[l];  P[4] = nd->nY[l];
			P[2] = 0.;        P[5] = 0.;
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);
			
			hx = nd->get_param(0, l);
			hy = nd->get_param(1, l);
			p3 = nd->get_param(2, l);
			f  = nd->get_param(3, l);

			if (p3 == MIN_HIT || p3 == 2.) {
			  hy  = nd->nY[l]*hx;
			  hx *= nd->nX[l];
			}
			else 
			if (p3 == NUMS_BND) { //...специальный случай -- одноосное растяжение;
				double nju = get_param(NUM_SHEAR+1), G1 = .5/get_param(NUM_SHEAR), 
						 AAA = -nju*G1, BBB = (1.-nju)*G1;
				hx = nd->X[l]*BBB;
				hy = nd->Y[l]*AAA; 
			}
			else 
			if (p3 == 0.) hx = hy = 0.;

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < solver.n; num++)
				memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));

/////////////////////////////////////////////
//...jump of all neaded displacement moments;
			B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
			B[i].shape->parametrization_hess(P, 1);

			hessian_deriv_N(1, P, i);
			jump1_common_x (P, i, 0);
			jump1_common_y (P, i, 1);

			jump4_common_x (P, i, 2);
			jump4_common_y (P, i, 3);

			B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
			B[i].shape->parametrization_hess(1, P, 1);

			hessian_deriv_N (1, P, i);
			jump1_compress_x(P, i, 0); solver.admittance(i, 0, G1);
			jump1_compress_y(P, i, 1); solver.admittance(i, 1, G1);

			jump4_compress_x(P, i, 2); 
			jump4_compress_y(P, i, 3); 

			jump_make_common(i, 0);
			jump_make_common(i, 2);

////////////////////////////////////////////////////////////////////////////////
//...граничные условия методом наименьших квадратов (классическая составляющая);
			if (p3 == MIN_HIT || p3 == NUMS_BND || p3 == MAX_HIT) {
				solver.admittance(i, 4, 0., 0, 1.); jump_admittance(1, i, 4, 0.);
				solver.admittance(i, 5, 0., 1, 1.); jump_admittance(1, i, 5, 0.); 

				if (fabs(hx) > EE)
				solver.to_equationHH(i, 0, solver.hh[i][0][m+4], hx*G1*f);
				solver.to_equationDD(i, solver.hh[i][0][m+4], solver.hh[i][0][m+4], f);

				if (fabs(hy) > EE)
				solver.to_equationHH(i, 0, solver.hh[i][0][m+5], hy*G1*f);
				solver.to_equationDD(i, solver.hh[i][0][m+5], solver.hh[i][0][m+5], f);
			}

///////////////////////////////////////////////////
//...обнуление когезионного поля (на всей границе);
			solver.admittance(i, 4, 0., 0, 1.); jump_admittance(0, i, 4, 0.);
			solver.admittance(i, 5, 0., 1, 1.); jump_admittance(0, i, 5, 0.); 
			solver.to_equationDD(i, solver.hh[i][0][m+4], solver.hh[i][0][m+4], f);
			solver.to_equationDD(i, solver.hh[i][0][m+5], solver.hh[i][0][m+5], f);

/////////////////////////////////////////////////
//...энергетические слагаемые (нормированные G1);
			solver.to_equationEE(i, solver.hh[i][0][m+4], solver.hh[i][0][m+2], -f); jump_admittance(1, i, 2, 0.);
			solver.to_equationEE(i, solver.hh[i][0][m+5], solver.hh[i][0][m+3], -f); jump_admittance(1, i, 3, 0.); 

			solver.to_equationEE(i, solver.hh[i][0][m],	 solver.hh[i][0][m+2], -f);
			solver.to_equationEE(i, solver.hh[i][0][m+1], solver.hh[i][0][m+3], -f);
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

///////////////////////////////////////////////////////////////
//...формирование матриц перехода с учетом функционала энергии;
Num_State CCohes2D::transfer4(CGrid * nd, int i, int k, int id_local)
{
	if (nd && 0 <= i && i < solver.N) {
      double G1 = get_param(NUM_SHEAR), f, P[6];
      int m = solver.id_norm, shift;

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int j = 0; j < solver.JR[i][0]; j++) if (k == solver.JR[i][j+solver.JR_SHIFT]) {
			for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
				P[0] = nd->X[l];  P[3] = nd->nX[l];
				P[1] = nd->Y[l];  P[4] = nd->nY[l];
				P[2] = 0.;        P[5] = 0.;
				B[i].shape->make_local(P);
				B[i].shape->norm_local(P+3);
				f = nd->get_param(0, l);

				for (int num = m; num < solver.n; num++) {
					 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}

////////////////////////////////////////////////////////////
//...jump of all neaded displacement moments for this block;
				shift = ( -B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, i);
				jump1_common_x (P, i, 0);
				jump1_common_y (P, i, 1);

				jump4_common_x (P, i, 2);
				jump4_common_y (P, i, 3);

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump1_compress_x(P, i, 0); solver.admittance(i, 0, G1);
				jump1_compress_y(P, i, 1); solver.admittance(i, 1, G1);

				jump4_compress_x(P, i, 2); 
				jump4_compress_y(P, i, 3); 

				jump_make_common(i, 0);
				jump_make_common(i, 2);

////////////////////////////////////////////////////////////////////
//...jump of all neaded displacement moments for neighbouring block;
				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);
				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				shift = ( -B[k].link[NUM_PHASE]-1)*NUM_SHIFT;
				B[k].shape->set_shape(1, B[k].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[k].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, k);
				jump1_common_x (P, k, 0);
				jump1_common_y (P, k, 1);

				jump4_common_x (P, k, 2);
				jump4_common_y (P, k, 3);

				B[k].shape->set_shape(1, B[k].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[k].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, k);
				jump1_compress_x(P, k, 0); solver.admittance(k, 0, G1);
				jump1_compress_y(P, k, 1); solver.admittance(k, 1, G1);

				jump4_compress_x(P, k, 2); 
				jump4_compress_y(P, k, 3); 

				jump_make_common(k, 0);
				jump_make_common(k, 2);

/////////////////////////////////////////////////////////////////////////////
//...сшивка функций методом наименьших квадратов (классическая составляющая);
				solver.admittance(i, 4, 0., 0, 1.); jump_admittance(1, i, 4, 0.);
				solver.admittance(k, 4, 0., 0, 1.); jump_admittance(1, k, 4, 0.); 
				solver.to_transferTR(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);

				solver.admittance(i, 5, 0., 1, 1.); jump_admittance(1, i, 5, 0.);
				solver.admittance(k, 5, 0., 1, 1.); jump_admittance(1, k, 5, 0.); 
				solver.to_transferTR(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);

////////////////////////////////////////////////////////////////////
//...сшивка функций методом наименьших квадратов (когезионное поле);
				solver.admittance(i, 4, 0., 0, 1.); jump_admittance(0, i, 4, 0.);
				solver.admittance(k, 4, 0., 0, 1.); jump_admittance(0, k, 4, 0.); 
				solver.to_transferTR(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);

				solver.admittance(i, 5, 0., 1, 1.); jump_admittance(0, i, 5, 0.);
				solver.admittance(k, 5, 0., 1, 1.); jump_admittance(0, k, 5, 0.); 
				solver.to_transferTR(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);

/////////////////////////////////////////////////
//...энергетические слагаемые (нормированные G1);
				solver.to_equationER(i, solver.hh[i][0][m+4], solver.hh[i][0][m+2], -f); jump_admittance(1, i, 2, 0.);
				solver.to_equationER(i, solver.hh[i][0][m+5], solver.hh[i][0][m+3], -f); jump_admittance(1, i, 3, 0.); 

				solver.to_equationER(i, solver.hh[i][0][m],	 solver.hh[i][0][m+2], -f);
				solver.to_equationER(i, solver.hh[i][0][m+1], solver.hh[i][0][m+3], -f);

				solver.to_equationEL(k, solver.hh[k][0][m+4], solver.hh[k][0][m+2],  f); jump_admittance(1, k, 2, 0.);
				solver.to_equationEL(k, solver.hh[k][0][m+5], solver.hh[k][0][m+3],  f); jump_admittance(1, k, 3, 0.); 

				solver.to_equationEL(k, solver.hh[k][0][m],	 solver.hh[k][0][m+2],  f);
				solver.to_equationEL(k, solver.hh[k][0][m+1], solver.hh[k][0][m+3],  f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

//////////////////////////////////////////////////////////////////////////
//...интегрирование НДС на заданном наборе узлов для периодической задачи;
Num_State CCohes2D::rigidy1(CGrid * nd, int i, double * K)
{
	if (nd) {
      int l, shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, m = solver.id_norm;
      double alpha = get_param(NUM_SHEAR+shift+1)/(.5-get_param(NUM_SHEAR+shift+1)), 
					 G0 = get_param(NUM_SHEAR+shift), f, P[6], UX, UY, RX, RY;

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		nd->TestGrid("nodes.bln", 0.0007, 0., 0., 0., AXIS_Z, 1);
		for ( l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = nd->nX[l];
			P[1] = nd->Y[l];  P[4] = nd->nY[l];
			P[2] = 0.;        P[5] = 0.;
			B[i].shape->make_local(P);
			B[i].shape->norm_local(P+3);
			f = nd->get_param(0, l);

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < solver.n; num++)
				memset(solver.hh[0][0][num], 0, solver.dim[0]*sizeof(double));

///////////////////////////////////////////////////////////////////////////////////////////////
//...формирование интеграла от полных деформаций и классических напряжений (два состояния НДС);
			B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
 			B[i].shape->parametrization_grad(0, P, 1);
			B[i].shape->parametrization_grad(1, P, 2);

			jump1_common_x(P, i, 0); 
			jump1_common_y(P, i, 1); 

			B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
			B[i].shape->parametrization_grad(1, P, 2);

			jump1_compress_x(P, i, 0);
			jump1_compress_y(P, i, 1);
			jump_make_common(i, 0); //...переходим к общей системе координат;

//////////////////////////////////
//...вычисляем полные перемещения;
			RX = B[i].shape->potential(solver.hh[i][0][m],   0);
			RY = B[i].shape->potential(solver.hh[i][0][m+1], 0);
			
			K[3] -= RX*nd->nX[l]*f;
			K[4] -=(RX*nd->nY[l]+RY*nd->nX[l])*f*.5;
			K[5] -= RY*nd->nY[l]*f;

			RX = B[i].shape->potential(solver.hh[i][0][m],   1);
			RY = B[i].shape->potential(solver.hh[i][0][m+1], 1);
			
			K[9]  -= RX*nd->nX[l]*f;
			K[10] -=(RX*nd->nY[l]+RY*nd->nX[l])*f*.5;
			K[11] -= RY*nd->nY[l]*f;

////////////////////////////////////////
//...вычисляем классические перемещения;
			jump_admittance(1, i, 0, 0.);
			jump_admittance(1, i, 1, 0.);

			UX = B[i].shape->potential(solver.hh[i][0][m],   0);
			UY = B[i].shape->potential(solver.hh[i][0][m+1], 0);
			
			K[0] -= G0*(UX*nd->nX[l]*2.+alpha*(UX*nd->nX[l]+UY*nd->nY[l]))*f;
			K[1] -= G0*(UX*nd->nY[l]+UY*nd->nX[l])*f;
			K[2] -= G0*(UY*nd->nY[l]*2.+alpha*(UX*nd->nX[l]+UY*nd->nY[l]))*f;

			UX = B[i].shape->potential(solver.hh[i][0][m],   1);
			UY = B[i].shape->potential(solver.hh[i][0][m+1], 1);
			
			K[6]  -= G0*(UX*nd->nX[l]*2.+alpha*(UX*nd->nX[l]+UY*nd->nY[l]))*f;
			K[7]  -= G0*(UX*nd->nY[l]+UY*nd->nX[l])*f;
			K[8]  -= G0*(UY*nd->nY[l]*2.+alpha*(UX*nd->nX[l]+UY*nd->nY[l]))*f;
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////
//...auxiliary function for parameters of doubly connected media;
void CCohes2D::set_fasa_hmg(double nju1, double nju2, double G1, double G2, double C1, double C2)
{
	if (size_of_param() > NUM_SHEAR+3+NUM_SHIFT) {
		param[NUM_SHEAR-1] = C1;
		param[NUM_SHEAR-1+NUM_SHIFT] = C2;
		param[NUM_SHEAR] = G1;
		param[NUM_SHEAR+NUM_SHIFT] = G2;
		param[NUM_SHEAR+1] = nju1;
		param[NUM_SHEAR+1+NUM_SHIFT] = nju2;
		param[NUM_SHEAR+2] = sqrt((.5-nju1)/(1.-nju1)*C1/G1);
		param[NUM_SHEAR+2+NUM_SHIFT] = sqrt((.5-nju2)/(1.-nju2)*C2/G2);
		param[NUM_SHEAR+3] = sqrt(C1/G1);
		param[NUM_SHEAR+3+NUM_SHIFT] = sqrt(C2/G2);
	}
}

void CCohes2D::set_fasa_hmg(double nju1, double nju2, double nju3, double G1, double G2, double G3, double C1, double C2, double C3)
{
	if (size_of_param() > NUM_SHEAR+3+NUM_SHIFT*2) {
		param[NUM_SHEAR-1] = C1;
		param[NUM_SHEAR-1+NUM_SHIFT] = C2;
		param[NUM_SHEAR-1+NUM_SHIFT*2] = C3;
		param[NUM_SHEAR] = G1;
		param[NUM_SHEAR+NUM_SHIFT] = G2;
		param[NUM_SHEAR+NUM_SHIFT*2] = G3;
		param[NUM_SHEAR+1] = nju1;
		param[NUM_SHEAR+1+NUM_SHIFT] = nju2;
		param[NUM_SHEAR+1+NUM_SHIFT*2] = nju3;
		param[NUM_SHEAR+2] = sqrt((.5-nju1)/(1.-nju1)*C1/G1);
		param[NUM_SHEAR+2+NUM_SHIFT] = sqrt((.5-nju2)/(1.-nju2)*C2/G2);
		param[NUM_SHEAR+2+NUM_SHIFT*2] = sqrt((.5-nju3)/(1.-nju3)*C3/G3);
		param[NUM_SHEAR+3] = sqrt(C1/G1);
		param[NUM_SHEAR+3+NUM_SHIFT] = sqrt(C2/G2);
		param[NUM_SHEAR+3+NUM_SHIFT*2] = sqrt(C3/G3);
	}
}

//////////////////////////////////////////////////////
//...counting kernel for solving double plane problem;
Num_State CCohes2D::computing_header(Num_Comput Num)
{
	int i, j, k, elem, id_dir, n_rhs = 2;
	char msg[201];

	if (! solver.mode(NO_MESSAGE)) {
		Message(" ");
		sprintf(msg, "CCohes2D sample: N_sm = %d, N_mpl = %d, C1 = %g, C2 = %g", N,
				  UnPackInts(get_param(NUM_MPLS)), get_param(NUM_SHEAR-1), get_param(NUM_SHEAR-1+NUM_SHIFT));
		Message(msg);

		Message(" ");
		Message("Junction counting...");

		switch (Num){
			case   BASIC_COMPUT: Message("Analytical Blocks..."); break;
			case MAPPING_COMPUT: Message("FEM Blocks...");			break;
		}
		Message(" ");
	}

//////////////////////////////////////////////
//...устанавливаем мультиполи межфазного слоя;
	if (NUM_PHASE > 0)
	for (k = 0; k < N; k++) if (B[k].link[0] >= NUM_PHASE)
	for (j = B[k].link[0]; j > 0; j--) 
	if ((i = B[k].link[j]) >= 0 && B[k].link[NUM_PHASE] != B[i].link[NUM_PHASE])
	if (/*B[k].link[NUM_PHASE] == -1 &&*/ B[k].mp) { 
		int m, m1, m2;
		m  = stru.geom[(i = stru.geom_ptr[k])+1]-2;
		m1 = stru.geom[i+4+(j-1)%m];
		m2 = stru.geom[i+4+j%m];
		double X0 = stru.X[m1], Y0 = stru.Y[m1], 
				 X1 = stru.X[m2], Y1 = stru.Y[m2];
		B[k].type = ELLI_BLOCK;
		B[k].mp[4] = arg0(comp(X1-X0, Y1-Y0))-M_PI_2;
		B[k].mp[5] = 0.;
		B[k].mp[6] = 0.;
	}

//////////////////////////////////
//...определяем блочную структуру;
	solver.set_blocks(N, n_rhs); //<==== number of saved potentials!!!
	solver.n += 12;//<==== number of additional auxilliary arrays!!!
	for (k = 0; k < solver.N;  k++)
		solver.set_links(k, B[k].link);

	shapes_init(INITIAL_STATE);
	shapes_init(NULL_STATE);

////////////////////////////////////////////////////////////////////
//...добавляем блоки на периодических связях и на границе включений;
	if (solv%ENERGY_SOLVING == PERIODIC_SOLVING) { 
		double par[6]; SetGeomBounding(par); d_parm[0] = par[1]-par[0]; d_parm[1] = par[3]-par[2];
		for (k = 0; k < N; k++) SkeletonBounding(B[k], par);
		for (k = 0; k < N; k++) if (B[k].link && B[k].link[NUM_PHASE] == -1) {
			i = 0; while ((elem =  geom_plink_2D(B[k], i, id_dir, par)) >= 0)			  
			solver.add_link(k, elem);
			i = j = 0; while ((elem = block_plink_2D(B[k], i, j, id_dir, par)) >= 0)			  
			solver.add_link(k, elem);
		}
	}
	LinkPhase2D(MAX_PHASE);

/////////////////////////////////////////////////////////
//...делаем перенумерацию структуры и задаем размерность;
	if (! solver.struct_permutat(solver.id_change == EXTERN_STATE ? NULL_STATE : /*NULL*/OK_STATE) || 
		 ! solver.inverse_index()) {
		return ERR_STATE;
	}
	for (k = 0; k < solver.N; k++)
		solver.set_dimension(k, freedom_block(k));

	solver.struct_init();

	if (solver.mode(FULLY_MODE)) { 
		solver.test_struct("1", 0);
		solver.test_struct("2", 1);
	}
	return OK_STATE;
}

//////////////////////////////////////////////////////////////////////////
//...calculation of function values (in common coordinate system) on grid;
void CCohes2D::GetFuncAllValues(double X, double Y, double Z, double * F, int i, Num_Value id_F, int id_variant, int iparam)
{
	if (! F) return;
	double P[6]  = { X, Y, Z, 1., 0., 0.};

/////////////////////////////////////
//...operation with all input points;
	if ( 0 <= i && i < N && B[i].shape && B[i].mp && B[i].link[0] >= NUM_PHASE) {
		int shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, m = solver.id_norm;

//////////////////////////////////////////////////////
//...reset auxilliary arrays and calculation function;
		for (int num = m; num < solver.n; num++)
			memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));

		B[i].shape->make_local(P);
		B[i].shape->norm_local(P+3);
		switch (id_F) {
			case DISPL_VALUE: {
//////////////////////////
//...common displacements;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				jump1_common_x(P, i, 0); 
				jump1_common_y(P, i, 1); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				jump1_compress_x(P, i, 0);
				jump1_compress_y(P, i, 1);
				//jump_make_common(i, 0); //...вместо B[i].shape->norm_common2D(F);

/////////////////////////////////////
//...calculation common displacement;
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common2D(F);

				//if ((solv%ENERGY_SOLVING) && 1) F[(solv%ENERGY_SOLVING)-1] -=X;
				F[1] = F[0]-X;
			}  break;
			case DISPL_GRAD_VALUE: {
///////////////////////////////////////////////////////////////
//...common displacements with independent cohesion potentials;
				B[i].shape->parametrization_hess(P, 1);

				jump1_common_x(P, i, 0); 
				jump1_common_y(P, i, 1); 

/////////////////////////////////////
//...calculation common displacement;
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common2D(F);
			}  break;
			case DISPL_CLASSIC_VALUE: {
//////////////////////////////////
//...displacement (classic field);
				B[i].shape->parametrization_hess(P, 1);

				jump1_classic_x(P, i, 0); 
				jump1_classic_y(P, i, 1); 

/////////////////////////////////////
//...calculation common displacement;
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common2D(F);
			}  break;
			case NORMAL_R_CLASSIC_VALUE: {
				double RR = sqrt(sqr(P[0])+sqr(P[1])+sqr(P[2]));
				if (RR != 0.) {
					P[3] = P[0]/RR; 
					P[4] = P[1]/RR;
					P[5] = P[2]/RR;
				}
				else {
					P[3] = 1.; P[4] = P[5] = 0.;
				}
				B[i].shape->norm_local(P+3);
//////////////////////////////////////////
//...classic derivatives of displacements;
				B[i].shape->parametrization_hess(P, 1);

				jump2_classic_x (P, i, 0); 
				jump2_classic_y (P, i, 1); 

/////////////////////////////////////////////////
//...calculation classic derivatives (ux and uy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common2D(F);
			}  break;
			case NORMAL_X_CLASSIC_VALUE: {
				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
//////////////////////////////////////////
//...classic derivatives of displacements;
				B[i].shape->parametrization_hess(P, 1);

				jump2_classic_x (P, i, 0); 
				jump2_classic_y (P, i, 1); 

/////////////////////////////////////////////////
//...calculation classic derivatives (ux and uy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common2D(F);
			}  break;
			case NORMAL_Y_CLASSIC_VALUE: {
				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
//////////////////////////////////////////
//...classic derivatives of displacements;
				B[i].shape->parametrization_hess(P, 1);

				jump2_classic_x (P, i, 0); 
				jump2_classic_y (P, i, 1); 

/////////////////////////////////////////////////
//...calculation classic derivatives (ux and uy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common2D(F);
			}  break;
			case DISPL_COHESION_VALUE: {
///////////////////////////////////
//...displacement (cohesion field);
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				jump1_common_x (P, i, 0); jump_admittance(0, i, 0, 0.);
				jump1_common_y (P, i, 1); jump_admittance(0, i, 1, 0.);

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				jump1_compress_x(P, i, 0);
				jump1_compress_y(P, i, 1);

///////////////////////////////////////
//...calculation cohesion displacement;
				F[0] = -B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = -B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common2D(F);
			}  break;
			case DISPL_GRAD_COHESION_VALUE: {
////////////////////////////////////////////////////////
//...cohesion displacements with independent potentials;
				B[i].shape->parametrization_hess(P, 1);

				jump1_cohesion_x(P, i, 0); 
				jump1_cohesion_y(P, i, 1); 

///////////////////////////////////////
//...calculation cohesion displacement;
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common2D(F);
			}  break;
			case STRESS_X_VALUE: {
				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
/////////////////////
//...common stresses;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, i);
				jump4_common_x (P, i, 0); 
				jump4_common_y (P, i, 1); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump4_compress_x(P, i, 0);
				jump4_compress_y(P, i, 1);

////////////////////////////////////////////////////
//...calculation common stress tensor (txx and txy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common2D(F);
			}  break;
			case NORMAL_X_GRAD_COHESION_VALUE: {
				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
/////////////////////////////////////////
//...common derivatives of displacements;
				B[i].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, i); hessian_deriv_N(4, P, i);
				jump2_cohesion_x (P, i, 0); 
				jump2_cohesion_y (P, i, 1); 

////////////////////////////////////////////////////
//...calculation common stress tensor (txx and txy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common2D(F);
			}  break;
			case NORMAL_Y_GRAD_COHESION_VALUE: {
				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
/////////////////////////////////////////
//...common derivatives of displacements;
				B[i].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, i); hessian_deriv_N(4, P, i);
				jump2_cohesion_x (P, i, 0); 
				jump2_cohesion_y (P, i, 1); 

////////////////////////////////////////////////////
//...calculation common stress tensor (txx and txy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common2D(F);
			}  break;
			case STRESS_Y_VALUE: {
				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
/////////////////////
//...common stresses;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, i);
				jump4_common_x (P, i, 0); 
				jump4_common_y (P, i, 1); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump4_compress_x(P, i, 0);
				jump4_compress_y(P, i, 1);

////////////////////////////////////////////////////
//...calculation common stress tensor (txy and tyy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common2D(F);
			}  break;
			case STRESS_X_GRAD_VALUE: {
				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
/////////////////////
//...common stresses;
				B[i].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, i); hessian_deriv_N(4, P, i);
				jump4_common_x (P, i, 0); 
				jump4_common_y (P, i, 1); 

////////////////////////////////////////////////////
//...calculation common stress tensor (txx and txy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common2D(F);
			}  break;
			case STRESS_Y_GRAD_VALUE: {
				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
/////////////////////
//...common stresses;
				B[i].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, i); hessian_deriv_N(4, P, i);
				jump4_common_x (P, i, 0); 
				jump4_common_y (P, i, 1); 

////////////////////////////////////////////////////
//...calculation common stress tensor (txy and tyy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common2D(F);
			}  break;
			case STRESS_X_CLASSIC_VALUE: {
				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
/////////////////////
//...common stresses;
				B[i].shape->parametrization_hess(P, 1);

				jump4_classic_x (P, i, 0); 
				jump4_classic_y (P, i, 1); 

////////////////////////////////////////////////////
//...calculation common stress tensor (txx and txy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common2D(F);
			}  break;
			case STRESS_Y_CLASSIC_VALUE: {
				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
/////////////////////
//...common stresses;
				B[i].shape->parametrization_hess(P, 1);

				jump4_classic_x (P, i, 0); 
				jump4_classic_y (P, i, 1); 

////////////////////////////////////////////////////
//...calculation common stress tensor (txy and tyy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common2D(F);
			}  break;
			case STRESS_X_COHESION_VALUE: {
				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
///////////////////////////////////////
//...surface stresses (cohesion field);
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N(1, P, i);
				jump4_common_x (P, i, 0); 
				jump4_common_y (P, i, 1); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump4_compress_x(P, i, 0); jump_admittance(0, i, 0, 0.);
				jump4_compress_y(P, i, 1); jump_admittance(0, i, 1, 0.);

/////////////////////////////////////////////
//...calculation stress tensor (txx and txy);
				F[0] = -B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = -B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common2D(F);
			}  break;
			case STRESS_Y_COHESION_VALUE: {
				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
///////////////////////////////////////
//...surface stresses (cohesion field);
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N(1, P, i);
				jump4_common_x (P, i, 0); 
				jump4_common_y (P, i, 1); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump4_compress_x(P, i, 0); jump_admittance(0, i, 0, 0.);
				jump4_compress_y(P, i, 1); jump_admittance(0, i, 1, 0.);

/////////////////////////////////////////////
//...calculation stress tensor (txy and tyy);
				F[0] = -B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = -B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common2D(F);
			}  break;
			case NORMAL_X_VALUE: {
				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
/////////////////////////////////////////
//...common derivatives of displacements;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, i);
				jump2_common_x (P, i, 0); 
				jump2_common_y (P, i, 1); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump2_compress_x(P, i, 0);
				jump2_compress_y(P, i, 1);

////////////////////////////////////////////////////
//...calculation common stress tensor (txx and txy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common2D(F);
			}  break;
			case NORMAL_Y_VALUE: {
				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
/////////////////////////////////////////
//...common derivatives of displacements;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, i);
				jump2_common_x (P, i, 0); 
				jump2_common_y (P, i, 1); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump2_compress_x(P, i, 0);
				jump2_compress_y(P, i, 1);

////////////////////////////////////////////////////
//...calculation common stress tensor (txx and txy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common2D(F);
			}  break;
			case NORMAL_R_GRAD_VALUE: {
				double RR = sqrt(sqr(P[0])+sqr(P[1])+sqr(P[2]));
				if (RR != 0.) {
					P[3] = P[0]/RR; 
					P[4] = P[1]/RR;
					P[5] = P[2]/RR;
				}
				else {
					P[3] = 1.; P[4] = P[5] = 0.;
				}
				B[i].shape->norm_local(P+3);
/////////////////////////////////////////
//...common derivatives of displacements;
				B[i].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, i); hessian_deriv_N(4, P, i);
				jump2_common_x (P, i, 0); 
				jump2_common_y (P, i, 1); 

////////////////////////////////////////////////////
//...calculation common stress tensor (txx and txy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common2D(F);
			}  break;
			case NORMAL_X_GRAD_VALUE: {
				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
/////////////////////////////////////////
//...common derivatives of displacements;
				B[i].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, i); hessian_deriv_N(4, P, i);
				jump2_common_x (P, i, 0); 
				jump2_common_y (P, i, 1); 

////////////////////////////////////////////////////
//...calculation common stress tensor (txx and txy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common2D(F);
			}  break;
			case NORMAL_Y_GRAD_VALUE: {
				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
/////////////////////////////////////////
//...common derivatives of displacements;
				B[i].shape->parametrization_hess(P, 1);

				hessian_deriv_N(1, P, i); hessian_deriv_N(4, P, i);
				jump2_common_x (P, i, 0); 
				jump2_common_y (P, i, 1); 

////////////////////////////////////////////////////
//...calculation common stress tensor (txx and txy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common2D(F);
			}  break;
			case DILAT_VALUE: {//...divergency of displacements;
				double G = (.5-get_param(NUM_SHEAR+1+shift))/((1.-get_param(NUM_SHEAR+1+shift))*get_param(NUM_SHEAR+shift));
///////////////////////////////////////
//...classical and pressure potentials;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
 				B[i].shape->parametrization_grad(0, P, 1);
				B[i].shape->parametrization_grad(1, P, 2);

				B[i].shape->adm_x(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 0), 1.); 
				B[i].shape->adm_y(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), 1.); 
				B[i].shape->adm_x(1, B[i].shape->FULL(solver.hh[i][0][m], 1, 0), -1.); 
				B[i].shape->adm_y(1, B[i].shape->FULL(solver.hh[i][0][m], 1, 1), -1.); 

/////////////////////////////////
//...divergency of displacements;
				F[0] = B[i].shape->potential(solver.hh[i][0][m], id_variant)*G;
			}	break;
			case ROTATION_VALUE: {//...rotation of displacements;
				double G = .5/get_param(NUM_SHEAR+shift);
////////////////////////////////////
//...classical and shear potentials;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
 				B[i].shape->parametrization_grad(0, P, 1);
				B[i].shape->parametrization_grad(1, P, 2);

				B[i].shape->adm_y(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 0), -1.); 
				B[i].shape->adm_x(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1), 1.); 
				B[i].shape->adm_y(1, B[i].shape->FULL(solver.hh[i][0][m], 1, 0), 1.); 
				B[i].shape->adm_x(1, B[i].shape->FULL(solver.hh[i][0][m], 1, 1), -1.); 

///////////////////////////////
//...rotation of displacements;
				F[0] = B[i].shape->potential(solver.hh[i][0][m], id_variant)*G;
			}	break;
			case ENERGY_VALUE: {//...local energy in cell;
				double e11, e12, e22, theta, u1, u2, G1 = get_param(NUM_SHEAR+shift), 
						 G2 = get_param(NUM_SHEAR+1+shift)*G1/(.5-get_param(NUM_SHEAR+1+shift)), C1 = get_param(NUM_SHEAR-1+shift);
/////////////////////////////////////////
//...common derivatives of displacements;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);

				hessian_deriv_N(1, P, i);
				jump2_common_x (P, i, 0); 
				jump2_common_y (P, i, 1); 

				jump1_common_x (P, i, 4); jump_admittance(0, i, 4, 0.);
				jump1_common_y (P, i, 5); jump_admittance(0, i, 5, 0.);

				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);

				hessian_deriv_N(1, P, i);
				jump2_common_x (P, i, 2); 
				jump2_common_y (P, i, 3); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump2_compress_x(P, i, 2);
				jump2_compress_y(P, i, 3);

				jump1_compress_x(P, i, 4); 
				jump1_compress_y(P, i, 5); 

				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);

				hessian_deriv_N (1, P, i);
				jump2_compress_x(P, i, 0);
				jump2_compress_y(P, i, 1);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//...вычисляем дилатацию, компоненты тензора деформации и когезионное поле (с учетом поворота компонент);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common2D(F);

				e11 =	F[0];
				e12 =	F[1];

				F[0] = B[i].shape->potential(solver.hh[i][0][m+2], id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+3], id_variant);
				B[i].shape->norm_common2D(F);

				e22 = F[1];
				e12 =(F[0]+e12)*.5;
				theta = e11+e22;

				u1 = B[i].shape->potential(solver.hh[i][0][m+4], id_variant);
				u2 = B[i].shape->potential(solver.hh[i][0][m+5], id_variant);

//////////////////
//...local energy;
				F[0] = G1*(e11*e11+2.*e12*e12+e22*e22)+(G2*theta*theta+C1*(u1*u1+u2*u2))*.5;
			}	break;
			case STRAIN_X_VALUE: {
/////////////////////////////////////////
//...common derivatives of displacements;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);

				hessian_deriv_N(1, P, i);
				jump2_common_x (P, i, 0); 
				jump2_common_y (P, i, 1); 

				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);

				hessian_deriv_N(1, P, i);
				jump2_common_x (P, i, 2); 
				jump2_common_y (P, i, 3); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump2_compress_x(P, i, 2);
				jump2_compress_y(P, i, 3);

				jump1_compress_x(P, i, 4); 
				jump1_compress_y(P, i, 5); 

				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);

				hessian_deriv_N (1, P, i);
				jump2_compress_x(P, i, 0);
				jump2_compress_y(P, i, 1);

////////////////////////////////////////////////////////
//...calculation common dilatation tensor (exx and exy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m+2], id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+3], id_variant);
				B[i].shape->norm_common2D(F); double exy = F[0];

				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common2D(F);

				F[1] = (F[1]+exy)*.5;
			}  break;
			case STRAIN_Y_VALUE: {
/////////////////////////////////////////
//...common derivatives of displacements;
				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+3+shift)); 
				B[i].shape->parametrization_hess(P, 1);

				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);

				hessian_deriv_N(1, P, i);
				jump2_common_x (P, i, 0); 
				jump2_common_y (P, i, 1); 

				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);

				hessian_deriv_N(1, P, i);
				jump2_common_x (P, i, 2); 
				jump2_common_y (P, i, 3); 

				B[i].shape->set_shape(1, B[i].shape->get_R(1), get_param(NUM_SHEAR+2+shift)); 
				B[i].shape->parametrization_hess(1, P, 1);

				hessian_deriv_N (1, P, i);
				jump2_compress_x(P, i, 2);
				jump2_compress_y(P, i, 3);

				jump1_compress_x(P, i, 4); 
				jump1_compress_y(P, i, 5); 

				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);

				hessian_deriv_N (1, P, i);
				jump2_compress_x(P, i, 0);
				jump2_compress_y(P, i, 1);

////////////////////////////////////////////////////////
//...calculation common dilatation tensor (exx and exy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common2D(F); double exy = F[1];

				F[0] = B[i].shape->potential(solver.hh[i][0][m+2], id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+3], id_variant);
				B[i].shape->norm_common2D(F);

				F[0] = (F[0]+exy)*.5;
			}  break;
			default : F[0] = i; F[1] = 0.;
     }
  }
}

///////////////////////////////////////////
//...calculaion energy by block functional;
void CCohes2D::GetEnergyValue(int k, double * energy)
{
   if (solver.dim && solver.dim[k] && solver.JR && solver.JR[k] && solver.hh && solver.hh[k][0].GetMatrix() && 
							solver.TL && solver.TL[k]) {
		double ** TL = solver.TL[k][solver.JR[k][solver.JR_DIAG]-solver.JR_SHIFT].GetMatrix(), * h = solver.hh[k][0][solver.id_norm];
		int i, j, NN = solver.dim[k];

///////////////////////
//...вычисляем энергию;
		block_shape_init(B[k], NO_STATE);
		if (h && TL) { 
			double measure = .5/get_param(NUM_SHEAR);
			for (i = 0; i < NN; i++)
			for (j = 0; j < NN; j++)
				energy[-B[k].link[NUM_PHASE]-1] += TL[i][j]*h[i]*h[j]*measure;
		}
	}
}

/////////////////////////////////////////////////////////////////////
//...одномерная модель межфазного слоя (симметричная слоистая среда);
void CCohes2D::TakeLayerModel(double L, double H, double l, double nju1, double nju2, double G1, double G2, double l1, double l2, double A)
{
	double P1[3] = {-L, -H*.5, 0.},
			 P2[3] = {-l, -H*.5, 0.},
			 P3[3] = {-l,  H*.5, 0.},
			 P4[3] = {-L,  H*.5, 0.},
			 P5[3] = {-l, -H*.5, 0.},
			 P6[3] = { l, -H*.5, 0.},
			 P7[3] = { l,  H*.5, 0.},
			 P8[3] = {-l,  H*.5, 0.};

/////////////////////////////////////
//...образуем образец из трех блоков;
	CCells * ce; 
   init_blocks(3);

	B[0].type = ELLI_BLOCK;
	B[0].bar	 = new CCells;
	ce = new CCells;
	ce->get_quad_facet(P1, P2, P3, P4, ZERO_STATE);
	B[0].bar->bar_add(ce);
	set_block3D(B[0], SPHERE_GENUS, L-l, 0, -L);
	set_link (B[0], NUM_PHASE+1);
	B[0].link[2] = -1+SRF_STATE;
	B[0].link[NUM_PHASE+1] = 1;
	B[0].mp[7] = 0.; //...нормируем центр на краю ячейки;

	B[1].type = ELLI_BLOCK;
	B[1].bar	 = new CCells;
	ce = new CCells;
	ce->get_quad_facet(P5, P6, P7, P8, ZERO_STATE);
	B[1].bar->bar_add(ce);
	set_block3D(B[1], SPHERE_GENUS, l);
	set_link (B[1], NUM_PHASE+2);
	B[1].link[2] = -2+SRF_STATE;
	B[1].link[4] = SRF_STATE;
	B[1].link[NUM_PHASE]	 = -2;
	B[1].link[NUM_PHASE+1] = 2;
	B[1].link[NUM_PHASE+2] = 0;
	B[1].mp[7] = 0.; //...нормируем центр в центре ячейки;

	B[2].type = ELLI_BLOCK;
	B[2].bar	 = new CCells;
	ce = new CCells; P1[0] = -P1[0]; P2[0] = -P2[0]; P3[0] = -P3[0]; P4[0] = -P4[0];
	ce->get_quad_facet(P2, P1, P4, P3, ZERO_STATE);
	B[2].bar->bar_add(ce);
	set_block3D(B[2], SPHERE_GENUS, L-l, 0, L);
	set_link (B[2], NUM_PHASE+1);
	B[2].link[4] = -1+SRF_STATE;
	B[2].link[NUM_PHASE+1] = 1;
	B[2].mp[7] = 0.; //...нормируем центр на краю ячейки;

////////////////////////////////////
//...устанавливаем параметры задачи;
	double C1, C2, k1 = (1.-nju1)/(.5-nju1)*G1, k2 = (1.-nju2)/(.5-nju2)*G2;
	set_mpls(PackInts(1, 1)); //...multipoles degree;
	set_quad(PackInts(4, 2)); //...quadrature degree;
	set_normaliz(1.);
	set_param(NUM_ADHES, A);
	set_lagrange(1e5);		  //...Lagrange corfficient for LSM;
	set_fasa_hmg(nju1, nju2, G1, G2, C1 = l1 ? G1/sqr(l1) : 1e33, C2 = l2 ? G2/sqr(l2) : 1e33);

/////////////////////////////////////////////////
//...инициализируем мультиполи и блочную матрицу;
	solver.set_blocks(N, 1); //<==== number of saved potentials !!!
	solver.n += 10;//<==== number of additional auxilliary arrays!!!
   shapes_init(INITIAL_STATE);
	B[0].shape->degree_init1(1, 0, solver.id_norm, 2);//...в слое задаем степень поменьше;
	B[1].shape->degree_init1(1, 0, solver.id_norm, 2);
	B[2].shape->degree_init1(1, 0, solver.id_norm, 2);
	shapes_init(NULL_STATE);
	int k;
	for (k = 0; k < solver.N;  k++)
	solver.set_dimension(k, freedom_block(k));
   solver.struct_init();

   double par[6] = {MAX_HIT, MIN_HIT, MAX_HIT, MIN_HIT, MAX_HIT, MIN_HIT};
	for (k = 0; k < N; k++)
		SkeletonBounding(B[k], par);
	d_parm[0] = par[1] - par[0];
	d_parm[1] = par[3] - par[2];

////////////////////////////////////////////////
//...вводим коэффициенты аналитического решения;
	double kk1 = sqrt(C1/k1), t1 = exp(-kk1*(L-l)), tt1 = sqr(t1), al1 = kk1*(1.+tt1)/(1.-tt1),
			 kk2 = sqrt(C2/k2), t2 = exp(-kk2*l),		tt2 = sqr(t2), al2 = kk2*(1.+tt2)/(1.-tt2),
			 AA = k2*al1+k1*al2-A*al1*al2,
			 BB = k2-k1-A*al2,
			 CC = k2-k1+A*al1,
			 QQ = AA*(k2*(L-l)+k1*l+A)-BB*CC,
			 HM = k2*L*AA/QQ,
			 HD = k1*L*AA/QQ,
			 DM = k2*L*BB/QQ*2.*t1/(1.-tt1),
			 DD = k1*L*CC/QQ*2.*t2/(tt2-1.);

////////////////////////////////////////////////////////////
//...заносим коэффициенты в представление Папковича-Нейбера;
	B[0].shape->set_R(1.);
	B[0].shape->A[0][0] = -L*4.*G1*(1.-nju1)/(3.-4.*nju1);
	B[0].shape->A[0][2] = HM*2.*G1*(1.-nju1)/(1.-2.*nju1);
	B[0].shape->A[0][7] = DM*C1/kk1;
	B[1].shape->set_R(1.);
	B[1].shape->A[0][2] = HD*2.*G2*(1.-nju2)/(1.-2.*nju2);
	B[1].shape->A[0][7] = DD*C2/kk2;
	B[2].shape->set_R(1.);
	B[2].shape->A[0][0] =  L*4.*G1*(1.-nju1)/(3.-4.*nju1);
	B[2].shape->A[0][2] = HM*2.*G1*(1.-nju1)/(1.-2.*nju1);
	B[2].shape->A[0][7] = DM*C1/kk1;
}

//////////////////////////////////////////////////////////////////////////
//...эффективный модуль Юнга на основе модели симметричной слоистой среды;
double CCohes2D::TakeLayer_E1(double ff)
{
	double nj1 = get_param(NUM_SHEAR+1), 
			 nj2 = get_param(NUM_SHEAR+1+NUM_SHIFT),
			 G1 = get_param(NUM_SHEAR),
			 G2 = get_param(NUM_SHEAR+NUM_SHIFT),
			 k1 = (1.-nj1)/(.5-nj1)*G1,
			 k2 = (1.-nj2)/(.5-nj2)*G2,
			 C1 = get_param(NUM_SHEAR-1),
			 C2 = get_param(NUM_SHEAR-1+NUM_SHIFT),
			 A  = get_param(NUM_ADHES),
			 B  = get_param(NUM_ADHES+1);
////////////////////////////////////////////////////////////
//...вывод эффективного модуля одноосного растяжения/сжатия;
	double kk1 = sqrt(C1/k1), t1 = ff ? exp(-kk1*(1./ff-1.)) : 0., tt1 = sqr(t1), al1 = kk1*(1.+tt1)/(1.-tt1),
			 kk2 = sqrt(C2/k2), t2 = exp(-kk2),	tt2 = sqr(t2), al2 = kk2*(1.+tt2)/(1.-tt2),
			 AA = k2*al1+k1*al2-A*al1*al2,
			 BB = k2-k1-A*al2,
			 CC = k2-k1+A*al1,
			 QQ = AA*(k2*(1.-ff)+k1*ff+A*ff)-BB*CC*ff,
			 kk_eff = k1*k2*AA/QQ, GG_eff;
////////////////////////////////////////////////////
//...повторяем вывод для эффективного модуля сдвига;
			 kk1 = sqrt(C1/G1); t1 = ff ? exp(-kk1*(1./ff-1.)) : 0.; tt1 = sqr(t1); al1 = kk1*(1.+tt1)/(1.-tt1);
			 kk2 = sqrt(C2/G2); t2 = exp(-kk2);	tt2 = sqr(t2); al2 = kk2*(1.+tt2)/(1.-tt2);
			 AA = G2*al1+G1*al2-B*al1*al2;
			 BB = G2-G1-B*al2;
			 CC = G2-G1+B*al1;
			 QQ = AA*(G2*(1.-ff)+G1*ff+B*ff)-BB*CC*ff;
			 GG_eff = G1*G2*AA/QQ;
	return(3.*kk_eff-4.*GG_eff)*GG_eff/(kk_eff-GG_eff);
}


////////////////////////////////////////////////////////////////////////////
//...эффективный модуль сдвига на основе модели симметричной слоистой среды;
double CCohes2D::TakeLayer_G1(double ff)
{
	double nj1 = get_param(NUM_SHEAR+1), 
			 nj2 = get_param(NUM_SHEAR+1+NUM_SHIFT),
			 G1 = get_param(NUM_SHEAR),
			 G2 = get_param(NUM_SHEAR+NUM_SHIFT),
			 C1 = get_param(NUM_SHEAR-1),
			 C2 = get_param(NUM_SHEAR-1+NUM_SHIFT),
			 A  = get_param(NUM_ADHES),
			 B  = get_param(NUM_ADHES+1);
////////////////////////////////////////////////////
//...повторяем вывод для эффективного модуля сдвига;
	double kk1 = sqrt(C1/G1), t1 = ff ? exp(-kk1*(1./ff-1.)) : 0., tt1 = sqr(t1), al1 = kk1*(1.+tt1)/(1.-tt1),
			 kk2 = sqrt(C2/G2), t2 = exp(-kk2),	tt2 = sqr(t2), al2 = kk2*(1.+tt2)/(1.-tt2),
			 AA = G2*al1+G1*al2-B*al1*al2,
			 BB = G2-G1-B*al2,
			 CC = G2-G1+B*al1,
			 QQ = AA*(G2*(1.-ff)+G1*ff+B*ff)-BB*CC*ff,
			 GG_eff = G1*G2*AA/QQ;
	return(GG_eff);
}

/////////////////////////////////////////////////
//...реализация регулярных масштабных множителей;
void scale_regul(double R, double kappa, double hh[3], double eps = 1e-14, int k_limit = 200)
{
	int k = 0;
	double dd = sqr(.5*R*kappa), f = 1., g = 1., d0 = sqr(kappa)*.5; 
	for (int k, m = 0; m < 2; m++, g *= (d0/m)) { 
		k = 0; hh[m] = (f = g);
		do {
			k += 1; hh[m] += (f *= dd/(k*(k+m)));
		} 
		while (fabs(f) > EE);
	}
	hh[2] = (sqr(kappa)*hh[0]-2.*hh[1]);
}

/////////////////////////////////////////////////
//...реализация сингулярных масштабных множителей;
void scale_irreg(double R, double kappa, double hh[3], double eps = 1e-14, int k_limit = 200)
{
	double dd = sqr(.5*R*kappa), dd_inv = 1./dd, g_euler = 0.577215664901532860606512, f = 1., g = 1., 
			 d0 = sqr(kappa)*.5, h_m = 0., h;
	for ( int k, m = 0; m < 2; m++, h_m += 1./m, g *= (d0/m)) { 
		for (hh[m]  = (f = -.5*dd_inv*m*g), k = 1; k < m; k++) 
			  hh[m] += (f *= -k*(m-k)*dd_inv);

		k = 0;  hh[m] += (h = -g_euler+(h_m-log(dd))*.5)*(f = g);
		do {
			k += 1; hh[m] += (h += (k+m*.5)/(k*(k+m)))*(f *= dd/(k*(k+m)));
		} 
		while (f*(1.+fabs(h)) > EE);
	}
	hh[2] = (sqr(kappa)*hh[0]-2.*hh[1]);
}

//////////////////////////////////////////////////////////////////////////
//...четырехфазная модель для цилиндрического включения (прямой алгоритм);
double CCohes2D::TakeEshelby_volm(double ff, double ff_l)
{
	double mu_I = get_param(NUM_SHEAR+NUM_SHIFT), nu_I = get_param(NUM_SHEAR+1+NUM_SHIFT), C_I = get_param(NUM_SHEAR-1+NUM_SHIFT),
			 mu_M = get_param(NUM_SHEAR), nu_M = get_param(NUM_SHEAR+1), C_M = get_param(NUM_SHEAR-1),
			 mu_L = get_param(NUM_SHEAR+NUM_SHIFT*2), nu_L = get_param(NUM_SHEAR+1+NUM_SHIFT*2), C_L = get_param(NUM_SHEAR-1+NUM_SHIFT*2),
			 K_I  = mu_I/(1.-2.*nu_I), K_M = mu_M/(1.-2.*nu_M), K_L = mu_L/(1.-2.*nu_L), c0 = ff+ff_l, c1 = ff/c0, DD,
			 ku_I = (1.-nu_I)/(.5-nu_I)*mu_I, ku_M = (1.-nu_M)/(.5-nu_M)*mu_M, ku_L = (1.-nu_L)/(.5-nu_L)*mu_L, AA = get_param(NUM_ADHES);
	if (C_I == 0. || C_M == 0. || C_L == 0.) { //...классическая четырехфазная модель;
		DD = ((K_I-K_M)-(1.-c1)*(K_I-K_L)*(1.+(K_M-K_L)/ku_L))/(1.+(1.-c1)*(K_I-K_L)/ku_L);
		return (K_M+c0*DD/(1.+(1.-c0)*DD/ku_M));
	}
	double kk_I = sqrt(C_I/ku_I), kk_M = sqrt(C_M/ku_M), kk_L = sqrt(C_L/ku_L), RR1 = 1./sqrt(c1), RR2 = 1./sqrt(ff), hh[3];
	scale_regul(1., kk_I, hh);
	double HH1 = hh[1], HH2 = hh[2];
	scale_regul(1., kk_L, hh);
	double HHL1 = hh[1], HHL2 = hh[2];
	scale_irreg(1., kk_L, hh);
	double JJL1 = hh[1], JJL2 = hh[2], HHM1, JJM1, HHM2, JJM2;
	scale_regul(RR1, kk_L, hh); HHM1 = hh[1]; HHM2 = hh[2];
	scale_irreg(RR1, kk_L, hh); JJM1 = hh[1]; JJM2 = hh[2];
	double HHK1, JJK1, HHK2, JJK2;
	scale_regul(RR1, kk_M, hh); HHK1 = hh[1]; HHK2 = hh[2];
	scale_irreg(RR1, kk_M, hh); JJK1 = hh[1]; JJK2 = hh[2];
	double HHP1, JJP1, HHP2, JJP2;
	scale_regul(RR2, kk_M, hh); HHP1 = hh[1]; HHP2 = hh[2];
	scale_irreg(RR2, kk_M, hh); JJP1 = hh[1]; JJP2 = hh[2];
	double matr[11][12] = {
		{ 1.,		  -HH1, -1., -1.,		  HHL1,		 JJL1, 0., 0., 0., 0., 0., 0. }, //...равенство функций на границе включения;
		{ 0.,			HH2,  0., -2.,		 -HHL2,		-JJL2, 0., 0., 0., 0., 0., 0. }, //...равенство нормальных производных на границе включения;
		{ 0., -ku_I*HH1,  0.,  0., ku_L*HHL1, ku_L*JJL1, 0., 0., 0., 0., 0., 0. }, //...равенство моментов на границе включения;
		{ K_I+mu_L, (mu_I-mu_L+ku_I*.5)*HH1, -ku_L, 0., -ku_L*.5*HHL1, -ku_L*.5*JJL1, 0., 0., 0., 0., 0., 0.}, //...поверхностные силы;
		{ 0.,	0., 1.,     c1, -HHM1, -JJM1, -1., -1.,  HHK1,  JJK1, 0., 0. }, //...равенство функций на границе слоя;
		{ 0.,	0., 0.,  2.*c1,  HHM2,  JJM2,  0., -2., -HHK2, -JJK2, 0., 0. }, //...равенство нормальных производных на границе слоя;
		{ 0., 0., AA, -AA*c1, -ku_L*HHM1-AA*(HHM1+HHM2), -ku_L*JJM1-AA*(JJM1+JJM2), 0., 0., ku_M*HHK1, ku_M*JJK1, 0., 0. }, //...равенство моментов на границе слоя;
		{ 0., 0., K_L, -mu_L*c1, (mu_L+ku_L*.5)*HHM1, (mu_L+ku_L*.5)*JJM1, -K_M, mu_M, -(mu_M+ku_M*.5)*HHK1, -(mu_M+ku_M*.5)*JJK1, 0., 0. }, //...поверхностные силы;
		{ 0., 0., 0., 0., 0., 0., 1., c0, -HHP1, -JJP1, 0., 1. }, //...перемещения на границе эффективной области;
#ifdef __NORMAL_DERIV_SECOND__
//		{ 0., 0., 0., 0., 0., 0., 1.,   -c0, -HHP1-HHP2, -JJP1-JJP2, 0., 0. }, 
		{ 0., 0., 0., 0., 0., 0., 0., 2.*c0,  HHP2, JJP2, 0., 1. }, //...равенство нулю нормальной производной на границе матрицы;
#else
		{ 0., 0., 0., 0., 0., 0., 0., 0.,  HHP1,  JJP1, 0., 0. }, //...равенство нулю когезионного поля на границе матрицы;
#endif
//		{ 0., 0., 0., 0., 0., 0., K_M, -mu_M*c0, (mu_M+ku_M*.5)*HHP1, (mu_M+ku_M*.5)*JJP1, -1., 0.},
		{ 0., 0., 0., 0., 0., 0., ku_M,  0., ku_M*.5*HHP1, ku_M*.5*JJP1, -1., mu_M }, //...эффективный модуль;
	};

////////////////////////////////////////////////////////////////////////////////////////////
//...решаем систему линейных уравнений  A0, A0*, A1, A1^, A1*, A1*^, A2, A2^, A2*, A2*^, KH;
	int dim_N = 11, ii[11] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, i, k, l, k0, l0;
	for (i = 0; i < dim_N; i++) {
		double f = 0.;
///////////////////////////////////////
//...look for position maximal element;
		for (k = 0; k < dim_N; k++)
			if (ii[k] != 1) 
				for (l = 0; l < dim_N; l++) 
					if (! ii[l]) {
						if (fabs(matr[k][l]) >= f) f = fabs(matr[k0 = k][l0 = l]); 
					}
					else if (ii[l] > 1) return(0.);
		++(ii[l0]);
///////////////////////////////////////////////////////////
//...swapping row for diagonal position of maximal element;
		if (k0 != l0) 
			for (l = 0; l <= dim_N; l++) {
				f = matr[k0][l]; matr[k0][l] = matr[l0][l]; matr[l0][l] = f; 
			}
		if (matr[l0][l0] == 0.) return(0.);
////////////////////////////////
//...diagonal row normalization;
		double finv = 1./matr[l0][l0]; matr[l0][l0] = 1.;
		for (l = 0; l <= dim_N; l++) matr[l0][l] *= finv;
/////////////////////////////////
//...elimination all outher rows;
		for (k = 0; k < dim_N; k++)
			if ( k != l0) {
				finv = matr[k][l0]; matr[k][l0] = 0.;
				for (l = 0; l <= dim_N; l++) matr[k][l] -= matr[l0][l]*finv;
			}
	}
	for (l = 0; l < dim_N; l++) KK[l] = matr[l][dim_N];
	return (matr[dim_N-1][dim_N]);
}

////////////////////////////////////////////////////////////////////
//...трехфазная модель для цилиндрического включения (прямой алгоритм);
double CCohes2D::TakeEshelby_volm_two(double ff)
{
	double mu_I = get_param(NUM_SHEAR+NUM_SHIFT), nu_I = get_param(NUM_SHEAR+1+NUM_SHIFT), C_I = get_param(NUM_SHEAR-1+NUM_SHIFT),
			 mu_M = get_param(NUM_SHEAR), nu_M = get_param(NUM_SHEAR+1), C_M = get_param(NUM_SHEAR-1), AA = get_param(NUM_ADHES),
			 K_I = mu_I/(1.-2.*nu_I), K_M = mu_M/(1.-2.*nu_M);
	if (C_I == 0. || C_M == 0.) { //...классическая трехфазная модель;
		return(K_M+ff*(K_I-K_M)/(1.+(1.-ff)*(K_I-K_M)/(K_M+mu_M)));
	}
	double ku_I = (1.-nu_I)/(.5-nu_I)*mu_I, kk_I = sqrt(C_I/ku_I),
			 ku_M = (1.-nu_M)/(.5-nu_M)*mu_M, kk_M = sqrt(C_M/ku_M), RR1 = 1./sqrt(ff), hh[3];
	scale_regul(1., kk_I, hh);
	double HH1 = hh[1], HH2 = hh[2];
	scale_regul(1., kk_M, hh);
	double HHL1 = hh[1], HHL2 = hh[2];
	scale_irreg(1., kk_M, hh);
	double JJL1 = hh[1], JJL2 = hh[2], HHM1, JJM1, HHM2, JJM2;
	scale_regul(RR1, kk_M, hh); HHM1 = hh[1]; HHM2 = hh[2];
	scale_irreg(RR1, kk_M, hh); JJM1 = hh[1]; JJM2 = hh[2];
	double matr[7][8] = {
		{ 1., -HH1, -1., -1.,  HHL1,  JJL1, 0., 0. }, //...равенство функций;
		{ 0.,  HH2,  0., -2., -HHL2, -JJL2, 0., 0. }, //...равенство нормальных производных;
		{ AA, -ku_I*HH1-AA*(HH1+HH2), 0.,0., ku_M*HHL1, ku_M*JJL1, 0., 0.}, //...равенство моментов;
		{ K_I+mu_M, (mu_I-mu_M+ku_I*.5)*HH1, -ku_M, 0., -ku_M*.5*HHL1, -ku_M*.5*JJL1, 0., 0.}, //...поверхностные силы;
		{ 0., 0., 1.,    ff, -HHM1, -JJM1, 0., 1. }, //...перемещения на границе эффективной области;
#ifdef __NORMAL_DERIV_SECOND__
//		{ 0., 0., 1.,   -ff, -HHM1-HHM2, -JJM1-JJM2, 0., 0. }, 
		{ 0., 0., 0., 2.*ff,  HHM2,  JJM2, 0., 1. }, //...равенство нулю нормальной производной на границе матрицы;
#else
		{ 0., 0., 0.,    0.,  HHM1,  JJM1, 0., 0. }, //...равенство нулю когезионного поля на границе матрицы;
#endif
//		{ 0., 0., K_M, -mu_M*ff, (mu_M+ku_M*.5)*HHM1, (mu_M+ku_M*.5)*JJM1, -1., 0.},
		{ 0., 0., ku_M,  0., ku_M*.5*HHM1, ku_M*.5*JJM1, -1., mu_M }, //...эффективный модуль;
	};

////////////////////////////////////////////////////////////////////////
//...решаем систему линейных уравнений  A0, A0*, A1, A1^, A1*, A1*^, KH;
	int dim_N = 7, ii[7] = {0, 0, 0, 0, 0, 0, 0}, i, k, l, k0, l0;
	for (i = 0; i < dim_N; i++) {
		double f = 0.;
///////////////////////////////////////
//...look for position maximal element;
		for (k = 0; k < dim_N; k++)
			if (ii[k] != 1) 
				for (l = 0; l < dim_N; l++) 
					if (! ii[l]) {
						if (fabs(matr[k][l]) >= f) f = fabs(matr[k0 = k][l0 = l]); 
					}
					else if (ii[l] > 1) return(0.);
		++(ii[l0]);
///////////////////////////////////////////////////////////
//...swapping row for diagonal position of maximal element;
		if (k0 != l0) 
			for (l = 0; l <= dim_N; l++) {
				f = matr[k0][l]; matr[k0][l] = matr[l0][l]; matr[l0][l] = f; 
			}
		if (matr[l0][l0] == 0.) return(0.);
////////////////////////////////
//...diagonal row normalization;
		double finv = 1./matr[l0][l0]; matr[l0][l0] = 1.;
		for (l = 0; l <= dim_N; l++) matr[l0][l] *= finv;
/////////////////////////////////
//...elimination all outher rows;
		for (k = 0; k < dim_N; k++)
			if ( k != l0) {
				finv = matr[k][l0]; matr[k][l0] = 0.;
				for (l = 0; l <= dim_N; l++) matr[k][l] -= matr[l0][l]*finv;
			}
	}
	for (l = 0; l < dim_N; l++) KK[l] = matr[l][dim_N];
	return (matr[dim_N-1][dim_N]);
}

//////////////////////////////////////////////////////////////////////////////
//...алгоритмическое решение системы уравнений Эшелби для многослойной модели;
double take_system_cyl(double ** matrix, int * ii, int dim_N, double mu_H, double nu_H);

////////////////////////////////////////////////////////////////////////////////
//...четырехфазная модель для цилиндрического включения (итерационный алгоритм);
double CCohes2D::TakeEshelby_shear(double ff, double ff_l, double eps, int max_iter)
{
	double mu_I = get_param(NUM_SHEAR+NUM_SHIFT), nu_I = get_param(NUM_SHEAR+1+NUM_SHIFT), C_I = get_param(NUM_SHEAR-1+NUM_SHIFT),
			 mu_M = get_param(NUM_SHEAR), nu_M = get_param(NUM_SHEAR+1), C_M = get_param(NUM_SHEAR-1),
			 mu_L = get_param(NUM_SHEAR+NUM_SHIFT*2), nu_L = get_param(NUM_SHEAR+1+NUM_SHIFT*2), C_L = get_param(NUM_SHEAR-1+NUM_SHIFT*2),
			 K_I  = mu_I/(1.-2.*nu_I), K_M = mu_M/(1.-2.*nu_M), K_L = mu_L/(1.-2.*nu_L), c0 = ff+ff_l, c1 = ff/c0,
			 ku_I = (1.-nu_I)/(.5-nu_I)*mu_I, kk_I = sqrt(C_I/ku_I), kp_I = sqrt(C_I/mu_I), AA = get_param(NUM_ADHES), 
			 ku_M = (1.-nu_M)/(.5-nu_M)*mu_M, kk_M = sqrt(C_M/ku_M), kp_M = sqrt(C_M/mu_M), BB = get_param(NUM_ADHES+1), 
			 ku_L = (1.-nu_L)/(.5-nu_L)*mu_L, kk_L = sqrt(C_L/ku_L), kp_L = sqrt(C_L/mu_L), RR1 = 1./sqrt(c1), RR2 = 1./sqrt(ff), hh[3];
	scale_regul(1., kk_I, hh);
	double HH1 = hh[1], HH2 = hh[2];
	scale_regul(1., kk_L, hh);
	double HHL1 = hh[1], HHL2 = hh[2];
	scale_irreg(1., kk_L, hh);
	double JJL1 = hh[1], JJL2 = hh[2], HHM1, JJM1, HHM2, JJM2;
	scale_regul(RR1, kk_L, hh); HHM1 = hh[1]; HHM2 = hh[2];
	scale_irreg(RR1, kk_L, hh); JJM1 = hh[1]; JJM2 = hh[2];
	double HHK1, JJK1, HHK2, JJK2;
	scale_regul(RR1, kk_M, hh); HHK1 = hh[1]; HHK2 = hh[2];
	scale_irreg(RR1, kk_M, hh); JJK1 = hh[1]; JJK2 = hh[2];
	double HHP1, JJP1, HHP2, JJP2;
	scale_regul(RR2, kk_M, hh); HHP1 = hh[1]; HHP2 = hh[2];
	scale_irreg(RR2, kk_M, hh); JJP1 = hh[1]; JJP2 = hh[2];
/////////////////////////////////////////////////////////////////////////
//...дополнительные масштабные функции для сдвигового масштабного модуля;
	scale_regul(1., kp_I, hh);
	double HP1 = hh[1], HP2 = hh[2];
	scale_regul(1., kp_L, hh);
	double HPL1 = hh[1], HPL2 = hh[2];
	scale_irreg(1., kp_L, hh);
	double JPL1 = hh[1], JPL2 = hh[2], HPM1, JPM1, HPM2, JPM2;
	scale_regul(RR1, kp_L, hh); HPM1 = hh[1]; HPM2 = hh[2];
	scale_irreg(RR1, kp_L, hh); JPM1 = hh[1]; JPM2 = hh[2];
	double HPK1, JPK1, HPK2, JPK2;
	scale_regul(RR1, kp_M, hh); HPK1 = hh[1]; HPK2 = hh[2];
	scale_irreg(RR1, kp_M, hh); JPK1 = hh[1]; JPK2 = hh[2];
	double HPP1, JPP1, HPP2, JPP2;
	scale_regul(RR2, kp_M, hh); HPP1 = hh[1]; HPP2 = hh[2];
	scale_irreg(RR2, kp_M, hh); JPP1 = hh[1]; JPP2 = hh[2];

////////////////////////////////////////////////////////////////////////////
//...заполняем матрицу для определения модуля сдвига в поперечной плоскости;
	double matr[22][23] ={
//...равенство функций на границе включения;
		{ 1., -4.*nu_I, -HH1/ku_I+2.*HH2/C_I, -2.*HP2/C_I, -1., 4.*nu_L, nu_L-1., -2., HHL1/ku_L-2.*HHL2/C_L, 2.*HPL2/C_L, JJL1/ku_L-2.*JJL2/C_L, 2.*JPL2/C_L, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. },
		{ 1., -2.*(3.-2.*nu_I), -2.*HH2/C_I, -HP1/mu_I+2.*HP2/C_I, -1., 2.*(3.-2.*nu_L), nu_L-.5, 2., 2.*HHL2/C_L, HPL1/mu_L-2.*HPL2/C_L, 2.*JJL2/C_L, JPL1/mu_L-2.*JPL2/C_L, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. },
//...поверхностные силы на границе включения;
		{ mu_I, 0., (mu_I/ku_I+.5)*HH1-(4.*mu_I+ku_I)*HH2/C_I, -HP1+(4.*mu_I+ku_I)*HP2/C_I, -mu_L, 0., mu_L, 6.*mu_L, -(mu_L/ku_L+.5)*HHL1+(4.*mu_L+ku_L)*HHL2/C_L, HPL1-(4.*mu_L+ku_L)*HPL2/C_L, -(mu_L/ku_L+.5)*JJL1+(4.*mu_L+ku_L)*JJL2/C_L, JPL1-(4.*mu_L+ku_L)*JPL2/C_L, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. },
		{ mu_I, -6.*mu_I, (1.-2.*mu_I/ku_I)*HH1+(7.*mu_I-2.*ku_I)*HH2/C_I, 1.5*HP1-(7.*mu_I-2.*ku_I)*HP2/C_I, -mu_L, 6.*mu_L, -mu_L*.5, -6.*mu_L, -(1.-2.*mu_L/ku_L)*HHL1-(7.*mu_L-2.*ku_L)*HHL2/C_L, -1.5*HPL1+(7.*mu_L-2.*ku_L)*HPL2/C_L, -(1.-2.*mu_L/ku_L)*JJL1-(7.*mu_L-2.*ku_L)*JJL2/C_L, -1.5*JPL1+(7.*mu_L-2.*ku_L)*JPL2/C_L, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. },
//...равенство моментов на границе включения;
		{ 0., 0., -ku_I*(HH1/ku_I-2.*HH2/C_I), -2.*ku_I*HP2/C_I, 0., 0., 0., 0., ku_L*(HHL1/ku_L-2.*HHL2/C_L), 2.*ku_L*HPL2/C_L, ku_L*(JJL1/ku_L-2.*JJL2/C_L), 2.*ku_L*JPL2/C_L, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. },
		{ 0., 0., -2.*mu_I*HH2/C_I, -mu_I*(HP1/mu_I-2.*HP2/C_I), 0., 0., 0., 0., 2.*mu_L*HHL2/C_L, mu_L*(HPL1/mu_L-2.*HPL2/C_L), 2.*mu_L*JJL2/C_L, mu_L*(JPL1/mu_L-2.*JPL2/C_L), 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. },
//...равенство нормальных производных на границе включения;
		{ 1., -12.*nu_I, (HH1-(1.+6.*ku_I/C_I)*HH2)/ku_I, -2.*HP1/mu_I+6.*HP2/C_I, -1., 12.*nu_L, 1.-nu_L, 6., -(HHL1-(1.+6.*ku_L/C_L)*HHL2)/ku_L, 2.*HPL1/mu_L-6.*HPL2/C_L, -(JJL1-(1.+6.*ku_L/C_L)*JJL2)/ku_L, 2.*JPL1/mu_L-6.*JPL2/C_L, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. },
		{ 1., -6.*(3.-2.*nu_I), -2.*HH1/ku_I+6.*HH2/C_I, (HP1-(1.+6.*mu_I/C_I)*HP2)/mu_I, -1., 6.*(3.-2.*nu_L), .5-nu_L, -6., 2.*HHL1/ku_L-6.*HHL2/C_L, -(HPL1-(1.+6.*mu_L/C_L)*HPL2)/mu_L, 2.*JJL1/ku_L-6.*JJL2/C_L, -(JPL1-(1.+6.*mu_L/C_L)*JPL2)/mu_L, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. },
//...равенство функций на границе слоя;
		{ 0., 0., 0., 0., 1., -4.*nu_L/c1, (1.-nu_L)*c1, 2.*sqr(c1), -HHM1/ku_L+2.*c1*HHM2/C_L, -2.*c1*HPM2/C_L, -JJM1/ku_L+2.*c1*JJM2/C_L, -2.*c1*JPM2/C_L, -1., 4.*nu_M, nu_M-1., -2., HHK1/ku_M-2.*c1*HHK2/C_M, 2.*c1*HPK2/C_M, JJK1/ku_M-2.*c1*JJK2/C_M, 2.*c1*JPK2/C_M, 0., 0., 0. },
		{ 0., 0., 0., 0., 1., -2.*(3.-2.*nu_L)/c1, (.5-nu_L)*c1, -2.*sqr(c1), -2.*c1*HHM2/C_L, -HPM1/mu_L+2.*c1*HPM2/C_L, -2.*c1*JJM2/C_L, -JPM1/mu_L+2.*c1*JPM2/C_L, -1., 2.*(3.-2.*nu_M), nu_M-.5, 2., 2.*c1*HHK2/C_M, HPK1/mu_M-2.*c1*HPK2/C_M, 2.*c1*JJK2/C_M, JPK1/mu_M-2.*c1*JPK2/C_M, 0., 0., 0. },
//...поверхностные силы на границе слоя;
		{ 0., 0., 0., 0., mu_L, 0., -mu_L*c1, -6.*mu_L*sqr(c1), (mu_L/ku_L+.5)*HHM1-(4.*mu_L+ku_L)*c1*HHM2/C_L, -HPM1+(4.*mu_L+ku_L)*c1*HPM2/C_L, (mu_L/ku_L+.5)*JJM1-(4.*mu_L+ku_L)*c1*JJM2/C_L, -JPM1+(4.*mu_L+ku_L)*c1*JPM2/C_L, -mu_M, 0., mu_M, 6.*mu_M, -(mu_M/ku_M+.5)*HHK1+(4.*mu_M+ku_M)*c1*HHK2/C_M, HPK1-(4.*mu_M+ku_M)*c1*HPK2/C_M, -(mu_M/ku_M+.5)*JJK1+(4.*mu_M+ku_M)*c1*JJK2/C_M, JPK1-(4.*mu_M+ku_M)*c1*JPK2/C_M, 0., 0., 0. },
		{ 0., 0., 0., 0., mu_L, -6.*mu_L/c1, mu_L*.5*c1, 6.*mu_L*sqr(c1), (1.-2.*mu_L/ku_L)*HHM1+(7.*mu_L-2.*ku_L)*c1*HHM2/C_L, 1.5*HPM1-(7.*mu_L-2.*ku_L)*c1*HPM2/C_L, (1.-2.*mu_L/ku_L)*JJM1+(7.*mu_L-2.*ku_L)*c1*JJM2/C_L, 1.5*JPM1-(7.*mu_L-2.*ku_L)*c1*JPM2/C_L, -mu_M, 6.*mu_M, -mu_M*.5, -6.*mu_M, -(1.-2.*mu_M/ku_M)*HHK1-(7.*mu_M-2.*ku_M)*c1*HHK2/C_M, -1.5*HPK1+(7.*mu_M-2.*ku_M)*c1*HPK2/C_M, -(1.-2.*mu_M/ku_M)*JJK1-(7.*mu_M-2.*ku_M)*c1*JJK2/C_M, -1.5*JPK1+(7.*mu_M-2.*ku_M)*c1*JPK2/C_M, 0., 0., 0. },
//...равенство моментов на границе слоя;
		{ 0., 0., 0., 0., AA, -12.*AA*nu_L/c1, AA*(nu_L-1.)*c1, -6.*AA*sqr(c1), -ku_L*(HHM1/ku_L-2.*c1*HHM2/C_L)+AA*(HHM1-(1.+6.*c1*ku_L/C_L)*HHM2)/ku_L, -2.*ku_L*c1*HPM2/C_L-AA*(2.*HPM1/mu_L-6.*c1*HPM2/C_L), -ku_L*(JJM1/ku_L-2.*c1*JJM2/C_L)+AA*(JJM1-(1.+6.*c1*ku_L/C_L)*JJM2)/ku_L, -2.*ku_L*c1*JPM2/C_L-AA*(2.*JPM1/mu_L-6.*c1*JPM2/C_L), 0., 0., 0., 0., ku_M*(HHK1/ku_M-2.*c1*HHK2/C_M), 2.*ku_M*c1*HPK2/C_M, ku_M*(JJK1/ku_M-2.*c1*JJK2/C_M), 2.*ku_M*c1*JPK2/C_M, 0., 0., 0. },
		{ 0., 0., 0., 0., BB, -6.*BB*(3.-2.*nu_L)/c1, BB*(nu_L-.5)*c1, 6.*BB*sqr(c1), -2.*mu_L*c1*HHM2/C_L-BB*(2.*HHM1/ku_L-6.*c1*HHM2/C_L), -mu_L*(HPM1/mu_L-2.*c1*HPM2/C_L)+BB*(HPM1-(1.+6.*c1*mu_L/C_L)*HPM2)/mu_L, -2.*mu_L*c1*JJM2/C_L-BB*(2.*JJM1/ku_L-6.*c1*JJM2/C_L), -mu_L*(JPM1/mu_L-2.*c1*JPM2/C_L)+BB*(JPM1-(1.+6.*c1*mu_L/C_L)*JPM2)/mu_L, 0., 0., 0., 0., 2.*mu_M*c1*HHK2/C_M, mu_M*(HPK1/mu_M-2.*c1*HPK2/C_M), 2.*mu_M*c1*JJK2/C_M, mu_M*(JPK1/mu_M-2.*c1*JPK2/C_M), 0., 0., 0. },
//...равенство нормальных производных на границе слоя;
		{ 0., 0., 0., 0., 1., -12.*nu_L/c1, (nu_L-1.)*c1, -6.*sqr(c1), (HHM1-(1.+6.*c1*ku_L/C_L)*HHM2)/ku_L, -2.*HPM1/mu_L+6.*c1*HPM2/C_L, (JJM1-(1.+6.*c1*ku_L/C_L)*JJM2)/ku_L, -2.*JPM1/mu_L+6.*c1*JPM2/C_L, -1., 12.*nu_M, 1.-nu_M, 6., -(HHK1-(1.+6.*c1*ku_M/C_M)*HHK2)/ku_M, 2.*HPK1/mu_M-6.*c1*HPK2/C_M, -(JJK1-(1.+6.*c1*ku_M/C_M)*JJK2)/ku_M, 2.*JPK1/mu_M-6.*c1*JPK2/C_M, 0., 0., 0. },
		{ 0., 0., 0., 0., 1., -6.*(3.-2.*nu_L)/c1, (nu_L-.5)*c1, 6.*sqr(c1), -2.*HHM1/ku_L+6.*c1*HHM2/C_L, (HPM1-(1.+6.*c1*mu_L/C_L)*HPM2)/mu_L, -2.*JJM1/ku_L+6.*c1*JJM2/C_L, (JPM1-(1.+6.*c1*mu_L/C_L)*JPM2)/mu_L, -1., 6.*(3.-2.*nu_M), .5-nu_M, -6., 2.*HHK1/ku_M-6.*c1*HHK2/C_M, -(HPK1-(1.+6.*c1*mu_M/C_M)*HPK2)/mu_M, 2.*JJK1/ku_M-6.*c1*JJK2/C_M, -(JPK1-(1.+6.*c1*mu_M/C_M)*JPK2)/mu_M, 0., 0., 0. },
#ifdef __NORMAL_DERIV_SECOND__ //...равенство нулю нормальной производной на границе слоя;
		{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., -12.*nu_M/c0, (nu_M-1.)*c0, -6.*sqr(c0), (HHP1-(1.+6.*c0*ku_M/C_M)*HHP2)/ku_M, -2.*HPP1/mu_M+6.*c0*HPP2/C_M, (JJP1-(1.+6.*c0*ku_M/C_M)*JJP2)/ku_M, -2.*JPP1/mu_M+6.*c0*JPP2/C_M, 0., 0., 0. },
		{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., -6.*(3.-2.*nu_M)/c0, (nu_M-.5)*c0, 6.*sqr(c0), -2.*HHP1/ku_M+6.*c0*HHP2/C_M, (HPP1-(1.+6.*c0*mu_M/C_M)*HPP2)/mu_M, -2.*JJP1/ku_M+6.*c0*JJP2/C_M, (JPP1-(1.+6.*c0*mu_M/C_M)*JPP2)/mu_M, 0., 0., 0. },
#else //...равенство нулю когезионного поля на границе слоя;
		{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., HHP1/ku_M-2.*c0*HHP2/C_M, 2.*c0*HPP2/C_M, JJP1/ku_M-2.*c0*JJP2/C_M, 2.*c0*JPP2/C_M, 0., 0., 0. },
		{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 2.*c0*HHP2/C_M, HPP1/mu_M-2.*c0*HPP2/C_M, 2.*c0*JJP2/C_M, JPP1/mu_M-2.*c0*JPP2/C_M, 0., 0., 0. },
#endif
//...перемещения на границе эффективной области;
		{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., -4.*nu_M/c0, (1.-nu_M)*c0, 2.*sqr(c0), -HHP1/ku_M+2.*c0*HHP2/C_M, -2.*c0*HPP2/C_M, -JJP1/ku_M+2.*c0*JJP2/C_M, -2.*c0*JPP2/C_M, -2., -1. , 1. },
		{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., -2.*(3.-2.*nu_M)/c0, (.5-nu_M)*c0, -2.*sqr(c0), -2.*c0*HHP2/C_M, -HPP1/mu_M+2.*c0*HPP2/C_M, -2.*c0*JJP2/C_M, -JPP1/mu_M+2.*c0*JPP2/C_M, 2., -1. , 1.},
//...поверхностные силы на границе эффективной области;
		{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., mu_M, 0., -mu_M*c0, -6.*mu_M*sqr(c0), (mu_M/ku_M+.5)*HHP1-(4.*mu_M+ku_M)*c0*HHP2/C_M, -HPP1+(4.*mu_M+ku_M)*c0*HPP2/C_M, (mu_M/ku_M+.5)*JJP1-(4.*mu_M+ku_M)*c0*JJP2/C_M, -JPP1+(4.*mu_M+ku_M)*c0*JPP2/C_M, 6., 1., 1. },
		{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., mu_M, -6.*mu_M/c0, mu_M*.5*c0, 6.*mu_M*sqr(c0), (1.-2.*mu_M/ku_M)*HHP1+(7.*mu_M-2.*ku_M)*c0*HHP2/C_M, 1.5*HPP1-(7.*mu_M-2.*ku_M)*c0*HPP2/C_M, (1.-2.*mu_M/ku_M)*JJP1+(7.*mu_M-2.*ku_M)*c0*JJP2/C_M, 1.5*JPP1-(7.*mu_M-2.*ku_M)*c0*JPP2/C_M, -6., -0.5, 1. },
	};

/////////////////////////////////////////////////////////////////////////////
//...итерационный алгоритм вычисления модуля сдвига в поперечном направлении;
	double optim, sgn0, sgn1, nu_H, mu_H, mu_H0, mu_H1/*, KH = TakeEshelby_volm_two(ff)*/,
		** matrix = NULL; set_matrix(matrix, 22, 23); int * ii = new_struct<int>(22), k_iter = 0, k, l, m = 22;

/////////////////////////////////////////////
//...цикл по определению сдвиговой жесткости;
	mu_H1 = mu_I*ff+mu_M*(1.-ff); mu_H1 *= 10.;
	nu_H = /*0.5*(1.-mu_H1/KH)*/0.3;
	for (k = 0; k < m; k++)
		for (l = 0; l < m+1; l++) matrix[k][l] = matr[k][l];
	optim = take_system_cyl(matrix, ii, m, mu_H1, nu_H);
	if (optim > 0.) sgn1 = 1.; else sgn1 = -1.;

	mu_H0 = mu_H1;
	do {
		mu_H0 /= 2.; k_iter++;
		nu_H = /*0.5*(1.-mu_H0/KH)*/0.3;
		for (k = 0; k < m; k++)
			for (l = 0; l < m+1; l++) matrix[k][l] = matr[k][l];
		optim = take_system_cyl(matrix, ii, m, mu_H0, nu_H);
	} while (optim*sgn1 > 0. && k_iter < max_iter);
	if (optim > 0.) sgn0 = 1.; else sgn0 = -1.; k_iter = 0;

	do {
		mu_H = (mu_H0+mu_H1)*.5; k_iter++;
		nu_H = /*0.5*(1.-mu_H/KH)*/0.3;
		for (k = 0; k < m; k++)
			for (l = 0; l < m+1; l++)	matrix[k][l] = matr[k][l];
		optim = take_system_cyl(matrix, ii, m, mu_H, nu_H);
		if (optim*sgn0 > 0.) mu_H0 = mu_H; else
			if (optim*sgn0 < 0.) mu_H1 = mu_H; else mu_H0 = mu_H1 = mu_H;
	} while (fabs(mu_H1-mu_H0) > eps && k_iter < max_iter);
	mu_H = (mu_H0+mu_H1)*.5;

	//////////////////////////////////
	//...вычисляем эффективные модули;
	if (sgn0*sgn1 > 0. || fabs(optim) > 1e-6) mu_H = -mu_H;
	for (l = 0; l < m; l++) GG[l] = matrix[l][m];
	delete_struct(matrix); delete_struct(ii);

	for (optim = -1., l = 0; l < m; l++) optim += GG[l]*matr[m-4][l];
	for (optim = -1., l = 0; l < m; l++) optim += GG[l]*matr[m-3][l];

	return(mu_H);
}

/////////////////////////////////////////////////////////////////////////////
//...трехфазная модель для цилиндрического включения (итерационный алгоритм);
double CCohes2D::TakeEshelby_shear_two(double ff, double & det, double eps, int max_iter)
{
	double mu_M = get_param(NUM_SHEAR), nu_M = get_param(NUM_SHEAR+1), C_M = get_param(NUM_SHEAR-1),
		mu_I = get_param(NUM_SHEAR+NUM_SHIFT), nu_I = get_param(NUM_SHEAR+1+NUM_SHIFT), C_I = get_param(NUM_SHEAR-1+NUM_SHIFT),
		K_M = mu_M/(1.-2.*nu_M), K_I = mu_I/(1.-2.*nu_I), AA = get_param(NUM_ADHES), BB = get_param(NUM_ADHES+1), c0 = ff,
		ku_M = (1.-nu_M)/(.5-nu_M)*mu_M, lm_M = ku_M-2.*mu_M, kk_M = sqrt(C_M/ku_M), kp_M = sqrt(C_M/mu_M),
		ku_I = (1.-nu_I)/(.5-nu_I)*mu_I, lm_I = ku_I-2.*mu_I, kk_I = sqrt(C_I/ku_I), kp_I = sqrt(C_I/mu_I), RR1 = 1./sqrt(c0), hh[3];
	scale_regul(1., kk_I, hh);
	double HH1 = hh[1], HH2 = hh[2];
	scale_regul(1., kk_M, hh);
	double HHL1 = hh[1], HHL2 = hh[2];
	scale_irreg(1., kk_M, hh);
	double JJL1 = hh[1], JJL2 = hh[2], HHM1, JJM1, HHM2, JJM2;
	scale_regul(RR1, kk_M, hh); HHM1 = hh[1]; HHM2 = hh[2];
	scale_irreg(RR1, kk_M, hh); JJM1 = hh[1]; JJM2 = hh[2];
/////////////////////////////////////////////////////////////////////////
//...дополнительные масштабные функции для сдвигового масштабного модуля;
	scale_regul(1., kp_I, hh);
	double HP1 = hh[1], HP2 = hh[2];
	scale_regul(1., kp_M, hh);
	double HPL1 = hh[1], HPL2 = hh[2];
	scale_irreg(1., kp_M, hh);
	double JPL1 = hh[1], JPL2 = hh[2], HPM1, JPM1, HPM2, JPM2;
	scale_regul(RR1, kp_M, hh); HPM1 = hh[1]; HPM2 = hh[2];
	scale_irreg(RR1, kp_M, hh); JPM1 = hh[1]; JPM2 = hh[2];

////////////////////////////////////////////////////////////////////////////
//...заполняем матрицу для определения модуля сдвига в поперечной плоскости;
	double matr[14][15] ={
//...равенство функций;
		{ 1., -4.*nu_I, -HH1/ku_I+2.*HH2/C_I, -2.*HP2/C_I, -1., 4.*nu_M, nu_M-1., -2., HHL1/ku_M-2.*HHL2/C_M, 2.*HPL2/C_M, JJL1/ku_M-2.*JJL2/C_M, 2.*JPL2/C_M, 0., 0., 0. },
		{ 1., -2.*(3.-2.*nu_I), -2.*HH2/C_I, -HP1/mu_I+2.*HP2/C_I, -1., 2.*(3.-2.*nu_M), nu_M-.5, 2., 2.*HHL2/C_M, HPL1/mu_M-2.*HPL2/C_M, 2.*JJL2/C_M, JPL1/mu_M-2.*JPL2/C_M, 0., 0., 0. },
//...поверхностные силы;
		{ mu_I, 0., (mu_I/ku_I+.5)*HH1-(4.*mu_I+ku_I)*HH2/C_I, -HP1+(4.*mu_I+ku_I)*HP2/C_I, -mu_M, 0., mu_M, 6.*mu_M, -(mu_M/ku_M+.5)*HHL1+(4.*mu_M+ku_M)*HHL2/C_M, HPL1-(4.*mu_M+ku_M)*HPL2/C_M, -(mu_M/ku_M+.5)*JJL1+(4.*mu_M+ku_M)*JJL2/C_M, JPL1-(4.*mu_M+ku_M)*JPL2/C_M, 0., 0., 0. },
		{ mu_I, -6.*mu_I, (1.-2.*mu_I/ku_I)*HH1+(7.*mu_I-2.*ku_I)*HH2/C_I, 1.5*HP1-(7.*mu_I-2.*ku_I)*HP2/C_I, -mu_M, 6.*mu_M, -mu_M*.5, -6.*mu_M, -(1.-2.*mu_M/ku_M)*HHL1-(7.*mu_M-2.*ku_M)*HHL2/C_M, -1.5*HPL1+(7.*mu_M-2.*ku_M)*HPL2/C_M, -(1.-2.*mu_M/ku_M)*JJL1-(7.*mu_M-2.*ku_M)*JJL2/C_M, -1.5*JPL1+(7.*mu_M-2.*ku_M)*JPL2/C_M, 0., 0., 0. },
//...равенство моментов;
		{ AA, -12.*AA*nu_I, -ku_I*(HH1/ku_I-2.*HH2/C_I)+AA*(HH1-(1.+6.*ku_I/C_I)*HH2)/ku_I, -2.*ku_I*HP2/C_I-2.*AA*HP1/mu_I+6.*HP2/C_I, 0., 0., 0., 0., ku_M*(HHL1/ku_M-2.*HHL2/C_M), 2.*ku_M*HPL2/C_M, ku_M*(JJL1/ku_M-2.*JJL2/C_M), 2.*ku_M*JPL2/C_M, 0., 0., 0. },
		{ BB, -6.*BB*(3.-2.*nu_I), -2.*mu_I*HH2/C_I-2.*BB*HH1/ku_I+6.*HH2/C_I, -mu_I*(HP1/mu_I-2.*HP2/C_I)+BB*(HP1-(1.+6.*mu_I/C_I)*HP2)/mu_I, 0., 0., 0., 0., 2.*mu_M*HHL2/C_M, mu_M*(HPL1/mu_M-2.*HPL2/C_M), 2.*mu_M*JJL2/C_M, mu_M*(JPL1/mu_M-2.*JPL2/C_M), 0., 0., 0. },
//...равенство нормальных производных;
		{ 1., -12.*nu_I, (HH1-(1.+6.*ku_I/C_I)*HH2)/ku_I, -2.*HP1/mu_I+6.*HP2/C_I, -1., 12.*nu_M, 1.-nu_M, 6., -(HHL1-(1.+6.*ku_M/C_M)*HHL2)/ku_M, 2.*HPL1/mu_M-6.*HPL2/C_M, -(JJL1-(1.+6.*ku_M/C_M)*JJL2)/ku_M, 2.*JPL1/mu_M-6.*JPL2/C_M, 0., 0., 0. },
		{ 1., -6.*(3.-2.*nu_I), -2.*HH1/ku_I+6.*HH2/C_I, (HP1-(1.+6.*mu_I/C_I)*HP2)/mu_I, -1., 6.*(3.-2.*nu_M), .5-nu_M, -6., 2.*HHL1/ku_M-6.*HHL2/C_M, -(HPL1-(1.+6.*mu_M/C_M)*HPL2)/mu_M, 2.*JJL1/ku_M-6.*JJL2/C_M, -(JPL1-(1.+6.*mu_M/C_M)*JPL2)/mu_M, 0., 0., 0. },
#ifdef __NORMAL_DERIV_SECOND__ //...равенство нулю нормальной производной на границе матрицы;
		{ 0., 0., 0., 0., 1., -12.*nu_M/c0, (nu_M-1.)*c0, -6.*sqr(c0), (HHM1-(1.+6.*c0*ku_M/C_M)*HHM2)/ku_M, -2.*HPM1/mu_M+6.*c0*HPM2/C_M, (JJM1-(1.+6.*c0*ku_M/C_M)*JJM2)/ku_M, -2.*JPM1/mu_M+6.*c0*JPM2/C_M, 0., 0., 0. },
		{ 0., 0., 0., 0., 1., -6.*(3.-2.*nu_M)/c0, (nu_M-.5)*c0, 6.*sqr(c0), -2.*HHM1/ku_M+6.*c0*HHM2/C_M, (HPM1-(1.+6.*c0*mu_M/C_M)*HPM2)/mu_M, -2.*JJM1/ku_M+6.*c0*JJM2/C_M, (JPM1-(1.+6.*c0*mu_M/C_M)*JPM2)/mu_M, 0., 0., 0. },
#else //...равенство нулю когезионного поля на границе матрицы;
		{ 0., 0., 0., 0., 0., 0., 0., 0., HHM1/ku_M-2.*c0*HHM2/C_M, 2.*c0*HPM2/C_M, JJM1/ku_M-2.*c0*JJM2/C_M, 2.*c0*JPM2/C_M, 0. },
		{ 0., 0., 0., 0., 0., 0., 0., 0., 2.*c0*HHM2/C_M, HPM1/mu_M-2.*c0*HPM2/C_M, 2.*c0*JJM2/C_M, JPM1/mu_M-2.*c0*JPM2/C_M, 0. },
#endif
//...перемещения на границе эффективной области;
		{ 0., 0., 0., 0., 1., -4.*nu_M/c0, (1.-nu_M)*c0, 2.*sqr(c0), -HHM1/ku_M+2.*c0*HHM2/C_M, -2.*c0*HPM2/C_M, -JJM1/ku_M+2.*c0*JJM2/C_M, -2.*c0*JPM2/C_M, -2., -1. , 1. },
		{ 0., 0., 0., 0., 1., -2.*(3.-2.*nu_M)/c0, (-nu_M+.5)*c0, -2.*sqr(c0), -2.*c0*HHM2/C_M, -HPM1/mu_M+2.*c0*HPM2/C_M, -2.*c0*JJM2/C_M, -JPM1/mu_M+2.*c0*JPM2/C_M, 2., -1. , 1. },
//...поверхностные силы на границе эффективной области;
		{ 0., 0., 0., 0., mu_M, 0., -mu_M*c0, -6.*mu_M*sqr(c0), (mu_M/ku_M+.5)*HHM1-(4.*mu_M+ku_M)*c0*HHM2/C_M, -HPM1+(4.*mu_M+ku_M)*c0*HPM2/C_M, (mu_M/ku_M+.5)*JJM1-(4.*mu_M+ku_M)*c0*JJM2/C_M, -JPM1+(4.*mu_M+ku_M)*c0*JPM2/C_M, 6., 1., 1. },
		{ 0., 0., 0., 0., mu_M, -6.*mu_M/c0, mu_M*.5*c0, 6.*mu_M*sqr(c0), (1.-2.*mu_M/ku_M)*HHL1+(7.*mu_M-2.*ku_M)*c0*HHL2/C_M, 1.5*HPL1-(7.*mu_M-2.*ku_M)*c0*HPL2/C_M, (1.-2.*mu_M/ku_M)*JJL1+(7.*mu_M-2.*ku_M)*c0*JJL2/C_M, 1.5*JPL1-(7.*mu_M-2.*ku_M)*c0*JPL2/C_M, -6., -0.5, 1. },
	};

/////////////////////////////////////////////////////////////////////////////
//...итерационный алгоритм вычисления модуля сдвига в поперечном направлении;
	double optim, sgn0, sgn1, nu_H, mu_H, mu_H0, mu_H1/*, KH = TakeEshelby_volm_two(ff)*/,
		** matrix = NULL; set_matrix(matrix, 14, 15); int * ii = new_struct<int>(14), k_iter = 0, k, l, m = 14;

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//...вычисляем определитель главной матрицы (старший коэффициент квадратного уравнения для модуля сдвига);
	for (k = 0; k < m; k++)
		for (l = 0; l < m+1; l++) matrix[k][l] = matr[k][l];
	GaussJ(matrix, NULL, m-2, det);

/////////////////////////////////////////////
//...цикл по определению сдвиговой жесткости;
	mu_H1 = mu_I*ff+mu_M*(1.-ff); mu_H1 *= 10.;
	nu_H = /*0.5*(1.-mu_H1/KH)*/0.3;
	for (k = 0; k < m; k++)
		for (l = 0; l < m+1; l++) matrix[k][l] = matr[k][l];
	optim = take_system_cyl(matrix, ii, m, mu_H1, nu_H);
	if (optim > 0.) sgn1 = 1.; else sgn1 = -1.;

	mu_H0 = mu_H1;
	do {
		mu_H0 /= 2.; k_iter++;
		nu_H = /*0.5*(1.-mu_H0/KH)*/0.3;
		for (k = 0; k < m; k++)
			for (l = 0; l < m+1; l++) matrix[k][l] = matr[k][l];
		optim = take_system_cyl(matrix, ii, m, mu_H0, nu_H);
	} while (optim*sgn1 > 0. && k_iter < max_iter);
	if (optim > 0.) sgn0 = 1.; else sgn0 = -1.; k_iter = 0;

	do {
		mu_H = (mu_H0+mu_H1)*.5; k_iter++;
		nu_H = /*0.5*(1.-mu_H/KH)*/0.3;
		for (k = 0; k < m; k++)
			for (l = 0; l < m+1; l++)	matrix[k][l] = matr[k][l];
		optim = take_system_cyl(matrix, ii, m, mu_H, nu_H);
		if (optim*sgn0 > 0.) mu_H0 = mu_H; else
			if (optim*sgn0 < 0.) mu_H1 = mu_H; else mu_H0 = mu_H1 = mu_H;
	} while (fabs(mu_H1-mu_H0) > eps && k_iter < max_iter);
	mu_H = (mu_H0+mu_H1)*.5;

	//////////////////////////////////
	//...вычисляем эффективные модули;
	if (sgn0*sgn1 > 0. || fabs(optim) > 1e-6) mu_H = -mu_H;
	for (l = 0; l < m; l++) GG[l] = matrix[l][m];
	delete_struct(matrix); delete_struct(ii);

	return(mu_H);
}

///////////////////////////////////////////////////////////////////////////////////
//...модуль Юнга вдоль волокон для цилиндрического включения в градиентной матрице;
double CCohes2D::TakeGradMatrix_k0(double ff)
{
	double mu_M = get_param(NUM_SHEAR), nu_M = get_param(NUM_SHEAR+1),
		mu_I = get_param(NUM_SHEAR+NUM_SHIFT), nu_I = get_param(NUM_SHEAR+1+NUM_SHIFT),
		E_M = mu_M*2.*(1.+nu_M), E_I = mu_I*2.*(1.+nu_I), sum = 0.;
	//////////////////////////////////////////////////////////
	//...вычисляем эффективный модуль Юнга (по формуле смеси);
	sum = E_I*ff+E_M*(1.-ff);
	return(sum);
}

//////////////////////////////////////////////////////////////////////////////////
//...плоский объемный модуль для цилиндрического включения с градиентной матрицей;
double CCohes2D::TakeGradMatrix_k1(double ff)
{
	double mu_M = get_param(NUM_SHEAR), nu_M = get_param(NUM_SHEAR+1), C_M = get_param(NUM_SHEAR-1),
		mu_I = get_param(NUM_SHEAR+NUM_SHIFT), nu_I = get_param(NUM_SHEAR+1+NUM_SHIFT),
		K_M = mu_M / (1.-2.*nu_M), K_I = mu_I/(1.-2.*nu_I), ku_M = (1.-nu_M)/(.5-nu_M)*mu_M, c0 = ff;
	if (C_M == 0.) { //...классическая трехфазная модель;
		return(K_M+c0*(K_I-K_M)/(1.+(1.-c0)*(K_I-K_M)/ku_M));
	}
	double kk_M = sqrt(C_M/ku_M), RR1 = 1./sqrt(c0), hh[3];
	scale_regul(1., kk_M, hh);
	double HHL1 = hh[1], HHL2 = hh[2];
	scale_irreg(1., kk_M, hh);
	double JJL1 = hh[1], JJL2 = hh[2], HHM1, JJM1, HHM2, JJM2;
	scale_regul(RR1, kk_M, hh); HHM1 = hh[1]; HHM2 = hh[2];
	scale_irreg(RR1, kk_M, hh); JJM1 = hh[1]; JJM2 = hh[2];
	double matr[6][7] = {
		{ 1., -1., -1.,  HHL1,  JJL1, 0., 0. }, //...равенство функций;
#ifdef __NORMAL_DERIV_FIRST0__
//		{ 0., 1., -1., -HHL1-HHL2, -JJL1-JJL2, 0., 0. },
		{ 1., 0., -2., -HHL2, -JJL2, 0., 0. }, //...равенство нулю нормальной производной на границе включения;
#else
		{ 0., 0.,  0.,  HHL1,  JJL1, 0., 0. }, //...равенство нулю когезионного поля на границе включения;
#endif
		{ K_I+mu_M, -ku_M, 0., -ku_M*.5*HHL1, -ku_M*.5*JJL1, 0., 0. }, //...поверхностные силы;
		{ 0.,  1.,    c0, -HHM1, -JJM1, 0., 1. }, //...перемещения на границе эффективной области;
#ifdef __NORMAL_DERIV_SECOND__
//		{ 0.,  1.,   -c0, -HHM1-HHM2, -JJM1-JJM2, 0., 0. }, 
		{ 0.,  0., 2.*c0,  HHM2,  JJM2, 0., 1. }, //...равенство нулю нормальной производной на границе матрицы;
#else
		{ 0.,  0.,    0.,  HHM1,  JJM1, 0., 0. }, //...равенство нулю когезионного поля на границе матрицы;
#endif
		{ 0., ku_M, 0., ku_M*.5*HHM1, ku_M*.5*JJM1, -1., mu_M }, //...эффективный модуль;
};

//////////////////////////////////////////////////////////////////
//...решаем систему линейных уравнений A0, A1, A1^, A1*, A1*^, KH;
	int dim_N = 6, ii[6] = { 0, 0, 0, 0, 0, 0 }, i, k, l, k0, l0;
	for (i = 0; i < dim_N; i++) {
		double f = 0.;
///////////////////////////////////////
//...look for position maximal element;
		for (k = 0; k < dim_N; k++)
			if (ii[k] != 1)
				for (l = 0; l < dim_N; l++)
					if (!ii[l]) {
						if (fabs(matr[k][l]) >= f) f = fabs(matr[k0 = k][l0 = l]);
					}
					else if (ii[l] > 1) return(0.);
					++(ii[l0]);
///////////////////////////////////////////////////////////
//...swapping row for diagonal position of maximal element;
					if (k0 != l0)
						for (l = 0; l <= dim_N; l++) {
							f = matr[k0][l]; matr[k0][l] = matr[l0][l]; matr[l0][l] = f;
						}
					if (matr[l0][l0] == 0.) return(0.);
////////////////////////////////
//...diagonal row normalization;
					double finv = 1./matr[l0][l0]; matr[l0][l0] = 1.;
					for (l = 0; l <= dim_N; l++) matr[l0][l] *= finv;
/////////////////////////////////
//...elimination all outher rows;
					for (k = 0; k < dim_N; k++)
						if (k != l0) {
							finv = matr[k][l0]; matr[k][l0] = 0.;
							for (l = 0; l <= dim_N; l++) matr[k][l] -= matr[l0][l]*finv;
						}
	}
	for (l = 0; l < dim_N; l++) KK[l] = matr[l][dim_N];
	return(matr[dim_N-1][dim_N]);
}

///////////////////////////////////////////////////////////////////////////////////////////////
//...модуль сдвига в поперечной плоскости для цилиндрического включения с градиентной матрицей;
double CCohes2D::TakeGradMatrix_G1(double ff, double & det, double eps, int max_iter)
{
	double mu_M = get_param(NUM_SHEAR), nu_M = get_param(NUM_SHEAR+1), C_M = get_param(NUM_SHEAR-1),
		mu_I = get_param(NUM_SHEAR+NUM_SHIFT), nu_I = get_param(NUM_SHEAR+1+NUM_SHIFT),
		K_M = mu_M/(1.-2.*nu_M), K_I = mu_I/(1.-2.*nu_I), c0 = ff,
		ku_M = (1.-nu_M)/(.5-nu_M)*mu_M, lm_M = ku_M-2.*mu_M, kk_M = sqrt(C_M/ku_M), kp_M = sqrt(C_M/mu_M),
		ku_I = (1.-nu_I)/(.5-nu_I)*mu_I, lm_I = ku_I-2.*mu_I, RR1 = 1./sqrt(c0), hh[3];
	scale_regul(1., kk_M, hh);
	double HHL1 = hh[1], HHL2 = hh[2];
	scale_irreg(1., kk_M, hh);
	double JJL1 = hh[1], JJL2 = hh[2], HHM1, JJM1, HHM2, JJM2;
	scale_regul(RR1, kk_M, hh); HHM1 = hh[1]; HHM2 = hh[2];
	scale_irreg(RR1, kk_M, hh); JJM1 = hh[1]; JJM2 = hh[2];
/////////////////////////////////////////////////////////////////////////
//...дополнительные масштабные функции для сдвигового масштабного модуля;
	scale_regul(1., kp_M, hh);
	double HPL1 = hh[1], HPL2 = hh[2];
	scale_irreg(1., kp_M, hh);
	double JPL1 = hh[1], JPL2 = hh[2], HPM1, JPM1, HPM2, JPM2;
	scale_regul(RR1, kp_M, hh); HPM1 = hh[1]; HPM2 = hh[2];
	scale_irreg(RR1, kp_M, hh); JPM1 = hh[1]; JPM2 = hh[2];

////////////////////////////////////////////////////////////////////////////
//...заполняем матрицу для определения модуля сдвига в поперечной плоскости;
	double matr[12][13] ={
//...равенство функций;
		{ 1., -4.*nu_I, -1., 4.*nu_M, nu_M-1., -2., HHL1/ku_M-2.*HHL2/C_M, 2.*HPL2/C_M, JJL1/ku_M-2.*JJL2/C_M, 2.*JPL2/C_M, 0., 0., 0. },
		{ 1., -2.*(3.-2.*nu_I), -1., 2.*(3.-2.*nu_M), nu_M-.5, 2., 2.*HHL2/C_M, HPL1/mu_M-2.*HPL2/C_M, 2.*JJL2/C_M, JPL1/mu_M-2.*JPL2/C_M, 0., 0., 0. },
//...поверхностные силы;
		{ mu_I, 0., -mu_M, 0., mu_M, 6.*mu_M, -(mu_M/ku_M+.5)*HHL1+(4.*mu_M+ku_M)*HHL2/C_M, HPL1-(4.*mu_M+ku_M)*HPL2/C_M, -(mu_M/ku_M+.5)*JJL1+(4.*mu_M+ku_M)*JJL2/C_M, JPL1-(4.*mu_M+ku_M)*JPL2/C_M, 0., 0., 0. },
		{ mu_I, -6.*mu_I, -mu_M, 6.*mu_M, -mu_M*.5, -6.*mu_M, -(1.-2.*mu_M/ku_M)*HHL1-(7.*mu_M-2.*ku_M)*HHL2/C_M, -1.5*HPL1+(7.*mu_M-2.*ku_M)*HPL2/C_M, -(1.-2.*mu_M/ku_M)*JJL1-(7.*mu_M-2.*ku_M)*JJL2/C_M, -1.5*JPL1+(7.*mu_M-2.*ku_M)*JPL2/C_M, 0., 0., 0. },
#ifdef __NORMAL_DERIV_FIRST0__ //...равенство нулю нормальной производной на границе включения;
		{ 0., 0., 1., -12.*nu_M, nu_M-1., -6., (HHL1-(1.+6.*ku_M/C_M)*HHL2)/ku_M, -2.*HPL1/mu_M+6.*HPL2/C_M, (JJL1-(1.+6.*ku_M/C_M)*JJL2)/ku_M, -2.*JPL1/mu_M+6.*JPL2/C_M, 0., 0., 0. },
		{ 0., 0., 1., -6.*(3.-2.*nu_M), nu_M-.5, 6., -2.*HHL1/ku_M+6.*HHL2/C_M, (HPL1-(1.+6.*mu_M/C_M)*HPL2)/mu_M, -2.*JJL1/ku_M+6.*JJL2/C_M, (JPL1-(1.+6.*mu_M/C_M)*JPL2)/mu_M, 0., 0., 0. },
#else //...равенство нулю когезионного поля на границе включения;
		{ 0., 0., 0., 0., 0., 0., HHL1/ku_M-2.*HHL2/C_M, 2.*HPL2/C_M, JJL1/ku_M-2.*JJL2/C_M, 2.*JPL2/C_M, 0. },
		{ 0., 0., 0., 0., 0., 0., 2.*HHL2/C_M, HPL1/mu_M-2.*HPL2/C_M, 2.*JJL2/C_M, JPL1/mu_M-2.*JPL2/C_M, 0. },
#endif
#ifdef __NORMAL_DERIV_SECOND__ //...равенство нулю нормальной производной на границе матрицы;
		{ 0., 0., 1., -12.*nu_M/c0, (nu_M-1.)*c0, -6.*sqr(c0), (HHM1-(1.+6.*c0*ku_M/C_M)*HHM2)/ku_M, -2.*HPM1/mu_M+6.*c0*HPM2/C_M, (JJM1-(1.+6.*c0*ku_M/C_M)*JJM2)/ku_M, -2.*JPM1/mu_M+6.*c0*JPM2/C_M, 0., 0., 0. },
		{ 0., 0., 1., -6.*(3.-2.*nu_M)/c0, (nu_M-.5)*c0, 6.*sqr(c0), -2.*HHM1/ku_M+6.*c0*HHM2/C_M, (HPM1-(1.+6.*c0*mu_M/C_M)*HPM2)/mu_M, -2.*JJM1/ku_M+6.*c0*JJM2/C_M, (JPM1-(1.+6.*c0*mu_M/C_M)*JPM2)/mu_M, 0., 0., 0. },
#else //...равенство нулю когезионного поля на границе матрицы;
		{ 0., 0., 0., 0., 0., 0., HHM1/ku_M-2.*c0*HHM2/C_M, 2.*c0*HPM2/C_M, JJM1/ku_M-2.*c0*JJM2/C_M, 2.*c0*JPM2/C_M, 0. },
		{ 0., 0., 0., 0., 0., 0., 2.*c0*HHM2/C_M, HPM1/mu_M-2.*c0*HPM2/C_M, 2.*c0*JJM2/C_M, JPM1/mu_M-2.*c0*JPM2/C_M, 0. },
#endif
//...перемещения на границе эффективной области;
		{ 0., 0., 1., -4.*nu_M/c0, (1.-nu_M)*c0, 2.*sqr(c0), -HHM1/ku_M+2.*c0*HHM2/C_M, -2.*c0*HPM2/C_M, -JJM1/ku_M+2.*c0*JJM2/C_M, -2.*c0*JPM2/C_M, -2., -1. , 1. },
		{ 0., 0., 1., -2.*(3.-2.*nu_M)/c0, (-nu_M+.5)*c0, -2.*sqr(c0), -2.*c0*HHM2/C_M, -HPM1/mu_M+2.*c0*HPM2/C_M, -2.*c0*JJM2/C_M, -JPM1/mu_M+2.*c0*JPM2/C_M, 2., -1. , 1. },
//...поверхностные силы на границе эффективной области;
		{ 0., 0., mu_M, 0., -mu_M*c0, -6.*mu_M*sqr(c0), (mu_M/ku_M+.5)*HHM1-(4.*mu_M+ku_M)*c0*HHM2/C_M, -HPM1+(4.*mu_M+ku_M)*c0*HPM2/C_M, (mu_M/ku_M+.5)*JJM1-(4.*mu_M+ku_M)*c0*JJM2/C_M, -JPM1+(4.*mu_M+ku_M)*c0*JPM2/C_M, 6., 1., 1. },
		{ 0., 0., mu_M, -6.*mu_M/c0, mu_M*.5*c0, 6.*mu_M*sqr(c0), (1.-2.*mu_M/ku_M)*HHL1+(7.*mu_M-2.*ku_M)*c0*HHL2/C_M, 1.5*HPL1-(7.*mu_M-2.*ku_M)*c0*HPL2/C_M, (1.-2.*mu_M/ku_M)*JJL1+(7.*mu_M-2.*ku_M)*c0*JJL2/C_M, 1.5*JPL1-(7.*mu_M-2.*ku_M)*c0*JPL2/C_M, -6., -0.5, 1. },
	};

/////////////////////////////////////////////////////////////////////////////
//...итерационный алгоритм вычисления модуля сдвига в поперечном направлении;
	double optim, sgn0, sgn1, nu_H, mu_H, mu_H0, mu_H1/*, KH = TakeGradMatrix_k1(ff)*/,
		** matrix = NULL; set_matrix(matrix, 12, 13); int * ii = new_struct<int>(12), k_iter = 0, k, l, m = 12;

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//...вычисляем определитель главной матрицы (старший коэффициент квадратного уравнения для модуля сдвига);
	for (k = 0; k < m; k++)
		for (l = 0; l < m+1; l++) matrix[k][l] = matr[k][l];
	GaussJ(matrix, NULL, m-2, det);

/////////////////////////////////////////////
//...цикл по определению сдвиговой жесткости;
	mu_H1 = mu_I*ff+mu_M*(1.-ff); mu_H1 *= 10.;
	nu_H = /*0.5*(1.-mu_H1/KH)*/0.3;
	for (k = 0; k < m; k++)
		for (l = 0; l < m+1; l++) matrix[k][l] = matr[k][l];
	optim = take_system_cyl(matrix, ii, m, mu_H1, nu_H);
	if (optim > 0.) sgn1 = 1.; else sgn1 = -1.;

	mu_H0 = mu_H1;
	do {
		mu_H0 /= 2.; k_iter++;
		nu_H = /*0.5*(1.-mu_H0/KH)*/0.3;
		for (k = 0; k < m; k++)
			for (l = 0; l < m+1; l++) matrix[k][l] = matr[k][l];
		optim = take_system_cyl(matrix, ii, m, mu_H0, nu_H);
	} while (optim*sgn1 > 0. && k_iter < max_iter);
	if (optim > 0.) sgn0 = 1.; else sgn0 = -1.; k_iter = 0;

	do {
		mu_H = (mu_H0+mu_H1)*.5; k_iter++;
		nu_H = /*0.5*(1.-mu_H/KH)*/0.3;
		for (k = 0; k < m; k++)
			for (l = 0; l < m+1; l++)	matrix[k][l] = matr[k][l];
		optim = take_system_cyl(matrix, ii, m, mu_H, nu_H);
		if (optim*sgn0 > 0.) mu_H0 = mu_H; else
			if (optim*sgn0 < 0.) mu_H1 = mu_H; else mu_H0 = mu_H1 = mu_H;
	} while (fabs(mu_H1-mu_H0) > eps && k_iter < max_iter);
	mu_H = (mu_H0+mu_H1)*.5;

	//////////////////////////////////
	//...вычисляем эффективные модули;
	if (sgn0*sgn1 > 0. || fabs(optim) > 1e-6) mu_H = -mu_H;
	for (l = 0; l < m; l++) GG[l] = matrix[l][m];
	delete_struct(matrix); delete_struct(ii);

	return(mu_H);
}

////////////////////////////////////////////////////////////////////////////////////////////////
//...модуль сдвига в поперечной плоскости для цилиндрического включения с классической матрицей;
double CCohes2D::TakeGradMatrix_G1_classic(double ff, double & det, double eps, int max_iter)
{
	double mu_M = get_param(NUM_SHEAR), nu_M = get_param(NUM_SHEAR+1),
		mu_I = get_param(NUM_SHEAR+NUM_SHIFT), nu_I = get_param(NUM_SHEAR+1+NUM_SHIFT), c0 = ff, RR1 = 1./sqrt(c0);

////////////////////////////////////////////////////////////////////////////
//...заполняем матрицу для определения модуля сдвига в поперечной плоскости;
	double matr[8][9] ={
//...равенство функций;
		{ 1., -4.*nu_I, -1., 4.*nu_M, nu_M-1., -2., 0., 0., 0. },
		{ 1., -2.*(3.-2.*nu_I), -1., 2.*(3.-2.*nu_M), nu_M-.5, 2., 0., 0., 0. },
//...поверхностные силы;
		{ mu_I, 0., -mu_M, 0., mu_M, 6.*mu_M, 0., 0., 0. },
		{ mu_I, -6.*mu_I, -mu_M, 6.*mu_M, -mu_M*.5, -6.*mu_M, 0., 0., 0. },
//...перемещения на границе эффективной области;
		{ 0., 0., 1., -4.*nu_M/c0, (1.-nu_M)*c0, 2.*sqr(c0), -2., -1. , 1. },
		{ 0., 0., 1., -2.*(3.-2.*nu_M)/c0, (-nu_M+.5)*c0, -2.*sqr(c0), 2., -1. , 1. },
//...поверхностные силы на границе эффективной области;
		{ 0., 0., mu_M, 0., -mu_M*c0, -6.*mu_M*sqr(c0), 6., 1., 1. },
		{ 0., 0., mu_M, -6.*mu_M/c0, mu_M*.5*c0, 6.*mu_M*sqr(c0), -6., -0.5, 1. },
	};

/////////////////////////////////////////////////////////////////////////////
//...итерационный алгоритм вычисления модуля сдвига в поперечном направлении;
	double optim, sgn0, sgn1, nu_H, mu_H, mu_H0, mu_H1/*, KH = TakeLayer_kk(N, ff, kp, mu)*/,
		** matrix = NULL; set_matrix(matrix, 8, 9); int * ii = new_struct<int>(8), k_iter = 0, k, l, m = 8;

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//...вычисляем определитель главной матрицы (старший коэффициент квадратного уравнения для модуля сдвига);
	for (k = 0; k < m; k++)
		for (l = 0; l < m+1; l++) matrix[k][l] = matr[k][l];
	GaussJ(matrix, NULL, m-2, det);

/////////////////////////////////////////////
//...цикл по определению сдвиговой жесткости;
	mu_H1 = mu_I*ff+mu_M*(1.-ff); mu_H1 *= 10.;
	nu_H = /*0.5*(1.-mu_H1/KH)*/0.3;
	for (k = 0; k < m; k++)
		for (l = 0; l < m+1; l++) matrix[k][l] = matr[k][l];
	optim = take_system_cyl(matrix, ii, m, mu_H1, nu_H);
	if (optim > 0.) sgn1 = 1.; else sgn1 = -1.;

	mu_H0 = mu_H1;
	do {
		mu_H0 /= 2.; k_iter++;
		nu_H = /*0.5*(1.-mu_H0/KH)*/0.3;
		for (k = 0; k < m; k++)
			for (l = 0; l < m+1; l++) matrix[k][l] = matr[k][l];
		optim = take_system_cyl(matrix, ii, m, mu_H0, nu_H);
	} while (optim*sgn1 > 0. && k_iter < max_iter);
	if (optim > 0.) sgn0 = 1.; else sgn0 = -1.; k_iter = 0;

	do {
		mu_H = (mu_H0+mu_H1)*.5; k_iter++;
		nu_H = /*0.5*(1.-mu_H/KH)*/0.3;
		for (k = 0; k < m; k++)
			for (l = 0; l < m+1; l++)	matrix[k][l] = matr[k][l];
		optim = take_system_cyl(matrix, ii, m, mu_H, nu_H);
		if (optim*sgn0 > 0.) mu_H0 = mu_H; else
			if (optim*sgn0 < 0.) mu_H1 = mu_H; else mu_H0 = mu_H1 = mu_H;
	} while (fabs(mu_H1-mu_H0) > eps && k_iter < max_iter);
	mu_H = (mu_H0+mu_H1)*.5;

	//////////////////////////////////
	//...вычисляем эффективные модули;
	if (sgn0*sgn1 > 0. || fabs(optim) > 1e-6) mu_H = -mu_H;
	delete_struct(matrix); delete_struct(ii);

	return(mu_H);
}

//////////////////////////////////////////////////////////////////////////////////
//...модуль Юнга вдоль волокон  для цилиндрического включения с градиентным слоем;
double CCohes2D::TakeGradLayer_k0(double ff, double ff_l)
{
	double mu_M = get_param(NUM_SHEAR), nu_M = get_param(NUM_SHEAR+1),
		mu_I = get_param(NUM_SHEAR+NUM_SHIFT), nu_I = get_param(NUM_SHEAR+1+NUM_SHIFT),
		mu_L = get_param(NUM_SHEAR+NUM_SHIFT*2), nu_L = get_param(NUM_SHEAR+1+NUM_SHIFT*2),
		E_M = mu_M*2.*(1.+nu_M), E_I = mu_I*2.*(1.+nu_I), E_L = mu_L*2.*(1.+nu_L), sum = 0.;
//////////////////////////////////////////////////////////
//...вычисляем эффективный модуль Юнга (по формуле смеси);
	sum = E_I*ff+E_L*ff_l+E_M*(1.-ff-ff_l);
	return(sum);
}

///////////////////////////////////////////////////////////////////////////////
//...плоский объемный модуль для цилиндрического включения с градиентным слоем;
double CCohes2D::TakeGradLayer_k1(double ff, double ff_l)
{
	double mu_M = get_param(NUM_SHEAR), nu_M = get_param(NUM_SHEAR+1),
		mu_I = get_param(NUM_SHEAR+NUM_SHIFT), nu_I = get_param(NUM_SHEAR+1+NUM_SHIFT),
		mu_L = get_param(NUM_SHEAR+NUM_SHIFT*2), nu_L = get_param(NUM_SHEAR+1+NUM_SHIFT*2), C_L = get_param(NUM_SHEAR-1+NUM_SHIFT*2),
		K_M  = mu_M / (1.-2.*nu_M), K_I = mu_I/(1.-2.*nu_I), K_L = mu_L/(1.-2.*nu_L), DD, c0 = ff+ff_l, c1 = ff/c0,
		ku_M = (1.-nu_M)/(.5-nu_M)*mu_M,	ku_L = (1.-nu_L)/(.5-nu_L)*mu_L;
	if (C_L == 0.) { //...классическая четырехтрехфазная модель;
		DD = ((K_I-K_M)-(1.-c1)*(K_I-K_L)*(1.+(K_M-K_L)/ku_L))/(1.+(1.-c1)*(K_I-K_L)/ku_L);
		return(K_M+c0*DD/(1.+(1.-c0)*DD/ku_M));
	}
	double kk_L = sqrt(C_L/ku_L), RR1 = 1./sqrt(c1), hh[3];
	scale_regul(1., kk_L, hh);
	double HHL1 = hh[1], HHL2 = hh[2];
	scale_irreg(1., kk_L, hh);
	double JJL1 = hh[1], JJL2 = hh[2], HHM1, JJM1, HHM2, JJM2;
	scale_regul(RR1, kk_L, hh); HHM1 = hh[1]; HHM2 = hh[2];
	scale_irreg(RR1, kk_L, hh); JJM1 = hh[1]; JJM2 = hh[2];
	double matr[8][9] = {
		{ 1., -1., -1.,  HHL1,  JJL1, 0., 0., 0., 0. },	//...равенство функций;
#ifdef __NORMAL_DERIV_FIRST0__
//		{ 0.,  1., -1., -HHL1-HHL2, -JJL1-JJL2,  0., 0., 0., 0. },
		{ 1.,  0., -2., -HHL2, -JJL2, 0., 0., 0., 0. },	//...равенство нулю нормальных производных на границе включения;
#else
		{ 0.,  0.,  0.,  HHL1,  JJL1, 0., 0., 0., 0. }, //...равенство нулю когезионного поля на границе включения;
#endif
		{ K_I+mu_L, -ku_L, 0., -ku_L*.5*HHL1, -ku_L*.5*JJL1, 0., 0., 0., 0. }, //...поверхностные силы;
		{ 0., 1.,    c1, -HHM1, -JJM1, -1., -1., 0., 0. },	//...равенство функций на границе слоя;
#ifdef __NORMAL_DERIV_SECOND__
//		{ 0., 1., -c1,   -HHM1-HHM2, -JJM1-JJM2, 0., 0., 0., 0. }, 
		{ 0., 0., 2.*c1,  HHM2,  JJM2, -1., -1., 0., 0. },	//...равенство нулю нормальных производных на границе слоя;
#else
		{ 0., 0.,    0.,  HHM1,  JJM1,  0.,  0., 0., 0. }, //...равенство нулю когезионного поля на границе слоя;
#endif
		{ 0., ku_L, 0., ku_L*.5*HHM1, ku_L*.5*JJM1, -K_M-mu_L, mu_M-mu_L, 0., 0. }, //...поверхностные силы на границе слоя;
		{ 0., 0., 0.,  0., 0.,   1.,   c0,  0.,   1. },	//...перемещения на границе эффективной области;
		{ 0., 0., 0.,  0., 0., ku_M,   0., -1., mu_M },	//...эффективный модуль;
	};

///////////////////////////////////////////////////////////////////////////
//...решаем систему линейных уравнений A0, A1, A1^, A1*, A1*^, A2, A2^, KH;
	int dim_N = 8, ii[8] = { 0, 0, 0, 0, 0, 0, 0, 0 }, i, k, l, k0, l0;
	for (i = 0; i < dim_N; i++) {
		double f = 0.;
///////////////////////////////////////
//...look for position maximal element;
		for (k = 0; k < dim_N; k++)
			if (ii[k] != 1)
				for (l = 0; l < dim_N; l++)
					if (!ii[l]) {
						if (fabs(matr[k][l]) >= f) f = fabs(matr[k0 = k][l0 = l]);
					}
					else if (ii[l] > 1) return(0.);
					++(ii[l0]);
///////////////////////////////////////////////////////////
//...swapping row for diagonal position of maximal element;
					if (k0 != l0)
						for (l = 0; l <= dim_N; l++) {
							f = matr[k0][l]; matr[k0][l] = matr[l0][l]; matr[l0][l] = f;
						}
					if (matr[l0][l0] == 0.) return(0.);
////////////////////////////////
//...diagonal row normalization;
					double finv = 1./matr[l0][l0]; matr[l0][l0] = 1.;
					for (l = 0; l <= dim_N; l++) matr[l0][l] *= finv;
/////////////////////////////////
//...elimination all outher rows;
					for (k = 0; k < dim_N; k++)
						if (k != l0) {
							finv = matr[k][l0]; matr[k][l0] = 0.;
							for (l = 0; l <= dim_N; l++) matr[k][l] -= matr[l0][l]*finv;
						}
	}
	for (l = 0; l < dim_N; l++) KK[l] = matr[l][dim_N];
	return(matr[dim_N-1][dim_N]);
}

////////////////////////////////////////////////////////////////////////////////////////////
//...модуль сдвига в поперечной плоскости для цилиндрического включения с градиентным слоем;
double CCohes2D::TakeGradLayer_G1(double ff, double ff_l, double eps, int max_iter)
{
	double mu_M = get_param(NUM_SHEAR), nu_M = get_param(NUM_SHEAR+1),
		mu_I = get_param(NUM_SHEAR+NUM_SHIFT), nu_I = get_param(NUM_SHEAR+1+NUM_SHIFT),
		mu_L = get_param(NUM_SHEAR+NUM_SHIFT*2), nu_L = get_param(NUM_SHEAR+1+NUM_SHIFT*2), C_L = get_param(NUM_SHEAR-1+NUM_SHIFT*2),
		K_M = mu_M/(1.-2.*nu_M), K_I = mu_I/(1.-2.*nu_I), K_L = mu_L/(1.-2.*nu_L), c0 = ff+ff_l, c1 = ff/c0,
		ku_M = (1.-nu_M)/(.5-nu_M)*mu_M, lm_M = ku_M-2.*mu_M,
		ku_I = (1.-nu_I)/(.5-nu_I)*mu_I, lm_I = ku_I-2.*mu_I,
		ku_L = (1.-nu_L)/(.5-nu_L)*mu_L, lm_L = ku_L-2.*mu_L, kk_L = sqrt(C_L/ku_L), kp_L = sqrt(C_L/mu_L), RR1 = 1./sqrt(c1), hh[3];
	scale_regul(1., kk_L, hh);
	double HHL1 = hh[1], HHL2 = hh[2];
	scale_irreg(1., kk_L, hh);
	double JJL1 = hh[1], JJL2 = hh[2], HHM1, JJM1, HHM2, JJM2;
	scale_regul(RR1, kk_L, hh); HHM1 = hh[1]; HHM2 = hh[2];
	scale_irreg(RR1, kk_L, hh); JJM1 = hh[1]; JJM2 = hh[2];
/////////////////////////////////////////////////////////////////////////
//...дополнительные масштабные функции для сдвигового масштабного модуля;
	scale_regul(1., kp_L, hh);
	double HPL1 = hh[1], HPL2 = hh[2];
	scale_irreg(1., kp_L, hh);
	double JPL1 = hh[1], JPL2 = hh[2], HPM1, JPM1, HPM2, JPM2;
	scale_regul(RR1, kp_L, hh); HPM1 = hh[1]; HPM2 = hh[2];
	scale_irreg(RR1, kp_L, hh); JPM1 = hh[1]; JPM2 = hh[2];
////////////////////////////////////////////////////////////////////////////
//...заполняем матрицу для определения модуля сдвига в поперечной плоскости;
	double matr[16][17] ={
//...равенство функций;
		{ 1., -4.*nu_I, -1., 4.*nu_L, nu_L-1., -2., HHL1/ku_L-2.*HHL2/C_L, 2.*HPL2/C_L, JJL1/ku_L-2.*JJL2/C_L, 2.*JPL2/C_L, 0., 0., 0., 0., 0., 0., 0. },
		{ 1., -2.*(3.-2.*nu_I), -1., 2.*(3.-2.*nu_L), nu_L-.5, 2., 2.*HHL2/C_L, HPL1/mu_L-2.*HPL2/C_L, 2.*JJL2/C_L, JPL1/mu_L-2.*JPL2/C_L, 0., 0., 0., 0., 0., 0., 0. },
//...поверхностные силы;
		{ mu_I, 0., -mu_L, 0., mu_L, 6.*mu_L, -(mu_L/ku_L+.5)*HHL1+(4.*mu_L+ku_L)*HHL2/C_L, HPL1-(4.*mu_L+ku_L)*HPL2/C_L, -(mu_L/ku_L+.5)*JJL1+(4.*mu_L+ku_L)*JJL2/C_L, JPL1-(4.*mu_L+ku_L)*JPL2/C_L, 0., 0., 0., 0., 0., 0., 0. },
		{ mu_I, -6.*mu_I, -mu_L, 6.*mu_L, -mu_L*.5, -6.*mu_L, -(1.-mu_L/ku_L)*HHL1-(7.*mu_L-2.*ku_L)*HHL2/C_L, -1.5*HPL1+(7.*mu_L-2.*ku_L)*HPL2/C_L, -(1.-mu_L/ku_L)*JJL1-(7.*mu_L-2.*ku_L)*JJL2/C_L, -1.5*JPL1+(7.*mu_L-2.*ku_L)*JPL2/C_L, 0., 0., 0., 0., 0., 0., 0. },
#ifdef __NORMAL_DERIV_FIRST0__ //...равенство нулю нормальных производных на границе включения;
		{ 0., 0., 1., -12.*nu_L, nu_L-1., -6., (HHL1-(1.+6.*ku_L/C_L)*HHL2)/ku_L, -2.*HPL1/mu_L+6.*HPL2/C_L, (JJL1-(1.+6.*ku_L/C_L)*JJL2)/ku_L, -2.*JPL1/mu_L+6.*JPL2/C_L, 0., 0., 0., 0., 0., 0., 0. },
		{ 0., 0., 1., -6.*(3.-2.*nu_L), nu_L-.5, 6., -2.*HHL1/ku_L+6.*HHL2/C_L, (HPL1-(1.+6.*mu_L/C_L)*HPL2)/mu_L, -2.*JJL1/ku_L+6.*JJL2/C_L, (JPL1-(1.+6.*mu_L/C_L)*JPL2)/mu_L, 0., 0., 0., 0., 0., 0., 0. },
#else //...равенство нулю когезионного поля на границе включения;
		{ 0., 0., 0., 0., 0., 0., HHL1/ku_L-2.*HHL2/C_L, 2.*HPL2/C_L, JJL1/ku_L-2.*JJL2/C_L, 2.*JPL2/C_L, 0., 0., 0., 0., 0. },
		{ 0., 0., 0., 0., 0., 0., 2.*HHL2/C_L, HPL1/mu_L-2.*HPL2/C_L, 2.*JJL2/C_L, JPL1/mu_L-2.*JPL2/C_L, 0., 0., 0., 0., 0. },
#endif
//...равенство функций на границе слоя;
		{ 0., 0., 1., -4.*nu_L/c1, (1.-nu_L)*c1, 2.*sqr(c1), -HHM1/ku_L+2.*c1*HHM2/C_L, -2.*c1*HPM2/C_L, -JJM1/ku_L+2.*c1*JJM2/C_L, -2.*c1*JPM2/C_L, -1., 4.*nu_M, -1.+nu_M, -2., 0., 0., 0. },
		{ 0., 0., 1., -2.*(3.-2.*nu_L)/c1, (-nu_L+.5)*c1, -2.*sqr(c1), -2.*c1*HHM2/C_L, -HPM1/mu_L+2.*c1*HPM2/C_L, -2.*c1*JJM2/C_L, -JPM1/mu_L+2.*c1*JPM2/C_L, -1., 2.*(3.-2.*nu_M), nu_M-.5, 2., 0., 0., 0. },
//...поверхностные силы на границе слоя;
		{ 0., 0., mu_L, 0., -mu_L*c1, -6.*mu_L*sqr(c1), (mu_L/ku_L+.5)*HHM1-(4.*mu_L+ku_L)*c1*HHM2/C_L, -HPM1+(4.*mu_L+ku_L)*c1*HPM2/C_L, (mu_L/ku_L+.5)*JJM1-(4.*mu_L+ku_L)*c1*JJM2/C_L, -JPM1+(4.*mu_L+ku_L)*c1*JPM2/C_L, -mu_M, 0., mu_M, 6.*mu_M, 0., 0., 0. },
		{ 0., 0., mu_L, -6.*mu_L/c1, mu_L*.5*c1, 6.*mu_L*sqr(c1), (1.-2.*mu_L/ku_L)*HHL1+(7.*mu_L-2.*ku_L)*c1*HHL2/C_L, 1.5*HPL1-(7.*mu_L-2.*ku_L)*c1*HPL2/C_L, (1.-2.*mu_L/ku_L)*JJL1+(7.*mu_L-2.*ku_L)*c1*JJL2/C_L, 1.5*JPL1-(7.*mu_L-2.*ku_L)*c1*JPL2/C_L, -mu_M, 6.*mu_M, -mu_M*.5, -6.*mu_M, 0., 0., 0. },
#ifdef __NORMAL_DERIV_SECOND__ //...равенство нулю нормальных производных на границе слоя;
		{ 0., 0., 1., -12.*nu_L/c1, (nu_L-1.)*c1, -6.*sqr(c1), (HHM1-(1.+6.*c1*ku_L/C_L)*HHM2)/ku_L, -2.*HPM1/mu_L+6.*c1*HPM2/C_L, (JJM1-(1.+6.*c1*ku_L/C_L)*JJM2)/ku_L, -2.*JPM1/mu_L+6.*c1*JPM2/C_L, 0., 0., 0., 0., 0., 0., 0. },
		{ 0., 0., 1., -6.*(3.-2.*nu_L)/c1, (nu_L-.5)*c1, 6.*sqr(c1), -2.*HHM1/ku_L+6.*c1*HHM2/C_L, (HPM1-(1.+6.*c1*mu_L/C_L)*HPM2)/mu_L, -2.*JJM1/ku_L+6.*c1*JJM2/C_L, (JPM1-(1.+6.*c1*mu_L/C_L)*JPM2)/mu_L, 0., 0., 0., 0., 0., 0., 0. },
#else //...равенство нулю когезионного поля на границе слоя;
		{ 0., 0., 0., 0., 0., 0., HHM1/ku_L-2.*c1*HHM2/C_L, 2.*c1*HPM2/C_L, JJM1/ku_L-2.*c1*JJM2/C_L, 2.*c1*JPM2/C_L, 0., 0., 0., 0., 0. },
		{ 0., 0., 0., 0., 0., 0., 2.*c1*HHM2/C_L, HPM1/mu_L-2.*c1*HPM2/C_L, 2.*c1*JJM2/C_L, JPM1/mu_L-2.*c1*JPM2/C_L, 0., 0., 0., 0., 0. },
#endif
//...перемещения на границе эффективной области;
		{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., -4.*nu_M/c0, (1.-nu_M)*c0, 2.*sqr(c0), -2., -1. , 1. },
		{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., -2.*(3.-2.*nu_M)/c0, (.5-nu_M)*c0, -2.*sqr(c0), 2., -1. , 1. },
//...поверхностные силы на границе эффективной области;
		{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., mu_M, 0., -mu_M*c0, -6.*mu_M*sqr(c0), 6., 1., 1. },
		{ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., mu_M, -6.*mu_M*c0, mu_M*.5*c0, 6.*mu_M*sqr(c0), -6., -0.5, 1. },
	};

/////////////////////////////////////////////////////////////////////////////
//...итерационный алгоритм вычисления модуля сдвига в поперечном направлении;
	double optim, sgn0, sgn1, nu_H, mu_H, mu_H0, mu_H1/*, KH = TakeGradLayer_k1(ff, ff_l)*/,
		** matrix = NULL; set_matrix(matrix, 16, 17); int * ii = new_struct<int>(16), k_iter = 0, k, l, m = 16;

/////////////////////////////////////////////
//...цикл по определению сдвиговой жесткости;
	mu_H1 = mu_I*ff+mu_L*ff_l+mu_M*(1.-ff-ff_l); mu_H1 *= 10.;
	nu_H = /*0.5*(1.-mu_H1/KH)*/0.3;
	for (k = 0; k < m; k++)
		for (l = 0; l < m+1; l++) matrix[k][l] = matr[k][l];
	optim = take_system_cyl(matrix, ii, m, mu_H1, nu_H);
	if (optim > 0.) sgn1 = 1.; else sgn1 = -1.;

	mu_H0 = mu_H1;
	do {
		mu_H0 /= 2.; k_iter++;
		nu_H = /*0.5*(1.-mu_H0/KH)*/0.3;
		for (k = 0; k < m; k++)
			for (l = 0; l < m+1; l++) matrix[k][l] = matr[k][l];
		optim = take_system_cyl(matrix, ii, m, mu_H0, nu_H);
	} while (optim*sgn1 > 0. && k_iter < max_iter);
	if (optim > 0.) sgn0 = 1.; else sgn0 = -1.; k_iter = 0;

	do {
		mu_H = (mu_H0+mu_H1)*.5; k_iter++;
		nu_H = /*0.5*(1.-mu_H/KH)*/0.3;
		for (k = 0; k < m; k++)
			for (l = 0; l < m+1; l++)	matrix[k][l] = matr[k][l];
		optim = take_system_cyl(matrix, ii, m, mu_H, nu_H);
		if (optim*sgn0 > 0.) mu_H0 = mu_H; else
			if (optim*sgn0 < 0.) mu_H1 = mu_H; else mu_H0 = mu_H1 = mu_H;
	} while (fabs(mu_H1-mu_H0) > eps && k_iter < max_iter);
	mu_H = (mu_H0+mu_H1)*.5;

	//////////////////////////////////
	//...вычисляем эффективные модули;
	if (sgn0*sgn1 > 0. || fabs(optim) > 1e-6) mu_H = -mu_H;
	for (l = 0; l < m; l++) GG[l] = matrix[l][m];
	delete_struct(matrix); delete_struct(ii);

	return(mu_H);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//...обжатие при растяжении вдоль волокон (второй коэффициент Пуассона) для цилиндрического включения с градиентным слоем;
double CCohes2D::TakeGradLayer_k2(double ff, double ff_l)
{
	double mu_M = get_param(NUM_SHEAR), nu_M = get_param(NUM_SHEAR+1),
		mu_I = get_param(NUM_SHEAR+NUM_SHIFT), nu_I = get_param(NUM_SHEAR+1+NUM_SHIFT),
		mu_L = get_param(NUM_SHEAR+NUM_SHIFT*2), nu_L = get_param(NUM_SHEAR+1+NUM_SHIFT*2), C_L = get_param(NUM_SHEAR-1+NUM_SHIFT*2),
		K_M = mu_M/(1.-2.*nu_M), K_I = mu_I/(1.-2.*nu_I), K_L = mu_L/(1.-2.*nu_L), c0 = ff+ff_l, c1 = ff/c0/*, DD*/,
		ku_M = (1.-nu_M)/(.5-nu_M)*mu_M, lm_M = ku_M-2.*mu_M,
		ku_I = (1.-nu_I)/(.5-nu_I)*mu_I, lm_I = ku_I-2.*mu_I,
		ku_L = (1.-nu_L)/(.5-nu_L)*mu_L, lm_L = ku_L-2.*mu_L, kk_L = sqrt(C_L/ku_L), RR1 = 1./sqrt(c1), hh[3];
	//if (C_L == 0.) { //...классическая четырехтрехфазная модель;
	//	DD = 1./(1.+(K_I-K_L)/ku_L*(1.-c1*(c0+(1.-c0)*(mu_M-mu_L)/ku_M)/(1.+(1.-c0)*(K_L-K_M)/ku_M)));
	//	return(lm_M+DD*c0*(1.+(1.-c1)*(K_I-K_L)/ku_L)/(1.+(1.-c0)*(K_L-K_M)/ku_M)*(lm_L-lm_M+c1*(lm_I-lm_L)/(1.+(1.-c1)*(K_I-K_L)/ku_L)));
	//}
	scale_regul(1., kk_L, hh);
	double HHL1 = hh[1], HHL2 = hh[2];
	scale_irreg(1., kk_L, hh);
	double JJL1 = hh[1], JJL2 = hh[2], HHM1, JJM1;
	scale_regul(RR1, kk_L, hh); HHM1 = hh[1]; 
	scale_irreg(RR1, kk_L, hh); JJM1 = hh[1];
	double matr[8][9] ={
		{ 1., -1., -1.,  HHL1,  JJL1,  0.,  0., 0., 0. },  //...равенство функций;
		{ 1., -2.,  0., -HHL2, -JJL2,  0.,  0., 0., 0. },  //...нормальные производные;
		{ K_I+K_L*(2.-3.*nu_L), -ku_L*1.5, -ku_L*.5, 0., 0., 0., 0., 0., lm_L-lm_I }, //...поверхностные силы;
		{ 0.,  1.,  c1,    0.,    0., -1., -1., 0., 0. },	//...равенство функций на второй границе;
		{ 0.,  0.,  0.,  HHM1,  JJM1,  0.,  0., 0., 0. },	//...равенство нулю когезионного поля;
		{ 0., K_L, -c1*mu_L, 0., 0., -K_M, mu_M,  0., lm_M-lm_L }, //...поверхностные силы на второй границе;
		{ 0.,  0., 0., 0., 0., ku_M,   0., -1., mu_M-lm_M },	//...эффективный модуль;
		{ 0.,  0., 0., 0., 0.,   1.,   c0,  0.,   0. },		//...перемещения на границе эффективной области;
	};

	//////////////////////////////////////////////////////////////////
	//...решаем систему линейных уравнений A0, C0, A1, B1, C1, D1, KH;
	int dim_N = 8, ii[8] ={ 0, 0, 0, 0, 0, 0, 0, 0 }, i, k, l, k0, l0;
	for (i = 0; i < dim_N; i++) {
		double f = 0.;
		///////////////////////////////////////
		//...look for position maximal element;
		for (k = 0; k < dim_N; k++)
			if (ii[k] != 1)
				for (l = 0; l < dim_N; l++)
					if (!ii[l]) {
						if (fabs(matr[k][l]) >= f) f = fabs(matr[k0 = k][l0 = l]);
					}
					else if (ii[l] > 1) return(0.);
					++(ii[l0]);
					///////////////////////////////////////////////////////////
					//...swapping row for diagonal position of maximal element;
					if (k0 != l0)
						for (l = 0; l <= dim_N; l++) {
							f = matr[k0][l]; matr[k0][l] = matr[l0][l]; matr[l0][l] = f;
						}
					if (matr[l0][l0] == 0.) return(0.);
					////////////////////////////////
					//...diagonal row normalization;
					double finv = 1./matr[l0][l0]; matr[l0][l0] = 1.;
					for (l = 0; l <= dim_N; l++) matr[l0][l] *= finv;
					/////////////////////////////////
					//...elimination all outher rows;
					for (k = 0; k < dim_N; k++)
						if (k != l0) {
							finv = matr[k][l0]; matr[k][l0] = 0.;
							for (l = 0; l <= dim_N; l++) matr[k][l] -= matr[l0][l]*finv;
						}
	}
	return(matr[7][8]);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//...модуль сдвига при растяжении вдоль волокон для цилиндрического включения с градиентным слоем;
double CCohes2D::TakeGradLayer_G2(double ff, double ff_l)
{
	double mu_M = get_param(NUM_SHEAR), nu_M = get_param(NUM_SHEAR+1),
		mu_I = get_param(NUM_SHEAR+NUM_SHIFT), nu_I = get_param(NUM_SHEAR+1+NUM_SHIFT),
		mu_L = get_param(NUM_SHEAR+NUM_SHIFT*2), nu_L = get_param(NUM_SHEAR+1+NUM_SHIFT*2), C_L = get_param(NUM_SHEAR-1+NUM_SHIFT*2),
		K_M = mu_M/(1.-2.*nu_M), K_I = mu_I/(1.-2.*nu_I), K_L = mu_L/(1.-2.*nu_L), c0 = ff+ff_l, c1 = ff/c0/*, DD*/,
		ku_M = (1.-nu_M)/(.5-nu_M)*mu_M,
		ku_L = (1.-nu_L)/(.5-nu_L)*mu_L, kk_L = sqrt(C_L/ku_L), RR1 = 1./sqrt(c1), hh[3];
	//if (C_L == 0.) { //...классическая четырехтрехфазная модель;
	//	DD = 1./(1.+(K_I-K_L)/ku_L*(1.-c1*(c0+(1.-c0)*(mu_M-mu_L)/ku_M)/(1.+(1.-c0)*(K_L-K_M)/ku_M)));
	//	return(mu_M+DD*c0*(1.+(1.-c1)*(K_I-K_L)/ku_L)/(1.+(1.-c0)*(K_L-K_M)/ku_M)*(mu_L-mu_M+c1*(mu_I-mu_L)/(1.+(1.-c1)*(K_I-K_L)/ku_L)));
	//}
	scale_regul(1., kk_L, hh);
	double HHL1 = hh[1], HHL2 = hh[2];
	scale_irreg(1., kk_L, hh);
	double JJL1 = hh[1], JJL2 = hh[2], HHM1, JJM1;
	scale_regul(RR1, kk_L, hh); HHM1 = hh[1];
	scale_irreg(RR1, kk_L, hh); JJM1 = hh[1];
	double matr[8][9] ={
		{ 1., -1., -1.,  HHL1,  JJL1,  0.,  0., 0., 0. },  //...равенство функций;
		{ 1., -2.,  0., -HHL2, -JJL2,  0.,  0., 0., 0. },  //...нормальные производные;
		{ K_I+K_L*(2.-3.*nu_L), -ku_L*1.5, -ku_L*.5, 0., 0., 0., 0., 0., mu_L-mu_I }, //...поверхностные силы;
		{ 0.,  1.,  c1,    0.,    0., -1., -1., 0., 0. },	//...равенство функций на второй границе;
		{ 0.,  0.,  0.,  HHM1,  JJM1,  0.,  0., 0., 0. },	//...равенство нулю когезионного поля;
		{ 0., K_L, -c1*mu_L, 0., 0., -K_M, mu_M,  0., mu_M-mu_L }, //...поверхностные силы на второй границе;
		{ 0.,  0., 0., 0., 0., ku_M,   0., -1., 0. }, //...эффективный модуль;
		{ 0.,  0., 0., 0., 0.,   1.,   c0,  0., 0. }, //...перемещения на границе эффективной области;
	};

	//////////////////////////////////////////////////////////////////
	//...решаем систему линейных уравнений A0, C0, A1, B1, C1, D1, KH;
	int dim_N = 8, ii[8] ={ 0, 0, 0, 0, 0, 0, 0, 0 }, i, k, l, k0, l0;
	for (i = 0; i < dim_N; i++) {
		double f = 0.;
		///////////////////////////////////////
		//...look for position maximal element;
		for (k = 0; k < dim_N; k++)
			if (ii[k] != 1)
				for (l = 0; l < dim_N; l++)
					if (!ii[l]) {
						if (fabs(matr[k][l]) >= f) f = fabs(matr[k0 = k][l0 = l]);
					}
					else if (ii[l] > 1) return(0.);
					++(ii[l0]);
					///////////////////////////////////////////////////////////
					//...swapping row for diagonal position of maximal element;
					if (k0 != l0)
						for (l = 0; l <= dim_N; l++) {
							f = matr[k0][l]; matr[k0][l] = matr[l0][l]; matr[l0][l] = f;
						}
					if (matr[l0][l0] == 0.) return(0.);
					////////////////////////////////
					//...diagonal row normalization;
					double finv = 1./matr[l0][l0]; matr[l0][l0] = 1.;
					for (l = 0; l <= dim_N; l++) matr[l0][l] *= finv;
					/////////////////////////////////
					//...elimination all outher rows;
					for (k = 0; k < dim_N; k++)
						if (k != l0) {
							finv = matr[k][l0]; matr[k][l0] = 0.;
							for (l = 0; l <= dim_N; l++) matr[k][l] -= matr[l0][l]*finv;
						}
	}
	return(matr[7][8]);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//...продольный модуль сдвига на основе простого сдвига для цилиндрического включения с градиентным слоем;
double CCohes2D::TakeGradLayer_G2_simple(double ff, double ff_l)
{
	double mu_M = get_param(NUM_SHEAR), nu_M = get_param(NUM_SHEAR+1),
		mu_I = get_param(NUM_SHEAR+NUM_SHIFT), nu_I = get_param(NUM_SHEAR+1+NUM_SHIFT),
		mu_L = get_param(NUM_SHEAR+NUM_SHIFT*2), nu_L = get_param(NUM_SHEAR+1+NUM_SHIFT*2), C_L = get_param(NUM_SHEAR-1+NUM_SHIFT*2),	
		c0 = ff+ff_l, c1 = ff/c0,	kp_L = sqrt(C_L/mu_L), RR1 = 1./sqrt(c1), DD, hh[3];
	if (C_L == 0.) { //...классическая четырехтрехфазная модель;
		DD = ((mu_I-mu_M)-(1.-c1)*(mu_I-mu_L)*(1.+(mu_M-mu_L)/(2.*mu_L)))/(1.+(1.-c1)*(mu_I-mu_L)/(2.*mu_L));
		return(mu_M+c0*DD/(1.+(1.-c0)*DD/(2.*mu_M)));
	}
	scale_regul(1., kp_L, hh);
	double HHL1 = hh[1], HHL2 = hh[2];
	scale_irreg(1., kp_L, hh);
	double JJL1 = hh[1], JJL2 = hh[2], HHM1, JJM1, HHM2, JJM2;
	scale_regul(RR1, kp_L, hh); HHM1 = hh[1]; HHM2 = hh[2];
	scale_irreg(RR1, kp_L, hh); JJM1 = hh[1]; JJM2 = hh[2];
	double matr[8][9] ={
		{ 1., -1., -1.,  HHL1,  JJL1,  0.,  0., 0., 0. },  //...равенство функций;
#ifdef __NORMAL_DERIV_FIRST0__
		{ 1.,  0., -2., -HHL2, -JJL2, 0., 0., 0., 0. },	//...равенство нулю нормальных производных на границе включения;
#else
		{ 0.,  0.,  0.,  HHL1,  JJL1, 0., 0., 0., 0. }, //...равенство нулю когезионного поля на границе включения;
#endif
		{ mu_I, -mu_L,  mu_L, 0., 0., 0.,  0., 0., 0. },     //...поверхностные силы;
		{ 0., 1.,  c1,    0.,    0., -1., -1., 0., 0. },	  //...равенство функций на второй границе;
#ifdef __NORMAL_DERIV_SECOND__
		{ 0., 0., 2.*c1,  HHM2,  JJM2, -1., -1., 0., 0. },	  //...равенство нулю нормальных производных на границе матрицы;
#else
		{ 0., 0.,    0.,  HHM1,  JJM1,  0.,  0., 0., 0. },   //...равенство нулю когезионного поля на границе матрицы;
#endif
		{ 0., mu_L, -c1*mu_L, 0., 0., -mu_M, mu_M, 0., 0. }, //...поверхностные силы на второй границе;
		{ 0.,  0., 0., 0., 0.,    1., c0,       0., 1. },	  //...перемещения на границе эффективной области;
		{ 0.,  0., 0., 0., 0., mu_M, -c0*mu_M, -1., 0. },	  //...эффективный модуль;
	};

	//////////////////////////////////////////////////////////////////
	//...решаем систему линейных уравнений A0, C0, A1, B1, C1, D1, KH;
	int dim_N = 8, ii[8] ={ 0, 0, 0, 0, 0, 0, 0, 0 }, i, k, l, k0, l0;
	for (i = 0; i < dim_N; i++) {
		double f = 0.;
///////////////////////////////////////
//...look for position maximal element;
		for (k = 0; k < dim_N; k++)
			if (ii[k] != 1)
				for (l = 0; l < dim_N; l++)
					if (!ii[l]) {
						if (fabs(matr[k][l]) >= f) f = fabs(matr[k0 = k][l0 = l]);
					}
					else if (ii[l] > 1) return(0.);
					++(ii[l0]);
///////////////////////////////////////////////////////////
//...swapping row for diagonal position of maximal element;
					if (k0 != l0)
						for (l = 0; l <= dim_N; l++) {
							f = matr[k0][l]; matr[k0][l] = matr[l0][l]; matr[l0][l] = f;
						}
					if (matr[l0][l0] == 0.) return(0.);
////////////////////////////////
//...diagonal row normalization;
					double finv = 1./matr[l0][l0]; matr[l0][l0] = 1.;
					for (l = 0; l <= dim_N; l++) matr[l0][l] *= finv;
/////////////////////////////////
//...elimination all outher rows;
					for (k = 0; k < dim_N; k++)
						if (k != l0) {
							finv = matr[k][l0]; matr[k][l0] = 0.;
							for (l = 0; l <= dim_N; l++) matr[k][l] -= matr[l0][l]*finv;
						}
	}
	return(matr[7][8]);
}

//////////////////////////////////////////////////////////////////////////
//...четырехфазная градиентная модель Эщелби для цилиндрических включений;
void CCohes2D::TakeEshelbyModel(double ff, double ff_l, double fK, double fG)
{
	double KH = TakeEshelby_volm(ff, ff_l), GH = TakeEshelby_shear(ff, ff_l), 
			 EH = fabs(GH)*(3*KH-fabs(GH))/KH, nuH = (KH-GH)/(2.*KH), kuH = (1.-nuH)/(.5-nuH)*GH,
			mu_M = get_param(NUM_SHEAR), nu_M = get_param(NUM_SHEAR+1), C_M = get_param(NUM_SHEAR-1), kuM = (1.-nu_M)/(.5-nu_M)*mu_M, K_M = mu_M/(1.-2.*nu_M),
			mu_I = get_param(NUM_SHEAR+NUM_SHIFT), nu_I = get_param(NUM_SHEAR+1+NUM_SHIFT), C_I = get_param(NUM_SHEAR-1+NUM_SHIFT), kuI = (1.-nu_I)/(.5-nu_I)*mu_I, K_I = mu_I/(1.-2.*nu_I),
			mu_L = get_param(NUM_SHEAR+NUM_SHIFT*2), nu_L = get_param(NUM_SHEAR+1+NUM_SHIFT*2), C_L = get_param(NUM_SHEAR-1+NUM_SHIFT*2), kuL = (1.-nu_L)/(.5-nu_L)*mu_L, K_L = mu_L/(1.-2.*nu_L);
	fK *= kuH;
	fG *= GH;

////////////////////////////////////////////////////
//...образуем образец из четырех сферических блоков;
	double rad0 = 1., rad1 = 1./sqrt(ff/(ff+ff_l)), rad2 = 1./sqrt(ff), AA = 4.*rad2;
	GetCircleQuadStruct2(AA, AA, rad0, rad1-rad0, rad2-rad1);
	B[0].type = GRAD_POLY_BLOCK;
	B[2].type = B[3].type = GRAD_ZOOM_BLOCK;

////////////////////////////////////
//...устанавливаем параметры задачи;
	set_mpls(PackInts(3, 3)); //...multipoles degree;
	set_quad(PackInts(4, 2)); //...quadrature degree;
	set_normaliz(1.);			  //...normalization coeffitient;
	set_lagrange(1.);			  //...Lagrange corfficient for LSM;
	if (size_of_param() > NUM_SHEAR+3+NUM_SHIFT*3) {
		param[NUM_SHEAR-1+NUM_SHIFT*3] = param[NUM_SHEAR-1];
		param[NUM_SHEAR+NUM_SHIFT*3] = param[NUM_SHEAR];
		param[NUM_SHEAR+1+NUM_SHIFT*3] = param[NUM_SHEAR+1];
		param[NUM_SHEAR+2+NUM_SHIFT*3] = param[NUM_SHEAR+2];
		param[NUM_SHEAR+3+NUM_SHIFT*3] = param[NUM_SHEAR+3];
		param[NUM_SHEAR-1] = 0.;
		param[NUM_SHEAR] = GH;
		param[NUM_SHEAR+1] = nuH;
		param[NUM_SHEAR+2] = 0.;
		param[NUM_SHEAR+3] = 0.;
	}

//////////////////////////////////
//...определяем блочную структуру;
	solver.set_blocks(N, 3); //<==== number of saved potentials !!!
	solver.n += 12;//<==== number of additional auxilliary arrays!!!
	for (int k = 0; k < solver.N;  k++)
		  solver.set_links(k, B[k].link);

	shapes_init(INITIAL_STATE);
	shapes_init(NULL_STATE);
	LinkPhase2D(MAX_PHASE);

	for (int k = 0; k < solver.N;  k++)
	solver.set_dimension(k, freedom_block(k));
   solver.struct_init();

//////////////////////
//...проверка решения;
//double dd, kk_I = sqrt(C_I/kuI), kp_I = sqrt(C_I/mu_I), kk_L = sqrt(C_L/kuL), kp_L = sqrt(C_L/mu_L),
//		kk_M = sqrt(C_M/kuM), kp_M = sqrt(C_M/mu_M), c0 = ff+ff_l, hh[3];
//		scale_regul(1., kk_I, hh);
//		double HH1 = hh[1], HH2 = hh[2];
//		scale_regul(1., kk_L, hh);
//		double HHL1 = hh[1], HHL2 = hh[2];
//		scale_irreg(1., kk_L, hh);
//		double JJL1 = hh[1], JJL2 = hh[2], HHM1, JJM1, HHM2, JJM2;
//		scale_regul(rad1, kk_L, hh); HHM1 = hh[1]; HHM2 = hh[2];
//		scale_irreg(rad1, kk_L, hh); JJM1 = hh[1]; JJM2 = hh[2];
//		double HHK1, JJK1, HHK2, JJK2;
//		scale_regul(rad1, kk_M, hh); HHK1 = hh[1]; HHK2 = hh[2];
//		scale_irreg(rad1, kk_M, hh); JJK1 = hh[1]; JJK2 = hh[2];
//		double HHP1, JJP1, HHP2, JJP2;
//		scale_regul(rad2, kk_M, hh); HHP1 = hh[1]; HHP2 = hh[2];
//		scale_irreg(rad2, kk_M, hh); JJP1 = hh[1]; JJP2 = hh[2];
///////////////////////////////////////////////////////////////////////////
////...дополнительные масштабные функции для сдвигового масштабного модуля;
//		scale_regul(1., kp_I, hh);
//		double HP1 = hh[1], HP2 = hh[2];
//		scale_regul(1., kp_L, hh);
//		double HPL1 = hh[1], HPL2 = hh[2];
//		scale_irreg(1., kp_L, hh);
//		double JPL1 = hh[1], JPL2 = hh[2], HPM1, JPM1, HPM2, JPM2;
//		scale_regul(rad1, kp_L, hh); HPM1 = hh[1]; HPM2 = hh[2];
//		scale_irreg(rad1, kp_L, hh); JPM1 = hh[1]; JPM2 = hh[2];
//		double HPK1, JPK1, HPK2, JPK2;
//		scale_regul(rad1, kp_M, hh); HPK1 = hh[1]; HPK2 = hh[2];
//		scale_irreg(rad1, kp_M, hh); JPK1 = hh[1]; JPK2 = hh[2];
//		double HPP1, JPP1, HPP2, JPP2;
//		scale_regul(rad2, kp_M, hh); HPP1 = hh[1]; HPP2 = hh[2];
//		scale_irreg(rad2, kp_M, hh); JPP1 = hh[1]; JPP2 = hh[2];
////...перемещения на границе эффективной области;
//		dd = GG[12]-4.*nu_M/c0*GG[13]+(1.-nu_M)*c0*GG[14]+2.*sqr(c0)*GG[15]-(HHP1/kuM-2.*c0*HHP2/C_M)*GG[16]-2.*c0*HPP2/C_M*GG[17]-(JJP1/kuM-2.*c0*JJP2/C_M)*GG[18]-2.*c0*JPP2/C_M*GG[19]-2.*GG[20]-(1.-nuH)*GG[21]-1.;
//		dd = GG[12]-2.*(3.-2.*nu_M)/c0*GG[13]+(.5-nu_M)*c0*GG[14]-2.*sqr(c0)*GG[15]-2.*c0*HHP2/C_M*GG[16]-(HPP1/mu_M-2.*c0*HPP2/C_M)*GG[17]-2.*c0*JJP2/C_M*GG[18]-(JPP1/mu_M-2.*c0*JPP2/C_M)*GG[19]+2.*GG[20]-(.5-nuH)*GG[21]-1.;

/////////////////////////////////////////////////////////////////////////////////////////
//...заносим коэффициенты в представление Папковича-Нейбера, плоское всестороннее сжатие;
	B[0].shape->set_R(1.);
	B[0].shape->set_shape(1, 1., get_param(NUM_SHEAR+3+(-B[0].link[NUM_PHASE]-1)*NUM_SHIFT));
	B[0].shape->set_shape(2, 1., get_param(NUM_SHEAR+2+(-B[0].link[NUM_PHASE]-1)*NUM_SHIFT));
	B[0].shape->A[0][2] = KK[0]*kuI/kuH;
	B[0].shape->A[0][8] = KK[0]*kuI/kuH;
	B[0].shape->A[0][30] = KK[1]*.5*C_I/kuH;
	B[0].shape->A[0][36] = KK[1]*.5*C_I/kuH;
	B[1].shape->set_R(1.);
	B[1].shape->A[0][2] = 1.;
	B[1].shape->A[0][8] = 1.;
	B[2].shape->set_R(1.);
	B[2].shape->set_shape(1, 1., get_param(NUM_SHEAR+3+(-B[2].link[NUM_PHASE]-1)*NUM_SHIFT));
	B[2].shape->set_shape(2, 1., get_param(NUM_SHEAR+2+(-B[2].link[NUM_PHASE]-1)*NUM_SHIFT));
	B[2].shape->set_shape(4, 1., get_param(NUM_SHEAR+3+(-B[2].link[NUM_PHASE]-1)*NUM_SHIFT));
	B[2].shape->set_shape(5, 1., get_param(NUM_SHEAR+2+(-B[2].link[NUM_PHASE]-1)*NUM_SHIFT));
	B[2].shape->A[0][2] = KK[2]*kuL/kuH;
	B[2].shape->A[0][8] = KK[2]*kuL/kuH;
	B[2].shape->A[0][30] = KK[4]*.5*C_L/kuH;
	B[2].shape->A[0][36] = KK[4]*.5*C_L/kuH;
	B[2].shape->A[0][44] = KK[3]*mu_L/kuH;
	B[2].shape->A[0][50] = KK[3]*mu_L/kuH;
	B[2].shape->A[0][72] = KK[5]*.5*C_L/kuH;
	B[2].shape->A[0][78] = KK[5]*.5*C_L/kuH;
	B[3].shape->set_R(1.);
	B[3].shape->set_shape(1, 1., get_param(NUM_SHEAR+3+(-B[3].link[NUM_PHASE]-1)*NUM_SHIFT));
	B[3].shape->set_shape(2, 1., get_param(NUM_SHEAR+2+(-B[3].link[NUM_PHASE]-1)*NUM_SHIFT));
	B[3].shape->set_shape(4, 1., get_param(NUM_SHEAR+3+(-B[3].link[NUM_PHASE]-1)*NUM_SHIFT));
	B[3].shape->set_shape(5, 1., get_param(NUM_SHEAR+2+(-B[3].link[NUM_PHASE]-1)*NUM_SHIFT));
	B[3].shape->A[0][2] = KK[6]*kuM/kuH;
	B[3].shape->A[0][8] = KK[6]*kuM/kuH;
	B[3].shape->A[0][30] = KK[8]*.5*C_M/kuH;
	B[3].shape->A[0][36] = KK[8]*.5*C_M/kuH;
	B[3].shape->A[0][44] = KK[7]*mu_M/kuH*(1.+ff_l/ff);
	B[3].shape->A[0][50] = KK[7]*mu_M/kuH*(1.+ff_l/ff);
	B[3].shape->A[0][72] = KK[9]*.5*C_M/kuH;
	B[3].shape->A[0][78] = KK[9]*.5*C_M/kuH;

/////////////////////////////////
////...деформация чистого сдвига;
	B[0].shape->A[1][2] = GG[0]*kuI/kuH;
	B[0].shape->A[1][8] = -GG[0]*kuI/kuH;
	B[0].shape->A[1][6] = 4.*GG[1]*mu_I*(1.-nu_I)/kuH;
	B[0].shape->A[1][12] = 4.*GG[1]*mu_I*(1.-nu_I)/kuH;
	B[0].shape->A[1][16] = GG[3]*.5*C_I/(mu_I*kuH);
	B[0].shape->A[1][22] = -GG[3]*.5*C_I/(mu_I*kuH);
	B[0].shape->A[1][30] = GG[2]*.5*C_I/(kuI*kuH);
	B[0].shape->A[1][36] = -GG[2]*.5*C_I/(kuI*kuH);
	B[1].shape->A[1][2] = 1.;
	B[1].shape->A[1][8] = -1.;
	B[1].shape->A[1][16] = GG[21]*(.5-nuH)/ff;
	B[1].shape->A[1][22] = -GG[21]*(.5-nuH)/ff;
	B[1].shape->A[1][20] = 2.*GG[20]*(.5-nuH)/(1.5-nuH)/sqr(ff);
	B[1].shape->A[1][26] = 2.*GG[20]*(.5-nuH)/(1.5-nuH)/sqr(ff);
	B[2].shape->A[1][2] = GG[4]*kuL/kuH;
	B[2].shape->A[1][8] = -GG[4]*kuL/kuH;
	B[2].shape->A[1][6] = 4.*GG[5]*mu_L*(1.-nu_L)/kuH;
	B[2].shape->A[1][12] = 4.*GG[5]*mu_L*(1.-nu_L)/kuH;
	B[2].shape->A[1][16] = GG[9]*.5*C_L/(mu_L*kuH);
	B[2].shape->A[1][22] = -GG[9]*.5*C_L/(mu_L*kuH);
	B[2].shape->A[1][30] = GG[8]*.5*C_L/(kuL*kuH);
	B[2].shape->A[1][36] = -GG[8]*.5*C_L/(kuL*kuH);
	B[2].shape->A[1][44] = GG[6]*mu_L*(1.-nu_L)/kuH;
	B[2].shape->A[1][50] = -GG[6]*mu_L*(1.-nu_L)/kuH;
	B[2].shape->A[1][48] = 4.*GG[7]*mu_L*(1.-nu_L)/((3.-2.*nu_L)*kuH);
	B[2].shape->A[1][54] = 4.*GG[7]*mu_L*(1.-nu_L)/((3.-2.*nu_L)*kuH);
	B[2].shape->A[1][58] = GG[11]*.5*C_L/(mu_L*kuH);
	B[2].shape->A[1][64] = -GG[11]*.5*C_L/(mu_L*kuH);
	B[2].shape->A[1][72] = GG[10]*.5*C_L/(kuL*kuH);
	B[2].shape->A[1][78] = -GG[10]*.5*C_L/(kuL*kuH);
	B[3].shape->A[1][2] = GG[12]*kuM/kuH;
	B[3].shape->A[1][8] = -GG[12]*kuM/kuH;
	B[3].shape->A[1][6] = 4.*GG[13]*mu_M*(1.-nu_M)/(kuH*(1.+ff_l/ff));
	B[3].shape->A[1][12] = 4.*GG[13]*mu_M*(1.-nu_M)/(kuH*(1.+ff_l/ff));
	B[3].shape->A[1][16] = GG[17]*.5*C_M/(mu_M*kuH);
	B[3].shape->A[1][22] = -GG[17]*.5*C_M/(mu_M*kuH);
	B[3].shape->A[1][30] = GG[16]*.5*C_M/(kuM*kuH);
	B[3].shape->A[1][36] = -GG[16]*.5*C_M/(kuM*kuH);
	B[3].shape->A[1][44] = GG[14]*mu_M*(1.-nu_M)/kuH*(1.+ff_l/ff);
	B[3].shape->A[1][50] = -GG[14]*mu_M*(1.-nu_M)/kuH*(1.+ff_l/ff);
	B[3].shape->A[1][48] = 4.*GG[15]*mu_M*(1.-nu_M)/((3.-2.*nu_M)*kuH)*sqr(1.+ff_l/ff);
	B[3].shape->A[1][54] = 4.*GG[15]*mu_M*(1.-nu_M)/((3.-2.*nu_M)*kuH)*sqr(1.+ff_l/ff);
	B[3].shape->A[1][58] = GG[19]*.5*C_M/(mu_M*kuH);
	B[3].shape->A[1][64] = -GG[19]*.5*C_M/(mu_M*kuH);
	B[3].shape->A[1][72] = GG[18]*.5*C_M/(kuM*kuH);
	B[3].shape->A[1][78] = -GG[18]*.5*C_M/(kuM*kuH);

//////////////////////////////////
////...комбинированная деформация;
	B[0].shape->A[2][2] = B[0].shape->A[1][2]*fG+B[0].shape->A[0][2]*fK;
	B[0].shape->A[2][8] = B[0].shape->A[1][8]*fG+B[0].shape->A[0][8]*fK;
	B[0].shape->A[2][6] = B[0].shape->A[1][6]*fG;
	B[0].shape->A[2][12] = B[0].shape->A[1][12]*fG;
	B[0].shape->A[2][16] = B[0].shape->A[1][16]*fG;
	B[0].shape->A[2][22] = B[0].shape->A[1][22]*fG;
	B[0].shape->A[2][30] = B[0].shape->A[1][30]*fG+B[0].shape->A[0][30]*fK;
	B[0].shape->A[2][36] = B[0].shape->A[1][36]*fG+B[0].shape->A[0][36]*fK;
	B[1].shape->A[2][2] = B[1].shape->A[1][2]*fG+B[1].shape->A[0][2]*fK;
	B[1].shape->A[2][8] = B[1].shape->A[1][8]*fG+B[1].shape->A[0][8]*fK;
	B[1].shape->A[2][16] = B[1].shape->A[1][16]*fG;
	B[1].shape->A[2][22] = B[1].shape->A[1][22]*fG;
	B[1].shape->A[2][20] = B[1].shape->A[1][20]*fG;
	B[1].shape->A[2][26] = B[1].shape->A[1][26]*fG;
	B[2].shape->A[2][2] = B[2].shape->A[1][2]*fG+B[2].shape->A[0][2]*fK;
	B[2].shape->A[2][8] = B[2].shape->A[1][8]*fG+B[2].shape->A[0][8]*fK;
	B[2].shape->A[2][6] = B[2].shape->A[1][6]*fG;
	B[2].shape->A[2][12] = B[2].shape->A[1][12]*fG;
	B[2].shape->A[2][16] = B[2].shape->A[1][16]*fG;
	B[2].shape->A[2][22] = B[2].shape->A[1][22]*fG;
	B[2].shape->A[2][30] = B[2].shape->A[1][30]*fG+B[2].shape->A[0][30]*fK;
	B[2].shape->A[2][36] = B[2].shape->A[1][36]*fG+B[2].shape->A[0][36]*fK;
	B[2].shape->A[2][44] = B[2].shape->A[1][44]*fG+B[2].shape->A[0][44]*fK;
	B[2].shape->A[2][50] = B[2].shape->A[1][50]*fG+B[2].shape->A[0][50]*fK;
	B[2].shape->A[2][48] = B[2].shape->A[1][48]*fG;
	B[2].shape->A[2][54] = B[2].shape->A[1][54]*fG;
	B[2].shape->A[2][58] = B[2].shape->A[1][58]*fG;
	B[2].shape->A[2][64] = B[2].shape->A[1][64]*fG;
	B[2].shape->A[2][72] = B[2].shape->A[1][72]*fG+B[2].shape->A[0][72]*fK;
	B[2].shape->A[2][78] = B[2].shape->A[1][78]*fG+B[2].shape->A[0][78]*fK;
	B[3].shape->A[2][2] = B[3].shape->A[1][2]*fG+B[3].shape->A[0][2]*fK;
	B[3].shape->A[2][8] = B[3].shape->A[1][8]*fG+B[3].shape->A[0][8]*fK;
	B[3].shape->A[2][6] = B[3].shape->A[1][6]*fG;
	B[3].shape->A[2][12] = B[3].shape->A[1][12]*fG;
	B[3].shape->A[2][16] = B[3].shape->A[1][16]*fG+B[3].shape->A[0][16]*fK;
	B[3].shape->A[2][22] = B[3].shape->A[1][22]*fG+B[3].shape->A[0][22]*fK;
	B[3].shape->A[2][20] = B[3].shape->A[1][20]*fG;
	B[3].shape->A[2][26] = B[3].shape->A[1][26]*fG;

	return;
}

///////////////////////////////////////////////////////////////////////
//...трехфазная градиентная модель Эщелби для цилиндрических включений;
void CCohes2D::TakeEshelbyModel_two(double ff, double fK, double fG)
{
	double det, KH = TakeEshelby_volm_two(ff), GH = TakeEshelby_shear_two(ff, det), 
			 EH = fabs(GH)*(3*KH-fabs(GH))/KH, nuH = (KH-GH)/(2.*KH), kuH = (1.-nuH)/(.5-nuH)*GH,
			mu_M = get_param(NUM_SHEAR), nu_M = get_param(NUM_SHEAR+1), C_M = get_param(NUM_SHEAR-1), kuM = (1.-nu_M)/(.5-nu_M)*mu_M, K_M = mu_M/(1.-2.*nu_M),
			mu_I = get_param(NUM_SHEAR+NUM_SHIFT), nu_I = get_param(NUM_SHEAR+1+NUM_SHIFT), C_I = get_param(NUM_SHEAR-1+NUM_SHIFT), kuI = (1.-nu_I)/(.5-nu_I)*mu_I, K_I = mu_I/(1.-2.*nu_I);
	fK *= kuH;
	fG *= GH;

/////////////////////////////////////////////////
//...образуем образец из трех сферических блоков;
	double rad0 = 1., rad1 = 1./sqrt(ff), AA = 4.*rad1;
	GetCircleQuadStruct(AA, AA, rad0, rad1-rad0);
	B[0].type = GRAD_POLY_BLOCK;
	B[2].type = GRAD_ZOOM_BLOCK;

////////////////////////////////////
//...устанавливаем параметры задачи;
	set_mpls(PackInts(3, 3)); //...multipoles degree;
	set_quad(PackInts(4, 2)); //...quadrature degree;
	set_normaliz(1.);			  //...normalization coeffitient;
	set_lagrange(1.);			  //...Lagrange corfficient for LSM;
	set_geometry(rad0, rad1-rad0);
	if (size_of_param() > NUM_SHEAR+1+NUM_SHIFT*2) {
		param[NUM_SHEAR-1+NUM_SHIFT*2] = param[NUM_SHEAR-1];
		param[NUM_SHEAR+NUM_SHIFT*2] = param[NUM_SHEAR];
		param[NUM_SHEAR+1+NUM_SHIFT*2] = param[NUM_SHEAR+1];
		param[NUM_SHEAR+2+NUM_SHIFT*2] = param[NUM_SHEAR+2];
		param[NUM_SHEAR+3+NUM_SHIFT*2] = param[NUM_SHEAR+3];
		param[NUM_SHEAR-1] = 0.;
		param[NUM_SHEAR] = GH;
		param[NUM_SHEAR+1] = nuH;
		param[NUM_SHEAR+2] = 0.;
		param[NUM_SHEAR+3] = 0.;
	}

//////////////////////////////////
//...определяем блочную структуру;
	solver.set_blocks(N, 3); //<==== number of saved potentials !!!
	solver.n += 12;//<==== number of additional auxilliary arrays!!!
	for (int k = 0; k < solver.N;  k++)
		  solver.set_links(k, B[k].link);

	shapes_init(INITIAL_STATE);
	shapes_init(NULL_STATE);
	LinkPhase2D(MAX_PHASE);

	for (int k = 0; k < solver.N;  k++)
	solver.set_dimension(k, freedom_block(k));
   solver.struct_init();


/////////////////////////////////////////////////////////////////////////////////////////
//...заносим коэффициенты в представление Папковича-Нейбера, плоское всестороннее сжатие;
	B[0].shape->set_R(1.);
	B[0].shape->set_shape(1, 1., get_param(NUM_SHEAR+3+(-B[0].link[NUM_PHASE]-1)*NUM_SHIFT));
	B[0].shape->set_shape(2, 1., get_param(NUM_SHEAR+2+(-B[0].link[NUM_PHASE]-1)*NUM_SHIFT));
	B[0].shape->A[0][2] = KK[0]*kuI/kuH;
	B[0].shape->A[0][8] = KK[0]*kuI/kuH;
	B[0].shape->A[0][30] = KK[1]*.5*C_I/kuH;
	B[0].shape->A[0][36] = KK[1]*.5*C_I/kuH;
	B[1].shape->set_R(1.);
	B[1].shape->A[0][2] = 1.;
	B[1].shape->A[0][8] = 1.;
	B[2].shape->set_R(1.);
	B[2].shape->set_shape(1, 1., get_param(NUM_SHEAR+3+(-B[2].link[NUM_PHASE]-1)*NUM_SHIFT));
	B[2].shape->set_shape(2, 1., get_param(NUM_SHEAR+2+(-B[2].link[NUM_PHASE]-1)*NUM_SHIFT));
	B[2].shape->set_shape(4, 1., get_param(NUM_SHEAR+3+(-B[2].link[NUM_PHASE]-1)*NUM_SHIFT));
	B[2].shape->set_shape(5, 1., get_param(NUM_SHEAR+2+(-B[2].link[NUM_PHASE]-1)*NUM_SHIFT));
	B[2].shape->A[0][2] = KK[2]*kuM/kuH;
	B[2].shape->A[0][8] = KK[2]*kuM/kuH;
	B[2].shape->A[0][30] = KK[4]*.5*C_M/kuH;
	B[2].shape->A[0][36] = KK[4]*.5*C_M/kuH;
	B[2].shape->A[0][44] = KK[3]*mu_M/kuH;
	B[2].shape->A[0][50] = KK[3]*mu_M/kuH;
	B[2].shape->A[0][72] = KK[5]*.5*C_M/kuH;
	B[2].shape->A[0][78] = KK[5]*.5*C_M/kuH;

/////////////////////////////////
////...деформация чистого сдвига;
	B[0].shape->A[1][2] = GG[0]*kuI/kuH;
	B[0].shape->A[1][8] = -GG[0]*kuI/kuH;
	B[0].shape->A[1][6] = 4.*GG[1]*mu_I*(1.-nu_I)/kuH;
	B[0].shape->A[1][12] = 4.*GG[1]*mu_I*(1.-nu_I)/kuH;
	B[0].shape->A[1][16] = GG[3]*.5*C_I/(mu_I*kuH);
	B[0].shape->A[1][22] = -GG[3]*.5*C_I/(mu_I*kuH);
	B[0].shape->A[1][30] = GG[2]*.5*C_I/(kuI*kuH);
	B[0].shape->A[1][36] = -GG[2]*.5*C_I/(kuI*kuH);
	B[1].shape->A[1][2] = 1.;
	B[1].shape->A[1][8] = -1.;
	B[1].shape->A[1][16] = GG[13]*(.5-nuH)/ff;
	B[1].shape->A[1][22] = -GG[13]*(.5-nuH)/ff;
	B[1].shape->A[1][20] = 2.*GG[12]*(.5-nuH)/(1.5-nuH)/sqr(ff);
	B[1].shape->A[1][26] = 2.*GG[12]*(.5-nuH)/(1.5-nuH)/sqr(ff);
	B[2].shape->A[1][2] = GG[4]*kuM/kuH;
	B[2].shape->A[1][8] = -GG[4]*kuM/kuH;
	B[2].shape->A[1][6] = 4.*GG[5]*mu_M*(1.-nu_M)/kuH;
	B[2].shape->A[1][12] = 4.*GG[5]*mu_M*(1.-nu_M)/kuH;
	B[2].shape->A[1][16] = GG[9]*.5*C_M/(mu_M*kuH);
	B[2].shape->A[1][22] = -GG[9]*.5*C_M/(mu_M*kuH);
	B[2].shape->A[1][30] = GG[8]*.5*C_M/(kuM*kuH);
	B[2].shape->A[1][36] = -GG[8]*.5*C_M/(kuM*kuH);
	B[2].shape->A[1][44] = GG[6]*mu_M*(1.-nu_M)/kuH;
	B[2].shape->A[1][50] = -GG[6]*mu_M*(1.-nu_M)/kuH;
	B[2].shape->A[1][48] = 4.*GG[7]*mu_M*(1.-nu_M)/((3.-2.*nu_M)*kuH);
	B[2].shape->A[1][54] = 4.*GG[7]*mu_M*(1.-nu_M)/((3.-2.*nu_M)*kuH);
	B[2].shape->A[1][58] = GG[11]*.5*C_M/(mu_M*kuH);
	B[2].shape->A[1][64] = -GG[11]*.5*C_M/(mu_M*kuH);
	B[2].shape->A[1][72] = GG[10]*.5*C_M/(kuM*kuH);
	B[2].shape->A[1][78] = -GG[10]*.5*C_M/(kuM*kuH);

//////////////////////////////////
////...комбинированная деформация;
	B[0].shape->A[2][2] = B[0].shape->A[1][2]*fG+B[0].shape->A[0][2]*fK;
	B[0].shape->A[2][8] = B[0].shape->A[1][8]*fG+B[0].shape->A[0][8]*fK;
	B[0].shape->A[2][6] = B[0].shape->A[1][6]*fG;
	B[0].shape->A[2][12] = B[0].shape->A[1][12]*fG;
	B[0].shape->A[2][16] = B[0].shape->A[1][16]*fG;
	B[0].shape->A[2][22] = B[0].shape->A[1][22]*fG;
	B[0].shape->A[2][30] = B[0].shape->A[1][30]*fG+B[0].shape->A[0][30]*fK;
	B[0].shape->A[2][36] = B[0].shape->A[1][36]*fG+B[0].shape->A[0][36]*fK;
	B[1].shape->A[2][2] = B[1].shape->A[1][2]*fG+B[1].shape->A[0][2]*fK;
	B[1].shape->A[2][8] = B[1].shape->A[1][8]*fG+B[1].shape->A[0][8]*fK;
	B[1].shape->A[2][16] = B[1].shape->A[1][16]*fG;
	B[1].shape->A[2][22] = B[1].shape->A[1][22]*fG;
	B[1].shape->A[2][20] = B[1].shape->A[1][20]*fG;
	B[1].shape->A[2][26] = B[1].shape->A[1][26]*fG;
	B[2].shape->A[2][2] = B[2].shape->A[1][2]*fG+B[2].shape->A[0][2]*fK;
	B[2].shape->A[2][8] = B[2].shape->A[1][8]*fG+B[2].shape->A[0][8]*fK;
	B[2].shape->A[2][6] = B[2].shape->A[1][6]*fG;
	B[2].shape->A[2][12] = B[2].shape->A[1][12]*fG;
	B[2].shape->A[2][16] = B[2].shape->A[1][16]*fG;
	B[2].shape->A[2][22] = B[2].shape->A[1][22]*fG;
	B[2].shape->A[2][30] = B[2].shape->A[1][30]*fG+B[2].shape->A[0][30]*fK;
	B[2].shape->A[2][36] = B[2].shape->A[1][36]*fG+B[2].shape->A[0][36]*fK;
	B[2].shape->A[2][44] = B[2].shape->A[1][44]*fG+B[2].shape->A[0][44]*fK;
	B[2].shape->A[2][50] = B[2].shape->A[1][50]*fG+B[2].shape->A[0][50]*fK;
	B[2].shape->A[2][48] = B[2].shape->A[1][48]*fG;
	B[2].shape->A[2][54] = B[2].shape->A[1][54]*fG;
	B[2].shape->A[2][58] = B[2].shape->A[1][58]*fG;
	B[2].shape->A[2][64] = B[2].shape->A[1][64]*fG;
	B[2].shape->A[2][72] = B[2].shape->A[1][72]*fG+B[2].shape->A[0][72]*fK;
	B[2].shape->A[2][78] = B[2].shape->A[1][78]*fG+B[2].shape->A[0][78]*fK;

	return;
}

//////////////////////////////////////////////////////////////////////////////////
//...четырехфазная модель Эщелби для цилиндрических включений с градиентным слоем;
void CCohes2D::TakeEshelbyGradModel(double ff, double ff_l, double fK, double fG)
{
	double KH = TakeGradLayer_k1(ff, ff_l), GH = TakeGradLayer_G1(ff, ff_l), 
		   EH = fabs(GH)*(3*KH-fabs(GH))/KH, nuH = (KH-GH)/(2.*KH), kuH = (1.-nuH)/(.5-nuH)*GH,
			mu_M = get_param(NUM_SHEAR), nu_M = get_param(NUM_SHEAR+1), kuM = (1.-nu_M)/(.5-nu_M)*mu_M, K_M = mu_M/(1.-2.*nu_M),
			mu_I = get_param(NUM_SHEAR+NUM_SHIFT), nu_I = get_param(NUM_SHEAR+1+NUM_SHIFT), kuI = (1.-nu_I)/(.5-nu_I)*mu_I, K_I = mu_I/(1.-2.*nu_I), 
			mu_L = get_param(NUM_SHEAR+NUM_SHIFT*2), nu_L = get_param(NUM_SHEAR+1+NUM_SHIFT*2), C_L = get_param(NUM_SHEAR-1+NUM_SHIFT*2), kuL = (1.-nu_L)/(.5-nu_L)*mu_L, K_L = mu_L/(1.-2.*nu_L);
	fK *= kuH;
	fG *= GH;

////////////////////////////////////////////////////
//...образуем образец из четырех сферических блоков;
	double rad0 = 1., rad1 = 1./sqrt(ff/(ff+ff_l)), rad2 = 1./sqrt(ff), AA = 4.*rad2;
	GetCircleQuadStruct2(AA, AA, rad0, rad1-rad0, rad2-rad1);
	B[2].type = GRAD_ZOOM_BLOCK;

////////////////////////////////////
//...устанавливаем параметры задачи;
	set_mpls(PackInts(3, 3)); //...multipoles degree;
	set_quad(PackInts(4, 2)); //...quadrature degree;
	set_normaliz(1.);			  //...normalization coeffitient;
	set_lagrange(1.);			  //...Lagrange corfficient for LSM;
	if (size_of_param() > NUM_SHEAR+3+NUM_SHIFT*3) {
		param[NUM_SHEAR-1+NUM_SHIFT*3] = param[NUM_SHEAR-1];
		param[NUM_SHEAR+NUM_SHIFT*3] = param[NUM_SHEAR];
		param[NUM_SHEAR+1+NUM_SHIFT*3] = param[NUM_SHEAR+1];
		param[NUM_SHEAR+2+NUM_SHIFT*3] = param[NUM_SHEAR+2];
		param[NUM_SHEAR+3+NUM_SHIFT*3] = param[NUM_SHEAR+3];
		param[NUM_SHEAR-1] = 0.;
		param[NUM_SHEAR] = GH;
		param[NUM_SHEAR+1] = nuH;
		param[NUM_SHEAR+2] = 0.;
		param[NUM_SHEAR+3] = 0.;
	}

//////////////////////////////////
//...определяем блочную структуру;
	solver.set_blocks(N, 3); //<==== number of saved potentials !!!
	solver.n += 12;//<==== number of additional auxilliary arrays!!!
	for (int k = 0; k < solver.N;  k++)
		  solver.set_links(k, B[k].link);

	shapes_init(INITIAL_STATE);
	shapes_init(NULL_STATE);
	LinkPhase2D(MAX_PHASE);

	for (int k = 0; k < solver.N;  k++)
	solver.set_dimension(k, freedom_block(k));
   solver.struct_init();

/////////////////////////////////////////////////////////////////////////////////////////
//...заносим коэффициенты в представление Папковича-Нейбера, плоское всестороннее сжатие;
	B[0].shape->set_R(1.);
	B[0].shape->A[0][2] = KK[0]*kuI/kuH;
	B[0].shape->A[0][8] = KK[0]*kuI/kuH;
	B[1].shape->set_R(1.);
	B[1].shape->A[0][2] = 1.;
	B[1].shape->A[0][8] = 1.;
	B[2].shape->set_R(1.);
	B[2].shape->set_shape(1, 1., get_param(NUM_SHEAR+3+(-B[2].link[NUM_PHASE]-1)*NUM_SHIFT));
	B[2].shape->set_shape(2, 1., get_param(NUM_SHEAR+2+(-B[2].link[NUM_PHASE]-1)*NUM_SHIFT));
	B[2].shape->set_shape(4, 1., get_param(NUM_SHEAR+3+(-B[2].link[NUM_PHASE]-1)*NUM_SHIFT));
	B[2].shape->set_shape(5, 1., get_param(NUM_SHEAR+2+(-B[2].link[NUM_PHASE]-1)*NUM_SHIFT));
	B[2].shape->A[0][2] = KK[1]*kuL/kuH;
	B[2].shape->A[0][8] = KK[1]*kuL/kuH;
	B[2].shape->A[0][30] = KK[3]*.5*C_L/kuH;
	B[2].shape->A[0][36] = KK[3]*.5*C_L/kuH;
	B[2].shape->A[0][44] = KK[2]*mu_L/kuH;
	B[2].shape->A[0][50] = KK[2]*mu_L/kuH;
	B[2].shape->A[0][72] = KK[4]*.5*C_L/kuH;
	B[2].shape->A[0][78] = KK[4]*.5*C_L/kuH;
	B[3].shape->set_R(1.);
	B[3].shape->A[0][2] = KK[5]*kuM/kuH;
	B[3].shape->A[0][8] = KK[5]*kuM/kuH;
	B[3].shape->A[0][16] = KK[6]*mu_M/kuH*(1.+ff_l/ff);
	B[3].shape->A[0][22] = KK[6]*mu_M/kuH*(1.+ff_l/ff);

///////////////////////////////
//...деформация чистого сдвига;
	B[0].shape->A[1][2] = GG[0]*kuI/kuH;
	B[0].shape->A[1][8] = -GG[0]*kuI/kuH;
	B[0].shape->A[1][6] = 4.*GG[1]*mu_I*(1.-nu_I)/kuH;
	B[0].shape->A[1][12] = 4.*GG[1]*mu_I*(1.-nu_I)/kuH;
	B[1].shape->A[1][2] = 1.;
	B[1].shape->A[1][8] = -1.;
	B[1].shape->A[1][16] = GG[15]*(.5-nuH)/ff;
	B[1].shape->A[1][22] = -GG[15]*(.5-nuH)/ff;
	B[1].shape->A[1][20] = 2.*GG[14]*(.5-nuH)/(1.5-nuH)/sqr(ff);
	B[1].shape->A[1][26] = 2.*GG[14]*(.5-nuH)/(1.5-nuH)/sqr(ff);
	B[2].shape->A[1][2] = GG[2]*kuL/kuH;
	B[2].shape->A[1][8] = -GG[2]*kuL/kuH;
	B[2].shape->A[1][6] = 4.*GG[3]*mu_L*(1.-nu_L)/kuH;
	B[2].shape->A[1][12] = 4.*GG[3]*mu_L*(1.-nu_L)/kuH;
	B[2].shape->A[1][16] = GG[7]*.5*C_L/(mu_L*kuH);
	B[2].shape->A[1][22] = -GG[7]*.5*C_L/(mu_L*kuH);
	B[2].shape->A[1][30] = GG[6]*.5*C_L/(kuL*kuH);
	B[2].shape->A[1][36] = -GG[6]*.5*C_L/(kuL*kuH);
	B[2].shape->A[1][44] = GG[4]*mu_L*(1.-nu_L)/kuH;
	B[2].shape->A[1][50] = -GG[4]*mu_L*(1.-nu_L)/kuH;
	B[2].shape->A[1][48] = 4.*GG[5]*mu_L*(1.-nu_L)/((3.-2.*nu_L)*kuH);
	B[2].shape->A[1][54] = 4.*GG[5]*mu_L*(1.-nu_L)/((3.-2.*nu_L)*kuH);
	B[2].shape->A[1][58] = GG[9]*.5*C_L/(mu_L*kuH);
	B[2].shape->A[1][64] = -GG[9]*.5*C_L/(mu_L*kuH);
	B[2].shape->A[1][72] = GG[8]*.5*C_L/(kuL*kuH);
	B[2].shape->A[1][78] = -GG[8]*.5*C_L/(kuL*kuH);
	B[3].shape->A[1][2] = GG[10]*kuM/kuH;
	B[3].shape->A[1][8] = -GG[10]*kuM/kuH;
	B[3].shape->A[1][6] = 4.*GG[11]*mu_M*(1.-nu_M)/(kuH*(1.+ff_l/ff));
	B[3].shape->A[1][12] = 4.*GG[11]*mu_M*(1.-nu_M)/(kuH*(1.+ff_l/ff));
	B[3].shape->A[1][16] = GG[12]*mu_M*(1.-nu_M)/kuH*(1.+ff_l/ff);
	B[3].shape->A[1][22] = -GG[12]*mu_M*(1.-nu_M)/kuH*(1.+ff_l/ff);
	B[3].shape->A[1][20] = 4.*GG[13]*mu_M*(1.-nu_M)/((3.-2.*nu_M)*kuH)*sqr(1.+ff_l/ff);
	B[3].shape->A[1][26] = 4.*GG[13]*mu_M*(1.-nu_M)/((3.-2.*nu_M)*kuH)*sqr(1.+ff_l/ff);

////////////////////////////////
//...комбинированная деформация;
	B[0].shape->A[2][2] = B[0].shape->A[1][2]*fG+B[0].shape->A[0][2]*fK;
	B[0].shape->A[2][8] = B[0].shape->A[1][8]*fG+B[0].shape->A[0][8]*fK;
	B[0].shape->A[2][6] = B[0].shape->A[1][6]*fG;
	B[0].shape->A[2][12] = B[0].shape->A[1][12]*fG;
	B[1].shape->A[2][2] = B[1].shape->A[1][2]*fG+B[1].shape->A[0][2]*fK;
	B[1].shape->A[2][8] = B[1].shape->A[1][8]*fG+B[1].shape->A[0][8]*fK;
	B[1].shape->A[2][16] = B[1].shape->A[1][16]*fG;
	B[1].shape->A[2][22] = B[1].shape->A[1][22]*fG;
	B[1].shape->A[2][20] = B[1].shape->A[1][20]*fG;
	B[1].shape->A[2][26] = B[1].shape->A[1][26]*fG;
	B[2].shape->A[2][2] = B[2].shape->A[1][2]*fG+B[2].shape->A[0][2]*fK;
	B[2].shape->A[2][8] = B[2].shape->A[1][8]*fG+B[2].shape->A[0][8]*fK;
	B[2].shape->A[2][6] = B[2].shape->A[1][6]*fG;
	B[2].shape->A[2][12] = B[2].shape->A[1][12]*fG;
	B[2].shape->A[2][16] = B[2].shape->A[1][16]*fG;
	B[2].shape->A[2][22] = B[2].shape->A[1][22]*fG;
	B[2].shape->A[2][30] = B[2].shape->A[1][30]*fG+B[2].shape->A[0][30]*fK;
	B[2].shape->A[2][36] = B[2].shape->A[1][36]*fG+B[2].shape->A[0][36]*fK;
	B[2].shape->A[2][44] = B[2].shape->A[1][44]*fG+B[2].shape->A[0][44]*fK;
	B[2].shape->A[2][50] = B[2].shape->A[1][50]*fG+B[2].shape->A[0][50]*fK;
	B[2].shape->A[2][48] = B[2].shape->A[1][48]*fG;
	B[2].shape->A[2][54] = B[2].shape->A[1][54]*fG;
	B[2].shape->A[2][58] = B[2].shape->A[1][58]*fG;
	B[2].shape->A[2][64] = B[2].shape->A[1][64]*fG;
	B[2].shape->A[2][72] = B[2].shape->A[1][72]*fG+B[2].shape->A[0][72]*fK;
	B[2].shape->A[2][78] = B[2].shape->A[1][78]*fG+B[2].shape->A[0][78]*fK;
	B[3].shape->A[2][2] = B[3].shape->A[1][2]*fG+B[3].shape->A[0][2]*fK;
	B[3].shape->A[2][8] = B[3].shape->A[1][8]*fG+B[3].shape->A[0][8]*fK;
	B[3].shape->A[2][6] = B[3].shape->A[1][6]*fG;
	B[3].shape->A[2][12] = B[3].shape->A[1][12]*fG;
	B[3].shape->A[2][16] = B[3].shape->A[1][16]*fG+B[3].shape->A[0][16]*fK;
	B[3].shape->A[2][22] = B[3].shape->A[1][22]*fG+B[3].shape->A[0][22]*fK;
	B[3].shape->A[2][20] = B[3].shape->A[1][20]*fG;
	B[3].shape->A[2][26] = B[3].shape->A[1][26]*fG;

	return;
}

/////////////////////////////////////////////////////////////////////////////////
//...трехфазная модель Эщелби для цилиндрических включений в градиентной матрице;
void CCohes2D::TakeEshelbyGradModel_two(double ff, double fK, double fG)
{
	double det, KH = TakeGradMatrix_k1(ff), GH = TakeGradMatrix_G1(ff, det), 
			 EH = fabs(GH)*(3*KH-fabs(GH))/KH, nuH = (KH-GH)/(2.*KH), kuH = (1.-nuH)/(.5-nuH)*GH,
			mu_M = get_param(NUM_SHEAR), nu_M = get_param(NUM_SHEAR+1), C_M = get_param(NUM_SHEAR-1), kuM = (1.-nu_M)/(.5-nu_M)*mu_M, K_M = mu_M/(1.-2.*nu_M),
			mu_I = get_param(NUM_SHEAR+NUM_SHIFT), nu_I = get_param(NUM_SHEAR+1+NUM_SHIFT), kuI = (1.-nu_I)/(.5-nu_I)*mu_I, K_I = mu_I/(1.-2.*nu_I);
	fK *= kuH;
	fG *= GH;

/////////////////////////////////////////////////
//...образуем образец из трех сферических блоков;
	double rad0 = 1., rad1 = 1./sqrt(ff), AA = 4.*rad1;
	GetCircleQuadStruct(AA, AA, rad0, rad1-rad0);
	B[2].type = GRAD_ZOOM_BLOCK;

////////////////////////////////////
//...устанавливаем параметры задачи;
	set_mpls(PackInts(3, 3)); //...multipoles degree;
	set_quad(PackInts(4, 2)); //...quadrature degree;
	set_normaliz(1.);			  //...normalization coeffitient;
	set_lagrange(1.);			  //...Lagrange corfficient for LSM;
	set_geometry(rad0, rad1-rad0);
	if (size_of_param() > NUM_SHEAR+1+NUM_SHIFT*2) {
		param[NUM_SHEAR-1+NUM_SHIFT*2] = param[NUM_SHEAR-1];
		param[NUM_SHEAR+NUM_SHIFT*2] = param[NUM_SHEAR];
		param[NUM_SHEAR+1+NUM_SHIFT*2] = param[NUM_SHEAR+1];
		param[NUM_SHEAR+2+NUM_SHIFT*2] = param[NUM_SHEAR+2];
		param[NUM_SHEAR+3+NUM_SHIFT*2] = param[NUM_SHEAR+3];
		param[NUM_SHEAR-1] = 0.;
		param[NUM_SHEAR] = GH;
		param[NUM_SHEAR+1] = nuH;
		param[NUM_SHEAR+2] = 0.;
		param[NUM_SHEAR+3] = 0.;
	}

//////////////////////////////////
//...определяем блочную структуру;
	solver.set_blocks(N, 3); //<==== number of saved potentials !!!
	solver.n += 12;//<==== number of additional auxilliary arrays!!!
	for (int k = 0; k < solver.N;  k++)
		  solver.set_links(k, B[k].link);

	shapes_init(INITIAL_STATE);
	shapes_init(NULL_STATE);
	LinkPhase2D(MAX_PHASE);

	for (int k = 0; k < solver.N;  k++)
	solver.set_dimension(k, freedom_block(k));
   solver.struct_init();

//////////////////////
//...проверка решения;
//double dd, kk_M = sqrt(C_M/kuM), kp_M = sqrt(C_M/mu_M), hh[3];
//		scale_regul(rad0, kk_M, hh);
//		double HHL1 = hh[1], HHL2 = hh[2];
//		scale_irreg(rad0, kk_M, hh);
//		double JJL1 = hh[1], JJL2 = hh[2], HHM1, JJM1, HHM2, JJM2;
//		scale_regul(rad1, kk_M, hh); HHM1 = hh[1]; HHM2 = hh[2];
//		scale_irreg(rad1, kk_M, hh); JJM1 = hh[1]; JJM2 = hh[2];
//		dd = KK[0]-KK[1]-KK[2]+HHL1*KK[3]+JJL1*KK[4];
//		dd = KK[1]-KK[2]-(HHL1+HHL2)*KK[3]-(JJL1+JJL2)*KK[4];
//		//dd = HHL1*KK[3]+JJL1*KK[4];
//		dd = KK[0]*(K_I+mu_M)-KK[1]*kuM-KK[3]*kuM*.5*HHL1-KK[4]*kuM*.5*JJL1;
//		dd = KK[1]+KK[2]*ff-KK[3]*HHM1-KK[4]*JJM1-1.;
//		//dd = 2.*KK[1]+KK[3]*HHM2+KK[4]*JJM2-1.;
//		dd = HHM1*KK[3]+JJM1*KK[4];
//		dd = KK[1]*kuM+KK[3]*kuM*.5*HHM1+KK[4]*kuM*.5*JJM1-KH-mu_M;
///////////////////////////////////////////////////////////////////////////
////...дополнительные масштабные функции для сдвигового масштабного модуля;
//		scale_regul(1., kp_M, hh);
//		double HPL1 = hh[1], HPL2 = hh[2];
//		scale_irreg(1., kp_M, hh);
//		double JPL1 = hh[1], JPL2 = hh[2], HPM1, JPM1, HPM2, JPM2;
//		scale_regul(rad1, kp_M, hh); HPM1 = hh[1]; HPM2 = hh[2];
//		scale_irreg(rad1, kp_M, hh); JPM1 = hh[1]; JPM2 = hh[2];
//		dd = GG[2]-GG[3]*4.*nu_M/ff+GG[4]*(1.-nu_M)*ff+GG[5]*2.*sqr(ff)-GG[6]*(HHM1/kuM-2.*ff*HHM2/C_M)-GG[7]*2.*ff*HPM2/C_M-GG[8]*(JJM1/kuM-2.*ff*JJM2/C_M)-GG[9]*2.*ff*JPM2/C_M-2.*GG[10]-GG[1]-1.;
//		dd = GG[2]-GG[3]*2.*(3.-2.*nu_M)/ff+GG[4]*(.5-nu_M)*ff-GG[5]*2.*sqr(ff)-GG[6]*2.*ff*HHM2/C_M-GG[7]*(HPM1/mu_M-2.*ff*HPM2/C_M)-GG[8]*2.*ff*JJM2/C_M-GG[9]*(JPM1/mu_M-2.*ff*JPM2/C_M)+2.*GG[10]-GG[1]-1.;

/////////////////////////////////////////////////////////////////////////////////////////
//...заносим коэффициенты в представление Папковича-Нейбера, плоское всестороннее сжатие;
	B[0].shape->set_R(1.);
	B[0].shape->A[0][2] = KK[0]*kuI/kuH;
	B[0].shape->A[0][8] = KK[0]*kuI/kuH;
	B[1].shape->set_R(1.);
	B[1].shape->A[0][2] = 1.;
	B[1].shape->A[0][8] = 1.;
	B[2].shape->set_R(1.);
	B[2].shape->set_shape(1, 1., get_param(NUM_SHEAR+3+(-B[2].link[NUM_PHASE]-1)*NUM_SHIFT));
	B[2].shape->set_shape(2, 1., get_param(NUM_SHEAR+2+(-B[2].link[NUM_PHASE]-1)*NUM_SHIFT));
	B[2].shape->set_shape(4, 1., get_param(NUM_SHEAR+3+(-B[2].link[NUM_PHASE]-1)*NUM_SHIFT));
	B[2].shape->set_shape(5, 1., get_param(NUM_SHEAR+2+(-B[2].link[NUM_PHASE]-1)*NUM_SHIFT));
	B[2].shape->A[0][2] = KK[1]*kuM/kuH;
	B[2].shape->A[0][8] = KK[1]*kuM/kuH;
	B[2].shape->A[0][30] = KK[3]*.5*C_M/kuH;
	B[2].shape->A[0][36] = KK[3]*.5*C_M/kuH;
	B[2].shape->A[0][44] = KK[2]*mu_M/kuH;
	B[2].shape->A[0][50] = KK[2]*mu_M/kuH;
	B[2].shape->A[0][72] = KK[4]*.5*C_M/kuH;
	B[2].shape->A[0][78] = KK[4]*.5*C_M/kuH;

/////////////////////////////////
////...деформация чистого сдвига;
	B[0].shape->A[1][2] = GG[0]*kuI/kuH;
	B[0].shape->A[1][8] = -GG[0]*kuI/kuH;
	B[0].shape->A[1][6] = 4.*GG[1]*mu_I*(1.-nu_I)/kuH;
	B[0].shape->A[1][12] = 4.*GG[1]*mu_I*(1.-nu_I)/kuH;
	B[1].shape->A[1][2] = 1.;
	B[1].shape->A[1][8] = -1.;
	B[1].shape->A[1][16] = GG[11]*(.5-nuH)/ff;
	B[1].shape->A[1][22] = -GG[11]*(.5-nuH)/ff;
	B[1].shape->A[1][20] = 2.*GG[10]*(.5-nuH)/(1.5-nuH)/sqr(ff);
	B[1].shape->A[1][26] = 2.*GG[10]*(.5-nuH)/(1.5-nuH)/sqr(ff);
	B[2].shape->A[1][2] = GG[2]*kuM/kuH;
	B[2].shape->A[1][8] = -GG[2]*kuM/kuH;
	B[2].shape->A[1][6] = 4.*GG[3]*mu_M*(1.-nu_M)/kuH;
	B[2].shape->A[1][12] = 4.*GG[3]*mu_M*(1.-nu_M)/kuH;
	B[2].shape->A[1][16] = GG[7]*.5*C_M/(mu_M*kuH);
	B[2].shape->A[1][22] = -GG[7]*.5*C_M/(mu_M*kuH);
	B[2].shape->A[1][30] = GG[6]*.5*C_M/(kuM*kuH);
	B[2].shape->A[1][36] = -GG[6]*.5*C_M/(kuM*kuH);
	B[2].shape->A[1][44] = GG[4]*mu_M*(1.-nu_M)/kuH;
	B[2].shape->A[1][50] = -GG[4]*mu_M*(1.-nu_M)/kuH;
	B[2].shape->A[1][48] = 4.*GG[5]*mu_M*(1.-nu_M)/((3.-2.*nu_M)*kuH);
	B[2].shape->A[1][54] = 4.*GG[5]*mu_M*(1.-nu_M)/((3.-2.*nu_M)*kuH);
	B[2].shape->A[1][58] = GG[9]*.5*C_M/(mu_M*kuH);
	B[2].shape->A[1][64] = -GG[9]*.5*C_M/(mu_M*kuH);
	B[2].shape->A[1][72] = GG[8]*.5*C_M/(kuM*kuH);
	B[2].shape->A[1][78] = -GG[8]*.5*C_M/(kuM*kuH);

//////////////////////////////////
////...комбинированная деформация;
	B[0].shape->A[2][2] = B[0].shape->A[1][2]*fG+B[0].shape->A[0][2]*fK;
	B[0].shape->A[2][8] = B[0].shape->A[1][8]*fG+B[0].shape->A[0][8]*fK;
	B[0].shape->A[2][6] = B[0].shape->A[1][6]*fG;
	B[0].shape->A[2][12] = B[0].shape->A[1][12]*fG;
	B[1].shape->A[2][2] = B[1].shape->A[1][2]*fG+B[1].shape->A[0][2]*fK;
	B[1].shape->A[2][8] = B[1].shape->A[1][8]*fG+B[1].shape->A[0][8]*fK;
	B[1].shape->A[2][16] = B[1].shape->A[1][16]*fG;
	B[1].shape->A[2][22] = B[1].shape->A[1][22]*fG;
	B[1].shape->A[2][20] = B[1].shape->A[1][20]*fG;
	B[1].shape->A[2][26] = B[1].shape->A[1][26]*fG;
	B[2].shape->A[2][2] = B[2].shape->A[1][2]*fG+B[2].shape->A[0][2]*fK;
	B[2].shape->A[2][8] = B[2].shape->A[1][8]*fG+B[2].shape->A[0][8]*fK;
	B[2].shape->A[2][6] = B[2].shape->A[1][6]*fG;
	B[2].shape->A[2][12] = B[2].shape->A[1][12]*fG;
	B[2].shape->A[2][16] = B[2].shape->A[1][16]*fG;
	B[2].shape->A[2][22] = B[2].shape->A[1][22]*fG;
	B[2].shape->A[2][30] = B[2].shape->A[1][30]*fG+B[2].shape->A[0][30]*fK;
	B[2].shape->A[2][36] = B[2].shape->A[1][36]*fG+B[2].shape->A[0][36]*fK;
	B[2].shape->A[2][44] = B[2].shape->A[1][44]*fG+B[2].shape->A[0][44]*fK;
	B[2].shape->A[2][50] = B[2].shape->A[1][50]*fG+B[2].shape->A[0][50]*fK;
	B[2].shape->A[2][48] = B[2].shape->A[1][48]*fG;
	B[2].shape->A[2][54] = B[2].shape->A[1][54]*fG;
	B[2].shape->A[2][58] = B[2].shape->A[1][58]*fG;
	B[2].shape->A[2][64] = B[2].shape->A[1][64]*fG;
	B[2].shape->A[2][72] = B[2].shape->A[1][72]*fG+B[2].shape->A[0][72]*fK;
	B[2].shape->A[2][78] = B[2].shape->A[1][78]*fG+B[2].shape->A[0][78]*fK;

	return;
}

//////////////////////////////////////////////////////////////
//...схема вычисления градиентной модели для конечной трещины;
double CCohes2D::TakeCrack_two(double X, double Y, double kappa, double * AA, double * BB, double * AD, double * BD, int max_iter)
{
	if (X < 0) X = -X; 
	if (AA && BB && AD && BD) {
		complex w = comp(X, Y), z = w*w-1., PP, QQ, PD, QD, FF, FF_dop;

//////////////////////////////////////////////////////
//...вычисляем отображение (прямое и вспомогательное);
		double fi = arg2(2.*X*Y, X*X-Y*Y-1.)*.5, rr = sqrt(sqrt(sqr(X*X-Y*Y-1)+sqr(2.*X*Y))),
			arch_re = log(sqr(X+rr*cos(fi))+sqr(Y+rr*sin(fi))), fi2 = arg2(Y+rr*sin(fi), X+rr*cos(fi)), f_dop, RR,
			f = kappa, sum = 0.;
		FF = comp(rr*cos(fi), rr*sin(fi));

		AA[0] = BD[0] = 1.;
		if (0)
		for (int j = 0; j <= max_iter; j++) {
			PP = QQ = PD = QD = 0.;
////////////////////////////////////////
//...схема Горнера вычисления полиномов;
			PP = (AA[j] /= (2.*(j+1.)))+PP*z;
			for (int k = j-1; k >= 0; k--) {
				QQ = (BB[k+1] = BB[k]/(2.*(k+1.)))+QQ*z;
				PP = (AA[k] = (AA[k]-BB[k+1]-(2.*k+3.)*AA[k+1])/(2.*(k+1.)))+PP*z;
			}
			QQ = (BB[0] = -AA[0])+QQ*z;
			PP = PP*w;

///////////////////////////////////////////////////////
//...схема Горнера вычисления дополнительных полиномов;
			QD = (BD[j] /= (2.*j+1.))+QD*z;
			for (int k = j-1; k >= 0; k--) {
				QD = (BD[k] = (BD[k]-2.*(k+1.)*BD[k+1])/(2.*k+1.))+QD*z;
				PD = (AD[k+1] = (AD[k]-BD[k+1])/(2.*k+3.))+PD*z;
			}
			PD = (AD[0] = -BD[0])+PD*z;
			QD = QD*w;

///////////////////////////////
//...вычисляем вклад в решение;
			sum += f*imag(conj(QD)*PP*FF+QQ*conj(PD*FF)+arch_re*conj(QD)*QQ);
			//sum += imag(f*conj(QD)*QQ*FF);
			f *= kappa;

////////////////////////////////////////
//...следующий шаг вычисления полиномов;
			PP = QQ = PD = QD = 0.;

////////////////////////////////////////
//...схема Горнера вычисления полиномов;
			PP = (AA[j+1] = AA[j]/(2.*j+3.))+PP*z;
			QQ = (BB[j] /= (2.*j+1.))+QQ*z;
			for (int k = j-1; k >= 0; k--) {
				PP = (AA[k+1] = (AA[k]-BB[k+1])/(2.*k+3.))+PP*z;
				QQ = (BB[k] = (BB[k]-2.*(k+1.)*BB[k+1])/(2.*k+1.))+QQ*z;
			}
			PP = (AA[0] = -BB[0])+PP*z;
			QQ = QQ*w;

///////////////////////////////////////////////////////
//...схема Горнера вычисления дополнительных полиномов;
			QD = (BD[j+1] = BD[j]/(2.*(j+1.)))+QD*z;
			PD = (AD[j] = (AD[j]-BD[j+1])/(2.*(j+1.)))+PD*z;
			for (int k = j-1; k >= 0; k--) {
				QD = (BD[k+1] = BD[k]/(2.*(k+1.)))+QD*z;
				PD = (AD[k] = (AD[k]-BD[k+1]-(2.*k+3.)*AD[k+1])/(2.*(k+1.)))+PD*z;
			}
			QD = (BD[0] = -AD[0])+QD*z;
			PD = PD*w;

/////////////////////////////// 
//...вычисляем вклад в решение;
			sum += f*imag(conj(QD)*PP*FF+QQ*conj(PD*FF)+arch_re*conj(QD)*QQ);
			//sum += imag(f*conj(QD)*QQ*FF);
			f *= kappa;
		}

//////////////////////////////////////////////////////
//...реализация решения на основе функции Макдональда;
		double eps = 1e-18;
		fi = arg2(Y, X)*.5; rr = sqrt(RR = sqr(X)+sqr(Y));
		FF_dop = comp(rr*cos(fi), -rr*sin(fi));
		if (rr) FF = comp(1./rr*cos(fi), -1./rr*sin(fi));
		else	  FF = comp(0.);
		sum += imag(FF)-imag(FF_dop); 
		f = f_dop = 1.; RR *= .25;
		int k = 0;
		do {
			k += 1;
			f *= RR/(k*(-.5+k-1.));
			f_dop *= RR/(k*(.5+k-1.));
			sum += imag(FF)*f-imag(FF_dop)*f_dop;
		} 
		while (fabs(f) >= eps && fabs(f_dop) >= eps && k < max_iter);
		return(sum);
	}
	return(0.);
}
#undef __NORMAL_DERIV_FIRST0__
#undef __NORMAL_DERIV_SECOND__
#undef  Message
