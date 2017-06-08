#include "stdafx.h"

#include "cshapes.h"
#include "cporosity2d.h"

#define  Message(Msg)   { printf("%s", Msg);  printf("\n");}

int CPorosity2D::NUM_SHEAR = 7;

//////////////////////////////////
//...initialization of the blocks;
int CPorosity2D::block_shape_init(Block<double> & B, Num_State id_free)
{
	int k, m;
	if (  B.shape && id_free == INITIAL_STATE) delete_shapes(B.shape);
   if (! B.shape && B.mp) {
		B.shape = new CShapeMixer<double>;
		B.shape->add_shape(CreateShape<double>(MP2D_POLY_SHAPE));
		if (get_param(NUM_LOCAL) != 0.)	B.shape->add_shape(CreateShape<double>(SK2D_POLY_SHAPE));

////////////////////////
//...setting parameters;
      B.shape->degree_init1(0, UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));
		B.shape->degree_init1(1, UnPackInts(get_param(NUM_MPLS)), solver.id_norm, 1); //...scalar function of gardient porosity;
      B.shape->set_shape(0, get_param(NUM_MPLS+1)*fabs(B.mp[7])); 
      B.shape->set_shape(1, get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_SHEAR+2), 1., 1.); 

//////////////////////////////////////////////////////////////////
//...local system of coordinates and initializing parametrization;
      B.shape->set_local(B.mp+1);
      B.shape->release  ();
   }

///////////////////////////////////////////////
//...setting cohesion parameter and potentials;
   if (B.shape && id_free != INITIAL_STATE) {
		if (id_free == SPECIAL_STATE) { //...переустановка радиуса и центра мультиполей;
         B.shape->set_local(B.mp+1);
			B.shape->set_shape(0, get_param(NUM_MPLS+1)*fabs(B.mp[7])); 
			B.shape->set_shape(1, get_param(NUM_MPLS+1)*fabs(B.mp[7]), get_param(NUM_SHEAR+2), 1., 1.); 
		}
		else
			if (id_free == OK_STATE && solver.hh[k = (int)(&B-B.B)][0].GetMatrix())
			for (int m = 0; m < solver.id_norm; m++)
				B.shape->set_potential(solver.hh[k][0][m], m);
		else
		if (id_free == NO_STATE && solver.hh[k = (int)(&B-B.B)][0].GetMatrix())
			for (int m = 0; m < solver.id_norm; m++)
				B.shape->get_potential(solver.hh[k][0][solver.id_norm+m], m);
		else
		if (id_free == NULL_STATE) //...переустановка потенциалов (в случае перемены степени, например);
				B.shape->init_potential();
   }                    
   return(B.shape != NULL);
}

//////////////////////////////////////////////
//...realization of common displacements (Rx);
void CPorosity2D::jump1_common_x(double * P, int i, int m)
{
	double G1   = 1./get_param(NUM_SHEAR), beta = get_param(NUM_SHEAR+3)/get_param(NUM_SHEAR-1), 
			alpha = .25/(get_param(NUM_SHEAR+1)-1.)*get_param(NUM_SHEAR+4), * ptr; m += solver.id_norm;
	B[i].shape->cpy_x     (0, B[i].shape->deriv);
	B[i].shape->cpy       (0, B[i].shape->p_cpy);
	B[i].shape->admittance(0, B[i].shape->p_cpy, B[i].shape->deriv, (alpha+1.)*G1, P[0]*alpha*G1);
	B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[1]*alpha*G1, 0.);

	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));

	if (get_param(NUM_LOCAL) != 0.)	{
		B[i].shape->admittance(1, ptr = B[i].shape->FULL(solver.hh[i][0][m], 1), NULL, 0., 0.);
		B[i].shape->adm_x		 (1, ptr, beta);
	}
}

//////////////////////////////////////////////
//...realization of common displacements (Ry);
void CPorosity2D::jump1_common_y(double * P, int i, int m)
{
	double G1   = 1./get_param(NUM_SHEAR), beta = get_param(NUM_SHEAR+3)/get_param(NUM_SHEAR-1), 
			alpha = .25/(get_param(NUM_SHEAR+1)-1.)*get_param(NUM_SHEAR+4), * ptr; m += solver.id_norm;
	B[i].shape->cpy_y     (0, B[i].shape->deriv);
	B[i].shape->cpy       (0, B[i].shape->p_cpy);
	B[i].shape->admittance(0, B[i].shape->p_cpy, B[i].shape->deriv, (alpha+1.)*G1, P[1]*alpha*G1);
	B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[0]*alpha*G1, 0.);

	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));

	if (get_param(NUM_LOCAL) != 0.)	{
		B[i].shape->admittance(1, ptr = B[i].shape->FULL(solver.hh[i][0][m], 1), NULL, 0., 0.);
		B[i].shape->adm_y		 (1, ptr, beta);
	}
}

///////////////////////////////////////////
//...realization of common fracture (Teta);
void CPorosity2D::jump1_common_h(double * P, int i, int m)
{
	double beta = get_param(NUM_SHEAR+3), * ptr; m += solver.id_norm;
	B[i].shape->admittance(0, ptr = B[i].shape->FULL(solver.hh[i][0][m], 0), NULL, 0., 0.);
	B[i].shape->adm_x(0, ptr, beta);
	B[i].shape->admittance(0, ptr = B[i].shape->FULL(solver.hh[i][0][m], 0, 1), NULL, 0., 0.);
	B[i].shape->adm_y(0, ptr, beta);
	B[i].shape->cpy  (1, B[i].shape->FULL(solver.hh[i][0][m], 1));
}

///////////////////////////////////////////////////////////////////////
//...realization of normal derivative of common displacements (dRx/dn);
void CPorosity2D::jump2_common_x(double * P, int i, int m)
{
	double G1   = 1./get_param(NUM_SHEAR), beta = get_param(NUM_SHEAR+3)/get_param(NUM_SHEAR-1), 
			alpha = .25/(get_param(NUM_SHEAR+1)-1.)*get_param(NUM_SHEAR+4), * ptr; m += solver.id_norm;
	B[i].shape->cpy_xx	 (0, B[i].shape->deriv);
	B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[3], 0.);
	B[i].shape->adm_xy    (0, B[i].shape->deriv, P[4]);

	B[i].shape->cpy       (0, B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[0]*alpha*G1, 0.);
	B[i].shape->adm_x     (0, B[i].shape->deriv, P[3]*(2.*alpha+1.)*G1);
	B[i].shape->adm_y     (0, B[i].shape->deriv, P[4]*(1.+alpha)*G1);
	B[i].shape->admittance(0, B[i].shape->p_cpy, NULL, P[1]*alpha*G1, 0.);
	B[i].shape->adm_x     (0, B[i].shape->p_cpy, P[4]*alpha*G1);

	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));

	if (get_param(NUM_LOCAL) != 0.)	{
		B[i].shape->admittance(1, ptr = B[i].shape->FULL(solver.hh[i][0][m], 1), NULL, 0., 0.);
		B[i].shape->adm_xx	 (1, ptr, P[3]*beta);
		B[i].shape->adm_xy	 (1, ptr, P[4]*beta);
	}
}

///////////////////////////////////////////////////////////////////
//...realization of normal derivative common displacements (dRy/dn);
void CPorosity2D::jump2_common_y(double * P, int i, int m)
{
	double G1   = 1./get_param(NUM_SHEAR), beta = get_param(NUM_SHEAR+3)/get_param(NUM_SHEAR-1), 
			alpha = .25/(get_param(NUM_SHEAR+1)-1.)*get_param(NUM_SHEAR+4), * ptr; m += solver.id_norm;
	B[i].shape->cpy_yy	 (0, B[i].shape->deriv);
	B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[4], 0.);
	B[i].shape->adm_xy    (0, B[i].shape->deriv, P[3]);

	B[i].shape->cpy       (0, B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[1]*alpha*G1, 0.);
	B[i].shape->adm_x     (0, B[i].shape->deriv, P[3]*(1.+alpha)*G1);
	B[i].shape->adm_y     (0, B[i].shape->deriv, P[4]*(2.*alpha+1.)*G1);
	B[i].shape->admittance(0, B[i].shape->p_cpy, NULL, P[0]*alpha*G1, 0.);
	B[i].shape->adm_y     (0, B[i].shape->p_cpy, P[3]*alpha*G1);

	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0));

	if (get_param(NUM_LOCAL) != 0.)	{
		B[i].shape->admittance(1, ptr = B[i].shape->FULL(solver.hh[i][0][m], 1), NULL, 0., 0.);
		B[i].shape->adm_xy	 (1, ptr, P[3]*beta);
		B[i].shape->adm_yy	 (1, ptr, P[4]*beta);
	}
}

////////////////////////////////////////////////////////////////////
//...realization of normal derivative of common fracture (dTeta/dn);
void CPorosity2D::jump2_common_h(double * P, int i, int m)
{
	double beta = get_param(NUM_SHEAR+3), * ptr; m += solver.id_norm;
	B[i].shape->admittance (0, ptr = B[i].shape->FULL(solver.hh[i][0][m], 0), NULL, 0., 0.);
	B[i].shape->adm_xx(0, ptr, P[3]*beta);
	B[i].shape->adm_xy(0, ptr, P[4]*beta);

	B[i].shape->admittance (0, ptr = B[i].shape->FULL(solver.hh[i][0][m], 0, 1), NULL, 0., 0.);
	B[i].shape->adm_xy(0, ptr, P[3]*beta);
	B[i].shape->adm_yy(0, ptr, P[4]*beta);

	if (get_param(NUM_LOCAL) != 0.)	{
		B[i].shape->admittance(1, ptr = B[i].shape->FULL(solver.hh[i][0][m], 1), NULL, 0., 0.);
		B[i].shape->adm_x(1, ptr, P[3]);
		B[i].shape->adm_y(1, ptr, P[4]);
	}
}

/////////////////////////////////////////////////////////////////
//...realization of surface forces for common displacements (Px);
void CPorosity2D::jump4_common_x(double * P, int i, int m)
{
	double alpha = .5/(get_param(NUM_SHEAR+1)-1.), beta = get_param(NUM_SHEAR-3)*get_param(NUM_SHEAR+3)+1.,
			 gamma = get_param(NUM_SHEAR-3)*get_param(NUM_SHEAR+3)+alpha*beta, 
			 beta1 = get_param(NUM_SHEAR+3)/get_param(NUM_SHEAR-1)*get_param(NUM_SHEAR)*2.,
			 beta2 = get_param(NUM_SHEAR-3)*(-2.*alpha-1.), * ptr;
	m += solver.id_norm;
	B[i].shape->cpy_xx	 (0, B[i].shape->deriv);
	B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[3], 0.);
	B[i].shape->adm_xy    (0, B[i].shape->deriv, P[4]);
	B[i].shape->admittance(0, B[i].shape->deriv, NULL, gamma, 0.);
	B[i].shape->cpy       (0, B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[0], 0.);
	B[i].shape->adm_x     (0, B[i].shape->deriv, P[3]*beta);
	B[i].shape->adm_y     (0, B[i].shape->deriv, P[4]*beta*(1.+alpha));
	B[i].shape->admittance(0, B[i].shape->p_cpy, NULL, P[1], 0.);
	B[i].shape->adm_y     (0, B[i].shape->p_cpy, P[3]*beta*(-2.*alpha-1.));
	B[i].shape->adm_x     (0, B[i].shape->p_cpy, P[4]*beta*( 1.+alpha));

	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
	
	if (get_param(NUM_LOCAL) != 0.)	{
		B[i].shape->admittance (1, ptr = B[i].shape->FULL(solver.hh[i][0][m], 1), NULL, 0., 0.);
		B[i].shape->adm_xx(1, ptr, P[3]*beta1);
		B[i].shape->adm_xy(1, ptr, P[4]*beta1);
		B[i].shape->adm	(1, ptr, P[3]*beta2);
	}
}

/////////////////////////////////////////////////////////////////
//...realization of surface forces for common displacements (Py);
void CPorosity2D::jump4_common_y(double * P, int i, int m)
{
	double alpha = .5/(get_param(NUM_SHEAR+1)-1.), beta = get_param(NUM_SHEAR-3)*get_param(NUM_SHEAR+3)+1.,
			 gamma = get_param(NUM_SHEAR-3)*get_param(NUM_SHEAR+3)+alpha*beta, 
			 beta1 = get_param(NUM_SHEAR+3)/get_param(NUM_SHEAR-1)*get_param(NUM_SHEAR)*2.,
			 beta2 = get_param(NUM_SHEAR-3)*(-2.*alpha-1.), * ptr;
	m += solver.id_norm;
	B[i].shape->cpy_yy	 (0, B[i].shape->deriv);
	B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[4], 0.);
	B[i].shape->adm_xy    (0, B[i].shape->deriv, P[3]);
	B[i].shape->admittance(0, B[i].shape->deriv, NULL, gamma, 0.);
	B[i].shape->cpy       (0, B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(0, B[i].shape->deriv, NULL, P[1], 0.);
	B[i].shape->adm_x     (0, B[i].shape->deriv, P[3]*beta*(1.+alpha));
	B[i].shape->adm_y     (0, B[i].shape->deriv, P[4]*beta);
	B[i].shape->admittance(0, B[i].shape->p_cpy, NULL, P[0], 0.);
	B[i].shape->adm_x     (0, B[i].shape->p_cpy, P[4]*beta*(-2.*alpha-1));
	B[i].shape->adm_y     (0, B[i].shape->p_cpy, P[3]*beta*( 1.+alpha));

	B[i].shape->cpy(0, B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(0, B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));

	if (get_param(NUM_LOCAL) != 0.)	{
		B[i].shape->admittance (1, ptr = B[i].shape->FULL(solver.hh[i][0][m], 1), NULL, 0., 0.);
		B[i].shape->adm_xy(1, ptr, P[3]*beta1);
		B[i].shape->adm_yy(1, ptr, P[4]*beta1);
		B[i].shape->adm	(1, ptr, P[4]*beta2);
	}
}

//////////////////////////////////////////////////////////////////////////
//...realization of natural condition for common fracture (1/C0*dTeta/dn);
void CPorosity2D::jump4_common_h(double * P, int i, int m)
{
	double alpha = 1./get_param(NUM_SHEAR-1), beta = alpha*get_param(NUM_SHEAR+3), * ptr; m += solver.id_norm;
	B[i].shape->admittance (0, ptr = B[i].shape->FULL(solver.hh[i][0][m], 0), NULL, 0., 0.);
	B[i].shape->adm_xx(0, ptr, P[3]*beta);
	B[i].shape->adm_xy(0, ptr, P[4]*beta);

	B[i].shape->admittance (0, ptr = B[i].shape->FULL(solver.hh[i][0][m], 0, 1), NULL, 0., 0.);
	B[i].shape->adm_xy(0, ptr, P[3]*beta);
	B[i].shape->adm_yy(0, ptr, P[4]*beta);

	if (get_param(NUM_LOCAL) != 0.)	{
		B[i].shape->admittance(1, ptr = B[i].shape->FULL(solver.hh[i][0][m], 1), NULL, 0., 0.);
		B[i].shape->adm_x(1, ptr, P[3]*alpha);
		B[i].shape->adm_y(1, ptr, P[4]*alpha);
	}
}

//////////////////////////////////////////////
//...transformation of the collocation vector;
void CPorosity2D::jump_make_local(int i, int m)
{
	if (! solver.dim || ! solver.hh[i][0].GetMatrix()) return;
	m  += solver.id_norm;
	for (int j = 0; j < solver.dim[i]; j++) {
		double P[] = { solver.hh[i][0][m  ][j], 
							solver.hh[i][0][m+1][j], 0.};
		B[i].shape->norm_local(P);
		solver.hh[i][0][m  ][j] = P[0];
		solver.hh[i][0][m+1][j] = P[1];
	}
}

void CPorosity2D::jump_make_common(int i, int m)
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

///////////////////////////////////////////////////////////////////////////
//...inclusion of the boundary condition data to the solver for all blocks;
Num_State CPorosity2D::gram1(CGrid * nd, int i, int id_local)
{
	if (nd && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double G1 = get_param(NUM_SHEAR), hx, hy, p3, f, P[6];
		int 	 m  = solver.id_norm;

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
			B[i].shape->parametrization_hess(P);

			if (p3 == MIN_HIT || p3 == NUMS_BND || p3 == MAX_HIT) {
				jump1_common_x(P, i, 0); solver.admittance(i, 0, G1);
				jump1_common_y(P, i, 1); solver.admittance(i, 1, G1);
				jump1_common_h(P, i, 2); solver.admittance(i, 2, G1);
			}
			else
			if (0. <= p3 && p3 <= (double)(NUMS_BND-SPECIAL_BND)) {
				jump4_common_x(P, i, 0); 
				jump4_common_y(P, i, 1); 
				jump4_common_h(P, i, 2); 
			}
			jump_make_common(i, 0);

////////////////////////////
//...composition functional;
			solver.to_equationDD(i, solver.hh[i][0][m],	 solver.hh[i][0][m],	  f);
			solver.to_equationDD(i, solver.hh[i][0][m+1], solver.hh[i][0][m+1], f);
			solver.to_equationDD(i, solver.hh[i][0][m+2], solver.hh[i][0][m+2], f);
		
			if (p3 == MIN_HIT || p3 == (double)(NORMS_BND-SPECIAL_BND) || p3 == NUMS_BND || p3 == MAX_HIT) f *= G1;
			if (fabs(hx) > EE)
				solver.to_equationHH(i, 0, solver.hh[i][0][m], hx*f);
			if (fabs(hy) > EE)
				solver.to_equationHH(i, 0, solver.hh[i][0][m+1], hy*f);
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

////////////////////////////////////////////////////////////////
//...inclusion of the joining data to the solver for all blocks;
Num_State CPorosity2D::transfer1(CGrid * nd, int i, int k, int id_local)
{
	if (nd && 0 <= i && i < solver.N) {
      double G1 = get_param(NUM_SHEAR), f, P[6], 
				 f1 = 1., g1 = G1*.5, g2 = G1*.5, g0 = -G1;
      int m = solver.id_norm;

////////////////////////////////////////////////
//...switch on isolated norm on stitching sides;
		int id_isolated = 0;
		if (id_isolated) {
			g0 = g1 = G1;
			f1 = g2 = 0.;
		}

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
				B[i].shape->parametrization_hess(P);

				jump1_common_x(P, i, 0); 
				jump2_common_x(P, i, 3); 
				solver.admittance(i, 0, g1, 3, g2); 
				solver.admittance(i, 3, g0, 0, f1); 

				jump1_common_y(P, i, 1); 
				jump2_common_y(P, i, 4); 
				solver.admittance(i, 1, g1, 4, g2); 
				solver.admittance(i, 4, g0, 1, f1); 

				jump1_common_h(P, i, 2); 
				jump2_common_h(P, i, 5); 
				solver.admittance(i, 2, G1);
				solver.admittance(i, 5, G1);

				jump_make_common(i, 0);
				jump_make_common(i, 3);

////////////////////////////////////////////////////////////////////
//...jump of all neaded displacement moments for neighbouring block;
				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);
				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);
				B[k].shape->parametrization_hess(P);

				jump1_common_x(P, k, 0); 
				jump2_common_x(P, k, 3); 
				solver.admittance(k, 0, g1, 3, g2); 
				solver.admittance(k, 3, g0, 0, f1); 

				jump1_common_y(P, k, 1); 
				jump2_common_y(P, k, 4); 
				solver.admittance(k, 1, g1, 4, g2); 
				solver.admittance(k, 4, g0, 1, f1); 

				jump1_common_h(P, k, 2); 
				jump2_common_h(P, k, 5); 
				solver.admittance(k, 2, G1);
				solver.admittance(k, 5, G1);

				jump_make_common(k, 0);
				jump_make_common(k, 3);

////////////////////////////
//...composition functional;
				solver.to_transferTR(i, j, solver.hh[i][0][m],   solver.hh[k][0][m],   f);
				solver.to_transferDD(i, j, solver.hh[i][0][m],   solver.hh[k][0][m+3], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+4], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+5], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

///////////////////////////////////////////////////////////////////////////
//...inclusion of the boundary condition data to the solver for all blocks;
Num_State CPorosity2D::gram4(CGrid * nd, int i, int id_local)
{
	if (nd && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double G1 = get_param(NUM_SHEAR), hx, hy, p3, f, P[6];
		int 	 m  = solver.id_norm;

////////////////////////////////////////////////////////////////////
//...realization of pattern cell boundary condition for Gram matrix;
		nd->TestGrid("nodes.bln", 0.0005, 0., 0., 0., AXIS_Z, 1);
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
			B[i].shape->parametrization_hess(P);

			jump1_common_x(P, i, 0); 
			jump1_common_y(P, i, 1); 

			jump4_common_x(P, i, 3); 
			jump4_common_y(P, i, 4); 

			if (1 || get_param(NUM_LOCAL)) {//...using gradient defectness;
				jump1_common_h(P, i, 2); 
				jump4_common_h(P, i, 5); 
			}
			jump_make_common(i, 0);
			jump_make_common(i, 3);

////////////////////////////////////
//...composition collocation vector;
			if (p3 == MIN_HIT || p3 == NUMS_BND || p3 == MAX_HIT) {
				solver.admittance(i, 6, 0., 0, G1);
				solver.admittance(i, 7, 0., 1, G1);
				solver.admittance(i, 8, 0., 2, G1);
			}

////////////////////////////////////////////////////
//...граничные условия методом наименьших квадратов;
			if (p3 == MIN_HIT || p3 == NUMS_BND || p3 == MAX_HIT) {
				solver.to_equationDD(i, solver.hh[i][0][m+6], solver.hh[i][0][m+6], f);
				solver.to_equationDD(i, solver.hh[i][0][m+7], solver.hh[i][0][m+7], f);
				solver.to_equationDD(i, solver.hh[i][0][m+8], solver.hh[i][0][m+8], f);

				if (fabs(hx) > EE)
					solver.to_equationHH(i, 0, solver.hh[i][0][m+6], hx*G1*f);
				if (fabs(hy) > EE)
					solver.to_equationHH(i, 0, solver.hh[i][0][m+7], hy*G1*f);
			}

/////////////////////////////////////////////////
//...энергетические слагаемые (нормированные G1);
			solver.to_equationEE(i, solver.hh[i][0][m],	 solver.hh[i][0][m+3], -(f *= G1));
			solver.to_equationEE(i, solver.hh[i][0][m+1], solver.hh[i][0][m+4], -f);
			solver.to_equationEE(i, solver.hh[i][0][m+2], solver.hh[i][0][m+5], -f);

			if (1. <= p3 && p3 <= (double)(NUMS_BND-SPECIAL_BND)) {
				if (fabs(hx) > EE)
					solver.to_equationEH(i, 0, solver.hh[i][0][m], -hx*f);
				if (fabs(hy) > EE)
					solver.to_equationEH(i, 0, solver.hh[i][0][m+1], -hy*f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

///////////////////////////////////////////////////////////////
//...формирование матриц перехода с учетом функционала энергии;
Num_State CPorosity2D::transfer4(CGrid * nd, int i, int k, int id_local)
{
	if (nd && 0 <= i && i < solver.N) {
      double G1 = get_param(NUM_SHEAR), f, P[6];
      int m = solver.id_norm;

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
				B[i].shape->parametrization_hess(P);

				jump1_common_x(P, i, 0); solver.admittance(i, 0, G1); 
				jump4_common_x(P, i, 3); 

				jump1_common_y(P, i, 1); solver.admittance(i, 1, G1); 
				jump4_common_y(P, i, 4); 

				if (get_param(NUM_LOCAL)) {//...using gradient defectness;
					jump1_common_h(P, i, 2); solver.admittance(i, 2, G1); 
					jump4_common_h(P, i, 5); 
				}
				jump_make_common(i, 0);
				jump_make_common(i, 3);

////////////////////////////////////////////////////////////////////
//...jump of all neaded displacement moments for neighbouring block;
				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);
				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);
				B[k].shape->parametrization_hess(P);

				jump1_common_x(P, k, 0); solver.admittance(k, 0, G1); 
				jump4_common_x(P, k, 3); 

				jump1_common_y(P, k, 1); solver.admittance(k, 1, G1); 
				jump4_common_y(P, k, 4); 

				if (1 || get_param(NUM_LOCAL)) {//...using gradient defectness;
					jump1_common_h(P, k, 2); solver.admittance(k, 2, G1); 
					jump4_common_h(P, k, 5); 
				}
				jump_make_common(k, 0);
				jump_make_common(k, 3);

/////////////////////////////////////////////////
//...сшивка функций методом наименьших квадратов;
				solver.to_transferTR(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);

/////////////////////////////////////////////////
//...энергетические слагаемые (нормированные G1);
				solver.to_equationER(i, solver.hh[i][0][m], solver.hh[i][0][m+3], -f);
				solver.to_equationEL(k, solver.hh[k][0][m], solver.hh[k][0][m+3],  f);

				solver.to_equationER(i, solver.hh[i][0][m+1], solver.hh[i][0][m+4], -f);
				solver.to_equationEL(k, solver.hh[k][0][m+1], solver.hh[k][0][m+4],  f);

				solver.to_equationER(i, solver.hh[i][0][m+2], solver.hh[i][0][m+5], -f);
				solver.to_equationEL(k, solver.hh[k][0][m+2], solver.hh[k][0][m+5],  f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////
//...auxiliary function for parameters of doubly connected media;
void CPorosity2D::set_fasa_hmg(double nju, double G1, double C0, double alpha, double gamma)
{
	if (size_of_param() > NUM_SHEAR+1) {
		param[NUM_SHEAR-3] = alpha;
		param[NUM_SHEAR-2] = gamma;
		param[NUM_SHEAR-1] = C0;
		param[NUM_SHEAR]	 = G1;
		param[NUM_SHEAR+1] = nju;
		param[NUM_SHEAR+2] = (.5-nju)*alpha/((1.-nju)*gamma*G1);
		param[NUM_SHEAR+3] = param[NUM_SHEAR+2]/(1.-param[NUM_SHEAR+2]*alpha);
		param[NUM_SHEAR+4] = 1.-param[NUM_SHEAR+3]*(1.-2.*nju)*alpha;
		param[NUM_SHEAR+2] = sqrt(C0*gamma*(1.-param[NUM_SHEAR+2]*alpha));
	}
}

//////////////////////////////////////////////////////////
//...counting header for solving plane elasticity problem;
Num_State CPorosity2D::computing_header(Num_Comput Num)
{
	int k, n_rhs = 2;
	char msg[201];

//	solver.set_mode(REDUCED_MESSAGE); //...отключаем подробную печать;
	if (! solver.mode(NO_MESSAGE)) {
		Message(" ");
		sprintf(msg, "CPorosity2D sample: N_sm = %d, N_mpl = %d, G1 = %g, C0 = %g", N,
				  UnPackInts(get_param(NUM_MPLS)), get_param(NUM_SHEAR), get_param(NUM_SHEAR-1));
		Message(msg);

		Message(" ");
		Message("Junction counting...");

		switch (Num){
			case   BASIC_COMPUT: Message("Analytical Blocks..."); break;
			case MAPPING_COMPUT: Message("FEM Blocks...");			break;
		}
		Message(" ");
	}

///////////////////////////////////
//...устанавливаем параметры среды;
	set_fasa_hmg(param[NUM_SHEAR+1], param[NUM_SHEAR], param[NUM_SHEAR-1], param[NUM_SHEAR-3], param[NUM_SHEAR-2]);
	solver.set_blocks(N, n_rhs); //<==== number of saved potentials !!!
	solver.n += 9;//<==== number of additional auxilliary arrays!!!
	for (k = 0; k < solver.N;  k++)
		  solver.set_links(k, B[k].link);

	shapes_init(INITIAL_STATE);
	shapes_init(NULL_STATE);

/////////////////////////////////////////////////////////
//...делаем перенумерацию структуры и задаем размерность;
	if (! solver.struct_permutat(solver.id_change == EXTERN_STATE ? NULL_STATE : OK_STATE) || 
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
void CPorosity2D::GetFuncAllValues(double X, double Y, double Z, double * F, int i, Num_Value id_F, int id_variant, int iparam)
{
	if (! F) return;
	double P[6]  = { X, Y, Z, 1., 0., 0.};

/////////////////////////////////////
//...operation with all input points;
	if (0 <= i && i < N && B[i].shape && B[i].mp) {
		int m = solver.id_norm;

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
				B[i].shape->parametrization_hess(P);
				jump1_common_x(P, i, 0); 
				jump1_common_y(P, i, 1); 

/////////////////////////////////////
//...calculation common displacement;
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F);

				if ((solv%ENERGY_SOLVING) && 0) F[solv%ENERGY_SOLVING-1] -= X;
			}  break;
			case FRACT_VALUE: {
/////////////////////
//...common fracture;
				B[i].shape->parametrization_hess(P);
				jump1_common_h(P, i, 0); 

/////////////////////////////////
//...calculation common fracture;
				F[0] = B[i].shape->potential(solver.hh[i][0][m], id_variant);
			}  break;
			case PUREFRACT_VALUE: {
////////////////////////////////////
//...potential psi (gradient field);
				B[i].shape->parametrization(1, P, 2);//...calculation with hessian;
				B[i].shape->cpy(1, B[i].shape->FULL(solver.hh[i][0][m], 1));

/////////////////////
//...calculation psi;
				F[0] = B[i].shape->potential(solver.hh[i][0][m], id_variant);
			}     break;
			case STRESS_X_VALUE: {
				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
///////////////////
//...common forces;
				B[i].shape->parametrization_hess(P);
				jump4_common_x(P, i, 0); 
				jump4_common_y(P, i, 1); 

/////////////////////////////////////////////
//...calculation stress tensor (txx and txy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F);
			}  break;
			case STRESS_Y_VALUE: {
				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
///////////////////
//...common forces;
				B[i].shape->parametrization_hess(P);
				jump4_common_x(P, i, 0); 
				jump4_common_y(P, i, 1); 

/////////////////////////////////////////////
//...calculation stress tensor (txy and tyy);
				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
				B[i].shape->norm_common(F);
        }   break;
        case DILAT_VALUE: { //...common dilatation;
				double alpha = (.5/(get_param(NUM_SHEAR+1)-1)+1.)/get_param(NUM_SHEAR), beta = get_param(NUM_SHEAR+3);
///////////////////////////////////////
//...dilatation of classical potential;
				B[i].shape->parametrization_hess(P);
				B[i].shape->cpy_x(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 0));
				B[i].shape->cpy_y(0, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
/////////////////////
//...common fracture;
				jump1_common_h(P, i, 1); 

////////////////////////////////////////
//...common dilatation of displacements;
             F[0] = (B[i].shape->potential(solver.hh[i][0][m],   id_variant)+
							B[i].shape->potential(solver.hh[i][0][m+1], id_variant))*alpha;
        }     break;
        default : F[0] = i; F[1] = 0.;
     }
  }
}
#undef  Message
