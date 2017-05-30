#include "stdafx.h"

#include "cshapes.h"
#include "cvisco2d.h"

#define  Message(Msg)   { printf("%s", Msg);  printf("\n");}


int CVisco2D::NUM_SHIFT = 4;
int CVisco2D::NUM_SHEAR = 6;
int CVisco2D::NUM_BASIC = 8;

//////////////////////////////////
//...initialization of the blocks;
int  CVisco2D::block_shape_init(Block<complex> & B, Num_State id_free)
{
	int k, m;
   if (  B.shape && id_free == INITIAL_STATE) delete_shapes(B.shape);
   if (! B.shape && B.mp) {
		B.shape = new CShapeMixer<complex>;
		if ((B.type & ERR_CODE) == ZOOM_BLOCK && B.mp[0] == ID_MAP(1, SPHEROID_GENUS)) {
			B.shape->add_shape(CreateShape<complex>(MP2D_POLY_SHAPE));
			B.shape->add_shape(CreateShape<complex>(MP2D_ZOOM_SHAPE));
////////////////////////
//...setting parameters;
			B.shape->degree_init1(0, UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));
			B.shape->set_shape(0, get_param(NUM_MPLS+1)*fabs(B.mp[7]));

			B.shape->degree_init1(1, UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));
			B.shape->set_shape(1, fabs(B.mp[8]));
		}
		else {
			B.shape->add_shape(CreateShape<complex>(MP2D_POLY_SHAPE));

////////////////////////
//...setting parameters;
			B.shape->set_shape(get_param(NUM_MPLS+1)*fabs(B.mp[7]));
			B.shape->degree_init1(UnPackInts(get_param(NUM_MPLS)), solver.id_norm, draft_dim(type()));

			if (B.link[NUM_PHASE] == -2) //...another degree of multipoles for inclusion!!!
			B.shape->degree_init1(UnPackInts(get_param(NUM_MPLS), 1), solver.id_norm,  draft_dim(type()));
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
			if ((B.type & ERR_CODE) == ZOOM_BLOCK && B.mp[0] == ID_MAP(1, SPHEROID_GENUS)) {
				B.shape->set_shape(0, get_param(NUM_MPLS+1)*fabs(B.mp[7]));
				B.shape->set_shape(1, fabs(B.mp[8]));
			}
			else B.shape->set_shape(get_param(NUM_MPLS+1)*fabs(B.mp[7]));		}
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

/////////////////////////////////////////////////
//...realization of classical displacements (Ux);
void CVisco2D::jump1_classic_x(double * P, int i, int m)
{
  int  shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
  complex G0 = comp(get_param(NUM_BASIC+2+shift), get_param(NUM_BASIC+3+shift)), G1 =  0.5/G0, 
			 K0 = comp(get_param(NUM_BASIC+0+shift), get_param(NUM_BASIC+1+shift)), K1 = G1/(3.*K0+4.*G0),  
			 A0 = (3.*K0+7.*G0)*K1, A1 = (3.*K0+G0)*K1; m += solver.id_norm;
	B[i].shape->cpy_x     (B[i].shape->deriv);
	B[i].shape->cpy       (B[i].shape->p_cpy);
	//B[i].shape->admittance(B[i].shape->p_cpy, B[i].shape->deriv, A0, -(complex)P[0]*A1);
	//B[i].shape->admittance(B[i].shape->deriv, NULL, -P[1]*A1*G1, 0.);

	B[i].shape->cpy(B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
}

/////////////////////////////////////////////////
//...realization of classical displacements (Uy);
void CVisco2D::jump1_classic_y(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
	complex G1   = 1./get_param(NUM_SHEAR+shift), 
			alpha  = .25/(1.-get_param(NUM_SHEAR+1+shift)); m += solver.id_norm;
	B[i].shape->cpy_y     (B[i].shape->deriv);
	B[i].shape->cpy       (B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->p_cpy, B[i].shape->deriv, (1.-alpha)*G1, -P[1]*alpha*G1);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[0]*alpha*G1, 0.);

	B[i].shape->cpy(B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
}

///////////////////////////////////////////////////////////////////
//...realization of normal derivative classical displacements (Ux);
void CVisco2D::jump2_classic_x(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
	complex G1   = 1./get_param(NUM_SHEAR+shift), 
			alpha  = .25/(1.-get_param(NUM_SHEAR+1+shift)); m += solver.id_norm;
	B[i].shape->cpy_x		();
	B[i].shape->admittance(B[i].shape->deriv, NULL, 0., 0.);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, 0., 0.);
	B[i].shape->deriv_X	(B[i].shape->deriv, P[3]);
	B[i].shape->deriv_Y	(B[i].shape->p_cpy, P[4]);
	B[i].shape->admittance(B[i].shape->deriv, B[i].shape->p_cpy, 1., 1.);
	B[i].shape->cpy_x		();

	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[0]*alpha*G1, 0.);
	B[i].shape->adm_x     (B[i].shape->deriv, P[3]*(1.-2.*alpha)*G1);
	B[i].shape->adm_y     (B[i].shape->deriv, P[4]*(1.-alpha)*G1);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, -P[1]*alpha*G1, 0.);
	B[i].shape->adm_x     (B[i].shape->p_cpy, -P[4]*alpha*G1);

	B[i].shape->cpy(B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
}

///////////////////////////////////////////////////////////////////
//...realization of normal derivative classical displacements (Uy);
void CVisco2D::jump2_classic_y(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
	complex G1   = 1./get_param(NUM_SHEAR+shift), 
			alpha  = .25/(1.-get_param(NUM_SHEAR+1+shift)); m += solver.id_norm;
	B[i].shape->cpy_y		();
	B[i].shape->admittance(B[i].shape->deriv, NULL, 0., 0.);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, 0., 0.);
	B[i].shape->deriv_X	(B[i].shape->deriv, P[3]);
	B[i].shape->deriv_Y	(B[i].shape->p_cpy, P[4]);
	B[i].shape->admittance(B[i].shape->deriv, B[i].shape->p_cpy, 1., 1.);
	B[i].shape->cpy_y		();

	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[0]*alpha*G1, 0.);
	B[i].shape->adm_y     (B[i].shape->deriv, -P[3]*alpha*G1);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, -P[1]*alpha*G1, 0.);
	B[i].shape->adm_x     (B[i].shape->p_cpy, P[3]*(1.-alpha)*G1);
	B[i].shape->adm_y     (B[i].shape->p_cpy, P[4]*(1.-2.*alpha)*G1);

	B[i].shape->cpy(B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
}

//////////////////////////////////////////////////////////////////
//...realization shear derivative of classical displacements (Ux);
void CVisco2D::jump3_classic_x(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
	complex G1   = 1./get_param(NUM_SHEAR+shift), 
			alpha  = .25/(1.-get_param(NUM_SHEAR+1+shift)); m += solver.id_norm;
	B[i].shape->cpy_x		();
	B[i].shape->admittance(B[i].shape->deriv, NULL, 0., 0.);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, 0., 0.);
	B[i].shape->deriv_X	(B[i].shape->deriv, -P[4]);
	B[i].shape->deriv_Y	(B[i].shape->p_cpy, P[3]);
	B[i].shape->admittance(B[i].shape->deriv, B[i].shape->p_cpy, 1., 1.);
	B[i].shape->cpy_x		();

	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[0]*alpha*G1, 0.);
	B[i].shape->adm_x     (B[i].shape->deriv,-P[4]*(1.-2.*alpha)*G1);
	B[i].shape->adm_y     (B[i].shape->deriv, P[3]*(1.-alpha)*G1);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, -P[1]*alpha*G1, 0.);
	B[i].shape->adm_x     (B[i].shape->p_cpy, -P[3]*alpha*G1);

	B[i].shape->cpy(B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
}

//////////////////////////////////////////////////////////////////
//...realization shear derivative of classical displacements (Uy);
void CVisco2D::jump3_classic_y(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
	complex G1   = 1./get_param(NUM_SHEAR+shift), 
			alpha  = .25/(1.-get_param(NUM_SHEAR+1+shift)); m += solver.id_norm;
	B[i].shape->cpy_y		();
	B[i].shape->admittance(B[i].shape->deriv, NULL, 0., 0.);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, 0., 0.);
	B[i].shape->deriv_X	(B[i].shape->deriv, -P[4]);
	B[i].shape->deriv_Y	(B[i].shape->p_cpy, P[3]);
	B[i].shape->admittance(B[i].shape->deriv, B[i].shape->p_cpy, 1., 1.);
	B[i].shape->cpy_y		();

	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[0]*alpha*G1, 0.);
	B[i].shape->adm_y     (B[i].shape->deriv, P[4]*alpha*G1);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, -P[1]*alpha*G1, 0.);
	B[i].shape->adm_x     (B[i].shape->p_cpy,-P[4]*(1.-alpha)*G1);
	B[i].shape->adm_y     (B[i].shape->p_cpy, P[3]*(1.-2.*alpha)*G1);

	B[i].shape->cpy(B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
}

////////////////////////////////////////////////////////////////////
//...realization of surface forces for classical displacements (Px);
void CVisco2D::jump4_classic_x(double * P, int i, int m)
{
	int     shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
	complex alpha = .5/(1.-get_param(NUM_SHEAR+1+shift)); m += solver.id_norm;
	B[i].shape->cpy_x  ();
	B[i].shape->deriv_N();
	B[i].shape->cpy_x  ();
	B[i].shape->admittance(B[i].shape->deriv, NULL, -alpha, 0.);
	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[0], 0.);
	B[i].shape->adm_x     (B[i].shape->deriv, P[3]);
	B[i].shape->adm_y     (B[i].shape->deriv, P[4]*(1.-alpha));
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, P[1], 0.);
	B[i].shape->adm_y     (B[i].shape->p_cpy, P[3]*(2.*alpha-1.));
	B[i].shape->adm_x     (B[i].shape->p_cpy, P[4]*(1.-alpha));

	B[i].shape->cpy(B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
}

////////////////////////////////////////////////////////////////////
//...realization of surface forces for classical displacements (Py);
void CVisco2D::jump4_classic_y(double * P, int i, int m)
{
	int     shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
	complex alpha = .5/(1.-get_param(NUM_SHEAR+1+shift)); m += solver.id_norm;
	B[i].shape->cpy_y  ();
	B[i].shape->deriv_N();
	B[i].shape->cpy_y  ();
	B[i].shape->admittance(B[i].shape->deriv, NULL, -alpha, 0.);
	B[i].shape->cpy       (B[i].shape->deriv, B[i].shape->p_cpy);
	B[i].shape->admittance(B[i].shape->deriv, NULL, P[0], 0.);
	B[i].shape->adm_x     (B[i].shape->deriv, P[4]*(2.*alpha-1.));
	B[i].shape->adm_y     (B[i].shape->deriv, P[3]*(1.-alpha));
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, P[1], 0.);
	B[i].shape->adm_x     (B[i].shape->p_cpy, P[3]*(1.-alpha));
	B[i].shape->adm_y     (B[i].shape->p_cpy, P[4]);

	B[i].shape->cpy(B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
}

/////////////////////////////////////////////////////////////////////////
//...realization double shear derivative of classical displacements (Ux);
void CVisco2D::jump5_classic_x(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
	complex G1   = 1./get_param(NUM_SHEAR+shift), 
			alpha  = .25/(1.-get_param(NUM_SHEAR+1+shift)); m += solver.id_norm;
	B[i].shape->admittance(B[i].shape->deriv, NULL, 0., 0.);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, 0., 0.);
	B[i].shape->cpy_xx	 ();
	B[i].shape->deriv_X	 (B[i].shape->deriv, sqr(P[4]));
	B[i].shape->adm		 (B[i].shape->p_cpy, sqr(P[4]));
	B[i].shape->cpy_xx	 ();

	B[i].shape->cpy_yy	 ();
	B[i].shape->deriv_X	 (B[i].shape->deriv, sqr(P[3]));
	B[i].shape->adm		 (B[i].shape->p_cpy, sqr(P[3]));
	B[i].shape->cpy_yy	 ();

	B[i].shape->cpy_xy	 ();
	B[i].shape->deriv_X	 (B[i].shape->deriv, -P[4]*P[3]*2.);
	B[i].shape->adm		 (B[i].shape->p_cpy, -P[4]*P[3]*2.);
	B[i].shape->cpy_xy	 ();

	B[i].shape->admittance(B[i].shape->p_cpy, B[i].shape->deriv, (1.-alpha)*G1, -P[0]*alpha*G1);
	B[i].shape->adm_xx    (B[i].shape->p_cpy, -sqr(P[4])*2.*alpha*G1);
	B[i].shape->adm_xy    (B[i].shape->p_cpy,  P[4]*P[3]*2.*alpha*G1);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[1]*alpha*G1, 0.);
	B[i].shape->adm_xy    (B[i].shape->deriv, -sqr(P[3])*2.*alpha*G1);
	B[i].shape->adm_xx    (B[i].shape->deriv,  P[4]*P[3]*2.*alpha*G1);

	B[i].shape->cpy(B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
}

/////////////////////////////////////////////////////////////////////////
//...realization double shear derivative of classical displacements (Uy);
void CVisco2D::jump5_classic_y(double * P, int i, int m)
{
	int    shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT;
	complex G1   = 1./get_param(NUM_SHEAR+shift), 
			alpha  = .25/(1.-get_param(NUM_SHEAR+1+shift)); m += solver.id_norm;
	B[i].shape->admittance(B[i].shape->deriv, NULL, 0., 0.);
	B[i].shape->admittance(B[i].shape->p_cpy, NULL, 0., 0.);
	B[i].shape->cpy_xx	 ();
	B[i].shape->deriv_Y	 (B[i].shape->deriv, sqr(P[4]));
	B[i].shape->adm		 (B[i].shape->p_cpy, sqr(P[4]));
	B[i].shape->cpy_xx	 ();

	B[i].shape->cpy_yy	 ();
	B[i].shape->deriv_Y	 (B[i].shape->deriv, sqr(P[3]));
	B[i].shape->adm		 (B[i].shape->p_cpy, sqr(P[3]));
	B[i].shape->cpy_yy	 ();

	B[i].shape->cpy_xy	 ();
	B[i].shape->deriv_Y	 (B[i].shape->deriv, -P[4]*P[3]*2.);
	B[i].shape->adm		 (B[i].shape->p_cpy, -P[4]*P[3]*2.);
	B[i].shape->cpy_xy	 ();

	B[i].shape->admittance(B[i].shape->p_cpy, B[i].shape->deriv, (1.-alpha)*G1, -P[1]*alpha*G1);
	B[i].shape->adm_yy    (B[i].shape->p_cpy, -sqr(P[3])*2.*alpha*G1);
	B[i].shape->adm_xy    (B[i].shape->p_cpy,  P[4]*P[3]*2.*alpha*G1);
	B[i].shape->admittance(B[i].shape->deriv, NULL, -P[0]*alpha*G1, 0.);
	B[i].shape->adm_xy    (B[i].shape->deriv, -sqr(P[4])*2.*alpha*G1);
	B[i].shape->adm_yy    (B[i].shape->deriv,  P[4]*P[3]*2.*alpha*G1);

	B[i].shape->cpy(B[i].shape->deriv, B[i].shape->FULL(solver.hh[i][0][m], 0));
	B[i].shape->cpy(B[i].shape->p_cpy, B[i].shape->FULL(solver.hh[i][0][m], 0, 1));
}

//////////////////////////////////////////////
//...transformation of the collocation vector;
void CVisco2D::jump_make_local(int i, int m)
{
	if (! solver.dim || ! solver.hh[i][0].GetMatrix()) return;
	m  += solver.id_norm;
	for (int j = 0; j < solver.dim[i]; j++) {
		complex P[] = {solver.hh[i][0][m  ][j], 
							solver.hh[i][0][m+1][j], 0.};
		B[i].shape->norm_local_T(P);
		solver.hh[i][0][m  ][j] = P[0];
		solver.hh[i][0][m+1][j] = P[1];
	}
}

void CVisco2D::jump_make_common(int i, int m)
{
	if (! solver.dim || ! solver.hh[i][0].GetMatrix()) return;
	m  += solver.id_norm;
	for (int j = 0; j < solver.dim[i]; j++) {
		complex P[] = {solver.hh[i][0][m  ][j], 
							solver.hh[i][0][m+1][j], 0.};
		B[i].shape->norm_common_T(P);
		solver.hh[i][0][m  ][j] = P[0];
		solver.hh[i][0][m+1][j] = P[1];
	}
}

///////////////////////////////////////////////////////////////////////////
//...inclusion of the boundary condition data to the solver for all blocks;
Num_State CVisco2D::gram1(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double G1 = get_param(NUM_SHEAR), hx, hy, p3, f, P[6];
		int 	 m  = solver.id_norm;

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(FULLY_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////////////////////////
//...realization of pattern cell boundary condition for Gram matrix;
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

//////////////////////////////////////
//...jump of all displacement moments;
			B[i].shape->parametrization_grad(P, 2);

			if (p3 == MIN_HIT || p3 == NUMS_BND || p3 == MAX_HIT) {
				jump1_classic_x(P, i, 0); solver.admittance(i, 0, G1);
				jump1_classic_y(P, i, 1); solver.admittance(i, 1, G1);
			}
			else
			if (0. <= p3 && p3 <= (double)(NUMS_BND-SPECIAL_BND)) {
				jump4_classic_x(P, i, 0); 
				jump4_classic_y(P, i, 1); 
			}
			else
			if (p3 == (double)(NORMS_BND-SPECIAL_BND)) {
				jump2_classic_x(P, i, 0); solver.admittance(i, 0, G1);
				jump2_classic_y(P, i, 1); solver.admittance(i, 1, G1);
			}
			jump_make_common(i, 0);

////////////////////////////
//...composition functional;
			solver.to_equationDD(i, solver.hh[i][0][m],	   solver.hh[i][0][m], f);
			solver.to_equationDD(i, solver.hh[i][0][m+1], solver.hh[i][0][m+1], f);
		
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

/////////////////////////////////////////////////////////////////////
//...junction data of the periodic boundary condition for all blocks;
Num_State CVisco2D::gram2(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double G1 = get_param(NUM_SHEAR), AX, AY, f, P[6], TX, TY, hx, 
				 g1 = G1*.5, f1 = 1., g2 = G1*.5, g0 = -G1,  requl = get_param(NUM_GEOMT);
      int id_isolated = 0, m  = solver.id_norm, id_dir, k, j, first = 1, k0, j0;
		if (id_isolated) {
			g0 = g1 = G1;
			f1 = g2 = 0.;
		}
/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(FULLY_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

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
					case 1: TX =  AX; hx = -g1*AX; break;
					case 2: TX = -AX; hx =  g1*AX; break;
					case 3: TY =  AY; break;
					case 4: TY = -AY; break;
				}
				B[k].mp[1] -= TX;
				B[k].mp[2] -= TY; B[k].shape->set_local_P0(B[k].mp+1);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m+3; num >= m; num--) {
					 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}
				if (first && solver.mode(REGULARIZATION) && (id_dir == 1 || id_dir == 2)) {
					for (int num = m+4; num < solver.n; num--) {
						 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
						 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
					}
					first = 0; k0 = k; j0 = j;
				}

//////////////////////////////////////
//...jump of all displacement moments;
				B[i].shape->parametrization_grad(P, 2); 

				jump1_classic_x(P, i, 0); jump1_classic_y(P, i, 1); 
				jump2_classic_x(P, i, 2); jump2_classic_y(P, i, 3);
				jump_make_common(i, 0);	  jump_make_common(i, 2);

				if (! first) { //...интегрирование перемещений;
					solver.admittance (i, 4, 1., 0, f*G1); 
					solver.admittance (i, 5, 1., 1, f*G1);
				}
				else if (requl) {
					solver.admittance (i, 4, 0., 0, G1); 
					solver.admittance (i, 5, 0., 1, G1);
				}
				else {
					solver.admittance (i, 4, 0., 1, G1); 
					solver.admittance (i, 5, 0., 0, G1);
				}
				solver.admittance(i, 0, g1, 2, g2); solver.admittance(i, 2, g0, 0, f1); 
				solver.admittance(i, 1, g1, 3, g2); solver.admittance(i, 3, g0, 1, f1); 

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);

				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->parametrization_grad(P, 2);

				jump1_classic_x(P, k, 0); jump1_classic_y(P, k, 1); 
				jump2_classic_x(P, k, 2); jump2_classic_y(P, k, 3);
				jump_make_common(k, 0);	  jump_make_common(k, 2);

				if (! first) { //...интегрирование перемещений;
					solver.admittance (k, 4, 1., 0, -f*G1); 
					solver.admittance (k, 5, 1., 1, -f*G1);
				}
				else if (requl) {
					solver.admittance (k, 4, 0., 0, -G1); 
					solver.admittance (k, 5, 0., 1, -G1);
				}
				else {
					solver.admittance (k, 4, 0., 1, -G1); 
					solver.admittance (k, 5, 0., 0, -G1);
				}
				solver.admittance(k, 0, g1, 2, g2); solver.admittance(k, 2, g0, 0, f1); 
				solver.admittance(k, 1, g1, 3, g2); solver.admittance(k, 3, g0, 1, f1); 

////////////////////////////
//...composition functional;
				solver.to_transferTR(i, j, solver.hh[i][0][m],   solver.hh[k][0][m],   f);
				solver.to_transferDD(i, j, solver.hh[i][0][m],   solver.hh[k][0][m+2], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+3], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);

				if (fabs(hx) > EE) {
				  solver.to_equationHH(i, 0, solver.hh[i][0][m],   hx*f);
				  solver.to_equationHH(i, 1, solver.hh[i][0][m+1], hx*f);

				  solver.to_equationHL(k, 0, solver.hh[k][0][m+2], -hx*f);
				  solver.to_equationHL(k, 1, solver.hh[k][0][m+3], -hx*f);
				}
				if (first && solver.mode(REGUL_BOUNDARY)) {//...регуляризация матрицы через граничное условие;
					if (id_dir == 1 || id_dir == 2) {
						solver.to_transferTR(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
						solver.to_transferDD(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
						solver.to_transferTL(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
					}
					if (id_dir == 3 || id_dir == 4) {
						solver.to_transferTR(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
						solver.to_transferDD(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
						solver.to_transferTL(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
					}
				}
				B[k].mp[1] += TX;
				B[k].mp[2] += TY; B[k].shape->set_local_P0(B[k].mp+1);
			}
		}
		if (! first) {//...регуляризация матрицы по интегралу первого блока;
			solver.clean_mode(REGULARIZATION);
			solver.to_transferTR(i, j0, solver.hh[i][0][m+4], solver.hh[k0][0][m+4], requl);
			solver.to_transferDD(i, j0, solver.hh[i][0][m+4], solver.hh[k0][0][m+4], requl);
			solver.to_transferTL(i, j0, solver.hh[i][0][m+4], solver.hh[k0][0][m+4], requl);

			solver.to_transferTR(i, j0, solver.hh[i][0][m+5], solver.hh[k0][0][m+5], requl);
			solver.to_transferDD(i, j0, solver.hh[i][0][m+5], solver.hh[k0][0][m+5], requl);
			solver.to_transferTL(i, j0, solver.hh[i][0][m+5], solver.hh[k0][0][m+5], requl);
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

////////////////////////////////////////////////////////////////
//...inclusion of the joining data to the solver for all blocks;
Num_State CVisco2D::transfer1(CGrid * nd, int i, int k, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N) {
      double G1 = get_param(NUM_SHEAR), f, P[6], 
				 f1 = 1., g1 = G1*.5, g2 = G1*.5, g0 = -G1;
      int id_isolated = 0, m = solver.id_norm;
		if (id_isolated) {
			g0 = g1 = G1;
			f1 = g2 = 0.;
		}
/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(FULLY_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

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

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m; num < solver.n; num++) {
					 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}

//////////////////////////////////////
//...jump of all displacement moments;
				B[i].shape->parametrization_grad(P, 2);

				jump1_classic_x(P, i, 0); 
				jump2_classic_x(P, i, 2); 
				solver.admittance (i, 0, g1, 2, g2); solver.admittance(i, 2, g0, 0, f1); 

				jump1_classic_y(P, i, 1); 
				jump2_classic_y(P, i, 3); 
				solver.admittance (i, 1, g1, 3, g2); solver.admittance(i, 3, g0, 1, f1); 

				jump_make_common(i, 0);
				jump_make_common(i, 2);

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);

				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->parametrization_grad(P, 2);

				jump1_classic_x(P, k, 0); 
				jump2_classic_x(P, k, 2); 
				solver.admittance (k, 0, g1, 2, g2); solver.admittance(k, 2, g0, 0, f1); 

				jump1_classic_y(P, k, 1); 
				jump2_classic_y(P, k, 3); 
				solver.admittance (k, 1, g1, 3, g2); solver.admittance(k, 3, g0, 1, f1); 
				
				jump_make_common(k, 0);
				jump_make_common(k, 2);

////////////////////////////
//...composition functional;
				solver.to_transferTR(i, j, solver.hh[i][0][m],   solver.hh[k][0][m],   f);
				solver.to_transferDD(i, j, solver.hh[i][0][m],   solver.hh[k][0][m+2], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+3], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////
//...inclusion conjunction data to the solver for all blocks;
Num_State CVisco2D::transfer2(CGrid * nd, int i, int k, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N && B && B[i].link && B[i].link[0] > NUM_PHASE) {
      double G1 = get_param(NUM_SHEAR), f, P[6], 
				 f0 = -1., f1 = 1., f2 = .5, g1 = G1*.5;
      int id_isolated = 1, id_flag = 1, m = solver.id_norm;
		if (id_isolated) {
			f1 = f2 = 0.;
			g1 = G1; f0 = 1.;	
			id_flag = B[i].link[NUM_PHASE] == -1;
		}
/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(FULLY_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

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

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m; num < solver.n; num++) {
					 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}

//////////////////////////////////////
//...jump of all displacement moments;
				B[i].shape->parametrization_grad(P, 2);

				jump1_classic_x(P, i, 0); 
				jump4_classic_x(P, i, 2); 
				solver.admittance (i, 0, g1, 2, f2); solver.admittance(i, 2, f0, 0, f1); 

				jump1_classic_y(P, i, 1); 
				jump4_classic_y(P, i, 3); 
				solver.admittance (i, 1, g1, 3, f2); solver.admittance(i, 3, f0, 1, f1); 

				jump_make_common(i, 0);
				jump_make_common(i, 2);

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);

				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->parametrization_grad(P, 2);

				jump1_classic_x(P, k, 0); 
				jump4_classic_x(P, k, 2); 
				solver.admittance (k, 0, g1, 2, f2); solver.admittance(k, 2, f0, 0, f1); 

				jump1_classic_y(P, k, 1); 
				jump4_classic_y(P, k, 3); 
				solver.admittance (k, 1, g1, 3, f2); solver.admittance(k, 3, f0, 1, f1); 

				jump_make_common(k, 0);
				jump_make_common(k, 2);

///////////////////////////
//...composition functional;
				if (id_flag) {
					solver.to_transferTR(i, j, solver.hh[i][0][m],   solver.hh[k][0][m],   f);
					solver.to_transferDD(i, j, solver.hh[i][0][m],   solver.hh[k][0][m+2], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);

					solver.to_transferTR(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
					solver.to_transferDD(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+3], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);
				}
				else {
					solver.to_transferTR(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m+2], f);
					solver.to_transferDD(i, j, solver.hh[i][0][m+2], solver.hh[k][0][m],   f);
					solver.to_transferTL(i, j, solver.hh[i][0][m],	 solver.hh[k][0][m],	  f);

					solver.to_transferTR(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+3], f);
					solver.to_transferDD(i, j, solver.hh[i][0][m+3], solver.hh[k][0][m+1], f);
					solver.to_transferTL(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				}
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////
//...формирование матрицы Грама с учетом функционала энергии;
Num_State CVisco2D::gram3(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double G1 = get_param(NUM_SHEAR), AX, AY, f, P[6], TX, TY, hx,  requl = get_param(NUM_GEOMT);
      int	 m  = solver.id_norm, id_dir, k, j, first = 1, k0, j0;

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(FULLY_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

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
				for (int num = m+3; num >= m; num--) {
					 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}
				if (first && solver.mode(REGULARIZATION) && (id_dir == 1 || id_dir == 2)) {
					for (int num = m+4; num < solver.n; num--) {
						 memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
						 memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
					}
					first = 0; k0 = k; j0 = j;
				}

//////////////////////////////////////
//...jump of all displacement moments;
				B[i].shape->parametrization_grad(P, 2);

				jump1_classic_x(P, i, 0); jump1_classic_y(P, i, 1);
				jump4_classic_x(P, i, 2); jump4_classic_y(P, i, 3);
				jump_make_common(i, 0);	  jump_make_common(i, 2);

				if (! first) { //...интегрирование перемещений;
					solver.admittance (i, 4, 1., 0, f*G1); 
					solver.admittance (i, 5, 1., 1, f*G1);
				}
				else if (requl) {
					solver.admittance (i, 4, 0., 0, G1); 
					solver.admittance (i, 5, 0., 1, G1);
				}
				else {
					solver.admittance (i, 4, 0., 1, G1); 
					solver.admittance (i, 5, 0., 0, G1);
				} 
				solver.admittance(i, 0, G1); solver.admittance(i, 1, G1);

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);

				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->parametrization_grad(P, 2);

				jump1_classic_x(P, k, 0); jump1_classic_y(P, k, 1); 
				jump4_classic_x(P, k, 2); jump4_classic_y(P, k, 3);
				jump_make_common(k, 0);	  jump_make_common(k, 2);

				if (! first) { //...интегрирование перемещений;
					solver.admittance (k, 4, 1., 0, -f*G1); 
					solver.admittance (k, 5, 1., 1, -f*G1);
				}
				else if (requl) {
					solver.admittance (k, 4, 0., 0, -G1); 
					solver.admittance (k, 5, 0., 1, -G1);
				}
				else {
					solver.admittance (k, 4, 0., 1, -G1); 
					solver.admittance (k, 5, 0., 0, -G1);
				}
				solver.admittance(k, 0, G1); solver.admittance(k, 1, G1);

/////////////////////////////////////////////////
//...сшивка функций методом наименьших квадратов;
				solver.to_transferTR(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
			
				if (fabs(hx) > EE) {
				  solver.to_equationHH(i, 0, solver.hh[i][0][m],   hx*f);
				  solver.to_equationHH(i, 1, solver.hh[i][0][m+1], hx*f);

				  solver.to_equationHL(k, 0, solver.hh[k][0][m],   -hx*f);
				  solver.to_equationHL(k, 1, solver.hh[k][0][m+1], -hx*f);
				}
				if (first && solver.mode(REGUL_BOUNDARY)) {//...регуляризация матрицы через граничное условие;
					if (id_dir == 1 || id_dir == 2) {
						solver.to_transferTR(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
						solver.to_transferDD(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
						solver.to_transferTL(i, j, solver.hh[i][0][m+4], solver.hh[k][0][m+4], f);
					}
					if (id_dir == 3 || id_dir == 4) {
						solver.to_transferTR(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
						solver.to_transferDD(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
						solver.to_transferTL(i, j, solver.hh[i][0][m+5], solver.hh[k][0][m+5], f);
					}
				}
				
/////////////////////////////////////////////////
//...энергетические слагаемые (нормированные G1);
				solver.to_equationER(i, solver.hh[i][0][m], solver.hh[i][0][m+2], -f);
				solver.to_equationEL(k, solver.hh[k][0][m], solver.hh[k][0][m+2],  f);

				solver.to_equationER(i, solver.hh[i][0][m+1], solver.hh[i][0][m+3], -f);
				solver.to_equationEL(k, solver.hh[k][0][m+1], solver.hh[k][0][m+3],  f);

				B[k].mp[1] += TX;
				B[k].mp[2] += TY; B[k].shape->set_local_P0(B[k].mp+1);
			}
		}
		if (! first) {//...регуляризация матрицы по интегралу первого блока;
			solver.clean_mode(REGULARIZATION);
			solver.to_transferTR(i, j0, solver.hh[i][0][m+4], solver.hh[k0][0][m+4], requl);
			solver.to_transferDD(i, j0, solver.hh[i][0][m+4], solver.hh[k0][0][m+4], requl);
			solver.to_transferTL(i, j0, solver.hh[i][0][m+4], solver.hh[k0][0][m+4], requl);

			solver.to_transferTR(i, j0, solver.hh[i][0][m+5], solver.hh[k0][0][m+5], requl);
			solver.to_transferDD(i, j0, solver.hh[i][0][m+5], solver.hh[k0][0][m+5], requl);
			solver.to_transferTL(i, j0, solver.hh[i][0][m+5], solver.hh[k0][0][m+5], requl);
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

///////////////////////////////////////////////////////////////
//...формирование матриц перехода с учетом функционала энергии;
Num_State CVisco2D::transfer3(CGrid * nd, int i, int k, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N) {
      double G1 = get_param(NUM_SHEAR), f, P[6];
      int m = solver.id_norm, shift;

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(FULLY_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

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

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m; num < solver.n; num++) {
					memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));
					memset(solver.hh[k][0][num], 0, solver.dim[k]*sizeof(double));
				}

//////////////////////////////////////
//...jump of all displacement moments;
				B[i].shape->parametrization_grad(P, 2);

				jump1_classic_x(P, i, 0); jump1_classic_y(P, i, 1); 
				jump4_classic_x(P, i, 2); jump4_classic_y(P, i, 3); 
				jump_make_common(i, 0);	  jump_make_common(i, 2);
				solver.admittance(i, 0, G1); 
				solver.admittance(i, 1, G1); 

				B[i].shape->make_common(P);
				B[i].shape->norm_common(P+3);

				B[k].shape->make_local(P);
				B[k].shape->norm_local(P+3);

				B[k].shape->parametrization_grad(P, 2);

				jump1_classic_x(P, k, 0); jump1_classic_y(P, k, 1); 
				jump4_classic_x(P, k, 2); jump4_classic_y(P, k, 3); 
				jump_make_common(k, 0);   jump_make_common(k, 2);
				solver.admittance(k, 0, G1); 
				solver.admittance(k, 1, G1); 

/////////////////////////////////////////////////
//...сшивка функций методом наименьших квадратов;
				solver.to_transferTR(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m], solver.hh[k][0][m], f);

				solver.to_transferTR(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferDD(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);
				solver.to_transferTL(i, j, solver.hh[i][0][m+1], solver.hh[k][0][m+1], f);

/////////////////////////////////////////////////
//...энергетические слагаемые (нормированные G1);
				if (B[i].link[NUM_PHASE] != B[k].link[NUM_PHASE] && B[i].link[NUM_PHASE] < -1)  f = -f;
				solver.to_equationER(i, solver.hh[i][0][m], solver.hh[i][0][m+2], -f);
				solver.to_equationEL(k, solver.hh[k][0][m], solver.hh[k][0][m+2],  f);

				solver.to_equationER(i, solver.hh[i][0][m+1], solver.hh[i][0][m+3], -f);
				solver.to_equationEL(k, solver.hh[k][0][m+1], solver.hh[k][0][m+3],  f);

////////////////////////////////////////////////////////////////
//...добавляем поверхностную энергию адгезии (нормировка ?????);
				if (B[i].link[NUM_PHASE] != B[k].link[NUM_PHASE] && 0) {
					double K1, K2, d1, d2;
					K1 = get_param(NUM_SHEAR-2+(shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT));
					d1 = get_param(NUM_SHEAR-1+ shift); 
					K2 = get_param(NUM_SHEAR-2+(shift = (-B[k].link[NUM_PHASE]-1)*NUM_SHIFT));
					d2 = get_param(NUM_SHEAR-1+ shift);

					jump3_classic_x(P, k, 4); 
					jump3_classic_y(P, k, 5); 
					solver.admittance(k, 6,  0.,   4, P[3]);
					solver.admittance(k, 4, -P[4], 5, P[3]);
					solver.admittance(k, 5,  P[4], 6, 1.);

					B[k].shape->make_local(P);
					B[k].shape->norm_local(P+3);

					B[i].shape->make_common(P);
					B[i].shape->norm_common(P+3);

					jump3_classic_x(P, i, 4); 
					jump3_classic_y(P, i, 5); 
					solver.admittance (i, 6,  0.,   4, P[3]);
					solver.admittance (i, 4, -P[4], 5, P[3]);
					solver.admittance (i, 5,  P[4], 6, 1.);
					
					solver.to_equationER(i, solver.hh[i][0][m+4], solver.hh[i][0][m+4],  f*K1);
					solver.to_equationEL(k, solver.hh[k][0][m+4], solver.hh[k][0][m+4], -f*K2);

					solver.to_equationER(i, solver.hh[i][0][m+5], solver.hh[i][0][m+5],  f*d1);
					solver.to_equationEL(k, solver.hh[k][0][m+5], solver.hh[k][0][m+5], -f*d2);
				}
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

///////////////////////////////////////////////////////////////////////////
//...inclusion of the boundary condition data to the solver for all blocks;
Num_State CVisco2D::gram4(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < solver.N && B[i].shape && B[i].mp) {
		double G1 = get_param(NUM_SHEAR), hx, hy, p3, f, P[6];
		int 	 m  = solver.id_norm;

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(FULLY_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////////////////////////
//...realization of pattern cell boundary condition for Gram matrix;
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

//////////////////////////////////////
//...jump of all displacement moments;
			B[i].shape->parametrization_grad(P, 2);

			jump1_classic_x(P, i, 0); jump1_classic_y(P, i, 1); 
			jump4_classic_x(P, i, 2); jump4_classic_y(P, i, 3); 	
			jump_make_common(i, 0);	  jump_make_common(i, 2);
			solver.admittance(i, 0, G1); 
			solver.admittance(i, 1, G1); 

////////////////////////////////////
//...composition collocation vector;
			if (p3 == MIN_HIT || p3 == NUMS_BND || p3 == MAX_HIT) {
				solver.admittance(i, 4, 0., 0, 1.);
				solver.admittance(i, 5, 0., 1, 1.);
			}
			else
			if (p3 == (double)(NORMS_BND-SPECIAL_BND)) {
				jump2_classic_x(P, i, 4); solver.admittance(i, 4, G1);
				jump2_classic_y(P, i, 5); solver.admittance(i, 5, G1);
				jump_make_common(i, 4);
			}

////////////////////////////////////////////////////
//...граничные условия методом наименьших квадратов;
			if (p3 == MIN_HIT || p3 == (double)(NORMS_BND-SPECIAL_BND) || p3 == NUMS_BND || p3 == MAX_HIT) {
				if (fabs(hx) > EE)
				solver.to_equationHH(i, 0, solver.hh[i][0][m+4], hx*G1*f);
				solver.to_equationDD(i, solver.hh[i][0][m+4], solver.hh[i][0][m+4], f);
				
				if (fabs(hy) > EE)
				solver.to_equationHH(i, 0, solver.hh[i][0][m+5], hy*G1*f);
				solver.to_equationDD(i, solver.hh[i][0][m+5], solver.hh[i][0][m+5], f);
			}

/////////////////////////////////////////////////
//...энергетические слагаемые (нормированные G1);
			solver.to_equationEE(i, solver.hh[i][0][m],	 solver.hh[i][0][m+2], -f);
			solver.to_equationEE(i, solver.hh[i][0][m+1], solver.hh[i][0][m+3], -f);

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

//////////////////////////////////////////////////////////////////////////
//...интегрирование НДС на заданном наборе узлов для периодической задачи;
Num_State CVisco2D::rigidy1(CGrid * nd, int i, complex * K)
{
	if (nd) {
      int l, shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, m = solver.id_norm;
      double alpha = get_param(NUM_SHEAR+shift+1)/(.5-get_param(NUM_SHEAR+shift+1)), 
					 G0 = get_param(NUM_SHEAR+shift), f, P[6];

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (solver.mode(FULLY_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
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

/////////////////////////////////////////////////////////////////////////////
//...вычисляем функции для формирования интеграла от деформаций и напряжений;
			B[i].shape->parametrization_grad(P, 2);
			jump1_classic_x(P, i, 0);
			jump1_classic_y(P, i, 1);

			//UX = B[i].shape->potential(solver.hh[i][0][m],   0);
			//UY = B[i].shape->potential(solver.hh[i][0][m+1], 0);
			//
			//K[0] -= G0*(UX*P[3]*2.+alpha*(UX*P[3]+UY*P[4]))*f;
			//K[1] -= G0*(UX*P[4]+UY*P[3])*f;
			//K[2] -= G0*(UY*P[4]*2.+alpha*(UX*P[3]+UY*P[4]))*f;
			//K[3] -= UX*P[3]*f;
			//K[4] -=(UX*P[4]+UY*P[3])*f*.5;
			//K[5] -= UY*P[4]*f;

			//UX = B[i].shape->potential(solver.hh[i][0][m],   1);
			//UY = B[i].shape->potential(solver.hh[i][0][m+1], 1);

			//K[6]  -= G0*(UX*P[3]*2.+alpha*(UX*P[3]+UY*P[4]))*f;
			//K[7]  -= G0*(UX*P[4]+UY*P[3])*f;
			//K[8]  -= G0*(UY*P[4]*2.+alpha*(UX*P[3]+UY*P[4]))*f;
			//K[9]  -= UX*P[3]*f;
			//K[10] -=(UX*P[4]+UY*P[3])*f*.5;
			//K[11] -= UY*P[4]*f;
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////
//...auxiliary function for parameters of doubly connected media;
void CVisco2D::set_fasa_hmg(double GM_re, double GM_im, double KM_re, double KM_im, double GI_re, double GI_im, double KI_re, double KI_im, double GL_re, double GL_im, double KL_re, double KL_im)
{
	if (size_of_param() > NUM_SHEAR+3+NUM_SHIFT*2) {
		param[NUM_SHEAR] = GM_re;
		param[NUM_SHEAR+NUM_SHIFT] = GI_re;
		param[NUM_SHEAR+NUM_SHIFT*2] = GL_re;
		param[NUM_SHEAR+1] = GM_im;
		param[NUM_SHEAR+1+NUM_SHIFT] = GI_im;
		param[NUM_SHEAR+1+NUM_SHIFT*2] = GL_im;
		param[NUM_SHEAR+2] = KM_re;
		param[NUM_SHEAR+2+NUM_SHIFT] = KI_re;
		param[NUM_SHEAR+2+NUM_SHIFT*2] = KL_re;
		param[NUM_SHEAR+3] = KM_im;
		param[NUM_SHEAR+3+NUM_SHIFT] = KI_im;
		param[NUM_SHEAR+3+NUM_SHIFT*2] = KL_im;
	}
}

//////////////////////////////////////////////////////////
//...counting header for solving plane elasticity problem;
Num_State CVisco2D::counting_header(int Num_comput)
{
	int i, j, k, elem, id_dir, n_rhs = 2;
	char msg[201];

//	solver.set_mode(REDUCED_MESSAGE); //...отключаем подробную печать;
	if (! solver.mode(NO_MESSAGE)) {
		Message(" ");
		sprintf(msg, "CVisco2D_draft: N_sm = %d, N_mpl = %d, G1 = %g, G2 = %g", N,
				  UnPackInts(get_param(0)), get_param(NUM_SHEAR), get_param(NUM_SHEAR+NUM_SHIFT));
		Message(msg);

		Message(" ");
		Message("Junction counting...");

		switch (Num_comput){
			case   BASIC_COMPUT: Message("Analytical Blocks..."); break;
			case MAPPING_COMPUT: Message("FEM Blocks...");			break;
		}
		Message(" ");
	}

//////////////////////////////////
//...определяем блочную структуру;
	solver.set_blocks(N, n_rhs); //<==== number of saved potentials !!!
	solver.n += 9;//<==== number of additional auxilliary arrays!!!
	for (k = 0; k < solver.N;  k++)
		  solver.set_links(k, B[k].link);

	shapes_init(INITIAL_STATE);
	shapes_init(ZERO_STATE);

////////////////////////////////////////////////////////////////////
//...добавляем блоки на периодических связях и на границе включений;
	if (solv%ENERGY_SOLVING == PERIODIC_SOLVING) { 
		double par[6]; SetGeomBounding(par);
		for (k = 0; k < N; k++) SkeletonBounding(B[k], par);
		for (k = 0; k < N; k++) if (B[k].link) {
			j = 0; while ((elem = geom_plink_2D(B[k], j, id_dir, par)) >= 0)			  
			solver.add_link(k, elem);
			i = j = 0; while ((elem = block_plink_2D(B[k], i, j, id_dir, par)) >= 0)			  
			solver.add_link(k, elem);
		}
	}
	LinkPhase2D(MAX_PHASE);

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
void CVisco2D::GetFuncAllValues(double X, double Y, double Z, complex * F, int i, Num_Value id_F, int id_variant, int iparam)
{
	if (! F) return;
	double P[6]  = { X, Y, Z, 1., 0., 0.};

/////////////////////////////////////
//...operation with all input points;
	if ( 0 <= i && i < N && B[i].shape && B[i].mp) {
		int shift = (-B[i].link[NUM_PHASE]-1)*NUM_SHIFT, m = solver.id_norm;

//////////////////////////////////////////////////////
//...reset auxilliary arrays and calculation function;
		for (int num = m; num < solver.n; num++)
			memset(solver.hh[i][0][num], 0, solver.dim[i]*sizeof(double));

		B[i].shape->make_local(P);
		B[i].shape->norm_local(P+3);
		switch (id_F) {
			case DISPL_VALUE: {
///////////////////////////////
//...calculation displacements;
				B[i].shape->parametrization_grad(P, 2);

				jump1_classic_x(P, i, 0); 
				jump1_classic_y(P, i, 1); 
				jump_make_common(i, 0);

//				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
//				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);

				if ((solv%ENERGY_SOLVING) && 0) F[solv%ENERGY_SOLVING-1] -= X;
			}  break;
			case STRESS_X_VALUE: {
/////////////////////////////////////////////
//...calculation stress tensor (txx and txy);
				P[3] = 1.; P[4] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
				B[i].shape->parametrization_grad(P, 2);

				jump4_classic_x(P, i, 0); 
				jump4_classic_y(P, i, 1); 
				jump_make_common(i, 0);

//				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
//				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
        }   break;
        case STRESS_Y_VALUE: {
/////////////////////////////////////////////
//...calculation stress tensor (txy and tyy);
				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
				B[i].shape->parametrization_grad(P, 2);

				jump4_classic_x(P, i, 0); 
				jump4_classic_y(P, i, 1); 
				jump_make_common(i, 0);

//				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
//				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
        }   break;
        case DILAT_VALUE: {
/////////////////////////////////
//...dilatation of displacements;
				double alpha = .25/(1.-get_param(NUM_SHEAR+1+shift)), G = 1./get_param(NUM_SHEAR+shift);
				B[i].shape->parametrization_hess(P, 1);

//				B[i].shape->cpy_x(B[i].shape->FULL(solver.hh[i][0][m+1], 0));
//				B[i].shape->cpy_y(B[i].shape->FULL(solver.hh[i][0][m+1], 0, 1));

//				F[0] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant)*(1.-2*alpha)*G;
        }   break;
        case ROTATION_VALUE: {
///////////////////////////////
//...rotation of displacements;
				double G = .5/get_param(NUM_SHEAR+shift);
				B[i].shape->parametrization_hess(P, 1);

//				B[i].shape->adm_y(B[i].shape->FULL(solver.hh[i][0][m+1], 0), -1.);
//				B[i].shape->adm_x(B[i].shape->FULL(solver.hh[i][0][m+1], 0, 1), 1.);

//				F[0] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant)*G;
        }   break;
        case NORMAL_Y_VALUE: {
////////////////////////
//...normal derivatives;
				P[4] = 1.; P[3] = P[5] = 0.;
				B[i].shape->norm_local(P+3);
				B[i].shape->parametrization_grad(P, 2);

				jump2_classic_x(P, i, 0); 
				jump2_classic_y(P, i, 1); 
				jump_make_common(i, 0);

//				F[0] = B[i].shape->potential(solver.hh[i][0][m],   id_variant);
//				F[1] = B[i].shape->potential(solver.hh[i][0][m+1], id_variant);
        }   break;
        case POTENTIAL_VALUE: {
/////////////////////////
//...potential fX and fY;
              B[i].shape->parametrization(P, 2);//...calculation with grad;
 //             F[0]  = B[i].shape->potential(0, B[i].shape->p_cpy, id_variant);
 //             F[1]  = B[i].shape->potential(1, B[i].shape->p_cpy-B[i].shape->get_NN(), id_variant);
        }     break;
        default : F[0] = i; F[1] = 0.;
     }
  }
}

///////////////////////////////////////////////////
//...calculaion rigidity matrix for homogenization;
void CVisco2D::GetRigidy(complex * K, Num_Comput Num_comput)
{
   int N_elem = UnPackInts(get_param(3)), N_max = UnPackInts(get_param(3), 1), j, k, l, m, mm, id_dir;

////////////////////////////////////////
//...auxilliary grids for discrete norm;
	CGrid * bnd = CreateNodes(GRID_EL_NODES);

	CGrid * block_bnd = CreateNodes();
			  block_bnd->add_params(1);

	CGrid * gauss_bnd = CreateNodes(GRID_EL_NODES);
			  gauss_bnd->add_params(1);

////////////////////////////////////////////////////////////
//...discrete norm on inclusion and jump boundary condition;
	double pp[1], Po[6], par[6] = {MAX_HIT, MIN_HIT, MAX_HIT, MIN_HIT, MAX_HIT, MIN_HIT};
	for (k = 0; k < N; k++) SkeletonBounding(B[k], par);

	for (k = 0; k < N; k++) if (B[k].bar && B[k].link) {
		solver.struct_init(k, NULL_STATE);

		for (int i = 0; i < B[k].bar->graph[0]; i++) 
		if (B[k].bar->ce[i]->cells_dim() == 2 && 
			 B[k].bar->ce[i]->mp[size_of_map (B[k].bar->ce[i]->mp)] == FACET_CELL) {

/////////////////////////////////
//...накапливаем граничные точки;
			for (j = 0; j < B[k].bar->ce[i]->graph[1]; j++) {
					  
				bnd->release();

				B[k].bar->ce[B[k].bar->ce[i]->graph[j+2]]->line_correct();
				B[k].bar->ce[B[k].bar->ce[i]->graph[j+2]]->grid_cells1(bnd, 0., N_max);
			
				for (mm = 1, m = NUM_PHASE; mm && m < B[k].link[0]; m++) 
				if  (B[k].link[j+1] == -B[k].link[m+1]+SRF_STATE) mm = 0;
				if (mm && B[k].link[j+1] < 0)	block_plink_2D(B[k], l = i, m = j, id_dir, par);
	
				if  (bnd->geom) 
				for (l = 0; l < bnd->geom[0]; l++) {
					int num = bnd->geom_element(l), num_n = bnd->geom[num+1],
					  num_f = num_n+num;

					if (bnd->geom[num] == GL_LINE_STRIP) {
						for (; num < num_f-1; num++) {
							Po[0] = bnd->X[bnd->geom[num+2]];
							Po[1] = bnd->Y[bnd->geom[num+2]];
							Po[2] = 0.;
							Po[3] = bnd->X[bnd->geom[num+3]];
							Po[4] = bnd->Y[bnd->geom[num+3]];
							Po[5] = 0.;

							gauss_bnd->facet_QG(Po, N_elem, SPECIAL_STATE);
							if ((! mm || mm && B[k].link[j+1] < 0 && ! id_dir) && bar && bar->graph && bar->graph[0]) //...коррекция квадратур;
							gauss_bnd->QG_curve(bar->ce[0]->mp);

							for (int lp = 0; lp < gauss_bnd->N; lp++) {
								if (! mm && B[k].link[NUM_PHASE] < -1) { //... переворачиваем нормаль на включении;
									gauss_bnd->nX[lp] = -gauss_bnd->nX[lp];
									gauss_bnd->nY[lp] = -gauss_bnd->nY[lp];
									gauss_bnd->nZ[lp] = -gauss_bnd->nZ[lp];
								}
								pp[0] = gauss_bnd->get_param(0, lp);
								block_bnd->add_new_point(gauss_bnd->X[lp],  gauss_bnd->Y[lp],  gauss_bnd->Z[lp], 
																gauss_bnd->nX[lp], gauss_bnd->nY[lp], gauss_bnd->nZ[lp], pp); 
							}
							gauss_bnd->add_buffer(gauss_bnd->N);
						}
					}
				}
			}
		}
		RigidyAll(block_bnd, k, K, Num_comput);
		block_bnd->add_buffer(block_bnd->N);
	}
	delete block_bnd;
	delete gauss_bnd;
}

//////////////////////////////////////////
//...calculaion energy by domain integral;
void CVisco2D::GetEnergy(complex * energy, Num_Value _FMF)
{
   if (! solver.mode(NO_MESSAGE)) Message("Energy calculation...");
	int N_elem = UnPackInts(get_param(3));

	CGrid * gauss_bnd = CreateNodes(GRID_EL_NODES);
	gauss_bnd->add_params(1);

	int i, j, k, m, l, N_ini = 4, id_first;
	double pm[3*4], f;
	complex F[3] = {0., 0., 0.};

   for (k = 0; k < N; k++) if (B[k].bar && B[k].link) {
       gauss_bnd->release();

       for (i = 0; i < B[k].bar->graph[0]; i++) 
       if  (B[k].bar->ce[i]->cells_dim() == 2 && 
            B[k].bar->ce[i]->mp[size_of_map (B[k].bar->ce[i]->mp)] == FACET_CELL) {

//////////////////////////////////////////////////////////////
//...снимаем данные о контуpе в массив, представляющий фасету;
				int N_arc, arc, prev = B[k].bar->ce[i]->graph[(N_arc = B[k].bar->ce[i]->graph[1])+1], 
					 N_min = min(N_arc, N_ini), m1, m2;
				for (j = 1; j <= N_min; j++, prev = arc) {
					  arc = B[k].bar->ce[i]->graph[j+1];
					  m1  = get_num(B[k].bar->ce[i]->ce[arc]->graph, 0),
					  m2  = get_num(B[k].bar->ce[i]->ce[arc]->graph, 1);
					  if (! B[k].bar->ce[i]->ce[prev]->topo_id(m1)) swap(m1, m2);

					  pm[3*(j-1)]   = B[k].bar->ce[i]->ce[m1]->mp[1];
					  pm[3*(j-1)+1] = B[k].bar->ce[i]->ce[m1]->mp[2];
					  pm[3*(j-1)+2] = B[k].bar->ce[i]->ce[m1]->mp[3];
				}

///////////////////////////////////////////////////////////////////////////
//...определяем положение криволинейной границы и строим квадратуру Гаусса;
            if (B[k].link[0] > NUM_PHASE) { //...криволинейный блок;
					for (j = id_first = 0; ! id_first && j < B[k].bar->ce[i]->graph[1]; j++) 
					for (m = NUM_PHASE;    ! id_first && m < B[k].link[0]; m++) 
					if  (B[k].link[j+1] == -B[k].link[m+1]+SRF_STATE) id_first = 1;
					if (id_first && N_arc == 3) {
						if (j == 1) { //...переворачиваем описание контура;
							swap(pm[0], pm[3]); swap(pm[1], pm[4]); swap(pm[2], pm[5]);
							swap(pm[0], pm[6]); swap(pm[1], pm[7]); swap(pm[2], pm[8]);
						}
						if (j == 3) {
							swap(pm[0], pm[6]); swap(pm[1], pm[7]); swap(pm[2], pm[8]);
							swap(pm[0], pm[3]); swap(pm[1], pm[4]); swap(pm[2], pm[5]);
						}
						gauss_bnd->facet_QG(pm, N_elem, NULL_STATE, NULL_STATE);
						if (bar && bar->graph && bar->graph[0]) //...коррекция квадратур;
						gauss_bnd->QG_tria_curve(bar->ce[0]->mp, pm);
					}
					else
					if (id_first && N_arc == 4) {
						if (j == 1) { //...переворачиваем описание контура;
							swap(pm[0], pm[3]); swap(pm[1], pm[4]);  swap(pm[2], pm[5]);
							swap(pm[0], pm[6]); swap(pm[1], pm[7]);  swap(pm[2], pm[8]);
							swap(pm[3], pm[9]); swap(pm[4], pm[10]); swap(pm[5], pm[11]);
						}
						if (j == 2) {
							swap(pm[0], pm[9]); swap(pm[1], pm[10]); swap(pm[2], pm[11]);
							swap(pm[6], pm[9]); swap(pm[7], pm[10]); swap(pm[8], pm[11]);
							swap(pm[6], pm[3]); swap(pm[7], pm[4]);  swap(pm[8], pm[5]);
							swap(pm[0], pm[3]); swap(pm[1], pm[4]);  swap(pm[2], pm[5]);
						}
						if (j == 3) {
							swap(pm[0], pm[3]); swap(pm[1], pm[4]);  swap(pm[2], pm[5]);
						}
						if (j == 4) {
							swap(pm[6], pm[3]); swap(pm[7], pm[4]);  swap(pm[8], pm[5]);
							swap(pm[6], pm[9]); swap(pm[7], pm[10]); swap(pm[8], pm[11]);
							swap(pm[0], pm[9]); swap(pm[1], pm[10]); swap(pm[2], pm[11]);
							swap(pm[0], pm[3]); swap(pm[1], pm[4]);  swap(pm[2], pm[5]);
						}
						gauss_bnd->facet_QG(pm, N_elem, OK_STATE, NULL_STATE);
						if (bar && bar->graph && bar->graph[0]) //...коррекция квадратур;
						gauss_bnd->QG_quad_curve(NULL, bar->ce[0]->mp, pm);
					}
				}
            if (B[k].link[0] == NUM_PHASE ) { //...прямолинейный блок;
					if (N_arc == 3)
						gauss_bnd->facet_QG(pm, N_elem, NULL_STATE, NULL_STATE);
					else
					if (N_arc == 4) {
						swap(pm[0], pm[3]); swap(pm[1], pm[4]);  swap(pm[2], pm[5]);
						gauss_bnd->facet_QG(pm, N_elem, OK_STATE, NULL_STATE);
					}
				}

/////////////////////////
//...интегрируем энергию;
				for (l = 0; l < gauss_bnd->N; l++) {
					f = gauss_bnd->get_param(0, l);
					GetFuncAllValues(gauss_bnd->X[l], gauss_bnd->Y[l], gauss_bnd->Z[l], F, k, _FMF);

					energy[-B[k].link[NUM_PHASE]*3-3] += F[0]*f;
					energy[-B[k].link[NUM_PHASE]*3-2] += F[1]*f;
					energy[-B[k].link[NUM_PHASE]*3-1] += F[2]*f;
				}
		 }
	}
	delete gauss_bnd;
   if (! solver.mode(NO_MESSAGE)) Message("O'K");
}

///////////////////////////////////////////
//...calculaion energy by block functional;
void CVisco2D::GetEnergyValue(int k, double * energy)
{
   if (solver.dim && solver.dim[k] && solver.JR && solver.JR[k] && solver.hh && solver.hh[k][0].GetMatrix() && 
							solver.TL && solver.TL[k]) {
		complex ** TL = solver.TL[k][solver.JR[k][solver.JR_DIAG]-solver.JR_SHIFT].GetMatrix(), * h = solver.hh[k][0][solver.id_norm];
		int  NN = solver.dim[k];

///////////////////////
//...вычисляем энергию;
		block_shape_init(B[k], NO_STATE);
		if (h && TL) { 
			double measure = 1./get_param(NUM_SHEAR);
//			for (i = 0; i < NN; i++)
//			for (j = 0; j < NN; j++)
//				energy[-B[k].link[NUM_PHASE]*3-3] += TL[i][j]*h[i]*h[j]*measure;
		}
	}
}

/////////////////////////////////////////////////////////
//...вязкоупругая модель многих сферических тел;
complex CVisco2D::TakeSphereEshelby(int nn, double * RR, complex * KK, complex * MU)
{
	complex KH;
	if (! nn) { //...классическая трехфазная модель (комплексная);
		double c0 = RR[0]/RR[1]; c0 *= sqr(c0);
		KH = KK[1]*(1.+c0/((1.-c0)/(1.+(4./3.)*MU[1]/KK[1])+KK[1]/(KK[0]-KK[1])));
		return(KH);
	}
	int k, dim_N = 4+2*nn;

	complex ** matrix = NULL, ** hh = NULL;
	set_matrix(matrix, dim_N, dim_N);
	set_matrix(hh,  2, dim_N);

//////////////////////////////////////////////////////////////////////////////////
//...заполняем и решаем систему линейных уравнений A0...Ak, Bk...Ann+1, Bnn+1, KH;
	matrix[0][0] =  RR[0];
	matrix[0][1] = -RR[0];
	matrix[0][2] = -1./sqr(RR[0]);
//
	matrix[1][0] =  3.*KK[0]*RR[0];
	matrix[1][1] = -3.*KK[1]*RR[0];
	matrix[1][2] =  4.*MU[1]/sqr(RR[0]);
//
	for (k = 1; k <= nn; k++) { 
		matrix[2*k][2*k-1] =  RR[k];
		matrix[2*k][2*k]	 =  1./sqr(RR[k]);
		matrix[2*k][2*k+1] = -RR[k];
		matrix[2*k][2*k+2] = -1./sqr(RR[k]);
//
		matrix[2*k+1][2*k-1] =  3.*KK[k]*RR[k];
		matrix[2*k+1][2*k]	= -4.*MU[k]/sqr(RR[k]);
		matrix[2*k+1][2*k+1] = -3.*KK[k+1]*RR[k];
		matrix[2*k+1][2*k+2] =  4.*MU[k+1]/sqr(RR[k]);
	}
//
	matrix[2*k][2*k-1] = RR[k];
	matrix[2*k][2*k]	 = 1./sqr(RR[k]);
	hh [1][2*k] = RR[k];
//
	matrix[2*k+1][2*k-1] =  3.*KK[k]*RR[k];
	matrix[2*k+1][2*k]	= -4.*MU[k]/sqr(RR[k]);
	matrix[2*k+1][2*k+1] = -3.*RR[k];
//
	solver.pivot_init    (dim_N);
	solver.GaussJ(matrix, dim_N);
	int i, l;
	for (               i = 0; i < dim_N; i++)
	for (hh[0][i] = 0., l = 0; l < dim_N; l++) hh[0][i] += matrix[i][l]*hh[1][l];

	solver.pivot_init(dim_N);
	KH = hh[0][2*k+1];

	delete_struct(matrix);
	delete_struct(hh);

	return(KH);
}

/////////////////////////////////////////////////////////////
//...комплексная трехфазна модель для сферического включения;
complex CVisco2D::TakeSphereEshelby_grad (double R, double c0, double * par)
{
	complex  k1 = comp(par[0]+4./3.*par[2], par[1]+4./3.*par[3]), K1 = comp(par[0], par[1]),
				k2 = comp(par[4]+4./3.*par[6], par[5]+4./3.*par[7]), K2 = comp(par[4], par[5]),
				C2 = comp(par[8]),
				C1 = comp(par[9]), KH = comp(0.);
	double	R1 = R, R2 = R/pow(c0, 1./3.), r1 = abs(C1/k1), fi1 = arg2(C1/k1), r2 = abs(C2/k2), fi2 = arg2(C2/k2);
	complex	kk1 = polar(sqrt(r1), fi1*.5), t1 = comp(1.), tp1 = polar(exp(-2.*real(kk1)*R1), -2.*imag(kk1)*R1), sh1 = (t1-tp1)*.5, ch1 = (t1+tp1)*.5,
				kk2 = polar(sqrt(r1), fi1*.5), t2 = comp(1.), tp2 = polar(exp(-2.*real(kk2)*R1), -2.*imag(kk2)*R1), sh2 = (t2-tp2)*.5, ch2 = (t2+tp2)*.5,
				t3  = polar(exp( real(kk2)*(R2-R1)),  imag(kk2)*(R2-R1)), 
				tp3 = polar(exp(-real(kk2)*(R2+R1)), -imag(kk2)*(R2-R1)), sh3 = (t3-tp3)*.5, ch3 = (t3+tp3)*.5;
	int dim_N = 7;

	if (C1 == 0. || C2 == 0.) {
		KH = K2+c0*(K1-K2)/(1.+(1.-c0)*(K1-K2)/k2);
		return(KH);
	}
	complex **  matrix = NULL, ** hh = NULL;
	set_matrix(matrix, dim_N, dim_N);
	set_matrix(hh, 2, dim_N);

//////////////////////////////////////////////////////////////////////////////
//...заполняем и решаем систему линейных уравнений A0, C0, A1, B1, C1, D1, KH;
	matrix[0][0] = R1;
	matrix[0][1] = (sh1-kk1*R1*ch1)/sqr(R1);
	matrix[0][2] = -R1;
	matrix[0][3] = -1./sqr(R1);
	matrix[0][4] = -(sh2-kk2*R1*ch2)/sqr(R1);
	matrix[0][5] = -(1.+kk2*R1)/(sqr(R1)*t2);
//
	matrix[1][2] = R2;
	matrix[1][3] = 1./sqr(R2);
	matrix[1][4] = (sh3-kk2*R2*ch3)/sqr(R2);
	matrix[1][5] = (1.+kk2*R2)/(sqr(R2)*t3);
	hh [1][1] = R2;
//
	matrix[2][0] = 1.;
	matrix[2][1] = -(2.*(sh1-kk1*R1*ch1)+(kk1*R1*kk1*R1)*sh1)/(R1*sqr(R1));
	matrix[2][2] = -1.;
	matrix[2][3] = 2./(R1*sqr(R1));
	matrix[2][4] = (2.*(sh2-kk2*R1*ch2)+(kk2*R1*kk2*R1)*sh2)/(R1*sqr(R1));
	matrix[2][5] = (2.*(1.+kk2*R1)+(kk2*R1*kk2*R1))/(R1*sqr(R1)*t2);
//
	matrix[3][2] = 1.;
	matrix[3][3] = -2./(R2*sqr(R2));
	matrix[3][4] = -(2.*(sh3-kk2*R2*ch3)+(kk2*R2*kk2*R2)*sh3)/(R2*sqr(R2));
	matrix[3][5] = -(2.*(1.+kk2*R2)+(kk2*R2*kk2*R2))/(R2*sqr(R2)*t3);
	hh [1][3] = 1.;
//
	matrix[4][1] = k1*(sh1-kk1*R1*ch1)/sqr(R1);
	matrix[4][4] = -k2*(sh2-kk2*R1*ch2)/sqr(R1);
	matrix[4][5] = -k2*(1.+kk2*R1)/(sqr(R1)*t2);
//
	matrix[5][0] = K1;
	matrix[5][2] = -K2;
	matrix[5][3] = (k2-K2)/(R1*sqr(R1));
//
	matrix[6][2] = K2;
	matrix[6][3] = -(k2-K2)/(R2*sqr(R2));
	matrix[6][6] = -1.;
//
	solver.pivot_init    (dim_N);
	solver.GaussJ(matrix, dim_N);
	int i, l;
	for (               i = 0; i < dim_N; i++)
	for (hh[0][i] = 0., l = 0; l < dim_N; l++) hh[0][i] += matrix[i][l]*hh[1][l];

	solver.pivot_init(dim_N);
	KH = hh[0][6];

	delete_struct(matrix);
	delete_struct(hh);

	return(KH);
}

////////////////////////////////////
//...классическая трехфазная модель;
complex CVisco2D::TakeEshelby_visco_volm_two(double ff)
{
	complex mu_M = comp(get_param(NUM_SHEAR), get_param(NUM_SHEAR+1)), K_M = comp(get_param(NUM_SHEAR+2), get_param(NUM_SHEAR+3)), ku_M = K_M+mu_M,
		mu_I = comp(get_param(NUM_SHEAR+NUM_SHIFT), get_param(NUM_SHEAR+1+NUM_SHIFT)), K_I = comp(get_param(NUM_SHEAR+2+NUM_SHIFT), get_param(NUM_SHEAR+3+NUM_SHIFT)), ku_I = K_I+mu_I;
	complex matr[4][5] = {
		{ 1.,  -1., -1.,  0.,   0. },	//...равенство функций;
		{ K_I+mu_M, -ku_M,  0., 0., 0. }, //...поверхностные силы;
		{ 0., ku_M,	 0., -1., mu_M },	//...эффективный модуль;
		{ 0.,   1.,  ff,  0.,   1. },	//...перемещения на границе эффективной области;
};

///////////////////////////////////////////////////////
//...решаем систему линейных уравнений A0, A1, A^1, KH;
	int dim_N = 4, ii[4] = { 0, 0, 0, 0 }, i, k, l, k0, l0;
	for (i = 0; i < dim_N; i++) {
		double f = 0.;
///////////////////////////////////////
//...look for position maximal element;
		for (k = 0; k < dim_N; k++)
			if (ii[k] != 1)
				for (l = 0; l < dim_N; l++)
					if (!ii[l]) {
						if (abs(matr[k][l]) >= f) f = abs(matr[k0 = k][l0 = l]);
					}
					else if (ii[l] > 1) return(0.);
					++(ii[l0]);
///////////////////////////////////////////////////////////
//...swapping row for diagonal position of maximal element;
					if (k0 != l0)
						for (l = 0; l <= dim_N; l++) swap(matr[k0][l], matr[l0][l]);
					if (matr[l0][l0] == comp(0.)) return(0.);
////////////////////////////////
//...diagonal row normalization;
					complex finv = comp(1.)/matr[l0][l0]; matr[l0][l0] = comp(1.);
					for (l = 0; l <= dim_N; l++) matr[l0][l] *= finv;
/////////////////////////////////
//...elimination all outher rows;
					for (k = 0; k < dim_N; k++)
						if (k != l0) {
							finv = matr[k][l0]; matr[k][l0] = comp(0.);
							for (l = 0; l <= dim_N; l++) matr[k][l] -= matr[l0][l]*finv;
						}
	}
	return(matr[3][4]);
}

/////////////////////////////////////////////////////////////////////////////////////
//...вычисление эффективного модуля сдвига через определитель (алгоритм Кристенсена);
inline complex DETER3(complex matr[8][8], int i1, int i2, int i3, int j1, int j2, int j3)
{	return matr[i1-1][j1-1]*matr[i2-1][j2-1]*matr[i3-1][j3-1]+
			 matr[i1-1][j2-1]*matr[i2-1][j3-1]*matr[i3-1][j1-1]+ 
			 matr[i2-1][j1-1]*matr[i3-1][j2-1]*matr[i1-1][j3-1]- 
			 matr[i3-1][j1-1]*matr[i2-1][j2-1]*matr[i1-1][j3-1]- 
			 matr[i2-1][j1-1]*matr[i1-1][j2-1]*matr[i3-1][j3-1]- 
			 matr[i3-1][j2-1]*matr[i2-1][j3-1]*matr[i1-1][j1-1];
}
////////////////////////////////////////////////////////////////////////////////////////////////
//...модуль сдвига в поперечной плоскости для цилиндрического включения с классической матрицей;
complex CVisco2D::TakeEshelby_visco_shear_two(double ff, double eps, int max_iter)
{
	complex mu_M = comp(get_param(NUM_SHEAR), get_param(NUM_SHEAR+1)), 
				K_M = comp(get_param(NUM_SHEAR+2), get_param(NUM_SHEAR+3)), ku_M = K_M+mu_M,
			  mu_I = comp(get_param(NUM_SHEAR+NUM_SHIFT), get_param(NUM_SHEAR+1+NUM_SHIFT)), 
				K_I = comp(get_param(NUM_SHEAR+2+NUM_SHIFT), get_param(NUM_SHEAR+3+NUM_SHIFT)), ku_I = K_I+mu_I;
////////////////////////////////////////////////////////////////////////////
//...заполняем матрицу для определения модуля сдвига в поперечной плоскости;
	complex matr[8][8] ={
//...равенство функций;
		{ 1., -2.*(K_I-mu_I)/K_I, -1., 4.*(K_M-mu_M)/K_M, -(K_M+mu_M)/(2.*K_M), -2., 0., 0. },
		{ 1., -2.*(2.*K_I+mu_I)/K_I, -1., 2.*(2.*K_M+mu_M)/K_M, -mu_M/(2.*K_M),  2., 0., 0. },
//...поверхностные силы;
		{ mu_I, 0., -mu_M, 0., mu_M, 6.*mu_M, 0., 0. },
		{ mu_I, -6.*mu_I, -mu_M, 6.*mu_M, -mu_M*.5, -6.*mu_M, 0., 0. },
//...перемещения на границе эффективной области;
		{ 0., 0., 1., -4.*(K_M-mu_M)/(K_M*ff), (K_M+mu_M)/(2.*K_M)*ff, 2.*sqr(ff), -2., -1. },
		{ 0., 0., 1., -2.*(2.*K_M+mu_M)/(K_M*ff),   mu_M/(2.*K_M)*ff, -2.*sqr(ff),  2., -1. },
//...поверхностные силы на границе эффективной области;
		{ 0., 0., mu_M, 0., -mu_M*ff, -6.*mu_M*sqr(ff), 6., -1. },
		{ 0., 0., mu_M, -6.*mu_M/ff, mu_M*.5*ff, 6.*mu_M*sqr(ff), -6., -1. },
////...другая нормировка - равенство функций;
//		{ 1., -2.*(K_I-mu_I)/K_I*ff, -1., 4.*(K_M-mu_M)/K_M*ff, -(K_M+mu_M)/(2.*K_M*ff), -2./sqr(ff), 0., 0. },
//		{ 1., -2.*(2.*K_I+mu_I)/K_I*ff, -1., 2.*(2.*K_M+mu_M)/K_M*ff, -mu_M/(2.*K_M*ff),  2./sqr(ff), 0., 0. },
////...поверхностные силы;
//		{ mu_I, 0., -mu_M, 0., mu_M/ff, 6.*mu_M/sqr(ff), 0., 0. },
//		{ mu_I, -6.*mu_I*ff, -mu_M, 6.*mu_M*ff, -mu_M*.5/ff, -6.*mu_M/sqr(ff), 0., 0. },
////...перемещения на границе эффективной области;
//		{ 0., 0., 1., -4.*(K_M-mu_M)/K_M, (K_M+mu_M)/(2.*K_M), 2., -2., -1. },
//		{ 0., 0., 1., -2.*(2.*K_M+mu_M)/K_M,   mu_M/(2.*K_M), -2.,  2., -1. },
////...поверхностные силы на границе эффективной области;
//		{ 0., 0., mu_M, 0., -mu_M, -6.*mu_M, 6., -1. },
//		{ 0., 0., mu_M, -6.*mu_M, mu_M*.5, 6.*mu_M, -6., -1. },
	};
/////////////////////////////
//...вычисляем олпределители;
	complex mu_H1, mu_H2, mu_H, A, B, C, D;
	A = (DETER3(matr, 1, 2, 3, 1, 2, 3)*DETER3(matr, 4, 7, 8, 4, 5, 6)-
		  DETER3(matr, 1, 2, 4, 1, 2, 3)*DETER3(matr, 3, 7, 8, 4, 5, 6)+
		  DETER3(matr, 1, 2, 7, 1, 2, 3)*DETER3(matr, 3, 4, 8, 4, 5, 6)-
		  DETER3(matr, 1, 2, 8, 1, 2, 3)*DETER3(matr, 3, 4, 7, 4, 5, 6)+
		  DETER3(matr, 1, 3, 4, 1, 2, 3)*DETER3(matr, 2, 7, 8, 4, 5, 6)-
		  DETER3(matr, 1, 3, 7, 1, 2, 3)*DETER3(matr, 2, 4, 8, 4, 5, 6)+
		  DETER3(matr, 1, 3, 8, 1, 2, 3)*DETER3(matr, 2, 4, 7, 4, 5, 6)+
		  DETER3(matr, 1, 4, 7, 1, 2, 3)*DETER3(matr, 2, 3, 8, 4, 5, 6)-
		  DETER3(matr, 1, 4, 8, 1, 2, 3)*DETER3(matr, 2, 3, 7, 4, 5, 6)-
		  DETER3(matr, 2, 3, 4, 1, 2, 3)*DETER3(matr, 1, 7, 8, 4, 5, 6)+
		  DETER3(matr, 2, 3, 7, 1, 2, 3)*DETER3(matr, 1, 4, 8, 4, 5, 6)-
		  DETER3(matr, 2, 3, 8, 1, 2, 3)*DETER3(matr, 1, 4, 7, 4, 5, 6)-
		  DETER3(matr, 2, 4, 7, 1, 2, 3)*DETER3(matr, 1, 3, 8, 4, 5, 6)+
		  DETER3(matr, 2, 4, 8, 1, 2, 3)*DETER3(matr, 1, 3, 7, 4, 5, 6)+
		  DETER3(matr, 3, 4, 7, 1, 2, 3)*DETER3(matr, 1, 2, 8, 4, 5, 6)-
		  DETER3(matr, 3, 4, 8, 1, 2, 3)*DETER3(matr, 1, 2, 7, 4, 5, 6));
	B = (DETER3(matr, 1, 2, 3, 1, 2, 3)*DETER3(matr, 4, 6, 8, 4, 5, 6)-
		  DETER3(matr, 1, 2, 4, 1, 2, 3)*DETER3(matr, 3, 6, 8, 4, 5, 6)+
		  DETER3(matr, 1, 2, 6, 1, 2, 3)*DETER3(matr, 3, 4, 8, 4, 5, 6)-
		  DETER3(matr, 1, 2, 8, 1, 2, 3)*DETER3(matr, 3, 4, 6, 4, 5, 6)+
		  DETER3(matr, 1, 3, 4, 1, 2, 3)*DETER3(matr, 2, 6, 8, 4, 5, 6)-
		  DETER3(matr, 1, 3, 6, 1, 2, 3)*DETER3(matr, 2, 4, 8, 4, 5, 6)+
		  DETER3(matr, 1, 3, 8, 1, 2, 3)*DETER3(matr, 2, 4, 6, 4, 5, 6)+
		  DETER3(matr, 1, 4, 6, 1, 2, 3)*DETER3(matr, 2, 3, 8, 4, 5, 6)-
		  DETER3(matr, 1, 4, 8, 1, 2, 3)*DETER3(matr, 2, 3, 6, 4, 5, 6)-
		  DETER3(matr, 2, 3, 4, 1, 2, 3)*DETER3(matr, 1, 6, 8, 4, 5, 6)+
		  DETER3(matr, 2, 3, 6, 1, 2, 3)*DETER3(matr, 1, 4, 8, 4, 5, 6)-
		  DETER3(matr, 2, 3, 8, 1, 2, 3)*DETER3(matr, 1, 4, 6, 4, 5, 6)-
		  DETER3(matr, 2, 4, 6, 1, 2, 3)*DETER3(matr, 1, 3, 8, 4, 5, 6)+
		  DETER3(matr, 2, 4, 8, 1, 2, 3)*DETER3(matr, 1, 3, 6, 4, 5, 6)+
		  DETER3(matr, 3, 4, 6, 1, 2, 3)*DETER3(matr, 1, 2, 8, 4, 5, 6)-
		  DETER3(matr, 3, 4, 8, 1, 2, 3)*DETER3(matr, 1, 2, 6, 4, 5, 6))*(-1.)+
		 (DETER3(matr, 1, 2, 3, 1, 2, 3)*DETER3(matr, 4, 6, 7, 4, 5, 6)-
		  DETER3(matr, 1, 2, 4, 1, 2, 3)*DETER3(matr, 3, 6, 7, 4, 5, 6)+
		  DETER3(matr, 1, 2, 6, 1, 2, 3)*DETER3(matr, 3, 4, 7, 4, 5, 6)-
		  DETER3(matr, 1, 2, 7, 1, 2, 3)*DETER3(matr, 3, 4, 6, 4, 5, 6)+
		  DETER3(matr, 1, 3, 4, 1, 2, 3)*DETER3(matr, 2, 6, 7, 4, 5, 6)-
		  DETER3(matr, 1, 3, 6, 1, 2, 3)*DETER3(matr, 2, 4, 7, 4, 5, 6)+
		  DETER3(matr, 1, 3, 7, 1, 2, 3)*DETER3(matr, 2, 4, 6, 4, 5, 6)+
		  DETER3(matr, 1, 4, 6, 1, 2, 3)*DETER3(matr, 2, 3, 7, 4, 5, 6)-
		  DETER3(matr, 1, 4, 7, 1, 2, 3)*DETER3(matr, 2, 3, 6, 4, 5, 6)-
		  DETER3(matr, 2, 3, 4, 1, 2, 3)*DETER3(matr, 1, 6, 7, 4, 5, 6)+
		  DETER3(matr, 2, 3, 6, 1, 2, 3)*DETER3(matr, 1, 4, 7, 4, 5, 6)-
		  DETER3(matr, 2, 3, 7, 1, 2, 3)*DETER3(matr, 1, 4, 6, 4, 5, 6)-
		  DETER3(matr, 2, 4, 6, 1, 2, 3)*DETER3(matr, 1, 3, 7, 4, 5, 6)+
		  DETER3(matr, 2, 4, 7, 1, 2, 3)*DETER3(matr, 1, 3, 6, 4, 5, 6)+
		  DETER3(matr, 3, 4, 6, 1, 2, 3)*DETER3(matr, 1, 2, 7, 4, 5, 6)-
		  DETER3(matr, 3, 4, 7, 1, 2, 3)*DETER3(matr, 1, 2, 6, 4, 5, 6))*(-5.)+
		 (DETER3(matr, 1, 2, 3, 1, 2, 3)*DETER3(matr, 4, 5, 8, 4, 5, 6)-
		  DETER3(matr, 1, 2, 4, 1, 2, 3)*DETER3(matr, 3, 5, 8, 4, 5, 6)+
		  DETER3(matr, 1, 2, 5, 1, 2, 3)*DETER3(matr, 3, 4, 8, 4, 5, 6)-
		  DETER3(matr, 1, 2, 8, 1, 2, 3)*DETER3(matr, 3, 4, 5, 4, 5, 6)+
		  DETER3(matr, 1, 3, 4, 1, 2, 3)*DETER3(matr, 2, 5, 8, 4, 5, 6)-
		  DETER3(matr, 1, 3, 5, 1, 2, 3)*DETER3(matr, 2, 4, 8, 4, 5, 6)+
		  DETER3(matr, 1, 3, 8, 1, 2, 3)*DETER3(matr, 2, 4, 5, 4, 5, 6)+
		  DETER3(matr, 1, 4, 5, 1, 2, 3)*DETER3(matr, 2, 3, 8, 4, 5, 6)-
		  DETER3(matr, 1, 4, 8, 1, 2, 3)*DETER3(matr, 2, 3, 5, 4, 5, 6)-
		  DETER3(matr, 2, 3, 4, 1, 2, 3)*DETER3(matr, 1, 5, 8, 4, 5, 6)+
		  DETER3(matr, 2, 3, 5, 1, 2, 3)*DETER3(matr, 1, 4, 8, 4, 5, 6)-
		  DETER3(matr, 2, 3, 8, 1, 2, 3)*DETER3(matr, 1, 4, 5, 4, 5, 6)-
		  DETER3(matr, 2, 4, 5, 1, 2, 3)*DETER3(matr, 1, 3, 8, 4, 5, 6)+
		  DETER3(matr, 2, 4, 8, 1, 2, 3)*DETER3(matr, 1, 3, 5, 4, 5, 6)+
		  DETER3(matr, 3, 4, 5, 1, 2, 3)*DETER3(matr, 1, 2, 8, 4, 5, 6)-
		  DETER3(matr, 3, 4, 8, 1, 2, 3)*DETER3(matr, 1, 2, 5, 4, 5, 6))*10.+
		 (DETER3(matr, 1, 2, 3, 1, 2, 3)*DETER3(matr, 4, 5, 7, 4, 5, 6)-
		  DETER3(matr, 1, 2, 4, 1, 2, 3)*DETER3(matr, 3, 5, 7, 4, 5, 6)+
		  DETER3(matr, 1, 2, 5, 1, 2, 3)*DETER3(matr, 3, 4, 7, 4, 5, 6)-
		  DETER3(matr, 1, 2, 7, 1, 2, 3)*DETER3(matr, 3, 4, 5, 4, 5, 6)+
		  DETER3(matr, 1, 3, 4, 1, 2, 3)*DETER3(matr, 2, 5, 7, 4, 5, 6)-
		  DETER3(matr, 1, 3, 5, 1, 2, 3)*DETER3(matr, 2, 4, 7, 4, 5, 6)+
		  DETER3(matr, 1, 3, 7, 1, 2, 3)*DETER3(matr, 2, 4, 5, 4, 5, 6)+
		  DETER3(matr, 1, 4, 5, 1, 2, 3)*DETER3(matr, 2, 3, 7, 4, 5, 6)-
		  DETER3(matr, 1, 4, 7, 1, 2, 3)*DETER3(matr, 2, 3, 5, 4, 5, 6)-
		  DETER3(matr, 2, 3, 4, 1, 2, 3)*DETER3(matr, 1, 5, 7, 4, 5, 6)+
		  DETER3(matr, 2, 3, 5, 1, 2, 3)*DETER3(matr, 1, 4, 7, 4, 5, 6)-
		  DETER3(matr, 2, 3, 7, 1, 2, 3)*DETER3(matr, 1, 4, 5, 4, 5, 6)-
		  DETER3(matr, 2, 4, 5, 1, 2, 3)*DETER3(matr, 1, 3, 7, 4, 5, 6)+
		  DETER3(matr, 2, 4, 7, 1, 2, 3)*DETER3(matr, 1, 3, 5, 4, 5, 6)+
		  DETER3(matr, 3, 4, 5, 1, 2, 3)*DETER3(matr, 1, 2, 7, 4, 5, 6)-
		  DETER3(matr, 3, 4, 7, 1, 2, 3)*DETER3(matr, 1, 2, 5, 4, 5, 6))*10.;
	C = (DETER3(matr, 1, 2, 3, 1, 2, 3)*DETER3(matr, 4, 5, 6, 4, 5, 6)-
		  DETER3(matr, 1, 2, 4, 1, 2, 3)*DETER3(matr, 3, 5, 6, 4, 5, 6)+
		  DETER3(matr, 1, 2, 5, 1, 2, 3)*DETER3(matr, 3, 4, 6, 4, 5, 6)-
		  DETER3(matr, 1, 2, 6, 1, 2, 3)*DETER3(matr, 3, 4, 5, 4, 5, 6)+
		  DETER3(matr, 1, 3, 4, 1, 2, 3)*DETER3(matr, 2, 5, 6, 4, 5, 6)-
		  DETER3(matr, 1, 3, 5, 1, 2, 3)*DETER3(matr, 2, 4, 6, 4, 5, 6)+
		  DETER3(matr, 1, 3, 6, 1, 2, 3)*DETER3(matr, 2, 4, 5, 4, 5, 6)+
		  DETER3(matr, 1, 4, 5, 1, 2, 3)*DETER3(matr, 2, 3, 6, 4, 5, 6)-
		  DETER3(matr, 1, 4, 6, 1, 2, 3)*DETER3(matr, 2, 3, 5, 4, 5, 6)-
		  DETER3(matr, 2, 3, 4, 1, 2, 3)*DETER3(matr, 1, 5, 6, 4, 5, 6)+
		  DETER3(matr, 2, 3, 5, 1, 2, 3)*DETER3(matr, 1, 4, 6, 4, 5, 6)-
		  DETER3(matr, 2, 3, 6, 1, 2, 3)*DETER3(matr, 1, 4, 5, 4, 5, 6)-
		  DETER3(matr, 2, 4, 5, 1, 2, 3)*DETER3(matr, 1, 3, 6, 4, 5, 6)+
		  DETER3(matr, 2, 4, 6, 1, 2, 3)*DETER3(matr, 1, 3, 5, 4, 5, 6)+
		  DETER3(matr, 3, 4, 5, 1, 2, 3)*DETER3(matr, 1, 2, 6, 4, 5, 6)-
		  DETER3(matr, 3, 4, 6, 1, 2, 3)*DETER3(matr, 1, 2, 5, 4, 5, 6))*(-20.);
	D = B*B-4.*A*C;

///////////////////////////////////////////////////////
//...определяем модуль сдвига из квадратного уравнения;
	double RR = abs(D), fi = arg2(D);
	mu_H1 = (-B-polar(RR, .5*fi))/(2.*A);
	mu_H2 = (-B+polar(RR, .5*fi))/(2.*A);

	if (real(mu_H1) > 0.) mu_H = mu_H1; else mu_H = mu_H2;
	return  (mu_H);
}
#undef  Message
