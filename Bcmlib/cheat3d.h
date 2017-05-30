/*==========================================*/
/*                  CHEAT3D                 */
/*==========================================*/
#ifndef ___CHeat3D___
#define ___CHeat3D___

#include "ccomput3d.h"

#define  Message(Msg)   { printf("%s", Msg);  printf("\n");}

//////////////////////////////////////////////////////////////
//...class of blocks partition for space electrstatic problem;
template <typename T> 
class CHeat3D : public CComput3D<T> {
public:
		int size_of_param() { return(11);}
public:
		void init() {
			delete_struct(this->param);
			this->param = new_struct<Param>(size_of_param());
		}
public:
//...constructor;
		CHeat3D (int num_phase = 8) {
			this->NUM_PHASE = num_phase;
			init();
		}
protected:
static int NUM_BASIC, NUM_SHIFT, NUM_GEOMT, MAX_PHASE;
		 int block_shape_init(Block<T> & B, Num_State id_free);
//...auxilliary operations with block matrix;
		void jump1(double * P, int i, int m);
		void jump2(double * P, int i, int m);
		void jump3(double * P, int i, int m);
		void jump4(double * P, int i, int m);
		void jump4_compos (double * P, int i, int m);
		void jump1_classic(double * P, int i, int m);
//...forming block matrix elements;
		Num_State gram1    (CGrid * nd, int i, int id_local);
		Num_State gram2    (CGrid * nd, int i, int id_local);
		Num_State gram2_old(CGrid * nd, int i, int id_local);
		Num_State gram2peri(CGrid * nd, int i, int id_local);
		Num_State gram3    (CGrid * nd, int i, int id_local);
		Num_State gram3_old(CGrid * nd, int i, int id_local);
		Num_State gram4    (CGrid * nd, int i, int id_local);
		Num_State transfer1(CGrid * nd, int i, int k, int id_local);
		Num_State transfer2(CGrid * nd, int i, int k, int id_local);
		Num_State trans_esh(CGrid * nd, int i, int k, int id_local);
		Num_State transfer3(CGrid * nd, int i, int k, int id_local) { return transfer4(nd, i, k, id_local);}
		Num_State transfer4(CGrid * nd, int i, int k, int id_local);
		Num_State rigidy1  (CGrid * nd, int i, T * K);
		Num_State rigidy2  (CGrid * nd, int i, T * K);
		Num_State rigidy5  (CGrid * nd, int i, T * K);
		Num_State computing_header(Num_Comput Num);
public:
//...параметры задачи;
		void set_fasa_hmg(double K1, double K2, double K3) { set_fasa_hmg (0., 0., K1, K2, K3, 0.);}
		void set_fasa_hmg(double K1, double K2, double K3, double C) { set_fasa_hmg (0., 0., K1, K2, K3, C);}
		void set_fasa_hmg(double R1, double R2, double K3, double K1, double K2) { set_fasa_hmg (R1, R2, K3, K1, K2, 0.);}
		void set_fasa_hmg(double R1, double R2, double K3, double K1, double K2, double C);
//...результаты решения задачи;
		void GetFuncAllValues(double X, double Y, double Z, T * F, int id_block, Num_Value id_F, int id_variant = 0, int iparam = 0);
//...аналитические модели;
		double TakeEshelby(double ff, double ff_l);
		double TakeEshelby_two (double ff);
		double TakeEshelby_grad(double ff);
		double TakeEllipsoidEshelby(double ff, double eps, double phi_stream = 0., int NX = 20, int NY = 20);
};

/////////////////////////////////////////////////
//...parametrization of the sample (preliminary);
#define DRAFT_N                     2    //...degree of multipoles;
#undef  DRAFT_N                          //...param(0);
#define DRAFT_Q                     0.92 //...normalization coefficient;
#undef  DRAFT_Q                          //...param(1);
#define DRAFT_N_quad                0.   //...parameters of quadrature;
#undef  DRAFT_N_quad                     //...param(2);
#define DRAFT_Q_facet               0.   //...normalization coefficient in facet subdivision;
#undef  DRAFT_Q_facet                    //...param(3);
#define DRAFT_C							0.   //...параметр когезионного поля;
#undef  DRAFT_C									//...param(4);
#define DRAFT_R1							1.	//...geometry of inclusion;
#undef  DRAFT_R1									//...param(5);
#define DRAFT_R2							1.	//...geometry of intermediate layer;
#undef  DRAFT_R2									//...param(6);
#define DRAFT_K3							1.   //...heat conduction in matrix;
#undef  DRAFT_K3									//...param(7);
#define DRAFT_K1							1.   //...heat conduction in inclusion;
#undef  DRAFT_K1									//...param(8);
#define DRAFT_K2							1.   //...heat conduction in intermediate;
#undef  DRAFT_K2									//...param(9);
#define DRAFT_lagrange              0.   //...Lagrange coefficient for LSM in block functional;
#undef  DRAFT_lagrange                   //...param(10);

/////////////////////////////////////////////////////
//        TEMPLATE VARIANT OF CHeat3D CLASS        //
/////////////////////////////////////////////////////
template <typename T> int CHeat3D<T>::NUM_GEOMT = 5;
template <typename T> int CHeat3D<T>::NUM_BASIC = 7;
template <typename T> int CHeat3D<T>::NUM_SHIFT = 1;
template <typename T> int CHeat3D<T>::MAX_PHASE = 3;

//////////////////////////////////
//...initialization of the blocks;
template <typename T>
int CHeat3D<T>::block_shape_init(Block<T> & B, Num_State id_free)
{
	int m;
	if (  B.shape && id_free == INITIAL_STATE) delete_shapes(B.shape);
   if (! B.shape && B.mp) {
		B.shape = new CShapeMixer<T>;
		if ((B.type & ERR_CODE) == ZOOM_BLOCK && B.mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
			B.shape->add_shape(CreateShape<T>(MP3D_POLY_SHAPE), 1);
			B.shape->add_shape(CreateShape<T>(MP3D_ZOOM_SHAPE), 1);

			B.shape->degree_init1(0, UnPackInts(this->get_param(this->NUM_MPLS)), this->solver.id_norm, draft_dim(this->type()));
			B.shape->set_shape(0, this->get_param(this->NUM_MPLS+1)*fabs(B.mp[7]));

			B.shape->degree_init1(1, UnPackInts(this->get_param(this->NUM_MPLS)), this->solver.id_norm, draft_dim(this->type()));
			B.shape->set_shape(1, fabs(B.mp[8]));

			if (this->get_param(this->NUM_GEOMT-1)) { //...подключаем внешнюю систему функций;
				B.shape->add_shape(CreateShape<T>(SK3D_ZOOM_SHAPE, 1));

				B.shape->degree_init1(2, UnPackInts(this->get_param(this->NUM_MPLS), 1), this->solver.id_norm, draft_dim(this->type()));
				B.shape->set_shape(2, fabs(B.mp[8]), sqrt(fabs(this->get_param(this->NUM_GEOMT-1))));
			}
		}
		else
		if ((B.type & ERR_CODE) == ZOOM_BLOCK && B.mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
			B.shape->add_shape(CreateShape<T>(MP3D_POLY_SHAPE), 1);
			B.shape->add_shape(CreateShape<T>(MP3D_ZOOM_SHAPE), 1);

			B.shape->degree_init1(0, UnPackInts(this->get_param(this->NUM_MPLS)), this->solver.id_norm, draft_dim(this->type()));
			B.shape->set_shape(0, fabs(B.mp[7]));

			B.shape->degree_init1(1, UnPackInts(this->get_param(this->NUM_MPLS)), this->solver.id_norm, draft_dim(this->type()));
			B.shape->set_shape(1, fabs(B.mp[8]));

			if (this->get_param(this->NUM_GEOMT-1)) {
				B.shape->add_shape(CreateShape<T>(SK3D_ZOOM_SHAPE, 0)); //...внутренняя система функций;
				B.shape->add_shape(CreateShape<T>(SK3D_ZOOM_SHAPE, 1)); //...подключаем внешнюю систему функций;

				B.shape->degree_init1(2, UnPackInts(this->get_param(this->NUM_MPLS)), this->solver.id_norm, draft_dim(this->type()));
				B.shape->set_shape(2, fabs(B.mp[7]), sqrt(fabs(this->get_param(this->NUM_GEOMT-1))));

				B.shape->degree_init1(3, UnPackInts(this->get_param(this->NUM_MPLS), 1), this->solver.id_norm, draft_dim(this->type()));
				B.shape->set_shape(3, fabs(B.mp[8]), sqrt(fabs(this->get_param(this->NUM_GEOMT-1))));
			}
		}
		else {
			if ((B.type & ERR_CODE) ==		CLAYER_BLOCK) B.shape->add_shape(CreateShape<T>(MP3D_CLAYER_SHAPE), 1); else
			if ((B.type & ERR_CODE) ==		  ELLI_BLOCK) B.shape->add_shape(CreateShape<T>(MP3D_ELLI_SHAPE),	1); else
			if ((B.type & ERR_CODE) ==		  ESHE_BLOCK) B.shape->add_shape(CreateShape<T>(MP3D_SPHEROID_SHAPE), 1); else
			if ((B.type & ERR_CODE) == ESHE_ZOOM_BLOCK) B.shape->add_shape(CreateShape<T>(MP3D_ZOOM_SHAPE), 1);
			else													  B.shape->add_shape(CreateShape<T>(MP3D_POLY_SHAPE), 1);
			
			B.shape->degree_init1(0, UnPackInts(this->get_param(this->NUM_MPLS)), this->solver.id_norm, draft_dim(this->type()));
			if ((B.type & ERR_CODE) == CLAYER_BLOCK) 
				B.shape->set_shape(0, fabs(B.mp[7]), this->get_param(NUM_BASIC+NUM_SHIFT), this->get_param(NUM_BASIC+NUM_SHIFT*2), this->get_param(NUM_BASIC), this->get_param(this->NUM_GEOMT), this->get_param(this->NUM_GEOMT+1)); else 
			if ((B.type & ERR_CODE) == ESHE_BLOCK) B.shape->set_shape(fabs(B.mp[8]), this->get_param(this->NUM_GEOMT), this->get_param(this->NUM_GEOMT+1)/this->get_param(this->NUM_GEOMT)); else
			if ((B.type & ERR_CODE) == ESHE_ZOOM_BLOCK) B.shape->set_shape(fabs(B.mp[8])); else
				B.shape->set_shape(0, fabs(B.mp[7]));

			if (this->get_param(this->NUM_GEOMT-1)) { //...внутренняя система функций;
				B.shape->add_shape(CreateShape<T>(SK3D_ZOOM_SHAPE, 0));

				B.shape->degree_init1(1, UnPackInts(this->get_param(this->NUM_MPLS)), this->solver.id_norm, draft_dim(this->type()));
				B.shape->set_shape(1, fabs(B.mp[7]), sqrt(fabs(this->get_param(this->NUM_GEOMT-1))));
			}
		}
 
////////////////////////////////////////////////////
//...local system of coordinate and parametrization;
      B.shape->set_local(B.mp+1);
      B.shape->release  ();
   }

///////////////////////////////////////
//...setting parameters and potentials;
   if (B.shape && id_free != INITIAL_STATE) {
		if (id_free == SPECIAL_STATE) { //...переустановка радиуса и центра мультиполей;
			B.shape->set_local(B.mp+1);
			if ((B.type & ERR_CODE) == ZOOM_BLOCK && B.mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				B.shape->set_shape(0, this->get_param(this->NUM_MPLS+1)*fabs(B.mp[7]));
				B.shape->set_shape(1, fabs(B.mp[8]));
				B.shape->set_shape(2, fabs(B.mp[8]), sqrt(fabs(this->get_param(this->NUM_GEOMT-1))));
			}
			else
			if ((B.type & ERR_CODE) == ZOOM_BLOCK && B.mp[0] == ID_MAP(2, SPHEROID_GENUS)) {
				B.shape->set_shape(0, fabs(B.mp[7]));
				B.shape->set_shape(1, fabs(B.mp[8]));
				B.shape->set_shape(2, fabs(B.mp[7]), sqrt(fabs(this->get_param(this->NUM_GEOMT-1))));
				B.shape->set_shape(3, fabs(B.mp[8]), sqrt(fabs(this->get_param(this->NUM_GEOMT-1))));
			}
			else {
				if ((B.type & ERR_CODE) == CLAYER_BLOCK) 
					B.shape->set_shape(0, fabs(B.mp[7]), this->get_param(NUM_BASIC+NUM_SHIFT), this->get_param(NUM_BASIC+NUM_SHIFT*2), this->get_param(NUM_BASIC), this->get_param(this->NUM_GEOMT), this->get_param(this->NUM_GEOMT+1)); else 
				if ((B.type & ERR_CODE) == ESHE_BLOCK) B.shape->set_shape(fabs(B.mp[8]), this->get_param(this->NUM_GEOMT), this->get_param(this->NUM_GEOMT+1)/this->get_param(this->NUM_GEOMT)); else
				if ((B.type & ERR_CODE) == ESHE_ZOOM_BLOCK) B.shape->set_shape(fabs(B.mp[8])); else
					B.shape->set_shape(0, fabs(B.mp[7]));
				B.shape->set_shape(1, fabs(B.mp[7]), sqrt(fabs(this->get_param(this->NUM_GEOMT-1))));
			}
		}
		else
		if (id_free == OK_STATE) {
			for (m = 0; m < this->solver.id_norm; m++) {
				B.shape->set_potential(this->solver.hh[(int)(&B-B.B)][0][m], m);
			}
		}
		else
		if (id_free == NO_STATE) {
			for (m = 0; m < this->solver.id_norm; m++) {
				B.shape->get_potential(this->solver.hh[(int)(&B-B.B)][0][this->solver.id_norm+m], m);
			}
		}
		else
		if (id_free == NULL_STATE) //...переустановка потенциалов (в случае перемены степени, например);
				B.shape->init_potential();
   }                    
   return(B.shape != NULL);
}

//////////////////////////////
//...realization of potential;
template <typename T>
void CHeat3D<T>::jump1(double * P, int i, int m)
{
	m += this->solver.id_norm;
	this->B[i].shape->cpy(this->solver.hh[i][0][m]);
}

template <typename T>
void CHeat3D<T>::jump1_classic(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(); m += this->solver.id_norm;
	this->B[i].shape->cpy(this->solver.hh[i][0][m]);
	for (l = nm; l < this->B[i].shape->num_shape(); l++) //...обнуляем когезионное поле;
		this->B[i].shape->admittance(l, this->solver.hh[i][0][m], NULL, 0., 0.);
}

//////////////////////////////////////
//...realization of normal derivative;
template <typename T>
void CHeat3D<T>::jump2(double * P, int i, int m)
{
	m += this->solver.id_norm;
	this->B[i].shape->admittance(this->solver.hh[i][0][m], NULL, 0., 0.);
	this->B[i].shape->adm_x(this->solver.hh[i][0][m], P[3]);
	this->B[i].shape->adm_y(this->solver.hh[i][0][m], P[4]);
	this->B[i].shape->adm_z(this->solver.hh[i][0][m], P[5]);
}

//////////////////////////////////
//...realization of cohesion term;
template <typename T>
void CHeat3D<T>::jump3(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(), shift = (-this->B[i].link [this->NUM_PHASE]-1)*NUM_SHIFT;
	double lambda = this->get_param(NUM_BASIC+shift); m += this->solver.id_norm;

	this->B[i].shape->cpy		 (this->solver.hh[i][0][m]);
	this->B[i].shape->admittance(this->solver.hh[i][0][m], NULL, lambda, 0.);
	for (l = 0; l < nm; l++) //...обнуляем классическое поле;
		this->B[i].shape->admittance(l, this->solver.hh[i][0][m], NULL, 0., 0.);
}

///////////////////////////////////
//...realization of heat intensity;
template <typename T>
void CHeat3D<T>::jump4(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(), shift = (-this->B[i].link [this->NUM_PHASE]-1)*NUM_SHIFT;
	double lambda = this->get_param(NUM_BASIC+shift); m += this->solver.id_norm;

	this->B[i].shape->admittance(this->solver.hh[i][0][m], NULL, 0., 0.);
	this->B[i].shape->adm_x(this->solver.hh[i][0][m], -P[3]*lambda);
	this->B[i].shape->adm_y(this->solver.hh[i][0][m], -P[4]*lambda);
	this->B[i].shape->adm_z(this->solver.hh[i][0][m], -P[5]*lambda);
	for (l = nm; l < this->B[i].shape->num_shape(); l++) //...обнуляем когезионное поле;
		this->B[i].shape->admittance(l, this->solver.hh[i][0][m], NULL, 0., 0.);
}

////////////////////////////////////////////////////////
//...realization of heat intensity for clayer functions;
template <typename T>
void CHeat3D<T>::jump4_compos(double * P, int i, int m)
{
	int l, nm = this->B[i].shape->num_usual(), shift = (-this->B[i].link [this->NUM_PHASE]-1)*NUM_SHIFT;
	double lambda = 0.;	m += this->solver.id_norm;

	if (this->B[i].shape->inverse() > 0) lambda = this->B[i].shape->get_param(2); else  
	if (! this->B[i].shape->inverse())   lambda = this->B[i].shape->get_param(0);	else lambda = this->B[i].shape->get_param(1);

	this->B[i].shape->admittance(this->solver.hh[i][0][m], NULL, 0., 0.);
	this->B[i].shape->adm_x(this->solver.hh[i][0][m], -P[3]*lambda);
	this->B[i].shape->adm_y(this->solver.hh[i][0][m], -P[4]*lambda);
	this->B[i].shape->adm_z(this->solver.hh[i][0][m], -P[5]*lambda);
	for (l = nm; l < this->B[i].shape->num_shape(); l++) //...обнуляем когезионное поле;
		this->B[i].shape->admittance(l, this->solver.hh[i][0][m], NULL, 0., 0.);
}

////////////////////////////////////////////////////////////////////////////
//...realization of common jump boundary condition for matrix and inclusion;
template <typename T>
Num_State CHeat3D<T>::gram1(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < this->solver.N && this->B[i].shape && this->B[i].mp) {
		double hh, p4, f, P[6];
		int m = this->solver.id_norm;

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (this->solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////////////////////////
//...realization of pattern cell boundary condition for Gram matrix;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = nd->nX[l];
			P[1] = nd->Y[l];  P[4] = nd->nY[l];
			P[2] = nd->Z[l];  P[5] = nd->nZ[l];
			this->B[i].shape->make_local(P);
			this->B[i].shape->norm_local(P+3);

			hh = nd->get_param(0, l); 
			p4 = nd->get_param(3, l); if (p4 == 0.) hh = 0.;
			f  = nd->get_param(4, l);

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < this->solver.n; num++)
				memset(this->solver.hh[i][0][num], 0, this->solver.dim[i]*sizeof(T));

/////////////////////////
//...jump of all moments;
			this->B[i].shape->parametrization_grad(P);

			if (p4 == MIN_HIT || p4 == MAX_HIT) {
				jump1(P, i, 0);
			}
			else
			if (0. <= p4 && p4 <= (double)(NUMS_BND-SPECIAL_BND)) {
				jump4(P, i, 0);
			}

////////////////////////////
//...composition functional;
			this->solver.to_equationDD(i, this->solver.hh[i][0][m], this->solver.hh[i][0][m], f);
			if (fabs(hh) > EE)
			  this->solver.to_equationHH(i, 0, this->solver.hh[i][0][m], hh*f);
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////////
//...junction data of the periodic boundary condition for all blocks;
template <typename T>
Num_State CHeat3D<T>::gram2(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < this->solver.N && this->B[i].shape && this->B[i].mp) {
		double AX, AY, AZ, f, P[6], 
				 g1 = .5, f1 = 1., g2 = .5, g0 = -1., TX, TY, TZ, hh;
      int id_isolated = 0, m = this->solver.id_norm, id_dir, k, j;
		if (id_isolated) {
			g0 = g1 = 1.;
			f1 = g2 = 0.;
		}

/////////////////////////////////////
//...тестовая печать множества узлов;
		 if (this->solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			AX = nd->get_param(0, l); 
			AY = nd->get_param(1, l);
			AZ = nd->get_param(2, l);
			f  = nd->get_param(4, l);

			for (k  = (int)nd->get_param(5, l), j = 0; j < this->solver.JR[i][0]; j++) 
			if ( k == this->solver.JR[i][j+this->solver.JR_SHIFT]) {
				P[0] = nd->X[l]; P[3] = nd->nX[l];
				P[1] = nd->Y[l]; P[4] = nd->nY[l];
				P[2] = nd->Z[l]; P[5] = nd->nZ[l];
				this->B[i].shape->make_local(P);
				this->B[i].shape->norm_local(P+3);

////////////////////////////////////////////////////////////////////////////////////////
//...вычисляем граничные условия периодического скачка и сдвигаем соответственные блоки;
				TX = TY = TZ = hh = 0.;
				switch (abs(id_dir = (int)nd->get_param(3, l))) {
					case 1: TX =  AX; break;
					case 2: TX = -AX; break;
					case 3: TY =  AY; break;
					case 4: TY = -AY; break;
					case 5: TZ =  AZ; hh = -AZ; break;
					case 6: TZ = -AZ; hh =  AZ; break;
				}
				this->B[k].mp[1] -= TX;
				this->B[k].mp[2] -= TY;
				this->B[k].mp[3] -= TZ; this->B[k].shape->set_local_P0(this->B[k].mp+1);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m; num < this->solver.n; num++) {
					 memset(this->solver.hh[i][0][num], 0, this->solver.dim[i]*sizeof(T));
					 memset(this->solver.hh[k][0][num], 0, this->solver.dim[k]*sizeof(T));
				}

/////////////////////////
//...jump of all moments;
				this->B[i].shape->parametrization_grad(P);
				jump1(P, i, 0); this->solver.admittance (i, 2, 0., 0, 1.);
				jump4(P, i, 1); 
				this->solver.admittance(i, 0, g1, 1, g2); 
				this->solver.admittance(i, 1, g0, 0, f1); 

				this->B[i].shape->make_common(P);
				this->B[i].shape->norm_common(P+3);

				this->B[k].shape->make_local(P);
				this->B[k].shape->norm_local(P+3);

				this->B[k].shape->parametrization_grad(P);
				jump1(P, k, 0); this->solver.admittance (k, 2, 0., 0, 1.);
				jump4(P, k, 1); 
				this->solver.admittance(k, 0, g1, 1, g2); 
				this->solver.admittance(k, 1, g0, 0, f1); 

////////////////////////////
//...composition functional;
				this->solver.to_transferTR(i, j, this->solver.hh[i][0][m],   this->solver.hh[k][0][m], f);
				this->solver.to_transferDD(i, j, this->solver.hh[i][0][m],   this->solver.hh[k][0][m+1], f);
				this->solver.to_transferTL(i, j, this->solver.hh[i][0][m+1], this->solver.hh[k][0][m+1], f);

				if (fabs(hh) > EE) {
				  this->solver.to_equationHH(i, 0, this->solver.hh[i][0][m  ],  g1*hh*f);
				  this->solver.to_equationHL(k, 0, this->solver.hh[k][0][m+1], -g2*hh*f);
				}
				if (this->solver.mode(REGUL_BOUNDARY) && (id_dir == 5 || id_dir == 6)) {//...регуляризация матрицы через граничное условие;
					this->solver.to_transferDD(i, j, this->solver.hh[i][0][m+2], this->solver.hh[k][0][m+2], f);
					if (fabs(hh) > EE) {
					  if (id_dir == 6) this->solver.to_equationHH(i, 0, this->solver.hh[i][0][m+2],  hh*f);
					  if (id_dir == 5) this->solver.to_equationHL(k, 0, this->solver.hh[k][0][m+2], -hh*f);
					}
				}
				this->B[k].mp[1] += TX;
				this->B[k].mp[2] += TY;
				this->B[k].mp[3] += TZ; this->B[k].shape->set_local_P0(this->B[k].mp+1);
 			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}
template <typename T>
Num_State CHeat3D<T>::gram2_old(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < this->solver.N && this->B[i].shape && this->B[i].mp) {
		double AX, AY, AZ, f, P[6], 
				 g1 = .5, f1 = 1., g2 = .5, g0 = -1., TX, TY, TZ, hh;
      int id_isolated = 0, m = this->solver.id_norm, id_dir, k, j, first = 1, k0, j0;
		if (id_isolated) {
			g0 = g1 = 1.;
			f1 = g2 = 0.;
		}

/////////////////////////////////////
//...тестовая печать множества узлов;
		 if (this->solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			AX = nd->get_param(0, l); 
			AY = nd->get_param(1, l);
			AZ = nd->get_param(2, l);
			f  = nd->get_param(4, l);

			for (k  = (int)nd->get_param(5, l), j = 0; j < this->solver.JR[i][0]; j++) 
			if ( k == this->solver.JR[i][j+this->solver.JR_SHIFT]) {
				P[0] = nd->X[l]; P[3] = nd->nX[l];
				P[1] = nd->Y[l]; P[4] = nd->nY[l];
				P[2] = nd->Z[l]; P[5] = nd->nZ[l];
				this->B[i].shape->make_local(P);
				this->B[i].shape->norm_local(P+3);

////////////////////////////////////////////////////////////////////////////////////////
//...вычисляем граничные условия периодического скачка и сдвигаем соответственные блоки;
				TX = TY = TZ = hh = 0.;
				switch (abs(id_dir = (int)nd->get_param(3, l))) {
					case 1: TX =  AX; break;
					case 2: TX = -AX; break;
					case 3: TY =  AY; break;
					case 4: TY = -AY; break;
					case 5: TZ =  AZ; hh = -AZ; break;
					case 6: TZ = -AZ; hh =  AZ; break;
				}
				this->B[k].mp[1] -= TX;
				this->B[k].mp[2] -= TY;
				this->B[k].mp[3] -= TZ; this->B[k].shape->set_local_P0(this->B[k].mp+1);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m+1; num >= m; num--) {
					 memset(this->solver.hh[i][0][num], 0, this->solver.dim[i]*sizeof(T));
					 memset(this->solver.hh[k][0][num], 0, this->solver.dim[k]*sizeof(T));
				}
				if (first && this->solver.mode(REGULARIZATION) && (id_dir == 5 || id_dir == 6)) {
					 first = 0; k0 = k; j0 = j;
					 memset(this->solver.hh[i][0][m+2], 0, this->solver.dim[i]*sizeof(T));
					 memset(this->solver.hh[k][0][m+2], 0, this->solver.dim[k]*sizeof(T));
				}

/////////////////////////
//...jump of all moments;
				this->B[i].shape->parametrization_grad(P);
				jump1(P, i, 0); if (! first) this->solver.admittance (i, 2, 1., 0, f); else this->solver.admittance (i, 2, 0., 0, 1.);
				jump4(P, i, 1); 
				this->solver.admittance(i, 0, g1, 1, g2); 
				this->solver.admittance(i, 1, g0, 0, f1); 

				this->B[i].shape->make_common(P);
				this->B[i].shape->norm_common(P+3);

				this->B[k].shape->make_local(P);
				this->B[k].shape->norm_local(P+3);

				this->B[k].shape->parametrization_grad(P);
				jump1(P, k, 0); if (! first) this->solver.admittance (k, 2, 1., 0, -f); else this->solver.admittance (k, 2, 0., 0, -1.);
				jump4(P, k, 1); 
				this->solver.admittance(k, 0, g1, 1, g2); 
				this->solver.admittance(k, 1, g0, 0, f1); 

////////////////////////////
//...composition functional;
				this->solver.to_transferTR(i, j, this->solver.hh[i][0][m],   this->solver.hh[k][0][m], f);
				this->solver.to_transferDD(i, j, this->solver.hh[i][0][m],   this->solver.hh[k][0][m+1], f);
				this->solver.to_transferTL(i, j, this->solver.hh[i][0][m+1], this->solver.hh[k][0][m+1], f);

				if (fabs(hh) > EE) {
				  this->solver.to_equationHH(i, 0, this->solver.hh[i][0][m  ],  g1*hh*f);
				  this->solver.to_equationHL(k, 0, this->solver.hh[k][0][m+1], -g2*hh*f);
				}
				if (first && this->solver.mode(REGUL_BOUNDARY) && (id_dir == 5 || id_dir == 6)) {//...регуляризация матрицы через граничное условие;
					this->solver.to_transferTR(i, j, this->solver.hh[i][0][m+2], this->solver.hh[k][0][m+2], f);
					this->solver.to_transferDD(i, j, this->solver.hh[i][0][m+2], this->solver.hh[k][0][m+2], f);
					this->solver.to_transferTL(i, j, this->solver.hh[i][0][m+2], this->solver.hh[k][0][m+2], f);
				}
				this->B[k].mp[1] += TX;
				this->B[k].mp[2] += TY;
				this->B[k].mp[3] += TZ; this->B[k].shape->set_local_P0(this->B[k].mp+1);
 			}
		}
		if (! first) {//...регуляризация матрицы по интегралу первого блока;
			this->solver.clean_mode(REGULARIZATION);
			this->solver.to_transferTR(i, j0, this->solver.hh[i][0][m+2], this->solver.hh[k0][0][m+2], 1.);
			this->solver.to_transferDD(i, j0, this->solver.hh[i][0][m+2], this->solver.hh[k0][0][m+2], 1.);
			this->solver.to_transferTL(i, j0, this->solver.hh[i][0][m+2], this->solver.hh[k0][0][m+2], 1.);
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////////
//...формирование матрицы Грама для периодической задачи (один блок);
template <typename T>
Num_State CHeat3D<T>::gram2peri(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < this->solver.N && this->B[i].shape && this->B[i].mp) {
		double AX, AY, AZ, f, P[6], g1 = .5, f1 = 1., g2 = .5, g0 = -1.;
      int id_isolated = 0, m = this->solver.id_norm, id_dir;
		if (id_isolated) {
			g0 = g1 = 1.;
			f1 = g2 = 0.;
		}

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (this->solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

///////////////////////////////////
//...inclusion data in gram matrix;
		for ( int l = 0; l < nd->N; l++) if (nd->hit[l] && 
			((id_dir = (int)nd->get_param(3, l)) == 1 || id_dir == 3 || id_dir == 5)) {
			P[0] = nd->X[l]; P[3] = nd->nX[l];
			P[1] = nd->Y[l]; P[4] = nd->nY[l];
			P[2] = nd->Z[l]; P[5] = nd->nZ[l];
			this->B[i].shape->make_local(P);
			this->B[i].shape->norm_local(P+3);

			AX = nd->get_param(0, l); 
			AY = nd->get_param(1, l);
			AZ = nd->get_param(2, l);
			f  = nd->get_param(4, l);

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < this->solver.n; num++)
				 memset(this->solver.hh[i][0][num], 0, this->solver.dim[i]*sizeof(T));

////////////////////////////////////////////////////
//...вычисляем функции для формирования функционала;
			this->B[i].shape->parametrization_grad(P); 
			jump1(P, i, 0); 
			jump4(P, i, 1); 
			if (this->get_param(this->NUM_GEOMT-1)) {
				jump2(P, i, 4); 
				jump3(P, i, 5); 
			}
			this->B[i].shape->make_common(P);
			this->B[i].shape->norm_common(P+3);

			if (id_dir == 1) this->B[i].mp[1] -= AX; else
			if (id_dir == 3) this->B[i].mp[2] -= AY; else
			if (id_dir == 5) this->B[i].mp[3] -= AZ; 
			this->B[i].shape->set_local_P0(this->B[i].mp+1);
			this->B[i].shape->make_local(P);
			this->B[i].shape->norm_local(P+3);

			this->B[i].shape->parametrization_grad(P);
			jump1(P, i, 2); 
			this->solver.admittance (i, 0, 1., 2, -1.); this->solver.admittance (i, 2, 2., 0, 1.);
			jump4(P, i, 3); 
			this->solver.admittance (i, 1, 1., 3, -1.);
			if (this->get_param(this->NUM_GEOMT-1)) {
				jump2(P, i, 6); this->solver.admittance (i, 4, 1., 6, -1.);
				jump3(P, i, 7); this->solver.admittance (i, 5, 1., 7, -1.);
			}
			this->solver.admittance(i, 0, g1, 1, g2); 
			this->solver.admittance(i, 1, g0, 0, f1); 
			this->solver.admittance(i, 4, g1, 5, g2); 
			this->solver.admittance(i, 5, g0, 4, f1); 

///////////////////////////////////////////////////////////////
//...сшивка периодических условий методом наименьших квадратов;
			this->solver.to_equationDD(i, this->solver.hh[i][0][m],   this->solver.hh[i][0][m], f);
			this->solver.to_equationDD(i, this->solver.hh[i][0][m+1], this->solver.hh[i][0][m+1], f);
			if (this->get_param(this->NUM_GEOMT-1)) {
				this->solver.to_equationDD(i, this->solver.hh[i][0][m+4], this->solver.hh[i][0][m+4], f);
				this->solver.to_equationDD(i, this->solver.hh[i][0][m+5], this->solver.hh[i][0][m+5], f);
			}
			if (id_dir == 5) {//...регуляризация матрицы (всегда!) и правая часть;
				this->solver.to_equationDD(i, this->solver.hh[i][0][m+2], this->solver.hh[i][0][m+2], f);
				this->solver.to_equationHH(i, 0,	this->solver.hh[i][0][m], -g1*AZ*f);
				this->solver.to_equationHH(i, 0, this->solver.hh[i][0][m+1], -g2*AZ*f);
			}
			if (id_dir == 1) this->B[i].mp[1] += AX; else
			if (id_dir == 3) this->B[i].mp[2] += AY; else
			if (id_dir == 5) this->B[i].mp[3] += AZ; 
			this->B[i].shape->set_local_P0(this->B[i].mp+1);
      }
		return(OK_STATE);
	}
	return(ERR_STATE);
}

//////////////////////////////////////////////////////////////////
//...inclusion of the stitching data to the this->solver for all blocks;
template <typename T>
Num_State CHeat3D<T>::transfer1(CGrid * nd, int i, int k, int id_local)
{
	if (nd && nd->N && 0 <= i && i < this->solver.N) {
      double f, P[6], f0 = 1., g1 = .5, g2 = .5, g0 = -1.;
      int id_isolated = 0, m = this->solver.id_norm, j, l;
		if (id_isolated) {
			g0 = g1 = 1.;
			f0 = g2 = 0.;
		}

/////////////////////////////////////
//...тестовая печать множества узлов;
		 if (this->solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (j = 0; j < this->solver.JR[i][0]; j++) if (k == this->solver.JR[i][j+this->solver.JR_SHIFT]) {
			for (l = 0; l < nd->N; l++) if (nd->hit[l]) {
				P[0] = nd->X[l];  P[3] = nd->nX[l];
				P[1] = nd->Y[l];  P[4] = nd->nY[l];
				P[2] = nd->Z[l];  P[5] = nd->nZ[l];
				this->B[i].shape->make_local(P);
				this->B[i].shape->norm_local(P+3);
				f = nd->get_param(0, l);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m; num < this->solver.n; num++) {
					memset(this->solver.hh[i][0][num], 0, this->solver.dim[i]*sizeof(T));
					memset(this->solver.hh[k][0][num], 0, this->solver.dim[k]*sizeof(T));
				}

/////////////////////////
//...jump of all moments;
				this->B[i].shape->parametrization_grad(P);
				jump1(P, i, 0); 
				jump4(P, i, 1); 
				this->solver.admittance(i, 0, g1, 1, g2); 
				this->solver.admittance(i, 1, g0, 0, f0);

				this->B[i].shape->make_common(P);
				this->B[i].shape->norm_common(P+3);

				this->B[k].shape->make_local(P);
				this->B[k].shape->norm_local(P+3);

				this->B[k].shape->parametrization_grad(P);
				jump1(P, k, 0); 
				jump4(P, k, 1); 
				this->solver.admittance(k, 0, g1, 1, g2); 
				this->solver.admittance(k, 1, g0, 0, f0); 

////////////////////////////
//...composition functional;
				this->solver.to_transferTR(i, j, this->solver.hh[i][0][m],   this->solver.hh[k][0][m], f);
				this->solver.to_transferDD(i, j, this->solver.hh[i][0][m],   this->solver.hh[k][0][m+1], f);
				this->solver.to_transferTL(i, j, this->solver.hh[i][0][m+1], this->solver.hh[k][0][m+1], f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////
//...inclusion conjunction data to the this->solver for all blocks;
template <typename T>
Num_State CHeat3D<T>::transfer2(CGrid * nd, int i, int k, int id_local)
{
	if (nd && nd->N && 0 <= i && i < this->solver.N && this->B && this->B[i].link && this->B[i].link[0] > this->NUM_PHASE) {
      double f, P[6], f0 = -1., f1 = 1., f2 = .5, g1 = .5;
      int id_isolated = 1, id_flag = 1, m = this->solver.id_norm, j, l;
		if (id_isolated) {
			f0 = g1 = 1.;
			f1 = f2 = 0.;
			id_flag = this->B[i].link[this->NUM_PHASE] == -2;
		}

/////////////////////////////////////
//...тестовая печать множества узлов;
		 if (this->solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (j = 0; j < this->solver.JR[i][0]; j++) if (k == this->solver.JR[i][j+this->solver.JR_SHIFT]) {
			for (l = 0; l < nd->N; l++) if (nd->hit[l]) {
				P[0] = nd->X[l];  P[3] = nd->nX[l];
				P[1] = nd->Y[l];  P[4] = nd->nY[l];
				P[2] = nd->Z[l];  P[5] = nd->nZ[l];
				this->B[i].shape->make_local(P);
				this->B[i].shape->norm_local(P+3);
				f = nd->get_param(0, l);

				for (int num = m; num < this->solver.n; num++) {
					memset(this->solver.hh[i][0][num], 0, this->solver.dim[i]*sizeof(T));
					memset(this->solver.hh[k][0][num], 0, this->solver.dim[k]*sizeof(T));
				}

/////////////////////////
//...jump of all moments;
				this->B[i].shape->parametrization_grad(P);
				jump1(P, i, 0); 
				jump4(P, i, 1); 
				if (this->get_param(this->NUM_GEOMT-1)) {
					jump2(P, i, 2); 
					jump3(P, i, 3); 
				}
				this->solver.admittance(i, 0, g1, 1, f2); 
				this->solver.admittance(i, 1, f0, 0, f1); 
				this->solver.admittance(i, 2, g1, 3, f2); 
				this->solver.admittance(i, 3, f0, 2, f1); 

				this->B[i].shape->make_common(P);
				this->B[i].shape->norm_common(P+3);

				this->B[k].shape->make_local(P);
				this->B[k].shape->norm_local(P+3);

				this->B[k].shape->parametrization_grad(P);
				jump1(P, k, 0); 
				jump4(P, k, 1); 
				if (this->get_param(this->NUM_GEOMT-1)) {
					jump2(P, k, 2); 
					jump3(P, k, 3); 
				}
				this->solver.admittance(k, 0, g1, 1, f2); 
				this->solver.admittance(k, 1, f0, 0, f1); 
				this->solver.admittance(k, 2, g1, 3, f2); 
				this->solver.admittance(k, 3, f0, 2, f1); 

///////////////////////////
//...composition functional;
				if (id_flag) {
					this->solver.to_transferTR(i, j, this->solver.hh[i][0][m],   this->solver.hh[k][0][m], f);
					this->solver.to_transferDD(i, j, this->solver.hh[i][0][m],   this->solver.hh[k][0][m+1], f);
					this->solver.to_transferTL(i, j, this->solver.hh[i][0][m+1], this->solver.hh[k][0][m+1], f);
					if (this->get_param(this->NUM_GEOMT-1)) {
						this->solver.to_transferTR(i, j, this->solver.hh[i][0][m+2], this->solver.hh[k][0][m+2], f);
						this->solver.to_transferDD(i, j, this->solver.hh[i][0][m+2], this->solver.hh[k][0][m+3], f);
						this->solver.to_transferTL(i, j, this->solver.hh[i][0][m+3], this->solver.hh[k][0][m+3], f);
					}
				}
				else {
					this->solver.to_transferTR(i, j, this->solver.hh[i][0][m+1], this->solver.hh[k][0][m+1], f);
					this->solver.to_transferDD(i, j, this->solver.hh[i][0][m+1], this->solver.hh[k][0][m], f);
					this->solver.to_transferTL(i, j, this->solver.hh[i][0][m],	 this->solver.hh[k][0][m], f);
					if (this->get_param(this->NUM_GEOMT-1)) {
						this->solver.to_transferTR(i, j, this->solver.hh[i][0][m+3], this->solver.hh[k][0][m+3], f);
						this->solver.to_transferDD(i, j, this->solver.hh[i][0][m+3], this->solver.hh[k][0][m+2], f);
						this->solver.to_transferTL(i, j, this->solver.hh[i][0][m+2], this->solver.hh[k][0][m+2], f);
					}
				}
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////
//...inclusion conjunction data to the this->solver for Eselby problem;
template <typename T>
Num_State CHeat3D<T>::trans_esh(CGrid * nd, int i, int k, int id_local)
{
	if (nd && nd->N && 0 <= i && i < this->solver.N && this->B && this->B[i].link && this->B[i].link[0] > this->NUM_PHASE) {
      double f, P[6], f0 = -1., f1 = 1., f2 = .5, g1 = .5, hh, pp;
      int id_isolated = 0, id_flag = 1, m = this->solver.id_norm, j, l;
		if (id_isolated) {
			f0 = g1 = 1.;
			f1 = f2 = 0.;
			id_flag = this->B[i].link[this->NUM_PHASE] == -2;
		}

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (j = 0; j < this->solver.JR[i][0]; j++) if (k == this->solver.JR[i][j+this->solver.JR_SHIFT]) {
			for (l = 0; l < nd->N; l++) if (nd->hit[l]) {
				P[0] = nd->X[l];  P[3] = nd->nX[l];
				P[1] = nd->Y[l];  P[4] = nd->nY[l];
				P[2] = nd->Z[l];  P[5] = nd->nZ[l];
				this->B[i].shape->make_local(P);
				this->B[i].shape->norm_local(P+3);
				f = nd->get_param(0, l);

				for (int num = m; num < this->solver.n; num++) {
					memset(this->solver.hh[i][0][num], 0, this->solver.dim[i]*sizeof(T));
					memset(this->solver.hh[k][0][num], 0, this->solver.dim[k]*sizeof(T));
				}

/////////////////////////
//...jump of all moments;
				this->B[i].shape->parametrization_grad(P);
				jump1(P, i, 0); 
				jump4(P, i, 1); 
				this->solver.admittance(i, 0, g1, 1, f2); 
				this->solver.admittance(i, 1, f0, 0, f1); 

				this->B[i].shape->make_common(P);
				this->B[i].shape->norm_common(P+3);

				this->B[k].shape->make_local(P);
				this->B[k].shape->norm_local(P+3);

				this->B[k].shape->parametrization_grad(P);
				jump1(P, k, 0); 
				jump4(P, k, 1); 
				this->solver.admittance(k, 0, g1, 1, f2); 
				this->solver.admittance(k, 1, f0, 0, f1); 

///////////////////////////
//...composition functional;
				double lambda = this->get_param(NUM_BASIC);
				if (id_flag) {
					this->solver.to_transferTR(i, j, this->solver.hh[i][0][m],   this->solver.hh[k][0][m], f);
					this->solver.to_transferDD(i, j, this->solver.hh[i][0][m],   this->solver.hh[k][0][m+1], f);
					this->solver.to_transferTL(i, j, this->solver.hh[i][0][m+1], this->solver.hh[k][0][m+1], f);
				}
				else {
					this->solver.to_transferTR(i, j, this->solver.hh[i][0][m+1], this->solver.hh[k][0][m+1], f);
					this->solver.to_transferDD(i, j, this->solver.hh[i][0][m+1], this->solver.hh[k][0][m], f);
					this->solver.to_transferTL(i, j, this->solver.hh[i][0][m],	 this->solver.hh[k][0][m], f);
				}
/////////////////////////////////////
//...однородное поле -- правая часть;
				if (this->B[i].link[this->NUM_PHASE] == -2) f = -f;
				hh = -nd->Z[l]/lambda; pp = nd->nZ[l]; hh = hh*g1+pp*f2; pp = pp*f0+hh*f1;
				this->solver.to_equationHH(i, 0, this->solver.hh[i][0][m],  -hh*f);
				this->solver.to_equationHH(k, 0, this->solver.hh[k][0][m+1], pp*f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////
//...формирование матрицы Грама с учетом функционала энергии;
template <typename T>
Num_State CHeat3D<T>::gram3(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < this->solver.N && this->B[i].shape && this->B[i].mp) {
		double AX, AY, AZ, f, P[6], TX, TY, TZ, hh;
      int m = this->solver.id_norm, id_dir, k, j;

/////////////////////////////////////
//...тестовая печать множества узлов;
		 if (this->solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			AX = nd->get_param(0, l); 
			AY = nd->get_param(1, l);
			AZ = nd->get_param(2, l);
			f  = nd->get_param(4, l);

			for (k  = (int)nd->get_param(5, l), j = 0; j < this->solver.JR[i][0]; j++) 
			if ( k == this->solver.JR[i][j+this->solver.JR_SHIFT]) {
				P[0] = nd->X[l]; P[3] = nd->nX[l];
				P[1] = nd->Y[l]; P[4] = nd->nY[l];
				P[2] = nd->Z[l]; P[5] = nd->nZ[l];
				this->B[i].shape->make_local(P);
				this->B[i].shape->norm_local(P+3);

////////////////////////////////////////////////////////////////////////////////////////
//...вычисляем граничные условия периодического скачка и сдвигаем соответственные блоки;
				TX = TY = TZ = hh = 0.;
				switch (abs(id_dir = (int)nd->get_param(3, l))) {
					case 1: TX =  AX; break;
					case 2: TX = -AX; break;
					case 3: TY =  AY; break;
					case 4: TY = -AY; break;
					case 5: TZ =  AZ; hh = -AZ; break;
					case 6: TZ = -AZ; hh =  AZ; break;
				}
				this->B[k].mp[1] -= TX;
				this->B[k].mp[2] -= TY;
				this->B[k].mp[3] -= TZ; this->B[k].shape->set_local_P0(this->B[k].mp+1);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m; num < this->solver.n; num++) {
					 memset(this->solver.hh[i][0][num], 0, this->solver.dim[i]*sizeof(T));
					 memset(this->solver.hh[k][0][num], 0, this->solver.dim[k]*sizeof(T));
				}

///////////////////////////////////////////////////////////////
//...вычисляем все необходимые моменты коллокационного вектора;
				this->B[i].shape->parametrization_grad(P);
				jump1(P, i, 0);
				jump4(P, i, 1); 

				this->B[i].shape->make_common(P);
				this->B[i].shape->norm_common(P+3);

				this->B[k].shape->make_local(P);
				this->B[k].shape->norm_local(P+3);

				this->B[k].shape->parametrization_grad(P);
				jump1(P, k, 0);
				jump4(P, k, 1); 

//////////////////////////////////////////////////////////////////
//...сшивка функций и условие скачка методом наименьших квадратов;
				this->solver.to_transferTR(i, j, this->solver.hh[i][0][m], this->solver.hh[k][0][m], f);
				this->solver.to_transferDD(i, j, this->solver.hh[i][0][m], this->solver.hh[k][0][m], f);
				this->solver.to_transferTL(i, j, this->solver.hh[i][0][m], this->solver.hh[k][0][m], f);

				if (fabs(hh) > EE) {
				  this->solver.to_equationHH(i, 0, this->solver.hh[i][0][m],  hh*f);
				  this->solver.to_equationHL(k, 0, this->solver.hh[k][0][m], -hh*f);
				}
				if (this->solver.mode(REGUL_BOUNDARY) && (id_dir == 5 || id_dir == 6)) {//...регуляризация матрицы через граничное условие;
					this->solver.to_transferDD(i, j, this->solver.hh[i][0][m], this->solver.hh[k][0][m], f);
					if (fabs(hh) > EE) {
					  if (id_dir == 6) this->solver.to_equationHH(i, 0, this->solver.hh[i][0][m],  hh*f);
					  if (id_dir == 5) this->solver.to_equationHL(k, 0, this->solver.hh[k][0][m], -hh*f);
					}
				}

/////////////////////////////
//...энергетиеские слагаемые;
				this->solver.to_equationER(i, this->solver.hh[i][0][m], this->solver.hh[i][0][m+1], -f); //...поменяли знак (нормаль внутрь);
				this->solver.to_equationEL(k, this->solver.hh[k][0][m], this->solver.hh[k][0][m+1],  f);
				
				this->B[k].mp[1] += TX;
				this->B[k].mp[2] += TY;
				this->B[k].mp[3] += TZ; this->B[k].shape->set_local_P0(this->B[k].mp+1);
			}
      }
		return(OK_STATE);
	}
	return(ERR_STATE);
}
template <typename T>
Num_State CHeat3D<T>::gram3_old(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < this->solver.N && this->B[i].shape && this->B[i].mp) {
		double AX, AY, AZ, f, P[6], TX, TY, TZ, hh;
      int m = this->solver.id_norm, id_dir, k, j, first = 1, k0, j0;

/////////////////////////////////////
//...тестовая печать множества узлов;
		 if (this->solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			AX = nd->get_param(0, l); 
			AY = nd->get_param(1, l);
			AZ = nd->get_param(2, l);
			f  = nd->get_param(4, l);

			for (k  = (int)nd->get_param(5, l), j = 0; j < this->solver.JR[i][0]; j++) 
			if ( k == this->solver.JR[i][j+this->solver.JR_SHIFT]) {
				P[0] = nd->X[l]; P[3] = nd->nX[l];
				P[1] = nd->Y[l]; P[4] = nd->nY[l];
				P[2] = nd->Z[l]; P[5] = nd->nZ[l];
				this->B[i].shape->make_local(P);
				this->B[i].shape->norm_local(P+3);

////////////////////////////////////////////////////////////////////////////////////////
//...вычисляем граничные условия периодического скачка и сдвигаем соответственные блоки;
				TX = TY = TZ = hh = 0.;
				switch (abs(id_dir = (int)nd->get_param(3, l))) {
					case 1: TX =  AX; break;
					case 2: TX = -AX; break;
					case 3: TY =  AY; break;
					case 4: TY = -AY; break;
					case 5: TZ =  AZ; hh = -AZ; break;
					case 6: TZ = -AZ; hh =  AZ; break;
				}
				this->B[k].mp[1] -= TX;
				this->B[k].mp[2] -= TY;
				this->B[k].mp[3] -= TZ; this->B[k].shape->set_local_P0(this->B[k].mp+1);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m+1; num >= m; num--) {
					 memset(this->solver.hh[i][0][num], 0, this->solver.dim[i]*sizeof(T));
					 memset(this->solver.hh[k][0][num], 0, this->solver.dim[k]*sizeof(T));
				}
				if (first && this->solver.mode(REGULARIZATION) && (id_dir == 5 || id_dir == 6)) {
					 first = 0; k0 = k; j0 = j;
					 memset(this->solver.hh[i][0][m+2], 0, this->solver.dim[i]*sizeof(T));
					 memset(this->solver.hh[k][0][m+2], 0, this->solver.dim[k]*sizeof(T));
				}

///////////////////////////////////////////////////////////////
//...вычисляем все необходимые моменты коллокационного вектора;
				this->B[i].shape->parametrization_grad(P);
				jump1(P, i, 0); if (! first) this->solver.admittance (i, 2, 1., 0, f); else this->solver.admittance (i, 2, 0., 0, 1.);
				jump4(P, i, 1); 

				this->B[i].shape->make_common(P);
				this->B[i].shape->norm_common(P+3);

				this->B[k].shape->make_local(P);
				this->B[k].shape->norm_local(P+3);

				this->B[k].shape->parametrization_grad(P);
				jump1(P, k, 0); if (! first) this->solver.admittance (k, 2, 1., 0, -f); else this->solver.admittance (k, 2, 0., 0, -1.);
				jump4(P, k, 1); 

//////////////////////////////////////////////////////////////////
//...сшивка функций и условие скачка методом наименьших квадратов;
				this->solver.to_transferTR(i, j, this->solver.hh[i][0][m], this->solver.hh[k][0][m], f);
				this->solver.to_transferDD(i, j, this->solver.hh[i][0][m], this->solver.hh[k][0][m], f);
				this->solver.to_transferTL(i, j, this->solver.hh[i][0][m], this->solver.hh[k][0][m], f);

				if (fabs(hh) > EE) {
				  this->solver.to_equationHH(i, 0, this->solver.hh[i][0][m],  hh*f);
				  this->solver.to_equationHL(k, 0, this->solver.hh[k][0][m], -hh*f);
				}
				if (first && this->solver.mode(REGUL_BOUNDARY) && (id_dir == 5 || id_dir == 6)) {//...регуляризация матрицы через граничное условие;
					this->solver.to_transferTR(i, j, this->solver.hh[i][0][m+2], this->solver.hh[k][0][m+2], f);
					this->solver.to_transferDD(i, j, this->solver.hh[i][0][m+2], this->solver.hh[k][0][m+2], f);
					this->solver.to_transferTL(i, j, this->solver.hh[i][0][m+2], this->solver.hh[k][0][m+2], f);
				}

/////////////////////////////
//...энергетиеские слагаемые;
				this->solver.to_equationER(i, this->solver.hh[i][0][m], this->solver.hh[i][0][m+1],  f);
				this->solver.to_equationEL(k, this->solver.hh[k][0][m], this->solver.hh[k][0][m+1], -f);
				
				this->B[k].mp[1] += TX;
				this->B[k].mp[2] += TY;
				this->B[k].mp[3] += TZ; this->B[k].shape->set_local_P0(this->B[k].mp+1);
			}
      }
		if (! first) {//...регуляризация матрицы по интегралу первого блока;
			this->solver.clean_mode(REGULARIZATION);
			this->solver.to_transferTR(i, j0, this->solver.hh[i][0][m+2], this->solver.hh[k0][0][m+2], 1.);
			this->solver.to_transferDD(i, j0, this->solver.hh[i][0][m+2], this->solver.hh[k0][0][m+2], 1.);
			this->solver.to_transferTL(i, j0, this->solver.hh[i][0][m+2], this->solver.hh[k0][0][m+2], 1.);
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////
//...формирование матрицы Грама с учетом функционала энергии;
template <typename T>
Num_State CHeat3D<T>::gram4(CGrid * nd, int i, int id_local)
{
	if (nd && nd->N && 0 <= i && i < this->solver.N && this->B[i].shape && this->B[i].mp) {
		double hh, p4, f, P[6];
		int m = this->solver.id_norm;

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (this->solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.002, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////////////////////////
//...realization of pattern cell boundary condition for Gram matrix;
		for (int l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = nd->nX[l];
			P[1] = nd->Y[l];  P[4] = nd->nY[l];
			P[2] = nd->Z[l];  P[5] = nd->nZ[l];
			this->B[i].shape->make_local(P);
			this->B[i].shape->norm_local(P+3);

			hh = nd->get_param(0, l); 
			p4 = nd->get_param(3, l); if (p4 == 0.) hh = 0.;
			f  = nd->get_param(4, l);

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < this->solver.n; num++)
				memset(this->solver.hh[i][0][num], 0, this->solver.dim[i]*sizeof(T));

///////////////////////////////////////////////////////////////
//...вычисляем все необходимые моменты коллокационного вектора;
			this->B[i].shape->parametrization_grad(P);
			jump1(P, i, 0); 
			jump4(P, i, 1);

////////////////////////////////////////////////////
//...граничное условие методом наименьших квадратов;
			if (p4 == MIN_HIT || p4 >= NUMS_BND && p4 < NUMS_BND+20 || p4 == MAX_HIT) {
				this->solver.to_equationDD(i, this->solver.hh[i][0][m], this->solver.hh[i][0][m], f);
				if (fabs(hh) > EE)
					this->solver.to_equationHH(i, 0, this->solver.hh[i][0][m], hh*f);
			}

//////////////////////////////
//...энергетические слагаемые;
			this->solver.to_equationEE(i, this->solver.hh[i][0][m], this->solver.hh[i][0][m+1], f);
			if (1. <= p4 && p4 <= (double)(NUMS_BND-SPECIAL_BND) && fabs(hh) > EE)
				this->solver.to_equationEH(i, 0, this->solver.hh[i][0][m], hh*f);
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

///////////////////////////////////////////////////////////////
//...формирование матриц перехода с учетом функционала энергии;
template <typename T>
Num_State CHeat3D<T>::transfer4(CGrid * nd, int i, int k, int id_local)
{
	if (nd && nd->N && 0 <= i && i < this->solver.N) {
      double f, P[6];
      int m = this->solver.id_norm, j, l;

/////////////////////////////////////
//...тестовая печать множества узлов;
		 if (this->solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.02, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for (j = 0; j < this->solver.JR[i][0]; j++) if (k == this->solver.JR[i][j+this->solver.JR_SHIFT]) {
			for (l = 0; l < nd->N; l++) if (nd->hit[l]) {
				P[0] = nd->X[l];  P[3] = nd->nX[l];
				P[1] = nd->Y[l];  P[4] = nd->nY[l];
				P[2] = nd->Z[l];  P[5] = nd->nZ[l];
				this->B[i].shape->make_local(P);
				this->B[i].shape->norm_local(P+3);
				f = nd->get_param(0, l);

/////////////////////////////
//...reset auxilliary arrays;
				for (int num = m; num < this->solver.n; num++) {
					memset(this->solver.hh[i][0][num], 0, this->solver.dim[i]*sizeof(T));
					memset(this->solver.hh[k][0][num], 0, this->solver.dim[k]*sizeof(T));
				}

///////////////////////////////////////////////////////////////
//...вычисляем все необходимые моменты коллокационного вектора;
				this->B[i].shape->parametrization_grad(P);
				jump1(P, i, 0); 
				jump4(P, i, 1); 

				this->B[i].shape->make_common(P);
				this->B[i].shape->norm_common(P+3);

				this->B[k].shape->make_local(P);
				this->B[k].shape->norm_local(P+3);

				this->B[k].shape->parametrization_grad(P);
				jump1(P, k, 0); 
				jump4(P, k, 1); 

/////////////////////////////////////////////////
//...сшивка функций методом наименьших квадратов;
				this->solver.to_transferTR(i, j, this->solver.hh[i][0][m], this->solver.hh[k][0][m], f);
				this->solver.to_transferDD(i, j, this->solver.hh[i][0][m], this->solver.hh[k][0][m], f);
				this->solver.to_transferTL(i, j, this->solver.hh[i][0][m], this->solver.hh[k][0][m], f);

/////////////////////////////
//...энергетиеские слагаемые;
				this->solver.to_equationER(i, this->solver.hh[i][0][m], this->solver.hh[i][0][m+1], -f); //...поменяли знак (нормаль внутрь);
				this->solver.to_equationEL(k, this->solver.hh[k][0][m], this->solver.hh[k][0][m+1],  f);
			}
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

////////////////////////////////////////////////////////
//...интегрирование параметров потока по границе блоков;
template <typename T>
Num_State CHeat3D<T>::rigidy1(CGrid * nd, int i, T * K)
{
	if (nd) {
      int shift = (-this->B[i].link[this->NUM_PHASE]-1)*NUM_SHIFT, m = this->solver.id_norm, l;
      double kk = this->get_param(NUM_BASIC+shift), f, P[6]; T TH;

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (this->solver.mode(PRINT_MODE)) 
			nd->TestGrid("nodes.bln", 0.002, 10., 20., 30., AXIS_Z, 1);

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for ( l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = nd->nX[l];
			P[1] = nd->Y[l];  P[4] = nd->nY[l];
			P[2] = nd->Z[l];  P[5] = nd->nZ[l];
			this->B[i].shape->make_local(P);
			this->B[i].shape->norm_local(P+3);
			f = nd->get_param(0, l);

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < this->solver.n; num++)
				memset(this->solver.hh[i][0][num], 0, this->solver.dim[i]*sizeof(double));

/////////////////////////////////////////////////////////////////////////
//...вычисляем функцию температуры для формирования интеграла от потоков;
			this->B[i].shape->parametrization(P, 1);
			jump1_classic(P, i, 0); 

			TH = this->B[i].shape->potential(this->solver.hh[i][0][m], 0);
			if ((this->B[i].type & ERR_CODE) == CLAYER_BLOCK) { 
				if (this->B[i].shape->inverse() > 0) kk = this->B[i].shape->get_param(2); else  
				if (! this->B[i].shape->inverse())   kk = this->B[i].shape->get_param(0); else kk = this->B[i].shape->get_param(1);
			}		
			K[0] += kk*TH*P[3]*f;
			K[1] += kk*TH*P[4]*f;
			K[2] += kk*TH*P[5]*f;
			K[3+(-this->B[i].link[this->NUM_PHASE]-1)] += P[2]*P[5]*f;
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

///////////////////////////////////////////////////////
//...интегрирование параметров потока по объему блоков;
template <typename T>
Num_State CHeat3D<T>::rigidy2(CGrid * nd, int i, T * K)
{
	if (nd) {
      int shift = (-this->B[i].link[this->NUM_PHASE]-1)*NUM_SHIFT, m = this->solver.id_norm, l;
      double f, P[6];T PX, PY, PZ;

////////////////////////////////////////////////
//...inclusion data in gram and transfer matrix;
		for ( l = 0; l < nd->N; l++) if (nd->hit[l]) {
			P[0] = nd->X[l];  P[3] = -1;
			P[1] = nd->Y[l];  P[4] = 0.;
			P[2] = nd->Z[l];  P[5] = 0.;
			this->B[i].shape->make_local(P);
			this->B[i].shape->norm_local(P+3);
			f = nd->get_param(0, l);

/////////////////////////////
//...reset auxilliary arrays;
			for (int num = m; num < this->solver.n; num++)
				memset(this->solver.hh[i][0][num], 0, this->solver.dim[i]*sizeof(T));

/////////////////////////////////////////////////////////////////////////
//...вычисляем функцию температуры для формирования интеграла от потоков;
			this->B[i].shape->parametrization_grad(P);
			jump4(P, i, 0);

			PX = this->B[i].shape->potential(this->solver.hh[i][0][m], 0);

			P[4] = -1.; P[3] = P[5] = 0.;
			this->B[i].shape->norm_local(P+3);
			jump4(P, i, 0);

			PY = this->B[i].shape->potential(this->solver.hh[i][0][m], 0);

			P[5] = -1.; P[3] = P[4] = 0.;
			this->B[i].shape->norm_local(P+3);
			jump4(P, i, 0);

			PZ = this->B[i].shape->potential(this->solver.hh[i][0][m], 0);
			
			K[0] += PX*f;
			K[1] += PY*f;
			K[2] += PZ*f;
			K[3+(-this->B[i].link[this->NUM_PHASE]-1)] += f;
		}
		return(OK_STATE);
	}
	return(ERR_STATE);
}

/////////////////////////////////////////////////////////////////////////////////////////
//...дополнительное интегрирование параметров потока для блоков с усложненными функциями;
template <typename T>
Num_State CHeat3D<T>::rigidy5(CGrid * nd, int i, T * K)
{
   int shift = (-this->B[i].link[this->NUM_PHASE]-1)*NUM_SHIFT, 
		N_elem = UnPackInts(this->get_param(this->NUM_QUAD)), N_max = UnPackInts(this->get_param(this->NUM_QUAD), 1), m, k, l, l0;
	double kk = this->get_param(NUM_BASIC+shift), K1 = this->get_param(NUM_BASIC+NUM_SHIFT), K2 = this->get_param(NUM_BASIC+NUM_SHIFT*2), 
			 R1 = this->get_param(this->NUM_GEOMT), R2 = this->get_param(this->NUM_GEOMT+1), L0 = to_double(K[6]), f, P[6];	T TH, sum[4];

	if (this->B[i].shape->type() == MP3D_CLAYER_SHAPE && (R1 || R2))	{
		CGrid * bound_bnd = CreateNodes(GRID_QG_NODES);
				  bound_bnd->add_params(1);

//////////////////////////////////////////
//...образуем временную поверхность сферы;
		CCells * ce = new(CCells);
		ce->init(1, 2, (l = size_of_map(2, SPHERE_GENUS))+1+size_of_dop(SPH_SEGMENT));
		ce->mp[0] = (CMap)ID_MAP(2, SPHERE_GENUS);
		ce->mp[1] = this->B[i].mp[1];
		ce->mp[2] = this->B[i].mp[2];
		ce->mp[3] = this->B[i].mp[3];
		ce->mp[7] = R1;
		ce->mp[l] = (CMap)SPH_SEGMENT;
      ce->mp[++l] = 2.*M_PI/N_max;
      ce->mp[++l] = 0.;
      ce->mp[++l] = M_PI/N_max; l0 = l;

/////////////////////////////////
//...накапливаем граничные точки;
		for (m = 0; m < N_max; m++, ce->mp[l0-1] = ce->mp[l0], ce->mp[l0] *= (m+1.)/m)
		for (k = 0; k < N_max; k++) {
			ce->mp[6] = ce->mp[l0-2]*k;
			if (L0 < R1) bound_bnd->sphere_intrusion_QG(ce->mp, N_elem, L0);
			else	ce->segms_QG(bound_bnd, N_elem);
		}

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (this->solver.mode(PRINT_MODE)) 
			bound_bnd->TestGrid("nodes1.bln", 0.002, 10., 20., 30., AXIS_X, 1);

////////////////////////////////
//...интегрирование температуры;
		for (memset(sum, 0, 4*sizeof(T)), l = 0; l < bound_bnd->N; l++) {
			P[0] = bound_bnd->X[l]; P[3] = bound_bnd->nX[l];
			P[1] = bound_bnd->Y[l]; P[4] = bound_bnd->nY[l];
			P[2] = bound_bnd->Z[l]; P[5] = bound_bnd->nZ[l];
			this->B[i].shape->make_local(P);
			this->B[i].shape->norm_local(P+3);

			memset(this->solver.hh[i][0][this->solver.id_norm], 0, this->solver.dim[i]*sizeof(T));
			this->B[i].shape->parametrization(P, 1);
			jump1(P, i, 0); 

			TH = this->B[i].shape->potential(this->solver.hh[i][0][this->solver.id_norm], 0);

			sum[0] += TH*P[3]*(f = bound_bnd->get_param(0, l));
			sum[1] += TH*P[4]* f;
			sum[2] += TH*P[5]* f;
			sum[3] += (P[0]*P[3]+P[1]*P[4]+P[2]*P[5])/3.*f;
		}
		if (to_double(sum[3]) < 0.) { //...правим знак в случае ошибки в направлении нормали;
			sum[0] = -sum[0];
			sum[1] = -sum[1];
			sum[2] = -sum[2];
			sum[3] = -sum[3];
		}
		if (L0 < R1) //...правим площадь промежуточного слоя в случае взаимопроникающих включений;
			sum[3] += (R1*R1-L0*L0)*2.*M_PI*L0;

		K[0] += (K1-K2)*sum[0];
		K[1] += (K1-K2)*sum[1];
		K[2] += (K1-K2)*sum[2];
		K[4] += sum[3];
		K[5] -= sum[3];
		if (! this->solver.mode(ACCUMULATION)) bound_bnd->add_buffer(bound_bnd->N);

//////////////////////////////////////////////////////
//...повторяем граничные точки для второй поверхности;
		ce->mp[7] = R2;
      ce->mp[l0-1] = 0.;
      ce->mp[l0] = M_PI/N_max;
		
		for (m = 0; m < N_max; m++, ce->mp[l0-1] = ce->mp[l0], ce->mp[l0] *= (m+1.)/m)
		for (k = 0; k < N_max; k++) {
			ce->mp[6] = ce->mp[l0-2]*k;
			if (L0 < R2) bound_bnd->sphere_intrusion_QG(ce->mp, N_elem, L0);
			else	ce->segms_QG(bound_bnd, N_elem);
		}

/////////////////////////////////////
//...тестовая печать множества узлов;
		if (this->solver.mode(PRINT_MODE)) 
			bound_bnd->TestGrid("nodes2.bln", 0.002, 10., 20., 30., AXIS_X, 1);

////////////////////////////////
//...интегрирование температуры;
		for (memset(sum, 0, 4*sizeof(T)), l = 0; l < bound_bnd->N; l++) {
			P[0] = bound_bnd->X[l]; P[3] = bound_bnd->nX[l];
			P[1] = bound_bnd->Y[l]; P[4] = bound_bnd->nY[l];
			P[2] = bound_bnd->Z[l]; P[5] = bound_bnd->nZ[l];
			this->B[i].shape->make_local(P);
			this->B[i].shape->norm_local(P+3);

			memset(this->solver.hh[i][0][this->solver.id_norm], 0, this->solver.dim[i]*sizeof(T));
			this->B[i].shape->parametrization(P, 1);
			jump1(P, i, 0); 

			TH = this->B[i].shape->potential(this->solver.hh[i][0][this->solver.id_norm], 0);

			sum[0] += TH*P[3]*(f = bound_bnd->get_param(0, l));
			sum[1] += TH*P[4]* f;
			sum[2] += TH*P[5]* f;
			sum[3] += (P[0]*P[3]+P[1]*P[4]+P[2]*P[5])/3.*f;
		}
		if (to_double(sum[3]) < 0.) { //...правим знак в случае ошибки в направлении нормали;
			sum[0] = -sum[0];
			sum[1] = -sum[1];
			sum[2] = -sum[2];
			sum[3] = -sum[3];
		}
		if (L0 < R2) //...правим площадь промежуточного слоя в случае взаимопроникающих включений;
			sum[3] += (R2*R2-L0*L0)*2.*M_PI*L0;

		K[0] += (K2-kk)*sum[0];
		K[1] += (K2-kk)*sum[1];
		K[2] += (K2-kk)*sum[2];
		K[5] += sum[3];
		K[3+(-this->B[i].link[this->NUM_PHASE]-1)] -= sum[3];
		if (! this->solver.mode(ACCUMULATION)) bound_bnd->add_buffer(bound_bnd->N);

		delete ce;
		delete bound_bnd;
	}
	return(OK_STATE);
}

///////////////////////////////////////
//...auxiliary function for parameters;
template <typename T>
void CHeat3D<T>::set_fasa_hmg(double R1, double R2, double K3, double K1, double K2, double C)
{
	if (size_of_param() > this->NUM_GEOMT+2) {
		this->param[this->NUM_GEOMT-1] = C;
		this->param[this->NUM_GEOMT]   = R1;
		this->param[this->NUM_GEOMT+1] = R2;
	}
	if (size_of_param() > NUM_BASIC+NUM_SHIFT*2) {
		this->param[NUM_BASIC] = K3;
		this->param[NUM_BASIC+NUM_SHIFT] = K1;
		this->param[NUM_BASIC+NUM_SHIFT*2] = K2;
	}
}

///////////////////////////////////////////////////////
//...counting header for solving electrostatic problem;
template <typename T>
Num_State CHeat3D<T>::computing_header(Num_Comput Num)
{
	int N_elem = UnPackInts(this->get_param(this->NUM_QUAD)), i, k, elem, id_dir, n_rhs = 2;
	char msg[201];
	
	if (! this->solver.mode(NO_MESSAGE)) {
		Message(" ");

		sprintf(msg, "CHeat3D sample: N_sm = %d, N_mpl = %d, N_elem = %d, Q_facet = %g", this->N, UnPackInts(this->get_param(this->NUM_MPLS)), N_elem, this->get_param(this->NUM_QUAD+1));
		Message(msg);

		Message(" ");
		switch (Num){
			case	 BASIC_COMPUT: Message("Junction Blocks..."); break;
			case MAPPING_COMPUT: Message("Mapping Blocks..."); break;
			case  PERIOD_COMPUT: Message("Periodic Block..."); break;
		}
		Message(" ");
	}

//////////////////////////////////
//...определяем блочную структуру;
	this->solver.set_blocks(this->N, n_rhs); //<==== number of saved potentials !!!
	this->solver.n += 8;//<==== number of additional auxilliary arrays!!!
	for (k = 0; k < this->solver.N;  k++)
		  this->solver.set_links(k, this->B[k].link);

	this->shapes_init(INITIAL_STATE);
	this->shapes_init(NULL_STATE);

////////////////////////////////////////////////////////////////////
//...добавляем блоки на периодических связях и на границе включений;
	if (this->solv%ENERGY_SOLVING == PERIODIC_SOLVING) { 
		double par[6]; this->SetBounding(par);
		for (k = 0; k < this->N; k++) if (this->B[k].link && this->B[k].link[this->NUM_PHASE] == -1) {
			i = 0; while ((elem =  this->geom_plink_3D(this->B[k], i, id_dir, par)) >= 0)			  
			this->solver.add_link(k, elem);
			i = 0; while ((elem = this->block_plink_3D(this->B[k], i, id_dir, par)) >= 0)			  
			this->solver.add_link(k, elem);
		}
	}
	this->LinkPhase3D(MAX_PHASE);

/////////////////////////////////////////////////////////
//...делаем перенумерацию структуры и задаем размерность;
	if (! this->solver.struct_permutat(this->solver.id_change == EXTERN_STATE ? NULL_STATE : /*NULL*/OK_STATE) || 
		 ! this->solver.inverse_index()) {
		return ERR_STATE;
	}
	for (k = 0; k < this->solver.N; k++)
		this->solver.set_dimension(k, this->freedom_block(k));

	this->solver.struct_init();

	if (this->solver.mode(FULLY_MODE)) { 
		this->solver.test_struct("1", 0);
		this->solver.test_struct("2", 1);
	}
   return OK_STATE;
}

//////////////////////////////////////////////////////////////////////////
//...calculation of function values (in common coordinate system) on grid;
template <typename T>
void CHeat3D<T>::GetFuncAllValues(double X, double Y, double Z, T * F, int i, Num_Value id_F, int id_variant, int iparam)
{
	if (! F) return;
	double P[6]  = { X, Y, Z, 1., 0., 0.};

/////////////////////////////////////
//...operation with all input points;
	if ( 0 <= i && i < this->N && this->B[i].shape && this->B[i].mp) {
		int m = this->solver.id_norm;

//////////////////////////////////////////////////////
//...reset auxilliary arrays and calculation function;
		for (int num = m; num < this->solver.n; num++)
			memset(this->solver.hh[i][0][num], 0, this->solver.dim[i]*sizeof(T));

		this->B[i].shape->make_local(P);
		this->B[i].shape->norm_local(P+3);
		switch (id_F) {
			case HEAT_VALUE: {
///////////////////////////
//...calculation potential;
				this->B[i].shape->parametrization(P, 1);
				jump1(P, i, 0); 

				F[0] = F[1] = this->B[i].shape->potential(this->solver.hh[i][0][m], id_variant); 
				if (this->solv == PERIODIC_SOLVING || this->solv == E_PERIODIC_SOLVING) F[1] -= Z;
			}	break;
			case HEAT_ESHE_VALUE: {
///////////////////////////
//...calculation potential;
				this->B[i].shape->parametrization(P, 1);
				jump1(P, i, 0); 

				F[0] = F[1] = this->B[i].shape->potential(this->solver.hh[i][0][m], id_variant); 
				if (this->B[i].link[this->NUM_PHASE] == -1) {
					double lambda = this->get_param(NUM_BASIC);
					F[0] = (F[1] -= Z/lambda);
				}
				if (this->solv ==	PERIODIC_SOLVING) F[1] -= Z;
			}	break;
			case FLUX_VALUE: {
///////////////////////////
//...calculation heat flux;
				P[3] = -1.; P[4] = P[5] = 0.;
				this->B[i].shape->norm_local(P+3);
				this->B[i].shape->parametrization_grad(P);
				jump4(P, i, 0);

				F[0] = this->B[i].shape->potential(this->solver.hh[i][0][m], id_variant);

				P[4] = -1.; P[3] = P[5] = 0.;
				this->B[i].shape->norm_local(P+3);
				jump4(P, i, 0);

				F[1] = this->B[i].shape->potential(this->solver.hh[i][0][m], id_variant);

				P[5] = -1.; P[3] = P[4] = 0.;
				this->B[i].shape->norm_local(P+3);
				jump4(P, i, 0);

				F[2] = this->B[i].shape->potential(this->solver.hh[i][0][m], id_variant);
			}  break;
			case FLUX_ESHE_VALUE: {
///////////////////////////
//...calculation heat flux;
				P[3] = -1.; P[4] = P[5] = 0.;
				this->B[i].shape->norm_local(P+3);
				this->B[i].shape->parametrization_grad(P);
				jump4(P, i, 0);

				F[0] = this->B[i].shape->potential(this->solver.hh[i][0][m], id_variant);

				P[4] = -1.; P[3] = P[5] = 0.;
				this->B[i].shape->norm_local(P+3);
				jump4(P, i, 0);

				F[1] = this->B[i].shape->potential(this->solver.hh[i][0][m], id_variant);

				P[5] = -1.; P[3] = P[4] = 0.;
				this->B[i].shape->norm_local(P+3);
				jump4(P, i, 0);

				F[2] = this->B[i].shape->potential(this->solver.hh[i][0][m], id_variant);
				if (this->B[i].link[this->NUM_PHASE] == -1)
					F[2] += P[5];
			}  break;
			case FLUX_COMPOS_VALUE: {
///////////////////////////
//...calculation heat flux;
				P[3] = -1.; P[4] = P[5] = 0.;
				this->B[i].shape->norm_local(P+3);
				this->B[i].shape->parametrization_grad(P);

				if (this->B[i].shape->type() == MP3D_CLAYER_SHAPE)	
					jump4_compos(P, i, 0); else jump4(P, i, 0);

				F[0] = this->B[i].shape->potential(this->solver.hh[i][0][m], id_variant);

				P[4] = -1.; P[3] = P[5] = 0.;
				this->B[i].shape->norm_local(P+3);

				if (this->B[i].shape->type() == MP3D_CLAYER_SHAPE)	
					jump4_compos(P, i, 0); else jump4(P, i, 0);

				F[1] = this->B[i].shape->potential(this->solver.hh[i][0][m], id_variant);

				P[5] = -1.; P[3] = P[4] = 0.;
				this->B[i].shape->norm_local(P+3);

				if (this->B[i].shape->type() == MP3D_CLAYER_SHAPE)	
					jump4_compos(P, i, 0); else jump4(P, i, 0);

				F[2] = this->B[i].shape->potential(this->solver.hh[i][0][m], id_variant);
			}  break;
        default : F[0] = T(i); F[1] = F[2] = 0.;
     }
  }
}

/////////////////////////////////////////////////
//...трехфазная модель для сферических включений;
template <typename T>
double CHeat3D<T>::TakeEshelby_two(double ff)
{
	double K1 = this->get_param(NUM_BASIC+NUM_SHIFT), 
			 K3 = this->get_param(NUM_BASIC), 
			 KH = K3*(1.+ff/((1.-ff)/3.+K3/(K1-K3)));
	return(KH);
}

/////////////////////////////////////////////////////////////
//...четырехфазная модель для сферических включений со слоем;
template <typename T>
double CHeat3D<T>::TakeEshelby(double ff, double ff_l)
{
	double K1 = this->get_param(NUM_BASIC+NUM_SHIFT), 
			 K2 = this->get_param(NUM_BASIC+NUM_SHIFT*2),
			 K3 = this->get_param(NUM_BASIC), c0 = ff+ff_l, c1 = ff/c0, KH;
	KH = K3*(((1.-c0)*(2.+c1)+K2/K3*(1.+2.*c0)*(1.-c1))+K1/K2*((1.-c0)*(1.-c1)+K2/K3*(1.+2.*c0)*(.5+c1)))/
			  (((1.+c0*.5)*(2.+c1)+K2/K3*(1.-c0)*(1.-c1))+K1/K2*((1.+c0*.5)*(1.-c1)+K2/K3*(1.-c0)*(.5+c1)));
	return(KH);
}

////////////////////////////////////////////////////////////////////
//...трехфазная модель для сферического включения (прямой алгоритм);
template <typename T>
double CHeat3D<T>::TakeEshelby_grad(double ff)
{
	double K_I = this->get_param(NUM_BASIC+NUM_SHIFT), K_M = this->get_param(NUM_BASIC),
			 C_I = this->get_param(this->NUM_GEOMT-1), C_M = C_I;
	if (C_I == 0. || C_M == 0.) { //...классическая трехфазная модель;
		return(K_M*(1.+ff/((1.-ff)/3.+K_M/(K_I-K_M))));
	}
	double RR1 = pow(ff, 1./3.), fR3 = 1./ff, 
			 kk_I = sqrt(C_I), tt_I = exp(-2.*kk_I), HH1 = ((1.+tt_I)*kk_I-(1.-tt_I))*fR3*.5, HH2 = ((1.-tt_I)*(sqr(kk_I)+2.)-2.*(1.+tt_I)*kk_I)*fR3*.5,
			 kk_M = sqrt(C_M), tt_D = exp(kk_M*(1./RR1-1.)),
			 JJP1 = (kk_M-1.)*fR3, JJP2 = (sqr(kk_M)+2.*(1.-kk_M))*fR3, JJM1 = -(kk_M+1.)*fR3, JJM2 = (sqr(kk_M)+2.*(1.+kk_M))*fR3, 
			 matr[7][8] = {
					{ 1., HH1, -1.,   -fR3, -JJP1, -JJM1, 0., 0.},
					{ 1., HH2, -1., 2.*fR3, -JJP2, -JJM2, 0., 0.},
					{ 0., K_I*HH1, 0.,0., -K_M*JJP1, -K_M*JJM1, 0., 0.},
					{ K_I, 0., -K_M, 2.*K_M*fR3, 0., 0., 0., 0.},
					{ 0., 0., 0., 0., (kk_M/RR1-1.)*tt_D, -(kk_M/RR1+1.)/tt_D, 0., 0.},
					{ 0., 0., K_M, -2.*K_M, 0., 0., -1., 0.},
					{ 0., 0., 1., 1., 0., 0., 0., 1.}, //...равенство температур дает единицу в правую часть;
	};

//////////////////////////////////////////////////////////////////
//...решаем систему линейных уравнений A0, C0, A1, B1, C1, D1, KH;
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

/////////////////////////////////////////////////////////
//...распределение температурного поля в среднем сечении;
	int id_visual = 0;
	if (id_visual) {
		FILE * TST = fopen("CCohes3D_sphere_eshelby.dat", "w");
		double A = .7, X, Y = 0., Z, RR, func; int M = 400;
		for (i = 0; i <= M; i++)
		for (l = 0; l <= M; l++) {
			X = (-1.+i*2./M)*A; Z = (-1.+l*2./M)*A; RR = sqrt(X*X+Y*Y+Z*Z); func = Z;
			if (RR < RR1) func *= (matr[6][0]+ matr[6][1]*(kk_I*RR/RR1*cosh(kk_I*RR/RR1)-sinh(kk_I*RR/RR1))/(RR*RR*RR)); else
			if (RR < 1.0) func *= (matr[6][2]+(matr[6][3]+matr[6][4]*(kk_M*RR/RR1-1.)*exp(kk_M*RR/RR1)-matr[6][5]*(kk_M*RR/RR1+1.)*exp(-kk_M*RR/RR1))/(RR*RR*RR)); 
			fprintf(TST, "%g   %g   %g\n", Z, X, func-Z);
		}
		fclose(TST);
		if (NULL_STATE) { //...проверка сшивки в одной точке;
			X = 0.; Z = RR1; RR = sqrt(X*X+Y*Y+Z*Z); double func1 = Z, func2 = Z;
			func1 *= (matr[6][0]+ matr[6][1]*(kk_I*RR/RR1*cosh(kk_I*RR/RR1)-sinh(kk_I*RR/RR1))/(RR*RR*RR));
			func2 *= (matr[6][2]+(matr[6][3]+matr[6][4]*(kk_M*RR/RR1-1.)*exp(kk_M*RR/RR1)-matr[6][5]*(kk_M*RR/RR1+1.)*exp(-kk_M*RR/RR1))/(RR*RR*RR)); 

			X = 0.; Z = 1.; RR = sqrt(X*X+Y*Y+Z*Z); func1 = Z; func2 = Z;
			func2 *= (matr[6][2]+(matr[6][3]+matr[6][4]*(kk_M*RR/RR1-1.)*exp(kk_M*RR/RR1)-matr[6][5]*(kk_M*RR/RR1+1.)*exp(-kk_M*RR/RR1))/(RR*RR*RR)); 
		}
	}
	return(matr[6][7]);
}

//////////////////////////////////////////////////////
//...двухфазная модель для эллипсоидального включения;
template <typename T>
double CHeat3D<T>::TakeEllipsoidEshelby(double ff, double eps, double phi_stream, int NX, int NY)
{
	double KI = this->get_param(NUM_BASIC+NUM_SHIFT), 
			 KM = this->get_param(NUM_BASIC), KK = KM/KI, KH = 0., z0, ch_alpha, sh_alpha, h1, h2, B1, B2, A1, A2;
	if (0. < eps && eps < 1.) {
		ch_alpha = 1./sqrt(1.-eps*eps); sh_alpha = eps*ch_alpha;
		h1 =  1./ch_alpha+(h2 = 0.5*log((ch_alpha-1.)/(ch_alpha+1.))); h2 += 1./(eps*sh_alpha);
		B1 = (1.-KK)/(KK*ch_alpha*sqr(eps-1./eps)-(1.-KK)*h1); A1 = 1.+h1*B1;
		B2 = (KK-1.)/(KK*ch_alpha*sqr(eps-1./eps)*2.+(1.-KK)*h2); A2 = 1.+h2*B2;
		z0 = ch_alpha;

//////////////////////////////////////////
//...вычисление температуры и теплопотока;
 		int id_visual = 1;
		if (id_visual) {
			CGrid * nd = CreateNodes();
			double cos_phi = 1., sin_phi = 0., R0 = 1., f0 = R0/z0, 
				cos_stream = cos(phi_stream), sin_stream = sin(phi_stream), A = 2.5, B = 2.5;

///////////////////////////////////////////////////////////
//...разбиваем плоскость сечения в пространстве параметров;
			nd->grid_box2D(NX, NY); 
			for (int i = 0; i < nd->N; i++) {
				nd->Z[i] = (nd->X[i]-.5)*A;
				nd->X[i] = (nd->Y[i]-.5)*B;
				nd->Y[i] = 0.; 
			}

////////////////////////////////////////////////////////////////////
//...распределение поля температуры и теплопотока в среднем сечении;
			FILE * TST = fopen("CHeat3D_spheroid.dat", "w");
			double func, stream1, stream2, X, Y, Z, f, ff, f1, f2, FN2;
			for (int i = 0; i < nd->N; i++) {
				Z = nd->Z[i]*cos_stream+nd->X[i]*sin_stream;
				X = nd->X[i]*cos_stream-nd->Z[i]*sin_stream;
				Y = nd->Y[i]; 
				func  = Z*cos_stream-X*sin_stream;
				stream1 = KM*cos_stream;
				stream2 = -KM*sin_stream;
				f = sqr(X)+sqr(f0);	f1 = sqr(Z);
				if ((ff = sqr(Z)+sqr(X/eps)) < sqr(R0)) {
					func = A1*Z*cos_stream-A2*X*sin_stream;
					stream1 =  KI*A1*cos_stream;            
					stream2 = -KI*A2*sin_stream;
				}
				else {
					f2 = (f1+2.*f-4.*sqr(f0))/sqr(f);
					if (f1 < 1e-5) 
						ch_alpha = (1.-f2*f*.5*(1.-f1*f2*.25))*.5;
					else
						ch_alpha = (1.-f/f1*(sqrt(1.+f1*f2)-1.))*.5;
					ch_alpha = 1./sqrt(ch_alpha); 
					sh_alpha = sqrt(sqr(ch_alpha)-1.);
					FN2 = sqr(f0*ch_alpha)-sqr(Z/ch_alpha);
					h1 = 1./ch_alpha+(h2 = 0.5*log((ch_alpha-1.)/(ch_alpha+1.))); h2 += ch_alpha/(sqr(ch_alpha)-1.);
					func += B1*h1*Z*cos_stream-B2*h2*X*sin_stream;
					stream1 = KM*( B1*h1*cos_stream+Z/(f0*FN2)*(B1*Z/(ch_alpha*sqr(ch_alpha))*cos_stream+2.*B2*X/(ch_alpha*sqr(sh_alpha))*sin_stream));
					stream2 = KM*(-B2*h2*sin_stream+X/(f0*FN2)*(B1*Z/(ch_alpha*sqr(sh_alpha))*cos_stream+2.*B2*X*ch_alpha/sqr(sqr(sh_alpha))*sin_stream));
				}
				fprintf(TST, "%g   %g   %g   %g   %g   %g   %g\n", nd->X[i], nd->Y[i], nd->Z[i],	func, func-nd->Z[i], 
					stream1*cos_stream-stream2*sin_stream, stream1*sin_stream+stream2*cos_stream);
			}
			fclose(TST);
		}
	}
	return(KH);
}
#undef  Message
#endif
