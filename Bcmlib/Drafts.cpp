#include "stdafx.h"
#include "drafts.h"

int CBase::NUM_MPLS  = 0;
int CBase::NUM_QUAD  = 2;
int CBase::NUM_TIME  = 4;
int CBase::NUM_LOCAL = 3;
int CBase::NUM_ADHES = 4;
int CBase::NUM_VIBRO = 4;
int CBase::NUM_GEOMT = 4;
int CBase::NUM_PHASE = 0;
int CBase::MAX_PHASE = 3;
int CBase::BOX_LINK_PERIOD = 0;

//////////////////////////////
//...construction real drafts;
CDraft<double> * CreateDraftR(Num_Draft id_DRAFT, int id_dop)
{
	return CreateDraft<double>(id_DRAFT, id_dop);
}

/////////////////////////////////
//...construction complex drafts;
CDraft<complex> * CreateDraftC(Num_Draft id_DRAFT, int id_dop)
{
	return CreateDraft<complex>(id_DRAFT, id_dop);
}

///////////////////////////////////////
//...construction double-double drafts;
CDraft<dd_real> * CreateDraftD(Num_Draft id_DRAFT, int id_dop)
{
	return CreateDraft<dd_real>(id_DRAFT, id_dop);
}

/////////////////////////////////////
//...construction quad-double drafts;
CDraft<qd_real> * CreateDraftQ(Num_Draft id_DRAFT, int id_dop)
{
	return CreateDraft<qd_real>(id_DRAFT, id_dop);
}

////////////////////////////
//...drafts identification;
Num_Draft draft_method(CBase * draft) 
{
	if (! draft) return NULL_DRAFT;
	if (typeid(* draft) == typeid(CBase)) return  BASIC_DRAFT; else
	if (typeid(* draft) == typeid(CDraft<double>)  ||
		 typeid(* draft) == typeid(CDraft<complex>) || 
		 typeid(* draft) == typeid(CDraft<dd_real>) ||
		 typeid(* draft) == typeid(CDraft<qd_real>)) return STRUCT_DRAFT; else
	if (typeid(* draft) == typeid(CLame2D)) return LAME2D_DRAFT; else
	if (typeid(* draft) == typeid(CLame3D)) return LAME3D_DRAFT; else
	if (typeid(* draft) == typeid(CHeat2D<double>)  ||
		 typeid(* draft) == typeid(CHeat2D<complex>) ||
		 typeid(* draft) == typeid(CHeat2D<dd_real>) ||
		 typeid(* draft) == typeid(CHeat2D<qd_real>)) return HEAT2D_DRAFT; else
	if (typeid(* draft) == typeid(CHeat3D<double>)  ||
		 typeid(* draft) == typeid(CHeat3D<complex>) ||
		 typeid(* draft) == typeid(CHeat3D<dd_real>) ||
		 typeid(* draft) == typeid(CHeat3D<qd_real>)) return HEAT3D_DRAFT; else
	if (typeid(* draft) == typeid(CHydro3D<double>)  ||
		 typeid(* draft) == typeid(CHydro3D<complex>) ||
		 typeid(* draft) == typeid(CHydro3D<dd_real>) ||
		 typeid(* draft) == typeid(CHydro3D<qd_real>)) return HYDRO3D_DRAFT; else
	if (typeid(* draft) == typeid(CCohes2D))		 return COHES2D_DRAFT; else
	if (typeid(* draft) == typeid(CCohes3D))		 return COHES3D_DRAFT; else
	if (typeid(* draft) == typeid(CMindl2D))		 return MINDL2D_DRAFT; else
	if (typeid(* draft) == typeid(CMindl3D))		 return MINDL3D_DRAFT; else
	if (typeid(* draft) == typeid(CVisco2D))		 return VISCO2D_DRAFT; else
	if (typeid(* draft) == typeid(CVisco2D_grad)) return VISCO2D_GRAD_DRAFT; else 
	if (typeid(* draft) == typeid(CPorosity2D))	 return	 POROSITY2D_DRAFT; else return NULL_DRAFT;
}

////////////////////////////////////////////////////////////////
//          INTERFACE FUNCTIONS FOR ESHELBY PROBLEMS          //
////////////////////////////////////////////////////////////////
double TakeLayer_EH(int N, double * ff, double * kk, double * ll)
{
	CCohes3D draft;
	CLame3D ldraft; double EH;
	if (ll == NULL) EH = ldraft.TakeLayer_kk(N, ff, kk);
	else				 EH =  draft.TakeLayer_kk(N, ff, kk, ll);
	return(EH);
}

double TakeCylinder_EH(int N, double * ff, double * kk)
{
	CLame2D ldraft; double EH;
	EH = ldraft.TakeLayer_k0(N, ff, kk);
	return(EH);
}

double TakeCylinder_KH1(int N, double * ff, double * kp, double * mu)
{
	CLame2D ldraft; double KH;
	KH = ldraft.TakeLayer_k1(N, ff, kp, mu);
	return(KH);
}

double TakeCylinder_GH1(int N, double * ff, double * kp, double * mu, double * nj)
{
	CLame2D ldraft; double GH;
	GH = ldraft.TakeLayer_G1(N, ff, kp, mu, nj);
	return(GH);
}

double TakeCylinder_KH2(int N, double * ff, double * lm, double * mu, double * nj)
{
	CLame2D ldraft; double LH;
	LH = ldraft.TakeLayer_k2(N, ff, lm, mu, nj);
	return(LH);
}

double TakeCylinder_GH2(int N, double * ff, double * mu, double * nj)
{
	CLame2D ldraft; double GH;
	GH = ldraft.TakeLayer_G2(N, ff, mu, nj);
	return(GH);
}

double TakeCylinder_GH2_simple(int N, double * ff, double * mu)
{
	CLame2D ldraft; double GH;
	GH = ldraft.TakeLayer_G2_simple(N, ff, mu);
	return(GH);
}

double TakeCylinder_SH(int N, double * ff, double * kp, double * mu, double * nj)
{
	CLame2D ldraft; double GH;
	GH = ldraft.TakeLayer_sh(N, ff, kp, mu, nj);
	return(GH);
}

double TakeSphere_KH(int N, double * ff, double * kv, double * mu)
{
	CLame3D ldraft; double KH;
	KH = ldraft.TakeLayer_k1(N, ff, kv, mu);
	return(KH);
}

double TakeSphere_GH(int N, double * ff, double * kv, double * mu, double * nj)
{
	CLame3D ldraft; double GH;
	GH = ldraft.TakeLayer_G1(N, ff, kv, mu, nj);
	return(GH);
}

double TakeGradMatrix_EH(double cI, double EI, double EM)
{
	double EH = EI*cI+EM*(1.-cI);
	return(EH);
}

double TakeGradMatrix_KH1(double cI, double EI, double EM, double nI, double nM, double lM)
{
	CCohes2D draft; double KH;
	draft.set_fasa_hmg(nM, nI, EM*.5/(1.+nM), EI*.5/(1.+nI), (lM?EM*(1.-nM)/(sqr(lM)*(1.+nM)*(1.-2.*nM)):0.), 0.);
	KH = draft.TakeGradMatrix_k1(cI);
	return(KH);
}

double TakeGradMatrix_GH1(double cI, double EI, double EM, double nI, double nM, double lM, double & det)
{
	CCohes2D draft; double GH;
	if (fabs(lM) != 0.) {
		draft.set_fasa_hmg(nM, nI, EM*.5/(1.+nM), EI*.5/(1.+nI), EM*(1.-nM)/(sqr(lM)*(1.+nM)*(1.-2.*nM)), 0.);
		GH = draft.TakeGradMatrix_G1(cI, det);
	}
	else {
		draft.set_fasa_hmg(nM, nI, EM*.5/(1.+nM), EI*.5/(1.+nI), 0., 0.);
		GH = draft.TakeGradMatrix_G1_classic(cI, det);
	}
	return(GH);
}

double TakeGradLayer_EH(double cI, double cL, double EI, double EL, double EM)
{
	double EH = EI*cI+EL*cL+EM*(1.-cI-cL);
	return(EH);
}

double TakeGradLayer_KH1(double cI, double cL, double EI, double EL, double EM, double nI, double nL, double nM, double lL)
{
	CCohes2D draft; double KH;
	draft.set_fasa_hmg(nM, nI, nL, EM*.5/(1.+nM), EI*.5/(1.+nI), EL*.5/(1.+nL), 0., 0., (lL?EL*(1.-nL)/(sqr(lL)*(1.+nL)*(1.-2.*nL)):0.));
	KH = draft.TakeGradLayer_k1(cI, cL);
	return(KH);
}

double TakeGradLayer_GH1(double cI, double cL, double EI, double EL, double EM, double nI, double nL, double nM, double lL)
{
	CLame2D ldraft;
	CCohes2D draft; double GH;
	if (fabs(lL) != 0.) {
		draft.set_fasa_hmg(nM, nI, nL, EM*.5/(1.+nM), EI*.5/(1.+nI), EL*.5/(1.+nL), 0., 0., EL*(1.-nL)/(sqr(lL)*(1.+nL)*(1.-2.*nL)));
		GH = draft.TakeGradLayer_G1(cI, cL);
	}
	else {
		const int N = 3;
		double ff[N] = {cI, cL, 1.-cI-cL}, mu[N] ={EI*.5/(1.+nI), EL*.5/(1.+nL), EM*.5/(1.+nM)},
				 kp[N] = {mu[0]/(1.-2.*nI), mu[1]/(1.-2.*nL), mu[2]/(1.-2.*nM)}, nj[N] ={nI, nL, nM};
		GH = ldraft.TakeLayer_G1(N, ff, kp, mu, nj);
	}
	return(GH);
}

double TakeGradLayer_KH2(double cI, double cL, double EI, double EL, double EM, double nI, double nL, double nM, double lL)
{
	CLame2D ldraft;
	CCohes2D draft; double LH;
	if (fabs(lL) != 0.) {
		draft.set_fasa_hmg(nM, nI, nL, EM*.5/(1.+nM), EI*.5/(1.+nI), EL*.5/(1.+nL), 0., 0., EL*(1.-nL)/(sqr(lL)*(1.+nL)*(1.-2.*nL)));
		LH = draft.TakeGradLayer_k2(cI, cL);
	}
	else {
		const int N = 3;
		double ff[N] = {cI, cL, 1.-cI-cL}, mu[N] ={EI*.5/(1.+nI), EL*.5/(1.+nL), EM*.5/(1.+nM)},
				 lm[N] = {mu[0]*2.*nI/(1.-2.*nI), mu[1]*2.*nL/(1.-2.*nL), mu[2]*2.*nM/(1.-2.*nM)}, nj[N] ={nI, nL, nM};
		LH = ldraft.TakeLayer_k2(N, ff, lm, mu, nj);
	}
	return(LH);
}

double TakeGradLayer_GH2(double cI, double cL, double EI, double EL, double EM, double nI, double nL, double nM, double lL)
{
	CLame2D ldraft;
	CCohes2D draft; double GH;
	if (fabs(lL) != 0.) { 
		draft.set_fasa_hmg(nM, nI, nL, EM*.5/(1.+nM), EI*.5/(1.+nI), EL*.5/(1.+nL), 0., 0., EL*(1.-nL)/(sqr(lL)*(1.+nL)*(1.-2.*nL)));
		GH = draft.TakeGradLayer_G2(cI, cL);
	}
	else {
		const int N = 3;
		double ff[N] = {cI, cL, 1.-cI-cL},
				 mu[N] = {EI*.5/(1.+nI), EL*.5/(1.+nL), EM*.5/(1.+nM)}, nj[N] = {nI, nL, nM};
		GH = ldraft.TakeLayer_G2(N, ff, mu, nj);
	}
	return(GH);
}

double TakeGradLayer_GH2_simple(double cI, double cL, double EI, double EL, double EM, double nI, double nL, double nM, double lL)
{
	CLame2D ldraft;
	CCohes2D draft; double GH;
	draft.set_fasa_hmg(nM, nI, nL, EM*.5/(1.+nM), EI*.5/(1.+nI), EL*.5/(1.+nL), 0., 0., (lL?EL*(1.-nL)/(sqr(lL)*(1.+nL)*(1.-2.*nL)):0.));
	GH = draft.TakeGradLayer_G2_simple(cI, cL);
	return(GH);
}

double TakeSphere_GH_det(double c0, double nju1, double nju2, double E1, double E2, double alpha)
{
	CLame3D ldraft; double GH;
	ldraft.set_fasa_hmg(nju2, nju1, E2/(1.+nju2)*.5, E1/(1.+nju1)*.5);
	GH = ldraft.TakeEshelby_shear_det(c0, alpha);
	return(GH);
}

double TakeSphere_volm_two(double c0, double nju1, double nju2, double E1, double E2, double l1, double l2, double AA)
{
	CCohes3D draft; 
	CLame3D ldraft; double KH;
	if (l1 == 0. || l2 == 0.) { 
		ldraft.set_fasa_hmg(nju1, nju2, E1/(1.+nju1)*.5, E2/(1.+nju2)*.5);
		KH = ldraft.TakeEshelby_volm_two(c0);
	}
	else {
		draft.set_fasa_hmg(nju1, nju2, E1/(1.+nju1)*.5, E2/(1.+nju2)*.5, E1/((1.+nju1)*sqr(l1))*.5, E2/((1.+nju2)*sqr(l2))*.5);
		draft.set_adhesion(AA, 0.);
		KH = draft.TakeEshelby_volm_two(c0);
	}
	return(KH);
}

double TakeSphere_shear_two(double c0, double nju1, double nju2, double E1, double E2, double l1, double l2, double AA, double BB)
{
	CCohes3D draft; 
	CLame3D ldraft; double GH;
	if (l1 == 0. || l2 == 0.) { 
		ldraft.set_fasa_hmg(nju1, nju2, E1/(1.+nju1)*.5, E2/(1.+nju2)*.5);
		GH = ldraft.TakeEshelby_shear_two(c0);
	}
	else {
		draft.set_fasa_hmg(nju1, nju2, E1/(1.+nju1)*.5, E2/(1.+nju2)*.5, E1/((1.+nju1)*sqr(l1))*.5, E2/((1.+nju2)*sqr(l2))*.5);
		draft.set_adhesion(AA, BB);
		GH = draft.TakeEshelby_shear_two(c0);
	}
	return(GH);
}

double TakeSphere_volm_sym(double c0, double nju1, double nju2, double E1, double E2, double l1, double l2)
{
	CCohes3D draft; 
	CLame3D ldraft; double KH;
	if (l1 == 0. || l2 == 0.) { 
		ldraft.set_fasa_hmg(nju1, nju2, E1/(1.+nju1)*.5, E2/(1.+nju2)*.5);
		KH = ldraft.TakeEshelby_volm_two(c0);
	}
	else {
		draft.set_fasa_hmg(nju1, nju2, E1/(1.+nju1)*.5, E2/(1.+nju2)*.5, E1/((1.+nju1)*sqr(l1))*.5, E2/((1.+nju2)*sqr(l2))*.5);
		KH = draft.TakeEshelby_volm_sym(c0);
	}
	return(KH);
}

double TakeSphere_shear_sym(double c0, double nju1, double nju2, double E1, double E2, double l1, double l2)
{
	CCohes3D draft; 
	CLame3D ldraft; double GH;
	if (c0 < 1e-3) return  (GH = E1/(1.+nju1)*.5);
	if (fabs(l1) <= 0.01 && fabs(l2) <= 0.01) { 
		GH = ldraft.TakeEshelby_shear(c0, nju2, nju1, E2, E1);
	}
	else {
		if (l2 < 0.) l2 = 0.03;
		GH = draft.TakeEshelby_shear_old(c0, nju2, nju1, E2, E1, l2, l1);
	}
	return(GH);
}

double TakeSphere_volm_grad(double c0, double nju1, double nju2, double E1, double E2, double ll)
{
	CCohes3D draft; 
	CLame3D ldraft; double KH;
	if (ll == 0. ) { 
		ldraft.set_fasa_hmg(nju1, nju2, E1/(1.+nju1)*.5, E2/(1.+nju2)*.5);
		KH = ldraft.TakeEshelby_volm_two(c0);
	}
	else {
		draft.set_fasa_hmg(nju1, nju2, E1/(1.+nju1)*.5, E2/(1.+nju2)*.5, E1/((1.+nju1)*sqr(ll))*.5, 0.);
		KH = draft.TakeEshelby_volm_grad(c0);
	}
	return(KH);
}

double TakeSphere_volm_rigd(double c0, double nju1, double E1, double ll)
{
	CCohes3D draft; double KH;
	draft.set_fasa_hmg(nju1, 0., E1/(1.+nju1)*.5, 0., (! ll ? 0.: E1/((1.+nju1)*sqr(ll))*.5), 0.);
	KH = draft.TakeEshelby_volm_rigd(c0);
	return(KH);
}

double TakeSphere_shear_grd(double c0, double nju1, double nju2, double E1, double E2, double ll)
{
	CCohes3D draft; 
	CLame3D ldraft; double GH;
	if (ll == 0.) { 
		ldraft.set_fasa_hmg(nju1, nju2, E1/(1.+nju1)*.5, E2/(1.+nju2)*.5);
		GH = ldraft.TakeEshelby_shear_two(c0);
	}
	else {
		draft.set_fasa_hmg(nju1, nju2, E1/(1.+nju1)*.5, E2/(1.+nju2)*.5, E1/((1.+nju1)*sqr(ll))*.5, 0.);
		GH = draft.TakeEshelby_shear_grd(c0);
	}
	return(GH);
}

double TakeSphere_shear_rig(double c0, double nju1, double E1, double ll)
{
	CCohes3D draft; 
	CLame3D ldraft; double GH;
	if (ll == 0.) { 
		ldraft.set_fasa_hmg(nju1, 0., E1/(1.+nju1)*.5, 0.);
		GH = ldraft.TakeEshelby_shear_rig(c0);
	}
	else {
		draft.set_fasa_hmg(nju1, 0., E1/(1.+nju1)*.5, 0., E1/((1.+nju1)*sqr(ll))*.5, 0.);
		GH = draft.TakeEshelby_shear_rig(c0);
	}
	return(GH);
}

double TakeSphere_shear(double c0, double nju1, double nju2, double E1, double E2, double l1, double l2)
{
	CCohes3D draft; 
	CLame3D ldraft; double GH;
	if (c0 < 2e-3) return  (GH = E1/(1.+nju1)*.5);
	if (fabs(l1) < 0.01 && fabs(l2) < 0.01) { 
		GH = ldraft.TakeEshelby_shear(c0, nju2, nju1, E2, E1);
	}
	else {
		if (l2 < 0.) l2 = 0.03;
		GH = draft.TakeEshelby_shear(c0, nju2, nju1, E2, E1, l2, l1);
	}
	return(GH);
}

/////////////////////////////////////////////////////////////////////////////////////
//...realization gradient model with option of displacement and stresses calculation;
void TakeEshelbyModel(double cI, double cL, double EI, double EL, double EM, double nI, double nL, double nM, double lI, double lL, double lM,
	double & EH, double & GH, double & nH, double rad, double fi, double displ[2], double sigma[3], double fK, double fG)
{
	CDraft<double> * sm = (CDraft<double> *)CreateDraftR(COHES2D_DRAFT, 8);
	sm->set_fasa_hmg(nM, nI, nL, EM*.5/(1.+nM), EI*.5/(1.+nI), EL*.5/(1.+nL), (lM ? EM*.5/((1.+nM)*sqr(lM)) : 0.), (lI ? EI*.5/((1.+nI)*sqr(lI)) : 0.), (lL ? EL*(1.-nL)/((1.+nL)*sqr(lL)*(1.-2.*nL)) : 0.));
	if (! lL || ! lM) sm->TakeEshelbyGradModel(cI, cL, fK, fG); else sm->TakeEshelbyModel(cI, cL, fK, fG);

	double rad0 = 1., rad1 = 1./sqrt(cI/(cI+cL)), rad2 = 1./sqrt(cI), F[3], X, Y;
  	EH = sm->get_Young(0);
	GH = sm->get_Shear(0);
	nH = sm->get_Poisn(0);

	int hit = 1;
	if (rad < rad0) hit = 0; else  
	if (rad < rad1) hit = 2; else 
	if (rad < rad2) hit = 3;	

	sm->GetFuncAllValues(X = rad*cos(fi), Y = rad*sin(fi), 0., F, hit, DISPL_GRAD_VALUE); displ[0] = F[0]; displ[1]= F[1];
	sm->GetFuncAllValues(X, Y, 0., F, hit,  STRESS_X_GRAD_VALUE); sigma[0] = F[0]; sigma[1]= F[1];
	sm->GetFuncAllValues(X, Y, 0., F, hit,  STRESS_Y_GRAD_VALUE); sigma[2] = F[1];
	
	delete sm;
}

/////////////////////////////////////////////////////////
//            INTERPHASE LAYER APPROXIMATION           //
/////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//...function of approximation with ll = 1./sqrt(C0+sqr(t)*(C1+t*(C2+...+t*CN)...), coef[N+1], data[2*N+2];
int strata_approx(int N, double * data, double * coef)
{
	int  * ii = new_struct<int>(N+1), i, k, l, k0, l0;
	double f, ** matrix = 0; set_matrix(matrix, N+1, N+1);

/////////////////////////////////////////////////
//...transferring given data of interphase layer;
	for (i = 0; i <= N; i++) coef[i] = 1./sqr(data[2*i+1]); 

/////////////////////////////////////////////////////////////////////////
//...filling system of linear algebraiq equations with Wandermode matrix;
	for (i = 0; i <= N; i++)
	for (f = 1., k = 0; k <= N; k++, f *= data[2*i], f *= (k == 1 ? data[2*i] : 1.)) matrix[i][k] = f; 

////////////////////////////////////////////////////////////////////////////////////////
//...defining coefficients by solving system of linear equation wwith Wandermode matrix;
	for (i = 0; i <= N; i++) {
		for (f = 0., k = 0; k <= N; k++) //...look for position maximal element;
			if (ii[k] != 1) 
				for (l = 0; l <= N; l++) 
					if (! ii[l]) {
						if (fabs(matrix[k][l]) >= f) f = fabs(matrix[k0 = k][l0 = l]); 
					}
					else if (ii[l] > 1) {
						delete_struct(ii); delete_struct(matrix);	return(-1);
					}
		++(ii[l0]);
///////////////////////////////////////////////////////////
//...swapping row for diagonal position of maximal element;
		if (k0 != l0)
			for (swap(coef[k0], coef[l0]), l = 0; l <= N; l++) swap(matrix[k0][l], matrix[l0][l]);
		if (matrix[l0][l0] == 0.) {
			delete_struct(ii); delete_struct(matrix);	return(0);
		}
////////////////////////////////
//...diagonal row normalization;
		f = 1./matrix[l0][l0]; matrix[l0][l0] = 1.;
		for (coef[l0] *= f, l = 0; l <= N; l++) matrix[l0][l] *= f;
////////////////////////////////
//...elimination all other rows;
		for (k = 0; k <= N; k++)
			if ( k != l0) {
				f = matrix[k][l0]; matrix[k][l0] = 0.;
				for (coef[k] -= coef[l0]*f, l = 0; l <= N; l++) matrix[k][l] -= matrix[l0][l]*f;
			}
	}
	delete_struct(ii); delete_struct(matrix);
	return(1);
}

////////////////////////////////////////////////////////////////
//...approximation of aglomeration radius by Lagrange expansion;
double strata_aglom(double t, int N, double * coef)
{
	double sum = 0;
	for (int i = N; i > 0; i--) sum = coef[i]+t*sum; 
	return 1./sqrt(sum = coef[0]+sqr(t)*sum);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//...approximation of intarphase layer by least square method ll = 1./sqrt(ll0**(-2)+B*(t**2-t0**2)), data[2*N+2];
double strata_approx(int N, double * data, double ll0)
{
	double sum = 0., fnorm = 0.;

/////////////////////////////////////
//...defining scalar produc and norm;
	for (int i = N; i <= N; i++) {
		sum += (sqr(data[2*i])-sqr(data[0]))*(1./sqr(data[2*i+1])-1./sqr(ll0));
		fnorm += sqr(sqr(data[2*i])-sqr(data[0]));
	}
	return(sum = sum/fnorm);
}

/////////////////////////////////////////////////////////////////
//...approximation of aglomeration radius by least square method;
double strata_aglom(double t, double t0, double B, double ll0)
{
	return 1./sqrt(1./sqr(ll0)+B*(sqr(t)-sqr(t0)));
}

/////////////////////////////////////////////////////////
//...approximatiob of aglomeration radius by exponential;
double aglom_approx(double t, double alpha, double r_max, double r_min, double t_min)
{
	return max(r_min, r_max-(r_max-r_min)*exp(-alpha*(t-t_min)));
}

/////////////////////////////////////////////////////////
//...approximatiob of aglomeration radius by exponential;
double rigid_approx(double t, double E1_max, double E2_max, double E_min, double t_min, double E1, double E2, double E3, double t1, double t2, double t3)
{
	double a1 = -log((E1_max-E1)/(E1_max-E_min)), b1 = -log((E1_max-E2)/(E1_max-E_min)), c1 = -log((E2_max-E3)/(E2_max-E_min)),
			 alpha = (b1*sqr(t1)*sqr(t1)-a1*sqr(t2)*sqr(t2))/(sqr(t1)*sqr(t2)*(sqr(t1)-sqr(t2))),
			 alpha2 = -(b1*sqr(t1)-a1*sqr(t2))/(sqr(t1)*sqr(t2)*(sqr(t1)-sqr(t2))), beta = c1/sqr(t3);

		FILE *  TST = fopen("WSMP_grad2.dat", "w");
		fprintf(TST, " %g, %g, %g\n", alpha/sqr(t_min*100.), alpha2/sqr(sqr(t_min*100.)), t_min);
		fclose (TST);

	return (t < t_min ? E1_max-(E1_max-E_min)*exp(-sqr(t/t_min-1.)*(alpha+alpha2*sqr(t/t_min-1.))) : E2_max-(E2_max-E_min)*exp(-beta*sqr(t/t_min-1.)));
}
double rigid_approx(double t, double E1_max, double E2_max, double E_min, double t_min, double E1, double E2, double E3, double E4, double t1, double t2, double t3, double t4)
{
	double a1 = -log((E1_max-E1)/(E1_max-E_min)), b1 = -log((E1_max-E2)/(E1_max-E_min)), 
			 c1 = -log((E2_max-E3)/(E2_max-E_min)), d1 = -log((E2_max-E4)/(E2_max-E_min)),
			 alpha = (b1*sqr(t1)*sqr(t1)-a1*sqr(t2)*sqr(t2))/(sqr(t1)*sqr(t2)*(sqr(t1)-sqr(t2))),
			 alpha2 = -(b1*sqr(t1)-a1*sqr(t2))/(sqr(t1)*sqr(t2)*(sqr(t1)-sqr(t2))), 
			 beta = (d1*sqr(t3)*sqr(t3)-c1*sqr(t4)*sqr(t4))/(sqr(t3)*sqr(t4)*(sqr(t3)-sqr(t4))),
			 beta2 = -(d1*sqr(t3)-c1*sqr(t4))/(sqr(t3)*sqr(t4)*(sqr(t3)-sqr(t4)));
	return (t < t_min ? E1_max-(E1_max-E_min)*exp(-sqr(t/t_min-1.)*(alpha+alpha2*sqr(t/t_min-1.))) : 
							  E2_max-(E2_max-E_min)*exp(-sqr(t/t_min-1.)*(beta+beta2*sqr(t/t_min-1.))));
}

/////////////////////////////////////////////////////////////////////////////
//...approximatiob of aglomeration radius by exponential and power functions;
double rigid_approx(double t, double E1_max, double E2_max, double E_min, double alpha, double beta, double t_min)
{
	return (t < t_min ? E1_max-(E1_max-E_min)*exp(-alpha*sqr(t/t_min-1.)) : E2_max-(E2_max-E_min)*exp(-beta*sqr(t/t_min-1.)));
	//return (t < t_min ? (61.138+27.877*(t/t_min-1.))*sqr(t/t_min-1.)+60. : E2_max-(E2_max-E_min)*exp(-beta*sqr(t/t_min-1.)));
}

/////////////////////////////////////////////////////////////////
//            APPROXIMATION OF ELASTOMERIC MATERIALS           //
/////////////////////////////////////////////////////////////////
void approx_all(double & E1, double & E2, double E3, double & delta, double muni[][2], int N_muni, Num_Approx approx = EXPONENT_APPROX)
{
	double A1 = 0., A2 = 0., A3 = 0., h1 = 0., h2 = 0., sigma = 0., F1, F2;
	for (int i = 0; i < N_muni; i++) {
		double lambda = muni[i][0];
		//double lambda = 1.+muni[i][0]/100.; //...исходное - деформация в %;
		double invariant = sqrt(sqr(lambda)+2./lambda-3.), jacobian = lambda-1./sqr(lambda), eps = lambda-1.;
		if (eps < 0.0001) jacobian = sqrt(3.)*(1.-eps*.5*(1.-eps*(1.-(31./27.)*eps)));
		else jacobian /= invariant;  
/////////////////////////////////////////
//...acuumulates norm and scala products;
		sigma += sqr(muni[i][1]);
		F1 = invariant;
		F2 = sqr(F1);
		switch (approx) {
			case LOGARITHM1_APPROX: F2 = log(abs(1.+E3*invariant))/E3; F1 -= F2; break;
			case LOGARITHM2_APPROX: F2 = 1./(1.+E3*invariant)-1.; break;
			case FRACTIONAL_APPROX: F2 = E3*.5*(1./sqr(1.+E3*invariant)-1.); break;
			case EXPONENT_APPROX: F2 = (1.-exp(-E3*invariant))/E3; break;
			case POWER_APPROX: F2 = (pow(1.+invariant, 1.-E3)-1.)/(1.-E3); break;
		}
		F1 *= jacobian;
		F2 *= jacobian;
		A1 += F1*F1;
		A2 += F2*F2;
		A3 += F1*F2;
		h1 += F1*muni[i][1];
		h2 += F2*muni[i][1];
	}
	E1 = (A3*h2-A2*h1)/(A3*A3-A1*A2);
	E2 = (A3*h1-A1*h2)/(A3*A3-A1*A2);
	delta = sqrt(1.-(h1*E1+h2*E2)/sigma);
}

////////////////////////////////////////////////////
//...realization of deflection of gradient mambrane;
double w0_membrane(double R, double eps, int k_limit) 
{
	int k = 0;
	double h = 0., dd = sqr(.5*R), f = 1., I0 = 1., K0 = h, w0 = 0.;
	do {
		k +=1;
		f *= dd/sqr(k);
		h += 1./k;
		I0 += f;
		K0 += f*h;
	}
	while (f >= eps && f*h >= eps && k < k_limit);
	w0 = K0/I0;
	return(w0);
}

///////////////////////////////////////////////////////////////////////
//            FUNDAMENTAL SOLUTIONS IN GRADIENT ELASTICITY           //
///////////////////////////////////////////////////////////////////////
/////////////////////////////
//...generalized Grin tensor;
double Grin_tensor(int i, int k, double X, double Y, double Z, double mu, double nju, double s)
{
	double r = sqrt(sqr(X)+sqr(Y)+sqr(Z)), h0 = s ? exp(-r/s) : 0., d = 16.*M_PI*mu*(1.-nju),	P[3] = {X, Y, Z}, Grin = 0.;
	if (s) {
		if (r < 0.0001) {
			if (i == k)	Grin += (10./3.-4.*nju)/s-(7./4.-2.*nju)*r/sqr(s);
			Grin += P[i-1]*P[k-1]/(4.*sqr(s));
		}
		else {
			if (i == k) {
				Grin += (4.*(1.-nju)*(1.-h0)-1.)/r;
				Grin -= (h0*(1.+r/s)-1.)*2.*sqr(s)/(r*sqr(r));
			}
			Grin += P[i-1]*P[k-1]*(sqr(r)-2.*sqr(s)*(3.-h0*(sqr(r/s)+3.*r/s+3.)))/(r*sqr(r)*sqr(r));
		}
	}
	else {
		if (i == k) Grin += (3.-4.*nju)/r;
		Grin += P[i-1]*P[k-1]/(r*sqr(r));
	}
	return(Grin /= d);
}

/////////////////////////////
//...gradient Flaman problem;
void Flaman_displ(double X, double Y, double & RX, double & RY, double mu, double nju, double s)
{
	double r = sqrt(sqr(X)+sqr(Y)), h0 = r ? log(r) : 0., fi = arg2(Y, X), d = 2.*M_PI*mu, K0 = 0., K2 = 0.;

	if (s) {
		int k = 0, k_limit = 1000;
		double h1 = -0.5772156649, h2 = 1.5+h1, dd = sqr(.5*r/s), f1 = 1., f2 = dd/2., h3 = 0.;
		K0 = log(2.*s)+h1;
		K2 = 0.5+(h3 = h0-log(2.*s)-0.5*(h1+h2))*f2; h1 = K0-h0; 
		do {
			k +=1;
			f1 *= dd/sqr(k);
			f2 *= dd/(k*(k+2.));
			h1 += 1./k;
			h3 -= .5/k+.5/(k+2); 
			K0 += f1*h1;
			K2 += f2*h3;
		}
		while (f1*fabs(h1) >= EE && f2*fabs(h3) >= EE && k < k_limit);

		RX = 2.*(1.-nju)*K0-sqr(cos(fi))+(sqr(cos(fi))-sqr(sin(fi)))*K2;
		RY = (1.-2.*nju)*fi+cos(fi)*sin(fi)*(2.*K2-1.);
	}
	else {
		RX = 2.*(1.-nju)*h0-sqr(cos(fi));
		RY = (1.-2.*nju)*fi-cos(fi)*sin(fi);
	}
	RX /= d;
	RY /= d;
}

//////////////////////////////////
//...gradient Bussinesque problem;
void Boussinesq_displ(double X, double Y, double Z, double & RX, double & RY, double & RZ, double mu, double nju, double s = 0.)
{
	double r = sqrt(sqr(X)+sqr(Y)+sqr(Z)), h0 = r ? 1./r : 0., h1 = r+Z ? 1./(r+Z) : 0., d = 4.*M_PI*mu, 
			d3 = (1.-2.*nju)*h1, d4 = Z*sqr(h0), d5 = sqr(Z*h0), d1 = (d3-d4)*h0, d2 = (2.*(1.-nju)+d5)*h0, d0;
	int M_iter = 1000, k;
	if (s) {
		if (r < 10.*s || 0) {
			for (d1 = d2 = 0., d0 = 1., k = 0; k <= M_iter; k++, d0 *= -r/(s*k)) {
				d1 -= (-(1.-2.*nju)*h1/(k+1.)+Z*h0/s*2./((k+2.)*(k+4)))*d0/s;
				d2 += (((3.-2.*nju)-2./(k+3.))/(k+1.)+sqr(Z)*h0/s*2./((k+2.)*(k+4)))*d0/s;
			}
		}
		else {
			d1 -= (2.*d4*(1.+3.*s*h0)+d3)*(h1 = h0*exp(-r/s)); d1 += 6.*sqr(s*h0)*d4*(h0-h1);
			d2 += (2.*d5*(1.+3.*s*h0)-(3.-2.*nju)-2.*s*h0)*h1; d2 -= 2.*sqr(s*h0)*(3.*d5-1.)*(h0-h1);
		}
	}
	RX = -(d1 /= d)*X;
	RY = -d1*Y;
	RZ = (d2 /= d);
}

void Boussinesq_sigma(double X, double Y, double Z, double & SX, double & SY, double & SZ, double s = 0.)
{
	double r = sqrt(sqr(X)+sqr(Y)+sqr(Z)), h0 = r ? 1./r : 0., h1 = -h0, d = 2.*M_PI, 
			dX = X*sqr(h0), dY = Y*sqr(h0), dZ = Z*sqr(h0), d3 = sqr(Z*h0), d2 = -3.*d3*h0, d1 = d2, d0;
	int M_iter = 1000, k;
	if (s) {
		if (r < 10.*s || 0) {
			for (d1 = d2 = 0., d0 = r, k = 0; k <= M_iter; k++, d0 *= -r/(s*k)) {
				d1 -= ((2.*(d3-1.)/(k+4.)+1.)/(k+2.)+sqr(Z)*h0/s*2./((k+3.)*(k+5.)))*d0/sqr(s);
				d2 -= ((2.*(d3-3.)/(k+4.)+3.)/(k+2.)+sqr(Z)*h0/s*2./((k+3.)*(k+5.)))*d0/sqr(s);
			}
		}
		else {
			d1 += ((2.*d3-1.)*(1.+r/s)+(5.*d3-1.)*(3.+r/s)*2.*s*h0)*(h1 *= exp(-r/s)); d1 += 6.*sqr(s*h0)*(5.*d3-1.)*(h0+h1);
			d2 += ((2.*d3-3.)*(1.+r/s)+(5.*d3-3.)*(3.+r/s)*2.*s*h0)*h1; d2 += 6.*sqr(s*h0)*(5.*d3-3.)*(h0+h1);
		}
	}
	SX = (d1 /= d)*dX;
	SY =  d1*dY;
	SZ = (d2 /= d)*dZ;
}

//////////////////////////////////////////////
//...skew Bussinesque problem (not ready !!!);
void Boussiskew_displ(double X, double Y, double Z, double & RX, double & RY, double & RZ, double mu, double nju, double s)
{
	double r = sqrt(sqr(X)+sqr(Y)+sqr(Z)), h0 = r ? 1./r : 0., h1 = r+Z ? 1./(r+Z) : 0., d = 4.*M_PI*mu, 
			d3 = (1.-2.*nju)*h1, d4 = Z*sqr(h0), d5 = sqr(Z*h0), d1 = (d3-d4)*h0, d2 = (2.*(1.-nju)+d5)*h0;
	if (s) {
		d1 -= (2.*d4*(1.+s*h0)+d3)*(h1 = h0*exp(-r/s));
		d2 += (2.*d5*(1.+s*h0*3.)-(3.-2.*nju)-s*h0)*h1;
	}
	RX = (d1 /= d)*X;
	RY =  d1*Y;
	RZ = (d2 /= d);
}

void Boussiskew_sigma(double X, double Y, double Z, double & SX, double & SY, double & SZ, double s)
{
	double r = sqrt(sqr(X)+sqr(Y)+sqr(Z)), h0 = r ? 1./r : 0., h1 = -h0, d = 2.*M_PI, 
			dX = X*sqr(h0), dY = Y*sqr(h0), dZ = Z*sqr(h0), d3 = sqr(Z*h0), d2 = -3.*d3*h0, d1 = d2;
	if (s) {
		d1 += (2.*d3*(1.+s*h0*15.)-1.-s*h0*6.)*(h1 *= exp(-r/s)*(1.+r/s));
		d2 += (2.*d3*(1.+s*h0*15.)-3.-s*h0*18.)*h1+sqr(s*h0)*h0*24.;
	}
	SX = (d1 /= d)*dX;
	SY =  d1*dY;
	SZ = (d2 /= d)*dZ;
}
