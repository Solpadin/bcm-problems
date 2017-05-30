#include "stdafx.h"
#include "kernel.h"

///////////////////////////////////////////////////
//...removing all internal structure of space cell;
void release(Cells * ce)
{
	if (ce) {
		for (int l, m = ce->graph ? ce->graph[0] : -1, k = 0; k <= m; k++)
		if (ce->ce && ce->ce[k]) {
			delete_struct(ce->ce[k]->mp); l = ce->ce[k]->pm && ce->ce[k]->graph ? ce->ce[k]->graph[1] : 0;
			delete_struct(ce->ce[k]->graph);

			for (l--; l >= 0; l--)
				delete_struct(ce->ce[k]->pm[l]);
				delete_struct(ce->ce[k]->pm);

			if (k != m) delete_struct(ce->ce[k]);
		}
		delete_struct(ce->ce);
	}
	return;
}

///////////////////////////////////
//...deleting compose body surface;
void delete_bar(int N, Bar ** & bar)
{
	for (N--; bar && N >= 0; N--)
		delete_bar(bar[N]);
		delete_struct(bar);
	return;
}

///////////////////////////////////////////////////////////////////////////
//...reading geometrical card and topology in ASCII code (working program);
CMap * map_in(char * id_MAP, unsigned long & count, unsigned long upper_limit)
{
	CMap * mp = NULL;
	if (id_MAP && count < upper_limit) {
		char buf[BUF_SIZE+1];
		CMap  id_mp;
		int  m1, m2, m3, m4, m5, k, m, i;
		if ( user_Read (buf, id_MAP, count, upper_limit)) id_mp = user_strtod(buf); else goto err;
		m  = size_of_map(i = map_dim(& id_mp), k = map_genus(& id_mp))+1;

		if ((mp = get_map(i, k)) != NULL) for (i = k = 1; k < m; k++)
		if ( user_Read(buf, id_MAP, count, upper_limit)) mp[i++] = user_strtod(buf); else goto err;

////////////////////////////////////////////////
//...defining and writing additional parameters;
      if (! add_new_maps(mp, m, k = size_of_dop(mp))) goto err;
      switch ((int)mp[i-1]) {
			case      NULL_CELL: goto err;
			case  B1SPLINE_CELL: {
					m1 = (int)mp[7]/2;
					m2 = (int)mp[8]+2;
					m3 = (int)mp[9];

					for (k = 0; k < m2; k++)
					if  (user_Read(buf, id_MAP, count, upper_limit)) mp[i++] = user_strtod(buf); else goto err;
					for (k = 0; k < m3; k++)
					for (m = 0; m < m1; m++)
					if  (user_Read(buf, id_MAP, count, upper_limit)) mp[i++] = user_strtod(buf); else goto err;
			}		break;
			case  B2SPLINE_CELL: {
					m1 = (int)mp[7]/4;
					m2 = (int)mp[8]+2;
					m3 = (int)mp[9]+2;
					m4 = (int)mp[10];
					m5 = (int)mp[11];

					for (k = 0; k < m2; k++)
					if  (user_Read(buf, id_MAP, count, upper_limit)) mp[i++] = user_strtod(buf); else goto err;
					for (k = 0; k < m3; k++)
					if  (user_Read(buf, id_MAP, count, upper_limit)) mp[i++] = user_strtod(buf); else goto err;
					for (k = 0; k < m4*m5; k++)
					for (m = 0; m < m1;    m++)
					if  (user_Read(buf, id_MAP, count, upper_limit)) mp[i++] = user_strtod(buf); else goto err;
         }		break;
         case  B1POLY_CELL: {
					m1 = (int)mp[7]/2;
					m3 = (int)mp[9];

					for (k = 0; k < m3; k++)
					for (m = 0; m < m1; m++)
					if  (user_Read(buf, id_MAP, count, upper_limit)) mp[i++] = user_strtod(buf); else goto err;
         }		break;
         default: {
				for (m = 0; m < k; m++)
				if  (user_Read(buf, id_MAP, count, upper_limit)) mp[i++] = user_strtod(buf); else goto err;
			}
		}
	}
err:
	return mp;
}

Topo * topo_in(char * id_GRAF, unsigned long & count, unsigned long upper_limit, int N)
{
	Topo * graph = NULL;
	if (id_GRAF && count < upper_limit) {
		char buf[BUF_SIZE+1];
		int  m0, m1, k;
		if ( user_Read(buf, id_GRAF, count, upper_limit)) m0 = atoi(buf); else goto err;
		if ( user_Read(buf, id_GRAF, count, upper_limit)) m1 = atoi(buf); else goto err;

		if ((graph    = new_struct<Topo>(m1+2)) != NULL)
		for (graph[0] = m0, graph[1] = m1, k = 2; k < m1+2; k++)
		if ( user_Read(buf, id_GRAF, count, upper_limit)) graph[k] = atoi(buf)-N; else goto err;
	}
err:
	return graph;
}

///////////////////////////////////////////////////////////////////////////
//...writing geometrical card and topology in ASCII code (working program);
void map_out(FILE * id_MAP, CMap * mp, int id_long)
{
	if (mp) {
		int  k = size_of_map(mp)+1,
			  m = size_of_dop(mp)+k, i = 0, m1, m2, m3, m4, m5;
		const char * ch_format = id_long ? "%0.15lg " : "%lg ";

      for (; k > 0; k--) fprintf(id_MAP, ch_format, mp[i++]);
      fprintf(id_MAP, "\n");

      switch ((int)mp[i-1]) {
			case      NULL_CELL: return;
			case  B1SPLINE_CELL: {
					m1 = (int)mp[7]/2;
					m2 = (int)mp[8]+2;
					m3 = (int)mp[9];

					for (fprintf(id_MAP, "    "), k = 0; k < m2; k++) fprintf(id_MAP, ch_format, mp[i++]);
					fprintf(id_MAP, "\n");

					for (k = 0; k < m3; k++) {
						for (fprintf(id_MAP, "        "), m = 0; m < m1; m++) fprintf(id_MAP, ch_format, mp[i++]);
						fprintf(id_MAP, "\n");
					}
			}		break;
			case  B2SPLINE_CELL: {
					m1 = (int)mp[7]/4;
					m2 = (int)mp[8]+2;
					m3 = (int)mp[9]+2;
					m4 = (int)mp[10];
					m5 = (int)mp[11];

					for (fprintf(id_MAP, "    "), k = 0; k < m2; k++) fprintf(id_MAP, ch_format, mp[i++]);
					fprintf(id_MAP, "\n");

					for (fprintf(id_MAP, "    "), k = 0; k < m3; k++) fprintf(id_MAP, ch_format, mp[i++]);
					fprintf(id_MAP, "\n");

					for (k = 0; k < m4*m5; k++) {
						for (fprintf(id_MAP, "        "), m = 0; m < m1; m++) fprintf(id_MAP, ch_format, mp[i++]);
						fprintf(id_MAP, "\n");
					}
         }		break;
         case  B1POLY_CELL: {
					m1 = (int)mp[7]/2;
					m3 = (int)mp[9];

					for (k = 0; k < m3; k++) {
						for (fprintf(id_MAP, "        "), m = 0; m < m1; m++) fprintf(id_MAP, ch_format, mp[i++]);
						fprintf(id_MAP, "\n");
					}
         }		break;
         default: {
					for (fprintf(id_MAP, "    "); i < m; i++) fprintf(id_MAP, ch_format, mp[i]);
					fprintf(id_MAP, "\n");
			}
		}
	}
	return;
}

void topo_out(FILE * id_GRAF, Topo * graph)
{
	if (graph) {
		for (int k = graph[1]+2; k > 0; k--, graph++)
		fprintf(id_GRAF, " %d", graph[0]);
		fprintf(id_GRAF, "\n");
	}
	return;
}

////////////////////////////////////////////////////////////
//...exchange data of Bar structure with disk in ASCII code;
void bar_out(char * ch_BAR, Bar * bar)
{
	if (ch_BAR && bar && ! bar->mp) {
		FILE  *  id_BAR = fopen(ch_BAR, "w");
      fprintf (id_BAR, "0:| # Bar %d\n", bar->graph[0]);

		for (int k = 0; k < bar->graph[0]; k++) {
			fprintf (id_BAR, "%d:| ", k);  map_out(id_BAR, bar->ce[k]->mp);
			fprintf (id_BAR, "%d:| ", k); topo_out(id_BAR, bar->ce[k]->graph);

/////////////////////////////////////
//...writing cell in parameter space;
			for (int m = bar->ce[k]->pm && bar->ce[k]->graph ? bar->ce[k]->graph[1] : 0,
						l = 0; l < m; l++) {
				fprintf(id_BAR, "%d:| ", k);  map_out(id_BAR, bar->ce[k]->pm[l]);
			}
		}
		fclose(id_BAR);
	}
	return;
}

///////////////////////////////
//...dimention of body surface;
int bar_dim(Bar * bar)
{
	int dim = ERR_DIM;
	for (int k = 0; k < bar->graph[0]; k++)
		dim = max(dim, cells_dim((Cells *)bar->ce[k]));
	return dim;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//...normalization of geometrical card (i.e. extracting reference point and tree coordinate axis);
void map_normalizat(CMap * mp, double * P, double * ort)
{
	if (mp && ort) {
		ort[0] = ort[1] = ort[3] = ort[5] = ort[7] = ort[8] = 0.;
		ort[2] = ort[4] = ort[6] = 1.;

		if (ID_MAP(0, NULL_GENUS) != mp[0]) {
			double CZ = cos(mp[4]), SZ = sin(mp[4]), CY = cos(mp[5]), SY = sin(mp[5]),
					 CX = cos(mp[6]), SX = sin(mp[6]);
			point_iso<double>(ort,   NULL, CZ, SZ, CY, SY, CX, SX);
			point_iso<double>(ort+3, NULL, CZ, SZ, CY, SY, CX, SX);
			point_iso<double>(ort+6, NULL, CZ, SZ, CY, SY, CX, SX);
		}
		if (P) {
			P[0]  = mp[1];
			P[1]  = mp[2];
			P[2]  = mp[3];
			P[3]  = P[0]*ort[0]+P[1]*ort[1]+P[2]*ort[2];
			P[0] -= P[3]*ort[0];
			P[1] -= P[3]*ort[1];
			P[2] -= P[3]*ort[2];
		}
	}
	return;
}

///////////////////////////////////////////////////////
//...normalization of geometrical card for plane facet;
void map_normalizat_facet(CMap * mp, double * P, double * ort)
{
	if (mp && ort) {
		ort[0] = ort[1] = ort[3] = ort[5] = ort[7] = ort[8] = 0.;
		ort[2] = ort[4] = ort[6] = 1.;

		if (ID_MAP(0, NULL_GENUS) != mp[0]) {
			ort[0] = mp[4];
			ort[1] = mp[5];
			ort[2] = mp[6];
		}
		if (P) {
			P[0]  = mp[1];
			P[1]  = mp[2];
			P[2]  = mp[3];
			P[3]  = P[0]*ort[0]+P[1]*ort[1]+P[2]*ort[2];
			P[0] -= P[3]*ort[0];
			P[1] -= P[3]*ort[1];
			P[2] -= P[3]*ort[2];
		}
	}
	return;
}

///////////////////////////////////////////////////////////////////////////////////////
//...displacement and rotate of geometrical card (and calculation of new Euler angles);
void map_iso(CMap * mp, double * P, double & CZ, double & SZ,
                                    double & CY, double & SY, double & CX, double & SX)
{
	if (mp) {
		point_iso(mp+1, P, CZ, SZ, CY, SY, CX, SX);

		if (ID_MAP(0, NULL_GENUS) != mp[0]) {
			double ort[9];
			map_normalizat(mp, NULL, ort);
			point_iso<double>(ort, NULL, CZ, SZ, CY, SY, CX, SX);

			mp[4] = arg0(comp(ort[0], ort[1]));
			mp[5] = arg0(comp(ort[2], sqrt(ort[0]*ort[0]+ort[1]*ort[1])));

			point_iso<double>(ort+6, NULL, CZ, SZ, CY, SY, CX, SX);
			point_iso<double>(ort+6, NULL, 1., 0., cos(mp[5]), sin(-mp[5]), cos(mp[4]), sin(-mp[4]));

			mp[6] = arg0(comp(ort[6], ort[7]));
		}
	}
	return;
}

///////////////////////////////////////////////////////////////////////
//...displacement and rotate of geometrical card with angles in radian;
void map_iso(CMap * mp, double * P, double fi, double theta, double fX)
{
	if (mp) {
		double CZ = cos(fi), SZ = sin(fi), CY = cos(theta), SY = sin(theta),
				 CX = cos(fX), SX = sin(fX);
		map_iso(mp, P, CZ, SZ, CY, SY, CX, SX);
	}
	return;
}

///////////////////////////////////////////////////////////////
//...displacement and rotate of the cell (recurrent procedure);
void cells_iso(Cells * ce, double * P, double & CZ, double & SZ,
                                       double & CY, double & SY, double & CX, double & SX)
{
	map_iso(ce->mp, P, CZ, SZ, CY, SY, CX, SX);
	for (int k = ce->graph[1]+1; k > 1; k--)
		cells_iso((Cells *)ce->ce[ce->graph[k]], P, CZ, SZ, CY, SY, CX, SX);
	return;
}

/////////////////////////////////////////
//...displacement and rotate of the cell;
void cells_iso(Cells * ce, double * P, double fi, double theta, double fX)
{
	if (ce && ce->mp) {
		double CZ = cos(fi), SZ = sin(fi), CY = cos(theta), SY = sin(theta),
				 CX = cos(fX), SX = sin(fX);
		for (int k = 0; k <= ce->graph[0]; k++) map_iso(ce->ce[k]->mp, P, CZ, SZ, CY, SY, CX, SX);
	}
	return;
}

////////////////////////////////////////////
//...displacement and rotate of the surface;
void bar_iso(Bar * bar, double * P, double fi, double theta, double fX)
{
	if (bar) {
		double CZ = cos(fi), SZ = sin(fi), CY = cos(theta), SY = sin(theta),
				 CX = cos(fX), SX = sin(fX);
		for (int k = 0; k < bar->graph[0]; k++) map_iso(bar->ce[k]->mp, P, CZ, SZ, CY, SY, CX, SX);
	}
	return;
}

////////////////////////////////////////////////////////////////////////
//...transformation to standerd position, when nulling geometrical card;
void cells_to_origine(Cells * ce)
{
	if (ce && ce->mp) {
		double P[3] = {-ce->mp[1], -ce->mp[2], -ce->mp[3]};

		cells_iso(ce, P);
		if (ID_MAP(0, NULL_GENUS) != ce->mp[0])
		cells_iso(ce, NULL, -ce->mp[6], -ce->mp[5], -ce->mp[4]);
	}
	return;
}

/////////////////////////////////////////////////////
//...transformation point in local coordinate system;
void make_local(double * P, CMap * mp, int id_shift)
{
	if (P && mp) {
		if (id_shift) point_shift(P, mp+1, OK_STATE);
		if (ID_MAP (0, NULL_GENUS) != mp[0])
		point_iso  (P, NULL, -mp[6], -mp[5], -mp[4]);
	}
	return;
}

///////////////////////////////////////////////////////////////////
//...vise verse - transformation point in common coordinate system;
void make_common(double * P, CMap * mp, int id_shift)
{
	if (P && mp) {
		if (ID_MAP (0, NULL_GENUS) != mp[0])
		point_iso  (P, NULL, mp[4], mp[5], mp[6]);
		if (id_shift) point_shift(P, mp+1);
	}
	return;
}

////////////////////////////////////////////////////////////////
//...transformation geometrical card in local coordinate system;
void map_local(CMap * mp, CMap * ext_mp)
{
	if (mp && ext_mp) {
		point_shift(mp+1, ext_mp+1, OK_STATE);
		if (ID_MAP (0,  NULL_GENUS) != ext_mp[0])
		map_iso    (mp, NULL, -ext_mp[6], -ext_mp[5], -ext_mp[4]);
	}
	return;
}

/////////////////////////////////////////////////////////////////
//...transformation geometrical card in common coordinate system;
void map_common(CMap * mp, CMap * ext_mp)
{
	if (mp && ext_mp) {
		if (ID_MAP (0,  NULL_GENUS) != ext_mp[0])
		map_iso    (mp, NULL, ext_mp[4], ext_mp[5], ext_mp[6]);
		point_shift(mp+1, ext_mp+1);
	}
	return;
}

////////////////////////////////////////////////////
//...transformation cell in local coordinate system;
void cells_local(Cells * ce, CMap * ext_mp)
{
	if (ce && ext_mp) {
		double P[3] = {-ext_mp[1], -ext_mp[2], -ext_mp[3]};

		cells_iso (ce, P);
		if (ID_MAP(0,  NULL_GENUS) != ext_mp[0])
		cells_iso (ce, NULL, -ext_mp[6], -ext_mp[5], -ext_mp[4]);
	}
	return;
}

/////////////////////////////////////////////////////
//...transformation cell in common coordinate system;
void cells_common(Cells * ce, CMap * ext_mp)
{
	if (ce && ext_mp) {
		double P[3] = {ext_mp[1], ext_mp[2], ext_mp[3]};

		if (ID_MAP(0,  NULL_GENUS) != ext_mp[0])
		cells_iso (ce, NULL, ext_mp[4], ext_mp[5], ext_mp[6]);
		cells_iso (ce, P);
	}
	return;
}

///////////////////////////////////////////////////////
//...transformation surface in local coordinate system;
void bar_local(Bar * bar, CMap * ext_mp)
{
	if (bar && ext_mp) {
		double P[3] = {-ext_mp[1], -ext_mp[2], -ext_mp[3]};

		bar_iso  (bar, P);
		if (ID_MAP(0,  NULL_GENUS) != ext_mp[0])
		bar_iso  (bar, NULL, -ext_mp[6], -ext_mp[5], -ext_mp[4]);
	}
	return;
}

////////////////////////////////////////////////////////
//...transformation surface in common coordinate system;
void bar_common(Bar * bar, CMap * ext_mp)
{
	if (bar && ext_mp) {
		double P[3] = {ext_mp[1], ext_mp[2], ext_mp[3]};

		if (ID_MAP(0,  NULL_GENUS) != ext_mp[0])
		bar_iso  (bar, NULL, ext_mp[4], ext_mp[5], ext_mp[6]);
		bar_iso  (bar, P);
	}
	return;
}

////////////////////////////////////////////////////////////////////////
//...auxilliary procedure - transformation geometrical card to new kind;
void map_to(CMap *& mp, int dim, int genus)
{
	int m = size_of_map(mp), mmm = size_of_map(dim, genus);
	if (m < mmm) {
		CMap * new_map = new_struct<CMap>(mmm+1);
		if ( ! new_map) return;
		else {
			memcpy(new_map, mp, m*sizeof(CMap));
			delete_struct(mp); mp = new_map;
		}
	}
	if (mp) {
		mp[0]   = (CMap)ID_MAP(dim, genus);
		mp[mmm] = (CMap)NULL_CELL;
	}
	return;
}

///////////////////////////////////////////////////////////////////////////
//...auxilliary function - identification of the point in geometrical card;
int map_into(CMap * mp, double X, double Y, double Z, int id_special) //...for blocks intersection;
{
	int k;
	switch ((int)mp[k = size_of_map(mp)]*(id_special == 1 ? 1 : 0)) {
		case SPECIAL_SHEET : return(fabs(X+mp[1]-mp[k+3]) < mp[k+1]+EE_ker && fabs(Y+mp[2]-mp[k+4]) < mp[k+2]+EE_ker);
		case SPECIAL_BOX   : return(mp[7] > 0. &&  fabs(X+mp[1]-mp[k+4]) < mp[k+1]+EE_ker && fabs(Y+mp[2]-mp[k+5]) < mp[k+2]+EE_ker &&
															  fabs(Z+mp[3]-mp[k+6]) < mp[k+3]+EE_ker ||
										  mp[7] < 0. && (fabs(X+mp[1]-mp[k+4]) > mp[k+1]-EE_ker || fabs(Y+mp[2]-mp[k+5]) > mp[k+2]-EE_ker ||
															  fabs(Z+mp[3]-mp[k+6]) > mp[k+3]-EE_ker));
		default            : return(ID_MAP(2, SPHERE_GENUS) == mp[0] && (mp[7] > 0. && X*X+Y*Y+Z*Z < mp[7]*mp[7]+EE_ker  ||
																							mp[7] < 0. && X*X+Y*Y+Z*Z > mp[7]*mp[7]-EE_ker) ||
										  ID_MAP(2,    CYL_GENUS) == mp[0] &&                X*X+Y*Y     < mp[7]*mp[7]+EE_ker  &&
																			 (mp[8]-EE_ker < Z && Z   < mp[9]+EE_ker       ||
																			  X*X+Y*Y+(Z-mp[8])*(Z-mp[8]) < mp[7]*mp[7]+EE_ker ||
																			  X*X+Y*Y+(Z-mp[9])*(Z-mp[9]) < mp[7]*mp[7]+EE_ker));
	}
}

/////////////////////////////////////////////////////////////////////////////////
//...point identification in geometrical card without ending connecting elements;
int map1into(CMap * mp, double X, double Y, double Z) //...for blocks visualization;
{
	int k;
	switch ((int)mp[k = size_of_map(mp)]) {
		case   BASIS_SHEET :
		case SPECIAL_SHEET : return(fabs(X+mp[1]-mp[k+3]) < mp[k+1]+EE_ker && fabs(Y+mp[2]-mp[k+4]) < mp[k+2]+EE_ker);
		case SPECIAL_BOX   :
		case   BASIS_BOX   : return(mp[7] > 0. &&  fabs(X+mp[1]-mp[k+4]) < mp[k+1]+EE_ker && fabs(Y+mp[2]-mp[k+5]) < mp[k+2]+EE_ker &&
															  fabs(Z+mp[3]-mp[k+6]) < mp[k+3]+EE_ker ||
										  mp[7] < 0. && (fabs(X+mp[1]-mp[k+4]) > mp[k+1]-EE_ker || fabs(Y+mp[2]-mp[k+5]) > mp[k+2]-EE_ker ||
															  fabs(Z+mp[3]-mp[k+6]) > mp[k+3]-EE_ker));
		default            : return(ID_MAP(2, SPHERE_GENUS) == mp[0] && (mp[7] > 0. && X*X+Y*Y+Z*Z < mp[7]*mp[7]+EE_ker  ||
																							mp[7] < 0. && X*X+Y*Y+Z*Z > mp[7]*mp[7]+EE_ker) ||
										  ID_MAP(2,    CYL_GENUS) == mp[0] &&                X*X+Y*Y     < mp[7]*mp[7]+EE_ker  &&
										  mp[8]-EE_ker < Z && Z < mp[9]+EE_ker);
	}
}

///////////////////////////////////////////////////////////////////
//...identification of the cards (for inverse orientation id = -1);
int map_id(CMap * mp, CMap * ext_mp)
{
	if (! mp || ! ext_mp) return(0);
	int id = mp[0] == ext_mp[0], m, k;

///////////////////////////////////////////////////////////
//...compare orientation and position of geometrical cards;
	if (id) {
		double P[4], ext_P[4], ort[9], ext_ort[9];
		if (mp[size_of_map(mp)] == (CMap)FACET_CELL) map_normalizat_facet(mp, P, ort);
		else                                         map_normalizat      (mp, P, ort);

		if (ext_mp[size_of_map(ext_mp)] == (CMap)FACET_CELL) map_normalizat_facet(ext_mp, ext_P, ext_ort);
		else                                                 map_normalizat      (ext_mp, ext_P, ext_ort);

///////////////////////////////////
//...compare orientation of Z axis;
     if (abs(id) && ID_MAP(2,  SPHERE_GENUS) != mp[0]) //...сфера;
     id *= (max(fabs(ort[0]-ext_ort[0]),
            max(fabs(ort[1]-ext_ort[1]), fabs(ort[2]-ext_ort[2]))
           ) < EE_ker)-
           (max(fabs(ort[0]+ext_ort[0]),
            max(fabs(ort[1]+ext_ort[1]), fabs(ort[2]+ext_ort[2]))
           ) < EE_ker);

     if (abs(id) && (ID_MAP(2,       CONE_GENUS) == mp[0] ||  //...cone;
                     ID_MAP(2,   ELL_CONE_GENUS) == mp[0] ||  //...elliptice cone;
                     ID_MAP(2,      PARAB_GENUS) == mp[0] ||  //...paraboloid;
                     ID_MAP(2,  ELL_PARAB_GENUS) == mp[0] ||  //...elliptic paraboloid;
                     ID_MAP(2,     BI_HYP_GENUS) == mp[0] ||  //...two-sheeted hyperboloid;
                     ID_MAP(2, ELL_BI_HYP_GENUS) == mp[0] ||   //...elliptic two-sheeted hyperboloid;
                                    NURBS_GENUS  == map_genus(mp) ||   //...Nurbs;
                                   BEZIER_GENUS  == map_genus(mp)))    //...Bezier;
     id *= id > 0;

////////////////////////////////////////
//...compare position ofreference point;
     if (abs(id) && ID_MAP(1,    NULL_GENUS) != mp[0]  //...line;
                 && ID_MAP(2,     CYL_GENUS) != mp[0]  //...cylinder;
                 && ID_MAP(2, ELL_CYL_GENUS) != mp[0]) //...elliptic cylinder;
     id *= (fabs(P[3]-id*ext_P[3]) < EE_ker);

     if (abs(id) && ID_MAP(2,    NULL_GENUS) != mp[0]) //...plane;
     id *= max(fabs(P[0]-ext_P[0]),
           max(fabs(P[1]-ext_P[1]), fabs(P[2]-ext_P[2]))
           ) < EE_ker;

/////////////////////////////////////////
//...compare orientation of the axis Y  and axis X;
     if (abs(id) && ID_MAP(1,        CYL_GENUS) == mp[0] ||  //...ellipse;
                    ID_MAP(1,       CONE_GENUS) == mp[0] ||  //...hyperbola;
                    ID_MAP(1,   ELL_CONE_GENUS) == mp[0] ||  //...parabola;
                    ID_MAP(2,  ELLIPSOID_GENUS) == mp[0] ||  //...ellipsoid;
                    ID_MAP(2,  ELL_PARAB_GENUS) == mp[0] ||  //...elliptic paraboloid;
                    ID_MAP(2, ELL_HYPERB_GENUS) == mp[0] ||  //...elliptic hyperboloid;
                    ID_MAP(2, ELL_BI_HYP_GENUS) == mp[0]) {  //...elliptic two-sheeted hyperboloid;

         id *= (max(fabs(ort[3]-ext_ort[3]),
                max(fabs(ort[4]-ext_ort[4]), fabs(ort[5]-ext_ort[5]))
               ) < EE_ker)-
               (max(fabs(ort[3]+ext_ort[3]),
                max(fabs(ort[4]+ext_ort[4]), fabs(ort[5]+ext_ort[5]))
               ) < EE_ker);

         if (abs(id) && (ID_MAP(1,       CONE_GENUS) == mp[0] ||  //...hyperbola;
                         ID_MAP(1,   ELL_CONE_GENUS) == mp[0] ||  //...parabola;
                         ID_MAP(2,  ELL_PARAB_GENUS) == mp[0] ||  //...elliptic paraboloid;
                         ID_MAP(2, ELL_BI_HYP_GENUS) == mp[0] ||   //...elliptic two-sheeted hyperboloid;
                                        NURBS_GENUS  == map_genus(mp) ||   //...Nurbs;
                                       BEZIER_GENUS  == map_genus(mp)))    //...Bezier;
         id *= id > 0;
     }

//////////////////////////////////////////////////////////////////////
//...compare additional parameters in common part of geometrical card;
     if (abs(id) && BEZIER_GENUS != map_genus(mp))  //...Bezier;
         for (k = 7; id && k < (m = size_of_map(map_dim(mp), map_genus(mp))); k++)
         id *= fabs(fabs(mp[k])-fabs(ext_mp[k])) < EE_ker; //...по абсолютной величине!

//////////////////////////////////////////////////////////////////////////
//...more strictly verification for circle (verify coincidence of X axis);
     if (abs(id) && ID_MAP(1,  SPHERE_GENUS) == mp[0]) //...окружность;
     id *= (max(fabs(ort[6]-ext_ort[6]),
            max(fabs(ort[7]-ext_ort[7]), fabs(ort[8]-ext_ort[8]))
           ) < EE_ker);

////////////////////////////////////////////////////////////////////////
//...compare additional parameters for dimension cells (for Nurbs also);
     if (abs(id) && mp[m]  ==  ext_mp[m] && mp[m] != (CMap)FACET_CELL) {
         for (k = size_of_dop((int)mp[m], mp)+m; id && k > m; k--)
         id *= fabs(mp[k]-ext_mp[k]) < EE_ker;
     }

//////////////////////////////////////////////
//...compare additional parameters for Bezier;
     if (abs(id) && BEZIER_GENUS == map_genus(mp))  //...Bezier;
     id *= (int)mp[7] == (int)ext_mp[7];
  }
  return id;
}

///////////////////////////////////////////////////////////
//...cells identification (for reversing cells -- id = -1);
int cells_id(Cells * ce, Cells * ext_ce, int id_num_correct)
{//...rule: cekk with negative boundary in the structure -- identifying,
 //...with negative boundary out of the structure -- not identifying;
	if (! ext_ce	  || ! ce)		return(1);
	if (! ext_ce->mp || ! ce->mp) return(0);

/////////////////////////////
//...check of the cell boundary;
	int id = ce->graph[1] == ext_ce->graph[1], i, k, l, m_ce;
	for (k = ce->graph[1]+1; id && k > 1; k--) if (ext_ce->graph[k] <= -ce->graph[0]-1) id = 0;

////////////////////////////////////////
//...check geometrical card of the cell;
	if (id)
      id *= map_id (ce->mp, ext_ce->mp);

/////////////////////////////////////////////////////////////////////////
//...compare elements of cell boundary and marking of coincided elements;
	if (id) {
		for (l = ext_ce->graph[1]+1; id && l > 1; l--) {
			for (k = ce->graph[1]+1; k > 1; k--) if (ce->graph[k] >= 0)
				if ((m_ce  = ext_ce->graph[l]) < 0 && (i = (m_ce+ce->graph[k]+1) == 0) != 0 ||
					  m_ce >= 0 && (i = cells_id(ce->ce[ce->graph[k]], ext_ce->ce[m_ce])) != 0) break;

			id *= abs(i);
			if (k > 1) { //...marking boundary and coincided elements;
				ce->graph[k] = -ce->graph[k]-1;
				if (id_num_correct && m_ce >= 0) {
					topo_correct(ext_ce->ce[ext_ce->graph[0]], m_ce, ce->graph[k], NULL_STATE);
					if (m_ce >= 0) ext_ce->graph[l] = m_ce;
				}
			}
		}

/////////////////////////////////////
//...reconstruction marking elements;
		for (k = ce->graph[1]+1; k > 1; k--) if (ce->graph[k] < 0)
					ce->graph[k] = -ce->graph[k]-1;
	}
	return id;
}

///////////////////////////////////////////////
//...cells identification in the cell boundary;
int common_id(Cells * ce, Cells * ext_ce)
{
  if (! ext_ce	|| ! ce)							  return(1);
  if (ce->ce != ext_ce->ce || ext_ce == ce) return(0);
  int id, k;

/////////////////////////////////////////////////////
//...comparison of the elements in the cell boundary;
  for (   id = 0,  k = ce->graph[1]+1; ! id && k > 1; k--)
  for (int l = ext_ce->graph[1]+1; ! id && l > 1; l--) id = ce->graph[k] == ext_ce->graph[l];
  return id*k;
}

///////////////////////////////////////////////////////////////////////////////////////////
//...identification of the cell number into the structure (with opposite sign for id = -1);
int cells_in_struct(Cells * ce, Cells * ext_ce, int id_num_correct)
{//...rule: null cell -- identifying automatically;

	if (! ext_ce || ! ce) return(1);
	if (! ext_ce->mp)		 return(0);

	int N, i;
	if (! ce->mp)
	for ( N = i = 0; i < ce->graph[0]; i++) {
			N = i; if ((++N *= cells_id(ce->ce[i], ext_ce, id_num_correct)) != 0) break;
	}
	else {
		for ( N = i = 0; i < ce->graph[1]; i++) {
				N = ce->graph[i+2]; if ((++N *= cells_id(ce->ce[ce->graph[i+2]], ext_ce, id_num_correct)) != 0) break;
		}
		if (! N)
			for (i = 0; i < ce->graph[1]; i++)
				if ((N = cells_in_struct(ce->ce[ce->graph[i+2]], ext_ce, id_num_correct)) != 0) break;
	}
	return N;
}

//////////////////////////////////////////////////////////////////////////////////
//...indentification cell number into the surface (with negativ sign for id = -1);
int cells_in_bar(Bar * bar, Cells * ce)
{
	if (! bar || bar->mp || ! ce || ! ce->mp) return(0);

	int N, i;
	for (N = 0, i = 0; i < bar->graph[0]; i++) {
		 N = i; if (++N *= cells_id(bar->ce[i], ce)) break;
	}
	return N;
}

///////////////////////////////////////////////////////////////////////
//...auxilliary function -- construction list of substructure elements;
void search_cells_element(Cells * ce, int *& id, int & N_buf, int & buf_size, int buf_delta)
{
	if (ce && ce->graph)
	for ( int k = 0; k < ce->graph[1]; k++) {
			if (buf_size == N_buf) {
				int * new_id = new_struct<int>(buf_size += buf_delta);
				memcpy(new_id, id, N_buf*sizeof(int));
				delete_struct(id); id = new_id;
			}
			id[N_buf++] = ce->graph[k+2];
			search_cells_element(ce->ce[ce->graph[k+2]], id, N_buf, buf_size, buf_delta);
	}
}

//////////////////////////////////////////////
//...arc length of one-dimensional space cell;
double cells_length(Cells * ce)
{
	double length = 0.;
	if (! ce || ! ce->mp) return(length);

	if (ce->mp[0] == ID_MAP(1, NULL_GENUS)) {
		length = abs_point(ce->ce[get_num(ce->graph, 0)]->mp+1,
								 ce->ce[get_num(ce->graph, 1)]->mp+1);
	}
	else
	if (ce->mp[0] == ID_MAP(1, SPHERE_GENUS)) {
		length = fabs(ce->mp[7]*2.*ce->mp[8]);
	}
	else
	if (ce->mp[0] == ID_MAP(1, CYL_GENUS)) {
		length = 0.;
	}
	else
	if (ce->mp[0] == ID_MAP(1, CONE_GENUS)) {
		length = 0.;
	}
	else
	if (ce->mp[0] == ID_MAP(1, ELL_CONE_GENUS)) {
		length = 0.;
	}
	return(length);
}

////////////////////////////////////////////////////////////////
//...correction of local coordinate system;
void edge_correct(Cells * ce)
{
	if (ce && ce->mp && ce->graph[1] > 0)
	if (ce->mp[0] == ID_MAP(1, NULL_GENUS)) {
		ce->mp[1] = (CMap)(ce->ce[get_num(ce->graph, 0)]->mp[1]+
								 ce->ce[get_num(ce->graph, 1)]->mp[1])*.5;
		ce->mp[2] = (CMap)(ce->ce[get_num(ce->graph, 0)]->mp[2]+
								 ce->ce[get_num(ce->graph, 1)]->mp[2])*.5;
		ce->mp[3] = (CMap)(ce->ce[get_num(ce->graph, 0)]->mp[3]+
								 ce->ce[get_num(ce->graph, 1)]->mp[3])*.5;
	}
	else
	if (ce->mp[0] == ID_MAP(1, SPHERE_GENUS)) {
		int     m1    = get_num(ce->graph, 0),
				  m2    = get_num(ce->graph, 1);
		double  CZ    = cos(ce->mp[6]), SZ = -sin(ce->mp[6]),
				  CY    = cos(ce->mp[5]), SY = -sin(ce->mp[5]),
				  CX    = cos(ce->mp[4]), SX = -sin(ce->mp[4]), f1, f2,
				  p1[3] = { ce->ce[m1]->mp[1]-ce->mp[1],
								ce->ce[m1]->mp[2]-ce->mp[2],
								ce->ce[m1]->mp[3]-ce->mp[3] },
				  p2[3] = { ce->ce[m2]->mp[1]-ce->mp[1],
								ce->ce[m2]->mp[2]-ce->mp[2],
								ce->ce[m2]->mp[3]-ce->mp[3] };
		complex z1, z2;

		point_iso<double>(p1, NULL, CZ, SZ, CY, SY, CX, SX); f1 = arg2(z1 = comp(p1[0], p1[1]));
		point_iso<double>(p2, NULL, CZ, SZ, CY, SY, CX, SX); f2 = arg2(z2 = comp(p2[0], p2[1]));
		if (z1 == z2 || f2+EE_ker < f1) f2 += 2.*M_PI;

		ce->mp[6] += .5*(f1+f2);
		ce->mp[8]  = .5*fabs(f2-f1);
	}
	else
   if (ce->mp[0] == ID_MAP(1, CYL_GENUS)) {
	}
	else
   if (ce->mp[0] == ID_MAP(1, CONE_GENUS)) {
	}
	else
   if (ce->mp[0] == ID_MAP(1, ELL_CONE_GENUS)) {
   }
	return;
}

////////////////////////////////
//...ordering graph of topology;
void topo_ord(Topo * graph)
{
	for (int k = graph[1]+1; k > 1; k--)
	for (int l = k-1; l > 1; l--) if (graph[k] < graph[l]) swap(graph[k], graph[l]);
	return;
}

/////////////////////////////////////////////////////////
//...increasing size of topologe graph (working program);
int topo_pad(Topo *& graph, int N_pad)
{
	if (N_pad > 0) {
		int m = graph[1]+2;

		Topo * new_graph = new_struct<Topo>(m+N_pad);
		if (! new_graph) return(0);

		while (--m >= 2) new_graph[m+N_pad] = graph[m]; ++m;
		while (--m >= 0) new_graph[m] = graph[m];

		delete_struct(graph); graph = new_graph;
		graph[1] += N_pad;
	}
	return(1);
}

////////////////////////////////////////////////////////////////
//...addition new element in graph topology (working procedure);
int topo_insert(Topo *& graph, Topo new_element)
{
	if (topo_pad(graph, 1)) graph[2] = new_element;
	else                    return(0);
	return(1);
}

//////////////////////////////////////////////////////////////////////////
//...identification of element in one list of linking (working procwdure);
int topo_id(Topo * graph, Topo some_element)
{
	int j = graph[1]+1;

	while (1 < j && some_element != graph[j]) j--;
	return(1 < j);
}

////////////////////////////////////////////////////////////////////
//...deleting elements from one list of linking (working procedure);
void topo_exc(Topo * graph, Topo some_element)
{
	int j = graph[1]+1;

	while (1 < j && some_element != graph[j]) j--;
	if    (1 < j) for (graph[0]--, graph[1]--; j < graph[1]+2; j++) graph[j] = graph[j+1];

	return;
}

//////////////////////////////////////////////////////////
//...identification elementin list of linking of the cell;
int topo_id(Cells * ce, Topo some_element, Topo exc_element)
{
	if (! ce || ! ce->mp) return(0);
	int id = topo_id(ce->graph, some_element), k;

///////////////////////////////////////
//...compare      lists of subelements;
	for (k = ce->graph[1]+1; ! id && k > 1; k--)
	if (ce->graph[k] != exc_element) id = topo_id(ce->ce[ce->graph[k]], some_element, exc_element);
	return id;
}

////////////////////////////////////////
//...correction links in topology graph;
void topo_correct(Cells * ce, int N, int N_new, int id_list)
{
	int k, i;
	if (! ce) return;
	if (id_list || ! ce->mp) { //...correcting or installingnumber by common list;
		for (k = ce->graph[0]; k >= 0; k--) if (ce->ce[k])
		for (i = ce->ce[k]->graph[1]+1; i > 1; i--)
			if (N == -1 && N_new  == -1 && ce->ce[k]->graph[i] < 0) ++(ce->ce[k]->graph[i]) *= -1; else
			if (N ==  ce->ce[k]->graph[i]) ce->ce[k]->graph[i] = N_new;
	}
	else
		if (N == -1 && N_new == -1) { //...correcting number (for cell) by structure;
			for (i = 0; i < ce->graph[1]; i++)
				if (ce->graph[i+2] < 0) ++(ce->graph[i+2]) *= -1;
				else   topo_correct(ce->ce[ce->graph[i+2]], N, N_new);
		}
		else //...installing number (for cell) by structure;
			for (i = 0; i < ce->graph[1]; i++)
			if (ce->graph[i+2] >= 0) {
				if (ce->graph[i+2] == N) ce->graph[i+2] = N_new;
				else topo_correct(ce->ce[ce->graph[i+2]], N, N_new);
			}
}

////////////////////////////////////////////////////////////////////////////
//...permutation of two element in common list of surface (working program);
void bar_mutat(Bar * bar, int k, int m)
{
	if (! bar || bar->mp || bar->graph[0] <= k || k < 0 || bar->graph[0] <= m || m < 0) return;

	(bar->graph[0])--;

	topo_correct((Cells *)bar, k, -m-1);
	topo_correct((Cells *)bar, m, -k-1);
	topo_correct((Cells *)bar);

	bar->graph[0]++;

	Bar * temp_bar = bar->ce[k]; bar->ce[k] = bar->ce[m]; bar->ce[m] = temp_bar;
	return;
}

///////////////////////////////////////////////////////////////////////////////
//...permutation of two element in common list of the cell (working procedure);
void cells_mutat(Cells * ce, int k, int m)
{
	if (! ce || ! ce->mp || ce->graph[0] <= k || k < 0 || ce->graph[0] <= m || m < 0) return;

	topo_correct(ce, k, -m-1);
	topo_correct(ce, m, -k-1);
	topo_correct(ce);

	Cells * temp_ce = ce->ce[k]; ce->ce[k] = ce->ce[m]; ce->ce[m] = temp_ce;
	return;
}

/////////////////////////////////////////////////////////////////////////////////
//...ordering elemnts in common list of surface and correction lists of topology;
void bar_ord(Bar * bar)
{
	if (bar && ! bar->mp)
		for (int  j = 0, k, dim = bar_dim(bar); dim >= 0; dim--)
			for (k = j; k < bar->graph[0]; k++)
				if (dim == cells_dim(bar->ce[k])) bar_mutat(bar, k, j++);
	return;
}

///////////////////////////////////////////////////////////////////////////////////////////////
//...ordering elements of common list (not including in the list) and correction link topology;
void cells_ord(Cells * ce)
{
	if (ce && ce->mp) {
		int j, k, m, dim;
		for (     j = ce->graph[0], dim = cells_dim(ce)-1; dim >= 0; dim--)
			for (k = j-1; k >= 0; k--)
				if (dim == cells_dim(ce->ce[k])) cells_mutat(ce, k, --j);

////////////////////////////////////////////
//...mirror reflection of boundary elements;
		m = 0; while (cells_dim(ce->ce[m]) < cells_dim(ce)-1) m++;
		for (k = 0; k < (ce->graph[0]-m)/2; k++) cells_mutat(ce, m+k, ce->graph[0]-1-k);
	}
	return;
}

///////////////////////////////////////////////////////////////////////////////
//...ordering of elements by 0 dimension in common list of cell by 2 dimension;
void cells_ord0(Cells * ce, int id_ord)
{
	int  k, l, j = 0, arc, prev;
	if (id_ord) cells_ord(ce);
	if (cells_dim(ce) == 2 && ce->graph) {
		prev = ce->graph[ce->graph[1]+1];

///////////////////////////////////////
//...cycle on all elements of boundary;
		for (l = 0; l < ce->graph[1]; l++, prev = arc)
			if (ce->ce[arc = ce->graph[l+2]]->graph[1] == 2) {
			if (! topo_id(ce->ce[prev], k = ce->ce[arc]->graph[2]))
												 k = ce->ce[arc]->graph[3];
			cells_mutat(ce, k, j++);
		}
	}
	return;
}

/////////////////////////////////////////////////////
//...mirror reflection  of boundary elements of cell;
void bar_invers(Bar * bar)
{
	int dim, k, m;
	if (bar && ! bar->mp) {
		bar_ord(bar);

		dim = bar_dim(bar); m = bar->graph[0]-1;
		while (m > 0 && cells_dim(bar->ce[m]) < dim) m--;

		for (k = 0; k < (m+1)/2; k++) bar_mutat(bar, k, m-k);
	}
	return;
}

////////////////////////////////////////////////
//...change of one-dimensional cell orientation;
void cells_invers1(Cells * ce)
{
	if (ce && 1 == cells_dim(ce) && ce->graph) {
		int N = ce->graph[1], k;

////////////////////////////////////////////
//...mirror reflection of boundary elements;
		for (k = 0; k < N/2; k++) swap(ce->graph[k+2], ce->graph[ce->graph[1]+1-k]);

////////////////////////////////////////
//...transformation of geometrical card;
		if (ce->mp[0] == ID_MAP(1, SPHERE_GENUS) ||
			 ce->mp[0] == ID_MAP(1,   NULL_GENUS)) { //...transformation leaves X axis at place;
			 ce->mp[6]  = M_PI - ce->mp[6];
			 ce->mp[5] += M_PI;
		}
	}
	return;
}

////////////////////////////////////////////////
//...change orientation of two-dimensional cell;
void cells_invers2(Cells * ce)
{
	if (ce && 2 == cells_dim(ce) && ce->graph) {
		int N = ce->graph[1], k;

////////////////////////////////////////////////////////////////
//...mirror reflecion  and cyclic shifting of boundary elements;
		for (k = 0; k < N/2; k++) swap(ce->graph[k+2], ce->graph[N+1-k]);
		//for (k = 0; k < N-1; k++) swap(ce->graph[(k*(N-1))%N+2], ce->graph[((k+1)*(N-1))%N+2]); //...way shift of elements on one position???

/////////////////////////////////////////////
//...преобразуем геометрическую карту ячейки;
		if (ce->mp[0] == ID_MAP(2, NULL_GENUS)) { //...transformation leaves X axis at place;
			 ce->mp[6]  = M_PI - ce->mp[6];
			 ce->mp[5] += M_PI;
		}
		else
		if (ce->mp[0] == ID_MAP(2, SPHERE_GENUS) ||
			 ce->mp[0] == ID_MAP(2,    CYL_GENUS) ||
			 ce->mp[0] == ID_MAP(2,   CONE_GENUS) ||
			 ce->mp[0] == ID_MAP(2,  TORUS_GENUS)) { //...change sign of parameter, replyed to normal;
			 ce->mp[7]  = -ce->mp[7];
		}
	}
	return;
}

/////////////////////////////////////////
//...changing orientation of common cell;
void cells_invers(Cells * ce, int id_sub_element)
{
	if (ce && ce->mp && ce->graph) {
		switch(cells_dim(ce)) {
			case 1: cells_invers1(ce); break;
			case 2: cells_invers2(ce); break;
		}
		if (id_sub_element)
		for (int i = 0; i < ce->graph[1]; i++) cells_invers(ce->ce[ce->graph[i+2]]);
	}
	return;
}

//////////////////////////////////////////////////////////
//...cyclic shift of boundary elements of the common cell;
void cells_cyclic_shift(Cells * ce, int m)
{
	if (ce && ce->graph) {
		int   k, l, p, N = ce->graph[1];

///////////////////////////////////////
//...normalize modulus of cyclic shift;
		while ((m %= N) <= -1) m += N;
		while ( m       >=  N) m -= N;

////////////////////////////////////////////////
//...cyclic shift of the cell boundary elements;
		for (l = 0;                    m && l < abs(m); l++)
		for (k = 0, p = max(N/abs(m)-1, 1); k < p;      k++)
			swap(ce->graph[(k*(N-m)+l)%N+2], ce->graph[((k+1)*(N-m)+l)%N+2]);
	}
	return;
}

//////////////////////////////////////////////////////////////
//...auxilliary function -- removing of substructure elements;
void delete_cells_in_struct(Cells * ce)
{
	if (ce)
	for ( int k = 0; k < ce->graph[1]; k++) if (ce->graph[k+2] >= 0 && ce->ce[ce->graph[k+2]]) {
			int l = ce->pm ? ce->graph[1] : 0;
			delete_struct(ce->mp);

			for (l--; l >= 0; l--)
				delete_struct(ce->pm[l]);
				delete_struct(ce->pm);

			delete_cells_in_struct(ce->ce[ce->graph[k+2]]);

			delete_struct(ce->ce[ce->graph[k+2]]->graph);
			delete_cells (ce->ce[ce->graph[k+2]], NULL_STATE);
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//...identification substructure elements, adding new elements in list and removing repeting elements;
void search_cells_in_struct(Cells * ce, Cells * ext_ce, int & i, int id_cell, Cells **& dop_ce, int & N_buf, int & buf_size, int buf_delta)
{
	int k, m, l, m_ce;
	if (ce && ext_ce)
	for ( k = 0; k < ce->graph[1]; k++) if ((m_ce = ce->graph[k+2]) >= 0 && ce->ce[m_ce]) {
			m = id_cell != NULL_STATE  ? cells_in_struct(ext_ce, ce->ce[m_ce], OK_STATE) : 0;
			if (! m) { //...storage new number in surface structure;
				if (buf_size == N_buf) {
					Cells ** new_ce = new_struct<Cells *>(buf_size += buf_delta);
					memcpy  (new_ce, dop_ce, N_buf*sizeof(Cells *)); delete_struct(dop_ce); dop_ce = new_ce;
				}
				search_cells_in_struct(ce->ce[m_ce], ext_ce, i, id_cell, dop_ce, N_buf, buf_size, buf_delta);
				topo_correct(ce->ce[ce->graph[0]], m_ce, -(++i), NULL_STATE); dop_ce[N_buf++] = ce->ce[m_ce];
			}
			else if (ce->ce[m_ce]) { //...removing repeting elements of structure;
				l = ce->ce[m_ce]->pm ? ce->ce[m_ce]->graph[1] : 0;
				delete_struct(ce->ce[m_ce]->mp);

				for (l--; l >= 0; l--)
					delete_struct(ce->ce[m_ce]->pm[l]);
					delete_struct(ce->ce[m_ce]->pm);

				delete_cells_in_struct(ce->ce[m_ce]);
				topo_correct(ce->ce[ce->graph[0]], m_ce, -abs(m), NULL_STATE);

				delete_struct(ce->ce[m_ce]->graph);
				delete_cells (ce->ce[m_ce], NULL_STATE);
			}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//...including elemnts of structure in common list on the base of analysis of elements substructure;
int bar_add(Bar * bar, Cells *& ext_ce, int id_cell, int buf_delta)
{
////////////////////////////////////////////////////
//...initial checking and parameters initialization;
	if (! ext_ce)												  return(1);
	if (! bar || bar->mp || ! ext_ce || ! ext_ce->mp) return(0);

	int  i = bar->graph[0], k, l, m, N, N_buf = 0, buf_size = 0;
	Cells ** dop_ce = new_struct<Cells *>(buf_size += buf_delta);

///////////////////////////////////////////////////////////
//...identification of cell as a hole and links correction;
	m = id_cell != NULL_STATE ? cells_in_struct(bar, ext_ce, OK_STATE) : 0;

/////////////////////////////////////////////
//...storage new number in surface structure;
	if (! m) {
		dop_ce[N_buf++] = ext_ce; N = ++i;
		search_cells_in_struct(ext_ce, bar, i, id_cell, dop_ce, N_buf, buf_size, buf_delta);
		topo_correct(ext_ce);
	}

/////////////////////////////////////////////
//...removing repeting elements of structure;
	else {
		N = m; l = ext_ce->pm ? ext_ce->graph[1] : 0;
		delete_struct(ext_ce->mp);

		for (l--; l >= 0; l--)
			delete_struct(ext_ce->pm[l]);
			delete_struct(ext_ce->pm);

		delete_cells_in_struct(ext_ce);

		delete_struct(ext_ce->graph);
		delete_cells (ext_ce, NULL_STATE);
	}

//////////////////////////////////////
//...forming new common list of cells;
	if (i > bar->graph[0]) {
		Cells ** new_ce = new_struct<Cells *>(i+1);
		for (k = 0; k < bar->graph[0]; k++) {
			new_ce[k] =  bar->ce[k]; new_ce[k]->ce = new_ce;
		}
		for (m = 0; m < N_buf; m++) {
			new_ce[k] = dop_ce[m]; new_ce[k++]->ce = new_ce;
		}
		delete_struct(bar->ce); bar->ce = new_ce;
		bar->ce[bar->graph[0] = i] = bar;
	}
	install_struct(bar);

/////////////////////////////
//...return from the program;
	delete_struct(dop_ce);
	return abs(N);
}

////////////////////////////////////////////////////////////////////////////////////////////////
//...including of the cell insurface description (return number of included cell or error code);
int bar_add(Bar * bar, Cells *& ext_ce, int id_cell)
{
	if (! ext_ce)												  return(1);
	if (! bar || bar->mp || ! ext_ce || ! ext_ce->mp) return(0);

	Cells ** ce_ce = ext_ce->ce;
	int i, k, m, l, N, dim, dim_N = cells_dim(ext_ce);

////////////////////////////////////////////////////////////////////////////////////////////////////
//...identification of elements of geometrical strucutures and correction of link lists of the cell;
	for ( i = bar->graph[0], dim = dim_N; dim >= 0; dim--)
	for ( k = ext_ce->graph[0]; k >= 0; k--) if (dim == cells_dim(ce_ce[k])) {
			m = id_cell ? cells_in_bar(bar, ce_ce[k]) : 0;

////////////////////////////////////////////////////////
//...storing new numberof the cell in surface structure;
		if (k == ext_ce->graph[0]) N = m ? m : i+1;
		if (! m) topo_correct(ext_ce, k, -(++i));
		else {
			l = ce_ce[k]->pm && ce_ce[k]->graph ? ce_ce[k]->graph[1] : 0;
			delete_struct(ce_ce[k]->mp);

			for (l--; l >= 0; l--)
				delete_struct(ce_ce[k]->pm[l]);
				delete_struct(ce_ce[k]->pm);

			topo_correct(ext_ce, k, -abs(m));
		}
	}
	topo_correct(ext_ce);

/////////////////////////////////////////////////
//...forming new common list and install surface;
	if (i > bar->graph[0]) {
		Cells ** new_ce = new_struct<Cells *>(i+1);
      for (k = 0; k < bar->graph[0]; k++) {
			new_ce[k] = bar->ce[k]; bar->ce[k]->ce = new_ce;
      }
		for (dim = dim_N;  dim >= 0; dim--)
		for (m = ext_ce->graph[0]; m >= 0; m--) if (dim == cells_dim(ce_ce[m])) {
			new_ce[k] = ce_ce[m]; new_ce[k++]->ce = new_ce;
		}
		delete_struct(bar->ce); bar->ce = new_ce;
		bar->ce[bar->graph[0] = i] = bar;
	}
	install_struct(bar);

/////////////////////////////////////////////////////////
//...removing the used elementsand return of the program;
	for (k = 0; k <= ext_ce->graph[0]; k++) if (! ce_ce[k]->mp) {
		delete_struct(ce_ce[k]->graph); delete_struct(ce_ce[k]);
	}
	delete_struct(ce_ce); ext_ce = NULL;

	return abs(N);
}

//////////////////////////////////////////////////////////////
//...function for addition parametric description of the body;
void trim_add(Cells * ce, Cells *& trim, int id_mp)
{
	if (ce && trim && ce->mp && trim->mp && ! ce->pm && ce->graph && trim->graph && (ce->graph[1] <= trim->graph[1])) {
		int k, N_pad;

//////////////////////
//...justify topology;
		if  (topo_pad(ce->graph, N_pad = trim->graph[1]-ce->graph[1]))
		for (k = 0; k < N_pad; k++) ce->graph[k+2] = ce->graph[N_pad+2]/*ERR_GRAPH*/;

///////////////////////////////////////////////
//...distribute and transfer parametric curves;
		for (ce->pm = new_struct<CMap *>(ce->graph[1]), k = 0; ce->pm && k < ce->graph[1]; k++)
			swap(ce->pm[k], trim->ce[trim->graph[k+2]]->mp);

////////////////////////////////////////////////////////
//...transfer parametrization and release rest elements;
		if (id_mp) {
			delete_struct(ce->mp); swap(ce->mp, trim->mp);
		}
		delete_cells(trim);
	}
}

///////////////////////////////////////////////////////////
//...spanning of the geometrical card on surface structure;
int bar_span(Bar * bar, CMap *& ext_mp)
{
	if (! bar || ! ext_mp || bar->mp) return(1);

	int  k, j = 1, m = map_dim(ext_mp);
	if ((k = bar_dim(bar)) != m-1 && k != ERR_DIM) j = 0;

///////////////////////////////////////////////////////////////////////////////
//...search all elements by dimention minus one and including its in link list;
	for (k = 0; j && k < bar->graph[0]; k++)
		if (cells_dim(bar->ce[k]) == m-1 && ! topo_insert(bar->graph, k)) j = 0;

////////////////////////////////////////////////////////////////////////////
//...connecting geometrical card, ordering elements and return from program;
	bar->mp = ext_mp; ext_mp = NULL; cells_ord(bar);
	return(j);
}

/////////////////////////////////////////////////
//...excluding cell from the surface description;
void bar_exc(Bar * bar, int N)
{
	if (! bar ||  --N > bar->graph[0] || N < 0 || bar->mp) return;
	int i, j, k, l;

/////////////////////////////////////////////////////////////
//...look over common list and removing not needded elements;
	for (k = 0; k < bar->graph[0]; k++) if (topo_id(bar->ce[N], k)) {
		for (j = i = 0; ! j && i < bar->graph[0]; i++)
		if (i != N) j = topo_id(bar->ce[i], k, N);

////////////////////////////////////////////////////////////////
//...remove elemnt, which out of all cells and not a bar->ce[N];
		if (! j && k != N) {
			l = bar->ce[k]->pm && bar->ce[k]->graph ? bar->ce[k]->graph[1] : 0;
			delete_struct(bar->ce[k]->mp);

			for (l--; l >= 0; l--)
				delete_struct(bar->ce[k]->pm[l]);
				delete_struct(bar->ce[k]->pm);

			delete_struct(bar->ce[k]->graph);
			delete_struct(bar->ce[k]);

///////////////////////////////////////////////////////////
//...change numeration and remove element from common list;
			for (i = k; i < bar->graph[0];  i++) bar->ce[i] = bar->ce[i+1]; bar->graph[0] -= 1;
			for (i = 0; i < bar->graph[0];  i++) topo_exc(bar->ce[i]->graph, k);
			for (i = k; i < bar->graph[0]+1; i++) topo_correct((Cells *)bar, i+1, i);
		}
	}

///////////////////////////////
//...remove element bar->ce[N];
	l = bar->ce[N]->pm && bar->ce[N]->graph ? bar->ce[N]->graph[1] : 0;
	delete_struct(bar->ce[N]->mp);

	for (l--; l >= 0; l--)
		delete_struct(bar->ce[N]->pm[l]);
		delete_struct(bar->ce[N]->pm);

	delete_struct(bar->ce[N]->graph);
	delete_struct(bar->ce[N]);

/////////////////////////////////////////////////////////////////////////
//...change numbering and remove element bar->ce[N] from the common list;
	for (i = N; i < bar->graph[0];  i++) bar->ce[i] = bar->ce[i+1]; bar->graph[0] -= 1;
	for (i = 0; i < bar->graph[0];  i++) topo_exc(bar->ce[i]->graph, N);
	for (i = N; i < bar->graph[0]+1; i++) topo_correct((Cells *)bar, i+1, i);

///////////////////////////////////////////////
//...install structure and return from program;
	install_struct(bar);
	return;
}

/////////////////////////////////////////////
//...copy geometrical card (working program);
CMap * map_cpy(CMap * mp, int id_dop)
{
	if (mp) {
		int      N_mp = size_of_map(map_dim(mp), map_genus(mp)), k = N_mp;
		CMap * new_map = new_struct<CMap>(N_mp += (id_dop ? size_of_dop((int)mp[N_mp], mp) : 0)+1);
		if (! new_map) return(new_map);

///////////////////////////////////////////////////////////////////////////////
//...prepare new geometrical card of the selected cell and return from program;
		memcpy(new_map, mp, N_mp*sizeof(CMap));
		if  (! id_dop) new_map[k] = (CMap)NULL_CELL;

		return(new_map);
	}
	return(NULL);
}

///////////////////////////////////////////////////////////////////////////////////
//...copy description of the selected cell in parametric space (working procedure);
CMap ** pm_cpy(CMap ** pm, int l)
{
	if (pm) {
		CMap ** new_pm = new_struct<CMap *>(l);
		for (l--; new_pm && l >= 0; l--) new_pm[l] = map_cpy(pm[l]);
		return(new_pm);
	}
	return(NULL);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//...extractiong geometrical card of the selected cell from the surface description (working program);
CMap * map_cpy(Bar * bar, int N)
{
	if (bar && 0 <= N && N < bar->graph[0]) return(map_cpy(bar->ce[N]->mp));
	else                                    return(NULL);
}

/////////////////////////////////////////////////////////////////////////////////////////
//...extractiong of the selected cell description in parametric space (workin procedure);
CMap ** pm_cpy(Bar * bar, int N)
{
	if (bar && 0 <= N && N < bar->graph[0] && bar->ce[N] && bar->ce[N]->graph && bar->ce[N]->pm)
		return(pm_cpy(bar->ce[N]->pm, bar->ce[N]->graph[1])); else
		return(NULL);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
//...extracting topology graph of the selected cell from the surface description (working procedure);
Topo * graph_cpy(Bar * bar, int N)
{
	if (! bar || N >= bar->graph[0] || N < 0) return(NULL);
	Topo * new_graph = new_struct<Topo>(bar->ce[N]->graph[1]+2);
	if ( ! new_graph) return(new_graph);
	memcpy(new_graph, bar->ce[N]->graph, (bar->ce[N]->graph[1]+2)*sizeof(Topo));
	return(new_graph);
}

////////////////////////////////////////////////////////
//...copy cell from the surface description (by number);
Cells * bar_cpy(Bar * bar, int N, int id_origine)
{
	if (! bar || --N >= bar->graph[0] || N < 0) return(NULL);

	Cells * bar_N = bar->ce[N];
	Topo  * id;
	int N_dim = cells_dim(bar_N), i, k, m = 1;

/////////////////////////////////////////////
//...forming workpiece of the cell structure;
	Cells *	 ce = new_cells(); if (! ce) return(ce);
	init(ce, bar_N->graph[0]+1, 0, 0);

	id = new_struct<Topo>(bar_N->graph[0]+1);
	if (ce->ce && id) {
		ce->graph = graph_cpy(bar, N);
		ce->mp    = map_cpy  (bar, N);
		ce->pm    = pm_cpy   (bar, N);
		if (! ce->graph || ! ce->mp) m = 0;

/////////////////////////////////////////////////////////////////////////////////////
//...forming all lists, link list and geometrical card in structure of selected cell;
		for (i = k = 0; m && k < bar->graph[0]; k++)
		if ( i < ce->graph[0] && cells_dim(bar->ce[k]) < N_dim && topo_id(bar_N, k)) {
			if (! (ce->ce[i] = new_cells())) m = 0;
			else {
				ce->ce[i]->graph = graph_cpy(bar, k);
				ce->ce[i]->mp    = map_cpy  (bar, k);
				ce->ce[i]->pm    = pm_cpy   (bar, k);
				ce->ce[i]->ce    = ce->ce;
				if (! ce->ce[i]->graph || ! ce->ce[i]->mp) m = 0;
			}
			id[i++] = k;
		}
		if (! m || i != ce->graph[0]) goto err;

//////////////////////////////////////////////////////////////////////////////
//...reordering of the link lists and transforms cell to the initial position;
		for (i = 0; i < ce->graph[0]; i++) topo_correct(ce, id[i], -(i+1));
		topo_correct(ce);
		cells_ord	(ce);

		if (id_origine) cells_to_origine(ce);

/////////////////////////////
//...return from the program;
		delete_struct(id); return(ce);
	}
err:
	delete_struct(id); delete_cells(ce); return(NULL);
}

//////////////////////////////////////////////////////////////
//...extraction cell from the surface description (by number);
Cells * bar_sub(Bar * bar, int N, int id_origine)
{
	Cells * ce = bar_cpy(bar, N, id_origine);
	if (ce) bar_exc(bar, N);
	return(ce);
}

////////////////////////////////////////////////////////////////////
//...install dimension of the common structure for the all elements;
void install_struct(Bar * bar)
{
	if (bar && ! bar->mp)
	for ( int i = 0; i < bar->graph[0]; i++)
		if (bar->ce[i]) bar->ce[i]->graph[0] = bar->graph[0];
}

//////////////////////////////////////////////////////////////////////
//...removing degenerate elements from the surface (by one dimension);
void bar_generate(Bar * bar, double eps)
{
	for (int m, l, k = bar->graph[0]-1; bar && k >= 0; k--)
	if (1 == cells_dim(bar->ce[k]) &&
		(1 == bar->ce[k]->graph[1] ||
		 2 == bar->ce[k]->graph[1] && bar->ce[k]->graph[2] == bar->ce[k]->graph[3]))
	if (bar->ce[k]->mp[0] == ID_MAP(1, SPHERE_GENUS)) {
		if (fabs(bar->ce[k]->mp[7]) < eps ||
			 fabs(bar->ce[k]->mp[8]) < eps) bar_exc(bar, k+1);
	}
	else
	if (bar->ce[k]->mp[0] == ID_MAP(1, NULL_GENUS)) {
		 bar_exc(bar, k+1);
	}
/////////////////////////////////////////////////////////////////////////
//...removing from the surface two dimwnsionalelements (remov repetings);
	else
	if (2 == cells_dim(bar->ce[k]))
	for (l = bar->ce[k]->graph[1]-1; l > 0; l--) if (bar->ce[k]->graph[2+l] == bar->ce[k]->graph[1+l]) {
		for (m = l; m < bar->ce[k]->graph[1]; m++) {
			if (bar->ce[k]->pm)
			bar->ce[k]->pm[m-1]    = bar->ce[k]->pm[m];
			bar->ce[k]->graph[1+m] = bar->ce[k]->graph[2+m];
		}
		bar->ce[k]->graph[1]--;
	}
	return;
}

////////////////////////////////////////////////////////////////////////////////////
//...excluded double nbering of the point in the description of the surface patches;
void bar_correct_double_point(Bar * bar, double eps)
{
	if (bar)
	for (int k = bar->graph[0]-1, l, m; k >= 0; k--)
	if (1 == cells_dim(bar->ce[k]) && (
		 2 == bar->ce[k]->graph[1]  && (l = bar->ce[k]->graph[2]) != (m = bar->ce[k]->graph[3])
											 && max_point(bar->ce[l]->mp+1, bar->ce[m]->mp+1) < eps)) {

		(bar->graph[0])--; topo_correct((Cells *)bar, m, l);
		(bar->graph[0])++; bar_exc     ((Cells *)bar, m+1);

	}
	return;
}

//////////////////////////////////////////////////////////////////////////////////////
//...excluded double numberingof the point in description of the one dimensional cell;
void cell_correct_double_point(Cells * ce, double eps)
{
	int l, m;
	if (1 == cells_dim(ce) && (
		 2 == ce->graph[1]  && (l = ce->graph[2]) != (m = ce->graph[3])
								  && max_point(ce->ce[l]->mp+1, ce->ce[m]->mp+1) < eps)) {

		for (topo_correct(ce, m, l); m < ce->graph[0]; m++) {
			Cells * temp_ce = ce->ce[m]; ce->ce[m] = ce->ce[m+1]; ce->ce[m+1] = temp_ce;
		}
		ce->graph[0]--;

		delete_struct(ce->ce[m]->mp);
		delete_struct(ce->ce[m]->graph);
		delete_struct(ce->ce[m]);

		install_struct(ce);
	}
}

/*====================================================================================*/
/*                  AUXILLIARY TOOLS FOR LINES AND CIRCLE ARCS                        */
/*====================================================================================*/
//////////////////////////////////////////////////////////////
//...extractiong of the circle arc centre (for plane element);
complex get_arc_center(Cells * ce)
{
	if (ce && ce->mp && ID_MAP(1, SPHERE_GENUS) == ce->mp[0]) return comp(ce->mp[1], ce->mp[2]);
	else return comp(0.);
}

////////////////////////////////////////////////////////////////////////////////////////
//...extracting of the beginning and ending point of boundary chain (for plane element);
complex get_arc_point(Cells * ce, int m)
{
	if (1 == cells_dim(ce)) {
		if (m >= ce->graph[1]) m = ce->graph[1]-1;
		if (m <  0)            m = 0;
		return comp(ce->ce[ce->graph[2+m]]->mp[1], ce->ce[ce->graph[2+m]]->mp[2]);
	}
	return comp(0.);
}

////////////////////////////////////////////////////////////////////////////////
//...extracting of the beginning and ending point of boundary chain with return;
complex get_arc_point(Cells * ce, int m, double & CZ, double & SZ, double & CY, double & SY,
                                         double & CX, double & SX)
{
	if (1 == cells_dim(ce)) {
		if (m >= ce->graph[1]) m = ce->graph[1]-1;
      if (m <  0)            m = 0;
      double p[3] = {ce->ce[ce->graph[2+m]]->mp[1],
							ce->ce[ce->graph[2+m]]->mp[2],
							ce->ce[ce->graph[2+m]]->mp[3]};

      point_iso<double>(p, NULL, CZ, SZ, CY, SY, CX, SX);
      return comp(p[0], p[1]);
	}
	return comp(0.);
}

//////////////////////////////////////////////////////////////////////////////////////////
//...auxilliary function for in_bar method - increment of argument along the boundary arc;
double arc_delta(complex z0, complex z1, complex z2, int topo)
{
  double r, f0;
  if (z0 == z1) {
      if (abs(z2-z1) >= EE) {
          f0 = arg2(z2)-arg2(z1);
          if (f0 < -M_PI) f0 += 2.*M_PI;
          if (f0 >  M_PI) f0 -= 2.*M_PI;
          if (imag(conj(z2-z1)*(-z1)) < -EE_dop) f0 = -fabs(f0); else f0 = fabs(f0);
          return f0;
      }
  } else {
      if ((r = abs(z1-z0)) >= EE) {
          f0 = arg2(z2)-arg2(z1);
          if (f0 < -M_PI) f0 += 2.*M_PI;
          if (f0 >  M_PI) f0 -= 2.*M_PI;
          if (topo) {
             if (imag(conj(z2-z1)*(-z1)) < -EE_dop) f0 = -fabs(f0); else f0 = fabs(f0);
             if (f0 >= 0. && abs(z0)-r   < -EE_dop) f0 = f0-2.*M_PI;
          }
          if (! topo) {
             if (imag(conj(z2-z1)*(-z1)) < EE_dop) f0 = -fabs(f0); else f0 = fabs(f0);
             if (f0 <= 0. && abs(z0)-r   < EE_dop) f0 = f0+2.*M_PI;
          }
          return f0;
      }
  }
  return(0.);
}

//////////////////////////////////////////////////////////////////////
//...auxilliary function for the mesh - length arc (in complex fomat);
double arc_length(complex z0, complex z1, complex z2, int topo)
{
  if (z0 == z1) return abs(z2-z1);
  else {
     double f = (arg2(z2-z0)-arg2(z1-z0))*(1.-2.*topo);
     return abs(z2-z0)*(f+EE_dop < 0. ? 2.*M_PI+f : f);
  }
}

///////////////////////////////////////////////////////////////////////////////////////
//...installing of the beginning and ending point of the boundary chain (on the plane);
void set_arc_point(Cells * ce, int m, complex z)
{
	if (1 == cells_dim(ce)) {
		if (m >= ce->graph[1]) m = ce->graph[1]-1;
		if (m <  0)            m = 0;
		ce->ce[ce->graph[2+m]]->mp[1] = real(z);
		ce->ce[ce->graph[2+m]]->mp[2] = imag(z);
		ce->ce[ce->graph[2+m]]->mp[3] = 0.;
	}
	return;
}

////////////////////////////////////////////////////////
//...installing centre of the circle arc (on the plane);
void set_arc_center(Cells * ce, complex z)
{
	if (ce && ce->mp && ID_MAP(1, SPHERE_GENUS) == ce->mp[0]) {
		ce->mp[1] = real(z);
		ce->mp[2] = imag(z);
		ce->mp[3] = 0.;
	}
	return;
}

//////////////////////////////////////////////////////////////////////////////////
//...correction line description (installing local coordinate system and tangent);
void line_correct(Cells * ce, int id_fast)
{
	if (ce && ce->mp && ID_MAP(1, NULL_GENUS) == ce->mp[0] && ce->graph[1] > 0) {
		double X, Y, Z, fi,
			P1[3] = { X = ce->ce[ce->graph[2]]->mp[1],
						 Y = ce->ce[ce->graph[2]]->mp[2],
						 Z = ce->ce[ce->graph[2]]->mp[3]},
			P2[3] = { X, Y, Z };

		if (ce->graph[1] != 1) {
			P2[0] = ce->ce[ce->graph[3]]->mp[1];
			P2[1] = ce->ce[ce->graph[3]]->mp[2];
			P2[2] = ce->ce[ce->graph[3]]->mp[3];
		}
		X -= P2[0];
		Y -= P2[1];
		Z -= P2[2];
/////////////////////////////////
//...geometrical card correction;
		ce->mp[1] = (CMap)(P1[0]+P2[0])*.5;
		ce->mp[2] = (CMap)(P1[1]+P2[1])*.5;
		ce->mp[3] = (CMap)(P1[2]+P2[2])*.5;
		ce->mp[4] = (CMap)0.;
		ce->mp[5] = (CMap)0.;
		ce->mp[6] = (CMap)0.;
		if (id_fast == NULL_STATE) {
			ce->mp[4] = (CMap)(arg0(comp(-X, -Y)));
			ce->mp[5] = (CMap)(arg0(comp(-Z, sqrt(X*X+Y*Y))));
			ce->mp[6] = (CMap)(0.);
		}
		else {
			ce->mp[4] = 0.;
			ce->mp[5] = 0.;
			ce->mp[6] = 1.;
			if ((fi = sqrt(sqr(X)+sqr(Y)+sqr(Z))) > EE) {
				ce->mp[4] = (CMap)(-X*(fi = 1./fi));
				ce->mp[5] = (CMap)(-Y* fi);
				ce->mp[6] = (CMap)(-Z* fi);
			}
		}
	}
	return;
}

////////////////////////////////////////////////////////////////////////////////////////
//...correction circle description on the plane XY (installing local coordinate system);
void circ_correct(Cells * ce, double R, int topo, int points_inverse)
{
	if (ce && ce->mp && ID_MAP(1, SPHERE_GENUS) == ce->mp[0] && ce->graph[1] > 0) {
		if (! topo) {
			ce->mp[0]  = ID_MAP(1, NULL_GENUS); line_correct(ce);
		}
		else {
			double X, Y, Z, fi, theta, beta, sg, f0, f,
                P1[3] = { X = ce->ce[ce->graph[2]]->mp[1],
                          Y = ce->ce[ce->graph[2]]->mp[2],
                          Z = 0.  },
                P2[3] = { X, Y, Z };

			if (ce->graph[1] != 1) {
             P2[0] = ce->ce[ce->graph[3]]->mp[1];
             P2[1] = ce->ce[ce->graph[3]]->mp[2];
             P2[2] = Z;
         }
			if (points_inverse) {
				swap(P1[0], P2[0]); X = P1[0];
				swap(P1[1], P2[1]); Y = P1[1];
				swap(P1[2], P2[2]); Z = P1[2];
			}
         X -= P2[0];
         Y -= P2[1];
         Z -= P2[2];
         R  = (f = sqrt(X*X+Y*Y+Z*Z)) < EE_ker ? 0. : fabs(R);
         f0 = R == 0. ? 0. : R*R-f*f*.25;

         if ((f0 = filtr_(f0, EE_ker)) < 0.) {
             ce->mp[0] = ID_MAP(1, NULL_GENUS); line_correct(ce);
             return;
         }
         beta  = arg0(comp(sqrt(f0), f*.5)),
         fi    = arg0(comp(-X, -Y))-M_PI_2,
         sg    = topo > 0 ? -1. : 1.,
         theta = topo > 0 ?  0. : M_PI; if (2 == abs(topo)) beta = M_PI-beta;

         ce->mp[1] = (CMap)(P1[0]+sg*R*cos(fi+sg*beta));
         ce->mp[2] = (CMap)(P1[1]+sg*R*sin(fi+sg*beta));
         ce->mp[5] = (CMap)theta;
         ce->mp[6] = (CMap)(cos(theta)*fi); //...устанавливаем "доворот" в ce->mp[6];
         ce->mp[7] = (CMap)R;
         ce->mp[8] = (CMap)beta;
         ce->mp[3] =
         ce->mp[4] = (CMap)(fi = 0.);
		}
	}
	return;
}
