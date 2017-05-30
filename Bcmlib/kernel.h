/*============================================*/
/*                  KERNEL                    */
/*============================================*/
#ifndef ___KERNEL___
#define ___KERNEL___

#include "utils.h"

#define  MAX_DIM               4
#define  ERR_DIM              -1
#define  ERR_GRAPH            -1
#define  MAX_SPLINE            9
#define  ID_MAP(map_dim, map_genus)  ((int)map_dim+(int)map_genus*MAX_DIM)

/////////////////////////////////////////////////
//...list of element codes in geometrical kernel;
//enum Num_Genus {
const int
     ERR_GENUS          = -1,
     NULL_GENUS         =  0,
     CYL_GENUS          =  1,
     ELL_CYL_GENUS      =  2,
     SPHERE_GENUS       =  3,
     SPHEROID_GENUS     =  4,
     ELLIPSOID_GENUS    =  5,
     CONE_GENUS         =  6,
     ELL_CONE_GENUS     =  7,
     TORUS_GENUS        =  8,
     PARAB_GENUS        =  9,
     ELL_PARAB_GENUS    = 10,
     BI_HYP_GENUS       = 11,
     ELL_BI_HYP_GENUS   = 12,
     HYPERB_GENUS       = 13,
     ELL_HYPERB_GENUS   = 14,
     NURBS_GENUS        = 15,
     NURBS_NULL_GENUS   = 16,
     NURBS_CYL_GENUS    = 17,
     BEZIER_GENUS       = 18,
     WEDGE_GENUS        = 19;
//};
typedef int    Topo;
typedef double CMap;

////////////////////////////////////////////////////
//...mail object of geometrical kernel - space cell;
struct Cells {
			Cells ** ce;    //...general list of all elementsof cell (self including);
			Topo  *  graph; //...list of connections in cell (topology);
			CMap  *  mp;    //...map of cell geometry;
			CMap  ** pm;    //...map of cell geometry in parametric space;
};
typedef  Cells Bar;		 //...body surface (cell with mp == NULL, number of surface elements = graph[0]);

///////////////////////////////////
//...identification of dimension cells (geometry with known parameters);
enum Num_Cell {
     NULL_CELL				= -1,
/////////////////////////
//...some standard cells;
     SHEET_CELL			= -2,
     BOX_CELL				= -3,
     UGOLOK_CELL			= -4,
     FACET_CELL			= -5,
     SPH_INTRUSION_CELL = -6,
     SHT_INTRUSION_CELL	= -7,
     CURVE1_TRIA_CELL	= -8,
     CURVE2_TRIA_CELL	= -9,
     CURVE1_QUAD_CELL	= -10,
     CURVE2_QUAD_CELL	= -11,
     CURVE11QUAD_CELL	= -12,
     CURVE3_QUAD_CELL	= -13,
     CURVE4_QUAD_CELL	= -14,
     B1SPLINE_CELL		= -15,
     B2SPLINE_CELL		= -16,
     B1POLY_CELL			= -17,
////////////////////////////////////////
//...nodes elements with boundary values;
     BOUNDARY_NODE		= -50,
////////////////////
//...special curves;
     ELLIPT_ARC			= -101,
/////////////////////////////////////////////////////////////////////
//...cells in isometric coordinates of orthogonal coordinates system;
     RING_SEGMENT			= -102,
     CYL_SEGMENT			= -103,
     SPH_SEGMENT			= -104,
     SPR_SEGMENT			= -105,
     CONE_SEGMENT			= -106,
     TORUS_SEGMENT		= -107,
///////////////
//...cap cells;
     SPH_BLEND				= -108,
     SPH_BLEND_DOP		= -109,
     SHEET_BLEND			= -110,
     SHEET_BLEND_DOP		= -111,
/////////////////////////////////////////////////
//...special cells for specifying block geometry;
     SPECIAL_SHEET		= -112,
     SPECIAL_BOX			= -113,
     BASIS_SHEET			= -114,
     BASIS_BOX				= -115,
};

///////////////////////////////////////////
//...extraction genus from geometrical map;
inline int map_genus(CMap * mp)
{
	return mp ? (int)mp[0] / MAX_DIM : ERR_GENUS;
}

////////////////////////////////////////////////
//...extraction dimension from geometrical card;
inline int map_dim(CMap * mp)
{
	return mp ? (int)mp[0] % MAX_DIM : ERR_DIM;
}

///////////////////////////////////////////////////////////////////////////////////////////
//...nable of sizes of geometrical card of space cell depending on the dimension and genus;
inline int size_of_map(int dim, int genus)
{
	switch (genus) {
		case         NULL_GENUS: return 0 == dim ? 4  : 7;  //...line;
		case          CYL_GENUS: return 1 == dim ? 9  : 10; //...ellipse;
		case       SPHERE_GENUS:                            //...circle;
		case         CONE_GENUS: return 1 == dim ? 9  : 8;  //...hyperbola;
		case     ELL_CONE_GENUS: return 1 == dim ? 8  : 10; //...parabola;
		case        PARAB_GENUS: return 8;
		case        TORUS_GENUS:
		case       BI_HYP_GENUS:
		case       HYPERB_GENUS:
		case      ELL_CYL_GENUS:
		case     SPHEROID_GENUS:
		case    ELL_PARAB_GENUS: return 9;
		case    ELLIPSOID_GENUS:
		case   ELL_BI_HYP_GENUS:
		case   ELL_HYPERB_GENUS:
		case   NURBS_NULL_GENUS:                            //...cylindrical surface;
		case    NURBS_CYL_GENUS: return 10;                 //...axisymmetric surface;
		case        NURBS_GENUS: return 1 == dim ? 10 : 12; //...B-curve or B-surface;
		case       BEZIER_GENUS: return 8;
		case        WEDGE_GENUS: return 9;
		default                : return 0;                  //...unknown geometrical card;
	}
}

inline int size_of_map(CMap * mp)
{
	return size_of_map(map_dim(mp), map_genus(mp));
}

////////////////////////////////////////////////////////////////////////////////
//...additional numer of parameters for dimension cell with parametric geometry;
inline int size_of_dop(int id_CELL, CMap * mp = NULL)
{
	switch (id_CELL) {
		case			  NULL_CELL : return 0;
		case			 SHEET_CELL : return 2;
		case SHT_INTRUSION_CELL :
		case				BOX_CELL : return 3;
		case			UGOLOK_CELL : return 7;
		case			 FACET_CELL :
		case	 CURVE1_TRIA_CELL :
		case	 CURVE2_TRIA_CELL : return 5; //...square and 4 parameters of boundary condition;
		case SPH_INTRUSION_CELL :
		case		 BOUNDARY_NODE : return 4; //...4 parameters of boundary condition;
		case		 B1SPLINE_CELL : return mp ? ((int)mp[7]/2)*(int)mp[9] +(int)mp[8] +2 : 0;
		case		 B2SPLINE_CELL : return mp ? ((int)mp[7]/4)*(int)mp[10]*(int)mp[11]+(int)mp[8]+(int)mp[9]+4 : 0;
		case			B1POLY_CELL : return mp ? ((int)mp[7]/2)*(int)mp[9] : 0;
		case			 ELLIPT_ARC : return 2;
		case			CYL_SEGMENT : return 2;
		case		  RING_SEGMENT :
		case			SPH_SEGMENT :
     case			SPR_SEGMENT :
		case		  CONE_SEGMENT :
		case		 TORUS_SEGMENT :
		case			SHEET_BLEND : return 3;
		case			  SPH_BLEND : return 4;
		case	  SHEET_BLEND_DOP : return 2;
		case		 SPH_BLEND_DOP : return 3;
		case			BASIS_SHEET :
		case		 SPECIAL_SHEET : return 4;
		case			  BASIS_BOX :
		case			SPECIAL_BOX : return 6;
		default					   : return 0;  //...unknown additional part of gemetrical card;
	}
}

inline int size_of_dop(CMap * mp)
{
	return size_of_dop((int)mp[size_of_map(mp)], mp);
}

/////////////////////////
//...genus of space cell;
inline int cells_genus(Cells * ce)
{
	return ce ? map_genus(ce->mp) : ERR_GENUS;
}

/////////////////////////////
//...dimansion of space cell;
inline int cells_dim(Cells * ce)
{
	return ce ? map_dim(ce->mp) : ERR_DIM;
}

///////////////////////////////////////////////////////
//...distribution of empty cell (with null parameters);
inline Cells * new_cells()
{
  return new_struct<Cells>();
}

////////////////////////////////////////////////////////
//...destruction of internal structure of abstract cell;
void release(Cells * ce);

///////////////////////////////////
//...distribution of abstarct cell;
inline void init (Cells * ce, int N_ce, int N_graph, int N_mp)
{
	if (ce) {
		release(ce);

		ce->ce    = new_struct<Cells *>(N_ce);
		ce->graph = new_struct<Topo>(N_graph);
		ce->mp    = new_struct<CMap>(N_mp);
		ce->pm    = NULL;
		if (ce->mp) ce->mp[N_mp-1] = (CMap)NULL_CELL;
		if (! ce->ce    && N_ce    > 0 ||
			 ! ce->graph && N_graph > 0 ||
			 ! ce->mp    && N_mp    > 0) release(ce);
		if (ce->ce)
			 ce->ce[N_ce-1] = ce;
	}
	return;
}

///////////////////////////////////////////
//...distribution of abstarct body surface;
inline Bar * get_bar(int N_ce)
{
	Cells *  bar; init(bar = new_cells(), N_ce+1, 2, 0);
	if (bar) bar->graph[0] = N_ce;
	return(Bar *)bar;
}

/////////////////////////////////////////////
//...full deleting abstarct cell and surface;
inline void delete_cells(Cells *& ce, int id_zero = OK_STATE)
{
	if (id_zero == OK_STATE) release(ce); delete_struct(ce);
	return;
}

inline void delete_bar (Bar *& bar)
{
	Cells * ce; delete_cells(ce = (Cells *)bar); bar = NULL;
	return;
}

//////////////////////////////////////
//...full deleting of compose surface;
void delete_bar(int N, Bar ** & bar);

///////////////////////////////
//...dimension of body surface;
int bar_dim(Bar * bar);

///////////////////////////////////////////////////////////////////////////
//...writing geometrical card and topology in ASCII code (working program);
CMap *  map_in(char * id_MAP,  unsigned long & count, unsigned long upper_limit);
Topo * topo_in(char * id_GRAF, unsigned long & count, unsigned long upper_limit, int N = 0);

void map_out (FILE * id_MAP,  CMap * mp, int id_long = NULL_STATE);
void topo_out(FILE * id_GRAF, Topo * graph);
void bar_out (char * ch_BAR,  Bar  * bar);

///////////////////////////////////////////////////
//...addition cell description in parametric space;
inline void add_pm(Cells * ce, int k, CMap *& mp)
{
	if (ce && ce->graph && 0 <= k && k < ce->graph[1] && ce->pm) {
		ce->pm[k] = mp; mp = NULL;
	}
	return;
}

/////////////////////////////////////////////////////////////////////////////////////////////
//...isometrical transformation of the point (two rotation about Z and one rotation about Y);
template <typename T>
inline void point_rotY(T * P, double C, double S)
{
	T		 X = P[0]*C+P[2]*S, Z = P[2]*C-P[0]*S;
	P[0] = X;  P[2] = Z;
}
template <typename T>
inline void point_rotZ(T * P, double C, double S)
{
   T		 X = P[0]*C-P[1]*S, Y = P[1]*C+P[0]*S;
   P[0] = X;  P[1] = Y;
}
template <typename T>
inline void point_shift(T * P, T * P_shift, int id_reverse = NULL_STATE)
{
	if (id_reverse) {
		P[0] -= P_shift[0];
		P[1] -= P_shift[1];
		P[2] -= P_shift[2];
	}
	else {
		P[0] += P_shift[0];
		P[1] += P_shift[1];
		P[2] += P_shift[2];
	}
}
template <typename T>
inline void point_iso(T * P, T * P_shift, double CZ, double SZ,
							 double CY, double SY,  double CX, double SX)
{
   point_rotZ(P, CX, SX); //...Euler angles in inverse order !!!
   point_rotY(P, CY, SY);
   point_rotZ(P, CZ, SZ);
   if (P && P_shift) point_shift(P, P_shift);
}

///////////////////////////////////////////////////////////////
//...relative displacement and rotation of the point in radian;
inline void point_iso(double * P, double * P0, double fi = 0., double theta = 0., double fX = 0.)
{
	if (P) {
		double CZ = cos(fi), SZ = sin(fi), CY = cos(theta), SY = sin(theta),
				 CX = cos(fX), SX = sin(fX);
		point_iso(P, P0, CZ, SZ, CY, SY, CX, SX);
	}
	return;
}

//////////////////////////////////////////////////////////////////////
//...isometrical transformation of geometrical card, cell and surface;
void  map_normalizat      (CMap * mp, double * P, double * ort);
void  map_normalizat_facet(CMap * mp, double * P, double * ort);
void  map_iso  (CMap  * mp, double * P, double & CZ, double & SZ, double & CY, double & SY, double & CX, double & SX);
void  map_iso  (CMap  * mp, double * P, double fi = 0., double theta = 0., double fX = 0.);
void  cells_iso(Cells * ce, double * P, double & CZ, double & SZ, double & CY, double & SY, double & CX, double & SX);
void  cells_iso(Cells * ce, double * P, double fi = 0., double theta = 0., double fX = 0.);
void  bar_iso  (Bar  * bar, double * P, double fi = 0., double theta = 0., double fX = 0.);
void  cells_to_origine(Cells * ce);

/////////////////////////////////////////////////////////////////////////////////
//...transformation of gemetrical elements in local and common coordinate system;
void  make_local  (double * P, CMap * mp, int id_shift = OK_STATE);
void  make_common (double * P, CMap * mp, int id_shift = OK_STATE);
void  map_local   (CMap   * mp,  CMap * ext_mp);
void  map_common  (CMap   * mp,  CMap * ext_mp);
void  cells_local (Cells * ce,  CMap * ext_mp);
void  cells_common(Cells * ce,  CMap * ext_mp);
void  bar_local   (Bar   * bar, CMap * ext_mp);
void  bar_common  (Bar   * bar, CMap * ext_mp);

//////////////////////////////////////////////////////////////////////////////////
//...arc length of one-dimensional cell and correction of local coordinate system;
inline int get_num(Topo * graph, int m)
{
	if (graph) {
		if (m >= graph[1]) m = graph[1]-1; else
		if (m <  0)        m = 0;
		if (m >= 0) return graph[2+m];
	}
	return(-1/*0 -- до 11.03.2010*/);
}
double cells_length(Cells * ce);
void   edge_correct(Cells * ce);

/////////////////////////////////////////////////////////
//...generation of geometrical card with null parameters;
inline CMap * get_map(int dim, int genus, int id_cell = NULL_CELL)
{
	int     k = size_of_map(dim, genus), id_dop = size_of_dop(id_cell);
	CMap * mp = k && id_dop >= 0 ? new_struct<CMap>(k+1+id_dop) : NULL;
	if (mp) {
		mp[0] = ID_MAP(dim, genus);
		mp[k] = (CMap)id_cell;
	}
	return mp;
}

inline int add_new_maps(CMap *& mp, int N_map, int N_new_maps)
{
	CMap * new_map = N_new_maps > 0 ? new_struct<CMap>(N_map+N_new_maps) : NULL;
	if (new_map) {
		if (mp) memcpy(new_map, mp, N_map*sizeof(CMap));
		delete_struct(mp); mp = new_map; return(1);
	}
	return(0);
}

inline int add_new_cell_maps(CMap *& mp, int id_CELL)
{
	int m = size_of_map(mp);
	CMap * new_map = new_struct<CMap>(m+1+size_of_dop(id_CELL, mp));
	if (new_map) {
		memcpy(new_map, mp, m*sizeof(CMap)); new_map[m] = id_CELL;
		delete_struct(mp); mp = new_map;
		return(1);
	}
	return(0);
}

/////////////////////////////////////////////////////////
//...identification of geometrical card/ cell and surface;
void map_to(CMap *& mp, int dim, int genus);
inline void cells_to(Cells * ce, int genus) { if (ce) map_to(ce->mp, map_dim(ce->mp), genus);}
int  map_id    (CMap  * mp, CMap   * ext_mp);
int  map_into  (CMap  * mp, double X, double Y, double Z, int id_special = OK_STATE);
int  map1into  (CMap  * mp, double X, double Y, double Z);
int  cells_id  (Cells * ce, Cells * ext_ce, int id_num_correct = NULL_STATE);
int  common_id (Cells * ce, Cells * ext_ce);
int  cells_in_struct(Cells * ce, Cells * ext_ce, int id_num_correct);
int  cells_in_bar	  (Bar   * bar, Cells * ce);
void search_cells_element(Cells * ce, int *& id, int & N_buf, int & buf_size, int buf_delta = 20);

///////////////////////////////////////////////
//...addition new element in graph of topology;
int   topo_insert	(Topo *& graph, Topo new_element);
void  topo_exc		(Topo  * graph, Topo some_element);
void	topo_correct(Cells * ce, int N = -1, int N_new = -1, int id_list = OK_STATE);

//////////////////////////////////////////////////////
//...identification of topology and elements ordering;
void  topo_ord (Topo  * graph);
int   topo_pad (Topo *& graph, int N_pad);
int   topo_id  (Topo  * graph, Topo some_element);
int   topo_id  (Cells * ce,    Topo some_element, Topo exc_element = -1);
void  bar_mutat  (Bar * bar,  int k, int m);
void  cells_mutat(Cells * ce, int k, int m);
void  bar_ord    (Bar   * bar);
void  cells_ord  (Cells * ce);
void  cells_ord0 (Cells * ce, int id_ord = 0);
void  bar_invers (Bar   * bar);

//////////////////////////////////////////////////////////////////////////////
//...element orientation and cyclic shift of boundary elements in common cell;
void cells_invers      (Cells * ce, int id_all = NULL_STATE);
void cells_cyclic_shift(Cells * ce, int      m = 1);

////////////////////////////////////////////////////////////////////////////////////
//...extractiong of geometrical card of specified cell from the surface description;
CMap * map_cpy	 (CMap  * mp, int id_dop = OK_STATE);
CMap ** pm_cpy	 (CMap ** pm, int l);
CMap * map_cpy	 (Bar * bar, int N);
CMap ** pm_cpy	 (Bar * bar, int N);
Topo * graph_cpy(Bar * bar, int N);

///////////////////////////////////////////////////
//...composing body surface from separate elements;
void delete_cells_in_struct(Cells * ce);
void search_cells_in_struct(Cells * ce, Cells * ext_ce, int & i, int id_cell, Cells **& dop_ce, int & N_buf, int & buf_size, int buf_delta = 20);
int     bar_add (Bar * bar, Cells *& ext_ce, int id_cell, int buf_delta);
int     bar_add (Bar * bar, Cells *& ext_ce, int id_cell = OK_STATE);
int     bar_span(Bar * bar, CMap  *& mp);
void    trim_add(Cells * ce, Cells *& trim, int id_mp = OK_STATE);
void    bar_exc (Bar * bar, int N);
Cells * bar_cpy (Bar * bar, int N, int id_origine = NULL_STATE);
Cells * bar_sub (Bar * bar, int N, int id_origine = NULL_STATE);

///////////////////////////////////////////////////////////////////////////////////////////////////
//...deleting degenerate elementsof one dimention and excluding of double numeration of the points;
void  install_struct				(Bar * bar);
void  bar_generate            (Bar * bar,  double eps = EE_ker);
void  bar_correct_double_point(Bar * bar,  double eps = EE_ker);
void cell_correct_double_point(Cells * ce, double eps = EE_ker);

/////////////////////////////////////////////
//...operation with points of plane elements;
complex get_arc_center(Cells * ce);
complex get_arc_point (Cells * ce, int m);
complex get_arc_point (Cells * ce, int m, double & CZ, double & SZ, double & CY, double & SY, double & CX, double & SX);
double  arc_delta		 (complex z0, complex z1, complex z2, int topo);
double  arc_length	 (complex z0, complex z1, complex z2, int topo);
void set_arc_point    (Cells * ce, int m, complex z);
void set_arc_center   (Cells * ce, complex z);

////////////////////////////////////////////////////////////////////////////////////////////
//...correction of the distribution of the line and circle (set of local coordinate system);
void line_correct(Cells * ce, int id_fast= NULL_STATE);
void circ_correct(Cells * ce, double R, int topo, int points_inverse = NULL_STATE);

////////////////////////////////////////////////////////
//...description of the common part of geometrical card;
//   mp[0]               - card identification;
//   mp[1], mp[2], mp[3] - coordinates of the reference point;
//   mp[4], mp[5], mp[6] - Euler angles of the local coordinate system;

//////////////////////////////////////////////////////////////////////////////
//...description of the additional part of parameters of the geometrical card;
//   ID_MAP(1,       NULL_GENUS): absent;
//
//   ID_MAP(2,       NULL_GENUS): absent;
//
//   ID_MAP(3,       NULL_GENUS): absent;
//
//   ID_MAP(1,      NURBS_GENUS): mp[7] - control points dimention + curve periodicity;
//                                mp[8] - nodes number;
//                                mp[9] - control points number;
//
//   ID_MAP(2,      NURBS_GENUS): mp[7]  - control points dimention + periodicity in U- and V-direction;
//                                mp[8]  - nodes number in U-direction;
//                                mp[9]  - nodes number in V-direction;
//                                mp[10] - control points number in U-direction;
//                                mp[11] - control points number in V-direction;
//
//   ID_MAP(2, NURBS_NULL_GENUS): mp[7]  - control points dimention + curve periodicity;
//                                mp[8]  - nodes number;
//                                mp[9]  - control points number;
//
//   ID_MAP(2,  NURBS_CYL_GENUS): mp[7]  - control points dimention + curve periodicity;
//                                mp[8]  - nodes number;
//                                mp[9]  - control points number;
//
//   ID_MAP(1,     SPHERE_GENUS): mp[7] - circle arc radius R >= 0;
//                                mp[8] - semilength of the circle arc in radian;
//
//   ID_MAP(2,     SPHERE_GENUS): mp[7] - sphere radius (R < 0 - internal normal);
//
//   ID_MAP(2,   SPHEROID_GENUS): mp[7] - spheroid major semiaxis;
//                                mp[8] - spheroid minor semiaxis (may be > major semiaxis);
//
//   ID_MAP(2,  ELLIPSOID_GENUS): mp[7] -
//                                mp[8] -
//                                mp[9] - spheroid semiaxises;
//
//   ID_MAP(1,        CYL_GENUS): mp[7] - ellipse major semiaxis;
//                                mp[8] - ellipse minor semiaxis;
//
//   ID_MAP(2,        CYL_GENUS): mp[7] - cylindrical radius (R < 0 - internal normal);
//                                mp[8] -
//                                mp[9] - beginning and ending point of cylinder axis (if exist);
//
//   ID_MAP(2,       CONE_GENUS): mp[7] - semiangle of cone in radian  (mp[7] < 0 - internal normal);
//
//   ID_MAP(2,      TORUS_GENUS): mp[7] - major radius R of toroidal coordnate (R < 0 - internal normal);
//                                mp[8] - minor radius r
//                                (parameter alpha > 0 of toroidal coordinates z+ir = i th((alpha+ibeta)/2 is absent);

/////////////////////////////////////////
//...description of geometrical elements;
//       NULL_GENUS : 0 - point; 1 - line;       2 - plane;      3 - space;
//     SPHERE_GENUS :            1 - circle;     2 - sphere;
//        CYL_GENUS :            1 - ellipse;    2 - cylinder;
//    ELL_CYL_GENUS :                            2 - elliptic cylinder;
//       CONE_GENUS :            1 - hyperbola;  2 - cone;
//   ELL_CONE_GENUS :            1 - parabola;   2 - elliptic cone;
//      TORUS_GENUS :                            2 - torus;
//   SPHEROID_GENUS :                            2 - spheroid;
//  ELLIPSOID_GENUS :                            2 - ellipsoid;
//      PARAB_GENUS :                            2 - paraboloid;
//  ELL_PARAB_GENUS :                            2 - elliptic paraboloid;
//     BI_HYP_GENUS :                            2 - two-sheeted hyperboloid;
// ELL_BI_HYP_GENUS :                            2 - elliptic two-sheeted hyperboloid;
//     HYPERB_GENUS :                            2 - hyperboloid;
// ELL_HYPERB_GENUS :                            2 - elliptic hyperboloid;
//      NURBS_GENUS :            1 - B-curve;    2 - B-surface;
//  NURBS_CYL_GENUS :                            2 - B-axisymmetrical surface;
// NURBS_NULL_GENUS :                            2 - B-cylindrical surface;
#endif
