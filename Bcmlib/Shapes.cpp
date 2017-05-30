#include "stdafx.h"
#include "shapes.h"

///////////////////////////
//...shapes identification;
Num_Shape shape_method(CBase * shape) 
{
	if (! shape) return NULL_SHAPE;
	if (typeid(* shape) == typeid(CMapi2DPoly<double>)  ||
		 typeid(* shape) == typeid(CMapi2DPoly<complex>) || 
		 typeid(* shape) == typeid(CMapi2DPoly<dd_real>) ||
		 typeid(* shape) == typeid(CMapi2DPoly<qd_real>)) return MP2D_POLY_SHAPE; else
	if (typeid(* shape) == typeid(CMapi3DPoly<double>)  ||
		 typeid(* shape) == typeid(CMapi3DPoly<complex>) || 
		 typeid(* shape) == typeid(CMapi3DPoly<dd_real>) ||
		 typeid(* shape) == typeid(CMapi3DPoly<qd_real>)) return MP3D_POLY_SHAPE; else
	if (typeid(* shape) == typeid(CMapi2DZoom<double>)  ||
		 typeid(* shape) == typeid(CMapi2DZoom<complex>) || 
		 typeid(* shape) == typeid(CMapi2DZoom<dd_real>) ||
		 typeid(* shape) == typeid(CMapi2DZoom<qd_real>)) return MP2D_ZOOM_SHAPE; else
	if (typeid(* shape) == typeid(CMapi3DZoom<double>)  ||
		 typeid(* shape) == typeid(CMapi3DZoom<complex>) || 
		 typeid(* shape) == typeid(CMapi3DZoom<dd_real>) ||
		 typeid(* shape) == typeid(CMapi3DZoom<qd_real>)) return MP3D_ZOOM_SHAPE; else
	if (typeid(* shape) == typeid(CMapi3DEll<double>)  ||
		 typeid(* shape) == typeid(CMapi3DEll<complex>) || 
		 typeid(* shape) == typeid(CMapi3DEll<dd_real>) ||
		 typeid(* shape) == typeid(CMapi3DEll<qd_real>)) return MP3D_ELLI_SHAPE; else
	if (typeid(* shape) == typeid(CBeamShape<double>)  ||
		 typeid(* shape) == typeid(CBeamShape<complex>) || 
		 typeid(* shape) == typeid(CBeamShape<dd_real>) ||
		 typeid(* shape) == typeid(CBeamShape<qd_real>)) return BEAM_POLY_SHAPE; else
	if (typeid(* shape) == typeid(CMapi2DCorner<double>)  ||
		 typeid(* shape) == typeid(CMapi2DCorner<complex>) || 
		 typeid(* shape) == typeid(CMapi2DCorner<dd_real>) ||
		 typeid(* shape) == typeid(CMapi2DCorner<qd_real>)) return MP2D_CORNER_SHAPE; else
	if (typeid(* shape) == typeid(CMapi2DClayer<double>)  ||
		 typeid(* shape) == typeid(CMapi2DClayer<complex>) || 
		 typeid(* shape) == typeid(CMapi2DClayer<dd_real>) ||
		 typeid(* shape) == typeid(CMapi2DClayer<qd_real>)) return MP2D_CLAYER_SHAPE; else
	if (typeid(* shape) == typeid(CMapi2DCircle<double>)  ||
		 typeid(* shape) == typeid(CMapi2DCircle<complex>) || 
		 typeid(* shape) == typeid(CMapi2DCircle<dd_real>) ||
		 typeid(* shape) == typeid(CMapi2DCircle<qd_real>)) return MP2D_CIRCLE_SHAPE; else
	if (typeid(* shape) == typeid(CMapi3DClayer<double>)  ||
		 typeid(* shape) == typeid(CMapi3DClayer<complex>) || 
		 typeid(* shape) == typeid(CMapi3DClayer<dd_real>) ||
		 typeid(* shape) == typeid(CMapi3DClayer<qd_real>)) return MP3D_CLAYER_SHAPE; else
	if (typeid(* shape) == typeid(CMapi3DSphere<double>)  ||
		 typeid(* shape) == typeid(CMapi3DSphere<complex>) || 
		 typeid(* shape) == typeid(CMapi3DSphere<dd_real>) ||
		 typeid(* shape) == typeid(CMapi3DSphere<qd_real>)) return MP3D_SPHERE_SHAPE; else
	if (typeid(* shape) == typeid(CMapi3DSpheroid))			 return MP3D_SPHEROID_SHAPE; else
	if (typeid(* shape) == typeid(CMapi3DSpheroidFull))	 return MP3D_SPHEROID_FULL_SHAPE; else
	if (typeid(* shape) == typeid(CHeat1DPoly<double>)  ||
		 typeid(* shape) == typeid(CHeat1DPoly<complex>) || 
		 typeid(* shape) == typeid(CHeat1DPoly<dd_real>) ||
		 typeid(* shape) == typeid(CHeat1DPoly<qd_real>))	 return HT1D_POLY_SHAPE; else
	if (typeid(* shape) == typeid(CHeat2DPoly<double>)  ||
		 typeid(* shape) == typeid(CHeat2DPoly<complex>) || 
		 typeid(* shape) == typeid(CHeat2DPoly<dd_real>) ||
		 typeid(* shape) == typeid(CHeat2DPoly<qd_real>))	 return HT2D_POLY_SHAPE; else
	if (typeid(* shape) == typeid(CHeat3DPoly<double>)  ||
		 typeid(* shape) == typeid(CHeat3DPoly<complex>) || 
		 typeid(* shape) == typeid(CHeat3DPoly<dd_real>) ||
		 typeid(* shape) == typeid(CHeat3DPoly<qd_real>))	 return HT3D_POLY_SHAPE; else
	if (typeid(* shape) == typeid(CHeat2DStreep<double>)  ||
		 typeid(* shape) == typeid(CHeat2DStreep<complex>) || 
		 typeid(* shape) == typeid(CHeat2DStreep<dd_real>) ||
		 typeid(* shape) == typeid(CHeat2DStreep<qd_real>)) return HT2D_STREEP_SHAPE; else
	if (typeid(* shape) == typeid(CHeat2DCircle<double>)  ||
		 typeid(* shape) == typeid(CHeat2DCircle<complex>) || 
		 typeid(* shape) == typeid(CHeat2DCircle<dd_real>) ||
		 typeid(* shape) == typeid(CHeat2DCircle<qd_real>)) return HT2D_CIRCLE_SHAPE; else
	if (typeid(* shape) == typeid(CHeat3DSphere<double>)  ||
		 typeid(* shape) == typeid(CHeat3DSphere<complex>) || 
		 typeid(* shape) == typeid(CHeat3DSphere<dd_real>) ||
		 typeid(* shape) == typeid(CHeat3DSphere<qd_real>)) return HT3D_SPHERE_SHAPE; else
	if (typeid(* shape) == typeid(CSkin2DBeamZ<double>)  ||
		 typeid(* shape) == typeid(CSkin2DBeamZ<complex>) || 
		 typeid(* shape) == typeid(CSkin2DBeamZ<dd_real>) ||
		 typeid(* shape) == typeid(CSkin2DBeamZ<qd_real>))  return SK2D_BEAMZSHAPE; else
	if (typeid(* shape) == typeid(CSkin3DBeam<double>)  ||
		 typeid(* shape) == typeid(CSkin3DBeam<complex>) || 
		 typeid(* shape) == typeid(CSkin3DBeam<dd_real>) ||
		 typeid(* shape) == typeid(CSkin3DBeam<qd_real>))	 return SK3D_BEAM_SHAPE; else
	if (typeid(* shape) == typeid(CSkin2DEll<double>)  ||
		 typeid(* shape) == typeid(CSkin2DEll<complex>) || 
		 typeid(* shape) == typeid(CSkin2DEll<dd_real>) ||
		 typeid(* shape) == typeid(CSkin2DEll<qd_real>))	 return SK2D_ELLI_SHAPE; else
	if (typeid(* shape) == typeid(CSkin2DPoly<double>)  ||
		 typeid(* shape) == typeid(CSkin2DPoly<complex>) || 
		 typeid(* shape) == typeid(CSkin2DPoly<dd_real>) ||
		 typeid(* shape) == typeid(CSkin2DPoly<qd_real>))	 return SK2D_POLY_SHAPE; else
	if (typeid(* shape) == typeid(CSkin3DPoly<double>)  ||
		 typeid(* shape) == typeid(CSkin3DPoly<complex>) || 
		 typeid(* shape) == typeid(CSkin3DPoly<dd_real>) ||
		 typeid(* shape) == typeid(CSkin3DPoly<qd_real>))	 return SK3D_POLY_SHAPE; else
	if (typeid(* shape) == typeid(CSkin3DZoom<double>)  ||
		 typeid(* shape) == typeid(CSkin3DZoom<complex>) || 
		 typeid(* shape) == typeid(CSkin3DZoom<dd_real>) ||
		 typeid(* shape) == typeid(CSkin3DZoom<qd_real>))	 return SK3D_ZOOM_SHAPE; else
	if (typeid(* shape) == typeid(CSkin3DExpp<double>)  ||
		 typeid(* shape) == typeid(CSkin3DExpp<complex>) || 
		 typeid(* shape) == typeid(CSkin3DExpp<dd_real>) ||
		 typeid(* shape) == typeid(CSkin3DExpp<qd_real>))	 return SK3D_EXPP_SHAPE; else
	if (typeid(* shape) == typeid(CWave2DPoly))				 return WV2D_POLY_SHAPE; else
	if (typeid(* shape) == typeid(CWave3DPoly))				 return WV3D_POLY_SHAPE; else
	if (typeid(* shape) == typeid(CAcou3DEll))				 return AU3D_ELLI_SHAPE; else
	if (typeid(* shape) == typeid(CAcou3DZoom))				 return AU3D_ZOOM_SHAPE; else
	if (typeid(* shape) == typeid(CAcou3DWave))				 return AU3D_WAVE_SHAPE; else
	if (typeid(* shape) == typeid(CAcou3DBeam))				 return AU3D_BEAM_SHAPE; else
	if (typeid(* shape) == typeid(CAcou3DBeamZ))				 return AU3D_BEAMZSHAPE; else
	if (typeid(* shape) == typeid(CAcou3DPoly))				 return AU3D_POLY_SHAPE; else
	if (typeid(* shape) == typeid(CSkin2DRadiiPoly<double>)  ||
		 typeid(* shape) == typeid(CSkin2DRadiiPoly<complex>) || 
		 typeid(* shape) == typeid(CSkin2DRadiiPoly<dd_real>) ||
		 typeid(* shape) == typeid(CSkin2DRadiiPoly<qd_real>)) return SK2D_RADII_POLY_SHAPE; else
	if (typeid(* shape) == typeid(CSkin2DRadiiZoom<double>)  ||
		 typeid(* shape) == typeid(CSkin2DRadiiZoom<complex>) || 
		 typeid(* shape) == typeid(CSkin2DRadiiZoom<dd_real>) ||
		 typeid(* shape) == typeid(CSkin2DRadiiZoom<qd_real>)) return SK2D_RADII_ZOOM_SHAPE; else return NULL_SHAPE;
}
