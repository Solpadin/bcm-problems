/*=========================================*/
/*                  DRAFTS                 */
/*=========================================*/
#ifndef ___DRAFTS___
#define ___DRAFTS___

#include <typeinfo>

#include "clame2d.h"
#include "cheat2d.h"
#include "cheat3d.h"
#include "ccohes2d.h"
#include "ccohes3d.h"
#include "cmindl2d.h"
#include "cmindl3d.h"
#include "chydro3d.h"
#include "cvisco2d.h"
#include "cvisco2d_grad.h"
#include "cporosity2d.h"

///////////////////////////////////////////////////////////////////////
//              ALL EXISTING TYPES OF DRAFT PROBLEMS                 //
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
//...global functions for construction of all existing types of drafts;
template <typename T> 
CDraft<T> * CreateDraft(Num_Draft id_DRAFT = NULL_DRAFT, int id_dop = 8)
{
	switch (id_DRAFT) {
      case			 NULL_DRAFT: return NULL;
      case		   BASIC_DRAFT: return (CDraft<T> *)(new CBase);
      case		  STRUCT_DRAFT: return					(new CDraft<T>);
      case		  LAME2D_DRAFT: if (typeid(T) == typeid(double))  return (CDraft<T> *)(new CLame2D(id_dop));
      case		  LAME3D_DRAFT: if (typeid(T) == typeid(double))  return (CDraft<T> *)(new CLame3D(id_dop));
      case		  HEAT2D_DRAFT: if (typeid(T) != typeid(complex)) return (CDraft<T> *)(new CHeat2D<T>(id_dop));
      case		  HEAT3D_DRAFT: if (typeid(T) != typeid(complex)) return (CDraft<T> *)(new CHeat3D<T>(id_dop));
      case		 HYDRO3D_DRAFT: if (typeid(T) != typeid(complex)) return (CDraft<T> *)(new CHydro3D<T>(id_dop));
      case		 COHES2D_DRAFT: if (typeid(T) == typeid(double))  return (CDraft<T> *)(new CCohes2D(id_dop));
      case		 COHES3D_DRAFT: if (typeid(T) == typeid(double))  return (CDraft<T> *)(new CCohes3D(id_dop));
      case		 MINDL2D_DRAFT: if (typeid(T) == typeid(double))  return (CDraft<T> *)(new CMindl2D(id_dop));
      case		 MINDL3D_DRAFT: if (typeid(T) == typeid(double))  return (CDraft<T> *)(new CMindl3D(id_dop));
      case		 VISCO2D_DRAFT: if (typeid(T) == typeid(complex)) return (CDraft<T> *)(new CVisco2D(id_dop));
      case VISCO2D_GRAD_DRAFT: if (typeid(T) == typeid(complex)) return (CDraft<T> *)(new CVisco2D_grad(id_dop));
      case	 POROSITY2D_DRAFT: if (typeid(T) == typeid(double))  return (CDraft<T> *)(new CPorosity2D(id_dop));
	}
   return NULL;
}
#endif
