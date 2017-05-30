/*========================================*/
/*                  NODES                 */
/*========================================*/
#ifndef ___NODES___
#define ___NODES___

#include <typeinfo>

#include "cgrid_el.h"
#include "cgrid_qg.h"

//////////////////////////////////////////////////////////////////////////
//                  ALL EXISTING TYPES OF GRID NODES                    //
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//...global function for construction of all existing types of grid nodes;
CGrid * CreateNodes(Num_Nodes id_NODES)
{
	switch (id_NODES) {
     case GRID_IN_NODES: return new CGrid;
     case GRID_EL_NODES: return new CGrid_el;
     case GRID_QG_NODES: return new CGrid_QG;
	}
   return NULL;
}
#endif
