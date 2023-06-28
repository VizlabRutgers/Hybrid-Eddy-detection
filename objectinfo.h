#ifndef OBJECTINFO_H
#define OBJECTINFO_H

// ************************************************************************* //
//  File: objectinfo.h
// ************************************************************************* //
#include "nodeinfo.h"
#include "cellinfo.h"
#include "input.h"


struct stcellIndex;
/*struct stcellTetra;
struct stcellPrism;*/
struct stCellHex;


int AddCellToObj(stObject* list, stCellIndex *cell, char *mergeArray);
int AddCellListToObj(stObject *objPtr, stCellIndex* cellList,char *mergeArray);
stObject* CreateObject(void);






/*int outputAttrib(objListObj *obj, char *outFile);
int outputTrak( objListObj *obj, char *fileBaseName,int timeFrame);

void HANDLE_CELL_IND(inpObject<objTetra> *inPtr, objTetra *cell, long* vertex_ind,long* speedingtable);

void HANDLE_CELL_IND(inpObject<objPrism> *inPtr, objPrism *cell, long* vertex_ind,long* speedingtable);

void HANDLE_CELL_IND(inpObject<objHex> *inPtr, objHex *cell, long* vertex_ind,long* speedingtable);

void setoutfield(OMobj_id out, int nnodes, int ncells, int nspace, \
		 int celltype, float* coords, float* nodedata, int *connects);

*/


#endif
