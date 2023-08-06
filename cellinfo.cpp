// ************************************************************************* //
//  File: cellinfo.cpp
// ************************************************************************* //
#include "cellinfo.h"
#include "nodeinfo.h"
// What HANDLE macro does:
// curNodePtr is the pointer to the a th element of the Nnodelist which is the in the form of stNodePos 
// thus it gives the j th node's stNodePos info adress
// then check if that node's flag is signed as above the threshold(s)
//if so call the AddCellToNodeList method with the Nodelist index (a), the cell index (j) and the nodelist (which is the full stNodepos array in the Inptr )
// In the AddCellToNodeList method:
// add the current CellID to the inptr->pnode->list
#define HANDLE(a) { curNodePtr = GetCurNode( nodeList, a); \
                    if ( (curNodePtr->flag & F_THR) == F_THR ) {  \
                	    if(!AddCellToNodeList( nodeList, a, j))  \
                    	     return METHOD_FAILURE;} }
void SetCellNodes(unsigned long pos[8],stCellHex *cell)
{
   pos[0]=  cell->point0;
   pos[1]=  cell->point1; 
   pos[2]=  cell->point2; 
   pos[3]=  cell->point3;
   pos[4]=  cell->point4; 
   pos[5]=  cell->point5;
   pos[6]=  cell->point6; 
   pos[7]=  cell->point7;
}

int ProcessCells(InpObject<stCellHex> *inPtr)
{
  register unsigned long j = 0;
  
  stNodePos    *nodeList, *curNodePtr;
  stCellHex    *cellList, *curCellPtr;
  
  nodeList  = inPtr->pnode;
  cellList = inPtr->pcell; // cellist is a type of stCellHex (a struct)

  for (j=0;j<inPtr->numCell;j++)
  {  
       curCellPtr = GetCurCell<stCellHex>(cellList, j); // return the address of the j th cell = (in the form of stCellHex)
       HANDLE(curCellPtr->point0);
       HANDLE(curCellPtr->point1);
       HANDLE(curCellPtr->point2);
       HANDLE(curCellPtr->point3);
       HANDLE(curCellPtr->point4);    	
       HANDLE(curCellPtr->point5);
       HANDLE(curCellPtr->point6);    	
       HANDLE(curCellPtr->point7);
  }
  return METHOD_SUCCESS;
}
/*I am not sure if the template functions in the cellinfo.h file should be in .h file or here*/

