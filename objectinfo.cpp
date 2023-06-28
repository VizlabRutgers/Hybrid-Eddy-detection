#include "objectinfo.h"
#include "objSegmentUtil.h"
stObject* CreateObject(void)
{
     stObject *ptr = new stObject;
     ptr->cellPtr = 0x00;
     ptr->next     = 0x00;
     ptr->endPtr   = 0x00;
     return (ptr);
}


int AddCellListToObj(stObject *objPtr, stCellIndex* cellList,char *mergeArray)
{
     int rc;
     stCellIndex *ptr= cellList, *indexPtr;
     
     while (ptr)
     {
          indexPtr = new stCellIndex;
          indexPtr->index=ptr->index;
          rc = AddCellToObj(objPtr, indexPtr,mergeArray);     
          if(rc)
               delete indexPtr; 
     // rc=1 means this objIndex has been added before addCellToObj(),
     // so the space allocated can be freed.
     
          ptr=ptr->next;      /* move to next tetra on tetralist */
     }

     return METHOD_SUCCESS;
}

int AddCellToObj(stObject* list, stCellIndex *cell, char *mergeArray)
{
     stCellIndex *objPtr;
     objPtr = list->cellPtr;
     if(!objPtr) 
     {
         list->cellPtr= cell;
         list->endPtr  = cell;
         cell->next = 0;
         mergeArray[cell->index] = 1; 
            // if mergeArray[i]=1, means the cell i has been added to an object.
            // if mergeArray[i]=0, means the cell 0 has not been added to any object.
         return 0;
     }
  
     if(mergeArray[cell->index] > 0) return 1; 
   // return 1, means objIndex cell has already been added before this function is called. 

     cell->next = NULL;
     list->endPtr->next = cell;
     list->endPtr = cell;

     mergeArray[cell->index] = 1;

     return 0;
}
