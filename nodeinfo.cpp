#include<iostream>
#include<assert.h>
#include "nodeinfo.h"
#include "objSegmentUtil.h"
using std::cout;
stNodePos* CreateNodeArray( unsigned long numItems)
{
     stNodePos *ptr;
     unsigned long  i;
//     if (numItems < 0)
//         return (0);
     ptr = new stNodePos[numItems];
     if(ptr == NULL)
     {
          cout << "cannot allocate space for node Array\n";
     }
     
     for(i=0;i<numItems;i++)
     {
          ptr[i].numOfAdjPos=0;
          ptr[i].adjPosList=NULL;
     }
     return ptr;
}

stNodePos* GetCurNode( stNodePos nodeList[], unsigned long index )
{
//  if ( index < 0 )
//      return (0);

  return ( &nodeList[index]);
}

stNodeData* CreateDataArray(unsigned long numItems)
{
//     if (numItems < 0)  return(0);
     stNodeData *ptr;
     ptr = new stNodeData[numItems];
     assert(ptr);
     return ptr;
}

stNodeData* GetCurNodeData ( stNodeData dataList[], unsigned long index )
{
     return ( &dataList[index]);
}



int AddCellToNodeList( stNodePos nodeList[], unsigned long index, unsigned long cellId)
{
     stNodePos    *curNodePtr;
     stCellIndex  *indexPtr, *cellListPtr;
     if((curNodePtr = &nodeList[index])==NULL)
          cout << "in addCellToPosList, posPtr is NULL\n";
     cellListPtr = curNodePtr->list;
     indexPtr = new stCellIndex;
     if(indexPtr==NULL)
     {
          cout << "in addCellToPosList, new error\n";
          return 0;
     }
     indexPtr->index = cellId;
     indexPtr->next = cellListPtr;
     curNodePtr->list = indexPtr;
     return METHOD_SUCCESS;
 }

stCellIndex* GetIncidentCells( stNodePos nodeList[], unsigned long index )
{
     stNodePos * ptr;
     ptr = &(nodeList[index]);
     return (ptr->list);
}
