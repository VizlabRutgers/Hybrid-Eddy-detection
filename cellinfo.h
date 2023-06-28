#ifndef CELLINFO_H
#define CELLINFO_H
// ************************************************************************* //
//  File: cellinfo.h
// ************************************************************************* //

//#include "objectinfo.h"
#include "input.h"
//#include "objectinfo.h"
#include "nodeinfo.h"
#include "objectinfo.h"


extern int AddCellToObj(stObject* list, stCellIndex *cell, char *mergeArray);
extern void SetArray(unsigned long index,char *mergeArray);
extern void ResetArray(unsigned long nnodes, char *mergeArray);// ****************************************************************************
//  Struct: stcellTetra
//
//  Purpose: holds the information of a Tetragonal cell
//      
//  Programmer: Rohini Pangrikar 
//  Creation:   Wed May 28 2008
// ****************************************************************************//

struct stcellTetra 
{
  unsigned long point0;
  unsigned long point1;
  unsigned long point2;
  unsigned long point3;
  unsigned long objNum; /* this cell is a part Object with number 'objNum' in global object list */
  unsigned long flag;
};
// ****************************************************************************
//  Struct: stcellPrism
//
//  Purpose: holds the information of a Prism type cell
//      
//  Programmer: Rohini Pangrikar 
//  Creation:   Wed May 28 2008
// ****************************************************************************//
struct stcellPrism 
{
  unsigned long point0;
  unsigned long point1;
  unsigned long point2;
  unsigned long point3;
  unsigned long point4;
  unsigned long point5;
  unsigned long objNum;/* this cell is a part Object with number 'objNum' in global object list */
  unsigned long flag;
};

// ****************************************************************************
//  Struct: stCellHex
//
//  Purpose: holds the information of a Hexagonal cell
//      
//  Programmer: Rohini Pangrikar 
//  Creation:   Wed May 28 2008
// ****************extern int AddCellToObj(stObject* list, stCellIndex *cell, char *mergeArray);************************************************************//
struct stCellHex 
{
  unsigned long point0;
  unsigned long point1;
  unsigned long point2;
  unsigned long point3;
  unsigned long point4;
  unsigned long point5;
  unsigned long point6;
  unsigned long point7;
  unsigned long objNum;/* this cell is a part Object with number 'objNum' in global object list */
  unsigned long flag;
};

template <class objType>  // objType can be objPrism, objTetra, objHex
objType* CreateCellArray(unsigned long numItems)
{
     objType *ptr;
     ptr = new objType[numItems];
     if(ptr == NULL)
     {
          cout << "cannot allocate space for Cell Array\n";
     }
   
     return ptr;
}

template <class objType> // objType can be objPrism, objTetra, objHex
objType *GetCurCell( objType cellList[], unsigned long index )
{
  return( & cellList[index] );
}

template <class T>
int MarkCellList(T *cellList, stObject *objPtr)
{
     unsigned long objNum = objPtr->objNum ;
     stCellIndex *ptr;
     T *cell;
     ptr = objPtr->cellPtr;
     while (ptr != 0 )
    {
         cell=GetCurCell(cellList, ptr->index);
         cell->objNum=objNum;
         ptr = ptr->next;
    }
    return 1;
}

template <class UCDtype>
class InpObject;

template <class T>
int MergeCellList( InpObject<T> *inPtr, stObject *objPtr, long objNum, T *cellList, char *mergeArray)
{
     stObject  *ptr;
     stCellIndex    *indexPtr;
     ptr = inPtr->pobject;
     while(ptr)
     {
          if (ptr->objNum == objNum )
               break;
          ptr=ptr->next;
     }

     if (ptr == 0 )  return (-1);

/* merge the 2 lists */

     ResetArray(inPtr->numCell,mergeArray);
     indexPtr = ptr->cellPtr;
     while( indexPtr )
     {
          //SetArray(indexPtr->index,mergeArray);
          mergeArray[indexPtr->index] = 1;
          indexPtr = indexPtr->next;
     }
     indexPtr = objPtr->cellPtr;
     while(indexPtr != 0 )
     {
          AddCellToObj(ptr, indexPtr,mergeArray);
          indexPtr = indexPtr->next;
     }
  
     MarkCellList(cellList, ptr);
     return 1;
}


//from seg.h
//void SetCellNodes(long pos[8],stCellTetra *cell);
//void SetCellNodes(long pos[8],stCellPrism *cell);





// ************************************************************************* //
//  END: cellinfo.h
// ************************************************************************* //
#endif
