#ifndef INPUT_H
#define INPUT_H
// ************************************************************************* //
//  File: input.h
// ************************************************************************* //
#include <string>
#include "stobject1.h"
//#include "objectinfo.h"
#include <vtkRectilinearGrid.h>
#include <vtkDataSet.h>
#include "nodeinfo.h"
#include "cellinfo.h"
//#include "objectinfo.h"
#include"objSegmentUtil.h"

/*struct stNodePos;
struct stNodeData;
struct stObject;
struct stCellIndex;*/
//struct stObject;

// ****************************************************************************
//  Class: InpObject
//
//  Purpose: holds the information of the input data
//      
//  Programmer: Rohini Pangrikar 
//  Creation:   Wed May 28 2008
// ****************************************************************************//
template <class CELLtype>
class InpObject 
{ 
 public:
  static int celltype;	
  char fileName[256];
  unsigned long numObj;
  unsigned long numCell;
  unsigned long numNodes;
  int  colorMethod;
  float minData;
  float maxData;
  float minData1;
  float maxData1;
  float minData2;
  float maxData2;  
  float thrVal;
  float thrVal1;
  float thrVal2;  
  float thresh_deltax;
  float thresh_deltay;
  float thresh_deltaz;
  int   thrType;
  float threshPercent_FromTop;
  float threshPercent_FromBottom;
  int nncomp;
  long MaxPacketnumber; //TotalPacket number in the current frame
  CELLtype  *pcell; // simple free
  stNodePos    *pnode; // complex free
  stNodeData   *pnodeData; //simple free
  stObject *pobject; // complex free

  stPacket *packets;    // contains the packets 

  InpObject() {pcell=0; pnode=0; pnodeData=0; pobject=0; packets=0;}
  ~InpObject(); 
  int freeCells(stCellIndex*);
  int freePacketLists(ObjIndex*);
};

//*** +++


template<class T>
InpObject<T> :: ~InpObject() 
{
     if(pcell)
          delete [] pcell;
     if(pnodeData)
          delete [] pnodeData;  
     if(packets)
          delete [] packets;      

     stObject *objPtr, *nextPtr; // freeObjects
     stCellIndex   *cellPtr;
     objPtr=pobject;
       

    while(objPtr)
    {
         nextPtr=objPtr->next;
  	 cellPtr=objPtr->cellPtr;
  	 if(cellPtr)
              freeCells(cellPtr);
  	 delete objPtr;
	 objPtr=nextPtr;  	
    }
    stCellIndex *list, *nextlist; // freePos
    register int i;
    for(i=0;i<numNodes;i++)
    {
         delete pnode[i].adjPosList;
         list=pnode[i].list;
         while(list)
         {
              nextlist=list->next;
              delete list;
              list = nextlist;
         }
    }
    delete [] pnode;
/*/--------------------------
    ObjIndex *ObjList, *nextlist1; // freePos
    //register int i;
    for(i=0;i<MaxPacketnumber;i++)
    {
         ObjList=packets[i].ObjList;
         while(ObjList)
         {
              nextlist1=ObjList->next;
              delete ObjList;
              ObjList = nextlist1;
         }
    }
    delete [] packets;
*/ //--------------------------
/*     stPacket *objPtr1, *nextPtr1; // freeObjects
     ObjIndex   *cellPtr1;
     objPtr1=packets;
       

    while(objPtr1)
    {
         nextPtr1=objPtr1->next;
  	 cellPtr1=objPtr1->ObjList;
  	 if(cellPtr1)
              freePacketLists(cellPtr1);
  	 delete objPtr1;
	 objPtr1=nextPtr1;  	
    }
*/    
    
}

template<class T>
int InpObject<T>::freePacketLists(ObjIndex *ObjPtr )
{
     ObjIndex    *nextPtr;    
     while(ObjPtr)
     {
          nextPtr=ObjPtr->next;
          delete ObjPtr;
          ObjPtr=nextPtr;	
     } 
  return METHOD_SUCCESS; 
}



template<class T>
int InpObject<T>::freeCells( stCellIndex *cellPtr )
{
     stCellIndex    *nextPtr;    
     while(cellPtr)
     {
          nextPtr=cellPtr->next;
          delete cellPtr;
          cellPtr=nextPtr;	
     } 
  return METHOD_SUCCESS; 
}





//*** ---
template <class UCDtype>
class InpObject;



template <class T>
int AddObjToObjList( InpObject<T> *inPtr, stObject *obj, long vol,int minObjSize)
{
     stObject *ptr;
     
     if(vol<minObjSize)          /* filter test */
          return METHOD_SUCCESS; 
     else                        /* changed 11/97  */
     {   
          ptr = inPtr->pobject;
      
         if (ptr == NULL)          /* no prior object */
         {
             inPtr->pobject = obj;   /* link first object */
	     inPtr->numObj=0;
	     obj->objNum=0;
	     obj->next=NULL;
	     inPtr->numObj++;
       	     return METHOD_SUCCESS;
         }
      
         else                      /* append obj to linked list */
	 {
	      while( ptr->next != NULL )
	      {
	           ptr=ptr->next;
	      }
	      ptr->next=obj;
	      obj->next = NULL;
	      obj->objNum=inPtr->numObj;
	      inPtr->numObj++;  /* increment num of objects */
             return  METHOD_SUCCESS;
	 }
    }  
    return  METHOD_SUCCESS;
}


//extern int ReadVtkData(vtkRectilinearGrid *in_ds,int nspace,int cell_type, int &nnodes,int &ncells,int &cellpoints,float** coord_array,int** node_conn_array,float** node_data,int &nncomp);

//extern int  TrackHex(vtkRectilinearGrid *in_ds,vtkDataSet **outDS,int nspace,int cell_type, int currentTime, float threshold,std::string listfile, char* polyfile, int smallest_vol,int precision,bool dovolren,int color_Method, float firstthresh, float secondthresh, float deltax_val, float deltay_val, float deltaz_val);

//extern void resetArray(int);
//extern void setArray(int);
// ************************************************************************* //
//  END: inputdata.h
// ************************************************************************* //

#endif
