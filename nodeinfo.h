#ifndef NODEINFO_H
#define NODEINFO_H

// ************************************************************************* //
//  File: nodeinfo.h
// ************************************************************************* //

// ****************************************************************************
//  Struct: stCellIndex
//
//  Purpose: holds the list of the cells of which the node is part of
//      
//  Programmer: Rohini Pangrikar 
//  Creation:   Wed May 28 2008
// ****************************************************************************//
struct stCellIndex 
{ 
  unsigned long index;
  struct stCellIndex *next;
};

// ****************************************************************************
//  Struct: ObjIndex
//
//  Purpose: holds the list of the cells of which the node is part of
//      
//  Programmer: Sedat OZER
//  Creation:   July 20 2010
// ****************************************************************************//
struct ObjIndex 
{ 
  unsigned long objNumber;
  struct ObjIndex *next;
};

// ****************************************************************************
//  Struct: stPacket
//
//  Purpose: holds the list of the cells of which the node is part of
//      
//  Programmer: Sedat OZER
//  Creation:   July 20 2010
// ****************************************************************************//

struct stPacket
{
  ObjIndex  *ObjList;        // this is the list for the objects in the current packet
  unsigned long       P_ID;              // this is the current packet's ID
  unsigned long       numObjs;           // total number of objects
  float      PacketCx;          // whole packet's centroid X
  float      PacketCy;          // whole packet's centroid Y
  float      PacketCz;          // whole packet's centroid Z
  float      PacketMass;          // whole packet's centroid Z
  long       PacketVolume;       // whole packet's Volume
  float      minX;               // packet extends: min X,Y,Z and Max X,Y,Z values
  float      minY;
  float      minZ;
  float      maxX;
  float      maxY;
  float      maxZ;  

  struct stPacket *next;

};



// ****************************************************************************
//  Struct: stNodePos
//
//  Purpose: holds the information of a node of the input data
//      
//  Programmer: Rohini Pangrikar 
//  Creation:   Wed May 28 2008
// ****************************************************************************//

struct stNodePos
{
  float x;
  float y;
  float z;
  
  long flag;
  stCellIndex  *list;  //loop delete
  unsigned long *adjPosList; // store only the points who have the different sign other then this point. It should better be array because binary search and quick sort will be used on it.
  unsigned long numOfAdjPos;
};
// ****************************************************************************
//  Struct: stNodeData
//
//  Purpose: holds the data value at the node of input data
//      
//  Programmer: Rohini Pangrikar 
//  Creation:   Wed May 28 2008
//
// ****************************************************************************//
struct stNodeData
{
  float  val0;
  float  val1;
  float  val2;
};

stNodePos* CreateNodeArray( unsigned long numItems);
stNodePos* GetCurNode( stNodePos nodeList[], unsigned long index );

int AddCellToNodeList( stNodePos nodeList[], unsigned long index, unsigned long cellId);
stCellIndex* GetIncidentCells( stNodePos nodeList[], unsigned long index );

stNodeData* CreateDataArray(unsigned long numItems);
stNodeData* GetCurNodeData(stNodeData dataList[],unsigned long index);



// ************************************************************************* //
//  END: nodeinfo.h
// ************************************************************************* //

#endif