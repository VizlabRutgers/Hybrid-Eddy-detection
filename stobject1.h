#include "nodeinfo.h"
#include "stRGB.h"

// ****************************************************************************
//  Struct:      stObject
//
//  Purpose:     holds the information of the segmented object  
//
//  Programmer : Rohini Pangrikar 
//  Creation:    Wed May 28 2008
// ****************************************************************************//

struct stObject
{
  stCellIndex *cellPtr; 
  stObject    *next;
  stCellIndex *endPtr;
  long        objNum;
  long        numCell; 
  long        maxNode;
  long        minNode;
  float       objMax;
  float       objMin;
  float       minX;
  float       minY;
  float       minZ;
  float       maxX;
  float       maxY;
  float       maxZ;
  stRGB       color;
  long        volume;
  float       mass;
  float       massSqr;
  float       centroidX;
  float       centroidY;
  float       centroidZ;
  float       Ixx;
  float       Iyy;
  float       Izz;
  float       Ixy;
  float       Iyz;
  float       Izx;

  float       LowerLeft_Corner_x;
  float       LowerLeft_Corner_y;
  float       LowerLeft_Corner_z;
  float       UpperRight_Corner_x;
  float       UpperRight_Corner_y;
  float       UpperRight_Corner_z;
  long        Packet_ID;
  bool        packetFlag;
};

