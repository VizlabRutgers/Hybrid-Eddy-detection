#ifndef OBJSEGMENT_UTIL_H
#define OBJSEGMENT_UTIL_H

#define METHOD_SUCCESS 1
#define METHOD_FAILURE 0
#define NOT_FOUND  0
#define YES_FOUND 1

const int _PRISM_=0;
const int _HEX_=1;
const int _QUAD_=2;
//const int _UNIFORM_=0;
const int _NONUNIFORM_=1;
const int THRESH_PERCENT = 0;
const int INVALID_OBJ = -1;
const int INTERVAL = 1000;
const int ADAPTIVE = 1;
const int FIXED = 0;

#define  F_UNUSED  0x00000000
#define  F_THR     0x00000001
#define  F_USED    0x00000002
#define  F_UCD     0x00000004

#endif