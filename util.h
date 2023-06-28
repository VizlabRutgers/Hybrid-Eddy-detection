#ifndef UTIL_H
#define UTIL_H
// ************************************************************************* //
//  File: util.h
// ************************************************************************* //
#define METHOD_SUCCESS 1
#define METHOD_FAILURE 0
#define NOT_FOUND  0
#define YES_FOUND 1

//Pinakin ::added new externed constants converting all the const ints defined in os.cxx
extern const int _NNODES_;

extern const int _COORDS_;
extern const int _NODEDATA_;
extern const int _CONNECTS_;
extern const int _MIN_NODE_LIST_;
extern const int _MAX_NODE_LIST_;

extern const int MAX_CHAR;
extern const int THRESH_PERCENT;

extern const int FIXED;
extern const int ADAPTIVE;

extern const int _TETRA_;
extern const int _PRISM_;
extern const int _HEX_;

extern const int INVALID_OBJ;
extern const int INTERVAL;
//extern const int _UNIFRM_;
//extern const int _NONUNIFORM;

#define prt_rtn1(A) { printf(A); return METHOD_FAILURE; }
#define prt_rtn2(A, B) {printf(A,B); return METHOD_FAILURE;}
#define prt_rtn3(A, B, C) {printf(A,B,C); return METHOD_FAILURE;}

#define ARRFREE(a) { if(a) {ARRfree(a); a=0;} };

#define  F_UNUSED  0x00000000
#define  F_THR     0x00000001
#define  F_USED    0x00000002
#define  F_UCD     0x00000004

#define _DEBUG_LEVEL_  0 /* 0 : shut up; 1: report ; 2: ???*/
#define DEBUG_NETWORK

int ucdcompare(int, int);
void ucdsort(int *, int );
string precision_time(int const time,int const precision);
bool genfldfile(int* in_field_dims,char const * label);

//from interf/util.h


#endif
