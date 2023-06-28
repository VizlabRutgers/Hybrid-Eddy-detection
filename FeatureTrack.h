#ifndef FEATURE_TRACK_H
#define FEATURE_TRACK_H
#include <string>
#include "FeatureTrackUtil.h"



//these are to be shifted to another file life FeatureTrack.h
//struct Consts;


extern void OverlapTest(Frame& t1,Frame& t2,vector<vector<int> > &OverlapTable,vector<vector<int> >& ScoreBoard,vector<int>& tag1,vector<int>& tag2);
extern int ComputeScore(DIRECTION const direct,Frame& t1,Frame& t2,vector<vector<int> >& ScoreBoard,vector<vector<int> >& OverlapTable);
extern int ComputeScore2(DIRECTION const direct,Frame& t1,Frame& t2,vector<vector<int> >& ScoreBoard,vector<vector<int> >& OverlapTable);
extern void FindOverlap(int const obj,vector<int>& Overlaps, DIRECTION const direct,Frame& t1,Frame& t2,vector<vector<int> > &OverlapTable );
extern void GenCombination(vector<int>& comb,int i);
extern int Intersect(vector<int> const & Comb,int const NumOvlp,vector<int> const & Overlaps,int const obj, DIRECTION const direct,vector<vector<int> > &OverlapTable);
extern float Intersect2(int objectIndex,vector<int> const & Overlaps,int const obj, DIRECTION const direct,vector<vector<int> > &OverlapTable);
extern float GeomMean(vector<int> const & Comb,int const NumOvlp,vector<int> const & Overlaps,int const obj, DIRECTION const direct, Frame& t1, Frame& t2);
extern float GeomMean3(int objectIndex,vector<int> const & Overlaps,int const obj, DIRECTION const direct, Frame& t1, Frame& t2);
extern float GeomMean2(int objectIndex,vector<int> const & Overlaps,int const obj, DIRECTION const direct, Frame& t1, Frame& t2);
extern float Union(int objectIndex,vector<int> const & Overlaps,int const obj, DIRECTION const direct, Frame& t1, Frame& t2);
extern bool RewritePolyFiles(string const trakTable, int const step,vector<string>& time_polyfile);
extern bool ReadDagFile(string const Filename,int const step, vector<FRAME>& FrameList);
extern void ParseTrackInfo(string const buffer,TMPFR &tmpfr0,TMPFR &tmpfr1);
extern bool AddNodeToTmpFr(int const which,int const obj,vector<int> & tmpObj,TMPFR &tmpfr0,TMPFR &tmpfr1);
extern void CopyTmpFrToList(int const which,TMPFR &tmpfr,vector<FRAME>& FrameList);
extern bool ReadNumNodes(vector<FRAME>& FrameList,int step,vector<string>& time_polyfile);
extern string path_core(string const s);
extern bool Read1stPolyFile(string const polyfile,vector<FRAME>& FrameList);
extern void Colorize(vector<FRAME>& FrameList);
extern void ColorizeSuc(vector<FT_NODE>::iterator node, vector<FRAME>::iterator nextfr);
extern bool UniformColorPoly(int const step,vector<FRAME>& FrameList,vector<string>& time_polyfile);
extern void TrackSplit_Merge(int const step, SPLIT_MERGE const sm,Frame& t1,Frame &t2, FILE* outfile,vector<vector<TrackObject> >&  objs,vector<int>& tag1,vector<int>& tag2,vector<vector<int> >& ScoreBoard, double &tolerance);
extern void TrackSplit_Merge2(int const step, SPLIT_MERGE const sm,Frame& t1,Frame &t2, FILE* outfile,vector<vector<TrackObject> >&  objs,vector<int>& tag1,vector<int>& tag2,vector<vector<int> >& ScoreBoard, double tolerance);
extern void TrackContinuous(int step,Frame& t1,Frame &t2,FILE* outfile,vector<vector<TrackObject> >&  objs,vector<int>& tag1,vector<int>& tag2,vector<vector<int> >& ScoreBoard, double& tolerance);
extern void TrackNew_Dissipate(int const step, NEW_DISSIPATE const nd,Frame& t1,Frame &t2,FILE* outfile,vector<vector<TrackObject> >&  objs,vector<int>& tag1,vector<int>& tag2);
extern void NumNonZeros(int const row,vector<int> & list, ROW_COL const rc,Frame& t1,Frame &t2 ,vector<vector<int> >& ScoreBoard,vector<int>& tag1,vector<int>& tag2);
extern int GetSameScore(int const row,int const ind,int const nlist,vector<int> & list,vector<int> & split, ROW_COL const rc,vector<vector<int> >& ScoreBoard,vector<int>& tag1,vector<int>& tag2, double &tolerance);
extern int CheckSM_Condition(int const row,int const nlist,vector<int> & list,vector<int> & split, ROW_COL const rc,vector<vector<int> >& ScoreBoard,vector<int>& tag1,vector<int>& tag2, double &tolerance, int objVolume);
extern bool IsMax(int const row,int const col, ROW_COL const rc,vector<vector<int> >& ScoreBoard,vector<int>& tag1,vector<int>& tag2,Frame& t1, Frame& t2, double &tolerance);
extern void GetInts(string const buffer, vector<int> & line);
#endif
