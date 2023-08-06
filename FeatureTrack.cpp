/*! \file tracking.cxx
  \brief tracking functions
*/
#ifndef FEATURE_TRACK_CPP
#define FEATURE_TRACK_CPP
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <iterator>
#include <cmath>
#include <string.h>
using namespace std;
//#include<boost/tokenizer.hpp>
//using namespace boost;
#include <limits>
#include <sys/stat.h>
#include "interfaceutil.h"
//#include "./Ftrack/FeatureTrackUtil.h"
#include "FeatureTrackUtil.h"
#include "FeatureTrack.h"

extern double bounds[6];

bool TrackObjects(string const basename,int const step,int curtime,Frame& f1,Frame& f2, vector<string>& time_polyfile,vector<vector<TrackObject> > & objs) 
{ // objs are actually not used after assignment. It might be used in the future.
     vector<vector<int> > OverlapTable; // no_ext+rand_wr, init in OverlapTest()  
     vector<vector<int> > ScoreBoard; // no_ext+rand_wr, init in OverlapTest() 
     vector<int> tag1, tag2; // no_ext+rand_wr, init in OverlapTest()
  // tag1, tag2 is initialized in OverlapTest(). They are shared by TrackSplit_Merge(),TrackContinuous(),TrackNew_Dissipate()
  // They tag the elements have been handled in Frame1 and Frame2 respectively.
     double tolerance = 5.0;

//     int test=0;
//     if(curtime==5)
//         test=1;

     OverlapTest(f1,f2,OverlapTable,ScoreBoard,tag1,tag2);
  // ScoreBoard must be initialized in OverlapTest(), not in ComputScore(), otherwise, it will be init twice!
     //#cout<<"TrackObject:after OverlapTest"<<endl;
     //#cout<<"OverlapTable:\n";
     /*//#for(vector<vector<int> >::const_iterator it1=OverlapTable.begin();it1!=OverlapTable.end();++it1)
     {
          for(vector<int>::const_iterator it2=it1->begin();it2!=it1->end();++it2) 
          {
   	       cout<<*it2<<" ";
   	  }
   	       cout<<"\n";
     }*/
     
     ComputeScore2(_FORWARD_,f1,f2,ScoreBoard,OverlapTable);
     //#cout<<"TrackObject:after ComputeScore(_FORWARD_)\n";
     ComputeScore2(_BACKWARD_,f1,f2,ScoreBoard,OverlapTable);
     //#cout<<"TrackObject:after ComputeScoreBoard(_BACKWARD_)\n";
     //#cout<<"ScoreBoard:\n";
     /*//#for(vector<vector<int> >::const_iterator it1=ScoreBoard.begin();it1!=ScoreBoard.end();++it1)
     {
          for(vector<int>::const_iterator it2=it1->begin();it2!=it1->end();++it2)
          {
   		cout<<*it2<<" ";
   	  }
   	  cout<<"\n";
     }*/
     string trakTable=basename+".trakTable";
     FILE *outfile;  
     if(step==1) // step==1 means the second frame, because index starts from 0
     {
          //#cout<<"TrackObjects:open trakTable to write:"<<trakTable<<"\n";
          outfile=fopen(trakTable.c_str(),"w");  
      // If just begin tracking, then create a new trakTable file
      // If there is already a trakTable file, then truncate it to an empty file.
      // cannot be only ios_base::trunc to create a new file
     }
     else
     {
  	// else, append new tracking info to it.

  	// Don't worry, this segment is to double check whether the current Frame has been tracked
  	// If yes, return silently. Otherwise, track.
          //string const tmpstr1("Frame #"+TtoS<int>(step+1)); #rohini: write actual timestep number
	  string const tmpstr1("Frame #"+TtoS<int>(curtime));
          if(!(outfile=fopen(trakTable.c_str(), "r"))) 
               cout << "cannot open trakTable to read\n";
          while(!feof(outfile))
          {
               char tmpstr[Consts::TRAKTABLE_MAXLEN];  // hard limit, might result segment fault for large dataset! I cannot use fstream to circumvent this problem.(I tried this for a long time.) But ReadDagFile() use fstream to read in trakTable file. it is self-discrepant. I don't know why. 
               fscanf(outfile,"%[^\n]\n",tmpstr);
               if(!strcmp(tmpstr,tmpstr1.c_str()))
               {
	            fclose(outfile);
	            RewritePolyFiles(trakTable,step,time_polyfile);
	//You must rewrite polyfile even you need not append track info to the trakTable file. Because the os module will always generate a polyfile. If do not modify the colors, the polyfile will have untracked color and the rendering will be wrong.
	// I don't know why the animfilename always repeat the last frame for twice. Strange.
	            return true;
               }
          }
          fclose(outfile);
          outfile=fopen(trakTable.c_str(),"a+");
     }
     //fprintf(outfile,"Frame #%d\n",step+1); # rohini: write actual timestep
     fprintf(outfile,"Frame #%d\n",curtime); 
     //#cout<<"TrackObject:after preparation of trakTable file\n";
     
  // The three Tracks are executed in order. 
  // "In the Score_board, entries of maximum equal scores in a row indicate taht an object in t(i) break up into two
  // or more objects in t(i+1). Similarly, entries of maximum equal scores ina column, indicates that objects in t(i) merges
  //into an object in t(i+1). The entries with maximum score which are not chosen in the above steps are continuous. 
  // The remaining objects in the search space are tagged creation (for t(i+1)), and dissipation (for t(i)).

  // Modified by Weiping Hua
//     TrackSplit_Merge(step,_SPLIT_,f1,f2,outfile,objs,tag1,tag2,ScoreBoard,tolerance);
     TrackSplit_Merge2(step,_SPLIT_,f1,f2,outfile,objs,tag1,tag2,ScoreBoard,tolerance);
     //#cout<<"TrackObject:after TrackSplit\n";
//     TrackSplit_Merge2(step,_MERGE_,f1,f2,outfile,objs,tag1,tag2,ScoreBoard,tolerance);
     TrackSplit_Merge2(step,_MERGE_,f1,f2,outfile,objs,tag1,tag2,ScoreBoard,tolerance);
     //#cout<<"TrackObject:after TrackMerge\n";
     TrackContinuous(step,f1,f2,outfile,objs,tag1,tag2,ScoreBoard,tolerance);
     //#cout<<"TrackObject:after TrackContinuous\n";
     TrackNew_Dissipate(step, _NEW_,f1,f2,outfile,objs,tag1,tag2);
     //#cout<<"TrackObject:after TrackNew\n";
     TrackNew_Dissipate(step,_DISSPATE_,f1,f2,outfile,objs,tag1,tag2);
     //#cout<<"TrackObject:after TrackDissipate\n";
     fclose(outfile);
     //only rewrite the polyfile on time step.
     RewritePolyFiles(trakTable,step,time_polyfile);
     //#cout<<"TrackObject:after RewritePolyFiles\n";
     return true;
}

/*! \fn void OverlapTest(Frame& t1,Frame& t2,vector<vector<int> > &OverlapTable,vector<vector<int> >& ScoreBoard,vector<int>& tag1,vector<int>& tag2)
  \param t1 the first frame
  \param t2 the second frame
  \param OverlapTable the 2d overlap table
  \param ScoreBoard  the 2d ScoreBoard table,ScoreBoard is empty and allocated space in this function. But ScoreBoard will not be filled with numbers in this function.
  \param tag1  the 1d array 
  \param tag2  the 1d array
  \returns none
 */
void OverlapTest(Frame& t1,Frame& t2,vector<vector<int> > &OverlapTable,vector<vector<int> >& ScoreBoard,vector<int>& tag1,vector<int>& tag2)
{
    // initialize OverlapTable[numObjs1][numObjs2]
     int numObjs1=t1.objVols.size();
     int numObjs2=t2.objVols.size();
     ScoreBoard=OverlapTable = vector<vector<int> >(numObjs1, vector<int>(numObjs2, 0));
     tag1=vector<int>(numObjs1,0);
     tag2=vector<int>(numObjs2,0);
     // the total volume of each frame. (volume is actually nodes number) 
     int numNodes1=t1.nodes.size();
     int numNodes2=t2.nodes.size();
     double shift_thresh = std::pow((t1.xMax-t1.xMin)*0.001,2)+std::pow((t1.yMax-t1.yMin)*0.001,2); //use 9 for best result so far..
 //    int temp_shift_thresh = 0;
     // the while loop is extemely time consuming. I use [] intead of at() when accessing vectors.
     int i(0), j(0);
     //#cout<<"OverlapTest:numNodes1="<<numNodes1<<" numNodes2="<<numNodes2<<endl;
     while (i<numNodes1 || j<numNodes2)
     {

/*          if (shift_thresh > j)
          { 
            temp_shift_thresh = (int)(j/2);
          }
          else
          { 
            temp_shift_thresh = shift_thresh;
          } 
*/

// Edited by SEDAT


          //(t1.nodes[i].NodeID == t2.nodes[j].NodeID)   /* overlap */    //     (t1.nodes[i].NodeID == t2.nodes[j].NodeID)


         //          int tempval =  (t2.nodes[j].NodeID) - (t1.nodes[i].NodeID);
         // Edit by Weiping Hua
         // It's insane to use NodeID. Now we use the coordinates.
          double tempval =  std::pow(t2.nodes[j].xCoord-t1.nodes[i].xCoord,2)+std::pow(t2.nodes[j].yCoord-t1.nodes[i].yCoord,2)+std::pow(t2.nodes[j].zCoord-t1.nodes[i].zCoord,2);


          //          if ( tempval <= shift_thresh)
//          if ( (tempval <= shift_thresh) && ( tempval >= 0) ) //(t1.nodes[i].NodeID == t2.nodes[j].NodeID) 
//          if (t1.nodes[i].NodeID == t2.nodes[j].NodeID)  //(t1.nodes[i].NodeID == t2.nodes[j].NodeID) 
	  //if ( (tempval <= shift_thresh) && ( tempval >= 0) ) //(t1.nodes[i].NodeID == t2.nodes[j].NodeID)
          //if (t1.nodes[i].NodeID == t2.nodes[j].NodeID)
         if ( (tempval <= shift_thresh) && ( tempval >= 0) )
          {
//             ucdNode test1;
//             ucdNode test2;
//             int test=0;
//             if(t1.nodes[i].ObjID ==27 && t2.nodes[j].ObjID==31){
//                 test=1;
//                 test1 = t1.nodes[i];
//                 test2 = t2.nodes[j];
//             }
               OverlapTable[t1.nodes[i++].ObjID][t2.nodes[j++].ObjID]++;
                 //cout<<"An entry entered to the OVERLAPTABLE !!"<< endl;
          }
          else
          {
	       if(t1.nodes[i].NodeID>t2.nodes[j].NodeID) 
	       {
                    if (j<numNodes2)
                         j++;
               }
	       else if (i<numNodes1)
                    i++;
          }
      
          if (i==numNodes1 || j==numNodes2)
          { /* one frame finished, then terminate search */
	       i = numNodes1;
	       j = numNodes2;
          }

     }

}

int ComputeScore(DIRECTION const direct,Frame& t1,Frame& t2,vector<vector<int> >& ScoreBoard,vector<vector<int> >& OverlapTable)
{
  vector<int> Comb; // no_ext+rand_wr, init in GenCombination()
  vector<int> Overlaps; //no_ext+rand_wr, init in FindOverlap()

/* Overlaps looks like:
    3 5 6 8 
   It is constructed with push_back(). Each element is the index of obj it is overlapping.
   */
/* OverlapTable looks like:
          obj2
         0      1       2
      ___________________
      0  20   1000    204
      1  36    485     56
obj1  2  504    56      7
      3  23     3      13
      __________________
   Each element in this 2d array is the overlap volume (node number) of two objs.
   OverlapTable is not symmetric.
   */
/* ScoreBoard looks like OverlapTable. Each element in this 2d array is the score of two objs.
   It is not symmetric. OverlapTable[i][j]!=0 <=> ScoreBoard[i][j]. ScoreBoard is actually an overlapTable after
   processing. An element in OverlapTable is the raw overlapping volume between two objects in different frames.
   ScoreBoard's elements are obtained after overall considerated. A big element in OverlapTable does not mean 
   there is a corresponding big element in the ScoreBoard.
   */
   
     int numObjs1=t1.objVols.size();
     int numObjs2=t2.objVols.size();
     int numObjs=(direct==_FORWARD_) ? numObjs1:numObjs2;
     for (register int obj1 = 0; obj1<numObjs; obj1++) 
     {
          Overlaps.clear();
          FindOverlap(obj1, Overlaps, direct,t1,t2,OverlapTable); // store the overlapping object indexes into Overlaps vector
          int NumOvlp(Overlaps.size()); 
          int NumCom( (int) pow(2.0, NumOvlp) );
          Comb.resize(NumOvlp);
          for (register int i=1; i<NumCom; i++)
          {
              // By Weiping Hua
              // This combination is to calculate all possible combination (00/01/10/11) and pick up the biggest score. I guess it is because the previous
              // method would regard this as a merge/split by checking if the minimum difference is lower than the tolerence. However, in most cases, every
              // feature needs to pass the tolerence instead of finding the highest score.
               GenCombination(Comb, i); 
               float cost(  (float) Intersect(Comb, NumOvlp, Overlaps, obj1, direct,OverlapTable)/GeomMean(Comb, NumOvlp, Overlaps, obj1, direct,t1,t2) );
               int Score((int)(1000*cost));
               for (register int j=0; j<NumOvlp; j++)
               {
	            if (Comb.at(j))
                    {
	                 int obj2 = Overlaps.at(j);
	                 if(direct==_FORWARD_)
                         {
    	                      if ((Score>ScoreBoard.at(obj1).at(obj2)) && (cost>Consts::DEFAULT_TOLERANCE)) 
	                           ScoreBoard.at(obj1).at(obj2) = Score;
	  	         }
                         else
                         {
	  	              if ((Score>ScoreBoard.at(obj2).at(obj1)) && (cost>Consts::DEFAULT_TOLERANCE)) 
	                           ScoreBoard.at(obj2).at(obj1) = Score;
	  		 }
    	            }
               }
          }
     }
     return 1;
}

int ComputeScore2(DIRECTION const direct,Frame& t1,Frame& t2,vector<vector<int> >& ScoreBoard,vector<vector<int> >& OverlapTable)
{
  vector<int> Comb; // no_ext+rand_wr, init in GenCombination()
  vector<int> Overlaps; //no_ext+rand_wr, init in FindOverlap()


/* Modified by Weiping Hua*/

/* Overlaps looks like:
    3 5 6 8
   It is constructed with push_back(). Each element is the index of obj it is overlapping.
   */
/* OverlapTable looks like:
          obj2
         0      1       2
      ___________________
      0  20   1000    204
      1  36    485     56
obj1  2  504    56      7
      3  23     3      13
      __________________
   Each element in this 2d array is the overlap volume (node number) of two objs.
   OverlapTable is not symmetric.
   */
/* ScoreBoard looks like OverlapTable. Each element in this 2d array is the score of two objs.
   It is not symmetric. OverlapTable[i][j]!=0 <=> ScoreBoard[i][j]. ScoreBoard is actually an overlapTable after
   processing. An element in OverlapTable is the raw overlapping volume between two objects in different frames.
   ScoreBoard's elements are obtained after overall considerated. A big element in OverlapTable does not mean
   there is a corresponding big element in the ScoreBoard.
   */

    int numObjs1=t1.objVols.size();
    int numObjs2=t2.objVols.size();
    int numObjs=(direct==_FORWARD_) ? numObjs1:numObjs2;
    for (int obj1 = 0; obj1<numObjs; obj1++){
        Overlaps.clear();
        FindOverlap(obj1, Overlaps, direct,t1,t2,OverlapTable); // store the overlapping object indexes into Overlaps vector
        int NumOvlp(Overlaps.size());
        int NumCom( (int) pow(2.0, NumOvlp) );
        Comb.resize(NumOvlp);
        for (int i=0; i<NumOvlp; i++){
            //               GenCombination(Comb, i);
            float cost(Intersect2(i, Overlaps, obj1, direct,OverlapTable)/GeomMean3(i, Overlaps, obj1, direct,t1,t2) );
            int Score((int)(1000*cost));
            int obj2 = Overlaps.at(i);
            if(direct==_FORWARD_){
               if((Score>ScoreBoard.at(obj1).at(obj2)) && (cost>Consts::DEFAULT_TOLERANCE))
                   ScoreBoard.at(obj1).at(obj2) = Score;
            }
            else{
                if((Score>ScoreBoard.at(obj2).at(obj1)) && (cost>Consts::DEFAULT_TOLERANCE))
                    ScoreBoard.at(obj2).at(obj1) = Score;
            }
        }
    }
    return 1;
}


/*! \fn void FindOverlap(int const obj,vector<int>& Overlaps, DIRECTION const direct,Frame& t1,Frame& t2,vector<vector<int> > &OverlapTable )
  \brief With OverlapTable, find the overlaps for an object.
  \param obj the index of an object
  \param Overlaps the 1d array storing the overlaps of an object
  \param direct  can only be _FORWARD_ or _BACKWARD_
  \param t1    the first frame
  \param t2    the second frame
  \param OverlapTable the 2d table computed in OverlapTest()
  \returns none
 */
void FindOverlap(int const obj,vector<int>& Overlaps, DIRECTION const direct,Frame& t1,Frame& t2,vector<vector<int> > &OverlapTable )
{
     int numObjs;
     if(direct==_FORWARD_)
          numObjs=t2.objVols.size();
     else
          numObjs=t1.objVols.size();
     for (register int i(0); i<numObjs; i++)
     {
          if(direct==_FORWARD_)
          {
	       if (OverlapTable.at(obj).at(i)) 
	       Overlaps.push_back(i);
          }
          else if(OverlapTable.at(i).at(obj))
	       Overlaps.push_back(i);
     }
}


void GenCombination(vector<int>& comb,int i)
{
     for(vector<int>::iterator it=comb.begin();it!=comb.end();++it) 
          *it=0;
     int j(0);
     while (i > 0) 
     {
          comb.at(j++)=(i % 2); 
          i /= 2;
     }
}

/*! \fn int Intersect(vector<int> const & Comb,int const NumOvlp,vector<int> const & Overlaps,int const obj, DIRECTION const direct,vector<vector<int> > &OverlapTable)
  \brief Get the total overlapping volume between obj in this frame and objs in the other frame. 
  \param Comb the 1d combination table computed in GenCombination()
  \param NumOvlp the size of Comb
  \param Overlaps the 1d array storing Overlaps 
  \param obj the index of object to compute overlap
  \param direct  can only be _FORWARD_ or _BACKWARD_
  \return the intersecting volume between the object in this frame and objects in the other frame.
 */
// "total" with Comb ???
int Intersect(vector<int> const & Comb,int const NumOvlp,vector<int> const & Overlaps,int const obj, DIRECTION const direct,vector<vector<int> > &OverlapTable)
{
     int InterVol(0);
     for (int i=0; i<NumOvlp; i++)
     {
          if (Comb.at(i))
          {
	       int obj1 = Overlaps.at(i);
	       if(direct==_FORWARD_)
	            InterVol += OverlapTable.at(obj).at(obj1);
	       else
                    InterVol += OverlapTable.at(obj1).at(obj);
	  }
    }
    return InterVol;
}
// Modified by Weiping Hua
float Intersect2(int objectIndex,vector<int> const & Overlaps,int const obj, DIRECTION const direct,vector<vector<int> > &OverlapTable)
{
    float InterVol(0);
    int obj1 = Overlaps.at(objectIndex);
    if(direct==_FORWARD_)
        InterVol += OverlapTable.at(obj).at(obj1);
    else
        InterVol += OverlapTable.at(obj1).at(obj);
    return InterVol;
}
//  
/*! \fn float GeomMean(vector<int> const & Comb,int const NumOvlp,vector<int> const & Overlaps,int const obj, DIRECTION const direct, Frame& t1, Frame& t2)
  \brief compute the geometry mean for one object.
  \param Comb the 1d combination table computed in GenCombination
  \param NumOvlp the size of Comb
  \param Overlaps the 1d array storing the overlaps of this obj.
  \param obj   the current object
  \param direct can only be _FORWARD_ or _BACKWARD_
  \param t1    the first frame
  \param t2    the second frame
  \returns the mean
*/   
float GeomMean(vector<int> const & Comb,int const NumOvlp,vector<int> const & Overlaps,int const obj, DIRECTION const direct, Frame& t1, Frame& t2)
{
     long vol(0);
     for (register int i=0; i<NumOvlp; i++)
     {
          if (Comb.at(i))
          {
           int obj1 = Overlaps.at(i);
	       if(direct==_FORWARD_)
	            vol += t2.objVols.at(obj1).ObjVol;    
	       else
                    vol += t1.objVols.at(obj1).ObjVol; 
          }
    }
    float mean;
    if(direct==_FORWARD_)
         mean = static_cast<float>(sqrt(double(t1.objVols.at(obj).ObjVol*vol)));
    else
         mean = static_cast<float>(sqrt(double(t2.objVols.at(obj).ObjVol*vol)));
    return mean;
}
float GeomMean2(int objectIndex,vector<int> const & Overlaps,int const obj, DIRECTION const direct, Frame& t1, Frame& t2)
{
    long vol(0);

    int obj1 = Overlaps.at(objectIndex);
    if(direct==_FORWARD_)
        vol += t2.objVols.at(obj1).objSurfVol;
    else
        vol += t1.objVols.at(obj1).objSurfVol;
    float mean;
    if(direct==_FORWARD_)
         mean = static_cast<float>(sqrt(double(t1.objVols.at(obj).objSurfVol*vol)));
    else
         mean = static_cast<float>(sqrt(double(t2.objVols.at(obj).objSurfVol*vol)));
    return mean;
}

float GeomMean3(int objectIndex,vector<int> const & Overlaps,int const obj, DIRECTION const direct, Frame& t1, Frame& t2)
{
    long vol(0);

    int obj1 = Overlaps.at(objectIndex);
    if(direct==_FORWARD_)
        vol += t2.objVols.at(obj1).ObjVol;
    else
        vol += t1.objVols.at(obj1).ObjVol;
    float mean;
    if(direct==_FORWARD_)
         mean = static_cast<float>(sqrt(double(t1.objVols.at(obj).ObjVol*vol)));
    else
         mean = static_cast<float>(sqrt(double(t2.objVols.at(obj).ObjVol*vol)));
    return mean;
}


float Union(int objectIndex,vector<int> const & Overlaps,int const obj, DIRECTION const direct, Frame& t1, Frame& t2)
{
    long vol(0);

    int obj1 = Overlaps.at(objectIndex);
    if(direct==_FORWARD_)
        vol += t2.objVols.at(obj1).objSurfVol;
    else
        vol += t1.objVols.at(obj1).objSurfVol;
    float mean;
    if(direct==_FORWARD_)
         mean = static_cast<float>(sqrt(double(t1.objVols.at(obj).objSurfVol*vol)));
    else
         mean = static_cast<float>(sqrt(double(t2.objVols.at(obj).objSurfVol*vol)));
    return mean;
}




/*! bool RewritePolyFiles(string const trakTable, int const step,vector<string>& time_polyfile)
  \brief rewrite only current polyfile which just has been done with objseg 
*/
bool RewritePolyFiles(string const trakTable, int const step,vector<string>& time_polyfile)
{
     vector<FRAME> FrameList; 
     // only read in two frames: current tracked frame and the previous frame.
     ReadDagFile(trakTable, step, FrameList);
     //#cout<<"RewritePolyFiles:after ReadDagFiles\n";
     // Read object volumes for the two frames.
     // Used to decide the color when two objects merge.
     ReadNumNodes(FrameList,step,time_polyfile);
     //#cout<<"RewritePolyFiles:after ReadNumNodes\n";
     // read color for the first frame.
    if(!Read1stPolyFile(time_polyfile.at(step-1),FrameList))
          cout << "Error reading poly\n";
     //#cout<<"RewritePolyFiles:after Read1stPolyFiles\n";
     Colorize(FrameList);
     //#cout<<"RewritePolyFiles:after Colorize\n";
     UniformColorPoly(step,FrameList,time_polyfile);
     //#cout<<"RewritePolyFiles:after UniformColorPoly\n";
     return true;
}

/*! \fn bool ReadDagFile(string const Filename,int const step, vector<FRAME>& FrameList)
  \brief read one trakTable file and build a FrameList which only has two frames: the current frame and the previous one. This function is different from BuildNodeRelationGraph(), which is also a reading-Dag function. This function only read in the tracking relation of the two frames just tracked, to recolor the polygon file, while miss other tracking info. In the other hand, BuildNodeRelationGraph() reads in all tracking info to build a graph. */
bool ReadDagFile(string const Filename,int const step, vector<FRAME>& FrameList)
{ // before this function, FrameList should be empty
     ifstream fp(Filename.c_str());
     if(!fp.is_open())
          cout << "cannot open Dag file\n";
  
     string buffer;
     bool tmpfr0isNull(true);
     int numFrames(0);
     TMPFR tmpfr0;
     TMPFR tmpfr1;
     while(1)
     {
         getline(fp, buffer,'\n');
         if(!fp.good())
              break;
         if (buffer[0] == 'F')
     // indicates start of a new frame
         {
	      if(!tmpfr0isNull) 
	      if(numFrames==step)
                  break;
	
	      numFrames += 1;
	      tmpfr0isNull = false;
	
	      continue;
          }
          if(numFrames==step)
          ParseTrackInfo(buffer,tmpfr0,tmpfr1); //<<<<<<<
     } //end of while loop
  
  // process the last frame
     FrameList.push_back(FRAME(numFrames-1));
     CopyTmpFrToList(0,tmpfr0,FrameList); 
     FrameList.push_back(FRAME(numFrames));
     CopyTmpFrToList(1,tmpfr1,FrameList);
     fp.close();
      //verify FrameList
     //#cout<<"ReadDagFile():verify FrameList:\n";
    /*//# for(vector<FRAME>::const_iterator it=FrameList.begin();it!=FrameList.end();++it) 
     {
          //#cout<<"Frame index:"<<it->index<<" "; 
          //#cout<<"this frame has "<<it->NodeArray.size()<<" nodes\n";
          int i(1);
          for(vector<FT_NODE>::const_iterator it1=it->NodeArray.begin();it1!=it->NodeArray.end();++it1,++i) 
          {
               cout<<"node "<<i<<":"<<it1->getnumSuc()<<" Sucs(";
               for(vector<int>::const_iterator it2=it1->Suc.begin();*it2!=-1;++it2)
	            cout<<*it2<<" ";
               cout<<") "<<it1->getnumPre()<<" Pres(";
               for(vector<int>::const_iterator it2=it1->Pre.begin();*it2!=-1;++it2)
	            cout<<*it2<<" ";
               cout<<")\n";
          }
     }*/
     return true;
}

void ParseTrackInfo(string const buffer,TMPFR &tmpfr0,TMPFR &tmpfr1)
{
     bool minusFlag(false);
     vector<int> tmpObj(Consts::MAXSPLIT,-1); // ext+rand_wr, extend with resize(n,-1)
  //tmpObj stores the frame1's(frame1, frame2 are being tracked) objects occuring in this buffer. It helps AddNodeToTmpFr().
     vector<int> line;
     GetInts(buffer,line);
     //#cout<<"ParseTrackInfo:in ParseTrackInfo:verify line:";
     //#cout<<"(";
     /*//#for(vector<int>::const_iterator it=line.begin();it!=line.end();++it)
          cout<<*it<<",";*/
     //#cout<<")\n";
     int i(0);
     for(vector<int>::const_iterator it=line.begin(); it!=line.end(); ++it)
     {
          if (*it == -1)
               minusFlag = true;
          else if(minusFlag == true) 
               AddNodeToTmpFr(1,*it-1,tmpObj,tmpfr0,tmpfr1);
          else
          {
               AddNodeToTmpFr(0,*it-1,tmpObj,tmpfr0,tmpfr1); // The frame into which the Node will be added has been initialized. AddNodeToTmpFr(0...) does not use tmpObj.
       
       // add the number before -1 into tmpObj
               vector<int>::const_iterator it1=tmpObj.begin();
               for(it1+=i;it1!=tmpObj.end()&&*it1!=-1; ++it1,++i);
	            bound_check<int>(tmpObj,i,Consts::MAXSPLIT,-1);
                    tmpObj.at(i) = *it -1;
          }
     }
}

/*! bool AddNodeToTmpFr(int const which,int const obj,vector<int> & tmpObj,TMPFR &tmpfr0,TMPFR &tmpfr1)
  \brief add an object node to tmpfr0 or tmpfr1.
 */
bool AddNodeToTmpFr(int const which,int const obj,vector<int> & tmpObj,TMPFR &tmpfr0,TMPFR &tmpfr1)
{
  // 1. add the obj into tmpfr1 array
  // 2. add the obj into the corresponding obj's successor array in tmpfr0 array
  // obj starts from 0

  // tmpObj should like "2 5 3 6 1 -1 -1 -1 -1 -1"
     if(!which)
     {
          bound_check<TMPNODE>(tmpfr0.Obj,obj,Consts::MAXOBJS);
          tmpfr0.Obj.at(obj).ObjInd = obj;
          return true;
     }
     // which==1
     bound_check<TMPNODE>(tmpfr1.Obj,obj,Consts::MAXOBJS);
     tmpfr1.Obj.at(obj).ObjInd = obj;
  // for tmpfr1, we need not process tmpfr1.Obj[i].Suc[].

  //(1)add obj as a successor into the suc list of every object in tmpfr0
  //(2)add objs in tmpObj as the preceeds into the preceeds list of obj in tmpfr1.
  //NOTE: (2) is used in BldNodeRelateGraph(), but not useful in tracking.cxx.

  // here tmpObj has been set and fixed.
     for(vector<int>::iterator it=tmpObj.begin(); it!=tmpObj.end()&&*it!=-1; ++it)
     {
          int i= tmpfr0.Obj.at(*it).getnumSuc();  
      // correct. See Feature.h about Suc.
          bound_check<int>(tmpfr0.Obj.at(*it).Suc,i,Consts::MAXSPLIT,-1);
          tmpfr0.Obj.at(*it).Suc.at(i) = obj;
          i=tmpfr1.Obj.at(obj).getnumPre();
          bound_check<int>(tmpfr1.Obj.at(obj).Pre,i,Consts::MAXSPLIT,-1);
          tmpfr1.Obj.at(obj).Pre.at(i)=*it;
     }
     return true;
 }
void CopyTmpFrToList(int const which,TMPFR &tmpfr,vector<FRAME>& FrameList)
{  // which = 0, copy tmpfr0;  1 copy tmpfr1;
     int t(tmpfr.getnumObjs());
     FrameList.back().NodeArray.resize(t);
     register int i = 0;
     for(i=0; i<t; i++)
     {
          FrameList.back().NodeArray.at(i).FrameInd = FrameList.back().index;
     }

  // Because I only track color in forward direction, only Suc is considered.
  // Thus, only which==0 is handled.
  // Although FT_NODE has Pre member, it is not set here.
     if( which == 0) 
     {
          for( i=0; i<t; i++)
          {
               FrameList.back().NodeArray.at(i).Suc.swap(tmpfr.Obj.at(i).Suc);
               cout<<"CopyTmpFrToList:frame"<<FrameList.size()-1<<":obj"<<i<<":"<<FrameList.back().NodeArray.at(i).getnumSuc()<<"suc:";
      //    copy(FrameList.back().NodeArray.at(i).Suc.begin(), FrameList.back().NodeArray.at(i).Suc.end(),ostream_iterator<int>(cout," "));
               for(vector<int>::iterator it=FrameList.back().NodeArray.at(i).Suc.begin();it!=FrameList.back().NodeArray.at(i).Suc.end() && *it!=-1;++it)          cout<<*it<<" ";
               cout<<"\n";

          }
     }
     else
     { // which==1
          for( i=0; i<t; i++)
          {
               FrameList.back().NodeArray.at(i).Pre.swap(tmpfr.Obj.at(i).Pre);
               cout<<"CopyTmpFrToList:frame"<<FrameList.size()-1<<":obj"<<i<<":"<<FrameList.back().NodeArray.at(i).getnumPre()<<"pre:";
      //    copy(FrameList.back().NodeArray.at(i).Suc.begin(), FrameList.back().NodeArray.at(i).Suc.end(),ostream_iterator<int>(cout," "));
               for(vector<int>::iterator it=FrameList.back().NodeArray.at(i).Pre.begin();it!=FrameList.back().NodeArray.at(i).Pre.end() && *it!=-1;++it)          cout<<*it<<" ";
               cout<<"\n";
          }
     }
}
/*! \fn bool ReadNumNodes(vector<FRAME>& FrameList,int step,vector<string>& time_polyfile)
  \brief Find volume(num of nodes) for each object node in each frame. Actually FrameList has only two frames(current and previous). 
*/
bool ReadNumNodes(vector<FRAME>& FrameList,int step,vector<string>& time_polyfile)
{
     register int i(0);
     int dummy1 = 0, dummy2 = 0, nNode = 0;
     ifstream fp;
     string attributefile, tmpstr, buffer;
     for(vector<FRAME>::iterator it=FrameList.begin(); it!=FrameList.end(); ++it, ++i)
     {
          string tmpstr(time_polyfile.at(step-1+i));
          attributefile=path_core(tmpstr)+".attr";
          ifstream fp(attributefile.c_str());
          if(!fp.is_open())
          {
	       cout<<"cannot open attr file:"<<attributefile<<endl;
	       return false;
          }
          for(vector<FT_NODE>::iterator itr=it->NodeArray.begin(); itr!=it->NodeArray.end();++itr)
	  { // take easy, NodeArray is made sure to be full
	       tmpstr.clear();
	       while(tmpstr!="Volume:")
	       {
	            getline(fp,buffer,'\n');
	            if(!fp.good())
                         break;
	            tmpstr=buffer.substr(0,buffer.find(' '));
	       }
	       tmpstr=buffer.substr(buffer.find(' ')+1);
	       itr->numNodes=StoT<int>(tmpstr);
    	  }
     }
     //#cout<<"ReadNumNodes:verify FrameList numNodes:\n";
     i=0;
    /*//# for(vector<FRAME>::iterator it=FrameList.begin(); it!=FrameList.end(); ++it, ++i)
     {
          cout<<"Frame "<<i<<" numNodes are:";
          int j(0);
         for(vector<FT_NODE>::iterator itr=it->NodeArray.begin(); itr!=it->NodeArray.end();++itr,++j)
              cout<<"(node "<<j<<","<<itr->numNodes<<")";
         cout<<"\n";
     } */
     return true;
}

string path_core(string const s)
{
  string::size_type i=s.rfind(".");
  if(i==string::npos) return s;
  else return s.substr(0,i);
}

/*! \fn bool Read1stPolyFile(string const polyfile,vector<FRAME>& FrameList)
  \brief read color from the first polyfile to assign color to the first frame. 
 */
bool Read1stPolyFile(string const polyfile,vector<FRAME>& FrameList)
{
     FILE *fp;
     char buffer[255];
     int  color255[3];
     int CurObjConnNum = 0, CurObjNodeNum = 0;
     int i, j=0;
     fp = fopen(polyfile.c_str(), "r");
     if(!fp)
         return 0;
     while(!feof(fp) && j<FrameList.begin()->NodeArray.size())
     {
          CurObjNodeNum =0; 
          CurObjConnNum =0;
          fscanf(fp, "%3d %3d %3d\n", &color255[0], &color255[1], &color255[2]);
          FrameList.begin()->NodeArray[j].Color[0]=color255[0];
          FrameList.begin()->NodeArray[j].Color[1]=color255[1];
          FrameList.begin()->NodeArray[j++].Color[2]=color255[2];
          fscanf(fp, "%d\n", &CurObjNodeNum);
          for(i=0; i<CurObjNodeNum; i++)
               fgets(buffer, 50, fp);
          fscanf(fp, "%d\n", &CurObjConnNum);
          for(i=0; i<CurObjConnNum; i++)
               fgets(buffer, 50, fp);
          fgets(buffer, 10, fp);
          fgets(buffer, 10, fp);
     }
     fclose(fp);
     return true;
}

/*! \fn void Colorize(vector<FRAME>& FrameList)
  \brief  Color the objects in the frames from the first frame 
*/
void Colorize(vector<FRAME>& FrameList)
{ // I need not set the first frame color, which has been set in Read1stPolyFile()
// if you want to apply any color policy, do in Read1stPolyFile()
  
     int k(0),siz(0);
   //Pinakin :: Checking number of frames in list
     //#cout<<"\nTHe number of frames contained in the FrameList are :"<<FrameList.size();
   
     for(vector<FT_NODE>::iterator it=FrameList.begin()->NodeArray.begin(); \
  	it!=FrameList.begin()->NodeArray.end(); ++it)
     {
          if(it->getnumSuc() != 0)
          {
      //Pinakin :: I want to display successors here for debugging
     // cout<<"\nCurrent Object:Obj"<<k<<" Suc are :";
     // for(vector<int>::const_iterator iter=it->Suc.begin();*iter!=-1;++iter)
       //  cout<<*iter<<" ";
               vector<FRAME>::iterator nextfr = FrameList.begin();
               //#cout<<"\nObject No.:"<<k<<"Color :"<<it->Color[0]<<" "<<it->Color[1]<<" "<<it->Color[2];
               ColorizeSuc(it, ++nextfr);
               k++;
          }
     }
     //Pinakin :: Check if frame-to-be-tracked has complete color information 
     //#cout<<"\nDisplaying just filled color info for frame 2:";
     FRAME tempo=FrameList.back();
     register int ll=0;
     while(ll < tempo.NodeArray.size())
     {
          //#cout<<"\nOBject No:"<<ll<<endl;
          //#cout<<tempo.NodeArray.at(ll).Color[0]<<" "<<tempo.NodeArray.at(ll).Color[1]<<" "<<tempo.NodeArray.at(ll).Color[2];
          ll++;    
     }
}

/*! \fn bool UniformColorPoly(int const step,vector<FRAME>& FrameList,vector<string>& time_polyfile)
  \brief uniform the polyfile of time "step". Because it is called after Colorize(), FrameList has already the color info. Each pass of avsexp network, only one polyfile is recolored with this function. FrameList only has two frames.
*/
bool UniformColorPoly(int const step,vector<FRAME>& FrameList,vector<string>& time_polyfile)
{ // PROBLEM!!
     register int i =0, j =0,k=0;
     char buffer[255];
     char tmpfile[100];
     FILE *fp1=NULL, *fp2=NULL;
     int CurObjConnNum = 0, CurObjNodeNum = 0;
     int color255[3];
     struct stat stat_buf;
     if(stat(time_polyfile.at(step).c_str(),&stat_buf)==-1)                   cout << "error when stat()\n";
     if(!(fp1 = fopen(time_polyfile.at(step).c_str(), "r")))
         cout << "polyfile open error\n";
     strcpy(tmpfile, time_polyfile.at(step).c_str());
     strcat(tmpfile, "tmp");
     if(!(fp2 = fopen(tmpfile, "w")))
          cout << "polytmpfile open error\n";
     j=0;
     FRAME last=FrameList.back();
  //while(!feof(fp1) && j<last.NodeArray.size())  {
     while(j<last.NodeArray.size())
     {
          CurObjNodeNum =0;
          CurObjConnNum =0;
          fscanf(fp1, "%3d %3d %3d\n", &color255[0], &color255[1], &color255[2]);
          if(last.NodeArray.at(j).Color[0]!=-1 ||
          last.NodeArray.at(j).Color[1]!=-1 ||
          last.NodeArray.at(j).Color[2]!=-1) 
          {
               color255[0] =last.NodeArray.at(j).Color[0];
               color255[1] =last.NodeArray.at(j).Color[1];
               color255[2] =last.NodeArray.at(j).Color[2];
          }
          fprintf(fp2, "%3d %3d %3d\n", color255[0], color255[1], color255[2]);
          j++;
    // break;
          fscanf(fp1, "%d\n", &CurObjNodeNum);
          fprintf(fp2, "%d\n", CurObjNodeNum);
          for(k=0; k<CurObjNodeNum; k++)
          {
      //initbuffer(buffer, 255);
               fgets(buffer, 255, fp1);
               fputs(buffer, fp2);
          }
          fscanf(fp1, "%d\n", &CurObjConnNum);
          fprintf(fp2, "%d\n", CurObjConnNum);
          for(k=0; k<CurObjConnNum; k++)
          {
      //initbuffer(buffer, 255);
               fgets(buffer, 50, fp1);
               fputs(buffer, fp2);
          }
    //initbuffer(buffer, 255);
          fgets(buffer, 10, fp1);
          fputs(buffer, fp2);
    //initbuffer(buffer, 255);
          fgets(buffer, 10, fp1);
          fputs(buffer, fp2);
     }
     fclose(fp1);
     fclose(fp2);
     if(stat(tmpfile,&stat_buf)==-1)
         cout << "error when stat(tmpfile)\n";
 //Pinakin :: Save originally written polyfile here for checking
    // char tempstore[100]; *Rohini
   //  strcpy(tempstore,tmpfile); *Rohini
   //  strcat(tempstore,"555"); *Rohini
   //  rename(time_polyfile.at(step).c_str(),tempstore); *Rohini
     ////remove(time_polyfile.at(step).c_str()); 
     // Rohini has commented these lines. It was just storing the temp poly file which is waste of disk space
     rename(tmpfile, time_polyfile.at(step).c_str());
     if(stat(time_polyfile.at(step).c_str(),&stat_buf)==-1)                   cout << "error when stat(time_polyfile.at(step))\n";
     return true;
}


/*! \fn void ColorizeSuc(vector<FT_NODE>::iterator node, vector<FRAME>::iterator nextfr)
  \brief a recursion function assigning tracked color to the related object nodes.
 */
void ColorizeSuc(vector<FT_NODE>::iterator node, vector<FRAME>::iterator nextfr)
{
  // the node in parameter is the node in Predecessor frame
  // the nextfr points to the Successor frame

  assert(node->getnumSuc() != 0);
  for(vector<int>::const_iterator it=node->Suc.begin(); *it!=-1; ++it) { 
  	// node.Suc.begin() must not be -1, this has been checked in Colorize()
  	// Because PredMaxnumNodes is initialized to 0, this if branch will always be gone into.
     if((*node).numNodes>nextfr->NodeArray.at(*it).PredMaxnumNodes) 
      { // numNodes is actually the volume of this object
    	nextfr->NodeArray.at(*it).Color[0] = node->Color[0];
     	nextfr->NodeArray.at(*it).Color[1] = node->Color[1];
	nextfr->NodeArray.at(*it).Color[2] = node->Color[2];
	nextfr->NodeArray.at(*it).PredMaxnumNodes = node->numNodes;
	  
      }
      //Pinakin :: To debug I am displaying nodes/succ and their colors	
	/*//#cout<<"\nSuccessor:"<<*it;
	cout<<"\nColorR :"<<nextfr->NodeArray.at(*it).Color[0]<<" ";
	cout<<"ColorG :"<<nextfr->NodeArray.at(*it).Color[1]<<" ";
	cout<<"ColorB :"<<nextfr->NodeArray.at(*it).Color[2]<<" ";*/
      
    
     if(nextfr->NodeArray.at(*it).getnumSuc() != 0) {
       vector<FRAME>::iterator FramePtr = nextfr+1;
       vector<FT_NODE>::iterator itr=nextfr->NodeArray.begin();
       itr+=*it;
       ColorizeSuc(itr, FramePtr);
    }
  }
}

/*! void TrackSplit_Merge(int const step, SPLIT_MERGE const sm,Frame& t1,Frame &t2, FILE* outfile,vector<vector<Obj> >&  objs,vector<int>& tag1,vector<int>& tag2,vector<vector<int> >& ScoreBoard)
  \brief track split or merge events
*/
void TrackSplit_Merge(int const step, SPLIT_MERGE const sm,Frame& t1,Frame &t2, FILE* outfile,vector<vector<TrackObject> >&  objs,vector<int>& tag1,vector<int>& tag2,vector<vector<int> >& ScoreBoard, double& tolerance)
{
     int numObjs, nSplits(0); // nSplits counts the times of splits/merges, for debug purpose.
     if(sm==_SPLIT_) 
          numObjs=t1.objVols.size(); 
     else 
          numObjs=t2.objVols.size(); 
     string funcname;
     if(sm==_SPLIT_)
          funcname="TrackSplit_Merge(SPLIT)";
     else
          funcname="TrackSplit_Merge(MERGE)";


     /*
     // section added by Naveen to find max and min values in  
     // scoreboard-------------------------------------
     // finding the maximum and minimum values in the scoreboard and calculate tolerance which is some %  
     // of difference between max and min
     // this tolerance is passed as an argument to GetSameScore function . the objects are considered to 
     // be split or merged if objects in same row or column
     // are equal for lie between this tolerance level. 
     */

     int max = 0;
     int min = 1000;
 
     for(int i = 0; i< ScoreBoard.size() ; i++)
     {
       for(int j = 0; j<ScoreBoard[i].size(); j++)
       { 
         if(ScoreBoard[i][j] > max)
            max = ScoreBoard[i][j];
         if(ScoreBoard[i][j] < min)
	   min = ScoreBoard[i][j];
       }
     }
 
     tolerance = (max - min)*0.15;  // tolerance is 15% of difference between max and min values in 
     //tolerance= 0;                               // scoreboard

     // cout<<"\n debug Max , Min and tolerance values are "<<max<<"\t"<<min<<"\t"<<tolerance<<endl;   
     // end of section added by Naveen 
     // ------------------------------------------------------------------	


     for (register int obj=0; obj<numObjs; obj++)
     {
          vector<int> list; //ext+seq_wr+no_init
          if(sm==_MERGE_ && tag2.at(obj))
               continue;
          if(sm==_SPLIT_)
               NumNonZeros(obj, list, _ROW_,t1,t2,ScoreBoard,tag1,tag2);
          else 
               NumNonZeros(obj, list, _COL_,t1,t2,ScoreBoard,tag1,tag2);
          int nzero = list.size(); // nzero is actually the number of objects in the other frame overlapping with this object.
          //#cout<<funcname<<":obj"<<obj<<" has "<<nzero<<" overlapping objs in the other frame\n";
          if (nzero>1)
          {
               for (register int i=0; i<nzero; i++)
               {
	            vector<int>::iterator mytag1;
	            vector<int>::iterator mytag2;
	            mytag1= (sm==_SPLIT_) ? tag1.begin():tag2.begin();
	            mytag2= (sm==_SPLIT_) ? tag2.begin():tag1.begin();
	            if (*(mytag2+list.at(i))==0)
                    {
	                 vector<int> SMList; //ext+seq_wr+no_init
	                 int nSM;  
	                 if(sm==_SPLIT_)
	                      nSM=GetSameScore(obj, i, nzero, list, SMList,_ROW_,ScoreBoard,tag1,tag2, tolerance);
	                 else
	                      nSM=GetSameScore(obj, i, nzero, list, SMList,_COL_,ScoreBoard,tag1,tag2, tolerance);

	                 if (nSM>1)
                         { // means there is split/merge

	                     //# cout<<funcname<<":obj "<<obj<<" has split/merge with "<<nSM<<" objs in the other frames\n";
                              *(mytag1+obj) = 1;    /* tag the object in the list1 */ 
	                      ++nSplits;	
	                      if(sm==_SPLIT_)
                              {
	                            objs_bcheck(step-1,obj,objs);
	                            objs.at(step-1).at(obj).issplit = Consts::YES;
	                            fprintf(outfile,"%d\t-1\t", obj+1);
	                      }
	                      else
                              {
	                           objs_bcheck(step,obj,objs);
	                           objs.at(step).at(obj).ismerge=Consts::YES;
	                      }
	                      for(register int j=0; j<nSM; j++)
                              {
	                           int obj2(SMList.at(j));
	                           *(mytag2+obj2)=1;
	                           fprintf(outfile,"%d\t", obj2+1);
 	                           if(sm==_SPLIT_)
                                   {
                        	        objs_bcheck(step-1,obj,objs);
		                        objs.at(step-1).at(obj).children.push_back(obj2);
		                        objs_bcheck(step,obj2,objs);
		                        objs.at(step).at(obj2).parents.push_back(obj);
	                           }
                                   else
                                   {
		                        objs_bcheck(step-1,obj2,objs);
		                        objs.at(step-1).at(obj2).children.push_back(obj);
		                        objs_bcheck(step,obj,objs);
		                        objs.at(step).at(obj).parents.push_back(obj2);
	                           }
	                      }
                              if(sm==_SPLIT_)
                                  fprintf(outfile,"\n");
	                      else
                                   fprintf(outfile,"-1\t%d\n",obj+1);

	                      //# cout<<funcname<<":after handling obj"<<obj<<" who has split/merge\n";
                        }
	            }
               }
          }
     }
    //# cout<<"TrackSplit_Merge:";
     if(sm==_SPLIT_) 
          cout<<nSplits<<" splits\n";
     else 
          cout<<nSplits<<" merges\n";
}

// Modified by Weiping
// Changes: THe difference between every feature in the second frame with the previous feature should be lower than the tolerence
void TrackSplit_Merge2(int const step, SPLIT_MERGE const sm,Frame& t1,Frame &t2, FILE* outfile,vector<vector<TrackObject> >&  objs,vector<int>& tag1,vector<int>& tag2,vector<vector<int> >& ScoreBoard, double tolerance)
{
     int numObjs, nSplits(0), objVolume; // nSplits counts the times of splits/merges, for debug purpose.
     if(sm==_SPLIT_)
          numObjs=t1.objVols.size();
     else
          numObjs=t2.objVols.size();
     string funcname;
     if(sm==_SPLIT_)
          funcname="TrackSplit_Merge(SPLIT)";
     else
          funcname="TrackSplit_Merge(MERGE)";


     /*
     // section added by Naveen to find max and min values in
     // scoreboard-------------------------------------
     // finding the maximum and minimum values in the scoreboard and calculate tolerance which is some %
     // of difference between max and min
     // this tolerance is passed as an argument to GetSameScore function . the objects are considered to
     // be split or merged if objects in same row or column
     // are equal for lie between this tolerance level.

     // Modified by Weiping Hua
     // Min/Max are fixed to [0, 1000] in Score Computing
     // I doubt why we need the Get SameScore here
     // Thus I simplify the code to check if the Score is bigger than the score.
     */

     int max = __DBL_MIN__;
     int min = __DBL_MAX__;

//     for(int i = 0; i< ScoreBoard.size() ; i++)
//     {
//       for(int j = 0; j<ScoreBoard[i].size(); j++)
//       {
//         if(ScoreBoard[i][j] > max)
//            max = ScoreBoard[i][j];
//         if(ScoreBoard[i][j] < min)
//       min = ScoreBoard[i][j];
//       }
//     }

     //Defined by the ScoreBoard
     max=1000;
     min=0;
     tolerance = (max - min)*0.75;  // tolerance is 15% of difference between max and min values in
     //tolerance= 0;                               // scoreboard

     // cout<<"\n debug Max , Min and tolerance values are "<<max<<"\t"<<min<<"\t"<<tolerance<<endl;
     // end of section added by Naveen
     // ------------------------------------------------------------------


     for (int obj=0; obj<numObjs; obj++)
     {
          vector<int> list; //ext+seq_wr+no_init
          vector<int> SMList; //ext+seq_wr+no_init
          if(sm==_MERGE_ && tag2.at(obj))
               continue;
          if(sm==_SPLIT_){
               NumNonZeros(obj, list, _ROW_,t1,t2,ScoreBoard,tag1,tag2);
               objVolume=t1.objVols.at(obj).ObjVol;
          }
          else {
               NumNonZeros(obj, list, _COL_,t1,t2,ScoreBoard,tag1,tag2);
               objVolume=t2.objVols.at(obj).ObjVol;
          }
          int nzero = list.size(); // nzero is actually the number of objects in the other frame overlapping with this object.
          //#cout<<funcname<<":obj"<<obj<<" has "<<nzero<<" overlapping objs in the other frame\n";
          if (nzero>1)
          {
//               for (int i=0; i<nzero; i++)
//               {
                vector<int>::iterator mytag1;
                vector<int>::iterator mytag2;
                mytag1= (sm==_SPLIT_) ? tag1.begin():tag2.begin();
                mytag2= (sm==_SPLIT_) ? tag2.begin():tag1.begin();
//                if (*(mytag2+list.at(i))==0)
//                    {
//                     vector<int> SMList; //ext+seq_wr+no_init
                     int nSM;
                     if(sm==_SPLIT_)
//	                      nSM=GetSameScore(obj, i, nzero, list, SMList,_ROW_,ScoreBoard,tag1,tag2, tolerance);
                          nSM=CheckSM_Condition(obj, nzero, list, SMList,_ROW_,ScoreBoard,tag1,tag2, tolerance,objVolume);
                     else
//	                      nSM=GetSameScore(obj, i, nzero, list, SMList,_COL_,ScoreBoard,tag1,tag2, tolerance);
                          nSM=CheckSM_Condition(obj, nzero, list, SMList,_COL_,ScoreBoard,tag1,tag2, tolerance,objVolume);

                     if (nSM>1)
                         { // means there is split/merge

                         //# cout<<funcname<<":obj "<<obj<<" has split/merge with "<<nSM<<" objs in the other frames\n";
                              *(mytag1+obj) = 1;    /* tag the object in the list1 */
                          ++nSplits;
                          if(sm==_SPLIT_)
                              {
                                objs_bcheck(step-1,obj,objs);
                                objs.at(step-1).at(obj).issplit = Consts::YES;
                                fprintf(outfile,"%d\t-1\t", obj+1);
                          }
                          else
                              {
                               objs_bcheck(step,obj,objs);
                               objs.at(step).at(obj).ismerge=Consts::YES;
                          }
                          for(int j=0; j<nSM; j++)
                              {
                               int obj2(SMList.at(j));
                               *(mytag2+obj2)=1;
                               fprintf(outfile,"%d\t", obj2+1);
                               if(sm==_SPLIT_)
                                   {
                                    objs_bcheck(step-1,obj,objs);
                                objs.at(step-1).at(obj).children.push_back(obj2);
                                objs_bcheck(step,obj2,objs);
                                objs.at(step).at(obj2).parents.push_back(obj);
                               }
                                   else
                                   {
                                objs_bcheck(step-1,obj2,objs);
                                objs.at(step-1).at(obj2).children.push_back(obj);
                                objs_bcheck(step,obj,objs);
                                objs.at(step).at(obj).parents.push_back(obj2);
                               }
                          }
                              if(sm==_SPLIT_)
                                  fprintf(outfile,"\n");
                          else
                                   fprintf(outfile,"-1\t%d\n",obj+1);

                          //# cout<<funcname<<":after handling obj"<<obj<<" who has split/merge\n";
                        }
//                }
//               }
          }
     }
    //# cout<<"TrackSplit_Merge:";
     if(sm==_SPLIT_)
          cout<<nSplits<<" splits\n";
     else
          cout<<nSplits<<" merges\n";
}

/*! \fn void TrackContinuous(int step,Frame& t1,Frame &t2,FILE* outfile,vector<vector<Obj> >&  objs,vector<int>& tag1,vector<int>& tag2,vector<vector<int> >& ScoreBoard)
  \brief track continun events
 */
void TrackContinuous(int step,Frame& t1,Frame &t2,FILE* outfile,vector<vector<TrackObject> >&  objs,vector<int>& tag1,vector<int>& tag2,vector<vector<int> >& ScoreBoard, double& tolerance)
{
     int obj1, obj2, nCont(0); //nCont counts the times of continuns, for debug purpose.
    
     int numObjs1=t1.objVols.size(); 
     int numObjs2=t2.objVols.size(); 
     for (obj1=0; obj1<numObjs1; obj1++)
     {
          for (obj2=0; obj2<numObjs2; obj2++)
          {
           if (ScoreBoard.at(obj1).at(obj2)>250 && (!tag1.at(obj1)) && (!tag2.at(obj2)))
               {
	            if (IsMax(obj1, obj2,_ROW_,ScoreBoard,tag1,tag2,t1,t2, tolerance) && IsMax(obj1, obj2,_COL_,ScoreBoard,tag1,tag2,t1,t2, tolerance))
                    {
		         tag1.at(obj1)=1;
		         tag2.at(obj2)=1;
		         ++nCont;
		         fprintf(outfile,"%d\t-1\t%d\n", obj1+1,obj2+1);
 		         objs_bcheck(step-1,obj1,objs);
		         objs.at(step-1).at(obj1).children.push_back(obj2);
		         objs_bcheck(step,obj2,objs);
		         objs.at(step).at(obj2).parents.push_back(obj1);
		         objs.at(step).at(obj2).ismerge = Consts::NO;
		         objs.at(step).at(obj2).issplit = Consts::NO;
		         objs.at(step).at(obj2).isnew = Consts::NO;	  
		    }
	       }
	  }
     }
     //#cout<<"TrackContinuous:"<<nCont<<" Continuns\n";
}

/*! \fn TrackNew_Dissipate(int const step, NEW_DISSIPATE const nd,Frame& t1,Frame &t2,FILE* outfile,vector<vector<Obj> >&  objs,vector<int>& tag1,vector<int>& tag2)
  \brief track new or dissipating events
 */
void TrackNew_Dissipate(int const step, NEW_DISSIPATE const nd,Frame& t1,Frame &t2,FILE* outfile,vector<vector<TrackObject> >&  objs,vector<int>& tag1,vector<int>& tag2)
{
     int numObjs, nNew(0); //nNew counts the times of New/Dissipate, for debug purpose.
     if(nd==_NEW_) 
          numObjs=t2.objVols.size();
     else
          numObjs=t1.objVols.size();
     for (int obj=0; obj<numObjs; obj++)
     {
          vector<int>::iterator mytag;
          if(nd==_NEW_)
               mytag=tag2.begin();
          else 
               mytag=tag1.begin();
          if (!*(mytag+obj))
          {
               *(mytag+obj) = 1;
               ++nNew;
               if(nd==_NEW_)
               {
	            objs_bcheck(step,obj,objs);
	            objs.at(step).at(obj).isnew=Consts::YES;
	            fprintf(outfile,"-1\t%d\n",obj+1);
               }
               else
               {
	            objs_bcheck(step-1,obj,objs);
	            objs.at(step-1).at(obj).issplit=Consts::NO;
	            objs.at(step-1).at(obj).ismerge=Consts::NO;
	            fprintf(outfile,"%d\t-1\n",obj+1);
               }
          }
     }
     //#cout<<"TrackNew_Dissipate:";
     if(nd==_NEW_)
          cout<<nNew<<" News\n";
     else 
          cout<<nNew<<" Dissipates\n";
}
/*! \fn void NumNonZeros(int const row,vector<int> & list, ROW_COL const rc,Frame& t1,Frame &t2 ,vector<vector<int> >& ScoreBoard,vector<int>& tag1,vector<int>& tag2)
  \brief  Used in TrackSplit_Merge(). Get the nonzero elements in Scoreboard, put the nonzeror obj id(overlapping obj) into list.
 */
void NumNonZeros(int const row,vector<int> & list, ROW_COL const rc,Frame& t1,Frame &t2 ,vector<vector<int> >& ScoreBoard,vector<int>& tag1,vector<int>& tag2)
{
     int numObjs;
     if(rc==_ROW_)
          numObjs=t2.objVols.size();
    else
          numObjs=t1.objVols.size();
    for (register int i(0); i<numObjs; i++)
    {
         if(rc==_ROW_)
         {
              if ((ScoreBoard.at(row).at(i)) && (!tag2.at(i))) 
    	           list.push_back(i);
         }
         else if((ScoreBoard.at(i).at(row)) && (!tag1.at(i)))
    	      list.push_back(i);
     }
}

/*! \fn int GetSameScore(int const row,int const ind,int const nlist,vector<int> & list,vector<int> & split, ROW_COL const rc,vector<vector<int> >& ScoreBoard,vector<int>& tag1,vector<int>& tag2)
  \brief  Get the objects who have the same score as obj ind. It returns the number of such objects (including ind itself).
  \param row can be _ROW_ or _COL_
  \param ind the object index
  \param nlist the ???????????????????????????????????????
  list is an input parameter which has been known when this function is called. split is the output parameter where the objects found are stored.
nlist is an input which is the number of elements used to find the same scores.

*/
int GetSameScore(int const row,int const ind,int const nlist,vector<int> & list,vector<int> & split, ROW_COL const rc,vector<vector<int> >& ScoreBoard,vector<int>& tag1,vector<int>& tag2, double &tolerance)
{
     split.push_back(list.at(ind));
     for (register int i=ind+1; i<nlist; i++)
     {
          if(rc==_ROW_)
          {
               if (abs(ScoreBoard.at(row).at(list.at(i)) - ScoreBoard.at(row).at(list.at(ind))) < (int)tolerance && (tag2.at(list.at(i))==0)) 
                    split.push_back(list.at(i));        
          }
          else
          {
               if(abs(ScoreBoard.at(list.at(i)).at(row) - ScoreBoard.at(list.at(ind)).at(row)) < (int)tolerance && (tag1.at(list.at(i))==0)) 
	            split.push_back(list.at(i)); 
          }
     }
     return split.size();
}

// Modified by Weiping Hua
// Change the Split/Merge check
int CheckSM_Condition(int const row,int const nlist,vector<int> & list,vector<int> & split, ROW_COL const rc,vector<vector<int> >& ScoreBoard,vector<int>& tag1,vector<int>& tag2, double &tolerance, int objVolume)
{
    for (int i=0; i<nlist; i++)
    {
        if(rc==_ROW_)
        {
//            int test = ScoreBoard.at(row).at(list.at(ind));
           if (ScoreBoard.at(row).at(list.at(i)) > (int)tolerance && (tag2.at(list.at(i))==0))
                split.push_back(list.at(i));
        }
        else
        {
           if(ScoreBoard.at(list.at(i)).at(row) > (int)tolerance && (tag1.at(list.at(i))==0))
                split.push_back(list.at(i));
        }
    }
    return split.size();
}


/*! \fn bool IsMax(int const row,int const col, ROW_COL const rc,vector<vector<int> >& ScoreBoard,vector<int>& tag1,vector<int>& tag2,Frame& t1, Frame& t2)
  \brief see whether element <row,col> in scoreboard is the max in the row (with rc=_ROW_) or in the col (with rc=_COL_)
 */
bool IsMax(int const row,int const col, ROW_COL const rc,vector<vector<int> >& ScoreBoard,vector<int>& tag1,vector<int>& tag2,Frame& t1, Frame& t2, double &tolerance)
{
     int numObjs;
     if(rc==_ROW_)
          numObjs=t2.objVols.size();    
     else 
          numObjs=t1.objVols.size();
     for (register int i=0; i<numObjs; i++)
     {
          if(rc==_ROW_)
          {
               if (((ScoreBoard.at(row).at(i) - ScoreBoard.at(row).at(col)) > (int)tolerance) && (tag2.at(i)==0))
	            return false;
          }
          else
          {
              if (((ScoreBoard.at(i).at(col) - ScoreBoard.at(row).at(col)) > (int)tolerance) && (tag1.at(i)==0))
	            return false;
          }
     }
     return true;    
}
/*! \fn void GetInts(string const buffer, vector<int> & line)
  \brief convert a line including numbers into a number array and store them into a 1d array -- line.
 */
void GetInts(string const buffer, vector<int> & line)
{
 /* typedef boost::tokenizer<boost::char_separator<char> >tokenizer;
  boost::char_separator<char> sep("\t");
  tokenizer tok(buffer, sep);
  for(tokenizer::iterator it=tok.begin();it!=tok.end();++it)
    line.push_back(StoT<int>(*it));*/
     /* boost dependency is to be removed; so I have written this piece of code*/
//      char *pch;
//      char *buf = new char[buffer.length()];(t1.nodes[i].NodeID == t2.nodes[j].NodeID)
//      buffer.copy(buf,buffer.length(),0);
//      pch = strtok(buf,"\t");
//      while(pch!=NULL)
//      {
//           int k = atoi(pch);
//           if(k!=0)
//           {
//                line.push_back(k);
//           }
//           pch = strtok(NULL,"\t");
//      }
     string delimiters = "\t";
    vector<string> tokens;
    // Skip delimiters at beginning.
    string::size_type lastPos = buffer.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = buffer.find_first_of(delimiters, lastPos);
    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(buffer.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = buffer.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = buffer.find_first_of(delimiters, lastPos);
    }
    // now add the integers to line vector
    for(int i = 0; i < tokens.size();i++)
    {
         line.push_back(StoT<int>(tokens[i]));
    }
}
#endif
