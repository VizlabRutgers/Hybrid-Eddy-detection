#ifndef OBJECTSEGMENT_H
#define OBJECTSEGMENT_H
#include <limits.h>
//#include <vtkDataSet.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataArray.h>
#include <vtkContourFilter.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkIdList.h>
#include <vtkType.h>
#include <assert.h>
#include <string>
#include "input.h"
#include "nodeinfo.h"
#include "cellinfo.h"
#include "objectinfo.h"
#include "objSegmentUtil.h"
#include "classifier.h"
#define PI 3.14159265


using namespace std;
extern string getlabel(std::string listfile);
extern string precision_time(int const time,int const precision);
extern void SetCellNodes(unsigned long pos[8],stCellHex *cell);
extern int ProcessCells(InpObject<stCellHex> *inPtr);
extern void ResetArray(unsigned long nnodes, char *mergeArray);


/*
 extern string getlabel(std::string listfile);
 extern string precision_time(int const time,int const precision);
 extern void SetCellNodes(long pos[8],stCellHex *cell);
 extern int ProcessCells(InpObject<stCellHex> *inPtr);
 int Get_Conn(InpObject<stCellHex> *inPtr, int *node_conn_array);
 void ResetArray(int nnodes, char *mergeArray);
 void GetPoint(stCellHex *cell, long& point, int ptNumber);
 void SETARRAY(stCellHex *cell,char *mergeArray);
 void SetArray(int index, char *mergeArray);
 void HANDLE_CELL_IND(InpObject<stCellHex> *inPtr, stCellHex *cell, long* vertex_ind,long* speedingtable,char *mergeArray);
 void SetOutField(vtkUnstructuredGrid *out_ds, int nnodes, int ncells, int nspace,int cellpoints, float* coords, float* nodedata, int *connects);
 string getlabel(std::string listfile);
 int OutputAttribute(stObject *obj, char *outFile);
 int OutputTrak( stObject *obj, char *fileBaseName,int timeFrame);
 */
//extern int MergeCellList(InpObject<stCellHex> *inptr, stObject*objPtr, long objNum, stCellHex *cellList, char *mergeArray);

#define HANDLE_(a)     { if ( ((nodeList[a].flag & F_THR) == F_THR ) && \
((nodeList[a].flag & F_USED) != F_USED) )  {  \
cellPtr=GetIncidentCells( nodeList, a); \
AddCellListToObj(objPtr, cellPtr,mergeArray); \
nodeList[a].flag |= F_USED; vol++; }}

#define HANDLE_CELL_INDEX(a, b)  {	n_interval=(int)(b/INTERVAL); \
int t2, t1=n_interval*INTERVAL;  \
for(t2=0, i=t1;i<b;i++) if(mergeArray[i]) t2++;  \
cout<<"b=["<<b<<"] and t2 =["<<t2<<"]"<<endl;  \
if(n_interval) vertex_ind[a]=speedingtable[n_interval-1]+t2; \
else vertex_ind[a]=t2;  \
}

//extern void SetOutField(vtkUnstructuredGrid *out_ds, int nnodes, int ncells, int nspace,int cellpoints, float* coords, float* nodedata, int *connects);





string getlabel(std::string listfile)
{
    ifstream fp;
    fp.open(listfile.c_str());
    if (!fp.is_open())
    {
        cout<<"getlabel: cannot open listfile\n";
        string tmp="";
        return tmp;
    }
    string label;
    fp>>label;
    string file;
    fp>>file;
    fp.close();
    int s;
    for (s=file.size()-1;s>-1;--s) {
        if (file[s]>'9'||file[s]<'0')
            break;
    }
    label+=file.substr(0,s+1);
    return label;
}






//
//
//
//template<class T>
//int PopulateInputData(InpObject<T> *inPtr, unsigned long *node_conn_array, float *coord_array, float *node_data)
//{
//    //#cout << "Now in PopulateInputData \n";
//    register unsigned long i=0, j=0;
//    stNodePos    *nodePtr,   *curNode;
//    stNodeData   *nodeDataPtr,  *curDataPtr;
//    float max_node_data0=-99999.0, min_node_data0=99999.0;
//    float max_node_data1=-99999.0, min_node_data1=99999.0;
//    float max_node_data2=-99999.0, min_node_data2=99999.0;
//    int aaaa = Get_Conn(inPtr, node_conn_array);
//    //get coords
//    nodePtr= CreateNodeArray(inPtr->numNodes);
//    //assert(nodePtr);
//    cout<<"inPtr->numNodes: "<<inPtr->numNodes<<endl;
//    //#cout << "post creating node list\n";
//    for (j=0; j<inPtr->numNodes; j++)
//    {
//        i=3*j;
//        curNode = GetCurNode( nodePtr, j);
//        curNode->x=coord_array[i];
//        curNode->y=coord_array[i+1];
//        curNode->z=coord_array[i+2];
//        curNode->list=NULL;
//        curNode->flag=-1;
//    }
//    //get node_data
//    nodeDataPtr = CreateDataArray( inPtr->numNodes);
//    //assert(nodeDataPtr);
//    //#cout << "post creating data list\n";
//    
//    for (j=0; j<inPtr->numNodes; j++)
//    {
//        i=j*inPtr->nncomp;
//        curDataPtr = GetCurNodeData( nodeDataPtr, j);
//        curDataPtr->val0=node_data[i];
//        curDataPtr->val1=node_data[i+1];
//        curDataPtr->val2=node_data[i+2];
//        if(max_node_data0<node_data[i])
//            max_node_data0=node_data[i];
//        if(min_node_data0>node_data[i])
//            min_node_data0=node_data[i];
//        if(max_node_data1<node_data[i+1])
//            max_node_data1=node_data[i+1];
//        if(min_node_data1>node_data[i+1])
//            min_node_data1=node_data[i+1];
//        if(max_node_data2<node_data[i+2])
//            max_node_data2=node_data[i+2];
//        if(min_node_data2>node_data[i+2])
//            min_node_data2=node_data[i+2];
//    }
//    
//    
//    /*
//     for (j=0;j<inPtr->numNodes;j++)
//     {
//     curDataPtr = GetCurNodeData( nodeDataPtr, j);
//     curDataPtr->val0=node_data[j];
//     if(max_node_data<node_data[j])
//     max_node_data=node_data[j];
//     if(min_node_data>node_data[j])
//     min_node_data=node_data[j];
//     }
//     */
//    inPtr->pnode   = nodePtr;
//    inPtr->pnodeData  = nodeDataPtr;
//    inPtr->minData = min_node_data0;
//    inPtr->maxData = max_node_data0;
//    inPtr->minData1 = min_node_data1;
//    inPtr->maxData1 = max_node_data1;
//    inPtr->minData2 = min_node_data2;
//    inPtr->maxData2 = max_node_data2;
//    cout<<"Inside the populateInputdata: inPtr->thresh_deltaz: "<<inPtr->thresh_deltaz<<endl;
//    
//    return METHOD_SUCCESS;
//}
//




int Get_Conn(InpObject<stCellHex> *inPtr, unsigned long *node_conn_array)
{
    //#cout << "In Get_Conn\n";
    register long i = 0, j = 0;
    stCellHex *cellArrPtr, *curCell;
    cellArrPtr=CreateCellArray<stCellHex>(inPtr->numCell);
    // assert(cellArrPtr);
    if(cellArrPtr!=NULL)
        cout << "post cell array creation\n";
    else
        cout << "array null";
    for(j=0;j<inPtr->numCell;j++)
    {
        i=8*j;
        curCell = GetCurCell( cellArrPtr, j);
        curCell->point0=node_conn_array[i];
        curCell->point1=node_conn_array[i+1];
        curCell->point2=node_conn_array[i+2];
        curCell->point3=node_conn_array[i+3];
        curCell->point4=node_conn_array[i+4];
        curCell->point5=node_conn_array[i+5];
        curCell->point6=node_conn_array[i+6];
        curCell->point7=node_conn_array[i+7];
        curCell->objNum=-1; // the object No. to which this tetra belongs.
        curCell->flag = -1;
    }
    inPtr->pcell= cellArrPtr;
    return METHOD_SUCCESS;
}






//marking all nodes according to whether their node data is less than threshold or not.
template <class T>
int ProcessThresholds(InpObject<T> *inPtr, double threshold_value, int thrValType)
{
    cout << "inside processthreshold \n";
    cout<<"inPtr->thrVal:["<<inPtr->thrVal<<"] and inPtr->nncomp:["<<inPtr->nncomp<<"] "<<endl;
    
    register unsigned long j=0;
    stNodePos *nodeList, *curNodePtr;
    stNodeData *dataList, *curDataPtr;
    static unsigned long  threshCnt = 0;
    bool comparisonValue = 0;
    //int numUnused = 0; // counting the number of nodes whose node data is less than threshold
    threshCnt = 0;  // counting the number of nodes whose node data is greater than threshold
    //  cout<<"Inside ProcessThresholds!!: inPtr->nncomp"<<inPtr->nncomp<<endl;
    
    if(thrValType==THRESH_PERCENT)
    {
        inPtr->threshPercent_FromTop = threshold_value;
        inPtr->threshPercent_FromBottom = 0.0;
    }
    else
        cout << "cannot process threshold type other than PERCENTAGE\n";
    nodeList = inPtr->pnode;
    dataList= inPtr->pnodeData;
    for( j=0;j<inPtr->numNodes;j++ )
    {
        curDataPtr = GetCurNodeData(dataList, j );
        curNodePtr  = GetCurNode(nodeList, j);
        if(inPtr->thrType==THRESH_PERCENT)
        {
            if (inPtr->nncomp == 1)
            {
                comparisonValue = (curDataPtr->val0 >= inPtr->thrVal);
            }
            if (inPtr->nncomp == 3)
            {
                comparisonValue = ((curDataPtr->val0 >= inPtr->thrVal) && (curDataPtr->val1 >= inPtr->thrVal1));
            }
            if (comparisonValue)
            {
                curNodePtr->flag = F_THR;
                threshCnt++;
                //numUnused++;
                //curNodePtr->flag = F_UNUSED;
            }
            else
            {
                //curNodePtr->flag = F_THR;
                //threshCnt++;
                curNodePtr->flag = F_UNUSED;
                
            }
        }
        else
            cout << "cannot process threshold type other than PERCENTAGE\n";
    }
    cout<<"leaving processthreshold"<<endl;
    return METHOD_SUCCESS;
}


template<class T>
int SegmentObjects(InpObject<T> *inPtr, int minObjSize,char *mergeArray)
{
    //register int i =0;//, l=0;
    stCellIndex  *indexPtr = NULL, *nextPtr,  *index1Ptr, *cellPtr, *swapPtr;
    stObject *objPtr;
    stNodeData *dataList;
    stNodePos  *nodeList;
    T          *cellList, *cell;
    unsigned long pos[8],vol;
    
    dataList = inPtr->pnodeData;
    nodeList = inPtr->pnode;
    cellList = inPtr->pcell;
    
    indexPtr=NULL;
    //switch (inPtr->thrType)
    //{
    //case 0:    // upwards thresholding
    while(1)
    {
        //free indexPtr
        while(indexPtr)
        {
            nextPtr=indexPtr->next;
            delete indexPtr;
            indexPtr=nextPtr;
        }
        
        // find list of remaining unmarked nodes with greatest val
        indexPtr=GetMaxDataValue(inPtr);
        if (indexPtr == 0)  break; // stop if list is empty
      	
        // go through all nodes with val=maxVal in list
        while (indexPtr)
        {
            if((nodeList[indexPtr->index].flag & F_USED) == F_USED)
            {
                // if max node marked as used - advance to next one
                // used means this max node has been handled
                swapPtr = indexPtr->next;
                delete indexPtr;  // free temp mem held by node on max list
                indexPtr = swapPtr;
                continue;
            }
            
            ResetArray(inPtr->numCell,mergeArray);
            // before grow a new object, reset MergeArray.
            
            nodeList[indexPtr->index].flag |= F_USED; // mark node used
            objPtr = (stObject*)CreateObject();
            cellPtr = GetIncidentCells(nodeList, indexPtr->index );
            // get cells incident with this point.
            AddCellListToObj(objPtr, cellPtr,mergeArray);
            index1Ptr=objPtr->cellPtr;
            vol=0; // number of nodes(not cells) in the object
            while (index1Ptr != 0 )
            {
                cell = GetCurCell(cellList, index1Ptr->index);
                SetCellNodes(pos,cell);
                
                HANDLE_(pos[0]);
                HANDLE_(pos[1]);
                HANDLE_(pos[2]);
                HANDLE_(pos[3]);
                HANDLE_(pos[4]);
                HANDLE_(pos[5]);
                HANDLE_(pos[6]);
                HANDLE_(pos[7]);
                index1Ptr = index1Ptr->next;
            }
            
            
            AddObjToObjList( inPtr, objPtr, vol, minObjSize);
            index1Ptr=objPtr->cellPtr;
            while (index1Ptr)
            {
                cell = GetCurCell(cellList, index1Ptr->index);
                if (cell->objNum != INVALID_OBJ)
                    break;
                // problem: how about the new object is touching two or more objects???
                index1Ptr=index1Ptr->next;
            }
            if(index1Ptr == 0 )
                MarkCellList(cellList, objPtr);
            else
            {
                MergeCellList(inPtr, objPtr, cell->objNum, cellList,mergeArray);
                objPtr->cellPtr = 0x00;
                inPtr->pobject->objNum--;
            }
            swapPtr = indexPtr->next;
            delete indexPtr;  // free temp mem held by node on max list
            indexPtr = swapPtr;
        } //  while (indexPtr != 0 )
        
    }
    //   break;
    //}
    return METHOD_SUCCESS;
}

template <class T>
stCellIndex* GetMaxDataValue(InpObject<T> *inObj)
{
    unsigned long cnt;
    float val, maxVal;
    stCellIndex *tmpPtr, *maxretPtr;
    stNodeData  *dataList=inObj->pnodeData;
    stNodePos   *nodeList=inObj->pnode;
    
    maxretPtr=tmpPtr=0;
    cnt=inObj->numNodes;
    maxVal = LONG_MIN;
    register unsigned long i = 0;
    for (i=0;i<cnt;i++)
    {
        if ( nodeList[i].flag == F_UNUSED ) // means this cell is not in the threshold
            continue;
        if ( nodeList[i].flag & F_USED ) // this cell has been browsed
            continue;
        
        val = dataList[i].val0;
        maxVal = ( val > maxVal ) ? val : maxVal;
    }
    
    // SEDAT NOTES to SEDAT: try to find the min point (minVal) above as well.
    // And then instead of using the condition (maxVal == val ) below,
    // try to use (maxVal-val)/(maxVal - minVal)<0.01
    // also check and warn if (maxVal == minVal)
    
    /* make a list of all points with the exteme value */
    for ( i=0;i<cnt;i++ )
    {
        val = dataList[i].val0;
        if ( maxVal == val )
        {
            tmpPtr = new stCellIndex;
            //assert(tmpPtr);
            tmpPtr->index=i;
            tmpPtr->next = maxretPtr;
            maxretPtr = tmpPtr;
        }
    }
    return maxretPtr;
};
void ResetArray(unsigned long nnodes, char *mergeArray)
{
    
    for(register unsigned long i=0;i<nnodes;i++)
        mergeArray[i] = 0;
    
}

template<class T>
int SetObjMax( InpObject<T>  *inPtr,int cellpoints )
{
    stCellIndex   *cellPtr;
    stObject      *objPtr;
    stNodeData    *dataList;
    stNodePos     *nodeList;
    float          val;
    int            ptNumber;
    unsigned long           point;
    T              *cell, *cellList;
    int nspace = 3;
    
    float min, max, mass, massSqr; // setObjMax
    
    
    max = min = -100000.0;
    objPtr = inPtr->pobject;
    nodeList = inPtr->pnode;
    cellList = inPtr->pcell;
    dataList = inPtr->pnodeData;
    register unsigned long i = 0;
    for(i=0;i<inPtr->numNodes;i++)
        nodeList[i].flag &= (~F_USED);
    
    while( objPtr )
    {
        mass=0.0;
        massSqr=0.0;
        cellPtr = objPtr->cellPtr;
        objPtr->objMax = -100000.00;
        objPtr->objMin = 100000.00;
        objPtr->LowerLeft_Corner_x = 99999999999.99;
        objPtr->LowerLeft_Corner_y = 99999999999.99;
        objPtr->LowerLeft_Corner_z = 99999999999.99;
        objPtr->UpperRight_Corner_x = -999999999.99;
        objPtr->UpperRight_Corner_y = -999999999.99;
        objPtr->UpperRight_Corner_z = -999999999.99;
        
        
        while( cellPtr )
        {
            cell=GetCurCell(cellList, cellPtr->index);
            for(ptNumber=0; ptNumber<cellpoints; ptNumber++)
            {
                GetPoint(cell,point,ptNumber);
                // if point not visited already in object mark it and
                //  add its value to the object mass
                val=dataList[point].val0;
	            if ((((nodeList[point].flag) & F_USED)!=F_USED) &&
                    (((nodeList[point].flag) & F_THR)==F_THR))
                {
                    nodeList[point].flag|=F_USED;
                    mass += val;
                    massSqr += val*val;
                    
                    
                    // Find the box corners here for each objects
                    for(int j = 0; j <nspace;j++)
                    {
                        
                        if ( nodeList[point].x < objPtr->LowerLeft_Corner_x )
                            objPtr->LowerLeft_Corner_x = nodeList[point].x;
                        else if ( nodeList[point].x  > objPtr->UpperRight_Corner_x )
                            objPtr->UpperRight_Corner_x = nodeList[point].x;
                        
                        if ( nodeList[point].y < objPtr->LowerLeft_Corner_y )
                            objPtr->LowerLeft_Corner_y = nodeList[point].y;
                        else if ( nodeList[point].y  > objPtr->UpperRight_Corner_y )
                            objPtr->UpperRight_Corner_y = nodeList[point].y;
                        
                        
                        if ( nodeList[point].z < objPtr->LowerLeft_Corner_z )
                            objPtr->LowerLeft_Corner_z = nodeList[point].z;
                        else if ( nodeList[point].z  > objPtr->UpperRight_Corner_z )
                            objPtr->UpperRight_Corner_z = nodeList[point].z;
                        
                    }
                    
                    
                    
                    //
                    
                    if(val > objPtr->objMax )  // if value exceeds max update max
                    {
                        objPtr->objMax=val;
                        objPtr->maxX=nodeList[point].x;
                        objPtr->maxY=nodeList[point].y;
                        objPtr->maxZ=nodeList[point].z;
                        objPtr->maxNode=point;
                        
                    }
                    if(val < objPtr->objMin )  // if value exceeds min update min
                    {
                        objPtr->objMin=val;
                        objPtr->minX=nodeList[point].x;
                        objPtr->minY=nodeList[point].y;
                        objPtr->minZ=nodeList[point].z;
                        objPtr->minNode=point;
                    }
	            }
            }
            cellPtr = cellPtr->next;
        }
        
        objPtr->mass = mass;
        objPtr->massSqr = massSqr;
        
        objPtr = objPtr->next;
    }
    return METHOD_SUCCESS;
}

void GetPoint(stCellHex *cell, unsigned long& point, int ptNumber)
{
    switch (ptNumber)
    {
        case 0: point=cell->point0;  break;
        case 1: point=cell->point1;  break;
        case 2: point=cell->point2;  break;
        case 3: point=cell->point3;  break;
        case 4: point=cell->point4;  break;
        case 5: point=cell->point5;  break;
        case 6: point=cell->point6;  break;
        case 7: point=cell->point7;  break;
        default: break;
    }
}




void SETARRAY(stCellHex *cell,char *mergeArray)
{
  	SetArray(cell->point0,mergeArray);
   	SetArray(cell->point1,mergeArray);
   	SetArray(cell->point2,mergeArray);
   	SetArray(cell->point3,mergeArray);
   	SetArray(cell->point4,mergeArray);
   	SetArray(cell->point5,mergeArray);
   	SetArray(cell->point6,mergeArray);
   	SetArray(cell->point7,mergeArray);
    //sedat	cout<<"cell->point7 :"<<cell->point7<<endl;
	
}
void SetArray(unsigned long index, char *mergeArray)
{
    mergeArray[index] = 1;
}

void HANDLE_CELL_IND(InpObject<stCellHex> *inPtr, stCellHex *cell, unsigned long* vertex_ind,unsigned long* speedingtable,char *mergeArray)
{
    int n_interval=0, i = 0;
    HANDLE_CELL_INDEX(0,cell->point0);
    HANDLE_CELL_INDEX(1,cell->point1);
    HANDLE_CELL_INDEX(2,cell->point2);
    HANDLE_CELL_INDEX(3,cell->point3);
    HANDLE_CELL_INDEX(4,cell->point4);
    HANDLE_CELL_INDEX(5,cell->point5);
    HANDLE_CELL_INDEX(6,cell->point6);
    HANDLE_CELL_INDEX(7,cell->point7);
}


/************** Find object Centroid,volume & tetra num.************/
template <class T>
int FindCentAndVol( InpObject<T> *inPtr,int cellpoints)
{        /* this function can not be called prior to mass calc func. */
    stCellIndex   *cellPtr;
    stObject      *objPtr;
    stNodeData    *dataList;
    stNodePos     *nodeList;
    float         val=0.0, Cx, Cy, Cz;
    int           ptNumber=0,  cellNum=0;
    unsigned long          pointCount=0, point=0;  // pointCount is used to set volume
    T             *cell, *cellList;
    
    objPtr = inPtr->pobject;
    nodeList = inPtr->pnode;
    register unsigned long i = 0;
    for (i=0;i<inPtr->numNodes;i++)       // reset all point's Used flag
        nodeList[i].flag &= (~F_USED);
    nodeList = inPtr->pnode;
    cellList = inPtr->pcell;
    dataList = inPtr->pnodeData;
    while ( objPtr != NULL )
    {
        pointCount = 0;
        Cx=Cy=Cz=0.0;
        cellPtr = objPtr->cellPtr;
        cellNum=0;   // reset number of tetrahedra in object
        while(cellPtr)
        {
            cellNum++;
            cell=GetCurCell(cellList, cellPtr->index);
            for(ptNumber=0; ptNumber<cellpoints; ptNumber++)
            {
                GetPoint(cell,point,ptNumber);
                // if point not visited already in object mark it and count it for the object volume
                val=dataList[point].val0;
                if ((((nodeList[point].flag) & F_USED)!=F_USED) &&
                    (((nodeList[point].flag) & F_THR)==F_THR))
                {
                    nodeList[point].flag|=F_USED;
	                pointCount++;
                    Cx += (val*nodeList[point].x)/((float)(objPtr->mass));
                    Cy += (val*nodeList[point].y)/((float)(objPtr->mass));
	                Cz += (val*nodeList[point].z)/((float)(objPtr->mass));
                }
            }
            cellPtr = cellPtr->next;
        }
        objPtr->numCell = cellNum;
        objPtr->volume = pointCount;
        objPtr->centroidX = Cx;
        objPtr->centroidY = Cy;
        objPtr->centroidZ = Cz;
        objPtr = objPtr->next;
    }
    return 1;
}



// SetOutField(out_ds,cursor,objPtr->numCell,nspace,cellpoints,coords,nodevals,connects);

void SetOutField(vtkUnstructuredGrid *out_ds, unsigned long nnodes, unsigned long ncells, int nspace,int cellpoints, float* coords, double* nodedata, unsigned long *connects)
{
    //watch<float>(coords, 20, _COORDS_);
    //#cout << "In SetOutField \n";
    vtkFloatArray *pcoords = vtkFloatArray::New();
    pcoords->SetNumberOfComponents(nspace);
    pcoords->SetNumberOfTuples(nnodes);
    
    unsigned long n_connect = 0;
    unsigned long n_coords = 0;
    float *temp_coord = new float[3];
    //float temp_coord[3] = {0.0};
    //cout<<" ---------Inside setoutfield.. nnodes["<<nnodes<<"] ncells=["<< ncells  <<"]------------------ "<<endl;
    
    
    for(unsigned long i = 0; i < /*ncell*cellpoints*/ nnodes; i ++)
    {
        for(int ii = 0 ; ii < nspace;ii++)
        {
            //cout<<" Inside setoutfield.. coords["<<n_coords<<"]= "<< coords[n_coords]  <<endl;
            temp_coord[ii] = coords[n_coords++];

        }
        
        //cout<<"settuple3: i=["<<i<<"] temp_coord[0]=["<<temp_coord[0]  <<"] temp_coord[1]=["<<temp_coord[1]<<"] temp_coord[2]=["<<temp_coord[2]<<"]"<<endl;
        pcoords->SetTuple3(i,(double)temp_coord[0],(double)temp_coord[1],(double)temp_coord[2]);
        // pcoords->SetTuple(i,temp_coord);
        
    }
    
    out_ds->Allocate(ncells*cellpoints);// this step is important
    vtkPoints *points = vtkPoints::New();
    //points->SetNumberOfPoints(nnodes);
    
    points->SetData(pcoords);
    // ### just to test the SetPoint function++
    
    // ### just to test the SetPoint function --
    //#cout << "post writing coords to vtk\n";
    
    //     out_ds->GetPointData()->SetNumberOfComponents(1);
    // vtkPointData *vtk_pt_data = vtkPointData::New();
    //vtk_pt_data->GetScalars()->SetNumberOfTuples(nnodes);
    vtkFloatArray *pdata = vtkFloatArray::New();
    pdata->SetNumberOfTuples((vtkIdType)nnodes);
    //cout << " Num Tuples in data : " <<  vtk_pt_data->GetScalars()->GetNumberOfTuples() << endl;
    
    //vtkDataArray *vtk_node_data = vtk_pt_data->GetScalars();
    double *temp_data = new double;
    for(unsigned long i = 0; i < nnodes; i++)
    {
        //cout<<"inside the for loop. nodedata["<<i<<"]="<<nodedata[i]<<endl;
        *temp_data = nodedata[i];
        pdata->SetTuple1(i,(double)*temp_data);
    }
    //#cout << "post writing node_data to vtk\n";
    // int *temp_connect = new int[cellpoints];
    
    
    vtkIdType *temp_connect = new vtkIdType[cellpoints];     //int n_connect = 0;
    
    //cout<<"---------Before InsertNextCell --------- "<<endl;
    
    for(unsigned long i = 0; i < ncells;i++)
    {
        for(int j = 0; j < cellpoints;j++)
        {
            temp_connect[j] = connects[n_connect++];
            //cout<<"temp_connect["<<j<<"]="<<temp_connect[j]<<endl;
        }
        out_ds->InsertNextCell(VTK_VOXEL, (vtkIdType) cellpoints,temp_connect);
    }
   // cout<<"---------AFter InsertNextCell --------- "<<endl;
    
    out_ds->SetPoints(points);
    points->Delete();
    pcoords->Delete();
    out_ds->GetPointData()->SetScalars(pdata);
    pdata->Delete();
    //cout << "End of SetOutField\n";
}

/*
template<class T>
int ExtractSurf(vtkDataSet *in_ds,InpObject<T> *inPtr,vtkDataSet **outDS, char* mypolyfile,int nspace,int cellpoints)
{
   // cout<<"------  Inside ExtractSurf -------  [ 1 ]   ------- "<<endl;
    
    T             *cell,*cellList;
    stCellIndex   *cellPtr;
    stNodePos     *nodeList;
    stObject      *objPtr;
    stNodeData    *dataList;
    objPtr = inPtr->pobject;
    cellList = inPtr->pcell;
    nodeList = inPtr->pnode;
    dataList = inPtr->pnodeData;
    register unsigned long i =0,j=0;
    char *mergeArray;
    for (i=0;i<inPtr->numNodes;i++)
        nodeList[i].flag &= (~F_USED);
    
    //
    
    
    // Now I use mergeArray to process the coordinates
    // of the object.
    mergeArray=new char[inPtr->numNodes];
    assert(mergeArray);
    //#cout << "Post creating merge\n";
    bool firstobj=true;
    long *min_node_list, *max_node_list;
    long maxVol=0,minVol=9999999,minVolobj,maxVolobj,minAreaobj,maxAreaobj;
    double maxAr=0.0,minAr=9999999.99;
    ofstream afile;
    
    cout<<"\nmypolyfile has value"<<mypolyfile;
     std::string temp_file(mypolyfile);
     int dpos=temp_file.find('.');
     temp_file=temp_file.substr(0,dpos+1);
     cout<<"\ntemp_file has value"<<temp_file;
     temp_file.append("addl");
   // cout<<" objPtr value[ "<<objPtr<<" ] "<<endl;
    
    while (objPtr)
    {
        // build a field acceptable by UTILisosurf() in avs/dv_util.h
        cellPtr = objPtr->cellPtr;
        ResetArray(inPtr->numNodes,mergeArray);
        while(cellPtr)  //for each object.. set all the cell nodes =1 in mergeArray array.
        {
            cell=GetCurCell(cellList,cellPtr->index);
            SETARRAY(cell,mergeArray);
            cellPtr=cellPtr->next;
        }
        
        
        unsigned long *speedingtable;
        speedingtable=new unsigned long[(unsigned long)(inPtr->numNodes/INTERVAL)];
        assert(speedingtable);
        //#cout << "Post creating speeding table\n";
        unsigned long cursor=0;
        for(i=0;i<inPtr->numNodes;i++)
        {
            if(mergeArray[i])
                cursor++;
            if(i>0 && (i+1)%INTERVAL==0)
	            speedingtable[(unsigned long)((i+1)/INTERVAL)-1]=cursor;
        }
        
        float *coords=0;
        double *nodevals=0;
        coords=new float[cursor*nspace];
        assert(coords);
        nodevals=new double[cursor];
        assert(nodevals);
        for(i=0, j=0;i<inPtr->numNodes;i++) // set coords and node values
            if(mergeArray[i])
            {
                coords[j++]=nodeList[i].x;
                coords[j++]=nodeList[i].y;
                coords[j++]=nodeList[i].z;
                nodevals[(unsigned long)(j/nspace)-1]=dataList[i].val0;
            }
        //cout << "------ cursor value = ["<< cursor<<"\n";
        
        unsigned long *connects=0;
        //#cout << "objPtr->numCell " << objPtr->numCell <<endl;
        //#cout << "cellpoints " << cellpoints <<endl;
        connects=new unsigned long[objPtr->numCell*cellpoints];
        assert(connects);
        cellPtr=objPtr->cellPtr;
        j=0;
        unsigned long *vertex_ind;
        vertex_ind=new unsigned long[cellpoints];
        assert(vertex_ind);
        while(cellPtr)
        { //set connects
            cell=GetCurCell(cellList,cellPtr->index);
            HANDLE_CELL_IND(inPtr,cell,vertex_ind, speedingtable,mergeArray);
            for(i=0;i<cellpoints;i++)
            {
                connects[j++]=vertex_ind[i];
                
            }
            
            cellPtr=cellPtr->next;
        }
        delete [] vertex_ind;
        vertex_ind = NULL;
        //#cout << "Post filling the connects\n";
        vtkUnstructuredGrid *out_ds = vtkUnstructuredGrid::New();
        
        
        
        
        
        SetOutField(out_ds,cursor,objPtr->numCell,nspace,cellpoints,coords,nodevals,connects);
        //#cout << "GetNumberOfCells() : " << out_ds->GetNumberOfCells() << endl;
        // now call the isosurface module
        
       // cout<<"------  Inside ExtractSurf -------  [ 2 ]   ------- "<<endl;
        
        vtkContourFilter *isosurface = vtkContourFilter::New();
        isosurface->SetInputConnection( in_ds );
        isosurface->SetNumberOfContours(1);
        //SEDAT look at here, thresholds.....
        isosurface->SetValue(0,inPtr->thrVal); // from Sedat: some web on internet says that this is double type sometimes it says this is float. which one? :). Currently it is double.
        //      cout<<" isosurface->SetValue(0,inPtr->thrVal); : "<<inPtr->thrVal<<endl;
       // cout<<"------  Inside ExtractSurf -------  [ 3]   ------- "<<endl;
        
        //#cout << "isosurface->GetValue() "  << isosurface->GetValue(0) << endl;
        //*outDS = (vtkDataSet *)isosurface->GetOutput();
        vtkPolyData *outpoly = isosurface->GetOutput();
       // cout<<"------  Inside ExtractSurf -------  [ 3.5 ]   ------- "<<endl;
        //cout<<"test 1"<<endl;
        outpoly->Update();
       // cout<<"------  Inside ExtractSurf -------  [ 4 ]   ------- "<<endl;
      //  cout<<"outpoly->GetNumberOfPoints():["<<outpoly->GetNumberOfPoints()<<"]"<<endl;
        //cout<<"test 2"<<endl;
        outpoly->Register(NULL);
        // cout<<"test 3"<<endl;
        *outDS = outpoly;
        //  cout<<"test 4"<<endl;
        out_ds->Delete();
        out_ds = NULL;
        isosurface->Delete();
        isosurface= NULL;
        double p = 0;
        //#cout << "the poly data is : \n";
        ofstream pfile;
        if(firstobj)
            pfile.open(mypolyfile,ofstream::out | ofstream::trunc);
        else
            pfile.open(mypolyfile,ofstream::out | ofstream::app);
        if(!pfile.is_open())
        {
            cout<<"ExtractSurf:cannot open "<<mypolyfile<<endl;
            return 0;
        }
       // cout<<"------  Inside ExtractSurf -------  [ 5 ]   ------- "<<endl;
        
        //#cout << "Num cells post contour is : " << outpoly->GetNumberOfCells()<< endl;
        vtkIdList *listconnect = vtkIdList::New();
        
        //#cout << "Num points is : " << outpoly->GetNumberOfPoints()<< endl;
        //#cout << "writing to poly file\n";
        //#cout << "mypolyfie is : " << mypolyfile << endl;
        
        pfile<<objPtr->color.r<<" "<<objPtr->color.g<<" "<<objPtr->color.b<<endl;
        pfile<<outpoly->GetNumberOfPoints()<<endl;
        //Sedat      cout<<" outpoly->GetNumberOfPoints() : "<<outpoly->GetNumberOfPoints()<<endl;
        pfile.setf(ios::fixed,ios::floatfield);
        pfile<<setprecision(6);
        double *temp_coord_array = new double[3];
        for(i=0;i<outpoly->GetNumberOfPoints();i++)
        {
            temp_coord_array = outpoly->GetPoint(i);
            for(int j = 0; j <nspace;j++)
            {
                pfile<<temp_coord_array[j] << " ";
            }
            pfile<<endl;
        }
       // cout<<"------  Inside ExtractSurf -------  [ 6]   ------- "<<endl;
        
        pfile<<outpoly->GetNumberOfCells()<<endl;
        for(i=0;i<outpoly->GetNumberOfCells();i++)
        {
            outpoly->GetCellPoints(i,listconnect);
            pfile<<"3 ";
            for(j = 0; j < listconnect->GetNumberOfIds();j++)
            {
                unsigned long temp_id = (unsigned long)listconnect->GetId((const vtkIdType) j) ;
                pfile<< temp_id << " ";
            }
            pfile << endl;
        }
        pfile<<0<<endl<<endl;
        //#cout << "finished writing poly file\n";
        pfile.close();
        listconnect->Delete();
        listconnect = NULL;
        
        delete [] coords;
        coords = NULL;
        delete [] connects;
        connects = NULL;
        delete [] nodevals;
        nodevals = NULL;
        delete [] speedingtable;
        speedingtable = NULL;
        objPtr = objPtr->next;
        firstobj=false;
        //return 1;
        
    }
  //  cout<<"------  Inside ExtractSurf -------  [ 7 ]   ------- "<<endl;
    
    delete [] mergeArray;
    //#cout <<"End of Extract Surface\n";
    return METHOD_SUCCESS;
}

*/

template <class T>
int FindMoment( InpObject<T> *inPtr,int cellpoints )
{        // this function can not be called prior to Centroid calc func.
    stCellIndex   *cellPtr;
    stObject      *objPtr;
    stNodeData    *dataList;
    stNodePos     *nodeList;
    float         val =0.0, Ixx=0.0, Iyy=0.0, Izz=0.0, Iyz=0.0, Ixy=0.0, Izx=0.0;
    float         mx = 0.0, my = 0.0, mz = 0.0;
    unsigned long           /*rc =0, */ptNumber=0, objNum=0;
    register unsigned long i =0, point =0;
    T             *cell, *cellList;
    
    objPtr = inPtr->pobject;
    nodeList = inPtr->pnode;
    for (i=0;i<inPtr->numNodes;i++)       // reset all point's Used flag
        nodeList[i].flag &= (~F_USED);
    nodeList = inPtr->pnode;
    cellList = inPtr->pcell;
    dataList = inPtr->pnodeData;
    while( objPtr )
    {
        Ixx = Iyy = Izz = Ixy = Iyz = Izx = 0.0;
        cellPtr = objPtr->cellPtr;
        while(cellPtr)
        {
            cell = GetCurCell(cellList, cellPtr->index);
            for(ptNumber=0; ptNumber<cellpoints; ptNumber++)
            {
                GetPoint(cell,point,ptNumber);
                // if point not visited already in object mark it and
                // sum its weight into object moments
	            val=dataList[point].val0;
	            if ((((nodeList[point].flag) & F_USED)!=F_USED) &&
                    (((nodeList[point].flag) & F_THR)==F_THR))
                {
                    nodeList[point].flag|=F_USED;
                    mx = (nodeList[point].x-objPtr->centroidX);
                    my = (nodeList[point].y-objPtr->centroidY);
                    mz = (nodeList[point].z-objPtr->centroidZ);
                    Ixx += (val*pow(mx,2)/((float)(objPtr->mass)));
                    Iyy += (val*pow(my,2)/((float)(objPtr->mass)));
                    Izz += (val*pow(mz,2)/((float)(objPtr->mass)));
                    Iyz += (val*my*mz/((float)(objPtr->mass)));
                    Ixy += (val*mx*my/((float)(objPtr->mass)));
                    Izx += (val*mz*mx/((float)(objPtr->mass)));
                }
            }
            cellPtr = cellPtr->next;
        }
        objPtr->Ixx = Ixx;
        objPtr->Iyy = Iyy;
        objPtr->Izz = Izz;
        objPtr->Iyz = Iyz;
        objPtr->Ixy = Ixy;
        objPtr->Izx = Izx;
        objPtr = objPtr->next;
    }
    //rc=1;
    return METHOD_SUCCESS;
}

int OutputAttribute(stObject *obj, char *outFile)
{
    FILE *fpout;
    //int rc;
    //char buffer[256];
    stObject *objPtr = obj;
    char outAttr[256];
    char extention[6]=".attr";
    strcpy(outAttr,outFile);
    strcat(outAttr,extention);
    cout<<"outAttr is: "<<outAttr<<"!!!!!!!!"<<endl;
    if((fpout = fopen(outAttr, "w"))==NULL)
        cout << "cannot open outAttr file to write\n";
    
    while( objPtr != NULL )
    {
        cout<<"INSIDE THE WHILE LOOP"<<endl;
        fprintf(fpout,"------------------------------------------------\n");
        fprintf(fpout,"object %ld attributes:\n", objPtr->objNum);
        fprintf(fpout,"Max position: (%f, %f, %f) with value: %f\n",
                objPtr->maxX,objPtr->maxY,objPtr->maxZ,objPtr->objMax);
        //BySedat---
        cout<<"X, Y, Z and MaxValue: "<< objPtr->maxX<<" "<<objPtr->maxY<<" "<<objPtr->maxZ<<" "<<objPtr->objMax<<endl;
        //End of by sedat.... Just debugging purpose
        fprintf(fpout,"Node #: %ld\n",objPtr->maxNode);
        fprintf(fpout,"Min position: (%f, %f, %f) with value: %f\n",
                objPtr->minX,objPtr->minY,objPtr->minZ,objPtr->objMin);
        fprintf(fpout,"Node #: %ld\n",objPtr->minNode);
        fprintf(fpout,"Integrated content: %f\n",objPtr->mass);
        fprintf(fpout,"Sum of squared content values: %f\n",objPtr->massSqr);
        fprintf(fpout,"Volume: %ld\n",objPtr->volume);
        fprintf(fpout,"Centroid: (%f, %f, %f)\n",objPtr->centroidX,
                objPtr->centroidY, objPtr->centroidZ);
        fprintf(fpout,"Moment: Ixx = %f\nIyy = %f\nIzz = %f\n",
                objPtr->Ixx, objPtr->Iyy, objPtr->Izz);
        fprintf(fpout,"Ixy = %f\nIyz = %f\nIzx = %f\n",
                objPtr->Ixy, objPtr->Iyz, objPtr->Izx);
        objPtr = objPtr->next;
    }
    fclose(fpout);
    return METHOD_SUCCESS;
}


/********************** Output track file ***********************/
int OutputTrak( stObject *obj, char *fileBaseName,int timeFrame)
{
    FILE *fpout;
    stObject *objPtr;
    char extention[6]=".trak";
    char outTrak[256];
    float time;
    time=(float)(timeFrame);
    strcpy(outTrak,fileBaseName);
    strcat(outTrak,extention);
    objPtr = obj;
    if((fpout = fopen(outTrak, "w"))==NULL)
        cout << "cannot open outTrak file to write\n";
    while( objPtr != NULL )
    {
        fprintf(fpout,"%s  ",fileBaseName);
        fprintf(fpout,"%f  ",time);
        fprintf(fpout,"%14.6f   %5.0ld   %10.6f  %10.6f  %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %5.0ld\n",
                objPtr->mass,objPtr->volume,objPtr->centroidX,
                objPtr->centroidY, objPtr->centroidZ,objPtr->LowerLeft_Corner_x ,objPtr->LowerLeft_Corner_y ,objPtr->LowerLeft_Corner_z ,objPtr->UpperRight_Corner_x ,objPtr->UpperRight_Corner_y , objPtr->UpperRight_Corner_z, objPtr->Packet_ID);
        objPtr = objPtr->next;
    }
    fclose(fpout);
    return METHOD_SUCCESS;
}

/********************** Output uocd file *************************/
template <class T>
int OutputOcd( InpObject<T> *inPtr, char *outFile,int timeFrame,bool dovolren,int cellpoints/*string volfile,*//*OMobj_id in_field*/)
{ // outFile is the basename
    FILE           *fpout;
    stObject       *objPtr;
    char           extention[6]=".uocd";
    char           outOcd[256];
    float          time = 0.0, val = 0.0;
    unsigned long           objCount = 0;
    stCellIndex    *cellPtr;
    stNodeData     *dataList;
    stNodePos      *nodeList;
    int            ptNumber = 0;
    unsigned long           i = 0, point =0;
    T              *cell, *cellList;
    
    //volume rendering -- part 1
    /*int type,size;
     char *tmp_node_data; // 3. TO BE FREED WITH ARRFREE()
     if(FLDget_node_data(in_field,0,&type,&tmp_node_data,&size,OM_GET_ARRAY_RD)!=1)
     prt_rtn1("readavsdata:cannot get node_data\n");
     // type=DTYPE_FLOAT,size=140464
     
     
     float min, max, tempmin, tempmax;
     if(FLDget_node_data_minmax(in_field, 0,(char *)&tempmin,(char *)&tempmax)!= 1 )
     prt_rtn1("readavsdata:Failed to get field minmax values.\n");
     
     UTILtype_to_float(&min, (char *)&tempmin, type);
     UTILtype_to_float(&max, (char *)&tempmax, type);
     float delta=255/(max-min);
     
     #if _DEBUG_LEVEL_>0
     cout<<"readavsdata:min="<<min<<" max="<<max<<endl;
     #endif
     
     char *vol;
     FILE *volfp;
     int in_field_ndim=0, *in_field_dims,nnodes;
     if(dovolren) {
     cout<<"output volume file...\n";
     //get coordinate values
     if(FLDget_nnodes(in_field,&nnodes)!=1){
     //NOTE: in_field_id instead of cell_set!!!
     prt_rtn1("cannot get nnodes\n");}
     
     
     if(FLDget_dims(in_field, &in_field_dims, &in_field_ndim)!=1)
     prt_rtn1("in_field_dims not found\n");
     
     vol=(char*)malloc(nnodes*sizeof(char));
     if(!vol) prt_rtn1("readavsdata: allocate error\n");
     for(register int i=0;i<nnodes;++i) vol[i]=0;
     
     volfp=fopen(volfile.c_str(),"wb");
     if(!volfp) {
     cout<<"readavsdata:cannot open "<<volfile<<" to write\n";
     return false;
     }
     }*/
    //--- end of part 1
    
    time=(float)(timeFrame);
    strcpy(outOcd,outFile);
    strcat(outOcd,extention);
    objPtr = inPtr->pobject;
    if((fpout = fopen(outOcd, "w"))==NULL)
        cout << "cannot open outOcd file to write\n";
    nodeList = inPtr->pnode;
    for (i=0;i<inPtr->numNodes;i++)       // reset all point's Used flag
        nodeList[i].flag &= (~F_USED);
    nodeList = inPtr->pnode;
    cellList = inPtr->pcell;
    dataList = inPtr->pnodeData;
    fprintf(fpout,"%f\n",time);
    fprintf(fpout,"%ld\n", inPtr->numObj);
    while ( objPtr != NULL )
    {
        fprintf(fpout,"%ld\n", objPtr->objNum);
        fprintf(fpout,"%ld %f %f %f %f\n",objPtr->volume,
                objPtr->mass,objPtr->centroidX,objPtr->centroidY,
                objPtr->centroidZ);
        fprintf(fpout,"%f %f %f %f %f %f\n",objPtr->Ixx,
                objPtr->Iyy,objPtr->Izz,objPtr->Ixy,objPtr->Iyz,
                objPtr->Izx);
        cellPtr = objPtr->cellPtr;
        
        while ( cellPtr  != NULL )
        {
            cell=GetCurCell(cellList, cellPtr->index);
            for(ptNumber=0; ptNumber<cellpoints; ptNumber++)
            {
	            GetPoint(cell,point,ptNumber);
                // if point not visited already in object mark it and print it
	            val=dataList[point].val0;
	            if ((((nodeList[point].flag) & F_USED)!=F_USED) &&
                    (((nodeList[point].flag) & F_THR)==F_THR))
                {
                    nodeList[point].flag|=F_USED;
                    fprintf(fpout,"%6ld %9.6f %9.6f %9.6f  %f\n",
                            point, nodeList[point].x,nodeList[point].y,
                            nodeList[point].z,val);
                    /*if(dovolren)
                     {
                     float  node_data;    // 4.TO BE FREED with free()
                     char *node_data_elem;
                     // int xy_slice_size=in_field_dims[0]*in_field_dims[1];
                     // int index=(int)(nodeList[point].x+nodeList[point].y*in_field_dims[0]+nodeList[point].z*xy_slice_size);
                     // I assume the dims is xdimXydimXzdim
                     // when node data is stored, the change on z axis is the slowest, change on x axis is the fastest.
                     
                     // I need not compute the index because point is just the index I want.
                     node_data_elem=tmp_node_data+point*type; // tmp_node_data is also *char
                     // type=4 for float
                     assert(UTILtype_to_float(&node_data, node_data_elem, type)==1);
                     //int UTILtype_to_float( float *dst, const void *src, int dtype );
                     vol[point]=(char)((node_data-min)*delta);
	                 }*/
	            }
            }
            cellPtr = cellPtr->next;
        }
        objPtr = objPtr->next;
    }
    fclose(fpout);
    
    //volume rendering --part2
    /*if(dovolren)
     {
     ARRFREE(tmp_node_data); // 3 is ARRFREED
     //debug
     if(fwrite(vol,sizeof(char),nnodes,volfp)!=nnodes)
     prt_rtn1("readavsdata:cannot create vol file\n");
     fclose(volfp);
     string vol_templ_file;
     vol_templ_file=volfile.substr(0,volfile.rfind('_')+1);
     vol_templ_file+="templ.vol";
     volfp=0;
     volfp=fopen(vol_templ_file.c_str(),"wb");
     assert(volfp);
     for(register int i=0;i<nnodes;++i)
     if(vol[i]!=0)
     vol[i]=50;
     if(fwrite(vol,sizeof(char),nnodes,volfp)!=nnodes)
     prt_rtn1("readavsdata:cannot create vol file\n");
     fclose(volfp);
     free(vol);
     }*/
    // --end of part 2
    return METHOD_SUCCESS;
}

template <class T>
int SetColVol(InpObject<T> *inPtr, float *absMin, float *absMax,int rangeType)
{
    //return 0;
    //#cout << "Inside SetColVol \n";
    stObject *objPtr, *firstObj;
    long volume = 0, maxVol = 0, minVol = 0;
    float scalar = 0.0;
    objPtr = inPtr->pobject;
    firstObj = objPtr;
    maxVol = 0;
    minVol = 9999999;
    
    if (rangeType==ADAPTIVE)
        // if colormap is to strech adaptively to current min/max values
    {
        //#cout << "colormap type adaptive\n";
        while (objPtr != NULL) // find max volume of objects for normalization
        {
            volume=objPtr->volume;
            if (maxVol<volume)
                maxVol = volume;
            if (minVol>volume)
                minVol = volume;
            objPtr = objPtr->next;
        }
    }
    //#cout << "min value is : " << minVol << endl;
    //#cout << "max value is : " << maxVol << endl;
    if (rangeType==FIXED)  /* Range of values in colormap is user defined */
    {
        maxVol=(long)*absMax;
        minVol=(long)*absMin;
    }
    objPtr = firstObj;
    while (objPtr != NULL)
    {
        volume=objPtr->volume;
        /* normilize object volume over range */
        if(minVol!=maxVol)
        {
            scalar=(float)(volume-minVol)/(float)(maxVol-minVol);
            //#cout << "sclar is : " << scalar << endl;
        }
        else
            scalar=1.0;
        if(scalar > 1.0)        // object above colormap fixed top value
        {
            objPtr->color.r = 255;
            objPtr->color.g = 0;
            objPtr->color.b = 0;
        }
        else if(scalar < 0.0)   // object under colormap fixed bottom value
        {
            objPtr->color.r = 0;
            objPtr->color.g = 0;
            objPtr->color.b = 255;
        }
        else
            if(scalar <= 0.5)
            {
                objPtr->color.r = 0;
                objPtr->color.g = (int)((2.0 * scalar)*255);
                objPtr->color.b = (int)((1.0 - 2.0*scalar)*255);
            }
            else if(scalar <= 1.0)
            {
                objPtr->color.r = (int)((2.0 * (scalar - 0.5))*255);
                objPtr->color.g = (int)((1.0 - 2.0*(scalar - 0.5))*255);
                objPtr->color.b = 0;
            }
        objPtr = objPtr->next;
    }
    
    return(0);
}

template <class T>
int SetColMass(InpObject<T> *inPtr, float *absMin, float *absMax,int rangeType)
{
    
    //return 0;
    stObject *objPtr, *firstObj;
    float mass = 0.0, maxMass=0.0, minMass=0.0;
    float scalar = 0.0;
    
    objPtr = inPtr->pobject;
    firstObj = objPtr;
    maxMass = 0.0;
    minMass = 9999999.0;
    if (rangeType==ADAPTIVE)
        // if colormap is to strech adaptively to current min/max values
    {
        while (objPtr != NULL) // find max volume of objects for normalization
        {
            mass=objPtr->mass;
            if (maxMass<mass)
                maxMass = mass;
            if (minMass>mass)
                minMass = mass;
            objPtr = objPtr->next;
        }
    }
    if (rangeType==FIXED)  /* Range of values in colormap is user defined */
    {
        maxMass=*absMax;
        minMass=*absMin;
    }
    
    //printf("maxVol  %d   minVol %d\n", maxVol, minVol);
    objPtr = firstObj;
    while (objPtr != NULL)
    {
        mass=objPtr->mass;
        /* normilize object mass over range */
        if(minMass!=maxMass)
            scalar=(float)(mass-minMass)/(maxMass-minMass);
        else
            scalar=1.0;
        if(scalar > 1.0)        // object above colormap fixed top value
        {
            objPtr->color.r = 255;
            objPtr->color.g = 0;
            objPtr->color.b = 0;
        }
        else if(scalar < 0.0)   // object under colormap fixed bottom value
        {
            objPtr->color.r = 0;
            objPtr->color.g = 0;
            objPtr->color.b = 255;
        }
        else
            if(scalar <= 0.5)
            {
                objPtr->color.r = 0;
                objPtr->color.g = (int)((2.0 * scalar)*255);
                objPtr->color.b = (int)((1.0 - 2.0*scalar)*255);
            }
            else if(scalar <= 1.0)
            {
                objPtr->color.r = (int)((2.0 * (scalar - 0.5))*255);
                objPtr->color.g = (int)((1.0 - 2.0*(scalar - 0.5))*255);
                objPtr->color.b = 0;
            }
        
        objPtr = objPtr->next;
    }
    return(0);
}

//
//
//template <class T>
//int FindPackets( InpObject<T> *inPtr)
//{
//    float min_delta_x, delta_x, delta_y, delta_z, current_angle, current_angle1,thrs_delta_x,thrs_delta_y,thrs_delta_z;
//    long packetID(0),found_something(0);
//    int currentObjNum1(0),currentObjNum2(0),currentObjNum3(0),BackwardCloserObj(0),ForwardCloserObj(0);
//    
//    stObject *objPtr1;
//    stObject *objPtr2;
//    stObject *minobjPtr;
//    stObject *initial_obj;
//    stObject *current_obj;
//    
//    //stPacket *Packs;
//    
//    thrs_delta_x = inPtr->thresh_deltax;
//    thrs_delta_y = inPtr->thresh_deltay;
//    thrs_delta_z = inPtr->thresh_deltaz;
//    objPtr1 = inPtr->pobject;
//    
//    while( objPtr1 != NULL )
//    {
//        // cout<<" inside the initial setting loop !!!"<<endl;
//        objPtr1->packetFlag = 0;
//        objPtr1 = objPtr1->next;
//    }
//    
//    objPtr1 = inPtr->pobject;
//    //Packs = inPtr->packets; // should I use this... I should think about this
//    
//    //object by object go and search.
//    while( objPtr1 != NULL )
//    {
//        currentObjNum1++;
//        //cout<<" inside the first loop !!!"<<endl;
//        if (objPtr1->packetFlag == 0)
//        {
//            //cout<<" inside the first loop if statement!!!"<<endl;
//            packetID = packetID + 1;
//            objPtr1->packetFlag = 1;
//            objPtr1->Packet_ID = packetID;
//            current_obj = objPtr1;
//            initial_obj = objPtr1;
//            
//        }
//        else // should I delete this else statement? mmm let me think on that...
//        {
//            objPtr1 = objPtr1->next;
//            continue;
//        }
//        
//        //min_delta_x = 100000;
//        found_something = 1;
//        while( found_something ) // look for positive side objects
//        {
//            //cout<<" inside the first loop found SOMETHINGGGGGG!!!"<<endl;
//            found_something = 0;
//            objPtr2 = NULL;
//            objPtr2 = inPtr->pobject;
//            minobjPtr = NULL;
//            currentObjNum2 = 0;
//            min_delta_x = 100000;
//            while(objPtr2 != NULL )
//            {
//                currentObjNum2++;
//                
//                if (objPtr2->packetFlag == 0)
//                {
//                    delta_x = objPtr2->centroidX - current_obj->centroidX;
//                    delta_y = objPtr2->centroidY - current_obj->centroidY;
//                    delta_z = objPtr2->centroidZ - current_obj->centroidZ;
//                    
//                    if ((delta_x<0) || (delta_z<0))
//                    {
//                        objPtr2 = objPtr2->next;
//                        continue;
//                        
//                    }
//                    if (delta_y < 0)
//                    {
//                        delta_y = (-1)*delta_y;  // take the absolute value of delta_y, since it may be on both sides.
//                    }
//                    current_angle = atan (delta_z / delta_x) * 180 / PI; // angle in terms of degree
//                    current_angle1 = atan (delta_y / delta_x) * 180 / PI; // angle in terms of degree
//                    //	  cout<<"thrs_delta_x thrs_delta_y and z are: "<<thrs_delta_x<<" "<<thrs_delta_y<<" "<<thrs_delta_z<<endl;
//                    //	  cout<<"deltaX/Y/Z and current angle/1 are: "<<delta_x<<" "<<delta_y<<" "<<delta_z<<" "<<current_angle<< " "<<current_angle1<< endl;
//                    if ( (delta_y < thrs_delta_y) && (delta_x < thrs_delta_x) &&  (delta_z < thrs_delta_z ) && ( min_delta_x > delta_x ) && (current_angle < 80) && (current_angle > -1) && (current_angle1 < 35) && (current_angle1 > -10) )
//                    {
//                        //            cout<<" inside the LOOP: objPTR2 and if conditions!!!"<<endl;
//                        // cout<<"FORWARD search found! "<<"deltaX/Y/Z and current angle/1 are: "<<delta_x<<" "<<delta_y<<" "<<delta_z<<" "<<current_angle<< " "<<current_angle1<< endl;
//                        min_delta_x = delta_x;
//                        found_something = 1;
//                        minobjPtr =  objPtr2;
//                        ForwardCloserObj = currentObjNum2;
//                    }
//                }
//                objPtr2 = objPtr2->next;
//            } //end of while (objPtr2)
//            
//            if (minobjPtr)
//            {
//                minobjPtr->Packet_ID = current_obj->Packet_ID;
//                minobjPtr->packetFlag = 1;
//                current_obj = minobjPtr;
//                //	cout<<"FORWARD search found! Objectpair is: "<<currentObjNum1<<" & "<<ForwardCloserObj <<endl;//" deltaX/Y/Z and current angle/1 are: "<<delta_x<<" "<<delta_y<<" "<<delta_z<<" "<<current_angle<< " "<<current_angle1<< endl;
//            }
//        } // end of while found_something
//        //cout<<" inside the loop after the first positive side loop !!!"<<endl;
//        
//        //-------------------------------- left side is done, now look for the right side on x axis...
//        //cout<<"current_obj & initial_obj are: "<< current_obj<<" and "<<initial_obj<<endl;
//        current_obj = initial_obj;
//        found_something = 1;
//        while( found_something ) // look for negative side of the current objectobjects
//        {
//            found_something = 0;
//            objPtr2 = NULL;
//            objPtr2 = inPtr->pobject;
//            minobjPtr = NULL;
//            currentObjNum3 = 0;
//            min_delta_x = 100000;
//            while(objPtr2 != NULL )
//            {
//                currentObjNum3++;
//                if (objPtr2->packetFlag == 0)
//                {
//                    delta_x = current_obj->centroidX - objPtr2->centroidX;
//                    delta_y = current_obj->centroidY - objPtr2->centroidY;
//                    delta_z = current_obj->centroidZ - objPtr2->centroidZ;
//                    if ((delta_x<0) || (delta_z<0))
//                    {
//                        objPtr2 = objPtr2->next;
//                        continue;
//                        
//                    }
//                    if (delta_y < 0)
//                    {
//                        delta_y = (-1)*delta_y;  // take the absolute value of delta_y, since it may be on both sides.
//                    }
//                    current_angle = atan (delta_z / delta_x) * 180 / PI; // angle in terms of degree
//                    current_angle1 = atan (delta_y / delta_x) * 180 / PI; // angle in terms of degree
//                    if ( (delta_y < thrs_delta_y) && (delta_x < thrs_delta_x) &&  (delta_z < thrs_delta_z ) && ( min_delta_x > delta_x ) && (current_angle < 80) && (current_angle > -1) && (current_angle1 < 35) && (current_angle1 > -10) )
//                    {
//                        //    cout<<"BACKWARD search found something too!!! "<<endl;
//                        min_delta_x = delta_x;
//                        found_something = 1;
//                        minobjPtr =  objPtr2;
//                        BackwardCloserObj = currentObjNum3;
//                    }
//                }
//                objPtr2 = objPtr2->next;
//            } //end of while (objPtr2)
//            
//            if (minobjPtr)
//            {
//                minobjPtr->Packet_ID = current_obj->Packet_ID;
//                minobjPtr->packetFlag = 1;
//                current_obj = minobjPtr;
//                //	cout<<"Backward search found! Objectpair is: "<<currentObjNum1<<" & "<<BackwardCloserObj <<endl;//" deltaX/Y/Z and current angle/1 are: "<<delta_x<<" "<<delta_y<<" "<<delta_z<<" "<<current_angle<< " "<<current_angle1<< endl;
//            }
//            
//            
//        } // end of while found_something
//        
//        objPtr1 = objPtr1->next;
//    } // end of while
//    inPtr->MaxPacketnumber = packetID;
//    cout<<"the MAX Packet NUMBER in the current frame is: "<<packetID<<endl;
//    return METHOD_SUCCESS;
//}
//
//

template <class T>
int CreatePacketStats( InpObject<T> *inPtr, char *fileBaseName,int timeFrame)
{
    
    //-- file writing blog
    FILE *fpout;
    char extention[9]=".packet";
    char outPacketF[256];
    float time;
    time=(float)(timeFrame);
    strcpy(outPacketF,fileBaseName);
    strcat(outPacketF,extention);
    if((fpout = fopen(outPacketF, "w")) == NULL)
        cout << "cannot open outPacketF file to write\n";
    //-- file writing blog ends here
    
    
    long CurrentobjPacket_ID(0),Max_packetNum(0),pcounter1(0),pcounter(0),totalobjnumber(0);
    stObject *objPtr1;
    stPacket *CurrentPacket, *NextPacket;
    ObjIndex *CurObjlist, *Nextlist;
    
    objPtr1 = inPtr->pobject; // first segmented object
    
    Max_packetNum = inPtr->MaxPacketnumber;
    pcounter = 1;
    
    while( pcounter < Max_packetNum+1 )
    {
        if (inPtr->packets ==NULL)
            NextPacket = NULL;
        else
            NextPacket = inPtr->packets;//  ->next = CurrentPacket;
        CurrentPacket = new stPacket;
        CurrentPacket->next = NextPacket;
        inPtr->packets = CurrentPacket;
        //initial values for the current packet
        CurrentPacket->P_ID = pcounter;
        CurrentPacket->numObjs = 0;
        CurrentPacket->PacketCx = 0;
        CurrentPacket->PacketCy = 0;
        CurrentPacket->PacketCz = 0;
        CurrentPacket->PacketVolume = 0;
        CurrentPacket->PacketMass = 0;
        CurrentPacket->minX = 99999999;
        CurrentPacket->minY = 99999999;
        CurrentPacket->minZ = 99999999;
        CurrentPacket->maxX = -9999999;
        CurrentPacket->maxY = -9999999;
        CurrentPacket->maxZ = -9999999;
        
        objPtr1 = inPtr->pobject; // first segmented object
        
        pcounter1=0;
        CurrentPacket->ObjList = NULL;
        while( objPtr1 != NULL )
        {
            if (pcounter == objPtr1->Packet_ID)
            {
                pcounter1++;
                if (pcounter1==1) // if its the first object..
                {
                    CurrentPacket->minX = 9999999;
                    CurrentPacket->minY = 9999999;
                    CurrentPacket->minZ = 9999999;
                    CurrentPacket->maxX = -9999999;
                    CurrentPacket->maxY = -9999999;
                    CurrentPacket->maxZ = -9999999;
                }
                //	if (pcounter1==1) // if its the first object..
                //CurObjlist = NULL;
                //else
                CurObjlist = CurrentPacket->ObjList;
                Nextlist = new ObjIndex;
                //assert(Nextlist);
                Nextlist->objNumber = objPtr1->objNum;
                Nextlist->next = CurObjlist;
                CurrentPacket->ObjList = Nextlist;
                CurrentPacket->numObjs++;
                CurrentPacket->PacketVolume = (CurrentPacket->PacketVolume) + (objPtr1->volume);
                
                CurrentPacket->PacketMass = CurrentPacket->PacketMass + objPtr1->mass;
                CurrentPacket->PacketCx = CurrentPacket->PacketCx + (objPtr1->mass ) * (objPtr1->centroidX);
                CurrentPacket->PacketCy = CurrentPacket->PacketCy + (objPtr1->mass ) * (objPtr1->centroidY);
                CurrentPacket->PacketCz = CurrentPacket->PacketCz + (objPtr1->mass ) * (objPtr1->centroidZ);
                
                // find the extends..
                
                if (objPtr1->LowerLeft_Corner_x < CurrentPacket->minX )
                    CurrentPacket->minX = objPtr1->LowerLeft_Corner_x;
                if (objPtr1->LowerLeft_Corner_y < CurrentPacket->minY )
                    CurrentPacket->minY = objPtr1->LowerLeft_Corner_y;
                if (objPtr1->LowerLeft_Corner_z < CurrentPacket->minZ )
                    CurrentPacket->minZ = objPtr1->LowerLeft_Corner_z;
                if (objPtr1->UpperRight_Corner_x > CurrentPacket->maxX )
                    CurrentPacket->maxX = objPtr1->UpperRight_Corner_x;
                if (objPtr1->UpperRight_Corner_y > CurrentPacket->maxY )
                    CurrentPacket->maxY = objPtr1->UpperRight_Corner_y;
                if (objPtr1->UpperRight_Corner_z > CurrentPacket->maxZ )
                    CurrentPacket->maxZ = objPtr1->UpperRight_Corner_z;
                
            }
            objPtr1 = objPtr1->next;
        }
        CurrentPacket->PacketCx = CurrentPacket->PacketCx /CurrentPacket->PacketMass;
        CurrentPacket->PacketCy = CurrentPacket->PacketCy /CurrentPacket->PacketMass;
        CurrentPacket->PacketCz = CurrentPacket->PacketCz /CurrentPacket->PacketMass;
        
        //-- file writing blog
        fprintf(fpout,"%5.0ld %5.0ld %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %5.0ld ",
                pcounter,CurrentPacket->PacketVolume,CurrentPacket->PacketMass,
                CurrentPacket->PacketCx, CurrentPacket->PacketCy, CurrentPacket->PacketCz,
                CurrentPacket->minX, CurrentPacket->minY, CurrentPacket->minZ, CurrentPacket->maxX, CurrentPacket->maxY, CurrentPacket->maxZ,
                CurrentPacket->numObjs);
        CurObjlist = inPtr->packets->ObjList;
        while(CurObjlist!=NULL)
        {
            if (CurObjlist->objNumber == 0)
                fprintf(fpout,"     0 ");
            else
                fprintf(fpout," %5.0ld ",CurObjlist->objNumber);
            CurObjlist = CurObjlist->next;
        }
        fprintf(fpout," \n");
        // file writing blog ends
        /*
         //---for debugging purpose only
         //  cout<<"Within the WHILE loop!! "<<endl;
         CurObjlist = inPtr->packets->ObjList;
         cout<<"the Packet Cx/Cy/Cz ARE "<<CurrentPacket->PacketCx<<" "<<CurrentPacket->PacketCy<<" "<<CurrentPacket->PacketCz<<endl;
         cout<<"Theobjects within the Packet ID: "<<pcounter<<" ARE ";
         while(CurObjlist!=NULL)
         {
         //      if (CurObjlist->objNumber==NULL)
         //	cout<<"CurObjlist->objNumber==NULL  !!!!!!!!! "<<endl;
         cout<<CurObjlist->objNumber<<" ";
         CurObjlist = CurObjlist->next;
         }
         cout<<" "<<endl;
         cout<<"Total Object number in teh current Packet is: "<< CurrentPacket->numObjs<<endl;
         //--debugging ends here
         */
        pcounter++;// This packet is done, go to the next packet
    }
    fclose(fpout);
    return METHOD_SUCCESS;
}

/* // this method is not complete and needs to be edited.
 //But nevermind that, since I combined it with the CreatePacketStats method.
 int OutputPacket( stPacket *Packet, char *fileBaseName,int timeFrame)
 {
 FILE *fpout;
 stPacket *objPtr;
 ObjIndex *CurObjlist;
 char extention[9]=".packet";
 char outPacket[256];
 float time;
 time=(float)(timeFrame);
 strcpy(outPacket,fileBaseName);
 strcat(outPacket,extention);
 objPtr = Packet;
 CurObjlist= Packet->ObjList;
 if((fpout = fopen(outPacket, "w"))==NULL)
 cout << "cannot open outPacket file to write\n";
 while( objPtr != NULL )
 {
 fprintf(fpout,"%f  ",time);
 fprintf(fpout,"%5.0d %10.6f %10.6f %10.6f %10.6f %5.0d ",
 objPtr->P_ID,objPtr->PacketMass,
 objPtr->PacketCx,objPtr->PacketCy,objPtr->PacketCz,
 objPtr->numObjs);
 CurObjlist= objPtr->ObjList;
 
 while(CurObjlist!=NULL)
 {
 fprintf(fpout,"%5.0d ",CurObjlist->objNumber);
 CurObjlist = CurObjlist->next;
 }
 fprintf(fpout," \n",CurObjlist->objNumber);
 //           fprintf(fpout,"%14.6f   %5.0d   %10.6f  %10.6f  %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %5.0d\n",
 //                   objPtr->mass,objPtr->volume,objPtr->centroidX,
 //                   objPtr->centroidY, objPtr->centroidZ,objPtr->LowerLeft_Corner_x ,objPtr->LowerLeft_Corner_y ,objPtr->LowerLeft_Corner_z ,objPtr->UpperRight_Corner_x ,objPtr->UpperRight_Corner_y , objPtr->UpperRight_Corner_z, objPtr->Packet_ID);
 objPtr = objPtr->next;
 }
 fclose(fpout);
 return METHOD_SUCCESS;
 }
 
 */




template<class T>
int PopulateInputData(InpObject<T> *inPtr, unsigned long *node_conn_array, float *coord_array, double *node_data)
{
    cout << "inside PopulateInputData inPtr->nncomp: "<<inPtr->nncomp<<endl;
    
    //#cout << "Now in PopulateInputData \n";
    register unsigned long i=0, j=0;
    stNodePos    *nodePtr,   *curNode;
    stNodeData   *nodeDataPtr,  *curDataPtr;
    double max_node_data0=-99999.0, min_node_data0=99999.0;
    double max_node_data1=-99999.0, min_node_data1=99999.0;
    double max_node_data2=-99999.0, min_node_data2=99999.0;
    //assert(Get_Conn(inPtr, node_conn_array));
    Get_Conn(inPtr, node_conn_array);
    //get coords
    cout << "inside PopulateInputData 1.5 "<<endl;
    
    nodePtr= CreateNodeArray(inPtr->numNodes);
    assert(nodePtr);
    //#cout << "post creating node list\n";
    cout << "inside PopulateInputData 1.6 "<<endl;
    for (j=0; j < inPtr->numNodes; j++)
    {
        i=3*j; // 3 is the dimension number (x, y, z) --------> i.e. the space is 3 dimensional
        curNode = GetCurNode( nodePtr, j);
        curNode->x = coord_array[i];
        curNode->y = coord_array[i+1];
        curNode->z = coord_array[i+2];
        curNode->list = NULL;
        curNode->flag = -1;
    }
    //get node_data
    cout << "inside PopulateInputData 1.7 "<<endl;
    nodeDataPtr = CreateDataArray( inPtr->numNodes);
    assert(nodeDataPtr);
    //#cout << "post creating data list\n";
    //cout << "inside PopulateInputData 2 "<<endl;
    cout << "inside PopulateInputData 1.8 "<<endl;
    
    
    // from SEDAT to SEDAT: in the future these nncomp==1 and ==3 should be combined for a generic case!!!
    if (inPtr->nncomp == 1)
    {
        for (j = 0; j < inPtr->numNodes; j++)
        {
            curDataPtr = GetCurNodeData( nodeDataPtr, j);
            curDataPtr->val0 = node_data[j];
            if( max_node_data0 < node_data[j] )
                max_node_data0 = node_data[j];
            if( min_node_data0 > node_data[j] )
                min_node_data0 = node_data[j];
            inPtr->pnode   = nodePtr;
            inPtr->pnodeData  = nodeDataPtr;
            inPtr->minData = min_node_data0;
            inPtr->maxData = max_node_data0;
            
            
        }
    }
    if (inPtr->nncomp == 3)
    {
        for (j = 0; j < inPtr->numNodes; j++)
        {
            i = j*inPtr->nncomp;
            curDataPtr = GetCurNodeData( nodeDataPtr, j);
            curDataPtr->val0 = node_data[i];
            curDataPtr->val1 = node_data[i+1];
            curDataPtr->val2 = node_data[i+2];
            if( max_node_data0 < node_data[i] )
                max_node_data0 = node_data[i];
            if( min_node_data0 > node_data[i] )
                min_node_data0 = node_data[i];
            if( max_node_data1 < node_data[i+1] )
                max_node_data1 = node_data[i+1];
            if( min_node_data1 > node_data[i+1])
                min_node_data1 = node_data[i+1];
            if(max_node_data2 < node_data[i+2])
                max_node_data2 = node_data[i+2];
            if(min_node_data2 > node_data[i+2])
                min_node_data2 = node_data[i+2];
            inPtr->pnode   = nodePtr;
            inPtr->pnodeData  = nodeDataPtr;
            inPtr->minData = min_node_data0;
            inPtr->maxData = max_node_data0;
            inPtr->minData1 = min_node_data1;
            inPtr->maxData1 = max_node_data1;
            inPtr->minData2 = min_node_data2;
            inPtr->maxData2 = max_node_data2;
        }
    }
    cout << "inside PopulateInputData 1.9 "<<endl;
    /*
     for (j=0;j<inPtr->numNodes;j++)
     {
     curDataPtr = GetCurNodeData( nodeDataPtr, j);
     curDataPtr->val0=node_data[j];
     if(max_node_data<node_data[j])
     max_node_data=node_data[j];
     if(min_node_data>node_data[j])
     min_node_data=node_data[j];
     }
     */
    
    
    return METHOD_SUCCESS;
}





template <class T>
void FindClusters( InpObject<T> *inPtr)  // this is to use a generic way for grouping...
{
//    int n = 3;
//    int k = 2;
//    int numberofTotalObjects = 0; // total number of objects for the current frame
//    stObject * current_obj;
//    stObject * object_forCentroid;
//    stObject * object_toWriteItsPacketId;
//    
//    current_obj = inPtr->pobject;
//    
//    while( current_obj != NULL ) // initialize all the packetFlag's to 0.
//    {
//        // cout<<" inside the initial setting loop !!!"<<endl;
//        current_obj->packetFlag = 0;
//        current_obj = current_obj->next;
//        numberofTotalObjects++;
//    }
//    Point3D test_data[numberofTotalObjects];
//    int counter1 = 0;
//    int centroidX(0), centroidY(0),centroidZ(0);
//    //change centroid to max location...
//    
//    object_forCentroid = inPtr->pobject;
//    
//    while( object_forCentroid != NULL )
//    {
//        Point3D currentCentroid;
//        currentCentroid.x  = object_forCentroid->centroidX;
//        currentCentroid.y  = object_forCentroid->centroidY;
//        currentCentroid.z  = object_forCentroid->centroidZ;
//        test_data[counter1] = currentCentroid;
//        //cout<< "currentCentroid.x "<<currentCentroid.x << "currentCentroid.y "<<currentCentroid.y << "currentCentroid.z "<<currentCentroid.z << "counter1 "<< counter1<<endl;
//        object_forCentroid = object_forCentroid->next;
//        counter1++;
//    }
//    cout<<"Before the Classf function!! The  numberofTotalObjects: "<<  numberofTotalObjects<< endl;
//    classifier classf(test_data, numberofTotalObjects, numberofTotalObjects/2);
//    
//    cout<<"between the two clustering function !"<<endl;
//    
//    bool dumyval = classf.run_classification();
//    
//    cout<<"after the classification"<<endl;
//    
//    //  for(int ii=0; ii<numberofTotalObjects  ; ii++)
//    //    cout<<"the "<<ii<<"th objects PacketID: " <<classf.m_Index[ii]<<endl;
//    //Point3D a = classf.Points[0]; examples
//    //int group_index = classf;.m_Index[0]; examples
//    int MaxPacketnumber1 =  0;
//    object_toWriteItsPacketId = inPtr->pobject;
//    counter1 = 0;
//    while( object_toWriteItsPacketId != NULL ) // initialize all the packetFlag's to 0.
//    {
//        object_toWriteItsPacketId->packetFlag = 1;
//        object_toWriteItsPacketId->Packet_ID = classf.m_Index[counter1];
//        if (classf.m_Index[counter1]  > MaxPacketnumber1)
//            MaxPacketnumber1 = classf.m_Index[counter1];
//        object_toWriteItsPacketId = object_toWriteItsPacketId->next;
//        counter1++;
//    }
//    inPtr->MaxPacketnumber = MaxPacketnumber1;
//    cout<<"end of findcluster function !!!!!"<<endl;
}









/*

template <class T>
int Segmentation(vtkDataSet *in_ds,vtkDataSet **outDS, InpObject<T> *inPtr,int nspace, unsigned long nnodes, unsigned long ncells,int cellpoints, float *coord_array,unsigned long *node_conn_array, double *node_data,double thresh,std::string listfile, int timeFrame, char* mypolyfile, int smallest_vol,int precision,bool dovolren string volfile ,int color_Method,int nncomp)
{
    cout<<"inside Segmentation  "<<endl;
    cout<< "SEDAT1: mypolyfile is: "<< mypolyfile << endl;
    //#cout << "Now in Segmentation \n";
    inPtr->numCell=ncells;
    inPtr->nncomp=nncomp;
    inPtr->numNodes=nnodes;
    //assert(PopulateInputData(inPtr, node_conn_array, coord_array, node_data));
    PopulateInputData(inPtr, node_conn_array, coord_array, node_data);     //# cout << " post PopulateInputData\n";
    // put the vtk data into an inpObject structure
    delete [] node_conn_array;
    node_conn_array = NULL;
    delete [] node_data;
    node_data = NULL;
    delete [] coord_array;
    coord_array = NULL;
    inPtr->colorMethod = color_Method; //as per user preference
    inPtr->thrType = THRESH_PERCENT;
    //assert( ProcessThresholds(inPtr, thresh, inPtr->thrType ) );
    cout<<"before precessthresholds"<<endl;
    ProcessThresholds(inPtr, thresh, inPtr->thrType);
    cout << "after ProcessThresholds \n";
    ProcessCells(inPtr);
    //assert( ProcessCells(inPtr));
    cout << "after ProcessCells\n";
    
    // build a list of cells connected to a given point.
    
    char *mergeArray = new char[ncells];
    assert( mergeArray );
    cout<<"Before SegmentObjects"<<endl;
    //assert( SegmentObjects( inPtr, smallest_vol, mergeArray));
    SegmentObjects( inPtr, smallest_vol, mergeArray);
    cout << "after SegmentObjects\n";
    
    //assert( SetObjMax(inPtr,  cellpoints)); // find max and integrated content
    SetObjMax(inPtr,  cellpoints);     //#cout << "Post SetObjMax\n";
    FindCentAndVol( inPtr, cellpoints);
    //assert( FindCentAndVol( inPtr, cellpoints ));
    
    
    
    float absMin, absMax;
    
    //#cout<<"\nColor Method has value : "<<inPtr->colorMethod;
    
    switch ( inPtr->colorMethod )
    {
            //case 1:// rc = setColor(inPtr, inPtr->listObj, inPtr->thrType);
            //  break;
        case 0: SetColVol( inPtr, &absMin , &absMax, ADAPTIVE); // Because ADAPTIVE is used,
            // the two middle parameters doesn't matter
            break;
        case 1: SetColMass( inPtr, &absMin, &absMax, ADAPTIVE);
            break;
        case 2: /*setColRandom(inPtr,&absMin,&absMax,ADAPTIVE);
            break;
    }
    string label;
    label=getlabel(listfile);
    //#cout << "timeFrame: " <<  timeFrame<< endl;
    if(!label.size())
        cout << "Segmentation : cannot get label\n";
    char outFile[256];
    /* if(timeFrame==0)
     sprintf(mypolyfile, "%s",label.c_str());
     else
     {
    string buf=label+precision_time( timeFrame, precision);
    sprintf(mypolyfile, buf.c_str());
    //1}
    strcpy(outFile,mypolyfile);
    sprintf(mypolyfile, "%s.poly", mypolyfile);
    cout<< "SEDAT2: mypolyfile is: "<< mypolyfile << endl;
    delete [] mergeArray;
    //#cout << "No. of objects found is : " << inPtr->numObj;
    cout<<"Before the Extract Surf"<<endl;
    ExtractSurf(in_ds, inPtr,outDS, mypolyfile,nspace,cellpoints);
    //assert(ExtractSurf(inPtr,outDS, mypolyfile,nspace,cellpoints));
    cout<<"After the Extract Surf"<<endl;
    //#cout << "Post ExtractSurf\n";
    
    //assert(FindMoment(inPtr,cellpoints));
    FindMoment(inPtr,cellpoints);
    
    //SEDAT added
    
    cout<<" BEFORE THE PACKET FINDING !!!"<<endl;
    //assert(FindPackets( inPtr));
    //assert(FindClusters( inPtr));
    FindClusters( inPtr);
    
    CreatePacketStats( inPtr,outFile, timeFrame);
    //assert(CreatePacketStats( inPtr,outFile, timeFrame));
    cout<<" AFTER THE PACKET FINDING !!!"<<endl;
    
    OutputAttribute(inPtr->pobject, outFile);
    OutputTrak(inPtr->pobject, outFile, timeFrame);
    OutputOcd(inPtr, outFile, timeFrame,dovolren,cellpoints/*volfile,in_field);
    //assert(OutputAttribute(inPtr->pobject, outFile));
    //assert(OutputTrak(inPtr->pobject, outFile, timeFrame));
    //assert(OutputOcd(inPtr, outFile, timeFrame,dovolren,cellpoints/*volfile,in_field));
    //cout<<"Before OUTPUTPACKET method"<<endl;
    //     assert( OutputPacket( inPtr->packets, outFile, timeFrame));
    //cout<<"After OUTPUTPACKET method"<<endl;
    
    return METHOD_SUCCESS;
}
*/



#endif



