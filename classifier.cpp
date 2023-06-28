#include "classifier.h"
#include <stdio.h>
// constructor

classifier::classifier(Point3D* pPoints, int n, int k, float deltaxval, float deltayval, float deltazval)
{
 	// get number of classes to classify
	std::cout<<"classifier created with k = " << k << " and n = " << n << endl; 
	K = k;
	N = n;
    Parameter1 = deltaxval;
    Parameter2 = deltayval;
    Parameter3 = deltazval;
    Points.clear();
    for(int i = 0; i < n; i++){
        Points.push_back(pPoints[i]);
    }
    
    m_Index = new int[N];
}

classifier::~classifier(void)
{
  if(m_Index != NULL){
    delete []m_Index;
    m_Index = NULL;
  }
}


double classifier::FindDistance(Point3D a, Point3D b)
{
 //  Finds Distance and returns the inverse of the distance as a similarity metric 
 //  if the Distance is smaller than the given threshold

    double threshold = 0.2;
    double delta_x= (a.x - b.x);
    double delta_y= (a.y - b.y);
    double delta_z= (a.z - b.z);
    double A = sqrt(delta_x*delta_x + delta_y*delta_y + delta_z*delta_z) ;
    if (A < threshold && A != 0)
      A = 1/A;
    else
      A = 0;
  
    return A;
}


double classifier::FindSimilarity(Point3D a, Point3D b)
{
    // ThisDataType =1 means data is Hairpin and grouping is based on the packet information; 
    // ThisDataType =2 means data is something and the similarity metric is based on the distance only; 
    int ThisDataType = 2;
    
    double Thresh_x = (double) Parameter1;
    double Thresh_y = (double) Parameter2;
    double Thresh_z = (double) Parameter3;
    double A;
    
 if (ThisDataType == 2)
 {
     double delta_x= (a.x - b.x);
     double delta_y= (a.y - b.y);
     double delta_z= (a.z - b.z);
     
     double Dx = Thresh_x - fabs(delta_x);
     double Dy = Thresh_y - fabs(delta_y);
     double Dz = Thresh_z - fabs(delta_z);
     
     A = (sign(Dx)+1)*(sign(Dy)+1)*(sign(Dz)+1)/8;
 
 }
    
    
//       double Thresh_x = 0.4;
//     double Thresh_y = 0.1;
//     double Thresh_z = 0.25;
//-----------------------------------
//       double Thresh_x = 0.25;
//     double Thresh_y = 0.1;
//     double Thresh_z = 0.2;
//     double MaxAngle = 3.1415926/4; 

    
    
    
    if (ThisDataType == 1)
    {

        double MaxAngle = 3.1415926/4;
        //45 degree
        //180*atan(norm(a(2)-b(2))/norm(a(1)-b(1)))/pi;
        double delta_x= (a.x - b.x);
        double delta_y= (a.y - b.y);
        double delta_z= (a.z - b.z);
        double Dx = Thresh_x - fabs(delta_x);
        double Dy = Thresh_y - fabs(delta_y);
        double Dz = Thresh_z - fabs(delta_z);
        
        double currentAngle= atan(delta_z/delta_x);
        double Anglesign = sign(currentAngle);
        double higher = (sign(delta_x*delta_z)+1)/2;
        double WithinBox = (sign(Dx)+1)*(sign(Dy)+1)*(sign(Dz)+1)/8;
        double AngleSmallerThanMax;
        if (Anglesign > 0)
            AngleSmallerThanMax = (sign(MaxAngle - currentAngle)+1)/2;
        else
            AngleSmallerThanMax = 0;
        
        double alltheconditionsTrue = higher * WithinBox * AngleSmallerThanMax;
        if (delta_x !=0)
            A = alltheconditionsTrue/fabs(delta_x);
        else if (delta_z !=0)
            A = alltheconditionsTrue/fabs(delta_z);
        else
            A = 0;
    }
    
    return A;
}

void classifier::classification()
{
	// fill distances matrix C, then I and P
	A.clear();
	C.clear();
	II.clear();
	P.clear();
	int i, j;
	C.resize(N);
	P.resize(N);
	II.resize(N);
	A.resize(N);
	// complexity is: O(N*N*log(N))
	std::cout<<".. classification running .. " << endl << "N = " << N << endl;
 	for(i=0; i<N; ++i)
	{
		std::vector<similarities> V_temp;
		similarities D_temp;
		V_temp.clear();
		std::multiset<similarities, Cmp> Q_temp;
		Q_temp.clear();
		V_temp.resize(N);
		for(j=0; j<N; ++j)
		{
		    D_temp.simi = FindSimilarity(Points[i], Points[j]);
		    //D_temp.simi = FindDistance(Points[i], Points[j]);
			D_temp.index = j;
			V_temp[j]=D_temp;
			if(j!=i)
			{
				Q_temp.insert(D_temp);
			}
			//cout<<"Index "<<i<<" with Index "<<j<<"'s similarity is "<<D_temp.simi<<endl;
		}
		C[i]=V_temp;
		II[i]=1;
		P[i]=Q_temp;
        
		std::vector<int> A_i;
		A_i.push_back(i);
		A[i]=A_i;
	}
    
	// next: build clusters until K clusters left
	// complexity is: O(N*N*log(N))
	for(int n=0; n<N; ++n)
	{
	  
	  //cout<<"The iteration index is "<<n<<endl;
	  //cout<<"The N is: "<<N<< " The K value is: " << K <<endl;
		double max_simi = MAX_SIMI;
		int max_index = 0;
		double sum_simi = 0;
		for(int k=0; k<N-1; ++k)
		{
			if(II[k]==1)
			{
				if(P[k].begin()->simi > max_simi)
				{
					max_simi = P[k].begin()->simi;
					max_index = P[k].begin()->index;
				}
			  sum_simi += P[k].begin()->simi;
			}

		}
		
		if(sum_simi == 0.0)
		{
// 		  for(int k = 0; k < N -1; ++k){
// 		    if(II[k]==1)
// 			{
// 			  cout<<"The maximum simi for each point "<<P[k].begin()->simi<<endl;
// 			  
// 			}
//		  }
		  break;
		}
		//cout<<"maximum similarity "<<n<<endl;
		// we have minimum distance
		// k1, k2 - indexes of most nearest clusters
		int k1 = max_index;
		int k2 = P[k1].begin()->index;
		
		int N_k1 = A[k1].size();
		int N_k2 = A[k2].size();
		
		P[k1].clear();
		// add cluster k2 to A[k1]
		for(int l=0;l<A[k2].size();++l)
		{
			A[k1].push_back(A[k2][l]);
		}
		//cout<<"cluster reorganization "<<n<<endl;
		// clear the second cluster
		II[k2] = 0;
		// O(N*log(N))
		for(int m=0; m<N; ++m)
		{
			// O(log(N)): insert, erase operations
			if((II[m]!=0)&&(m!=k1))
			{
				P[m].erase(C[m][k1]);
				P[m].erase(C[m][k2]);
                
                C[m][k1].simi = C[m][k1].simi > C[m][k2].simi ? C[m][k1].simi : C[m][k2].simi;
                C[k1][m].simi = C[m][k1].simi;
				
				P[m].insert(C[m][k1]);					
				P[k1].insert(C[k1][m]);
			}
		}
		//cout<<"similarity update"<<n<<endl;
	}
}


// printing classes
void classifier::print_classes()
{
	int class_num = 0;
	for(int i=0; i<N; ++i)
	{
		if(II[i]==1)
		{
			++class_num;
			//cout<<std::endl<<"Class number: "<<class_num<<std::endl<<std::endl;
			std::vector<int>::iterator it_a(A[i].begin());
			while(it_a!=A[i].end())
			{
				int index = (*it_a);
				//cout<<Points[index].x<<" "<<Points[index].y<<" "<<Points[index].z<<std::endl;
				++it_a;
				m_Index[index] = class_num; 
			}
		}
	}
}


// run classification
bool classifier::run_classification(void)
{
	// check K
	if((K<1)||(K>N))
	{
		std::cout<<"Wrong K (number of clusters). Should be in [1, N] but K = " << K <<std::endl;	
		
		return false;
	}

	classification();
    
    print_classes();
    
	return true;
}
