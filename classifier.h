#include "iostream"
#include "vector"
#include "set"
#include "math.h"

using namespace std; 

struct Point3D 
{
    double x;
    double y;
    double z;
   // struct Point3D * next;  Sedat added this line here just as a thought.
};

struct similarities
{
    double simi;
    int index;
};

struct Cmp{
    bool operator()(const similarities s1, const similarities s2)
    const
    {
        if(s1.simi == s2.simi)
        {
            return s1.index < s2.index;
        }
        return s1.simi > s2.simi;
    }
};

const int MAX_SIMI=0;


class classifier
{
private:
    
    inline double sign(double x){
        if(x > 0){
            return 1.0;
        }
        else if(x == 0){
            return 0.0;
        }
        else{
            return -1.0;
        }
    }

    // similarity between 3D points
    double FindSimilarity(Point3D a, Point3D b);
    double FindDistance(Point3D a, Point3D b);
	// multiset for storing sorted similarities
	std::vector< std::multiset<similarities, Cmp> > P;
	// 2d vector for storing similarity matrix
	std::vector< std::vector<similarities> > C;
	// vector for storing flags for marking currently active clusters
	std::vector<int> II;
	// 2d vector for storing lists of data points in clusters
	std::vector< std::vector<int> > A;

	// classification with provided linkage_criteria and distance_function
	void classification();
	
public:
	// number of classes to classify
	int K;
	// number of titles
	int N;
    float Parameter1, Parameter2, Parameter3;

	std::vector<Point3D> Points;
	int * m_Index;

	// total classification (using all linkage criteries and distances)
	bool run_classification(void);
	void print_classes();
	// constructor
	classifier(Point3D* pPoints, int n, int k, float deltaxval, float deltayval, float deltazval);
	// destructor
	~classifier(void);
};

