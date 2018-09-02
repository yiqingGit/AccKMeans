#include<stdio.h>
#include<iostream>
#include<fstream>
#include<time.h>
#include<math.h>
#include<string>

#ifndef MY_KMEANS
#define MY_KMEANS

#define M	31696					//表示数据的个数上界
#define	N	512				//表示维度的上界
#define	K	320					//表示聚类的个数上界
#define INF 9999999


using namespace std;

string infile;
string outfile;

double **dataSet;
int clusterNum;							// cluster number
int dim;						// dimension of dataSet
int num;		//the number of dataSet

double **distC2C;			// dist between centroids,the minimum distance store in diagonal 
double *u;					//  upper bound
double **centroid;			// centroid 
int *newCentroidNum;	//the count point of centroid include
double **newCentroid;	// new centroid after iter
int *indexD2C;		// index of centroid which data belong to 

/*
string infile;
string outfile;

double dataSet[M][N];
int clusterNum;							// cluster number
int dim;						// dimension of dataSet
int num;		//the number of dataSet

double distC2C[K][K];			// dist between centroids,the minimum distance store in diagonal 
double u[M];					//  upper bound
double centroid[K][N];			// centroid 
int newCentroidNum[K];	//the count point of centroid include
double newCentroid[K][N];	// new centroid after iter
int indexD2C[M];		// index of centroid which data belong to 
*/

int iterCount;	// count of iterations
bool clusterChange;

bool *r;


void alloc_All();
void init_All();
void free_All();
bool readFile(string filename);
void writeFile(string filename);
inline double power2(double x);
double distK_X(int k,int x);
double distK_K(int k1,int k2);
double distK_NewK(int k1,int k2);
void minC2C();
void add_NewCentroid(int k, int x);
void update_NewCentroid();
void update_Centroid();
void iterCal();
void accKmeans();
bool programArgs(int argc, char *argv[]);
void initCentroid();



bool centroidDataSame(int n);
void initRandomCentroid();
void importCentroid(string filename);
void initSequenceCentroid();

#endif