#include "AccKmeans.h"


/////////////////////////////////////////////////////////////////////
bool centroidDataSame(int n){
	for(int i=0;i<n;i++){
		int sameSum=0;
		for(int j=0;j<dim;j++){
			if(centroid[n][j]!=centroid[i][j]){		
				break;
			}else{
				sameSum ++;
			}
		}
		if(sameSum==dim){
			return true;
		}
	}	
	return false;
}
// init centroid by select random data from dataSet
void initRandomCentroid(){
	srand((unsigned)time(NULL));		//initialize random number generator
	int n=0;
	while(n < clusterNum){
		int point = rand() % num;
		cout<<"point:"<<point<<endl;
		for(int j=0;j<dim;j++){
			centroid[n][j] = dataSet[point][j];
		}
		while(centroidDataSame(n)){
			point = rand() % num;
			for(int j=0;j<dim;j++){
				centroid[n][j] = dataSet[point][j];
			}
		}
		n++;
	}
	cout<<"total select centroid:"<< n <<endl;
}
//init centroid by import file 
void importCentroid(string filename){
	ifstream inStream;
	inStream.open(filename,ifstream::in);
	if(!inStream){
		cerr<<"error,unable to open the file"<<inStream<<endl;
		return ;
	}
	int point=0;
	while(!inStream.eof()){
		 int row=0;
		 inStream>>row;
		 for(int j=0;j<dim;j++){
			centroid[point][j]=dataSet[row][j];
		 }
		 point++;
	}
	cout<<"total import centroid(start by 0) : "<<(--point)<<endl;

}

// init centroid by sequence
void initSequenceCentroid(){
	for(int point=0;point< clusterNum;point++){								
		for(int d=0;d<dim;d++){
			centroid[point][d]=dataSet[point][d];
		}
	}
}
//////////////////////////////////////////////////////////////

void alloc_All(){
	dataSet = new double *[num];
	for(int i=0;i<num;i++)
		dataSet[i] = new double[dim];
	
	distC2C = new double*[clusterNum];
	for(int i=0;i<clusterNum;i++)
		distC2C[i]= new double[clusterNum];

	u = new double[num];

	centroid= new double*[clusterNum];
	for(int i=0;i<clusterNum;i++)
		centroid[i]= new double[dim];

	newCentroidNum = new int[clusterNum];

	newCentroid= new double *[clusterNum];
	for(int i=0;i<clusterNum;i++)
		newCentroid[i]= new double[dim];
	
	indexD2C = new int[num];
}

void init_All(){
	
	for(int i=0;i<num;i++)
		for(int j=0;j<dim;j++)
			dataSet[i][j] = 0;

	for(int i=0;i<clusterNum;i++)
		for(int j=0;j<clusterNum;j++)
			distC2C[i][j] = 0;

	for(int i=0;i<num;i++)
		u[i] = 0;

	for(int i=0;i<clusterNum;i++)
		for(int j=0;j<dim;j++)
			centroid[i][j] = 0;

	for(int i=0;i<clusterNum;i++)
		newCentroidNum[i] = 0;

	for(int i=0;i<clusterNum;i++)
		for(int j=0;j<dim;j++)
			newCentroid[i][j] = 0;

	for(int i=0;i<num;i++)
		indexD2C[i] = 0;
	
}
void free_All(){
	for(int i=0;i<num;i++)
		delete[] dataSet[i];
	delete []dataSet;

	for(int i=0;i<clusterNum;i++)
		delete[] distC2C[i];
	delete []distC2C;

	delete []u;

	for(int i=0;i<clusterNum;i++)
		delete[] centroid[i];
	delete []centroid;

	delete []newCentroidNum;

	for(int i=0;i<clusterNum;i++)
		delete[] newCentroid[i];
	delete []newCentroid;

	delete []indexD2C;
}

bool readFile(string filename){
	ifstream inStream;
	inStream.open(filename , ifstream::in);
	if(!inStream){
		cerr<<"error,unable to open the file"<<inStream<<endl;
		return false;
	}
	int x=0;
	while(!inStream.eof()){
		for(int i=0;i<dim;i++){
			inStream>>dataSet[x][i];
		}
		x ++;
	}
	inStream.close();
	return true;
}
void writeFile(string filename){
	ofstream outStream;
	outStream.open(filename,ofstream::out);
	int i=0;
	while(i<num){
		//outStream<<dataSet[i][0]<<"\t"<<dataSet[i][1]<<"\t"<<dataSet[i][2]<<"\t"<<dataSet[i][3]<<"\t";
		outStream<<(indexD2C[i]+1)<<endl;		// plus 1 means to convert index start with 1
		i++;
	}
	outStream.flush();
	outStream.close();
}

inline double power2(double x){
	return x * x;
}

double distK_X(int k,int x){
	double dist =0 ;
	for(int j=0;j < dim ;j++ ){
		dist += power2(dataSet[x][j] - centroid[k][j]);
	}
	return sqrt(dist);
}
double distK_K(int k1,int k2){
	double dist =0 ;
	for(int j=0;j < dim ;j++ ){
		dist += power2(centroid[k1][j] - centroid[k2][j]);
	}
	return sqrt(dist);
}

double distK_NewK(int k1,int k2){
	double dist =0 ;
	for(int j=0;j < dim ;j++ ){
		dist += power2(centroid[k1][j] - newCentroid[k2][j]);
	}
	return sqrt(dist);
}

void minC2C(){	
	int k1=0,k2=0;
	for(k1=0;k1< clusterNum;k1++){				
		for(k2=k1+1;k2< clusterNum;k2++){
			double dist=0;
			dist=distK_K(k1,k2);
			distC2C[k1][k2]=dist;
			distC2C[k2][k1]=dist;	
		}
		double min=INF;
		for(int t=0;t< clusterNum;t++){
			if(distC2C[k1][t] < min && k1!=t){
				min=distC2C[k1][t];
			}
		}
		distC2C[k1][k1]= min;
	}
}
void add_NewCentroid(int k, int x){
	newCentroidNum[k] ++;
	for(int i=0;i<dim;i++){
		newCentroid[k][i] +=dataSet[x][i];
	}
}
void update_NewCentroid(){
	for(int k=0;k< clusterNum;k++){
		for(int i=0;i<dim;i++){
			newCentroid[k][i]= newCentroid[k][i] / newCentroidNum[k];
		}
	}
}
void update_Centroid(){
	for(int i=0;i< clusterNum;i++){
		for(int j=0;j<dim;j++){					
			centroid[i][j]=newCentroid[i][j];
			newCentroid[i][j]=0;		
		}
		newCentroidNum[i]=0;
	}
}

void iterCal(){
	clusterChange = false;
	for(int x = 0;x < num ; x ++){
		//cout<<".";
		int cx = indexD2C[x];
		if(u[x] <= 0.5 * distC2C[cx][cx]){
			add_NewCentroid(cx , x);
			continue ;
		}
		double dcx = distK_X(cx,x);
		for(int c=0;c < clusterNum ; c++){
			r[x] = (c!=cx) && (u[x]>dcx) &&( u[x]>0.5* distC2C[cx][c]);
			if(r[x]){
				u[x] = dcx ;
				r[x] = false;
			}
			dcx = u[x];
			if(dcx > 0.5 * distC2C[cx][c]){
				double dc = distK_X(c , x);
				if(dc < dcx){
					cx = c;
					u[x] = dc;
				}
			}
		}
		if(cx != indexD2C[x]){
			indexD2C[x] = cx;
			clusterChange = true;
		}		
		add_NewCentroid(cx , x);
	}
	update_NewCentroid();
	for(int x = 0;x < num ;x ++){
		int k = indexD2C[x];
		u[x] = u[x] + distK_NewK(k,k);
		r[x] = true ;
	}	
}
void accKmeans(){
	clusterChange = true;
	iterCount = 0;
	r= new bool[num];
	while(clusterChange && iterCount < 1000){
		cout<<"iter : "<< iterCount <<endl;
		minC2C();
		iterCal();
		update_Centroid();			
		iterCount ++;
	}
	delete []r;
}
bool programArgs(int argc, char *argv[]){
	dim= 0;
	clusterNum = 0;
	num =0;
	if(argc>=5){		
		infile =  argv[1];
		outfile = argv[2];
		num = atoi(argv[3]);
		dim = atoi(argv[4]);
		clusterNum = atoi(argv[5]);
		cout<<"inputFilePath: "<<infile<<endl;
		cout<<"outFilePath: "<<outfile<<endl;
		cout<<"num :"<<num<<endl;
		cout<<"dimisions :"<<dim<<endl;
		cout<<"cluster number :"<<clusterNum<<endl;
		if(num<=0 || dim <=0 || clusterNum <=0){
			cout<<endl<<">> bad parameter !"<<endl;
			return false;
		}
		if(clusterNum >=num){
			cout<<endl<<">> cluster number is should not larger than the samples !"<<endl;
			return false;
		}		
		return true;
	}else{
		cout<<"[--help] \n"<<endl;
		cout<<"用法： AccKmeans 源文件（inputFilePath） 类文件（ouputFilePath）样本数(num) 维度（dim）类数(clusterNum)"<<endl;
		return false;
	}
}
void initCentroid(){
	initRandomCentroid();				//init random centroid
	//initSequenceCentroid();			// init sequence centroid
	for(int x=0;x<num;x++){		
		int cx=0;
		indexD2C[x]=cx;						//init every samples belongs to first centroid
		u[x]=distK_X(cx,x);			//
//		cout<<u[x]<<endl;
	}
}

int main(int argc, char *argv[])
{
	bool b = programArgs(argc,argv);
	if(!b){
		return 0;
	}
	alloc_All();
	init_All();
	clock_t t;	
	if(!readFile(infile)){
		return 0;
	}
	cout<< num <<" samples in total import from dataSet"<<endl;
	initCentroid();	
	t=clock();			// time clock
	accKmeans();
	t=clock()-t;		// time clock
	cout<<"total run "<<(double)t<<" milliseconds <"<<(double)t/CLOCKS_PER_SEC<<" seconds >"<<endl;
	writeFile(outfile); 
	free_All();	
	return 0;
}
