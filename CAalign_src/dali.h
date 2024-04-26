/*
typedef struct{
	int *real;//real number
	double **area,*exp;//order based contact area and exposed area
	int N,Max;
	double Tarea,TareaSq;//Total Area
}DMTX;
*/
typedef struct{
	unsigned int *id;
	double *bound,*suma,*sco;
	int N;
	unsigned int **stbl;
	int *Ns;
	double *single_sco;
}PLIST;

typedef struct{
	int *a1,*a2;//Alignment
	double sco,rate;
}ALIDATA;

typedef struct{
	//int **a1,**a2;//Alignment
	double *sco;
	double ave,std;
	int N,len1,len2;
	ALIDATA *ali;
}ALI;

typedef struct{
	int r[3];//triangle
	double area;//sum of contact area
	double upb;//upper bound
}TRIA;


int FindMaxResNum(char *);
int MallocDmtx(DMTX *,int);
int ReadDmtx(DMTX *,char *,bool);
int CalPairlist(PLIST *,DMTX *,DMTX *,double *,int,double,double,double,bool);
int GeneSeeds(ALI *,PLIST *,DMTX *,DMTX *,int,double,bool);
int GenePath(ALI *,int,int,int **,int *,int,DMTX *,DMTX *,double,bool);
int GeneTripletsAndAli(ALI *,int,int,int *,int,int,DMTX *,DMTX *);
int FilterSeeds(ALI *,DMTX *,DMTX *);
int Tri2Ali(int *,int *,int,int,int,int,int,int,int,int,int);
double Ali2Sco(int *,DMTX *,DMTX *,double *);
double Ali2ScoSep(int *,DMTX *,DMTX *,double *,int);
double Ali2Sco_sq(int *,DMTX *,DMTX *,double *);//New
double AliTbl2Sco(int *,int *,int,DMTX *,DMTX *,double *);
double AliTbl2Sco_sq(int *,int *,int,DMTX *,DMTX *,double *);//New
double BetweenFrag(int,int,int,int,int,DMTX *,DMTX *,double *);//New
int RefineAli(ALI *,DMTX *,DMTX *,double,double,int,PDB *,PDB *);
int IterDP(ALI *,DMTX *,DMTX *,double *,double,double,int,double,double);
//int IterDPfromInitAli(ALI *,DMTX *,DMTX *,double,double);

int IterDPfromInitAli(ALI *,DMTX *,DMTX *,double,double,int,double *,double *,ALIDATA *);
int IterDPfromInitAliSepW(ALI *,DMTX *,DMTX *,double,double,int,double *,double *,ALIDATA *,double);
int IterDPfromInitAliSubOpt(ALI *,DMTX *,DMTX *,double,double,int,double *,double *,ALIDATA *);

int IterDPfromInitAliFrag(ALI *,DMTX *,DMTX *,double,double,int);
double cad_bnd(double,double);
int frag2alitbl(int *,int *,int,int,int *,int *);
//Triangle ver
int Tri(DMTX *,TRIA **,int);
int TriPair(TRIA *,int,TRIA *,int,DMTX *,DMTX *);
int UpperLowerBound(TRIA *,int,TRIA *,int,DMTX *,DMTX *);
int TreeAlign(TRIA *,int,TRIA *,int,DMTX *,DMTX *,ALI *);
int NrAlignments(ALI *);
int NrAlignmentsRate(ALI *,double);
int LocalStrSim(PDB *,PDB *,int,double *);
int LocalStrSimWindow(PDB *,PDB *,int,double *);//Window ver

int ShowResults(ALI *,DMTX *,DMTX *,double,double,int,PDB *,PDB *);
int ShowResultsSepW(ALI *,DMTX *,DMTX *,double,double,int,PDB *,PDB *,double,bool);
int CopyAliSingle(int *,int *,int);
int CopyAli(int *,int *,int *,int *,int,int);
int ShowAli(int *, int);
int ShowPymolScript(int *,int *,int,int);
double CalAliRmsd(int *,int,PDB *,PDB *,int *);
int ShowAliCode(int *,int *,int,int,PDB *,PDB *);
int ShowLocalQ(int *,int *,int,int,DMTX *,DMTX *);
int Malloc_ALI(ALI *);
int LocalPair(int *,double *,DMTX *,DMTX *,int);
int SetFmtx(DMTX *,DMTX *,double ,double *);


//Gaussian Filter 
static int PASCAL_TRI[7][7] = {
{	1, 	0, 	0,	0,	0,	0,	0}, // ala
{ 	2, 	1, 	0,	0,	0,	0,	0}, 
{ 	6, 	4, 	1,	0,	0,	0,	0}, 
{ 	20, 	15, 	6,	1,	0,	0,	0}, 
{ 	70, 	56, 	28,	8,	1,	0,	0},
{ 	252, 	210, 	120,	45,	10,	1,	0}
}; 

static double G_FILTER[7][7] = {
{	1.00, 	0, 	0,	0,	0,	0,	0}, // ala
{ 	0.50, 	0.25, 	0,	0,	0,	0,	0}, 
{ 	0.375, 	0.25, 	0.0625,	0,	0,	0,	0}, 
{ 	0.3125,	0.234375,	0.09375,0.015625,	0,	0,	0}, 
{ 	0.2734375,	0.21875,	0.109375,	0.03125,	0.00390625,	0,	0},
{ 	0.24609375,	0.205078125,	0.1171875,	0.043945313,	0.009765625,	0.000976563,	0},
{	0.225585938,	0.193359375,	0.120849609,	0.053710938,	0.016113281,	0.002929688,	0.000244141 }
}; 
