/*
typedef struct{
	int *real;//real number
	double **area,*exp;//order based contact area and exposed area
	int N,Max;
	double Tarea,TareaSq;//Total Area
}DMTX;
*/

int setup_mammoth(double **,int,double **,int,int,double *);
double mammoth(int *,int *,double **,int,double **,int,int,double,double,double *,double *);//a1,a2,CAxyz1,len1,CAxyz2,len22,Nv,Gopen,Gext
double mammoth_subopt(int **,int **,int,int,int,int,double,double,double *,double *);//a1,a2,N,CAxyz1,len1,CAxyz2,len22,Nv,Gopen,Gext
double mammoth_subopt_fwd(int **,int **,int,int,int,int,double,double,double *,double *,DMTX *,DMTX *,double);//a1,a2,N,CAxyz1,len1,CAxyz2,len22,Nv,Gopen,Gext


