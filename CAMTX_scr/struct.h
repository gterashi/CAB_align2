
#define PI 3.141592
#define ATOM 50000
#define RES 5000
#define AMI 20
#define LIN 256
//grid
#define GRINO 20000


/*def �δĶ�����*/
#define INC 15/*ligand��ž���ٿ�*/
#define MOV 1/*ligand��ư���ٿ�*/
#define START_RT 0 /*ligand��ž��start*/
#define START_MV 0 /*ligand��ư��start*/
#define FIN_RT 360
#define FIN_MV 0
#define ON 1
#define OFF 0
#define TRUE 0
#define FALSE -1
/*#define MATFILE "result030106"*/ /*�������ޥȥꥯ���Υե�����*/
#define MATSIZE 1500
//#define MATLEN 600
#define MATAMI 24 /*�ޥȥꥯ���Υ��ߥλ�����*/

#define CSHNUM 20 /*���ͻĴ��������*/

/*�ޥ�������*/
#define X(a) (a)*(a)
#define L(a,b,c,d,e,f) (sqrt((d-a)*(d-a)+(e-b)*(e-b)+(f-c)*(f-c)))/*�٥��ȥ�Ĺ*/
#define RAS(a) (2.000000*PI*a/360.000000)/*��->�饸����*/
#define RAD(a) (2.000000*PI*a/360.000000)/*��->�饸����*/
#define ANG(a) (360.000000*a/(2.000000*PI))/*�饸����->��*/

#define GAUSS(a,b) 1.00/(sqrt(2.00*PI)*b)*exp(-(a*a)/(2.0*b*b))


/*triangle*/

#define MAXMTX   4

/*��̽���*/
#define TOP 3000

#define NOT_GAUSS 1
#define USE_GAUSS 0

#define VALIABLE 0
#define CONSTANT 1

/*======================*/
/*�����ѿ�*/
/*
int cb,debug,sumresi ;
double in,mid,out,e1,e2,collision ;

int rev_ch;

double enemat_r;

double inc_rt,inc_mv,start_rt,start_mv,fin_rt,fin_mv;
double cog_dis;
double co_k;

double jusin_1[3],jusin_2[3];
double jusin_3[3];
double vector[3];
char chain;
double res_x,res_y,res_z;
double ang_xy,ang_yz,ang_zx;
*/
	/*chain�Ǥ��ֹ��̤�*/
//int FINAL_ATNo;
//int FINAL_RENo;
typedef struct{
	        //double x,y,z;
	        float x,y,z;
}COORD;

//in file structure
//Target & Reference
typedef struct{
	int Num_T,Num_R;
	char ID_T[100],ID_R[100];
	char AA_T[RES*2],AA_R[RES*2];
}IN_DATA;

typedef struct{
	int NumOfRes;
        int AA2int_data[RES];
	double Fp[RES],Fb[RES];
	double phi[RES],psi[RES];
	int SS[RES];
	int hb[RES];
	double score[RES];
	double mesh[RES][RES];
	double area_ex[RES],area_all[RES],area_bur[RES];
	char Chain[ATOM][2],RealNum[RES][5];
}PROF;
/*PDB-file��data��¤*/
typedef struct{
	char fname[LIN];
	int NumOfAtom, NumOfRes; 
	//int ResNnum[RES];
	//int SS[RES];
	//int AA2int_data[RES];
	//int AA2int_data_real[RES];
	//New!! from sakai typeatm.h
	int TypeAtomOder[RES][17];
	//int AtomOnRes[ATOM], SosuiAtom[ATOM], ConservedAtom[ATOM]; 
	//char TypeAtom[ATOM][4], TypeRes[RES][4], Chain[ATOM][2],RealNum[RES][5]; 
	float *Charge; 
	COORD *coord;
	COORD CAcd[RES];
	COORD *CBcd;
	COORD *Cen;
	COORD *Intra;//interaction
	int NumOfIntra;
	COORD *Nonin;//non intra
	int NumOfNonin;
	float *phi,*psi;

	//HETATM
	char **HET_TypeAtom;
	COORD *HET_coord;
	int NumOfHet;
	int NumOfReal;
	int RealResNum[RES];
	//New 2012.10.29-----------------
	double **xyz;
	int *TypeAtomId,*TypeResId;
	char **TypeAtom, **TypeRes, *Chain,**RealNum;
	int *AtomOnRes, *SosuiAtom, *ConservedAtom;
	int *ResOnAtom;
	int *ResNum,*AtomNum;
	double MaxXyz[3],MinXyz[3];
	//For Complex
	int NumOfChain;
	int *Cid;
	int NumOfAtomPerChain[10];
	double **CAxyz;
} PDB;

/*vec_SCOMAT�ι�¤*/
typedef struct{
	char ATOM1[4],ATOM2[4];
        char AA1[MATSIZE][4];
        char AA2[MATSIZE][4];
        float SC[MATSIZE][100];
	float total_sc[MATSIZE];
	float total_AAsc[20];
	float logdata[MATSIZE];
        int NumOfMat;
} MAT;

typedef struct{
	//PDB *model;
	int Rido,Rkeido,Rtr_x,Rtr_y,Rtr_z;
	float RR;
	double Rsco;
	double Re_sco;
	double Rtotal_sco;
	double Sf_sco;
	int clast_flag;
} RESULT;

typedef struct{
	char filename[LIN];
	char list_tg[LIN];
	int rs_mtxmode,ws_mtxmode;
	char s_mtxname[LIN];
	int mtxmode;
	int no_g_mode;
	char mtxname[LIN];
	float len;//fraction distance
	float ang_dis;//phi psi distance
	int chi_dif;//chi difference cut off
	char chk_name[LIN];
	int chk_mode;
	double r,z,m;
	int cut,ext,Ncon;


	int ans;
	int w;
	float relative,Ca_dif;

	//FFT
	int threads;
	int m_mode,rank;
	char p1[LIN]; char p2[LIN];//pdb file name
	char ddi[LIN];
	int ref;
	//ShibyaA1
	char q[LIN];
	char sdb[LIN];
	char list[LIN];
	int LenOfSeg;
	int mkdb_mode;
	double rms,rate;
	bool jmode,fmode,IgBroken,A3mode,NVmode;

	//HPMC
	double grid,dif,ang,std;
	double g_step;
	bool MtxMode,PlyMode,DemoMode,ChkMissMode;
	int ShowResNum,Color;
	double MetaC;
		

} CMD;
typedef struct{ 
	char AAname[4]; 
	float sco; 
} POP;

typedef struct{
 int point[3];
}TRI;

typedef double MTX[MAXMTX][MAXMTX];
typedef float eMTX[20][20];

#define PATCH_GROUP_NUM 10

typedef struct{
	COORD coord[GRINO];
	float sco[GRINO];
	unsigned int NumOfGri;
	unsigned int NumOfSf;
	unsigned int NumOfAllSf;
	COORD sf[GRINO];
	COORD ALLsf[GRINO];
	//int BURI[GRINO];
	COORD bump[GRINO];
	//int patch[GRINO][4];
	//COORD patch_cd[GRINO][4];
	COORD patch_cd[GRINO][3*PATCH_GROUP_NUM];
	int patch_num[GRINO];
	int chk_flag[GRINO];
	int buri[GRINO];
} GRID;

typedef struct{
	double grisco;
	double popsco;
	float R;
} SCO;

//�ǽ��ǡ���

#define DATA_NUM 5000
typedef struct{
	MTX mtx[DATA_NUM];
	double con_e[DATA_NUM];
	double con_Z[DATA_NUM];
	double sha_e[DATA_NUM];
	double sha_Z[DATA_NUM];
	double pop_e[DATA_NUM];
	double pop_Z[DATA_NUM];
	double gri_e[DATA_NUM];
	double gri_Z[DATA_NUM];
	int NumOfData;
	double worst;
	int worst_num;
	int clust_no[DATA_NUM];
	unsigned int NumOfCluster;
	float rmsd[DATA_NUM];
} DATA;

#define LIST_MAX 10000
typedef struct{
 char REC[LIST_MAX][LIN];
 char LIG[LIST_MAX][LIN];
 //float R_AA[LIST_MAX][20];//AA area ratio
 //float L_AA[LIST_MAX][20];//AA area ratio
 unsigned int total_list_no;
} LIST;

#define AA_MAX_NUM 20 //���ߥλ�����
#define BURIED_CLASS 200
#define ENE_CLASS 200 //���ƥ����
#define MAX_ENE 11.00
#define MIN_ENE -10.00
#define TOR_CLS 6

//AAtype: Buried: Hydropathy
typedef struct{
        int b_cls,e_cls,p_cls;
	int fb_flag,fp_flag;
        float max_ene,min_ene;
	float r;
        //float AA_k[AA_MAX_NUM][3];
        float AA_k[AA_MAX_NUM];
        float sco[AA_MAX_NUM][BURIED_CLASS][ENE_CLASS*3];
}eVF3D_MTX;

typedef struct{
        int b_cls,p_cls,e_cls,ss_cls;
	int fb_flag,fp_flag;
        float max_ene,min_ene;
	float max_p,min_p;
	float max_b,min_b;
	float r,ang_r;
	//float AA_k[AA_MAX_NUM][3];
	float AA_k[AA_MAX_NUM*27*16];
        //unsigned int VFdata[AA_MAX_NUM][BURIED_CLASS][ENE_CLASS];
}eVF3D_DATA;


static char *int_to_AA[20]={
"ALA",
"VAL",
"PHE",
"PRO",
"MET",
"ILE",
"LEU",
"ASP", 
"GLU", 
"LYS",
"ARG", 
"SER",        
"THR", 
"TYR", 
"HIS", 
"CYS", 
"ASN", 
"TRP", 
"GLN",
"GLY"
};
#define  MAX_LEN 10000

typedef struct{
        int     chk_len;
        char    seq_AA[MAX_LEN];
        double  chk_mat[MAX_LEN][21];
} CHK;

#define MAX_MESH 15
typedef struct{
        int ID,bur,exp,tot;
	int m,AA,SS;
        int id[MAX_MESH],aa[MAX_MESH],area[MAX_MESH];
	int ss[MAX_MESH];
}MDB;

#define MAX_GLEN 5
#define MAX_FRAG 15
typedef struct{
	//fregment table
	int frag[MAX_FRAG][MAX_GLEN];

        int ID,bur,exp,tot;
        short int m,AA,SS;
        int id[MAX_GLEN*MAX_FRAG],fp[MAX_GLEN*MAX_FRAG][MAX_FRAG];
}GRAPH;

typedef struct{
 //short int aa,ss,bur,exp,area,res,phi,psi;
 short int aa,area,res,phi,psi;
 short int chi[5];
 unsigned int ed;
}NODE;

typedef struct{
 unsigned int id[2];
 unsigned short int area;
 COORD k_1to2,k_2to1;
 float dis;
 short int torsion;
}EDGE;
#define MAX_FRA 3000

typedef unsigned int UINT;
typedef unsigned short USHORT;

#define MAX_SP 200
typedef struct{
 short int aa;
 int id[MAX_SP],edge_num,main_num;
 double *w;
 //double dis[RES];
 double sub;
 int ss;
 COORD cd;
 COORD vec;
}CONT;

typedef struct{
 int a,b,c;//shift data
 double e[3];//Eular
 double s,z;
 MTX mtx;
} SLIST;

//DP pointer
#define NON -1
#define UP 0
#define LEFT 1
#define DIA 2
#define GAP -1

typedef struct{ 
	double sco; 
	int poi;
} DPMTX;
typedef struct{
 int id[2];
 double dist;
} DLIST;


//shib_search
typedef struct{
 int f,f1,f2;//F(U)
 //new!!
 int f15[15];
 int id;//Model ID
 int pos;//position
 bool miss;
} FTABLE;

typedef struct{
 int *pos,*n;
 char *ami;
 float *x,*y,*z;
 //double **cd;
 int Nst;
 int total;
} CAPDB;

typedef struct{
 FTABLE *f;
 double rms;
} SHIBSET;

typedef struct{
 FTABLE ***f;
 double **rms;
 SHIBSET **shib;
 int *n;
} SHIBRESULTS;

