int chkcmdline(int, char **,CMD *);

int readpdb(PDB *, char *,int);//include malloc
int MallocPdb(PDB *,int);//include malloc
int SetCaCoords(PDB *);//include malloc

void add_phipsi(PDB *);
void add_cbvec(PDB *);
int readgri(GRID *, char *);

int add_SS(PDB *,char *);

//cd_func
void cp_cd(COORD *,double ,double ,double );
void copy_cd(COORD *,COORD *);
void add_cd(COORD *,COORD *,COORD *);
void diff_cd(COORD *,COORD *,COORD *);
void se_cd(COORD *,COORD *,COORD *);
double dis_cd(COORD *,COORD *);
double len_cd(COORD *);
void cl_cd(COORD *);
void cog_cd(COORD *,COORD *,COORD *,COORD *);
double in_pro_cd(COORD *,COORD *);
double inpro(COORD *cd1,COORD *cd2);
void print_cd(COORD *);
void times_cd(COORD *,COORD *,double);
int  misschk_cd(COORD *);
void unit_vec(COORD *);
double vec2area(COORD *,COORD *);
double tortion(COORD *,COORD *,COORD *,COORD *);

double cd_fit(COORD *,COORD *,int *,int,int,int,double, MTX);

//void counter(PDB *,MAT *,CMD *);
void out_pro_cd(COORD *,COORD *,COORD *);//gaiseki

void new_cd_get(COORD *,COORD *,COORD *,COORD *,COORD *,COORD *);//new,X,Y,Z,cog,old

void cl_cmd(CMD *);
int readPROF(PROF *,char *);
int read_eMTX(eMTX,eMTX,char *);
//void COG(PDB *,double []);//center of gravity
void COG(PDB *,COORD *);//center of gravity
void VEC(double [],double [],double []);//vector
void MV(PDB *,COORD *);
void MV_Gri(GRID *,COORD *);
//void TOZ(PDB *,double []);
void TOZ(PDB *,COORD *);
void TOZ_Gri(GRID *,COORD *);

void cp_pdb(PDB *,PDB *);
void cp_Gri(GRID *,GRID *);


double in_pro(double,double,double,double,double,double);/*����*/


//void sch(PDB *,PDB *,MAT *,PDB *,CMD *);/*���������󥸥�*/
void sch(PDB *,PDB *,GRID *,GRID *,PDB *,CMD *);/*���������󥸥�*/
double sco(PDB *);/*�������׻�*/
double Mcounter(PDB *,PDB *,MAT *,CMD *);/*�ޥȥꥯ������μ�ͳ���ͥ륮���׻�*/
double Gcounter(PDB *,PDB *,GRID *,GRID *,CMD *);//get score from grid data

int MVandSCORE(SCO *,PDB *,PDB *,PDB *,GRID *,GRID *,CMD *,MAT *,float,int,int,int,int,int);//��ư�ȥ������׻�

void result(PDB *,char);/*���*/
float vdw_aa(char *);

int Main_coll_chk(PDB *,PDB *,CMD *);//collision check

//random-----------------------
void sch_rand(PDB *,PDB *,GRID *,GRID *,CMD *,DATA *,eMTX,eMTX,PDB *);
double Random();
double range_rand(double,double);
int range_rand_i(double,double);
int flip(double);

void input_data(DATA *,MTX,double);
double local_opti(MTX,MTX,PDB *,PDB *,eMTX,eMTX,CMD *);
void high_reso_doc(PDB *,PDB *,GRID *,GRID *,CMD *,DATA *);
double mtx_ene(PDB *,PDB *,eMTX,eMTX,CMD *);
double mtx_shape(PDB *,PDB *,eMTX,eMTX,CMD *);

double fast_shape_comp(GRID *,GRID *,CMD *);
void trn_ALLgri(GRID *,GRID *,MTX);


//lsfit
double lsfit(COORD *,COORD *,int,MTX);
double lsfit2(COORD *,COORD *,int,MTX);
double fast_lsfit(COORD *,COORD *,int,MTX);
void zero_mtx(MTX );
double rmsd(COORD *,COORD *,int );
void cal_g(COORD *,COORD *,int);
void mlt_vec(MTX ,COORD *,COORD *);
void cp_mtx(MTX ,MTX );
void mlt_mtx(MTX ,MTX ,MTX ,double );
void mlt_mtx2(MTX,MTX,MTX,double,int);
double det(MTX );

#define EPSIRON  1.0E-14
#define EPSIRON2 1.0E-37
//#define EPSIRON2 1.0E-100
//#define EPSIRON  1.0E-50
#define NTMAX    30

double enephobic_ca(PDB *,PDB *,PDB *);

void CACB_trn(PDB *,PDB *,MTX );
void CACen_trn(PDB *,PDB *,MTX );
int CACen_colli(PDB *,PDB *,CMD *);
void show_ca(PDB *);

void cluster(DATA *,PDB *);
void CAcoord_trn(COORD *,PDB *,MTX);

void final_score(DATA *,PDB *,PDB *,GRID *,GRID *,PDB *,eMTX,eMTX,CMD *);

void show_mtx(MTX);

void trn_pdb(PDB *,PDB *,MTX);
void trn_gri(GRID *,GRID *,MTX);
void COORDtrn(COORD *,COORD *,MTX);

int AA2int(char *);
int A2int(char);

double SDcal(double *,int);
double AVEcal(double *,int);

int side_det(PDB *,int);
float next_d(float,float,float,int,int,int);

//int select_tet(int *,GRID *,float, int);
int select_tet(COORD *,PDB *,GRID *,float, int);

double eVF3D(PDB *,PDB *,eVF3D_MTX *,float len,CMD *);

void eVF3D_for_db_make_dim(eVF3D_DATA *,PDB *,PDB *,CMD *);
void eVF3D_for_db_make_mon(eVF3D_DATA *,PDB *,CMD *);

void DVS(COORD *,COORD *);
void vec_separation(COORD *,COORD *,COORD *);
double ave(double *,UINT);
double std(double *,UINT);

void show_node(NODE **,int);
double log_convert(CONT **,double *,int,int *,int,CMD *);
void FreeRotate(COORD *,COORD *,COORD *,COORD *,int,double);

double one2one(CONT **,CONT **,int,int,CMD *,int **);
double one2one_prof(CONT **,CONT **,int,int,CMD *,double*);
double hist_ave(int *,int);
double hist_std(int *,int);

//pdbFFT
int pdbcog(PDB *,COORD *);
int pdbcogFull(PDB *,COORD *);
int shift2cd(COORD,COORD *,PDB *);
int shift2cdFull(COORD,COORD *,PDB *);
double max_dist(COORD,COORD *,int);
void Euler2mtx(double,double, double,MTX);

//dp Sij, gap open,extention,align
double dp(double *,double,double,int *,int,int *,int,int *,int *);
double dp_afp(double *,double,double,int *,int,int *,int,int,int *,int *);
double dp_afp_fwd(double *,double,double,int *,int,int *,int,int,int *,int *,DMTX *,DMTX *,double);
double TraceBack(DMTX *,DMTX *,int,int,int,int,int);

//shib.c
int LowerBond(char **,int,int,int,FTABLE *,int *);
int LowerBond_perstr(CAPDB *,FTABLE *);

int LowerBond_A3(char **,int,int,int,FTABLE *,int *);
int write_sdb(char **,int *,int,int,int,FTABLE *,char *);
int read_sdb(int *,int *,int *,int,FTABLE **,CAPDB *,char *);
int comp_sdb(int,FTABLE *,int,FTABLE *,int ,double,CAPDB *,CAPDB *,SHIBRESULTS *);
int comp_sdb_A3(int,FTABLE *,int,FTABLE *,int ,double,CAPDB *,CAPDB *,SHIBRESULTS *);
int comp_sdb_NV(int,FTABLE *,int,FTABLE *,int ,double,CAPDB *,CAPDB *,SHIBRESULTS *);
int read_ca(CAPDB *,char *);
int read_ca_from_list(CAPDB *,char **,int);
int show_shibresults(int,char **,FTABLE *,SHIBRESULTS);

//clustring
int set_dmtx(CAPDB *,FTABLE *,short int **,int *,int *,int,double);
int clustering_ward(short int **,int *,int *,int ,double,double);
int find_rmodel(short int **,int *,int *,int,int,char **);

double Ave(double *,int);
double Std(double *,int);
