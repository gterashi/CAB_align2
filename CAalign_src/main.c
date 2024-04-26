/*
Distance matrix allignment
*/

//#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "struct.h"
#include "func.h"
//#include "rmsd.h"
//#include "mc.h"
#include "dali.h"
#include "tm.h"

#define PDB_STRLEN 55
#define MARGIN 5.0
#define SPHERE_R 1.40
#define FRAG_LEN 6

void malloc_error(char *a){
 fprintf(stderr,"malloc error in %s\n",a);
 exit(0);
}
double gettimeofday_sec() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + (double)tv.tv_usec*1e-6;
}

int readlist(char *fname,char **list){
 int num=0;
 FILE *fp;
 int len;
 if((fp=fopen(fname,"r"))==NULL)
  return FALSE;
 while(fgets(list[num],LIN,fp)!=NULL){
  len=strlen(list[num]);
  list[num][len-1]='\0';//ignore terminal \n
  num++;
 }
 fclose(fp);
 return TRUE;
}

int line_num(char *fname){
 int num=0;
 FILE *fp;
 char line[LIN];
 if((fp=fopen(fname,"r"))==NULL)
  return FALSE;
 while(fgets(line,LIN,fp)!=NULL){
  num++;
 }
 fclose(fp);
 return num;
}

int cmp_ftbl(const void *c1, const void *c2){
 FTABLE a1=*(FTABLE *)c1;
 FTABLE a2=*(FTABLE *)c2;

 return(a1.f-a2.f);
}

int cmp_alidata(const void *c1, const void *c2){
 ALIDATA a1=*(ALIDATA *)c1;
 ALIDATA a2=*(ALIDATA *)c2;
 if(a1.sco>a2.sco)
  return 0;
 //return(a1.sco-a2.sco);
 return 1;
}


double Ave(double *d,int n){
 double sum=0;
 int i;
 for(i=0;i<n;i++){
  sum+=d[i];
 }
 return (sum/n);
}
double Max(double *d,int n){
 double max=-999999999;
 int i;
 for(i=0;i<n;i++){
  if(max<d[i])
   max=d[i];
 }
 return (max);
}


double Std(double *d,int n){
 double sum=0;
 double ave=Ave(d,n);
 int i;
 for(i=0;i<n;i++){
  sum+=(d[i]-ave)*(d[i]-ave);
 }
 return sqrt(sum/n);
}


double *pkdata;//for Zscore
int pk_cnt;

PDB pdb1,pdb2;
CMD cmd;

//for DP
DPMTX *dmtx;

char **dblist;

int count_ca(char **,int);
int CountAtom(char *);
//double mammoth(int *,PDB *,PDB *,int,double,double);//ali,p1,p2,Nv,Gopen,Gext

int main(int argc, char **argv)
{
 //CPU time
 double t1=gettimeofday_sec();
 double t2,t3;
 DMTX d1,d2;
 int i,j,k;
 //--------

 if(chkcmdline(argc,argv,&cmd)==FALSE)
  return(0);

 int LenOfFrag=cmd.len;

 int N1=FindMaxResNum(cmd.d1);
 int N2=FindMaxResNum(cmd.d2);

 //Malloc
 printf("#d1= %s %d\n",cmd.d1,N1);
 if(MallocDmtx(&d1,N1)==FALSE){
  printf("Malloc error in d1\n");
  return 0;
 }
 printf("#d2= %s %d\n",cmd.d2,N2);
 if(MallocDmtx(&d2,N2)==FALSE){
  printf("Malloc error in d2\n");
  return 0;
 }
 //END malloc

 //Read file
 ReadDmtx(&d1,cmd.d1,cmd.IgNeighbor);
 ReadDmtx(&d2,cmd.d2,cmd.IgNeighbor);

 printf("#SepWeight= %.3f\n",cmd.SepWeight);

 d1.Max=N1;
 d2.Max=N2;
 printf("#Tarea1= %.3f\n",d1.Tarea);
 printf("#Tarea1_sep4= %.3f\n",d1.Tarea4);
 if(d1.Tarea4==0)
  d1.Tarea4=0.001;
 printf("#Tarea2= %.3f\n",d2.Tarea);
 printf("#Tarea2_sep4= %.3f\n",d2.Tarea4);

 if(d2.Tarea4==0)
  d2.Tarea4=0.001;

 //reading pdb file
 //count Num of Atom from PDB
 int Natom=CountAtom(cmd.p1);
 if(Natom < 1){
  printf("Cannot find any ATOM recoords in PDB1 %s\n",cmd.p1);
  return 0;
 }
 PDB p1,p2;
 if(MallocPdb(&p1,Natom)==-1)
  return 0;
 if(readpdb(&p1,cmd.p1,Natom)==FALSE){
  printf("##ERROR in readpdb###\n");
  return 0;
 }
 Natom=CountAtom(cmd.p2);
 if(Natom < 1){
  printf("Cannot find any ATOM recoords in PDB2 %s\n",cmd.p2);
  return 0;
 }
 if(MallocPdb(&p2,Natom)==-1)
  return 0;
 if(readpdb(&p2,cmd.p2,Natom)==FALSE){
  printf("##ERROR in readpdb###\n");
  return 0;
 }

 //Set CA coords
 if(SetCaCoords(&p1)==0){
  printf("##ERROR in SetCaCoords###\n");
  return 0;
 }
 if(SetCaCoords(&p2)==0){
  printf("##ERROR in SetCaCoords###\n");
  return 0;
 }

 
/*
 double *LoSim;
 if((LoSim=(double *)malloc(sizeof(double)*p1.NumOfRes*p2.NumOfRes))==NULL)
  malloc_error("Losim");
*/
 //Resnumber check
 if(p1.NumOfRes != d1.Max || p2.NumOfRes != d2.Max){
  printf("Wrong PDB file?? Np1= %d Np2= %d Nd1= %d Nd2= %d\n",p1.NumOfRes,p2.NumOfRes,d1.Max,d2.Max);
  return 0;
 }

 //Generate Seed alignments
 ALI ali;
 double Gopen,Gext,DPsco,Msco=0;
 int DPcnt=0;
 int Nsub=5;//Number of sub optimal alignment 5?

 //param of initial mammoth
 double Go1,Go2,Go3,Ge1,Ge2,Ge3;
/* v05
 Go1=0.00;Go2=10;Go3=1.0;
 Ge1=0.00;Ge2=10.0;Ge3=1.0;
*/

 Go1=0.00;Go2=50.0;Go3=5.0;
 Ge1=0.00;Ge2=10.0;Ge3=2.0;

 //gap less matching trim range
 //int range[6]={p1.NumOfRes,p1.NumOfRes/2,p1.NumOfRes/4,p2.NumOfRes,p2.NumOfRes/2,p2.NumOfRes/4};
 //int range[6]={p1.NumOfRes,p1.NumOfRes/2,p1.NumOfRes/4,p2.NumOfRes,p2.NumOfRes/2,p2.NumOfRes/4};
 //int range[4]={p1.NumOfRes,p1.NumOfRes/2,p2.NumOfRes,p2.NumOfRes/2};
 //int Flen;
 //-----------------------
 ali.len1=d1.Max;
 ali.len2=d2.Max;
 ali.N=0;

 for(Gopen=Go1;Gopen <= Go2;Gopen +=Go3)
  for(Gext=Ge1;Gext <= Ge2;Gext +=Ge3)
   if(Gext<=Gopen)
   ali.N++;
 
 //for tm
 ali.N*=2*Nsub;
 //add Gap less
 //ali.N+=4;

 //malloc ali
 if((Malloc_ALI(&ali))==-1)
  malloc_error("ALIDATA");
 puts("#Fin malloc");

 if((dmtx=(DPMTX *)malloc(sizeof(DPMTX)*(d1.Max+1)*(d2.Max+1)))==NULL)
  return 0;

 double *smtx,*tmp_smtx,*fmtx;
 if((smtx=(double *)malloc(sizeof(double)*(d1.Max+1)*(d2.Max+1)))==NULL)
  return 0;
 
 //for suboptimal DP
 if((tmp_smtx=(double *)malloc(sizeof(double)*(d1.Max+1)*(d2.Max+1)))==NULL)
  return 0;

/*
no use....
 //New Filter Out Mode!!
 double *amtx,*bmtx,fsco;
 int tmpa1[RES],tmpa2[RES],ng,gal[RES*2];
 if(cmd.FilterOut>0){
  t2=gettimeofday_sec();
  printf("#Enter Filter Out Mode TIME= %f\n",t2-t1);
  if((fmtx=(double *)malloc(sizeof(double)*(d1.Max+1)*(d2.Max+1)))==NULL)
  return 0;

  if(SetFmtx(&d1,&d2,cmd.SepWeight,fmtx)==-1)
   return 0;

  fsco=dp(fmtx,0.0,0.0,tmpa1,d1.Max,tmpa2,d2.Max,gal,&ng); 
  printf("Fsco= %f, EstMaxSupRecSco= %f\n",fsco,fsco/pow(d1.Tarea,0.7)+fsco/pow(d2.Tarea,0.7));
  ShowAliCode(tmpa1,tmpa2,d1.Max,d2.Max,&p1,&p2);
  t2=gettimeofday_sec();
  printf("#Finished Filter Out Mode TIME= %f\n",t2-t1);
 }
*/
 if(setup_mammoth(p1.CAxyz,p1.NumOfRes,p2.CAxyz,p2.NumOfRes,FRAG_LEN,smtx)==-1)
  return 0;

 //printf("0:0 = %f 12:13= %f\n",smtx[0],smtx[13+d2.Max*12]);
 //printf("0:0 = %f 13:12= %f\n",smtx[0],smtx[12+d2.Max*13]);


 int **sub_a1,**sub_a2;
 sub_a1=(int **)malloc(sizeof(int *)*Nsub);
 sub_a2=(int **)malloc(sizeof(int *)*Nsub);
 for(i=0;i<Nsub;i++){
  sub_a1[i]=(int *)malloc(sizeof(int )*p1.NumOfRes);
  sub_a2[i]=(int *)malloc(sizeof(int )*p2.NumOfRes);
 }
 //LocalPair(sub_a1,smtx,&d1,&d2,FRAG_LEN);

 //int Nseed;
 DPcnt=0;
 for(Gopen=Go1;Gopen <= Go2;Gopen +=Go3){
  for(Gext=Ge1;Gext <= Ge2;Gext +=Ge3){
   if(Gext > Gopen)
    continue;
   //fragment based
   //printf("#Go= %.2f Ge= %.2f\n",Gopen,Gext);
   	//Copy smtx
	for(i=0;i<(d1.Max)*(d2.Max);i++)
	 tmp_smtx[i]=smtx[i];
	//---------

   	//mammoth_subopt(ali.ali[DPcnt].a1,ali.ali[DPcnt].a2,p1.CAxyz,p1.NumOfRes,p2.CAxyz,p2.NumOfRes,FRAG_LEN,Gopen,Gext,tmp_smtx,&DPsco);//ali,p1,p2,Nv,Gopen,Gext
   	mammoth_subopt(sub_a1,sub_a2,Nsub,p1.NumOfRes,p2.NumOfRes,FRAG_LEN,Gopen,Gext,tmp_smtx,&DPsco);//ali,p1,p2,Nv,Gopen,Gext
   	//mammoth_subopt_fwd(sub_a1,sub_a2,Nsub,p1.NumOfRes,p2.NumOfRes,FRAG_LEN,Gopen,Gext,tmp_smtx,&DPsco,&d1,&d2,1.0);//ali,p1,p2,Nv,Gopen,Gext

	//Copy Alignment
	for(j=0;j<Nsub;j++){
   	 //ShowAliCode(sub_a1[j],sub_a2[j],ali.len1,ali.len2,&p1,&p2);
	 CopyAli(sub_a1[j],sub_a2[j],ali.ali[DPcnt].a1,ali.ali[DPcnt].a2,p1.NumOfRes,p2.NumOfRes);
	 DPcnt++;
	}
	//puts("//");

   //printf("Mammoth %f\n",DPsco);

   //ShowAli(ali.ali[DPcnt].a1,p1.NumOfRes);
   //ShowAliCode(ali.ali[DPcnt].a1,ali.ali[DPcnt].a2,ali.len1,ali.len2,&p1,&p2);
   //puts("//");
   //ShowAli(ali.ali[DPcnt].a1,p1.NumOfRes);
   //ali.ali[DPcnt].sco= DPsco;
   //DPcnt++;
  }
 }
 //return 0;

 double *stbl;
 ali.N=DPcnt;

 printf("# of Seeds = %d\n",ali.N);
 //NR alignments
 NrAlignments(&ali);
 //NrAlignmentsRate(&ali,0.90);
 printf("# of Nr ALI = %d\n",ali.N);

 //Add TMalign like alignment
 for(i=0;i<ali.N;i++){
  //printf("Impose %d\n",i);
  //ShowAli(ali.ali[i].a1,p1.NumOfRes);
  DPsco=tm(ali.ali[i].a1,ali.ali[i].a2,ali.ali[i+ali.N].a1,ali.ali[i+ali.N].a2,p1.CAxyz,p1.NumOfRes,p2.CAxyz,p2.NumOfRes,0.6,0,smtx);
  //printf("mammoth DPscore= %f TMscore= %f\n",ali.ali[i].sco,DPsco);
  //ShowAliCode(ali.ali[i].a1,ali.ali[i].a2,ali.len1,ali.len2,&p1,&p2);
  //ShowAliCode(ali.ali[i+ali.N].a1,ali.ali[i+ali.N].a2,ali.len1,ali.len2,&p1,&p2);
  //puts("//");
  //ShowAli(ali.ali[i].a1,p1.NumOfRes);
  //printf("DPsco= %.3f\n",DPsco);
 }
 ali.N*=2;


/*
 //Gap less mode
 int tmp_ali1[RES],tmp_ali2[RES];
 for(i=0;i<4;i++){
  Flen=range[i];
  if(Flen<7 || Flen >p1.NumOfRes||Flen >p2.NumOfRes)
   continue;
  //find max tm gapless alignment
  //printf("#Gap less p1 %d ->p2\n",Flen);
  GapLessAlign(tmp_ali1,tmp_ali2,p1.CAxyz,p1.NumOfRes,p2.CAxyz,p2.NumOfRes,Flen);
  //ShowAliCode(tmp_ali1,tmp_ali2,ali.len1,ali.len2,&p1,&p2);
  //TMalign
  DPsco=tm(tmp_ali1,tmp_ali2,ali.ali[ali.N].a1,ali.ali[ali.N].a2,p1.CAxyz,p1.NumOfRes,p2.CAxyz,p2.NumOfRes,0.6,0,smtx);
  printf("#Gap less p1 %d ->p2 %f\n",Flen,DPsco);
  ShowAliCode(ali.ali[ali.N].a1,ali.ali[ali.N].a2,ali.len1,ali.len2,&p1,&p2);
  ali.N++;
 }
*/


 qsort(ali.ali,ali.N,sizeof(ALIDATA),cmp_alidata);

 double max=ali.ali[0].sco;
 //printf("Ave= %f Std= %f Max %.3f MaxZ= %.3f N= %d\n",ali.ave,ali.std,max,(max-ali.ave)/ali.std,ali.N);
 t2=gettimeofday_sec();
 printf("#INITIAL ALIGNMENTS TIME= %f\n",t2-t1);
 printf("# of NR = %d\n",ali.N);

 //fragment based
 //IterDPfromInitAliFrag(&ali,&d1,&d2,cmd.GapOpen,cmd.GapExt,FRAG_LEN+1);
 //return 0;
 //residue based

 int go=cmd.GapOpen;

/*
 while(go>0){
  if(ali.N<1)
   ali.N=1;
  if(go<10)
   go=0;
  IterDPfromInitAli(&ali,&d1,&d2,go,cmd.GapExt);
  printf("finished DP N= %d\n",ali.N);
  NrAlignments(&ali);
  qsort(ali.ali,ali.N,sizeof(ALIDATA),cmp_alidata);
  printf("Final N= %d Max= %.3f\n",ali.N,ali.ali[0].sco);
  if(go==0)
   break;
  go/=2.00;
  ali.N*=0.5;
 }
*/
  //IterDPfromInitAli(&ali,&d1,&d2,go,cmd.GapExt,3);
 NrAlignments(&ali);
 printf("#finished TMalign N= %d\n",ali.N);

 //Setup stock
 ALIDATA *stock;
 if((stock=(ALIDATA *)malloc(sizeof(ALIDATA)*ali.N*10))==NULL)
  return -1;
 for(i=0;i<ali.N*10;i++)
  if((stock[i].a1=(int *)malloc(sizeof(int)*p1.NumOfRes))==NULL)
   return -1;


 int window=5;
 while(window>=1){
  if(ali.N<1)
   ali.N=1;
  //IterDPfromInitAli(&ali,&d1,&d2,go*(2.00*window+1.00),0.00,window,smtx,tmp_smtx,stock);
  IterDPfromInitAliSepW(&ali,&d1,&d2,go*(2.00*window+1.00),0.00,window,smtx,tmp_smtx,stock,cmd.SepWeight);
  qsort(ali.ali,ali.N,sizeof(ALIDATA),cmp_alidata);
  NrAlignments(&ali);
  //NrAlignmentsRate(&ali,0.90);

  printf("#Window %d finished DP N= %d\n",window*2+1,ali.N);
  printf("#Final N= %d Max= %.3f\n",ali.N,ali.ali[0].sco);

  //ShowAliCode(ali.ali[0].a1,ali.ali[0].a2,ali.len1,ali.len2,&p1,&p2);
  //ShowAliCode(ali.ali[1].a1,ali.ali[1].a2,ali.len1,ali.len2,&p1,&p2);
  window-=2;
  //ali.N*=0.5;
 }

/*
 for(i=0;i<ali.N;i++){
  DPsco=tm(ali.ali[i].a1,ali.ali[i].a2,ali.ali[i+ali.N].a1,ali.ali[i+ali.N].a2,p1.CAxyz,p1.NumOfRes,p2.CAxyz,p2.NumOfRes,0.6,0,smtx);
 }
 ali.N*=2;
 NrAlignments(&ali);
 IterDPfromInitAli(&ali,&d1,&d2,go,cmd.GapExt,window,smtx,tmp_smtx,stock);
 NrAlignments(&ali);
*/
 //Show Results
 //ShowResults(&ali,&d1,&d2,cmd.zmc,cmd.ExtRate,LenOfFrag,&p1,&p2);
 ShowResultsSepW(&ali,&d1,&d2,cmd.zmc,cmd.ExtRate,LenOfFrag,&p1,&p2,cmd.SepWeight,cmd.LQmode);

 t2=gettimeofday_sec();
 printf("#FINISHED TOTAL TIME= %f\n",t2-t1);
 return 0;

}



