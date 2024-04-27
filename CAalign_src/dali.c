#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "struct.h"
#include "func.h"
#include "dali.h"
#include "rmsd.h"
#define DIST(a,b,c,d,e,f) sqrt((a-b)*(a-b)+(c-d)*(c-d)+(e-f)*(e-f))
#define L(a,b,c,d,e,f) (sqrt((d-a)*(d-a)+(e-b)*(e-b)+(f-c)*(f-c)))
#define SEP_LEN 5

/*
//Random from cd_func.c
double Random(){//0-1 random 
        double r;
        int i;
        i=random();
        if(i != 0)
                i--;
        //r=(double)i/0x7fff; 
        r=(double)i/RAND_MAX;
        return(r);
}

int range_rand_i(double min,double max){ //min-max random
        double Random();
        return((max-min)*Random()+min);
}

double range_rand(double min,double max){ //min-max random
        double Random();
        return((max-min)*Random()+min);
}

int flip(double probability){ //Íð¿ô¤Ë¤è¤ë³ÎÎ¨È½Äê
        double Random();
        if(probability == 1.0)
                return(0);
        else if(Random() <= probability)
                return(0);
        else
                return(-1);
}
*/


int FindMaxResNum(char *filename){
 FILE *fpin; 
 char line[LIN], buf[LIN]; 
	
 if((fpin=fopen(filename,"r")) == NULL){ 
	fprintf(stderr,"Can't open %s\n",filename); 
	return(-1); 
 }

 int i1,i2,i3,i4;
 int max=0;
 float a;
 while(fgets(line,LIN,fpin)){
  if(!strncmp(line,"#",1))
   continue;
  //printf("%s",line);
  sscanf(line,"%d\t%d\t%d\t%d\t%f\n",&i1,&i2,&i3,&i4,&a);
  if(max<i1)//order base
   max=i1;
  if(max<i3)
   max=i3;
 }
 fclose(fpin);
 return (max+1);
}

int MallocDmtx(DMTX *d,int n){

 d->N=n;

 if((d->area=(double **)malloc(sizeof(double *)*n))==NULL)
  return FALSE;

 int i;//2D matrix
 for(i=0;i<n;i++)
  if((d->area[i]=(double *)calloc(n,sizeof(double)))==NULL)
   return FALSE;

 if((d->exp=(double *)malloc(sizeof(double)*n))==NULL)
  return FALSE;
 if((d->real=(int *)malloc(sizeof(int)*n))==NULL)
  return FALSE;

 return TRUE;
}

int ReadDmtx(DMTX *d,char *filename, bool flag){
 FILE *fpin; 
 char line[LIN], buf[LIN]; 
	
 if((fpin=fopen(filename,"r")) == NULL){ 
	fprintf(stderr,"Can't open %s\n",filename); 
	return(FALSE); 
 }

 int i1,i2,i3,i4;
 float a;
 d->Tarea=0;//init
 d->TareaSq=0;//init
 d->Tarea4=0;//init
 while(fgets(line,LIN,fpin)){
  if(!strncmp(line,"#",1))
   continue;
  //printf("%s",line);
  sscanf(line,"%d\t%d\t%d\t%d\t%f\n",&i1,&i2,&i3,&i4,&a);
  if(i3<0){//exposed
   d->exp[i1]=a;
  }else{
 	//if(flag==true && abs(i2-i4)<=1)//ignore next residue
 	if(flag==true && abs(i2-i4)<=1)//ignore next residue
	 continue;
   d->area[i1][i3]=a;
   d->real[i1]=i2;
   d->Tarea+=a;
   if(abs(i2-i4)>4)
    d->Tarea4+=a;
   //d->TareaSq+=sqrt(a); nouse
  }
 }
 fclose(fpin);
 return TRUE;
}

double AveArea(DMTX *d,int i,int j,int len){
 double sum=0;
 int i1,i2;
 for(i1=i;i1-i<len;i1++){
 	for(i2=j;i2-j<len;i2++){
	 sum+=d->area[i1][i2];
	 sum+=d->area[i2][i1];
	}
 }
 return(sum/(len+len));
}

//double PairScore(DMTX *d1,int p1,int p2,DMTX *d2,int q1,int q2,int len,double *b,double *t){
double SingleScore(DMTX *d1,int p1,DMTX *d2,int q1,int len){
 int i1,i2;
 int a1,a2,b1,b2;
 double sum_b,sum_t;
 sum_b=sum_t=0;
 for(i1=0;i1<len;i1++){
  a1=p1+i1;
  b1=q1+i1;
	for(i2=0;i2<len;i2++){
	 if(i1==i2)
	  continue;
	 a2=p1+i2; b2=q1+i2;
	 sum_b+=cad_bnd(d1->area[a1][a2],d2->area[b1][b2]);
	 sum_t+=d1->area[a1][a2];
	}
 }
 //return(sum_t-sum_b);
 return(1-sum_b/sum_t);
}

double PairScore(DMTX *d1,int p1,int p2,DMTX *d2,int q1,int q2,int len){
 double sum=0;
 int i1,i2;
 int a1,a2,b1,b2;
 double M,bounded;
 double sum_b,sum_t;
 sum_b=sum_t=0;
 double bmd;
/*
|p1|-----|p2|
|q1|-----|q2|
*/
 //p1,q1
 for(i1=0;i1<len;i1++){
  a1=p1+i1;
  b1=q1+i1;
	for(i2=0;i2<len;i2++){
	 if(i1==i2)
	  continue;
	 a2=p1+i2; b2=q1+i2;
	 if(d1->area[a1][a2]==0||d1->area[a2][a1]==0||d2->area[b1][b2]==0||d2->area[b2][b1]==0)
	  continue;

	 sum_b+=cad_bnd(d1->area[a1][a2],d2->area[b1][b2]);
	 //sum_b+=cad_bnd(d1->area[a2][a1],d2->area[b2][b1]);

/*
	 M =fabs(d1->area[a1][a2] - d2->area[b1][b2]);
	 bounded=d1->area[a1][a2];
	 if(bounded > M)
	  bounded=M;

	 sum_b+=bounded;
*/
	 sum_t+=d1->area[a1][a2];
	}
 }
 //p2,q2
 for(i1=0;i1<len;i1++){
  a1=p2+i1;
  b1=q2+i1;
	for(i2=0;i2<len;i2++){
	 if(i1==i2)
	  continue;
	 a2=p2+i2; b2=q2+i2;
	 if(d1->area[a1][a2]==0||d1->area[a2][a1]==0||d2->area[b1][b2]==0||d2->area[b2][b1]==0)
	  continue;

	 sum_b+=cad_bnd(d1->area[a1][a2] , d2->area[b1][b2]);
	 //sum_b+=cad_bnd(d1->area[a2][a1] , d2->area[b2][b1]);
/*
	 M =fabs(d1->area[a1][a2] - d2->area[b1][b2]);
	 bounded=d1->area[a1][a2];
	 if(bounded > M)
	  bounded=M;

	 sum_b+=bounded;
*/
	 sum_t+=d1->area[a1][a2];
	}
 }

 
 for(i1=0;i1<len;i1++){
  a1=p1+i1;
  b1=q1+i1;
 	for(i2=0;i2<len;i2++){
	 a2=p2+i2; b2=q2+i2;
	 if(d1->area[a1][a2]==0||d1->area[a2][a1]==0||d2->area[b1][b2]==0||d2->area[b2][b1]==0)
	  continue;

	 //a1-----  ->  a2------
	 //b1-----  ->  b2------
	 //scoring function
	 //CAD score need modefied
/*
	 M =fabs(d1->area[a1][a2] - d2->area[b1][b2]);
	 bounded=d1->area[a1][a2];
	 if(bounded > M)
	  bounded=M;

	 sum_b+=bounded;
*/
	 sum_b+cad_bnd(d1->area[a1][a2],d2->area[b1][b2]);
	 sum_b+cad_bnd(d1->area[a2][a1],d2->area[b2][b1]);
	 sum_t+=d1->area[a1][a2];
/*
	 M =fabs(d1->area[a2][a1] - d2->area[b2][b1]);
	 bounded=d1->area[a2][a1];
	 if(bounded > M)
	  bounded=M;

	 sum_b+=bounded;
*/
	 sum_t+=d1->area[a2][a1];
	}
 }
 //double ans=1-sum_b/sum_t;
 //if(sum_t==0)
 // ans=0;

 //*b=sum_b;
 //*t=sum_t;
 //return(ans);
 //return(sum_t-sum_b);
 return(1.00-sum_b/sum_t);
}

int CalPairlist(PLIST *p,DMTX *d1,DMTX *d2,double *sim,int len,double acut,double rcut,double ecut,bool cflag){
 int i1,i2;
 double ave;
 int *list1,*list2;
 double *sdata;//singlet list
 int n1,n2;
 double sco,rate;
 int a1[1000],a2[1000],al1[1000],al2[1000];
 int acnt=0;
 n1=n2=0;

 //if((sdata=(double *)malloc(sizeof(double)*(d1->N)*(d2->N)))==NULL)
 ///if((p->single_sco=(double *)malloc(sizeof(double)*(d1->N)*(d2->N)))==NULL)
 // return FALSE;

//check local similarity
/*
 for(i1=0;i1<d1->N-len;i1++){
  for(i2=0;i2<d2->N-len;i2++){
   //sdata[i1*(d2->N)+i2]=SingleScore(d1,i1,d2,i2,len);
   p->single_sco[i1*(d2->N)+i2]=SingleScore(d1,i1,d2,i2,len);
   //printf("%d %d %f\n",i1,i2,p->single_sco[i1*(d2->N)+i2]);
  }
 }
*/

 if((list1=(int *)malloc(sizeof(int)*(d1->N-len*2)*(d1->N-len)*2))==NULL)
  return FALSE;
 if((list2=(int *)malloc(sizeof(int)*(d2->N-len*2)*(d2->N-len)*2))==NULL)
  return FALSE;

 for(i1=0;i1<d1->N-len*2;i1++){
 	for(i2=i1+len;i2<d1->N-len;i2++){
	 ave=AveArea(d1,i1,i2,len);
	 //remove no contact data
	 if(ave<acut)
	  continue;
	 //printf("%d-%d (%f)\n",i1,i2,ave);
	 list1[n1*2]=i1;
	 list1[n1*2+1]=i2;
	 n1++;
	}
 }
 for(i1=0;i1<d2->N-len*2;i1++){
 	for(i2=i1+len;i2<d2->N-len;i2++){
	 ave=AveArea(d2,i1,i2,len);
	 //remove no contact data
	 if(ave<acut)
	  continue;
	 //printf("%d-%d (%f)\n",i1,i2,ave);
	 list2[n2*2]=i1;
	 list2[n2*2+1]=i2;
	 n2++;
	}
 }
 //malloc pair list
 if((p->id=(unsigned int *)malloc(sizeof(unsigned int)*n1*n2*4*2))==NULL)
  return FALSE;
 if((p->sco=(double *)malloc(sizeof(double)*n1*n2*2))==NULL)
  return FALSE;

 printf("#%d vs %d *2\n",n1,n2);

 int num=0;
 unsigned int tag;
 double bnd,t;
 double tmp=0;
 for(i1=0;i1<n1;i1++){
 	for(i2=0;i2<n2;i2++){
	 //cal score
	 //if( p->single_sco[list1[2*i1]*d2->N +list2[2*i2]] > rcut 
	 //&& p->single_sco[list1[2*i1+1]*d2->N +list2[2*i2+1]]>rcut){
	 if( sim[list1[2*i1]*d2->N +list2[2*i2]] > rcut 
	 && sim[list1[2*i1+1]*d2->N +list2[2*i2+1]]>rcut){
/*
	 a1[0]=list1[2*i1];
	 a1[1]=list1[2*i1+1];
	 a2[0]=list2[2*i2];
	 a2[1]=list2[2*i2+1];
*/
	 //printf(">>%d %d %d %d\n",a1[0],a1[1],a2[0],a2[1]);
	 //just between 2 fragments?????
	 sco=BetweenFrag(list1[2*i1],list1[2*i1+1],list2[2*i2],list2[2*i2+1],len,d1,d2,&rate);
	 //frag2alitbl(a1,a2,2,len,al1,al2);
	 //sco=AliTbl2Sco(al1,al2,len*2,d1,d2,&rate);
	 //sco=AliTbl2Sco_sq(al1,al2,len*2,d1,d2,&rate);
	 	if(rate > ecut){
	 	 p->id[num*4]  =list1[2*i1];
	 	 p->id[num*4+1]=list1[2*i1+1];
	 	 p->id[num*4+2]=list2[2*i2];
	 	 p->id[num*4+3]=list2[2*i2+1];
	 	 p->sco[num]=rate;
		 //tmp+=rate;
	 	 num++;
	 	}
	 }
	}
 }
 for(i1=0;cflag==true && i1<n1;i1++){
 	for(i2=0;i2<n2;i2++){
	 //crossover
	 if(p->single_sco[list1[2*i1]*d2->N +list2[2*i2+1]]>rcut
	 && p->single_sco[list1[2*i1+1]*d2->N +list2[2*i2]]>rcut){
/*
	 a1[0]=list1[2*i1];
	 a1[1]=list1[2*i1+1];
	 a2[0]=list2[2*i2+1];
	 a2[1]=list2[2*i2];
*/
	 //frag2alitbl(a1,a2,2,len,al1,al2);
	 //sco=AliTbl2Sco(al1,al2,len*2,d1,d2,&rate);
	 //sco=AliTbl2Sco_sq(al1,al2,len*2,d1,d2,&rate);
	 sco=BetweenFrag(list1[2*i1],list1[2*i1+1],list2[2*i2+1],list2[2*i2],len,d1,d2,&rate);

	 	if(rate > ecut){//extended rate
	 	 p->id[num*4]=list1[2*i1];
	 	 p->id[num*4+1]=list1[2*i1+1];
	 	 p->id[num*4+2]=list2[2*i2+1];
	 	 p->id[num*4+3]=list2[2*i2];
	 	 //p->sco[num]=sco;
	 	 p->sco[num]=rate;
		 //tmp+=rate;
	 	 num++;
	 	}
   	 }
	}
 }
 //printf("tmp= %f\n",tmp);
 free(list1);
 free(list2);
 p->N=num;
 return num;
}

int GeneSeeds(ALI *ali,PLIST *p,DMTX *d1,DMTX *d2,int len,double cut,bool vflag){
 int i1,i2,i3;
 int n=0;
 int id1,id2,Maxid=0;
 int *CntTbl;

/*
A-------B
|       |
a-------b
singlet (A,a) and (B,b)
*/

 //make singlet table
 if((CntTbl=(int *)calloc((d1->Max*d2->Max),sizeof(int)))==NULL)
  return FALSE;
 //consider All singlet
 for(i1=0;i1<p->N;i1++){
  //if(p->sco[i1]<cut)
  // continue;
  //printf("Add %d - %d & %d - %d %.3f\n",p->id[i1*4],p->id[i1*4+1],p->id[i1*4+2],p->id[i1*4+3],p->sco[i1]);

  id1=d2->Max * p->id[i1*4] + p->id[i1*4+2];//A-a
  id2=d2->Max * p->id[i1*4+1] + p->id[i1*4+3];//B-b

  CntTbl[id1]++;
  CntTbl[id2]++;
  //printf("ID1= %d/%d ID2= %d/%d\n",id1,CntTbl[id1],id2,CntTbl[id2]);
  if(id1>Maxid)
   Maxid=id1;
  if(id2>Maxid)
   Maxid=id2;
   
  n++;
 }
 printf("#Passed pair= %d Maxid= %d\n",n,Maxid);

 //Malloc Pair index
 int **Ptbl;
 Maxid++;
 if((Ptbl=(int **)malloc(sizeof(int*)*Maxid))==NULL)
  return FALSE;


/*
 int *SingleetTbl;
 if((SingletTbl=(int *)malloc(sizeof(int)*Maxid))==NULL)
  return FALSE;
*/

 for(i1=0;i1<Maxid;i1++){
  if(CntTbl[i1]<2)
   continue;
  if((Ptbl[i1]=(int *)malloc(sizeof(int)*CntTbl[i1]*2))==NULL)
  return FALSE;
 }

 int *CntTbl2;
 if((CntTbl2=(int *)calloc((d1->Max*d2->Max),sizeof(int)))==NULL)
  return FALSE;

 //input again
 for(i1=0;i1<p->N;i1++){
  //if(p->sco[i1]<ecut)
  // continue;

  id1=d2->Max * p->id[i1*4] + p->id[i1*4+2];
  id2=d2->Max * p->id[i1*4+1] + p->id[i1*4+3];

  //printf("Add %d - %d & %d - %d %.3f\n",p->id[i1*4],p->id[i1*4+1],p->id[i1*4+2],p->id[i1*4+3],p->sco[i1]);
  //printf("id1= %d(%d) id2= %d(%d)\n",id1,CntTbl[id1],id2,CntTbl[id2]);

  if(CntTbl[id1]>1){
   Ptbl[id1][CntTbl2[id1]*2  ]=p->id[i1*4+1];
   Ptbl[id1][CntTbl2[id1]*2+1]=p->id[i1*4+3];
  }
  if(CntTbl[id2]>1){
   Ptbl[id2][CntTbl2[id2]*2  ]=p->id[i1*4];
   Ptbl[id2][CntTbl2[id2]*2+1]=p->id[i1*4+2];
  }
  CntTbl2[id1]++;
  CntTbl2[id2]++;
 }
 //generate initial path
 //foreach singlet....
 int Ntri=0;
 int TNtri=0;
/*
 if((ali->a1=(int **)malloc(sizeof(int*)*d1->Max*d2->Max))==NULL)
  return 0;
 if((ali->a2=(int **)malloc(sizeof(int*)*d1->Max*d2->Max))==NULL)
  return 0;
 if((ali->sco=(double *)malloc(sizeof(double)*d1->Max*d2->Max))==NULL)
  return 0;
*/
//new ALIDATA
 if((ali->ali=(ALIDATA *)malloc(sizeof(ALIDATA)*d1->Max*d2->Max))==NULL)
  return 0;
 ali->N=0;

 for(i1=0;i1<d1->Max;i1++){
  for(i2=0;i2<d2->Max;i2++){
   id1=d2->Max * i1 + i2;
   if(CntTbl[id1]<2)
    continue;
   //Seed.......
   Ntri=GenePath(ali,i1,i2,Ptbl,CntTbl,len,d1,d2,cut,vflag);
   //printf("Start singlet %d:%d N=%d Ntri= %d\n",i1,i2,CntTbl[id1],Ntri);
   TNtri+=Ntri;
  }
 }
 printf("Total Num of Path = %d\n",TNtri);
 //ali->N=TNtri;
 return ali->N;
}

int check_frag_ali(int *a1,int *a2,int n,int id1,int id2,int len,bool flag){
 int i;
 for(i=0;i<n;i++){
  if(abs(a1[i]-id1)<len||abs(a2[i]-id2)<len)
   return 1;
  if(flag==false && (a1[i]-id1)*(a2[i]-id2)<0)
   return 1;
 }
 return 0;
}

int frag2ali(int *a1,int *a2,int n,int len,int *al1,int *al2,int max1,int max2){
 int i,j;
 //init
 for(i=0;i<max1;i++) al1[i]=-1;
 for(i=0;i<max2;i++) al2[i]=-1;

 for(i=0;i<n;i++){
	for(j=0;j<len;j++){
	 al1[a1[i]+j]=a2[i]+j;
	 al2[a2[i]+j]=a1[i]+j;
	}
 }
 return 0;
}
int frag2alitbl(int *a1,int *a2,int n,int len,int *al1,int *al2){
 int i,j;
 int cnt=0;
 for(i=0;i<n;i++){
	for(j=0;j<len;j++){
	 al1[cnt]=a1[i]+j;
	 al2[cnt]=a2[i]+j;
	 cnt++;
	}
 }
 return cnt;
}
double frag2newsco(int *a1,int *a2,int n,int id1,int id2,int len,DMTX *d1,DMTX *d2,double *rate){
 int i,j,k,l;
 int p11,p12,p21,p22;
 double bnd=0;
 double t=0;
 //p11------------p12
 //p21------------p22

 for(i=0;i<n;i++){
 	for(j=0;j<len;j++){
	 p11=a1[i]+j;
	 p21=a2[i]+j;
		for(k=0;k<len;k++){
		 p12=id1+k;
		 p22=id2+k;
		 bnd+=cad_bnd(d1->area[p11][p12],d2->area[p21][p22]);
		 bnd+=cad_bnd(d1->area[p12][p11],d2->area[p22][p21]);
		 t+=d1->area[p11][p12]+d1->area[p12][p11];
		}
 	}
 }
 //itself important!!!
 for(j=0;j<len;j++){
  p11=id1+j;
  p21=id2+j;
	for(k=j+1;k<len;k++){
	 p12=id1+k;
	 p22=id2+k;
	 bnd+=cad_bnd(d1->area[p11][p12],d2->area[p21][p22]);
	 bnd+=cad_bnd(d1->area[p12][p11],d2->area[p22][p21]);
	 t+=d1->area[p11][p12]+d1->area[p12][p11];
	}
 }

 *rate=1-bnd/t;
 return (t-bnd);
}

double frag2newsco_sq(int *a1,int *a2,int n,int id1,int id2,int len,DMTX *d1,DMTX *d2,double *rate){
 int i,j,k,l;
 int p11,p12,p21,p22;
 double bnd=0;
 double t=0;
 double sco,sq_t;
 //p11------------p12
 //p21------------p22
 sco=0;
 for(i=0;i<n;i++){
 	for(j=0;j<len;j++){
	 p11=a1[i]+j;
	 p21=a2[i]+j;
		for(k=0;k<len;k++){
		 p12=id1+k;
		 p22=id2+k;
		 bnd=cad_bnd(d1->area[p11][p12],d2->area[p21][p22]);
		 bnd+=cad_bnd(d1->area[p12][p11],d2->area[p22][p21]);
		 sq_t=d1->area[p11][p12]+d1->area[p12][p11];
		 sco+=sqrt(sq_t-bnd);
		 t+=sqrt(sq_t);
		}
 	}
 }
 //itself important!!!
 for(j=0;j<len;j++){
  p11=id1+j;
  p21=id2+j;
	for(k=j+1;k<len;k++){
	 p12=id1+k;
	 p22=id2+k;
	 bnd=cad_bnd(d1->area[p11][p12],d2->area[p21][p22]);
	 bnd+=cad_bnd(d1->area[p12][p11],d2->area[p22][p21]);
	 sq_t=d1->area[p11][p12]+d1->area[p12][p11];
	 sco+=sqrt(sq_t-bnd);
	 t+=sqrt(sq_t);
	}
 }

 *rate=sco/t;
 return (sco);
}

int trimfrag(int *in1,int *in2,int n,int *out1,int *out2,int pos){
 int i;
 for(i=0;i<n;i++){
  if(i<pos){
   out1[i]=in1[i];
   out2[i]=in2[i];
  }else{//shift
   out1[i]=in1[i+1];
   out2[i]=in2[i+1];
  }
 }
 return 0;
}

int GenePath(ALI *a,int id11,int id12,int **tbl,int *cnt_tbl,int len,DMTX *d1,DMTX *d2,double rcut,bool vflag){
 int n1=d1->Max;
 int n2=d2->Max;
 int i1,i2;
 int id21,id22,id31,id32;
 int cnt=0;
 int a1[1000],a2[1000];//fragment
 int ali1[1000],ali2[1000];//all alignment
 int tmpa1[1000],tmpa2[1000];//all alignment
 int acnt=0;
 double sco;
 int id,maxid[2];
 double max=0;
 double pre=0;
 double rate,trate,brate,bsco,crate;
 int bestt,backflag=0;
 //start
 a1[acnt]=id11;
 a2[acnt]=id12;
 acnt++;

 //find next fragment best score....
 while(1){
  max=0;
	//add new regions
 	for(i1=0;i1<acnt;i1++){
 	 id=n2*a1[i1]+a2[i1];
 	 if(cnt_tbl[id]<2)
 	  continue;
  		for(i2=0;i2<cnt_tbl[id];i2++){
		 //check available alignment?
		 if(check_frag_ali(a1,a2,acnt,tbl[id][2*i2],tbl[id][2*i2+1],len,vflag)!=0)
		  continue;
		 //fast cal frag -> score
		 sco=frag2newsco(a1,a2,acnt,tbl[id][2*i2],tbl[id][2*i2+1],len,d1,d2,&crate);
		 //sco=frag2newsco_sq(a1,a2,acnt,tbl[id][2*i2],tbl[id][2*i2+1],len,d1,d2,&crate);
		 if(crate<rcut)
		  continue;
		 //printf("cand %d:%d sco= %.3f rate= %.3f\n",tbl[id][2*i2],tbl[id][2*i2+1],sco,crate);
		 //if(max<crate){
		 if(max<sco){
		  max=sco;
		  //max=crate;
		  maxid[0]=tbl[id][2*i2];
		  maxid[1]=tbl[id][2*i2+1];
		 }
		}
 	}
 //no improvement
 if(max==0)
  break;
  //rate check
 a1[acnt]=maxid[0];
 a2[acnt]=maxid[1];
 frag2alitbl(a1,a2,acnt+1,len,ali1,ali2);

 sco=AliTbl2Sco(ali1,ali2,(acnt+1)*len,d1,d2,&rate);
 //sco=AliTbl2Sco_sq(ali1,ali2,(acnt+1)*len,d1,d2,&rate);

 if(rate <rcut)
  break;

 //no improvement
 if(sco <= pre)
  break;
  pre=sco;
  acnt++;
  backflag++;
  //printf("----ID %d:%d -> %d:%d sco= %.3f ALEN= %d CAD= %.3f RATE= %.3f\n",id11,id12,maxid[0],maxid[1],sco,acnt,sco/d1->Tarea,rate);

 }
 //input alignment
/*
 if((a->a1[a->N]=(int *)malloc(sizeof(int)*d1->Max))==NULL)
  return 0;
 if((a->a2[a->N]=(int *)malloc(sizeof(int)*d2->Max))==NULL)
  return 0;
*/
 if((a->ali[a->N].a1=(int *)malloc(sizeof(int)*d1->Max))==NULL)
  return 0;
 if((a->ali[a->N].a2=(int *)malloc(sizeof(int)*d2->Max))==NULL)
  return 0;

 //frag2ali(a1,a2,acnt,len,a->a1[a->N],a->a2[a->N],n1,n2);
 frag2ali(a1,a2,acnt,len,a->ali[a->N].a1,a->ali[a->N].a2,n1,n2);
 //a->sco[a->N]=pre;
 a->ali[a->N].sco=pre;
 a->N++;
 return 1;
}

int GeneTripletsAndAli(ALI *a,int id11,int id12,int *tbl,int n,int len,DMTX *d1,DMTX *d2){
 int n1=d1->Max;
 int n2=d2->Max;
 int i1,i2;
 int id21,id22,id31,id32;
 int cnt=0;
 double sco;
 for(i1=0;i1<n;i1++){
  id21=tbl[2*i1];
  id22=tbl[2*i1+1];
  //ignore overlap
  if(abs(id11-id21)<len)
   continue;
  if(abs(id12-id22)<len)
   continue;
 	for(i2=i1+1;i2<n;i2++){
	 id31=tbl[2*i2];
	 id32=tbl[2*i2+1];
	 //remove overlap
	 if(abs(id11-id31)<len||abs(id21-id31)< len)
  	  continue;
	 if(abs(id12-id32)<len||abs(id22-id32)< len)
  	  continue;
	 //printf("%d) %d-%d %d-%d %d-%d\n",*Cnt,id11,id12,id21,id22,id31,id32);
	 //Tri2Ali(a->a1[a->N],a->a2[a->N],id11,id12,id21,id22,id31,id32,n1,n2,len);
	 //sco=Ali2Sco(a1,a2,d1,d2);
	 //printf("Cnt= %d %d\n",*Cnt,a1[(*Cnt)][0]);
	 cnt++;
	 a->N++;
	}
 }

 return cnt;
}

int Tri2Ali(int *a1,int *a2,int i11,int i12,int i21,int i22,int i31,int i32,int len1,int len2,int seg_len){
 int i;
 //init
 for(i=0;i<len1;i++)
  a1[i]=-1;
 for(i=0;i<len2;i++)
  a2[i]=-1;
 //input
 for(i=0;i<seg_len;i++){
  //seg1
  a1[i11+i]=i12+i;
  a2[i12+i]=i11+i;
  //seg2
  a1[i21+i]=i22+i;
  a2[i22+i]=i21+i;
  //seg3
  a1[i31+i]=i32+i;
  a2[i32+i]=i31+i;
  //printf("%d - %d %d\n",i11+i,i12+i,a1[i11+i]);
 }
 return 0;
}

int FilterSeeds(ALI *a,DMTX *d1,DMTX *d2){
 int i,j;
 int n=a->N;
 double rate;
 for(i=0;i<n;i++){
  //printf("#Seeds %d\n",i);
  //scoring....
  //a->sco[i]=Ali2Sco(a->a1[i],d1,d2,&rate);
  //printf("%d %f\n",i,a->sco[i]);
 }
 a->ave=Ave(a->sco,n);
 a->std=Std(a->sco,n);
 //printf("Ave= %f Std= %f\n",ave,std);
}


//Scoring function*********
double Ali2Sco(int *a1,DMTX *d1,DMTX *d2,double *rate){
 int i,j,p1,p2; 
 double Tarea,TareaAlign,bd,cad;
 //double miss=0;
 //double wc=0;
 double A1,A2;
 Tarea=TareaAlign=bd=0;
 for(i=0;i<d1->Max;i++){
  p1=a1[i];
  if(p1<0)
   continue;
 	for(j=i+1;j<d1->Max;j++){
	 p2=a1[j];
	 //no aligned res
	 if(p2 < 0)
	  continue;
	 
	 A1=d1->area[i][j]+d1->area[j][i];
	 A2=d2->area[p1][p2]+d2->area[p2][p1];

	 bd+=cad_bnd(A1,A2);
	 bd+=cad_bnd(A2,A1);
	 Tarea+=A1+A2;
	 TareaAlign+=A1;
/* asymmetry
	 //i -> j
	 bd+=cad_bnd(d1->area[i][j],d2->area[p1][p2]); 
	 //j -> i	 
	 bd+=cad_bnd(d1->area[j][i],d2->area[p2][p1]);
	 Tarea+=d1->area[i][j]+d1->area[j][i];
*/
	 //Symmetry
	 

	}
 }
 *rate=(Tarea-bd)/(TareaAlign*2.00);
 //printf("Sco= %.3f Rate= %.3f = %.1f/%.1f ali= %.1f\n",Tarea-bd,*rate,Tarea-bd,TareaAlign*2,TareaAlign);
 return(Tarea-bd);
 //return(Tarea-bd-wc);
}

//ignore neighber residues +-4:helix
double Ali2ScoSep(int *a1,DMTX *d1,DMTX *d2,double *rate,int L){
 int i,j,p1,p2; 
 double Tarea,TareaAlign,bd,cad;
 //double miss=0;
 //double wc=0;
 double A1,A2;
 int chk_p2;
 Tarea=TareaAlign=bd=0;
 for(i=0;i<d1->Max;i++){
  p1=a1[i];
  if(p1<0)
   continue;
  chk_p2=p1+L+1;
 	for(j=i+L+1;j<d1->Max;j++){
	 p2=a1[j];
	 //no aligned res
	 if(p2 < 0)
	  continue;
	 if(p2<chk_p2)
	  continue;

	 
	 A1=d1->area[i][j]+d1->area[j][i];
	 A2=d2->area[p1][p2]+d2->area[p2][p1];

	 bd+=cad_bnd(A1,A2);
	 bd+=cad_bnd(A2,A1);
	 Tarea+=A1+A2;
	 TareaAlign+=A1;
	}
 }
 *rate=(Tarea-bd)/(TareaAlign*2.00);
 if(TareaAlign==0)
  *rate=0;
 return(Tarea-bd);
}

//Weight
double Ali2ScoSepW(int *a1,DMTX *d1,DMTX *d2,double *rate,double w){
 int i,j,p1,p2; 
 double Tarea,TareaAlign,bd,cad;
 //double miss=0;
 //double wc=0;
 double A1,A2;
 int chk_p2;
 int near1,near2;
 Tarea=TareaAlign=bd=0;
 for(i=0;i<d1->Max;i++){
  if(a1[i]<0)
   continue;
  p1=a1[i];

 	for(j=i+1;j<d1->Max;j++){
	 //no aligned res
	 if(a1[j] < 0)
	  continue;
	 p2=a1[j];

	 A1=d1->area[i][j]+d1->area[j][i];
	 A2=d2->area[p1][p2]+d2->area[p2][p1];

	 //near
	 if(j-i < SEP_LEN || (p1-p2 < SEP_LEN && p2-p1 < SEP_LEN)){
	  A1*=w;
	  A2*=w;
	 }


	 bd+=cad_bnd(A1,A2);
	 bd+=cad_bnd(A2,A1);
	 Tarea+=A1+A2;
	 TareaAlign+=A1;
	}
 }
 *rate=(Tarea-bd)/(TareaAlign*2.00);
 if(TareaAlign==0)
  *rate=0;
 return(Tarea-bd);
}



//Scoring function*********
double Ali2Sco_sq(int *a1,DMTX *d1,DMTX *d2,double *rate){
 int i,j,p1,p2; 
 double Tarea,bd,cad,Total,sco;
 double miss=0;
 Tarea=bd=Total=sco=0;
 //Tarea=100;
 for(i=0;i<d1->Max;i++){
 	for(j=i+1;j<d1->Max;j++){
	 if(d1->area[i][j]==0||d1->area[j][i]==0)
	  continue;

	 p1=a1[i];
	 p2=a1[j];
	 //no aligned res
	 if(p1<0 || p2 < 0){
	  //miss+=d1->area[i][j],d1->area[j][i];
	  continue;
	 }
	 //i -> j
	 bd=cad_bnd(d1->area[i][j],d2->area[p1][p2]); 
	 //j -> i	 
	 bd+=cad_bnd(d1->area[j][i],d2->area[p2][p1]);
	 Tarea=d1->area[i][j]+d1->area[j][i];
	 sco+=sqrt(Tarea-bd);
	 Total+=sqrt(Tarea);
	}
 }
 *rate=sco/Total;
 return(sco);
}

double AliTbl2Sco(int *a1,int *a2,int n,DMTX *d1,DMTX *d2,double *rate){
 int i,j,p11,p12,p21,p22; 
 double Tarea,bd,cad,wc;
 Tarea=bd=wc=0;
 //Tarea=100;
 for(i=0;i<n;i++){
  p11=a1[i]; p21=a2[i];
 	for(j=i+1;j<n;j++){
  	 p12=a1[j]; p22=a2[j];
	 //if(d1->area[p11][p12]==0||d1->area[p12][p11]==0)
	 // continue;

	 //i -> j
	 bd+=cad_bnd(d1->area[p11][p12],d2->area[p21][p22]); 
	 //j -> i	 
	 bd+=cad_bnd(d1->area[p12][p11],d2->area[p22][p21]); 
	 //bd+=cad_bnd(d1->area[j][i],d2->area[p2][p1]);
	 //Tarea+=d1->area[i][j]+d1->area[j][i];
	 Tarea+=d1->area[p11][p12]+d1->area[p12][p11];
 //printf("%d:%d Bd= %.3f Tarea= %.3f/%.3f Sco= %.3f\n",p11,p12,bd,Tarea,d1->Tarea,1-bd/Tarea);
	 //if(d1->area[p11][p12]==0||d1->area[p12][p11]==0)
	 // wc+=d2->area[p21][p22]+d2->area[p22][p21];
	}
 }
 //printf("Bd= %.3f Tarea= %.3f/%.3f Sco= %.3f\n",bd,Tarea,d1->Tarea,1-bd/Tarea);
 *rate=1-bd/Tarea;
 return(Tarea-bd);
 //return(Tarea-bd-wc);
}

double AliTbl2Sco_sq(int *a1,int *a2,int n,DMTX *d1,DMTX *d2,double *rate){
 int i,j,p11,p12,p21,p22; 
 double Tarea,bd,cad,sco,Total;
 Tarea=bd=sco=Total=0;
 //Tarea=100;
 for(i=0;i<n;i++){
  p11=a1[i]; p21=a2[i];
 	for(j=i+1;j<n;j++){
  	 p12=a1[j]; p22=a2[j];

	 //i -> j
	 bd=cad_bnd(d1->area[p11][p12],d2->area[p21][p22]); 
	 //j -> i	 
	 bd+=cad_bnd(d1->area[p12][p11],d2->area[p22][p21]); 
	 Tarea=d1->area[p11][p12]+d1->area[p12][p11];
	 sco+=sqrt(Tarea-bd);
	 Total+=sqrt(Tarea);
	}
 }

 *rate=sco/Total;
 return(sco);
}

//Cal only delta score
double BetweenFrag(int a11,int a12,int a21,int a22,int len,DMTX *d1,DMTX *d2,double *rate){
 int i,j,p11,p12,p21,p22; 
 double Tarea,bd,cad,sco,Total;
 Tarea=bd=sco=Total=0;
 //Tarea=100;
 for(i=0;i<len;i++){
  p11=a11+i; p21=a21+i;
 	for(j=0;j<len;j++){
  	 p12=a12+j; p22=a22+j;

	 //|-a11-| <=> |-a12-|
	 //|-a21-| <=> |-a22-|
	 //i -> j
	 //bd=cad_bnd(d1->area[p11][p12],d2->area[p21][p22]); 
	 bd+=cad_bnd(d1->area[p11][p12],d2->area[p21][p22]); 
	 //j -> i	 
	 bd+=cad_bnd(d1->area[p12][p11],d2->area[p22][p21]); 
	 //Tarea=d1->area[p11][p12]+d1->area[p12][p11];
	 Tarea+=d1->area[p11][p12]+d1->area[p12][p11];

	 //sco+=sqrt(Tarea-bd);
	 //Total+=sqrt(Tarea);
	}
 }
 if(Tarea==0){//no contact
  *rate=0;
  return 0;
 }
  
 //*rate=sco/Total;
 *rate=1-bd/Tarea;//Simple ver
 //return(sco);
 return(Tarea-bd);
}



//Basic Movement
int CopyAli(int *a1,int *a2,int *b1,int *b2,int n1,int n2){
 int i;
 for(i=0;i<n1;i++)
  b1[i]=a1[i];
 for(i=0;i<n2;i++)
  b2[i]=a2[i];
 return 0;
}
int CopyAliSingle(int *a1,int *b1,int n1){
 int i;
 for(i=0;i<n1;i++)
  b1[i]=a1[i];
 return 0;
}

//Move and scoring
double ExtAli(int *a1,int *a2,DMTX *d1,DMTX *d2,double sco,double rcut){
 int i,j;
 int new;
 int cnt=0;
 int list1[1000],list2[1000];
 //search candidate position
 for(i=0;i<d1->Max;i++){
  if(a1[i]!=-1)
   continue;
  if(cnt >1000)
   continue;
  //case1   AAAAAA*---
  if(i>0 && a1[i-1] != -1 && a2[a1[i-1]+1]==-1){
   list1[cnt]=i;
   list2[cnt]=a1[i-1]+1;
   cnt++;
  }
  //case2   ----*AAAAA
  if(i<d1->Max-1 && a1[i+1] !=-1 && a2[a1[i+1]-1]==-1){
   list1[cnt]=i;
   list2[cnt]=a1[i+1]-1;
   cnt++;
  }
 }
 //no candidate position
 if(cnt==0)
  return 0;

 //check only delta
 double delta=0;
 double bnd,cad,t=0;

 //new=range_rand_i(0,(double)cnt);
 //search best extention

 //int p1=list1[new];
 //int p2=list2[new];
 int p1,p2,maxid=-1;
 double max=0,maxt=0;
 for(i=0;i<cnt;i++){
  bnd=0;
  p1=list1[i];
  p2=list2[i];
  t=0;
 	for(j=0;j<d1->Max;j++){
 	 if(p1==j||a1[j]==-1)
 	  continue;
 	 bnd+=cad_bnd(d1->area[p1][j],d2->area[p2][a1[j]]);
 	 bnd+=cad_bnd(d1->area[j][p1],d2->area[a1[j]][p2]);
 	 t+=(d1->area[j][p1]+d1->area[p1][j]);
	}
	if(1-bnd/t < rcut)
	 continue;
	if(max<t-bnd){
	 max=t-bnd;
	 maxid=i;
	 maxt=t;
	}
 }
 if(maxid==-1)
  return 0;
 a1[list1[maxid]]=list2[maxid];
 a2[list2[maxid]]=list1[maxid];
 //printf("bnd= %f\n",bnd);
 return(max);
}

double TrimAli(int *a1,int *a2,DMTX *d1,DMTX *d2){

}

int ShiftAli(int *a1,int *a2,int n1,int n2,int p,int tml){
 int i,j,k;//p1-----tml

 for(i=p;i<tml;i++){//p1-----tml
  a1[i]=a1[i+1];
  a1[i+1]=-1;
  if(a1[i]>-1)
   a2[a1[i]]=i;
 }
 for(i=p;i>tml;i--){//tml---p1
  a1[i]=a1[i-1];
  a1[i-1]=-1;
  if(a1[i]>-1)
   a2[a1[i]]=i;
 }

 return 0;
}

double ShiftAliSco(int *a1,int *a2,DMTX *d1,DMTX *d2,double sco,double rcut,int len){
 int i,j;
 int new;
 int cntR=0;
 int cntL=0;
 int listR[1000];//shift Right
 int listL[1000];//shift Left
 int tmpa1[1000],tmpa2[1000];
 int besta1[1000],besta2[1000];
 double tmpsco,rate;

 //copy
 CopyAli(a1,a2,tmpa1,tmpa2,d1->Max,d2->Max);

 //search candidate position
 for(i=0;i<d1->Max;i++){
  if(a1[i]!=-1)
   continue;
  if(cntR >1000||cntL>1000)
   continue;
  //case1   AAAAAA*--- gap was shift to left
  if(i>0 && a1[i-1] != -1){
   listL[cntL]=i;
   cntL++;
  }
  //case2   ----*AAAAA to right
  if(i<d1->Max-1 && a1[i+1] !=-1){
   listR[cntR]=i;
   cntR++;
  }
 }

 int p1,p2,maxid=-1;
 double max=sco,maxt=0,bnd,t;
 int s;
 //shift Left
 for(i=0;i<cntL;i++){
  bnd=0;
  p1=listL[i];
  t=0;
  //printf("ShiftL %d %d\n",i,p1);
  //for(p2=0;p2<d1->Max;p2++)
  // printf("%4d",a1[p2]);
  //puts("");
	for(s=1;s<=len && p1-s >0;s++){
	 //printf("%d - %d\n",p1,p1-s);
 	 CopyAli(a1,a2,tmpa1,tmpa2,d1->Max,d2->Max);
	 ShiftAli(tmpa1,tmpa2,d1->Max,d2->Max,p1,p1-s);//gap shift left
	 tmpsco=Ali2Sco(tmpa1,d1,d2,&rate);
	 if(max<tmpsco){
	  //printf("shiftL %d sco= %.3f\n",s,tmpsco);
	  max=tmpsco;
 	  CopyAli(tmpa1,tmpa2,besta1,besta2,d1->Max,d2->Max);
	  
	 }
		//for(p2=0;p2<d2->Max;p2++)
		// printf("%4d",tmpa2[p2]);
		//puts("");
	}
 }
 //shift Right
 for(i=0;i<cntR;i++){
  bnd=0;
  p1=listR[i];
  t=0;
  //printf("ShiftR %d %d\n",i,p1);
	for(s=1;s<=len && p1+s <d1->Max;s++){
	 //printf("%d - %d\n",p1,p1+s);
 	 CopyAli(a1,a2,tmpa1,tmpa2,d1->Max,d2->Max);
	 ShiftAli(tmpa1,tmpa2,d1->Max,d2->Max,p1,p1+s);//gap shift left
	 tmpsco=Ali2Sco(tmpa1,d1,d2,&rate);
		//for(p2=0;p2<d2->Max;p2++)
		// printf("%4d",tmpa2[p2]);
		//puts("");
	 if(max<tmpsco){
	  max=tmpsco;
 	  CopyAli(tmpa1,tmpa2,besta1,besta2,d1->Max,d2->Max);
	  //printf("shiftR %d sco= %.3f\n",s,tmpsco);
	 }
	}
 }
 if(max-sco==0)
  return 0;
 //renew a1 and a2
 CopyAli(besta1,besta2,a1,a2,d1->Max,d2->Max);
 return(max-sco);
}

//trim & extension
int RefineAli(ALI *a,DMTX *d1,DMTX *d2,double z,double rcut,int len,PDB *p1,PDB *p2){
 int i,j;
 int n=a->N;
 double df1,df2,df3;
 double asco,rate,orisco;
 double pre;
 int tmp[1000];

 for(i=0;i<n && i < 5;i++){
  pre=a->ali[i].sco;
  //asco=Ali2Sco_sq(a->ali[i].a1,d1,d2,&rate);
  printf("No%d sco= %.3f pre= %.3f rate= %.3f CAD= %.3f\n",i+1,a->ali[i].sco,pre,rate,a->ali[i].sco/d1->Tarea);
  //Filtering
 
 }

}

int cmp_tria(const void *c1, const void *c2){
 TRIA a1=*(TRIA *)c1;
 TRIA a2=*(TRIA *)c2;
 if(a1.area<a2.area)//Big -> small
  return 1;
 else
  return 0;
}


int Tri(DMTX *d1,TRIA **t,int AllowOv){
 int i1,i2,i3;
 int cnt1=0;
 int cnt2=0;
 double acut=10;
 double a;
 for(i1=0;i1<d1->Max;i1++){
  for(i2=i1+1;i2<d1->Max;i2++){
   if(d1->area[i1][i2]==0||d1->area[i2][i1]==0)
    continue;
   for(i3=i2+1;i3<d1->Max;i3++){
    if(d1->area[i1][i3]==0||d1->area[i3][i1]==0)
     continue;
    if(d1->area[i2][i3]==0||d1->area[i3][i2]==0)
     continue;
    a=d1->area[i1][i2]+d1->area[i1][i3]+d1->area[i2][i3];
    a+=d1->area[i2][i1]+d1->area[i3][i1]+d1->area[i3][i2];
    if(a<acut)
     continue;
    //printf("1:%d) %d %d %d %.3f\n",cnt1,i1,i2,i3,a);
    cnt1++;
   }
  }
 }
 //int *tlist1,*tlist2;
 //tlist1=(int *)malloc(sizeof(int)*3*cnt1);
 TRIA *t1,*t2;
 if((t1=(TRIA *)malloc(sizeof(TRIA)*cnt1))==NULL)
  return 0;
/*
 if((t1=(TRIA **)malloc(sizeof(TRIA*)*cnt1))==NULL)
  return 0;
 for(i1=0;i1<cnt1;i1++)
  if((t1[i1]=(TRIA *)malloc(sizeof(TRIA)))==NULL)
   return 0;
*/
 cnt1=0;
 for(i1=0;i1<d1->Max;i1++){
  for(i2=i1+1;i2<d1->Max;i2++){
   if(d1->area[i1][i2]==0||d1->area[i2][i1]==0)
    continue;
   for(i3=i2+1;i3<d1->Max;i3++){
    if(d1->area[i1][i3]==0||d1->area[i3][i1]==0)
     continue;
    if(d1->area[i2][i3]==0||d1->area[i3][i2]==0)
     continue;
    a=d1->area[i1][i2]+d1->area[i1][i3]+d1->area[i2][i3];
    a+=d1->area[i2][i1]+d1->area[i3][i1]+d1->area[i3][i2];
    if(a<acut)
     continue;
    //tlist1[cnt1*3]=i1; tlist1[cnt1*3+1]=i2; tlist1[cnt1*3+2]=i3;
    t1[cnt1].r[0]=i1;
    t1[cnt1].r[1]=i2;
    t1[cnt1].r[2]=i3;
    t1[cnt1].area=a;
    cnt1++;
   }
  }
 }
/*
 for(i1=0;i1<d2->Max;i1++){
  for(i2=i1+1;i2<d2->Max;i2++){
   if(d2->area[i1][i2]==0)
    continue;
   for(i3=i2+1;i3<d2->Max;i3++){
    if(d2->area[i1][i3]==0)
     continue;
    if(d2->area[i2][i3]==0)
     continue;
    a=d2->area[i1][i2]+d2->area[i1][i3]+d2->area[i2][i3];
    a+=d2->area[i2][i1]+d2->area[i3][i1]+d2->area[i3][i2];
    if(a<acut)
     continue;
    //printf("2:%d) %d %d %d %.3f\n",cnt2,i1,i2,i3,a);
    cnt2++;
   }
  }
 }
 //tlist2=(int *)malloc(sizeof(int)*3*cnt2);
 if((t2=(TRIA *)malloc(sizeof(TRIA)*cnt2))==NULL)
  return 0;
 cnt2=0;
 for(i1=0;i1<d2->Max;i1++){
  for(i2=i1+1;i2<d2->Max;i2++){
   if(d2->area[i1][i2]==0)
    continue;
   for(i3=i2+1;i3<d2->Max;i3++){
    if(d2->area[i1][i3]==0)
     continue;
    if(d2->area[i2][i3]==0)
     continue;
    a=d2->area[i1][i2]+d2->area[i1][i3]+d2->area[i2][i3];
    a+=d2->area[i2][i1]+d2->area[i3][i1]+d2->area[i3][i2];
    if(a<acut)
     continue;
    //tlist2[cnt2*3]=i1; tlist2[cnt2*3+1]=i2; tlist2[cnt2*3+2]=i3;
    t2[cnt2].r[0]=i1;
    t2[cnt2].r[1]=i2;
    t2[cnt2].r[2]=i3;
    t2[cnt2].area=a;
    cnt2++;
   }
  }
 }
*/
 //printf("Ntri= %d * %d *3 = %d\n",cnt1,cnt2,cnt1*cnt2*3);

 //sort and redundant
 qsort(t1,cnt1,sizeof(TRIA),cmp_tria);
 //qsort(t2,cnt2,sizeof(TRIA),cmp_tria);

 int *tag,NumNr=0;
 if((tag=(int *)malloc((d1->Max)*sizeof(int)))==NULL)
  return 0;
 for(i1=0;i1<d1->Max;i1++)
  tag[i1]=-1;

 int flag=0;
 for(i1=0;i1<cnt1;i1++){
  //check
  flag=0;
  if(tag[t1[i1].r[0]]!=-1) flag++;
  if(tag[t1[i1].r[1]]!=-1) flag++;
  if(tag[t1[i1].r[2]]!=-1) flag++;
  if(flag>AllowOv)
   continue;
  tag[t1[i1].r[0]]=tag[t1[i1].r[1]]=tag[t1[i1].r[2]]=1;
  //Shift
  t1[NumNr].r[0]=t1[i1].r[0];
  t1[NumNr].r[1]=t1[i1].r[1];
  t1[NumNr].r[2]=t1[i1].r[2];
  t1[NumNr].area=t1[i1].area;
  NumNr++;
 }

 *t=t1;
 free(tag);

 return NumNr;
}

#define my_abs(x) ((x) >= 0 ? (x) : -(x))

double cad_bnd(double a1,double a2){
 //double cad=fabs(a1-a2);
 double cad=a1-a2;
 if(cad<0)
  cad*=-1.00;
 //double cad=my_abs(a1-a2);

// if(a1>a2)
//  return cad;
// else
//  return 0;

 if(cad<a1)
  return cad;
 return a1;


/*
 if(a1>a2)
  return cad;
 if(cad*0.5<a1)
  return 0.5*cad;
 return 0;
*/

}

double cad_bnd_sym(double a1,double a2){
 //double cad=fabs(a1-a2);
 double cad=a1-a2;
 double sco=0;
 if(cad<0)
  cad*=-1.00;

 if(cad<a1)
  sco=a1-cad;
 if(cad<a2)
  sco+=a2-cad;

 return sco;
}


int check_near(int *i1,int *i2, DMTX *d){
 int i,j;
 for(i=0;i<3;i++)
  for(j=0;j<3;j++)
   if(d->area[i1[i]][i2[j]]>0||d->area[i1[j]][i2[i]]>0||i1[i]==i2[j])
    return 1;

 return 0;
}
int TriPair(TRIA *t1,int n1,TRIA *t2,int n2,DMTX *d1,DMTX *d2){
 int i,j,k,l,m,s;
 int cnt1=0;
 int cnt2=0;
 for(i=0;i<n1;i++){
  for(j=i+1;j<n1;j++){
   //near?
   if(check_near(t1[i].r,t1[j].r,d1)==0)
    continue;
   //printf("Pair %d %d %d %d\n",i,j,t1[i].r[0],t1[j].r[0]);
   cnt1++;
  }
 }
 printf("cnt1= %d cnt2= %d\n",cnt1,cnt2);
}


int UpperLowerBound(TRIA *t1,int n1,TRIA *t2,int n2,DMTX *d1,DMTX *d2){
 int i,j,k,l,m,s;
 int p1,p2,p3,q1,q2,q3;
 double cad,bnd,maxbnd,minbnd;
 double o_max,o_min,sum_max,sum_min,sum_a;
 for(i=0;i<n1;i++){

  p1=t1[i].r[0]; p2=t1[i].r[1]; p3=t1[i].r[2];
  minbnd=100000;
  maxbnd=0;

	for(j=0;j<n2;j++){

	 //six pattern
	 //pattern 1 0:0, 1:1, 2:2
  	 q1=t2[j].r[0]; q2=t2[j].r[1]; q3=t2[j].r[2];
	 bnd =cad_bnd(d1->area[p1][p2],d2->area[q1][q2]);
	 bnd+=cad_bnd(d1->area[p2][p1],d2->area[q2][q1]);
	 bnd+=cad_bnd(d1->area[p1][p3],d2->area[q1][q3]);
	 bnd+=cad_bnd(d1->area[p3][p1],d2->area[q3][q1]);
	 bnd+=cad_bnd(d1->area[p3][p2],d2->area[q3][q2]);
	 bnd+=cad_bnd(d1->area[p2][p3],d2->area[q2][q3]);

	 if(bnd > maxbnd)
	  maxbnd=bnd;
	 if(bnd < minbnd)
	  minbnd=bnd;

	 //pattern 2 0:1, 1:2, 2:0
	 q1=t2[j].r[1]; q2=t2[j].r[2]; q3=t2[j].r[0];
	 bnd =cad_bnd(d1->area[p1][p2],d2->area[q1][q2]);
	 bnd+=cad_bnd(d1->area[p2][p1],d2->area[q2][q1]);
	 bnd+=cad_bnd(d1->area[p1][p3],d2->area[q1][q3]);
	 bnd+=cad_bnd(d1->area[p3][p1],d2->area[q3][q1]);
	 bnd+=cad_bnd(d1->area[p3][p2],d2->area[q3][q2]);
	 bnd+=cad_bnd(d1->area[p2][p3],d2->area[q2][q3]);

	 if(bnd > maxbnd)
	  maxbnd=bnd;
	 if(bnd < minbnd)
	  minbnd=bnd;

	 //pattern 3 0:2, 1:0, 2:1
	 q1=t2[j].r[2]; q2=t2[j].r[0]; q3=t2[j].r[1];
	 bnd =cad_bnd(d1->area[p1][p2],d2->area[q1][q2]);
	 bnd+=cad_bnd(d1->area[p2][p1],d2->area[q2][q1]);
	 bnd+=cad_bnd(d1->area[p1][p3],d2->area[q1][q3]);
	 bnd+=cad_bnd(d1->area[p3][p1],d2->area[q3][q1]);
	 bnd+=cad_bnd(d1->area[p3][p2],d2->area[q3][q2]);
	 bnd+=cad_bnd(d1->area[p2][p3],d2->area[q2][q3]);

	 if(bnd > maxbnd)
	  maxbnd=bnd;
	 if(bnd < minbnd)
	  minbnd=bnd;
	
	 //pattern 4 0:0, 1:2, 2:1
  	 q1=t2[j].r[0]; q2=t2[j].r[2]; q3=t2[j].r[1];
	 bnd =cad_bnd(d1->area[p1][p2],d2->area[q1][q2]);
	 bnd+=cad_bnd(d1->area[p2][p1],d2->area[q2][q1]);
	 bnd+=cad_bnd(d1->area[p1][p3],d2->area[q1][q3]);
	 bnd+=cad_bnd(d1->area[p3][p1],d2->area[q3][q1]);
	 bnd+=cad_bnd(d1->area[p3][p2],d2->area[q3][q2]);
	 bnd+=cad_bnd(d1->area[p2][p3],d2->area[q2][q3]);

	 if(bnd > maxbnd)
	  maxbnd=bnd;
	 if(bnd < minbnd)
	  minbnd=bnd;

	 //pattern 5 0:1, 1:0, 2:2
	 q1=t2[j].r[1]; q2=t2[j].r[0]; q3=t2[j].r[2];
	 bnd =cad_bnd(d1->area[p1][p2],d2->area[q1][q2]);
	 bnd+=cad_bnd(d1->area[p2][p1],d2->area[q2][q1]);
	 bnd+=cad_bnd(d1->area[p1][p3],d2->area[q1][q3]);
	 bnd+=cad_bnd(d1->area[p3][p1],d2->area[q3][q1]);
	 bnd+=cad_bnd(d1->area[p3][p2],d2->area[q3][q2]);
	 bnd+=cad_bnd(d1->area[p2][p3],d2->area[q2][q3]);

	 if(bnd > maxbnd)
	  maxbnd=bnd;
	 if(bnd < minbnd)
	  minbnd=bnd;

	 //pattern 6 0:2, 1:1, 2:0
	 q1=t2[j].r[2]; q2=t2[j].r[1]; q3=t2[j].r[0];
	 bnd =cad_bnd(d1->area[p1][p2],d2->area[q1][q2]);
	 bnd+=cad_bnd(d1->area[p2][p1],d2->area[q2][q1]);
	 bnd+=cad_bnd(d1->area[p1][p3],d2->area[q1][q3]);
	 bnd+=cad_bnd(d1->area[p3][p1],d2->area[q3][q1]);
	 bnd+=cad_bnd(d1->area[p3][p2],d2->area[q3][q2]);
	 bnd+=cad_bnd(d1->area[p2][p3],d2->area[q2][q3]);

	 if(bnd > maxbnd)
	  maxbnd=bnd;
	 if(bnd < minbnd)
	  minbnd=bnd;

	}
	//other contact
	sum_min=sum_max=sum_a=0;

	for(s=0;s<3;s++){
	 p1=t1[i].r[s];
	for(j=0;j<d1->Max;j++){
	 if(d1->area[p1][j]==0||d1->area[j][p1]==0)
	  continue;
	 if(j==t1[i].r[0]||j==t1[i].r[1]||j==t1[i].r[2])
	  continue;
	 o_max=0;o_min=100000;
		for(k=0;k<d2->Max;k++){
		 for(l=k+1;l<d2->Max;l++){
	 	  if(d2->area[k][l]==0||d2->area[l][k]==0)
		   continue;

		  bnd =cad_bnd(d1->area[p1][j],d2->area[k][l]);
		  bnd+=cad_bnd(d1->area[j][p1],d2->area[l][k]);
		  if(o_max < bnd) o_max=bnd;
		  if(o_min > bnd) o_min=bnd;
		 }
		}
	 sum_min+=o_min;
	 sum_max+=o_max;
	 sum_a+=d1->area[p1][j]+d1->area[j][p1];
	}
  	}
  //printf("maxbnd= %.3f minbnd= %.3f /%.3f\n",maxbnd,minbnd,t1[i].area);
  t1[i].upb=(sum_a+t1[i].area)-(minbnd+sum_min);
  printf("%d sum_max= %.3f sum_min= %.3f /%.3f UPB= %.3f\n",i,sum_max,sum_min,sum_a,t1[i].upb);
 }
}
int TreeAlign(TRIA *t1,int n1,TRIA *t2,int n2,DMTX *d1,DMTX *d2,ALI *a){
 int i,pos;
 int tali[1000];

 //Random Alignment

 //init
 for(i=0;i<n1;i++){
  tali[i]=-1;//no alignment
  //printf("%.3f\n",t1[i].upb);
 }
 pos=n1-1;//position
 //return 0;
 while(1){
  tali[pos]++;

  for(i=0;i<=pos;i++)
   printf("%4d",tali[i]);
  puts("");

  if(tali[pos]==n2){
   tali[pos]=-1;
   pos--;
  }else if(pos<n1-1){
   pos++;
  }

 }

}

int NrAlignments(ALI *a){
 int i,j,k;
 bool *tag,flag1=false;
 int *sum;
 bool flag2;
 if((tag=(bool*)malloc(sizeof(bool)*a->N))==NULL)
  return FALSE;
 if((sum=(int*)calloc(a->N,sizeof(int)))==NULL)
  return FALSE;

 //init
 for(i=0;i<a->N;i++){
  tag[i]=true;
  for(j=0;j<a->len1;j++){
   //sum[i]+=a->a1[i][j];
   sum[i]+=a->ali[i].a1[j];
  }
 }

 for(i=0;i<a->N;i++){
  //flag1=false;
  for(j=i+1;j<a->N;j++){
   if(sum[i]!=sum[j])
    continue;
   if(tag[j]==false)
    continue;
   flag2=false;
	for(k=0;k<a->len1;k++){
	 //if(a->a1[i][k]!=a->a1[j][k]){
	 if(a->ali[i].a1[k]!=a->ali[j].a1[k]){
	  flag2=true;
	  break;
	 }
	}
   if(flag2==false){
    //printf("%d %d same\n",i,j);
    tag[j]=false;
   }
  }
 }
 //shift
 int n=0;
 for(i=0;i<a->N;i++){
  if(tag[i]==false)
   continue;
  //printf("%d <- %d\n",n,i);
  if(n!=i){
/*
   a->a1[n]=a->a1[i];
   a->a2[n]=a->a2[i];
   a->sco[n]=a->sco[i];
   a->ali[n]=a->ali[i];//new
*/
   a->ali[n]=a->ali[i];
   //printf("NR %d <- %d %f\n",n,i,a->ali[i].sco);
  }
  n++;
 }

 free(tag);
 free(sum);
 a->N=n;
 return n;
}

int NrAlignmentsRate(ALI *a,double rate){
 int i,j,k;
 bool *tag,flag1=false;
 bool flag2;
 if((tag=(bool*)malloc(sizeof(bool)*a->N))==NULL)
  return FALSE;
 int cnt;
 int check=(int)((1.00-rate)*a->len1);

 //init
 for(i=0;i<a->N;i++){
  tag[i]=true;//keep
 }

 for(i=0;i<a->N;i++){
  for(j=i+1;j<a->N;j++){
   if(tag[j]==false)
    continue;

	cnt=0;
	for(k=0;k<a->len1;k++){
	 if(a->ali[i].a1[k]!=a->ali[j].a1[k]){
	  cnt++;
	 }
	 if(cnt > check)
	  break;
	}
    //printf("Same %d %d\n",cnt,check);
   if(cnt <= check){
    tag[j]=false;//remove
   }
  }
 }
 //shift
 int n=0;
 for(i=0;i<a->N;i++){
  if(tag[i]==false)
   continue;
  //printf("%d <- %d\n",n,i);
  if(n!=i){
/*
   a->a1[n]=a->a1[i];
   a->a2[n]=a->a2[i];
   a->sco[n]=a->sco[i];
   a->ali[n]=a->ali[i];//new
*/
   a->ali[n]=a->ali[i];
   //printf("NR %d <- %d %f\n",n,i,a->ali[i].sco);
  }
  n++;
 }

 free(tag);
 a->N=n;
 return n;
}





//from dali elastic similarity score
#define THETA 0.20
#define ALPHA 20.0
int LocalStrSim(PDB *p1,PDB *p2,int len,double *out){
 int i,n1,n2;
 int j,k,i1,i2;
 double d1,d2,tmp,tmp2;
 //for(i=0;i<n1-len;i++){
 for(i=0;i<p1->NumOfRes-len;i++){
 	//for(j=0;j<n2-len;j++){
 	for(j=0;j<p2->NumOfRes-len;j++){
 	 tmp=0;
	 //printf("%d vs %d\n",i,j);
	 for(i1=0;i1<len;i1++){
	  for(i2=0;i2<len;i2++){
	   if(i1==i2){
	    tmp+=THETA;
	    continue;
	   }

	   //d1=L(cd1[i+i1][0],cd1[i+i1][1],cd1[i+i1][2],cd1[i+i2][0],cd1[i+i2][1],cd1[i+i2][2]);
	   d1=L(p1->CAxyz[i+i1][0],p1->CAxyz[i+i1][1],p1->CAxyz[i+i1][2],p1->CAxyz[i+i2][0],p1->CAxyz[i+i2][1],p1->CAxyz[i+i2][2]);
	   d2=L(p2->CAxyz[j+i1][0],p2->CAxyz[j+i1][1],p2->CAxyz[j+i1][2],p2->CAxyz[j+i2][0],p2->CAxyz[j+i2][1],p2->CAxyz[j+i2][2]);

	   tmp2=fabs(d1-d2);

	   tmp+=(THETA - tmp2/(0.5*(d1+d2)))*exp(-((0.5*(d1+d2)*(0.5*(d1+d2)))/(ALPHA*ALPHA)));
	   //printf("%f %f %f\n",cd1[i+i1][0],cd1[i+i1][1],cd1[i+i1][2]);
	   //printf("%f %f %f\n",cd2[j+i1][0],cd2[j+i1][1],cd2[j+i1][2]);
	   //printf("d1= %.3f %.3f\n",d1,d2);
	  }
	 }
	 //printf("SIM= %.3f %.3f\n",tmp,tmp/(double)(len*len));
	 out[p2->NumOfRes * i + j]=tmp/(double)(len*len);
	}
 }

}

int LocalStrSimWindow(PDB *p1,PDB *p2,int len,double *out){
 int i,n;
 int j,k,i1,i2;
 double d1,d2,tmp,tmp2;
 for(i=0;i<p1->NumOfRes;i++){
 	for(j=0;j<p2->NumOfRes;j++){
 	 tmp=0;
 	 n=0;
	 for(i1=-len;i1<len;i1++){

	  if(i+i1<0 || i+i1 >= p1->NumOfRes)
	   continue;
	  if(j+i1<0 || j+i1 >= p2->NumOfRes)
	   continue;
	  

	  for(i2=-len;i2<len;i2++){
	   if(i+i2<0 || i+i2 >= p1->NumOfRes)
	    continue;
	   if(j+i2<0 || j+i2 >= p2->NumOfRes)
	    continue;
	   n++;	   

	   if(i1==i2){
	    tmp+=THETA;
	    continue;
	   }

	   d1=L(p1->CAxyz[i+i1][0],p1->CAxyz[i+i1][1],p1->CAxyz[i+i1][2],p1->CAxyz[i+i2][0],p1->CAxyz[i+i2][1],p1->CAxyz[i+i2][2]);
	   d2=L(p2->CAxyz[j+i1][0],p2->CAxyz[j+i1][1],p2->CAxyz[j+i1][2],p2->CAxyz[j+i2][0],p2->CAxyz[j+i2][1],p2->CAxyz[j+i2][2]);

	   tmp2=fabs(d1-d2);

	   tmp+=(THETA - tmp2/(0.5*(d1+d2)))*exp(-((0.5*(d1+d2)*(0.5*(d1+d2)))/(ALPHA*ALPHA)));
	  }
	 }
	 out[p2->NumOfRes * i + j]=tmp/(double)(n);
	 //printf("SimWin %d = %.3f\n",p2->NumOfRes * i + j,out[p2->NumOfRes * i + j]);
	}
 }

 return 0;
}
/*
//Ultra Fast!!
int SetSmtx(int *A1,int *A2,DMTX *d1,DMTX *d2,double *smtx){
 int i,j,k,m,id;
 double a1,a2,sco,bnd,t,tmp_bnd,tmp_t;
 double fac=1.0;
 double wc;//wrong contact
 //double base=100;
 //init
 for(i=0;i<d1->Max*d2->Max;i++){
   smtx[i]=0;
 }
 //input
 for(i=0;i<d1->Max;i++){
  for(k=0;k<d1->Max;k++){//aligned positions
   if(d1->area[i][k]==0||A1[k]==-1||i==k)
    continue;
 	for(j=0;j<d2->Max;j++){
		 //remove conflicting position

		 //if((A1[k]>j && k < i) || (A1[k] < j && k > i) )
		 // continue;
		 if((i-k)*(j-A1[k])<=0)
		  continue;
		 if(j==A1[k])
		  continue;
// Target <-> Model symmetry score

		 //M->T
		 a1=d1->area[i][k]+d1->area[k][i];
		 a2=d2->area[j][A1[k]]+d2->area[A1[k]][j];
		 id=j+d2->Max*i;
		 smtx[id]-=cad_bnd(a1,a2);
		 smtx[id]+=a1;
		 //M<-T
		 smtx[id]-=cad_bnd(a2,a1);
		 smtx[id]+=a2;
		 
	}
  }
 }
 return 0;
}
*/
//Ultra Fast!!
int SetSmtx(int *A1,int *A2,DMTX *d1,DMTX *d2,double *smtx){
 int i,j,k,m,id;
 double a1,a2,sco,bnd,t,tmp_bnd,tmp_t;
 //double fac=1.0;
 double wc;//wrong contact
 //double base=100;
 //init
 for(i=0;i<d1->Max*d2->Max;i++){
   smtx[i]=0;
 }
 //input
 for(i=0;i<d1->Max;i++){
  m=d2->Max*i;
  for(k=0;k<d1->Max;k++){//aligned positions
   if(d1->area[i][k]==0||A1[k]==-1||i==k)
    continue;

   a1=d1->area[i][k]+d1->area[k][i];

 	for(j=0;j<d2->Max;j++){
		 //remove conflicting position
/*
i++++------k++++     k++++-----i++++
j++++------A1[k]  or A1[k]-----j++++
*/

		 //if((A1[k]>j && k < i) || (A1[k] < j && k > i) )
		 // continue;
		 if((i-k)*(j-A1[k])<=0)
		  continue;
		 //if(j==A1[k])
		 // continue;
// Target <-> Model symmetry score

		 //M->T
		 //a1=d1->area[i][k]+d1->area[k][i];
		 a2=d2->area[j][A1[k]]+d2->area[A1[k]][j];
		 id=j+m;
		 smtx[id]+=cad_bnd_sym(a1,a2);
	}
	//  smtx[j+d2->Max*i]=t-bnd;
  }
 }
 return 0;
}

int SetSmtxSepWeight(int *A1,int *A2,DMTX *d1,DMTX *d2,double *smtx,double w){
 int i,j,k,m,id;
 double a1,a2,sco,bnd,t,tmp_bnd,tmp_t;
 double wc;//wrong contact
 double near1,near2;
 //init
 for(i=0;i<d1->Max*d2->Max;i++){
   smtx[i]=0;
 }
 //input
 for(i=0;i<d1->Max;i++){
  m=d2->Max*i;
  for(k=0;k<d1->Max;k++){//aligned positions
   if(d1->area[i][k]==0||A1[k]==-1||i==k)
    continue;

   a1=d1->area[i][k]+d1->area[k][i];

   if(i-k < SEP_LEN && k-i < SEP_LEN)
    near1=0;//near
   else
    near1=1;//not near

 	for(j=0;j<d2->Max;j++){
		 //remove conflicting position
		 if((i-k)*(j-A1[k])<1)
		  continue;
		 a2=d2->area[j][A1[k]]+d2->area[A1[k]][j];

		 if(a2==0)
		  continue;
		 id=j+m;

		 if(j-A1[k]<SEP_LEN && A1[k]-j<SEP_LEN)
		  near2=0;//near
		 else
		  near2=1;//not near


		 if(near1*near2==1){//not near & not near
		  smtx[id]+=cad_bnd_sym(a1,a2);
		 }else{//including near
		  smtx[id]+=w*cad_bnd_sym(a1,a2);
		  //printf("Not near %d %d %d %d %f %f %f\n",i,k,j,A1[k],a1,a2,smtx[id]);
		 }
	}
  }
 }
 return 0;
}



//Use pre filter
int SmtxFilter(double *in,double *out,int n1,int n2,int L){
 int i,j,k,m,id;
 double fac=(double)(L*2+1);
 //double r[100];
 //make filter
 //printf("filter %d\n",L);

 //convert
 for(i=0;i<n1;i++){
  m=n2*i;
  for(j=0;j<n2;j++){
   id=j+m;
   //out[id]=in[id]*G_FILTER[L][0];
   //out[id]=in[id]/fac;
   out[id]=in[id];
   for(k=1;k<=L && i+k < n1 && j+k<n2 ; k++){
    //out[id]+=in[(j+k)+n2*(i+k)]*G_FILTER[L][k];
    //out[id]+=in[(j+k)+n2*(i+k)]/fac;
    out[id]+=in[(j+k)+n2*(i+k)];
    //printf(">>%d %d %f\n",L,k,G_FILTER[L][k]);
   }
    //out[id]+=in[(j+k)+n2*(i+k)]/fac;
   for(k=1;k<=L && i-k > -1 && j-k > -1 ; k++)
    //out[id]+=in[(j-k)+n2*(i-k)]*G_FILTER[L][k];
    //out[id]+=in[(j-k)+n2*(i-k)]/fac;
    out[id]+=in[(j-k)+n2*(i-k)];
  }
 }
}

int SetSmtx_old(int *A1,int *A2,DMTX *d1,DMTX *d2,double *smtx){
 int i,j,k;
 double a1,a2,sco,bnd,t,tmp_bnd,tmp_t;
 double fac=1.0;
 double wc;//wrong contact
 //double base=100;
 for(i=0;i<d1->Max;i++){
 	for(j=0;j<d2->Max;j++){
	 bnd=t=sco=0;
 	 wc=0;
///*
		for(k=0;k<i;k++){
		 if(A1[k]==-1||d1->area[i][k]==0||d2->area[j][A1[k]]==0)
		  continue;
		 //remove conflicting position
		 if(A1[k]>=j)
		  break;
		 //i->k
		 a1=d1->area[i][k];
		 a2=d2->area[j][A1[k]];
		 bnd+=cad_bnd(a1,a2);
		 //if(a1==0)wc+=a2;
		 t+=a1;
		 //k->i
		 a1=d1->area[k][i];
		 a2=d2->area[A1[k]][j];
		 bnd+=cad_bnd(a1,a2);
		 t+=a1;
		}
		for(k=d1->Max-1;k>i;k--){
		 if(A1[k]==-1||d1->area[i][k]==0||d2->area[j][A1[k]]==0)
		  continue;
		 //remove conflicting position
		 if(A1[k]<=j)
		  break;
		 //i->k
		 a1=d1->area[i][k];
		 a2=d2->area[j][A1[k]];
		 bnd+=cad_bnd(a1,a2);
		 //if(a1==0)wc+=a2;
		 t+=a1;
		 //k->i
		 a1=d1->area[k][i];
		 a2=d2->area[A1[k]][j];
		 bnd+=cad_bnd(a1,a2);
		 t+=a1;
		}
//*/
/*
	 	for(k=0;k<d1->Max;k++){//aligned positions
		 if(A1[k]==-1||i==k||j==A1[k])
		  continue;
		 //remove conflicting position
		 if((A1[k]>j && k < i) || (A1[k] < j && k > i) )
		  continue;
		 //i->k
		 a1=d1->area[i][k];
		 a2=d2->area[j][A1[k]];
		 bnd+=cad_bnd(a1,a2);
		 //if(a1==0)wc+=a2;
		 t+=a1;
		 //k->i
		 a1=d1->area[k][i];
		 a2=d2->area[A1[k]][j];
		 bnd+=cad_bnd(a1,a2);
		 t+=a1;
		 //if(a1==0)wc+=a2;
		}
*/
	  smtx[j+d2->Max*i]=t-bnd;
	}
 }
 return 0;
}


int SetSmtx_sq(int *A1,int *A2,DMTX *d1,DMTX *d2,double *smtx){
 int i,j,k;
 double a1,a2,sco,bnd,t,tmp_bnd,tmp_t;
 for(i=0;i<d1->Max;i++){
 	for(j=0;j<d2->Max;j++){
	 bnd=t=sco=0;
	 	for(k=0;k<d1->Max;k++){//aligned positions
		 
		 if(A1[k]==-1||i==k||j==A1[k])
		  continue;
		 //remove conflicting position
		 if((A1[k]>j && k < i) || (A1[k] < j && k > i) )
		  continue;
		 //i->k
		 a1=d1->area[i][k];
		 a2=d2->area[j][A1[k]];
		 //bnd+=cad_bnd(a1,a2);
		 tmp_bnd=cad_bnd(a1,a2);
		 //bnd+=cad_bnd(a2,a1);//invert
		 //t+=a1;
		 tmp_t=a1;
		 //t+=a2;//inv
		 //k->i
		 a1=d1->area[k][i];
		 a2=d2->area[A1[k]][j];
		 //bnd+=cad_bnd(a1,a2);
		 tmp_bnd+=cad_bnd(a1,a2);
		 //bnd+=cad_bnd(a2,a1);//inv
		 //t+=a1;
		 tmp_t+=a1;
		 //t+=a2;//inv

		 //new
		 sco+=sqrt(tmp_t-tmp_bnd);
		}
	 //sco=t-bnd;//matched area based
	 //sco=1-bnd/t;//matched area rate based
	 //sco=sqrt(t-bnd);
	 //sco=(t-bnd)*(t-bnd);
	 //printf("%d vs %d = %.3f rate= %.3f\n",i,j,sco,1-bnd/t);
	 //if(sco> 0 && 1-bnd/t < 0.2)
	 // smtx[j+d2->Max*i]=0;
	 //else
	  smtx[j+d2->Max*i]=sco;
	}
 }
 return 0;
}

//Fragment Based Smtx
int SetSmtxFrag(int *A1,int *A2,DMTX *d1,DMTX *d2,int n,double *smtx){
 int i,j,k,l,m;
 double a1,a2,sco,bnd,t,tmp_bnd,tmp_t;
 double fac=1.0;
 double wc;//wrong contact
 int len1,len2;
 len1=d1->Max-n;
 len2=d2->Max-n;
 //double base=100;
 for(i=0;i<len1;i++){
 	for(j=0;j<len2;j++){
	 bnd=t=sco=0;
 	 wc=0;

		//contact of inter fragment
	 	for(k=0;k<len1;k++){//aligned positions
		 if(A1[k]==-1||i==k||j==A1[k])
		  continue;
		 //remove conflicting position
		 if((A1[k]>j && k < i) || (A1[k] < j && k > i) )
		  continue;

			//for(l=0;l)
		 //symmetry
		 a1=d1->area[i][k]+d1->area[k][i];
		 a2=d2->area[j][A1[k]]+d2->area[A1[k]][j];
		 bnd+=cad_bnd(a1,a2);
		 bnd+=cad_bnd(a2,a1);
		 t+=a1+a2;
/*
		 a1=d1->area[i][k];
		 a2=d2->area[j][A1[k]];
		 bnd+=cad_bnd(a1,a2);
		 t+=a1;
		 //k->i
		 a1=d1->area[k][i];
		 a2=d2->area[A1[k]][j];
		 bnd+=cad_bnd(a1,a2);
		 t+=a1;
*/
		}
	  smtx[j+len2*i]=t-bnd;
	}
 }
 return 0;
}


//int TrimBadRegions(int *a,DMTX *d1, DMTX *d2,double *LoSim, double rcut){
int TrimBadRegions(int *a,DMTX *d1, DMTX *d2, double rcut){
 int i,j;
 int tid;
 int tmp[1000];
 double prate,rate,sco,min=1000;
 double psco=Ali2Sco(a,d1,d2,&prate);
 double Tarea,bd;
 int p1,p2;
 int tlist[1000],tcnt=0;

/*
 //Local Structure similarity
 for(i=0;i<d1->Max;i++){
  if(a[i]<0)
   continue;
  if(LoSim[d2->Max*i +a[i]]<=0.0)
   a[i]=-1;
 }
*/
 //Area
 for(i=0;i<d1->Max;i++){
  if(a[i]<0)
   continue;
  Tarea=bd=0;
 	for(j=0;j<d1->Max;j++){
	 if(i==j||d1->area[i][j]==0||d1->area[j][i]==0)
	  continue;
	 if(a[j]<0)
 	  continue;

	 p1=a[i];
	 p2=a[j];
	 //i -> j
	 bd+=cad_bnd(d1->area[i][j],d2->area[p1][p2]); 
	 //j -> i	 
	 //bd+=cad_bnd(d1->area[j][i],d2->area[p2][p1]);
	 //Tarea+=d1->area[i][j]+d1->area[j][i];
	 Tarea+=d1->area[i][j];
	}
  if(Tarea==0||bd==0||1-bd/Tarea < rcut){
  //if(Tarea==0||bd==0){
   tlist[tcnt]=i;
   tcnt++;
   //printf("%d Area= %.3f rate= %.3f\n",i,Tarea-bd,1-bd/Tarea);
  }
 }
 //Remove
 for(i=0;i<tcnt;i++)
  a[tlist[i]]=-1;
 return 0;
}

int ShowAli(int *a, int n){
 int i;
 for(i=0;i<n;i++){
  if(a[i]==-1){
   //printf("%d:%d ",i+1,-1);
   continue;
  }else{
   printf("%d:%d ",i+1,a[i]+1);
  }
 }
 puts("");
 return 0;
}

int ShowPymolScript(int *a,int *a2,int n,int n2){
 int j;
 printf("select align1, ** and resi ");

 for(j=0;j<n;j++){
  if(a[j]!=-1 && (j==0 || (j>0 && a[j-1]==-1)))
   printf("%d",j+1); 
  else if(a[j]!=-1 && (j==n-1||a[j+1]==-1))
   printf("-%d+",j+1); 
 } puts("");

 printf("select align2, ** and resi ");
 for(j=0;j<n2;j++){
  if(a2[j]!=-1 && (j==0 || (j>0 && a2[j-1]==-1)))
   printf("%d",j+1); 
  else if(a2[j]!=-1 && (j==n2-1||a2[j+1]==-1))
   printf("-%d+",j+1);
 }
 puts("");
 return 0;
}

int StockCheck(ALIDATA *a,int *b,int len,int n){
 int i,j;
 int flag=0;
 for(i=n-1;i>=0;i--){
  flag=0;
 	for(j=0;j<len;j++){
	 //printf("%d:%d ",a[i].a1[j], b[j]);
	 if(a[i].a1[j] != b[j]){
   	  //printf("Diff %4d vs %4d j=%d i=%d\n",a[i].a1[j],b[j],j,i);
	  flag=1;
	  break;
	 }
	}
  if(flag==0){//find same alignment
   //printf("%4d vs %4d i=%d\n",a[i].a1[0],b[0],i);
   return 1;
  }
 }
 return 0;
}
int StockCheckSco(ALIDATA *a,int *b,double s,int len,int n){
 int i,j;
 int flag=0;
 for(i=n-1;i>=0;i--){
  flag=0;
  if(a[i].sco != s)
   continue;
 	for(j=0;j<len;j++){
	 //printf("%d:%d ",a[i].a1[j], b[j]);
	 if(a[i].a1[j] != b[j]){
   	  //printf("Diff %4d vs %4d j=%d i=%d\n",a[i].a1[j],b[j],j,i);
	  flag=1;
	  break;
	 }
	}
  if(flag==0){//find same alignment
   //printf("%4d vs %4d i=%d s= %.1f %.1f\n",a[i].a1[0],b[0],i,s,a[i].sco);
   return 1;
  }
 }
 return 0;
}


int IterDP(ALI *a,DMTX *d1,DMTX *d2,double *LoSim,double z,double rcut,int len,double gopen,double gext){
 int i,j;
 int n=a->N;
 double df1,df2,df3;
 double asco,rate,orisco;
 double pre;
 int tmp[1000];
 int tmpa1[1000],tmpa2[1000],tmpa3[1000];
 int gal[2000],ng,gsco;
 int besta1[1000],besta2[1000];
 double *smtx,bestrate,bestsco=0;//score matrix
 int cnt=0;
 int bcnt=0;
 int Niter=10;

 ALIDATA *stock;
 int Nstock=0;
 if((stock=(ALIDATA *)malloc(sizeof(ALIDATA)*n*Niter))==NULL)
  return -1;
 for(i=0;i<n*Niter;i++)
  if((stock[i].a1=(int *)malloc(sizeof(int)*d1->Max))==NULL)
   return -1;
   

 if((smtx=(double *)malloc(sizeof(double)*(d1->Max+1)*(d2->Max+1)))==NULL)
  return -1;

 for(i=0;i<n;i++){
  if((a->ali[i].sco-a->ave)/a->std<z)//z-score cut off
   break;
  SetSmtx(a->ali[i].a1,a->ali[i].a2,d1,d2,smtx);
  //full length DP----------------------
  cnt=0;
  bestsco=a->ali[i].sco;
  CopyAli(a->ali[i].a1,a->ali[i].a2,besta1,besta2,d1->Max,d2->Max);

///*
	while(cnt<Niter){
  	 dp(smtx,gopen,gext,tmpa1,d1->Max,tmpa2,d2->Max,gal,&ng);

//
	 //stockcheck
	 if(StockCheck(stock,tmpa1,d1->Max,Nstock)==0){
	  //Add to Stock
	  CopyAliSingle(tmpa1,stock[Nstock].a1,d1->Max);
	  Nstock++;
	 }else{
	  //printf("Find Same Alignments!! %d\n",i);
	  break;
	 }
//
  	 asco=Ali2Sco(tmpa1,d1,d2,&rate);
  	 //asco=Ali2Sco_sq(tmpa1,d1,d2,&rate);
  	 printf("Pre= %.3f Asco= %.3f Rate= %.3f Best= %.3f\n",a->ali[i].sco,asco,rate,bestsco);
	 //keep best alignment
	 if(bestsco<asco){
	  //copy
	  CopyAli(tmpa1,tmpa2,besta1,besta2,d1->Max,d2->Max);
	  bestsco=asco;
	  bestrate=rate;
	 }else{
	  //break;
	 }
	 //Trim bad regions
	 //TrimBadRegions(tmpa1,d1,d2,rcut);
	 //TrimBadRegions(tmpa1,d1,d2,LoSim,0.01);

  	 SetSmtx(tmpa1,tmpa2,d1,d2,smtx);
	 cnt++;
 	}


  //input
  //CopyAli(besta1,besta2,a->a1[bcnt],a->a2[bcnt],d1->Max,d2->Max);
  CopyAli(besta1,besta2,a->ali[i].a1,a->ali[i].a2,d1->Max,d2->Max);
  //a->sco[bcnt]=bestsco;

    
  //printf("BestScore= %.3f Zsco= %.3f ",bestsco,(bestsco-a->ave)/a->std);
  //printf("Asco= %.3f Rate= %.3f CADscore= %.3f\n",a->ali[i].sco,rate,asco/d1->Tarea);
  a->ali[i].sco=bestsco;
  a->ali[i].rate=bestrate;
  //ShowAli(besta1,d1->Max);
  //ShowPymolScript(a->a1[i],a->a2[i],d1->Max,d2->Max);
  bcnt++;
 }
 a->N=bcnt;

 //free(stock);
 return 0;
}

int IterDPfromInitAli(ALI *a,DMTX *d1,DMTX *d2,double gopen,double gext,int Nf,double *smtx,double *tmp_smtx,ALIDATA *stock){
 int i,j;
 int n=a->N;
 double df1,df2,df3;
 double asco,rate,orisco;
 double pre;
 int tmp[RES];
 int tmpa1[RES],tmpa2[RES],tmpa3[RES];
 int gal[RES],ng,gsco;
 int besta1[RES],besta2[RES];
 //double *smtx,*tmp_smtx,bestrate,bestsco=0;//score matrix
 double bestrate,bestsco=0;//score matrix
 double pre_sco;
 int cnt=0;
 int bcnt=0;
 int Niter=10;
 int Nstock=0;

 pre_sco=a->ali[0].sco;
 for(i=0;i<n;i++){

  SetSmtx(a->ali[i].a1,a->ali[i].a2,d1,d2,tmp_smtx);
  SmtxFilter(tmp_smtx,smtx,d1->Max,d2->Max,Nf);

  //full length DP----------------------
  asco=Ali2Sco(a->ali[i].a1,d1,d2,&rate);
  bestsco=asco;
  bestrate=rate;
  CopyAli(a->ali[i].a1,a->ali[i].a2,besta1,besta2,d1->Max,d2->Max);
  cnt=0;
	while(cnt<Niter){
	 
  	 dp(smtx,gopen,gext,tmpa1,d1->Max,tmpa2,d2->Max,gal,&ng);
  	 asco=Ali2Sco(tmpa1,d1,d2,&rate);
	 //stockcheck
	 if(StockCheckSco(stock,tmpa1,asco,d1->Max,Nstock)==0){
	  //Add to Stock
	  CopyAliSingle(tmpa1,stock[Nstock].a1,d1->Max);
	  stock[Nstock].sco=asco;
	  Nstock++;
	 }else{
	  //printf("Find Same Alignments!! %d\n",i);
	  break;
	 }

  	 //asco=Ali2Sco_sq(tmpa1,d1,d2,&rate);
  	 //printf("%d:Asco= %.3f Rate= %.3f Best= %.3f\n",cnt,asco,rate,bestsco);
	 //keep best alignment
	 if(bestsco<asco){
	  //copy
	  CopyAli(tmpa1,tmpa2,besta1,besta2,d1->Max,d2->Max);
	  bestsco=asco;
	  bestrate=rate;
	 }else{
	  //break;
	 }
	 //Trim bad regions
	 //TrimBadRegions(tmpa1,d1,d2,0.1);
	 cnt++;
	 if(cnt>=Niter)
	  break;

  	 SetSmtx(tmpa1,tmpa2,d1,d2,tmp_smtx);
  	 SmtxFilter(tmp_smtx,smtx,d1->Max,d2->Max,Nf);

 	}


  //input or not
  //CopyAli(besta1,besta2,a->a1[bcnt],a->a2[bcnt],d1->Max,d2->Max);
  CopyAli(besta1,besta2,a->ali[i].a1,a->ali[i].a2,d1->Max,d2->Max);
  //a->sco[bcnt]=bestsco;

    
  //printf("BestScore= %.3f Zsco= %.3f ",bestsco,(bestsco-a->ave)/a->std);
  //printf("Asco= %.3f Rate= %.3f CADscore= %.3f\n",a->ali[i].sco,rate,asco/d1->Tarea);
  a->ali[i].sco=bestsco;
  a->ali[i].rate=bestrate;
  //ShowAli(besta1,d1->Max);
  //ShowPymolScript(a->a1[i],a->a2[i],d1->Max,d2->Max);
  bcnt++;
 }
 a->N=bcnt;

 //free(smtx);
 //free(tmp_smtx);
 return 0;
}

int IterDPfromInitAliSepW(ALI *a,DMTX *d1,DMTX *d2,double gopen,double gext,int Nf,double *smtx,double *tmp_smtx,ALIDATA *stock,double w){
 int i,j;
 int n=a->N;
 double df1,df2,df3;
 double asco,rate,orisco;
 double pre;
 int tmp[RES];
 int tmpa1[RES],tmpa2[RES],tmpa3[RES];
 int gal[RES],ng,gsco;
 int besta1[RES],besta2[RES];
 double bestrate,bestsco=0;//score matrix
 double pre_sco;
 int cnt=0;
 int bcnt=0;
 int Niter=10;
 int Nstock=0;

 pre_sco=a->ali[0].sco;
 for(i=0;i<n;i++){
  //printf("#Start %d\n",i);
  //SetSmtx(a->ali[i].a1,a->ali[i].a2,d1,d2,tmp_smtx);
  SetSmtxSepWeight(a->ali[i].a1,a->ali[i].a2,d1,d2,tmp_smtx,w);
  SmtxFilter(tmp_smtx,smtx,d1->Max,d2->Max,Nf);

  //full length DP----------------------
  //asco=Ali2Sco(a->ali[i].a1,d1,d2,&rate);
  asco=Ali2ScoSepW(a->ali[i].a1,d1,d2,&rate,w);
  bestsco=asco;
  bestrate=rate;
  CopyAli(a->ali[i].a1,a->ali[i].a2,besta1,besta2,d1->Max,d2->Max);
  cnt=0;
	while(cnt<Niter){
  	 dp(smtx,gopen,gext,tmpa1,d1->Max,tmpa2,d2->Max,gal,&ng);
	 //ShowAli(tmpa1,d1->Max);
  	 //asco=Ali2Sco(tmpa1,d1,d2,&rate);
  	 asco=Ali2ScoSepW(tmpa1,d1,d2,&rate,w);
	 //printf("Asco= %f\n",asco);
	 //stockcheck
	 if(StockCheckSco(stock,tmpa1,asco,d1->Max,Nstock)==0){
	  //Add to Stock
	  CopyAliSingle(tmpa1,stock[Nstock].a1,d1->Max);
	  stock[Nstock].sco=asco;
	  Nstock++;
	 }else{
	  //printf("Find Same Alignments!! %d\n",i);
	  break;
	 }

	 //keep best alignment
	 if(bestsco<asco){
	  //copy
	  CopyAli(tmpa1,tmpa2,besta1,besta2,d1->Max,d2->Max);
	  bestsco=asco;
	  bestrate=rate;
	 }else{
	  //break;
	 }
	 cnt++;
	 if(cnt>=Niter)
	  break;

  	 //SetSmtx(tmpa1,tmpa2,d1,d2,tmp_smtx);
  	 SetSmtxSepWeight(tmpa1,tmpa2,d1,d2,tmp_smtx,w);
  	 SmtxFilter(tmp_smtx,smtx,d1->Max,d2->Max,Nf);

 	}


  //input or not
  CopyAli(besta1,besta2,a->ali[i].a1,a->ali[i].a2,d1->Max,d2->Max);

    
  a->ali[i].sco=bestsco;
  a->ali[i].rate=bestrate;
  bcnt++;
 }
 a->N=bcnt;

 return 0;
}



int IterDPfromInitAliSubOpt(ALI *a,DMTX *d1,DMTX *d2,double gopen,double gext,int Nf,double *smtx,double *tmp_smtx,ALIDATA *stock){
 int i,j,k;
 int n=a->N;
 double df1,df2,df3;
 double asco,rate,orisco;
 double pre;
 int tmp[RES];
 int tmpa1[RES],tmpa2[RES],tmpa3[RES];
 int gal[RES],ng,gsco;
 int besta1[RES],besta2[RES];
 double bestrate,bestsco=0;//score matrix
 double pre_sco;
 int cnt=0;
 int bcnt=0;
 int Niter=10;
 int Nstock=0;

 //subopt mode
 double sub_best;
 int sub_a1[RES],sub_a2[RES];

 pre_sco=a->ali[0].sco;
 for(i=0;i<n;i++){

  SetSmtx(a->ali[i].a1,a->ali[i].a2,d1,d2,tmp_smtx);
  SmtxFilter(tmp_smtx,smtx,d1->Max,d2->Max,Nf);

  //full length DP----------------------
  asco=Ali2Sco(a->ali[i].a1,d1,d2,&rate);
  bestsco=asco;
  bestrate=rate;
  CopyAli(a->ali[i].a1,a->ali[i].a2,besta1,besta2,d1->Max,d2->Max);
  cnt=0;
	while(cnt<Niter){

	 //subopt mode
	 sub_best=0;
	 for(j=0;j<10;j++){	 
  	  dp(smtx,gopen,gext,sub_a1,d1->Max,sub_a2,d2->Max,gal,&ng);
  	  asco=Ali2Sco(sub_a1,d1,d2,&rate);
	  printf("Cnt= %d %d -> %f\n",cnt,j,asco);

	  if(sub_best<asco){
	   CopyAliSingle(sub_a1,tmpa1,d1->Max);
	   sub_best=asco;
	  }

	  //remove aligned smtx
	  	for(k=0;k<d1->Max;k++){
		 if(sub_a1[k]==-1)
		  continue;
		 //clear
		 smtx[sub_a1[k]+d2->Max*k]=0;
	  	}
	  
	 }
	 asco=sub_best;
	 //-----------------------
	 //stockcheck
	 if(StockCheckSco(stock,tmpa1,asco,d1->Max,Nstock)==0){
	  //Add to Stock
	  CopyAliSingle(tmpa1,stock[Nstock].a1,d1->Max);
	  stock[Nstock].sco=asco;
	  Nstock++;
	 }else{
	  printf("Find Same Alignments!! %d\n",i);
	  break;
	 }

	 //keep best alignment
	 if(bestsco<asco){
	  //copy
	  CopyAli(tmpa1,tmpa2,besta1,besta2,d1->Max,d2->Max);
	  bestsco=asco;
	  bestrate=rate;
	 }else{
	  //break;
	 }
	 cnt++;
	 if(cnt>=Niter)
	  break;

  	 SetSmtx(tmpa1,tmpa2,d1,d2,tmp_smtx);
  	 SmtxFilter(tmp_smtx,smtx,d1->Max,d2->Max,Nf);

 	}


  //input or not
  //CopyAli(besta1,besta2,a->a1[bcnt],a->a2[bcnt],d1->Max,d2->Max);
  CopyAli(besta1,besta2,a->ali[i].a1,a->ali[i].a2,d1->Max,d2->Max);
  //a->sco[bcnt]=bestsco;

    
  //printf("BestScore= %.3f Zsco= %.3f ",bestsco,(bestsco-a->ave)/a->std);
  //printf("Asco= %.3f Rate= %.3f CADscore= %.3f\n",a->ali[i].sco,rate,asco/d1->Tarea);
  a->ali[i].sco=bestsco;
  a->ali[i].rate=bestrate;
  //ShowAli(besta1,d1->Max);
  //ShowPymolScript(a->a1[i],a->a2[i],d1->Max,d2->Max);
  bcnt++;
 }
 a->N=bcnt;

 return 0;
}



int IterDPfromInitAliFrag(ALI *a,DMTX *d1,DMTX *d2,double gopen,double gext,int Nf){
 int i,j;
 int n=a->N;
 double df1,df2,df3;
 double asco,rate,orisco;
 double pre;
 int tmp[1000];
 int tmpa1[1000],tmpa2[1000],tmpa3[1000];
 int gal[2000],ng,gsco;
 int besta1[1000],besta2[1000];
 double *smtx,bestrate,bestsco=0;//score matrix
 int cnt=0;
 int bcnt=0;
 int Niter=10;

 ALIDATA *stock;
 int Nstock=0;
 if((stock=(ALIDATA *)malloc(sizeof(ALIDATA)*n*Niter))==NULL)
  return -1;
 for(i=0;i<n*Niter;i++)
  if((stock[i].a1=(int *)malloc(sizeof(int)*d1->Max))==NULL)
   return -1;
   
 if((smtx=(double *)malloc(sizeof(double)*(d1->Max+1)*(d2->Max+1)))==NULL)
  return -1;

 for(i=0;i<n;i++){
  SetSmtxFrag(a->ali[i].a1,a->ali[i].a2,d1,d2,Nf,smtx);
  //SetSmtx(a->ali[i].a1,a->ali[i].a2,d1,d2,smtx);
  //full length DP----------------------
  cnt=0;
  bestsco=a->ali[i].sco;
  CopyAli(a->ali[i].a1,a->ali[i].a2,besta1,besta2,d1->Max,d2->Max);

	while(cnt<Niter){
  	 dp(smtx,gopen,gext,tmpa1,d1->Max,tmpa2,d2->Max,gal,&ng);

	 //stockcheck
	 if(StockCheck(stock,tmpa1,d1->Max,Nstock)==0){
	  //Add to Stock
	  CopyAliSingle(tmpa1,stock[Nstock].a1,d1->Max);
	  Nstock++;
	 }else{
	  //printf("Find Same Alignments!! %d\n",i);
	  break;
	 }
  	 asco=Ali2Sco(tmpa1,d1,d2,&rate);
  	 printf("Pre= %.3f Asco= %.3f Rate= %.3f Best= %.3f\n",a->ali[i].sco,asco,rate,bestsco);
	 //keep best alignment
	 if(bestsco<asco){
	  //copy
	  CopyAli(tmpa1,tmpa2,besta1,besta2,d1->Max,d2->Max);
	  bestsco=asco;
	  bestrate=rate;
	 }else{
	  //break;
	 }
	 //Trim bad regions
	 //TrimBadRegions(tmpa1,d1,d2,0.1);
  	 SetSmtx(tmpa1,tmpa2,d1,d2,smtx);
	 cnt++;
 	}


  //input
  //CopyAli(besta1,besta2,a->a1[bcnt],a->a2[bcnt],d1->Max,d2->Max);
  CopyAli(besta1,besta2,a->ali[i].a1,a->ali[i].a2,d1->Max,d2->Max);
  //a->sco[bcnt]=bestsco;

    
  //printf("BestScore= %.3f Zsco= %.3f ",bestsco,(bestsco-a->ave)/a->std);
  //printf("Asco= %.3f Rate= %.3f CADscore= %.3f\n",a->ali[i].sco,rate,asco/d1->Tarea);
  a->ali[i].sco=bestsco;
  a->ali[i].rate=bestrate;
  //ShowAli(besta1,d1->Max);
  //ShowPymolScript(a->a1[i],a->a2[i],d1->Max,d2->Max);
  bcnt++;
 }
 a->N=bcnt;

 //free(stock);
 return 0;
}


int ShowResults(ALI *a,DMTX *d1,DMTX *d2,double z,double rcut,int len,PDB *p1,PDB *p2){
 int i,j;
 int n=a->N;
 double df1,df2,df3;
 double asco,asco4,asco_all,rate,rate4,rate_all,orisco,rms;
 double pre;
 int tmp[1000];
 int alen;

 for(i=0;i<n && i < 1;i++){
 //for(i=0;i<n ;i++){
  asco=Ali2Sco(a->ali[i].a1,d1,d2,&rate);
  asco=Ali2ScoSep(a->ali[i].a1,d1,d2,&rate,4);
  //printf("No%d sco= %.3f rate= %.3f CAD= %.3f ",i+1,a->ali[i].sco,a->ali[i].rate,0.5*a->ali[i].sco/d1->Tarea);
  printf("No%d sco= %.3f rate= %.3f CAD= %.3f ",i+1,a->ali[i].sco,a->ali[i].rate,0.5*a->ali[i].sco/d1->Tarea);
  printf("Sep>4 sco= %.3f rate= %.3f CAD= %.3f\n",asco,rate,(asco/d1->Tarea4)*0.50);
  rms=CalAliRmsd(a->ali[i].a1,a->len1,p1,p2,&alen);
  printf("RMSD= %.3f (%d Res)\n",rms,alen);

  //ShowAli(a->ali[i].a1,a->len1);
  //resiude coord based alignments
  ShowAliCode(a->ali[i].a1,a->ali[i].a2,a->len1,a->len2,p1,p2);
 }

}

int ShowResultsSepW(ALI *a,DMTX *d1,DMTX *d2,double z,double rcut,int len,PDB *p1,PDB *p2,double w,bool LQmode){
 int i,j;
 int n=a->N;
 double df1,df2,df3;
 double asco,asco4,asco_all,rate,rate4,rate_all,orisco,rms;
 double pre;
 int tmp[1000];
 int alen;

 for(i=0;i<n && i < 1;i++){
 //for(i=0;i<n ;i++){
  asco=Ali2Sco(a->ali[i].a1,d1,d2,&rate);
  asco4=Ali2ScoSep(a->ali[i].a1,d1,d2,&rate4,SEP_LEN-1);
  asco_all=Ali2ScoSepW(a->ali[i].a1,d1,d2,&rate_all,w);
  //printf("No%d sco= %.3f rate= %.3f CAD= %.3f ",i+1,a->ali[i].sco,a->ali[i].rate,0.5*a->ali[i].sco/d1->Tarea);
  printf("No%d sco= %.3f rate= %.3f CAD= %.3f ",i+1,asco,rate,0.5*asco/d1->Tarea);
  printf("Sep>4 sco= %.3f rate= %.3f CAD= %.3f ",asco4,rate4,(asco4/d1->Tarea4)*0.50);
  printf("Weighted sco= %.3f rate= %.3f CAD= %.3f ",asco_all,rate_all,(asco_all/(w*(d1->Tarea-d1->Tarea4) + d1->Tarea4))*0.50);
  printf("SupRec= %.3f ",asco_all/pow(d1->Tarea,0.7) + asco_all/pow(d2->Tarea,0.7));
  rms=CalAliRmsd(a->ali[i].a1,a->len1,p1,p2,&alen);
  printf("RMSD= %.3f (%d Res)\n",rms,alen);

  //ShowAli(a->ali[i].a1,a->len1);
  //resiude coord based alignments
  ShowAliCode(a->ali[i].a1,a->ali[i].a2,a->len1,a->len2,p1,p2);
  //Show Local Quality mode
  if(LQmode==true)
   ShowLocalQ(a->ali[i].a1,a->ali[i].a2,a->len1,a->len2,d1,d2);
 }

}
int ShowLocalQ(int *a1,int *a2,int n1,int n2,DMTX *d1,DMTX *d2){
 int i,j,k;
 int p1,p2; 
 int g[RES],Ng=0;
 double base1[RES],bnd1[RES];
 double base2[RES],bnd2[RES];
 //generate global alignment
 i=j=0;
 int now=0;

 double Tarea,TareaAlign,bd,cad;
 double A1,A2;
 Tarea=TareaAlign=bd=0;
 //init
 for(i=0;i<d1->Max;i++){
  base1[i]=0;
  bnd1[i]=0;
 }
/*
 for(i=0;i<d2->Max;i++){
  base2[i]=0;
  bnd2[i]=0;
 }
*/
 for(i=0;i<d1->Max;i++){
  p1=a1[i];
  if(p1<0)
   continue;
 	for(j=i+1;j<d1->Max;j++){
	 p2=a1[j];
	 //no aligned res
	 if(p2 < 0)
	  continue;
	 
	 A1=d1->area[i][j]+d1->area[j][i];
	 A2=d2->area[p1][p2]+d2->area[p2][p1];

	 //bd+=cad_bnd(A1,A2);
	 //bd+=cad_bnd(A2,A1);
	 base1[i]+=A1+A2;
	 base1[j]+=A1+A2;
	 bnd1[i]+=cad_bnd(A1,A2);
	 bnd1[j]+=cad_bnd(A2,A1);
/*
	 base2[p1]+=A1+A2;
	 base2[p2]+=A1+A2;
	 bnd2[p1]+=cad_bnd(A1,A2);
	 bnd2[p2]+=cad_bnd(A2,A1);
*/
	}
 }
 //show digit data
 for(i=0;i<n1;i++){
  if(a1[i]==-1 ){
   g[Ng*2  ]=i;
   g[Ng*2+1]=-1;
   Ng++;
  }
  if(a1[i]!=-1){
	for(j=now;j<a1[i];j++){
	 g[Ng*2  ]=-1;
  	 g[Ng*2+1]=j;
	 Ng++;
	}
   g[Ng*2  ]=i;
   g[Ng*2+1]=a1[i];
   Ng++;
   now=a1[i]+1;
  }
 }
 if(now<n2){
  	for(j=now;j<n2;j++){
	 g[Ng*2  ]=-1;
  	 g[Ng*2+1]=j;
   	 //printf("*%3d %4d%4d\n",Ng,g[Ng*2],g[Ng*2+1]);
	 Ng++;
	}
 }

 int s;
 for(i=0;i<Ng;i++){
   //printf("X",g[i*2]);
  if(g[i*2]!= -1 && g[i*2+1]!= -1){
   if(base1[g[i*2]]==0)
    s=0;
   else
    s=(int)(10.00*(base1[g[i*2]]-bnd1[g[i*2]])/base1[g[i*2]]);
   if(s==10)
    s=9;
   //printf("%f %f \n",base1[g[i*2]],bnd1[g[i*2]]);
   printf("%d",s);
  }else{
   printf("-");
  }
 }
 puts("");
/*
 for(i=0;i<Ng;i++){
  if(g[i*2]!= -1 && g[i*2+1]!= -1){
   if(base2[g[i*2+1]]==0)
    s=0;
   else
    s=(int)(10.00*(base2[g[i*2+1]]-bnd2[g[i*2+1]])/base2[g[i*2+1]]);
   if(s==10)
    s=9;
   //printf("%f %f \n",base1[g[i*2]],bnd1[g[i*2]]);
   printf("%d",s);
  }else{
   printf("-");
  }
 }
 puts("");
*/
 return 0;
}

//Show basic 1 char amino acid code
int ShowAliCode(int *a1,int *a2,int n1,int n2,PDB *p1,PDB *p2){
 int i,j,k;
 int g[RES],Ng=0;;
 //generate global alignment
 i=j=0;
 int now=0;
 for(i=0;i<n1;i++){
  if(a1[i]==-1 ){
   g[Ng*2  ]=i;
   g[Ng*2+1]=-1;
   //printf("%3d %4d%4d\n",Ng,g[Ng*2],g[Ng*2+1]);
   Ng++;
  }
  if(a1[i]!=-1){
   //printf("now= %d here is %d\n",now,a1[i]);
	for(j=now;j<a1[i];j++){
	 g[Ng*2  ]=-1;
  	 g[Ng*2+1]=j;
   	 //printf("*%3d %4d%4d\n",Ng,g[Ng*2],g[Ng*2+1]);
	 Ng++;
	}
   g[Ng*2  ]=i;
   g[Ng*2+1]=a1[i];
   //printf("@@%3d %4d%4d\n",Ng,g[Ng*2],g[Ng*2+1]);
   Ng++;
   now=a1[i]+1;
  }
 }
 if(now<n2){
  	for(j=now;j<n2;j++){
	 g[Ng*2  ]=-1;
  	 g[Ng*2+1]=j;
   	 //printf("*%3d %4d%4d\n",Ng,g[Ng*2],g[Ng*2+1]);
	 Ng++;
	}
 }
  //printf("%3d %4d%4d %d %d\n",i+1,g[i*2],g[i*2+1],p1->TypeResId[g[i*2]],p2->TypeResId[g[i*2+1]]);
 for(i=0;i<Ng;i++){
  if(g[i*2]!= -1)
   printf("%c",AACODE[p1->TypeResId[g[i*2]]]);
  else
   printf("-",AACODE[p1->TypeResId[g[i*2]]]);
 }
 puts("");
 for(i=0;i<Ng;i++){
  if(g[i*2+1]!= -1)
   printf("%c",AACODE[p2->TypeResId[g[i*2+1]]]);
  else
   printf("-",AACODE[p2->TypeResId[g[i*2+1]]]);
 }
 puts("");
 return 0;
}

double CalAliRmsd(int *a,int n,PDB *p1,PDB *p2,int *alen){
 int i,cnt=0;
 double xyz1[RES][3],xyz2[RES][3];

 //input CA xyz coords
 for(i=0;i<n;i++){
  if(a[i]==-1)
   continue;
  xyz1[cnt][0]=p1->CAxyz[i][0];
  xyz1[cnt][1]=p1->CAxyz[i][1];
  xyz1[cnt][2]=p1->CAxyz[i][2];

  xyz2[cnt][0]=p2->CAxyz[a[i]][0];
  xyz2[cnt][1]=p2->CAxyz[a[i]][1];
  xyz2[cnt][2]=p2->CAxyz[a[i]][2];

  cnt++;
 }
 *alen=cnt;
 double rms;
 fast_rmsd(xyz1,xyz2,cnt,&rms);
 return rms;
}

int Malloc_ALI(ALI *a){
 int i;
 int N=a->N;
 if((a->ali=(ALIDATA *)malloc(sizeof(ALIDATA)*N))==NULL)
  return -1;

 for(i=0;i<N;i++){
  if((a->ali[i].a1=(int *)malloc(sizeof(int)*a->len1))==NULL)
   return -1;
  if((a->ali[i].a2=(int *)malloc(sizeof(int)*a->len2))==NULL)
   return -1;
 }
 return 0;
}


int LocalPair(int *a,double *mmtx,DMTX *d1,DMTX *d2,int N){
 int i,j,k,l,id1,id2;
 int m,o;
 bool f;
 double a1,a2,bd,Tarea,sco,max=0;
/*
i+++--------j+++++
k+++--------l+++++
*/

 int cnt=0;
 for(i=0;i<d1->Max-N+1;i++){
  for(k=0;k<d2->Max-N+1;k++){
   if(mmtx[k+d2->Max*i]<0.1)
    continue;
	for(j=i+N+1;j<d1->Max-N+1;j++){

   	 f=false;
	 for(m=0;m<N;m++){
	  for(o=0;o<N;o++){
	   if(d1->area[i+m][j+o]>0){
	    f=true;
	    break;
	   }
	   if(f==true)
	    break;
	  }
	 }

	 if(f==false)
	  continue;

  	 for(l=k+N+1;l<d2->Max-N+1;l++){
   	  if(mmtx[l+d2->Max*j]<0.1)
	   continue;

	  bd=Tarea=0;
	  	for(m=0;m<N;m++){
	  	 for(o=0;o<N;o++){
	  	  if(d2->area[k+m][l+o]==0)
		   continue;
		  a1=d1->area[i+m][j+o]+d1->area[j+o][i+m];
	 	  a2=d2->area[k+m][l+o]+d2->area[l+o][k+m];

		  bd+=cad_bnd(a1,a2);
		  bd+=cad_bnd(a2,a1);
		  Tarea+=a1+a2;

	  	  
	  	 }
	 	}

	  sco=Tarea-bd;


	  cnt++;
	  if(max<sco){
	   max=sco;
	  }
	  if(sco>300)
	   printf("Cmp %d %d %d %d %f\n",i,j,k,l,sco);
	 }
	}
  }
 }


 printf("cnt= %d\n",cnt);
 return 0;
}

int SetFmtx(DMTX *d1,DMTX *d2,double w,double *smtx){
 int i,j,k;
 double a1,a2,sco,bnd,t,tmp_bnd,tmp_t;
 double fac=1.0;
 double wc;//wrong contact
 double af1[RES],af2[RES],bf1[RES],bf2[RES];

 //init
 for(i=0;i<d1->Max;i++) af1[i]=bf1[i]=0;
 for(i=0;i<d2->Max;i++) af2[i]=bf2[i]=0;

 //input
 for(i=0;i<d1->Max;i++){
  for(j=i+1;j<d1->Max;j++){
   if(j-i < SEP_LEN){
    //af1[i]+=(d1->area[i][j]+d1->area[j][i])*w;
    af1[i]+=(d1->area[i][j])*w;
   }else{
    //af1[i]+=(d1->area[i][j]+d1->area[j][i]);
    af1[i]+=(d1->area[i][j]);
   }
  }
  for(j=i-1;j>=0;j--){
   if(i-j < SEP_LEN){
    //bf1[i]+=(d1->area[i][j]+d1->area[j][i])*w;
    bf1[i]+=(d1->area[i][j])*w;
   }else{
    //bf1[i]+=(d1->area[i][j]+d1->area[j][i]);
    bf1[i]+=(d1->area[i][j]);
   }
  }
 }

 for(i=0;i<d2->Max;i++){
  for(j=i+1;j<d2->Max;j++){
   if(j-i < SEP_LEN){
    //af2[i]+=(d2->area[i][j]+d2->area[j][i])*w;
    af2[i]+=(d2->area[i][j])*w;
   }else{
    //af2[i]+=(d2->area[i][j]+d2->area[j][i]);
    af2[i]+=(d2->area[i][j]);
   }
  }
  for(j=i-1;j>=0;j--){
   if(i-j < SEP_LEN){
    //bf2[i]+=(d2->area[i][j]+d2->area[j][i])*w;
    bf2[i]+=(d2->area[i][j])*w;
   }else{
    //bf2[i]+=(d2->area[i][j]+d2->area[j][i]);
    bf2[i]+=(d2->area[i][j]);
   }
  }
 }

 //input smtx
 int id=0;
 for(i=0;i<d1->Max;i++){
 	for(j=0;j<d2->Max;j++){
	 id=j+d2->Max*i;
	 if(bf1[i] < bf2[j])
 	  smtx[id]=bf1[i];
	 else
 	  smtx[id]=bf2[j];

	 if(af1[i] < af2[j])
 	  smtx[id]+=af1[i];
	 else
 	  smtx[id]+=af2[j];

	 //printf("%d %d (%.1f %.1f) (%.1f %.1f)\n",i,j,bf1[i],af1[i],bf2[j],bf2[j]);
	}
 }
 return 0;
}



