#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
//#include <bits/nan.h>
#include <time.h>
#include <sys/time.h>
#include "struct.h"
#include "func.h"
#include "mammoth.h"
#include "rmsd.h"
#define DIST(a,b,c,d,e,f) sqrt((a-b)*(a-b)+(c-d)*(c-d)+(e-f)*(e-f))
#define L(a,b,c,d,e,f) (sqrt((d-a)*(d-a)+(e-b)*(e-b)+(f-c)*(f-c)))

double mammoth(int *a1,int *a2,double **cd1,int n1,double **cd2,int n2,int Nv,double Gopen,double Gext,double *smtx,double *rsco){
 double sco;
 //double Uv1[RES][Nv][3];
 //double Uv2[RES][Nv][3];
 int i1,i2;
 int len1=n1-Nv;
 int len2=n2-Nv;
 double l,rmsd,rmsd2;
 //DP
 int gal[RES],ng,i;
 int A1[RES],A2[RES];
 int cent=(int)(Nv/2);
 //sco=dp(smtx,Gopen,Gext,a1,len1,a2,len2,gal,&ng);
 //printf("1 sco= %f %d\n",sco, RES);
 //*rsco=dp(smtx,Gopen,Gext,A1,len1,A2,len2,gal,&ng);
 //*rsco=dp_afp(smtx,Gopen,Gext,A1,len1,A2,len2,Nv+1,gal,&ng);//non over-lapped DP

 *rsco=dp_afp(smtx,Gopen,Gext,A1,n1,A2,n2,Nv+1,gal,&ng);//non over-lapped DP

 //suboptimal dp mode
 

 //convert to residue based alignment
/*
if Nv=6

i->i+1->i+2->i+3->i+4->i+5->i+6
             ***
             Center
*/
/*
 //init
 for(i=0;i<n1;i++) a1[i]=-1;
 for(i=0;i<n2;i++) a2[i]=-1;
*/
 for(i=0;i<n1;i++) a1[i]=A1[i];
 for(i=0;i<n2;i++) a2[i]=A2[i];
 
 
/*
 for(i=0;i<len1;i++){
  if(A1[i]<0)
   continue;
  a1[i+cent]=A1[i]+cent;
 }
 
 for(i=0;i<len2;i++){
  if(A2[i]<0)
   continue;
  a2[i+cent]=A2[i]+cent;
 }
*/

 //printf("2 sco= %f\n",*rsco);
 return (sco);
}

double mammoth_subopt(int **a1,int **a2,int N,int n1,int n2,int Nv,double Gopen,double Gext,double *smtx,double *rsco){
 double sco;
 int i1,i2;
 //int len1=n1-Nv;
 //int len2=n2-Nv;
 double l,rmsd,rmsd2;
 //DP
 int gal[RES],ng,i,j;
 int A1[RES],A2[RES];
 int cent=(int)(Nv/2);

 //suboptimal dp mode
 for(i=0;i<N;i++){
  //remove smtx where aligned positon
  if(i>0){
   for(j=0;j<n1;j++){
    if(a1[i-1][j]<0)
     continue;
    smtx[a1[i-1][j]+n2*j]*=0.9;
    //j+=Nv;
   }
  }
  //*rsco=dp_afp(smtx,Gopen,Gext,a1[i],n1,a2[i],n2,Nv+1,gal,&ng);//non over-lapped DP
  *rsco=dp(smtx,Gopen,Gext,a1[i],n1,a2[i],n2,gal,&ng);//allow over-lapped DP
  //printf("Round %d sco= %f\n",i+1,*rsco);
 }

 //convert
 for(i=0;i<N;i++){
  for(j=0;j<n1;j++)
   A1[j]=-1;
  for(j=0;j<n1-Nv;j++){
   if(a1[i][j]<0)
    continue;
   else
    A1[j+cent]=a1[i][j]+cent;
  }
  for(j=0;j<n1;j++){
   //printf("%d <- %d\n",a1[i][j],A1[j]);
   a1[i][j]=A1[j];
  }

  for(j=0;j<n2;j++)
   A2[j]=-1;
  for(j=0;j<n2-Nv;j++){
   if(a2[i][j]<0)
    continue;
   else
    A2[j+cent]=a2[i][j]+cent;
  }
  for(j=0;j<n2;j++)
   a2[i][j]=A2[j];
 }

 return (sco);
}

double mammoth_subopt_fwd(int **a1,int **a2,int N,
		int n1,int n2,int Nv,
		double Gopen,double Gext,double *smtx,double *rsco,
		DMTX *d1,DMTX *d2, double weight){
 double sco;
 int i1,i2;
 //int len1=n1-Nv;
 //int len2=n2-Nv;
 double l,rmsd,rmsd2;
 //DP
 int gal[RES],ng,i,j;
 int A1[RES],A2[RES];

 //suboptimal dp mode
 for(i=0;i<N;i++){
  //remove smtx where aligned positon
  if(i>0){
   for(j=0;j<n1;j++){
    if(a1[i-1][j]<0)
     continue;
    smtx[a1[i-1][j]+n2*j]=0.0;
    //j+=Nv;
   }
  }
  //*rsco=dp_afp(smtx,Gopen,Gext,a1[i],n1,a2[i],n2,Nv+1,gal,&ng);//non over-lapped DP
  *rsco=dp_afp_fwd(smtx,Gopen,Gext,a1[i],n1,a2[i],n2,Nv+1,gal,&ng,d1,d2,weight);//non over-lapped DP
  //printf("Round %d sco= %f\n",i+1,*rsco);
 }
 return (sco);
}



int setup_mammoth(double **cd1,int n1,double **cd2,int n2,int Nv,double *smtx){
 double sco=0;
 double Uv1[RES][Nv][3];
 double Uv2[RES][Nv][3];
 int i1,i2,tag;
 //int len1=n1-Nv;
 int len1=n1;
 //int len2=n2-Nv;
 int len2=n2;
 double l,rmsd,rmsd2;
 //6 vectors
 //set unit vector
 for(i1=0;i1<n1-Nv;i1++){
  for(i2=0;i2<Nv;i2++){
   tag=i1+i2;
   l=L(cd1[tag][0],cd1[tag][1],cd1[tag][2],cd1[tag+1][0],cd1[tag+1][1],cd1[tag+1][2]);
   Uv1[i1][i2][0]=(cd1[tag+1][0]-cd1[tag][0])/l;
   Uv1[i1][i2][1]=(cd1[tag+1][1]-cd1[tag][1])/l;
   Uv1[i1][i2][2]=(cd1[tag+1][2]-cd1[tag][2])/l;
  }
 }
 for(i1=0;i1<n2-Nv;i1++){
  for(i2=0;i2<Nv;i2++){
   tag=i1+i2;
   l=L(cd2[tag][0],cd2[tag][1],cd2[tag][2],cd2[tag+1][0],cd2[tag+1][1],cd2[tag+1][2]);
   Uv2[i1][i2][0]=(cd2[tag+1][0]-cd2[tag][0])/l;
   Uv2[i1][i2][1]=(cd2[tag+1][1]-cd2[tag][1])/l;
   Uv2[i1][i2][2]=(cd2[tag+1][2]-cd2[tag][2])/l;
  }
 }
 //set URMS(r)
 double urms_r=sqrt(2.0-2.84/sqrt(Nv));
 double urms;
 int m;
 //double *smtx;
 //if((smtx=(double *)malloc(sizeof(double)*len1*len2))==NULL)
 // return -1;
 //printf("URMS(r)= %.3f\n",urms_r);
 //set score matrix
 for(i1=0;i1<n1-Nv;i1++){
  m=i1*len2;
  for(i2=0;i2<n2-Nv;i2++){
   fast_rmsd_noshift(Uv1[i1],Uv2[i2],Nv,&urms);
   //fast_rmsd(Uv1[i1],Uv2[i2],Nv,&rmsd);
   //printf("%d %d rmsd= %.3f\n",i1,i2,urms);
   //smtx[j+d2->Max*i]=t-bnd;
   if(urms_r > urms)
    smtx[i2+m]=(urms_r-urms)*10.00/urms_r;
   else
    smtx[i2+m]=0.00;
   //if(i1==i2)
   //printf("%d vs %d noshift_rmsd= %f s= %.2f\n",i1,i2,urms,smtx[i2+i1*len2]);
  }
 }
 puts("#Fin setup mammoth");
 return 0;
}
