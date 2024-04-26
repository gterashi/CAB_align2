#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "struct.h"
#include "func.h"
#include "tm.h"
#include "rmsd.h"
#include "dali.h"
#define DIST(a,b,c,d,e,f) sqrt((a-b)*(a-b)+(c-d)*(c-d)+(e-f)*(e-f))
#define L(a,b,c,d,e,f) (sqrt((d-a)*(d-a)+(e-b)*(e-b)+(f-c)*(f-c)))
#define L2(a,b,c,d,e,f) ((d-a)*(d-a)+(e-b)*(e-b)+(f-c)*(f-c))

static void cross(double a[3], double b[3], double c[3])
{
  a[0] = b[1]*c[2] - b[2]*c[1];
  a[1] = b[2]*c[0] - b[0]*c[2];
  a[2] = b[0]*c[1] - b[1]*c[0];
}

double TMscore(double cd1[][3],double cd2[][3],int n,double d0,double mov_com[3],double mov_to_ref[3],double rmtx[3][3]){
 double sco=0,pre_sco,max_sco=0;
 double xyz1[RES][3],xyz2[RES][3];
 double f1[RES][3],f2[RES][3],tmp_rmtx[3][3];

 int L[4]={n,n/2,n/4,4};
 int cnt,i,i1,i2;

 int p,len,pos;
 double dist,rmsd,tmp_m_c[3],tmp_m_r[3];
 double d02=d0*d0;

 for(p=0;p<4;p++){
  len=L[p];
  //Shift window From N to C
  for(pos=0;pos+len-1<n;pos++){
   //printf("L= %d pos= %d\n",len,pos);
 	//input CA xyz coords
 	for(i=0;i<len;i++){
	 i1=pos+i;
  	 f1[i][0]=cd1[i1][0];
  	 f1[i][1]=cd1[i1][1];
  	 f1[i][2]=cd1[i1][2];
  
	 f2[i][0]=cd2[i1][0];
	 f2[i][1]=cd2[i1][1];
	 f2[i][2]=cd2[i1][2];

 	}
 	//rotate with first alignment
 	calculate_rotation_rmsd(f1,f2,len,tmp_m_c,tmp_m_r,tmp_rmtx,&rmsd);
	
 	//cal moved cd
 	rotate_shift_cd(cd2,xyz2,n,tmp_m_c,tmp_m_r,tmp_rmtx);

	pre_sco=0;
	//iterative
	while(1){
	 pre_sco=sco;
	 sco=0;
	 cnt=0;
		

   		for(i=0;i<n;i++){
   		 dist=L2(cd1[i][0],cd1[i][1],cd1[i][2],xyz2[i][0],xyz2[i][1],xyz2[i][2]);
		 if(dist<d02){
		  f1[cnt][0]=cd1[i][0];
  		  f1[cnt][1]=cd1[i][1];
  		  f1[cnt][2]=cd1[i][2];
  
		  f2[cnt][0]=cd2[i][0];
		  f2[cnt][1]=cd2[i][1];
		  f2[cnt][2]=cd2[i][2];
		  //printf("%d %f %f\n",i,f1[cnt][0],f2[cnt][0]);
   		  //sco+=1.00/(1.00+(dist/d02));
		  cnt++;
		 }
   		 sco+=1.00/(1.00+(dist/d02));
   		}
  	 //printf("Len= %d D0= %.3f rmsd= %.3f sco= %.3f cnt= %d\n",n,d0,rmsd,sco,cnt);
	 if(sco <= pre_sco || cnt < 4)
	  break;
	 calculate_rotation_rmsd(f1,f2,cnt,tmp_m_c,tmp_m_r,tmp_rmtx,&rmsd);
	 rotate_shift_cd(cd2,xyz2,n,tmp_m_c,tmp_m_r,tmp_rmtx);
	}
   if(max_sco<sco){
    max_sco=sco;
    //printf("max_sco= %f\n",max_sco);
     mov_com[0]=tmp_m_c[0];
     mov_com[1]=tmp_m_c[1];
     mov_com[2]=tmp_m_c[2];

     mov_to_ref[0]=tmp_m_r[0];
     mov_to_ref[1]=tmp_m_r[1];
     mov_to_ref[2]=tmp_m_r[2];
/*
    for(i1=0;i1<3;i1++){
     //mov_com[i1]=tmp_m_c[i1];
     //mov_to_ref[i1]=tmp_m_r[i1];

     //for(i2=0;i2<3;i2++){
//	rmtx[i1][i2]=tmp_rmtx[i1][i2];
//     }

	rmtx[i1][0]=tmp_rmtx[i1][0];
	rmtx[i1][1]=tmp_rmtx[i1][1];
	rmtx[i1][2]=tmp_rmtx[i1][2];
	
    }*/
	//fast ver
	rmtx[0][0]=tmp_rmtx[0][0];
	rmtx[0][1]=tmp_rmtx[0][1];
	rmtx[0][2]=tmp_rmtx[0][2];
	rmtx[1][0]=tmp_rmtx[1][0];
	rmtx[1][1]=tmp_rmtx[1][1];
	rmtx[1][2]=tmp_rmtx[1][2];
	rmtx[2][0]=tmp_rmtx[2][0];
	rmtx[2][1]=tmp_rmtx[2][1];
	rmtx[2][2]=tmp_rmtx[2][2];


   }
  }
 }
 //printf("TMscore= %f\n",max_sco);
 return max_sco;
}


double tm(int *in1,int *in2,
	  int *out1,int *out2,
	  double **cd1,int n1,double **cd2,int n2,
	  double Gopen,double Gext,double *smtx){

 double sco=0;
 int i,i1,i2;
 double l,rmsd,rmsd2;
 //DP
 int gal[RES],ng,gsco;
 int A1[RES],A2[RES];
 //int a1[RES],a2[RES];
 double tmp_cd1[RES][3],tmp_cd2[RES][3],xyz1[RES][3],xyz2[RES][3];
 
 //set d0
 int Lmin=n1;
 if(Lmin>n2)
  Lmin=n2;
 double d0=1.24*pow(Lmin-15,1.000/3.000)-1.8;

 int cnt=0;
 double mov_com[3],mov_to_ref[3],rmtx[3][3];
 double shift[3],pre_sco,dpsco;

 for(i=0;i<n1;i++){
  tmp_cd1[i][0]=cd1[i][0];
  tmp_cd1[i][1]=cd1[i][1];
  tmp_cd1[i][2]=cd1[i][2];
 }
 for(i=0;i<n2;i++){
  tmp_cd2[i][0]=cd2[i][0];
  tmp_cd2[i][1]=cd2[i][1];
  tmp_cd2[i][2]=cd2[i][2];
 }

 CopyAli(in1,in2,A1,A2,n1,n2);
 pre_sco=0;
 while(1){
	cnt=0;
 	for(i=0;i<n1;i++){
 	 if(A1[i]==-1)
 	  continue;

 	 xyz1[cnt][0]=cd1[i][0];
   	 xyz1[cnt][1]=cd1[i][1];
   	 xyz1[cnt][2]=cd1[i][2];

   	 xyz2[cnt][0]=cd2[A1[i]][0];
   	 xyz2[cnt][1]=cd2[A1[i]][1];
   	 xyz2[cnt][2]=cd2[A1[i]][2];
   	 //printf("*%d:%d %.3f %.3f %.3f|",i,in1[i], xyz1[cnt][0],xyz1[cnt][1],xyz1[cnt][2]);
	 //printf("%.3f %.3f %.3f\n",xyz2[cnt][0],xyz2[cnt][1],xyz2[cnt][2]);
   	 cnt++;
 	}
  //printf("cnt=%d\n",cnt);
  sco=TMscore(xyz1,xyz2,cnt,d0,mov_com,mov_to_ref,rmtx);
  //rotate all cd2
  //printf("%f %f %f\n",mov_to_ref[0],mov_to_ref[1],mov_to_ref[2]);
  rotate_shift_cd(tmp_cd2,xyz2,n2,mov_com,mov_to_ref,rmtx);
  //set Smtx
  SetSmtxTm(tmp_cd1,n1,xyz2,n2,d0,smtx);
  //DP
  //printf("sco= %f pre_sco= %f %d %d\n",sco,pre_sco,n1,n2);
  dpsco=dp(smtx,Gopen,Gext,A1,n1,A2,n2,gal,&ng);
  //ShowAli(A1,n1);
  //printf("sco= %f pre_sco= %f dpsco= %f N= %d\n",sco,pre_sco,dpsco,ng);
  if(pre_sco >=sco)
   break;
  CopyAli(A1,A2,out1,out2,n1,n2);
  pre_sco=sco;
 }
 //CopyAli(A1,A2,out1,out2,n1,n2);
 return pre_sco;
}


void rotate_shift_cd(double in[][3],double out[][3],int n,double mov_com[3],double mov_to_ref[3],double rmtx[3][3]){
 int i,j,k;
 double tmp[3];
 for(i=0;i<n;i++){
  out[i][0]=out[i][1]=out[i][2]=0;

/* faster!
  for (j=0; j<3; j++)
   tmp[j]=in[i][j]-mov_com[j];
*/
   tmp[0]=in[i][0]-mov_com[0];
   tmp[1]=in[i][1]-mov_com[1];
   tmp[2]=in[i][2]-mov_com[2];
  

  //for (j=0; j<3; j++){
    //for (k=0; k<3; k++){
    // out[i][j] += rmtx[j][k] * tmp[k];
    //}
     out[i][0] += rmtx[0][0] * tmp[0];
     out[i][0] += rmtx[0][1] * tmp[1];
     out[i][0] += rmtx[0][2] * tmp[2];
     out[i][1] += rmtx[1][0] * tmp[0];
     out[i][1] += rmtx[1][1] * tmp[1];
     out[i][1] += rmtx[1][2] * tmp[2];
     out[i][2] += rmtx[2][0] * tmp[0];
     out[i][2] += rmtx[2][1] * tmp[1];
     out[i][2] += rmtx[2][2] * tmp[2];
  //}
  //for (j=0; j<3; j++)
  // out[i][j] +=mov_com[j]+mov_to_ref[j];
   out[i][0] +=mov_com[0]+mov_to_ref[0];
   out[i][1] +=mov_com[1]+mov_to_ref[1];
   out[i][2] +=mov_com[2]+mov_to_ref[2];
 }
}

void SetSmtxTm(double cd1[][3],int n1,double cd2[][3],int n2,double d0,double *smtx){
 int i,j,k;
 double dist,d02=d0*d0;
 for(i=0;i<n1;i++){
  k=i*n2;
 	for(j=0;j<n2;j++){
	 dist=L2(cd1[i][0],cd1[i][1],cd1[i][2],cd2[j][0],cd2[j][1],cd2[j][2]);
	  smtx[j+k]=1.00/(1.00+(dist)/(d02));
	  //smtx[j+i*n2]=0.00;
	 //printf("smtx %d %d %f %f %f\n",i,j,dist,cd1[i][0],cd2[j][0]);
	}
 }
}

double GapLessAlign(int *out1,int *out2,
	  double **cd1,int n1,double **cd2,int n2, int L){
 int i,j,cnt=0;
 int k,Lmin=n1;
 int max_i=0,max_j=0;
 double xyz1[RES][3],xyz2[RES][3];

 if(n1>n2)
  Lmin=n2;
 //Lmin=L;

 double d0=1.24*pow(Lmin-15,1.000/3.000)-1.8;
 double mov_com[3],mov_to_ref[3],rmtx[3][3];

 double sco,max=0;
 double rms;
 //position
 for(i=0;i+L-1<n1;i++){
  //copy coords
  for(k=0;k<L;k++){
   xyz1[k][0]=cd1[i+k][0];
   xyz1[k][1]=cd1[i+k][1];
   xyz1[k][2]=cd1[i+k][2];
  }
	for(j=0;j+L-1<n2;j++){
	 //copy coords
	 for(k=0;k<L;k++){
	  xyz2[k][0]=cd2[j+k][0];
   	  xyz2[k][1]=cd2[j+k][1];
   	  xyz2[k][2]=cd2[j+k][2];
	 }
	 fast_rmsd(xyz1,xyz2,L,&rms);
	 if(rms>10)
	  continue;
	 sco=TMscore(xyz1,xyz2,L,d0,mov_com,mov_to_ref,rmtx);
	 //printf("%d %d (%d) sco= %f rms= %f\n",i,j,L,sco,rms);
	 if(max<sco){
	  max=sco;
	  max_i=i;max_j=j;
	 }
	 cnt++;
	}
 }
 printf("#Gapless= %d max= %f i= %d j= %d\n",cnt,max,max_i,max_j);
 //input alignment
 //init
 for(i=0;i<n1;i++) out1[i]=-1;
 for(i=0;i<n2;i++) out2[i]=-1;

 for(i=0;i<L;i++){
  out1[max_i+i]=max_j+i;
  out2[max_j+i]=max_i+i;
 }

 return max;
}
