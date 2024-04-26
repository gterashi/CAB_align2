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
#define DIST(a,b,c,d,e,f) sqrt((a-b)*(a-b)+(c-d)*(c-d)+(e-f)*(e-f))
#define L(a,b,c,d,e,f) (sqrt((d-a)*(d-a)+(e-b)*(e-b)+(f-c)*(f-c)))

static void cross(double a[3], double b[3], double c[3])
{
  a[0] = b[1]*c[2] - b[2]*c[1];
  a[1] = b[2]*c[0] - b[0]*c[2];
  a[2] = b[0]*c[1] - b[1]*c[0];
}

double TMscore(double cd1[][3],double cd2[][3],int n,double d0,double shift[3],double rmtx[3][3]){
 double sco;
 double xyz1[RES][3],xyz2[RES][3];
 double f1[RES][3],f2[RES][3];

 int L[4]={n,n/2,n/4,4};
 int cnt,i,i1,i2;

 int p,len,pos;
 double rmsd,mov_com[3],mov_to_ref[3];
/*
 //for move cd2
 for(i=0;i<n;i++){
  xyz2[i][0]=cd2[i1][0];
  xyz2[i][1]=cd2[i1][1];
  xyz2[i][2]=cd2[i1][2];
 }
*/

 for(p=0;p<4;p++){
  len=L[p];
  //Shift window From N to C
  for(pos=0;pos+len-1<n;pos++){
   printf("L= %d pos= %d\n",len,pos);
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
 	calculate_rotation_rmsd(f1,f2,len,mov_com,mov_to_ref,rmtx,&rmsd);
	
 	//cal moved cd
 	rotate_shift_cd(cd2,xyz2,n,mov_com,mov_to_ref,rmtx);
	for(i=0;i<len;i++){
	 printf("%d %.3f %.3f %.3f|",i, cd1[pos+i][0],cd1[pos+i][1],cd1[pos+i][2]);
	 printf("%.3f %.3f %.3f = %.3f \n",cd2[pos+i][0],cd2[pos+i][1],cd2[pos+i][2],L(f1[pos+i][0],f1[pos+i][1],f1[pos+i][2],cd2[pos+i][0],cd2[pos+i][1],cd2[pos+i][2]));
	}
 	printf("Len= %d D0= %.3f rmsd= %.3f\n",len,d0,rmsd);

  }
 }

 return sco;
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
 int a1[RES],a2[RES];
 double xyz1[RES][3],xyz2[RES][3];
 
 //set d0
 int Lmin=n1;
 if(Lmin>n2)
  Lmin=n2;
 double d0=1.24*pow(Lmin-15,1.000/3.000)-1.8;

 int cnt=0;
 double mov_com[3],mov_to_ref[3],rmtx[3][3];
 double shift[3];
 for(i=0;i<n1;i++){
  if(in1[i]==-1)
   continue;
   xyz1[cnt][0]=cd1[i][0];
   xyz1[cnt][1]=cd1[i][1];
   xyz1[cnt][2]=cd1[i][2];

   xyz2[cnt][0]=cd2[in1[i]][0];
   xyz2[cnt][1]=cd2[in1[i]][1];
   xyz2[cnt][2]=cd2[in1[i]][2];
   //printf("*%d:%d %.3f %.3f %.3f|",i,in1[i], xyz1[cnt][0],xyz1[cnt][1],xyz1[cnt][2]);
	//printf("%.3f %.3f %.3f\n",xyz2[cnt][0],xyz2[cnt][1],xyz2[cnt][2]);
   cnt++;
 }
 TMscore(xyz1,xyz2,cnt,d0,shift,rmtx);

/*
 //Copy Alignment
 for(i=0;i<n1;i++)
  a1[i]=in1[i];

 //iterative
 int L[4]={n1,n1/2,n1/4,4};

 int p,len,pos;
 for(p=0;p<4;p++){
  len=L[p];
  //From N to C
  for(pos=0;pos+len-1<n1;pos++){
   printf("L= %d pos= %d\n",len,pos);
 //input CA xyz coords
  	cnt=0;
 	for(i=0;cnt<len;i++){
	 i1=pos+i;
 	 if(in1[i1]==-1)
 	  continue;
  	 xyz1[cnt][0]=cd1[i1][0];
  	 xyz1[cnt][1]=cd1[i1][1];
  	 xyz1[cnt][2]=cd1[i1][2];

  	 xyz2[cnt][0]=cd2[in1[i1]][0];
  	 xyz2[cnt][1]=cd2[in1[i1]][1];
  	 xyz2[cnt][2]=cd2[in1[i1]][2];

  	 cnt++;
 	}
 	//rotate with first alignment
 	calculate_rotation_rmsd(xyz1,xyz2,cnt,mov_com,mov_to_ref,rmtx,&rmsd);
	
 	//cal moved cd
 	rotate_shift_cd(cd2,xyz2,n2,mov_com,mov_to_ref,rmtx);
 	printf("Lmin= %d D0= %.3f rmsd= %.3f\n",Lmin,d0,rmsd);

  }
 }

 //set Smtx
 SetSmtxTm(cd1,n1,xyz2,n2,d0,smtx);
 //DP
 sco=dp(smtx,Gopen,Gext,A1,n1,A2,n2,gal,&ng);

*/
 return 0;
}


void rotate_shift_cd(double in[][3],double out[][3],int n,double mov_com[3],double mov_to_ref[3],double rmtx[3][3]){
 int i,j,k;
 double tmp[3];
 for(i=0;i<n;i++){
  out[i][0]=out[i][1]=out[i][2]=0;

  for (j=0; j<3; j++)
   tmp[j]=in[i][j]-mov_com[j];

  for (j=0; j<3; j++){
    for (k=0; k<3; k++){
     out[i][j] += rmtx[j][k] * tmp[k];
    }
  }
  for (j=0; j<3; j++)
   out[i][j] +=mov_com[j]+mov_to_ref[j];
 }
}

void SetSmtxTm(double cd1[][3],int n1,double cd2[][3],int n2,double d0,double *smtx){
 int i,j;
 double dist;
 for(i=0;i<n1;i++){
 	for(j=0;j<n2;j++){
	 dist=L(cd1[i][0],cd1[i][1],cd1[i][2],cd2[i][0],cd2[i][1],cd2[i][2]);
	 smtx[j+i*n2]=1.00/(1.00+(dist*dist)/(d0*d0));
	}
 }
}
