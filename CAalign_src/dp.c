//Dynamic Programming
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "struct.h"
#include "func.h"
#include "dali.h"
#define GRID(N,a,b,c) c+N*(b+N*a)
#define GRID3D(N1,N2,a,b,c) c+N2*(b+N1*a)
#define GRID2D(N1,N2,a,b) (b)+(N2)*(a)

/*
//Smith & Waterman Local Alignment
//pointer
#define NON -1
#define UP 0
#define LEFT 1
#define DIA 2
#define GAP -1

typedef struct{ 
	double sco; 
	int poi;
} DPMTX;
*/
//Semi-global alignment
//Gaps of N and C-terminal are not penalized
/*
!!No gap penalty
AAAAAAAAAAAAAA
---BBBBBBB----

*/

//malloced in main
extern DPMTX *dmtx;

double dp(double *Smtx,double GapOpen,double GapExt,int *al1,int n1,int *al2,int n2,int *gal,int *glen){
 //puts("#start DP");
 //DPMTX dmtx[RES*RES];
 //DPMTX dmtx[n1*n2];
 int i1,i2;
 int N1,N2;
 double dia_sco,up_sco,left_sco;
 int gridid;
 int id2d,m;
 //Warning!!
 /*
 0:0 is not corresponding with Smtx[GRID2D(n1,n2,0,0)]
 0:x and x:1 is gap data of N-terminal
 use GRID2D(N1,N2,a+1,b+1)
 */


 N1=n1+1;
 N2=n2+1;

 dmtx[0].sco=0;
 dmtx[0].poi=NON;
 //puts("init");
 //for(i1=1;i1<=n1;i1++)
 for(i1=0;i1<n1;i1++)
  al1[i1]=-1;
 //for(i1=1;i1<=n2;i1++)
 for(i1=0;i1<n2;i1++)
  al2[i1]=-1;
 
 //init N-terminal DPMTX data SW
 for(i1=1;i1<=n1;i1++){
  dmtx[GRID2D(N1,N2,i1,0)].sco=0;
  dmtx[GRID2D(N1,N2,i1,0)].poi=NON;
 }
 for(i2=1;i2<=n2;i2++){
  dmtx[GRID2D(N1,N2,0,i2)].sco=0;
  dmtx[GRID2D(N1,N2,0,i2)].poi=NON;
 }
 //fill
 for(i1=1;i1<=n1;i1++){
  m=N2*i1;
 	for(i2=1;i2<=n2;i2++){
/*
	 if(i1==i2){
	 printf("%d-%d SCO=%.2f\n",i1,i2,Smtx[GRID2D(n1,n2,i1-1,i2-1)]);
	 printf("%d*%d %d %d =%d\n",n1,n2,i1-1,i2-1,GRID2D(n1,n2,i1-1,i2-1));
	 }
*/
	 //diagonal score
	 dia_sco=dmtx[GRID2D(N1,N2,i1-1,i2-1)].sco+Smtx[GRID2D(n1,n2,i1-1,i2-1)];
	 //gap score
	 //up score
	 up_sco=dmtx[GRID2D(N1,N2,i1-1,i2)].sco;
	 //gap open? or extention?
	 //open
	 if(dmtx[GRID2D(N1,N2,i1-1,i2)].poi==DIA){
	  up_sco-=GapOpen;
	 }else{
	  up_sco-=GapExt;
	 }
	 //left score
	 left_sco=dmtx[GRID2D(N1,N2,i1,i2-1)].sco;
	 //gap open? or extention?
	 //open
	 if(dmtx[GRID2D(N1,N2,i1,i2-1)].poi==DIA){
	  left_sco-=GapOpen;
	 }else{
	  left_sco-=GapExt;
	 }
	 //choose max score
	 //id2d=GRID2D(N1,N2,i1,i2);
	 id2d=i2+m;

	 if(dia_sco<=0 && up_sco <=0 && left_sco <=0){//NEW SW method
	  dmtx[id2d].poi=NON;
	  dmtx[id2d].sco=0;
	 }else if(dia_sco>=up_sco){
		if(dia_sco>=left_sco){
		 dmtx[id2d].poi=DIA;
		 dmtx[id2d].sco=dia_sco;
		}else{
		 dmtx[id2d].poi=LEFT;
		 dmtx[id2d].sco=left_sco;
		}
	 }else{
		if(up_sco>=left_sco){
		 dmtx[id2d].poi=UP;
		 dmtx[id2d].sco=up_sco;
		}else{
		 dmtx[id2d].poi=LEFT;
		 dmtx[id2d].sco=left_sco;
		}
	 }

	 //printf("%d:%d dia:%.3f up%.3f left:%.3f sco:%.3f poi %d\n",i1,i2,dia_sco,up_sco,left_sco,dmtx[GRID2D(N1,N2,i1,i2)].sco,dmtx[GRID2D(N1,N2,i1,i2)].poi);

	}
 }
 //puts("Trace back");
 //trace back
 //start from highest score bottom or end of columms Smith-Waterman
 //end zero gaps
 double highest=0;
 int i=0,j=0;
/*
 for(i1=0;i1<=n1;i1++){
  if(highest < dmtx[GRID2D(N1,N2,i1,n2)].sco){
   highest=dmtx[GRID2D(N1,N2,i1,n2)].sco;
   i=i1;
   j=n2;
  }
 }
 for(i2=0;i2<=n2;i2++){
  if(highest < dmtx[GRID2D(N1,N2,n1,i2)].sco){
   highest=dmtx[GRID2D(N1,N2,n1,i2)].sco;
   i=n1;
   j=i2;
  }
 }
*/

 //New Smith-Waterman
 for(i1=1;i1<=n1;i1++){
  m=N2*i1;
  for(i2=1;i2<=n2;i2++){
   id2d=i2+m;
   //if(highest < dmtx[GRID2D(N1,N2,i1,i2)].sco){
   if(highest < dmtx[id2d].sco){
   //if(highest <= dmtx[id2d].sco){ //New 2014.5.22!!!<-More Coverage, less quality
    //highest=dmtx[GRID2D(N1,N2,i1,i2)].sco;
    highest=dmtx[id2d].sco;
    i=i1;
    j=i2;
   }
  }
 }

 //printf("#START FROM %d %d sco=%.3f\n",i,j,highest);
 //printf("START FROM %d %d sco=%.3f\n",i+1,j+1,dmtx[GRID2D(N1,N2,i+1,j+1)].sco);
 //printf("START FROM %d %d smtxsco=%.3f\n",i+1,j+1,Smtx[GRID2D(n1,n2,i-1,j-1)]);
 //puts("#TRACE BACK");
 int gpos=0;//reverse position for global alignment
 //Becareful!! i&j is start from 1 -> fixed start from 0
 while(1){
 	gridid=GRID2D(N1,N2,i,j);
	//printf("%d-%d %.3f\n",i,j,dmtx[GRID2D(N1,N2,i,j)].sco);
	//if(dmtx[gridid].poi==NON){//end when 0:x or x:0
	 //al1[i-1]=GAP; al2[j-1]=GAP;
	 //printf("DPEND in %d %d\n",i-1,j-1);
	// break;
	//}
	if(dmtx[gridid].poi==DIA||dmtx[gridid].poi==NON){
	 //al1[i]=j; al2[j]=i;
	 al1[i-1]=j-1; al2[j-1]=i-1;
	 //gal[gpos*2]=i; gal[gpos*2+1]=j; gpos++;
	 i--;j--;
	}else if(dmtx[gridid].poi==LEFT){
	 //al1[i]=j; al2[j]=GAP;
	 al1[i-1]=j-1; al2[j-1]=GAP;
	 //gal[gpos*2]=GAP; gal[gpos*2+1]=j; gpos++;
	 j--;
	}else if(dmtx[gridid].poi==UP){
	 //al1[i]=-1; al2[j]=i;
	 al1[i-1]=-1; al2[j-1]=i-1;
	 //gal[gpos*2]=i; gal[gpos*2+1]=GAP; gpos++;
	 i--;
	}
	if(dmtx[gridid].poi==NON){
	 //printf("DPEND in %d %d\n",i-1,j-1);
	 break;
	}

 }
 *glen=gpos;
 return highest;
}

double dp_afp(double *Smtx,double GapOpen,double GapExt,int *al1,int n1,int *al2,int n2,int N,int *gal,int *glen){
 int i1,i2;
 int N1,N2;
 double dia_sco,up_sco,left_sco;
 int gridid;
 int id2d;
 //Warning!!
 /*
 0:0 is not corresponding with Smtx[GRID2D(n1,n2,0,0)]
 0:x and x:1 is gap data of N-terminal
 use GRID2D(N1,N2,a+1,b+1)
 */

 //Best One AFP = 10.00

 N1=n1+1;
 N2=n2+1;

 dmtx[0].sco=0;
 dmtx[0].poi=NON;
 double max=0;
 //puts("init");
 for(i1=0;i1<n1;i1++)
  al1[i1]=-1;
 for(i1=0;i1<n2;i1++)
  al2[i1]=-1;
 
 //init N-terminal DPMTX data
 for(i1=1;i1<=n1-N+1;i1++){
  dmtx[GRID2D(N1,N2,i1,0)].sco=0;
  dmtx[GRID2D(N1,N2,i1,0)].poi=NON;
 }
 for(i2=1;i2<=n2-N+1;i2++){
  dmtx[GRID2D(N1,N2,0,i2)].sco=0;
  dmtx[GRID2D(N1,N2,0,i2)].poi=NON;
 }
 //fill N-termilan+N on al1
 for(i1=1;i1<=n1-N+1;i1++){
  for(i2=1;i2<=N;i2++){
   dmtx[GRID2D(N1,N2,i1,i2)].sco=Smtx[GRID2D(n1,n2,i1-1,i2-1)];
   dmtx[GRID2D(N1,N2,i1,i2)].poi=NON;
   //if(max<dmtx[GRID2D(N1,N2,i1,i2)].sco) max=dmtx[GRID2D(N1,N2,i1,i2)].sco;
  }
 }
 //fill N-termilan+N on al2
 for(i2=N+1;i2<=n2-N+1;i2++){
  for(i1=1;i1<=N;i1++){
   dmtx[GRID2D(N1,N2,i1,i2)].sco=Smtx[GRID2D(n1,n2,i1-1,i2-1)];
   dmtx[GRID2D(N1,N2,i1,i2)].poi=NON;
   //if(max<dmtx[GRID2D(N1,N2,i1,i2)].sco) max=dmtx[GRID2D(N1,N2,i1,i2)].sco;
  }
 }

 double est_later;
 //fill
 for(i1=N+1;i1<=n1-N+1;i1++){
 	for(i2=N+1;i2<=n2-N+1;i2++){
/*
	 if(i1==i2){
	 printf("%d-%d SCO=%.2f\n",i1,i2,Smtx[GRID2D(n1,n2,i1-1,i2-1)]);
	 printf("%d*%d %d %d =%d\n",n1,n2,i1-1,i2-1,GRID2D(n1,n2,i1-1,i2-1));
	 }
*/

/*
	//New---------------
	 if(dmtx[GRID2D(N1,N2,i1-1,i2)].poi==dmtx[GRID2D(N1,N2,i1,i2-1)].poi==dmtx[GRID2D(N1,N2,i1-1,i2-1)].poi==STOP){
	  dmtx[id2d].poi=STOP;
	  dmtx[id2d].sco=0;
	  continue;
	 }
	//END new-----------
*/

	 //diagonal score
	 dia_sco=dmtx[GRID2D(N1,N2,i1-N,i2-N)].sco+Smtx[GRID2D(n1,n2,i1-1,i2-1)];
	 //gap score
	 //up score
	 up_sco=dmtx[GRID2D(N1,N2,i1-1,i2)].sco;
	 //gap open? or extention?
	 //open
	 if(dmtx[GRID2D(N1,N2,i1-1,i2)].poi==DIA){
	  up_sco-=GapOpen;
	 }else{
	  up_sco-=GapExt;
	 }
	 //left score
	 left_sco=dmtx[GRID2D(N1,N2,i1,i2-1)].sco;
	 //gap open? or extention?
	 //open
	 if(dmtx[GRID2D(N1,N2,i1,i2-1)].poi==DIA){
	  left_sco-=GapOpen;
	 }else{
	  left_sco-=GapExt;
	 }
/*
	 //New---------------
	 est_later=((n1-i1)-(n2-i2))*10.00;
	 if(est_later<0)
	  est_later*=-1.00;

	 if(est_later+dia_sco<max && est_later+up_sco<max && est_later+left_sco<max){
	  dmtx[id2d].poi=STOP;
	  dmtx[id2d].sco=0;
	  continue;
	 }
	 if(max<dia_sco) max=dia_sco;
	 if(max<left_sco) max=left_sco;
	 if(max<up_sco) max=up_sco;
	 //END New----------
*/
	 //choose max score
	 id2d=GRID2D(N1,N2,i1,i2);
	 if(dia_sco<=0 && up_sco <=0 && left_sco <=0){//NEW SW method
	  dmtx[id2d].poi=NON;
	  dmtx[id2d].sco=0;
	 }else if(dia_sco>=up_sco){
		if(dia_sco>=left_sco){
		 dmtx[id2d].poi=DIA;
		 dmtx[id2d].sco=dia_sco;
		}else{
		 dmtx[id2d].poi=LEFT;
		 dmtx[id2d].sco=left_sco;
		}
	 }else{
		if(up_sco>=left_sco){
		 dmtx[id2d].poi=UP;
		 dmtx[id2d].sco=up_sco;
		}else{
		 dmtx[id2d].poi=LEFT;
		 dmtx[id2d].sco=left_sco;
		}
	 }

	 //printf("%d:%d dia:%.3f up%.3f left:%.3f sco:%.3f poi %d\n",i1,i2,dia_sco,up_sco,left_sco,dmtx[GRID2D(N1,N2,i1,i2)].sco,dmtx[GRID2D(N1,N2,i1,i2)].poi);

	}
 }
 //Trace Back---------------
 //puts("Trace back");
 //trace back
 //start from highest score bottom or end of columms
 //end zero gaps
 double highest=0;
 int i=0,j=0;

 //New Smith-Waterman
 for(i1=1;i1<=n1-N+1;i1++){
  for(i2=1;i2<=n2-N+1;i2++){
   if(highest < dmtx[GRID2D(N1,N2,i1,i2)].sco){
    highest=dmtx[GRID2D(N1,N2,i1,i2)].sco;
    i=i1;
    j=i2;
   }
  }
 }

 //printf("#START FROM %d %d sco=%.3f\n",i,j,highest);
 //printf("START FROM %d %d sco=%.3f\n",n1,n2,dmtx[GRID2D(N1,N2,n1-1,n2-1)].sco);
 //puts("#TRACE BACK");
 int gpos=0;//reverse position for global alignment
 //Becareful!! i&j is start from 1 -> fixed start from 0
 while(1){
 	gridid=GRID2D(N1,N2,i,j);
	//printf("%d-%d %.3f poi= %d\n",i,j,dmtx[GRID2D(N1,N2,i,j)].sco,dmtx[gridid].poi);
/*
	if(dmtx[gridid].poi==NON){//end when 0:x or x:0
	 printf("DPEND in %d %d\n",i-1,j-1);
	 break;
	}
*/
	if(dmtx[gridid].poi==DIA||dmtx[gridid].poi==NON){
	 	//AFP
		for(i1=0;i1<N;i1++){
		 al1[i-1+i1]=j-1+i1; 
		 al2[j-1+i1]=i-1+i1;
		 //gal[gpos*2]=i+i1; gal[gpos*2+1]=j+i1;
		 gpos++;
		}
	 //back Fragment length
	 i-=N;j-=N;
/*
	 al1[i-1]=j-1; al2[j-1]=i-1;
	 gal[gpos*2]=i; gal[gpos*2+1]=j; 
	 gpos++;
	 i--;j--;
*/
	}else if(dmtx[gridid].poi==LEFT){
	 al1[i-1]=j-1; al2[j-1]=GAP;
	 //gal[gpos*2]=GAP; gal[gpos*2+1]=j; gpos++;
	 j--;
	}else if(dmtx[gridid].poi==UP){
	 al1[i-1]=-1; al2[j-1]=i-1;
	 //gal[gpos*2]=i; gal[gpos*2+1]=GAP; gpos++;
	 i--;
	}
	if(dmtx[gridid].poi==NON){//end when 0:x or x:0
	 //printf("DPEND in %d %d\n",i-1,j-1);
	 break;
	}
 }
 *glen=gpos;
 return highest;
}



double dp_afp_fwd(double *Smtx,double GapOpen,double GapExt,
		int *al1,int n1,int *al2,int n2,int N,
		int *gal,int *glen,
		DMTX *d1,DMTX *d2,double w){
 int i1,i2;
 int N1,N2;
 double dia_sco,up_sco,left_sco;
 int gridid;
 int id2d;
 //Warning!!
 /*
 0:0 is not corresponding with Smtx[GRID2D(n1,n2,0,0)]
 0:x and x:1 is gap data of N-terminal
 use GRID2D(N1,N2,a+1,b+1)
 */


 N1=n1+1;
 N2=n2+1;

 dmtx[0].sco=0;
 dmtx[0].poi=NON;
 //puts("init");
 for(i1=0;i1<n1;i1++)
  al1[i1]=-1;
 for(i1=0;i1<n2;i1++)
  al2[i1]=-1;
 
 //init N-terminal DPMTX data
 for(i1=1;i1<=n1-N+1;i1++){
  dmtx[GRID2D(N1,N2,i1,0)].sco=0;
  dmtx[GRID2D(N1,N2,i1,0)].poi=NON;
 }
 for(i2=1;i2<=n2-N+1;i2++){
  dmtx[GRID2D(N1,N2,0,i2)].sco=0;
  dmtx[GRID2D(N1,N2,0,i2)].poi=NON;
 }
 //fill N-termilan+N on al1
 for(i1=1;i1<=n1-N+1;i1++){
  for(i2=1;i2<=N;i2++){
   dmtx[GRID2D(N1,N2,i1,i2)].sco=Smtx[GRID2D(n1,n2,i1-1,i2-1)];
   dmtx[GRID2D(N1,N2,i1,i2)].poi=NON;
  }
 }
 //fill N-termilan+N on al2
 for(i2=N+1;i2<=n2-N+1;i2++){
  for(i1=1;i1<=N;i1++){
   dmtx[GRID2D(N1,N2,i1,i2)].sco=Smtx[GRID2D(n1,n2,i1-1,i2-1)];
   dmtx[GRID2D(N1,N2,i1,i2)].poi=NON;
  }
 }

 double Scon=0;
 //fill
 for(i1=N+1;i1<=n1-N+1;i1++){
 	for(i2=N+1;i2<=n2-N+1;i2++){
	 //diagonal score
	 dia_sco=dmtx[GRID2D(N1,N2,i1-N,i2-N)].sco+Smtx[GRID2D(n1,n2,i1-1,i2-1)];
	 Scon=TraceBack(d1,d2,N1,N2,N,i1,i2);//i1&i2 start from 1 (dpmtx base)
	 //printf("dia= %f Scom= %f\n",dia_sco,Scon*w);
	 dia_sco+=w*Scon;
	 //gap score
	 //up score
	 up_sco=dmtx[GRID2D(N1,N2,i1-1,i2)].sco;
	 //gap open? or extention?
	 //open
	 if(dmtx[GRID2D(N1,N2,i1-1,i2)].poi==DIA){
	  up_sco-=GapOpen;
	 }else{
	  up_sco-=GapExt;
	 }
	 //left score
	 left_sco=dmtx[GRID2D(N1,N2,i1,i2-1)].sco;
	 //gap open? or extention?
	 //open
	 if(dmtx[GRID2D(N1,N2,i1,i2-1)].poi==DIA){
	  left_sco-=GapOpen;
	 }else{
	  left_sco-=GapExt;
	 }
	 //choose max score
	 id2d=GRID2D(N1,N2,i1,i2);
	 if(dia_sco<=0 && up_sco <=0 && left_sco <=0){//NEW SW method
	  dmtx[id2d].poi=NON;
	  dmtx[id2d].sco=0;
	 }else if(dia_sco>=up_sco){
		if(dia_sco>=left_sco){
		 dmtx[id2d].poi=DIA;
		 dmtx[id2d].sco=dia_sco;
		}else{
		 dmtx[id2d].poi=LEFT;
		 dmtx[id2d].sco=left_sco;
		}
	 }else{
		if(up_sco>=left_sco){
		 dmtx[id2d].poi=UP;
		 dmtx[id2d].sco=up_sco;
		}else{
		 dmtx[id2d].poi=LEFT;
		 dmtx[id2d].sco=left_sco;
		}
	 }

	 //printf("%d:%d dia:%.3f up%.3f left:%.3f sco:%.3f poi %d\n",i1,i2,dia_sco,up_sco,left_sco,dmtx[GRID2D(N1,N2,i1,i2)].sco,dmtx[GRID2D(N1,N2,i1,i2)].poi);

	}
 }
 //Trace Back---------------
 //puts("Trace back");
 //trace back
 //start from highest score bottom or end of columms
 //end zero gaps
 double highest=0;
 int i=0,j=0;

 //New Smith-Waterman
 for(i1=1;i1<=n1-N+1;i1++){
  for(i2=1;i2<=n2-N+1;i2++){
   if(highest < dmtx[GRID2D(N1,N2,i1,i2)].sco){
    highest=dmtx[GRID2D(N1,N2,i1,i2)].sco;
    i=i1;
    j=i2;
   }
  }
 }

 //printf("#START FROM %d %d sco=%.3f\n",i,j,highest);
 //printf("START FROM %d %d sco=%.3f\n",n1,n2,dmtx[GRID2D(N1,N2,n1-1,n2-1)].sco);
 //puts("#TRACE BACK");
 int gpos=0;//reverse position for global alignment
 //Becareful!! i&j is start from 1 -> fixed start from 0
 while(1){
 	gridid=GRID2D(N1,N2,i,j);
	//printf("%d-%d %.3f poi= %d\n",i,j,dmtx[GRID2D(N1,N2,i,j)].sco,dmtx[gridid].poi);
	if(dmtx[gridid].poi==DIA||dmtx[gridid].poi==NON){
	 	//AFP
		for(i1=0;i1<N;i1++){
		 al1[i-1+i1]=j-1+i1; 
		 al2[j-1+i1]=i-1+i1;
		 gal[gpos*2]=i+i1; gal[gpos*2+1]=j+i1;
		 gpos++;
		}
	 //back Fragment length
	 i-=N;j-=N;
	}else if(dmtx[gridid].poi==LEFT){
	 al1[i-1]=j-1; al2[j-1]=GAP;
	 gal[gpos*2]=GAP; gal[gpos*2+1]=j; gpos++;
	 j--;
	}else if(dmtx[gridid].poi==UP){
	 al1[i-1]=-1; al2[j-1]=i-1;
	 gal[gpos*2]=i; gal[gpos*2+1]=GAP; gpos++;
	 i--;
	}
	if(dmtx[gridid].poi==NON){//end when 0:x or x:0
	 //printf("DPEND in %d %d\n",i-1,j-1);
	 break;
	}
 }
 *glen=gpos;
 return highest;
}

double TraceBack(DMTX *d1,DMTX *d2,int N1,int N2,int N,int I,int J){
 int gridid;
 int i1,i2,p11,p12,p21,p22;
 int i=I-N;//DPMTX based, start from 1
 int j=J-N;
 double sco=0;
 //printf("TraceBack %d %d\n",I,J);
 double bnd=0,t=0,a1,a2;
 while(1){
 	gridid=GRID2D(N1,N2,i,j);
	//printf("%d-%d %.3f poi= %d\n",i,j,dmtx[GRID2D(N1,N2,i,j)].sco,dmtx[gridid].poi);
	if(dmtx[gridid].poi==DIA||dmtx[gridid].poi==NON){
	 	//AFP, check dmtx score
		/*
		i+i1-1            I+i2-1
		p11++++++---------p12++++++
		p21++++++---------p22++++++
		j+i1-1            J+i2-1
		*/
		for(i1=0;i1<N;i1++){
		 p11=i+i1-1;
		 p21=j+i1-1;
		 for(i2=0;i2<N;i2++){
		  p12=I+i2-1;
		  p22=J+i2-1;

		  //printf("%d %d %d %d\n",p11,p12,p21,p22);


		  //a1=d1->area[0][0]+d1->area[0][0];
		  a1=d1->area[p11][p12]+d1->area[p12][p11];
		  a2=d2->area[p21][p22]+d2->area[p22][p21];
		  //printf("%d %d %d %d %f %f\n",p11,p12,p21,p22,a1,a2);

		  bnd+=cad_bnd(a1,a2);
		  bnd+=cad_bnd(a2,a1);
		  t+=a1+a2;

		 //al1[i-1+i1]=j-1+i1; 
		 //al2[j-1+i1]=i-1+i1;
		 //gal[gpos*2]=i+i1; gal[gpos*2+1]=j+i1;
		 //gpos++;
		 }
		}
	 //back Fragment length
	 i-=N;j-=N;
	}else if(dmtx[gridid].poi==LEFT){
	 //al1[i-1]=j-1; al2[j-1]=GAP;
	 //gal[gpos*2]=GAP; gal[gpos*2+1]=j; gpos++;
	 j--;
	}else if(dmtx[gridid].poi==UP){
	 //al1[i-1]=-1; al2[j-1]=i-1;
	 //gal[gpos*2]=i; gal[gpos*2+1]=GAP; gpos++;
	 i--;
	}
	if(dmtx[gridid].poi==NON){//end when 0:x or x:0
	 //printf("DPEND in %d %d\n",i-1,j-1);
	 break;
	}
 }
 sco=(t-bnd)/t;
 return sco;
}

