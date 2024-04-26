#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "struct.h"
#include "func.h"
#define RMS_CUTOFF 2.00 //cluster形成のrmsd閾値
#define MIN_CLS 1 //最低形成クラスター箇数

//data内のマトリクスからクラスター解析を行ないます。

int clustering(double **r,int n,double rms_cut){
 UINT **v,**p;
 //DATA *d;
 PDB *lig;
 //float **r;
 UINT i,j,l,m;
 int *noOfG;
 int *gno;//group number
 UINT no;
 UINT res_no;
 COORD *tmp1,*tmp2;
 int total_cls,cnt;
 UINT top_m,top_l;
 double top;
 float min=0;
 int flag,new_id;
 //no=30000;
 //no=d->NumOfData;
 res_no=lig->NumOfRes;

 //Memory get
 if((noOfG=(int *)malloc(n*sizeof(int)))==NULL)
  malloc_error("noOfG");

 if((gno=(int *)malloc(n*sizeof(int)))==NULL)
  malloc_error("gno");
 printf("##FIN Memory get stage in clustering##\n");
 printf("##START Clast Calc##\n");

 	//init group number
        for(l=0;l<n;l++){
                gno[l]=l;
                noOfG[l]=1;
        }
	total_cls=n;
 	printf("##Finish set distance mtx##\n");
        while(total_cls > MIN_CLS){
		//クラスター間の最小距離
                cnt=0;
		top=rms_cut+1;
                for(l=0;l<n;l++){
                 if(gno[l] != l)
		  continue;
               		for(m=l+1;m<n;m++){
                     	 if(gno[m] != m)
			  continue;
                         //if(top>r[l][m] || cnt==0){
                         if(top>r[l][m]){
                          top=r[l][m];
                          top_m=m;
                          top_l=l;
			  //return FALSE;
                          cnt++;
                         }
                        }
                }
		//printf("mindis=%.3f %d(%d %d)\n",top,cnt,top_l,top_m);
		if(top > rms_cut )
		 break;

		//renew dmtx
		for(l=0;l<n;l++){
		 if(gno[l] != l)
		  continue;
	  	 for(m=l+1;m<n;m++){
		  //if(gno[m] == m)????
		  if(gno[m] != m)
		   continue;
		  if(l!=top_m&&m!=top_m){ 
		   if(l==top_l){
		    //r[l][m]=next_d(r[top_l][m],r[top_m][m],r[top_l][top_m],noOfG[m],noOfG[top_l],noOfG[top_m]); 
		    r[l][m]=next_d(r[l][m],r[top_m][m],r[l][top_m],noOfG[m],noOfG[l],noOfG[top_m]); 
		    r[m][l]=r[l][m];
		   } else if(m==top_l){ 
		    //r[l][m]=next_d(r[l][top_l],r[l][top_m],r[top_l][top_m],noOfG[l],noOfG[top_l],noOfG[top_m]); 
		    r[l][m]=next_d(r[l][m],r[l][top_m],r[m][top_m],noOfG[l],noOfG[m],noOfG[top_m]); 
		    r[m][l]=r[l][m]; 
		   } 
		  } 
		 } 
		}
		//renew group number 
		 noOfG[top_l]+=noOfG[top_m];
		 noOfG[top_m]=0;
		 for(l=0;l<n;l++){
                        if(gno[l]==top_m)
                                gno[l]=top_l;
                 }
                --total_cls;
		//printf("mindis=%.3f %d\n",top,total_cls);
	}
	puts("#FIN clustering");
	//結果の表示
	//IDのふりなおし
	new_id=0;
	for(l=0;l<n;l++){
		if(noOfG[l] <1)
		 continue;
		
		printf("ClusterNo:%5d(%4d)|",l,noOfG[l]);
		flag=0;
		for(m=0;m<n;m++){
		 	if(gno[m]==l){
			 flag++;
		  	 printf("%6d",m);
			 //d->clust_no[m]=new_id;
			 if(flag == noOfG[l])
			  break;
			}
		}
		printf(" =>new id is %d\n",new_id);
		if(flag >0)
		 new_id++;
	}

}

float next_d(float d1,float d2,float d3,int n1,int n2,int n3){
        float ans;
	//new fast mode!!
	if(d1<0 || d2<0 || d3< 0)
	 return -1.00;
        //ans=d1*(n1+n3)/(n1+n2+n3)+d2*(n2+n3)/(n1+n2+n3)+d3*(n3)/(n1+n2+n3);
        //ans=X(d1)*(n1+n3)/(n1+n2+n3)+X(d2)*(n2+n3)/(n1+n2+n3)+X(d3)*(n3)/(n1+n2+n3);
        ans=(X(d1)*(double)(n1+n3)+X(d2)*(double)(n2+n3)-X(d3)*(double)(n3))/(double)(n1+n2+n3);
        //return(ans);
        return(sqrt(ans));
}

int clustering_ward(short int **d,int *cid,int *csize, int n,double rms_cut,double rate_cut){

 UINT i,total_cls;
 UINT cnt,l,m;
 UINT top_m,top_l;
 bool flag;
 //double top;
 short int top;
 rms_cut*=100;//100 times
 //init group number
 for(i=0;i<n;i++){
  cid[i]=i;
  csize[i]=1;
 }
 total_cls=n;

 while(total_cls > 1){
  flag=false;
  top=(short int)(rms_cut+1.00);
	//search min dmtx
	for(l=0;l<n;l++){
      	 if(cid[l] != l)
	  continue;
       		for(m=l+1;m<n;m++){
              	 if(cid[m] != m)
		  continue;
                 if(top>d[l][m] && d[l][m] >= 0){
                  top=d[l][m];
                  top_m=m;
                  top_l=l;
                  flag=true;
                 }
         	}
  	}
	if(flag==false){
	 printf("#No new cluster\n");
	 break;
	}
	//check the seize of new cluster!!
	if((double)(csize[top_l]+csize[top_m])/(double)n > rate_cut){
	 printf("#New cluster size reach --> %.3f > %.3f\n",(double)(csize[top_l]+csize[top_m])/(double)n,rate_cut);
	 break;
	}
	//printf("mindis=%.3f (%d %d) size= %d = %d+%d\n",(double)top/100,top_l,top_m,csize[top_l]+csize[top_m],csize[top_l],csize[top_m]);
	if(top > rms_cut )
	 break;

	//renew dmtx
	for(l=0;l<n;l++){
	 if(cid[l] != l)
	  continue;
	  	 for(m=l+1;m<n;m++){
		  if(cid[m] != m)
		   continue;
		  if(l!=top_m&&m!=top_m){ 
		   if(l==top_l){
		    d[l][m]=next_d(d[l][m],d[top_m][m],d[l][top_m],csize[m],csize[l],csize[top_m]); 
		    d[m][l]=d[l][m];
		   } else if(m==top_l){ 
		    d[l][m]=next_d(d[l][m],d[l][top_m],d[m][top_m],csize[l],csize[m],csize[top_m]); 
		    d[m][l]=d[l][m]; 
		   } 
		  } 
		 } 
		}
		//renew group number 
		 csize[top_l]+=csize[top_m];
		 csize[top_m]=0;
		 for(l=0;l<n;l++){
                        if(cid[l]==top_m)
                                cid[l]=top_l;
                 }
                --total_cls;
		//printf("mindis=%.3f %d\n",(double)top,total_cls);
	}
	puts("#FIN clustering");
 return total_cls;
}

//find representative model for each cluster
int find_rmodel(short int **d,int *cid, int *csize,int n,int Ncls,char **dblist){
 UINT i,j,k;
 UINT cnt,Ccnt=0;
 int *s;//distance within the same cluster
 int *member;//Member list of the cluster
 int minid;
 int min;

 //malloc zone!
 if((s=(int *)malloc(sizeof(int)*n))==NULL)
  malloc_error("s in find_rmodel");
 if((member=(int *)malloc(sizeof(int)*n))==NULL)
  malloc_error("member in find_rmodel");

 
 for(i=0;i<n;i++){
  cnt=0;
  if(i < cid[i])
   printf("%d:cid=%d\n",i,cid[i]);
  if(i!=cid[i])
   continue;
  	//Find member
	for(j=i;j<n;j++){
	 if(cid[j] != i)
	  continue;
	 //printf("CID %d TAG= %d %d/%d\n",i,j,cnt,csize[i]);
	 member[cnt]=j;
	 cnt++;
	 //no new member
	 if(cnt==csize[i])
	  break;
	}
	//init
	minid=0;
	min=0;
	//init
	for(j=0;j<cnt;j++)
	 s[j]=0;

	//find cluster center
	for(j=0;j<cnt;j++){	
	 for(k=0;k<cnt;k++){
	  if(j==k)
	   continue;
	  //printf("*j%d %d %d %d %d\n",j,member[j],member[k],d[member[j]][member[k]],s[j]);
	  s[j]+=d[member[j]][member[k]];
	  //if(s[j]< 0)
	  //printf("j%d %d %d %d %d\n",j,member[j],member[k],d[member[j]][member[k]],s[j]);
	 }
	}
	//find min
	for(j=0;j<cnt;j++){
	 if(j==0||s[j]<min){
	   minid=j;
	   min=s[j];
  //printf("#%d %d = %d %.3f\n",j,minid,min,(double)s[minid]/(cnt*100));
	 }
	}
  printf("CID= %d CSIZE= %d avg=%.3f MODEL= %s\n",Ccnt,csize[i],(double)min/(cnt*100),dblist[member[minid]],(double)s[minid]/(cnt*100));
  Ccnt++;
  if(Ccnt==Ncls)//No new cluster
   break;  
  
 }
}
