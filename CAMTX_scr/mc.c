#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "mc.h"
#include "struct.h"


//#define POS2ID(a,b,c,d,e,f) ( (a) + (b)*(d) + (c)*(d)*(e) )
#define POS2ID(a,b,c,d,e,f) ( (a)*(e)*(f) + (b)*(f) + (c))

//area of triangle
#define HERON(a,b,c) 0.25*sqrt(((a)*(a)+(b)*(b)+(c)*(c))*((a)*(a)+(b)*(b)+(c)*(c))-2*((a)*(a)*(a)*(a)+(b)*(b)*(b)*(b)+(c)*(c)*(c)*(c)))

double heron(double a,double b,double c){
 double tmp=sqrt(X(a*a+b*b+c*c)-2*(a*a*a*a+b*b*b*b+c*c*c*c));
 return (0.25*tmp);
}

//Set voxel size
//2^x size
int SetupVoxelSize(UINT *size,double *min,double *max,double step){
 UINT a;
 double b[3];
 double c;

 b[0]=max[0]-min[0]; 
 b[1]=max[1]-min[1]; 
 b[2]=max[2]-min[2];

 c=b[0];

 if(b[1]>c)
  c=b[1];
 if(b[2]>c)
  c=b[2];

 UINT tmp=0;
 UINT i;

 if(step <=0)
  return FALSE;

 tmp=(UINT)(c/step);

/*
 //set minimum 2^x size
 for(i=2;;i*=2){
  if(i>tmp)
   break;
 }
*/
 //printf("#2^x = %d\n",i);

 //size[0]=size[1]=size[2]=i;


 //reduce memory space
 size[0]=(UINT)(b[0]/step)+1;
 size[1]=(UINT)(b[1]/step)+1;
 size[2]=(UINT)(b[2]/step)+1;

 return TRUE;
}

int MallocVoxel(VOXEL *v,double *b,UINT *s,double step){

 UINT size=s[0]*s[1]*s[2];
 UINT size2=(s[0]+1)*(s[1]+1)*(s[2]+1);
 UINT i;
 //printf("#NumOfVox = %d\n",size);
 if(size <= 0)
  return FALSE;

 //step
 v->step=step;

 //set base data
 v->base[0]=b[0];
 v->base[1]=b[1];
 v->base[2]=b[2];

 //set size
 v->size[0]=s[0];
 v->size[1]=s[1];
 v->size[2]=s[2];

 //class
 if((v->cl=(int *)malloc(sizeof(int)*size))==NULL)
  return FALSE;

 //corner 01
 if((v->corner=(unsigned int *)malloc(sizeof(unsigned int)*size2))==NULL)
  return FALSE;
 //corner value calloc
 //if((v->vertex=(float *)malloc(sizeof(float)*size*8))==NULL)
 if((v->cv=(double *)calloc(size2,sizeof(double)))==NULL)
  return FALSE;

 return TRUE;
}

int ClearVox(VOXEL *v){
 UINT size=v->size[0]*v->size[1]*v->size[2];
 UINT size2=(v->size[0]+1)*(v->size[1]+1)*(v->size[2]+1);
 UINT i;
 for(i=0;i<size;i++) v->cl[i]=0;
 for(i=0;i<size2;i++) v->corner[i]=0;
 for(i=0;i<size2;i++) v->cv[i]=0;

}

int PdbOnVoxel(PDB *p, VOXEL *v,double C){
 //set Atom radius custamizable
 //0:N 1:C 2:O 3:S 4:P 5:H 6:unknown
 
 //from Atom coords
 int i;
 UINT a[3],b[3],vpos[3];
 double tmp1[3],tmp2[3],vcd[3];
 int i1,i2,i3;
 UINT id;
 double dist;
 
 for(i=0;i<p->NumOfAtom;i++){
  //printf("#%d(Type %d) %.3f %.3f %.3f\n",i,p->TypeAtomId[i],p->xyz[i][0],p->xyz[i][1],p->xyz[i][2]);

  //ignore hydrogen and unknown atom
  if(p->TypeAtomId[i] > 4)
   continue;

  //CHECK HERE FOR METABALL!!

  tmp1[0]=p->xyz[i][0] - AtomRadius[p->TypeAtomId[i]]*C;
  tmp1[1]=p->xyz[i][1] - AtomRadius[p->TypeAtomId[i]]*C;
  tmp1[2]=p->xyz[i][2] - AtomRadius[p->TypeAtomId[i]]*C;

  tmp2[0]=p->xyz[i][0] + AtomRadius[p->TypeAtomId[i]]*C;
  tmp2[1]=p->xyz[i][1] + AtomRadius[p->TypeAtomId[i]]*C;
  tmp2[2]=p->xyz[i][2] + AtomRadius[p->TypeAtomId[i]]*C;


	//summation for near corner
	//FindNearBase(b,p->xyz[i],v->base,v->step);
	FindNearBase(a,b,tmp1,tmp2,v->base,v->step);
	//printf("BASE= %d %d %d\n",a[0],a[1],a[2]);
	//printf("TOP = %d %d %d\n",b[0],b[1],b[2]);
	for(vpos[0]=a[0];vpos[0]<=b[0];vpos[0]++){
	 for(vpos[1]=a[1];vpos[1]<=b[1];vpos[1]++){
	  for(vpos[2]=a[2];vpos[2]<=b[2];vpos[2]++){
	   //corner id
	   //set id **Corner Size is (n+1)*(n+1)*(n+1) !!!
	   id=POS2ID(vpos[0],vpos[1],vpos[2],v->size[0]+1,v->size[1]+1,v->size[2]+1);
	   CornerCoords(vcd,vpos,v->base,v->step);
	   //distance
	   dist=XyzDist(p->xyz[i],vcd);
	   //sum on corner value
	   v->cv[id]+=MetaBallFunc(dist,AtomRadius[p->TypeAtomId[i]],C);
	   //printf("%d:r=%.2f dist=%.3f %.2f %.2f %.2f <=> %.2f %.2f %.2f = %f\n",id,AtomRadius[p->TypeAtomId[i]],dist,p->xyz[i][0],p->xyz[i][1],p->xyz[i][2],vcd[0],vcd[1],vcd[2],v->cv[id]);

	  }
	 }
	}

 }
 return TRUE;
}

int PdbResOnVoxel(PDB *p,int r, VOXEL *v,double C){
 //set Atom radius custamizable
 //0:N 1:C 2:O 3:S 4:P 5:H 6:unknown
 
 //from Atom coords
 int i;
 UINT a[3],b[3],vpos[3];
 double tmp1[3],tmp2[3],vcd[3];
 int i1,i2,i3;
 UINT id;
 double dist;
 
 //for(i=0;i<p->NumOfAtom;i++){
 for(i=p->ResOnAtom[r];i<p->ResOnAtom[r+1];i++){
  //if(p->AtomOnRes[i]!=r)
  // break;
  //printf("#%d(Type %d) %.3f %.3f %.3f\n",i,p->TypeAtomId[i],p->xyz[i][0],p->xyz[i][1],p->xyz[i][2]);

  //ignore hydrogen and unknown atom
  if(p->TypeAtomId[i] > 4)
   continue;

  //CHECK HERE FOR METABALL!!

  tmp1[0]=p->xyz[i][0] - AtomRadius[p->TypeAtomId[i]]*C;
  tmp1[1]=p->xyz[i][1] - AtomRadius[p->TypeAtomId[i]]*C;
  tmp1[2]=p->xyz[i][2] - AtomRadius[p->TypeAtomId[i]]*C;

  tmp2[0]=p->xyz[i][0] + AtomRadius[p->TypeAtomId[i]]*C;
  tmp2[1]=p->xyz[i][1] + AtomRadius[p->TypeAtomId[i]]*C;
  tmp2[2]=p->xyz[i][2] + AtomRadius[p->TypeAtomId[i]]*C;


	//summation for near corner
	//FindNearBase(b,p->xyz[i],v->base,v->step);
	FindNearBase(a,b,tmp1,tmp2,v->base,v->step);
	//printf("BASE= %d %d %d\n",a[0],a[1],a[2]);
	//printf("TOP = %d %d %d\n",b[0],b[1],b[2]);
	for(vpos[0]=a[0];vpos[0]<=b[0];vpos[0]++){
	 i1=(v->size[1]+1)*(v->size[2]+1)*vpos[0];
	 for(vpos[1]=a[1];vpos[1]<=b[1];vpos[1]++){
	  i2=i1+(v->size[2]+1)*vpos[1];
	  for(vpos[2]=a[2];vpos[2]<=b[2];vpos[2]++){
	   //corner id
	   //set id **Corner Size is (n+1)*(n+1)*(n+1) !!!
	   //id=POS2ID(vpos[0],vpos[1],vpos[2],v->size[0]+1,v->size[1]+1,v->size[2]+1);
	   id=i2+vpos[2];
	   CornerCoords(vcd,vpos,v->base,v->step);
	   //distance
	   dist=XyzDist(p->xyz[i],vcd);
	   //sum on corner value
	   v->cv[id]+=MetaBallFunc(dist,AtomRadius[p->TypeAtomId[i]],C);
	   //printf("%d:r=%.2f dist=%.3f %.2f %.2f %.2f <=> %.2f %.2f %.2f = %f\n",id,AtomRadius[p->TypeAtomId[i]],dist,p->xyz[i][0],p->xyz[i][1],p->xyz[i][2],vcd[0],vcd[1],vcd[2],v->cv[id]);

	  }
	 }
	}

 }
 return TRUE;
}



double MetaBallFunc(double d,double r, double C){
 double val;
 double cr=C*r;

 if(cr <d)
  return 0;
 //set your own function here!!
 val=((cr-d)*(cr-d))/((cr - r)*(cr - r));
 //printf("d= %f r= %f c= %f v=%f\n",d,r,C,val);
 //============================

 return val;
}

int CornerCoords(double *cd,UINT *pos,double *base,double step){
 cd[0]=base[0] + step * pos[0];
 cd[1]=base[1] + step * pos[1];
 cd[2]=base[2] + step * pos[2];
 return TRUE;
}

int FindNearBase(UINT *a,UINT *b,double *cd1,double *cd2,double *base,double step){

 a[0]=(UINT)((cd1[0]-base[0])/step);
 a[1]=(UINT)((cd1[1]-base[1])/step);
 a[2]=(UINT)((cd1[2]-base[2])/step);

 b[0]=(UINT)((cd2[0]-base[0])/step);
 b[1]=(UINT)((cd2[1]-base[1])/step);
 b[2]=(UINT)((cd2[2]-base[2])/step);

 return TRUE;
}

double XyzDist(double *a,double *b){
 return sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2]));
}


int DetVoxClass(VOXEL *v){
 unsigned int class=0;
 UINT i1,i2,i3,i4;
 UINT id1,id2;
 UINT add,cnt=0;
 //detect inside(true) or outside(false)
 for(i1=0;i1<(v->size[0]+1)*(v->size[1]+1)*(v->size[2]+1);i1++)
  if(v->cv[i1] < 1)
   v->corner[i1]=0;
  else
   v->corner[i1]=1;


 //class
 for(i1=0;i1<v->size[0];i1++){
  for(i2=0;i2<v->size[1];i2++){
   for(i3=0;i3<v->size[2];i3++){
    class=0;
    id1=POS2ID(i1,i2,i3,v->size[0],v->size[1],v->size[2]); //voxel id
    add=1;
	for(i4=0;i4<8;i4++){
	 //corner id
	 id2=POS2ID(i1+CornerPosTable[i4][0],i2+CornerPosTable[i4][1],i3+CornerPosTable[i4][2],v->size[0]+1,v->size[1]+1,v->size[2]+1);
	 class+=add*v->corner[id2];
	 add*=2;
	//printf("%d %d %d id=%d id2=%d/%d class= %d add= %d\n",i1,i2,i3,id1,id2,v->corner[id2],class,add);
	}
    v->cl[id1]=class;
    if(class>0 && class<255){
	//printf("%d class= %d\n",id1,class);
     cnt++;
    }else if(class>300){
	//printf("%d %d %d id=%d id2=%d/%d class= %d add= %d\n",i1,i2,i3,id1,id2,(v->size[0]+1)*(v->size[1]+1)*(v->size[2]+1),class,add);
     //return 0;

    }
   }
  }
 }

 //printf("#Total Active Voxel = %d / %d %f\n",cnt,v->size[0]*v->size[1]*v->size[2],(double)cnt/(double)(v->size[0]*v->size[1]*v->size[2]));


 //puts("#FIN DetClass");
 return cnt;
}

//int SetupHp(VOXEL *v,UINT **hp,UINT *hptbl){
int SetupHp(VOXEL *v,HP *h){
 UINT i,cnt=0,n=0;
 UINT tmp[20];
 UINT Nv= v->size[0]* v->size[1]*v->size[2];

 while(Nv >0){
  cnt+=Nv;
  tmp[n]=Nv;
  printf("#R%d %d %d\n",n,cnt,Nv);
  Nv=Nv/8;
  n++;
 }
 if((h->hptbl=(UINT *)malloc(sizeof(UINT)*n))==NULL){
  return -1;
 }
 if((h->hp=(UINT **)malloc(sizeof(UINT *)*n))==NULL)
  return -1;
 for(i=0;i<n;i++){
  h->hptbl[i]=tmp[i];
 	if((h->hp[i]=(UINT *)calloc(h->hptbl[i],sizeof(UINT)))==NULL){
 	 return -1;
 	}
 }

 //set base data
 UINT a=0;
 for(i=0;i<h->hptbl[0];i++){
  //h->hp[0][i]=NumberOfTriOnClass[v->cl[i]];//triangle base
  h->hp[0][i]=NumberOfTriOnClass[v->cl[i]];//voxel base
  //a+=hp[0][i];
 }
 //reduce data
 //base -> 2nd layer -> 3rd -> 4th ...
 UINT s=v->size[0];
 UINT x,y,z,p,id,NewId;
 for(i=0;i<n-1;i++){
  printf("L%d s=%d\n",i,s);
	for(x=0;x*2<s;x++){
	 for(y=0;y*2<s;y++){
	  for(z=0;z*2<s;z++){
	   for(p=0;p<8;p++){
	    id=POS2ID(2*x+CornerPosTable[p][0],2*y+CornerPosTable[p][1],2*z+CornerPosTable[p][2],s,s,s);
	    NewId=POS2ID(x,y,z,s/2,s/2,s/2);
	    h->hp[i+1][NewId]+=h->hp[i][id];
	    //printf("%d/%d %d/%d -> %d\n",id,hptbl[i],NewId,hptbl[i+1],hp[i][id]);
	   }
	  }
	 }
	}
  s/=2;
 }
 //printf("#Ntri on TopLayer= %d / %d\n",hp[n-1][0],a);
 h->Ntri=h->hp[n-1][0];
 h->n=n;
 return n;
}

int AssignVertex(HP *h,VOXEL *v){
 int i1,i2,i3;
 int tmp;
 printf("#Ntri= %d\n",h->Ntri);
 for(i1=0;i1<h->Ntri;i1++){
  tmp=i1;
	//search
  	for(i2=h->n-1;i2 >=0;i2--){
	 printf("%d %d\n",i1,i2);
	 printf("%d :",h->hp[i2][0]);
	 printf("%d :",h->hp[i2][1]);
	 printf("%d :",h->hp[i2][2]);
	 printf("%d :",h->hp[i2][3]);
	 printf("%d :",h->hp[i2][4]);
	 printf("%d :",h->hp[i2][5]);
	 printf("%d :",h->hp[i2][6]);
	 printf("%d\n",h->hp[i2][7]);
		for(i3=0;i3<8;i3++){
		 printf("#Hp %d %d\n",i2,h->hp[i2][i3]);
		 tmp-=h->hp[i2][i3];
		 if(tmp<0){
		  printf("#Find s= %d p= %d\n",i2,i3);
		  tmp+=h->hp[i2][i3];
		  break;
		 }
		}
	}
 }
 return TRUE;
}

int SetupVoxelTbl(unsigned int *tbl,unsigned int *Nt,VOXEL *v){
 UINT cnt=0;
 UINT i1,i2,i3;
 UINT id,tmp=0;
 Nt[0]=0;
 for(i1=0;i1<v->size[0];i1++){
  for(i2=0;i2<v->size[1];i2++){
   for(i3=0;i3<v->size[2];i3++){
    id=POS2ID(i1,i2,i3,v->size[0],v->size[1],v->size[2]);
    if(v->cl[id] !=0 && v->cl[id] !=255){
     tbl[cnt*3]=i1;
     tbl[cnt*3+1]=i2;
     tbl[cnt*3+2]=i3;
	//printf("tmp= %d id= %d\n",tmp,v->cl[id]);
     Nt[cnt+1]=tmp+NrNumberOfVertexOnClass[v->cl[id]];
     tmp=Nt[cnt+1];
     //printf("class=%d %d\n",NrNumberOfVertexOnClass[v->cl[id]],cnt*3+2);
     cnt++;
    }
   }
  }
 }
 return cnt;
}


//very simple
double FindVertexValue(double a, double b,double cutoff){
 //base
 //a-------|(val)------->b

 double val;
 val=(cutoff-a)/(b-a);
 return val;
}


int SetupVertexValue(unsigned int *t,unsigned int  n,VOXEL *v){
 UINT i1,i2,i3;
 int id,class,cid1,cid2;
 int p1,p2;

 //if((v->vertex=(double *)malloc(sizeof(double)*12*n))==NULL)
 // return FALSE;//
 //if((v->vcd=(double *)malloc(sizeof(double)*36*n))==NULL)
 // return FALSE;//

 for(i1=0;i1<n;i1++){
  id=POS2ID(t[3*i1],t[3*i1+1],t[3*i1+2],v->size[0],v->size[1],v->size[2]);

  //printf("%d %d %d %d ",i1,t[3*i1],t[3*i1+1],t[3*i1+2]);
  //printf("Class= %d\n",v->cl[id]);
  class=v->cl[id];
  //Set vertex
  for(i2=0;i2<NrNumberOfVertexOnClass[class];i2++){
   //reduce cal-------
   //if(NrVertexOnClass[class][i2] != 0 && NrVertexOnClass[class][i2] != 3 && NrVertexOnClass[class][i2] !=8)//edge 0,3,8 only
   // continue;
   //-----------------
   p1=VertexOnEdge[NrVertexOnClass[class][i2]][0];
   p2=VertexOnEdge[NrVertexOnClass[class][i2]][1];
   //printf("%d => %d-%d ",NrVertexOnClass[class][i2],p1,p2);
   //corner id
   cid1=POS2ID(t[3*i1]+CornerPosTable[p1][0],t[3*i1+1]+CornerPosTable[p1][1],t[3*i1+2]+CornerPosTable[p1][2],v->size[0]+1,v->size[1]+1,v->size[2]+1);
   cid2=POS2ID(t[3*i1]+CornerPosTable[p2][0],t[3*i1+1]+CornerPosTable[p2][1],t[3*i1+2]+CornerPosTable[p2][2],v->size[0]+1,v->size[1]+1,v->size[2]+1);

   //show value
   v->vertex[i1*12+NrVertexOnClass[class][i2]]=FindVertexValue(v->cv[cid1],v->cv[cid2],1.00);
   //v->vertex[i1*12+NrVertexOnClass[class][i2]]=0.5;
   //printf("%f - %f (%f)\n",v->cv[cid1],v->cv[cid2],v->vertex[i1*12+i2]);
  }	
 }

 return 0;
}

//a->b
void UnitVector(double *a,double *b,double *c){
 double d=XyzDist(a,b);
 c[0]=(b[0]-a[0])/d;
 c[1]=(b[1]-a[1])/d;
 c[2]=(b[2]-a[2])/d;
}

int RenderTo3D(unsigned int *t,unsigned int *Nt,unsigned int  n,VOXEL *v){
 UINT i1,i2,i3;
 int id,class,cid1,cid2;
 int p1,p2;
 double shift;
 double u[3],zero[3];

 for(i1=0;i1<n;i1++){
  id=POS2ID(t[3*i1],t[3*i1+1],t[3*i1+2],v->size[0],v->size[2],v->size[2]);
  class=v->cl[id];
  zero[0]=t[3*i1  ]*v->step + v->base[0];
  zero[1]=t[3*i1+1]*v->step + v->base[1];
  zero[2]=t[3*i1+2]*v->step + v->base[2];

  //show vertex value
   	//convert xyz coords
	for(i2=0;i2<NrNumberOfVertexOnClass[class];i2++){
  	 p1=VertexOnEdge[NrVertexOnClass[class][i2]][0];
  	 p2=VertexOnEdge[NrVertexOnClass[class][i2]][1];

	 shift=v->step*v->vertex[i1*12+NrVertexOnClass[class][i2]];

	 u[0]=(double)(CornerPosTable[p2][0]-CornerPosTable[p1][0]);
	 u[1]=(double)(CornerPosTable[p2][1]-CornerPosTable[p1][1]);
	 u[2]=(double)(CornerPosTable[p2][2]-CornerPosTable[p1][2]);
	 //printf("%f * (%f %f %f) + edge %d\n",shift,u[0],u[1],u[2],NrVertexOnClass[class][i2]);
	 //Xyz!!
	 v->vcd[i1*36 +  NrVertexOnClass[class][i2]*3 + 0]= u[0]*shift + zero[0] + CornerPosTable[p1][0]*v->step;
	 v->vcd[i1*36 +  NrVertexOnClass[class][i2]*3 + 1]= u[1]*shift + zero[1] + CornerPosTable[p1][1]*v->step;
	 v->vcd[i1*36 +  NrVertexOnClass[class][i2]*3 + 2]= u[2]*shift + zero[2] + CornerPosTable[p1][2]*v->step;
	 //printf("=> (%f %f %f)\n",v->vcd[i1*36 +  NrVertexOnClass[class][i2] + 0]);

	}
 }

 return 0;
}

//malloc and init
int InitPly(PLY *ply,int Nr){
/*
typedef struct{
 double **vcd;//coords of vertex
 int **res,**atm;
 int **tri;//triangle
 int **near;//near atom or residue for tri
 int *Nrti,*Nv;
 unsigned int TotalNtri,TotalNv;
} PLY;

*/
 if((ply->vcd=(double **)malloc(sizeof(double *)*Nr))==NULL){
  printf("#cannnot malloc for vcdl\n");
  free(ply->vcd);
  return FALSE;
 }
 if((ply->tri=(int **)malloc(sizeof(int *)*Nr))==NULL){
  printf("#cannnot malloc for tri\n");
  free(ply->tri);
  return FALSE;
 }
 if((ply->near=(int **)malloc(sizeof(int *)*Nr))==NULL){
  printf("#cannnot malloc for near\n");
  free(ply->near);
  return FALSE;
 }
 if((ply->Nv=(int *)malloc(sizeof(int)*Nr))==NULL){
  printf("#cannnot malloc for Nv\n");
  free(ply->Nv);
  return FALSE;
 }
 if((ply->Ntri=(int *)malloc(sizeof(int)*Nr))==NULL){
  printf("#cannnot malloc for Ntri\n");
  free(ply->Ntri);
  return FALSE;
 }

 ply->TotalNtri=0;
 ply->TotalNv=0;
 return TRUE;

}

int RenderToPly(unsigned int r,unsigned int *t,unsigned int *Nt,unsigned int  n,VOXEL *v,PLY *ply){
 UINT i1,i2,i3;
 int id,class,cid1,cid2;
 int p1,p2;
 double shift;
 double u[3],zero[3];

 UINT Nv=0,Ntri=0;
 UINT cnt=0;
 //unsigned int *Vtbl;

 //face(triangle) info
 int p[3];
 UINT NewID,SearchKey;


 for(i1=0;i1<n;i1++){
  id=POS2ID(t[3*i1],t[3*i1+1],t[3*i1+2],v->size[0],v->size[1],v->size[2]);
  class=v->cl[id];
  zero[0]=t[3*i1  ]*v->step + v->base[0];
  zero[1]=t[3*i1+1]*v->step + v->base[1];
  zero[2]=t[3*i1+2]*v->step + v->base[2];
  
  //reduce data=========================
  Nv+=NrNumberOf038VertexOnClass[class];
  Ntri+=NumberOfTriOnClass[class];

  	//set vertex value
   	//convert xyz coords
	for(i2=0;i2<NrNumberOfVertexOnClass[class];i2++){
  	 p1=VertexOnEdge[NrVertexOnClass[class][i2]][0];
  	 p2=VertexOnEdge[NrVertexOnClass[class][i2]][1];

	 shift=v->step*v->vertex[i1*12+NrVertexOnClass[class][i2]];

	 u[0]=(double)(CornerPosTable[p2][0]-CornerPosTable[p1][0]);
	 u[1]=(double)(CornerPosTable[p2][1]-CornerPosTable[p1][1]);
	 u[2]=(double)(CornerPosTable[p2][2]-CornerPosTable[p1][2]);
	 //printf("%f * (%f %f %f) + edge %d\n",shift,u[0],u[1],u[2],NrVertexOnClass[class][i2]);
	 //Xyz!!
	 v->vcd[i1*36 +  NrVertexOnClass[class][i2]*3 + 0]= u[0]*shift + zero[0] + CornerPosTable[p1][0]*v->step;
	 v->vcd[i1*36 +  NrVertexOnClass[class][i2]*3 + 1]= u[1]*shift + zero[1] + CornerPosTable[p1][1]*v->step;
	 v->vcd[i1*36 +  NrVertexOnClass[class][i2]*3 + 2]= u[2]*shift + zero[2] + CornerPosTable[p1][2]*v->step;
	 //printf("=> (%f %f %f)\n",v->vcd[i1*36 +  NrVertexOnClass[class][i2] + 0]);

	}
	for(i2=0;i2<NrNumberOf038VertexOnClass[class];i2++){
  	 v->Vtbl[3*id+Nr038VertexOnClass[class][i2]]=cnt;
  	 cnt++;
  	}
 }
 //malloc
/*
typedef struct{
 double **vcd;//coords of vertex
 int **res,**atm;
 int **tri;//triangle
 int **near;//near atom or residue for tri
 int *Nrti,*Nv;
 unsigned int TotalNtri,TotalNv;
} PLY;

*/
 if((ply->vcd[r]=(double *)malloc(sizeof(double)*3*Nv))==NULL){
  printf("#cannnot malloc for Vtbl\n");
  free(ply->vcd[r]);
  return FALSE;
 }
 if((ply->tri[r]=(int *)malloc(sizeof(int)*3*Ntri))==NULL){
  printf("#cannnot malloc for tri\n");
  free(ply->tri[r]);
  return FALSE;
 }
 //END malloc

 ply->Ntri[r]=Ntri;
 ply->Nv[r]=Nv;
 ply->TotalNtri+=Ntri;
 ply->TotalNv+=Nv;
 //printf("Nv= %d Ntri= %d\n",Nv,Ntri);
 cnt=0;
 UINT cnt2=0;
 //INPUT DATA TO PLY
 for(i1=0;i1<n;i1++){
  id=POS2ID(t[3*i1],t[3*i1+1],t[3*i1+2],v->size[0],v->size[1],v->size[2]);
  class=v->cl[id];
	for(i2=0;i2<NrNumberOfVertexOnClass[class];i2++){
	 if(NrVertexOnClass[class][i2]!=0 && NrVertexOnClass[class][i2]!=3 && NrVertexOnClass[class][i2]!=8 )//reduced vertex
  	  continue;
	 //for meshlab
	 ply->vcd[r][3*cnt   ]=v->vcd[i1*36 +  NrVertexOnClass[class][i2]*3   ];
	 ply->vcd[r][3*cnt+1]=v->vcd[i1*36 +  NrVertexOnClass[class][i2]*3 + 1];
	 ply->vcd[r][3*cnt+2]=v->vcd[i1*36 +  NrVertexOnClass[class][i2]*3 + 2];

  	 //printf("%.3f %.3f %.3f",v->vcd[i1*36 +  NrVertexOnClass[class][i2]*3 + 2],v->vcd[i1*36 + NrVertexOnClass[class][i2]*3 +1],v->vcd[i1*36 + NrVertexOnClass[class][i2]*3 + 0]);
	 //color alpha
 	 //printf(" %d %d %d %d\n",255*t[3*i1]/v->size[0],255*t[3*i1+1]/v->size[1],255*t[3*i1+2]/v->size[2],255);
	 cnt++;
	}

	for(i2=0;i2<NumberOfTriOnClass[class];i2++){
	 //printf("Start= %d ",Nt[i1]);

	 //printf("3");
	 //reduced vertex!! <1>
	 p[0]=ReducedVertexOnEdge[triTable[class][i2*3]][0];
	 p[1]=ReducedVertexOnEdge[triTable[class][i2*3]][1];
	 p[2]=ReducedVertexOnEdge[triTable[class][i2*3]][2];

	 NewID=POS2ID(t[3*i1]+p[0],t[3*i1+1]+p[1],t[3*i1+2]+p[2],v->size[0],v->size[1],v->size[2]);
	 SearchKey=3*NewID + ReducedVertexOnEdge[triTable[class][i2*3]][3]/3;
	 
	 //printf(" %d",v->Vtbl[SearchKey]);
	 ply->tri[r][cnt2*3+2]=v->Vtbl[SearchKey];//please check converted 2<->0

	 //reduced vertex!! <2>
	 p[0]=ReducedVertexOnEdge[triTable[class][i2*3+1]][0];
	 p[1]=ReducedVertexOnEdge[triTable[class][i2*3+1]][1];
	 p[2]=ReducedVertexOnEdge[triTable[class][i2*3+1]][2];

	 NewID=POS2ID(t[3*i1]+p[0],t[3*i1+1]+p[1],t[3*i1+2]+p[2],v->size[0],v->size[1],v->size[2]);
	 SearchKey=3*NewID + ReducedVertexOnEdge[triTable[class][i2*3+1]][3]/3;
	 
	 ply->tri[r][cnt2*3+1]=v->Vtbl[SearchKey];
	 //printf(" %d",v->Vtbl[SearchKey]);

	 //reduced vertex!! <3>
	 p[0]=ReducedVertexOnEdge[triTable[class][i2*3+2]][0];
	 p[1]=ReducedVertexOnEdge[triTable[class][i2*3+2]][1];
	 p[2]=ReducedVertexOnEdge[triTable[class][i2*3+2]][2];

	 NewID=POS2ID(t[3*i1]+p[0],t[3*i1+1]+p[1],t[3*i1+2]+p[2],v->size[0],v->size[1],v->size[2]);
	 SearchKey=3*NewID + ReducedVertexOnEdge[triTable[class][i2*3+2]][3]/3;
	 
	 ply->tri[r][cnt2*3  ]=v->Vtbl[SearchKey];//Please check converted 0<->2
	 //printf(" %d\n",v->Vtbl[SearchKey]);
	 cnt2++;
	}
 }
 if(cnt != Nv){
  printf("##ERROR in Vertex coords###\n");
  return 0;
 }
 if(cnt2!= Ntri){
  printf("##ERROR in triangle data###\n");
  return 0;
 }
 //printf("cnt=%d\n",cnt);
 return 0;
}



int cmp_uint(const void *c1, const void *c2){
 UINT a1=*(UINT *)c1;
 UINT a2=*(UINT *)c2;
 return(a1-a2);
}


int ShowPlyFormat(unsigned int *t,unsigned int *Nt,unsigned int  n,VOXEL *v){
 UINT i1,i2,i3;
 int id,class,cid1,cid2;
 int p1,p2;
 double shift;
 double u[3],zero[3];
 UINT Nv=0,Ntri=0;
 UINT cnt=0;
 unsigned int *Vtbl;

 if((Vtbl=(UINT *)malloc(sizeof(UINT)*(v->size[0])*(v->size[1])*(v->size[2])*3))==NULL){
  printf("#cannnot malloc for Vtbl\n");
  free(Vtbl);
  return FALSE;
 }

 //if((v->key=(UINT *)malloc(sizeof(UINT)*n*3))==NULL)
 // return FALSE;//

 //check triangle and element
 for(i1=0;i1<n;i1++){
  id=POS2ID(t[3*i1],t[3*i1+1],t[3*i1+2],v->size[0],v->size[1],v->size[2]);
  class=v->cl[id];
  //reduce data=========================
  Nv+=NrNumberOf038VertexOnClass[class];
  //************************************
  //Nv+=NrNumberOfVertexOnClass[class];

  Ntri+=NumberOfTriOnClass[class];
  //Set vertex order for data reduce
  for(i2=0;i2<NrNumberOf038VertexOnClass[class];i2++){
   //v->key[cnt]=id*3+i2;
   //printf("class= %d id= %d key= %d %d i2= %d\n",class,id,cnt,v->key[cnt],i2);
   Vtbl[3*id+Nr038VertexOnClass[class][i2]]=cnt;
   cnt++;
  }
 }
 

 //Header
 puts("ply\nformat ascii 1.0");
 puts("comment HPMC generated");
 printf("element vertex %d\n",Nv);
 puts("property float x\nproperty float y\nproperty float z");
 puts("property uchar red\nproperty uchar green\nproperty uchar blue\nproperty uchar alpha");
 printf("element face %d\n",Ntri);
 puts("property list uchar int vertex_indices");
 puts("end_header");

 //vertex info
 for(i1=0;i1<n;i1++){
  id=POS2ID(t[3*i1],t[3*i1+1],t[3*i1+2],v->size[0],v->size[1],v->size[2]);
  class=v->cl[id];
	for(i2=0;i2<NrNumberOfVertexOnClass[class];i2++){
	 if(NrVertexOnClass[class][i2]!=0 && NrVertexOnClass[class][i2]!=3 && NrVertexOnClass[class][i2]!=8 )//reduced vertex
  	  continue;
	 //for meshlab
  	 printf("%.3f %.3f %.3f",v->vcd[i1*36 +  NrVertexOnClass[class][i2]*3 + 2],v->vcd[i1*36 + NrVertexOnClass[class][i2]*3 +1],v->vcd[i1*36 + NrVertexOnClass[class][i2]*3 + 0]);
	 //color alpha
 	 printf(" %d %d %d %d\n",255*t[3*i1]/v->size[0],255*t[3*i1+1]/v->size[1],255*t[3*i1+2]/v->size[2],255);
	}
 }

 //face(triangle) info
 int p[3];
 UINT NewID,SearchKey;
 UINT *point;

 for(i1=0;i1<n;i1++){
  id=POS2ID(t[3*i1],t[3*i1+1],t[3*i1+2],v->size[0],v->size[1],v->size[2]);
  class=v->cl[id];
	for(i2=0;i2<NumberOfTriOnClass[class];i2++){
	 //printf("Start= %d ",Nt[i1]);

	 printf("3");
	 //reduced vertex!! <1>
	 p[0]=ReducedVertexOnEdge[triTable[class][i2*3]][0];
	 p[1]=ReducedVertexOnEdge[triTable[class][i2*3]][1];
	 p[2]=ReducedVertexOnEdge[triTable[class][i2*3]][2];

	 NewID=POS2ID(t[3*i1]+p[0],t[3*i1+1]+p[1],t[3*i1+2]+p[2],v->size[0],v->size[1],v->size[2]);
	 SearchKey=3*NewID + ReducedVertexOnEdge[triTable[class][i2*3]][3]/3;
	 
	 printf(" %d",Vtbl[SearchKey]);

	 //reduced vertex!! <2>
	 p[0]=ReducedVertexOnEdge[triTable[class][i2*3+1]][0];
	 p[1]=ReducedVertexOnEdge[triTable[class][i2*3+1]][1];
	 p[2]=ReducedVertexOnEdge[triTable[class][i2*3+1]][2];

	 NewID=POS2ID(t[3*i1]+p[0],t[3*i1+1]+p[1],t[3*i1+2]+p[2],v->size[0],v->size[1],v->size[2]);
	 SearchKey=3*NewID + ReducedVertexOnEdge[triTable[class][i2*3+1]][3]/3;
	 
	 printf(" %d",Vtbl[SearchKey]);

	 //reduced vertex!! <3>
	 p[0]=ReducedVertexOnEdge[triTable[class][i2*3+2]][0];
	 p[1]=ReducedVertexOnEdge[triTable[class][i2*3+2]][1];
	 p[2]=ReducedVertexOnEdge[triTable[class][i2*3+2]][2];

	 NewID=POS2ID(t[3*i1]+p[0],t[3*i1+1]+p[1],t[3*i1+2]+p[2],v->size[0],v->size[1],v->size[2]);
	 SearchKey=3*NewID + ReducedVertexOnEdge[triTable[class][i2*3+2]][3]/3;
	 
	 printf(" %d\n",Vtbl[SearchKey]);



	 //printf("class= %d Skey=%d id= %d *3 + %d ",class,SearchKey,NewID,ReducedVertexOnEdge[triTable[class][i2*3]][3]);
	 //point=(UINT *)bsearch(&SearchKey,v->key,cnt,sizeof(UINT),cmp_uint);
	 //printf("point= %d, %d\n",(point-(v->key)),point);
	 //printf("3 %d %d %d\n",Nt[i1]+NrtriTable[class][i2*3],Nt[i1]+NrtriTable[class][i2*3+1],Nt[i1]+NrtriTable[class][i2*3+2]);
	}
	//puts("//");
 }
 free(Vtbl);
 return TRUE;
}
void Val2Rainbow(int c[3],double val,double max){
 //Blue -> cyan -> green -> yellow -> orange -> red
 double s=val/max;
 //printf("s= %f %f / %f\n",s,val,max);
 if(s>=0 && s<1.00/4.00){
  //blue->cyan
  c[0]=0;
  c[1]=(int)255*(s*4.0);
  c[2]=255;
 }else if(s <2.00/4.00){
//cyan->green
  c[0]=0;
  c[1]=(int)255;
  c[2]=255-(int)255*((s-1.0/4.0)*4.0);
 }else if(s <3.00/5.00){
//green->yellow
  c[0]=(int)255*((s-2.0/4.0)*4.0);
  c[1]=(int)255;
  c[2]=0;
 }else if(s <=1.00){
//yellow->red
  c[0]=255;
  c[1]=(int)255-(int)255*((s-3.0/4.0)*4.0);
  c[2]=0;
 }else{
  //white
  c[0]=c[1]=c[2]=255;
 }
}

int OutPly(int *tbl,PLY *p,int c){
 UINT i1,i2;
 puts("ply\nformat ascii 1.0");
 puts("comment HPMC generated");
 printf("element vertex %d\n",p->TotalNv);
 puts("property float x\nproperty float y\nproperty float z");
 //puts("property uchar red\nproperty uchar green\nproperty uchar blue\nproperty uchar alpha");
 printf("element face %d\n",p->TotalNtri);
 puts("property list uchar int vertex_indices");
 puts("property uchar red");
 puts("property uchar green");
 puts("property uchar blue");
 puts("property uchar alpha");
 puts("end_header");
 int color[3];
 //Vertex data
 for(i1=0;i1<p->Nr;i1++){
  for(i2=0;i2<p->Nv[i1];i2++){
   printf("%.3f %.3f %.3f",p->vcd[i1][i2*3],p->vcd[i1][i2*3+1],p->vcd[i1][i2*3+2]);
   //color data
   //printf(" %d %d %d %d",color[0],color[1],color[2],255);
   //printf(" %d %d %d %d",255,255,255,255);
   printf("\n");
  }
 }
 //tri data
 UINT tmp=0;
 for(i1=0;i1<p->Nr;i1++){
  for(i2=0;i2<p->Ntri[i1];i2++){
   printf("3 %d %d %d",tmp+p->tri[i1][i2*3],tmp+p->tri[i1][i2*3+1],tmp+p->tri[i1][i2*3+2]);
   /*
   //0:N 1:C 2:O 3:S 4:P
   color[0]=color[1]=color[2]=0;
   if(p->atm[i1][i2]==0)//N blue
    color[2]=255;
   if(p->atm[i1][i2]==1)//C green
    color[1]=255;
   if(p->atm[i1][i2]==2)//O red
    color[0]=255;
   if(p->atm[i1][i2]==3){//S red
    color[0]=255;
    color[1]=255;
   }
*/
   if(c==2 && p->near[i1][i2]>=0){
    Val2Rainbow(color,(double)tbl[p->near[i1][i2]],(double)p->Nr);
   }else if(c==1){
    Val2Rainbow(color,(double)i1,(double)p->Nr);
   }else{
    color[0]=color[1]=color[2]=255;
   }

   printf(" %d %d %d %d",color[0],color[1],color[2],255);
   printf("\n");
  }
  tmp+=p->Nv[i1];
 }
}



/*

void Val2Rainbow(int c[3],double val,double max){
 //Blue -> cyan -> green -> yellow -> orange -> red
 double s=val/max;
 //printf("s= %f %f / %f\n",s,val,max);
 if(s>=0 && s<1.00/6.00){
//red 255,0,0 ->yellow 255,255,0
  c[0]=255;
  c[1]=(int)255*(s*6.0);
  c[2]=0;
 }else if(s <2.00/6.00){
//yellow->green
  c[0]=255-(int)255*((s-1.0/6.0)*6.0);
  c[1]=(int)255;
  c[2]=0;
 }else if(s <3.00/6.00){
//green->cyan
  c[0]=0;
  c[1]=(int)255;
  c[2]=(int)255*((s-2.0/6.0)*6.0);
 }else if(s <4.00/6.00){
//cyan->blue
  c[0]=0;
  c[1]=(int)255-(int)255*((s-3.0/6.0)*6.0);
  c[2]=(int)255;
 }else if(s <5.00/6.00){
//blue->magenta
  c[0]=(int)255*((s-4.0/6.0)*6.0);
  c[1]=0;
  c[2]=(int)255;
 }else if(s <=1.00){
 //magenta->red
  c[0]=(int)255;
  c[1]=0;
  c[2]=(int)255-255*((s-5.0/6.0)*6.0);
  //printf("s= %f c2= %d\n",s,c[2]);
 }else{
  //white
  c[0]=c[1]=c[2]=255;
 }
}
*/


int OutPlyN(int *tbl,PLY *p,int n,int c){
 UINT i1,i2;
 puts("ply\nformat ascii 1.0");
 puts("comment HPMC generated");
 //printf("element vertex %d\n",p->TotalNv);
 printf("element vertex %d\n",p->Nv[n]);
 puts("property float x\nproperty float y\nproperty float z");
 //puts("property uchar red\nproperty uchar green\nproperty uchar blue\nproperty uchar alpha");
 //printf("element face %d\n",p->TotalNtri);
 printf("element face %d\n",p->Ntri[n]);
 puts("property list uchar int vertex_indices");
 puts("property uchar red");
 puts("property uchar green");
 puts("property uchar blue");
 puts("property uchar alpha");
 puts("end_header");
 int color[3];
 //Vertex data
 i1=n;
  for(i2=0;i2<p->Nv[i1];i2++){
   printf("%.3f %.3f %.3f",p->vcd[i1][i2*3],p->vcd[i1][i2*3+1],p->vcd[i1][i2*3+2]);
   printf("\n");
  }
 //tri data
 UINT tmp=0;
 i1=n;
  for(i2=0;i2<p->Ntri[i1];i2++){
   printf("3 %d %d %d",tmp+p->tri[i1][i2*3],tmp+p->tri[i1][i2*3+1],tmp+p->tri[i1][i2*3+2]);
   /*
   //0:N 1:C 2:O 3:S 4:P
   color[0]=color[1]=color[2]=0;
   if(p->atm[i1][i2]==0)//N blue
    color[2]=255;
   if(p->atm[i1][i2]==1)//C green
    color[1]=255;
   if(p->atm[i1][i2]==2)//O red
    color[0]=255;
   if(p->atm[i1][i2]==3){//S red
    color[0]=255;
    color[1]=255;
   }
*/
   if(c==2 && p->near[n][i2]>=0){
    Val2Rainbow(color,(double)tbl[p->near[n][i2]],(double)p->Nr);
   }else if(c==1){
    Val2Rainbow(color,(double)n,(double)p->Nr);
   }else{
    color[0]=color[1]=color[2]=255;
   }
   //printf(" %d %d %d %d",255,10,255,255);
   printf(" %d %d %d %d",color[0],color[1],color[2],255);
   printf("\n");
  }

}

int FindMaxSizeRes(PDB *p,double *size, double *base){
 int i;
 double min[3],max[3];

 max[0]=max[1]=max[2]=-99999;
 min[0]=min[1]=min[2]= 99999;

 for(i=0;i<p->NumOfAtom;i++){

  if(max[0]<p->xyz[i][0]) max[0]=p->xyz[i][0];
  if(max[1]<p->xyz[i][1]) max[1]=p->xyz[i][1];
  if(max[2]<p->xyz[i][2]) max[2]=p->xyz[i][2];

  if(min[0]>p->xyz[i][0]) min[0]=p->xyz[i][0];
  if(min[1]>p->xyz[i][1]) min[1]=p->xyz[i][1];
  if(min[2]>p->xyz[i][2]) min[2]=p->xyz[i][2];


   //printf("Check %d %d\n",i,p->AtomOnRes[i]);
  if((p->AtomOnRes[i] != p->AtomOnRes[i+1])||i==p->NumOfAtom-1){

   //re-new size
   if(max[0]-min[0]>size[0]) size[0]=max[0]-min[0];
   if(max[1]-min[1]>size[1]) size[1]=max[1]-min[1];
   if(max[2]-min[2]>size[2]) size[2]=max[2]-min[2];

   //printf("Newsize= %.3f %.3f %.3f\n",size[0],size[1],size[2]);
   base[p->AtomOnRes[i]*3  ]=min[0];
   base[p->AtomOnRes[i]*3+1]=min[1];
   base[p->AtomOnRes[i]*3+2]=min[2];

   max[0]=max[1]=max[2]=-99999;
   min[0]=min[1]=min[2]= 99999;
  }
 }
 return 0;
}

int ResResCont(PDB *p,int r1,int r2,double r){
 int i1,i2;
 for(i1=p->ResOnAtom[r1];i1<p->ResOnAtom[r1+1];i1++){
 	for(i2=p->ResOnAtom[r2];i2<p->ResOnAtom[r2+1];i2++){
	 if(r>XyzDist(p->xyz[i1],p->xyz[i2]))
	  return TRUE;
	}
 }
 return FALSE;
}


int FindNearResidues(PDB *p,int **tbl,int *ntbl,double r){
 int i1,i2,i3,i4,minid;
 double c[3];
 int v[3];
 double d,min;
 //init
 for(i1=0;i1<p->NumOfRes-1;i1++)
  ntbl[i1]=0;

 for(i1=0;i1<p->NumOfRes-1;i1++){
	for(i2=i1+1;i2<p->NumOfRes;i2++){
	 if(ResResCont(p,i1,i2,r)==TRUE){
	  //printf("Contact!! %d %d\n",i1,i2);
	  tbl[i1][ntbl[i1]]=i2;
	  ntbl[i1]++;
	  tbl[i2][ntbl[i2]]=i1;
	  ntbl[i2]++;

	 }
	}
 }
}

int FindContact(PLY *ply,PDB *p,int **tbl,int *ntbl,double sphr){
 int i1,i2,i3,minid;
 int r1;
 double c[3];
 int v[3];
 double d,min;
 //malloc for matrix
 double area;
 int NearRes;
 double SP=sphr*2;

 //printf("Nr= %d Res= %d Atom= %d\n",ply->Nr,p->NumOfRes,p->NumOfAtom);

 if((ply->area=(double **)malloc(sizeof(double *)*ply->Nr))==NULL)
  return FALSE;
 if((ply->near=(int **)malloc(sizeof(int *)*ply->Nr))==NULL)
  return FALSE;
 if((ply->atm=(int **)malloc(sizeof(int *)*ply->Nr))==NULL)
  return FALSE;

 for(i1=0;i1<ply->Nr;i1++){
  if((ply->area[i1]=(double *)calloc(ply->Nr+1,sizeof(double)))==NULL)
   return FALSE;
  if((ply->near[i1]=(int *)malloc(sizeof(double)*ply->Ntri[i1]))==NULL)
   return FALSE;

  if((ply->atm[i1]=(int *)malloc(sizeof(double)*ply->Ntri[i1]))==NULL)
   return FALSE;
  
  for(i2=0;i2<ply->Ntri[i1];i2++){
   //printf("%d- %d\n",i1,i2);
   //set coords
   //vertec id
   v[0]=ply->tri[i1][i2*3];
   v[1]=ply->tri[i1][i2*3+1];
   v[2]=ply->tri[i1][i2*3+2];

   c[0]=ply->vcd[i1][v[0]*3  ] + ply->vcd[i1][v[1]*3  ] + ply->vcd[i1][v[2]*3  ];//X
   c[1]=ply->vcd[i1][v[0]*3+1] + ply->vcd[i1][v[1]*3+1] + ply->vcd[i1][v[2]*3+1];//Y
   c[2]=ply->vcd[i1][v[0]*3+2] + ply->vcd[i1][v[1]*3+2] + ply->vcd[i1][v[2]*3+2];//Z

   c[0]/=3.00;
   c[1]/=3.00;
   c[2]/=3.00;

   //printf("tri= %d %d %d\n",v[0],v[1],v[2]);
   //printf("v1=%.3f %.3f %.3f\n",ply->vcd[i1][v[0]*3  ],ply->vcd[i1][v[0]*3+1],ply->vcd[i1][v[0]*3+2]);
   //printf("c =%.3f %.3f %.3f\n",c[0],c[1],c[2]);

	//search based atom
	minid=-1;//no contact area
	min=1000;
	for(i3=p->ResOnAtom[i1];i3<p->ResOnAtom[i1+1];i3++){
	 if(p->TypeAtomId[i3]>4)
	  continue;
	 d=XyzDist(c,p->xyz[i3]); 
	 if(min > d){
	 //if(min > d){
	  min=d;
	  minid=i3;
	 }
	}
   	//printf("%d min= %f minid= %d res= %d\n",i1,min,minid,p->AtomOnRes[minid]);
	
	//ply->atm[i1][i2]=p->TypeAtomId[minid];//Type Atom only
	ply->atm[i1][i2]=minid;//atom number
	//search near atom from pdb
	minid=-1;//no contact area
	min=0;

	//fast! consider only near residues
	for(r1=0;r1<ntbl[i1];r1++){
	 NearRes=tbl[i1][r1];
	//for(i3=0;i3<p->NumOfAtom;i3++){
	for(i3=p->ResOnAtom[NearRes];i3<p->ResOnAtom[NearRes+1];i3++){
	 if(p->AtomOnRes[i3]==i1)//ignore same residue
	  continue;
	 if(p->TypeAtomId[i3]>4)
	  continue;
   	 //printf("c =%.3f %.3f %.3f\n",c[0],c[1],c[2]);
   	 //printf("p =%.3f %.3f %.3f\n",p->xyz[i3][0],p->xyz[i3][1],p->xyz[i3][2]);
	 d=XyzDist(c,p->xyz[i3]); 
	 //d-=(sphr*2 + AtomRadius[p->TypeAtomId[i3]]);
	 d-=(SP + AtomRadius[p->TypeAtomId[i3]]);
	 if(d<0 && min > d){
	 //if(min > d){
	  min=d;
	  minid=i3;
	 }
	}
	}
   //printf("%d min= %f minid= %d res= %d\n",i1,min,minid,p->AtomOnRes[minid]);
   ply->near[i1][i2]=minid;
	//if(minid>=0)
   	// printf("p =%.3f %.3f %.3f\n",p->xyz[minid][0],p->xyz[minid][1],p->xyz[minid][2]);
   //3 length
   //0-1
   c[0]=sqrt(X(ply->vcd[i1][v[0]*3  ] - ply->vcd[i1][v[1]*3  ]) + X(ply->vcd[i1][v[0]*3+1] - ply->vcd[i1][v[1]*3+1])+X(ply->vcd[i1][v[0]*3+2] - ply->vcd[i1][v[1]*3+2]));
   //0-2
   c[1]=sqrt(X(ply->vcd[i1][v[0]*3  ] - ply->vcd[i1][v[2]*3  ]) + X(ply->vcd[i1][v[0]*3+1] - ply->vcd[i1][v[2]*3+1])+X(ply->vcd[i1][v[0]*3+2] - ply->vcd[i1][v[2]*3+2]));
   //2-1
   c[2]=sqrt(X(ply->vcd[i1][v[2]*3  ] - ply->vcd[i1][v[1]*3  ]) + X(ply->vcd[i1][v[2]*3+1] - ply->vcd[i1][v[1]*3+1])+X(ply->vcd[i1][v[2]*3+2] - ply->vcd[i1][v[1]*3+2]));

   area=HERON(c[0],c[1],c[2]);
    //ply->area[i1][ply->Nr]+=area; //For Sum mode
   if(minid<0){
    ply->area[i1][ply->Nr]+=area;
    //printf("area= %f %f c= %.3f %.3f %.3f\n",ply->area[i1][ply->Nr],area,c[0],c[1],c[2]);
   }
   else{
    ply->area[i1][p->AtomOnRes[minid]]+=area;
   //printf("area= %f %f c= %.3f %.3f %.3f\n",ply->area[i1][p->AtomOnRes[minid]],area,c[0],c[1],c[2]);
   }
  }
 }
 return 0;
}

int OutCAMtx(PLY *ply,PDB *p){
 int i1,i2;
 puts("#order1\treal1\torder2\treal2\tarea\tAA1\tAA2");
 for(i1=0;i1<ply->Nr;i1++){
  for(i2=i1+1;i2<ply->Nr;i2++){
   if(ply->area[i1][i2]==0||ply->area[i2][i1]==0)
    continue;
   printf("%d\t%d\t%d\t%d\t%.3f\t%3s\t%3s\n",i1,p->ResNum[i1],i2,p->ResNum[i2],ply->area[i1][i2],p->TypeRes[i1],p->TypeRes[i2]);
   printf("%d\t%d\t%d\t%d\t%.3f\t%3s\t%3s\n",i2,p->ResNum[i2],i1,p->ResNum[i1],ply->area[i2][i1],p->TypeRes[i2],p->TypeRes[i1]);
  }
  //exposed area
  printf("%d\t%d\t%d\t%d\t%.3f\t%3s\tEXP\n",i1,p->ResNum[i1],-1,-1,ply->area[i1][ply->Nr],p->TypeRes[i1]);//exposed
 }
 return 0;
}

//New ver
int OutCAMtxCACD(PLY *ply,PDB *p){
 int i1,i2;
 printf("#Number of Res= %d\n",p->NumOfRes);
 printf("#Number of Chain= %d\n",p->NumOfChain);
 puts("#order1\treal1\torder2\treal2\tarea\tAA1\tAA2\tCid1\tCid2");
 for(i1=0;i1<ply->Nr;i1++){
  for(i2=i1+1;i2<ply->Nr;i2++){
   if(ply->area[i1][i2]==0||ply->area[i2][i1]==0)
    continue;
   printf("%d\t%d\t%d\t%d\t%.3f\t%3s\t%3s\t",i1,p->ResNum[i1],i2,p->ResNum[i2],ply->area[i1][i2],p->TypeRes[i1],p->TypeRes[i2]);
   printf("%d\t%d\n",p->Cid[i1],p->Cid[i2]);
   printf("%d\t%d\t%d\t%d\t%.3f\t%3s\t%3s\t",i2,p->ResNum[i2],i1,p->ResNum[i1],ply->area[i2][i1],p->TypeRes[i2],p->TypeRes[i1]);
   printf("%d\t%d\n",p->Cid[i2],p->Cid[i1]);
  }
  //exposed area
  printf("%d\t%d\t%d\t%d\t%.3f\t%3s\tEXP\t",i1,p->ResNum[i1],-1,-1,ply->area[i1][ply->Nr],p->TypeRes[i1]);//exposed
   printf("%d\t%d\n",p->Cid[i1],-1);
 }
 //CAcoords
 puts("#CACD");
 for(i1=0;i1<p->NumOfRes;i1++){
  printf("CA\t%d\t%d\t%d\t%3s\t",i1,p->ResNum[i1],p->Cid[i1],p->TypeRes[i1]);
  printf("%.3f\t%.3f\t%.3f\n",p->CAxyz[i1][0],p->CAxyz[i1][1],p->CAxyz[i1][2]);
 }
 return 0;
}

