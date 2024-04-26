/*
Shibuya A3+alpha Algorithm
2011.2.14 Fast-next-skip

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
#include "rmsd.h"
#include "mc.h"

#define PDB_STRLEN 55
#define MARGIN 5.0
#define SPHERE_R 1.40

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



double *pkdata;//for Zscore
int pk_cnt;

PDB pdb1,pdb2;
CMD cmd;


char **dblist;

int count_ca(char **,int);
int CountAtom(char *);
int MallocPdb(PDB *,int);
int PdbOnVoxel(PDB *, VOXEL *,double);
int PdbResOnVoxel(PDB *,int , VOXEL *,double);
int FindMaxSizeRes(PDB *,double *,double *);
int FindNearResidues(PDB *,int **,int *,double);
int FindContact(PLY *,PDB *,int **,int *,double);
int OutCAMtx(PLY *,PDB *);
int OutCAMtxCACD(PLY *,PDB *);//<--New

int main(int argc, char **argv)
{
 //CPU time
 double t1=gettimeofday_sec();
 double t2,t3;
 //--------

 int i,j;

 if(chkcmdline(argc,argv,&cmd)==FALSE)
  return(0);

 //count Num of Atom from PDB
 int Natom=CountAtom(cmd.p1);
 //printf("#NumOfAgtom= %d\n",Natom);
 if(Natom < 1){
  printf("Cannot find any ATOM recoords in %s\n",cmd.p1);
  return 0;
 }
 //malloc and read PDB file
 PDB p;
 if(MallocPdb(&p,Natom)==-1)
  return 0;
 if(readpdb(&p,cmd.p1,Natom)==FALSE){
  printf("##ERROR in readpdb###\n");
  return 0;
 }
 //Ca coord check!! 01/19/2015
 bool *MissTbl;
 MissTbl=(bool *)malloc(sizeof(bool)*p.NumOfRes);
 ChkMissRes(&p,MissTbl);
 //Shift PDB ATOMs
 ShiftMissRes(&p,MissTbl);
 //Set Ca coords 01/19/2015
 if(SetCaCoords(&p)==-1){
  printf("##ERROR in readpdb###\n");
  return 0;
 }


 //Set up size of voxel here!!
 VOXEL vox;
 double BaseXyz[3],MaxXyz[3];
 double MaxSizeXyz[3];
 double *BaseResXyz;
 UINT SizeOfVox[3];

 if((BaseResXyz=(double *)malloc(sizeof(double)*3*p.NumOfRes))==NULL){
  return 0;
 }
 //find maximum size of voxel
 FindMaxSizeRes(&p,MaxSizeXyz,BaseResXyz);
 //printf("MaxSizeXyz= %.3f %.3f %.3f\n",MaxSizeXyz[0],MaxSizeXyz[1],MaxSizeXyz[2]);
 //printf("BaseResXyz0= %.3f %.3f %.3f\n",BaseResXyz[0],BaseResXyz[1],BaseResXyz[2]);

/*
 //modefied
 BaseXyz[0]=p.MinXyz[0] -  (cmd.MetaC*1.8)-cmd.g_step-5;
 BaseXyz[1]=p.MinXyz[1] -  (cmd.MetaC*1.8)-cmd.g_step-5;
 BaseXyz[2]=p.MinXyz[2] -  (cmd.MetaC*1.8)-cmd.g_step-5;
 MaxXyz[0] =p.MaxXyz[0] +  (cmd.MetaC*1.8)+cmd.g_step+5;
 MaxXyz[1] =p.MaxXyz[1] +  (cmd.MetaC*1.8)+cmd.g_step+5;
 MaxXyz[2] =p.MaxXyz[2] +  (cmd.MetaC*1.8)+cmd.g_step+5;

 if(SetupVoxelSize(SizeOfVox,BaseXyz,MaxXyz,cmd.g_step)==FALSE)
  return 0;
*/

 //printf("#BaseXyz= %.3f %.3f %.3f\n",BaseXyz[0],BaseXyz[1],BaseXyz[2]);

//reduce vox size
 SizeOfVox[0]=(UINT)((MaxSizeXyz[0] + MARGIN*2*cmd.MetaC)/(double)cmd.g_step)+1;
 SizeOfVox[1]=(UINT)((MaxSizeXyz[1] + MARGIN*2*cmd.MetaC)/(double)cmd.g_step)+1;
 SizeOfVox[2]=(UINT)((MaxSizeXyz[2] + MARGIN*2*cmd.MetaC)/(double)cmd.g_step)+1;

 //printf("#Size   = %d %d %d\n",SizeOfVox[0],SizeOfVox[1],SizeOfVox[2]);

 //Malloc Voxcel
 if(MallocVoxel(&vox,BaseXyz,SizeOfVox,cmd.g_step)==-1){
  printf("###Malloc Voxel Error###");
  return 0;
 }

 int NumOfActVox;
 unsigned int *VoxTbl;
 UINT *NvPerVox;

 //Assign size from res0
 //PdbResOnVoxel(&p,0,&vox,cmd.MetaC);
 //Detect class for each voxel
 //NumOfActVox=DetVoxClass(&vox);
 //NumOfActVox=NumOfActVox*30;//for safe

 NumOfActVox=SizeOfVox[0]*SizeOfVox[1]*SizeOfVox[2];
 //printf("%d\n",NumOfActVox*3);
 if((VoxTbl=(unsigned int *)malloc(sizeof(unsigned int)*NumOfActVox*3))==NULL)
  return 0;
 if((NvPerVox=(unsigned int *)malloc(sizeof(unsigned int)*(NumOfActVox+1)))==NULL)
  return 0;
 //malloc voxel
 if((vox.vertex=(double *)malloc(sizeof(double)*NumOfActVox*12))==NULL)
  return 0;
 if((vox.vcd=(double *)malloc(sizeof(double)*NumOfActVox*36))==NULL)
  return 0;
 
 if((vox.Vtbl=(UINT *)malloc(sizeof(UINT)*(vox.size[0])*(vox.size[1])*(vox.size[2])*3))==NULL){
  printf("#cannnot malloc for Vtbl\n");
  free(vox.Vtbl);
  return FALSE;
 }
 
 PLY ply;//polygon data
 if(InitPly(&ply,p.NumOfRes)==FALSE)
  return 0;
 
 //puts("START!!");
 //each residues !!
 for(i=0;i<p.NumOfRes;i++){
	//Clear voxel data
	ClearVox(&vox);
  	//printf("##FOR %d\n",i);

	//Set Basecoords in vox
	vox.base[0]=BaseResXyz[i*3  ] - MARGIN*cmd.MetaC;
	vox.base[1]=BaseResXyz[i*3+1] - MARGIN*cmd.MetaC;
	vox.base[2]=BaseResXyz[i*3+2] - MARGIN*cmd.MetaC;

	//printf("BASE= %.3f %.3f %.3f\n",vox.base[0],vox.base[1],vox.base[2]);

 	//Detect Grid value
 	PdbResOnVoxel(&p,i,&vox,cmd.MetaC);
 	//Detect class for each voxel
 	NumOfActVox=DetVoxClass(&vox);
	//printf("#NumOfAvtiveVox= %d\n",NumOfActVox);
 	if(NumOfActVox==0){
 	 printf("###No Active Voxel###\n");
 	 return 0;
 	}
	//SetupVoxelTable
	SetupVoxelTbl(VoxTbl,NvPerVox,&vox);
	SetupVertexValue(VoxTbl,NumOfActVox,&vox);
	//rednering
 	RenderToPly(i,VoxTbl,NvPerVox,NumOfActVox,&vox,&ply);
	//Find near atom
 }
 ply.Nr=p.NumOfRes;

 //Find near atom for each faces
 //Set near residue table
 int **tbl,*ntbl;
 double dcut=1.4*2 + 2*1.80*cmd.MetaC;

 if((tbl=(int **)malloc(sizeof(int *)*p.NumOfRes))==NULL)
  return 0;
 for(i=0;i<p.NumOfRes;i++)
  if((tbl[i]=(int *)malloc(sizeof(int)*p.NumOfRes))==NULL)
   return 0;

 if((ntbl=(int *)malloc(sizeof(int)*p.NumOfRes))==NULL)
  return 0;
 
 //contact distance = 1.4*2 + 2*1.80*cmd.MetaC
 FindNearResidues(&p,tbl,ntbl,dcut);
 FindContact(&ply,&p,tbl,ntbl,SPHERE_R);

 if(cmd.DemoMode==true)
  OutPlyN(p.AtomOnRes,&ply,cmd.ShowResNum,cmd.Color);

 if(cmd.PlyMode==true)
 OutPly(p.AtomOnRes,&ply,cmd.Color);
 
 if(cmd.MtxMode==true)
  OutCAMtxCACD(&ply,&p);

 t2=gettimeofday_sec();
 printf("#FINISHED TOTAL TIME= %f\n",t2-t1);
 return 0;

}



