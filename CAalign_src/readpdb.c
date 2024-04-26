#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "struct.h"
#include "func.h"

//New ver 2012.10.29 by terashi
//Custmized for marching cubes!!!

int Atom2Int(char);

int readpdb(PDB *pdb, char *filename, int n){

	int i,j; 
	FILE *fpin; 
	char line[LIN], buf[LIN]; 
	double x,y,z;
	//dock
	int no1=0,no2=0;
	int real;
	COORD sum_coord;
   	int sum_num=0,nmr=0;
	x=y=z;

	pdb->NumOfAtom=n;
	
	if((fpin=fopen(filename,"r")) == NULL){ 
		fprintf(stderr,"Can't open %s\n",filename); 
		return(-1); 
	}
	strcpy(pdb->fname,filename);//New!
	//printf("******%s\n",pdb->fname);
	//printf("#filename: %s\n",filename); 
	i=0;j=-1;
	int Nca=0;
	int PreResNum=-99999;
	int ResNum=-1;
	pdb->NumOfReal=0;
	while(fgets(line,LIN,fpin)){ 
		if(!strncmp(line,"ATOM",4) 
		 &&(!strncmp(&line[16],"A",1) 
		 ||!strncmp(&line[16]," ",1))){ 
/*
			if(!strncmp(&line[13],"N  ",3)){ 
				j++; 
				//strncpy(pdb->TypeRes[j],&line[17],3); 
				//pdb->TypeRes[j][3]='\0';
				pdb->ResNnum[j]=i;
				//initial CA CB
				//cp_cd(&pdb->CAcd[j],-999.999,-999.999,-999.999);
				//cp_cd(&pdb->CAcd[j],-999.999,-999.999,-999.999);
				pdb->AA2int_data[j]=AA2int(pdb->TypeRes[j]);
				//printf("%s %d\n",pdb->TypeRes[j],pdb->AA2int_data[j]);
			}
*/
			//Res num order
			strncpy(buf,&line[22],4);
			buf[4]='\0';
			ResNum=atoi(buf);
			//For each residue========================
			if(PreResNum != ResNum){
			 j++;
			 PreResNum=ResNum;
			 pdb->ResNum[j]=ResNum;//Real number
			 //Residue Type
			 strncpy(pdb->TypeRes[j],&line[17],3);
			 pdb->TypeRes[j][3]='\0';
			 pdb->TypeResId[j]=AA2int(pdb->TypeRes[j]);
			 pdb->ResOnAtom[j]=i;
			}
			//========================================

			//For each Atoms
			pdb->AtomOnRes[i]=j;//Order base

			//chain id
			pdb->Chain[i]=line[21];

			//Real Atom number
			strncpy(buf,&line[4],7); buf[7]='\0';
			pdb->AtomNum[i]=atoi(buf);
			strncpy(pdb->TypeAtom[i],&line[13],3);
			pdb->TypeAtom[i][3]='\0'; 

			pdb->TypeAtomId[i]=Atom2Int(pdb->TypeAtom[i][0]);			

			//xyz coords
			strncpy(buf,&line[30],8); buf[8]='\0'; 
			x=atof(buf);
			strncpy(buf,&line[38],8); buf[8]='\0'; 
			y=atof(buf);
			strncpy(buf,&line[46],8); buf[8]='\0'; 
			z=atof(buf);
			pdb->xyz[i][0]=x;
			pdb->xyz[i][1]=y;
			pdb->xyz[i][2]=z;

			//Max and Min
			if(pdb->MaxXyz[0] < x||i==0) pdb->MaxXyz[0]=x;
			if(pdb->MaxXyz[1] < y||i==0) pdb->MaxXyz[1]=y;
			if(pdb->MaxXyz[2] < z||i==0) pdb->MaxXyz[2]=z;
			if(pdb->MinXyz[0] > x||i==0) pdb->MinXyz[0]=x;
			if(pdb->MinXyz[1] > y||i==0) pdb->MinXyz[1]=y;
			if(pdb->MinXyz[2] > z||i==0) pdb->MinXyz[2]=z;

			i++;
			if(i>n)
			 return -1;
		}else if(!strncmp(line,"MODEL",5)){//NMR
    	 	 nmr++;
		 if(nmr>1)//first model only
		  break; 
		}

	}
	pdb->ResOnAtom[j+1]=i;
	pdb->AtomOnRes[i]=j+1;
	//pdb->ResNnum[j+1]=i; 
	fclose(fpin) ; 
	pdb->NumOfAtom=i; 
	pdb->NumOfRes=j+1;//New!!
	//printf("#NumOfRes = %d\n",pdb->NumOfRes);
	//printf("#NumOfAtom= %d\n",pdb->NumOfAtom);
	//printf("#Xrange %.3f - %.3f\n",pdb->MinXyz[0],pdb->MaxXyz[0]);
	//printf("#Yrange %.3f - %.3f\n",pdb->MinXyz[1],pdb->MaxXyz[1]);
	//printf("#Zrange %.3f - %.3f\n",pdb->MinXyz[2],pdb->MaxXyz[2]);
	//debug
/*
	for(i=0;i<n;i++){
	 //printf("%c",pdb->Chain[i]);
	 printf("%d(%d),",pdb->TypeAtomId[i],pdb->TypeResId[pdb->AtomOnRes[i]]);
	}
	puts("");
*/
	return(0);
}


int MallocPdb(PDB *p,int n){
 int i;
	//malloc here!!
	//xyz coord data double n*3 format
	if((p->xyz=(double **)malloc(sizeof(double *)*n))==NULL){
	 free(p->xyz);
	 return -1;
	}
	for(i=0;i<n;i++){
	 if((p->xyz[i]=(double *)malloc(sizeof(double)*3))==NULL){
	  free(p->xyz);
	  return -1;
	 }
	}
	//Residue number on each Atom
	if((p->AtomOnRes=(int *)malloc(sizeof(int *)*n+1))==NULL){
	 free(p->AtomOnRes);
	 return -1;
	}
	//Atom number on each residue
	if((p->ResOnAtom=(int *)malloc(sizeof(int *)*n+1))==NULL){
	 free(p->ResOnAtom);
	 return -1;
	}
	//Type of Atom 3char
	if((p->TypeAtom=(char **)malloc(sizeof(char*)*n))==NULL){
	 return -1;
	}
	for(i=0;i<n;i++){
	 if((p->TypeAtom[i]=(char *)malloc(sizeof(char)*4))==NULL){
	  return -1;
	 }
	}
	if((p->TypeAtomId=(int *)malloc(sizeof(int)*n))==NULL){
	 return -1;
	}
	
	//Other data
	if((p->TypeAtomId=(int*)malloc(sizeof(int)*n))==NULL){
	 free(p->TypeAtomId);
	 return -1;
	}
	if((p->TypeResId=(int*)malloc(sizeof(int)*n))==NULL){
	 free(p->TypeResId);
	 return -1;
	}

	//Residue data
	if((p->TypeRes=(char **)malloc(sizeof(char*)*n))==NULL){
	 return -1;
	}
	for(i=0;i<n;i++){
	 if((p->TypeRes[i]=(char *)malloc(sizeof(char)*4))==NULL){
	  return -1;
	 }
	}
	//Real Res number
	if((p->ResNum=(int *)malloc(sizeof(int)*n))==NULL){
	 return -1;
	}
	//Real Atom Number
	if((p->AtomNum=(int *)malloc(sizeof(int)*n))==NULL){
	 return -1;
	}
	if((p->Chain=(char *)malloc(sizeof(char)*n))==NULL){
	 return -1;
	}
	
 //puts("#fin malloc");
 return 0;
}

int SetCaCoords(PDB *p){
 int i,n1,n2;

//MALLOC
 if((p->CAxyz=(double **) malloc(sizeof(double *)*p->NumOfRes))==NULL)
  return 0;
 for(i=0;i<p->NumOfRes;i++)
  if((p->CAxyz[i]=(double *) malloc(sizeof(double)*3))==NULL)
   return 0;
//END MALLOC

 n1=0;
 for(i=0;i<p->NumOfAtom;i++){
  if(!strncmp(p->TypeAtom[i],"CA ",3)){
   p->CAxyz[n1][0]=p->xyz[i][0];
   p->CAxyz[n1][1]=p->xyz[i][1];
   p->CAxyz[n1][2]=p->xyz[i][2];
 	//printf("%d -> %f\n",n1,p->CAxyz[n1][0]);
   n1++;
  }
 }
 if(n1 != p->NumOfRes)
  return 0;
 return n1;
}

int CountAtom(char *filename){
	int i,j; 
	FILE *fpin; 
	char line[LIN], buf[LIN]; 
	double x,y,z;
	//dock
   	int nmr=0;
	x=y=z;
	
	if((fpin=fopen(filename,"r")) == NULL){ 
		fprintf(stderr,"Can't open %s\n",filename); 
		return(-1); 
	}

	//printf("#filename: %s\n",filename); 
	i=0;j=-1;
	int Natom=0;
	while(fgets(line,LIN,fpin)){ 
		if(!strncmp(line,"ATOM",4) 
		 &&(!strncmp(&line[16],"A",1) 
		 ||!strncmp(&line[16]," ",1))){ 
		 Natom++;
		}else if(!strncmp(line,"MODEL",5)){//NMR
    	 	 nmr++;
		 if(nmr>1)//first model only
		  break; 
		}
	}
 return Natom;
}




int readgri(GRID *gri,char *filename){ 
	int i,j;
        FILE *fpin; 
	char line[LIN], buf[LIN]; 
	double x,y,z; 
	x=y=z=0;
	if((fpin=fopen(filename,"r")) == NULL){ 
		fprintf(stderr,"Can't open %s\n",filename); 
		return(-1); 
	} 
	printf("grid filename: %s\n",filename); 
	i=0;
	while(fgets(line,LIN,fpin)){ 
		if(!strncmp(line,"HETATM",6)){ 
			strncpy(buf,&line[30],8); buf[8]='\0'; 
			x=atof(buf); 
			strncpy(buf,&line[38],8); buf[8]='\0'; 
			y=atof(buf); 
			strncpy(buf,&line[46],8); buf[8]='\0'; 
			z=atof(buf); 
			cp_cd(&gri->coord[i],x,y,z); 
			strncpy(buf,&line[60],6); buf[6]='\0';
			gri->sco[i]=atof(buf);
			i++; 
		} 
	} 
	printf("# of grid number:   %7d\n",i);
	gri->NumOfGri=i;
	return(0);
}

//side chain & asa det
int side_det(PDB *pdb,int num){
 if(!strncmp(pdb->TypeAtom[num],"N  ",3))
  return(1);
 if(!strncmp(pdb->TypeAtom[num],"CA ",3)
  && strncmp(pdb->TypeRes[pdb->AtomOnRes[num]],"GLY",3))
  return(1);
 if(!strncmp(pdb->TypeAtom[num],"C  ",3))
  return(1);
 if(!strncmp(pdb->TypeAtom[num],"O  ",3))
  return(1);
 if(!strncmp(pdb->TypeAtom[num],"OXT",3))
  return(1);

 return(0);
}
int AA2int(char *TypeAtom){
        if(!strcmp(TypeAtom,"ALA")){return 0;}
        if(!strcmp(TypeAtom,"VAL")){return 1;}
        if(!strcmp(TypeAtom,"PHE")){return 2;}
        if(!strcmp(TypeAtom,"PRO")){return 3;}
        if(!strcmp(TypeAtom,"MET")){return 4;}
        if(!strcmp(TypeAtom,"ILE")){return 5;}
        if(!strcmp(TypeAtom,"LEU")){return 6;}
        if(!strcmp(TypeAtom,"ASP")){return 7;}
        if(!strcmp(TypeAtom,"GLU")){return 8;}
        if(!strcmp(TypeAtom,"LYS")){return 9;}
        if(!strcmp(TypeAtom,"ARG")){return 10;}
        if(!strcmp(TypeAtom,"SER")){return 11;}
        if(!strcmp(TypeAtom,"THR")){return 12;}
        if(!strcmp(TypeAtom,"TYR")){return 13;}
        if(!strcmp(TypeAtom,"HIS")){return 14;}
        if(!strcmp(TypeAtom,"CYS")){return 15;}
        if(!strcmp(TypeAtom,"ASN")){return 16;}
        if(!strcmp(TypeAtom,"TRP")){return 17;}
        if(!strcmp(TypeAtom,"GLN")){return 18;}
        if(!strcmp(TypeAtom,"GLY")){return 19;}
        return -1;
}

int A2int(char TypeAtom){
        switch(TypeAtom){
                case 'A':return 0;break;
                case 'V':return 1;break;
                case 'F':return 2;break;
                case 'P':return 3;break;
                case 'M':return 4;break;
                case 'I':return 5;break;
                case 'L':return 6;break;
                case 'D':return 7;break;
                case 'E':return 8;break;
                case 'K':return 9;break;
                case 'R':return 10;break;
                case 'S':return 11;break;
                case 'T':return 12;break;
                case 'Y':return 13;break;
                case 'H':return 14;break;
                case 'C':return 15;break;
                case 'N':return 16;break;
                case 'W':return 17;break;
                case 'Q':return 18;break;
                case 'G':return 19;break;
                case '-':return 20;break;
                default: return -1;
        }
}

int Atom2Int(char TypeAtom){
	switch(TypeAtom){
                case 'N':return 0;break;
                case 'C':return 1;break;
                case 'O':return 2;break;
                case 'S':return 3;break;
                case 'P':return 4;break;
                case 'H':return 5;break;
                default: return 6;
        }
}
