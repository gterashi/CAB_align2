#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "struct.h"
#include "func.h"

#define PDB_STRLEN 55
#define CACADIST 4.00
#define HASH_FUNC(a,b) a%b
#define DIST(a,b,c,d,e,f) sqrt((a-b)*(a-b)+(c-d)*(c-d)+(e-f)*(e-f))
//0 :0-1 //1 :0-2 //2 :0-3 //3 :0-4 //4 :0-5
//5 :1-2 //6 :1-3 //7 :1-4 //8 :1-5
//9 :2-3 //10:2-4 //11:2-5
//12:3-4 //13:3-5
//14:4-5
static int pat[15][2]={
{0,1}, {0,2}, {0,3}, {0,4}, {0,5},
{1,2}, {1,3}, {1,4}, {1,5},
{2,3}, {2,4}, {2,5},
{3,4}, {3,5},
{4,5}
};
static int pat15[15][3]={
{0,9,14}, //0-1:2-3:4-5
{0,10,13},//0-1:2-4:3-5
{0,11,12},//0-1:2-5:3-4
{1,6,14}, //0-2:1-3:4-5
{1,7,13}, //0-2:1-4:3-5
{1,8,12}, //0-2:1-5:3-4
{2,5,14}, //0-3:1-2:4-5
{2,7,11}, //0-3:1-4:2-5
{2,8,10}, //0-3:1-5:2-4
{3,5,13}, //0-4:1-2:3-5
{3,6,11}, //0-4:1-3:2-5
{3,8,9 }, //0-4:1-5:2-3
{4,5,12}, //0-5:1-2:3-4
{4,6,10}, //0-5:1-3:2-4
{4,7,9 }, //0-5:1-4:2-3
};

extern CMD cmd;

int cmp_shib_rms(const void *c1, const void *c2){
 SHIBSET a1=*(SHIBSET *)c1;
 SHIBSET a2=*(SHIBSET *)c2;
 if(a1.rms>a2.rms)
  return 1;
 else
  return 0;
}


int input_hash(int hashid,int tag,int *tbl,int *ktbl,int n){
 int id=hashid;
 tag++;
 while(1){
  if(ktbl[id]!=tag){
   tbl[id]=hashid;
   ktbl[id]=tag;
   //printf("input %d tag=%d in hash(%d)\n",hashid,tag,id);
   break;
  }else{
   //printf("hash(%d)=%d\n",hashid,ktbl[id]);
   //printf("CANNOTinput %d tag=%d in hash(%d)\n",hashid,tag,id);
  }
  id++;
  if(id>=n)
   id=0;
 }
 return TRUE;
}
int search_hash(int hashid,int tag,int *tbl,int *ktbl,int n){
 int id=hashid;
 while(1){
  if(ktbl[id]==tag){
	if(tbl[id]==hashid)
	 return TRUE;
	else{
	 continue;
	 id++;
	}
  }else{
   return FALSE;
  }
  id++;
  if(id>=n)
   id=0;
  if(id==hashid)
   return FALSE;
 }
}

double rmsd_capdb(double (*)[3],CAPDB *,int,int,int);

void pg_bar(int n){
 if(n%100==0)
   printf(".");
 if(n%1000==0)
   printf("+");
 if(n%5000==0)
  printf("|%d\n",n);
}

int count_ca(char **list,int N){
 FILE *fp;
 char line[LIN],buf[9];
 int num=0,i;
 setbuf(stdout,NULL);
 for(i=0;i<N;i++){
  if((fp=fopen(list[i],"r"))==NULL)
   continue;
 	while(fgets(line,PDB_STRLEN,fp)!=NULL){
		if(!strncmp(&line[13],"CA ",3) && !strncmp(line,"ATOM",4) ){
		 //Ca only
 		 num++;
		}
		//NMR
		if(!strncmp(&line[0],"ENDMDL",6))
		 break;
 	}
  	fclose(fp);
	pg_bar(i);
 }
 return num;
}


//For A5 Algorithm!!
int LowerBonds(char **list,int N,int len,int LenOfSeg,FTABLE *ft,int *res_list){
 
 int i,id=0,pos,j,k,fpat;
 int pN,pC;//position N,C
 int num=0;
 FILE *fp;
 char line[LIN],buf[9];
 //int *Miss;//Missing residue flag
 double x,y,z,*gx,*gy,*gz;
 double *sx,*sy,*sz;
 //int len=m;
 int total=0,start_tot;

 setbuf(stdout,NULL);

 gx=(double *)malloc(sizeof(double)*RES);
 gy=(double *)malloc(sizeof(double)*RES);
 gz=(double *)malloc(sizeof(double)*RES);

 sx=(double *)malloc(sizeof(double)*RES);
 sy=(double *)malloc(sizeof(double)*RES);
 sz=(double *)malloc(sizeof(double)*RES);

 for(i=0;i<N;i++){
  if((fp=fopen(list[i],"r"))==NULL)
   continue;
  //printf("#Reading..%s\n",list[i]);
  num=0;
 	while(fgets(line,PDB_STRLEN,fp)!=NULL){
		if(!strncmp(&line[13],"CA ",3) && !strncmp(line,"ATOM",4) ){
		 //Ca only
		 strncpy(buf,&line[30],8); buf[8]='\0'; 
		 sx[num]=atof(buf);
		 strncpy(buf,&line[38],8); buf[8]='\0'; 
		 sy[num]=atof(buf);
		 strncpy(buf,&line[46],8); buf[8]='\0'; 
		 sz[num]=atof(buf);
 		 num++;
		}
		//NMR
		if(!strncmp(&line[0],"ENDMDL",6))
		 break;
 	}
  	fclose(fp);
	if(num<len*6){
	 pg_bar(i);
	 continue;
	}
	res_list[i]=num;
	printf("%d num=%d\n",i,num);

	//Liner cal F(S)=|G(S_left)-G(S_rignt)|/2
	//not real cog, only sum of coords
	gx[0]=gy[0]=gz[0]=0;
	for(pos=0;pos<len;pos++){
	 gx[0]+=sx[pos];
	 gy[0]+=sy[pos];
	 gz[0]+=sz[pos];
	}
	for(pos=1;pos<=num-len;pos++){
	 gx[pos]=gx[pos-1]+(sx[pos+len-1]-sx[pos-1]);
	 gy[pos]=gy[pos-1]+(sy[pos+len-1]-sy[pos-1]);
	 gz[pos]=gz[pos-1]+(sz[pos+len-1]-sz[pos-1]);
	}
	//convert to real cog
	for(pos=0;pos<=num-len;pos++){
	 gx[pos]/=(double)len;
	 gy[pos]/=(double)len;
	 gz[pos]/=(double)len;
	}
	//input ft
	start_tot=total;
	//f15[0] 0-1
	/*
        pos
        0     1      2      3      4      5
	|*****++*****++-----++-----++-----++-----|
    	   len    len    len    len    len    len      
	*/
	for(pos=0;pos<=num-2*len;pos++){
	 ft[start_tot+pos].id=i;
	 ft[start_tot+pos].pos=pos;
	 //
	 ft[start_tot+pos].f15[0]=(int)(500*sqrt((X(gx[pos]-gx[pos+len])
	 +X(gy[pos]-gy[pos+len]) +X(gz[pos]-gz[pos+len]))));

	 ft[start_tot+pos].miss=false;
	}
	//f15[1] 0-2
	/*
        pos
        0     1      2      3      4      5
	|*****++-----++*****++-----++-----++-----|
    	   len    len    len    len    len    len      
	*/

	for(pos=0;pos<=num-3*len;pos++){
	 ft[start_tot+pos].f15[1]=(int)(500*sqrt((X(gx[pos]-gx[pos+2*len])
	 +X(gy[pos]-gy[pos+2*len]) +X(gz[pos]-gz[pos+2*len]))));
	}
	//f15[2] 0-3
	/*
        pos
        0     1      2      3      4      5
	|*****++-----++-----++*****++-----++-----|
    	   len    len    len    len    len    len      
	*/
	for(pos=0;pos<=num-4*len;pos++){
	 //ft[start_tot+pos].f15[2]=(int)50*sqrt((X(gx[pos]-gx[pos+3*len])
	 //+X(gy[pos]-gy[pos+3*len]) +X(gz[pos]-gz[pos+3*len]))/(len*len));
	 ft[start_tot+pos].f15[2]=(int)(500*sqrt((X(gx[pos]-gx[pos+3*len])
	 +X(gy[pos]-gy[pos+3*len]) +X(gz[pos]-gz[pos+3*len]))));
	}
	//f15[3] 0-4
	for(pos=0;pos<=num-5*len;pos++){
	 //ft[start_tot+pos].f15[3]=(int)50*sqrt((X(gx[pos]-gx[pos+4*len])
	 //+X(gy[pos]-gy[pos+4*len]) +X(gz[pos]-gz[pos+4*len]))/(len*len));
	 ft[start_tot+pos].f15[3]=(int)(500*sqrt((X(gx[pos]-gx[pos+4*len])
	 +X(gy[pos]-gy[pos+4*len]) +X(gz[pos]-gz[pos+4*len]))));
	}
	//f15[4] 0-5
	for(pos=0;pos<=num-6*len;pos++){
	 //ft[start_tot+pos].f15[4]=(int)50*sqrt((X(gx[pos]-gx[pos+5*len])
	 //+X(gy[pos]-gy[pos+5*len]) +X(gz[pos]-gz[pos+5*len]))/(len*len));
	 ft[start_tot+pos].f15[4]=(int)(500*sqrt((X(gx[pos]-gx[pos+5*len])
	 +X(gy[pos]-gy[pos+5*len]) +X(gz[pos]-gz[pos+5*len]))));

	 ft[start_tot+pos].f=ft[start_tot+pos].f15[4];//sort key******
	}
	/*
        pos
        f            f1            f2
	|-----++-----++-----++-----++-----++-----|
    	   len    len    len    len    len    len      
	*/
	//input f1 and f2
	//for(pos=0;pos<num-6*len;pos++){
	for(pos=0;pos<=num-LenOfSeg;pos++){
	 //ft[start_tot+pos].f1=ft[start_tot+pos+2*len].f;
	 //ft[start_tot+pos].f2=ft[start_tot+pos+4*len].f;
	 //New 0-1
	 ft[start_tot+pos].f15[5] =ft[start_tot+pos+  len].f15[0];//1-2
	 ft[start_tot+pos].f15[9] =ft[start_tot+pos+2*len].f15[0];//2-3
	 ft[start_tot+pos].f15[12]=ft[start_tot+pos+3*len].f15[0];//3-4
	 ft[start_tot+pos].f15[14]=ft[start_tot+pos+4*len].f15[0];//4-5
	 //0-2
	 ft[start_tot+pos].f15[6] =ft[start_tot+pos+   len].f15[1];//1-3
	 ft[start_tot+pos].f15[10] =ft[start_tot+pos+2*len].f15[1];//2-4
	 ft[start_tot+pos].f15[13] =ft[start_tot+pos+3*len].f15[1];//3-5
	 //0-3
	 ft[start_tot+pos].f15[7] =ft[start_tot+pos+   len].f15[2];//1-4
	 ft[start_tot+pos].f15[11] =ft[start_tot+pos+2*len].f15[2];//2-5
	 //0-4
	 ft[start_tot+pos].f15[8] =ft[start_tot+pos+  len].f15[3];//1-5
	}

	total=start_tot+pos;
  pg_bar(i);
 }
 return(total);
}

//one lowerbound per one str
int LowerBonds_perstr(CAPDB *db,FTABLE *ft){
 char **list;
 int LenOfSeg=db->n[0];
 int *res_list;
 int len=(int)(LenOfSeg/6.00);
 int i,id=0,pos,j,k,fpat,N;
 int pN,pC;//position N,C
 int num=0;
 int seg;
 FILE *fp;
 char line[LIN],buf[9];
 //int *Miss;//Missing residue flag
 double x,y,z,*gx,*gy,*gz;
 double *sx,*sy,*sz;
 //int len=m;
 int total=0,start_tot;
 int stpos;
 int p,p1,p2;

 printf("len= %d LenOfSeg= %d\n",len,LenOfSeg);
 setbuf(stdout,NULL);
 N=db->Nst;
 //use only 6 segments
 //gx=(double *)malloc(sizeof(double)*RES);
 gx=(double *)malloc(sizeof(double)*6);
 gy=(double *)malloc(sizeof(double)*6);
 gz=(double *)malloc(sizeof(double)*6);

/*
 //sx=(double *)malloc(sizeof(double)*RES);
 sx=(double *)malloc(sizeof(double)*6);
 sy=(double *)malloc(sizeof(double)*6);
 sz=(double *)malloc(sizeof(double)*6);
*/
 for(i=0;i<N;i++){
  //printf("#Reading..%s\n",list[i]);
  num=0;
	//res_list[i]=num;
	//printf("%d num=%d\n",i,num);

	//Liner cal F(S)=|G(S_left)-G(S_rignt)|/2
	//not real cog, only sum of coords
	//init
	for(seg=0;seg<6;seg++)
	 gx[seg]=gy[seg]=gz[seg]=0;
	
	//gx[0]=gy[0]=gz[0]=0;
	stpos=db->pos[i];
	for(pos=0;pos<len;pos++){
	 for(seg=0;seg<6;seg++){
	  gx[seg]+=db->x[stpos+seg*len+pos];
	  gy[seg]+=db->y[stpos+seg*len+pos];
	  gz[seg]+=db->z[stpos+seg*len+pos];
	 }
	}

	//convert to real cog
	for(seg=0;seg<6;seg++){
	 gx[seg]/=(double)len;
	 gy[seg]/=(double)len;
	 gz[seg]/=(double)len;
	}
	//input ft
	start_tot=total;
	for(p=0;p<15;p++){
	 p1=pat[p][0];
	 p2=pat[p][1];
	 ft[i].f15[p]=(int)(500*sqrt((X(gx[p1]-gx[p2])+X(gy[p1]-gy[p2]) +X(gz[p1]-gz[p2]))));
	}
	 ft[i].id=i;
	 ft[i].pos=0;
	 ft[i].miss=false;
	 ft[i].f=ft[i].f15[4];//sort key******
	//total=start_tot+pos;
  //printf("%d f0= %d %d\n",i,ft[i].f,db->pos[i]);
//	for(p=0;p<15;p++)
//	  printf("%d,",ft[i].f15[p]);
//	puts("");
  pg_bar(i);
 }
 return(i);
}

//A3 mode simple
//Bug fixed -> pos<=num-len
int LowerBonds_A3(char **list,int N,int len,int LenOfSeg,FTABLE *ft,int *res_list){
 //Use only f15[0],[9],[14] data
 int i,id=0,pos,j,k,fpat;
 int pN,pC;//position N,C
 int num=0;
 FILE *fp;
 char line[LIN],buf[9];
 //int *Miss;//Missing residue flag
 double x,y,z,*gx,*gy,*gz;
 double *sx,*sy,*sz;
 //int len=m;
 int total=0,start_tot;

 setbuf(stdout,NULL);

 gx=(double *)malloc(sizeof(double)*RES);
 gy=(double *)malloc(sizeof(double)*RES);
 gz=(double *)malloc(sizeof(double)*RES);

 sx=(double *)malloc(sizeof(double)*RES);
 sy=(double *)malloc(sizeof(double)*RES);
 sz=(double *)malloc(sizeof(double)*RES);

 for(i=0;i<N;i++){
  if((fp=fopen(list[i],"r"))==NULL)
   continue;
  //printf("#Reading..%s\n",list[i]);
  num=0;
 	while(fgets(line,PDB_STRLEN,fp)!=NULL){
		if(!strncmp(&line[13],"CA ",3) && !strncmp(line,"ATOM",4) ){
		 //Ca only
		 strncpy(buf,&line[30],8); buf[8]='\0'; 
		 sx[num]=atof(buf);
		 strncpy(buf,&line[38],8); buf[8]='\0'; 
		 sy[num]=atof(buf);
		 strncpy(buf,&line[46],8); buf[8]='\0'; 
		 sz[num]=atof(buf);
 		 num++;
		}
		//NMR
		if(!strncmp(&line[0],"ENDMDL",6))
		 break;
 	}
  	fclose(fp);
	if(num<len*6){
	 pg_bar(i);
	 continue;
	}
	res_list[i]=num;
	//printf("num=%d\n",num);


	//Liner cal F(S)=|G(S_left)-G(S_rignt)|/2
	//not real cog, only sum of coords
	gx[0]=gy[0]=gz[0]=0;
	for(pos=0;pos<len;pos++){
	 gx[0]+=sx[pos];
	 gy[0]+=sy[pos];
	 gz[0]+=sz[pos];
	}
	for(pos=1;pos<=num-len;pos++){
	 gx[pos]=gx[pos-1]+(sx[pos+len-1]-sx[pos-1]);
	 gy[pos]=gy[pos-1]+(sy[pos+len-1]-sy[pos-1]);
	 gz[pos]=gz[pos-1]+(sz[pos+len-1]-sz[pos-1]);
	}
	//convert to real cog
	for(pos=0;pos<=num-len;pos++){
	 gx[pos]/=(double)len;
	 gy[pos]/=(double)len;
	 gz[pos]/=(double)len;
	}
	//input ft
	start_tot=total;
	//f15[0] 0-1
	/*
        pos
        0     1      2      3      4      5
	|*****++*****++-----++-----++-----++-----|
    	   len    len    len    len    len    len      
	*/
	for(pos=0;pos<=num-2*len;pos++){
	 ft[start_tot+pos].id=i;
	 ft[start_tot+pos].pos=pos;
	 //more simple
	 //convert to int data *1000
	 //(int)(1000*sqrt((X(gx[pos]-gx[pos+len])
	 //+X(gy[pos]-gy[pos+len]) +X(gz[pos]-gz[pos+len]))/2.000));
	 // is
	 //ft[start_tot+pos].f=(int)(500*sqrt((X(gx[pos]-gx[pos+len])
	 //+X(gy[pos]-gy[pos+len]) +X(gz[pos]-gz[pos+len]))));

	 //ft[start_tot+pos].f15[0]=ft[start_tot+pos].f;
	 ft[start_tot+pos].f15[0]=(int)(500*sqrt((X(gx[pos]-gx[pos+len])
	 +X(gy[pos]-gy[pos+len]) +X(gz[pos]-gz[pos+len]))));
	 ft[start_tot+pos].f=ft[start_tot+pos].f15[0];//sort key******
	 ft[start_tot+pos].miss=false;
	}

	/*
        pos
        f            f1            f2
	|-----++-----++-----++-----++-----++-----|
    	   len    len    len    len    len    len      
	*/
	for(pos=0;pos<=num-LenOfSeg;pos++){
	 //New 0-1
	 //ft[start_tot+pos].f15[5] =ft[start_tot+pos+  len].f15[0];//1-2
	 ft[start_tot+pos].f15[9] =ft[start_tot+pos+2*len].f15[0];//2-3
	 //ft[start_tot+pos].f15[12]=ft[start_tot+pos+3*len].f15[0];//3-4
	 ft[start_tot+pos].f15[14]=ft[start_tot+pos+4*len].f15[0];//4-5
	}
	
	//Check Missing and broken residue
	for(k=1;k<num;k++){
	 if(sqrt(X(sx[k]-sx[k-1])+X(sy[k]-sy[k-1])+X(sz[k]-sz[k-1]))>CACADIST)
	  //printf("Broken in id=%d pos=%d\n",i,pos);
	 	//input true data
	 	for(j=1;k-j>=0 && j<LenOfSeg;j++){
		 ft[start_tot+k-j].miss=true;
		}
	}

	total=start_tot+pos;
  pg_bar(i);
 }
 return(total);
}


int write_sdb(char **list,int *res_list,int N,int TotalSubstr,int m,FTABLE *ftbl,char *file){
 int outsize;
 FILE *fo,*fp;
 int i;
 float x,y,z;
 char line[LIN],buf[9];

	if((fo=fopen(file,"wb"))==NULL){
	 puts("ERROR when writing sdb");
	 return FALSE;
	}
	outsize=fwrite(&TotalSubstr,sizeof(int),1,fo);
	outsize=fwrite(&N,sizeof(int),1,fo);
	outsize=fwrite(&m,sizeof(int),1,fo);
	outsize=fwrite(ftbl,sizeof(FTABLE),TotalSubstr,fo);

	//MODEL list....
	int len;
	for(i=0;i<N;i++){
	 len=strlen(list[i]);
	 outsize=fwrite(&i,sizeof(int),1,fo);
	 outsize=fwrite(&len,sizeof(int),1,fo);
	 outsize=fwrite(&res_list[i],sizeof(int),1,fo);
	 outsize=fwrite(list[i],sizeof(char),len,fo);
	 //printf("%d %s %d res= %d\n",i,list[i],len,res_list[i]);
	}
	//Aloso Ca coords!!??
	for(i=0;i<N;i++){
	 if(res_list[i]==0)
	  continue;
	 //set check id
	 outsize=fwrite(&i,sizeof(int),1,fo);
	 if((fp=fopen(list[i],"r"))==NULL)
  	  continue;
 		while(fgets(line,PDB_STRLEN,fp)!=NULL){
		 if(!strncmp(&line[13],"CA ",3) && !strncmp(line,"ATOM",4) ){
		  //Ca only
		  strncpy(buf,&line[30],8); buf[8]='\0'; 
		  x=atof(buf);
	 	  outsize=fwrite(&x,sizeof(float),1,fo);
		  strncpy(buf,&line[38],8); buf[8]='\0'; 
		  y=atof(buf);
	 	  outsize=fwrite(&y,sizeof(float),1,fo);
		  strncpy(buf,&line[46],8); buf[8]='\0'; 
		  z=atof(buf);
	 	  outsize=fwrite(&z,sizeof(float),1,fo);
		  //Amino3char
		  strncpy(buf,&line[17],3);
	 	  outsize=fwrite(buf,sizeof(char),3,fo);
		 }
		 //NMR
		 if(!strncmp(&line[0],"ENDMDL",6))
		  break;
		}
  	 fclose(fp);
 	}

	//for(i=0;i<10;i++)
	// printf("%d %d %d\n",ftbl[i].id,ftbl[i].pos,ftbl[i].f);

	printf("wrote %d substr, %d model\n",TotalSubstr,N);
	fclose(fo);
 return TRUE;
}

extern char **dblist;

int read_sdb(int *res_list,int *N,int *TotalSubstr,int m,FTABLE **ftbl,CAPDB *ca,char *file){
 int outsize;
 FILE *fo;
 int i,M,id;
 //FTABLE *ftbl;

	if((fo=fopen(file,"rb"))==NULL){
	 puts("ERROR when reading sdb");
	 return FALSE;
	}
	outsize=fread(TotalSubstr,sizeof(int),1,fo);
	outsize=fread(N,sizeof(int),1,fo);
	outsize=fread(&M,sizeof(int),1,fo);
	if(m != M){
	 printf("wrong substr length... %d(-n option) != %d(database)\n",m,M);
	 return FALSE;
	}
	//malloc
	if((res_list=(int *)malloc(sizeof(int)*(*N)))==NULL)
	 return FALSE;

	if((*ftbl=(FTABLE *)malloc(sizeof(FTABLE)*(*TotalSubstr)))==NULL)
	 return FALSE;

	if((dblist=(char **)malloc(sizeof(char *)*(*N)))==NULL)
	 return FALSE;
	//position
	if((ca->pos=(int *)malloc(sizeof(int)*(*N)))==NULL)
	 return FALSE;
	//num
	if((ca->n=(int *)calloc((*N),sizeof(int)))==NULL)
	 return FALSE;
	

	outsize=fread(*ftbl,sizeof(FTABLE),*TotalSubstr,fo);
	ca->Nst=*N;
	printf("TotalSubStr=%d Model=%d len=%d\n",*TotalSubstr,*N,M);

	//list....
	int len,tmp=0;
	for(i=0;i<(*N);i++){
	 outsize=fread(&id,sizeof(int),1,fo);
	 if(i!=id){
	  puts("sdb format error");
	  return FALSE;
	 }
	 outsize=fread(&len,sizeof(int),1,fo);
	 outsize=fread(&res_list[i],sizeof(int),1,fo);
	 tmp+=res_list[i];
	 if((dblist[i]=(char *)malloc(sizeof(char)*len+1))==NULL)
	  return FALSE;
	 outsize=fread(dblist[i],sizeof(char),len,fo);
	 dblist[i][len]='\0';

	 //printf("**%d|%s|%d|%d\n",i,list[i],len,res_list[i]);
	}
	printf("TOTAL RES= %d\n",tmp);


	if((ca->x=(float *)malloc(sizeof(float)*(tmp)))==NULL)
	 return FALSE;
	if((ca->y=(float *)malloc(sizeof(float)*(tmp)))==NULL)
	 return FALSE;
	if((ca->z=(float *)malloc(sizeof(float)*(tmp)))==NULL)
	 return FALSE;
/*
	if((ca->cd=(double **)malloc(sizeof(double*)*(tmp)))==NULL)
	 return FALSE;
	for(i=0;i<tmp;i++)
	 if((ca->cd[i]=(double *)malloc(sizeof(double)*(3)))==NULL)
	  return FALSE;
*/
	if((ca->ami=(char *)malloc(sizeof(char)*3*(tmp)))==NULL)
	 return FALSE;

	//Read Ca Coords And Amino3char!!
	float x,y,z;
	int total=0;
	int num;
	char ami[3];
	for(i=0;i<(*N);i++){
	 ca->n[i]=res_list[i];
	 if(res_list[i]==0)
	  continue;
	 //malloc for ca
	 //printf("%d>>%d\n",i,res_list[i]);
	 ca->pos[i]=total;

	 outsize=fread(&id,sizeof(int),1,fo);
	 if(i!=id){
	  puts("sdb format error");
	  return FALSE;
	 }
	 num=0;
		while(num<res_list[i]){
		 outsize=fread(&x,sizeof(float),1,fo);
		 outsize=fread(&y,sizeof(float),1,fo);
		 outsize=fread(&z,sizeof(float),1,fo);
		 outsize=fread(&(ca->ami[total*3]),sizeof(char),3,fo);

		 ca->x[total]=x; ca->y[total]=y; ca->z[total]=z;
		 //ca->cd[total][0]=x;
		 //ca->cd[total][1]=y;
		 //ca->cd[total][2]=z;
		 //printf("%c%c%c\n",ca->ami[total*3],ca->ami[total*3+1],ca->ami[total*3+2]);

		 total++;
		 num++;
		}
	}

	fclose(fo);
	ca->total=total;
	//for(i=0;i<10;i++)
	// printf("%d %d %d\n",ftbl[i]->id,ftbl[i]->pos,ftbl[i]->f);
 return TRUE;
}

int cmp_f_ftbl(const void *c1, const void *c2){
 int a1=*(int *)c1;
 FTABLE a2=*(FTABLE *)c2;
 return(a1-a2.f);
}

int pat15check(FTABLE *q,FTABLE *ftbl,int c){
 int i,p0,p1,p2;
 double dist;
 for(i=14;i>=0;i--){
  p0=pat15[i][0];
  p1=pat15[i][1];
  p2=pat15[i][2];
  dist=DIST(q->f15[p0],ftbl->f15[p0],q->f15[p1],ftbl->f15[p1],q->f15[p2],ftbl->f15[p2]);
  if(dist>(double)c)
   return FALSE;
 }
 return TRUE;
}

int comp_sdb(int m,FTABLE *q,int Nq,FTABLE *ftbl,int Ndb,\
double rms,CAPDB *qca,CAPDB *db,SHIBRESULTS *res){
 int iq,idb;

 double k=1.00;
 FTABLE *a,*fwd,*end;
 double r;
 int num,i,j;
 //COORD *Qcd;
 double Qcd[RES][3];
 int seg_len=(int)((double)m/3.000);
 int m_dash=(int)((double)m/3.000);

 FTABLE **Stock,**Fstock;
 if((Stock=(FTABLE **)malloc(sizeof(FTABLE *)*Ndb))==NULL)
  return FALSE;
 if((Fstock=(FTABLE **)malloc(sizeof(FTABLE *)*Ndb))==NULL)
  return FALSE;

 double *RMSstock;
 if((RMSstock=(double *)malloc(sizeof(double)*Ndb))==NULL)
  return FALSE;

 //malloc for results data
 if((res->f=(FTABLE ***)malloc(sizeof(FTABLE **)*Nq))==NULL)
  return FALSE;
 if((res->rms=(double **)malloc(sizeof(double *)*Nq))==NULL)
  return FALSE;
 if((res->n=(int *)malloc(sizeof(int)*Nq))==NULL)
  return FALSE;
 if((res->shib=(SHIBSET **)malloc(sizeof(SHIBSET *)*Nq))==NULL)
  return FALSE;

 //malloc flgs
 int **flgs;
 if(cmd.fmode==true){
 	if((flgs=(int **)malloc(sizeof(int *)*db->Nst))==NULL)
 	 return FALSE;
 	for(i=0;i<db->Nst;i++){
 	 if((flgs[i]=(int *)malloc(db->n[i]*sizeof(int)))==NULL)
 	  return FALSE;
		//init
		for(j=0;j<db->n[i];j++)
		 flgs[i][j]=-1;
	 }
 }

 //odds or even?
 if(m_dash%2==1)
  k=sqrt((double)(m_dash-1)/(double)m_dash);

 //RMSD>=sqrt(m_dash/m)*sqrt{D(P1,Q1)^2+D(P2,Q2)^2+D(P3,Q3)^2}
 //D(P1,Q1)=k*|F(P1)-F(Q1)|
 //sqrt{D(P1,Q1)^2+D(P2,Q2)^2+D(P3,Q3)^2}=DIST(P,Q)
 //RMSD/k*sqrt(m/m_dash)
 //
 int c=(int)(sqrt((double)m/(double)m_dash)*1000*rms/k);

 printf("#k=%.3f c=%.1f*k*sqrt(m=%d/m'=%d)=%d\n",k,rms*1000,m,m_dash,c);

 int Fs,Fs2,Ft;
 int Fs0S,Fs0E,Fs1S,Fs1E,Fs2S,Fs2E;
 int chk=0;
 //int **HASHTBL,**KTBL,hashid;
 int hash_pos;
 puts("#START Compare db");


 // RMSD(S,T)>=sqrt(m'/m)*k*|DIST(S,T)|
 // F(S)-rms/k <= F(T) <= F(S)+rms/k

	for(iq=0;iq<Nq;iq++){//positions!!
	 Fs0S=q[iq].f - c;
	 if(Fs0S<0)
	  Fs0S=ftbl[0].f;
	 //puts("S0");
	 Fs0E=q[iq].f + c;

	 num=0;
	 chk=0;
	 res->n[iq]=0;//init
	 //New check broken
	 if(cmd.IgBroken==true && q[iq].miss==true){
	  printf("Ignore Broken Segments\n");
	  res->n[iq]=0;
	  continue;
	 }
	 	//binary search
	 	while(Fs0S<ftbl[Ndb-1].f){
	 	 a=(FTABLE *)bsearch(&Fs0S,ftbl,Ndb,sizeof(FTABLE),cmp_f_ftbl);
	 	 if(a==NULL)
		  Fs0S++;
	 	 else
		  break;
		}
		//Forward
		fwd=a;
		fwd--;
		while(fwd-ftbl>=0 && fwd->f>=Fs0S){
		 //check f1, f2 ->no need
		 //if(DIST(q[iq].f,fwd->f,q[iq].f1,fwd->f1,q[iq].f2,fwd->f2) <=(double)c ){
		  if(pat15check(&q[iq],fwd,c)==TRUE){
		   Stock[chk]=fwd;
		   chk++;
		  }
		 //}
		 fwd--;
		}
		//Check D(S,T)
		while(a-ftbl<Ndb){
		 //if(DIST(q[iq].f,a->f,q[iq].f1,a->f1,q[iq].f2,a->f2) <=(double)c ){
		  if(pat15check(&q[iq],a,c)==TRUE){
		   Stock[chk]=a;
		   chk++;
		  }
		 //}
		 if(a->f > Fs0E)
		  break;
		 a++;
		}
	 	printf("#Qpos=%4d FsRange=%d-%d-%d ",iq,Fs0S,q[iq].f,Fs0E);
	 	printf("chkRate= %.3f (%d)\n",(double)chk/(double)Ndb,chk);

		if(chk==0)
		 continue;
		//Copy Coords
                for(i=0;i<m;i++){
                 //Qcd[i].x=qca->x[iq+i];
                 Qcd[i][0]=qca->x[iq+i];
                 //Qcd[i].y=qca->y[iq+i];
                 Qcd[i][1]=qca->y[iq+i];
                 //Qcd[i].z=qca->z[iq+i];
                 Qcd[i][2]=qca->z[iq+i];
                }
		//check Real RMSD
		num=0;
		for(i=0;i<chk;i++){
		 //check broken?
		 if(cmd.IgBroken==true && Stock[i]->miss==true){
                  //printf("#BROK %d pos=%d mid=%d\n"\
		 ,iq,Stock[i]->pos+j,Stock[i]->id);
		  //printf("b");
		  continue;
		 }
		 //check flgs and skip
		 if(cmd.fmode==true && flgs[Stock[i]->id][Stock[i]->pos]==iq){
		  //printf("s");
		  continue;
		 }

		 //RMSD HERE!!
		 r=rmsd_capdb(Qcd,db,Stock[i]->pos,Stock[i]->id,m);

		 if(cmd.jmode==true)
		  continue;

		 	//check next positions
			//fast-skip-next mode
		 if(cmd.fmode==true){
			for(j=1;Stock[i]->pos+j<db->n[Stock[i]->id];j++){
			 if(r>rms+CACADIST*j){
                  	  //printf("#FLAG %d pos=%d mid=%d\n"\
		 	  ,iq,Stock[i]->pos+j,Stock[i]->id);
			  flgs[Stock[i]->id][Stock[i]->pos+j]=iq;
			 } else{
			  break;
			 }
			}
			for(j=1;Stock[i]->pos-j>=0;j++){
			 if(r>rms+CACADIST*j){
                  	  //printf("#FLAG*%d pos=%d mid=%d\n"\
		 	  ,iq,Stock[i]->pos-j,Stock[i]->id);
			  flgs[Stock[i]->id][Stock[i]->pos-j]=iq;
			 } else{
			  break;
			 }
			}
		 }
		 if(r<=rms){
			//Zero bug check
			//if(r==0 && DIST(q[iq].f,Stock[i]->f,q[iq].f1,Stock[i]->f1,q[iq].f2,Stock[i]->f2)!=0){
			if(r==0 && DIST(q[iq].f15[4],Stock[i]->f15[4],q[iq].f15[7]\
			,Stock[i]->f15[7],q[iq].f15[9],Stock[i]->f15[9])!=0){
			 continue;
			}
                  //printf("#Hit %d pos=%d mid=%d rmsd= %.3f D(Q,T)= %.3f\n"\
		  ,iq,Stock[i]->pos,Stock[i]->id,r,\
		  DIST(q[iq].f15[4],Stock[i]->f15[4],q[iq].f15[7],Stock[i]->f15[7]\
		  ,q[iq].f15[9],Stock[i]->f15[9]));
		  Fstock[num]=Stock[i];
		  RMSstock[num]=r;
		  num++;
		 }
		}
	 //input data!!
	 res->n[iq]=num;
	 if(num==0)
	  continue;
	 //malloc AREA
	 if((res->f[iq]=(FTABLE **)malloc(sizeof(FTABLE *)*num))==NULL){
	  puts("Memory Over Flow IN SHIBRESULTS");
	  return FALSE;
	 }
	 if((res->rms[iq]=(double *)malloc(sizeof(double)*num))==NULL){
	  puts("Memory Over Flow IN SHIBRESULTS");
	  return FALSE;
	 }
	  if((res->shib[iq]=(SHIBSET *)malloc(sizeof(SHIBSET)*num))==NULL){
	  puts("Memory Over Flow IN SHIBRESULTS");
	  return FALSE;
	 }
	 //END malloc

	 for(i=0;i<num;i++){
	  res->rms[iq][i]=RMSstock[i];
	  res->f[iq][i]=Fstock[i];
	  res->shib[iq][i].f=Fstock[i];
	  res->shib[iq][i].rms=RMSstock[i];
	 }
	 //INPUT END
	}
 return TRUE;
}


int comp_sdb_A3(int m,FTABLE *q,int Nq,FTABLE *ftbl,int Ndb,\
double rms,CAPDB *qca,CAPDB *db,SHIBRESULTS *res){
 int iq,idb;

 double k=1.00;
 FTABLE *a,*fwd,*end;
 double r;
 int num,i,j;
 //COORD *Qcd;
 double Qcd[RES][3];
 int seg_len=(int)((double)m/3.000);
 int m_dash=(int)((double)m/3.000);

 FTABLE **Stock,**Fstock;
 if((Stock=(FTABLE **)malloc(sizeof(FTABLE *)*Ndb))==NULL)
  return FALSE;
 if((Fstock=(FTABLE **)malloc(sizeof(FTABLE *)*Ndb))==NULL)
  return FALSE;

 double *RMSstock;
 if((RMSstock=(double *)malloc(sizeof(double)*Ndb))==NULL)
  return FALSE;

 //malloc for results data
 if((res->f=(FTABLE ***)malloc(sizeof(FTABLE **)*Nq))==NULL)
  return FALSE;
 if((res->rms=(double **)malloc(sizeof(double *)*Nq))==NULL)
  return FALSE;
 if((res->n=(int *)malloc(sizeof(int)*Nq))==NULL)
  return FALSE;
 if((res->shib=(SHIBSET **)malloc(sizeof(SHIBSET *)*Nq))==NULL)
  return FALSE;

 //malloc flgs
 int **flgs;
/*
 if(cmd.fmode==true){
 	if((flgs=(int **)malloc(sizeof(int *)*db->Nst))==NULL)
 	 return FALSE;
 	for(i=0;i<db->Nst;i++){
 	 if((flgs[i]=(int *)malloc(db->n[i]*sizeof(int)))==NULL)
 	  return FALSE;
		//init
		for(j=0;j<db->n[i];j++)
		 flgs[i][j]=-1;
	 }
 }
*/
 //odds or even?
 if(m_dash%2==1)
  k=sqrt((double)(m_dash-1)/(double)m_dash);

 //RMSD>=sqrt(m_dash/m)*sqrt{D(P1,Q1)^2+D(P2,Q2)^2+D(P3,Q3)^2}
 //D(P1,Q1)=k*|F(P1)-F(Q1)|
 //sqrt{D(P1,Q1)^2+D(P2,Q2)^2+D(P3,Q3)^2}=DIST(P,Q)
 //RMSD/k*sqrt(m/m_dash)
 //
 int c=(int)(sqrt((double)m/(double)m_dash)*1000*rms/k);

 printf("#k=%.3f c=%.1f*k*sqrt(m=%d/m'=%d)=%d\n",k,rms*1000,m,m_dash,c);

 int Fs,Fs2,Ft;
 int Fs0S,Fs0E,Fs1S,Fs1E,Fs2S,Fs2E;
 int chk=0;
 int hash_pos;
 puts("#START Compare db");


 // RMSD(S,T)>=sqrt(m'/m)*k*|DIST(S,T)|
 // F(S)-rms/k <= F(T) <= F(S)+rms/k

	for(iq=0;iq<Nq;iq++){//positions!!
	 Fs0S=q[iq].f - c;
	 if(Fs0S<0)
	  Fs0S=ftbl[0].f;
	 //puts("S0");
	 Fs0E=q[iq].f + c;

	 num=0;
	 chk=0;
	 //New check broken
	 if(cmd.IgBroken==true && q[iq].miss==true){
	  printf("Ignore Broken Segments\n");
	  res->n[iq]=0;
	  continue;
	 }
	 	//binary search
	 	while(Fs0S<ftbl[Ndb-1].f){
	 	 a=(FTABLE *)bsearch(&Fs0S,ftbl,Ndb,sizeof(FTABLE),cmp_f_ftbl);
	 	 if(a==NULL)
		  Fs0S++;
	 	 else
		  break;
		}
		//Forward
		fwd=a;
		fwd--;
		while(fwd-ftbl>=0 && fwd->f>=Fs0S){
		  //original A3
		  if(DIST(q[iq].f15[0],fwd->f15[0],q[iq].f15[9],fwd->f15[9],q[iq].f15[14],fwd->f15[14]) <=(double)c ){
		   Stock[chk]=fwd;
		   chk++;
		  }
		 //}
		 fwd--;
		}
		//Check D(S,T)
		while(a-ftbl<Ndb){
		  if(DIST(q[iq].f15[0],a->f15[0],q[iq].f15[9],a->f15[9],q[iq].f15[14],a->f15[14]) <=(double)c ){
		   Stock[chk]=a;
		   chk++;
		  }
		 //}
		 if(a->f > Fs0E)
		  break;
		 a++;
		}
	 	printf("#Qpos=%4d FsRange=%d-%d-%d ",iq,Fs0S,q[iq].f,Fs0E);
	 	printf("chkRate= %.3f (%d)\n",(double)chk/(double)Ndb,chk);

		if(chk==0)
		 continue;
		//Copy Coords
                for(i=0;i<m;i++){
                 Qcd[i][0]=qca->x[iq+i];
                 Qcd[i][1]=qca->y[iq+i];
                 Qcd[i][2]=qca->z[iq+i];
                }
		//check Real RMSD
		num=0;
		for(i=0;i<chk;i++){
		 //check broken?
		 if(cmd.IgBroken==true && Stock[i]->miss==true){
		  continue;
		 }
		 //check flgs and skip
		 if(cmd.fmode==true && flgs[Stock[i]->id][Stock[i]->pos]==iq){
		  continue;
		 }
		 //just calculate check rate
		 if(cmd.jmode==true)
		  continue;

		 //RMSD HERE!!
		 r=rmsd_capdb(Qcd,db,Stock[i]->pos,Stock[i]->id,m);


		 	//check next positions
			//fast-skip-next mode
/*
		 if(cmd.fmode==true){
			for(j=1;Stock[i]->pos+j<db->n[Stock[i]->id];j++){
			 if(r>rms+CACADIST*j){
                  	  //printf("#FLAG %d pos=%d mid=%d\n"\
		 	  ,iq,Stock[i]->pos+j,Stock[i]->id);
			  flgs[Stock[i]->id][Stock[i]->pos+j]=iq;
			 } else{
			  break;
			 }
			}
			for(j=1;Stock[i]->pos-j>=0;j++){
			 if(r>rms+CACADIST*j){
                  	  //printf("#FLAG*%d pos=%d mid=%d\n"\
		 	  ,iq,Stock[i]->pos-j,Stock[i]->id);
			  flgs[Stock[i]->id][Stock[i]->pos-j]=iq;
			 } else{
			  break;
			 }
			}
		 }
*/
		 if(r<=rms){
			//Zero bug check
			if(r==0 && (DIST(q[iq].f15[0],Stock[i]->f15[0],q[iq].f15[9]\
			,Stock[i]->f15[9],q[iq].f15[14],Stock[i]->f15[14])!=0\ 
			|| isnan(DIST(q[iq].f15[0],Stock[i]->f15[0],q[iq].f15[9]\
                        ,Stock[i]->f15[9],q[iq].f15[14],Stock[i]->f15[14])))){
			 continue;
			}
                  printf("#Hit %d pos=%d mid=%d rmsd= %.3f D(Q,T)= %.3f\n"\
		  ,iq,Stock[i]->pos,Stock[i]->id,r,\
		  DIST(q[iq].f15[0],Stock[i]->f15[0],q[iq].f15[9]\
			,Stock[i]->f15[9],q[iq].f15[14],Stock[i]->f15[14]));
		  Fstock[num]=Stock[i];
		  RMSstock[num]=r;
		  num++;
		 }
		}
	 //input data!!
	 res->n[iq]=num;
	 if(num==0)
	  continue;
	 //malloc AREA
	 if((res->f[iq]=(FTABLE **)malloc(sizeof(FTABLE *)*num))==NULL){
	  puts("Memory Over Flow IN SHIBRESULTS");
	  return FALSE;
	 }
	 if((res->rms[iq]=(double *)malloc(sizeof(double)*num))==NULL){
	  puts("Memory Over Flow IN SHIBRESULTS");
	  return FALSE;
	 }
	  if((res->shib[iq]=(SHIBSET *)malloc(sizeof(SHIBSET)*num))==NULL){
	  puts("Memory Over Flow IN SHIBRESULTS");
	  return FALSE;
	 }
	 //END malloc

	 for(i=0;i<num;i++){
	  res->rms[iq][i]=RMSstock[i];
	  res->f[iq][i]=Fstock[i];
	  res->shib[iq][i].f=Fstock[i];
	  res->shib[iq][i].rms=RMSstock[i];
	 }
	 //INPUT END
	}
 return TRUE;
}


int comp_sdb_NV(int m,FTABLE *q,int Nq,FTABLE *ftbl,int Ndb,\
double rms,CAPDB *qca,CAPDB *db,SHIBRESULTS *res){
 int iq,idb;

 double k=1.00;
 FTABLE *a,*fwd,*end;
 double r;
 int num,i,j;
 //COORD *Qcd;
 double Qcd[RES][3];
 int seg_len=(int)((double)m/3.000);
 int m_dash=(int)((double)m/3.000);

 FTABLE **Stock,**Fstock;
 if((Stock=(FTABLE **)malloc(sizeof(FTABLE *)*Ndb))==NULL)
  return FALSE;
 if((Fstock=(FTABLE **)malloc(sizeof(FTABLE *)*Ndb))==NULL)
  return FALSE;

 double *RMSstock;
 if((RMSstock=(double *)malloc(sizeof(double)*Ndb))==NULL)
  return FALSE;

 //malloc for results data
 if((res->f=(FTABLE ***)malloc(sizeof(FTABLE **)*Nq))==NULL)
  return FALSE;
 if((res->rms=(double **)malloc(sizeof(double *)*Nq))==NULL)
  return FALSE;
 if((res->n=(int *)malloc(sizeof(int)*Nq))==NULL)
  return FALSE;
 if((res->shib=(SHIBSET **)malloc(sizeof(SHIBSET *)*Nq))==NULL)
  return FALSE;

 //malloc flgs
 int **flgs;
 //odds or even?
 if(m_dash%2==1)
  k=sqrt((double)(m_dash-1)/(double)m_dash);

 int c=(int)(sqrt((double)m/(double)m_dash)*1000*rms/k);

 printf("#k=%.3f c=%.1f*k*sqrt(m=%d/m'=%d)=%d\n",k,rms*1000,m,m_dash,c);

 int Fs,Fs2,Ft;
 int Fs0S,Fs0E,Fs1S,Fs1E,Fs2S,Fs2E;
 int chk=0;
 int hash_pos;
 puts("#START Compare db");

	for(iq=0;iq<Nq;iq++){//positions!!

	 num=0;
	 chk=0;
	 //New check broken
	 if(cmd.IgBroken==true && q[iq].miss==true){
	  printf("Ignore Broken Segments\n");
	  res->n[iq]=0;
	  continue;
	 }
	 	printf("#Qpos=%4d FsRange=%d-%d-%d ",iq,Fs0S,q[iq].f,Fs0E);
	 	printf("chkRate= %.3f (%d)\n",(double)chk/(double)Ndb,chk);

		//Copy Coords
                for(i=0;i<m;i++){
                 Qcd[i][0]=qca->x[iq+i];
                 Qcd[i][1]=qca->y[iq+i];
                 Qcd[i][2]=qca->z[iq+i];
                }
		//check Real RMSD
		num=0;
		chk=Ndb;
		for(i=0;i<chk;i++){
		 //check broken?
		 Stock[i]=&ftbl[i];
		 if(cmd.IgBroken==true && Stock[i]->miss==true){
		  continue;
		 }
		 //check flgs and skip
		 if(cmd.fmode==true && flgs[Stock[i]->id][Stock[i]->pos]==iq){
		  continue;
		 }
		 if(cmd.jmode==true)
		  continue;

		 //RMSD HERE!!
		 r=rmsd_capdb(Qcd,db,Stock[i]->pos,Stock[i]->id,m);


		 	//check next positions
			//fast-skip-next mode
		 if(r<=rms){
			//Zero bug check
			if(r==0 && DIST(q[iq].f15[0],Stock[i]->f15[0],q[iq].f15[9]\
			,Stock[i]->f15[9],q[iq].f15[14],Stock[i]->f15[14])!=0){
			 continue;
			}
                  printf("#Hit %d pos=%d mid=%d rmsd= %.3f D(Q,T)= %.3f\n"\
		  ,iq,Stock[i]->pos,Stock[i]->id,r,\
		  DIST(q[iq].f15[0],Stock[i]->f15[0],q[iq].f15[9]\
			,Stock[i]->f15[9],q[iq].f15[14],Stock[i]->f15[14]));
		  Fstock[num]=Stock[i];
		  RMSstock[num]=r;
		  num++;
		 }
		}
	 //input data!!
	 res->n[iq]=num;
	 if(num==0)
	  continue;
	 //malloc AREA
	 if((res->f[iq]=(FTABLE **)malloc(sizeof(FTABLE *)*num))==NULL){
	  puts("Memory Over Flow IN SHIBRESULTS");
	  return FALSE;
	 }
	 if((res->rms[iq]=(double *)malloc(sizeof(double)*num))==NULL){
	  puts("Memory Over Flow IN SHIBRESULTS");
	  return FALSE;
	 }
	  if((res->shib[iq]=(SHIBSET *)malloc(sizeof(SHIBSET)*num))==NULL){
	  puts("Memory Over Flow IN SHIBRESULTS");
	  return FALSE;
	 }
	 //END malloc

	 for(i=0;i<num;i++){
	  res->rms[iq][i]=RMSstock[i];
	  res->f[iq][i]=Fstock[i];
	  res->shib[iq][i].f=Fstock[i];
	  res->shib[iq][i].rms=RMSstock[i];
	 }
	 //INPUT END
	}
 return TRUE;
}

int read_ca(CAPDB *ca,char *file){

 char line[LIN],buf[LIN];
 FILE *fp;
 int num=0;
 //malloc
 ca->pos=(int *)malloc(sizeof(int)*RES);
 ca->ami=(char *)malloc(sizeof(char)*RES*3);
 ca->x=(float *)malloc(sizeof(float)*RES);
 ca->y=(float *)malloc(sizeof(float)*RES);
 ca->z=(float *)malloc(sizeof(float)*RES);

	if((fp=fopen(file,"r"))==NULL)
  	 return FALSE;

	while(fgets(line,PDB_STRLEN,fp)!=NULL){
	 if(!strncmp(&line[13],"CA ",3) && !strncmp(line,"ATOM",4) ){
	  //Ca only
	  strncpy(buf,&line[30],8); buf[8]='\0'; 
	  ca->x[num]=atof(buf);
	  strncpy(buf,&line[38],8); buf[8]='\0'; 
	  ca->y[num]=atof(buf);
	  strncpy(buf,&line[46],8); buf[8]='\0'; 
	  ca->z[num]=atof(buf);
	  //Amino3char
	  strncpy(&(ca->ami[num*3]),&line[17],3);
	  //printf("%f\n",ca->x[num]);
	  num++;
	 }
	 //NMR
	 if(!strncmp(&line[0],"ENDMDL",6))
	  break;
	}
  	 fclose(fp);
 ca->total=num;
 return TRUE;
}

int read_ca_from_list(CAPDB *ca,char **list,int N){

 char line[LIN],buf[LIN];
 FILE *fp;
 int num=0,i;
 char *file;
 setbuf(stdout,NULL);
 for(i=0;i<N;i++){
  //check the Number of Res
  if((fp=fopen(list[i],"r"))==NULL)
   continue;
 	while(fgets(line,PDB_STRLEN,fp)!=NULL){
		if(!strncmp(&line[13],"CA ",3) && !strncmp(line,"ATOM",4) ){
		 //Ca only
 		 num++;
		}
		//NMR
		if(!strncmp(&line[0],"ENDMDL",6))
		 break;
 	}
  	fclose(fp);
	pg_bar(i);
 }

 //malloc
 //ca->pos=(int *)malloc(sizeof(int)*RES);
 //ca->ami=(char *)malloc(sizeof(char)*RES*3);
 ca->x=(float *)malloc(sizeof(float)*num);
 ca->y=(float *)malloc(sizeof(float)*num);
 ca->z=(float *)malloc(sizeof(float)*num);
 ca->pos=(int *)malloc(sizeof(int)*N);
 ca->n=(int *)malloc(sizeof(int)*N);
 ca->Nst=N;
 num=0;//init
 for(i=0;i<N;i++){
  ca->pos[i]=num;
	if((fp=fopen(list[i],"r"))==NULL)
  	 continue;

	while(fgets(line,PDB_STRLEN,fp)!=NULL){
	 if(!strncmp(&line[13],"CA ",3) && !strncmp(line,"ATOM",4) ){
	  //Ca only
	  strncpy(buf,&line[30],8); buf[8]='\0'; 
	  ca->x[num]=atof(buf);
	  strncpy(buf,&line[38],8); buf[8]='\0'; 
	  ca->y[num]=atof(buf);
	  strncpy(buf,&line[46],8); buf[8]='\0'; 
	  ca->z[num]=atof(buf);
	  //Amino3char no need
	  //strncpy(&(ca->ami[num*3]),&line[17],3);
	  //printf("%f\n",ca->x[num]);
	  num++;
	 }
	 //NMR
	 if(!strncmp(&line[0],"ENDMDL",6))
	  break;
	}
  	 fclose(fp);
  ca->n[i]=num-ca->pos[i];
  //printf("Reading...%s pos:%d n:%d\n",list[i],ca->pos[i],ca->n[i]);
 }
 ca->total=num;
 return num;
}


//double rmsd_capdb(COORD *Qcd,CAPDB *t,int tpos,int Nt,int m){
double rmsd_capdb(double Qcd[][3],CAPDB *t,int tpos,int Nt,int m){
 double rmsd=0,rmsd2;
 //COORD dbcd[RES];
 double f[RES][3];
 //double q[RES][3];
 int i;
 int zero_pos=t->pos[Nt]+tpos;
 //MTX mtx;
 //copy db cd
 
///*
//puts("S1");
	for(i=0;i<m;i++){
	 //fast_rmsd
	 f[i][0]=t->x[zero_pos+i];
	 f[i][1]=t->y[zero_pos+i];
	 f[i][2]=t->z[zero_pos+i];
	}
	//i=m-1;
	//for(i=0;i<m;i++)
	// printf("ZERO%d:%.3f %.3f %.3f %.3f %.3f %.3f %f\n"\
	//,i,f[i][0],f[i][1],f[i][2],Qcd[i][0],Qcd[i][1],Qcd[i][2],rmsd);
//*/
	//rms

	//rmsd=fast_lsfit(dbcd,Qcd,m,mtx);
	fast_rmsd(f,Qcd,m,&rmsd);
//puts("S2");
	//printf("RMSD=%f\n",rmsd);
/*
	if(rmsd<2.0){
	 for(i=0;i<m;i++)
	 printf("%.3f %.3f %.3f %.3f %.3f %.3f %f\n"\
	,f[i][0],f[i][1],f[i][2],Qcd[i][0],Qcd[i][1],Qcd[i][2],rmsd);
	}
*/
 return rmsd;
}


int show_shibresults(int Nq,char **list,FTABLE *f,SHIBRESULTS res){
 int i,m,id;
 puts("RESULTS");
 for(i=0;i<Nq;i++){
  printf("POS %5d NumOfHits= %d\n",i,res.n[i]);
  if(res.n[i]==0)
   printf("NoHit\n");

   qsort(res.shib[i],res.n[i],sizeof(SHIBSET),cmp_shib_rms);

	for(m=0;m<res.n[i];m++){
	 //id=res.f[i][m]->id;
	 id=res.shib[i][m].f->id;
  	 //printf("SUBSTR RMSD= %.3f ID %5d POS %5d FILE %s\n",res.rms[i][m],id,res.f[i][m]->pos,list[id]);
  	 //printf("SUBSTR RMSD= %.3f ID %5d POS %5d FILE %s\n",res.shib[i][m].rms,id,res.shib[i][m].f->pos,list[id]);
  	 //printf("SUBSTR RMSD= %.3f ID %5d POS %5d FILE %s\n",res.rms[i][m],id,res.f[i][m]->pos);
	}
  puts("//");
 }
 return TRUE;
}

int set_dmtx(CAPDB *db,FTABLE *f,short int **d, int *cid,int *csize,int N,double r){
 int i,i1,i2;
 int LenOfSeg=db->n[0];
 double rmsd=0;
 int p1,p2;
 double Tcd[RES][3],Qcd[RES][3];

 int m_dash=(int)((double)LenOfSeg/3.000);
 double k=1.00;
 if(m_dash%2==1)
  k=sqrt((double)(m_dash-1)/(double)m_dash);
 int c=(int)(sqrt((double)LenOfSeg/(double)m_dash)*1000*r/k);

 printf("#k= %f c=%d m=%d m'=%d\n",k,c,LenOfSeg,m_dash);
 int igcnt=0;
/*
 //int cid
 for(i1=0;i1<N;i1++){
  cid[i1]=i1;
  csize[i1]=1;
 }
*/
 for(i1=0;i1<N;i1++){
  //copy cd
  p1=db->pos[i1];
  for(i=0;i<db->n[i1];i++){
   Tcd[i][0]=db->x[i+p1];
   Tcd[i][1]=db->y[i+p1];
   Tcd[i][2]=db->z[i+p1];
  }
 	for(i2=i1+1;i2<N;i2++){
	 if(pat15check(&f[i1],&f[i2],c)==FALSE){
	  //printf("Ignore");
	  d[i1][i2]=-1;
	  d[i2][i1]=-1;
	  igcnt++;
	  continue;
	 }
	 //cal rms
	 p2=db->pos[i2];
	 for(i=0;i<db->n[i2];i++){
  	  Qcd[i][0]=db->x[i+p2];
  	  Qcd[i][1]=db->y[i+p2];
  	  Qcd[i][2]=db->z[i+p2];
  	 }
	 fast_rmsd(Tcd,Qcd,db->n[i1],&rmsd);
	 d[i1][i2]=(short int)(rmsd*100.00);
	 d[i2][i1]=(short int)(rmsd*100.00);
	 //printf("%d vs %d = %.3f\n",i1,i2,rmsd);
	}
 }
 printf("#IG cnt= %d rate= %f\n",igcnt,(double)igcnt*2/(double)(N*N-N));
 return 0;
}
