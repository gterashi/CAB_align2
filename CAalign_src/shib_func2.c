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
#define HASH_FUNC(a,b) a%b
#define DIST(a,b,c,d,e,f) sqrt((a-b)*(a-b)+(c-d)*(c-d)+(e-f)*(e-f))


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

//For A3 Algorithm!!
int LowerBonds(char **list,int N,int len,FTABLE *ft,int *res_list){
 
 int i,id=0,pos;
 int pN,pC;//position N,C
 int num=0;
 FILE *fp;
 char line[LIN],buf[9];
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
 	}
  	fclose(fp);
	if(num<len*6){
	 pg_bar(i);
	 continue;
	}
	res_list[i]=num;


	//Liner cal F(S)=|G(S_left)-G(S_rignt)|
	gx[0]=gy[0]=gz[0];
	for(pos=0;pos<len;pos++){
	 gx[0]+=sx[pos];
	 gy[0]+=sy[pos];
	 gz[0]+=sz[pos];
	}
	for(pos=1;pos<num-len;pos++){
	 gx[pos]=gx[pos-1]+(sx[pos+len-1]-sx[pos-1]);
	 gy[pos]=gy[pos-1]+(sy[pos+len-1]-sy[pos-1]);
	 gz[pos]=gz[pos-1]+(sz[pos+len-1]-sz[pos-1]);
	}

	//input ft
	start_tot=total;
	for(pos=0;pos<num-2*len;pos++){
	 //ft[total].id=i;
	 ft[start_tot+pos].id=i;
	 //ft[total].pos=pos;
	 ft[start_tot+pos].pos=pos;
	 //ft[total].f=(int)50*sqrt(((gx[pos]-gx[pos+len])*(gx[pos]-gx[pos+len])
	 ft[start_tot+pos].f=(int)50*sqrt(((gx[pos]-gx[pos+len])*(gx[pos]-gx[pos+len])
	 +(gy[pos]-gy[pos+len])*(gy[pos]-gy[pos+len])
	 +(gz[pos]-gz[pos+len])*(gz[pos]-gz[pos+len]))/(len*len));
	 //printf("f=%d\n",ft[total].f);
	 //total++;
	}
	/*
        pos
        f            f1            f2
	|-----++-----++-----++-----++-----++-----|
    	   len    len    len    len    len    len      
	*/
	//input f1 and f2
	for(pos=0;pos<num-6*len;pos++){
	 ft[start_tot+pos].f1=ft[start_tot+pos+2*len].f;
	 ft[start_tot+pos].f2=ft[start_tot+pos+4*len].f;
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
	 //printf("%d %s %d res=%d\n",i,list[i],len,res_list[i]);
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
		}
  	 fclose(fp);
 	}

	//for(i=0;i<10;i++)
	// printf("%d %d %d\n",ftbl[i].id,ftbl[i].pos,ftbl[i].f);

	printf("wrote %d substr, %d model\n",TotalSubstr,N);
	fclose(fo);
 return TRUE;
}

int read_sdb(char **list,int *res_list,int *N,int *TotalSubstr,int m,FTABLE **ftbl,CAPDB *ca,char *file){
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

	if((list=(char **)malloc(sizeof(char *)*(*N)))==NULL)
	 return FALSE;

	if((ca->pos=(int *)malloc(sizeof(int)*(*N)))==NULL)
	 return FALSE;
	

	outsize=fread(*ftbl,sizeof(FTABLE),*TotalSubstr,fo);

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
	 list[i]=(char *)malloc(sizeof(char)*len+1);
	 outsize=fread(list[i],sizeof(char),len,fo);
	 list[i][len]='\0';

	 //printf("**%d|%s|%d|%d\n",i,list[i],len,res_list[i]);
	}
	printf("TOTAL RES= %d\n",tmp);


	if((ca->x=(float *)malloc(sizeof(float)*(tmp)))==NULL)
	 return FALSE;
	if((ca->y=(float *)malloc(sizeof(float)*(tmp)))==NULL)
	 return FALSE;
	if((ca->z=(float *)malloc(sizeof(float)*(tmp)))==NULL)
	 return FALSE;
	if((ca->ami=(char *)malloc(sizeof(char)*3*(tmp)))==NULL)
	 return FALSE;

	//Read Ca Coords And Amino3char!!
	float x,y,z;
	int total=0;
	int num;
	char ami[3];
	for(i=0;i<(*N);i++){
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

int comp_sdb(int m,FTABLE *q,int Nq,FTABLE *ftbl,int Ndb,double rms,CAPDB *qca,CAPDB *db){
 int iq,idb;

 double k=1.00;
 FTABLE *a,*fwd,*end;
 double r;
 int num,i;
 //COORD *Qcd;
 double Qcd[RES][3];
 int seg_len=(int)((double)m/3.000);

 FTABLE **Stock;

 if((Stock=(FTABLE **)malloc(sizeof(FTABLE *)*Ndb))==NULL)
  return FALSE;

 //if((Qcd=(double **)malloc(sizeof(double)*m))==NULL)
 // return FALSE;

 //odds or even?
 if(seg_len%2==1)
  k=sqrt((double)(seg_len-1)/(double)seg_len);
 int c=(int)(sqrt((double)m/(double)seg_len)*100*rms/k);
 printf("#k=%.3f c=%.1f*k*sqrt(%d/%d)=%d\n",k,rms*100,m,seg_len,c);

 int Fs,Fs2,Ft;
 int Fs0S,Fs0E,Fs1S,Fs1E,Fs2S,Fs2E;
 int chk=0;
 int **HASHTBL,**KTBL,hashid;
 int hash_pos;
 puts("#START Compare db");


 // RMSD(S,T)>=sqrt(m'/m)*k*|DIST(S,T)|
 // F(S)-rms/k <= F(T) <= F(S)+rms/k

	for(iq=0;iq<Nq;iq++){//positions!!
	 Fs0S=q[iq].f - c;
	 if(Fs0S<0)
	  Fs0S=ftbl[0].f;
	 Fs0E=q[iq].f + c;
	 //Fs1S=q[iq].f1 - c;
	 //Fs1E=q[iq].f1 + c;
	 //Fs2S=q[iq].f2 - c;
	 //Fs2E=q[iq].f2 + c;

	 num=0;
	 chk=0;
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
		while(fwd->f>=Fs0S && fwd-ftbl>=0){
		 
		 //check f1, f2
		 //if(Fs1S<fwd->f1 && Fs1E > fwd->f1 && Fs2S<fwd->f2 && Fs2E > fwd->f1){
		 if(DIST(q[iq].f,fwd->f,q[iq].f1,fwd->f1,q[iq].f2,fwd->f2) <=(double)c ){
		 //printf("%.1f <=C %d\n",DIST(q[iq].f,fwd->f,q[iq].f1,fwd->f1,q[iq].f2,fwd->f2),c);
		 //printf("FWD%d:%d)%d %d %d vs %d %d %d\n",fwd->id,fwd->pos,q[iq].f,q[iq].f1,q[iq].f2,fwd->f,fwd->f1,fwd->f2);
		  Stock[chk]=fwd;
		  chk++;
		 }
		 fwd--;
		}
		//Check D(S,T)
		while(a-ftbl<Ndb){
		 if(DIST(q[iq].f,a->f,q[iq].f1,a->f1,q[iq].f2,a->f2) <=c ){
		  //printf("%.1f <=C %d\n",DIST(q[iq].f,a->f,q[iq].f1,a->f1,q[iq].f2,a->f2),c);
		  //printf("AFT%d:%d)%d %d %d vs %d %d %d\n",a->id,a->pos,q[iq].f,q[iq].f1,q[iq].f2,a->f,a->f1,a->f2);
		  Stock[chk]=a;
		  chk++;
		 }
		 if(a->f > Fs0E)
		  break;
		 a++;
		}
	 	printf("#Qpos=%4d FsRange=%d %d-%d-%d ",iq,Fs0S,q[iq].f,Fs0E);
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
		//check RMSD
		num=0;
		for(i=0;i<chk;i++){
		 
		 //printf("#CHK %7d%5d\n",Stock[i]->id,Stock[i]->pos);
		 r=rmsd_capdb(Qcd,db,Stock[i]->pos,Stock[i]->id,m);
		 //printf("CHK%d:%d)%d %d %d vs %d %d %d\n",Stock[i]->id,Stock[i]->pos,q[iq].f,q[iq].f1,q[iq].f2,Stock[i]->f,Stock[i]->f1,Stock[i]->f2);
                  //printf("#%d(pos=%d mid=%d) rmsd= %.3f diff=%d\n",iq,Stock[i]->pos,Stock[i]->id,r,Stock[i]->f-q[iq].f);
		 if(r<=rms){
                  printf("#Hit %d pos=%d mid=%d) rmsd= %.3f D(Q,T)=%d\n",iq,Stock[i]->pos,Stock[i]->id,r,Stock[i]->f-q[iq].f);
		  num++;
		 }
		}
	 	//printf("Find %d\n",num);
	}

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
	}
  	 fclose(fp);
 ca->total=num;
 return TRUE;
}

//double rmsd_capdb(COORD *Qcd,CAPDB *t,int tpos,int Nt,int m){
double rmsd_capdb(double Qcd[][3],CAPDB *t,int tpos,int Nt,int m){
 double rmsd=0,rmsd2;
 //COORD dbcd[RES];
 double f[RES][3];
 //double q[RES][3];
 int i;
 int zero_pos=t->pos[Nt]+tpos;
 MTX mtx;
 //copy db cd
 
	for(i=0;i<m;i++){
	 //dbcd[i].x=t->x[zero_pos+i];
	 //dbcd[i].y=t->y[zero_pos+i];
	 //dbcd[i].z=t->z[zero_pos+i];
	 //print_cd(&dbcd[i]);puts("**");
	 //print_cd(&Qcd[i]);puts("**11");
	 //fast_rmsd
	 f[i][0]=t->x[zero_pos+i];
	 f[i][1]=t->y[zero_pos+i];
	 f[i][2]=t->z[zero_pos+i];
	 //q[i][0]=Qcd[i].x;
	 //q[i][1]=Qcd[i].y;
	 //q[i][2]=Qcd[i].z;
	}
	//rms

	//rmsd=fast_lsfit(dbcd,Qcd,m,mtx);
	fast_rmsd(f,Qcd,m,&rmsd);
	//printf("RMSD=%f %f\n",rmsd,rmsd2);
 return rmsd;
}
