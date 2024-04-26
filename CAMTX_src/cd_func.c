#include <stdio.h> 
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "struct.h"
#include "func.h"

double dis_cd(COORD *cd1,COORD *cd2){
	        double len; 
		len=L(cd1->x,cd1->y,cd1->z,cd2->x,cd2->y,cd2->z); 
		return len;
}

void cl_cd(COORD *cd){ 
	cd->x=0.000000; 
	cd->y=0.000000; 
	cd->z=0.000000;
}


double in_pro(double x1,double y1,double z1,double x2,double y2,double z2) { 
	double ans; 
	ans=acos((x1*x2+y1*y2+z1*z2)/(L(x1,y1,z1,0,0,0)*L(x2,y2,z2,0,0,0))); 
	return(ans);
}

double in_pro_cd(COORD *cd1,COORD *cd2) { 
	double ans; 
	COORD cd0; 
	cl_cd(&cd0); 
	ans=acos(((cd1->x*cd2->x)+(cd1->y*cd2->y)+(cd1->z*cd2->z))/((dis_cd(cd1,&cd0))*(dis_cd(cd2,&cd0)))); 
	return(ans);
}

double inpro(COORD *cd1,COORD *cd2){
	double ans;
	ans=(cd1->x*cd2->x)+(cd1->y*cd2->y)+(cd1->z*cd2->z);
	return(ans);
}

void COORDtrn(COORD *af,COORD *bf,MTX mtx){
 af->x=mtx[0][0]*bf->x+mtx[0][1]*bf->y+mtx[0][2]*bf->z+mtx[0][3];
 af->y=mtx[1][0]*bf->x+mtx[1][1]*bf->y+mtx[1][2]*bf->z+mtx[1][3];
 af->z=mtx[2][0]*bf->x+mtx[2][1]*bf->y+mtx[2][2]*bf->z+mtx[2][3];
}

void cp_cd(COORD *cd,double x,double y,double z){ 
	cd->x=x; 
	cd->y=y; 
	cd->z=z;
}

void copy_cd(COORD *cd1,COORD *cd2){ 
	cd1->x=cd2->x; 
	cd1->y=cd2->y; 
	cd1->z=cd2->z;
}
void add_cd(COORD *cd1,COORD *cd2,COORD *cd3) { 
	cd1->x=cd2->x+cd3->x; 
	cd1->y=cd2->y+cd3->y; 
	cd1->z=cd2->z+cd3->z;
}

//cd3->cd2
void diff_cd(COORD *cd1,COORD *cd2,COORD *cd3) { 
	cd1->x=cd2->x-cd3->x; 
	cd1->y=cd2->y-cd3->y; 
	cd1->z=cd2->z-cd3->z;
}

void se_cd(COORD *cd1,COORD *cd2,COORD *cd3) { 
	cd1->x=(cd2->x+cd3->x)/2.000000; 
	cd1->y=(cd2->y+cd3->y)/2.000000; 
	cd1->z=(cd2->z+cd3->z)/2.000000;
}

void unit_vec(COORD *cd){
 	double len;
	COORD init;
        cl_cd(&init);
        len=dis_cd(cd,&init);
	cd->x=cd->x/len;
	cd->y=cd->y/len;
	cd->z=cd->z/len;
}

double len_cd(COORD *cd){
	double len;
 	COORD init;
	cl_cd(&init);
	len=dis_cd(cd,&init);
	return len;
}

void cog_cd(COORD *cd,COORD *cd1,COORD *cd2,COORD *cd3) { 
	cd->x=((cd1->x)+(cd2->x)+(cd3->x))/3.000000; 
	cd->y=((cd1->y)+(cd2->y)+(cd3->y))/3.000000; 
	cd->z=((cd1->z)+(cd2->z)+(cd3->z))/3.000000;
}


void print_cd(COORD *cd){ 
	printf("%8.3f%8.3f%8.3f",cd->x,cd->y,cd->z);
}

//(a~b).x = a.y*b.z - a.z*b.y
////(a~b).y = a.z*b.x - a.x*b.z
////(a~b).z = a.x*b.y - a.y*b.x

void out_pro_cd(COORD *cd0,COORD *cd1,COORD *cd2){ 
	cd0->x=(cd1->y)*(cd2->z)-(cd2->y)*(cd1->z); 
	cd0->y=(cd1->z)*(cd2->x)-(cd2->z)*(cd1->x); 
	cd0->z=(cd1->x)*(cd2->y)-(cd2->x)*(cd1->y);
}
void new_cd_get(COORD *new,COORD *cdX,COORD *cdY,COORD *cdZ,COORD *cog,COORD *old){ 
	COORD cd0; 
	diff_cd(&cd0,old,cog); 
	new->x=dis_cd(cog,old)*(cos(in_pro_cd(cdX,&cd0))); 
	new->y=dis_cd(cog,old)*(cos(in_pro_cd(cdY,&cd0))); 
	new->z=dis_cd(cog,old)*(cos(in_pro_cd(cdZ,&cd0)));
}

void times_cd(COORD *new,COORD *old,double time){
	new->x=old->x*time;
	new->y=old->y*time;
	new->z=old->z*time;
}

int misschk_cd(COORD *cd){
 	if(cd->x == -999.999 && cd->y == -999.999 &&cd->z == -999.999)
	 return(-1);
 return(0);
}

//»°³Ñ·Á¤ÎÌÌÀÑ
double vec2area(COORD *cd1,COORD *cd2){
 double area=0;
 double ip=0;
 double len1,len2;

 len1=len_cd(cd1);
 len2=len_cd(cd2);
 ip=inpro(cd1,cd2);
 area=0.500*sqrt(X(len1)*X(len2)-X(ip));
 return(area);
}
//double tortion(COORD *,COODR *,COORD *,COORD *);
double tortion(COORD *cd1,COORD *cd2,COORD *cd3,COORD *cd4){
 COORD cd12,cd21,cd23,cd32,cd34,cd43;
 COORD out1,out2;
 double ang=0;

 //check
 if(misschk_cd(cd1)==-1
  ||misschk_cd(cd2)==-1
  ||misschk_cd(cd3)==-1
  ||misschk_cd(cd4)==-1)
  return(-999.99);
  //return(360.000);

 diff_cd(&cd12,cd2,cd1);
 diff_cd(&cd21,cd1,cd2);
 diff_cd(&cd23,cd3,cd2);
 diff_cd(&cd32,cd2,cd3);
 diff_cd(&cd34,cd4,cd3);
 diff_cd(&cd43,cd3,cd4);

 out_pro_cd(&out1,&cd23,&cd21);
 out_pro_cd(&out2,&cd34,&cd32);
 ang=in_pro_cd(&out1,&out2);
 if(inpro(&cd12,&out2)<0)
  ang*=-1.000;
 return(ANG(ang));
}

//dual vector space
void DVS(COORD *cd1,COORD *cd2){
 double k=0;
 COORD tmp[3];

 //e^1 = (e_2 x e_3)/(e_1 * (e_2 x e_3))
 out_pro_cd(&tmp[0],&cd2[1],&cd2[2]);
 k=inpro(&cd2[0],&tmp[0]);
 times_cd(&cd1[0],&tmp[0],1.00/k);
 
 //e^2 = (e_3 x e_1)/(e_2 * (e_3 x e_1))
 out_pro_cd(&tmp[0],&cd2[2],&cd2[0]);
 k=inpro(&cd2[1],&tmp[0]);
 times_cd(&cd1[1],&tmp[0],1.00/k);

 //e^3 = (e_1 x e_2)/(e_3 * (e_1 x e_2))
 out_pro_cd(&tmp[0],&cd2[0],&cd2[1]);
 k=inpro(&cd2[2],&tmp[0]);
 times_cd(&cd1[2],&tmp[0],1.00/k);
}

void vec_separation(COORD *k,COORD *v,COORD *p){
 COORD dual[3];

 DVS(dual,p);

  k->x=inpro(v,&dual[0]);
  k->y=inpro(v,&dual[1]);
  k->z=inpro(v,&dual[2]);

}

double ave(double *d,UINT n){
 double sum=0,av=0;
 int i;
 for(i=0;i<n;i++)
  sum+=d[i];
 av=sum/(double)n;
 return(av);
}

double std(double *d,UINT n){
 double a,sum=0,std;
 int i;
 a=ave(d,n);
 for(i=0;i<n;i++)
  sum+=X(d[i]-a);

 std=sqrt(sum/(double)n);
 return(std);
}

//histogram -> AVE
double hist_ave(int *d,int n){
 double sum=0,av=0;
 int i,num=0;
 for(i=0;i<n;i++){
  sum+=((double)i+0.500)*(double)d[i];
  num+=d[i];
 }
 av=sum/(double)num;
 return(av);
}
//histogram -> STD
double hist_std(int *d,int n){
 double a,sum=0,std;
 int i,num=0;
 a=hist_ave(d,n);
 //printf("ace:%f\n",a);
 for(i=0;i<n;i++){
  if(d[i]==0)
   continue;
  sum+=X(((double)i+0.500)-a)*(double)d[i];
  num+=d[i];
 }
 std=sqrt(sum/(double)num);
 return(std);
}


void show_node(NODE **n,int i1){
//show node
 int i2=0;
	 //printf("#N%4d%3d",i1,n[i1]->aa);
	 printf("%3d",n[i1]->aa);
  	 printf("%4d",n[i1]->phi);
  	 printf("%4d",n[i1]->psi);
  	 //chi1-5
	 for(i2=0;i2<5;i2++)
	  printf("%4d",n[i1]->chi[i2]);
  	 //printf("\n");
}

void freerotate(COORD *, double r ,MTX);
//2 point free Rotation
void FreeRotate(COORD *cd1,COORD *cd2,COORD *bf,COORD *af,int n,double r){
 int i1,i2,i3;
 COORD tan,shift,tmp[n];
 MTX mtx;

 diff_cd(&tan,cd1,cd2);
 unit_vec(&tan);
 freerotate(&tan,r,mtx);

 times_cd(&shift,cd2,-1);

 //printf("%.2f,%.2f,%.2f\n",mtx[0][0],mtx[0][1],mtx[0][2]);

 //shift
 for(i1=0;i1<n;i1++){
  //puts("0->");print_cd(&bf[i1]);puts("");
  add_cd(&tmp[i1],&shift,&bf[i1]);  
  //puts("1->");print_cd(&tmp[i1]);puts("");
  COORDtrn(&af[i1],&tmp[i1],mtx);
  //puts("2->");print_cd(&af[i1]);puts("");
  add_cd(&af[i1],cd2,&af[i1]);
  //puts("3->");print_cd(&af[i1]);puts("");
 }

}

void freerotate(COORD *cd,double  r ,MTX mtx)
{
    float v[16];
    double w =(double)cos(RAD(r)/2.0);
    double w2 = w * w;
    double s = (double)sin( RAD(r)/2.0);
    double x = cd->x * s;
    double y = cd->y * s;
    double z = cd->z * s;
/*
    float x = n[0] * s;
    float y = n[1] * s;
    float z = n[2] * s;
*/
    double x2 = X(x);
    double y2 = X(y);
    double z2 = X(z);
/*
0 1 2 3
4 5 6 7
8 9 10 11
*/
    v[0] = w2 + x2 - y2 - z2;
    mtx[0][0] = w2 + x2 - y2 - z2;
    v[4] = 2 * ( ( x * y ) - ( w * z ) );
    mtx[1][0] = 2 * ( ( x * y ) - ( w * z ) );
    v[8] = 2 * ( ( x * z ) + ( w * y ) );
    mtx[2][0] = 2 * ( ( x * z ) + ( w * y ) );
    v[12] = 0.0f;

    v[1] = 2 * ( ( x * y ) + ( w * z ) );
    mtx[0][1] = 2 * ( ( x * y ) + ( w * z ) );
    v[5] = w2 - x2 + y2 - z2;
    mtx[1][1] = w2 - x2 + y2 - z2;
    v[9] = 2 * ( ( y * z ) - ( w * x ) );
    mtx[2][1] = 2 * ( ( y * z ) - ( w * x ) );
    v[13] = 0.0f;

    v[2] = 2 * ( ( x * z ) - ( w * y ) );
    mtx[0][2] = 2 * ( ( x * z ) - ( w * y ) );
    v[6] = 2 * ( ( y * z ) + ( w * x ) );
    mtx[1][2] = 2 * ( ( y * z ) + ( w * x ) );
    v[10] = w2 - x2 - y2 + z2;
    mtx[2][2] = w2 - x2 - y2 + z2;
    v[14] = 0.0f;

    mtx[0][3]=0;
    mtx[1][3]=0;
    mtx[2][3]=0;

    v[3] = 0.0f;
    v[7] = 0.0f;
    v[11] = 0.0f;
    v[15] = 1.0f;
}

/*
 *
 af->x=mtx[0][0]*bf->x+mtx[0][1]*bf->y+mtx[0][2]*bf->z+mtx[0][3];
 af->y=mtx[1][0]*bf->x+mtx[1][1]*bf->y+mtx[1][2]*bf->z+mtx[1][3];
 af->z=mtx[2][0]*bf->x+mtx[2][1]*bf->y+mtx[2][2]*bf->z+mtx[2][3];
*/
//ZYZ
void Euler2mtx_ZYZ(double a,double b, double c,MTX m){

 m[0][0]= cos(a)*cos(c)-sin(a)*cos(b)*sin(c);
 m[0][1]=-cos(a)*sin(c)-sin(a)*cos(b)*cos(c);
 m[0][2]= sin(a)*sin(b);

 m[1][0]= sin(a)*cos(c)+cos(a)*cos(b)*sin(c);
 m[1][1]=-sin(a)*sin(c)+cos(a)*cos(b)*cos(c);
 m[1][2]=-cos(a)*sin(b);

 m[2][0]= sin(b)*sin(c);
 m[2][1]= sin(b)*cos(c);
 m[2][2]= cos(b);

 m[0][3]=m[1][3]=m[2][3]=0.000;

}
//X->Y->Z
void Euler2mtx(double a,double b, double c,MTX m){

 m[0][0]= cos(b)*cos(c);
 m[0][1]= sin(a)*sin(b)*cos(c)-cos(a)*sin(c);
 m[0][2]= cos(a)*sin(b)*cos(c)+sin(a)*sin(c);

 m[1][0]= cos(b)*sin(c);
 m[1][1]= cos(a)*cos(c)+sin(a)*sin(b)*sin(c);
 m[1][2]= cos(a)*sin(b)*sin(c)-sin(a)*cos(c);

 m[2][0]=-sin(b);
 m[2][1]= sin(a)*cos(b);
 m[2][2]= cos(a)*cos(b);

 m[0][3]=m[1][3]=m[2][3]=0.000;

}



void show_mtx(MTX m){
 printf("     X2 = (%9.6f)*X1 + (%9.6f)*Y1 + (%9.6f)*Z1 + (%12.6f)\n",m[0][0],m[0][1],m[0][2],m[0][3]);
 printf("     Y2 = (%9.6f)*X1 + (%9.6f)*Y1 + (%9.6f)*Z1 + (%12.6f)\n",m[1][0],m[1][1],m[1][2],m[1][3]);
 printf("     Z2 = (%9.6f)*X1 + (%9.6f)*Y1 + (%9.6f)*Z1 + (%12.6f)\n",m[2][0],m[2][1],m[2][2],m[2][3]);
}

double Random(){//0-1 random 
        double r;
        int i;
        i=random();
        if(i != 0)
                i--;
        //r=(double)i/0x7fff; 
        r=(double)i/RAND_MAX;
        return(r);
}

int range_rand_i(double min,double max){ //min-max random
        double Random();
        return((max-min)*Random()+min);
}

double range_rand(double min,double max){ //min-max random
        double Random();
        return((max-min)*Random()+min);
}

int flip(double probability){ //Íð¿ô¤Ë¤è¤ë³ÎÎ¨È½Äê
        double Random();
        if(probability == 1.0)
                return(0);
        else if(Random() <= probability)
                return(0);
        else
                return(-1);
}

