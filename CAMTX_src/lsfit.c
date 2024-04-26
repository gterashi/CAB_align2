#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
//#include <omp.h>
//#include "data.h"
//#include "tbl.h"
#include "struct.h"
#include "func.h"

// char *gtenv(char *ENV);

//extern SW sw;

static void shft(COORD *,COORD *,int ,int );
static void gett(COORD *,COORD *,int ,MTX );
static double cal_sqr(COORD *,int ); 
static void cal_tensor(COORD *,COORD *,int ,MTX );
static void cal_diag1(MTX ,MTX ,MTX ,double ,double );
static void cal_diag2(MTX ,MTX ,double *,int ,int ,int );
static void cal_trns2(MTX ,MTX );
static void rot_cd(COORD *,MTX ,int );
static void st_crd(COORD *,COORD *,int );
static int diagnl(MTX,double *,double *);

double fast_lsfit(COORD *cd0,COORD *cd1,int cnt,MTX trns)
{
	//omp_set_num_threads(2);
	int i;
	double rms,d[MAXMTX];
	COORD g0,g1;

	COORD trans;

	zero_mtx(trns);
	trns[0][0] = trns[1][1] = trns[2][2] = 1.0;
	trns[0][3] = trns[1][3] = trns[2][3] = 0.0;
	//rms = rmsd(cd0,cd1,cnt);
//printf("1)%.3f\n",rms);
	//if(rms < 1.0E-04)
	//	return rms;
//#pragma omp parallel
//#pragma omp sections
//{
//#pragma omp section
	cal_g(cd0,&g0,cnt);
//#pragma omp section
	cal_g(cd1,&g1,cnt);
//}
//#pragma omp parallel
//#pragma omp sections
//{
//#pragma omp section
	shft(cd0,&g0,cnt,1);
//#pragma omp section
	shft(cd1,&g1,cnt,1);
//}
	//rms = rmsd(cd0,cd1,cnt);
//printf("2)%.3f\n",rms);
	//if(rms < 1.0E-04)
	//	return rms;

	gett(cd0,cd1,cnt,trns);
	rot_cd(cd0,trns,cnt);
	rms = rmsd(cd0,cd1,cnt);
	return rms;

	//add trans matrix
	rot_cd(&g0,trns,1);
	diff_cd(&trans,&g1,&g0);//g0->g1
	//rot_cd(&trans,trns,1);

	trns[0][3]=trans.x;
	trns[1][3]=trans.y;
	trns[2][3]=trans.z;
//printf("3)%.3f\n",rms);
/*
	shft(cd1,&g1,cnt,-1);
*/
	return rms;
}


double lsfit(COORD *cd0,COORD *cd1,int cnt,MTX trns)
{
	int i;
	double rms,d[MAXMTX];
	COORD g0,g1;

	COORD trans;

	zero_mtx(trns);
	trns[0][0] = trns[1][1] = trns[2][2] = 1.0;
	trns[0][3] = trns[1][3] = trns[2][3] = 0.0;
	//rms = rmsd(cd0,cd1,cnt);
//printf("1)%.3f\n",rms);
	//if(rms < 1.0E-04)
	//	return rms;
	cal_g(cd0,&g0,cnt);
	cal_g(cd1,&g1,cnt);
	shft(cd0,&g0,cnt,1);
	shft(cd1,&g1,cnt,1);
	rms = rmsd(cd0,cd1,cnt);
//printf("2)%.3f\n",rms);
	if(rms < 1.0E-04)
		return rms;
	gett(cd0,cd1,cnt,trns);
	rot_cd(cd0,trns,cnt);
	rms = rmsd(cd0,cd1,cnt);
	//add trans matrix
	rot_cd(&g0,trns,1);
	diff_cd(&trans,&g1,&g0);//g0->g1
	//rot_cd(&trans,trns,1);
	trns[0][3]=trans.x;
	trns[1][3]=trans.y;
	trns[2][3]=trans.z;
//printf("3)%.3f\n",rms);
/*
	shft(cd1,&g1,cnt,-1);
*/
	return rms;
}

double lsfit2(COORD *cd0,COORD *cd1,int cnt,MTX trns)
{
        int i;
        double rms,d[MAXMTX];
        COORD g0,g1;
        COORD trans;

        zero_mtx(trns);
        trns[0][0] = trns[1][1] = trns[2][2] = 1.0;
	trns[0][3] = trns[1][3] = trns[2][3] = 0.0;
        rms = rmsd(cd0,cd1,cnt);
        if(rms < 1.0E-04)
                return rms;
        //calculate gravity point
        //cal_g(cd0,&g0,cnt); 
        //cal_g(cd1,&g1,cnt);
		//Attention!!! this is changed points!!
        //2 original shift to cd0[0] and cd1[0]
        //g0.x=cd0[0].x; g0.y=cd0[0].y; g0.z=cd0[0].z;
        //g1.x=cd1[0].x; g1.y=cd1[0].y; g1.z=cd1[0].z;

		//no shift!!
        //shft(cd0,&g0,cnt,1); 
        //shft(cd1,&g1,cnt,1);
/*
        rms = rmsd(cd0,cd1,cnt);
        if(rms < 1.0E-04)
                return rms;
*/
        gett(cd0,cd1,cnt,trns);
        rot_cd(cd0,trns,cnt);
        rms = rmsd(cd0,cd1,cnt);
        //add trans matrix
/*
        rot_cd(&g0,trns,1);
        diff_cd(&trans,&g1,&g0);//g0->g1
        trns[0][3]=trans.x; trns[1][3]=trans.y; trns[2][3]=trans.z;
*/
        return rms;
}


void zero_mtx(MTX a)
{
        int i,j;

        for(i = 0; i < MAXMTX - 1; i++){
                for(j = 0; j < MAXMTX - 1; j++)
                        a[i][j] = 0.0;
        }
}


static void shft(COORD *cd,COORD *g,int a_cnt,int vec)
{
	int i;

	for(i = 0; i < a_cnt; i++,cd++){
		cd->x -= (double )vec*g->x;
		cd->y -= (double )vec*g->y;
		cd->z -= (double )vec*g->z;
	}
}

static void gett(COORD *tgt,COORD *ref,int a_cnt,MTX trns)
{
	int i,j,k,l;
	double rms,rmst,rmsr,dett,trs,
		e[MAXMTX],wk[MAXMTX],wsq[MAXMTX],sq[MAXMTX];
	MTX w,s,u,work,tr;



	rms = DBL_MAX;
	rmst = cal_sqr(tgt,a_cnt);
	rmsr = cal_sqr(ref,a_cnt);
	cal_tensor(tgt,ref,a_cnt,w);
	cal_trns2(w,s);
	
/*
	mlt_mtx(tp,p,v,1.0);
*/
	cp_mtx(s,u);
	(void)diagnl(u,e,wk);
//ATTENTION!!!! for triangle only!!
	//if(a_cnt<4 && e[0] < 1.0E-7)
	//	e[0]=1.0E-07;
	for(i = 0; i < MAXMTX - 1; i++){
		if(fabs(e[i]) < 1.0E-07) 
			e[i] = 0.0;
	}
	for(i = 0; i < MAXMTX - 1; i++)
		wsq[i] = sqrt(e[i]);
	for(i = 1; i > -2; i -= 2){
		sq[0] = wsq[0]*(double )i;
		if(fabs(wsq[0]) < 1.0E-07)
			continue;
		for(j = 1; j > -2; j -= 2){
			sq[1] = wsq[1]*(double )j;
			if(fabs(wsq[1]) < 1.0E-07)
				continue;
			for(k = 1; k > -2; k -= 2){
				sq[2] = wsq[2]*(double )k;
				if(fabs(wsq[2]) < 1.0E-07)
					continue;
				cal_diag2(work,u,wsq,i,j,k);
				mlt_mtx(w,work,tr,1.0);
				trs = sq[0] + sq[1] + sq[2];
				cal_diag1(work,tr,w,1.0,trs);
				diagnl(work,e,wk);
				dett = det(tr);
				for(l = 0; l < MAXMTX - 1; l++){
					if(fabs(e[l]) < 1.0E-13)
						e[l] = 0.0;
				}
				if((e[0] >= 0.0) && (e[1] >= 0.0)
				&& (e[2] >= 0.0) && (dett >= 0.0))
					cp_mtx(tr,trns);
			}
		}
	}
}

static double cal_sqr(COORD *cd,int a_cnt)
{
	int i;
	double rms;

	rms = 0.0;
	for(i = 0; i < a_cnt; i++)
		rms += cd->x*cd->x + cd->y*cd->y + cd->z*cd->z;
	return (double )rms/a_cnt;
}

static void cal_tensor(COORD *tgt,COORD *ref,int a_cnt,MTX w)
{
	int i,j;

	for(i = 0; i < MAXMTX - 1; i++){
		for(j = 0; j < MAXMTX - 1; j++)
			w[i][j] = 0.0;
	}
	for(i = 0; i < a_cnt; i++,tgt++,ref++){
		w[0][0] += ref->x*tgt->x;
		w[0][1] += ref->x*tgt->y;
		w[0][2] += ref->x*tgt->z;
		w[1][0] += ref->y*tgt->x;
		w[1][1] += ref->y*tgt->y;
		w[1][2] += ref->y*tgt->z;
		w[2][0] += ref->z*tgt->x;
		w[2][1] += ref->z*tgt->y;
		w[2][2] += ref->z*tgt->z;
	}
	for(i = 0; i < MAXMTX - 1; i++){
		for(j = 0; j < MAXMTX - 1; j++)
			w[i][j] /= (double )a_cnt;
	}
}

static void cal_trns2(MTX p,MTX tp)
{
	int i,j;

	for(i = 0; i < MAXMTX - 1; i++){
		for(j = 0; j <= i; j++){
			tp[i][j] = p[0][i]*p[0][j];
			tp[i][j] += p[1][i]*p[1][j];
			tp[i][j] += p[2][i]*p[2][j];
			if(i != j)
				tp[j][i] = tp[i][j];
		}
	}
}
static void cal_diag1(MTX v,MTX w1,MTX w2,double r,double trs)
{
	int i,j;

	for(i = 0; i < MAXMTX - 1; i++){
		for(j = 0; j <= i; j++){
			v[i][j] = -r*w1[0][i]*w2[0][j];
			v[i][j] += -r*w1[1][i]*w2[1][j];
			v[i][j] += -r*w1[2][i]*w2[2][j];
			if(i != j)
				v[j][i] = v[i][j];
			else
				v[i][j] += trs;
		}
	}
}

static void cal_diag2(MTX v,MTX w,double *sq,int n0,int n1,int n2)
{
	int i,j;

	for(i = 0; i < MAXMTX - 1; i++){
		for(j = 0; j < MAXMTX - 1; j++){
			v[i][j] = w[i][0]*w[j][0]/(*sq)*(double )n0;
			v[i][j] += w[i][1]*w[j][1]/(*(sq+1))*(double )n1;
			v[i][j] += w[i][2]*w[j][2]/(*(sq+2))*(double )n2;
		}
	}
}

static void rot_cd(COORD *cd,MTX trns,int a_cnt)
{
	int i;
	COORD w;

	for(i = 0 ; i < a_cnt; i++,cd++){
		w = *cd;
		mlt_vec(trns,&w,cd);
	}
}

static void st_crd(COORD *cd0,COORD *cd1,int a_cnt)
{
	int i;

	for(i = 0; i < a_cnt; i++,cd0++,cd1++)
		*cd0 = *cd1;
}

double rmsd(COORD *cd1,COORD *cd2,int cnt)
{
        int i;
        double r;

        r = 0.0;
        for(i = 0; i < cnt; i++,cd1++,cd2++){
                r += (cd1->x - cd2->x)*(cd1->x - cd2->x);
                r += (cd1->y - cd2->y)*(cd1->y - cd2->y);
                r += (cd1->z - cd2->z)*(cd1->z - cd2->z);
 	//printf("%5.1f%5.1f%5.1f\n",cd1->x,cd1->y,cd1->z );
 	//printf("%5.1f%5.1f%5.1f\n",cd2->x,cd2->y,cd2->z );
        //printf("--------------\n");
        }
        if(fabs(r) < 1.0E-06)
                return 0.0;
        return sqrt((double )r/cnt);
}
void cal_g(COORD *cd,COORD *g,int cnt)
{
        int i;

        g->x = g->y = g->z = 0.0;
        for(i = 0; i < cnt; i++,cd++){
                g->x += cd->x;
                g->y += cd->y;
                g->z += cd->z;
        }
        g->x /= cnt;
        g->y /= cnt;
        g->z /= cnt;
}

void mlt_vec(MTX a,COORD *pos,COORD *trns)
{
        trns->x = a[0][0]*pos->x + a[0][1]*pos->y + a[0][2]*pos->z;
        trns->y = a[1][0]*pos->x + a[1][1]*pos->y + a[1][2]*pos->z;
        trns->z = a[2][0]*pos->x + a[2][1]*pos->y + a[2][2]*pos->z;
}
void cp_mtx(MTX a,MTX b)
{
        int i,j;

        for(i = 0; i < MAXMTX; i++){
                for(j = 0; j < MAXMTX; j++)
                        b[i][j] = a[i][j];
        }
}

int diagnl(MTX z,double d[],double e[])
{
        double f,x,g,h,hh,p,r,s;
        int i,j,k,l,m,ip1;

        for(i = MAXMTX - 2; i > 0; i--){
                l = i-1;
                f = z[i][l];
                g = 0.0;
                for(k = 0; k < l; k++)
                        g += z[i][k]*z[i][k];
                h = g + f*f;
                if(g <= EPSIRON2){
                        e[i] = f;
                        h=0.0;
                }
                else{
                        if(h > 1.0E-06)
                                g = sqrt(h);
                        else
                                g = 1.0E-06;
                        if(f >= 0.0)
                                g = -g;
                        e[i] = g;
                        h -= f*g;
                        z[i][l] = f - g;
                        f = 0.0;
                        for(j = 0; j < i; j++){
                                z[j][i] = z[i][j]/h;
                                g = 0.0;
                                for(k = 0; k < j; k++)
                                        g += z[j][k]*z[i][k];
                                for(k = j; k < i; k++)
                                        g += z[k][j]*z[i][k];
                                e[j] = g/h;
                                f += z[j][i]*g;
                        }
                        hh = f/(2.0*h);
                        for(j = 0; j < i; j++){
                                f = z[i][j];
                                e[j] -= hh*f;
                                g = e[j];
                                for(k = 0; k <= j; k++)
                                        z[j][k] -= e[k]*f + z[i][k]*g;
                        }
                }
                d[i] = h;
        }
        d[0] = e[0] = 0.0;
        for(i = 0; i < MAXMTX - 1; i++){
                if(d[i] != 0.0){
                        for(j = 0; j < i; j++){
                                g = 0.0;
                                for(k = 0; k < i; k++)
                                g += z[i][k]*z[k][j];
                                for(k = 0; k < i; k++)
                                        z[k][j] -= z[k][i]*g;
                        }
                }
                d[i] = z[i][i];
                z[i][i] = 1.0;
                for(j = 0; j < i; j++)
                        z[i][j] = z[j][i] = 0.0;
        }
        for(i = 1; i < MAXMTX - 1; i++)
                e[i-1] = e[i];
        f = x = e[MAXMTX-2] = 0.0;
        for(l = 0; l < MAXMTX - 1; l++){
                j = 0;
                h = (fabs(d[l]) + fabs(e[l]))*EPSIRON;
                if(x < h)
                        x = h;
                for(m = l; m < MAXMTX - 1; m++){
                        if(fabs(e[m]) <= x)
                                break;
                }
                if(m != l){
                        do{
                                if(j++ >= NTMAX)
                                        return -1;
                                p = (d[l+1]-d[l])/(e[l]+e[l]);
                                if(p*p + 1.0 > 1.0E-06)
                                        r = sqrt(p*p+1.0);
                                else
                                        r = 0.0;
                                if(p < 0.0)
                                        h = d[l]-e[l]/(p-r);
                                else
                                        h = d[l]-e[l]/(p+r);
                                for(i = l; i < MAXMTX - 1; i++)
                                        d[i] -= h;
                                f += h;
                                p = d[m];
                                hh = 1.0;
                                s = 0.0;
                                for(ip1 = m; ip1 > l; ip1--){
                                        i = ip1-1;
                                        g = e[i]*hh;
                                        h = hh*p;
                                        if(fabs(p) >= fabs(e[i])){
                                                hh = e[i]/p;
                                                if(hh*hh + 1.0 > 1.0E-06)
                                                        r = sqrt(hh*hh+1.0);
                                                else
                                                        r = 0.0;
                                                e[ip1] = s*p*r;
                                                s = hh/r;
                                                hh = 1.0/r;
                                        }
                                        else{
                                                hh = p/e[i];
                                                if(hh*hh + 1.0 > 1.0E-06)
                                                        r = sqrt(hh*hh+1.0);
                                                else
                                                        r = 0.0;
                                                e[ip1] = e[i]*s*r;
                                                s = 1.0/r;
                                                hh /= r;
                                        }
                                        p = d[i]*hh - s*g;
                                        d[ip1] = (d[i]*s+hh*g)*s+h;
                                        for(k = 0; k < MAXMTX - 1; k++){
                                                h = z[k][ip1];
                                                z[k][ip1] = z[k][i]*s + hh*h;
                                                z[k][i] = z[k][i]*hh - s*h;
                                        }
                                }
                                e[l] = s*p;
                                d[l] = hh*p;
                        } while(fabs(e[l]) > x);
                }
                d[l] += f;
        }
        for(i = 0; i < MAXMTX - 1; i++){
                k = i;
                p = d[i];
                for(j = i+1; j < MAXMTX - 1; j++){
                        if(d[j] < p){
                                k = j;
                                p = d[j];
                        }
                }
                if( k != i){
                        d[k] = d[i];
                        d[i] = p;
                        for(j = 0; j < MAXMTX - 1; j++){
                                p = z[j][i];
                                z[j][i] = z[j][k];
                                z[j][k] = p;
                        }
               }
        }

        return 0;
}

void mlt_mtx(MTX a,MTX b,MTX c,double r)
{
        int i,j,k;
        double w;

        for(i = 0; i < MAXMTX-1; i++){
                for(j = 0; j < MAXMTX-1; j++){
                        w = 0.0;
                        for(k = 0; k < MAXMTX-1; k++)
                                w += a[i][k]*b[k][j];
                        c[i][j] = w*r;
                }
        }
}

double det(MTX a)
{
        double r;

        r =  a[0][0]*a[1][1]*a[2][2];
        r += a[0][1]*a[1][2]*a[2][0];
        r += a[0][2]*a[2][1]*a[1][0];
        r -= a[0][1]*a[1][0]*a[2][2];
        r -= a[0][2]*a[1][1]*a[2][0];
        r -= a[0][0]*a[1][2]*a[2][1];
        return r;
}

