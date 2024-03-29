#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>

int compare2(const void *x, const void *y)
{
  double xx = *(double*)x, yy = *(double*)y;
  if (xx < yy) return -1;
  if (xx > yy) return  1;
  return 0;
}


double scaleTau(int n, double*x, double c1, double c2, double correc) {
int i;
double median=0;
double sigma0=1;
double xbalt[n];

for(i = 0; i < n; i++) {
	xbalt[i]=x[i];
	}

qsort(x,n,sizeof(double),compare2);
if (n % 2==0) { 
	median=(x[n/2]+x[n/2-1])/2;
	}
if (n % 2==1) { 
	median=x[(n+1)/2-1];
	}

for(i = 0; i < n; i++) {
	x[i]=xbalt[i];
	}

double xb[n];
for(i = 0; i < n; i++) {
	xb[i]=fabs(x[i]-median);
	}

for(i = 0; i < n; i++) {
	xbalt[i]=xb[i];
	}
qsort(xb,n,sizeof(double),compare2);
if (n % 2==0) { 
	sigma0=(xb[n/2]+xb[n/2-1])/2;
	}
if (n % 2==1) { 
	sigma0=xb[(n+1)/2-1];
	}
double mu=0;
for(i = 0; i < n; i++) {
	xb[i]=xbalt[i];
	}
double w[n];
for(i = 0; i < n; i++) {
	xb[i] = xb[i]/(sigma0*c1);
	w[i]=1-xb[i]*xb[i];
	w[i]=((fabs(w[i])+w[i])/2)*((fabs(w[i])+w[i])/2);
	mu=mu+w[i]*x[i];
	}
double deno=0;
for(i = 0; i < n; i++) {
	deno=deno+w[i];
	}
mu=mu/deno;
double vari=0;
double rho;
for(i = 0; i < n; i++) {
	rho=(x[i]-mu)/sigma0*(x[i]-mu)/sigma0;
	if (rho>c2*c2){
		rho=c2*c2;
	}
	vari=vari+rho;
}
return sigma0*sqrt(vari/n/correc);
}


SEXP scaleTau2(SEXP x, SEXP c1, SEXP c2, SEXP correc) {
	int n=length(x);
	double a=scaleTau(n, REAL(x),asReal(c1),asReal(c2),asReal(correc));
	return ScalarReal(a);
}

static void smoothpsi(int n, double*x, double a, double b, double d, double e, double k, double l) {
int i;
double bx;
for(i = 0; i < n; i++) {
	bx =fabs(x[i]);
	if ((bx>=k)&&(bx<l)) {
		if (x[i]>0)
			{x[i]=a+b*bx+d*x[i]*x[i]+e*bx*bx*bx;}
		if (x[i]<0)
			{x[i]=-(a+b*bx+d*x[i]*x[i]+e*bx*bx*bx);}
		}
	if (bx>=l){
		x[i]=0;
		}
	}
}

double smoothpsi3(double x, double a, double b, double d, double e, double k, double l) {
double bx;
	bx =fabs(x);
	if ((bx>=k)&&(bx<l)) {
		if (x>0)
			{x=a+b*bx+d*x*x+e*bx*bx*bx;}
		if (x<0)
			{x=-(a+b*bx+d*x*x+e*bx*bx*bx);}
		}
	if (bx>=l){
		x=0;
		}
return x;
}

SEXP smoothpsi2(SEXP x, SEXP a, SEXP b, SEXP d, SEXP e, SEXP k, SEXP l) {
	int n=length(x);
	SEXP ans = duplicate(x);
	smoothpsi(n,REAL(ans),asReal(a),asReal(b),asReal(d),asReal(e),asReal(k),asReal(l));
	return ans;
}

static void filterinit(int n,double*timeseries,double autocor,double a, double b, double d, double e, double k, double l, double c1, double c2, double correc, double varpred) {

int i;
double xhat[n+1];
double ug[n];
double sigma,sigmau,P,shat,ugt,varest;
double m;
xhat[0]=0;
sigma=varpred;
sigmau=sigma*(1-autocor*autocor);
P=sigma;
for(i = 1; i < n+1; i++) {
	xhat[i]=autocor*xhat[i-1];
	ug[i-1]=timeseries[i-1]-xhat[i];
	m=autocor*autocor*P+sigmau;
	shat=sqrt(m);
	ugt=smoothpsi3(ug[i-1]/shat,a,b,d,e,k,l);
	P=m-ugt/(ug[i-1]/shat)*m;
	xhat[i]=xhat[i]+shat*ugt;
	}
varest=scaleTau(n,ug,c1,c2,correc);
for(i = 0; i < n; i++) {
	timeseries[i]=xhat[i+1];
	}
timeseries[n]=varest;
for(i = n+1; i < 2*n+1; i++) {
	timeseries[i]=ug[i-(n+1)];
	}
}

SEXP filterinit2(SEXP timeseries, SEXP autocor1, SEXP a, SEXP b, SEXP d, SEXP e, SEXP k, SEXP l,SEXP c1, SEXP c2, SEXP correc, SEXP varpred) {
	int n=(length(timeseries)-1)/2;
	SEXP ans = duplicate(timeseries);
	filterinit(n,REAL(ans),asReal(autocor1),asReal(a),asReal(b),asReal(d),asReal(e),asReal(k),asReal(l),asReal(c1),asReal(c2),asReal(correc),asReal(varpred));
	return ans;
}



static void filter(int n, int p, double*timeseries, double autocor1, double oldvar, double*oldmemory, double a, double b, double d, double e, double k, double l, double c1, double c2, double correc, double*acf) {

int i,j,k1,i2;
double xhat[n+p];
double Phi[p][p];
double P[p][p];
double M[p][p];
double Mz[p][p];
double ug[n];
double shat;
double mhatma[p][p];
double ugt;
double fak;

memset(Phi, 0, sizeof(double) * p * p);
memset(P, 0, sizeof(double) * p * p);
memset(M, 0, sizeof(double) * p * p);
memset(Mz, 0, sizeof(double) * p * p);
Phi[0][p-1]=autocor1;
for(i = 0; i < p-1; i++) {
	Phi[0][i]=oldmemory[i]-autocor1*oldmemory[p-i-2];
	}
for(i = 0; i < p-1; i++) {
	Phi[i+1][i]=1;
	}

for(i = 0; i < p; i++) {
	for (j = 0; j < p; j++) {
		P[i][j]=acf[abs(i-j)];
		}
	}

double sigmau;
sigmau=oldvar*oldvar*(1-autocor1*autocor1);
for(i = 0; i < n+p; i++) {
	xhat[i]=0;
	}

for (i = p; i < n+p; i++) {
    ug[i-p]=timeseries[i-p];
    for (j = 0; j < p; j++) {
        ug[i-p]=ug[i-p]-Phi[0][j]*xhat[i-j-1];
        }
	for (j = 0; j < p; j++) {
		xhat[i]=xhat[i]+xhat[i-j-1]*Phi[0][j];
		}

	
	memset(Mz, 0, sizeof(double) * p * p);
	for(i2 = 0; i2 < p; i2++) {
		for (j = 0; j < p; j++) {
			for (k1 = 0; k1 < p; k1++) {
				Mz[i2][j]=Mz[i2][j]+Phi[i2][k1]*P[k1][j];
				}
			}
		}

	memset(M, 0, sizeof(double) * p * p);
	for(i2 = 0; i2 < p; i2++) {
		for (j = 0; j < p; j++) {
			for (k1 = 0; k1 < p; k1++) {
				M[i2][j]=M[i2][j]+Mz[i2][k1]*Phi[j][k1];
				}
			}
		}

	M[0][0]=M[0][0]+sigmau;	
	for(i2 = 0; i2 < p; i2++) {
		for (j = 0; j < p; j++) {
			mhatma[i2][j]=M[j][0]*M[i2][0];
			}
		}
	shat=sqrt(M[0][0]);
	ugt=smoothpsi3(ug[i-p]/shat,a,b,d,e,k,l);
	fak=1/shat/shat*ugt/ug[i-p]*shat;

	for(i2 = 0; i2 < p; i2++) {
		for (j = 0; j < p; j++) {
			P[i2][j]=M[i2][j]-fak*mhatma[i2][j];
			}
		}

	for (j = 0; j < p; j++) {
		xhat[i-j]=xhat[i-j]+M[j][0]/shat*ugt;
		}
    timeseries[i-p]=xhat[i];
	}
timeseries[n]=scaleTau(n,ug,c1,c2,correc);
for (i = n+1; i < 2*n+1; i++) {
	timeseries[i]=ug[i-(n+1)];
	}
}


SEXP filter2(SEXP timeseries, SEXP autocor, SEXP oldvar, SEXP oldmemory, SEXP a, SEXP b, SEXP d, SEXP e, SEXP k, SEXP l,SEXP c1, SEXP c2, SEXP correc, SEXP acf) {
	int n=(length(timeseries)-1)/2;
	int p=length(oldmemory)+1;
	SEXP ans = duplicate(timeseries);
	filter(n,p,REAL(ans),asReal(autocor),asReal(oldvar),REAL(oldmemory),asReal(a),asReal(b),asReal(d),asReal(e),asReal(k),asReal(l),asReal(c1),asReal(c2),asReal(correc),REAL(acf));
	return ans;
}
 
static void filter3(int n, int p, int q, double*timeseries, double*ar, double*ma, double varpred, double a, double b, double dd, double e, double kk, double l, double c1, double c2, double correc, double*acf) {

int k=0;
if (p>=q+1) {k=p;}
if (q>=p-1) {k=q+1;}
int i,j,k1,i2;
double alpha[n+k];
double Phi[k][k];
double P[k][k];
double M[k][k];
double Mz[k][k];
double ug[n];
double shat;
double mhatma[k][k];
double ugt;
double fak;
double d[k];

memset(d, 0, sizeof(double) * k);
d[0]=1;
for (i = 1; i < q+1; i++) {
    d[i]= -ma[i-1];
    }

memset(Phi, 0, sizeof(double) * k * k);
memset(P, 0, sizeof(double) * k * k);
memset(M, 0, sizeof(double) * k * k);
memset(Mz, 0, sizeof(double) * k * k);

for(i = 0; i < p; i++) {
	Phi[0][i]=ar[i];
	}
for(i = 0; i < k-1; i++) {
	Phi[i+1][i]=1;
	}

for(i = 0; i < k; i++) {
	for (j = 0; j < k; j++) {
		P[i][j]=acf[abs(i-j)];
		}
	}


double sigmau=varpred;
for(i = 0; i < n+k; i++) {
	alpha[i]=0;
	}

for (i = k; i < n+k; i++) {

	for (j = 0; j < p; j++) {
		alpha[i]=alpha[i]+alpha[i-j-1]*Phi[0][j];
		}

    ug[i-k]=timeseries[i-k]-alpha[i];

	memset(Mz, 0, sizeof(double) * k * k);
	for(i2 = 0; i2 < k; i2++) {
		for (j = 0; j < k; j++) {
			for (k1 = 0; k1 < k; k1++) {
				Mz[i2][j]=Mz[i2][j]+Phi[i2][k1]*P[k1][j];
				}
			}
		}

	memset(M, 0, sizeof(double) * k * k);
	for(i2 = 0; i2 < k; i2++) {
		for (j = 0; j < k; j++) {
			for (k1 = 0; k1 < k; k1++) {
				M[i2][j]=M[i2][j]+Mz[i2][k1]*Phi[j][k1];
				}
			}
		}

    for(i2 = 0; i2 < k; i2++) {
        for (j = 0; j < k; j++) {
            M[i2][j]=M[i2][j]+sigmau*d[i2]*d[j];
            }        
        }

	for(i2 = 0; i2 < k; i2++) {
		for (j = 0; j < k; j++) {
			mhatma[i2][j]=M[j][0]*M[i2][0];
			}
		}
	shat=sqrt(M[0][0]);	

	ugt=smoothpsi3(ug[i-k]/shat,a,b,dd,e,kk,l);


	fak=1/shat/shat*ugt/ug[i-k]*shat;

	for(i2 = 0; i2 < k; i2++) {
		for (j = 0; j < k; j++) {
			P[i2][j]=M[i2][j]-fak*mhatma[i2][j];
			}
		}

	for (j = 0; j < k; j++) {
		alpha[i-j]=alpha[i-j]+M[j][0]/shat*ugt;
		}
    timeseries[i-k]=alpha[i];
	}
timeseries[n]=scaleTau(n,ug,c1,c2,correc);
for (i = n+1; i < 2*n+1; i++) {
	timeseries[i]=ug[i-(n+1)];
	}
}


SEXP filter4(SEXP timeseries, SEXP ar, SEXP ma, SEXP varpred, SEXP a, SEXP b, SEXP d, SEXP e, SEXP k, SEXP l,SEXP c1, SEXP c2, SEXP correc, SEXP acf) {
	int n=(length(timeseries)-1)/2;
	int p=length(ar);
    int q=length(ma);
	SEXP ans = duplicate(timeseries);
	filter3(n,p,q,REAL(ans),REAL(ar),REAL(ma),asReal(varpred),asReal(a),asReal(b),asReal(d),asReal(e),asReal(k),asReal(l),asReal(c1),asReal(c2),asReal(correc),REAL(acf));
	return ans;
}
