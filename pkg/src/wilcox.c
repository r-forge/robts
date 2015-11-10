#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdlib.h>

static void wilcoxsuk(int n, double *x)
{
int suk[n];
int i,j;
double test[n];


suk[0] = 0;
for(i = 1; i < n; i++) {
	if (x[0]>x[i]) {
		suk[0]++;
		}
	}
test[0] = sqrt(n-1)*(2*suk[0]-(n-1))/sqrt(n)/n/n;
for(i = 1; i < n-1; i++) {
	suk[i]=suk[i-1];
	for(j = i + 1;j < n;j++){
		if (x[i]>x[j]) {
			suk[i]++;
			}
		}
	for(j = 0;j < i; j++){
		if (x[i]<x[j]) {
			suk[i]--;
			}
		}
	test[i]=sqrt((i+1)*(n-(i+1)))*(2*suk[i]-(i+1)*(n-(i+1)))/sqrt(n)/n/n;
	}
for(i = 0; i < n-1; i++) {
	x[i] = test[i];
	}
}

SEXP wilcoxsukz(SEXP data)
{
    int n = LENGTH(data);
    SEXP ans = duplicate(data);
    wilcoxsuk(n,REAL(ans));
    return ans;
}

int compare(const void *x, const void *y)
{
  double xx = *(double*)x, yy = *(double*)y;
  if (xx < yy) return -1;
  if (xx > yy) return  1;
  return 0;
}

double meddiff2(int n1, int n2, double*x, double*y) {
int i,j;
double diff[n1*n2];

    for(i = 1; i < n1; i++) {
	for(j = 1; i < n2; j++) {
		diff[i*n2+j]=x[i]-y[j];
		}
	}
    qsort(diff,n1*n2,sizeof(double),compare);
    if (n1*n2 % 2==0) { 
	return diff[n1*n2/2];
	}
    if (n1*n2 % 2==0) { 
	return (diff[(n1*n2-1)/2]+diff[(n1*n2+1)/2])/2;
	}
}
double meddiff3(int n1, int n2, double*x, double*y) {
int i,j;
double diff[n1*n2];

    for(i = 0; i < n1; i++) {
	for(j = 0; j < n2; j++) {
		diff[i*n2+j]=x[i]-y[j];
		}
	}
	
qsort(diff,n1*n2,sizeof(double),compare);
return diff[n1*n2/2];
    }



SEXP meddiff(SEXP x, SEXP y)
{
double a;
    int n1 = LENGTH(x);
    int n2 = LENGTH(y);
	a=meddiff3(n1,n2,REAL(x),REAL(y));
	return ScalarReal(a);
}


static void meddiffneu2(int n, double*x, double*i1, double*i2,double*erg) {
int i;
int j;
int nn=n*(n+1)/2;
double x1;
int grenze;
for(j = 0; j < n; j++) {
	int k = 0;
	grenze = (j+1)*(n-j);
	for (i = 0; i < nn; i++) {
		if ((i1[i]<=(j+1))&&(i2[i]>j+1)) {k++;
			if (k==(grenze/2)) {
				x1=x[i];
				}
			if (k==((grenze+1)/2)){
				erg[j]=x[i];
				if ((grenze+1)%2==0) {break;}
				}
			if (k==(grenze/2+1)) {
				erg[j]=(x1+x[i])/2;
				break;
				}
			}
		}
	}
}

SEXP meddiffneu(SEXP x,SEXP i1, SEXP i2, SEXP erg) {
	int n=LENGTH(erg);
	SEXP ans = duplicate(erg);
	meddiffneu2(n,REAL(x),REAL(i1),REAL(i2),REAL(ans));
	return ans;
}


static void runmean2(int n, double *x, int l) {
int i;
double x2[n-l]; 
double start=0;
for (i = 0; i < l; i++) {
	start=start+x[i];
	}
x2[0]=start;
for (i = 0; i < n-l; i++) {
	x2[i+1]=x2[i]-x[i]+x[i+l];
	}
for (i = 0; i < n-l; i++) {
	x[i]=x2[i];
	}
}

SEXP runmean(SEXP x, SEXP l) {
	int n=LENGTH(x);
	SEXP ans = duplicate(x);
	runmean2(n,REAL(ans),asInteger(l));
	return ans;
}




