#include <math.h>
#include <stddef.h>

extern double	R_NaN;			/* IEEE NaN or -DBL_MAX */

static int nn;
static double alphai, etai, setai, cetai, yyi;

static double func1(double s){
  double sa;
  sa=pow(s,alphai);
  return((sin(yyi*s-sa*setai)/s)*exp(-sa*cetai));}

static double func2(double s){
  double sa;
  sa=pow(s,-alphai);
  return((sin(yyi/s-sa*setai)*s)*exp(-sa*cetai))/(s*s);}

void pstable(int *n, double *y, double *beta, double *alpha, double *eps, int *err, double *ffy)
{
  int i, j;
  double h, s;
  *err=0;
  nn=*n;
  for(i=0;i<*n;i++){
    ffy[i]=0.0;
    etai=beta[i]*(1.0-fabs(1.0-alpha[i]))*M_PI/2.0;
    setai=sin(etai);
    cetai=cos(etai);
    alphai=alpha[i];
    yyi=y[i];
    if(etai==0.&&yyi==0)
      ffy[i]=0.5;
    else {
      ffy[i]=romberg(func1, *eps)+romberg(func2, *eps);
      ffy[i]=0.5+ffy[i]/M_PI;}}}

