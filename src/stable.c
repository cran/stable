#include <math.h>
#include <stddef.h>

extern double	R_NaN;			/* IEEE NaN or -DBL_MAX */

static int nn;
static double alphai, yi, setai, cetai;

static double fcn1(double s){
  double sa;
  sa=pow(s,alphai);
  return(cos(-yi*s+sa*setai)*exp(-sa*cetai));}

static double fcn2(double s){
  double sa;
  sa=pow(s,-alphai);
  return(cos(-yi/s+sa*setai)*exp(-sa*cetai))/(s*s);}

static void interp(double a[], double fa[], int n, double *f, double *df)
{
  int i, j, ni=1;
  double diff1, diff2, tmp1, tmp2, lim1, lim2, *tab1, *tab2;
 
  tmp1=fabs(a[1]);
  tab1=(double*)malloc((size_t)((n+1)*sizeof(double)));
  tab2=(double*)malloc((size_t)((n+1)*sizeof(double)));
  if(!tab1||!tab2)return;
  for(i=1;i<=n;i++){
    tmp2=fabs(a[i]);
    if(tmp2<tmp1){
      ni=i;
      tmp1=tmp2;}
    tab1[i]=fa[i];
    tab2[i]=fa[i];}
  *f=fa[ni--];
  for(j=1;j<n;j++){
    for(i=1;i<=n-j;i++){
      lim1=a[i];
      lim2=a[i+j];
      diff1=tab1[i+1]-tab2[i];
      diff2=lim1-lim2;
      if(diff2==0.0)goto end;
      diff2=diff1/diff2;
      tab2[i]=lim2*diff2;
      tab1[i]=lim1*diff2;}
    *df=2*ni<(n-j)?tab1[ni+1]:tab2[ni--];
    *f+=*df;}
 end: free((char *)tab2);
  free((char *)tab1);}

#define FUNC(x) ((*fcn)(x))

static double mpr(double (*fcn)(double), int n)
{
  double x, nn, tmpsum, pnt1, pnt2;
  static double sum;
  int i,j;

  if (n==1){
    sum=FUNC(0.5);
    return(sum);}
  else {
    for(i=1,j=1;j<n-1;j++) i*=3;
    nn=i;
    pnt1=1/(3.0*nn);
    pnt2=2.0*pnt1;
    x=0.5*pnt1;
    tmpsum=0.0;
    for(j=1;j<=i;j++){
      tmpsum+=FUNC(x);
      x+=pnt2;
      tmpsum+=FUNC(x);
      x+=pnt1;}
    sum=(sum+tmpsum/nn)/3.0;
    return(sum);}}

double romberg(double (*fcn)(double), double eps)
{
  int j;
  double sum,errsum,x[17],fx[17];

  x[1]=1.0;
  for(j=1;j<=16;j++){
    fx[j]=mpr(fcn,j);
    if(j>=5){
      interp(&x[j-5],&fx[j-5],5,&sum,&errsum);
      if(fabs(errsum)<eps*fabs(sum))return(sum);}
    x[j+1]=x[j]/9.0;
    fx[j+1]=fx[j];}
  return(R_NaN);}

void stable(int *n, double *y, double *beta, double *alpha, int *npt, double *up, double *eps, int *type, int *err, double *ffy)
{
  int i, j;
  double h, s, *eta, *seta, *ceta, *sa;
  *err=0;
  eta=(double*)malloc((size_t)((*n+1)*sizeof(double)));
  seta=(double*)malloc((size_t)((*n+1)*sizeof(double)));
  ceta=(double*)malloc((size_t)((*n+1)*sizeof(double)));
  sa=(double*)malloc((size_t)((*n+1)*sizeof(double)));
  nn=*n;
  if(!eta||!seta||!ceta||!sa){
    *err=1;
    return;}
  for(i=0;i<*n;i++){
    ffy[i]=0.0;
    eta[i]=beta[i]*(1.0-fabs(1.0-alpha[i]))*M_PI/2.0;
    seta[i]=sin(eta[i]);
    ceta[i]=cos(eta[i]);}
  if(*type==1){
    *npt=*npt-*npt%2;
    h=*up/ *npt;
	   for(j=0;j<=*npt;j++){
	     s=(*npt-j)*h;
	     for(i=0;i<*n;i++){
	       sa[i]=pow(s,alpha[i]);
	       ffy[i]=ffy[i]+(4-2*(j%2==0)-(j==1||j==*npt))*cos(-y[i]*s+sa[i]*seta[i])*exp(-sa[i]*ceta[i]);}}
	   for(i=0;i<*n;i++)ffy[i]=ffy[i]*h/3.0/M_PI;}
  else {
    for(i=0;i<*n;i++){
      alphai=alpha[i];
      yi=y[i];
      setai=seta[i];
      cetai=ceta[i];
      ffy[i]=romberg(fcn1, *eps)+romberg(fcn2, *eps);}
    for(i=0;i<*n;i++)ffy[i]=ffy[i]/M_PI;}
  free((char *)sa);
  free((char *)ceta);
  free((char *)seta);
  free((char *)eta);}
