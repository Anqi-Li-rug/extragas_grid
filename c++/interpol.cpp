/* interpol.h

   Interpolation and extrapolation functions
*/

#ifndef INTERPOL_H
#define INTERPOL_H

#include<iostream>
#include<cstdlib>
#include<cmath> 
#include<cstdio> 
#include<stddef.h> 
#include "interpol.h"

#define NR_END 1
#define FREE_ARG char*
using std::cout;

float *vector(long nl,long nh) 
{ float *v;
  v=(float*)malloc((size_t)((nh-nl+1+NR_END)*sizeof(float))); 
  if(!v) cout << "allocation failure in vector ()"; 
  return v-nl+NR_END; 
} 

void free_vector(float *v,long nl,long nh) 
{ 
  free((FREE_ARG) (v+nl-NR_END)); 
} 

void polint(float xa[],float ya[],int n,float x,float *y,float *dy)
{
  int i,m,ns=1;
  float den,dif,dift,ho,hp,w;

  float *c,*d; 
  dif=fabs(x-xa[1]); 
  c=vector(1,n); 
  d=vector(1,n); 
  for(i=1;i<=n;i++){ 
    if((dift=fabs(x-xa[i]))<dif){ ns=i; dif=dift; } 
    c[i]=ya[i]; 
    d[i]=ya[i]; 
  }
  *y=ya[ns--]; 
  for(m=1;m<n;m++){ 
    for(i=1;i<=n-m;i++) { 
      ho=xa[i]-x; 
      hp=xa[i+m]-x; 
      w=c[i+1]-d[i]; 
      if((den=ho-hp)==0.0) cout << "Error in routine polint"; 
      den=w/den; 
      d[i]=hp*den; 
      c[i]=ho*den; 
    } 
    *y+=(*dy=(2*ns<(n-m)?c[ns+1]:d[ns--])); 
  } 
  free_vector(d,1,n); 
  free_vector(c,1,n);
}

void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]) 
{ 
  int i,k; 
  float p,qn,sig,un,*u; 
  u=vector(1,n-1); 
  if (yp1 > 0.99e30) 
    y2[1]=u[1]=0.0; 
  else { 
    y2[1] = -0.5; 
    u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1); 
  }
  for (i=2;i<=n-1;i++) { 
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]); 
    p=sig*y2[i-1]+2.0; 
    y2[i]=(sig-1.0)/p; 
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]); 
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p; 
  }
  if (ypn > 0.99e30) 
    qn=un=0.0; 
  else { 
    qn=0.5; 
    un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1])); 
  } 
  y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0); 
  for (k=n-1;k>=1;k--) 
    y2[k]=y2[k]*y2[k+1]+u[k]; 
  free_vector(u,1,n-1); 
}

void splint(float xa[], float ya[], float y2a[], int n, float x, float *y) 
{ 
  //  void nrerror(char error_text[]); 
  int klo,khi,k; 
  float h,b,a; 
  klo=1; 
  khi=n; 
  while (khi-klo > 1) { 
    k=(khi+klo) >> 1; 
    if (xa[k] >x) khi=k; 
    else klo=k; 
  } 
  h=xa[khi]-xa[klo]; 
  if (h ==0.0) cout << "Bad xa input to routine splint";
  a=(xa[khi]-x)/h; 
  b=(x-xa[klo])/h; 
  *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0; 
}

/*
void polin2(float x1a[],float x2a[],float **ya,int m,int n,float x1, float x2,float *y, float *dy) 
{
  void polint(float xa[],float ya[],int n,float x,float *y,float *dy);
  int j;
  float *ymtmp; 
  ymtmp=vector(1,m); 
  for (j=1;j<=m;j++){
  polint(x2a,ya[j],n,x2,&ymtmp[j],dy); 
  }
  polint(x1a,ymtmp,m,x1,y,dy); 
  free_vector(ymtmp,1,m); 
} 
*/


void polin2(float x1a[],float x2a[],float *psubarray, int m, int n, float x1, float x2, float *y, float *dy) 
{
  void polint(float xa[],float ya[],int n,float x,float *y,float *dy);
  int j,i,k;
  float *ymtmp; 
  ymtmp=vector(1,m);
  /*  float subarray[m][n];
      for (i=0;i<m;i++){
      for (j=0;j<n;j++){
      subarray[i][j]=*(psubarray+i*m+j);
      cout << i << "\t" << j << "\t" << subarray[i][j] << "\n";
      }
      }
  */
  float *subarray = new float[n];
  for (int i=0;i<m;i++){
    for (int j=0;j<n;j++){
      subarray[j]=*(psubarray+i*n+j);
    }
    polint(x2a,subarray,n-1,x2,&ymtmp[i],dy);
  }
  polint(x1a,ymtmp,m,x1,y,dy);
  free_vector(ymtmp,1,m);
  delete [] subarray;
}


float bilin(float subx1[], float subx2[], float *psubarray, float x1, float x2)
{
  // BI-LINEAR INTERPOLATION
  float t,u,y,subarray[2][2];

  subarray[0][0]=*(psubarray);
  subarray[1][0]=*(psubarray+2);
  subarray[1][1]=*(psubarray+3);
  subarray[0][1]=*(psubarray+1);
  
  t=(x1-subx1[0])/(subx1[1]-subx1[0]);
  u=(x2-subx2[0])/(subx2[1]-subx2[0]);
  y=(1-t)*(1-u)*subarray[0][0]+t*(1-u)*subarray[1][0]+t*u*subarray[1][1]+(1-t)*u*subarray[0][1];
  return y;
}

void splie2(float x1a[], float x2a[], float *ya, int m, int n, float *y2a)
{ 
  void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]); 
  float *subarray = new float[n];
  for (int i=0;i<m;i++){
    for (int j=0;j<n;j++){
      subarray[j]=*(ya+i*n+j);
    }
    spline(x2a,subarray,n,1.0e30,1.0e30,y2a);    
  }
    
}

void splin2(float x1a[], float x2a[], float *ya, float *y2a, int m, int n, 
float x1, float x2, float *y) 
{
  void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]); 
  void splint(float xa[], float ya[], float y2a[], int n, float x, float *y); 
  float *ytmp,*yytmp; 
  ytmp=vector(1,m); 
  yytmp=vector(1,m); 
  //  for (j=1;j<=m;j++) 
  //    splint(x2a,ya[j],y2a[j],n,x2,&yytmp[j]); 

  float *subarray = new float[n];
  for (int i=0;i<m;i++){
    for (int j=0;j<n;j++){
      subarray[j]=*(ya+i*n+j);
    }
    splint(x2a,subarray,y2a,n,x2,&yytmp[i]); 
  }
  spline(x1a,yytmp,m,1.0e30,1.0e30,ytmp); 
  free_vector(yytmp,1,m); 
  free_vector(ytmp,1,m); 
  delete [] subarray;
}
#endif
