/* numint.h
   -- numerical integrations via Romberg method --
 */

/* 
Note: when integration from a to infinite use 
routine qromo with 
- function midpnt
- midinf: if function decreases faster than 1/x^2
- midexp: if function decreases faster than exp(-x)
*/

#include"interpol.h"
#include<iostream> 
#include<cmath> 
#include<stdio.h> 
#include<stddef.h> 
#include<stdlib.h> 

#define EPS 1.0e-6
#define K 5
using namespace std;

float qromb(float(*func)(float),float a,float b) 

#define JMAX 20
#define JMAXP (JMAX+1)

{ 
  void polint(float xa[],float ya[],int n,float x,float *y,float *dy); 
  float trapzd(float(*func)(float),float a,float b,int n); 
  float ss,dss; 
  float s[JMAXP], h[JMAXP+1]; 
  int j;
  h[1]=1.0; 
  for(j=1;j<=JMAX;j++){
    s[j]=trapzd(func,a,b,j); 
    if(j>=K){
      polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss); 
      //      cout << j << "\t" << ss << "\t" << fabs(dss) << "\t" << EPS*fabs(ss) << "\n";
      if(fabs(dss)<=EPS*fabs(ss)) return ss; 
    } 
    h[j+1]=0.25*h[j]; 
  } 
  cout << "Too many steps in routine qromb"; 
  return 0.0;
} 

float qromo(float(*func)(float),float a,float b, char type, float acc) 

//#define EPSo 3.0e-7
#define JMAXo 20
#define JMAXPo (JMAXo+1)

{ 
  float EPSo=acc;
  void polint(float xa[],float ya[],int n,float x,float*y,float*dy); 
  float midinf(float(*func)(float),float aa,float bb,int n);
  float midpnt(float(*func)(float),float aa,float bb,int n);
  float midexp(float(*func)(float),float aa,float bb,int n);
  int j; 
  float ss,dss,h[JMAXPo+1], s[JMAXPo]; 
  h[1]=1.0; 
  for(j=1;j<=JMAXo;j++){ 
    if (type == 'p') {s[j]=midpnt(func,a,b,j);}
    if (type == 'i') {s[j]=midinf(func,a,b,j);}
    if (type == 'e') {s[j]=midexp(func,a,b,j);}
    if(j>=K){ 
      polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss); 
      //      cout << j << "\t" << ss << "\t" << fabs(dss) << "\t" << EPSo*fabs(ss) << "\n";
      if(fabs(dss)<=EPSo*fabs(ss)) {
	return ss;
      }
    } 
    h[j+1]=h[j]/9.0; 
  } 
  cout << "Too many steps in routing qromo"; 
  return 0.0; 
} 



float trapzd(float(*func)(float),float a,float b,int n) 
  
#define FUNC1(x)((*func)(x)) 

{ float x,tnm,sum,del; 
  static float s; 
  int it,j; 
  if(n==1){ 
    return(s=0.5*(b-a)*(FUNC1(a)+FUNC1(b))); 
  }else
    { for(it=1,j=1;j<n-1;j++)it<<=1; tnm=it; del=(b-a)/tnm;
      x=a+0.5*del; 
      for(sum=0.0,j=1;j<=it;j++,x+=del) {
	sum+=FUNC1(x);
      }
      s=0.5*(s+(b-a)*sum/tnm); 
      return s; 
    } 
} 


float midinf(float(*funk)(float),float aa,float bb,int n) 

#define FUNC2(x)((*funk)(1.0/(x))/((x)*(x))) 

{ 
  float x,tnm,sum,del,ddel,b,a; 
  static float s; 
  int it,j; 
  b=1.0/aa; 
  a=1.0/bb; 
  if(n==1){ 
    return (s=(b-a)*FUNC2(0.5*(a+b))); 
  }
  else { 
    for(it=1,j=1;j<n-1;j++) it*=3; 
    tnm=it; 
    del=(b-a)/(3.0*tnm); 
    ddel=del+del; 
    x=a+0.5*del; 
    sum=0.0; 
    for(j=1;j<=it;j++) {
      sum+=FUNC2(x); 
      x+=ddel; 
      sum+=FUNC2(x); 
      x+=del;
    } 
    return (s=(s+(b-a)*sum/tnm)/3.0); 
  } 
}



float midpnt(float(*func)(float),float a,float b,int n) 

#define FUNC3(x) ((*func)(x)) 

{ 
  float x,tnm,sum,del,ddel; 
  static float s; 
  int it,j; 

  if (n == 1){ 
    return (s=(b-a)*FUNC3(0.5*(a+b))); 
  }
  else { 
    for(it=1,j=1;j<n-1;j++) it*=3; 
    tnm=it; 
    del=(b-a)/(3.0*tnm); 
    ddel=del+del; 
    x=a+0.5*del; 
    sum=0.0; 
    for(j=1;j<=it;j++) {
      sum+=FUNC3(x); 
      x+=ddel; 
      sum+=FUNC3(x); 
      x+=del;
    } 
    return (s=(s+(b-a)*sum/tnm)/3.0); 
  } 
}



float midexp(float(*funk)(float),float aa,float bb,int n) 

#define FUNC4(x)((*funk)(-log(x))/(x)) 

{ float x,tnm,sum,del,ddel,b,a; 
  static float s; 
  int it,j; 
  b=exp(-aa); 
  a=0.; 
  if(n==1){ 
    return (s=(b-a)*FUNC4(0.5*(a+b))); 
  }
  else { 
    for(it=1,j=1;j<n-1;j++) it*=3; 
    tnm=it; 
    del=(b-a)/(3.0*tnm); 
    ddel=del+del; 
    x=a+0.5*del; 
    sum=0.0; 
    for(j=1;j<=it;j++) {
      sum+=FUNC4(x); 
      x+=ddel; 
      sum+=FUNC4(x); 
      x+=del;
    } 
    return (s=(s+(b-a)*sum/tnm)/3.0); 
  } 
}

