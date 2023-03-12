/* statist.h
   statistics routines, random number generators
 */

#include <cstdlib>
#include <cmath>
//#include "statist.h"

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define MASK 123459876

float ran0(long *idum)
{
  long k;
  float ans;

  *idum^=MASK;
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  ans=AM*(*idum);
  *idum^=MASK;
  return ans;
}

float ran_gau(long int &seed, float sigma)
{
  /* 
     positive and negative sides
  */
	
	if(sigma==0) {return 0;}
	
  float ran0(long*);
  float num;
  do {
    num=ran0(&seed);
  } while (num==0.0); 
  
  float t=sqrt(log(1/num/num));
  float c0=2.515517;
  float c1=0.802853;
  float c2=0.010328;
  float d1=1.432788;
  float d2=0.189269;
  float d3=0.001308;
  float ranval= sigma*(t-(c0+c1*t+c2*t*t)/(1+d1*t+d2*t*t+d3*t*t*t));
  
  return ranval;
}

double ran_gau(long int &seed, double sigma)
{
	if(sigma==0) {return 0;}
  /* 
     positive and negative sides
  */
  float ran0(long*);
  float num;
  do {
    num=ran0(&seed);
  } while (num==0.0); 
  
  float t=sqrt(log(1/num/num));
  float c0=2.515517;
  float c1=0.802853;
  float c2=0.010328;
  float d1=1.432788;
  float d2=0.189269;
  float d3=0.001308;
  float ranval= sigma*(t-(c0+c1*t+c2*t*t)/(1+d1*t+d2*t*t+d3*t*t*t));
  
  return ranval;
}

float ran_lin(long int &seed, float ab)
{
  /* 
     normalized
     f=2b/a*(1-b/a*x), ab=a/b
  */
  float ran0(long*);
  float num;
  do {
    num=ran0(&seed);
  } while (num==0.0); 
  float ranval = ab*(1-sqrt(1-num));

  return ranval;
}

float ran_pow(long int &seed, float alpha, float h, float x_cut)
{
  /*
    f=(alpha-1)/h*(x_cut/h+1)**(alpha-1)*(x/h+1)**(-alpha)
  */
  float ran0(long*);
  float num;
  do {
    num=ran0(&seed);
  } while (num==0.0); 
  float ranval=h*((x_cut/h+1)*pow(1-num,1/(1-alpha))-1);

  return ranval;
}

float ran_exp(long int &seed, float h, float x_cut)
{
  /*
    normalized
    f=1/h*exp(-((x-x_cut)/h))
   */
  float ran0(long*);
  float num;
  do {
    num=ran0(&seed);
  } while (num==0.0); 
  float ranval=-h*log(1-num)+x_cut; 

  return ranval;
}

float ran_cos(long int &seed)
{
  /*
   */
  float ran0(long*);
  float num;
  do {
    num=ran0(&seed);
  } while (num==0.0); 
  float ranval=asin(num);

  return ranval;
}
