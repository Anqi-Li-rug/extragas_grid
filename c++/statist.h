/* statist.h
   statistics routines, random number generators
 */

#ifndef STATIST_H
#define STATIST_H

float ran0(long *);
float ran_gau(long int &, float );
double ran_gau(long int &, double );
float ran_lin(long int &, float );
float ran_pow(long int &, float , float , float );
float ran_exp(long int &, float , float );
float ran_cos(long int &);

#endif
