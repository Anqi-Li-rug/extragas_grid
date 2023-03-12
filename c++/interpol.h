/* interpol.h

   Interpolation and extrapolation functions
*/

#ifndef INTERPOL_H
#define INTERPOL_H

float *vector(long ,long );
void free_vector(float *,long ,long );
void polint(float ,float ,int ,float ,float *,float *);
void spline(float , float , int , float , float , float );
void splint(float , float , float , int , float , float *);
void polin2(float ,float ,float *, int , int , float , float , float *, float *);
float bilin(float , float , float *, float , float );
void splie2(float , float , float *, int , int , float *);
void splin2(float , float , float *, float *, int , int , float , float , float *);
#endif
