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

#include "interpol.h"

float qromb(float(*)(float),float ,float );
float qromo(float(*)(float),float ,float , char , float );
float trapzd(float(*)(float),float ,float ,int );
float midinf(float(*)(float),float ,float ,int );
float midpnt(float(*)(float),float ,float ,int );
float midexp(float(*)(float),float ,float ,int );
