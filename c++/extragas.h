/*
  EXGAS.H

  header file of all the routines (excluded classes) used by exgas
*/

#ifndef EXTRAGAS_H
#define EXTRAGAS_H

#include "array.h"

ffarray::Array3D<float> extragas_model(double v_ki,double ion_i,double accre_i, double &,double &);
//using namespace Magick; 

#endif
