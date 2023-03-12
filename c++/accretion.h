/*
  DENSITY and VELOCITY PATTERNS for COLD ACCRETION
*/

#include "array.h"
using namespace ffarray;
#define inline

class Accretion
{
 private: 
  float r0, v_r0, v_theta0, v_phi0, sigma_v0;

 public:
  Accretion(float, float, float, float, float);

  void Patterns(Array2D<double>&,
		  Array2D<double>&,
		  Array2D<double>&,
		  Array2D<double>&,
		  Array1D<double>, 
		  Array1D<double>,
		  Array2D<double>&,
		  Array2D<double>&,
		  Array2D<double>&
		  );
  
  ~Accretion();
};
