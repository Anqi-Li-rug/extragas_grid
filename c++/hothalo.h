/*
  HOT HALO

*/

#include "array.h"
using namespace ffarray;
#define inline

class HotHalo
{
 private:
  float v0, T, Rmax, zmin, zmax, mu;
  double lum;
  //  Array2D<double> rho_ik1(), rho_ik2();
  float norm_lum, mass_halo, mass_halo2;
  double conv_X;

 public:
  HotHalo(float, float, float, float, float, float, double);

  void SetDensity(Array2D<double>&,
		  const Array1D<double>, 
		  const Array1D<double>,
		  const Array2D<double>&,
		  const Array2D<double>&,
		  Array2D<double>&, 
		  const Array1D<double>, 
		  const Array1D<double>,
		  const Array2D<double>&,
		  const Array2D<double>&);
  
  void SetDensity(Array2D<double>&,
		  const Array1D<double>, 
		  const Array1D<double>,
		  const Array2D<double>&,
		  const Array2D<double>&);

  ~HotHalo();
};
