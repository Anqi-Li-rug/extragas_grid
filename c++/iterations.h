#include "array.h"
using namespace std;
using namespace ffarray;

class iterations
{
private:
  char par1[40], par2[40], par3[40];
  int n_iter1, n_iter2, n_iter3, RAsize, DECsize, VELsize;
  float x_shift, y_shift, kpctograd, cdelt1, cdelt2;
  // Distance above the plane from which residuals are calculated 
  float res_above;
  // Radial size where residuals are calculated 
  float res_Rsize;
  // Vertical size where residuals are calculated 
  float res_zsize;
  // Range of channel maps where noise is calculated 
  int freechannels;
  /* Range of channel maps where residuals are calculated in target mode*/
  int target_chan1;
  int target_chan2;
  float par1_0, par2_0, par3_0, delta_par1, delta_par2, delta_par3;
  int above_pix, below_pix, Rsize_pix, zsize_pix;
  float absnoise, noise;
  Array2D<float> above_tot, below_tot, above_left_tot, below_left_tot,
    average_tot;
  Array2D<float> above, below, above_app, below_app, average, targeted;
public:
	iterations(char param1[], char param2[], char param3[],
		       int &niter1, int &niter2, int &niter3,
		       float &param1_0, float &delta_param1,
		       float &param2_0, float &delta_param2,
			   float &param3_0, float &delta_param3);
  void CalculateResiduals(int, int, Array3D<float>, Array2D<float>, char);
  void WriteResiduals();
  //~iterations();
};
