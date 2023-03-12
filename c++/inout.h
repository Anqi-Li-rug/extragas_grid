/*
  INOUT

  inout is a class that deals with input and output rates in exgas.
  
*/

#include "array.h"
#include <vector>
using namespace ffarray;
#define inline

class inout
{
 protected:
    ffarray::Array1D<float> dens_out, dens_in,angular_in, dens_accr, 
    dens_out_above, dens_in_above, 
    tmass_out, tmass_in, tmass_accr,
    tmass_out_above, tmass_in_above;
  float tau, n_tau, deltaR;
  double Rmax_sf, Rmin_in, Rmin_sf, R1_in, R2_in,angular_fin,angular_fout,m_fin,m_fout,v_fin,v_fout;
  int i_in, i_out, above, nR, ii;
  char single_event; 
  double total_energy_above, ctotal_energy_above;
  float Global_acc_rate;	
  bool found, found_above;
 public:
  inout(int, float );
  //void CalculateRates(double,double, double,double,float, float, float, float, int, float, float, float, float, float, float, int, double&);
void CalculateRates(double,float, float, float, float, int, float, float, float, float, float, float, int, double&);  
void NormAndPlot(float, float, float);
  float get_GlobalAR () {return Global_acc_rate;}; 
	void back_to_zero (); //all quantities go to zero.
  ~inout();
};

