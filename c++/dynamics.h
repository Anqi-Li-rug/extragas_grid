/*
  DYNAMICS
  
  Dynamics is a class that contains all the function to perform dynamical 
  integration (in cartisian coordinates) given:
  1) initial positions and velocities;
  2) 2 matrices for the forces;
  3) parameters for drag forces;
  4) 
*/

#include <iostream>
#include <cmath>
#include <ctime>
#include "array.h"
#include "statist.h"
#include "galFunctions.h"
using namespace ffarray;
#define inline

//template <class T>
class dynamics
{
 protected:
  ffarray::Array2D<double> F_R_ik, F_z_ik, pot_ik, rho_hot_ik;
  ffarray::Array1D<double> R_i, z_k;
  double x, y, z;
  double dt;
  double R, phi;
  double g_x, g_y, g_z;
  double v_accr_x, v_accr_y, v_accr_z, v_accr_R, v_accr_theta, v_accr_phi;
  double vrel_x, vrel_y, vrel_z, vrel_tot;
  long seed,seed2;
  float K_D, v0_hot;
  double Vrel, sigmaV;
  double M_cloud, R_cloud, rho_corona, drag_coeff;
  double alpha_accr;
  float n_dot, accr_z_in, accr_z_fin, dn_accr, v_esc;
  ffarray::Array2D<double> subF_R_ik, subF_z_ik, sub_rho_ik;
  ffarray::Array1D<double> subR_i, subz_k;
  int i0, j0, subdim1, subdim2;
  bool ArraysCreated;
  char effects[20], accr_type[20], v_accr_type[20];
 public:
  dynamics(Array2D<double>&, 
	   Array2D<double>&, 
	   Array2D<double>&, 
	   Array1D<double>, 
	   Array1D<double>, 
	   int, 
	   int, 
	   char[], 
	   float, 
	   float, 
	   Array2D<double>&,
	   double,  //relative AZIMUTHAL velocity between disk and corona (km/s)
	   double,	//sigmaV (km/s)
	   double,  //cloud mass (Mo)
	   double, 	//cloud radius (pc)
	   double,  //corona density (g/cm3)
	   double  //delay time (Myr)
		   );
  dynamics(Array2D<double>&, Array2D<double>&, Array1D<double>, Array1D<double>);
  //  dynamics(Array2D<float>&, Array2D<float>&, Array1D<float>, Array1D<float>);
  //  float, float, float);
	void set_accr_z_in(float value){accr_z_in = value;}
  void set_accr_z_fin(float value){accr_z_fin = value;}
  void rk4(double,
	   double, 
	   double&, 
	   double&, 
	   double&, 
	   double&, 
	   double&, 
	   double&, 
	   float);
  void getForce(double, double, double, double, double, double, double);
  void getForceMem(double, double, double, double, double, double);
  void DragForce(double, double, double, double, double, double, 
		 double&, double&, double&);
  void DragForceMem(double, double, double, double&, double&, double&);
  void AccrForce(double, double, double, double&, double&, double&);
  void AccrDragForce(double, double, double, double, double, double&, double&, double&);
  void AccrNoDragForce(double, double, double, double, double, double&, double&, double&);
  float accretion(double, double, double, double,
		  char[], char[], float, float);
  float accretion(double, double, double, double,
		  Array2D<double>&, 
		  Array2D<double>&, Array2D<double>&, Array2D<double>&, 
		  Array1D<double>, Array1D<double>,
		  float);
  ~dynamics();
};
