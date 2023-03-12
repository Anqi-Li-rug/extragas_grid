extern const double CONV1, PI, RAD;
extern char inputfile[];
const double CONV3=4.04742e-8; // Mo/kpc^3 -> H_atoms/cm^-3
const double K_B=1.38066e-16; // Boltzmann constant
const double M_P=1.672649e-24; // proton mass (g)
const double KEV=1.1604e7; // kev -> K
const double KPC=3.08567802e21; // kpc (cm)

#include <cmath>
#include "accretion.h"
#include "dynamics.h"
#include "gnuplot.h"
#include "findValue.h"
#include <string.h>
extern long int seed;

using namespace ffarray;
void plot();

Accretion::Accretion(float r,
		     float vr, 
		     float vphi,
		     float vtheta,
		     float sigmav0)
{
  r0=r;
  v_r0=vr;
  v_theta0=vtheta;
  v_phi0=vphi;
  sigma_v0=sigmav0;
}

Accretion::~Accretion()
{}

void Accretion::Patterns(Array2D<double>& rho_cold_ik,
			 Array2D<double>& cold_v_R_ik,
			 Array2D<double>& cold_v_phi_ik,
			 Array2D<double>& cold_v_z_ik,
			 Array1D<double> R_i1, 
			 Array1D<double> z_k1,
			 Array2D<double>& pot_ik1,
			 Array2D<double>& F_R1,
			 Array2D<double>& F_z1)
{
  /* 
     DENSITY AND VELOCITY PATTERN OF THE ACCRETING MATERIAL
  */
  FILE* f_accr_orbit=NULL;
  FILE* f_accretion_R=NULL;
  FILE* f_accretion_z=NULL;
  f_accretion_R=fopen("accretion_R.dat", "w");
  f_accretion_z=fopen("accretion_z.dat", "w");
  f_accr_orbit=fopen("accr_orbit.dat", "w");

  float Etot, Ltot, L_z;
  int nR1=R_i1.dim();
  int nz1=z_k1.dim();
  float stepR1=R_i1[2]-R_i1[1];
  float stepz1=z_k1[2]-z_k1[1];
  bool keep_integr, keep_going=false;

  int n_accr=1000;
  float t_max;
  float dt=0.3;
  int i=0;

  dynamics accretion(F_R1, F_z1, R_i1, z_k1);  
  Array2D<double> trho_cold_ik(nR1,nz1,0.f),
    tcold_v_R_ik(nR1,nz1,0.f), 
    tcold_v_phi_ik(nR1,nz1,0.f), 
    tcold_v_z_ik(nR1,nz1,0.f);
  char accr_int_end[10]; findValue(inputfile, accr_int_end, "accr_int_end");
  float accr_R_max; findValue(inputfile, accr_R_max, "accr_R_max");

  while (i < n_accr){

    /* settings */
    //    out_boundaries=false;
    long seed=time(NULL)*((i+1)*15+3)*(int)sqrt((i+1)*(i+1)*7+(i+1)^i);

    //    float theta=ran0(&seed)*PI/2.;
    float theta=ran_cos(seed);
    float ran_r=ran_gau(seed,r0*0.03);

    // Initial positions
    float R=(r0+ran_r)*cos(theta);
    float phi=0;
    float z=(r0+ran_r)*sin(theta);
    
    // Initial velocities
    float v_R=ran_gau(seed,sigma_v0);
    float v_phi=ran_gau(seed,sigma_v0)+v_phi0;
    float v_z=ran_gau(seed,sigma_v0);

    //    cout << ran_r << " " << R << " " << z << " " << v_R << " " << v_phi << " " << v_z << endl;

    double t=0., x1=R, x2=phi, x3=z, x4=v_R*CONV1, x5=v_phi*CONV1, x6=v_z*CONV1;

    if (strcmp(accr_int_end, "first0") == 0){
      t_max=1e3;
    } else {
      keep_going=true; // for integration up to t=t_max
      findValue(inputfile, t_max, "accr_int_end");
    }

    cyltocar(x1,x2,x3,x4,x5,x6);
    keep_integr=true;

    while (t <= t_max and keep_integr == true) {
      
      /* Check if the particle is inside the potential grid */
      if ((x1 < R_i1[nR1-3]) && (fabs(x3) <z_k1[nz1-3])){

	cyltocar(x1,x2,x3,x4,x5,x6);
	accretion.rk4(t,dt,x1, x2, x3, x4, x5, x6, 0);
	cartocyl(x1,x2,x3,x4,x5,x6);

	Etot=(.5*(x4*x4+x5*x5+x6*x6)
	      +getValue(pot_ik1,R_i1,z_k1,x1,fabs(x3),1));
	Ltot=sqrt(x3*x3*x5*x5+(x3*x4-x1*x6)
		  *(x3*x4-x1*x6)+x1*x1*x5*x5);
	L_z=(x1*x5);

	//	fprintf(f_accr_orbit,
	//		"%.4f %.4f %.4f %.4f %.3f %.3f %.3f %.9f %.5f %.9f \n",
	//		t/1000.,x1,x2,x3,x4/CONV1,x5/CONV1,x6/CONV1,
	//		Etot,Ltot,L_z);
	
	for (int ai=0; ai<nR1-1; ai++) 
	  if (R_i1[ai] <= x1 && x1 < R_i1[ai+1]) 
	    for (int ak=0; ak<nz1-1; ak++) 
	      if (z_k1[ak] <= fabs(x3) && fabs(x3) < z_k1[ak+1]) {
		trho_cold_ik[ai][ak]++;
		tcold_v_R_ik[ai][ak]+=x4/CONV1;
		tcold_v_phi_ik[ai][ak]+=x5/CONV1;
		tcold_v_z_ik[ai][ak]+=x6/CONV1;
	      }
      }
      else {
	keep_integr=false;
      }
      /* Conditions to keep on integrating */
      if ((x3 <= 0) && (x1 < accr_R_max) && keep_going == false)
	keep_integr=false;
      t+=dt;
    }
    i++;
  }
  fclose(f_accr_orbit);

  int smoo=3;

  for (int ai=0; ai<nR1-1; ai++) 
    for (int ak=0; ak<nz1-1; ak++) {
      if (trho_cold_ik[ai][ak] != 0) {
	tcold_v_R_ik[ai][ak]/=trho_cold_ik[ai][ak];
	tcold_v_phi_ik[ai][ak]/=trho_cold_ik[ai][ak];
	tcold_v_z_ik[ai][ak]/=trho_cold_ik[ai][ak];
      }
    }
  rho_cold_ik=trho_cold_ik.smooth2D_keepsize(smoo);
  cold_v_R_ik=tcold_v_R_ik.smooth2D_keepsize(smoo);
  cold_v_phi_ik=tcold_v_phi_ik.smooth2D_keepsize(smoo);
  cold_v_z_ik=tcold_v_z_ik.smooth2D_keepsize(smoo);
  
  rho_cold_ik=rho_cold_ik/total(rho_cold_ik);

  //  plot();
}

void plot()
{
  Gnuplot_Tmpfile gp;
  gp.begin();
  gp.commandln("set terminal postscript portrait enhanced 10");
  gp.commandln("set out 'accretion.ps'");
  gp.commandln("set multiplot");
  gp.commandln("set nokey");
  gp.commandln("set size 0.33,0.25");
  gp.commandln("set origin 0,0.75");
  gp.commandln("set xlabel 't (Gyr)'");
  gp.commandln("set ylabel 'R (kpc)'");
  gp.commandln("plot 'accr_orbit.dat' using 1:2 with poi 0");
  gp.commandln("set ylabel 'phi'");
  gp.commandln("set origin 0,0.5");
  gp.commandln("plot 'accr_orbit.dat' using 1:3 with poi 0");
  gp.commandln("set origin 0,0.25");
  gp.commandln("set ylabel 'z (kpc)'");
  gp.commandln("plot 'accr_orbit.dat' using 1:4 with poi 0");
  gp.commandln("set origin 0.33,0.75");
  gp.commandln("set ylabel 'v_R (km/s)'");
  gp.commandln("plot 'accr_orbit.dat' using 1:5 with poi 0");
  gp.commandln("set origin 0.33,0.5");
  gp.commandln("set ylabel 'v_{phi} (km/s)'");
  gp.commandln("plot 'accr_orbit.dat' using 1:6 with poi 0");
  gp.commandln("set origin 0.33,0.25");
  gp.commandln("set ylabel 'v_z (km/s)'");
  gp.commandln("plot 'accr_orbit.dat' using 1:7 with poi 0");
  gp.commandln("set origin 0.66,0.75");
  gp.commandln("set xlabel 'R (kpc)'");
  gp.commandln("set ylabel 'z (kpc)'");
  gp.commandln("plot 'accr_orbit.dat' using 2:4 with poi 0");
  gp.commandln("set xlabel 'R (kpc)'");
  gp.commandln("set ylabel 'v_R (km/s)'");
  gp.commandln("set origin 0.66,0.5");
  gp.commandln("plot 'accr_orbit.dat' using 2:5 with poi 0");
  gp.commandln("set xlabel 'z (kpc)'");
  gp.commandln("set ylabel 'v_z (km/s)'");
  gp.commandln("set origin 0.66,0.25");
  gp.commandln("plot 'accr_orbit.dat' using 4:7 with poi 0");
  //  gp.commandln("set yrange [0.98:1.02]");
  gp.commandln("set xlabel 't (Gyr)'");
  gp.commandln("set ylabel 'E_{tot} (kpc/Myr)^2'");
  gp.commandln("set origin 0,0");
  //  gp.commandln("set yrange [0.99:1.01]");
  gp.commandln("set autosca");
  gp.commandln("plot 'accr_orbit.dat' using 1:8 with poi 0");
  gp.commandln("set xlabel 't (Gyr)'");
  gp.commandln("set ylabel 'L_{tot} (kpc^2/Myr)'");
  gp.commandln("set origin 0.33,0");
  gp.commandln("set auto y");
  //  gp.commandln("set yrange [0.8:1.2]");
  gp.commandln("plot 'accr_orbit.dat' using 1:9 with poi 0");
  gp.commandln("set xlabel 't (Gyr)'");
  gp.commandln("set ylabel 'L_z (kpc^2/Myr)'");
  gp.commandln("set origin 0.66,0");
  gp.commandln("set auto y");
  //  gp.commandln("set yrange [0.9999:1.0001]");
  gp.commandln("plot 'accr_orbit.dat' using 1:10 with poi 0");
  // panel 2
  gp.commandln("set multiplot");
  gp.commandln("set nokey");
  gp.commandln("set size 0.5,0.25");
  gp.commandln("set origin 0,0.75");
  gp.commandln("set xlabel 't (Gyr)'");
  gp.commandln("set ylabel 'R (kpc)'");
  gp.commandln("plot 'accr_orbit.dat' using 1:2 with poi 0");
  gp.commandln("set ylabel 'phi'");
  gp.commandln("set origin 0,0.5");
  gp.commandln("plot 'accr_orbit.dat' using 1:3 with poi 0");
  gp.commandln("set origin 0,0.25");
  gp.commandln("set ylabel 'z (kpc)'");
  gp.commandln("plot 'accr_orbit.dat' using 1:4 with poi 0");
  gp.commandln("set origin 0.5,0.75");
  gp.commandln("set ylabel 'v_R (km/s)'");
  gp.commandln("plot 'accr_orbit.dat' using 1:5 with poi 0");
  gp.commandln("set origin 0.5,0.5");
  gp.commandln("set ylabel 'v_{phi} (km/s)'");
  gp.commandln("plot 'accr_orbit.dat' using 1:6 with poi 0");
  gp.commandln("set origin 0.5,0.25");
  gp.commandln("set ylabel 'v_z (km/s)'");
  gp.commandln("plot 'accr_orbit.dat' using 1:7 with poi 0");
  gp.commandln("set origin 0,0");
  gp.commandln("set xlabel 'R (kpc)'");
  gp.commandln("set ylabel 'z (kpc)'");
  gp.commandln("plot 'accr_orbit.dat' using 2:4 with poi 0");
  /*
  gp.commandln("set xlabel 'R (kpc)'");
  gp.commandln("set ylabel 'v_R (km/s)'");
  gp.commandln("set origin 0.66,0.5");
  gp.commandln("plot 'accr_orbit.dat' using 2:5 with poi 0");
  gp.commandln("set xlabel 'z (kpc)'");
  gp.commandln("set ylabel 'v_z (km/s)'");
  gp.commandln("set origin 0.66,0.25");
  gp.commandln("plot 'accr_orbit.dat' using 4:7 with poi 0");
  //  gp.commandln("set yrange [0.98:1.02]");
  gp.commandln("set xlabel 't (Gyr)'");
  gp.commandln("set ylabel 'E_{tot} (kpc/Myr)^2'");
  gp.commandln("set origin 0,0");
  //  gp.commandln("set yrange [0.99:1.01]");
  gp.commandln("set autosca");
  gp.commandln("plot 'accr_orbit.dat' using 1:8 with poi 0");
  gp.commandln("set xlabel 't (Gyr)'");
  gp.commandln("set ylabel 'L_{tot} (kpc^2/Myr)'");
  gp.commandln("set origin 0.33,0");
  gp.commandln("set auto y");
  //  gp.commandln("set yrange [0.8:1.2]");
  gp.commandln("plot 'accr_orbit.dat' using 1:9 with poi 0");
  gp.commandln("set xlabel 't (Gyr)'");
  gp.commandln("set ylabel 'L_z (kpc^2/Myr)'");
  gp.commandln("set origin 0.66,0");
  gp.commandln("set auto y");
  //  gp.commandln("set yrange [0.9999:1.0001]");
  gp.commandln("plot 'accr_orbit.dat' using 1:10 with poi 0");
  */
  gp.end();
}
