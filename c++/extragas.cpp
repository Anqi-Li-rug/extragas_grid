/* extragas.cpp
 
 NAME : EXTRAGAS
 
 PURPOSE : integration of orbits of particles (extraplanar gas clouds) 
 in a galactic potential and creation of a particle cube for comparison 
 with HI data
 
 VERSION : 2.0 
 
 STRUCTURE:
 
 exgas : main program
 
 headers: 
 extragas.h:
 array.h:
 dynamics.h:
 projections.h:
 iterations.h:
 numint.h, statist.h, interpol.h(in numint.h) : from NR
 fits_io.h:
 hothalo.h:
 inout.h:
 
 modules:
 dynamics.cpp:
 projections.cpp:
 extragas_functions.cpp:
 extragas_utils.cpp:
 extragas_plots.cpp:
 hothalo.cpp:
 
 
 OUTLINE:
 
 declarations of variables, files...
 declarations and reading of parameters
 declarations of output matrices
 reading forces and potentials
 
 >0 iteration cycle
 defining iteration parameters
 building the hot halo
 setting initial paramters
 
 >1 cycle in R
 settings
 
 >2  cycle in v_k
 settings
 initial conditions
 distribution of kick velocities
 initial energy and L
 
 >3   cycle in time
 settings
 selection of v_k above threshold
 check if the particle is inside the potential grid
 ONE STEP R-K INTEGRATION 
 check integration step
 energy and angular momenta
 building the cubes 
 
 >4    cycle in phi
 settings
 
 >5     cycle in z
 external and internal projections
 calculating the outflow rate
 <5     end z
 
 <4    end phi
 
 <3   end time
 outflow/inflow rate
 transfer of angular momentum
 plotting orbit
 
 <2  end v_k
 
 <1 end R
 
 */
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <stdio.h>
#include <time.h>
#include <cstring>
#include <sstream>
//#include<Magick++.h>
#include "array.h"
#include "arrays3d.h"
#include "numint.h"
#include "statist.h"
#include "findValue.h"
#include "galFunctions.h"
#include "projections.h"
#include "dynamics.h"
#include "fitsUtils.h"
#include "hothalo.h"
#include "accretion.h"
#include "inout.h"
#include "extragasPlots.h"
#include "extragasFits.h"
#include <sys/types.h>
#include <sys/resource.h>
#include "extragas.h"
#include "iterations.h" // keep it after fits_io.h
#include "fitsio.h" 
#include "milkyway.h"
#include "fitsutil.h"
#include "myutilities.h"
#include "profiles.h"
#include <vector>
#include <fstream>
#include <string>

#define TOKMS 1.e-3
//using namespace Magick;

/*
 DECLARATIONS: CONSTANTS
 */
const double PI=3.14159265359;
const double RAD=180./PI;
const double CONV1=1.02269012107e-3; // km/s -> kpc/Myr
const double CONV2=1.98892e43; // Mo*km^2/s^2 -> erg
const double YR=3.155691e7; // yr -> sec
const double MO=1.989e33; // Mo (g)
const double PC=3.08567802e18; // pc (cm)
const double DEGTORAD=0.01745329252;

inline double Rtod_logspiral (double R, double pitch_deg)
{
	return R/cos(atan(1./tan(pitch_deg*DEGTORAD)));
}
inline double dtoR_logspiral (double d, double pitch_deg)
{
	return d*cos(atan(1./tan(pitch_deg*DEGTORAD)));
}

/*
 GLOBAL VARIABLES
 */

//one seed to rule them all
long int seed = time(NULL);

extern char inputfile[];
extern bool verbose; 
Array1D<float> FunctionPar(10);


//initialise surface density HI and H2
std::vector<float> RHI,DHI,RH2,DH2,RSFR,DSFR;
//surface density threshold to smooth the kenninicut law (prop to 1-(threshold/gas)^2)
float Kennicut_T, gammaSF;

inline double gauss2D(float l, float b, float sl, float sb)
{
	return exp(-l*l/(2.*sl*sl))*exp(-b*b/(2.*sb*sb))/(2.*PI*sl*sb);
};

inline double R_logspiral(double theta_deg, double pitch_deg, double A)
{
	return A*exp(tan(pitch_deg*DEGTORAD)*theta_deg*DEGTORAD);
}

inline double Phi_logspiral(double R_kpc, double pitch_deg, double A)
{
	return log(R_kpc/A)/(tan(pitch_deg*DEGTORAD))/DEGTORAD;
}

inline double z0_warp(double R, double phi, bool create_warp)
{
	if(create_warp)
	{
		double W0 = R>15? -66+150*(R-15)-0.47*(R-15)*(R-15) : -66;
		double W1 = R>10?   9+197*(R-10) -3.1*(R-10)*(R-10) :   9;
		double W2 = R>15? -70+171*(R-15) -5.3*(R-15)*(R-15) : -70;
		double phi1 = R>10? 88.49 + pow(R-10.,10.94)*exp(-R/0.884425) : 88.49;
		double phi2 = R>15? 0.1217*R*R - 7.1470*R + 181.89 : 102.0675;
		phi1 -= 90; phi2 += 45;
		return 1e-3*(W0 + W1*sin(((phi-180)-phi1)*DEGTORAD) + W2*sin((2*(phi-180)-phi2)*DEGTORAD));
	}
	else {return 0;}
}

Array3D<float> extragas_model(double v_ki,double ion_i,double accre_i, double &HImass, double &accretion)
{
	//cout.precision(6);
    //fclose(stdout); //close printf
	/* 
	 DECLARATIONS: VARIABLES
	 */
	// Cylindircal Coordinates
	double z, phi, v_R, v_phi, v_z;
	// Absolute Coordinates
	double x1,x2,x3,x4,x5,x6,xx2; 
	// Old Coordinates (used for refining integration time)
	double x1_prev, x2_prev, x3_prev, x4_prev, x5_prev, x6_prev, x_prev, v_prev, delta_x, delta_v, n_refin;
	// Forces and potential
	double previous_fx1,fx1,fx3,pot,fx4,fx6;
	// Elements in inner and outer matrices for forces
	int nR, nz, nR1, nz1, nR2, nz2, cnR, cnz; 
	// Grids in inner and outer matrices
	double stepR1, stepz1, stepR2, stepz2; 
	// Radii
	double R1, R2, R_ave, Rmin_in, R1_in, R2_in;
	// Kick velocities
	float v_k, theta_k, phi_k, mean_v_k=0;
	// Conversion factors
	double pixtopc, pctoJy, ntomass, JytoK;
	// Boolean variables
	bool second_passage=false, integrate_until_tmax=false, keep_integr=false, below=false, out_boundaries=false, integr_decr=false;
	// Times
	float t_max, tau, dt_default, t, t_above, tau_above, n_tau_above;;
	int n_tau;
	double counts_HI_above; //number of HI clouds at |z|>0.2 kpc 
	double counts_halo_above; //number of HI+HII clouds at |z|>0.2 kpc
	double noise;
	double counts_accreted_above;
	// Hot halo
	int xc, yc, vc, np_vk, np_halo, cp_halo, cp_vk, np_phi, np_disc, cp_disc, i_in, i_out, xc_z, pxc, pyc, pvc, ppxc; 
	int already_above, above;
	float Lz_hot_tot;
	//
	float xi,yi,vrad,l,b,d,vlos;
	float vlos_max,vlos_min; // used for HVCs
	float R_out_above, R_in_above; 
	
	/* 
	 DECLARATIONS: OUTPUTS
	 */
	
	double M_halo_above, M_HI_above, M_accreted_above, np_fraction, counts_HI, counts_disc;
	// Total energy and angular momentum per unit mass
	double Etot, Ltot, L_z, dL_z;
	double Etot_M, Ltot_M, L_z_M, Etot_M0, Ltot_M0, L_z_M0;
	double Etot0, Ltot0, L_z0; 
	float cn_dot, cn_disc, cn_accr, cn_accr_prev;
	double total_energy_above, ctotal_energy_above;
	// Aitoff projection
	float aitoff_z, aitoff_l, aitoff_b;
	int aitoff_xc, aitoff_yc;
	float blue, lightblue, green, yellow;
	int R_c, z_c, v_phi_c, v_R_c;
	
	//VELOCITY DISTRIBUTION OF THE PRISTINE MATERIAL (and other stuff)
	double myRmin = 0;
	double myRmax = 20.;
	double delta_myR = 0.25;
	
	double vphi_min = 100;
	double vphi_max = 250;
	double delta_vphi = 5.;
	
	double vR_min = -60;
	double vR_max = 60;
	double delta_vR = 5.;
	
	const int nmyR = (myRmax - myRmin)/delta_myR + 1;
	const int nvphi = (vphi_max - vphi_min)/delta_vphi + 1;
	const int nvR = (vR_max - vR_min)/delta_vR + 1;
	
	double myR[nmyR], MVrot[nmyR], MVrot_ratio[nmyR], M[nmyR], npart[nmyR];
	double vR[nvR], m_vR[nvR];
	double vphi[nvphi], m_vphi[nvphi];
	
	myR[0]=myRmin;
	for(int i=1;i<nmyR;i++) myR[i]=myR[i-1]+delta_myR;
	for(int i=0;i<nmyR;i++) {MVrot[i]=0; M[i]=0; MVrot_ratio[i]=0; npart[i]=0;}
	vR[0]=vR_min;
	for(int i=1;i<nvR;i++) vR[i]=vR[i-1]+delta_vR;
	for(int i=0;i<nvR;i++) m_vR[i]=0.;
	vphi[0]=vphi_min;
	for(int i=1;i<nvphi;i++) vphi[i]=vphi[i-1]+delta_vphi;
	for(int i=0;i<nvphi;i++) m_vphi[i]=0.;
	
	/* 
	 DECLARATIONS: STRINGS
	 */
	char par1[40], par2[40], par3[40];
	
	///* FOR MOVIE
	 string imageName;
	 //list<Image> imageList; 
	 int movieIncr=0;
	 //*/
	
	/*
	 DECLARATIONS: FILES 
	 */
	FILE* f_z_profile1=NULL;
	FILE* f_z_profile2=NULL;
	FILE* pinput_file;
	FILE* f_orbit=NULL;
	FILE* f_veldistr=NULL;
	FILE* f_veldistr2=NULL;
	FILE* f_times=NULL;
    FILE* f_param=NULL;
	FILE* f_hvcs_0=NULL;
	FILE* f_hvcs_p1=NULL;
	FILE* f_hvcs_n1=NULL;
	FILE* f_hvcs_p2=NULL;
	FILE* f_hvcs_n2=NULL;
	FILE* f_hvcs_axes=NULL;
	FILE * f_matrix1;
	FILE * f_matrix2;
	FILE* f_xspec=NULL;
	
	if(verbose)
	{
		cout << "\n";
		cout << "EXTRAGAS\n";
		cout << "Integrating orbits of extraplanar gas particles in a galactic potential\n";
		cout << "\n";
		cout << "Reading parameters from file 'extragas.in'..." << endl;
	}
	
	/* 
	 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	 READING PARAMETERS from INPUT FILE 
	 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	 */
    
//if 'y', the code will largerly use random number to sample many things
//if 'n', the code will regolarize most parts.
    char use_random; findValue(inputfile, use_random, "use_random");
    
	/* Input Files */
	char filename1[90], filename2[90], surfdens_HI[90], surfdens_H2[90], sfr_file[90],datacube_file[90], labfile[90], resfile[90], compare_to[90], absorption_list[90],param_name[90];
	findValue(inputfile, filename1, "inputfile1");
	findValue(inputfile, filename2, "inputfile2");
	findValue(inputfile, resfile, "residuals_file");
	findValue(inputfile, surfdens_HI, "surfdens_HI");
	findValue(inputfile, surfdens_H2, "surfdens_H2");
    findValue(inputfile, sfr_file, "sfr_file");
	findValue(inputfile, datacube_file, "datacube");
	//findValue(inputfile, outputfilename, "outputfilename");
	findValue(inputfile, compare_to, "compare_to");
	findValue(inputfile, absorption_list, "absorption_list");
	int low_radii; findValue(inputfile, low_radii, "low_radii"); 
	int subdim1; findValue(inputfile, subdim1, "subdim1"); 
	int subdim2; findValue(inputfile, subdim2, "subdim2"); 
    findValue(inputfile,param_name,"param_name");
    std::ifstream  src(inputfile, std::ios::binary);
    std::ofstream  dst(param_name,   std::ios::binary);
    dst << src.rdbuf();
	/*
	 GALACTIC PARAMETERS
	 */
	
	//reading HI and H2 surface density profile from files
	int lenght;
	char* buffer;
	int number_of_raws=0;
        cout<<"begin"<<endl;	
	//HI
	//if(verbose)cout<<"Reading HI surface density from file: "<<surfdens_HI<<endl;
	ifstream HIdata(surfdens_HI);
	HIdata.seekg(0, ios::end);
	lenght = HIdata.tellg();
	HIdata.seekg(0, ios::beg);
	buffer = new char [lenght];
	HIdata.read (buffer, lenght);
	for(int i=0; i<lenght; i++) if(buffer[i]==0x0A) number_of_raws++;
	HIdata.seekg(0, ios::beg);
	RHI.resize(number_of_raws); DHI.resize(number_of_raws);
	for(int i=0; i<RHI.size();i++) HIdata>>RHI[i]>>DHI[i];
	HIdata.close();
	
	//for(int i=0; i<RHI.size();i++) cout<<RHI[i]<<" "<<DHI[i]<<endl;
	
	//H2
	//if(verbose)cout<<"Reading H2 surface density from file: "<<surfdens_H2<<endl;
	number_of_raws = 0;
	ifstream H2data(surfdens_H2);
	H2data.seekg(0, ios::end);
	lenght = H2data.tellg();
	H2data.seekg(0, ios::beg);
	buffer = new char [lenght];
	H2data.read (buffer, lenght);
	for(int i=0; i<lenght; i++) if(buffer[i]==0x0A) number_of_raws++;
	H2data.seekg(0, ios::beg);
	RH2.resize(number_of_raws); DH2.resize(number_of_raws);
	for(int i=0; i<RH2.size();i++) H2data>>RH2[i]>>DH2[i];
	H2data.close();
	//SFR
    number_of_raws = 0;
    ifstream SFRdata(sfr_file);
    SFRdata.seekg(0, ios::end);
    lenght = SFRdata.tellg();
    SFRdata.seekg(0, ios::beg);
    buffer = new char [lenght];
    SFRdata.read (buffer, lenght);
    for(int i=0; i<lenght; i++) if(buffer[i]==0x0A) number_of_raws++;
    SFRdata.seekg(0, ios::beg);
    RSFR.resize(number_of_raws); DSFR.resize(number_of_raws);
    for(int i=0; i<RSFR.size();i++) SFRdata>>RSFR[i]>>DSFR[i];
    SFRdata.close();

	/* threshold density (in Msun/pc^2) applied to kennicut law (prop to 1-(threshold/gas)^2) */
	findValue(inputfile, Kennicut_T, "Kennicut_T");
	/* parameters for HoleExp and HoleExpSchmidt Note: they do not change */
	findValue(inputfile, FunctionPar[0], "ReSF"); 
	findValue(inputfile, FunctionPar[1], "alphaSF"); 
	findValue(inputfile, FunctionPar[2], "gammaSF");
	gammaSF = FunctionPar[2];
	/*
	 ofstream mydata("/Users/antoninomarasco/Desktop/SFRDtest.txt");
	 float myR=0;
	 while (myR<20.)
	 {
	 mydata<<myR<<" "<<HoleExpSchmidt(myR)<<endl;
	 myR+=0.1;
	 }
	 mydata.close();
	 */ 
	/* Gas scale length to estimate SFR */
	float alphaSF; findValue(inputfile, alphaSF, "alphaSF"); 
	float ReSF; findValue(inputfile, ReSF, "ReSF"); 
	/* Power to gas density to estimate SFR */
	float gammaSF; findValue(inputfile, gammaSF, "gammaSF"); 
	/* Inclination */
	float incl; findValue(inputfile, incl, "incl");
    float PA; findValue(inputfile, PA, "PA");
	/* Mirror? */
	char mirror_ra; findValue(inputfile, mirror_ra, "mirror_ra"); 
	char mirror_dec; findValue(inputfile, mirror_dec, "mirror_dec"); 
	/* Position of the Sun (for MW projection) */
	float R_Sun; findValue(inputfile, R_Sun, "R_Sun"); 
	/* Cloud size (in pc) to decide whether they are resolved or not, inner view*/
	float cloudsize; findValue(inputfile, cloudsize, "cloudsize"); 
	float d_resolved;
	int nPixels;
	long pixsize1, pixsize2;
	/* Flat velocity (for logaritmic potential) */
	float v0_phi_HVCs; findValue(inputfile, v0_phi_HVCs, "v0_phi_HVCs"); 
	/* Rotation velocity at the sum */
	float v0_phi_Rsun; findValue(inputfile, v0_phi_Rsun, "v0_phi_Rsun"); 
	/* Distance from observer */
	float dist; findValue(inputfile, dist, "distance");
	/* conversion between kpc and degrees */
	double kpctograd=1./dist/1.e3*RAD;
	/* Radial range for fountain */
	double RminSF; findValue(inputfile, RminSF, "RminSF"); 
	double RmaxSF; findValue(inputfile, RmaxSF, "RmaxSF"); 
	double deltaR; findValue(inputfile, deltaR, "deltaR"); 
	double passage_Rmax; findValue(inputfile, passage_Rmax, "passage_Rmax"); 
	//azimuthal range for the fountain
	double PhiminSF=-1e10; double PhimaxSF=1e10;
	/* Starting value for R */
	float R=RminSF;
	/* Radial range for galactic disc */
	float Rmax_disc; findValue(inputfile, Rmax_disc, "Rmax_disk"); 
	float Rmin_disc; findValue(inputfile, Rmin_disc, "Rmin_disk"); 
	int nRSF=(int)((RmaxSF-RminSF)/deltaR)+1;
	/* Total HI mass */
	float M_HI_tot; findValue(inputfile, M_HI_tot, "M_HI_tot");
	/* Total mass of the HI halo */
	float M_HI_halo; findValue(inputfile, M_HI_halo, "M_HI_halo"); 
	/* Mass of the HI disc */
	float M_disc = M_HI_tot - M_HI_halo;
	/* Total HI mass corrected for the primary beam attenuation */
	float M_tot_pbc = M_HI_tot; 
	float pbc=M_tot_pbc/M_HI_tot;

	/* SPIRAL ARMS */
	//if 'y', it creates up to 4 separate logarithmic spiral arms in both disk and halo
	//NOTE: you can remove single arms by using a_overderns<=1
	char create_arms; findValue(inputfile, create_arms, "create_arms");
	//Starting radius (in kpc) of all the arms
	float Rstart_arms; findValue(inputfile, Rstart_arms, "Rstart_arms");
	//if 'y', implement the velocity gap between the fountain cloud and the spiral arm rotation
	char velocity_gap; findValue(inputfile, velocity_gap, "velocity_gap");
	//pattern speed of the spiral arms (in km/s/kpc)
	float pattern_speed; findValue(inputfile, pattern_speed, "pattern_speed");
	//amplitudes of the logarithmic spiral arms
	float a1_amp; findValue(inputfile, a1_amp, "a1_amp");
	float a2_amp; findValue(inputfile, a2_amp, "a2_amp");
	float a3_amp; findValue(inputfile, a3_amp, "a3_amp");
	float a4_amp; findValue(inputfile, a4_amp, "a4_amp");
	//pitch angles of the logarithmic spiral arms
	float a1_pitch; findValue(inputfile, a1_pitch, "a1_pitch");
	float a2_pitch; findValue(inputfile, a2_pitch, "a2_pitch");
	float a3_pitch; findValue(inputfile, a3_pitch, "a3_pitch");
	float a4_pitch; findValue(inputfile, a4_pitch, "a4_pitch");
	//local contrast density of the arms with respect to the axysimmetric disc
	float a1_overdens; findValue(inputfile, a1_overdens, "a1_overdens");
	float a2_overdens; findValue(inputfile, a2_overdens, "a2_overdens");
	float a3_overdens; findValue(inputfile, a3_overdens, "a3_overdens");
	float a4_overdens; findValue(inputfile, a4_overdens, "a4_overdens");
	if(a1_overdens<=1.) a1_overdens=1.; if(a2_overdens<=1.) a2_overdens=1.;if(a3_overdens<=1.) a3_overdens=1.;if(a4_overdens<=1.) a4_overdens=1.;
	//radial extent of the arms (FWHM, in kpc)
	float a1_thickness; findValue(inputfile, a1_thickness, "a1_thickness");
	float a2_thickness; findValue(inputfile, a2_thickness, "a2_thickness");
	float a3_thickness; findValue(inputfile, a3_thickness, "a3_thickness");
	float a4_thickness; findValue(inputfile, a4_thickness, "a4_thickness");
	a1_thickness /= 2.355; a2_thickness /= 2.355; a3_thickness /= 2.355; a4_thickness /= 2.355; //convert to rms
	//number of "passage" intercepted at a given theta
	float minpitch1 = std::min(a1_pitch,a2_pitch);
	float minpitch2 = std::min(a3_pitch,a4_pitch);
	float minpitch = std::min(minpitch1,minpitch2);
	int kk = 0;
	while(R_logspiral(kk*360, minpitch, 1.0)<=Rmax_disc) {kk++;}
	const int n_R_arms = kk+1;
	double R_overdensity1[n_R_arms],R_overdensity2[n_R_arms],R_overdensity3[n_R_arms],R_overdensity4[n_R_arms];
	double arms_weight;
	//weights for normalization
	double Weight_noarms = 0;
	double Weight_witharms = 0;

	//WARP (Levine+06)
	char include_warp; findValue(inputfile, include_warp, "create_warp");
	bool create_warp = false; 
	if(include_warp=='y') create_warp = true;
	
	/*
	 PARAMETERS OF THE CUBE
	 */
	
	char projectionType[10]; findValue(inputfile, projectionType, "projectionType");
	if (strcoll(projectionType,"external") !=0) kpctograd=1.;
	/* Spatial grid and size of the data cube */
	float cdelt1_data; findValue(inputfile, cdelt1_data, "cdelt1_data"); 
	// in degree
	float cdelt2_data; findValue(inputfile, cdelt2_data, "cdelt2_data"); 
	// in degree
	int RAsize_data; findValue(inputfile, RAsize_data, "RAsize_data"); 
	int DECsize_data; findValue(inputfile, DECsize_data, "DECsize_data"); 
	/* Scale of the grid in spatial axes */
	float scale; findValue(inputfile, scale, "scale"); 
	/* Size of the model cube (pixel) */
	int RAsize=(int)(RAsize_data/scale+0.5);
	int DECsize=(int)(DECsize_data/scale+0.5);
	/* Spatial grid of the model cube (degrees) */
	float cdelt1;
	// in degree
	float cdelt2;
	// in degree
	if (strcoll(projectionType,"external") ==0) {
		cdelt1=cdelt1_data*scale;
		cdelt2=cdelt2_data*scale;
	}
	else {
		cdelt1=-360./(RAsize-1); // N.B. in degrees
		cdelt2=180./(DECsize-1); // N.B. in degrees
	}
	/* Grid of the model cube for internal projection */
	float ldelt=cdelt1_data*scale;
	float bdelt=cdelt2_data*scale;
	/* Beam size (resolution) of the data */
	float beam_x; findValue(inputfile, beam_x, "beam_x"); 
	// in degree
	float beam_y; findValue(inputfile, beam_y, "beam_y"); 
	// in degree
	/* Shift of galactic centre with respect to the centre of the cube */
	float x_shift; findValue(inputfile, x_shift, "x_shift"); 
	float y_shift; findValue(inputfile, y_shift, "y_shift");
	/* Velocity grid of the data cube */
	float cdelt3_data; findValue(inputfile, cdelt3_data, "cdelt3_data"); 
	/* Velocity grid and size of the model cube */
	float cdelt3; findValue(inputfile, cdelt3, "chansep"); 
	// NNBB: in the *.in file it cannot be called cdelt3 because it would read
	//       cdelt3_data!
	int VELsize; findValue(inputfile, VELsize, "VELsize"); 
	float scale3=fabs(cdelt3/cdelt3_data);
	/* Channel number in the model that corresponds */
	float channel0; findValue(inputfile, channel0, "channel0"); 
	/* Systemic velocity and shift in velocity */
	float v_sys; findValue(inputfile, v_sys, "v_sys"); 
	double crval1; findValue(inputfile, crval1, "CRVAL1"); 
	double crval2; findValue(inputfile, crval2, "CRVAL2"); 
	double crval3; findValue(inputfile, crval3, "CRVAL3"); 
	crval3=crval3-cdelt3*(channel0);
	float v_shift=v_sys-crval3;
	/* Cut off distance for closeby HVCs */
	float distance_cut; findValue(inputfile, distance_cut, "dist_cut"); 
	/*  thing */
	char hanning; findValue(inputfile, hanning, "hanning");
	//  char smooth; findValue(inputfile, &smooth, "smooth");
	char smooth; findValue(inputfile, smooth, "smooth");
	//  Conversion from Jy/beam into K?
	char Jy2K; findValue(inputfile, Jy2K, "Jy2K");
	/* Adding disc component */
	char add_disk; findValue(inputfile, add_disk, "add_disk");
	//scale height of the not-flared part of the disc
	float h_disk; findValue(inputfile, h_disk, "h_disk"); 
	float ch_disk;
	//The flare stars at R>Rstart_flare
	float Rstart_flare; findValue(inputfile, Rstart_flare, "Rstart_flare"); 
	//scale-radius of the flare
	float Rs_flare; findValue(inputfile, Rs_flare, "Rs_flare"); 
	//consider only positive (p) or negative (n) latitudes if you wish
	char which_latitude; findValue(inputfile, which_latitude, "which_latitude");

	/* 
	 PLOTTING PARAMETERS
	 */
	/* Starting parameters for plotting a single orbit */
	float R_plot; findValue(inputfile, R_plot, "R_plot"); 
	float v_plot; findValue(inputfile, v_plot, "v_plot"); 
	float deltav_plot; findValue(inputfile, deltav_plot, "deltav_plot");
	float deltaR_plot; findValue(inputfile, deltaR_plot, "deltaR_plot");
	/* Position for plotting the vertical profiles */
	float R_z_plot1; findValue(inputfile, R_z_plot1, "R_z_plot1");
	float R_z_plot2; findValue(inputfile, R_z_plot2, "R_z_plot2");
	/* Grid for plotting rotation curves at different heights */
	float rotsurf_grid; findValue(inputfile, rotsurf_grid, "rotsurf_grid");
	/* Clip to calculated weithed mean velocity (in percentage of the peak) */
	float clip_wm; findValue(inputfile, clip_wm, "clip_wm");
	float Rmin_gradient; findValue(inputfile, Rmin_gradient, "Rmin_gradient");
	float Rmax_gradient; findValue(inputfile, Rmax_gradient, "Rmax_gradient");
	float zmax_gradient; findValue(inputfile, zmax_gradient, "zmax_gradient");
	/* For HVCs */
	char v_hvcs_plot[10]; findValue(inputfile, v_hvcs_plot, "v_hvcs_plot");
	float vlos_cut; findValue(inputfile, vlos_cut, "vlos_cut");
	float zmax_hvcs; findValue(inputfile, zmax_hvcs, "zmax_hvcs");
	float Rmax_hvcs; findValue(inputfile, Rmax_hvcs, "Rmax_hvcs");
	char make_fit; findValue(inputfile, make_fit, "make_fit");
	char make_movie; findValue(inputfile, make_movie, "make_movie");
	
	float delta_param1, delta_param2, delta_param3, param1_0, param2_0, param3_0;
	char iter_mode; findValue(inputfile, iter_mode, "iter_mode"); 
	char projAC; findValue(inputfile, projAC, "projAC");

	/* 
	 PARAMETERS for MINIMIZATION OF the RESIDUALS
	*/
	float threshold_res; findValue(inputfile, threshold_res, "threshold_res");//It defines the "residual calculation zone" in the Milky Way mode
	char disk_cut; findValue(inputfile, disk_cut, "disk_cut");//If you want to delete the disk emission using deviation velocity
	/* Calculate residuals? */
		char normalize_flux; findValue(inputfile, normalize_flux, "normalize_flux");//if 'y' the model flux is normalized to data flux before residuals calculation.	
	/* Distance above the plane above which compare with the data */
	float z_above_in; findValue(inputfile, z_above_in, "z_above_in"); //ejecting clouds
	float z_above_fin; findValue(inputfile, z_above_fin, "z_above_fin");//returning clouds
	char stop;

	/*
	 INITIAL CONDITIONS AND INTEGRATION PARAMETERS
	 */
	/* Number of particles per bin (before cutting low velocities) */
	int n_part; findValue(inputfile, n_part, "n_part"); 
	//  int nRSF=n_part; POSSIBLE?
	/* Lower threshold in kick velocities */
	float v_k_thres; findValue(inputfile, v_k_thres, "v_k_thres"); 
	/* Characteristic kick velocity */
	//float h_v_k; findValue(inputfile_temp, h_v_k, "h_v_k");
    float h_v_k=v_ki;
	/* Central kick velocity (for gaussian distribution) */
	float v_k_0; findValue(inputfile, v_k_0, "v_k_0"); 
	/* Exponent for power law distribution of kick velocities */
	float alpha_v_k; findValue(inputfile, alpha_v_k, "alpha_v_k"); 
	/* Tilt for linear distribution of kick velocities */
	float ab_v_k; findValue(inputfile, ab_v_k, "ab_v_k"); 
	/* Type of potential: 1) logarithmic 2) matrix */
	//  char type_pot[10]; findValue(inputfile, type_pot, "type_pot"); 
	/* Intergration step */

	float dt; findValue(inputfile, dt, "delta_t"); 
	/* Default integration step */
	dt_default=dt;
	/* Maximum number of integrations */
	int n_dt_max; findValue(inputfile, n_dt_max, "n_dt_max"); 
	/* Type of vertical distribution: gau, pow, lin... */
	char v_k_distr[10]; findValue(inputfile, v_k_distr, "v_k_distr"); 
	/* Maximum kick velocity */
	float v_k_max; findValue(inputfile, v_k_max, "v_k_max");
	/* Opening angle for kicks (in degree) */
	float kick_angle; findValue(inputfile, kick_angle, "kick_angle");

	/* Exponent of the opening angle for kick velocities (it decreases the particle velocity) */
	float alpha_cos;findValue(inputfile, alpha_cos, "alpha_cos");
	/* Random velocity dispersion (if dsigma_v =! 0 then this is the dispersion at Rmax/2) */
	float sigma_v;findValue(inputfile, sigma_v, "sigma_v"); 
	/* Change in velocity dispersion from the inner to the outer parts */
	float dsigma_v;findValue(inputfile, dsigma_v, "dsigma_v");
	float csigma_v; // current sigma, for each R (see below) 
	/* Vertical dispersion (disc scaleheight) */
	float sigma_z;findValue(inputfile, sigma_z, "sigma_z");

	/* End of integration: first passage, second or number of Myr */
	char integr_end[10]; findValue(inputfile, integr_end, "integr_end");
	/* Visible part of the orbit: all, from apocentre, from top za fraction of the rising orbit */
	char visible[10]; findValue(inputfile, visible, "visible_part"); // option are: all, up, down, R_max, partial
	//fraction of the rising orbit (assumed 1/2 of the full one) not visible. A linear run of the vertical velocity is assumed
	//float ion_frac; findValue(inputfile_temp, ion_frac, "ion_frac");
    float ion_frac=ion_i;
	/* Parameters for refining the integration */

	float delta_x_max; findValue(inputfile, delta_x_max, "delta_x_max"); 
	float delta_v_max; findValue(inputfile, delta_v_max, "delta_v_max"); 
	float delta_x_min; findValue(inputfile, delta_x_min, "delta_x_min"); 
	float delta_v_min; findValue(inputfile, delta_v_min, "delta_v_min"); 
	//  float delta_t_min; findValue(inputfile, delta_t_min, "delta_t_min"); 
	int n_refin_max; findValue(inputfile, n_refin_max, "n_refin_max");
	

	/*
	 ISOLATED EVENTS
	*/
	//allow for single (or double) burst in a precise location of an arm
	char single_event; findValue(inputfile, single_event, "single_event");
	//choose one of the four spiral arms 
	int which_arm; findValue(inputfile, which_arm, "which_arm");
	//location of the burst in R along that particular arm
	float R_burst; findValue(inputfile, R_burst, "R_burst"); 
	//spatial extent of the burst (in kpc along the arm around the defined location) 
	float arm_extent; findValue(inputfile, arm_extent, "arm_extent"); 
	double pitch = NAN;
	double amplitude = NAN;
	double thickness = NAN; 
	if(single_event=='y')
	{
		if(which_arm==1)      {a1_overdens=1e8;a2_overdens=1;a3_overdens=1;a4_overdens=1;pitch=a1_pitch;amplitude=a1_amp; thickness=a1_thickness;}
		else if (which_arm==2){a1_overdens=1;a2_overdens=1e8;a3_overdens=1;a4_overdens=1;pitch=a2_pitch;amplitude=a2_amp; thickness=a2_thickness;}
		else if (which_arm==3){a1_overdens=1;a2_overdens=1;a3_overdens=1e8;a4_overdens=1;pitch=a3_pitch;amplitude=a3_amp; thickness=a3_thickness;}
		else if (which_arm==4){a1_overdens=1;a2_overdens=1;a3_overdens=1;a4_overdens=1e8;pitch=a4_pitch;amplitude=a4_amp; thickness=a4_thickness;}
		RminSF = R_burst; 
		deltaR = dtoR_logspiral(Rtod_logspiral(R_burst,pitch)+0.5*arm_extent,pitch)-dtoR_logspiral(Rtod_logspiral(R_burst,pitch)-0.5*arm_extent,pitch);
		RmaxSF = RminSF+deltaR-1e-10;
		nRSF=(int)((RmaxSF-RminSF)/deltaR)+1;
	}
	//lookback time (in Myr) at the beginning of the burst
	float t_burst; findValue(inputfile, t_burst, "t_burst");
	//duration (in Myr) of the first burst. If equal to t_burst, it lasts till now
	float t_duration; findValue(inputfile, t_duration, "t_duration");
	float showt_min = t_burst - t_duration;
	float showt_max = t_burst;
	//second burst, if any
	char second_burst; findValue(inputfile, second_burst, "second_burst");
	float t2_burst; findValue(inputfile, t2_burst, "t2_burst");
	float t2_duration; findValue(inputfile, t2_duration, "t2_duration");
	float showt2_min = t2_burst - t2_duration;
	float showt2_max = t2_burst;

	/*
	 PARAMETERS FOR INTERACTIONS AND ACCRETION
	 */
	char effects[20]; findValue(inputfile, effects, "effects");
	/* Temperature of the hot gas */
	float T_hot; findValue(inputfile, T_hot, "T_hot");
	/* Metallicity */
	float mu;   findValue(inputfile, mu, "mu");
	/* Luminosity of the hot halo */
	double L_hot; findValue(inputfile, L_hot, "L_hot");
	/* Range where luminosity is calculated */
	float R_max_hot; findValue(inputfile, R_max_hot, "R_max_hot"); // radial extent of the hot halo
	float z_min_hot; findValue(inputfile, z_min_hot, "z_min_hot"); 
	float z_max_hot; findValue(inputfile, z_max_hot, "z_max_hot");
	/* Drag paramter K_D=C_D*pi*(D_cloud/2.)*(D_cloud/2.)/m_cloud*mu 
	 C_D= drag constant ~0.5
	 D_cloud= diameter in pc
	 m_cloud= Mass in Mo
	 mu= metallicity
	 */
	float K_D; findValue(inputfile, K_D, "K_D"); 
	/* Rotation velocity of the hot halo (flat part) */
	float v0_hot;   findValue(inputfile, v0_hot, "v0_hot"); 
	
	/* Parameters of accretion */
	//float alpha_accr;  findValue(inputfile_temp, alpha_accr, "alpha_accr");
    float alpha_accr=accre_i;
	float accr_z_in;   findValue(inputfile, accr_z_in, "accr_z_in");
	float accr_z_fin;   findValue(inputfile, accr_z_fin, "accr_z_fin");
	float accr_r0; findValue(inputfile, accr_r0, "accr_r");
	float accr_v_phi0; findValue(inputfile, accr_v_phi0, "accr_v_phi");
	float accr_sigma_v0; findValue(inputfile, accr_sigma_v0, "accr_sigma_v");
	
	/* Type of accretion */
	char accr_dens[20];   findValue(inputfile, accr_dens, "accr_dens");
	/* Velocity pattern for accreting material */
	char accr_vel[20];   findValue(inputfile, accr_vel, "accr_vel");//see dynamic.cpp for different values
	//lag between disk and environment rotation, if "accr_vel" is set to "slowrot" 
	double relative_vel; findValue(inputfile, relative_vel, "relative_vel");
		
	/* Drag parameters */
	double M_cloud; findValue(inputfile, M_cloud, "M_cloud");
	double R_cloud; findValue(inputfile, R_cloud, "R_cloud");
	double rho_corona; findValue(inputfile, rho_corona, "rho_corona");
		
	/* DEFINING OUTPUT MATRIXES */
	Array3D<float> outcube(RAsize,DECsize,VELsize,0.f);
	myarray::double3d datacube;
    myarray::double3d maskcube;
	double* datacube_ptr = NULL;
    double* maskcube_ptr = NULL;
	Array2D<float> outcube_tot(RAsize,DECsize,0.f);
	Array2D<float> maxax(RAsize,VELsize,0.f);
	Array2D<float> tmaxax(RAsize,VELsize,0.f);
	Array2D<float> hvc_upcut(RAsize,DECsize,0.f);
	Array2D<float> hvc_lowcut(RAsize,DECsize,0.f);
	
	//cube axes
	std::vector<double> RA, DEC, VEL;
	for(int i=0;i<RAsize;i++) RA.push_back(crval1+i*cdelt1);
	for(int i=0;i<DECsize;i++) DEC.push_back(crval2+i*cdelt2);
	for(int i=0;i<VELsize;i++) VEL.push_back(crval3+i*cdelt3);
	

	/*
	 INTEGRATION PARAMETERS
	 */
	/* integration ending (first or second passage through the disk) */
	if (strcmp(integr_end, "first0") == 0){
		second_passage=false;
		t_max=1e6;
	} 
	else {
		if (strcmp(integr_end, "second0") == 0){
			second_passage=true;
			t_max=1e6;
		} 
		else {
			integrate_until_tmax=true; // for integration up to t=t_max
			findValue(inputfile, t_max, "integr_end");
			n_dt_max=100000000;
		}
	}
	
	// HVCS ******TEMPORARY******
	
	if (strcoll(projectionType,"hvcs") ==0) {
		f_hvcs_0=fopen("extragas_hvcs_0.dat", "w");
		f_hvcs_p1=fopen("extragas_hvcs_p1.dat", "w");
		f_hvcs_n1=fopen("extragas_hvcs_n1.dat", "w");
		f_hvcs_p2=fopen("extragas_hvcs_p2.dat", "w");
		f_hvcs_n2=fopen("extragas_hvcs_n2.dat", "w");
		f_hvcs_axes=fopen("extragas_hvcs_axes.dat", "w");
		/* building boundary matrices */
		for (int i=0;i<RAsize;i++) {
			for (int j=0;j<DECsize;j++) {
				//	float tl=(i-RAsize/2)*ldelt+180;
				float tl=(i-RAsize/2)*ldelt+180; // 0 < l < 360
				if (tl >= 360) {tl=359.9;}
				if (tl <= 0) {tl=0.1;}
				float tb=(j-DECsize/2)*bdelt;
				/*	
				 for (float tl=-180;tl<180;tl++) {
				 for (float tb=-90;tb<=90;tb++) {
				 matrixcell_int(-tl, tb, 0);
				 if (xc < RAsize && xc >= 0 && yc < DECsize && yc >=0) {
				 */
				projections m;
				hvc_upcut[i][j]=m.modelgal(tl,tb,+1, R_Sun, v0_phi_HVCs, Rmax_hvcs, zmax_hvcs);
				hvc_lowcut[i][j]=m.modelgal(tl,tb,-1, R_Sun, v0_phi_HVCs, Rmax_hvcs, zmax_hvcs);
				//	guihvc[i][j]=hvc_upcut[i][j]+fabs(hvc_lowcut[i][j]);
				//	hvc_lowcut[RAsize-i][j]=-modelgal(tl,tb,-1, R_Sun, v0_phi, Rmax_hvcs, zmax_hvcs)-vlos_cut;
			}
		}
	}

	/*
	 READING FORCES AND POTENTIAL
	 */
	float v0_phi=0.;
	int n = 0, nFiles;
	
	f_matrix1=fopen(filename1,"r");
	if (strcmp(filename2,"none") == 0) {
		//if(verbose)cout << "Using potential from file: " << filename1 << endl;
		nFiles=1;}
	else {
		f_matrix2=fopen (filename2,"r");
		//if(verbose)cout << "Using potential from files: " << filename1 << " and " << filename2 << endl;
		nFiles=2;
	}
	/* DETERMINE INPUT GRID */
	/* FILE 1 (outer matrix) */
	int if_stop=0;
	fscanf (f_matrix1, "%lf %*lf %*lf %*lf %*lf ", &fx1);
	previous_fx1=fx1;
	nz1=1, nR1=1;
	while (!feof(f_matrix1)) {
		fscanf (f_matrix1, "%lf %*lf %*lf %*lf %*lf ", &fx1);
		if (fx1 == previous_fx1 && if_stop==0) {nz1++;}
		if (fx1 != previous_fx1) {
			nR1++;
			if_stop=1;
		}
		previous_fx1=fx1;
	}
	//if(verbose)printf ("Reading %d x %d matrix \n",nR1,nz1);
	inout io(nR1, deltaR);

	/* Defining Arrays for the extended potential */
	Array1D<double> R_i1(nR1), z_k1(nz1);
	Array2D<double> pot_ik1(nR1, nz1), F_R1(nR1, nz1), F_z1(nR1, nz1); 
	Array2D<double> rho_hot_ik1(nR1, nz1, 0.); 
	Array2D<float> dLz_hot_ik1(nR1, nz1, 0.), dvphi_hot_ik1(nR1, nz1, 0.),
	cdLz_hot_ik1(nR1, nz1, 0.), cdvphi_hot_ik1(nR1, nz1, 0.);
	
	/* reading first matrix */
	rewind (f_matrix1);
	for (cnR=0;cnR<nR1;cnR++){
		for (cnz=0;cnz<nz1;cnz++){
			fscanf (f_matrix1, "%lf %lf %lf %lf %lf ", &fx1,&fx3,&pot,&fx4,&fx6);
			R_i1[cnR]=fx1;
			z_k1[cnz]=fx3;
			pot_ik1[cnR][cnz]=pot; // (kpc/Myr)^2
			F_R1[cnR][cnz]=fx4; // kpc/Myr^2
			F_z1[cnR][cnz]=fx6; // kpc/Myr^2
		}
	}
	stepR1=R_i1[1]-R_i1[0];
	stepz1=z_k1[1]-z_k1[0];
	fclose (f_matrix1);
	
	/* FILE 2 (inner matrix) */
	if (nFiles == 2) {
		if_stop=0;		fscanf (f_matrix2, "%lf %lf %lf %lf %lf ", &fx1,&fx3,&pot,&fx4,&fx6);
		previous_fx1=fx1;
		nz2=1;
		nR2=1;
		while (!feof(f_matrix2)) {
			fscanf (f_matrix2, "%lf %*lf %*lf %*lf %*lf ", &fx1);
			if (fx1 == previous_fx1 && if_stop==0) {nz2++;}
			if (fx1 != previous_fx1) {
				nR2++;
				if_stop=1;
			}
			previous_fx1=fx1;
		}
		//if(verbose)printf ("Reading %d x %d matrix \n",nR2,nz2);
	}
	else {nR2=1, nz2=1;}
	
	/* Defining Arrays for the inner potential */
	Array1D<double> R_i2(nR2), z_k2(nz2);
	Array2D<double> pot_ik2(nR2, nz2), F_R2(nR2, nz2), F_z2(nR2, nz2);
	Array2D<double> rho_hot_ik2(nR2, nz2, 0.); 
	//  Array2D<double> dLz_hot_ik2(nR2, nz2, 0.), dvphi_hot_ik2(nR2, nz2, 0.),
	//    cdLz_hot_ik2(nR2, nz2, 0.), cdvphi_hot_ik2(nR2, nz2, 0.);
	
	if (nFiles == 2) {
		rewind (f_matrix2);
		for (cnR=0;cnR<nR2;cnR++){
			for (cnz=0;cnz<nz2;cnz++){
				fscanf (f_matrix2, "%lf %lf %lf %lf %lf ", &fx1,&fx3,&pot,&fx4,&fx6);
				R_i2[cnR]=fx1;
				z_k2[cnz]=fx3;
				pot_ik2[cnR][cnz]=pot;
				F_R2[cnR][cnz]=fx4;
				F_z2[cnR][cnz]=fx6;
			}
		}
		stepR2=R_i2[1]-R_i2[0];
		stepz2=z_k2[1]-z_k2[0];
		fclose (f_matrix2);
	}
	
	/*
	 INITIAL CONDITION
	 */
	x5=0.;
	
	/*
	 ACCRETION PATTERNS
	 */
	Array2D<double> cold_rho_ik(nR1,nz1, 0.f), 
	cold_v_R_ik(nR1,nz1), cold_v_phi_ik(nR1,nz1), cold_v_z_ik(nR1,nz1);
	
	if (strcoll(effects,"accretion")==0 || strcoll(effects,"accrdrag")==0 || strcoll(effects,"accrnodrag")==0){
		if (strcoll(accr_dens,"build") ==0){
			Accretion smooth_accretion(accr_r0, 0, accr_v_phi0, 0, accr_sigma_v0);
			if(verbose)printf("Building accretion pattern \n");
			smooth_accretion.Patterns(cold_rho_ik,
									  cold_v_R_ik,
									  cold_v_phi_ik,
									  cold_v_z_ik,
									  R_i1,z_k1,pot_ik1,F_R1,F_z1);
			
			printHeaderAccr(nR1, nz1, stepR1, stepz1);
			//writefits_2D_DP("accr_density.fits", &cold_rho_ik[0][0], nR1, nz1, "header.tab");
			printHeaderAccr(nR1, nz1, stepR1, stepz1);
			//writefits_2D_DP("accr_v_R.fits", &cold_v_R_ik[0][0], nR1, nz1, "header.tab");
			printHeaderAccr(nR1, nz1, stepR1, stepz1);
			//writefits_2D_DP("accr_v_phi.fits", &cold_v_phi_ik[0][0], nR1, nz1, "header.tab");
			printHeaderAccr(nR1, nz1, stepR1, stepz1);
			//writefits_2D_DP("accr_v_z.fits", &cold_v_z_ik[0][0], nR1, nz1, "header.tab");
		}
	}
	/*
	 ------------------
	 0) ITERATION CYCLE
	 ------------------
	 */
	
	
	stop=n;

	time_t start,end;
	time (&start);	

	int ii=0, jj=0, vkk=0;
	io.back_to_zero();
	/*
	 BUILDING THE HOT HALO
	 */
	if (strcoll(effects,"drag") ==0) {
		HotHalo halo(v0_hot,T_hot,mu,R_max_hot,z_min_hot, z_max_hot, L_hot);
		if (nFiles == 2) {
			halo.SetDensity(rho_hot_ik1,R_i1,z_k1,pot_ik1,F_R1,
							rho_hot_ik2,R_i2,z_k2,pot_ik2,F_R2);
		}
		else {
			halo.SetDensity(rho_hot_ik1,R_i1,z_k1,pot_ik1,F_R1);
		}
	}

	/* 
	 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	 INTEGRATION AND PROJECTIONS
	 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	 */

	/*
	 setting initial parameters
	 */      
	np_halo=0;
	n_tau=0;
	tau=0;
	tau_above = 0;
	n_tau_above = 0;
	R=RminSF;
	counts_HI=0;
	counts_HI_above=0;
	counts_halo_above=0;
	counts_accreted_above=0;
	double counts_halo_local = 0.;
	double counts_accreted_local = 0.;
	total_energy_above=0;
	double n_lost = 0;
	double n_total = 0;
	Lz_hot_tot=0;
	//      Array2D<double> dLz_hot_ik1(nR1, nz1, 0.), dvphi_hot_ik1(nR1, nz1, 0.), tdLz_hot_ik1(nR1, nz1, 0.), tdvphi_hot_ik1(nR1, nz1, 0.);
	/*
	 output files
	 */
	if (iter_mode == 'n'){//replace with n
        //f_orbit=fopen("orbit.dat", "w");
		//f_veldistr=fopen("extragas_veldistr.dat", "w");
		//f_veldistr2=fopen("extragas_veldistr2.dat", "w");
		//f_times=fopen("extragas_time.dat", "w");
		if ((RminSF-deltaR) > 0) {
			/*fprintf(f_times,"0 0 0 0 %f %.2e \n", RminSF-deltaR, 
					2*PI*(RminSF-deltaR)/sqrt
					((RminSF-deltaR)*
					 getValue(F_R1,R_i1,z_k1,(RminSF-deltaR),0.,1))*1e6);*/
		}
	}

	/* 
	 output velocities and densities 
	 */
	int rotsurf_max = 25; // in kpc;
	int rotsurf_nR=(int)(rotsurf_max/rotsurf_grid);
	int rotsurf_nz=(int)(rotsurf_max/rotsurf_grid);
	int rotsurf_nvel=500;
	Array3D<float> out_v_phi(rotsurf_nR,rotsurf_nz,rotsurf_nvel);
	Array3D<float> out_v_R(rotsurf_nR,rotsurf_nz,rotsurf_nvel);
	Array2D<float> out_v_phi_wm(rotsurf_nR,rotsurf_nz, 0.f);
	Array2D<float> out_dens(rotsurf_nR,rotsurf_nz);
	/*if(verbose)
	{
		
		cout<<"h_v_k        = "<<h_v_k<<endl;
		cout<<"ion_frac     = "<<ion_frac<<endl;
		cout<<"alpha_accr   = "<<alpha_accr<<endl;
	}*/
	double avg_R = 0;
	double avg_z = 0;
	double avg_weight = 0;
	/*
	 ---------------------------
	 1) MAIN CYCLE IN R (RADIUS)
	 ---------------------------
	 */
	//if(verbose)cout<<"Starting integration" << endl;
	//      K_D=pow(10,K_D);
	dynamics dyn1(F_R1, 
				  F_z1,
				  pot_ik1,
				  R_i1, 
				  z_k1, 
				  subdim1,
				  subdim2,
				  effects, 
				  K_D, 
				  v0_hot, 
				  rho_hot_ik1,
				  relative_vel,
				  accr_sigma_v0,
				  M_cloud,
				  R_cloud,
				  rho_corona,
				  alpha_accr
				  );
	dynamics dyn2(F_R2, 
				  F_z2, 
				  pot_ik2,
				  R_i2, 
				  z_k2, 
				  subdim1,
				  subdim2,
				  effects,
				  K_D, 
				  v0_hot, 
				  rho_hot_ik2,
				  relative_vel,
				  accr_sigma_v0,
				  M_cloud,
				  R_cloud,
				  rho_corona,
				  alpha_accr
				  );

	dyn1.set_accr_z_in(accr_z_in);
	dyn1.set_accr_z_fin(accr_z_fin);
	dyn2.set_accr_z_in(accr_z_in);
	dyn2.set_accr_z_fin(accr_z_fin);

	while (ii <= nRSF-1) {
		/* settings */
		cp_vk=0;
		cp_halo=0;
		ctotal_energy_above=0;
		//	R=RminSF+ran0(&seed)*(RmaxSF-RminSF);
		R=RminSF+deltaR*ii;
		if(verbose)cout <<"Working at Radius R= "<<R << endl;
		/*
		 calculate number of particles in this bin
		 */
		R1=R-deltaR/2.;
		R2=R+deltaR/2.;
		if(single_event=='y')
		{
			PhiminSF = Phi_logspiral(R1, pitch, amplitude); 
			PhimaxSF = Phi_logspiral(R2, pitch, amplitude);
			while(PhiminSF>360.) PhiminSF-=360.;
			while(PhimaxSF>360.) PhimaxSF-=360.;
			R1-=2.5*thickness; if(R1<0.05) R1=0.05;
			R2+=2.5*thickness; if(R2>Rmax_disc) R2=Rmax_disc;
		}
		R_ave=R;
		np_phi=(int)((2.*PI*R/deltaR+0.5)); //old version
        //np_phi = 3.0*float(RAsize)/beam_x; //3 points per beam
        if(single_event=='y') {np_vk= n_part*arm_extent + 0.5;}
		else {np_vk=n_part*(R2-R1)+0.5;}
		if(np_phi<1) np_phi = 1;
		if(np_vk<1) np_vk = 1;
		if(single_event=='y') {cn_dot=HoleExpSchmidt(R_ave)/qromo(HoleExpSchmidt,R1,R2,'p',1.e-6);}
		else 	{cn_dot=HoleExpSchmidt(R_ave)/qromo(HoleExpSchmidt,RminSF,RmaxSF,'p',1.e-6);} // normalized n_dot
        cn_dot=DSFR[ii];
		/*
		 ---------------------------
		 2) CYCLE IN KICK VELOCITIES
		 ---------------------------
                
		 */
		while (cp_vk < np_vk) {
			if (iter_mode == 'n'){//replace with n
				//freopen("orbit.dat","w",f_orbit);
			}
			/* settings */
			out_boundaries=false;
			already_above=0;
			above=0;
			R_out_above=0.;
			R_in_above=0.;
			t_above=0;
			cn_accr=0;
			
			/*
			 Initial positions
			 */
			/* 1) Initial radius */
            if(use_random=='y')
            {
                R=ran0(&seed)*(R2-R1)+R1;
                if (R <= stepR2){
                    if (low_radii == 1){
                        R=ran0(&seed)*(R2-stepR2)+stepR2;
                    }
                    else {
                        if (R <= R_i2[0]){
                            R=ran0(&seed)*(R2-R_i2[0]*2.)+R_i2[0]*2.;
                        }
                    }
                }
            }
            else {R = R1 + (R2-R1)*float(cp_vk)/float(np_vk);}
			
			if (R <= stepR2){
				if (low_radii == 1){
					R=ran0(&seed)*(R2-stepR2)+stepR2;
				} 
				else {
					if (R <= R_i2[0]){
						R=ran0(&seed)*(R2-R_i2[0]*2.)+R_i2[0]*2.;
					}
				}
			}
			const double R_kick = R;
			/* 2) azimuthal angle */
			phi=0;
			/* 3) vertical positon */
			if (sigma_z > 0)
				z=fabs(ran_exp(seed,sigma_z,0.));
			else
				z=0;
			/*
			 Initial velocities
			*/
			/* velocity dispersion */
			/* 4) radial velocity */
			//csigma_v = sigma_v+dsigma_v/2.-dsigma_v/Rmax_disc*R;
			
			if(R>R_Sun) {csigma_v = sigma_v;}
			else{csigma_v = sigma_v + dsigma_v*(1.-R/R_Sun);}
			
            if(use_random=='y') {v_R=ran_gau(seed,csigma_v);}
            else {v_R = 0;}
            
			/* 5) azimuthal velocity */
			v0_phi=sqrt(R_ave*getValue(F_R1,R_i1,z_k1,(R_ave),0.,1))/CONV1;
            if(use_random=='y') {v_phi=v0_phi+ran_gau(seed,csigma_v);}
            else {v_phi=v0_phi;}
			/* 6) vertical velocity */
            if(use_random=='y') {v_z=ran_gau(seed,csigma_v);}
            else {v_z=0;}
			
			/* 
			 Kick velocities
			*/
			/* direction distributions */
			phi_k=ran0(&seed)*2*PI;
            //cout<<"cp_vk"<<cp_vk<<"phi_k"<<phi_k/2/PI*180.<<endl;
            if(use_random=='y') {theta_k = kick_angle*ran0(&seed)*PI/180.;}
            else {theta_k = kick_angle*PI/180.;}
            
			/* modulus distribution */
			// *************** MOVE TO FUNCTIONS.CPP??? *****************
			v_k=v_k_max;
			int v_k_iter=1;
			while (v_k >= v_k_max) 
			{
				if (strcmp(v_k_distr, "gau") == 0) {v_k=(v_k_0+ran_gau(seed,h_v_k))*pow(cos(theta_k),alpha_cos);}
                
                else if (strcmp(v_k_distr, "exp") == 0){v_k=ran_exp(seed,h_v_k,v_k_thres)*pow(cos(theta_k),alpha_cos);}
				else if (strcmp(v_k_distr, "pow") == 0){v_k=ran_pow(seed,alpha_v_k,csigma_v,v_k_thres)*	pow(cos(theta_k),alpha_cos);}
				else if (strcmp(v_k_distr, "lin") == 0){v_k=ran_lin(seed,ab_v_k)*pow(cos(theta_k),alpha_cos);}
				else if (strcmp(v_k_distr, "con") == 0){v_k=h_v_k;}
				if (v_k > v_k_max)
					//cout << "Kick velocity above limit, R = " << R <<" v_k = "<<v_k <<endl;
				v_k_iter++;
			}
			//cout<<"Kick velocity = "<<v_k<<endl;
			//v_k = h_v_k - 100*(R-1.5);
			//if(v_k<70) v_k=70;
			
			/* 
			 Resulting initial velocity components 
			 */
			
			// cout<<"v_k="<<v_k<<", theta="<<theta_k*(180/PI)<<endl;
			
			v_R=v_R+v_k*sin(theta_k)*sin(phi_k+PI/2.);
			v_phi=v_phi+v_k*sin(theta_k)*cos(phi_k+PI/2.);
			v_z=v_z+v_k*cos(theta_k);
								
			double v_z_0 = v_z;
			if (iter_mode =='n') {//replace with n
				//fprintf(f_veldistr,"%f %f %f %f %f %f \n", theta_k*RAD, phi_k*RAD, v_k, v_R, v_phi, v_z);
				if (v_k > v_k_thres) {
					//fprintf(f_veldistr2,"%f %f %f %f %f %f \n", theta_k*RAD, phi_k*RAD, v_k, v_R, v_phi, v_z);
				}
			}
			
			/* 
			 Switch to new coordinates: positions [kpc], velocities [kpc/Myr]
			 */
			t=0., x1=R, x2=phi, x3=z, x4=v_R*CONV1, x5=v_phi*CONV1, x6=v_z*CONV1;
			
			/* 
			 Initial Energy and Angular Momentum of the particle
			 */
			if (nFiles == 2) {
				if ((x1 < R_i2[nR2-subdim1]) && (fabs(x3) <z_k2[nz2-subdim2])) {
					Etot0=.5*(x4*x4+x5*x5+x6*x6)
					+getValue(pot_ik2,R_i2,z_k2,x1,fabs(x3),1); //(kpc/myr)^2
					Etot_M0=Etot0*(cn_dot+cn_accr)/cn_dot; 
					// but here cn_accr is zero
				}
				else {
					Etot0=.5*(x4*x4+x5*x5+x6*x6)
					+getValue(pot_ik1,R_i1,z_k1,x1,fabs(x3),1); //(kpc/myr)^2
					Etot_M0=Etot0*(cn_dot+cn_accr)/cn_dot;  
				}
			} 
			else {
				Etot0=.5*(x4*x4+x5*x5+x6*x6)
				+getValue(pot_ik1,R_i1,z_k1,x1,fabs(x3),1); //(kpc/myr)^2
				Etot_M0=Etot0*(cn_dot+cn_accr)/cn_dot;  
			}
			Ltot0=sqrt(x3*x3*x5*x5+(x3*x4-x1*x6)*(x3*x4-x1*x6)+x1*x1*x5*x5);
			L_z0=(x1*x5);
			L_z=L_z0;
			Ltot_M0=Ltot0*(cn_dot+cn_accr)/cn_dot; 
			L_z_M0=L_z0*(cn_dot+cn_accr)/cn_dot; 
			if (iter_mode == 'n'){//replace with n
				if (R <= R_plot+deltaR_plot 
					&& R >= R_plot-deltaR_plot 
					&& v_k <= v_plot+deltav_plot 
					&& v_k >= v_plot-deltav_plot)
                                        double nonewew=133.3;
					/*fprintf(f_orbit,"%.1f %.4f %.4f %.4f %.3f %.3f %.3f %.9f %.5f %.9f %.5f %.5f %.9f %.9f  %.9f 1 1 1 1 1 1 \n",
							t,x1,x2,x3,x4/CONV1,x5/CONV1,x6/CONV1,
							Etot0,Ltot0,L_z0,
							getValue(pot_ik1,R_i1,z_k1,x1,fabs(x3),1), 
							(cn_dot+cn_accr)/cn_dot,
							Etot_M0,Ltot_M0,L_z_M0);*/
			}
//cout<<"cycle in time line 1329"<<endl;		 
 	
			/* 
			 ------------------------------
			 3) CYCLE IN TIME (INTEGRATION)
			 ------------------------------
			 */
			/* settings */
			below=false;
			keep_integr=true;
			int n_dt=0;
			/*
			 Selection of v_k above threshold
			 */
			if (v_k > v_k_thres) {
				np_halo+=1;
				cp_halo++;
				mean_v_k+=(v_k*v_k);
				//while (t>= showt_min && t< showt_max && t <= t_max && keep_integr == true){
				while (t <= t_max and keep_integr == true) {
					/* Check if the particle is inside the potential grid */
					if ((x1 < R_i1[nR1-subdim1-1]) && (fabs(x3) <z_k1[nz1-subdim2-1])){
						
						/* 
						 Transform to cartesian coordinates 
						 (ONLY FOR INTEGRATION)
						 */
						cyltocar(x1,x2,x3,x4,x5,x6);
						
						/* definitions for refining integration step */
						x_prev=sqrt(x1*x1+x2*x2+x3*x3);
						v_prev=sqrt(x4*x4+x5*x5+x6*x6);
						x1_prev=x1;
						x2_prev=x2;
						x3_prev=x3;
						x4_prev=x4;
						x5_prev=x5;
						x6_prev=x6;
						cn_accr_prev=cn_accr;
						
						/* 
						 ========================
						 ONE STEP R-K INTEGRATION
						 ========================
						 */
						
						if ((nFiles == 2)
							&&((x1*x1+x2*x2) < R_i2[nR2-subdim1-1]) 
							&& (fabs(x3) <z_k2[nz2-subdim2-1])) {
							if (strcoll(effects,"accretion") ==0 || strcoll(effects,"accrdrag")==0 || strcoll(effects,"accrnodrag")==0) {
								if (strcoll(accr_dens,"build") ==0)
									cn_accr+=dyn2.accretion(x1,x2,fabs(x3),x6,
															cold_rho_ik,
															cold_v_R_ik,
															cold_v_phi_ik,
															cold_v_z_ik,
															R_i1,
															z_k1,
															cn_dot+cn_accr
															)*dt;
								else {
									cn_accr+=dyn2.accretion(x1, x2, fabs(x3),x6,
															accr_dens,
															accr_vel,
															cn_dot,
															t)*dt;
								}
							}
							dyn2.rk4(t,dt,x1, x2, x3, x4, x5, x6,cn_dot+cn_accr);
							//		  dyn2.rk4(dt,x1, x2, x3, x4, x5, x6,cn_dot);
						}
						else { // nFiles == 1 or in File 1
							if (strcoll(effects,"accretion") ==0 || strcoll(effects,"accrdrag")==0 || strcoll(effects,"accrnodrag")==0) {
								if (strcoll(accr_dens,"build") ==0)
									cn_accr+=dyn1.accretion(x1,x2,fabs(x3),x6,
															cold_rho_ik,
															cold_v_R_ik,
															cold_v_phi_ik,
															cold_v_z_ik,
															R_i1,
															z_k1,
															cn_dot+cn_accr
															)*dt;
								else {
									cn_accr+=dyn1.accretion(x1, x2, fabs(x3),x6,
															accr_dens,
															accr_vel,
															cn_dot,
															t)*dt;
								}
							}
							// cout<<"z:"<<fabs(x3)<<endl;
							dyn1.rk4(t,dt,x1, x2, x3, x4, x5, x6,cn_dot+cn_accr);
							//		  dyn1.rk4(dt,x1, x2, x3, x4, x5, x6,cn_dot);
						}
						
						/*
						 CHECK INTEGRATION STEP
						 */
						delta_x=fabs(sqrt(x1*x1+x2*x2+x3*x3)-x_prev);
						delta_v=(fabs(sqrt(x4*x4+x5*x5+x6*x6)-v_prev))/CONV1;
						
						n_refin=0;
						while ((delta_x > delta_x_max 
								|| delta_v > delta_v_max)
							   || 
							   (delta_x < delta_x_min
								|| delta_v < delta_v_min)
							   && n_refin < n_refin_max) {
							if (integr_decr == false) {
								//		    cout << ".";
							}
							if (delta_x > delta_x_max || delta_v > delta_v_max)
								dt=dt/3.;
							else 
								dt=dt*3.;
							
							x1=x1_prev;
							x2=x2_prev;
							x3=x3_prev;
							x4=x4_prev;
							x5=x5_prev;
							x6=x6_prev;
							cn_accr=cn_accr_prev;
							
							if ((nFiles == 2) // refined
								&& ((x1*x1+x2*x2) < R_i2[nR2-subdim1-1]) 
								&& (fabs(x3) <z_k2[nz2-subdim2-1])) {
								if (strcoll(effects,"accretion") ==0 || strcoll(effects,"accrdrag")==0 || strcoll(effects,"accrnodrag")==0) {
									if (strcoll(accr_dens,"build") ==0)
										cn_accr+=dyn2.accretion(x1,x2,fabs(x3),x6,
																cold_rho_ik,
																cold_v_R_ik,
																cold_v_phi_ik,
																cold_v_z_ik,
																R_i1,
																z_k1,
																cn_dot+cn_accr
																)*dt;
									else {
										cn_accr+=dyn2.accretion(x1, x2, fabs(x3),x6,
																// here there was
																// accr_norm, 
																accr_dens,
																accr_vel,
																cn_dot,
																t)*dt;
									}
								}
								//cout << "refin " << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5 << " " << x6 << " " << cn_dot+cn_accr<< endl;
								dyn2.rk4(t,dt,x1, x2, x3, x4, x5, x6,cn_dot+cn_accr);
								//		    dyn2.rk4(dt,x1, x2, x3, x4, x5, x6,cn_dot);
							}
							else { // nFiles == 1, refined
								if (strcoll(effects,"accretion") ==0 || strcoll(effects,"accrdrag")==0 || strcoll(effects,"accrnodrag")==0) {
									if (strcoll(accr_dens,"build") ==0)
										cn_accr+=dyn1.accretion(x1,x2,fabs(x3),x6,
																cold_rho_ik,
																cold_v_R_ik,
																cold_v_phi_ik,
																cold_v_z_ik,
																R_i1,
																z_k1,
																cn_dot+cn_accr
																)*dt;
									else {
										cn_accr+=dyn1.accretion(x1, x2, fabs(x3),x6,
																// here there was
																// accr_norm,
																accr_dens,
																accr_vel,
																cn_dot,
																t)*dt;
									}
								}
								dyn1.rk4(t,dt,x1, x2, x3, x4, x5, x6,cn_dot+cn_accr);
								//		    dyn1.rk4(dt,x1, x2, x3, x4, x5, x6,cn_dot);
							}
							delta_x=fabs(sqrt(x1*x1+x2*x2+x3*x3)-x_prev);
							delta_v=fabs(sqrt(x4*x4+x5*x5+x6*x6)-v_prev)/CONV1;
							integr_decr=true;
							n_refin+=1;
						}
						/*
						 END INTEGRATION
						 */
						//		fclose("orbit.dat");
						
						/* Back to cylindric coordinates */
						cartocyl(x1,x2,x3,x4,x5,x6);
						
						//azimuthal gap between spiral arms and cloud after a time t from the launch
						double phigap = 0;
						if(velocity_gap=='y' && create_arms=='y') 
						{
							phigap = 0.06*t*pattern_speed-x2;
							//cout<<"time="<<t<<", vphi="<<x5/CONV1<<", x2="<<x2<<", phiarm="<<0.06*t*pattern_speed<<", gap="<<phigap<<endl;
						}
						//		cout << t << " " << x5/CONV1 << " " << x1
						//		  << " " << x3 << " " << x4/CONV1 << " " << x6/CONV1 << endl;
						
						t+=dt;
						
						/* Check if the new position is inside the potential grid */
						//cout<<"R="<<x1<<" z="<<fabs(x3)<<endl;
						if ((x1 < R_i1[nR1-subdim1-1]) && (fabs(x3) <z_k1[nz1-subdim2-1])){
							/* Calculate new Energy and Angular Momenta */
							if (nFiles == 2) {
								if ((x1 < R_i2[nR2-subdim1]) 
									&& (fabs(x3) <z_k2[nz2-subdim2])) {
									Etot=(.5*(x4*x4+x5*x5+x6*x6)
										  +getValue(pot_ik2,R_i2,z_k2,x1,fabs(x3),1));
									Etot_M=Etot*(cn_dot+cn_accr)/cn_dot; 
								}
								else {
									Etot=(.5*(x4*x4+x5*x5+x6*x6)
										  +getValue(pot_ik1,R_i1,z_k1,x1,fabs(x3),1));
									Etot_M=Etot*(cn_dot+cn_accr)/cn_dot;  
								}
							} 
							else {
								Etot=(.5*(x4*x4+x5*x5+x6*x6)
									  +getValue(pot_ik1,R_i1,z_k1,x1,fabs(x3),1));
								Etot_M=Etot0*(cn_dot+cn_accr)/cn_dot;  
							}
							Ltot=sqrt(x3*x3*x5*x5+(x3*x4-x1*x6)
									  *(x3*x4-x1*x6)+x1*x1*x5*x5);
							Ltot_M=Ltot*(cn_dot+cn_accr)/cn_dot;  
							dL_z=((x1*x5)-L_z)
							*(cn_dot+cn_accr)/cn_dot;  //kpc^2/Myr
							L_z=(x1*x5);
							L_z_M=L_z*(cn_dot+cn_accr)/cn_dot;  
							
							//		  cout << x3 << " " << cn_accr << " " << L_z << endl;
							
							/* Calculate transfer of angular momentum to the halo */
							if (strcoll(effects,"drag") ==0) {
								
								for (int hi=0;hi<nR1-1;hi++) 
									if (R_i1[hi] <= x1 && x1 < R_i1[hi+1]) 
										for (int hk=0;hk<nz1-1;hk++) 
											if (z_k1[hk] <= fabs(x3) && fabs(x3) < z_k1[hk+1]) {
												cdLz_hot_ik1[hi][hk]-=dL_z*dt/dt_default;
												cdvphi_hot_ik1[hi][hk]-=dL_z
												/(2.*PI*stepR1*stepz1*(x1-delta_x/2.)
												  *(x1-delta_x/2.)*rho_hot_ik1[hi][hk])
												*dt/dt_default; // kpc/Myr Mo^-1
											}
								//		    HotHalo.LTransfer(cdLz_hot_ik1,
								//		  halo.transfer
							}
						} 
						else {
							//cout << "Warning: particles outside the potential grid! \n";
							n_lost += PI*(R2*R2 - R1*R1);
							keep_integr=false;
						}
						/* 
						 =========================
						 CREATE ARTIFICIAL CUBES
						 =========================
						 */ 
						/* 
						 -------------------------------
						 4) CYCLE IN AZIMUTH (PHI ANGLE)
						 -------------------------------
						 */
						
						if (single_event=='n' || (single_event=='y' && t>showt_min && t<=showt_max) || (single_event=='y' && second_burst=='y' && t>showt2_min && t<=showt2_max)){
							
							avg_R += fabs(x1)*(cn_dot+cn_accr);
							avg_z += fabs(x3)*(cn_dot+cn_accr);
							avg_weight += (cn_dot+cn_accr);
							for (jj=0; jj<np_phi; jj++){
								int sign=1;
								/* 
								 -----------------------------------------
								 5) CYCLE IN Z (MIRROWING BELOW THE PLANE)
								 -----------------------------------------
								 */
								// **** IN PROGRESS ****
								
								
								for (int kk2=0; kk2<2; kk2++){
									
									/* 
									 PROJECTIONS
									 */

                                    if(use_random=='y') {xx2=ran0(&seed)*360.;}
                                    else{xx2=360*float(jj)/float(np_phi);}
                                    
									if (kk2 == 1) {sign=-1;}
									// x3 and x6 ARE MULTIPLIED BY sign
									
									
									//
									//   External projection 
									
									if (strcoll(projectionType,"external") ==0) {
										
										double phi_halo = xx2;
										if((phi_halo>PhiminSF && phi_halo<PhimaxSF) || (PhiminSF>PhimaxSF && (phi_halo<PhimaxSF || phi_halo>PhiminSF)))
										{
											projections proj(x1,
															 phi_halo+phigap,
															 x3*sign,
															 x4/CONV1,
															 x5/CONV1,
															 x6/CONV1*sign);
											proj.external_view(incl, PA, xi, yi, vrad);
                                            float xit;
                                            float yit;
                                            float cosPA;
                                            float sinPA;
                                            cosPA=cos(PA/RAD-PI*0.5);
                                            sinPA=sin(PA/RAD-PI*0.5);
                                            xit=xi;
                                            yit=yi;
                                            xi=-xi*cosPA+yit*sinPA;
                                            yi=-xit*sinPA-yit*cosPA;
                                            
											arms_weight = 0;
											//SPIRAL ARMS (in the halo)
											if(create_arms=='y' && R_kick>=Rstart_arms)
											{
												if(phi_halo>360) phi_halo-=360;
												for(int k=0;k<n_R_arms;k++)
												{
													R_overdensity1[k] = R_logspiral(phi_halo+k*360, a1_pitch, a1_amp);
													R_overdensity2[k] = R_logspiral(phi_halo+k*360, a2_pitch, a2_amp);
													R_overdensity3[k] = R_logspiral(phi_halo+k*360, a3_pitch, a3_amp);
													R_overdensity4[k] = R_logspiral(phi_halo+k*360, a4_pitch, a4_amp);
												}
												for(int k=0;k<n_R_arms;k++) 
												{
													arms_weight += (a1_overdens-1.)*exp(-(R_kick-R_overdensity1[k])*(R_kick-R_overdensity1[k])/(2.*a1_thickness*a1_thickness))
													+  (a2_overdens-1.)*exp(-(R_kick-R_overdensity2[k])*(R_kick-R_overdensity2[k])/(2.*a2_thickness*a2_thickness))
													+  (a3_overdens-1.)*exp(-(R_kick-R_overdensity3[k])*(R_kick-R_overdensity3[k])/(2.*a3_thickness*a3_thickness))
													+  (a4_overdens-1.)*exp(-(R_kick-R_overdensity4[k])*(R_kick-R_overdensity4[k])/(2.*a4_thickness*a4_thickness));
												}
											}
											
											// Building the final cube
											xc=int(RAsize/2+xi/fabs(cdelt1)*kpctograd+x_shift+0.5)-1; 
											// in this way the centre is 127
											yc=int(DECsize/2+yi/cdelt2*kpctograd+y_shift+0.5)-1;
											vc=int((vrad+v_shift)/cdelt3+0.5); 
											// NB this has been tested with a face-on galaxy
											if (xc < RAsize && xc >= 0 
												&& yc < DECsize && yc >=0 
												&& vc < VELsize && vc >= 0) {
												if (mirror_ra =='y') {xc=RAsize-xc-1;}
												if (mirror_dec =='y') {yc=DECsize-yc-1;}
												if ((x6>0 && fabs(x3)>z_above_in) || (x6<0 && fabs(x3)>z_above_fin)) 
												{
													counts_halo_above+=(cn_dot+cn_accr)
													*deltaR*deltaR*dt/dt_default;
													counts_accreted_above += cn_accr*deltaR*deltaR*dt/dt_default;
													
													above=1;
												}
												if (strcmp(visible, "all") == 0 
													|| 
													((strcmp(visible, "down") == 0 && x6 < 0) ||
													 (strcmp(visible, "down") == 0 && x3 < 0))
                                                    ||
                                                    ((strcmp(visible, "up") == 0 && x6 > 0) ||
                                                     (strcmp(visible, "up") == 0 && x3 < 0))
													|| 
													((strcmp(visible, "R_max")==0 && x4 < 0 && x6 < 0)
													 || (strcmp(visible, "R_max") == 0 && x3 < 0)) 
													|| 
													((strcmp(visible, "partial")==0 && x6/CONV1 < v_z_0*(1.-ion_frac))
													 || (strcmp(visible, "partial") == 0 && x3 < 0)))
												{
													
													counts_HI+=(1.+arms_weight)*(cn_dot+cn_accr)
													*deltaR*deltaR*dt/dt_default;
            
													outcube[xc][yc][vc]+=(1.+arms_weight)*(cn_dot+cn_accr)
													*deltaR*deltaR*dt/dt_default;
													outcube_tot[xc][yc]+=(1.+arms_weight)*(cn_dot+cn_accr)
													*deltaR*deltaR*dt/dt_default;
													//			  outcube[xc][DECsize-yc][vc]+=(cn_dot+cn_accr)
													//			    *deltaR*deltaR*dt/dt_default;
													//			  outcube_tot[xc][DECsize-yc]+=(cn_dot+cn_accr)
													//			    *deltaR*deltaR*dt/dt_default;
													Weight_noarms += cn_dot+cn_accr;
													Weight_witharms += (1.+arms_weight)*(cn_dot+cn_accr);
													if ((x6>0 && fabs(x3)>z_above_in) || (x6<0 && fabs(x3)>z_above_fin)) 
													{
														counts_HI_above+=(1.+arms_weight)*(cn_dot+cn_accr)
														*deltaR*deltaR*dt/dt_default;
														above=1;
													}
												}
												//
												//  NB. cn_dot is an outflow density in counts/kpc^2
												//  that is added np_phi times for each R
												//  the outflow "mass" contributed at radius R will be
												//  delta_mass_out[counts]= 2pi *R *cn_dot *deltaR
												//  = np_phi *deltaR *cn_dot *deltaR
												//  This "mass" has to be normalised below using the
												//  input paramter total mass of the gaseous halo
												//  (*ntomass).
												//
											} 
											else if (out_boundaries == false) {
												//cout << "Warning: particles outside the cube boundaries (xc=" << xc << "," << "yc=" << yc << "," << "vc=" << vc << ")" << endl;
												n_lost += PI*(R2*R2 - R1*R1);
												out_boundaries=true;
											}
										}
									}
									
									/* 
									 Internal projection 
									 */ 
									if (strcoll(projectionType,"internal") ==0) {
										/*
										 //previous values
										 projections proj_prev(x1_prev,
										 x2_prev+xx2,
										 x3_prev*sign,
										 x4_prev/CONV1,
										 x5_prev/CONV1,
										 
										 x6_prev/CONV1*sign);
										 float l_prev, b_prev, d_prev, vlos_prev;
										 proj_prev.internal_view(R_Sun,v0_phi_Rsun,l_prev,b_prev,d_prev,vlos_prev);
										 int xc_prev = int(RAsize/2+l_prev/fabs(cdelt1)+0.5);
										 int yc_prev = int(DECsize/2+b_prev/cdelt2+0.5);
										 int vc_prev = int(VELsize/2-vlos_prev/cdelt3+0.5); 
										 if(xc_prev == RAsize) xc_prev = 0;
										 */
										//last values
										double phi_halo = xx2;
										if((phi_halo>PhiminSF && phi_halo<PhimaxSF) || (PhiminSF>PhimaxSF && (phi_halo<PhimaxSF || phi_halo>PhiminSF)))
										{
											projections proj(x1,
															 phi_halo+phigap,
															 x3*sign,
															 x4/CONV1,
															 x5/CONV1,
															 x6/CONV1*sign);
											proj.internal_view(R_Sun,v0_phi_Rsun,l,b,d,vlos); // output: l(degree),b(degree),vlos;
											arms_weight = 0;
											//SPIRAL ARMS (in the halo)
											if(create_arms=='y' && R_kick>=Rstart_arms)
											{
												if(phi_halo>360) phi_halo-=360;
												for(int k=0;k<n_R_arms;k++)
												{
													R_overdensity1[k] = R_logspiral(phi_halo+k*360, a1_pitch, a1_amp);
													R_overdensity2[k] = R_logspiral(phi_halo+k*360, a2_pitch, a2_amp);
													R_overdensity3[k] = R_logspiral(phi_halo+k*360, a3_pitch, a3_amp);
													R_overdensity4[k] = R_logspiral(phi_halo+k*360, a4_pitch, a4_amp);
												}
												for(int k=0;k<n_R_arms;k++) 
												{
													arms_weight += (a1_overdens-1.)*exp(-(R_kick-R_overdensity1[k])*(R_kick-R_overdensity1[k])/(2.*a1_thickness*a1_thickness))
													+  (a2_overdens-1.)*exp(-(R_kick-R_overdensity2[k])*(R_kick-R_overdensity2[k])/(2.*a2_thickness*a2_thickness))
													+  (a3_overdens-1.)*exp(-(R_kick-R_overdensity3[k])*(R_kick-R_overdensity3[k])/(2.*a3_thickness*a3_thickness))
													+  (a4_overdens-1.)*exp(-(R_kick-R_overdensity4[k])*(R_kick-R_overdensity4[k])/(2.*a4_thickness*a4_thickness));
												}
											}
											
											//DELETE THIS!
											//if(phi_halo<-90 || phi_halo>-30) arms_weight = 0;
											
											xc=int(RAsize/2+l/fabs(cdelt1)+0.5);
											yc=int(DECsize/2+b/cdelt2+0.5);
											vc=int(VELsize/2-vlos/cdelt3+0.5);
											if (xc == RAsize) xc=0;
											//cout<<l<<" "<<xc<<endl;
											//if(fabs(b)>61)cout<<b<<" "<<yc<<endl;
											//cout<<vlos<<" "<<vc<<endl<<endl;
											
											/* Building the final cube */
											if (vc < VELsize && vc >= 0) {
												if ((x6>0 && fabs(x3)>z_above_in) || (x6<0 && fabs(x3)>z_above_fin)) 
												{
													above=1;
													counts_halo_above+=(1.+arms_weight)*(cn_dot+cn_accr)
													*deltaR*deltaR*dt/dt_default;
													counts_accreted_above += (1.+arms_weight)*cn_accr*deltaR*deltaR*dt/dt_default;
													
													//VELOCITY DISTRIBUTION OF THE ACCRETED MATERIAL
													//RADIAL
													for(int i=0; i<nvR;i++)
													{
														if(x4/CONV1>=(vR[i]-delta_vR/2.) && x4/CONV1<(vR[i]+delta_vR/2.))
														{
															m_vR[i]+=(1.+arms_weight)*cn_accr*deltaR*deltaR*dt/dt_default;
															break;
														}
													}
													//AZIMUTHAL
													for(int i=0; i<nvphi; i++)
													{
														if(x5/CONV1>=(vphi[i]-delta_vphi/2.) && x5/CONV1<(vphi[i]+delta_vphi/2.))
														{
															m_vphi[i]+=(1.+arms_weight)*cn_accr*deltaR*deltaR*dt/dt_default;
															break;
														}
													}
													
													//SPECIFIC ANGULAR MOMENTUM
													
													// if(cp_vk!=cp_vk_prev)
													
													if(x1>7. && x1<10.)
													{
														counts_halo_local+=(1.+arms_weight)*(cn_dot+cn_accr)
														*deltaR*deltaR*dt/dt_default;
														counts_accreted_local += (1.+arms_weight)*cn_accr*deltaR*deltaR*dt/dt_default;  
													}
													//above=0;
												}
												//if(x6/CONV1>v_z_0) cout<<"vz="<<x6/CONV1<<" while vk="<<v_k<<endl;
												if (strcmp(visible, "all") == 0 
													|| 
                                                    ((strcmp(visible, "down") == 0 && x6 < 0) ||
                                                     (strcmp(visible, "down") == 0 && x3 < 0))
                                                    ||
                                                    ((strcmp(visible, "up") == 0 && x6 > 0) ||
                                                     (strcmp(visible, "up") == 0 && x3 < 0))
													|| 
													((strcmp(visible, "R_max")==0 && x4 < 0 && x6 < 0)
													 || (strcmp(visible, "R_max") == 0 && x3 < 0)) 
													|| 
													((strcmp(visible, "partial")==0 && x6/CONV1 <= v_z_0*(1.-ion_frac))
													 || (strcmp(visible, "partial") == 0 && x3 < 0))) {
														
														counts_HI+=(1.+arms_weight)*(cn_dot+cn_accr)*deltaR*deltaR*dt/dt_default;
														Weight_noarms += cn_dot+cn_accr;
														Weight_witharms += (1.+arms_weight)*(cn_dot+cn_accr);
														if ((x6>0 && fabs(x3)>z_above_in) || (x6<0 && fabs(x3)>z_above_fin)) 
														{
															//cout<<"z="<<fabs(x3)<<", vz="<<x6/CONV1<<", vthresh="<<v_k*(1.-ion_frac)<<", vkick="<<v_k<<endl;
															counts_HI_above+=(1.+arms_weight)*(cn_dot+cn_accr)
															*deltaR*deltaR*dt/dt_default;
															//above=1;
														}
														
														/*
														 The direction of the Sun should be at phi=0
														 left side: 0-180 deg, receding
														 right side: 180-360 deg, approaching
														 */
														if (d>distance_cut){
															
															/*
															 Here cn_dot is an outflow density in counts/kpc^2
															 that is added np_phi times for each R
															 the outflow "mass" contributed at radius R will be
															 delta_mass_out[counts]= 2pi *R *cn_dot *deltaR
															 = np_phi *deltaR *cn_dot *deltaR
															 the flux contribution changes every dt because the
															 distance of the cloud changes.
															 This "mass" has to be normalised below using the
															 input paramter total mass of the gaseous halo
															 (*ntomass).
															 */
															pixtopc=fabs(cdelt1*cdelt2)/RAD/RAD*d*d*1.e6; 
															// pix -> pc^2
															/*
															 outcube[xc][yc][vc]+=(cn_dot+cn_accr)
															 *deltaR*deltaR*dt/dt_default
															 /pixtopc;
															 outcube_tot[xc][yc]+=(cn_dot+cn_accr)
															 *deltaR*deltaR*dt/dt_default
															 /pixtopc;
															 */
															
															d_resolved = RAD/sqrt(abs(cdelt1*cdelt2))*cloudsize/1000.; // distance at which a cloud of size = cloudsize is resolved with
															//d_resolved = -1.;//a little test															// the assumed spatial resolution
															if (d < d_resolved)
															{ // the cloud is resolved
																// emission is spread in all directions with gaussian distributions
																pixsize1=(int)((atan(cloudsize/d/1000.)*RAD/fabs(cdelt1)+.5)/cos(b/RAD));
																pixsize2=(int)(atan(cloudsize/d/1000.)*RAD/fabs(cdelt2)+.5);
																if(use_random=='y')
                                                                {
                                                                    if(pixsize1<1 || pixsize2<1) cout<<"Warning!!!"<<endl;
                                                                    for (int p=0; p<pixsize1*pixsize2; p++)
                                                                    {
                                                                        // randomization in each direction
                                                                        pxc=(int)(xc+ran_gau(seed,pixsize1/3.));
                                                                        pyc=(int)(yc+ran_gau(seed,pixsize2/3.));
                                                                        pvc=(int)(vc+ran_gau(seed,csigma_v/cdelt3));
                                                                        if (pxc >= RAsize) pxc=pxc-RAsize;
                                                                        if (pxc < 0) pxc=pxc+RAsize;
                                                                        if (pyc >= DECsize)
                                                                        {
                                                                            pyc = 2*DECsize - pyc;
                                                                            if(pxc<RAsize/2) {pxc+=RAsize/2;}
                                                                            else{pxc-=RAsize/2;}
                                                                        }
                                                                        if (pyc < 0)
                                                                        {
                                                                            pyc = abs(pyc);
                                                                            if(pxc<RAsize/2) {pxc+=RAsize/2;}
                                                                            else{pxc-=RAsize/2;}
                                                                        }
                                                                        if (pvc >= VELsize) pyc=VELsize;
                                                                        if (pvc < 0) pyc=0;
                                                                        // effect of curvature: number of pixels as a function of b
                                                                        // note that differently from below here the latitude is recalculate at the exact position of pyc
                                                                        
                                                                        nPixels = (int)(1. / cos( (pyc - DECsize/2) * cdelt2 /RAD) / fabs(cdelt1) + ran0(&seed));//it's ok!
                                                                        for (int pp=0; pp< nPixels; pp++)
                                                                        {
                                                                            ppxc = (int)(pxc - nPixels/2 + pp + 0.5);
                                                                            if (ppxc >= RAsize) {ppxc=ppxc-RAsize;}
                                                                            if (ppxc < 0) {ppxc=ppxc+RAsize;}
                                                                            if(ppxc>=0 && ppxc<RAsize)
                                                                            {outcube[ppxc][pyc][pvc]+=(cn_dot+cn_accr)
                                                                            *deltaR*deltaR*dt/dt_default/pixtopc
                                                                            /pixsize1/pixsize2
                                                                            *(ran_gau(seed,.3)+1) // WHY is this there? (randomization around 1)
                                                                            /nPixels;
                                                                            outcube_tot[ppxc][pyc]+=(cn_dot+cn_accr)
                                                                            *deltaR*deltaR*dt/dt_default/pixtopc
                                                                            *fabs(cdelt3)
                                                                            /pixsize1/pixsize2*(ran_gau(seed,.3)+1)/nPixels;}
                                                                            
                                                                        }
                                                                    }
                                                                }
                                                                else
                                                                {
                                                                if(pixsize1>RAsize/2.) pixsize1=RAsize/2.;
																if(pixsize1<1 || pixsize2<1) {cout<<"Warning!!!"<<endl;}
																int i1 = xc-pixsize1; int i2=xc+pixsize1;
																int j1 = yc-pixsize2; int j2 = yc+pixsize2;
																int k1 = vc-3*csigma_v/cdelt3; int k2 = vc+3*csigma_v/cdelt3;
																if (i1 >= RAsize) i1=i1-RAsize; if (i2 >= RAsize) i2=i2-RAsize;
																if (i1 < 0) i1=i1+RAsize; if (i2 < 0) i2=i2+RAsize; 
																if (j1 >= DECsize)
																{	
																	j1 = 2*DECsize - j1;
																	if(i1<RAsize/2) {i1+=RAsize/2;}
																	else{i1-=RAsize/2;}
																	if(i2<RAsize/2) {i2+=RAsize/2;}
																	else{i2-=RAsize/2;}
																}
																if (j2 >= DECsize)
																{	
																	j2 = 2*DECsize - j2;
																	if(i1<RAsize/2) {i1+=RAsize/2;}
																	else{i1-=RAsize/2;}
																	if(i2<RAsize/2) {i2+=RAsize/2;}
																	else{i2-=RAsize/2;}
																}
																if (j1 < 0)
																{
																	j1 = abs(j1);
																	if(i1<RAsize/2) {i1+=RAsize/2;}
																	else{i1-=RAsize/2;}
																	if(i2<RAsize/2) {i2+=RAsize/2;}
																	else{i2-=RAsize/2;}
																}
																if (j2 < 0)
																{
																	j2 = abs(j2);
																	if(i1<RAsize/2) {i1+=RAsize/2;}
																	else{i1-=RAsize/2;}
																	if(i2<RAsize/2) {i2+=RAsize/2;}
																	else{i2-=RAsize/2;}
																}
																if (k1 < 0) k1=0; if (k2 < 0) k2=0; 
																for(int i=i1;i<=i2;i++) for(int j=j1;j<=j2;j++) for(int k=k1;k<=k2;k++)	
																{
                                                                    float value = (1.+arms_weight)*(cn_dot+cn_accr)
																	*deltaR*deltaR*dt/dt_default/pixtopc
																	*myutil::gaussian(xc-i,pixsize1/3.)*myutil::gaussian(yc-j,pixsize2/3.)*myutil::gaussian(vc-k,csigma_v/cdelt3)
																	*fabs(cdelt1*cdelt2*cdelt3);
                                                                    outcube[i][j][k]+=value;
											
                                                                    }
                                                                }
															}
															else {
																// effect of curvature: number of pixels as a function of b
if(use_random=='y') {nPixels = (int)(1./cos(b/RAD)/fabs(cdelt1)+ran0(&seed));}
                                                                else {nPixels = max(1,(int)(1./cos(b/RAD)/fabs(cdelt1)+0.5));}
																//cout<<"nPixels="<<nPixels<<endl;
																for (int p=0; p< nPixels; p++)
																{
																	pxc = (int)(xc - nPixels/2 + p + 0.5);
																	if (pxc >= RAsize) pxc=pxc-RAsize;
																	if (pxc < 0) pxc=pxc+RAsize; 
																	if(pxc>=0 && pxc<RAsize)//is this ok?
																	{
																		float value = (1.+arms_weight)*(cn_dot+cn_accr)
																		*deltaR*deltaR*dt/dt_default/pixtopc
																		/nPixels;
																		outcube[pxc][yc][vc]+=value;
																		
																		outcube_tot[pxc][yc]+=(1.+arms_weight)*(cn_dot+cn_accr)
																		*deltaR*deltaR*dt/dt_default/pixtopc
																		/nPixels;
																		
																	}
																	else {//cout<<"Trouble! pxc="<<pxc<<" for v_k="<<v_k<<", R="<<int(R)<<endl;
																		}
																}
																
															}
															
															
															
														}
													}
											}
											else if (out_boundaries == false) {
												//cout << "Warning: particles outside the cube boundaries (xc="<<xc<<" yc="<<yc<<" vc=" << vc << ")" << endl;
												n_lost += PI*(R2*R2 - R1*R1);
												out_boundaries=true;
											}
										}
									}
									
									/* 
									 Internal projection (HVCS)
									 */
									if (strcoll(projectionType,"hvcs") ==0) {
										projections proj(x1,
														 x2+xx2,
														 x3*sign,
														 x4/CONV1,
														 x5/CONV1,
														 x6/CONV1*sign);
										proj.internal_view_hvcs(R_Sun,v0_phi_Rsun,Rmax_hvcs,zmax_hvcs,l,b,d,vlos); // output: l(degree),b(degree),vlos;
										xc=int(RAsize/2+l/fabs(cdelt1)+0.5);
										yc=int(DECsize/2+b/cdelt2+0.5);
										vc=int(VELsize/2-vlos/cdelt3+0.5);
										if (xc == 360) xc=0;
										if (vc < VELsize && vc >= 0) {
											/* Building the final cube */
											if (d>distance_cut){
												if (xc < RAsize/2) {
													xc+=RAsize/2;
												} 
												else {
													xc-=RAsize/2;
												}
												vlos_max=hvc_upcut[xc][yc];
												vlos_min=hvc_lowcut[xc][yc];
												if (vlos > vlos_max+vlos_cut || vlos < vlos_min-vlos_cut) {
													if (strcmp(v_hvcs_plot, "vdev") == 0){
														if (vlos > vlos_max) {vlos=vlos-vlos_max;}
														if (vlos < vlos_min) {vlos=vlos-vlos_min;}
													}
													xc=int(RAsize/2+l/fabs(cdelt1)+0.5);
													yc=int(DECsize/2+b/cdelt2+0.5);
													vc=int(VELsize/2-vlos/cdelt3+0.5);
													if (xc == 360) xc=0;
													if (vc < VELsize && vc >= 0) {
														if (xc < RAsize/2) {
															xc+=RAsize/2;
														} 
														else {
															xc-=RAsize/2;
														}
														pixtopc=fabs(cdelt1*cdelt2)/RAD/RAD*
														d*d*1.e6; // pix -> pc^2
														outcube[xc][yc][vc]+=(1.+arms_weight)*(cn_dot+cn_accr)
														*deltaR*deltaR*dt/dt_default;
														outcube_tot[xc][yc]+=(1.+arms_weight)*(cn_dot+cn_accr)
														*deltaR*deltaR*dt/dt_default;
														if ((x6>0 && fabs(x3)>z_above_in) || (x6<0 && fabs(x3)>z_above_fin)) {
															counts_HI_above+=(1.+arms_weight)*(cn_dot+cn_accr)
															*deltaR*deltaR*dt/dt_default;
															counts_accreted_above += (1.+arms_weight)*cn_accr*deltaR*deltaR*dt/dt_default;
															//above=1;
														}
													}
												}
											}
										} 
										else if (out_boundaries == false) {
											//cout << "Warning: particles outside the cube boundaries (xc="<<xc<<" yc="<<yc<<" vc=" << vc << ")" << endl;
											n_lost += PI*(R2*R2 - R1*R1);
											out_boundaries=true;
										}
									}
								}
								/* END CYCLE IN z */
								
								
								// Rotsurfs and density at different heights
								R_c=int(fabs(x1/rotsurf_grid)+0.5);
								z_c=int(fabs(x3/rotsurf_grid)+0.5);
								// rotsurf_nvel=500;
								// we set a velocity grid of 1 km/s
								// v_phi_c = v_phi 
								v_phi_c=int(x5/CONV1+0.5);
								// we set a velocity grid of 1 km/s
								// with v_R=0 at rotsurf_nvel/2
								v_R_c=int(x4/CONV1+0.5)+rotsurf_nvel/2;
								if (R_c < rotsurf_nR && R_c >= 0 
									&& z_c < rotsurf_nz && z_c >=0 
									&& v_phi_c < rotsurf_nvel && v_phi_c >= 0) 
								{
									if (strcmp(visible, "all") == 0 
										|| 
                                        ((strcmp(visible, "down") == 0 && x6 < 0) ||
                                         (strcmp(visible, "down") == 0 && x3 < 0))
                                        ||
                                        ((strcmp(visible, "up") == 0 && x6 > 0) ||
                                         (strcmp(visible, "up") == 0 && x3 < 0))
										|| 
										((strcmp(visible, "R_max")==0 && x4 < 0 && x6 < 0)
										 || (strcmp(visible, "R_max") == 0 && x3 < 0)) 
										|| 
										((strcmp(visible, "partial")==0 && x6/CONV1 < v_k*(1.-ion_frac))
										 || (strcmp(visible, "partial") == 0 && x3 < 0)))
									{
										out_v_phi[R_c][z_c][v_phi_c]
										+=(1.+arms_weight)*(cn_dot+cn_accr)*rotsurf_grid*rotsurf_grid
										/np_phi*dt/dt_default;
										out_v_R[R_c][z_c][v_R_c]
										+=(1.+arms_weight)*(cn_dot+cn_accr)*rotsurf_grid*rotsurf_grid
										/np_phi*dt/dt_default;
										out_dens[R_c][z_c]
										+=(1.+arms_weight)*(cn_dot+cn_accr)*rotsurf_grid*rotsurf_grid
										/np_phi*dt/dt_default;
										//cout<<"density in R="<<R_c<<",z="<<z_c<<" = "<<out_dens[R_c][z_c]<<endl;
									}
								}
							}
							/* END CYCLE in AZIMUTH */
						}
						
						/*
						 Outflow and Inflow radii above z_above
						 */
						if ((x6>0 && fabs(x3)>z_above_in) || (x6<0 && fabs(x3)>z_above_fin)) {
							t_above+=dt;
							if (already_above==0) {
								R_out_above=x1;
								already_above=1;
							}
						} 
						else {
							if (already_above==1) {
								R_in_above=x1;
								already_above=0;
							}	  
						}
						
						/* 
						 Writing output file for orbit plot
						 (t, R, phi, z, 
						 */
						if (iter_mode == 'n'){//replace with n
							if (R <= R_plot+deltaR_plot 
								&& R >= R_plot-deltaR_plot 
								&& v_k <= v_plot+deltav_plot 
								&& v_k >= v_plot-deltav_plot)
                            {
								/*fprintf(f_orbit,"%.1f %.4f %.4f %.4f %.3f %.3f %.3f %.9f %.5f %.9f %.5f %.5f %.9f %.5f %.9f %.9f %.9f %.9f %.9f %.9f %.9f \n",
										t,x1,x2,x3,x4/CONV1,x5/CONV1,x6/CONV1,
										Etot,Ltot,L_z, 
										getValue(pot_ik1,R_i1,z_k1,x1,fabs(x3),1), 
										(cn_dot+cn_accr)/cn_dot,
										Etot_M,Ltot_M,L_z_M,
										Etot/Etot0, Ltot/Ltot0, L_z/L_z0,
										Etot_M/Etot_M0, Ltot_M/Ltot_M0, L_z_M/L_z_M0);
						*/		
								//cyltocar(x1,x2,x3,x4,x5,x6);
								//		    fprintf(f_orbit,"%.5f %.5f %.5f \n", getforce(type_pot,1,x1,x2,x3,0,0,0),getforce(type_pot,2,x1,x2,x3,0,0,0),getforce(type_pot,3,x1,x2,x3,0,0,0));
								//		    cartocyl(x1,x2,x3,x4,x5,x6);
							}
						}
					} 
					else {
						//cout << "Warning2: particles outside the potential grid!" << endl;
						n_lost += PI*(R2*R2 - R1*R1);
						keep_integr=false;
					}
					
					/* Conditions to keep on integrating */
					if (integrate_until_tmax == false) {
						if (below == false && x3 < 0) {
							if (second_passage == false 
								&& x1 <= passage_Rmax)
								keep_integr=false;
							else 
								below=true; // keep_integr is already true
						}
						if (below == true && x3 > 0) {
							keep_integr=false;
							below=false;
						}
					}
					n_dt+=1;	    
					if (n_dt > n_dt_max) {
						//printf("Maximum number of integrations exceeded : %i \n",n_dt-1);
						keep_integr=false;
					}
					dt=dt_default;
					/*  
					 if(check_up && fabs(x3)>=z_above && x6>0)
					 {
					 sLz_in = x1*x5/CONV1;
					 Lz_in = sLz*(cn_dot+cn_accr);
					 check_up = false;
					 }
					 if(check_down && fabs(x3)<=z_above && x6<0)
					 {
					 sLz_fin = x1*x5/CONV1;
					 Lz_fin = sLz*(cn_dot+cn_accr);
					 check_down = false;
					 }
					 */
				}
				
				/* END INTEGRATION CYCLE */
				//	    cout << cn_dot << " " << cn_accr << " " << cn_accr/cn_dot << endl;
				
				
				
				///***** MOVIE
				 if (iter_mode == 'n' && make_movie == 'y'){
				 // Initialize image (ImageMagick) //
				 int IMsize=RAsize*DECsize*3;
				 float pixels[IMsize];
				 //Image image(RAsize, DECsize, "RGB", FloatPixel, pixels);
				 float r, g, b, norma=max(outcube_tot);
				 for (int i=0; i<RAsize; i++){
				 for (int j=0; j<DECsize; j++){
				 r=0;
				 g=outcube_tot[i][j]/norma*sqrt(fabs(j-DECsize/2)+1);
				 b=0;
				 if (g > 1) g=1;
				 //image.pixelColor(i,j, ColorRGB(r,g,b));
				 }
				 }
				 movieIncr++;
				 imageName="test/movie/image";
				 stringstream oss;
				 oss << imageName << movieIncr;
				 string output(oss.str());
				 imageName=output+".gif";
				 int zoomRAsize=RAsize+movieIncr;
				 int zoomDECsize=DECsize+movieIncr;
				 //	    image.zoom(Geometry(zoomRAsize,zoomDECsize));
				 //	    image.rotate(100-movieIncr/10.);
				 //image.write(imageName);
				 //readImages( &imageList, imageName );
				 }
				 //*****/
				/*
				 OUTFLOW AND INFLOW RATES
				 */
				n_tau+=1;
				tau+=t;
				if(above)
				{
					n_tau_above+=cn_dot*(R2*R2-R1*R1);
					tau_above+=t*cn_dot*(R2*R2-R1*R1);
				}
				io.CalculateRates(x1,
								  t,
								  t_above,
								  dt,
								  tau,
								  n_tau,
								  R1,
								  R2,
								  cn_dot,
								  cn_accr,
								  v_k,
								  R_in_above,
								  above,
								  ctotal_energy_above);
				
				// TRANSFER OF ANGULAR MOMENTUM TO THE HALO
				for (int hi=0;hi<nR1-1;hi++) {
					for (int hk=0;hk<nz1-1;hk++) {
						cdLz_hot_ik1[hi][hk]=cdLz_hot_ik1[hi][hk]*2.*PI*R_ave*deltaR
						*cn_dot/dt; /* mass kpc^2/Myr Myr^-1 */
						// there was another factor 2 ????
						/* 
						 np_phi*deltaR=2*pi*R_ave
						 */
						cdvphi_hot_ik1[hi][hk]=cdvphi_hot_ik1[hi][hk]*2.*PI*R_ave*deltaR*cn_dot/dt; /* mass kpc/Myr Mo^-1 Myr^-1*/
						dLz_hot_ik1[hi][hk]+=cdLz_hot_ik1[hi][hk];
						dvphi_hot_ik1[hi][hk]+=cdvphi_hot_ik1[hi][hk];
						
						Lz_hot_tot+=cdLz_hot_ik1[hi][hk]*2.;
						cdLz_hot_ik1[hi][hk]=0;
						cdvphi_hot_ik1[hi][hk]=0;
					}
				}
				
				/* 
				 Plotting orbit
				 */
				if (iter_mode == 'n') {//replace with n
					if (R <= R_plot+deltaR_plot and R >= R_plot-deltaR_plot and v_k <= v_plot+deltav_plot and v_k >= v_plot-deltav_plot) {
						//fclose(f_orbit);
                        //string psfile="set out 'long_term/orbit1/extragas_";
                        //psfile+=std::to_string(ii);
                        //psfile+=".ps'";
                        //char psfile_array[31];
                        //strcpy(psfile_array, psfile.c_str());
                        //char psfile[30]="set out 'movie1/extra_"+to_string(ii)+".ps'";
						//plot1(alpha_accr,psfile_array);
						// these plotting routines are only for Fraternali&Binney2008
						//		plot0();
						//		plot0col();
						//if(verbose)cout << "writing output plot file extragas.ps\n";
						//		cout << '\a';
					}
					//fprintf(f_times,"%i %f %f %.2e %f %.2e \n", n_tau, R, v_k, t*1.e6, R_ave, 2*PI*R_ave/v0_phi/CONV1*1e6);
				}
			}
			cp_vk++;
			integr_decr=false;
			//cout<<x1<<" "<<x4/CONV1<<" "<<x5/CONV1<<" "<<x6/CONV1<<" "<<above<<endl;;
			for(int i=0; i<nmyR; i++)
			{
				//if(x1>=(myR[i]-delta_myR/2.) && x1<(myR[i]+delta_myR/2.) && above==1)
				if(x1>=(myR[i]-delta_myR/2.) && x1<(myR[i]+delta_myR/2.))
				{
					//cout<<"R="<<x1<<" z="<<x3<<" Vz="<<x6/CONV1<<endl;
					npart[i] += 1;
					MVrot[i] += (cn_dot+cn_accr)*x5/CONV1;
					MVrot_ratio[i] += (cn_dot+cn_accr)*x5/sqrt(x1*getValue(F_R1,R_i1,z_k1,x1,0.,1));
					M[i] += (cn_dot+cn_accr);
					break;
				}
			}
		}
                
                
		/* END CYCLE in KICK VELOCITIES */
		ii++;
		//      if (integr_decr == 1) {cout << "\n";}
		n_total += cp_halo*PI*(R2*R2 - R1*R1);	
		if(verbose) cout << " number of sample particles : " << cp_halo << endl;
		//	cout<<" counts_HI : "<<counts_HI<<endl;
		total_energy_above+=(ctotal_energy_above);
		
		
	}
	//      fclose(f_orbit);
	/* END CYCLE IN R */
	
	//mass-weighted height
	if(single_event=='y' && verbose) cout<<"Mass-weighted radius : "<<avg_R/avg_weight<<" kpc"<<endl;
	if(single_event=='y' && verbose) cout<<"Mass-weighted height : "<<avg_z/avg_weight<<" kpc"<<endl;
	
	
	/* 
	 %%%%%%%%%%%%%%%%%%%%%%%%
	 OUTPUTS AND PLOTTING
	 %%%%%%%%%%%%%%%%%%%%%%%%
	 */
	/* 
	 Plotting kick velocity distribution and travel times
	 */
	if (iter_mode == 'n'){//replace with n
		//fclose(f_veldistr);
		//fclose(f_veldistr2);
		//fclose(f_times);
		//plot3();
		//plot4();
	}

	np_fraction=1.*np_halo/nRSF/np_vk;
	mean_v_k=sqrt(mean_v_k/np_halo);
	/*if(verbose)
	{
		cout << "Fraction of particles above velocity threshold : " << np_fraction << endl; 
		cout << "Mean travel time : " << tau/n_tau << " Myrs"<< endl;
		cout << "Mean travel time for clouds crossing z="<<z_above_in<<" kpc : " << tau_above/n_tau_above << " Myrs"<< endl;
		cout << "Mean kick velocity : " << mean_v_k << endl;
		cout << "Fraction of particles lost during the integration: " << 100.*n_lost/n_total<<" percent"<<endl; 
	}*/
	/*
	 NORMALIZATION OF THE TOTAL FLUX
	 */
	//    cout << "original counts halo : " << counts_HI << "\n";
	//  Convertions formulae
	//  N_H = 1.25 * 10^21 * B[mJ/beam] * dv [km/s] / beam["^2]
	//  Mo/pc^2 = 10 * B[mJ/beam] * dv [km/s] / beam["^2]
	//  B[my/beam] = 0.1 / dv[km/s] * beam["^2] * B[Mo/pc^2]
	//  T_B[K] = 685.25 / beam["^2] * B[mJy/beam]
	//  T_B[K] = 68.525 * B[Mo/pc^2] / dv[km/s]  
	if (smooth == 'n') {
		beam_x=fabs(cdelt1)*3600.;
		beam_y=cdelt2*3600.;
	}
	pctoJy=0.1*(1.13*beam_x*beam_y)/fabs(cdelt3)/1000; // Mo/pc^2 -> Jy/Beam
	if (Jy2K == 'y')
		JytoK=685.25/(1.13*beam_x*beam_y)*1000.;
	else
		JytoK=1.;

	if(verbose)cout<<"JytoK factor: "<<JytoK<<endl;

	//double tot = 0;
	//for (int i=0;i<RAsize;i++) for (int j=0;j<DECsize;j++) for(int k=0;k<VELsize;k++) tot+=outcube[i][j][k];
	//cout<<" TOTAL NOW = "<<tot<<endl;
	
	/* External galaxy */
	if (strcoll(projectionType,"external") ==0) {
		// for external galaxy counts_HI can be calculated here because each
		// cloud can be considered at the same distance
		counts_HI=0.;
		for (xc=0;xc<RAsize;xc++)
			for (yc=0;yc<DECsize;yc++)
				for (vc=0;vc<VELsize;vc++)
					counts_HI+=outcube[xc][yc][vc];	
		ntomass=M_HI_halo/counts_HI; // mass -> Mo
        

        
        //rat=M_HI_halo/2.36e5/D**2/dv/tot;
		pixtopc=fabs(cdelt1*cdelt2)/RAD/RAD*dist*dist*1.e12;
		//	cout << ntomass << " " << pixtopc << " " << pctoJy << endl;
		// pix -> pc^2
        float counts_Jy_HI=0;
		for (xc=0;xc<RAsize;xc++){
			for (yc=0;yc<DECsize;yc++){
				outcube_tot[xc][yc]=outcube_tot[xc][yc]*ntomass/pixtopc*pctoJy*JytoK*fabs(cdelt3);
				for (vc=0;vc<VELsize;vc++){
					// NB. this flux will be the same as the data only if scale3=1
					// because it is in mJy/beam and not mJy/beam*km/s!
					//outcube[xc][yc][vc]=outcube[xc][yc][vc]*ntomass/pixtopc*pctoJy*JytoK;
                    outcube[xc][yc][vc]=outcube[xc][yc][vc]/pixtopc*pctoJy*JytoK;
                    counts_Jy_HI+=outcube[xc][yc][vc];
                    
				}
			}
		}
        float ntomassr1=M_HI_halo/2.36e5/dist/dist/fabs(cdelt3)/4/log(2)*PI*8.12;
        float rrrr=ntomassr1/(counts_Jy_HI);
        for (xc=0;xc<RAsize;xc++)
        for (yc=0;yc<DECsize;yc++)
            for (vc=0;vc<VELsize;vc++)
                outcube[xc][yc][vc]=outcube[xc][yc][vc]*rrrr;
	}

	else {
		/* Milky Way */
		// counts_HI has already been calculated
		ntomass=M_HI_halo/counts_HI; // counts -> Mo
		//cout<<"counts_HI = "<<counts_HI<<endl;
		for (xc=0;xc<RAsize;xc++){
			for (yc=0;yc<DECsize;yc++){
				outcube_tot[xc][yc]=outcube_tot[xc][yc]*ntomass*pctoJy*JytoK*fabs(cdelt3);
				for (vc=0;vc<VELsize;vc++){
					outcube[xc][yc][vc]=outcube[xc][yc][vc]*ntomass*pctoJy*JytoK;
				}
			}
		}
	}
	
	// density and velocity maps
	for (xc=0;xc<rotsurf_nR;xc++){
		for (yc=0;yc<rotsurf_nz; yc++){
			out_dens[xc][yc]=out_dens[xc][yc]*ntomass/rotsurf_grid*pctoJy*JytoK;
			for (vc=0;vc<rotsurf_nvel; vc++){
				out_v_phi[xc][yc][vc]=out_v_phi[xc][yc][vc]*ntomass
				/rotsurf_grid*pctoJy*JytoK;
				out_v_R[xc][yc][vc]=out_v_R[xc][yc][vc]*ntomass
				/rotsurf_grid*pctoJy*JytoK;
			}
		}
	}

	M_HI_above = counts_HI_above*ntomass;	//HI mass above .2 kpc		
	M_halo_above = counts_halo_above*ntomass; //HI+HII mass above .2 kpc
	M_accreted_above = counts_accreted_above*ntomass; //pristine gas accreted above .2 kpc

	//ESTIMATE THE FRACTION OF THE PRISTINE MATERIAL WITH A CERTAIN VELOCITY			
	for(int i=0; i<nvR; i++) {m_vR[i]/=counts_accreted_above;}
	for(int i=0; i<nvphi; i++) {m_vphi[i]/=counts_accreted_above;}
	//ofstream myvR("vR_distrib.txt");
	//for(int i=0; i<nvR; i++) myvR<<vR[i]<<" "<<m_vR[i]<<endl;
	//myvR.close();
	//ofstream myvphi("vphi_distrib.txt");
	//for(int i=0; i<nvphi; i++) myvphi<<vphi[i]<<" "<<m_vphi[i]<<endl;
	//myvphi.close();

	//COMPARISON BETWEEN THE AZIMUTHAL VELOCITY OF THE INFALLING MATERIAL AND THAT OF THE DISC			
	//ofstream LZ("Lz_comparison.txt");
	//LZ<<"# Radius   mean velocity   vhalo/vdisk   n_cloud"<<endl;  	
	//for(int i=0; i<nmyR; i++) if(npart[i]>0) LZ<<myR[i]<<" "<<MVrot[i]/M[i]<<" "<<MVrot_ratio[i]/M[i]<<" "<<npart[i]<<endl;
	//LZ.close();

	/*    questo no, perche' il tipo di calcolo dipende se external o
	 internal projection
	 counts_HI=0.;
	 for (xc=0;xc<RAsize;xc++){
	 for (yc=0;yc<DECsize;yc++){
	 counts_HI+=outcube_tot[xc][yc];
	 }
	 }
	 */
	//if(verbose)cout<<endl;
	//if(verbose)printf("HI mass above %.1f kpc = %.2e Mo\n",z_above_fin, M_HI_above);
	//if(verbose)printf("HI+HII mass above %.1f kpc = %.2e Mo\n",z_above_fin, M_halo_above);
	//if(verbose)printf("Accreted mass above %.1f kpc = %.2e Mo\n",z_above_fin, M_accreted_above);
	//if(verbose)cout<<"Percentage of pristine gas in the halo = "<<1e2*M_accreted_above/M_halo_above<<" percent"<<endl;
	//cout<<"Percentage of pristine gas in the local halo = "<<1e2*counts_accreted_local/counts_halo_local<<" percent"<<endl;
	//if(verbose)cout<<endl;
	total_energy_above=total_energy_above*ntomass*CONV2/YR;
	//if(verbose)printf("Total energy input to shoot particles above %.1f kpc: %.2e erg/s = %.1e SN/yr\n", z_above_fin, total_energy_above, total_energy_above*YR/1.e51);


	/*
	 Normalization of the outflow/inflow rates 
	 */
	//if (iter_mode == 'y')
	io.NormAndPlot(ntomass*Weight_witharms/Weight_noarms,
				   mean_v_k,
				   M_tot_pbc);
	
	accretion = io.get_GlobalAR();
	HImass    = M_HI_halo;

	for (int hi=0;hi<nR1-1;hi++) {
		for (int hk=0;hk<nz1-1;hk++) {
			dLz_hot_ik1[hi][hk]=dLz_hot_ik1[hi][hk]*ntomass; //kpc^2/Myr Myr^-1
			dvphi_hot_ik1[hi][hk]=dvphi_hot_ik1[hi][hk]*ntomass/CONV1; //km/s Myr^-1
			if (dvphi_hot_ik1[hi][hk] > 0) {
				//	  cout << R_i1[hi] << " " << z_k1[hk] << " " << dLz_hot_ik1[hi][hk] << " " << dvphi_hot_ik1[hi][hk] << "\n";
			}
		}
	}
	//    print(dLz_hot_ik1);
	if (strcoll(effects,"drag") ==0) {
		Lz_hot_tot=Lz_hot_tot*ntomass;
		cout << "Angular momentum transfer (to the halo) rate " << Lz_hot_tot 
		<< " Mo kpc^2 Myr^-2" << endl;
		printHeaderDeltaV(cdelt1, cdelt2, kpctograd);
		//writefits_2D("Delta_v_phi.fits", &dvphi_hot_ik1[0][0], nR1,nz1, "header.tab");
		printHeaderDeltaL(cdelt1, cdelt2, kpctograd);
		//writefits_2D("Delta_L.fits", &dLz_hot_ik1[0][0], nR1,nz1, "header.tab");
	}

	/* 
	 ADDING DISC COMPONENT
	 
	 ********* MOVE OUTSIDE MAIN PROGRAM? ********
	 */

	if ((add_disk== 'y') && (strcoll(projectionType,"hvcs") !=0)){
		if(verbose)cout << "Adding a disc component with mass = " << M_disc << " Mo"<<endl;
		if(verbose)cout << "Note: mass of HI disc computed from HI distribution is "<<qromo(HIsurf_R,Rmin_disc,Rmax_disc,'p',1.e-6)*1e6<<" Msun"<<endl;;  
		R=Rmin_disc;
		deltaR = 0.5;
		int intR = int(Rmin_disc);
		ii=0;
		
		//additional cube and matrix for the disc. It's needed if you are using spiral arms!
		Array3D<float> outcube_disc(RAsize,DECsize,VELsize,0.f);
		Array2D<float> outcube_tot_disc(RAsize,DECsize,0.f);
		
		/* 
		 Cycle in R (for disc component) 
		 */
		findValue(inputfile, FunctionPar[0], "h_disk"); 
		findValue(inputfile, FunctionPar[1], "Rstart_flare"); 
		double weight_noarms = 0;
		double weight_witharms = 0;
		
		while (R < Rmax_disc-deltaR){
			R=Rmin_disc+deltaR*ii;
			if(intR!=int(R))
			{
				intR = int(R);
				if(verbose)cout<<"Working at R = "<<intR<<" kpc"<<endl;
			}
			R1=R;
			R2=R+deltaR;
			R_ave=(R1+R2)/2.;
			np_phi=(int)(2.*PI*R_ave/deltaR+0.5);
		
			np_disc=n_part*deltaR;
			if (np_disc > 1000) {np_disc=250;}
			if (np_disc <= 100) {np_disc=250;}
			//np_disc = 5000; //change if you want
			//cloudsize = 25; //change if you want
			if(R_ave > Rstart_flare)
			{
				ch_disk = h_disk*exp((R_ave-Rstart_flare)/Rs_flare);
			}
			else
			{
				ch_disk = h_disk;
			}
			cn_disc=HIsurf(R_ave)/(2.*ch_disk);
			cp_disc=0;
			//counts_disc=qromo(HIsurf_R,Rmin_disc,Rmax_disc,'p',1.e-6);
			counts_disc=qromo(HIsurf_R_hs,Rmin_disc,Rmax_disc,'p',1.e-6);
			ntomass=M_disc/counts_disc;
			
			/* Cycle in particles (for disc component) */
			while (cp_disc < np_disc) {
				/* 1) radius */
				R=ran0(&seed)*(R2-R1)+R1;
				//csigma_v = sigma_v+dsigma_v/2.-dsigma_v/Rmax_disc*R;
				if(R>R_Sun) {csigma_v = sigma_v;}
				else{csigma_v = sigma_v + dsigma_v*(1.-R/R_Sun);}
				/* 2) azimuthal angle */
				phi=0;
				/* 3) vertical position */
				z=fabs(ran_exp(seed,ch_disk,0.)); //cout<<z<<endl;
				/* 4) radial velocity */
				v_R=ran_gau(seed ,csigma_v);
				/* 5) azimuthal velocity */
				v0_phi=sqrt(R_ave*getValue(F_R1,R_i1,z_k1,(R_ave),0.,1))/CONV1;
				v_phi=v0_phi+ran_gau(seed ,csigma_v);
				//			cout << R_ave << " " << ch_disk << endl;
				/* 6) vertical velocity */
				v_z=ran_gau(seed ,csigma_v);
				
				/* Cycle in azimuth (for disc component) */
				//Assuming a constant overdensity of 5 for the arms in the disc 
				a1_overdens = 5;
				a2_overdens = 5;
				a3_overdens = 5;
				a4_overdens = 5;
				
				for (jj=0; jj<np_phi; jj++){
					int sign=1;
					/* Mirrowing below the plane */
					for (int kk2=0; kk2<2; kk2++){
						
						/* Projections (for disc) */
						if (kk2 == 1) {sign=-1;}
						double phi_disc = phi+ran0(&seed)*360.;
						if(phi_disc>360) phi_disc-=360;
						double zwarp = z0_warp(R, phi_disc, create_warp);  
						
						projections proj_disc(R,
											  phi_disc,
											  z*sign+zwarp,
											  v_R,
											  v_phi,
											  v_z*sign);
						arms_weight = 0;
						//SPIRAL ARMS (in the disc)
						if(create_arms=='y' && R>=Rstart_arms)
						{
							//if(phi_disc>360) phi_disc-=360;
							for(int k=0;k<n_R_arms;k++)
							{
								R_overdensity1[k] = R_logspiral(phi_disc+k*360, a1_pitch, a1_amp);
								R_overdensity2[k] = R_logspiral(phi_disc+k*360, a2_pitch, a2_amp);
								R_overdensity3[k] = R_logspiral(phi_disc+k*360, a3_pitch, a3_amp);
								R_overdensity4[k] = R_logspiral(phi_disc+k*360, a4_pitch, a4_amp);
							}
							for(int k=0;k<n_R_arms;k++) 
							{
								arms_weight += (a1_overdens-1.)*exp(-(R-R_overdensity1[k])*(R-R_overdensity1[k])/(2.*a1_thickness*a1_thickness))
								+  (a2_overdens-1.)*exp(-(R-R_overdensity2[k])*(R-R_overdensity2[k])/(2.*a2_thickness*a2_thickness))
								+  (a3_overdens-1.)*exp(-(R-R_overdensity3[k])*(R-R_overdensity3[k])/(2.*a3_thickness*a3_thickness))
								+  (a4_overdens-1.)*exp(-(R-R_overdensity4[k])*(R-R_overdensity4[k])/(2.*a4_thickness*a4_thickness));
							}
						}
						
						/*  External projection */
						if (strcoll(projectionType,"external") ==0) {		
							proj_disc.external_view(incl,PA,xi,yi,vrad);
							/* Building the final cube (disc) */
							xc=int(RAsize/2+xi/fabs(cdelt1)*kpctograd+x_shift+0.5)-1; 
							// in this way the centre is 127
							yc=int(DECsize/2+yi/cdelt2*kpctograd+y_shift+0.5)-1;
							vc=int((vrad+v_shift)/cdelt3+0.5); 
							if (xc < RAsize && xc >= 0 && 
								yc < DECsize && yc >=0 && 
								vc < VELsize && vc >= 0) {
								if (mirror_ra =='y') {xc=RAsize-xc-1;}
								if (mirror_dec =='y') {yc=DECsize-yc-1;}
								outcube_disc[xc][yc][vc]+=
								(1.+arms_weight)*cn_disc/np_disc*deltaR*deltaR*ntomass/pixtopc*pctoJy*JytoK/2.;
								outcube_tot_disc[xc][yc]+=
								(1.+arms_weight)*cn_disc/np_disc*deltaR*deltaR*ntomass/pixtopc*pctoJy*JytoK/2.*fabs(cdelt3);
								weight_noarms += cn_disc;
								weight_witharms += cn_disc*(1.+arms_weight);
							}
						}
						/* Internal projection (DISC) */
						if ((strcoll(projectionType,"internal") ==0) || 
							(strcoll(projectionType,"hvcs") ==0)) { 
							proj_disc.internal_view(R_Sun,v0_phi_Rsun,l,b,d,vlos); // output: l(degree),b(degree),vlos;
							//if (R>=20) {cout<<"(R,phi_disc,l) = ("<<R<<","<<phi_disc<<","<<l<<")"<<endl;}
							xc=int(RAsize/2+l/fabs(cdelt1)+0.5);
							yc=int(DECsize/2+b/cdelt2+0.5);
							vc=int(VELsize/2-vlos/cdelt3+0.5);
							
							//		cout << z <<" "<< b <<" "<< d <<" "<< endl;
							
							
							/* Building the final cube */
							if (xc < RAsize && xc >= 0 
								&& yc < DECsize && yc >=0 
								&& vc < VELsize && vc >= 0) {
								if (d>distance_cut){
									pixtopc=fabs(cdelt1*cdelt2)/RAD/RAD*d*d*1.e6; 
									// pix -> pc^2
									if (d > RAD/sqrt(abs(cdelt1*cdelt2))*cloudsize/1000.){
										//		      vc=(int)(vc+ran_gau(seed,sigma_v));
										nPixels = (int)(1./cos(b/RAD)/fabs(cdelt1)+ran0(&seed)); 
										for (int p=0; p< nPixels; p++){
											pxc = (int)(xc - nPixels/2 + p + 0.5);
											if (pxc >= RAsize) pxc=pxc-RAsize;
											if (pxc < 0) pxc=pxc+RAsize; 
											outcube_disc[pxc][yc][vc]+=
											(1.+arms_weight)*cn_disc/np_disc*deltaR*deltaR*ntomass/pixtopc
											*pctoJy/2*JytoK/nPixels;
											outcube_tot_disc[pxc][yc]+=
											(1.+arms_weight)*cn_disc/np_disc*deltaR*deltaR*ntomass/pixtopc
											*pctoJy/2*JytoK*fabs(cdelt3)/nPixels;
											weight_noarms += cn_disc;
											weight_witharms += cn_disc*(1.+arms_weight);
										}
									}
									else { // resolved clouds
										pixsize1=(int)(atan(cloudsize/d/1000.)*RAD/fabs(cdelt1)+.5);
										pixsize2=(int)(atan(cloudsize/d/1000.)*RAD/cdelt2+.5);
										for (int p=0; p<pixsize1*pixsize2; p++) {
											//			cout << pixsize1 << " " << xc << " " << yc << endl;
											// randomization in each direction
											pxc=(int)(xc+ran_gau(seed,pixsize1/3.));
											pyc=(int)(yc+ran_gau(seed ,pixsize2/3.));
											pvc=(int)(vc+ran_gau(seed ,csigma_v/cdelt3));
											//			cout << xc << " " << yc << endl;
											if (pxc >= RAsize) pxc=pxc-RAsize;
											if (pxc < 0) pxc=pxc+RAsize; 
											if (pyc >= DECsize)
											{	
												pyc = 2*DECsize - pyc;
												if(pxc<RAsize/2) {pxc+=RAsize/2;}
												else{pxc-=RAsize/2;}
											}
											if (pyc < 0)
											{
												pyc = abs(pyc);
												if(pxc<RAsize/2) {pxc+=RAsize/2;}
												else{pxc-=RAsize/2;}
											} 
											if (pvc >= VELsize) pyc=VELsize;
											if (pvc < 0) pyc=0; 
											//		      cout << d << " " << cloudsize << " " 
											//			   << pixsize1 << endl;
											// effect of curvature: number of pixels as a function of b
											// note that differently from below here the latitude is recalculate at the exact position of pyc
											nPixels = (int)(1. / cos( (pyc - DECsize/2) * cdelt2 /RAD) / fabs(cdelt1) + ran0(&seed));
											for (int pp=0; pp< nPixels; pp++){
												ppxc = (int)(pxc - nPixels/2 + pp + 0.5);
												if (ppxc >= RAsize) ppxc=ppxc-RAsize;
												if (ppxc < 0) ppxc=ppxc+RAsize; 
												outcube_disc[ppxc][pyc][pvc]+=
												(1.+arms_weight)*cn_disc/np_disc*deltaR*deltaR*ntomass/pixtopc
												*pctoJy/2*JytoK
												/pixsize1/pixsize2*(ran_gau(seed,.3)+1)/nPixels;
												outcube_tot_disc[ppxc][pyc]+=
												(1.+arms_weight)*cn_disc/np_disc*deltaR*deltaR*ntomass/pixtopc
												*pctoJy/2*JytoK*fabs(cdelt3)
												/pixsize1/pixsize2*(ran_gau(seed,.3)+1)/nPixels;
												weight_noarms += cn_disc;
												weight_witharms += cn_disc*(1.+arms_weight);
											}
										}
									}
								}
							}
						}
					}
				}
				// One could leave out the contribution to out_XXX by the disk 
				// (negligible above 0.6 kpc)
				R_c=int(fabs(R/rotsurf_grid)+0.5);
				z_c=int(fabs(z/rotsurf_grid)+0.5);
				v_phi_c=int(v_phi+0.5);
				v_R_c=int(v_R+0.5)+rotsurf_nvel/2;
				if (R_c < rotsurf_nR && R_c >= 0 
					&& z_c < rotsurf_nz && z_c >=0 
					&& v_phi_c < rotsurf_nvel && v_phi_c >= 0) {
					out_v_phi[R_c][z_c][v_phi_c]+=
					(1.+arms_weight)*cn_disc/np_disc*rotsurf_grid*ntomass*pctoJy*JytoK/2.;
					out_v_R[R_c][z_c][v_R_c]+=
					(1.+arms_weight)*cn_disc/np_disc*rotsurf_grid*ntomass*pctoJy*JytoK/2.;
					out_dens[R_c][z_c]+=
					(1.+arms_weight)*cn_disc/np_disc*rotsurf_grid*ntomass*pctoJy*JytoK/2.;
				}
				
				cp_disc++;
			}
			ii++;
		}
		//normalize the disc and add it to the HI cube
		for(int i=0;i<RAsize;i++) for(int j=0;j<DECsize; j++) for(int k=0;k<VELsize; k++)
		{
			outcube_disc[i][j][k]*=weight_noarms/weight_witharms;
			outcube[i][j][k]+=outcube_disc[i][j][k];
		}
		//the same for the total map
		for(int i=0;i<RAsize;i++) for(int j=0;j<DECsize; j++)
		{
			outcube_tot_disc[i][j]*=weight_noarms/weight_witharms;
			outcube_tot[i][j]+=outcube_tot_disc[i][j];
		}
	}
	/* END ADDING DISC */

	PlotCurves(out_v_phi,
			   out_v_phi_wm,
			   out_v_R,
			   out_dens,
			   kpctograd,
			   rotsurf_grid,
			   rotsurf_nvel,
			   rotsurf_max,
			   clip_wm,
			   Rmin_gradient,
			   Rmax_gradient,
			   zmax_gradient,
			   make_fit,
			   incl);

	/*    printHeader(rotsurf_nR, rotsurf_nz, rotsurf_nvel, 
	 rotsurf_grid, rotsurf_grid, 2, 
	 0, 0, 0,
	 projectionType, "first", Jy2K);
	 writefits_3D("guip.fits", &out_v_phi[0][0][0], 
	 rotsurf_nR, rotsurf_nz, rotsurf_nvel, "header.tab");
	 */
         
	/* total counts */
	if (Jy2K == 'y')
	{
		if(verbose) printf("Total flux = %f K km/s \n", total(outcube));
		//	printf("Total flux = %f K km/s \n", total(outcube_tot));
	}
	else    
		if(verbose) printf("Total flux = %f Jy km/s \n", 
			   total(outcube_tot)*fabs(cdelt1*cdelt2*3600*3600)/beam_x/beam_y);
                 
		
	/*
	 HANNING SMOOTHING (do it BEFORE the cut!)
	 */			
	if (hanning== 'y') 
	{
		outcube.Hanning(3);
	}

	if(projAC=='y' && strcoll(projectionType,"internal") ==0)
	{
		double temp[RAsize];
		myarray::double3d tempcube(RAsize,DECsize,VELsize);
		for(int k=0; k<VELsize; k++) for(int j=0; j<DECsize; j++) for(int i=0;i<RAsize;i++) tempcube(i,j,k)=outcube[i][j][k];
		for (int k=0; k<VELsize; k++)
		{
			for (int j=0; j<DECsize; j++)
			{
				for (int i=0; i<=RAsize/2; i++) temp[i]=tempcube(i+RAsize/2,j,k);
				for (int i=RAsize/2 +1; i<RAsize; i++) temp[i]=tempcube(i-RAsize/2,j,k);
				for (int i=0;i<RAsize;i++) outcube[i][j][k] = temp[i]; 
			}
		}
		crval1 = 360;
	}

	/* 
	 SMOOTHING 
	 */
	Array3D<float> smoothedOutcube(RAsize,DECsize,VELsize,0.f);
	Array2D<float> smoothedOutcubeTot(RAsize,DECsize,0.f);

	if (smooth == 'y')
	{
		if (strcoll(projectionType,"external") ==0)	
		{
			myarray::double3d tempcube(RAsize,DECsize,VELsize);
			for(int i=0;i<RAsize;i++) for(int j=0;j<DECsize;j++) for(int k=0;k<VELsize;k++) tempcube(i,j,k) = outcube[i][j][k];
			fitsutil::smooth_cube(tempcube,RA,DEC,VEL,fabs(cdelt1)*3600.,fabs(cdelt2)*3600.,beam_x,beam_y,5);
			for(int i=0;i<RAsize;i++) for(int j=0;j<DECsize;j++) for(int k=0;k<VELsize;k++) smoothedOutcube[i][j][k] = tempcube(i,j,k);
			
			//smoothedOutcube=outcube.smooth2D(fabs(cdelt1)*3600.,fabs(cdelt2)*3600.,beam_x,beam_y);
			smoothedOutcubeTot=outcube_tot.smooth2D(fabs(cdelt1)*3600.,fabs(cdelt2)*3600.,beam_x,beam_y);
		}
		else 
		{
			double temp[RAsize];
			myarray::double3d tempcube(RAsize,DECsize,VELsize);
			for(int k=0; k<VELsize; k++) for(int j=0; j<DECsize; j++) for(int i=0;i<RAsize;i++) tempcube(i,j,k)=outcube[i][j][k];
			if(crval1!=360)
			{
				for (int k=0; k<VELsize; k++)
				{
					for (int j=0; j<DECsize; j++)
					{
						for (int i=0; i<=RAsize/2; i++) temp[i]=tempcube(i+RAsize/2,j,k);
						for (int i=RAsize/2 +1; i<RAsize; i++) temp[i]=tempcube(i-RAsize/2,j,k);
						for (int i=0;i<RAsize;i++) outcube[i][j][k] = temp[i]; 
					}
				} 
				
			}
			//smoothing
			for(int i=0;i<RAsize;i++) for(int j=0;j<DECsize;j++) for(int k=0;k<VELsize;k++) tempcube(i,j,k) = outcube[i][j][k];
			fitsutil::smooth_allsky_cube(tempcube,RA,DEC,VEL,fabs(cdelt1),fabs(cdelt2),beam_x,beam_y,3);
			for(int i=0;i<RAsize;i++) for(int j=0;j<DECsize;j++) for(int k=0;k<VELsize;k++) smoothedOutcube[i][j][k] = tempcube(i,j,k);
			
			if(projAC=='n')
			{
				//back to the original projection (longitude axis from 180 to -180)
				for(int k=0; k<VELsize; k++) for(int j=0; j<DECsize; j++) for(int i=0;i<RAsize;i++) tempcube(i,j,k)=smoothedOutcube[i][j][k];
				for (int k=0; k<VELsize; k++)
				{
					for (int j=0; j<DECsize; j++)
					{
						for (int i=0; i<=RAsize/2; i++) temp[i]=tempcube(i+RAsize/2,j,k);
						for (int i=RAsize/2 +1; i<RAsize; i++) temp[i]=tempcube(i-RAsize/2,j,k);
						for (int i=0;i<RAsize;i++) smoothedOutcube[i][j][k] = temp[i]; 
					}
				}
			}
			else{crval1=360;}
			//smoothedOutcube=outcube.smooth2D(fabs(cdelt1)*3600.,fabs(cdelt2)*3600.,beam_x,beam_y);
			smoothedOutcubeTot=outcube_tot.smooth2D(fabs(cdelt1),fabs(cdelt2),beam_x,beam_y);
		}
	}
	else
	{
		if(strcoll(projectionType,"internal")==0 && projAC=='y')
		{
			crval1 = 360;
			//reorganizing the cube
			double temp[RAsize];
			myarray::double3d tempcube(RAsize,DECsize,VELsize);
			for(int k=0; k<VELsize; k++) for(int j=0; j<DECsize; j++) for(int i=0;i<RAsize;i++) tempcube(i,j,k)=outcube[i][j][k];
			for (int k=0; k<VELsize; k++)
			{
				for (int j=0; j<DECsize; j++)
				{
					for (int i=0; i<=RAsize/2; i++) temp[i]=tempcube(i+RAsize/2,j,k);
					for (int i=RAsize/2 +1; i<RAsize; i++) temp[i]=tempcube(i-RAsize/2,j,k);
					for (int i=0;i<RAsize;i++) outcube[i][j][k] = temp[i]; 
				}
			} 
		}
	}
	
	//simulate a velocity dispersion
	/*
	int kmin, kmax;
	double halo_vdisp = kick_angle;
	if(verbose) cout<<"Including a velocity dispersion of "<<halo_vdisp<<" km/s."<<endl;
	double vdisp_convol = sqrt(halo_vdisp*halo_vdisp - cdelt3*cdelt3);
	if(verbose) cout<<"Convolving the cube with a gaussian of rms "<<vdisp_convol<<" km/s..."<<endl;
	double HIline[VELsize], myvalue;
	for(int i=0;i<RAsize;i++) for(int j=0;j<DECsize;j++)
	{
		for(int k=0;k<VELsize;k++) HIline[k] = smoothedOutcube[i][j][k];
		for(int k=0;k<VELsize;k++) 
		{
			kmin = k-3*vdisp_convol/cdelt3-0.5;
			kmax = k+3*vdisp_convol/cdelt3+0.5;
			if(kmin<0) kmin = 0;
			if(kmax>VELsize-1) kmax=VELsize-1;
			myvalue = 0;
			for(int kk=kmin; kk<=kmax; kk++) myvalue+=HIline[kk]*myutil::gaussian(VEL[k]-VEL[kk],vdisp_convol);
			smoothedOutcube[i][j][k] = myvalue*cdelt3; 
		}
	}
	*/
	
	//this defeat the purpose of calculating the total map BEFORE, right?
	if (smooth == 'y') 
	{
		for(int i=0; i<RAsize;i++) for(int j=0;j<DECsize;j++)
		{
			double value = 0;
			for(int k=0;k<VELsize;k++) if(!isnan(smoothedOutcube[i][j][k])) value += smoothedOutcube[i][j][k];
			smoothedOutcubeTot[i][j] = value;
		}
	}
	else
	{
		for(int i=0; i<RAsize;i++) for(int j=0;j<DECsize;j++)
		{
			double value = 0;
			for(int k=0;k<VELsize;k++) if(!isnan(outcube[i][j][k])) value += outcube[i][j][k];
			outcube_tot[i][j] = value;
		}
	}
		
	/*
	 WRITING TOTAL MAP
	 */
	if (strcoll(projectionType,"external") ==0)
		printHeaderTot(RAsize, DECsize, cdelt1, cdelt2, 
					   crval1,crval2,projectionType, Jy2K);
	if (strcoll(projectionType,"internal") ==0)
		printHeader2D_AllSky(RAsize, DECsize, 
							 cdelt1, cdelt2, 
							 1, 1, 
							 crval1,crval2,"K*km/s");

	if(iter_mode=='n')
	{
		//if (smooth == 'y')
			//writefits_2D("outcube_tot.fits", &smoothedOutcubeTot[0][0], RAsize, DECsize, "header.tab");
		//else
			//writefits_2D("outcube_tot.fits", &outcube_tot[0][0], RAsize, DECsize, "header.tab");  
	}

	//      P-V ALONG MAJOR AXIS
	for (xc=0;xc<RAsize;xc++)
		for (vc=0;vc<VELsize;vc++)
			maxax[RAsize-1-xc][vc]=outcube[xc][DECsize/2][vc];

	/*
	 Writing fits file of the P-V
	 */
	/*
	if (strcoll(projectionType,"external") == 0){
		printHeaderMax(RAsize, VELsize, cdelt1, cdelt3,crval1,crval3);//NB was crval2!!
		writefits_2D("maxax.fits", &maxax[0][0], RAsize, VELsize, "header.tab");
	}
	 */
	/* 
	 WRITING DATA CUBE
	 */
	if (iter_mode=='n') {//replace with n
		if(verbose)
		{
			cout << "Writing final cube" <<endl;
			cout<<"X_SIZE="<<RAsize<<endl;
			cout<<"Y_SIZE="<<DECsize<<endl;
			cout<<"Z_SIZE="<<VELsize<<endl;
			cout<<"CRVAL1 = "<<crval1<<endl;
			cout<<"CRVAL2 = "<<crval2<<endl;
			cout<<"CRVAL3 = "<<crval3<<endl;
			cout<<"CDELT1 ="<<cdelt1<<endl;
			cout<<"CDELT2 ="<<cdelt2<<endl;
			cout<<"CDELT3 ="<<cdelt3<<endl;
		}
		
        //cout<<projectionType<<endl;
        
        
        if (strcoll(projectionType,"external") ==0)
			printHeader(RAsize, DECsize, VELsize, 
						cdelt1, cdelt2, cdelt3, 
						crval1, crval2, crval3, 
						projectionType, "middle", Jy2K);
		if (strcoll(projectionType,"internal") ==0)
			printHeader3D_AllSky(RAsize, DECsize, VELsize,
								 cdelt1, cdelt2, cdelt3,
								 1, 1, 1,
								 crval1, crval2, crval3, "K");
		
		//if (smooth == 'y')
			//writefits_3D(outputfilename, &smoothedOutcube[0][0][0], RAsize, DECsize, VELsize, "header.tab");
		//else
			//writefits_3D(outputfilename, &outcube[0][0][0], RAsize, DECsize, VELsize, "header.tab");
	}

	if(which_latitude=='p')
	{
		for(int i=0;i<RAsize;i++) for(int j=0;j<DECsize/2;j++) for(int k=0;k<VELsize;k++)
		{
			smoothedOutcube[i][j][k] = 0;
			outcube[i][j][k] = 0;
		}
	}
	if(which_latitude=='n')
	{
		for(int i=0;i<RAsize;i++) for(int j=DECsize/2;j<DECsize;j++) for(int k=0;k<VELsize;k++)
		{
			smoothedOutcube[i][j][k] = 0;
			outcube[i][j][k] = 0;
		}
	}
    
    //cout<<"fdsfas"<<sizeof(outcube)<<endl;
    //return outcube;
    if(smooth=='y') {return smoothedOutcube;}
	else {return outcube;}
   
    
}
