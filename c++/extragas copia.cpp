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

#include "extragas.h"
#include "iterations.h" // keep it after fits_io.h
#include "arrays3d.h"
#include "fitsio.h" 
#include "milkyway.h"
#include "fitsutil.h"
#include "myutilities.h"
#include "profiles.h"
#include "nr3.h"
#include <vector>
#include <fstream>
#define TOKMS 1.e-3

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

extern char inputfile[];
Array1D<float> FunctionPar(10);

//initialise surface density HI and H2
std::vector<float> RHI,DHI,RH2,DH2;
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

int extragas_model()
{
	cout.precision(6);
	
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
	long int seed, seed2;
	float Lz_hot_tot;
	//
	float xi,yi,vrad,l,b,d,vlos;
	float vlos_max,vlos_min; // used for HVCs
	float R_out_above, R_in_above; 
	
	/* 
	 DECLARATIONS: OUTPUTS
	 */
	
	int n_iter1, n_iter2, n_iter3;
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
	//residual
	double res_absolute, res_squared;
	
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
	
	/* FOR MOVIE 
	 string imageName;
	 list<Image> imageList; 
	 int movieIncr=0;
	 */	
	
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
	FILE* f_hvcs_0=NULL;
	FILE* f_hvcs_p1=NULL;
	FILE* f_hvcs_n1=NULL;
	FILE* f_hvcs_p2=NULL;
	FILE* f_hvcs_n2=NULL;
	FILE* f_hvcs_axes=NULL;
	FILE * f_matrix1;
	FILE * f_matrix2;
	FILE* f_xspec=NULL;
	
	cout << "\n";
	cout << "EXTRAGAS\n";
	cout << "Integrating orbits of extraplanar gas particles in a galactic potential\n";
	cout << "\n";
	cout << "Reading parameters from file 'extragas.in'..." << endl;
	
	/* 
	 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	 READING PARAMETERS from INPUT FILE 
	 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	 */
	/* Input Files */
	char filename1[60], filename2[60], surfdens_HI[60], surfdens_H2[60], datacube_file[90], outputfilename[50], labfile[60], resfile[60], compare_to[50], absorption_list[50]; 
	findValue(inputfile, filename1, "inputfile1");
	findValue(inputfile, filename2, "inputfile2");
	findValue(inputfile, resfile, "residuals_file");
	findValue(inputfile, surfdens_HI, "surfdens_HI");
	findValue(inputfile, surfdens_H2, "surfdens_H2");
	findValue(inputfile, datacube_file, "datacube");
	findValue(inputfile, outputfilename, "outputfilename");
	findValue(inputfile, compare_to, "compare_to");
	findValue(inputfile, absorption_list, "absorption_list");
	int low_radii; findValue(inputfile, low_radii, "low_radii"); 
	int subdim1; findValue(inputfile, subdim1, "subdim1"); 
	int subdim2; findValue(inputfile, subdim2, "subdim2"); 
		
	
	/*
	 GALACTIC PARAMETERS
	 */
	
	//reading HI and H2 surface density profile from files
	int lenght;
	char* buffer;
	int number_of_raws=0;
	
	//HI
	cout<<"Reading HI surface density from file: "<<surfdens_HI<<endl;
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
	cout<<"Reading H2 surface density from file: "<<surfdens_H2<<endl;
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
	float crval1; findValue(inputfile, crval1, "CRVAL1"); 
	float crval2; findValue(inputfile, crval2, "CRVAL2"); 
	float crval3; findValue(inputfile, crval3, "CRVAL3"); 
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
	//datacube containing the distance of the emission with respect to the Sun. Useful in the inner projection!
	myarray::array3d<double> distcube; double* distcube_ptr = NULL;
	//datacube containing the life-time of HI clouds  
	myarray::array3d<double> timecube; double* timecube_ptr = NULL;
	//datacube containing the kick radius of HI clouds  
	myarray::array3d<double> rkickcube; double* rkickcube_ptr = NULL;
	//datacube containing the pristine fraction associated to each voxel 
	myarray::array3d<double> metcube; double* metcube_ptr = NULL;
	//datacube containing the ejection radius of 
	char dcube; findValue(inputfile, dcube, "dcube");
	char tcube; findValue(inputfile, tcube, "tcube");
	char rcube; findValue(inputfile, rcube, "rcube");
	char Zcube; findValue(inputfile, Zcube, "Zcube");
	double Z_corona; findValue(inputfile, Z_corona, "Z_corona");
	double Z_disk; findValue(inputfile, Z_disk, "Z_disk");
	if(dcube=='y') 
	{
		distcube_ptr = new double [RAsize*DECsize*VELsize];
		distcube.init(distcube_ptr,RAsize,DECsize,VELsize);
		for(int i=0;i<RAsize;i++) for(int j=0;j<DECsize;j++) for(int k=0;k<VELsize;k++) distcube(i,j,k)=0; 
	}
	if(tcube=='y') 
	{
		timecube_ptr = new double [RAsize*DECsize*VELsize];
		timecube.init(timecube_ptr,RAsize,DECsize,VELsize);
		for(int i=0;i<RAsize;i++) for(int j=0;j<DECsize;j++) for(int k=0;k<VELsize;k++) timecube(i,j,k)=0; 
	}
	if(rcube=='y') 
	{
		rkickcube_ptr = new double [RAsize*DECsize*VELsize];
		rkickcube.init(rkickcube_ptr,RAsize,DECsize,VELsize);
		for(int i=0;i<RAsize;i++) for(int j=0;j<DECsize;j++) for(int k=0;k<VELsize;k++) rkickcube(i,j,k)=0; 
	}
	if(Zcube=='y') 
	{
		metcube_ptr = new double [RAsize*DECsize*VELsize];
		metcube.init(metcube_ptr,RAsize,DECsize,VELsize);
		for(int i=0;i<RAsize;i++) for(int j=0;j<DECsize;j++) for(int k=0;k<VELsize;k++) metcube(i,j,k)=0; 
	}

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
	if(projAC=='y' && iter_mode =='y' && (strcoll(projectionType,"internal") ==0))
	{
		cout<<"WARNING! ARE YOU SURE THAT THE DATACUBE USES THE ANTI-CENTER PROJECTION?"<<endl;
	}
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
	float h_v_k; findValue(inputfile, h_v_k, "h_v_k"); 
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
	char visible[10]; findValue(inputfile, visible, "visible_part"); // option are: all, z_max, R_max, partial
	//fraction of the rising orbit (assumed 1/2 of the full one) not visible. A linear run of the vertical velocity is assumed
	float ion_frac; findValue(inputfile, ion_frac, "ion_frac");
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
	if(single_event=='y')
	{
		double pitch;
		if(which_arm==1)      {a1_overdens=1e8;a2_overdens=1;a3_overdens=1;a4_overdens=1;pitch=a1_pitch;}
		else if (which_arm==2){a1_overdens=1;a2_overdens=1e8;a3_overdens=1;a4_overdens=1;pitch=a2_pitch;}
		else if (which_arm==3){a1_overdens=1;a2_overdens=1;a3_overdens=1e8;a4_overdens=1;pitch=a3_pitch;}
		else if (which_arm==4){a1_overdens=1;a2_overdens=1;a3_overdens=1;a4_overdens=1e8;pitch=a4_pitch;}
		RminSF = R_burst; RmaxSF = RminSF+0.0001;
		deltaR = dtoR_logspiral(Rtod_logspiral(R_burst,pitch)+0.5*arm_extent,pitch)-dtoR_logspiral(Rtod_logspiral(R_burst,pitch)-0.5*arm_extent,pitch);
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
	/* Normalization of accretion at vrel=75 km/s */
	float accr_norm;  findValue(inputfile, accr_norm, "accr_norm_75"); 
	accr_norm*=1e-3;
	
	/* Parameters of accretion */
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
	
	cout<<"ACCRETION PATTERN = "<<accr_vel<<endl;
	
	/* Drag timescale at vrel=75 km/s (in Myrs) */
	float t_drag; findValue(inputfile, t_drag, "t_drag_75"); 
	
	/* DEFINING OUTPUT MATRIXES */
	Array3D<float> outcube(RAsize,DECsize,VELsize,0.f);
	myarray::double3d datacube; 
	double* datacube_ptr = NULL;
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
	
	int dim_3 = 3; 
	
	//  float guihvc[RAsize][DECsize];
	
	/*
	 READING ITERATIVE PARAMETERS
	 */
	findValue(inputfile, par1, "iter_par1");
	findValue(inputfile, par2, "iter_par2");
	findValue(inputfile, par3, "iter_par3");
	char resume; //allow to continue a previous iteration, using "resfile" for previous data 
	findValue(inputfile, resume, "resume");
	char fullResiduals; findValue(inputfile, fullResiduals, "fullResiduals"); // if y=normal residual fitting, if n=minimizes only some sight lines
	
	iterations Iter(par1,par2,par3,
					n_iter1,n_iter2,n_iter3,
					param1_0, delta_param1,
					param2_0, delta_param2,
					param3_0, delta_param3);
	
	if (iter_mode == 'y')
	{
		printf("Number of iterations %i x %i x %i = %i \n", n_iter1, n_iter2, n_iter3, n_iter1*n_iter2*n_iter3);
	}
	else 
	{
		n_iter1=1;
		n_iter2=1;
		n_iter3=1;
	}
	
	//ITERATIVE MODE
	int n_start1 = 0, n_start2 = 0, n_start3 = 0;
	double dataflux = 0.;
	ofstream myresidual;
	std::vector<double> l_abs, b_abs, v_abs;
	std::vector<char> type_abs;
	std::vector<string> id_abs, ref_abs;
	if(iter_mode=='y')
	{
		if(strcmp(compare_to, "cube") == 0)
		{
			//READING THE DATACUBE
			cout<<"Reading the datacube "<<datacube_file<<endl;
			fitsutil::read_cube(datacube_file,RA,DEC,VEL,datacube_ptr);
			datacube.init(datacube_ptr, RA.size(), DEC.size(), VEL.size());
			if(RAsize!=RA.size() || DECsize!=DEC.size() || VELsize!=VEL.size())
			{
				cout<<"FATAL ERROR! The datacube and the modelcube have different sizes"<<endl;
				cout<<"Sizes of the datacube  : "<<RA.size()<<","<<DEC.size()<<","<<VEL.size()<<endl;
				cout<<"Sizes of the modelcube : "<<RAsize<<","<<DECsize<<","<<VELsize<<endl;
				cout<<"Seriously man, you can do better"<<endl; 
				return EXIT_FAILURE;
			}
			if (strcoll(projectionType,"external") ==0)
			{
				for(int i=0;i<RAsize;i++) for(int j=0;j<DECsize;j++) for(int k=0;k<VELsize;k++) dataflux += datacube(i,j,k);
			}
			else if (strcoll(projectionType,"internal") ==0)
			{
				double cosb;
				for(int i=0;i<RAsize;i++) for(int j=0;j<DECsize;j++) for(int k=0;k<VELsize;k++)
				{
					if(!isnan(datacube(i,j,k)))
					{
						cosb = cos(DEC[j]*DEGTORAD);
						dataflux += datacube(i,j,k)*cosb;
					}
				}
			}
			else {dataflux = NAN;}
			dataflux *= fabs(cdelt1*cdelt2*cdelt3); 
			cout<<"Data flux = "<<dataflux<<endl;
		}
		else if(strcmp(compare_to, "list") == 0)
		{
			cout<<"Reading the absorption list "<<absorption_list<<endl;
			int n_abs=myutil::nlines_in(absorption_list);
			cout<<"    "<<n_abs<<" absorptions found"<<endl;
			l_abs.resize(n_abs); b_abs.resize(n_abs); v_abs.resize(n_abs); type_abs.resize(n_abs); id_abs.resize(n_abs); ref_abs.resize(n_abs);
			ifstream absdata(absorption_list);
			for(int i=0;i<n_abs;i++) absdata>>id_abs[i]>>l_abs[i]>>b_abs[i]>>v_abs[i]>>type_abs[i]>>ref_abs[i];
			absdata.close();
			//for(int i=0;i<n_abs;i++) cout<<id_abs[i]<<" "<<l_abs[i]<<" "<<b_abs[i]<<" "<<v_abs[i]<<" "<<type_abs[i]<<" "<<ref_abs[i]<<endl;
		}
		else
		{
			cout<<"ERROR! 'compare_to' keyword can be only 'cube' or 'list'!"<<endl;
			return EXIT_FAILURE;
		}
		if(resume == 'n')
		{
			//Opening residual file (only if confirmed)
			char choice;
			cout<<"Do you want to open a new file "<<resfile<<" ?"<<endl;
			cout<<"This will delete any previous version of this file"<<endl;
			cout<<"Continue? (y/n) ";
			cin>>choice;
			if(choice != 'y') {cout<<endl<<"Ok, bye bye!"<<endl; return EXIT_SUCCESS;}
			else {myresidual.open(resfile);}
		}
		else
		{
			//resuming residual file
			cout<<"Opening the already existing file "<<resfile<<endl;
			
			int lenght;
			char* buffer;
			int number_of_raws=0;
			ifstream data(resfile);
			data.seekg(0, ios::end);
			lenght = data.tellg();
			data.seekg(0, ios::beg);
			buffer = new char [lenght];
			data.read (buffer, lenght);
			for(int i=0; i<lenght; i++) if(buffer[i]==0x0A) number_of_raws++;
			data.close();
			cout<<"Number of iteration already done: "<<number_of_raws<<endl;
			n_start1 = number_of_raws/(n_iter2*n_iter3);
			n_start2 = (number_of_raws%(n_iter2*n_iter3))/n_iter3;
			n_start3 = number_of_raws%n_iter3;
			myresidual.open(resfile, ios_base::app);
		}
	}
	
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
		cout << "Using potential from file: " << filename1 << endl;
		nFiles=1;}
	else {
		f_matrix2=fopen (filename2,"r");
		cout << "Using potential from files: " << filename1 << " and " << filename2 << endl;
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
	printf ("Reading %d x %d matrix \n",nR1,nz1);
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
		if_stop=0;
		fscanf (f_matrix2, "%lf %lf %lf %lf %lf ", &fx1,&fx3,&pot,&fx4,&fx6);
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
		printf ("Reading %d x %d matrix \n",nR2,nz2);
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
	
	if (strcoll(effects,"accretion")==0 || strcoll(effects,"accrdrag")==0){
		if (strcoll(accr_dens,"build") ==0){
			Accretion smooth_accretion(accr_r0, 0, accr_v_phi0, 0, accr_sigma_v0);
			printf("Building accretion pattern \n");
			smooth_accretion.Patterns(cold_rho_ik,
									  cold_v_R_ik,
									  cold_v_phi_ik,
									  cold_v_z_ik,
									  R_i1,z_k1,pot_ik1,F_R1,F_z1);
			
			printHeaderAccr(nR1, nz1, stepR1, stepz1);
			writefits_2D_DP("accr_density.fits", &cold_rho_ik[0][0], nR1, nz1, "header.tab");
			printHeaderAccr(nR1, nz1, stepR1, stepz1);
			writefits_2D_DP("accr_v_R.fits", &cold_v_R_ik[0][0], nR1, nz1, "header.tab");
			printHeaderAccr(nR1, nz1, stepR1, stepz1);
			writefits_2D_DP("accr_v_phi.fits", &cold_v_phi_ik[0][0], nR1, nz1, "header.tab");
			printHeaderAccr(nR1, nz1, stepR1, stepz1);
			writefits_2D_DP("accr_v_z.fits", &cold_v_z_ik[0][0], nR1, nz1, "header.tab");
		}
	}
	/*
	 ------------------
	 0) ITERATION CYCLE
	 ------------------
	 */
	
	int global_iter= n_start1*n_iter2*n_iter3 + n_start2*n_iter3 + n_start3;
	float mypar1, mypar2, mypar3;
	stop=n;
	
	for (int iter1=n_start1 ; iter1<n_iter1; iter1++){
		for (int iter2=n_start2 ; iter2<n_iter2; iter2++){	
			for (int iter3=n_start3 ; iter3<n_iter3; iter3++){
			
			time_t start,end;
			time (&start);	
			
			int ii=0, jj=0, vkk=0;
			io.back_to_zero();
			if (iter_mode=='y'){
				global_iter+=1;
				
				if (strcoll(par1,"M_HI_halo") ==0) {M_HI_halo=param1_0+delta_param1*iter1; mypar1=M_HI_halo;}
				if (strcoll(par2,"M_HI_halo") ==0) {M_HI_halo=param2_0+delta_param2*iter2; mypar2=M_HI_halo;}
				if (strcoll(par3,"M_HI_halo") ==0) {M_HI_halo=param3_0+delta_param3*iter3; mypar3=M_HI_halo;}
				M_disc = M_HI_tot - M_HI_halo;
				if (strcoll(par1,"h_v_k") ==0) {h_v_k=param1_0+delta_param1*iter1;mypar1=h_v_k;}
				if (strcoll(par2,"h_v_k") ==0) {h_v_k=param2_0+delta_param2*iter2;mypar2=h_v_k;}
				if (strcoll(par3,"h_v_k") ==0) {h_v_k=param3_0+delta_param3*iter3;mypar3=h_v_k;}
				if (strcoll(par1,"ion_frac") ==0) {ion_frac=param1_0+delta_param1*iter1;mypar1=ion_frac;}
				if (strcoll(par2,"ion_frac") ==0) {ion_frac=param2_0+delta_param2*iter2;mypar2=ion_frac;}
				if (strcoll(par3,"ion_frac") ==0) {ion_frac=param3_0+delta_param3*iter3;mypar3=ion_frac;}
				if (strcoll(par1,"n_part") ==0) {n_part=param1_0+delta_param1*iter1;mypar1=n_part;}
				if (strcoll(par2,"n_part") ==0) {n_part=param2_0+delta_param2*iter2;mypar2=n_part;}
				if (strcoll(par3,"n_part") ==0) {n_part=param3_0+delta_param3*iter3;mypar3=n_part;}
				if (strcoll(par1,"showt_min") ==0) {showt_min=param1_0+delta_param1*iter1;mypar1=showt_min;}
				if (strcoll(par2,"showt_min") ==0) {showt_min=param2_0+delta_param2*iter2;mypar2=showt_min;}
				if (strcoll(par3,"showt_min") ==0) {showt_min=param3_0+delta_param3*iter3;mypar3=showt_min;}
				if (strcoll(par1,"showt_max") ==0) {showt_max=param1_0+delta_param1*iter1;mypar1=showt_max;}
				if (strcoll(par2,"showt_max") ==0) {showt_max=param2_0+delta_param2*iter2;mypar2=showt_max;}
				if (strcoll(par3,"showt_max") ==0) {showt_max=param3_0+delta_param3*iter3;mypar3=showt_max;}
				if (strcoll(par1,"showt2_min") ==0) {showt2_min=param1_0+delta_param1*iter1; mypar1=showt2_min; showt2_max=showt2_min + showt_max - showt_min;}
				if (strcoll(par2,"showt2_min") ==0) {showt2_min=param2_0+delta_param2*iter2; mypar2=showt2_min; showt2_max=showt2_min + showt_max - showt_min;}
				if (strcoll(par3,"showt2_min") ==0) {showt2_min=param3_0+delta_param3*iter3; mypar3=showt2_min; showt2_max=showt2_min + showt_max - showt_min;}
				if (strcoll(par1,"v_k_thres") ==0) 
				{v_k_thres=param1_0+delta_param1*iter1;mypar1=v_k_thres;}
				if (strcoll(par2,"v_k_thres") ==0) 
				{v_k_thres=param2_0+delta_param2*iter2;mypar2=v_k_thres;}
				if (strcoll(par3,"v_k_thres") ==0) 
				{v_k_thres=param3_0+delta_param3*iter3;mypar3=v_k_thres;}
				if (strcoll(par1,"RmaxSF") ==0) {RmaxSF=param1_0+delta_param1*iter1; RminSF=RmaxSF-delta_param1; mypar1=RmaxSF;}
				if (strcoll(par2,"RmaxSF") ==0) {RmaxSF=param2_0+delta_param2*iter2; RminSF=RmaxSF-delta_param2; mypar2=RmaxSF;}
				if (strcoll(par3,"RmaxSF") ==0) {RmaxSF=param3_0+delta_param3*iter3; RminSF=RmaxSF-delta_param3; mypar3=RmaxSF;}
				if (strcoll(par1,"arm_section") ==0) {RminSF=param1_0+delta_param1*iter1; RmaxSF=RminSF+0.01; deltaR=delta_param1; nRSF=(int)((RmaxSF-RminSF)/deltaR)+1; mypar1=RminSF;}
				if (strcoll(par2,"arm_section") ==0) {RminSF=param2_0+delta_param2*iter2; RmaxSF=RminSF+0.01; deltaR=delta_param2; nRSF=(int)((RmaxSF-RminSF)/deltaR)+1; mypar2=RminSF;}
				if (strcoll(par3,"arm_section") ==0) {RminSF=param3_0+delta_param3*iter3; RmaxSF=RminSF+0.01; deltaR=delta_param3; nRSF=(int)((RmaxSF-RminSF)/deltaR)+1; mypar3=RminSF;}
				if (strcoll(par1,"disk_section") ==0) {RminSF=param1_0+delta_param1*iter1; RmaxSF=RminSF+0.01; deltaR=delta_param1; nRSF=(int)((RmaxSF-RminSF)/deltaR)+1; mypar1=RminSF;}
				if (strcoll(par2,"disk_section") ==0) {RminSF=param2_0+delta_param2*iter2; RmaxSF=RminSF+0.01; deltaR=delta_param2; nRSF=(int)((RmaxSF-RminSF)/deltaR)+1; mypar2=RminSF;}
				if (strcoll(par3,"disk_section") ==0) {RminSF=param3_0+delta_param3*iter3; RmaxSF=RminSF+0.01; deltaR=delta_param3; nRSF=(int)((RmaxSF-RminSF)/deltaR)+1; mypar3=RminSF;}
				if (strcoll(par1,"disk_section")==0 || strcoll(par2,"disk_section")==0 || strcoll(par3,"disk_section")==0) create_arms='n';
				if (strcoll(par1,"arm_section")==0 || strcoll(par2,"arm_section")==0 || strcoll(par3,"arm_section")==0) 
			
				if (strcoll(par1,"alpha_v_k") ==0)
				{alpha_v_k=param1_0+delta_param1*iter1;mypar1=alpha_v_k;}
				if (strcoll(par2,"alpha_v_k") ==0)
				{alpha_v_k=param2_0+delta_param2*iter2;mypar2=alpha_v_k;}
				if (strcoll(par3,"alpha_v_k") ==0)
				{alpha_v_k=param3_0+delta_param3*iter3;mypar3=alpha_v_k;}
				if (strcoll(par1,"v_k_max") ==0) {v_k_max=param1_0+delta_param1*iter1;mypar1=v_k_max;}
				if (strcoll(par2,"v_k_max") ==0) {v_k_max=param2_0+delta_param2*iter2;mypar2=v_k_max;}
				if (strcoll(par3,"v_k_max") ==0) {v_k_max=param3_0+delta_param3*iter3;mypar3=v_k_max;}
				if (strcoll(par1,"gammaSF") ==0)
				{gammaSF=param1_0+delta_param1*iter1;mypar1=gammaSF;}
				if (strcoll(par2,"gammaSF") ==0)
				{gammaSF=param2_0+delta_param2*iter2;mypar2=gammaSF;}
				if (strcoll(par3,"gammaSF") ==0)
				{gammaSF=param3_0+delta_param3*iter3;mypar3=gammaSF;}
				if (strcoll(par1,"v0_hot") ==0) {v0_hot=param1_0+delta_param1*iter1;mypar1=v0_hot;}
				if (strcoll(par2,"v0_hot") ==0) {v0_hot=param2_0+delta_param2*iter2;mypar2=v0_hot;}
				if (strcoll(par3,"v0_hot") ==0) {v0_hot=param3_0+delta_param3*iter3;mypar3=v0_hot;}
				if (strcoll(par1,"K_D") ==0) {K_D=param1_0+delta_param1*iter1;mypar1=K_D;}
				if (strcoll(par2,"K_D") ==0) {K_D=param2_0+delta_param2*iter2;mypar2=K_D;}
				if (strcoll(par3,"K_D") ==0) {K_D=param3_0+delta_param3*iter3;mypar3=K_D;}
				if (strcoll(par1,"accr_norm") ==0) {accr_norm=param1_0+delta_param1*iter1;mypar1=accr_norm; accr_norm*=1e-3;}
				if (strcoll(par2,"accr_norm") ==0) {accr_norm=param2_0+delta_param2*iter2;mypar2=accr_norm; accr_norm*=1e-3;}
				if (strcoll(par3,"accr_norm") ==0) {accr_norm=param3_0+delta_param3*iter3;mypar3=accr_norm; accr_norm*=1e-3;}
				//accr_norm=pow(10,accr_norm);
				if (strcoll(par1,"accr_z") ==0) {accr_z_fin=param1_0+delta_param1*iter1;mypar1=accr_z_fin;}
				if (strcoll(par2,"accr_z") ==0) {accr_z_fin=param2_0+delta_param2*iter2;mypar2=accr_z_fin;}
				if (strcoll(par3,"accr_z") ==0) {accr_z_fin=param3_0+delta_param3*iter3;mypar3=accr_z_fin;}
				if (strcoll(par1,"t_drag") ==0) {t_drag=param1_0+delta_param1*iter1;mypar1=t_drag;}
				if (strcoll(par2,"t_drag") ==0) {t_drag=param2_0+delta_param2*iter2;mypar2=t_drag;}
				if (strcoll(par3,"t_drag") ==0) {t_drag=param3_0+delta_param3*iter3;mypar3=t_drag;}
				if (strcoll(par1,"kick_angle") ==0) {kick_angle=param1_0+delta_param1*iter1;mypar1=kick_angle;}
				if (strcoll(par2,"kick_angle") ==0) {kick_angle=param2_0+delta_param2*iter2;mypar2=kick_angle;}
				if (strcoll(par3,"kick_angle") ==0) {kick_angle=param3_0+delta_param3*iter3;mypar3=kick_angle;}
				
				cout<<endl<<endl;
				cout<<"ITERATION NUMBER: "<<global_iter<<" / "<<n_iter1*n_iter2*n_iter3<<endl;
				cout<<"Values of parameters:"<<endl;
				cout<<"   "<<par1<<" = "<<param1_0+delta_param1*iter1<<endl;   	
				cout<<"   "<<par2<<" = "<<param2_0+delta_param2*iter2<<endl;
				cout<<"   "<<par3<<" = "<<param3_0+delta_param3*iter3<<endl;
			}
			
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
				f_orbit=fopen("orbit.dat", "w");
				f_veldistr=fopen("extragas_veldistr.dat", "w");
				f_veldistr2=fopen("extragas_veldistr2.dat", "w");
				f_times=fopen("extragas_time.dat", "w");
				if ((RminSF-deltaR) > 0) {
					fprintf(f_times,"0 0 0 0 %f %.2e \n", RminSF-deltaR, 
							2*PI*(RminSF-deltaR)/sqrt
							((RminSF-deltaR)*
							 getValue(F_R1,R_i1,z_k1,(RminSF-deltaR),0.,1))*1e6);
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
			
			/*
			 ---------------------------
			 1) MAIN CYCLE IN R (RADIUS)
			 ---------------------------
			 */
			cout <<endl<< "Starting integration" << endl;
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
						  t_drag,
						  accr_norm
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
						  t_drag,
						  accr_norm
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
				seed=time(NULL)*((ii+1)*(ii+1)*654136);
				//	R=RminSF+ran0(&seed)*(RmaxSF-RminSF);
				R=RminSF+deltaR*ii;
				cout <<"Working at Radius R= "<<R << endl;
				/*
				 calculate number of particles in this bin
				 */
				R1=R-deltaR/2.;
				R2=R+deltaR/2.;
				R_ave=R;
				np_phi=(int)((2.*PI*R/deltaR+0.5));
				np_vk=n_part;
				cn_dot=HoleExpSchmidt(R_ave)/
				qromo(HoleExpSchmidt,RminSF,RmaxSF,'p',1.e-6); /* normalized n_dot*/ 
				
				/* 
				 ---------------------------
				 2) CYCLE IN KICK VELOCITIES
				 ---------------------------
				 */
				while (cp_vk < np_vk) {
					if (iter_mode == 'n'){//replace with n
						freopen("orbit.dat","w",f_orbit);
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
					seed=time(NULL)*(cp_vk+1)*(ii+2);
					/* 1) Initial radius */
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
					
					v_R=ran_gau(seed,csigma_v);
					/* 5) azimuthal velocity */
					v0_phi=sqrt(R_ave*getValue(F_R1,R_i1,z_k1,(R_ave),0.,1))/CONV1;
					v_phi=v0_phi+ran_gau(seed/2,csigma_v);
					/* 6) vertical velocity */
					v_z=ran_gau(seed/3,csigma_v); // N.B. here v_z was = to fabs(ran_gau()), why?
					
					/* 
					 Kick velocities
					 */
					/* direction distributions */
					phi_k=ran0(&seed)*2*PI;
					theta_k=ran0(&seed)*kick_angle*PI/180.; // linear in theta
					//theta_k=kick_angle*PI/180.;
					
					/* modulus distribution */
					// *************** MOVE TO FUNCTIONS.CPP??? *****************
					v_k=v_k_max;
					int v_k_iter=1;
					while (v_k >= v_k_max) {
						seed=time(NULL)*(cp_vk+1)*(ii+2)*v_k_iter;
						if (strcmp(v_k_distr, "gau") == 0){
							v_k=(v_k_0+ran_gau(seed/4,h_v_k))*pow(cos(theta_k),alpha_cos);
						} 
						else {
							if (strcmp(v_k_distr, "exp") == 0){
								v_k=
								ran_exp(seed/4,h_v_k,v_k_thres)*pow(cos(theta_k),alpha_cos);
							} 
							else {
								if (strcmp(v_k_distr, "pow") == 0){
									// NB. here I substituted v_k_cut with v_k_thres
									v_k=
									ran_pow(seed/4,alpha_v_k,csigma_v,v_k_thres)*
									pow(cos(theta_k),alpha_cos);
								} 
								else {
									if (strcmp(v_k_distr, "lin") == 0){
										v_k=ran_lin(seed/4,ab_v_k)*pow(cos(theta_k),alpha_cos);
									} 
									else {
										if (strcmp(v_k_distr, "con") == 0){
											v_k=h_v_k;
										}
									}
								}
							}
						}
						if (v_k > v_k_max)
							cout << "Kick velocity above limit, R = " << R <<" v_k = "<<v_k <<endl;
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
						fprintf(f_veldistr,"%f %f %f %f %f %f \n", theta_k*RAD, phi_k*RAD, v_k, v_R, v_phi, v_z);
						if (v_k > v_k_thres) {
							fprintf(f_veldistr2,"%f %f %f %f %f %f \n", theta_k*RAD, phi_k*RAD, v_k, v_R, v_phi, v_z);
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
							fprintf(f_orbit,"%.1f %.4f %.4f %.4f %.3f %.3f %.3f %.9f %.5f %.9f %.5f %.5f %.9f %.9f  %.9f 1 1 1 1 1 1 \n",
									t,x1,x2,x3,x4/CONV1,x5/CONV1,x6/CONV1,
									Etot0,Ltot0,L_z0,
									getValue(pot_ik1,R_i1,z_k1,x1,fabs(x3),1), 
									(cn_dot+cn_accr)/cn_dot,
									Etot_M0,Ltot_M0,L_z_M0);
					}
					
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
									if (strcoll(effects,"accretion") ==0 || strcoll(effects,"accrdrag")==0) {
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
																	cn_dot+cn_accr
																	)*dt;
										}
									}
									dyn2.rk4(t,dt,x1, x2, x3, x4, x5, x6,cn_dot+cn_accr);
									//		  dyn2.rk4(dt,x1, x2, x3, x4, x5, x6,cn_dot);
								}
								else { // nFiles == 1 or in File 1
									if (strcoll(effects,"accretion") ==0 || strcoll(effects,"accrdrag")==0) {
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
																	cn_dot+cn_accr
																	)*dt;
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
										if (strcoll(effects,"accretion") ==0 || strcoll(effects,"accrdrag")==0) {
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
																		cn_dot+cn_accr
																		)*dt;
											}
										}
										//cout << "refin " << x1 << " " << x2 << " " << x3 << " " << x4 << " " << x5 << " " << x6 << " " << cn_dot+cn_accr<< endl;
										dyn2.rk4(t,dt,x1, x2, x3, x4, x5, x6,cn_dot+cn_accr);
										//		    dyn2.rk4(dt,x1, x2, x3, x4, x5, x6,cn_dot);
									}
									else { // nFiles == 1, refined
										if (strcoll(effects,"accretion") ==0 || strcoll(effects,"accrdrag")==0) {
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
																		cn_dot+cn_accr
																		)*dt;
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
								if(velocity_gap=='y' && create_arms=='y') phigap = 0.06*t*pattern_speed - x2;
								
								//		cout << t << " " << x5/CONV1 << " " << x1
								//		  << " " << x3 << " " << x4/CONV1 << " " << x6/CONV1 << endl;
								
								t+=dt;
								
								/* Check if the new position is inside the potential grid */
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
									cout << "Warning: particles outside the potential grid! \n";
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
								
								if ((single_event=='y' && t>showt_min && t<=showt_max) || (single_event=='y' && second_burst=='y' && t>showt2_min && t<=showt2_max)){
									
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

											xx2=ran0(&seed2)*360.;
											if (kk2 == 1) {sign=-1;}
											// x3 and x6 ARE MULTIPLIED BY sign
											
											seed2=time(NULL)*(jj+1)*(ii+1)*(cp_vk+1);
											
											//
											//   External projection 
											
											if (strcoll(projectionType,"external") ==0) {
												
												double phi_halo = x2+xx2;
												projections proj(x1,
																 x2+xx2+phigap,
																 x3*sign,
																 x4/CONV1,
																 x5/CONV1,
																 x6/CONV1*sign);
												proj.external_view(incl, xi, yi, vrad);
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
														((strcmp(visible, "z_max") == 0 && x6 < 0) || 
														 (strcmp(visible, "z_max") == 0 && x3 < 0)) 
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
													cout << "Warning: particles outside the cube boundaries (xc=" << xc << "," << "yc=" << yc << "," << "vc=" << vc << ")" << endl;
													n_lost += PI*(R2*R2 - R1*R1);
													out_boundaries=true;
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
												double phi_halo = x2+xx2;
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
														((strcmp(visible, "z_max") == 0 && x6 < 0) || 
														 (strcmp(visible, "z_max") == 0 && x3 < 0)) 
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
																	pixsize1=(int)(atan(cloudsize/d/1000.)*RAD/fabs(cdelt1)+.5);
																	pixsize2=(int)(atan(cloudsize/d/1000.)*RAD/fabs(cdelt2)+.5);
																	
																	if(pixsize1<1 || pixsize2<1) cout<<"Warning!!!"<<endl;
																	for (int p=0; p<pixsize1*pixsize2; p++) 
																	{
																		seed=time(NULL)*((jj+1)*(p+1)*432143);
																		// randomization in each direction
																		pxc=(int)(xc+ran_gau(seed,pixsize1/3.));
																		pyc=(int)(yc+ran_gau(seed/3,pixsize2/3.));
																		pvc=(int)(vc+ran_gau(seed/2,csigma_v/cdelt3));
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
																		//cout<<"ran0="<<ran0(&seed)<<endl;
																		for (int pp=0; pp< nPixels; pp++)
																		{
																			ppxc = (int)(pxc - nPixels/2 + pp + 0.5);
																			if (ppxc >= RAsize) ppxc=ppxc-RAsize;
																			if (ppxc < 0) ppxc=ppxc+RAsize; 
																			float myran = (ran_gau(seed,.3)+1);
																			float value = (1.+arms_weight)*(cn_dot+cn_accr)	
																			*deltaR*deltaR*dt/dt_default/pixtopc
																			/pixsize1/pixsize2
																			*myran // WHY is this there? (randomization around 1)
																			/nPixels;
																			outcube[ppxc][pyc][pvc]+=value;
																			if(dcube=='y') distcube(ppxc,pyc,pvc)+=value*d;
																			if(tcube=='y') timecube(ppxc,pyc,pvc)+=value*t;
																			if(rcube=='y') rkickcube(ppxc,pyc,pvc)+=value*R_kick;
																			if(Zcube=='y') metcube(ppxc,pyc,pvc)+=(1.+arms_weight)*cn_accr*deltaR*deltaR*dt/dt_default/pixtopc/pixsize1/pixsize2*myran/nPixels;
																			outcube_tot[ppxc][pyc]+=(1.+arms_weight)*(cn_dot+cn_accr)
																			*deltaR*deltaR*dt/dt_default/pixtopc
																			*fabs(cdelt3)
																			/pixsize1/pixsize2*(ran_gau(seed,.3)+1)/nPixels;
																		}
																	}
																}
																else {
																	// effect of curvature: number of pixels as a function of b
																	seed=time(NULL)*((jj+1)*xc*432143);
																	nPixels = (int)(1./cos(b/RAD)/fabs(cdelt1)+ran0(&seed)); 
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
																			if(dcube=='y') distcube(pxc,yc,vc)+=value*d;
																			if(tcube=='y') timecube(pxc,yc,vc)+=value*t;
																			if(rcube=='y') rkickcube(pxc,yc,vc)+=value*R_kick;
																			if(Zcube=='y') metcube(pxc,yc,vc)+=(1.+arms_weight)*cn_accr
																				*deltaR*deltaR*dt/dt_default/pixtopc
																				/nPixels;
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
													cout << "Warning: particles outside the cube boundaries (xc="<<xc<<" yc="<<yc<<" vc=" << vc << ")" << endl;
													n_lost += PI*(R2*R2 - R1*R1);
													out_boundaries=true;
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
													cout << "Warning: particles outside the cube boundaries (xc="<<xc<<" yc="<<yc<<" vc=" << vc << ")" << endl;
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
												((strcmp(visible, "z_max") == 0 && x6 < 0) || 
												 (strcmp(visible, "z_max") == 0 && x3 < 0)) 
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
										fprintf(f_orbit,"%.1f %.4f %.4f %.4f %.3f %.3f %.3f %.9f %.5f %.9f %.5f %.5f %.9f %.5f %.9f %.9f %.9f %.9f %.9f %.9f %.9f \n",
												t,x1,x2,x3,x4/CONV1,x5/CONV1,x6/CONV1,
												Etot,Ltot,L_z, 
												getValue(pot_ik1,R_i1,z_k1,x1,fabs(x3),1), 
												(cn_dot+cn_accr)/cn_dot,
												Etot_M,Ltot_M,L_z_M,
												Etot/Etot0, Ltot/Ltot0, L_z/L_z0,
												Etot_M/Etot_M0, Ltot_M/Ltot_M0, L_z_M/L_z_M0);
										
										//cyltocar(x1,x2,x3,x4,x5,x6);
										//		    fprintf(f_orbit,"%.5f %.5f %.5f \n", getforce(type_pot,1,x1,x2,x3,0,0,0),getforce(type_pot,2,x1,x2,x3,0,0,0),getforce(type_pot,3,x1,x2,x3,0,0,0));
										//		    cartocyl(x1,x2,x3,x4,x5,x6);
									}
								}
							} 
							else {
								cout << "Warning: particles outside the potential grid!" << endl;
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
								printf("Maximum number of integrations exceeded : %i \n",n_dt-1);
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
						
						
						
						/***** MOVIE 
						 if (iter_mode == 'n' && make_movie == 'y'){
						 // Initialize image (ImageMagick) //
						 int IMsize=RAsize*DECsize*3;
						 float pixels[IMsize];
						 Image image(RAsize, DECsize, "RGB", FloatPixel, pixels);
						 float r, g, b, norma=max(outcube_tot);
						 for (int i=0; i<RAsize; i++){
						 for (int j=0; j<DECsize; j++){
						 r=0;
						 g=outcube_tot[i][j]/norma*sqrt(fabs(j-DECsize/2)+1);
						 b=0;
						 if (g > 1) g=1;
						 image.pixelColor(i,j, ColorRGB(r,g,b));
						 }
						 }
						 movieIncr++;
						 imageName="movie/image";
						 stringstream oss;
						 oss << imageName << movieIncr;
						 string output(oss.str());
						 imageName=output+".gif";
						 int zoomRAsize=RAsize+movieIncr;
						 int zoomDECsize=DECsize+movieIncr;
						 //	    image.zoom(Geometry(zoomRAsize,zoomDECsize));
						 //	    image.rotate(100-movieIncr/10.);
						 image.write(imageName);
						 readImages( &imageList, imageName );
						 }
						 *****/
						
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
								fclose(f_orbit);
								plot1(accr_norm);
								// these plotting routines are only for Fraternali&Binney2008
								//		plot0();
								//		plot0col();
								cout << "writing output plot file extragas.ps\n";
								//		cout << '\a';
							}
							fprintf(f_times,"%i %f %f %.2e %f %.2e \n", n_tau, R, v_k, t*1.e6, R_ave, 2*PI*R_ave/v0_phi/CONV1*1e6);
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
				cout << " number of sample particles : " << cp_halo << endl;
				//	cout<<" counts_HI : "<<counts_HI<<endl;
				total_energy_above+=(ctotal_energy_above);
				
				
			}
			//      fclose(f_orbit);
			/* END CYCLE IN R */
			//flux-weighted distance
			if(dcube=='y' && strcoll(projectionType,"internal") ==0)
			{
				for(int i=0;i<RAsize;i++) for(int j=0;j<DECsize;j++) for(int k=0;k<VELsize;k++) if(outcube[i][j][k]>0) {distcube(i,j,k)/=outcube[i][j][k];}
				for(int i=0;i<RAsize;i++) for(int j=0;j<DECsize;j++) for(int k=0;k<VELsize;k++) if(distcube(i,j,k)<=0 || outcube[i][j][k]<threshold_res) {distcube(i,j,k)=NAN;}
			}
			//flux-weighted lifetime
			if(tcube=='y' && strcoll(projectionType,"internal") ==0)
			{
				for(int i=0;i<RAsize;i++) for(int j=0;j<DECsize;j++) for(int k=0;k<VELsize;k++) if(outcube[i][j][k]>0) {timecube(i,j,k)/=outcube[i][j][k];}
				for(int i=0;i<RAsize;i++) for(int j=0;j<DECsize;j++) for(int k=0;k<VELsize;k++) if(timecube(i,j,k)<=0 || outcube[i][j][k]<threshold_res) {timecube(i,j,k)=NAN;}
			}
			//flux-weighted kick radii
			if(rcube=='y' && strcoll(projectionType,"internal") ==0)
			{
				for(int i=0;i<RAsize;i++) for(int j=0;j<DECsize;j++) for(int k=0;k<VELsize;k++) if(outcube[i][j][k]>0) {rkickcube(i,j,k)/=outcube[i][j][k];}
				for(int i=0;i<RAsize;i++) for(int j=0;j<DECsize;j++) for(int k=0;k<VELsize;k++) if(rkickcube(i,j,k)<=0 || outcube[i][j][k]<threshold_res) {rkickcube(i,j,k)=NAN;}
			}
			//metallicity 
			if(Zcube=='y' && strcoll(projectionType,"internal") ==0)
			{
				//pristine fraction
				for(int i=0;i<RAsize;i++) for(int j=0;j<DECsize;j++) for(int k=0;k<VELsize;k++) if(outcube[i][j][k]>0) {metcube(i,j,k)/=outcube[i][j][k];}
				for(int i=0;i<RAsize;i++) for(int j=0;j<DECsize;j++) for(int k=0;k<VELsize;k++) if(metcube(i,j,k)<0 || outcube[i][j][k]<threshold_res) {metcube(i,j,k)=NAN;}
				//convert to metallicity
				for(int i=0;i<RAsize;i++) for(int j=0;j<DECsize;j++) for(int k=0;k<VELsize;k++) {metcube(i,j,k) = Z_corona*metcube(i,j,k) + Z_disk*(1.-metcube(i,j,k));}
			}
				
				
			/* 
			 %%%%%%%%%%%%%%%%%%%%%%%%
			 OUTPUTS AND PLOTTING
			 %%%%%%%%%%%%%%%%%%%%%%%%
			 */
			/* 
			 Plotting kick velocity distribution and travel times
			 */
			if (iter_mode == 'n'){//replace with n
				fclose(f_veldistr);
				fclose(f_veldistr2);
				fclose(f_times);
				plot3();
				plot4();
			}
			
			np_fraction=1.*np_halo/nRSF/np_vk;
			mean_v_k=sqrt(mean_v_k/np_halo);
			cout << "Fraction of particles above velocity threshold : " << np_fraction << endl; 
			cout << "Mean travel time : " << tau/n_tau << " Myrs"<< endl;
			cout << "Mean travel time for clouds crossing z="<<z_above_in<<" kpc : " << tau_above/n_tau_above << " Myrs"<< endl;
			cout << "Mean kick velocity : " << mean_v_k << endl;
			cout << "Fraction of particles lost during the integration: " << 100.*n_lost/n_total<<" percent"<<endl; 
			
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
			
			cout<<"JytoK factor: "<<JytoK<<endl;
			
			
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
				pixtopc=fabs(cdelt1*cdelt2)/RAD/RAD*dist*dist*1.e12; 
				//	cout << ntomass << " " << pixtopc << " " << pctoJy << endl;
				// pix -> pc^2
				for (xc=0;xc<RAsize;xc++){
					for (yc=0;yc<DECsize;yc++){
						outcube_tot[xc][yc]=outcube_tot[xc][yc]*ntomass/pixtopc*pctoJy*JytoK*fabs(cdelt3);
						for (vc=0;vc<VELsize;vc++){
							// NB. this flux will be the same as the data only if scale3=1
							// because it is in mJy/beam and not mJy/beam*km/s!
							outcube[xc][yc][vc]=outcube[xc][yc][vc]*ntomass/pixtopc*pctoJy*JytoK;
						}
					}
				}
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
			ofstream myvR("vR_distrib.txt");
			for(int i=0; i<nvR; i++) myvR<<vR[i]<<" "<<m_vR[i]<<endl;
			myvR.close();
			ofstream myvphi("vphi_distrib.txt");
			for(int i=0; i<nvphi; i++) myvphi<<vphi[i]<<" "<<m_vphi[i]<<endl;
			myvphi.close();
			
			//COMPARISON BETWEEN THE AZIMUTHAL VELOCITY OF THE INFALLING MATERIAL AND THAT OF THE DISC			
			ofstream LZ("Lz_comparison.txt");
			LZ<<"# Radius   mean velocity   vhalo/vdisk   n_cloud"<<endl;  	
			for(int i=0; i<nmyR; i++) if(npart[i]>0) LZ<<myR[i]<<" "<<MVrot[i]/M[i]<<" "<<MVrot_ratio[i]/M[i]<<" "<<npart[i]<<endl;
			LZ.close();
			
			/*    questo no, perche' il tipo di calcolo dipende se external o
			 internal projection
			 counts_HI=0.;
			 for (xc=0;xc<RAsize;xc++){
			 for (yc=0;yc<DECsize;yc++){
			 counts_HI+=outcube_tot[xc][yc];
			 }
			 }
			 */
			cout<<endl;
			printf("HI mass above %.1f kpc = %.2e Mo\n",z_above_fin, M_HI_above);
			printf("HI+HII mass above %.1f kpc = %.2e Mo\n",z_above_fin, M_halo_above);
			printf("Accreted mass above %.1f kpc = %.2e Mo\n",z_above_fin, M_accreted_above);
			cout<<"Percentage of pristine gas in the halo = "<<1e2*M_accreted_above/M_halo_above<<" percent"<<endl;
			//cout<<"Percentage of pristine gas in the local halo = "<<1e2*counts_accreted_local/counts_halo_local<<" percent"<<endl;
			cout<<endl;
			total_energy_above=total_energy_above*ntomass*CONV2/YR;
			printf("Total energy input to shoot particles above %.1f kpc: %.2e erg/s = %.1e SN/yr\n", z_above_fin, total_energy_above, total_energy_above*YR/1.e51);

			
			/*
			 Normalization of the outflow/inflow rates 
			 */
			//if (iter_mode == 'y')
			io.NormAndPlot(ntomass*Weight_witharms/Weight_noarms,
						   mean_v_k,
						   M_tot_pbc);
			
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
				writefits_2D("Delta_v_phi.fits", &dvphi_hot_ik1[0][0], nR1,nz1, "header.tab");
				printHeaderDeltaL(cdelt1, cdelt2, kpctograd);
				writefits_2D("Delta_L.fits", &dLz_hot_ik1[0][0], nR1,nz1, "header.tab");
			}
			
			/* 
			 ADDING DISC COMPONENT
			 
			 ********* MOVE OUTSIDE MAIN PROGRAM? ********
			 */
			
			if ((add_disk== 'y') && (strcoll(projectionType,"hvcs") !=0)){
				cout << "Adding a disc component with mass = " << M_disc << " Mo"<<endl;
				cout << "Note: mass of HI disc computed from HI distribution is "<<qromo(HIsurf_R,Rmin_disc,Rmax_disc,'p',1.e-6)*1e6<<" Msun"<<endl;;  
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
						cout<<"Working at R = "<<intR<<" kpc"<<endl;
					}
					R1=R;
					R2=R+deltaR;
					R_ave=(R1+R2)/2.;
					np_phi=(int)(2.*PI*R_ave/deltaR+0.5);
				
					np_disc=n_part;
					if (n_part > 1000) {np_disc=250;}
					if (n_part <= 100) {np_disc=250;}
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
						seed=time(NULL)*(cp_disc+1)*(ii+1);
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
						v_R=ran_gau(seed/2,csigma_v);
						/* 5) azimuthal velocity */
						v0_phi=sqrt(R_ave*getValue(F_R1,R_i1,z_k1,(R_ave),0.,1))/CONV1;
						v_phi=v0_phi+ran_gau(seed/3,csigma_v);
						//			cout << R_ave << " " << ch_disk << endl;
						/* 6) vertical velocity */
						v_z=ran_gau(seed/4,csigma_v);
						
						/* Cycle in azimuth (for disc component) */
						//Assuming a constant overdensity of 2 for the arms in the disc 
						a1_overdens = 2;
						a2_overdens = 2;
						a3_overdens = 2;
						a4_overdens = 2;
						
						for (jj=0; jj<np_phi; jj++){
							int sign=1;
							/* Mirrowing below the plane */
							for (int kk2=0; kk2<2; kk2++){
								
								/* Projections (for disc) */
								if (kk2 == 1) {sign=-1;}
								seed2=time(NULL)*(ii+1)*(jj+1)*(cp_disc+1);
								double phi_disc = phi+ran0(&seed2)*360.;
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
									proj_disc.external_view(incl,xi,yi,vrad);
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
												seed=time(NULL)*((jj+1)*pxc*10001);
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
													seed=time(NULL)*((jj+1)*(p+1)*432143);
													//			cout << pixsize1 << " " << xc << " " << yc << endl;
													// randomization in each direction
													pxc=(int)(xc+ran_gau(seed,pixsize1/3.));
													pyc=(int)(yc+ran_gau(seed/3,pixsize2/3.));
													pvc=(int)(vc+ran_gau(seed/2,csigma_v/cdelt3));
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
				printf("Total flux = %f K km/s \n", total(outcube));
				//	printf("Total flux = %f K km/s \n", total(outcube_tot));
			}
			else 
				printf("Total flux = %f Jy km/s \n", 
					   total(outcube_tot)*fabs(cdelt1*cdelt2*3600*3600)/beam_x/beam_y);
			
				
			/*
			 HANNING SMOOTHING (do it BEFORE the cut!)
			 */			
			if (hanning== 'y') 
			{
				outcube.Hanning(3);
			}
			
			/*
			 CUT THE DISK EMISSION + EXTERNAL GALAXIES
			 */
			
			if (iter_mode=='y' && disk_cut == 'y' && strcmp(compare_to, "cube")==0)
			{
				if(projAC=='y')
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
				cout<<"Blanking the region of the disc (works only if data are blanked)"<<endl;
				for(int i=0; i<RAsize; i++)
				{
					for(int j=0; j<DECsize; j++)
					{
						for (int k=0; k<VELsize; k++)
						{
							if(isnan(datacube(i,j,k))) {outcube[i][j][k] = NAN;} //cout<<"NAN in ("<<i<<","<<j<<","<<k<<")"<<endl;}
						}
					}
				}
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
					fitsutil::smooth_cube(tempcube,RA,DEC,VEL,fabs(cdelt1)*3600.,fabs(cdelt2)*3600.,beam_x,beam_y,5,false);
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
					fitsutil::smooth_allsky_cube(tempcube,RA,DEC,VEL,fabs(cdelt1),fabs(cdelt2),beam_x,beam_y,3,false);
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
			if(strcoll(projectionType,"internal")==0 && projAC=='y' && dcube=='y')
			{
				crval1 = 360;
				//reorganizing the distance cube
				double temp[RAsize];
				myarray::double3d tempcube(RAsize,DECsize,VELsize);
				for(int k=0; k<VELsize; k++) for(int j=0; j<DECsize; j++) for(int i=0;i<RAsize;i++) tempcube(i,j,k)=distcube(i,j,k);
				for (int k=0; k<VELsize; k++)
				{
					for (int j=0; j<DECsize; j++)
					{
						for (int i=0; i<=RAsize/2; i++) temp[i]=tempcube(i+RAsize/2,j,k);
						for (int i=RAsize/2 +1; i<RAsize; i++) temp[i]=tempcube(i-RAsize/2,j,k);
						for (int i=0;i<RAsize;i++) distcube(i,j,k) = temp[i]; 
					}
				}
			}
			if(strcoll(projectionType,"internal")==0 && projAC=='y' && tcube=='y')
			{
				crval1 = 360;
				//reorganizing the distance cube
				double temp[RAsize];
				myarray::double3d tempcube(RAsize,DECsize,VELsize);
				for(int k=0; k<VELsize; k++) for(int j=0; j<DECsize; j++) for(int i=0;i<RAsize;i++) tempcube(i,j,k)=timecube(i,j,k);
				for (int k=0; k<VELsize; k++)
				{
					for (int j=0; j<DECsize; j++)
					{
						for (int i=0; i<=RAsize/2; i++) temp[i]=tempcube(i+RAsize/2,j,k);
						for (int i=RAsize/2 +1; i<RAsize; i++) temp[i]=tempcube(i-RAsize/2,j,k);
						for (int i=0;i<RAsize;i++) timecube(i,j,k) = temp[i]; 
					}
				}
			}
			if(strcoll(projectionType,"internal")==0 && projAC=='y' && rcube=='y')
			{
				crval1 = 360;
				//reorganizing the distance cube
				double temp[RAsize];
				myarray::double3d tempcube(RAsize,DECsize,VELsize);
				for(int k=0; k<VELsize; k++) for(int j=0; j<DECsize; j++) for(int i=0;i<RAsize;i++) tempcube(i,j,k)=rkickcube(i,j,k);
				for (int k=0; k<VELsize; k++)
				{
					for (int j=0; j<DECsize; j++)
					{
						for (int i=0; i<=RAsize/2; i++) temp[i]=tempcube(i+RAsize/2,j,k);
						for (int i=RAsize/2 +1; i<RAsize; i++) temp[i]=tempcube(i-RAsize/2,j,k);
						for (int i=0;i<RAsize;i++) rkickcube(i,j,k) = temp[i]; 
					}
				}
			}
			if(strcoll(projectionType,"internal")==0 && projAC=='y' && Zcube=='y')
			{
				crval1 = 360;
				//reorganizing the distance cube
				double temp[RAsize];
				myarray::double3d tempcube(RAsize,DECsize,VELsize);
				for(int k=0; k<VELsize; k++) for(int j=0; j<DECsize; j++) for(int i=0;i<RAsize;i++) tempcube(i,j,k)=metcube(i,j,k);
				for (int k=0; k<VELsize; k++)
				{
					for (int j=0; j<DECsize; j++)
					{
						for (int i=0; i<=RAsize/2; i++) temp[i]=tempcube(i+RAsize/2,j,k);
						for (int i=RAsize/2 +1; i<RAsize; i++) temp[i]=tempcube(i-RAsize/2,j,k);
						for (int i=0;i<RAsize;i++) metcube(i,j,k) = temp[i]; 
					}
				}
			}
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
			
			if (smooth == 'y')
				writefits_2D("outcube_tot.fits", &smoothedOutcubeTot[0][0], RAsize, DECsize, "header.tab");
			else
				writefits_2D("outcube_tot.fits", &outcube_tot[0][0], RAsize, DECsize, "header.tab");  
			
			
			//      P-V ALONG MAJOR AXIS
			for (xc=0;xc<RAsize;xc++)
				for (vc=0;vc<VELsize;vc++)
					maxax[RAsize-1-xc][vc]=outcube[xc][DECsize/2][vc];
			
			/*
			 Writing fits file of the P-V
			 */
			if (strcoll(projectionType,"external") == 0){
				printHeaderMax(RAsize, VELsize, cdelt1, cdelt3,crval1,crval3);//NB was crval2!!
				writefits_2D("maxax.fits", &maxax[0][0], RAsize, VELsize, "header.tab");
			}
			
			/* 
			 WRITING DATA CUBE
			 */
			if (iter_mode=='n') {//replace with n
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
				
				if (smooth == 'y')
					writefits_3D(outputfilename, &smoothedOutcube[0][0][0], RAsize, DECsize, VELsize, "header.tab");
				else
					writefits_3D(outputfilename, &outcube[0][0][0], RAsize, DECsize, VELsize, "header.tab");
				if(dcube=='y' && strcoll(projectionType,"internal") ==0) 
				{
					std::vector<double> RAinv(RAsize);
					if(projAC=='y') {for(int i=0;i<RAsize;i++) RAinv[i]=RA[i]+180;}
					else{for(int i=0;i<RAsize;i++) RAinv[i]=RA[i];}
					fitsutil::print_cube("mydistance.fits",RAinv,DEC,VEL,distcube_ptr,"GLON-CAR","GLAT-CAR","VELO-LSR","DEG","DEG","km/s","kpc");
					delete[] distcube_ptr; distcube_ptr = NULL;
				}
				if(tcube=='y' && strcoll(projectionType,"internal") ==0) 
				{
					std::vector<double> RAinv(RAsize);
					if(projAC=='y') {for(int i=0;i<RAsize;i++) RAinv[i]=RA[i]+180;}
					else{for(int i=0;i<RAsize;i++) RAinv[i]=RA[i];}
					fitsutil::print_cube("mytime.fits",RAinv,DEC,VEL,timecube_ptr,"GLON-CAR","GLAT-CAR","VELO-LSR","DEG","DEG","km/s","Myr");
					delete[] timecube_ptr; timecube_ptr = NULL;
				}
				if(rcube=='y' && strcoll(projectionType,"internal") ==0) 
				{
					std::vector<double> RAinv(RAsize);
					if(projAC=='y') {for(int i=0;i<RAsize;i++) RAinv[i]=RA[i]+180;}
					else{for(int i=0;i<RAsize;i++) RAinv[i]=RA[i];}
					fitsutil::print_cube("myRkick.fits",RAinv,DEC,VEL,rkickcube_ptr,"GLON-CAR","GLAT-CAR","VELO-LSR","DEG","DEG","km/s","kpc");
					delete[] rkickcube_ptr; rkickcube_ptr = NULL;
				}
				if(Zcube=='y' && strcoll(projectionType,"internal") ==0) 
				{
					std::vector<double> RAinv(RAsize);
					if(projAC=='y') {for(int i=0;i<RAsize;i++) RAinv[i]=RA[i]+180;}
					else{for(int i=0;i<RAsize;i++) RAinv[i]=RA[i];}
					fitsutil::print_cube("mymetallicity.fits",RAinv,DEC,VEL,metcube_ptr,"GLON-CAR","GLAT-CAR","VELO-LSR","DEG","DEG","km/s","Z0");
					delete[] metcube_ptr; metcube_ptr = NULL;
				}
			} 		
			/*
			 LONGITUDE GAS DISTRIBUTION
			*/
			/*	
			if(strcmp(projectionType, "internal")==0)
			{
				std::vector<double> gas_lon(RAsize);
				for(int i=0;i<RAsize;i++) gas_lon[i]=0;
				if (smooth == 'y')
				{
					for(int i=0;i<RAsize;i++)
					{
						double value = 0;
						for(int j=0;j<DECsize;j++) for(int k=0;k<VELsize;k++) value+=smoothedOutcube[i][j][k]*cos(DEC[j]*DEGTORAD);
						gas_lon[i]=value*fabs(cdelt1*cdelt2*cdelt3);
					}
				}
				else
				{
					for(int i=0;i<RAsize;i++)
					{
						double value = 0;
						for(int j=0;j<DECsize;j++) for(int k=0;k<VELsize;k++) value+=outcube[i][j][k]*cos(DEC[j]*DEGTORAD);
						gas_lon[i]=value*fabs(cdelt1*cdelt2*cdelt3);
					}
				}	
				double max_gas_lon = myutil::max(gas_lon);
				for(int i=0;i<RAsize;i++) gas_lon[i]/=max_gas_lon;
				ofstream data_lon("coldens_longitude.txt");
				for(int i=0;i<RAsize;i++) data_lon<<RA[i]<<setw(12)<<gas_lon[i]<<endl;
				data_lon.close();
			}
			*/
			/*
			 COMPARISON WITH OBSERVATIONS 
			 */
			Array2D<float> rescube_tot(RAsize,DECsize, 0.f);
			Array3D<float> rescube(RAsize,DECsize,VELsize, 0.f);
			double fluxratio,eta;
			int nabs_model; 
			if (iter_mode=='y')
			{
				double res_1 = 0., res_2 = 0., res_3 = 0.;
				
				//COMPARISON WITH A DATACUBE
				if(strcmp(compare_to, "cube") == 0)
				{
					if(normalize_flux == 'y')
					{
						double modelflux = 0; double cosb;
						if (strcoll(projectionType,"external") ==0)
						{
							//normalizing the total flux of the model to the same of the data
							if(smooth == 'y')
							{
								for(int i=0;i<RAsize;i++) for(int j=0;j<DECsize;j++) for(int k=0;k<VELsize;k++) if(!isnan(smoothedOutcube[i][j][k])) modelflux += smoothedOutcube[i][j][k];
							}
							else
							{
								for(int i=0;i<RAsize;i++) for(int j=0;j<DECsize;j++) for(int k=0;k<VELsize;k++) if(!isnan(outcube[i][j][k])) modelflux += outcube[i][j][k];
							}
						}
						else if(strcoll(projectionType,"internal") ==0)
						{
							//normalizing the total flux of the model to the same of the data
							if(smooth == 'y')
							{
								for(int i=0;i<RAsize;i++) for(int j=0;j<DECsize;j++) for(int k=0;k<VELsize;k++) 
								{
									if(!isnan(smoothedOutcube[i][j][k]))
									{
										cosb = cos(DEC[j]*DEGTORAD);
										modelflux += smoothedOutcube[i][j][k]*cosb;
									}
								}
							}
							else
							{
								for(int i=0;i<RAsize;i++) for(int j=0;j<DECsize;j++) for(int k=0;k<VELsize;k++) 
								{
									if(!isnan(outcube[i][j][k]))
									{
										cosb = cos(DEC[j]*DEGTORAD);
										modelflux += outcube[i][j][k]*cosb;
									}
								}
							}
						} 
						modelflux *= fabs(cdelt1*cdelt2*cdelt3);
						cout<<"Model flux = "<<modelflux<<endl;
						fluxratio  = dataflux/modelflux;
						cout<<"DATA vs MODEL HI flux ratio: "<<fluxratio<<endl;
						if(smooth == 'y')
						{
							for(int i=0;i<RAsize;i++) for(int j=0;j<DECsize;j++) for(int k=0;k<VELsize;k++) smoothedOutcube[i][j][k] *= fluxratio;
						}
						else
						{
							for(int i=0;i<RAsize;i++) for(int j=0;j<DECsize;j++) for(int k=0;k<VELsize;k++) outcube[i][j][k] *= fluxratio;
						}
						cout<<"MODEL HI flux normalized to the datacube"<<endl;
					}
					else {cout<<"Warning: DATA and MODEL HI flux may differ!"<<endl;}
					
					//evaluating residuals
					cout<<"Evaluating residuals..."<<endl;
					double cosb = 1;
					if(smooth=='y')
					{
						for(int j=0;j<DECsize;j++) 
						{
							if(strcmp(projectionType, "internal")==0) cosb = cos(DEC[j]*0.0174532925);
							for(int i=0;i<RAsize;i++) for(int k=0;k<VELsize;k++)
							{
								if((datacube(i,j,k)>threshold_res || smoothedOutcube[i][j][k]>threshold_res) && !isnan(datacube(i,j,k)*smoothedOutcube[i][j][k]))
								{
									// R1 = |dat-mod|
									res_1 += cosb*fabs(datacube(i,j,k) - smoothedOutcube[i][j][k]);
									
									// R2 = (dat-mod)^2  
									res_2 += cosb*powf(datacube(i,j,k) - smoothedOutcube[i][j][k],2);
									
									// R3 = |dat-mod|/max(dat,mod)
									res_3 += cosb*fabs(datacube(i,j,k) - smoothedOutcube[i][j][k])/max(datacube(i,j,k),double(smoothedOutcube[i][j][k]));
								}
							}
						}
					}
					else
					{
						for(int j=0;j<DECsize;j++) 
						{
							if(strcmp(projectionType, "internal")==0) cosb = cos(DEC[j]*0.0174532925);
							for(int i=0;i<RAsize;i++) for(int k=0;k<VELsize;k++)
							{
								if((datacube(i,j,k)>threshold_res || outcube[i][j][k]>threshold_res) && !isnan(datacube(i,j,k)*smoothedOutcube[i][j][k]))
								{
									// R1 = |dat-mod|
									res_1 += cosb*fabs(datacube(i,j,k) - outcube[i][j][k]);
									
									// R2 = (dat-mod)^2  
									res_2 += cosb*powf(datacube(i,j,k) - outcube[i][j][k],2);
									
									// R3 = |dat-mod|/max(dat,mod)
									res_3 += cosb*fabs(datacube(i,j,k) - outcube[i][j][k])/max(datacube(i,j,k),double(outcube[i][j][k]));
								}
							}
						}
					}
					if(isnan(res_1)) res_1 = 1e10;
					if(isnan(res_2)) res_2 = 1e10;
					if(isnan(res_3)) res_3 = 1e10;
				}
				//COMPARISON WITH AN ABSORPTION LIST
				
				else if (strcmp(compare_to, "list") == 0)
				{
					ofstream reproducibility("reproducibility.txt");
					string profileplot = "profiles.pdf";
					profiles_macro myprofiles(profileplot);
					std::vector<double> profile(VELsize);
					nabs_model = 0; 
					double vel_peak, vel_max, vel_med;
					for(int n=0;n<l_abs.size();n++)
					{
						int i = (l_abs[n]-crval1)/cdelt1+0.5;
						int j = (b_abs[n]-crval2)/cdelt2+0.5;
						if(smooth=='y') {for(int k=0;k<VELsize;k++) profile[k] = smoothedOutcube[i][j][k];}
						else {for(int k=0;k<VELsize;k++) profile[k] = outcube[i][j][k];}
						vel_peak = myutil::closest_peak(VEL,profile,v_abs[n],0.05);
						vel_max = myutil::max(VEL,profile);
						vel_med = myutil::median(VEL,profile); 
						
						if(!isnan(vel_peak)) 
						{
							cout<<"Absorption line "<<id_abs[n]<<" found in the model"<<endl;
							cout<<"		velocity observed:  "<<v_abs[n]<<" km/s"<<endl;
							cout<<"		velocity predicted: "<<vel_peak<<" km/s (peak)"<<endl;
							cout<<"		                    "<<vel_max<<" km/s (max)"<<endl;
							cout<<"		                    "<<vel_med<<" km/s (med)"<<endl;
							
							myprofiles.add_profile(id_abs[n],ref_abs[n],type_abs[n],VEL,profile,l_abs[n],b_abs[n],v_abs[n],vel_peak,vel_max,vel_med);
							
							nabs_model+=1.;
							
							// R1: closest peak
							res_1 += pow(vel_peak - v_abs[n],2);
						
							// R2: highest peak 
							res_2 += pow(vel_max - v_abs[n],2);
						
							// R3: median 
							res_3 += pow(vel_med - v_abs[n],2);
							
							reproducibility<<l_abs[n]<<"  "<<pow(vel_peak - v_abs[n],2)<<endl;
						}
					}
					reproducibility.close();
					eta = static_cast<double>(nabs_model)/static_cast<double>(l_abs.size());
					cout<<"Percentage of absorptions having a counterpart in the model: "<<eta*100<<" percent"<<endl;
					if(nabs_model>0)
					{
						res_1/= nabs_model*nabs_model;
						res_2/= nabs_model*nabs_model;
						res_3/= nabs_model*nabs_model;
					}
					else
					{
						res_1 = NAN;
						res_2 = NAN;
						res_3 = NAN;
					}
					//close the profile macro
					myprofiles.close();
				}
				cout<<"done"<<endl;
				
				if(strcoll(effects,"accrdrag")==0)
				{
					if(strcmp(compare_to, "cube")==0)
					{
						if(normalize_flux=='y') {myresidual<<mypar1<<setw(12)<<mypar2<<setw(12)<<mypar3<<setw(12)<<fluxratio*M_HI_halo<<setw(12)<<io.get_GlobalAR()*fluxratio<<setw(12)<<res_1<<setw(12)<<res_2<<setw(12)<<res_3<<endl;}
						else {myresidual<<mypar1<<setw(12)<<mypar2<<setw(12)<<mypar3<<setw(12)<<M_halo_above<<setw(12)<<io.get_GlobalAR()<<setw(12)<<res_1<<setw(12)<<res_2<<setw(12)<<res_3<<endl;}
					}
					else if(strcmp(compare_to, "list")==0)
					{
						myresidual<<mypar1<<setw(12)<<mypar2<<setw(12)<<mypar3<<setw(12)<<eta*100<<setw(12)<<io.get_GlobalAR()<<setw(12)<<res_1<<setw(12)<<res_2<<setw(12)<<res_3<<endl;
					}
					
				}
				//erasing outcube and smoothedOutcube for the next iteration	
				for(int i=0; i<RAsize; i++)
				{
					for(int j=0; j<DECsize; j++)
					{
						for (int k=0; k<VELsize; k++) {outcube[i][j][k] = 0.; smoothedOutcube[i][j][k] = 0.;}
					}
				}
			}
			/* END cycle n_iter2 */
			time(&end);
			double time_min = difftime(end, start)/60.;
			if(time_min>=1.) {cout<<"Computational time for this model: "<<time_min<<" minutes"<<endl;}
			else {cout<<"Computational time for this model: "<<time_min*60.<<" seconds"<<endl;}
			if(iter_mode=='y')
			{
				double time_min_remaining = time_min*(n_iter1*n_iter2*n_iter3 - global_iter);
				cout<<"APPROX. REMAINING COMPUTATIONAL TIME: ";
				if(time_min_remaining<60.) {cout<<time_min_remaining<<" minutes"<<endl;}
				else if(time_min_remaining<1440) {cout<<time_min_remaining/60.<<" hours"<<endl;} 
				else {cout<<time_min_remaining/1440.<<" days"<<endl;}
			}
				
		}
			n_start3 = 0;
		}
		n_start2 = 0;
	}
	
	
	if (iter_mode=='y' && fullResiduals=='y') 
	{
		cout<<endl;
		cout<<"BUILDING THE RESIDUAL CUBES"<<endl;
		
		//temporary vectors
		const int res_dim = n_iter1*n_iter2*n_iter3;
		std::vector<double> temp_abs(res_dim), temp_sqr(res_dim), temp_rel(res_dim), temp_tot(res_dim);
		
		//read the residual file
		ifstream data(resfile);
		string temp;
		for(int i=0;i<res_dim;i++)
		{
			data>>temp>>temp>>temp>>temp>>temp>>temp_abs[i]>>temp_sqr[i]>>temp_rel[i];
		}
		data.close();
		//delete the 0s 
		for(int i=0;i<res_dim;i++) if(temp_abs[i]<=0) {temp_abs[i]=NAN; temp_sqr[i]=NAN; temp_rel[i]=NAN;}
		
		//normalize residuals to the minimum
		double min_abs = myutil::min(temp_abs);
		double min_sqr = myutil::min(temp_sqr);
		double min_rel = myutil::min(temp_rel);
		for(int i=0;i<res_dim;i++) 
		{
			temp_abs[i]/=min_abs; temp_sqr[i]/=min_sqr; temp_rel[i]/=min_rel; 
		}
		
		//total residual
		for(int i=0;i<res_dim;i++) temp_tot[i] = temp_abs[i] + temp_sqr[i] + temp_rel[i];
		double min_tot = myutil::min(temp_tot);
		for(int i=0;i<res_dim;i++) temp_tot[i] /= min_tot;
		
		//initialize cubes for residuals
		double *res_abs_ptr = new double [res_dim]; 
		double *res_sqr_ptr = new double [res_dim];
		double *res_rel_ptr = new double [res_dim];
		double *res_tot_ptr = new double [res_dim];
		myarray::double3d res_abs, res_sqr, res_rel, res_tot;
		res_abs.init(res_abs_ptr, n_iter1, n_iter2, n_iter3);
		res_sqr.init(res_sqr_ptr, n_iter1, n_iter2, n_iter3);
		res_rel.init(res_rel_ptr, n_iter1, n_iter2, n_iter3);
		res_tot.init(res_tot_ptr, n_iter1, n_iter2, n_iter3);
		for(int i=0;i<n_iter1;i++) for(int j=0;j<n_iter2;j++) for(int k=0; k<n_iter3; k++) 
		{
			res_abs(i,j,k) = temp_abs[i*n_iter2*n_iter3 + j*n_iter3 + k];
			res_sqr(i,j,k) = temp_sqr[i*n_iter2*n_iter3 + j*n_iter3 + k];
			res_rel(i,j,k) = temp_rel[i*n_iter2*n_iter3 + j*n_iter3 + k];
			res_tot(i,j,k) = temp_tot[i*n_iter2*n_iter3 + j*n_iter3 + k];
		}
		 
		//cube axis
		std::vector<double> axis1, axis2, axis3;
		for(int i=0;i<n_iter1;i++) axis1.push_back(param1_0 + i*delta_param1);
		for(int i=0;i<n_iter2;i++) axis2.push_back(param2_0 + i*delta_param2);
		for(int i=0;i<n_iter3;i++) axis3.push_back(param3_0 + i*delta_param3);
		
		char smooth_res; findValue(inputfile, smooth_res, "smooth_res");
		if(smooth_res=='y')
		{
			cout<<"Smoothing the residuals to the desired resolution: "<<endl;
			
			cout<<"Old beam : ("<<fabs(delta_param1)<<","<<fabs(delta_param2)<<","<<fabs(delta_param3)<<")"<<endl;   
			double beam_par1; findValue(inputfile, beam_par1, "beam_par1");
			double beam_par2; findValue(inputfile, beam_par2, "beam_par2");
			double beam_par3; findValue(inputfile, beam_par3, "beam_par3");
			cout<<"New beam : ("<<beam_par1<<","<<beam_par2<<","<<beam_par3<<")"<<endl;
			fitsutil::smooth3D(res_abs,axis1,axis2,axis3,beam_par1,beam_par2,beam_par3,4,false);
			fitsutil::smooth3D(res_sqr,axis1,axis2,axis3,beam_par1,beam_par2,beam_par3,4,false);
			fitsutil::smooth3D(res_rel,axis1,axis2,axis3,beam_par1,beam_par2,beam_par3,4,false);
			fitsutil::smooth3D(res_tot,axis1,axis2,axis3,beam_par1,beam_par2,beam_par3,4,false);
			min_abs = myutil::min(res_abs); 
			min_sqr = myutil::min(res_sqr); 
			min_rel = myutil::min(res_rel); 
			min_tot = myutil::min(res_tot);
			for(int i=0;i<n_iter1;i++) for(int j=0;j<n_iter2;j++) for(int k=0;k<n_iter3;k++)
			{
				res_abs(i,j,k)/=min_abs;
				res_sqr(i,j,k)/=min_sqr;
				res_rel(i,j,k)/=min_rel;
				res_tot(i,j,k)/=min_tot;
			}
		}
		
		//write cubes
		if(strcmp(compare_to, "cube")==0)
		{
			fitsutil::print_cube("residuals_abs.fits",axis1,axis2,axis3,res_abs_ptr,par1,par2,par3,"units","units","units","units");
			fitsutil::print_cube("residuals_sqr.fits",axis1,axis2,axis3,res_sqr_ptr,par1,par2,par3,"units","units","units","units");
			fitsutil::print_cube("residuals_rel.fits",axis1,axis2,axis3,res_rel_ptr,par1,par2,par3,"units","units","units","units");
			fitsutil::print_cube("residuals_tot.fits",axis1,axis2,axis3,res_tot_ptr,par1,par2,par3,"units","units","units","units");
		}
		else if(strcmp(compare_to, "list")==0)
		{
			fitsutil::print_cube("residuals_peak.fits",axis1,axis2,axis3,res_abs_ptr,par1,par2,par3,"units","units","units","units");
			fitsutil::print_cube("residuals_max.fits",axis1,axis2,axis3,res_sqr_ptr,par1,par2,par3,"units","units","units","units");
			fitsutil::print_cube("residuals_med.fits",axis1,axis2,axis3,res_rel_ptr,par1,par2,par3,"units","units","units","units");
			fitsutil::print_cube("residuals_tot.fits",axis1,axis2,axis3,res_tot_ptr,par1,par2,par3,"units","units","units","units");
		}
	}
	
	return 0;
}
