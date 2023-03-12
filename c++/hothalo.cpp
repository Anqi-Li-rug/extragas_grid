#include <cmath>
#include "hothalo.h"
#include "extragasPlots.h"
#include <stdio.h>
extern const double CONV1, PI, RAD;
const double CONV3=4.04742e-8; // Mo/kpc^3 -> H_atoms/cm^-3
const double K_B=1.38066e-16; // Boltzmann constant
const double M_P=1.672649e-24; // proton mass (g)
const double KEV=1.1604e7; // kev -> K
const double KPC=3.08567802e21; // kpc (cm)

using namespace ffarray;

HotHalo::HotHalo(float v0_hot, 
		 float T_hot, 
		 float mu_hot,
		 float Rmax_hot, 
		 float zmin_hot, 
		 float zmax_hot, 
		 double lum_hot)
{
  v0=v0_hot;
  T=T_hot;
  mu=mu_hot;
  Rmax=Rmax_hot;
  zmin=zmin_hot;
  zmax=zmax_hot;
  lum=lum_hot;
  /* CALCULATE HALO EMISSIVITY AND LUMINOSITY */
  conv_X=-2.1106e-23
    +3.82959e-22*(T/KEV)
    -1.30922e-21*(T/KEV)*(T/KEV)
    +2.11345e-21*(T/KEV)*(T/KEV)*(T/KEV)
    -1.28213e-21*(T/KEV)*(T/KEV)*(T/KEV)*(T/KEV); 
  // only valid for 0.08 < T < 0.6 keV = 9.3e5 < T < 7.0e6 K
  /* 
     polinomial fit to xspec results (MEKAL model)
     flux in the range 0.3-2 keV 
     units is cm^6 
     (calculated from xspec with norm=1e-14 and np=1 cm^-3 assumed
     already multiplied by 4pi*d^2)
  */
  // conv_X=-8.3395e-24+1.67826e-22*(T_hot/keV)
  // -2.44989e-22*(T_hot/keV)*(T_hot/keV)
  // +1.05846e-22*(T_hot/keV)*(T_hot/keV)*(T_hot/keV); // for .08<T<1 keV
}

HotHalo::~HotHalo()
{}

void HotHalo::SetDensity(Array2D<double>& rho_hot_ik1,
			 const Array1D<double> R_i1, 
			 const Array1D<double> z_k1,
			 const Array2D<double>& pot_ik1,
			 const Array2D<double>& F_R1,
			 Array2D<double>& rho_hot_ik2,
			 const Array1D<double> R_i2, 
			 const Array1D<double> z_k2,
			 const Array2D<double>& pot_ik2,
			 const Array2D<double>& F_R2)
{
  /* 
     CONSTRUCTION OF THE HOT HALO
  */
  FILE* f_hot_halo_R=NULL;
  FILE* f_hot_halo_z=NULL;
  int hi, hk;
  float mass_halo, mass_halo2;
  int nR1=R_i1.dim();
  int nz1=z_k1.dim();
  int nR2=R_i2.dim();
  int nz2=z_k2.dim();
  float stepR1=R_i1[1]-R_i1[0];
  float stepz1=z_k1[1]-z_k1[0];
  Array2D<float> emissivity_ik(nR1,nz1), brightness_ik;
 
  // writing outer matrix of density (not normalized)
  for (hi=0; hi<nR1; hi++) {
    for (hk=0; hk<nz1; hk++) {
      rho_hot_ik1[hi][hk]=exp(-(mu*M_P/K_B/T)*
			      ((pot_ik1[hi][hk])/CONV1/CONV1-
			       v0*v0/2.)*1.e10); // kpc^-3
    }
  }
  
  // writing inner matrix of density (not normalized)
  for (hi=0; hi<nR2; hi++) {
    for (hk=0; hk<nz2; hk++) {
      rho_hot_ik2[hi][hk]=exp(-(mu*M_P/K_B/T)*
			      (pot_ik2[hi][hk]/CONV1/CONV1-
			       v0*v0/2.)*1.e10); // kpc^-3
    }
  }
 
  norm_lum=0.;
  for (hi=0; hi<nR1; hi++) {
    for (hk=0; hk<nz1; hk++) {
      //      cout << hi << " " << hk << " " << nR1 << " " << nz1 << endl;
      emissivity_ik[hi][hk]=
				conv_X*rho_hot_ik1[hi][hk]*rho_hot_ik1[hi][hk]*CONV3*CONV3; 
      // erg s^-1 cm^-3
      if (R_i1[hi] <= Rmax && z_k1[hk] 
	  >= zmin && z_k1[hk] <= zmax) {
	norm_lum+=
	  4*PI*R_i1[hi]*stepR1*stepz1*emissivity_ik[hi][hk]*KPC*KPC*KPC;
      }
    }
  }
  
  f_hot_halo_R=fopen("extragas_halo_R.dat", "w");
  f_hot_halo_z=fopen("extragas_halo_z.dat", "w"); 
  
  for (hi=0; hi<nR2; hi++) {
    for (hk=0; hk<nz2; hk++) {
      /* 
	 NORMALIZATION inner matrix 
      */
      rho_hot_ik2[hi][hk]=
				rho_hot_ik2[hi][hk]*sqrt(lum/norm_lum); //Mo/kpc^3
      //	  vphi_hot_ik2[hi][hk]=v0_hot;
      // Lz_hot_ik2[hi][hk]=2.*pi*R_i2[hi]*R_i2[hi]*stepR2*stepz2*rho_hot_ik2[hi][hk]*v0_hot*conv1;
      if (hi == 0) {
				fprintf(f_hot_halo_z,"%f %.3e %.3e \n",z_k2[hk], rho_hot_ik2[hi][hk], rho_hot_ik2[hi][hk]*CONV3);
      }
    }
    // write files
    if (hi == 0) {
      fprintf(f_hot_halo_R,"%f %.3e %.3e \n",z_k2[0], rho_hot_ik2[0][0], rho_hot_ik2[0][0]*CONV3);
    } 
    else {
      fprintf(f_hot_halo_R,"%f %.3e %.3e \n",R_i2[hi], rho_hot_ik2[hi][0], rho_hot_ik2[hi][0]*CONV3);
    }
  }
  for (hi=0; hi<nR1; hi++) {
    for (hk=0; hk<nz1; hk++) {
      /* 
	 NORMALIZATION outer matrix 
      */
      rho_hot_ik1[hi][hk]=
				rho_hot_ik1[hi][hk]*sqrt(lum/norm_lum); //Mo/kpc^3
      emissivity_ik[hi][hk]=emissivity_ik[hi][hk]*(lum/norm_lum);
      //	  vphi_hot_ik1[hi][hk]=v0_hot;
      //	  vphi_hot_ik[hi][hk]=v0_hot;
      //	  Lz_hot_ik1[hi][hk]=2.*pi*R_i1[hi]*R_i1[hi]*stepR1*stepz1*rho_hot_ik1[hi][hk]*v0_hot*conv1;
      mass_halo+=4*PI*R_i1[hi]*rho_hot_ik1[hi][hk]*stepR1*stepz1;
      if (R_i1[hi] <= Rmax && z_k1[hk]
	  >= zmin && z_k1[hk] <= zmax) {
	mass_halo2+=4*PI*R_i1[hi]*rho_hot_ik1[hi][hk]*stepR1*stepz1;
      }
      if (hi == 0) {
	if (z_k1[hk] > z_k2[nz2-1]){
	  fprintf(f_hot_halo_z,"%f %.3e %.3e \n",z_k1[hk], rho_hot_ik1[hi][hk], rho_hot_ik1[hi][hk]*CONV3);
	}
      }
      // PB values along R_i[0] too LOW!!!!!!!!!!
    }
    if (R_i1[hi] > R_i2[nR2-1]){
      fprintf(f_hot_halo_R,"%f %.3e %.3e \n",R_i1[hi], rho_hot_ik1[hi][0], rho_hot_ik1[hi][0]*CONV3);
      //	rmass_hot_halo+=4*pi*R_i1[hi]*R_i1[hi]*rho_hot_ik1[hi][0]*(R_i1[1]-R_i1[0]);
    }
  }
  //  print(rho_hot_ik1);
  
  fclose(f_hot_halo_R);
  fclose(f_hot_halo_z);
  plot8();
  
  printf("Central density of the hot halo = %.9e Mo/kpc^3\n",rho_hot_ik2[0][0]/CONV1/CONV1);
  printf("Total mass of the hot halo = %.2e Mo\n",mass_halo);
  printf("Mass between z>%.1f kpc and z<%.1f kpc and within R<%.1f kpc = %.2e Mo\n",zmin,zmax,Rmax,mass_halo2);
  
}

void HotHalo::SetDensity(Array2D<double>& rho_hot_ik1,
			 const Array1D<double> R_i1, 
			 const Array1D<double> z_k1,
			 const Array2D<double>& pot_ik1,
			 const Array2D<double>& F_R1)
{
  /* 
     CONSTRUCTION OF THE HOT HALO
  */
  FILE* f_hot_halo_R=NULL;
  FILE* f_hot_halo_z=NULL;
  int hi, hk;
  float mass_halo, mass_halo2;
  int nR1=R_i1.dim();
  int nz1=z_k1.dim();
  float stepR1=R_i1[1]-R_i1[0];
  float stepz1=z_k1[1]-z_k1[0];
  Array2D<float> emissivity_ik(nR1,nz1), brightness_ik;
  
  // writing outer matrix of density (not normalized)
  for (hi=0; hi<nR1; hi++) {
    for (hk=0; hk<nz1; hk++) {
      rho_hot_ik1[hi][hk]=exp(-(mu*M_P/K_B/T)*
			      ((pot_ik1[hi][hk])/CONV1/CONV1-
			       v0*v0/2.)*1.e10); // kpc^-3
    }
  }
  
  norm_lum=0.;
  for (hi=0; hi<nR1; hi++) {
    for (hk=0; hk<nz1; hk++) {
      emissivity_ik[hi][hk]=
	conv_X*rho_hot_ik1[hi][hk]*rho_hot_ik1[hi][hk]*CONV3*CONV3; 
      // erg s^-1 cm^-3
      if (R_i1[hi] <= Rmax && z_k1[hk] 
	  >= zmin && z_k1[hk] <= zmax) {
	norm_lum+=
	  4*PI*R_i1[hi]*stepR1*stepz1*emissivity_ik[hi][hk]*KPC*KPC*KPC;
      }
    }
  }
  
  for (hi=0; hi<nR1; hi++) {
    for (hk=0; hk<nz1; hk++) {
      /* 
	 NORMALIZATION outer matrix 
      */
      rho_hot_ik1[hi][hk]=
	rho_hot_ik1[hi][hk]*sqrt(lum/norm_lum); //Mo/kpc^3
      emissivity_ik[hi][hk]=emissivity_ik[hi][hk]*(lum/norm_lum);
      //	  vphi_hot_ik1[hi][hk]=v0_hot;
      //	  vphi_hot_ik[hi][hk]=v0_hot;
      //	  Lz_hot_ik1[hi][hk]=2.*pi*R_i1[hi]*R_i1[hi]*stepR1*stepz1*rho_hot_ik1[hi][hk]*v0_hot*conv1;
      mass_halo+=4*PI*R_i1[hi]*rho_hot_ik1[hi][hk]*stepR1*stepz1;
      if (R_i1[hi] <= Rmax && z_k1[hk]
	  >= zmin && z_k1[hk] <= zmax) {
	mass_halo2+=4*PI*R_i1[hi]*rho_hot_ik1[hi][hk]*stepR1*stepz1;
      }
      if (hi == 0) {
	//	if (z_k1[hk] > z_k2[nz2-1]){
	  fprintf(f_hot_halo_z,"%f %.3e %.3e \n",z_k1[hk], rho_hot_ik1[hi][hk], rho_hot_ik1[hi][hk]*CONV3);
	  //	}
      }
      // PB values along R_i[0] too LOW!!!!!!!!!!
    }
    //    if (R_i1[hi] > R_i2[nR2-1]){
      fprintf(f_hot_halo_R,"%f %.3e %.3e \n",R_i1[hi], rho_hot_ik1[hi][0], rho_hot_ik1[hi][0]*CONV3);
      //	rmass_hot_halo+=4*pi*R_i1[hi]*R_i1[hi]*rho_hot_ik1[hi][0]*(R_i1[1]-R_i1[0]);
      //    }
  }

  fclose(f_hot_halo_R);
  fclose(f_hot_halo_z);
  plot8();
  
  printf("Central density of the hot halo = %.9e Mo/kpc^3\n",rho_hot_ik1[0][0]/CONV1/CONV1);
  printf("Total mass of the hot halo = %.2e Mo\n",mass_halo);
  printf("Mass between z>%.1f kpc and z<%.1f kpc and within R<%.1f kpc = %.2e Mo\n",zmin,zmax,Rmax,mass_halo2);
  
}



/*
  FITS IMAGES OF THE HOT HALO

// HALO DENSITY 
Array2D<float> hot_halo(RAsize,DECsize); 
float ttr, ttz, tv;
Array2D<float> norm_hot_halo1(RAsize,DECsize);
Array2D<float> norm_hot_halo2(RAsize,DECsize);
projections halo(ttr,90.,ttz,0.,0.,0.); // output: xi,yi,vrad;
for (int ti=0; ti<RAsize; ti++) {
  for (int tk=0; tk<DECsize; tk++) {
    hot_halo[ti][tk]=0;
    norm_hot_halo1[ti][tk]=0;
    norm_hot_halo2[ti][tk]=0;
  }
}
for (hi=0; hi<nR1; hi++) {
  for (hk=0; hk<nz1; hk++) {
    ttr=R_i1[hi];
    ttz=z_k1[hk];
    halo.external_view(incl,xi,yi,vrad); // output: xi,yi,vrad;
    matrixcell_ext(xi, yi, vrad,xc, yc, vc);
    if (xc < RAsize && xc >= 0 && yc < DECsize && yc >=0) {
      tv=fabs(rho_hot_ik1[hi][hk]);
      hot_halo[xc][yc]+=tv;
      hot_halo[RAsize-xc-1][yc]+=tv;
      hot_halo[xc][DECsize-yc-1]+=tv;
      hot_halo[RAsize-xc-1][DECsize-yc-1]+=tv;
      norm_hot_halo1[xc][yc]+=1;
      norm_hot_halo1[RAsize-xc-1][yc]+=1;
      norm_hot_halo1[xc][DECsize-yc-1]+=1;
      norm_hot_halo1[RAsize-xc-1][DECsize-yc-1]+=1;
    }
  }
}
for (xc=0; xc<RAsize; xc++) {
  for (yc=0; yc<DECsize; yc++) {
    if (norm_hot_halo1[xc][yc] != 0){
      hot_halo[xc][yc]=hot_halo[xc][yc]/norm_hot_halo1[xc][yc];
    }
  }
}
for (hi=0; hi<nR2; hi++) {
  for (hk=0; hk<nz2; hk++) {
    ttr=R_i2[hi];
    ttz=z_k2[hk];
    halo.external_view(incl,xi,yi,vrad); // output: xi,yi,vrad;
    matrixcell_ext(xi, yi, vrad,xc, yc, vc);
    if (xc < RAsize && xc >= 0 && yc < DECsize && yc >=0) {
      tv=fabs(rho_hot_ik2[hi][hk]);
      if (norm_hot_halo2[xc][yc] == 0) {
	hot_halo[xc][yc]=tv;
	hot_halo[RAsize-xc-1][yc]=tv;
	hot_halo[xc][DECsize-yc-1]=tv;
	hot_halo[RAsize-xc-1][DECsize-yc-1]=tv;
      } else {
	hot_halo[xc][yc]+=tv;
	hot_halo[RAsize-xc-1][yc]+=tv;
		hot_halo[xc][DECsize-yc-1]+=tv;
		hot_halo[RAsize-xc-1][DECsize-yc-1]+=tv;
	      }
	      norm_hot_halo2[xc][yc]+=1;
	      norm_hot_halo2[RAsize-xc-1][yc]+=1;
	      norm_hot_halo2[xc][DECsize-yc-1]+=1;
	      norm_hot_halo2[RAsize-xc-1][DECsize-yc-1]+=1;
	    }
	  }
	}
	for (xc=0; xc<RAsize; xc++) {
	  for (yc=0; yc<DECsize; yc++) {
	    if (norm_hot_halo2[xc][yc] != 0){
	      hot_halo[xc][yc]=hot_halo[xc][yc]/norm_hot_halo2[xc][yc];
	    }
	  }
	}
	printHeaderTot(RAsize, DECsize, cdelt1, cdelt2, proj_type);
	writefits_2D("hhalo_density.fits", &hot_halo[0][0], RAsize, DECsize, "header.tab");  
	
	// HALO SURFACE BRIGHTNESS 
	float subR[2], subz[2], subemissivity_ik[2][2];
	for (hi=0; hi<nR1; hi++) {
	  for (hk=0; hk<nz1; hk++) {
	    brightness_ik[hi][hk]=emissivity_ik[hi][hk]*KPC;
	    for (int li=1; li<nR1; li++) {
	      if ((int)sqrt(hi*hi+li*li) < nR1){
		brightness_ik[hi][hk]+=2*emissivity_ik[(int)sqrt(hi*hi+li*li)][hk]*KPC;
	      }
	    }
	  }
	}
	
	for (int ti=0; ti<RAsize; ti++) {
	  for (int tk=0; tk<DECsize; tk++) {
	    hot_halo[ti][tk]=0;
	    norm_hot_halo1[ti][tk]=0;
	  }
	}
	for (hi=0; hi<nR1; hi++) {
	  for (hk=0; hk<nz1; hk++) {
	    ttr=R_i1[hi];
	    ttz=z_k1[hk];
	    halo.external_view(incl,xi,yi,vrad); // output: xi,yi,vrad;
	    matrixcell_ext(xi, yi, vrad,xc, yc, vc);
	    if (xc < RAsize && xc >= 0 && yc < DECsize && yc >=0) {
	      tv=fabs(brightness_ik[hi][hk]);
	      hot_halo[xc][yc]+=tv;
	      hot_halo[RAsize-xc-1][yc]+=tv;
	      hot_halo[xc][DECsize-yc-1]+=tv;
	      hot_halo[RAsize-xc-1][DECsize-yc-1]+=tv;
	      norm_hot_halo1[xc][yc]+=1;
	      norm_hot_halo1[RAsize-xc-1][yc]+=1;
	      norm_hot_halo1[xc][DECsize-yc-1]+=1;
	      norm_hot_halo1[RAsize-xc-1][DECsize-yc-1]+=1;
	    }
	  }
	}
	for (xc=0; xc<RAsize; xc++) {
	  for (yc=0; yc<DECsize; yc++) {
	    if (norm_hot_halo1[xc][yc] != 0){
	      hot_halo[xc][yc]=hot_halo[xc][yc]/norm_hot_halo1[xc][yc];
	    }
	  }
	}
	printHeaderTot(RAsize, DECsize, cdelt1, cdelt2, proj_type);
	writefits_2D("hhalo_brightness.fits", &hot_halo[0][0], RAsize, DECsize, "header.tab");  
*/      
