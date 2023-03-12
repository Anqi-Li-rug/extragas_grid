/*
  Utilities
  
  This file contains general utilities to be used by exgas
*/

#include <iostream>
#include <cstring>
#include <cmath>
#include <stdio.h>
#include "extragasFits.h"

/* required by ever y program that uses CFITSIO  */
#include "fitsio.h"  
//#include "/scisoft/i386/Packages/cfitsio/include/fitsio.h"

using namespace std;

const double pig=3.14159265359;


void printHeader(int RAsize, int DECsize, int VELsize,
		 float cdelt1, float cdelt2, float cdelt3,
		 float crval1, float crval2, float crval3,
		 char proj_type[], char crpix[], char Jy2K)
{
  FILE* fitsheader=NULL;
  fitsheader=fopen("header.tab", "w");

  if (strcoll(proj_type,"hvcs") ==0)
    crval1+=180;
  /*
    if (strcoll(proj_type,"external") ==0)
    crval2=atof(value);//+cdelt2/scale;//is needed or not(???) 
  if (strcoll(proj_type,"internal") ==0 || strcoll(proj_type,"hvcs") ==0)
    crval2=atof(value);
  */
  
  fprintf(fitsheader,"CDELT1 = %f \n",cdelt1);
  fprintf(fitsheader,"CDELT2 = %f \n",cdelt2);
  fprintf(fitsheader,"CDELT3 = %f \n",cdelt3);
  fprintf(fitsheader,"CUNIT1 = DEGREE \n");
  fprintf(fitsheader,"CUNIT2 = DEGREE \n");
  fprintf(fitsheader,"CUNIT3 = KM/S \n");
  if (strcoll(proj_type,"external") ==0) {
    fprintf(fitsheader,"CTYPE1 = RA---NCP \n");
    fprintf(fitsheader,"CTYPE2 = DEC--NCP \n");
  }
  if (strcoll(proj_type,"internal") ==0 || strcoll(proj_type,"hvcs") ==0) {
    fprintf(fitsheader,"CTYPE1 = RA \n");
    fprintf(fitsheader,"CTYPE2 = DEC \n");
  }
  fprintf(fitsheader,"CTYPE3 = VELO-HEL \n");
  int crpix1, crpix2;
  if (strcoll(crpix,"first")==0) {
    crpix1=1;
    crpix2=1;
  } 
  else {
    if (strcoll(crpix,"middle")==0) {
      crpix1=RAsize/2;
      crpix2=DECsize/2;
    }
    else {
      crpix1=1;
      crpix2=1;
    }
  }
  fprintf(fitsheader,"CRPIX1 = %i \n",crpix1);
  fprintf(fitsheader,"CRPIX2 = %i \n",crpix2);
  int crpix3=1;
  fprintf(fitsheader,"CRPIX3 = %i \n",crpix3);
  fprintf(fitsheader,"CRVAL1 = %f \n",crval1);
  fprintf(fitsheader,"CRVAL2 = %f \n",crval2);
  fprintf(fitsheader,"CRVAL3 = %f \n",crval3);
  char bunit[10];
  if (Jy2K == 'y') 
    strcpy(bunit,"K");
  else
    strcpy(bunit,"JY/BEAM");
  fprintf(fitsheader,"BUNIT  = %s \n",bunit);
  char btype[]="intensity";
  fprintf(fitsheader,"BTYPE  = %s \n",btype);
  fclose(fitsheader);
}

void printHeaderTot(int RAsize, int DECsize, 
		    float cdelt1, float cdelt2,
		    float crval1, float crval2,
		    char proj_type[], char Jy2K
		    )
{
  FILE* fitsheader_tot=NULL;
  fitsheader_tot=fopen("header.tab", "w");

  fprintf(fitsheader_tot,"CDELT1 = %f \n",cdelt1);
  fprintf(fitsheader_tot,"CDELT2 = %f \n",cdelt2);
  fprintf(fitsheader_tot,"CUNIT1 = DEGREE \n");
  fprintf(fitsheader_tot,"CUNIT2 = DEGREE \n");
  if (strcoll(proj_type,"external") ==0) {
    fprintf(fitsheader_tot,"CTYPE1 = RA---NCP \n");
    fprintf(fitsheader_tot,"CTYPE2 = DEC--NCP \n");
  }
  if (strcoll(proj_type,"internal") ==0 || strcoll(proj_type,"hvcs") ==0) {
    fprintf(fitsheader_tot,"CTYPE1 = RA---NCP \n");
    fprintf(fitsheader_tot,"CTYPE2 = DEC--NCP \n");
  }
  int crpix1=RAsize/2;
  fprintf(fitsheader_tot,"CRPIX1 = %i \n",crpix1);
  int crpix2=DECsize/2;
  fprintf(fitsheader_tot,"CRPIX2 = %i \n",crpix2);
  fprintf(fitsheader_tot,"CRVAL1 = %f \n",crval1);
  fprintf(fitsheader_tot,"CRVAL2 = %f \n",crval2);
  if (Jy2K == 'y') 
    fprintf(fitsheader_tot,"BUNIT  = K*km/s \n");
  else
    fprintf(fitsheader_tot,"BUNIT  = JY/BEAM*km/s \n");
  char btype[]="intensity";
  fprintf(fitsheader_tot,"BTYPE  = %s \n",btype);
  fclose(fitsheader_tot);
}

void printHeaderDeltaV(float cdelt1, float cdelt2, float kpctograd)
{
  FILE* fitsheader_tot=NULL;
  fitsheader_tot=fopen("header.tab", "w");

  fprintf(fitsheader_tot,"CDELT1 = %f \n",-cdelt1/kpctograd);
  fprintf(fitsheader_tot,"CDELT2 = %f \n",cdelt2/kpctograd);
  fprintf(fitsheader_tot,"CUNIT1 = KPC \n");
  fprintf(fitsheader_tot,"CUNIT2 = KPC \n");
  fprintf(fitsheader_tot,"CTYPE1 = R(kpc) \n");
  fprintf(fitsheader_tot,"CTYPE2 = z(kpc) \n");
  fprintf(fitsheader_tot,"CRPIX1 = %i \n",3);
  fprintf(fitsheader_tot,"CRPIX2 = %i \n",3);
  fprintf(fitsheader_tot,"CRVAL1 = %f \n",0);
  fprintf(fitsheader_tot,"CRVAL2 = %f \n",0);
  char bunit[]="km/s/Myr";
  fprintf(fitsheader_tot,"BUNIT  = %s \n",bunit);
  char btype[]="intensity";
  fprintf(fitsheader_tot,"BTYPE  = %s \n",btype);
  fclose(fitsheader_tot);
}

void printHeaderDeltaL(float cdelt1, float cdelt2, float kpctograd)
{
  FILE* fitsheader_tot=NULL;
  fitsheader_tot=fopen("header.tab", "w");
  fprintf(fitsheader_tot,"CDELT1 = %f \n",-cdelt1/kpctograd);
  fprintf(fitsheader_tot,"CDELT2 = %f \n",cdelt2/kpctograd);
  fprintf(fitsheader_tot,"CUNIT1 = KPC \n");
  fprintf(fitsheader_tot,"CUNIT2 = KPC \n");
  fprintf(fitsheader_tot,"CTYPE1 = R(kpc) \n");
  fprintf(fitsheader_tot,"CTYPE2 = z(kpc) \n");
  fprintf(fitsheader_tot,"CRPIX1 = %i \n",3);
  fprintf(fitsheader_tot,"CRPIX2 = %i \n",3);
  fprintf(fitsheader_tot,"CRVAL1 = %f \n",0);
  fprintf(fitsheader_tot,"CRVAL2 = %f \n",0);
  char bunit[]="Mokpc^2/Myr^2";
  fprintf(fitsheader_tot,"BUNIT  = %s \n",bunit);
  char btype[]="intensity";
  fprintf(fitsheader_tot,"BTYPE  = %s \n",btype);
  fclose(fitsheader_tot);
}

void printHeaderIter(char par1[], char par2[], 
		     float param1_0, float param2_0,
		     float delta_param1, float delta_param2)
{
  FILE* fitsheader_iter=NULL;
  fitsheader_iter=fopen("header.tab", "w+");

  float icdelt1, icdelt2, crval1, crval2;

  if (param1_0 > 1e4) {
    crval1=param1_0/pow(10.,(int)log10(param1_0));
    icdelt1=delta_param1/pow(10.,(int)log10(param1_0));
  } else {
    icdelt1=delta_param1;
    crval1=param1_0;
  }
  if (param2_0 > 1e4) {
    crval2=param2_0/pow(10.,(int)log10(param2_0));
    icdelt2=delta_param2/pow(10.,(int)log10(param2_0));
  } else {
    icdelt2=delta_param2;
    crval2=param2_0;
  }
  fprintf(fitsheader_iter,"CDELT1 = %f \n",icdelt1);
  fprintf(fitsheader_iter,"CDELT2 = %f \n",icdelt2);
  fprintf(fitsheader_iter,"CTYPE1 = %s \n", par1);
  fprintf(fitsheader_iter,"CTYPE2 = %s \n", par2);
  fprintf(fitsheader_iter,"CRPIX1 = 0 \n");
  fprintf(fitsheader_iter,"CRPIX2 = 0 \n");
  fprintf(fitsheader_iter,"CRVAL1 = %f \n",crval1);
  fprintf(fitsheader_iter,"CRVAL2 = %f \n",crval2);
  fclose(fitsheader_iter);
  
}

void printHeaderMax(int RAsize, int VELsize, 
		    float cdelt1, float cdelt3,
		    float crval1, float crval2
		    )
{
  FILE* fitsheader=NULL;
  fitsheader=fopen("header.tab", "w");


  fprintf(fitsheader,"CDELT1 = %f \n",cdelt1);
  fprintf(fitsheader,"CDELT2 = %f \n",-cdelt3);
  fprintf(fitsheader,"CUNIT1 = DEGREE \n");
  fprintf(fitsheader,"CUNIT2 = KM/S \n");
  fprintf(fitsheader,"CTYPE1 = RA---NCP \n");
  fprintf(fitsheader,"CTYPE2 = VELO-HEL \n");
  int crpix1=RAsize/2;
  fprintf(fitsheader,"CRPIX1 = %i \n",crpix1);
  int crpix2=1;
  fprintf(fitsheader,"CRPIX2 = %i \n",crpix2);
  fprintf(fitsheader,"CRVAL1 = %f \n",crval1);
  fprintf(fitsheader,"CRVAL2 = %f \n",crval2);
  char bunit[]="JY/BEAM";
  fprintf(fitsheader,"BUNIT  = %s \n",bunit);
  char btype[]="intensity";
  fprintf(fitsheader,"BTYPE  = %s \n",btype);
  fclose(fitsheader);
}

void printHeaderAccr(int RAsize, int VELsize, 
		     float cdelt1, float cdelt2
		     )
{
  FILE* fitsheader=NULL;
  fitsheader=fopen("header.tab", "w");

  fprintf(fitsheader,"CDELT1 = %f \n",cdelt1);
  fprintf(fitsheader,"CDELT2 = %f \n",cdelt2);
  fprintf(fitsheader,"CUNIT1 = KPC \n");
  fprintf(fitsheader,"CUNIT2 = KPC \n");
  fprintf(fitsheader,"CTYPE1 = R(kpc) \n");
  fprintf(fitsheader,"CTYPE2 = z(kpc) \n");
  fprintf(fitsheader,"CRPIX1 = 0 \n");
  fprintf(fitsheader,"CRPIX2 = 0 \n");
  fprintf(fitsheader,"CRVAL1 = 0 \n");
  fprintf(fitsheader,"CRVAL2 = 0 \n");
  char bunit[]="JY/BEAM";
  fprintf(fitsheader,"BUNIT  = %s \n",bunit);
  char btype[]="intensity";
  fprintf(fitsheader,"BTYPE  = %s \n",btype);
  fclose(fitsheader);
}

void printHeaderSimple(int RAsize, int VELsize, 
		       float cdelt1, float cdelt2
		       )
{
  FILE* fitsheader=NULL;
  fitsheader=fopen("header.tab", "w");

  fprintf(fitsheader,"CDELT1 = %f \n",cdelt1);
  fprintf(fitsheader,"CDELT2 = %f \n",cdelt2);
  fprintf(fitsheader,"CUNIT1 = arcsec \n");
  fprintf(fitsheader,"CUNIT2 = arcsec \n");
  fprintf(fitsheader,"CTYPE1 = PAR \n");
  fprintf(fitsheader,"CTYPE2 = PAR \n");
  fprintf(fitsheader,"CRPIX1 = 0 \n");
  fprintf(fitsheader,"CRPIX2 = 0 \n");
  fprintf(fitsheader,"CRVAL1 = 0 \n");
  fprintf(fitsheader,"CRVAL2 = 0 \n");
  char bunit[]="JY/BEAM";
  fprintf(fitsheader,"BUNIT  = %s \n",bunit);
  char btype[]="intensity";
  fprintf(fitsheader,"BTYPE  = %s \n",btype);
  fclose(fitsheader);
}

