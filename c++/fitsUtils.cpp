/*
 *  fitsUtils.cpp
 *  Extragas
 *
 *  Created by Antonino Marasco on 07/10/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#include <iostream>
#include <cstring>
#include <cmath>
#include "fitsio.h"  
#include "fitsUtils.h"

using namespace std;
/*
void writefits_3D(const char *filename, float *outcube, const int xsize, const int ysize, const int zsize, const char *header_file)
{
  fitsfile *fptr;        //pointer to the FITS file; defined in fitsio.h 
  int status, ii, jj, kk;
  long  fpixel = 1, naxis = 3, nelements;
  //  long naxes[3] = {RAsize, DECsize, VELsize};
  //    float array[RAsize][DECsize][VELsize]; 
  long dnaxes[3] = {xsize, ysize, zsize};
  long naxes[3] = {zsize, ysize, xsize};
  float array[zsize][ysize][xsize];
  
  remove(filename);               // Delete old file if it already exists 
  status = 0;    // initialize status before calling fitsio routines 
  fits_create_file(&fptr, filename, &status);   // create new file
  
  // Create the primary array image (32-bit floating pixels 
  fits_create_img(fptr, FLOAT_IMG, naxis, dnaxes, &status);
    
  FILE* fitsheader=NULL;
  fitsheader=fopen(header_file, "r");
  char key[7];
  char eq[3];
  char cval[10];
  char c;
  while (c!=EOF){
    fscanf(fitsheader, "%s %s %s",key,eq,cval);
    float val=atof(cval);
    c=fgetc(fitsheader);
    if (val == 0){
          fits_update_key(fptr, TSTRING, key, &cval, " ", &status);
    }
    if (val != 0){
          fits_update_key(fptr, TFLOAT, key, &val, " ", &status);
    }
  }
  fclose(fitsheader);

  printf("Writing %i x %i x %i fits file \n", naxes[0], naxes[1], naxes[2]);

  // Initialize the values in the image with a linear ramp function 
  for (ii = 0; ii < naxes[0]; ii++)
    {
      for (jj = 0; jj < naxes[1]; jj++)
	{
	  for (kk = 0; kk < naxes[2]; kk++)
	    {
	       // This is the way the array is seen
		  //array[ii][jj][kk] = outcube[ii*naxes[1]*naxes[2]+jj*naxes[2]+kk];  
		  //In the declaration to write the fits file the axes are 
		  //swapped, as in fortran, ii->kk, kk->ii
		
	      array[ii][jj][kk] = outcube[kk*naxes[1]*naxes[0]+jj*naxes[0]+ii];
	    }
	}
    }
  
  nelements = naxes[0] * naxes[1] * naxes[2]; // number of pixels to write 
  
  fits_write_img(fptr, TFLOAT, fpixel, nelements, array[0], &status);

  fits_close_file(fptr, &status);            // close the file
  
  fits_report_error(stderr, status);  //print out any error messages
};
*/

void writefits_3D(const char *filename, float *outcube, const int xsize, const int ysize, const int zsize, const char *header_file)
{
	fitsfile *fptr;        //pointer to the FITS file; defined in fitsio.h 
  int status, i, j, k;
  long  fpixel = 1, naxis = 3, nelements;
  //  long naxes[3] = {RAsize, DECsize, VELsize};
  //    float array[RAsize][DECsize][VELsize]; 
  long dnaxes[3] = {xsize, ysize, zsize};
  long naxes[3] = {zsize, ysize, xsize};
  float* array_3d = NULL; //array dinamico
  array_3d = new float[xsize*ysize*zsize];
  
  remove(filename);               // Delete old file if it already exists 
  status = 0;    // initialize status before calling fitsio routines 
  fits_create_file(&fptr, filename, &status);   // create new file
  
  // Create the primary array image (32-bit floating pixels 
  fits_create_img(fptr, FLOAT_IMG, naxis, dnaxes, &status);
    
  FILE* fitsheader=NULL;
  fitsheader=fopen(header_file, "r");
  char key[20];
  char eq[20];
  char cval[20];
  char c;
  while (c!=EOF){
    fscanf(fitsheader, "%s %s %s",key,eq,cval);
    float val=atof(cval);
    c=fgetc(fitsheader);
    if (val == 0){fits_update_key(fptr, TSTRING, key, &cval, " ", &status);}
    if (val != 0){fits_update_key(fptr, TFLOAT, key, &val, " ", &status);}
  }
  fclose(fitsheader);

  printf("Writing %i x %i x %i fits file \n", naxes[0], naxes[1], naxes[2]);

  // Initialize the values in the image with a linear ramp function 
  for (k=0; k<zsize; k++)
  {
	for (j=0; j<ysize; j++)
	{
		for (i=0; i<xsize; i++) array_3d[i+j*xsize+k*xsize*ysize] = outcube[k+j*zsize+i*ysize*zsize];
	}
  }
  
  nelements = naxes[0] * naxes[1] * naxes[2]; // number of pixels to write 
  
  fits_write_img(fptr, TFLOAT, fpixel, nelements, array_3d, &status);

  fits_close_file(fptr, &status);            // close the file
  
  fits_report_error(stderr, status);  //print out any error messages
  delete[] array_3d;
  array_3d = NULL;
};

void writefits_2D(const char *filename, float *outimage, const int xsize, const int ysize, const char *header_file) 

{
	fitsfile *fptr;       /* pointer to the FITS file; defined in fitsio.h */
  int status, i, j;
  long  fpixel = 1, naxis = 2, nelements;
  long dnaxes[2] = {xsize, ysize};
  long naxes[2] = {xsize, ysize};
  float* array_2d = NULL; //array dinamico
	array_2d = new float[xsize*ysize];
  
  remove(filename);               /* Delete old file if it already exists */
  status = 0;    /* initialize status before calling fitsio routines */
  fits_create_file(&fptr, filename, &status);   /* create new file */
  
  /* Create the primary array image (32-bit floating pixels */
  fits_create_img(fptr, FLOAT_IMG, naxis, dnaxes, &status);
    
  FILE* fitsheader=NULL;
  fitsheader=fopen(header_file, "r");
  char key[20];
  char eq[20];
  char cval[20];
  char c;
  while (c!=EOF){
    fscanf(fitsheader, "%s %s %s",key,eq,cval);
    float val=atof(cval);
    c=fgetc(fitsheader);
    if (val == 0){
          fits_update_key(fptr, TSTRING, key, &cval, " ", &status);
    }
    if (val != 0){
          fits_update_key(fptr, TFLOAT, key, &val, " ", &status);
    }
  }
  fclose(fitsheader);

  printf("Writing %i x %i fits file \n", naxes[0], naxes[1]);

  /* Initialize the values in the image with a linear ramp function */
  for (j = 0; j < ysize; j++)
  {
      for (i = 0; i < xsize; i++)
	  {
		  array_2d[i+j*xsize] = outimage[j+i*ysize];
	  }
  }
  
  nelements = naxes[0] * naxes[1]; /* number of pixels to write */
  
  fits_write_img(fptr, TFLOAT, fpixel, nelements, array_2d, &status);

  fits_close_file(fptr, &status);            /* close the file */
  
  fits_report_error(stderr, status);  /* print out any error messages */
	delete[] array_2d;
	array_2d = NULL;
};

void writefits_2D_DP(const char *filename, double *outimage, const int xsize, const int ysize, const char *header_file) 

{
  fitsfile *fptr;       /* pointer to the FITS file; defined in fitsio.h */
  int status, ii, jj;
  long  fpixel = 1, naxis = 2, nelements;
  long dnaxes[2] = {xsize, ysize};
  long naxes[2] = {xsize, ysize};
  double array[ysize][xsize];
  
  remove(filename);               /* Delete old file if it already exists */
  status = 0;    /* initialize status before calling fitsio routines */
  fits_create_file(&fptr, filename, &status);   /* create new file */
  
  /* Create the primary array image (32-bit floating pixels */
  fits_create_img(fptr, FLOAT_IMG, naxis, dnaxes, &status);
    
  FILE* fitsheader=NULL;
  fitsheader=fopen(header_file, "r");
  char key[7];
  char eq[3];
  char cval[10];
  char c;
  while (c!=EOF){
    fscanf(fitsheader, "%s %s %s",key,eq,cval);
    double val=atof(cval);
    c=fgetc(fitsheader);
    if (val == 0){
          fits_update_key(fptr, TSTRING, key, &cval, " ", &status);
    }
    if (val != 0){
          fits_update_key(fptr, TDOUBLE, key, &val, " ", &status);
    }
  }
  fclose(fitsheader);

  printf("Writing %i x %i fits file \n", naxes[0], naxes[1]);

  /* Initialize the values in the image with a linear ramp function */
  for (ii = 0; ii < naxes[0]; ii++)
    {
      for (jj = 0; jj < naxes[1]; jj++)
	{
	  //*(AutoArray + 2 * naxes[1] + 3)  *(*(DynArray+3)+4)
	  array[jj][ii] = outimage[ii*naxes[1]+jj];
	}
    }
  
  nelements = naxes[0] * naxes[1]; /* number of pixels to write */
  
  fits_write_img(fptr, TDOUBLE, fpixel, nelements, array[0], &status);

  fits_close_file(fptr, &status);            /* close the file */
  
  fits_report_error(stderr, status);  /* print out any error messages */
};


void readfits_2D(const char *filename, float *image, int xsize, int ysize) 
{
    fitsfile *fptr2;       /* pointer to the FITS file, defined in fitsio.h */
    int status,  nfound, anynull;
    long naxes[2], fpixel, nelements, npixels;

#define imagesize xsize*ysize
    float nullval;
    float timage[ysize][xsize];
  
    status = 0;

    if ( fits_open_file(&fptr2, filename, READONLY, &status) );
    fits_report_error(stderr, status);  /* print out any error messages */
    
    /* read the NAXIS1 and NAXIS2 keyword to get image size */
    if ( fits_read_keys_lng(fptr2, "NAXIS", 1, 2, naxes, &nfound, &status) );
    fits_report_error(stderr, status);  /* print out any error messages */

    npixels  = naxes[0] * naxes[1];         /* number of pixels in the image */
    fpixel   = 1;
    nullval  = 0;                /* don't check for null values in the image */


    while (npixels > 0) {
      nelements = npixels;
      if (npixels > imagesize)
	nelements = imagesize;     /* read as many pixels as will fit in buffer */
      /* Note that even though the FITS images contains unsigned integer */
      /* pixel values (or more accurately, signed integer pixels with    */
      /* a bias of 32768),  this routine is reading the values into a    */
      /* float array.   Cfitsio automatically performs the datatype      */
      /* conversion in cases like this.                                  */
      
      if (fits_read_img(fptr2, TFLOAT, fpixel, nelements, &nullval,
			timage, &anynull, &status) );
      fits_report_error(stderr, status);  /* print out any error messages */
      
      npixels -= nelements;    /* increment remaining number of pixels */
      fpixel  += nelements;    /* next pixel to be read in image */
    }
    for (int jj = 0; jj < naxes[1]; jj++) {
      for (int ii = 0; ii < naxes[0]; ii++) {
       	image[ii*naxes[1]+jj]= timage[jj][ii];
      }
    }
    printf("Reading %i x %i fits file \n", naxes[0], naxes[1]);

    if ( fits_close_file(fptr2, &status) );
    fits_report_error(stderr, status);  /* print out any error messages */
    
    return;
};


void readfits_3D(const char *filename, float *cube, int xsize, int ysize, int zsize)

{
    //streambuf* orig_buf = cout.rdbuf();
    //cout.rdbuf(NULL);
    fitsfile *fptr3;
    int status,  nfound, anynull;
    long naxes[3], fpixel, nelements, npixels;

    #define cubesize xsize*ysize*zsize
    float nullval;
    float tcube[zsize][ysize][xsize];
    //float tcube[xsize][ysize][zsize];

    status = 0;

    if ( fits_open_file(&fptr3, filename, READONLY, &status) );
    fits_report_error(stderr, status);


    if ( fits_read_keys_lng(fptr3, "NAXIS", 1, 3, naxes, &nfound, &status) );
    fits_report_error(stderr, status);

    npixels  = naxes[0] * naxes[1] * naxes[2];   
    fpixel   = 1;
    nullval  = 0;

    while (npixels > 0)
      {
	nelements = npixels;
	if (npixels > cubesize)
	  nelements = cubesize;
	
	if ( fits_read_img(fptr3, TFLOAT, fpixel, nelements, &nullval,
			   tcube, &anynull, &status) );
	fits_report_error(stderr, status); 
	
	npixels -= nelements; 
	fpixel  += nelements; 
      }
    
    for (int kk = 0; kk < naxes[2]; kk++) {
      for (int jj = 0; jj < naxes[1]; jj++) {
	for (int ii = 0; ii < naxes[0]; ii++) {
	  /* This is the way the array is seen
	     array[ii][jj][kk] = outcube[ii*naxes[1]*naxes[2]+jj*naxes[2]+kk];  
	     In the declaration to write the fits file the axes are 
	     swapped, as in fortran, ii->kk, kk->ii
	  */
	  //	  cube[kk*naxes[1]*naxes[0]+jj*naxes[0]+ii]= tcube[ii][jj][kk];
	  cube[ii*naxes[2]*naxes[1]+jj*naxes[2]+kk]= tcube[kk][jj][ii];
	}
      }
    }
    streambuf* orig_buf1 = cout.rdbuf();
    cout.rdbuf(NULL);
    printf("Reading %i x %i x %i fits file \n", naxes[0], naxes[1], naxes[2]);

    if ( fits_close_file(fptr3, &status) );
    fits_report_error(stderr, status);

    return;
};


void printHeader3D(int RAsize, int DECsize, int VELsize,
		   float cdelt1, float cdelt2, float cdelt3,
		   int crpix1, int crpix2, int crpix3,
		   float crval1, float crval2, float crval3)
{
  FILE* fitsheader=NULL;
  fitsheader=fopen("header.tab", "w");

  fprintf(fitsheader,"CDELT1 = %f \n",cdelt1);
  fprintf(fitsheader,"CDELT2 = %f \n",cdelt2);
  fprintf(fitsheader,"CDELT3 = %f \n",cdelt3);
  fprintf(fitsheader,"CUNIT1 = DEGREE \n");
  fprintf(fitsheader,"CUNIT2 = DEGREE \n");
  fprintf(fitsheader,"CUNIT3 = KM/S \n");
  fprintf(fitsheader,"CTYPE1 = RA---NCP \n");
  fprintf(fitsheader,"CTYPE2 = DEC--NCP \n");
  fprintf(fitsheader,"CTYPE3 = VELO-HEL\n");
  fprintf(fitsheader,"CRPIX1 = %i \n",crpix1);
  fprintf(fitsheader,"CRPIX2 = %i \n",crpix2);
  fprintf(fitsheader,"CRPIX3 = %i \n",crpix3);
  fprintf(fitsheader,"CRVAL1 = %f \n",crval1);
  fprintf(fitsheader,"CRVAL2 = %f \n",crval2);
  fprintf(fitsheader,"CRVAL3 = %f \n",crval3);
  char bunit[]="JY/BEAM";
  fprintf(fitsheader,"BUNIT  = %s \n",bunit);
  char btype[]="intensity";
  fprintf(fitsheader,"BTYPE  = %s \n",btype);
  fclose(fitsheader);
};

void printHeader3D_AllSky(int RAsize, int DECsize, int VELsize,
			  float cdelt1, float cdelt2, float cdelt3,
			  int crpix1, int crpix2, int crpix3,
			  float crval1, float crval2, float crval3,
			  char bunit[])
{
  FILE* fitsheader=NULL;
  fitsheader=fopen("header.tab", "w");

  fprintf(fitsheader,"CDELT1 = %f \n",cdelt1);
  fprintf(fitsheader,"CDELT2 = %f \n",cdelt2);
  fprintf(fitsheader,"CDELT3 = %f \n",cdelt3);
  fprintf(fitsheader,"CUNIT1 = DEGREE \n");
  fprintf(fitsheader,"CUNIT2 = DEGREE \n");
  fprintf(fitsheader,"CUNIT3 = KM/S \n");
  fprintf(fitsheader,"CTYPE1 = RA---NCP \n");
  fprintf(fitsheader,"CTYPE2 = DEC--NCP \n");
  fprintf(fitsheader,"CTYPE3 = VELO-HEL \n");
  fprintf(fitsheader,"CRPIX1 = %i \n",crpix1);
  fprintf(fitsheader,"CRPIX2 = %i \n",crpix2);
  fprintf(fitsheader,"CRPIX3 = %i \n",crpix3);
  fprintf(fitsheader,"CRVAL1 = %f \n",crval1);
  fprintf(fitsheader,"CRVAL2 = %f \n",crval2);
  fprintf(fitsheader,"CRVAL3 = %f \n",crval3);
  fprintf(fitsheader,"BUNIT  = %s \n",bunit);
  char btype[]="intensity";
  fprintf(fitsheader,"BTYPE  = %s \n",btype);
  fclose(fitsheader);
};

void printHeader2D(int RAsize, int DECsize,
		   float cdelt1, float cdelt2,
		   int crpix1, int crpix2,
		   float crval1, float crval2)
{
  FILE* fitsheader=NULL;
  fitsheader=fopen("header.tab", "w");

  fprintf(fitsheader,"CDELT1 = %f \n",cdelt1);
  fprintf(fitsheader,"CDELT2 = %f \n",cdelt2);
  fprintf(fitsheader,"CUNIT1 = DEGREE \n");
  fprintf(fitsheader,"CUNIT2 = DEGREE \n");
  fprintf(fitsheader,"CTYPE1 = RA---NCP \n");
  fprintf(fitsheader,"CTYPE2 = DEC--NCP \n");
  fprintf(fitsheader,"CRPIX1 = %i \n",crpix1);
  fprintf(fitsheader,"CRPIX2 = %i \n",crpix2);
  fprintf(fitsheader,"CRVAL1 = %f \n",crval1);
  fprintf(fitsheader,"CRVAL2 = %f \n",crval2);
  char bunit[]="JY/BEAM";
  fprintf(fitsheader,"BUNIT  = %s \n",bunit);
  char btype[]="intensity";
  fprintf(fitsheader,"BTYPE  = %s \n",btype);
  fclose(fitsheader);
};

void printHeader2D_AllSky(int RAsize, int DECsize,
			  float cdelt1, float cdelt2,
			  int crpix1, int crpix2,
			  float crval1, float crval2,
			  char bunit[])
{
  FILE* fitsheader=NULL;
  fitsheader=fopen("header.tab", "w");

  fprintf(fitsheader,"CDELT1 = %f \n",cdelt1);
  fprintf(fitsheader,"CDELT2 = %f \n",cdelt2);
  fprintf(fitsheader,"CUNIT1 = DEGREE \n");
  fprintf(fitsheader,"CUNIT2 = DEGREE \n");
  fprintf(fitsheader,"CTYPE1 = RA---NCP \n");
  fprintf(fitsheader,"CTYPE2 = DEC--NCP \n");
  fprintf(fitsheader,"CRPIX1 = %i \n",crpix1);
  fprintf(fitsheader,"CRPIX2 = %i \n",crpix2);
  fprintf(fitsheader,"CRVAL1 = %f \n",crval1);
  fprintf(fitsheader,"CRVAL2 = %f \n",crval2);
  fprintf(fitsheader,"BUNIT  = %s \n",bunit);
  char btype[]="intensity";
  fprintf(fitsheader,"BTYPE  = %s \n",btype);
  fclose(fitsheader);
};
