/*
 *  fitsUtils.h
 *  Extragas
 *
 *  Created by Antonino Marasco on 07/10/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef FRATERNALI_FITSUTIL
#define FRATERNALI_FITSUTIL

using namespace std;

void writefits_3D(const char *filename, float *outcube, const int xsize, const int ysize, const int zsize, const char *header_file);
void writefits_2D(const char *filename, float *outimage, const int xsize, const int ysize, const char *header_file);
void writefits_2D_DP(const char *filename, double *outimage, const int xsize, const int ysize, const char *header_file);
void readfits_2D(const char *filename, float *image, int xsize, int ysize);
void readfits_3D(const char *filename, float *cube, int xsize, int ysize, int zsize);
void printHeader3D(int RAsize, int DECsize, int VELsize, float cdelt1, float cdelt2, float cdelt3, int crpix1, int crpix2, int crpix3, float crval1, float crval2, float crval3);
void printHeader3D_AllSky(int RAsize, int DECsize, int VELsize, float cdelt1, float cdelt2, float cdelt3, int crpix1, int crpix2, int crpix3, float crval1, float crval2, float crval3, char bunit[]);
void printHeader2D(int RAsize, int DECsize, float cdelt1, float cdelt2, int crpix1, int crpix2, float crval1, float crval2);
void printHeader2D_AllSky(int RAsize, int DECsize, float cdelt1, float cdelt2, int crpix1, int crpix2, float crval1, float crval2, char bunit[]);

#endif