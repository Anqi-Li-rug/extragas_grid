/*
 *  fitsutil.h
 *  Bulge
 *
 *  Created by Antonino Marasco on 16/04/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef FITSUTIL_H
#define FITSUTIL_H

#include <iostream>
#include <vector>
#include "arrays3d.h"
#include "fitsio.h"

using namespace std;

namespace fitsutil
{
	//read a 2D fit file
    void read_map(char* input_fits, std::vector<double> &x, std::vector<double> &y, double* &map_ptr);	

	//print a 2D fit file
	void print_map(char* output_fits, std::vector<double> &x, std::vector<double> &y, double* &map_ptr, char* ctype1, char* ctype2, char* cunit1, char* cunit2, char* bunit);
	
	//read a 3D fit file
	void read_cube(char* input_fits, std::vector<double> &x, std::vector<double> &y, std::vector<double> &z, int &crpix1, int &crpix2, int &crpix3, double* &cube_ptr);
	
	//print a 3D fit file
	void print_cube(char* output_fits, std::vector<double> &x, std::vector<double> &y, std::vector<double> &z, int crpix1, int crpix2, int crpix3, double* &cube_ptr, char* ctype1, char* ctype2, char* ctype3, char* cunit1, char* cunit2, char* cunit3, char* bunit);
	
	//smooth a 3D fit file
	void smooth_cube(myarray::double3d &cube, std::vector<double> &x, std::vector<double> &y, std::vector<double> &z, double oldx, double oldy, double newx, double newy, double sigma_accuracy);
	
	//smooth an all-sky 3D fit file. It's basically for the Milky Way only.
	void smooth_allsky_cube(myarray::double3d &cube, std::vector<double> &x, std::vector<double> &y, std::vector<double> &z, double oldx, double oldy, double newx, double newy, double sigma_accuracy);

	//smooth a cube in all the directions (rarely useful)
	void smooth3D(myarray::array3d<double> &cube, std::vector<double> &x, std::vector<double> &y, std::vector<double> &z, double beamx, double beamy, double beamz, double sigma_accuracy);
}
#endif