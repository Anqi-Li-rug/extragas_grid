/*
 *  extragas_fit.h
 *  Extragas
 *
 *  Created by Antonino Marasco on 27/11/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef EXTRAGAS_FIT_H
#define EXTRAGAS_FIT_H

#include <fstream>
#include <iostream>
#include <vector>
#include "arrays3d.h"
#include "array2d.h"
#include "nr3.h"

//using namespace std;

class extragas_fit
{
public:
	//constructor. It loads the data-cube and evaluate it's flux 
	extragas_fit(char* datacube_file, char* maskcube_file);
    ~extragas_fit(void);
	//create a single model based on param_file (and calculate the residual)
	//use the file extragas.in
	//use the file free_parameters.txt (only initial estimate)
    double modelmaker(char* param_file, bool normalize, bool printcube,char* outputfilename,double v_k1,double f_k1,double accre_k1);
	//the same, but it keeps on doing the same model and save that with the best residual
	void modelmaker_loop(char* param_file);
	
	//Downhill simplex algorithm for finding the best-fit model
	//use the file extragas.in
	//use the file free_parameters.txt (both initial estimate and initial displacements)
	void bestmodel_finder(char* param_file);
	
	//Markov-chain Monte Carlo for best-parameter finding
	double mcmc(char* param_file, long ncycle, double pixel_per_beam, double h, double v_k, double accre);
    
	
private:
	//return the residual |model-data|
	//double tominimize (VecDoub freepar);
	//write the free-parameters in a temporary input file that is read by extragas
    double tominimize (VecDoub freepar);
	void assign();
	
	//cube variables	
	//axis
	std::vector<double> RA, DEC, VEL;
	int X_SIZE, Y_SIZE, Z_SIZE, crpix1, crpix2, crpix3;
	//cube
	myarray::double3d datacube; double* datacube_ptr;
    myarray::double3d maskcube; double* maskcube_ptr;
	double dataflux, modelflux, fluxratio;
	//miscellaneous
	bool write_cube, use_blanks;
    //projection type
	char projectionType[10];
	//HI halo mass and accretion rate
	double mass_HIhalo, accretionrate;
	
	//parameter variables
	int nparam;
	std::vector<double> p, dp;	
	//flux-density threshold
	double residual_above;
	//to check if residual are getting better	
	bool getting_better;
	double residual_last;
	
	////////////////////////////////////////////////////////////	
	//VARIABLES for the DOWNHILL SIMPLEX ALGORITHM (from nr3) //
	////////////////////////////////////////////////////////////
	//some stuff (are they really useful?)
	Doub ftol; Int nfunc; Int mpts; Int ndim; Doub fmin; VecDoub y; MatDoub pnow;
	
	//MINIMISATION ROUTINES
	//Generic implementation, you must pass the coordinates of all N+1 points
	VecDoub downhill(MatDoub_I &p0);
	//pass 1 point only, the others are displaced by different amounts
	VecDoub downhill(VecDoub_I &point, VecDoub_I &dels);
	//pass 1 point only, the others are displaced by the same amounts
	VecDoub downhill(VecDoub_I &point, const Doub del);
	
	//utility 1	
	void get_psum(MatDoub_I &p, VecDoub_O &psum);
	//utility 2
	Doub amotry(MatDoub_IO &p, VecDoub_O &y, VecDoub_IO &psum, const Int ihi, const Doub fac);
	//logfile
	void logfile_update(char* file, int step, double myresidual);
	
	//////////////////////////
	//FUNCTIONS USED IN MCMC//
	//////////////////////////
	double assumed_error(VecDoub freepar);
	double prior(VecDoub xparam);
};

#endif
