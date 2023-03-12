/*
 *  extragas_fit.cpp
 *  Extragas
 *
 *  Created by Antonino Marasco on 27/11/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *
 */

//#include "extragas_fit.h"
#include "extragasFits.h"
#include "fitsutil.h" 
#include "fitsUtils.h"
#include "findValue.h"
#include "array.h"
#include "arrays3d.h"
#include "extragas.h"
#include "statist.h"
#include "myutilities.h"
#include <fstream>
#include "zap.h"
#include "extragas_fit.h"
extern bool verbose;
extern long int seed;
using namespace std;

//GLOBAL VARIABLES
char* tempfile = "/net/cuby3.apertif/data/users/li/extragas/files/temp_extragas.in";
extern char inputfile[];
extern "C"
{
extragas_fit::extragas_fit(char* datacube_file,char* maskcube_file)
{
    //streambuf* orig_buf = cout.rdbuf();
    //cout.rdbuf(NULL);
	//NUMBER OF FREE-PARAMETERS
	nparam = 3;
   // int model_n;
    //string prefix="/net/cuby3.apertif/data/users/li/extragas/cubes/model_n";
    //stringstream oss1;
    //oss1 <<prefix<<model_n<<".fits";
    //string outputfilename(oss1.str());
	//cout<<"Reading the datacube file "<<datacube_file<<endl;
	fitsutil::read_cube(maskcube_file,RA,DEC,VEL,crpix1,crpix2,crpix3,maskcube_ptr);
	fitsutil::read_cube(datacube_file,RA,DEC,VEL,crpix1,crpix2,crpix3,datacube_ptr);
	X_SIZE = RA.size(); Y_SIZE = DEC.size(); Z_SIZE = VEL.size();
	datacube.init(datacube_ptr, X_SIZE, Y_SIZE, Z_SIZE);
    maskcube.init(maskcube_ptr, X_SIZE, Y_SIZE, Z_SIZE);
	
	//evaluate the totalflux of the datacube
	findValue(inputfile, projectionType, "projectionType");
	dataflux = 0;
	if (strcoll(projectionType,"external") ==0)
	{
        for(int i=0;i<X_SIZE; i++) for(int j=0;j<Y_SIZE; j++) for(int k=0;k<Z_SIZE; k++) if((!isnan(datacube(i,j,k)))&((!isnan(maskcube(i,j,k))))&(maskcube(i,j,k)==1)) {dataflux+=datacube(i,j,k);}
	}
	else if (strcoll(projectionType,"internal") ==0)
	{
		double cosb;
		for(int i=0;i<X_SIZE;i++) for(int j=0;j<Y_SIZE;j++) for(int k=0;k<Z_SIZE;k++)
		{
			if(!isnan(datacube(i,j,k)))
			{
				cosb = cos(DEC[j]*0.01745329252);
				dataflux += datacube(i,j,k)*cosb;
			}
		}
	}
	//cout<<"Total datacube flux = "<<dataflux<<endl<<endl;
	residual_last = 1e20; getting_better = false;
}

void extragas_fit::assign()
{	
	//ofstream temp(tempfile);
	//temp<<"h_v_k  "<<p[0]<<endl;
	//temp<<"ion_frac  "<<p[1]<<endl;
	//temp<<"alpha_accr  "<<p[2]<<endl;
	//temp.close();
}

double extragas_fit::modelmaker(char* param_file, bool normalize, bool printcube,char* outputfilename,double v_k1,double f_k1,double accre_k1)
{
	p.resize(nparam); dp.resize(nparam);
        double rms_noise=0.00019;
	//if(verbose)cout<<endl<<"CREATING THE MODELCUBE"<<endl<<endl;
	//if(verbose)cout<<endl<<"Reading the parameter file "<<param_file<<endl;
	findValue(param_file, p[0], "vkick");
	findValue(param_file, p[1], "fion");
	findValue(param_file, p[2], "alpha_accr");
    p[0]=v_k1;
    p[1]=f_k1;
    p[2]=accre_k1;
    double ratio_k1=pow(10.0,accre_k1);
    ratio_k1=1.0;
    //char outputfilename[60];
    //findValue(inputfile,outputfilename,"outputfilename");
    double mass_HIhalo=0.1;
    double accretionrate=0.1;
    //cout<<mass_HIhalo<<accretionrate<<endl;
	assign();
    //cout<<mass_HIhalo<<endl;
    //cout<<accretionrate<<endl;
	ffarray::Array3D<float> modelcube = extragas_model(p[0],p[1],p[2], mass_HIhalo, accretionrate);
	//normalize the model to the same flux of the data
    float ptt1=0;
	if(normalize)
    {
        if(verbose)cout<<"Normalizing the modelcube to the same flux of the datacube"<<endl;
        modelflux = 0;
        if (strcoll(projectionType,"external") ==0)
        {
            for(int i=0;i<X_SIZE; i++) for(int j=0;j<Y_SIZE; j++) for(int k=0;k<Z_SIZE; k++) if((!isnan(modelcube[i][j][k]))&((!isnan(maskcube(i,j,k))))&(maskcube(i,j,k)==1)) {modelflux+=modelcube[i][j][k];}
        }
        
        else if (strcoll(projectionType,"internal") ==0)
        {
            double cosb;
            for(int i=0;i<X_SIZE;i++) for(int j=0;j<Y_SIZE;j++) for(int k=0;k<Z_SIZE;k++)
            {
                if(!isnan(datacube(i,j,k)))
                {
                    cosb = cos(DEC[j]*0.01745329252);
                    modelflux += modelcube[i][j][k]*cosb;
                }
            }
        }
        fluxratio = dataflux/modelflux;
        mass_HIhalo*=fluxratio;
        accretionrate*=fluxratio;
        cout<<"nonormdata to model flux ratio = "<<fluxratio<<endl;
        for(int i=0;i<X_SIZE;i++) for(int j=0;j<Y_SIZE;j++) for(int k=0;k<Z_SIZE;k++) modelcube[i][j][k]*=1;
    }
	
	double* cube_ptr = new double [X_SIZE*Y_SIZE*Z_SIZE];
	myarray::double3d cube; cube.init(cube_ptr, X_SIZE, Y_SIZE, Z_SIZE);
	for(int i=0;i<X_SIZE;i++) for(int j=0;j<Y_SIZE;j++) for(int k=0;k<Z_SIZE;k++) cube(i,j,k)=modelcube[i][j][k];
	//evaluate the residual
	findValue(param_file, residual_above, "residual_above");
	findValue(param_file, use_blanks, "use_blanks");
	
	//evaluate the residual
	if(verbose)cout<<"Evaluating the residuals..."<<endl;
	
	double residual = 0, weight, cosb, beam=2;
	for(int i=0;i<X_SIZE;i++) for(int j=0;j<Y_SIZE;j++) for(int k=0;k<Z_SIZE;k++)
	{
		cosb = cos(DEC[j]*0.01745329252);
		weight = cosb;
                residual_above=-100000000;
		if((datacube(i,j,k)>residual_above || modelcube[i][j][k]>residual_above))
		{
			if((!isnan(modelcube[i][j][k]))&(!isnan(maskcube(i,j,k)))&((modelcube[i][j][k]*ratio_k1>2.*rms_noise)|(datacube(i,j,k)>2.*rms_noise))&(maskcube(i,j,k)>0)) {residual += pow(fabs(datacube(i,j,k) - ratio_k1*modelcube[i][j][k]),1);}
			else 
			{
				if(use_blanks) {residual += weight*modelcube[i][j][k];}
				else {residual += 0;}
			}
		}
	}
	if(isnan(residual) || isinf(residual)) residual=1e20;
	if(verbose)cout<<"Residual = "<<residual<<endl;
		if(printcube)
        {
            if (strcoll(projectionType,"external") ==0)
            {
            fitsutil::print_cube(outputfilename,RA,DEC,VEL,crpix1,crpix2,crpix3,cube_ptr,"RA---NCP","DEC--NCP","VELO-HEL","DEGREE","DEGREE","KM/S","JY/BEAM");
            }
            else if (strcoll(projectionType,"internal") ==0)
            {
                for(int i=0;i<X_SIZE;i++) RA[i]+=180.;
                fitsutil::print_cube("/net/cuby3.apertif/data/users/li/extragas/cubes/modelcube_halo.fits",RA,DEC,VEL,crpix1,crpix2,crpix3,cube_ptr,"GLON-CAR","GLAT-CAR","VELO-LSR","DEGREE","DEGREE","M/S","K");
            }
        }
	zaparr(cube_ptr);
    return residual;
}

double extragas_fit::tominimize (VecDoub freepar)
{
	getting_better = false;
	
	//read the parameters
	int i,j = 0;
	for(i=0;i<p.size();i++)
	{
		if(dp[i]!=0)
		{
			p[i] = freepar[j];
			j++;
		}
	}
	assign();
	
	//make the model
	if(verbose) cout<<"Creating the modelcube"<<endl;
	ffarray::Array3D<float> modelcube = extragas_model(p[0],p[1],p[2], mass_HIhalo, accretionrate);
	
	//check that model and data have the same size
	if(X_SIZE!=modelcube.dim1() || Y_SIZE!=modelcube.dim2() || Z_SIZE!=modelcube.dim3())
	{
		cout<<"ERROR! DATACUBE AND MODELCUBE HAVE DIFFERENT SIZES!"<<endl;
		return NAN;
	}	
	
	//normalize the model to the same flux of the data
	if(verbose)cout<<"Normalizing the modelcube to the same flux of the datacube"<<endl;
	modelflux = 0;
	if (strcoll(projectionType,"external") ==0)
	{
		for(int i=0;i<X_SIZE; i++) for(int j=0;j<Y_SIZE; j++) for(int k=0;k<Z_SIZE; k++)
            if((!isnan(modelcube[i][j][k]))&((!isnan(maskcube(i,j,k))))&(maskcube(i,j,k)==1))
            modelflux+=modelcube[i][j][k];
	}
	else if (strcoll(projectionType,"internal") ==0)
	{
		double cosb;
		for(int i=0;i<X_SIZE;i++) for(int j=0;j<Y_SIZE;j++) for(int k=0;k<Z_SIZE;k++)
		{
			if(!isnan(datacube(i,j,k)))
			{
				cosb = cos(DEC[j]*0.01745329252);
				modelflux += modelcube[i][j][k]*cosb;
			}
		}
	}
	fluxratio = dataflux/modelflux;
	mass_HIhalo*=fluxratio;
	accretionrate*=fluxratio;
	cout<<"		data to model flux ratio = "<<fluxratio<<endl; 
	for(int i=0;i<X_SIZE;i++) for(int j=0;j<Y_SIZE;j++) for(int k=0;k<Z_SIZE;k++) modelcube[i][j][k]*=fluxratio;
	
	//evaluate the residual
	if(verbose)cout<<"Evaluating the residuals..."<<endl;
	
	double residual = 0, weight, cosb;
	for(int k=0;k<Z_SIZE;k++) 
	{
		for(int i=0;i<X_SIZE;i++) for(int j=0;j<Y_SIZE;j++) 
		{
			cosb = cos(DEC[j]*0.01745329252);
			//weight = 1./(1.4*beam*beam/cosb);//1.4 for hanning smoothing
			weight = cosb;
            double maskncount=0;
			if((datacube(i,j,k)>residual_above || modelcube[i][j][k]>residual_above))
			{
                if((!isnan(modelcube[i][j][k]))&((!isnan(maskcube(i,j,k))))&(maskcube(i,j,k)==1)) {residual += weight*pow(fabs(datacube(i,j,k) - modelcube[i][j][k]),1);maskncount++;}
				else 
				{
					if(use_blanks) {residual += weight*modelcube[i][j][k];}
					else {residual += 0;}
				}
			}
		}
	}
	if(isnan(residual) || isinf(residual)) residual=1e20;
	if(verbose)cout<<"Residual = "<<residual<<endl;
	
	//update residuals and print the model if they are improving. In this way you always have the best model in real time!
	if(residual<residual_last) {residual_last = residual; getting_better = true;}
	if(getting_better)
	{
		if(dp[0]!=0)  cout<<"vkick			[km/s]  = "<<p[0]<<endl;
		if(dp[1]!=0)  cout<<"fion			        = "<<p[1]<<endl;
		if(dp[2]!=0)  cout<<"alpha_accr	    [1/Gyr] = "<<p[2]<<endl;
					  cout<<"HI halo mass   [Mo]    = "<<mass_HIhalo<<endl;
					  cout<<"accretion rate [Mo/yr] = "<<accretionrate<<endl;
                bool write_cube=false;
		if(write_cube)
		{
			double* cube_ptr = new double [X_SIZE*Y_SIZE*Z_SIZE];
			myarray::double3d cube; cube.init(cube_ptr, X_SIZE, Y_SIZE, Z_SIZE);
			for(int i=0;i<X_SIZE;i++) for(int j=0;j<Y_SIZE;j++) for(int k=0;k<Z_SIZE;k++) cube(i,j,k)=modelcube[i][j][k];
			fitsutil::print_cube("modelcube_temp.fits",RA,DEC,VEL,crpix1,crpix2,crpix3,cube_ptr,"RA---NCP","DEC--NCP","VELO-HEL","DEGREE","DEGREE","KM/S","JY/BEAM");
			zaparr(cube_ptr);
		}
	}
	return residual;
}


void extragas_fit::modelmaker_loop(char* param_file)
{
	p.resize(nparam); dp.resize(nparam);
    char outputfilename[60];
	if(verbose)cout<<endl<<"CREATING THE MODELCUBE"<<endl<<endl;
	if(verbose)cout<<endl<<"Reading the parameter file "<<param_file<<endl;
	findValue(param_file, p[0], "vkick");
	findValue(param_file, p[1], "fion");
	findValue(param_file, p[2], "alpha_accr");
    findValue(inputfile,outputfilename,"outputfilename");
	assign();
	
	double residual_best = 1e20;
	
	while(1)
	{
		ffarray::Array3D<float> modelcube = extragas_model(p[0],p[1],p[2], mass_HIhalo, accretionrate);
		//normalize the model to the same flux of the data
		if(verbose)cout<<"Normalizing the modelcube to the same flux of the datacube"<<endl;
		modelflux = 0;
		if (strcoll(projectionType,"external") ==0)
		{
			for(int i=0;i<X_SIZE; i++) for(int j=0;j<Y_SIZE; j++) for(int k=0;k<Z_SIZE; k++) if(!isnan(datacube(i,j,k))) modelflux+=modelcube[i][j][k];
		}
		else if (strcoll(projectionType,"internal") ==0)
		{
			double cosb;
			for(int i=0;i<X_SIZE;i++) for(int j=0;j<Y_SIZE;j++) for(int k=0;k<Z_SIZE;k++)
			{
				if(!isnan(datacube(i,j,k)))
				{
					cosb = cos(DEC[j]*0.01745329252);
					modelflux += modelcube[i][j][k]*cosb;
				}
			}
		}
		
		fluxratio = dataflux/modelflux;
		mass_HIhalo*=fluxratio;
		accretionrate*=fluxratio;
		if(verbose)cout<<"		data to model flux ratio = "<<fluxratio<<endl; 
		for(int i=0;i<X_SIZE;i++) for(int j=0;j<Y_SIZE;j++) for(int k=0;k<Z_SIZE;k++) modelcube[i][j][k]*=fluxratio;
		
		double* cube_ptr = new double [X_SIZE*Y_SIZE*Z_SIZE];
		myarray::double3d cube; cube.init(cube_ptr, X_SIZE, Y_SIZE, Z_SIZE);
		for(int i=0;i<X_SIZE;i++) for(int j=0;j<Y_SIZE;j++) for(int k=0;k<Z_SIZE;k++) cube(i,j,k)=modelcube[i][j][k];
		
		//evaluate the residual
		findValue(param_file, residual_above, "residual_above");
		findValue(param_file, use_blanks, "use_blanks");
		
		//evaluate the residual
		if(verbose)cout<<"Evaluating the residuals..."<<endl;
		
		double residual = 0, weight, cosb;
		for(int i=0;i<X_SIZE;i++) for(int j=0;j<Y_SIZE;j++) for(int k=0;k<Z_SIZE;k++)
		{
			cosb = cos(DEC[j]*0.01745329252);
			weight = cosb;
			if((datacube(i,j,k)>residual_above || modelcube[i][j][k]>residual_above))
			{
				if(!isnan(datacube(i,j,k))) {residual += weight*pow(fabs(datacube(i,j,k) - modelcube[i][j][k]),1);}
				else 
				{
					if(use_blanks) {residual += weight*modelcube[i][j][k];}
					else {residual += 0;}
				}
			}
		}
		if(isnan(residual) || isinf(residual)) residual=1e20;		
		if(residual<residual_best)
		{
			residual_best = residual;
			cout<<"A better model has been found, with residual "<<residual_best<<endl;
			//fitsutil::print_cube(outputfilename,RA,DEC,VEL,crpix1,crpix2,crpix3,cube_ptr,"RA---NCP","DEC--NCP","VELO-HEL","DEGREE","DEGREE","m/s","K");
		}
		zaparr(cube_ptr);
	}
}


double extragas_fit::assumed_error(VecDoub freepar)
{	
	//read the parameters
	int i,j = 0;
	for(i=0;i<p.size();i++)
	{
		if(dp[i]!=0)
		{
			p[i] = freepar[j];
			j++;
		}
	}
	assign();
	
	//make the model
	if(verbose) cout<<"Creating the modelcube"<<endl;
	ffarray::Array3D<float> modelcube = extragas_model(p[0],p[1],p[2], mass_HIhalo, accretionrate);
	
	//check that model and data have the same size
	if(X_SIZE!=modelcube.dim1() || Y_SIZE!=modelcube.dim2() || Z_SIZE!=modelcube.dim3())
	{
		cout<<"ERROR! DATACUBE AND MODELCUBE HAVE DIFFERENT SIZES!"<<endl;
		return NAN;
	}	
	
	//normalize the model to the same flux of the data
	if(verbose)cout<<"Normalizing the modelcube to the same flux of the datacube"<<endl;
	modelflux = 0;
	if (strcoll(projectionType,"external") ==0)
	{
		for(int i=0;i<X_SIZE; i++) for(int j=0;j<Y_SIZE; j++) for(int k=0;k<Z_SIZE; k++) if(!isnan(datacube(i,j,k))) modelflux+=modelcube[i][j][k];
	}
	else if (strcoll(projectionType,"internal") ==0)
	{
		double cosb;
		for(int i=0;i<X_SIZE;i++) for(int j=0;j<Y_SIZE;j++) for(int k=0;k<Z_SIZE;k++)
		{
			if(!isnan(datacube(i,j,k)))
			{
				cosb = cos(DEC[j]*0.01745329252);
				modelflux += modelcube[i][j][k]*cosb;
			}
		}
	}
	
	fluxratio = dataflux/modelflux;
	mass_HIhalo*=fluxratio;
	accretionrate*=fluxratio;
	//if(verbose)cout<<"= "<<fluxratio<<endl; 
	for(int i=0;i<X_SIZE;i++) for(int j=0;j<Y_SIZE;j++) for(int k=0;k<Z_SIZE;k++) modelcube[i][j][k]*=fluxratio;
	
	//evaluate the residual	
	double residual = 0;
	long npoints = 0;
	for(int k=0;k<Z_SIZE;k++) 
	{
		for(int i=0;i<X_SIZE;i++) for(int j=0;j<Y_SIZE;j++) 
		{
			if((datacube(i,j,k)>residual_above || modelcube[i][j][k]>residual_above))
			{
				npoints++;
				if(!isnan(datacube(i,j,k))) {residual += pow(fabs(datacube(i,j,k) - modelcube[i][j][k]),1);}
				else 
				{
					if(use_blanks) {residual += modelcube[i][j][k];}
					else {residual += 0;}
				}
			}
		}
	}
	double Ndof = npoints - freepar.size();
	return residual/Ndof;
}

double extragas_fit::prior(VecDoub xparam)
{
	//demanding that the parameters make physically sense 
	double prior_1 = 1;
	int dim = xparam.size();
	VecDoub probability(dim);
	for(int i=0;i<dim;i++)				probability[i] = 1;
	if (xparam[0]<0 || xparam[0]>600)	probability[0] = 0;
	if (xparam[1]<0 || xparam[1]>1)		probability[1] = 0;
	if (xparam[2]<0)                    probability[2] = 0;
	for(int i=0;i<dim;i++) prior_1*=probability[i];
	
	return prior_1;
}

double extragas_fit::mcmc(char* param_file, long ncycle, double pixel_per_beam,double v_k, double h, double accre)
{
    cout<<"------------"<<endl;
    cout<<"MCMC AT WORK"<<endl;
    cout<<"------------"<<endl;
    //some variables used by mcmc
    double residual, residual_old, prob_ratio, acceptance;
    const double residual_weight = 1./(pixel_per_beam*pixel_per_beam);
    VecDoub xparam(3);
    xparam[0]=v_k;
    xparam[1]=h;
    xparam[2]=accre;
    residual = tominimize(xparam);
    return residual;
    cout<<"    residual : "<<residual<<endl;
    cout<<"    evaluating the assumed error on the data..."<<endl;
    double assumederror = assumed_error(xparam);
    cout<<"    assumed error : "<<assumederror<<" K"<<endl;
    double finalerror = assumederror;
    cout<<"Final error used in the residual calculation: "<<finalerror<<endl;
    double prob = exp(residual_weight*(-residual)/finalerror);
    double lnprob = log(prior(xparam))+residual_weight*(-residual)/finalerror;
    cout<<lnprob<<endl;
    return lnprob;
    cout<<"Number of iterations: "<<ncycle<<endl;
    
    p.resize(nparam); dp.resize(nparam);
    
    //cout<<endl<<"    Reading the parameter file "<<param_file<<endl;
    findValue(param_file, p[0], "vkick");
    findValue(param_file, p[1], "fion");
    findValue(param_file, p[2], "alpha_accr");
    
    findValue(param_file, dp[0], "vkick_delta");
    findValue(param_file, dp[1], "fion_delta");
    findValue(param_file, dp[2], "alpha_accr_delta");
    
    findValue(param_file, residual_above, "residual_above");
    findValue(param_file, ftol, "ftol");
    findValue(param_file, write_cube, "write_cube");
    findValue(param_file, use_blanks, "use_blanks");
    
    //const double residual_weight = 1./(pixel_per_beam*pixel_per_beam);
    
    //number of free parameters
    int nfreeparam = 0;
    for(int i=0; i<nparam; i++) if(dp[i]!=0) nfreeparam+=1;
    if(nfreeparam<=0)
    {
        cout<<"No free parameters. Nothing to do, apparently!"<<endl;
        return -INFINITY;
    }
    else {cout<<"    number of free parameters: "<<nfreeparam<<endl;}
    
    //initialise the initial points
    VecDoub  xparam_old;
    VecDoub xdelta(nfreeparam);
    int j=0;
    for(int i=0;i<nparam;i++)
    {
        if(dp[i]!=0)
        {
            xparam[j]=p[i];
            xdelta[j]=dp[i];
            j++;
        }
    }
    myarray::double2d mcmc_all(nfreeparam,ncycle+1);
    for(int i=0;i<nfreeparam;i++) for(int j=0;j<ncycle+1;j++) mcmc_all(i,j)=NAN;
    for(int i=0;i<nfreeparam;i++) mcmc_all(i,0) = xparam[i];
    
    //residual of the best-fit model
    cout<<"    evaluating the residual for the fiducial model..."<<endl;
    //residual = tominimize(xparam);
    //cout<<"    residual : "<<residual<<endl;
    //cout<<"    evaluating the assumed error on the data..."<<endl;
    //double assumederror = assumed_error(xparam);
    //cout<<"    assumed error : "<<assumederror<<" K"<<endl;
    //double finalerror = assumederror;
    //cout<<"Final error used in the residual calculation: "<<finalerror<<endl;
    //double prob = exp(residual_weight*(-residual)/finalerror);
    //double lnprob = log(prob*prior(xparam));
    //return lnprob;
    /*
    //evaluating proper values for xdelta of each parameter
    cout<<endl<<"Evaluating proper values for xdelta"<<endl;
    VecDoub xparam_bestmodel = xparam;
    const double bestfit_residual = residual;
    double acceptance_old;
    for(int i=0;i<nfreeparam;i++)
    {
        cout<<"Working on parameter "<<i+1<<endl;
        double xparam_ini = xparam[i];
        double acceptance = 1;
        do
        {
            acceptance_old = acceptance;
            xparam[i] += xdelta[i];
            cout<<"    testing "<<xparam_ini<<" -> "<<xparam[i]<<endl;
            residual = tominimize(xparam);
            prob_ratio = exp(residual_weight*(bestfit_residual - residual)/finalerror);
            acceptance = min(1.,prob_ratio*prior(xparam)/prior(xparam_bestmodel));
            cout<<"    acceptance probability = "<<acceptance*100<<" percent"<<endl;
        }
        while (acceptance>0.68);
        double x1 =  xparam[i]-xdelta[i];
        double x2 =     xparam[i];
        double y1 =  acceptance_old;
        double y2 =  acceptance;
        xdelta[i] = fabs(x1+(0.68-y1)/(y2-y1)*(x2-x1)-xparam_ini);
        cout<<"   xdelta for this parameter fixed at "<<xdelta[i]<<endl;
        xparam[i] = xparam_ini;
    }
    */
    //the cycle begins!
    ofstream results("mcmc_results.txt");
    results<<"0";
    for(int i=0;i<nfreeparam;i++) results<<"  "<<xparam[i];
    results<<endl;
    int newper, oldper = 0;
    for(int n=0;n<ncycle;n++)
    {
        newper = 100*(double(n)/double(ncycle));
        if(newper>oldper)
        {
            cout<<"   progress = "<<newper<<" percent"<<endl;
            oldper = newper;
        }
        
        //copy the values of the previous iteration
        xparam_old   = xparam;
        residual_old = residual;
        
        //try to upgrade the point using the proposal function (multivariate gaussian)
        for(int i=0;i<nfreeparam;i++)
        {
            xparam[i] = xparam_old[i] + ran_gau(seed,xdelta[i]);
        }
        
        cout<<endl<<"ITERATION NUMBER "<<n+1<<endl;
        cout<<"New proposed model:"<<endl;
        for(int i=0;i<nfreeparam;i++) cout<<xparam_old[i]<<" -> "<<xparam[i]<<endl;
        
        if(prior(xparam)>0)
        {
            residual = tominimize(xparam);
            
            //acceptance probability
            cout<<"   Residual (previous model) : "<<residual_old<<endl;
            cout<<"   Prior    (previous model) : "<<prior(xparam_old)<<endl;
            cout<<"   Residual (this model)     : "<<residual<<endl;
            cout<<"   Prior    (this model)     : "<<prior(xparam)<<endl;
            prob_ratio = exp(residual_weight*(residual_old - residual)/finalerror);
            acceptance = min(1.,prob_ratio*prior(xparam)/prior(xparam_old));
        }
        else
        {
            cout<<"   Prior is 0 here!"<<endl;
            acceptance = 0;
        }
        
        //if the new point is rejected, just stick on the old one
        cout<<"Acceptance probability: "<<acceptance*100<<" percent";
        if(ran0(&seed)>acceptance)
        {
            cout<<" (REJECTED)"<<endl;
            xparam   = xparam_old;
            residual = residual_old;
        }
        else {cout<<" (ACCEPTED)"<<endl;}
        for(int i=0;i<nfreeparam;i++) mcmc_all(i,n+1) = xparam[i];
        
        //trick: adjust the proposal function for the next iteration, according to how "acceptance" behaves
        if(acceptance<0.1) for(int i=0;i<nfreeparam;i++) xdelta[i]*=0.75;
        if(acceptance>0.4) for(int i=0;i<nfreeparam;i++) xdelta[i]*=1.5;
        
        //print results
        results<<n+1;
        for(int i=0;i<nfreeparam;i++) results<<"  "<<xparam[i];
        results<<endl;
        
        //evaluate average and stdv of each variable at each step
        for(int i=0;i<nfreeparam;i++)
        {
            vector<double> temp;
            for(int j=0;j<ncycle+1;j++) if(!isnan(mcmc_all(i,j))) temp.push_back(mcmc_all(i,j));
            cout<<"    parameter ["<<i+1<<"] = "<<myutil::average(temp)<<" +/- "<<myutil::stdev(temp,myutil::average(temp))<<endl;
        }
    }
    results.close();
}

void extragas_fit::bestmodel_finder(char* param_file)
{
	p.resize(nparam); dp.resize(nparam);
	
	//cout<<endl<<"ITERATIVE SEARCH FOR THE BEST MODEL (DOWNHILL ALGORITHM)"<<endl<<endl;
	//cout<<endl<<"Reading the parameter file "<<param_file<<endl;
	findValue(param_file, p[0], "vkick");
	findValue(param_file, p[1], "fion");
	findValue(param_file, p[2], "alpha_accr");
	
	findValue(param_file, dp[0], "vkick_delta");
	findValue(param_file, dp[1], "fion_delta");
	findValue(param_file, dp[2], "alpha_accr_delta");
	
	findValue(param_file, residual_above, "residual_above");
	findValue(param_file, ftol, "ftol");
	findValue(param_file, write_cube, "write_cube");
	findValue(param_file, use_blanks, "use_blanks");
	
	//number of free parameters
	int nfreeparam = 0;
	for(int i=0; i<nparam; i++) if(dp[i]!=0) nfreeparam+=1;
	if(nfreeparam<=0)
	{
		cout<<"No free parameters. No minimization is required"<<endl;
		return;
	}
	else {cout<<"NUMBER OF FREE PARAMETERS: "<<nfreeparam<<endl;}
	
	//initialise the initial points
	VecDoub xparam(nfreeparam);
	VecDoub xdelta(nfreeparam);
	int j=0;
	for(int i=0;i<nparam;i++)
	{
		if(dp[i]!=0) 
		{
			xparam[j]=p[i];
			xdelta[j]=dp[i];
			j++;
		}
	}
	VecDoub best_param;	 Doub res_value;
	//minimization using the downhill algorithm
	best_param = downhill(xparam,xdelta);
	res_value  = tominimize(best_param); 
	
	cout<<endl;
	cout<<"---------------------------------------------------"<<endl;
	cout<<"THE ALGORITHM CONVERGED TO THE FOLLOWING BEST MODEL"<<endl;
	if(dp[0]!=0)  cout<<"vkick			[km/s]  = "<<p[0]<<endl;
	if(dp[1]!=0)  cout<<"fion			        = "<<p[1]<<endl;
	if(dp[2]!=0)  cout<<"alpha_accr	    [1/Gyr] = "<<p[2]<<endl;
	cout<<"HI halo mass   [Mo]    = "<<mass_HIhalo<<endl;
	cout<<"accretion rate [Mo/yr] = "<<accretionrate<<endl;
	
	cout<<endl;
	cout<<"BEST RESIDUAL = "<<res_value<<endl;
}


void extragas_fit::logfile_update(char* file, int step, double myresidual)
{
	ofstream logfile(file, ios_base::app);
	logfile<<"STEP "<<step<<endl<<endl;
	if(dp[0]!=0)  logfile<<"vkick			[km/s]  = "<<p[0]<<endl;
	if(dp[1]!=0)  logfile<<"fion			        = "<<p[1]<<endl;
	if(dp[2]!=0)  logfile<<"alpha_accr	    [1/Gyr] = "<<p[2]<<endl;
	logfile<<"HI halo mass   [Mo]    = "<<mass_HIhalo<<endl;
	logfile<<"accretion rate [Mo/yr] = "<<accretionrate<<endl;
	logfile<<endl;
	logfile<<"RESIDUAL = "<<myresidual<<endl;
	logfile<<"------------------"<<endl<<endl;
	logfile.close();
}


///////////////////////////////////////////////////////////////
//STUFF FOR THE RESIDUAL MINIMIZATION VIA THE DOWNHILL METHOD//
///////////////////////////////////////////////////////////////

inline void extragas_fit::get_psum(MatDoub_I &p, VecDoub_O &psum)
{
	for (Int j=0;j<ndim;j++) 
	{
		Doub sum=0.0;
		for (Int i=0;i<mpts;i++)
			sum += p[i][j];
		psum[j]=sum;
	}
}

Doub extragas_fit::amotry(MatDoub_IO &p, VecDoub_O &y, VecDoub_IO &psum, const Int ihi, const Doub fac)
{
	VecDoub ptry(ndim);
	Doub fac1=(1.0-fac)/ndim;
	Doub fac2=fac1-fac;
	for (Int j=0;j<ndim;j++)
		ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
	Doub ytry=tominimize(ptry);
	if (ytry < y[ihi]) 
	{
		y[ihi]=ytry;
		for (Int j=0;j<ndim;j++) 
		{
			psum[j] += ptry[j]-p[ihi][j];
			p[ihi][j]=ptry[j];
		}
	}
	return ytry;
}

}
VecDoub extragas_fit::downhill(MatDoub_I &p0)
{
	char* logfile_name = "downhill.log";
	ofstream logfile(logfile_name);
	
	const Int NMAX=10000;
	const Doub TINY=1.0e-10;
	Int ihi,ilo,inhi;
	mpts=p0.nrows();
	ndim=p0.ncols();
	VecDoub psum(ndim),pmin(ndim),x(ndim);
	pnow=p0;
	y.resize(mpts);
	
	cout<<"Evaluating the initial "<<mpts<<" vertices"<<endl; 
	for (Int i=0;i<mpts;i++) 
	{
		for (Int j=0;j<ndim;j++)
			x[j]=pnow[i][j];
		y[i]=tominimize(x);
		cout<<"Vertex "<<i+1<<" completed. Here the residual is : "<<y[i]<<endl;
	}
	cout<<endl;
	nfunc=0;
	get_psum(pnow,psum);
	int step_number = 0;
	for (;;) 
	{
		step_number ++;
		cout<<"STEP "<<step_number<<" IN PROGRESS"<<endl;
		ilo=0;
		ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
		for (Int i=0;i<mpts;i++) 
		{
			if (y[i] <= y[ilo]) ilo=i;
			if (y[i] > y[ihi]) 
			{
				inhi=ihi;
				ihi=i;
			} else if (y[i] > y[inhi] && i != ihi) inhi=i;
		}
		Doub rtol=2.0*abs(y[ihi]-y[ilo])/(abs(y[ihi])+abs(y[ilo])+TINY);
		cout<<"		Residual = "<<y[ilo]<<". Ratio = "<<rtol<<"/"<<ftol<<endl;
		if (rtol < ftol) 
		{
			SWAP(y[0],y[ilo]);
			for (Int i=0;i<ndim;i++) 
			{
				SWAP(pnow[0][i],pnow[ilo][i]);
				pmin[i]=pnow[0][i];
			}
			fmin=y[0];
			return pmin;
		}
		if (nfunc >= NMAX) 
		{	
			cout<<endl<<"WARNING!"<<endl;
			cout<<"Maximum value of fuction evaluations ("<<NMAX<<") has been reached"<<endl;
			cout<<"Returning the last best set of parameters"<<endl;
			SWAP(y[0],y[ilo]);
			for (Int i=0;i<ndim;i++) 
			{
				SWAP(pnow[0][i],pnow[ilo][i]);
				pmin[i]=pnow[0][i];
			}
			fmin=y[0];
			return pmin;
		}
		nfunc += 2;
		Doub ytry=amotry(pnow,y,psum,ihi,-1.0);
		if (ytry <= y[ilo])
			ytry=amotry(pnow,y,psum,ihi,2.0);
		else if (ytry >= y[inhi]) 
		{
			Doub ysave=y[ihi];
			ytry=amotry(pnow,y,psum,ihi,0.5);
			if (ytry >= ysave) 
			{
				for (Int i=0;i<mpts;i++) 
				{
					if (i != ilo) 
					{
						for (Int j=0;j<ndim;j++)
							pnow[i][j]=psum[j]=0.5*(pnow[i][j]+pnow[ilo][j]);
						y[i]=tominimize(psum);
					}
				}
				nfunc += ndim;
				get_psum(pnow,psum);
			}
		} 
		else --nfunc;
		
		//print each passage in the logfile
		logfile_update(logfile_name,step_number,y[ilo]);
	}
}

VecDoub extragas_fit::downhill(VecDoub_I &point, VecDoub_I &dels)
{
	Int ndim=point.size();
	MatDoub p0(ndim+1,ndim);
	for (Int i=0;i<ndim+1;i++) {
		for (Int j=0;j<ndim;j++)
			p0[i][j]=point[j];
		if (i !=0 ) p0[i][i-1] += dels[i-1];
	}
	return downhill(p0);
}

VecDoub extragas_fit::downhill(VecDoub_I &point, const Doub del)
{
	VecDoub dels(point.size(),del);
	return downhill(point,dels);
}
extragas_fit::~extragas_fit(void)
{
    delete[] datacube_ptr; datacube_ptr=NULL;
}
