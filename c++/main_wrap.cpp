//
//  pythonconvert.cpp
//  Extragas
//
//  Created by Anqi Li on 26/08/2021.
//


#include <iostream>
#include <vector>
#include "extragas_fit.h"
#include "myutilities.h"
#include "array2d.h"
using namespace std;
bool verbose = true;
char inputfile[] = "/net/cuby3.apertif/data/users/li/extragas_grid/files/extragas_n2403.in";
extern "C"
{
    double grid(float v_k, float h, float accre,char* outfits)
{
        char* datacube_file  = "/net/cuby3.apertif/data/users/li/extragas/cubes/NGC2403.fits";
        char* maskcube_file  = "/net/cuby3.apertif/data/users/li/extragas/cubes/NGC2403_mask.fits";
          //list of free parameters and their range
        char* parameter_file = "/net/cuby3.apertif/data/users/li/extragas/files/free_parameters.txt";
          //file containing the mcmc chain
        char* chain_file = "/net/cuby3.apertif/data/users/li/extragas/files/mcmc_results.txt";
        //initializing object
        extragas_fit paperino(datacube_file,maskcube_file);
        //produce a single modelcube using the parameters of parameter_file
        double result=0;
        result=paperino.modelmaker(parameter_file, true, true,outfits,v_k,h,accre);
        try{result=paperino.modelmaker(parameter_file, true, true,outfits,v_k,h,accre);}
        catch( const std::invalid_argument){
        //catch(std::exception& e) {
        return 0;
     }
        //stdout_on();
        return result;
        

}
}
