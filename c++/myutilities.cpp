/*
 *  myfits_utilities.cpp
 *  TestCfits
 *
 *  Created by Antonino Marasco on 14/04/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */


#include "myutilities.h"
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <sstream>
#include "fitsio.h"
#include "fitsutil.h"
#include "statist.h"
#include "array2d.h"
#define PI 3.14159265359 

using namespace std;

//long int seed = time(NULL);

inline double myutil::sgn (double value)
{
	if (value>0.){return 1.;}
	else if (value<0.){return -1;}
	else {return 0.;}
};

double myutil::max(std::vector<double> v)
{
	double mymax = -1e20;
	if(v.size()==1) {return mymax;}
	else
	{
		for(int i=0; i<v.size(); i++) if(v[i]>mymax) {mymax=v[i];}
	}
	return mymax;
};

double myutil::max(std::vector<double> x, std::vector<double> y)
{
	if(x.size()==0) return NAN;
	if(x.size()!=y.size())
	{
		cout<<"Error in function myutil::max. Vectors x and y have different sizes"<<endl;
		return NAN;
	}
	double xmax = x[0]; double ymax = -1e30;
	if(y.size()==1) {return xmax;}
	else
	{
		for(int i=0; i<y.size(); i++) if(y[i]>ymax) 
		{
			ymax = y[i];
			xmax = x[i];
		}
	}
	return xmax;
};

double myutil::min(std::vector<double> v)
{
	double mymin = 1e20;
	if(v.size()==0) return NAN;
	else
	{
		for(int i=0; i<v.size(); i++) if(v[i]<mymin) {mymin=v[i];}
	}
	return mymin;
};

double myutil::max(myarray::double3d cube)
{
	if(cube.xsize()*cube.ysize()*cube.zsize()<=0) return NAN;
	double mymax = cube(0,0,0);
	if(cube.xsize()*cube.ysize()*cube.zsize()==1) {return mymax;}
	else
	{
		for(int i=0;i<cube.xsize();i++) for(int j=0; j<cube.ysize();j++) for(int k=0;k<cube.zsize();k++)
		{
			if(cube(i,j,k)>mymax) mymax=cube(i,j,k);
		}
	}
	return mymax;
};

double myutil::min(myarray::double3d cube)
{
	if(cube.xsize()*cube.ysize()*cube.zsize()<=0) return NAN;
	double mymin = cube(0,0,0);
	if(cube.xsize()*cube.ysize()*cube.zsize()==1) {return mymin;}
	else
	{
		for(int i=0;i<cube.xsize();i++) for(int j=0; j<cube.ysize();j++) for(int k=0;k<cube.zsize();k++)
		{
			if(cube(i,j,k)<mymin) mymin=cube(i,j,k);
		}
	}
	return mymin;
};

int myutil::nlines_in(std::string file)
{
	ifstream data(file.c_str());
	char* buffer;
	int number_of_raws=0, lenght;
	data.seekg(0, ios::end);
	lenght = data.tellg();
	data.seekg(0, ios::beg);
	buffer = new char [lenght];
	data.read (buffer, lenght);
	for(int i=0; i<lenght; i++) if(buffer[i]==0x0A) number_of_raws++;
	data.seekg(0, ios::beg);
	delete[] buffer;
	data.close();
	return number_of_raws;
}

int myutil::nlines_in(char* file)
{
	ifstream data(file);
	char* buffer;
	int number_of_raws=0, lenght;
	data.seekg(0, ios::end);
	lenght = data.tellg();
	data.seekg(0, ios::beg);
	buffer = new char [lenght];
	data.read (buffer, lenght);
	for(int i=0; i<lenght; i++) if(buffer[i]==0x0A) number_of_raws++;
	data.seekg(0, ios::beg);
	delete[] buffer;
	data.close();
	return number_of_raws;
}

double myutil::interpol_lin(double X, std::vector<double> x, std::vector<double> y)
{
	if(X<=x[0]) {return y[0];}
	else if(X>=x[x.size()-1]) {return y[y.size()-1];}
	else
	{
		int i = (X-x[0])/(x[1]-x[0]);
		return y[i] + (y[i+1]-y[i])*(X-x[i])/(x[i+1]-x[i]);
	}
}

double bisection (double(*f)(double), double value, double xmin, double xmax, double precision,int n_iterations)
{
	double x1, x2, xtemp, v1, v2, vtemp;
	int counter = 0;
	x1 = xmin; x2 = xmax; 
	v1 = f(x1)-value;
	v2 = f(x2)-value;
	if(myutil::sgn(v1) == myutil::sgn(v2))
	{
		cout<<"ERROR IN BISECTION FUNCTION. F(X1) and F(X2) HAVE THE SAME SIGN!"<<endl;
		return x1;
	}
	else
	{	
		while(counter<n_iterations)
		{
			xtemp = (x2+x1)/2.;
			vtemp = f(xtemp)-value;
			if (abs(vtemp)<=precision) {return xtemp;}
			else
			{
				if(myutil::sgn(vtemp)==myutil::sgn(f(x1)-value)) {x1=xtemp;}
				else {x2=xtemp;}
				counter += 1;
			}
		}
		cout<<"MAX NUMBER OF ITERATION IN BISECTION FCT"<<endl;
		return xtemp;
	}
};

double myutil::gaussian (double x, double Sigma)
{
	return (1./(Sigma*sqrt(2.*PI)) * exp (-x*x/(2.*Sigma*Sigma)));
};

double myutil::aitoffX (double lambda, double phi)//in radianti!
{
	double alpha = acos(cos(phi)*cos(lambda/2.));
	return (2.*cos(phi)*sin(lambda/2.)/(sin(alpha)/alpha));
};

double myutil::aitoffY (double lambda, double phi)//in radianti!
{
	double alpha = acos(cos(phi)*cos(lambda/2.));
	return (sin(phi)/(sin(alpha)/alpha));
};

double myutil::hammerX (double lambda, double phi)//in radianti!
{
	return (2.*sqrt(2)*cos(phi)*sin(lambda/2.)/sqrt(1.+cos(phi)*cos(lambda/2.)));
};

double myutil::hammerY (double lambda, double phi)//in radianti!
{
	return (sqrt(2)*sin(phi)/sqrt(1.+cos(phi)*cos(lambda/2.)));
};
 
//calcola l'integrale col metodo del trapezio
double myutil::integral(double(*f)(double x), double x_ini, double x_fin, double dx)
{	
	double sum = 0.;
	double x1, x2;
	x1 = x_ini;
	x2 = x1+dx;
	while (x2<=x_fin)
	{
		sum += (f(x1)+f(x2))*dx/2.;
		x1=x2;
		x2 += dx;
	}
	return sum;	
};

double myutil::average(vector<double> x)
{
	double sum = 0;
	for(int i=0;i<x.size();i++) sum+=x[i];
	return sum/x.size();
}

double myutil::stdev(vector<double> x, double M)
{
	vector<double> newx(x.size());
	for(int i=0;i<x.size();i++) newx[i]=(x[i]-M)*(x[i]-M);
	return sqrt(myutil::average(newx));
}

double myutil::median(vector<double> x)
{
	sort(x.begin(),x.end());
	int dim = x.size();
	if(dim%2==1) {return x[(dim-1)/2];}
	else {return 0.5*(x[dim/2]+x[dim/2-1]);}
}

double myutil::median(vector<double> x, vector<double> y)
{
	double integral = 0;
	for(int i=0;i<x.size()-1;i++) integral += 0.5*(x[i+1]-x[i])*(y[i]+y[i+1]); 
	double intH = 0.5*integral;
	double int1 = 0, int2 = 0;
	for(int i=0;i<x.size()-1;i++)
	{
		int2 += 0.5*(x[i+1]-x[i])*(y[i]+y[i+1]);
		if(int2>=intH) 
		{
			return x[i]+(intH-int1)*(x[i+1]-x[i])/(int2-int1);//linear interpolation
		}
		int1 = int2;
	}
	return NAN;
}

double myutil::closest_peak(std::vector<double> &x, std::vector<double> &y, double value, double threshold_fraction)
{
	//normalize the profile
	double ymax = myutil::max(y);
	
	int peak_index = -10;
	double peak_diff = 1e30;
	for(int n=2; n<x.size()-2;n++)
	{
		//if(y[n-2]<y[n-1]<y[n]>y[n+1]>y[n+2] && y[n]>threshold_fraction*ymax && fabs(x[n]-value)<=peak_diff) //sure about this?
		if(y[n-2]<y[n-1] && y[n-1]<=y[n] && y[n]>=y[n+1] && y[n+1]>y[n+2] && y[n]>threshold_fraction*ymax && fabs(x[n]-value)<=peak_diff)
		{
			peak_index = n;
			peak_diff = fabs(x[peak_index]-value);
		}
	}
	if(peak_index<0) {return NAN;}
	else {return x[peak_index];}
}

double myutil::mad(vector<double> x, vector<double> y, double M)
{
	vector<double> tempx = x; 
	for(int i=0;i<tempx.size();i++) tempx[i] = fabs(x[i]-M);
	double newx_max = std::max(tempx[0],tempx[tempx.size()-1]);
	double dx = fabs(x[1]-x[0]);
	int dim = newx_max/(x[1]-x[0])+1;
	vector<double> newx(dim), newy(dim);
	for(int i=0;i<dim;i++) newx[i] = i*dx;
	for(int i=0;i<dim;i++)
	{
		for(int n=0;n<tempx.size();n++) if(tempx[n]>=newx[i]-0.5*dx && tempx[n]<newx[i]+0.5*dx) newy[i]+=y[n];
	}
	newy[0]*=2.;
	return myutil::median(newx,newy);
}

double myutil::gmad(vector<double> x, vector<double> y, double M)
{
	vector<double> tempx = x; 
	for(int i=0;i<tempx.size();i++) tempx[i] = fabs(x[i]-M);
	double newx_max = std::max(tempx[0],tempx[tempx.size()-1]);
	double dx = fabs(x[1]-x[0]);
	int dim = newx_max/(x[1]-x[0])+1;
	vector<double> newx(dim), newy(dim);
	for(int i=0;i<dim;i++) newx[i] = i*dx;
	for(int i=0;i<dim;i++)
	{
		for(int n=0;n<tempx.size();n++) if(tempx[n]>=newx[i]-0.5*dx && tempx[n]<newx[i]+0.5*dx) newy[i]+=y[n];
	}
	newy[0]*=2.;
	return 1.4826*myutil::median(newx,newy);
}

string myutil::convertInt(int number)
{
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

