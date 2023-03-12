/*
 *  milkyway.cpp
 *  TestCfits
 *
 *  Created by Antonino Marasco on 20/11/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
#include <iostream>
#include <fstream>
#include <cmath>
#include "milkyway.h"
#define DEGTORAD 0.01745329252
using namespace std;

milkyway::milkyway (float max_radius, float sun_radius, float buldge_radius, float zin, float zfin, float velocity)
{
	r_max = max_radius; z1 = zin; z2 = zfin; r_sun = sun_radius; v_sun = velocity; r_buldge = buldge_radius;
};

milkyway::milkyway (float max_radius, float sun_radius, float Hsun, float RS_FLARE, float velocity)
{
	r_max = max_radius; r_sun = sun_radius; hsun = Hsun; Rs_flare = RS_FLARE; v_sun = velocity;
};

milkyway::milkyway (float sun_radius, float velocity)//Fissa i parametri a quelli di Kalberla
{
	r_sun = sun_radius;
	v_sun = velocity;
	hsun = 0.15;
	Rs_flare = 9.8;
	Rs_disk = 3.15;
	n_disk_sun = 0.9;
};

//dal modello di Wakker 91
inline float milkyway::v (float radius)
{
	return v_sun;
	//if(radius<0.5) {return (v_sun*radius/0.5);}
	//else {return v_sun;}
};
inline float milkyway::hs_disk(float radius)
{	
	if(radius>5.) {return (hsun*exp((radius-r_sun)/Rs_flare));}
	else {return (hsun*exp((5.-r_sun)/Rs_flare));}
	//return (hsun*exp((radius-r_sun)/Rs_flare));
}
inline float milkyway::midplane_disk(float radius)
{
	if (radius>7.) {return (n_disk_sun*exp(-(radius-r_sun)/Rs_disk));}
	else {return (n_disk_sun*exp(-(7.-r_sun)/Rs_disk));}
	//return (n_disk_sun*exp(-(radius-r_sun)/Rs_disk));
};
inline float milkyway::n_disk(float radius, float height)
{		
//exponential
	//return (midplane_disk(radius)*exp(-abs(height)/hs_disk(radius)));
//gaussian
	return (midplane_disk(radius)*exp(-log(2.)*height*height/(hs_disk(radius)*hs_disk(radius))));
//sech^2
	//return (midplane_disk(radius)*exp(-abs(height)/hs_disk(radius)));
};
void milkyway::plot_n_disk(float z, float Rmin, float Rmax, float dR, char* namefile)
{
	ofstream data;
	data.open(namefile);
	float R=Rmin;
	while(R<=Rmax)
	{
		data<<R<<" "<<n_disk(R,z)<<endl;
		R += dR;
	}
	data.close();
};
void milkyway::plot_n_hscale(float Rmax, float dR, int n_height, char* namefile)
{
	ofstream data;
	data.open(namefile);
	float R=0;
	float hs;
	while(R<=Rmax)
	{
		hs=hs_disk(R);
		data<<R<<" "<<n_height*hs<<endl;
		data<<R<<" "<<-n_height*hs<<endl;
		data<<-R<<" "<<n_height*hs<<endl;
		data<<-R<<" "<<-n_height*hs<<endl;
		R += dR;
	}
	data.close();
};
void milkyway::plot_isodensity(float Rmax, float dR, float zmax, float dz, float densvalue, float denstoll,char* namefile)
{
	ofstream data;
	data.open(namefile);
	data<<8.5<<" "<<0<<endl;
	float R,z;
	float density;
	z=0;
	while(z<=zmax)
	{
		R = 0.;
		while(R<=Rmax)
		{
			density = n_disk(R,z);
			if(density<=(densvalue+denstoll) && density>=(densvalue-denstoll))
			{
				data<<R<<" "<<z<<endl;
				data<<-R<<" "<<z<<endl;
				data<<R<<" "<<-z<<endl;
				data<<-R<<" "<<-z<<endl;
			}
			R+=dR;
		}
		z+=dz;
	}
	data.close();
}

//da Kalberla 2008
inline float milkyway::zmax_Kalberla(float radius)
{
	if(radius>5.) {return (hsun*exp((radius-r_sun)/Rs_flare));}
	else {return (hsun*exp((5.-r_sun)/Rs_flare));}
};
//Da Wakker, High Velocity Clouds (book)
inline float milkyway::zmax_Wakker (float radius)
{
	if (radius<=r_sun) {return z1;}
	else {return (z1 + (z2-z1)*(radius/r_sun - 1.)*(radius/r_sun - 1.)/4.);}
};
void milkyway::envelope_Kalberla_phys (float lon, float lat, float step, float Rmax, float zmax, float density_limit)
{
	l = lon*DEGTORAD;
	b =	lat*DEGTORAD;
	vmin = 1e20; vmax = -1e20;
	float d = 0, R = r_sun, z = 0, vlsr = 0;
	while (R<=Rmax && abs(z)<=zmax)
	{
		if (n_disk(R,z)>=density_limit)
		{
			vlsr = (r_sun*v(R)/R - v_sun)*sin(l)*cos(b);
			if (vlsr>vmax) vmax = vlsr;
			if (vlsr<vmin) vmin = vlsr;
		}
		d+=step;
		z = d*sin(b);
		R = r_sun*sqrt(pow(cos(b)*d/r_sun,2) - 2.*cos(b)*cos(l)*d/r_sun + 1.);
	}
};

void milkyway::envelope_Kalberla_geom (float lon, float lat, float step, int number_of_scaleheight)
{
	l = lon*DEGTORAD;
	b =	lat*DEGTORAD;
	float Hmax;
	vmin = 1e20; vmax = -1e20;
	float d = 0, R = r_sun, z = 0, vlsr = 0;
	while (R<=r_max && abs(z)<=r_max/2)
	{
		Hmax = zmax_Kalberla(R) * number_of_scaleheight; 
		if (abs(z)<=Hmax)
		{
			vlsr = (r_sun*v(R)/R - v_sun)*sin(l)*cos(b);
			if (vlsr>vmax) vmax = vlsr;
			if (vlsr<vmin) vmin = vlsr;
		}
		d+=step;
		z = d*sin(b);
		R = r_sun*sqrt(pow(cos(b)*d/r_sun,2) - 2.*cos(b)*cos(l)*d/r_sun + 1.);
	}
};

void milkyway::envelope_Wakker (float lon, float lat, float step)
{
	l = lon*DEGTORAD;
	b =	lat*DEGTORAD;
	vmin = 0.; vmax = 0.;
	float d = 0, R = 0, z = 0, vlsr = 0;
	while (R<=r_max && abs(z)<=r_max/2)
	{
		z = d*sin(b);
		R = r_sun*sqrt(pow(cos(b)*d/r_sun,2) - 2.*cos(b)*cos(l)*d/r_sun + 1.);
		if (abs(z)<=zmax_Wakker(R))
		{
			vlsr = (r_sun*v(R)/R - v_sun)*sin(l)*cos(b);
			if (vlsr>vmax) vmax = vlsr;
			if (vlsr<vmin) vmin = vlsr;
			d+=step;
		}
		else 
		{
			d+=step;
			continue;
		}
	}
};


//Da de Heij et al. 2002
float milkyway::zmax_deHeij (float radius, float sigma)
{
//standard
	if (radius<=11.5) {return (sigma*z1);}
	else {return (sigma*(z1+(radius-11.5)*z2));}
//nostro
	//if (radius<=11.5) {return z1;}
	//else {return (z1+sigma*(radius-11.5)*z2);}
};
void milkyway::envelope_deHeij (float lon, float lat, float step)
{
	l = lon*DEGTORAD;
	b =	lat*DEGTORAD;
	vmin = 0.; vmax = 0.;
	
	float mysigma = 5.;
	float alpha = 67;
	alpha*=DEGTORAD;
	
	float d = 0, R = 0, z = 0, vlsr = 0;
	float sin_phi, cos_phi, mysin, mycos, z0;
	while (R<=r_max && abs(z)<=15.)
	{
		z = d*sin(b);
		R = r_sun*sqrt(pow(cos(b)*d/r_sun,2) - 2.*cos(b)*cos(l)*d/r_sun + 1.);
		if (R<=11.5) {z0 = 0.;}
		else
		{
			sin_phi = (d/R)*sin(l)*cos(b);
			cos_phi = -(-d*cos(b)*cos(l)+r_sun)/R;
			mysin = sin_phi*cos(alpha)+cos_phi*sin(alpha);
			mycos = cos_phi*cos(alpha)-sin_phi*sin(alpha);
			z0 = (R-11.5)*mysin/6. + 0.3*(1.-2.*mycos)*pow((R-11.5)/6.,2);
		}
		if (z<=(z0+zmax_deHeij(R,mysigma)/2.) && z>=(z0-zmax_deHeij(R,mysigma)/2.))
		{
			vlsr = (r_sun*v(R)/R - v_sun)*sin(l)*cos(b);
			if (vlsr>vmax) vmax = vlsr;
			if (vlsr<vmin) vmin = vlsr;
			d+=step;
		}
		else
		{
			d+=step;
			continue;
		}
	}
};

float milkyway::vel_min () {return vmin;};
float milkyway::vel_max () {return vmax;};
