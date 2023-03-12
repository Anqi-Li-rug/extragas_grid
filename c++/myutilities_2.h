/*
 *  utilities.h
 *  TestCfits
 *
 *  Created by Antonino Marasco on 27/11/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef MYUTILITIES_H
#define MYUTILITIES_H
#include <vector>
#include <cmath>
#include <string>

using namespace std;

namespace myutil
{
//convert an int into a string
	string convertInt(int number);
//restituisce il segno di value
	double sgn (double value);
//trova il massimo elemento del vettore v
	double max (std::vector<double> v);
//trova il minimo elemento del vettore v
	double min (std::vector<double> v);
//linear interpolation of the binned y(x) at position X	
	double interpol_lin(double X, std::vector<double> x, std::vector<double> y);
//return the number of lines inside the file	
	int nlines_in(std::string file);
	int nlines_in(char* file);
//trova la x tale che f(x)=value, ricercandone il valore tra xmin e xmax, con precisione precisione, e facendo un numero massimo di step n_iterations
	double bisection (double(*f)(double), double value, double xmin, double xmax, double precision,int n_iterations);
//ritorna il valore della funzione normale (gaussiana normalizzata)	
	double gaussian (double x, double Sigma);
//proiezioni
	double aitoffX (double lambda, double phi);//in radianti!
	double aitoffY (double lambda, double phi);//in radianti!
	double hammerX (double lambda, double phi);//in radianti!
	double hammerY (double lambda, double phi);//in radianti!
//calcola l'integrale col metodo del trapezio
	double integral(double(*f)(double x), double x_ini, double x_fin, double dx);
//calcola media e deviazione standard	
	double average(std::vector<double> x);
	double stdev(std::vector<double> x, double M);
//evaluate median and the mad from a distribution y(x)
 	double median(std::vector<double> x, std::vector<double> y);
	double median(std::vector<double> x);
	double mad(std::vector<double> x, std::vector<double> y, double M);
	double gmad(std::vector<double> x, std::vector<double> y, double M);//the same as mad, but multiplied by 1.4826 to be consistent with 1sigma from a gaussian distribution
//DATACUBE MANIPULATION
//create a total map and a (median) velocity field map starting from a datacube
	void create_maps(char* inputcube, char* totalmap_file, char* velocityfield_file, double density_threshold);
//Point in polygon (PIP) routine: check if the point (testx,testy) is inside the polygon	
	int pnpoly(int nvert, double *vertx, double *verty, double testx, double testy);
}
#endif
 
