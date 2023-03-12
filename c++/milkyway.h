/*
 *  milkyway.h
 *  TestCfits
 *
 *  Created by Antonino Marasco on 13/11/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef MILKYWAY_H
#define MILKYWAY_H

using namespace std;

//classe per la creazione del modello di vel(l,b) in LSR per la nostra galassia
class milkyway
{
	public:
	milkyway (float max_radius, float sun_radius, float buldge_radius, float zin, float zfin, float velocity); //costruttore per modello a spessore variabile (Wakker)
	milkyway (float max_radius, float sun_radius, float Hsun, float RS_FLARE, float velocity);//costruttore modello di disco con flare esponenziale. Da usare anche per il modello di Kalberla, in taglio geometrico
	milkyway (float sun_radius, float velocity);//costruttore per modello di disco di Kalberla (2008) - taglio fisico. I parametri del disco sono impostati automaticamente
	void envelope_deHeij (float lon, float lat, float step); //calcola il massimo e il minimo dell'inviluppo V(l,b), i dati in input sono in gradi
	void envelope_Wakker (float lon, float lat, float step);
	void envelope_Kalberla_geom (float lon, float lat, float step, int number_of_scaleheight); //valido per il disco di Kalberla (2008), la linea di vista si interrompe ad un taglio geometrico pari a number_of_scaleheight*altezzascala
	void envelope_Kalberla_phys (float lon, float lat, float step, float Rmax, float zmax, float density_limit);//valido per il disco di Kalberla (2008), la linea di vista si interrompe quando la densità diventa inferiore a density limit (in cm^-3)
	void plot_n_disk(float z, float Rmin, float Rmax, float dR, char* namefile); //genera un file namefile per graifcare n_disk
	void plot_n_hscale(float Rmax, float dR, int n_height, char* namefile);
	void plot_isodensity(float Rmax, float dR, float zmax, float dz, float densvalue, float denstoll,char* namefile);
	float vel_min (); //ritorna il minimo dell'inviluppo V(l,b)
	float vel_max ();//ritorna il massimo dell'inviluppo V(l,b)
//intervalli da utilizzare:
//logitudine: [0,360]
//latitudine: [-90,90]	
	private:
	float r_max, z1, z2, r_buldge, r_sun, v_sun, vmin, vmax;
	float hsun, Rs_flare, Rs_disk, n_disk_sun;//parametri validi per il disco di Kalberla (2008)
	float l, b; 
	
//funzioni private
	float v (float radius);
	float zmax_Wakker (float radius);
	float zmax_deHeij (float radius, float sigma);
	//funzioni che gestiscono il disco di Kalberla (2008)
	float zmax_Kalberla(float radius);
	float hs_disk(float radius);
	float midplane_disk(float radius);
	float n_disk(float radius, float height);
	
};
#endif