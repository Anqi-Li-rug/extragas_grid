/*
 *  profiles.h
 *  Extragas
 *
 *  Created by Antonino Marasco on 29/07/13.
 *  Copyright 2013 __MyCompanyName__. All rights reserved.
 *	A very simple class that generate a gnuplot macro for the line profiles
 *
*/


#ifndef PROFILES_H
#define PROFILES_H

#include <vector>
#include <fstream>
#include "myutilities.h"

using namespace std;

 class profiles_macro
 {
 public:
	 profiles_macro(string plotfile_name)
	 {
		 plotfile = plotfile_name;
		 datafile.open("MACRO_profiles.txt");
		 datafile<<"reset"<<endl;
		 datafile<<"set term pdf enhanced colour fsize 8"<<endl;
		 datafile<<"set output '"<<plotfile<<"'"<<endl;
		 datafile<<"set style data lines"<<endl;
		 datafile<<"set xrange [-450:450]"<<endl;
		 datafile<<"set yrange [-0.05:1.05]"<<endl;
		 datafile<<"set xlabel 'velocity [km/s]'"<<endl;
		 datafile<<"set ylabel ''"<<endl;
		 datafile<<"set parametric"<<endl;
		 datafile<<"set trange [-10:10]"<<endl;
	}
	 
	 void add_profile(string id, string ref, char type, std::vector<double> vel, std::vector<double> profile, double l, double b, double vobs, double vpeak, double vmax, double vmed)
	 {
		 //normalize profile
		 double mymax = myutil::max(profile);
		 for(int i=0;i<profile.size();i++) profile[i]/=mymax;
		 //write the profile
		 string format = ".txt";
		 string id_complete = id+format;
		 ofstream datatemp(id_complete.c_str());
		 for(int i=0;i<vel.size();i++) datatemp<<vel[i]<<" "<<profile[i]<<endl;
		 datatemp.close();
		 //label some info
		 datafile<<"set label '"<<id<<"' left font 'Helvetica,6' at graph 0.02, 0.95"<<endl;
		 datafile<<"set label '(l,b):("<<l<<","<<b<<")' left font 'Helvetica,6' at graph 0.02, 0.9"<<endl;
		 datafile<<"set label 'type : "<<type<<"' left font 'Helvetica,6' at graph 0.02, 0.85"<<endl;
		 datafile<<"set label 'ref  : "<<ref<<"' left font 'Helvetica,6' at graph 0.02, 0.8"<<endl;
		 //add the plot to the macro
		 datafile<<"vobs = "<<vobs<<endl;
		 datafile<<"vpeak = "<<vpeak<<endl;
		 datafile<<"vmax = "<<vmax<<endl;
		 datafile<<"vmed = "<<vmed<<endl;
		 datafile<<"plot '"<<id_complete<<"' lw 4 lc rgb 'black' title '', vobs,t lw 3 lc rgb 'red' title 'observed', vpeak,t lw 3 lc rgb '#00008B' title 'peak', vmax,t lw 3 lc rgb '#1E90FF' title 'max', vmed,t lw 3 lc rgb '#00FFFF' title 'median'"<<endl;
		 datafile<<"unset label"<<endl;
	 }
	 
	 void close()
	 {
		 datafile<<"unset output"<<endl;
		 datafile.close();
	 };
 
	private:
	 ofstream datafile;
	 string plotfile;
 };

#endif