/*
  Plotting routines for exgas

  This file contains the plotting routines used by exgas.
  These routines make use of the gnuplot library.

*/

#include "gnuplot.h"
#include "array.h"
//#include "extragas.h"
#include <cmath>

using namespace ffarray;

extern int RAsize, DECsize, VELsize;
extern const double CONV1;
extern const double RAD;


void plot0col()
{
  Gnuplot_Tmpfile gp;
  gp.begin();
  // gp.commandln("set terminal x11");
  gp.commandln("set terminal postscript eps color enhanced 14 solid");
  gp.commandln("set out 'extragas0col.ps'");
  // page 1
  // ********************
  // ** THE INPUT FILE FOR THESE PLOTS IS CALLED extragas_2403_orbitplots.in **
  // ********************
  /*
  gp.commandln("set multiplot");
  gp.commandln("set nokey");
  gp.commandln("set size 0.5,0.5");
  gp.commandln("set origin 0,0");
  gp.commandln("set xlabel 'R (kpc)'");
  gp.commandln("set ylabel 'z (kpc)'");
  gp.commandln("set xrange [4.89:5.14]");
  gp.commandln("set yrange [0:2]");
  gp.commandln("set label 'Fountain + Accretion' at 4.95, 1.8 font 'Helvetica,16'");
  gp.commandln("plot 'orbit.dat' using 2:4 with lines 3");
  gp.commandln("set origin 0.5,0");
  gp.commandln("set xlabel 't (Myr)'");
  gp.commandln("set ylabel 'v (km/s)'");
  gp.commandln("set xrange [0:120]");
  gp.commandln("set yrange [-60:140]");
  gp.commandln("unset label");
  gp.commandln("set label 'v_R' at 75, 5 font 'Helvetica,18'");
  gp.commandln("set label 'v_z' at 15, 60 font 'Helvetica,18'");
  gp.commandln("set label 'v_{/Symbol f}' at 20, 105 font 'Helvetica,18'");
  gp.commandln("plot 'orbit.dat' using 1:7 with lines, 'orbit.dat' using 1:5 with lines, 'orbit.dat' using 1:6 with lines");
  */
  // ACCRETION RATE
  gp.commandln("set auto y");
  gp.commandln("set multiplot");
  gp.commandln("set nokey");
  gp.commandln("set size 0.5,0.5");
  gp.commandln("set origin 0,0");
  gp.commandln("set xrange [0:120]");
  gp.commandln("set yrange [1:1.18]");
  gp.commandln("set xlabel 't (Myr)'");
  gp.commandln("set ylabel 'Mass of the cloud'");
  gp.commandln("plot 'orbit.dat' using 1:12 with lines 3");
  gp.commandln("set origin 0.5,0");
  gp.commandln("set yrange [0.85:1]");
  gp.commandln("set ylabel 'L_z'");
  gp.commandln("plot 'orbit.dat' using 1:18 with lines 3");
  gp.end();
}

void plot0()
{
  Gnuplot_Tmpfile gp;
  gp.begin();
  // gp.commandln("set terminal x11");
  gp.commandln("set terminal postscript eps enhanced 14");
  gp.commandln("set out 'extragas0.ps'");
  // page 1
  /*
  gp.commandln("set multiplot");
  gp.commandln("set nokey");
  gp.commandln("set size 0.5,0.5");
  gp.commandln("set origin 0,0");
  gp.commandln("set xlabel 'R (kpc)'");
  gp.commandln("set ylabel 'z (kpc)'");
  gp.commandln("set xrange [4.89:5.14]");
  gp.commandln("set yrange [0:2]");
  gp.commandln("set label 'Fountain + Accretion' at 4.95, 1.8 font 'Helvetica,16'");
  gp.commandln("plot 'orbit.dat' using 2:4 with lines");
  gp.commandln("set origin 0.5,0");
  gp.commandln("set xlabel 't (Myr)'");
  gp.commandln("set ylabel 'v (km/s)'");
  gp.commandln("set xrange [0:120]");
  gp.commandln("set yrange [-60:140]");
  gp.commandln("unset label");
  gp.commandln("set label 'v_R' at 75, 5 font 'Helvetica,18'");
  gp.commandln("set label 'v_z' at 15, 60 font 'Helvetica,18'");
  gp.commandln("set label 'v_{/Symbol f}' at 20, 105 font 'Helvetica,18'");
  gp.commandln("plot 'orbit.dat' using 1:7 with lines, 'orbit.dat' using 1:5 with lines, 'orbit.dat' using 1:6 with lines");
  */
  // ACCRETION RATE
  gp.commandln("set auto y");
  gp.commandln("set multiplot");
  gp.commandln("set nokey");
  gp.commandln("set size 0.5,0.5");
  gp.commandln("set origin 0,0");
  gp.commandln("set xrange [0:120]");
  gp.commandln("set yrange [1:1.18]");
  gp.commandln("set xlabel 't (Myr)'");
  gp.commandln("set ylabel 'Mass of the cloud'");
  gp.commandln("plot 'orbit.dat' using 1:12 with lines");
  gp.commandln("set origin 0.5,0");
  gp.commandln("set yrange [0.85:1]");
  gp.commandln("set ylabel 'L_z'");
  gp.commandln("plot 'orbit.dat' using 1:18 with lines");
  gp.end();
  gp.end();
}

void plot1(float accr_norm,char psfile[30])
{
  Gnuplot_Tmpfile gp;
  gp.begin();
  // gp.commandln("set terminal x11");
  // gp.commandln("set terminal postscript eps color");
  // page 1
  gp.commandln("set terminal postscript portrait enhanced 10");
  gp.commandln(psfile);
  gp.commandln("set multiplot");
  gp.commandln("set nokey");
  gp.commandln("set size 0.5,0.25");
  gp.commandln("set origin 0,0.75");
  gp.commandln("set xlabel 'R (kpc)'");
  gp.commandln("set ylabel 'z (kpc)'");
  gp.commandln("plot 'orbit.dat' using 2:4 with lines");
  gp.commandln("set origin 0.5,0.75");
  gp.commandln("set xlabel 't (Myr)'");
  gp.commandln("set ylabel 'v (km/s)'");
  gp.commandln("plot 'orbit.dat' using 1:7 with lines, 'orbit.dat' using 1:5 with lines, 'orbit.dat' using 1:6 with lines");
  // page 2
  gp.commandln("set multiplot");
  gp.commandln("unset label");
  gp.commandln("set nokey");
  gp.commandln("set auto");
  gp.commandln("set size 0.33,0.25");
  gp.commandln("set origin 0,0.75");
  gp.commandln("set xlabel 't (Myr)'");
  gp.commandln("set ylabel 'R (kpc)'");
  gp.commandln("plot 'orbit.dat' using 1:2 with lines");
  gp.commandln("set ylabel 'phi'");
  gp.commandln("set origin 0,0.5");
  gp.commandln("plot 'orbit.dat' using 1:3 with lines");
  gp.commandln("set origin 0,0.25");
  gp.commandln("set ylabel 'z (kpc)'");
  gp.commandln("plot 'orbit.dat' using 1:4 with lines");
  gp.commandln("set origin 0.33,0.75");
  gp.commandln("set ylabel 'v_R (km/s)'");
  gp.commandln("plot 'orbit.dat' using 1:5 with lines");
  gp.commandln("set origin 0.33,0.5");
  gp.commandln("set ylabel 'v_{phi} (km/s)'");
  gp.commandln("plot 'orbit.dat' using 1:6 with lines");
  gp.commandln("set origin 0.33,0.25");
  gp.commandln("set ylabel 'v_z (km/s)'");
  gp.commandln("plot 'orbit.dat' using 1:7 with lines");
  gp.commandln("set origin 0.66,0.75");
  gp.commandln("set xlabel 'R (kpc)'");
  gp.commandln("set ylabel 'z (kpc)'");
  gp.commandln("plot 'orbit.dat' using 2:4 with lines");
  gp.commandln("set xlabel 'R (kpc)'");
  gp.commandln("set ylabel 'v_R (km/s)'");
  gp.commandln("set origin 0.66,0.5");
  gp.commandln("plot 'orbit.dat' using 2:5 with lines");
  gp.commandln("set xlabel 'z (kpc)'");
  gp.commandln("set ylabel 'v_z (km/s)'");
  gp.commandln("set origin 0.66,0.25");
  gp.commandln("plot 'orbit.dat' using 4:7 with lines");
  //  gp.commandln("set yrange [0.98:1.02]");
  gp.commandln("set xlabel 't (Myr)'");
  gp.commandln("set ylabel 'E_{tot} (kpc/Myr)^2'");
  gp.commandln("set origin 0,0");
  //  gp.commandln("set yrange [0.99:1.01]");
  gp.commandln("set autosca");
  gp.commandln("plot 'orbit.dat' using 1:8 with lines");
  gp.commandln("set xlabel 't (Myr)'");
  gp.commandln("set ylabel 'L_{tot} (kpc^2/Myr)'");
  gp.commandln("set origin 0.33,0");
  gp.commandln("set auto y");
  //  gp.commandln("set yrange [0.8:1.2]");
  gp.commandln("plot 'orbit.dat' using 1:9 with lines");
  gp.commandln("set xlabel 't (Myr)'");
  gp.commandln("set ylabel 'L_z (kpc^2/Myr)'");
  gp.commandln("set origin 0.66,0");
  gp.commandln("set auto y");
  //  gp.commandln("set yrange [0.9999:1.0001]");
  gp.commandln("plot 'orbit.dat' using 1:10 with lines");
  // page 3
  gp.commandln("set auto y");
  gp.commandln("set multiplot");
  gp.commandln("set nokey");
  gp.commandln("set size 0.33,0.25");
  gp.commandln("set xlabel 't (Myr)'");
  gp.commandln("set ylabel 'Mass of cloud'");
  gp.commandln("set origin 0,0.75");
  gp.commandln("plot 'orbit.dat' using 1:12 with lines");

  gp.commandln("set yrange [10:-50]");
  gp.commandln("set xlabel 't (Myr)'");
  gp.commandln("set ylabel 'Gradient'");
  gp.commandln("set origin 0.33,0.75");
  gp.commandln("plot 'orbit.dat' using ($1):(($6-220.)/$4) with lines");

  gp.commandln("set auto y");
  gp.commandln("set xlabel 't (Myr)'");
  gp.commandln("set ylabel 'Energy'");
  gp.commandln("set origin 0,0.5");
  gp.commandln("plot 'orbit.dat' using 1:13 with lines");
  gp.commandln("set xlabel 't (Myr)'");
  gp.commandln("set ylabel 'L'");
  gp.commandln("set origin 0.33,0.5");
  gp.commandln("plot 'orbit.dat' using 1:14 with lines");
  gp.commandln("set xlabel 't (Myr)'");
  gp.commandln("set ylabel 'L_z'");
  gp.commandln("set origin 0.66,0.5");
  gp.commandln("plot 'orbit.dat' using 1:15 with lines");

  gp.commandln("set yrange [0.9:1.1]");
  gp.commandln("set xlabel 't (Myr)'");
  gp.commandln("set ylabel 'Energy'");
  gp.commandln("set origin 0,0.25");
  gp.commandln("plot 'orbit.dat' using 1:16 with lines");
  gp.commandln("set xlabel 't (Myr)'");
  gp.commandln("set ylabel 'L'");
  gp.commandln("set origin 0.33,0.25");
  gp.commandln("plot 'orbit.dat' using 1:17 with lines");
  gp.commandln("set xlabel 't (Myr)'");
  gp.commandln("set ylabel 'L_z'");
  gp.commandln("set origin 0.66,0.25");
  gp.commandln("plot 'orbit.dat' using 1:18 with lines");

  gp.commandln("set xlabel 't (Myr)'");
  gp.commandln("set ylabel 'Energy'");
  gp.commandln("set origin 0,0");
  gp.commandln("plot 'orbit.dat' using 1:19 with lines");
  gp.commandln("set xlabel 't (Myr)'");
  gp.commandln("set ylabel 'L'");
  gp.commandln("set origin 0.33,0.");
  gp.commandln("plot 'orbit.dat' using 1:20 with lines");
  gp.commandln("set xlabel 't (Myr)'");
  gp.commandln("set ylabel 'L_z'");
  gp.commandln("set origin 0.66,0.");
  gp.commandln("plot 'orbit.dat' using 1:21 with lines");
  // page 4
  gp.commandln("set auto y");
  gp.commandln("set multiplot");
  gp.commandln("set nokey");
  gp.commandln("set size 0.5,0.25");
  gp.commandln("set origin 0,0.75");
  gp.commandln("set xlabel 't (Myr)'");
  gp.commandln("set ylabel 'Mass of the cloud'");
  gp.commandln("set origin 0,0.75");
  //  gp.commandln("plot 'orbit.dat' using 1:12 with lines, 1+x*0.001");
  //  gp.commandln("set origin 0,0.5");   
  //  gp.commandln("set log x");
  gp.commandln("plot 'orbit.dat' using 1:12 with lines");

  gp.commandln("unset log");
  gp.commandln("set origin 0.5,0.75");
  gp.commandln("set ylabel 'L_z'");
  gp.commandln("plot 'orbit.dat' using 1:18 with lines");
  
 // gp.commandln("set size 1,0.7");
  // gp.commandln("set origin 0,0");
  // gp.commandln("plot 'hvcs_out3.dat' using 4:5 2, 'hvcs_out2.dat' using 4:5 7 'hvcs_out4.dat' using 4:5 5,'hvcs_out1.dat' using 4:5 1, 'hvcs_out5.dat' using  4:5 3");
  // gp.commandln("plot [-180:180] [-90:90] 'hvcs_out3.dat' using 7:8 2, 'hvcs_out2.dat' using 7:8 6, 'hvcs_out4.dat' using 7:8 5,'hvcs_out1.dat' using 7:8 1, 'hvcs_out5.dat' using 7:8 3");
  // gp.commandln("splot 'hvcs_out.dat' using 4:5:6");
  // gp.commandln("pause -1");
  gp.end();
    cout<<"gnuppppppppplot"<<endl;

}

void plot2()
{ cout<<"plot2222222"<<endl;
  Gnuplot_Tmpfile gp;
  gp.begin();
  // panel 2
  gp.commandln("set terminal postscript portrait enhanced 16");
  gp.commandln("set out 'extragas2.ps'");
  gp.commandln("set multiplot");
  gp.commandln("set nokey");
  gp.commandln("set size 0.5,0.25");
  gp.commandln("set origin 0,.5");
  gp.commandln("set log y");
  gp.commandln("set xrange [0:30]");
  gp.commandln("set yrange [0.0001:.3]");
  //  gp.commandln("set yrange [0.00001:.03]");
  gp.commandln("set xlabel 'R (kpc)'");
  gp.commandln("set ylabel 'Density flow (M_o yr^{-1} kpc^{-2})'");
  gp.commandln("plot 'extragas_inout.dat' using 15:16 with boxes, 'extragas_inout.dat' using 15:17 with boxes");
  gp.commandln("set origin 0.5,0.5");
  gp.commandln("unset log y");
  gp.commandln("set autoscale");
  gp.commandln("set xrange [0:30]");
  gp.commandln("set yrange [0:35]");
  //  gp.commandln("set yrange [0:3.5]");
  gp.commandln("set xlabel 'R (kpc)'");
  gp.commandln("set ylabel 'Comunlative flows (M_o yr^{-1})'");
  gp.commandln("plot 'extragas_inout.dat' using 1:6 with l, 'extragas_inout.dat' using 1:7 with l");
  // panel 2
  gp.commandln("set multiplot");
  gp.commandln("set nokey");
  gp.commandln("set style line 1 lt 3 lw 2");
  gp.commandln("set style line 2 lt 1 lw 1");
  gp.commandln("set size 1,0.33");
  gp.commandln("set origin 0,.66");
  gp.commandln("set log y");
  gp.commandln("set xrange [0:30]");
  gp.commandln("set yrange [0.0001:1]");
  gp.commandln("unset xlabel");
  gp.commandln("set ylabel 'Density outflow (M_o yr^{-1} kpc^{-2})'");
  gp.commandln("plot 'extragas_inout.dat' using 1:2 with boxes ls 1, 'extragas_inout.dat' using 1:18 with boxes ls 2");
  gp.commandln("set origin 0,0.33");
  gp.commandln("set xrange [0:30]");
  gp.commandln("set yrange [0.0001:1]");
  gp.commandln("unset xlabel");
  gp.commandln("set ylabel 'Density inflow (M_o yr^{-1} kpc^{-2})'");
  gp.commandln("plot 'extragas_inout.dat' using 1:3 with boxes ls 1, 'extragas_inout.dat' using 1:19 with boxes ls 2");
  gp.commandln("set origin 0,0");
  gp.commandln("unset log y");
  gp.commandln("set xrange [0:30]");
  gp.commandln("set yrange [-2.5:2.5]");
  gp.commandln("set xlabel 'R (kpc)'");
  gp.commandln("set ylabel 'Net outflow (M_o yr^{-1})'");
  gp.commandln("plot 'extragas_inout.dat' using 1:8 with l ls 1, 'extragas_inout.dat' using 1:22 with l ls 2");
  // panel 3
  gp.commandln("set multiplot");
  gp.commandln("set nokey");
  gp.commandln("set size 0.5,0.25");
  gp.commandln("set origin 0,.5");
  gp.commandln("set autoscale");
  gp.commandln("set xrange [0:30]");
  gp.commandln("set xlabel 'R (kpc)'");
  gp.commandln("set ylabel 'Density flow (M_o/yr/kpc^2)'");
  gp.commandln("plot 'extragas_inout.dat' using 15:16 with boxes, 'extragas_inout.dat' using 15:17 with boxes");
  gp.commandln("set origin 0.5,0.5");
  gp.commandln("unset log y");
  gp.commandln("set autoscale");
  gp.commandln("set xrange [0:30]");
  gp.commandln("set xlabel 'R (kpc)'");
  gp.commandln("set ylabel 'Comunlative flows (M_o/yr)'");
  gp.commandln("plot 'extragas_inout.dat' using 1:6 with l, 'extragas_inout.dat' using 1:7 with l");
  // panel 4
  gp.commandln("set multiplot");
  gp.commandln("set nokey");
  gp.commandln("set size 0.5,0.25");
  gp.commandln("set origin 0,.75");
  gp.commandln("set xrange [0:30]");
  gp.commandln("set xlabel 'R (kpc)'");
  gp.commandln("set ylabel 'Outflow (M_O/yr/kpc^2)'");
  gp.commandln("plot 'extragas_inout.dat' using 1:2 with boxes");
  gp.commandln("set origin .5,.75");
  gp.commandln("set xlabel 'R (kpc)'");
  gp.commandln("set ylabel 'Inflow (M_O/yr/kpc^2)'");
  gp.commandln("plot 'extragas_inout.dat' using 1:3 with boxes");
  gp.commandln("set origin 0.5,.5");
  gp.commandln("set xlabel 'R (kpc)'");
  gp.commandln("set ylabel 'Total flow (M_O)'");
  gp.commandln("plot 'extragas_inout.dat' using 1:5 with boxes");
  gp.commandln("set origin 0,0.5");
  gp.commandln("set xlabel 'R (kpc)'");
  gp.commandln("set ylabel 'Comunlative flows (M_O/yr)'");
  gp.commandln("plot 'extragas_inout.dat' using 1:6 with l, 'extragas_inout.dat' using 1:7 with l");
  // panel 5
  gp.commandln("set multiplot");
  gp.commandln("set size 0.5,0.25");
  gp.commandln("set origin 0,0.75");
  gp.commandln("set xlabel 'R (kpc)'");
  gp.commandln("set ylabel 'Net flow (M_O/yr)'");
  gp.commandln("plot 'extragas_inout.dat' using 1:4 with boxes");
  gp.commandln("set origin 0.5,.75");
  gp.commandln("set xlabel 'R (kpc)'");
  gp.commandln("set ylabel 'Total outflow (M_O/yr)'");
  gp.commandln("plot 'extragas_inout.dat' using 1:8 with boxes");
  gp.commandln("set origin 0,0.5");
  gp.commandln("set xlabel 'R (kpc)'");
  gp.commandln("set ylabel 'Cumulative Energy outflow (erg/s)'");
  gp.commandln("plot 'extragas_inout.dat' using 1:10 with l");
  gp.commandln("set origin 0.5,0.5");
  gp.commandln("set xlabel 'R (kpc)'");
  gp.commandln("set ylabel 'Surface density (M_O/kpc^2)'");
  gp.commandln("plot 'extragas_inout.dat' using 1:12 w l, 'extragas_inout.dat' using 1:11 w l");
  gp.commandln("set origin 0,0.25");
  gp.commandln("set xlabel 'R (kpc)'");
  gp.commandln("set ylabel 'Surface density (M_O/kpc^2)'");
  gp.commandln("plot 'extragas_inout.dat' using 1:13 with l, 'extragas_inout.dat' using 1:11 with l");
  gp.commandln("set origin 0.5,0.25");
  gp.commandln("set xlabel 'R (kpc)'");
  gp.commandln("set ylabel 'Surface density (M_O/kpc^2)'");
  gp.commandln("plot 'extragas_inout.dat' using 1:14 with l, 'extragas_inout.dat' using 1:11 with l");
  // panel 6
  gp.commandln("set multiplot");
  gp.commandln("set nokey");
  gp.commandln("set style line 1 lt 3 lw 2");
  gp.commandln("set style line 2 lt 1 lw 1");
  gp.commandln("set size 0.5,0.25");
  gp.commandln("set origin 0,0.75");
  gp.commandln("set log y");
  gp.commandln("set xrange [0:30]");
  gp.commandln("set yrange [0.0001:1]");
  gp.commandln("set xlabel 'R (kpc)'");
  gp.commandln("set ylabel 'Density inflow (M_o yr^{-1} kpc^{-2})'");
  gp.commandln("plot 'extragas_inout.dat' using 1:23 with boxes ls 1");
  gp.commandln("set origin 0.5,0.75");
  gp.commandln("set xrange [0:30]");
  gp.commandln("set auto y");
  gp.commandln("set ylabel 'Density inflow (M_o yr^{-1} kpc^{-2})'");
  gp.commandln("plot 'extragas_inout.dat' using 15:24 with boxes ls 1");
  gp.commandln("set origin 0,0.5");
  gp.commandln("unset log y");
  gp.commandln("set autoscale y");
  gp.commandln("set xrange [0:30]");
  gp.commandln("set ylabel 'Comulative flow (M_o/yr)'");
  gp.commandln("plot 'extragas_inout.dat' using 1:25 with l");
  //  gp.commandln("set xlabel 'R (kpc)'");
  //  gp.commandln("set ylabel 'Comulative inflow (Mo)'");
  //  gp.commandln("plot 'extragas_inout.dat' using 1:7 with lin");
  // gp.commandln("pause -1");
  gp.end();
}

void plot3()
{
  Gnuplot_Tmpfile gp;
  gp.begin();
  gp.commandln("set terminal postscript portrait enhanced 16");
  gp.commandln("set out 'extragas3.ps'");
  gp.commandln("set nokey");
  gp.commandln("set size 1,.5");
  gp.commandln("set origin 0,0");
  gp.commandln("set xrange [0:90]");
  gp.commandln("set yrange [0:250]");
  gp.commandln("set xlabel '{/Symbol q} (degree)'");
  gp.commandln("set ylabel 'v_{kick} (km s^{-1})'");
  gp.commandln("plot 'extragas_veldistr.dat' using 1:3 with p pt 7, 25 w l");
  // panel 2
  gp.commandln("plot 'extragas_veldistr2.dat' using 1:3 with p pt 7");
  // panel 3
  gp.commandln("set multiplot");
  gp.commandln("set nokey");
  gp.commandln("set size 0.5,0.25");
  gp.commandln("set origin 0,0.25");
  gp.commandln("set xrange [0:90]");
  gp.commandln("set yrange [0:250]");
  gp.commandln("set xlabel 'theta (deg)'");
  gp.commandln("set ylabel 'v_{kick} (km/s)'");
  gp.commandln("plot 'extragas_veldistr.dat' using 1:3 with poi");
  gp.commandln("set origin .5,.25");
  gp.commandln("set xlabel 'theta (deg)'");
  gp.commandln("set ylabel 'v_R (km/s)'");
  gp.commandln("plot 'extragas_veldistr.dat' using 1:4 with poi");
  gp.commandln("set origin 0,0");
  gp.commandln("set xlabel 'theta (deg)'");
  gp.commandln("set ylabel 'v_{phi} (km/s)'");
  gp.commandln("set yrange [100:350]");
  gp.commandln("plot 'extragas_veldistr.dat' using 1:5 with poi");
  gp.commandln("set origin 0.5,0");
  gp.commandln("set xlabel 'theta (deg)'");
  gp.commandln("set ylabel 'v_z (km/s)'");
  gp.commandln("set yrange [0:250]");
  gp.commandln("plot 'extragas_veldistr.dat' using 1:6 with poi");
  // gp.commandln("pause -1");
  gp.end();
}

void plot4() 
{
	// TIME
  Gnuplot_Tmpfile gp;
  gp.begin();
  gp.commandln("set terminal postscript portrait enhanced 12");
  gp.commandln("set out 'extragas4.ps'");
  gp.commandln("set multiplot");
  gp.commandln("set nokey");
  gp.commandln("set log y");
  gp.commandln("set size 0.5,0.25");
  gp.commandln("set origin 0,.5");
  gp.commandln("set xrange [0:20]");
  gp.commandln("set yrange [1.e7:5.e9]");
  gp.commandln("set xlabel 'R (kpc)'");
  gp.commandln("set ylabel 'Travel time (yr)'");
  gp.commandln("plot 'extragas_time.dat' using 2:4 with poi pt 1, 'extragas_time.dat' u 5:6 with l");
  gp.commandln("set origin .5,.5");
  gp.commandln("set auto x");
  gp.commandln("set yrange [1.e7:5.e9]");
  gp.commandln("set xlabel 'V_{kick} (km/s)'");
  gp.commandln("plot 'extragas_time.dat' using 3:4 with poi pt 1");
  // gp.commandln("pause -1");
  gp.end();
}

void plot5()
{
  Gnuplot_Tmpfile gp;
  gp.begin();
  gp.commandln("set terminal postscript portrait enhanced 12");
  gp.commandln("set out 'extragas5.ps'");
  gp.commandln("set multiplot");
  gp.commandln("set nokey");
  gp.commandln("set log y");
  gp.commandln("set size 0.5,0.2");
  gp.commandln("set origin 0,.5");
  gp.commandln("set xrange [-13:13]");
  gp.commandln("set yrange [0.05:200]");
  gp.commandln("set xlabel 'z (kpc)'");
  gp.commandln("set ylabel 'Column density'");
  gp.commandln("set style line 1 lt 1 lw 1.6");
  gp.commandln("set style line 2 lt 2 lw 1");
  gp.commandln("plot 'z_profile1.dat' using 2:3 with poi 6, 'z_profile1.dat' using 2:4 w l ls 1, 0.2 w l ls 2");
  gp.commandln("set origin .5,.5");
  gp.commandln("plot 'z_profile2.dat' using 2:3 with poi 6, 'z_profile2.dat' using 2:4 w l ls 1, 0.2 w l ls 2");
  // gp.commandln("pause -1");
  gp.end();
}

void plot6()
{
  Gnuplot_Tmpfile gp;
  gp.begin();
  gp.commandln("set terminal postscript color enhanced 12");
  gp.commandln("set out 'extragas6.ps'");
  gp.commandln("set nokey");
  gp.commandln("set size 1,0.75");
  gp.commandln("set xrange [-180:180]");
  gp.commandln("set yrange [-90:90]");
  gp.commandln("set style line 1 lt 0 lw 1.6");
  gp.commandln("set style line 1 lt 2 lw 1");
  gp.commandln("unset border");
  gp.commandln("unset xtics");
  gp.commandln("unset x2tics");
  gp.commandln("unset ytics");
  gp.commandln("unset y2tics");
  gp.commandln("plot 'extragas_hvcs_n2.dat' u 3:4 with poi 3, 'extragas_hvcs_n1.dat' u 3:4 with poi 5, 'extragas_hvcs_0.dat' u 3:4 with poi 2, 'extragas_hvcs_p1.dat' u 3:4 with poi 6, 'extragas_hvcs_p2.dat' u 3:4 with poi 8, 'extragas_hvcs_axes.dat' u 3:4 with p 0");
  gp.commandln("");
  gp.commandln("");
  gp.commandln("");
  // gp.commandln("pause -1");
  gp.end();
}

void plot7(char make_fit)
{
  /* ROTATION VELOCITIES */
  Gnuplot_Tmpfile gp;
  gp.begin();
  gp.commandln("set terminal postscript enhanced 18");
  gp.commandln("set out 'extragas7.ps'");
  gp.commandln("set nokey");
  // panel 1
  gp.commandln("set size 1,1");
  gp.commandln("set xrange [0:18]");
  gp.commandln("set yrange [0:250]");
  gp.commandln("set xlabel 'R (kpc)'");
  gp.commandln("set ylabel 'Rotation velocity (km/s)'");
  gp.commandln("set style line 1 lt 1 lw 2");
  gp.commandln("set style line 2 lt 0 lw 2");
  gp.commandln("set style line 3 lt 4 lw 2");
  gp.commandln("set style line 4 lt 5 lw 2");
  gp.commandln("set multiplot");
  gp.commandln("set size 0.5,0.5");
  gp.commandln("set origin 0,0.5");
  gp.commandln("plot 'extragas_rotsurf.dat' u 1:2 w l ls 1, '../../Astro/plots/rotsurf_et_1.3.dat' u ($1):($2-11):($3) w yerr 5, 'extragas_rotsurf.dat' u 1:5 w l ls 3");
  gp.commandln("set origin .5,.5");
  gp.commandln("plot 'extragas_rotsurf.dat' u 1:2 w l ls 1, '../../Astro/plots/rotsurf_et_2.6.dat' u ($1):($2-11):($3) w yerr 5, 'extragas_rotsurf.dat' u 1:9 w l ls 3");
  gp.commandln("set origin 0,0");
  gp.commandln("plot 'extragas_rotsurf.dat' u 1:2 w l ls 1, '../../Astro/plots/rotsurf_et_3.9.dat' u ($1):($2-11):($3) w yerr 5, 'extragas_rotsurf.dat' u 1:13 w l ls 3");
  gp.commandln("set origin 0.5,0");
  gp.commandln("plot 'extragas_rotsurf.dat' u 1:2 w l ls 1, '../../Astro/plots/rotsurf_et_5.2.dat' u ($1):($2-11):($3) w yerr 5, 'extragas_rotsurf.dat' u 1:17 w l ls 3");
  // panel 2
  gp.commandln("set multiplot");
  gp.commandln("set origin 0,0.5");
  gp.commandln("plot 'extragas_rotsurf.dat' u 1:2 w l ls 1, '../../Astro/plots/rotsurf_et_1.3.dat' u ($1):($2-11):($3) w yerr 5, 'extragas_rotsurf.dat' u 1:6 w l ls 3");
  gp.commandln("set origin .5,.5");
  gp.commandln("plot 'extragas_rotsurf.dat' u 1:2 w l ls 1, '../../Astro/plots/rotsurf_et_2.6.dat' u ($1):($2-11):($3) w yerr 5, 'extragas_rotsurf.dat' u 1:10 w l ls 3");
  gp.commandln("set origin 0,0");
  gp.commandln("plot 'extragas_rotsurf.dat' u 1:2 w l ls 1, '../../Astro/plots/rotsurf_et_3.9.dat' u ($1):($2-11):($3) w yerr 5, 'extragas_rotsurf.dat' u 1:14 w l ls 3");
  gp.commandln("set origin 0.5,0");
  gp.commandln("plot 'extragas_rotsurf.dat' u 1:2 w l ls 1, '../../Astro/plots/rotsurf_et_5.2.dat' u ($1):($2-11):($3) w yerr 5, 'extragas_rotsurf.dat' u 1:18 w l ls 3");
  // panel 3 GRADIENT
  gp.commandln("unset multiplot");
  gp.commandln("set size 1,1");
  gp.commandln("set origin 0,0");
  gp.commandln("set xrange [0:6]");
  gp.commandln("set auto y");
  gp.commandln("set xlabel 'z (kpc)'");
  gp.commandln("set ylabel 'Rotation velocity (km/s)'");
  gp.commandln("vdisk=228");
  gp.commandln("dvel=-3.");
  gp.commandln("g(x)=vdisk+dvel*x");
  if (make_fit == 'y')
    gp.commandln("fit g(x) 'extragas_gradient2.dat' via dvel,vdisk");
  gp.commandln("plot 'extragas_gradient.dat' u 1:2 w p 2, 'extragas_gradient2.dat' u 1:2 w p 5, 'extragas_gradient2.dat' u 1:2:3 w yerr 4, g(x)");
  // panel 4
  gp.commandln("set xrange [0:18]");
  gp.commandln("set yrange [0:250]");
  gp.commandln("plot 'extragas_rotsurf.dat' u 1:2 w l ls 1, '../../Astro/plots/rotsurf_et_2.6.dat' u 1:2:3 w yerr 5, 'extragas_rotsurf.dat' u 1:9 w l ls 2, '../../Astro/plots/rotsurf_et_5.2.dat' u 1:2:3 w yerr 4, 'extragas_rotsurf.dat' u 1:17 w l ls 3, 'extragas_rotsurf.dat' u 1:25 w l ls 3");
  // panel 5
  gp.commandln("plot 'extragas_rotsurf.dat' u 1:2 w l ls 1, '../../Astro/plots/rotsurf_et_2.6.dat' u 1:2:3 w yerr 5, 'extragas_rotsurf.dat' u 1:10 w l ls 2, '../../Astro/plots/rotsurf_et_5.2.dat' u 1:2:3 w yerr 4, 'extragas_rotsurf.dat' u 1:18 w l ls 3, 'extragas_rotsurf.dat' u 1:26 w l ls 3");
  //  gp.commandln("plot 'extragas_rotsurf.dat' u 1:2 w l ls 1, 'extragas_rotsurf.dat' u 1:6 w l ls 2, 'extragas_rotsurf.dat' u 1:10 w l ls 3, 'extragas_rotsurf.dat' u 1:14 w l ls 2, 'extragas_rotsurf.dat' u 1:18 w l ls 3");
  //  gp.commandln("set multiplot");
  //  gp.commandln("plot 'extragas_rotsurf.dat' u 1:2 w l, 'extragas_rotsurf.dat' u 1:3 w l, 'extragas_rotsurf.dat' u 1:4 w p, 'extragas_rotsurf.dat' u 1:6 w p, 'extragas_rotsurf.dat' u 1:8 w p, 'extragas_rotsurf.dat' u 1:10 w p, 'extragas_rotsurf.dat' u 1:12 w p");
  // panel 6 RADIAL VELOCITIES
  gp.commandln("set multiplot");
  gp.commandln("set size 0.5,0.5");
  gp.commandln("set auto x");
  gp.commandln("set yrange [-50:50]");
  gp.commandln("set xlabel 'R (kpc)'");
  gp.commandln("set ylabel 'Radial velocity (km/s)'");
  gp.commandln("set origin 0,0.5");
  gp.commandln("plot 'extragas_radsurf.dat' u 1:2 w l ls 1, 'extragas_radsurf.dat' u 1:5 w l ls 3");
  gp.commandln("set origin .5,.5");
  gp.commandln("plot 'extragas_radsurf.dat' u 1:2 w l ls 1, 'extragas_radsurf.dat' u 1:9 w l ls 3");
  gp.commandln("set origin 0,0");
  gp.commandln("plot 'extragas_radsurf.dat' u 1:2 w l ls 1, 'extragas_radsurf.dat' u 1:13 w l ls 3");
  gp.commandln("set origin 0.5,0");
  gp.commandln("plot 'extragas_radsurf.dat' u 1:2 w l ls 1, 'extragas_radsurf.dat' u 1:17 w l ls 3");
  // panel 7
  gp.commandln("set multiplot");
  gp.commandln("set size 0.5,0.5");
  gp.commandln("set origin 0,0.5");
  gp.commandln("plot 'extragas_radsurf.dat' u 1:2 w l ls 1, 'extragas_radsurf.dat' u 1:6 w l ls 3");
  gp.commandln("set origin .5,.5");
  gp.commandln("plot 'extragas_radsurf.dat' u 1:2 w l ls 1, 'extragas_radsurf.dat' u 1:10 w l ls 3");
  gp.commandln("set origin 0,0");
  gp.commandln("plot 'extragas_radsurf.dat' u 1:2 w l ls 1, 'extragas_radsurf.dat' u 1:14 w l ls 3");
  gp.commandln("set origin 0.5,0");
  gp.commandln("plot 'extragas_radsurf.dat' u 1:2 w l ls 1, 'extragas_radsurf.dat' u 1:18 w l ls 3");
  // panel 8
  gp.commandln("unset multiplot");
  gp.commandln("set origin 0,0");
  gp.commandln("set size 1,1");
  gp.commandln("set xrange [0:18]");
  gp.commandln("set yrange [0:250]");
  gp.commandln("set auto y");
  gp.commandln("plot 'extragas_dens.dat' u 1:2 w l, 'extragas_dens.dat' u 1:4 w l, 'extragas_dens.dat' u 1:6 w l, 'extragas_dens.dat' u 1:8 w l");
  // panel 9
  gp.commandln("set multiplot");
  gp.commandln("set size 0.5,0.5");
  gp.commandln("set origin 0,.5");
  gp.commandln("plot 'extragas_dens.dat' u 1:2 w p");
  gp.commandln("set origin .5,.5");
  gp.commandln("plot 'extragas_dens.dat' u 1:4 w p");
  gp.commandln("set origin 0,0");
  gp.commandln("plot 'extragas_dens.dat' u 1:6 w p");
  gp.commandln("set origin .5,0");
  gp.commandln("plot 'extragas_dens.dat' u 1:8 w p");
  // gp.commandln("pause -1");
  gp.end();
}

void plot8()
{
  /* HOT HALO */
  Gnuplot_Tmpfile gp;
  gp.begin();
  gp.commandln("set terminal postscript portrait enhanced 16");
  gp.commandln("set out 'extragas8.ps'");
  gp.commandln("set nokey");
  gp.commandln("set size 1,.5");
  gp.commandln("set xrange [0.1:100]");
  gp.commandln("set yrange [.000005:5]");
  gp.commandln("set log");
  gp.commandln("set style line 1 lt 1 lw 2");
  gp.commandln("set style line 2 lt 2 lw 2");
  gp.commandln("set xlab 'Distance from the centre (kpc)' ");
  gp.commandln("set ylab 'Volume density (atoms cm^{-3})' ");
  gp.commandln("plot 'extragas_halo_R.dat' u 1:3 w l ls 1, 'extragas_halo_z.dat' u 1:3 w l ls 2");
  gp.commandln("unset log");
  gp.commandln("set xrange [0:30]");
  gp.commandln("plot 'extragas_halo_R.dat' u 1:3 w p 3, 'extragas_halo_z.dat' u 1:3 w p 4");
  gp.commandln("set xrange [0.1:30]");
  gp.commandln("set log");
  gp.commandln("plot 'extragas_halo_R.dat' u 1:3 w p 3, 'extragas_halo_z.dat' u 1:3 w p 4");
  // gp.commandln("pause -1");
  gp.end();
}

void plot9()
{
  /* RESIDUALS */
  Gnuplot_Tmpfile gp;
  gp.begin();
  gp.commandln("set terminal postscript color enhanced 12");
  gp.commandln("set out 'extragas9.ps'");
  gp.commandln("set key");
  gp.commandln("set size 1,1");
  //  gp.commandln("set xrange [0:30]");
  //  gp.commandln("set yrange [150:240]");
  gp.commandln("set style line 1 lt 0 lw 1.6");
  gp.commandln("set style line 1 lt 2 lw 1");
  gp.commandln("set xlab 'Par 1' ");
  gp.commandln("set ylab 'Residuals' ");
  gp.commandln("plot 'residuals.dat' u 1:3 w p, 'residuals.dat' u 1:4 w p, 'residuals.dat' u 1:5 w p, 'residuals.dat' u 1:6 w p");
  gp.commandln("plot 'residuals.dat' u 2:3 w p, 'residuals.dat' u 2:4 w p, 'residuals.dat' u 2:5 w p, 'residuals.dat' u 2:6 w p");
  gp.commandln("plot 'residuals.dat' u 1:7 w p, 'residuals.dat' u 1:8 w p, 'residuals.dat' u 1:9 w p, 'residuals.dat' u 1:10 w p");
  gp.commandln("plot 'residuals.dat' u 2:7 w p, 'residuals.dat' u 2:8 w p, 'residuals.dat' u 2:9 w p, 'residuals.dat' u 2:10 w p");
  gp.commandln("splot 'residuals.dat' using 1:2:5");
  gp.commandln("splot 'residuals.dat' using 1:2:6");
  // gp.commandln("pause -1");
  gp.end();
}

void PlotCurves(Array3D<float> out_v_phi,
		Array2D<float>& out_v_phi_wm,
		Array3D<float> out_v_R,
		Array2D<float> out_dens,
		float kpctograd,
		float rotsurf_grid,
		int rotsurf_nvel,
		float rotsurf_max,
		float clip_wm,
		float Rmin_gradient,
		float Rmax_gradient,
		float zmax_gradient,
		char make_fit,
		float incl
		)
{
  /*
    Plotting rotation curves at different heights
  */
  FILE* f_v_R;
  FILE* f_v_phi;
  FILE* f_dens;
  FILE* f_grad; 
  FILE* f_grad2;
  //f_v_phi=fopen("extragas_rotsurf.dat", "w");
  //f_v_R=fopen("extragas_radsurf.dat", "w");
  float norm_v_phi, norm_v_R;
  int rotsurf_nR=(int)(rotsurf_max/rotsurf_grid);
  int rotsurf_nz=(int)(rotsurf_max/rotsurf_grid);
  int i=0, imin, imax, jmax;
  double tR, tv, tz;

  // WEIGHTED MEAN in v_phi
  clip_wm*=max(out_v_phi);
  for (int i=0; i<rotsurf_nR; i++){
    for (int j=0; j<rotsurf_nz; j++){
      norm_v_phi=0.;
      for (int k=0; k<rotsurf_nvel; k++){
	if (out_v_phi[i][j][k] > clip_wm) {
	  // k is the velocity because we use a grid of 1 km/s
	  out_v_phi_wm[i][j]+=k*out_v_phi[i][j][k];
	  norm_v_phi+=out_v_phi[i][j][k];
	}
      }
      if (norm_v_phi != 0)
	out_v_phi_wm[i][j]=out_v_phi_wm[i][j]/norm_v_phi;
    }
  }
  Array2D<float> smoOut_v_phi_wm(rotsurf_nR,rotsurf_nz, 0.f);
  //  smoOut_v_phi_wm=out_v_phi_wm.smooth2D(1,1,3,1);
  smoOut_v_phi_wm=out_v_phi_wm.smooth1D_noZero('x',1,3);

  // WEIGHTED MEAN in v_R
  Array2D<float> out_v_R_wm(rotsurf_nR,rotsurf_nz, 0.f);
  for (int i=0; i<rotsurf_nR; i++){
    for (int j=0; j<rotsurf_nz; j++){
      norm_v_R=0.;
      for (int k=0; k<rotsurf_nvel; k++){
	if (out_v_R[i][j][k] > clip_wm) {
	  // we use a grid of 1 km/s
	  // centred in rotsurf_nvel/2
	  // velocities below rotsurf_nvel/2 are negative
	  out_v_R_wm[i][j]+=(k-rotsurf_nvel/2)*out_v_R[i][j][k];
	  norm_v_R+=out_v_R[i][j][k];
	}
      }
      if (norm_v_R != 0) {
	out_v_R_wm[i][j]=out_v_R_wm[i][j]/norm_v_R;
      }
    }
  }
  Array2D<float> smoOut_v_R_wm(rotsurf_nR,rotsurf_nz, 0.f);
  //  smoOut_v_phi_wm=out_v_phi_wm.smooth2D(1,1,3,1);
  smoOut_v_R_wm=out_v_R_wm.smooth1D_noZero('x',1,3);

  // WRITING ROTATION CURVES IN and ABOVE THE PLANE
  while (i < rotsurf_nR-1){
    // columns:
    // radius, vrot0, vrot[i], smooVrot[i]... (vrot0 is not smoothed)
    tR=i*rotsurf_grid;
    //fprintf(f_v_phi,"%.2f %.2f ", tR, out_v_phi_wm[i][0]);
    for (int j=1; j<rotsurf_nz; j++){
      // cout << i << " " << j << " " << tR  << "\n";
      if (out_v_phi_wm[i][j] != 0) {
	//fprintf(f_v_phi,"%.2f %.2f ", out_v_phi_wm[i][j], 
		//smoOut_v_phi_wm[i][j]);
      }
      else {
	//fprintf(f_v_phi,"- - ");
      }
    }
    //fprintf(f_v_phi,"\n");
    i++;
  }
  //fclose(f_v_phi);

  // WRITING RADIAL VELOCITIES IN and ABOVE THE PLANE
  i=0;
  while (i < rotsurf_nR-1){
    // columns:
    // radius, vR0, vR[i], smooVR[i]... (vR0 is not smoothed)
    tR=i*rotsurf_grid;
    //fprintf(f_v_R,"%.2f %.2f ", tR, out_v_R_wm[i][0]);
    // WRITING ROTATION CURVES ABOVE THE PLANE
    for (int j=1; j<rotsurf_nz; j++){
      //cout << i << " " << j << " " << tR  << " " << out_v_R_wm[i][0] << "\n";
      if (out_v_R_wm[i][j] != 0) {
	//fprintf(f_v_R,"%.2f %.2f ", out_v_R_wm[i][j], 
		//smoOut_v_R_wm[i][j]);
      }
      else {
	//fprintf(f_v_R,"- - ");
      }
    }
    //fprintf(f_v_R,"\n");
    i++;
  }
  //fclose(f_v_R);

  //    Gradient
  //f_grad=fopen("extragas_gradient.dat", "w");
  //f_grad2=fopen("extragas_gradient2.dat", "w");
  imin=(int)(Rmin_gradient/rotsurf_grid)+1;
  imax=(int)(Rmax_gradient/rotsurf_grid)+1;
  jmax=(int)(zmax_gradient/rotsurf_grid)+1;
  Array1D<double> aveVel(jmax-1,0.f), sigVel(jmax-1,0.f);
  int nVel;
  for (int j=1; j<jmax; j++){
    nVel=0;
    for (int i=imin; i<imax; i++){
      tz=j*rotsurf_grid;
      if (out_v_phi_wm[i][j] != 0) {
	//fprintf(f_grad,"%.2f %.2f \n", tz, out_v_phi_wm[i][j]);
	aveVel[j-1]+=out_v_phi_wm[i][j];
	sigVel[j-1]+=(out_v_phi_wm[i][j]*out_v_phi_wm[i][j]);
	nVel+=1;
      }
    }
    aveVel[j-1]/=nVel;
    sigVel[j-1]=sqrt(1./(nVel-1)*(sigVel[j-1]-nVel*aveVel[j-1]*aveVel[j-1]));
  }
  for (int j=1; j<jmax; j++){
    tz=j*rotsurf_grid;
    //fprintf(f_grad2,"%.2f %.2f %.2f \n", tz, aveVel[j-1], sigVel[j-1]);
  }
  //fclose(f_grad);
  //fclose(f_grad2);

  // Plotting velocity field  
  for (int i=0; i<rotsurf_nR; i++){
    for (int j=0; j<rotsurf_nz; j++){
      out_v_phi_wm[i][j]=out_v_phi_wm[i][(int)(j*cos(incl/RAD))];
    }
  }

  /*
    Plotting density
  */
  //f_dens=fopen("extragas_dens.dat", "w");
  for (int i=0; i<rotsurf_nR; i++){
    tR=i*rotsurf_grid;
    //fprintf(f_dens,"%.2f ", tR);
    for (int j=0; j<rotsurf_nz; j++){
      //      tz=j*rotsurf_grid;
      //fprintf(f_dens,"%.5f ", out_dens[i][j]);
    }
    //fprintf(f_dens,"\n");
  }
  //fclose(f_dens);
//  plot7(make_fit);
}
