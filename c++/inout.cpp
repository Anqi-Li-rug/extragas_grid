/*
  INOUT

  inout is a class that deals with input and output rates in extragas.
  
*/
#include "inout.h"
#include "findValue.h"
#include "numint.h"
#include "galFunctions.h"
#include "extragasPlots.h"
#include <stdio.h>
const double YR=3.155691e7; // yr -> sec
const double CONV2=1.98892e43; // Mo*km^2/s^2 -> erg
extern char inputfile[];
extern bool verbose;
extern std::vector<float> DSFR;

using namespace std;

inout::inout(int nR1, float dR)
{
  nR=nR1;
  deltaR=dR;
  Array1D<float> inout1(nR, 0.f);
  dens_out=inout1;
  Array1D<float> inout2(nR, 0.f);
  dens_in=inout2; 
  Array1D<float> inout3(nR, 0.f);
  dens_out_above=inout3; 
  Array1D<float> inout4(nR, 0.f);
  dens_in_above=inout4; 
  Array1D<float> inout5(nR, 0.f);
  tmass_out=inout5; 
  Array1D<float> inout6(nR, 0.f);
  tmass_in=inout6; 
  Array1D<float> inout7(nR, 0.f);
  tmass_out_above=inout7; 
  Array1D<float> inout8(nR, 0.f);
  tmass_in_above=inout8; 
  Array1D<float> inout9(nR, 0.f);
  dens_accr=inout9; 
  Array1D<float> inout10(nR, 0.f);
  tmass_accr=inout10; 
  Global_acc_rate = 0.;
}

void inout::back_to_zero()
{
	for (int i=0;i<nR;i++)
	{
		dens_out[i] = 0.;
		dens_in[i] = 0.;
		dens_out_above[i] = 0.;
		dens_in_above[i] = 0.;
		tmass_out[i] = 0.;
		tmass_in[i] = 0.;
		tmass_out_above[i] = 0.;
		tmass_in_above[i] = 0.;
		dens_accr[i] = 0.;
		tmass_accr[i] = 0.;
	}
	Global_acc_rate = 0.;
}
void inout::CalculateRates(double R_in,
			   float t,
			   float t_above,
			   float dt,
			   float tau,
			   int n_tau,
			   float R1,
			   float R2,
			   float cn_dot,
			   float cn_accr,
			   float v_k,
			   float R_in_above,
			   int above,
			   double& ctotal_energy_above)
{  // cout<<"R"<<R_in<<"t_above"<<t_above/t<<"R_above"<<R_in_above<<endl;
    findValue(inputfile, single_event, "single_event");
  float RmaxSF, RminSF; 
        //float qromo(HoleExpSchmidt,RminSF,RmaxSF,'p',1.e-6)=1.0;	
	if(single_event=='y') {RminSF=R1; RmaxSF=R2;}
	else
	{
		findValue(inputfile, RmaxSF, "RmaxSF");
		findValue(inputfile, RminSF, "RminSF");
	}
 
  // OUTPUT DENSITY
  i_out=(int)((R1+R2)/2./deltaR); /* +0.5 ??? */
  //dens_out[i_out]+=cn_dot*(t/dt)/(t*1.e6)*2.; /* mass s^-1 */
  dens_out[i_out]+=cn_dot*(t/dt)/(t*1.e6)*2.;
  dens_out_above[i_out]+=cn_dot/dt/1.e6*2.*t_above/t;

  //cout << i_out << " " << R1 << " " << R2  << " " << t << " " << dt << " "  << dens_out[i_out] << " " <<cn_accr/cn_dot<< endl;

  // OUTPUT MASS
  if (single_event == 'n') {
    tmass_out[i_out]+=qromo(HoleExpSchmidt_R,R1,R2,'p',1.e-6)/qromo(HoleExpSchmidt,RminSF,RmaxSF,'p',1.e-6)*(t/dt)/(t*1.e6)*2.;
    tmass_out_above[i_out]+=qromo(HoleExpSchmidt_R,R1,R2,'p',1.e-6)/qromo(HoleExpSchmidt,RminSF,RmaxSF,'p',1.e-6)/dt/1.e6*2.*t_above/t;
  }
  else {
    tmass_out[i_out]+=qromo(HoleExpSchmidt,R1,R2,'p',1.e-6)*deltaR/qromo(HoleExpSchmidt,RminSF,RmaxSF,'p',1.e-6)*(t/dt)/(t*1.e6)*2.;
    tmass_out_above[i_out]+=qromo(HoleExpSchmidt,R1,R2,'p',1.e-6)*deltaR/qromo(HoleExpSchmidt,RminSF,RmaxSF,'p',1.e-6)/dt/1.e6*2.*t_above/t;
  }
  i_in=0;

  // ************** CHECK THE VALUE OF THE ENERGY!!!! ***********

  // OUTPUT ENERGY
  if (above == 1) {
    if (single_event == 'n')
      ctotal_energy_above+=.5*qromo(HoleExpSchmidt_R,R1,R2,'p',1.e-6)/qromo(HoleExpSchmidt,RminSF,RmaxSF,'p',1.e-6)/dt/1.e6*2.*t_above/t*v_k*v_k;
    else
      ctotal_energy_above+=.5*qromo(HoleExpSchmidt,R1,R2,'p',1.e-6)/qromo(HoleExpSchmidt,RminSF,RmaxSF,'p',1.e-6)*deltaR/dt/1.e6*2.*t_above/t*v_k*v_k;
  }

  Rmin_in=(RminSF/deltaR-(int)(RminSF/deltaR))*deltaR;
  found=false;
  found_above=false;

  // INPUT DENSITY AND MASS
  while (i_in < nR) { 
    // *************** INCLUDE THIS *****************
    //&& found == false && found_above == false) {
    R1_in=Rmin_in+deltaR*i_in-deltaR/2.;
    R2_in=Rmin_in+deltaR*i_in+deltaR/2.;
    if (R1_in < 0.) {R1_in=0.;}
    //    cout << R_in << " " << R2_in << " " << R1_in << endl;   
    //if (R_in<=Rmin_in-deltaR/2.){dens_accr[0]+=cn_accr/dt/1.e6*2;dens_in[0]+=HoleExpSchmidt((R2+R1)/2.)/qromo(HoleExpSchmidt,RminSF,RmaxSF,'p',1.e-6)/dt/1.e6*2.+cn_accr/dt/1.e6*2;}
    if (R_in <= R2_in and R_in >= R1_in) {
      dens_accr[i_in]+=cn_accr/dt/1.e6*2;
      // remember:
      //     cn_dot=HoleExpSchmidt(R_ave)/
      //	    qromo(HoleExpSchmidt,RminSF,RmaxSF,'p',1.e-6);
      dens_in[i_in]+=HoleExpSchmidt((R2+R1)/2.)/qromo(HoleExpSchmidt,RminSF,RmaxSF,'p',1.e-6)/dt/1.e6*2.*(R2+R1)/(2.*R_in)+cn_accr/dt/1.e6*2;
        //dens_in[i_in]+=cn_dot/dt/1.e6*2.+cn_accr/dt/1.e6*2;
      if (single_event == 'n') {
	tmass_accr[i_in]+=2.*PI*R_in*cn_accr*deltaR/dt/1.e6*2;
	tmass_in[i_in]+=qromo(HoleExpSchmidt_R,R1,R2,'p',1.e-6)
	  /qromo(HoleExpSchmidt,RminSF,RmaxSF,'p',1.e-6)/dt/1.e6*2.
	  +2.*PI*R_in*cn_accr*deltaR/dt/1.e6*2;
      }
      else {
	tmass_in[i_in]+=qromo(HoleExpSchmidt,R1,R2,'p',1.e-6)*deltaR
	  /qromo(HoleExpSchmidt,RminSF,RmaxSF,'p',1.e-6)/dt/1.e6*2.;
      }
      found=true;
    }
    // ABOVE
    if (R_in_above <= R2_in and R_in_above >= R1_in) {
      dens_in_above[i_in]+=
	HoleExpSchmidt((R2+R1)/2.)/qromo(HoleExpSchmidt,RminSF,RmaxSF,'p',1.e-6)
	/dt/1.e6*2.*t_above/t*(R2+R1)/(2.*R_in);
      if (single_event == 'n')
	tmass_in_above[i_in]+=qromo(HoleExpSchmidt_R,R1,R2,'p',1.e-6)
	  /qromo(HoleExpSchmidt,RminSF,RmaxSF,'p',1.e-6)/dt/1.e6*2.*t_above/t;
      else 
	tmass_in_above[i_in]+=qromo(HoleExpSchmidt,R1,R2,'p',1.e-6)*deltaR
	  /qromo(HoleExpSchmidt,RminSF,RmaxSF,'p',1.e-6)/dt/1.e6*2.*t_above/t;
      found_above=true;
    }
    i_in++;
  }
}

void inout::NormAndPlot(float ntomass,
			float mean_v_k,
			float M_tot_pbc
			)
{
  float Rmax_disc; findValue(inputfile, Rmax_disc, "Rmax_disk"); 
  float Rmin_disc; findValue(inputfile, Rmin_disc, "Rmin_disk"); 
  FILE* f_inout=NULL;
  f_inout=fopen("extragas_inout.dat", "w");

  Array1D<float> sR(nR,0.f),sdens_in(nR,0.f),sdens_out(nR,0.f),
    sdens_accr(nR,0.f);   // smoothed
  Array1D<float> cum_mass_out(nR,0.f),cum_mass_in(nR,0.f),
    cum_mass_accr(nR,0.f),
    cum_mass_out_above(nR,0.f),cum_mass_in_above(nR,0.f);
  Array1D<float> dens_now(nR,0.f), dens0(nR,0.f),
    dens0_2(nR,0.f), dens0_3(nR,0.f);
  Array1D<float> total_mass(nR,0.f);
  Array1D<double> energy_out(nR,0.f), cum_energy_out(nR,0.f);

  for (int iii=0; iii<nR/2; iii++){
    ii=1+iii*2;
    sR[iii]=Rmin_in+(deltaR*ii*2+deltaR)/2.;
    sdens_in[iii]=(dens_in[ii-1]+dens_in[ii])/2.;
    sdens_out[iii]=(dens_out[ii-1]+dens_out[ii])/2.;
    sdens_accr[iii]=(dens_accr[ii-1]+dens_accr[ii])/2.;
  }
  dens_out=dens_out*ntomass;
  dens_in=dens_in*ntomass;
  dens_accr=dens_accr*ntomass;
  //sdens_in=sdens_in*ntomass;
  //sdens_out=sdens_out*ntomass;
sdens_in=tmass_in*ntomass;
sdens_out=tmass_out*ntomass;  
sdens_accr=sdens_accr*ntomass;
  dens_in_above=dens_in_above*ntomass;
  dens_out_above=dens_out_above*ntomass;
  total_mass=(tmass_out-tmass_in)*(ntomass*1.e10);

  for (ii=0; ii<nR; ii++){
    for (int ii2=0; ii2<ii; ii2++){
      cum_mass_out[ii]+=tmass_out[ii2]*ntomass;
      cum_mass_in[ii]+=tmass_in[ii2]*ntomass;
      cum_mass_accr[ii]+=tmass_accr[ii2]*ntomass;
      cum_mass_out_above[ii]+=tmass_out_above[ii2]*ntomass;
      cum_mass_in_above[ii]+=tmass_in_above[ii2]*ntomass;
    }
    cum_energy_out[ii]=.5*cum_mass_out[ii]*mean_v_k*mean_v_k*CONV2/YR;
    energy_out[ii]=dens_out[ii]*(.5*mean_v_k*mean_v_k*CONV2);
  
    //Change in radial density distribution 
    dens_now[ii]=
      HoleExp(Rmin_in+deltaR*ii)
      /qromo(HoleExp_R,Rmin_disc,Rmax_disc,'p',1.e-6)*M_tot_pbc;
    dens0[ii]=
      HoleExp(Rmin_in+deltaR*ii)/qromo(HoleExp_R,Rmin_disc,Rmax_disc,'p',1.e-6)
      *M_tot_pbc+(dens_out[ii]-dens_in[ii])*tau/n_tau*1.e6;
    dens0_2[ii]=
      HoleExp(Rmin_in+deltaR*ii)/qromo(HoleExp_R,Rmin_disc,Rmax_disc,'p',1.e-6)
      *M_tot_pbc+(dens_out[ii]-dens_in[ii])*10.e9;
    dens0_3[ii]=
      HoleExp(Rmin_in+deltaR*ii)/qromo(HoleExp_R,Rmin_disc,Rmax_disc,'p',1.e-6)
      *M_tot_pbc+(dens_out[ii]-dens_in[ii])*10.e10;
    fprintf (f_inout,"%.2f %.5f %.5f %.5f %.2e %.2f %.2f %.2f %.2e %.2e %.2e %.2e %.2e %.2f %.2f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f \n", 
	     Rmin_in+deltaR*ii+deltaR/2., //#1
	     dens_out[ii], //#2
	     dens_in[ii], //#3
	     dens_out[ii]-dens_in[ii], //#4
	     total_mass[ii], //#5
	     cum_mass_out[ii], //#6
	     cum_mass_in[ii], //#7
	     cum_mass_out[ii]-cum_mass_in[ii], //#8
	     energy_out[ii], //#9
	     cum_energy_out[ii], //#10
	     dens_now[ii], //#11
	     dens0[ii], //#12
	     dens0_2[ii], //#13
	     dens0_3[ii], //#14
	     sR[ii], //#15
	     sdens_out[ii], //#16
	     sdens_in[ii], //#17
	     dens_out_above[ii], //#18
	     dens_in_above[ii], //#19
	     cum_mass_out_above[ii], //#20
	     cum_mass_in_above[ii], //#21
	     cum_mass_out_above[ii]-cum_mass_in_above[ii],//#22
	     dens_accr[ii],//#23
	     sdens_accr[ii],//#24
	     cum_mass_accr[ii]);//#25
  }
  if(verbose) printf ("Global accretion rate = %.3f Mo/yr \n", cum_mass_accr[nR-1]);
  Global_acc_rate = cum_mass_accr[nR-1];	
  fclose(f_inout);
}

inout::~inout()
{
}
