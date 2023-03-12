/*
  FUNCTIONS

  Mathematical and geometriacal functions to be use by exgas
*/

#include <cmath>
#include <iostream>
#include <vector>
#include "array.h"
using namespace std;
using namespace ffarray;

extern const double PI;//=3.14159265359;
const double RAD=180./PI;
/*
float HoleExpSchmidt(float R) 
{
  extern Array1D<float> FunctionPar;
  float Re=FunctionPar[0], alpha=FunctionPar[1], gammaSF=FunctionPar[2]; 
  float nDot_val=0.3*pow(R/Re+1,alpha*gammaSF)*exp(-gammaSF*(R/Re));
  return nDot_val;
}

float HoleExpSchmidt_R(float R)
{
  extern Array1D<float> FunctionPar;
  float Re=FunctionPar[0], alpha=FunctionPar[1], gammaSF=FunctionPar[2]; 
  float nDotR_val=0.3*2.*PI*R*pow(R/Re+1,alpha*gammaSF)*exp(-gammaSF*(R/Re));
  return nDotR_val;
}
*/
/*float HoleExpSchmidt(float R) 
{
	extern vector<float> RHI;
	extern vector<float> DHI;
	extern vector<float> RH2;
	extern vector<float> DH2;
	extern float Kennicut_T;
	extern float gammaSF;
	
	float valueHI, valueH2, step;
	
	
	int sizeH2 = RH2.size();
	
	if(R<=RH2[0]) {valueH2 = DH2[0];}
	else if (R>RH2[sizeH2-1]) {valueH2 = 0.;}
	else
	{
		for(int i=0; i<sizeH2-1; i++)
		{
			step = fabs(RH2[i+1]-RH2[i]);
			if(fabs(R-RH2[i])<=step && fabs(R-RH2[i+1])<=step)
			{
				valueH2 = DH2[i] + (DH2[i+1]-DH2[i])*(R-RH2[i])/(RH2[i+1]-RH2[i]);
				break;
			}
		}
	}
	return(pow(valueH2, gammaSF));
}
*/

 float HoleExpSchmidt(float R)
 {
     extern std::vector<float> RSFR;
     extern std::vector<float> DSFR;
     float r_sfr;
     int len_r=RSFR.size();
     if (R<RSFR[0]) {r_sfr=DSFR[0]; return r_sfr;}
     if (R>RSFR[len_r-1]) {r_sfr=0; return r_sfr;}
     for(int i=0; i<len_r-1; i++){
     if ((R>=RSFR[i])&(R<=RSFR[i+1])){
        r_sfr=(DSFR[i+1]-DSFR[i])*(R-RSFR[i])/(RSFR[i+1]-RSFR[i])+DSFR[i];
        return r_sfr;}}
 }


float HoleExpSchmidt_R(float R)
 {
     extern std::vector<float> RSFR;
     extern std::vector<float> DSFR;
     float r_sfr;
     int len_r=RSFR.size();
     if (R<RSFR[0]) {r_sfr=DSFR[0]; return 2.*PI*R*r_sfr;}
     if (R>RSFR[len_r-1]) {r_sfr=0; return 2.*PI*R*r_sfr;}
     for(int i=0; i<len_r-1; i++){
     if ((R>=RSFR[i])&(R<=RSFR[i+1])){
        r_sfr=(DSFR[i+1]-DSFR[i])*(R-RSFR[i])/(RSFR[i+1]-RSFR[i])+DSFR[i];
        return 2.*PI*R*r_sfr;}}
 }
 
/*float HoleExpSchmidt_R(float R) 
{
	extern vector<float> RHI;
	extern vector<float> DHI;
	extern vector<float> RH2;
	extern vector<float> DH2;
	extern float Kennicut_T;
	extern float gammaSF;
	
	float valueHI, valueH2, step;
	/*
	int sizeHI = RHI.size();
	
	if(R<=RHI[0]) valueHI = DHI[0];
	else if (R>RHI[sizeHI-1]) valueHI = 0.;
	else
	{
		for(int i=0; i<sizeHI-1; i++)
		{
			step = fabs(RHI[i+1]-RHI[i]);
			if(fabs(R-RHI[i])<=step && fabs(R-RHI[i+1])<=step)
			{
				valueHI = DHI[i] + (DHI[i+1]-DHI[i])*(R-RHI[i])/(RHI[i+1]-RHI[i]);
				break;
			}
		}
	}
	*/
	/*int sizeH2 = RH2.size();
	
	if(R<=RH2[0]) valueH2 = DH2[0];
	else if (R>RH2[sizeH2-1]) valueH2 = 0.;
	else
	{
		for(int i=0; i<sizeH2-1; i++)
		{
			step = fabs(RH2[i+1]-RH2[i]);
			if(fabs(R-RH2[i])<=step && fabs(R-RH2[i+1])<=step)
			{
				valueH2 = DH2[i] + (DH2[i+1]-DH2[i])*(R-RH2[i])/(RH2[i+1]-RH2[i]);
				break;
			}
		}
	}
	return 2.*PI*R*pow(valueH2, gammaSF); 
}*/



/*
// Funzioni a gradino per test
 float HoleExpSchmidt(float R) 
 {
	if (R<5.) {return 1.;}
	else {return 0.;}
 }
 
 float HoleExpSchmidt_R(float R)
 {
	if (R<5.) {return 2.*PI*R;}
	else {return 0.;}
 }
*/

float ExpDisc(float R)
{
  extern Array1D<float> FunctionPar;
  /* Exponential disk */
  float Exp_val=exp(-R/FunctionPar[0]);
  return Exp_val;
}

float ExpDisc_R(float R)
{
  extern Array1D<float> FunctionPar;
  /* Exponential disk */
  float ExpR_val=2.*PI*R*exp(-R/FunctionPar[0]);
  return ExpR_val;
}

float ExpDisc0(float R)
{
  extern Array1D<float> FunctionPar;
  /* Exponential disk */
  float Exp_val=FunctionPar[0]*exp(-R/FunctionPar[1]);
  return Exp_val;
}

float ExpDisc0_R(float R)
{
  extern Array1D<float> FunctionPar;
  /* Exponential disk */
  float Exp_val=2*PI*R*FunctionPar[0]*exp(-R/FunctionPar[1]);
  return Exp_val;
}

float HoleExp(float R)
{
  extern Array1D<float> FunctionPar;
  /* Exponential disk with inner depression */
  float HoleExp_val=pow(R/FunctionPar[0]+1,FunctionPar[1])
    *exp(-R/FunctionPar[0]);
  return HoleExp_val;
}

float HoleExp_R(float R)
{
  extern Array1D<float> FunctionPar;
  /* Exponential disk with inner depression */
  float HoleExp_val=2*PI*R*pow(R/FunctionPar[0]+1,FunctionPar[1])
    *exp(-R/FunctionPar[0]);
  return HoleExp_val;
}

float HIsurf(float R)
{
	extern vector<float> RHI;
	extern vector<float> DHI;
	
	float valueHI, step;
	
	int sizeHI = RHI.size();
	
	if(R<=RHI[0]) {valueHI = DHI[0];}
	else if (R>RHI[sizeHI-1]) {valueHI = 0.;}
	else
	{
		for(int i=0; i<sizeHI-1; i++)
		{
			step = fabs(RHI[i+1]-RHI[i]);
			if(fabs(R-RHI[i])<=step && fabs(R-RHI[i+1])<=step)
			{
				valueHI = DHI[i] + (DHI[i+1]-DHI[i])*(R-RHI[i])/(RHI[i+1]-RHI[i]);
				break;
			}
		}
	}
	return valueHI;
}

float HIsurf_R(float R)
{
	extern vector<float> RHI;
	extern vector<float> DHI;
	
	float valueHI, step;
	
	int sizeHI = RHI.size();
	
	if(R<=RHI[0]) {valueHI = DHI[0];}
	else if (R>RHI[sizeHI-1]) {valueHI = 0.;}
	else
	{
		for(int i=0; i<sizeHI-1; i++)
		{
			step = fabs(RHI[i+1]-RHI[i]);
			if(fabs(R-RHI[i])<=step && fabs(R-RHI[i+1])<=step)
			{
				valueHI = DHI[i] + (DHI[i+1]-DHI[i])*(R-RHI[i])/(RHI[i+1]-RHI[i]);
				break;
			}
		}
	}
	return 2.*PI*R*valueHI;
}

float HIsurf_R_hs(float R)
{
	extern vector<float> RHI;
	extern vector<float> DHI;
	extern Array1D<float> FunctionPar;
	
	
	float valueHI, step;
	
	int sizeHI = RHI.size();
	
	if(R<=RHI[0]) {valueHI = DHI[0];}
	else if (R>RHI[sizeHI-1]) {valueHI = 0.;}
	else
	{
		for(int i=0; i<sizeHI-1; i++)
		{
			step = fabs(RHI[i+1]-RHI[i]);
			if(fabs(R-RHI[i])<=step && fabs(R-RHI[i+1])<=step)
			{
				valueHI = DHI[i] + (DHI[i+1]-DHI[i])*(R-RHI[i])/(RHI[i+1]-RHI[i]);
				break;
			}
		}
	}
	
	float hs;
	if(R>5.) {hs = FunctionPar[0]*exp((R-8.5)/FunctionPar[1]);}
	else{hs = FunctionPar[0]*exp((5.-8.5)/FunctionPar[1]);}
	
	return PI*R*valueHI/hs;
}

float PowerLaw(float r)
{
  extern Array1D<float> FunctionPar;
  /* Power law */
  float PL_val=pow(r/FunctionPar[0],FunctionPar[1]);
  return PL_val;
}

float PowerLaw0(float r)
{
  extern Array1D<float> FunctionPar;
  /* Power law */
  float PL_val=FunctionPar[0]*pow(r/FunctionPar[1],FunctionPar[2]);
  return PL_val;
}

float PowerLaw_r(float r)
{
  extern Array1D<float> FunctionPar;
  /* Power law */
  float PLr_val=4*PI*r*r*pow(r/FunctionPar[0],FunctionPar[1]);
  return PLr_val;
}

float PowerLaw0_r(float r)
{
  extern Array1D<float> FunctionPar;
  /* Power law */
  float PLr_val=4*PI*r*r*FunctionPar[0]*pow(r/FunctionPar[1],FunctionPar[2]);
  return PLr_val;
}

float DoublePowerLaw(float r)
{
  extern Array1D<float> FunctionPar;
  /* Power law */
  float DPL_val=pow(r/FunctionPar[0],FunctionPar[1]);
  return DPL_val;
}

float DoublePowerLaw0(float r)
{
  extern Array1D<float> FunctionPar;
  /* Power law */
  float DPL_val=FunctionPar[0]*pow(r/FunctionPar[1],FunctionPar[2]);
  return DPL_val;
}

float DoublePowerLaw_r(float r)
{
  extern Array1D<float> FunctionPar;
  /* Power law */
  float DPLr_val=4*PI*r*r*pow(r/FunctionPar[1],FunctionPar[2]);
  return DPLr_val;
}

float DoublePowerLaw0_r(float r)
{
  extern Array1D<float> FunctionPar;
  /* Power law */
  float DPLr_val=4*PI*r*r*FunctionPar[0]*pow(r/FunctionPar[1],FunctionPar[2]);
  return DPLr_val;
}

void cyltocar (double& x1, double& x2, double& x3, double& x4, double& x5, double& x6)
// READS PHASE SPACE CILINDRICAL COORDINATES AND RETURNS CARTESIANS
{
  double px2=x2, px4=x4;
  x2=x1*sin(x2/RAD);
  x1=x1*cos(px2/RAD);
  //  x3=;
  x4=x4*cos(px2/RAD)-x5*sin(px2/RAD);
  x5=px4*sin(px2/RAD)+x5*cos(px2/RAD);
  // x6=; 
}

void cartocyl (double& x1, double& x2, double& x3, double& x4, double& x5, double& x6)
// READS PHASE SPACE CARTESIAN COORDINATES AND RETURNS CILINDRICALS
{
  double px1=x1, px2=x2, px4=x4;
  x1=sqrt(x1*x1+x2*x2);
  x2=atan2(px2,px1)*RAD;
  if (x2 <= 0) {x2=x2+360;}
  //  z=;
  x4=(px1*x4+px2*x5)/x1;
  x5=(px1*x5-px2*px4)/x1;
  //  v_z=;
}

void sphertocar (double& x1, double& x2, double& x3, double& x4, double& x5, double& x6)
// READS PHASE SPACE SPHERICAL COORDINATES AND RETURNS CARTESIANS
{
  double px1=x1, px4=x4;
  x4=x4*cos(x2)*sin(x3)-x5*sin(x2)+x6*cos(x2)*cos(x3);
  x5=px4*sin(x2)*sin(x3)+x5*cos(x2)+x6*sin(x2)*cos(x3);
  x6=px4*cos(x3)-x6*sin(x3); //

  x1=x1*cos(x2)*sin(x3);
  x2=px1*sin(x2)*sin(x3);//
  x3=px1*cos(x3);//
}

void cartospher (double& x1, double& x2, double& x3, double& x4, double& x5, double& x6)
// READS PHASE SPACE CARTESIAN COORDINATES AND RETURNS SPHERICAL
{
  double px1=x1, px3=x3, px4=x4;
  x1=sqrt(x1*x1+x2*x2+x3*x3); //
  x3=acos(x3/x1);//

  x4=(px1*x4+x2*x5+px3*x6)/x1;//
  x5=x1*sin(x3)*(px1*x5-x2*px4)/(px1*px1+x2*x2);
  x6=(px3*x4-x1*x6)/sqrt(px1*px1+x2*x2); //

  x2=atan2(x2,px1);
  if (x2 <= 0) {x2=x2+2*PI;}
}

