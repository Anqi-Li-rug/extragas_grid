/* los.cpp
   PROJECTIONS ALONG THE LINE OF SIGHT 
   
   - INPUT: SPACE CILINDRIC COORDINATES (R,phi,z,v_R,v_phi,v_z) of a particle
   - OUTPUT: ON-SKY PROJECTED COORDINATES (xi,yi,v_los)
   
   2 possible projections:
   1) EXTERNAL: seen from outside at an inclination (incl) with respect to the plane
   2) INTERNAL: as seen from inside in the plane at a distance (R0) from the centre
 */

#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <cstdio>
#include "projections.h"
using namespace std;

// DECLARATIONS

extern const double PI;//=3.14159265359;
const double RAD=180./PI;

// CONSTRUCTORS

projections::projections(double cR,
			 double cphi, 
			 double cz, 
			 double cv_R, 
			 double cv_phi, 
			 double cv_z)
{
  R=cR, phi=cphi, z=cz, v_R=cv_R, v_phi=cv_phi, v_z=cv_z;
}

projections::projections() {}

projections::~projections() {}

// MEMBER FUNCTIONS


/*
  - EXTERNAL VIEW: 

  The observer is outside the system, the line of nodes is along Y, 
  the inclination (incl) of the symmetry axis and the line of sight is a 
  free parameter.
*/
void projections::external_view(float incl, float PA, float& xi, float& yi, float& vrad)
{
  xi=x_proj(R,phi);
  yi=y_proj(z,incl,R,phi);
  vrad=v_proj(v_R,v_phi,v_z,incl,phi);
   // float cosPA;
    //float sinPA;
    //float temp;
    //cosPA=cos(PA/RAD-PI*0.5);
    //sinPA=sin(PA/RAD-PI*0.5);
    //temp = R*cos(phi/RAD)*cos(PI-incl/RAD)-z*sin(PI-incl/RAD);
    //xi=R*sin(phi/RAD)*cosPA-temp*sinPA;
    //yi=R*sin(phi/RAD)*sinPA+temp*cosPA;
    //vrad=v_proj(v_R,v_phi,v_z,incl,phi);
}

/*
   - INTERNAL VIEW:

   The observer is inside the system, on the plane of symmetry, at a 
   distance R0 from the centre and moving at a velocity v0_phi_R0. The 
   azimuthal (l) coordinate is set at 0 towards the centre of the system.
*/
void projections::internal_view(float R0,float v0_phi_R0, float& l, float& b, float& d, float& vlos)
{
  l=longitude(R,phi,R0);
  l=l*RAD;
  b=latitude(R,phi,z,l,R0);
  b=b*RAD;
  d=distance(R,z,l,R0,phi); // in kpc
  vlos=v_los(v_R, v_phi, v_z, l, b, phi);
  vlos=vlos-v0_phi_R0*sin(l/RAD)*cos(b/RAD);
}

void projections::internal_view_hvcs(float R0,float v0_phi_R0, float Rmax_hvcs, float zmax_hvcs, float& l, float& b, float& d, float& vlos) 
{
  l=longitude(R,phi,R0);
  l=l*RAD;
  b=latitude(R,phi,z,l,R0);
  b=b*RAD;
  d=distance(R,z,l,R0,phi);
  vlos=v_los(v_R, v_phi, v_z, l, b, phi);
  vlos=vlos-v0_phi_R0*sin(l/RAD)*cos(b/RAD);
  hvlos=hvcvel(vlos, l, b, R0, v0_phi_R0, Rmax_hvcs, zmax_hvcs);
}


// OTHER FUNCTIONS

/* EXTERNAL VIEW
   phi=0 along the line of sight, then increases to the right
   when phi<270 phi>90 near side of the galaxy
 */
float projections::x_proj(float R, float phi)
{
  float x;
    x=R*sin(phi/RAD);
  return x;
}
float projections::y_proj(float z,float incl,float R,float phi)
{
  float y;
  y=z*sin((180-incl)/RAD)+R*cos((180-incl)/RAD)*cos(phi/RAD);
  return y;
}
float projections::v_proj(float v_R, float v_phi, float v_z, float incl, float phi)
{
  float v;
  v=-v_R*sin(incl/RAD)*cos((phi+180.)/RAD)+v_phi*sin(incl/RAD)*sin((phi)/RAD)-v_z*cos((180-incl)/RAD);
  return v;
}


/* INTERNAL VIEW
 */
float projections::longitude(float R, float phi, float R0)
{
  if (phi > 360){phi=phi-360;}
  float l=atan(sin(phi/RAD)/(R0/R-cos(phi/RAD)));
  if (phi < 180 && l < 0) {
    l=l+180/RAD;
  }
  if (phi > 180 && l > 0) {
    l=l-180/RAD;
  }
  if (phi == 180) {l=0;}
  if (phi == 0 && R <= R0) {l=0;}
  if (phi == 0 && R >= R0) {l=360/RAD;}
  return l;
}

float projections::latitude(float R, float phi, float z, float l, float R0)
{
  float b=atan(z*sin(l/RAD)/(R*sin(phi/RAD)));
  if (phi == 180) {b=atan(z/(R+R0));}
  if (phi == 0) {b=atan(z/fabs(R0-R));}
  return b;
}

float projections::distance(float R, float z, float l, float R0, float phi) 
{
  float d=0.;
  float cphi=phi;
	while (cphi>360) {cphi-=360;}
	while (cphi<0) {cphi+=360;}
  if (R < R0) 
  {
    if (cphi >= 0 && cphi < 90-asin(R/R0)*RAD){
      d=R0*(cos(l/RAD)-sqrt(fabs(R*R/R0/R0-sin(l/RAD)*sin(l/RAD))));
      // abs is there otherwise error if sqrt(-0.00000);
      d=sqrt(d*d+z*z);
    }
    if (cphi >= 90-asin(R/R0)*RAD && cphi <= 270+asin(R/R0)*RAD){ 
      d=R0*(cos(l/RAD)+sqrt(fabs(R*R/R0/R0-sin(l/RAD)*sin(l/RAD))));
      d=sqrt(d*d+z*z);
    }
    if (cphi > 270+asin(R/R0)*RAD && cphi <= 360){  
      d=R0*(cos(l/RAD)-sqrt(fabs(R*R/R0/R0-sin(l/RAD)*sin(l/RAD))));
      d=sqrt(d*d+z*z);
    }
  }
  if (R >= R0) {
    d=R0*(cos(l/RAD)+sqrt(fabs(R*R/R0/R0-sin(l/RAD)*sin(l/RAD))));
    d=sqrt(d*d+z*z);
  }
  if (d > 0) {
    return d;
  } else {
	  cout << "Warning: zero distance! - (R,phi,z) = ("<<R<<","<<phi<<","<<z<<")"<<endl;
    return 0;
  }
}

float projections::v_los(float v_R, float v_phi, float v_z, float l, float b, float phi)
{
  float vlos;
  vlos=v_R*cos((l+phi)/RAD)*cos(b/RAD)+v_phi*sin((l+phi)/RAD)*cos(b/RAD)-v_z*sin(b/RAD);
  return vlos;
}

/*
void projections::BuildInternalCube(Array3D<float>& cube, Array2D<float>& image,
			    int RAsize, int DECsize, int VELsize, 
			    float cdelt1, float cdelt2, float cdelt3)
{
  // Position of the Sun (for MW projection)
  float R0; findValue(inputfile, R0, "R0"); 
  // Rotation velocity at the sum
  float v0_phi_R0; findValue(inputfile, v0_phi_R0, "v0_phi_R0"); 
  // Cut off distance for closeby HVCs
  float distance_cut; findValue(inputfile, distance_cut, "dist_cut"); 
  float pixtopc;

  projections MW;
  
  MW.internal_view(R0,v0_phi_R0,l,b,d,vlos); // output: l(degree),b(degree),vlos;
  
  int xc=int(RAsize/2+l/fabs(cdelt1)+0.5);
  int yc=int(DECsize/2+b/cdelt2+0.5);
  int vc=int(VELsize/2-vlos/cdelt3+0.5);
  
  if (xc == RAsize) xc=0;
  //		      cout << l << " " << xc << " " <<endl;
  
  // Building the final cube
  if (vc < VELsize && vc >= 0) {
    if (d>distance_cut){
      pixtopc=fabs(cdelt1*cdelt2)/RAD/RAD*d*d*1.e6; 
      // pix -> pc^2
      outcube[xc][yc][vc]+=(cn_dot+cn_accr)
	*deltaR*deltaR*dt/dt_default/pixtopc;
      outcube_tot[xc][yc]+=(cn_dot+cn_accr)
	*deltaR*deltaR*dt/dt_default/pixtopc;
      // if the projected cloudsize is smaller than the 
      // pixel size
      if (d > RAD/sqrt(abs(cdelt1*cdelt2))*cloudsize/1000.){
	outcube[xc][yc][vc]+=(cn_dot+cn_accr)
	  *deltaR*deltaR*dt/dt_default/pixtopc;
	outcube_tot[xc][yc]+=(cn_dot+cn_accr)
	  *deltaR*deltaR*dt/dt_default/pixtopc;
	counts_HI+=(cn_dot+cn_accr)
	  *deltaR*deltaR*dt/dt_default;
	if (fabs(x3) > z_above) {
	  counts_HI_above+=(cn_dot+cn_accr)
	    *deltaR*deltaR*dt/dt_default;
	  above=1;
			    }
      }
      else {
	pixsize1=(int)(asin(cloudsize/d/1000.)*RAD/fabs(cdelt1)+.5);
	pixsize2=(int)(asin(cloudsize/d/1000.)*RAD/cdelt2+.5);
	for (int p=0; p<pixsize1*pixsize2; p++) {
	  seed=time(NULL)*((jj+1)*(p+1)*432143);
	  // randomization in each direction
	  pxc=(int)(xc+ran_gau(seed,pixsize1/3.));
	  pyc=(int)(yc+ran_gau(seed/3,pixsize2/3.));
	  pvc=(int)(vc+ran_gau(seed/2,sigma_v/cdelt3));
	  if (pxc >= RAsize) pxc=pxc-RAsize;
	  if (pxc < 0) pxc=pxc+RAsize; 
	  if (pyc >= DECsize) pyc=pyc-DECsize;
	  if (pyc < 0) pyc=pyc+DECsize; 
	  if (pvc >= VELsize) pyc=VELsize;
	  if (pvc < 0) pyc=0; 
	  outcube[pxc][pyc][pvc]+=(cn_dot+cn_accr)	
	    *deltaR*deltaR*dt/dt_default/pixtopc
	    /pixsize1/pixsize2*(ran_gau(seed,.3)+1);
	  outcube_tot[pxc][pyc]+=(cn_dot+cn_accr)
	    *deltaR*deltaR*dt/dt_default/pixtopc
	    *fabs(cdelt3)
	    /pixsize1/pixsize2*(ran_gau(seed,.3)+1);
	}
      }
    }
  }
} 
*/

/* HVCS VIEW
   Transforms the velocity along the line of sight (vlos) in velocity
   difference from a simple model of galactic differential rotation 
   given by MODELGAL (see Wakker 1991).   
*/
float projections::hvcvel(float tvlos, float tl, float tb, float R0, float v0_phi_R0, float Rmax_hvcs, float zmax_hvcs)
{
  if (tl < 0) {tl=tl+360;}
  float vlosmin=modelgal(tl,tb,-1, R0, v0_phi_R0, Rmax_hvcs, zmax_hvcs);
  float vlosmax=modelgal(tl,tb,+1, R0, v0_phi_R0, Rmax_hvcs, zmax_hvcs);
  //cout << vlosmin << " " << tvlos << " " << vlosmax << "\n";
  if (tvlos < vlosmax && tvlos > vlosmin) 
    {tvlos=0;
    }
  else
    {
      if (tvlos > 0) {tvlos=tvlos-vlosmax;}
      if (tvlos < 0) {tvlos=tvlos-vlosmin;}
    }
  return tvlos;
}


float projections::modelgal(float tl, float tb, int minmax, float R0, float v0, float Rmax_hvcs, float zmax_hvcs) {
/* MODELGAL
   Calculates a model galaxy (MW type) in differential rotation with a flat
   rotation curve of v0_phi_R0, maximum R and z of 26 and 4 kpc respectively.
   The galaxy is seen from the position of the Sun at distance R0 from the
   centre.
   Return vlosmin and vlosmax (depending on value of parameter minmax: -1
   or +1) that are the minumum and maximum velocities permitted in such a 
   model. Clouds belove and above this boundary are anomalous or high 
   velocity clouds.
*/
  tb=fabs(tb);
  
  float phimax1=0.0;
  phimax1=fabs(asin(R0/Rmax_hvcs*sin((tl-180)/RAD))*RAD-tl+180);
  // why is tl-180 necessary???
  //  phimax1=fabs(asin(R0/Rmax_hvcs*sin((tl)/RAD))*RAD-tl+180);

  float bstar;
  //  bstar=fabs(atan(zmax_hvcs/Rmax_hvcs*sin(phimax1/RAD)/sin((tl)/RAD)))*RAD;
  //  bstar=fabs(atan(zmax_hvcs/(R0*cos(tl/RAD)+sqrt(Rmax_hvcs*Rmax_hvcs-R0*R0*sin(tl/RAD)*sin(tl/RAD)))))*RAD;
  bstar=fabs(atan2((double) zmax_hvcs,(R0*cos(tl/RAD)+sqrt(Rmax_hvcs*Rmax_hvcs-R0*R0*sin(tl/RAD)*sin(tl/RAD)))))*RAD;

  float phimax=0.0;
  float phimax2=0.0;
  if (fabs(tb) <= bstar) {
    phimax=phimax1;
  } else {
    phimax2=atan(sin(tl/RAD)/((R0/zmax_hvcs)*tan(tb/RAD)-cos(tl/RAD)))*RAD;
    if (tl > 180) {phimax2=-phimax2;}
    if (phimax2 < 0) {phimax2=phimax2+180;}
    phimax=phimax2;
  }  
  
  float vlosmin=0;
  float vlosmax=0;
  float tphi=0.0;

  if (fabs(phimax) < 1) {tphi=phimax;}

  /*
    if (tl <= 90) {
    vlosmin=v0*cos(tb/RAD);//*(sin((tl+phimax)/RAD)-sin(tl/RAD));
    vlosmax=0.;
    } else {
    if (tl <= 180) {
    vlosmin=0.;
    vlosmax=v0*cos(tb/RAD);//*(1+sin(tl/RAD));
    } else {
    if (tl <= 270) {
    vlosmin=v0*cos(tb/RAD);//*(1-sin(tl/RAD));
    vlosmax=0.;
    } else {
    if (tl <= 360) {
    vlosmin=0.;
    vlosmax=v0*cos(tb/RAD);//*(sin((tl+phimax)/RAD)-sin(tl/RAD));
    }
    }
    }
    }
  */

  while (tphi < phimax && tphi <= 180.0) {
    float vlos=0.0;
    if (tl < 180) {
      vlos=v0*cos(tb/RAD)*(sin((tl+tphi)/RAD)-sin(tl/RAD));
    }
    if (tl > 180) {
      vlos=-v0*cos(tb/RAD)*(sin((-tl+tphi)/RAD)-sin((-tl)/RAD));
    }
    if (vlos < vlosmin) {vlosmin=vlos;}
    if (vlos > vlosmax) {vlosmax=vlos;}
    tphi=tphi+1.0;
  }

  float retvlos;
  if (minmax == -1) {retvlos=vlosmin;}
  if (minmax == +1) {retvlos=vlosmax;}
  return retvlos;
}
