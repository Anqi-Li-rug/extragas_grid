#include "dynamics.h"
#include <string.h>
extern const double PI;//=3.14159265359;
extern const double RAD=180./PI;
extern const double CONV1=1.02269012107e-3; // km/s -> kpc/Myr
extern long int seed;

//template <class T>(vz>0 && fabs(z)>accr_z_in) || (vz<0 fabs(z)>accr_z_fin)
dynamics::dynamics(Array2D<double>& F_R, 
		   Array2D<double>& F_z,
		   Array1D<double> Ri, 
		   Array1D<double> zk)
{
  F_R_ik=F_R;
  F_z_ik=F_z;
  R_i=Ri;
  z_k=zk;
}

dynamics::dynamics(Array2D<double>& F_R, 
		   Array2D<double>& F_z,
		   Array2D<double>& pot,
		   Array1D<double> Ri, 
		   Array1D<double> zk,
		   int dim1,
		   int dim2,
		   char eff[],
		   float KD,
		   float v0,
		   Array2D<double>& rho_hot,
		   double relative_v,
		   double sigma_accr,
		   double cloud_mass,
		   double cloud_radius,
		   double corona_density,
		   double alpha_input	
			)
{
  F_R_ik=F_R;
  F_z_ik=F_z;
  pot_ik=pot;
  R_i=Ri;
  z_k=zk;
  subdim1=dim1;
  subdim2=dim2;
  strcpy(effects,eff);
  rho_hot_ik=rho_hot;
  K_D=KD;
  v0_hot=v0;
  Vrel = relative_v;
  sigmaV = sigma_accr;
  M_cloud = cloud_mass;
  R_cloud = cloud_radius;
  rho_corona = corona_density;
  alpha_accr = alpha_input*1e-3; //in Myr^-1	
  drag_coeff = 1.51e22*PI*R_cloud*R_cloud*rho_corona/M_cloud;
}

/*
dynamics::dynamics(Array2D<float>& F_R, 
		   Array2D<float>& F_z,
		   Array1D<float> Ri, 
		   Array1D<float> zi)
{
}
*/

//template <class T>
dynamics::~dynamics() {}

//template <class T>
void dynamics::rk4(double time,
		   double d_t,
		   double& x, 
		   double& y, 
		   double& z, 
		   double& v_x, 
		   double& v_y, 
		   double& v_z,
		   float n)
  //  PERFORMS RUNGE-KUTTA INTEGRATION to 4th order accurancy
  //  INPUTS: dt=increment
  //  READS coordinates at time t, RETURNS coordinates at time t+dt
{	

  dt=d_t;
  n_dot=n;
  double dx1,dx2,dx3,dx4,dy1,dy2,dy3,dy4,dz1,dz2,dz3,dz4;
  double dvx1,dvx2,dvx3,dvx4,dvy1,dvy2,dvy3,dvy4,dvz1,dvz2,dvz3,dvz4;

  //  printf("%.1f %.1f %.3f %.1f %.1f %.1f %.1f %.1f \n", x, y, z, v_x/CONV1, v_y/CONV1, v_z/CONV1, vR/CONV1, vphi/CONV1);
  //  printf("%.1f %.1f %.3f %.1f %.1f %.1f %.1f %.1f \n", x, y, z, v_x/CONV1, v_y/CONV1, v_z/CONV1, dt, n);

  dx1=v_x*dt;
  dy1=v_y*dt;
  dz1=v_z*dt;
  getForce(time,x,y,z,v_x,v_y,v_z);
  dvx1=-g_x*dt;
  dvy1=-g_y*dt;
  dvz1=-g_z*dt;

  //  printf("%.1f %.1f %.3f %.1f %.1f %.1f %.1f %.1f \n", x, y, z, v_x/CONV1, v_y/CONV1, v_z/CONV1, dt, n);

  dx2=(v_x+dvx1/2)*dt;
  dy2=(v_y+dvy1/2)*dt;
  dz2=(v_z+dvz1/2)*dt;
  getForce(time,x+dx1/2, y+dy1/2, z+dz1/2, v_x+dvx1/2, v_y+dvy1/2, v_z+dvz1/2);
  dvx2=-g_x*dt;
  dvy2=-g_y*dt;
  dvz2=-g_z*dt;

  //  printf("%.1f %.1f %.3f %.1f %.1f %.1f %.1f %.1f \n", x, y, z, v_x/CONV1, v_y/CONV1, v_z/CONV1, dt, n);

  dx3=(v_x+dvx2/2)*dt;
  dy3=(v_y+dvy2/2)*dt;
  dz3=(v_z+dvz2/2)*dt;
  getForce(time,x+dx2/2,y+dy2/2,z+dz2/2,v_x+dvx2/2, v_y+dvy2/2, v_z+dvz2/2);
  dvx3=-g_x*dt;
  dvy3=-g_y*dt;
  dvz3=-g_z*dt;

  //  printf("%.1f %.1f %.3f %.1f %.1f %.1f %.1f %.1f \n", x, y, z, v_x/CONV1, v_y/CONV1, v_z/CONV1, dt, n);

  dx4=(v_x+dvx3)*dt;
  dy4=(v_y+dvy3)*dt;
  dz4=(v_z+dvz3)*dt;
  getForce(time,x+dx3,y+dy3,z+dz3,v_x+dvx3, v_y+dvy3, v_z+dvz3);
  dvx4=-g_x*dt;
  dvy4=-g_y*dt;
  dvz4=-g_z*dt;

  //  printf("%.1f %.1f %.3f %.1f %.1f %.1f %.1f %.1f \n", x, y, z, v_x/CONV1, v_y/CONV1, v_z/CONV1, dt, n);

  x=x+(dx1+2*dx2+2*dx3+dx4)/6.;
  y=y+(dy1+2*dy2+2*dy3+dy4)/6.;
  z=z+(dz1+2*dz2+2*dz3+dz4)/6.;
  v_x=v_x+(dvx1+2*dvx2+2*dvx3+dvx4)/6.;
  v_y=v_y+(dvy1+2*dvy2+2*dvy3+dvy4)/6.;
  v_z=v_z+(dvz1+2*dvz2+2*dvz3+dvz4)/6.;
	
  //  printf("%.1f %.1f %.3f %.1f %.1f %.1f %.1f %.1f \n", x, y, z, v_x/CONV1, v_y/CONV1, v_z/CONV1, dt, n);
  //  cout << endl;

  //  printf("%.1f %.1f %.3f %.1f \n", dvx1/CONV1, dvx2/CONV1, dvx3/CONV1, dvx4/CONV1);
  //  printf("%.1f %.1f %.3f %.1f %.1f %.1f \n", x, y, z, v_x/CONV1, v_y/CONV1, v_z/CONV1);
}

void dynamics::getForce(double time, double x, double y, double z,
			double vx, double vy, double vz)
  // GIVES THE VALUE OF THE ACCELERATION IN (CARTESIAN) POSITION (x,y,z)
  // N.B.
  // g_z > 0 means attraction towards the disc
{
  R=sqrt(x*x+y*y);
  phi=atan2(y,x);
  if (phi <= 0) {phi=phi+2*PI;}
  
  g_x=getValue(F_R_ik,R_i,z_k,R,fabs(z),1);
  g_y=g_x*sin(phi);
  g_x*=cos(phi);
  if (z == fabs(z)) {
    g_z=getValue(F_z_ik,R_i,z_k,R,z,1);
  } else {
    g_z=-getValue(F_z_ik,R_i,z_k,R,fabs(z),1);
  }
  
  //  cout << R << " " << fabs(z) << " " << g_x << " " << g_z << endl;

  if (strcmp(effects,"drag")==0)
    DragForce(x,y,z,vx,vy,vz,g_x,g_y,g_z);

  //  cout << R << " " << g_x << " " << g_y << " " << g_z << endl;

  if (strcmp(effects,"accretion")==0)
    AccrForce(vx,vy,vz,g_x,g_y,g_z);
	//  cout << R << " " << g_x << " " << g_y << " " << g_z << endl;
	
  if (strcmp(effects,"accrdrag")==0)
    AccrDragForce(time,z,vx,vy,vz,g_x,g_y,g_z);
  if (strcmp(effects,"accrnodrag")==0)
    AccrNoDragForce(time,z,vx,vy,vz,g_x,g_y,g_z);
    //cout << R << " " << g_x << " " << g_y << " " << g_z << endl;
}


/*
//template <class T>
void dynamics::getForceMem(double x, double y, double z, 
			   double vx, double vy, double vz)
  // GIVES THE VALUE OF THE ACCELERATION IN (CARTESIAN) POSITION (x,y,z)
  // remembers the previous subarray and avoids to extract it each time
{
  R=sqrt(x*x+y*y);
  phi=atan2(y,x);
  if (phi <= 0) {phi=phi+2*PI;}

  // N.B. *** REPLACE -1 with SIZE OF SUBARRAY ***

  if ((R < R_i[R_i.dim()-subdim1]) && (fabs(z) <z_k[z_k.dim()-subdim2])) {
    if (i0 >= 0 && i0 < R_i.dim()-subdim1 &&
	j0 >= 0 && j0 < z_k.dim()-subdim2 &&
	ArraysCreated) {
      if (R >= R_i[i0] && R < R_i[i0+1] &&
	  fabs(z) >= z_k[j0] && fabs(z) < z_k[j0+1] &&
	  ArraysCreated) {
	g_x=getValueSubarray(subF_R_ik,subR_i,subz_k,R,fabs(z));
	g_y=g_x*sin(phi);
	g_x*=cos(phi);
	if (z == fabs(z)) {
	  g_z=getValueSubarray(subF_z_ik,subR_i,subz_k,R,z);
	} else {
	  g_z=-getValueSubarray(subF_z_ik,subR_i,subz_k,R,fabs(z));
	}
	//      cout << R << " " << z << " " <<  R_i[i0] << " " << z_k[j0] << " " << i0 << " " << j0 << " " << endl;
	//      print(subF_R_ik);
	if (strcmp(effects,"drag")==0) {
	  DragForceMem(vx,vy,vz,g_x,g_y,g_z);
	}
	if (strcmp(effects,"accretion")==0)
	  AccrForce(vx,vy,vz,g_x,g_y,g_z);
      }
      else {
	g_x=getValue(F_R_ik,R_i,z_k,R,fabs(z),1);
	g_y=g_x*sin(phi);
	g_x*=cos(phi);
	if (z == fabs(z)) {
	  g_z=getValue(F_z_ik,R_i,z_k,R,z,1);
	} else {
	  g_z=-getValue(F_z_ik,R_i,z_k,R,fabs(z),1);
	}
	subF_R_ik=subarray(F_R_ik, R_i, z_k, R,fabs(z), subdim1, subdim2);
	subF_z_ik=subarray(F_z_ik, R_i, z_k, R,fabs(z), subdim1, subdim2);
	subR_i=subarray(R_i, R, 2);
	subz_k=subarray(z_k, fabs(z), 2);
	ArraysCreated=true;
	i0=prev_i0;
	j0=prev_j0;
	//      cout << "IN 1" << endl;
	//      cout << R << " " << z << " " <<  R_i[i0] << " " << z_k[j0] << " " << i0 << " " << j0 << " " << endl;
	//      print(subF_R_ik);
	if (strcmp(effects,"drag")==0) {
	  DragForce(x,y,z,vx,vy,vz,g_x,g_y,g_z);
	  sub_rho_ik=subarray(rho_hot_ik, R_i, z_k, R,fabs(z), subdim1, subdim2);
	}
	if (strcmp(effects,"accretion")==0)
	  AccrForce(vx,vy,vz,g_x,g_y,g_z);
      }
    }
    else {
      g_x=getValue(F_R_ik,R_i,z_k,R,fabs(z),1);
      g_y=g_x*sin(phi);
      g_x*=cos(phi);
      if (z == fabs(z)) {
	g_z=getValue(F_z_ik,R_i,z_k,R,z,1);
      } else {
	g_z=-getValue(F_z_ik,R_i,z_k,R,fabs(z),1);
      }
      subF_R_ik=subarray(F_R_ik, R_i, z_k, R,fabs(z), subdim1, subdim2);
      subF_z_ik=subarray(F_z_ik, R_i, z_k, R,fabs(z), subdim2, subdim2);
      subR_i=subarray(R_i, R, 2);
      subz_k=subarray(z_k, fabs(z), 2);
      i0=prev_i0;
      j0=prev_j0;
      ArraysCreated=true;
      if (strcmp(effects,"drag")==0) {
      	DragForce(x,y,z,vx,vy,vz,g_x,g_y,g_z);
      	sub_rho_ik=subarray(rho_hot_ik, R_i, z_k, R,fabs(z),subdim1, subdim2);
      }
      if (strcmp(effects,"accretion")==0)
	AccrForce(vx,vy,vz,g_x,g_y,g_z);
    }
  }
  else {
    cout << "***Warning: Particle outside (inner?) matrix... skipped" << endl;
  }

}
*/
 
void dynamics::DragForce(double x, double y, double z, 
			 double vx, 
			 double vy, 
			 double vz, 
			 double& g_x, 
			 double& g_y, 
			 double& g_z)
{
  double vrel, K_drag, crho_hot;
  R=sqrt(x*x+y*y);


  vrel_x=vx+v0_hot*sin(phi/RAD)*CONV1;
  vrel_y=vy-v0_hot*cos(phi/RAD)*CONV1;
  vrel_z=vz;
  vrel=sqrt(vrel_x*vrel_x+vrel_y*vrel_y+vrel_z*vrel_z);
  
  K_drag=.5*PI*K_D/1.e6; // kpc^2/Mo
  //      K_drag=.5*K0_drag*pi*(D_cloud/2.)*(D_cloud/2.)/m_cloud/1e6; // kpc^2/Mo

  if (R > (R_i[1]-R_i[0])) {
    crho_hot=getValue(rho_hot_ik,R_i,z_k,R,fabs(z),1);
  }
  else {
    // to avoid very large densities
    crho_hot=getValue(rho_hot_ik,R_i,z_k,(R_i[1]-R_i[0]),fabs(z),1);
  }
  
  g_x+=K_drag*crho_hot*vrel_x*vrel;
  g_y+=K_drag*crho_hot*vrel_y*vrel;
  g_z+=K_drag*crho_hot*vrel_z*vrel;  

  //  cout << R << " " << fabs(z) << " " << K_drag << " " << crho_hot << " " << vrel_x << " " << g_x << " " << g_z << endl;
  //      printf("%.1f %.1f %.3f %.4e %.1f %.1f %.1f %.1f %.1f %.1f %f %f %f %f %f %f \n", cR, cphi, sign*cz, crho_hot, vR/conv1, vphi/conv1, vz/conv1, vrel_R/conv1, vrel_phi/conv1, vrel_z/conv1, g_x, g_y, g_z, cF_R*cos(cphi/rad), cF_R*sin(cphi/rad), cF_z);

}

void dynamics::DragForceMem(double vx, 
			    double vy, 
			    double vz, 
			    double& g_x, 
			    double& g_y, 
			    double& g_z)
{
  vrel_x=vx+v0_hot*sin(phi/RAD)*CONV1;
  vrel_y=vy-v0_hot*cos(phi/RAD)*CONV1;
  vrel_z=vz;
  
  double K_drag=.5*PI*K_D/1.e6; // kpc^2/Mo

  double crho_hot=getValueSubarray(sub_rho_ik,subR_i,subz_k,R,fabs(z));
  g_x+=K_drag*crho_hot*vrel_x*fabs(vrel_x);
  g_y+=K_drag*crho_hot*vrel_y*fabs(vrel_y);
  g_z+=K_drag*crho_hot*vrel_z*fabs(vrel_z);

  //  cout << K_drag << " " << crho_hot << " " << vrel_x << " " << g_x << endl;
}

void dynamics::AccrForce(double vx, 
			 double vy, 
			 double vz, 
			 double& g_x, 
			 double& g_y, 
			 double& g_z)
{
  // IN THE FRAME WHERE THE CLOUD IS AT REST
  // CORRECTIONS TO THE GRAVITATIONAL ACCELERATIONS
    g_x+=-(dn_accr)/(n_dot-dn_accr*dt)*(v_accr_x-vx);
	g_y+=-(dn_accr)/(n_dot-dn_accr*dt)*(v_accr_y-vy);
	g_z+=-(dn_accr)/(n_dot-dn_accr*dt)*(v_accr_z-vz);
  /*  
      cout << " " << dt << " " << v_accr_x << " " 
      << vx << " " 
      << g_x << " " 
      << dn_accr << " " 
      << n_dot << endl;
      
      cout << dn_accr << " " << (dn_accr)/(n_dot+dn_accr) << " " 
      << (v_accr_x-vx)/CONV1 << " " 
      << (dn_accr)/(n_dot+dn_accr)*(v_accr_x-vx)/CONV1
      << endl;
      cout << dn_accr << " " << (dn_accr)/(n_dot+dn_accr) << " " 
      << (v_accr_y-vy)/CONV1 << " " 
      << (dn_accr)/(n_dot+dn_accr)*(v_accr_y-vy)/CONV1
      << endl;
      cout << "a " << g_x << " " << g_y << " " << g_z << endl;
  */
}

void dynamics::AccrDragForce(
						double time,
						double z,
						double vx, 
						double vy, 
						double vz, 
						double& g_x, 
						double& g_y, 
						double& g_z
							 )
{
	// IN THE FRAME WHERE THE CLOUD IS AT REST
	// CORRECTIONS TO THE GRAVITATIONAL ACCELERATIONS
	if((vz>0 && fabs(z)>accr_z_in) || (vz<0 && fabs(z)>accr_z_fin))
	{
		//NEW VERSION		
		vrel_x = (v_accr_x-vx);
		vrel_y = (v_accr_y-vy);
		vrel_z = (v_accr_z-vz);
		vrel_tot = sqrt(vrel_x*vrel_x + vrel_y*vrel_y + vrel_z*vrel_z)/CONV1;
		
	//deceleration produced by accretion
		
		//ad-hoc simulation-like accretion
		//double g_accr = d_accr_sim(time+tdelay)/(1.+accr_sim(time+tdelay)-accr_sim(tdelay));
		
		//exp accretion
		double g_accr  = alpha_accr;		
		//power-law accretion
		//double g_accr = (beta/alpha)*pow(time/alpha,beta-1.)/(1.+pow(time/alpha,beta));
		
		//deceleration produced by drag
		double g_drag  = (vrel_tot-Vrel)*drag_coeff;
		if (g_drag<0){
			g_drag=0;
		}
		double g_hydro = g_accr+g_drag;
		g_x+=-g_hydro*vrel_x;
		g_y+=-g_hydro*vrel_y;
		g_z+=-g_hydro*vrel_z;
	}
}



void dynamics::AccrNoDragForce(
						double time,
						double z,
						double vx, 
						double vy, 
						double vz, 
						double& g_x, 
						double& g_y, 
						double& g_z
							 )
{
	// IN THE FRAME WHERE THE CLOUD IS AT REST
	// CORRECTIONS TO THE GRAVITATIONAL ACCELERATIONS
	if((vz>0 && fabs(z)>accr_z_in) || (vz<0 && fabs(z)>accr_z_fin))
	{
		//NEW VERSION		
		vrel_x = (v_accr_x-vx);
		vrel_y = (v_accr_y-vy);
		vrel_z = (v_accr_z-vz);
		vrel_tot = sqrt(vrel_x*vrel_x + vrel_y*vrel_y + vrel_z*vrel_z)/CONV1;
		
	//deceleration produced by accretion
		
		//ad-hoc simulation-like accretion
		//double g_accr = d_accr_sim(time+tdelay)/(1.+accr_sim(time+tdelay)-accr_sim(tdelay));
		
		//exp accretion
		double g_accr  = alpha_accr;		
		//power-law accretion
		//double g_accr = (beta/alpha)*pow(time/alpha,beta-1.)/(1.+pow(time/alpha,beta));
		
		//deceleration produced by drag
		//double g_drag  = vrel_tot*drag_coeff;
		double g_drag=0.;
                double g_hydro = g_accr+g_drag;
				
		g_x+=-g_hydro*vrel_x;
		g_y+=-g_hydro*vrel_y;
		g_z+=-g_hydro*vrel_z;
	}
}





float dynamics::accretion(double x, double y, double z, double vz,
						  char accr_dens[],
						  char accr_vel[],
						  float m0,
						  float time)
{
	//comment this! 
	//accr_norm = dyn_accr_norm;
	R=sqrt(x*x+y*y);
	v_esc=sqrt(fabs(2.*getValue(pot_ik,R_i,z_k,R,fabs(z),1)))/CONV1;
	int i_ran=1;
	double x_1=x,x_2=y,x_3=z,x_4=v_esc,x_5=v_esc,x_6=v_esc;
	
	if (strcoll(accr_dens,"constant")==0) {
		
		//ad-hoc simulation-like 
		//dn_accr = m0*d_accr_sim(time+tdelay);
		
		//exponential
		dn_accr=alpha_accr*m0*exp(alpha_accr*time);
		
		//power law
		//dn_accr=(beta/alpha)*m0*pow(time/alpha,beta-1.);
	}
	if (strcoll(accr_dens,"linear_z")==0) {
		dn_accr=alpha_accr*m0*exp(alpha_accr*time)*z;
	}
	if (strcoll(accr_dens,"squareR")==0) {
		dn_accr=alpha_accr*m0*exp(alpha_accr*time)/R/R*m0;
	}
	if (strcoll(accr_dens,"custom")==0) {
		dn_accr=alpha_accr*m0*exp(alpha_accr*time)/(R/5.)*m0;
	}
	if (strcoll(accr_vel,"static")==0) {
		
		cartospher(x_1, x_2, x_3, x_4, x_5, x_6);      
		while (sqrt(x_4*x_4+x_5*x_5+x_6*x_6) >= v_esc) {
			x_4=ran_gau(seed,sigmaV);
			
			x_5=ran_gau(seed,sigmaV);
			
			x_6=ran_gau(seed,sigmaV);
			i_ran+=1;
		}
		sphertocar(x_1,x_2,x_3,x_4,x_5,x_6);
	}
	
	if (strcoll(accr_vel,"radial")==0) {
		
		cartospher(x_1, x_2, x_3, x_4, x_5, x_6);      
		while (sqrt(x_4*x_4+x_5*x_5+x_6*x_6) >= v_esc) {
			x_4=-(fabs(ran0(&seed)*v_esc));
			
			x_5=ran_gau(seed,sigmaV);
			
			x_6=ran_gau(seed,sigmaV);
			i_ran+=1;
		}
		sphertocar(x_1,x_2,x_3,x_4,x_5,x_6);
	}
	
	if (strcoll(accr_vel,"radial_low")==0) {
		
		cartospher(x_1, x_2, x_3, x_4, x_5, x_6);
		while (sqrt(x_4*x_4+x_5*x_5+x_6*x_6) >= v_esc) {
			x_4=-(fabs(ran_gau(seed,v_esc/2)));
			
			x_5=ran_gau(seed,sigmaV);
			
			x_6=ran_gau(seed,sigmaV);
			i_ran+=1;
		}
		sphertocar(x_1,x_2,x_3,x_4,x_5,x_6);
	}
	
	if (strcoll(accr_vel,"polar")==0) {
		while (sqrt(x_4*x_4+x_5*x_5+x_6*x_6) >= v_esc) {
			x_4=ran_gau(seed,sigmaV);
			
			x_5=ran_gau(seed,sigmaV);
			
			x_6=-(fabs(ran0(&seed)*v_esc));
			//	x_6=-(fabs(ran_gau(seed,v_esc/2)));
			i_ran+=1;
		}
	}
	
	if (strcoll(accr_vel,"random")==0) {
		while (sqrt(x_4*x_4+x_5*x_5+x_6*x_6) >= v_esc) {
			x_4=ran_gau(seed,v_esc/3);
			
			x_5=ran_gau(seed,v_esc/3);
			
			x_6=ran_gau(seed,v_esc/3);
			i_ran+=1;
		}
	}
	
	if (strcoll(accr_vel,"rotating")==0) {
		cartocyl(x_1, x_2, x_3, x_4, x_5, x_6);      
		while (sqrt(x_4*x_4+x_5*x_5+x_6*x_6) >= v_esc) {
			x_4=ran_gau(seed,sigmaV);
			
			x_5=(fabs(ran0(&seed)*sqrt(R*getValue(F_R_ik,R_i,z_k,R,0.,1))/CONV1));
			
			x_6=ran_gau(seed,sigmaV);
			i_ran+=1;
		}
		cyltocar(x_1, x_2, x_3, x_4, x_5, x_6);      
	}
	
	//rotating with lag Vrel with respect to the disk	  
	if (strcoll(accr_vel,"slowrot")==0) {
		cartocyl(x_1, x_2, x_3, x_4, x_5, x_6);      
		while (sqrt(x_4*x_4+x_5*x_5+x_6*x_6) >= v_esc) {
			x_4=0;
			x_5=(sqrt(R*getValue(F_R_ik,R_i,z_k,R,fabs(x_3),1))/CONV1) - Vrel;
			if(x_5<0){x_5=0;}
			x_6=0;
		}
		cyltocar(x_1, x_2, x_3, x_4, x_5, x_6);      
	}
	
	//rotating with lag Vrel up to |z|=height_cut, then static  
	if (strcoll(accr_vel,"slowrot_cut")==0) {
		cartocyl(x_1, x_2, x_3, x_4, x_5, x_6);      
		while (sqrt(x_4*x_4+x_5*x_5+x_6*x_6) >= v_esc) {
			x_4=ran_gau(seed,sigmaV);
			
			if(fabs(x_3)<3.5) {x_5=(sqrt(R*getValue(F_R_ik,R_i,z_k,R,fabs(x_3),1))/CONV1) - Vrel;}
			else {x_5 = ran_gau(seed,sigmaV);}
			if(x_5<0){x_5=0;}//just in case
			
			x_6=ran_gau(seed,sigmaV);
			
			i_ran+=1;
		}
		//cout<<"vR="<<x_4<<"  vphi="<<x_5<<"  vz="<<x_6<<endl;
		cyltocar(x_1, x_2, x_3, x_4, x_5, x_6);      
	}
	
	// polar + angular momentum (as the disk)
	if (strcoll(accr_vel,"rotpol")==0) {
		cartocyl(x_1, x_2, x_3, x_4, x_5, x_6);      
		while (sqrt(x_4*x_4+x_5*x_5+x_6*x_6) >= v_esc) {
			x_4=ran_gau(seed,sigmaV);
			
			x_5=fabs(sqrt(R*getValue(F_R_ik,R_i,z_k,R,0.,1))/CONV1)/2.;
			//	cout << x_5 << endl;
			
			x_6=-(fabs(ran0(&seed)*v_esc));
			i_ran+=1;
		}
		cyltocar(x_1, x_2, x_3, x_4, x_5, x_6);      
	}
	
	v_accr_x=x_4*CONV1;
	v_accr_y=x_5*CONV1;
	v_accr_z=x_6*CONV1;
	//    dn_accr=dn_accr*n_dot_dot;
	
	if ((vz>0 && fabs(z)>accr_z_in) || (vz<0 && fabs(z)>accr_z_fin)) {return dn_accr;}
	else
	{
		dn_accr = 0;
		return 0;
	}
}

float dynamics::accretion(double x, double y, double z, double vz,
			  Array2D<double>& cold_rho_ik,
			  Array2D<double>& cold_v_R_ik,
			  Array2D<double>& cold_v_phi_ik,
			  Array2D<double>& cold_v_z_ik,
			  Array1D<double> R_i,
			  Array1D<double> z_k,
			  float n_dot_dot)
{
	//comment this! 
	//accr_norm = dyn_accr_norm;
	
  double x_1=x,x_2=y,x_3=z,x_4=0,x_5=0,x_6=0;

  cartocyl(x_1, x_2, x_3, x_4, x_5, x_6);

  int i=0, k=0;
  
  if ((vz>0 && fabs(z)>accr_z_in) || (vz<0 && fabs(z)>accr_z_fin)) {
    while (x_1 > R_i[i])
      i++;
    while (x_3 > z_k[k])
      k++;

    x_4=cold_v_R_ik[i][k];
    x_5=cold_v_phi_ik[i][k];
    x_6=cold_v_z_ik[i][k];
    
    cyltocar(x_1, x_2, x_3, x_4, x_5, x_6);

    v_accr_x=x_4*CONV1;
    v_accr_y=x_5*CONV1;
    v_accr_z=x_6*CONV1;
    dn_accr=cold_rho_ik[i][k];
    //    cout << i << " " << k << " " << x << " " 
    //	 << R_i[i] << " " << y << " " << z_k[k] << endl;
    
    dn_accr=dn_accr*n_dot;
    return dn_accr;
  }
  else {
    dn_accr=0;
    return 0;
  }
}
