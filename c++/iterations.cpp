#include "iterations.h"
#include "galFunctions.h"
#include "fitsUtils.h"
#include "findValue.h"
#include "extragasPlots.h"
#include "extragasFits.h"
#include <string.h>
#include <stdio.h>
const double DEGTORAD=0.01745329252;

extern char inputfile[];
extern const double RAD;

inline double Rtod_logspiral (double R, double pitch_deg)
{
	return R/cos(atan(1./tan(pitch_deg*DEGTORAD)));
}
inline double dtoR_logspiral (double d, double pitch_deg)
{
	return d*cos(atan(1./tan(pitch_deg*DEGTORAD)));
}

iterations::iterations(char param1[], char param2[], char param3[],
		       int &niter1, int &niter2, int &niter3,
		       float &param1_0, float &delta_param1,
		       float &param2_0, float &delta_param2,
			   float &param3_0, float &delta_param3)
{
  /*
    DECLARATIONS
  */
  findValue(inputfile, RAsize, "RAsize"); 
  findValue(inputfile, DECsize, "DECsize"); 
  findValue(inputfile, VELsize, "VELsize"); 
  findValue(inputfile, x_shift, "x_shift"); 
  findValue(inputfile, y_shift, "y_shift"); 
  findValue(inputfile, cdelt1, "cdelt1"); 
  findValue(inputfile, cdelt2, "cdelt2"); 
  float distance; findValue(inputfile, distance, "distance"); 
  /* conversion between kpc and degrees */
  kpctograd=1./distance/1.e3*RAD;
  float h_v_k1; findValue(inputfile, h_v_k1, "h_v_k1"); 
  float h_v_k2; findValue(inputfile, h_v_k2, "h_v_k2"); 
  float delta_h_v_k; findValue(inputfile, delta_h_v_k, "delta_h_v_k"); 
  float M_halo1; findValue(inputfile, M_halo1, "M_halo1");  
  float M_halo2; findValue(inputfile, M_halo2, "M_halo2"); 
  float delta_M; findValue(inputfile, delta_M, "delta_M"); 
  float gammaSF1; findValue(inputfile, gammaSF1, "gammaSF1"); 
  float gammaSF2; findValue(inputfile, gammaSF2, "gammaSF2"); 
  float delta_gSF; findValue(inputfile, delta_gSF, "delta_gSF");
  float showt_min1; findValue(inputfile, showt_min1, "showt_min1"); 
  float showt_min2; findValue(inputfile, showt_min2, "showt_min2");
  float delta_showt_min; findValue(inputfile, delta_showt_min, "delta_showt_min");
  float showt2_min1; findValue(inputfile, showt2_min1, "showt2_min1"); 
  float showt2_min2; findValue(inputfile, showt2_min2, "showt2_min2");
  float delta_showt2_min; findValue(inputfile, delta_showt2_min, "delta_showt2_min");  
  float showt_max1; findValue(inputfile, showt_max1, "showt_max1"); 
  float showt_max2; findValue(inputfile, showt_max2, "showt_max2");
  float delta_showt_max; findValue(inputfile, delta_showt_max, "delta_showt_max");
  float kick_angle1; findValue(inputfile, kick_angle1, "kick_angle1"); 
  float kick_angle2; findValue(inputfile, kick_angle2, "kick_angle2");
  float delta_kick_angle; findValue(inputfile, delta_kick_angle, "delta_kick_angle");
  float K_D1; findValue(inputfile, K_D1, "K_D1"); 
  float K_D2; findValue(inputfile, K_D2, "K_D2"); 
  float delta_K_D; findValue(inputfile, delta_K_D, "delta_K_D"); 
  float v0_hot1; findValue(inputfile, v0_hot1, "v0_hot1"); 
  float v0_hot2; findValue(inputfile, v0_hot2, "v0_hot2"); 
  float delta_v0_hot; findValue(inputfile, delta_v0_hot, "delta_v0_hot"); 
  float alpha_v_k1; findValue(inputfile, alpha_v_k1, "alpha_v_k1"); 
  float alpha_v_k2; findValue(inputfile, alpha_v_k2, "alpha_v_k2"); 
  float delta_a_v_k; findValue(inputfile, delta_a_v_k, "delta_a_v_k"); 
  float RmaxSF1; findValue(inputfile, RmaxSF1, "RmaxSF1"); 
  float RmaxSF2; findValue(inputfile, RmaxSF2, "RmaxSF2"); 
  float delta_RmaxSF; findValue(inputfile, delta_RmaxSF, "delta_RmaxSF"); 
  float v_k_thres1; findValue(inputfile, v_k_thres1, "v_k_thres1"); 
  float v_k_thres2; findValue(inputfile, v_k_thres2, "v_k_thres2"); 
  float delta_v_k_t; findValue(inputfile, delta_v_k_t, "delta_v_k_t"); 
  float v_k_max1; findValue(inputfile, v_k_max1, "v_k_max1"); 
  float v_k_max2; findValue(inputfile, v_k_max2, "v_k_max2"); 
  float delta_v_k_max; findValue(inputfile, delta_v_k_max, "delta_v_k_max"); 
  float accr_norm1; findValue(inputfile, accr_norm1, "accr_norm1"); 
  float accr_norm2; findValue(inputfile, accr_norm2, "accr_norm2"); 
  float delta_accr; findValue(inputfile, delta_accr, "delta_accr"); 
  float accr_z1; findValue(inputfile, accr_z1, "accr_z1"); 
  float accr_z2; findValue(inputfile, accr_z2, "accr_z2"); 
  float d_accr_z; findValue(inputfile, d_accr_z, "d_accr_z");
  float ion_frac1; findValue(inputfile, ion_frac1, "ion_frac1");
  float ion_frac2; findValue(inputfile, ion_frac2, "ion_frac2");
  float delta_ion_frac; findValue(inputfile, delta_ion_frac, "delta_ion_frac");
  float t_drag1; findValue(inputfile, t_drag1, "t_drag1");
  float t_drag2; findValue(inputfile, t_drag2, "t_drag2");
  float delta_t_drag; findValue(inputfile, delta_t_drag, "delta_t_drag");
  int which_arm; findValue(inputfile,which_arm,"which_arm");
  float pitch;		
	if(which_arm==1) {findValue(inputfile,pitch,"a1_pitch");}
	else if(which_arm==2) {findValue(inputfile,pitch,"a2_pitch");}
	else if(which_arm==3) {findValue(inputfile,pitch,"a3_pitch");}
	else if(which_arm==4) {findValue(inputfile,pitch,"a4_pitch");}
	else {cout<<"ERROR: which_arm must be between 1 and 4"<<endl;}
  float delta_arm; findValue(inputfile, delta_arm, "delta_arm");
  float delta_disk; findValue(inputfile, delta_disk, "delta_disk");
  float RminSF; findValue(inputfile, RminSF, "RminSF");	
  float RmaxSF; findValue(inputfile, RmaxSF, "RmaxSF");		
  int n_part1; findValue(inputfile, n_part1, "n_part1");
  int n_part2; findValue(inputfile, n_part2, "n_part2");
  int delta_npart; findValue(inputfile, delta_npart, "delta_npart");

  strcpy(par1,param1);
  strcpy(par2,param2);
  strcpy(par3,param3);
	
  if (strcoll(par1,"M_HI_halo") ==0){
    param1_0=M_halo1;
    delta_param1=delta_M;
    n_iter1=(int)((M_halo2-M_halo1)/delta_M+0.5)+1;
  }
  else {
    if (strcoll(par1,"h_v_k") ==0){
      param1_0=h_v_k1;
      delta_param1=delta_h_v_k;
      n_iter1=(int)((h_v_k2-h_v_k1)/delta_h_v_k+0.5)+1;
    }
    else {
      if (strcoll(par1,"K_D") ==0){
	param1_0=K_D1;
	delta_param1=delta_K_D;
	n_iter1=(int)((K_D2-K_D1)/delta_K_D+0.5)+1;
      } 
      else {
	if (strcoll(par1,"v0_hot") ==0){
	  param1_0=v0_hot1;
	  delta_param1=delta_v0_hot;
	  n_iter1=(int)((v0_hot2-v0_hot1)/delta_v0_hot+0.5)+1;
	} 
	else {
	  if (strcoll(par1,"gammaSF") ==0){
	    param1_0=gammaSF1;
	    delta_param1=delta_gSF;
	    n_iter1=(int)((gammaSF2-gammaSF1)/delta_gSF+0.5)+1;
	  }
	  else {
	    if (strcoll(par1,"RmaxSF") ==0){
	      param1_0=RmaxSF1;
	      delta_param1=delta_RmaxSF;
	      n_iter1=(int)((RmaxSF2-RmaxSF1)/delta_RmaxSF+0.5)+1;
	    }
	    else {
	      if (strcoll(par1,"v_k_thres") ==0){
		param1_0=v_k_thres1;
		delta_param1=delta_v_k_t;
		n_iter1=(int)((v_k_thres2-v_k_thres1)/delta_v_k_t+0.5)+1;
	      } 
	      else {
		if (strcoll(par1,"v_k_max") ==0){
		  param1_0=v_k_max1;
		  delta_param1=delta_v_k_max;
		  n_iter1=(int)((v_k_max2-v_k_max1)/delta_v_k_max+0.5)+1;
		} 
		else {
		  if (strcoll(par1,"alpha_v_k") ==0){
		    param1_0=alpha_v_k1;
		    delta_param1=delta_a_v_k;
		    n_iter1=(int)((alpha_v_k2-alpha_v_k1)/delta_a_v_k+0.5)+1;
		  }
		  else {
		    if (strcoll(par1,"accr_norm") ==0){
		      param1_0=accr_norm1;
		      delta_param1=delta_accr;
		      n_iter1=(int)((accr_norm2-accr_norm1)/delta_accr+0.5)+1;
		    }
		    else {
		      if (strcoll(par1,"accr_z") ==0){
			    param1_0=accr_z1;
			    delta_param1=d_accr_z;
			    n_iter1=(int)((accr_z2-accr_z1)/d_accr_z+0.5)+1;
		      }
			  else {
				if (strcoll(par1,"ion_frac") ==0){
					param1_0=ion_frac1;
					delta_param1=delta_ion_frac;
					n_iter1=(int)((ion_frac2-ion_frac1)/delta_ion_frac+0.5)+1;
				}
				else {
					if (strcoll(par1,"n_part") ==0){
						param1_0=n_part1;
						delta_param1=delta_npart;
						n_iter1=(int)((n_part2-n_part1)/delta_npart+0.5)+1;
					}
					else {
						if (strcoll(par1,"showt_min") ==0){
							param1_0=showt_min1;
							delta_param1=delta_showt_min;
							n_iter1=(int)((showt_min2-showt_min1)/delta_showt_min+0.5)+1;
						}
						else {
							if (strcoll(par1,"showt_max") ==0){
								param1_0=showt_max1;
								delta_param1=delta_showt_max;
								n_iter1=(int)((showt_max2-showt_max1)/delta_showt_max+0.5)+1;
							}
							else {
								if (strcoll(par1,"t_drag") ==0){
									param1_0=t_drag1;
									delta_param1=delta_t_drag;
									n_iter1=(int)((t_drag2-t_drag1)/delta_t_drag+0.5)+1;
								}
								else {
									if (strcoll(par1,"showt2_min") ==0){
										param1_0=showt2_min1;
										delta_param1=delta_showt2_min;
										n_iter1=(int)((showt2_min2-showt2_min1)/delta_showt2_min+0.5)+1;
									}
									else {
										if (strcoll(par1,"kick_angle") ==0){
											param1_0=kick_angle1;
											delta_param1=delta_kick_angle;
											n_iter1=(int)((kick_angle2-kick_angle1)/delta_kick_angle+0.5)+1;
										}
										else {
											if(strcoll(par1,"arm_section")==0){
												   param1_0=RminSF;
												   delta_param1 = dtoR_logspiral(10+delta_arm,pitch)-dtoR_logspiral(10,pitch);
												   n_iter1=(int)((RmaxSF-RminSF)/delta_param1+0.5)+1;
										}
										else
										{
											if(strcoll(par1,"disk_section")==0){
												param1_0=RminSF;
												delta_param1 = delta_disk;
												n_iter1=(int)((RmaxSF-RminSF)/delta_param1+0.5)+1;
											}
										}
									}
								}
							}
						}
					}
				  }	
				}
			  }
			}
		  }
		}
	  }
	}
  }
}
    }
  }
}
  
  if (strcoll(par2,"M_HI_halo") ==0){
    param2_0=M_halo1;
    delta_param2=delta_M;
    n_iter2=(int)((M_halo2-M_halo1)/delta_M+0.5)+1;
  }
  else {
    if (strcoll(par2,"h_v_k") ==0){
      param2_0=h_v_k1;
      delta_param2=delta_h_v_k;
      n_iter2=(int)((h_v_k2-h_v_k1)/delta_h_v_k+0.5)+1;
    }
    else {
      if (strcoll(par2,"K_D") ==0){
	param2_0=K_D1;
	delta_param2=delta_K_D;
	n_iter2=(int)((K_D2-K_D1)/delta_K_D+0.5)+1;
      }
      else {
	if (strcoll(par2,"v0_hot") ==0){
	  param2_0=v0_hot1;
	  delta_param2=delta_v0_hot;
	  n_iter2=(int)((v0_hot2-v0_hot1)/delta_v0_hot+0.5)+1;
	} 
	else {
	  if (strcoll(par2,"gammaSF") ==0){
	    param2_0=gammaSF1;
	    delta_param2=delta_gSF;
	    n_iter2=(int)((gammaSF2-gammaSF1)/delta_gSF+0.5)+1;
	  }
	  else {
	    if (strcoll(par2,"RmaxSF") ==0){
	      param2_0=RmaxSF1;
	      delta_param2=delta_RmaxSF;
	      n_iter2=(int)((RmaxSF2-RmaxSF1)/delta_RmaxSF+0.5)+1;
	    } 
	    else {
	      if (strcoll(par2,"v_k_thres") ==0){
		param2_0=v_k_thres1;
		delta_param2=delta_v_k_t;
		n_iter2=(int)((v_k_thres2-v_k_thres1)/delta_v_k_t+0.5)+1;
	      } 
	      else {
		if (strcoll(par2,"v_k_max") ==0){
		  param2_0=v_k_max1;
		  delta_param2=delta_v_k_max;
		  n_iter2=(int)((v_k_max2-v_k_max1)/delta_v_k_max+0.5)+1;
		} 
		else {
		  if (strcoll(par2,"alpha_v_k") ==0){
		    param2_0=alpha_v_k1;
		    delta_param2=delta_a_v_k;
		    n_iter2=(int)((alpha_v_k2-alpha_v_k1)/delta_a_v_k+0.5)+1;
		  } 
		  else {
		    if (strcoll(par2,"accr_norm") ==0){
		      param2_0=accr_norm1;
		      delta_param2=delta_accr;
		      n_iter2=(int)((accr_norm2-accr_norm1)/delta_accr+0.5)+1;
		    }
		    else {
		      if (strcoll(par2,"accr_z") ==0){
				param2_0=accr_z1;
				delta_param2=d_accr_z;
				n_iter2=(int)((accr_z2-accr_z1)/d_accr_z+0.5)+1;
			//			cout << accr_z1 << " " << accr_z2 << " " << d_accr_z << endl;
			  }
			  else {
				  if (strcoll(par2,"ion_frac") ==0){
					  param2_0=ion_frac1;
					  delta_param2=delta_ion_frac;
					  n_iter2=(int)((ion_frac2-ion_frac1)/delta_ion_frac+0.5)+1;
				}
				else {
					if (strcoll(par2,"n_part") ==0){
						param2_0=n_part1;
						delta_param2=delta_npart;
						n_iter2=(int)((n_part2-n_part1)/delta_npart+0.5)+1;
					}
					else {
						if (strcoll(par2,"showt_min") ==0){
							param2_0=showt_min1;
							delta_param2=delta_showt_min;
							n_iter2=(int)((showt_min2-showt_min1)/delta_showt_min+0.5)+1;
						}
						else {
							if (strcoll(par2,"showt_max") ==0){
								param2_0=showt_max1;
								delta_param2=delta_showt_max;
								n_iter2=(int)((showt_max2-showt_max1)/delta_showt_max+0.5)+1;
							}
							else {
								if (strcoll(par2,"t_drag") ==0){
									param2_0=t_drag1;
									delta_param2=delta_t_drag;
									n_iter2=(int)((t_drag2-t_drag1)/delta_t_drag+0.5)+1;
								}
								else {
									if (strcoll(par2,"showt2_min") ==0){
										param2_0=showt2_min1;
										delta_param2=delta_showt2_min;
										n_iter2=(int)((showt2_min2-showt2_min1)/delta_showt2_min+0.5)+1;
									}
									else {
										if (strcoll(par2,"kick_angle") ==0){
											param2_0=kick_angle1;
											delta_param2=delta_kick_angle;
											n_iter2=(int)((kick_angle2-kick_angle1)/delta_kick_angle+0.5)+1;
										}
										else {
										if(strcoll(par2,"arm_section")==0)
										   {
											   param2_0=RminSF;
											   delta_param2 = dtoR_logspiral(10+delta_arm,pitch)-dtoR_logspiral(10,pitch);
											   n_iter2=(int)((RmaxSF-RminSF)/delta_param2+0.5)+1;
										   }
										else
										{
											if(strcoll(par2,"disk_section")==0){
												param2_0=RminSF;
												delta_param2 = delta_disk;
												n_iter2=(int)((RmaxSF-RminSF)/delta_param2+0.5)+1;
											}
										}
										   }
									}
								}
						}
					}
				}
			  }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
	}
  }
	if (strcoll(par3,"M_HI_halo") ==0){
		param3_0=M_halo1;
		delta_param3=delta_M;
		n_iter3=(int)((M_halo2-M_halo1)/delta_M+0.5)+1;
	}
	else {
		if (strcoll(par3,"h_v_k") ==0){
			param3_0=h_v_k1;
			delta_param3=delta_h_v_k;
			n_iter3=(int)((h_v_k2-h_v_k1)/delta_h_v_k+0.5)+1;
		}
		else {
			if (strcoll(par3,"K_D") ==0){
				param3_0=K_D1;
				delta_param3=delta_K_D;
				n_iter3=(int)((K_D2-K_D1)/delta_K_D+0.5)+1;
			} 
			else {
				if (strcoll(par3,"v0_hot") ==0){
					param3_0=v0_hot1;
					delta_param3=delta_v0_hot;
					n_iter3=(int)((v0_hot2-v0_hot1)/delta_v0_hot+0.5)+1;
				} 
				else {
					if (strcoll(par3,"gammaSF") ==0){
						param3_0=gammaSF1;
						delta_param3=delta_gSF;
						n_iter3=(int)((gammaSF2-gammaSF1)/delta_gSF+0.5)+1;
					}
					else {
						if (strcoll(par3,"RmaxSF") ==0){
							param3_0=RmaxSF1;
							delta_param3=delta_RmaxSF;
							n_iter3=(int)((RmaxSF2-RmaxSF1)/delta_RmaxSF+0.5)+1;
						}
						else {
							if (strcoll(par3,"v_k_thres") ==0){
								param3_0=v_k_thres1;
								delta_param3=delta_v_k_t;
								n_iter3=(int)((v_k_thres2-v_k_thres1)/delta_v_k_t+0.5)+1;
							} 
							else {
								if (strcoll(par3,"v_k_max") ==0){
									param3_0=v_k_max1;
									delta_param3=delta_v_k_max;
									n_iter3=(int)((v_k_max2-v_k_max1)/delta_v_k_max+0.5)+1;
								} 
								else {
									if (strcoll(par3,"alpha_v_k") ==0){
										param3_0=alpha_v_k1;
										delta_param3=delta_a_v_k;
										n_iter3=(int)((alpha_v_k2-alpha_v_k1)/delta_a_v_k+0.5)+1;
									}
									else {
										if (strcoll(par3,"accr_norm") ==0){
											param3_0=accr_norm1;
											delta_param3=delta_accr;
											n_iter3=(int)((accr_norm2-accr_norm1)/delta_accr+0.5)+1;
										}
										else {
											if (strcoll(par3,"accr_z") ==0){
												param3_0=accr_z1;
												delta_param3=d_accr_z;
												n_iter3=(int)((accr_z2-accr_z1)/d_accr_z+0.5)+1;
											}
											else {
												if (strcoll(par3,"ion_frac") ==0){
													param3_0=ion_frac1;
													delta_param3=delta_ion_frac;
													n_iter3=(int)((ion_frac2-ion_frac1)/delta_ion_frac+0.5)+1;
												}
												else {
													if (strcoll(par3,"n_part") ==0){
														param3_0=n_part1;
														delta_param3=delta_npart;
														n_iter3=(int)((n_part2-n_part1)/delta_npart+0.5)+1;
													}
													else {
														if (strcoll(par3,"showt_min") ==0){
															param3_0=showt_min1;
															delta_param3=delta_showt_min;
															n_iter3=(int)((showt_min2-showt_min1)/delta_showt_min+0.5)+1;
														}
														else {
															if (strcoll(par3,"showt_max") ==0){
																param3_0=showt_max1;
																delta_param3=delta_showt_max;
																n_iter3=(int)((showt_max2-showt_max1)/delta_showt_max+0.5)+1;
															}
															else {
																if (strcoll(par3,"t_drag") ==0){
																	param3_0=t_drag1;
																	delta_param3=delta_t_drag;
																	n_iter3=(int)((t_drag2-t_drag1)/delta_t_drag+0.5)+1;
																}
																else {
																	if (strcoll(par3,"showt2_min") ==0){
																		param3_0=showt2_min1;
																		delta_param3=delta_showt2_min;
																		n_iter3=(int)((showt2_min2-showt2_min1)/delta_showt2_min+0.5)+1;
																	}
																	else {
																		if (strcoll(par3,"kick_angle") ==0){
																			param3_0=kick_angle1;
																			delta_param3=delta_kick_angle;
																			n_iter3=(int)((kick_angle2-kick_angle1)/delta_kick_angle+0.5)+1;
																		}
																		else{
																			if(strcoll(par3,"arm_section")==0)
																			   {
																				   param3_0=RminSF;
																				   delta_param3 = dtoR_logspiral(15+delta_arm,pitch)-dtoR_logspiral(15,pitch);
																				   n_iter3=(int)((RmaxSF-RminSF)/delta_param3+0.5)+1;
																			   }
																			else
																			{
																				if(strcoll(par3,"disk_section")==0){
																					param3_0=RminSF;
																					delta_param3 = delta_disk;
																					n_iter3=(int)((RmaxSF-RminSF)/delta_param3+0.5)+1;
																				}
																			}
																			}
																	}
																}
														}
													}
												}	
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
 }
	
	
  par1_0=param1_0;
  par2_0=param2_0;
  par3_0=param3_0;  
  delta_par1=delta_param1;
  delta_par2=delta_param2;
  delta_par3=delta_param3;
  niter1=n_iter1;
  niter2=n_iter2;
  niter3=n_iter3;
  findValue(inputfile, res_above, "res_above");
  findValue(inputfile, res_Rsize, "res_Rsize");
  findValue(inputfile, res_zsize, "res_zsize");
  findValue(inputfile, freechannels, "freechannels");
  findValue(inputfile, target_chan1, "target_chan1");
  findValue(inputfile, target_chan2, "target_chan2");
}

void iterations::CalculateResiduals(int iter1, int iter2,
				    Array3D<float> rescube,
				    Array2D<float> rescube_tot,
				    char stop)
{
  static Array2D<float> counts_above_tot(n_iter1,n_iter2,0.f),
    counts_below_tot(n_iter1,n_iter2,0.f), 
    counts_above_left_tot(n_iter1,n_iter2,0.f), 
    counts_below_left_tot(n_iter1,n_iter2,0.f), 
    counts_average_tot(n_iter1,n_iter2,0.f);
  static Array2D<float> counts_above(n_iter1,n_iter2,0.f), 
    counts_below(n_iter1,n_iter2,0.f), 
    counts_above_app(n_iter1,n_iter2,0.f), 
    counts_below_app(n_iter1,n_iter2,0.f), 
    counts_average(n_iter1,n_iter2,0.f), 
    counts_targeted(n_iter1,n_iter2,0.f);
  
  int xc, yc, vc;
  above_pix=int(res_above*(kpctograd/fabs(cdelt2))+y_shift+0.5);
  below_pix=int(res_above*(kpctograd/fabs(cdelt1))-y_shift+0.5);

  if (res_Rsize == 0)
    Rsize_pix=RAsize/2-1;
  else
    Rsize_pix=int(res_Rsize*(kpctograd/fabs(cdelt1))+0.5)+1;
  if (res_zsize == 0)
    zsize_pix=DECsize/2-1;
  else 
    zsize_pix=int(res_zsize*(kpctograd/fabs(cdelt2))+0.5)+1;

  for (xc=RAsize/2-Rsize_pix;xc<RAsize/2+Rsize_pix;xc++){
    for (yc=DECsize/2-zsize_pix;yc<DECsize/2+zsize_pix;yc++){
      for (vc=0;vc<freechannels;vc++)
	absnoise+=fabs(rescube[xc][yc][vc]);
      for (vc=VELsize-freechannels;vc<VELsize;vc++)
	absnoise+=fabs(rescube[xc][yc][vc]);
    }
  }

  int n_above=0, n_below=0, n_above_app=0, n_below_app=0, n_targeted=0;
  // average absolute value of the data in mJy/Beam 
  absnoise=absnoise/((2*Rsize_pix)*(2*zsize_pix))/(2*freechannels); 
  noise=absnoise/sqrt(2./PI);
  if (iter1 == 0 && iter2==0) 
    cout << "Estimated noise in the data: " << noise*1000. << " mJy/Beam\n";

  for (vc=0;vc<VELsize;vc++)
    for (xc=RAsize/2-Rsize_pix;xc<RAsize/2+Rsize_pix;xc++)
      for (yc=DECsize/2+above_pix;yc<DECsize/2+zsize_pix;yc++)
	if (rescube[xc][yc][vc] > noise) {
	  counts_above[iter1][iter2]+=fabs(rescube[xc][yc][vc]);
	  n_above++;
	}
  for (vc=0;vc<VELsize;vc++)
    for (xc=RAsize/2-Rsize_pix;xc<RAsize/2+Rsize_pix;xc++)
      for (yc=DECsize/2-zsize_pix;yc<DECsize/2-below_pix;yc++)
	if (rescube[xc][yc][vc] > noise){
	  counts_below[iter1][iter2]+=fabs(rescube[xc][yc][vc]);
	  n_below++;
	}
  
  for (xc=RAsize/2-Rsize_pix;xc<RAsize/2+Rsize_pix;xc++)
    for (yc=DECsize/2+above_pix;yc<DECsize/2+zsize_pix;yc++)
      for (vc=0;vc<VELsize/2;vc++)
	if (rescube[xc][yc][vc] > noise) {
	  counts_above_app[iter1][iter2]+=fabs(rescube[xc][yc][vc]);
	  n_above_app++;
	}
  
  for (vc=0;vc<VELsize/2;vc++)
    for (xc=RAsize/2-Rsize_pix;xc<RAsize/2+Rsize_pix;xc++)
      for (yc=DECsize/2-zsize_pix;yc<DECsize/2-below_pix;yc++)
	if (rescube[xc][yc][vc] > noise) {
	  counts_below_app[iter1][iter2]+=fabs(rescube[xc][yc][vc]);
	  n_below_app++;
	}

  for (vc=target_chan1;vc<target_chan2;vc++)
    for (xc=RAsize/2-Rsize_pix;xc<RAsize/2;xc++)
      for (yc=DECsize/2+above_pix;yc<DECsize/2+zsize_pix;yc++)
	if (rescube[xc][yc][vc] > noise) {
	  counts_targeted[iter1][iter2]+=fabs(rescube[xc][yc][vc]);
	  n_targeted++;
	}

  for (vc=target_chan1;vc<target_chan2;vc++)
    for (xc=RAsize/2-Rsize_pix;xc<RAsize/2;xc++)
      for (yc=DECsize/2-zsize_pix;yc<DECsize/2-below_pix;yc++)
	if (rescube[xc][yc][vc] > noise) {
	  counts_targeted[iter1][iter2]+=fabs(rescube[xc][yc][vc]);
	  n_targeted++;
	}

  for (xc=RAsize/2-Rsize_pix;xc<RAsize/2+Rsize_pix;xc++)
    for (yc=DECsize/2+above_pix;yc<DECsize/2+zsize_pix;yc++)
      counts_above_tot[iter1][iter2]+=fabs(rescube_tot[xc][yc]);

  for (xc=RAsize/2-Rsize_pix;xc<RAsize/2+Rsize_pix;xc++)
    for (yc=DECsize/2-zsize_pix;yc<DECsize/2-below_pix;yc++)
      counts_below_tot[iter1][iter2]+=fabs(rescube_tot[xc][yc]);

  for (xc=RAsize/2-Rsize_pix;xc<RAsize/2+int(x_shift+.5);xc++)
    for (yc=DECsize/2+above_pix;yc<DECsize/2+zsize_pix;yc++)
      counts_above_left_tot[iter1][iter2]+=fabs(rescube_tot[xc][yc]);

  for (xc=RAsize/2-Rsize_pix;xc<RAsize/2+int(x_shift+.5);xc++)
    for (yc=DECsize/2-zsize_pix;yc<DECsize/2-below_pix;yc++)
      counts_below_left_tot[iter1][iter2]+=fabs(rescube_tot[xc][yc]);
  
  counts_above[iter1][iter2]=fabs(counts_above[iter1][iter2]/
    (2*Rsize_pix*(zsize_pix-above_pix))/VELsize-absnoise);

  counts_below[iter1][iter2]=fabs(counts_below[iter1][iter2]/
    (2*Rsize_pix*(zsize_pix-below_pix))/VELsize-absnoise);

  counts_above_app[iter1][iter2]=fabs(counts_above_app[iter1][iter2]/
    (Rsize_pix*(zsize_pix-above_pix))/(VELsize/2)-absnoise);

  counts_below_app[iter1][iter2]=fabs(counts_below_app[iter1][iter2]/
    (Rsize_pix*(zsize_pix-below_pix))/(VELsize/2)-absnoise);

  counts_above_tot[iter1][iter2]=counts_above_tot[iter1][iter2]/
    (2*Rsize_pix*(zsize_pix-above_pix));

  counts_below_tot[iter1][iter2]=counts_below_tot[iter1][iter2]/
    (2*Rsize_pix*(zsize_pix-below_pix));

  counts_above_left_tot[iter1][iter2]=counts_above_left_tot[iter1][iter2]/
    (Rsize_pix*(zsize_pix-above_pix));

  counts_below_left_tot[iter1][iter2]=counts_below_left_tot[iter1][iter2]/
    (Rsize_pix*(zsize_pix-below_pix));
  
  counts_average[iter1][iter2]=(counts_above[iter1][iter2]+
				counts_below[iter1][iter2]+
				counts_above_app[iter1][iter2]+
				counts_below_app[iter1][iter2])/4./noise;
  
  counts_targeted[iter1][iter2]=(counts_targeted[iter1][iter2]
				 /(Rsize_pix*(2*zsize_pix-above_pix-below_pix))
				 /fabs(target_chan2-target_chan1)
				 -absnoise)/noise;
  
  counts_average_tot[iter1][iter2]=(counts_above_tot[iter1][iter2]+
				    counts_below_tot[iter1][iter2]+
				    counts_above_left_tot[iter1][iter2]+
				    counts_below_left_tot[iter1][iter2])/4.;

  //  cout << noise << " "  << " " << n_iter1-1 << " " << n_iter2-1 << endl;
  if ((iter1 == (n_iter1-1) && iter2 == (n_iter2-1))
      || stop == 'n'){
    above_tot=counts_above_tot;
    below_tot=counts_below_tot;
    above_left_tot=counts_above_left_tot;
    below_left_tot=counts_below_left_tot;
    average_tot=counts_average_tot;
    above=counts_above;
    below=counts_below;
    above_app=counts_above_app;
    below_app=counts_below_app;
    average=counts_average;
    targeted=counts_targeted;
  }
}

void iterations::WriteResiduals()
{
  //  void writefits_2D(const char*, float*, const int, const int, const char*); 

  FILE* f_residuals=NULL;
  f_residuals=fopen("residuals.dat", "w");
  
  for (int i1=0; i1<n_iter1; i1++){
    for (int i2=0; i2<n_iter2; i2++){
      fprintf(f_residuals,"%.2f %.2f %f %f %f %f %f %f %f %f \n", 
	      par1_0+delta_par1*i1, par2_0+delta_par2*i2, 
	      above[i1][i2], below[i1][i2], 
	      above_app[i1][i2], below_app[i1][i2], 
	      above_tot[i1][i2], below_tot[i1][i2],
	      above_left_tot[i1][i2], below_left_tot[i1][i2]);
    }
  }
  fclose(f_residuals);
  plot9();
  printHeaderIter(par1, par2, par1_0, par2_0, delta_par1, delta_par2);
  writefits_2D("above.fits", &above[0][0], n_iter1, n_iter2, "header.tab");
  printHeaderIter(par1, par2, par1_0, par2_0, delta_par1, delta_par2);
  writefits_2D("below.fits", &below[0][0], n_iter1, n_iter2, "header.tab");
  printHeaderIter(par1, par2, par1_0, par2_0, delta_par1, delta_par2);
  writefits_2D("above_app.fits", &above_app[0][0], n_iter1, n_iter2, "header.tab");
  printHeaderIter(par1, par2, par1_0, par2_0, delta_par1, delta_par2);    
  writefits_2D("below_app.fits", &below_app[0][0], n_iter1, n_iter2, "header.tab");
  printHeaderIter(par1, par2, par1_0, par2_0, delta_par1, delta_par2);
  writefits_2D("average.fits", &average[0][0], n_iter1, n_iter2, "header.tab");
  printHeaderIter(par1, par2, par1_0, par2_0, delta_par1, delta_par2);
  writefits_2D("target.fits", &targeted[0][0], n_iter1, n_iter2, "header.tab");
  printHeaderIter(par1, par2, par1_0, par2_0, delta_par1, delta_par2);
  writefits_2D("average_tot.fits", &average_tot[0][0], n_iter1, n_iter2, "header.tab");
}
