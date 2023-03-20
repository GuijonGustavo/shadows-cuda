//
///*Loop function*/
////void remited();
//
///*funtions that we will be use*/
//
// /*some variables*/
__device__ double r_a = 0.0;
// double epsilon(double,double); 
// double delta(double);
//
// /*4D metric*/
// double g_00(double,double);
// double g_RR(double,double);
// double g_0phi(double,double);
// double g_thetatheta(double,double);
// double g_phiphi(double,double);
//
// /*refraction index*/
// double n(double,double,double);
//
// /*derivates of the 4D metric*/
// double dg_00_R(double,double);
// double dg_00_theta(double,double);
// double dg_0phi_R(double,double);
// double dg_0phi_theta(double,double);
//
// /*3D metric*/
// double gamma_RR(double,double);
// double gamma_thetatheta(double,double);
// double gamma_phiphi(double,double);
//
// /*derivates of the 3D metric*/
// double dgamma_RR_R(double,double);
// double dgamma_RR_theta(double,double);
// double dgamma_thetatheta_R(double,double);
// double dgamma_thetatheta_theta(double,double);
// double dgamma_phiphi_R(double,double);
// double dgamma_phiphi_theta(double,double);
//
// /*evolution equations of e*/
// double common_term(double,double,double,double,double,double);
// double de_R(double,double,double,double,double,double);
// double de_theta(double,double,double,double,double,double);
// double de_phi(double,double,double,double,double,double);

//int main ()
//{
//  loop();
//  return 0;
//}
//
__device__
double epsilon(double R,double theta)
  {
    return R*R + r_a*r_a*cosf(theta)*cosf(theta);
  }

__device__
double delta(double R)
  {
    return R*R - R + r_a*r_a;
  }

__device__
double g_00(double R, double theta)
  {
    return 1.0 - R/epsilon(R,theta);
  }

__device__
double g_RR(double R,double theta)
  {
    return -1.0*epsilon(R,theta)/delta(R);
  }

__device__
double g_0phi(double R,double theta)
  {
    return 2.0*r_a*R*sinf(theta)*sinf(theta)/epsilon(R,theta);
  }

__device__
double g_thetatheta(double R,double theta)
  {
    return -1.0*epsilon(R,theta);
  }

__device__
double g_phiphi(double R,double theta)
  {
    return -1.0*( R*R + r_a*r_a + r_a*r_a*R*sinf(theta)*sinf(theta)/epsilon(R,theta) )*sinf(theta)*sinf(theta);
  }

__device__
double n(double R,double theta,double e_phi)
  {
    return epsilon(R,theta)*( sqrt(1.0-R/epsilon(R,theta)) - 2.0*r_a*R*sinf(theta)*sinf(theta)*e_phi/epsilon(R,theta) )/(epsilon(R,theta)-R);
  }

__device__
double dg_00_R(double R,double theta)
  {
    return 2.0*R*R/(epsilon(R,theta)*epsilon(R,theta)) - 1.0/epsilon(R,theta);
  }

__device__
double dg_00_theta(double R,double theta)
  {
    return -2.0*r_a*r_a*R*cosf(theta)*sinf(theta)/(epsilon(R,theta)*epsilon(R,theta));
  }

__device__
double dg_0phi_R(double R,double theta)
  {
    return 2.0*r_a*sinf(theta)*sinf(theta)*( r_a*r_a*cosf(theta)*cosf(theta) - R*R )/(epsilon(R,theta)*epsilon(R,theta));
  }

__device__
double dg_0phi_theta(double R,double theta)
  {
    return 4.0*r_a*R*sinf(theta)*cosf(theta)*( R*R + r_a*r_a )/(epsilon(R,theta)*epsilon(R,theta)); 
  }

__device__
double gamma_RR(double R,double theta)
  {
    return epsilon(R,theta) / delta(R);
  }

__device__
double gamma_thetatheta(double R,double theta)
  {
    return epsilon(R,theta);
  }

__device__
double gamma_phiphi(double R,double theta)
  {
    return ( R*R + r_a*r_a + r_a*r_a*R*sinf(theta)*sinf(theta)/epsilon(R,theta) + 2.0*r_a*R/(epsilon(R,theta)-R)) * sinf(theta)*sinf(theta);
  }

__device__
double dgamma_RR_R(double R,double theta)
  {
    return (2.0*R*r_a*r_a*sinf(theta)*sinf(theta) + r_a*r_a*cosf(theta)*cosf(theta) - R*R) / (delta(R)*delta(R));
  }

__device__
double dgamma_RR_theta(double R,double theta)
  {
    return -2.0*r_a*r_a*sinf(theta)*cosf(theta)/delta(R);
  }

__device__
double dgamma_thetatheta_R(double R,double theta)
  {
    return 2.0*R;
  }

__device__
double dgamma_thetatheta_theta(double R,double theta)
  {
    return -2.0*r_a*r_a*sinf(theta)*cosf(theta);
  }

__device__
double dgamma_phiphi_R(double R,double theta)
  {
    return ( 2.0*R + r_a*sinf(theta)*sinf(theta)*(r_a*r_a*cosf(theta)*cosf(theta)-R*R)*(r_a/(epsilon(R,theta)*epsilon(R,theta)) + 2.0/((epsilon(R,theta)-R)*(epsilon(R,theta)-R))) ) * sinf(theta)*sinf(theta);
  }

__device__
double dgamma_phiphi_theta(double R,double theta)
  {
    return  sinf(theta)*sinf(theta)*( 2.0*r_a*r_a*R*sinf(theta)*cosf(theta)*(R*R+r_a*r_a) / (epsilon(R,theta)*epsilon(R,theta)) \
		    		    + 4.0*r_a*r_a*r_a*R*sinf(theta)*cosf(theta)/( (epsilon(R,theta)-R)*(epsilon(R,theta)-R) ) )\
	    + 2.0*sinf(theta)*cosf(theta)*( R*R + r_a*r_a + r_a*r_a*R*sinf(theta)*sinf(theta)/epsilon(R,theta) + 2.0*r_a*R/(epsilon(R,theta)-R));
  }

__device__
double common_term(double R,double theta,double phi,double e_r,double e_theta,double e_phi)
  {
    return ( e_r*( (0.5*dg_00_R(R,theta)/sqrt(g_00(R,theta))-e_phi*dg_0phi_R(R,theta))/g_00(R,theta) \
	            - dg_00_R(R,theta)*(sqrt(g_00(R,theta))-g_0phi(R,theta)*e_phi)/(g_00(R,theta)*g_00(R,theta)) ) \
	     + e_theta*( (0.5*dg_00_theta(R,theta)/sqrt(g_00(R,theta))-e_phi*dg_0phi_theta(R,theta))/g_00(R,theta) \
	                  - dg_00_theta(R,theta)*(sqrt(g_00(R,theta))-g_0phi(R,theta)*e_phi)/(g_00(R,theta)*g_00(R,theta)) ) ) \
	   / n(R,theta,e_phi);
  }

__device__
double de_r(double R,double theta,double phi,double e_r,double e_theta,double e_phi)
  {
    return ( ( (0.5*dg_00_R(R,theta)/sqrt(g_00(R,theta))-e_phi*dg_0phi_R(R,theta))/g_00(R,theta) \
	     - dg_00_R(R,theta)*(sqrt(g_00(R,theta))-g_0phi(R,theta)*e_phi)/(g_00(R,theta)*g_00(R,theta)) ) ) \
	   /(n(R,theta,e_phi)*gamma_RR(R,theta)) \
	   - e_r*common_term(R,theta,phi,e_r,e_theta,e_phi) \
	   + ( 0.5*( dgamma_RR_R(R,theta)*e_r*e_r + dgamma_thetatheta_R(R,theta)*e_theta*e_theta \
	            + dgamma_phiphi_R(R,theta)*e_phi*e_phi ) \
	      - (dgamma_RR_R(R,theta)*e_r*e_r + dgamma_RR_theta(R,theta)*e_theta*e_r) ) / gamma_RR(R,theta);
  }

__device__
double de_theta(double R,double theta,double phi,double e_r,double e_theta,double e_phi)
  {
    return ( ( (0.5*dg_00_theta(R,theta)/sqrt(g_00(R,theta))-e_phi*dg_0phi_theta(R,theta))/g_00(R,theta) \
	     - dg_00_theta(R,theta)*(sqrt(g_00(R,theta))-g_0phi(R,theta)*e_phi)/(g_00(R,theta)*g_00(R,theta)) ) ) \
	   /(n(R,theta,e_phi)*gamma_thetatheta(R,theta)) \
	   - e_theta*common_term(R,theta,phi,e_r,e_theta,e_phi) \
	   + ( 0.5*( dgamma_RR_theta(R,theta)*e_r*e_r + dgamma_thetatheta_theta(R,theta)*e_theta*e_theta \
	            + dgamma_phiphi_theta(R,theta)*e_phi*e_phi ) \
              - (dgamma_thetatheta_R(R,theta)*e_r*e_theta + dgamma_thetatheta_theta(R,theta)*e_theta*e_theta) ) / gamma_thetatheta(R,theta);
  }

__device__
double de_phi(double R,double theta,double phi,double e_r,double e_theta,double e_phi)
  {
    return -1.0*(dgamma_phiphi_R(R,theta)*e_r*e_phi + dgamma_phiphi_theta(R,theta)*e_theta*e_phi)/gamma_phiphi(R,theta) \
	   - e_phi*common_term(R,theta,phi,e_r,e_theta,e_phi); 
  }


