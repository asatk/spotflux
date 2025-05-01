flow_u0 = 12;   // maximum meridional flow speed m/s
const double diffr0 = 2.894e-6;   // equatorial angular velocity, rad/s
const double diffr_a2 = -0.1264;
const double diffr_a4 = -0.1591;

double flow_q = 7;
double flow_w = 8;
double flow_v = 2.0;
double flow_n = 1;

double *calc_flow(int ntheta);
double *calc_difr(int ntheta);
double *calc_flow_grad(int ntheta);
