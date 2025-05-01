static const double flow_u0 = 12;   // maximum meridional flow speed m/s
static const double difr0 = 2.894e-6;   // equatorial angular velocity, rad/s
static const double difr_a2 = -0.1264;
static const double difr_a4 = -0.1591;

static const double flow_q = 7;
static const double flow_w = 8;
static const double flow_v = 2.0;
static const double flow_n = 1;

double *calc_flow(int ntheta);
double *calc_difr(int ntheta);
double *calc_flow_grad(int ntheta);
