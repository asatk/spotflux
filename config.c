/**
 * This file defines all of the constants for this 2D stellar surface flux
 * transport model.
 */

#define PI 3.14159265358979323846

////////////////////////////////////////////////////////////////////////////////

/** SIMULATION **/

const int ntheta = 128;
const int nphi = 256;
int nt = 20001;
double dt = 3e3;
const int freq = 1000;
const char *fname = "bfld.dat";

const double dth = PI / (ntheta - 2);
const double dph = 2 * PI / (nphi - 1);

// UPDATE
double alpha = 1.0;

// RANDOM
const unsigned long long seed = 0x2025;

////////////////////////////////////////////////////////////////////////////////

/** FIELD AND STAR **/

// Lemerle+ 15
const int bprof = 1;
const double field_b0 = 8.5;
const double field_eta = 3.5e12;
const double field_tau = 32 * 3.15e7;
const double field_rad = 6.957e10;

////////////////////////////////////////////////////////////////////////////////

/** FLUID FLOWS **/

const double flow_u0 = 1200;
const double difr0 = 2.894e-6;
const double difr_a2 = -0.1264;
const double difr_a4 = -0.1591;

// Lemerle+ 15 eqn 3 profile parameters
const double flow_q = 7;
const double flow_w = 8;
const double flow_v = 2.0;
const double flow_n = 1;

////////////////////////////////////////////////////////////////////////////////

/** BIPOLAR MAGNETIC REGION **/

const double theta_avg = (90 - 17.5) / 180 * PI;
const double inc_avg = 4.2 / 180 * PI;
const double spotsep = 1.8e9 / field_rad / 2;

//// SCHRIJVER+ 01 EMPIRICAL CYCLE
double activity = 3.15e7 / 52;
const double p_bmr = 1.9;
const double p_eph = 2.9;
const double A0 = 1.0;
const double beta = 1.0;
const double tcycle = 11 * 3.15e7;
const double a0 = 8.0 / 86400 * 180 * 180 / PI / PI;
const double a1 = 8.0 / 86400 * 180 * 180 / PI / PI;

//// SCHRIJVER+ 01 EMPIRICAL EMERGENCE
// FLUX
const double fluxmax = 1.5e22;
const double fluxmin = 1.2e19;
const double avgfluxd = 180;
// SIZE
const double smax = fluxmax / avgfluxd;
const double smin = fluxmin / avgfluxd;
// COLATITUDE
const double th0 = 25 / 180 * PI;
const double th1 = 4 / 180 * PI;
const double mu_th = theta_avg;
const double flux_th = 5e20;
// INCLINATION
const double i0 = PI / 2;
const double i1 = PI / 10;
const double mu_i = inc_avg;
const double flux_i = 8e19;

//// NAIVE EMERGENCE
const double bmr_b0 = 100.0;
const double bmr_sigma = 4. / 180 * PI;
const double bmr_th = theta_avg;
const double bmr_ph = PI;
const double bmr_i = inc_avg;

////////////////////////////////////////////////////////////////////////////////
