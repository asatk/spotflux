/**
 * This file defines all of the constants for this 2D stellar surface flux
 * transport model.
 */

#include "math.h"

////////////////////////////////////////////////////////////////////////////////

/** SIMULATION **/

const int ntheta = 128;
const int nphi = 256;
const int nt = 2001;
const double dt = 3e3;
const int freq = 20;
const char *fname = "bfld.dat";

// UPDATE
double alpha = 1.0;

// RANDOM
const unsigned long long seed = 0x2025LL;

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

#define THETA_AVG (90 - 17.5) / 180 * M_PI
#define INC_AVG 4.2 / 180 * M_PI

//// SCHRIJVER+ 01 EMPIRICAL CYCLE
double activity = 3.15e7 / 52;
const double p_bmr = 1.9;
const double p_eph = 2.9;
const double A0 = 1.0;
const double beta = 1.0;
const double tcycle = 11 * 3.15e7;
const double a0 = 8.0 / 86400 * 180 * 180 / M_PI / M_PI;
const double a1 = 8.0 / 86400 * 180 * 180 / M_PI / M_PI;

//// SCHRIJVER+ 01 EMPIRICAL EMERGENCE
// FLUX
const double fluxmax = 1.5e22;
const double fluxmin = 1.2e19;
const double avgfluxd = 180;
// SIZE
const double smax = fluxmax / avgfluxd;
const double smin = fluxmin / avgfluxd;
// COLATITUDE
const double th0 = 25 / 180 * M_PI;
const double th1 = 4 / 180 * M_PI;
const double mu_th = THETA_AVG;
const double flux_th = 5e20;
// INCLINATION
const double i0 = M_PI / 2;
const double i1 = M_PI / 10;
const double mu_i = INC_AVG;
const double flux_i = 8e19;

//// NAIVE EMERGENCE
const double bmr_b0 = 100.0;
const double bmr_sigma = 4. / 180 * M_PI;
const double bmr_th = THETA_AVG;
const double bmr_ph = M_PI;
const double bmr_i = INC_AVG;
const double bmr_sep = 1.8e9 / 2;

////////////////////////////////////////////////////////////////////////////////
