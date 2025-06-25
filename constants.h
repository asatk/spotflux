/**
 * This file defines all of the constants for this 2D stellar surface flux
 * transport model.
 */

#ifndef CONSTANTS_H
#define CONSTANTS_H

////////////////////////////////////////////////////////////////////////////////

/** BIPOLAR MAGNETIC REGION **/

//#define THETA_AVG (90 - 17.5) / 180 * M_PI
//#define INC_AVG 4.2 / 180 * M_PI

//// SCHRIJVER+ 01 EMPIRICAL CYCLE
extern const double p_bmr;          // power law exponent for active regions
extern const double p_eph;          // power law exponent for ephemeral regions
extern const double A0;             // activity factor (rel. to Sun at max)
extern const double beta;           // normalization of activity variation
extern const double tcycle;         // solar cycle 11yrs in s
extern const double a0;             // power law coef (rad^-2 s^-1 per hemi)
extern const double a1;             // power law coef (s)

//// SCHRIJVER+ 01 EMPIRICAL EMERGENCE
// FLUX
extern const double fluxmax;        // max BMR flux (Mx)
extern const double fluxmin;        // min BMR flux (Mx)
extern const double avgfluxd;       // avg flux density (field) in G
// SIZE
extern const double smax;           // in cm^2
extern const double smin;           // in cm^2
extern double rangefactor;
// COLATITUDE
extern const double th0;
extern const double th1;
extern const double mu_th;
extern const double flux_th;        // flux normalization for colat in Mx
// INCLINATION
extern const double i0;
extern const double i1;
extern const double mu_i;
extern const double flux_i;         // flux normalization for inc in Mx

//// NAIVE EMERGENCE
extern const double bmr_b0;         // strength of BMR (G)
extern const double bmr_sigma;      // spread of BMR
extern const double bmr_th;
extern const double bmr_ph;
extern const double bmr_i;
extern const double bmr_sep;

////////////////////////////////////////////////////////////////////////////////

/** FIELD AND STAR **/

// Lemerle+ 15
extern const double field_b0;       // field strength (G)
extern const double field_eta;      // diffusivity (cm^2/s)
extern const double field_tau;      // sink time (s)
extern const double field_rad;      // stellar radius (cm)

////////////////////////////////////////////////////////////////////////////////

/** FLUID FLOWS **/

extern const double flow_u0;        // maximum meridional flow speed cm/s
extern const double difr0;          // equatorial angular velocity, rad/s
extern const double difr_a2;        // difr quadratic term coef
extern const double difr_a4;        // difr quartic term coef

// Lemerle+ 15 eqn 3 profile parameters
extern const double flow_q;
extern const double flow_w;
extern const double flow_v;
extern const double flow_n;

////////////////////////////////////////////////////////////////////////////////
#endif  // CONSTANTS_H
