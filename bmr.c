#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "random.h"
#include "field.h"
#include "bmr.h"

// return the sign of the hemisphere for now
static char sample_hemi() {
    return 2 * (ranq1() % 2) - 1;
}

static double sample_size_bmr() {
    return ranq1pwr(-p_bmr, smin, smax);
}

// Sample the colatitude of a BMR from gaussian with variable width.
// Schrijver 2001 ApJ 547 475 -- p487
static double sample_th(double flux) {
    double sig_th;
    sig_th = th0 * exp(-flux / flux_th) + th1;
    return randnorm(mu_th, sig_th);
}

static double sample_ph() {
    return ranq1dbl() * 2 * M_PI;
}

// Sample the inclination of a BMR from gaussian with variable width
// Schrijver 2001 ApJ 547 475 -- p487
static double sample_i(double flux) {
    double sig_i;
    sig_i = i0 * exp(-flux / flux_i) + i1;
    return randnorm(mu_i, sig_i);
}

// TODO: check order of ops for mod/%
static double cycle_activity(double time) {
    double t;
    t = time / tcycle;
    return A0 * beta * sin(2 * M_PI * t) * fmod(2 * t, 1) * exp(-5 * fmod(2 * t * t, 1));
}

static double activity = -A0;

void set_activity(double act) {
    activity = act;
}

/**
static double mean_field(double **grid, int ntheta, int nphi) {
    
    int i, j;
    double field;

    for ( i = 1 ; i < ntheta - 1 ; i++ ) {
        for ( j = 0 ; j < nphi ; j++ ) {
            field = fabs(grid[i][j]);
        }
    }

    field /= (ntheta - 2) * (nphi - 1);
    return field;
}

static double mean_flux(double **grid, int ntheta, int nphi) {
    
    int i, j;
    double dth, dph, th, dS, flux;

    dth = M_PI / (ntheta - 2);
    dph = 2 * M_PI / (ntheta - 1);

    for ( i = 1 ; i < ntheta - 1 ; i++ ) {
        th = (i - 0.5) * dth;
        dS = dth * dph * sin(th);
        for ( j = 0 ; j < nphi ; j++ ) {
            flux = fabs(grid[i][j] * dS);
        }
    }

    flux /= (ntheta - 2) * (nphi - 1);
    return flux;
}
**/

/**
 * Emergence of BMRs according to the prescription of Schrijver+ 01a.
 */
void schrijver(double **grid, int ntheta, int nphi, double time, double dt) {
    
    int i, j, k, N;
    double n, A, dth, dph, th, ph, u, v, size, flux, hemi, sep, colat, lng, inc, dist, sig, val;
    double *fluxes, *colats, *longs, *sigs;

    int t;


    double rangefactor = (pow(smax / field_rad / field_rad * 180 * 180 / M_PI / M_PI, 1-p_bmr) - pow(smin / field_rad / field_rad * 180 * 180 / M_PI / M_PI, 1-p_bmr)) / (1 - p_bmr);

    dth = M_PI / (ntheta - 2);
    dph = 2 * M_PI / (ntheta - 1);

//    avgflux = mean_flux(grid, ntheta, nphi);
//    avgfield = mean_field(grid, ntheta, nphi);

    // sample flux from power law
    A = ( activity < 0.0 ) ? -activity : cycle_activity(time);
    
    // number of spots
    n = a0 * A * dt * rangefactor;

    // fractional number of spots
    v = ranq1dbl();
    u = (v < (n - trunc(n))) ? 1 : 0;
    (void)t;
    /**
    if ( u > 0.0 ) {
        t = round(time/dt);

        printf("[%d] v %.3f n %.3f trunc(n) %d sunspot+1\n",t, v, n, (int) trunc(n));
    }
    **/

    N = (int) trunc(n + u);
    if ( N < 1.0 ) return;

    // separation btwn spots -- 18000km = 1.8e9cm
    sep = 1.8e9 / field_rad / 2;

    fluxes = (double *) malloc(2 * N * sizeof(double));
    colats = (double *) malloc(2 * N * sizeof(double));
    longs = (double *) malloc(2 * N * sizeof(double));
    sigs = (double *) malloc(2 * N * sizeof(double));

    for ( k = 0 ; k < N ; k++ ) {
        size = sample_size_bmr();
        flux = size * avgfluxd;
        hemi = sample_hemi();
        colat = sample_th(flux);
        lng = sample_ph();
        inc = sample_i(flux);

        // gaussian profile where FWHM is diameter of 'size'
        sig = sqrt(size / 2 / M_PI / M_LN2);

        // leading spot
        sigs[2*k] = sig;
        fluxes[2*k] = flux * hemi;
        colats[2*k] = colat - hemi * sep * sin(inc);
        longs[2*k] = lng + sep * cos(inc);
        
        // trailing spot
        sigs[2*k+1] = sig;
        fluxes[2*k+1] = -flux * hemi;
        colats[2*k+1] = colat + hemi * sep * sin(inc);
        longs[2*k+1] = lng - sep * cos(inc);

        printf("size = %.3e | sig = %.3e | flux = %.3e | colat = %.3e | lng = %.3e\n",
                size, sig, flux, colat, lng);
    }

    //exit(1);

    
    for ( i = 1 ; i < ntheta - 1 ; i++ ) {
        th = (i - 0.5) * dth;
        for ( j = 0 ; j < nphi ; j++ ) {
            ph = j * dph;
            for ( k = 0 ; k < 2 * N ; k++ ) {
                dist = sqrt(pow(th - colats[k], 2) + pow(ph - longs[k], 2)) / sigs[k];
                if ( dist >= 5 ) continue;
                val = fluxes[k] * exp(-pow(dist, 2) / 2);
//                printf("%.3e + %.3e\n", grid[i][j], val);
                grid[i][j] += val;
            }
        }
    }
}

/**
 * Emergence of BMRs according to the prescription of Lemerle+ 15.
 */
void lemerle(double **grid, int ntheta, int nphi, double time, double dt) {
    (void) grid;
    (void) ntheta;
    (void) nphi;
    (void) time;
    (void) dt;
}

/**
 * Naive emergence of a single bipolar magnetic region
 */
void naive(double **grid, int ntheta, int nphi, double time, double dt) {
    int i, j;
    double dist, distp, distn, val, th, ph, dth, dph;

    (void) time;
    (void) ph;

    dth = M_PI / (ntheta - 2);
    dph = 2 * M_PI / (nphi - 1);

    for ( i = 1 ; i < ntheta - 1 ; i++ ) {
        th = (i - 0.5) * dth;
        for ( j = 0 ; j < nphi ; j++ ) {
            ph = j * dph;
            dist = sqrt(pow(th - bmr_th, 2) + pow(ph - bmr_ph, 2)) / bmr_sigma;
            if ( dist >= 5 ) continue;
            distp = sqrt(pow(th - bmrp_th, 2) + pow(ph - bmrp_ph, 2));
            distn = sqrt(pow(th - bmrn_th, 2) + pow(ph - bmrn_ph, 2));
            val = exp(-pow(distp / bmr_sigma, 2) / 2) - exp(-pow(distn / bmr_sigma, 2) / 2);
            val *= bmr_b0;
            grid[i][j] += val;
        }
    }
}
