#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "bmr.h"
#include "constants.h"
#include "field.h"
#include "random.h"

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))

static double *fluxes;
static double *colats;
static double *longs;
static double *sigs;
static double rangefactor;
static int Nmax;


// TODO separate emergence into diff submodules

void init_schrijver(double activity) {

    double A;

    rangefactor = (pow(smax / field_rad / field_rad * 180 * 180 / M_PI / M_PI, 1-p_bmr) - pow(smin / field_rad / field_rad * 180 * 180 / M_PI / M_PI, 1-p_bmr)) / (1 - p_bmr);

    // 0.6 is larger than the max of the cycle function sin * mod * exp(mod)
    A = ( activity < 0.0 ) ? -activity : A0 * beta * 0.6;
    Nmax = (int) ceil(a0 * A * dt * rangefactor);

    // perhaps init this once with the max N could possibly be?
    fluxes = (double *) malloc(2 * Nmax * sizeof(double));
    colats = (double *) malloc(2 * Nmax * sizeof(double));
    longs = (double *) malloc(2 * Nmax * sizeof(double));
    sigs = (double *) malloc(2 * Nmax * sizeof(double));

}

// return the sign of the hemisphere for now
static char sample_hemi() {
    return 2 * (ranq1() % 2) - 1;
}

static double sample_size_bmr() {
    //static const double smax = fluxmax / avgfluxd; // in cm^2
    //static const double smin = fluxmin / avgfluxd; // in cm^2
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

static double cycle_activity(double time) {
    double t;
    t = time / tcycle;
    return A0 * beta * sin(2 * M_PI * t) * fmod(2 * t, 1) * exp(-5 * fmod(2 * t * t, 1));
}

void set_activity(double act) {
    activity = act;
}

/**
 * Emergence of BMRs according to the prescription of Schrijver+ 01a.
 * 'activity' param: curly A in Schrijver formulae (01a eq A1, 01b eq 1,2)
 */
void schrijver(double **grid, double time) {
    
    int i, ipole, imin, imax, j, jwrap, jmin, jmax, k, N;
    double n, A, th, ph, u, v, size, flux, hemi, colat, lng, inc, dist, distth, distph, sig, dS, omega;

    // sample flux from power law
    A = ( activity < 0.0 ) ? -activity : cycle_activity(time);
    
    // number of spots
    n = a0 * A * dt * rangefactor;

    // fractional number of spots
    v = ranq1dbl();
    u = (v < (n - trunc(n))) ? 1 : 0;
    N = (int) trunc(n + u); // TODO trunc redundant here?
    if ( N < 1 ) return;
    if ( N > Nmax ) {
        printf("Number of spots N exceeds calculated maximum: %d vs. Nmax %d\n",
                N, Nmax);
        return;
    }

    // TODO pre-initialize a bunch of spots in a grid and either sample or index
    // them with a total count. if not sampling but just making, create new once
    // you hit the size of the thing.

    for ( k = 0 ; k < N ; k++ ) {
        //size = sample_size_bmr();
        size = sample_size_bmr();
        omega = size / field_rad / field_rad;   // solid angle in rad^2
        flux = size * avgfluxd; // TODO check if flux should be divided by R^2
        hemi = sample_hemi();
        colat = sample_th(flux);
        if ( hemi < 0.0 ) colat = M_PI - colat;
        lng = sample_ph();
        inc = sample_i(flux);

        // gaussian profile where FWHM is diameter of 'size'
        sig = sqrt(omega / 2 / M_PI / M_LN2);

        // leading spot
        sigs[2*k] = sig;
        fluxes[2*k] = -flux * hemi;
        colats[2*k] = colat - hemi * spotsep * sin(inc);
        longs[2*k] = lng + spotsep * cos(inc);
        
        // trailing spot
        sigs[2*k+1] = sig;
        fluxes[2*k+1] = flux * hemi;
        colats[2*k+1] = colat + hemi * spotsep * sin(inc);
        longs[2*k+1] = lng - spotsep * cos(inc);

        //printf("size = %.3e | sig = %.3e | flux = %.3e | colat = %.3e | lng = %.3e\n",
        //        size, sig, flux, colat, lng);
    }

    // TODO this gets bad for N >= 20 spots... wow
    // use Bresenham's line algorithm / midpoint circle algorithm to draw each
    // level of the circle
    
    // TODO why does this impl produce diff results from prev? can't just be
    // the different order that floating pt arithmetic is performed...
    for ( k = 0 ; k < 2 * N ; k++ ) {

        // gaussian profile out to 5 sigma
        imin = (int) ((colats[k] - 5 * sigs[k] * sin(colats[k])) / dth);
        imax = (int) ((colats[k] + 5 * sigs[k] * sin(colats[k])) / dth);
        jmin = (int) ((longs[k] - 5 * sigs[k]) / dph);
        jmax = (int) ((longs[k] + 5 * sigs[k]) / dph);

        imin = MAX(1, imin);
        imax = MIN(ntheta - 1, imax);
        jmin = MAX(0, jmin);
        jmax = MIN(nphi, jmax);
        
        // TODO make polar spots wrap to ph + 180 if it goes over...
        for ( i = imin ; i < imax ; i++ ) {

            // north polar spot crosses pole
            /**
            if ( i > ntheta - 1 ) {
                ipole = ntheta - 1

            }
            **/
            ipole = i;

            th = (i - 0.5) * dth;
            dS = dth * dph * sin(th) * field_rad * field_rad;

            distth = th - colats[k];

            for ( j = jmin ; j < jmax ; j++ ) {

                // wrap around phi
                jwrap = j % jmax;

                ph = jwrap * dph;

                // TODO document this and change expr to trunc using bitwise ops
                distph = 2 * M_PI * ((ph - longs[k]) / 2 / M_PI - floor((ph - longs[k]) / 2 / M_PI + 1. / 2));

                // squared distance
                dist = (distth * distth + distph * distph) / sigs[k] / sigs[k];
                if ( dist <= 25 ) {
                    grid[ipole][jwrap] += fluxes[k] / dS * exp(-dist / 2);
                }
            }
        }
    }

    /**
    for ( i = 1 ; i < ntheta - 1 ; i++ ) {
        th = (i - 0.5) * dth;
        dS = dth * dph * sin(th) * field_rad * field_rad;
        for ( j = 0 ; j < nphi ; j++ ) {
            ph = j * dph;
            for ( k = 0 ; k < 2 * N ; k++ ) {
                distph = 2 * M_PI * ((ph - longs[k]) / 2 / M_PI - floor((ph - longs[k]) / 2 / M_PI + 1. / 2));
                dist = sqrt(pow(th - colats[k], 2) + pow(distph, 2)) / sigs[k];
                // phi wraps around like triangular wave form (dropped fabs b/c
                // squared)
                if ( dist >= 5 ) continue;
                val = fluxes[k] / dS * exp(-pow(dist, 2) / 2);
                grid[i][j] += val;
            }
        }
    }
    **/
}

/**
 * Emergence of BMRs according to the prescription of Lemerle+ 15.
 */
void lemerle(double **grid, double time) {
    (void) grid;
    (void) ntheta;
    (void) nphi;
    (void) time;
    (void) dt;
}

/**
 * Naive emergence of a single bipolar magnetic region
 * 'activity' param: interval of time btwn spot emergence.
 * sign of 'activity' flips inclination (anti-Hale).
 */
void naive(double **grid, double time) {
    int i, j, step;
    double dist, distph, val, th, ph, dth, dph, spot_th, spot_ph;

    step = (int) trunc(time / dt);
    double bmr_freq = fabs(activity);
    int hale = (0 < activity) - (activity < 0); // activity sign
    if ( step % ((int) trunc(bmr_freq / dt)) != 0 )
        return;

    dth = M_PI / (ntheta - 2);
    dph = 2 * M_PI / (nphi - 1);

    for ( i = 1 ; i < ntheta - 1 ; i++ ) {
        th = (i - 0.5) * dth;
        for ( j = 0 ; j < nphi ; j++ ) {
            ph = j * dph;
            val = 0.0;

            // leading spot northern hemisphere
            spot_th = bmr_th + spotsep * sin(bmr_i);
            spot_ph = bmr_ph + spotsep * cos(bmr_i);
            distph = 2 * M_PI * ((ph - spot_ph) / 2 / M_PI - floor((ph - spot_ph) / 2 / M_PI + 1. / 2));
            dist = sqrt(pow(th - spot_th, 2) + pow(distph, 2)) / bmr_sigma;
            if ( dist < 5 )
                val += hale * bmr_b0 * exp(-pow(dist, 2) / 2);

            //printf("bmr: %.2f | spot: %.2f | sep: %.2f | bmr_i %.2f\n", bmr_th, spot_th, spotsep * sin(bmr_i), bmr_i);

            // trailing spot northern hemisphere
            spot_th = bmr_th - spotsep * sin(bmr_i);
            spot_ph = bmr_ph - spotsep * cos(bmr_i);
            distph = 2 * M_PI * ((ph - spot_ph) / 2 / M_PI - floor((ph - spot_ph) / 2 / M_PI + 1. / 2));
            dist = sqrt(pow(th - spot_th, 2) + pow(distph, 2)) / bmr_sigma;
            if ( dist < 5 )
                val -= hale * bmr_b0 * exp(-pow(dist, 2) / 2);

            // leading spot southern hemisphere
            spot_th = M_PI - bmr_th - spotsep * sin(bmr_i);
            spot_ph = bmr_ph + spotsep * cos(bmr_i);
            distph = 2 * M_PI * ((ph - spot_ph) / 2 / M_PI - floor((ph - spot_ph) / 2 / M_PI + 1. / 2));
            dist = sqrt(pow(th - spot_th, 2) + pow(distph, 2)) / bmr_sigma;
            if ( dist < 5 )
                val -= hale * bmr_b0 * exp(-pow(dist, 2) / 2);

            // trailing spot southern hemisphere
            spot_th = M_PI - bmr_th + spotsep * sin(bmr_i);
            spot_ph = bmr_ph - spotsep * cos(bmr_i);
            distph = 2 * M_PI * ((ph - spot_ph) / 2 / M_PI - floor((ph - spot_ph) / 2 / M_PI + 1. / 2));
            dist = sqrt(pow(th - spot_th, 2) + pow(distph, 2)) / bmr_sigma;
            if ( dist < 5 )
                val += hale * bmr_b0 * exp(-pow(dist, 2) / 2);
            
            grid[i][j] += val;
        }
    }
}

void none(double **grid, double time) {}
