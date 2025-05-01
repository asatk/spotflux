#include <math.h>
#include <stdlib.h>

#include "field.h"

/**
 * Initialize 2D grid for magnetic field and apply initial condition.
 */
double **init_grid(int ntheta, int nphi, char bprof) {
    
    int i, j;
    double th, ph;
    double **grid;
    initb_t init_field;

    grid = (double **) malloc(ntheta * nphi * sizeof(double));
    init_field = init_fields[(int) bprof];

    for ( i = 0 ; i < ntheta ; i++ ) {
        th = M_PI * (2 * i - 1) / 2 / (ntheta - 2);
        for ( j = 0 ; j < nphi ; j++ ) {
            ph = 2 * M_PI * i / (nphi - 1);
            grid[i][j] = init_field(th, ph);
        }
    }

    return grid;

}

/**
 * Create theta axis which ranges from [-dth/2, pi + dth/2]
 */
double *init_theta(int ntheta) {

    int i;
    double dth, *theta;

    theta = (double *) malloc(ntheta * sizeof(double));
    dth = M_PI / (ntheta - 2);

    for (i = 0 ; i < ntheta ; i++ ) {
        theta[i] = (i - 1/2) * dth;
    }

    return theta;

}

/**
 * Create phi axis which ranges from [0, 2pi]
 */
double *init_phi(int nphi) {
    int i;
    double dph, *phi;

    phi = (double *) malloc(nphi * sizeof(double));
    dph = 2 * M_PI / (nphi - 1);

    for (i = 0 ; i < nphi ; i++ ) {
        phi[i] = i * dph;
    }

    return phi;

}
