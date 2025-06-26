#include <math.h>
#include <stdlib.h>

#include <stdio.h>

#include "constants.h"
#include "field.h"

/**
 * Initialize 2D grid for magnetic field and apply initial condition.
 */
double **init_grid(int bprof) {
    
    int i, j;
    double th, ph;
    double **grid;
    initb_t init_field;

    // TODO refactor to be one contiguous location in memory
    grid = (double **) malloc(ntheta * sizeof(double *));
    init_field = init_fields[bprof];

    for ( i = 0 ; i < ntheta ; i++ ) {
        grid[i] = (double *) malloc(nphi * sizeof(double));
        th = dth * (i - 0.5);
        for ( j = 0 ; j < nphi ; j++ ) {
            ph = dph * j;
            grid[i][j] = init_field(th, ph);
        }
    }

    return grid;

}

/**
 * Create theta axis which ranges from [-dth/2, pi + dth/2]
 */
double *init_theta() {

    int i;
    double *theta;

    theta = (double *) malloc(ntheta * sizeof(double));

    for (i = 0 ; i < ntheta ; i++ ) {
        theta[i] = (i - 1/2) * dth;
    }

    return theta;

}

/**
 * Create phi axis which ranges from [0, 2pi]
 */
double *init_phi() {
    int i;
    double *phi;

    phi = (double *) malloc(nphi * sizeof(double));

    for (i = 0 ; i < nphi ; i++ ) {
        phi[i] = i * dph;
    }

    return phi;

}
