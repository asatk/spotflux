#include <stdio.h>
#include <stdlib.h>

#include "constants.h"
#include "linalg.h"

static double *gam_th;
static double *gam_ph;
static double *bb;
static double *u;
static double *z;

/**
 * Initialize the solvers of tridiagonal systems.
 */
void init_solvers() {

    gam_th = (double *) malloc(ntheta * sizeof(double));
    gam_ph = (double *) malloc(nphi * sizeof(double));

    bb = (double *) malloc(nphi * sizeof(double));
    u = (double *) malloc(nphi * sizeof(double));
    z = (double *) malloc(nphi * sizeof(double));
}


/**
 * Solves a tridiagonal system of linear equations A . x = r
 * Source: NR (Press 2007)
 */
void tridiag(double *a, double *b, double *c, double *r, double *x, int n) {
    int i;
    double beta, *gam;

    if ( b[0] == 0.0 ) {
        printf("tridiag -- first element on primary diagonal b[0] must be nonzero.\n");
        exit(1);
    }

    gam = (n == nphi) ? gam_ph : gam_th;
    
    beta = b[0];
    x[0] = r[0] / beta;

    // decomposition and forward substitution
    for ( i = 1 ; i < n ; i++ ) {
        gam[i] = c[i-1] / beta;
        beta = b[i] - a[i] * gam[i];
        if (beta == 0.0) {
            printf("tridiag -- beta is zero (zero pivot)\n");
            exit(1);
        }
        x[i] = (r[i] - a[i] * x[i-1]) / beta;
    }

    // back-substitution
    for ( i = n - 2 ; i >= 0 ; i-- )
        x[i] -= gam[i+1] * x[i+1];
}

/**
 * Solves a cyclic tridiagonal system of linear equations A . x = r where A is a
 * tridiagonal matrix with the lower, primary, and upper diagonal values given
 * by the vectors a, b, and c. The corner elements are given by alpha in the
 * first row and beta in the first column.
 * Assumes bottom left element ((N-1,0) aka alpha) is cN-1
 * Assumes top right element ((0,N-1) aka beta) is a0
 */
void tridiag_c(double *a, double *b, double *c, double *r, double *x, int n) {
    int i;
    double factor, alpha, beta, gamma;

    alpha = c[n-1];
    beta = a[0];

    if ( n <= 2 ) {
        printf("tridiag_c -- n too small for cyclic tridiagonal solver");
        exit(1);
    }

    gamma = -b[0];
    bb[0] = b[0] - gamma;
    bb[n-1] = b[n-1] - alpha * beta / gamma;

    for ( i = 1 ; i < n - 1 ; i++ )
        bb[i] = b[i];

    // solve A . x = r
    tridiag(a, bb, c, r, x, n);

    u[0] = gamma;
    u[n-1] = alpha;

    for ( i = 1 ; i < n - 1 ; i++ )
        u[i] = 0.0;

    // solve A . z = u
    tridiag(a, bb, c, u, z, n);

    // factor from Sherman-Morrison formula v . z / (1 + v . z)
    factor = (x[0] + beta * x[n-1] / gamma) / 
        (1.0 + z[0] + beta * z[n-1] / gamma);
    for ( i = 0 ; i < n ; i++ )
        x[i] -= factor * z[i];
}
