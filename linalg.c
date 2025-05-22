#include <stdio.h>
#include <stdlib.h>

#include "linalg.h"

//void init_tridiag_c(int *n);

/**
 * Solves a tridiagonal system of linear equations A . x = r
 * Source: NR (Press 2007)
 */
void tridiag(double *a, double *b, double *c, double *r, double *x, int n) {
    int i;
    double beta, *gamma;

    if ( b[0] == 0.0 ) {
        printf("tridiag -- first element on primary diagonal b[0] must be nonzero.\n");
        exit(1);
    }
    
    beta = b[0];
    x[0] = r[0] / beta;

    // decomposition and forward substitution
    for ( i = 1 ; i < n ; i++ ) {
        gamma[i] = c[i-1] / beta;
        beta = b[i] - a[i] * gamma[i];
        if (beta == 0.0) {
            printf("tridiag -- beta is zero (zero pivot)\n");
            exit(1);
        }
        x[i] = (r[i] - a[i] * x[i-1]) / beta;
    }

    // back-substitution
    for ( i = n - 2 ; i >= 0 ; i-- )
        x[i] -= gamma[i+1] * x[i+1];
}

/**
 * Solves a cyclic tridiagonal system of linear equations A . x = r where A is a
 * tridiagonal matrix with the lower, primary, and upper diagonal values given
 * by the vectors a, b, and c. The corner elements are given by alpha in the
 * first row and beta in the first column.
 */
void tridiag_c(double *a, double *b, double *c, double *r, double *x, int n,
        const double beta, const double alpha) {
    int i;
    double factor, gamma, *bb, *u, *z;

    if ( n <= 2 ) {
        printf("tridiag_c -- n too small for cyclic tridiagonal solver");
        exit(1);
    }

    bb = (double *) malloc(n * sizeof(double));
    u = (double *) malloc(n * sizeof(double));
    z = (double *) malloc(n * sizeof(double));

    gamma = -b[0];
    bb[0] = b[0] - gamma;
    bb[n-1] = b[n-1] - alpha * beta / gamma;

    for ( i = 1 ; i < n - 1 ; i++ )
        bb[i] = b[i];

    // solve A . x = r
    tridiag(a, bb, c, r, x);

    u[0] = gamma;
    u[n-1] = alpha;

    for ( i = 1 ; i < n - 1 ; i++ )
        u[i] = 0.0;

    // solve A . z = u
    tridiag(a, bb, c, u, z);

    // factor from Sherman-Morrison formula v . z / (1 + v . z)
    factor = (x[0] + beta * x[n-1] / gamma) / 
        (1.0 + z[0] + beta * z[n-1] / gamma);
    for ( i = 0 ; i < n ; i++ )
        x[i] -= factor * z[i];
}
