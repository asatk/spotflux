#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "flow.h"
#include "field.h"
#include "methods.h"

char verbose = 0;

static double *coefa;
static double *coefb;
static double coefc;
static double *coefd;
static double *coefe;

// TODO make all theta vecs ntheta but set end elmts to 0
void init_coef(double *flow, double *grad, double *difr,
        double ntheta, double nphi, double dt) {

    coefa = (double *) malloc((ntheta - 1) * sizeof(double));
    coefb = (double *) malloc((ntheta - 1) * sizeof(double));
    coefd = (double *) malloc((ntheta - 1) * sizeof(double));
    coefe = (double *) malloc((ntheta - 1) * sizeof(double));

    int i;
    double val, th, dth, dph, cos_th, sin_th;

    dth = M_PI / (ntheta - 2);
    dph = 2 * M_PI / (nphi - 1);

    coefa[0] = 0;
    coefb[0] = 0;
    coefc = dt * field_eta / pow(field_rad * dth, 2);
    coefd[0] = 0;
    coefe[0] = 0;

    for ( i = 1 ; i < ntheta - 1 ; i++ ) {
        th = (i - 0.5) * dth;
        cos_th = cos(th);
        sin_th = sin(th);

        // linear term
        val = (grad[i] + flow[i] * cos_th / sin_th) / field_rad;
        val = dt / 2 * (val - 1 / field_tau);
        coefa[i] = val;

        // first deriv theta
        val = flow[i] + field_eta * cos_th / sin_th / field_rad;
        val *= dt / 2 / field_rad / dth;
        coefb[i] = val;

        // first deriv phi
        coefd[i] = - dt * difr[i] / 2 / dph;

        // second deriv phi
        coefe[i] = dt * field_eta / pow(field_rad * sin_th * dph, 2);

    }
}

void ftcs(double **grid, double **newgrid, double *flow, double *grad, double *difr,
        int ntheta, int nphi, double dt) {

    int i, j;
    double val;

    for ( i = 1 ; i < ntheta - 1; i++ ) {
        for ( j = 0 ; j < nphi ; j++ ) {

            val = (1 + (2 * (coefa[i] - coefc - coefe[i]))) * grid[i][j];
            val += (coefc + coefb[i]) * grid[i+1][j];
            val += (coefc - coefb[i]) * grid[i-1][j];
            val += (coefe[i] + coefd[i]) * grid[i][(j+1)%nphi];
            val += (coefe[i] - coefd[i]) * grid[i][(j-1+nphi)%nphi];
            newgrid[i][j] = val;
        }
    }

    // Neumann boundary conditions - zero flux at poles
    for ( j = 0 ; j < nphi ; j++ ) {
        newgrid[0][j] = newgrid[1][j];
        newgrid[ntheta-1][j] = newgrid[ntheta-2][j];
    }
}

/**
void ftcs_tri(double **grid, double **newgrid, double *flow, double *grad, double *difr,
        int ntheta, int nphi, double dt) {

    int i, j;
    double val, th, dth, dph;


    dth = M_PI / (ntheta - 2);
    dph = 2 * M_PI / (nphi - 1);

    th = (i - 0.5) * dth;
    ph = j * dph;

}
 */
