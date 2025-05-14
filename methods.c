#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "flow.h"
#include "field.h"

char verbose = 0;

static double *coefa;
static double *coefb;
static double coefc;
static double *coefd;
static double *coefe;

void init_coef(double *flow, double *grad, double *difr,
        double ntheta, double nphi, double dt) {

    coefa = (double *) malloc(ntheta * sizeof(double));
    coefb = (double *) malloc(ntheta * sizeof(double));
    coefd = (double *) malloc(ntheta * sizeof(double));
    coefe = (double *) malloc(ntheta * sizeof(double));

    int i;
    double val, th, dth, dph, cos_th, sin_th;

    dth = M_PI / (ntheta - 2);
    dph = 2 * M_PI / (nphi - 1);

    coefc = dt * field_eta / pow(field_rad * dth, 2);

    for ( i = 0 ; i < ntheta ; i++ ) {
        th = (i - 0.5) * dt;
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

double **ftcs(double **grid, double **newgrid, double *flow, double *grad, double *difr,
        int ntheta, int nphi, double dt) {

    int i, j;
    double val, term, th, dth, ph, dph;
    double bij, bipj, bimj, bijp, bijm;

    (void)ph;

    dth = M_PI / (ntheta - 2);
    dph = 2 * M_PI / (nphi - 1);

    for ( i = 1 ; i < ntheta - 1; i++ ) {
        th = (i - 0.5) * dth;
        for ( j = 0 ; j < nphi ; j++ ) {
            ph = j * dph;

            // store these values to avoid repeated computations
            bij = grid[i][j];
            bipj = grid[i+1][j];
            bimj = grid[i-1][j];
            bijp = grid[i][(j+1)%nphi];
            bijm = grid[i][(j-1+nphi)%nphi];

            // new value of field at (theta[i], phi[j])
            // term 0
            val = bij / dt;

            // term 1
            term = grad[i] * bij;
            term += flow[i] * (bipj - bimj) / 2 / dth;
            term += cos(th) / sin(th) * flow[i] * bij;
            term /= -field_rad;
            val += term;

            // term 2
            term = bijp - bijm;
            term *= -difr[i] / 2 / dph;
            val += term;

            // term 3
            term = (bipj + bimj - 2 * bij) / dth;
            term += cos(th) / sin(th) * (bipj - bimj) / 2;
            term *= field_eta / pow(field_rad,2) / dth;
            val += term;

            // term 4
            term = bijp + bijm - 2 * bij;
            term *= field_eta / pow(field_rad * sin(th) * dph,2);
            val += term;

            // term 5
            term = -bij / field_tau;
            val += term;

            // update new grid with new value
            newgrid[i][j] = val * dt;
        }
    }

    // Neumann boundary conditions - zero flux at poles
    for ( j = 0 ; j < nphi ; j++ ) {
        newgrid[0][j] = newgrid[1][j];
        newgrid[ntheta-1][j] = newgrid[ntheta-2][j];
    }

    return newgrid;
}

double **ftcs_tri(double **grid, dobule **newgrid, double *flow, double *grad, double *difr,
        int ntheta, int nphi, double dt) {

    int i, j;
    double val, term, th, dth, ph, dph;
    double bij, bipj, bimj, bijp, bijm;

    (void)ph;

    dth = M_PI / (ntheta - 2);
    dph = 2 * M_PI / (nphi - 1);

    th = (i - 0.5) * dth;
    ph = j * dph;

    return NULL;
}

