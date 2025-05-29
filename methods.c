#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "linalg.h"
#include "flow.h"
#include "field.h"
#include "methods.h"

char verbose = 0;

// Intermediate grid -- index: phi (j), theta (i)
static double **grids1;

// STAGE 1 tridiag system
static double *msub;
static double *mpri;
static double *msup;
static double *q1;
static double *q2;
static double *q3;
static double *q4;
static double *q5;
static double *rt;

// STAGE 2 tridiag system
static double **asub;
static double **apri;
static double **asup;
static double *b1;
static double *b2;
static double *b3;
static double *b4;
static double *b5;
static double *rp;

void init_ftcs(double *flow, double *grad, double *difr,
        double ntheta, double nphi, double dt, double alpha) {

    int i, j;
    double coefa, coefb, coefc, coefd, coefe, th, dth, dph, cos_th, sin_th;

    // intermediate grid going from stage 1 to stage 2
    grids1 = (double **) malloc(nphi * sizeof(double));

    dth = M_PI / (ntheta - 2);
    dph = 2 * M_PI / (nphi - 1);

    // second deriv theta
    coefc = dt * field_eta / pow(field_rad * dth, 2);

    // STAGE 1 INITIALIZATION
    msub = (double *) malloc(ntheta * sizeof(double));
    mpri = (double *) malloc(ntheta * sizeof(double));
    msup = (double *) malloc(ntheta * sizeof(double));
    q1 = (double *) malloc(ntheta * sizeof(double));
    q2 = (double *) malloc(ntheta * sizeof(double));
    q3 = (double *) malloc(ntheta * sizeof(double));
    q4 = (double *) malloc(ntheta * sizeof(double));
    q5 = (double *) malloc(ntheta * sizeof(double));
    rt = (double *) malloc(ntheta * sizeof(double));

    // STAGE 2 INITIALIZATION
    asub = (double **) malloc(ntheta * sizeof(double));
    apri = (double **) malloc(ntheta * sizeof(double));
    asup = (double **) malloc(ntheta * sizeof(double));
    b1 = (double *) malloc(nphi * sizeof(double));
    b2 = (double *) malloc(nphi * sizeof(double));
    b3 = (double *) malloc(nphi * sizeof(double));
    b4 = (double *) malloc(nphi * sizeof(double));
    b5 = (double *) malloc(nphi * sizeof(double));
    rp = (double *) malloc(nphi * sizeof(double));

    for ( i = 1 ; i < ntheta - 1 ; i++ ) {
        th = (i - 0.5) * dth;
        cos_th = cos(th);
        sin_th = sin(th);

        // linear term
        coefa = (grad[i] + flow[i] * cos_th / sin_th) / field_rad;
        coefa = dt / 2 * (coefa - 1 / field_tau);

        // first deriv theta
        coefb = flow[i] + field_eta * cos_th / sin_th / field_rad;
        coefb *= dt / 2 / field_rad / dth;

        // first deriv phi
        coefd = - dt * difr[i] / 2 / dph;

        // second deriv phi
        coefe = dt * field_eta / pow(field_rad * sin_th * dph, 2);

        // construct matrix M
        msub[i] = -alpha * (coefc - coefb);
        mpri[i] = 1 - alpha * (coefa - 2 * coefc); // TODO check 2 on coefa
        msup[i] = -alpha * (coefc + coefb);

        // construct vector q
        q1[i] = (1 - alpha) * (coefc - coefb);
        q2[i] = 1 + coefa + (1 - alpha) * (coefa - 2 * coefc) - 2 * coefe;
        q3[i] = (1 - alpha) * (coefc + coefb);
        q4[i] = coefe - coefd;
        q5[i] = coefe + coefd;

        // construct matrix A
        asub[i] = (double *) malloc(nphi * sizeof(double));
        apri[i] = (double *) malloc(nphi * sizeof(double));
        asup[i] = (double *) malloc(nphi * sizeof(double));

        for ( j = 0 ; j < nphi ; j++ ) {

            // intermediate grid
            grids1[j] = (double *) malloc(ntheta * sizeof(double));

            asub[i][j] = -alpha * (coefe - coefd);
            apri[i][j] = 1 - alpha * (coefa - 2 * coefe);
            asup[i][j] = -alpha * (coefe + coefd);
        }

        // construct vector b
        b1[i] = (1 - alpha) * (coefe - coefd);
        b2[i] = 1 + coefa + (1 - alpha) * (coefa - 2 * coefe) - 2 * coefc;
        b3[i] = (1 - alpha) * (coefe + coefd);
        b4[i] = coefc - coefb;
        b5[i] = coefc + coefb;

    }
}

void ftcs(double **grid, double **newgrid, double *flow, double *grad, double *difr,
        int ntheta, int nphi, double dt) {

    int i, j;
    double val, c1, c2, c3, c4, c5;

    for ( i = 1 ; i < ntheta - 1; i++ ) {
        
        //c1 = 1 + (2 * (coefa[i] - coefc - coefe[i]));
        // TODO check math for c1
        c1 = (q2[i] + b2[i]) / 2;
        c2 = b5[i];
        c3 = b4[i];
        c4 = q5[i];
        c5 = q4[i];

        for ( j = 0 ; j < nphi ; j++ ) {
        
            val = c1 * grid[i][j];
            val += c2 * grid[i+1][j];
            val += c3 * grid[i-1][j];
            val += c4 * grid[i][(j+1)%nphi];
            val += c5 * grid[i][(j-1+nphi)%nphi];
            newgrid[i][j] = val;
        }
    }

    // Neumann boundary conditions - zero flux at poles
    for ( j = 0 ; j < nphi ; j++ ) {
        newgrid[0][j] = newgrid[1][j];
        newgrid[ntheta-1][j] = newgrid[ntheta-2][j];
    }
}

void ftcs_tri(double **grid, double **newgrid, double *flow, double *grad, double *difr,
        int ntheta, int nphi, double dt) {

    int i, j;

    // STAGE 1 -- theta
    for ( j = 0 ; j < nphi ; j++ ) {
        for ( i = 1 ; i < ntheta - 1 ; i++ ) {
            rt[i] = q1[i] * grid[i-1][j] +
                    q2[i] * grid[i][j] +
                    q3[i] * grid[i+1][j] +
                    q4[i] * grid[i][(j-1+nphi)%nphi] +
                    q5[i] * grid[i][(j+1)%nphi];
        }

        // just inline the whole tridiag stuff...
        // eh maybe not
        tridiag(msub+1, mpri+1, msup+1, rt+1, grids1[j]+1, ntheta-2);
    }

    // Neumann BCs for theta
    for ( j = 0 ; j < nphi ; j++ ) {
        grids1[0][j] = grids1[1][j];
        grids1[ntheta-1][j] = grids1[ntheta-2][j];
    }

    // STAGE 2 -- phi
    for ( i = 1 ; i < ntheta - 1 ; i++ ) {
        for ( j = 0 ; j < nphi ; j++ ) {

            rp[j] = b1[i] * grids1[j][i-1] +
                    b2[i] * grids1[j][i] +
                    b3[i] * grids1[j][i+1] +
                    b4[i] * grids1[(j-1+nphi)%nphi][i] +
                    b5[i] * grids1[(j+1)%nphi][i];
        }
        // just inline the whole tridiag stuff...
        // eh maybe not
        tridiag_c(asub[i], apri[i], asup[i], rp, newgrid[i], nphi);
        // when we place new values into the grid, we can't use old ones right?
        // or we can at least have a separate scheme that allows but for now sep
    }

    // Neumann BCs for theta
    for ( j = 0 ; j < nphi ; j++ ) {
        newgrid[0][j] = newgrid[1][j];
        newgrid[ntheta-1][j] = newgrid[ntheta-2][j];
    }

}
