#include <math.h>
#include <stdlib.h>

#include "flow.h"
#include "field.h"

double **ftcs(double **grid, double *flow, double *grad, double *difr,
        double ntheta, double nphi, double dt) {

    int i, j;
    double val, term, **newgrid, th, dth, ph, dph;
    double bij, bipj, bimj, bijp, bijm;

    (void)ph;

    // hows abouts you flip flop btwn two mem areas??!! no more mallocing
    newgrid = (double**) malloc(ntheta * nphi * sizeof(double));

    dth = M_PI / (ntheta - 2);
    dph = 2 * M_PI / (nphi - 1);

    for ( i = 0 ; i < ntheta ; i++ ) {
        th = (i - 0.5) * dth;
        for ( j = 0 ; j < nphi ; j++ ) {
            ph = i * dph;
           
            // store these values to avoid repeated computations
            bij = grid[i][j];
            bipj = grid[i+1][j];
            bimj = grid[i-1][j];
            bijp = grid[i][j+1];
            bijm = grid[i][j-1];

            // new value of field at (theta[i], phi[j])
            // term 0
            val = bij / dt;

            // term 1
            term = grad[i] * bij;
            term += flow[i] * (bipj - bimj) / 2 / dth;
            term += cos(th) / sin(th) * flow[i] * bij;
            term /= -rad;
            val += term;

            // term 2
            term = bijp - bijm;
            term /= -difr[i] / 2 / dph;
            val += term;

            // term 3
            term = (bipj + bimj - 2 * bij) / dph;
            term += cos(th) / sin(th) * (bipj - bimj) / 2;
            term *= field_eta / pow(rad,2) / dph;
            val += term;

            // term 4
            term = bijp + bijm - 2 * bij;
            term *= field_eta / pow(rad * sin(th) * dph,2);
            val += term;

            // term 5
            val += -bij / field_tau;
            
            // term 6 - active regions
            val += 0.0;

            // update new grid with new value
            newgrid[i][j] = val / dt;
        }
    }
    return newgrid;
}
