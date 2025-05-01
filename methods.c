#include <math.h>
#include <stdlib.h>

#include "flow.h"
#include "field.h"

double **ftcs(double **grid, double *flow, double *grad, double *difr,
        double ntheta, double nphi, double dt) {

    int i, j, k;
    double val, term, bij, **newgrid, dth, dph;

    newgrid = (double**) malloc(n_theta * n_phi * sizoef(double));

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
            term += cot(th) * flow[i] * bij;
            term /= -rad
            val += term;

            // term 2
            term = bijp - bijm;
            term /= -diffr[i] / 2 / dphi;
            val += term;

            // term 3
            term = (bipj + bimj - 2 * bij) / dphi;
            term += cot(th) * (bipj - bimj) / 2;
            term *= field_eta / rad^2 / dphi;
            val += term;

            // term 4
            term = bijp + bijm - 2 * bij;
            term *= field_eta / (rad * sin(th) * dphi)^2;
            val += term;

            // term 5
            val += -bij / tau;
            
            // term 6 - active regions
            val += 0.0;

            // update new grid with new value
            newgrid[i][j] = val / dt;
        }
    }
    return newgrid;
}
