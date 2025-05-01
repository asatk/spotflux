#include <math.h>
#include <stdlib.h>

#include "flow.h"
#include "field.h"

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
            ph = i * dph;

            // store these values to avoid repeated computations
            bij = grid[i][j];
            bipj = grid[i+1][j];
            bimj = grid[i-1][j];
            bijp = grid[i][(j+1)%nphi];
            bijm = grid[i][(j-1)%nphi];

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
            term /= -difr[i] / 2 / dph;
            val += term;

            // term 3
            term = (bipj + bimj - 2 * bij) / dph;
            term += cos(th) / sin(th) * (bipj - bimj) / 2;
            term *= field_eta / pow(field_rad,2) / dph;
            val += term;

            // term 4
            term = bijp + bijm - 2 * bij;
            term *= field_eta / pow(field_rad * sin(th) * dph,2);
            val += term;

            // term 5
            val += -bij / field_tau;
            
            // term 6 - active regions
            val += 0.0;

            // update new grid with new value
            newgrid[i][j] = val / dt;
        }
    }

    // Neumann boundary conditions - zero flux at poles
    for ( j = 0 ; j < nphi ; j++ ) {
        newgrid[0][j] = newgrid[1][j];
        newgrid[ntheta-1][j] = newgrid[ntheta-2][j];
    }

    return newgrid;
}
