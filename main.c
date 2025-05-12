#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "flow.h"
#include "field.h"
#include "bmr.h"
#include "init.h"
#include "methods.h"
#include "io.h"

/* GLOBAL VARIABLES */
int ntheta = 128;
int nphi = 256;
int nt = 1000001;
double dt = 3.15e7 / 500000;   // 5000 time steps per yr
char bprof = 1;
char rartype = 0;
int nrar = 0;
method_t update = ftcs;
int freq = 5000;
int bmr_step = 50000;

char *fname = "bfld.dat";


int main(int argc, char **argv) {

    int t;
    double **grid, **newgrid, **tempgrid, *flow, *grad, *difr, dth, dph, cfl;
    FILE *f;
//    args a;

    // open file for saving data
    f = fopen(fname, "w");

    // Initialize grid
    grid = init_grid(ntheta, nphi, bprof);
    newgrid = init_grid(ntheta, nphi, bprof);

    dth = M_PI / (ntheta - 2);
    dph = 2 * M_PI / (nphi - 1);
    cfl = field_eta * dt * (1 / dth / dth + 1 / dph / dph);
    printf("CFL est: %.4le\n", cfl);

    // pre-calculate motion of charges on surface - steady, axisymmetric flow.
    flow = calc_flow(ntheta);
    grad = calc_flow_grad(ntheta);
    difr = calc_difr(ntheta);

    // evolve surface magnetic field over time
    for ( t = 0 ; t < nt ; t++ ) {


        // inject active region
        if ( t % bmr_step == 0 ) {
            grid = inject_bmr(grid, ntheta, nphi);
        }

        // save snapshot of solution
        if ( t % freq == 0 ) {
            printf("Storing solution for step %*d\n", 8, t);
            store(grid, f, ntheta, nphi);
        }

        tempgrid = update(grid, newgrid, flow, grad, difr, ntheta, nphi, dt);
        newgrid = grid;
        grid = tempgrid;

    }

    printf("Simulation terminated.\n");
    fclose(f);

}
