#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "random.h"
#include "linalg.h"
#include "bmr.h"
#include "flow.h"
#include "field.h"
#include "init.h"
#include "methods.h"
#include "io.h"

/* GLOBAL VARIABLES */
int ntheta = 128;
int nphi = 256;
int nt = 10001;
double dt = 3e1;   // 1e6 steps per year
char bprof = 1;
int freq = 100;

// UPDATE
method_t update = ftcs;
double alpha = 0.0;

// RANDOM
unsigned long long seed = 0x2025LL;

// BMR
emerge_t emerge = naive;
double activity = -1.0;
double bmr_freq = 3.15e7 / 365; // 1 bmr every week

char *fname = "bfld.dat";


int main(int argc, char **argv) {

    int t;
    double **grid, **newgrid, **tempgrid;
    double *flow, *grad, *difr;
    double dth, dph, cfl, time;

    FILE *f;

    // open file for saving data
    f = fopen(fname, "w");

    // Initialize PRNG
    set_seed(seed);
    set_activity(activity);
    set_bmr_freq(bmr_freq);

    // Initialize tridiagonal solvers
    init_solvers(ntheta, nphi);

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

    if ( alpha == 0.0 )
        update = ftcs;

    if ( update == ftcs ) {
        alpha = 0.0;
        init_ftcs(flow, grad, difr, ntheta, nphi, dt, alpha);
    } else if ( update == ftcs_tri)
        init_ftcs(flow, grad, difr, ntheta, nphi, dt, alpha);

    // evolve surface magnetic field over time
    time = 0;
    for ( t = 0 ; t < nt ; t++ ) {

        // inject active region
        emerge(grid, ntheta, nphi, time, dt);

        // save snapshot of solution
        if ( t % freq == 0 ) {
            printf("Storing solution for step %*d\n", 8, t);
            store(grid, f, ntheta, nphi);
        }

        update(grid, newgrid, flow, grad, difr, ntheta, nphi, dt);
        tempgrid = newgrid;
        newgrid = grid;
        grid = tempgrid;
        time += dt;

    }

    printf("Simulation terminated.\n");
    fclose(f);

}
