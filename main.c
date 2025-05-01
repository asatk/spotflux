#include <stdio.h>
#include <stdlib.h>

#include "flow.h"
#include "field.h"
#include "init.h"
#include "methods.h"
#include "io.h"

/* GLOBAL VARIABLES */
int ntheta = 128;
int nphi = 256;
int nt = 1000;
int dt = 1 / 100 * 3.15e7;   // 100 time steps per yr
char bprof = 0;
char rartype = 0;
int nrar = 0;
method_t update = ftcs;
int freq = 50;

char *fname = "bfld.dat";


int main(int argc, char **argv) {

    int t;
    double **grid, **newgrid, **tempgrid, *flow, *grad, *difr;
    FILE *f;
    method_t update;
//    args a;

    // open file for saving data
    f = fopen(fname, "w");

    // Initialize grid
    grid = init_grid(ntheta, nphi, bprof);
    newgrid = init_grid(ntheta, nphi, bprof);


    // pre-calculate motion of charges on surface - steady, axisymmetric flow.
    flow = calc_flow(ntheta);
    grad = calc_flow_grad(ntheta);
    difr = calc_difr(ntheta);

    update = ftcs;

    // evolve surface magnetic field over time
    for ( t = 0 ; t < nt ; t++ ) {

        // save snapshot of solution
        if ( t % freq == 0 ) {
            printf("Storing solution for step %*d\n", 4, t);
            store(grid, f, ntheta, nphi);
        }

        /**
        // inject active region
        // TODO linked list
        for ( k = 0 ; k < nrar ; k++ ) {
            if (inject_t[k] == t) {
                grid = inject(grid, a.rartype);
            }
        }
        */

        tempgrid = update(grid, newgrid, flow, grad, difr, ntheta, nphi, dt);
        newgrid = grid;
        grid = tempgrid;

    }

    printf("Simulation terminated.\n");
    fclose(f);

}
