#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "src/constants.h"
#include "src/random.h"
#include "src/linalg.h"
#include "src/bmr.h"
#include "src/flow.h"
#include "src/field.h"
#include "src/init.h"
#include "src/methods.h"
#include "src/io.h"

method_t update = ftcs_tri;
emerge_t emerge = schrijver;

int main(int argc, char **argv) {

    int t;
    double **grid, **newgrid, **tempgrid;
    double *flow, *grad, *difr;
    double cfl, time;

    FILE *f;

    parseargs(argc, argv);

    // open file for saving data
    f = fopen(fname, "w");

    // Initialize PRNG
    set_seed(seed);
    set_activity(activity);

    // Initialize tridiagonal solvers
    init_solvers();

    // Initialize grid
    grid = init_grid(bprof);
    newgrid = init_grid(bprof);

    // TODO separate emergence into diff submodules
    // init spot emergence
    if ( emerge == schrijver )
        init_schrijver(activity);

    cfl = field_eta * dt * (1 / dth / dth + 1 / dph / dph);
    printf("CFL est: %.4le\n", cfl);

    // pre-calculate motion of charges on surface - steady, axisymmetric flow.
    flow = calc_flow();
    grad = calc_flow_grad();
    difr = calc_difr();

    if ( alpha == 0.0 )
        update = ftcs;

    if ( update == ftcs ) {
        alpha = 0.0;
        init_ftcs(flow, grad, difr);
    } else if ( update == ftcs_tri)
        init_ftcs(flow, grad, difr);

    printf("Simulation initialized.\n");
    // evolve surface magnetic field over time
    time = 0;
    for ( t = 0 ; t < nt ; t++ ) {

        // inject active region
        emerge(grid, time);

        // save snapshot of solution
        if ( t % freq == 0 ) {
            printf("Storing solution for step %*d\n", 8, t);
            store(grid, f);
        }

        update(grid, newgrid, flow, grad, difr);
        tempgrid = newgrid;
        newgrid = grid;
        grid = tempgrid;
        time += dt;

    }

    printf("Simulation terminated.\n");
    fclose(f);

}
