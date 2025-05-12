#include <math.h>

#include "bmr.h"

double **inject_bmr(double **grid, int ntheta, int nphi) {
    int i, j;
    double distp, distn, val, th, ph, dth, dph;

    (void)ph;

    dth = M_PI / (ntheta - 2);
    dph = 2 * M_PI / (nphi - 1);

    for ( i = 0 ; i < ntheta ; i++ ) {
        th = (i - 0.5) * dth;
        for ( j = 0 ; j < nphi ; j++ ) {
            ph = j * dph;
            distp = sqrt(pow(th - bmrp_th, 2) + pow(ph - bmrp_ph, 2));
            distn = sqrt(pow(th - bmrn_th, 2) + pow(ph - bmrn_ph, 2));
            val = exp(-pow(distp / bmr_sigma, 2) / 2) - exp(-pow(distn / bmr_sigma, 2) / 2);
            val *= bmr_b0;
            grid[i][j] += val;
        }
    }
    return grid;
}
