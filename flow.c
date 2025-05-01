#include <math.h>
#include <stdlib.h>

#include "flow.h"

/**
 * Meridional flow
 */
double *calc_flow(int ntheta) {
    
    int i;
    double th, dth, cos_th, sin_th, val, *flow;

    flow = (double *) malloc(ntheta * sizeof(double));

    dth = M_PI / (ntheta - 2);

    for ( i = 0 ; i < ntheta ; i++ ) {
        th = (i - 0.5) * dth;
        cos_th = cos(th);
        sin_th = sin(th);

        val = pow(erf(flow_v * sin_th), flow_q);
        val += pow(erf(flow_w * cos_th), flow_n);
        val *= -flow_u0;

        flow[i] = val;
    }

    return flow;

}

/**
 * Differential rotation
 */
double *calc_difr(int ntheta) {
    
    int i;
    double th, dth, cos_th, val, *difr;

    difr = (double *) malloc(ntheta * sizeof(double));

    dth = M_PI / (ntheta - 2);

    for ( i = 0 ; i < ntheta ; i++ ) {
        th = (i - 0.5) * dth;
        cos_th = cos(th);

        val = 1 + difr_a2 * pow(cos_th,2) + difr_a4 * pow(cos_th,4);
        val *= difr0;

        difr[i] = val;
    }

    return difr;

}

double *calc_flow_grad(int ntheta) {

    int i;
    double *grad, vsin, wcos, th, dth, val;

    grad = calc_flow(ntheta);

    dth = M_PI / (ntheta - 2);

    for ( i = 0 ; i < ntheta ; i++ ) {
        th = (i - 0.5) * dth;
        vsin = flow_v * sin(th);
        wcos = flow_w * cos(th);

        val = exp(-pow(vsin,2)) * flow_q * flow_v * cos(th) / erf(vsin);
        val -= exp(-pow(wcos,2)) * flow_n * flow_w * sin(th) / erf(wcos);
        val *= 2 / sqrt(M_PI);
        
        grad[i] = val;
    }

    return grad;

}
