typedef double **(*method_t)(double **, double *, double *, double *,
        double, double, double);

double **ftcs(double **grid, double *flow, double *grad, double *difr, double ntheta, double nphi, double dt);

method_t methods[] = {
    ftcs
};
