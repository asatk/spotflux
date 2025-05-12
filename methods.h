typedef double **(*method_t)(double **, double **, double *, double *, double *,
        int, int, double);

double **ftcs(double **grid, double **newgrid, double *flow, double *grad, double *difr, int ntheta, int nphi, double dt);
double **lax(double **grid, double **newgrid, double *flow, double *grad, double *difr, int ntheta, int nphi, double dt);

method_t methods[] = {
    ftcs
};
