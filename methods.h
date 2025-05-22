typedef void(*method_t)(double **, double **, double *, double *, double *,
        int, int, double);

void ftcs(double **grid, double **newgrid, double *flow, double *grad, double *difr, int ntheta, int nphi, double dt);
void init_coef(double *flow, double *grad, double *difr,
        double ntheta, double nphi, double dt);

static const method_t methods[] = {
    ftcs
};
