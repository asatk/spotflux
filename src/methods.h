typedef void(*method_t)(double **, double **, double *, double *, double *);

void ftcs(double **grid, double **newgrid,
        double *flow, double *grad, double *difr);
void ftcs_tri(double **grid, double **newgrid,
        double *flow, double *grad, double *difr);

void init_ftcs(double *flow, double *grad, double *difr);

static const method_t methods[] = {
    ftcs,
    ftcs_tri
};
