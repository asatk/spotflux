typedef void(*method_t)(double **, double **); 

void ftcs(double **grid, double **newgrid);
void ftcs_tri(double **grid, double **newgrid);

void init_ftcs(double *flow, double *grad, double *difr);

static const method_t methods[] = {
    ftcs,
    ftcs_tri
};
