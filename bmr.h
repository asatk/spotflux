typedef void(*emerge_t)(double **, int, int, double, double);

void none(double **grid, int ntheta, int nphi, double time, double dt);
void naive(double **grid, int ntheta, int nphi, double time, double dt);
void lemerle(double **grid, int ntheta, int nphi, double time, double dt);
void schrijver(double **grid, int ntheta, int nphi, double time, double dt);
void set_activity(double act);

static const emerge_t emerge_modes[] = {
    none,
    naive,
    schrijver,
    lemerle
};
