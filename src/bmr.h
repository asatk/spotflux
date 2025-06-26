typedef void(*emerge_t)(double **, double);

void none(double **grid, double time);
void naive(double **grid, double time);
void lemerle(double **grid, double time);
void schrijver(double **grid, double time);
void set_activity(double act);

static const emerge_t emerge_modes[] = {
    none,
    naive,
    schrijver,
    lemerle
};
