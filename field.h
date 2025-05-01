typedef double (*initb_t)(double, double);

double init_field_1(double th, double ph);
double init_field_2(double th, double ph);

static const initb_t init_fields[] = {
    init_field_1,
    init_field_2
};

static const double field_b0 = 8.5;    // field strength (G)
static const double field_eta = 350 * 10^6;    // diffusivity (m^2/s)
static const double field_tau = 32 * 3.15e7;   // sink time (s)

static const double field_rad = 6.957e8;     // stellar radius (m)
