typedef void(*emerge_t)(double **, int, int, double, double);

void none(double **grid, int ntheta, int nphi, double time, double dt);
void naive(double **grid, int ntheta, int nphi, double time, double dt);
void lemerle(double **grid, int ntheta, int nphi, double time, double dt);
void schrijver(double **grid, int ntheta, int nphi, double time, double dt);
void set_activity(double act);
void set_bmr_freq(double bmr_f);

static const emerge_t emerge_modes[] = {
    none,
    naive,
    schrijver,
    lemerle
};

////// SCHRIJVER
static const double p_bmr = 1.9;
static const double p_eph = 2.9;
static const double A0 = 1.0;   // activity factor (rel. to Sun at max)
static const double beta = 1.0; // normalization of activity variation
static const double tcycle = 11 * 3.15e7;   // solar cycle 11yrs in s
static const double a0 = 8.0 / 86400 * 180 * 180 / M_PI / M_PI;       // power law const (rad^-2 s^-1 per hemisphere)
static const double a1 = 8.0 / 86400 * 180 * 180 / M_PI / M_PI;       // power law const (s)
static const double fluxmax = 1.5e22;       // max BMR flux (Mx)
static const double fluxmin = 1.2e19;       // min BMR flux (Mx)
static const double avgfluxd = 180;         // avg flux density (field) in G
static const double smax = fluxmax / avgfluxd; // in cm^2
static const double smin = fluxmin / avgfluxd; // in cm^2

// COLATITUDE
static const double th0 = (90 - 25) / 180 * M_PI;
static const double th1 = 4 / 180 * M_PI;
static const double mu_th = (90 - 17.5) / 180 * M_PI;
static const double flux_th = 5e20;

// INCLINATION
static const double i0 = M_PI / 2;
static const double i1 = M_PI / 10;
static const double mu_i = 4.2 / 180 * M_PI;
static const double flux_i = 8e19;



////// NAIVE
/**
 * Hard-coded location of pos and neg naive BMR
 */
static const double bmr_b0 = 10.0;  // strength of BMR (G)
static const double bmr_sigma = 4. / 180 * M_PI;  // spread of BMR
static const double bmr_th = M_PI / 3;
static const double bmr_ph = M_PI;
static const double bmr_i = 4.2 / 180 * M_PI;
static const double bmr_sep = 1.8e9 / 2;
