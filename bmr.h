typedef void(*emerge_t)(double **, int, int, double, double);

void naive(double **grid, int ntheta, int nphi, double time, double dt);
void lemerle(double **grid, int ntheta, int nphi, double time, double dt);
void schrijver(double **grid, int ntheta, int nphi, double time, double dt);
void set_activity(double act);

static const emerge_t emerge_modes[] = {
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
static const double a0 = 8.0 / 86400;       // power law const (s)
static const double a1 = 8.0 / 86400;       // power law const (s)
static const double fluxmax = 1.5e22;       // max BMR flux
static const double fluxmin = 1.2e19;       // min BMR flux
static const double avgfluxd = 180;         // avg flux density (field) in G
static const double smax = fluxmax / avgfluxd;
static const double smin = fluxmin / avgfluxd;

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
static const double bmrp_th = bmr_th + bmr_sigma;
static const double bmrn_th = bmr_th - bmr_sigma;
static const double bmrp_ph = bmr_ph + bmr_sigma;
static const double bmrn_ph = bmr_ph - bmr_sigma;
