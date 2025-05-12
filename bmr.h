double **inject_bmr(double **grid, int ntheta, int nphi);

static const double bmr_b0 = 10.0e-4;  // strength of BMR (T)
static const double bmr_sigma = 4. / 180 * M_PI;  // spread of BMR
/**
 * Hard-coded location of pos and neg BMR
 */
static const double bmrp_th = 2 * M_PI / 3 + bmr_sigma;
static const double bmrn_th = 2 * M_PI / 3 - bmr_sigma;
static const double bmrp_ph = M_PI + bmr_sigma;
static const double bmrn_ph = M_PI - bmr_sigma;
