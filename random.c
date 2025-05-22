#include <math.h>
#include <stdlib.h>

#include "random.h"

/**
 * Efficient and good PRNG from Numerical Recipes
 */

unsigned long long state = 4101842887655102017LL;

void set_seed(unsigned long long seed) {
    state ^= seed;
    state = ranq1();
}

inline unsigned long long ranq1() {
    state ^= state >> 21;
    state ^= state << 35;
    state ^= state >> 4;
    return state * 2685821657736338717LL;
}

/**
 * Use ranq1 PRNG to sample a uniform variate from [0.0, 1.0]
 */
inline double ranq1dbl() {
    return 5.42101086242752217E-20 * ranq1();
}

inline double ranq1pwr(double p, double x0, double x1) {
    double y;
    y = ranq1dbl();
    return pow(((pow(x1, p+1) - pow(x0, p+1)) * y + pow(x0, p+1)), 1/(p+1));
}

double randnorm(double mu, double sig) {

    double u, v, x, y, q;

    do {
        u = ranq1dbl();
        v = 1.7156 * (ranq1dbl() - 0.5);
        x = u - 0.449871;
        y = abs(v) + 0.386595;
        q = sqrt(x) + y * (0.19600 * y - 0.25472 * x);
    } while ( q > 0.27597 &&
            (q > 0.27846 || sqrt(v) > -4 * log(u) * sqrt(u)));
    return mu + sig * v / u;

}
