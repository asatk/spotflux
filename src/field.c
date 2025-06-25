#include <math.h>
#include <stdlib.h>

#include "constants.h"
#include "field.h"


/**
 * Initial surface field profile 1
 */
inline double init_field_1(double th, double ph) {
    (void)ph;
    return field_b0 * pow(fabs(cos(th)),7) * cos(th);
}

/**
 * Initial surface field profile 2
 */
inline double init_field_2(double th, double ph) {
    (void)ph;
    return field_b0 * erf(8 * pow(fabs(cos(th)),11) * cos(th) / M_PI);
}
