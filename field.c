#include <math.h>
#include <stdlib.h>

#include "field.h"


/**
 * Initial radial field profile 1
 */
inline double init_field_1(double th, double ph) {
    (void)ph;
    return field_b0 * abs(pow(cos(th),7)) * cos(th);
}

/**
 * Initial radial field profile 2
 */
inline double init_field_2(double th, double ph) {
    (void)ph;
    return field_b0 * erf(8 * abs(pow(cos(th),11)) * cos(th) / M_PI);
}
