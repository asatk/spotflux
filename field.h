typedef double (*initb_t)(double, double);

double init_field_1(double th, double ph);
double init_field_2(double th, double ph);

static const initb_t init_fields[] = {
    init_field_1,
    init_field_2
};
