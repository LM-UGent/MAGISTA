#ifndef DEF_RANDGEN
#define DEF_RANDGEN

#include "stdint.h"

typedef struct randgen_s randgen_s;

randgen_s* randgen_alloc();
void randgen_free(randgen_s* generator);
void randgen_unseed(randgen_s* generator);
void randgen_seed(randgen_s* generator, uint64_t value);
int64_t randgen_uniform_i64(randgen_s* generator, int64_t minvalue, int64_t maxvalue);
double randgen_uniform_f64(randgen_s* generator, double minvalue, double maxvalue);

#endif

