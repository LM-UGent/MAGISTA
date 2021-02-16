#include "randgen.h"
#include <stdlib.h>
#include <time.h>


struct randgen_s {
    uint64_t state_1, state_2;
    int mode;
};

randgen_s* randgen_alloc() {
    randgen_s* result;
    result = (randgen_s*)calloc(1, sizeof(randgen_s));
    result->state_2 = (((uint64_t)time(NULL)) << 1) | 1; /*must be an odd number to avoid collapsing to 0*/
    result->mode = 1;
    return result;
}
void randgen_free(randgen_s* generator){
    free(generator);
}


void randgen_unseed(randgen_s* generator) {
    generator->mode = 0;
    srand((unsigned int)time(NULL)*1000 + (unsigned int)clock());
}

void randgen_seed(randgen_s* generator, uint64_t value) {
    generator->mode = 1;
    generator->state_2 = (value << 1) | 1; /*must be an odd number to avoid collapsing to 0*/
}

static inline uint64_t _multu64_rem(uint64_t A, uint64_t B, uint64_t* p_carry) {
    uint64_t Alow, Ahigh, Blow, Bhigh;
    uint64_t Rlow, Rhigh, Rmid1, Rmid2 ,tmp;
    const uint64_t u32_count = (((uint64_t)1) << 32);
    Alow = A % u32_count;
    Ahigh = (A >> 32);
    Blow = B % u32_count;
    Bhigh = (B >> 32);

    tmp = Alow * Blow;
    Rlow = tmp % u32_count;
    tmp = Ahigh* Blow + (tmp >> 32);
    Rmid1 = tmp % u32_count;
    Rmid2 = (tmp >> 32);
    tmp = Rmid1 + Alow * Bhigh;
    Rmid1 = tmp % u32_count;
    tmp = Rmid2 + Ahigh * Bhigh + (tmp >> 32);
    Rmid2 = tmp % u32_count;
    Rhigh = (tmp >> 32);

    *p_carry = (Rhigh << 32) | Rmid2;
    return (Rmid1 << 32) | Rlow;
}
static inline void _multu128(uint64_t Ahigh, uint64_t Alow, uint64_t Bhigh, uint64_t Blow, uint64_t* p_RHigh, uint64_t* p_Rlow) {
    /*
    this version assumes that all operation work for 64-bit integers, but not 128-bit integers
    If such operations exist, they should replace the code below (or above) through preprocessor conditionals
    */
    uint64_t Rlow, Rhigh;
    Rlow = _multu64_rem(Alow, Blow, &Rhigh);
    Rhigh += Blow*Ahigh + Bhigh*Alow;
    *p_RHigh = Rhigh;
    *p_Rlow = Rlow;
}

static uint64_t randgen_generate(randgen_s* generator) {
    /*
    a = 25096281518912105342191851917838718629
    = 1360472147205615982 * 2^64 + 3346542527535191717
    M = 2^128
    next = (previous*a) % M
    see: TABLES OF LINEAR CONGRUENTIAL GENERATORS OF DIFFERENT SIZES AND GOOD LATTICE STRUCTURE, L'ECUYER, 1999
    */
    uint64_t tmp1, tmp2;
    if (generator->mode == 1) {
        _multu128(generator->state_1, generator->state_2, 1360472147205615982, 3346542527535191717, &tmp1, &tmp2);
        generator->state_1 = tmp1;
        generator->state_2 = tmp2;
        return generator->state_1;
    }
    else {
        /* hardware random number generator should go here */
        return (((uint64_t)rand()) << 32) + rand();
    }
}

int64_t randgen_uniform_i64(randgen_s* generator, int64_t minvalue, int64_t maxvalue) {
    if (minvalue >= maxvalue)return minvalue;
    if (maxvalue == INT64_MAX && minvalue==INT64_MIN) {
        return (int64_t)randgen_generate(generator);
    }
    /* consider doing something more involved if you want to avoid modulo bias*/
    return randgen_generate(generator) % (maxvalue - minvalue + 1) + minvalue;
}
double randgen_uniform_f64(randgen_s* generator, double minvalue, double maxvalue) {
    return minvalue + (maxvalue-minvalue)*((double)(randgen_generate(generator)>>12))/((double)(UINT64_MAX>>12));
}