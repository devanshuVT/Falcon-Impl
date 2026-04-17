#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define N 512

typedef struct {
    int coeffs[N];
} poly;

static const int q = 12289;
static const double sigma_min = 1.277833697;
static const double sigma_max = 1.8205;

uint8_t get_random_byte(void) {
    return (uint8_t)arc4random_uniform(256);
}

//rough 72-bit container: 64 high bits + 8 low bits
typedef struct {
    uint64_t hi;
    uint8_t lo;
} uint72_t;

// compare two 72-bit values: return 1 if a < b else 0
int compare_uint72_less(uint72_t a, uint72_t b) {
    if (a.hi < b.hi) return 1;
    if (a.hi > b.hi) return 0;
    if (a.lo < b.lo) return 1;
    return 0;
}

// make one 72-bit random value
uint72_t get_random_72bits(void) {
    uint72_t u;
    u.hi = 0;
    for (int i = 0; i < 8; i++) {
        u.hi = (u.hi << 8) | get_random_byte();
    }
    u.lo = get_random_byte();
    return u;
}

static const uint72_t RCDT[18] = {
    {11815180629386737944ULL,   2},
    { 6112276501875359195ULL, 130},
    { 2485368865086253505ULL, 255},
    {  779533143144634698ULL, 228},
    {  186200561932255007ULL, 111},
    {   33577742212363453ULL,  95},
    {    4544132645877610ULL, 228},
    {     459595263094115ULL, 218},
    {      34638249229156ULL,  40},
    {       1941286552588ULL, 105},
    {         80784707633ULL, 251},
    {          2493483785ULL,  31},
    {            57040297ULL, 152},
    {              966510ULL, 187},
    {               12125ULL, 126},
    {                 112ULL, 152},
    {                   0ULL, 198},
    {                   0ULL,   1}
};

uint64_t ApproxExp(double x, double ccs)
{
    // Constants C copied from Falcon spec Algorithm 13; ULL added for C uint64_t literals
    static const uint64_t C[13] = {
        0x00000004741183A3ULL,
        0x00000036548CFC06ULL,
        0x0000024FDCBF140AULL,
        0x0000171D939DE045ULL,
        0x0000D00CF58F6F84ULL,
        0x000680681CF796E3ULL,
        0x002D82D8305B0FEAULL,
        0x011111110E066FD0ULL,
        0x0555555555070F00ULL,
        0x155555555581FF00ULL,
        0x400000000002B400ULL,
        0x7FFFFFFFFFFF4800ULL,
        0x8000000000000000ULL
    };
    uint64_t y = C[0];
    uint64_t z = (uint64_t)floor(ldexp(x, 63));
    for (int u = 1; u <= 12; u++) {
        __uint128_t prod = (__uint128_t)z * y;
        y = C[u] - (prod >> 63);
    }
    z = (uint64_t)floor(ldexp(ccs, 63));
    __uint128_t prod2 = (__uint128_t)z * y;
    y = (uint64_t)(prod2 >> 63);
    return y;
}

int BerExp(double x, double ccs)
{
    double ln2 = log(2.0); //computes ln⁡(2) because the algorithm uses ln(2) in both Step 1 and Step 2
    int s = (int)floor(x / ln2);
    double r = x - s * ln2;
    if (s > 63) {
        s = 63;
    }
    uint64_t z = (2 * ApproxExp(r, ccs) - 1) >> s;
    int i = 64;
    int w;
    do {
        i = i - 8;
        w = get_random_byte() - ((z >> i) & 0xFF);
    } while ((w == 0) && (i > 0));
    return (w < 0);
}

int BaseSampler(void) {
    uint72_t u = get_random_72bits();
    int z0 = 0;
    for (int i = 0; i <= 17; i++) {
        z0 = z0 + compare_uint72_less(u, RCDT[i]);
    }
    return z0;
}

int SamplerZ(double mu, double sigma_prime)
{
    double r = mu - floor(mu);
    double ccs = sigma_min / sigma_prime;
    while (1) {
        int z0 = BaseSampler();
        int b = get_random_byte() & 0x1;
        int z = b + (2 * b - 1) * z0;
        double x = ((z - r) * (z - r)) / (2.0 * sigma_prime * sigma_prime) - ((double)z0 * (double)z0) / (2.0 * sigma_max * sigma_max);
        if (BerExp(x, ccs) == 1) {
            return z + (int)floor(mu);
        }
    }
}

int main(void)
{
    double mu = 0.3;
    double sigma_prime = 1.4;

    for (int t = 0; t < 20; t++) {
        int z = SamplerZ(mu, sigma_prime);
        printf("z = %d\n", z);
    }

    return 0;
}