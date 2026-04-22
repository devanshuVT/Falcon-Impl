#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "poly.h"
#include "xpoly.h"

static const int q = 12289;
static const double sigma_min = 1.277833697;
static const double sigma_max = 1.8205;

static const uint8_t *replay_bytes = NULL;
static size_t replay_len = 0;
static size_t replay_pos = 0;

void set_replay_bytes(const uint8_t *buf, size_t len) {
    replay_bytes = buf;
    replay_len = len;
    replay_pos = 0;
}

void clear_replay_bytes(void) {
    replay_bytes = NULL;
    replay_len = 0;
    replay_pos = 0;
}

uint8_t get_random_byte(void) {
    if (replay_bytes != NULL) {
        if (replay_pos >= replay_len) {
            fprintf(stderr, "Replay byte stream exhausted\n");
            exit(1);
        }
        return replay_bytes[replay_pos++];
    }
    return (uint8_t)arc4random_uniform(256);
}

typedef struct {
    uint64_t hi;
    uint8_t lo;
} uint72_t;

int compare_uint72_less(uint72_t a, uint72_t b) {
    if (a.hi < b.hi) return 1;
    if (a.hi > b.hi) return 0;
    if (a.lo < b.lo) return 1;
    return 0;
}

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
    double ln2 = log(2.0);
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
static int BaseSampler_from_u(uint72_t u)
{
    int z0 = 0;
    for (int i = 0; i <= 17; i++) {
        z0 += compare_uint72_less(u, RCDT[i]);
    }
    return z0;
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
void Reduce(poly f, poly g, poly *F, poly *G, int n)
{
    poly k;
    int iter = 0;

    do {
        k = compute_k(f, g, *F, *G, n);

        poly kf = poly_mul_mod_phi_n(k, f, n);
        poly kg = poly_mul_mod_phi_n(k, g, n);

        *F = poly_sub_n(*F, kf, n);
        *G = poly_sub_n(*G, kg, n);

        iter++;
        if (iter > 40) {
            printf("Reduce stopped after 20 iterations\n");
            break;
        }
    } while (!poly_is_zero_n(k, n));
}

int extended_gcd(int a, int b, int *u, int *v)
{
    if (b == 0) {
        *u = 1;
        *v = 0;
        return a;
    }
    int u1, v1;
    int d = extended_gcd(b, a % b, &u1, &v1);
    *u = v1;
    *v = u1 - (a / b) * v1;
    return d;
}

int sample_fg_coeff(void)
{
    int z = 0;
    double sigma_star = 1.43300980528773;  /* spec Eq. (3.29) for Falcon-512 */

    for (int i = 0; i < 8; i++) {
        z += SamplerZ(0.0, sigma_star);
    }

    return z;
}

int NTRUSolve(poly f, poly g, poly *F, poly *G, int n)
{
    printf("Entering NTRUSolve with n = %d\n", n);

    if (n == 1) {
        int f0 = (int)llround(f.coeffs[0]);
        int g0 = (int)llround(g.coeffs[0]);
        int u, v;
        int d = extended_gcd(f0, g0, &u, &v);
        int v_spec = -v;

        if (d != 1 && d != -1) {
            printf("Base case failed at n = %d, f0 = %d, g0 = %d, d = %d\n", n, f0, g0, d);
            return 0;
        }

        if (d == -1) {
            u = -u;
            v_spec = -v_spec;
        }

        *F = poly_zero();
        *G = poly_zero();

        F->coeffs[0] = (double)(v_spec * q);
        G->coeffs[0] = (double)(u * q);

        return 1;
    } else {
        poly fp = poly_field_norm(f, n);
        poly gp = poly_field_norm(g, n);

        if (n <= 8) {
            printf("After norm, n = %d\n", n);

            printf("fp: ");
            for (int i = 0; i < n / 2; i++) {
                printf("%.0f ", fp.coeffs[i]);
            }
            printf("\n");

            printf("gp: ");
            for (int i = 0; i < n / 2; i++) {
                printf("%.0f ", gp.coeffs[i]);
            }
            printf("\n");
        }

        poly Fp = poly_zero();
        poly Gp = poly_zero();

        if (!NTRUSolve(fp, gp, &Fp, &Gp, n / 2)) {
            printf("Recursive child failed at n = %d\n", n);
            return 0;
        }

        poly g_negx = poly_substitute_neg_x(g, n);
        poly f_negx = poly_substitute_neg_x(f, n);

        poly Fp_x2 = poly_substitute_x2(Fp, n);
        poly Gp_x2 = poly_substitute_x2(Gp, n);

        poly F_x2 = poly_mul_mod_phi_n(Fp_x2, g_negx, n);
        poly G_x2 = poly_mul_mod_phi_n(Gp_x2, f_negx, n);

        *F = F_x2;
        *G = G_x2;

        Reduce(f, g, F, G, n);

        return 1;
    }
}

void dump_poly_to_file(const char *filename, poly p, int n)
{
    FILE *fp = fopen(filename, "w");
    if (fp == NULL) {
        printf("Could not open %s for writing\n", filename);
        return;
    }

    for (int i = 0; i < n; i++) {
        fprintf(fp, "%.0f\n", p.coeffs[i]);
    }

    fclose(fp);
}

int NTRUGen(poly *f, poly *g, poly *F, poly *G)
{
    int attempts = 0;

    while (1) {
        attempts++;

        *f = poly_zero();
        *g = poly_zero();
        *F = poly_zero();
        *G = poly_zero();

        for (int i = 0; i < N; i++) {
            f->coeffs[i] = (double)sample_fg_coeff();
            g->coeffs[i] = (double)sample_fg_coeff();
        }

        if (ntt_has_zero_coefficient(*f)) {
            continue;
        }

        double gamma = compute_gamma(*f, *g, q);
        double gamma_threshold = 1.17 * sqrt((double)q);

        if (gamma > gamma_threshold) {
            printf("Attempt %d failed: gamma = %.6f\n", attempts, gamma);
            continue;
        }

        printf("Attempt %d passed gamma: gamma = %.6f\n", attempts, gamma);

        if (!NTRUSolve(*f, *g, F, G, N)) {
            printf("Attempt %d failed: NTRUSolve\n", attempts);
            dump_poly_to_file("f_fail.txt", *f, N);
            dump_poly_to_file("g_fail.txt", *g, N);
            printf("Saved full failing sample to f_fail.txt and g_fail.txt\n");
            return 0;
        }

        printf("Success after %d attempt(s)\n", attempts);
        return 1;
    }
}


// int main(void)
// {
//     poly f = poly_zero();
//     poly g = poly_zero();
//     poly F = poly_zero();
//     poly G = poly_zero();

//     FILE *ff = fopen("f_fail.txt", "r");
//     FILE *fg = fopen("g_fail.txt", "r");

//     if (ff == NULL || fg == NULL) {
//         printf("Could not open failing sample files\n");
//         return 1;
//     }

//     for (int i = 0; i < N; i++) {
//         fscanf(ff, "%lf", &f.coeffs[i]);
//         fscanf(fg, "%lf", &g.coeffs[i]);
//     }

//     fclose(ff);
//     fclose(fg);

//     int ok = NTRUSolve(f, g, &F, &G, N);

//     printf("NTRUSolve ok = %d\n", ok);

//     return 0;
// }


static void load_fail_sample_as_xpoly(xpoly *f, xpoly *g) {
    FILE *ff = fopen("f_fail.txt", "r");
    FILE *fg = fopen("g_fail.txt", "r");

    if (ff == NULL || fg == NULL) {
        fprintf(stderr, "Could not open failing sample files\n");
        exit(1);
    }

    for (int i = 0; i < N; i++) {
        long fv, gv;
        if (fscanf(ff, "%ld", &fv) != 1 || fscanf(fg, "%ld", &gv) != 1) {
            fprintf(stderr, "Failed reading coefficient %d\n", i);
            exit(1);
        }
        mpz_set_si(f->coeffs[i], fv);
        mpz_set_si(g->coeffs[i], gv);
    }

    fclose(ff);
    fclose(fg);
}
static poly sample_small_poly(void) {
    poly p = poly_zero();
    for (int i = 0; i < N; i++) {
        p.coeffs[i] = (double)sample_fg_coeff();
    }
    return p;
}

static void print_mpz_bits(const xpoly *p, int n, const char *name) {
    size_t max_bits = 0;
    for (int i = 0; i < n; i++) {
        size_t b = mpz_sizeinbase(p->coeffs[i], 2);
        if (mpz_sgn(p->coeffs[i]) == 0) b = 0;
        if (b > max_bits) max_bits = b;
    }
    printf("%s max_bits = %zu\n", name, max_bits);
}

static int sample_and_solve_exact_once(int trial, int q_local) {
    poly f_d = sample_small_poly();
    poly g_d = sample_small_poly();

    if (ntt_has_zero_coefficient(f_d) || ntt_has_zero_coefficient(g_d)) {
        printf("[trial %d] skipped: NTT precheck failed\n", trial);
        return -1;
    }

    xpoly f, g, F, G;
    xpoly_init(&f);
    xpoly_init(&g);
    xpoly_init(&F);
    xpoly_init(&G);

    xpoly_zero(&f, N);
    xpoly_zero(&g, N);
    xpoly_zero(&F, N);
    xpoly_zero(&G, N);

    xpoly_from_int_poly(&f, &f_d, N);
    xpoly_from_int_poly(&g, &g_d, N);

    int ok = NTRUSolve_exact(&f, &g, &F, &G, N, q_local);
    if (!ok) {
        printf("[trial %d] exact solve failed\n", trial);
        xpoly_clear(&f);
        xpoly_clear(&g);
        xpoly_clear(&F);
        xpoly_clear(&G);
        return 0;
    }

    int verify_ok = xpoly_verify_ntru(&f, &g, &F, &G, N, q_local);
    printf("[trial %d] solve=%d verify=%d\n", trial, ok, verify_ok);

    if (verify_ok) {
        print_mpz_bits(&F, N, "F");
        print_mpz_bits(&G, N, "G");
    }

    xpoly_clear(&f);
    xpoly_clear(&g);
    xpoly_clear(&F);
    xpoly_clear(&G);

    return verify_ok ? 1 : 0;
}

static void batch_test_exact_solver(int trials) {
    const int q_local = 12289;
    int attempted = 0;
    int skipped = 0;
    int solved = 0;
    int verified = 0;

    for (int t = 1; t <= trials; t++) {
        int r = sample_and_solve_exact_once(t, q_local);
        if (r == -1) {
            skipped++;
            continue;
        }
        attempted++;
        if (r == 1) {
            solved++;
            verified++;
        }
    }

    printf("\n=== exact batch summary ===\n");
    printf("requested = %d\n", trials);
    printf("attempted = %d\n", attempted);
    printf("skipped   = %d\n", skipped);
    printf("solved    = %d\n", solved);
    printf("verified  = %d\n", verified);
}
static uint72_t make_uint72_from_hex(const char *hex18)
{
    uint72_t u;
    char hi_str[17];
    char lo_str[3];
    for (int i = 0; i < 16; i++) {
        hi_str[i] = hex18[i];
    }
    hi_str[16] = '\0';

    lo_str[0] = hex18[16];
    lo_str[1] = hex18[17];
    lo_str[2] = '\0';

    u.hi = strtoull(hi_str, NULL, 16);
    u.lo = (uint8_t)strtoul(lo_str, NULL, 16);

    return u;
}


static void test_basesampler_vectors(void)
{
    struct {
        const char *u_hex;
        int expected_z0;
    } tests[] = {
        {"0FC5442FF043D66E91", 3},
        {"CAC64EA5450A22941E", 0},
        {"F4DA0F8D8444D1A772", 0},
        {"2C1DED26CD52AED770", 2},
    };

    int num_tests = (int)(sizeof(tests) / sizeof(tests[0]));
    int pass = 1;

    for (int i = 0; i < num_tests; i++) {
        uint72_t u = make_uint72_from_hex(tests[i].u_hex);
        int got = BaseSampler_from_u(u);

        printf("[BaseSampler test %d] u=0x%s expected=%d got=%d %s\n",
               i + 1,
               tests[i].u_hex,
               tests[i].expected_z0,
               got,
               (got == tests[i].expected_z0) ? "PASS" : "FAIL");

        if (got != tests[i].expected_z0) {
            pass = 0;
        }
    }

    if (pass) {
        printf("BaseSampler vector tests: PASS\n");
    } else {
        printf("BaseSampler vector tests: FAIL\n");
    }
}
static void test_BerExp_vector_1(void)
{
    double x = 0x1.44B03ABD4E2FAp+0;
    double ccs = 0x1.7FFECD162AE20p-1;
    uint8_t bytes[] = {0xEA};

    set_replay_bytes(bytes, sizeof(bytes));
    int got = BerExp(x, ccs);
    clear_replay_bytes();

    printf("[BerExp test 1] expected=0 got=%d %s\n",
           got,
           (got == 0) ? "PASS" : "FAIL");
}
static void test_BerExp_vector_2(void)
{
    double x = 0x1.99F8C1F3DD125p-10;
    double ccs = 0x1.7FFECD162AE20p-1;

    uint8_t bytes[] = {0x6C};

    set_replay_bytes(bytes, sizeof(bytes));
    int got = BerExp(x, ccs);
    clear_replay_bytes();

    printf("[BerExp test 2] expected=1 got=%d %s\n",
           got,
           (got == 1) ? "PASS" : "FAIL");
}
static void test_SamplerZ_vector_1(void)
{
    double mu = -0x1.6F9E6CB3119A4p+6;
    double sigma_prime = 1.0 / (0x1.2C8142A489B3Cp-1);

    uint8_t bytes[] = {
        0x0F, 0xC5, 0x44, 0x2F, 0xF0, 0x43, 0xD6, 0x6E, 0x91,
        0xD1,
        0xEA,

        0xCA, 0xC6, 0x4E, 0xA5, 0x45, 0x0A, 0x22, 0x94, 0x1E,
        0xDC,
        0x6C
    };

    set_replay_bytes(bytes, sizeof(bytes));
    int got = SamplerZ(mu, sigma_prime);
    clear_replay_bytes();

    printf("[SamplerZ test 1] expected=-92 got=%d %s\n",
           got,
           (got == -92) ? "PASS" : "FAIL");
}
int NTRUGen_exact(poly *f_out, poly *g_out, xpoly *F_out, xpoly *G_out) {
    const int q_local = 12289;
    int attempts = 0;

    while (1) {
        attempts++;

        poly f = poly_zero();
        poly g = poly_zero();

        for (int i = 0; i < N; i++) {
            f.coeffs[i] = (double)sample_fg_coeff();
            g.coeffs[i] = (double)sample_fg_coeff();
        }

        if (ntt_has_zero_coefficient(f) || ntt_has_zero_coefficient(g)) {
            printf("[attempt %d] reject: NTT precheck failed\n", attempts);
            continue;
        }

        double gamma = compute_gamma(f, g, q_local);
        if (gamma > 1.17 * sqrt((double)q_local)) {
            printf("[attempt %d] reject: gamma too large (gamma=%f)\n", attempts, gamma);
            continue;
        }

        xpoly fx, gx, Fx, Gx;
        xpoly_init(&fx);
        xpoly_init(&gx);
        xpoly_init(&Fx);
        xpoly_init(&Gx);

        xpoly_zero(&fx, N);
        xpoly_zero(&gx, N);
        xpoly_zero(&Fx, N);
        xpoly_zero(&Gx, N);

        xpoly_from_int_poly(&fx, &f, N);
        xpoly_from_int_poly(&gx, &g, N);

        int ok = NTRUSolve_exact(&fx, &gx, &Fx, &Gx, N, q_local);
        if (!ok) {
            printf("[attempt %d] reject: NTRUSolve_exact failed\n", attempts);
            xpoly_clear(&fx);
            xpoly_clear(&gx);
            xpoly_clear(&Fx);
            xpoly_clear(&Gx);
            continue;
        }

        if (!xpoly_verify_ntru(&fx, &gx, &Fx, &Gx, N, q_local)) {
            printf("[attempt %d] reject: verification failed\n", attempts);
            xpoly_clear(&fx);
            xpoly_clear(&gx);
            xpoly_clear(&Fx);
            xpoly_clear(&Gx);
            continue;
        }

        *f_out = f;
        *g_out = g;
        xpoly_set(F_out, &Fx, N);
        xpoly_set(G_out, &Gx, N);

        xpoly_clear(&fx);
        xpoly_clear(&gx);
        xpoly_clear(&Fx);
        xpoly_clear(&Gx);

        printf("[attempt %d] success\n", attempts);
        return 1;
    }
}
int main(void) {
    poly f, g;
    xpoly F, G;

    xpoly_init(&F);
    xpoly_init(&G);
    xpoly_zero(&F, N);
    xpoly_zero(&G, N);

    int ok = NTRUGen_exact(&f, &g, &F, &G);
    printf("NTRUGen_exact ok = %d\n", ok);

    if (ok) {
        xpoly_print(&F, 8, "F first 8");
        xpoly_print(&G, 8, "G first 8");
    }

    xpoly_clear(&F);
    xpoly_clear(&G);
    return 0;
}