#ifndef POLY_H
#define POLY_H
#include <stdint.h>
#define N 512

typedef struct {
    double coeffs[N];
} poly;

typedef struct {
    double real[N];
    double imag[N];
} poly_fft;

poly poly_zero(void);
int poly_is_zero(poly p);
poly poly_add(poly a, poly b);
poly poly_sub(poly a, poly b);
poly poly_scalar_mul(poly a, double c);
poly poly_adjoint(poly a);
poly poly_mul_mod_phi(poly a, poly b);
poly poly_round_coeffs(poly a);
double compute_gamma(poly f, poly g, int q);
double poly_inner(poly a, poly b);
double poly_norm(poly a);
poly_fft fft(poly a);
poly ifft(poly_fft a);
poly poly_div_fft(poly num, poly den);
poly compute_k(poly f, poly g, poly F, poly G, int n);
double poly_pair_norm(poly a, poly b);
poly poly_substitute_neg_x(poly a, int n);
int mod_q(int64_t x);
int mod_add(int a, int b);
int mod_mul(int a, int b);
int mod_pow(int a, int e);
int ntt_has_zero_coefficient(poly f);

poly poly_mul_mod_phi_n(poly a, poly b, int n);
poly poly_field_norm(poly a, int n);

poly poly_substitute_x2(poly a, int n);

poly poly_adjoint_n(poly a, int n);
poly poly_sub_n(poly a, poly b, int n);
poly poly_round_coeffs_n(poly a, int n);
int poly_is_zero_n(poly p, int n);
poly poly_div_fft_n(poly num, poly den, int n);
poly compute_k(poly f, poly g, poly F, poly G, int n);

#endif