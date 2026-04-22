#ifndef XPOLY_H
#define XPOLY_H

#include <gmp.h>
#include "poly.h"   // for N

typedef struct {
    mpz_t coeffs[N];
} xpoly;

/* lifecycle */
void xpoly_init(xpoly *p);
void xpoly_clear(xpoly *p);
void xpoly_zero(xpoly *p, int n);
void xpoly_set(xpoly *dst, const xpoly *src, int n);

/* conversions / debug */
void xpoly_from_int_poly(xpoly *dst, const poly *src, int n);
void xpoly_print(const xpoly *p, int n, const char *name);

/* arithmetic */
void xpoly_add_n(xpoly *res, const xpoly *a, const xpoly *b, int n);
void xpoly_sub_n(xpoly *res, const xpoly *a, const xpoly *b, int n);

/* exact negacyclic multiplication modulo x^n + 1 */
void xpoly_mul_mod_phi_n(xpoly *res, const xpoly *a, const xpoly *b, int n);

void xpoly_substitute_neg_x(xpoly *res, const xpoly *a, int n);
void xpoly_substitute_x2(xpoly *res, const xpoly *a, int n);

void xpoly_field_norm(xpoly *res, const xpoly *a, int n);
int xpoly_is_zero_n(const xpoly *a, int n);
long xpoly_get_si_checked(const xpoly *a, int idx);
int NTRUSolve_exact(const xpoly *f, const xpoly *g, xpoly *F, xpoly *G, int n, int q);

int extended_gcd(int a, int b, int *u, int *v);
void xpoly_scalar_mul_si(xpoly *res, const xpoly *a, long c, int n);
int xpoly_equals_monomial_q(const xpoly *a, int n, int q);
int xpoly_verify_ntru(const xpoly *f, const xpoly *g, const xpoly *F, const xpoly *G, int n, int q);
#endif
