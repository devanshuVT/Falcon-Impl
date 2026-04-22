#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <limits.h>
#include "xpoly.h"
void xpoly_init(xpoly *p) {
    for (int i = 0; i < N; i++) {
        mpz_init(p->coeffs[i]);
    }
}

void xpoly_clear(xpoly *p) {
    for (int i = 0; i < N; i++) {
        mpz_clear(p->coeffs[i]);
    }
}

void xpoly_zero(xpoly *p, int n) {
    for (int i = 0; i < n; i++) {
        mpz_set_ui(p->coeffs[i], 0);
    }
}

void xpoly_set(xpoly *dst, const xpoly *src, int n) {
    for (int i = 0; i < n; i++) {
        mpz_set(dst->coeffs[i], src->coeffs[i]);
    }
}

void xpoly_from_int_poly(xpoly *dst, const poly *src, int n) {
    for (int i = 0; i < n; i++) {
        long v = (long)(src->coeffs[i]);
        mpz_set_si(dst->coeffs[i], v);
    }
}

void xpoly_print(const xpoly *p, int n, const char *name) {
    printf("%s = [", name);
    for (int i = 0; i < n; i++) {
        gmp_printf("%Zd", p->coeffs[i]);
        if (i != n - 1) printf(", ");
    }
    printf("]\n");
}

void xpoly_add_n(xpoly *res, const xpoly *a, const xpoly *b, int n) {
    for (int i = 0; i < n; i++) {
        mpz_add(res->coeffs[i], a->coeffs[i], b->coeffs[i]);
    }
}

void xpoly_sub_n(xpoly *res, const xpoly *a, const xpoly *b, int n) {
    for (int i = 0; i < n; i++) {
        mpz_sub(res->coeffs[i], a->coeffs[i], b->coeffs[i]);
    }
}

void xpoly_mul_mod_phi_n(xpoly *res, const xpoly *a, const xpoly *b, int n) {
    mpz_t *tmp;
    mpz_t prod;
    int deg;

    if (n <= 0 || n > N) {
        fprintf(stderr, "xpoly_mul_mod_phi_n: invalid n=%d\n", n);
        exit(1);
    }

    deg = 2 * n - 1;

    tmp = malloc((size_t)deg * sizeof(mpz_t));
    if (tmp == NULL) {
        fprintf(stderr, "xpoly_mul_mod_phi_n: malloc failed\n");
        exit(1);
    }

    for (int i = 0; i < deg; i++) {
        mpz_init(tmp[i]);
        mpz_set_ui(tmp[i], 0);
    }

    mpz_init(prod);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            mpz_mul(prod, a->coeffs[i], b->coeffs[j]);
            mpz_add(tmp[i + j], tmp[i + j], prod);
        }
    }
    for (int k = deg - 1; k >= n; k--) {
        if (mpz_sgn(tmp[k]) != 0) {
            mpz_sub(tmp[k - n], tmp[k - n], tmp[k]);
        }
    }

    /* write result */
    for (int i = 0; i < n; i++) {
        mpz_set(res->coeffs[i], tmp[i]);
    }

    mpz_clear(prod);

    for (int i = 0; i < deg; i++) {
        mpz_clear(tmp[i]);
    }
    free(tmp);
}


void xpoly_substitute_neg_x(xpoly *res, const xpoly *a, int n) {
    if (n <= 0 || n > N) {
        fprintf(stderr, "xpoly_substitute_neg_x: invalid n=%d\n", n);
        exit(1);
    }

    for (int i = 0; i < n; i++) {
        if ((i & 1) == 0) {
            mpz_set(res->coeffs[i], a->coeffs[i]);
        } else {
            mpz_neg(res->coeffs[i], a->coeffs[i]);
        }
    }
}

void xpoly_substitute_x2(xpoly *res, const xpoly *a, int n) {
    if (n <= 0 || n > N || (n & 1) != 0) {
        fprintf(stderr, "xpoly_substitute_x2: invalid n=%d (must be positive, <= N, and even)\n", n);
        exit(1);
    }

    int m = n / 2;

    for (int i = 0; i < n; i++) {
        mpz_set_ui(res->coeffs[i], 0);
    }

    for (int i = 0; i < m; i++) {
        mpz_set(res->coeffs[2 * i], a->coeffs[i]);
    }
}

int xpoly_is_zero_n(const xpoly *a, int n) {
    if (n <= 0 || n > N) {
        fprintf(stderr, "xpoly_is_zero_n: invalid n=%d\n", n);
        exit(1);
    }

    for (int i = 0; i < n; i++) {
        if (mpz_sgn(a->coeffs[i]) != 0) {
            return 0;
        }
    }
    return 1;
}

long xpoly_get_si_checked(const xpoly *a, int idx) {
    if (idx < 0 || idx >= N) {
        fprintf(stderr, "xpoly_get_si_checked: invalid idx=%d\n", idx);
        exit(1);
    }

    if (!mpz_fits_slong_p(a->coeffs[idx])) {
        fprintf(stderr, "xpoly_get_si_checked: coefficient does not fit in signed long\n");
        exit(1);
    }

    return mpz_get_si(a->coeffs[idx]);
}
void xpoly_field_norm(xpoly *res, const xpoly *a, int n) {
    if (n <= 0 || n > N || (n & 1) != 0) {
        fprintf(stderr, "xpoly_field_norm: invalid n=%d (must be positive, <= N, and even)\n", n);
        exit(1);
    }

    int m = n / 2;
    xpoly a_neg;
    xpoly prod;

    xpoly_init(&a_neg);
    xpoly_init(&prod);

    xpoly_zero(&a_neg, n);
    xpoly_zero(&prod, n);
    xpoly_zero(res, m);

    xpoly_substitute_neg_x(&a_neg, a, n);

    xpoly_mul_mod_phi_n(&prod, a, &a_neg, n);
    for (int i = 0; i < m; i++) {
        mpz_set(res->coeffs[i], prod.coeffs[2 * i]);
    }
    for (int i = 1; i < n; i += 2) {
        if (mpz_sgn(prod.coeffs[i]) != 0) {
            fprintf(stderr, "xpoly_field_norm: odd coefficient at index %d is nonzero\n", i);
            exit(1);
        }
    }

    xpoly_clear(&a_neg);
    xpoly_clear(&prod);
}
int NTRUSolve_exact(const xpoly *f, const xpoly *g, xpoly *F, xpoly *G, int n, int q) {
    if (n <= 0 || n > N) {
        fprintf(stderr, "NTRUSolve_exact: invalid n=%d\n", n);
        exit(1);
    }
if (n == 1) {
    mpz_t d, u, v, tmp;
    mpz_init(d);
    mpz_init(u);
    mpz_init(v);
    mpz_init(tmp);

    mpz_gcdext(d, u, v, f->coeffs[0], g->coeffs[0]);

    if (mpz_cmp_si(d, 1) != 0 && mpz_cmp_si(d, -1) != 0) {
        mpz_clear(d);
        mpz_clear(u);
        mpz_clear(v);
        mpz_clear(tmp);
        return 0;
    }
    if (mpz_cmp_si(d, -1) == 0) {
        mpz_neg(u, u);
        mpz_neg(v, v);
    }

    xpoly_zero(F, 1);
    xpoly_zero(G, 1);

    mpz_mul_si(tmp, v, q);
    mpz_neg(F->coeffs[0], tmp);
    mpz_mul_si(G->coeffs[0], u, q);
    mpz_clear(d);
    mpz_clear(u);
    mpz_clear(v);
    mpz_clear(tmp);

    return 1;
}

    if ((n & 1) != 0) {
        fprintf(stderr, "NTRUSolve_exact: n=%d must be a power-of-two style even size above base case\n", n);
        exit(1);
    }

    int m = n / 2;

    xpoly fp, gp;
    xpoly Fp, Gp;

    xpoly gx_neg, fx_neg;
    xpoly Fp_x2, Gp_x2;

    xpoly Frec, Grec;

    xpoly_init(&fp);
    xpoly_init(&gp);
    xpoly_init(&Fp);
    xpoly_init(&Gp);
    xpoly_init(&gx_neg);
    xpoly_init(&fx_neg);
    xpoly_init(&Fp_x2);
    xpoly_init(&Gp_x2);
    xpoly_init(&Frec);
    xpoly_init(&Grec);

    xpoly_zero(&fp, m);
    xpoly_zero(&gp, m);
    xpoly_zero(&Fp, m);
    xpoly_zero(&Gp, m);
    xpoly_zero(&gx_neg, n);
    xpoly_zero(&fx_neg, n);
    xpoly_zero(&Fp_x2, n);
    xpoly_zero(&Gp_x2, n);
    xpoly_zero(&Frec, n);
    xpoly_zero(&Grec, n);

    xpoly_field_norm(&fp, f, n);
    xpoly_field_norm(&gp, g, n);

    /* recursive solve on size n/2 */
    if (!NTRUSolve_exact(&fp, &gp, &Fp, &Gp, m, q)) {
        xpoly_clear(&fp);
        xpoly_clear(&gp);
        xpoly_clear(&Fp);
        xpoly_clear(&Gp);
        xpoly_clear(&gx_neg);
        xpoly_clear(&fx_neg);
        xpoly_clear(&Fp_x2);
        xpoly_clear(&Gp_x2);
        xpoly_clear(&Frec);
        xpoly_clear(&Grec);
        return 0;
    }

    /* reconstruct:
       F = F'(x^2) * g(-x)
       G = G'(x^2) * f(-x)
    */
    xpoly_substitute_neg_x(&gx_neg, g, n);
    xpoly_substitute_neg_x(&fx_neg, f, n);

    xpoly_substitute_x2(&Fp_x2, &Fp, n);
    xpoly_substitute_x2(&Gp_x2, &Gp, n);

    xpoly_mul_mod_phi_n(&Frec, &Fp_x2, &gx_neg, n);
    xpoly_mul_mod_phi_n(&Grec, &Gp_x2, &fx_neg, n);

    xpoly_set(F, &Frec, n);
    xpoly_set(G, &Grec, n);

    xpoly_clear(&fp);
    xpoly_clear(&gp);
    xpoly_clear(&Fp);
    xpoly_clear(&Gp);
    xpoly_clear(&gx_neg);
    xpoly_clear(&fx_neg);
    xpoly_clear(&Fp_x2);
    xpoly_clear(&Gp_x2);
    xpoly_clear(&Frec);
    xpoly_clear(&Grec);

    return 1;
}

int xpoly_equals_monomial_q(const xpoly *a, int n, int q) {
    if (n <= 0 || n > N) {
        fprintf(stderr, "xpoly_equals_monomial_q: invalid n=%d\n", n);
        exit(1);
    }

    if (mpz_cmp_si(a->coeffs[0], q) != 0) {
        return 0;
    }

    for (int i = 1; i < n; i++) {
        if (mpz_sgn(a->coeffs[i]) != 0) {
            return 0;
        }
    }

    return 1;
}

int xpoly_verify_ntru(const xpoly *f, const xpoly *g, const xpoly *F, const xpoly *G, int n, int q) {
    if (n <= 0 || n > N) {
        fprintf(stderr, "xpoly_verify_ntru: invalid n=%d\n", n);
        exit(1);
    }

    xpoly t1, t2, lhs;
    xpoly_init(&t1);
    xpoly_init(&t2);
    xpoly_init(&lhs);

    xpoly_zero(&t1, n);
    xpoly_zero(&t2, n);
    xpoly_zero(&lhs, n);

    xpoly_mul_mod_phi_n(&t1, f, G, n);

    xpoly_mul_mod_phi_n(&t2, g, F, n);

    xpoly_sub_n(&lhs, &t1, &t2, n);

    int ok = xpoly_equals_monomial_q(&lhs, n, q);

    if (!ok) {
        fprintf(stderr, "xpoly_verify_ntru: verification failed\n");
        xpoly_print(f, n, "f");
        xpoly_print(g, n, "g");
        xpoly_print(F, n, "F");
        xpoly_print(G, n, "G");
        xpoly_print(&t1, n, "f*G");
        xpoly_print(&t2, n, "g*F");
        xpoly_print(&lhs, n, "fG-gF");
    }

    xpoly_clear(&t1);
    xpoly_clear(&t2);
    xpoly_clear(&lhs);

    return ok;
}

void test_NTRUSolve_exact_small(void) {
    const int q = 12289;

    xpoly f, g, F, G;
    xpoly_init(&f);
    xpoly_init(&g);
    xpoly_init(&F);
    xpoly_init(&G);

    xpoly_zero(&f, 1);
    xpoly_zero(&g, 1);
    xpoly_zero(&F, 1);
    xpoly_zero(&G, 1);

    mpz_set_si(f.coeffs[0], 1);
    mpz_set_si(g.coeffs[0], 1);

    if (!NTRUSolve_exact(&f, &g, &F, &G, 1, q)) {
        printf("NTRUSolve_exact failed at n=1\n");
    } else if (!xpoly_verify_ntru(&f, &g, &F, &G, 1, q)) {
        printf("Verification failed at n=1\n");
    } else {
        printf("PASS n=1\n");
    }
    xpoly_zero(&f, 2);
    xpoly_zero(&g, 2);
    xpoly_zero(&F, 2);
    xpoly_zero(&G, 2);

    /* f = 1 + x, g = 1 */
    mpz_set_si(f.coeffs[0], 1);
    mpz_set_si(f.coeffs[1], 1);
    mpz_set_si(g.coeffs[0], 1);
    mpz_set_si(g.coeffs[1], 0);

    if (!NTRUSolve_exact(&f, &g, &F, &G, 2, q)) {
        printf("NTRUSolve_exact failed at n=2\n");
    } else if (!xpoly_verify_ntru(&f, &g, &F, &G, 2, q)) {
        printf("Verification failed at n=2\n");
    } else {
        printf("PASS n=2\n");
    }
    xpoly_zero(&f, 4);
    xpoly_zero(&g, 4);
    xpoly_zero(&F, 4);
    xpoly_zero(&G, 4);

    mpz_set_si(f.coeffs[0], 1);
    mpz_set_si(f.coeffs[1], 1);
    mpz_set_si(f.coeffs[2], 0);
    mpz_set_si(f.coeffs[3], 0);

    mpz_set_si(g.coeffs[0], 1);
    mpz_set_si(g.coeffs[1], 0);
    mpz_set_si(g.coeffs[2], 1);
    mpz_set_si(g.coeffs[3], 0);

    if (!NTRUSolve_exact(&f, &g, &F, &G, 4, q)) {
        printf("NTRUSolve_exact failed at n=4\n");
    } else if (!xpoly_verify_ntru(&f, &g, &F, &G, 4, q)) {
        printf("Verification failed at n=4\n");
    } else {
        printf("PASS n=4\n");
    }

    xpoly_clear(&f);
    xpoly_clear(&g);
    xpoly_clear(&F);
    xpoly_clear(&G);
}
void xpoly_scalar_mul_si(xpoly *res, const xpoly *a, long c, int n) {
    if (n <= 0 || n > N) {
        fprintf(stderr, "xpoly_scalar_mul_si: invalid n=%d\n", n);
        exit(1);
    }

    for (int i = 0; i < n; i++) {
        mpz_mul_si(res->coeffs[i], a->coeffs[i], c);
    }
}
