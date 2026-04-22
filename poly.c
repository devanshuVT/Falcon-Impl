#include "poly.h"
#include <math.h>
#include "kiss_fft.h"
#include <stdint.h>
#include <math.h>

static const int Q = 12289;
static const int OMEGA = 49;

int mod_q(int64_t x) {
    int r = (int)(x % Q);
    if (r < 0) r += Q;
    return r;
}

int mod_add(int a, int b) {
    return mod_q((int64_t)a + b);
}

int mod_mul(int a, int b) {
    return mod_q((int64_t)a * b);
}

int mod_pow(int a, int e) {
    int result = 1;
    int base = mod_q(a);

    while (e > 0) {
        if (e & 1) {
            result = mod_mul(result, base);
        }
        base = mod_mul(base, base);
        e >>= 1;
    }
    return result;
}

poly poly_zero(void) {
    poly p;
    for (int i = 0; i < N; i++) {
        p.coeffs[i] = 0.0;
    }
    return p;
}

int poly_is_zero(poly p) {
    for (int i = 0; i < N; i++) {
        if (p.coeffs[i] != 0.0) {
            return 0;
        }
    }
    return 1;
}

poly poly_add(poly a, poly b) {
    poly out;
    for (int i = 0; i < N; i++) {
        out.coeffs[i] = a.coeffs[i] + b.coeffs[i];
    }
    return out;
}

poly poly_sub(poly a, poly b) {
    poly out;
    for (int i = 0; i < N; i++) {
        out.coeffs[i] = a.coeffs[i] - b.coeffs[i];
    }
    return out;
}

poly poly_scalar_mul(poly a, double c) {
    poly out;
    for (int i = 0; i < N; i++) {
        out.coeffs[i] = a.coeffs[i] * c;
    }
    return out;
}

poly poly_adjoint(poly a) {
    poly out = poly_zero();
    out.coeffs[0] = a.coeffs[0];
    for (int i = 1; i < N; i++) {
        out.coeffs[N - i] = -a.coeffs[i];
    }
    return out;
}

poly poly_mul_mod_phi(poly a, poly b) {
    poly out = poly_zero();

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int deg = i + j;
            double term = a.coeffs[i] * b.coeffs[j];

            if (deg < N) {
                out.coeffs[deg] += term;
            } else {
                out.coeffs[deg - N] -= term;
            }
        }
    }

    return out;
}

poly poly_round_coeffs(poly a) {
    poly out;
    for (int i = 0; i < N; i++) {
        out.coeffs[i] = round(a.coeffs[i]);
    }
    return out;
}

double poly_inner(poly a, poly b) {
    double sum = 0.0;
    for (int i = 0; i < N; i++) {
        sum += a.coeffs[i] * b.coeffs[i];
    }
    return sum;
}

double poly_norm(poly a) {
    return sqrt(poly_inner(a, a));
}

poly_fft fft(poly a) {
    poly_fft out;
    kiss_fft_cpx in[N];
    kiss_fft_cpx out_cpx[N];
    kiss_fft_cfg cfg = kiss_fft_alloc(N, 0, NULL, NULL);
    for (int i = 0; i < N; i++) {
        in[i].r = a.coeffs[i];
        in[i].i = 0.0;
    }
    kiss_fft(cfg, in, out_cpx);
    for (int i = 0; i < N; i++) {
        out.real[i] = out_cpx[i].r;
        out.imag[i] = out_cpx[i].i;
    }
    free(cfg);
    return out;
}

poly ifft(poly_fft a) {
    poly out = poly_zero();
    kiss_fft_cpx in[N];
    kiss_fft_cpx out_cpx[N];
    kiss_fft_cfg cfg = kiss_fft_alloc(N, 1, NULL, NULL);
    for (int i = 0; i < N; i++) {
        in[i].r = a.real[i];
        in[i].i = a.imag[i];
    }
    kiss_fft(cfg, in, out_cpx);
    for (int i = 0; i < N; i++) {
        out.coeffs[i] = out_cpx[i].r / N;
    }
    free(cfg);
    return out;
}

poly poly_div_fft(poly num, poly den) {
    return poly_div_fft_n(num, den, N);
}

poly compute_k(poly f, poly g, poly F, poly G, int n)
{
    poly f_star = poly_adjoint_n(f, n);
    poly g_star = poly_adjoint_n(g, n);

    poly num1 = poly_mul_mod_phi_n(F, f_star, n);
    poly num2 = poly_mul_mod_phi_n(G, g_star, n);
    poly num = poly_zero();

    poly den1 = poly_mul_mod_phi_n(f, f_star, n);
    poly den2 = poly_mul_mod_phi_n(g, g_star, n);
    poly den = poly_zero();

    for (int i = 0; i < n; i++) {
        num.coeffs[i] = num1.coeffs[i] + num2.coeffs[i];
        den.coeffs[i] = den1.coeffs[i] + den2.coeffs[i];
    }

    poly k = poly_div_fft_n(num, den, n);
    return poly_round_coeffs_n(k, n);
}

double poly_pair_norm(poly a, poly b) {
    double na = poly_norm(a);
    double nb = poly_norm(b);
    return sqrt(na * na + nb * nb);
}

double compute_gamma(poly f, poly g, int q) {
    poly minus_f = poly_scalar_mul(f, -1.0);
    double term1 = poly_pair_norm(g, minus_f);

    poly f_star = poly_adjoint(f);
    poly g_star = poly_adjoint(g);

    poly ff_star = poly_mul_mod_phi(f, f_star);
    poly gg_star = poly_mul_mod_phi(g, g_star);
    poly den = poly_add(ff_star, gg_star);

    poly qf_star = poly_scalar_mul(f_star, (double)q);
    poly qg_star = poly_scalar_mul(g_star, (double)q);

    poly part1 = poly_div_fft(qf_star, den);
    poly recon = poly_mul_mod_phi(part1, den);
    poly part2 = poly_div_fft(qg_star, den);

    double term2 = poly_pair_norm(part1, part2);
    return (term1 > term2) ? term1 : term2;
}

static int poly_eval_mod_q(poly f, int x) {
    int acc = 0;

    for (int i = N - 1; i >= 0; i--) {
        int coeff = mod_q((int64_t)llround(f.coeffs[i]));
        acc = mod_add(mod_mul(acc, x), coeff);
    }

    return acc;
}

int ntt_has_zero_coefficient(poly f) {
    for (int j = 0; j < N; j++) {
        int root = mod_pow(OMEGA, 2 * j + 1);   /* roots of x^N + 1 */
        int value = poly_eval_mod_q(f, root);

        if (value == 0) {
            return 1;
        }
    }
    return 0;
}
poly poly_mul_mod_phi_n(poly a, poly b, int n) {
    poly out = poly_zero();

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int deg = i + j;
            double term = a.coeffs[i] * b.coeffs[j];

            if (deg < n) {
                out.coeffs[deg] += term;
            } else {
                out.coeffs[deg - n] -= term;
            }
        }
    }

    return out;
}
poly poly_field_norm(poly a, int n) {
    poly even = poly_zero();
    poly odd = poly_zero();
    poly odd_sq;
    poly even_sq;
    poly out = poly_zero();
    int half = n / 2;
    for (int i = 0; i < half; i++) {
        even.coeffs[i] = a.coeffs[2 * i];
        odd.coeffs[i] = a.coeffs[2 * i + 1];
    }
    even_sq = poly_mul_mod_phi_n(even, even, half);
    odd_sq = poly_mul_mod_phi_n(odd, odd, half);
    out.coeffs[0] = even_sq.coeffs[0] + odd_sq.coeffs[half - 1];
    for (int i = 1; i < half; i++) {
        out.coeffs[i] = even_sq.coeffs[i] - odd_sq.coeffs[i - 1];
    }
    return out;
}

poly poly_substitute_x2(poly a, int n) {
    poly out = poly_zero();
    int half = n / 2;
    for (int i = 0; i < half; i++) {
        out.coeffs[2 * i] = a.coeffs[i];
    }
    return out;
}

poly poly_substitute_neg_x(poly a, int n) {
    poly out = poly_zero();

    for (int i = 0; i < n; i++) {
        if (i % 2 == 0) {
            out.coeffs[i] = a.coeffs[i];
        } else {
            out.coeffs[i] = -a.coeffs[i];
        }
    }
    return out;
}

poly poly_adjoint_n(poly a, int n) {
    poly out = poly_zero();
    out.coeffs[0] = a.coeffs[0];

    for (int i = 1; i < n; i++) {
        out.coeffs[n - i] = -a.coeffs[i];
    }

    return out;
}

poly poly_sub_n(poly a, poly b, int n) {
    poly out = poly_zero();

    for (int i = 0; i < n; i++) {
        out.coeffs[i] = a.coeffs[i] - b.coeffs[i];
    }

    return out;
}

poly poly_round_coeffs_n(poly a, int n) {
    poly out = poly_zero();

    for (int i = 0; i < n; i++) {
        out.coeffs[i] = round(a.coeffs[i]);
    }

    return out;
}

int poly_is_zero_n(poly p, int n) {
    for (int i = 0; i < n; i++) {
        if (p.coeffs[i] != 0.0) {
            return 0;
        }
    }
    return 1;
}
poly poly_div_fft_n(poly num, poly den, int n) {
    poly out = poly_zero();

    kiss_fft_cpx in_num[n];
    kiss_fft_cpx in_den[n];
    kiss_fft_cpx out_num[n];
    kiss_fft_cpx out_den[n];
    kiss_fft_cpx out_quo[n];
    kiss_fft_cpx time_quo[n];

    kiss_fft_cfg cfg_fwd = kiss_fft_alloc(n, 0, NULL, NULL);
    kiss_fft_cfg cfg_inv = kiss_fft_alloc(n, 1, NULL, NULL);

    for (int i = 0; i < n; i++) {
        double angle = M_PI * (double)i / (double)n;
        double c = cos(angle);
        double s = sin(angle);

        in_num[i].r = num.coeffs[i] * c;
        in_num[i].i = num.coeffs[i] * s;

        in_den[i].r = den.coeffs[i] * c;
        in_den[i].i = den.coeffs[i] * s;
    }

    kiss_fft(cfg_fwd, in_num, out_num);
    kiss_fft(cfg_fwd, in_den, out_den);

    for (int i = 0; i < n; i++) {
        double ar = out_num[i].r;
        double ai = out_num[i].i;
        double br = out_den[i].r;
        double bi = out_den[i].i;

        double denom_mag = br * br + bi * bi;

        if (denom_mag == 0.0) {
            out_quo[i].r = 0.0;
            out_quo[i].i = 0.0;
        } else {
            out_quo[i].r = (ar * br + ai * bi) / denom_mag;
            out_quo[i].i = (ai * br - ar * bi) / denom_mag;
        }
    }

    kiss_fft(cfg_inv, out_quo, time_quo);

    for (int i = 0; i < n; i++) {
        double angle = M_PI * (double)i / (double)n;
        double c = cos(angle);
        double s = sin(angle);

        double rr = time_quo[i].r / (double)n;
        double ii = time_quo[i].i / (double)n;
        out.coeffs[i] = rr * c + ii * s;
    }

    free(cfg_fwd);
    free(cfg_inv);

    return out;
}