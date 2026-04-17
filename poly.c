#include "poly.h"
#include <math.h>

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