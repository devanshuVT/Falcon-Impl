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