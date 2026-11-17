#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <cmath>

/* Real roots of quadratics and cubics, used by the binomial ADMM solvers. */

inline int solve_quadratic(double A, double B, double C, double &r1, double &r2) {
    const double EPS = 1e-14;
    double D = B*B - 4*A*C;
    if (std::fabs(D) < EPS) {
        r1 = r2 = -B/(2*A);
        return 1;
    } else if (D > 0) {
        double s = std::sqrt(D);
        r1 = (-B + s)/(2*A);
        r2 = (-B - s)/(2*A);
        return 2;
    }
    return 0;
};

inline void solve_cubic(double a, double b, double c, double d, double &x1, double &x2, double &x3) {
    const double EPS = 1e-14;
    if (std::fabs(d) < EPS) {
        x1 = 0;
        int n = solve_quadratic(a, b, c, x2, x3);
        if (n == 0) x2 = x3 = 0;
        return;
    }

    double A = b/a, B = c/a, C = d/a;
    double p = B - (A*A)/3, q = (2*A*A*A)/27 - (A*B)/3 + C;
    double D = (q*q)/4 + (p*p*p)/27;
    double shift = -A/3;
    if (D > EPS) {
        double s = std::sqrt(D);
        double alpha = -q/2 + s, beta = -q/2 - s;
        double u = std::cbrt(alpha), v = std::cbrt(beta);
        x1 = u + v + shift;
        x2 = x3 = x1;
    } else if (std::fabs(D) < EPS) {
        double u = std::cbrt(-q/2);
        x1 = 2*u + shift;
        x2 = x3 = -u + shift;
    } else {
        double r = 2*std::sqrt(-p/3);
        double phi = std::acos((-q/2)/std::sqrt(-(p*p*p)/27))/3;
        x1 = r*std::cos(phi) + shift;
        x2 = r*std::cos(phi + 2*M_PI/3) + shift;
        x3 = r*std::cos(phi + 4*M_PI/3) + shift;
    }
}

#endif
