#include <cmath>

#include "helm_parameterization.h"

namespace PandaXDM {

    double j1 (double x)
    {
        return sin(x)/(x*x) - cos(x)/x;
    }

    double HelmParameterization::formFactor2 (double e)
    {
        c = 1.23 * pow(A, 1./3) - 0.60; // fm
        double rn = sqrt(c*c + 7./3*M_PI*M_PI*a*a - 5 * s * s);
        double q = 6.92e-3 * sqrt(A*e);
        double qrn = q*rn;
        double f = 3.* j1(qrn)/qrn * exp(-q*q*s*s/2.);
        return f*f;
    }
}
