tast
#include <cmath>
#include <iostream>

#include "idm_math.h"

using namespace std;

namespace PandaXDM {

    IDMMath::IDMMath()
    {
        helm.A = A;
    }

    double IDMMath::k0()
    {
        return pow(M_PI*v_0*v_0, 1.5);
    }

    double IDMMath::k1()
    {
        return k0() * (erf(v_esc/v_0) - 2/sqrt(M_PI)*v_esc/v_0*exp(-v_esc*v_esc/v_0/v_0));
    }

    double IDMMath::get_form_factor_square (double e)
    {
        return helm.formFactor2(e);
    }

    double IDMMath::v_min (double e)
    {
        double m_target = get_target_mass_GeV();
        double m_reduce = target_reduce_mass();
        double v = (e * m_target /m_reduce + delta_m_) * 1e-6
            / sqrt(2 * m_target * e) * v_c;
        return v;
    }

    double IDMMath::target_reduce_mass ()
    {
        double m_target = get_target_mass_GeV ();
        return (m_dm * m_target) / (m_dm + m_target);
    }

    double IDMMath::proton_reduce_mass ()
    {
        return (m_dm * m_proton) / (m_dm + m_proton);
    }

    double IDMMath::kg_to_GeV(double m)
    {
        return m/1.783e-27;
    }

    double IDMMath::GeV_to_kg(double m)
    {
        return m*1.783e-27;
    }

    double IDMMath::get_target_mass_GeV ()
    {
        return kg_to_GeV(A/1000/N_0);
    }

    double IDMMath::v_dis_integration_analytical(double e)
    {
        double vv_min = v_min(e);
        if (vv_min > v_esc + v_E) {
            return 0;
        }
        if (vv_min > v_esc - v_E) {
            return 0.5*M_PI*v_0*v_0/v_E/k1()
                * (-2. * exp(-v_esc*v_esc/v_0/v_0)*(v_E+v_esc-vv_min)
                   + sqrt(M_PI)*v_0*(erf(v_esc/v_0) + erf((v_E-vv_min)/v_0)));
        }
        return (0.5*M_PI*v_0*v_0*(-4*exp(-v_esc*v_esc/v_0/v_0)
                                 + sqrt(M_PI)*v_0*(erf((2*v_E-v_esc)/v_0) + erf(v_esc/v_0))/v_E)
            + 0.5*pow(M_PI, 1.5)*v_0*v_0*v_0/v_E*(-erf(v_esc/v_0)
                                                  + erf((-2*v_E+v_esc)/v_0)
                                                  + erf((v_E-vv_min)/v_0)
                                                  + erf((v_E+vv_min)/v_0)))/k1();
    }

    double IDMMath::dRdE (double e)
    {
        double mu_p = proton_reduce_mass();
        double res = 0.5 * rho_dm * sigma_0 * A * A * helm.formFactor2(e)
            * v_dis_integration_analytical(e)
            / m_dm / mu_p / mu_p;
        res = res / 1e7 / GeV_to_kg(1) * 86400 * v_c * v_c;
        return res;
    }

    double IDMMath::getFormFactor2(double e)
    {
        return helm.formFactor2(e);
    }

    double IDMMath::dRdE_r (double * e, double *) {
        return dRdE(*e);
    }

    double IDMMath::dRdE_r_p (double * e, double *p) {
        set_darkmatter_mass(p[0]);
        set_cross_section(p[1]);
        set_mass_splitting(p[2]);
        return dRdE(*e);
    }

    double IDMMath::getFormFactor2_r(double * e, double *) {
        return getFormFactor2(*e);
    }

    double IDMMath::v_dis_integration_analytical_r(double *e, double *) {
        return v_dis_integration_analytical(*e);
    }
}
