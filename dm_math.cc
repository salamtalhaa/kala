#include <cmath>
#include <iostream>

#include <TF2.h>
#include "dm_math.h"

using namespace std;

namespace PandaXDM {

    double f_v(double *x, double *par)
    {
        double v = x[0];
        double cos_theta = x[1];
        double v_0 = par[0];
        double v_min = par[1];
        double v_E = par[2];
        double v_esc = par[3];
        double k1 = par[4];
        double v_dm = sqrt(v*v + v_E*v_E + 2*v*v_E*cos_theta);
        if (v < v_min || v_dm>v_esc) {
            return 0;
        }
        double ratio = v_dm/v_0;
        return 2*M_PI*exp(-ratio*ratio)*v/k1;
    }

    DMMath::DMMath()
    {
        helm.A = A;
        double x[37]={1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,25,30,40,50,60,70,80,90,100};
        double y[37]={0.00165961,0.00251345,0.00388933,0.00546004,0.00729435,0.00972375,0.013181,0.0163768,0.0209521,0.0252034,0.105437,0.238027,0.378688,0.510054,0.619531,0.703173,0.761702,0.809053,0.834818,0.856425,0.870151,0.881208,0.886329,0.888934,0.889893,0.892343,0.889077,0.883818,0.801981,0.63864,0.285158,0.114345,0.0494859,0.0221751,0.0129725,0.00827964,0.00438017};
        nr_eff_curve_ = new TGraph(37, x, y);
    }

    double DMMath::k0()
    {
        return pow(M_PI*v_0*v_0, 1.5);
    }

    double DMMath::k1()
    {
        return k0() * (erf(v_esc/v_0) - 2/sqrt(M_PI)*v_esc/v_0*exp(-v_esc*v_esc/v_0/v_0));
    }

    double DMMath::get_form_factor_square (double e)
    {
        return helm.formFactor2(e);
    }

    double DMMath::v_min (double e)
    {
        double m_target = get_target_mass_GeV();
        double m_reduce = target_reduce_mass();
        // cout << "dark matter mass: " << m_dm << " GeV/c2" << endl;
        // cout << "target mass: " << m_target << " GeV/c2" << endl;
        // cout << "reduced mass: " << m_reduce << " GeV/c2" << endl;
        double v = sqrt(0.5 * m_target * e * 1e-6)/m_reduce * v_c;
        return v*1e-3;
    }

    double DMMath::target_reduce_mass ()
    {
        double m_target = get_target_mass_GeV ();
        return (m_dm * m_target) / (m_dm + m_target);
    }

    double DMMath::proton_reduce_mass ()
    {
        return (m_dm * m_proton) / (m_dm + m_proton);
    }

    double DMMath::kg_to_GeV(double m)
    {
        return m/1.783e-27;
    }

    double DMMath::GeV_to_kg(double m)
    {
        return m*1.783e-27;
    }

    double DMMath::get_target_mass_GeV ()
    {
        return kg_to_GeV(A/1000/N_0);
    }

    double DMMath::v_dis_integration(double e)
    {
        TF2 f("f", &f_v, 0, 1000, -1, 1, 5);
        f.SetParameters(v_0, v_min(e), v_E, v_esc, k1());
        double res = f.Integral(0, 1000, -1, 1);
        return res;
    }

    double DMMath::v_dis_integration_analytical(double e)
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

    double DMMath::dRdE (double e)
    {
        double mu_p = proton_reduce_mass();
        double res = 0.5 * rho_dm * sigma_0 * A * A * helm.formFactor2(e)
            * v_dis_integration_analytical(e)
            / m_dm / mu_p / mu_p;
        res = res / 1e7 / GeV_to_kg(1) * 86400 * v_c * v_c;
        return res;
    }

    double DMMath::getFormFactor2(double e)
    {
        return helm.formFactor2(e);
    }

    double DMMath::dRdE_r (double * e, double *) {
        return dRdE(*e);
    }

    double DMMath::dRdE_r_p (double * e, double *p) {
        set_darkmatter_mass(p[0]);
        set_cross_section(p[1]);
        return dRdE(*e);
    }

    double DMMath::getFormFactor2_r(double * e, double *) {
        return getFormFactor2(*e);
    }

    double DMMath::v_dis_integration_r(double *e, double *) {
        return v_dis_integration(*e);
    }

    double DMMath::v_dis_integration_analytical_r(double *e, double *) {
        return v_dis_integration_analytical(*e);
    }
}
