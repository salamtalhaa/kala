#include "quanta_nest.h"

#include <TMath.h>
#include <iostream>

QuantaNest::QuantaNest (){
    // set default values
    energy = 5;                   // keV
    field = 1.03;                 // V/cm F_0
    type = 0;                     // 0: NR; 1: ER.
    n_electron = -1;
    n_photon = -1;
    light_pE = -1;
    charge_pE = -1;
    lindhard_factor = -1;
    epsilon = 11.5 * energy * TMath::Power(54, -7./3.);
    gamma = -1;
    kappa = -1;
    density = 2.9;
    resolution = (0.12724-0.032152*density-0.0013492*TMath::Power(density,2.))*1.5;
    zeta = 0.01385 * TMath::Power(field, -0.062);
    n_quanta = -1;
    ratio = -1;
    light_quenching = 1./(1 + 3.3 * TMath::Power(epsilon, 1.14));
    dpe_fraction = 0.215;
    pmt_resolution = 0.6;
    pde = 0.1;
    eee = 0.5;
    gas_gain = 20;
    recomb_fluctuation = false;
    enable_electron_life_time = false;
    CalculateF0F1F2();
    tr.SetSeed(12345);
}

void QuantaNest::calculate()
{
    if (type==0) {
        get_lindhard_factor();
    }
    CalculateNQuanta();
    CalculateNexNi();
    CalculateRecombinationRate();
    CalculatePhotonElectron();
    CalculateCharge();
}

int QuantaNest::get_num_electron()
{
    if (n_electron<0) {
        calculate();
    }
    return n_electron;
}

int QuantaNest::get_num_photon()
{
    if (n_photon<0) {
        calculate();
    }
    return n_photon;
}

double QuantaNest::get_light_in_pe()
{
    if (light_pE<0) {
        calculate();
    }
    return light_pE;
}

double QuantaNest::get_charge_in_pe()
{
    if (charge_pE<0) {
        calculate();
    }
    return charge_pE;
}

void QuantaNest::set_energy (double e)
{
    energy = e;
    epsilon = 11.5 * energy * TMath::Power(54, -7./3.);
    light_quenching = 1./(1 + 3.3 * TMath::Power(epsilon, 1.14));
}

void QuantaNest::set_field (double f)
{
    field = f;
    zeta = 0.01385 * TMath::Power(field, -0.062);
}

void QuantaNest::set_type (int t)
{
    type = t;
}

void QuantaNest::set_density (double d)
{
    density = d;
    resolution = (0.12724-0.032152*density-0.0013492*TMath::Power(density,2.))*1.5;
}

double QuantaNest::get_lindhard_factor ()
{
    gamma = 3 * TMath::Power(epsilon, 0.15) + 0.7 * TMath::Power(epsilon, 0.6) + epsilon;
    kappa = 0.1394;
    lindhard_factor = kappa * gamma / (1 + kappa * gamma);
    return lindhard_factor;
}

int QuantaNest::get_num_quanta ()
{
    return n_quanta;
}

int QuantaNest::get_num_ionization ()
{
    return n_ionization;
}

int QuantaNest::get_num_excitation ()
{
    return n_excitation;
}

double QuantaNest::get_ratio ()
{
    return ratio;
}

double QuantaNest::get_recombination_rate ()
{
    return recombination_rate;
}

double QuantaNest::get_recombination_rate_t ()
{
    return recombination_rate_t;
}

void QuantaNest::CalculateNQuanta ()
{
    double mean_quanta = energy*1000/kW;
    double sigma = TMath::Sqrt(resolution*mean_quanta);
    n_quanta = (int) (tr.Gaus(mean_quanta, sigma));
    if (type == 0) {
        double smeared_lF = tr.Gaus(lindhard_factor, 0.25*lindhard_factor);
        n_quanta = tr.Binomial(n_quanta, smeared_lF);
    }
    if (n_quanta < 0) {
        n_quanta = 0;
    }
}

void QuantaNest::CalculateNexNi ()
{
    if (n_quanta<=0) {
        n_ionization = 0;
        n_excitation = 0;
        return;
    }
    CalculateNexNiRatio();
    n_excitation = tr.Binomial(n_quanta, ratio/(1+ratio));
    n_ionization = n_quanta - n_excitation;
}

void QuantaNest::CalculateNexNiRatio ()
{
    // only for nr.
    if (type==0) {
        ratio = 1.24 * TMath::Power(field, -0.0472)*(1 - TMath::Exp(-239*epsilon)) * 1.26;
    } else if (type==1) {
        ratio = 0.059813 + 0.031228 * density;
    }
}

void QuantaNest::CalculateRecombinationRate ()
{
    if (n_ionization == 0) {
        recombination_rate = 0;
    } else {
        if (type == 0) {
            recombination_rate = 1 - TMath::Log(1+n_ionization*zeta)/(n_ionization*zeta);
        } else {
            CalculateErRecombinationRate();
        }
        recombination_rate_t = recombination_rate;
        if (recomb_fluctuation) {
            while (1) {
                double fac = 0.1;
                if (type == 1) {
                    fac = 0.127;
                }
                recombination_rate = tr.Gaus(recombination_rate_t, fac * TMath::Sqrt(recombination_rate_t));
                if (recombination_rate<1) {
//          std::cout << recombination_rate << std::endl;
                    break;
                }
            }
        }
        if (recombination_rate<0) {
            recombination_rate = 0;
        }
        //    std::cout << recombination_rate << std::endl;
    }
}

void QuantaNest::CalculateErRecombinationRate ()
{
    double tib_f = 0.6347 * TMath::Exp(-0.00014 * field);
    double tib_e = -0.373 * TMath::Exp(-(field*0.001)/tib_f) + 1.5;
    double tib_curl_a = 10 * TMath::Power(field, -0.04) * TMath::Exp(18/field);
    double tib_curl_z = 1 - TMath::Power(field, 0.2147) + 3;
    double tomas_imel = tib_f * TMath::Power(energy, -tib_e) * (1 - TMath::Exp(-TMath::Power((energy - tib_curl_z)/tib_curl_a, 0.188 * TMath::Power(field, 1./3)))) * TMath::Power(density/2.888, 0.3);
    recombination_rate = 1 - TMath::Log(1+tomas_imel/4.*n_ionization)/(tomas_imel/4.*n_ionization);
}

void QuantaNest::CalculatePhotonElectron ()
{
    n_photon = n_excitation + tr.Binomial(n_ionization, recombination_rate);
    n_electron = n_quanta - n_photon;
    if (type==0) {
        n_photon = tr.Binomial(n_photon, light_quenching);
    }
}

void QuantaNest::CalculateCharge ()
{
    // calculate s1
    // 1. detected photon
    int i_photon = tr.Binomial(n_photon, 1-f0);
    if (i_photon>0) {
        // 2. generated pe
        int i_photon2 = tr.Binomial(i_photon, dpe_fraction);
        i_photon += i_photon2;
        // 3. smearing to non-negative value
        light_pE = -1;
        while (light_pE <= 0) {
            light_pE = tr.Gaus(i_photon, pmt_resolution * TMath::Sqrt(i_photon));
        }
    } else {
        light_pE = 0;
    }

    // calculate s2
    // 1. extracted electrons
    int i_electron;
    if (enable_electron_life_time) {
        double survive_rate = TMath::Exp(-drift_time/electron_life_time);
        int ns_electron = tr.Binomial(n_electron, survive_rate);
        i_electron = tr.Binomial(ns_electron, eee);
    } else {
        i_electron = tr.Binomial(n_electron, eee);
    }
    if (i_electron > 0) {
        double g_var = -1;
        while (g_var <= 0) {
            g_var = tr.Gaus(1, 0.165);
        }
        charge_pE = -1;
        while (charge_pE <=0 ) {
            charge_pE = tr.Gaus(i_electron * gas_gain * g_var, pmt_resolution * TMath::Sqrt(i_electron * gas_gain * g_var));
        }
    } else {
        charge_pE = 0;
    }
}

void QuantaNest::set_pmt_resolution (double r)
{
    pmt_resolution = r;
}

void QuantaNest::set_dpe_fraction (double f)
{
    dpe_fraction = f;
    CalculateF0F1F2();
}

void QuantaNest::set_pde (double p)
{
    pde = p;
    CalculateF0F1F2();
}

void QuantaNest::set_eee (double e)
{
    eee = e;
}

void QuantaNest::set_gas_gain (double g)
{
    gas_gain = g;
}

void QuantaNest::set_recomb_fluctuation (bool t)
{
    recomb_fluctuation = t;
}

void QuantaNest::set_electron_life_time (double t)
{
    if (t>0) {
        electron_life_time = t;
        enable_electron_life_time = true;
    }
}

void QuantaNest::set_drift_time (double t)
{
    if (t>0) {
        drift_time = t;
    }
}

void QuantaNest::CalculateF0F1F2 ()
{
    f0 = 1 - pde/(1 + dpe_fraction);
}
