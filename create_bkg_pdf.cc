// this program is used to create the histograms for different background source

#include "quanta_nest.h"
#include "nest_er_model.h"

#include <TFile.h>
#include <TH2F.h>
#include <TRandom3.h>
#include <TMath.h>
#include <TGraph.h>

#include <iostream>
#include <sstream>
#include <string>
#include <memory>
#include <vector>

using std::cout;
using std::endl;
using std::cerr;
using std::stringstream;
using std::string;
using std::vector;

template<typename T>
void parse_parameter(T& p, const char * c) {
    std::stringstream ss(c);
    ss >> p;
}

void ER_sim2 (double e_comb, double & s1, double & logs2s1);
void ER_sim (double & s1, double & logs2s1, double &eff);

struct Params {
    double pde;
    double eee;
    double seg;
    double e_life;
    double s2_thres;
    double s2_w;
    double s1_thres;
    double s1_w;
};

namespace util {
    TRandom3 tr3;
    Params params;
}

using namespace util;

int main (int argc, char * argv[])
{
    if (argc!=11) {
        cerr << "Usage: " << argv [0] << " pde eee seg e_life s2_thres s2_w s1_thres s1_w num_events output_file" << endl;
        return 1;
    }
    int num_events;
    parse_parameter(params.pde, argv[1]);
    parse_parameter(params.eee, argv[2]);
    parse_parameter(params.seg, argv[3]);
    parse_parameter(params.e_life, argv[4]);
    parse_parameter(params.s2_thres, argv[5]);
    parse_parameter(params.s2_w, argv[6]);
    parse_parameter(params.s1_thres, argv[7]);
    parse_parameter(params.s1_w, argv[8]);
    parse_parameter(num_events, argv[9]);
    string output_file_name{argv[10]};

    // quanta nest
    QuantaNest nest;
    nest.set_field(400);
    nest.set_gas_gain(params.seg);
    nest.set_pde(params.pde);
    nest.set_eee(params.eee);
    nest.set_recomb_fluctuation(false);

    // create the output file
    TFile fout(output_file_name.c_str(), "RECREATE");

    vector<string> bkg_names{"xe127", "kr85", "other_ER", "accidental", "nbkg"};

    // generate the bkg histogram one by one
    vector<TH2F *> hPDF_bkg;

    for (size_t i = 0; i<bkg_names.size(); ++i) {
        hPDF_bkg.push_back(new TH2F(bkg_names[i].c_str(), "", 200, 0, 100, 200, 0, 4));
    }

    // initialize the random generator
    tr3.SetSeed(0);
    gRandom->SetSeed(0); // required for histogram random number

    // xe127
    // we need to sample from two peaks: 5.2 keV and 33.2 keV
    // the ratio is 13.1% over 83.4%
    constexpr double ratio = 13.1/(13.1 + 83.4);
    for (int i = 0; i < num_events; ++i) {
        double v = tr3.Uniform(0, 1);
        double e_comb;
        if (v < ratio)
            e_comb = tr3.Gaus(5.2, 0.5);
        else
            e_comb = tr3.Gaus(33.2, 1);
        double s1, logs2s1;
        ER_sim2(e_comb, s1, logs2s1);
        if (s1 > 0 && s1 < 100)
            hPDF_bkg[0]->Fill(s1, logs2s1);
    }

    vector<double> er_edep, er_s1;
    for (auto p: e_s1_map) {
        er_edep.push_back(p.first);
        er_s1.push_back(p.second);
    }
    TGraph er_nest_graph(er_edep.size(), er_edep.data(), er_s1.data());

    // Kr85
    for (int i = 0; i < num_events; ++i) {
        double e_comb = tr3.Uniform(0, 50);
        double s1 = er_nest_graph.Eval(e_comb);
        double logs2s1, eff;
        ER_sim(s1, logs2s1, eff);
        hPDF_bkg[1]->Fill(s1, logs2s1, eff);
    }

    // other ER
    for (int i=0; i < num_events; ++i) {
        double e_comb = tr3.Uniform(0, 50);
        double s1 = er_nest_graph.Eval(e_comb);
        double logs2s1, eff;
        ER_sim(s1, logs2s1, eff);
        hPDF_bkg[2]->Fill(s1, logs2s1, eff);
    }

    // accidental bkg
    TFile f_acc("accidental.root", "READ");
    TH2F * hh = static_cast<TH2F *>(f_acc.Get("hh2"));
    for (int i = 0; i < num_events; ++i) {
        double s1, logs2s1;
        hh->GetRandom2(s1, logs2s1);
        hPDF_bkg[3]->Fill(s1, logs2s1);
    }
    f_acc.Close();

    // neutron bkg with nest
    TFile f_neutron("nr_bkg_spectrum_total_norm.root", "READ");
    TH1F * h_Enr = static_cast<TH1F *> (f_neutron.Get("spectot"));
    nest.set_type(0);
    TH2F neutron_bkg("neutron_bkg", "", 200, 0, 100, 200, 0, 4);
    constexpr double e_cut = 1.1;
    for (int i = 0; i < num_events; ++i) {
        double e_nr = h_Enr->GetRandom();
        nest.set_energy(e_nr);
        nest.calculate();
        double s1 = nest.get_light_in_pe();
        double s2 = nest.get_charge_in_pe();
        double dt = tr3.Uniform(18, 310);
        double s2_raw = s2 * TMath::Exp(- dt / params.e_life);
        double eff_s1 = 1 / (TMath::Exp(- (s1 - params.s1_thres) / params.s1_w) + 1);
        double eff_s2 = 1 / (TMath::Exp(- (s2_raw - params.s2_thres) / params.s2_w) + 1);
        double eff = eff_s1 * eff_s2 * 0.94;
        double logs2s1 = TMath::Log10(s2/s1);
        if ( e_nr < e_cut || s2_raw < 100 || s2 > 12000
             || logs2s1 < 1.074 - 0.595 * TMath::Exp(-s1/6.38)) {
            neutron_bkg.Fill(-1, -1);
        } else {
            neutron_bkg.Fill(s1, logs2s1, eff);
        }
    }
    for (int i = 0; i < num_events; ++i) {
        double s1, logs2s1;
        neutron_bkg.GetRandom2(s1, logs2s1);
        hPDF_bkg[4]->Fill(s1, logs2s1);
    }
    f_neutron.Close();

    fout.cd();
    // write the histograms into the file
    for (size_t i = 0; i < hPDF_bkg.size(); ++i) {
        hPDF_bkg[i]->Write();
    }

    fout.Close();
}

void ER_sim2 (double e_comb, double & s1, double & logs2s1)
{
    double band_center = -0.2012 + 0.917 * TMath::Exp(-e_comb/2.34);
    double sigma = 0.116 + 0.097 * TMath::Exp(-e_comb/7);
    double log2 = tr3.Gaus(band_center, sigma);

    s1 = e_comb * params.pde / 0.0137 / (1 + TMath::Power(10, log2));
    double s2 = e_comb * params.eee * params.seg / 0.0137 * TMath::Power(10, log2) / (1 + TMath::Power(10, log2));

    logs2s1 = TMath::Log10(s2/s1);
}

void ER_sim (double & s1, double & logs2s1, double &eff)
{
    double band_center = 1.79 + 0.81 * TMath::Exp(-s1/10.65);
    constexpr double sigma = 0.142;

    logs2s1 = tr3.Gaus(band_center, sigma);
    double s2 = TMath::Power(10, logs2s1) * s1;
    double dt = tr3.Uniform(18, 310);
    double s2_raw = s2 * TMath::Exp(- dt / params.e_life);
    double eff_s1 = 1 / (TMath::Exp(-(s1 - params.s1_thres) / params.s1_w) + 1);
    double eff_s2 = 1 / (TMath::Exp(-(s2_raw - params.s2_thres) / params.s2_w) + 1);
    eff = eff_s1 * eff_s2;
    // add BDT efficiency
    eff *= 0.9716 / (TMath::Exp(-(s1 - 2.567)/0.9268) + 1) * (1 - 0.0005674 * s1);

}
