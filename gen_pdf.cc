// this program is used to generate the s1 and s2
// distribution for a given wimp mass m and cross section

#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <cstdio>

#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TRandom.h>

#include "dm_math.h"
#include "quanta_nest.h"

int main(int argc, char * argv[])
{
    if (argc!=3) {
        std::cerr << "Usage: " << argv[0] << " mass cross_section" << std::endl;
        return 1;
    }
    double mass;
    double cross_section;
    {
        std::stringstream ss(argv[1]);
        ss >> mass;
    }
    {
        std::stringstream ss(argv[2]);
        ss >> cross_section;
    }

    using namespace PandaXDM;
    DMMath dmm;
    dmm.set_target_mass(131.293);
    TF1 * f_drde = new TF1("f_drde", &dmm, &DMMath::dRdE_r_p, 0, 100, 2, "DMMath", "dRdE_r_p");
    f_drde->SetParameter(0, mass);
    f_drde->SetParameter(1, cross_section);

    gRandom->SetSeed(0);

    QuantaNest nest;
    nest.set_field(400);
    nest.set_type(0);
    nest.set_gas_gain(24.4);
    nest.set_eee(0.4604);
    nest.set_pde(0.1176);
    nest.set_recomb_fluctuation(false);

    std::string output_file_name;
    {
        std::unique_ptr<char[]> buf (new char[128]);
        sprintf(buf.get(), "s_dis_mass_%.2f_%.2e.root", mass, cross_section);
        output_file_name = buf.get();
    }

    TFile f(output_file_name.c_str(), "RECREATE");
    TTree * signal_tree = new TTree("signal", "Signal");
    double e, s1, s2;
    signal_tree->Branch("Energy", &e, "Energy/D");
    signal_tree->Branch("S1", &s1, "S1/D");
    signal_tree->Branch("S2", &s2, "S2/D");
    double sum_rate = f_drde->Integral(0.1, 80);
    double sum_weight{0};
    int cycle = 100000;
    for (auto i=0; i<cycle; ++i) {
        e = f_drde->GetRandom(0.1, 80);
        nest.set_energy(e);
        nest.calculate();
        s1 = nest.get_light_in_pe();
        s2 = nest.get_charge_in_pe();
        if (s1>0&&s2>0) {
            signal_tree->Fill();
            // following is the search window
            if ((s1>3 && s1<45)
                && (log10(s2/s1) < 1.46+0.55*exp(-s1/12.45))    // nr median
                && (log10(s2/s1) > 1.074-0.595*exp(-s1/6.38))   // nr 999
                && s2>100) {
                sum_weight += 0.94
                    / (exp(-(s2-76)/19) + 1) // s2 efficiency
                    / (exp(-(s1-3.2)/0.44) + 1) // s1 efficiency
                    * 0.9716/(exp(-(s1-2.567)/0.9268) + 1); // BDT efficiency
            }
        }
    }
    signal_tree->Write();
    f.Close();
    std::cout << "N (expected counts/day/kg) = " << sum_weight/cycle*sum_rate << std::endl;
    std::cout << "sum weight (cycle): " << sum_weight << " (" << cycle << ")" << std::endl;
    std::cout << "sum_rate: " << sum_rate << std::endl;
}
