// this program is used to demonstrate the functions in DMMath

#include <iostream>
#include <vector>
#include "dm_math.h"

#include <TF1.h>
#include <TF2.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TRandom.h>
#include <TError.h>

using std::cout;
using std::endl;

int main () {
    using namespace PandaXDM;
    DMMath dmm;

    dmm.set_darkmatter_mass(50);
    dmm.set_cross_section(1e-45);
    dmm.set_target_mass(131.293);

    gErrorIgnoreLevel = kWarning;

    cout << "event rate" << endl;
    std::vector<double> ev = {1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    for (auto & e: ev) {
        cout << "minimal velocity (" << e << " keV): " << dmm.v_min(e) << endl;
        cout << "average inverse velocity (" << e << " keV): " << dmm.v_dis_integration(e) << endl;
        cout << "form factor (" << e << " keV): " << dmm.getFormFactor2(e) << endl;
        cout << e << " ==> " << dmm.dRdE(e) << endl;
        cout << " == " << endl;
    }

    // Let's plot the dR/dE distribution of xenon
    TF1 * fr = new TF1("fr", &dmm, &DMMath::dRdE_r, 0, 100, 0, "DMMath", "dRdE_r");
    TF1 * ffr = new TF1("ffr", &dmm, &DMMath::getFormFactor2_r, 0, 100, 0, "DMMath", "getFormFactor2_r");
    TF1 * fvr = new TF1("fvr", &dmm, &DMMath::v_dis_integration_analytical_r, 0, 100, 0, "DMMath", "v_dis_integration_analytical_r");
    TCanvas c1("c1");
    c1.SetLogy(1);
    fr->Draw();
    c1.Print("rate_xenon.png");
    ffr->Draw();
    c1.Print("form_xenon.png");
    fvr->Draw();
    c1.Print("v_xenon.png");

    // plot the dr/dE distribution for different mass of dark matter
    TF1 * fr2 = new TF1("fr2", &dmm, &DMMath::dRdE_r_p, 0, 100, 2, "DMMath", "dRdE_r_p");

    // m_dm = 100 GeV
    fr2->SetParameter(0, 100);
    fr2->SetParameter(1, 1e-45);
    fr2->Draw();
    cout << fr2->Integral(0.1, 60) << endl;
    c1.Print("drde_xenon_100.png");

    // m_dm = 50 GeV
    fr2->SetParameter(0, 50);
    fr2->SetParameter(1, 1e-45);
    fr2->Draw();
    cout << fr2->Integral(0.1, 60) << endl;
    c1.Print("drde_xenon_50.png");

    // plot dr/de for different WIMP mass
    // and different cross sections
    std::vector<double> masses = {3, 4, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 500, 1000, 2000, 5000, 10000};
    std::vector<double> c_sections = {1e-42, 1e-43, 1e-44, 1e-45, 1e-46};
    for (auto &m: masses) {
        fr2->SetParameter(0, m);
        for (auto &sec: c_sections) {
            fr2->SetParameter(1, sec);
            fr2->Draw();
            fr2->SetTitle(Form("dR/dE WIMP Mass %.2f GeV @%.2e cm^{2}", m, sec));
            if ((&m == &masses.front()) && (&sec == &c_sections.front())) {
                c1.Print("drde_xenon.pdf(");
            } else if ((&m == &masses.back()) && (&sec == &c_sections.back())) {
                c1.Print("drde_xenon.pdf)");
            } else {
                c1.Print("drde_xenon.pdf");
            }
        }
    }

    // example of sampling from the dr/de distribution
    TH1D * edis = new TH1D("edis", "Recoil Energy;Energy (keV); Count/keV", 60, 0, 60);
    gRandom->SetSeed(0);
    for (int i=0; i<100000; ++i) {
        edis->Fill(fr2->GetRandom(1, 60));
    }
    edis->Scale(1./100000);
    edis->Scale(fr2->Integral(1, 60));
    edis->Draw();
    c1.Print("edis.png");
}
