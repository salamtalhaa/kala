// this program is used to demonstrate the functions in DMMath

#include <iostream>
#include <vector>

#include <TF1.h>
#include <TCanvas.h>

#include "idm_math.h"

using std::cout;
using std::endl;

int main () {

    using namespace PandaXDM;
    IDMMath idmm;

    idmm.set_darkmatter_mass(1000);
    idmm.set_cross_section(1e-40);
    idmm.set_target_mass(131.293);
    idmm.set_mass_splitting(100);

    std::vector<double> ev = {1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    for (auto & e: ev) {
        cout << "minimal velocity (" << e << " keV): " << idmm.v_min(e) << endl;
    }

    std::vector<double> masses = {1000, 10000};
    std::vector<double> dms = {0, 100, 200, 300, 400};

    TF1 * f_drde = new TF1("f_drde", &idmm, &IDMMath::dRdE_r_p, 0, 600, 3, "IDMMath", "dRdE_r_p");
    TCanvas c1("c1", "c1");
    c1.SetLogy(1);
    for (auto &m: masses) {
        f_drde->SetParameter(0, m);
        f_drde->SetParameter(1, 1e-40);
        for (auto &delta : dms) {
            f_drde->SetParameter(2, delta);
            f_drde->Draw();
            f_drde->SetTitle(Form("Inelastic WIMP dR/dE mass %.2f GeV #delta %.0f keV", m, delta));
            if ((&m == &masses.front()) && (&delta == &dms.front())) {
                c1.Print("idm_drde_xenon.pdf(");
            } else if ((&m == &masses.back()) && (&delta == &dms.back())) {
                c1.Print("idm_drde_xenon.pdf)");
            } else {
                c1.Print("idm_drde_xenon.pdf");
            }
        }
    }
}
