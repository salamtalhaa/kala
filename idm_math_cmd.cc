#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdio>

#include <TF1.h>

#include "idm_math.h"

using std::cin;
using std::cout;
using std::endl;

int main () {

    using namespace PandaXDM;
    IDMMath idmm;

    idmm.set_darkmatter_mass(1000);
    idmm.set_cross_section(1e-40);
    idmm.set_target_mass(131.293);
    idmm.set_mass_splitting(100);

    std::string s;
    cout << "Please input mass splitting and cross section, type 'q' to exit." << endl;
    cout << "Example: 100 1e-34" << endl;
    double dm, xsec;

    TF1 * f_drde = new TF1("f_drde", &idmm, &IDMMath::dRdE_r_p, 0, 300, 3, "IDMMath", "dRdE_r_p");
    f_drde->SetNpx(10000);

    while (1) {
        std::getline (cin, s);
        auto n = s.find("q");
        if (n != std::string::npos) {
            break;
        }
        {
            std::stringstream ss(s);
            ss >> dm >> xsec;
            f_drde->SetParameter(0, 1000);
            f_drde->SetParameter(1, xsec);
            f_drde->SetParameter(2, dm);

            cout << "rate (cts/day/kg) = " << f_drde->Integral(1, 100) << "." << endl;
            cout << "expected in pandax ii 98.7 days: " << f_drde->Integral(1, 100) * 3.3e4 << "." << endl;
            cout << "----------" << endl;
        }
    }
    cout << "bye bye!" << endl;
}
