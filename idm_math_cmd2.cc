#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdio>

#include <TF1.h>

#include "idm_math.h"

using std::cin;
using std::cerr;
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
    cerr << "Please input mass splitting, type 'q' to exit." << endl;
    cerr << "Example: 100" << endl;
    double dm;

    TF1 * f_drde = new TF1("f_drde", &idmm, &IDMMath::dRdE_r_p, 0, 300, 3, "IDMMath", "dRdE_r_p");
    f_drde->SetNpx(10000);

    while (1) {
        std::getline (cin, s);
        if (cin.eof())
            return 0;
        auto n = s.find("q");
        if (n != std::string::npos) {
            cerr << "eof" << endl;
            break;
        }
        {
            std::stringstream ss(s);
            ss >> dm;
            f_drde->SetParameter(0, 1000);
            f_drde->SetParameter(1, 1e-45);
            f_drde->SetParameter(2, dm);
            double n_exp = f_drde->Integral(1, 100) * 3.3e4;
            double sigma = 6/n_exp*1e-45;
            for (auto i=0; i<10; ++i) {
                cout << dm << " " << sigma << endl;
                sigma /= 3;
            }
        }
    }
    cerr << "bye bye!" << endl;
}
