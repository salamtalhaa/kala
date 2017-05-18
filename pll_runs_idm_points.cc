#include "max_likelihood_runs.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cstdio>
#include <vector>
#include <utility>

#include <TFile.h>
#include <TTree.h>

template<typename T>
void parse_parameter(T& p, const char * c) {
    std::stringstream ss(c);
    ss >> p;
}

struct FitResult {
    int flag;
    int status;
    double log_likelihood;
    double dm_mass;
    double delta_mass;
    double xsec;
    int is_fixed;
};

using namespace std;

using ResultsVec = std::vector<FitResult>;

namespace var {
    int flag, status, is_fixed;
    double delta_mass, delta_mass_old, xsec, q_mu;
    double dm_mass = 1000;
}

using namespace var;

void fill_results_to_tree (ResultsVec & results_vec, TTree * tree );
int main(int argc, char * argv[])
{
    if (argc!=5) {
        cerr << "Usage: " << argv[0] << " config_file points_file output_root_file output_q_file"<< endl;
        return 1;
    }

    // read points from file
    std::ifstream fin(argv[2]);
    std::vector<std::pair<double, double>> points;
    while(fin.good()) {
        double d_mass, xsec;
        fin >> d_mass >> xsec;
        if (fin.eof())
            break;
        points.emplace_back(d_mass, xsec);
    }
    fin.close();
    std::vector<double> unique_delta_masses;
    std::ofstream fout(argv[4]);
    for (const auto & point: points) {
        std::cout << "mass = " << point.first << ", xsec = " << point.second << std::endl;
    }
// //    string output_q_file_name(argv[10]);

    TFile output_root_file(argv[3], "RECREATE");
    TTree * tree = new TTree("pll_tree", "tree to store xsec scan results");
    tree->Branch("dm_mass", &dm_mass, "dm_mass/D");
    tree->Branch("delta_mass", &delta_mass_old, "delta_mass/D");
    tree->Branch("xsec", &xsec, "xsec/D");
    tree->Branch("q_mu", &q_mu, "q_mu/D");
    tree->Branch("flag", &flag, "flag/I");
    tree->Branch("status", &status, "status/I");
    tree->Branch("is_fixed", &is_fixed, "is_fixed/I");
    MaxLikelihoodRuns mlrs(argv[1]);

//     double delta_mass_delta = (delta_m_high - delta_m_low) / delta_m_steps;

    std::cout << "trying to do the fit" << std::endl;

    std::vector<FitResult> results_vec;
    // loop over all mass steps
    double global_cross_section_min, max_log_likelihood;
    for (const auto & point: points) {
        delta_mass = point.first;
        if (std::find(unique_delta_masses.begin(),
                      unique_delta_masses.end(),
                      delta_mass) == unique_delta_masses.end()) {
            // first occurrence of the delta mass
            if (results_vec.size()>0) {
                fill_results_to_tree(results_vec, tree);
                // force the best fit cross section and fit again,
                // then dump the best fit histogram to file
                mlrs.MaximumLikelihood(results_vec[0].xsec, max_log_likelihood, false);
                mlrs.dumpBestFitResultsToFile(&output_root_file);
                // clear the results vec
                results_vec.clear();
            }
            unique_delta_masses.push_back(delta_mass);
            mlrs.SetIDMProperties(dm_mass, delta_mass, 1);
            mlrs.computeAverageEfficiencies();
            global_cross_section_min = 0;
            max_log_likelihood = 0;
            // fit by floating the cross section
            std::cout << "fitting with floating cross section for delta mass " << delta_mass << " keV" << std::endl;
            status = mlrs.MaximumLikelihood(global_cross_section_min, max_log_likelihood, true);
            if (status == 0) {
                results_vec.insert(results_vec.end(), {0, status, max_log_likelihood, dm_mass, delta_mass, global_cross_section_min, 0});
            }
            // fit with fixed 0 cross section
            std::cout << "fitting with fixed cross section 0 for delta mass " << delta_mass << " keV" << std::endl;
            xsec = 0;
            status = mlrs.MaximumLikelihood(xsec, max_log_likelihood, false);
            if (status == 0) {
                results_vec.insert(results_vec.end(), {0, status, max_log_likelihood, dm_mass, delta_mass, xsec, 1});
            }
        }
        // fit with fixed non-zero cross section
        xsec = point.second;
        std::cout << "fitting delta mass " << delta_mass << " keV with cross section " << xsec << std::endl;
        status = mlrs.MaximumLikelihood(xsec, max_log_likelihood, false);
        if (status == 0) {
            results_vec.insert(results_vec.end(), {0, status, max_log_likelihood, dm_mass, delta_mass, xsec, 1});
        }
        delta_mass_old = delta_mass;
    }
    // fill the last energy point
    if (results_vec.size()>0) {
        fill_results_to_tree(results_vec, tree);
        // force best fit and fit again, dump results to file
        mlrs.MaximumLikelihood(results_vec[0].xsec, max_log_likelihood, false);
        mlrs.dumpBestFitResultsToFile(&output_root_file);
    }
    output_root_file.cd();
    tree->Write();
    output_root_file.Close();
}

void fill_results_to_tree (ResultsVec & results_vec, TTree * tree )
{
    std::sort(results_vec.begin(), results_vec.end(),
              [](FitResult & a, FitResult &b) {
                  return a.log_likelihood < b.log_likelihood;
              });
    for (auto & res: results_vec) {
        xsec = res.xsec;
        if (&res == &results_vec[0]) {
            flag = 0;
            q_mu = res.log_likelihood;
        } else {
            flag = 1;
            if (res.xsec <= results_vec[0].xsec) {
                q_mu = 0;
            } else {
                q_mu = res.log_likelihood - results_vec[0].log_likelihood;
            }
        }
        status = res.status;
        is_fixed = res.is_fixed;
        tree->Fill();
    }
}
