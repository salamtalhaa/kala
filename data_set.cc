#include "data_set.h"

#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <regex>
#include <cmath>
#include <memory>

#include <TKey.h>
#include <TROOT.h>
#include <TClass.h>
#include <TMath.h>
#include <TRandom.h>

// the energy window
constexpr double energy_window = 25;

double acceptance(double s1, double log_s2_s1)
{
    double s2 = std::pow(10, log_s2_s1) * s1;

    // modified for extended s1 and s2 region
    if (s1<3 || s1>100)
        return 0;

    if (s2<100 || s2 > 12000)
        return 0;

    if (log_s2_s1 < 1.074 - 0.595 * std::exp(-s1/6.38))
        return 0;

    return 1;
}

double GetBinCellSize (const TH2 * h)
{
    double x = h->GetXaxis()->GetBinUpEdge(1) - h->GetXaxis()->GetBinLowEdge(1);
    double y = h->GetYaxis()->GetBinUpEdge(1) - h->GetYaxis()->GetBinLowEdge(1);
    return x * y;
}

double LookupPDF (TH2F * h, double x, double y)
{
    double v = h->GetBinContent(h->FindBin(x, y));
    v = v/h->GetEntries()/GetBinCellSize(h);
    return v;
}


DataSet::DataSet()
    : days{0}, target_mass{100}, n_observed{0},
      data_file{""}, dm_pdf_file{""}, bkg_pdf_file{""},
      average_dm_efficiency{1}
{
}

DataSet::DataSet(const DataSet & set)
    : days{set.days}, target_mass{set.target_mass}, n_observed{set.n_observed},
      bkg_mdru{set.bkg_mdru},
      data_file{set.data_file}, dm_pdf_file{set.dm_pdf_file}, bkg_pdf_file{set.bkg_pdf_file},
      dm_pdf{set.dm_pdf}, bkg_pdf{set.bkg_pdf},
      x_{set.x_}, y_{set.y_},
      average_dm_efficiency{set.average_dm_efficiency},
      average_bkg_efficiencies{set.average_bkg_efficiencies}
{
    std::cout << "copy constructor!!!" << std::endl;
}

DataSet::~DataSet()
{
}

DataSet & DataSet::operator= (const DataSet& set)
{
    days = set.days;
    target_mass = set.target_mass;
    n_observed = set.n_observed;
    bkg_mdru = set.bkg_mdru;
    data_file = set.data_file;
    dm_pdf_file = set.dm_pdf_file;
    bkg_pdf_file = set.bkg_pdf_file;
    dm_pdf = set.dm_pdf;
    bkg_pdf = set.bkg_pdf;
    x_ = set.x_;
    y_ = set.y_;
    average_dm_efficiency = set.average_dm_efficiency;
    average_bkg_efficiencies = set.average_bkg_efficiencies;
    std::cout << " operator = " << std::endl;
    return * this;
}

DataSet & DataSet::operator= (DataSet&& set)
{
    days = std::move(set.days);
    target_mass = std::move(set.target_mass);
    n_observed = std::move(set.n_observed);
    bkg_mdru = std::move(set.bkg_mdru);
    data_file = std::move(set.data_file);
    dm_pdf_file = std::move(set.dm_pdf_file);
    bkg_pdf_file = std::move(set.bkg_pdf_file);
    dm_pdf = std::move(set.dm_pdf);
    bkg_pdf = std::move(set.bkg_pdf);
    x_ = std::move(set.x_);
    y_ = std::move(set.y_);
    average_dm_efficiency = std::move(set.average_dm_efficiency);
    average_bkg_efficiencies = std::move(set.average_bkg_efficiencies);
    std::cout << " operator = && " << std::endl;
    return * this;
}

void DataSet::LoadConfig (const string &data_description_file)
{
    std::string line;
    std::ifstream file(data_description_file);
    std::cout << "loading config for data set " << data_description_file << "..." << std::endl;
    set_name = data_description_file;
    while (file.good()) {
        std::getline(file, line);
        if (file.eof()) {
            break;
        }
        // find the first occurance of space or tab.
        size_t p = std::find_if(line.begin(), line.end(), [] (auto & c) {
                return (c == ' ' || c == '\t');
            }) - line.begin();
        if (p==line.size()) {
            continue;
        }
        auto key = line.substr(0, p);
        if (key == "days") {
            days = std::stod(line.data()+p+1);
        } else if (key == "target_mass") {
            target_mass = std::stod(line.data()+p+1);
        } else if (key == "data_file") {
            data_file = line.data()+p+1;
            data_file.erase(std::remove(data_file.begin(), data_file.end(), ' '), data_file.end());
        } else if (key == "dm_pdf") {
            dm_pdf_file = line.data()+p+1;
            dm_pdf_file.erase(std::remove(dm_pdf_file.begin(), dm_pdf_file.end(), ' '), dm_pdf_file.end());
        } else if (key == "bkg_pdf") {
            bkg_pdf_file = line.data()+p+1;
            bkg_pdf_file.erase(std::remove(bkg_pdf_file.begin(), bkg_pdf_file.end(), ' '), bkg_pdf_file.end());
        } else if (key == "bkg_mdru") {
            std::stringstream ss(line.data()+p+1);
            double v;
            while ( ss >> v) {
                bkg_mdru.push_back(v);
            }
        } else {
            std::cout << key << " is not used!" << std::endl;
        }
    }
    file.close();
}

void DataSet::PrintParameters()
{
    using std::cout;
    using std::endl;
    cout << " ///////////////////////////" << endl;
    cout << endl;
    cout << "    Data Set Parameters:" << endl << endl;
    cout << "    Days: " << days << endl;
    cout << "    Mass: " << target_mass << " kg" << endl;
    cout << "    Data File: " << data_file << endl;
    cout << "    DM PDF File: " << dm_pdf_file << endl;
    cout << "    BKG PDF File: " << bkg_pdf_file << endl;
    cout << "    BKG MDRU: ";
    for (auto & bm : bkg_mdru) {
        cout << bm << " ";
    }
    cout << endl << endl;
    cout << " ///////////////////////////" << endl;
}

void DataSet::ReadPDFs(const stringvec & bkg_components)
{
    TFile fdm(dm_pdf_file.data(), "READ");
    TIter next(fdm.GetListOfKeys());
    TKey * key;
    while((key = (TKey *) next())) {
        if (!(gROOT->GetClass(key->GetClassName())->InheritsFrom("TH2"))) {
            continue;
        }
        TH2F * histo = (TH2F *)key->ReadObj()->Clone();
        histo->SetDirectory(0);
        dm_pdf.push_back(histo);
    }
    fdm.Close();

    // get all the background components
    TFile fbkg(bkg_pdf_file.data(), "READ");
    for (auto & bc : bkg_components) {
        TH2F * histo = (TH2F*) fbkg.Get(bc.data())->Clone();
        histo->SetDirectory(0);
        bkg_pdf.push_back(histo);
    }
    fbkg.Close();
}

void DataSet::ReadData()
{
    std::string line;
    std::ifstream fin(data_file);
    std::regex data_format("([0-9]+\\.[0-9]+) \\* ([0-9]+\\.[0-9]+)", std::regex::extended);
    std::smatch results;
    while (fin.good()) {
        std::getline(fin, line);
        if (fin.eof())
            break;
        std::regex_search(line, results, data_format);
        if (!results.empty()) {
            {
                std::stringstream ss (results[1]);
                double v;
                ss >> v;
                x_.push_back(v);
            }
            {
                std::stringstream ss (results[2]);
                double v;
                ss >> v;
                y_.push_back(v);
            }
        }
    }
    fin.close();
    n_observed = x_.size();
}

void DataSet::ReadToyDataFromTree(TTree * tree, int round_number)
{
    std::string * pset{new std::string()};
    int round{0}, nEvents{0};
    double s1[10000], logs2s1[10000];
    tree->SetBranchAddress("round", &round);
    tree->SetBranchAddress("set", &pset);
    tree->SetBranchAddress("nEvents", &nEvents);
    tree->SetBranchAddress("s1", &s1);
    tree->SetBranchAddress("logs2s1", &logs2s1);
    // loop over the trees to find out the round and set
    long nEntries = tree->GetEntries();
    x_.clear();
    y_.clear();
    std::cout << "set name = " << set_name << " and round = " << round_number << std::endl;
    for (long iEntry = 0; iEntry<nEntries; ++iEntry) {
        tree->GetEntry(iEntry);
        if (round!=round_number || set_name != *pset) {
            pset->clear();
            continue;
        }
        x_.insert(x_.begin(), s1, s1+nEvents);
        y_.insert(y_.begin(), logs2s1, logs2s1+nEvents);
        n_observed = nEvents;
        pset->clear();
        break;
    }
    std::cout << "load " << x_.size() << " points from tree." << std::endl;
    // don't forget to delete the pointer
    delete pset;
}

double DataSet::averageEfficiency(TH2F *pdf)
{
    double sum = 0;
    double totalEntries = pdf->GetEntries();
    for (int i=1; i<pdf->GetNbinsX(); ++i) {
        for (int j=1; j<pdf->GetNbinsY(); ++j) {
            double weight = pdf->GetBinContent(i, j)/totalEntries;
            double x = pdf->GetXaxis()->GetBinCenter(i);
            double y = pdf->GetYaxis()->GetBinCenter(j);
            sum += weight * acceptance(x, y);
        }
    }
    return sum;
}

TH2F * DataSet::FindMatchDMPdfByDeltaM(double dm) const
{
    TH2F * pdf = 0;
    std::regex pattern("dm([0-9\\.]+)");
    std::smatch results;
    double ddm = 10000;
    for (auto p: dm_pdf) {
        std::string name = p->GetName();
        std::regex_search(name, results, pattern);
        if (!results.empty()) {
            double dm_p = std::stod(results[1]);
            if (std::abs(dm - dm_p) < ddm) {
                ddm = std::abs(dm - dm_p);
                pdf = p;
            }
        }
    }
    return pdf;
}

void DataSet::computeAverageEfficiency(double delta_mass)
{
    TH2F * match_dm_pdf = FindMatchDMPdfByDeltaM(delta_mass);
    average_dm_efficiency = averageEfficiency(match_dm_pdf);
    computeAverageBkgEfficiencies();
//    std::cout << "average dm efficency for delta mass " << delta_mass << " is " << average_dm_efficiency << std::endl;
}

void DataSet::computeAverageBkgEfficiencies()
{
    average_bkg_efficiencies.clear();
    for (auto & bp: bkg_pdf) {
        average_bkg_efficiencies.push_back(averageEfficiency(bp));
    }
}

double DataSet::LogLikelihood (double delta_mass,
                               double signal_count_per_day_per_kg,
                               double delta_dm,
                               const std::vector<double> & delta_bkg)
{
    double dm_raw = signal_count_per_day_per_kg * (1 + delta_dm) * days * target_mass;
    double dm_detected = dm_raw * average_dm_efficiency;

    std::vector<double> bkg_raw, bkg_detected;

    auto n_expected = dm_detected;
//    std::cout << "signal count per day per kg: "  << signal_count_per_day_per_kg << " "  << average_dm_efficiency << " " << delta_dm << " " << dm_raw << " " << dm_detected << std::endl;
    for (size_t i=0; i<delta_bkg.size(); ++i) {
        bkg_raw.emplace_back(bkg_mdru[i] * (1 + delta_bkg[i]) * days * target_mass * energy_window / 1e3);
//        std::cout << "bkg info: " << bkg_mdru[i] << ", " << delta_bkg[i] << ", " << average_bkg_efficiencies[i] << std::endl;
        bkg_detected.emplace_back(bkg_raw[i] * average_bkg_efficiencies[i]);
        n_expected += bkg_detected[i];
    }

    double log_likelihood = 0;
    auto match_dm_pdf = FindMatchDMPdfByDeltaM(delta_mass);
    for (int i =0; i<n_observed; ++i) {
        double prob = 0;
        double pdf = LookupPDF(match_dm_pdf, x_[i], y_[i]);
        double factor = acceptance(x_[i], y_[i]);
        prob = dm_raw * pdf * factor / n_expected;
        for (size_t j=0; j<bkg_pdf.size(); ++j) {
            pdf = LookupPDF(bkg_pdf[j], x_[i], y_[i]);
            prob += bkg_raw[j] * pdf * factor / n_expected;;
        }
        if (prob <= 0)
            prob = 1e-31;
//        std::cout << "prob = " << prob << std::endl;
        log_likelihood += std::log(prob);
    }
    log_likelihood += std::log(TMath::Poisson(n_observed, n_expected));
//    std::cout << "log likelihood: " << log_likelihood << " (" << n_expected << ", " << n_observed << ")" << std::endl;
    return log_likelihood;
}

void DataSet::NormalizePdf(TH2F * pdf)
{
    std::unique_ptr<TH2F> tmp((TH2F*)pdf->Clone());
    tmp->Reset();

    for (int i=1; i<=pdf->GetNbinsX(); ++i) {
        for (int j=1; j<=pdf->GetNbinsY(); ++j) {
            double x = pdf->GetXaxis()->GetBinCenter(i);
            double y = pdf->GetYaxis()->GetBinCenter(j);
            double weight = pdf->GetBinContent(i, j);
            if (acceptance(x ,y) > 0) {
                tmp->SetBinContent(i, j, weight);
            }
        }
    }
    tmp->Scale(1/tmp->Integral());
    pdf->Reset();
    for (int i=1; i<=pdf->GetNbinsX(); ++i) {
        for (int j=1; j<=pdf->GetNbinsY(); ++j) {
            pdf->SetBinContent(i, j, tmp->GetBinContent(i, j));
        }
    }
}

void DataSet::dumpBestFit(double delta_mass,
                          const std::vector<double> parameters,
                          TFile *file,
                          const std::string tag)
{
    double signal_count_per_day_per_kg = parameters[0];
    double delta_dm = parameters[1];

    double dm_raw = signal_count_per_day_per_kg * (1 + delta_dm) * days * target_mass;
    double dm_detected = dm_raw * average_dm_efficiency;
    std::unique_ptr<TH2F> hist_dm_bestfit ((TH2F *)(FindMatchDMPdfByDeltaM(delta_mass)->Clone()));
    hist_dm_bestfit->AddDirectory(kFALSE);
    hist_dm_bestfit->SetName((tag + hist_dm_bestfit->GetName()).c_str());
    NormalizePdf(hist_dm_bestfit.get());
    hist_dm_bestfit->Scale(dm_detected);
    file->cd();
    hist_dm_bestfit->Write();

    std::unique_ptr<TH2F> hist_data((TH2F*)(hist_dm_bestfit->Clone()));
    hist_data->AddDirectory(kFALSE);
    hist_data->SetName((tag + "data").c_str());
    hist_data->Reset();
    for (int i=0; i<n_observed; ++i) {
        hist_data->Fill(x_[i], y_[i]);
    }
    file->cd();
    hist_data->Write();

    for (size_t i = 0; i < bkg_pdf.size(); ++i) {
        double bkg_raw = bkg_mdru[i] * (1 + parameters[i+2]) * days * target_mass * energy_window / 1e3;
        double bkg_detected = bkg_raw * average_bkg_efficiencies[i];
        std::unique_ptr<TH2F> hist_bkg_bestfit ((TH2F *)(bkg_pdf[i]->Clone()));
        hist_bkg_bestfit->AddDirectory(kFALSE);
        hist_bkg_bestfit->SetName((tag + bkg_pdf[i]->GetName()).c_str());
        NormalizePdf(hist_bkg_bestfit.get());
        hist_bkg_bestfit->Scale(bkg_detected);
        file->cd();
        hist_bkg_bestfit->Write();
    }
}

int DataSet::GetMCDetectedBkgCount(int index, double delta) const
{
    if ((unsigned int)index >= bkg_mdru.size()) {
        return 0;
    }
    delta = gRandom->Gaus(0, delta);
    double nbkg_raw = bkg_mdru[index] * (1 + delta) * days * target_mass * energy_window / 1e3;
    double nbkg_exp = nbkg_raw * average_bkg_efficiencies[index];
    double nbkg_det = gRandom->Poisson(nbkg_exp);
    return (int)nbkg_det;
}

int DataSet::GetMCDetectedDMCount(double ref_count, double delta) const
{
    delta = gRandom->Gaus(0, delta);
    double ndm_raw = ref_count * (1 + delta) * days * target_mass;
    double ndm_exp = ndm_raw * average_dm_efficiency;
    double ndm_det = gRandom->Poisson(ndm_exp);
    return (int)ndm_det;
}

void DataSet::SampleFromPDF(TH2F *pdf, double &x, double &y)
{
    do {
        pdf->GetRandom2(x, y);
    } while(acceptance(x, y)<1);
}
