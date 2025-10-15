#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <regex>
#include <numeric>
#include <filesystem>
#include <map>
#include <set>
#include <algorithm>
#include <cmath>
#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TLatex.h"
#include "TFile.h"

namespace fs = std::filesystem;

const int nDiscriminators = 31;
const int nASICs = 8;
const int nChannels = 128;
const int totalChannels = nASICs * nChannels;

using DiscriminatorVector = std::vector<float>;
using ChannelVector = std::vector<DiscriminatorVector>;
using ASICVector = std::vector<ChannelVector>;

// Mapping for physical ASIC order
std::map<int,int> sortMapFEB_A = {{0,7},{1,6},{2,5},{3,4},{4,3},{5,2},{6,1},{7,0}};
std::map<int,int> sortMapFEB_B = {{0,1},{1,0},{2,3},{3,2},{4,5},{5,4},{6,7},{7,6}};

// Struct to store header info
struct FileHeaderInfo {
    int vref_p = -1;
    int vref_n = -1;
    int thr2_glb = -1;
    int vref_t = -1;
    int vref_t_range = -1;
    int pol = -1;
    std::pair<int,int> adc_range = {-1,-1};
    int fast_disc = -1;
};

void load_adc_data(ASICVector& adc_data_elect,
                   ASICVector& adc_data_holes,
                   std::map<std::string, FileHeaderInfo>& header_info,
                   std::set<int>& skip_elect,
                   std::set<int>& skip_holes,
                   const std::string& folderPath = "./",
                   const std::string& feb_side = "A") 
{
    const float CUT_MIN = 60.0f;
    const float CUT_MAX = 180.0f;

    std::regex re_hw("HW_([0-7])");
    std::regex re_side("(elect|holes)");
    std::smatch match;

    // Header regexes
    std::regex re_vref_p("Vref_p\\s+(\\d+)");
    std::regex re_vref_n("Vref_n\\s+(\\d+)");
    std::regex re_thr2("Thr2_glb\\s+(\\d+)");
    std::regex re_vref_t("Vref_t\\s+(\\d+)");
    std::regex re_vref_t_range("Vref_t_range\\s+(\\d+)");
    std::regex re_pol("Pol\\s+(\\d+)");
    std::regex re_adc_range("ADC range\\s+(\\d+)-(\\d+)");
    std::regex re_fast_disc("FAST disc\\s+(\\d+)");

    // Count files per hw/side
    std::map<int,int> ecounts, hcounts;
    for (const auto& entry : fs::directory_iterator(folderPath)) {
        std::string fname = entry.path().string();
        if (fname.find(".txt") == std::string::npos) continue;

        int hw=-1;
        std::string side;
        if (std::regex_search(fname, match, re_hw)) hw = std::stoi(match[1]);
        if (std::regex_search(fname, match, re_side)) side = match[1];
        if (hw<0 || side.empty()) continue;

        if (side=="elect") ecounts[hw]++;
        else if (side=="holes") hcounts[hw]++;
    }

    // Determine skips per side using HW IDs
    for(int hw=0; hw<nASICs; ++hw) {
        if(ecounts[hw]!=1) {
            std::cout << "Skipping HW_" << hw << ": elect files = " << ecounts[hw] << " (need 1)\n";
            skip_elect.insert(hw);
        }
        if(hcounts[hw]!=1) {
            std::cout << "Skipping HW_" << hw << ": holes files = " << hcounts[hw] << " (need 1)\n";
            skip_holes.insert(hw);
        }
    }

    // Load allowed files
    for (const auto& entry : fs::directory_iterator(folderPath)) {
        std::string fname = entry.path().string();
        if (fname.find(".txt") == std::string::npos) continue;

        int hw=-1;
        std::string side;
        if (std::regex_search(fname, match, re_hw)) hw = std::stoi(match[1]);
        if (std::regex_search(fname, match, re_side)) side = match[1];
        if (hw<0 || side.empty()) continue;

        if((side=="elect" && skip_elect.count(hw)) ||
           (side=="holes" && skip_holes.count(hw))) continue;

        std::ifstream file(fname);
        if(!file.is_open()) { std::cerr << "Cannot open " << fname << std::endl; continue; }

        // Parse header
        FileHeaderInfo header;
        std::string line;
        while(std::getline(file,line)) {
            if(line.find("ch:")==0) break;
            std::smatch m;
            if(std::regex_search(line,m,re_vref_p)) header.vref_p = std::stoi(m[1]);
            else if(std::regex_search(line,m,re_vref_n)) header.vref_n = std::stoi(m[1]);
            else if(std::regex_search(line,m,re_thr2)) header.thr2_glb = std::stoi(m[1]);
            else if(std::regex_search(line,m,re_vref_t)) header.vref_t = std::stoi(m[1]);
            else if(std::regex_search(line,m,re_vref_t_range)) header.vref_t_range = std::stoi(m[1]);
            else if(std::regex_search(line,m,re_pol)) header.pol = std::stoi(m[1]);
            else if(std::regex_search(line,m,re_adc_range)) header.adc_range = {std::stoi(m[1]), std::stoi(m[2])};
            else if(std::regex_search(line,m,re_fast_disc)) header.fast_disc = std::stoi(m[1]);
        }

        std::string key = "HW_" + std::to_string(hw) + "_" + side;
        header_info[key] = header;

        // Rewind for channel data
        file.clear();
        file.seekg(0);
        while(std::getline(file,line) && line.find("ch:")!=0) {}

        do {
            if(line.find("ch:")!=0) continue;
            std::istringstream iss(line.substr(3));
            int ch; iss >> ch;
            if(ch<0||ch>=nChannels) continue;

            DiscriminatorVector vals(nDiscriminators,-1.f);
            for(int i=0;i<nDiscriminators;++i) {
                iss >> vals[i];
                // Apply cut: mark as invalid if outside [CUT_MIN, CUT_MAX]
                if (vals[i] < CUT_MIN || vals[i] > CUT_MAX) vals[i] = -1.f;
            }

            if(side=="elect") adc_data_elect[hw][ch] = vals;
            else adc_data_holes[hw][ch] = vals;

        } while(std::getline(file,line));
    }
}

ASICVector reorder_asics(const ASICVector& data, const std::map<int,int>& sortMap) {
    ASICVector reordered(nASICs, ChannelVector(nChannels, DiscriminatorVector(nDiscriminators,-1.f)));
    for(int hw=0; hw<nASICs; ++hw) {
        int sortedIndex = sortMap.at(hw);
        reordered[sortedIndex] = data[hw];
    }
    return reordered;
}

void MODULE_NAME_statistics(const ASICVector& adc_data_elect, 
                       const ASICVector& adc_data_holes, 
                       const std::set<int>& skip_elect_sorted,
                       const std::set<int>& skip_holes_sorted,
                       const std::map<std::string, FileHeaderInfo>& header_info,
                       const std::string& feb_side)
{
    TCanvas* c = new TCanvas("asic_mean_stddev", "ASIC Mean and StdDev LSB", 1400, 700);

    // Determine P-side and N-side mapping
    const auto& pMap = (feb_side == "A") ? sortMapFEB_B : sortMapFEB_A; // P-side
    const auto& nMap = (feb_side == "A") ? sortMapFEB_A : sortMapFEB_B; // N-side

    // Vectors for TGraphErrors
    std::vector<double> xElect, meanElect, stdElect;
    std::vector<double> xHoles, meanHoles, stdHoles;
    std::vector<double> allElectrons, allHoles;

    auto compute_stats = [](const DiscriminatorVector& v, double& mean, double& stddev) {
        std::vector<float> filtered;
        for (float val : v) if (val > 0) filtered.push_back(val);
        if (filtered.empty()) { mean = stddev = NAN; return; }
        double sum = std::accumulate(filtered.begin(), filtered.end(), 0.0);
        mean = sum / filtered.size();
        if (filtered.size() > 1) {
            double sq_sum = 0;
            for (float val : filtered) sq_sum += (val - mean)*(val - mean);
            stddev = std::sqrt(sq_sum / (filtered.size() - 1)) / std::sqrt(filtered.size());
        } else stddev = 0;
    };

    // Loop over sorted ASICs
    for (int sorted_hw = 0; sorted_hw < nASICs; ++sorted_hw) {

        // Electrons (P-side)
        int hw_p = -1;
        for (const auto& [hw, idx] : pMap) { if (idx == sorted_hw) { hw_p = hw; break; } }
        if (hw_p != -1 && !skip_elect_sorted.count(sorted_hw)) {
            for (int ch = 0; ch < nChannels; ++ch) {
                const auto& ve = adc_data_elect[sorted_hw][ch];
                double mean, stddev;
                compute_stats(ve, mean, stddev);
                if (!std::isnan(mean) && mean > 0) {
                    int idx = sorted_hw * nChannels + ch;
                    xElect.push_back(idx);
                    meanElect.push_back(mean);
                    stdElect.push_back(stddev);
                    for (float val : ve) if (val > 0) allElectrons.push_back(val);
                }
            }
        }

        // Holes (N-side)
        int hw_n = -1;
        for (const auto& [hw, idx] : nMap) { if (idx == sorted_hw) { hw_n = hw; break; } }
        if (hw_n != -1 && !skip_holes_sorted.count(sorted_hw)) {
            for (int ch = 0; ch < nChannels; ++ch) {
                const auto& vh = adc_data_holes[sorted_hw][ch];
                double mean, stddev;
                compute_stats(vh, mean, stddev);
                if (!std::isnan(mean) && mean > 0) {
                    int idx = sorted_hw * nChannels + ch;
                    xHoles.push_back(idx);
                    meanHoles.push_back(mean);
                    stdHoles.push_back(stddev);
                    for (float val : vh) if (val > 0) allHoles.push_back(val);
                }
            }
        }
    }

    // Draw TGraphErrors
    TGraphErrors* gElect = new TGraphErrors(xElect.size(), xElect.data(), meanElect.data(), nullptr, stdElect.data());
    TGraphErrors* gHoles = new TGraphErrors(xHoles.size(), xHoles.data(), meanHoles.data(), nullptr, stdHoles.data());

    gElect->SetMarkerStyle(20); gElect->SetMarkerColor(kRed); gElect->SetLineColor(kRed);
    gHoles->SetMarkerStyle(24); gHoles->SetMarkerColor(kBlue); gHoles->SetLineColor(kBlue);

    gElect->SetTitle("MODULE_NAME Mean Trim Statistics;Channel;Mean [LSB]");
    gElect->GetXaxis()->SetLimits(0, totalChannels - 1);
    gElect->SetMinimum(105); gElect->SetMaximum(140);
    gElect->Draw("AP");
    gHoles->Draw("P SAME");
    c->Update();  // updates pad coordinates

    // Get pad limits
    double ymin = gPad->GetUymin();
    double ymax = gPad->GetUymax();

    // Instantiate TLatex and TLine
    TLatex tex; tex.SetTextAngle(0); tex.SetTextSize(0.025); tex.SetTextAlign(22);
    TLine line; line.SetLineStyle(2); line.SetLineColor(kGray+1);

    // Draw dividers and Vref labels
    for (int sorted_hw = 0; sorted_hw < nASICs; ++sorted_hw) {
        double xpos = (sorted_hw + 1) * nChannels - 0.5;
        line.DrawLine(xpos, ymin, xpos, ymax);

        // Find HW IDs
        int hw_p = -1, hw_n = -1;
        for (const auto& [hw, idx] : pMap) if (idx == sorted_hw) { hw_p = hw; break; }
        for (const auto& [hw, idx] : nMap) if (idx == sorted_hw) { hw_n = hw; break; }

        double text_x = xpos - nChannels/2.0;

        // Electron Vref
        if (hw_p != -1) {
            std::string key = "HW_" + std::to_string(hw_p) + "_elect";
            if (header_info.count(key)) {
                auto& info = header_info.at(key);
                tex.SetTextColor(kRed);
                tex.DrawLatex(text_x, ymin + 0.05*(ymax-ymin), Form("Vp=%d Vn=%d", info.vref_p, info.vref_n));
            }
        }

        // Hole Vref
        if (hw_n != -1) {
            std::string key = "HW_" + std::to_string(hw_n) + "_holes";
            if (header_info.count(key)) {
                auto& info = header_info.at(key);
                tex.SetTextColor(kBlue);
                tex.DrawLatex(text_x, ymin + 0.02*(ymax-ymin), Form("Vp=%d Vn=%d", info.vref_p, info.vref_n));
            }
        }
    }

    // Compute overall stats for TLegend
    auto getStats = [](const std::vector<double>& vec, double& mean, double& stddev) {
        if (vec.empty()) { mean = 0; stddev = 0; return; }
        double sum = std::accumulate(vec.begin(), vec.end(), 0.0);
        mean = sum / vec.size();
        double sqsum = 0;
        for (auto v : vec) sqsum += (v - mean)*(v - mean);
        stddev = std::sqrt(sqsum / vec.size());
    };

    double meanE=0, stdE=0, meanH=0, stdH=0;
    getStats(allElectrons, meanE, stdE);
    getStats(allHoles, meanH, stdH);

    // Draw legend
    TLegend* legend = new TLegend(0.12, 0.80, 0.45, 0.89); 
    legend->SetTextSize(0.025);
    legend->AddEntry(gElect, "Electrons", "lp"); 
    legend->AddEntry(gHoles, "Holes", "lp");
    legend->AddEntry((TObject*)nullptr, Form("Electrons: %.2f +/- %.2f [LSB]", meanE,stdE),"");
    legend->AddEntry((TObject*)nullptr, Form("Holes:     %.2f +/- %.2f [LSB]", meanH,stdH),"");
    legend->Draw();

    // Save plots
    c->SaveAs("MODULE_NAME_statistics.root");
    c->SaveAs("MODULE_NAME_statistics.pdf");
}

void mean_adc_discriminators(const ASICVector& adc_data_elect,
                             const ASICVector& adc_data_holes,
                             const std::set<int>& skip_elect,
                             const std::set<int>& skip_holes) 
{
    std::vector<double> allElectronsValues;
    std::vector<double> allHolesValues;

    std::vector<double> x_vals, mean_e, mean_h, err_e, err_h;

    auto compute_stats = [](const std::vector<double>& values, double& mean, double& stddev, double& sem) {
        if (values.empty()) { mean = stddev = sem = 0; return; }
        double sum = std::accumulate(values.begin(), values.end(), 0.0);
        mean = sum / values.size();
        if (values.size() > 1) {
            double sqsum = 0;
            for (auto v : values) sqsum += (v - mean) * (v - mean);
            stddev = std::sqrt(sqsum / (values.size() - 1));
            sem = stddev / std::sqrt(values.size());
        } else { stddev = sem = 0; }
    };

    for (int d = 0; d < nDiscriminators; ++d) {
        std::vector<double> values_e, values_h;

        for (int hw = 0; hw < nASICs; ++hw) {
            bool skip_e = skip_elect.count(hw);
            bool skip_h = skip_holes.count(hw);
            
            for (int ch = 0; ch < nChannels; ++ch) {
                if (!skip_e) {
                    float ve = adc_data_elect[hw][ch][d];
                    if (ve != -1) { values_e.push_back(ve); allElectronsValues.push_back(ve); }
                }
                if (!skip_h) {
                    float vh = adc_data_holes[hw][ch][d];
                    if (vh != -1) { values_h.push_back(vh); allHolesValues.push_back(vh); }
                }
            }
        }

        double m_e=0, s_e=0, se_e=0, m_h=0, s_h=0, se_h=0;
        compute_stats(values_e, m_e, s_e, se_e);
        compute_stats(values_h, m_h, s_h, se_h);

        x_vals.push_back(d);
        mean_e.push_back(m_e); err_e.push_back(se_e);
        mean_h.push_back(m_h); err_h.push_back(se_h);
    }

    // Convert to TGraphErrors
    auto g_elect = new TGraphErrors(nDiscriminators, x_vals.data(), mean_e.data(), nullptr, err_e.data());
    auto g_holes = new TGraphErrors(nDiscriminators, x_vals.data(), mean_h.data(), nullptr, err_h.data());

    g_elect->SetLineColor(kRed);   g_elect->SetMarkerColor(kRed);   g_elect->SetMarkerStyle(20);
    g_holes->SetLineColor(kBlue);  g_holes->SetMarkerColor(kBlue);  g_holes->SetMarkerStyle(24);

    TCanvas* c = new TCanvas("mean_adc_per_discriminator", "MODULE_NAME Mean LSB per Discriminator", 1000, 600);
    g_elect->GetXaxis()->SetTitle("Discriminator");
    g_elect->GetYaxis()->SetTitle("Mean [LSB]");
    g_elect->GetYaxis()->SetRangeUser(105, 140);
    g_elect->SetTitle("MODULE_NAME Mean LSB per Discriminator;Discriminator;Mean [LSB]");
    g_holes->SetTitle("MODULE_NAME Mean LSB per Discriminator;Discriminator;Mean [LSB]");

    g_elect->Draw("AP");
    g_holes->Draw("P SAME");

    auto legend = new TLegend(0.7,0.8,0.9,0.9);
    legend->AddEntry(g_elect,"Electrons","lp");
    legend->AddEntry(g_holes,"Holes","lp");
    legend->Draw();

    c->SaveAs("MODULE_NAME_mean_adc_per_discriminator.root");
    c->SaveAs("MODULE_NAME_mean_adc_per_discriminator.pdf");
}


void mean_adc_discriminators_odd_even(const ASICVector& adc_data_elect,
                                      const ASICVector& adc_data_holes,
                                      const std::set<int>& skip_elect,
                                      const std::set<int>& skip_holes) 
{
    std::vector<double> x_vals, e_even, e_odd, h_even, h_odd;
    std::vector<double> err_e_even, err_e_odd, err_h_even, err_h_odd;

    auto compute_stats = [](const std::vector<double>& values, double& mean, double& stddev, double& sem) {
        if (values.empty()) { mean = stddev = sem = 0; return; }
        double sum = std::accumulate(values.begin(), values.end(), 0.0);
        mean = sum / values.size();
        if (values.size() > 1) {
            double sqsum = 0;
            for (auto v : values) sqsum += (v - mean) * (v - mean);
            stddev = std::sqrt(sqsum / (values.size() - 1));
            sem = stddev / std::sqrt(values.size());
        } else { stddev = sem = 0; }
    };

    for (int d = 0; d < nDiscriminators; ++d) {
        std::vector<double> v_e_even, v_e_odd, v_h_even, v_h_odd;

        for (int hw = 0; hw < nASICs; ++hw) {
            bool skip_e = skip_elect.count(hw);
            bool skip_h = skip_holes.count(hw);
            
            for (int ch = 0; ch < nChannels; ++ch) {
                if (!skip_e) {
                    float ve = adc_data_elect[hw][ch][d];
                    if (ve != -1) ((ch % 2 == 0) ? v_e_even : v_e_odd).push_back(ve);
                }
                if (!skip_h) {
                    float vh = adc_data_holes[hw][ch][d];
                    if (vh != -1) ((ch % 2 == 0) ? v_h_even : v_h_odd).push_back(vh);
                }
            }
        }

        double m_e_e=0,s_e_e=0,se_e_e=0, m_e_o=0,s_e_o=0,se_e_o=0;
        double m_h_e=0,s_h_e=0,se_h_e=0, m_h_o=0,s_h_o=0,se_h_o=0;
        compute_stats(v_e_even, m_e_e, s_e_e, se_e_e);
        compute_stats(v_e_odd,  m_e_o, s_e_o, se_e_o);
        compute_stats(v_h_even, m_h_e, s_h_e, se_h_e);
        compute_stats(v_h_odd,  m_h_o, s_h_o, se_h_o);

        x_vals.push_back(d);
        e_even.push_back(m_e_e); err_e_even.push_back(se_e_e);
        e_odd.push_back(m_e_o);  err_e_odd.push_back(se_e_o);
        h_even.push_back(m_h_e); err_h_even.push_back(se_h_e);
        h_odd.push_back(m_h_o);  err_h_odd.push_back(se_h_o);
    }

    auto g_e_even = new TGraphErrors(nDiscriminators, x_vals.data(), e_even.data(), nullptr, err_e_even.data());
    auto g_e_odd  = new TGraphErrors(nDiscriminators, x_vals.data(), e_odd.data(),  nullptr, err_e_odd.data());
    auto g_h_even = new TGraphErrors(nDiscriminators, x_vals.data(), h_even.data(), nullptr, err_h_even.data());
    auto g_h_odd  = new TGraphErrors(nDiscriminators, x_vals.data(), h_odd.data(),  nullptr, err_h_odd.data());

    g_e_even->SetLineColor(kRed);      g_e_even->SetMarkerColor(kRed);     g_e_even->SetMarkerStyle(20);
    g_e_odd->SetLineColor(kMagenta);   g_e_odd->SetMarkerColor(kMagenta);  g_e_odd->SetMarkerStyle(21);
    g_h_even->SetLineColor(kBlue);     g_h_even->SetMarkerColor(kBlue);    g_h_even->SetMarkerStyle(24);
    g_h_odd->SetLineColor(kGreen+2);   g_h_odd->SetMarkerColor(kGreen+2);  g_h_odd->SetMarkerStyle(25);

    TCanvas* c = new TCanvas("mean_adc_odd_even", "MODULE_NAME Mean ADC per Discriminator (Odd/Even)", 1200, 700);
    g_e_even->GetXaxis()->SetTitle("Discriminator");
    g_e_even->GetYaxis()->SetTitle("Mean [LSB]");
    g_e_even->GetYaxis()->SetRangeUser(105, 140);
    g_e_even->SetTitle("MODULE_NAME Mean LSB per Discriminator;Discriminator;Mean [LSB]");
    g_h_even->SetTitle("MODULE_NAME Mean LSB per Discriminator;Discriminator;Mean [LSB]");
    
    g_e_even->Draw("AP");
    g_e_odd->Draw("P SAME");
    g_h_even->Draw("P SAME");
    g_h_odd->Draw("P SAME");

    auto legend = new TLegend(0.7,0.8,0.9,0.9);
    legend->AddEntry(g_e_even,"Electrons Even Ch","lp");
    legend->AddEntry(g_e_odd,"Electrons Odd Ch","lp");
    legend->AddEntry(g_h_even,"Holes Even Ch","lp");
    legend->AddEntry(g_h_odd,"Holes Odd Ch","lp");
    legend->Draw();

    c->SaveAs("MODULE_NAME_mean_adc_odd_even.root");
    c->SaveAs("MODULE_NAME_mean_adc_odd_even.pdf");
}

void mean_adc_discriminators_perASIC_oddeven(const ASICVector& adc_data_elect,
                                             const ASICVector& adc_data_holes,
                                             const std::set<int>& skip_elect,
                                             const std::set<int>& skip_holes)
{
    auto compute_stats = [](const std::vector<double>& values, double& mean, double& sem) {
        if (values.empty()) { mean = sem = 0; return; }
        mean = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
        if (values.size() > 1) {
            double sqsum = 0;
            for (auto v : values) sqsum += (v - mean) * (v - mean);
            double stddev = std::sqrt(sqsum / (values.size() - 1));
            sem = stddev / std::sqrt(values.size());
        } else sem = 0;
    };

    // x-axis = discriminators
    std::vector<double> x_vals(nDiscriminators);
    std::iota(x_vals.begin(), x_vals.end(), 0);

    // Storage for graphs
    std::vector<TGraphErrors*> gElect_even(nASICs), gElect_odd(nASICs);
    std::vector<TGraphErrors*> gHoles_even(nASICs), gHoles_odd(nASICs);

    for (int hw = 0; hw < nASICs; ++hw) {
        if (!skip_elect.count(hw)) {
            std::vector<double> mean_even, err_even, mean_odd, err_odd;
            for (int d = 0; d < nDiscriminators; ++d) {
                std::vector<double> vals_even, vals_odd;
                for (int ch = 0; ch < nChannels; ++ch) {
                    float v = adc_data_elect[hw][ch][d];
                    if (v != -1) ((ch % 2 == 0) ? vals_even : vals_odd).push_back(v);
                }
                double m_even, se_even, m_odd, se_odd;
                compute_stats(vals_even, m_even, se_even);
                compute_stats(vals_odd,  m_odd,  se_odd);
                mean_even.push_back(m_even); err_even.push_back(se_even);
                mean_odd.push_back(m_odd);   err_odd.push_back(se_odd);
            }
            gElect_even[hw] = new TGraphErrors(nDiscriminators, x_vals.data(), mean_even.data(), nullptr, err_even.data());
            gElect_odd[hw]  = new TGraphErrors(nDiscriminators, x_vals.data(), mean_odd.data(),  nullptr, err_odd.data());

            gElect_even[hw]->SetLineColor(kRed);    gElect_even[hw]->SetMarkerColor(kRed);    gElect_even[hw]->SetMarkerStyle(20);
            gElect_odd[hw]->SetLineColor(kMagenta); gElect_odd[hw]->SetMarkerColor(kMagenta); gElect_odd[hw]->SetMarkerStyle(21);
        }

        if (!skip_holes.count(hw)) {
            std::vector<double> mean_even, err_even, mean_odd, err_odd;
            for (int d = 0; d < nDiscriminators; ++d) {
                std::vector<double> vals_even, vals_odd;
                for (int ch = 0; ch < nChannels; ++ch) {
                    float v = adc_data_holes[hw][ch][d];
                    if (v != -1) ((ch % 2 == 0) ? vals_even : vals_odd).push_back(v);
                }
                double m_even, se_even, m_odd, se_odd;
                compute_stats(vals_even, m_even, se_even);
                compute_stats(vals_odd,  m_odd,  se_odd);
                mean_even.push_back(m_even); err_even.push_back(se_even);
                mean_odd.push_back(m_odd);   err_odd.push_back(se_odd);
            }
            gHoles_even[hw] = new TGraphErrors(nDiscriminators, x_vals.data(), mean_even.data(), nullptr, err_even.data());
            gHoles_odd[hw]  = new TGraphErrors(nDiscriminators, x_vals.data(), mean_odd.data(),  nullptr, err_odd.data());

            gHoles_even[hw]->SetLineColor(kBlue);   gHoles_even[hw]->SetMarkerColor(kBlue);   gHoles_even[hw]->SetMarkerStyle(24);
            gHoles_odd[hw]->SetLineColor(kGreen+2); gHoles_odd[hw]->SetMarkerColor(kGreen+2); gHoles_odd[hw]->SetMarkerStyle(25);
        }
    }

    // Canvas for Electrons
    TCanvas* c1 = new TCanvas("c_elect_oddeven","MODULE_NAME Mean ADC per Discriminator per ASIC - Electrons",1200,800);
    c1->Divide(4,2);
    for (int hw = 0; hw < nASICs; ++hw) {
        c1->cd(hw+1);
        if (gElect_even[hw] && gElect_odd[hw]) {
            gElect_even[hw]->SetTitle(Form("Electrons ASIC %d;Discriminator;Mean [LSB]", hw));
            gElect_even[hw]->GetYaxis()->SetRangeUser(100,160);
            gElect_even[hw]->Draw("AP");
            gElect_odd[hw]->Draw("P SAME");
            auto legend = new TLegend(0.6,0.8,0.9,0.9);
            legend->AddEntry(gElect_even[hw],"Even ch","lp");
            legend->AddEntry(gElect_odd[hw],"Odd ch","lp");
            legend->Draw();
        }
    }
    c1->SaveAs("MODULE_NAME_mean_adc_perASIC_oddeven_electrons.root");
    c1->SaveAs("MODULE_NAME_mean_adc_perASIC_oddeven_electrons.pdf");

    // Canvas for Holes
    TCanvas* c2 = new TCanvas("c_holes_oddeven","MODULE_NAME Mean ADC per Discriminator per ASIC - Holes",1200,800);
    c2->Divide(4,2);
    for (int hw = 0; hw < nASICs; ++hw) {
        c2->cd(hw+1);
        if (gHoles_even[hw] && gHoles_odd[hw]) {
            gHoles_even[hw]->SetTitle(Form("Holes ASIC %d;Discriminator;Mean [LSB]", hw));
            gHoles_even[hw]->GetYaxis()->SetRangeUser(100,160);
            gHoles_even[hw]->Draw("AP");
            gHoles_odd[hw]->Draw("P SAME");
            auto legend = new TLegend(0.6,0.8,0.9,0.9);
            legend->AddEntry(gHoles_even[hw],"Even ch","lp");
            legend->AddEntry(gHoles_odd[hw],"Odd ch","lp");
            legend->Draw();
        }
    }
    c2->SaveAs("MODULE_NAME_mean_adc_perASIC_oddeven_holes.root");
    c2->SaveAs("MODULE_NAME_mean_adc_perASIC_oddeven_holes.pdf");
}


void discriminator_raw_adc_distributions(const ASICVector& adc_data_elect, 
                                         const ASICVector& adc_data_holes, 
                                         const std::set<int>& skip_elect,
                                         const std::set<int>& skip_holes) {
    TCanvas* c = new TCanvas("discriminator_adc_distributions", 
                             "MODULE_NAME_Discriminator Distributions", 
                             1000, 600);
    c->Divide(6, 6);

    const int nDiscriminators = adc_data_elect[0][0].size();
    TGraph* g_skew_elect = new TGraph();
    TGraph* g_skew_holes = new TGraph();

    for (int d = 0; d < nDiscriminators; ++d) {
        c->cd(d + 1);
        TH1F* h_elect = new TH1F(Form("h_disc_e_%d", d), 
                                 Form("MODULE_NAME Raw Discriminator Distributions %d;LSB;Counts", d), 
                                 100, 0, 500);
        TH1F* h_holes = new TH1F(Form("h_disc_h_%d", d), 
                                 Form("MODULE_NAME Raw Discriminator Distributions %d;LSB;Counts", d), 
                                 100, 0, 500);

        // --- skewness variables ---
        double mean_e=0, m2_e=0, m3_e=0; int count_e=0;
        double mean_h=0, m2_h=0, m3_h=0; int count_h=0;

        for (int hw = 0; hw < nASICs; ++hw) {
            bool skip_e = skip_elect.count(hw);
            bool skip_h = skip_holes.count(hw);
            
            for (int ch = 0; ch < nChannels; ++ch) {
                if (!skip_e) {
                    float ve = adc_data_elect[hw][ch][d];
                    if (ve != -1) {
                        h_elect->Fill(ve);

                        // --- update skewness on-the-fly ---
                        count_e++;
                        double delta = ve - mean_e;
                        mean_e += delta / count_e;
                        m2_e += delta * (ve - mean_e);
                        m3_e += delta * delta * (ve - mean_e);
                    }
                }
                if (!skip_h) {
                    float vh = adc_data_holes[hw][ch][d];
                    if (vh != -1) {
                        h_holes->Fill(vh);

                        // --- update skewness on-the-fly ---
                        count_h++;
                        double delta = vh - mean_h;
                        mean_h += delta / count_h;
                        m2_h += delta * (vh - mean_h);
                        m3_h += delta * delta * (vh - mean_h);
                    }
                }
            }
        }

        // --- compute skewness ---
        double skew_e = (count_e > 1) ? (sqrt(count_e) * m3_e) / pow(m2_e, 1.5) : 0;
        double skew_h = (count_h > 1) ? (sqrt(count_h) * m3_h) / pow(m2_h, 1.5) : 0;

        std::cout << "Discriminator " << d 
                  << " | n-side skew: " << skew_e 
                  << " | p-side skew: " << skew_h << std::endl;

        g_skew_elect->SetPoint(d, d, skew_e);
        g_skew_holes->SetPoint(d, d, skew_h);

        h_elect->SetLineColor(kRed);
        h_holes->SetLineColor(kBlue);
        h_elect->SetStats(kFALSE);
        h_holes->SetStats(kFALSE);
        h_elect->Draw();
        h_holes->Draw("SAME");
    }

    // --- save canvas with distributions ---
    c->SaveAs("MODULE_NAME_discriminator_distributions.root");
    c->SaveAs("MODULE_NAME_discriminator_distributions.pdf");

    // --- plot skewness TGraph ---
    TCanvas* c_skew = new TCanvas("skewness_canvas","Skewness per Discriminator",1000,600);
    g_skew_elect->SetMarkerStyle(20);
    g_skew_elect->SetMarkerColor(kRed);

    g_skew_holes->SetMarkerStyle(21);
    g_skew_holes->SetMarkerColor(kBlue);

    g_skew_elect->SetTitle("Skewness per Discriminator;Discriminator Index;Skewness");
    g_skew_elect->Draw("AP");
    g_skew_holes->Draw("P SAME");

    // --- zero reference line ---
    TLine *line = new TLine(0, 0, nDiscriminators - 1, 0);
    line->SetLineStyle(2);
    line->SetLineColor(kBlack);
    line->Draw();

    // --- fix Y axis range ---
    double minY = std::min(g_skew_elect->GetY()[TMath::LocMin(g_skew_elect->GetN(), g_skew_elect->GetY())],
                           g_skew_holes->GetY()[TMath::LocMin(g_skew_holes->GetN(), g_skew_holes->GetY())]);
    double maxY = std::max(g_skew_elect->GetY()[TMath::LocMax(g_skew_elect->GetN(), g_skew_elect->GetY())],
                           g_skew_holes->GetY()[TMath::LocMax(g_skew_holes->GetN(), g_skew_holes->GetY())]);
    gPad->Update();
    gPad->GetUymin(); // Ensure pad is updated before range setting
    gPad->GetUxmax(); // Just in case

    g_skew_elect->GetYaxis()->SetRangeUser(minY - 0.5, maxY + 0.5);

    // --- legend ---
    TLegend *leg = new TLegend(0.15, 0.85, 0.35, 0.95);
    leg->AddEntry(g_skew_elect, "n-side (electrons)", "p");
    leg->AddEntry(g_skew_holes, "p-side (holes)", "p");
    leg->Draw();

    c_skew->SaveAs("MODULE_NAME_skewness_per_discriminator.pdf");
}



void uncalibrated_discriminators_per_chn(const ASICVector& adc_data_elect, 
                                         const ASICVector& adc_data_holes, 
                                         const std::set<int>& skip_elect,
                                         const std::set<int>& skip_holes) 
{
    std::vector<int> uncalib_elect(totalChannels, 0);
    std::vector<int> uncalib_holes(totalChannels, 0);

    for (int hw = 0; hw < nASICs; ++hw) {
        bool skip_e = skip_elect.count(hw);
        bool skip_h = skip_holes.count(hw);
        
        for (int ch = 0; ch < nChannels; ++ch) {
            int globalChn = hw * nChannels + ch;

            if (!skip_e) {
                const auto& elect_vals = adc_data_elect[hw][ch];
                for (float val : elect_vals) {
                    if (val == -1) uncalib_elect[globalChn]++;
                }
            }
            
            if (!skip_h) {
                const auto& holes_vals = adc_data_holes[hw][ch];
                for (float val : holes_vals) {
                    if (val == -1) uncalib_holes[globalChn]++;
                }
            }
        }
    }

    TH1F* hElect = new TH1F("hUncalibElect", "MODULE_NAME Uncalibrated Discriminators;Channel;# of -1 values", 
                           totalChannels, -0.5, totalChannels - 0.5);
    TH1F* hHoles = new TH1F("hUncalibHoles", "MODULE_NAME Uncalibrated Discriminators;Channel;# of -1 values", 
                           totalChannels, -0.5, totalChannels - 0.5);

    for (int i = 0; i < totalChannels; ++i) {
        hElect->SetBinContent(i + 1, uncalib_elect[i]);
        hHoles->SetBinContent(i + 1, uncalib_holes[i]);
    }

    TCanvas* c = new TCanvas("cUncalibrated", "MODULE_NAME Uncalibrated Discriminators per Channel", 1200, 600);
    c->SetGrid();
    hElect->SetLineColor(kRed);
    hHoles->SetLineColor(kBlue);
    hElect->SetStats(kFALSE);
    hHoles->SetStats(kFALSE);
    hElect->Draw("HIST");
    hHoles->Draw("HIST SAME");

    TLegend* legend = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend->AddEntry(hElect, "Electrons", "l");
    legend->AddEntry(hHoles, "Holes", "l");
    legend->Draw();

    c->SaveAs("MODULE_NAME_uncalibrated_discriminators_per_chn.pdf");
    TFile* outfile = new TFile("MODULE_NAME_uncalibrated_discriminators_per_chn.root", "RECREATE");
    hElect->Write();
    hHoles->Write();
    c->Write();
    outfile->Close();
}

void plot_missing_per_discriminator(const std::vector<int>& missing_elect,
                                    const std::vector<int>& missing_holes,
                                    const int totalChannels = 1024) 
{
    const int nDisc = 31;
    TH1F* h_missing_elect = new TH1F("h_missing_elect", "MODULE_NAME Missing Values per Discriminator;Discriminator;Avg Missing per Channel", nDisc, -0.5, 30.5);
    TH1F* h_missing_holes = new TH1F("h_missing_holes", "MODULE_NAME Missing Values per Discriminator;Discriminator;Avg Missing per Channel", nDisc, -0.5, 30.5);

    for (int d = 0; d < nDisc; ++d) {
        float p_e = missing_elect[d] / float(totalChannels);
        float p_h = missing_holes[d] / float(totalChannels);

        h_missing_elect->SetBinContent(d + 1, p_e);
        h_missing_holes->SetBinContent(d + 1, p_h);

        float sem_e = std::sqrt(p_e * (1 - p_e) / totalChannels);
        float sem_h = std::sqrt(p_h * (1 - p_h) / totalChannels);

        h_missing_elect->SetBinError(d + 1, sem_e);
        h_missing_holes->SetBinError(d + 1, sem_h);
    }

    TCanvas* c2 = new TCanvas("c2", "MODULE_NAME Missing Values per Discriminator", 1000, 600);
    h_missing_elect->SetLineColor(kRed); h_missing_elect->SetMarkerStyle(20); h_missing_elect->SetMarkerColor(kRed);
    h_missing_holes->SetLineColor(kBlue); h_missing_holes->SetMarkerStyle(24); h_missing_holes->SetMarkerColor(kBlue);
    h_missing_elect->SetStats(kFALSE); h_missing_holes->SetStats(kFALSE);

    h_missing_elect->Draw("E1 P");
    h_missing_holes->Draw("SAME E1 P");

    TLegend* legend2 = new TLegend(0.7, 0.8, 0.9, 0.9);
    legend2->AddEntry(h_missing_elect, "Electrons", "lp");
    legend2->AddEntry(h_missing_holes, "Holes", "lp");
    legend2->Draw();

    TFile* outfile = new TFile("MODULE_NAME_missing_per_discriminator.root", "RECREATE");
    h_missing_elect->Write();
    h_missing_holes->Write();
    c2->Write();
    outfile->Close();
}

void heatmaps_mean_skewness(const ASICVector& adc_data_elect,
                            const ASICVector& adc_data_holes,
                            const std::set<int>& skip_elect_sorted,
                            const std::set<int>& skip_holes_sorted,
                            const std::map<int,int>& pMap,
                            const std::map<int,int>& nMap,
                            double mean_min = -1, double mean_max = -1,
                            double skew_min = -2, double skew_max = 2)
{
    gStyle->SetOptStat(0);
    gStyle->SetPaintTextFormat("4.0f"); // integers, no decimals

    // ---------- Helper: compute skewness ----------
    auto compute_skewness = [](const std::vector<double>& v) -> double {
        if (v.size() < 2) return 0.0;
        double mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
        double m2 = 0.0, m3 = 0.0;
        for (auto x : v) {
            double d = x - mean;
            m2 += d*d;
            m3 += d*d*d;
        }
        m2 /= v.size();
        m3 /= v.size();
        double sd = std::sqrt(m2);
        if (sd == 0) return 0.0;
        return m3 / (sd*sd*sd);
    };

    // ---------- GLOBAL HISTOGRAMS ----------
    TH2F* h_mean_elect = new TH2F("h_mean_elect", "Mean (Electrons);Discriminator;ASIC(sorted)",
                                  nDiscriminators, -0.5, nDiscriminators-0.5,
                                  nASICs, -0.5, nASICs-0.5);
    TH2F* h_skew_elect = new TH2F("h_skew_elect", "Skewness (Electrons);Discriminator;ASIC(sorted)",
                                  nDiscriminators, -0.5, nDiscriminators-0.5,
                                  nASICs, -0.5, nASICs-0.5);
    TH2F* h_mean_holes = new TH2F("h_mean_holes", "Mean (Holes);Discriminator;ASIC(sorted)",
                                  nDiscriminators, -0.5, nDiscriminators-0.5,
                                  nASICs, -0.5, nASICs-0.5);
    TH2F* h_skew_holes = new TH2F("h_skew_holes", "Skewness (Holes);Discriminator;ASIC(sorted)",
                                  nDiscriminators, -0.5, nDiscriminators-0.5,
                                  nASICs, -0.5, nASICs-0.5);

    for (int sorted_hw = 0; sorted_hw < nASICs; ++sorted_hw) {

        // ----- ELECTRONS -----
        int hw_p = -1;
        for (const auto& [hw, idx] : pMap) { if (idx == sorted_hw) { hw_p = hw; break; } }
        if (hw_p != -1 && !skip_elect_sorted.count(sorted_hw)) {
            for (int d = 0; d < nDiscriminators; ++d) {
                std::vector<double> vals;
                for (int ch = 0; ch < nChannels; ++ch) {
                    float v = adc_data_elect[sorted_hw][ch][d];
                    if (v != -1) vals.push_back(v);
                }
                if (!vals.empty()) {
                    double mean = std::accumulate(vals.begin(), vals.end(), 0.0)/vals.size();
                    h_mean_elect->SetBinContent(d+1, sorted_hw+1, mean);
                    h_skew_elect->SetBinContent(d+1, sorted_hw+1, compute_skewness(vals));
                }
            }
        }

        // ----- HOLES -----
        int hw_n = -1;
        for (const auto& [hw, idx] : nMap) { if (idx == sorted_hw) { hw_n = hw; break; } }
        if (hw_n != -1 && !skip_holes_sorted.count(sorted_hw)) {
            for (int d = 0; d < nDiscriminators; ++d) {
                std::vector<double> vals;
                for (int ch = 0; ch < nChannels; ++ch) {
                    float v = adc_data_holes[sorted_hw][ch][d];
                    if (v != -1) vals.push_back(v);
                }
                if (!vals.empty()) {
                    double mean = std::accumulate(vals.begin(), vals.end(), 0.0)/vals.size();
                    h_mean_holes->SetBinContent(d+1, sorted_hw+1, mean);
                    h_skew_holes->SetBinContent(d+1, sorted_hw+1, compute_skewness(vals));
                }
            }
        }
    }

    if (mean_min >= 0 && mean_max > mean_min) {
        h_mean_elect->SetMinimum(mean_min); h_mean_elect->SetMaximum(mean_max);
        h_mean_holes->SetMinimum(mean_min); h_mean_holes->SetMaximum(mean_max);
    }
    h_skew_elect->SetMinimum(skew_min); h_skew_elect->SetMaximum(skew_max);
    h_skew_holes->SetMinimum(skew_min); h_skew_holes->SetMaximum(skew_max);

    auto draw_save = [&](TH2F* h, const char* fname){
        TCanvas* c = new TCanvas(Form("c_%s", h->GetName()), h->GetTitle(), 900, 600);
        h->Draw("COLZ TEXT");
        h->SetMarkerSize(1.0);
        c->Update();
        c->SaveAs(fname);
    };

    draw_save(h_mean_elect, "heatmap_mean_electrons_sorted.pdf");
    draw_save(h_mean_holes, "heatmap_mean_holes_sorted.pdf");
    draw_save(h_skew_elect, "heatmap_skew_electrons_sorted.pdf");
    draw_save(h_skew_holes, "heatmap_skew_holes_sorted.pdf");

    // ---------- NEW: PER-ASIC MEAN HEATMAPS IN ONE CANVAS ----------
    auto make_mean_per_side = [&](const ASICVector& data,
                                  const std::set<int>& skip,
                                  const std::map<int,int>& hwMap,
                                  const char* side)
    {
        TCanvas* c = new TCanvas(Form("c_mean_%s_all", side),
                                 Form("%s (all ASICs)", side), 2400, 1200);
        c->Divide(4,2); // 8 ASICs in 2 rows, 4 cols

        int pad = 1;
        for (int sorted_hw = 0; sorted_hw < nASICs; ++sorted_hw) {
            int hw = -1;
            for (const auto& [raw, idx] : hwMap) {
                if (idx == sorted_hw) { hw = raw; break; }
            }
            if (hw == -1 || skip.count(sorted_hw)) continue;

            TH2F* h_mean = new TH2F(Form("h_mean_%s_%d", side, sorted_hw),
                                    Form("%s ASIC %d;Channel;Discriminator", side, sorted_hw),
                                    nChannels, -0.5, nChannels-0.5,
                                    nDiscriminators, -0.5, nDiscriminators-0.5);

            for (int ch = 0; ch < nChannels; ++ch) {
                for (int d = 0; d < nDiscriminators; ++d) {
                    float v = data[sorted_hw][ch][d];
                    if (v == -1) continue;
                    h_mean->SetBinContent(ch+1, d+1, v);
                }
            }

            if (mean_min >= 0 && mean_max > mean_min) {
                h_mean->SetMinimum(mean_min);
                h_mean->SetMaximum(mean_max);
            }

            c->cd(pad++);
            h_mean->Draw("COLZ");
        }

        c->SaveAs(Form("heatmaps_mean_%s_allASICs.pdf", side));
    };

    make_mean_per_side(adc_data_elect, skip_elect_sorted, pMap, "elect");
    make_mean_per_side(adc_data_holes, skip_holes_sorted, nMap, "holes");
}

void channels_per_asic(const ASICVector& adc_data_elect,
                       const ASICVector& adc_data_holes,
                       const std::set<int>& skip_elect,
                       const std::set<int>& skip_holes,
                       int ch_start) 
{
    const int nDisc = nDiscriminators; // 31 discriminators
    const int nChn = 4;                // number of consecutive channels to plot
    const int colors[4] = {kRed, kBlue, kGreen+2, kMagenta};

    // =======================
    // Electrons
    // =======================
    TCanvas* cElect = new TCanvas("cElect", "Electrons - 4 channels per ASIC", 1600, 1000);
    cElect->Divide(4,2); // 8 ASICs

    for (int asic = 0; asic < nASICs; ++asic) {
        if (skip_elect.count(asic)) continue;
        cElect->cd(asic+1);

        TMultiGraph* mg = new TMultiGraph();
        TLegend* leg = new TLegend(0.65, 0.70, 0.90, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);

        for (int ch = 0; ch < nChn; ++ch) {
            int channel = ch_start + ch;
            std::vector<double> x, y;
            for (int d = 0; d < nDisc; ++d) {
                float val = adc_data_elect[asic][channel][d];
                if (val < 0) continue;
                x.push_back(d);
                y.push_back(val);
            }
            if (x.empty()) continue;

            TGraph* gr = new TGraph(x.size(), x.data(), y.data());
            gr->SetLineColor(colors[ch]);
            gr->SetMarkerColor(colors[ch]);
            gr->SetMarkerStyle(20+ch);
            gr->SetLineWidth(2);

            mg->Add(gr, "LP");
            leg->AddEntry(gr, Form("Channel %d", channel), "lp");
        }

        mg->SetTitle(Form("Electron ASIC %d;Discriminator;ADC Value", asic));
        mg->Draw("A");
        leg->Draw();
    }

    // Save electron canvas
    cElect->Update();
    cElect->SaveAs(Form("channels_perASIC_electrons_chStart%d.png", ch_start));
    cElect->SaveAs(Form("channels_perASIC_electrons_chStart%d.root", ch_start));


    // =======================
    // Holes
    // =======================
    TCanvas* cHoles = new TCanvas("cHoles", "Holes - 4 channels per ASIC", 1600, 1000);
    cHoles->Divide(4,2); // 8 ASICs

    for (int asic = 0; asic < nASICs; ++asic) {
        if (skip_holes.count(asic)) continue;
        cHoles->cd(asic+1);

        TMultiGraph* mg = new TMultiGraph();
        TLegend* leg = new TLegend(0.65, 0.70, 0.90, 0.90);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);

        for (int ch = 0; ch < nChn; ++ch) {
            int channel = ch_start + ch;
            std::vector<double> x, y;
            for (int d = 0; d < nDisc; ++d) {
                float val = adc_data_holes[asic][channel][d];
                if (val < 0) continue;
                x.push_back(d);
                y.push_back(val);
            }
            if (x.empty()) continue;

            TGraph* gr = new TGraph(x.size(), x.data(), y.data());
            gr->SetLineColor(colors[ch]);
            gr->SetMarkerColor(colors[ch]);
            gr->SetMarkerStyle(20+ch);
            gr->SetLineWidth(2);

            mg->Add(gr, "LP");
            leg->AddEntry(gr, Form("Channel %d", channel), "lp");
        }

        mg->SetTitle(Form("Hole ASIC %d;Discriminator;ADC Value", asic));
        mg->Draw("A");
        leg->Draw();
    }

    // Save hole canvas
    cHoles->Update();
    cHoles->SaveAs(Form("channels_perASIC_holes_chStart%d.png", ch_start));
    cHoles->SaveAs(Form("channels_perASIC_holes_chStart%d.root", ch_start));
}


void plot_groups(const ASICVector& adc_data,
                             const std::string& side,
                             const std::map<int,int>& sortMap) {
    int nASICs = adc_data.size();
    int nDiscriminators = adc_data[0][0].size();

    TCanvas* c = new TCanvas(Form("c_%s_groups", side.c_str()),
                             Form("%s side: ASIC group means", side.c_str()), 1600, 1200);
    c->Divide(4,2); // 8 ASICs

    int colors[4] = {kRed, kBlue, kGreen+2, kMagenta};

    for (int hw_sorted = 0; hw_sorted < nASICs; ++hw_sorted) {
        // map sorted index -> HW ID
        int hw = -1;
        for (const auto& kv : sortMap) {
            if (kv.second == hw_sorted) { hw = kv.first; break; }
        }
        if (hw == -1) continue;

        c->cd(hw_sorted+1);
        TMultiGraph* mg = new TMultiGraph();
        TLegend* leg = new TLegend(0.65,0.7,0.9,0.9);

        // 4 groups
        for (int g = 0; g < 4; ++g) {
            std::vector<double> x, y, ey;

            for (int d = 0; d < nDiscriminators; ++d) {
                double sum = 0, sum2 = 0;
                int nchan = 0;
                for (int ch = g; ch < 128; ch += 4) {
                    double val = adc_data[hw][ch][d];
                    if (val >= 0) { // skip missing -1
                        sum += val;
                        sum2 += val*val;
                        nchan++;
                    }
                }
                if (nchan>0) {
                    double mean = sum / nchan;
                    double sem = sqrt((sum2 - sum*sum/nchan)/(nchan-1)) / sqrt(nchan); // SEM
                    x.push_back(d);
                    y.push_back(mean);
                    ey.push_back(sem);
                }
            }

            TGraphErrors* gr = new TGraphErrors(x.size(), &x[0], &y[0], nullptr, &ey[0]);
            gr->SetLineColor(colors[g]);
            gr->SetLineWidth(2);
            gr->SetMarkerStyle(20);
            gr->SetMarkerColor(colors[g]);
            mg->Add(gr);
            leg->AddEntry(gr, Form("Group %d", g), "lp");
        }

        mg->SetTitle(Form("ASIC %d", hw));
        mg->Draw("APL");
        mg->GetXaxis()->SetTitle("Discriminator index");
        mg->GetYaxis()->SetTitle("Mean ADC ± SEM");
        leg->Draw();
    }

    c->Update();
    c->SaveAs(Form("groups_%s_allASICs.pdf", side.c_str()));
    c->SaveAs(Form("groups_%s_allASICs.root", side.c_str()));
}

void plot_groups_allASICs(const ASICVector& adc_elect, 
                          const ASICVector& adc_holes,
                          const std::map<int,int>& pMap,
                          const std::map<int,int>& nMap) {
    int nASICs = adc_elect.size();
    int nDiscriminators = adc_elect[0][0].size();

    TCanvas* c = new TCanvas("c_groups_allASICs", "4-Channel Groups: Electrons & Holes", 1600, 1000);
    TMultiGraph* mg = new TMultiGraph();
    TLegend* leg = new TLegend(0.65,0.7,0.9,0.9);

    int colors[4] = {kRed, kBlue, kGreen+2, kMagenta};

    // Loop over 4 groups
    for(int g=0; g<4; ++g) {
        // Electrons
        {
            std::vector<double> x, y, ey;
            for(int d=0; d<nDiscriminators; ++d) {
                std::vector<double> vals;
                for(int hw=0; hw<nASICs; ++hw) {
                    int hw_sorted = nMap.at(hw);
                    for(int ch = g; ch < 128; ch += 4) {
                        double val = adc_elect[hw_sorted][ch][d];
                        if(val >= 0) vals.push_back(val);
                    }
                }
                if(vals.size() > 0) {
                    double sum=0, sum2=0, N = vals.size();
                    for(double v : vals) { sum += v; sum2 += v*v; }
                    double mean = sum/N;
                    double sem = sqrt((sum2 - sum*sum/N)/(N-1)) / sqrt(N);
                    x.push_back(d); y.push_back(mean); ey.push_back(sem);
                }
            }
            TGraphErrors* gr = new TGraphErrors(x.size(), &x[0], &y[0], nullptr, &ey[0]);
            gr->SetLineColor(colors[g]);
            gr->SetMarkerColor(colors[g]);
            gr->SetMarkerStyle(20);
            gr->SetLineWidth(2);
            mg->Add(gr);
            leg->AddEntry(gr, Form("Electrons Group %d", g), "lp");
        }

        // Holes
        {
            std::vector<double> x, y, ey;
            for(int d=0; d<nDiscriminators; ++d) {
                std::vector<double> vals;
                for(int hw=0; hw<nASICs; ++hw) {
                    int hw_sorted = pMap.at(hw);
                    for(int ch = g; ch < 128; ch += 4) {
                        double val = adc_holes[hw_sorted][ch][d];
                        if(val >= 0) vals.push_back(val);
                    }
                }
                if(vals.size() > 0) {
                    double sum=0, sum2=0, N = vals.size();
                    for(double v : vals) { sum += v; sum2 += v*v; }
                    double mean = sum/N;
                    double sem = sqrt((sum2 - sum*sum/N)/(N-1)) / sqrt(N);
                    x.push_back(d); y.push_back(mean); ey.push_back(sem);
                }
            }
            TGraphErrors* gr = new TGraphErrors(x.size(), &x[0], &y[0], nullptr, &ey[0]);
            gr->SetLineColor(colors[g]);
            gr->SetMarkerColor(colors[g]);
            gr->SetMarkerStyle(24); // different marker for holes
            gr->SetLineStyle(2);    // dashed for holes
            gr->SetLineWidth(2);
            mg->Add(gr);
            leg->AddEntry(gr, Form("Holes Group %d", g), "lp");
        }
    }

    mg->Draw("APL");
    mg->GetXaxis()->SetTitle("Discriminator index");
    mg->GetYaxis()->SetTitle("Mean ADC ± SEM");
    mg->SetTitle("4-Channel Groups: Electrons & Holes (All ASICs)");
    leg->Draw();
    c->Update();
    c->SaveAs("groups_allASICs_electrons_vs_holes.pdf");
    c->SaveAs("groups_allASICs_electrons_vs_holes.root");
}

//---------------------------CORRELATION-----------------------------------------------------------

void correlation_map_channels(const ASICVector& adc_elect,
                                      const ASICVector& adc_holes,
                                      const std::map<int,int>& nMap,
                                      const std::map<int,int>& pMap) {
    int nASICs = adc_elect.size();
    int nChannels = 128;
    int nDiscriminators = adc_elect[0][0].size();

    // Canvas for electrons
    TCanvas* c_e = new TCanvas("c_corr_channels_elect", "Channel Correlation Electrons", 1600, 1200);
    c_e->Divide(4,2); // 8 pads: 4 columns × 2 rows
    for(int hw_sorted=0; hw_sorted<nASICs; ++hw_sorted){
        int hw = -1;
        for(const auto& kv : nMap) if(kv.second==hw_sorted){ hw=kv.first; break; }

        c_e->cd(hw_sorted+1);
        TH2F* h = new TH2F(Form("h_corr_ch_elec_%d", hw), Form("Electrons ASIC %d", hw),
                           nChannels,0,nChannels, nChannels,0,nChannels);

        for(int i=0;i<nChannels;i++){
            for(int j=0;j<nChannels;j++){
                std::vector<double> vi,vj;
                for(int d=0;d<nDiscriminators;d++){
                    double val_i = adc_elect[hw][i][d];
                    double val_j = adc_elect[hw][j][d];
                    if(val_i>=0 && val_j>=0){ vi.push_back(val_i); vj.push_back(val_j);}
                }
                double corr=0;
                if(vi.size()>1){
                    double mean_i=0, mean_j=0;
                    for(size_t k=0;k<vi.size();k++){ mean_i+=vi[k]; mean_j+=vj[k]; }
                    mean_i/=vi.size(); mean_j/=vi.size();
                    double num=0, denom_i=0, denom_j=0;
                    for(size_t k=0;k<vi.size();k++){
                        num += (vi[k]-mean_i)*(vj[k]-mean_j);
                        denom_i += (vi[k]-mean_i)*(vi[k]-mean_i);
                        denom_j += (vj[k]-mean_j)*(vj[k]-mean_j);
                    }
                    corr = num / sqrt(denom_i*denom_j);
                }
                h->SetBinContent(i+1,j+1,corr);
            }
        }
        h->Draw("COLZ");
        c_e->SaveAs("correlation_channels_electrons_allASICs.pdf");
	c_e->SaveAs("correlation_channels_electrons_allASICs.root");
    }
    c_e->Update();

    // Canvas for holes
    TCanvas* c_h = new TCanvas("c_corr_channels_holes", "Channel Correlation Holes", 1600, 1200);
    c_h->Divide(4,2); // 8 pads
    for(int hw_sorted=0; hw_sorted<nASICs; ++hw_sorted){
        int hw = -1;
        for(const auto& kv : pMap) if(kv.second==hw_sorted){ hw=kv.first; break; }

        c_h->cd(hw_sorted+1);
        TH2F* h = new TH2F(Form("h_corr_ch_holes_%d", hw), Form("Holes ASIC %d", hw),
                           nChannels,0,nChannels, nChannels,0,nChannels);

        for(int i=0;i<nChannels;i++){
            for(int j=0;j<nChannels;j++){
                std::vector<double> vi,vj;
                for(int d=0;d<nDiscriminators;d++){
                    double val_i = adc_holes[hw][i][d];
                    double val_j = adc_holes[hw][j][d];
                    if(val_i>=0 && val_j>=0){ vi.push_back(val_i); vj.push_back(val_j);}
                }
                double corr=0;
                if(vi.size()>1){
                    double mean_i=0, mean_j=0;
                    for(size_t k=0;k<vi.size();k++){ mean_i+=vi[k]; mean_j+=vj[k]; }
                    mean_i/=vi.size(); mean_j/=vi.size();
                    double num=0, denom_i=0, denom_j=0;
                    for(size_t k=0;k<vi.size();k++){
                        num += (vi[k]-mean_i)*(vj[k]-mean_j);
                        denom_i += (vi[k]-mean_i)*(vi[k]-mean_i);
                        denom_j += (vj[k]-mean_j)*(vj[k]-mean_j);
                    }
                    corr = num / sqrt(denom_i*denom_j);
                }
                h->SetBinContent(i+1,j+1,corr);
            }
        }
        h->Draw("COLZ");
        c_h->SaveAs("correlation_channels_holes_allASICs.pdf");
	c_h->SaveAs("correlation_channels_holes_allASICs.root");
    }
    c_h->Update();
}


void correlation_map_ASICs(const ASICVector& adc_elect,
                                    const ASICVector& adc_holes,
                                    const std::map<int,int>& nMap,
                                    const std::map<int,int>& pMap) {
    int nASICs = adc_elect.size();
    int nDiscriminators = adc_elect[0][0].size();

    TCanvas* c = new TCanvas("c_corr_ASICs_combined", "ASIC Correlation: Electrons & Holes", 1200, 600);
    c->Divide(2,1); // 2 pads: left=Electrons, right=Holes

    // Electrons
    c->cd(1);
    TH2F* h_e = new TH2F("h_corr_ASICs_e", "Electrons ASIC Correlation;ASIC i;ASIC j", nASICs,0,nASICs, nASICs,0,nASICs);
    for(int i=0;i<nASICs;i++){
        int hw_i = -1; for(const auto& kv: nMap) if(kv.second==i){hw_i=kv.first; break;}
        for(int j=0;j<nASICs;j++){
            int hw_j = -1; for(const auto& kv: nMap) if(kv.second==j){hw_j=kv.first; break;}
            std::vector<double> vi,vj;
            for(int ch=0;ch<128;ch++) for(int d=0;d<nDiscriminators;d++){
                double val_i = adc_elect[hw_i][ch][d];
                double val_j = adc_elect[hw_j][ch][d];
                if(val_i>=0 && val_j>=0){ vi.push_back(val_i); vj.push_back(val_j);}
            }
            double corr=0;
            if(vi.size()>1){
                double mean_i=0, mean_j=0; for(size_t k=0;k<vi.size();k++){mean_i+=vi[k]; mean_j+=vj[k];} mean_i/=vi.size(); mean_j/=vi.size();
                double num=0, denom_i=0, denom_j=0; for(size_t k=0;k<vi.size();k++){num+=(vi[k]-mean_i)*(vj[k]-mean_j); denom_i+=(vi[k]-mean_i)*(vi[k]-mean_i); denom_j+=(vj[k]-mean_j)*(vj[k]-mean_j);}
                corr=num/sqrt(denom_i*denom_j);
            }
            h_e->SetBinContent(i+1,j+1,corr);
        }

    }
    h_e->Draw("COLZ");
    h_e->SaveAs("correlation_ASICs_electrons.pdf");
    h_e->SaveAs("correlation_ASICs_electrons.root");

    // Holes
    c->cd(2);
    TH2F* h_h = new TH2F("h_corr_ASICs_h", "Holes ASIC Correlation;ASIC i;ASIC j", nASICs,0,nASICs, nASICs,0,nASICs);
    for(int i=0;i<nASICs;i++){
        int hw_i=-1; for(const auto& kv: pMap) if(kv.second==i){hw_i=kv.first; break;}
        for(int j=0;j<nASICs;j++){
            int hw_j=-1; for(const auto& kv: pMap) if(kv.second==j){hw_j=kv.first; break;}
            std::vector<double> vi,vj;
            for(int ch=0;ch<128;ch++) for(int d=0;d<nDiscriminators;d++){
                double val_i = adc_holes[hw_i][ch][d];
                double val_j = adc_holes[hw_j][ch][d];
                if(val_i>=0 && val_j>=0){ vi.push_back(val_i); vj.push_back(val_j);}
            }
            double corr=0;
            if(vi.size()>1){
                double mean_i=0, mean_j=0; for(size_t k=0;k<vi.size();k++){mean_i+=vi[k]; mean_j+=vj[k];} mean_i/=vi.size(); mean_j/=vi.size();
                double num=0, denom_i=0, denom_j=0; for(size_t k=0;k<vi.size();k++){num+=(vi[k]-mean_i)*(vj[k]-mean_j); denom_i+=(vi[k]-mean_i)*(vi[k]-mean_i); denom_j+=(vj[k]-mean_j)*(vj[k]-mean_j);}
                corr=num/sqrt(denom_i*denom_j);
            }
            h_h->SetBinContent(i+1,j+1,corr);
        }
    }
    h_h->Draw("COLZ");
    h_h->SaveAs("correlation_ASICs_holes.pdf");
    h_h->SaveAs("correlation_ASICs_holes.root");
    c->Update();   
}

void correlation_map_groups(const ASICVector& adc_elect,
                            const ASICVector& adc_holes,
                            const std::map<int,int>& nMap,
                            const std::map<int,int>& pMap) {

    int nASICs = adc_elect.size();
    int nChannels = 128;
    int nDiscriminators = adc_elect[0][0].size();
    const int nGroups = 4; // 4 groups

    // Define channel indices per group (example: group 0 = 0,3,7,...)
    std::vector<std::vector<int>> groups(nGroups);
    for(int ch=0; ch<nChannels; ch++){
        int g = ch % nGroups;
        groups[g].push_back(ch);
    }

    // --- Electrons ---
    TCanvas* c_e = new TCanvas("c_corr_groups_elect", "Group Correlation Electrons", 1600, 1200);
    c_e->Divide(4,2); // 8 pads
    for(int hw_sorted=0; hw_sorted<nASICs; ++hw_sorted){
        int hw=-1;
        for(const auto& kv:nMap) if(kv.second==hw_sorted){ hw=kv.first; break; }

        c_e->cd(hw_sorted+1);
        TH2F* h = new TH2F(Form("h_corr_grp_elec_%d", hw),
                           Form("Electrons ASIC %d", hw),
                           nGroups,0,nGroups, nGroups,0,nGroups);

        // Compute group correlation
        for(int g1=0; g1<nGroups; g1++){
            for(int g2=0; g2<nGroups; g2++){
                std::vector<double> v1, v2;
                for(int d=0; d<nDiscriminators; d++){
                    double val1=0, val2=0;
                    int count1=0, count2=0;
                    for(int ch : groups[g1]){
                        double v = adc_elect[hw][ch][d];
                        if(v>=0){ val1+=v; count1++; }
                    }
                    for(int ch : groups[g2]){
                        double v = adc_elect[hw][ch][d];
                        if(v>=0){ val2+=v; count2++; }
                    }
                    if(count1>0 && count2>0){ v1.push_back(val1/count1); v2.push_back(val2/count2); }
                }

                double corr=0;
                if(v1.size()>1){
                    double mean1=0, mean2=0;
                    for(size_t k=0;k<v1.size();k++){ mean1+=v1[k]; mean2+=v2[k]; }
                    mean1/=v1.size(); mean2/=v1.size();
                    double num=0, den1=0, den2=0;
                    for(size_t k=0;k<v1.size();k++){
                        num += (v1[k]-mean1)*(v2[k]-mean2);
                        den1 += (v1[k]-mean1)*(v1[k]-mean1);
                        den2 += (v2[k]-mean2)*(v2[k]-mean2);
                    }
                    corr = num / sqrt(den1*den2);
                }

                h->SetBinContent(g1+1,g2+1,corr);
            }
        }
        h->Draw("COLZ");
        c_e->SaveAs(Form("correlation_groups_electrons_ASIC%d.pdf", hw));
	c_e->SaveAs(Form("correlation_groups_electrons_ASIC%d.root", hw));
    }
    c_e->Update();

    // --- Holes (same logic) ---
    TCanvas* c_h = new TCanvas("c_corr_groups_holes", "Group Correlation Holes", 1600, 1200);
    c_h->Divide(4,2);
    for(int hw_sorted=0; hw_sorted<nASICs; ++hw_sorted){
        int hw=-1;
        for(const auto& kv:pMap) if(kv.second==hw_sorted){ hw=kv.first; break; }

        c_h->cd(hw_sorted+1);
        TH2F* h = new TH2F(Form("h_corr_grp_holes_%d", hw),
                           Form("Holes ASIC %d", hw),
                           nGroups,0,nGroups, nGroups,0,nGroups);

        for(int g1=0; g1<nGroups; g1++){
            for(int g2=0; g2<nGroups; g2++){
                std::vector<double> v1, v2;
                for(int d=0; d<nDiscriminators; d++){
                    double val1=0, val2=0;
                    int count1=0, count2=0;
                    for(int ch : groups[g1]){
                        double v = adc_holes[hw][ch][d];
                        if(v>=0){ val1+=v; count1++; }
                    }
                    for(int ch : groups[g2]){
                        double v = adc_holes[hw][ch][d];
                        if(v>=0){ val2+=v; count2++; }
                    }
                    if(count1>0 && count2>0){ v1.push_back(val1/count1); v2.push_back(val2/count2); }
                }

                double corr=0;
                if(v1.size()>1){
                    double mean1=0, mean2=0;
                    for(size_t k=0;k<v1.size();k++){ mean1+=v1[k]; mean2+=v2[k]; }
                    mean1/=v1.size(); mean2/=v1.size();
                    double num=0, den1=0, den2=0;
                    for(size_t k=0;k<v1.size();k++){
                        num += (v1[k]-mean1)*(v2[k]-mean2);
                        den1 += (v1[k]-mean1)*(v1[k]-mean1);
                        den2 += (v2[k]-mean2)*(v2[k]-mean2);
                    }
                    corr = num / sqrt(den1*den2);
                }

                h->SetBinContent(g1+1,g2+1,corr);
            }
        }
        h->Draw("COLZ");
        c_h->SaveAs(Form("correlation_groups_holes_ASIC%d.pdf", hw));
	c_h->SaveAs(Form("correlation_groups_holes_ASIC%d.root", hw));
    }
    c_h->Update();
       
}

void group_discriminator_correlation(const ASICVector& adc_data,
                                     const std::map<int,int>& hwMap,
                                     const std::string& name_prefix = "electrons") {
    int nASICs = adc_data.size();
    int nChannels = 128;
    int nDiscriminators = adc_data[0][0].size();
    const int nGroups = 4;

    // Define channel groups
    std::vector<std::vector<int>> groups(nGroups);
    for(int ch=0; ch<nChannels; ch++) groups[ch % nGroups].push_back(ch);

    for(int hw_sorted=0; hw_sorted<nASICs; ++hw_sorted){
        int hw=-1;
        for(const auto& kv : hwMap) if(kv.second==hw_sorted){ hw=kv.first; break; }

        // One canvas per ASIC
        TCanvas* c = new TCanvas(Form("c_grp_disc_corr_%s_%d", name_prefix.c_str(), hw),
                                 Form("Group vs Discriminator Correlation %s ASIC %d", name_prefix.c_str(), hw),
                                 1600, 1200);
        c->Divide(2,2); // 4 groups

        for(int g=0; g<nGroups; g++){
            c->cd(g+1);
            TH2F* h = new TH2F(Form("h_grp_disc_corr_%s_%d_g%d", name_prefix.c_str(), hw, g),
                               Form("%s ASIC %d Group %d", name_prefix.c_str(), hw, g),
                               nDiscriminators,0,nDiscriminators,
                               nDiscriminators,0,nDiscriminators);

            // Build vectors for each discriminator in this group
            std::vector<std::vector<double>> adc_vectors(nDiscriminators);
            for(int d=0; d<nDiscriminators; d++){
                for(int ch : groups[g]){
                    double v = adc_data[hw][ch][d];
                    if(v>=0) adc_vectors[d].push_back(v);
                }
            }

            // Compute correlation matrix
            for(int d1=0; d1<nDiscriminators; d1++){
                for(int d2=0; d2<nDiscriminators; d2++){
                    double corr = 0;
                    const auto &v1 = adc_vectors[d1];
                    const auto &v2 = adc_vectors[d2];

                    if(v1.size()>1 && v2.size()>1){
                        double mean1=0, mean2=0;
                        for(double val : v1) mean1+=val;
                        for(double val : v2) mean2+=val;
                        mean1/=v1.size(); mean2/=v2.size();

                        double num=0, den1=0, den2=0;
                        size_t n = std::min(v1.size(), v2.size());
                        for(size_t k=0;k<n;k++){
                            num += (v1[k]-mean1)*(v2[k]-mean2);
                            den1 += (v1[k]-mean1)*(v1[k]-mean1);
                            den2 += (v2[k]-mean2)*(v2[k]-mean2);
                        }
                        if(den1>0 && den2>0) corr = num / sqrt(den1*den2);
                    }
                    h->SetBinContent(d1+1,d2+1,corr);
                }
            }
            h->Draw("COLZ");
        }
        c->Update();
                // Save results
        c->SaveAs(Form("group_disc_corr_%s_ASIC%d.pdf", name_prefix.c_str(), hw));
        c->SaveAs(Form("group_disc_corr_%s_ASIC%d.root", name_prefix.c_str(), hw));
    }
}


void trim_QA() {
    // 1) Allocate raw data vectors
    ASICVector adc_data_elect(nASICs, ChannelVector(nChannels, DiscriminatorVector(nDiscriminators, -1.f)));
    ASICVector adc_data_holes(nASICs, ChannelVector(nChannels, DiscriminatorVector(nDiscriminators, -1.f)));

    // 2) Store header information
    std::map<std::string, FileHeaderInfo> header_info;

    // 3) Ask user for FEB side
    std::string feb_side;
    std::cout << "Enter p-side FEB (A or B): ";
    std::cin >> feb_side;

    // 4) Prepare raw skip sets (using HW IDs)
    std::set<int> skip_elect_hw;
    std::set<int> skip_holes_hw;

    // 5) Load ADC data
    load_adc_data(adc_data_elect, adc_data_holes, header_info, skip_elect_hw, skip_holes_hw, "./", feb_side);

    std::cout << "Skipped electron ASICs (HW IDs): ";
    for (int hw : skip_elect_hw) std::cout << hw << " ";
    std::cout << "\nSkipped hole ASICs (HW IDs): ";
    for (int hw : skip_holes_hw) std::cout << hw << " ";
    std::cout << std::endl;

    // 6) Decide mapping for p-side and n-side
    const auto& nMap = (feb_side == "A") ? sortMapFEB_B : sortMapFEB_A;
    const auto& pMap = (feb_side == "A") ? sortMapFEB_A : sortMapFEB_B;

    // 7) Reorder ADC data
    ASICVector adc_data_elect_sorted = reorder_asics(adc_data_elect, nMap);
    ASICVector adc_data_holes_sorted = reorder_asics(adc_data_holes, pMap);

    // 8) Convert skip sets from HW IDs to sorted indices
    std::set<int> skip_elect_sorted;
    std::set<int> skip_holes_sorted;
    
    for (int hw : skip_elect_hw) skip_elect_sorted.insert(nMap.at(hw));
    for (int hw : skip_holes_hw) skip_holes_sorted.insert(pMap.at(hw));

    // 9) Call existing analysis functions
    MODULE_NAME_statistics(adc_data_elect_sorted, adc_data_holes_sorted, 
                              skip_elect_sorted, skip_holes_sorted, 
                              header_info, feb_side);
    mean_adc_discriminators(adc_data_elect_sorted, adc_data_holes_sorted, skip_elect_sorted, skip_holes_sorted);
    mean_adc_discriminators_odd_even(adc_data_elect_sorted, adc_data_holes_sorted, skip_elect_sorted, skip_holes_sorted);
    discriminator_raw_adc_distributions(adc_data_elect_sorted, adc_data_holes_sorted, skip_elect_sorted, skip_holes_sorted);
    uncalibrated_discriminators_per_chn(adc_data_elect_sorted, adc_data_holes_sorted, skip_elect_sorted, skip_holes_sorted);

    // 10) Compute missing discriminator counts per side
    std::vector<int> missing_elect(nDiscriminators, 0);
    std::vector<int> missing_holes(nDiscriminators, 0);
    for (int hw = 0; hw < nASICs; ++hw) {
        if (skip_elect_sorted.count(hw) && skip_holes_sorted.count(hw)) continue;
        for (int ch = 0; ch < nChannels; ++ch) {
            const auto& ve = adc_data_elect_sorted[hw][ch];
            const auto& vh = adc_data_holes_sorted[hw][ch];
            for (int d = 0; d < nDiscriminators; ++d) {
                if (!skip_elect_sorted.count(hw) && ve[d] == -1) missing_elect[d]++;
                if (!skip_holes_sorted.count(hw) && vh[d] == -1) missing_holes[d]++;
            }
        }
    }
    plot_missing_per_discriminator(missing_elect, missing_holes);
    
heatmaps_mean_skewness(adc_data_elect_sorted, adc_data_holes_sorted,
                       skip_elect_sorted, skip_holes_sorted,
                       pMap, nMap,
                       100, 150,   // mean min/max
                       -1.5, 1.5);     // skew min/max

channels_per_asic(adc_data_elect_sorted, adc_data_holes_sorted,
                       skip_elect_sorted, skip_holes_sorted, 
                       0); // for channels 0,1,2,3
                       
mean_adc_discriminators_perASIC_oddeven(adc_data_elect_sorted, adc_data_holes_sorted, skip_elect_sorted, skip_holes_sorted);
plot_groups(adc_data_elect_sorted, "Electrons", nMap);
plot_groups(adc_data_holes_sorted, "Holes", pMap);
plot_groups_allASICs(adc_data_elect_sorted, adc_data_holes_sorted, pMap, nMap);
correlation_map_channels(adc_data_elect_sorted, adc_data_holes_sorted, nMap, pMap);
correlation_map_ASICs(adc_data_elect_sorted, adc_data_holes_sorted, nMap, pMap);
correlation_map_groups(adc_data_elect_sorted, adc_data_holes_sorted, nMap, pMap);
group_discriminator_correlation(adc_data_elect_sorted, nMap, "electrons");
group_discriminator_correlation(adc_data_holes_sorted, pMap, "holes");
}
