#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPad.h"
#include "TString.h"
#include "TF1.h"
#include "TGraph.h"
#include "TSpline.h"
#include "TStyle.h"
#include <tuple>
#include <vector>
#include <fstream>


double ModBox(double dEdx);
double Recombination(double dEdx);


void recomb_finder(){
   
   /*
   // Get lookup curve
   TFile *f  = new TFile("TEST.root","READ");
   // Get TGraph
   TGraph *g = (TGraph*)f->Get("ESTAR_energy_lookup_curve_lin"); 

   const int n_points = g->GetN();
*/

    //read in estar data
    std::ifstream estar_data("estar_final.txt");
    double energy, CSDA;
    std::vector<double> dEdx_vec;
    std::vector<double> range_vec;
    std::vector<double> energy_vec;
    // calculate dEdx
    while(estar_data >> energy >> CSDA){
        range_vec.push_back(CSDA/1.38434);
        energy_vec.push_back(energy);
        dEdx_vec.push_back(energy/(CSDA/1.38434));
    }

    std::vector<double> dQdx;
    std::vector<double> dQ;
    std::vector<double> recomb;
    for(auto i = 0; i < dEdx_vec.size(); i++){
        dQdx.push_back(ModBox(dEdx_vec[i]));
        //std::cout << (ModBox3(dEdx_vec[i]) * CSDA_vec[i])/(ModBox2(dEdx_vec[i]) * CSDA_vec[i]) << std::endl; 
        recomb.push_back(Recombination(dEdx_vec[i]));
        dQ.push_back(dQdx[i] * range_vec[i]);
        std::cout << "energy:  " << energy_vec[i] << "  dEdx  " << dEdx_vec[i] << "  dQdx  " << dQdx[i] << "   range  " << range_vec[i] <<  "  recomb: " << Recombination(dEdx_vec[i]) <<  "  recomb2:  " << (23.6E-6) * dQdx[i] / dEdx_vec[i] << std::endl;
    }
   
    // line at y = 0.64
    TLine *line_064 = new TLine(0,0.64,10000,0.64);
    line_064->SetLineWidth(2);
    line_064->SetLineColor(kRed);

    // line at x = 
    TLine *xline = new TLine(0.29, 0, 0.29, 10000);
    xline->SetLineWidth(2);
    xline->SetLineColor(kAzure+1);

    TLine *xline2 = new TLine(0.2, 0, 0.2, 10000);
    xline2->SetLineWidth(2);
    xline2->SetLineStyle(9);
    xline2->SetLineColor(kAzure+1);

    TGraph *g_dEdx_recomb   = new TGraph(dEdx_vec.size(), &dEdx_vec[0], &recomb[0]);
    TGraph *g_energy_recomb = new TGraph(energy_vec.size(), &energy_vec[0], &recomb[0]);
    TGraph *g_energy_dEdx   = new TGraph(energy_vec.size(), &energy_vec[0], &dEdx_vec[0]);
    TGraph *g_charge_recomb = new TGraph(dQ.size(), &dQ[0], &recomb[0]);
    TGraph *g_range_recomb  = new TGraph(range_vec.size(), &range_vec[0], &recomb[0]);
    
    TCanvas *c = new TCanvas();
    g_dEdx_recomb->GetXaxis()->SetTitle("dE/dx [MeV/cm]");
    g_dEdx_recomb->GetYaxis()->SetTitle("Recombination");
    g_dEdx_recomb->SetLineWidth(4);
    g_dEdx_recomb->Draw();
    line_064->Draw("same");
    gPad->SetLogx();

    TCanvas *c1 = new TCanvas();
    g_energy_recomb->GetXaxis()->SetTitle("Energy [MeV]");
    g_energy_recomb->GetYaxis()->SetTitle("Recombination");
    g_energy_recomb->SetLineWidth(4);
    g_energy_recomb->Draw();
    line_064->Draw("same");
    xline->Draw("same");
    //xline2->Draw("same");

    TCanvas *c2 = new TCanvas();
    g_charge_recomb->GetXaxis()->SetTitle("Charge [nelectrons]");
    g_charge_recomb->GetYaxis()->SetTitle("Recombination");
    g_charge_recomb->Draw();

    TCanvas *c3 = new TCanvas();
    g_range_recomb->GetXaxis()->SetTitle("Range [cm]");
    g_range_recomb->GetYaxis()->SetTitle("Recombination");
    g_range_recomb->Draw();

    TFile *f = TFile::Open("dedx_recomb.root", "RECREATE");      
    //g_energy_recomb->Write("g_energy_vs_recomb");
    //g_charge_recomb->Write("g_charge_vs_recomb");
    //g_range_recomb->Write("g_range_vs_recomb");
    g_dEdx_recomb->Write("g_dEdx_vs_recomb");

}

// ModBox function to calculate recomobination

double ModBox(const double dEdx){
    double rho   = 1.38434;
    double Alpha = 0.93;
    double Beta  = 0.212/(rho * 0.5);
    double Wion = 23.6E-6;
    double dQdx;
    dQdx = std::log(Alpha + Beta * dEdx)/(Beta * Wion);

    return dQdx; 
}

double Recombination(const double dEdx){
    double rho   = 1.38434;
    double Alpha = 0.93;
    double Beta  = 0.212/(rho * 0.5);
    double R;
    R = std::log(Alpha + Beta * dEdx)/(Beta * dEdx);

    return R; 
}

