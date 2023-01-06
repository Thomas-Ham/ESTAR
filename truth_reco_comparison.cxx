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


int ModBox(double dEdx, double efield);
double ModBox2(double dEdx);
int ModBox3(double dEdx);

void truth_reco_comparison(){
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

    /*
    for (int i =0; i < range_vec.size(); i++){
        std::cout << "range vec: " << range_vec[i] << std::endl;
    }
    */

    // min max Efield
    double max_efield = 0.55;
    double min_efield = 0.45;
    double efield_interval = (max_efield - min_efield)/dEdx_vec.size();
    std::vector<double> efield_vec;
    double efield_temp = min_efield;
    for(auto i = 0; i < dEdx_vec.size(); i++){
        efield_vec.push_back(efield_temp);
        efield_temp += efield_interval;
    }


    // Get corresponding dQdx
    std::vector<double> Q;
    // loop over each value of the efield and get Q for each dE/dx
    for(auto j = efield_vec.begin(); j != efield_vec.end(); j++){;
        for(auto i = 0; i < dEdx_vec.size(); i++){
            Q.push_back(ModBox(dEdx_vec[i], *j) * range_vec[i]);
        }
    }       

    std::vector<double> Q2;
    for(auto i = 0; i < dEdx_vec.size(); i++){
        Q2.push_back(ModBox2(dEdx_vec[i]) * range_vec[i]);
        //std::cout << (ModBox3(dEdx_vec[i]) * CSDA_vec[i])/(ModBox2(dEdx_vec[i]) * CSDA_vec[i]) << std::endl; 
    }

    // If we have N values for the energy and E-Field, we'll have N^2 values for Q
    // since we calculate Q for each combination of energy and E-field.
    // So we need to 'extend' the energy and E-field vector appropriately so the TGraph is made correctly.
    // i.e. The E-field is at the min value N times while we scan over all the energy values,
    // then we go to the next E-field value N times and repeat.
    // Note; This seems like a really stupid way to do things, but it isn't slow so don't care
    std::vector<double> efield_vec2;
    std::vector<double> energy_vec2;
    for(auto i = 0; i < dEdx_vec.size(); i++){
        for(auto j = 0; j < dEdx_vec.size(); j++){
            efield_vec2.push_back(efield_vec[i]);
            energy_vec2.push_back(energy_vec[j]);
        }
    }


    TGraph2D *g = new TGraph2D(efield_vec2.size(), &Q[0], &efield_vec2[0] ,&energy_vec2[0]);
    
    TGraph *g2 = new TGraph(dEdx_vec.size(), &Q2[0], &energy_vec[0]);
     
   
    g->GetHistogram()->GetXaxis()->SetTitle("Collected Charge [No. of electrons]");
    g->GetHistogram()->GetXaxis()->SetTitleOffset(2);
    g->GetHistogram()->GetXaxis()->SetRangeUser(0,150E6);
//    g->GetHistogram()->GetYaxis()->SetRangeUser(0,5);
    g->GetHistogram()->GetYaxis()->SetTitle("E-Field [kV/cm]");
    g->GetHistogram()->GetYaxis()->SetNdivisions(505);
    g->GetHistogram()->GetYaxis()->SetTitleOffset(2);
    g->GetHistogram()->GetZaxis()->SetTitle("Energy [MeV]");
    g->GetHistogram()->GetZaxis()->SetTitleOffset(1.3);
   
//    g->SetNpx(101);
//    g->SetNpy(100);
    //g->SetMarkerStyle(20);
    g->Draw("surf1");


    TF1 *QvsE_fit = new TF1("QvsE_fit","pol1");                                      
    //QvsE_fit->SetParameters(0,0);
    //QvsE_fit->FixParameter(0,0);
    g2->Fit(QvsE_fit); //, "rob=0.95");
     

    g2->SetMarkerStyle(20);
    //g2->GetXaxis()->SetRangeUser(0, 151E3);
    g2->GetYaxis()->SetRangeUser(0, 5);
    //g2->Fit(QvsE_fit);
    g2->Draw("AP");


    TFile *f = TFile::Open("TEST.root", "RECREATE");      
    g->Write("ESTAR_energy_lookup_curve");
    g2->Write("ESTAR_energy_lookup_curve_lin");

    std::cout << "eval: " << g->Interpolate(8.05E6,0.467) << std::endl;
    std::cout << "eval: " << g->Interpolate(8.05E6,0.51) << std::endl;
                                                               
    TCanvas* c2 = new TCanvas();
    g2->Draw();

}

// ModBox function to calculate recomobination
int ModBox(const double dEdx, const double efield){
    double rho   = 1.38;
    double Alpha = 0.93;
    double Beta  = 0.212/(rho * efield);
    double Wion = 23.6E-6;
    double dQdx;
    dQdx = std::log(Alpha + Beta * dEdx)/(Beta * Wion);

    return dQdx; 
}

double ModBox2(const double dEdx){
    double rho   = 1.38434;
    double Alpha = 0.93;
    double Beta  = 0.212/(rho * 0.5);
    double Wion = 23.6E-6;
    double dQdx;
    dQdx = std::log(Alpha + Beta * dEdx)/(Beta * Wion);

    return dQdx; 
}

int ModBox3(const double dEdx){
    double rho   = 1.38;
    double Alpha = 0.93;
    double Beta  = 0.212/(rho * 0.52);
    double Wion = 23.6E-6;
    double dQdx;
    dQdx = std::log(Alpha + Beta * dEdx)/(Beta * Wion);

    return dQdx; 
}

