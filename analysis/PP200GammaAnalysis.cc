
/*
 * Photon spectrum analysis for a single final-state hadrons file.
 *
 * - Reads /scratch/user/cameron.parker/projects/runs/pp-gamma/test_out_final_state_hadrons.dat
 * - Fills photon pT spectrum (pid==22)
 * - Loads Table 7 from ../external_packages/gammaLoss/data/pp2760PhotonSpectrum.root
 * - Plots both spectra for visual comparison
 * 
 * Usage: Just run, no arguments needed.
 * Output: photon_spectrum.root (contains spectrum), comparison.png (overlay plot)
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <dirent.h>
#include <cmath>
#include <memory>
#include <cstdlib>  // for exit()
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"
#include "TAxis.h"
#include "TColor.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TPaveText.h"

// Organization-specific includes (Pythia8 and reader)
#include "JetScapeReaderFinalStateHadrons.h"
#include "JetScapeBanner.h" // if needed
#include "JetScapeParticles.h"
#include "Pythia8/Pythia.h"

// Helper function to get Table 7 TGraphAsymmErrors from a ROOT file
TGraphAsymmErrors* LoadTableGraph(const std::string& rootFilePath, string tableName) {
    TFile *f = TFile::Open(rootFilePath.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: Cannot open ROOT file: " << rootFilePath << std::endl;
        return nullptr;
    }
    TDirectory* tab7 = (TDirectory*)f->Get(tableName.c_str());
    if (!tab7) {
        std::cerr << "ERROR: Table not found in " << rootFilePath << std::endl;
        return nullptr;
    }
    TGraphAsymmErrors* g = (TGraphAsymmErrors*)tab7->Get("Graph1D_y2");
    if (!g) {
        std::cerr << "ERROR: Could not retrieve Graph1D from ROOT file." << std::endl;
        return nullptr;
    }
    return g;
}

void ScaleByCenter(TH1D* hist) {
    // Loop over all bins (ROOT uses bin numbering starting at 1)
    for (int i = 1; i <= hist->GetNbinsX(); i++) {
        double binCenter = hist->GetBinCenter(i);
        // Avoid division by zero.
        if (binCenter != 0) {
            double content = hist->GetBinContent(i);
            double error = hist->GetBinError(i);
            hist->SetBinContent(i, content / binCenter);
            hist->SetBinError(i, error / binCenter);
        }
    }
}

int main(int argc, char** argv){
    if(argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <analysisDir>" << std::endl;
        return 1;
    }
    
    std::string analysisDir = argv[1];

    // --- File paths
    const std::string hadronFile = analysisDir + "test_out_final_state_hadrons.dat";
    const std::string rootDataFile = "../external_packages/gammaLoss/data/pp200PhotonSpectrum.root";
    const std::string outputRootFile = analysisDir + "photon_spectrum.root";

    // Define pT bins to match Table 7 as closely as possible (adjust these to match binning of Table 7 if known)
    double etaCut = 0.35;
    std::vector<double> ptBins = {
        0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0,
        3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 16.0
    };
    TH1D* hPhotonPt = new TH1D("photonPt", "Prompt photon p_{T} spectrum; p_{T} [GeV]; dN/dp_{T}", ptBins.size()-1, ptBins.data());
    TH1D* hDecayPhotonPt = new TH1D("decayphotonPt", "Decay photon p_{T} spectrum; p_{T} [GeV]; dN/dp_{T}", ptBins.size()-1, ptBins.data());

    long nPhotons = 0;
    long nEvents = 0;

    auto reader = make_shared<Jetscape::JetScapeReaderFinalStateHadrons>(hadronFile);
    
    while(!reader->Finished())
    {
        try{
            reader->Next();
        }
        catch(...){
            break;
        }
        nEvents++;
        std::vector<std::shared_ptr<Jetscape::Hadron>> particles = reader->GetHadrons();
        if(particles.empty())
            break;
        
        // Loop through particles in event.
        for(auto& particle : particles)
        {
            if(std::fabs(particle.get()->eta()) > etaCut)
                continue;

            if(particle.get()->pid() == 22)
            {
                nPhotons++;
                if(particle.get()->pstat() == 23 or particle.get()->pstat() == 22)
                    hPhotonPt->Fill(particle.get()->pt());
                    
                if(particle.get()->pstat() == 0)
                    hDecayPhotonPt->Fill(particle.get()->pt());
            }
        } // end particle loop
    } // end event loop
    double sigmaGen = reader->GetSigmaGen();
    if(sigmaGen < 0) sigmaGen = 41.9734;
    cout << "SigmaGen: " << sigmaGen << endl;

    // --- Save histogram to ROOT file
    TFile* outf = TFile::Open(outputRootFile.c_str(), "RECREATE");
    hPhotonPt->Scale(1.0/(2.0*M_PI*2.0*etaCut*nEvents), "width");
    hDecayPhotonPt->Scale(1.0e9*sigmaGen/(2.0*M_PI*2.0*etaCut*nEvents), "width");
    ScaleByCenter(hPhotonPt);;
    ScaleByCenter(hDecayPhotonPt);
    hPhotonPt->Write("dNdpt");
    hPhotonPt->Scale(1.0e9*sigmaGen);
    hPhotonPt->Write("dsigma/dpt");
    hDecayPhotonPt->Write("decay photons");

    std::cout << "Found " << nPhotons << " photons in file.\n";
    std::cout << "Photon pT spectrum histogram written to " << outputRootFile << std::endl;
    
    // --- Load experimental spectrum for comparison
    TGraphAsymmErrors* expGraph = LoadTableGraph(rootDataFile, "Figure 2");
    if (!expGraph) {
        std::cerr << "ERROR: Unable to retrieve experimental Table 2 graph for overlay!" << std::endl;
        return 2;
    }

    // --- Create comparison plot
    TCanvas* c1 = new TCanvas("c1", "Photon spectrum comparison", 800, 600);
    c1->SetLogy();
    c1->SetLogx();

    hPhotonPt->SetLineColor(kBlue+1);
    hPhotonPt->SetMarkerColor(kBlue+1);
    hPhotonPt->SetFillColorAlpha(kBlue+1, 0.35);
    hPhotonPt->SetLineWidth(2);
    auto photonSigmaGraph = new TGraphErrors(hPhotonPt);

    expGraph->SetMarkerStyle(24);
    expGraph->SetLineWidth(2);
    expGraph->Scale(1.0e9);

    auto fa1 = new TF1("parameterization","2*6495*200/((0.8+x)^5)",.1,15); //https://arxiv.org/pdf/hep-ph/0308085 eq 17

    auto mg = new TMultiGraph;
    mg->SetTitle("Photon Yields;p_{T} [GeV]; E d^{3}#sigma/dp^{3} [1/GeV^{2}]");
    mg->Add(photonSigmaGraph,"l3");
    mg->Add(expGraph);
    mg->SetMinimum(0.1);
    mg->Draw("AP");
    fa1->Draw("SAME");
    c1->BuildLegend();

    //plot details
    TPaveText *pt = new TPaveText(.6, .65, .85, .85, "NDC");
    pt->AddText("p-p 200 GeV");
    pt->AddText("-0.35 < #eta_{#gamma} < 0.35");
    pt->Draw("SAME");

    c1->SaveAs((analysisDir + "comparison.png").c_str());

    // Ratio Comparisons
    hPhotonPt->Add(hDecayPhotonPt); //making total
    hDecayPhotonPt->Divide(hPhotonPt); //getting relative amount
    auto photonRatioGraph = new TGraphErrors(hDecayPhotonPt);
    outf->cd();
    photonRatioGraph->Write("ratio");

    TCanvas* c2 = new TCanvas("c2", "Photon ratio comparison", 800, 600);
    c2->SetLogy();
    c2->SetLogx();
    photonRatioGraph->Draw("l3");
    c1->BuildLegend();
    c2->SaveAs((analysisDir + "ratio.png").c_str());

    outf->Close();
    return 0;
}
