
/*
 * Photon spectrum analysis for a single final-state hadrons file.
 *
 * - Reads /scratch/user/cameron.parker/projects/runs/pp-gamma/test_out_final_state_hadrons.dat
 * - Fills photon pT spectrum (pid==22)
 * - Loads Table 7 from ../external_packages/gammaLoss/data/pp2760PhotonSpectrum.root
 * - Plots both spectra for visual comparison
 * 
 * Usage: Just run, no arguments needed.
 * Output: photon_spectrum.root (contains spectrum), comparison.pdf (overlay plot)
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
        std::cerr << "ERROR: Table 7 not found in " << rootFilePath << std::endl;
        return nullptr;
    }
    TGraphAsymmErrors* g = (TGraphAsymmErrors*)tab7->Get("Graph1D_y1");
    if (!g) {
        std::cerr << "ERROR: Could not retrieve Table 7 Graph1D_y1 from ROOT file." << std::endl;
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
    const std::string rootDataFile = "../external_packages/gammaLoss/data/pp2760PhotonSpectrum.root";
    const std::string outputRootFile = analysisDir + "photon_spectrum.root";
    
    // Initialize Pythia8 instance (internal organization version).
    Pythia8::Pythia pythia;

    // Define pT bins to match Table 7 as closely as possible (adjust these to match binning of Table 7 if known)
    double etaCut = 0.9;
    std::vector<double> ptBins = {
        0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0,
        3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0, 16.0, 20.0
    };
    TH1D* hPhotonPt = new TH1D("photonPt", "Prompt photon p_{T} spectrum; p_{T} [GeV]; dN/dp_{T}", ptBins.size()-1, ptBins.data());

    // Data to compare to
    TFile* hadDataFile = TFile::Open("../external_packages/gammaLoss/data/pp2760hadrons.root");
    TDirectory* hadDataDir = (TDirectory*)hadDataFile->Get("Table 1");
    TGraphErrors* hadDataGraph =  (TGraphErrors*)hadDataDir->Get("Graph1D_y1");
    
    // Hadron histograms.
    TFile* hadFile = TFile::Open("../external_packages/gammaLoss/data/PbPb2760Hadrons.root");
    TDirectory* tableObj2 = (TDirectory*)hadFile->Get("Table 4");
    TH1D* hadDataSpecHist =  (TH1D*)tableObj2->Get("Hist1D_y1");
    TH1D* hadSpecHist = (TH1D*)hadDataSpecHist->Clone("HadronSpec");
    TH1D* etaHist = new TH1D("EtaDist", "Rapidity Distribution; #eta; Photons Out (a.u.)", 20, -2.5, 2.5);
    TH1D* phiHist = new TH1D("PhiDist", "Azimuthal Distribution; #phi; Photons Out (a.u.)", 20, 0, 2.*M_PI);
    hadSpecHist->Reset();
    double hadEtaCut = 1.0;

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

        //if(nEvents > 1000000) break; 
        
        // Loop through particles in event.
        for(auto& particle : particles)
        {
            // Check charge using the internal Pythia particleData.
            if(pythia.particleData.charge(particle.get()->pid()) != 0 && fabs(particle.get()->eta()) <= hadEtaCut)
            {
                hadSpecHist->Fill(particle.get()->pt());
                continue;
            }

            if(std::fabs(particle.get()->eta()) > etaCut)
                continue;

            if(particle.get()->pid() == 22)
            {
                nPhotons++;
                hPhotonPt->Fill(particle.get()->pt());
                etaHist->Fill(particle.get()->eta());
                phiHist->Fill(particle.get()->phi());
                continue;
            }
        } // end particle loop
    } // end event loop
    double sigmaGen = reader->GetSigmaGen();
    if(sigmaGen < 0) sigmaGen = 43.061;
    cout << "SigmaGen: " << sigmaGen << endl;
    cout << nEvents << " Events" << endl;

    // --- Save histogram to ROOT file
    TFile* outf = TFile::Open(outputRootFile.c_str(), "RECREATE");
    hPhotonPt->Scale(1.0/(2.0*M_PI*2.0*etaCut*nEvents), "width");
    ScaleByCenter(hPhotonPt);
    hPhotonPt->Write();

    // Spatial distributions
    etaHist->Write("eta dist");
    phiHist->Write("phi hist");

    TCanvas* cPhi = new TCanvas("cPhi", "Photon Azimuthal Distribution", 800, 600);
    phiHist->Scale(1./nEvents, "width");
    phiHist->SetMinimum(0.);
    phiHist->SetStats(0);
    phiHist->Draw("B");
    cPhi->SaveAs((analysisDir + "phi.pdf").c_str());
    
    TCanvas* cEta = new TCanvas("cPhi", "Photon Rapidity Distribution", 800, 600);
    etaHist->SetStats(0);
    etaHist->Scale(1./nEvents, "width");
    etaHist->Draw("B");
    cEta->SaveAs((analysisDir + "eta.pdf").c_str());

    // Process hadron spectra.
    hadSpecHist->Scale(1.0/(2.0*M_PI*2.0*etaCut*nEvents), "width");
    ScaleByCenter(hadSpecHist);
    TGraphErrors* hadronGraph = new TGraphErrors(hadSpecHist);
    hadSpecHist->Scale(sigmaGen);
    hadSpecHist->Write("pp hadron baseline");

    TH1D* photonSigma = (TH1D*)hPhotonPt->Clone("xsec");
    photonSigma->Scale(sigmaGen*1.0e9*1.0);
    photonSigma->Write();

    std::cout << "Found " << nPhotons << " photons in file.\n";
    std::cout << "Photon pT spectrum histogram written to " << outputRootFile << std::endl;
    
    // --- Load experimental spectrum for comparison
    TGraphAsymmErrors* expGraph = LoadTableGraph(rootDataFile, "Table 9");
    if (!expGraph) {
        std::cerr << "ERROR: Unable to retrieve experimental Table 7 graph for overlay!" << std::endl;
        return 2;
    }
    expGraph->MovePoints(-10,0);

    TGraphAsymmErrors* expSigmaGraph = LoadTableGraph(rootDataFile, "Table 5");
    if (!expSigmaGraph) {
        std::cerr << "ERROR: Unable to retrieve experimental Table 5 graph for overlay!" << std::endl;
        return 2;
    }
    expSigmaGraph->MovePoints(-10,0);
    expSigmaGraph->SetTitle("ALICE");
    expSigmaGraph->SetMarkerStyle(24);

    //making confidence interval graphic
    std::vector<double> sigmabins = {
        0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.8, 
        3.25, 3.75, 4.5, 5.5, 9.0
    };

    std::vector<double> sigmaLimits = {
        7.68e9, 9.22e8, 3.82e8, 1.46e8, 4.37e7, 1.6e7, 1.1e7, 6.72e6, 5.25e6,
        3.04e6, 1.26e6, 8.21e5, 5.98e5, 1.79e5, 9.85e4, 3.69e4, 2.92e3
    };

    std::vector<double> sigmaLimitsLowers = {
        7.68e9, 9.22e8, 3.82e8, 1.46e8, 4.37e7, 1.6e7, 1.1e7, 6.72e6, 5.25e6,
        3.04e6, 1.26e6, 8.21e5, 5.98e5, 1.79e5, 9.85e4, 3.69e4, 2.92e3
    };

    std::vector<double> sigmaLimitsUppers = {
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0
    };
    
    std::vector<double> sigmaLimitsWidths = {
        0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
        0.1, 0.1, 0.2, 0.25, 0.25, 0.5, 0.5, 1.
    };

    TGraphAsymmErrors* limits = new TGraphAsymmErrors(sigmabins.size(),
        sigmabins.data(), sigmaLimits.data(), sigmaLimitsWidths.data(), sigmaLimitsWidths.data(),
        sigmaLimitsLowers.data(), sigmaLimitsUppers.data());


    // --- Create comparison plot
    TCanvas* c1 = new TCanvas("c1", "Photon spectrum comparison", 800, 600);
    c1->SetLogy();
    c1->SetLogx();

    hPhotonPt->SetLineColor(kBlue+1);
    hPhotonPt->SetMarkerStyle(20);
    hPhotonPt->SetMarkerColor(kBlue+1);
    hPhotonPt->SetFillColorAlpha(kBlue+1, 0.35);
    hPhotonPt->SetLineWidth(2);
    auto hPhotonPtGraph = new TGraphErrors(hPhotonPt);

    expGraph->SetMarkerStyle(24);
    expGraph->SetLineWidth(2);

    auto mg1 = new TMultiGraph;
    mg1->SetTitle("Photon Yields;p_{T}; dN/dp_{T}");
    mg1->Add(expGraph);
    mg1->Add(hPhotonPtGraph,"l3");
    mg1->Draw("AP");
    c1->BuildLegend();

    c1->SaveAs((analysisDir + "yields.pdf").c_str());

    // --- Create sigma plot
    TCanvas* c2 = new TCanvas("c2", "Photon sigma comparison", 800, 600);
    c2->SetLogy();
    c2->SetLogx();

    auto fa1 = new TF1("parameterization","2*6495*2760/((0.8+x)^5)",.1,15); //https://arxiv.org/pdf/hep-ph/0308085 eq 17
    fa1->SetTitle("Parameterization");

    photonSigma->SetLineColor(kBlue+1);
    photonSigma->SetMarkerStyle(20);
    photonSigma->SetMarkerColor(kBlue+1);
    photonSigma->SetFillColorAlpha(kBlue+1, 0.35);
    photonSigma->SetLineWidth(2);
    auto photonSigmaGraph = new TGraphErrors(photonSigma);
    photonSigmaGraph->SetTitle("X-SCAPE");

    photonSigma->Scale(5.0);
    photonSigma->SetLineColor(kBlue+3);
    photonSigma->SetMarkerColor(kBlue+3);
    photonSigma->SetFillColorAlpha(kBlue+3, 0.35);
    auto kfactorgraph = new TGraphErrors(photonSigma);
    kfactorgraph->SetTitle("K factor = 5");

    auto mg2 = new TMultiGraph;
    mg2->SetTitle("p-p Direct Photon Yields;p_{T} [GeV]; E d^{3}#sigma/dp^{3} [pb/GeV^{2}]");
    mg2->Add(expSigmaGraph);
    mg2->Add(photonSigmaGraph,"l3");
    mg2->Add(kfactorgraph,"l3");
    mg2->Add(limits);
    mg2->SetMinimum(1e1);
    mg2->Draw("AP");
    fa1->Draw("SAME");
    
    auto legend = new TLegend();
    legend->AddEntry(expSigmaGraph);
    legend->AddEntry(photonSigmaGraph, "", "lf");
    legend->AddEntry(kfactorgraph, "", "lf");
    legend->AddEntry(fa1);
    legend->Draw();

    //plot details
    TPaveText *pt = new TPaveText(.6, .65, .85, .85, "NDC");
    pt->AddText("p-p 2.76 TeV");
    pt->AddText("-0.9 < #eta_{#gamma} < 0.9");
    pt->SetFillColor(kWhite);
    pt->Draw("SAME NDC");

    c2->SaveAs((analysisDir + "sigma.pdf").c_str());

    // --- Create hadron plot
    TCanvas* c3 = new TCanvas("c3", "Hadron sigma comparison", 800, 600);
    c3->SetLogy();
    c3->SetLogx();

    hadronGraph->SetLineColor(kBlue+1);
    hadronGraph->SetMarkerStyle(20);
    hadronGraph->SetMarkerColor(kBlue+1);
    hadronGraph->SetFillColorAlpha(kBlue+1, 0.35);
    hadronGraph->SetLineWidth(2);

    auto mg3 = new TMultiGraph;
    mg3->SetTitle("Hadron Yields;p_{T}; E d^{3}N/dp^{3}");
    mg3->Add(hadDataGraph);
    mg3->Add(hadronGraph,"l3");
    mg3->Draw("AP");
    mg3->SetMinimum(1.e-14);
    c3->BuildLegend();

    c3->SaveAs((analysisDir + "hadrons.pdf").c_str());

    outf->Close();
    return 0;
}
