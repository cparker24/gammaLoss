/*******************************************************************************
 * Copyright (c) The JETSCAPE Collaboration, 2018
 *
 * Modular, task-based framework for simulating all aspects of heavy-ion collisions
 *
 * For the list of contributors see AUTHORS.
 *
 * Report issues at https://github.com/JETSCAPE/JETSCAPE/issues
 *
 * or via email to bugs.jetscape@gmail.com
 *
 * Distributed under the GNU General Public License 3.0 (GPLv3 or later).
 * See COPYING for details.
 ******************************************************************************/

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
#include "TGraphAsymmErrors.h"
#include "TPaveText.h"

// Organization-specific includes (Pythia8 and reader)
#include "JetScapeReaderFinalStateHadrons.h"
#include "JetScapeBanner.h" // if needed
#include "JetScapeParticles.h"
#include "JetScape.h"
#include "JetScapeXML.h"
#include "Pythia8/Pythia.h"
#include "analysis.h"
#include "QVecV2.h"

using namespace Jetscape;

int main(int argc, char** argv)
{
    if(argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <analysisDir>" << std::endl;
        return 1;
    }

    bool abson = false;
    if(argc == 3){
        if(std::string(argv[2]) == "--abson")
            abson = true;
    }

    
    std::string analysisDir = argv[1];
    std::string pointsDir = analysisDir + "points/";
    std::string outDir = analysisDir + "analysis/";

    if(abson) outDir = analysisDir + "abs-analysis/";
    
    //run settings
    auto jetscape = make_shared<JetScape>();
    const char* mainXMLName = "../config/jetscape_main.xml";
    std::string  userXMLName = analysisDir + "/config.xml";
    jetscape->SetXMLMainFileName(mainXMLName);
    jetscape->SetXMLUserFileName(userXMLName.c_str());

    JetScapeXML::Instance()->OpenXMLMainFile(mainXMLName);
    JetScapeXML::Instance()->OpenXMLUserFile(userXMLName.c_str());
    double emissionRate = JetScapeXML::Instance()->GetElementDouble({"Photons", "gammaLoss", "emissionFactor"});
    bool absorption = JetScapeXML::Instance()->GetElementDouble({"Photons", "gammaLoss", "reabsorption"});
    // resampled photons per event: thermal (gammaLoss) and prompt (PhotonGun)
    double thermalResample = JetScapeXML::Instance()->GetElementDouble({"Photons", "gammaLoss", "resample"});
    double hardResample = JetScapeXML::Instance()->GetElementDouble({"Hard", "PhotonGun", "resample"});

    cout << "Emission Rate: " << emissionRate << endl;
    cout << "Reabsorption: " << absorption << endl;

    // moving to output dir
    chdir(outDir.c_str());
    std::filesystem::create_directory("fits");
    
    // Initialize Pythia8 instance (internal organization version).
    Pythia8::Pythia pythia;
    
    // Open data files for overlay histograms.
    TFile* specFile = TFile::Open("../external_packages/gammaLoss/data/PbPb2760Spectra.root");
    TFile* photonDatav2file = TFile::Open("../external_packages/gammaLoss/data/PbPb2760v2.root");
    TFile* hadDatav2file   = TFile::Open("../external_packages/gammaLoss/data/PbPb2760HadronV2.root");
    TFile* hadronDataSpecFile = TFile::Open("../external_packages/gammaLoss/data/PbPb2760Hadrons.root");
    // pp baseline from the analysis directory if present, otherwise the default run
    std::string ppBaselineFile = analysisDir + "/pp_photon_spectrum.root";
    if(!std::ifstream(ppBaselineFile)) ppBaselineFile = "/scratch/user/cameron.parker/projects/runs/pp/pp2760baseline/photon_spectrum.root";
    TFile* ppFile = TFile::Open(ppBaselineFile.c_str());
    
    // Retrieve graphs from files.
    bool central = true;
    if (analysisDir.find("20-40") != std::string::npos) 
        central = false;

    TGraphAsymmErrors* photonDatav2;
    TGraphAsymmErrors* specData;
    TGraphAsymmErrors* hadDatav2;
    TGraphAsymmErrors* hadDataSpec;
    TH1D* pThist;

    if(central){
        photonDatav2 = (TGraphAsymmErrors*)getFileObject(photonDatav2file, "Table 7", "Graph1D_y1");
        specData = (TGraphAsymmErrors*)getFileObject(specFile, "Table 1", "Graph1D_y1");
        pThist = replicateHistogram(specFile, "Table 1", "Hist1D_y1");
        hadDatav2 =  (TGraphAsymmErrors*)getFileObject(hadDatav2file, "Table 2", "Graph1D_y1");
        hadDataSpec =  (TGraphAsymmErrors*)getFileObject(hadronDataSpecFile, "Table 11", "Graph1D_y1");
    }
    else{
        photonDatav2 = (TGraphAsymmErrors*)getFileObject(photonDatav2file, "Table 8", "Graph1D_y1");
        specData = (TGraphAsymmErrors*)getFileObject(specFile, "Table 2", "Graph1D_y1");
        pThist = replicateHistogram(specFile, "Table 2", "Hist1D_y1");
        hadDatav2 =  (TGraphAsymmErrors*)getFileObject(hadDatav2file, "Table 2", "Graph1D_y3");
        hadDataSpec =  (TGraphAsymmErrors*)getFileObject(hadronDataSpecFile, "Table 12", "Graph1D_y1");
    }
    specData->SetTitle("ALICE Data; p_{T} [GeV]; 1/(2 #pi N_{ev}) d^{2} N/p_{T}dp_{T} d#eta [GeV^{-2}]");

    TH1D* photonDatav2Hist = (TH1D*)getFileObject(photonDatav2file, "Table 8", "Hist1D_y1");
    TH1D* hadDatav2Hist = (TH1D*)getFileObject(hadDatav2file, "Table 2", "Hist1D_y2");
    TH1D* hadSpecHist = replicateHistogram(hadronDataSpecFile, "Table 4", "Hist1D_y1");

    cout << "Data graphs acquired" << endl;
    
    // Define pT and phi bin edges.
    std::vector<double> hadpTs = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,3.0,3.3,3.7,4.1,4.6,5.4,6.2};
    std::vector<double> phis;
    int nPhiBins = 20;
    const double twoPi = 2 * M_PI;
    for(int i = 0; i <= nPhiBins; i++){
        phis.push_back(i * M_PI/10);
    }
    
    // Extended range
    double extendedBins[] = {0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,3.0,4.0,5.0,7.0,10.0};
    int nExtendedBins = 15;

    double hadBins[] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.5,3.0,3.5,4.0,5.0};
    int nHadBins = 18;

    double v2bins[] = {0.2,0.3,0.5,0.7,0.9,1.1,1.3,1.5,2.0,2.5,3.0,3.5,4.0,5.0};
    int nv2bins = 13;
    
    // Create histograms.
    TH2D* totalhist = new TH2D("PhotonRate", "Photon Counts; #phi; p_{T} [GeV]; Photons Out",
                               phis.size()-1, &phis[0],
                               photonDatav2Hist->GetNbinsX(), photonDatav2Hist->GetXaxis()->GetXbins()->GetArray());
    TH2D* exthist = new TH2D("ExtPhotonRate", "Photon Counts; #phi; pT [GeV]; Photons Out",
                               phis.size()-1, &phis[0],
                               nv2bins, v2bins);
    TH2D* QGPspread = new TH2D("QGPPhotonRate", "Total; phi; pT (GeV); Photons Out",
                               phis.size()-1, &phis[0],
                               nv2bins, v2bins);
    TH2D* HRGspread = new TH2D("HRGPhotonRate", "Total; phi; pT (GeV); Photons Out",
                               phis.size()-1, &phis[0],
                               nv2bins, v2bins);
    TH2D* Hardspread = new TH2D("HardPhotonRate", "Total; phi; pT (GeV); Photons Out",
                               phis.size()-1, &phis[0],
                               nv2bins, v2bins);
    TH1D* QGPHist = new TH1D("QGPRate", "QGP; pT (GeV); Photons Out",
                             15, extendedBins);
    TH1D* HRGHist = new TH1D("HRGRate", "Hadronic; pT (GeV); Photons Out",
                             15, extendedBins);
    TH1D* showerHist = new TH1D("ShowerRate", "Shower; pT (GeV); Photons Out",
                                15, extendedBins);
    TH1D* hardHist = new TH1D("HardProcessRate", "Prompt; pT (GeV); Photons Out",
                              15, extendedBins);
    TH2D* otherhist = new TH2D("OtherDistribution", "Hadron Counts; #phi; p_{T} (GeV); Hadrons Out",
                               phis.size()-1, &phis[0],
                               nHadBins, hadBins);
    TH1D* etaHist = new TH1D("EtaDist", "Photon Counts; eta; Photons Out", 20, -2.5, 2.5);
    TH1D* promptetaHist = new TH1D("PromptEtaDist", "Photon Counts; eta; Photons Out", 20, -2.5, 2.5);
    TH1D* QGPetaHist = new TH1D("QGPEtaDist", "Photon Counts; eta; Photons Out", 20, -2.5, 2.5);
    TH1D* HRGetaHist = new TH1D("HRGEtaDist", "Photon Counts; eta; Photons Out", 20, -2.5, 2.5);

    QVecV2* photonQv2 = new QVecV2("photonhist", "Direct Photon v_{2}; p_{T} [GeV]; v_{2}",
                             nv2bins, v2bins);

    QVecV2* hadronQv2 = new QVecV2("hadronhist", "Charged Hadron v_{2}; p_{T} [GeV]; v_{2}",
                             nHadBins, hadBins);

    cout << "Analysis graphs made" << endl;
    
    // Output ROOT file.
    TFile* outfile = TFile::Open((outDir + "analysis.root").c_str(), "RECREATE");
    
    // Get list of directories (each representing one analysis point).
    std::vector<std::string> directories = getDirs(pointsDir);
    
    int eventCount = 0;
    double etaCut = 0.8;
    double hardFactor = 0.0735818/61.7678/hardResample; //prompt photon over inelastic xsec, per PhotonGun photon
    double issFactor = 1.0/JetScapeXML::Instance()->GetElementDouble({"SoftParticlization", "iSS", "number_of_repeated_sampling"});
    double thermalFactor = 1.0/(emissionRate*thermalResample);
    double jetFactor = 1./100.;

    // Ncoll calculations
    double averageNcoll = 0.0;
    
    int nphotons = 0;
    int nabsphotons = 0;

    // Loop over each directory (analysis point).
    for(const auto& dirName : directories)
    {
        std::string filePath = pointsDir + dirName + "/test_out_final_state_hadrons.dat";
        // Create reader for this file.
        auto reader = make_shared<Jetscape::JetScapeReaderFinalStateHadrons>(filePath);
        double EPangle = -10.0;

        int Ncoll = GetNcoll(pointsDir + dirName);
        cout << "Event Ncoll: " << Ncoll << endl;
        averageNcoll += Ncoll;

        double fileEvents = reader->TotalEventCount();
        cout << "Events in file: " << fileEvents << endl;

        //if(fileEvents < 2) continue;
        jetFactor = 1.0/fileEvents;
        
        while(!reader->Finished())
        {
            try{
                cout << dirName << ": ";
                reader->Next();
            }
            catch(...){
                break;
            }
            eventCount++;
            std::vector<std::shared_ptr<Jetscape::Hadron>> particles = reader->GetHadrons();
            if(particles.empty())
                break;
            
            if(EPangle < -9)
            {
                EPangle = GetEpAngle(particles);
                std::cout << "Event Plane Angle: " << EPangle << std::endl;
            }
            
            // Loop through particles in event.
            for(auto& particle : particles)
            {
                double weight = 1.0;
                double delphi = particle->phi() - EPangle;
                if(delphi > twoPi)
                    delphi -= twoPi;
                if(delphi < 0)
                    delphi += twoPi;
                
                // Check for photon selection based on pid and status.
                if((particle.get()->pid() == 22 && (abs(particle.get()->pstat()) == 24 || abs(particle.get()->pstat()) == 25)))
                {
                    weight = thermalFactor; //photon production increase
                }
                
                if((particle.get()->pstat() == 11 || particle.get()->pstat() == 27) 
                    && std::fabs(particle.get()->eta()) < etaCut && pythia.particleData.charge(particle.get()->pid()) != 0) //smash tag
                {
                    weight = issFactor; //iSS resampler count
                    otherhist->Fill(delphi, particle.get()->pt(), weight);
                    hadronQv2->Fill(particle.get()->pt(), weight, delphi);
                }
                
                //starting photon section now
                if(std::fabs(particle.get()->eta()) > etaCut)
                    continue;

                // Check charge using the internal Pythia particleData.
                if(particle.get()->pstat() == 11|| particle.get()->pstat() == 27)
                {
                    if(pythia.particleData.charge(particle.get()->pid()) != 0)
                    {
                        hadSpecHist->Fill(particle.get()->pt(), weight);
                    }
                    continue;
                }

                // Adding shower hadrons to v2
                if(particle.get()->pstat() == 0)
                {
                    if(pythia.particleData.charge(particle.get()->pid()) != 0)
                    {
                        weight = jetFactor;
                        hadSpecHist->Fill(particle.get()->pt(), weight);
                        otherhist->Fill(delphi, particle.get()->pt(), weight);
                    }
                    continue;
                }
                if(particle.get()->pstat() == -1)
                {
                    if(pythia.particleData.charge(particle.get()->pid()) != 0)
                    {
                        weight = -1.*jetFactor;
                        hadSpecHist->Fill(particle.get()->pt(), weight);
                        otherhist->Fill(delphi, particle.get()->pt(), weight);
                    }
                    continue;
                }

                //photon section now 
                if(particle.get()->pid() == 22)
                {
                    if(particle.get()->pstat() < 0) nabsphotons++;
                    nphotons++;
                    //skipping absorbed ones in abs runs
                    if(abson == true and particle.get()->pstat() < 0){
                        continue;
                    }

                    if(abs(particle.get()->pstat()) == 22)
                    {
                        weight = hardFactor;
                        hardHist->Fill(particle.get()->pt(), weight);
                        exthist->Fill(delphi, particle.get()->pt(), weight);
                        Hardspread->Fill(delphi, particle.get()->pt(), weight);
                        if(particle.get()->pt() > 1.0) promptetaHist->Fill(particle.get()->eta(), weight);

                        if(absorption)
                            totalhist->Fill(delphi, particle.get()->pt(), weight);
                        else
                            smoothFill(totalhist, particle.get()->pt(), weight);
                    }

                    if(abs(particle.get()->pstat()) == 23)
                    {
                        weight = hardFactor;
                        showerHist->Fill(particle.get()->pt(), weight);
                        hardHist->Fill(particle.get()->pt(), weight);
                        Hardspread->Fill(delphi, particle.get()->pt(), weight);
                        totalhist->Fill(delphi, particle.get()->pt(), weight);
                        exthist->Fill(delphi, particle.get()->pt(), weight);
                        if(particle.get()->pt() > 1.0) promptetaHist->Fill(particle.get()->eta(), weight);
                    }

                    if(abs(particle.get()->pstat()) == 24){
                        QGPHist->Fill(particle.get()->pt(), weight);
                        QGPspread->Fill(delphi, particle.get()->pt(), weight);
                        totalhist->Fill(delphi, particle.get()->pt(), weight);
                        exthist->Fill(delphi, particle.get()->pt(), weight);
                        if(particle.get()->pt() > 1.0) QGPetaHist->Fill(particle.get()->eta(), weight);
                    }
                    
                    if(abs(particle.get()->pstat()) == 25){
                        HRGHist->Fill(particle.get()->pt(), weight);
                        HRGspread->Fill(delphi, particle.get()->pt(), weight);
                        totalhist->Fill(delphi, particle.get()->pt(), weight);
                        exthist->Fill(delphi, particle.get()->pt(), weight);
                        if(particle.get()->pt() > 1.0) HRGetaHist->Fill(particle.get()->eta(), weight);
                    }
                    pThist->Fill(particle.get()->pt(), weight);
                    photonQv2->Fill(particle.get()->pt(), weight, delphi);
                    
                    if(particle.get()->pt() > 1.0) etaHist->Fill(particle.get()->eta(), weight);
                }
            } // end particle loop
        } // end event loop
        reader->Close();
    } // end directory loop

    // Scale and write histograms.
    if(eventCount > 0)
    {
        // Final Ncoll
        averageNcoll /= directories.size();
        cout << "Average Ncoll: " << averageNcoll << endl;

        cout << "Total photons: " << nphotons << endl;
        cout << "Absorbed photons: " << nabsphotons << endl;
        
        totalhist->Scale(1.0/eventCount);
        totalhist->Write();
        
        otherhist->Scale(1.0/eventCount);
        otherhist->Write();

        exthist->Scale(1.0/eventCount, "width");
        exthist->Write();
        
        QGPHist->SetLineColor(kPink);
        QGPHist->SetFillColorAlpha(kPink, 0.35);
        QGPHist->SetLineWidth(2);
        QGPHist->Scale(1.0/(2.0*M_PI*2.0*etaCut*eventCount), "width");
        ScaleByCenter(QGPHist);
        QGPHist->Write();
        ExtractTemp(QGPHist);
        TGraphAsymmErrors* QGPGraph = new TGraphAsymmErrors(QGPHist);
        
        HRGHist->SetLineColor(kViolet);
        HRGHist->SetFillColorAlpha(kViolet, 0.35);
        HRGHist->SetLineWidth(2);
        HRGHist->Scale(1.0/(2.0*M_PI*2.0*etaCut*eventCount), "width");
        ScaleByCenter(HRGHist);
        HRGHist->Write();
        TGraphAsymmErrors* HRGGraph = new TGraphAsymmErrors(HRGHist);
        
        hardHist->SetLineColor(kGreen);
        hardHist->SetFillColorAlpha(kGreen, 0.35);
        hardHist->SetLineWidth(2);
        hardHist->Scale(1.0/(2.0*M_PI*2.0*etaCut*eventCount), "width");
        ScaleByCenter(hardHist);
        hardHist->Write();
        TGraphAsymmErrors* hardGraph = new TGraphAsymmErrors(hardHist);
        
        showerHist->SetLineColor(kOrange);
        showerHist->SetFillColorAlpha(kOrange, 0.35);
        showerHist->SetLineWidth(2);
        showerHist->Scale(1.0/(2.0*M_PI*2.0*etaCut*eventCount), "width");
        ScaleByCenter(showerHist);
        showerHist->Write();
        TGraphAsymmErrors* showerGraph = new TGraphAsymmErrors(showerHist);
        
        etaHist->Write();
        promptetaHist->Write();
        QGPetaHist->Write();
        HRGetaHist->Write();

        //plot details
        TPaveText *pt = new TPaveText(.15, .1, .4, .3, "NDC");
        pt->AddText("Pb-Pb 2.76 TeV");
        if(central) pt->AddText("Centrality 0-20%");
        else pt->AddText("Centrality 20-40%");
        pt->AddText("-0.8 < #eta_{#gamma} < 0.8");
        pt->SetFillColor(kWhite);

        
        TPaveText *pthad = new TPaveText(.15, .1, .4, .3, "NDC");
        pthad->AddText("Pb-Pb 2.76 TeV");
        pthad->AddText("Centrality 20-40%");
        pthad->AddText("-0.8 < #eta < 0.8");
        pthad->SetFillColor(kWhite);

        // Create pT spectrum projection.
        pThist->Scale(1.0/(2.0*M_PI*2.0*etaCut*eventCount), "width");
        pThist->SetFillColorAlpha(kBlue, 0.35);
        pThist->SetLineWidth(2);
        ScaleByCenter(pThist);
        pThist->Write();
        ExtractTemp(pThist);
        TGraphAsymmErrors* totalGraph = new TGraphAsymmErrors(pThist);
        
        // Process hadron spectra.
        hadSpecHist->Scale(1.0/(2.0*M_PI*2.0*etaCut*directories.size()), "width");
        ScaleByCenter(hadSpecHist);
        hadSpecHist->SetFillColorAlpha(kBlue, 0.35);
        hadSpecHist->SetLineColor(kBlue);
        hadSpecHist->SetLineWidth(3);
        TGraphAsymmErrors* hadSpecGraph = new TGraphAsymmErrors(hadSpecHist);
        hadSpecHist->Write("hadrons");
        
        // Draw hadron spectrum
        hadDataSpec->SetTitle("ALICE Data; p_{T} [GeV]; 1/(2 #pi N_{ev}) d^{2} N/p_{T}dp_{T} d#eta [GeV^{-2}]");
        myRatioPlot(hadDataSpec, hadSpecGraph, "Charged Hadron Spectrum", {}, true, true, pthad);

        // Make Raa plot
        TH1D* ppHads = (TH1D*)ppFile->Get("pp hadron baseline");
        TH1D* Raa = (TH1D*)hadSpecHist->Clone("Raa");
        Raa->Divide(ppHads);
        Raa->Scale(1./10.);
        Raa->Write("Raa charged hadrons");
        TGraphAsymmErrors* RaaGraph = new TGraphAsymmErrors(Raa);
        TGraphAsymmErrors* raaData = (TGraphAsymmErrors*)getFileObject(hadronDataSpecFile,"Table 20", "Graph1D_y1");

        TCanvas* raaCanvas = new TCanvas("Raa Canvas", "Hadron Ratio", 800, 600);
        raaCanvas->SetLogy();
        raaCanvas->SetLogx();
        raaData->Draw();
        RaaGraph->Draw("SAME l3");
        raaCanvas->Print((outDir + "Raa.png").c_str());
        
        //test new v2 method
        auto QVecv2 = photonQv2->GetV2Graph();
        QVecv2->Write("photon SPv2");
        QVecv2->SetTitle("SP v_{2}");

        auto hadronQVecV2 = hadronQv2->GetV2Graph();
        hadronQVecV2->Write("hadron SPv2");
        hadronQVecV2->SetTitle("SP v_{2}");

        // Drawing v2
        TGraphAsymmErrors* flowHadsFit = MakevnGraph(otherhist, "Fit Hadron v_{2}", 2, true);
        flowHadsFit->Write();
        flowHadsFit->SetLineColor(kRed);
        flowHadsFit->SetLineWidth(2);
        flowHadsFit->SetFillColorAlpha(kRed, 0.35);

        TGraphAsymmErrors* flowHads = MakevnGraph(otherhist, "Integral Hadron v_{2}");
        flowHads->SetLineColor(kGreen+2);
        flowHads->SetLineWidth(2);
        flowHads->SetFillColorAlpha(kGreen+2, 0.35);
        flowHads->Write();
        hadDatav2->SetTitle("ALICE Data; p_{T} [GeV]; v_{2}^{hadron}");
        myRatioPlot(hadDatav2, hadronQVecV2, "Charged Hadron v_{2}", {flowHadsFit, flowHads}, false, false, pthad);
        
        TGraphAsymmErrors* flowPhotons = MakevnGraph(totalhist, "Integral Photon v_{2}");
        flowPhotons->SetLineColor(kGreen+2);
        flowPhotons->SetLineWidth(2);
        flowPhotons->SetFillColorAlpha(kGreen+2, 0.35);
        flowPhotons->Write();
        
        TGraphAsymmErrors* fitPhotons = MakevnGraph(exthist, "Fit Photon v_{2}", 2, true);
        fitPhotons->Write();
        fitPhotons->SetLineColor(kRed);
        fitPhotons->SetLineWidth(2);
        fitPhotons->SetFillColorAlpha(kRed, 0.35);

        photonDatav2->SetTitle("ALICE Data; p_{T} [GeV]; v_{2}^{#gamma}");
        QVecv2->SetTitle("SP photon v_{2}");
        myRatioPlot(photonDatav2, QVecv2, "Photon v_{2}", {fitPhotons, flowPhotons}, false, false, pt);

        TGraphAsymmErrors* extPhotons = MakevnGraph(exthist, "Extended Photon v2");
        extPhotons->Write();
        
        TCanvas* specCanvas = new TCanvas("specCanvas", "Photon Spec", 800, 600);
        specCanvas->SetLogy();
        vector<TGraphAsymmErrors*> photonGraphs = {};
        // pp baseline per inelastic collision: dsigma/dpT [pb] / sigma_inel
        TH1D* ppSigma = (TH1D*)ppFile->Get("xsec");
        ppSigma->Scale(1.0/(61.7678*1.0e9));
        TGraphAsymmErrors* ppBase = new TGraphAsymmErrors(ppSigma);
        ppBase->SetTitle("p-p baseline scaled");
        ppBase->SetLineColor(kGreen-2);
        ppBase->SetFillColorAlpha(kGreen-2,0.35);
        ppBase->SetLineWidth(2);
        ppBase->Scale(averageNcoll*1.); //value from https://journals.aps.org/prc/abstract/10.1103/PhysRevC.93.044906 Table I
        photonGraphs.push_back(hardGraph);
        photonGraphs.push_back(QGPGraph);
        photonGraphs.push_back(HRGGraph);
        photonGraphs.push_back(ppBase);
        totalGraph->SetTitle("X-SCAPE");
        myRatioPlot(specData, totalGraph, "Direct Photon p_{T} Spectrum", photonGraphs, true, true, pt);
        totalGraph->Write("total");

        //other v2s
        TGraphAsymmErrors* QGPv2 = MakevnGraph(QGPspread, "QGP v2");
        QGPv2->Write("QGPv2");
        
        TGraphAsymmErrors* HRGv2 = MakevnGraph(HRGspread, "HRG v2");
        HRGv2->Write("HRGv2");
        
        TGraphAsymmErrors* Hardv2 = MakevnGraph(Hardspread, "Hard v2");
        Hardv2->Write("Hardv2");
    }
    else
    {
        std::cerr << "No events processed!" << std::endl;
    }
    
    outfile->Close();
    std::cout << "Analysis complete. Output written to " << outDir + "analysis.root" << std::endl;
    
    return 0;
}