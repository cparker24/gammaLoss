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
#include "TMultiGraph.h"
#include "TGraphAsymmErrors.h"
#include "TLine.h"
#include "TLegend.h"
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
    std::string pointsDir = analysisDir + "/points/";
    std::string outDir = analysisDir + "/analysis/";

    if(abson) outDir = analysisDir + "abs-analysis/";

    //run settings
    auto jetscape = make_shared<JetScape>();
    const char* mainXMLName = "../config/jetscape_main.xml";
    std::string  userXMLName = analysisDir + "/config.xml";
    jetscape->SetXMLMainFileName(mainXMLName);
    jetscape->SetXMLUserFileName(userXMLName.c_str());

    JetScapeXML::Instance()->OpenXMLMainFile(mainXMLName);
    JetScapeXML::Instance()->OpenXMLUserFile(userXMLName.c_str());
    double emissionRate = JetScapeXML::Instance()->GetElementDouble({"Eloss", "gammaLoss", "emissionFactor"});

    cout << "Emission Rate: " << emissionRate << endl;

    bool central = true;
    if (analysisDir.find("20-40") != std::string::npos) 
        central = false;
    
    // Initialize Pythia8 instance (internal organization version).
    Pythia8::Pythia pythia;
    
    // Open data files for overlay histograms.
    TFile* photonDatav2file = TFile::Open("../external_packages/gammaLoss/data/AuAu200photonv2.root");
    TFile* hadFile = TFile::Open("../external_packages/gammaLoss/data/AuAu200IDHads.root");
    TFile* hadv2File = TFile::Open("../external_packages/gammaLoss/data/AuAu200hadronv2.root");
    TFile* specFileSTAR = TFile::Open("../external_packages/gammaLoss/data/AuAu200photonsSTAR.root");
    TFile* specFile = TFile::Open("../external_packages/gammaLoss/data/AuAu200PhotonSpectrum.root");
    TFile* ppPhotonBaseline = TFile::Open("/scratch/user/cameron.parker/projects/runs/pp/pp200/photon_spectrum.root");

    // Retrieve graphs from files.
    TGraphAsymmErrors* photonDatav2;
    TGraphAsymmErrors* specData;
    TGraphAsymmErrors* specDataSTAR;
    TGraphAsymmErrors* hadDatav2;
    TGraphAsymmErrors* pionGraph;

    if(central){
        photonDatav2 = (TGraphAsymmErrors*)getFileObject(photonDatav2file, "Figure 2b", "Graph1D_y9");
        specData = (TGraphAsymmErrors*)getFileObject(specFile, "Figure 15", "Graph1D_y1");
        specDataSTAR = (TGraphAsymmErrors*)getFileObject(specFileSTAR, "Table 18", "Graph1D_y1");
        hadDatav2 =  (TGraphAsymmErrors*)getFileObject(hadv2File, "Figure 3a", "Graph1D_y2");
        pionGraph =  (TGraphAsymmErrors*)getFileObject(hadFile, "Figure 1", "Graph1D_y17");
    }
    else{
        photonDatav2 = (TGraphAsymmErrors*)getFileObject(photonDatav2file, "Figure 2b", "Graph1D_y11");
        specData = (TGraphAsymmErrors*)getFileObject(specFile, "Figure 15", "Graph1D_y2");
        specDataSTAR = (TGraphAsymmErrors*)getFileObject(specFileSTAR, "Table 17", "Graph1D_y1");
        hadDatav2 =  (TGraphAsymmErrors*)getFileObject(hadv2File, "Figure 3g", "Graph1D_y2");
        pionGraph =  (TGraphAsymmErrors*)getFileObject(hadFile, "Figure 1", "Graph1D_y18");
    }

    specDataSTAR->SetTitle("STAR data");
    specData->SetTitle("PHENIX data");
    photonDatav2->SetPointError(3,0.25,0.25,0.2,0.2);
    pionGraph->RemovePoint(15); pionGraph->RemovePoint(15); //deleting random duplicates

    TGraphAsymmErrors* ppbase = new TGraphAsymmErrors((TH1D*)ppPhotonBaseline->Get("dNdpt"));
    ppbase->SetTitle("N_{coll} scaled pp baseline");
    ppbase->SetLineColor(kGreen-2);
    ppbase->SetFillColorAlpha(kGreen-2,0.35);
    ppbase->SetLineWidth(2);

    cout << "Data graphs acquired" << endl;

    // Hadron histograms.
    TH1D* hadSpecHist = replicateHistogram(hadFile, "Figure 1", "Hist1D_y18");
    TH1D* hadDatav2Hist = (TH1D*)getFileObject(hadv2File, "Figure 3g", "Hist1D_y2");
    TH1D* pThist = replicateHistogram(specFile, "Figure 15", "Hist1D_y2");
    TH1D* photonDatav2Hist = (TH1D*)getFileObject(photonDatav2file, "Figure 2b", "Hist1D_y11");

    // Define pT and phi bin edges.
    std::vector<double> phis;
    int nPhiBins = 20;
    const double twoPi = 2 * M_PI;
    for(int i = 0; i <= nPhiBins; i++){
        phis.push_back(i * M_PI/10);
    }

    // Extended range
    double extendedBins[] = {0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,3.0,4.0,5.0,6.0,8.0,10.0,12.0};
    int nExtendedBins = 17;
    
    double v2bins[] = {0.2,0.3,0.5,0.7,0.9,1.1,1.3,1.5,2.0,2.5,3.0,3.5,4.0,5.0};
    int nv2bins = 13;
    
    // Create histograms.
    TH2D* totalhist = new TH2D("PhotonRate", "Total; phi; pT (GeV); Photons Out",
                               phis.size()-1, &phis[0],
                               photonDatav2Hist->GetNbinsX(), photonDatav2Hist->GetXaxis()->GetXbins()->GetArray());
    TH2D* exthist = new TH2D("ExtPhotonRate", "Total; phi; pT (GeV); Photons Out",
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
    TH1D* extendedHist = new TH1D("total", "total; pT (GeV); Photons Out",
                             nExtendedBins, extendedBins);
    TH1D* QGPHist = new TH1D("QGPRate", "QGP; pT (GeV); Photons Out",
                             nExtendedBins, extendedBins);
    TH1D* HRGHist = new TH1D("HRGRate", "Hadronic; pT (GeV); Photons Out",
                             nExtendedBins, extendedBins);
    TH1D* showerHist = new TH1D("ShowerRate", "Shower; pT (GeV); Photons Out",
                                nExtendedBins, extendedBins);
    TH1D* hardHist = new TH1D("HardProcessRate", "Prompt; pT (GeV); Photons Out",
                              nExtendedBins, extendedBins);
    TH2D* otherhist = new TH2D("OtherDistribution", "Hadron Counts; phi; pT (GeV); Hadrons Out",
                               phis.size()-1, &phis[0],
                               hadDatav2Hist->GetNbinsX(), hadDatav2Hist->GetXaxis()->GetXbins()->GetArray());
    TH1D* etaHist = new TH1D("EtaDist", "Photon Counts; eta; Photons Out", 20, -2.5, 2.5);
    TH1D* phiHist = new TH1D("PhiDist", "Photon Counts; #phi; Photons Out", 20, -3.14159, 2.*3.14159);

    QVecV2* photonQv2 = new QVecV2("photonhist", "Direct Photon v_{2}; p_{T} [GeV]; v_{2}",
                             nv2bins, v2bins);

    cout << "Analysis graphs made" << endl;

    // Output ROOT file.
    TFile* outfile = TFile::Open((outDir + "analysis.root").c_str(), "RECREATE");
    
    // Get list of directories (each representing one analysis point).
    std::vector<std::string> directories = getDirs(pointsDir);
    
    int eventCount = 0;
    double etaCut = 0.5;
    double hadEtaCut = 0.35;
    double hardBase = 1.0/41.9734;
    double hardFactor = 3.*0.115526/41.9734; //prompt photon over inelastic xsec
    double issFactor = 1.0/JetScapeXML::Instance()->GetElementDouble({"SoftParticlization", "iSS", "number_of_repeated_sampling"});
    double thermalFactor = 1.0/emissionRate;
    double jetFactor = 1./100.;
    
    // Ncoll calculations
    double averageNcoll = 0.0;
    vector<double> Ncolls = {};

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
        Ncolls.push_back(Ncoll);

        double fileEvents = reader->TotalEventCount();
        cout << "Events in file: " << fileEvents << endl;
        jetFactor = 1.0/fileEvents;

        double sigmaGen = reader->GetSigmaGen();
        //hardFactor = hardBase*sigmaGen;
        cout << "SigmaGen: " << sigmaGen << endl;
        cout << "Hard factor: " << hardFactor << endl;
        
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
                
                if((particle.get()->pstat() == 11 || particle.get()->pstat() == 27) && std::fabs(particle.get()->eta()) < hadEtaCut && pythia.particleData.charge(particle.get()->pid()) < 0) //smash tag
                {
                    if(reader->GetCurrentEvent() > 1) continue;
                    weight = issFactor; //iSS resampler count
                    otherhist->Fill(delphi, particle.get()->pt(), weight);
                }
                
                if(std::fabs(particle.get()->eta()) > etaCut)
                    continue;
                
                // Check charge using the internal Pythia particleData.
                if(particle.get()->pstat() == 11 || particle.get()->pstat() == 27)
                {
                    if(std::fabs(particle.get()->eta()) < 0.5)
                    {
                        weight = issFactor; //iSS resampler count
                        if(particle.get()->pid() == 211) hadSpecHist->Fill(particle.get()->pt(), weight);
                    }
                    continue;
                }

                //Adding shower hadrons to v2
                if(particle.get()->pstat() == 0)
                {
                    if(pythia.particleData.charge(particle.get()->pid()) != 0)
                    {
                        weight = jetFactor;
                        if(particle.get()->pid() == 211) hadSpecHist->Fill(particle.get()->pt(), weight);
                        otherhist->Fill(delphi, particle.get()->pt(), weight);
                    }
                    continue;
                }
                if(particle.get()->pstat() == -1)
                {
                    if(pythia.particleData.charge(particle.get()->pid()) != 0)
                    {
                        weight = -1.*jetFactor;
                        if(particle.get()->pid() == 211) hadSpecHist->Fill(particle.get()->pt(), weight);
                        otherhist->Fill(delphi, particle.get()->pt(), weight);
                    }
                    continue;
                }
                
                if(particle.get()->pid() == 22)
                {
                    //skipping absorbed ones in abs runs
                    if(abson == true and particle.get()->pstat() < 0)
                        continue;

                    if(abs(particle.get()->pstat()) == 22)
                    {
                        weight = hardFactor;
                        hardHist->Fill(particle.get()->pt(), weight);
                        Hardspread->Fill(delphi, particle.get()->pt(), weight);
                    }

                    if(abs(particle.get()->pstat()) == 23)
                    {
                        weight = hardFactor;
                        showerHist->Fill(particle.get()->pt(), weight);
                        hardHist->Fill(particle.get()->pt(), weight);
                        Hardspread->Fill(delphi, particle.get()->pt(), weight);
                    }

                    if((particle.get()->pstat()) == 24){
                        QGPHist->Fill(particle.get()->pt(), weight);
                        QGPspread->Fill(delphi, particle.get()->pt(), weight);
                    }
                    
                    if(abs(particle.get()->pstat()) == 25){
                        HRGHist->Fill(particle.get()->pt(), weight);
                        HRGspread->Fill(delphi, particle.get()->pt(), weight);
                    }

                    //if(particle.get()->pstat() == 27) //SMASH photons maybe include
                    //    HRGHist->Fill(particle.get()->pt(), 100.);
                    totalhist->Fill(delphi, particle.get()->pt(), weight);
                    exthist->Fill(delphi, particle.get()->pt(), weight);
                    pThist->Fill(particle.get()->pt(), weight);
                    extendedHist->Fill(particle.get()->pt(), weight);
                    photonQv2->Fill(particle.get()->pt(), weight, delphi);

                    if(particle.get()->pt() > 1.0) {
                        etaHist->Fill(particle.get()->eta());
                        phiHist->Fill(particle.get()->phi());
                    }
                }
            } // end particle loop
        } // end event loop
        reader->Close();
    } // end directory loop

    // Scale and write histograms.
    if(eventCount > 0)
    {
        // moving to output dir
        chdir(outDir.c_str());

        // Final Ncoll
        averageNcoll /= directories.size();
        //cout << "Average Ncoll: " << averageNcoll << endl;
        ppbase->Scale(averageNcoll); //value from https://journals.aps.org/prc/abstract/10.1103/PhysRevC.93.044906 Table I

        //plot info
        TPaveText *pt = new TPaveText(.15, .1, .4, .3, "NDC");
        pt->AddText("Au-Au 200 GeV");
        if(central) pt->AddText("Centrality 0-20%");
        else pt->AddText("Centrality 20-40%");
        pt->AddText("-0.35 < #eta_{#gamma} < 0.35");
        pt->SetFillColor(kWhite);

        totalhist->Scale(1.0/eventCount);
        totalhist->Write();
        
        otherhist->Scale(1.0/eventCount);
        otherhist->Write();
        
        extendedHist->Scale(1.0/(2.0*M_PI*2.0*etaCut*eventCount), "width");
        ScaleByCenter(extendedHist);
        
        vector<TGraphAsymmErrors*> photonGraphs = {};
        QGPHist->SetLineColor(kPink);
        QGPHist->SetFillColorAlpha(kPink, 0.35);
        QGPHist->SetLineWidth(2);
        QGPHist->Scale(1.0/(2.0*M_PI*2.0*etaCut*eventCount), "width");
        ScaleByCenter(QGPHist);
        QGPHist->Write();
        ExtractTemp(QGPHist);
        TGraphAsymmErrors* QGPGraph = new TGraphAsymmErrors(QGPHist);
        QGPGraph->SetDrawOption("l3"); photonGraphs.push_back(QGPGraph);
        
        HRGHist->SetLineColor(kViolet);
        HRGHist->SetFillColorAlpha(kViolet, 0.35);
        HRGHist->SetLineWidth(2);
        HRGHist->Scale(1.0/(2.0*M_PI*2.0*etaCut*eventCount), "width");
        ScaleByCenter(HRGHist);
        HRGHist->Write();
        TGraphAsymmErrors* HRGGraph = new TGraphAsymmErrors(HRGHist);
        HRGGraph->SetDrawOption("l3"); photonGraphs.push_back(HRGGraph);
        
        hardHist->SetLineColor(kGreen);
        hardHist->SetFillColorAlpha(kGreen, 0.35);
        hardHist->SetLineWidth(2);
        hardHist->Scale(1.0/(2.0*M_PI*2.0*etaCut*eventCount), "width");
        ScaleByCenter(hardHist);
        hardHist->Write();
        TGraphAsymmErrors* hardGraph = new TGraphAsymmErrors(hardHist);
        hardGraph->SetDrawOption("l3"); photonGraphs.push_back(hardGraph);
        
        showerHist->SetLineColor(kOrange);
        showerHist->SetFillColorAlpha(kOrange, 0.35);
        showerHist->SetLineWidth(2);
        showerHist->Scale(1.0/(2.0*M_PI*2.0*etaCut*eventCount), "width");
        ScaleByCenter(showerHist);
        showerHist->Write();
        TGraphAsymmErrors* showerGraph = new TGraphAsymmErrors(showerHist);
        showerGraph->SetDrawOption("l3"); //photonGraphs.push_back(showerGraph);
        
        etaHist->Write();
        phiHist->Write();
        
        // Create pT spectrum projection.
        pThist->Scale(1.0/(2.0*M_PI*2.0*etaCut*eventCount), "width");
        pThist->SetFillColorAlpha(kBlue, 0.35);
        pThist->SetLineWidth(2);
        ScaleByCenter(pThist);
        ExtractTemp(pThist);
        TGraphAsymmErrors* totalGraph = new TGraphAsymmErrors(pThist);
        //totalGraph->Write("total");
        
        // Process hadron spectra.
        hadSpecHist->Scale(1.0/(2.0*M_PI*2.0*0.5*directories.size()), "width");
        ScaleByCenter(hadSpecHist);
        hadSpecHist->SetFillColorAlpha(kBlue, 0.35);
        hadSpecHist->SetLineColor(kBlue);
        hadSpecHist->SetLineWidth(3);
        TGraphAsymmErrors* hadSpecGraph = new TGraphAsymmErrors(hadSpecHist);
        hadSpecHist->Write("hadrons");
        
        // Draw hadron spectrum.
        pionGraph->SetTitle("STAR Data; p_{T} [GeV]; 1/(2 #pi N_{ev}) d^{2} N/dp_{T} d#eta [1/GeV]");
        myRatioPlot(pionGraph, hadSpecGraph, "Charged Hadron Spectrum", {}, true, true);
        
        TGraphAsymmErrors* flowHads = MakevnGraph(otherhist, "Hadron v2");
        flowHads->Write();
        flowHads->SetTitle("STAR Data; p_{T} [GeV]; v_{2}");
        myRatioPlot(hadDatav2, flowHads, "Charged Hadron v2");
        
        // Photon v2 then spectrum
        TGraphAsymmErrors* flowPhotons = MakevnGraph(totalhist, "Photon v2");
        flowPhotons->Write();
        myRatioPlot(photonDatav2, flowPhotons, "Photon v2");

        //Extended photon v2
        exthist->Scale(1.0/eventCount, "width");
        TGraphAsymmErrors* extPhotons = MakevnGraph(exthist, "Extended Photon v2");
        extPhotons->Write();
        
        specDataSTAR->SetMarkerStyle(kFullSquare);
        specDataSTAR->SetLineWidth(3);
        specDataSTAR->SetDrawOption("AP");
        photonGraphs.push_back(specDataSTAR);
        photonGraphs.push_back(ppbase);
        totalGraph->SetTitle("X-SCAPE");
        specData->SetTitle("PHENIX Data; p_{T} [GeV]; 1/(2 #pi N_{ev}) d^{2} N/p_{T}dp_{T} dy [GeV^{-2}]");
        myRatioPlot(specData, totalGraph, "Direct Photon p_{T} Spectrum", photonGraphs, true, true, pt);

        //other v2s
        TGraphAsymmErrors* QGPv2 = MakevnGraph(QGPspread, "QGP v2");
        QGPv2->Write("QGPv2");
        
        TGraphAsymmErrors* HRGv2 = MakevnGraph(HRGspread, "HRG v2");
        HRGv2->Write("HRGv2");
        
        TGraphAsymmErrors* Hardv2 = MakevnGraph(Hardspread, "Hard v2");
        Hardv2->Write("Hardv2");

        //saving for other usage
        TGraphAsymmErrors* expandedTotalGraph = new TGraphAsymmErrors(extendedHist);
        expandedTotalGraph->Write("total");

        //test new v2 method
        auto QVecv2 = photonQv2->GetV2Graph();
        QVecv2->Write("photon SPv2");

        JSINFO << "Ncoll Mean: " << Mean(Ncolls);
        JSINFO << "Ncoll StdDev: " << StdDev(Ncolls);
    }
    else
    {
        std::cerr << "No events processed!" << std::endl;
    }
    
    outfile->Close();
    std::cout << "Analysis complete. Output written to " << outDir + "analysis.root" << std::endl;
    std::cout << "Total Events: " << eventCount << std::endl;
    
    return 0;
}