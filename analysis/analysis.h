#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <dirent.h>
#include <cmath>
#include <memory>
#include <cstdlib>  // for exit()
#include <numeric>
#include <filesystem>
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

int GetNcoll(string dirName){
    string slurm_path = dirName + "/slurm.out";
    std::ifstream infile(slurm_path);
    if (!infile) return -1;
    std::string line;
    while (std::getline(infile, line)) {
        // Try to find 'ncoll' with a value
        if (line.find("npart") != std::string::npos) {
            std::getline(infile, line); //incrementing to next line after trento header
            std::istringstream iss(line);
            string a, b, c, d, e;
            iss >> a >> b >> c >> d >> e;
            return stoi(e);
        }
    }
    return -1;
}

// Helper: list directory names inside a given folder.
std::vector<std::string> getDirs(const std::string& basePath)
{
    std::vector<std::string> directories;
    DIR *dir = opendir(basePath.c_str());
    if (!dir) {
        std::cerr << "Error: Cannot open directory " << basePath << std::endl;
        return directories;
    }
    struct dirent *entry;
    while ((entry = readdir(dir)) != nullptr) {
        // Skip current and parent directory entries
        if (std::string(entry->d_name) == "." || std::string(entry->d_name) == "..")
            continue;
        // Check if this entry is a directory.
        if (entry->d_type == DT_DIR) {
            directories.push_back(std::string(entry->d_name));
        }
    }
    closedir(dir);
    
    // Sort the directories vector numerically.
    std::sort(directories.begin(), directories.end(), [](const std::string& a, const std::string& b) {
        // Convert the directory names to integers.
        // It is assumed that all directory names are valid numeric strings.
        return std::atoi(a.c_str()) < std::atoi(b.c_str());
    });
    
    return directories;
}

//compact way to grab item from ROOT file
TObject* getFileObject(TFile* rootfile, string dir, string objname){
    TDirectory* directory = (TDirectory*)rootfile->Get(dir.c_str());
    return directory->Get(objname.c_str());
}

TH1D* replicateHistogram(TFile* rootfile, string dir, string objname){
    TH1D* dataHist = (TH1D*)getFileObject(rootfile, dir, objname);
    TH1D* newHist = (TH1D*)dataHist->Clone();
    newHist->Reset();
    return newHist;
}

// Ratio plot generator
void myRatioPlot(TGraphAsymmErrors* dataGraph, TGraphAsymmErrors* predictionHist, string title, vector<TGraphAsymmErrors*> others = {}, bool xlog = false, bool ylog = false, TPaveText* pt = nullptr){   
    //values for plots
    int bins = dataGraph->GetN();
    double data[bins], prediction[bins], predictionerrors[bins], xcoords[bins], dataErrors[bins], binWidths[bins], asym[bins], ratioError[bins], asymErrorsHigh[bins], asymErrorsLow[bins];
    double predictionwidths[bins] = {0.0};

    //title shenanigans
    string tempTitle = dataGraph->GetTitle();
    dataGraph->SetTitle(title.c_str());

    //cycling through bins, note first one is an overflow bin so it is skipped
    for(int i = 0; i < bins; i++){
        //reading values from histograms
        dataGraph->GetPoint(i,xcoords[i],data[i]);
        predictionHist->GetPoint(i,xcoords[i],prediction[i]);
        dataErrors[i] = abs(dataGraph->GetErrorY(i));
        binWidths[i] = abs(dataGraph->GetErrorX(i));
        predictionerrors[i] = abs(predictionHist->GetErrorY(i));
        asym[i] = abs(prediction[i]/data[i]);
        asymErrorsHigh[i] = 1 + abs(dataErrors[i]/data[i]);
        asymErrorsLow[i] = 1 - abs(dataErrors[i]/data[i]);
        ratioError[i] = abs(predictionerrors[i]/data[i]);

        //debugging
        //cout << xcoords[i] << " " << data[i] << " " << dataErrors[i] << " " << asym[i] << " " << ratioError[i] << endl;
        
        //correcting for negatives when doing log
        //if(data[i] < 0) data[i] = 0;
        //if(dataErrors[i] > data[i]) dataErrors[i] = data[i]*0.99;
    }

    //Drawing main plot
    TCanvas* c = new TCanvas("c1","c1",1000,800);
    TPad* upper = new TPad("plot","plot",0,0.4,1,1);
    upper->SetBottomMargin(0);
    upper->Draw();
    upper->cd();
    if(ylog) upper->SetLogy();
    if(xlog) upper->SetLogx();

    //Data graph
    TGraphAsymmErrors* dataPlot = (TGraphAsymmErrors*)dataGraph->Clone();
    dataPlot->SetMarkerStyle(kFullDotLarge);
    dataPlot->SetLineWidth(3);
    if(!ylog) dataPlot->GetHistogram()->GetYaxis()->SetRangeUser(-0.08,0.35);

    //prediction graph
    TGraphAsymmErrors* predictionPlot = (TGraphAsymmErrors*)predictionHist->Clone();
    predictionPlot->SetLineColor(kBlue);
    predictionPlot->SetFillColorAlpha(kBlue, 0.3);
    predictionPlot->SetLineWidth(3);

    //layered graph
    dataPlot->Draw("AP");
    if(pt != nullptr) pt->Draw("SAME NDC");
    dataPlot->GetXaxis()->SetLimits(xcoords[0]-binWidths[0],xcoords[bins-1]+binWidths[bins-1]);
    dataPlot->GetYaxis()->SetLabelSize(0.04);
    dataPlot->GetYaxis()->SetTitleSize(0.07);
    dataPlot->GetYaxis()->SetTitleOffset(0.7);
    predictionPlot->Draw("l3 SAME");

    //ratio plot
    TGraphAsymmErrors* comparisonPlot = new TGraphAsymmErrors(bins, xcoords, asym, binWidths, binWidths, ratioError, ratioError);
    comparisonPlot->SetTitle("");
    comparisonPlot->GetXaxis()->SetLimits(xcoords[0]-binWidths[0],xcoords[bins-1]+binWidths[bins-1]);
    comparisonPlot->GetXaxis()->SetTitle(dataPlot->GetHistogram()->GetXaxis()->GetTitle());
    comparisonPlot->GetXaxis()->SetTitleSize(0.08);
    comparisonPlot->GetXaxis()->SetLabelSize(0.06);
    comparisonPlot->GetYaxis()->SetTitle("X-SCAPE/data");
    comparisonPlot->GetYaxis()->CenterTitle(true);
    comparisonPlot->GetYaxis()->SetTitleSize(0.08);
    comparisonPlot->GetYaxis()->SetTitleOffset(0.7);
    comparisonPlot->GetYaxis()->SetLabelSize(0.06);
    comparisonPlot->GetHistogram()->GetYaxis()->SetRangeUser(0,2.1);  
    comparisonPlot->SetLineColor(kBlue);
    comparisonPlot->SetFillColorAlpha(kBlue, 0.3);
    comparisonPlot->SetLineWidth(3);

    //Legend
    TLegend leg(.6,.50,.9,.9);
    leg.AddEntry(predictionPlot,predictionHist->GetTitle(),"l");
    leg.AddEntry(dataPlot,tempTitle.c_str(),"ep");
    leg.SetTextSize(0.05);
    for(int i=0; i<others.size(); i++){
        string name = others[i]->GetTitle();
        if(name.find("data") != std::string::npos){
            others[i]->SetLineColor(kGray+2);
            others[i]->SetMarkerColor(kGray+2);
            others[i]->Draw("p SAME");
            leg.AddEntry(others[i],others[i]->GetHistogram()->GetTitle(),"ep");
        }
        else{
            others[i]->Draw("l3 SAME");
            leg.AddEntry(others[i],others[i]->GetHistogram()->GetTitle(),"l");
        }
    }
    leg.DrawClone("Same");

    //drawing ratio plot
    c->cd();
    TPad* lower = new TPad("plot","plot",0,0,1,0.4);
    lower->SetTopMargin(0);
    lower->SetBottomMargin(0.2);
    lower->Draw();
    lower->cd();
    if(xlog) lower->SetLogx();
	comparisonPlot->Draw("");
    TLine* line = new TLine(xcoords[0],1,xcoords[bins-1]+binWidths[bins-1],1);
    line->SetLineStyle(9);
    line->Draw("SAME");
    TGraph* asymLow = new TGraph(bins, xcoords, asymErrorsLow);
    asymLow->SetLineStyle(9);
    asymLow->Draw("SAME");
    TGraph* asymHigh = new TGraph(bins, xcoords, asymErrorsHigh);
    asymHigh->SetLineStyle(9);
    asymHigh->Draw("SAME");

    string filename = title + ".pdf";
    c->Print(filename.c_str());
    c->Close();
}

double GetEpAngle(vector<shared_ptr<Jetscape::Hadron>> hadrons){
  double Qx = 0.0;
  double Qy = 0.0;

  for (auto hadron : hadrons) {
    if ((hadron.get()->pstat() == 11 || hadron.get()->pstat() == 27) && hadron.get()->pt() > 0.5) {
      Qx += hadron.get()->pt() * std::cos(2.0 * hadron.get()->phi());
      Qy += hadron.get()->pt() * std::sin(2.0 * hadron.get()->phi());
    }
  }

  return std::atan2(Qy, Qx) / 2.0;
}

double Getv2(TH2D* hist, int bin)
{
    double integral = 0.0;
    double denom = 0.0;
    int nbins = hist->GetNbinsX();

    // Loop over all bins (bins in ROOT are 1-indexed)
    for (int i = 1; i <= nbins; i++) {
        // Retrieve bin content, center, and width from X-axis
        double content = hist->GetBinContent(i,bin);
        double binCenter = hist->GetXaxis()->GetBinCenter(i);
        double binWidth = hist->GetXaxis()->GetBinWidth(i);
        integral += content * std::cos(2.0 * binCenter) * binWidth;
        denom += content * binWidth;
    }

    if (denom <= 0.0)
        return 0.0;
    
    return integral / denom;
}

// Function to extract a temperature by fitting the histogram with an exponential.
// hist: Pointer to a TH1D histogram.
void ExtractTemp(TH1D* hist) {
    // Define the fit range using the x-axis minimum and maximum.
    double xmin = hist->GetXaxis()->GetXmin();
    double xmax = hist->GetXaxis()->GetXmax();

    // Create TF1 for the exponential fit.
    TF1* QGPfit = new TF1("thisfit", "[0]*exp(-x/[1])", 1.1, 4.1); //limits compatible with data
    QGPfit->SetParNames("A", "T");
    // Set the initial guess for the temperature parameter.
    QGPfit->SetParLimits(1, 0.1, 1);
    QGPfit->SetLineColor(kBlack);

    // Fit the histogram with the TF1 function.
    hist->Fit("thisfit", "Q R"); // quiet fit

    // Print the fitted temperature parameter multiplied by 1000 (converted to MeV).
    std::cout << "Fit temperature is " << static_cast<int>(QGPfit->GetParameter(1) * 1000.0)
              << " MeV" << std::endl;

    // Note: Consider proper memory management for QGPfit if necessary.
}

std::pair<double, double> GetvnFit(TH2D* hist, int bin, int n) {
    // Project the histogram on the X axis for the specified y-bin range.
    TH1D* phiHist = hist->ProjectionX("low pT distribution", bin);
    
    // Calculate the amplitude coefficient (integral) using bin contents and widths.
    double integral = 0.0;
    int nbinsx = phiHist->GetNbinsX();
    for (int i = 1; i <= nbinsx; ++i) {
        integral += phiHist->GetBinContent(i) * phiHist->GetBinWidth(i);
    }
    
    // If the integral is non-positive, no meaningful fit can be done.
    if (integral <= 0.0) {
        delete phiHist;
        return std::make_pair(0.0, 0.0);
    }
    
    // Create a TF1 fit function: f(x) = [0]*(1 + ([1]*2*cos(2*x)))
    TF1* thisfit = new TF1("thisfit", "[0]*(1+([1]*2*cos([2]*x)))", 0, 2 * TMath::Pi());
    thisfit->SetParameter(0, integral / (2 * TMath::Pi()));
    thisfit->SetParameter(1, 0.0);
    thisfit->FixParameter(2, n);
    thisfit->SetParNames("A", "v2", "n");
    
    // Fit the projection histogram quietly.
    phiHist->Fit(thisfit, "Q");
    
    // Extract the fitted v2 parameter.
    double v2 = thisfit->GetParameter(1);
    double errorV2 = thisfit->GetParError(1);

    TCanvas* c = new TCanvas();
    phiHist->SetStats(0);
    phiHist->Draw();
    string title = hist->GetTitle();
    string name = "fits/"+title+to_string(bin)+".pdf";
    c->Print(name.c_str());
    
    // Cleanup allocated memory.
    delete phiHist;
    delete thisfit;
    
    return std::make_pair(v2, errorV2);
}

// Function to calculate v2 and its error
std::pair<double, double> GetvnWithError(TH2D* hist, int bin, int n)
{
    double integral = 0.0;
    double denom = 0.0;
    double errorIntegral2 = 0.0;
    double errorDenom2 = 0.0;
    int nbins = hist->GetNbinsX();

    // Loop over all bins (bins in ROOT are 1-indexed)
    for (int i = 1; i <= nbins; i++) {
        // Retrieve bin content, center, and width from X-axis
        double content = hist->GetBinContent(i, bin);
        double error = hist->GetBinError(i, bin); // Get the error for the bin content
        double binCenter = hist->GetXaxis()->GetBinCenter(i);
        double binWidth = hist->GetXaxis()->GetBinWidth(i);
        
        integral += content * std::cos(double(n) * binCenter) * binWidth;
        denom += content * binWidth;
        
        // Error propagation for the integral
        errorIntegral2 += std::pow(error * std::cos(double(n) * binCenter) * binWidth, 2);
        errorDenom2 += std::pow(error * binWidth, 2);
    }

    if (denom <= 0.0)
        return std::make_pair(0.0, 0.0);

    double v2 = integral / denom;
    double errorV2 = std::sqrt((errorIntegral2 / std::pow(denom, 2))
                     + (errorDenom2 * std::pow(integral / (denom * denom), 2))); // Error propagation for v2

    return std::make_pair(v2, errorV2);
}

// New function: conversion of the Python Makev2Graph to C++
TGraphAsymmErrors* MakevnGraph(TH2D* hist, std::string name, int n = 2, bool fit = false)
{
    int nBins = hist->GetNbinsY();
    vector<Double_t> pTs;
    vector<Double_t> v2s;
    vector<Double_t> errors;
    vector<Double_t> widths;
    
    // Loop over all bins in the Y-axis (ROOT bins are 1-indexed)
    for (int i = 1; i <= nBins; ++i)
    {
        std::pair<double, double> v2;
        if(fit)
            v2 = GetvnFit(hist, i, n);
        else
            v2 = GetvnWithError(hist, i, n);

        double pT = hist->GetYaxis()->GetBinCenter(i);
        double width = hist->GetYaxis()->GetBinWidth(i)/2.0;
        pTs.push_back(pT);
        v2s.push_back(v2.first);
        errors.push_back(v2.second);
        widths.push_back(width);
    }
    
    // Create a TGraph using the data arrays from std::vector
    TGraphAsymmErrors* graph = new TGraphAsymmErrors(nBins, pTs.data(), v2s.data(), widths.data(), widths.data(), errors.data(), errors.data());
    graph->SetNameTitle(name.c_str(), (name + ";p_{T} (GeV); v_{"+to_string(n)+"}").c_str());
    graph->SetLineColor(kBlue);
    graph->SetFillColorAlpha(kBlue, 0.35);
    graph->SetLineWidth(2);
    return graph;
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

//spreading out pythia photons
void smoothFill(TH2D* hist, double pT, double weight){
    int nBinsPhi = hist->GetNbinsX();
    double delphi = 2.*M_PI/nBinsPhi;
    double newWeight = weight / nBinsPhi;

    double thisPhi = delphi/2.;
    for(int i = 1; i <= nBinsPhi; i++){
        hist->Fill(thisPhi, pT, newWeight);
        thisPhi += delphi;
    }
}

double Mean(vector<double> v){
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / v.size();

    return mean;
}

double StdDev(vector<double> v){
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / v.size();

    double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
    double stdev = std::sqrt(sq_sum / v.size() - mean * mean);

    return stdev;
}