#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <dirent.h>
#include <cmath>
#include <memory>
#include <cstdlib>  // for exit()
#include <numeric>
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

class QVecV2{
    public:
        QVecV2();

        QVecV2(const char *name, const char *title, Int_t nbinsx, const Double_t *xbins){
            QvecHist = new TH1D(name, title, nbinsx, xbins);
            QvecBase = new TH1D(name, title, nbinsx, xbins);
        };

        void Fill(double x, double weight, double delphi){
            QvecHist->Fill(x, weight*cos(2.0*delphi));
            QvecBase->Fill(x, weight);
        };

        TGraphAsymmErrors* GetV2Graph(){
            QvecHist->Divide(QvecBase);
            return new TGraphAsymmErrors(QvecHist);
        }

    private:
        TH1D* QvecHist;
        TH1D* QvecBase;
};