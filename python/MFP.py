import ROOT
import math
from array import array
from functions import *

MFPfile = ROOT.TFile("plots/MFP-old/MFP.root","r")

MPHist = ROOT.TH1D(MFPfile.Get("MPHist"))
totalHist = ROOT.TH1D(MFPfile.Get("totalHist"))
outHist = ROOT.TH1D(MFPfile.Get("outHist"))

MPHist.Divide(outHist)
outHist.Divide(totalHist)

c1 = ROOT.TCanvas()
outHist.SetTitle("Ratio of photons out")
outHist.GetXaxis().SetTitle("N_{out}/N_{total}")
outGraph = ROOT.TGraphErrors(outHist)
outGraph.Draw("")
c1.Print("plots/MFP-old/Ratio-out.png")
c1.Close()

c2 = ROOT.TCanvas()
MPHist.SetTitle("Mean Medium Path")
MPHist.Draw("L3")
c2.Print("plots/MFP-old/Average-path.png")
c2.Close()

MFPlist = []
MFPerrList = []
pList = []
pErrList = []


for i in range(1,MPHist.GetNbinsX()-4):
    try:
        L = MPHist.GetBinContent(i)
        delL = MPHist.GetBinError(i)
        R = outHist.GetBinContent(i)
        delR = outHist.GetBinError(i)
        print("Momentum:",outHist.GetBinCenter(i))
        print("Ratio out:",R)
        print("Average medium path:",L)
        L0 = -L / math.log(R)
        delL0 = math.sqrt(pow(delL/math.log(R),2) + pow(delR*L/(R*math.log(R)*math.log(R)),2))
    except:
        L0 = 0
        delL0 = 0
    print("MFP:",L0,"+/-",delL0)
    MFPlist.append(L0)
    MFPerrList.append(delL0)
    pList.append(outHist.GetBinCenter(i))
    pErrList.append(outHist.GetBinWidth(i)/2.)

n = len(MFPlist)
x = array('d',pList)
xerr = array('d',pErrList)
y = array('d',MFPlist)
yerr = array('d',MFPerrList)

gr = ROOT.TGraphErrors(n, x, y, xerr, yerr)
gr.SetTitle("Mean Free Path;p_{T} [GeV];L_{0} [fm]")
gr.SetMarkerStyle(21)
gr.SetMarkerColor(4)
gr.SetLineColor(4)
gr.SetLineWidth(3)
gr.SetFillColorAlpha(4,0.3)

c3 = ROOT.TCanvas("c1", "Canvas", 800, 600)
gr.Draw("AL3") # A = Axis, P = Points
c3.Update()
c3.Print("plots/MFP-old/MFP.pdf")