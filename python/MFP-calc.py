import ROOT
import math
from array import array
from functions import *

# setting directory for analysis
analysisDir = sys.argv[1]
directories = getDirs(analysisDir)

# getting initial histograms
MFP0file = ROOT.TFile(analysisDir+"points/0/MFP.root","r")
MPHist = ROOT.TH1D(MFP0file.Get("MPHist"))
totalHist = ROOT.TH1D(MFP0file.Get("totalHist"))
outHist = ROOT.TH1D(MFP0file.Get("outHist"))

# loop over files
for i,directory in enumerate(directories):
    filename = analysisDir+"points/"+directory+"/MFP.root"
    print("Reading filename:",filename)
    MFPfile = ROOT.TFile(filename,"r")

    if i == 0:
        continue
    else:
        tempMPHist = ROOT.TH1D(MFPfile.Get("MPHist"))
        temptotalHist = ROOT.TH1D(MFPfile.Get("totalHist"))
        tempoutHist = ROOT.TH1D(MFPfile.Get("outHist"))
        MPHist.Add(tempMPHist)
        totalHist.Add(temptotalHist)
        outHist.Add(tempoutHist)

outHist.Print()
totalHist.Print()

MPHist.Divide(outHist)
outHist.Divide(totalHist)

MPHist.Sumw2()
totalHist.Sumw2()
outHist.Sumw2()

c1 = ROOT.TCanvas()
outHist.SetTitle("Ratio of Photons Out;p_{T} [GeV];N_{out}/N_{total}")
outHist.Draw()
outHist.SetStats(0)
c1.Print("plots/MFP/Ratio-out.png")
c1.Close()

c2 = ROOT.TCanvas()
MPHist.SetTitle("Mean Medium Path")
MPHist.SetMinimum(0)
MPHist.Draw()
MPHist.SetStats(0)
c2.Print("plots/MFP/Average-path.png")
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
gr.SetTitle("Effective Mean Free Path;p_{T} [GeV];<L_{0}> [fm]")
gr.SetMarkerStyle(21)
gr.SetMarkerColor(4)
gr.SetLineColor(4)
gr.SetLineWidth(3)
gr.SetFillColorAlpha(4,0.3)
gr.SetMinimum(0)

c3 = ROOT.TCanvas("c1", "Canvas", 800, 600)
gr.Draw("AL3") # A = Axis, P = Points
c3.Update()
c3.Print("plots/MFP/MFP.pdf")