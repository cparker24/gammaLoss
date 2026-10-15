import ROOT
from fontTools.misc.testTools import getXML

def GetXSCAPEGraph(filepath, color, title):
    plotFile = ROOT.TFile.Open(filepath)
    keys = plotFile.GetListOfKeys()
    firstKey = keys.At(1)

    plot = firstKey.ReadObj()
    plot.SetTitle(title)
    plot.SetLineColor(color)
    plot.SetFillColorAlpha(color,0.35)
    plot.SetLineWidth(2)
    return ROOT.TGraphErrors(plot)

minBiasGraph = GetXSCAPEGraph("/scratch/user/cameron.parker/projects/runs/pp/pp200/photon_spectrum.root", ROOT.kGreen-2, "PYTHIA8 Min-bias")
promptGraph = GetXSCAPEGraph("/scratch/user/cameron.parker/projects/runs/pp/pp200prompt/photon_spectrum.root", ROOT.kBlue-2, "PYTHIA8 Prompt only")
#nuclearGraph = GetXSCAPEGraph("/scratch/user/cameron.parker/projects/runs/pp/pp2760nuclearprompt/photon_spectrum.root", ROOT.kRed-2, "PYTHIA8 nPDF Prompt only")
gunGraph = GetXSCAPEGraph("/scratch/user/cameron.parker/projects/runs/pp/pp200gunnew/photon_spectrum.root", ROOT.kRed, "pQCD Gun")

rootDataFile = ROOT.TFile.Open("/scratch/user/cameron.parker/projects/C-SCAPE-COMP-gamma/external_packages/gammaLoss/data/pp200PhotonSpectrum.root")
expSigmaGraph = rootDataFile.Get("Figure 2").Get("Graph1D_y2")
expSigmaGraph.SetTitle("PHENIX")
expSigmaGraph.SetMarkerStyle(24)
expSigmaGraph.Scale(1.0e9)

fa1 = ROOT.TF1("parameterization","2*6495*200/(((200/3450)+x)^5)",.1,15)
fa1.SetTitle("New Parameterization")
fa1.SetLineStyle(8)
fa1.SetLineColor(ROOT.kBlack)

fa2 = ROOT.TF1("parameterization","6495*200/((x)^5)",.1,15)
fa2.SetTitle("Old Parameterization")
fa2.SetLineStyle(9)
fa2.SetLineColor(ROOT.kGray+2)

#plot details
pt = ROOT.TPaveText(.65, .65, .85, .85, "NDC")
pt.SetFillColor(ROOT.kWhite)
pt.AddText("p-p @ 200 GeV")
pt.AddText("-0.35 < #eta_{#gamma} < 0.35")

c2 = ROOT.TCanvas("c2", "Photon sigma comparison", 800, 600)
c2.SetLogy()
#c2.SetLogx()

mg2 = ROOT.TMultiGraph()
mg2.SetTitle(" ;p_{T} [GeV]; E d^{3}#sigma/dp^{3} [pb/GeV^{2}]")
mg2.Add(expSigmaGraph)
mg2.Add(minBiasGraph,"l3")
#mg2.Add(nuclearGraph,"l3")
mg2.Add(promptGraph,"l3")
mg2.Add(gunGraph,"l3")
mg2.SetMinimum(1e-1)
mg2.Draw("AP")
fa1.Draw("SAME")
fa2.Draw("SAME")

legend = ROOT.TLegend()
legend.AddEntry(expSigmaGraph)
legend.AddEntry(minBiasGraph, "", "lf")
#legend.AddEntry(nuclearGraph, "", "lf")
legend.AddEntry(promptGraph, "", "lf")
legend.AddEntry(gunGraph, "", "lf")
legend.AddEntry(fa1)
legend.AddEntry(fa2)
legend.Draw("SAME")

pt.Draw("SAME")

c2.SaveAs("plots/others/pp-200-comparison.pdf")