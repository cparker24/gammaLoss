import ROOT
from fontTools.misc.testTools import getXML

def GetXSCAPEGraph(filepath, color, title):
    plotFile = ROOT.TFile.Open(filepath)
    plot = plotFile.Get("xsec")
    plot.SetTitle(title)
    plot.SetLineColor(color)
    plot.SetFillColorAlpha(color,0.35)
    plot.SetLineWidth(2)
    return ROOT.TGraphErrors(plot)

minBiasGraph = GetXSCAPEGraph("/scratch/user/cameron.parker/projects/runs/pp/pp2760baseline/photon_spectrum.root", ROOT.kGreen-2, "PYTHIA8 Min-bias")
promptGraph = GetXSCAPEGraph("/scratch/user/cameron.parker/projects/runs/pp/pp2760prompt/photon_spectrum.root", ROOT.kBlue-2, "PYTHIA8 Prompt only")
nuclearGraph = GetXSCAPEGraph("/scratch/user/cameron.parker/projects/runs/pp/pp2760nuclearprompt/photon_spectrum.root", ROOT.kRed-2, "PYTHIA8 nPDF Prompt only")
gunGraph = GetXSCAPEGraph("/scratch/user/cameron.parker/projects/runs/pp/pp2760gunnew/photon_spectrum.root", ROOT.kRed, "pQCD Gun")

rootDataFile = ROOT.TFile.Open("/scratch/user/cameron.parker/projects/C-SCAPE-COMP-gamma/external_packages/gammaLoss/data/pp2760PhotonSpectrum.root")
expSigmaGraph = rootDataFile.Get("Table 5").Get("Graph1D_y1")
expSigmaGraph.MovePoints(-10,0)
expSigmaGraph.SetTitle("ALICE")
expSigmaGraph.SetMarkerStyle(24)

sigmabins = ROOT.std.vector('double')()
for v in [0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.8, 3.25, 3.75, 4.5, 5.5, 9.0]: sigmabins.push_back(v)
sigmaLimits = ROOT.std.vector('double')()
for v in [7.68e9, 9.22e8, 3.82e8, 1.46e8, 4.37e7, 1.6e7, 1.1e7, 6.72e6, 5.25e6, 3.04e6, 1.26e6, 8.21e5, 5.98e5, 1.79e5, 9.85e4, 3.69e4, 2.92e3]: sigmaLimits.push_back(v)
sigmaLimitsLowers = ROOT.std.vector('double')()
for v in [7.68e9, 9.22e8, 3.82e8, 1.46e8, 4.37e7, 1.6e7, 1.1e7, 6.72e6, 5.25e6, 3.04e6, 1.26e6, 8.21e5, 5.98e5, 1.79e5, 9.85e4, 3.69e4, 2.92e3]: sigmaLimitsLowers.push_back(v)
sigmaLimitsUppers = ROOT.std.vector('double')()
for v in [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]: sigmaLimitsUppers.push_back(v)
sigmaLimitsWidths = ROOT.std.vector('double')()
for v in [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.25, 0.25, 0.5, 0.5, 1.0]: sigmaLimitsWidths.push_back(v)
limits = ROOT.TGraphAsymmErrors( sigmabins.size(), sigmabins.data(), sigmaLimits.data(), sigmaLimitsWidths.data(), sigmaLimitsWidths.data(), sigmaLimitsLowers.data(), sigmaLimitsUppers.data() )

fa1 = ROOT.TF1("parameterization","2*6495*2760/((0.8+x)^5)",.1,15)
fa1.SetTitle("New Parameterization")
fa1.SetLineStyle(8)
fa1.SetLineColor(ROOT.kBlack)

fa2 = ROOT.TF1("parameterization","6495*2760/((x)^5)",.1,15)
fa2.SetTitle("Old Parameterization")
fa2.SetLineStyle(9)
fa2.SetLineColor(ROOT.kGray+2)

#plot details
pt = ROOT.TPaveText(.65, .65, .85, .85, "NDC")
pt.SetFillColor(ROOT.kWhite)
pt.AddText("p-p @ 2.76 TeV")
pt.AddText("-0.9 < #eta_{#gamma} < 0.9")

c2 = ROOT.TCanvas("c2", "Photon sigma comparison", 800, 600)
c2.SetLogy()
c2.SetLogx()

mg2 = ROOT.TMultiGraph()
mg2.SetTitle(" ;p_{T} [GeV]; E d^{3}#sigma/dp^{3} [pb/GeV^{2}]")
mg2.Add(expSigmaGraph)
mg2.Add(minBiasGraph,"l3")
mg2.Add(nuclearGraph,"l3")
mg2.Add(promptGraph,"l3")
mg2.Add(gunGraph,"l3")
mg2.Add(limits)
mg2.SetMinimum(1e1)
mg2.Draw("AP")
fa1.Draw("SAME")
fa2.Draw("SAME")

legend = ROOT.TLegend()
legend.AddEntry(expSigmaGraph)
legend.AddEntry(minBiasGraph, "", "lf")
legend.AddEntry(nuclearGraph, "", "lf")
legend.AddEntry(promptGraph, "", "lf")
legend.AddEntry(gunGraph, "", "lf")
legend.AddEntry(fa1)
legend.AddEntry(fa2)
legend.Draw("SAME")

pt.Draw("SAME")

c2.SaveAs("plots/others/pp-comparison.pdf")