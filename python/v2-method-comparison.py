import ROOT

rootfile = ROOT.TFile("/scratch/user/cameron.parker/projects/runs/new/PbPb2760/20-40/PCE-narrow/analysis/analysis.root")

fitGraph = rootfile.Get("ExtendedPhoton v2 fit")
integralGraph = rootfile.Get("Extended Photon v2")
spGraph = rootfile.Get("photon SPv2")

canvas = ROOT.TCanvas("comparison", "v2 Comparison", 800, 600)

fitGraph.SetLineColor(ROOT.kRed)
fitGraph.SetLineWidth(2)
fitGraph.SetFillColorAlpha(ROOT.kRed, 0.3)

integralGraph.SetLineColor(ROOT.kBlue)
integralGraph.SetLineWidth(2)
integralGraph.SetFillColorAlpha(ROOT.kBlue, 0.3)

spGraph.SetLineColor(ROOT.kGreen)
spGraph.SetLineWidth(2)
spGraph.SetFillColorAlpha(ROOT.kGreen, 0.3)

mg = ROOT.TMultiGraph()
mg.Add(fitGraph, "l3")
mg.Add(integralGraph, "l3")
mg.Add(spGraph, "l3")

mg.Draw("AP")
mg.SetTitle("Photon v_{2} Method Comparison")
mg.GetXaxis().SetTitle("p_{T} (GeV)")
mg.GetYaxis().SetTitle("v_{2}")

legend = ROOT.TLegend(0.6, 0.7, 0.9, 0.9)
legend.AddEntry(fitGraph, "Fit", "l")
legend.AddEntry(integralGraph, "Integral", "l")
legend.AddEntry(spGraph, "SP", "l")
legend.Draw()

canvas.Update()

canvas.SaveAs("plots/others/v2_comparison.pdf")