import sys
import ROOT
import os
import numpy as np
from ctypes import *
from sys import argv

ROOT.gStyle.SetTitleSize(0.09,"X")
ROOT.gStyle.SetTitleSize(0.068,"Y")
ROOT.gStyle.SetTitleOffset(0.9,"Y")
ROOT.gStyle.SetTitleSize(0.1,"T")
ROOT.gStyle.SetLabelSize(0.1,"X")
ROOT.gStyle.SetLabelSize(0.07,"Y")
ROOT.gStyle.SetPadLeftMargin(0.15)
ROOT.gStyle.SetPadBottomMargin(0.90)
ROOT.gStyle.SetPadTopMargin(0.15)

def histToGraph(hist):
    if(not hist): return hist

    n = hist.GetNbinsX()
    x = []
    xerr = []
    y = []
    yerr = []

    for i in range(n):
        thisy = hist.GetBinContent(i+1)
        thisyerr = hist.GetBinError(i+1)

        if thisy <= 0:
            continue
        
        if thisyerr > thisy:
            thisyerr = 0

        x.append(hist.GetBinCenter(i+1))
        xerr.append(hist.GetBinWidth(i+1)/2.)
        y.append(thisy)
        yerr.append(thisyerr)

    graph = ROOT.TGraphErrors(len(x),np.array(x),np.array(y),np.array(xerr),np.array(yerr))
    graph.SetTitle(hist.GetTitle())
    graph.SetName(hist.GetName())

    return graph

def cleanGraph(inGraph):
    n = inGraph.GetN()
    x = []
    xerr = []
    y = []
    yerr = []

    gmin = 1e30
    gmax = 0
    maxindex = 0

    for i in range(n):
        maxindex = i
        thisx = c_double(0)
        thisy = c_double(0)
        inGraph.GetPoint(i,thisx,thisy)
        thisyerr = inGraph.GetErrorY(i)
        thisx = thisx.value
        thisy = thisy.value


        if thisy <= 0.0:
            thisy = 1e-30
        
        if thisyerr > thisy:
            thisyerr = thisy*0.99

        x.append(thisx)
        xerr.append(0.0)
        y.append(thisy)
        yerr.append(thisyerr)

        if thisy < gmin:
            gmin = abs(thisy - thisyerr)
        if thisy > gmax:
            gmax = thisy + thisyerr

    graph = inGraph.Clone() #ROOT.TGraphErrors(n,np.array(x),np.array(y),np.array(xerr),np.array(yerr))
    graph.SetTitle(inGraph.GetTitle())
    graph.SetName(inGraph.GetName())

    # getting x ranges
    thisx = c_double(0)
    thisy = c_double(0)
    inGraph.GetPoint(0,thisx,thisy)
    xmin = thisx.value - inGraph.GetErrorX(0)
    inGraph.GetPoint(maxindex,thisx,thisy)
    xmax = thisx.value + inGraph.GetErrorX(maxindex)

    return inGraph, gmin, gmax, xmin, xmax

def divideGraph(numerator, denominator):
    bins = int(numerator.GetN())

    xs = []
    ys = []
    xerrs = []
    yerrs = []
    for i in range(bins):
        x = c_double(0)
        y1 = c_double(0)
        y2 = c_double(0)
        numerator.GetPoint(i,x,y1)
        y2 = c_double(denominator.Eval(x,0,""))
        yerr = numerator.GetErrorY(i)
        xerr = numerator.GetErrorX(i)
    
        if y2.value == 0:
            y2 = c_double(1.0)
            y1 = c_double(0)
            yerr = 0

        xs.append(float(x.value))
        ys.append(float(y1.value/y2.value))
        xerrs.append(float(xerr))
        yerrs.append(float(yerr/y2.value))

    return ROOT.TGraphErrors(bins,np.array(xs),np.array(ys),np.array(xerrs),np.array(yerrs))

def makePlot(MCs, plotvars, plotDir, pt):
    # plot settings
    xlog = True if plotvars[2] == "log" else False
    ylog = True if plotvars[3] == "log" else False
    colors = [ROOT.kRed,ROOT.kBlue,ROOT.kBlue+3,ROOT.kAzure+10,ROOT.kBlue-10, ROOT.kMagenta, ROOT.kOrange, ROOT.kGreen]
    datacolors = [ROOT.kBlack, ROOT.kGray+2]
    dataIndex = 0
    title = plotvars[4]

    # Drawing main plot
    upperMG = ROOT.TMultiGraph()
    c = ROOT.TCanvas("c1","c1",1000,700)
    upper = ROOT.TPad("plot","plot",0,0.4,1,1)
    upper.SetBottomMargin(0)
    upper.Draw()
    upper.cd()

    # Data graph
    dataPlot, gmin, gmax, xmin, xmax = cleanGraph(MCs["Default"])
    del MCs["Default"]
    dataPlot.SetTitle("Default")    
    dataPlot.SetLineWidth(3)    
    dataPlot.SetLineColor(colors[0])
    dataPlot.SetFillColorAlpha(colors[0], 0.35)
    upperMG.Add(dataPlot,"l3")
    plotvars[1] = plotvars[1].replace('$', '')
    upperMG.GetYaxis().SetTitle(plotvars[1])
    upperMG.GetHistogram().SetTitle(title)

    if("v_{2}" in title):
        xmax = 6
        upperMG.SetMinimum(-0.07)
        upperMG.SetMaximum(0.4)  

    if("p_{T}" in title):
        xmax = 4  

    # prediction graph formatting
    for i, MC in enumerate(MCs):
        MCs[MC].SetLineWidth(3)
        MCs[MC].SetName(MC)
        MCs[MC].SetTitle(MC)
        if "data" in MC:
            MCs[MC].SetLineColor(datacolors[dataIndex])
            MCs[MC].SetMarkerColor(datacolors[dataIndex])
            MCs[MC].SetMarkerStyle(20+dataIndex)
            upperMG.Add(MCs[MC],"p")
            dataIndex = dataIndex + 1
        else: 
            MCs[MC].SetLineColor(colors[i+1])
            MCs[MC].SetFillColorAlpha(colors[i+1], 0.35)
            upperMG.Add(MCs[MC],"l3")
            MCs[MC].SetLineStyle(i+2)
    dataIndex = 0

    # drawing
    if(ylog): 
        upper.SetLogy()
        upperMG.SetMinimum(gmin/2)
        upperMG.SetMaximum(gmax*2)
    if(xlog): upper.SetLogx()
    upperMG.Draw("AP")

    if("AuAu200" in plotDir):
        xmax = 2

    if("v_{2}" in title): 
        xmin = 0
        xmax = 5
        if("AuAu200" in plotDir):
            xmax = 2
        upline = ROOT.TLine(xmin,0,xmax,0)
        upline.SetLineStyle(9)
        upline.Draw("SAME")
    else: 
        upperMG.SetMinimum(3e-5)

    upperMG.GetXaxis().SetLimits(xmin,xmax)

    pt.Draw("SAME NDC")

    legend = upper.BuildLegend()
    legend.SetY1(0.10)
    legend.SetY2(0.45)
    legend.SetX1(0.17)
    legend.SetX2(0.4)
    legend.SetTextSize(0.04)
    if("v_{2}" in title): 
        legend.GetListOfPrimitives().RemoveLast()
        legend.SetY1(0.5)
        legend.SetY2(0.85)

    # ratio plots
    lowerMG = ROOT.TMultiGraph()

    # Data graph
    defRatio = divideGraph(dataPlot,dataPlot)
    defRatio.SetFillColorAlpha(ROOT.kBlack, 0.15)
    lowerMG.Add(defRatio,"l3")
    #lowerMG.GetYaxis().SetLabelSize(0.06)
    #lowerMG.GetYaxis().SetTitleSize(0.06)
    #lowerMG.GetXaxis().SetLabelSize(0.06)
    lowerMG.GetYaxis().SetTitleSize(0.1)
    lowerMG.GetYaxis().SetTitleOffset(0.6)
    lowerMG.GetYaxis().SetTitle("Set/Default")
    lowerMG.GetYaxis().CenterTitle(True)
    lowerMG.GetXaxis().CenterTitle(True)
    plotvars[0] = plotvars[0].replace('$', '')
    lowerMG.GetXaxis().SetTitle(plotvars[0])

    MCratios = []
    for i, MC in enumerate(MCs):
        thisRatio = divideGraph(MCs[MC],dataPlot)
        MCratios.append(thisRatio.Clone())
        MCratios[i].SetLineWidth(3)
        
        if "data" in MC:
            MCratios[i].SetLineColor(datacolors[dataIndex])
            MCratios[i].SetMarkerColor(datacolors[dataIndex])
            MCratios[i].SetMarkerStyle(20+dataIndex)
            lowerMG.Add(MCratios[i],"p")
            dataIndex = dataIndex + 1
        else: 
            MCratios[i].SetLineColor(colors[i+1])
            MCratios[i].SetFillColorAlpha(colors[i+1], 0.35)
            lowerMG.Add(MCratios[i],"l3")
            MCratios[i].SetLineStyle(i+2)

    # drawing ratio plot
    c.cd()
    lower = ROOT.TPad("plot","plot",0,0,1,0.4)
    lower.SetTopMargin(0)
    lower.SetBottomMargin(0.2)
    lower.Draw()
    lower.cd()
    if(xlog): lower.SetLogx()
    lowerMG.SetMaximum(2.1)
    lowerMG.SetMinimum(0.)
    lowerMG.Draw("AP")
    lowerMG.GetXaxis().SetLimits(xmin,xmax)

    line = ROOT.TLine(xmin,1,xmax,1)
    line.SetLineStyle(9)
    line.Draw("SAME")

    filename  = plotDir + title + ".pdf"
    c.Print(filename)
    c.Close()

#############################################################################################
# start of main

# getting data 
dataDir = "/scratch/user/cameron.parker/projects/C-SCAPE-COMP-gamma/external_packages/gammaLoss/data/"

STARfile = ROOT.TFile(dataDir+"AuAu200IDHads.root","r")
STARv2file = ROOT.TFile(dataDir+"AuAu200hadronv2.root","r")
LHCspecfile = ROOT.TFile(dataDir+"PbPb2760Hadrons.root","r")
LHCv2File = ROOT.TFile(dataDir+"PbPb2760HadronV2.root","r")

# getting input
system = argv[1]

# plot details
pt = ROOT.TPaveText(.68, .12, .88, .32, "NDC")
pt.SetFillColor(ROOT.kWhite)

if system == "AuAu200-0-20":
    files = {"Normal freezeout": "/scratch/user/cameron.parker/projects/runs/new/AuAu200/0-20/hydro-hadrons/analysis",
             "Photon freezeout": "/scratch/user/cameron.parker/projects/runs/new/AuAu200/0-20/default/analysis",}
    
    STARspec = STARfile.Get("Figure 1").Get("Graph1D_y17")
    STARv2 = STARv2file.Get("Figure 3a").Get("Graph1D_y2")

    pt.AddText("Au-Au 200 GeV")
    pt.AddText("Centrality 0-20%")

if system == "AuAu200-20-40":
    files = {"Default": "/scratch/user/cameron.parker/projects/runs/new/AuAu200/20-40/default-norm/analysis",
             "Normal freezeout": "/scratch/user/cameron.parker/projects/runs/new/AuAu200/20-40/hydro-hadrons/analysis",
             "Photon freezeout": "/scratch/user/cameron.parker/projects/runs/new/AuAu200/20-40/default/analysis",}
    
    STARspec = STARfile.Get("Figure 1").Get("Graph1D_y18")
    STARv2 = STARv2file.Get("Figure 3g").Get("Graph1D_y2")

    pt.AddText("Au-Au 200 GeV")
    pt.AddText("Centrality 20-40%")

if system == "PbPb2760-0-20":
    files = {"Default": "/scratch/user/cameron.parker/projects/runs/new/PbPb2760/0-20/default/analysis",
             "Normal freezeout": "/scratch/user/cameron.parker/projects/runs/new/PbPb2760/0-20/default/abs-analysis",
             "Photon freezeout": "/scratch/user/cameron.parker/projects/runs/new/PbPb2760/0-20/thermal-K/analysis",
             "Thermal K factor+Abs": "/scratch/user/cameron.parker/projects/runs/new/PbPb2760/0-20/thermal-K/abs-analysis",
             "#PsiCE": "/scratch/user/cameron.parker/projects/runs/new/PbPb2760/0-20/PCE/analysis",
             "#PsiCE+Abs": "/scratch/user/cameron.parker/projects/runs/new/PbPb2760/0-20/PCE/abs-analysis",
             "#PsiCE-5": "/scratch/user/cameron.parker/projects/runs/new/PbPb2760/0-20/PCE-5/analysis",
             "#PsiCE-5+Abs": "/scratch/user/cameron.parker/projects/runs/new/PbPb2760/0-20/PCE-5/abs-analysis",}
        
    LHCspec = LHCspecfile.Get("Table 1").Get("Graph1D_y1")
    LHCv2 = LHCv2File.Get("Table 7").Get("Graph1D_y1")

    pt.AddText("Pb-Pb 2.76 TeV")
    pt.AddText("Centrality 20-40%")

if system == "PbPb2760-20-40":
    files = {"Default": "/scratch/user/cameron.parker/projects/runs/new/PbPb2760/20-40/default-norm/analysis",
             "Changed Normalization": "/scratch/user/cameron.parker/projects/runs/new/PbPb2760/20-40/hydro-hadrons/analysis",
             "Photon settings": "/scratch/user/cameron.parker/projects/runs/new/PbPb2760/20-40/default/analysis",}
    
    LHCspec = LHCspecfile.Get("Table 12").Get("Graph1D_y1")
    LHCv2 = LHCv2File.Get("Table 2").Get("Graph1D_y3")

    pt.AddText("Pb-Pb 2.76 TeV")
    pt.AddText("Centrality 20-40%")

roots = {option: ROOT.TFile(files[option]+"/analysis.root", "r") for option in files}
Hadrons = {option: ROOT.TGraphErrors(roots[option].Get("hadrons")) for option in roots}
Hadronv2s = {option: roots[option].Get("Hadron v2") for option in roots}

if "PbPb2760" in system:
    Hadrons.update({"ALICE data": LHCspec})
    Hadronv2s.update({"ALICE data": LHCv2})
    pt.AddText("-0.8 < #eta < 0.8")

v2pt = pt.Clone()

if "AuAu200" in system:
    Hadrons.update({"STAR data": STARspec})
    Hadronv2s.update({"PHENIX data": STARv2})
    v2pt.AddText("-0.35 < #eta < 0.35")
    pt.AddText("-0.5 < #eta < 0.5")


# making directory to put plots
plotDir = "plots/" + system + "/hadrons/"
try:
    os.mkdir(plotDir)
except:
    print("Output dir exists, will be overwritten")

# Hadron benchmark
if "AuAu200" in system:
    specVars = ["p_{T} [GeV]", "$#frac{1}{2 #pi N_{ev} p_{T}} #frac{d^{2} N}{dp_{T} d#eta} [GeV^{-2}]", "log", "log", "#pi^{+} p_{T}"]
    v2Vars = ["p_{T} [GeV]", "v_{2}", "lin", "lin", "h^{+} v_{2}"]
else:
    specVars = ["p_{T} [GeV]", "$#frac{1}{2 #pi N_{ev} p_{T}} #frac{d^{2} N}{dp_{T} d#eta} [GeV^{-2}]", "log", "log", "Charged Hadron p_{T}"]
    v2Vars = ["p_{T} [GeV]", "v_{2}", "lin", "lin", "Charged Hadron v_{2}"]

makePlot(Hadrons, specVars, plotDir, pt)
makePlot(Hadronv2s, v2Vars, plotDir, v2pt)