import matplotlib as plt
import numpy as np
import pandas as pd
import ROOT
from functions import *

temp = 0.2
alpha = 0.0073
alphaS = 0.3

fmtoGeVinv = 0.1973
volume = 5000/pow(fmtoGeVinv,4)
resample = 10000/1.0

df = pd.read_csv("../build/AMYcheck.out", header=None, delimiter=' ')
def RatePlot():
    df.columns = ["$p/T$", "NLO","LO","LL"]

    print(df)

    plot = df.plot(title="Emission Rate Check", x=0, y=[1,2,3], sharex=True)
    plot.set_ylabel(r"$p \frac{dR}{dp} / (\alpha \alpha_{s} T^{4})$")
    pltimage = plot.get_figure()
    pltimage.savefig("plots/others/AMYrate.png")

def GetHist(datfilepath):
    datfile = open(datfilepath,"r")
    datlines = datfile.readlines()
    totalhist = ROOT.TH1D('Photon Rate','Monte Carlo; p/T; p #frac{dR}{dp} / #alpha #alpha_{s} T^{4}',20,0,10)
    events = 0
    totalphotons = 0
    for line in datlines:
        if "Event" in line:
            events = events + 1
            print("Event:", events)

        if '#' in line:
            continue

        part = parton(line)
        if part.PID == 22 and part.status == 24:
            totalhist.Fill(part.E/temp, part.E/temp) #temp scaling
            totalphotons += 1

    print("Total photons:", totalphotons)

    totalhist.Scale(1.0/events, "width")
    totalhist.Scale(1./(alpha*alphaS*pow(temp,4)))
    totalhist.Scale(1./(volume*resample))
    totalhist.SetMarkerSize(2)
    totalhist.SetMarkerStyle(21)
    totalhist.SetStats(0)

    return totalhist

RatePlot()

# initializing histograms
analysisDir = "/scratch/user/cameron.parker/projects/runs/brick/Brick-emission/"
totalhist1 = GetHist("/scratch/user/cameron.parker/projects/runs/brick/Brick-emission/nlo_test_out_final_state_hadrons.dat")
totalhist2 = GetHist("/scratch/user/cameron.parker/projects/runs/brick/Brick-emission/lo_test_out_final_state_hadrons.dat")
totalhist3 = GetHist("/scratch/user/cameron.parker/projects/runs/brick/Brick-emission/ll_test_out_final_state_hadrons.dat")

x1s = df['$p/T$'].to_numpy()
y1s = df['NLO'].to_numpy()
graph1 = ROOT.TGraph(len(x1s),x1s,y1s)

y2s = df['LO'].to_numpy()
graph2 = ROOT.TGraph(len(x1s),x1s,y2s)

y3s = df['LL'].to_numpy()
graph3 = ROOT.TGraph(len(x1s),x1s,y3s)

c1 = ROOT.TCanvas("c1","c1",1400,1200)
ROOT.gStyle.SetOptTitle(0)
                   
totalhist1.SetTitle("NLO Monte Carlo")
totalhist1.GetYaxis().SetTitleOffset(0.9)
totalhist1.Draw("EX0 P")

t = ROOT.TPaveText(0.2, 0.9, 0.8, 1.0,'NDC')
t.AddText("Emission Rate Tests")
t.SetFillColor(ROOT.kWhite)
t.Draw("SAME")

totalhist2.SetMarkerColor(ROOT.kBlue)
totalhist2.SetTitle("LO Monte Carlo")
totalhist2.Draw("EX0 P SAME")

totalhist3.SetMarkerColor(ROOT.kRed)
totalhist3.SetTitle("LL Monte Carlo")
#totalhist3.Draw("EX0 P SAME")

graph1.SetLineWidth(4)
graph1.SetTitle("NLO")
graph1.Draw("L SAME")

graph2.SetLineColor(ROOT.kBlue)
graph2.SetLineWidth(4)
graph2.SetTitle("LO")
graph2.Draw("L SAME")

graph3.SetLineColor(ROOT.kRed)
graph3.SetLineWidth(4)
graph3.SetTitle("LL")
#graph3.Draw("L SAME")

c1.BuildLegend()
c1.SetLogy()
c1.Print(analysisDir+"data.png")
c1.Close()