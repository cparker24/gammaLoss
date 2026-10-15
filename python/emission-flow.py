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

noflowvec = ROOT.TVector3(0,0,0)
flowvec = ROOT.TVector3(0.3,0.3,-0.3)

def GetHist(datfilepath, medvec):
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
            momVec = ROOT.TLorentzVector(part.Px,part.Py,part.Pz,part.E)
            momVec.Boost(-medvec)
            totalhist.Fill(momVec.E()/temp, momVec.E()/temp) #temp scaling
            totalphotons += 1

    print("Total photons:", totalphotons)

    totalhist.Scale(1.0/events, "width")
    totalhist.Scale(1./(alpha*alphaS*pow(temp,4)))
    totalhist.Scale(1./(volume*resample))
    totalhist.SetMarkerSize(2)
    totalhist.SetMarkerStyle(21)
    totalhist.SetStats(0)

    return totalhist

# initializing histograms
analysisDir = "/scratch/user/cameron.parker/projects/runs/brick/Brick-emission/"
totalhist1 = GetHist("/scratch/user/cameron.parker/projects/runs/brick/Brick-emission/nlo_test_out_final_state_hadrons.dat", noflowvec)
totalhist2 = GetHist("/scratch/user/cameron.parker/projects/runs/brick/Brick-emission/boost_test_out_final_state_hadrons.dat", flowvec)

c1 = ROOT.TCanvas("c1","c1",1400,1200)
ROOT.gStyle.SetOptTitle(0)
                   
totalhist1.SetTitle("No Flow")
totalhist1.GetYaxis().SetTitleOffset(0.9)
totalhist1.Draw("EX0 P")

t = ROOT.TPaveText(0.2, 0.9, 0.8, 1.0,'NDC')
t.AddText("Medium Flow Emission Test")
t.SetFillColor(ROOT.kWhite)
t.Draw("SAME")

totalhist2.SetMarkerColor(ROOT.kBlue)
totalhist2.SetTitle("Flow")
totalhist2.Draw("EX0 P SAME")

c1.BuildLegend()
c1.SetLogy()
c1.Print(analysisDir+"flow.png")
c1.Close()