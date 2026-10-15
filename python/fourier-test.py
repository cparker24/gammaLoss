from array import array
from ctypes import Array
from math import pi
from reader import FinalHadronReader
from reader import Particle
from functions import *
import ROOT
import os
import sys
import pythia8

pTs = [0.9,1.1]
phis = [i*pi/10 for i in range(20+1)]
hist = ROOT.TH2D('Photon Rate','Photon Counts; phi; pT (GeV); Photons Out',
                      len(phis)-1,phis[0],phis[len(phis)-1],
                      len(pTs)-1,array('d',pTs))

for i in range(hist.GetNbinsX()):
    hist.SetBinContent(i+1, 1, 1 + cos(2 * hist.GetXaxis().GetBinCenter(i+1)))


canvas = ROOT.TCanvas()
hist.Draw("COLZ")
canvas.Print("fourier-test.png")
print(Getv2(hist,1))