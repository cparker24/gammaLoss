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

analysisDir = sys.argv[1]
directories = getDirs(analysisDir)
pythia = pythia8.Pythia()

#data openings
photonDatav2file = ROOT.TFile("data/PbPb2760v2.root")
photonDatav2 = photonDatav2file.Get("Table 8").Get("Graph1D_y1")
hadDatav2file = ROOT.TFile("data/PbPb2760HadronV2.root")
hadDatav2 = hadDatav2file.Get("Table 2").Get("Graph1D_y2")
hadronDataSpecFile = ROOT.TFile("data/PbPb2760Hadrons.root")
hadDataSpec1 = hadronDataSpecFile.Get("Table 3").Get("Graph1D_y1")
hadDataSpec2 = hadronDataSpecFile.Get("Table 4").Get("Graph1D_y1")

# histogram declarations
pTs = [0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,3.0,3.3,3.7,4.1,4.6,5.4,6.2]
phis = [i*pi/10 for i in range(20+1)]
totalhist = ROOT.TH2D('Photon Rate','Photon Counts; phi; pT (GeV); Photons Out',
                      len(phis)-1,phis[0],phis[len(phis)-1],
                      len(pTs)-1,array('d',pTs))
QGPHist = ROOT.TH1D('QGP Rate','QGP; pT (GeV); Photons Out',
                      len(pTs)-1,array('d',pTs))
HRGHist = ROOT.TH1D('HRG Rate','HRG; pT (GeV); Photons Out',
                      len(pTs)-1,array('d',pTs))
showerHist = ROOT.TH1D('Shower Rate','Shower; pT (GeV); Photons Out',
                      len(pTs)-1,array('d',pTs))
hardHist = ROOT.TH1D('Hard Process Rate','Hard; pT (GeV); Photons Out',
                      len(pTs)-1,array('d',pTs))
otherhist = ROOT.TH2D('Other Distribution','Hadron Counts; phi; pT (GeV); Hadrons Out',
                      len(phis)-1,phis[0],phis[len(phis)-1],
                      len(pTs)-1,array('d',pTs))
etaHist = ROOT.TH1D('Eta Dist','Photon Counts; eta; Photons Out',
                      20,-2.5,2.5)

# hadron histograms
hadFile = ROOT.TFile("data/PbPb2760Spectra.root")
hadSpecData = hadFile.Get("Table 4").Get("Graph1D_y1")
hadSpecHist = hadFile.Get("Table 4").Get("Hist1D_y1").Clone()
hadSpecHist.Reset()
hadSpecHist.SetNameTitle("Hadron Spec","Hadron Spec")
binArray = hadSpecHist.GetXaxis().GetXbins()
hadHistTotal = ROOT.TH2D("Hadron Total Spec", "Hadron Total Spec", 
                      len(phis)-1,phis[0],phis[len(phis)-1],
                      len(pTs)-1,array('d',pTs))

# output directory
outdir = analysisDir+"analysis/"
outfile = ROOT.TFile.Open(outdir+'analysis.root', "RECREATE")

eventCount = 0
etaCut = 0.8
# file loop
for directory in directories:
    XSCAPEfile = FinalHadronReader(analysisDir+"points/"+directory+"/test_out_final_state_hadrons.dat")
    EPangle = -10
    
    # event loop
    while XSCAPEfile.IsFinished() is False:
        eventCount += 1
        #EPangle = XSCAPEfile.EPangle
        particles = XSCAPEfile.GetHadrons()

        # particle loop
        if particles == None: break

        if EPangle < -9:
            EPangle = GetEPAngle(particles)
            print("Event Plane Angle:",EPangle)

        for particle in particles:
            weight = 1
            delphi = particle.phi()-EPangle
            if delphi > 2*pi: delphi = delphi-2*pi
            if delphi < 0: delphi = delphi+2*pi

            if (particle.PID == 22 and particle.status == 24) or (particle.PID == 22 and particle.status == 25):
                etaHist.Fill(particle.p.Eta())
                weight = 1/25
                
            if particle.status == 11:
                weight = 1/10
                otherhist.Fill(delphi, particle.pT, weight)

            if abs(particle.p.Eta()) > etaCut: continue

            if(pythia.particleData.charge(int(particle.PID)) and particle.status == 11):
                hadSpecHist.Fill(particle.pT)
                hadHistTotal.Fill(delphi, particle.pT, weight)

            if particle.PID == 22 and particle.status != 0:
                totalhist.Fill(delphi, particle.pT, weight)
                if particle.status == 22: hardHist.Fill(particle.pT, weight)
                if particle.status == 23: showerHist.Fill(particle.pT, weight)
                if particle.status == 24: QGPHist.Fill(particle.pT, weight)
                if particle.status == 25: HRGHist.Fill(particle.pT, weight)


totalhist.Scale(1/eventCount)
totalhist.Write()
otherhist.Scale(1/eventCount)
otherhist.Write()

QGPHist.SetLineColor(ROOT.kPink)
QGPHist.Scale(1.0/(2.0*pi*2.0*etaCut*eventCount),"width")
ScaleByCenter(QGPHist)
QGPHist.Write()
ExtractTemp(QGPHist)

HRGHist.SetLineColor(ROOT.kViolet)
HRGHist.Scale(1.0/(2.0*pi*2.0*etaCut*eventCount),"width")
ScaleByCenter(HRGHist)
HRGHist.Write()

hardHist.SetLineColor(ROOT.kGreen)
hardHist.Scale(1.0/(2.0*pi*2.0*etaCut*eventCount),"width")
ScaleByCenter(hardHist)
hardHist.Write()

showerHist.SetLineColor(ROOT.kOrange)
showerHist.Scale(1.0/(2.0*pi*2.0*etaCut*eventCount),"width")
ScaleByCenter(showerHist)
showerHist.Write()

etaHist.Write()

# making total spectrum
pThist = totalhist.ProjectionY("pT spectrum")
pThist.Scale(1.0/(2.0*pi*2.0*etaCut),"width")
ScaleByCenter(pThist)
pThist.Write()

# hadron writing
hadSpecHist.Scale(1.0/(2.0*pi*2.0*etaCut*eventCount),"width")
ScaleByCenter(hadSpecHist)
hadSpecHist.Write()
hadHistTotal.Write()

hadSpecCanvas = ROOT.TCanvas()
hadSpecCanvas.SetLogy()
hadSpecCanvas.SetLogx()
hadDataSpec1.Draw()
hadDataSpec2.Draw("SAME")
hadSpecHist.Draw("SAME")
hadSpecCanvas.Print(outdir+"hadronSpec.png")

# v2 calculations
print("Photon v2 calc: ", Getv2(totalhist,1))
print("Photon v2 fit: ", Getv2Fit(totalhist,1))
print("Test hadron v2 calc: ", Getv2(otherhist,2))

flowHads = Makev2Graph(otherhist, "Hadron v2")
flowHads.Write()
hadCanvas = ROOT.TCanvas()
flowHads.Draw()
hadDatav2.Draw("SAME")
hadCanvas.Print(outdir+"hadronV2.png")

flowPhotons = Makev2Graph(totalhist, "Photon v2")
flowPhotons.Write()
photonCanvas = ROOT.TCanvas()
flowPhotons.Draw()
photonDatav2.Draw("SAME")
photonCanvas.Print(outdir+"photonV2.png")

# getting total spectra
specFile = ROOT.TFile("data/PbPb2760Spectra.root")
specData = specFile.Get("Table 2").Get("Graph1D_y1")

specCanvas = ROOT.TCanvas()
specCanvas.SetLogy()
specData.Draw()
pThist.Draw("SAME")
hardHist.Draw("SAME")
showerHist.Draw("SAME")
QGPHist.Draw("SAME")
HRGHist.Draw("SAME")
specCanvas.BuildLegend()
specCanvas.Print(outdir+"photonSpec.png")

outfile.Close()