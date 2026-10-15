# makes set of brick systems for comparing photon absoroption
from multiprocessing import connection
import os
from datetime import date
from functions import *
import pandas as pd
from operator import itemgetter
import xml.etree.ElementTree as ET
import ROOT
from functions import parton

# setting directory for analysis
analysisDir = sys.argv[1]
directories = getDirs(analysisDir)
outdir = analysisDir+"analysis/"
outfile = ROOT.TFile.Open(outdir+'analysis.root', "RECREATE")

# reading design
design = pd.read_csv(outdir+"parameters.txt")
rawlist = design[design.columns[0]].to_list()
step = (rawlist[1]-rawlist[0])/2.0
paramarray = [point-step for point in rawlist]
paramarray.append(rawlist[-1]+step)
paramname = design.columns[0]

# initializing histograms
totalhist = ROOT.TH1D('Photon Rate','Photon Rate; param; Photons Out (%)',
                      len(paramarray)-1,np.array(paramarray))
eventHist = ROOT.TH1D('Event Count','Event Count; param; pT (GeV); Events',
                      len(paramarray)-1,np.array(paramarray))

# directory loop
for i,dir in enumerate(directories):
    xml = ET.parse(analysisDir+"points/"+dir+"/xml/brick_0.xml")
    root = xml.getroot()

    # reading length
    length = 0
    for temp in root.iter('brick_length'):
        length = float(temp.text)

    # reading pT
    thispT = 0
    for temp in root.iter('pT'):
        thispT = float(temp.text)

    # reading partons
    print(f'Reading {dir}')
    datfile = open(analysisDir+"points/"+dir+"/dat/Brick_Bin0_final_state_hadrons.dat","r")
    datlines = datfile.readlines()
    for line in datlines:
        if 'Event' in line:
            eventHist.Fill(rawlist[i])

        if '#' in line:
            continue

        part = parton(line)
        if part.PID == 22 and part.index == 0 and part.status < 0:
            totalhist.Fill(rawlist[i])

# fits
totalhist.Scale(100.0)
totalhist.Divide(eventHist)

totalhist = InvertHist(totalhist)

totalhist.GetYaxis().SetTitle("Photons Out (%)") 
totalhist.GetXaxis().SetTitle("v_{y}") 
totalhist.SetStats(0)
totalhist.SetLineWidth(3)
totalhist.SetMarkerStyle(21)

# output
c1 = ROOT.TCanvas("c1","c1",1400,1200)
totalhist.Draw("EX0 P")
c1.Print(outdir+"plot.png")
c1.Close()
totalhist.Write()

# end behavior
totalhist.Write()
outfile.Close()