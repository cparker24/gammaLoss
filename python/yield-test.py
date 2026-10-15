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
from zmq import Event
import matplotlib.pyplot as plt
import numpy

# setting directory for analysis
analysisDir = sys.argv[1]
directories = getDirs(analysisDir)
outdir = analysisDir+"analysis/"

photonCount = 0
eventCount = 0
energySum = 0

# histogram declarations
pTs = [0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,3.0,3.3,3.7,4.1,4.6,5.4,6.2]
QGPHist = ROOT.TH1D('QGP Rate','QGP; p (GeV); dR/pdp (GeV^3)',
                      len(pTs)-1,array('d',pTs))

# directory loop
for dir in directories:
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
            eventCount += 1

        if '#' in line:
            continue

        part = parton(line)
        if part.PID == 22 and part.status == 24:
            photonCount += 1
            energySum += part.E
            QGPHist.Fill(part.E)

    datfile.close()
    del datfile, datlines
    print("Finished")

print("Photons per event:", photonCount/eventCount)

# Calculation
alpha = 1/137

def alphaS(temp):
  return 6.0 * pi / (27.0 * math.log(temp / 0.022))

def gS(temp):
  return sqrt(4.0 * pi * alphaS(temp))

def B(temp):
   return 24.*alpha*alphaS(temp)*pow(temp,4)/(9.*pi)

def fermi(p, T):
    x = p / T
    expo = math.exp(-x)
    return expo / (1.0 + expo)

def f(x):
   return x/(1+math.exp(x))

def g(x):
   C22 = 0.041 / x - 0.3615 + 1.01 * math.exp(-1.35 * x)
   Cab = sqrt(1.5) * (0.548 / pow(x, 1.5) * math.log(12.28 + 1.0 / x) + 0.133 * x / sqrt(1.0 + x / 16.27))
  
   return (0.5 * math.log(2.0 * x) + C22 + Cab) * f(x)

def AMYdRd3p(p, T):
  x = p / T
  expo = math.exp(-x)
  fermi = expo / (1.0 + expo)
  
  prfph = (6./9.) * (alpha*alphaS(T)*pow(T,2))/(pow(pi,2)*p) * fermi
  
  C22 = 0.041 / x - 0.3615 + 1.01 * math.exp(-1.35 * x)
  Cab = sqrt(1.5) * (0.548 / pow(x, 1.5) * math.log(12.28 + 1.0 / x) + 0.133 * x / sqrt(1.0 + x / 16.27))
  
  Ctot = 0.5 * math.log(2.0 * x) + C22 + Cab
  
  dRd3p = prfph * (math.log(sqrt(3.0) / gS(T)) + Ctot)
  return dRd3p

def AMYdRdx(x, T):
   return B(T)*(math.log(sqrt(3)/gS(T))*f(x) + g(x))

def y(x, T):
   return (x*T) * AMYdRdx(x, T) / (alpha*alphaS(T)*pow(T,4))

temp = 0.2
integral1 = 0.384869
integral2 = 0.549816
volumeSampled = 73760
fmToGeVinv = 5.067728853174532

Bfactor = 24.*alpha*alphaS(temp)*pow(temp,4)/(9.*pi)
integralFactor = math.log(sqrt(3)/gS(temp))*integral1 + integral2
photonsProduced = Bfactor*integralFactor*volumeSampled*pow(fmToGeVinv,4)
print("Expected photons per event:", photonsProduced)

xs = np.linspace(0.01,10,100)
ys = [y(x,temp) for x in xs]

# plot the data
plt.figure(figsize=(5,5))
plt.plot(xs, ys, color='tab:blue')
plt.xlabel(r"$p/T$")
plt.ylabel(r"$(p*dR/dp)/(\alpha_{EM} \alpha_s)$")
plt.xlim(0,10)
plt.ylim(0,0.5)
plt.savefig('plots/energy-yield.png', dpi=1000)

# energy comparisons
start = 2
end = 100
count = 1000
width = (end-start)/count
xs = np.linspace(start,end,count)

energyPrediction = 0
for x in xs:
   energyPrediction += x*temp*AMYdRdx(x,temp)*width*volumeSampled*pow(fmToGeVinv,4)

print("Energy from X-SCAPE:", energySum/eventCount)
print("Energy prediction:", energyPrediction)

# plotting spectra
# x[0] = p, p[0] = T
def AMYpdRdp(x: numpy.ndarray, pars: numpy.ndarray) -> float:
   pgamma = x[0]
   temp = pars[0]
   return AMYdRdx(pgamma/temp,temp)/(temp*x[0])

AMYfunc = ROOT.TF1("AMY", AMYpdRdp, 0.5, 7, npar=1, ndim=1)
AMYfunc.SetParameter(0,temp)
QGPHist.SetLineColor(ROOT.kBlue)
QGPHist.Scale(1.0/(volumeSampled*pow(fmToGeVinv,4)*eventCount),"width")
ScaleByCenter(QGPHist)
ExtractTemp(QGPHist)
specCanvas = ROOT.TCanvas()
specCanvas.SetLogy()
QGPHist.Draw()
AMYfunc.Draw("SAME")
specCanvas.BuildLegend()
specCanvas.Print("plots/spec-comparison.png")