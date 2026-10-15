from math import cos, sin, sqrt
import os
import sys
import numpy as np
import random
import pandas as pd
import math
from reader import *
from array import array
import xml.etree.ElementTree as ET
import xml.dom.minidom as minidom

# make base directories for the analysis
def makeTotalDir(totaldir):
    try:
        os.makedirs(totaldir)
    except:
        pass

    try:
        os.makedirs(totaldir + "points/")
    except:
        pass

    try:
        os.makedirs(totaldir + "analysis/")
    except:
        pass

    try:
        os.makedirs(totaldir + "abs-analysis/")
    except:
        pass

# gets list of directories in the target directory
def getDirs(dir):
    startdir = os.getcwd()
    os.chdir(dir+"points")
    files = os.listdir(".")

    directories = []
    for item in files:
        if os.path.isdir(item) and item != 'analysis': directories.append(item)

    os.chdir(startdir)

    
    directories.sort()  # to make sure the directories are sorted in case it doesnt complete
    return directories

# to see what ran last if runs dont complete
def update(cmd):
    outfile = open('last-ran.txt','a')
    outfile.write(cmd+'\n')
    outfile.close()

    print(cmd)

# runs the xml in jetscape
def runxml(xmlname):
    runjetscapecmd = "./runJetscape " + xmlname
    update(runjetscapecmd)
    os.system(runjetscapecmd)

# writing set of xmls to a file to be run later
def writeXmls(xmls,baseDir):
    outfile = open(baseDir+"/xml/cmds.txt","w")
    outfile.writelines(xml+"\n" for xml in xmls)
    outfile.close()

# reading xmls to run
def readXmls(baseDir):
    infile = open(baseDir+"/xml/cmds.txt","r")
    lines = infile.readlines()
    return lines

# sorting list of strings by int value
def intSort(dirs):
    dirs = [int(x) for x in dirs]
    dirs.sort()
    dirs = [str(x) for x in dirs]
    return dirs

# gets list of directories in the target directory
def concatDats(dir):
    startdir = os.getcwd()
    os.chdir(dir+"/dat")
    dats = os.listdir(".")

    cmd = "cat "
    for dat in dats:
        cmd = cmd + dat + " "
    cmd = cmd + " > " + dats[0].rstrip(".dat.gz") + "0.dat.gz"
    os.system(cmd)
    
    for dat in dats:
        os.remove(dat)

    os.chdir(startdir)

class parton:
    def __init__(self, index, PID, status, E, Px, Py, Pz):
        self.index = index
        self.PID = PID
        self.status = status
        self.E = E
        self.Px = Px
        self.Py = Py
        self.Pz = Pz
        self.pT = sqrt(Px**2 + Py**2)

    def __init__(self, line):
        params = line.split(' ')
        if len(params) < 7:
            self.index = -1
            self.PID = -1
            self.status = -1
            self.E = -1
            self.Px = -1
            self.Py = -1
            self.Pz = -1
            self.pT = -1
            return

        self.index = int(params[0])
        self.PID = int(params[1])
        self.status = int(params[2])
        self.E = float(params[3])
        self.Px = float(params[4])
        self.Py = float(params[5])
        self.Pz = float(params[6])
        self.pT = sqrt(self.Px**2 + self.Py**2)

    def __str__(self):
        return f"{self.PID} {self.E} {self.pT}"
    
##chance for photon to be absorbed from https:##arxiv.org/abs/hep-ph/9405309
def absFactor1(pVec, T):
    alpha = 1./137.
    alphaS = 0.3
    p = pVec.P()
    return 2.0*(5.*math.pi/9.)*(alpha*alphaS*T*T/p)*math.log(0.2317*p/(alphaS*T))

def absFactor2(pVec, T):
    ## Constants
    alpha = 1.0 / 137.0 
    hbc = 0.1973 
    
    ## Calculate alpha_s at temperature 'temp'
    alphsT = 6.0 * math.pi / (27.0 * math.log(T / 0.022)) 
    gsT = math.sqrt(alphsT * 4.0 * math.pi) 
    
    ## Calculate x and exponential term
    p = pVec.P() 
    x = pVec.P() / T 
    expo = math.exp(-x) 
    fermi = expo / (1.0 + expo) 
    
    ## Calculate prfph
    prfph = alpha * alphsT * T**2 * (5.0 / 9.0) / (p * math.pi**2) 
    
    ## Calculate C22 and Cab
    C22 = 0.041 / x - 0.3615 + 1.01 * math.exp(-1.35 * x) 
    Cab = math.sqrt(1.5) * (0.548 / pow(x, 1.5) * math.log(12.28 + 1.0 / x) + 0.133 * x / math.sqrt(1.0 + x / 16.27)) 
    
    ## Calculate Ctot
    Ctot = 0.5 * math.log(2.0 * x) + C22 + Cab 
    
    ## Calculate dRd3p
    dRd3p = prfph * fermi * (math.log(math.sqrt(3.0) / gsT) + Ctot)
    return dRd3p * 4 * pow(math.pi,3) / expo

def GetEPAngle(hadronList):
    Qx = 0
    Qy = 0

    for hadron in hadronList:
        if hadron.status == 11 and hadron.pT > 1:
            Qx += hadron.pT * cos(2.0*hadron.phi())
            Qy += hadron.pT * sin(2.0*hadron.phi())

    return np.arctan2(Qy,Qx)/2.0

def Getv2(hist, bin):
    integral = 0
    denom = 0

    for i in range(hist.GetNbinsX()):
        integral += hist.GetBinContent(i+1,bin) * cos(2 * hist.GetXaxis().GetBinCenter(i+1)) * hist.GetXaxis().GetBinWidth(i+1)
        denom += hist.GetBinContent(i+1,bin) * hist.GetXaxis().GetBinWidth(i+1)

    if denom <= 0:
        return 0

    return integral/(denom)

def Getv2Fit(hist, bin):
    phiHist = hist.ProjectionX("low pT distribution",bin)

    # getting amplitude coefficient
    integral = 0
    for i in range(phiHist.GetNbinsX()):
        integral += phiHist.GetBinContent(i+1,bin) * phiHist.GetBinWidth(i+1)

    if integral <= 0:
        return 0

    thisfit = ROOT.TF1("thisfit", "([0])*(1 + ([1]*2*cos(2*x)))", 0, 2*pi)
    thisfit.SetParameter(0,integral/(2*pi))
    thisfit.SetParameter(1,0)
    thisfit.SetParNames("A","v2")
    phiHist.Fit("thisfit")

    return thisfit.GetParameter(1)

def Makev2Graph(hist,name):
    pTs = []
    v2s = []

    for i in range(hist.GetNbinsY()):
        x = i+1
        v2 = Getv2(hist,x)
        pT = hist.GetYaxis().GetBinCenter(x)
        pTs.append(pT)
        v2s.append(v2)

    graph = ROOT.TGraph(hist.GetNbinsY(),array('d',pTs), array('d',v2s))
    graph.SetNameTitle(name,name+";p_{T} (GeV); v_{2}")    
    graph.SetLineColor(ROOT.kBlue)
    return graph

def ScaleByCenter(hist):

    for i in range(hist.GetNbinsX()):
        x = i+1
        hist.SetBinContent(x,hist.GetBinContent(x)/hist.GetBinCenter(x))
        hist.SetBinError(x,hist.GetBinError(x)/hist.GetBinCenter(x))

def MultByCenter(hist):

    for i in range(hist.GetNbinsX()):
        x = i+1
        hist.SetBinContent(x,hist.GetBinContent(x)*hist.GetBinCenter(x))
        hist.SetBinError(x,hist.GetBinError(x)*hist.GetBinCenter(x))
    return hist
    

def ExtractTemp(hist):
    QGPfit = ROOT.TF1("thisfit", "[0]*exp(-x/[1])", hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmin())
    QGPfit.SetParNames("A","T")
    QGPfit.SetParameter(1,0.2)
    QGPfit.SetLineColor(ROOT.kBlack)
    hist.Fit("thisfit")
    print("Fit temperature is",int(QGPfit.GetParameter(1)*1000),"MeV")

def MakeXml(index, baseDir, xmlPath):
    # Parse
    tree = ET.parse(xmlPath)
    root = tree.getroot()

    # Find or create Random and seed elements
    random_el = root.find("Random")
    if random_el is None:
        random_el = ET.SubElement(root, "Random")
    seed_el = random_el.find("seed")
    if seed_el is None:
        seed_el = ET.SubElement(random_el, "seed")

    # Set the new seed
    seed_el.text = str(int(index+1000))

    # Pretty-print via minidom (ElementTree writes a single line by default)
    rough_xml = ET.tostring(root, encoding="utf-8")
    pretty_xml = minidom.parseString(rough_xml).toprettyxml(indent="  ", encoding="utf-8")

    # Decide output path
    out_path = baseDir+"/config.xml"

    # Write with XML declaration
    with open(out_path, "wb") as f:
        f.write(rough_xml)

    return str(out_path)

# flips hist that has percent of absorbed photons to percent of passed photons 
def InvertHist(hist):
    for i in range(hist.GetNbinsX()):
        x = i+1
        hist.SetBinContent(x,100.0-hist.GetBinContent(x))

    return hist