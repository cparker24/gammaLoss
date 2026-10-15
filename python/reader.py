from math import sqrt
from math import pi
from matplotlib.pylab import arctan
import ROOT

class Particle:
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
        self.p = ROOT.TLorentzVector(self.Px,self.Py,self.Pz,self.E)

    def __str__(self):
        return f"{self.PID} {self.E} {self.pT}"
    
    def phi(self):
        phi = self.p.Phi()
        if phi < 0:
            phi = 2*pi + phi
        
        return phi
    
class FinalHadronReader:
    def __init__(self, filename) -> None:
        self.file = open(filename,"r")
        self.lines = self.file.readlines()
        self.index = 2
        self.finished = False
        self.eventNumber = 0
        self.EPangle = float(self.lines[1].split("\t")[6])

    # Get the next set of hadrons from the event
    def GetHadrons(self):
        self.eventNumber += 1
        print("Event", self.eventNumber)
        particleList = []

        while self.index < len(self.lines):
            line = self.lines[self.index]
            self.index += 1

            if self.index >= len(self.lines):
                self.finished = True

            # ending at start of the next event
            if "Event" in line: 
                self.EPangle = float(line.split("\t")[6])
                return particleList
            
            
            # ending at end of the file
            if "sigmaGen" in line: 
                self.finished = True
                args = line.split('\t')
                self.sigmaGen = float(args[2])
                self.sigmaErr = float(args[4])
                self.weight = float(args[6])
                return particleList
            
            # skipping comments
            if "#" in line: 
                continue

            particleList.append(Particle(line))

    # making sure there is more to read
    def IsFinished(self):
        return self.finished