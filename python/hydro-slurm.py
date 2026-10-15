from sys import argv
from simple_slurm import Slurm
import yaml
import os
from datetime import date
from functions import *
import pandas as pd
from operator import itemgetter
import xml.etree.ElementTree as ET
import shutil

# option reading
runName = argv[1]
startdir = ""

npoints = 100
if len(argv) == 3:
    npoints = int(argv[2])

print("System:", runName)
print("Points to be made:", npoints)

# date and file initiation
userDir = "/scratch/user/cameron.parker/projects/"
todays_date = date.today()
totaldir = userDir + "runs/" + runName + "/"

# xml shenanigans
xmltemplate = userDir + "C-SCAPE-COMP-gamma/external_packages/gammaLoss/config/photon/"+runName+".xml"

# making directory for the runs
def makeDir(index):
    # setting names
    baseDir = totaldir + "points/" + str(index)

    try:
        os.makedirs(baseDir)
    except:
        pass

    #returning base directory name
    return baseDir

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Start of main method

# Creating total directory for use
makeTotalDir(totaldir)
shutil.copyfile(xmltemplate, totaldir+"/config.xml")

# number of hydro events to run
indices = list(range(npoints))

# Changing to build directory
os.chdir(userDir + "C-SCAPE-COMP-gamma/build")

for i in range(len(indices)):    

    # Making xmls
    baseDir = makeDir(i)
    os.chdir(baseDir)
    print("Making hydro event for",baseDir)
    shutil.copyfile("/scratch/user/cameron.parker/projects/C-SCAPE-COMP-gamma/build/music_input","music_input")
    thisXml = MakeXml(i,baseDir,totaldir+"/config.xml")

    # Running jetscape for them
    slurm = Slurm(**yaml.load(open('/scratch/user/cameron.parker/projects/C-SCAPE-COMP-gamma/external_packages/gammaLoss/python/default.yaml'),Loader=yaml.FullLoader))
    slurm.add_arguments(job_name= runName + "=" + str(i))
    slurm.add_arguments(mem="10G")
    slurm.add_cmd("source ~/modules.sh")
    slurm.add_cmd("cd " + baseDir)
    #slurm.sbatch("pwd")
    slurm.sbatch("/scratch/user/cameron.parker/projects/C-SCAPE-COMP-gamma/build/runJetscape " + thisXml)
