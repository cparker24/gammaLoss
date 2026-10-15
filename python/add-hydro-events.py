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

# starting options initiation
userDir = "/scratch/user/cameron.parker/projects/"
totaldir = argv[1]
xmltemplate = argv[1] + "config.xml"
addedEvents = int(argv[2])

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

# number of hydro events to run
startIndex = len(os.listdir(totaldir+"points"))
indices = list(range(startIndex, startIndex+addedEvents))
print(indices)

# Changing to build directory
os.chdir(userDir + "C-SCAPE-COMP-gamma/build")

for i in indices:    
    # Making xmls
    baseDir = makeDir(i)
    os.chdir(baseDir)
    print("Making hydro event for",baseDir)
    shutil.copyfile("/scratch/user/cameron.parker/projects/C-SCAPE-COMP-gamma/build/music_input","music_input")
    thisXml = MakeXml(i,baseDir,totaldir+"/config.xml")

    # Running jetscape for them
    slurm = Slurm(**yaml.load(open('/scratch/user/cameron.parker/projects/C-SCAPE-COMP-gamma/external_packages/gammaLoss/python/default.yaml'),Loader=yaml.FullLoader))
    slurm.add_arguments(job_name='hydro'+str(i))
    slurm.add_arguments(mem="10G")
    slurm.add_cmd("source ~/modules.sh")
    slurm.add_cmd("cd " + baseDir)
    #slurm.sbatch("pwd")
    slurm.sbatch("/scratch/user/cameron.parker/projects/C-SCAPE-COMP-gamma/build/runJetscape " + thisXml)
