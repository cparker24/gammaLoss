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
softOnly = True
startdir = ""

# date and file initiation
userDir = "/scratch/user/cameron.parker/projects/"
todays_date = date.today()
totaldir = userDir + "runs/Hydro-" + str(todays_date.month) + "-" + str(todays_date.day) + "/"

# xml shenanigans
xmltemplate = userDir + "C-SCAPE-COMP-gamma/test/sims_tune.xml"

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
indices = list(range(5))

# Changing to build directory
os.chdir(userDir + "C-SCAPE-COMP-gamma/build")

for i in range(len(indices)):    

    # Making xmls
    baseDir = makeDir(i)
    os.chdir(baseDir)
    print("Making hydro event for",baseDir)
    shutil.copyfile("/scratch/user/cameron.parker/projects/C-SCAPE-COMP-gamma/build/music_input","music_input")

    # Running jetscape for them
    slurm = Slurm(**yaml.load(open('/scratch/user/cameron.parker/projects/C-SCAPE-COMP-gamma/external_packages/gammaLoss/python/short-job.yaml'),Loader=yaml.FullLoader))
    slurm.add_arguments(job_name='hydro'+str(i))
    slurm.add_cmd("source ~/modules.sh")
    slurm.add_cmd("cd " + baseDir)
    #slurm.sbatch("pwd")
    slurm.sbatch("/scratch/user/cameron.parker/projects/C-SCAPE-COMP-gamma/build/runJetscape " + xmltemplate)
