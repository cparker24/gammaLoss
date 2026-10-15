#makes set of brick systems for comparing photon absoroption
import os
import multiprocessing as mp
from datetime import date
from functions import *
import pandas as pd
from operator import itemgetter
import xml.etree.ElementTree as ET

# option reading
softOnly = True
design = pd.read_csv(sys.argv[1])
startdir = ""

# date and file initiation
userDir = "/scratch/user/cameron.parker/projects/"
todays_date = date.today()
totaldir = userDir + "runs/brick/Brick-" + str(todays_date.month) + "-" + str(todays_date.day) + "/"

# xml shenanigans
xmltemplate = userDir + "C-SCAPE-COMP-gamma/external_packages/gammaLoss/config/photon/jetscape_user_gamma_abs_test.xml"

# function to make xmls
def makexml(index, parameters, baseDir):
    # appending option and filename
    filename = baseDir + "/dat/Brick_Bin" + str(index)

    # xml declarations
    newxml = ET.parse(xmltemplate)
    root = newxml.getroot()

    # building lines for xml file
    for temp in root.iter('outputFilename'):
        temp.text = str(filename)
    
    for key in parameters.columns:
        for temp in root.iter(key):
            temp.text = str(parameters.iloc[-1][key])

    xmlname = baseDir+"/xml/brick_"+str(index)+".xml"
    newxml.write(xmlname)

    output = {'xml': xmlname, 'dat': filename}
    return output

# making directory for the runs
def makeDir(index):
    # setting names
    baseDir = totaldir + "points/" + str(index)
    plotsDir = baseDir + "/plots"
    rootDir = baseDir + "/root"
    datDir = baseDir + "/dat"
    xmlDir = baseDir + "/xml"

    #making dirs
    try:
        os.makedirs(plotsDir)
    except:
        pass
    
    try:
        os.makedirs(datDir)
    except:
        pass

    try:
        os.makedirs(rootDir)
    except:
        pass

    try:
        os.makedirs(xmlDir)
    except:
        pass

    #returning base directory name
    return baseDir

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Start of main method

# Creating total directory for use
makeTotalDir(totaldir)

# number of processes to run
indices = list(range(1))

# Changing to build directory
os.chdir(userDir + "C-SCAPE-COMP-gamma/build")

pool = mp.Pool(len(indices))

# Loop over design points and make parameter file
design.to_csv(totaldir+'analysis/parameters.txt',index=False)
print(design)

for i in range(len(design)):    
    # Making xmls
    baseDir = makeDir(i)
    output = pool.starmap(makexml, [(index, design.loc[[i]], baseDir) for index in indices]) # star map to each set of bounds

    # Running jetscape for them
    xmls, dats = zip(*map(itemgetter('xml', 'dat'), output)) 
    writeXmls(xmls,baseDir)
    #continue
    pool.map(runxml, xmls)

pool.close()