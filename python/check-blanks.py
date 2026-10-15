from functions import *
import os

# setting directory for analysis
analysisDir = sys.argv[1]
directories = getDirs(analysisDir)

for dir in directories:
    datFile = analysisDir+"/points/"+dir+"/test_out_final_state_hadrons.dat"
    filesize = os.path.getsize(datFile)
    if(filesize == 0):
        print(dir)
