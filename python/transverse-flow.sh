#!/bin/bash

##Run with sbatch command

##NECESSARY JOB SPECIFICATIONS
#SBATCH --job-name=hydro-analysis       #Set the job name to "JobExample1"
#SBATCH --time=20:00:00              #Set the wall clock limit
#SBATCH --ntasks=1                  #Request n tasks
#SBATCH --output=/scratch/user/cameron.parker/slurmoutput/transverse-analysis.txt      #Send stdout/err to "Example1Out.[jobID]"
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G

##OPTIONAL JOB SPECIFICATIONS
#SBATCH --mail-type=ALL              #Send email on all job events
#SBATCH --mail-user=cameron.parker@tamu.edu    #Send all emails to email_address

#First Executable Line

source ~/modules.sh
echo "loaded modules"

cd /scratch/user/cameron.parker/projects/C-SCAPE-COMP-gamma/external_packages/gammaLoss/python/
python3 run-brick.py designs/transverse-flows.txt
python3 1D-analysis.py /scratch/user/cameron.parker/projects/runs/brick/Brick-5-27/