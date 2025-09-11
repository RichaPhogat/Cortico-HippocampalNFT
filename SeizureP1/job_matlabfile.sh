#!/bin/bash
 
#PBS -l ncpus=2
#PBS -l mem=190GB
#PBS -l jobfs=10GB
#PBS -q normal
#PBS -P hl36
#PBS -l walltime=47:50:00
#PBS -l storage=gdata/hl36+scratch/hl36
#PBS -l wd
#PBS -l ngpus=0
#PBS -N JobName
#PBS -o /home/575/rr8777/HippSeizureSlowP2I33/outputfile.out
#PBS -e /home/575/rr8777/HippSeizureSlowP2I33/errorfile.err

module load matlab_licence/newcastle
module load matlab/R2023b
nohup matlab -nodisplay -nosplash -nodesktop -r "NFTTwoSheetChecking; exit;"
exit
