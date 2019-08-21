#!/bin/bash 
#SBATCH -p general       # Partition to submit to 
#SBATCH -n 1                   # Number of cores 
#SBATCH --reservation=bioinformatics
#SBATCH -A informatics_workshop
#SBATCH -t 01:00:00               # Runtime in days-hours:minutes 
#SBATCH --mem 1500              # Memory in MB 
#SBATCH -J fastqc               # job name 
#SBATCH -o fastqc.%A.out        # File to which standard out will be written 
#SBATCH -e fastqc.%A.err        # File to which standard err will be written 

module purge
module load fastqc/0.11.8-fasrc01

fastqc --outdir `pwd` $1
