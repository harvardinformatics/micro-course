#!/bin/bash 
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p serial_requeue,shared
#SBATCH -A informatics_workshop
#SBATCH -e rcorrector_%A.e
#SBATCH -o rcorrector_%A.o
#SBATCH -J rcorrector
#SBATCH --mem=24000
#SBATCH --time=00:20:00

module purge
module load Rcorrector/20180919-fasrc01

# $1 = comma-separated list of R1 files
# $2 = comma-separated list of R2 files
# one could also specify an optional -od $(pwd) cmd to and supply a 3rd cmd line argument to specify 
#explicitly the output directory. NOTE: this wouldn't have to be the current working directory

perl /n/helmod/apps/centos7/Core/Rcorrector/20180919-fasrc01/bin/run_rcorrector.pl -t 8  -1 $1  -2 $2 

