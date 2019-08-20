#!/bin/bash 
#SBATCH -p shared,serial_requeue       # Partition to submit to 
#SBATCH -n 1                   # Number of cores 
$SBATCH -A informatics_workshop
#SBATCH -t 00:20:00               # Runtime in days-hours:minutes 
#SBATCH --mem 1500              # Memory in MB 
#SBATCH -J trinstats               # job name 
#SBATCH -o trinstats.%A.out        # File to which standard out will be written 
#SBATCH -e trinstats.%A.err        # File to which standard err will be written 

module purge
module load trinityrnaseq/2.8.5-fasrc01

/n/helmod/apps/centos7/Core/trinityrnaseq/2.8.5-fasrc01/util/TrinityStats.pl $1 > ${1}.assemblymetrics.txt


