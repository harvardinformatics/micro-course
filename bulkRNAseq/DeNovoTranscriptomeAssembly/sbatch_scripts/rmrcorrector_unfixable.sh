#!/bin/bash 
#SBATCH -p shared,serial_requeue       # Partition to submit to 
#SBATCH -n 1                   # Number of cores 
#SBATCH -t 01:00:00               # Runtime in days-hours:minutes 
#SBATCH --mem 1500              # Memory in MB 
#SBATCH -A informatics_workshop
#SBATCH -J rmunfix               # job name 
#SBATCH -o rmunfix.%A.out        # File to which standard out will be written 
#SBATCH -e rmunfix.%A.err        # File to which standard err will be written 

module purge
module load python/2.7.14-fasrc01
# $1 == rcorrector corrected R1 fastq
# $2 == rcorrector corrected R2 fastq
# $3 == sample name (to be included in metrics output log filename)

python /n/scratchlfs/informatics/nanocourse/rna-seq/denovo_assembly/python_scripts/FilterUncorrectabledPEfastq.py -1 $1 -2 $2 -s $3
