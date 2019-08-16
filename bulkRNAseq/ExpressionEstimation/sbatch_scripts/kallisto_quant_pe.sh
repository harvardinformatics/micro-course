#!/bin/bash 
#SBATCH -n 6
#SBATCH -p shared,serial_requeue,general
#SBATCH -e quant_%A.err
#SBATCH -o quant_%A.out
#SBATCH -J quant
#SBATCH --mem=16000
#SBATCH -t 06:00:00

# $1 == sample id
module purge
module load kallisto/0.45.1-fasrc01

kallisto quant -i /n/scratchlfs/informatics/nanocourse/rna-seq/expression/ref/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.transcripts.idx -b 100 -t 6 -o ${1}_quant ${1}_1_val_1.fq ${1}_2_val_2.fq

