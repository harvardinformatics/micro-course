#!/bin/bash 
#SBATCH -n 12
#SBATCH -p shared,serial_requeue,general
#SBATCH -e rsem_%A.err
#SBATCH -o rsem_%A.out
#SBATCH -J rsem
#SBATCH --mem=36000
#SBATCH -t 06:00:00

module purge
module load rsem/1.2.29-fasrc03

# $1 = R1
# $2 = R2
# $3 = ref genome without fasta/fa file extension
# $4 = outfile prefix
# single end mode requires you have an independent estimate of sequenced fragment size mean and sd

rsem-calculate-expression --bowtie2 -p 12 --time --paired-end $1 $2 $3 $4
# with a reference genome, we can also supply --output-genome-bam so the alignments will be in a format for easy visualization in IGV
