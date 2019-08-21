#!/bin/bash 
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -p general            
#SBATCH --reservation=bioinformatics
#SBATCH -A informatics_workshop
#SBATCH -e rs_pt1_%A.err        # File to which STDERR will be written
#SBATCH -o rs_pt1_%A.out        # File to which STDOUT will be written
#SBATCH -J rs_pt1               # Job name
#SBATCH --mem=1000                  # Memory requested
#SBATCH --time=00:10:00              # Runtime in HH:MM:SS

module purge
module load gcc/4.9.3-fasrc01
module load samtools/1.5-fasrc02
module load bowtie2/2.3.2-fasrc02

# $1 == path to trinity bowtie2 index
# $2 == R1 fastq
# $3 ==R2 fastq
bowtie2 -p 8 -q --no-unal -k 20 -x $1  -1 $2 -2 $3 2>align_stats.txt| samtools view -@10 -Sb -o bowtie2.bam

cat 2>&1 align_stats.txt
