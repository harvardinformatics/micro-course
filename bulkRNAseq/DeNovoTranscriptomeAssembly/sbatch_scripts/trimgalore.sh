#!/bin/sh
#SBATCH -n 1 
#SBATCH -t 05:00:00 #Runtime in minutes
#SBATCH -p general
#SBATCH --reservation=bioinformatics
#SBATCH -A informatics_workshop
#SBATCH -e tgalore_ERR1101637_%A.err
#SBATCH -o tgalore_ERR1101637_%A.out
#SBATCH --mem=3000 
#SBATCH -J tgalore
module purge
module load cutadapt/1.8.1-fasrc01

# $1 == R1
# $2 == R2
# stringency == number of bases at end of read that must match adapter sequence, 1 is extremely stringent
# -e is the error tolerance for the match to the adapter sequence
# even though we don't use the unpaired reads, i.e. if one is completely removed bec <36 bp post trimming,
# we opt to keep those reads in the output, in case there is some downstream application you'd want to use them in

/n/scratchlfs/informatics/nanocourse/rna-seq/denovo_assembly/util/TrimGalore/trim_galore --paired  --illumina --retain_unpaired  --phred33  --output_dir $(pwd) --length 36 -q 5 --stringency 1 -e 0.1 $1 $2
