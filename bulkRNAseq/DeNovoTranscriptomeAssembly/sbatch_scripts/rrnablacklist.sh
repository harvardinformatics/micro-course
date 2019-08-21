#!/bin/bash
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -A informatics_workshop
#SBATCH -t 00:20:00  
#SBATCH -p shared,serial_requeue 
#SBATCH --mem=24000  
#SBATCH -e blacklist_ERR1101637_%A.e
#SBATCH -o blacklist_ERR1101637_%A.o
module purge
module load bowtie2/2.3.2-fasrc02

# $1 == R1
# $2 == R2
# each output fq entry generates a pair with embedded "1" and "2" strings to indicate which read
# paired_unaligned are those you use moving forward, because neither read could be mapped to the rRNA database

bowtie2 --quiet --very-sensitive-local --phred33  -x /n/scratchlfs/informatics/nanocourse/rna-seq/denovo_assembly/databases/silva/DNA_SILVA_132_tax_SSUParc_LSUParcPlusCarolinensisReferenceRNA -1 $1  -2 $2 --threads 8 --met-file ERR1101637_bowtie2_metrics.txt --al-conc-gz paired_aligned_ERR1101637.fq.gz --un-conc-gz paired_unaligned_ERR1101637.fq.gz  --al-gz unpaired_aligned_ERR1101637.fq.gz --un-gz unpaired_unaligned_ERR1101637.fq.gz
