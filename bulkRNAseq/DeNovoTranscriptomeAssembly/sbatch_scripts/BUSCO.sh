#!/bin/bash
#SBATCH -J BUSCO
#SBATCH -N 1                   # Ensure that all cores are on one machine
#SBATCH -n 16                  # Use 16 cores for the job
#SBATCH -t 04:00:00              # Runtime in D-HH:MM
#SBATCH -p shared,serial_requeue      # Partition to submit to
#SBATCH --mem=32000            # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o BUSCO.%A.out  # File to which STDOUT will be written
#SBATCH -e BUSCO.%A.err  # File to which STDERR will be written

module purge
module load BUSCO/3.1.0-fasrc01
module load R/3.5.1-fasrc01

export BUSCO_CONFIG_FILE="/n/scratchlfs/informatics/nanocourse/rna-seq/denovo_assembly/util/BUSCO/config.ini"

# $1 input fasta file (your assembly, e.g. Trinty.fasta)
# $2 output directory (BUSCO will prepend run_)
# $3 lineage directory see lineages to choose from at: http://busco.ezlab.org/  

run_BUSCO.py  -c 16 -o $2 -i $1 -l $3 -m transcriptome 


#Example submission:   sbatch BUSCO.sh Trinity.fasta trinity_BUSCO /n/holylfs/EXTERNAL_REPOS/INFORMATICS/BUSCO/eukaryota_odb9/
