#!/bin/bash 
#SBATCH -N 1
#SBATCH -n 24
#SBATCH -p bigmem
#SBATCH -e trinity_%A.e
#SBATCH -o trinity_%A.o
#SBATCH --mail-type=END
#SBATCH --mail-type=BEGIN
#SBATCH -J trinity
#SBATCH --mem=250000 # for full assembly
#SBATCH -t 48:00:00

module purge
module load trinityrnaseq/2.8.5-fasrc01

Trinity --seqType fq --max_memory 240G --min_kmer_cov 1 --left paired_unaligned_ERR1101637.fq.1.gz --right paired_unaligned_ERR1101637.fq.2.gz --output mouse_Trinity_2019.02.21 --CPU 24

###########################################################
#add to trinity call to run parallel jobs on compute farm #
##########################################################
## --grid_exec "/n/home_rc/afreedman/software/HpcGridRunner-1.0.2/hpc_cmds_GridRunner.pl --grid_conf $(pwd)/grid.conf -c"  --grid_node_max_memory 5G 
