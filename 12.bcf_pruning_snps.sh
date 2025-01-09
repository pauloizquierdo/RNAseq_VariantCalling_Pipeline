#!/bin/bash
#SBATCH --job-name=prusnps                # Job name
#SBATCH --nodes=1                           # Run all processes on a single node       
#SBATCH --ntasks=1                          # Run a single task        
#SBATCH --cpus-per-task=10                 # Number of CPU cores per task
#SBATCH --mem=50GB                          # Total memory limit
#SBATCH --time=04:00:00                     # Time limit hrs:min:sec
#SBATCH -A glbrc
#SBATCH --output=snps.out
#SBATCH --error=snps.err
#SBATCH --constraint="[intel18|amd20|amd22]"

cd /mnt/research/glbrc_group/lowry_lab/RNASeq_cistrans/population

# Load modules
module purge
module load GCC/12.2.0
module load BCFtools/1.17

bcftools view -v snps population.vcf.gz -Ob -o population_snps.bcf

bcftools view -i 'QUAL>40 && INFO/DP>20 && F_MISSING<0.1' population_snps.bcf -Ov -o population_snps_Q40_DP20_FM01.vcf

scontrol show job > population_pruning_snps.txt