#!/bin/bash
#SBATCH --job-name=GVCF                # Job name
#SBATCH --nodes=1                           # Run all processes on a single node       
#SBATCH --ntasks=1                          # Run a single task        
#SBATCH --cpus-per-task=10                 # Number of CPU cores per task
#SBATCH --mem=50GB                          # Total memory limit
#SBATCH --time=04:00:00                     # Time limit hrs:min:sec
#SBATCH -A glbrc
#SBATCH --output=concatenate.out
#SBATCH --error=concatenate.err
#SBATCH --constraint="[intel18|amd20|amd22]"

cd /mnt/research/glbrc_group/lowry_lab/RNASeq_cistrans/population

# Load modules
module purge
module load GCC/7.3.0-2.30  OpenMPI/3.1.1
module load GATK/4.1.4.1-Python-3.6.6

gatk GatherVcfs \
  -I Chr01K_pop.vcf.gz \
  -I Chr01N_pop.vcf.gz \
  -I Chr02K_pop.vcf.gz \
  -I Chr02N_pop.vcf.gz \
  -I Chr03K_pop.vcf.gz \
  -I Chr03N_pop.vcf.gz \
  -I Chr04K_pop.vcf.gz \
  -I Chr04N_pop.vcf.gz \
  -I Chr05K_pop.vcf.gz \
  -I Chr05N_pop.vcf.gz \
  -I Chr06K_pop.vcf.gz \
  -I Chr06N_pop.vcf.gz \
  -I Chr07K_pop.vcf.gz \
  -I Chr07N_pop.vcf.gz \
  -I Chr08K_pop.vcf.gz \
  -I Chr08N_pop.vcf.gz \
  -I Chr09K_pop.vcf.gz \
  -I Chr09N_pop.vcf.gz \
  -I scaffold_174_pop.vcf.gz \
  -I scaffold_19_pop.vcf.gz \
  -I scaffold_196_pop.vcf.gz \
  -I scaffold_218_pop.vcf.gz \
  -I scaffold_260_pop.vcf.gz \
  -I scaffold_283_pop.vcf.gz \
  -I scaffold_319_pop.vcf.gz \
  -I scaffold_320_pop.vcf.gz \
  -I scaffold_337_pop.vcf.gz \
  -I scaffold_338_pop.vcf.gz \
  -I scaffold_34_pop.vcf.gz \
  -I scaffold_381_pop.vcf.gz \
  -I scaffold_390_pop.vcf.gz \
  -I scaffold_401_pop.vcf.gz \
  -I scaffold_413_pop.vcf.gz \
  -I scaffold_414_pop.vcf.gz \
  -I scaffold_428_pop.vcf.gz \
  -I scaffold_434_pop.vcf.gz \
  -I scaffold_448_pop.vcf.gz \
  -I scaffold_451_pop.vcf.gz \
  -I scaffold_467_pop.vcf.gz \
  -I scaffold_497_pop.vcf.gz \
  -I scaffold_509_pop.vcf.gz \
  -I scaffold_516_pop.vcf.gz \
  -I scaffold_546_pop.vcf.gz \
  -I scaffold_564_pop.vcf.gz \
  -I scaffold_574_pop.vcf.gz \
  -I scaffold_91_pop.vcf.gz \
  -O population.vcf.gz

gatk IndexFeatureFile -I population.vcf.gz