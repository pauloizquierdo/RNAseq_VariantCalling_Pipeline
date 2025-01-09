#!/bin/bash
#SBATCH --job-name=GVCF                # Job name
#SBATCH --nodes=1                           # Run all processes on a single node 
#SBATCH --array=1-46       
#SBATCH --ntasks=1                          # Run a single task        
#SBATCH --cpus-per-task=10                 # Number of CPU cores per task
#SBATCH --mem=50GB                          # Total memory limit
#SBATCH --time=04:00:00                     # Time limit hrs:min:sec
#SBATCH -A glbrc
#SBATCH --output=GVCF_%A_%a.out
#SBATCH --error=GVCF_%A_%a.err
#SBATCH --constraint="[intel18|amd20|amd22]"

cd /mnt/research/glbrc_group/lowry_lab/RNASeq_cistrans/STAR_mapping

TMP_DIR="/mnt/gs21/scratch/izquier7/tpm_gatk"
GDBI="GenomicsDB/"

CHR=$(awk "NR==${SLURM_ARRAY_TASK_ID}" ../../WGS/referenceGenome/Pvirgatumvar_AP13HAP1_772_v6.0.hardmasked.fa.fai | cut -f1)

CHROMOSOME_DB_PATH="${GDBI}/${CHR}"

# Load modules
module purge
module load GCC/7.3.0-2.30  OpenMPI/3.1.1
module load GATK/4.1.4.1-Python-3.6.6

gatk --java-options "-Xms40g -Xmx40g -XX:ParallelGCThreads=2" GenomicsDBImport \
  --sample-name-map sample_map.txt \
  --genomicsdb-workspace-path $CHROMOSOME_DB_PATH \
  --tmp-dir $TMP_DIR \
  -L $CHR \
  --reader-threads 8

scontrol show job > GenomicsDBImport_L_job_info.txt