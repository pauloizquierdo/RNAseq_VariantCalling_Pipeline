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

REFGEN_DIR="/mnt/research/glbrc_group/lowry_lab/WGS/referenceGenome" # reference genome directory
TMP_DIR="/mnt/gs21/scratch/izquier7/tpm_gatk"
POP_DIR="/mnt/research/glbrc_group/lowry_lab/RNASeq_cistrans/population"

CHR=$(awk "NR==${SLURM_ARRAY_TASK_ID}" ../../WGS/referenceGenome/Pvirgatumvar_AP13HAP1_772_v6.0.hardmasked.fa.fai | cut -f1)

# Load modules
module purge
module load GCC/7.3.0-2.30  OpenMPI/3.1.1
module load GATK/4.1.4.1-Python-3.6.6

gatk --java-options "-Xms35G -Xmx35G -XX:ParallelGCThreads=2" GenotypeGVCFs \
    -R $REFGEN_DIR/Pvirgatumvar_AP13HAP1_772_v6.0.hardmasked.fa \
    -V gendb://GenomicsDB/$CHR \
    -L $CHR \
    --tmp-dir $TMP_DIR \
    -O $POP_DIR/${CHR}_pop.vcf.gz

scontrol show job $SLURM_JOB_ID > ${CHR}_${SLURM_ARRAY_TASK_ID}_hardmasked.txt