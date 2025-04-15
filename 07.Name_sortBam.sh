#!/bin/bash
#SBATCH --job-name=index     # Job name
#SBATCH --nodes=1              # Run all processes on a single node        
#SBATCH --ntasks=1             # Run a single task 
#SBATCH --array=1-98           # Number of tasks       
#SBATCH --cpus-per-task=4     # Number of CPU cores per task
#SBATCH --mem=44GB             # Total memory limit
#SBATCH --time=02:00:00        # Time limit hrs:min:sec
#SBATCH -A data-machine
#SBATCH --output=index_%A_%a.out
#SBATCH --error=index_%A_%a.err
#SBATCH --constraint="[intel18|amd20|amd22]"

cd /mnt/research/glbrc_group/lowry_lab/RNASeq_cistrans

GENOME_DIR="/mnt/research/glbrc_group/lowry_lab/WGS/referenceGenome"
MAPPING="STAR_mapping/"
TMP_DIR="/mnt/gs21/scratch/izquier7/tpm_gatk"

SAMPLE_NAME=$(awk "NR==${SLURM_ARRAY_TASK_ID}" Sample_name_AP13_DAC6_F1.txt | cut -f2)

Load modules
module purge
module load picard/2.25.1-Java-11 

# GATK require read groups for proper variant calling and base recalibration, and after duplication the RGPL was misssing, replace with illumina
java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
      I=${MAPPING}/${SAMPLE_NAME}_2ndPass_Aligned.sortedByCoord.out.dedup.split.bam \
      O=${MAPPING}/${SAMPLE_NAME}_2ndPass_Aligned.sortedByCoord.out.dedup.split_rg.bam \
      RGID=${SAMPLE_NAME} \
      RGLB=${SAMPLE_NAME} \
      RGPL=ILLUMINA \
      RGPU=${SAMPLE_NAME} \
      RGSM=${SAMPLE_NAME}

# There is a conflict between picard and samtools modules, so we need to purge and load samtools after picard
module purge
module load GCC/11.3.0
module load SAMtools/1.16.1

samtools index ${MAPPING}/${SAMPLE_NAME}_2ndPass_Aligned.sortedByCoord.out.dedup.split_rg.bam

# Save job info
scontrol show job $SLURM_JOB_ID > ${MAPPING}/${SAMPLE_NAME}_job_info_index_${SLURM_ARRAY_TASK_ID}.txt
