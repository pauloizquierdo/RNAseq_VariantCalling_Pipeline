#!/bin/bash
#SBATCH --job-name=markdup     # Job name
#SBATCH --nodes=1              # Run all processes on a single node        
#SBATCH --ntasks=1             # Run a single task 
#SBATCH --array=1-98           # Number of tasks       
#SBATCH --cpus-per-task=10     # Number of CPU cores per task
#SBATCH --mem=22GB             # Total memory limit
#SBATCH --time=03:00:00        # Time limit hrs:min:sec
#SBATCH -A data-machine
#SBATCH --output=markdup_%A_%a.out
#SBATCH --error=markdup_%A_%a.err
#SBATCH --constraint="[intel18|amd20|amd22]"

cd /mnt/research/glbrc_group/lowry_lab/RNASeq_cistrans

READS="/mnt/gs21/scratch/seguraab/SW_Cold/00_raw_data"
GENOME_DIR="/mnt/research/glbrc_group/lowry_lab/RNASeq_cistrans/STAR_index"
MAPPING="STAR_mapping/"
TMP_DIR="/mnt/gs21/scratch/izquier7/tpm_gatk"

SAMPLE_NAME=$(awk "NR==${SLURM_ARRAY_TASK_ID}" Sample_name_AP13_DAC6_F1.txt | cut -f2)

# Load modules
module purge
module load picard/2.25.1-Java-11 


java -Xms20G -Xmx20G -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    I=${MAPPING}/${SAMPLE_NAME}_2ndPass_Aligned.sortedByCoord.out.bam \
    O=${MAPPING}/${SAMPLE_NAME}_2ndPass_Aligned.sortedByCoord.out.dedup.bam \
    M=${MAPPING}/${SAMPLE_NAME}_2ndPass_Aligned.sortedByCoord.out.dedup_metrics.txt 
    TMP_DIR=$TMP_DIR

# Save job info
scontrol show job $SLURM_JOB_ID > ${MAPPING}/${SAMPLE_NAME}_job_info_star_align_dedup_${SLURM_ARRAY_TASK_ID}.txt
