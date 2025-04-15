#!/bin/bash
#SBATCH --job-name=STAR_P2     # Job name
#SBATCH --nodes=1              # Run all processes on a single node        
#SBATCH --ntasks=1             # Run a single task 
#SBATCH --array=1-98           # Number of tasks       
#SBATCH --cpus-per-task=10     # Number of CPU cores per task
#SBATCH --mem=22GB             # Total memory limit
#SBATCH --time=03:00:00        # Time limit hrs:min:sec
#SBATCH -A data-machine
#SBATCH --output=star_align_pass2_%A_%a.out
#SBATCH --error=star_align_pass2_%A_%a.err
#SBATCH --constraint="[intel18|amd20|amd22]"

cd /mnt/research/glbrc_group/lowry_lab/RNASeq_cistrans

READS="/mnt/gs21/scratch/seguraab/SW_Cold/00_raw_data"
GENOME_DIR="/mnt/research/glbrc_group/lowry_lab/RNASeq_cistrans/STAR_index"
MAPPING="STAR_mapping/"

SAMPLE_FASTA=$(awk "NR==${SLURM_ARRAY_TASK_ID}" Sample_name_AP13_DAC6_F1.txt | cut -f1)
SAMPLE_NAME=$(awk "NR==${SLURM_ARRAY_TASK_ID}" Sample_name_AP13_DAC6_F1.txt | cut -f2)

module load GCC/11.3.0
module load STAR/2.7.10b

STAR --genomeDir ${GENOME_DIR} \
     --readFilesIn ${READS}/${SAMPLE_FASTA} \
     --runThreadN 10 \
     --outFileNamePrefix ${MAPPING}/${SAMPLE_NAME}_2ndPass_ \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --sjdbFileChrStartEnd ${MAPPING}/${SAMPLE_NAME}_SJ.out.tab \
     --limitBAMsortRAM 20000000000

# Save job info
scontrol show job $SLURM_JOB_ID > ${MAPPING}/${SAMPLE_NAME}_job_info_star_align_pass2_${SLURM_ARRAY_TASK_ID}.txt
