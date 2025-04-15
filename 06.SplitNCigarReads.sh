#!/bin/bash
#SBATCH --job-name=splitN     # Job name
#SBATCH --nodes=1              # Run all processes on a single node        
#SBATCH --ntasks=1             # Run a single task 
#SBATCH --array=1-98           # Number of tasks       
#SBATCH --cpus-per-task=10     # Number of CPU cores per task
#SBATCH --mem=22GB             # Total memory limit
#SBATCH --time=03:00:00        # Time limit hrs:min:sec
#SBATCH -A data-machine
#SBATCH --output=splitN_%A_%a.out
#SBATCH --error=splitN_%A_%a.err
#SBATCH --constraint="[intel18|amd20|amd22]"

cd /mnt/research/glbrc_group/lowry_lab/RNASeq_cistrans


GENOME_DIR="/mnt/research/glbrc_group/lowry_lab/WGS/referenceGenome"
MAPPING="STAR_mapping/"

SAMPLE_NAME=$(awk "NR==${SLURM_ARRAY_TASK_ID}" Sample_name_AP13_DAC6_F1.txt | cut -f2)

# Load modules
module purge
module load GCC/7.3.0-2.30  OpenMPI/3.1.1
module load GATK/4.1.4.1-Python-3.6.6

gatk SplitNCigarReads \
    -R ${GENOME_DIR}/Pvirgatumvar_AP13HAP1_772_v6.0.hardmasked.fa \
    -I ${MAPPING}/${SAMPLE_NAME}_2ndPass_Aligned.sortedByCoord.out.dedup.bam \
    -O ${MAPPING}/${SAMPLE_NAME}_2ndPass_Aligned.sortedByCoord.out.dedup.split.bam 

# Save job info
scontrol show job $SLURM_JOB_ID > ${MAPPING}/${SAMPLE_NAME}_job_info_star_align_dedup_split_${SLURM_ARRAY_TASK_ID}.txt
