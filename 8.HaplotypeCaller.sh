#!/bin/bash
#SBATCH --job-name=Haploty     # Job name
#SBATCH --nodes=1              # Run all processes on a single node        
#SBATCH --ntasks=1             # Run a single task 
#SBATCH --array=1-98           # Number of tasks       
#SBATCH --cpus-per-task=8     # Number of CPU cores per task
#SBATCH --mem=44GB             # Total memory limit
#SBATCH --time=12:00:00        # Time limit hrs:min:sec
#SBATCH -A data-machine
#SBATCH --output=Haplotychr1N-3N_%A_%a.out
#SBATCH --error=Haplotychr1N-3N_%A_%a.err
#SBATCH --constraint="[intel18|amd20|amd22]"

cd /mnt/research/glbrc_group/lowry_lab/RNASeq_cistrans

GENOME_DIR="/mnt/research/glbrc_group/lowry_lab/WGS/referenceGenome"
MAPPING="STAR_mapping/"
TMP_DIR="/mnt/gs21/scratch/izquier7/tpm_gatk"

SAMPLE_NAME=$(awk "NR==${SLURM_ARRAY_TASK_ID}" Sample_name_AP13_DAC6_F1.txt | cut -f2)

# Load modules
module purge
module load GCC/7.3.0-2.30  OpenMPI/3.1.1
module load GATK/4.1.4.1-Python-3.6.6

gatk --java-options "-Xms40g -Xmx40g -XX:ParallelGCThreads=2" HaplotypeCaller \
    -R ${GENOME_DIR}/Pvirgatumvar_AP13HAP1_772_v6.0.hardmasked.fa \
    -I ${MAPPING}/${SAMPLE_NAME}_2ndPass_Aligned.sortedByCoord.out.dedup.split_rg.bam  \
    -O ${MAPPING}/${SAMPLE_NAME}_2ndPass_Aligned.sortedByCoord.out.dedup.split_rg.vcf.gz \
    -ERC GVCF \
    --tmp-dir $TMP_DIR \

# Save job info
scontrol show job $SLURM_JOB_ID > ${MAPPING}/${SAMPLE_NAME}_job_info_haplotypecall${SLURM_ARRAY_TASK_ID}.txt

zcat STAR_mapping/DAC6_10AS_D1_2ndPass_Aligned.sortedByCoord.out.dedup.split_rg.vcf.gz | tail
zcat STAR_mapping/AP13_3AS_D3_2ndPass_Aligned.sortedByCoord.out.dedup.split_rg.vcf.gz | tail
zcat STAR_mapping/F1_3AS_D3_2ndPass_Aligned.sortedByCoord.out.dedup.split_rg.vcf.gz | tail

