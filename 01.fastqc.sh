#!/bin/bash
#SBATCH --job-name=fastqc                # Job name
#SBATCH --array=1-295                   # Array range
#SBATCH --nodes=1                           # Run all processes on a single node        
#SBATCH --ntasks=1                          # Run a single task        
#SBATCH --cpus-per-task=5                 # Number of CPU cores per task
#SBATCH --mem=20GB                         # Total memory limit
#SBATCH --time=04:00:00                     # Time limit hrs:min:sec
#SBATCH --output=fastqc_%A_%a.txt    # Standard output log
#SBATCH -A data-machine
#SBATCH --constraint="[intel18|amd20|amd22]"

# Load the FastQC module
module load FastQC/0.11.7-Java-1.8.0_162

# Directory containing  .fastq.gz files
FASTQ_DIR="/mnt/gs21/scratch/seguraab/SW_Cold/00_raw_data"

# Output directory for FastQC reports
OUTPUT_DIR="/mnt/research/glbrc_group/lowry_lab/RNASeq_cistrans/fastqc_files"

# Create an array of all .fastq.gz files in the directory
FASTQ_FILES=($FASTQ_DIR/*.fastq.gz)

# Use SLURM_ARRAY_TASK_ID to get the file for this task
FILE=${FASTQ_FILES[$SLURM_ARRAY_TASK_ID]}

# Run FastQC on the selected file
fastqc -o "$OUTPUT_DIR" -t 5 --noextract "$FILE"
