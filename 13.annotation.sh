#!/bin/bash
#SBATCH --job-name=annotation                # Job name
#SBATCH --nodes=1                           # Run all processes on a single node       
#SBATCH --ntasks=1                          # Run a single task        
#SBATCH --cpus-per-task=10                 # Number of CPU cores per task
#SBATCH --mem=50GB                          # Total memory limit
#SBATCH --time=5-00:00:00                     # Time limit hrs:min:sec
#SBATCH -A data-machine
#SBATCH --output=annotation_snp.out
#SBATCH --error=annotation_snp.err
#SBATCH --constraint="[intel18|amd20|amd22]"

cd /mnt/research/glbrc_group/lowry_lab/RNASeq_cistrans/population

java -jar /mnt/home/izquier7/Documents/software/NGSEPcore_4.3.1.jar VCFAnnotate \
    -i population_snps_Q40_DP20_FM01.vcf \
    -r /mnt/research/glbrc_group/lowry_lab/WGS/referenceGenome/Pvirgatumvar_AP13HAP1_772_v6.0.hardmasked.fa \
    -t /mnt/research/glbrc_group/lowry_lab/WGS/referenceGenome/Pvirgatumvar_AP13HAP1_772_v6.1.gene_exons.gff3 \
    -o population_snps_Q40_DP20_FM01_annotated.vcf

