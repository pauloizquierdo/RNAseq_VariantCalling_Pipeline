#!/bin/bash
#SBATCH --job-name=STAR_index                # Job name
#SBATCH --nodes=1                           # Run all processes on a single node        
#SBATCH --ntasks=1                          # Run a single task        
#SBATCH --cpus-per-task=10                 # Number of CPU cores per task
#SBATCH --mem=50GB                         # Total memory limit
#SBATCH --time=04:00:00                     # Time limit hrs:min:sec
#SBATCH --output=STAR_index_output.txt    # Standard output log
#SBATCH -A data-machine
#SBATCH --constraint="[intel18|amd20|amd22]"

#module load Cufflinks/2.2.1
module load GCC/11.3.0
module load STAR/2.7.10b

#gffread Pvirgatumvar_AP13HAP1_772_v6.1.gene_exons.gff3 -T -o Pvirgatumvar_AP13HAP1_772_v6.1.gene_exons.gtf

STAR --runThreadN 10 \
     --runMode genomeGenerate \
     --genomeDir STAR_index \
     --genomeFastaFiles ../WGS/referenceGenome/Pvirgatumvar_AP13HAP1_772_v6.0.hardmasked.fa \
     --sjdbGTFfile ../WGS/referenceGenome/Pvirgatumvar_AP13HAP1_772_v6.1.gene_exons.gtf \
     --sjdbOverhang 149 # sjdbOverhang should be set to the length of the RNA-seq read sequences minus one, most of reads length is 150.
