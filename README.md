Notes:
  - The FASTQ files used in this pipeline are located in the folder:
  /mnt/gs21/scratch/seguraab/SW_Cold/00_raw_data.

  - The FASTQ files have the suffix .filter, which suggests that the reads were filtered for quality control. However, I could not find additional information to confirm this.

The pipeline described in this repository follows the best practices outlined in the following link:
https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels
________________________________________________________  
  Scripts:
1.	fastqc: Quality Check
  -	FASTQC was run on all FASTQ files to assess their quality.
  -	All files had 0 sequences flagged as poor quality, and the sequence lengths ranged from 49–150 bp, suggesting that the reads were already clean.
2.	genome_index_star: Genome Index Creation
  - Pvirgatumvar_AP13HAP1_772_v6.0.hardmasked.fa
  -	This script creates the genome index required for the STAR aligner.
3.	mapping_star: RNA Read Mapping
  -	Maps RNA-seq reads to a reference genome using the STAR aligner.
4.	mapping_star_pass2: Two-Pass Mapping
  -	Implements STAR’s two-pass mode, as recommended by GATK, to improve alignment accuracy around novel splice junctions.
5.	markduplicates: Duplicate Read Identification
  -	Identifies and marks duplicate reads in the BAM file.
6.	SplitNCigarReads: Splitting Reads with N in CIGAR Strings
  -	Splits reads containing N in the CIGAR string into multiple supplementary alignments.
  -	Hard-clips mismatching overhangs and adjusts mapping qualities to follow DNA conventions.
7.	NamesortBam: Name-Sorting BAM Files
  -	The outputs from SplitNCigarReads are sorted and indexed for downstream processing.
8.	HaplotypeCaller: Variant Calling
  -	Identifies variants in the RNA-seq data using GATK’s HaplotypeCaller.
9.	GenomicsDBImport: Database Creation
  -	Imports single-sample GVCFs (per chromosome or scaffold) into GenomicsDB for joint genotyping.
10.	Genotype: Joint Genotyping
  -	Performs joint genotyping on pre-called samples to produce a multi-sample VCF file.
11.	concatenate: Merge Variants
  -	Combines variant files from all chromosomes or scaffolds into a single VCF file.
12.	bcf_pruning: Variant Filtering
  -	Filters SNPs and indels, retaining only variants with:
    -	Quality > 40
    -	Depth of coverage > 20
    -	Missing data < 0.1.
13.	annotation: Variant Annotation
  -	Annotates the final VCF file using the Pvirgatumvar_AP13HAP1_772_v6.0 GFF3 file.
