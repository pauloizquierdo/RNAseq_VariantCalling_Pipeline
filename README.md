# RNA-seq Variant Calling Pipeline (GATK Best Practices)

This repository contains a complete RNA-seq variant discovery pipeline based on the GATK Best Practices for SNP and Indel calling from short-read RNA-seq data.

🔗 Reference:  
[GATK RNA-seq Best Practices Guide](https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels)

---

## 📁 Input Data

The FASTQ files used in this pipeline are located at:
/mnt/gs21/scratch/seguraab/SW_Cold/00_raw_data

---

## 🔬 Pipeline Overview

The steps below follow GATK's best practices for RNA-seq variant discovery:

### 1. `fastqc`: Quality Control
- FASTQC was run on all FASTQ files to evaluate read quality.
- No sequences were flagged as poor quality.
- Sequence lengths ranged from 49–150 bp, indicating pre-cleaned reads.

### 2. `genome_index_star`: STAR Genome Indexing
- Creates a STAR genome index using:
Pvirgatumvar_AP13HAP1_772_v6.0.hardmasked.fa

### 3. `mapping_star`: First-Pass Read Mapping
- Aligns RNA-seq reads to the reference genome using STAR.

### 4. `mapping_star_pass2`: Two-Pass Mapping
- Implements STAR's two-pass mapping mode to improve alignment across novel splice junctions, as recommended by GATK.

### 5. `markduplicates`: Mark Duplicates
- Uses Picard to identify and mark PCR duplicates in BAM files.

### 6. `SplitNCigarReads`: Read Processing
- Splits reads with `N` operators in the CIGAR string into separate exonic segments.
- Hard-clips mismatching overhangs.
- Adjusts mapping qualities to conform to DNA variant calling conventions.

### 7. `namesortBam`: Name Sort and Index
- BAM files from `SplitNCigarReads` are name-sorted and indexed for use with GATK tools.

### 8. `HaplotypeCaller`: Variant Calling
- Calls SNPs and indels using GATK's `HaplotypeCaller` in GVCF mode.

### 9. `GenomicsDBImport`: GVCF Consolidation
- Combines individual sample GVCFs (by chromosome or scaffold) into a GenomicsDB database.

### 10. `Genotype`: Joint Genotyping
- Performs joint genotyping across all samples to generate a multi-sample VCF.

### 11. `concatenate`: Variant Merging
- Merges VCFs from all scaffolds or chromosomes into a single file.

### 12. `bcf_pruning`: Variant Filtering
- Applies filters to retain high-confidence variants:
- Quality > 40
- Depth > 20
- Missing data < 10%

### 13. `annotation`: Functional Annotation
- Annotates variants using the `Pvirgatumvar_AP13HAP1_772_v6.0.gff3` annotation file.

---

## 🧬 Notes

- All reference files used (e.g., reference genome, GFF3) are from the **P. virgatum AP13 v6.0** assembly.
- This pipeline is designed for **RNA-seq-based variant discovery** and may require modification for DNA-seq or other data types.

---

## ✅ Output

- Cleaned and aligned BAM files
- Multi-sample VCF with high-confidence SNPs and indels
- Annotated variants linked to gene features

---





