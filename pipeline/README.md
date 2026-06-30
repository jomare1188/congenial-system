# GBS Pipeline: FASTQ to VCF for Population Analysis

Complete pipeline for processing Genotyping-by-Sequencing (GBS) data from raw FASTQ files to VCF format suitable for Structure and other population genetics analyses.

## Overview

This pipeline processes single-ended GBS reads using:
- **Enzymes**: MspI and PstI
- **Quality Control**: FastQC and MultiQC
- **Trimming**: GBStrim
- **SNP Calling**: Stacks de novo pipeline
- **Output**: VCF, Structure, Genepop formats

## Directory Structure

```
/home/diegoj/TATIANA/
├── raw/test/              # Input: FASTQ.GZ files
├── 01_fastqc/             # Quality control reports
├── 02_trimmed/            # Trimmed FASTQ files
├── 03_stacks/             # Stacks intermediate files
└── 04_results/            # Final VCF and Structure files
```

## Quick Start

### 1. Check Dependencies

```bash
bash check_dependencies.sh
```

### 2. Run the Pipeline

```bash
# Make executable
chmod +x gbs_pipeline.sh

# Run the full pipeline
bash gbs_pipeline.sh
```

### 3. Monitor Progress

The pipeline will create log files in each step. Monitor with:

```bash
tail -f /home/diegoj/TATIANA/03_stacks/denovo_map.log
tail -f /home/diegoj/TATIANA/04_results/populations.log
```

## Pipeline Steps

### Step 1: Quality Control (FastQC/MultiQC)
- Runs FastQC on all raw FASTQ files in parallel
- Aggregates results with MultiQC
- **Output**: `01_fastqc/raw_multiqc_report.html`

### Step 2: Read Trimming (GBStrim)
- Removes adapter sequences and low-quality bases
- Uses MspI and PstI enzyme recognition sites
- Crops reads to 50bp
- **Parameters**:
  - `--enzyme1 mspi`: First restriction enzyme
  - `--enzyme2 psti`: Second restriction enzyme
  - `--croplength 50`: Final read length
  - `--read R1`: Single-end reads

### Step 3: Post-Trim QC
- Quality control on trimmed reads
- **Output**: `01_fastqc/trimmed/trimmed_multiqc_report.html`

### Step 4: De Novo Assembly (Stacks)
- Builds loci de novo without reference genome
- Creates catalog of loci across samples
- **Key Parameters**:
  - `-M 2`: Mismatches between stacks within individual
  - `-n 1`: Mismatches between catalog loci
  - `-m 3`: Minimum depth of coverage
  - `-T 100`: Number of threads

### Step 5: SNP Calling and Filtering (populations)
- Calls SNPs across populations
- Applies filters for downstream analysis
- Exports multiple formats
- **Key Filters**:
  - `-r 0.8`: Locus must be present in 80% of individuals
  - `-p 1`: Minimum number of populations
  - `--min-maf 0.05`: Minimum minor allele frequency of 5%
  - `--max-obs-het 0.7`: Maximum observed heterozygosity
  - `--write-single-snp`: One SNP per locus (reduces LD)

## Output Files

### Quality Control
- `01_fastqc/raw_multiqc_report.html` - Raw data QC summary
- `01_fastqc/trimmed/trimmed_multiqc_report.html` - Trimmed data QC summary

### Final Results (04_results/)
- **populations.snps.vcf** - VCF file for general use
- **populations.structure** - Input for Structure software
- **populations.genepop** - Input for Genepop software
- **populations.fst_*.tsv** - FST statistics between populations
- **populations.sumstats*.tsv** - Summary statistics per locus
- **pipeline_summary.txt** - Complete pipeline summary

## Parameter Tuning



### GBStrim Parameters

**Current settings are optimized for MspI/PstI combination:**
- `--croplength 50`: Adjust based on your read length and quality

### Stacks denovo_map.pl Parameters

**For higher stringency (fewer, higher quality SNPs):**
```bash
-M 1   # Allow fewer mismatches within individual
-n 0   # Require exact matches between catalog loci
-m 5   # Increase minimum depth
```

**For more permissive calling (more SNPs, potentially more error):**
```bash
-M 3   # Allow more mismatches
-n 2   # Allow more catalog mismatches
-m 2   # Lower minimum depth
```

### Stacks populations Parameters

**For more stringent filtering:**
```bash
-r 0.9          # Locus present in 90% of samples
--min-maf 0.1   # Higher MAF threshold
--max-obs-het 0.5  # Stricter heterozygosity filter
```

**For more permissive filtering:**
```bash
-r 0.5          # Locus present in 50% of samples
--min-maf 0.01  # Lower MAF threshold
--max-obs-het 0.9  # More permissive
```

## Population Map File

The pipeline creates a default population map (`03_stacks/popmap.txt`) with all samples in one population. 

**You MUST edit this file** to reflect your actual population structure:

```
sample1    pop1
sample2    pop1
sample3    pop2
sample4    pop2
```

Format: `sample_name<TAB>population_id`

## References

- **Stacks**: Rochette NC, et al. (2019) Stacks 2: Analytical framework for paired-end sequencing improves RADseq-based population genomics. Molecular Ecology 28:4737-4754
- **GBStrim**: https://bitbucket.org/jgarbe/gbstrim/
- **Structure**: Pritchard JK, et al. (2000) Genetics 155:945-959

## Support

- **Stacks Manual**: http://catchenlab.life.illinois.edu/stacks/manual/
- **Stacks Google Group**: https://groups.google.com/g/stacks-users

## Notes

- Pipeline designed for single-ended reads only
- Uses 100 threads - adjust based on your system
- All intermediate files are kept for troubleshooting
- Compressed FASTQ files (.gz) are handled natively
