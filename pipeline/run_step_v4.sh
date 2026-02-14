#!/bin/bash

################################################################################
# GBS Pipeline - Individual Step Runner
# Run specific steps of the pipeline independently
################################################################################

set -euo pipefail

# Configuration
RUN_NAME="run2"
THREADS=100
RAW_DIR="/dados01/jorge/TATIANA/raw"
BASE_DIR="/dados01/jorge/TATIANA"
QC_DIR="${BASE_DIR}/run2/01_fastqc"
TRIMMED_DIR="${BASE_DIR}/run2/02_trimmed"
RESYNCED_DIR="${BASE_DIR}/run2/02b_resynced"
STACKS_DIR="${BASE_DIR}/${RUN_NAME}/03_stacks"
GBSTRIM_DIR="/dados01/jorge/TATIANA/bin/gbstrim"
RESULTS_DIR="${BASE_DIR}/${RUN_NAME}/04_results"


# Function definitions
run_fastqc() {
    echo "Running FastQC on raw reads..."
    mkdir -p "${QC_DIR}"
    find "${RAW_DIR}" -name "*.fastq.gz" | \
        parallel -j "${THREADS}" \
        "fastqc {} -o ${QC_DIR} --quiet"
    multiqc "${QC_DIR}" -o "${QC_DIR}" -n raw_multiqc_report --quiet --force
    echo "✓ FastQC complete. Report: ${QC_DIR}/raw_multiqc_report.html"
}

run_gbstrim() {
    echo "Running GBStrim..."
    mkdir -p "${TRIMMED_DIR}"

    # Get the directory where this script is located
    SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

    find "${RAW_DIR}" -name "*_R1_001.fastq.gz" | \
        parallel -j "${THREADS}" \
        "${SCRIPT_DIR}/gbstrim_wrapper.sh {} ${TRIMMED_DIR} R1" >> logs_gbstrim.txt 2>&1

    echo "✓ GBStrim R1 complete. Output: ${TRIMMED_DIR}"

    find "${RAW_DIR}" -name "*_R2_001.fastq.gz" | \
        parallel -j "${THREADS}" \
        "${SCRIPT_DIR}/gbstrim_wrapper.sh {} ${TRIMMED_DIR} R2" >> logs_gbstrim.txt 2>&1

    echo "✓ GBStrim R2 complete. Output: ${TRIMMED_DIR}"
}

run_resync() {
    echo "Running resync.pl to synchronize paired-end reads..."
    mkdir -p "${RESYNCED_DIR}"

    # Get the directory where gbstrim tools are located
    SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

    # Export variables for parallel
    export TRIMMED_DIR
    export RESYNCED_DIR
    export SCRIPT_DIR

    # Create a wrapper function for parallel
    resync_pair() {
        R1_FILE="$1"
        R2_FILE="${R1_FILE/_R1_/_R2_}"

        if [ ! -f "${R2_FILE}" ]; then
            echo "⚠ WARNING: R2 file not found for ${R1_FILE}"
            return 1
        fi

        # Get base sample name
        SAMPLE=$(basename "${R1_FILE}" | sed 's/_R1_.*trimmed\.fastq$//')

        # Output files
        R1_SYNCED="${RESYNCED_DIR}/${SAMPLE}.1.fastq"
        R2_SYNCED="${RESYNCED_DIR}/${SAMPLE}.2.fastq"

        echo "Resyncing: ${SAMPLE}"

        # Run resync.pl (assuming it's in same directory as this script or in PATH)
        if command -v /dados01/jorge/TATIANA/bin/gbstrim/resync.pl &> /dev/null; then
            /dados01/jorge/TATIANA/bin/gbstrim/resync.pl "${R1_FILE}" "${R2_FILE}" "${R1_SYNCED}" "${R2_SYNCED}" > "${RESYNCED_DIR}/"${SAMPLE}"_resync.log" 2>&1
        elif [ -f "${SCRIPT_DIR}/resync.pl" ]; then
            "${SCRIPT_DIR}/resync.pl" "${R1_FILE}" "${R2_FILE}" "${R1_SYNCED}" "${R2_SYNCED}" 2>&1
        else
            echo "ERROR: resync.pl not found in PATH or ${SCRIPT_DIR}"
            return 1
        fi
    }

    # Export the function for parallel
    export -f resync_pair

    # Run parallel on all R1 files
    find "${TRIMMED_DIR}" -name "*_R1_*trimmed.fastq" -type f | \
        parallel -j "${THREADS}" resync_pair {} > "${RESYNCED_DIR}/resync.log" 2>&1

    # Count successful outputs
    RESYNC_COUNT=$(find "${RESYNCED_DIR}" -name "*.1.fastq" -type f | wc -l)

    if [ ${RESYNC_COUNT} -eq 0 ]; then
        echo "ERROR: No paired files found to resync!"
        echo "Check ${RESYNCED_DIR}/resync.log for details"
        exit 1
    fi

    echo "✓ Resync complete. Processed ${RESYNC_COUNT} sample pairs"
    echo "  Output: ${RESYNCED_DIR}"
    echo "  Log: ${RESYNCED_DIR}/resync.log"
}

run_fastqc_trimmed() {
    echo "Running FastQC on resynced reads..."
    QC_TRIMMED_DIR="${QC_DIR}/resynced"
    mkdir -p "${QC_TRIMMED_DIR}"
    find "${RESYNCED_DIR}" -name "*.fastq" | \
        parallel -j "${THREADS}" \
        "fastqc {} -o ${QC_TRIMMED_DIR} --quiet"
    multiqc "${QC_TRIMMED_DIR}" -o "${QC_TRIMMED_DIR}" -n resynced_multiqc_report --quiet --force
    echo "✓ FastQC (resynced) complete. Report: ${QC_TRIMMED_DIR}/resynced_multiqc_report.html"
}


run_stacks() {
    echo "Running Stacks de novo pipeline..."
    mkdir -p "${STACKS_DIR}"

    POPMAP="${STACKS_DIR}/popmap.txt"
    if [ ! -f "${POPMAP}" ]; then
        echo "Creating default population map..."
        # Use resynced R1 files for popmap
        for file in "${RESYNCED_DIR}"/*.1.fastq; do
            sample=$(basename "${file}" .1.fastq)
            echo -e "${sample}\t1" >> "${POPMAP}"
        done
        echo "⚠ WARNING: All samples assigned to population 1"
        echo "  Edit ${POPMAP} before proceeding if you have multiple populations"
    fi

    denovo_map.pl \
        -M 1 \
        -n 1 \
        -m 3 \
        -T "${THREADS}" \
        -o "${STACKS_DIR}" \
        --samples "${RESYNCED_DIR}" \
        --popmap "${POPMAP}" \
        --paired \
        2>&1 | tee "${STACKS_DIR}/denovo_map.log"

    echo "✓ Stacks complete. Output: ${STACKS_DIR}"
}

run_populations() {
    echo "Running populations to generate VCF..."
    mkdir -p "${RESULTS_DIR}"

    POPMAP="${STACKS_DIR}/popmap.txt"
    if [ ! -f "${POPMAP}" ]; then
        echo "ERROR: Population map not found: ${POPMAP}"
        echo "Run step 5 (stacks) first"
        exit 1
    fi

    populations \
        -P "${STACKS_DIR}" \
        -M "${POPMAP}" \
        -t "${THREADS}" \
        --min-mac 2 \
        -r 0.65 \
        -p 1 \
        --min-maf 0.05 \
        --max-obs-het 0.7 \
        --write-single-snp \
        --vcf \
        --structure \
        --genepop \
        --fstats \
        -O "${RESULTS_DIR}" \
        2>&1 | tee "${RESULTS_DIR}/populations.log"

    if [ -f "${RESULTS_DIR}/populations.snps.vcf" ]; then
        N_SNPS=$(grep -v "^#" "${RESULTS_DIR}/populations.snps.vcf" | wc -l)
        echo "✓ Populations complete. SNPs found: ${N_SNPS}"
        echo "  VCF: ${RESULTS_DIR}/populations.snps.vcf"
        echo "  Structure: ${RESULTS_DIR}/populations.structure"
    fi
}

run_variant_polish() {
    echo "Running VCF filtering and polishing..."
    POLISH_DIR="${BASE_DIR}/${RUN_NAME}/05_polish"
    mkdir -p "${POLISH_DIR}"

    VCF_INPUT="${RESULTS_DIR}/populations.snps.vcf"
    
    if [ ! -f "${VCF_INPUT}" ]; then
        echo "ERROR: VCF file not found: ${VCF_INPUT}"
        echo "Run step 6 (populations) first"
        exit 1
    fi

    # Calculate max-meanDP (2x the mean depth across all samples)
    echo "Calculating mean depth threshold for paralog detection..."
    MEAN_DP=$(bcftools query -f '%INFO/DP\n' "${VCF_INPUT}" | \
              awk '{sum+=$1; count++} END {print int(sum/count*2)}')
    
    echo "  Mean depth threshold (2x): ${MEAN_DP}"

    # Apply filters
    VCF_FILTERED="${POLISH_DIR}/populations.filtered.vcf"
    
    vcftools \
        --vcf "${VCF_INPUT}" \
        --minDP 8 \
        --maxDP 100 \
        --hwe 0.01 \
        --max-meanDP 100\
        --recode \
        --recode-INFO-all \
        --out "${POLISH_DIR}/populations.filtered" \
        2>&1 | tee "${POLISH_DIR}/vcftools.log"

    # Rename output file
    if [ -f "${POLISH_DIR}/populations.filtered.recode.vcf" ]; then
        mv "${POLISH_DIR}/populations.filtered.recode.vcf" "${VCF_FILTERED}"
        
        N_SNPS_BEFORE=$(grep -v "^#" "${VCF_INPUT}" | wc -l)
        N_SNPS_AFTER=$(grep -v "^#" "${VCF_FILTERED}" | wc -l)
        
        echo "✓ VCF filtering complete"
        echo "  SNPs before filtering: ${N_SNPS_BEFORE}"
        echo "  SNPs after filtering: ${N_SNPS_AFTER}"
        echo "  Filtered VCF: ${VCF_FILTERED}"
        echo "  Log: ${POLISH_DIR}/vcftools.log"
    else
        echo "ERROR: Filtering failed. Check ${POLISH_DIR}/vcftools.log"
        exit 1
    fi
}


# Show usage
show_usage() {
    cat <<EOF
GBS Pipeline - Individual Step Runner (with Resync)

Usage: $0 [step_number|step_name|all]

Steps:
  1 | fastqc          - Run FastQC on raw reads
  2 | gbstrim         - Trim reads with GBStrim
  3 | resync          - Resynchronize paired-end reads
  4 | fastqc_trimmed  - Run FastQC on resynced reads
  5 | stacks          - Run Stacks de novo pipeline (paired-end mode)
  6 | populations     - Generate VCF and other formats
  7 | variant polish  - Filter VCF (depth, HWE, paralogs)
  all                 - Run all steps in sequence

Examples:
  $0 1              # Run only FastQC
  $0 resync         # Run only resync
  $0 stacks         # Run only Stacks
  $0 all            # Run complete pipeline

Configuration:
  Threads: ${THREADS}
  Raw data: ${RAW_DIR}
  Output: ${BASE_DIR}
  Resynced reads: ${RESYNCED_DIR}

Note: Stacks now runs in paired-end mode using resynced reads
EOF
}

# Main execution
if [ $# -eq 0 ]; then
    show_usage
    exit 0
fi

case "$1" in
    1|fastqc)
        run_fastqc
        ;;
    2|gbstrim)
        run_gbstrim
        ;;
    3|resync)
        run_resync
        ;;
    4|fastqc_trimmed)
        run_fastqc_trimmed
        ;;
    5|stacks)
        run_stacks
        ;;
    6|populations)
        run_populations
        ;;
    7|variant_polish)
	run_variant_polish
	;;
    all)
        echo "Running complete pipeline..."
        run_fastqc
        run_gbstrim
        run_resync
        run_fastqc_trimmed
        run_stacks
        run_populations
	run_variant_polish
        echo ""
        echo "✓ Complete pipeline finished successfully!"
        ;;
    *)
        echo "ERROR: Unknown step '$1'"
        echo ""
        show_usage
        exit 1
        ;;
esac
