#!/bin/bash

set -euo pipefail

# Default values
BASE_DIR=""
OUTPUT_DIR=""

# Help message
usage() {
    echo "Usage: $0 -b BASE_DIR -o OUTPUT_DIR"
    echo ""
    echo "Options:"
    echo "  -b BASE_DIR      Path to base directory containing species subfolders with genome sequence and SRR data"
    echo "  -o OUTPUT_DIR    Path to directory where results will be written"
    echo "  -h               Show this help message"
    echo ""
    echo "Example:"
    echo "  ./hisat2_mapping.sh -b /vol/hisat2-mapping/ -o /vol/hisat2-mapping/output"
    exit 1
}

# Parse command-line options
while getopts "b:o:h" opt; do
    case $opt in
        b) BASE_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Check required arguments
if [[ -z "$BASE_DIR" || -z "$OUTPUT_DIR" ]]; then
    usage
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Loop over species directories
for SPECIES_DIR in "$BASE_DIR"/*/; do
    # Extract species name
    SPECIES_NAME=$(basename "$SPECIES_DIR")

    # Find reference genome (.fna, .fa, .fasta)
    REF_GENOME=$(find "$SPECIES_DIR" -maxdepth 1 \( -name "*.fna" -o -name "*.fa" -o -name "*.fasta" \) | head -n 1)

    if [[ -z "$REF_GENOME" ]]; then
        echo "[INFO] No reference genome found for $SPECIES_NAME. Skipping."
        continue
    fi

    # Build HISAT2 index if needed
    INDEX_PREFIX="$SPECIES_DIR/${SPECIES_NAME}_hisat2_index"
    if [[ ! -f "${INDEX_PREFIX}.1.ht2" ]]; then
        echo "[INFO] Building HISAT2 index for $SPECIES_NAME..."
        hisat2-build "$REF_GENOME" "$INDEX_PREFIX"
    fi

    # Create species output subdirectory
    SPECIES_OUTPUT_DIR="$OUTPUT_DIR/${SPECIES_NAME}_hisat2"
    mkdir -p "$SPECIES_OUTPUT_DIR"

    BAM_FILES=()

    # Loop over RNA-seq read directories
    for READS_DIR in "$SPECIES_DIR"/*/; do
        for SAMPLE_DIR in $(find "$READS_DIR" -mindepth 1 -maxdepth 4 -type d); do
            R1=$(find "$SAMPLE_DIR" \( -name "*_1.fastq" -o -name "*_1.fastq.gz" \) | head -n 1)
            R2=$(find "$SAMPLE_DIR" \( -name "*_2.fastq" -o -name "*_2.fastq.gz" \) | head -n 1)

            if [[ -z "$R1" || -z "$R2" ]]; then
                echo "[INFO] Paired-end reads not found in $SAMPLE_DIR. Skipping."
                continue
            fi

            SAMPLE_NAME=$(basename "$SAMPLE_DIR")
            TEMP_BAM="$SPECIES_OUTPUT_DIR/${SAMPLE_NAME}_${SPECIES_NAME}.bam"

            echo "[INFO] Aligning $SAMPLE_NAME for species $SPECIES_NAME..."
            hisat2 -x "$INDEX_PREFIX" -1 "$R1" -2 "$R2" -p 4 | samtools view -Sb - > "$TEMP_BAM"

            BAM_FILES+=("$TEMP_BAM")
        done
    done

    if [[ ${#BAM_FILES[@]} -eq 0 ]]; then
        echo "[INFO] No BAM files generated for $SPECIES_NAME. Skipping."
        continue
    fi

    MERGED_BAM="$SPECIES_OUTPUT_DIR/${SPECIES_NAME}_merged.bam"
    SORTED_BAM="$SPECIES_OUTPUT_DIR/${SPECIES_NAME}_sorted.bam"

    echo "[INFO] Merging BAM files for species $SPECIES_NAME..."
    samtools merge -@ 4 "$MERGED_BAM" "${BAM_FILES[@]}"

    echo "[INFO] Sorting BAM file for species $SPECIES_NAME..."
    samtools sort -@ 4 "$MERGED_BAM" -o "$SORTED_BAM"

    echo "[INFO] Indexing sorted BAM file for species $SPECIES_NAME..."
    samtools index "$SORTED_BAM"

    echo "[INFO] Cleaning up intermediate files..."
    rm "${BAM_FILES[@]}" "$MERGED_BAM"

    echo "[INFO] Processing of species $SPECIES_NAME complete."
done

echo "[INFO] All species processed."
