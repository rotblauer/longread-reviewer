#!/usr/bin/env bash
# Download and extract NA19240 long-read data from public sources for testing.
#
# Source: IGSR (International Genome Sample Resource) / Human Pangenome Reference Consortium
#   - Sample: NA19240 (Yoruba, 1000 Genomes)
#   - Technology: PacBio HiFi (CCS), aligned to GRCh38
#   - Reference: GRCh38 (hg38)
#   - Region: chr17:10958130-11017414
#
# Uses the HPRC NA19240 aligned BAM which is indexed on S3, allowing
# efficient random access to extract just the region of interest.
#
# Requirements: samtools >= 1.10
#
# Usage:
#   ./download_from_source.sh [output_dir]
#
# The current test data included in this repo is a small synthetic BAM that mimics
# the characteristics of real NA19240 HiFi data, suitable for unit and integration testing.
# Run this script to replace it with real data for manual inspection and validation.

set -euo pipefail

OUTDIR="${1:-$(dirname "$0")}"
CHROM="chr17"
REGION="${CHROM}:10958130-11017414"
SAMPLE="NA19240"

# HPRC NA19240 PacBio HiFi aligned to GRCh38 with winnowmap (has .bai index)
# Source: https://s3-us-west-2.amazonaws.com/human-pangenomics/ -> working/HPRC_PLUS/NA19240/
HPRC_BASE="https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/NA19240/analysis/aligned_reads/hifi/GRCh38"
HPRC_BAM="${HPRC_BASE}/NA19240_aligned_GRCh38_winnowmap.sorted.bam"
HPRC_BAI="${HPRC_BAM}.bai"

# Fallback: 1000 Genomes high-coverage Illumina data (has index, good for testing)
# Note: This is Illumina short-read data, not HiFi, but useful for testing basic functionality
FALLBACK_BAM="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/data/YRI/NA19240/alignment/NA19240.alt_bwamem_GRCh38DH.20150718.YRI.high_coverage.cram"
FALLBACK_REF="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"

OUTPUT_BAM="${OUTDIR}/${SAMPLE}.${CHROM}_fragment.bam"

echo "Downloading ${SAMPLE} reads for region ${REGION}..."
echo "Output: ${OUTPUT_BAM}"
echo ""

# Function to try extracting from an indexed remote BAM
try_indexed_extract() {
    local bam_url="$1"
    local bai_url="$2"
    local output="$3"
    local extra_args="${4:-}"

    echo "Attempting indexed extraction from:"
    echo "  BAM: ${bam_url}"
    echo "  Region: ${REGION}"

    # Try with explicit index specification using -X
    # This tells samtools to use the second file as the index
    if samtools view -b ${extra_args} -X -o "${output}" "${bam_url}" "${bai_url}" "${REGION}" 2>/dev/null; then
        return 0
    fi

    # Fallback: samtools may auto-find .bai if it exists at the expected location
    if samtools view -b ${extra_args} -o "${output}" "${bam_url}" "${REGION}" 2>/dev/null; then
        return 0
    fi

    return 1
}

# Verify output file has content
verify_bam() {
    local bam="$1"
    if [[ ! -f "${bam}" ]]; then
        return 1
    fi
    local count
    count=$(samtools view -c "${bam}" 2>/dev/null || echo "0")
    [[ "${count}" -gt 0 ]]
}

# Try the HPRC aligned BAM first
echo "Trying HPRC aligned BAM..."
if try_indexed_extract "${HPRC_BAM}" "${HPRC_BAI}" "${OUTPUT_BAM}" && verify_bam "${OUTPUT_BAM}"; then
    echo "Successfully extracted region from HPRC BAM"
else
    echo "HPRC BAM not available or not indexed. Trying 1000 Genomes high-coverage CRAM..."

    # The 1000G high-coverage data uses CRAM format and needs the reference
    if try_indexed_extract "${FALLBACK_BAM}" "${FALLBACK_BAM}.crai" "${OUTPUT_BAM}" "-T ${FALLBACK_REF}" && verify_bam "${OUTPUT_BAM}"; then
        echo "Successfully extracted region from 1000G CRAM"
    else
        echo ""
        echo "ERROR: Could not extract region from remote indexed files."
        echo ""
        echo "Manual fallback: Download a pre-aligned, indexed BAM locally."
        echo "You can find NA19240 HiFi data at:"
        echo "  - https://www.internationalgenome.org/data-portal/sample/NA19240"
        echo "  - https://github.com/human-pangenomics/HPP_Year1_Assemblies"
        echo ""
        echo "Or generate synthetic test data with:"
        echo "  # The existing test BAM in this directory should work for basic testing"
        exit 1
    fi
fi

# Index the output BAM
echo ""
echo "Indexing output BAM..."
samtools index "${OUTPUT_BAM}"

echo ""
echo "Done. Extracted reads:"
samtools idxstats "${OUTPUT_BAM}"
echo ""
echo "Read count in region:"
samtools view -c "${OUTPUT_BAM}"
