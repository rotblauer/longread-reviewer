#!/usr/bin/env bash
# Download and extract NA19240 long-read data from public sources for testing.
#
# Source: 1000 Genomes Project / Human Pangenome Reference Consortium (HPRC)
#   - Sample: NA19240 (Yoruba, 1000 Genomes)
#   - Technology: PacBio HiFi
#   - Reference: GRCh38 (hg38)
#   - Region: chr17:10958130-10971414
#
# Public data URLs:
#   HPRC: https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html
#   1KG ONT: https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1KG_ONT_VIENNA/
#   HPRC HiFi: https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HiFi/
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
REGION="${CHROM}:10958130-10971414"
SAMPLE="NA19240"

# HPRC HiFi aligned BAM (hg38) — update URL as needed
# See: https://github.com/human-pangenomics/HPP_Year1_Assemblies
BAM_URL="https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HiFi/alignments/${SAMPLE}/${SAMPLE}.hifi.bam"

echo "Downloading ${SAMPLE} HiFi reads for region ${REGION}..."
echo "Source: ${BAM_URL}"
echo "Output: ${OUTDIR}"

# Use samtools to extract just the region of interest
samtools view -b -h "${BAM_URL}" "${REGION}" -o "${OUTDIR}/${SAMPLE}.${CHROM}_fragment.bam"
samtools index "${OUTDIR}/${SAMPLE}.${CHROM}_fragment.bam"

echo "Done. Extracted reads:"
samtools idxstats "${OUTDIR}/${SAMPLE}.${CHROM}_fragment.bam"
