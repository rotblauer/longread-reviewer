#!/usr/bin/env bash
# Download the reference genome region for testing.
#
# Source: UCSC Genome Browser / NCBI GRCh38
#   - Reference: GRCh38 (hg38)
#   - Region: chr17:10958130-11017414 (with buffer)
#
# This script extracts the reference sequence for the NA19240 test region
# plus a configurable buffer on each side to accommodate read overhangs.
#
# Requirements: samtools >= 1.10
#
# Usage:
#   ./download_reference.sh [output_dir] [buffer_size]
#
# Arguments:
#   output_dir  - Directory to save the reference FASTA (default: script directory)
#   buffer_size - Base pairs to add on each side (default: 5000)

set -euo pipefail

OUTDIR="${1:-$(dirname "$0")}"
BUFFER="${2:-5000}"

# Region of interest (same as download_from_source.sh)
CHROM="chr17"
REGION_START=10958130
REGION_END=11017414

# Calculate buffered region
BUFFERED_START=$((REGION_START - BUFFER))
BUFFERED_END=$((REGION_END + BUFFER))

# Ensure start doesn't go below 1
if [[ ${BUFFERED_START} -lt 1 ]]; then
    BUFFERED_START=1
fi

REGION="${CHROM}:${BUFFERED_START}-${BUFFERED_END}"
OUTPUT_FA="${OUTDIR}/${CHROM}_fragment.fa"

echo "Downloading reference genome region..."
echo "  Original region: ${CHROM}:${REGION_START}-${REGION_END}"
echo "  Buffer: ${BUFFER} bp on each side"
echo "  Buffered region: ${REGION}"
echo "  Output: ${OUTPUT_FA}"
echo ""

# GRCh38 reference sources (in order of preference)
# 1. UCSC - has indexed FASTA accessible via samtools
UCSC_REF="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"

# 2. NCBI/1000 Genomes - full analysis set
NCBI_REF="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"

# 3. Ensembl GRCh38
ENSEMBL_REF="http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

# Function to extract region from indexed remote FASTA
try_samtools_faidx() {
    local ref_url="$1"
    local region="$2"
    local output="$3"

    echo "Trying: ${ref_url}"

    # samtools faidx can work with remote indexed FASTA files
    if samtools faidx "${ref_url}" "${region}" > "${output}" 2>/dev/null; then
        # Check if we got actual sequence (not just header or error)
        if [[ -s "${output}" ]] && grep -q "^[ACGTNacgtn]" "${output}"; then
            return 0
        fi
    fi

    return 1
}

# Function to use UCSC DAS server (always works, no index needed)
try_ucsc_das() {
    local chrom="$1"
    local start="$2"
    local end="$3"
    local output="$4"

    # UCSC DAS server - reliable but limited to ~10MB chunks
    local das_url="https://genome.ucsc.edu/cgi-bin/das/hg38/dna?segment=${chrom}:${start},${end}"

    echo "Trying UCSC DAS server..."

    # DAS returns XML, we need to extract the sequence
    local xml_response
    xml_response=$(curl -sL "${das_url}" 2>/dev/null) || return 1

    # Extract sequence from XML (between <DNA> tags)
    local sequence
    sequence=$(echo "${xml_response}" | sed -n '/<DNA/,/<\/DNA>/p' | grep -v '<' | tr -d ' \n\r')

    if [[ -n "${sequence}" ]]; then
        # Write as FASTA
        echo ">${chrom}:${start}-${end}" > "${output}"
        # Format to 60 characters per line
        echo "${sequence}" | fold -w 60 >> "${output}"
        return 0
    fi

    return 1
}

# Function to use Ensembl REST API
try_ensembl_rest() {
    local chrom="$1"
    local start="$2"
    local end="$3"
    local output="$4"

    # Ensembl uses chromosome names without 'chr' prefix
    local ensembl_chrom="${chrom#chr}"

    local api_url="https://rest.ensembl.org/sequence/region/human/${ensembl_chrom}:${start}..${end}?content-type=text/x-fasta"

    echo "Trying Ensembl REST API..."

    local response
    response=$(curl -sL "${api_url}" 2>/dev/null) || return 1

    if [[ "${response}" == ">"* ]]; then
        # Rename the header to match our naming convention
        echo ">${chrom}:${start}-${end}" > "${output}"
        echo "${response}" | tail -n +2 >> "${output}"
        return 0
    fi

    return 1
}

# Try methods in order of preference
success=false

# Method 1: Try UCSC DAS (most reliable, no index needed)
if try_ucsc_das "${CHROM}" "${BUFFERED_START}" "${BUFFERED_END}" "${OUTPUT_FA}"; then
    echo "Successfully downloaded from UCSC DAS server"
    success=true
fi

# Method 2: Try Ensembl REST API
if [[ "${success}" == "false" ]]; then
    if try_ensembl_rest "${CHROM}" "${BUFFERED_START}" "${BUFFERED_END}" "${OUTPUT_FA}"; then
        echo "Successfully downloaded from Ensembl REST API"
        success=true
    fi
fi

# Method 3: Try samtools faidx with remote indexed FASTA
if [[ "${success}" == "false" ]]; then
    if try_samtools_faidx "${NCBI_REF}" "${REGION}" "${OUTPUT_FA}"; then
        echo "Successfully extracted from NCBI reference"
        success=true
    fi
fi

if [[ "${success}" == "false" ]]; then
    echo ""
    echo "ERROR: Could not download reference sequence."
    echo ""
    echo "Manual fallback: Download GRCh38 reference locally and extract region:"
    echo "  samtools faidx hg38.fa ${REGION} > ${OUTPUT_FA}"
    exit 1
fi

# Index the output FASTA
echo ""
echo "Indexing output FASTA..."
samtools faidx "${OUTPUT_FA}"

# Summary
echo ""
echo "Done. Reference sequence:"
echo "  File: ${OUTPUT_FA}"
echo "  Index: ${OUTPUT_FA}.fai"
echo ""
cat "${OUTPUT_FA}.fai"
echo ""
echo "Sequence length: $(grep -v '^>' "${OUTPUT_FA}" | tr -d '\n' | wc -c | tr -d ' ') bp"

