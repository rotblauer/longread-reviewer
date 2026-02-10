use std::path::Path;

use anyhow::{Context, Result};
use noodles::sam::alignment::record::data::field::value::Value;

use crate::region::Region;

/// Simplified CIGAR operation for internal use.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarOp {
    /// Alignment match (can be sequence match or mismatch).
    Match(u32),
    /// Insertion to the reference.
    Insertion(u32),
    /// Deletion from the reference.
    Deletion(u32),
    /// Soft clip (bases present in read but not aligned).
    SoftClip(u32),
    /// Hard clip (bases not present in read).
    HardClip(u32),
}

impl CigarOp {
    /// Number of bases this operation consumes on the reference.
    pub fn ref_len(&self) -> u32 {
        match self {
            CigarOp::Match(n) | CigarOp::Deletion(n) => *n,
            CigarOp::Insertion(_) | CigarOp::SoftClip(_) | CigarOp::HardClip(_) => 0,
        }
    }

    /// Number of bases this operation consumes on the read.
    pub fn read_len(&self) -> u32 {
        match self {
            CigarOp::Match(n) | CigarOp::Insertion(n) | CigarOp::SoftClip(n) => *n,
            CigarOp::Deletion(_) | CigarOp::HardClip(_) => 0,
        }
    }
}

/// A single aligned read extracted from a BAM/CRAM file.
#[derive(Debug, Clone)]
pub struct AlignedRead {
    /// Read name / query name.
    pub name: String,
    /// Chromosome / reference name the read is aligned to.
    pub chrom: String,
    /// 1-based start position on the reference.
    pub start: u64,
    /// Mapping quality.
    pub mapq: u8,
    /// CIGAR operations describing the alignment.
    pub cigar: Vec<CigarOp>,
    /// Read sequence (ASCII bases).
    pub sequence: Vec<u8>,
    /// Per-base quality scores (Phred-scaled).
    pub qualities: Vec<u8>,
    /// Whether this read is on the reverse strand.
    pub is_reverse: bool,
    /// HP tag value if present (for haplotagged reads).
    pub haplotype_tag: Option<u8>,
    /// Alignment end on reference (1-based, inclusive).
    pub end: u64,
}

impl AlignedRead {
    /// Get the aligned portion of the sequence (excluding soft/hard clips).
    pub fn aligned_sequence(&self) -> Vec<u8> {
        let mut result = Vec::new();
        let mut read_pos = 0usize;
        for op in &self.cigar {
            match op {
                CigarOp::Match(n) => {
                    let n = *n as usize;
                    if read_pos + n <= self.sequence.len() {
                        result.extend_from_slice(&self.sequence[read_pos..read_pos + n]);
                    }
                    read_pos += n;
                }
                CigarOp::Insertion(n) => {
                    let n = *n as usize;
                    if read_pos + n <= self.sequence.len() {
                        result.extend_from_slice(&self.sequence[read_pos..read_pos + n]);
                    }
                    read_pos += n;
                }
                CigarOp::SoftClip(n) => {
                    read_pos += *n as usize;
                }
                CigarOp::Deletion(_) => {}
                CigarOp::HardClip(_) => {}
            }
        }
        result
    }

    /// Compute the aligned bases at each reference position within the read's span.
    /// Returns a Vec of (ref_pos_1based, read_base) tuples for match/mismatch operations.
    pub fn reference_aligned_bases(&self) -> Vec<(u64, u8)> {
        let mut result = Vec::new();
        let mut ref_pos = self.start;
        let mut read_pos = 0usize;

        for op in &self.cigar {
            match op {
                CigarOp::Match(n) => {
                    for i in 0..*n as usize {
                        if read_pos + i < self.sequence.len() {
                            result.push((ref_pos + i as u64, self.sequence[read_pos + i]));
                        }
                    }
                    ref_pos += *n as u64;
                    read_pos += *n as usize;
                }
                CigarOp::Insertion(n) => {
                    read_pos += *n as usize;
                }
                CigarOp::Deletion(n) => {
                    ref_pos += *n as u64;
                }
                CigarOp::SoftClip(n) => {
                    read_pos += *n as usize;
                }
                CigarOp::HardClip(_) => {}
            }
        }
        result
    }
}

/// Reader for BAM/CRAM alignment files.
pub struct AlignmentReader;

/// The HP (haplotype) tag used in haplotagged BAM files.
const HP_TAG: [u8; 2] = [b'H', b'P'];

impl AlignmentReader {
    /// Read alignments from a BAM file for a given region.
    ///
    /// The BAM file must be sorted and indexed (.bai).
    pub fn read_bam(path: &Path, region: &Region) -> Result<Vec<AlignedRead>> {
        let mut reader = noodles::bam::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .with_context(|| format!("failed to open BAM file: {}", path.display()))?;

        let header = reader.read_header().context("failed to read BAM header")?;

        let region_str = format!("{}:{}-{}", region.chrom, region.start, region.end);
        let query_region: noodles::core::Region = region_str
            .parse()
            .with_context(|| format!("failed to parse region: {region_str}"))?;

        let mut reads = Vec::new();
        let query = reader
            .query(&header, &query_region)
            .context("failed to query BAM region")?;

        for result in query.records() {
            let record = result.context("failed to read BAM record")?;

            if record.flags().is_unmapped() {
                continue;
            }

            if let Some(read) = Self::convert_bam_record(&record, &header)? {
                reads.push(read);
            }
        }

        Ok(reads)
    }

    /// Convert a noodles BAM record to our AlignedRead type.
    fn convert_bam_record(
        record: &noodles::bam::Record,
        header: &noodles::sam::Header,
    ) -> Result<Option<AlignedRead>> {
        let name = record
            .name()
            .map(|n| String::from_utf8_lossy(n).into_owned())
            .unwrap_or_else(|| "unknown".to_string());

        let flags = record.flags();
        let is_reverse = flags.is_reverse_complemented();

        let ref_seq_id = match record.reference_sequence_id() {
            Some(Ok(id)) => id,
            Some(Err(e)) => return Err(e).context("failed to read reference sequence ID"),
            None => return Ok(None),
        };

        let chrom = header
            .reference_sequences()
            .get_index(ref_seq_id)
            .map(|(name, _)| name.to_string())
            .unwrap_or_else(|| "unknown".to_string());

        let start = match record.alignment_start() {
            Some(Ok(p)) => p.get() as u64,
            Some(Err(e)) => return Err(e).context("failed to read alignment start"),
            None => 0,
        };

        let mapq = record.mapping_quality().map(|q| q.get()).unwrap_or(0);

        let cigar = Self::convert_cigar(record)?;
        let ref_consumed: u64 = cigar.iter().map(|op| op.ref_len() as u64).sum();
        let end = start + ref_consumed.saturating_sub(1);

        let sequence: Vec<u8> = (0..record.sequence().len())
            .filter_map(|i| record.sequence().get(i))
            .collect();

        let qualities: Vec<u8> = record
            .quality_scores()
            .as_ref().to_vec();

        // Try to read HP (haplotype) tag
        let haplotype_tag = record
            .data()
            .get(&HP_TAG)
            .and_then(|v| match v {
                Ok(Value::UInt8(n)) => Some(n),
                _ => None,
            });

        Ok(Some(AlignedRead {
            name,
            chrom,
            start,
            mapq,
            cigar,
            sequence,
            qualities,
            is_reverse,
            haplotype_tag,
            end,
        }))
    }

    /// Convert noodles CIGAR to our simplified CIGAR representation.
    fn convert_cigar(record: &noodles::bam::Record) -> Result<Vec<CigarOp>> {
        use noodles::sam::alignment::record::cigar::op::Kind;

        let mut ops = Vec::new();
        for result in record.cigar().iter() {
            let op = result.context("failed to read CIGAR operation")?;
            let len = op.len() as u32;
            let converted = match op.kind() {
                Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => CigarOp::Match(len),
                Kind::Insertion => CigarOp::Insertion(len),
                Kind::Deletion => CigarOp::Deletion(len),
                Kind::SoftClip => CigarOp::SoftClip(len),
                Kind::HardClip => CigarOp::HardClip(len),
                _ => continue, // skip padding, etc.
            };
            ops.push(converted);
        }
        Ok(ops)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cigar_op_lengths() {
        assert_eq!(CigarOp::Match(10).ref_len(), 10);
        assert_eq!(CigarOp::Match(10).read_len(), 10);
        assert_eq!(CigarOp::Insertion(5).ref_len(), 0);
        assert_eq!(CigarOp::Insertion(5).read_len(), 5);
        assert_eq!(CigarOp::Deletion(3).ref_len(), 3);
        assert_eq!(CigarOp::Deletion(3).read_len(), 0);
        assert_eq!(CigarOp::SoftClip(4).ref_len(), 0);
        assert_eq!(CigarOp::SoftClip(4).read_len(), 4);
    }

    #[test]
    fn test_aligned_read_basics() {
        let read = AlignedRead {
            name: "read1".to_string(),
            chrom: "chr1".to_string(),
            start: 100,
            mapq: 60,
            cigar: vec![CigarOp::Match(10)],
            sequence: b"ACGTACGTAC".to_vec(),
            qualities: vec![30; 10],
            is_reverse: false,
            haplotype_tag: None,
            end: 109,
        };

        assert_eq!(read.aligned_sequence(), b"ACGTACGTAC");
    }

    #[test]
    fn test_aligned_read_with_softclip() {
        let read = AlignedRead {
            name: "read2".to_string(),
            chrom: "chr1".to_string(),
            start: 100,
            mapq: 60,
            cigar: vec![
                CigarOp::SoftClip(2),
                CigarOp::Match(6),
                CigarOp::SoftClip(2),
            ],
            sequence: b"NNACGTACNN".to_vec(),
            qualities: vec![30; 10],
            is_reverse: false,
            haplotype_tag: None,
            end: 105,
        };

        // Aligned portion should exclude soft-clipped bases
        assert_eq!(read.aligned_sequence(), b"ACGTAC");
    }

    #[test]
    fn test_reference_aligned_bases() {
        let read = AlignedRead {
            name: "read3".to_string(),
            chrom: "chr1".to_string(),
            start: 100,
            mapq: 60,
            cigar: vec![CigarOp::Match(4), CigarOp::Deletion(2), CigarOp::Match(4)],
            sequence: b"ACGTACGT".to_vec(),
            qualities: vec![30; 8],
            is_reverse: false,
            haplotype_tag: None,
            end: 109,
        };

        let bases = read.reference_aligned_bases();
        assert_eq!(bases.len(), 8);
        // First 4 bases at positions 100-103
        assert_eq!(bases[0], (100, b'A'));
        assert_eq!(bases[3], (103, b'T'));
        // After deletion (skip 104,105), next bases at 106-109
        assert_eq!(bases[4], (106, b'A'));
        assert_eq!(bases[7], (109, b'T'));
    }

    #[test]
    fn test_cigar_op_hardclip() {
        assert_eq!(CigarOp::HardClip(10).ref_len(), 0);
        assert_eq!(CigarOp::HardClip(10).read_len(), 0);
    }

    #[test]
    fn test_aligned_read_with_insertion() {
        let read = AlignedRead {
            name: "ins_read".to_string(),
            chrom: "chr1".to_string(),
            start: 100,
            mapq: 60,
            cigar: vec![
                CigarOp::Match(4),
                CigarOp::Insertion(2),
                CigarOp::Match(4),
            ],
            // seq: ACGT + NN (insertion) + ACGT = 10 bases
            sequence: b"ACGTNNACGT".to_vec(),
            qualities: vec![30; 10],
            is_reverse: false,
            haplotype_tag: None,
            end: 107,
        };

        let aligned = read.aligned_sequence();
        // aligned_sequence includes insertion bases
        assert_eq!(aligned, b"ACGTNNACGT");

        let ref_bases = read.reference_aligned_bases();
        // Reference-aligned bases skip the insertion
        assert_eq!(ref_bases.len(), 8);
        assert_eq!(ref_bases[0], (100, b'A'));
        assert_eq!(ref_bases[3], (103, b'T'));
        // After insertion, ref pos continues from 104
        assert_eq!(ref_bases[4], (104, b'A'));
        assert_eq!(ref_bases[7], (107, b'T'));
    }

    #[test]
    fn test_aligned_read_with_deletion() {
        let read = AlignedRead {
            name: "del_read".to_string(),
            chrom: "chr1".to_string(),
            start: 100,
            mapq: 60,
            cigar: vec![
                CigarOp::Match(3),
                CigarOp::Deletion(5),
                CigarOp::Match(3),
            ],
            sequence: b"ACGACG".to_vec(),
            qualities: vec![30; 6],
            is_reverse: false,
            haplotype_tag: None,
            end: 110,
        };

        let aligned = read.aligned_sequence();
        assert_eq!(aligned, b"ACGACG");

        let ref_bases = read.reference_aligned_bases();
        assert_eq!(ref_bases.len(), 6);
        // First 3 at 100-102
        assert_eq!(ref_bases[0], (100, b'A'));
        assert_eq!(ref_bases[2], (102, b'G'));
        // Skip 5 positions (103-107), next at 108
        assert_eq!(ref_bases[3], (108, b'A'));
    }

    #[test]
    fn test_aligned_read_complex_cigar() {
        // S2 M3 I1 M2 D2 M3 S2
        let read = AlignedRead {
            name: "complex".to_string(),
            chrom: "chr1".to_string(),
            start: 100,
            mapq: 60,
            cigar: vec![
                CigarOp::SoftClip(2),
                CigarOp::Match(3),
                CigarOp::Insertion(1),
                CigarOp::Match(2),
                CigarOp::Deletion(2),
                CigarOp::Match(3),
                CigarOp::SoftClip(2),
            ],
            // 2 SC + 3 M + 1 I + 2 M + 3 M + 2 SC = 13
            sequence: b"NNACGNACGTTNN".to_vec(),
            qualities: vec![30; 13],
            is_reverse: false,
            haplotype_tag: None,
            end: 111,
        };

        let aligned = read.aligned_sequence();
        // Soft clips excluded, rest included
        assert_eq!(aligned, b"ACGNACGTT");

        let ref_bases = read.reference_aligned_bases();
        // M3 gives 3 bases at 100-102, M2 gives 2 at 103-104, skip D2 (105-106), M3 at 107-109
        assert_eq!(ref_bases.len(), 8);
    }

    #[test]
    fn test_aligned_read_only_softclips() {
        let read = AlignedRead {
            name: "all_sc".to_string(),
            chrom: "chr1".to_string(),
            start: 100,
            mapq: 0,
            cigar: vec![CigarOp::SoftClip(10)],
            sequence: b"AAAAAAAAAA".to_vec(),
            qualities: vec![30; 10],
            is_reverse: false,
            haplotype_tag: None,
            end: 99,
        };

        let aligned = read.aligned_sequence();
        assert!(aligned.is_empty());

        let ref_bases = read.reference_aligned_bases();
        assert!(ref_bases.is_empty());
    }

    #[test]
    fn test_aligned_read_empty_sequence() {
        let read = AlignedRead {
            name: "empty".to_string(),
            chrom: "chr1".to_string(),
            start: 100,
            mapq: 0,
            cigar: vec![],
            sequence: vec![],
            qualities: vec![],
            is_reverse: false,
            haplotype_tag: None,
            end: 100,
        };

        assert!(read.aligned_sequence().is_empty());
        assert!(read.reference_aligned_bases().is_empty());
    }

    #[test]
    fn test_aligned_read_reverse_strand() {
        let read = AlignedRead {
            name: "rev".to_string(),
            chrom: "chr1".to_string(),
            start: 100,
            mapq: 60,
            cigar: vec![CigarOp::Match(4)],
            sequence: b"TGCA".to_vec(),
            qualities: vec![30; 4],
            is_reverse: true,
            haplotype_tag: None,
            end: 103,
        };

        assert!(read.is_reverse);
        assert_eq!(read.aligned_sequence(), b"TGCA");
    }

    #[test]
    fn test_aligned_read_with_haplotype_tag() {
        let read = AlignedRead {
            name: "hp".to_string(),
            chrom: "chr1".to_string(),
            start: 100,
            mapq: 60,
            cigar: vec![CigarOp::Match(4)],
            sequence: b"ACGT".to_vec(),
            qualities: vec![30; 4],
            is_reverse: false,
            haplotype_tag: Some(1),
            end: 103,
        };

        assert_eq!(read.haplotype_tag, Some(1));
    }

    #[test]
    fn test_hardclip_does_not_affect_alignment() {
        let read = AlignedRead {
            name: "hc".to_string(),
            chrom: "chr1".to_string(),
            start: 100,
            mapq: 60,
            cigar: vec![
                CigarOp::HardClip(5),
                CigarOp::Match(4),
                CigarOp::HardClip(5),
            ],
            sequence: b"ACGT".to_vec(),
            qualities: vec![30; 4],
            is_reverse: false,
            haplotype_tag: None,
            end: 103,
        };

        assert_eq!(read.aligned_sequence(), b"ACGT");
        let ref_bases = read.reference_aligned_bases();
        assert_eq!(ref_bases.len(), 4);
    }

    #[test]
    fn test_cigar_ref_consumed_total() {
        let cigar = vec![
            CigarOp::SoftClip(5),
            CigarOp::Match(10),
            CigarOp::Insertion(3),
            CigarOp::Match(5),
            CigarOp::Deletion(2),
            CigarOp::Match(8),
            CigarOp::SoftClip(4),
        ];
        let ref_consumed: u32 = cigar.iter().map(|op| op.ref_len()).sum();
        assert_eq!(ref_consumed, 10 + 5 + 2 + 8); // 25

        let read_consumed: u32 = cigar.iter().map(|op| op.read_len()).sum();
        assert_eq!(read_consumed, 5 + 10 + 3 + 5 + 8 + 4); // 35
    }
}
