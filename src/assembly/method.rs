use std::collections::HashMap;

use anyhow::Result;

use crate::alignment::AlignedRead;

/// Result of an assembly over a region.
#[derive(Debug, Clone)]
pub struct AssemblyResult {
    /// The assembled consensus sequence.
    pub sequence: Vec<u8>,
    /// Per-base depth (number of reads supporting each base).
    pub depth: Vec<u32>,
    /// Per-base confidence score (0.0 to 1.0).
    pub confidence: Vec<f64>,
    /// Name/identifier for this assembly method.
    pub method_name: String,
}

/// Trait for pluggable assembly methods.
///
/// Users can implement this trait to experiment with different assembly strategies.
/// Each method takes a set of aligned reads and a reference sequence, and produces
/// an assembly result.
pub trait AssemblyMethod: Send + Sync {
    /// Name of this assembly method.
    fn name(&self) -> &str;

    /// Generate an assembly from the given reads against the reference.
    fn assemble(&self, reads: &[AlignedRead], reference: &[u8]) -> Result<AssemblyResult>;
}

/// Simple majority-vote consensus assembly.
///
/// At each reference position, the most common base among overlapping reads is chosen.
pub struct ConsensusAssembly;

impl AssemblyMethod for ConsensusAssembly {
    fn name(&self) -> &str {
        "consensus"
    }

    fn assemble(&self, reads: &[AlignedRead], reference: &[u8]) -> Result<AssemblyResult> {
        let ref_len = reference.len();
        if ref_len == 0 {
            return Ok(AssemblyResult {
                sequence: Vec::new(),
                depth: Vec::new(),
                confidence: Vec::new(),
                method_name: self.name().to_string(),
            });
        }

        // For each position, count base occurrences using quality-weighted voting.
        let mut base_counts: Vec<HashMap<u8, f64>> = vec![HashMap::new(); ref_len];
        let mut depth = vec![0u32; ref_len];

        // We assume reads are aligned to the same coordinate system as the reference slice.
        // The reference slice starts at position 0 internally.
        let ref_start = reads
            .iter()
            .filter(|r| !r.cigar.is_empty())
            .map(|r| r.start)
            .min()
            .unwrap_or(0);

        for read in reads {
            let aligned_bases = read.reference_aligned_bases();
            for (ref_pos, base) in aligned_bases {
                let idx = (ref_pos - ref_start) as usize;
                if idx < ref_len {
                    let quality_weight = read
                        .qualities
                        .first()
                        .copied()
                        .unwrap_or(30) as f64
                        / 60.0;
                    *base_counts[idx]
                        .entry(base.to_ascii_uppercase())
                        .or_insert(0.0) += 1.0 + quality_weight;
                    depth[idx] += 1;
                }
            }
        }

        let mut sequence = Vec::with_capacity(ref_len);
        let mut confidence = Vec::with_capacity(ref_len);

        for (i, counts) in base_counts.iter().enumerate() {
            if counts.is_empty() {
                // No coverage: fall back to reference base.
                sequence.push(reference[i]);
                confidence.push(0.0);
            } else {
                let total: f64 = counts.values().sum();
                let (best_base, best_count) = counts
                    .iter()
                    .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
                    .unwrap();
                sequence.push(*best_base);
                confidence.push(*best_count / total);
            }
        }

        Ok(AssemblyResult {
            sequence,
            depth,
            confidence,
            method_name: self.name().to_string(),
        })
    }
}

/// Windowed consensus assembly that processes the region in overlapping windows.
///
/// This allows experimenting with different window sizes and overlap strategies
/// to handle complex structural variants.
pub struct WindowConsensusAssembly {
    /// Size of each assembly window in bases.
    pub window_size: usize,
    /// Number of bases to overlap between windows.
    pub overlap: usize,
}

impl WindowConsensusAssembly {
    pub fn new(window_size: usize, overlap: usize) -> Self {
        Self {
            window_size,
            overlap,
        }
    }
}

impl AssemblyMethod for WindowConsensusAssembly {
    fn name(&self) -> &str {
        "window_consensus"
    }

    fn assemble(&self, reads: &[AlignedRead], reference: &[u8]) -> Result<AssemblyResult> {
        let ref_len = reference.len();
        if ref_len == 0 {
            return Ok(AssemblyResult {
                sequence: Vec::new(),
                depth: Vec::new(),
                confidence: Vec::new(),
                method_name: self.name().to_string(),
            });
        }

        let inner = ConsensusAssembly;
        let step = if self.window_size > self.overlap {
            self.window_size - self.overlap
        } else {
            1
        };

        // Accumulate weighted votes from each window.
        let mut base_votes: Vec<HashMap<u8, f64>> = vec![HashMap::new(); ref_len];
        let mut total_depth = vec![0u32; ref_len];

        let mut window_start = 0;
        while window_start < ref_len {
            let window_end = std::cmp::min(window_start + self.window_size, ref_len);
            let window_ref = &reference[window_start..window_end];

            // Filter reads that overlap this window.
            let ref_start = reads
                .iter()
                .filter(|r| !r.cigar.is_empty())
                .map(|r| r.start)
                .min()
                .unwrap_or(0);

            let window_reads: Vec<_> = reads
                .iter()
                .filter(|r| {
                    let r_start = (r.start - ref_start) as usize;
                    let r_end = (r.end - ref_start) as usize;
                    r_start < window_end && r_end >= window_start
                })
                .cloned()
                .collect();

            let result = inner.assemble(&window_reads, window_ref)?;

            for (i, (&base, &conf)) in result
                .sequence
                .iter()
                .zip(result.confidence.iter())
                .enumerate()
            {
                let global_idx = window_start + i;
                if global_idx < ref_len {
                    *base_votes[global_idx].entry(base).or_insert(0.0) += conf;
                    total_depth[global_idx] = total_depth[global_idx].max(result.depth[i]);
                }
            }

            window_start += step;
        }

        let mut sequence = Vec::with_capacity(ref_len);
        let mut confidence = Vec::with_capacity(ref_len);

        for (i, votes) in base_votes.iter().enumerate() {
            if votes.is_empty() {
                sequence.push(reference[i]);
                confidence.push(0.0);
            } else {
                let total: f64 = votes.values().sum();
                let (best_base, best_score) = votes
                    .iter()
                    .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
                    .unwrap();
                sequence.push(*best_base);
                confidence.push(if total > 0.0 { *best_score / total } else { 0.0 });
            }
        }

        Ok(AssemblyResult {
            sequence,
            depth: total_depth,
            confidence,
            method_name: self.name().to_string(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::CigarOp;

    fn make_read(name: &str, start: u64, seq: &[u8]) -> AlignedRead {
        let len = seq.len() as u64;
        AlignedRead {
            name: name.to_string(),
            chrom: "chr1".to_string(),
            start,
            mapq: 60,
            cigar: vec![CigarOp::Match(seq.len() as u32)],
            sequence: seq.to_vec(),
            qualities: vec![30; seq.len()],
            is_reverse: false,
            haplotype_tag: None,
            end: start + len - 1,
        }
    }

    #[test]
    fn test_consensus_assembly_simple() {
        let reference = b"ACGTACGT";
        let reads = vec![
            make_read("r1", 1, b"ACGTACGT"),
            make_read("r2", 1, b"ACGTACGT"),
            make_read("r3", 1, b"ACGTACGT"),
        ];

        let method = ConsensusAssembly;
        let result = method.assemble(&reads, reference).unwrap();

        assert_eq!(result.sequence, b"ACGTACGT");
        assert_eq!(result.depth, vec![3, 3, 3, 3, 3, 3, 3, 3]);
        assert!(result.confidence.iter().all(|&c| c == 1.0));
    }

    #[test]
    fn test_consensus_assembly_with_variant() {
        let reference = b"ACGTACGT";
        // Position 4 (index 3): 2 reads say T, 1 says G -> consensus should be T
        let reads = vec![
            make_read("r1", 1, b"ACGTACGT"),
            make_read("r2", 1, b"ACGTACGT"),
            make_read("r3", 1, b"ACGGACGT"),
        ];

        let method = ConsensusAssembly;
        let result = method.assemble(&reads, reference).unwrap();

        assert_eq!(result.sequence[3], b'T');
    }

    #[test]
    fn test_consensus_empty_reads() {
        let reference = b"ACGT";
        let reads: Vec<AlignedRead> = vec![];

        let method = ConsensusAssembly;
        let result = method.assemble(&reads, reference).unwrap();

        // Should fall back to reference
        assert_eq!(result.sequence, b"ACGT");
        assert!(result.confidence.iter().all(|&c| c == 0.0));
    }

    #[test]
    fn test_window_consensus_assembly() {
        let reference = b"ACGTACGTACGT";
        let reads = vec![
            make_read("r1", 1, b"ACGTACGTACGT"),
            make_read("r2", 1, b"ACGTACGTACGT"),
        ];

        let method = WindowConsensusAssembly::new(6, 2);
        let result = method.assemble(&reads, reference).unwrap();

        assert_eq!(result.sequence.len(), reference.len());
    }

    #[test]
    fn test_empty_reference() {
        let reads: Vec<AlignedRead> = vec![];

        let method = ConsensusAssembly;
        let result = method.assemble(&reads, b"").unwrap();
        assert!(result.sequence.is_empty());
    }
}
