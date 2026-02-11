use std::collections::HashMap;

use anyhow::Result;

use crate::alignment::AlignedRead;
use crate::haplotype::{HaplotypeAssigner, HaplotypeLabel};

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
    ///
    /// # Arguments
    /// * `reads` - Aligned reads overlapping the region
    /// * `reference` - Reference sequence for the region
    /// * `ref_start_pos` - The genomic position where the reference slice starts (1-based)
    fn assemble(&self, reads: &[AlignedRead], reference: &[u8], ref_start_pos: u64) -> Result<AssemblyResult>;
}

/// Simple majority-vote consensus assembly.
///
/// At each reference position, the most common base among overlapping reads is chosen.
pub struct ConsensusAssembly;

impl AssemblyMethod for ConsensusAssembly {
    fn name(&self) -> &str {
        "consensus"
    }

    fn assemble(&self, reads: &[AlignedRead], reference: &[u8], ref_start_pos: u64) -> Result<AssemblyResult> {
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

        // Use the provided ref_start_pos to map genomic coordinates to array indices
        let ref_start = ref_start_pos;

        for read in reads {
            let aligned_bases = read.reference_aligned_bases();
            for (ref_pos, base) in aligned_bases {
                // Skip bases outside our reference region
                if ref_pos < ref_start {
                    continue;
                }
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

    fn assemble(&self, reads: &[AlignedRead], reference: &[u8], ref_start_pos: u64) -> Result<AssemblyResult> {
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

            // Use the provided ref_start_pos for coordinate mapping
            let window_genomic_start = ref_start_pos + window_start as u64;
            let window_genomic_end = ref_start_pos + window_end as u64;

            let window_reads: Vec<_> = reads
                .iter()
                .filter(|r| {
                    // Filter reads that overlap this window in genomic coordinates
                    r.start < window_genomic_end && r.end > window_genomic_start
                })
                .cloned()
                .collect();

            let result = inner.assemble(&window_reads, window_ref, window_genomic_start)?;

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

/// Result of per-haplotype assembly, containing separate assemblies for each
/// detected haplotype and a divergence track highlighting positions where the
/// two haplotypes differ.
#[derive(Debug, Clone)]
pub struct HaplotypeAssemblyResult {
    /// Assembly from Haplotype 1 reads only.
    pub hap1_assembly: AssemblyResult,
    /// Assembly from Haplotype 2 reads only.
    pub hap2_assembly: AssemblyResult,
    /// Number of reads assigned to Haplotype 1.
    pub hap1_read_count: usize,
    /// Number of reads assigned to Haplotype 2.
    pub hap2_read_count: usize,
    /// Per-base divergence between the two haplotype assemblies (0.0 = same, 1.0 = different).
    pub divergence: Vec<f64>,
    /// Per-base SV likelihood score (0.0-1.0), derived from divergence, depth, and confidence.
    pub sv_likelihood: Vec<f64>,
}

/// Produces per-haplotype assemblies by splitting reads into haplotype groups
/// and assembling each independently using a given inner assembly method.
///
/// This is the recommended approach for reviewing structural variants: if an SV
/// is real, each haplotype assembly should show a clean, consistent signal.
/// Divergence between haplotypes pinpoints heterozygous events.
pub struct HaplotypeAwareAssembly {
    inner: ConsensusAssembly,
}

impl HaplotypeAwareAssembly {
    pub fn new() -> Self {
        Self {
            inner: ConsensusAssembly,
        }
    }

    /// Split reads by haplotype and assemble each group independently.
    pub fn assemble_by_haplotype(
        &self,
        reads: &[AlignedRead],
        reference: &[u8],
        ref_start_pos: u64,
    ) -> Result<HaplotypeAssemblyResult> {
        let assigner = HaplotypeAssigner::new();
        let assignments = assigner.assign(reads, reference);

        let assignment_map: HashMap<String, HaplotypeLabel> = assignments
            .iter()
            .map(|a| (a.read_name.clone(), a.haplotype))
            .collect();

        let hap1_reads: Vec<AlignedRead> = reads
            .iter()
            .filter(|r| assignment_map.get(&r.name) == Some(&HaplotypeLabel::Hap1))
            .cloned()
            .collect();

        let hap2_reads: Vec<AlignedRead> = reads
            .iter()
            .filter(|r| assignment_map.get(&r.name) == Some(&HaplotypeLabel::Hap2))
            .cloned()
            .collect();

        let hap1_count = hap1_reads.len();
        let hap2_count = hap2_reads.len();

        let mut hap1_assembly = self.inner.assemble(&hap1_reads, reference, ref_start_pos)?;
        hap1_assembly.method_name = "haplotype_1".to_string();

        let mut hap2_assembly = self.inner.assemble(&hap2_reads, reference, ref_start_pos)?;
        hap2_assembly.method_name = "haplotype_2".to_string();

        // Compute per-base divergence and SV likelihood
        let len = reference.len();
        let mut divergence = vec![0.0f64; len];
        let mut sv_likelihood = vec![0.0f64; len];

        for i in 0..len {
            let h1_base = if i < hap1_assembly.sequence.len() {
                hap1_assembly.sequence[i]
            } else {
                reference[i]
            };
            let h2_base = if i < hap2_assembly.sequence.len() {
                hap2_assembly.sequence[i]
            } else {
                reference[i]
            };

            // Divergence: do the two haplotypes disagree?
            let is_divergent = !h1_base.eq_ignore_ascii_case(&h2_base);
            divergence[i] = if is_divergent { 1.0 } else { 0.0 };

            // SV likelihood: combine divergence, depth, and confidence signals
            let h1_conf = if i < hap1_assembly.confidence.len() {
                hap1_assembly.confidence[i]
            } else {
                0.0
            };
            let h2_conf = if i < hap2_assembly.confidence.len() {
                hap2_assembly.confidence[i]
            } else {
                0.0
            };
            let h1_depth = if i < hap1_assembly.depth.len() {
                hap1_assembly.depth[i]
            } else {
                0
            };
            let h2_depth = if i < hap2_assembly.depth.len() {
                hap2_assembly.depth[i]
            } else {
                0
            };

            // Both haplotypes need adequate depth for a reliable call
            let depth_factor = ((h1_depth.min(h2_depth) as f64) / 5.0).min(1.0);
            // Both haplotypes should have good confidence in their respective calls
            let conf_factor = (h1_conf * h2_conf).sqrt();

            sv_likelihood[i] = divergence[i] * depth_factor * conf_factor;
        }

        Ok(HaplotypeAssemblyResult {
            hap1_assembly,
            hap2_assembly,
            hap1_read_count: hap1_count,
            hap2_read_count: hap2_count,
            divergence,
            sv_likelihood,
        })
    }
}

impl Default for HaplotypeAwareAssembly {
    fn default() -> Self {
        Self::new()
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
        let result = method.assemble(&reads, reference, 1).unwrap();

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
        let result = method.assemble(&reads, reference, 1).unwrap();

        assert_eq!(result.sequence[3], b'T');
    }

    #[test]
    fn test_consensus_empty_reads() {
        let reference = b"ACGT";
        let reads: Vec<AlignedRead> = vec![];

        let method = ConsensusAssembly;
        let result = method.assemble(&reads, reference, 1).unwrap();

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
        let result = method.assemble(&reads, reference, 1).unwrap();

        assert_eq!(result.sequence.len(), reference.len());
    }

    #[test]
    fn test_empty_reference() {
        let reads: Vec<AlignedRead> = vec![];

        let method = ConsensusAssembly;
        let result = method.assemble(&reads, b"", 1).unwrap();
        assert!(result.sequence.is_empty());
    }

    fn make_hap_read(name: &str, start: u64, seq: &[u8], hp: Option<u8>) -> AlignedRead {
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
            haplotype_tag: hp,
            end: start + len - 1,
        }
    }

    #[test]
    fn test_haplotype_aware_assembly_with_hp_tags() {
        let reference = b"ACGTACGT";
        let reads = vec![
            make_hap_read("r1", 1, b"ACGTACGT", Some(1)),
            make_hap_read("r2", 1, b"ACGTACGT", Some(1)),
            make_hap_read("r3", 1, b"ATGTACGT", Some(2)),
            make_hap_read("r4", 1, b"ATGTACGT", Some(2)),
        ];

        let method = HaplotypeAwareAssembly::new();
        let result = method.assemble_by_haplotype(&reads, reference, 1).unwrap();

        assert_eq!(result.hap1_read_count, 2);
        assert_eq!(result.hap2_read_count, 2);
        assert_eq!(result.hap1_assembly.sequence.len(), reference.len());
        assert_eq!(result.hap2_assembly.sequence.len(), reference.len());

        // Hap1 should match reference at position 2 (C), Hap2 should have T
        assert_eq!(result.hap1_assembly.sequence[1], b'C');
        assert_eq!(result.hap2_assembly.sequence[1], b'T');

        // Divergence should be 1.0 at position 2, 0.0 elsewhere
        assert_eq!(result.divergence[1], 1.0);
        assert_eq!(result.divergence[0], 0.0);
        assert_eq!(result.divergence[2], 0.0);
    }

    #[test]
    fn test_haplotype_aware_assembly_no_haplotypes() {
        let reference = b"ACGT";
        let reads = vec![
            make_hap_read("r1", 1, b"ACGT", None),
            make_hap_read("r2", 1, b"ACGT", None),
        ];

        let method = HaplotypeAwareAssembly::new();
        let result = method.assemble_by_haplotype(&reads, reference, 1).unwrap();

        // Without HP tags and no variants, all reads go to Unassigned
        assert_eq!(result.hap1_read_count, 0);
        assert_eq!(result.hap2_read_count, 0);
    }

    #[test]
    fn test_haplotype_aware_assembly_divergence_length() {
        let reference = b"ACGTACGT";
        let reads = vec![
            make_hap_read("r1", 1, b"ACGTACGT", Some(1)),
            make_hap_read("r2", 1, b"ACGTACGT", Some(2)),
        ];

        let method = HaplotypeAwareAssembly::new();
        let result = method.assemble_by_haplotype(&reads, reference, 1).unwrap();

        assert_eq!(result.divergence.len(), reference.len());
        assert_eq!(result.sv_likelihood.len(), reference.len());
    }
}
