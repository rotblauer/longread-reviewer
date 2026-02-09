use crate::alignment::AlignedRead;
use crate::assembly::method::AssemblyResult;

/// Per-base quality metric at a single reference position.
#[derive(Debug, Clone)]
pub struct BaseMetric {
    /// 1-based reference position.
    pub position: u64,
    /// Reference base at this position.
    pub ref_base: u8,
    /// Consensus/assembly base at this position.
    pub assembly_base: u8,
    /// Read depth at this position.
    pub depth: u32,
    /// Fraction of reads agreeing with the assembly base (0.0-1.0).
    pub agreement: f64,
    /// Average base quality of reads at this position.
    pub avg_quality: f64,
    /// Whether assembly differs from reference.
    pub is_variant: bool,
}

/// Overall fitness score for an assembly.
#[derive(Debug, Clone)]
pub struct FitnessScore {
    /// Per-base metrics.
    pub base_metrics: Vec<BaseMetric>,
    /// Average agreement across all positions (0.0-1.0).
    pub mean_agreement: f64,
    /// Average depth across the assembly.
    pub mean_depth: f64,
    /// Fraction of positions where assembly matches reference.
    pub reference_identity: f64,
    /// Overall fitness score (0.0-1.0), combining all metrics.
    pub overall: f64,
}

/// Calculator for per-base metrics and fitness scores.
pub struct MetricsCalculator {
    /// Weight for agreement score in overall fitness.
    pub agreement_weight: f64,
    /// Weight for depth score in overall fitness.
    pub depth_weight: f64,
    /// Weight for quality score in overall fitness.
    pub quality_weight: f64,
}

impl MetricsCalculator {
    pub fn new() -> Self {
        Self {
            agreement_weight: 0.5,
            depth_weight: 0.3,
            quality_weight: 0.2,
        }
    }

    pub fn with_weights(agreement: f64, depth: f64, quality: f64) -> Self {
        Self {
            agreement_weight: agreement,
            depth_weight: depth,
            quality_weight: quality,
        }
    }

    /// Compute per-base metrics and overall fitness for an assembly.
    pub fn compute_fitness(
        &self,
        assembly: &AssemblyResult,
        reads: &[AlignedRead],
        reference: &[u8],
    ) -> FitnessScore {
        let len = assembly.sequence.len();
        if len == 0 {
            return FitnessScore {
                base_metrics: Vec::new(),
                mean_agreement: 0.0,
                mean_depth: 0.0,
                reference_identity: 0.0,
                overall: 0.0,
            };
        }

        let ref_start = reads
            .iter()
            .filter(|r| !r.cigar.is_empty())
            .map(|r| r.start)
            .min()
            .unwrap_or(1);

        let mut base_metrics = Vec::with_capacity(len);

        for i in 0..len {
            let ref_base = if i < reference.len() {
                reference[i]
            } else {
                b'N'
            };
            let assembly_base = assembly.sequence[i];
            let depth = assembly.depth[i];
            let confidence = assembly.confidence[i];

            // Count reads agreeing with assembly base at this position.
            let ref_pos = ref_start + i as u64;
            let mut agree_count = 0u32;
            let mut total_count = 0u32;
            let mut quality_sum = 0.0f64;

            for read in reads {
                for &(rp, base) in &read.reference_aligned_bases() {
                    if rp == ref_pos {
                        total_count += 1;
                        if base.eq_ignore_ascii_case(&assembly_base) {
                            agree_count += 1;
                        }
                        // Use the average quality of the read as a proxy
                        let avg_q = if read.qualities.is_empty() {
                            30.0
                        } else {
                            read.qualities.iter().map(|&q| q as f64).sum::<f64>()
                                / read.qualities.len() as f64
                        };
                        quality_sum += avg_q;
                    }
                }
            }

            let agreement = if total_count > 0 {
                agree_count as f64 / total_count as f64
            } else {
                confidence
            };

            let avg_quality = if total_count > 0 {
                quality_sum / total_count as f64
            } else {
                0.0
            };

            let is_variant =
                !ref_base.eq_ignore_ascii_case(&assembly_base);

            base_metrics.push(BaseMetric {
                position: ref_pos,
                ref_base,
                assembly_base,
                depth,
                agreement,
                avg_quality,
                is_variant,
            });
        }

        let mean_agreement = if base_metrics.is_empty() {
            0.0
        } else {
            base_metrics.iter().map(|m| m.agreement).sum::<f64>() / base_metrics.len() as f64
        };

        let mean_depth = if base_metrics.is_empty() {
            0.0
        } else {
            base_metrics.iter().map(|m| m.depth as f64).sum::<f64>() / base_metrics.len() as f64
        };

        let reference_identity = if base_metrics.is_empty() {
            0.0
        } else {
            let matches = base_metrics.iter().filter(|m| !m.is_variant).count();
            matches as f64 / base_metrics.len() as f64
        };

        // Normalize depth score (saturates at 30x).
        let depth_score = (mean_depth / 30.0).min(1.0);

        // Normalize quality score (Phred 30 = 1.0).
        let mean_quality = if base_metrics.is_empty() {
            0.0
        } else {
            base_metrics.iter().map(|m| m.avg_quality).sum::<f64>() / base_metrics.len() as f64
        };
        let quality_score = (mean_quality / 30.0).min(1.0);

        let overall = self.agreement_weight * mean_agreement
            + self.depth_weight * depth_score
            + self.quality_weight * quality_score;

        FitnessScore {
            base_metrics,
            mean_agreement,
            mean_depth,
            reference_identity,
            overall,
        }
    }
}

impl Default for MetricsCalculator {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::CigarOp;
    use crate::assembly::method::AssemblyResult;

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
    fn test_perfect_assembly_fitness() {
        let reference = b"ACGT";
        let reads = vec![
            make_read("r1", 1, b"ACGT"),
            make_read("r2", 1, b"ACGT"),
            make_read("r3", 1, b"ACGT"),
        ];
        let assembly = AssemblyResult {
            sequence: b"ACGT".to_vec(),
            depth: vec![3, 3, 3, 3],
            confidence: vec![1.0, 1.0, 1.0, 1.0],
            method_name: "test".to_string(),
        };

        let calc = MetricsCalculator::new();
        let fitness = calc.compute_fitness(&assembly, &reads, reference);

        assert_eq!(fitness.mean_agreement, 1.0);
        assert_eq!(fitness.reference_identity, 1.0);
        assert!(fitness.overall > 0.5);
    }

    #[test]
    fn test_variant_detection() {
        let reference = b"ACGT";
        let reads = vec![
            make_read("r1", 1, b"ACGT"),
            make_read("r2", 1, b"ATGT"),
            make_read("r3", 1, b"ATGT"),
        ];
        // Assembly with a variant at position 2 (C->T)
        let assembly = AssemblyResult {
            sequence: b"ATGT".to_vec(),
            depth: vec![3, 3, 3, 3],
            confidence: vec![1.0, 0.67, 1.0, 1.0],
            method_name: "test".to_string(),
        };

        let calc = MetricsCalculator::new();
        let fitness = calc.compute_fitness(&assembly, &reads, reference);

        // Position 2 should be marked as variant
        assert!(fitness.base_metrics[1].is_variant);
        assert!(!fitness.base_metrics[0].is_variant);
        assert!(fitness.reference_identity < 1.0);
    }

    #[test]
    fn test_empty_assembly() {
        let calc = MetricsCalculator::new();
        let assembly = AssemblyResult {
            sequence: Vec::new(),
            depth: Vec::new(),
            confidence: Vec::new(),
            method_name: "test".to_string(),
        };

        let fitness = calc.compute_fitness(&assembly, &[], b"");
        assert_eq!(fitness.overall, 0.0);
    }

    #[test]
    fn test_custom_weights() {
        let reference = b"ACGT";
        let reads = vec![make_read("r1", 1, b"ACGT")];
        let assembly = AssemblyResult {
            sequence: b"ACGT".to_vec(),
            depth: vec![1, 1, 1, 1],
            confidence: vec![1.0, 1.0, 1.0, 1.0],
            method_name: "test".to_string(),
        };

        let calc1 = MetricsCalculator::with_weights(1.0, 0.0, 0.0);
        let calc2 = MetricsCalculator::with_weights(0.0, 1.0, 0.0);

        let f1 = calc1.compute_fitness(&assembly, &reads, reference);
        let f2 = calc2.compute_fitness(&assembly, &reads, reference);

        // With 100% agreement weight, score should be high (agreement = 1.0)
        assert!(f1.overall > f2.overall);
    }
}
