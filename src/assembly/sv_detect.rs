//! Assembly-based structural variant detection.
//!
//! Compares a local assembly to the reference sequence to detect structural
//! variants without relying on the original aligner's CIGAR strings.  This
//! provides an independent, assembly-first signal that is combined with
//! read-level evidence to produce a per-event confidence score.

use crate::alignment::AlignedRead;
use crate::assembly::method::{AssemblyResult, HaplotypeAssemblyResult};

/// A structural variant detected by comparing assembly to reference.
#[derive(Debug, Clone)]
pub struct AssemblySVEvent {
    /// Type of structural variant.
    pub sv_type: AssemblySVType,
    /// 1-based start position on the reference.
    pub start: u64,
    /// 1-based end position on the reference (inclusive).
    pub end: u64,
    /// Size of the event in base pairs (positive for insertions, negative for deletions).
    pub size: i64,
    /// Confidence score (0.0–1.0) combining assembly divergence, depth, and
    /// haplotype concordance.
    pub confidence: f64,
    /// Number of reads supporting this event (via CIGAR evidence).
    pub supporting_reads: usize,
    /// Source of detection.
    pub source: DetectionSource,
}

/// Structural variant type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AssemblySVType {
    Deletion,
    Insertion,
}

impl std::fmt::Display for AssemblySVType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AssemblySVType::Deletion => write!(f, "DEL"),
            AssemblySVType::Insertion => write!(f, "INS"),
        }
    }
}

/// How this SV was detected.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DetectionSource {
    /// Detected from assembly vs reference divergence.
    Assembly,
    /// Detected from haplotype divergence.
    HaplotypeDivergence,
    /// Detected from CIGAR-level read evidence.
    CigarReads,
    /// Multiple sources agree.
    Combined,
}

impl std::fmt::Display for DetectionSource {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DetectionSource::Assembly => write!(f, "Assembly"),
            DetectionSource::HaplotypeDivergence => write!(f, "HapDiv"),
            DetectionSource::CigarReads => write!(f, "CIGAR"),
            DetectionSource::Combined => write!(f, "Combined"),
        }
    }
}

/// Detects structural variants by comparing assemblies to the reference.
pub struct AssemblySVCaller {
    /// Minimum number of consecutive divergent bases to call an SV.
    pub min_divergent_run: usize,
    /// Minimum SV size (bp) to report.
    pub min_sv_size: usize,
    /// Minimum confidence to report.
    pub min_confidence: f64,
}

impl AssemblySVCaller {
    pub fn new() -> Self {
        Self {
            min_divergent_run: 10,
            min_sv_size: 30,
            min_confidence: 0.1,
        }
    }

    /// Detect SVs from assembly vs reference comparison.
    ///
    /// Scans for runs of low-confidence or divergent positions in the assembly
    /// and interprets them as candidate structural variants.
    pub fn detect_from_assembly(
        &self,
        assembly: &AssemblyResult,
        reference: &[u8],
        ref_start_pos: u64,
        reads: &[AlignedRead],
    ) -> Vec<AssemblySVEvent> {
        let mut events = Vec::new();
        let len = assembly.sequence.len().min(reference.len());
        if len == 0 {
            return events;
        }

        // Find runs of divergent positions (assembly ≠ reference)
        let mut run_start: Option<usize> = None;
        for i in 0..=len {
            let is_divergent = if i < len {
                !assembly.sequence[i].eq_ignore_ascii_case(&reference[i])
            } else {
                false // sentinel to close the last run
            };

            if is_divergent {
                if run_start.is_none() {
                    run_start = Some(i);
                }
            } else if let Some(start) = run_start {
                let run_len = i - start;
                if run_len >= self.min_divergent_run {
                    // Determine SV type from depth pattern
                    let avg_depth: f64 = assembly.depth[start..i]
                        .iter()
                        .map(|&d| d as f64)
                        .sum::<f64>()
                        / run_len as f64;

                    let avg_confidence: f64 = assembly.confidence[start..i]
                        .iter()
                        .sum::<f64>()
                        / run_len as f64;

                    // Count CIGAR-based supporting reads
                    let genomic_start = ref_start_pos + start as u64;
                    let genomic_end = ref_start_pos + i as u64;
                    let supporting = count_cigar_support(
                        reads,
                        genomic_start,
                        genomic_end,
                        self.min_sv_size,
                    );

                    // Classify: low depth in the divergent region suggests deletion
                    let sv_type = if avg_depth < 2.0 {
                        AssemblySVType::Deletion
                    } else {
                        AssemblySVType::Insertion
                    };

                    let size = if sv_type == AssemblySVType::Deletion {
                        -(run_len as i64)
                    } else {
                        run_len as i64
                    };

                    // Confidence: combine assembly signal and CIGAR support
                    let asm_conf = 1.0 - avg_confidence.min(1.0);
                    let read_conf = (supporting as f64 / reads.len().max(1) as f64).min(1.0);
                    let confidence = (asm_conf * 0.5 + read_conf * 0.5).min(1.0);

                    if size.unsigned_abs() as usize >= self.min_sv_size
                        && confidence >= self.min_confidence
                    {
                        events.push(AssemblySVEvent {
                            sv_type,
                            start: genomic_start,
                            end: genomic_end,
                            size,
                            confidence,
                            supporting_reads: supporting,
                            source: if supporting > 0 {
                                DetectionSource::Combined
                            } else {
                                DetectionSource::Assembly
                            },
                        });
                    }
                }
                run_start = None;
            }
        }

        events
    }

    /// Detect SVs from haplotype divergence tracks.
    ///
    /// Scans the divergence and SV-likelihood vectors for runs of high
    /// divergence, indicating heterozygous structural variants.
    pub fn detect_from_haplotype_assembly(
        &self,
        hap_result: &HaplotypeAssemblyResult,
        ref_start_pos: u64,
        reads: &[AlignedRead],
    ) -> Vec<AssemblySVEvent> {
        let mut events = Vec::new();
        let len = hap_result.divergence.len();
        if len == 0 {
            return events;
        }

        let mut run_start: Option<usize> = None;
        for i in 0..=len {
            let is_divergent = if i < len {
                hap_result.divergence[i] > 0.0
            } else {
                false
            };

            if is_divergent {
                if run_start.is_none() {
                    run_start = Some(i);
                }
            } else if let Some(start) = run_start {
                let run_len = i - start;
                if run_len >= self.min_divergent_run {
                    let avg_sv_likelihood: f64 = hap_result.sv_likelihood[start..i]
                        .iter()
                        .sum::<f64>()
                        / run_len as f64;

                    let genomic_start = ref_start_pos + start as u64;
                    let genomic_end = ref_start_pos + i as u64;

                    let supporting = count_cigar_support(
                        reads,
                        genomic_start,
                        genomic_end,
                        self.min_sv_size,
                    );

                    // Use depth difference between haplotypes to classify
                    let h1_depth: f64 = hap_result.hap1_assembly.depth[start..i]
                        .iter()
                        .map(|&d| d as f64)
                        .sum::<f64>()
                        / run_len as f64;
                    let h2_depth: f64 = hap_result.hap2_assembly.depth[start..i]
                        .iter()
                        .map(|&d| d as f64)
                        .sum::<f64>()
                        / run_len as f64;

                    let sv_type = if h1_depth < 1.0 || h2_depth < 1.0 {
                        AssemblySVType::Deletion
                    } else {
                        AssemblySVType::Insertion
                    };

                    let size = if sv_type == AssemblySVType::Deletion {
                        -(run_len as i64)
                    } else {
                        run_len as i64
                    };

                    let confidence = avg_sv_likelihood.min(1.0);

                    if size.unsigned_abs() as usize >= self.min_sv_size
                        && confidence >= self.min_confidence
                    {
                        events.push(AssemblySVEvent {
                            sv_type,
                            start: genomic_start,
                            end: genomic_end,
                            size,
                            confidence,
                            supporting_reads: supporting,
                            source: if supporting > 0 {
                                DetectionSource::Combined
                            } else {
                                DetectionSource::HaplotypeDivergence
                            },
                        });
                    }
                }
                run_start = None;
            }
        }

        events
    }
}

impl Default for AssemblySVCaller {
    fn default() -> Self {
        Self::new()
    }
}

/// Count reads with large CIGAR indels overlapping the given genomic range.
fn count_cigar_support(
    reads: &[AlignedRead],
    region_start: u64,
    region_end: u64,
    min_size: usize,
) -> usize {
    use crate::alignment::CigarOp;

    let mut count = 0;
    for read in reads {
        let mut ref_pos = read.start;
        for op in &read.cigar {
            match op {
                CigarOp::Match(n) => ref_pos += *n as u64,
                CigarOp::Insertion(n) if *n as usize >= min_size => {
                    if ref_pos >= region_start && ref_pos <= region_end {
                        count += 1;
                        break;
                    }
                }
                CigarOp::Deletion(n) if *n as usize >= min_size => {
                    let del_end = ref_pos + *n as u64;
                    if ref_pos <= region_end && del_end >= region_start {
                        count += 1;
                        break;
                    }
                    ref_pos += *n as u64;
                }
                CigarOp::Deletion(n) => ref_pos += *n as u64,
                _ => {}
            }
        }
    }
    count
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
    fn test_detect_from_assembly_no_divergence() {
        let reference = b"ACGTACGT";
        let assembly = AssemblyResult {
            sequence: b"ACGTACGT".to_vec(),
            depth: vec![10; 8],
            confidence: vec![1.0; 8],
            method_name: "test".to_string(),
        };

        let caller = AssemblySVCaller::new();
        let events = caller.detect_from_assembly(&assembly, reference, 1, &[]);
        assert!(events.is_empty(), "No SVs expected when assembly matches reference");
    }

    #[test]
    fn test_detect_from_assembly_with_divergent_run() {
        let ref_len = 100;
        let reference: Vec<u8> = (0..ref_len).map(|_| b'A').collect();

        // Assembly diverges for 40 bases in the middle (positions 30..70)
        let mut sequence = reference.clone();
        for i in 30..70 {
            sequence[i] = b'T';
        }
        let assembly = AssemblyResult {
            sequence,
            depth: vec![1; ref_len], // low depth suggests deletion
            confidence: vec![0.2; ref_len],
            method_name: "test".to_string(),
        };

        let reads = vec![make_read("r1", 1, &reference)];
        let caller = AssemblySVCaller::new();
        let events = caller.detect_from_assembly(&assembly, &reference, 1, &reads);

        assert!(
            !events.is_empty(),
            "Should detect at least one SV from divergent region"
        );
        assert!(events[0].size.unsigned_abs() >= 30);
    }

    #[test]
    fn test_assembly_sv_type_display() {
        assert_eq!(format!("{}", AssemblySVType::Deletion), "DEL");
        assert_eq!(format!("{}", AssemblySVType::Insertion), "INS");
    }

    #[test]
    fn test_detection_source_display() {
        assert_eq!(format!("{}", DetectionSource::Assembly), "Assembly");
        assert_eq!(format!("{}", DetectionSource::Combined), "Combined");
    }
}
