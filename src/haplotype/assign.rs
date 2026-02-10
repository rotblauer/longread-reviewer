use crate::alignment::AlignedRead;

/// Label for haplotype assignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum HaplotypeLabel {
    /// Haplotype 1.
    Hap1,
    /// Haplotype 2.
    Hap2,
    /// Unassigned (could not determine haplotype).
    Unassigned,
}

impl std::fmt::Display for HaplotypeLabel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            HaplotypeLabel::Hap1 => write!(f, "H1"),
            HaplotypeLabel::Hap2 => write!(f, "H2"),
            HaplotypeLabel::Unassigned => write!(f, "U"),
        }
    }
}

/// Assignment of a read to a haplotype.
#[derive(Debug, Clone)]
pub struct ReadAssignment {
    pub read_name: String,
    pub haplotype: HaplotypeLabel,
    /// Confidence of the assignment (0.0-1.0).
    pub confidence: f64,
}

/// Assigns reads to haplotypes.
///
/// Uses the HP tag if available (from haplotagged BAMs), otherwise
/// performs a simple allele-based clustering.
pub struct HaplotypeAssigner;

impl HaplotypeAssigner {
    pub fn new() -> Self {
        Self
    }

    /// Assign reads to haplotypes.
    ///
    /// Strategy:
    /// 1. If reads have HP tags, use those directly.
    /// 2. Otherwise, cluster reads by shared alleles at variant positions.
    pub fn assign(&self, reads: &[AlignedRead], reference: &[u8]) -> Vec<ReadAssignment> {
        // Check if any reads have HP tags
        let has_hp_tags = reads.iter().any(|r| r.haplotype_tag.is_some());

        if has_hp_tags {
            self.assign_from_tags(reads)
        } else {
            self.assign_from_alleles(reads, reference)
        }
    }

    /// Assign using HP tags from haplotagged reads.
    fn assign_from_tags(&self, reads: &[AlignedRead]) -> Vec<ReadAssignment> {
        reads
            .iter()
            .map(|read| {
                let haplotype = match read.haplotype_tag {
                    Some(1) => HaplotypeLabel::Hap1,
                    Some(2) => HaplotypeLabel::Hap2,
                    _ => HaplotypeLabel::Unassigned,
                };
                ReadAssignment {
                    read_name: read.name.clone(),
                    haplotype,
                    confidence: if read.haplotype_tag.is_some() {
                        1.0
                    } else {
                        0.0
                    },
                }
            })
            .collect()
    }

    /// Simple allele-based haplotype clustering.
    ///
    /// Identifies positions where reads disagree and clusters them
    /// into two groups based on which allele they carry.
    fn assign_from_alleles(
        &self,
        reads: &[AlignedRead],
        reference: &[u8],
    ) -> Vec<ReadAssignment> {
        if reads.is_empty() {
            return Vec::new();
        }

        let ref_start = reads
            .iter()
            .filter(|r| !r.cigar.is_empty())
            .map(|r| r.start)
            .min()
            .unwrap_or(1);

        // Find variant positions (where reads disagree with reference).
        let mut variant_positions: Vec<u64> = Vec::new();

        for (pos_idx, &ref_byte) in reference.iter().enumerate() {
            let ref_pos = ref_start + pos_idx as u64;
            let ref_base = ref_byte.to_ascii_uppercase();

            let mut has_alt = false;
            let mut has_ref = false;

            for read in reads {
                for &(rp, base) in &read.reference_aligned_bases() {
                    if rp == ref_pos {
                        if base.to_ascii_uppercase() == ref_base {
                            has_ref = true;
                        } else {
                            has_alt = true;
                        }
                    }
                }
            }

            if has_ref && has_alt {
                variant_positions.push(ref_pos);
            }
        }

        if variant_positions.is_empty() {
            // No variants found - all reads go to Hap1
            return reads
                .iter()
                .map(|r| ReadAssignment {
                    read_name: r.name.clone(),
                    haplotype: HaplotypeLabel::Unassigned,
                    confidence: 0.0,
                })
                .collect();
        }

        // For each read, determine if it carries ref or alt at variant positions.
        // Cluster based on majority allele pattern.
        reads
            .iter()
            .map(|read| {
                let aligned = read.reference_aligned_bases();
                let mut ref_count = 0;
                let mut alt_count = 0;

                for &var_pos in &variant_positions {
                    let pos_idx = (var_pos - ref_start) as usize;
                    let ref_base = if pos_idx < reference.len() {
                        reference[pos_idx].to_ascii_uppercase()
                    } else {
                        continue;
                    };

                    for &(rp, base) in &aligned {
                        if rp == var_pos {
                            if base.to_ascii_uppercase() == ref_base {
                                ref_count += 1;
                            } else {
                                alt_count += 1;
                            }
                        }
                    }
                }

                let total = ref_count + alt_count;
                if total == 0 {
                    ReadAssignment {
                        read_name: read.name.clone(),
                        haplotype: HaplotypeLabel::Unassigned,
                        confidence: 0.0,
                    }
                } else if ref_count > alt_count {
                    ReadAssignment {
                        read_name: read.name.clone(),
                        haplotype: HaplotypeLabel::Hap1,
                        confidence: ref_count as f64 / total as f64,
                    }
                } else {
                    ReadAssignment {
                        read_name: read.name.clone(),
                        haplotype: HaplotypeLabel::Hap2,
                        confidence: alt_count as f64 / total as f64,
                    }
                }
            })
            .collect()
    }
}

impl Default for HaplotypeAssigner {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::CigarOp;

    fn make_read(name: &str, start: u64, seq: &[u8], hp: Option<u8>) -> AlignedRead {
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
    fn test_assign_from_hp_tags() {
        let reads = vec![
            make_read("r1", 1, b"ACGT", Some(1)),
            make_read("r2", 1, b"ACGT", Some(2)),
            make_read("r3", 1, b"ACGT", None),
        ];

        let assigner = HaplotypeAssigner::new();
        let assignments = assigner.assign(&reads, b"ACGT");

        assert_eq!(assignments[0].haplotype, HaplotypeLabel::Hap1);
        assert_eq!(assignments[1].haplotype, HaplotypeLabel::Hap2);
        assert_eq!(assignments[2].haplotype, HaplotypeLabel::Unassigned);
    }

    #[test]
    fn test_assign_from_alleles() {
        let reference = b"ACGT";
        let reads = vec![
            make_read("r1", 1, b"ACGT", None), // matches reference
            make_read("r2", 1, b"ACGT", None), // matches reference
            make_read("r3", 1, b"ATGT", None), // variant at pos 2
            make_read("r4", 1, b"ATGT", None), // variant at pos 2
        ];

        let assigner = HaplotypeAssigner::new();
        let assignments = assigner.assign(&reads, reference);

        // r1, r2 should be one haplotype, r3, r4 the other
        assert_eq!(assignments[0].haplotype, assignments[1].haplotype);
        assert_eq!(assignments[2].haplotype, assignments[3].haplotype);
        assert_ne!(assignments[0].haplotype, assignments[2].haplotype);
    }

    #[test]
    fn test_assign_no_variants() {
        let reference = b"ACGT";
        let reads = vec![
            make_read("r1", 1, b"ACGT", None),
            make_read("r2", 1, b"ACGT", None),
        ];

        let assigner = HaplotypeAssigner::new();
        let assignments = assigner.assign(&reads, reference);

        // All unassigned when no variants
        assert!(assignments
            .iter()
            .all(|a| a.haplotype == HaplotypeLabel::Unassigned));
    }

    #[test]
    fn test_assign_empty() {
        let assigner = HaplotypeAssigner::new();
        let assignments = assigner.assign(&[], b"ACGT");
        assert!(assignments.is_empty());
    }

    #[test]
    fn test_assign_default() {
        let assigner = HaplotypeAssigner::default();
        let assignments = assigner.assign(&[], b"ACGT");
        assert!(assignments.is_empty());
    }

    #[test]
    fn test_assign_hp_tag_values() {
        let reads = vec![
            make_read("r1", 1, b"ACGT", Some(1)),
            make_read("r2", 1, b"ACGT", Some(2)),
            make_read("r3", 1, b"ACGT", Some(3)), // unknown HP value
        ];

        let assigner = HaplotypeAssigner::new();
        let assignments = assigner.assign(&reads, b"ACGT");

        assert_eq!(assignments[0].haplotype, HaplotypeLabel::Hap1);
        assert_eq!(assignments[0].confidence, 1.0);

        assert_eq!(assignments[1].haplotype, HaplotypeLabel::Hap2);
        assert_eq!(assignments[1].confidence, 1.0);

        // HP=3 is not recognized
        assert_eq!(assignments[2].haplotype, HaplotypeLabel::Unassigned);
    }

    #[test]
    fn test_assign_mixed_hp_and_no_hp() {
        // When some reads have HP tags, tag-based assignment is used for all
        let reads = vec![
            make_read("r1", 1, b"ACGT", Some(1)),
            make_read("r2", 1, b"ACGT", None),
        ];

        let assigner = HaplotypeAssigner::new();
        let assignments = assigner.assign(&reads, b"ACGT");

        assert_eq!(assignments[0].haplotype, HaplotypeLabel::Hap1);
        assert_eq!(assignments[1].haplotype, HaplotypeLabel::Unassigned);
        assert_eq!(assignments[1].confidence, 0.0);
    }

    #[test]
    fn test_assign_allele_confidence() {
        let reference = b"ACGT";
        let reads = vec![
            make_read("r1", 1, b"ACGT", None),
            make_read("r2", 1, b"ATGT", None),
        ];

        let assigner = HaplotypeAssigner::new();
        let assignments = assigner.assign(&reads, reference);

        // Both should have confidence of 1.0 (100% ref or 100% alt at variant)
        assert_eq!(assignments[0].confidence, 1.0);
        assert_eq!(assignments[1].confidence, 1.0);
    }

    #[test]
    fn test_assign_multiple_variant_positions() {
        let reference = b"AAAAA";
        // Reads with variants at multiple positions
        let reads = vec![
            make_read("r1", 1, b"AAAAA", None), // ref
            make_read("r2", 1, b"ATAAA", None), // variant at pos 2
            make_read("r3", 1, b"AAATA", None), // variant at pos 4
            make_read("r4", 1, b"ATATA", None), // variants at pos 2 and 4
        ];

        let assigner = HaplotypeAssigner::new();
        let assignments = assigner.assign(&reads, reference);

        // r1 should be Hap1 (reference alleles)
        assert_eq!(assignments[0].haplotype, HaplotypeLabel::Hap1);
        // r4 should be Hap2 (all alt alleles)
        assert_eq!(assignments[3].haplotype, HaplotypeLabel::Hap2);
    }

    #[test]
    fn test_haplotype_label_display() {
        assert_eq!(HaplotypeLabel::Hap1.to_string(), "H1");
        assert_eq!(HaplotypeLabel::Hap2.to_string(), "H2");
        assert_eq!(HaplotypeLabel::Unassigned.to_string(), "U");
    }

    #[test]
    fn test_read_assignment_fields() {
        let reads = vec![
            make_read("test_read", 1, b"ACGT", Some(1)),
        ];

        let assigner = HaplotypeAssigner::new();
        let assignments = assigner.assign(&reads, b"ACGT");

        assert_eq!(assignments[0].read_name, "test_read");
        assert_eq!(assignments[0].haplotype, HaplotypeLabel::Hap1);
        assert_eq!(assignments[0].confidence, 1.0);
    }

    #[test]
    fn test_assign_all_same_allele() {
        let reference = b"AAAA";
        // All reads have the same variant -> no variant positions detected
        let reads = vec![
            make_read("r1", 1, b"TTTT", None),
            make_read("r2", 1, b"TTTT", None),
        ];

        let assigner = HaplotypeAssigner::new();
        let assignments = assigner.assign(&reads, reference);

        // No variant positions (need both ref AND alt), all unassigned
        assert!(assignments.iter().all(|a| a.haplotype == HaplotypeLabel::Unassigned));
    }

    #[test]
    fn test_haplotype_label_equality() {
        assert_eq!(HaplotypeLabel::Hap1, HaplotypeLabel::Hap1);
        assert_ne!(HaplotypeLabel::Hap1, HaplotypeLabel::Hap2);
        assert_ne!(HaplotypeLabel::Hap1, HaplotypeLabel::Unassigned);
    }
}
