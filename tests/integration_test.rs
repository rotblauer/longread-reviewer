//! Integration tests using the example NA19240 BAM and reference FASTA.
//!
//! These tests exercise the full pipeline: BAM reading → assembly → haplotype → metrics,
//! using the synthetic test data included in `tests/data/`.

use std::path::PathBuf;

use longread_reviewer::alignment::{AlignedRead, AlignmentReader, CigarOp};
use longread_reviewer::assembly::engine::AssemblyEngine;
use longread_reviewer::assembly::method::{
    AssemblyMethod, ConsensusAssembly, WindowConsensusAssembly,
};
use longread_reviewer::haplotype::{HaplotypeAssigner, HaplotypeLabel};
use longread_reviewer::metrics::MetricsCalculator;
use longread_reviewer::reference::ReferenceGenome;
use longread_reviewer::region::Region;

// ---------------------------------------------------------------------------
// helpers
// ---------------------------------------------------------------------------

fn test_bam_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/NA19240.chr17_fragment.bam")
}

fn test_ref_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data/chr17_fragment.fa")
}

/// Region covering the full test fragment (1-based, inclusive).
fn full_region() -> Region {
    Region::new("chr17", 1, 2001).unwrap()
}

/// A smaller sub-region for targeted tests.
fn sub_region() -> Region {
    Region::new("chr17", 100, 500).unwrap()
}

/// Load reads + reference for the full region.
fn load_full() -> (Vec<AlignedRead>, Vec<u8>) {
    let region = full_region();
    let reads = AlignmentReader::read_bam(&test_bam_path(), &region)
        .expect("failed to read BAM");
    let reference = ReferenceGenome::from_file(&test_ref_path())
        .expect("failed to read reference");
    let ref_seq = reference.fetch(&region).expect("failed to fetch region");
    (reads, ref_seq)
}

/// Load reads + reference for the sub-region.
fn load_sub() -> (Vec<AlignedRead>, Vec<u8>) {
    let region = sub_region();
    let reads = AlignmentReader::read_bam(&test_bam_path(), &region)
        .expect("failed to read BAM");
    let reference = ReferenceGenome::from_file(&test_ref_path())
        .expect("failed to read reference");
    let ref_seq = reference.fetch(&region).expect("failed to fetch region");
    (reads, ref_seq)
}

// ===========================================================================
// 1. File & data loading tests
// ===========================================================================

#[test]
fn test_bam_file_exists() {
    assert!(
        test_bam_path().exists(),
        "Test BAM not found at {:?}",
        test_bam_path()
    );
}

#[test]
fn test_bam_index_exists() {
    let bai = test_bam_path().with_extension("bam.bai");
    assert!(
        bai.exists(),
        "BAM index not found at {:?}",
        bai
    );
}

#[test]
fn test_reference_file_exists() {
    assert!(
        test_ref_path().exists(),
        "Reference FASTA not found at {:?}",
        test_ref_path()
    );
}

#[test]
fn test_reference_index_exists() {
    let fai = PathBuf::from(format!("{}.fai", test_ref_path().display()));
    assert!(
        fai.exists(),
        "FASTA index not found at {:?}",
        fai
    );
}

// ===========================================================================
// 2. BAM reading tests
// ===========================================================================

#[test]
fn test_read_bam_full_region() {
    let region = full_region();
    let reads = AlignmentReader::read_bam(&test_bam_path(), &region).unwrap();
    assert!(
        reads.len() >= 10,
        "Expected at least 10 reads, got {}",
        reads.len()
    );
    assert!(
        reads.len() <= 30,
        "Expected at most 30 reads, got {}",
        reads.len()
    );
}

#[test]
fn test_read_bam_sub_region() {
    let region = sub_region();
    let reads = AlignmentReader::read_bam(&test_bam_path(), &region).unwrap();
    assert!(
        !reads.is_empty(),
        "Expected reads in sub-region, got none"
    );
}

#[test]
fn test_read_bam_empty_region() {
    // Query a region with no reads — well beyond our fragment
    let region = Region::new("chr17", 100_000_000, 100_001_000).unwrap();
    let reads = AlignmentReader::read_bam(&test_bam_path(), &region).unwrap();
    assert!(
        reads.is_empty(),
        "Expected no reads in empty region"
    );
}

#[test]
fn test_all_reads_on_correct_chrom() {
    let (reads, _) = load_full();
    for read in &reads {
        assert_eq!(
            read.chrom, "chr17",
            "Read {} on wrong chrom: {}",
            read.name, read.chrom
        );
    }
}

#[test]
fn test_read_names_are_nonempty() {
    let (reads, _) = load_full();
    for read in &reads {
        assert!(
            !read.name.is_empty(),
            "Read has empty name"
        );
    }
}

#[test]
fn test_read_names_look_like_pacbio_ccs() {
    let (reads, _) = load_full();
    for read in &reads {
        assert!(
            read.name.contains("/ccs"),
            "Read name does not look like PacBio CCS: {}",
            read.name
        );
    }
}

#[test]
fn test_reads_have_valid_coordinates() {
    let (reads, _) = load_full();
    for read in &reads {
        assert!(
            read.end >= read.start,
            "Read {} has end < start: {}-{}",
            read.name, read.start, read.end
        );
    }
}

#[test]
fn test_reads_have_cigar() {
    let (reads, _) = load_full();
    for read in &reads {
        assert!(
            !read.cigar.is_empty(),
            "Read {} has empty CIGAR",
            read.name
        );
    }
}

#[test]
fn test_reads_have_sequence() {
    let (reads, _) = load_full();
    for read in &reads {
        assert!(
            !read.sequence.is_empty(),
            "Read {} has empty sequence",
            read.name
        );
    }
}

#[test]
fn test_reads_have_qualities() {
    let (reads, _) = load_full();
    for read in &reads {
        assert!(
            !read.qualities.is_empty(),
            "Read {} has empty qualities",
            read.name
        );
        assert_eq!(
            read.qualities.len(),
            read.sequence.len(),
            "Read {} quality/sequence length mismatch: {} vs {}",
            read.name,
            read.qualities.len(),
            read.sequence.len()
        );
    }
}

#[test]
fn test_reads_have_valid_mapq() {
    let (reads, _) = load_full();
    for read in &reads {
        assert!(
            read.mapq > 0,
            "Read {} has zero MAPQ",
            read.name
        );
    }
}

#[test]
fn test_some_reads_have_hp_tags() {
    let (reads, _) = load_full();
    let tagged = reads.iter().filter(|r| r.haplotype_tag.is_some()).count();
    assert!(
        tagged > 0,
        "Expected some reads with HP tags, got 0"
    );
}

#[test]
fn test_hp_tags_are_valid() {
    let (reads, _) = load_full();
    for read in &reads {
        if let Some(hp) = read.haplotype_tag {
            assert!(
                hp == 1 || hp == 2,
                "Read {} has invalid HP tag: {}",
                read.name, hp
            );
        }
    }
}

#[test]
fn test_both_haplotypes_present() {
    let (reads, _) = load_full();
    let hap1 = reads.iter().filter(|r| r.haplotype_tag == Some(1)).count();
    let hap2 = reads.iter().filter(|r| r.haplotype_tag == Some(2)).count();
    assert!(hap1 > 0, "No haplotype 1 reads");
    assert!(hap2 > 0, "No haplotype 2 reads");
}

#[test]
fn test_some_reads_are_untagged() {
    let (reads, _) = load_full();
    let untagged = reads.iter().filter(|r| r.haplotype_tag.is_none()).count();
    assert!(untagged > 0, "Expected some untagged reads");
}

#[test]
fn test_cigar_ops_are_valid_types() {
    let (reads, _) = load_full();
    for read in &reads {
        for op in &read.cigar {
            match op {
                CigarOp::Match(n) | CigarOp::Insertion(n) | CigarOp::Deletion(n)
                | CigarOp::SoftClip(n) | CigarOp::HardClip(n) => {
                    assert!(*n > 0, "Read {} has zero-length CIGAR op", read.name);
                }
            }
        }
    }
}

#[test]
fn test_cigar_read_length_matches_sequence() {
    let (reads, _) = load_full();
    for read in &reads {
        let cigar_read_len: u32 = read.cigar.iter().map(|op| op.read_len()).sum();
        assert_eq!(
            cigar_read_len as usize,
            read.sequence.len(),
            "Read {} CIGAR read length ({}) != sequence length ({})",
            read.name,
            cigar_read_len,
            read.sequence.len()
        );
    }
}

#[test]
fn test_aligned_sequence_nonempty() {
    let (reads, _) = load_full();
    for read in &reads {
        let aligned = read.aligned_sequence();
        assert!(
            !aligned.is_empty(),
            "Read {} has empty aligned sequence",
            read.name
        );
    }
}

#[test]
fn test_reference_aligned_bases_nonempty() {
    let (reads, _) = load_full();
    for read in &reads {
        let bases = read.reference_aligned_bases();
        assert!(
            !bases.is_empty(),
            "Read {} produced no reference-aligned bases",
            read.name
        );
    }
}

#[test]
fn test_reference_aligned_bases_positions_ascending() {
    let (reads, _) = load_full();
    for read in &reads {
        let bases = read.reference_aligned_bases();
        for window in bases.windows(2) {
            assert!(
                window[0].0 <= window[1].0,
                "Read {} has non-ascending ref positions: {} > {}",
                read.name,
                window[0].0,
                window[1].0
            );
        }
    }
}

#[test]
fn test_reference_aligned_bases_are_valid_nucleotides() {
    let (reads, _) = load_full();
    let valid = b"ACGTNacgtn";
    for read in &reads {
        let bases = read.reference_aligned_bases();
        for &(pos, base) in &bases {
            assert!(
                valid.contains(&base),
                "Read {} has invalid base '{}' at position {}",
                read.name,
                base as char,
                pos
            );
        }
    }
}

// ===========================================================================
// 3. Reference loading tests
// ===========================================================================

#[test]
fn test_reference_loads() {
    let reference = ReferenceGenome::from_file(&test_ref_path()).unwrap();
    let chroms = reference.chromosomes();
    assert!(
        chroms.contains(&"chr17"),
        "Reference missing chr17"
    );
}

#[test]
fn test_reference_fetch_full_region() {
    let reference = ReferenceGenome::from_file(&test_ref_path()).unwrap();
    let region = full_region();
    let seq = reference.fetch(&region).unwrap();
    assert_eq!(
        seq.len(),
        2001,
        "Reference sequence length mismatch"
    );
}

#[test]
fn test_reference_fetch_sub_region() {
    let reference = ReferenceGenome::from_file(&test_ref_path()).unwrap();
    let region = sub_region();
    let seq = reference.fetch(&region).unwrap();
    assert_eq!(
        seq.len() as u64,
        region.len(),
        "Sub-region fetch length mismatch"
    );
}

#[test]
fn test_reference_bases_are_valid() {
    let reference = ReferenceGenome::from_file(&test_ref_path()).unwrap();
    let seq = reference.fetch(&full_region()).unwrap();
    let valid = b"ACGTN";
    for (i, &base) in seq.iter().enumerate() {
        assert!(
            valid.contains(&base),
            "Invalid reference base '{}' at position {}",
            base as char,
            i
        );
    }
}

#[test]
fn test_reference_fetch_wrong_chrom() {
    let reference = ReferenceGenome::from_file(&test_ref_path()).unwrap();
    let region = Region::new("chrZ", 1, 100).unwrap();
    assert!(reference.fetch(&region).is_err());
}

// ===========================================================================
// 4. Consensus assembly tests
// ===========================================================================

#[test]
fn test_consensus_assembly_produces_output() {
    let (reads, ref_seq) = load_full();
    let method = ConsensusAssembly;
    let result = method.assemble(&reads, &ref_seq).unwrap();

    assert_eq!(
        result.sequence.len(),
        ref_seq.len(),
        "Assembly length should match reference length"
    );
    assert_eq!(result.depth.len(), ref_seq.len());
    assert_eq!(result.confidence.len(), ref_seq.len());
    assert_eq!(result.method_name, "consensus");
}

#[test]
fn test_consensus_assembly_has_depth() {
    let (reads, ref_seq) = load_full();
    let method = ConsensusAssembly;
    let result = method.assemble(&reads, &ref_seq).unwrap();

    let mean_depth: f64 =
        result.depth.iter().map(|&d| d as f64).sum::<f64>() / result.depth.len() as f64;
    assert!(
        mean_depth > 1.0,
        "Expected mean depth > 1.0, got {:.1}",
        mean_depth
    );
}

#[test]
fn test_consensus_confidence_in_range() {
    let (reads, ref_seq) = load_full();
    let method = ConsensusAssembly;
    let result = method.assemble(&reads, &ref_seq).unwrap();

    for (i, &conf) in result.confidence.iter().enumerate() {
        assert!(
            (0.0..=1.0).contains(&conf),
            "Confidence at position {} out of range: {}",
            i, conf
        );
    }
}

#[test]
fn test_consensus_bases_are_valid() {
    let (reads, ref_seq) = load_full();
    let method = ConsensusAssembly;
    let result = method.assemble(&reads, &ref_seq).unwrap();

    let valid = b"ACGTNacgtn";
    for (i, &base) in result.sequence.iter().enumerate() {
        assert!(
            valid.contains(&base),
            "Invalid assembly base '{}' at position {}",
            base as char,
            i
        );
    }
}

#[test]
fn test_consensus_high_identity_to_reference() {
    let (reads, ref_seq) = load_full();
    let method = ConsensusAssembly;
    let result = method.assemble(&reads, &ref_seq).unwrap();

    let matches = result
        .sequence
        .iter()
        .zip(ref_seq.iter())
        .filter(|(a, b)| a.eq_ignore_ascii_case(b))
        .count();
    let identity = matches as f64 / ref_seq.len() as f64;
    assert!(
        identity > 0.9,
        "Expected >90% identity to reference, got {:.1}%",
        identity * 100.0
    );
}

// ===========================================================================
// 5. Window consensus assembly tests
// ===========================================================================

#[test]
fn test_window_consensus_assembly_produces_output() {
    let (reads, ref_seq) = load_full();
    let method = WindowConsensusAssembly::new(100, 20);
    let result = method.assemble(&reads, &ref_seq).unwrap();

    assert_eq!(result.sequence.len(), ref_seq.len());
    assert_eq!(result.method_name, "window_consensus");
}

#[test]
fn test_window_consensus_various_sizes() {
    let (reads, ref_seq) = load_full();

    for (window, overlap) in &[(50, 10), (100, 20), (200, 40), (500, 100)] {
        let method = WindowConsensusAssembly::new(*window, *overlap);
        let result = method.assemble(&reads, &ref_seq).unwrap();

        assert_eq!(
            result.sequence.len(),
            ref_seq.len(),
            "Window({},{}) produced wrong length",
            window, overlap
        );
        assert!(
            result.depth.iter().any(|&d| d > 0),
            "Window({},{}) has no depth",
            window, overlap
        );
    }
}

#[test]
fn test_window_consensus_confidence_in_range() {
    let (reads, ref_seq) = load_full();
    let method = WindowConsensusAssembly::new(100, 20);
    let result = method.assemble(&reads, &ref_seq).unwrap();

    for (i, &conf) in result.confidence.iter().enumerate() {
        assert!(
            (0.0..=1.0).contains(&conf),
            "Confidence at position {} out of range: {}",
            i, conf
        );
    }
}

// ===========================================================================
// 6. Assembly engine tests
// ===========================================================================

#[test]
fn test_engine_evaluate_all_with_test_data() {
    let (reads, ref_seq) = load_full();
    let region = full_region();

    let mut engine = AssemblyEngine::new();
    engine.add_method(Box::new(ConsensusAssembly));
    engine.add_method(Box::new(WindowConsensusAssembly::new(100, 20)));
    engine.add_method(Box::new(WindowConsensusAssembly::new(50, 10)));

    let results = engine.evaluate_all(&reads, &ref_seq, &region).unwrap();
    assert_eq!(results.len(), 3);

    // Results should be sorted by fitness (best first)
    for i in 1..results.len() {
        assert!(
            results[i - 1].fitness.overall >= results[i].fitness.overall,
            "Results not sorted: {} < {}",
            results[i - 1].fitness.overall,
            results[i].fitness.overall
        );
    }
}

#[test]
fn test_engine_run_specific_method_with_test_data() {
    let (reads, ref_seq) = load_full();

    let mut engine = AssemblyEngine::new();
    engine.add_method(Box::new(ConsensusAssembly));

    let result = engine.run_method("consensus", &reads, &ref_seq).unwrap();
    assert!(result.is_some());

    let assembly = result.unwrap();
    assert_eq!(assembly.method_name, "consensus");
    assert_eq!(assembly.sequence.len(), ref_seq.len());
}

#[test]
fn test_engine_nonexistent_method() {
    let (reads, ref_seq) = load_full();

    let mut engine = AssemblyEngine::new();
    engine.add_method(Box::new(ConsensusAssembly));

    let result = engine.run_method("nonexistent", &reads, &ref_seq).unwrap();
    assert!(result.is_none());
}

// ===========================================================================
// 7. Haplotype assignment tests
// ===========================================================================

#[test]
fn test_haplotype_assignment_from_tags() {
    let (reads, ref_seq) = load_full();
    let assigner = HaplotypeAssigner::new();
    let assignments = assigner.assign(&reads, &ref_seq);

    assert_eq!(
        assignments.len(),
        reads.len(),
        "Should have one assignment per read"
    );
}

#[test]
fn test_haplotype_assignments_contain_both_haplotypes() {
    let (reads, ref_seq) = load_full();
    let assigner = HaplotypeAssigner::new();
    let assignments = assigner.assign(&reads, &ref_seq);

    let hap1 = assignments
        .iter()
        .filter(|a| a.haplotype == HaplotypeLabel::Hap1)
        .count();
    let hap2 = assignments
        .iter()
        .filter(|a| a.haplotype == HaplotypeLabel::Hap2)
        .count();

    assert!(hap1 > 0, "No Hap1 assignments");
    assert!(hap2 > 0, "No Hap2 assignments");
}

#[test]
fn test_haplotype_assignments_have_names() {
    let (reads, ref_seq) = load_full();
    let assigner = HaplotypeAssigner::new();
    let assignments = assigner.assign(&reads, &ref_seq);

    for assignment in &assignments {
        assert!(!assignment.read_name.is_empty());
    }
}

#[test]
fn test_haplotype_assignment_confidence_in_range() {
    let (reads, ref_seq) = load_full();
    let assigner = HaplotypeAssigner::new();
    let assignments = assigner.assign(&reads, &ref_seq);

    for assignment in &assignments {
        assert!(
            (0.0..=1.0).contains(&assignment.confidence),
            "Assignment confidence out of range: {}",
            assignment.confidence
        );
    }
}

#[test]
fn test_haplotype_assignment_matches_hp_tags() {
    let (reads, ref_seq) = load_full();
    let assigner = HaplotypeAssigner::new();
    let assignments = assigner.assign(&reads, &ref_seq);

    // Since our data has HP tags, assignments should use them
    for (read, assignment) in reads.iter().zip(assignments.iter()) {
        match read.haplotype_tag {
            Some(1) => assert_eq!(
                assignment.haplotype,
                HaplotypeLabel::Hap1,
                "Read {} with HP:1 assigned to {:?}",
                read.name,
                assignment.haplotype
            ),
            Some(2) => assert_eq!(
                assignment.haplotype,
                HaplotypeLabel::Hap2,
                "Read {} with HP:2 assigned to {:?}",
                read.name,
                assignment.haplotype
            ),
            _ => {} // untagged reads can go either way
        }
    }
}

#[test]
fn test_haplotype_label_display() {
    assert_eq!(format!("{}", HaplotypeLabel::Hap1), "H1");
    assert_eq!(format!("{}", HaplotypeLabel::Hap2), "H2");
    assert_eq!(format!("{}", HaplotypeLabel::Unassigned), "U");
}

// ===========================================================================
// 8. Metrics / fitness tests
// ===========================================================================

#[test]
fn test_fitness_score_with_test_data() {
    let (reads, ref_seq) = load_full();
    let method = ConsensusAssembly;
    let assembly = method.assemble(&reads, &ref_seq).unwrap();

    let calc = MetricsCalculator::new();
    let fitness = calc.compute_fitness(&assembly, &reads, &ref_seq);

    assert!(
        fitness.overall > 0.0 && fitness.overall <= 1.0,
        "Fitness score out of range: {}",
        fitness.overall
    );
}

#[test]
fn test_fitness_mean_agreement() {
    let (reads, ref_seq) = load_full();
    let method = ConsensusAssembly;
    let assembly = method.assemble(&reads, &ref_seq).unwrap();

    let calc = MetricsCalculator::new();
    let fitness = calc.compute_fitness(&assembly, &reads, &ref_seq);

    assert!(
        fitness.mean_agreement > 0.5,
        "Mean agreement too low: {}",
        fitness.mean_agreement
    );
    assert!(
        fitness.mean_agreement <= 1.0,
        "Mean agreement above 1.0: {}",
        fitness.mean_agreement
    );
}

#[test]
fn test_fitness_mean_depth() {
    let (reads, ref_seq) = load_full();
    let method = ConsensusAssembly;
    let assembly = method.assemble(&reads, &ref_seq).unwrap();

    let calc = MetricsCalculator::new();
    let fitness = calc.compute_fitness(&assembly, &reads, &ref_seq);

    assert!(
        fitness.mean_depth > 1.0,
        "Mean depth too low: {}",
        fitness.mean_depth
    );
}

#[test]
fn test_fitness_reference_identity() {
    let (reads, ref_seq) = load_full();
    let method = ConsensusAssembly;
    let assembly = method.assemble(&reads, &ref_seq).unwrap();

    let calc = MetricsCalculator::new();
    let fitness = calc.compute_fitness(&assembly, &reads, &ref_seq);

    assert!(
        fitness.reference_identity > 0.5,
        "Reference identity too low: {}",
        fitness.reference_identity
    );
    assert!(
        fitness.reference_identity <= 1.0,
        "Reference identity above 1.0: {}",
        fitness.reference_identity
    );
}

#[test]
fn test_fitness_base_metrics_length() {
    let (reads, ref_seq) = load_full();
    let method = ConsensusAssembly;
    let assembly = method.assemble(&reads, &ref_seq).unwrap();

    let calc = MetricsCalculator::new();
    let fitness = calc.compute_fitness(&assembly, &reads, &ref_seq);

    assert_eq!(
        fitness.base_metrics.len(),
        ref_seq.len(),
        "Base metrics length should match reference"
    );
}

#[test]
fn test_fitness_base_metrics_agreement_in_range() {
    let (reads, ref_seq) = load_full();
    let method = ConsensusAssembly;
    let assembly = method.assemble(&reads, &ref_seq).unwrap();

    let calc = MetricsCalculator::new();
    let fitness = calc.compute_fitness(&assembly, &reads, &ref_seq);

    for metric in &fitness.base_metrics {
        assert!(
            (0.0..=1.0).contains(&metric.agreement),
            "Base metric agreement out of range at pos {}: {}",
            metric.position,
            metric.agreement
        );
    }
}

#[test]
fn test_fitness_base_metrics_detect_variants() {
    let (reads, ref_seq) = load_full();
    let method = ConsensusAssembly;
    let assembly = method.assemble(&reads, &ref_seq).unwrap();

    let calc = MetricsCalculator::new();
    let fitness = calc.compute_fitness(&assembly, &reads, &ref_seq);

    let variant_count = fitness.base_metrics.iter().filter(|m| m.is_variant).count();
    // With our test data, some variants should be detected (het sites)
    // But most positions should match the reference
    let total = fitness.base_metrics.len();
    assert!(
        variant_count < total / 2,
        "Too many variants: {}/{}",
        variant_count,
        total
    );
}

#[test]
fn test_fitness_custom_weights() {
    let (reads, ref_seq) = load_full();
    let method = ConsensusAssembly;
    let assembly = method.assemble(&reads, &ref_seq).unwrap();

    let calc_agreement = MetricsCalculator::with_weights(1.0, 0.0, 0.0);
    let calc_depth = MetricsCalculator::with_weights(0.0, 1.0, 0.0);
    let calc_quality = MetricsCalculator::with_weights(0.0, 0.0, 1.0);

    let f_agreement = calc_agreement.compute_fitness(&assembly, &reads, &ref_seq);
    let f_depth = calc_depth.compute_fitness(&assembly, &reads, &ref_seq);
    let f_quality = calc_quality.compute_fitness(&assembly, &reads, &ref_seq);

    // Each should produce a valid score
    assert!((0.0..=1.0).contains(&f_agreement.overall));
    assert!((0.0..=1.0).contains(&f_depth.overall));
    assert!((0.0..=1.0).contains(&f_quality.overall));

    // With 100% agreement weight and high agreement data, score should be high
    assert!(
        f_agreement.overall > 0.5,
        "Agreement-only fitness too low: {}",
        f_agreement.overall
    );
}

// ===========================================================================
// 9. Full pipeline integration tests
// ===========================================================================

#[test]
fn test_full_pipeline_assemble_and_evaluate() {
    let (reads, ref_seq) = load_full();
    let region = full_region();

    // 1. Assembly
    let method = ConsensusAssembly;
    let assembly = method.assemble(&reads, &ref_seq).unwrap();
    assert!(!assembly.sequence.is_empty());

    // 2. Metrics
    let calc = MetricsCalculator::new();
    let fitness = calc.compute_fitness(&assembly, &reads, &ref_seq);
    assert!(fitness.overall > 0.0);

    // 3. Haplotyping
    let assigner = HaplotypeAssigner::new();
    let assignments = assigner.assign(&reads, &ref_seq);
    assert_eq!(assignments.len(), reads.len());

    // 4. Engine evaluation
    let mut engine = AssemblyEngine::new();
    engine.add_method(Box::new(ConsensusAssembly));
    engine.add_method(Box::new(WindowConsensusAssembly::new(100, 20)));
    let results = engine.evaluate_all(&reads, &ref_seq, &region).unwrap();
    assert!(!results.is_empty());
}

#[test]
fn test_full_pipeline_sub_region() {
    let (reads, ref_seq) = load_sub();

    // Assembly should still work with subset of reads
    let method = ConsensusAssembly;
    let assembly = method.assemble(&reads, &ref_seq).unwrap();
    assert_eq!(assembly.sequence.len(), ref_seq.len());

    // Metrics
    let calc = MetricsCalculator::new();
    let fitness = calc.compute_fitness(&assembly, &reads, &ref_seq);
    assert!(fitness.overall > 0.0);
}

#[test]
fn test_different_assembly_methods_produce_different_results() {
    let (reads, ref_seq) = load_full();

    let consensus = ConsensusAssembly.assemble(&reads, &ref_seq).unwrap();
    let window = WindowConsensusAssembly::new(50, 10)
        .assemble(&reads, &ref_seq)
        .unwrap();

    // Both should produce same-length output
    assert_eq!(consensus.sequence.len(), window.sequence.len());

    // Method names should differ
    assert_ne!(consensus.method_name, window.method_name);
}

// ===========================================================================
// 10. App state tests
// ===========================================================================

#[test]
fn test_app_creation_with_test_data() {
    use longread_reviewer::viewer::App;

    let (reads, ref_seq) = load_full();
    let region = full_region();
    let app = App::new(region, ref_seq.clone(), reads.clone());

    assert_eq!(app.reads.len(), reads.len());
    assert_eq!(app.reference.len(), ref_seq.len());
    assert!(!app.should_quit);
    assert!(app.assembly.is_none());
}

#[test]
fn test_app_run_assembly_with_test_data() {
    use longread_reviewer::viewer::App;

    let (reads, ref_seq) = load_full();
    let region = full_region();
    let mut app = App::new(region, ref_seq, reads);

    app.run_assembly().unwrap();
    assert!(app.assembly.is_some());
    assert!(app.fitness.is_some());
}

#[test]
fn test_app_assign_haplotypes_with_test_data() {
    use longread_reviewer::viewer::App;

    let (reads, ref_seq) = load_full();
    let num_reads = reads.len();
    let region = full_region();
    let mut app = App::new(region, ref_seq, reads);

    app.assign_haplotypes();
    assert_eq!(app.haplotype_assignments.len(), num_reads);
}

#[test]
fn test_app_next_method_cycles() {
    use longread_reviewer::viewer::App;

    let (reads, ref_seq) = load_full();
    let region = full_region();
    let mut app = App::new(region, ref_seq, reads);

    let initial_method = app.current_method;
    app.next_method().unwrap();
    assert_ne!(app.current_method, initial_method);

    // Cycle back to start
    let num_methods = app.engine.method_names().len();
    for _ in 1..num_methods {
        app.next_method().unwrap();
    }
    assert_eq!(app.current_method, initial_method);
}

#[test]
fn test_app_evaluate_all_with_test_data() {
    use longread_reviewer::viewer::App;

    let (reads, ref_seq) = load_full();
    let region = full_region();
    let app = App::new(region, ref_seq, reads);

    let results = app.evaluate_all().unwrap();
    assert!(!results.is_empty());

    for (name, score) in &results {
        assert!(!name.is_empty());
        assert!(
            *score >= 0.0 && *score <= 1.0,
            "Score out of range for {}: {}",
            name, score
        );
    }
}

// ===========================================================================
// 11. Region parsing edge-case tests
// ===========================================================================

#[test]
fn test_region_chr17_parse() {
    let r: Region = "chr17:10958130-10971414".parse().unwrap();
    assert_eq!(r.chrom, "chr17");
    assert_eq!(r.start, 10958130);
    assert_eq!(r.end, 10971414);
}

#[test]
fn test_region_length() {
    let r = full_region();
    assert_eq!(r.len(), 2001);
}

#[test]
fn test_region_display_roundtrip() {
    let r = full_region();
    let s = r.to_string();
    let r2: Region = s.parse().unwrap();
    assert_eq!(r, r2);
}
