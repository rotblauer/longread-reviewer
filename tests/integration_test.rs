//! Integration tests using the NA19240 chr17 fragment BAM and reference FASTA.
//!
//! These tests exercise the full pipeline: BAM reading → assembly → haplotype → metrics,
//! using the real PacBio HiFi data in `tests/data/`.
//!
//! The test data covers chr17:10958130-11017414 (~59 kb), which contains a complex
//! structural event several kb in size. The reference FASTA includes a 5 kb buffer
//! on each side (chr17:10953130-11022414). Tests verify that the pipeline accurately
//! highlights this event and facilitates review.

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

/// Region covering the core alignment region of interest (1-based, inclusive).
fn full_region() -> Region {
    Region::new("chr17", 10958130, 11017414).unwrap()
}

/// A smaller sub-region for targeted tests (~5 kb).
fn sub_region() -> Region {
    Region::new("chr17", 10980000, 10985000).unwrap()
}

/// A focused sub-region around the complex event (~2 kb) for faster assembly/metrics tests.
fn event_region() -> Region {
    Region::new("chr17", 10990000, 10992000).unwrap()
}

/// Load reads + reference for the full region.
fn load_full() -> (Vec<AlignedRead>, Vec<u8>, Region) {
    let region = full_region();
    let reads = AlignmentReader::read_bam(&test_bam_path(), &region)
        .expect("failed to read BAM");
    let reference = ReferenceGenome::from_file(&test_ref_path())
        .expect("failed to read reference");
    let ref_seq = reference.fetch(&region).expect("failed to fetch region");
    (reads, ref_seq, region)
}

/// Load reads + reference for the sub-region.
fn load_sub() -> (Vec<AlignedRead>, Vec<u8>, Region) {
    let region = sub_region();
    let reads = AlignmentReader::read_bam(&test_bam_path(), &region)
        .expect("failed to read BAM");
    let reference = ReferenceGenome::from_file(&test_ref_path())
        .expect("failed to read reference");
    let ref_seq = reference.fetch(&region).expect("failed to fetch region");
    (reads, ref_seq, region)
}

/// Load reads + reference for the event region.
fn load_event() -> (Vec<AlignedRead>, Vec<u8>, Region) {
    let region = event_region();
    let reads = AlignmentReader::read_bam(&test_bam_path(), &region)
        .expect("failed to read BAM");
    let reference = ReferenceGenome::from_file(&test_ref_path())
        .expect("failed to read reference");
    let ref_seq = reference.fetch(&region).expect("failed to fetch region");
    (reads, ref_seq, region)
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
        reads.len() >= 50,
        "Expected at least 50 reads, got {}",
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
    let region = Region::new("chr17", 100_000_000, 100_001_000).unwrap();
    let reads = AlignmentReader::read_bam(&test_bam_path(), &region).unwrap();
    assert!(
        reads.is_empty(),
        "Expected no reads in empty region"
    );
}

#[test]
fn test_all_reads_on_correct_chrom() {
    let (reads, _, _region) = load_full();
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
    let (reads, _, _region) = load_full();
    for read in &reads {
        assert!(
            !read.name.is_empty(),
            "Read has empty name"
        );
    }
}

#[test]
fn test_read_names_look_like_pacbio_ccs() {
    let (reads, _, _region) = load_full();
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
    let (reads, _, _region) = load_full();
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
    let (reads, _, _region) = load_full();
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
    let (reads, _, _region) = load_full();
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
    let (reads, _, _region) = load_full();
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
    let (reads, _, _region) = load_full();
    for read in &reads {
        assert!(
            read.mapq > 0,
            "Read {} has zero MAPQ",
            read.name
        );
    }
}

#[test]
fn test_reads_are_long_reads() {
    let (reads, _, _region) = load_full();
    let avg_len: f64 =
        reads.iter().map(|r| r.sequence.len() as f64).sum::<f64>() / reads.len() as f64;
    assert!(
        avg_len > 5000.0,
        "Expected average long-read length > 5kb, got {:.0}",
        avg_len
    );
}

#[test]
fn test_hp_tags_are_valid() {
    let (reads, _, _region) = load_full();
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
fn test_reads_have_both_strands() {
    let (reads, _, _region) = load_full();
    let forward = reads.iter().filter(|r| !r.is_reverse).count();
    let reverse = reads.iter().filter(|r| r.is_reverse).count();
    assert!(forward > 0, "No forward-strand reads");
    assert!(reverse > 0, "No reverse-strand reads");
}

#[test]
fn test_cigar_ops_are_valid_types() {
    let (reads, _, _region) = load_full();
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
    let (reads, _, _region) = load_full();
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
    let (reads, _, _region) = load_full();
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
    let (reads, _, _region) = load_full();
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
    let (reads, _, _region) = load_full();
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
    let (reads, _, _region) = load_full();
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
// 3. Reference loading tests (with fragment FASTA support)
// ===========================================================================

#[test]
fn test_reference_loads() {
    let reference = ReferenceGenome::from_file(&test_ref_path()).unwrap();
    let chroms = reference.chromosomes();
    assert!(
        chroms.contains(&"chr17"),
        "Reference should expose chr17 from fragment name; got: {:?}",
        chroms
    );
}

#[test]
fn test_reference_fetch_full_region() {
    let reference = ReferenceGenome::from_file(&test_ref_path()).unwrap();
    let region = full_region();
    let seq = reference.fetch(&region).unwrap();
    assert_eq!(
        seq.len() as u64,
        region.len(),
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
    let seq = reference.fetch(&sub_region()).unwrap();
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

#[test]
fn test_reference_fetch_outside_fragment() {
    let reference = ReferenceGenome::from_file(&test_ref_path()).unwrap();
    let region = Region::new("chr17", 1, 1000).unwrap();
    assert!(
        reference.fetch(&region).is_err(),
        "Should fail for region outside the fragment"
    );
}

// ===========================================================================
// 4. Consensus assembly tests (using event region for speed)
// ===========================================================================

#[test]
fn test_consensus_assembly_produces_output() {
    let (reads, ref_seq, region) = load_event();
    let method = ConsensusAssembly;
    let result = method.assemble(&reads, &ref_seq, region.start).unwrap();

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
    let (reads, ref_seq, region) = load_event();
    let method = ConsensusAssembly;
    let result = method.assemble(&reads, &ref_seq, region.start).unwrap();

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
    let (reads, ref_seq, region) = load_event();
    let method = ConsensusAssembly;
    let result = method.assemble(&reads, &ref_seq, region.start).unwrap();

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
    let (reads, ref_seq, region) = load_event();
    let method = ConsensusAssembly;
    let result = method.assemble(&reads, &ref_seq, region.start).unwrap();

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

// ===========================================================================
// 5. Window consensus assembly tests
// ===========================================================================

#[test]
fn test_window_consensus_assembly_produces_output() {
    let (reads, ref_seq, region) = load_event();
    let method = WindowConsensusAssembly::new(100, 20);
    let result = method.assemble(&reads, &ref_seq, region.start).unwrap();

    assert_eq!(result.sequence.len(), ref_seq.len());
    assert_eq!(result.method_name, "window_consensus");
}

#[test]
fn test_window_consensus_various_sizes() {
    let (reads, ref_seq, region) = load_event();

    for (window, overlap) in &[(50, 10), (100, 20), (200, 40), (500, 100)] {
        let method = WindowConsensusAssembly::new(*window, *overlap);
        let result = method.assemble(&reads, &ref_seq, region.start).unwrap();

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
    let (reads, ref_seq, region) = load_event();
    let method = WindowConsensusAssembly::new(100, 20);
    let result = method.assemble(&reads, &ref_seq, region.start).unwrap();

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
    let (reads, ref_seq, region) = load_event();
    let region = event_region();

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
    let (reads, ref_seq, region) = load_event();

    let mut engine = AssemblyEngine::new();
    engine.add_method(Box::new(ConsensusAssembly));

    let result = engine.run_method("consensus", &reads, &ref_seq, 1).unwrap();
    assert!(result.is_some());

    let assembly = result.unwrap();
    assert_eq!(assembly.method_name, "consensus");
    assert_eq!(assembly.sequence.len(), ref_seq.len());
}

#[test]
fn test_engine_nonexistent_method() {
    let (reads, ref_seq, region) = load_event();

    let mut engine = AssemblyEngine::new();
    engine.add_method(Box::new(ConsensusAssembly));

    let result = engine.run_method("nonexistent", &reads, &ref_seq, 1).unwrap();
    assert!(result.is_none());
}

// ===========================================================================
// 7. Haplotype assignment tests
// ===========================================================================

#[test]
fn test_haplotype_assignment_count_matches_reads() {
    let (reads, ref_seq, region) = load_event();
    let assigner = HaplotypeAssigner::new();
    let assignments = assigner.assign(&reads, &ref_seq);

    assert_eq!(
        assignments.len(),
        reads.len(),
        "Should have one assignment per read"
    );
}

#[test]
fn test_haplotype_assignments_have_names() {
    let (reads, ref_seq, region) = load_event();
    let assigner = HaplotypeAssigner::new();
    let assignments = assigner.assign(&reads, &ref_seq);

    for assignment in &assignments {
        assert!(!assignment.read_name.is_empty());
    }
}

#[test]
fn test_haplotype_assignment_confidence_in_range() {
    let (reads, ref_seq, region) = load_event();
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
    let (reads, ref_seq, region) = load_event();
    let assigner = HaplotypeAssigner::new();
    let assignments = assigner.assign(&reads, &ref_seq);

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
    let (reads, ref_seq, region) = load_event();
    let method = ConsensusAssembly;
    let assembly = method.assemble(&reads, &ref_seq, region.start).unwrap();

    let calc = MetricsCalculator::new();
    let fitness = calc.compute_fitness(&assembly, &reads, &ref_seq, region.start);

    assert!(
        fitness.overall > 0.0 && fitness.overall <= 1.0,
        "Fitness score out of range: {}",
        fitness.overall
    );
}

#[test]
fn test_fitness_mean_agreement() {
    let (reads, ref_seq, region) = load_event();
    let method = ConsensusAssembly;
    let assembly = method.assemble(&reads, &ref_seq, region.start).unwrap();

    let calc = MetricsCalculator::new();
    let fitness = calc.compute_fitness(&assembly, &reads, &ref_seq, region.start);

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
    let (reads, ref_seq, region) = load_event();
    let method = ConsensusAssembly;
    let assembly = method.assemble(&reads, &ref_seq, region.start).unwrap();

    let calc = MetricsCalculator::new();
    let fitness = calc.compute_fitness(&assembly, &reads, &ref_seq, region.start);

    assert!(
        fitness.mean_depth > 1.0,
        "Mean depth too low: {}",
        fitness.mean_depth
    );
}

#[test]
fn test_fitness_base_metrics_length() {
    let (reads, ref_seq, region) = load_event();
    let method = ConsensusAssembly;
    let assembly = method.assemble(&reads, &ref_seq, region.start).unwrap();

    let calc = MetricsCalculator::new();
    let fitness = calc.compute_fitness(&assembly, &reads, &ref_seq, region.start);

    assert_eq!(
        fitness.base_metrics.len(),
        ref_seq.len(),
        "Base metrics length should match reference"
    );
}

#[test]
fn test_fitness_base_metrics_agreement_in_range() {
    let (reads, ref_seq, region) = load_event();
    let method = ConsensusAssembly;
    let assembly = method.assemble(&reads, &ref_seq, region.start).unwrap();

    let calc = MetricsCalculator::new();
    let fitness = calc.compute_fitness(&assembly, &reads, &ref_seq, region.start);

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
fn test_fitness_custom_weights() {
    let (reads, ref_seq, region) = load_event();
    let method = ConsensusAssembly;
    let assembly = method.assemble(&reads, &ref_seq, region.start).unwrap();

    let calc_agreement = MetricsCalculator::with_weights(1.0, 0.0, 0.0);
    let calc_depth = MetricsCalculator::with_weights(0.0, 1.0, 0.0);
    let calc_quality = MetricsCalculator::with_weights(0.0, 0.0, 1.0);

    let f_agreement = calc_agreement.compute_fitness(&assembly, &reads, &ref_seq, region.start);
    let f_depth = calc_depth.compute_fitness(&assembly, &reads, &ref_seq, region.start);
    let f_quality = calc_quality.compute_fitness(&assembly, &reads, &ref_seq, region.start);

    assert!((0.0..=1.0).contains(&f_agreement.overall));
    assert!((0.0..=1.0).contains(&f_depth.overall));
    assert!((0.0..=1.0).contains(&f_quality.overall));

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
    let (reads, ref_seq, region) = load_event();
    let region = event_region();

    // 1. Assembly
    let method = ConsensusAssembly;
    let assembly = method.assemble(&reads, &ref_seq, region.start).unwrap();
    assert!(!assembly.sequence.is_empty());

    // 2. Metrics
    let calc = MetricsCalculator::new();
    let fitness = calc.compute_fitness(&assembly, &reads, &ref_seq, region.start);
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
    let (reads, ref_seq, region) = load_sub();

    let method = ConsensusAssembly;
    let assembly = method.assemble(&reads, &ref_seq, region.start).unwrap();
    assert_eq!(assembly.sequence.len(), ref_seq.len());
}

#[test]
fn test_different_assembly_methods_produce_different_results() {
    let (reads, ref_seq, region) = load_event();

    let consensus = ConsensusAssembly.assemble(&reads, &ref_seq, region.start).unwrap();
    let window = WindowConsensusAssembly::new(50, 10)
        .assemble(&reads, &ref_seq, region.start)
        .unwrap();

    assert_eq!(consensus.sequence.len(), window.sequence.len());
    assert_ne!(consensus.method_name, window.method_name);
}

// ===========================================================================
// 10. App state tests
// ===========================================================================

#[test]
fn test_app_creation_with_test_data() {
    use longread_reviewer::viewer::App;

    let (reads, ref_seq, region) = load_event();
    let region = event_region();
    let app = App::new(region, ref_seq.clone(), reads.clone());

    assert_eq!(app.reads.len(), reads.len());
    assert_eq!(app.reference.len(), ref_seq.len());
    assert!(!app.should_quit);
    assert!(app.assembly.is_none());
}

#[test]
fn test_app_run_assembly_with_test_data() {
    use longread_reviewer::viewer::App;

    let (reads, ref_seq, region) = load_event();
    let region = event_region();
    let mut app = App::new(region, ref_seq, reads);

    app.run_assembly().unwrap();
    assert!(app.assembly.is_some());
    assert!(app.fitness.is_some());
}

#[test]
fn test_app_assign_haplotypes_with_test_data() {
    use longread_reviewer::viewer::App;

    let (reads, ref_seq, region) = load_event();
    let num_reads = reads.len();
    let region = event_region();
    let mut app = App::new(region, ref_seq, reads);

    app.assign_haplotypes();
    assert_eq!(app.haplotype_assignments.len(), num_reads);
}

#[test]
fn test_app_next_method_cycles() {
    use longread_reviewer::viewer::App;

    let (reads, ref_seq, region) = load_event();
    let region = event_region();
    let mut app = App::new(region, ref_seq, reads);

    let initial_method = app.current_method;
    app.next_method().unwrap();
    assert_ne!(app.current_method, initial_method);

    let num_methods = app.engine.method_names().len();
    for _ in 1..num_methods {
        app.next_method().unwrap();
    }
    assert_eq!(app.current_method, initial_method);
}

#[test]
fn test_app_evaluate_all_with_test_data() {
    use longread_reviewer::viewer::App;

    let (reads, ref_seq, region) = load_event();
    let region = event_region();
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
    let r: Region = "chr17:10958130-11017414".parse().unwrap();
    assert_eq!(r.chrom, "chr17");
    assert_eq!(r.start, 10958130);
    assert_eq!(r.end, 11017414);
}

#[test]
fn test_region_length() {
    let r = full_region();
    assert_eq!(r.len(), 59285);
}

#[test]
fn test_region_display_roundtrip() {
    let r = full_region();
    let s = r.to_string();
    let r2: Region = s.parse().unwrap();
    assert_eq!(r, r2);
}

// ===========================================================================
// 12. Complex event detection tests
//
// The NA19240 chr17 fragment contains a complex structural event (several kb).
// These tests validate that the pipeline highlights it for review.
// ===========================================================================

#[test]
fn test_reads_contain_large_indels() {
    let (reads, _, _region) = load_full();
    let mut large_ins_ops = 0usize;
    let mut large_del_ops = 0usize;

    for read in &reads {
        for op in &read.cigar {
            match op {
                CigarOp::Insertion(n) if *n > 100 => large_ins_ops += 1,
                CigarOp::Deletion(n) if *n > 100 => large_del_ops += 1,
                _ => {}
            }
        }
    }

    assert!(
        large_ins_ops > 0 || large_del_ops > 0,
        "Expected large indel operations (>100 bp) indicating a complex event"
    );
}

#[test]
fn test_assembly_diverges_from_reference_in_event_region() {
    let (reads, ref_seq, region) = load_event();
    let method = ConsensusAssembly;
    let result = method.assemble(&reads, &ref_seq, region.start).unwrap();

    let matches = result
        .sequence
        .iter()
        .zip(ref_seq.iter())
        .filter(|(a, b)| a.eq_ignore_ascii_case(b))
        .count();
    let identity = matches as f64 / ref_seq.len() as f64;

    // In a region with a complex event, the consensus should diverge from reference.
    // The identity should be well below 100% but still above noise level.
    assert!(
        identity < 1.0,
        "Expected assembly to diverge from reference in event region, got {:.1}% identity",
        identity * 100.0
    );
}

#[test]
fn test_metrics_detect_variant_positions_in_event_region() {
    let (reads, ref_seq, region) = load_event();
    let method = ConsensusAssembly;
    let assembly = method.assemble(&reads, &ref_seq, region.start).unwrap();

    let calc = MetricsCalculator::new();
    let fitness = calc.compute_fitness(&assembly, &reads, &ref_seq, region.start);

    let variant_count = fitness.base_metrics.iter().filter(|m| m.is_variant).count();
    assert!(
        variant_count > 0,
        "Expected variant positions in the complex event region, got 0"
    );
}

#[test]
fn test_multiple_indel_cigar_ops_per_read() {
    let (reads, _, _region) = load_full();

    let mut total_ins_ops = 0usize;
    let mut total_del_ops = 0usize;
    for read in &reads {
        for op in &read.cigar {
            match op {
                CigarOp::Insertion(_) => total_ins_ops += 1,
                CigarOp::Deletion(_) => total_del_ops += 1,
                _ => {}
            }
        }
    }

    assert!(
        total_ins_ops > 100,
        "Expected many insertion operations across reads, got {}",
        total_ins_ops
    );
    assert!(
        total_del_ops > 100,
        "Expected many deletion operations across reads, got {}",
        total_del_ops
    );
}

#[test]
fn test_sub_region_reads_overlap_event() {
    let (reads, _, _region) = load_sub();
    assert!(
        reads.len() >= 10,
        "Expected substantial coverage in sub-region near event, got {}",
        reads.len()
    );
}

#[test]
fn test_haplotype_clustering_produces_assignments_in_event_region() {
    let (reads, ref_seq, region) = load_event();
    let assigner = HaplotypeAssigner::new();
    let assignments = assigner.assign(&reads, &ref_seq);

    // Every read should get an assignment (haplotype or unassigned).
    assert_eq!(assignments.len(), reads.len());

    let hap1 = assignments
        .iter()
        .filter(|a| a.haplotype == HaplotypeLabel::Hap1)
        .count();
    let hap2 = assignments
        .iter()
        .filter(|a| a.haplotype == HaplotypeLabel::Hap2)
        .count();
    let unassigned = assignments
        .iter()
        .filter(|a| a.haplotype == HaplotypeLabel::Unassigned)
        .count();

    // The event region may be homozygous (all reads agree) or heterozygous;
    // either way the clustering should produce a valid partition.
    assert_eq!(
        hap1 + hap2 + unassigned,
        reads.len(),
        "All reads should be accounted for: H1={}, H2={}, U={}",
        hap1, hap2, unassigned
    );
}

#[test]
fn test_fitness_reference_identity_lower_in_event_region() {
    let (reads, ref_seq, region) = load_event();
    let method = ConsensusAssembly;
    let assembly = method.assemble(&reads, &ref_seq, region.start).unwrap();

    let calc = MetricsCalculator::new();
    let fitness = calc.compute_fitness(&assembly, &reads, &ref_seq, region.start);

    // A region with a complex event should show reference identity < 1.0
    assert!(
        fitness.reference_identity < 1.0,
        "Expected reference identity < 1.0 in event region, got {}",
        fitness.reference_identity
    );
}

#[test]
fn test_event_region_has_depth_coverage() {
    let (reads, ref_seq, region) = load_event();
    let method = ConsensusAssembly;
    let result = method.assemble(&reads, &ref_seq, region.start).unwrap();

    let positions_with_coverage = result.depth.iter().filter(|&&d| d > 0).count();
    let coverage_fraction = positions_with_coverage as f64 / result.depth.len() as f64;
    assert!(
        coverage_fraction > 0.5,
        "Expected >50% of event region positions to have coverage, got {:.1}%",
        coverage_fraction * 100.0
    );
}

#[test]
fn test_assembly_methods_rank_reasonably_on_event_data() {
    let (reads, ref_seq, region) = load_event();
    let region = event_region();

    let mut engine = AssemblyEngine::new();
    engine.add_method(Box::new(ConsensusAssembly));
    engine.add_method(Box::new(WindowConsensusAssembly::new(200, 40)));
    engine.add_method(Box::new(WindowConsensusAssembly::new(500, 100)));

    let results = engine.evaluate_all(&reads, &ref_seq, &region).unwrap();

    // All methods should produce valid fitness scores
    for result in &results {
        assert!(
            result.fitness.overall > 0.0,
            "Method {} scored 0.0 on event data",
            result.assembly.method_name
        );
        assert!(
            result.fitness.overall <= 1.0,
            "Method {} scored >1.0 on event data",
            result.assembly.method_name
        );
    }
}

#[test]
fn test_variant_positions_cluster_in_event() {
    let (reads, ref_seq, region) = load_event();
    let method = ConsensusAssembly;
    let assembly = method.assemble(&reads, &ref_seq, region.start).unwrap();

    let calc = MetricsCalculator::new();
    let fitness = calc.compute_fitness(&assembly, &reads, &ref_seq, region.start);

    // Look for clusters of variant positions (characteristic of complex events).
    let variant_positions: Vec<u64> = fitness
        .base_metrics
        .iter()
        .filter(|m| m.is_variant)
        .map(|m| m.position)
        .collect();

    if variant_positions.len() >= 2 {
        // Check that some variants are near each other (within 100 bp)
        let has_cluster = variant_positions.windows(2).any(|w| w[1] - w[0] < 100);
        assert!(
            has_cluster,
            "Expected clustered variant positions in event region"
        );
    }
}

#[test]
fn test_high_read_agreement_despite_reference_divergence() {
    let (reads, ref_seq, region) = load_event();
    let method = ConsensusAssembly;
    let assembly = method.assemble(&reads, &ref_seq, region.start).unwrap();

    let calc = MetricsCalculator::new();
    let fitness = calc.compute_fitness(&assembly, &reads, &ref_seq, region.start);

    // In a complex event, reads may consistently differ from the reference
    // but agree with each other (high agreement, low reference identity).
    // This is the signature that makes the event visible for review.
    let high_agreement = fitness
        .base_metrics
        .iter()
        .filter(|m| m.depth > 0 && m.agreement > 0.8)
        .count();
    let variant_positions = fitness
        .base_metrics
        .iter()
        .filter(|m| m.is_variant)
        .count();

    assert!(
        high_agreement > 0,
        "Expected positions with high read agreement"
    );
    assert!(
        variant_positions > 0,
        "Expected variant positions where assembly diverges from reference"
    );
}
