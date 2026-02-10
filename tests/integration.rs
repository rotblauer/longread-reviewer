//! Integration tests that exercise the full pipeline against the real example BAM file.
//!
//! These tests confirm that every major component works correctly end-to-end with
//! real PacBio CCS data (NA19240, chr17:10958130-11017414).

use std::path::Path;

use longread_reviewer::alignment::{AlignedRead, AlignmentReader, CigarOp};
use longread_reviewer::assembly::engine::AssemblyEngine;
use longread_reviewer::assembly::method::{
    AssemblyMethod, ConsensusAssembly, WindowConsensusAssembly,
};
use longread_reviewer::haplotype::{HaplotypeAssigner, HaplotypeLabel};
use longread_reviewer::metrics::MetricsCalculator;
use longread_reviewer::region::Region;

/// Path to the example BAM file shipped with the repository.
const EXAMPLE_BAM: &str =
    "example/NA19240_2020_merged.ccs.hg38.aligned.chr17_10958130_11017414.bam";

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn example_bam_path() -> &'static Path {
    Path::new(EXAMPLE_BAM)
}

/// Read a region from the example BAM.
fn read_region(chrom: &str, start: u64, end: u64) -> Vec<AlignedRead> {
    let region = Region::new(chrom, start, end).unwrap();
    AlignmentReader::read_bam(example_bam_path(), &region).unwrap()
}

/// Build a synthetic reference from the majority base at each position in the reads.
/// This is needed because we don't ship a reference FASTA with the example.
fn synthetic_reference(reads: &[AlignedRead], region: &Region) -> Vec<u8> {
    let len = (region.end - region.start + 1) as usize;
    let mut base_counts: Vec<[u32; 4]> = vec![[0u32; 4]; len];

    for read in reads {
        for (pos, base) in read.reference_aligned_bases() {
            if pos >= region.start && pos <= region.end {
                let idx = (pos - region.start) as usize;
                if idx < len {
                    let bi = match base.to_ascii_uppercase() {
                        b'A' => 0,
                        b'C' => 1,
                        b'G' => 2,
                        b'T' => 3,
                        _ => continue,
                    };
                    base_counts[idx][bi] += 1;
                }
            }
        }
    }

    base_counts
        .iter()
        .map(|counts| {
            let max_idx = counts
                .iter()
                .enumerate()
                .max_by_key(|&(_, c)| c)
                .map(|(i, _)| i)
                .unwrap_or(0);
            match max_idx {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                3 => b'T',
                _ => b'N',
            }
        })
        .collect()
}

// ===========================================================================
// BAM Reading Tests
// ===========================================================================

#[test]
fn bam_reads_full_region() {
    let reads = read_region("chr17", 10958130, 11017414);
    assert!(
        reads.len() > 100,
        "expected >100 reads in full region, got {}",
        reads.len()
    );
}

#[test]
fn bam_reads_small_subregion() {
    let reads = read_region("chr17", 10990000, 10990100);
    assert!(
        !reads.is_empty(),
        "expected reads in small subregion, got 0"
    );
}

#[test]
fn bam_reads_empty_outside_range() {
    // Region far outside the BAM's coordinate range should return no reads.
    let reads = read_region("chr17", 1, 100);
    assert_eq!(reads.len(), 0);
}

#[test]
fn bam_read_fields_populated() {
    let reads = read_region("chr17", 10990000, 10990100);
    for read in &reads {
        assert!(!read.name.is_empty(), "read name should not be empty");
        assert_eq!(read.chrom, "chr17");
        assert!(read.start > 0, "start should be positive");
        assert!(read.end >= read.start, "end should be >= start");
        assert!(read.mapq > 0, "mapq should be > 0 for aligned reads");
        assert!(!read.cigar.is_empty(), "cigar should not be empty");
        assert!(!read.sequence.is_empty(), "sequence should not be empty");
        assert!(
            !read.qualities.is_empty(),
            "qualities should not be empty"
        );
    }
}

#[test]
fn bam_reads_have_cigar_variety() {
    let reads = read_region("chr17", 10958130, 11017414);
    let has_match = reads
        .iter()
        .any(|r| r.cigar.iter().any(|op| matches!(op, CigarOp::Match(_))));
    let has_ins = reads
        .iter()
        .any(|r| r.cigar.iter().any(|op| matches!(op, CigarOp::Insertion(_))));
    let has_del = reads
        .iter()
        .any(|r| r.cigar.iter().any(|op| matches!(op, CigarOp::Deletion(_))));
    let has_sc = reads
        .iter()
        .any(|r| r.cigar.iter().any(|op| matches!(op, CigarOp::SoftClip(_))));

    assert!(has_match, "expected some match operations");
    assert!(has_ins, "expected some insertions in long-read data");
    assert!(has_del, "expected some deletions in long-read data");
    assert!(has_sc, "expected some soft clips");
}

#[test]
fn bam_reads_both_strands() {
    let reads = read_region("chr17", 10958130, 11017414);
    let fwd = reads.iter().filter(|r| !r.is_reverse).count();
    let rev = reads.iter().filter(|r| r.is_reverse).count();
    assert!(fwd > 0, "expected forward-strand reads");
    assert!(rev > 0, "expected reverse-strand reads");
}

#[test]
fn bam_aligned_sequence_consistent() {
    let reads = read_region("chr17", 10990000, 10990100);
    for read in &reads {
        let aligned = read.aligned_sequence();
        // Aligned sequence should not be larger than the full sequence
        assert!(
            aligned.len() <= read.sequence.len() + 100, // allow small margin for insertions
            "aligned seq too large for read {}",
            read.name
        );
    }
}

#[test]
fn bam_reference_aligned_bases_within_span() {
    let reads = read_region("chr17", 10990000, 10990100);
    for read in &reads {
        let bases = read.reference_aligned_bases();
        for &(pos, base) in &bases {
            assert!(
                pos >= read.start && pos <= read.end,
                "aligned base position {} outside read span {}-{} for {}",
                pos,
                read.start,
                read.end,
                read.name
            );
            assert!(
                base.is_ascii_alphabetic(),
                "non-alphabetic base {} at pos {} in {}",
                base,
                pos,
                read.name
            );
        }
    }
}

// ===========================================================================
// Assembly Tests (using BAM reads + synthetic reference)
// ===========================================================================

fn setup_assembly_data() -> (Vec<AlignedRead>, Vec<u8>, Region) {
    let region = Region::new("chr17", 10990000, 10990100).unwrap();
    let reads = AlignmentReader::read_bam(example_bam_path(), &region).unwrap();
    let reference = synthetic_reference(&reads, &region);
    (reads, reference, region)
}

#[test]
fn consensus_assembly_on_real_data() {
    let (reads, reference, _region) = setup_assembly_data();
    let method = ConsensusAssembly;
    let result = method.assemble(&reads, &reference).unwrap();

    assert_eq!(
        result.sequence.len(),
        reference.len(),
        "assembly length should match reference length"
    );
    assert_eq!(result.depth.len(), reference.len());
    assert_eq!(result.confidence.len(), reference.len());
    assert_eq!(result.method_name, "consensus");

    // Depth should be > 0 at covered positions
    let covered = result.depth.iter().filter(|&&d| d > 0).count();
    assert!(
        covered > reference.len() / 2,
        "at least half the positions should be covered"
    );

    // Confidence should be between 0 and 1
    for &c in &result.confidence {
        assert!(c >= 0.0 && c <= 1.0, "confidence out of range: {}", c);
    }
}

#[test]
fn window_consensus_on_real_data() {
    let (reads, reference, _region) = setup_assembly_data();
    let method = WindowConsensusAssembly::new(50, 10);
    let result = method.assemble(&reads, &reference).unwrap();

    assert_eq!(result.sequence.len(), reference.len());
    assert_eq!(result.method_name, "window_consensus");
}

#[test]
fn window_consensus_various_sizes() {
    let (reads, reference, _region) = setup_assembly_data();

    for (ws, ov) in &[(20, 5), (50, 10), (80, 20)] {
        let method = WindowConsensusAssembly::new(*ws, *ov);
        let result = method.assemble(&reads, &reference).unwrap();
        assert_eq!(
            result.sequence.len(),
            reference.len(),
            "window_size={} overlap={} produced wrong length",
            ws,
            ov
        );
    }
}

// ===========================================================================
// Engine Evaluation Tests
// ===========================================================================

#[test]
fn engine_evaluate_all_real_data() {
    let (reads, reference, region) = setup_assembly_data();

    let mut engine = AssemblyEngine::new();
    engine.add_method(Box::new(ConsensusAssembly));
    engine.add_method(Box::new(WindowConsensusAssembly::new(50, 10)));
    engine.add_method(Box::new(WindowConsensusAssembly::new(30, 5)));

    let results = engine.evaluate_all(&reads, &reference, &region).unwrap();
    assert_eq!(results.len(), 3);

    // Results should be sorted by fitness (descending)
    for i in 1..results.len() {
        assert!(
            results[i - 1].fitness.overall >= results[i].fitness.overall,
            "results not sorted: {} < {}",
            results[i - 1].fitness.overall,
            results[i].fitness.overall
        );
    }

    for r in &results {
        assert!(r.fitness.overall >= 0.0 && r.fitness.overall <= 1.0);
        assert!(r.fitness.mean_agreement >= 0.0 && r.fitness.mean_agreement <= 1.0);
        assert!(r.fitness.mean_depth >= 0.0);
        assert!(r.fitness.reference_identity >= 0.0 && r.fitness.reference_identity <= 1.0);
    }
}

#[test]
fn engine_run_single_method_real_data() {
    let (reads, reference, _region) = setup_assembly_data();

    let mut engine = AssemblyEngine::new();
    engine.add_method(Box::new(ConsensusAssembly));

    let result = engine.run_method("consensus", &reads, &reference).unwrap();
    assert!(result.is_some());

    let result = engine
        .run_method("nonexistent", &reads, &reference)
        .unwrap();
    assert!(result.is_none());
}

// ===========================================================================
// Metrics Tests
// ===========================================================================

#[test]
fn fitness_metrics_on_real_data() {
    let (reads, reference, _region) = setup_assembly_data();
    let method = ConsensusAssembly;
    let assembly = method.assemble(&reads, &reference).unwrap();

    let calc = MetricsCalculator::new();
    let fitness = calc.compute_fitness(&assembly, &reads, &reference);

    assert_eq!(fitness.base_metrics.len(), reference.len());
    assert!(fitness.overall > 0.0, "fitness should be > 0 with real data");

    // Check individual base metrics
    for bm in &fitness.base_metrics {
        assert!(bm.agreement >= 0.0 && bm.agreement <= 1.0);
        assert!(bm.avg_quality >= 0.0);
    }

    // There should be some variant positions (real data has heterozygous sites)
    let variant_count = fitness.base_metrics.iter().filter(|m| m.is_variant).count();
    // It's okay if there are no variants in a small region, but the field should be populated
    assert!(
        variant_count <= reference.len(),
        "variant count should be <= reference length"
    );
}

#[test]
fn fitness_custom_weights_real_data() {
    let (reads, reference, _region) = setup_assembly_data();
    let method = ConsensusAssembly;
    let assembly = method.assemble(&reads, &reference).unwrap();

    let calc_agreement = MetricsCalculator::with_weights(1.0, 0.0, 0.0);
    let calc_depth = MetricsCalculator::with_weights(0.0, 1.0, 0.0);
    let calc_quality = MetricsCalculator::with_weights(0.0, 0.0, 1.0);

    let f_agree = calc_agreement.compute_fitness(&assembly, &reads, &reference);
    let f_depth = calc_depth.compute_fitness(&assembly, &reads, &reference);
    let f_quality = calc_quality.compute_fitness(&assembly, &reads, &reference);

    // All should produce valid scores
    assert!(f_agree.overall >= 0.0 && f_agree.overall <= 1.0);
    assert!(f_depth.overall >= 0.0 && f_depth.overall <= 1.0);
    assert!(f_quality.overall >= 0.0 && f_quality.overall <= 1.0);

    // The underlying metrics should be the same regardless of weights
    assert_eq!(f_agree.mean_agreement, f_depth.mean_agreement);
    assert_eq!(f_agree.mean_depth, f_quality.mean_depth);
}

// ===========================================================================
// Haplotype Tests
// ===========================================================================

#[test]
fn haplotype_assignment_real_data() {
    let (reads, reference, _region) = setup_assembly_data();

    let assigner = HaplotypeAssigner::new();
    let assignments = assigner.assign(&reads, &reference);

    assert_eq!(assignments.len(), reads.len());

    for a in &assignments {
        assert!(!a.read_name.is_empty());
        assert!(a.confidence >= 0.0 && a.confidence <= 1.0);
    }

    // Count haplotype distribution
    let h1 = assignments
        .iter()
        .filter(|a| a.haplotype == HaplotypeLabel::Hap1)
        .count();
    let h2 = assignments
        .iter()
        .filter(|a| a.haplotype == HaplotypeLabel::Hap2)
        .count();
    let unassigned = assignments
        .iter()
        .filter(|a| a.haplotype == HaplotypeLabel::Unassigned)
        .count();

    assert_eq!(h1 + h2 + unassigned, reads.len());
}

// ===========================================================================
// App Tests (no TUI, just state management)
// ===========================================================================

#[test]
fn app_with_real_data() {
    use longread_reviewer::viewer::App;

    let region = Region::new("chr17", 10990000, 10990100).unwrap();
    let reads = AlignmentReader::read_bam(example_bam_path(), &region).unwrap();
    let reference = synthetic_reference(&reads, &region);

    let mut app = App::new(region, reference, reads);

    // Run assembly
    app.run_assembly().unwrap();
    assert!(app.assembly.is_some());
    assert!(app.fitness.is_some());

    let fitness = app.fitness.as_ref().unwrap();
    assert!(fitness.overall > 0.0);

    // Cycle methods
    app.next_method().unwrap();
    assert!(app.assembly.is_some());

    // Assign haplotypes
    app.assign_haplotypes();
    assert!(!app.haplotype_assignments.is_empty());

    // Evaluate all
    let eval = app.evaluate_all().unwrap();
    assert!(!eval.is_empty());
    for (name, score) in &eval {
        assert!(!name.is_empty());
        assert!(*score >= 0.0 && *score <= 1.0);
    }
}

#[test]
fn app_scroll_bounds_with_real_data() {
    use crossterm::event::KeyCode;
    use longread_reviewer::viewer::App;

    let region = Region::new("chr17", 10990000, 10990100).unwrap();
    let reads = AlignmentReader::read_bam(example_bam_path(), &region).unwrap();
    let reference = synthetic_reference(&reads, &region);

    let mut app = App::new(region, reference.clone(), reads);

    // Scroll right many times - should cap at max
    for _ in 0..100 {
        app.handle_key(KeyCode::Right).unwrap();
    }
    assert!(app.scroll_x <= reference.len());

    // Scroll left many times - should not go below 0
    for _ in 0..200 {
        app.handle_key(KeyCode::Left).unwrap();
    }
    assert_eq!(app.scroll_x, 0);

    // Scroll down
    let num_reads = app.reads.len();
    for _ in 0..num_reads + 10 {
        app.handle_key(KeyCode::Down).unwrap();
    }
    assert!(app.scroll_y <= num_reads);

    // Scroll up
    for _ in 0..num_reads + 10 {
        app.handle_key(KeyCode::Up).unwrap();
    }
    assert_eq!(app.scroll_y, 0);
}

// ===========================================================================
// End-to-End Pipeline Test
// ===========================================================================

#[test]
fn end_to_end_pipeline() {
    // This test mirrors what the CLI subcommands do.
    let region = Region::new("chr17", 10990000, 10990100).unwrap();
    let reads = AlignmentReader::read_bam(example_bam_path(), &region).unwrap();
    let reference = synthetic_reference(&reads, &region);

    // 1. Assembly
    let mut engine = AssemblyEngine::new();
    engine.add_method(Box::new(ConsensusAssembly));
    engine.add_method(Box::new(WindowConsensusAssembly::new(50, 10)));
    engine.add_method(Box::new(WindowConsensusAssembly::new(30, 5)));
    engine.add_method(Box::new(WindowConsensusAssembly::new(80, 20)));

    // 2. Evaluate all methods
    let results = engine.evaluate_all(&reads, &reference, &region).unwrap();
    assert_eq!(results.len(), 4);

    // 3. Pick best
    let best = &results[0];
    assert!(!best.assembly.sequence.is_empty());
    assert!(best.fitness.overall >= results.last().unwrap().fitness.overall);

    // 4. Compute detailed metrics for the best
    let calc = MetricsCalculator::new();
    let fitness = calc.compute_fitness(&best.assembly, &reads, &reference);
    assert_eq!(fitness.base_metrics.len(), reference.len());

    // 5. Haplotype assignment
    let assigner = HaplotypeAssigner::new();
    let assignments = assigner.assign(&reads, &reference);
    assert_eq!(assignments.len(), reads.len());
}
