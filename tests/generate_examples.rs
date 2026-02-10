//! Generate example output files from the example BAM file.
//!
//! Run with: cargo test --test generate_examples -- --nocapture --ignored
//!
//! This creates files under `example/output/` showing representative reviewer output.

use std::io::Write;
use std::path::Path;

use longread_reviewer::alignment::{AlignedRead, AlignmentReader, CigarOp};
use longread_reviewer::assembly::engine::AssemblyEngine;
use longread_reviewer::assembly::method::{
    AssemblyMethod, ConsensusAssembly, WindowConsensusAssembly,
};
use longread_reviewer::haplotype::{HaplotypeAssigner, HaplotypeLabel};
use longread_reviewer::metrics::MetricsCalculator;
use longread_reviewer::region::Region;

const EXAMPLE_BAM: &str =
    "example/NA19240_2020_merged.ccs.hg38.aligned.chr17_10958130_11017414.bam";

fn example_bam_path() -> &'static Path {
    Path::new(EXAMPLE_BAM)
}

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

#[test]
#[ignore]
fn generate_read_summary() {
    let region = Region::new("chr17", 10958130, 11017414).unwrap();
    let reads = AlignmentReader::read_bam(example_bam_path(), &region).unwrap();

    let out_path = Path::new("example/output/read_summary.txt");
    let mut f = std::fs::File::create(out_path).unwrap();

    writeln!(f, "=== Read Summary ===").unwrap();
    writeln!(
        f,
        "BAM file: NA19240_2020_merged.ccs.hg38.aligned.chr17_10958130_11017414.bam"
    )
    .unwrap();
    writeln!(f, "Region: {}", region).unwrap();
    writeln!(f, "Total reads: {}", reads.len()).unwrap();
    writeln!(f).unwrap();

    let min_start = reads.iter().map(|r| r.start).min().unwrap();
    let max_end = reads.iter().map(|r| r.end).max().unwrap();
    let avg_len =
        reads.iter().map(|r| r.sequence.len() as f64).sum::<f64>() / reads.len() as f64;
    let avg_mapq = reads.iter().map(|r| r.mapq as f64).sum::<f64>() / reads.len() as f64;
    let fwd = reads.iter().filter(|r| !r.is_reverse).count();
    let rev = reads.iter().filter(|r| r.is_reverse).count();

    let min_len = reads.iter().map(|r| r.sequence.len()).min().unwrap();
    let max_len = reads.iter().map(|r| r.sequence.len()).max().unwrap();

    writeln!(f, "--- Coordinate Range ---").unwrap();
    writeln!(f, "Leftmost alignment start: {}", min_start).unwrap();
    writeln!(f, "Rightmost alignment end:  {}", max_end).unwrap();
    writeln!(f).unwrap();

    writeln!(f, "--- Read Lengths ---").unwrap();
    writeln!(f, "Min read length: {}", min_len).unwrap();
    writeln!(f, "Max read length: {}", max_len).unwrap();
    writeln!(f, "Mean read length: {:.0}", avg_len).unwrap();
    writeln!(f).unwrap();

    writeln!(f, "--- Mapping Quality ---").unwrap();
    writeln!(f, "Mean MAPQ: {:.1}", avg_mapq).unwrap();
    writeln!(f).unwrap();

    writeln!(f, "--- Strand Distribution ---").unwrap();
    writeln!(f, "Forward: {} ({:.1}%)", fwd, fwd as f64 / reads.len() as f64 * 100.0).unwrap();
    writeln!(f, "Reverse: {} ({:.1}%)", rev, rev as f64 / reads.len() as f64 * 100.0).unwrap();
    writeln!(f).unwrap();

    // CIGAR operation statistics
    let mut match_count = 0u64;
    let mut ins_count = 0u64;
    let mut del_count = 0u64;
    let mut sc_count = 0u64;
    let mut match_bases = 0u64;
    let mut ins_bases = 0u64;
    let mut del_bases = 0u64;
    let mut sc_bases = 0u64;

    for read in &reads {
        for op in &read.cigar {
            match op {
                CigarOp::Match(n) => {
                    match_count += 1;
                    match_bases += *n as u64;
                }
                CigarOp::Insertion(n) => {
                    ins_count += 1;
                    ins_bases += *n as u64;
                }
                CigarOp::Deletion(n) => {
                    del_count += 1;
                    del_bases += *n as u64;
                }
                CigarOp::SoftClip(n) => {
                    sc_count += 1;
                    sc_bases += *n as u64;
                }
                CigarOp::HardClip(_) => {}
            }
        }
    }

    writeln!(f, "--- CIGAR Operations ---").unwrap();
    writeln!(f, "Match operations:     {} ({} bases)", match_count, match_bases).unwrap();
    writeln!(f, "Insertion operations: {} ({} bases)", ins_count, ins_bases).unwrap();
    writeln!(f, "Deletion operations:  {} ({} bases)", del_count, del_bases).unwrap();
    writeln!(f, "Soft-clip operations: {} ({} bases)", sc_count, sc_bases).unwrap();
    writeln!(f).unwrap();

    // Per-read table (first 20 reads)
    writeln!(f, "--- First 20 Reads ---").unwrap();
    writeln!(
        f,
        "{:<45} {:>10} {:>10} {:>6} {:>8} {:>7}",
        "Name", "Start", "End", "MAPQ", "SeqLen", "Strand"
    )
    .unwrap();
    writeln!(f, "{}", "-".repeat(90)).unwrap();
    for read in reads.iter().take(20) {
        writeln!(
            f,
            "{:<45} {:>10} {:>10} {:>6} {:>8} {:>7}",
            read.name,
            read.start,
            read.end,
            read.mapq,
            read.sequence.len(),
            if read.is_reverse { "REV" } else { "FWD" }
        )
        .unwrap();
    }

    println!("Generated: {}", out_path.display());
}

#[test]
#[ignore]
fn generate_assembly_output() {
    let region = Region::new("chr17", 10990000, 10990100).unwrap();
    let reads = AlignmentReader::read_bam(example_bam_path(), &region).unwrap();
    let reference = synthetic_reference(&reads, &region);

    let out_path = Path::new("example/output/assembly_results.txt");
    let mut f = std::fs::File::create(out_path).unwrap();

    writeln!(f, "=== Assembly Results ===").unwrap();
    writeln!(f, "Region: {}", region).unwrap();
    writeln!(f, "Reads used: {}", reads.len()).unwrap();
    writeln!(f, "Reference length: {}", reference.len()).unwrap();
    writeln!(f).unwrap();

    // Run consensus assembly
    let consensus = ConsensusAssembly;
    let result = consensus.assemble(&reads, &reference).unwrap();

    writeln!(f, "--- Consensus Assembly ---").unwrap();
    writeln!(f, "Method: {}", result.method_name).unwrap();
    writeln!(f, "Assembly length: {}", result.sequence.len()).unwrap();
    writeln!(
        f,
        "Sequence: {}",
        String::from_utf8_lossy(&result.sequence)
    )
    .unwrap();
    writeln!(f).unwrap();

    writeln!(f, "Per-base depth:").unwrap();
    let depth_str: String = result
        .depth
        .iter()
        .map(|d| format!("{:>3}", d))
        .collect::<Vec<_>>()
        .join(" ");
    writeln!(f, "  {}", depth_str).unwrap();
    writeln!(f).unwrap();

    writeln!(f, "Per-base confidence:").unwrap();
    let conf_str: String = result
        .confidence
        .iter()
        .map(|c| format!("{:.2}", c))
        .collect::<Vec<_>>()
        .join(" ");
    writeln!(f, "  {}", conf_str).unwrap();
    writeln!(f).unwrap();

    // Window consensus
    for (ws, ov) in &[(50, 10), (30, 5), (80, 20)] {
        let method = WindowConsensusAssembly::new(*ws, *ov);
        let result = method.assemble(&reads, &reference).unwrap();

        writeln!(
            f,
            "--- Window Consensus (window={}, overlap={}) ---",
            ws, ov
        )
        .unwrap();
        writeln!(f, "Assembly length: {}", result.sequence.len()).unwrap();
        writeln!(
            f,
            "Sequence: {}",
            String::from_utf8_lossy(&result.sequence)
        )
        .unwrap();

        let avg_depth =
            result.depth.iter().map(|&d| d as f64).sum::<f64>() / result.depth.len() as f64;
        let avg_conf = result.confidence.iter().sum::<f64>() / result.confidence.len() as f64;
        writeln!(f, "Average depth: {:.1}", avg_depth).unwrap();
        writeln!(f, "Average confidence: {:.3}", avg_conf).unwrap();
        writeln!(f).unwrap();
    }

    println!("Generated: {}", out_path.display());
}

#[test]
#[ignore]
fn generate_evaluation_output() {
    let region = Region::new("chr17", 10990000, 10990100).unwrap();
    let reads = AlignmentReader::read_bam(example_bam_path(), &region).unwrap();
    let reference = synthetic_reference(&reads, &region);

    let out_path = Path::new("example/output/evaluation_results.txt");
    let mut f = std::fs::File::create(out_path).unwrap();

    writeln!(f, "=== Assembly Method Evaluation ===").unwrap();
    writeln!(f, "Region: {}", region).unwrap();
    writeln!(f, "Reads used: {}", reads.len()).unwrap();
    writeln!(f).unwrap();

    let mut engine = AssemblyEngine::new();
    engine.add_method(Box::new(ConsensusAssembly));
    engine.add_method(Box::new(WindowConsensusAssembly::new(50, 10)));
    engine.add_method(Box::new(WindowConsensusAssembly::new(30, 5)));
    engine.add_method(Box::new(WindowConsensusAssembly::new(80, 20)));

    let results = engine.evaluate_all(&reads, &reference, &region).unwrap();

    writeln!(
        f,
        "{:<25} {:>10} {:>10} {:>10} {:>10}",
        "Method", "Agreement", "Depth", "RefIdent", "Fitness"
    )
    .unwrap();
    writeln!(f, "{}", "-".repeat(65)).unwrap();

    for result in &results {
        writeln!(
            f,
            "{:<25} {:>10.3} {:>10.1} {:>10.3} {:>10.3}",
            result.assembly.method_name,
            result.fitness.mean_agreement,
            result.fitness.mean_depth,
            result.fitness.reference_identity,
            result.fitness.overall,
        )
        .unwrap();
    }
    writeln!(f).unwrap();

    if let Some(best) = results.first() {
        writeln!(
            f,
            "Best method: {} (fitness: {:.3})",
            best.assembly.method_name, best.fitness.overall
        )
        .unwrap();
    }

    println!("Generated: {}", out_path.display());
}

#[test]
#[ignore]
fn generate_fitness_output() {
    let region = Region::new("chr17", 10990000, 10990100).unwrap();
    let reads = AlignmentReader::read_bam(example_bam_path(), &region).unwrap();
    let reference = synthetic_reference(&reads, &region);

    let out_path = Path::new("example/output/fitness_metrics.txt");
    let mut f = std::fs::File::create(out_path).unwrap();

    let method = ConsensusAssembly;
    let assembly = method.assemble(&reads, &reference).unwrap();
    let calc = MetricsCalculator::new();
    let fitness = calc.compute_fitness(&assembly, &reads, &reference);

    writeln!(f, "=== Per-Base Fitness Metrics ===").unwrap();
    writeln!(f, "Region: {}", region).unwrap();
    writeln!(f, "Method: consensus").unwrap();
    writeln!(f).unwrap();

    writeln!(f, "--- Summary ---").unwrap();
    writeln!(f, "Mean agreement:    {:.3}", fitness.mean_agreement).unwrap();
    writeln!(f, "Mean depth:        {:.1}", fitness.mean_depth).unwrap();
    writeln!(f, "Reference identity: {:.3}", fitness.reference_identity).unwrap();
    writeln!(f, "Overall fitness:   {:.3}", fitness.overall).unwrap();
    writeln!(f).unwrap();

    let variant_count = fitness.base_metrics.iter().filter(|m| m.is_variant).count();
    writeln!(
        f,
        "Variant positions: {} / {} ({:.1}%)",
        variant_count,
        fitness.base_metrics.len(),
        variant_count as f64 / fitness.base_metrics.len() as f64 * 100.0
    )
    .unwrap();
    writeln!(f).unwrap();

    writeln!(f, "--- Per-Base Detail ---").unwrap();
    writeln!(
        f,
        "{:>10} {:>5} {:>5} {:>6} {:>10} {:>10} {:>8}",
        "Position", "Ref", "Asm", "Depth", "Agreement", "AvgQual", "Variant"
    )
    .unwrap();
    writeln!(f, "{}", "-".repeat(60)).unwrap();

    for bm in &fitness.base_metrics {
        writeln!(
            f,
            "{:>10} {:>5} {:>5} {:>6} {:>10.3} {:>10.1} {:>8}",
            bm.position,
            bm.ref_base as char,
            bm.assembly_base as char,
            bm.depth,
            bm.agreement,
            bm.avg_quality,
            if bm.is_variant { "YES" } else { "." },
        )
        .unwrap();
    }

    println!("Generated: {}", out_path.display());
}

#[test]
#[ignore]
fn generate_haplotype_output() {
    let region = Region::new("chr17", 10990000, 10990100).unwrap();
    let reads = AlignmentReader::read_bam(example_bam_path(), &region).unwrap();
    let reference = synthetic_reference(&reads, &region);

    let out_path = Path::new("example/output/haplotype_assignments.txt");
    let mut f = std::fs::File::create(out_path).unwrap();

    let assigner = HaplotypeAssigner::new();
    let assignments = assigner.assign(&reads, &reference);

    writeln!(f, "=== Haplotype Assignments ===").unwrap();
    writeln!(f, "Region: {}", region).unwrap();
    writeln!(f, "Total reads: {}", reads.len()).unwrap();
    writeln!(f).unwrap();

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

    writeln!(f, "--- Summary ---").unwrap();
    writeln!(f, "Haplotype 1: {} reads", h1).unwrap();
    writeln!(f, "Haplotype 2: {} reads", h2).unwrap();
    writeln!(f, "Unassigned:  {} reads", unassigned).unwrap();
    writeln!(f).unwrap();

    writeln!(f, "--- Per-Read Assignments ---").unwrap();
    writeln!(
        f,
        "{:<45} {:>10} {:>12}",
        "Read Name", "Haplotype", "Confidence"
    )
    .unwrap();
    writeln!(f, "{}", "-".repeat(70)).unwrap();

    for a in &assignments {
        writeln!(
            f,
            "{:<45} {:>10} {:>12.3}",
            a.read_name,
            a.haplotype.to_string(),
            a.confidence,
        )
        .unwrap();
    }

    println!("Generated: {}", out_path.display());
}
