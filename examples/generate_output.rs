//! Generate example output from the test BAM data.
//!
//! Usage: cargo run --example generate_output

use std::path::PathBuf;

use longread_reviewer::alignment::AlignmentReader;
use longread_reviewer::assembly::engine::AssemblyEngine;
use longread_reviewer::assembly::method::{AssemblyMethod, ConsensusAssembly, WindowConsensusAssembly};
use longread_reviewer::haplotype::{HaplotypeAssigner, HaplotypeLabel};
use longread_reviewer::metrics::MetricsCalculator;
use longread_reviewer::reference::ReferenceGenome;
use longread_reviewer::region::Region;

fn main() {
    let manifest = env!("CARGO_MANIFEST_DIR");
    let bam_path = PathBuf::from(manifest).join("tests/data/NA19240.chr17_fragment.bam");
    let ref_path = PathBuf::from(manifest).join("tests/data/chr17_fragment.fa");

    let region = Region::new("chr17", 1, 2001).unwrap();

    // Load data
    let reference = ReferenceGenome::from_file(&ref_path).expect("failed to read reference");
    let ref_seq = reference.fetch(&region).expect("failed to fetch region");
    let reads = AlignmentReader::read_bam(&bam_path, &region).expect("failed to read BAM");

    println!("=== Longread Reviewer — Example Output ===");
    println!();
    println!("BAM:       {}", bam_path.display());
    println!("Reference: {}", ref_path.display());
    println!("Region:    {}", region);
    println!("Reads:     {}", reads.len());
    println!();

    // --- Read summary ---
    println!("--- Read Summary ---");
    println!(
        "{:<45} {:>6} {:>6} {:>8} {:>4} {:>3} {:>4}",
        "Name", "Start", "End", "Length", "MAPQ", "Rev", "HP"
    );
    println!("{}", "-".repeat(85));
    for read in &reads {
        let hp_str = match read.haplotype_tag {
            Some(h) => format!("{}", h),
            None => "-".to_string(),
        };
        println!(
            "{:<45} {:>6} {:>6} {:>8} {:>4} {:>3} {:>4}",
            read.name,
            read.start,
            read.end,
            read.sequence.len(),
            read.mapq,
            if read.is_reverse { "Y" } else { "N" },
            hp_str,
        );
    }
    println!();

    // --- Consensus assembly ---
    println!("--- Consensus Assembly ---");
    let method = ConsensusAssembly;
    let assembly = method
        .assemble(&reads, &ref_seq)
        .expect("assembly failed");

    let calc = MetricsCalculator::new();
    let fitness = calc.compute_fitness(&assembly, &reads, &ref_seq);

    println!("Method:           {}", assembly.method_name);
    println!("Sequence length:  {}", assembly.sequence.len());
    println!("Mean depth:       {:.1}", fitness.mean_depth);
    println!("Mean agreement:   {:.3}", fitness.mean_agreement);
    println!("Ref identity:     {:.3}", fitness.reference_identity);
    println!("Fitness score:    {:.3}", fitness.overall);
    println!(
        "Consensus (first 80bp): {}",
        String::from_utf8_lossy(&assembly.sequence[..80.min(assembly.sequence.len())])
    );
    println!();

    // --- Per-base metrics (first 20 positions) ---
    println!("--- Per-base Metrics (first 20 positions) ---");
    println!(
        "{:>6} {:>4} {:>4} {:>6} {:>10} {:>10} {:>8}",
        "Pos", "Ref", "Asm", "Depth", "Agreement", "AvgQual", "Variant"
    );
    println!("{}", "-".repeat(56));
    for metric in fitness.base_metrics.iter().take(20) {
        println!(
            "{:>6} {:>4} {:>4} {:>6} {:>10.3} {:>10.1} {:>8}",
            metric.position,
            metric.ref_base as char,
            metric.assembly_base as char,
            metric.depth,
            metric.agreement,
            metric.avg_quality,
            if metric.is_variant { "YES" } else { "" },
        );
    }
    println!();

    // --- Method evaluation ---
    println!("--- Assembly Method Evaluation ---");
    let mut engine = AssemblyEngine::new();
    engine.add_method(Box::new(ConsensusAssembly));
    engine.add_method(Box::new(WindowConsensusAssembly::new(100, 20)));
    engine.add_method(Box::new(WindowConsensusAssembly::new(50, 10)));
    engine.add_method(Box::new(WindowConsensusAssembly::new(200, 40)));

    let results = engine
        .evaluate_all(&reads, &ref_seq, &region)
        .expect("evaluation failed");

    println!(
        "{:<25} {:>10} {:>10} {:>10} {:>10}",
        "Method", "Agreement", "Depth", "RefIdent", "Fitness"
    );
    println!("{}", "-".repeat(65));
    for result in &results {
        println!(
            "{:<25} {:>10.3} {:>10.1} {:>10.3} {:>10.3}",
            result.assembly.method_name,
            result.fitness.mean_agreement,
            result.fitness.mean_depth,
            result.fitness.reference_identity,
            result.fitness.overall,
        );
    }

    if let Some(best) = results.first() {
        println!(
            "\nBest method: {} (fitness: {:.3})",
            best.assembly.method_name, best.fitness.overall
        );
    }
    println!();

    // --- Haplotype assignment ---
    println!("--- Haplotype Assignment ---");
    let assigner = HaplotypeAssigner::new();
    let assignments = assigner.assign(&reads, &ref_seq);

    let hap1: Vec<_> = assignments
        .iter()
        .filter(|a| a.haplotype == HaplotypeLabel::Hap1)
        .collect();
    let hap2: Vec<_> = assignments
        .iter()
        .filter(|a| a.haplotype == HaplotypeLabel::Hap2)
        .collect();
    let unassigned: Vec<_> = assignments
        .iter()
        .filter(|a| a.haplotype == HaplotypeLabel::Unassigned)
        .collect();

    println!("H1:         {} reads", hap1.len());
    println!("H2:         {} reads", hap2.len());
    println!("Unassigned: {} reads", unassigned.len());
    println!();

    println!(
        "{:<45} {:>10} {:>12}",
        "Read", "Haplotype", "Confidence"
    );
    println!("{}", "-".repeat(67));
    for assignment in &assignments {
        println!(
            "{:<45} {:>10} {:>12.3}",
            assignment.read_name, assignment.haplotype, assignment.confidence,
        );
    }
}
