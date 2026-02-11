use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::{Parser, Subcommand};

use longread_reviewer::alignment::AlignmentReader;
use longread_reviewer::assembly::engine::AssemblyEngine;
use longread_reviewer::assembly::method::{ConsensusAssembly, WindowConsensusAssembly};
use longread_reviewer::haplotype::HaplotypeAssigner;
use longread_reviewer::metrics::MetricsCalculator;
use longread_reviewer::reference::ReferenceGenome;
use longread_reviewer::region::Region;
use longread_reviewer::viewer::App;

#[derive(Parser)]
#[command(
    name = "longread-reviewer",
    about = "Terminal-based long-read alignment viewer with local assembly and fitness evaluation",
    version
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Open an interactive TUI viewer for a region
    View {
        /// Path to BAM file (must be sorted and indexed)
        #[arg(short, long)]
        bam: PathBuf,

        /// Path to reference FASTA file
        #[arg(short, long)]
        reference: PathBuf,

        /// Region to view (format: chr:start-end)
        #[arg(short = 'L', long)]
        region: String,
    },

    /// Generate a local assembly for a region and print results
    Assemble {
        /// Path to BAM file (must be sorted and indexed)
        #[arg(short, long)]
        bam: PathBuf,

        /// Path to reference FASTA file
        #[arg(short, long)]
        reference: PathBuf,

        /// Region to assemble (format: chr:start-end)
        #[arg(short = 'L', long)]
        region: String,

        /// Assembly method to use (consensus, window_consensus)
        #[arg(short, long, default_value = "consensus")]
        method: String,

        /// Window size for window-based methods
        #[arg(short, long, default_value = "100")]
        window_size: usize,

        /// Window overlap for window-based methods
        #[arg(short, long, default_value = "20")]
        overlap: usize,
    },

    /// Evaluate multiple assembly methods and rank by fitness
    Evaluate {
        /// Path to BAM file (must be sorted and indexed)
        #[arg(short, long)]
        bam: PathBuf,

        /// Path to reference FASTA file
        #[arg(short, long)]
        reference: PathBuf,

        /// Region to evaluate (format: chr:start-end)
        #[arg(short = 'L', long)]
        region: String,
    },
}

fn load_data(
    bam_path: &Path,
    ref_path: &Path,
    region_str: &str,
) -> Result<(Region, Vec<u8>, Vec<longread_reviewer::alignment::AlignedRead>)> {
    let region: Region = region_str
        .parse()
        .context("failed to parse region")?;

    let reference = ReferenceGenome::from_file(ref_path)?;
    let ref_seq = reference.fetch(&region)?;
    let reads = AlignmentReader::read_bam(bam_path, &region)?;

    Ok((region, ref_seq, reads))
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::View {
            bam,
            reference,
            region,
        } => {
            let (region, ref_seq, reads) = load_data(&bam, &reference, &region)?;
            println!(
                "Loaded {} reads in region {}",
                reads.len(),
                region
            );

            let mut app = App::new(region, ref_seq, reads);
            app.run_tui()?;
        }

        Commands::Assemble {
            bam,
            reference,
            region,
            method,
            window_size,
            overlap,
        } => {
            let (region, ref_seq, reads) = load_data(&bam, &reference, &region)?;
            println!(
                "Loaded {} reads in region {}",
                reads.len(),
                region
            );

            let mut engine = AssemblyEngine::new();
            match method.as_str() {
                "consensus" => engine.add_method(Box::new(ConsensusAssembly)),
                "window_consensus" => {
                    engine.add_method(Box::new(WindowConsensusAssembly::new(window_size, overlap)))
                }
                _ => anyhow::bail!("unknown assembly method: {method}"),
            };

            let result = engine
                .run_method(&method, &reads, &ref_seq, region.start)?
                .context("assembly method not found")?;

            let calc = MetricsCalculator::new();
            let fitness = calc.compute_fitness(&result, &reads, &ref_seq, region.start);

            println!("Assembly method: {}", result.method_name);
            println!(
                "Sequence:        {}",
                String::from_utf8_lossy(&result.sequence)
            );
            println!("Length:          {}", result.sequence.len());
            println!("Mean depth:      {:.1}", fitness.mean_depth);
            println!("Mean agreement:  {:.3}", fitness.mean_agreement);
            println!("Ref identity:    {:.3}", fitness.reference_identity);
            println!("Fitness score:   {:.3}", fitness.overall);

            // Haplotype analysis
            let assigner = HaplotypeAssigner::new();
            let assignments = assigner.assign(&reads, &ref_seq);
            let hap1_count = assignments
                .iter()
                .filter(|a| a.haplotype == longread_reviewer::haplotype::HaplotypeLabel::Hap1)
                .count();
            let hap2_count = assignments
                .iter()
                .filter(|a| a.haplotype == longread_reviewer::haplotype::HaplotypeLabel::Hap2)
                .count();
            let unassigned = assignments
                .iter()
                .filter(|a| {
                    a.haplotype == longread_reviewer::haplotype::HaplotypeLabel::Unassigned
                })
                .count();

            println!("\nHaplotype assignment:");
            println!("  H1: {hap1_count} reads");
            println!("  H2: {hap2_count} reads");
            println!("  Unassigned: {unassigned} reads");
        }

        Commands::Evaluate {
            bam,
            reference,
            region,
        } => {
            let (region, ref_seq, reads) = load_data(&bam, &reference, &region)?;
            println!(
                "Loaded {} reads in region {}",
                reads.len(),
                region
            );

            let mut engine = AssemblyEngine::new();
            engine.add_method(Box::new(ConsensusAssembly));
            engine.add_method(Box::new(WindowConsensusAssembly::new(100, 20)));
            engine.add_method(Box::new(WindowConsensusAssembly::new(50, 10)));
            engine.add_method(Box::new(WindowConsensusAssembly::new(200, 40)));

            let results = engine.evaluate_all(&reads, &ref_seq, &region)?;

            println!("\nAssembly method evaluation (ranked by fitness):");
            println!("{:<25} {:>10} {:>10} {:>10} {:>10}", "Method", "Agreement", "Depth", "RefIdent", "Fitness");
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
        }
    }

    Ok(())
}
