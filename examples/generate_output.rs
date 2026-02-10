//! Generate example output from the test BAM data.
//!
//! Produces `examples/example_output.txt` and `docs/index.html` automatically.
//!
//! Usage: cargo run --example generate_output

use std::fmt::Write as FmtWrite;
use std::fs;
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

    // Use a ~2 kb region around the complex structural event in the test data.
    let region = Region::new("chr17", 10990000, 10992000).unwrap();

    // Load data
    let reference = ReferenceGenome::from_file(&ref_path).expect("failed to read reference");
    let ref_seq = reference.fetch(&region).expect("failed to fetch region");
    let reads = AlignmentReader::read_bam(&bam_path, &region).expect("failed to read BAM");

    let mut out = String::new();

    writeln!(out, "=== Longread Reviewer — Example Output ===").unwrap();
    writeln!(out).unwrap();
    writeln!(out, "BAM:       tests/data/NA19240.chr17_fragment.bam").unwrap();
    writeln!(out, "Reference: tests/data/chr17_fragment.fa").unwrap();
    writeln!(out, "Region:    {}", region).unwrap();
    writeln!(out, "Reads:     {}", reads.len()).unwrap();
    writeln!(out).unwrap();

    // --- Read summary ---
    writeln!(out, "--- Read Summary ---").unwrap();
    writeln!(
        out,
        "{:<45} {:>6} {:>6} {:>8} {:>4} {:>3} {:>4}",
        "Name", "Start", "End", "Length", "MAPQ", "Rev", "HP"
    )
    .unwrap();
    writeln!(out, "{}", "-".repeat(85)).unwrap();
    for read in &reads {
        let hp_str = match read.haplotype_tag {
            Some(h) => format!("{}", h),
            None => "-".to_string(),
        };
        writeln!(
            out,
            "{:<45} {:>6} {:>6} {:>8} {:>4} {:>3} {:>4}",
            read.name,
            read.start,
            read.end,
            read.sequence.len(),
            read.mapq,
            if read.is_reverse { "Y" } else { "N" },
            hp_str,
        )
        .unwrap();
    }
    writeln!(out).unwrap();

    // --- Consensus assembly ---
    writeln!(out, "--- Consensus Assembly ---").unwrap();
    let method = ConsensusAssembly;
    let assembly = method
        .assemble(&reads, &ref_seq)
        .expect("assembly failed");

    let calc = MetricsCalculator::new();
    let fitness = calc.compute_fitness(&assembly, &reads, &ref_seq);

    writeln!(out, "Method:           {}", assembly.method_name).unwrap();
    writeln!(out, "Sequence length:  {}", assembly.sequence.len()).unwrap();
    writeln!(out, "Mean depth:       {:.1}", fitness.mean_depth).unwrap();
    writeln!(out, "Mean agreement:   {:.3}", fitness.mean_agreement).unwrap();
    writeln!(out, "Ref identity:     {:.3}", fitness.reference_identity).unwrap();
    writeln!(out, "Fitness score:    {:.3}", fitness.overall).unwrap();
    writeln!(
        out,
        "Consensus (first 80bp): {}",
        String::from_utf8_lossy(&assembly.sequence[..80.min(assembly.sequence.len())])
    )
    .unwrap();
    writeln!(out).unwrap();

    // --- Per-base metrics (first 20 positions) ---
    writeln!(out, "--- Per-base Metrics (first 20 positions) ---").unwrap();
    writeln!(
        out,
        "{:>6} {:>4} {:>4} {:>6} {:>10} {:>10} {:>8}",
        "Pos", "Ref", "Asm", "Depth", "Agreement", "AvgQual", "Variant"
    )
    .unwrap();
    writeln!(out, "{}", "-".repeat(56)).unwrap();
    for metric in fitness.base_metrics.iter().take(20) {
        writeln!(
            out,
            "{:>6} {:>4} {:>4} {:>6} {:>10.3} {:>10.1} {:>8}",
            metric.position,
            metric.ref_base as char,
            metric.assembly_base as char,
            metric.depth,
            metric.agreement,
            metric.avg_quality,
            if metric.is_variant { "YES" } else { "" },
        )
        .unwrap();
    }
    writeln!(out).unwrap();

    // --- Method evaluation ---
    writeln!(out, "--- Assembly Method Evaluation ---").unwrap();
    let mut engine = AssemblyEngine::new();
    engine.add_method(Box::new(ConsensusAssembly));
    engine.add_method(Box::new(WindowConsensusAssembly::new(100, 20)));
    engine.add_method(Box::new(WindowConsensusAssembly::new(50, 10)));
    engine.add_method(Box::new(WindowConsensusAssembly::new(200, 40)));

    let results = engine
        .evaluate_all(&reads, &ref_seq, &region)
        .expect("evaluation failed");

    writeln!(
        out,
        "{:<25} {:>10} {:>10} {:>10} {:>10}",
        "Method", "Agreement", "Depth", "RefIdent", "Fitness"
    )
    .unwrap();
    writeln!(out, "{}", "-".repeat(65)).unwrap();
    for result in &results {
        writeln!(
            out,
            "{:<25} {:>10.3} {:>10.1} {:>10.3} {:>10.3}",
            result.assembly.method_name,
            result.fitness.mean_agreement,
            result.fitness.mean_depth,
            result.fitness.reference_identity,
            result.fitness.overall,
        )
        .unwrap();
    }

    if let Some(best) = results.first() {
        writeln!(
            out,
            "\nBest method: {} (fitness: {:.3})",
            best.assembly.method_name, best.fitness.overall
        )
        .unwrap();
    }
    writeln!(out).unwrap();

    // --- Haplotype assignment ---
    writeln!(out, "--- Haplotype Assignment ---").unwrap();
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

    writeln!(out, "H1:         {} reads", hap1.len()).unwrap();
    writeln!(out, "H2:         {} reads", hap2.len()).unwrap();
    writeln!(out, "Unassigned: {} reads", unassigned.len()).unwrap();
    writeln!(out).unwrap();

    writeln!(
        out,
        "{:<45} {:>10} {:>12}",
        "Read", "Haplotype", "Confidence"
    )
    .unwrap();
    writeln!(out, "{}", "-".repeat(67)).unwrap();
    for assignment in &assignments {
        writeln!(
            out,
            "{:<45} {:>10} {:>12.3}",
            assignment.read_name, assignment.haplotype, assignment.confidence,
        )
        .unwrap();
    }

    // Print to stdout
    print!("{out}");

    // Write text output
    let txt_path = PathBuf::from(manifest).join("examples/example_output.txt");
    fs::write(&txt_path, &out).expect("failed to write example_output.txt");
    eprintln!("Wrote {}", txt_path.display());

    // Write docs/index.html
    let docs_dir = PathBuf::from(manifest).join("docs");
    fs::create_dir_all(&docs_dir).expect("failed to create docs directory");
    let html = generate_html(&out);
    let html_path = docs_dir.join("index.html");
    fs::write(&html_path, html).expect("failed to write docs/index.html");
    eprintln!("Wrote {}", html_path.display());
}

fn generate_html(text_output: &str) -> String {
    // Escape HTML special characters in the text output
    let escaped = text_output
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;");

    format!(
        r#"<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>longread-reviewer — Example Output</title>
<style>
  :root {{ --bg: #0d1117; --fg: #c9d1d9; --accent: #58a6ff; --border: #30363d; --code-bg: #161b22; }}
  * {{ margin: 0; padding: 0; box-sizing: border-box; }}
  body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif;
         background: var(--bg); color: var(--fg); line-height: 1.6; padding: 2rem; max-width: 960px; margin: 0 auto; }}
  h1 {{ color: var(--accent); margin-bottom: 0.5rem; }}
  p.subtitle {{ color: #8b949e; margin-bottom: 2rem; }}
  a {{ color: var(--accent); text-decoration: none; }}
  a:hover {{ text-decoration: underline; }}
  pre {{ background: var(--code-bg); border: 1px solid var(--border); border-radius: 6px;
         padding: 1rem; overflow-x: auto; font-size: 0.85rem; line-height: 1.5;
         font-family: 'SFMono-Regular', Consolas, 'Liberation Mono', Menlo, monospace; }}
  .note {{ background: #1c2128; border-left: 3px solid var(--accent); padding: 0.75rem 1rem;
           margin-bottom: 1.5rem; border-radius: 0 6px 6px 0; font-size: 0.9rem; }}
  footer {{ margin-top: 2rem; color: #484f58; font-size: 0.8rem; border-top: 1px solid var(--border); padding-top: 1rem; }}
</style>
</head>
<body>
<h1>longread-reviewer</h1>
<p class="subtitle">Terminal-based long-read alignment viewer with local assembly and fitness evaluation</p>
<div class="note">
  This page is <strong>auto-generated</strong> on every push by the
  <a href="https://github.com/rotblauer/longread-reviewer/actions">generate-examples</a>
  workflow using the bundled test BAM and reference FASTA.
</div>
<h2 style="margin-bottom:0.5rem;">Example Output</h2>
<pre>{escaped}</pre>
<footer>
  Generated by <code>cargo run --example generate_output</code> &mdash;
  <a href="https://github.com/rotblauer/longread-reviewer">source on GitHub</a>
</footer>
</body>
</html>
"#
    )
}
