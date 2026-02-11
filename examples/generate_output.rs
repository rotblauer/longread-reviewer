//! Generate example output from the test BAM data.
//!
//! Scans the full test region chr17:10958130-11017414 for complex events,
//! evaluates multiple assembly methods, and produces visualizations useful
//! for assessing local realignment/assembly quality.
//!
//! Produces `examples/example_output.txt` and `docs/index.html` automatically.
//!
//! Usage: cargo run --example generate_output

use std::collections::HashMap;
use std::fmt::Write as FmtWrite;
use std::fs;
use std::path::PathBuf;

use longread_reviewer::alignment::AlignmentReader;
use longread_reviewer::alignment::CigarOp;
use longread_reviewer::assembly::engine::AssemblyEngine;
use longread_reviewer::assembly::method::{AssemblyMethod, AssemblyResult, ConsensusAssembly, WindowConsensusAssembly};
use longread_reviewer::haplotype::{HaplotypeAssigner, HaplotypeLabel};
use longread_reviewer::metrics::{FitnessScore, MetricsCalculator};
use longread_reviewer::reference::ReferenceGenome;
use longread_reviewer::region::Region;
use longread_reviewer::alignment::AlignedRead;

/// A detected complex event region with metrics.
#[derive(Debug, Clone)]
struct ComplexEvent {
    start: u64,
    end: u64,
    variant_count: usize,
    mean_agreement: f64,
    depth_variance: f64,
    event_type: String,
}

/// Scan for regions with high variant density or coverage anomalies.
fn find_complex_events(
    fitness: &FitnessScore,
    window_size: usize,
    min_variants: usize,
    min_disagreement: f64,
) -> Vec<ComplexEvent> {
    let mut events = Vec::new();
    let metrics = &fitness.base_metrics;

    if metrics.len() < window_size {
        return events;
    }

    let mut i = 0;
    while i + window_size <= metrics.len() {
        let window = &metrics[i..i + window_size];

        let variant_count = window.iter().filter(|m| m.is_variant).count();
        let mean_agreement: f64 = window.iter().map(|m| m.agreement).sum::<f64>() / window_size as f64;
        let depths: Vec<f64> = window.iter().map(|m| m.depth as f64).collect();
        let mean_depth = depths.iter().sum::<f64>() / depths.len() as f64;
        let depth_variance = depths.iter().map(|d| (d - mean_depth).powi(2)).sum::<f64>() / depths.len() as f64;

        let is_complex = variant_count >= min_variants || mean_agreement < (1.0 - min_disagreement);

        if is_complex {
            let start = window.first().unwrap().position;
            let end = window.last().unwrap().position;

            // Classify event type
            let event_type = if variant_count >= min_variants * 2 {
                "HIGH_VARIANT_DENSITY"
            } else if depth_variance > mean_depth * 2.0 {
                "COVERAGE_ANOMALY"
            } else if mean_agreement < 0.7 {
                "LOW_AGREEMENT"
            } else {
                "VARIANT_CLUSTER"
            };

            // Extend to include adjacent complex windows
            let mut extended_end = end;
            let mut j = i + window_size;
            while j + window_size <= metrics.len() {
                let next_window = &metrics[j..j + window_size];
                let next_variants = next_window.iter().filter(|m| m.is_variant).count();
                let next_agreement: f64 = next_window.iter().map(|m| m.agreement).sum::<f64>() / window_size as f64;

                if next_variants >= min_variants / 2 || next_agreement < (1.0 - min_disagreement / 2.0) {
                    extended_end = next_window.last().unwrap().position;
                    j += window_size / 2;
                } else {
                    break;
                }
            }

            events.push(ComplexEvent {
                start,
                end: extended_end,
                variant_count,
                mean_agreement,
                depth_variance,
                event_type: event_type.to_string(),
            });

            // Skip past this event
            i = j;
        } else {
            i += window_size / 2;
        }
    }

    events
}

/// Generate a coverage histogram for a region.
fn coverage_histogram(fitness: &FitnessScore, bins: usize) -> Vec<(u32, usize)> {
    let mut depth_counts: HashMap<u32, usize> = HashMap::new();
    for metric in &fitness.base_metrics {
        *depth_counts.entry(metric.depth).or_insert(0) += 1;
    }

    let max_depth = *depth_counts.keys().max().unwrap_or(&0);
    let bin_size = (max_depth as usize / bins).max(1);

    let mut histogram: Vec<(u32, usize)> = Vec::new();
    for bin in 0..bins {
        let bin_start = (bin * bin_size) as u32;
        let bin_end = ((bin + 1) * bin_size) as u32;
        let count: usize = depth_counts
            .iter()
            .filter(|entry| *entry.0 >= bin_start && *entry.0 < bin_end)
            .map(|entry| *entry.1)
            .sum();
        histogram.push((bin_start, count));
    }
    histogram
}

/// Render ASCII coverage plot for a region.
fn render_coverage_track(fitness: &FitnessScore, width: usize, height: usize) -> Vec<String> {
    let metrics = &fitness.base_metrics;
    if metrics.is_empty() {
        return vec!["No data".to_string()];
    }

    let step = metrics.len().max(1) / width.max(1);
    let step = step.max(1);

    let mut sampled_depths: Vec<f64> = Vec::new();
    for i in (0..metrics.len()).step_by(step) {
        let end = (i + step).min(metrics.len());
        let avg_depth: f64 = metrics[i..end].iter().map(|m| m.depth as f64).sum::<f64>() / (end - i) as f64;
        sampled_depths.push(avg_depth);
    }
    sampled_depths.truncate(width);

    let max_depth = sampled_depths.iter().cloned().fold(1.0_f64, f64::max);

    let mut lines: Vec<String> = Vec::new();
    for row in (0..height).rev() {
        let threshold = (row as f64 / height as f64) * max_depth;
        let mut line = String::new();
        for &depth in &sampled_depths {
            if depth >= threshold {
                line.push('█');
            } else if depth >= threshold - max_depth / height as f64 / 2.0 {
                line.push('▄');
            } else {
                line.push(' ');
            }
        }
        lines.push(line);
    }

    // Add scale
    lines.push("─".repeat(width.min(sampled_depths.len())));
    lines.push(format!("Depth: 0{:>width$}", format!("{:.0}", max_depth), width = width - 2));

    lines
}

/// Render ASCII variant density track.
fn render_variant_track(fitness: &FitnessScore, width: usize) -> String {
    let metrics = &fitness.base_metrics;
    if metrics.is_empty() {
        return String::new();
    }

    let step = metrics.len().max(1) / width.max(1);
    let step = step.max(1);

    let mut track = String::new();
    for i in (0..metrics.len()).step_by(step) {
        let end = (i + step).min(metrics.len());
        let variant_count = metrics[i..end].iter().filter(|m| m.is_variant).count();
        let density = variant_count as f64 / (end - i) as f64;

        let ch = if density > 0.5 {
            '█'
        } else if density > 0.2 {
            '▓'
        } else if density > 0.1 {
            '▒'
        } else if density > 0.0 {
            '░'
        } else {
            '·'
        };
        track.push(ch);

        if track.len() >= width {
            break;
        }
    }
    track
}

/// Render ASCII agreement heatmap.
fn render_agreement_track(fitness: &FitnessScore, width: usize) -> String {
    let metrics = &fitness.base_metrics;
    if metrics.is_empty() {
        return String::new();
    }

    let step = metrics.len().max(1) / width.max(1);
    let step = step.max(1);

    let mut track = String::new();
    for i in (0..metrics.len()).step_by(step) {
        let end = (i + step).min(metrics.len());
        let mean_agreement: f64 = metrics[i..end].iter().map(|m| m.agreement).sum::<f64>() / (end - i) as f64;

        // Color scale: green (high agreement) -> yellow -> red (low agreement)
        let ch = if mean_agreement >= 0.95 {
            '█'
        } else if mean_agreement >= 0.9 {
            '▓'
        } else if mean_agreement >= 0.8 {
            '▒'
        } else if mean_agreement >= 0.7 {
            '░'
        } else {
            '·'
        };
        track.push(ch);

        if track.len() >= width {
            break;
        }
    }
    track
}

/// Compare assemblies at a specific event region.
fn compare_assemblies_at_event(
    event: &ComplexEvent,
    reads: &[AlignedRead],
    reference: &ReferenceGenome,
    region: &Region,
) -> Vec<(String, FitnessScore, AssemblyResult)> {
    let event_region = Region::new(&region.chrom, event.start, event.end).unwrap();
    let ref_seq = match reference.fetch(&event_region) {
        Ok(seq) => seq,
        Err(_) => return Vec::new(),
    };

    // Filter reads overlapping this event
    let event_reads: Vec<_> = reads
        .iter()
        .filter(|r| r.start <= event.end && r.end >= event.start)
        .cloned()
        .collect();

    if event_reads.is_empty() {
        return Vec::new();
    }

    let calc = MetricsCalculator::new();
    let methods: Vec<Box<dyn AssemblyMethod>> = vec![
        Box::new(ConsensusAssembly),
        Box::new(WindowConsensusAssembly::new(50, 10)),
        Box::new(WindowConsensusAssembly::new(100, 20)),
        Box::new(WindowConsensusAssembly::new(200, 40)),
    ];

    let mut results = Vec::new();
    for method in methods {
        if let Ok(assembly) = method.assemble(&event_reads, &ref_seq, event.start) {
            let fitness = calc.compute_fitness(&assembly, &event_reads, &ref_seq, event.start);
            results.push((method.name().to_string(), fitness, assembly));
        }
    }

    results.sort_by(|a, b| b.1.overall.partial_cmp(&a.1.overall).unwrap());
    results
}

fn main() {
    let manifest = env!("CARGO_MANIFEST_DIR");
    let bam_path = PathBuf::from(manifest).join("tests/data/NA19240.chr17_fragment.bam");
    let ref_path = PathBuf::from(manifest).join("tests/data/chr17_fragment.fa");

    // Full test region
    let region = Region::new("chr17", 10958130, 11017414).unwrap();

    // Load data
    let reference = ReferenceGenome::from_file(&ref_path).expect("failed to read reference");
    let ref_seq = reference.fetch(&region).expect("failed to fetch region");
    let reads = AlignmentReader::read_bam(&bam_path, &region).expect("failed to read BAM");

    let mut out = String::new();

    // === HEADER ===
    writeln!(out, "╔══════════════════════════════════════════════════════════════════════════════╗").unwrap();
    writeln!(out, "║            LONGREAD REVIEWER — Complex Event Analysis Report                 ║").unwrap();
    writeln!(out, "╚══════════════════════════════════════════════════════════════════════════════╝").unwrap();
    writeln!(out).unwrap();

    writeln!(out, "┌─ Input Data ────────────────────────────────────────────────────────────────┐").unwrap();
    writeln!(out, "│ BAM:       tests/data/NA19240.chr17_fragment.bam").unwrap();
    writeln!(out, "│ Reference: tests/data/chr17_fragment.fa").unwrap();
    writeln!(out, "│ Region:    {} ({} bp)", region, region.len()).unwrap();
    writeln!(out, "│ Reads:     {}", reads.len()).unwrap();
    writeln!(out, "└─────────────────────────────────────────────────────────────────────────────┘").unwrap();
    writeln!(out).unwrap();

    // === GLOBAL ASSEMBLY & METRICS ===
    writeln!(out, "┌─ Region-Wide Assembly ────────────────────────────────────────────────────────┐").unwrap();

    let method = ConsensusAssembly;
    let assembly = method.assemble(&reads, &ref_seq, region.start).expect("assembly failed");
    let calc = MetricsCalculator::new();
    let fitness = calc.compute_fitness(&assembly, &reads, &ref_seq, region.start);

    writeln!(out, "│").unwrap();
    writeln!(out, "│ Consensus Statistics:").unwrap();
    writeln!(out, "│   Assembly length:  {} bp", assembly.sequence.len()).unwrap();
    writeln!(out, "│   Mean depth:       {:.1}x", fitness.mean_depth).unwrap();
    writeln!(out, "│   Mean agreement:   {:.1}%", fitness.mean_agreement * 100.0).unwrap();
    writeln!(out, "│   Reference match:  {:.1}%", fitness.reference_identity * 100.0).unwrap();
    writeln!(out, "│   Fitness score:    {:.3}", fitness.overall).unwrap();

    let variant_count = fitness.base_metrics.iter().filter(|m| m.is_variant).count();
    let variant_rate = variant_count as f64 / fitness.base_metrics.len() as f64 * 100.0;
    writeln!(out, "│   Variants:         {} ({:.2}%)", variant_count, variant_rate).unwrap();
    writeln!(out, "│").unwrap();
    writeln!(out, "└─────────────────────────────────────────────────────────────────────────────┘").unwrap();
    writeln!(out).unwrap();

    // === VISUAL TRACKS ===
    let track_width = 78;
    writeln!(out, "┌─ Coverage & Variant Tracks ───────────────────────────────────────────────────┐").unwrap();
    writeln!(out, "│").unwrap();
    writeln!(out, "│ Coverage Profile ({}bp → {} chars):", region.len(), track_width).unwrap();
    for line in render_coverage_track(&fitness, track_width, 5) {
        writeln!(out, "│   {}", line).unwrap();
    }
    writeln!(out, "│").unwrap();
    writeln!(out, "│ Variant Density (█=high ▓=med ▒=low ░=rare ·=none):").unwrap();
    writeln!(out, "│   {}", render_variant_track(&fitness, track_width)).unwrap();
    writeln!(out, "│").unwrap();
    writeln!(out, "│ Read Agreement  (█=95%+ ▓=90%+ ▒=80%+ ░=70%+ ·=<70%):").unwrap();
    writeln!(out, "│   {}", render_agreement_track(&fitness, track_width)).unwrap();
    writeln!(out, "│").unwrap();
    writeln!(out, "│ Position scale:").unwrap();
    writeln!(out, "│   {:.<width$}", format!("{}", region.start), width = track_width).unwrap();
    writeln!(out, "│   {: <width$}{}", "", region.end, width = track_width - 10).unwrap();
    writeln!(out, "│").unwrap();
    writeln!(out, "└─────────────────────────────────────────────────────────────────────────────┘").unwrap();
    writeln!(out).unwrap();

    // === COMPLEX EVENT DETECTION ===
    writeln!(out, "┌─ Complex Event Detection ─────────────────────────────────────────────────────┐").unwrap();
    writeln!(out, "│").unwrap();

    let events = find_complex_events(&fitness, 100, 3, 0.15);

    if events.is_empty() {
        writeln!(out, "│ No complex events detected with current thresholds.").unwrap();
        writeln!(out, "│ (window=100bp, min_variants=3, min_disagreement=15%)").unwrap();
    } else {
        writeln!(out, "│ Found {} potential complex event(s):", events.len()).unwrap();
        writeln!(out, "│").unwrap();
        writeln!(out, "│ {:>4}  {:>12}  {:>12}  {:>6}  {:>8}  {:>10}  {:<20}",
            "#", "Start", "End", "Size", "Variants", "Agreement", "Type").unwrap();
        writeln!(out, "│ {}", "─".repeat(76)).unwrap();

        for (i, event) in events.iter().enumerate() {
            writeln!(out, "│ {:>4}  {:>12}  {:>12}  {:>6}  {:>8}  {:>9.1}%  {:<20}",
                i + 1,
                event.start,
                event.end,
                event.end - event.start + 1,
                event.variant_count,
                event.mean_agreement * 100.0,
                event.event_type
            ).unwrap();
        }
    }
    writeln!(out, "│").unwrap();
    writeln!(out, "└─────────────────────────────────────────────────────────────────────────────┘").unwrap();
    writeln!(out).unwrap();

    // === STRUCTURAL VARIANT DETECTION ===
    writeln!(out, "┌─ Structural Variant Detection (from CIGAR) ─────────────────────────────────┐").unwrap();
    writeln!(out, "│").unwrap();

    // Detect large indels from CIGAR strings
    #[derive(Debug, Clone)]
    struct IndividualSV {
        read_name: String,
        sv_type: String,
        start: u64,
        end: u64,
        size: i64,
    }

    let mut all_svs: Vec<IndividualSV> = Vec::new();

    for read in &reads {
        let mut ref_pos = read.start;
        for op in &read.cigar {
            match op {
                CigarOp::Match(n) => {
                    ref_pos += *n as u64;
                }
                CigarOp::Insertion(n) => {
                    if *n >= 50 {
                        all_svs.push(IndividualSV {
                            read_name: read.name.clone(),
                            sv_type: "INS".to_string(),
                            start: ref_pos,
                            end: ref_pos,
                            size: *n as i64,
                        });
                    }
                }
                CigarOp::Deletion(n) => {
                    if *n >= 50 {
                        all_svs.push(IndividualSV {
                            read_name: read.name.clone(),
                            sv_type: "DEL".to_string(),
                            start: ref_pos,
                            end: ref_pos + *n as u64,
                            size: -(*n as i64),
                        });
                    }
                    ref_pos += *n as u64;
                }
                _ => {}
            }
        }
    }

    // Group SVs by position (within 50bp tolerance) and type
    #[derive(Debug, Clone)]
    struct SVEvent {
        sv_type: String,
        start: u64,
        end: u64,
        size: i64,
        supporting_reads: Vec<String>,
    }

    let mut sv_events: Vec<SVEvent> = Vec::new();

    for sv in &all_svs {
        // Try to find matching event
        let mut found = false;
        for event in &mut sv_events {
            if event.sv_type == sv.sv_type &&
               ((event.start as i64 - sv.start as i64).abs() < 50) &&
               (event.size - sv.size).abs() < 50 {
                // Same event, different read
                if !event.supporting_reads.contains(&sv.read_name) {
                    event.supporting_reads.push(sv.read_name.clone());
                }
                found = true;
                break;
            }
        }
        if !found {
            sv_events.push(SVEvent {
                sv_type: sv.sv_type.clone(),
                start: sv.start,
                end: sv.end,
                size: sv.size,
                supporting_reads: vec![sv.read_name.clone()],
            });
        }
    }

    // Sort by position
    sv_events.sort_by_key(|e| e.start);

    if sv_events.is_empty() {
        writeln!(out, "│ No structural variants (>=50bp) detected in read alignments.").unwrap();
    } else {
        // Calculate totals
        let total_del: i64 = sv_events.iter()
            .filter(|e| e.sv_type == "DEL")
            .map(|e| e.size.abs())
            .sum();
        let total_ins: i64 = sv_events.iter()
            .filter(|e| e.sv_type == "INS")
            .map(|e| e.size.abs())
            .sum();
        let del_count = sv_events.iter().filter(|e| e.sv_type == "DEL").count();
        let ins_count = sv_events.iter().filter(|e| e.sv_type == "INS").count();

        writeln!(out, "│ ⚠ STRUCTURAL VARIANTS DETECTED!").unwrap();
        writeln!(out, "│").unwrap();
        writeln!(out, "│ Summary:").unwrap();
        writeln!(out, "│   Total deletions:   {} bp across {} events", total_del, del_count).unwrap();
        writeln!(out, "│   Total insertions:  {} bp across {} events", total_ins, ins_count).unwrap();
        writeln!(out, "│   Net effect:        {} bp", total_ins - total_del).unwrap();
        writeln!(out, "│").unwrap();

        writeln!(out, "│ {:>4}  {:>5}  {:>12}  {:>12}  {:>8}  {:>6}  {:<25}",
            "#", "Type", "Start", "End", "Size", "Reads", "Supporting Reads").unwrap();
        writeln!(out, "│ {}", "─".repeat(76)).unwrap();

        for (i, event) in sv_events.iter().enumerate() {
            let size_str = if event.size < 0 {
                format!("-{}", event.size.abs())
            } else {
                format!("+{}", event.size)
            };
            let sample_reads: String = event.supporting_reads.iter()
                .take(2)
                .map(|r| {
                    // Extract short read identifier
                    let parts: Vec<&str> = r.split('/').collect();
                    if parts.len() >= 2 {
                        parts[1].to_string()
                    } else {
                        r.chars().take(12).collect()
                    }
                })
                .collect::<Vec<_>>()
                .join(", ");

            writeln!(out, "│ {:>4}  {:>5}  {:>12}  {:>12}  {:>8}  {:>6}  {:<25}",
                i + 1,
                event.sv_type,
                event.start,
                event.end,
                size_str,
                event.supporting_reads.len(),
                sample_reads
            ).unwrap();
        }

        // Identify complex SV region (multiple deletions in close proximity)
        let del_events: Vec<_> = sv_events.iter().filter(|e| e.sv_type == "DEL").collect();
        if del_events.len() >= 3 {
            let first_del = del_events.first().unwrap();
            let last_del = del_events.last().unwrap();
            let span = last_del.end - first_del.start;

            writeln!(out, "│").unwrap();
            writeln!(out, "│ ★ COMPLEX DELETION REGION DETECTED").unwrap();
            writeln!(out, "│   Location:   {}:{}-{}", region.chrom, first_del.start, last_del.end).unwrap();
            writeln!(out, "│   Span:       {} bp", span).unwrap();
            writeln!(out, "│   Total deleted: {} bp across {} events", total_del, del_count).unwrap();
            writeln!(out, "│   Supported by {} unique reads",
                del_events.iter().flat_map(|e| &e.supporting_reads).collect::<std::collections::HashSet<_>>().len()
            ).unwrap();
        }
    }
    writeln!(out, "│").unwrap();
    writeln!(out, "└─────────────────────────────────────────────────────────────────────────────┘").unwrap();
    writeln!(out).unwrap();

    // === ASSEMBLY METHOD COMPARISON ===
    writeln!(out, "┌─ Assembly Method Comparison ──────────────────────────────────────────────────┐").unwrap();
    writeln!(out, "│").unwrap();
    writeln!(out, "│ Evaluating multiple assembly strategies across the full region:").unwrap();
    writeln!(out, "│").unwrap();

    let mut engine = AssemblyEngine::new();
    engine.add_method(Box::new(ConsensusAssembly));
    engine.add_method(Box::new(WindowConsensusAssembly::new(50, 10)));
    engine.add_method(Box::new(WindowConsensusAssembly::new(100, 20)));
    engine.add_method(Box::new(WindowConsensusAssembly::new(200, 40)));
    engine.add_method(Box::new(WindowConsensusAssembly::new(500, 100)));

    let results = engine.evaluate_all(&reads, &ref_seq, &region).expect("evaluation failed");

    writeln!(out, "│ {:<22} {:>10} {:>8} {:>10} {:>10} {:>8}",
        "Method", "Agreement", "Depth", "RefMatch", "Variants", "Fitness").unwrap();
    writeln!(out, "│ {}", "─".repeat(70)).unwrap();

    for result in &results {
        let var_count = result.fitness.base_metrics.iter().filter(|m| m.is_variant).count();
        let name = if result.assembly.method_name == "window_consensus" {
            // Include window parameters in name
            format!("window({})", result.assembly.sequence.len())
        } else {
            result.assembly.method_name.clone()
        };
        writeln!(out, "│ {:<22} {:>9.1}% {:>7.1}x {:>9.1}% {:>10} {:>8.3}",
            name,
            result.fitness.mean_agreement * 100.0,
            result.fitness.mean_depth,
            result.fitness.reference_identity * 100.0,
            var_count,
            result.fitness.overall,
        ).unwrap();
    }

    if let Some(best) = results.first() {
        writeln!(out, "│").unwrap();
        writeln!(out, "│ ★ Best method: {} (fitness: {:.3})", best.assembly.method_name, best.fitness.overall).unwrap();
    }
    writeln!(out, "│").unwrap();
    writeln!(out, "└─────────────────────────────────────────────────────────────────────────────┘").unwrap();
    writeln!(out).unwrap();

    // === DETAILED EVENT ANALYSIS ===
    if !events.is_empty() {
        writeln!(out, "┌─ Detailed Event Analysis ─────────────────────────────────────────────────────┐").unwrap();

        for (i, event) in events.iter().take(3).enumerate() {
            writeln!(out, "│").unwrap();
            writeln!(out, "│ Event #{}: {} ({}:{}-{})",
                i + 1, event.event_type, region.chrom, event.start, event.end).unwrap();
            writeln!(out, "│ {}", "─".repeat(70)).unwrap();

            let event_assemblies = compare_assemblies_at_event(event, &reads, &reference, &region);

            if event_assemblies.is_empty() {
                writeln!(out, "│   No assemblies generated (insufficient data)").unwrap();
            } else {
                writeln!(out, "│   Assembly comparison at this locus:").unwrap();
                writeln!(out, "│   {:<18} {:>10} {:>10} {:>12}", "Method", "Agreement", "RefMatch", "Fitness").unwrap();
                writeln!(out, "│   {}", "─".repeat(52)).unwrap();

                for (name, fit, _asm) in &event_assemblies {
                    writeln!(out, "│   {:<18} {:>9.1}% {:>9.1}% {:>12.3}",
                        name,
                        fit.mean_agreement * 100.0,
                        fit.reference_identity * 100.0,
                        fit.overall
                    ).unwrap();
                }

                // Show sequence differences for top 2 methods
                if event_assemblies.len() >= 2 {
                    let (name1, _, asm1) = &event_assemblies[0];
                    let (name2, _, asm2) = &event_assemblies[1];

                    let diff_count = asm1.sequence.iter()
                        .zip(asm2.sequence.iter())
                        .filter(|(a, b)| !a.eq_ignore_ascii_case(b))
                        .count();

                    writeln!(out, "│").unwrap();
                    writeln!(out, "│   Sequence differences between {} and {}: {} bp", name1, name2, diff_count).unwrap();
                }
            }
        }
        writeln!(out, "│").unwrap();
        writeln!(out, "└─────────────────────────────────────────────────────────────────────────────┘").unwrap();
        writeln!(out).unwrap();
    }

    // === HAPLOTYPE ANALYSIS ===
    writeln!(out, "┌─ Haplotype Analysis ──────────────────────────────────────────────────────────┐").unwrap();
    writeln!(out, "│").unwrap();

    let assigner = HaplotypeAssigner::new();
    let assignments = assigner.assign(&reads, &ref_seq);

    let hap1: Vec<_> = assignments.iter().filter(|a| a.haplotype == HaplotypeLabel::Hap1).collect();
    let hap2: Vec<_> = assignments.iter().filter(|a| a.haplotype == HaplotypeLabel::Hap2).collect();
    let unassigned: Vec<_> = assignments.iter().filter(|a| a.haplotype == HaplotypeLabel::Unassigned).collect();

    writeln!(out, "│ Haplotype Distribution:").unwrap();
    writeln!(out, "│   H1:         {:>3} reads ({:.1}%)", hap1.len(), hap1.len() as f64 / reads.len().max(1) as f64 * 100.0).unwrap();
    writeln!(out, "│   H2:         {:>3} reads ({:.1}%)", hap2.len(), hap2.len() as f64 / reads.len().max(1) as f64 * 100.0).unwrap();
    writeln!(out, "│   Unassigned: {:>3} reads ({:.1}%)", unassigned.len(), unassigned.len() as f64 / reads.len().max(1) as f64 * 100.0).unwrap();
    writeln!(out, "│").unwrap();

    // Show high-confidence assignments
    let mut high_conf: Vec<_> = assignments.iter().filter(|a| a.confidence > 0.7).collect();
    high_conf.sort_by(|a, b| b.confidence.partial_cmp(&a.confidence).unwrap());

    if !high_conf.is_empty() {
        writeln!(out, "│ High-confidence assignments (>70%):").unwrap();
        writeln!(out, "│   {:<40} {:>10} {:>12}", "Read", "Haplotype", "Confidence").unwrap();
        writeln!(out, "│   {}", "─".repeat(64)).unwrap();
        for a in high_conf.iter().take(10) {
            let short_name = if a.read_name.len() > 38 {
                format!("{}...", &a.read_name[..35])
            } else {
                a.read_name.clone()
            };
            writeln!(out, "│   {:<40} {:>10} {:>11.1}%", short_name, a.haplotype, a.confidence * 100.0).unwrap();
        }
        if high_conf.len() > 10 {
            writeln!(out, "│   ... and {} more", high_conf.len() - 10).unwrap();
        }
    }
    writeln!(out, "│").unwrap();
    writeln!(out, "└─────────────────────────────────────────────────────────────────────────────┘").unwrap();
    writeln!(out).unwrap();

    // === READ SUMMARY ===
    writeln!(out, "┌─ Read Summary ─────────────────────────────────────────────────────────────────┐").unwrap();
    writeln!(out, "│").unwrap();
    writeln!(out, "│ {:<36} {:>10} {:>10} {:>6} {:>4} {:>4}",
        "Name", "Start", "End", "Len", "MQ", "HP").unwrap();
    writeln!(out, "│ {}", "─".repeat(74)).unwrap();

    for read in reads.iter().take(20) {
        let hp_str = match read.haplotype_tag {
            Some(h) => format!("{}", h),
            None => "-".to_string(),
        };
        let short_name = if read.name.len() > 34 {
            format!("{}...", &read.name[..31])
        } else {
            read.name.clone()
        };
        writeln!(out, "│ {:<36} {:>10} {:>10} {:>6} {:>4} {:>4}",
            short_name, read.start, read.end, read.sequence.len(), read.mapq, hp_str
        ).unwrap();
    }
    if reads.len() > 20 {
        writeln!(out, "│ ... and {} more reads", reads.len() - 20).unwrap();
    }
    writeln!(out, "│").unwrap();
    writeln!(out, "└─────────────────────────────────────────────────────────────────────────────┘").unwrap();

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
