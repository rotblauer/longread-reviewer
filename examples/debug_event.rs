//! Debug script to analyze the complex event read
//!
//! Usage: cargo run --example debug_event

use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

use longread_reviewer::alignment::AlignmentReader;
use longread_reviewer::reference::ReferenceGenome;
use longread_reviewer::region::Region;

fn main() {
    let manifest = env!("CARGO_MANIFEST_DIR");
    let bam_path = PathBuf::from(manifest).join("tests/data/NA19240.chr17_fragment.bam");
    let ref_path = PathBuf::from(manifest).join("tests/data/chr17_fragment.fa");
    let output_path = PathBuf::from(manifest).join("examples/debug_event_output.txt");

    let mut out = File::create(&output_path).expect("failed to create output file");

    let region = Region::new("chr17", 10958130, 11017414).unwrap();

    let reference = ReferenceGenome::from_file(&ref_path).expect("failed to read reference");
    let ref_seq = reference.fetch(&region).expect("failed to fetch region");
    let reads = AlignmentReader::read_bam(&bam_path, &region).expect("failed to read BAM");

    writeln!(out, "Region: {} ({} bp)", region, region.len()).unwrap();
    writeln!(out, "Reference length: {} bp", ref_seq.len()).unwrap();
    writeln!(out, "Total reads: {}", reads.len()).unwrap();
    writeln!(out).unwrap();

    // Find the read of interest
    let target_read = "m64043_200201_000449/3932797/ccs";

    for read in &reads {
        if read.name == target_read || read.name.contains("3932797") {
            writeln!(out, "=== Found target read: {} ===", read.name).unwrap();
            writeln!(out, "  Start:    {}", read.start).unwrap();
            writeln!(out, "  End:      {}", read.end).unwrap();
            writeln!(out, "  Length:   {} bp", read.sequence.len()).unwrap();
            writeln!(out, "  Span:     {} bp (on reference)", read.end - read.start + 1).unwrap();
            writeln!(out, "  MAPQ:     {}", read.mapq).unwrap();
            writeln!(out, "  Reverse:  {}", read.is_reverse).unwrap();
            writeln!(out).unwrap();

            // Analyze CIGAR for large indels
            writeln!(out, "  CIGAR analysis:").unwrap();
            let mut ref_pos = read.start;
            let mut read_pos = 0u64;
            let mut total_deletions = 0u64;
            let mut total_insertions = 0u64;
            let mut large_events = Vec::new();

            for op in &read.cigar {
                match op {
                    longread_reviewer::alignment::CigarOp::Match(n) => {
                        ref_pos += *n as u64;
                        read_pos += *n as u64;
                    }
                    longread_reviewer::alignment::CigarOp::Insertion(n) => {
                        total_insertions += *n as u64;
                        if *n >= 50 {
                            large_events.push(format!(
                                "    INSERTION: {} bp at ref_pos {} (read_pos {})",
                                n, ref_pos, read_pos
                            ));
                        }
                        read_pos += *n as u64;
                    }
                    longread_reviewer::alignment::CigarOp::Deletion(n) => {
                        total_deletions += *n as u64;
                        if *n >= 50 {
                            large_events.push(format!(
                                "    DELETION:  {} bp at ref_pos {}-{} (read_pos {})",
                                n, ref_pos, ref_pos + *n as u64, read_pos
                            ));
                        }
                        ref_pos += *n as u64;
                    }
                    longread_reviewer::alignment::CigarOp::SoftClip(n) => {
                        read_pos += *n as u64;
                    }
                    longread_reviewer::alignment::CigarOp::HardClip(_) => {}
                }
            }

            writeln!(out, "    Total insertions: {} bp", total_insertions).unwrap();
            writeln!(out, "    Total deletions:  {} bp", total_deletions).unwrap();
            writeln!(out, "    Net indel:        {} bp", total_insertions as i64 - total_deletions as i64).unwrap();
            writeln!(out).unwrap();

            if !large_events.is_empty() {
                writeln!(out, "  Large events (>=50bp):").unwrap();
                for event in &large_events {
                    writeln!(out, "{}", event).unwrap();
                }
            } else {
                writeln!(out, "  No large events (>=50bp) found in CIGAR").unwrap();
            }
            writeln!(out).unwrap();
        }
    }

    // Also look for any reads with large indels
    writeln!(out, "=== All reads with large indels (>=100bp) ===").unwrap();
    for read in &reads {
        let mut has_large_indel = false;
        let mut events = Vec::new();
        let mut ref_pos = read.start;

        for op in &read.cigar {
            match op {
                longread_reviewer::alignment::CigarOp::Match(n) => {
                    ref_pos += *n as u64;
                }
                longread_reviewer::alignment::CigarOp::Insertion(n) => {
                    if *n >= 100 {
                        has_large_indel = true;
                        events.push(format!("INS:{}bp@{}", n, ref_pos));
                    }
                }
                longread_reviewer::alignment::CigarOp::Deletion(n) => {
                    if *n >= 100 {
                        has_large_indel = true;
                        events.push(format!("DEL:{}bp@{}-{}", n, ref_pos, ref_pos + *n as u64));
                    }
                    ref_pos += *n as u64;
                }
                longread_reviewer::alignment::CigarOp::SoftClip(_) => {}
                longread_reviewer::alignment::CigarOp::HardClip(_) => {}
            }
        }

        if has_large_indel {
            writeln!(out, "  {} ({}:{}-{}): {}",
                read.name, read.chrom, read.start, read.end,
                events.join(", ")
            ).unwrap();
        }
    }

    eprintln!("Wrote output to {}", output_path.display());
}


