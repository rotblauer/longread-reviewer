use ratatui::buffer::Buffer;
use ratatui::layout::Rect;
use ratatui::style::{Color, Modifier, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, Paragraph, Widget};

use crate::alignment::AlignedRead;
use crate::assembly::method::AssemblyResult;
use crate::haplotype::HaplotypeLabel;
use crate::metrics::BaseMetric;

/// Color a base by its nucleotide identity (IGV-style coloring).
pub fn base_color(base: u8) -> Color {
    match base.to_ascii_uppercase() {
        b'A' => Color::Green,
        b'C' => Color::Blue,
        b'G' => Color::Yellow,
        b'T' => Color::Red,
        b'N' => Color::DarkGray,
        _ => Color::White,
    }
}

/// Color for haplotype labels.
pub fn haplotype_color(hap: HaplotypeLabel) -> Color {
    match hap {
        HaplotypeLabel::Hap1 => Color::Cyan,
        HaplotypeLabel::Hap2 => Color::Magenta,
        HaplotypeLabel::Unassigned => Color::DarkGray,
    }
}

/// Render the reference track as a line of colored bases.
pub fn render_reference_line(reference: &[u8], offset: usize, width: usize) -> Line<'static> {
    if offset >= reference.len() {
        return Line::from(Vec::new());
    }
    let end = std::cmp::min(offset + width, reference.len());
    let spans: Vec<Span> = reference[offset..end]
        .iter()
        .map(|&b| {
            Span::styled(
                String::from(b as char),
                Style::default()
                    .fg(base_color(b))
                    .add_modifier(Modifier::BOLD),
            )
        })
        .collect();
    Line::from(spans)
}

/// Render a single read aligned against the reference coordinate system.
pub fn render_read_line(
    read: &AlignedRead,
    ref_start: u64,
    offset: usize,
    width: usize,
    haplotype: Option<HaplotypeLabel>,
    reference: &[u8],
) -> Line<'static> {
    let view_start = ref_start + offset as u64;
    let view_end = view_start + width as u64;

    let mut spans: Vec<Span> = Vec::with_capacity(width);
    let aligned_bases = read.reference_aligned_bases();

    for pos in view_start..view_end {
        if pos < read.start || pos > read.end {
            spans.push(Span::raw(" "));
        } else if let Some(&(_, base)) = aligned_bases.iter().find(|&&(rp, _)| rp == pos) {
            let ref_idx = (pos - ref_start) as usize;
            let ref_base = if ref_idx < reference.len() {
                reference[ref_idx]
            } else {
                b'N'
            };

            let style = if !base.eq_ignore_ascii_case(&ref_base) {
                // Mismatch: highlight with inverse colors
                Style::default()
                    .fg(Color::White)
                    .bg(base_color(base))
                    .add_modifier(Modifier::BOLD)
            } else {
                // Match: dim the base
                let fg = match haplotype {
                    Some(h) => haplotype_color(h),
                    None => Color::Gray,
                };
                Style::default().fg(fg)
            };

            spans.push(Span::styled(String::from(base as char), style));
        } else {
            // Deletion
            spans.push(Span::styled(
                "-".to_string(),
                Style::default().fg(Color::Red),
            ));
        }
    }

    Line::from(spans)
}

/// Render assembly consensus line with confidence coloring.
pub fn render_assembly_line(
    assembly: &AssemblyResult,
    offset: usize,
    width: usize,
    reference: &[u8],
) -> Line<'static> {
    let end = std::cmp::min(offset + width, assembly.sequence.len());
    let spans: Vec<Span> = (offset..end)
        .map(|i| {
            let base = assembly.sequence[i];
            let conf = assembly.confidence[i];
            let ref_base = if i < reference.len() {
                reference[i]
            } else {
                b'N'
            };

            let is_variant = !base.eq_ignore_ascii_case(&ref_base);

            let fg = if is_variant {
                Color::White
            } else if conf > 0.9 {
                base_color(base)
            } else if conf > 0.5 {
                Color::Yellow
            } else {
                Color::Red
            };

            let bg = if is_variant {
                Color::Magenta
            } else {
                Color::Reset
            };

            let style = Style::default().fg(fg).bg(bg);
            Span::styled(String::from(base as char), style)
        })
        .collect();
    Line::from(spans)
}

/// Render a metrics bar showing per-base agreement as a colored bar.
pub fn render_metrics_line(
    metrics: &[BaseMetric],
    offset: usize,
    width: usize,
) -> Line<'static> {
    let end = std::cmp::min(offset + width, metrics.len());
    let spans: Vec<Span> = metrics[offset..end]
        .iter()
        .map(|m| {
            let ch = if m.agreement > 0.9 {
                '█'
            } else if m.agreement > 0.7 {
                '▓'
            } else if m.agreement > 0.5 {
                '▒'
            } else if m.agreement > 0.3 {
                '░'
            } else {
                ' '
            };

            let color = if m.agreement > 0.9 {
                Color::Green
            } else if m.agreement > 0.7 {
                Color::Yellow
            } else if m.agreement > 0.5 {
                Color::Red
            } else {
                Color::DarkGray
            };

            Span::styled(ch.to_string(), Style::default().fg(color))
        })
        .collect();
    Line::from(spans)
}

/// A widget that renders the full pileup view (reference + reads + assembly + metrics).
pub struct PileupView<'a> {
    pub reference: &'a [u8],
    pub reads: &'a [AlignedRead],
    pub assembly: Option<&'a AssemblyResult>,
    pub metrics: Option<&'a [BaseMetric]>,
    pub haplotype_assignments: &'a [(String, HaplotypeLabel)],
    pub ref_start: u64,
    pub offset: usize,
    pub scroll_y: usize,
    pub region_label: &'a str,
}

impl Widget for PileupView<'_> {
    fn render(self, area: Rect, buf: &mut Buffer) {
        let block = Block::default()
            .title(format!(" {} ", self.region_label))
            .borders(Borders::ALL);
        let inner = block.inner(area);
        block.render(area, buf);

        let width = inner.width as usize;
        let mut y = inner.y;

        // Coordinate ruler
        if y < inner.y + inner.height {
            let ruler = render_ruler(self.ref_start, self.offset, width);
            buf.set_line(inner.x, y, &ruler, inner.width);
            y += 1;
        }

        // Reference track
        if y < inner.y + inner.height {
            let ref_line = render_reference_line(self.reference, self.offset, width);
            buf.set_line(inner.x, y, &ref_line, inner.width);
            y += 1;
        }

        // Assembly track (if present)
        if let Some(assembly) = self.assembly
            && y < inner.y + inner.height {
                let asm_line = render_assembly_line(assembly, self.offset, width, self.reference);
                buf.set_line(inner.x, y, &asm_line, inner.width);
                y += 1;
            }

        // Metrics bar (if present)
        if let Some(metrics) = self.metrics
            && y < inner.y + inner.height {
                let met_line = render_metrics_line(metrics, self.offset, width);
                buf.set_line(inner.x, y, &met_line, inner.width);
                y += 1;
            }

        // Separator
        if y < inner.y + inner.height {
            let sep = Line::from("─".repeat(width));
            buf.set_line(inner.x, y, &sep, inner.width);
            y += 1;
        }

        // Reads
        let reads_to_show = self
            .reads
            .iter()
            .skip(self.scroll_y)
            .take((inner.y + inner.height - y) as usize);

        for read in reads_to_show {
            if y >= inner.y + inner.height {
                break;
            }

            let hap = self
                .haplotype_assignments
                .iter()
                .find(|(name, _)| name == &read.name)
                .map(|(_, h)| *h);

            let read_line =
                render_read_line(read, self.ref_start, self.offset, width, hap, self.reference);
            buf.set_line(inner.x, y, &read_line, inner.width);
            y += 1;
        }
    }
}

/// Render a coordinate ruler line.
fn render_ruler(ref_start: u64, offset: usize, width: usize) -> Line<'static> {
    let start_pos = ref_start + offset as u64;
    let mut ruler = String::with_capacity(width);

    for i in 0..width {
        let pos = start_pos + i as u64;
        if pos.is_multiple_of(10) {
            let label = format!("{pos}");
            if i + label.len() <= width {
                ruler.push('|');
            } else {
                ruler.push(' ');
            }
        } else {
            let prev_tick = pos - (pos % 10);
            let label = format!("{prev_tick}");
            let label_offset = (pos - prev_tick) as usize;
            if label_offset < label.len() && label_offset > 0 {
                ruler.push(label.as_bytes()[label_offset] as char);
            } else {
                ruler.push('·');
            }
        }
    }

    Line::from(Span::styled(
        ruler,
        Style::default().fg(Color::DarkGray),
    ))
}

/// Render a status bar with summary information.
pub fn render_status_bar(
    region: &str,
    num_reads: usize,
    assembly_method: Option<&str>,
    fitness_score: Option<f64>,
) -> Paragraph<'static> {
    let mut parts = vec![
        Span::styled(
            format!(" Region: {region} "),
            Style::default()
                .fg(Color::White)
                .bg(Color::DarkGray),
        ),
        Span::styled(
            format!(" Reads: {num_reads} "),
            Style::default()
                .fg(Color::White)
                .bg(Color::DarkGray),
        ),
    ];

    if let Some(method) = assembly_method {
        parts.push(Span::styled(
            format!(" Assembly: {method} "),
            Style::default()
                .fg(Color::Cyan)
                .bg(Color::DarkGray),
        ));
    }

    if let Some(score) = fitness_score {
        let color = if score > 0.8 {
            Color::Green
        } else if score > 0.5 {
            Color::Yellow
        } else {
            Color::Red
        };
        parts.push(Span::styled(
            format!(" Fitness: {score:.3} "),
            Style::default().fg(color).bg(Color::DarkGray),
        ));
    }

    parts.push(Span::styled(
        " [q]uit [a]ssemble [h]aplotype [←→]scroll [↑↓]reads [n]ext method ".to_string(),
        Style::default()
            .fg(Color::Yellow)
            .bg(Color::DarkGray),
    ));

    Paragraph::new(Line::from(parts))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_base_color() {
        assert_eq!(base_color(b'A'), Color::Green);
        assert_eq!(base_color(b'C'), Color::Blue);
        assert_eq!(base_color(b'G'), Color::Yellow);
        assert_eq!(base_color(b'T'), Color::Red);
        assert_eq!(base_color(b'N'), Color::DarkGray);
    }

    #[test]
    fn test_render_reference_line() {
        let reference = b"ACGTACGT";
        let line = render_reference_line(reference, 0, 4);
        assert_eq!(line.spans.len(), 4);
    }

    #[test]
    fn test_render_reference_line_offset() {
        let reference = b"ACGTACGT";
        let line = render_reference_line(reference, 4, 4);
        assert_eq!(line.spans.len(), 4);
    }

    #[test]
    fn test_render_metrics_line() {
        let metrics = vec![
            BaseMetric {
                position: 1,
                ref_base: b'A',
                assembly_base: b'A',
                depth: 10,
                agreement: 1.0,
                avg_quality: 30.0,
                is_variant: false,
            },
            BaseMetric {
                position: 2,
                ref_base: b'C',
                assembly_base: b'T',
                depth: 10,
                agreement: 0.4,
                avg_quality: 30.0,
                is_variant: true,
            },
        ];

        let line = render_metrics_line(&metrics, 0, 2);
        assert_eq!(line.spans.len(), 2);
    }

    #[test]
    fn test_base_color_lowercase() {
        // Should handle lowercase bases too
        assert_eq!(base_color(b'a'), Color::Green);
        assert_eq!(base_color(b'c'), Color::Blue);
        assert_eq!(base_color(b'g'), Color::Yellow);
        assert_eq!(base_color(b't'), Color::Red);
    }

    #[test]
    fn test_base_color_unknown() {
        assert_eq!(base_color(b'X'), Color::White);
        assert_eq!(base_color(b'-'), Color::White);
    }

    #[test]
    fn test_haplotype_color() {
        assert_eq!(haplotype_color(HaplotypeLabel::Hap1), Color::Cyan);
        assert_eq!(haplotype_color(HaplotypeLabel::Hap2), Color::Magenta);
        assert_eq!(haplotype_color(HaplotypeLabel::Unassigned), Color::DarkGray);
    }

    #[test]
    fn test_render_reference_line_empty() {
        let reference: &[u8] = b"";
        let line = render_reference_line(reference, 0, 10);
        assert_eq!(line.spans.len(), 0);
    }

    #[test]
    fn test_render_reference_line_wider_than_ref() {
        let reference = b"ACGT";
        let line = render_reference_line(reference, 0, 100);
        assert_eq!(line.spans.len(), 4);
    }

    #[test]
    fn test_render_reference_line_offset_beyond() {
        let reference = b"ACGT";
        let line = render_reference_line(reference, 10, 4);
        assert_eq!(line.spans.len(), 0);
    }

    #[test]
    fn test_render_status_bar_minimal() {
        let status = render_status_bar("chr1:1-100", 10, None, None);
        // Should not panic
        let _ = status;
    }

    #[test]
    fn test_render_status_bar_full() {
        let status = render_status_bar("chr1:1-100", 42, Some("consensus"), Some(0.95));
        let _ = status;
    }

    #[test]
    fn test_render_status_bar_low_fitness() {
        let status = render_status_bar("chr1:1-100", 5, Some("window_consensus"), Some(0.3));
        let _ = status;
    }

    #[test]
    fn test_render_assembly_line() {
        let assembly = AssemblyResult {
            sequence: b"ACGT".to_vec(),
            depth: vec![10, 10, 10, 10],
            confidence: vec![1.0, 0.8, 0.4, 0.1],
            method_name: "test".to_string(),
        };
        let reference = b"ACGT";

        let line = render_assembly_line(&assembly, 0, 4, reference);
        assert_eq!(line.spans.len(), 4);
    }

    #[test]
    fn test_render_assembly_line_with_variant() {
        let assembly = AssemblyResult {
            sequence: b"ATGT".to_vec(), // C->T variant at position 2
            depth: vec![10, 10, 10, 10],
            confidence: vec![1.0, 0.8, 1.0, 1.0],
            method_name: "test".to_string(),
        };
        let reference = b"ACGT";

        let line = render_assembly_line(&assembly, 0, 4, reference);
        assert_eq!(line.spans.len(), 4);
    }

    #[test]
    fn test_render_assembly_line_offset() {
        let assembly = AssemblyResult {
            sequence: b"ACGTACGT".to_vec(),
            depth: vec![10; 8],
            confidence: vec![1.0; 8],
            method_name: "test".to_string(),
        };
        let reference = b"ACGTACGT";

        let line = render_assembly_line(&assembly, 4, 4, reference);
        assert_eq!(line.spans.len(), 4);
    }

    #[test]
    fn test_render_read_line_basic() {
        use crate::alignment::CigarOp;

        let read = AlignedRead {
            name: "read1".to_string(),
            chrom: "chr1".to_string(),
            start: 1,
            mapq: 60,
            cigar: vec![CigarOp::Match(8)],
            sequence: b"ACGTACGT".to_vec(),
            qualities: vec![30; 8],
            is_reverse: false,
            haplotype_tag: None,
            end: 8,
        };
        let reference = b"ACGTACGT";

        let line = render_read_line(&read, 1, 0, 8, None, reference);
        assert_eq!(line.spans.len(), 8);
    }

    #[test]
    fn test_render_read_line_with_haplotype() {
        use crate::alignment::CigarOp;

        let read = AlignedRead {
            name: "read1".to_string(),
            chrom: "chr1".to_string(),
            start: 1,
            mapq: 60,
            cigar: vec![CigarOp::Match(4)],
            sequence: b"ACGT".to_vec(),
            qualities: vec![30; 4],
            is_reverse: false,
            haplotype_tag: None,
            end: 4,
        };
        let reference = b"ACGT";

        let line = render_read_line(&read, 1, 0, 4, Some(HaplotypeLabel::Hap1), reference);
        assert_eq!(line.spans.len(), 4);
    }

    #[test]
    fn test_render_read_line_with_mismatch() {
        use crate::alignment::CigarOp;

        let read = AlignedRead {
            name: "read1".to_string(),
            chrom: "chr1".to_string(),
            start: 1,
            mapq: 60,
            cigar: vec![CigarOp::Match(4)],
            sequence: b"ATGT".to_vec(), // C->T mismatch at pos 2
            qualities: vec![30; 4],
            is_reverse: false,
            haplotype_tag: None,
            end: 4,
        };
        let reference = b"ACGT";

        let line = render_read_line(&read, 1, 0, 4, None, reference);
        assert_eq!(line.spans.len(), 4);
    }

    #[test]
    fn test_render_read_line_partial_overlap() {
        use crate::alignment::CigarOp;

        // Read starts at 5, but viewport starts at 1
        let read = AlignedRead {
            name: "read1".to_string(),
            chrom: "chr1".to_string(),
            start: 5,
            mapq: 60,
            cigar: vec![CigarOp::Match(4)],
            sequence: b"ACGT".to_vec(),
            qualities: vec![30; 4],
            is_reverse: false,
            haplotype_tag: None,
            end: 8,
        };
        let reference = b"ACGTACGT";

        let line = render_read_line(&read, 1, 0, 8, None, reference);
        assert_eq!(line.spans.len(), 8);
    }

    #[test]
    fn test_render_metrics_line_various_agreements() {
        let metrics = vec![
            BaseMetric {
                position: 1, ref_base: b'A', assembly_base: b'A',
                depth: 10, agreement: 0.95, avg_quality: 30.0, is_variant: false,
            },
            BaseMetric {
                position: 2, ref_base: b'C', assembly_base: b'C',
                depth: 10, agreement: 0.75, avg_quality: 30.0, is_variant: false,
            },
            BaseMetric {
                position: 3, ref_base: b'G', assembly_base: b'G',
                depth: 10, agreement: 0.55, avg_quality: 30.0, is_variant: false,
            },
            BaseMetric {
                position: 4, ref_base: b'T', assembly_base: b'T',
                depth: 10, agreement: 0.35, avg_quality: 30.0, is_variant: false,
            },
            BaseMetric {
                position: 5, ref_base: b'A', assembly_base: b'A',
                depth: 10, agreement: 0.1, avg_quality: 30.0, is_variant: false,
            },
        ];

        let line = render_metrics_line(&metrics, 0, 5);
        assert_eq!(line.spans.len(), 5);
    }

    #[test]
    fn test_render_metrics_line_with_offset() {
        let metrics = vec![
            BaseMetric {
                position: 1, ref_base: b'A', assembly_base: b'A',
                depth: 10, agreement: 1.0, avg_quality: 30.0, is_variant: false,
            },
            BaseMetric {
                position: 2, ref_base: b'C', assembly_base: b'C',
                depth: 10, agreement: 0.8, avg_quality: 30.0, is_variant: false,
            },
            BaseMetric {
                position: 3, ref_base: b'G', assembly_base: b'G',
                depth: 10, agreement: 0.6, avg_quality: 30.0, is_variant: false,
            },
        ];

        let line = render_metrics_line(&metrics, 1, 2);
        assert_eq!(line.spans.len(), 2);
    }
}
