//! Enhanced IGV-like viewer with structural variant detection and local assembly comparison.

use std::collections::HashMap;

use ratatui::buffer::Buffer;
use ratatui::layout::Rect;
use ratatui::style::{Color, Modifier, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, Widget};

use crate::alignment::{AlignedRead, CigarOp};
use crate::assembly::method::AssemblyResult;
use crate::haplotype::HaplotypeLabel;
use crate::metrics::BaseMetric;

/// Detected structural variant for display.
#[derive(Debug, Clone)]
pub struct SVEvent {
    pub sv_type: String,  // "DEL" or "INS"
    pub start: u64,
    pub end: u64,
    pub size: i64,
    pub supporting_reads: Vec<String>,
}

/// View mode for the enhanced viewer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ViewMode {
    /// Standard pileup view
    Pileup,
    /// Focus on a specific SV event
    SVDetail,
    /// Side-by-side haplotype comparison
    HaplotypeCompare,
    /// Read detail view
    ReadDetail,
}

/// Enhanced pileup view with SV awareness.
pub struct EnhancedPileupView<'a> {
    pub reference: &'a [u8],
    pub reads: &'a [AlignedRead],
    pub assembly: Option<&'a AssemblyResult>,
    pub metrics: Option<&'a [BaseMetric]>,
    pub haplotype_assignments: &'a [(String, HaplotypeLabel)],
    pub sv_events: &'a [SVEvent],
    pub ref_start: u64,
    pub offset: usize,
    pub scroll_y: usize,
    pub selected_read: Option<usize>,
    pub selected_sv: Option<usize>,
    pub view_mode: ViewMode,
    pub show_insertions: bool,
    pub show_soft_clips: bool,
    pub highlight_mismatches: bool,
    pub region_label: &'a str,
}

impl<'a> EnhancedPileupView<'a> {
    /// Render a single read with full CIGAR awareness (insertions, deletions, soft clips).
    fn render_read_with_cigar(
        &self,
        read: &AlignedRead,
        width: usize,
        haplotype: Option<HaplotypeLabel>,
    ) -> Vec<Span<'static>> {
        let view_start = self.ref_start + self.offset as u64;
        let view_end = view_start + width as u64;

        let mut spans: Vec<Span<'static>> = Vec::new();
        let mut ref_pos = read.start;
        let mut read_pos = 0usize;

        // Pre-compute what we'll show at each reference position
        let mut position_content: HashMap<u64, (u8, bool, bool)> = HashMap::new(); // pos -> (base, is_mismatch, is_insertion_before)
        let mut insertions_at: HashMap<u64, Vec<u8>> = HashMap::new();
        let mut deletions: Vec<(u64, u64)> = Vec::new(); // (start, end)

        for op in &read.cigar {
            match op {
                CigarOp::Match(n) => {
                    for _ in 0..*n {
                        if read_pos < read.sequence.len() {
                            let base = read.sequence[read_pos];
                            let ref_idx = (ref_pos.saturating_sub(self.ref_start)) as usize;
                            let ref_base = if ref_idx < self.reference.len() {
                                self.reference[ref_idx]
                            } else {
                                b'N'
                            };
                            let is_mismatch = !base.eq_ignore_ascii_case(&ref_base);
                            position_content.insert(ref_pos, (base, is_mismatch, false));
                        }
                        ref_pos += 1;
                        read_pos += 1;
                    }
                }
                CigarOp::Insertion(n) => {
                    let ins_bases: Vec<u8> = read.sequence[read_pos..read_pos + *n as usize].to_vec();
                    insertions_at.entry(ref_pos).or_default().extend(ins_bases);
                    read_pos += *n as usize;
                }
                CigarOp::Deletion(n) => {
                    deletions.push((ref_pos, ref_pos + *n as u64));
                    ref_pos += *n as u64;
                }
                CigarOp::SoftClip(n) => {
                    read_pos += *n as usize;
                }
                CigarOp::HardClip(_) => {}
            }
        }

        // Now render the view
        for pos in view_start..view_end {
            if pos < read.start || pos > read.end {
                spans.push(Span::raw(" "));
                continue;
            }

            // Check if there's an insertion to show before this position
            if self.show_insertions {
                if let Some(ins_bases) = insertions_at.get(&pos) {
                    // Show insertion marker
                    if ins_bases.len() <= 3 {
                        for &b in ins_bases {
                            spans.push(Span::styled(
                                String::from(b as char),
                                Style::default().fg(Color::Magenta).add_modifier(Modifier::BOLD),
                            ));
                        }
                    } else {
                        spans.push(Span::styled(
                            format!("[+{}]", ins_bases.len()),
                            Style::default().fg(Color::Magenta).add_modifier(Modifier::BOLD),
                        ));
                    }
                }
            }

            // Check if this position is in a deletion
            let in_deletion = deletions.iter().any(|(s, e)| pos >= *s && pos < *e);

            if in_deletion {
                spans.push(Span::styled(
                    "─".to_string(),
                    Style::default().fg(Color::Red).add_modifier(Modifier::DIM),
                ));
            } else if let Some((base, is_mismatch, _)) = position_content.get(&pos) {
                let style = if *is_mismatch && self.highlight_mismatches {
                    Style::default()
                        .fg(Color::White)
                        .bg(base_color(*base))
                        .add_modifier(Modifier::BOLD)
                } else {
                    let fg = match haplotype {
                        Some(HaplotypeLabel::Hap1) => Color::Cyan,
                        Some(HaplotypeLabel::Hap2) => Color::Magenta,
                        _ => Color::Gray,
                    };
                    Style::default().fg(fg)
                };
                spans.push(Span::styled(String::from(*base as char), style));
            } else {
                spans.push(Span::raw(" "));
            }
        }

        spans
    }

    /// Render the SV track showing locations of structural variants.
    fn render_sv_track(&self, width: usize) -> Line<'static> {
        let view_start = self.ref_start + self.offset as u64;
        let view_end = view_start + width as u64;

        let mut track = vec![' '; width];

        for (idx, sv) in self.sv_events.iter().enumerate() {
            let is_selected = self.selected_sv == Some(idx);

            if sv.end >= view_start && sv.start <= view_end {
                let start_idx = if sv.start > view_start {
                    (sv.start - view_start) as usize
                } else {
                    0
                };
                let end_idx = std::cmp::min(
                    if sv.end > view_start { (sv.end - view_start) as usize } else { 0 },
                    width
                );

                let ch = if sv.sv_type == "DEL" { '▼' } else { '▲' };

                for i in start_idx..end_idx {
                    if i < track.len() {
                        track[i] = ch;
                    }
                }

                // Mark the start with a different character
                if start_idx < track.len() {
                    track[start_idx] = if is_selected { '█' } else { '┃' };
                }
            }
        }

        let spans: Vec<Span> = track.iter().map(|&ch| {
            let color = match ch {
                '▼' => Color::Red,     // Deletion
                '▲' => Color::Green,   // Insertion
                '█' | '┃' => Color::Yellow,  // Start marker
                _ => Color::DarkGray,
            };
            Span::styled(ch.to_string(), Style::default().fg(color))
        }).collect();

        Line::from(spans)
    }

    /// Render the depth track as a sparkline.
    fn render_depth_track(&self, width: usize) -> Line<'static> {
        let view_start = self.ref_start + self.offset as u64;

        let depths: Vec<u32> = (0..width).map(|i| {
            let pos = view_start + i as u64;
            self.reads.iter()
                .filter(|r| r.start <= pos && r.end >= pos)
                .count() as u32
        }).collect();

        let max_depth = *depths.iter().max().unwrap_or(&1);
        let sparkline_chars = ['▁', '▂', '▃', '▄', '▅', '▆', '▇', '█'];

        let spans: Vec<Span> = depths.iter().map(|&d| {
            let level = if max_depth > 0 {
                ((d as f64 / max_depth as f64) * 7.0) as usize
            } else {
                0
            };
            let ch = sparkline_chars[level.min(7)];
            let color = if d < 5 {
                Color::Red
            } else if d < 20 {
                Color::Yellow
            } else {
                Color::Green
            };
            Span::styled(ch.to_string(), Style::default().fg(color))
        }).collect();

        Line::from(spans)
    }

    /// Render read name labels (left margin).
    fn render_read_label(&self, read: &AlignedRead, haplotype: Option<HaplotypeLabel>) -> Line<'static> {
        let short_name = if read.name.len() > 20 {
            format!("{}...", &read.name[..17])
        } else {
            read.name.clone()
        };

        let hp_marker = match haplotype {
            Some(HaplotypeLabel::Hap1) => "H1",
            Some(HaplotypeLabel::Hap2) => "H2",
            _ => "  ",
        };

        let color = match haplotype {
            Some(HaplotypeLabel::Hap1) => Color::Cyan,
            Some(HaplotypeLabel::Hap2) => Color::Magenta,
            _ => Color::Gray,
        };

        Line::from(vec![
            Span::styled(hp_marker.to_string(), Style::default().fg(color).add_modifier(Modifier::BOLD)),
            Span::raw(" "),
            Span::styled(short_name, Style::default().fg(Color::DarkGray)),
        ])
    }
}

impl Widget for EnhancedPileupView<'_> {
    fn render(self, area: Rect, buf: &mut Buffer) {
        let block = Block::default()
            .title(format!(" {} [{}] ", self.region_label,
                match self.view_mode {
                    ViewMode::Pileup => "Pileup",
                    ViewMode::SVDetail => "SV Detail",
                    ViewMode::HaplotypeCompare => "Haplotype Compare",
                    ViewMode::ReadDetail => "Read Detail",
                }
            ))
            .borders(Borders::ALL);
        let inner = block.inner(area);
        block.render(area, buf);

        // Reserve space for read names on the left
        let label_width = 25u16;
        let seq_width = inner.width.saturating_sub(label_width) as usize;
        let seq_x = inner.x + label_width;

        let mut y = inner.y;

        // Row 1: Coordinate ruler
        if y < inner.y + inner.height {
            let ruler = render_ruler(self.ref_start, self.offset, seq_width);
            buf.set_string(inner.x, y, "Position", Style::default().fg(Color::DarkGray));
            buf.set_line(seq_x, y, &ruler, seq_width as u16);
            y += 1;
        }

        // Row 2: Depth track
        if y < inner.y + inner.height {
            let depth = self.render_depth_track(seq_width);
            buf.set_string(inner.x, y, "Depth", Style::default().fg(Color::DarkGray));
            buf.set_line(seq_x, y, &depth, seq_width as u16);
            y += 1;
        }

        // Row 3: SV track
        if !self.sv_events.is_empty() && y < inner.y + inner.height {
            let sv_track = self.render_sv_track(seq_width);
            let sv_label = format!("SVs ({})", self.sv_events.len());
            buf.set_string(inner.x, y, &sv_label, Style::default().fg(Color::Yellow));
            buf.set_line(seq_x, y, &sv_track, seq_width as u16);
            y += 1;
        }

        // Row 4: Reference track
        if y < inner.y + inner.height {
            let ref_line = render_reference_line(self.reference, self.offset, seq_width);
            buf.set_string(inner.x, y, "Reference", Style::default().fg(Color::Green).add_modifier(Modifier::BOLD));
            buf.set_line(seq_x, y, &ref_line, seq_width as u16);
            y += 1;
        }

        // Row 5: Assembly track (if present)
        if let Some(assembly) = self.assembly {
            if y < inner.y + inner.height {
                let asm_line = render_assembly_line(assembly, self.offset, seq_width, self.reference);
                buf.set_string(inner.x, y, "Assembly", Style::default().fg(Color::Cyan).add_modifier(Modifier::BOLD));
                buf.set_line(seq_x, y, &asm_line, seq_width as u16);
                y += 1;
            }
        }

        // Row 6: Metrics bar (if present)
        if let Some(metrics) = self.metrics {
            if y < inner.y + inner.height {
                let met_line = render_metrics_line(metrics, self.offset, seq_width);
                buf.set_string(inner.x, y, "Agreement", Style::default().fg(Color::DarkGray));
                buf.set_line(seq_x, y, &met_line, seq_width as u16);
                y += 1;
            }
        }

        // Separator
        if y < inner.y + inner.height {
            let sep = "─".repeat(inner.width as usize);
            buf.set_string(inner.x, y, &sep, Style::default().fg(Color::DarkGray));
            y += 1;
        }

        // Reads pileup
        let reads_to_show = self.reads
            .iter()
            .enumerate()
            .skip(self.scroll_y)
            .take((inner.y + inner.height - y) as usize);

        for (idx, read) in reads_to_show {
            if y >= inner.y + inner.height {
                break;
            }

            let hap = self.haplotype_assignments
                .iter()
                .find(|(name, _)| name == &read.name)
                .map(|(_, h)| *h);

            let is_selected = self.selected_read == Some(idx);

            // Read label
            let label = self.render_read_label(read, hap);
            if is_selected {
                buf.set_string(inner.x, y, ">", Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD));
            }
            buf.set_line(inner.x + 1, y, &label, label_width - 1);

            // Read sequence with CIGAR-aware rendering
            let spans = self.render_read_with_cigar(read, seq_width, hap);
            let read_line = Line::from(spans);
            buf.set_line(seq_x, y, &read_line, seq_width as u16);

            y += 1;
        }
    }
}

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

/// Render the reference track as a line of colored bases.
pub fn render_reference_line(reference: &[u8], offset: usize, width: usize) -> Line<'static> {
    let end = std::cmp::min(offset + width, reference.len());
    if offset >= reference.len() {
        return Line::from(vec![Span::raw(" ".repeat(width))]);
    }
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

/// Render assembly consensus line with confidence coloring.
pub fn render_assembly_line(
    assembly: &AssemblyResult,
    offset: usize,
    width: usize,
    reference: &[u8],
) -> Line<'static> {
    let end = std::cmp::min(offset + width, assembly.sequence.len());
    if offset >= assembly.sequence.len() {
        return Line::from(vec![Span::raw(" ".repeat(width))]);
    }
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
    if offset >= metrics.len() {
        return Line::from(vec![Span::raw(" ".repeat(width))]);
    }
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

/// Render a coordinate ruler line.
pub fn render_ruler(ref_start: u64, offset: usize, width: usize) -> Line<'static> {
    let start_pos = ref_start + offset as u64;
    let mut ruler = String::with_capacity(width);

    for i in 0..width {
        let pos = start_pos + i as u64;
        if pos % 10 == 0 {
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

