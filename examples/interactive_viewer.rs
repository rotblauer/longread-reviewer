//! Interactive IGV-like viewer for longread alignment review.
//!
//! Usage: cargo run --example interactive_viewer
//!
//! Features:
//! - Navigate with arrow keys (Shift for faster)
//! - Jump between structural variants with [ and ]
//! - Toggle SV panel with 's'
//! - Run assembly with 'a', cycle methods with 'n'
//! - Assign haplotypes with 'h'
//! - Toggle insertions with 'i', mismatches with 'm'
//! - Press '?' for help

use std::io;
use std::path::PathBuf;
use std::collections::HashMap;

use anyhow::Result;
use crossterm::event::{self, Event, KeyCode, KeyEventKind, KeyModifiers};
use crossterm::terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen};
use crossterm::ExecutableCommand;
use ratatui::backend::CrosstermBackend;
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Color, Modifier, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, Clear, List, ListItem, ListState, Paragraph, Widget, Wrap};
use ratatui::buffer::Buffer;
use ratatui::Terminal;

use longread_reviewer::alignment::{AlignedRead, AlignmentReader, CigarOp};
use longread_reviewer::assembly::engine::AssemblyEngine;
use longread_reviewer::assembly::method::{AssemblyResult, ConsensusAssembly, WindowConsensusAssembly};
use longread_reviewer::haplotype::{HaplotypeAssigner, HaplotypeLabel};
use longread_reviewer::metrics::{FitnessScore, MetricsCalculator};
use longread_reviewer::reference::ReferenceGenome;
use longread_reviewer::region::Region;

/// Detected structural variant.
#[derive(Debug, Clone)]
struct SVEvent {
    sv_type: String,
    start: u64,
    end: u64,
    size: i64,
    supporting_reads: Vec<String>,
}

/// Application state.
struct App {
    region: Region,
    reference: Vec<u8>,
    reads: Vec<AlignedRead>,
    assembly: Option<AssemblyResult>,
    fitness: Option<FitnessScore>,
    haplotype_assignments: Vec<(String, HaplotypeLabel)>,
    sv_events: Vec<SVEvent>,
    scroll_x: usize,
    scroll_y: usize,
    engine: AssemblyEngine,
    current_method: usize,
    selected_sv: Option<usize>,
    show_insertions: bool,
    highlight_mismatches: bool,
    show_help: bool,
    show_sv_panel: bool,
    sv_list_state: ListState,
    should_quit: bool,
    status_message: Option<String>,
}

impl App {
    fn new(region: Region, reference: Vec<u8>, reads: Vec<AlignedRead>) -> Self {
        let mut engine = AssemblyEngine::new();
        engine.add_method(Box::new(ConsensusAssembly));
        engine.add_method(Box::new(WindowConsensusAssembly::new(100, 20)));
        engine.add_method(Box::new(WindowConsensusAssembly::new(50, 10)));

        let mut app = Self {
            region, reference, reads, assembly: None, fitness: None,
            haplotype_assignments: Vec::new(), sv_events: Vec::new(),
            scroll_x: 0, scroll_y: 0, engine, current_method: 0,
            selected_sv: None, show_insertions: true, highlight_mismatches: true,
            show_help: false, show_sv_panel: true, sv_list_state: ListState::default(),
            should_quit: false, status_message: None,
        };
        app.detect_svs();
        let _ = app.run_assembly();
        app.assign_haplotypes();

        // Jump to first SV if any
        if !app.sv_events.is_empty() {
            app.sv_list_state.select(Some(0));
            app.goto_selected_sv();
        }

        app
    }

    fn detect_svs(&mut self) {
        let mut all_svs: Vec<(String, String, u64, u64, i64)> = Vec::new();
        for read in &self.reads {
            let mut ref_pos = read.start;
            for op in &read.cigar {
                match op {
                    CigarOp::Match(n) => ref_pos += *n as u64,
                    CigarOp::Insertion(n) if *n >= 30 => {
                        all_svs.push((read.name.clone(), "INS".into(), ref_pos, ref_pos, *n as i64));
                    }
                    CigarOp::Deletion(n) if *n >= 30 => {
                        all_svs.push((read.name.clone(), "DEL".into(), ref_pos, ref_pos + *n as u64, -(*n as i64)));
                        ref_pos += *n as u64;
                    }
                    CigarOp::Deletion(n) => ref_pos += *n as u64,
                    _ => {}
                }
            }
        }

        // Group similar SVs
        let mut sv_events: Vec<SVEvent> = Vec::new();
        for (rn, st, s, e, sz) in &all_svs {
            let mut found = false;
            for ev in &mut sv_events {
                if ev.sv_type == *st && ((ev.start as i64 - *s as i64).abs() < 50) && (ev.size - sz).abs() < 50 {
                    if !ev.supporting_reads.contains(rn) { ev.supporting_reads.push(rn.clone()); }
                    found = true; break;
                }
            }
            if !found { sv_events.push(SVEvent { sv_type: st.clone(), start: *s, end: *e, size: *sz, supporting_reads: vec![rn.clone()] }); }
        }
        sv_events.sort_by_key(|e| e.start);
        self.sv_events = sv_events;
    }

    fn run_assembly(&mut self) -> Result<()> {
        let methods = self.engine.method_names();
        if methods.is_empty() { return Ok(()); }
        let mn = methods[self.current_method % methods.len()];
        if let Some(asm) = self.engine.run_method(mn, &self.reads, &self.reference, self.region.start)? {
            let calc = MetricsCalculator::new();
            let fit = calc.compute_fitness(&asm, &self.reads, &self.reference, self.region.start);
            self.status_message = Some(format!("Assembly: {} | Fitness: {:.3} | Agreement: {:.1}%", mn, fit.overall, fit.mean_agreement * 100.0));
            self.fitness = Some(fit); self.assembly = Some(asm);
        }
        Ok(())
    }

    fn assign_haplotypes(&mut self) {
        let assigner = HaplotypeAssigner::new();
        let asgn = assigner.assign(&self.reads, &self.reference);
        self.haplotype_assignments = asgn.into_iter().map(|a| (a.read_name, a.haplotype)).collect();
        let hm: HashMap<_, _> = self.haplotype_assignments.iter().cloned().collect();
        self.reads.sort_by(|a, b| {
            let ha = hm.get(&a.name).unwrap_or(&HaplotypeLabel::Unassigned);
            let hb = hm.get(&b.name).unwrap_or(&HaplotypeLabel::Unassigned);
            ha.to_string().cmp(&hb.to_string()).then(a.start.cmp(&b.start))
        });
    }

    fn goto_position(&mut self, pos: u64) {
        if pos >= self.region.start && pos <= self.region.end {
            // Center the position in view
            self.scroll_x = (pos.saturating_sub(self.region.start).saturating_sub(40)) as usize;
        }
    }

    fn goto_selected_sv(&mut self) {
        if let Some(idx) = self.sv_list_state.selected() {
            if idx < self.sv_events.len() {
                let sv_start = self.sv_events[idx].start;
                let sv_type = self.sv_events[idx].sv_type.clone();
                let sv_size = self.sv_events[idx].size;
                let sv_end = self.sv_events[idx].end;
                let sv_reads_len = self.sv_events[idx].supporting_reads.len();
                self.goto_position(sv_start);
                self.selected_sv = Some(idx);
                self.status_message = Some(format!("▶ {} {}bp @ chr17:{}-{} ({} reads)",
                    sv_type, sv_size.abs(), sv_start, sv_end, sv_reads_len));
            }
        }
    }

    fn next_sv(&mut self) {
        if self.sv_events.is_empty() { return; }
        let i = self.sv_list_state.selected().map(|i| (i + 1) % self.sv_events.len()).unwrap_or(0);
        self.sv_list_state.select(Some(i)); self.goto_selected_sv();
    }

    fn prev_sv(&mut self) {
        if self.sv_events.is_empty() { return; }
        let i = self.sv_list_state.selected().map(|i| if i == 0 { self.sv_events.len() - 1 } else { i - 1 }).unwrap_or(0);
        self.sv_list_state.select(Some(i)); self.goto_selected_sv();
    }

    fn handle_key(&mut self, code: KeyCode, mods: KeyModifiers) -> Result<()> {
        if self.show_help { self.show_help = false; return Ok(()); }
        match code {
            KeyCode::Char('q') | KeyCode::Esc => self.should_quit = true,
            KeyCode::Left => { let s = if mods.contains(KeyModifiers::SHIFT) { 100 } else { 10 }; self.scroll_x = self.scroll_x.saturating_sub(s); self.status_message = None; }
            KeyCode::Right => { let s = if mods.contains(KeyModifiers::SHIFT) { 100 } else { 10 }; self.scroll_x = (self.scroll_x + s).min(self.reference.len().saturating_sub(1)); self.status_message = None; }
            KeyCode::Up => { self.scroll_y = self.scroll_y.saturating_sub(1); self.status_message = None; }
            KeyCode::Down => { self.scroll_y = (self.scroll_y + 1).min(self.reads.len().saturating_sub(1)); self.status_message = None; }
            KeyCode::Home => { self.scroll_x = 0; self.status_message = Some("Start of region".into()); }
            KeyCode::End => { self.scroll_x = self.reference.len().saturating_sub(80); self.status_message = Some("End of region".into()); }
            KeyCode::PageUp => self.scroll_y = self.scroll_y.saturating_sub(20),
            KeyCode::PageDown => self.scroll_y = (self.scroll_y + 20).min(self.reads.len().saturating_sub(1)),
            KeyCode::Char('a') => { self.run_assembly()?; }
            KeyCode::Char('n') => { self.current_method = (self.current_method + 1) % self.engine.method_names().len().max(1); self.run_assembly()?; }
            KeyCode::Char('h') => { self.assign_haplotypes(); self.status_message = Some("Haplotypes assigned".into()); }
            KeyCode::Char('s') => self.show_sv_panel = !self.show_sv_panel,
            KeyCode::Char(']') | KeyCode::Tab => self.next_sv(),
            KeyCode::Char('[') | KeyCode::BackTab => self.prev_sv(),
            KeyCode::Enter => self.goto_selected_sv(),
            KeyCode::Char('i') => { self.show_insertions = !self.show_insertions; self.status_message = Some(format!("Insertions: {}", if self.show_insertions { "ON" } else { "OFF" })); }
            KeyCode::Char('m') => { self.highlight_mismatches = !self.highlight_mismatches; self.status_message = Some(format!("Mismatches: {}", if self.highlight_mismatches { "ON" } else { "OFF" })); }
            KeyCode::Char('?') | KeyCode::F(1) => self.show_help = true,
            _ => {}
        }
        Ok(())
    }

    fn run_tui(&mut self) -> Result<()> {
        enable_raw_mode()?;
        io::stdout().execute(EnterAlternateScreen)?;
        let mut term = Terminal::new(CrosstermBackend::new(io::stdout()))?;

        while !self.should_quit {
            term.draw(|f| self.render_ui(f))?;
            if event::poll(std::time::Duration::from_millis(50))? {
                if let Event::Key(k) = event::read()? {
                    if k.kind == KeyEventKind::Press { let _ = self.handle_key(k.code, k.modifiers); }
                }
            }
        }

        disable_raw_mode()?;
        io::stdout().execute(LeaveAlternateScreen)?;
        Ok(())
    }

    fn render_ui(&self, f: &mut ratatui::Frame) {
        let sz = f.area();
        let main = Layout::default().direction(Direction::Vertical).constraints([Constraint::Min(10), Constraint::Length(2)]).split(sz);

        if self.show_sv_panel && !self.sv_events.is_empty() {
            let h = Layout::default().direction(Direction::Horizontal).constraints([Constraint::Length(40), Constraint::Min(40)]).split(main[0]);
            self.render_sv_panel(f, h[0]);
            self.render_pileup(f, h[1]);
        } else {
            self.render_pileup(f, main[0]);
        }

        self.render_status_bar(f, main[1]);
        if self.show_help { self.render_help(f, sz); }
    }

    fn render_sv_panel(&self, f: &mut ratatui::Frame, area: Rect) {
        let total_del: i64 = self.sv_events.iter().filter(|e| e.sv_type == "DEL").map(|e| e.size.abs()).sum();
        let total_ins: i64 = self.sv_events.iter().filter(|e| e.sv_type == "INS").map(|e| e.size.abs()).sum();

        let items: Vec<ListItem> = self.sv_events.iter().enumerate().map(|(i, sv)| {
            let sym = if sv.sv_type == "DEL" { "▼DEL" } else { "▲INS" };
            let col = if sv.sv_type == "DEL" { Color::Red } else { Color::Green };
            let sel = self.sv_list_state.selected() == Some(i);
            let sty = if sel { Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD | Modifier::REVERSED) } else { Style::default().fg(col) };
            ListItem::new(Line::from(vec![
                Span::styled(format!("{} {:>5}bp ", sym, sv.size.abs()), sty),
                Span::styled(format!("@{} ", sv.start), Style::default().fg(Color::DarkGray)),
                Span::styled(format!("({})", sv.supporting_reads.len()), Style::default().fg(Color::Cyan)),
            ]))
        }).collect();

        let list = List::new(items)
            .block(Block::default()
                .title(format!(" {} SVs: -{}bp/+{}bp ", self.sv_events.len(), total_del, total_ins))
                .borders(Borders::ALL)
                .border_style(Style::default().fg(Color::Yellow)));
        f.render_stateful_widget(list, area, &mut self.sv_list_state.clone());
    }

    fn render_pileup(&self, f: &mut ratatui::Frame, area: Rect) {
        let pv = PileupWidget {
            reference: &self.reference,
            reads: &self.reads,
            assembly: self.assembly.as_ref(),
            fitness: self.fitness.as_ref(),
            haplotype_assignments: &self.haplotype_assignments,
            sv_events: &self.sv_events,
            ref_start: self.region.start,
            offset: self.scroll_x,
            scroll_y: self.scroll_y,
            selected_sv: self.selected_sv,
            show_insertions: self.show_insertions,
            highlight_mismatches: self.highlight_mismatches,
            region_label: &self.region.to_string(),
        };
        f.render_widget(pv, area);
    }

    fn render_status_bar(&self, f: &mut ratatui::Frame, area: Rect) {
        let pos = self.region.start + self.scroll_x as u64;
        let mut spans = vec![
            Span::styled(format!(" chr17:{} ", pos), Style::default().fg(Color::White).bg(Color::Blue)),
            Span::styled(format!(" {} reads ", self.reads.len()), Style::default().fg(Color::White).bg(Color::DarkGray)),
        ];
        if let Some(ref a) = self.assembly {
            spans.push(Span::styled(format!(" {} ", a.method_name), Style::default().fg(Color::Cyan).bg(Color::DarkGray)));
        }
        if let Some(ref fit) = self.fitness {
            let col = if fit.overall > 0.95 { Color::Green } else if fit.overall > 0.8 { Color::Yellow } else { Color::Red };
            spans.push(Span::styled(format!(" fit:{:.3} ", fit.overall), Style::default().fg(col).bg(Color::DarkGray)));
        }
        if !self.sv_events.is_empty() {
            spans.push(Span::styled(format!(" {} SVs ", self.sv_events.len()), Style::default().fg(Color::Yellow).bg(Color::DarkGray)));
        }
        f.render_widget(Paragraph::new(Line::from(spans)), Rect { height: 1, ..area });

        let msg = self.status_message.clone().unwrap_or_else(||
            "[?]Help [←→]Scroll [↑↓]Reads [Tab]NextSV [s]Panel [a]Asm [n]Method [h]Hap [i]Ins [m]Mis".into()
        );
        f.render_widget(Paragraph::new(Span::styled(format!(" {}", msg), Style::default().fg(Color::DarkGray))), Rect { y: area.y + 1, height: 1, ..area });
    }

    fn render_help(&self, f: &mut ratatui::Frame, area: Rect) {
        let txt = "\
NAVIGATION
  ←/→         Scroll 10bp (Shift: 100bp)
  ↑/↓         Scroll reads
  Home/End    Jump to start/end
  PageUp/Down Scroll 20 reads

STRUCTURAL VARIANTS
  Tab / ]     Next SV
  Shift-Tab   Previous SV
  [           Previous SV
  Enter       Jump to selected SV
  s           Toggle SV panel

ASSEMBLY & ANALYSIS
  a           Run assembly
  n           Next assembly method
  h           Assign haplotypes

DISPLAY
  i           Toggle insertions
  m           Toggle mismatch highlighting

OTHER
  ?/F1        This help
  q/Esc       Quit";

        let (w, h) = (48u16, 26u16);
        let pop = Rect::new((area.width.saturating_sub(w)) / 2, (area.height.saturating_sub(h)) / 2, w, h);
        f.render_widget(Clear, pop);
        f.render_widget(
            Paragraph::new(txt)
                .block(Block::default()
                    .title(" Help (any key to close) ")
                    .borders(Borders::ALL)
                    .border_style(Style::default().fg(Color::Yellow)))
                .wrap(Wrap { trim: false }),
            pop
        );
    }
}

/// Pileup rendering widget.
struct PileupWidget<'a> {
    reference: &'a [u8],
    reads: &'a [AlignedRead],
    assembly: Option<&'a AssemblyResult>,
    fitness: Option<&'a FitnessScore>,
    haplotype_assignments: &'a [(String, HaplotypeLabel)],
    sv_events: &'a [SVEvent],
    ref_start: u64,
    offset: usize,
    scroll_y: usize,
    selected_sv: Option<usize>,
    show_insertions: bool,
    highlight_mismatches: bool,
    region_label: &'a str,
}

impl Widget for PileupWidget<'_> {
    fn render(self, area: Rect, buf: &mut Buffer) {
        let block = Block::default()
            .title(format!(" {} ", self.region_label))
            .borders(Borders::ALL);
        let inner = block.inner(area);
        block.render(area, buf);

        let label_w = 28u16;
        let seq_w = inner.width.saturating_sub(label_w) as usize;
        let seq_x = inner.x + label_w;
        let mut y = inner.y;

        // Ruler
        if y < inner.y + inner.height {
            buf.set_string(inner.x, y, "Position", Style::default().fg(Color::DarkGray));
            let ruler = render_ruler(self.ref_start, self.offset, seq_w);
            buf.set_line(seq_x, y, &ruler, seq_w as u16);
            y += 1;
        }

        // Depth sparkline
        if y < inner.y + inner.height {
            buf.set_string(inner.x, y, "Depth", Style::default().fg(Color::DarkGray));
            let depth = self.render_depth(seq_w);
            buf.set_line(seq_x, y, &depth, seq_w as u16);
            y += 1;
        }

        // SV track
        if !self.sv_events.is_empty() && y < inner.y + inner.height {
            buf.set_string(inner.x, y, "SVs", Style::default().fg(Color::Yellow));
            let sv_track = self.render_sv_track(seq_w);
            buf.set_line(seq_x, y, &sv_track, seq_w as u16);
            y += 1;
        }

        // Reference
        if y < inner.y + inner.height {
            buf.set_string(inner.x, y, "Reference", Style::default().fg(Color::Green).add_modifier(Modifier::BOLD));
            let ref_line = self.render_reference(seq_w);
            buf.set_line(seq_x, y, &ref_line, seq_w as u16);
            y += 1;
        }

        // Assembly
        if let Some(asm) = self.assembly {
            if y < inner.y + inner.height {
                buf.set_string(inner.x, y, "Assembly", Style::default().fg(Color::Cyan).add_modifier(Modifier::BOLD));
                let asm_line = self.render_assembly(asm, seq_w);
                buf.set_line(seq_x, y, &asm_line, seq_w as u16);
                y += 1;
            }
        }

        // Separator
        if y < inner.y + inner.height {
            buf.set_string(inner.x, y, &"─".repeat(inner.width as usize), Style::default().fg(Color::DarkGray));
            y += 1;
        }

        // Reads
        for read in self.reads.iter().skip(self.scroll_y).take((inner.y + inner.height - y) as usize) {
            if y >= inner.y + inner.height { break; }

            let hap = self.haplotype_assignments.iter().find(|(n, _)| n == &read.name).map(|(_, h)| *h);

            // Label
            let label = self.render_read_label(read, hap);
            buf.set_line(inner.x, y, &label, label_w);

            // Sequence
            let seq_line = self.render_read(read, hap, seq_w);
            buf.set_line(seq_x, y, &seq_line, seq_w as u16);
            y += 1;
        }
    }
}

impl PileupWidget<'_> {
    fn render_depth(&self, width: usize) -> Line<'static> {
        let view_start = self.ref_start + self.offset as u64;
        let depths: Vec<u32> = (0..width).map(|i| {
            let pos = view_start + i as u64;
            self.reads.iter().filter(|r| r.start <= pos && r.end >= pos).count() as u32
        }).collect();
        let max_d = *depths.iter().max().unwrap_or(&1);
        let chars = ['▁', '▂', '▃', '▄', '▅', '▆', '▇', '█'];
        let spans: Vec<Span> = depths.iter().map(|&d| {
            let lvl = if max_d > 0 { ((d as f64 / max_d as f64) * 7.0) as usize } else { 0 };
            let col = if d < 5 { Color::Red } else if d < 20 { Color::Yellow } else { Color::Green };
            Span::styled(chars[lvl.min(7)].to_string(), Style::default().fg(col))
        }).collect();
        Line::from(spans)
    }

    fn render_sv_track(&self, width: usize) -> Line<'static> {
        let view_start = self.ref_start + self.offset as u64;
        let view_end = view_start + width as u64;
        let mut track = vec![(' ', Color::DarkGray); width];

        for (idx, sv) in self.sv_events.iter().enumerate() {
            if sv.end >= view_start && sv.start <= view_end {
                let s = if sv.start > view_start { (sv.start - view_start) as usize } else { 0 };
                let e = ((sv.end - view_start) as usize).min(width);
                let ch = if sv.sv_type == "DEL" { '▼' } else { '▲' };
                let col = if sv.sv_type == "DEL" { Color::Red } else { Color::Green };
                let is_sel = self.selected_sv == Some(idx);
                for i in s..e {
                    if i < track.len() {
                        track[i] = (ch, if is_sel { Color::Yellow } else { col });
                    }
                }
            }
        }
        Line::from(track.iter().map(|(c, col)| Span::styled(c.to_string(), Style::default().fg(*col))).collect::<Vec<_>>())
    }

    fn render_reference(&self, width: usize) -> Line<'static> {
        let end = (self.offset + width).min(self.reference.len());
        if self.offset >= self.reference.len() { return Line::from(" ".repeat(width)); }
        let spans: Vec<Span> = self.reference[self.offset..end].iter().map(|&b| {
            Span::styled((b as char).to_string(), Style::default().fg(base_color(b)).add_modifier(Modifier::BOLD))
        }).collect();
        Line::from(spans)
    }

    fn render_assembly(&self, asm: &AssemblyResult, width: usize) -> Line<'static> {
        let end = (self.offset + width).min(asm.sequence.len());
        if self.offset >= asm.sequence.len() { return Line::from(" ".repeat(width)); }
        let spans: Vec<Span> = (self.offset..end).map(|i| {
            let b = asm.sequence[i];
            let ref_b = if i < self.reference.len() { self.reference[i] } else { b'N' };
            let is_var = !b.eq_ignore_ascii_case(&ref_b);
            let (fg, bg) = if is_var { (Color::White, Color::Magenta) } else { (base_color(b), Color::Reset) };
            Span::styled((b as char).to_string(), Style::default().fg(fg).bg(bg))
        }).collect();
        Line::from(spans)
    }

    fn render_read_label(&self, read: &AlignedRead, hap: Option<HaplotypeLabel>) -> Line<'static> {
        let hp = match hap {
            Some(HaplotypeLabel::Hap1) => ("H1", Color::Cyan),
            Some(HaplotypeLabel::Hap2) => ("H2", Color::Magenta),
            _ => ("  ", Color::DarkGray),
        };
        let name = if read.name.len() > 20 { format!("{}...", &read.name[..17]) } else { read.name.clone() };
        Line::from(vec![
            Span::styled(hp.0, Style::default().fg(hp.1).add_modifier(Modifier::BOLD)),
            Span::raw(" "),
            Span::styled(name, Style::default().fg(Color::DarkGray)),
        ])
    }

    fn render_read(&self, read: &AlignedRead, hap: Option<HaplotypeLabel>, width: usize) -> Line<'static> {
        let view_start = self.ref_start + self.offset as u64;
        let view_end = view_start + width as u64;

        // Build position map from CIGAR
        let mut pos_bases: HashMap<u64, u8> = HashMap::new();
        let mut deletions: Vec<(u64, u64)> = Vec::new();
        let mut ref_pos = read.start;
        let mut read_pos = 0usize;

        for op in &read.cigar {
            match op {
                CigarOp::Match(n) => {
                    for _ in 0..*n {
                        if read_pos < read.sequence.len() {
                            pos_bases.insert(ref_pos, read.sequence[read_pos]);
                        }
                        ref_pos += 1; read_pos += 1;
                    }
                }
                CigarOp::Insertion(n) => { read_pos += *n as usize; }
                CigarOp::Deletion(n) => {
                    deletions.push((ref_pos, ref_pos + *n as u64));
                    ref_pos += *n as u64;
                }
                CigarOp::SoftClip(n) => { read_pos += *n as usize; }
                CigarOp::HardClip(_) => {}
            }
        }

        let hap_col = match hap {
            Some(HaplotypeLabel::Hap1) => Color::Cyan,
            Some(HaplotypeLabel::Hap2) => Color::Magenta,
            _ => Color::Gray,
        };

        let spans: Vec<Span> = (0..width).map(|i| {
            let pos = view_start + i as u64;
            if pos < read.start || pos > read.end { return Span::raw(" "); }

            let in_del = deletions.iter().any(|(s, e)| pos >= *s && pos < *e);
            if in_del {
                return Span::styled("─", Style::default().fg(Color::Red));
            }

            if let Some(&base) = pos_bases.get(&pos) {
                let ref_idx = (pos - self.ref_start) as usize;
                let ref_b = if ref_idx < self.reference.len() { self.reference[ref_idx] } else { b'N' };
                let is_mm = !base.eq_ignore_ascii_case(&ref_b);

                let sty = if is_mm && self.highlight_mismatches {
                    Style::default().fg(Color::White).bg(base_color(base)).add_modifier(Modifier::BOLD)
                } else {
                    Style::default().fg(hap_col)
                };
                Span::styled((base as char).to_string(), sty)
            } else {
                Span::raw(" ")
            }
        }).collect();

        Line::from(spans)
    }
}

fn base_color(b: u8) -> Color {
    match b.to_ascii_uppercase() {
        b'A' => Color::Green, b'C' => Color::Blue, b'G' => Color::Yellow, b'T' => Color::Red,
        _ => Color::DarkGray,
    }
}

fn render_ruler(ref_start: u64, offset: usize, width: usize) -> Line<'static> {
    let start = ref_start + offset as u64;
    let mut s = String::with_capacity(width);
    for i in 0..width {
        let pos = start + i as u64;
        if pos % 10 == 0 { s.push('|'); }
        else {
            let tick = pos - (pos % 10);
            let lbl = format!("{}", tick);
            let off = (pos - tick) as usize;
            if off > 0 && off < lbl.len() { s.push(lbl.as_bytes()[off] as char); }
            else { s.push('·'); }
        }
    }
    Line::from(Span::styled(s, Style::default().fg(Color::DarkGray)))
}

fn main() -> Result<()> {
    let manifest = env!("CARGO_MANIFEST_DIR");
    let bam_path = PathBuf::from(manifest).join("tests/data/NA19240.chr17_fragment.bam");
    let ref_path = PathBuf::from(manifest).join("tests/data/chr17_fragment.fa");

    let region = Region::new("chr17", 10958130, 11017414)?;

    eprintln!("Loading data...");
    let reference = ReferenceGenome::from_file(&ref_path)?;
    let ref_seq = reference.fetch(&region)?;
    let reads = AlignmentReader::read_bam(&bam_path, &region)?;

    eprintln!("Loaded {} reads for {}", reads.len(), region);
    eprintln!("Starting interactive viewer...");

    let mut app = App::new(region, ref_seq, reads);
    app.run_tui()?;

    Ok(())
}

