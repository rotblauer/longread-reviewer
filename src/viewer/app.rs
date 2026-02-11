use std::io;

use anyhow::Result;
use crossterm::event::{self, Event, KeyCode, KeyEventKind};
use crossterm::terminal::{disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen};
use crossterm::ExecutableCommand;
use ratatui::backend::CrosstermBackend;
use ratatui::layout::{Constraint, Direction, Layout};
use ratatui::Terminal;

use crate::alignment::AlignedRead;
use crate::assembly::engine::AssemblyEngine;
use crate::assembly::method::{AssemblyResult, ConsensusAssembly, HaplotypeAssemblyResult, WindowConsensusAssembly};
use crate::haplotype::{HaplotypeAssigner, HaplotypeLabel};
use crate::metrics::{FitnessScore, MetricsCalculator};
use crate::region::Region;
use crate::viewer::render::{render_status_bar, PileupView};

/// Application state for the TUI viewer.
pub struct App {
    /// Current genomic region being viewed.
    pub region: Region,
    /// Reference sequence for the region.
    pub reference: Vec<u8>,
    /// Aligned reads in the region.
    pub reads: Vec<AlignedRead>,
    /// Current assembly result.
    pub assembly: Option<AssemblyResult>,
    /// Current fitness score.
    pub fitness: Option<FitnessScore>,
    /// Haplotype assignments for reads.
    pub haplotype_assignments: Vec<(String, HaplotypeLabel)>,
    /// Per-haplotype assembly result with divergence and SV likelihood.
    pub haplotype_assembly: Option<HaplotypeAssemblyResult>,
    /// Horizontal scroll offset.
    pub scroll_x: usize,
    /// Vertical scroll offset for reads.
    pub scroll_y: usize,
    /// Assembly engine with registered methods.
    pub engine: AssemblyEngine,
    /// Index of current assembly method.
    pub current_method: usize,
    /// Whether the app should quit.
    pub should_quit: bool,
}

impl App {
    /// Create a new App with the given data.
    pub fn new(region: Region, reference: Vec<u8>, reads: Vec<AlignedRead>) -> Self {
        let mut engine = AssemblyEngine::new();
        engine.add_method(Box::new(ConsensusAssembly));
        engine.add_method(Box::new(WindowConsensusAssembly::new(100, 20)));
        engine.add_method(Box::new(WindowConsensusAssembly::new(50, 10)));

        Self {
            region,
            reference,
            reads,
            assembly: None,
            fitness: None,
            haplotype_assignments: Vec::new(),
            haplotype_assembly: None,
            scroll_x: 0,
            scroll_y: 0,
            engine,
            current_method: 0,
            should_quit: false,
        }
    }

    /// Run the assembly with the current method.
    pub fn run_assembly(&mut self) -> Result<()> {
        let methods = self.engine.method_names();
        if methods.is_empty() {
            return Ok(());
        }

        let method_name = methods[self.current_method % methods.len()];
        if let Some(assembly) = self
            .engine
            .run_method(method_name, &self.reads, &self.reference, self.region.start)?
        {
            let calc = MetricsCalculator::new();
            let fitness = calc.compute_fitness(&assembly, &self.reads, &self.reference, self.region.start);
            self.fitness = Some(fitness);
            self.assembly = Some(assembly);
        }
        Ok(())
    }

    /// Assign reads to haplotypes.
    pub fn assign_haplotypes(&mut self) {
        let assigner = HaplotypeAssigner::new();
        let assignments = assigner.assign(&self.reads, &self.reference);
        self.haplotype_assignments = assignments
            .into_iter()
            .map(|a| (a.read_name, a.haplotype))
            .collect();

        // Sort reads by haplotype for cleaner visualization
        let hap_map: std::collections::HashMap<String, HaplotypeLabel> = self
            .haplotype_assignments
            .iter()
            .cloned()
            .collect();

        self.reads.sort_by(|a, b| {
            let ha = hap_map.get(&a.name).unwrap_or(&HaplotypeLabel::Unassigned);
            let hb = hap_map.get(&b.name).unwrap_or(&HaplotypeLabel::Unassigned);
            ha.to_string().cmp(&hb.to_string()).then(a.start.cmp(&b.start))
        });
    }

    /// Run per-haplotype assembly, splitting reads by haplotype and assembling
    /// each group independently. Produces divergence and SV likelihood tracks.
    pub fn run_haplotype_assembly(&mut self) -> Result<()> {
        let result = self.engine.assemble_by_haplotype(
            &self.reads,
            &self.reference,
            self.region.start,
        )?;
        self.haplotype_assembly = Some(result);
        Ok(())
    }

    /// Cycle to the next assembly method and re-run.
    pub fn next_method(&mut self) -> Result<()> {
        let num_methods = self.engine.method_names().len();
        if num_methods > 0 {
            self.current_method = (self.current_method + 1) % num_methods;
            self.run_assembly()?;
        }
        Ok(())
    }

    /// Evaluate all assembly methods and report results.
    pub fn evaluate_all(&self) -> Result<Vec<(String, f64)>> {
        let results = self
            .engine
            .evaluate_all(&self.reads, &self.reference, &self.region)?;
        Ok(results
            .into_iter()
            .map(|r| (r.assembly.method_name, r.fitness.overall))
            .collect())
    }

    /// Handle a key event.
    pub fn handle_key(&mut self, code: KeyCode) -> Result<()> {
        match code {
            KeyCode::Char('q') | KeyCode::Esc => self.should_quit = true,
            KeyCode::Left => {
                self.scroll_x = self.scroll_x.saturating_sub(10);
            }
            KeyCode::Right => {
                let max = self.reference.len().saturating_sub(1);
                self.scroll_x = (self.scroll_x + 10).min(max);
            }
            KeyCode::Up => {
                self.scroll_y = self.scroll_y.saturating_sub(1);
            }
            KeyCode::Down => {
                let max = self.reads.len().saturating_sub(1);
                self.scroll_y = (self.scroll_y + 1).min(max);
            }
            KeyCode::Char('a') => {
                self.run_assembly()?;
            }
            KeyCode::Char('h') => {
                self.assign_haplotypes();
            }
            KeyCode::Char('H') => {
                self.run_haplotype_assembly()?;
            }
            KeyCode::Char('n') => {
                self.next_method()?;
            }
            KeyCode::Char('e') => {
                // Evaluate all methods (results shown in status)
                let _ = self.evaluate_all();
            }
            _ => {}
        }
        Ok(())
    }

    /// Run the TUI event loop.
    pub fn run_tui(&mut self) -> Result<()> {
        enable_raw_mode()?;
        io::stdout().execute(EnterAlternateScreen)?;
        let backend = CrosstermBackend::new(io::stdout());
        let mut terminal = Terminal::new(backend)?;

        while !self.should_quit {
            terminal.draw(|frame| {
                let chunks = Layout::default()
                    .direction(Direction::Vertical)
                    .constraints([Constraint::Min(5), Constraint::Length(1)])
                    .split(frame.area());

                let method_name = self.assembly.as_ref().map(|a| a.method_name.as_str());
                let fitness_score = self.fitness.as_ref().map(|f| f.overall);

                let status = render_status_bar(
                    &self.region.to_string(),
                    self.reads.len(),
                    method_name,
                    fitness_score,
                );
                frame.render_widget(status, chunks[1]);

                let pileup = PileupView {
                    reference: &self.reference,
                    reads: &self.reads,
                    assembly: self.assembly.as_ref(),
                    metrics: self
                        .fitness
                        .as_ref()
                        .map(|f| f.base_metrics.as_slice()),
                    haplotype_assignments: &self.haplotype_assignments,
                    ref_start: self.region.start,
                    offset: self.scroll_x,
                    scroll_y: self.scroll_y,
                    region_label: &self.region.to_string(),
                };
                frame.render_widget(pileup, chunks[0]);
            })?;

            if event::poll(std::time::Duration::from_millis(100))?
                && let Event::Key(key) = event::read()?
                    && key.kind == KeyEventKind::Press {
                        self.handle_key(key.code)?;
                    }
        }

        disable_raw_mode()?;
        io::stdout().execute(LeaveAlternateScreen)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::CigarOp;

    fn make_test_app() -> App {
        let region = Region::new("chr1", 1, 8).unwrap();
        let reference = b"ACGTACGT".to_vec();
        let reads = vec![
            AlignedRead {
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
            },
            AlignedRead {
                name: "read2".to_string(),
                chrom: "chr1".to_string(),
                start: 1,
                mapq: 60,
                cigar: vec![CigarOp::Match(8)],
                sequence: b"ACGTACGT".to_vec(),
                qualities: vec![30; 8],
                is_reverse: false,
                haplotype_tag: None,
                end: 8,
            },
        ];

        App::new(region, reference, reads)
    }

    #[test]
    fn test_app_creation() {
        let app = make_test_app();
        assert_eq!(app.reads.len(), 2);
        assert_eq!(app.reference, b"ACGTACGT");
        assert!(!app.should_quit);
    }

    #[test]
    fn test_run_assembly() {
        let mut app = make_test_app();
        app.run_assembly().unwrap();
        assert!(app.assembly.is_some());
        assert!(app.fitness.is_some());
    }

    #[test]
    fn test_assign_haplotypes() {
        let mut app = make_test_app();
        app.assign_haplotypes();
        assert_eq!(app.haplotype_assignments.len(), 2);
    }

    #[test]
    fn test_next_method() {
        let mut app = make_test_app();
        assert_eq!(app.current_method, 0);
        app.next_method().unwrap();
        assert_eq!(app.current_method, 1);
        assert!(app.assembly.is_some());
    }

    #[test]
    fn test_evaluate_all() {
        let app = make_test_app();
        let results = app.evaluate_all().unwrap();
        assert_eq!(results.len(), 3); // 3 default methods
        for (name, score) in &results {
            assert!(!name.is_empty());
            assert!(*score >= 0.0);
            assert!(*score <= 1.0);
        }
    }

    #[test]
    fn test_handle_key_quit() {
        let mut app = make_test_app();
        app.handle_key(KeyCode::Char('q')).unwrap();
        assert!(app.should_quit);
    }

    #[test]
    fn test_handle_key_scroll() {
        let mut app = make_test_app();
        app.handle_key(KeyCode::Right).unwrap();
        assert!(app.scroll_x > 0);
        app.handle_key(KeyCode::Left).unwrap();
        // May be back to 0 depending on step size
    }
}
