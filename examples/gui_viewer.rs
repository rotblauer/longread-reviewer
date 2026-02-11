//! GUI-based longread alignment viewer with zoomable multi-scale visualization.
//!
//! Usage: cargo run --example gui_viewer --release
//!
//! Features:
//! - Multi-scale zoomable view (from single-base to full-region overview)
//! - Structural variant highlighting with clear visual indicators
//! - Side-by-side comparison of multiple assembly methods
//! - Haplotype visualization with color coding
//! - Interactive navigation with pan/zoom
//! - Evidence quality metrics display

use std::collections::HashMap;
use std::path::PathBuf;

use anyhow::Result;
use eframe::egui::{self, Color32, FontId, Pos2, Rect, RichText, Sense, Stroke, Vec2};
use egui_plot::{Bar, BarChart, Legend, Plot};

use longread_reviewer::alignment::{AlignedRead, AlignmentReader, CigarOp};
use longread_reviewer::assembly::engine::AssemblyEngine;
use longread_reviewer::assembly::method::{AssemblyResult, ConsensusAssembly, HaplotypeAssemblyResult, WindowConsensusAssembly};
use longread_reviewer::haplotype::{HaplotypeAssigner, HaplotypeLabel};
use longread_reviewer::metrics::{FitnessScore, MetricsCalculator};
use longread_reviewer::reference::ReferenceGenome;
use longread_reviewer::region::Region;

// Color palette for consistent theming
const COLOR_HAP1: Color32 = Color32::from_rgb(100, 200, 255);
const COLOR_HAP2: Color32 = Color32::from_rgb(255, 150, 200);
const COLOR_UNASSIGNED: Color32 = Color32::from_rgb(180, 180, 180);
const COLOR_DELETION: Color32 = Color32::from_rgb(220, 60, 60);
const COLOR_INSERTION: Color32 = Color32::from_rgb(60, 180, 60);
const COLOR_MISMATCH: Color32 = Color32::from_rgb(255, 200, 60);
#[allow(dead_code)]
const COLOR_REFERENCE: Color32 = Color32::from_rgb(80, 80, 80);
#[allow(dead_code)]
const COLOR_ASSEMBLY: Color32 = Color32::from_rgb(100, 150, 255);

/// Structural variant event
#[derive(Debug, Clone)]
struct SVEvent {
    sv_type: SVType,
    start: u64,
    end: u64,
    size: i64,
    supporting_reads: Vec<String>,
    confidence: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum SVType {
    Deletion,
    Insertion,
}

impl SVEvent {
    fn color(&self) -> Color32 {
        match self.sv_type {
            SVType::Deletion => COLOR_DELETION,
            SVType::Insertion => COLOR_INSERTION,
        }
    }

    fn label(&self) -> &'static str {
        match self.sv_type {
            SVType::Deletion => "DEL",
            SVType::Insertion => "INS",
        }
    }
}

/// Zoom level presets
#[derive(Debug, Clone, Copy, PartialEq)]
enum ZoomLevel {
    BasePair,      // Single base view (~100bp visible)
    Fine,          // Fine detail (~500bp)
    Medium,        // Medium zoom (~2kb)
    Coarse,        // Overview (~10kb)
    Region,        // Full region view
}

impl ZoomLevel {
    fn bases_per_pixel(&self, region_len: u64) -> f64 {
        match self {
            ZoomLevel::BasePair => 0.1,
            ZoomLevel::Fine => 0.5,
            ZoomLevel::Medium => 2.0,
            ZoomLevel::Coarse => 10.0,
            ZoomLevel::Region => region_len as f64 / 1000.0,
        }
    }

    fn label(&self) -> &'static str {
        match self {
            ZoomLevel::BasePair => "Base",
            ZoomLevel::Fine => "Fine",
            ZoomLevel::Medium => "Medium",
            ZoomLevel::Coarse => "Overview",
            ZoomLevel::Region => "Full Region",
        }
    }

    fn all() -> &'static [ZoomLevel] {
        &[
            ZoomLevel::BasePair,
            ZoomLevel::Fine,
            ZoomLevel::Medium,
            ZoomLevel::Coarse,
            ZoomLevel::Region,
        ]
    }
}

/// Assembly method info
struct AssemblyInfo {
    name: String,
    result: AssemblyResult,
    fitness: FitnessScore,
}

/// Main application state
struct SVReviewerApp {
    // Data
    region: Region,
    reference: Vec<u8>,
    reads: Vec<AlignedRead>,
    sv_events: Vec<SVEvent>,
    haplotype_assignments: HashMap<String, HaplotypeLabel>,

    // Assembly state
    engine: AssemblyEngine,
    assemblies: Vec<AssemblyInfo>,
    selected_assembly: usize,
    haplotype_assembly: Option<HaplotypeAssemblyResult>,

    // View state
    zoom: ZoomLevel,
    custom_zoom: f64, // bases per pixel
    view_center: f64, // center position in genomic coordinates
    scroll_y: f32,

    // Selection state
    selected_sv: Option<usize>,
    hovered_sv: Option<usize>,
    show_sv_panel: bool,
    show_metrics_panel: bool,
    show_assembly_comparison: bool,

    // Display options
    show_insertions: bool,
    show_deletions: bool,
    highlight_mismatches: bool,
    color_by_haplotype: bool,
    show_depth_track: bool,
    show_confidence_track: bool,
    show_haplotype_assemblies: bool,

    // Computed values (cached)
    depth_by_position: Vec<u32>,
    max_depth: u32,
}

impl SVReviewerApp {
    fn new(region: Region, reference: Vec<u8>, reads: Vec<AlignedRead>) -> Self {
        let mut engine = AssemblyEngine::new();
        engine.add_method(Box::new(ConsensusAssembly));
        engine.add_method(Box::new(WindowConsensusAssembly::new(100, 20)));
        engine.add_method(Box::new(WindowConsensusAssembly::new(50, 10)));
        engine.add_method(Box::new(WindowConsensusAssembly::new(200, 40)));

        let view_center = (region.start + region.end) as f64 / 2.0;

        let mut app = Self {
            region: region.clone(),
            reference: reference.clone(),
            reads,
            sv_events: Vec::new(),
            haplotype_assignments: HashMap::new(),
            engine,
            assemblies: Vec::new(),
            selected_assembly: 0,
            haplotype_assembly: None,
            zoom: ZoomLevel::Medium,
            custom_zoom: 2.0,
            view_center,
            scroll_y: 0.0,
            selected_sv: None,
            hovered_sv: None,
            show_sv_panel: true,
            show_metrics_panel: true,
            show_assembly_comparison: false,
            show_insertions: true,
            show_deletions: true,
            highlight_mismatches: true,
            color_by_haplotype: true,
            show_depth_track: true,
            show_confidence_track: true,
            show_haplotype_assemblies: true,
            depth_by_position: Vec::new(),
            max_depth: 1,
        };

        app.detect_svs();
        app.assign_haplotypes();
        app.run_all_assemblies(&reference, region.start);
        app.run_haplotype_assembly(&reference, region.start);
        app.compute_depth_track();

        // Jump to first SV if any
        if !app.sv_events.is_empty() {
            app.selected_sv = Some(0);
            app.center_on_sv(0);
        }

        app
    }

    fn detect_svs(&mut self) {
        let mut all_svs: Vec<(String, SVType, u64, u64, i64)> = Vec::new();

        for read in &self.reads {
            let mut ref_pos = read.start;
            for op in &read.cigar {
                match op {
                    CigarOp::Match(n) => ref_pos += *n as u64,
                    CigarOp::Insertion(n) if *n >= 30 => {
                        all_svs.push((
                            read.name.clone(),
                            SVType::Insertion,
                            ref_pos,
                            ref_pos,
                            *n as i64,
                        ));
                    }
                    CigarOp::Deletion(n) if *n >= 30 => {
                        all_svs.push((
                            read.name.clone(),
                            SVType::Deletion,
                            ref_pos,
                            ref_pos + *n as u64,
                            -(*n as i64),
                        ));
                        ref_pos += *n as u64;
                    }
                    CigarOp::Deletion(n) => ref_pos += *n as u64,
                    _ => {}
                }
            }
        }

        // Group similar SVs with position tolerance
        let mut sv_events: Vec<SVEvent> = Vec::new();
        for (read_name, sv_type, start, end, size) in all_svs {
            let mut found = false;
            for ev in &mut sv_events {
                if ev.sv_type == sv_type
                    && (ev.start as i64 - start as i64).abs() < 50
                    && (ev.size - size).abs() < 50
                {
                    if !ev.supporting_reads.contains(&read_name) {
                        ev.supporting_reads.push(read_name.clone());
                    }
                    found = true;
                    break;
                }
            }
            if !found {
                sv_events.push(SVEvent {
                    sv_type,
                    start,
                    end,
                    size,
                    supporting_reads: vec![read_name],
                    confidence: 0.0,
                });
            }
        }

        // Calculate confidence based on read support
        let total_reads = self.reads.len() as f64;
        for sv in &mut sv_events {
            sv.confidence = sv.supporting_reads.len() as f64 / total_reads.max(1.0);
        }

        sv_events.sort_by_key(|e| e.start);
        self.sv_events = sv_events;
    }

    fn assign_haplotypes(&mut self) {
        let assigner = HaplotypeAssigner::new();
        let assignments = assigner.assign(&self.reads, &self.reference);
        self.haplotype_assignments = assignments
            .into_iter()
            .map(|a| (a.read_name, a.haplotype))
            .collect();

        // Sort reads by haplotype then by position
        let hm = &self.haplotype_assignments;
        self.reads.sort_by(|a, b| {
            let ha = hm.get(&a.name).unwrap_or(&HaplotypeLabel::Unassigned);
            let hb = hm.get(&b.name).unwrap_or(&HaplotypeLabel::Unassigned);
            ha.to_string()
                .cmp(&hb.to_string())
                .then(a.start.cmp(&b.start))
        });
    }

    fn run_all_assemblies(&mut self, reference: &[u8], ref_start: u64) {
        let calc = MetricsCalculator::new();
        self.assemblies.clear();

        for method_name in self.engine.method_names() {
            if let Ok(Some(result)) = self.engine.run_method(method_name, &self.reads, reference, ref_start) {
                let fitness = calc.compute_fitness(&result, &self.reads, reference, ref_start);
                self.assemblies.push(AssemblyInfo {
                    name: method_name.to_string(),
                    result,
                    fitness,
                });
            }
        }
    }

    fn run_haplotype_assembly(&mut self, reference: &[u8], ref_start: u64) {
        match self.engine.assemble_by_haplotype(&self.reads, reference, ref_start) {
            Ok(result) => {
                self.haplotype_assembly = Some(result);
            }
            Err(e) => {
                eprintln!("Warning: haplotype assembly failed: {}", e);
                self.haplotype_assembly = None;
            }
        }
    }

    fn compute_depth_track(&mut self) {
        let len = self.reference.len();
        let mut depth = vec![0u32; len];

        for read in &self.reads {
            let start_idx = read.start.saturating_sub(self.region.start) as usize;
            let end_idx = ((read.end - self.region.start) as usize).min(len);
            for i in start_idx..end_idx {
                if i < len {
                    depth[i] += 1;
                }
            }
        }

        self.max_depth = *depth.iter().max().unwrap_or(&1);
        self.depth_by_position = depth;
    }

    fn center_on_sv(&mut self, idx: usize) {
        if idx < self.sv_events.len() {
            let sv = &self.sv_events[idx];
            self.view_center = (sv.start + sv.end) as f64 / 2.0;
            // Adjust zoom based on SV size
            let sv_size = (sv.end - sv.start).max(sv.size.unsigned_abs());
            if sv_size < 100 {
                self.zoom = ZoomLevel::Fine;
            } else if sv_size < 1000 {
                self.zoom = ZoomLevel::Medium;
            } else {
                self.zoom = ZoomLevel::Coarse;
            }
        }
    }

    fn bases_per_pixel(&self) -> f64 {
        self.zoom.bases_per_pixel(self.region.len())
    }

    fn visible_range(&self, view_width: f32) -> (u64, u64) {
        let half_width = (view_width as f64 * self.bases_per_pixel()) / 2.0;
        let start = (self.view_center - half_width).max(self.region.start as f64) as u64;
        let end = (self.view_center + half_width).min(self.region.end as f64) as u64;
        (start, end)
    }

    fn pos_to_x(&self, pos: u64, view_rect: Rect) -> f32 {
        let (vis_start, vis_end) = self.visible_range(view_rect.width());
        let frac = (pos as f64 - vis_start as f64) / (vis_end - vis_start).max(1) as f64;
        view_rect.left() + (frac * view_rect.width() as f64) as f32
    }

    #[allow(dead_code)]
    fn x_to_pos(&self, x: f32, view_rect: Rect) -> u64 {
        let (vis_start, vis_end) = self.visible_range(view_rect.width());
        let frac = (x - view_rect.left()) / view_rect.width();
        (vis_start as f64 + frac as f64 * (vis_end - vis_start) as f64) as u64
    }
}

impl eframe::App for SVReviewerApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        // Top panel - toolbar
        egui::TopBottomPanel::top("toolbar").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.heading("SV Reviewer");
                ui.separator();

                // Zoom controls
                ui.label("Zoom:");
                for level in ZoomLevel::all() {
                    if ui.selectable_label(self.zoom == *level, level.label()).clicked() {
                        self.zoom = *level;
                    }
                }

                ui.separator();

                // Navigation
                if ui.button("⏮").on_hover_text("First SV").clicked() {
                    if !self.sv_events.is_empty() {
                        self.selected_sv = Some(0);
                        self.center_on_sv(0);
                    }
                }
                if ui.button("◀").on_hover_text("Previous SV").clicked() {
                    if let Some(idx) = self.selected_sv {
                        let new_idx = if idx > 0 { idx - 1 } else { self.sv_events.len().saturating_sub(1) };
                        self.selected_sv = Some(new_idx);
                        self.center_on_sv(new_idx);
                    }
                }
                if ui.button("▶").on_hover_text("Next SV").clicked() {
                    let new_idx = self.selected_sv.map(|i| (i + 1) % self.sv_events.len()).unwrap_or(0);
                    self.selected_sv = Some(new_idx);
                    self.center_on_sv(new_idx);
                }
                if ui.button("⏭").on_hover_text("Last SV").clicked() {
                    if !self.sv_events.is_empty() {
                        let last = self.sv_events.len() - 1;
                        self.selected_sv = Some(last);
                        self.center_on_sv(last);
                    }
                }

                ui.separator();

                // Position display
                let (vis_start, vis_end) = self.visible_range(800.0);
                ui.label(format!(
                    "{}:{}-{} ({} bp)",
                    self.region.chrom,
                    vis_start,
                    vis_end,
                    vis_end - vis_start
                ));
            });
        });

        // Left panel - SV list
        if self.show_sv_panel {
            egui::SidePanel::left("sv_panel")
                .default_width(280.0)
                .show(ctx, |ui| {
                    self.render_sv_panel(ui);
                });
        }

        // Right panel - metrics and assembly comparison
        if self.show_metrics_panel {
            egui::SidePanel::right("metrics_panel")
                .default_width(320.0)
                .show(ctx, |ui| {
                    self.render_metrics_panel(ui);
                });
        }

        // Bottom panel - controls and legend
        egui::TopBottomPanel::bottom("controls").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.checkbox(&mut self.show_sv_panel, "SV Panel");
                ui.checkbox(&mut self.show_metrics_panel, "Metrics");
                ui.separator();
                ui.checkbox(&mut self.show_insertions, "Insertions");
                ui.checkbox(&mut self.show_deletions, "Deletions");
                ui.checkbox(&mut self.highlight_mismatches, "Mismatches");
                ui.checkbox(&mut self.color_by_haplotype, "Haplotypes");
                ui.checkbox(&mut self.show_depth_track, "Depth");
                ui.checkbox(&mut self.show_confidence_track, "Confidence");
                ui.separator();
                let hap_label = if let Some(ref ha) = self.haplotype_assembly {
                    format!("Per-Haplotype Assembly (H1:{} H2:{})", ha.hap1_read_count, ha.hap2_read_count)
                } else {
                    "Per-Haplotype Assembly".to_string()
                };
                ui.checkbox(&mut self.show_haplotype_assemblies, hap_label);
            });
        });

        // Central panel - main visualization
        egui::CentralPanel::default().show(ctx, |ui| {
            self.render_main_view(ui);
        });

        // Handle keyboard shortcuts
        ctx.input(|i| {
            if i.key_pressed(egui::Key::ArrowLeft) {
                self.view_center -= 50.0 * self.bases_per_pixel();
            }
            if i.key_pressed(egui::Key::ArrowRight) {
                self.view_center += 50.0 * self.bases_per_pixel();
            }
            if i.key_pressed(egui::Key::Minus) || i.key_pressed(egui::Key::Num1) {
                self.zoom = match self.zoom {
                    ZoomLevel::BasePair => ZoomLevel::Fine,
                    ZoomLevel::Fine => ZoomLevel::Medium,
                    ZoomLevel::Medium => ZoomLevel::Coarse,
                    ZoomLevel::Coarse => ZoomLevel::Region,
                    ZoomLevel::Region => ZoomLevel::Region,
                };
            }
            if i.key_pressed(egui::Key::Equals) || i.key_pressed(egui::Key::Num2) {
                self.zoom = match self.zoom {
                    ZoomLevel::Region => ZoomLevel::Coarse,
                    ZoomLevel::Coarse => ZoomLevel::Medium,
                    ZoomLevel::Medium => ZoomLevel::Fine,
                    ZoomLevel::Fine => ZoomLevel::BasePair,
                    ZoomLevel::BasePair => ZoomLevel::BasePair,
                };
            }
            if i.key_pressed(egui::Key::N) {
                let new_idx = self.selected_sv.map(|i| (i + 1) % self.sv_events.len()).unwrap_or(0);
                self.selected_sv = Some(new_idx);
                self.center_on_sv(new_idx);
            }
            if i.key_pressed(egui::Key::P) {
                if let Some(idx) = self.selected_sv {
                    let new_idx = if idx > 0 { idx - 1 } else { self.sv_events.len().saturating_sub(1) };
                    self.selected_sv = Some(new_idx);
                    self.center_on_sv(new_idx);
                }
            }
        });
    }
}

impl SVReviewerApp {
    fn render_sv_panel(&mut self, ui: &mut egui::Ui) {
        ui.heading("Structural Variants");
        ui.separator();

        // Summary stats
        let del_count = self.sv_events.iter().filter(|e| e.sv_type == SVType::Deletion).count();
        let ins_count = self.sv_events.iter().filter(|e| e.sv_type == SVType::Insertion).count();
        let total_del: i64 = self.sv_events.iter()
            .filter(|e| e.sv_type == SVType::Deletion)
            .map(|e| e.size.abs())
            .sum();
        let total_ins: i64 = self.sv_events.iter()
            .filter(|e| e.sv_type == SVType::Insertion)
            .map(|e| e.size.abs())
            .sum();

        ui.horizontal(|ui| {
            ui.colored_label(COLOR_DELETION, format!("{} DEL (-{}bp)", del_count, total_del));
            ui.colored_label(COLOR_INSERTION, format!("{} INS (+{}bp)", ins_count, total_ins));
        });

        ui.separator();

        // SV list - collect click actions to process after iteration
        let mut clicked_sv: Option<usize> = None;
        let mut hovered_sv_new: Option<usize> = None;

        egui::ScrollArea::vertical().show(ui, |ui| {
            for (idx, sv) in self.sv_events.iter().enumerate() {
                let is_selected = self.selected_sv == Some(idx);
                let is_hovered = self.hovered_sv == Some(idx);

                let bg_color = if is_selected {
                    Color32::from_rgb(60, 60, 80)
                } else if is_hovered {
                    Color32::from_rgb(50, 50, 60)
                } else {
                    Color32::TRANSPARENT
                };

                egui::Frame::new()
                    .fill(bg_color)
                    .inner_margin(4.0)
                    .corner_radius(4.0)
                    .show(ui, |ui| {
                        let response = ui.horizontal(|ui| {
                            // SV type indicator
                            let symbol = match sv.sv_type {
                                SVType::Deletion => "▼",
                                SVType::Insertion => "▲",
                            };
                            ui.colored_label(sv.color(), RichText::new(symbol).size(20.0));

                            ui.vertical(|ui| {
                                ui.horizontal(|ui| {
                                    ui.colored_label(
                                        sv.color(),
                                        RichText::new(format!("{} {}bp", sv.label(), sv.size.abs())).strong(),
                                    );
                                });
                                ui.label(
                                    RichText::new(format!("chr17:{}-{}", sv.start, sv.end))
                                        .small()
                                        .color(Color32::GRAY),
                                );
                                ui.horizontal(|ui| {
                                    ui.label(
                                        RichText::new(format!("{} reads", sv.supporting_reads.len()))
                                            .small(),
                                    );
                                    // Confidence bar
                                    let conf_width = 60.0;
                                    let (rect, _) = ui.allocate_exact_size(Vec2::new(conf_width, 8.0), Sense::hover());
                                    let painter = ui.painter();
                                    painter.rect_filled(rect, 2.0, Color32::from_gray(40));
                                    let filled = Rect::from_min_size(
                                        rect.min,
                                        Vec2::new((conf_width * sv.confidence as f32).max(2.0), 8.0),
                                    );
                                    painter.rect_filled(filled, 2.0, sv.color());
                                });
                            });
                        });

                        if response.response.clicked() {
                            clicked_sv = Some(idx);
                        }
                        if response.response.hovered() {
                            hovered_sv_new = Some(idx);
                        }
                    });
            }
        });

        // Process click actions after iteration
        if let Some(idx) = clicked_sv {
            self.selected_sv = Some(idx);
            self.center_on_sv(idx);
        }
        if hovered_sv_new.is_some() {
            self.hovered_sv = hovered_sv_new;
        }
    }

    fn render_metrics_panel(&mut self, ui: &mut egui::Ui) {
        ui.heading("Assembly Metrics");
        ui.separator();

        // Assembly selector
        if !self.assemblies.is_empty() {
            ui.horizontal(|ui| {
                ui.label("Method:");
                egui::ComboBox::from_id_salt("assembly_selector")
                    .selected_text(&self.assemblies[self.selected_assembly].name)
                    .show_ui(ui, |ui| {
                        for (idx, asm) in self.assemblies.iter().enumerate() {
                            if ui.selectable_label(self.selected_assembly == idx, &asm.name).clicked() {
                                self.selected_assembly = idx;
                            }
                        }
                    });
            });

            if let Some(asm) = self.assemblies.get(self.selected_assembly) {
                ui.separator();

                // Fitness scores
                ui.heading("Fitness Scores");
                Self::render_metric_bar(ui, "Overall", asm.fitness.overall, Color32::from_rgb(100, 200, 100));
                Self::render_metric_bar(ui, "Agreement", asm.fitness.mean_agreement, Color32::from_rgb(100, 150, 255));
                Self::render_metric_bar(ui, "Ref Identity", asm.fitness.reference_identity, Color32::from_rgb(255, 200, 100));

                ui.separator();

                ui.label(format!("Mean Depth: {:.1}", asm.fitness.mean_depth));
                ui.label(format!("Assembly Length: {} bp", asm.result.sequence.len()));
            }
        }

        ui.separator();
        ui.heading("Assembly Comparison");
        ui.checkbox(&mut self.show_assembly_comparison, "Show all methods");

        if self.show_assembly_comparison && !self.assemblies.is_empty() {
            // Comparison chart
            let bars: Vec<Bar> = self.assemblies
                .iter()
                .enumerate()
                .map(|(i, asm)| {
                    Bar::new(i as f64, asm.fitness.overall)
                        .name(&asm.name)
                        .fill(if i == self.selected_assembly {
                            Color32::from_rgb(100, 200, 100)
                        } else {
                            Color32::from_rgb(80, 80, 120)
                        })
                })
                .collect();

            Plot::new("assembly_comparison")
                .height(120.0)
                .show_axes(true)
                .legend(Legend::default())
                .show(ui, |plot_ui| {
                    plot_ui.bar_chart(BarChart::new(bars).name("Fitness"));
                });
        }

        // Haplotype assembly summary
        if let Some(ref hap_asm) = self.haplotype_assembly {
            ui.separator();
            ui.heading("Haplotype Assembly");

            egui::Frame::new()
                .fill(Color32::from_rgb(35, 35, 45))
                .inner_margin(8.0)
                .corner_radius(4.0)
                .show(ui, |ui| {
                    ui.horizontal(|ui| {
                        ui.colored_label(COLOR_HAP1, RichText::new(format!("H1: {} reads", hap_asm.hap1_read_count)).strong());
                        ui.colored_label(COLOR_HAP2, RichText::new(format!("H2: {} reads", hap_asm.hap2_read_count)).strong());
                    });

                    let total_divergent: usize = hap_asm.divergence.iter().filter(|&&d| d > 0.0).count();
                    let total_bases = hap_asm.divergence.len();
                    let pct_divergent = if total_bases > 0 {
                        total_divergent as f64 / total_bases as f64 * 100.0
                    } else {
                        0.0
                    };
                    ui.label(format!("Divergent positions: {} ({:.2}%)", total_divergent, pct_divergent));

                    let high_sv: usize = hap_asm.sv_likelihood.iter().filter(|&&s| s > 0.5).count();
                    let sv_color = if high_sv > 100 {
                        Color32::from_rgb(255, 60, 60)
                    } else if high_sv > 0 {
                        Color32::from_rgb(255, 180, 60)
                    } else {
                        Color32::from_rgb(100, 200, 100)
                    };
                    ui.colored_label(
                        sv_color,
                        format!("High-confidence SV bases: {}", high_sv),
                    );
                });
        }

        // Selected SV details
        if let Some(idx) = self.selected_sv {
            if let Some(sv) = self.sv_events.get(idx) {
                ui.separator();
                ui.heading("Selected SV");

                egui::Frame::new()
                    .fill(Color32::from_rgb(40, 40, 50))
                    .inner_margin(8.0)
                    .corner_radius(4.0)
                    .show(ui, |ui| {
                        ui.colored_label(
                            sv.color(),
                            RichText::new(format!("{} {}bp", sv.label(), sv.size.abs())).heading(),
                        );
                        ui.label(format!("Position: chr17:{}-{}", sv.start, sv.end));
                        ui.label(format!("Supporting reads: {}", sv.supporting_reads.len()));
                        ui.label(format!("Confidence: {:.1}%", sv.confidence * 100.0));

                        ui.separator();
                        ui.label("Supporting reads:");
                        egui::ScrollArea::vertical().max_height(100.0).show(ui, |ui| {
                            for read_name in &sv.supporting_reads {
                                let short_name = if read_name.len() > 25 {
                                    format!("{}...", &read_name[..22])
                                } else {
                                    read_name.clone()
                                };
                                ui.label(RichText::new(short_name).small().monospace());
                            }
                        });
                    });
            }
        }
    }

    fn render_metric_bar(ui: &mut egui::Ui, label: &str, value: f64, color: Color32) {
        ui.horizontal(|ui| {
            ui.label(format!("{}: ", label));
            let width = 100.0;
            let height = 16.0;
            let (rect, _) = ui.allocate_exact_size(Vec2::new(width, height), Sense::hover());
            let painter = ui.painter();

            // Background
            painter.rect_filled(rect, 3.0, Color32::from_gray(40));

            // Filled portion
            let filled_width = (width * value as f32).max(1.0);
            let filled = Rect::from_min_size(rect.min, Vec2::new(filled_width, height));
            painter.rect_filled(filled, 3.0, color);

            // Text
            painter.text(
                rect.center(),
                egui::Align2::CENTER_CENTER,
                format!("{:.1}%", value * 100.0),
                FontId::proportional(12.0),
                Color32::WHITE,
            );
        });
    }

    fn render_main_view(&mut self, ui: &mut egui::Ui) {
        let available_rect = ui.available_rect_before_wrap();
        let (vis_start, vis_end) = self.visible_range(available_rect.width());

        // Handle pan/zoom with mouse
        let response = ui.allocate_rect(available_rect, Sense::click_and_drag());

        if response.dragged() {
            let delta = response.drag_delta();
            self.view_center -= delta.x as f64 * self.bases_per_pixel();
            self.scroll_y -= delta.y;
        }

        // Zoom with scroll wheel
        ui.input(|i| {
            let scroll = i.raw_scroll_delta.y;
            if scroll != 0.0 {
                let factor = if scroll > 0.0 { 0.9 } else { 1.1 };
                self.custom_zoom *= factor;
                self.custom_zoom = self.custom_zoom.clamp(0.05, 100.0);
            }
        });

        let painter = ui.painter();

        // Calculate track positions
        let mut y = available_rect.top() + 10.0;
        let track_height = 20.0;

        // 1. Position ruler
        self.render_ruler(painter, available_rect, y, vis_start, vis_end);
        y += track_height + 5.0;

        // 2. Depth track
        if self.show_depth_track {
            self.render_depth_track(painter, available_rect, y, vis_start, vis_end);
            y += track_height * 2.0 + 5.0;
        }

        // 3. SV overview track
        self.render_sv_overview(painter, available_rect, y, vis_start, vis_end);
        y += track_height + 10.0;

        // 4. Reference sequence
        self.render_reference_track(painter, available_rect, y, vis_start, vis_end);
        y += track_height + 5.0;

        // 5. Assembly track
        if let Some(asm) = self.assemblies.get(self.selected_assembly) {
            self.render_assembly_track(painter, available_rect, y, vis_start, vis_end, asm);
            y += track_height + 5.0;
        }

        // 6. Confidence track
        if self.show_confidence_track {
            if let Some(asm) = self.assemblies.get(self.selected_assembly) {
                self.render_confidence_track(painter, available_rect, y, vis_start, vis_end, asm);
                y += track_height + 5.0;
            }
        }

        // 7. Per-haplotype assembly tracks
        if self.show_haplotype_assemblies {
            if let Some(ref hap_asm) = self.haplotype_assembly {
                // Hap1 assembly track
                self.render_haplotype_assembly_track(
                    painter, available_rect, y, vis_start, vis_end,
                    &hap_asm.hap1_assembly, "H1", COLOR_HAP1, hap_asm.hap1_read_count,
                );
                y += track_height + 2.0;

                // Hap2 assembly track
                self.render_haplotype_assembly_track(
                    painter, available_rect, y, vis_start, vis_end,
                    &hap_asm.hap2_assembly, "H2", COLOR_HAP2, hap_asm.hap2_read_count,
                );
                y += track_height + 2.0;

                // Divergence / SV likelihood track
                self.render_divergence_track(
                    painter, available_rect, y, vis_start, vis_end,
                    &hap_asm.divergence, &hap_asm.sv_likelihood,
                );
                y += track_height + 5.0;
            }
        }

        // Separator
        painter.line_segment(
            [Pos2::new(available_rect.left(), y), Pos2::new(available_rect.right(), y)],
            Stroke::new(1.0, Color32::from_gray(60)),
        );
        y += 5.0;

        // 7. Read pileup
        let pileup_rect = Rect::from_min_max(
            Pos2::new(available_rect.left(), y),
            available_rect.max,
        );
        self.render_read_pileup(painter, pileup_rect, vis_start, vis_end);
    }

    fn render_ruler(&self, painter: &egui::Painter, rect: Rect, y: f32, vis_start: u64, vis_end: u64) {
        let view_rect = Rect::from_min_max(
            Pos2::new(rect.left(), y),
            Pos2::new(rect.right(), y + 20.0),
        );

        painter.rect_filled(view_rect, 0.0, Color32::from_rgb(30, 30, 40));

        // Determine tick spacing based on zoom
        let range = vis_end - vis_start;
        let tick_spacing = if range < 200 { 10 } else if range < 1000 { 50 } else if range < 5000 { 200 } else { 1000 };

        let first_tick = (vis_start / tick_spacing) * tick_spacing;
        let mut pos = first_tick;

        while pos <= vis_end {
            let x = self.pos_to_x(pos, view_rect);
            if x >= view_rect.left() && x <= view_rect.right() {
                painter.line_segment(
                    [Pos2::new(x, y), Pos2::new(x, y + 8.0)],
                    Stroke::new(1.0, Color32::GRAY),
                );
                if pos % (tick_spacing * 5) == 0 || tick_spacing >= 200 {
                    painter.text(
                        Pos2::new(x, y + 12.0),
                        egui::Align2::CENTER_TOP,
                        format!("{}", pos),
                        FontId::monospace(10.0),
                        Color32::GRAY,
                    );
                }
            }
            pos += tick_spacing;
        }
    }

    fn render_depth_track(&self, painter: &egui::Painter, rect: Rect, y: f32, vis_start: u64, vis_end: u64) {
        let height = 40.0;
        let view_rect = Rect::from_min_max(
            Pos2::new(rect.left(), y),
            Pos2::new(rect.right(), y + height),
        );

        painter.rect_filled(view_rect, 0.0, Color32::from_rgb(25, 25, 35));

        // Label
        painter.text(
            Pos2::new(rect.left() + 5.0, y + 5.0),
            egui::Align2::LEFT_TOP,
            "Depth",
            FontId::proportional(10.0),
            Color32::GRAY,
        );

        // Draw depth bars
        let num_bins = (rect.width() / 2.0) as usize;
        let bin_width = rect.width() / num_bins as f32;

        for i in 0..num_bins {
            let bin_start = vis_start + (i as u64 * (vis_end - vis_start) / num_bins as u64);
            let bin_end = vis_start + ((i + 1) as u64 * (vis_end - vis_start) / num_bins as u64);

            let start_idx = bin_start.saturating_sub(self.region.start) as usize;
            let end_idx = (bin_end.saturating_sub(self.region.start) as usize).min(self.depth_by_position.len());

            if start_idx < self.depth_by_position.len() && start_idx < end_idx {
                let avg_depth: f64 = self.depth_by_position[start_idx..end_idx]
                    .iter()
                    .map(|&d| d as f64)
                    .sum::<f64>()
                    / (end_idx - start_idx).max(1) as f64;

                let bar_height = (avg_depth / self.max_depth as f64 * (height as f64 - 15.0)) as f32;
                let color = if avg_depth < 5.0 {
                    Color32::from_rgb(200, 60, 60)
                } else if avg_depth < 15.0 {
                    Color32::from_rgb(200, 200, 60)
                } else {
                    Color32::from_rgb(60, 180, 60)
                };

                let bar_rect = Rect::from_min_max(
                    Pos2::new(rect.left() + i as f32 * bin_width, y + height - bar_height),
                    Pos2::new(rect.left() + (i + 1) as f32 * bin_width - 1.0, y + height),
                );
                painter.rect_filled(bar_rect, 0.0, color);
            }
        }
    }

    fn render_sv_overview(&self, painter: &egui::Painter, rect: Rect, y: f32, vis_start: u64, vis_end: u64) {
        let height = 24.0;
        let view_rect = Rect::from_min_max(
            Pos2::new(rect.left(), y),
            Pos2::new(rect.right(), y + height),
        );

        painter.rect_filled(view_rect, 0.0, Color32::from_rgb(30, 30, 35));

        // Label
        painter.text(
            Pos2::new(rect.left() + 5.0, y + 2.0),
            egui::Align2::LEFT_TOP,
            "SVs",
            FontId::proportional(10.0),
            Color32::GRAY,
        );

        for (idx, sv) in self.sv_events.iter().enumerate() {
            if sv.end >= vis_start && sv.start <= vis_end {
                let x_start = self.pos_to_x(sv.start, view_rect).max(view_rect.left());
                let x_end = self.pos_to_x(sv.end, view_rect).min(view_rect.right());
                let width = (x_end - x_start).max(8.0);

                let is_selected = self.selected_sv == Some(idx);
                let is_hovered = self.hovered_sv == Some(idx);

                let base_color = sv.color();
                let color = if is_selected {
                    Color32::from_rgb(255, 255, 100)
                } else if is_hovered {
                    Color32::from_rgb(
                        base_color.r().saturating_add(40),
                        base_color.g().saturating_add(40),
                        base_color.b().saturating_add(40),
                    )
                } else {
                    base_color
                };

                let sv_rect = Rect::from_min_max(
                    Pos2::new(x_start, y + 4.0),
                    Pos2::new(x_start + width, y + height - 4.0),
                );

                painter.rect_filled(sv_rect, 3.0, color);

                // Draw symbol
                let symbol = match sv.sv_type {
                    SVType::Deletion => "▼",
                    SVType::Insertion => "▲",
                };
                if width > 15.0 {
                    painter.text(
                        sv_rect.center(),
                        egui::Align2::CENTER_CENTER,
                        symbol,
                        FontId::proportional(12.0),
                        Color32::WHITE,
                    );
                }

                // Size label if room
                if width > 50.0 {
                    painter.text(
                        Pos2::new(sv_rect.right() + 3.0, sv_rect.center().y),
                        egui::Align2::LEFT_CENTER,
                        format!("{}bp", sv.size.abs()),
                        FontId::proportional(10.0),
                        color,
                    );
                }
            }
        }
    }

    fn render_reference_track(&self, painter: &egui::Painter, rect: Rect, y: f32, vis_start: u64, vis_end: u64) {
        let height = 20.0;
        let view_rect = Rect::from_min_max(
            Pos2::new(rect.left(), y),
            Pos2::new(rect.right(), y + height),
        );

        painter.rect_filled(view_rect, 0.0, Color32::from_rgb(20, 40, 20));

        // Label
        painter.text(
            Pos2::new(rect.left() + 5.0, y + 2.0),
            egui::Align2::LEFT_TOP,
            "Ref",
            FontId::proportional(10.0),
            Color32::from_rgb(100, 200, 100),
        );

        // Show bases if zoomed in enough
        let bases_per_pixel = self.bases_per_pixel();
        if bases_per_pixel < 0.5 {
            // Show individual bases
            for pos in vis_start..=vis_end {
                let idx = (pos - self.region.start) as usize;
                if idx < self.reference.len() {
                    let base = self.reference[idx] as char;
                    let x = self.pos_to_x(pos, view_rect);
                    let color = base_color(self.reference[idx]);
                    painter.text(
                        Pos2::new(x, y + height / 2.0),
                        egui::Align2::CENTER_CENTER,
                        base.to_string(),
                        FontId::monospace(12.0),
                        color,
                    );
                }
            }
        } else {
            // Show colored blocks
            let block_size = (bases_per_pixel * 2.0).max(1.0) as usize;
            let mut pos = vis_start;
            while pos < vis_end {
                let idx = (pos - self.region.start) as usize;
                if idx < self.reference.len() {
                    let x = self.pos_to_x(pos, view_rect);
                    let x_end = self.pos_to_x(pos + block_size as u64, view_rect);
                    let block_rect = Rect::from_min_max(
                        Pos2::new(x, y + 4.0),
                        Pos2::new(x_end, y + height - 4.0),
                    );
                    painter.rect_filled(block_rect, 1.0, Color32::from_rgb(60, 100, 60));
                }
                pos += block_size as u64;
            }
        }
    }

    fn render_assembly_track(
        &self,
        painter: &egui::Painter,
        rect: Rect,
        y: f32,
        vis_start: u64,
        vis_end: u64,
        asm: &AssemblyInfo,
    ) {
        let height = 20.0;
        let view_rect = Rect::from_min_max(
            Pos2::new(rect.left(), y),
            Pos2::new(rect.right(), y + height),
        );

        painter.rect_filled(view_rect, 0.0, Color32::from_rgb(20, 30, 50));

        // Label
        painter.text(
            Pos2::new(rect.left() + 5.0, y + 2.0),
            egui::Align2::LEFT_TOP,
            &asm.name,
            FontId::proportional(10.0),
            Color32::from_rgb(100, 150, 255),
        );

        let bases_per_pixel = self.bases_per_pixel();
        if bases_per_pixel < 0.5 {
            // Show individual bases with mismatch highlighting
            for pos in vis_start..=vis_end {
                let idx = (pos - self.region.start) as usize;
                if idx < asm.result.sequence.len() && idx < self.reference.len() {
                    let asm_base = asm.result.sequence[idx];
                    let ref_base = self.reference[idx];
                    let is_variant = !asm_base.eq_ignore_ascii_case(&ref_base);

                    let x = self.pos_to_x(pos, view_rect);

                    if is_variant && self.highlight_mismatches {
                        // Highlight variant
                        let bg_rect = Rect::from_center_size(
                            Pos2::new(x, y + height / 2.0),
                            Vec2::new(10.0, height - 4.0),
                        );
                        painter.rect_filled(bg_rect, 2.0, Color32::from_rgb(180, 60, 180));
                    }

                    let color = if is_variant { Color32::WHITE } else { base_color(asm_base) };
                    painter.text(
                        Pos2::new(x, y + height / 2.0),
                        egui::Align2::CENTER_CENTER,
                        (asm_base as char).to_string(),
                        FontId::monospace(12.0),
                        color,
                    );
                }
            }
        }
    }

    fn render_confidence_track(
        &self,
        painter: &egui::Painter,
        rect: Rect,
        y: f32,
        vis_start: u64,
        vis_end: u64,
        asm: &AssemblyInfo,
    ) {
        let height = 20.0;
        let view_rect = Rect::from_min_max(
            Pos2::new(rect.left(), y),
            Pos2::new(rect.right(), y + height),
        );

        painter.rect_filled(view_rect, 0.0, Color32::from_rgb(25, 25, 35));

        // Label
        painter.text(
            Pos2::new(rect.left() + 5.0, y + 2.0),
            egui::Align2::LEFT_TOP,
            "Conf",
            FontId::proportional(10.0),
            Color32::GRAY,
        );

        // Draw confidence as colored bars
        let num_bins = (rect.width() / 3.0) as usize;
        let bin_width = rect.width() / num_bins as f32;

        for i in 0..num_bins {
            let bin_start = vis_start + (i as u64 * (vis_end - vis_start) / num_bins as u64);
            let bin_end = vis_start + ((i + 1) as u64 * (vis_end - vis_start) / num_bins as u64);

            let start_idx = bin_start.saturating_sub(self.region.start) as usize;
            let end_idx = (bin_end.saturating_sub(self.region.start) as usize).min(asm.result.confidence.len());

            if start_idx < asm.result.confidence.len() && start_idx < end_idx {
                let avg_conf: f64 = asm.result.confidence[start_idx..end_idx]
                    .iter()
                    .sum::<f64>()
                    / (end_idx - start_idx).max(1) as f64;

                let bar_height = (avg_conf * (height - 6.0) as f64) as f32;
                let color = if avg_conf < 0.5 {
                    Color32::from_rgb(200, 60, 60)
                } else if avg_conf < 0.8 {
                    Color32::from_rgb(200, 200, 60)
                } else {
                    Color32::from_rgb(60, 180, 60)
                };

                let bar_rect = Rect::from_min_max(
                    Pos2::new(rect.left() + i as f32 * bin_width, y + height - bar_height - 3.0),
                    Pos2::new(rect.left() + (i + 1) as f32 * bin_width - 1.0, y + height - 3.0),
                );
                painter.rect_filled(bar_rect, 0.0, color);
            }
        }
    }

    fn render_read_pileup(&self, painter: &egui::Painter, rect: Rect, vis_start: u64, vis_end: u64) {
        let read_height = 12.0;
        let read_spacing = 2.0;

        painter.rect_filled(rect, 0.0, Color32::from_rgb(20, 20, 25));

        // Label
        painter.text(
            Pos2::new(rect.left() + 5.0, rect.top() + 2.0),
            egui::Align2::LEFT_TOP,
            format!("Reads ({})", self.reads.len()),
            FontId::proportional(10.0),
            Color32::GRAY,
        );

        let max_visible_reads = ((rect.height() - 20.0) / (read_height + read_spacing)) as usize;
        let scroll_offset = (self.scroll_y / (read_height + read_spacing)).max(0.0) as usize;

        let mut y = rect.top() + 20.0;

        for read in self.reads.iter().skip(scroll_offset).take(max_visible_reads) {
            // Skip reads outside visible range
            if read.end < vis_start || read.start > vis_end {
                y += read_height + read_spacing;
                continue;
            }

            // Get haplotype color
            let hap = self.haplotype_assignments.get(&read.name);
            let hap_color = if self.color_by_haplotype {
                match hap {
                    Some(HaplotypeLabel::Hap1) => COLOR_HAP1,
                    Some(HaplotypeLabel::Hap2) => COLOR_HAP2,
                    _ => COLOR_UNASSIGNED,
                }
            } else {
                Color32::from_gray(120)
            };

            // Draw read bar
            let x_start = self.pos_to_x(read.start.max(vis_start), rect).max(rect.left());
            let x_end = self.pos_to_x(read.end.min(vis_end), rect).min(rect.right());

            if x_end > x_start {
                let read_rect = Rect::from_min_max(
                    Pos2::new(x_start, y),
                    Pos2::new(x_end, y + read_height),
                );

                painter.rect_filled(read_rect, 2.0, hap_color.gamma_multiply(0.4));

                // Draw CIGAR features
                self.render_read_cigar(painter, read, rect, y, read_height, vis_start, vis_end, hap_color);
            }

            y += read_height + read_spacing;
        }
    }

    fn render_read_cigar(
        &self,
        painter: &egui::Painter,
        read: &AlignedRead,
        rect: Rect,
        y: f32,
        height: f32,
        vis_start: u64,
        vis_end: u64,
        _base_color: Color32,
    ) {
        let mut ref_pos = read.start;
        let mut read_pos = 0usize;
        let bases_per_pixel = self.bases_per_pixel();

        for op in &read.cigar {
            match op {
                CigarOp::Match(n) => {
                    // Draw matches with mismatch highlighting
                    if bases_per_pixel < 1.0 && self.highlight_mismatches {
                        for i in 0..*n as usize {
                            let pos = ref_pos + i as u64;
                            if pos >= vis_start && pos <= vis_end && read_pos + i < read.sequence.len() {
                                let ref_idx = (pos - self.region.start) as usize;
                                if ref_idx < self.reference.len() {
                                    let read_base = read.sequence[read_pos + i];
                                    let reference_base = self.reference[ref_idx];
                                    if !read_base.eq_ignore_ascii_case(&reference_base) {
                                        let x = self.pos_to_x(pos, rect);
                                        let mm_rect = Rect::from_center_size(
                                            Pos2::new(x, y + height / 2.0),
                                            Vec2::new(6.0, height - 2.0),
                                        );
                                        painter.rect_filled(mm_rect, 1.0, COLOR_MISMATCH);
                                    }
                                }
                            }
                        }
                    }
                    ref_pos += *n as u64;
                    read_pos += *n as usize;
                }
                CigarOp::Insertion(n) => {
                    if self.show_insertions && *n >= 10 {
                        if ref_pos >= vis_start && ref_pos <= vis_end {
                            let x = self.pos_to_x(ref_pos, rect);
                            // Draw insertion marker
                            painter.line_segment(
                                [Pos2::new(x, y), Pos2::new(x, y + height)],
                                Stroke::new(2.0, COLOR_INSERTION),
                            );
                            if bases_per_pixel < 2.0 {
                                painter.text(
                                    Pos2::new(x + 2.0, y),
                                    egui::Align2::LEFT_TOP,
                                    format!("+{}", n),
                                    FontId::proportional(8.0),
                                    COLOR_INSERTION,
                                );
                            }
                        }
                    }
                    read_pos += *n as usize;
                }
                CigarOp::Deletion(n) => {
                    if self.show_deletions {
                        let del_start = ref_pos.max(vis_start);
                        let del_end = (ref_pos + *n as u64).min(vis_end);
                        if del_end > del_start {
                            let x_start = self.pos_to_x(del_start, rect);
                            let x_end = self.pos_to_x(del_end, rect);

                            // Draw deletion as red line
                            let del_y = y + height / 2.0;
                            painter.line_segment(
                                [Pos2::new(x_start, del_y), Pos2::new(x_end, del_y)],
                                Stroke::new(2.0, COLOR_DELETION),
                            );

                            // Size label for large deletions
                            if *n >= 30 && x_end - x_start > 20.0 {
                                painter.text(
                                    Pos2::new((x_start + x_end) / 2.0, del_y - 6.0),
                                    egui::Align2::CENTER_BOTTOM,
                                    format!("-{}", n),
                                    FontId::proportional(8.0),
                                    COLOR_DELETION,
                                );
                            }
                        }
                    }
                    ref_pos += *n as u64;
                }
                CigarOp::SoftClip(n) => {
                    read_pos += *n as usize;
                }
                CigarOp::HardClip(_) => {}
            }
        }
    }

    fn render_haplotype_assembly_track(
        &self,
        painter: &egui::Painter,
        rect: Rect,
        y: f32,
        vis_start: u64,
        vis_end: u64,
        assembly: &AssemblyResult,
        label: &str,
        color: Color32,
        read_count: usize,
    ) {
        let height = 20.0;
        let view_rect = Rect::from_min_max(
            Pos2::new(rect.left(), y),
            Pos2::new(rect.right(), y + height),
        );

        // Background with haplotype-tinted color
        painter.rect_filled(
            view_rect,
            0.0,
            Color32::from_rgb(
                (color.r() as u16 * 30 / 255) as u8,
                (color.g() as u16 * 30 / 255) as u8,
                (color.b() as u16 * 30 / 255) as u8 + 20,
            ),
        );

        // Label
        painter.text(
            Pos2::new(rect.left() + 5.0, y + 2.0),
            egui::Align2::LEFT_TOP,
            format!("{} ({})", label, read_count),
            FontId::proportional(10.0),
            color,
        );

        let bases_per_pixel = self.bases_per_pixel();
        if bases_per_pixel < 0.5 {
            // Show individual bases with mismatch highlighting
            for pos in vis_start..=vis_end {
                let idx = (pos - self.region.start) as usize;
                if idx < assembly.sequence.len() && idx < self.reference.len() {
                    let asm_base = assembly.sequence[idx];
                    let ref_base = self.reference[idx];
                    let is_variant = !asm_base.eq_ignore_ascii_case(&ref_base);

                    let x = self.pos_to_x(pos, view_rect);

                    if is_variant && self.highlight_mismatches {
                        let bg_rect = Rect::from_center_size(
                            Pos2::new(x, y + height / 2.0),
                            Vec2::new(10.0, height - 4.0),
                        );
                        painter.rect_filled(bg_rect, 2.0, color.gamma_multiply(0.6));
                    }

                    let text_color = if is_variant { Color32::WHITE } else { base_color(asm_base) };
                    painter.text(
                        Pos2::new(x, y + height / 2.0),
                        egui::Align2::CENTER_CENTER,
                        (asm_base as char).to_string(),
                        FontId::monospace(12.0),
                        text_color,
                    );
                }
            }
        } else {
            // Show colored blocks indicating depth/confidence
            let num_bins = (rect.width() / 3.0).max(1.0) as usize;
            let bin_width = rect.width() / num_bins as f32;

            for i in 0..num_bins {
                let bin_start = vis_start + (i as u64 * (vis_end - vis_start) / num_bins as u64);
                let bin_end = vis_start + ((i + 1) as u64 * (vis_end - vis_start) / num_bins as u64);

                let start_idx = bin_start.saturating_sub(self.region.start) as usize;
                let end_idx = (bin_end.saturating_sub(self.region.start) as usize).min(assembly.confidence.len());

                if start_idx < assembly.confidence.len() && start_idx < end_idx {
                    let avg_conf: f64 = assembly.confidence[start_idx..end_idx]
                        .iter()
                        .sum::<f64>()
                        / (end_idx - start_idx).max(1) as f64;

                    // Check for variants in this bin
                    let has_variant = (start_idx..end_idx).any(|idx| {
                        idx < assembly.sequence.len()
                            && idx < self.reference.len()
                            && !assembly.sequence[idx].eq_ignore_ascii_case(&self.reference[idx])
                    });

                    let bar_height = (avg_conf * (height - 6.0) as f64) as f32;
                    let bar_color = if has_variant {
                        color
                    } else {
                        color.gamma_multiply(0.3)
                    };

                    let bar_rect = Rect::from_min_max(
                        Pos2::new(rect.left() + i as f32 * bin_width, y + height - bar_height - 3.0),
                        Pos2::new(rect.left() + (i + 1) as f32 * bin_width - 1.0, y + height - 3.0),
                    );
                    painter.rect_filled(bar_rect, 0.0, bar_color);
                }
            }
        }
    }

    fn render_divergence_track(
        &self,
        painter: &egui::Painter,
        rect: Rect,
        y: f32,
        vis_start: u64,
        vis_end: u64,
        divergence: &[f64],
        sv_likelihood: &[f64],
    ) {
        let height = 20.0;
        let view_rect = Rect::from_min_max(
            Pos2::new(rect.left(), y),
            Pos2::new(rect.right(), y + height),
        );

        painter.rect_filled(view_rect, 0.0, Color32::from_rgb(30, 25, 25));

        // Label
        painter.text(
            Pos2::new(rect.left() + 5.0, y + 2.0),
            egui::Align2::LEFT_TOP,
            "H1/H2 Div",
            FontId::proportional(10.0),
            Color32::from_rgb(255, 180, 100),
        );

        // Draw divergence/SV likelihood as colored bars
        let num_bins = (rect.width() / 3.0).max(1.0) as usize;
        let bin_width = rect.width() / num_bins as f32;

        for i in 0..num_bins {
            let bin_start = vis_start + (i as u64 * (vis_end - vis_start) / num_bins as u64);
            let bin_end = vis_start + ((i + 1) as u64 * (vis_end - vis_start) / num_bins as u64);

            let start_idx = bin_start.saturating_sub(self.region.start) as usize;
            let end_idx = (bin_end.saturating_sub(self.region.start) as usize).min(divergence.len());

            if start_idx < divergence.len() && start_idx < end_idx {
                let avg_div: f64 = divergence[start_idx..end_idx]
                    .iter()
                    .sum::<f64>()
                    / (end_idx - start_idx).max(1) as f64;

                let avg_sv: f64 = sv_likelihood[start_idx..end_idx.min(sv_likelihood.len())]
                    .iter()
                    .sum::<f64>()
                    / (end_idx - start_idx).max(1) as f64;

                if avg_div > 0.0 {
                    let bar_height = (avg_div * (height - 6.0) as f64) as f32;

                    // Color from yellow (low likelihood) to bright red (high SV likelihood)
                    let color = if avg_sv > 0.5 {
                        Color32::from_rgb(255, 60, 60) // High SV likelihood - bright red
                    } else if avg_sv > 0.1 {
                        Color32::from_rgb(255, 180, 60) // Moderate - orange
                    } else {
                        Color32::from_rgb(255, 255, 100) // Low - yellow (divergent but low confidence)
                    };

                    let bar_rect = Rect::from_min_max(
                        Pos2::new(rect.left() + i as f32 * bin_width, y + height - bar_height - 3.0),
                        Pos2::new(rect.left() + (i + 1) as f32 * bin_width - 1.0, y + height - 3.0),
                    );
                    painter.rect_filled(bar_rect, 0.0, color);
                }
            }
        }
    }
}

fn base_color(b: u8) -> Color32 {
    match b.to_ascii_uppercase() {
        b'A' => Color32::from_rgb(100, 200, 100),
        b'C' => Color32::from_rgb(100, 100, 255),
        b'G' => Color32::from_rgb(255, 200, 100),
        b'T' => Color32::from_rgb(255, 100, 100),
        _ => Color32::GRAY,
    }
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
    eprintln!("Starting GUI viewer...");

    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1400.0, 900.0])
            .with_title("SV Reviewer - Longread Alignment Viewer"),
        ..Default::default()
    };

    eframe::run_native(
        "SV Reviewer",
        options,
        Box::new(|_cc| Ok(Box::new(SVReviewerApp::new(region, ref_seq, reads)))),
    )
    .map_err(|e| anyhow::anyhow!("Failed to run GUI: {}", e))?;

    Ok(())
}

