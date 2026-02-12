# longread-reviewer

A terminal-based and GUI long-read alignment viewer with **structural variant detection**, **dynamic local assembly** (independent of the original aligner), haplotype-aware read visualization, and per-base fitness evaluation metrics.

> ⚠️ **Prototype generated via iterative AI prompting.** This repository (code, docs, and examples) was produced by iteratively prompting AI models as a rapid prototype. Review every component carefully and validate independently before relying on it for research, clinical, or production use.

## Prototype status

- Experimental, research-only prototype with no guarantees on accuracy, safety, or performance.
- Generated through iterative AI prompting; expect rough edges and simplifications typical of LLM-produced code.
- Please file issues or patches where you spot problems so others are not surprised.

## Features

- **Structural Variant Detection**: Detects large insertions and deletions from both CIGAR strings *and* assembly-vs-reference comparison for independent confirmation
- **Dynamic Local Assembly**: Multiple aligner-independent assembly methods including de Bruijn graph–based assembly, majority-vote consensus, and windowed consensus
- **Assembly-Based SV Calling**: Compares local assemblies to the reference and identifies divergent regions as candidate SVs — no reliance on the original aligner
- **Interactive IGV-like TUI Viewer**: Navigate through alignments, jump between SVs, and inspect sequences at base-level in the terminal
- **GUI Viewer**: Zoomable, pannable graphical interface with multi-track visualization, per-haplotype assembly comparison, and SV confidence scoring
- **BAM/CRAM Alignment Loading**: Efficiently loads aligned reads from indexed BAM files for any genomic region
- **Haplotype Assignment**: Automatically assign reads to haplotypes using HP tags or allele-based clustering
- **Per-Base Fitness Metrics**: Quantitative evaluation of assembly quality including agreement, depth, and quality scores
- **Assembly Method Evaluation**: Compare multiple assembly strategies and automatically rank them by fitness score

## Quick Start: Interactive TUI Viewer

```bash
cargo run --release --example interactive_viewer
```

This launches an IGV-like terminal interface with the bundled test data (NA19240 PacBio HiFi reads on chr17).

### TUI Controls

**Navigation:**
- `←` / `→` — Scroll 10bp (hold `Shift` for 100bp)
- `↑` / `↓` — Scroll through reads
- `Home` / `End` — Jump to start/end of region
- `PageUp` / `PageDown` — Scroll 20 reads

**Structural Variants:**
- `Tab` / `]` — Jump to next SV
- `Shift+Tab` / `[` — Jump to previous SV  
- `Enter` — Jump to selected SV
- `s` — Toggle SV panel

**Assembly & Analysis:**
- `a` — Run assembly for current region
- `n` — Cycle to next assembly method
- `h` — Assign reads to haplotypes

**Display:**
- `i` — Toggle insertion display
- `m` — Toggle mismatch highlighting
- `?` / `F1` — Show help
- `q` / `Esc` — Quit

## Quick Start: GUI Viewer

The GUI viewer provides a zoomable, pannable graphical interface with multi-track
visualization. It uses the bundled example data (NA19240 HiFi reads on chr17,
~59 kb region containing a complex structural event).

### 1. Launch the GUI

```bash
cargo run --release --example gui_viewer
```

The viewer opens automatically with the test data loaded. On startup it will:
- Load 206 PacBio HiFi reads for **chr17:10,958,130–11,017,414**
- Detect structural variants from CIGAR strings
- Run multiple assembly methods (consensus, windowed consensus, de Bruijn)
- Perform per-haplotype assembly and divergence analysis
- Run assembly-based SV detection (independent of the original aligner)

### 2. Navigating the GUI

**Zoom levels** (top toolbar or keyboard):
| Level | Description | Keyboard |
|-------|-------------|----------|
| Base | Single-base resolution (~100 bp visible) | `+` / `=` to zoom in |
| Fine | Fine detail (~500 bp) | |
| Medium | Standard view (~2 kb) | |
| Overview | Broad view (~10 kb) | `-` / `1` to zoom out |
| Full Region | Entire loaded region | |

**Pan and zoom:**
- **Click and drag** — Pan left/right and scroll through reads
- **Scroll wheel** — Fine zoom in/out
- **Arrow keys** — `←` / `→` to pan

**SV navigation** (top toolbar or keyboard):
- **⏮ / ⏭** — Jump to first / last SV
- **◀ / ▶** — Previous / next SV
- `N` — Next SV, `P` — Previous SV

### 3. Understanding the Tracks

The main visualization area contains several tracks from top to bottom:

| Track | Description |
|-------|-------------|
| **Position ruler** | Genomic coordinates with tick marks |
| **Depth** | Per-position read depth as a colored bar chart (red < 5×, yellow < 15×, green ≥ 15×) |
| **SVs** | Structural variant locations: ▼ = deletion (red), ▲ = insertion (green). Selected SV highlighted in yellow |
| **Ref** | Reference sequence (individual bases at high zoom, colored blocks at low zoom) |
| **Assembly** | Currently selected assembly method's consensus (variants highlighted in magenta) |
| **Conf** | Assembly confidence track (red < 50%, yellow < 80%, green ≥ 80%) |
| **H1 / H2** | Per-haplotype assembly tracks with variant highlighting |
| **H1/H2 Div** | Haplotype divergence track — bright red = high SV likelihood, orange = moderate, yellow = low confidence divergence |
| **Reads** | Read pileup colored by haplotype (blue = H1, pink = H2, gray = unassigned) with CIGAR-aware rendering |

### 4. Side Panels

**Left panel — Structural Variants:**
- Lists all detected SVs with type, size, position, supporting reads, and confidence bar
- Shows **assembly-based SV detection** summary with events detected independently from assembly–reference comparison
- Click any SV to navigate to it

**Right panel — Assembly Metrics:**
- Select and compare assembly methods (consensus, windowed consensus, de Bruijn)
- View fitness scores: overall, agreement, reference identity
- Toggle "Show all methods" for a comparison bar chart
- **Haplotype Assembly** section shows H1/H2 read counts, divergent positions, and high-confidence SV bases
- **Selected SV** detail shows:
  - ✓ **Assembly-confirmed** if the SV is corroborated by assembly-based detection
  - ⚠ **CIGAR-only** if the SV is only detected from read CIGAR strings

### 5. Display Toggles (bottom toolbar)

| Toggle | Description |
|--------|-------------|
| SV Panel | Show/hide the left SV panel |
| Metrics | Show/hide the right metrics panel |
| Insertions | Show insertion markers on reads |
| Deletions | Show deletion lines on reads |
| Mismatches | Highlight bases differing from reference |
| Haplotypes | Color reads by haplotype assignment |
| Depth | Show the depth track |
| Confidence | Show the assembly confidence track |
| Per-Haplotype Assembly | Show the H1, H2, and divergence tracks |

### 6. Interpreting SV Evidence

A **real structural variant** will show multiple converging signals:

1. **CIGAR evidence**: Multiple reads show large insertions or deletions at the same position
2. **Assembly divergence**: The local assembly differs from the reference in the same region
3. **Haplotype divergence**: The two haplotype assemblies differ (for heterozygous events)
4. **Depth anomaly**: Read depth drops at deletions or shows unusual patterns at insertions
5. **High confidence**: The assembly confidence is high despite differing from reference

When the **Selected SV** panel shows "✓ Assembly-confirmed", it means the SV was detected
both from read CIGARs and independently from the assembly comparison — this is strong
evidence that the event is real. Events marked "⚠ CIGAR-only" should be reviewed more
carefully.

## Command Line Usage

### View alignments interactively

```bash
longread-reviewer view \
  --bam alignments.bam \
  --reference genome.fa \
  --region chr1:1000-2000
```

### Generate a local assembly

```bash
longread-reviewer assemble \
  --bam alignments.bam \
  --reference genome.fa \
  --region chr1:1000-2000 \
  --method consensus
```

Available assembly methods:
- `consensus` — Simple majority-vote consensus
- `window_consensus` — Windowed consensus with configurable window size and overlap
- `debruijn` — De Bruijn graph–based local assembly (aligner-independent)

### Evaluate assembly methods

```bash
longread-reviewer evaluate \
  --bam alignments.bam \
  --reference genome.fa \
  --region chr1:1000-2000
```

Runs all built-in assembly methods (including de Bruijn) and ranks them by overall fitness score.

## Example Output

A comprehensive set of example results is auto-generated from the bundled test
data (PacBio HiFi reads for NA19240 chr17) on every push via
[GitHub Actions](.github/workflows/generate-examples.yml). You can view the
latest output in two ways:

- **Plain text:** [`examples/example_output.txt`](examples/example_output.txt)
- **HTML report:** [`docs/index.html`](docs/index.html) (also available via
  [GitHub Pages](https://rotblauer.github.io/longread-reviewer/))

To regenerate locally:

```bash
cargo run --release --example generate_output
```

This writes `examples/example_output.txt` and `docs/index.html`.

## Architecture

```
src/
├── alignment/        # BAM/CRAM reading, CIGAR parsing, read representation
├── assembly/         # Pluggable assembly engine and methods
│   ├── engine.rs     # Assembly runner and evaluator
│   ├── method.rs     # AssemblyMethod trait + implementations (consensus, windowed, de Bruijn)
│   └── sv_detect.rs  # Assembly-based SV detection (compares assembly vs reference)
├── haplotype/        # Read-to-haplotype assignment
├── metrics/          # Per-base fitness metrics and scoring
├── reference/        # Reference genome (FASTA) loading
├── region.rs         # Genomic region parsing (chr:start-end)
└── viewer/           # TUI visualization
    ├── app.rs        # Application state and event loop
    ├── render.rs     # Rendering widgets (pileup, ruler, status bar)
    └── enhanced_view.rs  # SV-aware CIGAR rendering
```

### Dynamic Local Assembly Methods

The assembly engine supports multiple methods, all of which operate on read
sequences directly rather than relying on the original aligner's coordinate
mapping:

| Method | Description |
|--------|-------------|
| `ConsensusAssembly` | Quality-weighted majority-vote at each reference position |
| `WindowConsensusAssembly` | Overlapping-window consensus for handling complex SVs |
| `DeBruijnAssembly` | De Bruijn graph–based assembly: builds a k-mer graph from reads, traverses the heaviest path, and aligns back to the reference |
| `HaplotypeAwareAssembly` | Per-haplotype assembly with divergence and SV likelihood scoring |

### Extending with custom assembly methods

Implement the `AssemblyMethod` trait to add new assembly strategies:

```rust
use longread_reviewer::assembly::method::{AssemblyMethod, AssemblyResult};
use longread_reviewer::alignment::AlignedRead;

struct MyAssembly;

impl AssemblyMethod for MyAssembly {
    fn name(&self) -> &str { "my_method" }

    fn assemble(&self, reads: &[AlignedRead], reference: &[u8], ref_start_pos: u64) -> anyhow::Result<AssemblyResult> {
        // Your assembly logic here
        todo!()
    }
}
```

Register it with the engine:

```rust
use longread_reviewer::assembly::engine::AssemblyEngine;

let mut engine = AssemblyEngine::new();
engine.add_method(Box::new(MyAssembly));
```

## AI review notes and common pitfalls

These are areas that warrant extra manual review because they reflect common LLM-generated shortcuts:

- **Strand handling in alignments** (`src/alignment/reader.rs`): helper methods ignore the `is_reverse` flag, so reverse-strand reads are not reverse-complemented before consensus/metrics, which can invert bases in pileups.
- **Quality and depth heuristics in consensus** (`src/assembly/method.rs`): consensus weighting reuses only the first quality score per read and counts depth without indel-aware adjustments, so confidence can be inflated relative to the true pileup.
- **Naïve haplotype clustering** (`src/haplotype/assign.rs`): allele-based grouping is O(N×L) over reads and bases, ignores indels/quality, and drops ties to `Unassigned`, making phasing brittle on sparse or noisy data.
- **Fitness scoring trust in assembly output** (`src/metrics/fitness.rs`): depth/confidence come from the assembly result instead of recalculating from raw reads, so scores may look better than the underlying pileup warrants if assemblies truncate or misplace bases.
- **Reference-biased de Bruijn assembly** (`src/assembly/method.rs`): the graph traversal seeds with the reference k-mer and caps length to the reference window, so real structural variants or low-coverage gaps can collapse back to the reference path with optimistic confidence.

## Building

```bash
cargo build --release
```

## Testing

```bash
cargo test
```

## License

MIT
