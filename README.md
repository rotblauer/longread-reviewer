# longread-reviewer

A terminal-based long-read alignment viewer with local assembly generation, haplotype-aware read visualization, and per-base fitness evaluation metrics.

## Features

- **BAM/CRAM Alignment Loading**: Efficiently loads aligned reads from indexed BAM files for any genomic region
- **Interactive TUI Viewer**: IGV-like terminal visualization with colored nucleotides, mismatch highlighting, and scrollable pileup view
- **Dynamic Local Assembly**: Generate consensus assemblies on-the-fly for any region using pluggable assembly methods
- **Haplotype Assignment**: Automatically assign reads to haplotypes using HP tags or allele-based clustering
- **Per-Base Fitness Metrics**: Quantitative evaluation of assembly quality including agreement, depth, and quality scores
- **Assembly Method Evaluation**: Compare multiple assembly strategies and automatically rank them by fitness score

## Usage

### View alignments interactively

```bash
longread-reviewer view \
  --bam alignments.bam \
  --reference genome.fa \
  --region chr1:1000-2000
```

**Keyboard controls:**
- `←` / `→` — Scroll horizontally
- `↑` / `↓` — Scroll through reads
- `a` — Run assembly for current region
- `h` — Assign reads to haplotypes
- `n` — Cycle to next assembly method
- `e` — Evaluate all assembly methods
- `q` / `Esc` — Quit

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

### Evaluate assembly methods

```bash
longread-reviewer evaluate \
  --bam alignments.bam \
  --reference genome.fa \
  --region chr1:1000-2000
```

Runs all built-in assembly methods and ranks them by overall fitness score.

## Architecture

```
src/
├── alignment/     # BAM/CRAM reading, CIGAR parsing, read representation
├── assembly/      # Pluggable assembly engine and methods
│   ├── engine.rs  # Assembly runner and evaluator
│   └── method.rs  # AssemblyMethod trait + implementations
├── haplotype/     # Read-to-haplotype assignment
├── metrics/       # Per-base fitness metrics and scoring
├── reference/     # Reference genome (FASTA) loading
├── region.rs      # Genomic region parsing (chr:start-end)
└── viewer/        # TUI visualization
    ├── app.rs     # Application state and event loop
    └── render.rs  # Rendering widgets (pileup, ruler, status bar)
```

### Extending with custom assembly methods

Implement the `AssemblyMethod` trait to add new assembly strategies:

```rust
use longread_reviewer::assembly::method::{AssemblyMethod, AssemblyResult};
use longread_reviewer::alignment::AlignedRead;

struct MyAssembly;

impl AssemblyMethod for MyAssembly {
    fn name(&self) -> &str { "my_method" }

    fn assemble(&self, reads: &[AlignedRead], reference: &[u8]) -> anyhow::Result<AssemblyResult> {
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