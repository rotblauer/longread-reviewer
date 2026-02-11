use anyhow::Result;

use crate::alignment::AlignedRead;
use crate::assembly::method::{AssemblyMethod, AssemblyResult};
use crate::metrics::fitness::{FitnessScore, MetricsCalculator};
use crate::region::Region;

/// Engine for running and evaluating multiple assembly methods on a region.
///
/// The engine takes a set of assembly methods, runs each one, computes
/// fitness metrics, and ranks the results. This allows automated evaluation
/// of different assembly strategies.
pub struct AssemblyEngine {
    methods: Vec<Box<dyn AssemblyMethod>>,
}

/// Result of evaluating a single assembly method.
#[derive(Debug)]
pub struct EvaluationResult {
    /// The assembly output.
    pub assembly: AssemblyResult,
    /// Fitness score for this assembly.
    pub fitness: FitnessScore,
}

impl AssemblyEngine {
    pub fn new() -> Self {
        Self {
            methods: Vec::new(),
        }
    }

    /// Register an assembly method.
    pub fn add_method(&mut self, method: Box<dyn AssemblyMethod>) {
        self.methods.push(method);
    }

    /// Run all registered assembly methods on the given reads and reference.
    /// Returns results sorted by fitness score (best first).
    pub fn evaluate_all(
        &self,
        reads: &[AlignedRead],
        reference: &[u8],
        region: &Region,
    ) -> Result<Vec<EvaluationResult>> {
        let calculator = MetricsCalculator::new();
        let mut results = Vec::new();

        for method in &self.methods {
            let assembly = method.assemble(reads, reference, region.start)?;
            let fitness = calculator.compute_fitness(&assembly, reads, reference, region.start);

            results.push(EvaluationResult { assembly, fitness });
        }

        // Sort by overall fitness score, best first
        results.sort_by(|a, b| {
            b.fitness
                .overall
                .partial_cmp(&a.fitness.overall)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        Ok(results)
    }

    /// Run a single named method.
    pub fn run_method(
        &self,
        method_name: &str,
        reads: &[AlignedRead],
        reference: &[u8],
        ref_start_pos: u64,
    ) -> Result<Option<AssemblyResult>> {
        for method in &self.methods {
            if method.name() == method_name {
                return method.assemble(reads, reference, ref_start_pos).map(Some);
            }
        }
        Ok(None)
    }

    /// List available method names.
    pub fn method_names(&self) -> Vec<&str> {
        self.methods.iter().map(|m| m.name()).collect()
    }
}

impl Default for AssemblyEngine {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::CigarOp;
    use crate::assembly::method::{ConsensusAssembly, WindowConsensusAssembly};

    fn make_read(name: &str, start: u64, seq: &[u8]) -> AlignedRead {
        let len = seq.len() as u64;
        AlignedRead {
            name: name.to_string(),
            chrom: "chr1".to_string(),
            start,
            mapq: 60,
            cigar: vec![CigarOp::Match(seq.len() as u32)],
            sequence: seq.to_vec(),
            qualities: vec![30; seq.len()],
            is_reverse: false,
            haplotype_tag: None,
            end: start + len - 1,
        }
    }

    #[test]
    fn test_engine_evaluate_all() {
        let mut engine = AssemblyEngine::new();
        engine.add_method(Box::new(ConsensusAssembly));
        engine.add_method(Box::new(WindowConsensusAssembly::new(4, 1)));

        let reference = b"ACGTACGT";
        let reads = vec![
            make_read("r1", 1, b"ACGTACGT"),
            make_read("r2", 1, b"ACGTACGT"),
        ];

        let region = Region::new("chr1", 1, 8).unwrap();
        let results = engine.evaluate_all(&reads, reference, &region).unwrap();

        assert_eq!(results.len(), 2);
        // Both should produce valid results
        for result in &results {
            assert!(!result.assembly.sequence.is_empty());
            assert!(result.fitness.overall >= 0.0);
            assert!(result.fitness.overall <= 1.0);
        }
    }

    #[test]
    fn test_engine_method_names() {
        let mut engine = AssemblyEngine::new();
        engine.add_method(Box::new(ConsensusAssembly));
        engine.add_method(Box::new(WindowConsensusAssembly::new(10, 2)));

        let names = engine.method_names();
        assert!(names.contains(&"consensus"));
        assert!(names.contains(&"window_consensus"));
    }

    #[test]
    fn test_engine_run_single_method() {
        let mut engine = AssemblyEngine::new();
        engine.add_method(Box::new(ConsensusAssembly));

        let reference = b"ACGT";
        let reads = vec![make_read("r1", 1, b"ACGT")];

        let result = engine.run_method("consensus", &reads, reference, 1).unwrap();
        assert!(result.is_some());

        let missing = engine.run_method("nonexistent", &reads, reference, 1).unwrap();
        assert!(missing.is_none());
    }
}
