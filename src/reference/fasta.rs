use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};

use crate::region::Region;

/// A loaded reference genome backed by an indexed FASTA file.
pub struct ReferenceGenome {
    /// Chromosome name -> full sequence (loaded on demand).
    sequences: HashMap<String, Vec<u8>>,
}

impl ReferenceGenome {
    /// Load all sequences from a FASTA file into memory.
    pub fn from_file(path: &Path) -> Result<Self> {
        use noodles::fasta;
        use std::io::BufReader;

        let file = std::fs::File::open(path)
            .with_context(|| format!("failed to open FASTA file: {}", path.display()))?;
        let mut reader = fasta::io::Reader::new(BufReader::new(file));
        let mut sequences = HashMap::new();

        for result in reader.records() {
            let record = result.context("failed to read FASTA record")?;
            let name = String::from_utf8_lossy(record.name()).into_owned();
            let seq: Vec<u8> = record.sequence().as_ref().to_vec();
            sequences.insert(name, seq);
        }

        Ok(Self { sequences })
    }

    /// Create a reference genome from in-memory sequences (useful for testing).
    pub fn from_sequences(seqs: HashMap<String, Vec<u8>>) -> Self {
        Self { sequences: seqs }
    }

    /// Fetch the reference sequence for the given region (1-based coordinates).
    /// Returns the bases as uppercase ASCII bytes.
    pub fn fetch(&self, region: &Region) -> Result<Vec<u8>> {
        let seq = self
            .sequences
            .get(&region.chrom)
            .with_context(|| format!("chromosome '{}' not found in reference", region.chrom))?;

        let start = (region.start as usize).saturating_sub(1); // convert 1-based to 0-based
        let end = std::cmp::min(region.end as usize, seq.len());

        if start >= seq.len() {
            anyhow::bail!(
                "region {} is beyond sequence length {}",
                region,
                seq.len()
            );
        }

        let bases: Vec<u8> = seq[start..end].iter().map(|b| b.to_ascii_uppercase()).collect();
        Ok(bases)
    }

    /// Get the list of available chromosome names.
    pub fn chromosomes(&self) -> Vec<&str> {
        self.sequences.keys().map(|s| s.as_str()).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_reference() -> ReferenceGenome {
        let mut seqs = HashMap::new();
        seqs.insert("chr1".to_string(), b"ACGTACGTACGTACGT".to_vec());
        seqs.insert("chr2".to_string(), b"TTTTAAAACCCCGGGG".to_vec());
        ReferenceGenome::from_sequences(seqs)
    }

    #[test]
    fn test_fetch_region() {
        let reference = test_reference();
        let region = Region::new("chr1", 1, 4).unwrap();
        let bases = reference.fetch(&region).unwrap();
        assert_eq!(bases, b"ACGT");
    }

    #[test]
    fn test_fetch_middle() {
        let reference = test_reference();
        let region = Region::new("chr1", 5, 8).unwrap();
        let bases = reference.fetch(&region).unwrap();
        assert_eq!(bases, b"ACGT");
    }

    #[test]
    fn test_fetch_missing_chrom() {
        let reference = test_reference();
        let region = Region::new("chrX", 1, 100).unwrap();
        assert!(reference.fetch(&region).is_err());
    }

    #[test]
    fn test_fetch_beyond_length() {
        let reference = test_reference();
        let region = Region::new("chr1", 100, 200).unwrap();
        assert!(reference.fetch(&region).is_err());
    }

    #[test]
    fn test_chromosomes() {
        let reference = test_reference();
        let chroms = reference.chromosomes();
        assert!(chroms.contains(&"chr1"));
        assert!(chroms.contains(&"chr2"));
    }

    #[test]
    fn test_from_fasta_file() {
        let dir = tempfile::tempdir().unwrap();
        let fasta_path = dir.path().join("test.fa");
        std::fs::write(&fasta_path, b">chr1\nACGTACGT\n>chr2\nTTTTAAAA\n").unwrap();

        let reference = ReferenceGenome::from_file(&fasta_path).unwrap();
        let region = Region::new("chr1", 1, 4).unwrap();
        let bases = reference.fetch(&region).unwrap();
        assert_eq!(bases, b"ACGT");
    }

    #[test]
    fn test_fetch_clamps_to_length() {
        let reference = test_reference();
        // Region extends beyond sequence length
        let region = Region::new("chr1", 14, 20).unwrap();
        let bases = reference.fetch(&region).unwrap();
        // Should return bases up to end of sequence (16 bases total, 0-indexed 13..15)
        assert!(!bases.is_empty());
        assert!(bases.len() <= 7); // Can't exceed 20-14+1 but capped at seq len
    }

    #[test]
    fn test_fetch_single_base() {
        let reference = test_reference();
        let region = Region::new("chr1", 1, 1).unwrap();
        let bases = reference.fetch(&region).unwrap();
        assert_eq!(bases, b"A");
    }

    #[test]
    fn test_fetch_entire_chromosome() {
        let reference = test_reference();
        let region = Region::new("chr2", 1, 16).unwrap();
        let bases = reference.fetch(&region).unwrap();
        assert_eq!(bases, b"TTTTAAAACCCCGGGG");
    }

    #[test]
    fn test_from_sequences() {
        let mut seqs = HashMap::new();
        seqs.insert("test".to_string(), b"GATTACA".to_vec());
        let reference = ReferenceGenome::from_sequences(seqs);

        let region = Region::new("test", 1, 7).unwrap();
        let bases = reference.fetch(&region).unwrap();
        assert_eq!(bases, b"GATTACA");
    }

    #[test]
    fn test_fetch_uppercase() {
        let mut seqs = HashMap::new();
        seqs.insert("chr1".to_string(), b"acgtacgt".to_vec()); // lowercase
        let reference = ReferenceGenome::from_sequences(seqs);

        let region = Region::new("chr1", 1, 4).unwrap();
        let bases = reference.fetch(&region).unwrap();
        assert_eq!(bases, b"ACGT"); // should be uppercase
    }

    #[test]
    fn test_from_file_missing_file() {
        let result = ReferenceGenome::from_file(std::path::Path::new("/nonexistent/file.fa"));
        assert!(result.is_err());
    }

    #[test]
    fn test_chromosomes_count() {
        let reference = test_reference();
        assert_eq!(reference.chromosomes().len(), 2);
    }

    #[test]
    fn test_from_fasta_file_multiline() {
        let dir = tempfile::tempdir().unwrap();
        let fasta_path = dir.path().join("test.fa");
        // Multi-line FASTA
        std::fs::write(&fasta_path, b">chr1\nACGT\nACGT\n>chr2\nTTTT\n").unwrap();

        let reference = ReferenceGenome::from_file(&fasta_path).unwrap();
        let region = Region::new("chr1", 1, 8).unwrap();
        let bases = reference.fetch(&region).unwrap();
        assert_eq!(bases, b"ACGTACGT");
    }
}
