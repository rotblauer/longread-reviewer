use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};

use crate::region::Region;

/// Metadata for a FASTA entry whose name encodes a genomic region (e.g. `chr17:100-200`).
#[derive(Debug, Clone)]
struct FragmentInfo {
    /// The base chromosome name (e.g. `chr17`).
    chrom: String,
    /// 1-based start position on the chromosome.
    start: u64,
    /// 1-based end position on the chromosome (inclusive).
    end: u64,
    /// Key into `sequences` for the actual data.
    seq_key: String,
}

/// A loaded reference genome backed by an indexed FASTA file.
///
/// Supports fragment FASTA files where sequence names encode genomic regions
/// (e.g. `>chr17:10953130-11022414`). When fetching by chromosome name with
/// genomic coordinates, matching fragments are used automatically.
pub struct ReferenceGenome {
    /// Chromosome name -> full sequence (loaded on demand).
    sequences: HashMap<String, Vec<u8>>,
    /// Fragment entries parsed from FASTA names like `chr:start-end`.
    fragments: Vec<FragmentInfo>,
}

impl ReferenceGenome {
    /// Try to parse a FASTA name like `chr17:10953130-11022414` into fragment info.
    fn parse_fragment_name(name: &str) -> Option<(String, u64, u64)> {
        let (chrom, rest) = name.split_once(':')?;
        let (start_str, end_str) = rest.split_once('-')?;
        let start: u64 = start_str.parse().ok()?;
        let end: u64 = end_str.parse().ok()?;
        if start > 0 && end >= start {
            Some((chrom.to_string(), start, end))
        } else {
            None
        }
    }

    /// Build fragment index from the current sequence names.
    fn build_fragments(sequences: &HashMap<String, Vec<u8>>) -> Vec<FragmentInfo> {
        let mut fragments = Vec::new();
        for key in sequences.keys() {
            if let Some((chrom, start, end)) = Self::parse_fragment_name(key) {
                fragments.push(FragmentInfo {
                    chrom,
                    start,
                    end,
                    seq_key: key.clone(),
                });
            }
        }
        fragments
    }

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

        let fragments = Self::build_fragments(&sequences);
        Ok(Self { sequences, fragments })
    }

    /// Create a reference genome from in-memory sequences (useful for testing).
    pub fn from_sequences(seqs: HashMap<String, Vec<u8>>) -> Self {
        let fragments = Self::build_fragments(&seqs);
        Self { sequences: seqs, fragments }
    }

    /// Fetch the reference sequence for the given region (1-based coordinates).
    /// Returns the bases as uppercase ASCII bytes.
    ///
    /// Supports both direct chromosome matches and fragment FASTA entries where
    /// the sequence name encodes a region (e.g. `chr17:10953130-11022414`).
    pub fn fetch(&self, region: &Region) -> Result<Vec<u8>> {
        // Direct match: chromosome name exists as-is.
        if let Some(seq) = self.sequences.get(&region.chrom) {
            let start = (region.start as usize).saturating_sub(1);
            let end = std::cmp::min(region.end as usize, seq.len());

            if start >= seq.len() {
                anyhow::bail!(
                    "region {} is beyond sequence length {}",
                    region,
                    seq.len()
                );
            }

            let bases: Vec<u8> = seq[start..end].iter().map(|b| b.to_ascii_uppercase()).collect();
            return Ok(bases);
        }

        // Fragment match: look for a fragment entry whose base chromosome matches
        // and whose coordinate range covers the requested region.
        for frag in &self.fragments {
            if frag.chrom == region.chrom
                && region.start >= frag.start
                && region.end <= frag.end
            {
                let seq = self.sequences.get(&frag.seq_key).unwrap();
                let offset_start = (region.start - frag.start) as usize;
                let offset_end = (region.end - frag.start + 1) as usize;
                let end = std::cmp::min(offset_end, seq.len());

                if offset_start >= seq.len() {
                    anyhow::bail!(
                        "region {} is beyond fragment {} length {}",
                        region,
                        frag.seq_key,
                        seq.len()
                    );
                }

                let bases: Vec<u8> = seq[offset_start..end]
                    .iter()
                    .map(|b| b.to_ascii_uppercase())
                    .collect();
                return Ok(bases);
            }
        }

        anyhow::bail!("chromosome '{}' not found in reference", region.chrom)
    }

    /// Get the list of available chromosome names.
    ///
    /// Includes both raw FASTA entry names and base chromosome names from fragments.
    pub fn chromosomes(&self) -> Vec<&str> {
        let mut seen = std::collections::HashSet::new();
        let mut chroms: Vec<&str> = Vec::new();
        for key in self.sequences.keys() {
            if seen.insert(key.as_str()) {
                chroms.push(key.as_str());
            }
        }
        for frag in &self.fragments {
            if seen.insert(frag.chrom.as_str()) {
                chroms.push(&frag.chrom);
            }
        }
        chroms
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
    fn test_fragment_fasta_fetch() {
        let mut seqs = HashMap::new();
        // Simulates a FASTA entry named "chr1:100-115" containing 16 bases
        seqs.insert("chr1:100-115".to_string(), b"ACGTACGTACGTACGT".to_vec());
        let reference = ReferenceGenome::from_sequences(seqs);

        // Fetch using the base chromosome name and original coordinates
        let region = Region::new("chr1", 100, 103).unwrap();
        let bases = reference.fetch(&region).unwrap();
        assert_eq!(bases, b"ACGT");

        // Fetch a middle slice
        let region = Region::new("chr1", 104, 107).unwrap();
        let bases = reference.fetch(&region).unwrap();
        assert_eq!(bases, b"ACGT");
    }

    #[test]
    fn test_fragment_chromosomes_includes_base_name() {
        let mut seqs = HashMap::new();
        seqs.insert("chr17:1000-2000".to_string(), b"AAAA".to_vec());
        let reference = ReferenceGenome::from_sequences(seqs);
        let chroms = reference.chromosomes();
        assert!(chroms.contains(&"chr17"));
        assert!(chroms.contains(&"chr17:1000-2000"));
    }

    #[test]
    fn test_fragment_out_of_range() {
        let mut seqs = HashMap::new();
        seqs.insert("chr1:100-115".to_string(), b"ACGTACGTACGTACGT".to_vec());
        let reference = ReferenceGenome::from_sequences(seqs);

        // Request outside fragment range
        let region = Region::new("chr1", 50, 60).unwrap();
        assert!(reference.fetch(&region).is_err());
    }

    #[test]
    fn test_fragment_fasta_file() {
        let dir = tempfile::tempdir().unwrap();
        let fasta_path = dir.path().join("test.fa");
        std::fs::write(&fasta_path, b">chr5:500-507\nACGTACGT\n").unwrap();

        let reference = ReferenceGenome::from_file(&fasta_path).unwrap();
        let region = Region::new("chr5", 500, 503).unwrap();
        let bases = reference.fetch(&region).unwrap();
        assert_eq!(bases, b"ACGT");
    }
}
