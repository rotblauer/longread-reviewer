use std::fmt;
use std::str::FromStr;

use thiserror::Error;

#[derive(Error, Debug)]
pub enum RegionError {
    #[error("invalid region format: expected 'chr:start-end', got '{0}'")]
    InvalidFormat(String),
    #[error("invalid coordinate: {0}")]
    InvalidCoordinate(#[from] std::num::ParseIntError),
    #[error("start ({start}) must be less than end ({end})")]
    InvalidRange { start: u64, end: u64 },
}

/// A genomic region specified as chromosome:start-end (1-based, inclusive).
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Region {
    pub chrom: String,
    pub start: u64,
    pub end: u64,
}

impl Region {
    pub fn new(chrom: impl Into<String>, start: u64, end: u64) -> Result<Self, RegionError> {
        if start > end {
            return Err(RegionError::InvalidRange { start, end });
        }
        Ok(Self {
            chrom: chrom.into(),
            start,
            end,
        })
    }

    /// Length of the region in bases.
    pub fn len(&self) -> u64 {
        self.end - self.start + 1
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl fmt::Display for Region {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}-{}", self.chrom, self.start, self.end)
    }
}

impl FromStr for Region {
    type Err = RegionError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (chrom, rest) = s
            .split_once(':')
            .ok_or_else(|| RegionError::InvalidFormat(s.to_string()))?;
        let (start_str, end_str) = rest
            .split_once('-')
            .ok_or_else(|| RegionError::InvalidFormat(s.to_string()))?;
        let start: u64 = start_str.parse()?;
        let end: u64 = end_str.parse()?;
        Region::new(chrom, start, end)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_valid_region() {
        let r: Region = "chr1:1000-2000".parse().unwrap();
        assert_eq!(r.chrom, "chr1");
        assert_eq!(r.start, 1000);
        assert_eq!(r.end, 2000);
        assert_eq!(r.len(), 1001);
    }

    #[test]
    fn test_parse_single_base() {
        let r: Region = "chrX:500-500".parse().unwrap();
        assert_eq!(r.len(), 1);
    }

    #[test]
    fn test_parse_invalid_format() {
        assert!("chr1".parse::<Region>().is_err());
        assert!("chr1:1000".parse::<Region>().is_err());
        assert!("chr1:abc-def".parse::<Region>().is_err());
    }

    #[test]
    fn test_invalid_range() {
        assert!(Region::new("chr1", 2000, 1000).is_err());
    }

    #[test]
    fn test_display() {
        let r = Region::new("chr1", 100, 200).unwrap();
        assert_eq!(r.to_string(), "chr1:100-200");
    }
}
