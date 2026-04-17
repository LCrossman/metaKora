//! # metakora
//!
//! `metakora` is a high-performance Rust tool designed to calculate **α-diversity metrics**
//! directly from frequency-of-frequencies (histogram) data.
//!
//! ## Overview
//!
//! In metagenomics and biomarker discovery, datasets often involve billions of k-mer observations.
//! Processing these as raw vectors is memory-intensive. `metakora` solves this by using a
//! histogram-based approach ($k$ abundance $\to$ $n_k$ features), allowing for constant-memory
//! complexity relative to the number of abundance classes.
//!
//! ## Supported Metrics
//!
//! This crate implements standard ecological indices, aligned with R packages like `abdiv` and `vegan`:
//!
//! * **Shannon Entropy ($H$):** Measures uncertainty and complexity using natural logarithms ($\ln$).
//! * **Chao1:** A richness estimator that predicts total unique features (including unobserved ones) using singletons and doubletons.
//! * **Pielou's Evenness ($J$):** Measures how evenly individuals are distributed across features.
//! * **Robbins Estimator:** Calculates the probability that the next sampled individual represents a previously unobserved feature.
//! * **Berger-Parker Index:** A simple measure of the numerical dominance of the most abundant feature.
//! * **Simpson Index:** Probability that two randomly selected individuals belong to different features ($1 - D$).
//!
//! ## Input File Format
//!
//! The tool expects a two-column, whitespace-separated text file:
//! ```text
//! 1  14502   # 14,502 distinct k-mers were seen exactly once (singletons)
//! 2  3200    # 3,200 distinct k-mers were seen exactly twice (doubletons)
//! 15 1       # 1 distinct k-mer was seen 15 times
//! ```
//!
//! ## Example Usage
//!
//! ```rust
//! use std::collections::HashMap;
//! use metakora::{chao1, compute_entropy_and_evenness};
//!
//! fn main() {
//!     let mut histogram = HashMap::new();
//!     histogram.insert(1, 10); // 10 singletons
//!     histogram.insert(2, 5);  // 5 doubletons
//!
//!     let richness = chao1(&histogram);
//!     let (shannon, _, _) = compute_entropy_and_evenness(&histogram);
//!
//!     println!("Estimated True Richness: {}", richness);
//! }
//! ```
#![allow(non_snake_case)]
use clap::Parser;
use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader};

///metaKora: α-Diversity Metrics from K-mer Histograms
///metaKora calculates biological diversity indices (Shannon, Chao1, etc.)
///from a two-column frequency distribution file.
#[derive(Parser, Debug)]
#[command(
    author = "LCrossman",
    version = "1.0",
    about = "Computes alpha-diversity for genomic biomarkers",
    long_about = "A high-performance estimator designed for large-scale k-mer histograms. \
                  Supports Shannon, Chao1, Pielou, Robbins, and Berger-Parker metrics."
)]
struct Arguments {
    //path to the input histogram file (format: abundance count)
    #[arg(short, long, value_name = "FILE")]
    filename: String,
}

//Error types for histogram error
#[derive(Debug)]
pub enum HistogramError {
    Io(std::io::Error),
    ParseError {
        line_number: usize,
        content: String,
        details: String,
    },
    InvalidData(String),
}
//implement display for clean error messages
impl fmt::Display for HistogramError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(e) => write!(f, "IO Error: {}", e),
            Self::ParseError {
                line_number,
                content,
                details,
            } => {
                write!(
                    f,
                    "Line {}: Could not parse '{}' - {}",
                    line_number, content, details
                )
            }
            Self::InvalidData(msg) => write!(f, "Invalid Data: {}", msg),
        }
    }
}
//allow automatic conversion of std::io::Error to HistogramError
impl From<std::io::Error> for HistogramError {
    fn from(err: std::io::Error) -> Self {
        HistogramError::Io(err)
    }
}

//computes shannnon, maximum entropy and pielou evenness
fn compute_entropy_and_evenness(histogram: &HashMap<u32, u32>) -> (f64, f64, f64) {
    // Total number of k-mers (including repeats)
    let total_kmers: f64 = histogram
        .iter()
        .map(|(&abundance, &n_k)| (abundance as f64) * (n_k as f64))
        .sum();

    // Number of distinct k-mers
    let distinct_kmers: f64 = histogram
        .iter()
        .filter(|&(abundance, _)| *abundance > 0)
        .map(|(_, &n_k)| n_k as f64)
        .sum();

    // Shannon entropy
    let mut Shannon = 0.0;
    for (&abundance, &n_k) in histogram {
        if abundance > 0 {
            let p = abundance as f64 / total_kmers;
            //using ln to match R default, TODO: provide option for log2
            Shannon -= (n_k as f64) * (p * p.ln());
        }
    }

    //maximum entropy H_max (if all distinct k-mers were equally likely)
    let H_max = if distinct_kmers > 0.0 {
        distinct_kmers.ln()
    } else {
        0.0
    };

    // Pielou’s evenness
    let Pielou = if H_max > 0.0 { Shannon / H_max } else { 0.0 };

    (Shannon, H_max, Pielou)
}
//computes the simpson index
pub fn simpson_index(histogram: &HashMap<u32, u32>) -> f64 {
    let mut total = 0.0;
    let mut numerator = 0.0;

    for (count, n_k) in histogram {
        let n_i = (*count as f64) * (*n_k as f64);
        numerator += n_i * (n_i - 1.0);
        total += n_i;
    }

    if total > 1.0 {
        1.0 - (numerator / (total * (total - 1.0)))
    } else {
        0.0
    }
}
//calculate the observed features from HashMap
pub fn observed_features(histogram: &HashMap<u32, u32>) -> u32 {
    histogram
        .iter()
        .filter(|&(abundance, _)| *abundance > 0)
        .map(|(_, &n_k)| n_k)
        .sum()
}
//calculate the chao1 from HashMap
pub fn chao1(histogram: &HashMap<u32, u32>) -> f64 {
    //s_obs is the number of species with abundance > 0
    let s_obs = histogram
        .iter()
        .filter(|&(abundance, _)| *abundance > 0)
        .map(|(_, &n_k)| n_k as f64)
        .sum::<f64>();

    //f1 and f2 are Rare species counts from the histogram
    let f1 = *histogram.get(&1).unwrap_or(&0) as f64;
    let f2 = *histogram.get(&2).unwrap_or(&0) as f64;
    //check for small or sparse samples
    if f2 > 0.0 {
        //classic chao1 formula - use if there are doubletons
        s_obs + (f1 * f1) / (2.0 * f2)
    } else {
        //bias-corrected formula for use when f2 is 0
        s_obs + (f1 * (f1 - 1.0)) / (2.0 * (f2 + 1.0))
    }
}
//calculate the robbins from HashMap
pub fn robbins(histogram: &HashMap<u32, u32>) -> f64 {
    let n1 = *histogram.get(&1).unwrap_or(&0) as f64;
    let total: f64 = histogram
        .iter()
        .map(|(&k, &n_k)| (k as f64) * (n_k as f64))
        .sum();
    if total > 0.0 {
        n1 / total
    } else {
        0.0
    }
}
//calculate the berger_parker dominance from HashMap
pub fn berger_parker_d(histogram: &HashMap<u32, u32>) -> f64 {
    let total: f64 = histogram
        .iter()
        .map(|(k, n_k)| (*k as f64) * (*n_k as f64))
        .sum();
    let mut max_class = 0.0;
    for (count, n_k) in histogram {
        let occ = (*count as f64) * (*n_k as f64);
        if occ > max_class {
            max_class = occ;
        }
    }
    if total > 0.0 {
        max_class / total
    } else {
        0.0
    }
}

pub fn read_histogram_from_file(path: &str) -> Result<HashMap<u32, u32>, HistogramError> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut histogram: HashMap<u32, u32> = HashMap::new();

    for (index, line) in reader.lines().enumerate() {
        let line_number = index + 1;
        let line = line?;
        //trim and skip empty lines
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        let parts: Vec<&str> = trimmed.split_whitespace().collect();

        // Skip empty or malformed lines
        if parts.len() != 2 {
            return Err(HistogramError::ParseError {
                line_number,
                content: line,
                details: "Expected exactly 2 columns".to_string(),
            });
        }
        //parse the abundance
        let count: u32 = parts[0].parse().map_err(|_| HistogramError::ParseError {
            line_number,
            content: line.clone(),
            details: format!("Invalid abundance key: {}", parts[0]),
        })?;
        let n_k: u32 = parts[1].parse().map_err(|_| HistogramError::ParseError {
            line_number,
            content: line.clone(),
            details: format!("Invalid frequency count: {}", parts[1]),
        })?;
        if count == 0 {
            return Err(HistogramError::InvalidData(format!(
                "Zero abundance key found at line {}",
                line_number
            )));
        }
        histogram.insert(count, n_k);
    }
    if histogram.is_empty() {
        return Err(HistogramError::InvalidData(
            "File was empty or contained no valid data".into(),
        ));
    }

    Ok(histogram)
}

fn main() -> Result<(), HistogramError> {
    //example histogram: k-mer count -> number of k-mers with that count
    let args = Arguments::parse();
    let histogram = read_histogram_from_file(&args.filename)?;

    let (Shannon, H_max, Pielou) = compute_entropy_and_evenness(&histogram);
    let c = chao1(&histogram);
    let o = observed_features(&histogram);
    let r = robbins(&histogram);
    let simp = simpson_index(&histogram);
    let bpd = berger_parker_d(&histogram);
    println!("Sample\tShannon\tH_max\tPielou\tChao1\tObserved\tRobbins\tBerger_Parker\tSimpson");
    println!(
        "{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
        &args.filename, Shannon, H_max, Pielou, c, o, r, bpd, simp
    );
    Ok(())
}
