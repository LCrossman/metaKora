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
//! * **Inv Berger-Parker Index:** The reciprocal measure of the numerical dominance of the most abundant feature, calculated after noise removal.
//! * **Simpson Index:** Probability that two randomly selected individuals belong to different features ($1 - D$).
//! * **Peak Area:** The area under the primary biological peak of the k-mer abundance distribution.
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
//!
//! To create the berger-parker index, we have used the inverse berger-parker distance `inverse_berger_parker_d`.
//! Berger-Parker index is a measure of dominance in the dataset.  Using kmers, this was found to be very sensitive to noise.
//! Here we have created an option to pass in a parameter for a minimum abundance threshold.
//!
//! If the parameter is not provided, a minimum abundance threshold is calculated automatically from the histogram using
//! a peak detection algorithm ignoring the first peak (which is often a background noise peak).
//! If you choose to pass a manual threshold, it will be used instead of the automatically calculated value, however,
//! the script will still warn you if the peak detection algorithm identifies that the manual threshold is still within the noise peak.
//! If you are unsure about the optimal threshold, you may be best to let the script calculate it automatically.
//! Note that if you run the script on multiple samples, the automatic threshold may vary between samples, since noise levels may differ.
//!
#![allow(non_snake_case)]
use clap::Parser;
use log::{LevelFilter, info};
use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use std::io::Write;
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
                  Supports Shannon, Chao1, Pielou, Robbins, Peak Area and Berger-Parker metrics."
)]
struct Arguments {
    //path to the input histogram file (format: abundance count)
    #[arg(short, long, value_name = "FILE")]
    filename: String,
    #[arg(short, long, value_name = "MIN_ABUNDANCE")]
    min_abundance: Option<u32>,
    #[arg(short, long, default_value = "metakora.log")]
    log_file: String,
}

//setup a log file to keep details of the analysis without printing to stdout
fn setup_logging(log_path: &str) -> Result<(), std::io::Error> {
    let target = Box::new(File::create(log_path)?);
    env_logger::Builder::new()
        .format(|buf, record| {
            writeln!(
                buf,
                "{} [{}] - {}",
                chrono::Local::now().format("%Y-%m-%d %H:%M:%S"),
                record.level(),
                record.args()
            )
        })
        .filter(None, LevelFilter::Info)
        .target(env_logger::Target::Pipe(target))
        .init();
    Ok(())
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
    if total > 0.0 { n1 / total } else { 0.0 }
}
//functions to calculate a real biological peak for dominance by removing noise peak at the start of the histogram
//we calculate the peak area dominance from HashMap
pub fn peak_area_dominance_integrated(
    histogram: &HashMap<u32, u32>,
    noise_floor: Option<u32>,
) -> f64 {
    //first calculate Total Mass (N) for the entire library
    let total_mass: f64 = histogram
        .iter()
        .map(|(&k, &v)| (k as f64) * (v as f64))
        .sum();

    if total_mass <= 0.0 {
        //if total mass is zero, return a value of zero dominance
        return 0.0;
    }
    let floor = noise_floor.unwrap_or(1);
    //gatekeeper function to find a real peak in the histogram rather than just the first peak
    let peak_x = find_highest_biological_peak(histogram, floor);
    match peak_x {
        Some(x) => {
            //calculate Area Under the peak (15% width)
            let width = (x as f64 * 0.15).max(2.0) as u32;
            let start = (x as u32).saturating_sub(width);
            let end = (x as u32) + width;
            let peak_area: f64 = (start..=end)
                .map(|abundance| {
                    (abundance as f64) * (*histogram.get(&abundance).unwrap_or(&0) as f64)
                })
                .sum();
            //calculate dominance of the peak
            peak_area / total_mass
        }
        None => 0.0, //no biological peak found, therefore no dominance
    }
}
//smooth the histogram with a moving average window to eliminate micro-fluctuations
//before trying to find peaks/valleys in noisy data
fn smooth_histogram(histogram: &HashMap<u32, u32>, window: u32) -> HashMap<u32, f64> {
    let max_x = *histogram.keys().max().unwrap_or(&0);
    let half = window / 2;
    let mut smoothed = HashMap::new();
    for x in 1..=max_x {
        let lo = x.saturating_sub(half);
        let hi = (x + half).min(max_x);
        let count = (hi - lo + 1) as f64;
        let sum: f64 = (lo..=hi)
            .map(|i| *histogram.get(&i).unwrap_or(&0) as f64)
            .sum();
        smoothed.insert(x, sum / count);
    }
    info!("smoothed histogram");
    smoothed
}
//find the noise valley (after the first local peak) in the smoothed histogram
fn find_noise_valley(histogram: &HashMap<u32, u32>, smooth_window: u32) -> u32 {
    let smoothed = smooth_histogram(histogram, smooth_window);
    let max_x = *histogram.keys().max().unwrap_or(&0);

    //find the FIRST local peak (the noise/error spike, usually at low x)
    let mut noise_peak_x: Option<u32> = None;
    for x in 2..max_x {
        let prev = smoothed.get(&(x - 1)).copied().unwrap_or(0.0);
        let curr = smoothed.get(&x).copied().unwrap_or(0.0);
        let next = smoothed.get(&(x + 1)).copied().unwrap_or(0.0);
        if curr > prev && curr > next {
            noise_peak_x = Some(x);
            break; //only care about the FIRST peak
        }
    }
    let noise_peak = match noise_peak_x {
        Some(p) => p,
        None => return 1, //fallback: no clear noise peak found
    };
    //walk forward from noise peak, find the first local MINIMUM (the valley)
    for x in (noise_peak + 1)..max_x {
        let prev = smoothed.get(&(x - 1)).copied().unwrap_or(0.0);
        let curr = smoothed.get(&x).copied().unwrap_or(0.0);
        let next = smoothed.get(&(x + 1)).copied().unwrap_or(0.0);
        if curr < prev && curr < next {
            return x; //this is the min_abundance figure for the berger-parker diversity calculation
        }
    }
    info!("noise peak: {}", noise_peak);
    noise_peak //fallback if no valley found after the noise peak
}

//full pipeline: auto-detect min_abundance, then calculate diversity
pub fn calculate_diversity_auto(histogram: &HashMap<u32, u32>, min_abundance: Option<u32>) -> f64 {
    //we fix a smoothing window at 5 for consistent noise reduction
    const SMOOTH_WINDOW: u32 = 5;
    let min_abundance =
        min_abundance.unwrap_or_else(|| find_noise_valley(histogram, SMOOTH_WINDOW));
    info!("threshold min_abundance: {}", min_abundance);
    //use the biological peak finder for everything above this min_abundance floor
    let bio_peak = find_highest_biological_peak(histogram, min_abundance);
    info!("Biological peak at: {:?}", bio_peak);
    inverse_berger_parker_d(histogram, min_abundance)
}

//create a kmer peak finder to look for the "Best" peak after noise peak
fn find_highest_biological_peak(histogram: &HashMap<u32, u32>, noise_floor: u32) -> Option<u32> {
    let mut best_peak_x = None;
    let mut highest_frequency = 0;

    let max_x = *histogram.keys().max().unwrap_or(&0);
    let start = noise_floor.max(noise_floor);
    //iterate to find local maxima
    for x in start..max_x {
        let prev = *histogram.get(&(x - 1)).unwrap_or(&0);
        let curr = *histogram.get(&x).unwrap_or(&0);
        let next = *histogram.get(&(x + 1)).unwrap_or(&0);
        //check if it's a peak (local maximum)
        if curr > prev && curr > next {
            if curr > highest_frequency {
                highest_frequency = curr;
                best_peak_x = Some(x);
            }
        }
    }
    info!("noise floor: {}", noise_floor);
    info!("best peak: {:?}", best_peak_x);
    best_peak_x
}

pub fn inverse_berger_parker_d(histogram: &HashMap<u32, u32>, min_abundance: u32) -> f64 {
    //filter noise to get the Signal bins
    let valid_bins: Vec<(u32, u32)> = histogram
        .iter()
        .filter(|&(&a, &f)| a > min_abundance && f > min_abundance)
        .map(|(&a, &f)| (a, f))
        .collect();
    //calculate total sum of all abundance * frequency
    let total: f64 = valid_bins
        .iter()
        .map(|&(a, f)| (a as f64) * (f as f64))
        .sum();
    if total <= 0.0 {
        return 0.0;
    }
    //find n_max: The area under the curve of the single most dominant peak
    let n_max = valid_bins
        .iter()
        .map(|&(a, f)| (a as f64) * (f as f64))
        .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
        .unwrap_or(0.0);
    let d = n_max / total;
    info!("Total sum in signal: {}", total);
    info!("N_max area sum: {}", n_max);
    info!("Dominance (d): {:.6}", d);
    if d == 0.0 { 0.0 } else { 1.0 / d }
}

//previous berger_parker_d calculation is not strictly berger_parker dominance
pub fn berger_parker_orig(histogram: &HashMap<u32, u32>) -> f64 {
    let total: f64 = histogram
        .iter()
        .map(|(&abundance, &n_k)| (abundance as f64) * (n_k as f64))
        .sum();
    info!("berger-parker histogram sum total is {:?}", &total);
    let n_max = histogram
        .iter()
        .filter(|&(_, &count)| count > 0)
        .map(|(&abundance, _)| abundance)
        .max()
        .unwrap_or(0) as f64;
    info!("berger-parker n_max is {:?}", &n_max);
    if total > 0.0 { n_max / total } else { 0.0 }
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
    info!("histogram read OK");
    Ok(histogram)
}

/// Diagnose why min_abundance sensitivity is high
pub fn abundance_sensitivity_report(histogram: &HashMap<u32, u32>, min_abundance: Option<u32>) {
    let valley_x = find_noise_valley(histogram, 5);
    let bio_peak = find_highest_biological_peak(histogram, valley_x).unwrap_or(0);

    // Check valley depth relative to peaks — shallow valley = unstable cutoff
    let noise_peak =
        find_highest_biological_peak(histogram, min_abundance.unwrap_or(1)).unwrap_or(0);
    let noise_height = *histogram.get(&(noise_peak as u32)).unwrap_or(&0) as f64;
    let valley_height = *histogram.get(&valley_x).unwrap_or(&0) as f64;
    let bio_height = *histogram.get(&bio_peak).unwrap_or(&0) as f64;

    let valley_depth_ratio = valley_height / noise_height.max(bio_height);
    let total_mass: f64 = histogram
        .iter()
        .map(|(&a, &f)| (a as f64) * (f as f64))
        .sum();
    if total_mass > 0.0 {
        // Look at the area ±5 abundance units around the valley
        let start = valley_x.saturating_sub(5);
        let end = valley_x + 5;

        let ambiguous_mass: f64 = (start..=end)
            .map(|x| {
                let f = *histogram.get(&x).unwrap_or(&0) as f64;
                (x as f64) * f
            })
            .sum();

        info!(
            "Mass in ±5 of valley (x={}-{}): {:.1}% of total",
            start,
            end,
            100.0 * ambiguous_mass / total_mass
        );
    }
    info!(
        "Noise peak:          x={}, freq={}",
        noise_peak, noise_height
    );
    info!("Bio peak:            x={}, freq={}", bio_peak, bio_height);
    info!(
        "Valley depth ratio:  {:.3}  (0=deep clean valley, 1=flat/no valley)",
        valley_depth_ratio
    );
    //warn if valley is shallow — means cutoff is unreliable
    if valley_depth_ratio > 0.3 {
        println!(
            "WARNING: kmer histogram noise valley depth ratio {:.3} — min_abundance cutoff is ambiguous",
            valley_depth_ratio
        );
        println!(
            "consider adjusting min_abundance cutoff to a higher value or inspecting the histogram"
        );
        info!(
            "WARNING: kmer histogram noise valley depth ratio {:.3} — min_abundance cutoff is ambiguous",
            valley_depth_ratio
        );
        info!(
            "consider adjusting min_abundance cutoff to a higher value or inspecting the histogram"
        );
    }
}

fn main() -> Result<(), HistogramError> {
    //example histogram: k-mer count -> number of k-mers with that count
    let args = Arguments::parse();
    //initialize the logger
    if let Err(e) = setup_logging(&args.log_file) {
        eprintln!("Failed to initialize logging: {}", e); //continue without logging
    }
    info!("Starting metakora analysis for file: {}", args.filename);
    let histogram = read_histogram_from_file(&args.filename)?;
    let min_abundance = args.min_abundance;
    let (Shannon, H_max, Pielou) = compute_entropy_and_evenness(&histogram);
    let c = chao1(&histogram);
    let o = observed_features(&histogram);
    let r = robbins(&histogram);
    let simp = simpson_index(&histogram);
    let peak_area = peak_area_dominance_integrated(&histogram, min_abundance);
    let bpd = calculate_diversity_auto(&histogram, min_abundance);
    //let bpd = inverse_berger_parker_d(&histogram, *min_abund);
    println!(
        "Sample\tShannon\tH_max\tPielou\tChao1\tObserved\tRobbins\tInv_Berger_Parker\tSimpson\tPeak_Area"
    );
    println!(
        "{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
        &args.filename, Shannon, H_max, Pielou, c, o, r, bpd, simp, peak_area
    );
    let _ = abundance_sensitivity_report(&histogram, min_abundance);
    Ok(())
}
