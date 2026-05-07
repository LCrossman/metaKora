## metaKora ## 

meta*Kora* is a metagenomics diversity program aimed at non-biased diversity counts including all of the data present.  
Many metagenomics programs call taxa first and then provide diversity statistics on the known component.
Here we provide a method to calculate diversity metrics on the whole dataset.

This is a blazing fast Rust tool for calculating alpha diversity metrics directly from **frequency-of-frequencies (histogram)** data. 

This crate is specifically optimized for large-scale metagenomic data (like k-mer counts) where traditional observation vectors are memory-intensive. Here we take a histogram approach ($k$ abundance $\to$ $n_k$ features), so that even datasets with billions of individuals can be processed quickly and in constant memory relative to the number of unique abundance classes.

You will first need to create a kmer count frequency histogram for each file with your prefered kmer count software, e.g. KMC3, Jellyfish
You can produce histograms from sample metagenomic assemblies or from fastq reads

## Installation 

System Requirements
OS: Linux, macOS, or Windows (via WSL2).

Rust: Version 1.70 or higher (required for compilation), (can be installed with conda/mamba), if you prefer NOT to install Rust, metaKora could be called directly from Python, please leave a GitHub issue to request this feature

Memory: Efficient O(1) memory usage relative to total individuals; performs well even on standard laptops with billion-count datasets.

## 1. Using Conda/Mamba
This is the best method for ensuring your environment is isolated and reproducible across different machines. 
Use the included environment.yml 

```
conda env create -n metakora_env -f environment.yml
conda activate metakora_env

cargo install --git [https://github.com/LCrossman/metaKora.git](https://github.com/LCrossman/metaKora.git) --root $CONDA_PREFIX
```

### 2. Build from source 

If you need to install Rust:
```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```
# Clone the repository
```
git clone https://github.com/LCrossman/metaKora.git
cd metakora
```
# Build the optimized release binary using Rust Cargo 
```
cargo build --release
```
# The binary will be located at:
```
./target/release/metakora --filename <PATH_TO_FILE>
```
### Quick Start 

```
./target/release/metakora --help
```
Tests will be added shortly

## Features

* **Memory Efficient:** Processes data as a "count of counts."
* **Validated Metrics:** Implements Shannon Entropy, Chao1, Pielou's Evenness, Robbins Estimator, Simpson Index, and Berger-Parker Dominance in a consistent manner according to the same metrics as R abdiv and vegan libraries.
* **Robust:** Detailed error reporting with line-number tracking for malformed input files.
* **Text Output Format** Output file is provided to stdout as a tab delimited text file with a single header line

## Diversity Metrics Supported

| Metric | Description 
| :--- | :--- 
| **Shannon Index ($H$)** | Measures uncertainty/diversity using natural logs. 
| **Pielou’s Evenness ($J$)** | Measures how close the community is to numerical equality. 
| **Observed** | Count of total number of distinct kmers.
| **Chao1** | Predicts total estimated richness/complexity of the sample considering total number of unique kmers. 
| **Robbins Estimator** | Probability that the next sample represents a new feature. A 0 indicates that every kmer has already been seen.
| **Inv Berger-Parker** | The reciprocal of the measure of dominance by the most abundant features, calculated after noise removal. 
| **Simpson Index** | Probability that two individuals belong to different features ($1 - D$). 

## Handling Noise Floor (min_abundance)
By default, the script uses a peak-detection algorithm to find the "valley" after the initial noise spike. However, you can override this:

Automatic (Default): The script finds the noise valley and biological peak automatically, logging the details .

Manual: Use --min-abundance <INT> to skip a specific number of abundance classes.  Check the log for any warnings regarding if noise is still detected at the specified threshold.

```bash
./target/release/metakora --filename data.txt --min-abundance 10
```

## Diagnostics & Log Files
Every run generates a log file (default: metakora.log). This file contains a Sensitivity Report which is vital for quality control. It includes:

Noise vs. Bio Peak frequencies.

Valley Depth Ratio: A value > 0.3 warns you if your noise is still ambiguous.

Ambiguous Mass: The % of total k-mer mass sitting near the cutoff point.

```bash
./target/release/metakora --filename data.txt --log-file sample_A.log
```

## Input Format

meta*Kora* expects a two-column, tab-separated file representing the kmer frequency histogram (no headers):
```text
1                  14502    # 14,502 kmers seen once (singletons)
2                  3200     # 3,200 kmers seen twice (doubletons)
3                  1643     # 1,643 kmers seen three times
4                  786      # 786 kmers seen 4 times...
```
**This type of file is outputted by kmer counting programs such as kmc3 and jellyfish.  You may want to consider filtering rare reads at this stage, or rarefaction by subsampling all the samples to the size of the smallest read file.**

You can count kmers from each separate sample metagenomic assembly and run meta*Kora* on each file combining the results in a table for visualization in R 

Please add any issues or requests for features in the GitHub issues.  One feature will be the ability to call this library directly from Python

The output file is a print to stdout of a text file with a single line header: 
```
Sample	Shannon	H_max	Pielou	Chao1	Observed	Robbins	Inv_Berger_Parker	Simpson	Peak_Area
test_input.txt	9.783776	9.910016	0.987261	52991.625625	20131	0.500500	1.000000	0.660022	0.891493
```
