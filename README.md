## metaKora ## 

meta*Kora* is a metagenomics diversity program aimed at non-biased diversity counts including all of the data present.  
Many metagenomics programs call taxa first and then provide diversity statistics on the known component.
Here we provide a method to calculate diversity metrics on the whole dataset.

This is a blazing fast Rust tool for calculating alpha diversity metrics directly from **frequency-of-frequencies (histogram)** data. 

This crate is specifically optimized for large-scale metagenomic data (like k-mer counts) where traditional observation vectors are memory-intensive. Here we take a histogram approach ($k$ abundance $\to$ $n_k$ features), so that even datasets with billions of individuals can be processed quickly and in constant memory relative to the number of unique abundance classes.

## Installation

## Features

* **Memory Efficient:** Processes data as a "count of counts."
* **Validated Metrics:** Implements Shannon Entropy, Chao1, Pielou's Evenness, Robbins Estimator, Simpson Index, and Berger-Parker Dominance in a consistent manner according to the same metrics as R abdiv and vegan libraries.
* **Robust:** Detailed error reporting with line-number tracking for malformed input files.
* **Text Output Format** Output file is provided to stdout as a tab delimited text file with a single header line

## Diversity Metrics Supported

| Metric | Function | Description |
| :--- | :--- | :--- |
| **Shannon Index ($H$)** | Measures uncertainty/diversity using natural logs. |
| **Pielou’s Evenness ($J$)** | Measures how close the community is to numerical equality. |
| **Chao1** | Predicts total richness including unobserved species (Bias-Corrected). |
| **Robbins Estimator** | Probability that the next sample represents a new feature. |
| **Berger-Parker** | Measure of dominance by the most abundant feature. |
| **Simpson Index** | Probability that two individuals belong to different species ($1 - D$). |

## Input Format

The tool expects a two-column, space-separated file representing the kmer frequency histogram (no headers):
```text
1                  14502    # 14,502 species seen once (singletons)
2                  3200     # 3,200 species seen twice (doubletons)
15                 1        # 1 species seen fifteen times
```
This type of file is outputted by kmer counting programs such as kmc3 and jellyfish 

You can count kmers from each separate sample metagenomic assembly and run meta*Kora* on each file combining the results in a table for visualization in R 

Please add any issues or requests for features in the GitHub issues.  One feature will be the ability to call this library directly from Python
