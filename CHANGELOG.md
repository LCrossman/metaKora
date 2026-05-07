# Changelog

All notable changes to the **metaKora** project will be documented in this file.

## [0.2.0] - 2026-05-07

### Added
- **Automated Noise Detection:** Implemented a "noise valley" detection algorithm that automatically identifies the transition point between sequencing errors (low-abundance noise) and true biological signal in k-mer histograms.
- **Allowed a manual threshold to be chosen, and a logged comment for a warning if the automatic detection suggests that the threshold has not passed the sequencing noise
- **Inverse Berger-Parker Dominance:** Added a refined dominance metric that calculates the reciprocal of the Berger-Parker index ($1/d$) after filtering for the detected noise floor.
- **Peak Area Metric:** Introduced a new metric that integrates the area under the primary biological peak (15% width) to assess sample dominance more robustly than single-point indices within the kmer data.
- **Logging System:** Integrated `env_logger` and `chrono` to generate persistent log files. This replaces internal println! debug statements with structured logging.
- **Sensitivity Reporting:** A logging suite that calculates:
    - **Valley Depth Ratio:** Quantifies the "cleanliness" of the noise-to-signal separation.
    - **Ambiguous Mass:** Measures the percentage of total library mass residing in the transition zone near the cutoff.
- **CLI Arguments:** Added `--log-file` (default: `metakora.log`) and `--min-abundance` flags to allow for manual override of the automatic detection parameters.

### Changed
- **Code Organization:** Moved error types into a dedicated `HistogramError` enum for better error propagation and more informative terminal output.
- **Output Format:** Updated the final `stdout` results table to include `Inv_Berger_Parker` and `Peak_Area` columns

### Improved
- **Signal Processing:** Optimized moving-average smoothing for the histogram using a fixed window size to prevent micro-fluctuations from triggering false peaks during auto-thresholding.

### Fixed
- Resolved a potential panic in peak detection by using `saturating_sub` for low-abundance calculations
---

## [0.1.0] - 2026-04-16
### Added
- Initial release of metaKora.
- Basic support for Shannon, Pielou, Chao1, observed, berger-parker_d and Robbins metrics.
- Support for whitespace-separated frequency-of-frequency input files.
