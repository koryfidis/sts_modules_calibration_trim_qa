# sts_modules_calibration_trim_qa
GSI ISSP 2025 Research Project: CBM STS Detector Electronics QA Pipeline. C++/ROOT script for analyzing and validating ADC trim data from SMX-v2.2 ASICs (1024 ch/side) on silicon microstrip detectors. The pipeline generates QA plots including Mean/SEM statistics, Odd/Even channel analysis, Skewness heatmaps, and Correlation matrices. This script has been used to assess the data from over 30 STS modules (one silicon module contains 8 ASICs per side, 16 in total).

# CBM STS Detector Electronics QA Pipeline (C++/ROOT)

**A quality assurance and calibration validation pipeline developed during the GSI Helmholtz Centre International Summer Student Programme (ISSP) for the CBM (Compressed Baryonic Matter) experiment's Silicon Tracking System (STS).**

### Project Goal

The primary objective of this script is to perform an in-depth **Quality Assurance (QA)** analysis and validate the **ADC trim data** acquired from double-sided silicon microstrip detector modules. It processes data from the **SMX-v2.2 readout ASICs** (8 ASICs, 1024 channels per detector side).

### Key Analysis Features

The C++ script utilizes the **ROOT framework** to execute a multi-step analysis and generate standardized QA plots:

* **Global Performance Metrics**: Calculates and plots the **Mean and Standard Deviation (SEM)** of trim ADC values across all 1024 channels for both P-side (holes) and N-side (electrons), visualizing performance stability.
* **Systematic Effects Check**: Discriminator index analysis plots the mean ADC value as a function of the 31 discriminator settings, including a split by **Odd/Even channels** to identify potential readout crosstalk or systematic noise.
* **Data Integrity & Robustness**:
    * Generates **Heatmaps** of the Mean and **Skewness** (a measure of noise symmetry) across the ASIC and Discriminator dimensions for high-level module QA.
    * Tracks the count of **missing/invalid data points** (values outside the accepted trim range) per channel and discriminator.
* **Correlation Studies**: Computes and visualizes **correlation matrices** for Channel-to-Channel, ASIC-to-ASIC, and Channel-Group-to-Group performance, crucial for diagnosing electronic coupling issues.
* **Channel Grouping Analysis**: Plots the average trim values for **4-channel groups** ($\text{ch } \pmod 4$) to diagnose common-mode noise effects inherent to the ASIC readout structure.
