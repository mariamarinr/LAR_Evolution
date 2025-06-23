
<!-- PROJECT SHIELDS -->
[![Python Version](https://img.shields.io/badge/Python-3.8%2B-blue)](https://www.python.org/)
[![R Version](https://img.shields.io/badge/R-4.0%2B-blue)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-MIT-green)](LICENSE)

# LAR_Evolution

<p align="center">
  <strong><font size="4"><h1>Evolutionary fate of the proanthocyanidin biosynthesis gene <em>LAR </em></h1></font></strong>
</p>

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="#overview">Overview</a></li>
    <li><a href="#repository-structure">Repository Structure</a></li>
    <li><a href="#scripts">Scripts</a></li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#requirements">Requirements</a></li>
    <li><a href="#contact">Contact</a></li>
  </ol>
</details>

## Overview
This repository contains scripts used in the study *Evolutionary fate of the proanthocyanidin biosynthesis gene LAR*. These scripts were utilized for data processing, analysis, and visualization.

## Repository Structure

- `data_preprocessing/` - Scripts for cleaning and preparing raw data.
- `expression_analysis/` - Scripts for analyzing gene expression data.
- `README.md` - This file.

## Scripts

### Data Preprocessing
- `clean_fasta.py` - Filters raw sequence data for quality control.
- `extract_kipes_cds.py` - Extract high-quality LAR/DFR/ANR protein sequences from KIPES output and map them to corresponding CDS sequences using a CSV table
  
### Expression Analysis
- `prep_expression_matrix.R` - Filters RNA-seq count data based on a list of candidate genes, transposes the expression matrix, and merges it with sample metadata to prepare a clean dataset for downstream analyses (e.g., visualization, statistical modeling).
- `LMM_expression_analysis.R` - Conducts linear mixed-effects modeling (LMM) analysis on RNA-seq data, performs post-hoc comparisons, and visualizes gene expression patterns with boxplots and violin plots.

## Usage
Each script contains comments and instructions for execution. To run a script, ensure you have the required dependencies installed. Example:

```bash
python3 clean_fasta.py input.fasta output.fasta
```

## Requirements
- Python (>=3.8)
- R (>=4.0)
- Required Python libraries: Biopython, NumPy, Pandas
- Required R packages: ggplot2, lme4, lmerTest, dplyr, tidyr, emmeans 

## Contact
For any questions or issues, please open an issue in this repository or contact Maria F Marin-Recinos at mafer.recinos92@gmail.com.

---

