# Monumentome Publication

This GitHub repository contains the Bash scripts used to process raw sequencing data and R Markdown (`.Rmd`) files used to generate publication figures from MetaPhlAn 4 results.

## 1. Bash Scripts: Sequence Processing and Taxonomic Profiling

The Bash scripts cover the preprocessing and microbial taxonomic profiling workflow:

1. **Adapter trimming:** Remove residual sequencing adapters using **Trim Galore**.
2. **Human read removal:** Align reads against the human reference genome (hg38) using **Bowtie2** and remove reads mapping to the human genome.
3. **Microbial profiling:** Process the cleaned, non-human reads using **MetaPhlAn 4** to identify and estimate the relative abundances of microbial taxa.

## 2. Publication Figure Generation

The repository includes two R Markdown (`.Rmd`) files:

- **Main Figures:** Generates figures included in the main manuscript.
- **Supplementary Figures:** Generates figures included in the supplementary materials.

The analysis uses the following inputs:

1. **MetaPhlAn 4 results:** Species-level abundance data generated from the sequencing analysis.
2. **Sample metadata:** Associated sample information used for downstream analyses and visualizations.
3. **Microbe Directory annotations:** Microbial characteristics and annotations available from [The Microbe Directory](https://docs.superhuman.com/@themicrobedirectory/home/the-microbe-directory-1).
4. **Literature review data:** A curated list of microorganisms previously reported in association with monuments, including organisms implicated in biodeterioration and degradation.

## 3. Running the Analysis

To reproduce the publication figures:

1. Download the required **Microbe Directory** dataset from the link above. If the dataset is unavailable, comment out the corresponding sections of the R Markdown files.
2. Ensure the MetaPhlAn 4 results, sample metadata, and literature review datasets are available at the paths specified in the scripts.
3. Open and execute the `.Rmd` files in **RStudio**, or use AI-assisted coding tools such as **Codex** or **Cursor** to help run and troubleshoot the analysis.

R Markdown files combine R code, explanatory text, and analysis steps in a single document, making the figure-generation workflow easier to follow, modify, and reproduce.

**Note:** Some figures may depend on external annotation files or specific R packages. Ensure the required dependencies and input files are available before running the scripts.
