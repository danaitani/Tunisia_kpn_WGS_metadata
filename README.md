**Tunisia *Klebsiella pneumoniae* WGS Data and Code Repository**

## Description

This repository contains de-identified metadata and analysis code associated with whole-genome sequencing (WGS) of third-generation cephalosporin-resistant *Klebsiella pneumoniae* isolates collected through the Tunisian Antimicrobial Resistance Surveillance System (TARSS).

These data were used for genomic epidemiological analyses presented in a PhD thesis at the London School of Hygiene & Tropical Medicine and in related research manuscripts.

## Sequence data

Raw sequencing reads are publicly available in the European Nucleotide Archive (ENA):

* **BioProject:** PRJEB98550
* https://www.ebi.ac.uk/ena/browser/view/PRJEB98550

The ENA archive contains raw sequencing reads and accession numbers for all isolates included in the study.

## Metadata availability

This repository provides curated, de-identified metadata used for analysis and reproducibility

To minimise the risk of re-identification, detailed patient-level information (e.g. exact dates, age, ward identifiers) is not publicly shared.

## Code availability

This repository also contains the analysis code used to generate the results and figures presented in the manuscript, including:

* R scripts for data cleaning, statistical analysis, and visualization
* Python scripts for gene detection and MLST calling

The code is provided to support reproducibility of the analyses.

## Requirements

### R packages

* dplyr, tidyr, ggplot2, patchwork, cowplot
* AMRgen, DescTools, igraph, ggraph
* readxl, writexl, stringi

### Python

* Python 3
* pandas, pysam, biopython

### External tools

* bowtie2
* samtools

## Citation

If you use these data or code, please cite the associated manuscript:

Itani D, Smaoui H, Thabet L, Zribi M, Dhraief S, Kanzari L, et al. Genomic Surveillance of Third-Generation Cephalosporin-Resistant Klebsiella pneumoniae in Tunisian AMR Surveillance System Hospitals [Internet]. medRxiv; 2026 [cited 2026 Apr 11]. p. 2026.04.08.26350452. Available from: https://www.medrxiv.org/content/10.64898/2026.04.08.26350452v1 doi:10.64898/2026.04.08.26350452

