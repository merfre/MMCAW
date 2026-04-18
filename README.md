# MMCAW — Metagenomic Microbiome Combination Analysis Workflow

## Overview

The Metagenomic Microbiome Combination Analysis Workflow (MMCAW) is a reproducible Snakemake-based pipeline for the analysis of Oxford Nanopore Technologies (ONT) metagenomic sequencing data. The workflow was developed for the analysis of skin and wound microbiome datasets as part of a PhD thesis at the University of Hull, but is applicable to other metagenomic nanopore datasets with appropriate configuration.

MMCAW integrates multiple taxonomic assigners, performs comparative and consensus-based classification, and generates community-level diversity analyses while tracking computational performance and resource usage.

## Scope of the workflow

MMCAW supports the following major analytical components:

- **Read preprocessing and quality control**
- **Optional host read removal**
- **Optional metagenomic assembly**
- **Taxonomic assignment using multiple tools:**
  - CAT  
  - Kraken2  
  - BLAST (with modified LCA implementation)  
- **Comparison and consensus across assigners**
- **Community composition and diversity analyses**
- **Built-in benchmarking and performance tracking**

All steps are implemented as modular Snakemake rules, allowing users to enable or disable components via the configuration file.

## Workflow structure (DAG)

The directed acyclic graph (DAG) below represents the structure of the workflow, including major processing steps and dependencies between rules.

<img width="1851" height="1499" alt="MMCAW_dag" src="https://github.com/user-attachments/assets/82b9ac4c-122a-4500-ba86-a7345aefe2b5" />

## Requirements

### Software
- Snakemake **v7.22.0**
- Conda (Python **v3.10.8**)

### Environment management
Software dependencies are managed via Conda using the environment specified in: `workflow/envs/environment.yml`

## Installation

Clone the repository:
```bash
git clone <repository-url>
cd MMCAW
```

Create and activate the Conda environment:
```bash
conda env create -f workflow/envs/environment.yml
conda activate mmcaw
```

## Usage
### Basic run (local)

Perform a dry run to check the workflow:
```bash
snakemake -n
```
Run the workflow using multiple cores:
```bash
snakemake --cores 10
```
Execute the workflow locally:
```bash
snakemake --printshellcmds --use-conda --cores 10
```
After successful execution, you can create a self-contained interactive HTML report with all results:
```bash
snakemake --report hmcw_final_report.html
```

### Running on an HPC / cluster

MMCAW was developed and tested on the University of Hull’s Viper HPC. If using a cluster, configure and run with an appropriate Snakemake profile (e.g., SLURM, PBS, etc.):
```bash
snakemake --profile <your-cluster-profile>
```
### Inputs (brief)
MMCAW expects the following inputs:

* Basecalled and demultiplexed FASTQ files (Guppy v6.0.6, high-accuracy model)
* Guppy sequence summary files
* Unblocked read ID lists
* A sample/run metadata file specified in config/config.yaml

Detailed input requirements are described in `resources/README.md`.

### Outputs (brief)

Key outputs include:

* Quality-filtered reads (fastp)
* Host-filtered reads (optional)
* Metagenomic assemblies (Flye, optional)
* Taxonomic classifications from CAT, Kraken2, and BLAST
* Consensus taxonomy tables across assigners
* Community composition and diversity analyses (alpha/beta diversity)
* Snakemake benchmarking reports (runtime and resource usage per rule)

### Repository structure

.
├── workflow/
│   ├── Snakefile
│   ├── rules/
│   ├── envs/
│   └── config/
├── resources/
│   └── databases/
└── config/
    └── config.yaml

Detailed descriptions of resources and configuration are provided in:

* `resources/README.md`
* `config/README.md`

## Reproducibility & benchmarking
* Workflow implemented in Snakemake v7.22.0
* All software dependencies are managed via Conda
* Snakemake’s built-in benchmarking is enabled by default to record:
  * Rule-level runtime
  * CPU and memory usage
  * Resource performance across datasets of varying size and complexity

This supports systematic evaluation of workflow efficiency and scalability.

### Data availability

Where feasible, raw sequencing data and associated bioinformatic workflows have been archived:

* Zenodo: doi: 10.5281/zenodo.17752049

## Citation / Thesis

If you use MMCAW in your work, please cite:

Merideth Naomi Freiheit (2025). Development of Reproducible Metagenomic Approaches for Skin and Wound Microbiome Analysis. University of Hull.
