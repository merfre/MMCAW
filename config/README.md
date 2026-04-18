# Configuration

This document describes the structure and parameters of `config/config.yaml` for MMCAW. Users should modify this file to control workflow behaviour, input locations, and analysis settings.

## Overview of config.yaml

The configuration file defines:
- Software environments
- Input metadata location
- Core workflow options (which analysis components to run)
- Database locations
- Tool-specific parameters

## Required metadata

### Sample/run metadata table

The metadata file is specified by: `metadata_file: "config/ONT_mock_cont_assem.txt"`

This file must include, at minimum:
- Sample identifier
- Run identifier
- Paths to FASTQ files
- Paths to sequence summary files
- Paths to unblocked read ID lists (if used)

Column names and exact format should be consistent with the workflow’s expectations in the Snakefile.

## Analysis options (workflow toggles)

These boolean parameters control which components of the workflow are executed:

| Parameter | Description |
|-----------|-------------|
| `include_db_creation` | Build databases within the workflow |
| `include_cat` | Enable CAT classification |
| `include_kraken2` | Enable Kraken2 classification |
| `include_blast` | Enable BLAST classification |
| `include_sourmash` | Enable Sourmash (optional) |
| `include_phylophlan` | Enable PhyloPhlAn phylogenetic analysis |
| `include_comparison` | Compare and merge assigner results (requires ≥2 assigners) |
| `include_rgi` | Identify resistance genes via CARD (requires CAT) |

## Database locations

Users must either place databases in the default locations below or update the paths accordingly:

```yaml
filtering_reference: "resources/databases/human_reference/GCF_000001405.40_GRCh38.p14_genomic.fna"
kraken_db: "~/Kraken2_Simple_Workflow/resources/databases/krakenstd_06_2023/kraken2_std_database"
cat_db: "resources/databases/20240422_CAT_nr/db"
cat_taxonomy: "resources/databases/20240422_CAT_nr/tax"
blast_db: "resources/databases/NCBI_blast_database/nt"
taxdump: "resources/databases/taxdump"
```

## Preprocessing parameters (fastp)

| Parameter                   | Description                                  |
| --------------------------- | -------------------------------------------- |
| `qualified_quality_phred`   | Minimum Phred score to count as qualified    |
| `unqualified_percent_limit` | Maximum % of unqualified bases allowed       |
| `average_qual`              | Minimum average read quality (0 = no filter) |
| `min_length`                | Minimum read length                          |
| `front_trim`                | Bases trimmed from 5’ end                    |
| `tail_trim`                 | Bases trimmed from 3’ end                    |

## Assembly parameters (Flye)

| Parameter         | Description                                         |
| ----------------- | --------------------------------------------------- |
| `read_type`       | `--nano-raw` (default) for uncorrected ONT reads    |
| `minimum_overlap` | Minimum overlap length between reads (default 1000) |

## Taxonomic assignment parameters

### Kraken2
| Parameter           | Description                                     |
| ------------------- | ----------------------------------------------- |
| `kraken_confidence` | Confidence threshold (0–1) for taxonomic labels |

### BLAST
| Parameter               | Description                        |
| ----------------------- | ---------------------------------- |
| `BLAST_min_perc_ident`  | Minimum percent identity           |
| `BLAST_min_evalue`      | Maximum e-value                    |
| `BLAST_max_target_seqs` | Maximum number of target sequences |

### Modified LCA (MLCA) parameters
| Parameter       | Description                |
| --------------- | -------------------------- |
| `MLCA_bitscore` | Minimum bitscore           |
| `MLCA_identity` | Minimum percent identity   |
| `MCLA_coverage` | Minimum alignment coverage |
| `MLCA_majority` | Majority threshold (%)     |
| `MLCA_hits`     | Minimum number of hits     |

### Assigner comparison and consensus

When `include_comparison: True`, MMCAW:

* Standardizes taxonomy outputs from CAT, Kraken2, and BLAST
* Assigns consensus taxonomy if at least two tools agree
* Labels contigs as no_agreement if all three disagree
* Reports percent agreement across taxonomy levels

## Plotting and reporting

| Parameter    | Description                                                        |
| ------------ | ------------------------------------------------------------------ |
| `prevalence` | Number of most prevalent species shown in final plots (default 25) |

## Benchmarking and resources

* `threads`: Number of cores available to Snakemake rules (default 10)
* Benchmarking is enabled by default and outputs rule-level performance metrics.

## Minimal example config (snippet)

```yaml
conda_envs: "workflow/envs/environment.yml"
metadata_file: "config/ONT_mock_cont_assem.txt"
threads: "10"

include_cat: True
include_kraken2: True
include_blast: True
include_comparison: True

filtering_reference: "resources/databases/human_reference/GCF_000001405.40_GRCh38.p14_genomic.fna"
```
