# Resources

## Purpose of this directory

This directory contains all external reference data and input resources required to run MMCAW, including sequencing inputs, reference genomes, and taxonomic databases. Correct placement and formatting of these resources is essential for reproducible execution of the workflow.

## Sequencing inputs

### FASTQ files
- **Description:** Basecalled and demultiplexed ONT reads generated using Guppy v6.0.6 (high-accuracy model).
- **Format:** `.fastq.gz`
- **Expected location:** User-defined; specified in the metadata file referenced in `config/config.yaml`.

### Guppy sequence summaries
- **Description:** Guppy summary files corresponding to each sequencing run.
- **Format:** `.txt`
- **Usage:** Used for run-level QC and tracking.

### Unblocked read ID lists
- **Description:** Lists of read IDs used to exclude blocked reads from downstream analysis.
- **Format:** Plain text files.

## Reference data

### Human reference genome (GRCh38.p14)
- **Description:** Used for optional host read filtering.
- **Default path in config:** resources/databases/human_reference/GCF_000001405.40_GRCh38.p14_genomic.fna

## Taxonomic databases

### Kraken2 database
- **Description:** Kraken2 standard database including bacterial, archaeal, viral, and human genomes.
- **Default path in config:** ~/Kraken2_Simple_Workflow/resources/databases/krakenstd_06_2023/kraken2_std_database

### BLAST (NT) database
- **Description:** NCBI nucleotide (NT) database used for BLAST-based taxonomic assignment.
- **Default path in config:** resources/databases/NCBI_blast_database/nt

### CAT databases
- **Description:** NCBI taxonomy and protein databases used by CAT for contig classification.
- **Default paths in config:** resources/databases/20240422_CAT_nr/db, resources/databases/20240422_CAT_nr/tax

### NCBI Taxonomy (taxdump)
- **Description:** Required for taxonomy name resolution and LCA assignment.
- **Default path in config:** resources/databases/taxdump

## Optional resources

Additional databases may be generated using the included database creation subworkflow by setting: include_db_creation: True in `config/config.yaml`.

## How to obtain and prepare resources

- Human reference genome: Download from NCBI (GRCh38.p14).
- Kraken2 database: Build using Kraken2 standard database scripts or provide an existing installation.
- BLAST NT database: Download via NCBI BLAST+ `update_blastdb.pl`.
- CAT databases: Download from the official CAT repository or build using provided scripts.

Ensure all paths in `config/config.yaml` correctly point to the corresponding files.

## Directory structure (recommended)

resources/

└── databases/

├── human_reference/

├── krakenstd_06_2023/

├── 20240422_CAT_nr/

├── NCBI_blast_database/

└── taxdump/
