Analysis steps and rules for the Metagenomic Microbiome Combination Analysis Workflow


This report contains the results for the samples located in this metadata file "{{ snakemake.config["metadata_file"] }}" and analysis was performed using the conda environment "{{ snakemake.config["conda_envs"] }}"



Analysis options included in this report (indicated by True/False in config):



Database creation subworkflow: {{ snakemake.config["include_db_creation"] }}


Contig annotation tool for taxonomic assignment: {{ snakemake.config["include_cat"] }}


Kraken2 for taxonomic assignment: {{ snakemake.config["include_kraken2"] }}


BLAST for taxonomic assignment: {{ snakemake.config["include_blast"] }}


Sourmash for taxonomic assignment: {{ snakemake.config["include_sourmash"] }}


PhyloPhlAn for phylogenetic analysis: {{ snakemake.config["include_phylophlan"] }}


Comparison of selected taxnomic assigners: {{ snakemake.config["include_comparison"] }}


Resistance gene identification with CARD: {{ snakemake.config["include_rgi"] }}