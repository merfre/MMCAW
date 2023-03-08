### Rule to assign taxonomy using Kraken2 ###

configfile: "config/config.yaml"

### Run kraken2 and create reports

rule kraken2:
  conda:
    "../envs/environment.yaml"
  input:
    "results/preprocessing/flye_results/{PATHS}/assembly.fasta"
  output:
    report = "results/kraken2/{PATHS}_kraken_report.txt",
    kraken = "results/kraken2/{PATHS}_kraken.krk"
    # Needs to be .krk for recentrifuge
  params:
    kraken_db = config['kraken_db'],
    threads = config['threads'],
    confidence = config['kraken_confidence']
  shell:
    """
    kraken2 --threads {params.threads} --db {params.kraken_db} {input} \
    --confidence {params.confidence} --report {output.report} --output {output.kraken}
    """
    # use taxonomy names, use threads from config file, path to database is in the config file too

### Add full taxonomy to blast results

rule taxonomy_to_kraken:
  conda:
    "../envs/environment.yaml"
  input:
    kraken = "results/kraken2/{PATHS}_kraken.krk",
    rankedlineage = "resources/databases/taxdump/rankedlineage.dmp",
    merged = "resources/databases/taxdump/merged.dmp"
  output:
    kraken_tax = "results/kraken2/{PATHS}_kraken_tax.tsv"
  script:
    "scripts/taxonomy_to_kraken.py"

### Count reads for each taxonomy

rule kraken_read_count:
  conda:
    "../envs/environment.yaml"
  input:
    "results/kraken2/{PATHS}_kraken_tax.tsv"
  params:
    path = "{PATHS}"
  output:
    "results/kraken2/kraken_read_counts/{PATHS}_kraken_counts.tsv"
  script:
    "scripts/kraken_read_count.py"

### Create list of taxonomies identified by kraken2 in all SAMPLES

rule create_kraken_lists:
# use unix to create lists for R script input
  conda:
    "../envs/environment.yaml"
  input:
    expand("results/kraken2/kraken_read_counts/{path}_kraken_counts.tsv", path=PATHS)
  output:
    file_path_list = "results/kraken2/file_path_list.tsv",
    taxonomy_list = "results/kraken2/taxonomy_list.tsv"
  shell:
    """
    ls {input} > {output.file_path_list};
    cut -f 2-8 {input} | sort | uniq | tee {output.taxonomy_list};
    """

### Merge Kraken2 results for all SAMPLES

rule combine_kraken_results:
  conda:
    "../envs/environment.yaml"
  input:
    file_path_list = "results/kraken2/file_path_list.tsv",
    taxonomy_list = "results/kraken2/taxonomy_list.tsv"
  output:
    "results/kraken2/kraken_merged_results.tsv"
  script:
    "scripts/merging_results.R"

### Create barplots of taxonomy prevalence at various levels

rule taxonomy_summary_barplots_kraken:
  conda:
    "../envs/environment.yaml"
  input:
    "results/kraken2/kraken_merged_results.tsv"
  output:
    kraken_species = report("results/kraken2/taxonomy_plots/kraken2_species_barplot.pdf", caption="report/kraken2_species_barplot.rst", category="Kraken2"),
    kraken_family = report("results/kraken2/taxonomy_plots/kraken2_family_barplot.pdf", caption="report/kraken2_family_barplot.rst", category="Kraken2"),
    kraken_class = report("results/kraken2/taxonomy_plots/kraken2_class_barplot.pdf", caption="report/kraken2_class_barplot.rst", category="Kraken2"),
    kraken_kingdom = report("results/kraken2/taxonomy_plots/kraken2_kingdom_barplot.pdf", caption="report/kraken2_kingdom_barplot.rst", category="Kraken2")
  script:
    "scripts/taxonomy_barplots_kraken.py"

### Separate the taxonomy assigner report by taxonomy level(s) for analysis

rule kraken_tax_levels:
# separating kraken2 reports by taxonomy levels using R
  conda:
    "../envs/environment.yaml"
  input:
    "results/kraken2/kraken_merged_results.tsv"
  output:
    species = "results/kraken2/kraken2_species.tsv",
  script:
    "scripts/separating_tax_levels_kraken.R"

### Create species heatmaps of kraken2 results

rule species_heatmap_kraken:
  conda:
    "../envs/environment.yaml"
  input:
    "results/kraken2/kraken2_species.tsv"
  params:
    prevalence = config['prevalence']
  output:
    report("results/kraken2/taxonomy_plots/kraken2_species_heatmap.pdf", caption="report/kraken2_species_heatmap.rst", category="Kraken2")
  script:
    "scripts/taxonomy_heatmaps_kraken.R"

### Create stacked taxonomy barplot of assigner outputs

rule taxonomy_plots_kraken:
  conda:
    "../envs/environment.yaml"
  input:
    "results/kraken2/kraken2_species.tsv"
  params:
    prevalence = config['prevalence']
  output:
    report("results/kraken2/taxonomy_plots/kraken2_taxonomy_plot.pdf", caption="report/kraken2_taxonomy_plot.rst", category="Kraken2")
  script:
    "scripts/stacked_taxonomy_barplot_kraken.R"

### List Kraken2 analysis outputs for target rule

rule kraken_plots:
  conda:
    "../envs/environment.yaml"
  input:
    tax_plot = "results/kraken2/taxonomy_plots/kraken2_taxonomy_plot.pdf",
    barplot = "results/kraken2/taxonomy_plots/kraken2_species_barplot.pdf",
    heatmap = "results/kraken2/taxonomy_plots/kraken2_species_heatmap.pdf"
  output:
    "results/kraken2/kraken2_plot_list.txt"
  shell:
    "ls {input.tax_plot} {input.barplot} {input.heatmap} > {output}"
