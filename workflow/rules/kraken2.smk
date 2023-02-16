### Rule to assign taxonomy using kraken2 ###

configfile: "config/config.yaml"

### Run kraken2 and create reports

rule kraken2:
  #conda:
    #"../environment.yml"
  input:
    "results/Preprocessing/flye_results/{PATHS}/assembly.fasta"
  output:
    report = "results/Kraken2/{PATHS}_kraken_report.txt",
    kraken = "results/Kraken2/{PATHS}_kraken.krk"
    # Needs to be .krk for recentrifuge
  params:
    kraken_db = config['kraken_db'],
    threads = config['threads'],
    confidence = config['kraken_confidence']
  shell:
    "kraken2 --threads {params.threads} --db {params.kraken_db} {input} \
    --confidence {params.confidence} --report {output.report} --output {output.kraken}"
    # use taxonomy names, use threads from config file, path to database is in the config file too

### Add full taxonomy to blast results

rule taxonomy_to_kraken:
  #conda:
    #"../environment.yml"
  input:
    kraken = "results/Kraken2/{PATHS}_kraken.krk",
    rankedlineage = "resources/databases/taxdump/rankedlineage.dmp"
  output:
    kraken_tax = "results/Kraken2/{PATHS}_kraken_tax.tsv"
  params:
    taxdump = config['taxdump']
  script:
    "scripts/taxonomy_to_kraken.py"

### Count reads for each taxonomy

rule kraken_read_count:
  #conda:
    #"../environment.yml"
  input:
    "results/Kraken2/{PATHS}_kraken_tax.tsv"
  params:
    path = "{PATHS}"
  output:
    "results/Kraken2/kraken_read_counts/{PATHS}_kraken_counts.tsv"
  script:
    "scripts/kraken_read_count.py"

### Create list of taxonomies identified by kraken2 in all SAMPLES

rule create_kraken_lists:
# use unix to create lists for R script input
  #conda:
  #    #"../environment.yml"
  input:
    expand("results/Kraken2/kraken_read_counts/{path}_kraken_counts.tsv", path=PATHS)
  output:
    file_path_list = "results/Kraken2/file_path_list.tsv",
    taxonomy_list = "results/Kraken2/taxonomy_list.tsv"
  shell:
    """
    ls {input} > {output.file_path_list};
    cut -f 2-8 {input} | sort | uniq | tee {output.taxonomy_list};
    """

### Merge Kraken2 results for all SAMPLES

rule combine_kraken_results:
  #conda:
  #    #"../environment.yml"
  input:
    file_path_list = "results/Kraken2/file_path_list.tsv",
    taxonomy_list = "results/Kraken2/taxonomy_list.tsv"
  output:
    "results/Kraken2/kraken_merged_results.tsv"
  script:
    "scripts/merging_results.R"

### Create barplots of taxonomy prevalence at various levels

rule taxonomy_summary_barplots_kraken:
  #conda:
   #"../environment.yml"
  input:
    "results/Kraken2/kraken_merged_results.tsv"
  output:
    kraken_species = report("results/Kraken2/taxonomy_plots/kraken2_species_barplot.pdf", caption="report/kraken2_species_barplot.rst", category="Kraken2"),
    kraken_family = report("results/Kraken2/taxonomy_plots/kraken2_family_barplot.pdf", caption="report/kraken2_family_barplot.rst", category="Kraken2"),
    kraken_class = report("results/Kraken2/taxonomy_plots/kraken2_class_barplot.pdf", caption="report/kraken2_class_barplot.rst", category="Kraken2"),
    kraken_kingdom = report("results/Kraken2/taxonomy_plots/kraken2_kingdom_barplot.pdf", caption="report/kraken2_kingdom_barplot.rst", category="Kraken2")
  script:
    "scripts/taxonomy_barplots_kraken.py"

### Separate the taxonomy assigner report by taxonomy level(s) for analysis

rule kraken_tax_levels:
# separating kraken2 reports by taxonomy levels using R
  #conda:
    #"../environment.yml"
  input:
    "results/Kraken2/kraken_merged_results.tsv"
  output:
    species = "results/Kraken2/kraken2_species.tsv",
  script:
    "scripts/separating_tax_levels_kraken.R"

### Create species heatmaps of kraken2 results

rule species_heatmap_kraken:
  #conda:
       #"../environment.yml"
  input:
    "results/Kraken2/kraken2_species.tsv"
  output:
    report("results/Kraken2/taxonomy_plots/kraken2_species_heatmap.pdf", caption="report/kraken2_species_heatmap.rst", category="Kraken2")
  script:
    "scripts/taxonomy_heatmaps_kraken.R"

### Create stacked taxonomy barplot of assigner outputs

rule taxonomy_plots_kraken:
  #conda:
       #"../environment.yml"
  input:
    "results/Kraken2/kraken2_species.tsv"
  output:
    report("results/Kraken2/taxonomy_plots/kraken2_taxonomy_plot.pdf", caption="report/kraken2_taxonomy_plot.rst", category="Kraken2")
  script:
    "scripts/stacked_taxonomy_barplot_kraken.R"

### List Kraken2 analysis outputs for target rule

rule kraken_plots:
  #conda:
    #"../environment.yml"
  input:
    tax_plot = "results/Kraken2/taxonomy_plots/kraken2_taxonomy_plot.pdf",
    barplot = "results/Kraken2/taxonomy_plots/kraken2_species_barplot.pdf",
    heatmap = "results/Kraken2/taxonomy_plots/kraken2_species_heatmap.pdf"
  output:
    "results/Kraken2/kraken2_plot_list.txt"
  shell:
    "ls {input.tax_plot} {input.barplot} {input.heatmap} > {output}"
