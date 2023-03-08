### Rules to assign taxonomy and annotate with CAT and produce a CAT report ###

configfile: "config/config.yaml"

### Run CAT

rule cat:
  conda:
    "../envs/environment.yaml"
  input:
    "results/preprocessing/flye_results/{PATHS}/assembly.fasta"
  output:
    next_rule_input = "results/cat/{PATHS}.contig2classification.txt"
    # output necessary to run rule because it is the input for the report rule
  params:
    prefix = "results/cat/{PATHS}",
    # prefix includes the temporary directory that cat is running in: ./results/cat/
    cat_db = config['cat_db'],
    cat_taxonomy = config['cat_taxonomy'],
    threads = config['threads']
  shell:
    "CAT contigs -n {params.threads} --force -c {input} -d {params.cat_db} -t {params.cat_taxonomy} -o {params.prefix}"
    # using CAT to produce taxonomy report using the NCBI NR database from the CAT database
    # which was created prior to running snakemake with: CAT prepare --fresh
    # force option overwrites files if already run
    # -n 10 flag is the number of cores to use when running

### Create CAT clean report

rule cat_report:
# this rule is to add names to the cat report and the input is the output from the contigs run
  conda:
    "../envs/environment.yaml"
  input:
    "results/cat/{PATHS}.contig2classification.txt"
  output:
    cat = "results/cat/{PATHS}_cat_names.txt"
  params:
    cat_taxonomy = config['cat_taxonomy']
  shell:
    """
    CAT add_names -i {input} -o {output.cat} -t {params.cat_taxonomy} --only_official --exclude_scores;
    sed -i 's/# //g' {output.cat};
    sed -i 's/superkingdom/kingdom/g' {output.cat};
    sed -i 's/contig/query/g' {output.cat};
    sed -i 's/no support/unidentified/g' {output.cat};
    sed -i 's/ /_/g' {output.cat};
    """
    # the exclude scores flag will remove the matching scores CAT produces and cleans up the report

### Count reads for each taxonomy

rule cat_read_count:
  conda:
    "../envs/environment.yaml"
  input:
    "results/cat/{PATHS}_cat_names.txt"
  params:
    path = "{PATHS}"
  output:
    "results/cat/cat_read_counts/{PATHS}_cat_counts.tsv"
  script:
    "scripts/cat_read_count.py"

### Create list of taxonomies identified by CAT in all SAMPLES

rule create_cat_lists:
# use unix to create lists for R script input
  conda:
    "../envs/environment.yaml"
  input:
    expand("results/cat/cat_read_counts/{path}_cat_counts.tsv", path=PATHS)
  output:
    file_path_list = "results/cat/file_path_list.tsv",
    taxonomy_list = "results/cat/taxonomy_list.tsv"
  shell:
    """
    ls {input} > {output.file_path_list};
    cut -f 2-8 {input} | sort | uniq | tee {output.taxonomy_list};
    """

### Merge CAT results for all SAMPLES

rule combine_cat_results:
# use R to merge cat reports by clade
  conda:
    "../envs/environment.yaml"
  input:
    file_path_list = "results/cat/file_path_list.tsv",
    taxonomy_list = "results/cat/taxonomy_list.tsv"
  output:
    "results/cat/cat_merged_results.tsv"
  script:
    "scripts/merging_results.R"

### Create barplots of taxonomy prevalence at various levels

rule taxonomy_summary_barplots_cat:
  conda:
    "../envs/environment.yaml"
  input:
    "results/cat/cat_merged_results.tsv"
  output:
    cat_species = report("results/cat/taxonomy_plots/cat_species_barplot.pdf", caption="report/cat_species_barplot.rst", category="CAT"),
    cat_family = report("results/cat/taxonomy_plots/cat_family_barplot.pdf", caption="report/cat_family_barplot.rst", category="CAT"),
    cat_class = report("results/cat/taxonomy_plots/cat_class_barplot.pdf", caption="report/cat_class_barplot.rst", category="CAT"),
    cat_kingdom = report("results/cat/taxonomy_plots/cat_kingdom_barplot.pdf", caption="report/cat_kingdom_barplot.rst", category="CAT")
  script:
    "scripts/taxonomy_barplots_cat.py"

### Separate the taxonomy assigner report by taxonomy level(s) for analysis

rule cat_tax_levels:
# separating CAT reports by taxonomy levels using grep
  conda:
    "../envs/environment.yaml"
  input:
    "results/cat/cat_merged_results.tsv"
  output:
    species = "results/cat/cat_species.tsv"
  script:
    "scripts/separating_tax_levels_cat.R"

### Create species heatmaps of cat results

rule species_heatmap_cat:
  conda:
    "../envs/environment.yaml"
  input:
    "results/cat/cat_species.tsv"
  params:
    prevalence = config['prevalence']
  output:
    report("results/cat/taxonomy_plots/cat_species_heatmap.pdf", caption="report/cat_species_heatmap.rst", category="CAT")
  script:
    "scripts/taxonomy_heatmaps_cat.R"

### Create stacked taxonomy barplot of assigner outputs

rule taxonomy_plots_cat:
  conda:
    "../envs/environment.yaml"
  input:
    "results/cat/cat_species.tsv"
  params:
    prevalence = config['prevalence']
  output:
    report("results/cat/taxonomy_plots/cat_taxonomy_plot.pdf", caption="report/cat_taxonomy_plot.rst", category="CAT")
  script:
    "scripts/stacked_taxonomy_barplot_cat.R"

### List CAT analysis outputs for target rule

rule cat_plots:
  conda:
    "../envs/environment.yaml"
  input:
    tax_plot = "results/cat/taxonomy_plots/cat_taxonomy_plot.pdf",
    barplot = "results/cat/taxonomy_plots/cat_species_barplot.pdf",
    heatmap = "results/cat/taxonomy_plots/cat_species_heatmap.pdf"
  output:
    "results/cat/cat_plot_list.txt"
  shell:
    "ls {input.tax_plot} {input.barplot} {input.heatmap} > {output}"
