### Rules to perform BLAST on fasta files and convert the results ###

configfile: "config/config.yaml"

### Run blast

rule blastn:
  conda:
    "../envs/environment.yaml"
  input:
    query = "results/preprocessing/flye_results/{PATHS}/assembly.fasta",
  output:
    blast = "results/blast/blast_results/{PATHS}_blast.tsv"
  params:
    outformat = "'6 qseqid stitle sacc staxids pident qcovs evalue bitscore'",
    blast_db = config['blast_db'],
    BLAST_min_perc_ident = config['BLAST_min_perc_ident'],
    BLAST_min_evalue = config['BLAST_min_evalue'],
    BLAST_max_target_seqs = config['BLAST_max_target_seqs'],
    threads = config['threads']
  shell:
    """
    blastn \
    -query {input.query} \
    -db {params.blast_db} \
    -num_threads {params.threads} \
    -outfmt {params.outformat} \
    -perc_identity {params.BLAST_min_perc_ident} \
    -evalue {params.BLAST_min_evalue} \
    -max_target_seqs {params.BLAST_max_target_seqs} \
    -out {output.blast}
    """

### Add full taxonomy to blast results

rule taxonomy_to_blast:
  conda:
    "../envs/environment.yaml"
  input:
    blast = "results/blast/blast_results/{PATHS}_blast.tsv",
    rankedlineage = "resources/databases/taxdump/rankedlineage.dmp"
  output:
    blast_tax = "results/blast/blast_tax/{PATHS}_blast_tax.tsv"
  params:
    taxdump = config['taxdump']
  script:
    "scripts/taxonomy_to_blast.py"

### Assign lowest common ancestor to blast results

rule mlca:
  conda:
    "../envs/environment.yaml"
  input:
    blast = "results/blast/blast_tax/{PATHS}_blast_tax.tsv"
  output:
    lca = "results/blast/mlca/{PATHS}_lca.tsv"
  params:
    bitscore = config['MLCA_bitscore'],
    identity = config['MLCA_identity'],
    coverage = config['MCLA_coverage'],
    majority = config['MLCA_majority'],
    min_hits = config['MLCA_hits']
  script:
    "scripts/mlca.py"

### Count reads for each taxonomy

rule mlca_read_count:
  conda:
    "../envs/environment.yaml"
  input:
    "results/blast/mlca/{PATHS}_lca.tsv"
  params:
    path = "{PATHS}"
  output:
    "results/blast/mlca_read_counts/{PATHS}_lca_counts.tsv"
  script:
    "scripts/mlca_read_count.py"

### Create list of taxonomies identified by blast in all SAMPLES

rule create_blast_lists:
# use unix to create lists for R script input
  conda:
    "../envs/environment.yaml"
  input:
    expand("results/blast/mlca_read_counts/{path}_lca_counts.tsv", path=PATHS)
  output:
    file_path_list = "results/blast/file_path_list.tsv",
    taxonomy_list = "results/blast/taxonomy_list.tsv"
  shell:
    """
    ls {input} > {output.file_path_list};
    cut -f 2-8 {input} | sort | uniq | tee {output.taxonomy_list};
    """
# results/blast/mlca_read_counts/lib/sample_lca_counts.tsv

### Merge BLAST results for all SAMPLES

rule combine_blast_results:
  conda:
    "../envs/environment.yaml"
  input:
    file_path_list = "results/blast/file_path_list.tsv",
    taxonomy_list = "results/blast/taxonomy_list.tsv"
  output:
    "results/blast/blast_merged_results.tsv"
  script:
    "scripts/merging_results.R"

### Create barplots of taxonomy prevalence at various levels

rule taxonomy_summary_barplots_blast:
  conda:
    "../envs/environment.yaml"
  input:
    "results/blast/blast_merged_results.tsv"
  output:
    blast_species = report("results/blast/taxonomy_plots/blast_species_barplot.pdf", caption="report/blast_species_barplot.rst", category="BLAST"),
    blast_family = report("results/blast/taxonomy_plots/blast_family_barplot.pdf", caption="report/blast_family_barplot.rst", category="BLAST"),
    blast_class = report("results/blast/taxonomy_plots/blast_class_barplot.pdf", caption="report/blast_class_barplot.rst", category="BLAST"),
    blast_kingdom = report("results/blast/taxonomy_plots/blast_kingdom_barplot.pdf", caption="report/blast_kingdom_barplot.rst", category="BLAST")
  script:
    "scripts/taxonomy_barplots_blast.py"

### Separate the taxonomy assigner report by taxonomy level(s) for analysis

rule blast_tax_levels:
# separating blast reports by taxonomy levels using R
  conda:
    "../envs/environment.yaml"
  input:
    "results/blast/blast_merged_results.tsv"
  output:
    species = "results/blast/blast_species.tsv",
  script:
    "scripts/separating_tax_levels_blast.R"

### Create species heatmaps of blast results

rule species_heatmap_blast:
  conda:
    "../envs/environment.yaml"
  input:
    "results/blast/blast_species.tsv"
  params:
    prevalence = config['prevalence']
  output:
    report("results/blast/taxonomy_plots/blast_species_heatmap.pdf", caption="report/blast_species_heatmap.rst", category="BLAST")
  script:
    "scripts/taxonomy_heatmaps_blast.R"

### Create stacked taxonomy barplot of assigner outputs

rule taxonomy_plots_blast:
  conda:
    "../envs/environment.yaml"
  input:
    "results/blast/blast_species.tsv"
  params:
    prevalence = config['prevalence']
  output:
    report("results/blast/taxonomy_plots/blast_taxonomy_plot.pdf", caption="report/blast_taxonomy_plot.rst", category="BLAST")
  script:
    "scripts/stacked_taxonomy_barplot_blast.R"

### List BLAST analysis outputs for target rule

rule blast_plots:
  conda:
    "../envs/environment.yaml"
  input:
    tax_plot = "results/blast/taxonomy_plots/blast_taxonomy_plot.pdf",
    barplot = "results/blast/taxonomy_plots/blast_species_barplot.pdf",
    heatmap = "results/blast/taxonomy_plots/blast_species_heatmap.pdf"
  output:
    "results/blast/blast_plot_list.txt"
  shell:
    "ls {input.tax_plot} {input.barplot} {input.heatmap} > {output}"
