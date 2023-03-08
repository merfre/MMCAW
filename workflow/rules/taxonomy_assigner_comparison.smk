### Rules to compare taxonomy assigner results ###

configfile: "config/config.yaml"

### Creating alpha diversity scores and plots for all taxonomic assigners

rule alpha_diversity:
  conda:
    "../envs/environment.yaml"
  input:
    cat = "results/cat/cat_species.tsv",
    kraken = "results/kraken2/kraken2_species.tsv",
    blast = "results/blast/blast_species.tsv"
  output:
    alpha_div_table = report("results/assigner_comparison/alpha_diversity/alpha_div_table.tsv", caption="report/alpha_div_table.rst", category="Taxonomy assigner comparison", subcategory="Alpha diversity"),
    shannon_blast_cat = report("results/assigner_comparison/alpha_diversity/shannon_scatterplot_blast_cat.pdf", caption="report/shannon_scatterplot_blast_cat.rst", category="Taxonomy assigner comparison", subcategory="Alpha diversity"),
    shannon_blast_kraken = report("results/assigner_comparison/alpha_diversity/shannon_scatterplot_blast_kraken.pdf", caption="report/shannon_scatterplot_blast_kraken.rst", category="Taxonomy assigner comparison", subcategory="Alpha diversity"),
    shannon_cat_kraken = report("results/assigner_comparison/alpha_diversity/shannon_scatterplot_cat_kraken.pdf", caption="report/shannon_scatterplot_cat_kraken.rst", category="Taxonomy assigner comparison", subcategory="Alpha diversity"),
    shannon_barplot = report("results/assigner_comparison/alpha_diversity/shannon_barplot.pdf", caption="report/shannon_barplot.rst", category="Taxonomy assigner comparison", subcategory="Alpha diversity")
  script:
    "scripts/alpha_diversity.R"

### Creating beta diversity PCA plot of all taxonomic assigners

rule beta_diversity_pca:
  conda:
    "../envs/environment.yaml"
  input:
    kraken2 = "results/kraken2/kraken2_species.tsv",
    cat = "results/cat/cat_species.tsv",
    blast = "results/blast/blast_species.tsv"
  output:
    report("results/assigner_comparison/beta_diversity/beta_diversity_pca.pdf", caption="report/beta_diversity_pca.rst", category="Taxonomy assigner comparison", subcategory="Beta diversity")
  script:
    "scripts/beta_div_pca.R"

### Creating summary table of taxonomic assigner results to compare

rule taxonomy_assigner_summary:
  conda:
    "../envs/environment.yaml"
  input:
    cat_results = "results/cat/cat_merged_results.tsv",
    kraken_results = "results/kraken2/kraken_merged_results.tsv",
    blast_results = "results/blast/blast_merged_results.tsv"
  output:
    report("results/assigner_comparison/taxonomy_assigner_summary.tsv", caption="report/taxonomy_assigner_summary.rst", category="Taxonomy assigner comparison", subcategory="Comparison summaries")
  params:
    seq_stats = "results/qc_reports/humrm_qc_report.tsv"
  script:
    "scripts/assigner_summary.py"

### Comparing taxonomy assignment for each read of all assigner_summary

rule taxonomy_comparison:
  conda:
    "../envs/environment.yaml"
  input:
    cat_results = "results/cat/{PATHS}_cat_names.txt",
    kraken_results = "results/kraken2/{PATHS}_kraken_tax.tsv",
    blast_results = "results/blast/mlca/{PATHS}_lca.tsv"
  output:
    sample_comparison = report("results/assigner_comparison/{PATHS}_comparison.tsv", caption="report/tax_read_comparison.rst", category="Taxonomy assigner comparison", subcategory="Direct taxonomy comparison"),
    sample_summary = "results/assigner_comparison/{PATHS}_comparison_summary.tsv"
  params:
    path = "{PATHS}"
  script:
    "scripts/read_tax_comparison.py"

### Create list of tax comparison summaries for all SAMPLES

rule create_summary_lists:
# use unix to create lists for R script input
  conda:
    "../envs/environment.yaml"
  input:
    expand("results/assigner_comparison/{path}_comparison_summary.tsv", path=PATHS)
  output:
    file_path_list = "results/assigner_comparison/file_path_list.tsv"
  shell:
    "ls {input} > {output.file_path_list}"

### Summary of taxonomy comparison of all samples

rule combine_tax_comparison:
  conda:
    "../envs/environment.yaml"
  input:
    file_path_list = "results/assigner_comparison/file_path_list.tsv"
  output:
    report("results/assigner_comparison/tax_comparison_summary.tsv", caption="report/tax_comparison_summary.rst", category="Taxonomy assigner comparison", subcategory="Comparison summaries")
  script:
    "scripts/merging_tax_comparison.R"

### List taxonomy assigner comparison outputs for target rule

rule comparison_outputs:
  conda:
    "../envs/environment.yaml"
  input:
    table = "results/assigner_comparison/taxonomy_assigner_summary.tsv",
    alpha_div = "results/assigner_comparison/alpha_diversity/alpha_div_table.tsv",
    beta_div = "results/assigner_comparison/beta_diversity/beta_diversity_pca.pdf",
    tax_compare = "results/assigner_comparison/tax_comparison_summary.tsv"
  output:
    "results/assigner_comparison/output_list.txt"
  shell:
    "ls {input.table} {input.alpha_div} {input.beta_div} {input.tax_compare} > {output}"
