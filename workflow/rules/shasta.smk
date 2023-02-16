### Rules to run shasta and assemble reads into contigs ###

### Assemble reads into contigs using Shasta

rule shasta:
  #conda:
    #"../environment.yml"
  input:
    "results/Preprocessing/fasta_converted/{PATHS}_trimmed_filtered_humrm.fasta"
  output:
    directory = directory("results/Preprocessing/shasta_results/{PATHS}"),
    fasta = "results/Preprocessing/shasta_results/{PATHS}/Assembly.fasta"
  params:
    config_file = config['shasta_config']
  shell:
    "shasta --input {input} --assemblyDirectory {output.directory} --config {params.config_file}"

### Post assembly reports

rule assembly_stat_report:
  #conda:
    #"../environment.yml"
  input:
    expand("results/Preprocessing/shasta_results/{path}/Assembly.fasta", path=PATHS)
  output:
    report = report("results/QC_reports/assembly_stat_reports/stats_report.tsv", caption="report/assembly_stat_reports.rst", category="QC reports")
  shell:
    "seqkit stats {input} -a -T > {output.report}"
    # flag -a denotes all statistics, including quartiles of seq length, sum_gap, N50
    # -T output in machine-friendly tabular format
