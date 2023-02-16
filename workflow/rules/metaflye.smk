### Rule to assemble nanopore reads with metaflye ###

### Assemble reads into contigs using metaflye

rule metaflye:
  #conda:
    #"../environment.yml"
  input:
    "results/Preprocessing/fasta_converted/{PATHS}_trimmed_filtered_humrm.fasta"
  output:
    directory = directory("results/Preprocessing/flye_results/{PATHS}"),
    fasta = "results/Preprocessing/flye_results/{PATHS}/assembly.fasta"
    # requires fasta in output to be input for next rules
  params:
    read_type = config['read_type'],
    min_overlap = config['minimum_overlap']
  shell:
    """
    flye {params.read_type} {input} --out-dir {output.directory} --meta \
    --min-overlap {params.min_overlap}
    """
    # output is assembly.fasta and assembly information in txt in directory

# add info to report and stats report

### Post assembly reports

rule assembly_stat_report:
  #conda:
    #"../environment.yml"
  input:
    expand("results/Preprocessing/flye_results/{path}/assembly.fasta", path=PATHS)
  output:
    report = report("results/QC_reports/assembly_stat_reports/stats_report.tsv", caption="report/assembly_stat_reports.rst", category="QC reports")
  shell:
    "seqkit stats {input} -a -T > {output.report}"
    # flag -a denotes all statistics, including quartiles of seq length, sum_gap, N50
    # -T output in machine-friendly tabular format
