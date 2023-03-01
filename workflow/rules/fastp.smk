### Rule to run fastp for QC and QC reports ###

configfile: "config/config.yaml"

### Pre-QC read report

rule preqc_stats:
  #conda:
    #"../workflow/envs/environment.yaml"
  input:
    fastq_file = expand("resources/{path}.fastq", path=PATHS)
  output:
    report = report("results/qc_reports/sample_pre_qc_report.tsv", caption="report/pre_qc_reports.rst", category="QC Reports")
  shell:
    "seqkit stats {input.fastq_file} -a -T > {output.report}"
    # flag -a denotes all statistics, including quartiles of seq length, sum_gap, N50
    # -T output in machine-friendly tabular format

### Run fastp

rule fastp:
  #conda:
    #"../workflow/envs/environment.yaml"
  input:
    fastq = "resources/{PATHS}.fastq",
    preqc_stats = "results/qc_reports/sample_pre_qc_report.tsv"
  output:
    reads_trimmed = "results/preprocessing/trimmed_filtered/{PATHS}_trimmed_filtered.fastq",
    html = report("results/qc_reports/fastp_reports/{PATHS}.html", caption="report/fastp_reports.rst", category="QC Reports")
  params:
    qualified_quality_phred = config['qualified_quality_phred'],
    unqualified_percent_limit = config['unqualified_percent_limit'],
    average_qual = config['average_qual'],
    min_length = config['min_length'],
    front_trim = config['front_trim'],
    tail_trim = config['tail_trim']
  shell:
    """
    fastp -i {input} -q {params.qualified_quality_phred} -u {params.unqualified_percent_limit} \
    -e {params.average_qual} -l {params.min_length} -f {params.front_trim} -t {params.tail_trim} \
    -y -o {output.reads_trimmed} -h {output.html}
    """
    # flag q denotes quality value that a base is qualified - Default 15 means phred quality >=Q15 is qualified.
    # flag u denotes how many percents of bases are allowed to be unqualified (0~100) - Default 40 means 40%
    # flag e denotes if one read's average quality score <avg_qual, then this read/pair is discarded - Default 0 means no requirement (int [=0])
    # flag l denotes the minimum length requirement
    # flag f denotes the number of bases to trim off the front of the read
    # flag t denotes the number of bases to trim off the tail of the read
    # flag y denotes the activation of low complexity filter
    # Adapter trimming is enabled by default and sequences can be automatically detected for both PE/SE data
    # Overrepresented sequence analysis is disabled by default, to enable --overrepresentation_analysis

### Post-QC read report

rule postqc_stats:
  #conda:
    #"../workflow/envs/environment.yaml"
  input:
    fastq_file = expand("results/preprocessing/trimmed_filtered/{path}_trimmed_filtered.fastq", path=PATHS)
  output:
    report = report("results/qc_reports/sample_post_qc_report.tsv", caption="report/post_qc_reports.rst", category="QC Reports")
  shell:
    "seqkit stats {input.fastq_file} -a -T > {output.report}"
    # flag -a denotes all statistics, including quartiles of seq length, sum_gap, N50
    # -T output in machine-friendly tabular format
