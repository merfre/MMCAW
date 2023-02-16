### Rule to run fastp for QC and QC reports ###

rule fastp:
  #conda:
    #"../environment.yml"
  input:
    "{PATHS}.fastq"
  output:
    reads_trimmed = "results/Preprocessing/trimmed_filtered/{PATHS}_trimmed_filtered.fastq",
    html = report("results/QC_reports/fastp_reports/{PATHS}.html", caption="report/fastp_reports.rst", category="QC reports")
  params:
    qualified_quality_phred = config['qualified_quality_phred'],
    unqualified_percent_limit = config['unqualified_percent_limit'],
    average_qual = config['average_qual'],
    min_length = config['min_length'],
    front_trim = config['front_trim'],
    tail_trim = config['tail_trim']
  shell:
    "fastp -i {input} -q {params.qualified_quality_phred} -u {params.unqualified_percent_limit} \
    -e {params.average_qual} -l {params.min_length} -f {params.front_trim} -t {params.tail_trim} \
    -o {output.reads_trimmed} -h {output.html}"
    # flag q denotes quality value that a base is qualified - Default 15 means phred quality >=Q15 is qualified.
    # flag u denotes how many percents of bases are allowed to be unqualified (0~100) - Default 40 means 40%
    # flag e denotes if one read's average quality score <avg_qual, then this read/pair is discarded - Default 0 means no requirement (int [=0])
    # flag l denotes the minimum length requirement
    # flag f denotes the number of bases to trim off the front of the read
    # flag t denotes the number of bases to trim off the tail of the read
    # Adapter trimming is enabled by default and sequences can be automatically detected for both PE/SE data
    # Overrepresented sequence analysis is disabled by default, to enable --overrepresentation_analysis

### Perform multiqc on all fastp reports of trimmed reads

#rule multiqc_post:
  #conda:
    #"../environment.yml"
  #input:
    #expand("results/QC_reports/fastqc_reports_post/{path}_fastqc.zip", path=PATHS)
  #output:
    #"results/QC_reports/multiqc_reports/{LIBRARIES}/multiqc_post_qc.html"
  #params:
    #""  # Optional: extra parameters for multiqc.
  #log:
    #"results/QC_reports/multiqc_reports/{LIBRARIES}/multiqc_post.log"
  #wrapper:
    #config['multiqc_wrapper']
