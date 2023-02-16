### Rules to remove human sequences from the nanopore data ###

configfile: "config/config.yaml"

### Align reads to reference human sequences with minimap2

rule minimap2_align_human:
  #conda:
    #"../environment.yml"
  input:
    "results/Preprocessing/trimmed_filtered/{PATHS}_trimmed_filtered.fastq"
  output:
    "results/Preprocessing/minimap2/{PATHS}_human_alignment.sam"
  params:
    filtering_reference = config['filtering_reference']
  shell:
    "minimap2 -ax map-ont {params.filtering_reference} {input} > {output}"
    # the options -ax map-ont instruct minimap2 to align the remove_human_sequences
    # and specifies that it is nanopore data (rather than pacbio)

### Use samtools to remove the reads that matched the reference

rule remove_human_sequences:
  #conda:
    #"../environment.yml"
  input:
    "results/Preprocessing/minimap2/{PATHS}_human_alignment.sam"
  output:
    "results/Preprocessing/trimmed_filtered_humrm/{PATHS}_trimmed_filtered_humrm.fastq"
  shell:
    "samtools view -buSh -f 4 {input} | samtools fastq - > {output}"
    # -f denotes extracting only reads that match that sam flag (aka un-mapped reads)
    # -b is bam output, -u is uncompressed, -s is sam input, -h is to include the header in the output

### Convert samples from fastq to fasta with seqkit

rule fasta_conversion:
  #conda:
    #"../environment.yml"
  input:
    "results/Preprocessing/trimmed_filtered_humrm/{PATHS}_trimmed_filtered_humrm.fastq"
  output:
    "results/Preprocessing/fasta_converted/{PATHS}_trimmed_filtered_humrm.fasta"
  shell:
    "seqkit fq2fa {input} -o {output}"

### Post human removed QC reports

rule qc_report_humrm:
  #conda:
    #"../environment.yml"
  input:
    expand("results/Preprocessing/fasta_converted/{path}_trimmed_filtered_humrm.fasta", path=PATHS)
  output:
    report = report("results/QC_reports/humrm_qc_reports/stats_report.tsv", caption="report/humrm_qc_reports.rst", category="QC reports")
  shell:
    "seqkit stats {input} -a -T > {output.report}"
    # flag -a denotes all statistics, including quartiles of seq length, sum_gap, N50
    # -T output in machine-friendly tabular format
