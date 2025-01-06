### Rules to remove human sequences from the nanopore data ###

configfile: "config/config.yaml"

### Align reads to reference human sequences with minimap2

rule minimap2_align_human:
  conda:
    "envs/environment.yaml"
  input:
    "results/preprocessing/trimmed_filtered/{PATHS}_trimmed_filtered.fastq"
  output:
    "results/preprocessing/minimap2/{PATHS}_human_alignment.sam"
  benchmark:
    "benchmarks/{PATHS}_minimap2_align_human.tsv"
  params:
    filtering_reference = config['filtering_reference']
  shell:
    "minimap2 -ax map-ont {params.filtering_reference} {input} > {output}"
    # the options -ax map-ont instruct minimap2 to align the remove_human_sequences
    # and specifies that it is nanopore data (rather than pacbio)

### Use samtools to remove the reads that matched the reference

rule remove_human_sequences:
  conda:
    "envs/environment.yaml"
  input:
    "results/preprocessing/minimap2/{PATHS}_human_alignment.sam"
  output:
    "results/preprocessing/trimmed_filtered_humrm/{PATHS}_trimmed_filtered_humrm.fastq"
  benchmark:
    "benchmarks/{PATHS}_remove_human_sequences.tsv"
  shell:
    "samtools view -buSh -f 4 {input} | samtools fastq - > {output}"
    # -f denotes extracting only reads that match that sam flag (aka un-mapped reads)
    # -b is bam output, -u is uncompressed, -s is sam input, -h is to include the header in the output

### Convert samples from fastq to fasta with seqkit

rule fasta_conversion:
  conda:
    "envs/environment.yaml"
  input:
    "results/preprocessing/trimmed_filtered_humrm/{PATHS}_trimmed_filtered_humrm.fastq"
  output:
    "results/preprocessing/fasta_converted/{PATHS}_trimmed_filtered_humrm.fasta"
  benchmark:
    "benchmarks/{PATHS}_fasta_conversion.tsv"
  shell:
    "seqkit fq2fa {input} -o {output}"

### Post human removed QC reports

rule qc_report_humrm:
  conda:
    "envs/environment.yaml"
  input:
    expand("results/preprocessing/fasta_converted/{path}_trimmed_filtered_humrm.fasta", path=PATHS)
  output:
    report = report("results/qc_reports/humrm_qc_report.tsv", caption="report/humrm_qc_reports.rst", category="QC reports")
  benchmark:
    "benchmarks/qc_report_humrm.tsv"
  shell:
    "seqkit stats {input} -a -T > {output.report}"
    # flag -a denotes all statistics, including quartiles of seq length, sum_gap, N50
    # -T output in machine-friendly tabular format
