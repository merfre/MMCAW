### Rule to run MetaBAT2 and create bins of reads ###

configfile: "config/config.yaml"

### Optional step prior to running MetaBAT2

# MetaBAT2 allows for the use of a summary bam file to be used when binning
# It is a file having mean and variance of base coverage depth (tab delimited;
# the first column should be contig names, and the first row will be
# considered as the header and be skipped)
# The effect of including or excluding this file is not described by MetaBAT2
# To produce the optional file: jgi_summarize_bam_contig_depths --outputDepth depth.txt *.bam

### Run metabat2

rule metabat2:
  #conda:
  #"../environment.yml"
  input:
    "results/Preprocessing/flye_results/{PATHS}/assembly.fasta"
  output:
    "results/Preprocessing/metabat2/{PATHS}_bins/bin"
  params:
    MetaBAT2_min_size = config['MetaBAT2_min_size'],
    MetaBAT2_max_perc = config['MetaBAT2_max_perc'],
    MetaBAT2_min_edge = config['MetaBAT2_min_edge'],
    MetaBAT2_max_edge = config['MetaBAT2_max_edge'],
    MetaBAT2_tnf_prob = config['MetaBAT2_tnf_prob'],
    MetaBAT2_min_coverage = config['MetaBAT2_min_coverage'],
    MetaBAT2_min_coverage_sum = config['MetaBAT2_min_coverage_sum'],
    MetaBAT2_min_bin_size = config['MetaBAT2_min_bin_size'],
    threads = config['threads']
  shell:
    "metabat2 -i {input} --minContig {params.MetaBAT2_min_size} \
    --maxP {params.MetaBAT2_max_perc} --minS {params.MetaBAT2_min_edge} \
    --maxEdges {params.MetaBAT2_max_edge} --pTNF {params.MetaBAT2_tnf_prob} \
    --minCV {params.MetaBAT2_min_coverage} --numThreads {params.threads} \
    --minCVSum {params.MetaBAT2_min_coverage_sum} \
    --minClsSize {params.MetaBAT2_min_bin_size} -o {output}"
    # Descriptions for all parameters are in the configuration file next to the parameters
    # input os the fasta format of filered and trimmed samples
    # output is a base file name and path for each bin. The default output is fasta format.
    # Use -l option to output only contig names.
    # Set --noAdd when added small or leftover contigs cause too much contamination.
    # if the optional summary file is being used the flag is -a depth.txt
