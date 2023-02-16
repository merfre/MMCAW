fastp report for {{ snakemake.wildcards["PATHS"] }} that contains sequence statistics pre and post initial quality control.

Parameters used for fastp were a quality phred value of {{ snakemake.config["qualified_quality_phred"] }}, unqualified percent limit of {{ snakemake.config["unqualified_percent_limit"] }},
average read quality requirement of {{ snakemake.config["average_qual"] }}, minimum read length of {{ snakemake.config["min_length"] }}, front trim of {{ snakemake.config["front_trim"] }}, and a tail trim of {{ snakemake.config["tail_trim"] }}.
