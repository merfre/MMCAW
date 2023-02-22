### Rules to run phylophlan and transform the outputs ###

rule phylophlan:
  #conda:
    #"../environment.yml"
  input:
  output:
  shell:

  # first step is to build config file - after databse is built
  # second is to run phylophlan_metagenomic with directory output and tsv output
  # third is to use phylophlan_draw_metagenomic for a heatmap
  # fourth is to create a tree?

#phylophlan_metagenomic \
#    -i input_metagenomic \
#    -o output_metagenomic \
#    --nproc 4 \
#    -n 1 \
#    -d SGB.Jan19 \
#    --verbose 2>&1 | tee logs/phylophlan_metagenomic.log

#making heatmap:

#phylophlan_draw_metagenomic \
#    -i output_metagenomic.tsv \
#    -o output_heatmap \
#    --map bin2meta.tsv \
#    --top 20 \
#    --verbose 2>&1 | tee logs/phylophlan_draw_metagenomic.log

  # automatically retrieving reference genomes and species-specific sets of UniRef90 proteins.
  # updated when the users specify the --database_update parameter.
  # Database files comprise the sets of precomputed species-specific UniRef90 proteins, the list of
  # available genomes from GenBank, and the SGB release.
# phylophlan_setup_database -i /home/639893/Subway_2015_Reproduction/data/databases/ncbi_database -d phylophlan_database -o resources/databases -t n --verbose
