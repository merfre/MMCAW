### Rules to search for ARGs in annotated genes using CARD RGI ###

configfile: "config/config.yaml"

### Run CARD resistance gene identifier

rule rgi:
  conda:
    "envs/environment.yaml"

# need to load the card db first (and add this to the db creation workflow):
# do bulk load all reference data from github:
# https://github.com/arpcard/rgi#install-rgi-from-conda
# *** note that the parameter version_number depends upon the versions of WildCARD data downloaded, please adjust accordingly
# *** note that the FASTA filenames plus the parameter version_number depend on the versions of CARD and WildCARD data downloaded, please adjust accordingly
# wget https://card.mcmaster.ca/latest/data
#tar -xvf data ./card.json
#wget -O wildcard_data.tar.bz2 https://card.mcmaster.ca/latest/variants
#mkdir -p wildcard
#tar -xjf wildcard_data.tar.bz2 -C wildcard
#gunzip wildcard/*.gz
#rgi card_annotation -i /path/to/card.json > card_annotation.log 2>&1
#rgi wildcard_annotation -i wildcard --card_json /path/to/card.json
#  -v version_number > wildcard_annotation.log 2>&1
#rgi load \
#  --card_json /path/to/card.json \
#  --debug --local \
#  --card_annotation card_database_v3.2.4.fasta \
#  --card_annotation_all_models card_database_v3.2.4_all.fasta \
#  --wildcard_annotation wildcard_database_v4.0.0.fasta \
#  --wildcard_annotation_all_models wildcard_database_v4.0.0_all.fasta \
#  --wildcard_index /path/to/wildcard/index-for-model-sequences.txt \
#  --wildcard_version 4.0.0 \
#  --amr_kmers /path/to/wildcard/all_amr_61mers.txt \
#  --kmer_database /path/to/wildcard/61_kmer_db.json \
#  --kmer_size 61


# running rgi:
# rgi main OR rgi bwt (metagenomic, beta status)
# for running main rgi:
# --low_quality         use for short contigs to predict partial genes (default: False)
# -i INPUT_SEQUENCE, --input_sequence INPUT_SEQUENCE
#                        input file must be in either FASTA (contig and
#                        protein) or gzip format! e.g myFile.fasta,
#                        myFasta.fasta.gz
# -o OUTPUT_FILE, --output_file OUTPUT_FILE
#                        output folder and base filename
#  -t {contig,protein}, --input_type {contig,protein}
#                        specify data input type (default = contig)
# -a {DIAMOND,BLAST}, --alignment_tool {DIAMOND,BLAST}
#                        specify alignment tool (default = BLAST)
# --include_loose       include loose hits in addition to strict and perfect
#                        hits (default: False)
#  --include_nudge       include hits nudged from loose to strict hits
#                        (default: False)
# --local               use local database (default: uses database in
#                        executable directory) CONFUSING BC HOW TO SPECIFY
# looks like bwt/meta requires short read sequencing with paried reads


# make heatmaps of rgi results as well: rgi heatmap
# -i INPUT, --input INPUT
#                        Directory containing the RGI .json files (REQUIRED)
#  -cat {drug_class,resistance_mechanism,gene_family}, --category {drug_class,resistance_mechanism,gene_family}
#                        The option to organize resistance genes based on a category.
#  -f, --frequency       Represent samples based on resistance profile.
#  -o OUTPUT, --output OUTPUT
#                        Name for the output EPS and PNG files.
#                        The number of files run will automatically
#                        be appended to the end of the file name.(default=RGI_heatmap)
# -clus {samples,genes,both}, --cluster {samples,genes,both}
#                        Option to use SciPy's hiearchical clustering algorithm to cluster rows (AMR genes) or columns (samples).
#  -d {plain,fill,text}, --display {plain,fill,text}
#                        Specify display options for categories (deafult=plain).