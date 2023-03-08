### Rules to search for ARGs in annotated genes using CARD RGI ###

configfile: "config/config.yaml"

### Run CARD resistance gene identifier

rule rgi:
  conda:
    "envs/environment.yaml"