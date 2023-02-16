### Rule to get tax database for recentrifuge from NCBI ###

configfile: "config/config.yaml"

# This is dave's script from Tapirs to retrieve the recntrifuge input:

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

### Retrieve taxdump directory from NCBI

rule get_taxdump:
  #conda:
    #"../environment.yml"
  input:
    FTP.remote("ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.zip")
  output:
    taxdump = directory("resources/databases/taxdump"),
    rankedlineage = "resources/databases/taxdump/rankedlineage.dmp",
    nodes = "resources/databases/taxdump/nodes.dmp",
    taxdump_version = "resources/databases/taxdump_version.txt"
  shell:
    """
    echo "Taxdump retrieved on `date`" > {output.taxdump_version};
    unzip {input} -d {output.taxdump};
    """
