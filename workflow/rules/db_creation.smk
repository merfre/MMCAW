### Rule to load the database creation workflow as a subworkflow ###

subworkflow db_creation:
  workdir:
    "rules/Database_Creation_Workflow"
  snakefile:
    "rules/Database_Creation_Workflow/workflow/Snakefile"
  configfile:
    "workflow/rules/Database_Creation_Workflow/config/config.yaml"


rule run_subworkflow:
  input:
    blast = db_creation("results/blast_db.done"),
    cat = db_creation("results/taxdump_version.txt"),
    taxdump= db_creation("results/cat_db.done"),
    reference = db_creation("results/reference.done")
  output:
    touch("results/database_creation.done")
  shell:
    "ls ./"
