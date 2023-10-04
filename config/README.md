Describe how to configure the workflow (using config.yaml and maybe additional files).
All of them need to be present with example entries inside of the config folder.

This file is split into four sections:


1. The first is at the top of the file and contains general information for the workflow, including the location of the desired libraries, the metadata file for the samples, and the software environment to use for analysis.


2. The second section is titled 'Analysis options' and allows certain steps in analysis to be toggled on or off. Each analysis option is described below its title and describes which software is used for that step. When the word 'True' is next to the analysis option it will be included in the next workflow run. When 'False' is next to the option it will be excluded and the software described in that step will not be used for the next run.


3. The next section is below the second and labeled 'Database locations'. In this section the location of the reference databases for analysis is specified. For this workflow, there is one reference required, a  Kraken2 database for taxonomy assignment.


4. The final section, located below the third, is titled 'Parameters'. This section contains the list of adjustable parameters used in this workflow. Each parameter has a title, a value assigned to it, and a description. The parameters are organized by the software that uses them. For instance, fastp is used for initial quality control and has six adjustable parameters for its performance.