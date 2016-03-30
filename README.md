#TEPIC
------
Annotation of genomic regions using Transcription factor (TF) binding sites and epigenetic data.

##Introduction
*TEPIC* segments the genome into user specified regions and annotates those with TF binding using TRAP (Roider et al., 2007). 
These predictions are aggregated to gene scores. 
Within this aggregation TEPIC offers exponential decay (Ouyang et al., 2009) and scaling of TF region scores using the signal of an open chromatin assay.

##Using TEPIC
To start TEPIC, run the script *TEPIC.sh*

    ./TEPIC.sh

The following parameters are required to run TEPIC:

* -g The reference genome.
* -b Regions the user want to be annotated.
* -o Prefix of the output files.
* -p File containing position weight matrices (PWMs) in Jaspar format.

The optional parameters are:

* -a Genome annotation file (gtf). All genes contained in this file will be annotated.
* -w Size of the window around the TSS of genes.
* -d Signal of the open chromatin assay in bg format. Used to compute the average per peak coverage within the regions specified in -b
* -e Deactivates the exponential decay
* -n Indicates that the file in -b contains the average signal in the peaks in the specified column (default 6). In this case the -d option is not required to obtain scaled TF affinities.
* -c Number of cores used within TRAP.

Depending on the used arguments, TEPIC produces the files containing:

* TF affinities for all user specified regions.
* Scaled TF affinities for all user specified regions.
* TF affinities for all genes contained in the annotation file.
* Scaled TF affinities for all genes contained in the annotation file.
* A file holding all regions which were annotated.
* A file containing the factors used to scale the original TF affinities.


Each run of TEPIC generates an *analysis meta datafile* containing all parameters, files, and outputs associated with the last run of TEPIC.
Together with the provided process xml file, the executed command lines  can be reconstructed (Ebert et al., 2015).
