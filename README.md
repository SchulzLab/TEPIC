# TEPIC (version 2.0)
-------
Annotation of genomic regions using Transcription factor (TF) binding sites and epigenetic data.

##News
08.06.2017: Version 2.0 of TEPIC is available.

## Introduction
*TEPIC* segments the genome into user specified regions and annotates those with TF binding using TRAP (1). 
These predictions are aggregated to gene scores. 
Within this aggregation TEPIC offers exponential decay (2) and scaling of TF region scores using the signal of an open chromatin assay.

With version 2 of TEPIC, we introduced new features:
* Affinities can now be normalised for peak length during TF-gene score computation.
* The length of the PSEMs can be considered in the normalisation.
* We introduced features for peak length and peak counts.
* Scaling can now be performed in two ways: The traditional way as proposed in the TEPIC manuscript by directly multiplying
peak signal and TF affinities or by generating a separate signal feature.

## Installing TEPIC
To run *TEPIC* the following packages/software must be installed:
* Python (minimum version 2.7)
* [bedtools](https://github.com/arq5x/bedtools2)
* A g++ compiler supporting openmp to use the parallel implementation of TRAP.

To compile the C++ version of TRAP execute the script
	[Code/compileTRAP.sh](Code/compileTRAP.sh).

## Position specific energy matrices
The position weight matrices used in the *TEPIC* manuscript are stored in the file
	[PWMs/pwm_vertebrates_jaspar_uniprobe_original.PSEM](PWMs/pwm_vertebrates_jaspar_uniprobe_original.PSEM).
An extended set of pwms is also available for human, mouse, rat, drosophila melanogaster, and Caenorhabditis elegans.
We collected motifs from *JASPAR* (4), *HOCOMOCO* (5), and the *Kellis Lab ENCODE Motif database* (6).
* The human set contains 515 *JASPAR Vertebrata* matrices, 81 *Hocomoco human* matrices, and 130 matrices from the *Kellis Lab database*.
* The mouse set contains 499 *JASPAR Vertebrata* matrices, 67 *Hocomoco mouse* matrices , and 121 matrcies from the *Kellis Lab database*.
* The rat set contains 489 *JASPAR Vertebrata* matrices, 67 *Hocomoco mouse* matrices, and 121 matrices from the *Kellis Lab database*.
* The drosophila melanogaster set contains 129 *JASPAR* matrices retrieved from *JASPAR Insecta*.
* The Caenorhabditis elegans set contains 26 *JASPAR* matrices retrieved from *JASPAR Nematoda*.

Additional position weight matrices can be transformed to a usable format using 
	[Code/PSCM_to_PSEM.cpp](Code/PSCM_to_PSEM.cpp).
This program converts matrices in TRANSFAC format to the energy format used by TRAP. 
Details on the parameters used for conversion can be found in the header of the provided files.

## Using TEPIC
To start TEPIC, run the script *TEPIC.sh*

    ./TEPIC.sh

The following parameters are required to run TEPIC:

* -g The reference genome in plain (uncompressed) FASTA format with Ensembl-style chromosome names (i.e., without "chr" prefix).
* -b Regions the user wants to be annotated; chromosome naming compatible to the reference genome file.
* -o Prefix of the output files.
* -p File containing position specific energy matrices (PSEM).

The optional parameters are:

* -a Genome annotation file (gtf). All genes contained in this file will be annotated. The file must have the original format provided by gencode. 
* -w Size of the window around the TSS of genes.
* -d Signal of the open chromatin assay in bg format. Used to compute the average per peak coverage within the regions specified in -b
* -e Deactivates the exponential decay
* -n Indicates that the file in -b contains the average signal in the peaks in the specified column. In this case the -d option is not required to obtain scaled TF affinities.
* -c Number of cores used within TRAP.
* -f A gtf file containing genes of interest. Only regions contained in the file specified by the -b option that are within the window specified by the -w option around these genes will be annotated.
* -y Flag indicating whether the entire gene body should be annotated with TF affinities. A window of half the size of the -w option will be additionaly considered upstream of the genes TSS.
* -l Flag to be set if affinities should not be normalised by peak length.
* -u Flag to be set if peak features for peak length and peak counts should not be generated.
* -x If -d or -n is used together with this flag, the original (Decay-)Scaling formulation of TEPIC is used to compute gene-TF scores.
* -m Path to a tab delimited file containing the length of the used PSEMs. This is incorporated in normalising peak length.
* -s Flag to be set if a sparse matrix representation should be used.
* -z Flag indicating that the output of TEPIC should be zipped.

Depending on the used arguments, TEPIC produces files containing:
* TF affinities for all user specified regions.
* Scaled TF affinities for all user specified regions.
* TF affinities for all genes contained in the annotation file.
* Scaled TF affinities for all genes contained in the annotation file.
* A file containing the factors used to scale the original TF affinities.
* Peak length, peak counts as well as the average signal within a peak. 

Each run of TEPIC generates an *analysis meta datafile (amd)* containing all parameters, files, and outputs associated with the last run of TEPIC.
Together with the provided process xml file, the executed command lines  can be reconstructed (3). We provide amd files in the folder
*MetaData*. These correspond to the gene scores of the *50kb* and *50kb-S* annotation introduced in the *TEPIC* manuscript.

Note that the input files **have to** have unix file endings. Using the scaled features is currently only supported for homo sapiens, mus musculus, and
rattus norvegicus.

## Example
To run a test trial of *TEPIC*, you can use the data provided in the *Test* folder. You can run it with the command

	./TEPIC.sh -g ../Example/example_sequence.fa -b ../Example/example_regions.bed -o TEPIC-Example -p ../PWMs/pwm_vertebrates_jaspar_uniprobe_original.PSEM -a ../Example/example_annotation.gtf -w 3000 -e

This will generate gene scores for the genes contained in *example_annotation.gtf*, using a window of size 3000bp, all pwms contained in *pwm_vertebrates_jaspar_uniprobe_converted.PSEM*, and without 
exponential decay. 

Additionally, we provide a script to test several annotation versions of TEPIC. Execute the script

	bash runTestCases.sh

to compute the trial cases. 

## Citation
If you are using TEPIC please cite:

**Combining transcription factor binding affinities with open-chromatin data for accurate gene expression prediction**
Schmidt et al., Nucleic Acids Research 2016; doi: 10.1093/nar/gkw1061 [full text](http://nar.oxfordjournals.org/content/early/2016/11/29/nar.gkw1061.full) 

Other works that have influenced ours:
> (1) Predicting transcription factor affinities to DNA from a biophysical model, Roider HG, et al., Bioinformatics, 2007.

> (2) ChIP-Seq of transcription factors predicts absolute and differential gene expression in embryonic stem cells, Ouyang Z, et al.,  PNAS, 2009.

> (3) A general concept for consistent documentation of computational analyses, Ebert P, et al.,  Database, 2015.

> (4) JASPAR: an open-access database for eukaryotic transcription factor binding profiles, Sandelin A., et al., Nucleic Acids Research, 2004.
 
> (5) HOCOMOCO: a comprehensive collection of human transcription factor binding sites models , Kulakovskiy Ivan V., et al., Nucleic Acids Research, 2013.

> (6) Systematic discovery and characterization of regulatory motifs in ENCODE TF binding experiments, Kheradpour P, and Kellis M, Nucleic Acids Research, 2013.
