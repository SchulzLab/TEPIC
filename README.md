# TEPIC (version 2.0)
-------
Annotation of genomic regions using Transcription factor (TF) binding sites and epigenetic data. Learning of key regulatory TFs in  individual cell types or learning of discriminatory TFs that show a difference in regulation between two cell types.

## News
17.08.2017: TEPIC TF-gene scores can now be binarisied using background regions provided by the user.

21.06.2017: TEPIC TF-gene scores can be binarisied using a TF and tissue specific affinity cut off. This can be combined with the dynamic networks learner [DREM](http://www.sb.cs.cmu.edu/drem/) to build gene regulatory networks that make appropriate use of time-specific epigenetics data. Further details on this new feature are available in the [description](docs/Description.pdf).

13.06.2017: [DYNAMITE](https://github.com/SchulzLab/TEPIC/tree/master/MachineLearningPipelines/DYNAMITE) our workflow for learning differential TF regulation is now included in the repository.

09.06.2017: Version 2.0 of TEPIC is available.
With version 2 of TEPIC, we introduced new features:
* We extended the set of PSEMs.
* TF affinities are computed using a C++ implementation of TRAP.
* Affinities can be normalised by peak length during TF-gene score computation.
* The length of the PSEMs can be considered in the normalisation.
* We introduced features for peak length and peak counts.
* Scaling can be performed in two ways: The original way as proposed in the TEPIC manuscript by directly multiplying
peak signal and TF affinities or by generating a separate signal feature.

Further, the repository now includes the code required to learn linear models from TF gene scores to predict gene expression.
For further details, please see the [INVOKE](https://github.com/SchulzLab/TEPIC/tree/master/MachineLearningPipelines/INVOKE) section.

## Introduction
*TEPIC* segments the genome into user specified regions and annotates those with TF binding using TRAP (1). 
These predictions are aggregated to gene scores. 
Within this aggregation TEPIC offers exponential decay (2) and scaling of TF region scores using the signal of an open chromatin assay.
While computing gene TF scores, TEPIC can perform normalisation for peak length (optionally correcting for the length of the binding motifs as well)
and produces separate features for peak length, peak count and/or peak signal. These can be used in downstream applications, e.g. to determine
the influence of chromatin accessiblity on gene expression, without considering detailed information on TF binding. 
In case a binary assignment of TF binding is required, TEPIC allows to apply a TF specific affinity threshold derived from random genomic sequences that
show similar characteristica (GC content, length) then the provided regions. 

Further details on TEPIC can be found in the provided [description](docs/Description.pdf).

## Installing TEPIC
To run *TEPIC* the following packages/software must be installed:
* Python (minimum version 2.7)
* [bedtools](https://github.com/arq5x/bedtools2)
* A g++ compiler supporting openmp to use the parallel implementation of TRAP.

To compile the C++ version of TRAP and to install possibly missing R-packages for downstream processing execute the script
	[Code/compile_TRAP_install_R_packages.sh](Code/compile_TRAP_install_R_packages.sh).

To use the script [findBackground](Code/findBackground.py), which is necessary to compute TF specific affinity thresholds, the following python libraries are required:
* numpy
* scipy
* twobitreader


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
Files holding the length of the provided PSEMs are provided too. 

Additional position weight matrices can be transformed to a usable format using 
	[Code/PSCM_to_PSEM.cpp](Code/PSCM_to_PSEM.cpp).
This program converts matrices in TRANSFAC format to the energy format used by TRAP. 
Details on the parameters used for conversion can be found in the header of the provided files.

In TRANSFAC format, a matrix has to have the following structure:
XX
ID <RunningNumber>        <TF-Name>
XX
P0        A         C         G         T
1         0         93        6         1
2         7         81        1         12
.
.
.
14        1         3         95        0
XX
//   


In the PSEM format, a matrix has to be in the following structure:
><RunningNumber>  <TF-Name> lnR0: <value>
1.56945   -0.108976 1.46047   0
5.06003   4.54008   0         4.06982
.
.
.
4.59839   4.07844   4.07844   0
5.11834   4.59839   0         5.11834



## Using TEPIC
To start TEPIC, run the script *TEPIC.sh*

    ./TEPIC.sh

The following parameters are required to run TEPIC:

* -g The reference genome in plain (uncompressed) FASTA format with Ensembl-style chromosome names (i.e., without "chr" prefix). If a "chr" prefix is present, use the -j option. 
* -b Regions the user wants to be annotated; chromosome naming compatible to the reference genome file.
* -o Prefix of the output files.
* -p File containing position specific energy matrices (PSEM).

The optional parameters are:

* -a Genome annotation file (gtf). All genes contained in this file will be annotated. The file must have the original format provided by gencode, gzipped files are not supported. 
* -w Size of the window around the TSS of genes.
* -d Signal of the open chromatin assay in bg format. Used to compute the average per peak coverage within the regions specified in -b.
* -e Boolean controlling exponential decay (default TRUE).
* -n Indicates that the file in -b contains the average signal in the peaks in the specified column. In this case the -d option is not required to obtain scaled TF affinities.
* -c Number of cores used within TRAP.
* -f A gtf file containing genes of interest. Only regions contained in the file specified by the -b option that are within the window specified by the -w option around these genes will be annotated. The file must have the original format provided by gencode, gzipped files are not supported.
* -y Flag indicating whether the entire gene body should be annotated with TF affinities. A window of half the size of the -w option will be additionaly considered upstream of the genes TSS.
* -l Flag to be set if affinities should not be normalised by peak length.
* -u Flag to be set if peak features for peak length and peak counts should not be generated.
* -x If -d or -n is used together with this flag, the original (Decay-)Scaling formulation of TEPIC is used to compute gene-TF scores.
* -m Path to a tab delimited file containing the length of the used PSEMs. This is incorporated in normalising peak length.
* -z Flag indicating that the output of TEPIC should be zipped.
* -k Path to a file containing background regions provided by the user. This option can not be used together with the -r option. 
* -r Path to a 2bit representation of the reference genome. This is required to compute a TF specific affinity threshold as well as a binary and sparse TF-gene interaction list. This can not be used together with the -k option. 
* -v p-value cut off used to determine a cut off to derive a binary score for TF binding (default 0.05).
* -i minutes that should be spend at most per chromosome to find matching random regions (default 3).
* -j Flag indicating that the reference genome contains a chr prefix. 

Depending on the used arguments, TEPIC produces files containing:
* TF affinities for all user specified regions.
* Scaled TF affinities for all user specified regions.
* TF affinities for all genes contained in the annotation file.
* Scaled TF affinities for all genes contained in the annotation file.
* A file containing the factors used to scale the original TF affinities.
* TF affinities along with features for peak length, peak counts and/or the average signal within a peak. 
* Thresholded TF affinities and TF-gene scores.
* A sparse representation that contains only those TF-gene interactions with affinities above a affinity threshold derived from random genomic sequences.

Each run of TEPIC generates an *analysis meta datafile (amd)* containing all parameters, files, and outputs associated with the last run of TEPIC.
Together with the provided process xml file, the executed command lines  can be reconstructed (3). We provide amd files in the folder
*MetaData*. These correspond to the gene scores of the *50kb* and *50kb-S* annotation introduced in the *TEPIC* manuscript.

Note that the input files **have to** have unix file endings. Using bed graph files to compute the signal in peaks is currently only supported for homo sapiens, mus musculus, and
rattus norvegicus.

## Example
To run a test trial of *TEPIC*, you can use the data provided in the *Test* folder. You can run it with the command

	./TEPIC.sh -g ../Test/example_sequence.fa -b ../Test/example_regions.bed -o TEPIC-Example -p ../PWMs/pwm_vertebrates_jaspar_uniprobe_original.PSEM -a ../Test/example_annotation.gtf -w 3000 -e FALSE

This will generate gene scores for the genes contained in *example_annotation.gtf*, using a window of size 3000bp, all pwms contained in *pwm_vertebrates_jaspar_uniprobe_converted.PSEM*, and without 
exponential decay. 

Additionally, we provide a script to test several annotation versions of TEPIC. Execute the script

	bash runTestCases.sh

to compute multiple trial cases.

## Citation
If you are using TEPIC and/or [INVOKE](https://github.com/SchulzLab/TEPIC/tree/master/MachineLearningPipelines/INVOKE) please cite:

**Combining transcription factor binding affinities with open-chromatin data for accurate gene expression prediction**
Schmidt et al., Nucleic Acids Research 2016; doi: 10.1093/nar/gkw1061 [full text](http://nar.oxfordjournals.org/content/early/2016/11/29/nar.gkw1061.full) 

If you are using [DYNAMITE](https://github.com/SchulzLab/TEPIC/tree/master/MachineLearningPipelines/DYNAMITE) please also cite:

**Epigenomic Profiling of Human CD4+ T Cells Supports a Linear Differentiation Model and Highlights Molecular Regulators of Memory Development**
Durek et al. Cell Immunity, Volume 45, Issue 5, 15 November 2016, [full text](http://www.cell.com/immunity/fulltext/S1074-7613(16)30433-2)


Other works that have influenced ours:
> (1) Predicting transcription factor affinities to DNA from a biophysical model, Roider HG, et al., Bioinformatics, 2007.

> (2) ChIP-Seq of transcription factors predicts absolute and differential gene expression in embryonic stem cells, Ouyang Z, et al.,  PNAS, 2009.

> (3) A general concept for consistent documentation of computational analyses, Ebert P, et al.,  Database, 2015.

> (4) JASPAR: an open-access database for eukaryotic transcription factor binding profiles, Sandelin A., et al., Nucleic Acids Research, 2004.
 
> (5) HOCOMOCO: a comprehensive collection of human transcription factor binding sites models , Kulakovskiy Ivan V., et al., Nucleic Acids Research, 2013.

> (6) Systematic discovery and characterization of regulatory motifs in ENCODE TF binding experiments, Kheradpour P, and Kellis M, Nucleic Acids Research, 2013.
