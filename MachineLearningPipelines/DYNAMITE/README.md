# DYNAMITE (version 1.0)
-------
Differential analysis to identify novel transcriptional regulators for differentially expressed genes

## News
13.06.2017: DYNAMITE is now included in the TEPIC repository

## Introduction
*DYNAMITE* offers logistic regression with elastic net regularisation to infer
potentially important transcriptional regulators by predicting up/down-regulation for differentially expressed genes from *TEPIC* TF-gene score ratios. 

Further details on DYNAMITE can be found in the provided [description](/docs/Description.pdf).

## Installing TEPIC
To run *DYNAMITE* the following packages/software must be installed additionally to TEPIC:
* R (minimum version 2.7.4)
* glmnet
* doMC
* optional: ggplot2

## Required Input:
*DYNAMITE* can be applied to several samples at once. For each sample, one tab-delimited file has to be provided
containing TF-gene scores and log2 gene expression ratios per gene. We provide  scripts to combine TF-gene scores computed
by TEPIC with log2 gene expression ratios to generate the correct input [integrateData.py](Scripts/integrateData.py) and [prepareForClassificiation.R](Scripts/prepareForClassificiation.R).

The gene expression data file is a two column, tab delimited file, including a header, with EnsembleGeneIDs in the the first column and
log2 gene expression ratios in the second column. An example gene expression file is provided as well: 
[Example](ExampleData/TN_vs_TEM0_000001.txt).

## Using DYNAMITE
DYNAMITE can be used in two different ways:

### Combined pipeline with TEPIC
We offer an automated pipeline that includes the computation of TF-gene scores. To run this pipeline, execute the command

    bash runDYNAMITE.sh ./DYNAMITE.cfg

The setting for this pipeline can be changed in [DYNAMITE.cfg](DYNAMITE.cfg).

The following options can be changed in the configuration file:

General parameters:
* open_regions_Group1 = Bed files with open/active regions for one more multiple samples belonging to group 1.
* open_regions_Group2 = Bed files with open/active regions for one more multiple samples belonging to group 2.
* differential_Gene_Expression = Log2 gene expression ratios computed between group 1 and group 2.
* outputDirectory = Desired output directory.
* existing_TEPIC_Results_Group1 = Path to precomputed TEPIC results of group 1.
* existing_TEPIC_Results_Group2 = Path to precomputed TEPIC results of group 2.

*TEPIC* parameters:
* path = Path to TEPIC folder.
* referenceGenome = Filename of the reference genome that should be used.
* pwm = Filename to the PSEM file that should be used.
* cores_TEPIC = Number of cores that can be used.
* geneAnnotation = Filename of the gene annotation file (can be obtained from GENCODE).
* window = Size of the window centered around the TSS of genes.
* peakFeatures = TRUE if peak features (peak length, peak count) should be computed, FALSE otherwise.
* decay = TRUE if exponential decay should be used, FALSE otherwise.
* coverage_Files_Group1 = bg files containing the signal of the open chromatin assays of group 1 (Note that the order must be the same as in open_regions_Group1).
* coverage_Files_Group2 = bg files containing the signal of the open chromatin assays of group 2 (Note that the order must be the same as in open_regions_Group2).
* coverage_Columns_Group1 = Column containing the signal of the assay in the open region stored in the files listed in open_regions_Group1.
* coverage_Columns_Group2 = Column containing the signal of the assay in the open region stored in the files listed in open_regions_Group2.

*DYNAMITE* parameters:
* rpath = Path to the *DYNAMITE* Rscript.
* cores_R = Number of cores to be used (suggestion: 1/alpha_Step_Size+1).
* alpha_Step_Size = Step size for optimising alpha in elastic net regularisation.
* innerCV = Number of inner cross validation folds
* outerCV = Number of outer cross validation folds
* performance = TRUE if model performance should be assessed using an outer cross validation, FALSE otherwise.
* balanced = TRUE if the data set should be balanced through downsampling for model training and evaluation, FALSE otherwise. 

For reasons of simplicity and ease of use this mode does not offer all options that are supported by TEPIC.
Note that the pipeline does support multiple samples per group. However, it is not possible to compare more than two groups at once. 
This can be done if the *DYNAMITE* Rscript is ran manually. 

### Running DYNAMITE manually
Running an *DYNAMITE* analysis manually allows full access to all options of *TEPIC* and *DYNAMITE*. 
There are three main steps the user needs to carry out:
(1) Running *TEPIC*: For details on how to run *TEPIC*, see the *TEPIC* README.
(2) Combine TF gene-scores with gene expression data. We provide a script for this task [integrateData.py](Scripts/integrateData.py).
It can be used with the command

	python Scripts/integrateData.py <TF gene scores> <Gene expression data> <Combined data>

(3) Discretise the log2 ratios. To do this use the script [prepareForClassificiation.R](Scripts/prepareForClassificiation.R) with the command:

	Rscript Scripts/prepareForClassificiation.R <Input> <Output>
	
(4) Run the actual model using the script [Scripts/DYNAMITE.R](Scripts/DYNAMITE.R). The basic command is

	Rscript Scripts/DYNAMITE.R --dataDir=<Data Path> --outDir=<Output Path> --response=Expression --regularisation=<E,L,R> --performance=<TRUE,FALSE>

## Outputs
*DYNAMITE* always produces the following output:
* A *txt* file with the regression coefficients learned on the entire data set.
* If ggplot2 is available, a bar plot is generated that shows all nonzero regression coefficients.

If the performance of the model is assessed, *DYNAMITE* additionally generates:
* A file *Performance_Overview.txt* that holds information on model performance: Accuracy on Training and Test data, F1 measures.
* A boxplot showing model performance.
* Confusion matrices per sample and outer cross validation run. 
* A heatmap showing at the top positive and negative features

## Example
To run a test trial of *DYNAMITE*, execute the script runDYNAMITE.sh. You can run it with the command

	bash runDYNAMITE.sh ./DYNAMITE.cfg

We provide precomputed TEPIC results to run the example. If a reference genome and a gene annotation file are provided by the user,
the first step of the pipeline, the run of TEPIC, can be carried out as well. Sample open chromatin regions, log2 gene expresssion ratios, and a reduced set of PSEMs (in the interest of runtime)
are already contained in the folder [ExampleData](ExampleData). Note that the genome version used to compute the TF gene-scores in this example is hg19. 

## Citation
If you are using DYNAMITE please cite:

**Epigenomic Profiling of Humand CD4+ T Cells Supports a Linear Differentiation Model and Highlights Molecular Regulators of Memory Development.**
Durek et al., Immunity 2016; doi: 10.1016/j.immuni.2016.10.022 [full text](http://www.sciencedirect.com/science/article/pii/S1074761316304332?via%3Dihub) 

Other works that have influenced ours:
> (1) Predicting transcription factor affinities to DNA from a biophysical model, Roider HG, et al., Bioinformatics, 2007.
> (2) ChIP-Seq of transcription factors predicts absolute and differential gene expression in embryonic stem cells, Ouyang Z, et al.,  PNAS, 2009.
> (3) A general concept for consistent documentation of computational analyses, Ebert P, et al.,  Database, 2015.
> (4) JASPAR: an open-access database for eukaryotic transcription factor binding profiles, Sandelin A., et al., Nucleic Acids Research, 2004.
> (5) HOCOMOCO: a comprehensive collection of human transcription factor binding sites models , Kulakovskiy Ivan V., et al., Nucleic Acids Research, 2013.
> (6) Systematic discovery and characterization of regulatory motifs in ENCODE TF binding experiments, Kheradpour P, and Kellis M, Nucleic Acids Research, 2013.
