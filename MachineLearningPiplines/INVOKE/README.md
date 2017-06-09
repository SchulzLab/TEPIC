# INVOKE (version 1.0)
-------
Identification of key transcriptional regulators using epigenetics data

## News
09.06.2017: INVOKE is now included in the TEPIC repository

## Introduction
*INVOKE* offers linear regression with various regularisation techniques (Lasso, Ridge, Elastic net) to infer
potentially important transcriptional regulators by predicting gene expression from *TEPIC* TF-gene scores. 

Further details on INVOKE can be found in the provided description [TEPIC-INVOKE.pdf](/docs/TEPIC-INVOKE.pdf).

## Installing TEPIC
To run *INVOKE* the following packages/software must be installed additionally to TEPIC:
* R (minimum version 2.7.4)
* glmnet
* doMC
* parallel
* optional: gplots
* optional: ggplot2

##Required Input:
*INVOKE* can be applied to several samples at once. For each sample, one tab-delimited file has to be provided
containing TF-gene scores and gene expression values per gene. We provide a script to combine TF-gene scores computed
by TEPIC with gene expression data to generate the correct input [integrateData.py](Scripts/integrateData.py).

The gene expression data file is a two column, tab delimited file, including a header, with EnsembleGeneIDs in the the first column and
gene expression values in the second column. An example gene expression file is provided as well: 
[ExampleData/S001S745_ERX616976_GRCh38_hotspot_peaks_20150709_chr1.bed](S001S745_ERX616976_GRCh38_hotspot_peaks_20150709_chr1.bed).

## Using INVOKE
INVOKE can be used in two different ways:

###Combined pipeline with TEPIC
We offer an automated pipeline that includes the computation of TF-gene scores. To run this pipeline, execute the command

    bash runInvokeAnalysis.sh

The setting for this pipeline can be changed in [invokeAnalysis.cfg](invokeAnalysis.cfg).

The following options can be changed in the configuration file:
General parameters:
* open_regions = Bed file containing open chromatin regions (an example is provided).
* gene_Expression_Data = Two column, tab delimited file containing gene expression data (an example is provided).
* outputDirectory = Desired output directory.
* preComputedTEPIC = Filename(s) of precomputed TEPIC results (an example is provided).

*TEPIC* parameters:
* path = Path to TEPIC folder.
* referenceGenome = Filename of the reference genome that should be used.
* pwm = Filename to the PSEM file that should be used.
* cores_TEPIC = Number of cores that can be used.
* geneAnnotation = Filename of the gene annotation file (can be obtained from GENCODE).
* window = Size of the window centered around the TSS of genes.
* peakFeatures = TRUE if peak features (peak length, peak count) should be computed, FALSE otherwise.
* decay = TRUE if exponential decay should be used, FALSE otherwise.
* coverageFile = Filename(s) of bed graph files used to comput the signal in open regions.
* coverageColumn = Number of the column in the file specified by the open_regions option, containing the average signal within the open region.

*INVOKE* parameters:
* rpath = Path to the *INVOKE* script.
* cores_R = Number of cores to be used (suggestion: 1/alpha_Step_Size+1).
* alpha_Step_Size = Step size for optimising alpha in case of elastic net regularisation.
* regularisation = E for Elastic net, L for Lasso, R for Ridge Regression.
* innerCV = Number of inner cross validation folds
* outerCV = Number of outer cross validation folds
* testsize = Size of the test set
* performance = TRUE if model performance should be assessed using an outer cross validation, FALSE otherwise.

For reasons of simplicity and ease of use this mode does not offer all options that are supported by TEPIC and INVOKE.

###Running INVOKE manually
Running an *INVOKE* analysis manually allows full access to all options of *TEPIC* and *INVOKE*. 
There are three main steps the user needs to carry out:
(1) Running *TEPIC*: For details on how to run *TEPIC*, see the *TEPIC* README.
(2) Combine TF gene-scores with gene expression data. We provide a script for this task [integrateData.py](Scripts/integrateData.py).
It can be used with the command

	python Scripts/integrateData.py <TF gene scores> <Gene expression data> <Combined data>

(3) Run the actual model using the script [Scripts/INVOKE.R](Scripts/INVOKE.R). The basic command is

	Rscript Scripts/INVOKE.R --dataDir=<Data Path> --outDir=<Output Path> --response=Expression --regularisation=<E,L,R> --performance=<TRUE,FALSE>

Overall, the following parameters can be specified:
* outDir = Output directory (will be created if it does not exist).
* dataDir = Directory containing the data.
* response = Name of the response variable.
* cores = Number of cores to be used by R (suggestion: 1/alpha_Step_Size+1).
* fixedAlpha = Use a fixed value for alpha instead of optimising alpha in elastic net regularisation.
* alpha = Stepsize for the optimisation of alpha.
* testsize = Size of the test set.
* regularisation = E for Elastic net, L for Lasso, R for Ridge Regression
* innerCV = Number of inner cross validation folds 
* outerCV = Number of outer cross validation folds
* performance = TRUE if model performance should be assessed using an outer cross validation, FALSE otherwise.
* seed = Specify a random SEED to allow reproducability (can be only used if 1 core is used).
* leaveOneOutCV = Flag indicating whether the models should be learned using leave one out cross validation
* asRData = Flag indicating whether feature coefficients should be stored as RData files. 


##Outputs
*INVOKE* always produces the following output:
* A *txt* file with the regression coefficients learned on the entire data set.
* If ggplot2 is available, a bar plot is generated that shows all coefficients with a mean absolute value > 0.025.

If the performance of the model is assessed, *INVOKE* additionally generates:
* A file *Performance_Overview.txt* that holds information on model performance: Pearson correlation, Spearman correlation, and MSE.
* A boxplot showing model performance as well as the input *txt* file to generate the Figure. 
* If gplots is available, a heatmap that shows at most the top 10 positive and the top 10 negative coefficients across the outer cross validation runs as well as a *txt* file with the coefficient values.
* Scatter Plots per sample and outer cross validation run showing the predicted versus the measured gene expression on test data.

## Example
To run a test trial of *INVOKE*, you can run the script runInvokeAnalysis.sh. You can run it with the command

	bash runInvokeAnalysis.sh

We provide precomputed TEPIC results to run the example. If a reference genome and a gene annotation file are provided by the user,
the first step of the pipeline, the run of TEPIC, can be computed as well. Sample open chromatin regions as well as gene expresssion data
are already contained in the folder [ExampleData](ExampleData).

## Citation
If you are using INVOKE please cite:

**Combining transcription factor binding affinities with open-chromatin data for accurate gene expression prediction**
Schmidt et al., Nucleic Acids Research 2016; doi: 10.1093/nar/gkw1061 [full text](http://nar.oxfordjournals.org/content/early/2016/11/29/nar.gkw1061.full) 

Other works that have influenced ours:
> (1) Predicting transcription factor affinities to DNA from a biophysical model, Roider HG, et al., Bioinformatics, 2007.
> (2) ChIP-Seq of transcription factors predicts absolute and differential gene expression in embryonic stem cells, Ouyang Z, et al.,  PNAS, 2009.
> (3) A general concept for consistent documentation of computational analyses, Ebert P, et al.,  Database, 2015.
> (4) JASPAR: an open-access database for eukaryotic transcription factor binding profiles, Sandelin A., et al., Nucleic Acids Research, 2004.
> (5) HOCOMOCO: a comprehensive collection of human transcription factor binding sites models , Kulakovskiy Ivan V., et al., Nucleic Acids Research, 2013.
> (6) Systematic discovery and characterization of regulatory motifs in ENCODE TF binding experiments, Kheradpour P, and Kellis M, Nucleic Acids Research, 2013.
