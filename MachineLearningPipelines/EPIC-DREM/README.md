# EPIC-DREM (version 1.0)
-------
Using epigenetics data in the dynamic regulatory events miner

## News
14.03.2018: We added further details to the repository on how to use EPIC-DREM
21.06.2017: TEPIC now supports to generate input for DREM

## Introduction
*EPIC-DREM* combines time-point specific TF-gene scores computed with *TEPIC* with time-series gene expression data to
identify TFs that can be associated to gene-expression changes occuring over time.

Further details on EPIC-DREM can be found in the provided [description](/docs/Description.pdf) as well as
in the manuscripts (+) and (0).

## Installing TEPIC
To run *EPIC-DREM* the following packages/software must be installed additionally to TEPIC:
* R (minimum version 2.7.4)
* DREM (available [here](http://www.sb.cs.cmu.edu/drem/))

Required external python packages:
* numpy
* numpy.random
* scipy.spatial
* twobitreader

## Required Input:
*INVOKE* can be applied to several samples at once. For each sample, one tab-delimited file has to be provided
containing TF-gene scores and gene expression values per gene. We provide a script to combine TF-gene scores computed
by TEPIC with gene expression data to generate the correct input [integrateData.py](Scripts/integrateData.py).

The gene expression data file is a two column, tab delimited file, including a header, with EnsembleGeneIDs in the the first column and
gene expression values in the second column. An example gene expression file is provided as well: 
[ExampleData/S001S745_ERX616976_GRCh38_hotspot_peaks_20150709_chr1.bed](S001S745_ERX616976_GRCh38_hotspot_peaks_20150709_chr1.bed).

## Using EPIC-DREM
ADD details

## Citation
If you are using EPIC-DREM please cite:

(0) **Temporal epigenomic profiling identifies AHR as dynamic super-enhancer controlled regulator of mesenchymal multipotency**,
Deborah Gerard, Florian Schmidt, Aurelien Ginolhac, Martine Schmitz, Rashi Halder, Peter Ebert, Marcel H. Schulz, Thomas Sauter, Lasse Sinkkonen, bioRxiv, 
[full text](https://www.biorxiv.org/content/early/2017/11/17/183988) 

Other works that have influenced ours:
> (1) Predicting transcription factor affinities to DNA from a biophysical model, Roider HG, et al., Bioinformatics, 2007.
> (2) ChIP-Seq of transcription factors predicts absolute and differential gene expression in embryonic stem cells, Ouyang Z, et al.,  PNAS, 2009.
> (3) A general concept for consistent documentation of computational analyses, Ebert P, et al.,  Database, 2015.
> (4) JASPAR: an open-access database for eukaryotic transcription factor binding profiles, Sandelin A., et al., Nucleic Acids Research, 2004.
> (5) HOCOMOCO: a comprehensive collection of human transcription factor binding sites models , Kulakovskiy Ivan V., et al., Nucleic Acids Research, 2013.
> (6) Systematic discovery and characterization of regulatory motifs in ENCODE TF binding experiments, Kheradpour P, and Kellis M, Nucleic Acids Research, 2013.
> (7) Reconstructing Dynamic Regulatory Maps, Ernst J, et al., Nature-EMBO Molecular Systems Biology, 3:74, 2007.
