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
in the manuscripts (1,2).

## Installing TEPIC
To run *EPIC-DREM* the following packages/software must be installed additionally to TEPIC:
* R (minimum version 2.7.4)
* DREM (available [here](http://www.sb.cs.cmu.edu/drem/))

Required external python packages are:
* numpy
* numpy.random
* scipy.spatial
* twobitreader

## Required options:
To run *EPIC-DREM*, we need to compute TF binding predictions for all time points of interest using *TEPIC*.
Due to memory limitations of the DREM software, thresholded TF affinities must be computed. This means that
all TF affinities below a certain threshold will be set to zero. The threshold is determined using TF affinities computed
on a random sequence set and a user defined p-value cut-off (*-v* option). Details to the thresholding are provided in the [documentation](/docs/Description.pdf).
These sequences can either be provided by the user (*-k* option), or be generated automatically by *TEPIC* (*-r* option).

Using the *-r* option, TEPIC will compute a set of background regions that matches length and GC content of the provided candidate binding sites. 
The options specifies the refererence genome in 2bit format. The *-i* option can be used to set the maximum time spent per chromosome, measured in minutes, to find matching background sequences.
Using the *-k* option, a set of user defined background sequences can be specified in bed format. The options *-k* and *-r* are mutually exclusive.

To run *EPIC-DREM* gene expression data must be available for each time-point as well.

## Input to DREM
DREM requires the following input:

* Gene expression data in a tab seperated file (including a header) with the structure:

				Timepoint1	Timepoint2 	...	Timepoint<N>
		GeneID	value		value		...	value

* TF-gene links in a tab seperated file (including a header) with the structure:

		TF	GeneID	score	Timepoint

Note that score can be set to 1 (e.g. using binary TF affinities, as done in (0)), or using the actual TF affinities.
*TEPIC* automatically provides columnes 1 to 3 of the *TF-gene links* file when either the *-k* or *-r* is used option. 
The user (1) has to aggreate these files and (2) has to add the timepoint label. Please also note that the GeneIDs used in the
*expression* file have to match the GeneIDs used in the *TF-gene link* file. 

##Output of EPIC-DREM
The *TEPIC* part of EPIC-DREM provides thresholded TF affinity files. In detail, these are in addition to the standard output

* thresholded TF affinities in all provided genomic regions,
* a file with binary TF-gene assignments which can be further processed as input for DREM,
* continous TF-gene scores using the tresholded affinity values.

the *DREM* software can be used to obtain:

* split tables per timepoint indicating essential regulators at specific splits and time-points,
* gene tables per path indicating coexpressed/coregulated genes,
* the network of splits itself.

Further details can be found in (2).

## Citation
If you are using EPIC-DREM please cite:

(1) **Temporal epigenomic profiling identifies AHR as dynamic super-enhancer controlled regulator of mesenchymal multipotency**,
Deborah Gerard, Florian Schmidt, Aurelien Ginolhac, Martine Schmitz, Rashi Halder, Peter Ebert, Marcel H. Schulz, Thomas Sauter, Lasse Sinkkonen; bioRxiv; 2017
[full text](https://www.biorxiv.org/content/early/2017/11/17/183988) 

(2) **DREM 2.0: Improved reconstruction of dynamic regulatory networks from time-series expression data **,
Marcel H. Schulz, William E Devanny, Anthony Gitter, Shan Zhong, Jason Ernst, Ziv Bar-Joseph; BMC Systems Biology, 2012,
[full text](https://bmcsystbiol.biomedcentral.com/track/pdf/10.1186/1752-0509-6-104)


