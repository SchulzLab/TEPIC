We offer two machine learning approaches as poosible downstream applications for TEPIC gene-TF scores:

*INVOKE* is a based on linear regression to identify key regulators in a sample of interest. All genes are treated
as they are regulated similiary, thus TFs that affect the expression of many genes in the considered sample can be revealed. 

*DYNAMITE* uses logistic regression to identify TFs that can classify genes as being up or downregulated between tissues/samples. 
This can be used to find TFs that might be related to differences between two distinct cell types or cell states. 

For details on the methods, please check the individual subfolders as well as the [documentation](/docs/Description.pdf).
