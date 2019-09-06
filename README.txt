The files in this repository are the code that formed the statistical analyses for the manuscript titled: "Effect of Ocrelizumab on B and T-cell immune repertoires in patients with relapsing MS". The code herein is provided to serve as a reference, so that researchers can see exactly what was done. As such, if imported, much of these scripts will not run for the following reasons: they have many dependencies, in some instances were designed to run on a computational cluster (using SGE), and may even have absolute paths hardcoded in them.

The code is divided into two classes: 

"workflow_scripts" contains the scripts that execute a given analysis. They will do simple things like parsing filepaths, and walk a dataset through the different steps of an analysis. They will typically invoke a given python class, and execute the methods of that class. 

"python_classes" contains the bulk of the code. It is as the name suggest, a directory that contains the python code for a variety of classes that were designed for different data types. For example, there is a class for immune time series data ("immune_time_series_class.py") which contains a variety of methods useful for immune repertoire sequence data that is of a longitudinal nature. 


Contents of workflow_scripts:

diversity_linear_mixed_model.R - This is an R script that runs the linear mixed model for modeling the Shannon Index diversity data.

get_clonal_freq_spectra.py - This script processes the clonal sequence data to get the frequency spectra across clonal clusters for each sample.

get_gene_usage.py - This script processes the clonal sequence data to get the expression level (i.e. number of reads that map to a given gene, normalized by the total number of reads) for each gene segment (V or J) across each time-point.

regression_based_geneByGene_test.R - This is an R script that runs the linear regression for each gene segment. 
