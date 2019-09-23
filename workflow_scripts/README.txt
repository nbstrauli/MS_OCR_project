Contents of workflow_scripts:

diversity_linear_mixed_model.R - This is an R script that runs the linear mixed model for modeling the Shannon Index diversity data.

get_clonal_freq_spectra.py - This script processes the clonal sequence data to get the frequency spectra across clonal clusters for each sample.

get_gene_usage.py - This script processes the clonal sequence data to get the expression level (i.e. number of reads that map to a given gene, normalized by the total number of reads) for each gene segment (V or J) across each time-point.

regression_based_geneByGene_test.R - This is an R script that runs the linear regression for each gene segment.