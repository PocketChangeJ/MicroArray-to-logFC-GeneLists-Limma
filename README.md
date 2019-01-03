# MicroArray-to-logFC-GeneLists-Limma
From gene expression microarray data to Limma and curation of resulting Gene Lists.

Input Files:

a. Gene Expression Data, tab-separated (colnames, rownames and data)

b. Classlabels, tab separated (0 or 1, last line empty)

c. Platform, tab separated (2 columns; probe ids and gene symbols)


Pre checking - do manually:

a. Set respective working directory (setwd("..."))

b. boxplot and decide on setting toNormalize = TRUE or FALSE


No console version due to manual checking that needs to be done during execution (i.e. check if normalized)
