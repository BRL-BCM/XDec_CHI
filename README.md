# Sliding Window Extension - performing cell type specific correlations

This method uses proportions estimated via a deconvolution algorithm and a matrix of feature x samples to perform sliding window deconvolution. The output is either correlations between two specific genes in a cell type specific manner or all genes which correlate with a single gene of interest in a cell type specific manner.

## Using this package

### Installation
You can install using devtools: require(devtools) install_github(“BRL-BCM/XDec-CHI”).

Look at the package Rmd vignette.
Located in /vignette/introduction.Rmd. It will take you through both functions in the package:

Identifying correlations and plotting cell type specific correlations between two genes
Identifying all genes which correlate to a gene of interest in a cell type specific manner
