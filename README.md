# ptbmethylation
Preterm Birth temporal DNA methylation data analysis pipeline.

This repository contains initial source codes to reproduce the temporal maternal peripheral blood DNA methylation data analysis performed in the paper by Das, et al. which is titled "Pattern of variation in DNA methylation during pregnancy among mothers who delivered preterm in the GARBH-Ini cohort".

## Table of Contents

* [R](https://www.r-project.org/)
* [Python](https://www.pitt.edu/~naraehan/python2/)

## R

For reading the IDAT files and necessary preprocessing and filtering, you would require the following packages in R:

* [ChAMP](https://www.bioconductor.org/packages/release/bioc/vignettes/ChAMP/inst/doc/ChAMP.html)
* [ENmix](https://www.bioconductor.org/packages/release/bioc/vignettes/ENmix/inst/doc/ENmix.pdf)
* [preprocessCore](https://bioconductor.riken.jp/packages/3.8/bioc/manuals/preprocessCore/man/preprocessCore.pdf)
* [dplyr](https://cran.r-project.org/web/packages/dplyr/vignettes/dplyr.html)
* [tibble](https://cran.r-project.org/web/packages/tibble/tibble.pdf)

The raw IDAT files for Infinium MethylationEPIC BeadChip (Illumina) and processed data is available in [GEO Submission GSE169338](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169338).

1. Load required packages on R

2. Read the IDAT files 
3. Remove SNPs, Multi Hit probes, Sex chromosomes, etc. from the data
4. Load rgSet data matrix for preprocessing by background correction, dye-bias correction, Quality control by identifying and excluding low quality probes
5. Load the pre-processed matrix for intra-array normalisation using BMIQ normalisation
6. Simultaneously perform batch effect correction (if necessary)
7. Load the batch effect corrected matrix for inter-array normalisation using Quantile Normalisation
8. Load the normalised data matrix for cell type normalisation, in case using a heterogeneous tissue sample for methylation analysis.

We used Illumina Infinium EPIC manifest file in which annotation is based on hg19 reference sequence. We considered probes in promoter specific regions of the genome, i.e. from 200 to 1500 bases upstream of transcription start sites and 5’ untranslated regions (TSS200, TSS1500, 5’UTR) for analysis.

And we want to remove the CpG sites with low beta values (intensity values), i.e. any value ≤0.2 in all samples at all time points. We then proceeded with the beta matrix for subsequent analysis.

## Python

For the analysis, I have used Python 2.7.18. Use the preprocessed, normalised and cell-type corrected beta matrix for this set of analysis. The goal here is to:

* Compute the temporal variance of beta values for each pregnant woman (each mother in our case were sampled through three time points through the period of gestation and at delivery), and find the mean of the temporal variances for each group (term / preterm)
* Compute the Standard deviation for each CpG site
* Compute the d-statistic
* Bootstrap (n=1000) d-statistic values with shuffled term and preterm status identity
* Find if the original d-statistic lies above 99th quantile of the 1000 bootstrapped d-statistic 

For performing the above analysis, you would require the following packages in python: sys, csv, random, math, statistics, and numpy.

1. Import necessary packages for the anaylysis (sys, csv, random, math, statistics, numpy)

2. Set input and output file parameters

3. Read the header for the input file

4. Draft output file header

5. For each row corresponding each CpG site in the input file, perform the following functions:

   a. Each CpG given as rownames for the input file

   b. Beta values clubbed in groups of 4, as for each mother, the data are arranged in 4 consecutive columns

   c. Calculate the overall variance for each CpG sites

   d, e. Calculate the mean of temporal variances and standard deviations for mothers who delivered at term

   f, g. Calculate the mean of temporal variances and standard deviations for mothers who delivered at preterm

   h. Find the standard deviation for each CpG site across all samples

   i. Compute the d-statistic, which is the difference in the mean of variances between preterm and term, divided by the pooled standard deviation for each CpG site

   j. Now, for each CpG site, we generate 1000 d-statistics by bootstrapping with shuffled term and preterm status identity

   k. We check now if the d-statistic that we obtained for each probe belong above the 99th quantile of the bootstrapped d-statistic

   l. If the d-statistic from for the CpG from the initial computation lies below the 99th quantile, we flag the CpG as "No", or else, we flag as "Yes"

   m. Write the computed values in the output file.

The CpG sites with largest d-statistic (top 5% of the list of CpG sites that were flagged  "Yes") were considered for downstream pathway analysis.