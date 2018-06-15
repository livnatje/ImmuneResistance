# **Mapping immune resistance in cancer with single-cell data**

This resource provides the code developed in the study of Jerby-Arnon _et al._ **_"Single-cell RNA-seq of melanoma ecosystems reveals sources of T cell exclusion linked to immunotherapy clinical outcomes"_**. It reproduces the key results of the study and can be applied to other single-cell cohorts to explore cell-cell interactions in cancer.

## **Requirements**

* R (tested in R version 3.4.0 (2017-04-21) -- "You Stupid Darkness").
* R libraries: scde, matrixStats, plotrix, plyr, ppcor, survival, ROCR, Hmisc, rms, mixtools, lme4, lmerTest

## **Data**

The data is provided in the [Single Cell Portal](https://portals.broadinstitute.org/single_cell/study/melanoma-immunotherapy-resistance#study-download) (_**ImmRes_Rfiles.zip**_).

In the Portal you will also find the processed single-cell [gene expression](https://portals.broadinstitute.org/single_cell/study/melanoma-immunotherapy-resistance#study-download) along with [interactive views](https://portals.broadinstitute.org/single_cell/study/melanoma-immunotherapy-resistance#study-visualize).

## **Quick start**

To reproduce the results reported in Jerby-Arnon _et al._ download _**ImmRes_Rfiles.zip**_ from the [Single Cell Portal](https://portals.broadinstitute.org/single_cell/study/melanoma-immunotherapy-resistance#study-download). Unzip the file and move the resulting _Data_ directory to the ImmuneResistance directory. 

In _R_ go to the _Code_ directory and run ```master.code()``` which is provided in _ImmRes_master.R_. The ```master.code()``` will walk you through the different stages of the study, divided into six main modules:

[**(1-2)** First](https://github.com/livnatje/ImmuneResistance/wiki/Mapping-immune-resistance-in-melanoma), analyzing the single-cell data to generate various gene signatures that characterize different cell subtypes and immune resistant cell states. For more information see [_Mapping immune resistance in melanoma_](https://github.com/livnatje/ImmuneResistance/wiki/Mapping-immune-resistance-in-melanoma).

[**(3-5)** Next](https://github.com/livnatje/ImmuneResistance/wiki/Predicting-immunotherapy-resistance), analyzing independent cohorts obtained from bulk melanoma tumors to explore and test the immune resistance program. For more information see [_Predicting immunotherapy resistance_](https://github.com/livnatje/ImmuneResistance/wiki/Predicting-immunotherapy-resistance).

[**(6)** Lastly](https://github.com/livnatje/ImmuneResistance/wiki/Repressing-the-immune-resistance-program), performing a pan-cancer analysis to identify drugs that could repress the immune resistance program in cancer cells. For more information see [_Repressing the immune resistance program_](https://github.com/livnatje/ImmuneResistance/wiki/Repressing-the-immune-resistance-program).

# General notes

The code provided in ```ImmRes_master.R``` reproduces the key results of the study. It also generates the study figures and table in the ```Output``` directory. The code follows the analyses that were performed in the study in a sequential manner. 

As the results are already provided in the ```Results``` directory, it is possible to run only some parts of the code and focus on specific analyses, or [apply the approach to other datasets](https://github.com/livnatje/ImmuneResistance/wiki/Applying-the-approach-to-other-datasets).

## Citation

Jerby-Arnon L _et al._ _**Single-cell RNA-seq of melanoma ecosystems reveals sources of T cell exclusion linked to immunotherapy clinical outcomes**_.
