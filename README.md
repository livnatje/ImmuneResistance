# **Mapping immune resistance in cancer with single-cell data**

The purpose of this resource is to enable you to: **(1)** reproduce the key results of our recent study of immune resistance in melanoma (Jerby-Arnon _et al._); and **(2)** apply our approach to other single-cell cohorts to explore cell-cell interactions in cancer.

Visit our [wiki](https://github.com/livnatje/ImmuneResistance/wiki) page for instructions and details.

## **Requirements**

* R (tested in R version 3.4.0 (2017-04-21) -- "You Stupid Darkness").
* R libraries: scde, matrixStats, plotrix, plyr, ppcor, survival, ROCR, Hmisc, rms, mixtools, lme4, lmerTest

## **Data**

The data is provided through the [Single Cell Portal](https://portals.broadinstitute.org/single_cell/study/melanoma-immunotherapy-resistance#study-download) (_**ImmRes_Rfiles.zip**_).

In the Portal you will also find the processed single-cell [gene expression](https://portals.broadinstitute.org/single_cell/study/melanoma-immunotherapy-resistance#study-download) of the clinical cohort and experimental data, as well as [interactive views](https://portals.broadinstitute.org/single_cell/study/melanoma-immunotherapy-resistance#study-visualize).

## **Quick start**

To reproduce the results reported in Jerby-Arnon _et al._ download _**ImmRes_Rfiles.zip**_ from the [Single Cell Portal](https://portals.broadinstitute.org/single_cell/study/melanoma-immunotherapy-resistance#study-download). Unzip the file and move the resulting _Data_ directory to the ImmuneResistance directory. 

In _R_ go to the _Code_ directory and run ```master.code()``` which is provided in _ImmRes_master.R_.

The ```master.code()``` will walk you through the different stages of the study, divided into six main analyses:

**(1-2)** First, analyzing the single-cell data to generate various gene signatures that characterize different cell subtypes and immune resistant cell states. For more information see [I. Mapping immune resistance in melanoma](https://github.com/livnatje/ImmuneResistance/wiki/I.-Mapping-immune-resistance-in-melanoma)

**(3-5)** Next, analyzing independent cohorts obtained from bulk melanoma tumors to explore and test the immune resistance program. For more details see [II. Predicting immunotherapy resistance](https://github.com/livnatje/ImmuneResistance/wiki/II.-Predicting-immunotherapy-resistance).

**(6)** Lastly, performing a pan-cancer analysis to identify drugs that could repress the immune resistance program in cancer cells.

## Citation

Please use the following citation:

Jerby-Arnon L _et al._ Single-cell RNA-seq of melanoma ecosystems reveals sources of T cell exclusion linked to immunotherapy clinical outcomes.
