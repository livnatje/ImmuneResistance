# ImmuneResistance
Code and additional processed data for manuscript "Single-cell RNA-seq of melanoma ecosystems reveals sources of T cells exclusion linked to immunotherapy clinical outcomes"

The ImmRes_master.R file reproduces the key results of the study.

#### Generating de-novo cell subtype signatures

First it downloads the annotated clinical single-cell RNA-seq (scRNA-seq) data and analyzes it to derive cell subytpe specific signatures, by running the code which is provided in ```GitHub1_denovoCellTypeSig.R```

```R
print("1. Generating de-novo cell-type signatures.")
# The code is provided in "GitHub1_denovoCellTypeSig.R"
cell.type.de<-get.cell.type.sig()
```

The resulting signatures will be written to tables ```Output/Tables/TableS3B_denovo.cell.subtype.sig.csv``` (reproducing **Table S3B** from Jerby-Arnon et al. 2018).

It will then compute the overall expression of each cell subtype signature across the cells, and plot the results on a two dimentional embedding (t-Distributed Stochastic Neighbor Embedding (t-SNE)). The resulting figures will be found in  ```Output/Figures/Fig1D_tSNE.nonmal.pdf``` and ```Output/Figures/FigS1FG_tSNE.nonmal.pdf```(reproducing **Figures 1D** and **S1F-G** from Jerby-Arnon et al. 2018).

![Fig1](/Images/tSNE_nonmal.pdf)
