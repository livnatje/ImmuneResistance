# ImmuneResistance
Code and additional processed data for manuscript "Single-cell RNA-seq of melanoma ecosystems reveals sources of T cells exclusion linked to immunotherapy clinical outcomes"

The ImmRes_master.R file reproduces the key results of the study.

#### 1. Generating de-novo cell subtype signatures

First it downloads the annotated clinical single-cell RNA-seq (scRNA-seq) data and analyzes it to derive cell subytpe specific signatures, by running the code which is provided in ```GitHub1_denovoCellTypeSig.R```

```R
cell.type.de<-get.cell.type.sig()
```

The resulting signatures will be written to tables ```Output/Tables/TableS3B_denovo.cell.subtype.sig.csv``` (reproducing **Table S3B** from Jerby-Arnon et al. 2018).

It will then compute the overall expression of each cell subtype signature across the cells, and plot the results on a two dimentional embedding (t-Distributed Stochastic Neighbor Embedding (t-SNE)). The resulting figures will be found in  ```Output/Figures/Fig1D_tSNE.nonmal.pdf``` and ```Output/Figures/FigS1FG_tSNE.nonmal.pdf```(reproducing **Figures 1D** and **S1F-G** from Jerby-Arnon et al. 2018).

![tSNE_nonmal_small](/Images/tSNE_nonmal_small.png)

Interactive tSNE plots of the clinical and experimental single cell data are provided in the [Single Cell Portal](https://portals.broadinstitute.org/single_cell/study/melanoma-immunotherapy-resistance).

#### 2. Identifying immune resistance programs in malignant melanoma cells

2.1. Generating the T cell exclusion signatures.
```R
exc.de<-mal.t.cell.exclusion(rB = r.tcga,r.sc = r.sc,cell.sig = cell.sig)
```
2.2. Generating the post-treatment signatures.
```R
trt.de<-get.post.trt.sig(r = r.sc,subs.trt = subs.trt)
```
2.3. Generating the functional signatures
```R
fnc.de<-get.fnc.res.sig(r = r.sc)
```
2.4. Combining the signatures into the immune resistance program.
```R
res.sig<-get.res.program()
```
