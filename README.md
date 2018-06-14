# ImmuneResistance
Code and additional processed data for manuscript "Single-cell RNA-seq of melanoma ecosystems reveals sources of T cells exclusion linked to immunotherapy clinical outcomes"

The ImmRes_master.R file reproduces the key results of the study.

#### 1. Generating de-novo cell subtype signatures

First it downloads the annotated clinical single-cell RNA-seq (scRNA-seq) data and analyzes it to derive cell subtype specific signatures, by running the code which is provided in ```GitHub1_denovoCellTypeSig.R```

```R
cell.type.de<-get.cell.type.sig()
```

The resulting signatures will be written to ```Output/Tables/TableS3B_denovo.cell.subtype.sig.csv``` (reproducing **Table S3B** from Jerby-Arnon et al. 2018).

It will then compute the overall expression of each cell subtype signature across the cells, and plot the results on a two dimensional embedding (t-Distributed Stochastic Neighbor Embedding (t-SNE)). The resulting figures will be found in  ```Output/Figures/Fig1D_tSNE.nonmal.pdf``` and ```Output/Figures/FigS1FG_tSNE.nonmal.pdf```(reproducing **Figures 1D** and **S1F-G** from Jerby-Arnon et al. 2018).

![tSNE_nonmal_small](/Images/tSNE_nonmal_small.png)


Interactive tSNE plots of the clinical and experimental single cell data are provided in the [Single Cell Portal](https://portals.broadinstitute.org/single_cell/study/melanoma-immunotherapy-resistance).

#### 2. Identifying immune resistance programs in malignant melanoma cells

Next, we will use the signatures that were generated for malignant and CD8 T cells (```cell.sig```) to characterize malignant cells in "cold" melanoma tumors (with low levels of T cell infiltration). Our data-driven approach will integrates single-cell (```r.sc```) and bulk RNA-Seq data (```rB```) to uncover the *T cell exclusion program* of malignant melanoma cells.

```R
exc.de<-mal.t.cell.exclusion(rB = r.tcga,r.sc = r.sc,cell.sig = cell.sig)
```
![Fig1A](/Images/Fig1A.png)

Next, we will compare malignant melanoma cells from post-immunotherapy (anti-PD1 and anti-CTLA4) resistant melanoma tumors to those from treatment naive tumors. We will identify differentially expressed genes and derive the the post-treatment signatures.

```R
trt.de<-get.post.trt.sig(r = r.sc,subs.trt = subs.trt)
```
We will then also identify genes which are co-regulated (positively) or anti-regulated (negatively) with genes whose inhibition desensitized melanoma cells to T cell mediated killing in functional screens [(Patel et al., 2017)](https://www.nature.com/articles/nature23477).

```R
fnc.de<-get.fnc.res.sig(r = r.sc)
```
Lastly, we will combine the different signatures into the immune resistance program.

```R
res.sig<-get.res.program()
```

The resulting signatures will be written to ```Output/Tables/TableS4A_resistance.program.csv``` (reproducing **Table S4B** from Jerby-Arnon et al. 2018).
