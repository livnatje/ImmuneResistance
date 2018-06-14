Testing and applying the immune resistance program to predict immunotherapy resistance in melanoma


### 3. Longitudinal analysis
Analyzing Validation Cohort 1.
The code is provided in ```ImmRes3_longitudinal.R```

```R
valCo1<-set.ValCo1()
valCo1.results<-test.ValCo1(r = valCo1)
```
Analyzing a MAPKi resistance cohort (Hugo et al. 2015).
```R
mapkiCo<-set.matched.MAPKi.Hugo()
mapki.results<-test.matched.MAPKi.Hugo(r = mapkiCo)
```
### 4. Predicting ICB responses in published cohorts
The code is provided in ```ImmRes4_predictICBresponses.R```
Use different signatures to predict ICB responses.
```R
r.tcga<-set.TCGA(r.tcga = r.tcga)
R<-set.public.ICB.cohorts()
```
Testing OS predictions (TCGA).
```R
tcga.OS.prf<-prd.TCGA.survival(r.tcga)
```
Testing ICB response predictions.
```R
  ICB.prf<-prd.public.ICB.response(R)
```
### 5. Predicting ICB responses in Validation Cohort 2.
The code is provided in ```ImmRes5_valCohort2.R```
5.1 Use different signatures to predict ICB responses.
```R  
r.pd1<-set.aPD1()
```
Testing ICB response predictions.
```R  
aPD1.val<-prd.aPD1()
generate.fig5() # Provided in "ImmRes_output.R"
```
### 6. Prioritizing drugs to repress the immune resistance program.") 
```R
find.sensitizing.drugs()
```
