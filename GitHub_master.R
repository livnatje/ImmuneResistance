git.master.code<-function(){
  ######################################################################################
  print("1. Generating de-novo cell-type signatures.")
  # The code is provided in "GitHub1_denovoCellTypeSig.R"
  cell.type.de<-git.get.cell.type.sig()
  
  ######################################################################################
  print("2. Identifying immune resistance programs in malignant melanoma cells.")
  # The code is provided in "GitHub2_immuneResistanceProgram.R"
  r.sc<-readRDS("/Volumes/MerckIT/GitHub/Data/scData/Mel.malignant.rds")
  r.tcga<-readRDS("/Volumes/MerckIT/GitHub/Data/PublicData/TCGA_SKCM.rds")
  load("/Volumes/MerckIT/GitHub/Results/CellTypes/cell.type.sig.RData")
  load("/Volumes/MerckIT/GitHub/Results/Resistance/PostTreatment/post.treatment.subsample.de.RData")
  
  print("2.1. Generating the T cell exclusion signatures.")
  exc.de<-git.mal.t.cell.exclusion(rB = r.tcga,r.sc = r.sc,cell.sig = cell.sig)
  print("2.2. Generating the post-treatment signatures.")
  trt.de<-git.get.post.trt.sig(r = r.sc,subs.trt = subs.trt)
  print("2.3. Generating the functional signatures.")
  fnc.de<-git.get.fnc.res.sig(r = r.sc)
  print("2.4. Combining the signatures into the immune resistance program.")
  res.sig<-git.get.res.program()
  iciA.sigs<-save.all.ici.predictors()
  
  ######################################################################################
  print("3. Longitudinal analysis.")
  # The code is provided in "GitHub3_longitudinal.R"
  print("3.1 Analyzing Validation Cohort 1.")
  valCo1<-git.set.ValCo1()
  valCo1.results<-git.test.ValCo1(r = valCo1)
  print("3.2 Analyzing a MAPKi resistance cohort (Hugo et al. 2015).")
  mapkiCo<-git.set.matched.MAPKi.Hugo()
  mapki.results<-git.test.matched.MAPKi.Hugo(r = mapkiCo)
  
  ######################################################################################
  print("4. Predicting ICB responses in published cohorts.")
  # The code is provided in "GitHub4_predictICBresponses.R"
  print("4.1 Use different signatures to predict ICB responses.")
  r.tcga<-git.set.TCGA(r.tcga = r.tcga)
  R<-git.set.public.ICB.cohorts()
  print("4.2 Testing OS predictions (TCGA).")
  tcga.OS.prf<-git.prd.TCGA.survival(r.tcga)
  print("4.3 Testing ICB response predictions.")
  ICB.prf<-git.prd.public.ICB.response(R)
  
  ######################################################################################
  print("5. Predicting ICB responses in Validation Cohort #2.")
  # The code is provided in "GitHub5_valCohort2.R"
  print("5.1 Use different signatures to predict ICB responses.")
  r.pd1<-git.set.aPD1()
  print("5.2 Testing ICB response predictions.")
  aPD1.val<-git.prd.aPD1()
  
  ######################################################################################
  print("6. Prioritizing drugs to repress the immune resistance program.")
  git.find.sensitizing.drugs()
  
  ######################################################################################
  print("7. Exploring the impact of CDK4/6 inhibition on melanoma cells.")
  
  
  return()
}

git.master.extra<-function(){
  ##########################************** EXTRA **************#########################
  
  
  ######################################################################################
  print("Co-regulation of the resistance program within and across tumors")
  # Code provided in Melanoma.functions.R
  ICR.coregulation()
  
  print("c-map analysis.")
  # Code provided in Melanoma.functions.R
  # Write the shorter version of the resistance program
  cmap.analysis()
  # Run c-map (from the website)
  # Analyze the results
  cmap.analysis()
  
  print("Comparison to other signatures.")
  
  ######################################################################################
  print("Spatial analyses.")
  # Code provided in GitHub_spatial.R
  spatial.master()
  
  ######################################################################################
  print("Intrinsic variation analyses (scRNA-seq of mice cancer cell lines).")
  # Code provided in Melanoma.functions.R
  # Load the mouse cell lines scRNA-seq data; identify DEG and compare to the resistance program.
  mouse.cl.rsl<-mouse.cell.lines.analysis()
  
}


