git.master.code<-function(prf.stages){
  if(is.element(1,prf.stages)){
    print("1. Generating de-novo cell-type signatures..")
    rF<-readRDS("/Volumes/MerckIT/GitHub/Data/scData/Mel.all.data.QC.rds")
    cell.type.de<-git.get.cell.type.sig(rF)
    cell.sig<-cell.type.de$sig
  }else{
    print("1. Loading de-novo cell-type signatures..")
    load("/Volumes/MerckIT/GitHub/Results/Signatures/cell.type.sig.RData") # load cell.sig
  }
  
  if(is.element(2,prf.stages)){
    print("2. Generating the T cell exclusion signatures.")
    r.mal<-readRDS("/Volumes/MerckIT/GitHub/Data/scData/Mel.malignant.rds")
    r.tcga<-readRDS("/Volumes/MerckIT/GitHub/Data/PublicData/TCGA_SKCM.rds")
    exc.de<-git.mal.t.cell.exclusion(r.tcga,r.mal,cell.sig)
  }
  
  if(is.element(3,prf.stages)){
    print("3. Generating the post-treatment signatures.")
    if(!exists("r.mal")){r.mal<-readRDS("/Volumes/MerckIT/GitHub/Data/scData/Mel.malignant.rds")}
    load("/Volumes/MerckIT/GitHub/Results/PostTreatment/post.treatment.subsample.de.RData")
    trt.de<-git.get.post.trt.sig(r.mal,subs.trt)
  }
  
  if(is.element(4,prf.stages)){
    if(!exists("r.mal")){r.mal<-readRDS("/Volumes/MerckIT/GitHub/Data/scData/Mel.malignant.rds")}
    print("4. Generating the functional signatures.")
    fnc.de<-git.get.fnc.res.sig(r.mal)
  }
  
  if(is.element(5,prf.stages)){
    print("5. Combining the signatures into the immune resistance program.")
    load("/Volumes/MerckIT/GitHub/Results/Signatures/Tcell.exclusion.sig.RData")
    load("/Volumes/MerckIT/GitHub/Results/Signatures/post.treatment.sig.RData")
    load("/Volumes/MerckIT/GitHub/Results/Signatures/functional.resistance.sig.RData")
    res.sig<-git.get.res.program(exc.sig = exc.sig,trt.sig = trt.sig,fnc.sig = fnc.sig)
    iciA.sigs<-save.all.ici.predictors()
  }else{
    load("/Volumes/MerckIT/GitHub/Results/Signatures/resistance.program.RData")
  }
  if(is.element(6,prf.stages)){
    print("6. Predicting ICB responses in external and new cohorts.")
    r.tcga<-git.set.TCGA()
    R<-git.set.public.ICB.cohorts()
    r.pd1<-git.set.aPD1()
  }
  if(is.element(7,prf.stages)){
    print("7. Testing ICB response predictions.")
    ICB.prf<-git.prd.public.ICB.response(R)
    aPD1.val<-git.prd.aPD1()
  }
  return()
}

git.get.normalized.gene.counts<-function(cell.type){
  fileName<-paste0("/Volumes/MerckIT/GitHub/Data/scData/Mel.",cell.type,".QC.rds")
  r<-readRDS(fileName)
  pagodaR<-pagoda.get.denovo.modules(r$cd,sampleName = r$data.name,work.dir = "MerckIT")
  r$varinfo<-pagodaR$varinfo
  r$knn<-pagodaR$knn
  r$norm<-r$varinfo$mat
  saveRDS(r,file = fileName)
}

git.get.cell.type.sig<-function(r,subtype.flag = F){
  if(is.null(r)){
    r<-readRDS("/Volumes/MerckIT/GitHub/Data/scData/Mel.all.data.QC.rds")
    cell.type.de<-git.get.cell.type.sig(r)
    return(cell.type.de)
  }
  cell.types.u<-sort(unique(r$cell.types))
  # Make all pairwise comparisons
  cell.type.pairs <- t(combn(unique(r$cell.types),2))
  de.zscores<-apply(cell.type.pairs,1,function(x){
    print(paste("Comparing",x[1],"to",x[2]))
    b1<-is.element(r$cell.types,x[1])
    b2<-is.element(r$cell.types,x[2])
    de<-git.apply.ttest(r$tpm[,b1|b2],b1[b1|b2],ranksum.flag = T)[,"zscores"]
    return(de)
  })
  rownames(de.zscores)<-r$genes
  colnames(de.zscores)<-paste(cell.type.pairs[,1],cell.type.pairs[,2],sep = "_")
  results<-list(cell.type.pairs = cell.type.pairs,cell.types = cell.types.u,de.zscores = de.zscores)
  results$de.sum<-lapply(cell.types.u, function(x){
    b1<-is.element(cell.type.pairs[,1],x)
    b2<-is.element(cell.type.pairs[,2],x)
    de<-cbind.data.frame(de.zscores[,b1],-de.zscores[,b2])
    colnames(de)<-c(cell.type.pairs[b1,2],cell.type.pairs[b2,1])
    de$Min<-rowMin(de)
    return(de)
  })
  names(results$de.sum)<-cell.types.u
  results$sig<-lapply(results$de.sum,function(x) sort(rownames(x)[x$Min>5]))
  print(summary(results$sig))
  if(!subtype.flag){results<-git.get.t.cell.sig(results)}
  
  results$sig<-remove.ribo(results$sig)
  
  cell.type.de<-results
  cell.sig<-results$sig
  if(!subtype.flag){
    save(cell.type.de,file = "/Volumes/MerckIT/GitHub/Results/CellTypes/cell.type.de.RData")
    save(cell.sig,file = "/Volumes/MerckIT/GitHub/Results/CellTypes/cell.type.sig.RData")
    save(cell.sig,file = "/Volumes/MerckIT/GitHub/Results/Signatures/cell.type.sig.RData")
  }
  return(cell.type.de)
}

git.get.t.cell.sig<-function(results){
  f<-function(m,cell.types){
    cell.types<-intersect(colnames(m),cell.types)
    m$Min_noT<-rowMin(m[,setdiff(cell.types,c("T.CD4","T.CD8","NK"))])
    print(colnames(m))
    m$Min_noNK<-rowMin(m[,setdiff(cell.types,"NK")])
    return(m)
  }
  results$de.sum$T.CD4<-f(results$de.sum$T.CD4,results$cell.types)
  results$de.sum$T.CD8<-f(results$de.sum$T.CD8,results$cell.types)
  genes<-rownames(results$de.zscores)
  results$sig$T.cell<-sort(genes[results$de.sum$T.CD4$Min_noT>5&results$de.sum$T.CD8$Min_noT>5])
  results$sig$T.CD8.NK<-sort(genes[results$de.sum$T.CD8$Min_noNK>5])
  results$sig$T.CD8<-get.top.elements(-results$de.sum$T.CD8,q = 50)[["Min"]]
  results$sig$T.CD4<-get.top.elements(-results$de.sum$T.CD4,q = 50)[["Min"]]
  print(summary(results$sig))
  return(results)
}

git.get.cell.supertype.sig<-function(de){
  f<-function(de,cell.types,supertype){
    other.types<-setdiff(de$cell.types,cell.types)
    m<-get.mat(rownames(de$de.zscores),m.cols = cell.types)
    for(i in cell.types){
      m[,i]<-rowMin(de$de.sum[[i]][,other.types])
    }
    m<-cbind.data.frame(m,min = rowMin(m))
    de$de.sum[[supertype]]<-m
    de$sig[[supertype]]<-get.top.elements(-m,100,min.ci = (-3))$min
    return(de)
  }
  supertypes<-list(stroma = c('CAF','Endo.'),
                   lymphocyte = c('B.cell','T.CD4','T.CD8'),
                   immune = c('B.cell','Macrophage','NK','T.CD4','T.CD8'))
  for(x in names(supertypes)){
    de<-f(de = de,cell.types = supertypes[[x]],supertype = x)
  }
  print(summary(de$sig))
  cell.type.de<-de
  cell.sig<-de$sig
  save(cell.type.de,file = "/Volumes/MerckIT/GitHub/Results/CellTypes/cell.type.de.full.RData")
  save(cell.sig,file = "/Volumes/MerckIT/GitHub/Results/CellTypes/cell.type.sig.full.RData")
  save(cell.sig,file = "/Volumes/MerckIT/GitHub/Results/Signatures/cell.type.sig.full.RData")
  return(de)
}

git.mal.t.cell.exclusion<-function(rB,r.sc,cell.sig){
  if(is.null(rB)){
    r.sc<-readRDS("/Volumes/MerckIT/GitHub/Data/scData/Mel.malignant.rds")
    rB<-readRDS("/Volumes/MerckIT/GitHub/Data/PublicData/TCGA_SKCM.rds")
    load("/Volumes/MerckIT/GitHub/Results/CellTypes/cell.type.sig.RData")
    results<-git.mal.t.cell.exclusion(rB,r.sc,cell.sig)
    return(results)
    
    load("/Volumes/MerckIT/Data/PublicResources/Harbst.mel.bulk.RData");r.har<-r;rm(r)
    results<-git.cell.cell.intr(r1,r.sc,cellA.markers = cell.sig$Mal,
                                cellB.markers = cell.sig$T.CD8,
                                cellA.name = "malignant",cellB.name = "T.CD8",
                                fileName = NULL, pval = 0.1)
    
  }
  
  results<-git.cell.cell.intr(rB,r.sc,cellA.markers = cell.sig$Mal,
                              cellB.markers = cell.sig$T.CD8,
                              cellA.name = "malignant",cellB.name = "T.CD8",
                              fileName = "FINAL_Tcell_Exclusion", pval = 0.1)
  exc.sig<-results$sig[c("exc","att")]
  exc.sig<-lapply(exc.sig,function(x) sort(unique(x)))
  names(exc.sig)<-c("exc.up","exc.down")
  print(summary(exc.sig))
  
  save(exc.sig,file = "/Volumes/MerckIT/GitHub/Results/Resistance/Exclusion/FINAL_Tcell_Exclusion.sig.RData")
  save(exc.sig,file = "/Volumes/MerckIT/GitHub/Results/Signatures/Tcell.exclusion.sig.RData")
  
  results<-git.test.seed.sig(r.sc,results,no.itr = 1000)
  save(results,file = "/Volumes/MerckIT/GitHub/Results/Resistance/Exclusion/FINAL_Tcell_Exclusion.full.RData")
  
  return(results)
}

git.cell.cell.intr.full<-function(rB,r.sc,cellA.markers,cellB.markers,
                                  cellA.name = "malignant",cellB.name = "T.cell",
                                  bulk.confounders = NULL,sc.confounders = NULL,
                                  fileName = NULL, pval = 0.1,no.itr = 50){
  B<-get.mat(r.sc$genes,1:no.itr)
  B.seed.exc<-B;B.seed.att<-B
  B.exc<-B;B.att<-B;
  B.excA<-B;B.attA<-B
  for(i in 1:no.itr){
    print(paste("Iteration no.",i))
    results<-git.cell.cell.intr(rB,r.sc,cellA.markers = cellA.markers,
                                cellB.markers = cellB.markers,
                                cellA.name = cellA.name,cellB.name = cellB.name,
                                bulk.confounders = bulk.confounders,
                                sc.confounders = sc.confounders,
                                fileName = NULL, pval = pval)
    B.seed.exc[,i]<-is.element(r.sc$genes,results$seed.sig$exc)
    B.seed.att[,i]<-is.element(r.sc$genes,results$seed.sig$att)
    B.exc[,i]<-is.element(r.sc$genes,results$sig$exc)
    B.att[,i]<-is.element(r.sc$genes,results$sig$att)
    B.excA[,i]<-is.element(r.sc$genes,results$sigA$exc)
    B.attA[,i]<-is.element(r.sc$genes,results$sigA$att)
    print(table(rowSums(B.seed.exc,na.rm = T)))
    print(table(rowSums(B.seed.att,na.rm = T)))
    print(table(rowSums(B.exc,na.rm = T)))
    print(table(rowSums(B.att,na.rm = T)))
    print(table(rowSums(B.excA,na.rm = T)))
    print(table(rowSums(B.attA,na.rm = T)))
  }
  results$seed.sig$exc.final<-r.sc$genes[rowSums(B.seed.exc,na.rm = T)==no.itr]
  results$seed.sig$att.final<-r.sc$genes[rowSums(B.seed.att,na.rm = T)==no.itr]
  results$sig$exc.final<-r.sc$genes[rowSums(B.exc,na.rm = T)==no.itr]
  results$sig$att.final<-r.sc$genes[rowSums(B.att,na.rm = T)==no.itr]
  results$sig$excA.final<-r.sc$genes[rowSums(B.excA,na.rm = T)==no.itr]
  results$sig$attA.final<-r.sc$genes[rowSums(B.attA,na.rm = T)==no.itr]
  results$B.exc<-B.exc;results$B.att<-B.att
  results$B.excA<-B.excA;results$B.attA<-B.attA
  print(summary(results$sig))
  if(!is.null(fileName)){
    save(results,file = paste0("/Volumes/MerckIT/GitHub/Results/Resistance/Exclusion/",fileName,".full.RData"))
  }
  return(results)
}

git.cell.cell.intr<-function(rB,r.sc,cellA.markers,cellB.markers,cellA.name = "malignant",
                             cellB.name = "T.cell",bulk.confounders = NULL,sc.confounders = NULL,
                             fileName = NULL, pval = 0.1){
  print(paste("Characterizing",cellA.name,"cells in tumors with low",cellB.name,"infiltration."))
  # print("Setting the seed");set.seed(1234)
  cellA.markers<-sort(intersect(rB$genes,cellA.markers))
  results<-list(cellA.markers = cellA.markers,
                cellB.markers = cellB.markers,
                bulk.confounders = bulk.confounders,
                sc.confounders = sc.confounders)
  
  print("1. Estimating cell type B abundance in bulk gene expression")
  results$bulk.cellB.abn <- round(get.sign.scores.bulk(rB,gene.sign = list(cellB.markers),num.rounds = 1000),2)
  rownames(results$bulk.cellB.abn)<-rB$samples
  
  print("2. Looking for genes correlated with cell type B abundance in bulk gene expression")
  if(is.null(bulk.confounders)){
    results$bulk.cor<-spearman.cor(t(rB$tpm),results$bulk.cellB.abn,method = "pearson")
  }else{
    bulk.confounders<-as.matrix(bulk.confounders)
    b<-!is.na(rowSums(bulk.confounders))
    print(paste("Using",sum(b),"bulk samples (with full data)."))
    bulk.confounders<-bulk.confounders[b,]
    results$bulk.cor<-pcor.mat(t(rB$tpm[,b]),results$bulk.cellB.abn[b],bulk.confounders,method = "pearson")
  }
  results$bulk.cor<-git.add.onesided.p(results$bulk.cor)
  
  print("3. Getting the seed signatures")
  results$seed.sig<-get.top.elements(results$bulk.cor[cellA.markers,c("p.pos","p.neg")],q = 20,min.ci = pval)
  names(results$seed.sig)<-c("att","exc")
  print(summary(results$seed.sig))
  
  print("4. Testing if the seed signatures are anti correlated")
  results$sc.seed.scores<-round(get.sign.scores(r.sc,results$seed.sig,num.rounds = 1000),2)
  cor.plot(results$sc.seed.scores,main = "Seed overall expression")
  
  print("5. Expanding the seed signatures")
  r.sc$q<-cbind(r.sc$comp,log(r.sc$comp.reads))
  if(!is.null(sc.confounders)){
    r.sc$q<-cbind(r.sc$q,sc.confounders)
    results$sc.seed.scores.rgr<-t(get.residuals(t(results$sc.seed.scores),sc.confounders))
    cor.plot(results$sc.seed.scores.rgr,main = "Seed residuals")
  }
  
  results$att.sc<-pcor.mat(t(r.sc$tpm),results$sc.seed.scores[,1],r.sc$q,method = "spearman")#"pearson")
  results$exc.sc<-pcor.mat(t(r.sc$tpm),results$sc.seed.scores[,2],r.sc$q,method = "spearman")#"pearson")
  results$att.sc<-git.add.onesided.p(results$att.sc)
  results$exc.sc<-git.add.onesided.p(results$exc.sc)
  
  print("6. Generating the final signatures")
  f<-function(s){
    names(s)<-gsub("p.","",names(s),fixed = T)
    s$exc<-intersect(s$att.neg,s$exc.pos)
    s$att<-intersect(s$att.pos,s$exc.neg)
    return(s)
  }
  b.rp<-!startsWith(r.sc$genes,"RP")
  s1<-c(get.top.elements(m = results$att.sc[,c("p.pos","p.neg")],q = 200,min.ci = 0.01,main = "att"),
        get.top.elements(m = results$exc.sc[,c("p.pos","p.neg")],q = 200,min.ci = 0.01,main = "exc"))
  s2<-c(get.top.elements(m = results$att.sc[b.rp,c("p.pos","p.neg")],q = 200,min.ci = 0.01,main = "att"),
        get.top.elements(m = results$exc.sc[b.rp,c("p.pos","p.neg")],q = 200,min.ci = 0.01,main = "exc"))
  s1<-f(s1);s2<-f(s2)
  results$sigA<-s1
  results$sig.no.rp<-s1
  results$sig<-lapply(names(s1),function(x) sort(unique(c(s1[[x]],s2[[x]]))))
  names(results$sig)<-names(s1)
  print(summary(results$sig[c("exc","att")]))
  if(!is.null(fileName)){
    save(results,file = paste0("/Volumes/MerckIT/GitHub/Results/Resistance/Exclusion/",fileName,".RData"))
  }
  return(results)
}

git.test.seed.sig<-function(r,results,no.itr = 1000){
  no.cellA.markers <- length(results$cellA.markers)
  seed.cor<-spearman.cor(results$sc.seed.scores,method = "pearson")
  rnd.seed.att<-list();rnd.seed.exc<-list()
  B.rnd.cellA.markers <- get.compatible.non.overlapping.sig(r = r,sig.genes = results$cellA.markers,num.rounds = no.itr)
  for(i in 1:no.itr){
    rnd.cellA.markers <- r$genes[B.rnd.cellA.markers[,i]]
    rnd.seed<-get.top.elements(results$bulk.cor[rnd.cellA.markers,c("p.pos","p.neg")],q = 20,min.ci = 0.1)
    rnd.seed.att[[i]]<-rnd.seed$p.pos
    rnd.seed.exc[[i]]<-rnd.seed$p.neg
  }
  scores.att<-round(get.sign.scores(r,rnd.seed.att,num.rounds = 100),2);scores.att[,1]<-results$sc.seed.scores[,"att"]
  scores.exc<-round(get.sign.scores(r,rnd.seed.exc,num.rounds = 100),2);scores.exc[,1]<-results$sc.seed.scores[,"exc"]
  results$test.seed<-list()
  results$test.seed$rand.seed.cor<-spearman.cor(scores.att,scores.exc,method = "pearson",match.flag = T)
  hist(results$test.seed$rand.seed.cor[,"R"],100);abline(v = results$test.seed$rand.seed.cor[1,"R"],col = "red")
  print(get.empirical.p.value(results$test.seed$rand.seed.cor[,"R"]))
  
  X<-10*((2^r$tpm)-1)
  r$genes.dist<-log2(rowMeans(X,na.rm = T)+1)
  r$genes.dist.q<-my.discretize(r$genes.dist,n.cat = 50)
  r$zscores<-center.matrix(r$tpm,dim = 1,sd.flag = F)
  raw.scores.att<-get.random.sig.scores(r = r,genes.dist.q = r$genes.dist.q,
                                        b.sign = is.element(r$genes,results$seed.sig$att),
                                        num.rounds = no.itr,full.flag = T)
  raw.scores.exc<-get.random.sig.scores(r = r,genes.dist.q = r$genes.dist.q,
                                        b.sign = is.element(r$genes,results$seed.sig$exc),
                                        num.rounds = no.itr,full.flag = T)
  scores.att<-t(laply(1:no.itr,function(i){raw.scores.att[,i] - rowMeans(raw.scores.att[,-i])}))
  scores.exc<-t(laply(1:no.itr,function(i){raw.scores.exc[,i] - rowMeans(raw.scores.exc[,-i])}))
  
  scores.att[,1]<-results$sc.seed.scores[,"att"]
  scores.exc[,1]<-results$sc.seed.scores[,"exc"]
  
  scores.att<-round(scores.att,2)
  scores.exc<-round(scores.exc,2)
  
  results$test.seed$rand.sig.cor<-spearman.cor(scores.att,scores.exc,method = "pearson",match.flag = T)
  hist(results$test.seed$rand.sig.cor[,"R"],100);abline(v = results$test.seed$rand.sig.cor[1,"R"],col = "red")
  results$test.seed$emp.p<-rbind(rand.seed = get.empirical.p.value(results$test.seed$rand.seed.cor[,"R"]),
                                 rand.sig = get.empirical.p.value(results$test.seed$rand.sig.cor[,"R"]))
  print(results$test.seed$emp.p)
  return(results)
}

git.get.post.trt.sig<-function(r,subs.trt){
  if(is.null(r)){
    load("/Volumes/MerckIT/GitHub/Results/Resistance/PostTreatment/post.treatment.subsample.de.RData")
    r<-readRDS("/Volumes/MerckIT/GitHub/Data/scData/Mel.malignant.rds")
    trt.de<-git.get.post.trt.sig(r,subs.trt)
    return(trt.de)
  }
  trt.de<-git.find.DEG(r,main = "post.treatment",mainP = "treated")
  genes<-rownames(trt.subsmp)
  sig.subsmp<-list(up = genes[trt.subsmp[,"F.up"]>quantile(trt.subsmp[,"F.up"],0.9,na.rm = T)],
                   down = genes[trt.subsmp[,"F.down"]>quantile(trt.subsmp[,"F.down"],0.9,na.rm = T)])
  trt.de$sig<-list(trt.up = intersect(sig.subsmp$up,trt.de$sig$up),
                   trt.down = intersect(sig.subsmp$down,trt.de$sig$down))
  trt.sig<-trt.de$sig
  print(summary(trt.sig))
  save(trt.de,file = "/Volumes/MerckIT/GitHub/Results/Resistance/PostTreatment/post.treatment.de.RData")
  save(trt.sig,file = "/Volumes/MerckIT/GitHub/Results/Resistance/PostTreatment/post.treatment.sig.RData")
  save(trt.sig,file = "/Volumes/MerckIT/GitHub/Results/Signatures/post.treatment.sig.RData")
  return(trt.de)
}

git.get.fnc.res.sig<-function(r){
  if(is.null(r)){
    r<-readRDS("/Volumes/MerckIT/GitHub/Data/scData/Mel.malignant.rds")
    fnc.res.de<-git.get.fnc.res.sig(r)
    return(fnc.res.de)
  }
  res.down.seed<-c('B2M','CD58','HLA-A','MLANA','SOX10','SRP54','TAF3','TAP2','TAPBP')
  res.down.seed<-intersect(res.down.seed,r$genes)
  print("Characterzing the transcriptional state of malignant cells that underexpress one of these genes:")
  print(res.down.seed)
  
  B<-apply(r$norm[res.down.seed,],1,function(x) x<quantile(x,0.01));print(table(rowSums(B)))
  r$res.cells<-rowSums(B)>0
  print(paste("Found",sum(r$res.cells),"resistant cells."))
  fnc.res.de<-git.find.DEG(r,main = "functional.resistance",mainP = "res.cells",q = 250)
  # fnc.res.de$sig<-git.get.scde.ttest.sig(fnc.res.de,ttest.flag = F,rank.flag = T,q = 250)
  names(fnc.res.de$sig)<-c("fnc.up","fnc.down")
  fnc.sig<-fnc.res.de$sig
  save(fnc.res.de,file = "/Volumes/MerckIT/GitHub/Results/Resistance/Functional/functional.resistance.de.RData")
  save(fnc.sig,file = "/Volumes/MerckIT/GitHub/Results/Resistance/Functional/functional.resistance.sig.RData")
  save(fnc.sig,file = "/Volumes/MerckIT/GitHub/Results/Signatures/functional.resistance.sig.RData")
  return(fnc.res.de)
}

git.get.res.program<-function(exc.sig,trt.sig,fnc.sig){
  if(is.null(exc.sig)){
    load("/Volumes/MerckIT/GitHub/Results/Signatures/Tcell.exclusion.sig.RData")
    load("/Volumes/MerckIT/GitHub/Results/Signatures/post.treatment.sig.RData")
    load("/Volumes/MerckIT/GitHub/Results/Signatures/functional.resistance.sig.RData")
  }
  res.sig<-c(exc.sig,trt.sig,fnc.sig,
             list(resi.up = intersect(trt.sig$trt.up,exc.sig$exc.up),
                  resi.down = intersect(trt.sig$trt.down,exc.sig$exc.down),
                  res.up = c(trt.sig$trt.up,exc.sig$exc.up),
                  res.down = c(trt.sig$trt.down,exc.sig$exc.down)))
  res.sig$res.final.up<-intersect(res.sig$res.up,fnc.sig$fnc.up)
  res.sig$res.final.down<-intersect(res.sig$res.down,fnc.sig$fnc.down)
  res.sig<-lapply(res.sig, function(x) sort(unique(x)))
  print(summary(res.sig))
  save(res.sig,file = "/Volumes/MerckIT/GitHub/Results/Signatures/resistance.program.RData")
  
  load("/Volumes/MerckIT/Data/PublicResources/CRISPR_screen/mouse_human_mapping.RData")
  res.sig.mouse<-lapply(res.sig, function(x) convert.to.mice(mouse.human = mouse.human,x))
  saveRDS(res.sig.mouse,"/Volumes/MerckIT/GitHub/Results/Signatures/resistance.program.mouse.rds")
  
  return(res.sig)
}

git.find.DEG<-function(r = NULL,main = "",mainP = "treated",batch = NULL,q = 200){
  if(!is.null(batch)){batch<-as.factor(batch)}
  print(paste("Running DE:",main))
  groups<-r[[mainP]]
  print(paste0(mainP,": Positive ",sum(groups),", negative ",sum(!groups)))
  de<-list(main = main)
  de$ttest<-git.apply.ttest(r$tpm,groups,ranksum.flag = F)
  de$scde <- scde.expression.difference(r$knn, r$cd, r$varinfo$prior, groups = factor(!groups), 
                                        n.randomizations  =  100, n.cores  =  1, verbose  =  1,
                                        batch = batch)
  if(!is.null(batch)){
    de$scdeF<-de$scde;de$scde<-de$scdeF$batch.adjusted
  }
  b<-rowSums(r$norm,na.rm = T)==0;r$norm[b,]<-r$tpm[b,]
  de$ranksum<-git.apply.ttest(r$norm,groups,ranksum.flag = T)
  
  cor.plot(de$scde[,"Z"],de$ttest[,"zscores"])
  de$sum<-cbind.data.frame(scde = de$scde[,c("Z","cZ")],
                           ttest = de$ttest[,c("BH.more","BH.less")],
                           ranksum = de$ranksum[,c("BH.more","BH.less")])
  print(spearman.cor(de$sum)$cor)
  de$sig<-git.get.scde.ttest.sig(de,q = q,scde.c = (-1.96),
                                 ttest.c = 0.01,rank.c = 0.01,
                                 ttest.flag = F,rank.flag = T)
  return(de)
}

git.get.scde.ttest.sig<-function(de,q = 200,scde.c = (-1.96),
                                 ttest.c = 0.01,rank.c = 0.01,
                                 ttest.flag = F,rank.flag = T){
  sig.scde<-get.top.cor(de$sum,q = q,min.ci = scde.c)[c("scde.cZ.up","scde.cZ.down")]
  sig<-sig.scde
  if(ttest.flag){
    sig.ttest<-get.top.elements(de$sum[,c("ttest.BH.more","ttest.BH.less")],
                                q = q,min.ci = ttest.c)
    sig<-sig.ttest
  }
  if(rank.flag){
    sig.rank<-get.top.elements(de$sum[,c("ranksum.BH.more","ranksum.BH.less")],
                               q = q,min.ci = rank.c)
    sig<-sig.rank
  }
  if(ttest.flag&rank.flag){
    sig<-list(up = intersect(sig.ttest$ttest.BH.more,sig.rank$ranksum.BH.more),
              down = intersect(sig.ttest$ttest.BH.less,sig.rank$ranksum.BH.less))
  }
  print(summary(sig))
  sig<-list(up = intersect(sig.scde$scde.cZ.up,sig[[1]]),
            down = intersect(sig.scde$scde.cZ.down,sig[[2]]));print(summary(sig))
  return(sig)
}

git.plot.cell.sig<-function(){
  load("/Volumes/MerckIT/GitHub/Data/Mel.all.data.QC.RData")
  load("/Volumes/MerckIT/GitHub/Results/cell.type.de_tpm.ttest.RData")
  load("/Volumes/MerckIT/Results/Prognostic/all.ici.sigs.RData")
  source("source.heatmap.R")
  
  idx<-order(r$cell.types)
  r<-get.expZ(r)
  
  cell.sig<-setdiff.lists(cell.type.de$sig,cell.type.de.ttest$sig)
  summary(cell.sig)
  l<-names(cell.sig)
  names(iciA.sigs)<-gsub(" ",".",names(iciA.sigs))
  l<-names(iciA.sigs)[29:30]
  
  toMatch<-c("Fluidgm","Ayers","Riaz","Hugo","mapki")
  l<-unique(grep(paste(toMatch,collapse="|"), names(iciA.sigs),value=TRUE))
  b<-rowSums(r$zscoresC>0)>0
  
  cell.types<-r$cell.types
  cell.types[r$cell.types=="Mal"]<-"1.Mal"
  b<-is.element(r$cell.types,c("CAF","Endo"))
  cell.types[b]<-paste0("2.",cell.types[b])
  table(cell.types)
  
  plot.list<-lapply(l, function(x){
    my.violin(r$scores[,x],cell.types,order.flag = F,main = x)
  })
  
  lapply(l, function(x){
    pdf(paste0("~/Desktop/CellRevisions/Figures/Heatmap_",x,"_sig_rank.pdf"))
    genes<-intersect(r$genes[b],iciA.sigs[[x]])
    my.heatmap.function(r$zscoresC[genes,idx],col.labels = r$cell.types[idx],
                        row.labels = startsWith(genes,"RP"),
                        cluster.flag = "both",cor.cluster = T,main = x)
    dev.off()
  })
  
}

git.prd.MAPKi.Hugo<-function(r,res.sig,cell.sig,mapki.sig,mal.genes,all.genes){
  if(is.null(r)){
    library(lme4);library(lmerTest);attach(mtcars)
    r<-readRDS("/Volumes/MerckIT/GitHub/Data/PublicData/MAPKi.Hugo.Cell.2015.rds")
    load("/Volumes/MerckIT/GitHub/Results/Signatures/resistance.program.RData")
    load("/Volumes/MerckIT/GitHub/Results/Signatures/cell.type.sig.RData")
    load("/Volumes/MerckIT/Data/PublicResources/MAPKi.matched/MAPKi.matched.sig.RData")
    load("/Volumes/MerckIT/GitHub/Data/Mel.malignant.genes.RData")
    all.genes<-readRDS("/Volumes/MerckIT/GitHub/Data/Mel.all.data.QC_genes.rds")
    r<-git.compute.samples.res.scores(r = r,res.sig = res.sig,cell.sig = cell.sig,
                                      residu.flag = F,cc.sig = F,num.rounds = 100,even.flag = T)
    r$mapki<-get.sign.scores.bulk(r,mapki.sig,num.rounds = 100)
    r$cc.scores<-get.sign.scores.bulk(r,go.env[c("Melanoma_cell_cycle_Tirosh","G1_S_Tirosh","G2_M_Tirosh")])
    results<-git.prd.MAPKi.Hugo(r,res.sig,cell.sig,mapki.sig,mal.genes,all.genes)
  }
  
  r$tme<-r$tme[,c("B.cell","CAF","Endo.","Macrophage","NK","T.cell")]
  r$y<-r$res[,"res.cmb"]
  f<-function(y){
    r$y<-y
    M1 <- with(r,lmer (y ~ prog + (1 | patient) + tme, mtcars))
    v<-summary(M1)$coefficient["progTRUE",c("Estimate","Pr(>|t|)")]
    return(v)
  }
  
  f.tme<-function(y){
    r$y<-y;M1 <- with(r,lmer (y ~ prog + (1 | patient), mtcars))
    v<-summary(M1)$coefficient["progTRUE",c("Estimate","Pr(>|t|)")]
    return(v)
  }
  
  results<-list(HLM.res = t(apply(r$res,2,f)),
                HLM.mapki = t(apply(r$mapki,2,f)),
                HLM.cc = t(apply(r$cc.scores,2,f)),
                HLM.tme = t(apply(r$tme,2,f.tme)),
                anova = t(anova.mat(cbind(r$res,r$cc.scores,r$mapki),r$patient)),
                HG.res = GO.enrichment.lapply(res.sig,mal.genes,mapki.sig),
                HG.tme = GO.enrichment.lapply(cell.sig,all.genes,mapki.sig))
  save(results,file = "/Volumes/MerckIT/GitHub/Results/Predictors/Prd.MAPKi.res.Hugo2015.RData")
  return(results)
}

git.prd.aPD1.Riaz<-function(rF,riaz.sig,res.sig,cell.sig,mal.genes,post.nivo = NULL){
  if(is.null(rF)){
    rF<-readRDS("/Volumes/MerckIT/GitHub/Data/PublicData/aPD1.Riaz.Cell.2017.rds")
    mal.genes<-readRDS("/Volumes/MerckIT/GitHub/Data/scData/Mel.malignant_genes.rds")
    load("/Volumes/MerckIT/Data/PublicResources/Riaz/Riaz.sig.RData")
    load("/Volumes/MerckIT/GitHub/Results/Signatures/resistance.program.RData")
    load("/Volumes/MerckIT/GitHub/Results/CellTypes/cell.type.sig.RData")
    load("/Volumes/MerckIT/Results/Prognostic/all.ici.sigs.RData")
    if(is.null(rF$res)){
      # rF<-git.compute.samples.res.scores(r = rF,res.sig = iciA.sigs,cell.sig = cell.sig[c("T.CD8","T.cell")])
      rF<-git.compute.samples.res.scores(r = rF,res.sig = res.sig,cell.sig = cell.sig[c("T.CD8","T.cell")])
      # residu.flag = T,cc.sig = go.env[c("G1_S_Tirosh","G2_M_Tirosh")])
    }
    results<-git.prd.aPD1.Riaz(rF,riaz.sig,res.sig,cell.sig,mal.genes,post.nivo = NULL)
    results.naive<-git.prd.aPD1.Riaz(rF,riaz.sig,res.sig,cell.sig,mal.genes,post.nivo = F)
    results.prog<-git.prd.aPD1.Riaz(rF,riaz.sig,res.sig,cell.sig,mal.genes,post.nivo = T)
    save(results,results.naive,results.prog,file = "/Volumes/MerckIT/Results/Prognostic/rsp_riaz.RData")
    save(results,results.naive,results.prog,file = "/Volumes/MerckIT/GitHub/Results/Predictors/rsp_riaz.RData")
    
    intersect(riaz.sig$on.rsp.up,res.sig$uICR.down)
    intersect(riaz.sig$on.rsp.down,res.sig$uICR.up)
    return(results)
  }
  rF<-set.list(rF,rF$response!="NE")
  r<-set.list(rF,rF$on==0);r<-set.list(r,!duplicated(r$patients))
  if(!is.null(post.nivo)){
    r<-set.list(r,r$post.nivo==post.nivo)
    rF<-set.list(rF,rF$post.nivo==post.nivo)
  }
  
  results<-list(data.type = ifelse(is.null(post.nivo),"All",unique(as.character(r$clin.info$Cohort))))
  # Survival
  results$OS<-cox.mat(t(r$res),r)
  # PFS
  r<-set.list(r,r$response!="NE")
  r$survival<-r$pfs
  results$PFS<-cox.mat(t(r$res),r)
  # CR vs. others
  # CR/PR vs. others
  # CR/PR vs. PD
  results$response<-cbind.data.frame(CR.vs.others = t.test.mat(t(r$res),r$response=="CR")[,"zscores"],
                                     CRPR.vs.others = t.test.mat(t(r$res),r$cb)[,"zscores"],
                                     CRPR.vs.PD = t.test.mat(t(r$res[r$response!="SD",]),
                                                             r$cb[r$response!="SD"])[,"zscores"])
  if(!is.null(post.nivo)){return(results)}
  # Compare to signatures:
  # (1) response, (2) molecular response,
  # (3) on-treatment, (4) on-treatment responders
  results$HG<-GO.enrichment.lapply(res.sig,mal.genes,riaz.sig)
  # HLM: (1) test if the resistance program is up/down regulated on-treatment
  #      (2) test if the resistance program is up/down regulated on-treatment specifically in responders
  
  rF<-set.list(rF,rF$response!="NE")
  rF$tmeHLM<-rF$tme[,c("T.CD8")]
  rF$on.cb<-rF$cb&rF$on==1
  rF$on.ncb<-!rF$cb&rF$on==1
  f1<-function(r1,y){
    r1$y<-y
    M1 <- with(r1,lmer (y ~ on + (1 | patients) , mtcars))
    M2 <- with(r1,lmer (y ~ on + (1 | patients) + post.nivo, mtcars))
    M3 <- with(r1,lmer (y ~ on.cr + (1 | patients) , mtcars))
    M4 <- with(r1,lmer (y ~ on.cb + (1 | patients) , mtcars))
    M5 <- with(r1,lmer (y ~ on.ncb + (1 | patients) , mtcars))
    f0<-function(M,x){
      p1<-summary(M)$coefficients
      p1<-get.cor.zscores(p1[,"Estimate"],p1[,"Pr(>|t|)"])
      p1<-p1[x]
      return(p1)
    }
    p<-c(f0(M1,"on"),f0(M2,"post.nivoTRUE"),f0(M3,"on.cr"),f0(M4,"on.cbTRUE"),f0(M5,"on.ncbTRUE"))
    print(p)
    return(p)
  }
  f2<-function(r1,y,x){
    r1$y<-y;r1$x<-x
    M1 <- with(r1,lmer (y ~ on + (1 | patients) + tmeHLM + x, mtcars))
    p1<-summary(M1)$coefficients
    p1<-get.cor.zscores(p1[,"Estimate"],p1[,"Pr(>|t|)"])
    p1<-p1["x"]
    return(p1)
  }
  rF$on.cr<-rF$on.response[,"CR"]
  results$HLM <- t(apply(rF$res,2,function(x) f1(rF,x)))
  HLM.response <- t(apply(rF$res,2,function(y) {
    v<-c(CR = f2(rF,y,rF$on.response[,"CR"]),
         PR = f2(rF,y,rF$on.response[,"PR"]),
         SD = f2(rF,y,rF$on.response[,"SD"]),
         PD = f2(rF,y,rF$on.response[,"PD"]),
         CB = f2(rF,y,ifelse(rF$cb&rF$on==1,1,0)))
    print(v)
    return(v)
  }))
  colnames(HLM.response)<-gsub(".x","",colnames(HLM.response),fixed = T)
  results$HLM <- cbind.data.frame(results$HLM,on = HLM.response)
  return(results)
  
}

git.set.public.ICB.cohorts<-function(compute.scores.flag = T){
  if(compute.scores.flag){
    load("/Volumes/MerckIT/GitHub/Results/Resistance/denovo.sig.names.RData")
    load("/Volumes/MerckIT/GitHub/Results/CellTypes/cell.type.sig.RData")
    load("/Volumes/MerckIT/Results/Prognostic/all.ici.sigs.RData")
    load("/Volumes/MerckIT/GitHub/Data/PublicData/public.ICB.datasets.RData")
    R$aPD1.riaz.pre<-set.list(R$aPD1.riaz,R$aPD1.riaz$on==0)
    R$aPD1.riaz.pre<-set.list(R$aPD1.riaz.pre,!duplicated(R$aPD1.riaz.pre$patients))
    R$aPD1.riaz.pre.prog<-set.list(R$aPD1.riaz.pre,R$aPD1.riaz.pre$post.nivo)
    R<-lapply(R,function(r){
      git.compute.samples.res.scores(r,iciA.sigs,cell.sig = cell.sig,num.rounds = 1000)
    })
    print("Saving ICB cohorts..")
    saveRDS(R,file = "/Volumes/MerckIT/GitHub/Data/PublicData/public.ICB.datasets_withScores.rds")
  }else{
    R<-readRDS("/Volumes/MerckIT/GitHub/Data/PublicData/public.ICB.datasets_withScores.rds")
  }
  return(R)
}

git.set.TCGA<-function(r.tcga = NULL,res.sig,cell.sig){
  if(is.null(r.tcga)){
    r.tcga<-readRDS("/Volumes/MerckIT/GitHub/Data/PublicData/TCGA_SKCM.rds")
    load("/Volumes/MerckIT/GitHub/Results/Signatures/resistance.program.RData")
    load("/Volumes/MerckIT/GitHub/Results/Signatures/cell.type.sig.RData")
  }
  r.tcga<-git.compute.samples.res.scores(r = r.tcga,
                                         res.sig = res.sig,
                                         cell.sig = cell.sig,
                                         cc.sig = NULL,
                                         residu.flag = F,num.rounds = 1000)
  print("Saving the TCGA melanoma cohort..")
  saveRDS(r.tcga,file = "/Volumes/MerckIT/GitHub/Data/PublicData/TCGA_SKCM_withScores.rds")
  return(r.tcga)
}

git.prd.public.ICB.response<-function(R){
  if(is.null(R)){
    R<-git.get.public.ICB.cohorts(compute.scores.flag = T)
  }
  
  test.sig<-function(r){
    r$cb<-is.element(r$response,c("CR","PR"))
    b<-is.element(r$response,c("CR","PR","PD"))
    bA<-!is.element(r$response,"NE")
    m<-cbind.data.frame(CR.vs.others = t.test.mat(t(r$res[b,]),r$response[b]=="CR")[,"zscores"],
                        CR.vs.PRPD = t.test.mat(t(r$res[bA,]),r$response[bA]=="CR")[,"zscores"],
                        CRPR.vs.others = t.test.mat(t(r$res[b,]),r$cb[b])[,"zscores"],
                        CRPR.vs.PD = t.test.mat(t(r$res[bA,]),r$cb[bA])[,"zscores"])
    return(m)
  }
  rearrange<-function(m,R1,comp){
    m1<-t(laply(m,function(x) x[[comp]]))
    rownames(m1)<-colnames(R1[[1]]$res)
    colnames(m1)<-names(R1)
    return(m1)
  }
  idx<-setdiff(names(R),c("mapki","aPD1.riaz"))
  ICB.prd1<-lapply(R[idx],test.sig)
  ICB.prd2<-lapply(colnames(ICB.prd1[[1]]), function(x) rearrange(ICB.prd1,R[idx],x))
  names(ICB.prd2)<-colnames(ICB.prd1[[1]])
  save(ICB.prd1,ICB.prd2,file = "/Volumes/MerckIT/GitHub/Results/Predictors/ExtVal.prf.RData")
  return(ICB.prd2)
}

git.compute.samples.res.scores<-function(r,res.sig,cell.sig,residu.flag = F,
                                         cc.sig = NULL,num.rounds = 1000){
  r$res.ori<-NULL;r$res<-NULL;r$tme<-NULL;r$X<-NULL
  names(res.sig)<-gsub(" ",".",names(res.sig))
  two.sided<-unique(gsub(".up","",gsub(".down","",names(res.sig))))
  b<-is.element(paste0(two.sided,".up"),names(res.sig))&
    is.element(paste0(two.sided,".down"),names(res.sig))
  two.sided<-two.sided[b]
  r$res<-get.sign.scores.bulk(r,res.sig,num.rounds = num.rounds)
  res2<-r$res[,paste0(two.sided,".up")]-r$res[,paste0(two.sided,".down")]
  colnames(res2)<-two.sided
  r$res<-cbind(res2,r$res)
  rownames(r$res)<-colnames(r$tpm)
  r$tme<-get.sign.scores.bulk(r,cell.sig,num.rounds = num.rounds)
  r<-git.cmb.res.scores(r)
  if(residu.flag){
    print("Computing residuales")
    r$res<-t(get.residuals(t(r$res),colSums(r$tpm>0)))
    r$tme<-t(get.residuals(t(r$tme),colSums(r$tpm>0)))
  }
  if(!is.null(cc.sig)){
    print("Filtering cell-cycle effects")
    if(is.null(r$cc.scores)){
      r$cc.scores<-get.sign.scores.bulk(r,cc.sig,num.rounds = num.rounds)
    }
    r$res.ori<-r$res
    r$res<-t(get.residuals(t(r$res),r$cc.scores))
  }
  if(is.element("res.cmb",colnames(r$res))){
    r$res<-cbind.data.frame(r$res,res.final.minus.TIL = r$res[,"res.final"]-r$tme[,"T.CD8"])
    if(!is.null(r$res.ori)){
      r$res.ori<-cbind.data.frame(r$res.ori,res.final.minus.TIL = r$res.ori[,"res.final"]-r$tme[,"T.CD8"]) 
    }
  }
  return(r)
}

git.cmb.res.scores<-function(r){
  if(!all(is.element(c("trt","exc","fnc"),colnames(r$res)))){
    return(r)
  }
  r$res<-cbind.data.frame(res.cmb2 = rowMeans(r$res[,c("trt","exc")]),
                          res.cmb2.up = rowMeans(r$res[,c("trt.up","exc.up")]),
                          res.cmb2.down = rowMeans(r$res[,c("trt.down","exc.down")]),
                          res.cmb = rowMeans(r$res[,c("trt","exc","fnc")]),
                          res.cmb.up = rowMeans(r$res[,c("trt.up","exc.up","fnc.up")]),
                          res.cmb.down = rowMeans(r$res[,c("trt.down","exc.down","fnc.down")]),
                          res.cmb.ue = (rowMeans(r$res[,c("trt","exc")])+r$res[,"fnc"])/2,
                          res.cmb.ue.up = (rowMeans(r$res[,c("trt.up","exc.up")])+r$res[,"fnc.up"])/2,
                          res.cmb.ue.down = (rowMeans(r$res[,c("trt.down","exc.down")])+r$res[,"fnc.down"])/2,
                          r$res)
  return(r)
}

exc.cut<-function(results){
  
  sig<-lapply(c(100,150,200,250,300), function(x) tmp.exc.cut(results,x))
  sig<-unlist(sig,recursive = F)
  summary(sig)
  v<-git.prd.aPD1.simple(r.pd1,sig)
  
}

tmp.exc.cut<-function(results,q = 200){
  print("6. Generating the final signatures")
  f<-function(s){
    names(s)<-gsub("p.","",names(s),fixed = T)
    s$exc<-intersect(s$att.neg,s$exc.pos)
    s$att<-intersect(s$att.pos,s$exc.neg)
    return(s)
  }
  b.rp<-!startsWith(r.sc$genes,"RP")
  s1<-c(get.top.elements(m = results$att.sc[,c("p.pos","p.neg")],q = q,min.ci = 0.01,main = "att"),
        get.top.elements(m = results$exc.sc[,c("p.pos","p.neg")],q = q,min.ci = 0.01,main = "exc"))
  s2<-c(get.top.elements(m = results$att.sc[b.rp,c("p.pos","p.neg")],q = q,min.ci = 0.01,main = "att"),
        get.top.elements(m = results$exc.sc[b.rp,c("p.pos","p.neg")],q = q,min.ci = 0.01,main = "exc"))
  s1<-f(s1);s2<-f(s2)
  results$sigA<-s1
  results$sig.no.rp<-s1
  results$sig<-lapply(names(s1),function(x) sort(unique(c(s1[[x]],s2[[x]]))))
  names(results$sig)<-names(s1)
  print(summary(results$sig[c("exc","att")]))
  sig<-results$sig[c("exc","att")]
  names(sig)<-c(paste0("exc.",q,".up"),
                paste0("exc.",q,".down"))
  return(sig)
}

call.get.t.subset.sig<-function(){
  
  load("../Data/all_gene_sets_0204_withSnS.RData")
  rF<-readRDS("/Volumes/MerckIT/GitHub/Data/scData/Mel.all.data.QC.rds")
  
  r4<-readRDS("/Volumes/MerckIT/GitHub/Data/scData/Mel.T.CD4.QC.rds")
  de4<-git.get.t.cell.subset.sig(rF,r4,go.env,T.type = "T.CD4")
  
  r8<-readRDS("/Volumes/MerckIT/GitHub/Data/scData/Mel.T.CD8.QC.rds")
  de8<-git.get.t.cell.subset.sig(rF,r8,go.env,T.type = "T.CD8")
  
  t.sig4<-readRDS("/Volumes/MerckIT/GitHub/Results/CellTypes/T.CD4.subtype.sig.rds")
  t.sig8<-readRDS("/Volumes/MerckIT/GitHub/Results/CellTypes/T.CD8.subtype.sig.rds")
  load("/Volumes/MerckIT/GitHub/Results/Signatures/cell.type.sig.full.RData")
  t.sig<-c(cell.sig[c("T.CD4","T.CD8","T.cell")],t.sig4,t.sig8)
  print(summary(t.sig))
  saveRDS(t.sig,file = "/Volumes/MerckIT/GitHub/Results/CellTypes/T.subtype.sig.rds")
}

git.get.t.cell.subset.sig<-function(rF,r,go.env,T.type = "T.CD8"){
  if(is.null(r)){
    load("../Data/all_gene_sets_0204_withSnS.RData")
    r<-readRDS(paste0("/Volumes/MerckIT/GitHub/Data/scData/Mel.",T.type,".QC.rds"))
    rF<-readRDS("/Volumes/MerckIT/GitHub/Data/scData/Mel.all.data.QC.rds")
    de<-git.get.t.cell.subset.sig(rF,r,go.env,T.type = ifelse(grepl("CD4",r$data.name),"T.CD4","T.CD8"))
  }
  
  print(paste("1. Assign",T.type,"to cell subtypes according to well-established markers."))
  t.sig<-go.env[c("Cytotoxic_T_cell","Inhibitory.receptors","Tirosh_Naive_Tcell")]
  names(t.sig)<-c("cytotoxic","exhausted","naive")
  r$scores<-get.sign.scores(r,t.sig)
  B<-discretize.mat.q(r$scores,q1 = 0.9)
  L<-get.mat(m.rows = r$cells,m.cols = names(t.sig),data = "")
  for(i in names(t.sig)){L[,i]<-ifelse(B[,i],i,"")}
  r$cell.types<-remove.redundant.dashs(paste(L[,1],L[,2],L[,3],sep = "_"))
  if(T.type=="T.CD4"){
    r<-git.find.Tregs(r)
    B<-cbind.data.frame(B,Treg = r$B.treg$Treg)
    r$cell.types[r$B.treg$Treg]<-paste0(r$cell.types[r$B.treg$Treg],"_Treg")
    r$cell.types<-remove.redundant.dashs(r$cell.types)
  }
  r$cell.types[r$cell.types==""]<-"other"
  r1<-set.list(r,rowSums(B)<2)
  print("Cell subtypes:");print(table(r$cell.types))
  print("Cell subtypes:");print(table(r1$cell.types))
  
  print("2. Get subtype signatures.")
  de<-git.get.cell.type.sig(r1,subtype.flag = T)
  de$sig<-remove.ribo(lapply(de$de.sum,function(x) sort(rownames(x)[x$Min>2])))
  print(summary(de$sig))
  
  print("3. Ensure the genes in the subtype signatures do not represent other cell types.")
  u.cell.types<-sort(setdiff(unique(r1$cell.types),"other"))
  de$other.cell.types<-list()
  de$sig.specific<-list()
  cell.types.no.nk<-setdiff(sort(unique(rF$cell.types)),c(T.type,"NK"))
  for(x in u.cell.types){
    genes<-de$sig[[x]]
    b<-is.element(rF$cells,r$cells[r$cell.types==x])&rF$cell.types==T.type
    print(table(b))
    m<-ranksum.t.test.multi.refs(rF$tpm[genes,],b,ref.groups = rF$cell.types)$up
    m<-cbind.data.frame(max.p.no.nk = rowMax(m[,cell.types.no.nk]),m)
    de$sig.specific[[x]]<-rownames(m)[m$max.p<0.05]
    de$sig.specific.no.nk[[x]]<-rownames(m)[m$max.p.no.nk<0.05]
    de$other.cell.types[[x]]<-m
  }
  print("Without NKs");print(summary(de$sig.specific.no.nk))
  print("With NKs");print(summary(de$sig.specific))
  t.sig<-de$sig.specific
  names(t.sig)<-paste0(T.type,".",names(t.sig))
  t.sig<-t.sig[laply(t.sig,length)>5]
  print(summary(t.sig))
  
  saveRDS(r,file = paste0("/Volumes/MerckIT/GitHub/Data/scData/Mel.",T.type,".QC.rds"))
  saveRDS(de,file = paste0("/Volumes/MerckIT/GitHub/Results/CellTypes/",T.type,".subtype.de.rds"))
  saveRDS(t.sig,file = paste0("/Volumes/MerckIT/GitHub/Results/CellTypes/",T.type,".subtype.sig.rds"))
  saveRDS(t.sig,file = paste0("/Volumes/MerckIT/GitHub/Results/Signatures/",T.type,".subtype.sig.rds"))
  return(de)
}

git.find.Tregs<-function(r){
  if(is.null(r)){r<-readRDS("/Volumes/MerckIT/GitHub/Data/scData/Mel.T.CD4.QC.rds")}
  D1<-get.seurat.dataset(r$genes,round(r$tpm,2),D1.name = r$data.name)
  Treg.genes<-c("FOXP3","IL2RA")# IL2RA = CD25
  D1<-AddImputedScore(D1,genes.fit = Treg.genes)
  if(!identical(r$cells,D1@cell.names)){return()}
  r$imputed<-as.matrix(D1@imputed)
  B.treg<-get.mat(r$cells,Treg.genes)
  for(x in Treg.genes){
    par(mfrow=c(2,2))
    hist(r$tpm[x,],100,main = x)
    hist(r$imputed[x,],100,main = paste(x,"imputed"))
    cor.plot(r$tpm[x,],r$imputed[x,])
    b1<-r$tpm[x,]>quantile(r$tpm[x,],0.75,na.rm = T)
    b2<-r$imputed[x,]>quantile(r$imputed[x,],0.75,na.rm = T)
    print(table(b1,b2))
    B.treg[,x]<-b1&b2
  }
  r$B.treg<-cbind.data.frame(Treg = rowSums(B.treg)==2,B.treg)
  return(r)
}

git.cmap.sig<-function(){
  load("/Volumes/MerckIT/GitHub/Results/Resistance/PostTreatment/post.treatment.subsample.de.RData")
  load("/Volumes/MerckIT/GitHub/Results/Resistance/Exclusion/FINAL_Tcell_Exclusion.RData")
  load("/Volumes/MerckIT/GitHub/Results/Resistance/PostTreatment/post.treatment.de.RData")
  f<-function(s){
    names(s)<-gsub("p.","",names(s),fixed = T)
    s$exc<-intersect(s$att.neg,s$exc.pos)
    s$att<-intersect(s$att.pos,s$exc.neg)
    return(s)
  }
  exc.sig<-c(get.top.elements(m = results$att.sc[,c("p.pos","p.neg")],q = 130,min.ci = 0.01,main = "att"),
             get.top.elements(m = results$exc.sc[,c("p.pos","p.neg")],q = 130,min.ci = 0.01,main = "exc"))
  exc.sig<-f(exc.sig)
  summary(exc.sig)
  
  genes<-rownames(trt.subsmp)
  sig.subsmp<-list(up = genes[trt.subsmp[,"F.up"]>quantile(trt.subsmp[,"F.up"],0.9,na.rm = T)],
                   down = genes[trt.subsmp[,"F.down"]>quantile(trt.subsmp[,"F.down"],0.9,na.rm = T)])
  trt.sig<-git.get.scde.ttest.sig(trt.de,q = 100,scde.c = (-1.96),
                                  ttest.c = 0.01,rank.c = 0.01,
                                  ttest.flag = F,rank.flag = T)
  trt.sig<-list(trt.up = intersect(sig.subsmp$up,trt.sig$up),
                trt.down = intersect(sig.subsmp$down,trt.sig$down))
  summary(trt.sig)
  res.sig<-list(up = unique(sort(c(trt.sig$trt.up,exc.sig$exc))),
                down = unique(sort(c(trt.sig$trt.down,exc.sig$att))))
  res.sig$up.no.rib <- remove.ribo(res.sig$up)
  summary(res.sig)
  summary(remove.ribo(res.sig))
  write.csv(list.2.mat(res.sig),file = "~/Desktop/CellRevisions/cmap.res.sig.csv")
  save(res.sig,file = "/Volumes/MerckIT/GitHub/Results/Resistance/cmap.res.sig.RData")
}

function(){
  trg.patients<-unique(r$patients[rowSums(r$B.target[,c("on","post")])>0])
  b.pre.trg<-rowSums(r$B.target[,c("on","post")])==0
  b.trg.pt<-is.element(r$patients,trg.patients)
  table(b.pre.trg,ifelse(b.trg.pt,"trg.patient","no.trg.patient"))
  table(r$B.ici[b.trg.pt,"on"],b.pre.trg[b.trg.pt]) # A sample from a patient that recieved targeted therapy and is now on ICI
  table(r$B.ici[b.trg.pt,"post"],b.pre.trg[b.trg.pt]) # A sample from a patient that recieved targeted therapy and is now on ICI
  b<-(r$B.ici[,"on"]|r$B.ici[,"post"])&b.trg.pt&b.pre.trg
  table(b)
  r$clin.info[b,]
}

git.trt.immune<-function(r = NULL, cell.type = ""){
  if(is.null(r)){
    fileName<-paste0("/Volumes/MerckIT/GitHub/Data/scData/Mel.",cell.type,".QC.RData")
    r<-readRDS(fileName)
  }
  cell.type<-my.gsub(c("Mel.",".QC"),"",r$data.name)
  print(paste("Identifying a post-treatment signature for",cell.type))
  de<-git.find.DEG(r,main = "post.treatment",mainP = "treated")
  de<-git.trt.HLM(r,de)
  fileName<-paste0("/Volumes/MerckIT/GitHub/Results/Resistance/PostTreatment/",
                   cell.type,"_post.treatment.de.RData")
  saveRDS(de,file = fileName)
  return(de)
}

git.trt.HLM<-function(r,de){
  if(max(r$comp.reads)>2){
    print("Transforming complexity")
    r$comp.reads<-log(r$comp.reads)
    r$comp.reads<-r$comp.reads/max(r$comp.reads)
    r$comp<-r$comp/max(r$comp)
  }
  
  library(lme4);library(lmerTest);attach(mtcars)
  f1<-function(r,x){
    r$x<-x
    M1 <- with(r, lmer (x ~ comp.reads + (1 | samples) + treated, mtcars))
    v<-summary(M1)$coefficients["treatedTRUE",c("Estimate","Pr(>|t|)")]
    return(v)
  }
  f2<-function(r,x){
    v<-tryCatch({f1(r,x)},
                error = function(err){return(c(NA,NA))})
    return(v)
  }
  b<-rowSums(r$norm,na.rm = T)==0;r$norm[b,]<-r$tpm[b,]
  P<-t(apply(r$norm[unlist(de$sig),],1,function(x) f2(r,x)))
  P<-cbind.data.frame(Z = get.cor.zscores(P[,"Estimate"],P[,"Pr(>|t|)"]),P)
  de$hlm<-P
  de$sig.hlm<-get.top.cor(P,min.ci = log10(0.1),idx = "Z",q = Inf)
  de$sig.hlm<-list(up = intersect(de$sig$up,de$sig.hlm$Z.up),
                   down = intersect(de$sig$down,de$sig.hlm$Z.down))
  
  sig<-c(de$sig,de$sig.hlm)
  names(sig)<-c("scde.rank.up","scde.rank.down","hlm.up","hlm.down")
  print(summary(sig))
  de$sig<-sig
  r$res<-compute.state.scores(r,sig)
  lapply(colnames(r$res[,1:2]), function(x){
    plot.immuno.de(r,r$res[,x],filpd = grepl("down",x),
                   f = colMeans,subfigureF = T,
                   main = x)
  })
  load("../Data/all_gene_sets_0204_withSnS.RData")
  de$go<-GO.enrichment.lapply(go.env,r$genes,de$sig)
  return(de)
}




