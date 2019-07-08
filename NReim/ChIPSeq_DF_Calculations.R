## I use DESeq2 to get pvalues
  ## Counts over genes
if(F){
  rm(list=ls())
  setwd("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr")
  files <- Sys.glob("*.gr")
  
  load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenes.RData")
  genes <- subset(genes,genes$verification=="Verified"  & genes$chr!="chrM")
  genes <- as(genes,"GRanges")
  
  res <- read.delim("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Normalization.factors.NataliaMayNov2017_CorrOrlando.txt",header = T)
  res$file <- paste0(res$sample,".gr")
  
  libsize=1e7
  count <- c()
  countsc <- c()
  counta <- c()
  countb <- c()
  for(f in files){
    cat(f,"\n")
    load(f)
    ## raw counts
    gr <-resize(gr,200)
    count <- cbind(count,countOverlaps(genes,gr,ignore.strand=T))
    cat(" scale", length(gr),' ',1e7," ",1e7/length(gr) ,"\n")
    
    ## spikein counts : multiplication
    x <- libsize/length(gr)
    x <- x*res[res$file==f,]$normF
    x <- round(x*length(gr))
    cat(" multiply", length(gr),' ',x," ",x/length(gr) ,"\n")
    rn <- sample(1:length(gr),size = x,replace = T)
    gr <- gr[rn]
    counta <- cbind(counta,countOverlaps(genes,gr,ignore.strand=T))
    
    ## spikein counts : division
    x <- libsize/length(gr)
    x <- x/res[res$file==f,]$normF
    x <- round(x*length(gr))
    cat(" div", length(gr),' ',x," ",x/length(gr) ,"\n")
    rn <- sample(1:length(gr),size = x,replace = T)
    gr <- gr[rn]
    countb <- cbind(countb,countOverlaps(genes,gr,ignore.strand=T))
    
    rn <- sample(1:length(gr),size = libsize,replace = T)
    gr <- gr[rn]
    countsc <- cbind(countsc,countOverlaps(genes,gr,ignore.strand=T))
    rm(gr,f,x,rn)
  }
  colnames(countsc) <- colnames(counta) <- colnames(countb) <- colnames(count) <- gsub(".gr","",files)
  cf <- cbind(as.data.frame(genes),as.data.frame(count))
  #save(cf,file="A:/work/WinstonLab/Natalia_Reim/DataFrames/RawCounts_Signal_genesSGD.RData")
  #write.table(cf,file="RawCountsperGene_ChIPAssay.txt",quote = F,sep = "\t",row.names = F)
  cf <- cbind(as.data.frame(genes),countsc)
  #save(cf,file="A:/work/WinstonLab/Natalia_Reim/DataFrames/LibSizeScaledCounts_Signal_genesSGD.RData")
  cf <- cbind(as.data.frame(genes),counta)
  save(cf,file="A:/work/WinstonLab/Natalia_Reim/DataFrames/SpikeinScaled_MultipCounts_Signal_genesSGD.RData")
  cf <- cbind(as.data.frame(genes),countb)
  save(cf,file="A:/work/WinstonLab/Natalia_Reim/DataFrames/SpikeinScaled_DivCounts_Signal_genesSGD.RData")
  
}
## Counts over first 500bp for H3
if(F){
  rm(list=ls())
  setwd("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr")
  files <- Sys.glob("*.gr")
  #files <- files[grep("H3_ChIP",files)]
  
  load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenes.RData")
  genes <- subset(genes,genes$verification=="Verified"  & genes$chr!="chrM")
  genes <- as(genes,"GRanges")
  genes <- flank(resize(genes,1,"start"),500,start = F)
  
  res <- read.delim("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Normalization.factors.NataliaMayNov2017_CorrOrlando.txt",header = T)
  res$file <- paste0(res$sample,".gr")
  
  libsize=1e7
  count <- c()
  countsc <- c()
  counta <- c()
  countb <- c()
  for(f in files){
    cat(f,"\n")
    load(f)
    ## raw counts
    gr <-resize(gr,200)
    count <- cbind(count,countOverlaps(genes,gr,ignore.strand=T))
    cat(" scale", length(gr),' ',1e7," ",1e7/length(gr) ,"\n")
    
    ## spikein counts : multiplication
    x <- libsize/length(gr)
    x <- x*res[res$file==f,]$normF
    x <- round(x*length(gr))
    cat(" multiply", length(gr),' ',x," ",x/length(gr) ,"\n")
    rn <- sample(1:length(gr),size = x,replace = T)
    gr <- gr[rn]
    counta <- cbind(counta,countOverlaps(genes,gr,ignore.strand=T))
    
    ## spikein counts : division
    x <- libsize/length(gr)
    x <- x/res[res$file==f,]$normF
    x <- round(x*length(gr))
    cat(" div", length(gr),' ',x," ",x/length(gr) ,"\n")
    rn <- sample(1:length(gr),size = x,replace = T)
    gr <- gr[rn]
    countb <- cbind(countb,countOverlaps(genes,gr,ignore.strand=T))
    
    rn <- sample(1:length(gr),size = libsize,replace = T)
    gr <- gr[rn]
    countsc <- cbind(countsc,countOverlaps(genes,gr,ignore.strand=T))
    rm(gr,f,x,rn)
  }
  colnames(countsc) <- colnames(counta) <- colnames(countb) <- colnames(count) <- gsub(".gr","",files)
  cf <- cbind(as.data.frame(genes),as.data.frame(count))
  cf <- cbind(as.data.frame(genes),countsc)
  cf <- cbind(as.data.frame(genes),counta)
  save(cf,file="A:/work/WinstonLab/Natalia_Reim/DataFrames/SpikeinScaled_MultipCounts_Signal_H3binSGD.RData")
  cf <- cbind(as.data.frame(genes),countb)
  #save(cf,file="A:/work/WinstonLab/Natalia_Reim/DataFrames/SpikeinScaled_DivCounts_Signal_genesSGD.RData")
  
}
## Counts over intron/exons 
if(T){
  rm(list = ls())
  genome.anno.path="A:/work/yeast_Annotations/Yeast_ref_Annotations_DJ.gff"
  d <- read.delim(genome.anno.path,header=F,comment.char = "#")
  d$tracking_id = gsub("gene_id ","",gsub("; transcript_id.*","",d$V9))
  d$exon.number <- as.numeric(gsub(";","",gsub("exon_number ", "", stringr::str_extract(d$V9,"exon_number \\S*?;"))))
  allIntrons <- read.xlsx("A:/work/WinstonLab/Natalia_Reim/Paper_Thesis_Figs/IntronsAndRPGenes.xlsx",sheetIndex = 1)
  d <- d[d$tracking_id%in%allIntrons$reference & d$V3=="exon",]
  two.introns <- names(sort(table(allIntrons$reference),decreasing = T)[1:8])
  allIntrons$intron.numer <- 1
  allIntrons[c(18,27,58,73,86,96,168,242),]$intron.numer <- 2
  allIntrons$feature <- "intron"
  allIntrons <- allIntrons[,c(2:5,1,6,14,13)]
  d <- d[,c(1,4,5,7,10,10,3,11)]
  names(d) <- names(allIntrons) <- c("chr","start","end","strand","tracking_id","gene","feature","number")
  genes <- rbind(d,allIntrons)
  genes <- as(genes,"GRanges")
  rm(allIntrons,d,two.introns)
  
  
  setwd("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr")
  files <- Sys.glob("*.gr")
  res <- read.delim("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Normalization.factors.NataliaMayNov2017_CorrOrlando.txt",header = T)
  res$file <- paste0(res$sample,".gr")
  res$file%in%files
 
  libsize=1e7
  count <- c()
  counta <- c()
  for(f in files){
    cat(f,"\n")
    load(f)
    ## raw counts
    gr <-resize(gr,200)
    count <- cbind(count,countOverlaps(genes,gr,ignore.strand=T))
    cat(" scale", length(gr),' ',1e7," ",1e7/length(gr) ,"\n")
    
    ## spikein counts : multiplication
    x <- libsize/length(gr)
    x <- x*res[res$file==f,]$normF
    x <- round(x*length(gr))
    cat(" multiply", length(gr),' ',x," ",x/length(gr) ,"\n")
    rn <- sample(1:length(gr),size = x,replace = T)
    gr <- gr[rn]
    counta <- cbind(counta,countOverlaps(genes,gr,ignore.strand=T))
    
    rm(gr,f,x,rn)
  }
  colnames(counta) <- colnames(count) <- gsub(".gr","",files)
  cf <- cbind(as.data.frame(genes),as.data.frame(count))
  save(cf,file="A:/work/WinstonLab/Natalia_Reim/DataFrames/RawCounts_Signal_IntronExonSGD.RData")
  #write.table(cf,file="RawCountsperGene_ChIPAssay.txt",quote = F,sep = "\t",row.names = F)
  cf <- cbind(as.data.frame(genes),counta)
  save(cf,file="A:/work/WinstonLab/Natalia_Reim/DataFrames/SpikeinScaled_MultipCounts_Signal_IntronExonSGD.RData")
  
  
}
  

## Spikein normalized Rel occupancy per gene
if(F){
  rm(list=ls())
  ## get my_relative_signal function from MultigeneHeatMeta.R 
  res <- read.delim("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Normalization.factors.NataliaMayNov2017_CorrOrlando.txt",header = T)
  
  protein <- list()
  protein[["H3K4me3"]] <- "H3"
  protein[["H3K36me2"]] <- "H3"
  protein[["H3K36me3"]] <- "H3 (round1)"
  m.rel.occup <- c()
  for(i in  1:length(protein)){
    f = paste("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_",names(protein)[i],".RData",sep="")
    fc = paste("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_",protein[[i]],".RData",sep="")
    
    if(file.exists(f) & file.exists(fc)){
      cat(names(protein)[i],"\n")
      cat(res[res$Factor==names(protein)[i],]$normF,"\n")
      cat(protein[[i]],"\n")
      cat(res[res$Factor==protein[[i]],]$normF,"\n")
      
      load(f)
      l1 <- l; rm(l);
      load(fc)
      l2 <- l; rm(l);
      normF1=res[res$Factor==names(protein)[[i]],]$normF
      normF2=res[res$Factor==protein[[i]],]$normF
      
      l <-my_relative_signal(l1,l2,normF1,normF2)
      rm(l1,l2)
      
      if(i==1){
        m.rel.occup <- data.frame(tracking_id=rownames(l[[1]]))
      }
      myfun <- function(m){ apply(m,1,sum) }
      cf <- data.frame(tracking_id=rownames(l[[1]]))
      for(j in 1:length(l)){
        cf  <- cbind(cf,myfun(l[[j]])) 
      }
      names(cf)[2:ncol(cf)] <- names(l)
      cf$x <- apply(cf[,2:3],1,function(x) sqrt(prod(x)))
      cf$y <- apply(cf[,4:5],1,function(x) sqrt(prod(x)))
      names(cf)[6:7] <- c(paste0(names(protein)[i],"_Depl."),paste0(names(protein)[i]))
      m.rel.occup <- cbind(m.rel.occup,cf[,6:7])
    }
  }
  save(m.rel.occup,file="A:/work/WinstonLab/Natalia_Reim/DataFrames/relOccupancyChIPSeq.RData")    
}


## Log2 factor occupancy changes over genes
if(T){
  
  ################## NOTE
  ## 1) Orlando et al suggest computing spikein normalization factors by incorporating % of mixed spikein chromatin
  ## 2) Usually for RNASeq data is it considered constant
  ## 3) For ChIPSeq data it cn be obtained from input sequencing
  ## 4) If DeSeq2 is used to calculate spikein size-factors, this contribution is lost
  ## 5) Howver, log2FC values for spikein genome center around 0 (as we want)
  ## 6) If I use Orlando et al spikein factors for DESeq2 computation, the log2FC values for spikein are thrown off in positive direction
  ## 7) It is unclear which is correct approach - as for metagene plots, I use earlier Orlando approach (which is more logical)
  ## 8) For computing log2FC values, i decide to use Deseq2 spikein size factor approach as this makes sure taht spikein log2FC values are centered around 0
  
  ## Approach:
  ## 1) Resample the ChIPseq libraries by adjusting for the spikein proportions per million
  ## 2) Perform DESeq2 analysis using unit size factors (i.e. 1,1,1,1) to preserve spikein proportions 
  
  rm(list=ls()) 
  mydeseq2fun <- function(m,type,sizeF){
    require(DESeq2)
    condition <- gsub("_[1234]","",type)
    colData <- data.frame(condition,type)
    row.names(colData) <- colData$type
    dds <- DESeqDataSetFromMatrix(countData = as.matrix(m),colData = colData,design = ~ condition)
    if(length(sizeF)==length(type)){
      cat("Using the provided size Factors\n")
      sizeFactors(dds) <- sizeF
    }else{
      cat("Computing the size factors from data\n")
    }
    dds <- DESeq(dds)
    res <- as.data.frame(results(dds))
    z <- round(as.data.frame(DESeq2::counts(dds,normalized=T)),2)
    res$tracking_id <- rownames(m)
    res <- cbind(res,z)
    res <- merge(res,genes,by="tracking_id",all.x=T)
    res <- res[,-c(4,5,6)]
    names(res)[4] <- "FDR"
    j = names(dds@rowRanges@elementMetadata@listData)[12]
    j <- gsub("condition_","",j)
    j = unlist(strsplit(j,"_vs_"))
    res$log2FoldChange <- res$log2FoldChange*(-1)
    names(res)[3] <- paste0("log2 (",j[2],"/",j[1],")")
    return(res)
  }
  
  load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenes.RData")
  genes <- subset(genes,genes$verification=="Verified" & genes$species=="Scer" & genes$chr!="chrM")
  load("A:/work/yeast_Annotations/refs.RData")
  genes <- merge(genes,unique(refs[,c("tracking_id","gene","link")]),by="tracking_id",all.x=T)
  res <- read.delim("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Normalization.factors.NataliaMayNov2017_CorrOrlando.txt",header=T)
  #res <- read.delim("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Normalization.factors.NataliaMayNov2017_revised.txt",header=T)
  normF <- res$normF
  names(normF) <- res$sample
  normF <- 1/normF
  
  
  ### (1) Spiken normalization
  load("A:/work/WinstonLab/Natalia_Reim/DataFrames/SpikeinScaled_MultipCounts_Signal_H3binSGD.RData")
  #load("A:/work/WinstonLab/Natalia_Reim/DataFrames/SpikeinScaled_MultipCounts_Signal_genesSGD.RData")
  #load("A:/work/WinstonLab/Natalia_Reim/DataFrames/LibSizeScaledCounts_Signal_genesSGD.RData")
  rownames(cf) <- cf$tracking_id
  m0 <- cf[cf$species=="Scer",16:79]
  mr <- cf[cf$species=="Spom",16:79]
  m0 <- m0[,names(normF)]
  mr <- mr[,names(normF)]
  l <- list()
  ln <- list()
  for(i in seq(1,45,4)){
    cat(names(m0)[i:(i+3)],"\n")
    cat(names(normF[i:(i+3)]),"\n")
    cat(normF[i:(i+3)],"\n")
    cat(paste(as.character(res$Factor[i])),"\n")
    l[[paste(as.character(res$Factor[i]))]] <- mydeseq2fun(m = round(mr[,i:(i+3)]),type=names(mr)[i:(i+3)],sizeF= c(1,1,1,1))  #normF[i:(i+3)]
    ln[[paste(as.character(res$Factor[i]))]] <- mydeseq2fun(m = round(m0[,i:(i+3)]),type=names(m0)[i:(i+3)],sizeF=c(1,1,1,1))  #normF[i:(i+3)]
    ## If libraries are only scaled for their sizes, following lines of codes can be used
    ## The results are highly correlative to the final results!
    #l[[paste(as.character(res$Factor[i]))]] <- mydeseq2fun(m = round(mr[,i:(i+3)]),type=names(mr)[i:(i+3)],sizeF= normF[i:(i+3)])  #normF[i:(i+3)]
    #ln[[paste(as.character(res$Factor[i]))]] <- mydeseq2fun(m = round(m0[,i:(i+3)]),type=names(m0)[i:(i+3)],sizeF=normF[i:(i+3)])  #normF[i:(i+3)]
  }
  rm(mr,m0,cf)
  
  ## without spikein normalization
  load("A:/work/WinstonLab/Natalia_Reim/DataFrames/RawCounts_ChIPSignal_genesSGD.RData")
  rownames(cf) <- cf$tracking_id
  m0 <- cf[cf$species=="Scer",16:79]
  mr <- cf[cf$species=="Spom",16:79]
  m0 <- m0[,names(normF)]
  mr <- mr[,names(normF)]
  w <- list()
  wn <- list()
  for(i in seq(1,45,4)){
    cat(names(m0)[i:(i+3)],"\n")
    cat(names(normF[i:(i+3)]),"\n")
    cat(normF[i:(i+3)],"\n")
    cat(paste(as.character(res$Factor[i])),"\n")
    w[[paste(as.character(res$Factor[i]))]] <- mydeseq2fun(m = round(mr[,i:(i+3)]),type=names(mr)[i:(i+3)],sizeF= 0)  #normF[i:(i+3)]
    wn[[paste(as.character(res$Factor[i]))]] <- mydeseq2fun(m = round(m0[,i:(i+3)]),type=names(m0)[i:(i+3)],sizeF=0)  #normF[i:(i+3)]
  }
  rm(mr,m0,cf)
  save(l,ln,w,wn,file="A:/work/WinstonLab/Natalia_Reim/DataFrames/Diff_Factor_Occupancy_DESeq2SGD.RData")
  
  
  #load(file="A:/work/WinstonLab/Natalia_Reim/DataFrames/Diff_Factor_Occupancy_DESeq2SGD.RData")
  pl <- list()
  pln <- list()
  for(i in 1:length(l)){
    pl[[paste(names(l)[i])]] <- l[[i]][,3]
    pln[[paste(names(ln)[i])]] <- ln[[i]][,3]
  }
  boxplot(pl,ylim=c(-3,3))
  boxplot(pln,ylim=c(-3,2))
  abline(h=0)
  if(F){
    for(i in 1:12){
      if(names(ln)[i]%in%c("S2P","S5P","RNAPII","RNAPII (round1)")){
        cat(names(ln)[i],": choosing non-spikein l2fc values\n")
        write.table(wn[[i]],file= paste("A:/work/WinstonLab/Natalia_Reim/DataFrames/",names(wn)[i],"ChIP_OccupancyLog2FC_DESeq2CorrOrlando.txt",sep=""),quote = F,sep = "\t",row.names = F)
      }else{
        write.table(ln[[i]],file= paste("A:/work/WinstonLab/Natalia_Reim/DataFrames/",names(ln)[i],"ChIP_OccupancyLog2FC_DESeq2CorrOrlando.txt",sep=""),quote = F,sep = "\t",row.names = F)
      }
    }
  }
  
  ### H3 signal in 500bp bin downstream of TSS
  if(F){
    for(i in 1:2){
      write.table(ln[[i]],file= paste("A:/work/WinstonLab/Natalia_Reim/DataFrames/",names(ln)[i],"ChIP_OccupancyLog2FC_DESeq2CorrOrlando_H3Bins.txt",sep=""),quote = F,sep = "\t",row.names = F)
    }
  }
  
  m.spike <- ln[[1]]$`log2 (Spn1.Depl_H3_ChIP/Spn1_H3_ChIP)`
  fdr.spike <- ln[[1]]$FDR
  m <- wn[[1]]$`log2 (Spn1.Depl_H3_ChIP/Spn1_H3_ChIP)`
  fdr <- wn[[1]]$FDR
  for(i in 2:length(ln)){
    m.spike <- cbind(m.spike,ln[[i]][,3])
    fdr.spike <- cbind(fdr.spike,ln[[i]][,4])
    m <- cbind(m,wn[[i]][,3])
    fdr <- cbind(fdr,wn[[i]][,4])
  }
  m.spike <- as.data.frame(m.spike)
  rownames(m.spike) <- ln[[1]]$tracking_id
  names(m.spike) <- names(ln)
  fdr.spike <- as.data.frame(fdr.spike)
  rownames(fdr.spike) <- ln[[1]]$tracking_id
  names(fdr.spike) <- names(ln)
  m <- as.data.frame(m)
  rownames(m) <- wn[[1]]$tracking_id
  names(m) <- names(wn)
  fdr <- as.data.frame(fdr)
  rownames(fdr) <- wn[[1]]$tracking_id
  names(fdr) <- names(wn)

  m.occup <- data.frame(tracking_id=wn[[1]]$tracking_id)
  m.spike.occup <- data.frame(tracking_id=ln[[1]]$tracking_id)
  for(i in 1:length(ln)){
    cat(names(ln)[i],"\n")
    cat(" ",paste0(names(ln[[i]][,5:8]),collapse = ","),"\n")
    qwe <- ln[[i]][,5:8]
    qwe$x <- apply(qwe[,1:2],1,function(x) sqrt(prod(x)))
    qwe$y <- apply(qwe[,3:4],1,function(x) sqrt(prod(x)))
    names(qwe)[5:6] <- names(qwe)[c(1,3)]
    qwe <- qwe[,5:6]
    names(qwe) <- gsub("_[1234567]$","",names(qwe))
    m.spike.occup <- cbind(m.spike.occup,qwe)
    rm(qwe)
    
    qwe <- wn[[i]][,5:8]
    cat(" ",paste0(names(wn[[i]][,5:8]),collapse = ","),"\n")
    qwe$x <- apply(qwe[,1:2],1,function(x) sqrt(prod(x)))
    qwe$y <- apply(qwe[,3:4],1,function(x) sqrt(prod(x)))
    names(qwe)[5:6] <- names(qwe)[c(1,3)]
    qwe <- qwe[,5:6]
    names(qwe) <- gsub("_[1234567]$","",names(qwe))
    m.occup <- cbind(m.occup,qwe)
    rm(qwe)
    cat("\n\n")
  }
  names(m.occup)[4:5] <- gsub("H3_ChIP","H3 (round1)_ChIP",names(m.occup)[4:5])
  names(m.spike.occup)[4:5] <- gsub("H3_ChIP","H3 (round1)_ChIP",names(m.spike.occup)[4:5])
  names(m.occup)[14:15] <- gsub("RNAPII","RNAPII (round1)",names(m.occup)[14:15])
  names(m.spike.occup)[14:15] <- gsub("RNAPII","RNAPII (round1)",names(m.spike.occup)[14:15])
  
  save(m,m.spike,m.occup,m.spike.occup,fdr.spike,fdr,file="A:/work/WinstonLab/Natalia_Reim/DataFrames/Log2FC values for the ChIPFactors_DESeq2.RData")

}

## relative log2FC values over genes
if(T){
  rm(list=ls())
  load("A:/work/WinstonLab/Natalia_Reim/DataFrames/Diff_Factor_Occupancy_DESeq2SGD.RData")
  myl <- list()
  myl[[ "H3K36me2"]] <- "H3"
  myl[[ "H3K4me3"]] <- "H3"
  myl[[ "H3K36me3"]] <- "H3 (round1)"
  myl[[ "S2P"]] <- "RNAPII"
  myl[[ "S5P"]] <- "RNAPII"
  
  for(i in 1:length(myl)){
    prot= names(myl)[i]
    control = myl[[i]]
    if(names(myl)[i]%in%c("S2P","S5P")){
      cat("choosing non-spikein l2fc values\n")
      m <- wn[[prot]]
      n <- wn[[control]]
    }else{
      m <- ln[[prot]]
      n <- ln[[control]]
    }
    z <- m
    z[,3] <- m[,3]-n[,3]
    z[,4] <- "*"
    z[,5:8] <- (m[,5:8]*1+1)/(n[,5:8]*1+1)
    write.table(z,file=paste0("A:/work/WinstonLab/Natalia_Reim/DataFrames/RelLog2FC_",prot,"_",control,"_CorrOrlando.txt"),quote = F,sep = "\t",row.names = F)
  }
}

## Log2 factor occupancy changes/relative log2FC values over Introns and exons
if(T){
  
  ################## NOTE
  ## 1) Orlando et al suggest computing spikein normalization factors by incorporating % of mixed spikein chromatin
  ## 2) Usually for RNASeq data is it considered constant
  ## 3) For ChIPSeq data it cn be obtained from input sequencing
  ## 4) If DeSeq2 is used to calculate spikein size-factors, this contribution is lost
  ## 5) Howver, log2FC values for spikein genome center around 0 (as we want)
  ## 6) If I use Orlando et al spikein factors for DESeq2 computation, the log2FC values for spikein are thrown off in positive direction
  ## 7) It is unclear which is correct approach - as for metagene plots, I use earlier Orlando approach (which is more logical)
  ## 8) For computing log2FC values, i decide to use Deseq2 spikein size factor approach as this makes sure taht spikein log2FC values are centered around 0
  
  ## Approach:
  ## 1) Resample the ChIPseq libraries by adjusting for the spikein proportions per million
  ## 2) Perform DESeq2 analysis using unit size factors (i.e. 1,1,1,1) to preserve spikein proportions 
  
  rm(list=ls()) 
  library(xlsx)
  mydeseq2fun <- function(m,genes,type,sizeF){
    require(DESeq2)
    condition <- gsub("_[1234]","",type)
    colData <- data.frame(condition,type)
    row.names(colData) <- colData$type
    dds <- DESeqDataSetFromMatrix(countData = as.matrix(m),colData = colData,design = ~ condition)
    if(length(sizeF)==length(type)){
      cat("Using the provided size Factors\n")
      sizeFactors(dds) <- sizeF
    }else{
      cat("Computing the size factors from data\n")
    }
    dds <- DESeq(dds)
    res <- as.data.frame(results(dds))
    z <- round(as.data.frame(DESeq2::counts(dds,normalized=T)),2)
    res$tracking_id <- rownames(m)
    res <- cbind(res,z)
    res <- merge(res,genes,by="tracking_id",all.x=T)
    res <- res[,-c(4,5,6)]
    names(res)[4] <- "FDR"
    j = names(dds@rowRanges@elementMetadata@listData)[12]
    j <- gsub("condition_","",j)
    j = unlist(strsplit(j,"_vs_"))
    res$log2FoldChange <- res$log2FoldChange*(-1)
    names(res)[3] <- paste0("log2 (",j[2],"/",j[1],")")
    return(res)
  }
  
  
  genome.anno.path="A:/work/yeast_Annotations/Yeast_ref_Annotations_DJ.gff"
  d <- read.delim(genome.anno.path,header=F,comment.char = "#")
  d$tracking_id = gsub("gene_id ","",gsub("; transcript_id.*","",d$V9))
  d$exon.number <- as.numeric(gsub(";","",gsub("exon_number ", "", stringr::str_extract(d$V9,"exon_number \\S*?;"))))
  allIntrons <- read.xlsx("A:/work/WinstonLab/Natalia_Reim/Paper_Thesis_Figs/IntronsAndRPGenes.xlsx",sheetIndex = 1)
  d <- d[d$tracking_id%in%allIntrons$reference & d$V3=="exon",]
  two.introns <- names(sort(table(allIntrons$reference),decreasing = T)[1:8])
  allIntrons$intron.numer <- 1
  allIntrons[c(18,27,58,73,86,96,168,242),]$intron.numer <- 2
  allIntrons$feature <- "intron"
  allIntrons <- allIntrons[,c(2:5,1,6,14,13)]
  d <- d[,c(1,4,5,7,10,10,3,11)]
  names(d) <- names(allIntrons) <- c("chr","start","end","strand","tracking_id","gene","feature","number")
  genes <- rbind(d,allIntrons)
  genes <- as(genes,"GRanges")
  rm(allIntrons,d,two.introns)
  
  
  res <- read.delim("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Normalization.factors.NataliaMayNov2017_CorrOrlando.txt",header=T)
  normF <- res$normF
  names(normF) <- res$sample
  normF <- 1/normF
  
  g <- as.data.frame(genes)
  g$tracking_id <- paste0(g$tracking_id,"_",g$feature,g$number)
  g$species <- "Scer"
  
  ### (1) Spiken normalization
  load("A:/work/WinstonLab/Natalia_Reim/DataFrames/SpikeinScaled_MultipCounts_Signal_IntronExonSGD.RData")
  cf$tracking_id <- paste0(cf$tracking_id,"_",cf$feature,cf$number)
  cf$species <- "Scer"
  rownames(cf) <- cf$tracking_id
  m0 <- cf[cf$species=="Scer",10:74]
  m0 <- m0[,names(normF)]
  ln <- list()
  for(i in seq(1,45,4)){
    cat(names(m0)[i:(i+3)],"\n")
    cat(names(normF[i:(i+3)]),"\n")
    cat(normF[i:(i+3)],"\n")
    cat(paste(as.character(res$Factor[i])),"\n")
    ln[[paste(as.character(res$Factor[i]))]] <- mydeseq2fun(m = round(m0[,i:(i+3)]),type=names(m0)[i:(i+3)],genes = g,sizeF=c(1,1,1,1))  #normF[i:(i+3)]
  }
  rm(m0,cf)
  
  ## without spikein normalization
  load("A:/work/WinstonLab/Natalia_Reim/DataFrames/RawCounts_Signal_IntronExonSGD.RData")
  cf$tracking_id <- paste0(cf$tracking_id,"_",cf$feature,cf$number)
  cf$species <- "Scer"
  rownames(cf) <- cf$tracking_id
  m0 <- cf[cf$species=="Scer",10:74]
  m0 <- m0[,names(normF)]
  wn <- list()
  for(i in seq(1,45,4)){
    cat(names(m0)[i:(i+3)],"\n")
    cat(names(normF[i:(i+3)]),"\n")
    cat(normF[i:(i+3)],"\n")
    cat(paste(as.character(res$Factor[i])),"\n")
    wn[[paste(as.character(res$Factor[i]))]] <- mydeseq2fun(m = round(m0[,i:(i+3)]),genes = g,type=names(m0)[i:(i+3)],sizeF=0)  #normF[i:(i+3)]
  }
  rm(m0,cf)
  save(ln,wn,file="A:/work/WinstonLab/Natalia_Reim/DataFrames/Diff_Factor_Occupancy_DESeq2SGD_IntronExons.RData")
  
  
  load(file="A:/work/WinstonLab/Natalia_Reim/DataFrames/Diff_Factor_Occupancy_DESeq2SGD_IntronExons.RData")
  pln <- list()
  for(i in 1:length(ln)){
    pln[[paste(names(ln)[i])]] <- ln[[i]][,3]
  }
  boxplot(pln,ylim=c(-3,2))
  abline(h=0)
  
  m.spike <- ln[[1]]$`log2 (Spn1.Depl_H3_ChIP/Spn1_H3_ChIP)`
  fdr.spike <- ln[[1]]$FDR
  m <- wn[[1]]$`log2 (Spn1.Depl_H3_ChIP/Spn1_H3_ChIP)`
  fdr <- wn[[1]]$FDR
  for(i in 2:length(ln)){
    m.spike <- cbind(m.spike,ln[[i]][,3])
    fdr.spike <- cbind(fdr.spike,ln[[i]][,4])
    m <- cbind(m,wn[[i]][,3])
    fdr <- cbind(fdr,wn[[i]][,4])
  }
  m.spike <- as.data.frame(m.spike)
  rownames(m.spike) <- ln[[1]]$tracking_id
  names(m.spike) <- names(ln)
  fdr.spike <- as.data.frame(fdr.spike)
  rownames(fdr.spike) <- ln[[1]]$tracking_id
  names(fdr.spike) <- names(ln)
  m <- as.data.frame(m)
  rownames(m) <- wn[[1]]$tracking_id
  names(m) <- names(wn)
  fdr <- as.data.frame(fdr)
  rownames(fdr) <- wn[[1]]$tracking_id
  names(fdr) <- names(wn)
  
  m.occup <- data.frame(tracking_id=wn[[1]]$tracking_id)
  m.spike.occup <- data.frame(tracking_id=ln[[1]]$tracking_id)
  for(i in 1:length(ln)){
    cat(names(ln)[i],"\n")
    cat(" ",paste0(names(ln[[i]][,5:8]),collapse = ","),"\n")
    qwe <- ln[[i]][,5:8]
    qwe$x <- apply(qwe[,1:2],1,function(x) sqrt(prod(x)))
    qwe$y <- apply(qwe[,3:4],1,function(x) sqrt(prod(x)))
    names(qwe)[5:6] <- names(qwe)[c(1,3)]
    qwe <- qwe[,5:6]
    names(qwe) <- gsub("_[1234567]$","",names(qwe))
    m.spike.occup <- cbind(m.spike.occup,qwe)
    rm(qwe)
    
    qwe <- wn[[i]][,5:8]
    cat(" ",paste0(names(wn[[i]][,5:8]),collapse = ","),"\n")
    qwe$x <- apply(qwe[,1:2],1,function(x) sqrt(prod(x)))
    qwe$y <- apply(qwe[,3:4],1,function(x) sqrt(prod(x)))
    names(qwe)[5:6] <- names(qwe)[c(1,3)]
    qwe <- qwe[,5:6]
    names(qwe) <- gsub("_[1234567]$","",names(qwe))
    m.occup <- cbind(m.occup,qwe)
    rm(qwe)
    cat("\n\n")
  }
  names(m.occup)[4:5] <- gsub("H3_ChIP","H3 (round1)_ChIP",names(m.occup)[4:5])
  names(m.spike.occup)[4:5] <- gsub("H3_ChIP","H3 (round1)_ChIP",names(m.spike.occup)[4:5])
  names(m.occup)[14:15] <- gsub("RNAPII","RNAPII (round1)",names(m.occup)[14:15])
  names(m.spike.occup)[14:15] <- gsub("RNAPII","RNAPII (round1)",names(m.spike.occup)[14:15])
  
  
  
  myl <- list()
  myl[[ "H3K36me2"]] <- "H3"
  myl[[ "H3K4me3"]] <- "H3"
  myl[[ "H3K36me3"]] <- "H3 (round1)"
  myl[[ "S2P"]] <- "RNAPII"
  myl[[ "S5P"]] <- "RNAPII"
  m.rel <- data.frame(tracking_id=wn[[1]]$tracking_id)
  for(i in 1:length(myl)){
    prot= names(myl)[i]
    control = myl[[i]]
    if(names(myl)[i]%in%c("S2P","S5P")){
      cat("choosing non-spikein l2fc values\n")
      mt <- wn[[prot]]
      nt <- wn[[control]]
    }else{
      mt <- ln[[prot]]
      nt <- ln[[control]]
    }
    z <- mt
    z[,3] <- mt[,3]-nt[,3]
    z[,4] <- "*"
    z[,5:8] <- (mt[,5:8]*1+1)/(nt[,5:8]*1+1)
    m.rel <- cbind(m.rel,z[,3])
  }
  names(m.rel)[2:6] <- names(myl)
  save(m.rel,genes,m,m.spike,m.occup,m.spike.occup,fdr.spike,fdr,file="A:/work/WinstonLab/Natalia_Reim/DataFrames/Log2FC values for the ChIPFactors_DESeq2_IntronsExons.RData")
}


### Generate bedgraphs for ChIPSeq with modified Orlando (merged replicates)
if(T){
  standardize.cov <- function(covs,seq.depth=1e7){
    total.cov <- sum(as.numeric(unlist(covs)))
    covss <- round(((covs * seq.depth) / total.cov), 4)
    return(covss)
  }
  spiked.chr <- c("chrXVII","chrXVIII","chrXIX")
  chr.flt <- c("chrI","chrII","chrIII","chrIV","chrIX","chrV","chrVI","chrVII","chrVIII","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
  res <- read.delim("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Normalization.factors.NataliaMayNov2017_CorrOrlando.txt",header = T)
  res$file <- paste0(res$sample,".gr")
  res <- res[res$Factor%in%c("RNAPII","S2P","S5P","Spt6","H3"),]
  res$spl <- paste0(res$condition,"_",res$Factor)
  res <- split(res,res$spl)
  
  for(i in 1:length(res)){
    for (j in 1:length(res[[i]]$file)){
      f= paste0("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/",as.character(res[[i]]$file[j]))
      normF = as.numeric(res[[i]]$normF[j])
      cat(f," ", normF,"\n")
      load(f)
      gr <- subset(gr,seqnames(gr)%in%c(chr.flt,spiked.chr))
      gr <- keepSeqlevels(gr, c(chr.flt,spiked.chr))
      cov1 <- standardize.cov(coverage(gr))
      cov1 <- cov1*normF
      if(j==1){
        cov = cov1
      }else{
        cov <- cov+cov1  
      }
    }
    cov <- round(cov/j,2)
    rm(cov1,gr)
    gr <- as(cov,'GRanges')
    gr <- gr[gr$score!=0]
    gr <- sort(gr)
    gr <- subset(gr,seqnames(gr)%in%c(chr.flt))
    gr <- keepSeqlevels(gr, c(chr.flt))
    gr <- as.data.frame(gr)
    gr$strand <- NULL
    gr$width <- NULL
    cat("written file: ",paste0("A:/work/WinstonLab/Natalia_Reim/MayNov17/vis/",names(res)[i],".merged.bedGraph"), "\n")
    write.table(gr,file = paste0("A:/work/WinstonLab/Natalia_Reim/MayNov17/vis/",names(res)[i],".merged.bedGraph"),sep = "\t",quote = F,row.names = F,col.names = F )
    rm(cov,gr)
  }
  
}







############################################################################################################
### DESeq2 older method that have been used by people in the literature:
if(T){
  
  rm(list=ls()) 
  load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenes.RData")
  genes <- subset(genes,genes$verification=="Verified" & genes$species=="Scer" & genes$chr!="chrM")
  
  load("A:/work/WinstonLab/Natalia_Reim/DataFrames/RawCounts_ChIPSignal_genesSGD.RData")
  load("A:/work/yeast_Annotations/refs.RData")
  genes <- merge(genes,unique(refs[,c("tracking_id","gene","link")]),by="tracking_id",all.x=T)
  rownames(cf) <- cf$tracking_id
  
  res <- read.delim("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Normalization.factors.NataliaMayNov2017_revised.txt",header=T)
  
  m0 <- cf[cf$species=="Scer",16:79]
  mr <- cf[cf$species=="Spom",16:79]
  
  normF <- res$normF
  names(normF) <- res$sample
  m0 <- m0[,names(normF)]
  mr <- mr[,names(normF)]
  
  type=names(mr)
  condition <- gsub("_[1234]","",type)
  colData <- data.frame(condition,type)
  row.names(colData) <- colData$type
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(mr),colData = colData,design = ~ condition)
  sizeF <- estimateSizeFactors(dds)
  normF = sizeF@colData@listData$sizeFactor
  m0 <- m0[,names(normF)]
  mr <- mr[,names(normF)]
  
  mydeseq2fun <- function(m,type,sizeF){
    require(DESeq2)
    condition <- gsub("_[1234]","",type)
    colData <- data.frame(condition,type)
    row.names(colData) <- colData$type
    dds <- DESeqDataSetFromMatrix(countData = as.matrix(m),colData = colData,design = ~ condition)
    sizeFactors(dds) <- sizeF
    dds <- DESeq(dds)
    res <- as.data.frame(results(dds))
    z <- round(as.data.frame(DESeq2::counts(dds,normalized=T)),2)
    res$tracking_id <- rownames(m)
    res <- cbind(res,z)
    res <- merge(res,genes,by="tracking_id",all.x=T)
    res <- res[,-c(4,5,6)]
    names(res)[4] <- "FDR"
    j = names(dds@rowRanges@elementMetadata@listData)[12]
    j <- gsub("condition_","",j)
    j = unlist(strsplit(j,"_vs_"))
    res$log2FoldChange <- res$log2FoldChange*(-1)
    names(res)[3] <- paste0("log2 (",j[2],"/",j[1],")")
    return(res)
  }
  #rownames(m0) <- cf$tracking_id
  #rownames(mr) <- cf$tracking_id
  l <- list()
  ln <- list()
  for(i in seq(1,45,4)){
    cat(names(m0)[i:(i+3)],"\n")
    cat(names(normF[i:(i+3)]),"\n")
    cat(normF[i:(i+3)],"\n")
    cat(paste(as.character(res$Factor[i])),"\n")
    l[[paste(as.character(res$Factor[i]))]] <- mydeseq2fun(m = round(mr[,i:(i+3)]),type=names(mr)[i:(i+3)],sizeF= normF[i:(i+3)])  #normF[i:(i+3)]
    ln[[paste(as.character(res$Factor[i]))]] <- mydeseq2fun(m = round(m0[,i:(i+3)]),type=names(m0)[i:(i+3)],sizeF=normF[i:(i+3)])  #normF[i:(i+3)]
    #m <- mydeseq2fun(m = round(m0[,i:(i+3)]),type=names(m0)[i:(i+3)],sizeF=normF[i:(i+3)])  #normF[i:(i+3)]
  }
  save(l,ln,file="A:/work/WinstonLab/Natalia_Reim/DataFrames/Diff_Factor_Occupancy_DESeq2SGD_published.RData")
  
  
  load(file="A:/work/WinstonLab/Natalia_Reim/DataFrames/Diff_Factor_Occupancy_DESeq2SGD_published.RData")
  pl <- list()
  pln <- list()
  for(i in 1:length(l)){
    pl[[paste(names(l)[i])]] <- l[[i]][,3]
    pln[[paste(names(ln)[i])]] <- ln[[i]][,3]
  }
  boxplot(pl,ylim=c(-3,3))
  boxplot(pln,ylim=c(-3,2))
  abline(h=0)
  
}
############################################################################################################
## Log2FC comparison between ChIP-seq and RNA-seq with DESeq2 (old; can be deleted)
if(F){
  rm(list=ls())
  load("Gene groups.RData")
  load("A:/work/WinstonLab/Natalia_Reim/DataFrames/Log2FC values for the ChIPFactors_DESeq2.RData")
  m.spike$tracking_id <- rownames(m.spike)
  load(file="A:/work/WinstonLab/Natalia_Reim/Analysis_Feb18/DE_Genes_Filtered.RData")
  
  m.spike <- merge(m.spike,de.genes[,c("tracking_id","log2FoldChange","padj")],by="tracking_id",all.x=T)
  
  myplotfun <- function(m.spike,i,j){
    plot(m.spike[,i],m.spike[,j],main=names(m.spike)[i],xlab="log2 FC, factor occupancy", ylab="log2 FC,RNA",cex.lab=1.5,cex.axis=1.5,cex=1.2,cex.main=1.5,xaxt='n')
    axis(1,at =seq(-5,5,0.5),labels = seq(-5,5,0.5),tick = T ,cex=1.5,cex.axis=1.5)
    abline(h=seq(-5,5,1),v = seq(-5,5,1),col="gray",lty=2)
    mtext(paste("Pearson correlation=",round(cor(m.spike[,i],m.spike[,j],use="complete.obs"),3),sep=""),3,-1.5,col="red",cex=1.25)
    r <- lm(m.spike[,j]~m.spike[,i],na.action = na.omit)
    abline(r,col="red",lty=1)
  }
  
  pdf("A:/work/WinstonLab/Natalia_Reim/DataFrames/FactorOccupancy and expression correlations_DESeq2.pdf",width = 12,height = 15)
  par(mfrow=c(4,3))
  for(i in 2:13){
    myplotfun(m.spike,i,14)
  }
  dev.off()
  
  pdf("A:/work/WinstonLab/Natalia_Reim/DataFrames/FactorOccupancy and expression correlations_RPgenes_DESeq2.pdf",width = 12,height = 15)
  par(mfrow=c(4,3))
  for(i in 2:13){
    myplotfun(m.spike[m.spike$tracking_id%in%rp$tracking_id,],i,14)
  }
  dev.off()
  
  de.genes <- de.genes[de.genes$padj<0.05 & de.genes$log2FoldChange<0,]
  pdf("A:/work/WinstonLab/Natalia_Reim/DataFrames/FactorOccupancy and expression correlations_downgenes_DESeq2.pdf",width = 12,height = 15)
  par(mfrow=c(4,3))
  for(i in 2:13){
    myplotfun(m.spike[m.spike$tracking_id%in%de.genes$tracking_id,],i,14)
  }
  dev.off()
  
  m.spike$padj <- -log10(m.spike$padj)
  m.spike$padj <- ifelse(is.na(m.spike$padj),0,m.spike$padj)
  m.spike$padj <- ifelse(m.spike$padj>20,20,m.spike$padj)
  myplotfun <- function(m.spike,i,j){
    plot(m.spike[,i],m.spike[,j],main=names(m.spike)[i],xlab="log2 FC, factor occupancy", ylab="Sensitivity to expression change",cex.lab=1.5,cex.axis=1.5,cex=1.2,cex.main=1.5,xaxt='n')
    axis(1,at =seq(-5,5,1),labels = seq(-5,5,1),tick = T ,cex=1.5,cex.axis=1.5)
    abline(h=seq(-5,5,1),v = seq(-5,5,1),col="gray",lty=2)
    mtext(paste("r=",round(cor(m.spike[,i],m.spike[,j],use="complete.obs"),3),sep=""),3,-1.5,col="red",cex=1.25)
    r <- lm(m.spike[,j]~m.spike[,i],na.action = na.omit)
    abline(r,col="red",lty=1)
  }
  
  pdf("A:/work/WinstonLab/Natalia_Reim/DataFrames/FactorOccupancy and expression sensitivity correlations_DESeq2.pdf",width = 12,height = 15)
  par(mfrow=c(4,3))
  for(i in 2:13){
    myplotfun(m.spike[m.spike$tracking_id%in%de.genes$tracking_id,],i,15)
  }
  dev.off()
  
}






