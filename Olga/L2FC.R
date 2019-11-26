rm(list=ls())
library(DESeq2)
library(GenomicRanges)
library(ggplot2)
chr.flt <- c("chrI","chrII","chrIII","chrIV","chrIX","chrV","chrVI","chrVII","chrVIII","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
spiked.chr <- c("chrXVII","chrXVIII","chrXIX")
mydeseq2fun <- function(m,sizeF){
  require(DESeq2)
  type=names(m)
  condition <- gsub("_[01234]","",type)
  condition <- gsub("_I|_II","",condition)
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
  DESeq2::plotMA(dds,main=paste0("log2 (",j[2],"/",j[1],")"))
  
  res$log2FoldChange <- res$log2FoldChange*(-1)
  names(res)[3] <- paste0("log2 (",j[2],"/",j[1],")")
  return(res)
}

##################################################################################
#### TSS-seq data
##################################################################################
OUTDIR="A:/work/WinstonLab/Olga/TSS_Seq/Steinmetz/"
## DE of peaks, promoter regions etc
if(T){
  load(paste0(OUTDIR,"TSS-SeqPeaks_AdaptiveIDR_SteinmetzAnno.RData"))
  g <- as(res,"GRanges")
  g <- reduce(g)
  range(width(g))
  
  ## annotae peaks
  load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenesFull.RData")
  genes <- genes[!is.na(genes$verification),]
  genes <- as(genes,"GRanges")
  prom <- flank(genes,width = 200)
  g <- my.overlap.annotation.function(genes,prom,g[,0])
  length(res)
  
  setwd("A:/work/WinstonLab/Olga/TSS_Seq/gr")
  files <- Sys.glob("*.gr")
  files <- files[-9] ## remove bad replicate
  
  ## counts along unified peaks
  ct <- c()
  for(f in files){
    cat(f,"\n")
    load(f)
    gr <- resize(gr,1,"start")
    ct <- cbind(ct,countOverlaps(g,gr,ignore.strand=F))
  }
  colnames(ct) <- gsub(".gr","",files)
  df1 <- cbind(as.data.frame(g),as.data.frame(ct))
  rm(ct,gr,g)
  save(df1,file=paste0(OUTDIR,"readcounts along unifiedPeaks.RData"))
  
  ## counts along promoters
  prom <- subset(prom,seqnames(prom)%in%c(spiked.chr,chr.flt))
  prom <- resize(prom,width(prom)+50,fix = "start")
  cnt <- count.grfilelist.genomicBins(files,prom,filtgr = NULL,gwidth = 1,fix = "start",ignore.strand = F)
  cntas <- count.grfilelist.genomicBins(files,prom,filtgr = NULL,gwidth = 1,fix = "start",ignore.strand = T)
  cntas[,19:30] <- cntas[,19:30] - cnt[,19:30]
  save(cnt,cntas,file=paste0(OUTDIR,"readcounts along promoters.RData"))
  
  ## spikein normalization factors
  require(DESeq2)
  w <- cnt[cnt$seqnames%in%spiked.chr,]
  w <- w[,c("WT_1", "WT_2", "spt6-1004_1", "spt6-1004_2", "spt6_1", "spt6_2", 
            "ts_Spt6_0", "ts_Spt6_1", "ts_WT_1", "ts_WT_2", "ts_spt6-1004_1", 
            "ts_spt6-1004_2")]
  w <- w[rowSums(w)>0,]
  type <- names(w)
  condition <- gsub("_[1230]$","",type)
  colData <- data.frame(condition,type)
  row.names(colData) <- colData$type
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(w),colData = colData,design = ~ condition)
  dds <- DESeq(dds)
  dd <- sizeFactors(dds)
  save(dd,file=paste0(OUTDIR,"TSS_Seq_SizeFactors.spikeinNormalization.RData"))
}

if(F){
  setwd("A:/work/WinstonLab/Olga/TSS_Seq/Steinmetz")
  #rm(list=ls())
  OUTDIR="A:/work/WinstonLab/Olga/TSS_Seq/Steinmetz/"
  load(paste0(OUTDIR,"readcounts along unifiedPeaks.RData"))
  load(paste0(OUTDIR,"TSS_Seq_SizeFactors.spikeinNormalization.RData"))
  load(paste0(OUTDIR,"readcounts along promoters.RData"))
  
  load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenes.RData")
  #genes <- subset(genes,genes$verification=="Verified" & genes$species=="Scer" & genes$chr!="chrM")
  load("A:/work/yeast_Annotations/refs.RData")
  genes <- merge(genes,unique(refs[,c("tracking_id","gene","link")]),by="tracking_id",all.x=T)
  
  rownames(cnt) <- cnt$tracking_id
  rownames(cntas) <- cntas$tracking_id
  
  l <- list()
  ######### Sense unified peaks pre-post spikein 
  pdf(file=paste0(OUTDIR,"PromoterSense_PostSpike_MAplots.pdf"),width = 8,height = 8)
  par(mfrow=c(2,2))
  l[["spt61004_30C"]] <- mydeseq2fun(cnt[,c(19:22)],sizeF=dd[c(1:4)])
  l[["spt6-YW_30C"]] <- mydeseq2fun(cnt[,c(19:20,23:24)],sizeF=dd[c(1,2,5,6)])
  l[["spt6-YW_37C"]] <- mydeseq2fun(cnt[,c(27:28,25:26)],sizeF=dd[c(9:10,7:8)])
  l[["spt61004_37C"]] <- mydeseq2fun(cnt[,c(27:30)],sizeF=dd[c(9:12)])
  dev.off()

  save(l,file="A:/work/WinstonLab/Olga/TSS_Seq/Steinmetz/DE_promoters_PostSpikein_L2FC.RData")
  
  
  ######### Sense unified peaks pre-post spikein 
  pdf(file=paste0(OUTDIR,"SenseUnifiedPeaksExpression_MAplots.pdf"),width = 8,height = 8)
  par(mfrow=c(2,2))
  spt6YW <- mydeseqfun(df1,ff = c(21,22,25,26))
  tsspt6YW <- mydeseqfun(df1,c(29,30,27,28))
  spt61004 <- mydeseqfun(df1,21:24)
  tsspt61004 <- mydeseqfun(df1,c(29:32))
  dev.off()
  
  pdf(file=paste0(OUTDIR,"SenseUnifiedPeaksExpression_MAplots_postspike.pdf"),width = 8,height = 8)
  par(mfrow=c(2,2))
  spt6YWsp <- mydeseqfun(df1,c(21,22,25,26),spike.norm = T,dd=dd[c(1,2,5,6)])
  tsspt6YWsp <- mydeseqfun(df1,c(29,30,27,28),spike.norm = T,dd =dd[c(9,10,7,8)] )
  spt61004sp <- mydeseqfun(df1,21:24,spike.norm = T,dd =dd[c(1:4)])
  tsspt61004sp <- mydeseqfun(df1,c(29:32),spike.norm = T,dd =dd[c(9:12)])
  dev.off()
  save.image(file="A:/work/WinstonLab/Olga/TSS_Seq/Steinmetz/DE UNified peaks.RData")
  
  
  
  
  rm(list=ls())
  load(file="A:/work/WinstonLab/Olga/TSS_Seq/Steinmetz/DE UNified peaks.RData")
  library(xlsx)
  write.xlsx(spt6YW[spt6YW$padj<0.05 & !is.na(spt6YW$padj),-c(23:25)],file="Log2FC Peak report.xlsx",sheetName = "spt6-YW @30°C",row.names = F,showNA = T)
  write.xlsx(tsspt6YW[tsspt6YW$padj<0.05 & !is.na(tsspt6YW$padj),-c(23:25)],file="Log2FC Peak report.xlsx",sheetName = "spt6-YW @37°C",row.names = F,showNA = T,append = T)
  write.xlsx(spt61004[spt61004$padj<0.05 & !is.na(spt61004$padj),-c(23:25)],file="Log2FC Peak report.xlsx",sheetName = "spt6-1004 @30°C",row.names = F,showNA = T,append = T)
  write.xlsx(tsspt61004[tsspt61004$padj<0.05 & !is.na(tsspt61004$padj),-c(16:18)],file="Log2FC Peak report.xlsx",sheetName = "spt6-1004 @37°C",row.names = F,showNA = T,append = T)
  
  write.table(spt6YWsp[spt6YWsp$padj<0.05 & !is.na(spt6YWsp$padj),-c(23:25)],file="temp.txt",quote = F,sep = "\t",row.names = F)
  write.table(tsspt6YWsp[tsspt6YWsp$padj<0.05 & !is.na(tsspt6YWsp$padj),-c(23:25)],file="temp.txt",quote = F,sep = "\t",row.names = F)
  write.table(spt61004sp[spt61004sp$padj<0.05 & !is.na(spt61004sp$padj),-c(23:25)],file="temp.txt",quote = F,sep = "\t",row.names = F)
  write.table(tsspt61004sp[tsspt61004sp$padj<0.05 & !is.na(tsspt61004sp$padj),-c(16:18)],file="temp.txt",quote = F,sep = "\t",row.names = F)
  
  
  ## heatmaps
  library(UpSetR)
  library(gplots)
  m <- data.frame(ifelse(spt6YW[,"padj"]<0.05 & !is.na(spt6YW[,"padj"]),spt6YW[,22],NA),
                  ifelse(tsspt6YW[,"padj"]<0.05 & !is.na(tsspt6YW[,"padj"]),tsspt6YW[,22],NA),
                  ifelse(spt61004[,"padj"]<0.05 & !is.na(spt61004[,"padj"]),spt61004[,22],NA),
                  ifelse(tsspt61004[,"padj"]<0.05 & !is.na(tsspt61004[,"padj"]),tsspt61004[,22],NA))
  names(m) <- c("spt6-YW @30°C","spt6-YW @37°C","spt6-1004 @30°C","spt6-1004 @37°C")
  rownames(m) <- paste0(df1$tracking_id,"_",df1$Category,"_",df1$start)
  m <- m[apply(m, 1, function(y) !all(is.na(y))),]
  z <- df1[,c(21:32)]
  z1 <- as.data.frame(sapply(seq(1, 11,2), function(x) rowSums(z[, c(x,x+1)])))
  names(z1) <- names(z)[seq(1,11,2)]
  names(z1) <- gsub("_[10]","",names(z1))
  z1 <- apply(z1,2,function(x) round((x*1e6)/sum(x),2))
  z1 <- as.data.frame(z1)
  z1$width <- df1$width
  z1[,1:6] <- z1[,1:6]*1000/z1[,7]
  z1$width <- NULL
  rownames(z1) <- paste0(df1$tracking_id,"_",df1$Category,"_",df1$start)
  z1 <- z1[rowSums(z1)>10,]
  z1 <- as.data.frame(t(apply(z1,1,function(x){
    round((x-mean(x))/sd(x),3)  
  })))
  names(z1) <- c("WT @30°C","spt6-1004 @30°C","spt6-YW @30°C","spt6-YW @37°C","WT @37°C","spt6-1004 @37°C")
  
  
  pdf(file=paste0(OUTDIR,"Heatplots related to peak expression.pdf"),width = 8,height = 8)
  for(qq in c("Promoter_S", "Promoter_AS", "GeneBody_S", "GeneBody_AS", "Intergenic")){
    cat(qq,"\n")
    ## plot1
    hmcols<-(colorRampPalette(c("yellow","orange","red4"))(256))
    z2 <- z1[grep(qq,rownames(z1)),]
    heatmap.2(as.matrix(z2),main=paste0("expression (",qq,")"),dendrogram = "column",key=T,density.info="none",key.title="row Z-score",trace="none",cexRow = 0.45,labRow = NA,symm = F,margins = c(6,16),cexCol = 0.7,col = hmcols,scale="none")
    ## plot2
    z2 <- subset(z1,rownames(z1)%in%rownames(m))
    z2 <- z2[grep(qq,rownames(z2)),]
    hmcols<-(colorRampPalette(c("yellow","orange","red4"))(256))
    x <- rownames(z2)
    x[-seq(1,length(x),200)] <- NA
    heatmap.2(as.matrix(z2),main=paste0("expression, DEpeaks (",qq,")"),dendrogram = "column",key=T,density.info="none",key.title="row Z-score",trace="none",cexRow = 0.45,labRow = NA,symm = F,margins = c(6,16),cexCol = 0.7,col = hmcols,scale="none")
    ## plot3
    m[is.na(m)] <- 0
    m <- m*(-1)
    z2 <- m[grep(qq,rownames(m)),]
    hmcols<-(colorRampPalette(c("blue4","white","red4"))(256))
    heatmap.2(as.matrix(m),key=T,main=paste0("log2FC (",qq,")"),dendrogram = "column",density.info="none",key.title="row Z-score",trace="none",cexRow = 0.45,labRow = NA,symm = F,margins = c(6,15),cexCol = 0.7,col = hmcols,scale="none")
    rm(z2)
  }
  dev.off()
  
  
  m <- data.frame(ifelse(spt6YW[,"padj"]<0.05 & !is.na(spt6YW[,"padj"]),spt6YW[,22],NA),
                  ifelse(tsspt6YW[,"padj"]<0.05 & !is.na(tsspt6YW[,"padj"]),tsspt6YW[,22],NA),
                  ifelse(spt61004[,"padj"]<0.05 & !is.na(spt61004[,"padj"]),spt61004[,22],NA),
                  ifelse(tsspt61004[,"padj"]<0.05 & !is.na(tsspt61004[,"padj"]),tsspt61004[,22],NA))
  names(m) <- c("spt6-YW @30°C","spt6-YW @37°C","spt6-1004 @30°C","spt6-1004 @37°C")
  rownames(m) <- paste0(df1$tracking_id,"_",df1$Category,"_",df1$start)
  m <- m[apply(m, 1, function(y) !all(is.na(y))),]
  m[is.na(m)]<- 0
  m <- m*(-1)
  mup <- mdown <- m
  mup <- ifelse(mup <= 0,0,1)
  mdown <- ifelse(mdown >= 0,0,1)
  colnames(mup) <- paste(colnames(mup),"(U)",sep=":")
  colnames(mdown) <- paste(colnames(mdown),"(D)",sep=":")
  mup<- as.data.frame(mup)
  mdown <- as.data.frame(mdown)
  m <- cbind(mup,mdown)
  pdf(fle=paste0(OUTDIR,"PeakOlaps_UpSetR.pdf"),width = 10,height = 6)
  print(upset(m, sets = names(m), mb.ratio = c(0.55, 0.45), order.by = "freq"))
  print(upset(mup, sets = names(mup), mb.ratio = c(0.55, 0.45), order.by = "freq"))
  print(upset(mdown, sets = names(mdown), mb.ratio = c(0.55, 0.45), order.by = "freq"))
  dev.off()
  
  
  ########################################################################################################### 
  m <- data.frame(ifelse(spt6YWsp[,"padj"]<0.05 & !is.na(spt6YWsp[,"padj"]),spt6YWsp[,22],NA),
                  ifelse(tsspt6YWsp[,"padj"]<0.05 & !is.na(tsspt6YWsp[,"padj"]),tsspt6YWsp[,22],NA),
                  ifelse(spt61004sp[,"padj"]<0.05 & !is.na(spt61004sp[,"padj"]),spt61004sp[,22],NA),
                  ifelse(tsspt61004sp[,"padj"]<0.05 & !is.na(tsspt61004sp[,"padj"]),tsspt61004sp[,22],NA))
  names(m) <- c("spt6-YW @30°C","spt6-YW @37°C","spt6-1004 @30°C","spt6-1004 @37°C")
  rownames(m) <- paste0(df1$tracking_id,"_",df1$Category,"_",df1$start)
  m <- m[apply(m, 1, function(y) !all(is.na(y))),]
  z <- df1[,c(21:32)]
  z <- t(t(z)/(dd)) # spikein normalization
  z1 <- as.data.frame(sapply(seq(1, 11,2), function(x) rowSums(z[, c(x,x+1)])))
  names(z1) <- names(z)[seq(1,11,2)]
  names(z1) <- gsub("_[10]","",names(z1))
  z1 <- as.data.frame(z1)
  z1$width <- df1$width
  z1[,1:6] <- z1[,1:6]*1000/z1[,7]
  z1$width <- NULL
  rownames(z1) <- paste0(df1$tracking_id,"_",df1$Category,"_",df1$start)
  z1 <- z1[rowSums(z1)>10,]
  z1 <- as.data.frame(t(apply(z1,1,function(x){
    round((x-mean(x))/sd(x),3)  
  })))
  names(z1) <- c("WT @30°C","spt6-1004 @30°C","spt6-YW @30°C","spt6-YW @37°C","WT @37°C","spt6-1004 @37°C")
  
  
  pdf(file=paste0(OUTDIR,"Heatplots related to peak expression_postSpikein.pdf"),width = 8,height = 8)
  for(qq in c("Promoter_S", "Promoter_AS", "GeneBody_S", "GeneBody_AS", "Intergenic")){
    cat(qq,"\n")
    ## plot1
    hmcols<-(colorRampPalette(c("yellow","orange","red4"))(256))
    z2 <- z1[grep(qq,rownames(z1)),]
    heatmap.2(as.matrix(z2),main=paste0("expression (",qq,")"),dendrogram = "column",key=T,density.info="none",key.title="row Z-score",trace="none",cexRow = 0.45,labRow = NA,symm = F,margins = c(6,16),cexCol = 0.7,col = hmcols,scale="none")
    ## plot2
    z2 <- subset(z1,rownames(z1)%in%rownames(m))
    z2 <- z2[grep(qq,rownames(z2)),]
    hmcols<-(colorRampPalette(c("yellow","orange","red4"))(256))
    x <- rownames(z2)
    x[-seq(1,length(x),200)] <- NA
    heatmap.2(as.matrix(z2),main=paste0("expression, DEpeaks (",qq,")"),dendrogram = "column",key=T,density.info="none",key.title="row Z-score",trace="none",cexRow = 0.45,labRow = NA,symm = F,margins = c(6,16),cexCol = 0.7,col = hmcols,scale="none")
    ## plot3
    m[is.na(m)] <- 0
    m <- m*(-1)
    z2 <- m[grep(qq,rownames(m)),]
    hmcols<-(colorRampPalette(c("blue4","white","red4"))(256))
    heatmap.2(as.matrix(m),key=T,main=paste0("log2FC (",qq,")"),dendrogram = "column",density.info="none",key.title="row Z-score",trace="none",cexRow = 0.45,labRow = NA,symm = F,margins = c(6,15),cexCol = 0.7,col = hmcols,scale="none")
    rm(z2)
  }
  dev.off()
  
  
  m <- data.frame(ifelse(spt6YWsp[,"padj"]<0.05 & !is.na(spt6YWsp[,"padj"]),spt6YWsp[,22],NA),
                  ifelse(tsspt6YWsp[,"padj"]<0.05 & !is.na(tsspt6YWsp[,"padj"]),tsspt6YWsp[,22],NA),
                  ifelse(spt61004sp[,"padj"]<0.05 & !is.na(spt61004sp[,"padj"]),spt61004sp[,22],NA),
                  ifelse(tsspt61004sp[,"padj"]<0.05 & !is.na(tsspt61004sp[,"padj"]),tsspt61004sp[,22],NA))
  names(m) <- c("spt6-YW @30°C","spt6-YW @37°C","spt6-1004 @30°C","spt6-1004 @37°C")
  rownames(m) <- paste0(df1$tracking_id,"_",df1$Category,"_",df1$start)
  m <- m[apply(m, 1, function(y) !all(is.na(y))),]
  m[is.na(m)]<- 0
  m <- m*(-1)
  mup <- mdown <- m
  mup <- ifelse(mup <= 0,0,1)
  mdown <- ifelse(mdown >= 0,0,1)
  colnames(mup) <- paste(colnames(mup),"(U)",sep=":")
  colnames(mdown) <- paste(colnames(mdown),"(D)",sep=":")
  mup<- as.data.frame(mup)
  mdown <- as.data.frame(mdown)
  m <- cbind(mup,mdown)
  pdf(file=paste0(OUTDIR,"PeakOlaps_UpSetR_postSpikein.pdf"),width = 10,height = 6)
  print(upset(m, sets = names(m), mb.ratio = c(0.55, 0.45), order.by = "freq"))
  print(upset(mup, sets = names(mup), mb.ratio = c(0.55, 0.45), order.by = "freq"))
  print(upset(mdown, sets = names(mdown), mb.ratio = c(0.55, 0.45), order.by = "freq"))
  dev.off()
  
}


##################################################################################
#### NET-seq data
##################################################################################
OUTDIR="A:/work/WinstonLab/Olga/Netseq/"
if(F){
  library(GenomicRanges)
  setwd("A:/work/WinstonLab/Olga/Netseq")
  load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenes.RData")
  #genes <- subset(genes,genes$verification=="Verified" & genes$species=="Scer" & genes$chr!="chrM")
  genes <- as(genes,"GRanges")
  files <- Sys.glob("*.gr")
  s <- count.grfilelist.genomicBins_covfiles(grfiles = files,genomicBins = genes,ignore.strand = F)
  as <- count.grfilelist.genomicBins_covfiles(grfiles = files,genomicBins = genes,ignore.strand = T)
  as[,16:ncol(as)] <- as[,16:ncol(as)]-s[,16:ncol(s)]
  
  genes <- resize(genes,500,fix = "start")
  s500 <- count.grfilelist.genomicBins(grfiles = files,genomicBins = genes,gwidth = NULL,ignore.strand = F)
  as500 <- count.grfilelist.genomicBins(grfiles = files,genomicBins = genes,gwidth = NULL,ignore.strand = T)
  as500[,16:ncol(as)] <- as500[,16:ncol(as)]-s500[,16:ncol(s)]
  save(s,as,s500,as500,file=paste0(OUTDIR,"/NETSeq_countspergenes_CountsPer500bp.RData"))

  
  ######### Sense unified peaks pre-post spikein 
  load(paste0(OUTDIR,"/NETSeq_countspergenes_CountsPer500bp.RData"))
  load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenes.RData")
  genes <- subset(genes,genes$verification=="Verified" & genes$species=="Scer" & genes$chr!="chrM")
  s <- s[s$tracking_id%in%genes$tracking_id,]
  as <- as[as$tracking_id%in%genes$tracking_id,]
  rownames(s) <- s$tracking_id
  rownames(as) <- as$tracking_id
  s[,16:27] <- round(s[,16:27])
  as[,16:27] <- round(as[,16:27])
  
  l <- list()
  pdf(file=paste0(OUTDIR,"PromoterSense_MAplots.pdf"),width = 8,height = 8)
  par(mfrow=c(2,2))
  l[["spt6-YW_30C"]] <- mydeseq2fun(s[,c(24,25,16,17)],sizeF = c(1,1,1,1))
  l[["spt61004_30C"]] <- mydeseq2fun(s[,c(24,25,20,21)],sizeF = c(1,1,1,1))
  l[["spt6-YW_37C"]] <- mydeseq2fun(s[,c(26,27,18,19)],sizeF = c(1,1,1,1))
  l[["spt61004_37C"]] <- mydeseq2fun(s[,c(26,27,22,23)],sizeF = c(1,1,1,1))
  dev.off()
  
  save(l,file=paste0(OUTDIR,"DE_promoters_PostSpikein_L2FC.RData"))

}



##################################################################################
#### ChIP-seq data
##################################################################################
OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/"
## Counts over genes
if(F){
  get_counts <- function(OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/"){
    setwd(OUTDIR)
    load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenes.RData")
    genes <- subset(genes,genes$verification=="Verified"  & genes$chr!="chrM" & !is.na(genes$verification))
    genes <- as(genes,"GRanges")
    
    res <- read.delim(paste0(OUTDIR,"NormalizationFactors.txt"),header = T,stringsAsFactors = F)
    res <- res[res$Factor !="Inp",]
    res$file <- paste0(res$sample,".gr")
    
    libsize=1e7
    count <- c()
    countsc <- c()
    counta <- c()
    for(f in res$file){
      cat(f,"\n")
      load(f)
      ## raw counts
      gr <- subset(gr,seqnames(gr)%in%c(chr.flt,spiked.chr))
      gr <- keepSeqlevels(gr,c(chr.flt,spiked.chr))
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
      #x <- libsize/length(gr)
      #x <- x/res[res$file==f,]$normF
      #x <- round(x*length(gr))
      #cat(" div", length(gr),' ',x," ",x/length(gr) ,"\n")
      #rn <- sample(1:length(gr),size = x,replace = T)
      #gr <- gr[rn]
      #countb <- cbind(countb,countOverlaps(genes,gr,ignore.strand=T))
      
      rn <- sample(1:length(gr),size = libsize,replace = T)
      gr <- gr[rn]
      countsc <- cbind(countsc,countOverlaps(genes,gr,ignore.strand=T))
      rm(gr,f,x,rn)
    }
    colnames(countsc) <- colnames(counta) <- colnames(count) <- gsub(".gr","",res$file)
    cf <- cbind(as.data.frame(genes),as.data.frame(count))
    save(cf,file=paste0(OUTDIR,"RawCounts_Signal_genesSGD.RData"))
    cf <- cbind(as.data.frame(genes),countsc)
    save(cf,file=paste0(OUTDIR,"LibSizeScaledCounts_Signal_genesSGD.RData"))
    cf <- cbind(as.data.frame(genes),counta)
    save(cf,file=paste0(OUTDIR,"SpikeinScaled_MultipCounts_Signal_genesSGD.RData"))
    return(1)
  }
  
  get_counts(OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/")
  get_counts(OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/")
  
}

## log2FC calculations
if(T){
  mydeseq2fun <- function(m,type,sizeF){
    require(DESeq2)
    condition <- gsub("_r[1234]","",type)
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
  
  ### (1) Spiken normalization
  l <- list()
  ln <- list()
  
  
  OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/"
  load(paste0(OUTDIR,"SpikeinScaled_MultipCounts_Signal_genesSGD.RData"))
  rownames(cf) <- cf$tracking_id
  names(cf) <- gsub("MycChIP","Spt16",names(cf))
  names(cf) <- gsub("8WG16ChIP","Pol2-July",names(cf))
  names(cf) <- gsub("wt","WT",names(cf))
  names(cf) <- gsub("spt6-YW-pob3","double",names(cf))
  m0 <- cf[cf$species=="Scer",16:ncol(cf)]
  mr <- cf[cf$species=="Spom",16:ncol(cf)]
  genotypes <- c("pob3","spt6-YW", "double") #,   "spt6-YW_Flag","WT_Flag",  "pob3_Flag",     "spt6-YW-pob3_Flag")
  chips <- c("Pol2-July", "Spt16") #,   "spt6-YW_Flag","WT_Flag",  "pob3_Flag",     "spt6-YW-pob3_Flag")
  for(tmp in c(30,37)){
    for(chip in chips){
      for(gt in genotypes){
        cntr <- names(m0)[grepl("WT_",names(m0),fixed = T) & grepl(chip,names(m0),fixed = T) & grepl(tmp,names(m0),fixed = T)] 
        exp <- names(m0)[grepl(gt,names(m0),fixed = T) & grepl(chip,names(m0),fixed = T) & grepl(chip,names(m0),fixed = T) & grepl(tmp,names(m0),fixed = T)] 
        cat(cntr,"\n",exp,"\n\n")
        z0 <- m0[,c(cntr,exp)]
        zr <- mr[,c(cntr,exp)]
        l[[paste(gt,chip,tmp,sep = "_")]] <- mydeseq2fun(m = zr,type=names(zr),sizeF= c(1,1,1,1))  #normF[i:(i+3)]
        ln[[paste(gt,chip,tmp,sep = "_")]] <- mydeseq2fun(m = z0,type=names(z0),sizeF=c(1,1,1,1))  #normF[i:(i+3)]
      }
    }
  }
  rm(mr,m0,cf,z0,zr,chip,tmp,gt)
  
  OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/"
  load(paste0(OUTDIR,"SpikeinScaled_MultipCounts_Signal_genesSGD.RData"))
  rownames(cf) <- cf$tracking_id
  names(cf) <- gsub("FlagChIP","Spt6",names(cf))
  names(cf) <- gsub("8WG16ChIP","Pol2-Aug",names(cf))
  names(cf) <- gsub("wt","WT",names(cf))
  names(cf) <- gsub("spt6-YW-pob3","double",names(cf))
  m0 <- cf[cf$species=="Scer",16:ncol(cf)]
  mr <- cf[cf$species=="Spom",16:ncol(cf)]
  genotypes <- c("pob3","spt6-YW", "double") #,   "spt6-YW_Flag","WT_Flag",  "pob3_Flag",     "spt6-YW-pob3_Flag")
  chips <- c("Pol2-Aug", "Spt6") #,   "spt6-YW_Flag","WT_Flag",  "pob3_Flag",     "spt6-YW-pob3_Flag")
  for(tmp in c(30,37)){
    for(chip in chips){
      for(gt in genotypes){
        cntr <- names(m0)[grepl("WT_",names(m0),fixed = T) & grepl(chip,names(m0),fixed = T) & grepl(tmp,names(m0),fixed = T)] 
        exp <- names(m0)[grepl(gt,names(m0),fixed = T) & grepl(chip,names(m0),fixed = T) & grepl(chip,names(m0),fixed = T) & grepl(tmp,names(m0),fixed = T)] 
        cat(cntr,"\n",exp,"\n\n")
        z0 <- m0[,c(cntr,exp)]
        zr <- mr[,c(cntr,exp)]
        l[[paste(gt,chip,tmp,sep = "_")]] <- mydeseq2fun(m = zr,type=names(zr),sizeF= c(1,1,1,1))  #normF[i:(i+3)]
        ln[[paste(gt,chip,tmp,sep = "_")]] <- mydeseq2fun(m = z0,type=names(z0),sizeF=c(1,1,1,1))  #normF[i:(i+3)]
      }
    }
  }
  rm(mr,m0,cf,z0,zr,chip,tmp,gt)
  
  ###(2) non spikein normalized
  w <- list()
  wn <- list()
  
  OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/"
  load(paste0(OUTDIR,"RawCounts_Signal_genesSGD.RData"))
  rownames(cf) <- cf$tracking_id
  names(cf) <- gsub("MycChIP","Spt16",names(cf))
  names(cf) <- gsub("8WG16ChIP","Pol2-July",names(cf))
  names(cf) <- gsub("wt","WT",names(cf))
  names(cf) <- gsub("spt6-YW-pob3","double",names(cf))
  m0 <- cf[cf$species=="Scer",16:ncol(cf)]
  mr <- cf[cf$species=="Spom",16:ncol(cf)]
  genotypes <- c("pob3","spt6-YW", "double") #,   "spt6-YW_Flag","WT_Flag",  "pob3_Flag",     "spt6-YW-pob3_Flag")
  chips <- c("Pol2-July", "Spt16") #,   "spt6-YW_Flag","WT_Flag",  "pob3_Flag",     "spt6-YW-pob3_Flag")
  for(tmp in c(30,37)){
    for(chip in chips){
      for(gt in genotypes){
        cntr <- names(m0)[grepl("WT_",names(m0),fixed = T) & grepl(chip,names(m0),fixed = T) & grepl(tmp,names(m0),fixed = T)] 
        exp <- names(m0)[grepl(gt,names(m0),fixed = T) & grepl(chip,names(m0),fixed = T) & grepl(chip,names(m0),fixed = T) & grepl(tmp,names(m0),fixed = T)] 
        cat(cntr,"\n",exp,"\n\n")
        z0 <- m0[,c(cntr,exp)]
        zr <- mr[,c(cntr,exp)]
        w[[paste(gt,chip,tmp,sep = "_")]] <- mydeseq2fun(m = zr,type=names(zr),sizeF= 0)  #normF[i:(i+3)]
        wn[[paste(gt,chip,tmp,sep = "_")]] <- mydeseq2fun(m = z0,type=names(z0),sizeF=0)  #normF[i:(i+3)]
      }
    }
  }
  rm(mr,m0,cf,z0,zr,chip,tmp,gt)
  
  
  OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/"
  load(paste0(OUTDIR,"RawCounts_Signal_genesSGD.RData"))
  rownames(cf) <- cf$tracking_id
  names(cf) <- gsub("FlagChIP","Spt6",names(cf))
  names(cf) <- gsub("8WG16ChIP","Pol2-Aug",names(cf))
  names(cf) <- gsub("wt","WT",names(cf))
  names(cf) <- gsub("spt6-YW-pob3","double",names(cf))
  m0 <- cf[cf$species=="Scer",16:ncol(cf)]
  mr <- cf[cf$species=="Spom",16:ncol(cf)]
  genotypes <- c("pob3","spt6-YW", "double") #,   "spt6-YW_Flag","WT_Flag",  "pob3_Flag",     "spt6-YW-pob3_Flag")
  chips <- c("Pol2-Aug", "Spt6") #,   "spt6-YW_Flag","WT_Flag",  "pob3_Flag",     "spt6-YW-pob3_Flag")
  for(tmp in c(30,37)){
    for(chip in chips){
      for(gt in genotypes){
        cntr <- names(m0)[grepl("WT_",names(m0),fixed = T) & grepl(chip,names(m0),fixed = T) & grepl(tmp,names(m0),fixed = T)] 
        exp <- names(m0)[grepl(gt,names(m0),fixed = T) & grepl(chip,names(m0),fixed = T) & grepl(chip,names(m0),fixed = T) & grepl(tmp,names(m0),fixed = T)] 
        cat(cntr,"\n",exp,"\n\n")
        z0 <- m0[,c(cntr,exp)]
        zr <- mr[,c(cntr,exp)]
        w[[paste(gt,chip,tmp,sep = "_")]] <- mydeseq2fun(m = zr,type=names(zr),sizeF= 0)  #normF[i:(i+3)]
        wn[[paste(gt,chip,tmp,sep = "_")]] <- mydeseq2fun(m = z0,type=names(z0),sizeF=0)  #normF[i:(i+3)]
      }
    }
  }
  rm(mr,m0,cf,z0,zr,chip,tmp,gt)
  save(l,ln,w,wn,file="A:/work/WinstonLab/Olga/ChIPSeq/Diff_Factor_Occupancy_DESeq2SGD.RData")
  
  
  
  ####### print sanity check boxplots for lo2FC values
  rm(list=ls())
  load(file="A:/work/WinstonLab/Olga/ChIPSeq/Diff_Factor_Occupancy_DESeq2SGD.RData")
  pl <- c()
  rl <- c()
  for(i in 1:length(l)){
    rl <- rbind(rl,
                data.frame(name=names(l)[i],
                           spikeinnorm_Spom=sum(l[[i]][,4] <0.05, na.rm = T ),
                           spikeinnorm_Scer=sum(ln[[i]][,4] <0.05, na.rm = T ),
                           LibNorm_Spom=sum(w[[i]][,4] <0.05, na.rm = T ),
                           LibNorm_Scer=sum(wn[[i]][,4] <0.05, na.rm = T ))
    )
    pl <- rbind(pl, data.frame(l2fc=l[[i]][,3],profile=names(l)[i],genotype="Spom",norm="Spikein"),
                    data.frame(l2fc=ln[[i]][,3],profile=names(l)[i],genotype="Scer",norm="Spikein"),
                    data.frame(l2fc=w[[i]][,3],profile=names(l)[i],genotype="Spom",norm="Library Size"),
                    data.frame(l2fc=wn[[i]][,3],profile=names(l)[i],genotype="Scer",norm="Library Size"))
  }
    q <- ggplot(pl,aes(x=profile,y=l2fc,fill=profile))+geom_boxplot()+
    facet_wrap(~genotype+norm,ncol = 2,scales = "free_x")+theme_bw()+
    xlab("")+
    theme(legend.position = "none",
          axis.text = element_text(size=14,color="black"),
          strip.text = element_text(size=16,color="black"))+
    coord_flip()+
    geom_hline(yintercept = 0,col="black")
  pdf(file="A:/work/WinstonLab/Olga/ChIPSeq/DESeq2_L2FC_boxplots_ScerSpom.pdf",width = 10,height = 12)
  print(q)
  grid.newpage()
  grid.table(rl)
  dev.off()
  
  
  ##### (2) Sanity check: RLE for Spom reads
  rm(list=ls())
  load(file="A:/work/WinstonLab/Olga/ChIPSeq/Diff_Factor_Occupancy_DESeq2SGD.RData")
  m.occup <- data.frame(tracking_id=w[[1]]$tracking_id)
  m.spike.occup <- data.frame(tracking_id=l[[1]]$tracking_id)
  for(i in 1:length(l)){
    cat(names(ln)[i],"\n")
    cat(" ",paste0(names(ln[[i]][,5:8]),collapse = ","),"\n")
    m.spike.occup <- cbind(m.spike.occup,l[[i]][,5:8])
    m.occup <- cbind(m.occup,w[[i]][,5:8])
    cat("\n\n")
  }
  rownames(m.occup) <- rownames(m.spike.occup) <- m.occup$tracking_id
  m.occup$tracking_id <- m.spike.occup$tracking_id <- NULL 
  m.occup <- m.occup[,unique(names(m.occup))]
  m.spike.occup <- m.spike.occup[,unique(names(m.spike.occup))]
  
  m.occup <- log2( (m.occup/apply(m.occup, 1, function(x) (median(x,na.rm=T)))))
  m.spike.occup <- log2( (m.spike.occup/apply(m.spike.occup, 1, function(x) (median(x,na.rm=T)))))
  
  chip="Pol2-July"
  tmp=30
  pl.list <- list()
  for (chip in c("Pol2-July","Pol2-Aug","Spt6","Spt16")){
    z <- m.occup[,grep(chip,names(m.occup))]
    z1<- m.spike.occup[,grep(chip,names(m.spike.occup))]
    for(tmp in c(30,37)){
      exp = (1:ncol(z))[grepl(tmp,names(z),fixed = T)]
      myname= paste0( tmp,"_",chip  )
      cat(myname,"\n")
      cat(names(z)[exp],"\n\n")
      
      zz <- reshape::melt(z[,exp])
      zz$condition <- 'LibSizeNormalized'
      zz1 <- reshape::melt(z1[,exp])
      zz1$condition <- 'SpikeinNormalized'
      zz <- rbind(zz,zz1)
      
      q <- ggplot(zz, aes(x=variable,y=value,fill=variable)) + geom_violin(col='black') +
        geom_boxplot(width=0.1,fill='white', color="black",outlier.shape = NA,outlier.size = 0)+
        geom_hline(yintercept = 0,col="black",size=0.7,linetype="dashed")+
        facet_wrap(~condition, nrow = 1,scales = "free_x")+
        ylim(-2,2)+ xlab("") + ylab("relative log occupancy")+
        coord_flip()+
        theme_bw()+theme(legend.position = "none",
                         axis.title = element_text(size=14,colour = "black"),
                         axis.text = element_text(size=14,color="black"))
      pl.list[[myname]] <- print(q)
      rm(zz,zz1,q)
    }
    rm(z,z1)
  }
  
  pdf("SpikeinNorm_effectAssessment_RLEBoxes_DESeq2NormReads.pdf",width = 8,height = 6)
  for(i in 1:length(pl.list)){
    print(pl.list[[i]])
  }
  dev.off()
  
  
  ##### (3) Sanity check: PCA analyses
  rm(list=ls())
  load(file="A:/work/WinstonLab/Olga/ChIPSeq/Diff_Factor_Occupancy_DESeq2SGD.RData")
  m.occup <- data.frame(tracking_id=w[[1]]$tracking_id)
  m.spike.occup <- data.frame(tracking_id=l[[1]]$tracking_id)
  for(i in 1:length(l)){
    cat(names(ln)[i],"\n")
    cat(" ",paste0(names(ln[[i]][,5:8]),collapse = ","),"\n")
    m.spike.occup <- cbind(m.spike.occup,l[[i]][,5:8])
    m.occup <- cbind(m.occup,w[[i]][,5:8])
    cat("\n\n")
  }
  rownames(m.occup) <- rownames(m.spike.occup) <- m.occup$tracking_id
  m.occup$tracking_id <- m.spike.occup$tracking_id <- NULL 
  m.occup <- m.occup[,unique(names(m.occup))]
  m.spike.occup <- m.spike.occup[,unique(names(m.spike.occup))]
  
  pl <- data.frame(id=names(m.occup))
  pl$profile <- "Pol2-July"
  pl[grep("_Spt16_",pl$id),]$profile <- "Spt16"
  pl[grep("_Spt6_",pl$id),]$profile <- "Spt6"
  pl[grep("_Pol2-Aug_",pl$id),]$profile <- "Pol2-Aug"
  pl$repl <- "r1_30"
  pl[grep("r1_37",pl$id),]$repl <- "r1_37"
  pl[grep("r2_37",pl$id),]$repl <- "r2_37"
  pl[grep("r2_30",pl$id),]$repl <- "r2_30"
  
  
  
  pdf(file = "A:/work/WinstonLab/Olga/ChIPSeq/SampleNormCountPCA_Scer.pdf",width = 12,height = 8)
  z <- m.occup
  pca_data=prcomp(t(z))
  pca_data_perc=round(100*pca_data$sdev^2/sum(pca_data$sdev^2),1)
  df_pca_data=data.frame(PC1 = pca_data$x[,1], PC2 = pca_data$x[,2], 
                         sample = names(z), condition=pl$profile,repl=pl$repl,label=pl$id)
  ggplot(df_pca_data, aes(PC1,PC2, color = condition,label=label))+geom_text_repel(box.padding = 0.1)+
    geom_point(size=8,aes(shape=repl))+
    labs(x=paste0("PC1 (",pca_data_perc[1],")"), y=paste0("PC2 (",pca_data_perc[2],")"))+
    scale_color_brewer(palette = "Dark2")+
    ggtitle("before spike-in normalization")+
    theme(legend.position = "bottom",
          plot.title = element_text(size=16,colour = "black",hjust = 0.5))
  
  z <- m.spike.occup
  pca_data=prcomp(t(z))
  pca_data_perc=round(100*pca_data$sdev^2/sum(pca_data$sdev^2),1)
  df_pca_data=data.frame(PC1 = pca_data$x[,1], PC2 = pca_data$x[,2], 
                         sample = names(z), condition=pl$profile,repl=pl$repl,label=pl$id)
  ggplot(df_pca_data, aes(PC1,PC2, color = condition,label=label))+geom_text_repel(box.padding = 0.1)+
    geom_point(size=8,aes(shape=repl))+
    labs(x=paste0("PC1 (",pca_data_perc[1],")"), y=paste0("PC2 (",pca_data_perc[2],")"))+
    scale_color_brewer(palette = "Dark2")+
    ggtitle("after spike-in normalization")+
    theme(legend.position = "bottom",
          plot.title = element_text(size=16,colour = "black",hjust = 0.5))
    theme(legend.position = "bottom")
  dev.off()
  
  
  ### extract log2FC values: use non-spikein values
  rm(list=ls())
  load(file="A:/work/WinstonLab/Olga/ChIPSeq/Diff_Factor_Occupancy_DESeq2SGD.RData")
  m.spike <- data.frame(tracking_id=ln[[1]]$tracking_id)
  m <- data.frame(tracking_id=ln[[1]]$tracking_id)
  for(i in 1:length(ln)){
    m.spike <- cbind(m.spike,ln[[i]][,3])
    m <- cbind(m,wn[[i]][,3])
  }
  names(m)[2:25] <- names(m.spike)[2:25] <- names(l)
  
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
    names(qwe) <- gsub("_r[1234567]","",names(qwe))
    m.spike.occup <- cbind(m.spike.occup,qwe)
    rm(qwe)
    
    qwe <- wn[[i]][,5:8]
    cat(" ",paste0(names(wn[[i]][,5:8]),collapse = ","),"\n")
    qwe$x <- apply(qwe[,1:2],1,function(x) sqrt(prod(x)))
    qwe$y <- apply(qwe[,3:4],1,function(x) sqrt(prod(x)))
    names(qwe)[5:6] <- names(qwe)[c(1,3)]
    qwe <- qwe[,5:6]
    names(qwe) <- gsub("_r[1234567]","",names(qwe))
    m.occup <- cbind(m.occup,qwe)
    rm(qwe)
    cat("\n\n")
  }
  save(m.occup,m.spike.occup,m,m.spike,file="A:/work/WinstonLab/Olga/ChIPSeq/DESeq2_ChIPSeq_L2FC_DFs.RData")
  
  
  
  pl <- data.frame(profile=names(m)[2:25],
                   raw=apply(m[,2:25],2,function(x) median(x,na.rm=T) ),
                   spikein=apply(m.spike[,2:25],2,function(x) median(x,na.rm=T) )
                   )
  pl <- melt(pl,measure.vars = c("raw","spikein"))
  pl$variable <- factor(pl$variable,levels=rev(c("raw","spikein")))
  
  pdf(file="A:/work/WinstonLab/Olga/ChIPSeq/GlobalEffectAfterSpikein_MedianL2FC.pdf",height = 6,width = 5)
  ggplot(pl,aes(x=profile,y=value,fill=variable))+geom_bar(stat="identity")+coord_flip()+
    theme_bw()+ylim(-1,1)
  dev.off()
  
}







