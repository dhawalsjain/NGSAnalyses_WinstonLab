#################### Natalia :thesis figures
rm(list=ls())
library(ggplot2)
library(data.table)
library(GenomicRanges)
library(Hmisc)
library(DESeq2)
library(ggpubr)
library(ggrepel)
library(gridExtra)
library(grid)
library(tidyr)
library(cowplot)


spiked.chr <- c("chrXVII","chrXVIII","chrXIX")
chr.flt <- c("chrI","chrII","chrIII","chrIV","chrIX","chrV","chrVI","chrVII","chrVIII","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
OUTDIR="A:/work/WinstonLab/Natalia_Reim/ThesisFigures_March2019"
OUTPREFIX="NR_March19"
cols = c("non-depleted"="#4477AA", 
         "Spn1-depleted"="#BB5566", 
         "WT_DMSO"="#DDAA33", 
         "WT_IAA"="#dd8033",
         "SPN1-AID_DMSO"="#4477AA", 
         "SPN1-AID_IAA"="#BB5566",
         "SPN1-AID_DMSO-1"="#4477AA", 
         "SPN1-AID_IAA-1"="#BB5566",
         "SPN1-AID_DMSO-2"="#4477AA", 
         "SPN1-AID_IAA-2"="#BB5566",
         "Spn1-AID_DMSO"="#4477AA", 
         "Spn1-AID_IAA"="#BB5566",
         "Spn1-AID_DMSO-1"="#4477AA", 
         "Spn1-AID_IAA-1"="#BB5566",
         "Spn1-AID_DMSO-2"="#4477AA", 
         "Spn1-AID_IAA-2"="#BB5566",
         "noAID_IAA-1"="#4477AA","noAID_IAA-2"="#4477AA",
         
         "Depleted (1)" = "#BB5566","Depleted (2)" = "red4",
         "Depleted (3)" = "#BB5566","Depleted (4)" = "red4",
         "Non-depleted (1)" = "#4477AA","Non-depleted (2)" ="slateblue4",
         "Non-depleted (3)" = "#4477AA","Non-depleted (4)" = "slateblue4")

gg <- theme_minimal()+
  theme(axis.text.x = element_text(size=16,colour = 'black'),
        axis.text.y = element_text(size=16,colour = 'black'),
        legend.position = "right",
        axis.title = element_text(size=18,color="black"),
        legend.text = element_text(size=15,color="black"),
        legend.title = element_text(size=16,color="black"),
        axis.line = element_line(colour = "black",size = 0.7),
        axis.ticks = element_line(colour = "black",size = 0.5)) 


### RT PCR bar plots
RTPCR_Barplot <- function(Inp.file="A:/work/WinstonLab/Natalia_Reim/ThesisFigures_March2019/RTPCR_verification_data_NR.txt")
{
  d <- read.delim(Inp.file,header =T)
  d$mean <- apply(d[,3:5],1,function(x) mean(x,na.rm=T))
  d$sd <- apply(d[,3:5],1,function(x) sd(x,na.rm=T))
  d$se <- apply(d[,3:5],1,function(x){sd(x,na.rm=T)/sqrt(length(x))})
  d$ci <- d$se*1.96
  
  d$gene <- factor(d$gene,levels=c("KAP120","GCV3","SRB4","SER3","UBI4","HSP82"))
  
  z <- d[,c("gene","condition","rep1","rep2","rep3")]
  z <- reshape::melt(z,measure.vars=c("rep1","rep2","rep3"))
  
  
  
  pl<-ggplot(z, aes(gene,value,fill = condition))+
    geom_bar(position = 'dodge', stat = 'summary', fun.y = 'mean') +
    geom_errorbar(stat = 'summary', width = 0.25,position=position_dodge(.9),
                  fun.ymin = function(x) mean(x) - sd(x),
                  fun.ymax = function(x) mean(x) + sd(x)) +
    geom_point(aes(x = gene), position = position_dodge(width = 0.9)) +gg+
    facet_wrap(~gene,nrow=1,scales = "free")+
    theme(axis.text.x = element_text(face = "italic"))+
    xlab("")+ylab("relative mRNA \nlevels")+theme(strip.text = element_blank())+
    scale_fill_manual(values = cols)+theme(legend.title = element_blank())+
    scale_y_continuous(expand = c(0,0))+
    theme(axis.ticks= element_line(size=0.75,color = "black"),
          axis.line = element_line(size = 0.5,colour = "black"))
    
  return(pl)
}

Splicing_Barplot <- function(Inp.file="A:/work/WinstonLab/Natalia_Reim/ThesisFigures_March2019/NR_March19_splicing.txt")
{
  d <- read.delim(Inp.file,header =T)
  d$mean <- apply(d[,4:5],1,mean)
  d$sd <- apply(d[,4:5],1,sd)
  d$se <- apply(d[,4:5],1,function(x){sd(x)/sqrt(length(x))})
  d$ci <- d$se*1.96
  
  d$gene <- factor(d$gene,levels=c("RPL14B","RPL34A","RPS27A","RPL43B","RPS21B","SFT1"))
  d$status <- factor(d$status,levels=c("unchanged","changed"))
  d$genotype <- gsub(" ","_",d$genotype)
  d$genotype <- factor(d$genotype,levels=c("WT_DMSO","SPN1-AID_IAA"))
  
  z <- d[,c("gene","genotype","status","rep1","rep2")]
  z <- reshape::melt(z,measure.vars=c("rep1","rep2"))
  z$gene <- factor(z$gene, levels=c("RPL43B","RPS21B","SFT1","RPL14B","RPL34A","RPS27A"))
  
  
  pl <- ggplot(z, aes(gene,value,fill = genotype))+
    geom_bar(position = 'dodge', stat = 'summary', fun.y = 'mean') +
    geom_errorbar(stat = 'summary', width = 0.25,position=position_dodge(.9),
                  fun.ymin = function(x) mean(x) - sd(x),
                  fun.ymax = function(x) mean(x) + sd(x)) +
    geom_point(aes(x = gene), position = position_dodge(width = 0.9)) +gg+
    facet_wrap(~gene,nrow=1,scales = "free")+
    xlab("")+ylab("unspliced/spliced \n transcripts (A.U.)")+
    scale_fill_manual(values = cols)+
    scale_y_continuous(expand = c(0,0))+
    theme(axis.text.x = element_text(face = "italic"),
          axis.ticks = element_line(size=1,color = "gray"),
          axis.line = element_line(size = 1,colour = "gray"),
          strip.text = element_blank(),
          legend.title = element_blank())
  
  return(pl)
}

#Splicing_Barplot()
#
batch_effects <- function(cntdata.path="A:/work/WinstonLab/Natalia_Reim/Analysis_Feb18/Sense_Antisese_ReadCounts_revised.RData", 
                          gene.anno.path="A:/work/yeast_Annotations/ScSp_SGD.RevisedGenesFull.RData",
                          threshold=0.1,
                          samples=c("SPN1-AID_DMSO-1", "SPN1-AID_DMSO-2", "SPN1-AID_IAA-1", "SPN1-AID_IAA-2"),
                          OUTDIR="A:/work/WinstonLab/Natalia_Reim/ThesisFigures_March2019",
                          OUTPREFIX="NR_March19",is.return=F){
  
  spiked.chr <- c("chrXVII","chrXVIII","chrXIX")
  chr.flt <- c("chrI","chrII","chrIII","chrIV","chrIX","chrV","chrVI","chrVII","chrVIII","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
  load(cntdata.path)
  
  mydeseq2func <- function(s,dd,chr.flt,samples,spikenorm=F,model.batch=F,toplotheat=F,
                           toplotpca=F,titl="myplot"){
    df1.e <- subset(s, s$chr%in%chr.flt)
    df1.e <-as.matrix( round (df1.e[,samples] ))
    type <- colnames(df1.e)
    condition <- sub("-[012]","",type)
    
    if(model.batch==T){
      batch <- factor(rep(1:2,6))
      cat("modeling batch effects\n")
      colData <- data.frame(batch,condition,type)
      row.names(colData) <- colData$type
      dds <- DESeqDataSetFromMatrix(countData = df1.e,colData = colData,design = ~ batch+condition)
    }else{
      colData <- data.frame(condition,type)
      row.names(colData) <- colData$type
      dds <- DESeqDataSetFromMatrix(countData = df1.e,colData = colData,design = ~ condition)
    }
    
    if(spikenorm==T){
      sizeFactors(dds) <- dd[samples]
    }
    dds <- DESeq(dds)
    rld <- rlog(dds)
    cat(resultsNames(dds),"\n")
    if(toplotheat==T){
      hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
      distsRL <- dist(t(assay(rld)))
      mat <- as.matrix(distsRL)
      rownames(mat) <- colnames(mat) <- with(colData(dds),paste(type))
      hc <- hclust(distsRL)
      heatmap(mat, Rowv=as.dendrogram(hc),symm=TRUE, trace="none",col = rev(hmcol), margin=c(10, 10),main = titl)
    }
    if(toplotpca == T){
      plotPCA(rld, intgroup=c("type"),main=titl)
    }
    
    #return(dds)
  }
  mydeseq2func(s,dd,c(chr.flt),samples = names(s)[1:12],
               spikenorm = T,model.batch = T,toplotheat = T,toplotpca = F,titl = "library size normalization")
  
  mydeseq2func(s,dd,c(chr.flt,spiked.chr),names(s)[1:12],spikenorm = T)
  mydeseq2func(s,dd,c(chr.flt),names(s)[1:12],spikenorm = F)
  mydeseq2func(s,dd,spiked.chr,names(s)[1:12],spikenorm = F)
  mydeseq2func(s,dd,spiked.chr,names(s)[1:12],spikenorm = T)
  mydeseq2func(s,dd,c(chr.flt),names(s)[1:12],spikenorm = T)
  
  
  #### sva
  mysvafun <- function(s){
    n <- s[,1:12]
    library(sva)
    library(limma)
    n <- as.matrix(n[rowSums(n)>0,])
    info <- data.frame(sampleID = colnames(n), sample=sub("-[1234]","",colnames(n)), 
                       Batch=rep(1:2,6))
    
    ## estimating numbers of significant variables  
    mod <- model.matrix(~as.factor(Batch),data = info)
    n.sv = num.sv(n,mod,method="leek")  ## too many!!
    ## Cosindering pipeline (includes experiment and count pipeline) as the major variable
    mod1 = model.matrix(~as.factor(info$Batch))
    mod0 = cbind(mod1[,1])
    svseq = svaseq(n,mod1,mod0,n.sv=1)$sv
    
    #pdf("Estimating batch effect using SVA.pdf",8,8)
    par(mar=c(10,4,2,2))
    plot(svseq,pch=19,col="red",xaxt="n",xlab="")
    axis(1,at=1:ncol(n),labels = colnames(n),las=2,tick = TRUE)
    abline(h = 0)
    #dev.off()
    
    batch = info$Batch
    modcombat = model.matrix(~1, data=info)
    m = as.data.frame(ComBat(dat=n, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=TRUE))
    m <- cbind(m,s[,13])
    names(m)[13] <- "chr"
    m
  }
  
  s <- s[rowSums(s[,1:12])>0,]
  
  sc <- mysvafun(s)
  sc[,1:12] <- round(sc[,1:12])
  
  
  
}

## Figure1 :
DE_expression_genes_tests <- function(cntdata.path="A:/work/WinstonLab/Natalia_Reim/Analysis_Feb18/Sense_Antisese_ReadCounts_revised.RData", 
                                gene.anno.path="A:/work/yeast_Annotations/ScSp_SGD.RevisedGenesFull.RData",
                                threshold=0.1,
                                samples=c("SPN1-AID_DMSO-1", "SPN1-AID_DMSO-2", "SPN1-AID_IAA-1", "SPN1-AID_IAA-2"),
                                OUTDIR="A:/work/WinstonLab/Natalia_Reim/ThesisFigures_March2019",
                                OUTPREFIX="NR_March19",is.return=F)
{
  library(DESeq2)
  library(GenomicRanges)
  
  load(cntdata.path) ## image contains sense, antisene count table, size factors for both using spombe
  load(gene.anno.path)
  s <- subset(s,rownames(s)%in%genes[genes$verification=="Verified",]$tracking_id)
  as <- subset(as,rownames(as)%in%genes[genes$verification=="Verified",]$tracking_id)
  s <- s[complete.cases(s),]
  
  mydeseq2func <- function(s,dd,chr.flt,samples,threshold=NULL,spikenorm=T,removebatch=F){
    df1.e <- subset(s, s$chr%in%chr.flt)
    df1.e <-as.matrix( round (df1.e[,samples] ))
    type <- colnames(df1.e)
    condition <- sub("-[0123456789]","",type)
    
    if(removebatch==T){
      cat("modeling batch\t")
      batch <- factor(gsub("\\S*-","",type))
      cat(batch,"\n")
      colData <- data.frame(batch,condition,type)
      row.names(colData) <- colData$type
      dds <- DESeqDataSetFromMatrix(countData = df1.e,colData = colData,design = ~ batch+condition)
    }else{
      colData <- data.frame(condition,type)
      row.names(colData) <- colData$type
      dds <- DESeqDataSetFromMatrix(countData = df1.e,colData = colData,design = ~ condition)
    }
    if(spikenorm==T){
      sizeFactors(dds) <- dd[samples]
    }
    dds <- DESeq(dds)
    cat(resultsNames(dds),"\n")
    rld <- rlog(dds)
    hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
    distsRL <- dist(t(assay(rld)))
    mat <- as.matrix(distsRL)
    rownames(mat) <- colnames(mat) <- with(colData(dds),paste(type))
    hc <- hclust(distsRL)
    heatmap(mat, Rowv=as.dendrogram(hc),symm=TRUE, trace="none",col = rev(hmcol), margin=c(10, 10))
    
    
    res <- as.data.frame(results(dds))
    cnt <- as.data.frame(DESeq2::counts(dds,normalized=T))
    res <- cbind(res,cnt)
    res$tracking_id <- rownames(res)
    res$comparison <- paste0("log2(",gsub("_vs_","/",gsub("condition_","",colnames(dds@rowRanges@elementMetadata)[12])),")")
    res <- merge(res,genes,by="tracking_id",all.x=T)
    res$prank <- trunc(rank(res$baseMean))/length(res$baseMean)
    if(!is.null(threshold)){
      res <- res[res$padj<threshold & !is.na(res$padj),]
      cat("upregulated genes: ",sum(res$log2FoldChange>0),"\n")
      cat("downregulated genes: ",sum(res$log2FoldChange<0),"\n")
    }
    return(res)
  }
  
  threshold =0.1
  res <- mydeseq2func(s,dd,chr.flt,samples = c("SPN1-AID_DMSO-1", "SPN1-AID_DMSO-2", "SPN1-AID_IAA-1", "SPN1-AID_IAA-2"))
  if(T){
    res<- mydeseq2func(s,dd,chr.flt,samples = c("noAID_IAA-1", "noAID_IAA-2", "SPN1-AID_IAA-1", "SPN1-AID_IAA-2"))
    save(res,file="A:/work/WinstonLab/Natalia_Reim/ThesisFigures_March2019/NR_March19_SenseDE_Spn1vsnoAID.Data.RData")
  }
  
  res1 <- mydeseq2func(s,dd,chr.flt,samples = c("WT_DMSO-1", "WT_DMSO-2","SPN1-AID_DMSO-1", "SPN1-AID_DMSO-2"),threshold=threshold)
  res1.1 <- mydeseq2func(s,dd,chr.flt,samples = c("WT_DMSO-1", "WT_DMSO-2","SPN1-AID_DMSO-1", "SPN1-AID_DMSO-2"),threshold=threshold,spikenorm = F)
  
  res2 <- mydeseq2func(s,dd,chr.flt,samples = c("WT_DMSO-1", "WT_DMSO-2","WT_IAA-1", "WT_IAA-2"),threshold=threshold)
  res2.1 <- mydeseq2func(s,dd,chr.flt,samples = c("WT_DMSO-1", "WT_DMSO-2","WT_IAA-1", "WT_IAA-2"))
  res3 <- mydeseq2func(s,dd,chr.flt,samples = c("WT_DMSO-1", "WT_DMSO-2","noAID_DMSO-1", "noAID_DMSO-2"),threshold=threshold)
  res4 <- mydeseq2func(s,dd,chr.flt,samples = c("WT_DMSO-1", "WT_DMSO-2","noAID_IAA-1", "noAID_IAA-2"),threshold=threshold)
  #res <- mydeseq2func(s,dd,chr.flt,samples = names(s)[1:4],spikenorm = F,removebatch = T,threshold = 0.1)
  
  res2.1$`WT_DMSO` <- apply(res2.1[,8:9],1,function(x){sqrt(prod(x))})
  write.table(res2.1,file=paste0(OUTDIR,"/RNASeq_WTComparison.txt"),sep="\t",quote = F,row.names = F)
  
  if(F){
    res5 <- mydeseq2func(s,dd,chr.flt,samples = c("WT_DMSO-1", "WT_DMSO-2","SPN1-AID_IAA-1", "SPN1-AID_IAA-2"),threshold=threshold)
    res6 <- mydeseq2func(s,dd,chr.flt,samples = c("WT_IAA-1", "WT_IAA-2","SPN1-AID_IAA-1", "SPN1-AID_IAA-2"),threshold=threshold)
    res7 <- mydeseq2func(s,dd,chr.flt,samples = c("noAID_DMSO-1", "noAID_DMSO-2","SPN1-AID_IAA-1", "SPN1-AID_IAA-2"),threshold=threshold)
    res8 <- mydeseq2func(s,dd,chr.flt,samples = c("noAID_IAA-1", "noAID_IAA-2","SPN1-AID_IAA-1", "SPN1-AID_IAA-2"),threshold=threshold)
    res81 <- mydeseq2func(s,dd,chr.flt,samples = c("noAID_IAA-1", "noAID_IAA-2","SPN1-AID_IAA-1", "SPN1-AID_IAA-2"),threshold=threshold,spikenorm = F)
    res9 <- mydeseq2func(s,dd,chr.flt,samples = c("SPN1-AID_DMSO-1", "SPN1-AID_DMSO-2", "SPN1-AID_IAA-1", "SPN1-AID_IAA-2"),threshold=threshold)
    res91 <- mydeseq2func(s,dd,chr.flt,samples = c("SPN1-AID_DMSO-1", "SPN1-AID_DMSO-2", "SPN1-AID_IAA-1", "SPN1-AID_IAA-2"),threshold=threshold,spikenorm = F)
    #res91 <- mydeseq2func(s,dd,chr.flt,samples = c("SPN1-AID_DMSO-1", "SPN1-AID_DMSO-2", "SPN1-AID_IAA-1", "SPN1-AID_IAA-2"),threshold=0.05,spikenorm = F)
    
    resx <- mydeseq2func(s,dd,chr.flt,samples = c("SPN1-AID_DMSO-1", "SPN1-AID_DMSO-2","noAID_DMSO-1", "noAID_DMSO-2"),threshold=threshold,spikenorm = T)
    resx <- mydeseq2func(s,dd,chr.flt,samples = c("SPN1-AID_DMSO-1", "SPN1-AID_DMSO-2","noAID_IAA-1", "noAID_IAA-2"),threshold=threshold,spikenorm = T)
    
    pdf(file=paste0(OUTDIR,"/Spn1Depl to noAIDIAA comparison wWo Spikein.pdf"),width = 12,height = 10)
    par(mfrow=c(2,2))
    Venn2(res81[res81$log2FoldChange<0,]$tracking_id,res91[res91$log2FoldChange<0,]$tracking_id,
          names=c("noAID-IAA","SPN1-AID_DMSO"))
    mtext("No spikein; Down in Spn1-depl",side = 3,cex = 2)
    Venn2(res81[res81$log2FoldChange>0,]$tracking_id,res91[res91$log2FoldChange>0,]$tracking_id,
          names=c("noAID-IAA","SPN1-AID_DMSO"))
    mtext("No spikein; Up in Spn1-depl",side = 3,cex = 2)
    
    Venn2(res8[res8$log2FoldChange<0,]$tracking_id,res9[res9$log2FoldChange<0,]$tracking_id,
          names=c("noAID-IAA","SPN1-AID_DMSO"))
    mtext("Spikein; Down in Spn1-depl",side = 3,cex = 2)
    Venn2(res8[res8$log2FoldChange>0,]$tracking_id,res9[res9$log2FoldChange>0,]$tracking_id,
          names=c("noAID-IAA","SPN1-AID_DMSO"))
    mtext("Spikein; Up in Spn1-depl",side = 3,cex = 2)
    dev.off()

    s1 <- s
    dd1 <- dd
    dd1[seq(5,11,2)+1] <- dd1[seq(5,11,2)]
    names(s1)[5:12] <- names(dd1)[5:12] <- c("SPN1-AID_DMSO-3","SPN1-AID_DMSO-4","SPN1-AID_DMSO-5","SPN1-AID_DMSO-6",
                         "SPN1-AID_DMSO-7","SPN1-AID_DMSO-8","SPN1-AID_DMSO-9","SPN1-AID_DMSO-0")
    
    
    res92 <- mydeseq2func(s1,dd1,chr.flt,samples = c("SPN1-AID_DMSO-1", "SPN1-AID_DMSO-2",
                                                     "SPN1-AID_DMSO-3","SPN1-AID_DMSO-4",
                                                     "SPN1-AID_DMSO-5","SPN1-AID_DMSO-6",
                                                     "SPN1-AID_DMSO-7","SPN1-AID_DMSO-8",
                                                     "SPN1-AID_DMSO-9","SPN1-AID_DMSO-0", "SPN1-AID_IAA-1", "SPN1-AID_IAA-2"),threshold=0.1)
    res93 <- mydeseq2func(s1,dd1,chr.flt,samples = c("SPN1-AID_DMSO-1", "SPN1-AID_DMSO-2",
                                                     "SPN1-AID_DMSO-3","SPN1-AID_DMSO-4",
                                                     "SPN1-AID_DMSO-5","SPN1-AID_DMSO-6",
                                                     "SPN1-AID_DMSO-7","SPN1-AID_DMSO-8",
                                                     "SPN1-AID_DMSO-9","SPN1-AID_DMSO-0", "SPN1-AID_IAA-1", "SPN1-AID_IAA-2"),threshold=0.05)
    s1 <- s
    dd1 <- dd
    names(dd1)[7:8] <- names(s1)[7:8] <- c("SPN1-AID_DMSO-3","SPN1-AID_DMSO-4")
    res92 <- mydeseq2func(s1,dd1,chr.flt,samples = c("SPN1-AID_DMSO-1", "SPN1-AID_DMSO-2","SPN1-AID_DMSO-3","SPN1-AID_DMSO-4", "SPN1-AID_IAA-1", "SPN1-AID_IAA-2"),threshold=0.05,spikenorm = T)
    res93 <- mydeseq2func(s1,dd1,chr.flt,samples = c("SPN1-AID_DMSO-1", "SPN1-AID_DMSO-2","SPN1-AID_DMSO-3","SPN1-AID_DMSO-4", "SPN1-AID_IAA-1", "SPN1-AID_IAA-2"),threshold=0.05,spikenorm = F)
    
    pdf(file=paste0(OUTDIR,"/PieChartComparison of DE_allControls genes.pdf"),width = 10,height = 5)
    par(mfrow=c(1,2))
    v <- c(5134-nrow(res92),sum(res92$log2FoldChange<0),sum(res92$log2FoldChange>0))
    names(v) <- c("Unchanged","Down","Up")
    pie(v,labels = paste0(names(v),"\n(",v,")"))
    mtext("10% FDR cutoff",3,cex=2)
    v <- c(5134-nrow(res93),sum(res93$log2FoldChange<0),sum(res93$log2FoldChange>0))
    names(v) <- c("Unchanged","Down","Up")
    pie(v,labels = paste0(names(v),"\n(",v,")"))
    mtext("5% FDR cutoff",3,cex=2)
    dev.off()
    
    write.xlsx(res92,file = paste0(OUTDIR,"/MulticontrolDE_analyses_allC.xlsx"),sheetName = "10%FDR",row.names = F,append = F)
    write.xlsx(res93,file = paste0(OUTDIR,"/MulticontrolDE_analyses_allC.xlsx"),sheetName = "5%FDR",row.names = F,append = T)
    
    pdf(file=paste0(OUTDIR,"/NR_March19_SPN1Depletion_compa2allcontrols_p10.pdf"),height =6,width = 7)
    if(T){
      library(Vennerable)
      library(UpSetR)
      dg <- list(WT_DMSO=res5[res5$log2FoldChange>0,]$tracking_id,
                 WT_IAA=res6[res6$log2FoldChange>0,]$tracking_id,
                 noAID_DMSO=res7[res7$log2FoldChange<0,]$tracking_id,
                 noAID_IAA=res8[res8$log2FoldChange<0,]$tracking_id,
                 spn1_dmso=res9[res9$log2FoldChange<0,]$tracking_id)
      dg <- Venn(dg)
      Tstem <- compute.Venn(dg)
      gp <- VennThemes(Tstem, colourAlgorithm = "sequential")
      gps <- gp[["Set"]]
      nSets <- length(gps)
      for (ix in 1:nSets) {
        gps[[ix]]$lwd <- nSets + 1 - ix
      }
      gp[["Set"]] <- gps
      plot(dg, type = "ChowRuskey", show = list(SetLabels = TRUE),gp = gp)
      
      dg <- list(WT_DMSO=res5[res5$log2FoldChange<0,]$tracking_id,
                 WT_IAA=res6[res6$log2FoldChange<0,]$tracking_id,
                 noAID_DMSO=res7[res7$log2FoldChange>0,]$tracking_id,
                 noAID_IAA=res8[res8$log2FoldChange>0,]$tracking_id,
                 spn1_dmso=res9[res9$log2FoldChange>0,]$tracking_id)
      dg <- Venn(dg)
      Tstem <- compute.Venn(dg)
      gp <- VennThemes(Tstem, colourAlgorithm = "sequential")
      gps <- gp[["Set"]]
      nSets <- length(gps)
      for (ix in 1:nSets) {
        gps[[ix]]$lwd <- nSets + 1 - ix
      }
      gp[["Set"]] <- gps
      plot(dg, type = "ChowRuskey", show = list(SetLabels = TRUE),gp = gp)
      
      dg <- list(WT_DMSO=res5[res5$log2FoldChange<0,]$tracking_id,
                 WT_IAA=res6[res6$log2FoldChange<0,]$tracking_id,
                 noAID_DMSO=res7[res7$log2FoldChange>0,]$tracking_id,
                 noAID_IAA=res8[res8$log2FoldChange>0,]$tracking_id,
                 spn1_dmso=res9[res9$log2FoldChange>0,]$tracking_id)
      m <- data.frame(tracking_id=res$tracking_id,
                      WT_DMSO=0,WT_IAA=0,noAID_DMSO=0,noAID_IAA=0,SPN1_DMSO=0)
      rownames(m) <- m$tracking_id
      m$WT_DMSO <- ifelse(m$tracking_id%in%dg$WT_DMSO,1,0)
      m$WT_IAA <- ifelse(m$tracking_id%in%dg$WT_IAA,1,0)
      m$noAID_DMSO <- ifelse(m$tracking_id%in%dg$noAID_DMSO,1,0)
      m$noAID_IAA <- ifelse(m$tracking_id%in%dg$noAID_IAA,1,0)
      m$SPN1_DMSO <- ifelse(m$tracking_id%in%dg$spn1_dmso,1,0)
      m$tracking_id <- NULL
      m <- m[rowSums(m)>0,]
      names(m) <- paste0(names(m)," (U)")
      print(upset(m, sets = names(m), mb.ratio = c(0.55, 0.45), order.by = "freq"))
      n <- m
      n$tracking_id <- rownames(n)
      n <- merge(n,genes,by="tracking_id",all.x=T)
      names(n)[2:6] <- paste0("log2 (SPN1-AID_IAA/",names(n)[2:6],")")
      
      dg <- list(WT_DMSO=res5[res5$log2FoldChange>0,]$tracking_id,
                 WT_IAA=res6[res6$log2FoldChange>0,]$tracking_id,
                 noAID_DMSO=res7[res7$log2FoldChange<0,]$tracking_id,
                 noAID_IAA=res8[res8$log2FoldChange<0,]$tracking_id,
                 spn1_dmso=res9[res9$log2FoldChange<0,]$tracking_id)
      m <- data.frame(tracking_id=res$tracking_id,
                      WT_DMSO=0,WT_IAA=0,noAID_DMSO=0,noAID_IAA=0,SPN1_DMSO=0)
      rownames(m) <- m$tracking_id
      m$WT_DMSO <- ifelse(m$tracking_id%in%dg$WT_DMSO,1,0)
      m$WT_IAA <- ifelse(m$tracking_id%in%dg$WT_IAA,1,0)
      m$noAID_DMSO <- ifelse(m$tracking_id%in%dg$noAID_DMSO,1,0)
      m$noAID_IAA <- ifelse(m$tracking_id%in%dg$noAID_IAA,1,0)
      m$SPN1_DMSO <- ifelse(m$tracking_id%in%dg$spn1_dmso,1,0)
      m$tracking_id <- NULL
      m <- m[rowSums(m)>0,]
      names(m) <- paste0(names(m)," (D)")
      print(upset(m, sets = names(m), mb.ratio = c(0.55, 0.45), order.by = "freq"))
      m$tracking_id <- rownames(m)
      m <- merge(m,genes,by="tracking_id",all.x=T)
      names(m)[2:6] <- paste0("log2 (SPN1-AID_IAA/",names(m)[2:6],")")
      
      
      ### Only noAID-IAA, SPN1-DMSO, SPN1-IAA
      dg <- list(noAID_IAA=res8[res8$log2FoldChange>0,]$tracking_id,
                 spn1_dmso=res9[res9$log2FoldChange>0,]$tracking_id)
      m <- data.frame(tracking_id=res$tracking_id,
                      noAID_IAA=0,SPN1_DMSO=0)
      rownames(m) <- m$tracking_id
      m$noAID_IAA <- ifelse(m$tracking_id%in%dg$noAID_IAA,1,0)
      m$SPN1_DMSO <- ifelse(m$tracking_id%in%dg$spn1_dmso,1,0)
      m$tracking_id <- NULL
      m <- m[rowSums(m)>0,]
      names(m) <- paste0(names(m)," (U)")
      print(upset(m, sets = names(m), mb.ratio = c(0.55, 0.45), order.by = "freq"))
      
      dg <- list(noAID_IAA=res8[res8$log2FoldChange<0,]$tracking_id,
                 spn1_dmso=res9[res9$log2FoldChange<0,]$tracking_id)
      m <- data.frame(tracking_id=res$tracking_id,
                      noAID_IAA=0,SPN1_DMSO=0)
      rownames(m) <- m$tracking_id
      m$noAID_IAA <- ifelse(m$tracking_id%in%dg$noAID_IAA,1,0)
      m$SPN1_DMSO <- ifelse(m$tracking_id%in%dg$spn1_dmso,1,0)
      m$tracking_id <- NULL
      m <- m[rowSums(m)>0,]
      names(m) <- paste0(names(m)," (D)")
      print(upset(m, sets = names(m), mb.ratio = c(0.55, 0.45), order.by = "freq"))
    }
    dev.off()
    
    
    library(xlsx)
    write.xlsx(res5,file = paste0(OUTDIR,"/NR_March19_Log2FC_Depletion2AllControls_10pFDR.xlsx"),sheetName = "WT_DMSO",row.names = F)
    write.xlsx(res6,file = paste0(OUTDIR,"/NR_March19_Log2FC_Depletion2AllControls_10pFDR.xlsx"),sheetName = "WT_IAA",row.names = F,append = T)
    write.xlsx(res7,file = paste0(OUTDIR,"/NR_March19_Log2FC_Depletion2AllControls_10pFDR.xlsx"),sheetName = "noAID_DMSO",row.names = F,append = T)
    write.xlsx(res8,file = paste0(OUTDIR,"/NR_March19_Log2FC_Depletion2AllControls_10pFDR.xlsx"),sheetName = "noAID_IAA",row.names = F,append = T)
    write.xlsx(res9,file = paste0(OUTDIR,"/NR_March19_Log2FC_Depletion2AllControls_10pFDR.xlsx"),sheetName = "SPN1-AID_DMSO",row.names = F,append = T)
    write.xlsx(m,file = paste0(OUTDIR,"/NR_March19_Log2FC_Depletion2AllControls_10pFDR.xlsx"),sheetName = "OverlapAnalysis",row.names = F,append = T)
    write.table(m,file = paste0(OUTDIR,"/temp.txt"),row.names = F,quote = F,sep = "\t")
    write.table(n,file = paste0(OUTDIR,"/temp.txt"),row.names = F,quote = F,sep = "\t")
    
  }
  
  res.cont <- rbind(res1,res2,res3,res4)
  de.genes <- res[!is.na(res$padj) & res$padj<threshold & res$log2FoldChange<0,]
  save(de.genes,file=paste0(OUTDIR,"/",OUTPREFIX,"_10pFDR_downGenes.RData"))
  
  ## print the report page
  write.xlsx(res,file=paste0(OUTDIR,"/",OUTPREFIX,"_SenseDE.xlsx"),sheetName = "senseAllGenes",row.names = F,append = F,showNA = T)
  write.xlsx(res[!is.na(res$padj) & res$padj < threshold & res$log2FoldChange>0,],file=paste0(OUTDIR,"/",OUTPREFIX,"_SenseDE.xlsx"),sheetName = "senseUP",row.names = F,append = T,showNA = T)
  write.xlsx(res[!is.na(res$padj) & res$padj < threshold & res$log2FoldChange<0,],file=paste0(OUTDIR,"/",OUTPREFIX,"_SenseDE.xlsx"),sheetName = "senseDown",row.names = F,append = T,showNA = T)
  write.xlsx(res.cont,file=paste0(OUTDIR,"/",OUTPREFIX,"_SenseDE.xlsx"),sheetName = "senseContaminant",row.names = F,append = T,showNA = T)
  
  save(res,file=paste0(OUTDIR,"/",OUTPREFIX,"_SenseDE.Data.RData"))
  cat("Wrote output files: \n")
  cat(" ",paste0(OUTDIR,"/",OUTPREFIX,"_SenseDE.xlsx"),"\n")
  cat(" ",paste0(OUTDIR,"/",OUTPREFIX,"_SenseDE.Data.RData"),"\n")
  if(is.return==T){
    return(res)
  }else{
    return(1)
  }
}

DE_expression_genes <- function(cntdata.path="A:/work/WinstonLab/Natalia_Reim/Analysis_Feb18/Sense_Antisese_ReadCounts_revised.RData", 
                                      gene.anno.path="A:/work/yeast_Annotations/ScSp_SGD.RevisedGenesFull.RData",
                                      threshold=0.1,
                                      samples=c("SPN1-AID_DMSO-1", "SPN1-AID_DMSO-2", "SPN1-AID_IAA-1", "SPN1-AID_IAA-2"),
                                      OUTDIR="A:/work/WinstonLab/Natalia_Reim/ThesisFigures_March2019",
                                      OUTPREFIX="NR_March19",is.return=F){
  library(DESeq2)
  library(GenomicRanges)
  
  load(cntdata.path) ## image contains sense, antisene count table, size factors for both using spombe
  load(gene.anno.path)
  s <- subset(s,rownames(s)%in%genes[genes$verification=="Verified",]$tracking_id)
  as <- subset(as,rownames(as)%in%genes[genes$verification=="Verified",]$tracking_id)
  s <- s[complete.cases(s),]
  
  mydeseq2func <- function(s,dd,chr.flt,samples,threshold=NULL,spikenorm=T,removebatch=F){
    df1.e <- subset(s, s$chr%in%chr.flt)
    df1.e <-as.matrix( round (df1.e[,samples] ))
    type <- colnames(df1.e)
    condition <- sub("-[0123456789]","",type)
    
    if(removebatch==T){
      cat("modeling batch\t")
      batch <- factor(gsub("\\S*-","",type))
      cat(batch,"\n")
      colData <- data.frame(batch,condition,type)
      row.names(colData) <- colData$type
      dds <- DESeqDataSetFromMatrix(countData = df1.e,colData = colData,design = ~ batch+condition)
    }else{
      colData <- data.frame(condition,type)
      row.names(colData) <- colData$type
      dds <- DESeqDataSetFromMatrix(countData = df1.e,colData = colData,design = ~ condition)
    }
    if(spikenorm==T){
      sizeFactors(dds) <- dd[samples]
    }
    dds <- DESeq(dds)
    cat(resultsNames(dds),"\n")
    rld <- rlog(dds)
    hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
    distsRL <- dist(t(assay(rld)))
    mat <- as.matrix(distsRL)
    rownames(mat) <- colnames(mat) <- with(colData(dds),paste(type))
    hc <- hclust(distsRL)
    heatmap(mat, Rowv=as.dendrogram(hc),symm=TRUE, trace="none",col = rev(hmcol), margin=c(10, 10))
    
    
    res <- as.data.frame(results(dds))
    cnt <- as.data.frame(DESeq2::counts(dds,normalized=T))
    res <- cbind(res,cnt)
    res$tracking_id <- rownames(res)
    res$comparison <- paste0("log2(",gsub("_vs_","/",gsub("condition_","",colnames(dds@rowRanges@elementMetadata)[12])),")")
    res <- merge(res,genes,by="tracking_id",all.x=T)
    res$prank <- trunc(rank(res$baseMean))/length(res$baseMean)
    if(!is.null(threshold)){
      res <- res[res$padj<threshold & !is.na(res$padj),]
      cat("upregulated genes: ",sum(res$log2FoldChange>0),"\n")
      cat("downregulated genes: ",sum(res$log2FoldChange<0),"\n")
    }
    return(res)
  }
  
  
  res <- mydeseq2func(s,dd,chr.flt,samples =samples,threshold=NULL,spikenorm = T)
  res$lfcSE <- res$stat <- res$pvalue <- NULL
  res$comparison <- "log2(SPN1.AID_IAA/SPN1.AID_DMSO_noAID-IAA)"
  names(res)[c(2,4)] <- c("average.expression","FDR")                       
  
  
  f = paste0(OUTDIR,"/",OUTPREFIX,"_SenseDE_2Controls.xlsx")
  write.xlsx(res,file=f,sheetName = "senseAllGenes",row.names = F,append = F,showNA = T)
  threshold=0.1
  write.xlsx(res[!is.na(res$FDR) & res$FDR < threshold & res$log2FoldChange>0,],file=f,sheetName = paste0("senseUP_",threshold,"FDR"),row.names = F,append = T,showNA = T)
  write.xlsx(res[!is.na(res$FDR) & res$FDR < threshold & res$log2FoldChange<0,],file=f,sheetName = paste0("senseDN_",threshold,"FDR"),row.names = F,append = T,showNA = T)
  write.table(res[!is.na(res$FDR) & res$FDR < threshold & res$log2FoldChange<0,],file = "temp",quote = F,sep = "\t",row.names = F)
  
  threshold=0.05
  write.xlsx(res[!is.na(res$FDR) & res$FDR < threshold & res$log2FoldChange>0,],file=f,sheetName = paste0("senseUP_",threshold,"FDR"),row.names = F,append = T,showNA = T)
  write.table(res[!is.na(res$FDR) & res$FDR < threshold & res$log2FoldChange>0,],file = "temp",quote = F,sep = "\t",row.names = F)
  
  write.xlsx(res[!is.na(res$FDR) & res$FDR < threshold & res$log2FoldChange<0,],file=f,sheetName = paste0("senseDN_",threshold,"FDR"),row.names = F,append = T,showNA = T)
  write.table(res[!is.na(res$FDR) & res$FDR < threshold & res$log2FoldChange<0,],file = "temp",quote = F,sep = "\t",row.names = F)
  
  
  write.xlsx(res.cont,file=paste0(OUTDIR,"/",OUTPREFIX,"_SenseDE.xlsx"),sheetName = "senseContaminant",row.names = F,append = T,showNA = T)
  
  save(res,file=paste0(OUTDIR,"/",OUTPREFIX,"_SenseDE_2Controls.Data.RData"))
  cat("Wrote output files: \n")
  cat(" ",paste0(OUTDIR,"/",OUTPREFIX,"_SenseDE.xlsx"),"\n")
  cat(" ",paste0(OUTDIR,"/",OUTPREFIX,"_SenseDE.Data.RData"),"\n")
  
  
  
}
#source("A:/work/scripts/Project_Winston/final/GeneCoveragePlot.R")

MAplot.PieChart.DEGenes <- function(DEExpression.filepath="A:/work/WinstonLab/Natalia_Reim/ThesisFigures_March2019/NR_March19_SenseDE_2Controls.Data.RData",
                          threshold=0.1,
                          samples=c("SPN1-AID_DMSO-1", "SPN1-AID_DMSO-2", "SPN1-AID_IAA-1", "SPN1-AID_IAA-2")
                          ){
  
  library(xlsx)
  require(ggplot2)
  plot.list <- list()
  
  load(DEExpression.filepath)
  #names(res)[c(2,4)] <- c("baseMean","padj")
  
  res1 <- res[res$padj<threshold,]
  res1 <- res1[!is.na(res1$padj),]
  res2 <- res1[res1$log2FoldChange>0,]
  res <- subset(res,!rownames(res)%in%rownames(res1))
  res1 <- res1[res1$log2FoldChange<0,]
  
  colnames(res1) <- make.unique(names(res1))
  colnames(res2) <- make.unique(names(res2))
  
  cols=c("Upregulated"="#90C987","Downregulated"="#F6C141","Unchanged"="gray80")
  ## MA-plot
  pl1 <- res[,c("tracking_id","baseMean","log2FoldChange","padj")]
  plot.list[["ma_plot"]] <- ggplot(pl1) + geom_point(aes(log2(baseMean),log2FoldChange),col="gray80",shape=21,size=3)+
    geom_point(data = res1,aes(log2(baseMean),log2FoldChange),col="#F6C141",shape=21,size=3)+
    geom_point(data = res2,aes(log2(baseMean),log2FoldChange),col="#90C987",shape=21,size=3)+
    xlab("Mean expression, log2") + ylab("log2 fold change")+
    gg+ theme(axis.line = element_line(size=0.75,color="black"),
              axis.ticks = element_line(size=1,colour = "black"))
  
  ## Piechart
  pl2 <- c(length(res2$tracking_id),length(res1$tracking_id),length(res$tracking_id))
  names(pl2) <- c("Upregulated","Downregulated","Unchanged")
  #names(pl2) <- paste(names(pl2)," (",pl2,")",sep="")
  pl2 <- as.data.frame(pl2)
  pl2$status <- rownames(pl2)
  names(pl2)[1] <- "numbers"
  pl2$Label <- paste0(pl2$status,"\n(",pl2$numbers,")")
  pl2$Label <- gsub("regulated","",pl2$Label)
  pl2$status <- factor(pl2$status,levels=c("Upregulated","Unchanged","Downregulated"))
  
  plot.list[["pie_chart"]] <- ggplot(pl2, aes(x="", y=numbers, fill=status))+geom_bar(width = 1, stat = "identity")+
    coord_polar("y", start=0)+scale_fill_manual(values = cols)+
    gg+theme_void()+ 
    geom_text_repel(aes(x=2,y=cumsum(pl2$numbers) - pl2$numbers / 2,label=Label),size=5,nudge_x = 0.1,segment.size = 0.6)+
      theme(legend.position = "none",
            axis.ticks.x = element_line(size = 2))
    
  
  return(plot.list)
}
#MAplot.PieChart.DEGenes(threshold = 0.05)


RNASeq.RLE.boxplot <- function(inp.file="A:/work/WinstonLab/Natalia_Reim/Analysis_Feb18/SpikeinNormalization_Assessment_RLE.RData"){
  if(F){
    rm(list=ls())
    load("A:/work/WinstonLab/Natalia_Reim/Analysis_Feb18/Sense_Antisese_ReadCounts_revised_01.RData")
    df1 <- s
    df1 <- df1[complete.cases(df1),]
    spiked.chr <- c("chrXVII","chrXVIII","chrXIX")
    chr.flt <- c("chrI","chrII","chrIII","chrIV","chrIX","chrV","chrVI","chrVII","chrVIII","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
    df1 <- df1[df1$chr%in%c(spiked.chr,chr.flt),]
    
    df1 <- df1[,c(13,1:12)]
    df2 <-df1
    df2[,2:13] <- apply(df2[,2:13],2,function(x){round(x*1000000/sum(x),2)})
    df1.c <- subset(df1, df1$chr%in%spiked.chr)
    df2.c <- subset(df2, df1$chr%in%spiked.chr)
    df1.e <- subset(df1, !df1$chr%in%spiked.chr)
    df2.e <- subset(df2, !df2$chr%in%spiked.chr)
    rm(df1,df2)
    
    # Normalization factor estimate from spikein using DESeq2
    library(DESeq2)
    df1.c <-as.matrix( round (df1.c[,2:13] ))
    type <- colnames(df1.c)
    condition <- sub("-[012]","",type)
    colData <- data.frame(condition,type)
    row.names(colData) <- colData$type
    dds <- DESeqDataSetFromMatrix(countData = df1.c,colData = colData,design = ~ condition)
    dds <- DESeq(dds)
    dd <- sizeFactors(dds)
    rm(dds)
    
    pl <- c()
    my.rle <- function(m){
      m <- t(apply(m, 1,function(x){round(x/median(x),2)}))
      m <- ifelse(m=="NaN",NA,m)
      m <- ifelse(m=="Inf",NA,m)
      m <- log(m)
      return(melt(m))
    }
    pl <- rbind(pl, cbind(my.rle(df1.c),category="Raw Sp tags"))
    pl <- rbind(pl, cbind(my.rle(df1.e[,2:13]),category="Raw Sc tags"))
    pl <- rbind(pl, cbind(my.rle(t(t(df1.c)/(dd))),category="SpikeNorm Sp tags"))
    pl <- rbind(pl, cbind(my.rle(t(t(df1.e[,2:13])/(dd))),category="SpikeNorm Sc tags"))
    save(pl,file="A:/work/WinstonLab/Natalia_Reim/Analysis_Feb18/SpikeinNormalization_Assessment_RLE.RData")
    
  }
  load(inp.file)
  names(pl)[2] <- "Var2"
  pl$category <- gsub("SpikeNorm","Normalized",pl$category)
  pl$category <- gsub("tags","reads",pl$category)
  pl$category <- gsub("Normalized","norm.",pl$category)
  pl$category <- gsub("Raw","raw",pl$category)
  pl$category <- factor(pl$category,levels=c("raw Sp reads","raw Sc reads","norm. Sp reads","norm. Sc reads"))
  pl$Var2 <- gsub("SPN","Spn",pl$Var2)
  
  q <- ggplot(pl, aes(x=Var2,y=value,fill=Var2)) #+ geom_vline(xintercept = c(0.365,0.4))
  q <- q + geom_boxplot(outlier.shape = NA) + ylab("Relative log2 expression") + facet_wrap(~category,nrow=2) +coord_flip()+xlab("")
  q <- q + ylim(-1,1)
  q <- q + geom_hline(yintercept = 0,col="black")
  q <- q + scale_color_manual(values = cols) +scale_fill_manual(values=cols)
  q <- q + theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = "grey80"),axis.line = element_line(size = 0.7, colour = "black"))
  q <- q + theme(axis.text.x = element_text(colour = 'black',size=14),
                 axis.text.y = element_text(colour = 'black',size=14),
                 axis.title = element_text(colour = 'black',size=15),
                 strip.text.x = element_text(colour = 'black',size=16),
                 strip.background = element_blank())
  q <- q + theme(legend.position = "none")
  q <- q + theme(panel.grid = element_blank())
  q
  
} 


RNASeq.SenseCoverage <- function(script.path="A:/work/scripts/Project_Winston/final/GeneCoveragePlot.R",
                                 gene.anno.path="A:/work/yeast_Annotations/ScSp_SGD.RevisedGenesFull.RData",
                                 genome.anno.path="A:/work/yeast_Annotations/Yeast_ref_Annotations_DJ.gff",
                                 gene.names =c("HSP82","SER3","KAP120","GCV3","UBI4","SRB4"),
                                 OUTDIR="A:/work/WinstonLab/Natalia_Reim/ThesisFigures_March2019",
                                 OUTPREFIX="NR_March19",which.cov="sense"){
  ## hardcoded sample names etc
  cols = c("non-depleted"="#4477AA", 
           "Spn1-depleted"="#BB5566", 
           "WT_DMSO"="#DDAA33", 
           "WT_IAA"="#dd8033",
           "antisense_SPN1-AID_DMSO"="#73A1CE", 
           "antisense_SPN1-AID_IAA"="#E27A8B", 
           "antisense_WT_DMSO"="#EFCD7C", 
           "sense_SPN1-AID_DMSO"="#4477AA", 
           "sense_SPN1-AID_IAA"="#BB5566", 
           "sense_WT_DMSO"="#DDAA33",
           "reverse"="gray80","forward"="gray80")
  #source(script.path)
  ## gene annotations
  load(gene.anno.path)
  genes$gene <- ifelse(is.na(genes$gene),as.character(genes$tracking_id),as.character(genes$gene))
  
  ## genome annotations
  d <- read.delim(genome.anno.path,header=F,comment.char = "#")
  d$tracking_id = gsub("gene_id ","",gsub("; transcript_id.*","",d$V9))
  d <- d[d$V3=="exon",]
  d <- d[,c(1:5,7,10)]
  d <- merge(d,genes[,c("tracking_id","start","end")],by="tracking_id",all.x=T)
  anno <- d
  anno <- anno[!is.na(anno$start),]
  names(anno) <- c("tracking_id","chr","biotype","feature","exon_start","exon_end","strand","gene_start","gene_end")
  anno$feature <- NULL
  anno$start <- anno$gene_start
  anno$end <- anno$gene_end
  anno <- merge(anno,genes[,c("tracking_id","gene")],by="tracking_id",all.x=T)
  anno <- as(anno,"GRanges")
  rm(d)
  
  ## gene list
  glist <- genes[genes$gene%in%gene.names,]$tracking_id
  glist=paste(glist,collapse = ",")
  
  ## create coverage plot list
  plotlist <- list()
  plotlist[["WT"]] <- gene.coverage.plot.df(path="A:/work/WinstonLab/Natalia_Reim/October2016/BedGraph/whole_read_replicates merged/WT_DMSO_Mergeq50WR.bedGraph",profile = "WT_DMSO",genes = genes,glist = glist,Margin = 500)
  plotlist[["SPN1_DMSO"]] <- gene.coverage.plot.df(path="A:/work/WinstonLab/Natalia_Reim/October2016/BedGraph/whole_read_replicates merged/SPN1_DMSO_Mergeq50WR.bedGraph",profile = "SPN1-AID_DMSO",genes = genes,glist = glist,Margin = 500)
  plotlist[["SPN1_IAA"]] <- gene.coverage.plot.df(path="A:/work/WinstonLab/Natalia_Reim/October2016/BedGraph/whole_read_replicates merged/SPN1_IAA_Mergeq50WR.bedGraph",profile = "SPN1-AID_IAA",genes = genes,glist = glist,Margin = 500)
  pl <- rbind(plotlist[[1]],plotlist[[2]],plotlist[[3]])
  pl$colguide <- paste0(pl$ori,"_",pl$profile)
  pl$profile <- factor(pl$profile,levels=c("WT_DMSO", "SPN1-AID_DMSO","SPN1-AID_IAA"))
  
  ## plot list
  if(which.cov=="sense"){
    p1 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,only.gene = T,
                             mygene = "KAP120",ori = "sense",cols = cols,anno = anno,Margin = 200)
    p2 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,only.gene = T,
                             mygene = "GCV3",ori = "sense",cols = cols,anno = anno,Margin = 200)
    p3 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,only.gene = T,
                             mygene = "SRB4",ori = "sense",cols = cols,anno = anno,Margin = 50)
    p4 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,stranded = T,
                             mygene = "SER3",ori = "sense",cols = cols,anno = anno,Margin = 200)
    p5 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,only.gene = T,
                             mygene = "UBI4",ori = "sense",cols = cols,anno = anno,Margin = 200)
    p6 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,only.gene = T,
                             mygene = "HSP82",ori = "sense",cols = cols,anno = anno,Margin = 200)
    pl.list <- list(p1[[1]],p2[[1]],p3[[1]],p1[[2]],p2[[2]],p3[[2]],
                    p4[[1]],p5[[1]],p6[[1]],p4[[2]],p5[[2]],p6[[2]])
  }
  
  if(which.cov=="splicing"){
    p1 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,only.gene = T,
                             mygene = "RPL43B",ori = "sense",cols = cols,anno = anno,Margin = 200)
    p2 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,only.gene = T,
                             mygene = "RPS21B",ori = "sense",cols = cols,anno = anno,Margin = 200)
    p3 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,only.gene = T,
                             mygene = "SFT1",ori = "sense",cols = cols,anno = anno,Margin = 200)
    p4 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,only.gene = T,
                             mygene = "RPL14B",ori = "sense",cols = cols,anno = anno,Margin = 200)
    p5 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,only.gene = T,
                             mygene = "RPL34A",ori = "sense",cols = cols,anno = anno,Margin = 200)
    p6 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,only.gene = T,
                             mygene = "RPS27A",ori = "sense",cols = cols,anno = anno,Margin = 200)
    pl.list <- list(p1[[1]],p2[[1]],p3[[1]],p1[[2]],p2[[2]],p3[[2]],
                    p4[[1]],p5[[1]],p6[[1]],p4[[2]],p5[[2]],p6[[2]])
    
  }
  
  if(which.cov=="antisense"){
    p1 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                       mygene = "CUS2",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 50)
    p2 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                       mygene = "RRP1",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 50)
    p3 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                       mygene = "MAF1",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 50)
    p4 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "FUN26",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 50)
    p5 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "GIM3",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 50)
    p6 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "PRP9",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 50)
    p7 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "SNU66",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 50)
    p8 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "VAM6",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 50)
    p9 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "STF1",ori = "",only.gene = F,cols = cols,anno = anno,Margin = 270)
    
    pl.list <- list(p1[[1]],p2[[1]],p3[[1]],p1[[2]],p2[[2]],p3[[2]],
                    p4[[1]],p5[[1]],p6[[1]],p4[[2]],p5[[2]],p6[[2]],
                    p7[[1]],p8[[1]],p9[[1]],p7[[2]],p8[[2]],p9[[2]])
  }
  
  if(which.cov=="antisense2"){
    p1 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "APP1",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 200)
    p2 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "ARF2",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 200)
    p3 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "BMT6",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 200)
    p4 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "PRS1",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 200)
    p5 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "PUF6",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 200)
    p6 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "RPL23A",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 200)
    p7 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "RPP1B",ori = "",only.gene = F,cols = cols,anno = anno,Margin = 800)
    p8 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "RPS1B",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 200)
    p9 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "RRI1",ori = "",only.gene = F,cols = cols,anno = anno,Margin = 200)
    p10 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "TFB4",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 200)
    p11 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "TRR1",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 200)
    p12 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "UMP1",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 200)
    p13 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "YCP4",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 200)
    p14 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "CYB5",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 200)
    p15 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "HMS1",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 200)
    p16 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "HOL1",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 200)
    p17 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "MOH1",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 100)
    p18 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "TMC1",ori = "",only.gene = F,cols = cols,anno = anno,Margin = 200)
    
    pl.list <- list(p1[[1]],p2[[1]],p3[[1]],p1[[2]],p2[[2]],p3[[2]],
                    p4[[1]],p5[[1]],p6[[1]],p4[[2]],p5[[2]],p6[[2]],
                    p7[[1]],p8[[1]],p9[[1]],p7[[2]],p8[[2]],p9[[2]],
                    p10[[1]],p11[[1]],p12[[1]],p10[[2]],p11[[2]],p12[[2]],
                    p13[[1]],p14[[1]],p15[[1]],p13[[2]],p14[[2]],p15[[2]],
                    p16[[1]],p17[[1]],p18[[1]],p16[[2]],p17[[2]],p18[[2]])
  }
  
  if(which.cov=="antisense3"){
    p1 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "CUS2",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 50)
    p2 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "RRP1",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 50)
    p3 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "MAF1",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 50)
    p4 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "ARF2",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 200)
    p5 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "BMT6",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 50)
    p6 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "PRS1",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 200)
    p7 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                              mygene = "TFB4",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 200)
    p8 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                             mygene = "STF1",ori = "",only.gene = F,cols = cols,anno = anno,Margin = 500)
    p9 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                              mygene = "TMC1",ori = "",only.gene = F,cols = cols,anno = anno,Margin = 200)
    p10 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                              mygene = "CYB5",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 200)
    p11 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                              mygene = "HMS1",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 200)
    p12 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                              mygene = "MOH1",ori = "",only.gene = T,cols = cols,anno = anno,Margin = 100)
    
    pl.list <- list(p1[[1]],p2[[1]],p3[[1]],p1[[2]],p2[[2]],p3[[2]],
                    p4[[1]],p5[[1]],p6[[1]],p4[[2]],p5[[2]],p6[[2]],
                    p7[[1]],p8[[1]],p9[[1]],p7[[2]],p8[[2]],p9[[2]],
                    p10[[1]],p11[[1]],p12[[1]],p10[[2]],p11[[2]],p12[[2]])
  }
  
  return(pl.list)
  ##
  #return(c(p1,p2,p3,p4,p5,p6)) 
}


ChIPSeq.Coverage <- function(script.path="A:/work/scripts/Project_Winston/final/GeneCoveragePlot.R",
                            gene.anno.path="A:/work/yeast_Annotations/ScSp_SGD.RevisedGenesFull.RData",
                            genome.anno.path="A:/work/yeast_Annotations/Yeast_ref_Annotations_DJ.gff",
                            gene.names =c("GCV3","UBI4"),
                            OUTDIR="A:/work/WinstonLab/Natalia_Reim/ThesisFigures_March2019",
                            OUTPREFIX="NR_March19"){
  ## hardcoded sample names etc
  cols = c("sense_Rpb1"="#4477AA",  "sense_Rpb1_Depl"="#BB5566", "sense_S2P"="#4477AA",  
           "sense_S2P_Depl"="#BB5566", 
           "sense_S5P"="#4477AA",  "sense_S5P_Depl"="#BB5566", "sense_Spt6"="#4477AA", 
           "sense_Spt6_Depl"="#BB5566","sense_H3"="#4477AA", 
           "sense_H3_Depl"="#BB5566")
  #source(script.path)
  ## gene annotations
  load(gene.anno.path)
  genes$gene <- ifelse(is.na(genes$gene),as.character(genes$tracking_id),as.character(genes$gene))
  
  ## genome annotations
  d <- read.delim(genome.anno.path,header=F,comment.char = "#")
  d$tracking_id = gsub("gene_id ","",gsub("; transcript_id.*","",d$V9))
  d <- d[d$V3=="exon",]
  d <- d[,c(1:5,7,10)]
  d <- merge(d,genes[,c("tracking_id","start","end")],by="tracking_id",all.x=T)
  anno <- d
  anno <- anno[!is.na(anno$start),]
  names(anno) <- c("tracking_id","chr","biotype","feature","exon_start","exon_end","strand","gene_start","gene_end")
  anno$feature <- NULL
  anno$start <- anno$gene_start
  anno$end <- anno$gene_end
  anno <- merge(anno,genes[,c("tracking_id","gene")],by="tracking_id",all.x=T)
  anno <- as(anno,"GRanges")
  rm(d)
  
  ## gene list
  glist <- genes[genes$gene%in%gene.names,]$tracking_id
  glist=paste(glist,collapse = ",")
  
  if(F){
    ## create coverage plot list
    plotlist <- list()
    plotlist[["Rpb1_Depl"]] <- gene.coverage.plot.ChIPSeq.df(path="A:/work/WinstonLab/Natalia_Reim/MayNov17/vis/Depleted_RNAPII.merged.bedGraph",profile = "Rpb1_Depl",genes = genes,glist = glist,Margin = 500)
    plotlist[["Rpb1"]] <- gene.coverage.plot.ChIPSeq.df(path="A:/work/WinstonLab/Natalia_Reim/MayNov17/vis/Non-depleted_RNAPII.merged.bedGraph",profile = "Rpb1",genes = genes,glist = glist,Margin = 500)
    plotlist[["S5P_Depl"]] <- gene.coverage.plot.ChIPSeq.df(path="A:/work/WinstonLab/Natalia_Reim/MayNov17/vis/Depleted_S5P.merged.bedGraph",profile = "S5P_Depl",genes = genes,glist = glist,Margin = 500)
    plotlist[["S5P"]] <- gene.coverage.plot.ChIPSeq.df(path="A:/work/WinstonLab/Natalia_Reim/MayNov17/vis/Non-depleted_S5P.merged.bedGraph",profile = "S5P",genes = genes,glist = glist,Margin = 500)
    plotlist[["S2P_Depl"]] <- gene.coverage.plot.ChIPSeq.df(path="A:/work/WinstonLab/Natalia_Reim/MayNov17/vis/Depleted_S2P.merged.bedGraph",profile = "S2P_Depl",genes = genes,glist = glist,Margin = 500)
    plotlist[["S2P"]] <- gene.coverage.plot.ChIPSeq.df(path="A:/work/WinstonLab/Natalia_Reim/MayNov17/vis/Non-depleted_S2P.merged.bedGraph",profile = "S2P",genes = genes,glist = glist,Margin = 500)
    plotlist[["Spt6_Depl"]] <- gene.coverage.plot.ChIPSeq.df(path="A:/work/WinstonLab/Natalia_Reim/MayNov17/vis/Depleted_Spt6.merged.bedGraph",profile = "Spt6_Depl",genes = genes,glist = glist,Margin = 500)
    plotlist[["Spt6"]] <- gene.coverage.plot.ChIPSeq.df(path="A:/work/WinstonLab/Natalia_Reim/MayNov17/vis/Non-depleted_Spt6.merged.bedGraph",profile = "Spt6",genes = genes,glist = glist,Margin = 500)
    save(plotlist,file = "A:/work/WinstonLab/Natalia_Reim/DataFrames/ChIPProfiles_CovPlotsPaper.RData")
  }
  if(F){
    plotlist <- list()
    plotlist[["H3_Depl"]] <- gene.coverage.plot.ChIPSeq.df(path="A:/work/WinstonLab/Natalia_Reim/MayNov17/vis/Depleted_H3.merged.bedGraph",profile = "H3_Depl",genes = genes,glist = glist,Margin = 500)
    plotlist[["H3"]] <- gene.coverage.plot.ChIPSeq.df(path="A:/work/WinstonLab/Natalia_Reim/MayNov17/vis/Non-depleted_H3.merged.bedGraph",profile = "H3",genes = genes,glist = glist,Margin = 500)
    save(plotlist,file = "A:/work/WinstonLab/Natalia_Reim/DataFrames/ChIPProfiles_H3CovPlotsPaper.RData")
  }
  #load(file = "A:/work/WinstonLab/Natalia_Reim/DataFrames/ChIPProfiles_CovPlotsPaper.RData")
  load(file = "A:/work/WinstonLab/Natalia_Reim/DataFrames/ChIPProfiles_H3CovPlotsPaper.RData")
  
  pl <- c()
  for(i in 1:length(plotlist)){
    pl <- rbind(pl, plotlist[[i]])
  }
  pl$colguide <- paste0(pl$ori,"_",pl$profile)
  pl$profile <- factor(pl$profile,levels=c("Rpb1","Rpb1_Depl","S5P","S5P_Depl","S2P","S2P_Depl","Spt6","Spt6_Depl","H3","H3_Depl"))
  
  ## plot list
  if(F){
    p1 <- gene.coverage.plot(pl=pl,samples = c("Rpb1","Rpb1_Depl"),genes=genes,plot.coverage = F,only.gene = T,
                             mygene = "GCV3",ori = "sense",cols = cols,anno = anno,Margin = 500)
    p2 <- gene.coverage.plot(pl=pl,samples = c("S5P","S5P_Depl"),genes=genes,plot.coverage = F,only.gene = T,
                             mygene = "GCV3",ori = "sense",cols = cols,anno = anno,Margin = 500)
    p3 <- gene.coverage.plot(pl=pl,samples = c("S2P","S2P_Depl"),genes=genes,plot.coverage = F,only.gene = T,
                             mygene = "GCV3",ori = "sense",cols = cols,anno = anno,Margin = 500)
    p4 <- gene.coverage.plot(pl=pl,samples = c("Spt6","Spt6_Depl"),genes=genes,plot.coverage = F,only.gene = T,
                             mygene = "GCV3",ori = "sense",cols = cols,anno = anno,Margin = 500)
    p5 <- gene.coverage.plot(pl=pl,samples = c("Rpb1","Rpb1_Depl"),genes=genes,plot.coverage = F,only.gene = T,
                             mygene = "UBI4",ori = "sense",cols = cols,anno = anno,Margin = 500)
    p6 <- gene.coverage.plot(pl=pl,samples = c("S5P","S5P_Depl"),genes=genes,plot.coverage = F,only.gene = T,
                             mygene = "UBI4",ori = "sense",cols = cols,anno = anno,Margin = 500)
    p7 <- gene.coverage.plot(pl=pl,samples = c("S2P","S2P_Depl"),genes=genes,plot.coverage = F,only.gene = T,
                             mygene = "UBI4",ori = "sense",cols = cols,anno = anno,Margin = 500)
    p8 <- gene.coverage.plot(pl=pl,samples = c("Spt6","Spt6_Depl"),genes=genes,plot.coverage = F,only.gene = T,
                             mygene = "UBI4",ori = "sense",cols = cols,anno = anno,Margin = 500)
    pl.list <- list(p1[[1]],p2[[1]],p3[[1]],p4[[1]],
                    p1[[2]],p2[[2]],p3[[2]],p4[[2]],
                    p5[[1]],p6[[1]],p7[[1]],p8[[1]],
                    p5[[2]],p6[[2]],p7[[2]],p8[[2]])
  }
  
  if(T){
    p1 <- gene.coverage.plot(pl=pl,samples = c("H3","H3_Depl"),genes=genes,plot.coverage = F,only.gene = T,
                             mygene = "SML1",ori = "sense",cols = cols,anno = anno,Margin = 500)
    p2 <- gene.coverage.plot(pl=pl,samples = c("H3","H3_Depl"),genes=genes,plot.coverage = F,only.gene = T,
                             mygene = "ILV5",ori = "sense",cols = cols,anno = anno,Margin = 500)
    p3 <- gene.coverage.plot(pl=pl,samples = c("H3","H3_Depl"),genes=genes,plot.coverage = F,only.gene = T,
                             mygene = "TSA1",ori = "sense",cols = cols,anno = anno,Margin = 500)
    p4 <- gene.coverage.plot(pl=pl,samples = c("H3","H3_Depl"),genes=genes,plot.coverage = F,only.gene = T,
                             mygene = "RPL10",ori = "sense",cols = cols,anno = anno,Margin = 500)
    p5 <- gene.coverage.plot(pl=pl,samples = c("H3","H3_Depl"),genes=genes,plot.coverage = F,only.gene = T,
                             mygene = "RPL8A",ori = "sense",cols = cols,anno = anno,Margin = 500)
    p6 <- gene.coverage.plot(pl=pl,samples = c("H3","H3_Depl"),genes=genes,plot.coverage = F,only.gene = T,
                             mygene = "RPL28",ori = "sense",cols = cols,anno = anno,Margin = 500)
    pl.list <- list(p1[[1]],p2[[1]],p3[[1]],
                    p1[[2]],p2[[2]],p3[[2]],
                    p4[[1]],p5[[1]],p6[[1]],
                    p4[[2]],p5[[2]],p6[[2]])
  }
  return(pl.list)
  ##
  #return(c(p1,p2,p3,p4,p5,p6)) 
}



##unfinished
CellGrowthSignature <- function(ODuibhir.table.paths=c("A:/work/yeast_Annotations/table1.RData","A:/work/yeast_Annotations/table3.RData"),
                                DEExpression.filepath="A:/work/WinstonLab/Natalia_Reim/ThesisFigures_March2019/NR_March19_SenseDE.Data.RData",
                                plot.test.examples=F,
                                threshold=0.1,
                                OUTDIR="A:/work/WinstonLab/Natalia_Reim/ThesisFigures_March2019",
                                OUTPREFIX="NR_March19"){
  require(ggpubr)
  load(DEExpression.filepath) ## res
  load(ODuibhir.table.paths[1]) #table1
  load(ODuibhir.table.paths[2]) #table3
  
  de.genes <- res[!is.na(res$padj) & res$padj<threshold,]
  de.genes$gene <- ifelse(is.na(de.genes$gene),as.character(de.genes$tracking_id),as.character(de.genes$gene))
  de.genes <- as.character(de.genes$gene)
  
  t = as.matrix(table1)
  x = as.matrix(table3)
  most=order(t[1,] / t[1,1])[1]
  PCA1 = t[,most]
  o1=order(cor(x,PCA1),decreasing=T)
  o2=order(t[,1] / t[1,1])
  t=t[o2,o1]
  x=x[o2,o1]
  PCA1 <- PCA1[rownames(x)]
  
  if(plot.test.examples==T){
    pdf(file=paste0(OUTDIR,"/",OUTPREFIX,"_Oduibhir_PosExamples.pdf"),width = 7,height = 6)
    par(mfrow=c(3,4))
    for (i in c(1, 136, 271, 405, 540, 675, 810, 945, 1080, 1214, 1349, 1484)){
      plot(PCA1,x[,i],xlab="PC1",ylab="log2 FC",pch=20,col=rgb(0,0,0,0.01), main=paste(colnames(x)[i]))
      points(PCA1[names(PCA1)%in%de.genes],x[rownames(x)%in%de.genes,i],pch=20,col=rgb(1,0,0,0.01))
      r1 <- lm(x[,i]~PCA1,na.action = na.omit)
      abline(v=0,h=0,col="gray10",lty=2)
      abline(r1,col="black")
      p <- round(cor(x[,i],PCA1,use = "complete.obs"),3)
      mtext(paste("r=",p,sep=""),3,0,col="black")
      p <- round(cor(PCA1[names(PCA1)%in%de.genes],x[rownames(x)%in%de.genes,i],use = "complete.obs"),3)
      mtext(paste("r=",p,sep=""),3,-1,col="red")
      rm(p,r1)
      #smoothScatter(PCA1$PCA1,x[,i],pch="",xlab="PC1",ylab="log2 FC",main=paste(colnames(x)[i]))
      #r1 <- lm(x[,i]~PCA1$PCA1,na.action = na.omit)
      #abline(r1,col="red")
      #abline(h = 0,v = 0,col="gray",lty=2)
      #mtext( paste("r=", round(cor(x[,i],PCA1$PCA1,use="complete.obs"),3),sep = ""),3,-1,col="red",cex=0.7)
    }
    dev.off()
  }
  
  PCA1 <- as.data.frame(PCA1)
  PCA1 <- cbind(PCA1,tracking_id=rownames(PCA1))
  names(PCA1)[2] <- "gene"
  res <- merge(res,PCA1,by="gene",all.x=T)
  
  
  res1 <- res[res$padj<threshold,]
  res1 <- res1[!is.na(res1$padj),]
  res2 <- res1[res1$log2FoldChange>0,]
  res <- subset(res,!rownames(res)%in%rownames(res1))
  res1 <- res1[res1$log2FoldChange<0,]
  
  
  cols=c("Upregulated"="#90C987","Downregulated"="#F6C141","Unchanged"="gray80")
  ## MA-plot
  pl1 <- res[,c("tracking_id","baseMean","log2FoldChange","padj")]
  corval <- cor(res1$log2FoldChange, res1$PCA1,use = "complete.obs")
  
  p<- ggplot(res) + geom_point(aes(PCA1,log2FoldChange),col="gray80",shape=21,size=3)+
    geom_point(data = res2,aes(PCA1,log2FoldChange),col="#90C987",shape=21,size=3)+
    geom_point(data = res1,aes(PCA1,log2FoldChange),col="#F6C141",shape=21,size=3)+
    xlab("1st PCA, O'Duibhir et al") + ylab("log2 fold change")+
    geom_hline(yintercept = 0,col="gray",linetype="dashed")+
    geom_vline(xintercept = 0,col="gray",linetype="dashed")+
    annotate(geom="text", x=2, y=2.5, label=paste0("r=",round(corval,3)),color="red",size=5,fontface="italic")+
    gg+ theme(axis.line = element_line(size=0.75,color="black"),
               axis.ticks = element_line(size=1,colour = "black"))
  p
}

GOAnalysis.figures <- function(DEExpression.filepath="A:/work/WinstonLab/Natalia_Reim/ThesisFigures_March2019/NR_March19_SenseDE_2Controls.Data.RData",
                              cntdata.path="A:/work/WinstonLab/Natalia_Reim/Analysis_Feb18/Sense_Antisese_ReadCounts_revised.RData",
                              threshold=0.05,
                              OUTDIR="A:/work/WinstonLab/Natalia_Reim/ThesisFigures_March2019",
                              OUTPREFIX="NR_March19"){
  
  if(F){
    load(cntdata.path)
    s <- s[complete.cases(s),]
    s1 <- s[s$chr%in%chr.flt,]
    s1$chr <- NULL
    s1 <- apply(s1,2,function(x){
      h <- quantile(x,0.25)
      return(ifelse(x>h,1,0))
    })
    l <- rowSums(s1)
    expressed.genes <- unique(rownames(s1)[which(l>0)])
    ## total 5521 expressed genes
    rm(s1,l)
    g <- genes[genes$species=="Scer" & genes$verification=="Verified" & !is.na(genes$verification),]
    expressed.genes <- subset(expressed.genes,expressed.genes%in%g$tracking_id)
    expressed.genes <- unique(expressed.genes)
    
    load(DEExpression.filepath)
    names(res)[c(2,4)] <- c("baseMean","padj")
    res1 <- res[res$log2FoldChange<0 & res$padj<threshold & !is.na(res$padj),]
    res2 <- res[res$log2FoldChange>0 & res$padj<threshold & !is.na(res$padj),]
    z <- GOStats.enrichment(genes = res1$tracking_id,universe = expressed.genes,p.val = 0.1,species = "Scer" )
    z1 <- GOStats.enrichment(genes = res2$tracking_id,universe = expressed.genes,p.val = 0.1 ,species = "Scer")
    write.table(z,file=paste0(OUTDIR,"/",OUTPREFIX,"_GOEnrich.txt"),sep="\t",append = F,quote = F,row.names = F)
    write.table(z1,file=paste0(OUTDIR,"/",OUTPREFIX,"_GOEnrichUp.txt"),sep="\t",append = F,quote = F,row.names = F)
    
  }
  
  ### plots
  library(ggplot2)
  library(ggrepel)
  one.data <- rbind(cbind(read.csv(paste0(OUTDIR,"/REVIGO.csv"),header=T),class="BP",status="Down"),
                    cbind(read.csv(paste0(OUTDIR,"/REVIGO (1).csv"),header=T),class="CC",status="Down"),
                    cbind(read.csv(paste0(OUTDIR,"/REVIGO (2).csv"),header=T),class="MF",status="Down"),
                    cbind(read.csv(paste0(OUTDIR,"/REVIGO (3).csv"),header=T),class="BP",status="UP"))
  
  one.data$plot_size <- log2(2.718281828459^one.data$plot_size)
  one.data$frequency <- gsub("%","",one.data$frequency)
  one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
  one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
  one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
  one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
  one.data$frequency <- as.numeric( as.character(one.data$frequency) );
  one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
  one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
  names(one.data)[7] <- "log10_p_value"
  
  cats <- c("gene expression",
            "DNA replication",
            "ncRNA metabolic process",
            "RNA processing",
            "nitrogen compound metabolic process",
            "establishment of ribosome localization",
            "nuclear transport",
            "cellular component organization or biogenesis",
            "ribonucleoprotein complex biogenesis",
            "ribosome biogenesis",
            "ribonucleoprotein complex localization",
            "nitrogen compound transport",
            "RNA modification")
  
  plot.list=list()
  pl <- one.data[one.data$class=="BP" & one.data$status=="Down",]
  names(pl)[6] <- "log2 (size)"
  names(pl)[7] <- "log10 (p-value)"
  p1 <- ggplot( data = pl );
  p1 <- p1 + geom_point( aes_string( x = "plot_X",y =  "plot_Y", colour = "`log10 (p-value)`", size = "`log2 (size)`"), alpha = I(0.4) ) + scale_size_area();
  p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "orange", "red"), limits = c( min(pl$`log10 (p-value)`), 0) );
  p1 <- p1 + scale_size( range=c(5, 20)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
  ex <- pl[pl$description%in%cats,]
  p1 <- p1 + geom_text_repel( data = ex, aes(plot_X, plot_Y, label = description), box.padding = unit(0.75, "lines"), point.padding=unit(0.75, "lines"), size = 6.1 );
  p1 <- p1 + labs (x = "semantic space x", y = "semantic space y");
  p1 <- p1 + theme(axis.text.y = element_text(color = "black",size=16),axis.title=element_text(size=20))
  p1 <- p1 + theme(axis.text.x = element_text(color = "black",size=16),axis.title=element_text(size=20))
  p1 <- p1 + theme(legend.key = element_blank(),legend.title = element_text(size=18,face="bold"),legend.text = element_text(colour="black", size = 16, face = "bold"))
  p1 <- p1 + ggtitle("Biological processes downregulated after Spn1 depletion") + theme(plot.title = element_text(size = 22,hjust=0.5))
  p1
  plot.list[["down_bp"]] <- p1
  
  pl <- one.data[one.data$class=="CC" & one.data$status=="Down",]
  names(pl)[6] <- "log2 (size)"
  names(pl)[7] <- "log10 (p-value)"
  p1 <- ggplot( data = pl );
  p1 <- p1 + geom_point( aes_string( x = "plot_X",y =  "plot_Y", colour = "`log10 (p-value)`", size = "`log2 (size)`"), alpha = I(0.4) ) + scale_size_area();
  p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "orange", "red"), limits = c( min(pl$`log10 (p-value)`), 0) );
  p1 <- p1 + scale_size( range=c(5, 20)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
  ex <- pl[pl$`log10 (p-value)`< -2,]
  p1 <- p1 + geom_text_repel( data = ex, aes(plot_X, plot_Y, label = description), box.padding = unit(0.75, "lines"), point.padding=unit(0.75, "lines"), size = 6.1 );
  p1 <- p1 + labs (x = "semantic space x", y = "semantic space y");
  p1 <- p1 + theme(axis.text.y = element_text(color = "black",size=16),axis.title=element_text(size=20))
  p1 <- p1 + theme(axis.text.x = element_text(color = "black",size=16),axis.title=element_text(size=20))
  p1 <- p1 + theme(legend.key = element_blank(),legend.title = element_text(size=18,face="bold"),legend.text = element_text(colour="black", size = 16, face = "bold"))
  p1 <- p1 + ggtitle("Cellular component downregulated after Spn1 depletion") + theme(plot.title = element_text(size = 22,hjust=0.5))
  p1
  plot.list[["down_cc"]] <- p1
  
  pl <- one.data[one.data$class=="MF" & one.data$status=="Down",]
  names(pl)[6] <- "log2 (size)"
  names(pl)[7] <- "log10 (p-value)"
  p1 <- ggplot( data = pl );
  p1 <- p1 + geom_point( aes_string( x = "plot_X",y =  "plot_Y", colour = "`log10 (p-value)`", size = "`log2 (size)`"), alpha = I(0.4) ) + scale_size_area();
  p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "orange", "red"), limits = c( min(pl$`log10 (p-value)`), 0) );
  p1 <- p1 + scale_size( range=c(5, 20)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
  ex <- pl[pl$`log10 (p-value)`< -2,]
  p1 <- p1 + geom_text_repel( data = ex, aes(plot_X, plot_Y, label = description), box.padding = unit(0.75, "lines"), point.padding=unit(0.75, "lines"), size = 6.1 );
  p1 <- p1 + labs (x = "semantic space x", y = "semantic space y");
  p1 <- p1 + theme(axis.text.y = element_text(color = "black",size=16),axis.title=element_text(size=20))
  p1 <- p1 + theme(axis.text.x = element_text(color = "black",size=16),axis.title=element_text(size=20))
  p1 <- p1 + theme(legend.key = element_blank(),legend.title = element_text(size=18,face="bold"),legend.text = element_text(colour="black", size = 16, face = "bold"))
  p1 <- p1 + ggtitle("Molecular function downregulated after Spn1 depletion") + theme(plot.title = element_text(size = 22,hjust=0.5))
  p1
  plot.list[["down_mf"]] <- p1
  
  pl <- one.data[one.data$class=="BP" & one.data$status=="UP",]
  names(pl)[6] <- "log2 (size)"
  names(pl)[7] <- "log10 (p-value)"
  p1 <- ggplot( data = pl );
  p1 <- p1 + geom_point( aes_string( x = "plot_X",y =  "plot_Y", colour = "`log10 (p-value)`", size = "`log2 (size)`"), alpha = I(0.4) ) + scale_size_area();
  p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "orange", "red"), limits = c( min(pl$`log10 (p-value)`), 0) );
  p1 <- p1 + scale_size( range=c(5, 20)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
  ex <- pl[pl$`log10 (p-value)`< -2,]
  p1 <- p1 + geom_text_repel( data = ex, aes(plot_X, plot_Y, label = description), box.padding = unit(0.75, "lines"), point.padding=unit(0.75, "lines"), size = 6.1 );
  p1 <- p1 + labs (x = "semantic space x", y = "semantic space y");
  p1 <- p1 + theme(axis.text.y = element_text(color = "black",size=16),axis.title=element_text(size=20))
  p1 <- p1 + theme(axis.text.x = element_text(color = "black",size=16),axis.title=element_text(size=20))
  p1 <- p1 + theme(legend.key = element_blank(),legend.title = element_text(size=18,face="bold"),legend.text = element_text(colour="black", size = 16, face = "bold"))
  p1 <- p1 + ggtitle("Biological processes upregulated after Spn1 depletion") + theme(plot.title = element_text(size = 22,hjust=0.5))
  p1
  plot.list[["up_bp"]] <- p1
  
  return(plot.list)
}


Spn1Occupancy_changes <- function(merged.metagene="A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/MergedMetagenes_ChIPSpike_final_unscaled.RData",
                                  viability.file="A:/work/WinstonLab/Natalia_Reim/ThesisFigures_March2019/NR_March19_Viability.txt"){
  load(merged.metagene)
  require(ggplot2)
  
  pl1 <- pl[pl$factor=="Spn1 (R2)",]
  pl1$condition <- factor(pl1$condition, levels=c("non-depleted","Spn1-depleted"))
  
  cols=c("Spn1-depleted" = "#BB5566","non-depleted" = "#4477AA",
         "DMSO"= "#4477AA","IAA"="#BB5566")
  
  q <- ggplot(pl1, aes(x=pos,y=signal,ymax=signal+se, ymin=signal-se,col=condition))
  q <- q + geom_vline(xintercept = c(50),colour="gray50", linetype = "dashed")
  q <- q + geom_linerange(col='gray') + geom_line(size=1.2) #+ ylim(-0.03,0.15) #ylim(-0.01,0.045) #ylim(-0.04,0.1)
  q <- q + scale_color_manual(values = cols)
  q <- q + scale_y_continuous(expand = c(0,0),limits=c(0,max(pl1$signal)+1))
  q <- q + scale_x_continuous(expand = c(0,0),breaks = c(0,50,150,250,350,450), 
                              labels=c("","TSS","1kb","2kb","3kb","4kb")) +
    ylab("normalized counts") + xlab("") 
  q <- q + guides(color=guide_legend(title="Spn1 levels")) 
  q <- q + gg
  q <- q + ggtitle("Spn1 occupancy by ChIP-seq")+
       theme(plot.title = element_text(color="black",size=20,vjust=0.5,hjust = 0),
          axis.line = element_line(colour = "black",size = 0.7),
          axis.ticks = element_line(colour = "black",size = 0.5),
          legend.position = c(1, 0), 
          legend.justification = c( 1,0),
          legend.title = element_blank(),
          legend.text = element_text(color="black",size=16),
          legend.background = element_blank())
  q
  
  z <- read.delim(viability.file,header=T)
  z$mean <- apply(z[,3:4],1,mean)
  z$sd <- apply(z[,3:4],1,sd)
  z$se <- apply(z[,3:4],1,function(x){sd(x)/sqrt(length(x))})
  z$ci <- z$se*1.96
  z$genotype <- gsub("SPN","Spn",z$genotype)
  z$genotype <- factor(z$genotype,levels=c("WT","noAID","Spn1-AID"))
  
  plot.list <- list()
   
  z <- z[,c("genotype","treatment","rep1","rep2")]  
  z <- reshape::melt(z,measure.vars=c("rep1","rep2"))
  
  plot.list[["viability"]]  <- ggplot(z, aes(genotype,value,fill = treatment))+
      geom_bar(position = 'dodge', stat = 'summary', fun.y = 'mean') +
      geom_errorbar(stat = 'summary', width = 0.25,position=position_dodge(.9),
                    fun.ymin = function(x) mean(x) - sd(x),
                    fun.ymax = function(x) mean(x) + sd(x)) +
      geom_point(aes(x = genotype), position = position_dodge(width = 0.9)) +gg+
      xlab("")+ylab("colonies per mL\n normalized by OD")+
      scale_fill_manual(values = cols)+
      theme(axis.line = element_line(colour = "black",size = 0.7),
            axis.ticks = element_line(colour = "black",size = 0.5),
            plot.title = element_text(color="black",size=20,vjust=0.5,hjust = 0),
            strip.text = element_blank(),
            legend.title = element_blank())+
      ggtitle("Viability after Spn1 depletion")
    
    
  plot.list[["chip_metagene"]] <- q
  return(plot.list)
}


RNASeq.sampleSimilarity.n.Spikein <- function(cntdata.path="A:/work/WinstonLab/Natalia_Reim/Analysis_Feb18/Sense_Antisese_ReadCounts_revised.RData", 
                                              gene.anno.path="A:/work/yeast_Annotations/ScSp_SGD.RevisedGenesFull.RData"){
  require(gplots)
  require(DESeq2)
  load(cntdata.path)
  load(gene.anno.path)
  s <- subset(s,rownames(s)%in%genes[genes$verification=="Verified",]$tracking_id)
  as <- subset(as,rownames(as)%in%genes[genes$verification=="Verified",]$tracking_id)
  s <- s[complete.cases(s),]
  
  df1.e <- subset(s, !s$chr%in%spiked.chr)
  df1.e <-as.matrix( round (df1.e[,1:12] ))
  type <- colnames(df1.e)
  condition <- sub("-[012]","",type)
  colData <- data.frame(condition,type)
  row.names(colData) <- colData$type
  dds <- DESeqDataSetFromMatrix(countData = df1.e,colData = colData,design = ~ condition)
  dds <- DESeq(dds)
  rld <- rlog(dds)
  
  dds1 <- DESeqDataSetFromMatrix(countData = df1.e,colData = colData,design = ~ condition)
  sizeFactors(dds1) <- dd[1:12]
  dds1 <- DESeq(dds1)
  rld1 <- rlog(dds1)
  
  
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  distsRL <- dist(t(assay(rld)))
  mat <- as.matrix(distsRL)
  distsRL1 <- dist(t(assay(rld)))
  mat1 <- as.matrix(distsRL1)
  
  rownames(mat) <- colnames(mat) <- rownames(mat1) <- colnames(mat1) <- with(colData(dds),paste(type))
  hc <- hclust(distsRL)
  #heatmap(mat, Rowv=as.dendrogram(hc),symm=TRUE, trace="none",col = rev(hmcol), margin=c(10, 10))
  
  
  load(file="A:/work/WinstonLab/Natalia_Reim/Analysis_Feb18/SpikeinNormalization_Assessment_RLE.RData")
  names(pl)[2] <- "Var2"
  pl$category <- gsub("SpikeNorm","normalized",pl$category)
  pl$category <- gsub("Raw","raw",pl$category)
  
  pl$category <- factor(pl$category,levels=c("raw Sp tags","raw Sc tags","normalized Sp tags","normalized Sc tags"))
  pl$Var2 <-  factor(pl$Var2,levels=rev(c("noAID_IAA-1", "noAID_IAA-2","SPN1-AID_DMSO-1", "SPN1-AID_DMSO-2", "SPN1-AID_IAA-1", "SPN1-AID_IAA-2"
  )))
  
  q <- ggplot(pl, aes(x=Var2,y=value,fill=Var2)) #+ geom_vline(xintercept = c(0.365,0.4))
  q <- q + geom_boxplot(outlier.shape = NA) + ylab("relative log2 expression") + 
    facet_wrap(~category,nrow=2) +coord_flip()+xlab("")
  q <- q + ylim(-1,1)
  q <- q + geom_hline(yintercept = 0,col="black")
  q <- q + scale_color_manual(values = cols) +scale_fill_manual(values=cols)
  q <- q + theme(panel.background = element_rect(fill = NA),
                 panel.grid.major = element_line(colour = "grey80"),
                 axis.line = element_line(size = 0.7, colour = "black"))
  q <- q 
  q <- q + theme(strip.text = element_text(color="black",size = 16))
  q <- q + theme(legend.position = "none")
  q <- q + theme(panel.grid = element_blank())
  q
  
  pdf(file=paste0(OUTDIR,"/",OUTPREFIX,"_Figure22.pdf"),width = 6,height = 6)
  heatmap.2(mat,notecol="black",notecex =1,density.info="none",trace="none",margins =c(10,10),col=rev(hmcol),sepwidth = c(0.01,0.01),key.title = "",key.xlab = "similarity",key.ylab = "") 
  heatmap.2(mat1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(10,10),col=rev(hmcol),sepwidth = c(0.01,0.01),key.title = "",key.xlab = "similarity",key.ylab = "") 
  dev.off()
  pdf(file=paste0(OUTDIR,"/",OUTPREFIX,"_Figure22_b1.pdf"),width = 9,height = 6)
  print(q)
  dev.off()
  
}



IntronRetention.plots <- function(inp.file="A:/work/WinstonLab/Natalia_Reim/Paper_Thesis_Figs/Natalia_IntronRetention_Jan4_2019.txt",
                                  gene.groups="A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Gene groups.RData",
                                  refs.rdata="A:/work/yeast_Annotations/refs.RData"){
  res <- read.delim(inp.file,header = T)
  res <- res[res$chr!="chrM",]
  if(F){
    write.table(res, file=paste0(OUTDIR,"/IntronRetentionQuantifications.txt"),sep = "\t",row.names = F,quote = F)
  }
  
  pl <- res
  pl$checksum <- pl$unspliced+pl$spliced
  p <- reshape2::dcast(pl,id~profile,value.var = "checksum")
  p[,2:13] <- apply(p[,2:13],2,function(x){
    ifelse(x<3,0,x)
  })
  p <- p[,c(1,6:9)]
  p$test <- p$`SPN1-AID_DMSO-1`*p$`SPN1-AID_DMSO-2`*p$`SPN1-AID_IAA-1`*p$`SPN1-AID_IAA-2`
  p <- p[p$test>0,]
  p <- p$id
  pl <- pl[pl$id%in%p,]
  
  
  #pl <- res[res$unspliced+res$spliced>=3,]
  #pl <- res
  #mergepl <- res
  mergepl <- pl 
  load(gene.groups)
  rp <- rp[rp$gene!="RPP1",]
  rp <- rp[-grep("MRP",rp$gene),]
  
  mergepl$profile <- gsub("-[12]","",mergepl$profile)
  mergepl <- data.table(mergepl)
  
  mergepl <- mergepl[,unspliced_5_m:=sum(FivePrimeUnspliced),by=list(chr,start,end,strand,reference,gene,profile,id,gene_width,intron_width,WT_RPKM,nondepleted_RPKM,depleted_RPKM)]
  mergepl <- mergepl[,unspliced_3_m:=sum(ThreePrimeUnspliced),by=list(chr,start,end,strand,reference,gene,profile,id,gene_width,intron_width,WT_RPKM,nondepleted_RPKM,depleted_RPKM)]
  mergepl <- mergepl[,spliced_m:=sum(spliced),by=list(chr,start,end,strand,reference,gene,profile,id,gene_width,intron_width,WT_RPKM,nondepleted_RPKM,depleted_RPKM)]
  mergepl <- mergepl[,internal_m:=sum(internal),by=list(chr,start,end,strand,reference,gene,profile,id,gene_width,intron_width,WT_RPKM,nondepleted_RPKM,depleted_RPKM)]
  mergepl <- mergepl[,gene.cov_m:=sum(gene.cov),by=list(chr,start,end,strand,reference,gene,profile,id,gene_width,intron_width,WT_RPKM,nondepleted_RPKM,depleted_RPKM)]
  mergepl <- mergepl[,intron.cov_m:=sum(intron.cov),by=list(chr,start,end,strand,reference,gene,profile,id,gene_width,intron_width,WT_RPKM,nondepleted_RPKM,depleted_RPKM)]
  mergepl <- mergepl[,unspliced_m:=sum(unspliced),by=list(chr,start,end,strand,reference,gene,profile,id,gene_width,intron_width,WT_RPKM,nondepleted_RPKM,depleted_RPKM)]
  mergepl<- as.data.frame(mergepl)
  
  mergepl <- unique(mergepl[,c("chr", "start", "end", "strand", "reference", "gene", "id", 
                               "gene_width", "intron_width", "profile", "WT_RPKM", "nondepleted_RPKM", 
                               "depleted_RPKM","unspliced_5_m", "unspliced_3_m", "unspliced_m", 
                               "spliced_m", "internal_m", "intron.cov_m", "gene.cov_m")])
  
  #### Calculate intron retention ratios
  ##  IRR1 = (unspliced+internal)/(unspliced+internal+spliced) 
  mergepl$IRR1 <- round((mergepl$unspliced_m+mergepl$internal_m)/(mergepl$unspliced_m+mergepl$internal_m+mergepl$spliced_m),2)
  ##  IRR2 = (unspliced/2+internal)/(unspliced/2+internal+spliced) 
  mergepl$IRR2 <- round(( (mergepl$unspliced_m/2)+mergepl$internal_m)/( (mergepl$unspliced_m/2)+mergepl$internal_m+mergepl$spliced_m),2)
  ##  IRR3 = (unspliced/2)/(unspliced/2+spliced) 
  mergepl$IRR3 <- round(( (mergepl$unspliced_m/2))/( (mergepl$unspliced_m/2)+mergepl$spliced_m),2)
  ##  IRR4 = len.norm.intron.cov /len.norm.gene.cov
  mergepl$IRR4 <- round ( ( (mergepl$intron.cov_m)/mergepl$intron_width) / (mergepl$gene.cov_m/mergepl$gene_width), 2)
  
  write.table(mergepl,file="A:/work/WinstonLab/Natalia_Reim/Paper_Thesis_Figs/Natalia_IntronRetention_MergedReplicates_May15_2019.txt",sep = "\t",quote = F,row.names = F)
  
  pl$profile <- factor(pl$profile,levels=c("SPN1-AID_DMSO-1", "SPN1-AID_DMSO-2", "SPN1-AID_IAA-1", "SPN1-AID_IAA-2",
                                           "noAID_DMSO-1", "noAID_DMSO-2", "noAID_IAA-1", "noAID_IAA-2", 
                                           "WT_DMSO-1", "WT_DMSO-2", "WT_IAA-1", "WT_IAA-2"))
  mergepl$profile <- factor(mergepl$profile,levels=c("SPN1-AID_DMSO", "SPN1-AID_IAA","noAID_DMSO", "noAID_IAA",  
                                                     "WT_DMSO", "WT_IAA"))
  load(refs.rdata)
  
  #pl <- pl[(pl$spliced+(pl$unspliced/2))>=3,]
  #mergepl <- mergepl[(mergepl$spliced_m+(mergepl$unspliced_m/2))>=3,]
  
  myIRR= "IRR1"
  #for(myIRR in c("IRR1")){
    cat(myIRR,"\n")
    m <- pl[,c("reference","id","profile",myIRR)]  
    n <- mergepl[,c("reference","id","profile",myIRR)]  
    names(m)[4] <- "IRR"
    names(n)[4] <- "IRR"
    if(myIRR=="IRR4"){
      m <- m[is.finite(m$IRR),]
      n <- n[is.finite(n$IRR),]
    }else{
      m$IRR <- ifelse(!is.finite(m$IRR),0,m$IRR)
      n$IRR <- ifelse(!is.finite(n$IRR),0,n$IRR)
    }
    
    m <- reshape2::dcast(m,id+reference~profile,value.var = "IRR")
    n <- reshape2::dcast(n,id+reference~profile,value.var = "IRR")
    mscatter = m
    nscatter = n
    n <- merge(n,refs[,c("tracking_id","gene","Gene.type")],by.x="reference",by.y="tracking_id",all.x=T)
    n$category <- "Non-RPL introns"
    n[grep("RP[SLP]\\S*",n$gene),]$category <- "RPL introns"
    n <- n[,c(1,2,9:11,4,3,5:8)]
    n$control_max <- apply(n[,7:11],1,max)
    
    
    nx <- n[,c("reference","category","SPN1-AID_IAA","SPN1-AID_DMSO")]
    nx$category <- "all introns"
    nx <- rbind(nx,n[,c("reference","category","SPN1-AID_IAA","SPN1-AID_DMSO")])
    nx <- nx[complete.cases(nx),]
    
    
    numintrons <- c()
    numgenes <- c()
    for(y in c("all introns","RPL introns", "Non-RPL introns")){
      y <- as.character(y)
      numintrons <- c(numintrons,sum(nx$category==y),
                      sum(nx[nx$reference%in%de.genes$tracking_id,]$category==y),
                      sum(nx[nx$`SPN1-AID_IAA`<nx$`SPN1-AID_DMSO` ,]$category==y),
                      sum(nx[nx$`SPN1-AID_IAA`>nx$`SPN1-AID_DMSO` ,]$category==y),
                      sum(nx[nx$`SPN1-AID_IAA`<nx$`SPN1-AID_DMSO` & nx$reference%in%de.genes$tracking_id,]$category==y),
                      sum(nx[nx$`SPN1-AID_IAA`>nx$`SPN1-AID_DMSO` & nx$reference%in%de.genes$tracking_id,]$category==y))
      
      numgenes <- c(numgenes,length(unique(nx[nx$category==y,]$reference)),
                    length(unique(nx[nx$reference%in%de.genes$tracking_id & nx$category==y,]$reference)),
                    length(unique(nx[nx$`SPN1-AID_IAA`<nx$`SPN1-AID_DMSO` & nx$category==y,]$reference)),
                    length(unique(nx[nx$`SPN1-AID_IAA`>nx$`SPN1-AID_DMSO` & nx$category==y,]$reference)),
                    length(unique(nx[nx$`SPN1-AID_IAA`<nx$`SPN1-AID_DMSO` & nx$reference%in%de.genes$tracking_id & nx$category==y,]$reference)),
                    length(unique(nx[nx$`SPN1-AID_IAA`>nx$`SPN1-AID_DMSO` & nx$reference%in%de.genes$tracking_id &nx$category==y,]$reference)))
      
    }
    resview <- data.frame(intron_type=c(rep("all introns",6),rep("RPL introns",6),rep("Non-RPL introns",6)),
                          filtered_for = rep(c("-","Downregulated genes","Decreased IRR","Increased IRR","Decreased IRR+Dwonregulated genes","Increased IRR+Downregulated genes"),3),
                          num_of_introns=numintrons,
                          num_of_genes = numgenes)
    
    
    for(i in levels(as.factor(nx$category))){
      #pval = ks.test(nx[nx$category==i,]$`SPN1-AID_IAA`,nx[nx$category==i,]$`SPN1-AID_DMSO`,exact = T)$p.value
      pval = wilcox.test(nx[nx$category==i,]$`SPN1-AID_IAA`, 
                         nx[nx$category==i,]$`SPN1-AID_DMSO`, 
                         paired = T, alternative = "greater")$p.value
      cat(i, " ", pval,"\n")
    }
    rr <- as.data.frame(table(as.character(nx$category)))
    names(rr) <- c("category","num")
    nx <- merge(nx,rr,by="category",all.x=T)
    nx$category <- paste0(nx$category," (",nx$num,")")
    nx$num <- NULL 
    nx$reference <- NULL
    
    nx$diff <- nx$`SPN1-AID_IAA`-nx$`SPN1-AID_DMSO`
    
    nx <- melt(nx)
    names(nx) <- c("category","profile","IRR")
    nx$profile <- gsub("SPN1-AID_DMSO","non-depleted",nx$profile)
    nx$profile <- gsub("SPN1-AID_IAA","Spn1-depleted",nx$profile)
    nx <- nx[nx$profile=="diff",]
    #col <- c("non-depleted"="#66CCEE","Spn1-depleted"="#BB5566")
    nx$category <- gsub("Non-RPL","non-RP",nx$category)
    nx$category <- gsub("RPL","RP",nx$category)
    
    my_comparisons <- list( c("all introns (144)","non-RP introns (53)"),
                            c("all introns (144)","RP introns (91)"),
                            c("non-RP introns (53)","RP introns (91)") )
    
    p3 <- ggplot(nx, aes(x=category,y=IRR,fill=profile)) + geom_violin(col='black') +
      geom_boxplot(width=0.1,fill='white', color="black",outlier.shape = NA,outlier.size = 0)+
      xlab("") + ylab("change in \nintron retention ratio") +
      scale_fill_manual(values = c("diff"="#dd8033"))+
      stat_compare_means(comparisons = my_comparisons,
                         label.x = 1.5,na.rm=T,method.args = list(alternative="two.sided"),
                         paired = F,aes(label=..p.signif..))+
      geom_hline(yintercept = 0,col="gray")+
      theme(panel.background = element_rect(fill = "white", colour = "black")) + 
      theme(axis.text = element_text(colour = 'black',size = 14),
            axis.text.x = element_text(angle=30,hjust=1),
            axis.title = element_text(color = "black",size=18),
            strip.background = element_blank(),
            strip.text = element_text(color = "black",size=18),
            legend.title = element_blank(),
            legend.background = element_rect(color = NA),
            legend.box.background = element_rect(color = NA),
            legend.text = element_text(size=14,color="black"))+
      theme(panel.grid = element_blank())+
      theme(legend.position = "none")
    
    pdf(file = paste0(OUTDIR,"/",OUTPREFIX,"_IntronRetentionRatioDifference.pdf"),height = 6,width = 6)
    p3
    dev.off()
    
  
}


shiftratio_ecdf <- function(shift.file="A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/RawCounts_ShiftRatioNucls_SGDFilt.RData",
                            key.file="A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/SampleKey.txt",
                            gene.path="A:/work/yeast_Annotations/ScSp_SGD.RevisedGenesFull.RData",
                            length_cutoff=1000){
  
  load(shift.file)
  keys <- read.delim(key.file)
  load(gene.path)
  
  genes <- genes[genes$chr%in%c(chr.flt),]
  genes <- genes[!is.na(genes$verification),]
  genes <- as(genes,"GRanges")
  d <- as.data.frame(findOverlaps(genes,select="all")) ## remove genes that are overlapping on the same strands
  d <- d[d$queryHits!=d$subjectHits,]
  d <- unique(c(d$queryHits,d$subjectHits))
  genes <- genes[-d]
  genes <- genes[genes$verification=="Verified"]
  
  genes <- subset(genes,width(genes)>=length_cutoff)
  
  m <- subset(m,m$tracking_id%in%as.character(genes$tracking_id))
  m <- subset(m,m$chr%in%c(chr.flt))
  m <- m[grep("H3",m$profile),]
  
  rm(d)
  d <- c()
  for(i in 1:length(keys$ChIP)){
    if(as.character(keys$ChIP[i])%in%m$profile){
      d <- rbind(d,
                 cbind(m[m$profile==as.character(keys$ChIP[i]),],
                       cond=paste(as.character(keys$Factor[i]),"",sep="")))
    }
  }
  d$replicate <- "Non-depleted"
  d[grep(".Depl",d$profile),]$replicate <- "Depleted"
  d$replicate <- paste(d$replicate," (",gsub("\\S*_","",d$profile),")",sep="")
  d$shift.ratio <- round(d$down.500/d$up.500,2)
  d$shift.ratio <- ifelse(!is.finite(d$shift.ratio),NA,d$shift.ratio)
  d$cond <- paste0(d$cond," (R2)")
  d$cond <- gsub("H3 [(]round1[)] [(]R2[)]","H3 (R1)",d$cond)
  d$cond <- gsub("H3K36me3 [(]R2[)]","H3K36me3 (R1)",d$cond)
  table(d$cond)
  d <- d[d$cond!="H3 (R1)",]
  d$cond <- factor(d$cond,levels=c("H3 (R2)", "H3K36me3 (R1)","H3K36me2 (R2)", 
                                    "H3K4me3 (R2)"))
  
  q <- ggplot(d, aes(shift.ratio,col=replicate))
  q <- q + stat_ecdf(geom = "line") + xlim(0,5)
  q <- q + scale_color_manual(values = cols)
  q <- q + facet_wrap(~cond, nrow=1,scales = "free")
  q <- q + ylab("% of genes") + xlab("3' Signal Shift")
  q <- q + gg
  q <- q + theme(strip.text.x = element_text(colour = 'black',size=18),
                 legend.position = "none")
  
  d$replicate <- factor(d$replicate,levels=c("Non-depleted (1)", "Non-depleted (2)", "Non-depleted (3)", "Non-depleted (4)",
                                             "Depleted (1)", "Depleted (2)", "Depleted (3)", "Depleted (4)"))
  
  p <- ggplot(d, aes(x=replicate,y=log2(shift.ratio),fill=replicate)) #+ geom_vline(xintercept = c(0.365,0.4))
  p <- p + geom_boxplot(outlier.shape = NA) + ylab("shift, log2") + xlab("")+
           facet_wrap(~cond,nrow=1,scales = "free")
  p <- p + geom_hline(yintercept = 0,col="black")
  p <- p + scale_color_manual(values = cols) +scale_fill_manual(values=cols)
  p <- p + ylim(-2.5,2.5)
  p <- p + gg
  p <- p + theme(strip.text.x = element_text(colour = 'black',size=18),
                 legend.position = "bottom",
                 axis.text.x = element_blank())
  p
  
  plist <- list()
  plist[["ECDF"]] <- q
  plist[["Box"]] <- p
  return(plist)
}

my.scatterplot.fun <- function(l2fc,xlab,ylab,x,y,main,side="right",pos="bottom",gg){
  q <-ggplot(l2fc, aes(x=l2fc[,x], y=l2fc[,y])) + 
    xlab(paste0(xlab))+ylab(paste0(ylab))+ggtitle(paste(main))+
    stat_bin_hex(geom="point", aes(color=log2(..count..)), binwidth=c(.05,.05), size=2, shape=16, stroke=0)+
    scale_color_viridis(option="inferno",direction = 1)+
    geom_hline(yintercept = 0,col="black")+geom_vline(xintercept = 0,col="black")+
    geom_smooth(method='lm',formula=y~x,col="red")+
    stat_cor(method = "pearson", size=5,col="red",label.sep = "\n",
             label.x.npc = side, label.y.npc = pos)+gg+
    scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+ylim(-2,3)
  return(q)
}


scatters_fig2_7 <- function(df.path="A:/work/WinstonLab/Natalia_Reim/DataFrames/NataliaChIPOccupShiftsExpression.RData"){
  load(df.path)
  gg1 <- theme_minimal() +theme(axis.text.x = element_text(size=14,colour = "black"),
                               axis.text.y = element_text(size=14,colour = "black"),
                               axis.title =  element_text(size=16,color="black"),
                               axis.line = element_line(color="black",size = 0.75),
                               axis.ticks = element_line(color="black",size = 1),
                               legend.position = "none",
                               plot.title = element_text(colour ="black",size = 18,hjust = 0.5))
  pl.list <- list()
  bs <- "basal_rna_exp"
  qs <- c("shift_H3 (R2)", "shift_H3K36me3 (R1)","shift_H3K36me2 (R2)", "shift_H3K4me3 (R2)")
  for(i in 1:length(qs)){
    xlab="basal expression, log2"
    ylab="shift,log2"
    x=bs
    y=qs[i]
    main=gsub("shift_","",qs[i])
    cat(main,"\n")
    pl.list[[paste0(main, "Expn")]] <- print(my.scatterplot.fun(pl,xlab=xlab,ylab=ylab,main=main,x=x,y=y,
                                                                side="right",pos = "bottom",gg1) )
    rm(xlab,ylab,main,x,y)
  }
  
  
  bs <-"occup_RNAPII (R2)"
  qs <- c("shift_H3 (R2)", "shift_H3K36me3 (R1)","shift_H3K36me2 (R2)", "shift_H3K4me3 (R2)")
  for(i in 1:length(qs)){
    xlab="Rpb1 occupancy, log2"
    ylab="shift,log2"
    x=bs
    y=qs[i]
    main=gsub("shift_","",qs[i])
    cat(main,"\n")
    pl.list[[paste0(main, "Occup")]] <- print(my.scatterplot.fun(pl,xlab=xlab,ylab=ylab,main=main,x=x,y=y,
                                                                 side="right",pos = "bottom",gg1) )
    rm(xlab,ylab,main,x,y)
  }
  
  return(pl.list)
}

scatters_fig2_3 <- function(df.path="A:/work/WinstonLab/Natalia_Reim/DataFrames/NataliaChIPOccupShiftsExpression.RData"){
  load(df.path)
  gg1 <- theme_minimal() +theme(axis.text.x = element_text(size=14,colour = "black"),
                                axis.text.y = element_text(size=14,colour = "black"),
                                axis.title =  element_text(size=16,color="black"),
                                axis.line = element_line(color="black",size = 0.75),
                                axis.ticks = element_line(color="black",size = 1),
                                legend.position = "none",
                                plot.title = element_text(colour ="black",size = 18,hjust = 0.5))
  pl.list <- list()
  bs <-"expn_l2fc"
  qs <- c("l2fc_RNAPII (R2)","l2fc_Ser5-P (R2)","l2fc_Ser2-P (R2)","l2fc_Spt6 (R1)")
  for(i in 1:length(qs)){
    xlab="log2FC, expression"
    ylab="log2FC, occupancy"
    x=bs
    y=qs[i]
    main=gsub("l2fc_","",qs[i])
    cat(main,"\n")
    pl.list[[paste0(main, "Occup")]] <- print(my.scatterplot.fun(pl,xlab=xlab,ylab=ylab,main=main,x=x,y=y,side="left",pos = "top",gg1) )
    rm(xlab,ylab,main,x,y)
  }
  return(pl.list)
}


## Figure1
if(F){
  p <- Spn1Occupancy_changes()
  lay <- rbind(c(1,2))
  
  pdf(file = paste0(OUTDIR,"/",OUTPREFIX,"_ThesisFig1.pdf"),height = 4,width = 10)
  grid.arrange(grobs = p, layout_matrix = lay,widths=c(5,4))
  dev.off()
}

### Figure2-2
if(F){
  p1 <- MAplot.PieChart.DEGenes(threshold = 0.1)
  p2 <- RTPCR_Barplot()
  p2 <- list(p2)
  #source("A:/work/scripts/Project_Winston/final/GeneCoveragePlot.R")
  ptmp <- RNASeq.SenseCoverage(script.path = "A:/work/scripts/Project_Winston/final/GeneCoveragePlot.R")
  gs <- c(p1,ptmp,p2)
  lay <- rbind(c(1,1,1,2,2,2),
               c(3,3,4,4,5,5),
               c(6,6,7,7,8,8),
               c(9,9,10,10,11,11),
               c(12,12,13,13,14,14),
               rep(15,6))
  pdf(file = paste0(OUTDIR,"/",OUTPREFIX,"_ThesisFig2.pdf"),height = 11,width = 9)
  grid.arrange(grobs = gs, heights=c(10,8,2.5,8,2.5,10), layout_matrix = lay)
  dev.off()
  #grid.arrange(grobs=p1,nrow=1,widths=c(5,3)) #,top= textGrob(paste0(mygene),gp=gpar(fontsize=22,font=3))
  #grid.arrange(grobs=p2,nrow=1,widths=c(4,4)) #,top= textGrob(paste0(mygene),gp=gpar(fontsize=22,font=3))
  #grid.arrange(grobs=ptmp,nrow=4,heights=c(8,1.5,8,1.5)) #,top= textGrob(paste0(mygene),gp=gpar(fontsize=22,font=3))
  save(gs,p1,p2,ptmp,file=paste0(OUTDIR,"/",OUTPREFIX,"_SenseRNASeqCoverage.PlotDFs.RData"))
  
}

## Figure 2.4
if(F){
  theme_set(theme_cowplot(font_size=20))
  p1 <- CellGrowthSignature()
  p1 <- list(p1)
  p2 <- GOAnalysis.figures()
  
  
  gs <- c(p1,p2[1])
  lay <- rbind(c(1,1,NA),
               c(2,2,2))
  pdf(file = paste0(OUTDIR,"/",OUTPREFIX,"_ThesisFig4.pdf"),height = 10,width = 8)
  grid.arrange(grobs = gs, heights=c(4,6),layout_matrix=lay)
  dev.off()
  
}

##Figure S2_3
if(F){
  ptmp <- RNASeq.SenseCoverage(gene.names =c("CUS2","MAF1","RRP1",
                                             "ARF2","BMT6","PRS1",
                                             "TFB4","STF1","TMC1",
                                             "CYB5","HMS1","MOH1"),which.cov = "antisense3")

  lay <- rbind(c(1,2,3),
               c(4,5,6),
               c(7,8,9),
               c(10,11,12),
               c(13,14,15),
               c(16,17,18),
               c(19,20,21),
               c(22,23,24))
  pdf(file = paste0(OUTDIR,"/",OUTPREFIX,"_ThesisFigS2_3.pdf"),height = 12,width = 12)
  grid.arrange(grobs = ptmp, heights=c(8,2,8,2,8,2,8,2), layout_matrix = lay)
  dev.off()
  
}

## Figure 2.6
if(F){
  p3 <- Splicing_Barplot()
  p1 <- IntronRetention.plots()
  p1 <- list(p1)
  p3 <- list(p3)
  
  #source("A:/work/scripts/Project_Winston/final/GeneCoveragePlot.R")
  ptmp <- RNASeq.SenseCoverage(gene.names =c("RPS21B","RPL43B","SFT1","RPL34A","RPL14B","RPS27A"),which.cov = "splicing")
  gs <- c(p1,ptmp,p3)
  lay <- rbind(c(1,1,1,1,1,1),
               c(3,3,4,4,5,5),
               c(6,6,7,7,8,8),
               c(9,9,10,10,11,11),
               c(12,12,13,13,14,14),
               c(NA,NA,NA,NA,NA,NA),
               c(15,15,15,15,15,15))
  pdf(file = paste0(OUTDIR,"/",OUTPREFIX,"_ThesisFig2_6.pdf"),height = 13,width = 10)
  grid.arrange(grobs = gs, heights=c(10,8,2.5,8,2.5,4,10), layout_matrix = lay)
  dev.off()
  #grid.arrange(grobs=p1,nrow=1,widths=c(5,3)) #,top= textGrob(paste0(mygene),gp=gpar(fontsize=22,font=3))
  #grid.arrange(grobs=p2,nrow=1,widths=c(4,4)) #,top= textGrob(paste0(mygene),gp=gpar(fontsize=22,font=3))
  #grid.arrange(grobs=ptmp,nrow=4,heights=c(8,1.5,8,1.5)) #,top= textGrob(paste0(mygene),gp=gpar(fontsize=22,font=3))
  save(gs,p1,p2,ptmp,file=paste0(OUTDIR,"/",OUTPREFIX,"_SenseRNASeqCoverage.PlotDFs.RData"))
  
}

## Figure 2.5
if(F){
  ptmp <- RNASeq.SenseCoverage(gene.names =c("CUS2","MAF1","RRP1","FUN26","GIM3","PRP9","SNU66","VAM6","STF1"),which.cov = "antisense")
  lay <- rbind(c(1,2,3),
               c(4,5,6),
               c(7,8,9),
               c(10,11,12),
               c(13,14,15),
               c(16,17,18))
  pdf(file = paste0(OUTDIR,"/",OUTPREFIX,"_ThesisFig2_5.pdf"),height = 12,width = 12)
  grid.arrange(grobs = ptmp, heights=c(8,2,8,2,8,2), layout_matrix = lay)
  dev.off()
  #grid.arrange(grobs=p1,nrow=1,widths=c(5,3)) #,top= textGrob(paste0(mygene),gp=gpar(fontsize=22,font=3))
  
  
  g0 <- c("CUS2","MAF1","RRP1","FUN26","GIM3","PRP9","SNU66","VAM6","STF1")
  g1 <- c("APP1","ARF2","BMT6","CUS2","FUN26","GIM3","MAF1","PRP9","PRS1","PUF6","RPL23A","RPP1B","RPS1B","RRI1","RRP1","SNU66","STF1","TFB4","TRR1","UMP1","VAM6","YCP4")
  g2 <- c("CYB5","HMS1","HOL1","MOH1","TMC1")
  g1 <- g1[!g1%in%g0]
  g2 <- g2[!g2%in%g0]
  g1 <- c(g1,g2)
  gene.names <- g1
  
  lay <- rbind(c(1,2,3),
               c(4,5,6),
               c(7,8,9),
               c(10,11,12),
               c(13,14,15),
               c(16,17,18))
  pdf(file = paste0(OUTDIR,"/",OUTPREFIX,"_ThesisFig2_5_additional1.pdf"),height = 12,width = 12)
  grid.arrange(grobs = ptmp[1:18], heights=c(8,2,8,2,8,2), layout_matrix = lay)
  dev.off()
  pdf(file = paste0(OUTDIR,"/",OUTPREFIX,"_ThesisFig2_5_additional2.pdf"),height = 12,width = 12)
  grid.arrange(grobs = ptmp[19:36], heights=c(8,2,8,2,8,2), layout_matrix = lay)
  dev.off()
  
  
}

## Figure2-3 suppliment
if(F){
  p1 <- MAplot.PieChart.DEGenes(DEExpression.filepath = "A:/work/WinstonLab/Natalia_Reim/ThesisFigures_March2019/NR_March19_SenseDE_Spn1vsnoAID.Data.RData")
  p2 <- CellGrowthSignature()
  p3 <- RNASeq.RLE.boxplot()
  
  gs <- c(list(p3),p1[1],list(p2))
  
  
  lay <- rbind(c(NA,1),
               c(2,3))
  pdf(file = paste0(OUTDIR,"/",OUTPREFIX,"_ThesisFig2_2supl.pdf"),height = 8,width = 10)
  grid.arrange(grobs = gs, layout_matrix = lay)
  dev.off()
  
  load("A:/work/WinstonLab/Natalia_Reim/ThesisFigures_March2019/NR_March19_SenseDE_Spn1vsnoAID.Data.RData")
  res$status <- "unchanged"
  res$status <- ifelse(res$log2FoldChange>0 & res$padj<0.1 & !is.na(res$padj), as.character("Up"),as.character(res$status))
  res$status <- ifelse(res$log2FoldChange<0 & res$padj<0.1 & !is.na(res$padj), as.character("Down"),as.character(res$status))
  table(res$status)
  write.table(res,file=paste0(OUTDIR,"/",OUTPREFIX,"Spn1IAA_vs_noAIDIAA_DEGenes.txt"),sep = "\t",quote = F,row.names = F)
   
}

## Figure 2_7
if(F){
  p1 <- shiftratio_ecdf()
  p2 <- scatters_fig2_7()
  
  gs <- c(p1[1],p1[2],p2[1],p2[2],p2[3],p2[4],p2[5],p2[6],p2[7],p2[8])
  lay <- rbind(c(1,1,1,1),
               c(2,2,2,2),
               c(3,4,5,6),
               c(7,8,9,10))
  
  pdf(file = paste0(OUTDIR,"/",OUTPREFIX,"_ThesisFig2_7.pdf"),height = 12,width = 12)
  grid.arrange(grobs = gs,heights=c(4,5,5,5), layout_matrix = lay)
  dev.off()
  
}

##Figure 2_3new
if(F){
  p1 <- ChIPSeq.Coverage(gene.names =c("UBI4","GCV3"))
  p2 <- scatters_fig2_3()
  
  gs <- c(p2[1],p2[2],p2[3],p2[4],p1)
  
          #p1[1],p1[2],p1[3],p1[4],p1[7],p1[8],p1[7],p1[8])
  
  
  lay <- rbind(c(1,2,3,4),
               c(5,6,7,8),
               c(9,10,11,12),
               c(13,14,15,16),
               c(17,18,19,20))
  pdf(file = paste0(OUTDIR,"/",OUTPREFIX,"_ThesisFig2_3new.pdf"),height = 8,width = 10)
  grid.arrange(grobs = gs, heights=c(10,8,2,8,2), layout_matrix = lay)
  dev.off()
  
}


##Figure 2_a1 #H3 coverages
if(F){
  p1 <- ChIPSeq.Coverage(gene.names =c("SML1","ILV5","TSA1","RPL10","RPL8A","RPL28"))
  lay <- rbind(c(1,2,3),
               c(4,5,6),
               c(7,8,9),
               c(10,11,12))
  pdf(file = paste0(OUTDIR,"/",OUTPREFIX,"_ThesisFig2_a1.pdf"),height = 6,width = 9)
  grid.arrange(grobs = p1, heights=c(8,2,8,2), layout_matrix = lay)
  dev.off()
}

