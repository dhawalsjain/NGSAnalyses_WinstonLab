rm(list=ls())
OUTDIR="A:/work/WinstonLab/Natalia_Reim/Paper_Thesis_Figs/"
cols = c("non-depleted"="#4477AA", "Spn1-depleted"="#BB5566", 
         "WT_DMSO"="#DDAA33", 
         "WT_IAA"="#dd8033",
         "noAID_DMSO"="#DDAA33", 
         "noAID_IAA"="#dd8033",
         "SPN1-AID_DMSO"="#4477AA","SPN1-AID_IAA"="#BB5566")
spiked.chr <- c("chrXVII","chrXVIII","chrXIX")
chr.flt <- c("chrI","chrII","chrIII","chrIV","chrIX","chrV","chrVI","chrVII","chrVIII","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")

## generate counts
#### Sense: HTSeq
#### Antisense: Custom script
if(T){
   setwd("A:/work/WinstonLab/Natalia_Reim/October2016/HTseq")
  files <- c("count_SPN1_DMSO_1.txt", 
             "count_SPN1_DMSO_2.txt",
             "count_SPN1_IAA_1.txt", "count_SPN1_IAA_2.txt", 
             "count_TIR4_DMSO_1.txt", 
             "count_TIR4_DMSO_2.txt", "count_TIR4_IAA_1.txt", 
             "count_TIR4_IAA_2.txt",  
             "count_WT_DMSO_1.txt", "count_WT_DMSO_2.txt", 
             "count_WT_IAA_1.txt", "count_WT_IAA_2.txt")
  htseq <- c()
  for (f in files) {
    d <- read.delim(f,header=F)
    d <- cbind(d, file=paste(f))
    htseq <- rbind(htseq,d)
    rm(d)
  }
  htseq$file <- sub(".txt","",htseq$file)
  htseq$file <- sub("count_","",htseq$file)
  names(htseq) <- c("tracking_id","counts","sample")
  htseq <- reshape2::dcast(htseq,tracking_id~sample,value.var = "counts")
  htseq <- htseq[-c(1:5),]
  write.table(htseq,file="Counttable.txt",quote = F,sep = "\t",row.names = F)
  
  files <- c("count_SPN1_DMSO_1.neg.txt", 
             "count_SPN1_DMSO_2.neg.txt",
             "count_SPN1_IAA_1.neg.txt", "count_SPN1_IAA_2.neg.txt", 
             "count_TIR4_DMSO_1.neg.txt", 
             "count_TIR4_DMSO_2.neg.txt", "count_TIR4_IAA_1.neg.txt", 
             "count_TIR4_IAA_2.neg.txt",  
             "count_WT_DMSO_1.neg.txt", "count_WT_DMSO_2.neg.txt", 
             "count_WT_IAA_1.neg.txt", "count_WT_IAA_2.neg.txt")
  htseq <- c()
  for (f in files) {
    d <- read.delim(f,header=F)
    d <- cbind(d, file=paste(f))
    htseq <- rbind(htseq,d)
    rm(d)
  }
  htseq$file <- sub(".txt","",htseq$file)
  htseq$file <- sub("count_","",htseq$file)
  names(htseq) <- c("tracking_id","counts","sample")
  htseq <- reshape2::dcast(htseq,tracking_id~sample,value.var = "counts")
  htseq <- htseq[-c(1:5),]
  write.table(htseq,file="ASCounttable.txt",quote = F,sep = "\t",row.names = F)
  
  
  
  load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenesFull.RData")
  htseq <- read.delim("A:/work/WinstonLab/Natalia_Reim/October2016/HTseq/Counttable.txt",header=T)
  htseq <- merge(htseq,unique(genes[,c("chr","tracking_id")]),by="tracking_id",all.x=T)
  rownames(htseq) <- htseq$tracking_id
  htseq$tracking_id <- NULL
  s <- htseq 
  names(s)[1:12] <- c("SPN1-AID_DMSO-1", "SPN1-AID_DMSO-2", "SPN1-AID_IAA-1", "SPN1-AID_IAA-2", "noAID_DMSO-1", 
                      "noAID_DMSO-2", "noAID_IAA-1", "noAID_IAA-2", "WT_DMSO-1", "WT_DMSO-2", 
                      "WT_IAA-1", "WT_IAA-2")
  htseq <- read.delim("A:/work/WinstonLab/Natalia_Reim/October2016/HTseq/ASCounttable.txt",header=T)
  htseq <- merge(htseq,unique(genes[,c("chr","tracking_id")]),by="tracking_id",all.x=T)
  rownames(htseq) <- htseq$tracking_id
  htseq$tracking_id <- NULL
  as <- htseq 
  names(as)[1:12] <- c("SPN1-AID_DMSO-1", "SPN1-AID_DMSO-2", "SPN1-AID_IAA-1", "SPN1-AID_IAA-2", "noAID_DMSO-1", 
                       "noAID_DMSO-2", "noAID_IAA-1", "noAID_IAA-2", "WT_DMSO-1", "WT_DMSO-2", 
                       "WT_IAA-1", "WT_IAA-2")
  
  
  check_spikeinValidity_sense <- function(s,return.norm=T){
    df1 <- s
    df1 <- df1[complete.cases(df1),]
    spiked.chr <- c("chrXVII","chrXVIII","chrXIX","chrMsp","chrM","Sp_chrAB325691")
    df1 <- df1[,c(13,1:12)]
    df2 <-df1
    df2[,2:13] <- apply(df2[,2:13],2,function(x){round(x*1000000/sum(x),2)})
    df1.c <- subset(df1, df1$chr%in%spiked.chr)
    df2.c <- subset(df2, df2$chr%in%spiked.chr)
    df1.e <- subset(df1, !df1$chr%in%spiked.chr)
    df2.e <- subset(df2, !df2$chr%in%spiked.chr)
    rm(df1,df2)
    
    # Normalization factor estimate from spikein using DESeq2
    library(DESeq2)
    df1.c <-as.matrix( round (df1.c[,2:13] ))
    type <- colnames(df1.c)
    condition <- sub("_[012]","",type)
    colData <- data.frame(condition,type)
    row.names(colData) <- colData$type
    dds <- DESeqDataSetFromMatrix(countData = df1.c,colData = colData,design = ~ condition)
    dds <- DESeq(dds)
    dd <- sizeFactors(dds)
    rm(dds)
    
    if(return.norm==T){
      return(dd)
      
    }else{
      # Normalization factor estimate from spikein using median expression value across SP genes
      m <- t(apply(df2.c[,2:13],2,median))
      n <- matrix(nrow=12,ncol=12)
      for (i in 1:12) {n[i,] <- round(m/m[i],2)}
      medsizeFactor <- as.vector(t(rowMeans(t(t(n)/n[1,]))))
      names(medsizeFactor) <- names(s)[1:12]
      rm(m,n)
      
      pdf("A:/work/WinstonLab/Natalia_Reim/Analysis_Feb18/SpikeinAssessment_HTSeqCounts.pdf",width = 11,height = 15)
      par(mfrow=c(4,2))
      par(mar=c(3,3,2,2))
      FP4.rle.boxplot(df1.c)         ## RLE distribution for spikein control
      mtext("Raw Sp tags",3,line=1)
      FP4.rle.boxplot(df1.e[,2:13])         ## RLE distribution for spikein control
      mtext("Raw Sc tags",3,line=1)
      
      FP4.rle.boxplot(df2.c[,2:13])         ## RLE distribution for spikein control
      mtext("Size-normalized Sp tags",3,line=1)
      FP4.rle.boxplot(df2.e[,2:13])         ## RLE distribution for spikein control
      mtext("Size-normalized Sc tags",3,line=1)
      
      FP4.rle.boxplot(t(t(df1.c)/(dd)))         ## RLE distribution for spikein control
      mtext("DESeq2 sizeFactor normalized raw Sp counts",3,line = 1)
      FP4.rle.boxplot(t(t(df1.e[,2:13])/(dd)))         ## RLE distribution for spikein control
      mtext("DESeq2 sizeFactor normalized raw Sc counts",3,line = 1)
      
      FP4.rle.boxplot(t(medsizeFactor*t(df2.c[,2:13])))         ## RLE distribution for spikein control
      mtext("Median sizeFactor and library size normalized Sc counts",3,line = 1)
      FP4.rle.boxplot(t(medsizeFactor*t(df2.e[,2:13])))         ## RLE distribution for spikein control
      mtext("Median sizeFactor and library size normalized Sp counts",3,line = 1)
      
      dev.off()
      rm(df1.c,df1.e,df2.c,df2.e,htseq,colData,type,condition)
      
    }
  }
  dd <- check_spikeinValidity_sense(s,return.norm = T)
  
  load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenesFull.RData")
  #genes <- genes[genes$species=="Scer",]
  genes <- genes[complete.cases(genes),]
  genes <- as(genes,"GRanges")
  length(genes[genes$species=="Scer" & seqnames(genes)!="chrM"])
  g <- disjoin(genes,ignore.strand=T)
  d <- as.data.frame(findOverlaps(g,genes,select="all",ignore.strand=T))
  d <- as.data.frame(table(d$queryHits))
  d <- d[d$Freq==1,]
  g <- g[as.numeric(d$Var1)]
  d <- findOverlaps(g,genes,select="all",ignore.strand=T)
  values(g) <- cbind(values(g),data.frame(values(genes[subjectHits(d)])))
  strand(g) <- strand(genes[subjectHits(d)])
  g$cov <- width(g)*100/width(genes[subjectHits(d)])
  #g <- g[g$cov>50]
  g[g$species=="Scer" & seqnames(g)!="chrM"]
  g[g$species=="Scer"]
  
  setwd("A:/work/WinstonLab/Natalia_Reim/October2016/gr")
  files <- Sys.glob("*.gr")
  as <- c()
  for (f in files) {
    cat(f,"\n")
    load(f)
    as <- cbind(as,countOverlaps(g,gr,minoverlap = 1,ignore.strand=F))# strand info in compl format in rnaseq bam
    rm(gr)
  }
  as <- as.data.frame(as)
  colnames(as) <- gsub(".gr","",files)
  g <- as.data.frame(g)
  as$tracking_id <- g$tracking_id
  as$chr <- g$seqnames
  as <- data.table(as)
  as <- as[,lapply(.SD,sum),by=list(tracking_id,chr)]
  as <- as.data.frame(as)
  rownames(as) <- as.character(as$tracking_id)
  as$tracking_id <- NULL 
  as <- as[,c(2:13,1)]
  
  ad <- check_spikeinValidity_sense(as,return.norm = T)
  
  load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenesFull.RData")
  names(as)[1:12] <- names(s)[1:12] <- names(ad) <- names(dd) <- c("SPN1-AID_DMSO-1", "SPN1-AID_DMSO-2", "SPN1-AID_IAA-1", "SPN1-AID_IAA-2", "noAID_DMSO-1", 
                                                                   "noAID_DMSO-2", "noAID_IAA-1", "noAID_IAA-2", "WT_DMSO-1", "WT_DMSO-2", 
                                                                   "WT_IAA-1", "WT_IAA-2")
  save(s,as,ad,dd,genes,file="A:/work/WinstonLab/Natalia_Reim/Analysis_Feb18/Sense_Antisese_ReadCounts_revised_01.RData")
}
## AS genes : non-overlapping, not within 100bp of each other on oppo strand
if(T){
  #rm(list=ls())
  load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenesFull.RData")
  #genes <- genes[genes$species=="Scer",]
  genes <- genes[genes$chr%in%chr.flt,]
  genes <- genes[!is.na(genes$verification),]
  
  genes <- as(genes,"GRanges")
  length(genes[genes$species=="Scer" & seqnames(genes)!="chrM"])
  g <- disjoin(genes,ignore.strand=T)
  d <- as.data.frame(findOverlaps(g,genes,select="all",ignore.strand=T))
  d <- as.data.frame(table(d$queryHits))
  d <- d[d$Freq==1,]
  g <- g[as.numeric(d$Var1)]
  d <- findOverlaps(g,genes,select="all",ignore.strand=T)
  values(g) <- cbind(values(g),data.frame(values(genes[subjectHits(d)])))
  strand(g) <- strand(genes[subjectHits(d)])
  g$cov <- width(g)*100/width(genes[subjectHits(d)])
  #g <- g[g$cov>50]
  
  
  ## AS 
  gas <- genes
  gas <- resize(gas,width(gas)+101,fix = "center")
  
  ol <- as.data.frame(findOverlaps(gas,ignore.strand=T,select="all"))
  ol <- ol[ol$queryHits!=ol$subjectHits,]
  ol <- ol[!duplicated(apply(ol,1,function(x) paste(sort(x),collapse=''))),]
  ol$chr1 <- as.character(seqnames(gas[ol$queryHits]))
  ol$tss1 <- as.numeric(gas[ol$queryHits]$TSS)
  ol$tts1 <- as.numeric(gas[ol$queryHits]$TTS)
  ol$str1 <- as.character(strand(gas[ol$queryHits]))
  ol$gene1 <- as.character(gas[ol$queryHits]$tracking_id)
  ol$chr2 <- as.character(seqnames(gas[ol$subjectHits]))
  ol$tss2 <- as.numeric(gas[ol$subjectHits]$TSS)
  ol$tts2 <- as.numeric(gas[ol$subjectHits]$TTS)
  ol$str2 <- as.character(strand(gas[ol$subjectHits]))
  ol$gene2 <- as.character(gas[ol$subjectHits]$tracking_id)
  
 
  ## completely enclosed genes (first 2), convergent genes within 100bp
  x1 <- ol[(ol$tss1 < ol$tss2 & ol$tts1 > ol$tts2)|
            (ol$tss2 < ol$tss1 & ol$tts2 > ol$tts1),]
  ol <- ol[!ol$queryHits%in%x1$queryHits,]
  x2 <- ol[(ol$tss2 < ol$tss1 & ol$tss2 > ol$tts1 & ol$str1!=ol$str2)|
           (ol$tss1 < ol$tss2 & ol$tss1 > ol$tts2 & ol$str1!=ol$str2),]
  ol <- ol[!ol$queryHits%in%x2$queryHits,]
  x3 <- ol[ol$str1!=ol$str2 & abs(ol$tss1-ol$tss2)>100, ]  
  ol <- ol[!ol$queryHits%in%x3$queryHits,]
  
  x <- rbind(x1,x2,x3)
  x <- unique(c(x$gene1,x$gene2))
  gas <- g[!g$tracking_id%in%x]
  gas <- as.character(gas$tracking_id)
  rm(g,x,x1,x2,ol,d,genes)
  
}

## Use my own counts in the antisense orientation
if(F){
  ## Figures produced using geomtric means of the signal across replicates
  
  load("A:/work/WinstonLab/Natalia_Reim/Analysis_Feb18/Sense_Antisese_ReadCounts_revised_01.RData")
  ## Sense_Antisese_ReadCounts_revised_01.RData = Sense counts using HTSeq and AScounts using my custom script
  ## AS counts along the gene body -including introns
  load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenesFull.RData")
  
  if(F){
    htseq <- read.delim("A:/work/WinstonLab/Natalia_Reim/October2016/HTseq/ASCounttable.txt",header=T)
    htseq <- merge(htseq,unique(genes[,c("chr","tracking_id")]),by="tracking_id",all.x=T)
    rownames(htseq) <- htseq$tracking_id
    htseq$tracking_id <- NULL
    as <- htseq 
    names(as)[1:12] <- c("SPN1-AID_DMSO-1", "SPN1-AID_DMSO-2", "SPN1-AID_IAA-1", "SPN1-AID_IAA-2", "noAID_DMSO-1", 
                         "noAID_DMSO-2", "noAID_IAA-1", "noAID_IAA-2", "WT_DMSO-1", "WT_DMSO-2", 
                         "WT_IAA-1", "WT_IAA-2")
    rm(htseq)
    #as <- as[rownames(as)%in%gas,]
    #s <- s[rownames(s)%in%gas,]
  }
  
  
  s <- subset(s,rownames(s)%in%genes[genes$verification=="Verified" & genes$species=="Scer",]$tracking_id)
  as <- subset(as,rownames(as)%in%genes[genes$verification=="Verified" & genes$species=="Scer",]$tracking_id)
  s <- s[s$chr%in%chr.flt,]
  as <- as[as$chr%in%chr.flt,]
  
  if(F){
    file=paste0(OUTDIR,"Natalia_proportion of AS to Sense.pdf")
    pdf(file,width = 7,height = 7)
    par(mar=c(4,14,5,2))
    barplot(apply(as[,1:12],2,sum)*100/apply(s[,1:12],2,sum),ylab="",main = "% antisense relative to sense",cex.lab=1.5,cex.axis=1.5,cex=1.5,cex.main=2,las=2,horiz=T)
    dev.off()
  }
  
  ## spikein normalize using size factors derived from Spombe
  s[,1:12] <- t(t(s[,1:12])/(dd[1:12]))
  as[,1:12] <- t(t(as[,1:12])/(dd[1:12]))
  s <- subset(s,rownames(s)%in%rownames(as))
  
  # genometric mean
  if(T){
    s1 <- c()
    as1 <- c()
    for( i in seq(1,11,2)){
      s1 <- cbind(s1, apply(s[,c(i,i+1)],1, function(x) sqrt(prod(x)) ))
      as1 <- cbind(as1, apply(as[,c(i,i+1)],1, function(x) sqrt(prod(x)) ))
    }
    s1 <- as.data.frame(s1)
    as1 <- as.data.frame(as1)
    names(s1) <- names(as1) <- gsub("-1","",names(s)[seq(1,11,2)])
  }
  
  if(T){
    zz <- as1
    zz$tracking_id <- rownames(zz)
    zz <- zz[zz$tracking_id%in%gas,]
    zz <- melt(zz)
    zz$value <- log2(zz$value+1)
    head(zz)
    zz$variable <- factor(zz$variable,levels=c("WT_DMSO", "WT_IAA","noAID_DMSO", "noAID_IAA","SPN1-AID_DMSO", "SPN1-AID_IAA"))
    p <- ggplot(zz, aes(x=variable,y=value,fill=variable))
    p <- p + geom_violin(color="black")
    p <- p + xlab("") + ylab("log2 (Antisense expression)")
    p <- p + geom_hline(yintercept = h,col="red",lty=2)
    p <- p + scale_fill_manual(values = cols)
    p <- p + geom_boxplot(width=0.1,fill='white', color="black",outlier.shape = NA)
    p <- p + theme(panel.background = element_rect(fill = "white", colour = "black")) 
    p <- p + theme(panel.grid = element_blank())
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=14,colour = "black"),
                   axis.text.y = element_text(size=14,colour = "black"),
                   axis.title =  element_text(size=15,color="black"))
    p <- p + theme(strip.background = element_blank(),
                   strip.text = element_text(size=16,colour = 'black',hjust = 0.5,face = "bold"))
    p <- p + theme(legend.position = "right",
                   legend.title = element_text(size=13, face="bold", color="black"),
                   legend.text = element_text(size=14, face="plain"),
                   legend.margin = margin(0,0,0,0),
                   legend.box.margin = margin(0,0,0,0))
    p <- p + theme(legend.position = "none")
    
    file=paste0(OUTDIR,"PF_ASExpression_log2expression.pdf")
    pdf(file,width = 6,height = 5)
    p
    dev.off()
    
  }
  
  ## average
  if(F){
    s1 <- round((s[,seq(1,11,by=2)] + s[,seq(2,12,by=2)])/2,2)
    as1 <- round((as[,seq(1,11,by=2)] + as[,seq(2,12,by=2)])/2,2)
    s1 <- as.data.frame(s1)
    as1 <- as.data.frame(as1)
    names(s1) <- names(as1) <- gsub("-1","",names(s1))
  }
  
  df.box <- log2(as1+1)-log2(s1+1)  ## for plotting boxes
  
  s1$tracking_id <- rownames(s1)
  as1$tracking_id <- rownames(as1)
  s1 <- s1[,c(1:2,3:4,5:7)]
  as1 <- as1[,c(1:2,3:4,5:7)]
  
  df.scatter <- as1  ## for plotting scatter
  plot(log2(as1$`SPN1-AID_DMSO`+1),log2(as1$`SPN1-AID_IAA`+1))
  plot(log2(as1$`SPN1-AID_DMSO`+1),log2(as1$WT_IAA+1))
  
  s1 <- reshape::melt.data.frame(s1,measure.vars = names(s1[1:6]))
  as1 <- reshape::melt.data.frame(as1,measure.vars = names(as1[1:6]))
  s1$ori <- "Sense"
  as1$ori <- "Antisense"
  df1 <- rbind(s1,as1)
  rm(s1,as1)
  
  df1 <- merge(df1,genes[,c("tracking_id","width")],by="tracking_id",all.x=T)
  df1$width <- df1$width/1000
  df1$value <- round(df1$value/df1$width,2)
  df1$value <- df1$value +0.1
  df1 <- df1[!df1$variable%in%c("noAID_IAA","noAID_DMSO"),]
  log2(range(df1$value))
  df1$class <- cut(log2(df1$value),breaks=-4:17)
  df1$variable <- gsub("-1","",df1$variable)
  #df1$variable <- gsub("SPN1-AID_DMSO","non-depleted",df1$variable)
  #df1$variable <- gsub("SPN1-AID_IAA","Spn1-depleted",df1$variable)
  names(df1)[2] <- "condition"
  
  p <- ggplot(df1,aes(class,fill=condition)) 
  p <- p + geom_bar(position="dodge") + geom_hline(yintercept =seq(0,2000,500),col="gray")
  p <- p + facet_wrap(~ori,nrow=2,scales = "free") + coord_cartesian() 
  p <- p + xlab("") + ylab("number of genes")
  p <- p + scale_fill_manual(values = cols)
  p <- p + theme(panel.background = element_rect(fill = "white", colour = "black")) 
  p <- p + theme(panel.grid = element_blank())
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1,size=14,colour = "black"),
                 axis.text.y = element_text(size=14,colour = "black"),
                 axis.title =  element_text(size=15,color="black")
                 )
  p <- p + theme(strip.background = element_blank(),
                 strip.text = element_text(size=16,colour = 'black',hjust = 0.5,face = "bold"))
  p <- p + theme(legend.position = "right",
                 legend.title = element_text(size=13, face="bold", color="black"),
                 legend.text = element_text(size=14, face="plain"),
                 legend.margin = margin(0,0,0,0),
                 legend.box.margin = margin(0,0,0,0))
  p <- p +scale_x_discrete(labels=c("(-4,-3]" = "0", "(-3,-2]" = "0.25", "(-2,-1]" = "0.5", "(-1,0]" = "1", "(0,1]" = "2", "(1,2]" = "4", "(2,3]" = "8", "(3,4]" = "16", "(4,5]" = "32", "(5,6]" = "64", "(6,7]" = "128", "(7,8]" = "256", "(8,9]" = "512", "(9,10]" = "1024", "(10,11]" = "2048", "(11,12]" = "4096", "(12,13]" = "8192", "(13,14]" = "16384", "(14,15]" = "32768", "(15,16]" = "65536", "(16,17]" = "131072","(17,18]" = "262144","(18,19]" = "524288","(19,20]" = "1048576"))
  p
  ######## BPKM plot
  file=paste0(OUTDIR,"Natalia_AntisenseBPKM plot_postspike_avgSignal.pdf")
  pdf(file,width = 8,height = 7)
  p
  dev.off()
  
  
  ####### Number of genes with 
  cf <- reshape2::dcast(df1,tracking_id+condition~ori,value.var = "value")
  cf$Antisense <- ifelse(is.na(cf$Antisense),0,cf$Antisense)
  cc <- cf[cf$condition=="WT_DMSO",]
  cc <- cc[order(cc$Sense,decreasing = F),]   
  cc$group <- rep(1:17,each=300)[1:length(cc$tracking_id)]  
  cf <- merge(cf,cc[,c("group","tracking_id")],by="tracking_id",all.x=T)  
  head(cf)
  names(cf)[2] <- "sample"
  
  z <- rbind(cbind(category="Sense, BPKM>1",as.data.frame(table(cf[cf$Sense>1,]$sample)*100/length(cc[,1]))),
             cbind(category="Antisense, BPKM>1",as.data.frame(table(cf[cf$Antisense>1,]$sample)*100/length(cc[,1]))),
             cbind(category="Sense, BPKM>3",as.data.frame(table(cf[cf$Sense>3,]$sample)*100/length(cc[,1]))),
             cbind(category="Antisense, BPKM>3",as.data.frame(table(cf[cf$Antisense>3,]$sample)*100/length(cc[,1]))),
             cbind(category="Sense, BPKM>10",as.data.frame(table(cf[cf$Sense>10,]$sample)*100/length(cc[,1]))),
             cbind(category="Antisense, BPKM>10",as.data.frame(table(cf[cf$Antisense>10,]$sample)*100/length(cc[,1]))),
             cbind(category="Sense, BPKM>50",as.data.frame(table(cf[cf$Sense>50,]$sample)*100/length(cc[,1]))),
             cbind(category="Antisense, BPKM>50",as.data.frame(table(cf[cf$Antisense>50,]$sample)*100/length(cc[,1]))))
  z$Freq <- round(z$Freq,2)
  z$Var1 <- gsub("-1","",z$Var1)
  z <- reshape2::dcast(z,category~Var1,value.var = "Freq")
  z
  
  library(gridExtra)
  library(grid)
  file=paste0(OUTDIR,"Natalia_PercentGenes_withDetectableExpression.pdf")
  pdf(file,width = 15,height = 5)
  grid.table(z)
  dev.off()
  
  
  ###### Log2 (AS/S) boxplots
  pl <- df.box[,c(1:2,5:6)]
  boxplot(pl)
  pl <- melt(pl)
  #pl$variable <- gsub("SPN1-AID_DMSO","non-depleted",pl$variable)
  #pl$variable <- gsub("SPN1-AID_IAA","Spn1-depleted",pl$variable)
  names(pl)[1] <- "condition"
  
  h <- median(pl[pl$condition=="SPN1-AID_DMSO",]$value)
  
  p <- ggplot(pl, aes(x=condition,y=value,fill=condition))
  p <- p + geom_violin(color="black")
  p <- p + xlab("") + ylab("log2 (Antisense/Sense)")
  p <- p + geom_hline(yintercept = h,col="red",lty=2)
  p <- p + scale_fill_manual(values = cols)
  p <- p + geom_boxplot(width=0.25,fill='white', color="black",outlier.shape = NA)
  p <- p + theme(panel.background = element_rect(fill = "white", colour = "black")) 
  p <- p + theme(panel.grid = element_blank())
  p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1,size=14,colour = "black"),
                 axis.text.y = element_text(size=14,colour = "black"),
                 axis.title =  element_text(size=15,color="black"))
  p <- p + theme(strip.background = element_blank(),
                 strip.text = element_text(size=16,colour = 'black',hjust = 0.5,face = "bold"))
  p <- p + theme(legend.position = "right",
                 legend.title = element_text(size=13, face="bold", color="black"),
                 legend.text = element_text(size=14, face="plain"),
                 legend.margin = margin(0,0,0,0),
                 legend.box.margin = margin(0,0,0,0))
  p <- p + theme(legend.position = "none")
  p
  
  file=paste0(OUTDIR,"PF_ASExpression_toSense expression.pdf")
  pdf(file,width = 6,height = 5)
  p
  dev.off()
  
  
  ######### Scatter plots of AS sginals
  pl <- df.scatter[,c(1,2,5,6)]
  pl$tracking_id <- rownames(pl)
  pl <- merge(pl,genes[,c("tracking_id","width")],all.x=T,by="tracking_id")
  pl$width <- pl$width/1000
  pl[,2:5] <- pl[,2:5]/pl$width
  pl$tracking_id <- NULL
  pl$width <- NULL
  names(pl) <- c("a","b","c","d")
  pscount <- 1
  
  g1 <-ggplot(pl, aes(x=log2(a+pscount), y=log2(b+pscount))) +
    stat_bin_hex(geom="point", aes(color=log2(..count..)), binwidth=c(.01,.01), size=2, shape=16, stroke=0)+
      scale_color_viridis(option="inferno",direction = 1)+
    theme(axis.text.x = element_text(size=14,colour = "black"),
          axis.text.y = element_text(size=14,colour = "black"),
          axis.title =  element_text(size=15,color="black"))+
    theme(legend.position = "none")+
    theme(panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid = element_blank())+
    geom_abline(slope = 1,intercept = 0,col="black")+
    xlab("non-depleted")+ylab("Spn1-depleted")

  g2 <- ggplot(pl, aes(x=log2(c+pscount), y=log2(a+pscount))) +
    stat_bin_hex(geom="point", aes(color=log2(..count..)), binwidth=c(.01,.01), size=2, shape=16, stroke=0)+
    scale_color_viridis(option="inferno",direction = 1)+
    theme(axis.text.x = element_text(size=14,colour = "black"),
          axis.text.y = element_text(size=14,colour = "black"),
          axis.title =  element_text(size=15,color="black"))+
    theme(legend.position = "none")+
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid = element_blank())+
    geom_abline(slope = 1,intercept = 0,col="black")+
    xlab("WT_DMSO")+ylab("non-depleted")
  
  file=paste0(OUTDIR,"PF_AS_NormalizedCountComparison.pdf")
  pdf(file,width = 8,height = 4)
  print(grid.arrange(g1,g2,ncol=2,top= textGrob("FPKM counts",gp=gpar(fontsize=20,font=2))))
  dev.off()
}

# which genes to choose
## finding candidates:Natlia
if(F){
  load("A:/work/WinstonLab/Natalia_Reim/Analysis_Feb18/Sense_Antisese_ReadCounts_revised_01.RData")
  load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenesFull.RData")
  load("A:/work/yeast_Annotations/refs.RData")
  
  if(F){
    htseq <- read.delim("A:/work/WinstonLab/Natalia_Reim/October2016/HTseq/ASCounttable.txt",header=T)
    htseq <- merge(htseq,unique(genes[,c("chr","tracking_id")]),by="tracking_id",all.x=T)
    rownames(htseq) <- htseq$tracking_id
    htseq$tracking_id <- NULL
    as <- htseq 
    names(as)[1:12] <- c("SPN1-AID_DMSO-1", "SPN1-AID_DMSO-2", "SPN1-AID_IAA-1", "SPN1-AID_IAA-2", "noAID_DMSO-1", 
                         "noAID_DMSO-2", "noAID_IAA-1", "noAID_IAA-2", "WT_DMSO-1", "WT_DMSO-2", 
                         "WT_IAA-1", "WT_IAA-2")
    rm(htseq)
  }
  
  
  s <- subset(s,rownames(s)%in%genes[genes$verification=="Verified",]$tracking_id)
  as <- subset(as,rownames(as)%in%genes[genes$verification=="Verified",]$tracking_id)
  as <- subset(as,as$chr%in%chr.flt)
  s <- subset(s,s$chr%in%chr.flt)
  
  ## need minimum 3 reads
  #x <- as[,1:12]
  #x[x<3] <- 0
  #as[,1:12] <- x
  ## spike in normalization
  as[,1:12] <- t(t(as[,1:12])/(dd[1:12]))
  #as[,1:12] <- apply(as[,1:12],2,function(x) x*1e6/sum(x))
  
  ast <- as
  ast$tracking_id <- rownames(as)
  ast <- merge(ast,genes[,c("tracking_id","width")],by="tracking_id",all.x=T)
  ast$width <- ast$width/1000
  ast[,2:13] <- ast[,2:13]/ast$width ## RPKM signal
  as1 <- c()
  for( i in seq(2,12,2)){
    as1 <- cbind(as1, apply(ast[,c(i,i+1)],1, function(x) sqrt(prod(x)) ))
  }
  as1 <- as.data.frame(as1)
  names(as1) <- gsub("-1","",names(ast)[seq(2,12,2)])
  rownames(as1) <- ast$tracking_id
  kz <- as1
  names(kz) <- paste0("GMean_",names(kz)) 
  kz <- cbind(kz,ast[,2:13])
  ast <- as1
  rm(as1)
  
  kz <- subset(kz,rownames(kz)%in%gas)
  kz$tracking_id <- rownames(kz)
  kz <- kz[,c(19,1:18)]
  kz <- merge(kz,genes[,c("tracking_id","gene")],by="tracking_id",all.x=T)
  #write.table(kz, "A:/work/WinstonLab/Natalia_Reim/ThesisFigures_March2019/AS.table.NataliaReqest.txt",quote = F,sep = "\t",row.names = F)
  #write.table(kz, "A:/work/WinstonLab/Natalia_Reim/ThesisFigures_March2019/AS.table.NataliaReqest_CustomAS.txt",quote = F,sep = "\t",row.names = F)
  
  
  myfun <- function(ast,genes,cutoff=5){
    #ast[ast<cutoff] <- 0 ## minimum RPKM of 3
    ast <- ast[,c(2,1,4,3,6,5)]
    gst <- ast
    gst <- gst[apply(gst[,1:2],1,max)>=cutoff,]
    gst <- gst*100+1
    gst$`SPN1-AID_IAA/SPN1-AID_DMSO` <- log2(gst$`SPN1-AID_IAA`/gst$`SPN1-AID_DMSO`)
    gst$`SPN1-AID_DMSO/WT_DMSO` <- log2(gst$`SPN1-AID_DMSO`/gst$WT_DMSO)
    gst$`noAID_IAA/WT_DMSO` <- log2(gst$noAID_IAA/gst$WT_DMSO)
    gst$`noAID_DMSO/WT_DMSO` <- log2(gst$noAID_DMSO/gst$WT_DMSO)
    gst$`WT_IAA/WT_DMSO` <- log2(gst$WT_IAA/gst$WT_DMSO)
    gst <- gst[,7:11]
    
    gst <- gst[abs(gst$`SPN1-AID_IAA/SPN1-AID_DMSO`)>1,]
    gst$absm <- apply(gst[,2:5],1,function(x){
      rr <- max(abs(x))
      if(rr == abs(min(x))){
        return(rr*(-1))
      }else{
        return(rr)
      }
    })
    gstf <- gst[gst$absm>1 & gst$`SPN1-AID_IAA/SPN1-AID_DMSO`*gst$absm>1,]
    gst <- subset(gst,!rownames(gst)%in%rownames(gstf))
    
    
    #ast$control_max <- apply(ast[,2:6],1,max)
    #ps=1
    #ast$l2fc <- log2(ast$`SPN1-AID_IAA`*100+ps) -log2(ast$`SPN1-AID_DMSO`*100+ps)
    #ast$l2fc_control <- log2(ast$`SPN1-AID_IAA`*100+ps) -log2(ast$control_max*100+ps)
    #ast$tracking_id <- rownames(ast)
    #ast <- merge(ast,refs[,c("tracking_id","gene")],by="tracking_id",all.x=T)
    #ast1 <- ast[abs(ast$l2fc)>1 & abs(ast$l2fc_control)>1,]
    #ast1 <- ast1[ast1$tracking_id%in%as.character(gas),]  
    #sum(ast1$l2fc>0)
    #sum(ast1$l2fc<0)
    #names(ast1)[9:10] <- c("Log2FC","Log2FC (SPN1_IAA/ControlMax)") 
    #ast1 <- merge(ast1,genes,by="tracking_id",all.x=T) 
    gst$tracking_id <- rownames(gst)
    gst <- gst[gst$tracking_id%in%as.character(gas),]  
    gst <- merge(gst,genes,by="tracking_id",all.x=T)
    names(gst)[2] <- "Log2FC"
    gst
  }
  ast_3rpkm <- myfun(ast,genes,cutoff = 3)
  ast_5rpkm <- myfun(ast,genes,cutoff = 5)
  ast_10rpkm <- myfun(ast,genes,cutoff = 10)
  
  sum(ast_3rpkm$Log2FC>0)
  sum(ast_3rpkm$Log2FC<0)
  sum(ast_5rpkm$Log2FC>0)
  sum(ast_5rpkm$Log2FC<0)
  sum(ast_10rpkm$Log2FC>0)
  sum(ast_10rpkm$Log2FC<0)
  
  sum(ast_3rpkm$Log2FC<0)/nrow(ast_3rpkm)
  sum(ast_5rpkm$Log2FC<0)/nrow(ast_5rpkm)
  sum(ast_10rpkm$Log2FC<0)/nrow(ast_10rpkm)
  
  
  library(xlsx)
  file=paste0(OUTDIR,"Log2FC_AS_transcription_April19.xlsx")
  write.xlsx(ast_3rpkm,file=paste0(file),sheetName = "3RPKM cutoff",row.names = F,append = F)
  write.xlsx(ast_5rpkm,file=paste0(file),sheetName = "5RPKM cutoff",row.names = F,append = T)
  write.xlsx(ast_10rpkm,file=paste0(file),sheetName = "10RPKM cutoff",row.names = F,append = T)
  #write.table(ast1,file,sep="\t",quote = F,row.names = F)
  
}




