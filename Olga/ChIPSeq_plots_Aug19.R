rm(list=ls())
require(GenomicRanges)
require(ggplot2)
require(factoextra)
require(Hmisc)
library(viridis)
library(EnrichedHeatmap)
library(circlize)
require(ape)
require(ggrepel)
library(ggdendro)
library(psych)

chr.flt <- c("chrI","chrII","chrIII","chrIV","chrIX","chrV","chrVI","chrVII","chrVIII","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
spiked.chr <- c("chrXVII","chrXVIII","chrXIX")
col=c( ## merged replicate names
      "WT_Myc"="black","WT_8WG16"="black", "WT_Flag"="black", 
      "pob3_Myc"="green3", "pob3_8WG16"="green3", "pob3_Flag"="green3", 
      "spt6-YW_Myc"="blue","spt6-YW_Flag"="blue", "spt6-YW_8WG16"="blue", 
      "spt6-YW-pob3_Myc"="red","spt6-YW-pob3_8WG16"="red", "spt6-YW-pob3_Flag"="red",
      
      ## individual replicates
      "pob3_8WG16-1"="green3", "pob3_8WG16-2"="green1", "pob3_Flag-1"="green3", 
      "pob3_Flag-2"="green1", "pob3_Myc-1"="green3", "pob3_Myc-2"="green1", 
      "spt6-YW-pob3_8WG16-1"="red", "spt6-YW-pob3_8WG16-2"="pink", "spt6-YW-pob3_Flag-1"="red",
      "spt6-YW-pob3_Flag-2"="pink", "spt6-YW-pob3_Myc-1"="red", "spt6-YW-pob3_Myc-2"="pink", 
      "spt6-YW_8WG16-1"="blue", "spt6-YW_8WG16-2"="lightblue", "spt6-YW_Flag-1"="blue", 
      "spt6-YW_Flag-2"="lightblue", "spt6-YW_Myc-1"="blue", "spt6-YW_Myc-2"="lightblue", 
      "WT_8WG16-1"="black", "WT_8WG16-2"="gray", "WT_Flag-1"="black", 
      "WT_Flag-2"="gray", "WT_Myc-1"="black", "WT_Myc-2"="gray",
      
      "pob3-1"="green3", "pob3-2"="green1", 
      "double-1"="red", "double-2"="pink", 
      "spt6-YW-1"="blue", "spt6-YW-2"="lightblue",
      "WT-1"="black", "WT-2"="gray",
      
      "pob3"="green3",
      "double"="red",
      "spt6-YW"="blue", 
      "WT"="black",
      
      ## temp
      "37"="slateblue4","30"="slateblue1"
)
ggscat <- function (data, columns = 1:ncol(data), color = NULL, alpha = 1, corMethod = "pearson",size=14) {
  data <- upgrade_scatmat_data(data)
  data.choose <- data[columns]
  dn <- data.choose[sapply(data.choose, is.numeric)]
  if (ncol(dn) == 0) {
    stop("All of your variables are factors. Need numeric variables to make scatterplot matrix.")
  }
  if (ncol(dn) < 2) {
    stop("Not enough numeric variables to make a scatter plot matrix")
  }
  a <- uppertriangle(data, columns = columns, color = color, 
                     corMethod = corMethod)
  if (is.null(color)) {
    plot <- scatmat(data, columns = columns, alpha = alpha) + 
      geom_text(data = a, aes_string(label = "r"), colour = "black",size=size)
  }
  else {
    plot <- scatmat(data, columns = columns, color = color, 
                    alpha = alpha) + geom_text(data = a, aes_string(label = "r", 
                                                                    color = "colorcolumn"),size=size) + labs(color = color)
  }
  factor <- data.choose[sapply(data.choose, is.factor)]
  if (ncol(factor) == 0) {
    return(plot)
  }
  else {
    warning("Factor variables are omitted in plot")
    return(plot)
  }
}
environment(ggscat) <- asNamespace('GGally')


### Manual handlig ##
### Sample normalization factors for July 2019 ChIPSeq
if(F){
  setwd("A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr")
  OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr"
  load("proportion.spikein.RData")
  head(res)
  write.table(res,file="spikein_percents.txt",sep = "\t",quote = F,row.names = F)
  res <- read.delim("NormalizationFactors.txt",header = T)
  res <- res[c(1,3,2,4,
               5,7,6,8,
               9,11,10,12,
               13,15,14,16,
               17,19,18,20,
               21,23,22,24,
               25,27,26,28,
               29,31,30,32,
               33,35,34,36,
               37,39,38,40,
               41,43,42,44,
               45,47,46,48),]
  write.table(res,file="NormalizationFactors.txt",sep="\t",quote = F,row.names = F)
}
### Manual handlig ##
### Sample normalization factors for August 2019 ChIPSeq
if(F){
  setwd("A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr")
  OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr"
  load("proportion.spikein.RData")
  head(res)
  write.table(res,file="spikein_percents.txt",sep = "\t",quote = F,row.names = F)
  res <- read.delim("NormalizationFactors.txt",header = T)
  res <- res[c(1,3,2,4,
               5,7,6,8,
               9,11,10,12,
               13,15,14,16,
               17,19,18,20,
               21,23,22,24,
               25,27,26,28,
               29,31,30,32,
               33,35,34,36,
               37,39,38,40,
               41,43,42,44,
               45,47,46,48),]
  write.table(res,file="NormalizationFactors.txt",sep="\t",quote = F,row.names = F)
}

## %spikein plots
if(F){
  ProportionSpikeInFig <- function(data.path="proportion.spikein.RData",norm.path="NormalizationFactors.txt"){
    load(paste0(data.path))
    head(res)
    d <- read.delim(paste0(norm.path),header = T,stringsAsFactors = F)
    res <- res[match(d$sample,res$sample),]
    res$sample <- factor(res$sample,levels=as.character(d$sample))
    res$Factor <- d$Factor
    res$normF <- d$normF
    res$condition <- d$condition
    head(res)
    res$sample <- as.character(res$sample)
    rownames(res) <- 1:nrow(res)
    res$Factor <- as.character(res$Factor)
    res$condition <- as.character(res$condition)
    res$condition[33:48] <- rep(c(30,30,37,37),4)
    
    res1 <- res
    res$fraction <- res$heterosomes/(res$autosomes+res$heterosomes)
    res$sample <- factor(res$sample,levels=res$sample)
    
    p1 <-ggplot(res,aes(x=sample, y=fraction, fill=condition)) + geom_bar(stat="identity",position="dodge")+ 
      facet_wrap(~Factor,nrow=4,scales = "free_x")+
      xlab("") + ylab("Fraction of total reads") + theme_minimal()+
      scale_fill_manual(values = col) + ylim(0,0.2) +
      theme(axis.text.x = element_text(angle = 70,size=14,color="black",hjust = 1),
            axis.text.y = element_text(size=16,color="black",hjust = 1),
            axis.title = element_text(size=18,color="black"),
            legend.position = "none")+
      theme(strip.background = element_blank(),
            strip.text = element_text(size=20,color="black"))+
      geom_hline(yintercept = 0.1,col="red",size=1,lty=2)
    p1
    
  }
  
  
  ####### % spikein signal
  out.path=paste0(OUTDIR,"/Spikein proportions.pdf")
  pdf(out.path,width = 10,height = 16)
  ProportionSpikeInFig(data.path=paste0(OUTDIR,"/proportion.spikein.RData"),
                       norm.path=paste0(OUTDIR,"/NormalizationFactors.txt"))
  dev.off()
  
  out.path=paste0(OUTDIR,"/Spikein proportions.pdf")
  pdf(out.path,width = 10,height = 16)
  ProportionSpikeInFig(data.path=paste0(OUTDIR,"/proportion.spikein.RData"),
                       norm.path=paste0(OUTDIR,"/NormalizationFactors.txt"))
  dev.off()
  
  
  
  
}


## ChIP/Input for Spombe and Scer :: estimation of IPing efficiency
if(T){
  ProportionSpikeInFig <- function(data.path="A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/ReadCounts.txt"){
    d <- read.delim(paste0(data.path),header = T,stringsAsFactors = F)
    d$tmp <- d[,4]+d[,5]
    d[,4:5] <- round(d[,4:5]*1e7/d$tmp)
    d$tmp <- d[,6]+d[,7]
    d[,6:7] <- round(d[,6:7]*1e7/d$tmp)
    d$Spombe <- d[,5]/d[,7]
    d$Scer <- d[,4]/d[,6]
    d$tmp <- gsub("\\S*_","",d$ChIP)
    d$'Scer/Spom' <- d$Spombe/d$Scer
    
    
    par(mfrow=c(2,2))
    par(mar=c(4,15,4,4))
    for(f in c("Flag","8WG16")){
      for(tmp in c(30,37)){
        v <- d[d$tmp==tmp,]$`Scer/Spom`
        if(length(v)==0){
          next
        }
        names(v) <- d[d$tmp==tmp,]$ChIP
        v <- v[grep(f,names(v))]
        h<- median(v[grep("WT|wt",names(v))])
        col <- rep("gray",length(v))
        col[grep("WT|wt",names(v))] <- "red"
        barplot(v,las=2,horiz=T,main=paste0(f," @",tmp),cex.main=2,col = col)
        abline(v=c(v[grep("WT|wt",names(v))]),col="red",lty=2)
        rm(v,h,col)
      }
    }
    
  }
  
  pdf("A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/GlobalChangeAssess_SpomtoScer_July.pdf",width = 10,height = 8)
  ProportionSpikeInFig(data.path="A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/ReadCounts.txt")
  dev.off()
  
  pdf("A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/GlobalChangeAssess_SpomtoScer_July.pdf",width = 10,height = 8)
  ProportionSpikeInFig(data.path="A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/ReadCounts.txt")
  dev.off()
  
}


### sample correlation plots
## Peason Correlations for the paper
if(F){
  #OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr"
  OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr"
  
  Inp.file=paste0(OUTDIR,"/CountMatrix_SignalRegions_5p.RData")
  norm.file=paste0(OUTDIR,"/NormalizationFactors.txt")
  out.file1=paste0(OUTDIR,"Sample correlations _5p signal regions_final.pdf")
  out.file2=paste0(OUTDIR,"Sample correlations _5p signal regions_subsets.pdf")
  
  source("A:/work/scripts/Project_Winston/final/winstonLab_functions.R")
  
  load(paste0(Inp.file)) 
  m <- cf[,6:ncol(cf)]
  normF <- read.delim(paste0(norm.file),header=T)
  normF$sample <- gsub("-",".",normF$sample)
  normF$sample%in%names(m)
  
  m <- m[,normF$sample]
  normF$sample==names(m)
  m1 <- matrix(NA,nrow = nrow(m),ncol = ncol(m))
  for(i in 1:ncol(m)){
    m1[,i] <- m[,i]*normF$normF[i]
  }
  m1 <- as.data.frame(m1)
  names(m1) <- names(m)
  
  ## print correlations
  pdf(paste0(out.file1),width = 14,height = 14)
  {
    require(gplots)
    ## All samples (Raw, SpikeNorm)
    z1 <- round(cor(log2(m+1),use="pairwise.complete.obs",method = "pearson"),2)
    
    hmcols <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
    heatmap.2(z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.0,0.0),sepcolor = "black",
              main="Pearson correlation (pre-spikein)",key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
    hmcols <- palette(brewer.pal(n = 11, name = "RdYlBu"))
    heatmap.2(z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.0,0.0),sepcolor = "black",
              main="Pearson correlation (pre-spikein)",key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
    hmcols <- (colorRampPalette(c("white", "black"))(256))
    heatmap.2(z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.0,0.0),sepcolor = "black",
              main="Pearson correlation (pre-spikein)",key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
    ## post spikein
    z1 <- round(cor(log2(m1+1),use="pairwise.complete.obs",method = "pearson"),2)
    hmcols <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
    heatmap.2(z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.0,0.0),sepcolor = "black",
              main="Pearson correlation (post-spikein)",key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
    hmcols <- palette(brewer.pal(n = 11, name = "RdYlBu"))
    heatmap.2(z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.0,0.0),sepcolor = "black",
              main="Pearson correlation (post-spikein)",key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
    hmcols <- (colorRampPalette(c("white", "black"))(256))
    heatmap.2(z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.0,0.0),sepcolor = "black",
              main="Pearson correlation (post-spikein)",key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
    
  }
  dev.off()
  
  hmcols <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  pdf(paste0(out.file1),width = 12,height = 12)
  {
    n <- names(m) ## total 16 inputs + 48 ChIPs
    n <- n[grep("8WG16ChIP",n)]
    z1 <- round(cor(log2(m[ ,names(m)%in%n]+1)),2)
    heatmap.2(main = "ChIP-Seq correlations (pre-normalization)",z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.01,0.01),sepcolor = "black",colsep=0:ncol(z1),rowsep=0:nrow(z1),
              key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
    z1 <- round(cor(log2(m1[ ,names(m1)%in%n]+1)),2)
    heatmap.2(main = "ChIP-Seq correlations (post-normalization)",z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.01,0.01),sepcolor = "black",colsep=0:ncol(z1),rowsep=0:nrow(z1),
              key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
    
    ## All Histone samples (Raw, SpikeNorm)
    n <- names(m) ## total 16 inputs + 48 ChIPs
    n <- n[grep("MycChIP|FlagChIP",n)]
    z1 <- round(cor(log2(m[ ,names(m)%in%n]+1)),2)
    heatmap.2(main = "ChIP-Seq correlations (pre-normalization)",z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.01,0.01),sepcolor = "black",colsep=0:ncol(z1),rowsep=0:nrow(z1),
              key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
    z1 <- round(cor(log2(m1[ ,names(m1)%in%n]+1)),2)
    heatmap.2(main = "ChIP-Seq correlations (post-normalization)",z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.01,0.01),sepcolor = "black",colsep=0:ncol(z1),rowsep=0:nrow(z1),
              key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
    
    
    ## All Histone samples (Raw, SpikeNorm)
    n <- names(m) ## total 16 inputs + 48 ChIPs
    n <- n[-grep("Inp",n)]
    z1 <- round(cor(log2(m[ ,names(m)%in%n]+1)),2)
    heatmap.2(main = "ChIP-Seq correlations (pre-normalization)",z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.01,0.01),sepcolor = "black",colsep=0:ncol(z1),rowsep=0:nrow(z1),
              key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
    z1 <- round(cor(log2(m1[ ,names(m1)%in%n]+1)),2)
    heatmap.2(main = "ChIP-Seq correlations (post-normalization)",z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.01,0.01),sepcolor = "black",colsep=0:ncol(z1),rowsep=0:nrow(z1),
              key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
  }
  dev.off()
  
  v <- c()
  v1 <-c()
  for(i in seq(1,47,2)){
    cat(names(m)[i],"\t",names(m)[i+1],"\n")
    v <- c(v,cor(m[,i],m[,(i+1)],method = "pearson",use = "complete.obs"))
    v1 <- c(v1,cor(m1[,i],m1[,(i+1)],method = "pearson",use = "complete.obs"))
  }
  names(v) <- names(v1) <- names(m)[seq(1,47,2)]
  names(v)[25:48] <- names(v1)[25:48] <- names(m)[seq(1,47,2)]
  
}

### Correlation plots :: scatters using SERs at5% for Scer 
if(T){
  myplot <- function(m,titl){
    q <- ggscat(m, columns = 1:ncol(m),color = NULL,size=10)+
      stat_bin_hex(geom="point", aes(color=(..count..)), binwidth=c(.05,.05), size=2, shape=16, stroke=0)+
      scale_color_viridis(option="inferno",direction = 1)+
      theme(axis.text.x = element_text(size=18,colour = "black"),
            axis.text.y = element_text(size=18,colour = "black"),
            axis.title =  element_text(size=20,color="black"),
            strip.text = element_text(color="black",size=20))+
      theme(legend.position = "none")+
      theme(panel.background = element_rect(fill = "white", colour = "black"),
            panel.grid = element_blank())+
      geom_smooth(method='lm',formula=y~x,col="red")+
      xlab("factor occupancy,log2")+ylab("factor occupancy, log2")+
      ggtitle(paste0(titl))+
      theme(plot.title = element_text(colour ="black",size = 18,hjust = 0.5))
    q
  }
  
  OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr"
  Inp.file=paste0(OUTDIR,"/CountMatrix_SignalRegions_5p.RData")
  norm.file=paste0(OUTDIR,"/NormalizationFactors.txt")
  load(paste0(Inp.file)) 
  m <- cf[!cf$seqnames%in%spiked.chr & cf$width<500,6:ncol(cf)]
  normF <- read.delim(paste0(norm.file),header=T)
  normF$sample <- gsub("-",".",normF$sample)
  normF$sample%in%names(m)
  m <- m[,normF$sample]
  normF$sample==names(m)
  m1 <- matrix(NA,nrow = nrow(m),ncol = ncol(m))
  for(i in 1:ncol(m)){
    m1[,i] <- m[,i]*normF$normF[i]
  }
  m1 <- as.data.frame(m1)
  h <- names(m)
  h <- gsub("8WG16ChIP","Pol2",h)
  h <- gsub("MycChIP","Spt16",h)
  h <- gsub("FlagChIP","Spt6",h)
  h <- gsub("spt6.YW","spt6-YW",h)
  h <- gsub("_r1_30",",30°C,R1",h)
  h <- gsub("_r2_30",",30°C,R2",h)
  h <- gsub("_r1_37",",37°C,R1",h)
  h <- gsub("_r2_37",",37°C,R2",h)
  h <- gsub("spt6-YW.pob3","double",h)
  names(m) <- names(m1) <- h
  for(f in c("Pol2","Inp")){
    for(tmp in c(30,37)){
      n <- names(m) ## total 16 inputs + 48 ChIPs
      n <- n[grep(f,n)]
      n <- n[grep(tmp,n)]
      z1 <- log2(m[ ,names(m)%in%n]+1)
      fout = paste0(f,"_SampleCorrelations_preSpikein_",tmp,"_July.png")
      cat(fout,"\n")
      cat(dim(z1),"\n")
      png(file = paste0(fout),width = 2000,height = 2000)
       print(myplot(z1,titl = f))
      dev.off()
      z1 <- log2(m1[ ,names(m1)%in%n]+1)
      fout = paste0(f,"_SampleCorrelations_postSpikein_",tmp,".png")
      cat(fout,"\n\n")
      cat(dim(z1),"\n")
      png(file = paste0(fout),width = 2000,height = 2000)
      print(myplot(z1,titl = f))
      dev.off()
    }    
  }
  
  
  OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr"
  Inp.file=paste0(OUTDIR,"/CountMatrix_SignalRegions_5p.RData")
  norm.file=paste0(OUTDIR,"/NormalizationFactors.txt")
  load(paste0(Inp.file)) 
  m <- cf[!cf$seqnames%in%spiked.chr & cf$width<500,6:ncol(cf)]
  normF <- read.delim(paste0(norm.file),header=T)
  normF$sample <- gsub("-",".",normF$sample)
  normF$sample%in%names(m)
  m <- m[,normF$sample]
  normF$sample==names(m)
  m1 <- matrix(NA,nrow = nrow(m),ncol = ncol(m))
  for(i in 1:ncol(m)){
    m1[,i] <- m[,i]*normF$normF[i]
  }
  m1 <- as.data.frame(m1)
  h <- names(m)
  h <- gsub("8WG16ChIP","Pol2",h)
  h <- gsub("MycChIP","Spt16",h)
  h <- gsub("FlagChIP","Spt6",h)
  h <- gsub("spt6.YW","spt6-YW",h)
  h <- gsub("_r1_30",",30°C,R1",h)
  h <- gsub("_r2_30",",30°C,R2",h)
  h <- gsub("_r1_37",",37°C,R1",h)
  h <- gsub("_r2_37",",37°C,R2",h)
  h <- gsub("spt6-YW.pob3","double",h)
  names(m) <- names(m1) <- h
   for(f in c("Spt6","Pol2","Inp")){
    for(tmp in c(30,37)){
      n <- names(m) ## total 16 inputs + 48 ChIPs
      n <- n[grep(f,n)]
      n <- n[grep(tmp,n)]
      z1 <- log2(m[ ,names(m)%in%n]+1)
      fout = paste0(f,"_SampleCorrelations_preSpikein_",tmp,".png")
      cat(fout,"\n")
      cat(dim(z1),"\n")
      png(file = paste0(fout),width = 2000,height = 2000)
      print(myplot(z1,titl = f))
      dev.off()
      z1 <- log2(m1[ ,names(m1)%in%n]+1)
      fout = paste0(f,"_SampleCorrelations_postSpikein_",tmp,".png")
      cat(fout,"\n\n")
      cat(dim(z1),"\n")
      png(file = paste0(fout),width = 2000,height = 2000)
      print(myplot(z1,titl = f))
      dev.off()
    }    
  }
}

### Correlation plots :: scatters using genes for Spom 
if(T){
  myplot <- function(m,titl){
    q <- ggscat(m, columns = 1:ncol(m),color = NULL,size=10)+
      stat_bin_hex(geom="point", aes(color=(..count..)), binwidth=c(.05,.05), size=2, shape=16, stroke=0)+
      scale_color_viridis(option="inferno",direction = 1)+
      theme(axis.text.x = element_text(size=18,colour = "black"),
            axis.text.y = element_text(size=18,colour = "black"),
            axis.title =  element_text(size=20,color="black"),
            strip.text = element_text(color="black",size=20))+
      theme(legend.position = "none")+
      theme(panel.background = element_rect(fill = "white", colour = "black"),
            panel.grid = element_blank())+
      geom_smooth(method='lm',formula=y~x,col="red")+
      xlab("factor occupancy,log2")+ylab("factor occupancy, log2")+
      ggtitle(paste0(titl))+
      theme(plot.title = element_text(colour ="black",size = 18,hjust = 0.5))
    q
  }
  
  OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr"
  Inp.file=paste0(OUTDIR,"/ChIPSignalAlongGenes.RData")
  norm.file=paste0(OUTDIR,"/NormalizationFactors.txt")
  load(paste0(Inp.file)) 
  #cf[,11:ncol(cf)] <- apply(cf[,11:ncol(cf)],2,function(x) x*1e6/sum(x))
  m <- cf[cf$seqnames%in%spiked.chr & cf$width>50,11:ncol(cf)]
  normF <- read.delim(paste0(norm.file),header=T,stringsAsFactors = F)
  #normF$sample <- gsub("-",".",normF$sample)
  normF$sample%in%names(m)
  m <- m[,normF$sample]
  normF$sample==names(m)
  m1 <- matrix(NA,nrow = nrow(m),ncol = ncol(m))
  for(i in 1:ncol(m)){
    m1[,i] <- m[,i]*normF$normF[i]
  }
  m1 <- as.data.frame(m1)
  h <- names(m)
  h <- gsub("8WG16ChIP","Pol2",h)
  h <- gsub("MycChIP","Spt16",h)
  h <- gsub("FlagChIP","Spt6",h)
  h <- gsub("spt6.YW","spt6-YW",h)
  h <- gsub("_r1_30",",30°C,R1",h)
  h <- gsub("_r2_30",",30°C,R2",h)
  h <- gsub("_r1_37",",37°C,R1",h)
  h <- gsub("_r2_37",",37°C,R2",h)
  h <- gsub("spt6-YW.pob3","double",h)
  names(m) <- names(m1) <- h
  for(f in c("Pol2","Spt16","Inp")){
    for(tmp in c(30,37)){
      n <- names(m) ## total 16 inputs + 48 ChIPs
      n <- n[grep(f,n)]
      n <- n[grep(tmp,n)]
      z1 <- log2(m[ ,names(m)%in%n]+1)
      fout = paste0("Spombe_",f,"_SampleCorrelations_preSpikein_",tmp,"_July.png")
      cat(fout,"\n")
      cat(dim(z1),"\n")
      png(file = paste0(fout),width = 2000,height = 2000)
      print(myplot(z1,titl = f))
      dev.off()
      z1 <- log2(m1[ ,names(m1)%in%n]+1)
      fout = paste0("Spombe_",f,"_SampleCorrelations_postSpikein_",tmp,".png")
      cat(fout,"\n\n")
      cat(dim(z1),"\n")
      png(file = paste0(fout),width = 2000,height = 2000)
      print(myplot(z1,titl = f))
      dev.off()
    }    
  }
  
  
  OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr"
  Inp.file=paste0(OUTDIR,"/ChIPSignalAlongGenes.RData")
  norm.file=paste0(OUTDIR,"/NormalizationFactors.txt")
  load(paste0(Inp.file)) 
  m <- cf[cf$seqnames%in%spiked.chr & cf$width>50,11:ncol(cf)]
  normF <- read.delim(paste0(norm.file),header=T,stringsAsFactors = F)
  #normF$sample <- gsub("-",".",normF$sample)
  normF$sample%in%names(m)
  m <- m[,normF$sample]
  normF$sample==names(m)
  m1 <- matrix(NA,nrow = nrow(m),ncol = ncol(m))
  for(i in 1:ncol(m)){
    m1[,i] <- m[,i]*normF$normF[i]
  }
  m1 <- as.data.frame(m1)
  h <- names(m)
  h <- gsub("8WG16ChIP","Pol2",h)
  h <- gsub("MycChIP","Spt16",h)
  h <- gsub("FlagChIP","Spt6",h)
  h <- gsub("spt6.YW","spt6-YW",h)
  h <- gsub("_r1_30",",30°C,R1",h)
  h <- gsub("_r2_30",",30°C,R2",h)
  h <- gsub("_r1_37",",37°C,R1",h)
  h <- gsub("_r2_37",",37°C,R2",h)
  h <- gsub("spt6-YW.pob3","double",h)
  names(m) <- names(m1) <- h
  for(f in c("Spt6","Pol2","Inp")){
    for(tmp in c(30,37)){
      n <- names(m) ## total 16 inputs + 48 ChIPs
      n <- n[grep(f,n)]
      n <- n[grep(tmp,n)]
      z1 <- log2(m[ ,names(m)%in%n]+1)
      fout = paste0("Spombe_",f,"_SampleCorrelations_preSpikein_",tmp,".png")
      cat(fout,"\n")
      cat(dim(z1),"\n")
      png(file = paste0(fout),width = 2000,height = 2000)
      print(myplot(z1,titl = f))
      dev.off()
      z1 <- log2(m1[ ,names(m1)%in%n]+1)
      fout = paste0("Spombe_",f,"_SampleCorrelations_postSpikein_",tmp,".png")
      cat(fout,"\n\n")
      cat(dim(z1),"\n")
      png(file = paste0(fout),width = 2000,height = 2000)
      print(myplot(z1,titl = f))
      dev.off()
    }    
  }
}

### CNV analyses using median signal
if(F){
  #OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr"
  OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr"
  Inp.file=paste0(OUTDIR,"/CountMatrix_SignalRegions_5p.RData")
  out.file=paste0(OUTDIR,"/CNV_MedianReadCountplots.pdf")
  
  load(pste0(Inp.file))
  cf[,6:ncol(cf)] <- cf[,6:ncol(cf)]*1000/cf$width
  cf[,6:ncol(cf)] <- apply(cf[,6:ncol(cf)],2,function(x){round(x*1e6/sum(x),3)  })
  
  mn <- c()
  for(chr in levels(as.factor(cf$seqnames))){
    mn <- rbind(mn,data.frame(chr=chr,
                              med=apply(cf[cf$seqnames==chr,6:ncol(cf)],2,median),
                              sample=names(cf)[6:ncol(cf)]))
  }
  rownames(mn) <- NULL
  
  
  mn$variable <- "autosomes"
  mn[mn$chr%in%spiked.chr,]$variable <- "heterosomes"
  n <- mn[grep("Inp",mn$sample),]
  n <- n[n$variable=="autosomes",]
  p <-ggplot(n,aes(x=sample, y=med, fill=variable)) + 
    geom_bar(stat="identity",position="dodge")+ facet_wrap(~chr,nrow=5)+
    xlab("") + ylab("median read count") + theme_minimal()+
    scale_fill_brewer(palette = "Dark2") +coord_flip()+
    theme(axis.text.x = element_text(angle = 0,size=14,color="black",hjust = 1),
          axis.text.y = element_text(size=14,color="black",hjust = 1),
          axis.title = element_text(size=16,color="black"),
          legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_text(size=18,color="black"))
  
  pdf(paste0(out.file),width = 8,height = 16)
  p
  dev.off()
}

## Read summary
if(F){
  rm(list=ls())
  d <- read.delim("readSummary.txt.txt",header=F)
  head(d)
  d$step <- ""
  d[grep('raw',d$V1),]$step <- "raw"
  d[grep('clean',d$V1),]$step <- "trimmed"
  d[grep('sort.bam',d$V1),]$step <- "InBamFile"
  d[grep('map.sort.bam',d$V1),]$step <- "mapped"
  table(d$step)
  
  d$V1 <- gsub("_raw.fq.gz","",d$V1)
  d$V1 <- gsub("_clean.fq.gz","",d$V1)
  d$V1 <- gsub(".map.sort.bam","",d$V1)
  d$V1 <- gsub(".sort.bam","",d$V1)
  d$V1 <- gsub(" ","",d$V1)
  
  d2 <- d[d$V3=="August",]
  d1 <- d[d$V3=="July",]
  
  require(reshape2)
  d1 <- reshape2::dcast(d1,V1~step,value.var = "V2")
  d2 <- reshape2::dcast(d2,V1~step,value.var = "V2")
  
  d1$raw <- d1$raw/4
  d1$trimmed <- d1$trimmed/4
  d2$raw <- d2$raw/4
  d2$trimmed <- d2$trimmed/4
  d1$month <- "July"
  d2$month <- 'August'
  
  d <- rbind(d1,d2)
  rm(d1,d2)
  head(d)
  
  d$rawseq <- paste0(round(d$raw/1e6,2),"M")
  d$percentofRawMapped <- paste0(round(d$mapped*100/d$raw,2),"% mapped")
  d$trimseq <- paste0(round(d$trimmed/1e6,2),"M")
  d$mapseq <- paste0(round(d$mapped/1e6,2),"M")
  names(d)[1] <- "sample"
  write.table(d,file="ReadSummary_Formatted.txt",sep="\t",quote = F,row.names = F)
  
}

### sample PCA analyses and correlations in groups
if(F){
  rm(list=ls())
  if(F){
    OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr"
    load(paste0(OUTDIR,"/ChIPSignalAlongGenes.RData")) 
    cf1 <- cf;rm(cf)
    load(paste0(OUTDIR,"/ChIPSignalAlongGenes_1st500bp.RData")) 
    cf[,11:ncol(cf)] <- round(cf[,11:ncol(cf)]/cf1[,11:ncol(cf1)],2)
    save(cf,file=paste0(OUTDIR,"/ChIPSignalAlongGenes_Ratio500byGL.RData")) 
    rm(cf,cf1)
  }
  
  OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr"
  Inp.file=paste0(OUTDIR,"/ChIPSignalAlongGenes.RData")
  norm.file=paste0(OUTDIR,"/NormalizationFactors.txt")
  load(paste0(Inp.file)) 
  cf <- cf[cf$species=="Scer",]
  m <- cf[cf$species=="Scer",11:ncol(cf)]
  normF <- read.delim(paste0(norm.file),header=T,stringsAsFactors = F)
  #normF$sample <- gsub("-",".",normF$sample)
  normF$sample%in%names(m)
  m <- m[,as.character(normF$sample)]
  normF$sample==names(m)
  m1 <- matrix(NA,nrow = nrow(m),ncol = ncol(m))
  m <- data.frame(apply(m,2,function(x) round(x*1e7/sum(x),2) ))
  for(i in 1:ncol(m)){
    m1[,i] <- m[,i]*normF$normF[i]
  }
  m1 <- as.data.frame(m1)
  names(m1) <- names(m)
  rownames(m) <- rownames(m1) <- cf$tracking_id
  rm(cf,normF,i,Inp.file,norm.file)
  
  OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr"
  Inp.file=paste0(OUTDIR,"/ChIPSignalAlongGenes.RData")
  norm.file=paste0(OUTDIR,"/NormalizationFactors.txt")
  load(paste0(Inp.file)) 
  cf <- cf[cf$species=="Scer",]
  n <- cf[cf$species=="Scer",11:ncol(cf)]
  normF <- read.delim(paste0(norm.file),header=T,stringsAsFactors = F)
  #normF$sample <- gsub("-",".",normF$sample)
  normF$sample%in%names(n)
  n <- n[,as.character(normF$sample)]
  normF$sample==names(n)
  n1 <- matrix(NA,nrow = nrow(n),ncol = ncol(n))
  n <- as.data.frame(apply(n,2,function(x) round(x*1e7/sum(x),2) ))
  for(i in 1:ncol(n)){
    n1[,i] <- n[,i]*normF$normF[i]
  }
  n1 <- as.data.frame(n1)
  names(n1) <- names(n)
  rownames(n) <- rownames(n1) <- cf$tracking_id
  rm(cf,normF,i,Inp.file,norm.file)
  
  names(m) <- paste0("Spt16_",names(m))
  names(n) <- paste0("Spt6_",names(n))
  names(m1) <- paste0("Spt16_",names(m1))
  names(n1) <- paste0("Spt6_",names(n1))
  m <- cbind(m,n)
  m1 <- cbind(m1,n1)
  rm(n,n1)
  
  ######### Sample correlations
  hmcols <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  pdf("ChIPSampleCorrelations_groupwise.pdf",width = 12,height = 12)
  {
    n <- names(m) ## total 16 inputs + 48 ChIPs
    n <- n[grep("8WG16ChIP",n)]
    n <- n[grep("Spt16_",n)]
    z1 <- round(cor(log2(m[ ,names(m)%in%n]+1)),2)
    rownames(z1) <- gsub('Spt16_','',rownames(z1))
    colnames(z1) <- gsub('Spt16_','',colnames(z1))
    heatmap.2(main = "Pol2 ChIP correlations (July)",z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.01,0.01),sepcolor = "black",colsep=0:ncol(z1),rowsep=0:nrow(z1),
              key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
    
    n <- names(m) ## total 16 inputs + 48 ChIPs
    n <- n[grep("8WG16ChIP",n)]
    n <- n[grep("Spt6_",n)]
    z1 <- round(cor(log2(m[ ,names(m)%in%n]+1)),2)
    rownames(z1) <- gsub('Spt6_','',rownames(z1))
    colnames(z1) <- gsub('Spt6_','',colnames(z1))
    heatmap.2(main = "Pol2 ChIP correlations (August)",z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.01,0.01),sepcolor = "black",colsep=0:ncol(z1),rowsep=0:nrow(z1),
              key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
    
    n <- names(m) ## total 16 inputs + 48 ChIPs
    n <- n[grep("FlagChIP",n)]
    n <- n[grep("Spt6_",n)]
    z1 <- round(cor(log2(m[ ,names(m)%in%n]+1)),2)
    rownames(z1) <- gsub('Spt6_','',rownames(z1))
    colnames(z1) <- gsub('Spt6_','',colnames(z1))
    heatmap.2(main = "FlagChIP correlations",z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.01,0.01),sepcolor = "black",colsep=0:ncol(z1),rowsep=0:nrow(z1),
              key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
    
    n <- names(m) ## total 16 inputs + 48 ChIPs
    n <- n[grep("MycChIP",n)]
    n <- n[grep("Spt16_",n)]
    z1 <- round(cor(log2(m[ ,names(m)%in%n]+1)),2)
    rownames(z1) <- gsub('Spt16_','',rownames(z1))
    colnames(z1) <- gsub('Spt16_','',colnames(z1))
    heatmap.2(main = "MycChIP correlations",z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.01,0.01),sepcolor = "black",colsep=0:ncol(z1),rowsep=0:nrow(z1),
              key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
    
    n <- names(m) ## total 16 inputs + 48 ChIPs
    n <- n[grep("Inp",n)]
    n <- n[grep("Spt16_",n)]
    z1 <- round(cor(log2(m[ ,names(m)%in%n]+1)),2)
    rownames(z1) <- gsub('Spt16_','',rownames(z1))
    colnames(z1) <- gsub('Spt16_','',colnames(z1))
    heatmap.2(main = "input correlations (July)",z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.01,0.01),sepcolor = "black",colsep=0:ncol(z1),rowsep=0:nrow(z1),
              key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
    
    
    n <- names(m) ## total 16 inputs + 48 ChIPs
    n <- n[grep("Inp",n)]
    n <- n[grep("Spt6_",n)]
    z1 <- round(cor(log2(m[ ,names(m)%in%n]+1)),2)
    rownames(z1) <- gsub('Spt6_','',rownames(z1))
    colnames(z1) <- gsub('Spt6_','',colnames(z1))
    heatmap.2(main = "input correlations (August)",z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.01,0.01),sepcolor = "black",colsep=0:ncol(z1),rowsep=0:nrow(z1),
              key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
    
  }
  dev.off()
  rm(z1,n,hmcols)
  
  ######### Sample clustering
  pdf("ChIPSampleClustering_groupwise.pdf",width = 8,height = 12)
  {
    n <- names(m) ## total 16 inputs + 48 ChIPs
    n <- n[grep("8WG16ChIP",n)]
    n <- n[grep("Spt16_",n)]
    z1 <- round(cor(log2(m[ ,names(m)%in%n]+1)),2)
    rownames(z1) <- gsub('Spt16_','',rownames(z1))
    colnames(z1) <- gsub('Spt16_','',colnames(z1))
    fit <- hclust(dist(z1, method = "euclidean"), method="ward.D") 
    plot(as.phylo(fit), type = "unrooted", cex = 0.8,no.margin = T,main = "Pol2 ChIP correlations (July)")
    
    n <- names(m) ## total 16 inputs + 48 ChIPs
    n <- n[grep("8WG16ChIP",n)]
    n <- n[grep("Spt6_",n)]
    z1 <- round(cor(log2(m[ ,names(m)%in%n]+1)),2)
    rownames(z1) <- gsub('Spt6_','',rownames(z1))
    colnames(z1) <- gsub('Spt6_','',colnames(z1))
    fit <- hclust(dist(z1, method = "euclidean"), method="ward.D") 
    plot(as.phylo(fit), type = "unrooted", cex = 0.8,no.margin = T,main = "Pol2 ChIP correlations (August)")
    
    n <- names(m) ## total 16 inputs + 48 ChIPs
    n <- n[grep("FlagChIP",n)]
    n <- n[grep("Spt6_",n)]
    z1 <- round(cor(log2(m[ ,names(m)%in%n]+1)),2)
    rownames(z1) <- gsub('Spt6_','',rownames(z1))
    colnames(z1) <- gsub('Spt6_','',colnames(z1))
    fit <- hclust(dist(z1, method = "euclidean"), method="ward.D") 
    plot(as.phylo(fit), type = "unrooted", cex = 0.8,no.margin = T,main = "FlagChIP correlations")
    
    n <- names(m) ## total 16 inputs + 48 ChIPs
    n <- n[grep("MycChIP",n)]
    n <- n[grep("Spt16_",n)]
    z1 <- round(cor(log2(m[ ,names(m)%in%n]+1)),2)
    rownames(z1) <- gsub('Spt16_','',rownames(z1))
    colnames(z1) <- gsub('Spt16_','',colnames(z1))
    fit <- hclust(dist(z1, method = "euclidean"), method="ward.D") 
    plot(as.phylo(fit), type = "unrooted", cex = 0.8,no.margin = T,main = "MycChIP correlations")
    
    n <- names(m) ## total 16 inputs + 48 ChIPs
    n <- n[grep("Inp",n)]
    n <- n[grep("Spt16_",n)]
    z1 <- round(cor(log2(m[ ,names(m)%in%n]+1)),2)
    rownames(z1) <- gsub('Spt16_','',rownames(z1))
    colnames(z1) <- gsub('Spt16_','',colnames(z1))
    fit <- hclust(dist(z1, method = "euclidean"), method="ward.D") 
    plot(as.phylo(fit), type = "unrooted", cex = 0.8,no.margin = T,main = "input correlations (July)")
    
    n <- names(m) ## total 16 inputs + 48 ChIPs
    n <- n[grep("Inp",n)]
    n <- n[grep("Spt6_",n)]
    z1 <- round(cor(log2(m[ ,names(m)%in%n]+1)),2)
    rownames(z1) <- gsub('Spt6_','',rownames(z1))
    colnames(z1) <- gsub('Spt6_','',colnames(z1))
    fit <- hclust(dist(z1, method = "euclidean"), method="ward.D") 
    plot(as.phylo(fit), type = "unrooted", cex = 0.8,no.margin = T,main = "input correlations (August)")
    
  }
  dev.off()
  rm(fit,z1,n)
  
  ############## PCA analyses
  z <- m
  z <- z[,-grep("Inp",names(z))]
  z <- as.matrix(z)
  z[!is.finite(z)] <- 0
  z <- as.data.frame(z)
  pl <- data.frame(id=names(z))
  pl$col <- "input (July2019)"
  pl[grep("Spt16_",pl$id),]$col <- "input (Aug2019)"
  pl[grep("8WG16ChIP",pl$id),]$col <- "Pol2 (Aug2019)"
  pl[c(9:12,17:20,25:28),]$col="Pol2 (July2019)"
  pl[grep("MycChIP",pl$id),]$col <- "Spt16"
  pl[grep("FlagChIP",pl$id),]$col <- "Spt6"
  pl <- merge(pl,x[,c("sample","lab")],by.x="id",by.y="sample",all.x=T)
  
  pca_data=prcomp(t(z))
  pca_data_perc=round(100*pca_data$sdev^2/sum(pca_data$sdev^2),1)
  df_pca_data=data.frame(PC1 = pca_data$x[,1], PC2 = pca_data$x[,2], 
                         sample = names(z), condition=pl$col,label=pl$id)
  df_pca_data$label <- gsub("Spt6_","",df_pca_data$label)
  df_pca_data$label <- gsub("Spt16_","",df_pca_data$label)
  
  
  #pdf(file = "A:/work/WinstonLab/Olga/ChIPSeq/SampleSVAplots_nonspikein.pdf",width = 10,height = 8)
  pdf(file = "A:/work/WinstonLab/Olga/ChIPSeq/SampleSVAplots_SigRatios.pdf",width = 12,height = 8)
  ggplot(df_pca_data, aes(PC1,PC2, color = condition,label=label))+geom_text_repel(box.padding = 0.1)+
    geom_point(size=8)+
    labs(x=paste0("PC1 (",pca_data_perc[1],")"), y=paste0("PC2 (",pca_data_perc[2],")"))+
    scale_color_brewer(palette = "Dark2")+
    theme(legend.position = "bottom")
  
  ggplot(df_pca_data, aes(PC1,PC2, color = condition,label=label))+
    geom_point(size=8)+
    labs(x=paste0("PC1 (",pca_data_perc[1],")"), y=paste0("PC2 (",pca_data_perc[2],")"))+
    scale_color_brewer(palette = "Dark2")+
    theme(legend.position = "bottom")
  dev.off()

   

  
  
  
}

### Mutation analysis
if(F){
  rm(list=ls())
  OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt6/"
  #OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt16/"
  
  d = rbind(read.delim(paste0(OUTDIR,"MutationAnalysis.sam"),header=F))
  d = d[!is.na(d$V6),]  ## maped reads
  head(d)
  
  d$V7 = gsub("M","",d$V7)
  d$V7 = as.numeric(d$V7)
  table(d$V7)
  d = d[d$V7>28,]  ## minimum seq len
  d$V5= as.numeric(as.vector(d$V5))
  
  gr = GRanges(c("SPT6_WT","SPT6_mutant","SPT5_WT","SPT5_mutant","POB3_WT","POB3_mutant","SPN1_WT","SPN1_mutant"),
               IRanges(start = c(42,42,66,66,99,99,56,56),end = c(53,53,66,66,99,99,56,56)),"+")
  d$end = d$V5+d$V7
  names(d)[4:8] = c("chr","start","score","width","strand")
  d = as(d,"GRanges")
  d = split(d,seqnames(d))
  gr = split(gr,seqnames(gr))
  
  for(i in 1:length(gr)){
    cat(names(gr)[i],"\n")
    z = subsetByOverlaps(d[[names(gr)[i]]], gr[[i]],ignore.strand=T)
    z$mut = as.numeric(start(gr[[i]]))
    d[[names(gr)[i]]] = z
  }
  
  d= unlist(d)
  names(d) = 1:length(d)
  table(d$V1, seqnames(d))
  
  d$md <- gsub("MD:Z:","",d$V14)
  d$nm <- as.numeric(gsub("NM:i:","",d$V15))
  d = as.data.frame(d) 
  
  nm      <- d$nm[1171]
  md      <- d$md[1171] #"0G15^GAC0T60T4^AA0C0"
  start   <- d$start[1171]
  mut <- d$mut[1171]
  
  d$res <- ""
  for(i in 1:length(d$seqnames)){
    nm=d$nm[i]
    md=d$md[i]
    start=d$start[i]
    mut=d$mut[i]
    
    cat(start,"\n")   
    # First crucial step is to split the md string at each substitution/deletion operation
    md.gsub <- gsub("([\\^]*[ACGT]+)[0]*", " \\1 ", md)
    # split each operation using strsplit
    md.spl  <- strsplit(md.gsub, "[ ]+")[[1]]
    this    <- as.integer()
    # since its a relatively time-consuming operation
    # use NM flag to calculate position of mismatches
    # only if there are mismatches
    if (nm != 0) {
      this <- lapply(md.spl, function(y) { 
        if (!is.na(as.numeric(y))) {
          o <- rep("M", as.numeric(y))
        } else if( length(grep("\\^", y)) > 0) {
          o <- rep("D", nchar(y) - 1)
        } else if (nchar(y) == 1) {
          o <- rep("MM", 1)
        }
      })
      this <- do.call(c, this)
      # after this step, we have a vector of length = 
      # read length and each position telling if its a 
      # match(M), mismatch(MM) or deletion (D)
      this <- which(this == "MM")
      this <- this+start
      if(mut%in%this){
        d$res[i] = paste0("TRUE|",mut,"|",paste(this,collapse=","))
      }else{
        d$res[i] = paste0("FALSE|",mut,"|",paste(this,collapse=","))
      }  
    }else{
      d$res[i] = paste0("FALSE|",mut,"|-")
    }
  }
  d1 = d[grep("TRUE",d$res),]
  d = d[-as.numeric(rownames(d1)),]
  
  mt <- as.data.frame(table(d$V1,d$seqnames))
  mt <- reshape2::dcast(mt,Var1~Var2,value.var = "Freq")
  
  mt$Var1 <- gsub("[{]","",mt$Var1)
  mt$Var1 <- gsub("[}]","",mt$Var1)
  mt$Var1 <- gsub("_clean","",mt$Var1)
  names(mt)[1] <- "sequence.file"
  mt$profile <- gsub("_[1234567]$","",mt$sequence.file)
  write.table(mt,file = paste0(OUTDIR,"MutationMappingReadSummary.txt"),quote = F,sep = "\t",row.names = F)
}

### TSS-centered metagene plots with/without merging the replicates 
 ###  DF
if(F){
  ## bring the functions from MultgenHeatMeta.R
  rm(list=ls())
  mymetagenefun_tss <- function(l,genes,facet=NULL,normF=c(1,1)){
    l[[1]] <- subset(l[[1]],rownames(l[[1]])%in%genes$tracking_id)
    l[[2]] <- subset(l[[2]],rownames(l[[2]])%in%genes$tracking_id)
    
    l[[1]] <- l[[1]]*normF[1]
    l[[2]] <- l[[2]]*normF[2]
    
    depl1 <- l[[1]]
    depl2 <- l[[2]]
    avgdepl <- round((depl1+depl2)/2,2)  ## average signal after spikein normalizationnnnnn
    
    if(!is.null(facet)){
      depl1 <-subset(depl1,rownames(depl1)%in%facet$tracking_id)
      depl2 <-subset(depl2,rownames(depl2)%in%facet$tracking_id)
    }
    get.avg <- function(m,nF){
      df_m <- cbind(1:ncol(m),
                    as.data.frame(apply(m,2,mean,na.rm=TRUE)),
                    as.data.frame(apply(m,2,sd,na.rm=TRUE)), 
                    as.data.frame(apply(m,2, function(x) sd(x,na.rm=TRUE)/sqrt(length(x))))
      )
      names(df_m) <- c("pos","signal", "sd","se")
      df_m$ci <- df_m$se*(1.96)
      return(df_m)
    }
    pl <- rbind(cbind(get.avg(depl1),condition="Repl-1"),
                cbind(get.avg(depl2),condition="Repl-2"),
                cbind(get.avg(avgdepl),condition="Average")
           )
    return(pl)
  }
  
  load(file="A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Gene groups.RData")
  rp <- rp[rp$gene!="RPP1",]
  rp <- rp[-grep("MRP",rp$gene),]
  
  genotypes <- c("WT",  "pob3","spt6-YW", "spt6-YW-pob3") #,   "spt6-YW_Flag","WT_Flag",  "pob3_Flag",     "spt6-YW-pob3_Flag")
  pl <- c()
  ## July
  res <- read.delim("A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/NormalizationFactors.txt",header=T,stringsAsFactors = F)
  chips <- c("8WG16",  "Myc") #,   "spt6-YW_Flag","WT_Flag",  "pob3_Flag",     "spt6-YW-pob3_Flag")
  for(tmp in c(30,37)){
    for(chip in chips){
      for(gt in genotypes){
        f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/MultigeneHeatmap_TSS_",gt,"_",chip,".RData",sep="")
        if(file.exists(f)){
          cat(gt,"\t",chip,"\t",tmp,"\n")
          cat(f,"\n")
          load(f)
          if(tmp==30){
            myname = paste0(gt,"_",chip,"_30°C")
            l <- l[1:2]
            normF = res[res$Factor==paste0(gt,"_",chip),]$normF[1:2]
          }else{
            myname = paste0(gt,"_",chip,"_37°C")
            l <- l[3:4]
            normF = res[res$Factor==paste0(gt,"_",chip),]$normF[3:4]
          }
          cat(normF,"\n")
          lx <- l
          lx[[1]] <- lx[[1]]*1e6/sum(lx[[1]])
          lx[[2]] <- lx[[2]]*1e6/sum(lx[[2]])
          m <- mymetagenefun_tss(lx,genes=genes)
          m$profile <- myname
          m$month="July"
          m$spike="LibSizeNorm"
          pl <- rbind(pl,m)
          m <- mymetagenefun_tss(l,genes=genes,normF=normF)
          m$profile <- myname
          m$month="July"
          m$spike="SpininNorm"
          pl <- rbind(pl,m)
          rm(m,l,f)
        }else{
          cat("no such file: ", f,"\n")
        }
      }
    }
  }
  ## Aug
  res <- read.delim("A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/NormalizationFactors.txt",header=T,stringsAsFactors = F)
  chips <- c("8WG16",  "Flag") #,   "spt6-YW_Flag","WT_Flag",  "pob3_Flag",     "spt6-YW-pob3_Flag")
  for(tmp in c(30,37)){
    for(chip in chips){
      for(gt in genotypes){
        f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/MultigeneHeatmap_TSS_",gt,"_",chip,".RData",sep="")
        if(file.exists(f)){
          cat(gt,"\t",chip,"\t",tmp,"\n")
          cat(f,"\n")
          load(f)
          if(tmp==30){
            myname = paste0(gt,"_",chip,"_30°C")
            l <- l[1:2]
            normF = res[res$Factor==paste0(gt,"_",chip),]$normF[1:2]
          }else{
            myname = paste0(gt,"_",chip,"_37°C")
            l <- l[3:4]
            normF = res[res$Factor==paste0(gt,"_",chip),]$normF[3:4]
          }
          cat(myname,"\n")
          cat(normF,"\n")
          lx <- l
          lx[[1]] <- lx[[1]]*1e6/sum(lx[[1]])
          lx[[2]] <- lx[[2]]*1e6/sum(lx[[2]])
          m <- mymetagenefun_tss(lx,genes=genes)
          m$profile <- myname
          m$month="Aug"
          m$spike="LibSizeNorm"
          pl <- rbind(pl,m)
          m <- mymetagenefun_tss(l,genes=genes,normF=normF)
          m$profile <- myname
          m$month="Aug"
          m$spike="SpininNorm"
          pl <- rbind(pl,m)
          rm(m,l,f)
        }else{
          cat("no such file: ", f,"\n")
        }
      }
    }
  }
  
  
  #########################
  myplot <- function(pl1,titl){
    ggplot(pl1, aes(x=pos,y=signal,ymax=signal+se, ymin=signal-se,col=as.factor(profile)))+
      facet_wrap(~spike,nrow = 5,scales = "free")+
      geom_vline(xintercept = c(50),colour="gray50", linetype = "dashed") +
      geom_linerange(col='gray') + geom_line(size=1.2) +
      #scale_color_manual(values =col) +
      scale_color_brewer(palette = "Dark2") +
      scale_y_continuous(expand = c(0,0))+
      scale_x_continuous(breaks =c(1,50,150,250,350,450), labels=c("","TSS","1kb","2kb","3kb","4kb"),expand = c(0,0)) +ylab("normalized signal") + xlab("") +
      theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = NA),axis.line = element_line(size = 0.7, colour = "black"))+
      theme(axis.text.x = element_text(colour = 'black',size=15),
            axis.text.y = element_text(colour = 'black',size=15), 
            axis.title.y = element_text(colour = 'black',size=18),
            strip.text.x = element_text(colour = 'black',size=16),
            strip.background = element_blank())+
      theme(legend.text =element_text(colour = 'black',size=14),
            legend.title = element_text(colour = 'black',size=15),
            legend.position = "bottom",
            plot.title = element_text(size=20,hjust = 0.5,vjust = 0.5))+
      ggtitle(titl)
    
  }
  
  pl$profile <- gsub("spt6-YW-pob3","double",pl$profile)
  pl$profile <- gsub("8WG16","Pol2",pl$profile)
  pl$profile <- gsub("Myc","Spt16",pl$profile)
  pl$profile <- gsub("Flag","Spt6",pl$profile)
  pl$spike <- gsub("SpininNorm","Spike-in scaled",pl$spike)
  pl$spike <- gsub("LibSizeNorm","Library size scaled",pl$spike)
  pdf("TSS-centered_Metageneplots_RawSignal.pdf",width = 8,height = 5)
  for(chip in c("Pol2","Spt16")){
    for(tmp in c("30°C","37°C")){
      pl1 <- pl[pl$condition=="Average" & pl$month=="July",]
      pl1 <- pl1[grepl(chip,pl1$profile,fixed = T) & grepl(tmp,pl1$profile,fixed = T),]
      pl1$profile <- gsub(paste0("_",chip,"_",tmp),"",pl1$profile)
      print(myplot(pl1,paste0(chip," ChIP (July),",tmp)))
      
      pl1 <- pl[pl$condition!="Average" & pl$month=="July",]
      pl1 <- pl1[grepl(chip,pl1$profile,fixed = T) & grepl(tmp,pl1$profile,fixed = T),]
      pl1$profile <- gsub(paste0("_",chip,"_",tmp),"",pl1$profile)
      pl1$condition <- gsub("Repl-","",pl1$condition)
      pl1$profile <- paste0(pl1$profile,"-",pl1$condition)
      print(myplot(pl1,paste0(chip," ChIP (July),",tmp)))
      
    }
  }
  for(chip in c("Pol2","Spt6")){
    for(tmp in c("30°C","37°C")){
      pl1 <- pl[pl$condition=="Average" & pl$month=="Aug",]
      pl1 <- pl1[grepl(chip,pl1$profile,fixed = T) & grepl(tmp,pl1$profile,fixed = T),]
      pl1$profile <- gsub(paste0("_",chip,"_",tmp),"",pl1$profile)
      print(myplot(pl1,paste0(chip," ChIP (August),",tmp)))
      
      pl1 <- pl[pl$condition!="Average" & pl$month=="Aug",]
      pl1 <- pl1[grepl(chip,pl1$profile,fixed = T) & grepl(tmp,pl1$profile,fixed = T),]
      pl1$profile <- gsub(paste0("_",chip,"_",tmp),"",pl1$profile)
      pl1$condition <- gsub("Repl-","",pl1$condition)
      pl1$profile <- paste0(pl1$profile,"-",pl1$condition)
      print(myplot(pl1,paste0(chip," ChIP (August),",tmp)))
    }
  }
  dev.off()
  
}  
if(T){
  rm(list=ls())
  load(file="A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Gene groups.RData")
  rp <- rp[rp$gene!="RPP1",]
  rp <- rp[-grep("MRP",rp$gene),]
  
  pl <- c()
  gl <- c()
  pl.repl <- c()
  gl.repl <- c()
  ## July samples
  protein<- c("WT_8WG16",  "WT_Myc",  "pob3_8WG16",  "pob3_Myc",  "spt6-YW_Myc",  "spt6-YW_8WG16",
              "spt6-YW-pob3_8WG16",  "spt6-YW-pob3_Myc")
  res <- read.delim("A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/NormalizationFactors.txt",header = T,stringsAsFactors = F)
  for( prot in protein){
    f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/MultigeneHeatmap_TSS_",prot,".RData",sep="")
    if(file.exists(f)){
      cat(prot,"\n")
      load(f)
      cat(names(l)[1:2],"\n")
      cat(names(l)[3:4],"\n")
      for(j in 1:length(l)){
        nn <- as.character(names(l)[j])
        cat(" ",nn,"\t",res[res$sample==nn,]$normF,"\n")
        l[[nn]] <- l[[nn]]*res[res$sample==nn,]$normF
      }
      gl <- rbind(gl, 
                  cbind(mymetagenefun_tss(l,genes,facet=NULL),factor=paste0(prot),category="spikenorm"))
      gl.repl <- rbind(gl.repl, 
                  cbind(mymetagenefun_tss(l,genes,facet=NULL,merge.repl=F),factor=paste0(prot),category="spikenorm"))
      load(f)
      for(j in 1:length(l)){
        l[[j]] <- l[[j]]*1e6/sum(l[[j]],na.rm=T)
      }
      
      gl <- rbind(gl, 
                  cbind(mymetagenefun_tss(l,genes,facet=NULL),factor=paste0(prot),category="non-spikenorm"))
      gl.repl <- rbind(gl.repl, 
                       cbind(mymetagenefun_tss(l,genes,facet=NULL,merge.repl=F),factor=paste0(prot),category="non-spikenorm"))
      rm(l,f,prot,normF)
    }else{
      cat("file corresponding to ",prot, " does not exists! \n")
    }
  }
  ## July samples
  protein<- c("WT_8WG16",  "WT_Flag",  "pob3_8WG16",  "pob3_Flag",  "spt6-YW_Flag",  "spt6-YW_8WG16",
              "spt6-YW-pob3_8WG16",  "spt6-YW-pob3_Flag")
  res <- read.delim("A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/NormalizationFactors.txt",header = T,stringsAsFactors = F)
  for( prot in protein){
    f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/MultigeneHeatmap_TSS_",prot,".RData",sep="")
    if(file.exists(f)){
      cat(prot,"\n")
      load(f)
      cat(names(l)[1:2],"\n")
      cat(names(l)[3:4],"\n")
      for(j in 1:length(l)){
        nn <- as.character(names(l)[j])
        cat(" ",nn,"\t",res[res$sample==nn,]$normF,"\n")
        l[[nn]] <- l[[nn]]*res[res$sample==nn,]$normF
      }
      pl <- rbind(pl, 
                  cbind(mymetagenefun_tss(l,genes,facet=NULL),factor=paste0(prot),category="spikenorm"))
      pl.repl <- rbind(pl.repl, 
                       cbind(mymetagenefun_tss(l,genes,facet=NULL,merge.repl=F),factor=paste0(prot),category="spikenorm"))
      load(f)
      for(j in 1:length(l)){
        l[[j]] <- l[[j]]*1e6/sum(l[[j]],na.rm=T)
      }
      pl <- rbind(pl, 
                  cbind(mymetagenefun_tss(l,genes,facet=NULL),factor=paste0(prot),category="non-spikenorm"))
      pl.repl <- rbind(pl.repl, 
                       cbind(mymetagenefun_tss(l,genes,facet=NULL,merge.repl=F),factor=paste0(prot),category="non-spikenorm"))
      rm(l,f,prot,normF)
    }else{
      cat("file corresponding to ",prot, " does not exists! \n")
    }
  }
  
  pl$condition <- gsub("Depleted","30C",pl$condition)
  pl$condition <- gsub("Non-depleted","37C",pl$condition)
  pl.repl$condition <- gsub("Depleted","30C",pl.repl$condition)
  pl.repl$condition <- gsub("Non-depleted","37C",pl.repl$condition)
  pl.repl$repl <- gl.repl$repl <-pl$repl <- gl$repl <- 1
  pl.repl[grep("C-2",pl.repl$condition),]$repl <- 2
  pl.repl$condition <- gsub("-1|-2","",pl.repl$condition)
  pl.repl$factor <- paste0(pl.repl$factor,'-',pl.repl$repl)
  
  gl$condition <- gsub("Depleted","30C",gl$condition)
  gl$condition <- gsub("Non-depleted","37C",gl$condition)
  gl.repl$condition <- gsub("Depleted","30C",gl.repl$condition)
  gl.repl$condition <- gsub("Non-depleted","37C",gl.repl$condition)
  gl.repl$repl <- 1
  gl.repl[grep("C-2",gl.repl$condition),]$repl <- 2
  gl.repl$condition <- gsub("-1|-2","",gl.repl$condition)
  gl.repl$factor <- paste0(gl.repl$factor,'-',gl.repl$repl)
  head(pl.repl)
  save(genes, pl,gl,pl.repl,gl.repl, file="A:/work/WinstonLab/Olga/ChIPSeq/TSS-Centered_MetageneDF_artScaled.RData")
}
 ## plots
if(F){
  OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq"
  require(ggplot2)
  require(ggrepel)
  require(gplots)
  myplot <- function(pl1,titl){
    ggplot(pl1, aes(x=pos,y=signal,ymax=signal+se, ymin=signal-se,col=as.factor(factor)))+
      facet_wrap(~condition,nrow = 5,scales = "free")+
      geom_vline(xintercept = c(50),colour="gray50", linetype = "dashed") +
      geom_linerange(col='gray') + geom_line(size=1.2) +
      scale_color_manual(values =col) +
      #scale_color_brewer(palette = "Dark2") +
      scale_y_continuous(expand = c(0,0))+
      scale_x_continuous(breaks =c(1,50,150,250,350,450), labels=c("","TSS","1kb","2kb","3kb","4kb"),expand = c(0,0)) +ylab("normalized signal") + xlab("") +
      theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = NA),axis.line = element_line(size = 0.7, colour = "black"))+
      theme(axis.text.x = element_text(colour = 'black',size=15),
            axis.text.y = element_text(colour = 'black',size=15), 
            axis.title.y = element_text(colour = 'black',size=18),
            strip.text.x = element_text(colour = 'black',size=16),
            strip.background = element_blank())+
      theme(legend.text =element_text(colour = 'black',size=14),
            legend.title = element_text(colour = 'black',size=15),
            legend.position = "bottom",
            plot.title = element_text(size=20,hjust = 0.5,vjust = 0.5))+
      ggtitle(titl)
    
  }
  
  ## pl serries August,
  ## gl series July
  
  load(file="A:/work/WinstonLab/Olga/ChIPSeq/TSS-Centered_MetageneDF_artScaled.RData")
  
  pdf(file = paste0(OUTDIR,"/Metagenplots_TSS_replicates.pdf"),height = 8,width = 10)
  {
    pl1 <- pl.repl[grep("_Flag",pl.repl$factor),]
    pl1 <- pl1[pl1$category=="spikenorm",]
    print(myplot(pl1,titl = "Flag ChIP, Spikein normalized"))
    
    pl1 <- pl.repl[grep("_Flag",pl.repl$factor),]
    pl1 <- pl1[pl1$category=="non-spikenorm",]
    print(myplot(pl1,titl = "Flag ChIP, non-spikein normalized"))
    
    pl1 <- pl.repl[grep("_8WG16",pl.repl$factor),]
    pl1 <- pl1[pl1$category=="spikenorm",]
    print(myplot(pl1,titl = "8WG16 ChIP, Spikein normalized"))
    
    pl1 <- pl.repl[grep("_8WG16",pl.repl$factor),]
    pl1 <- pl1[pl1$category=="non-spikenorm",]
    print(myplot(pl1,titl = "8WG16 ChIP, non-spikein normalized"))
    
    
    pl1 <- gl.repl[grep("_Myc",gl.repl$factor),]
    pl1 <- pl1[pl1$category=="spikenorm",]
    print(myplot(pl1,titl = "Myc ChIP, Spikein normalized"))
    
    pl1 <- gl.repl[grep("_Myc",gl.repl$factor),]
    pl1 <- pl1[pl1$category=="non-spikenorm",]
    print(myplot(pl1,titl = "Myc ChIP, non-spikein normalized"))
    
    pl1 <- gl.repl[grep("_8WG16",gl.repl$factor),]
    pl1 <- pl1[pl1$category=="spikenorm",]
    print(myplot(pl1,titl = "8WG16 ChIP, Spikein normalized (July)"))
    
    pl1 <- gl.repl[grep("_8WG16",gl.repl$factor),]
    pl1 <- pl1[pl1$category=="non-spikenorm",]
    print(myplot(pl1,titl = "8WG16 ChIP, non-spikein normalized (July)"))
  }
  dev.off() 
  
  pdf(file = paste0(OUTDIR,"/Metagenplots_TSS_merged.pdf"),height = 8,width = 10)
  {
    pl1 <- pl[grep("_Flag",pl$factor),]
    pl1 <- pl1[pl1$category=="spikenorm",]
    print(myplot(pl1,titl = "Flag ChIP, Spikein normalized"))
    
    pl1 <- pl[grep("_Flag",pl$factor),]
    pl1 <- pl1[pl1$category=="non-spikenorm",]
    print(myplot(pl1,titl = "Flag ChIP, non-spikein normalized"))
    
    pl1 <- pl[grep("_8WG16",pl$factor),]
    pl1 <- pl1[pl1$category=="spikenorm",]
    print(myplot(pl1,titl = "8WG16 ChIP, Spikein normalized"))
    
    pl1 <- pl[grep("_8WG16",pl$factor),]
    pl1 <- pl1[pl1$category=="non-spikenorm",]
    print(myplot(pl1,titl = "8WG16 ChIP, non-spikein normalized"))
    
    
    pl1 <- gl[grep("_Myc",gl$factor),]
    pl1 <- pl1[pl1$category=="spikenorm",]
    print(myplot(pl1,titl = "Myc ChIP, Spikein normalized"))
    
    pl1 <- gl[grep("_Myc",gl$factor),]
    pl1 <- pl1[pl1$category=="non-spikenorm",]
    print(myplot(pl1,titl = "Myc ChIP, non-spikein normalized"))
    
    pl1 <- gl[grep("_8WG16",gl$factor),]
    pl1 <- pl1[pl1$category=="spikenorm",]
    print(myplot(pl1,titl = "8WG16 ChIP, Spikein normalized (July)"))
    
    pl1 <- gl[grep("_8WG16",gl$factor),]
    pl1 <- pl1[pl1$category=="non-spikenorm",]
    print(myplot(pl1,titl = "8WG16 ChIP, non-spikein normalized (July)"))
  }
  dev.off() 
  
  
  pl.repl$month <- "August"
  gl.repl$month <- "July"
  z <- rbind(pl.repl[grep("8WG16",pl.repl$factor),], gl.repl[grep("8WG16",gl.repl$factor),],
             pl.repl[grep("Flag",pl.repl$factor),], gl.repl[grep("Flag",gl.repl$factor),],
             pl.repl[grep("Myc",pl.repl$factor),], gl.repl[grep("Myc",gl.repl$factor),])
  z$fact <- gsub("-[0123456]$","",z$factor)
  z$col <- paste0(z$month,"_",z$repl)
  
  myplot <- function(pl1,titl){
    ggplot(pl1, aes(x=pos,y=signal,ymax=signal+se, ymin=signal-se,col=as.factor(col)))+
      facet_wrap(~condition+fact,nrow = 5,scales = "free")+
      geom_vline(xintercept = c(50),colour="gray50", linetype = "dashed") +
      geom_linerange(col='gray') + geom_line(size=1.2) +
      #scale_color_manual(values =col) +
      scale_color_brewer(palette = "Dark2") +
      scale_y_continuous(expand = c(0,0))+
      scale_x_continuous(breaks =c(1,50,150,250,350,450), labels=c("","TSS","1kb","2kb","3kb","4kb"),expand = c(0,0)) +ylab("normalized signal") + xlab("") +
      theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = NA),axis.line = element_line(size = 0.7, colour = "black"))+
      theme(axis.text.x = element_text(colour = 'black',size=15),
            axis.text.y = element_text(colour = 'black',size=15), 
            axis.title.y = element_text(colour = 'black',size=18),
            strip.text.x = element_text(colour = 'black',size=16),
            strip.background = element_blank())+
      theme(legend.text =element_text(colour = 'black',size=14),
            legend.title = element_text(colour = 'black',size=15),
            legend.position = "bottom",
            plot.title = element_text(size=20,hjust = 0.5,vjust = 0.5))+
      ggtitle(titl)
    
  }
  pdf(file = paste0(OUTDIR,"/ChIP_JulyAugustComparisons.pdf"),height = 12,width = 10)
  {
    z1 <- z[grep("_8WG16",z$factor),]
    z1 <- z1[z1$category=="spikenorm",]
    print(myplot(z1,titl = "8WG16 ChIP, Spikein normalized"))
    z1 <- z[grep("_8WG16",z$factor),]
    z1 <- z1[z1$category=="non-spikenorm",]
    print(myplot(z1,titl = "8WG16 ChIP, Non-spikein normalized"))
    
    z1 <- z[grep("_Flag",z$factor),]
    z1 <- z1[z1$category=="spikenorm",]
    print(myplot(z1,titl = "Flag ChIP, Spikein normalized"))
    z1 <- z[grep("_Flag",z$factor),]
    z1 <- z1[z1$category=="non-spikenorm",]
    print(myplot(z1,titl = "Flag ChIP, Non-spikein normalized"))
    
    z1 <- z[grep("_Myc",z$factor),]
    z1 <- z1[z1$category=="spikenorm",]
    print(myplot(z1,titl = "Myc ChIP, Spikein normalized"))
    z1 <- z[grep("_Myc",z$factor),]
    z1 <- z1[z1$category=="non-spikenorm",]
    print(myplot(z1,titl = "Myc ChIP, Non-spikein normalized"))
    
  }
  dev.off()
  
  z$col <- paste0(z$col,"_",z$condition)
  z$col <- gsub("July_","",z$col)
  z$col <- gsub("August_","",z$col)
  z$fact <- paste0(z$month,"_",z$fact)
  myplot <- function(pl1,titl){
    ggplot(pl1, aes(x=pos,y=signal,ymax=signal+se, ymin=signal-se,col=as.factor(col)))+
      facet_wrap(~fact,nrow = 5,scales = "free")+
      geom_vline(xintercept = c(50),colour="gray50", linetype = "dashed") +
      geom_linerange(col='gray') + geom_line(size=1.2) +
      #scale_color_manual(values =col) +
      scale_color_brewer(palette = "Dark2") +
      scale_y_continuous(expand = c(0,0))+
      scale_x_continuous(breaks =c(1,50,150,250,350,450), labels=c("","TSS","1kb","2kb","3kb","4kb"),expand = c(0,0)) +ylab("normalized signal") + xlab("") +
      theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = NA),axis.line = element_line(size = 0.7, colour = "black"))+
      theme(axis.text.x = element_text(colour = 'black',size=15),
            axis.text.y = element_text(colour = 'black',size=15), 
            axis.title.y = element_text(colour = 'black',size=18),
            strip.text.x = element_text(colour = 'black',size=16),
            strip.background = element_blank())+
      theme(legend.text =element_text(colour = 'black',size=14),
            legend.title = element_text(colour = 'black',size=15),
            legend.position = "bottom",
            plot.title = element_text(size=20,hjust = 0.5,vjust = 0.5))+
      ggtitle(titl)
    
  }
  pdf(file = paste0(OUTDIR,"/ChIP_JulyAugust_TempComparisons.pdf"),height = 12,width = 10)
  {
    z1 <- z[grep("_8WG16",z$factor),]
    z1 <- z1[z1$category=="spikenorm",]
    print(myplot(z1,titl = "8WG16 ChIP, Spikein normalized"))
    z1 <- z[grep("_8WG16",z$factor),]
    z1 <- z1[z1$category=="non-spikenorm",]
    print(myplot(z1,titl = "8WG16 ChIP, Non-spikein normalized"))
    
    z1 <- z[grep("_Flag",z$factor),]
    z1 <- z1[z1$category=="spikenorm",]
    print(myplot(z1,titl = "Flag ChIP, Spikein normalized"))
    z1 <- z[grep("_Flag",z$factor),]
    z1 <- z1[z1$category=="non-spikenorm",]
    print(myplot(z1,titl = "Flag ChIP, Non-spikein normalized"))
    
    z1 <- z[grep("_Myc",z$factor),]
    z1 <- z1[z1$category=="spikenorm",]
    print(myplot(z1,titl = "Myc ChIP, Spikein normalized"))
    z1 <- z[grep("_Myc",z$factor),]
    z1 <- z1[z1$category=="non-spikenorm",]
    print(myplot(z1,titl = "Myc ChIP, Non-spikein normalized"))
    
  }
  dev.off()

}

## multigen heatmap in groups with sorting etc, 
if(F){
  setwd("A:/work/WinstonLab/Olga/ChIPSeq")
  rm(list=ls())
  generate.kmeans.change.plots <- function(k=3,fllist,norm=T,
                                           prots = c("WT_8WG16\n @30°C July", "pob3_8WG16\n @30°C July", "spt6-YW_8WG16\n @30°C July","spt6-YW-pob3_8WG16\n @30°C July"),
                                           control=c(), perc=1,
                                           genes=NULL,row_order=NULL, img="",lg=F,row.z.transform=F){
    
    plotm <- list()
    for(prot in c(prots)){
      cf <- subset(fllist,fllist$category==prot)
      load(unique(cf$path))
      if(norm==T){
        for(j in 1:nrow(cf)){
          l[[as.character(cf$sample[j])]] <- l[[as.character(cf$sample[j])]]*cf$normF[j]
        }
      }
      l <- subset(l,names(l)%in%as.character(cf$sample))
      m <- l[[1]]
      for(j in 2:length(l)){
        m <- m+l[[j]]
      }
      m <- round(m/length(l),2)
      ## cap top 1% signal
      if(!is.null(perc)){
        h <- quantile(m,probs=perc,na.rm=T)
        m[m>h] <- h
      }
      if(row.z.transform==T){
        m = as.data.frame(t(apply(m,1,function(x) {
          if(sum(x)==0){
            return(x)
          }
          #round(x*100/sum(x),2)}
          round( (x-mean(x,na.rm=T))/sd(x,na.rm=T),2 )}
        )))
      }
      plotm[[paste0(prot)]] <- m
      rm(m,j,cf,l)
    }
    
    if(length(control)==length(prots)){
      plotc <- list()
      counter=0
      for(prot in c(control)){
        cf <- subset(fllist,fllist$category==prot)
        load(unique(cf$path))
        counter=counter+1
        if(norm==T){
          for(j in 1:nrow(cf)){
            l[[as.character(cf$sample[j])]] <- l[[as.character(cf$sample[j])]]*cf$normF[j]
          }
        }
        l <- subset(l,names(l)%in%as.character(cf$sample))
        m <- l[[1]]
        for(j in 2:length(l)){
          m <- m+l[[j]]
        }
        m <- round(m/length(l),2)
        if(!is.null(perc)){
          h <- quantile(m,probs=perc,na.rm=T)
          m[m>h] <- h
        }
        if(row.z.transform==T){
          m = as.data.frame(t(apply(m,1,function(x) {
            if(sum(x)==0){
              return(x)
            }
            round( (x-mean(x,na.rm=T))/sd(x,na.rm=T),2 )}
            #round(x*100/sum(x),2)}
          )))
        }
        plotc[[paste0(prot,counter)]] <- m
        rm(m,j,cf,l)
      }
      
      for(j in 1:length(plotm)){
        if(row.z.transform==T){
          z <- as.matrix(plotm[[j]] - plotc[[j]])
        }else{
          z <- as.matrix(round((plotm[[j]]*100+1)/(plotc[[j]]*100+1),2))
        }
        z[!is.finite(z)]<- 0
        plotm[[j]] <- as.data.frame(z)
      }
      rm(plotc,j,z)
    }
    
    
    makematrix <- function(a,lg=lg){
      #a <- a[,1:100]
      a <- as.matrix(a)
      if(lg==T){
        a <-log2(a+1)
      }
      dim(a) = dim(a)
      attr(a, "upstream_index") = 1:50
      attr(a, "target_index") = 51
      attr(a, "downstream_index") = 51:150
      attr(a, "extend") = c(1,1)  # it must be a vector of length two
      class(a) = c("normalizedMatrix", "matrix")
      attr(a, "signal_name") = "Spt6"
      attr(a, "target_name") = "Insert"
      attr(a,"target_is_single_point") = TRUE
      return(a)
    }
    
    mx <- c()
    mn <- c()
    for(j in 1:length(plotm)){
      zz <- plotm[[j]]
      zz <- zz[,1:150]
      zz <- zz[match(genes$tracking_id,rownames(zz)),]
      h <- quantile(zz,probs=perc,na.rm=T)
      mx <- append(mx,h )
      mn <- append(mn,min(zz,na.rm=T))
      zz[zz>h] <- h
      plotm[[j]] <- zz
    }
    if(lg==T){
      mx <- log2(max(mx,na.rm=T)+1)
      mn <- log2(min(mn,na.rm=T)+1)
    }else{
      mx <- (max(mx,na.rm=T))
      mn <- (min(mn,na.rm=T))
    }
    if(!is.finite(mn)){
      mn=0
    }
    
    #col = colorRamp2(seq(mn,mx,length.out = 10), viridis(10,alpha = 0.2,direction = -1,option = "magma") )
    #col = colorRamp2(seq(mn,mx,length.out = 9), brewer.pal(n = 9, name = "Blues") )
    #col = colorRamp2(quantile( as.matrix(plotm[[1]]), probs=seq(0,1,length.out = 10)), viridis(10,alpha = 0.1,direction = -1,option = "cividis") )
    
    if(!is.null(row_order)){
      #col = colorRamp2(seq(mn,mx,length.out = 10), viridis(10,alpha = 0.2,direction = -1,option = "magma") )
      z <- makematrix(plotm[[1]],lg = lg)
      col = colorRamp2(quantile(z,seq(0,1,length.out = 10)), viridis(10,alpha = 0.2,direction = -1,option = "magma") )
      htlist <- EnrichedHeatmap(z,axis_name = c("","TSS",""),column_title_gp = gpar(fontsize = 30),axis_name_gp = gpar(fontsize = 25),row_order = row_order, col = col, name = "1",
                                top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:10),yaxis_facing = "left")), column_title = names(plotm)[1])
      for(j in 2:length(plotm)){
        z <- makematrix(plotm[[j]],lg = lg)
        col = colorRamp2(quantile(z,seq(0,1,length.out = 10)), viridis(10,alpha = 0.2,direction = -1,option = "magma") )
        htlist <- htlist + 
            EnrichedHeatmap(z,axis_name = c("","TSS",""),column_title_gp = gpar(fontsize = 30),axis_name_gp = gpar(fontsize = 25),row_order = row_order, col = col, name = as.character(j),
            top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:10),yaxis_facing = "left")), column_title = names(plotm)[j])
      }
      
    }else{
      partition = paste0("cluster", kmeans(plotm[[ length(plotm) ]][,1:100], centers = k)$cluster)
      htlist = Heatmap(partition, col = structure(2:10, names = paste0(levels(as.factor(partition)) )), name = "",
                        show_row_names = FALSE, width = unit(1, "mm"))
      for(j in 1:length(plotm)){
        z <- makematrix(plotm[[j]],lg = lg)
        col = colorRamp2(quantile(z,seq(0,1,length.out = 10)), viridis(10,alpha = 0.2,direction = -1,option = "magma") )
        htlist <- htlist+
        EnrichedHeatmap(z,axis_name = c("","TSS",""),column_title_gp = gpar(fontsize = 30),
                        axis_name_gp = gpar(fontsize = 25), col = col, name = as.character(j),
                        top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:6),yaxis_facing = "left")), 
                        column_title = paste0(names(plotm)[j]))
      }
    }
    
    if(!is.null(genes)){
      htlist = htlist+
        Heatmap(log2(genes$pol2+1), row_order = row_order,column_names_gp = gpar(fontsize = 25),col = c("white","pink", "red4"), name = "log2(RNAPII)", 
                show_row_names = FALSE, width = unit(10, "mm"))+
        Heatmap(log2(genes$tf2d+1), row_order = row_order,column_names_gp = gpar(fontsize = 25),col = c("white","pink", "red4"), name = "log2(TFIIB)", 
                show_row_names = FALSE, width = unit(10, "mm"))+
      Heatmap(log2(genes$summit.rpkm+1),column_names_gp = gpar(fontsize = 25),row_order = row_order, col = c("white","pink", "red4"), name = "log2(FPKM)", 
              show_row_names = FALSE, width = unit(10, "mm"))
        
      g = genes
    }

    
    if(!is.null(row_order)){
      png(paste0(img),width = 1400,height = 800)
      draw(htlist, heatmap_legend_side = "bottom", gap = unit(2, "mm"))
      dev.off()
    }else if(!is.null(genes)){
      g$cluster = partition
      png(paste0(img),width = 1400,height = 800)
      draw(htlist, split = partition, heatmap_legend_side = "bottom", gap = unit(2, "mm"))
      dev.off()
      if(is.null(row_order)){
        return(as.data.frame(g))
      }
    }
  }    

  
  ## get gene object
  if(T){
    Inp.file="A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/ChIPSignalAlongGenes.RData"
    load(paste0(Inp.file)) 
    d <- cf[,c(6,grep("wt_8WG16",names(cf)))]
    d <- d[,c(1,grep("_30",names(d)))]
    Inp.file="A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/ChIPSignalAlongGenes.RData"
    load(paste0(Inp.file)) 
    e <- cf[,c(6,grep("WT_8WG16",names(cf)))]
    e <- e[,c(1,grep("_30",names(e)))]
    d <- merge(d,e,by='tracking_id')
    rm(e)
    d <- merge(d,cf[,c("tracking_id",'width')],by="tracking_id")
    
    chr.flt <- c("chrI","chrII","chrIII","chrIV","chrIX","chrV","chrVI","chrVII","chrVIII","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
    load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenes.RData")
    genes = genes[genes$species=="Scer" & genes$verification=="Verified",]
    genes = genes[!is.na(genes$chr),]
    genes = genes[genes$chr%in%chr.flt,]
    d <- d[d$tracking_id%in%genes$tracking_id,]
    d[,2:5] <- apply(d[,2:5],2,function(x) round(x*1e6/sum(x),2)) 
    d$pol2 <- apply(d[,2:5],1,function(x) prod(x)^(1/length(x)) )
    d$pol2 <-d$pol2*1000/d$width
    genes <- merge(genes,d[,c("tracking_id","pol2")],by="tracking_id",all.x=T)
    genes$pol2 <- ifelse(is.na(genes$pol2),0,genes$pol2)
    genes <- as(genes,"GRanges")
    rm(d,cf)

    cnt <- c()
    load("A:/work/WinstonLab/Olga/ChIPNexus/wt_1.gr")
    cnt <- cbind(cnt,countOverlaps(genes,gr,ignore.strand=T))
    rm(gr)
    load("A:/work/WinstonLab/Olga/ChIPNexus/wt_2.gr")
    cnt <- cbind(cnt,countOverlaps(genes,gr,ignore.strand=T))
    rm(gr)
    genes$tf2d = apply(cnt,1,mean)
    rm(cnt)
    
    g = genes
    rm(genes)
    ## TSS-seq signal
    load("A:/work/WinstonLab/Olga/MNase_Seq_Aug/gr/Olga_MNaseSeq_plotmatrix_TSS_2kbcentered_TSSseq.RData")
    rm(plotm)
    genes <- as.data.frame(genes)
    genes <- genes[,c("tracking_id","summit.rpkm")]
    g <- as.data.frame(g)
    g <- merge(g,genes,by="tracking_id",all.x=T)
    genes <- g
    genes$summit.rpkm <- ifelse(is.na(genes$summit.rpkm),0,genes$summit.rpkm)
    genes = as(genes,"GRanges")
    #genes$width <- width(genes)
    genes$pol2 <- ifelse(genes$pol2>quantile(genes$pol2,0.99,na.rm=T),quantile(genes$pol2,0.99,na.rm=T),genes$pol2)
    genes$summit.rpkm <- ifelse(genes$summit.rpkm>quantile(genes$summit.rpkm,0.99,na.rm=T),quantile(genes$summit.rpkm,0.99,na.rm=T),genes$summit.rpkm)
    genes$tf2d <- ifelse(genes$tf2d>quantile(genes$tf2d,0.99,na.rm=T),quantile(genes$tf2d,0.99,na.rm=T),genes$tf2d)
    rm(g)
  }
  
  r1 <- read.delim("A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/NormalizationFactors.txt",header = T,stringsAsFactors = F)
  r1$path <- paste0("A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/MultigeneHeatmap_TSS_",r1$Factor,".RData")
  r1$month <- "July"
  r2 <- read.delim("A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/NormalizationFactors.txt",header = T,stringsAsFactors = F)
  r2$path <- paste0("A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/MultigeneHeatmap_TSS_",r2$Factor,".RData")
  r2$month <- "Aug"
  res <- rbind(r1,r2)
  rm(r1,r2)
  res <- res[-grep("Inp",res$sample),]
  res$category <- paste0(res$Factor,"\n @",res$condition,"°C ",res$month)
  fllist <- res
  rm(res)
  
  
  ## row order 
  generate.kmeans.change.plots(fllist = fllist,norm = F,prots = c("WT_8WG16\n @30°C July", "pob3_8WG16\n @30°C July", "spt6-YW_8WG16\n @30°C July","spt6-YW-pob3_8WG16\n @30°C July"),
                               perc=1,genes=genes,row_order = order(width(genes),decreasing = T), img = "Leng_8WG16_30C_July.png",lg=F,row.z.transform = T)
  generate.kmeans.change.plots(fllist = fllist,norm = F,prots = c("WT_8WG16\n @37°C July", "pob3_8WG16\n @37°C July", "spt6-YW_8WG16\n @37°C July","spt6-YW-pob3_8WG16\n @37°C July"),
                               perc=1,genes=genes,row_order = order(width(genes),decreasing = T), img = "Leng_8WG16_37C_July.png",lg=F,row.z.transform = T)
  generate.kmeans.change.plots(fllist = fllist,norm = F,prots = c("WT_8WG16\n @30°C Aug", "pob3_8WG16\n @30°C Aug", "spt6-YW_8WG16\n @30°C Aug","spt6-YW-pob3_8WG16\n @30°C Aug"),
                               perc=1,genes=genes,row_order = order(width(genes),decreasing = T), img = "Leng_8WG16_30C_Aug.png",lg=F,row.z.transform = T)
  
  generate.kmeans.change.plots(fllist = fllist,norm = F,prots = c("WT_8WG16\n @37°C Aug", "pob3_8WG16\n @37°C Aug", "spt6-YW_8WG16\n @37°C Aug","spt6-YW-pob3_8WG16\n @37°C Aug"),
                               perc=1,genes=genes,row_order = order(width(genes),decreasing = T), img = "Len_8WG16_37C_Aug.png",lg=F,row.z.transform = T)
  generate.kmeans.change.plots(fllist = fllist,norm = F,prots = c("WT_Flag\n @30°C Aug", "pob3_Flag\n @30°C Aug", "spt6-YW_Flag\n @30°C Aug","spt6-YW-pob3_Flag\n @30°C Aug"),
                               perc=1,genes=genes,row_order = order(width(genes),decreasing = T), img = "Len_Flag_30C_Aug.png",lg=F,row.z.transform = T)
  generate.kmeans.change.plots(fllist = fllist,norm = F,prots = c("WT_Flag\n @37°C Aug", "pob3_Flag\n @37°C Aug", "spt6-YW_Flag\n @37°C Aug","spt6-YW-pob3_Flag\n @37°C Aug"),
                               perc=1,genes=genes,row_order = order(width(genes),decreasing = T), img = "Len_Flag_37C_Aug.png",lg=F,row.z.transform = T)
  generate.kmeans.change.plots(fllist = fllist,norm = F,prots = c("WT_Myc\n @30°C July", "pob3_Myc\n @30°C July", "spt6-YW_Myc\n @30°C July","spt6-YW-pob3_Myc\n @30°C July"),
                               perc=1,genes=genes,row_order = order(width(genes),decreasing = T), img = "Len_Myc_30C_July.png",lg=F,row.z.transform = T)
  generate.kmeans.change.plots(fllist = fllist,norm = F,prots = c("WT_Myc\n @37°C July", "pob3_Myc\n @37°C July", "spt6-YW_Myc\n @37°C July","spt6-YW-pob3_Myc\n @37°C July"),
                               perc=1,genes=genes,row_order = order(width(genes),decreasing = T), img = "Len_Myc_37C_July.png",lg=F,row.z.transform = T)
  
  
  ## row order 
  generate.kmeans.change.plots(fllist = fllist,norm = F,prots = c("WT_8WG16\n @30°C July", "pob3_8WG16\n @30°C July", "spt6-YW_8WG16\n @30°C July","spt6-YW-pob3_8WG16\n @30°C July"),
    perc=1,genes=genes,row_order = order(genes$pol2,decreasing = T), img = "8WG16_30C_July.png",lg=F,row.z.transform = T)
  generate.kmeans.change.plots(fllist = fllist,norm = F,prots = c("WT_8WG16\n @37°C July", "pob3_8WG16\n @37°C July", "spt6-YW_8WG16\n @37°C July","spt6-YW-pob3_8WG16\n @37°C July"),
    perc=1,genes=genes,row_order = order(genes$pol2,decreasing = T), img = "8WG16_37C_July.png",lg=F,row.z.transform = T)
  generate.kmeans.change.plots(fllist = fllist,norm = F,prots = c("WT_8WG16\n @30°C Aug", "pob3_8WG16\n @30°C Aug", "spt6-YW_8WG16\n @30°C Aug","spt6-YW-pob3_8WG16\n @30°C Aug"),
    perc=1,genes=genes,row_order = order(genes$pol2,decreasing = T), img = "8WG16_30C_Aug.png",lg=F,row.z.transform = T)
  generate.kmeans.change.plots(fllist = fllist,norm = F,prots = c("WT_8WG16\n @37°C Aug", "pob3_8WG16\n @37°C Aug", "spt6-YW_8WG16\n @37°C Aug","spt6-YW-pob3_8WG16\n @37°C Aug"),
    perc=1,genes=genes,row_order = order(genes$pol2,decreasing = T), img = "8WG16_37C_Aug.png",lg=F,row.z.transform = T)
  generate.kmeans.change.plots(fllist = fllist,norm = F,prots = c("WT_Flag\n @30°C Aug", "pob3_Flag\n @30°C Aug", "spt6-YW_Flag\n @30°C Aug","spt6-YW-pob3_Flag\n @30°C Aug"),
    perc=1,genes=genes,row_order = order(genes$pol2,decreasing = T), img = "Flag_30C_Aug.png",lg=F,row.z.transform = T)
  generate.kmeans.change.plots(fllist = fllist,norm = F,prots = c("WT_Flag\n @37°C Aug", "pob3_Flag\n @37°C Aug", "spt6-YW_Flag\n @37°C Aug","spt6-YW-pob3_Flag\n @37°C Aug"),
    perc=1,genes=genes,row_order = order(genes$pol2,decreasing = T), img = "Flag_37C_Aug.png",lg=F,row.z.transform = T)
  generate.kmeans.change.plots(fllist = fllist,norm = F,prots = c("WT_Myc\n @30°C July", "pob3_Myc\n @30°C July", "spt6-YW_Myc\n @30°C July","spt6-YW-pob3_Myc\n @30°C July"),
    perc=1,genes=genes,row_order = order(genes$pol2,decreasing = T), img = "Myc_30C_July.png",lg=F,row.z.transform = T)
  generate.kmeans.change.plots(fllist = fllist,norm = F,prots = c("WT_Myc\n @37°C July", "pob3_Myc\n @37°C July", "spt6-YW_Myc\n @37°C July","spt6-YW-pob3_Myc\n @37°C July"),
    perc=1,genes=genes,row_order = order(genes$pol2,decreasing = T), img = "Myc_37C_July.png",lg=F,row.z.transform = T)
  
  ## k-means clusters 
  res <- list()
  res[["8WG16_30C_July"]] <- generate.kmeans.change.plots(k = 3,fllist = fllist,norm = F,prots = c("WT_8WG16\n @30°C July", "pob3_8WG16\n @30°C July", "spt6-YW_8WG16\n @30°C July","spt6-YW-pob3_8WG16\n @30°C July"),
    perc=1,genes=genes,img = "K3_8WG16_30C_July.png",lg=F,row.z.transform = T)
  res[["8WG16_37C_July"]] <- generate.kmeans.change.plots(k=3,fllist = fllist,norm = F,prots = c("WT_8WG16\n @37°C July", "pob3_8WG16\n @37°C July", "spt6-YW_8WG16\n @37°C July","spt6-YW-pob3_8WG16\n @37°C July"),
    perc=1,genes=genes,img = "K3_8WG16_37C_July.png",lg=F,row.z.transform = T)
  res[["8WG16_30C_Aug"]] <- generate.kmeans.change.plots(fllist = fllist,norm = F,prots = c("WT_8WG16\n @30°C Aug", "pob3_8WG16\n @30°C Aug", "spt6-YW_8WG16\n @30°C Aug","spt6-YW-pob3_8WG16\n @30°C Aug"),
    perc=1,genes=genes, img = "K3_8WG16_30C_Aug.png",lg=F,row.z.transform = T)
  res[["8WG16_37C_Aug"]] <- generate.kmeans.change.plots(fllist = fllist,norm = F,prots = c("WT_8WG16\n @37°C Aug", "pob3_8WG16\n @37°C Aug", "spt6-YW_8WG16\n @37°C Aug","spt6-YW-pob3_8WG16\n @37°C Aug"),
    perc=1,genes=genes,img = "K3_8WG16_37C_Aug.png",lg=F,row.z.transform = T)
  res[["Flag_30C_Aug"]] <- generate.kmeans.change.plots(fllist = fllist,norm = F,prots = c("WT_Flag\n @30°C Aug", "pob3_Flag\n @30°C Aug", "spt6-YW_Flag\n @30°C Aug","spt6-YW-pob3_Flag\n @30°C Aug"),
    perc=1,genes=genes,img = "K3_Flag_30C_Aug.png",lg=F,row.z.transform = T)
  res[["Flag_37C_Aug"]] <- generate.kmeans.change.plots(fllist = fllist,norm = F,prots = c("WT_Flag\n @37°C Aug", "pob3_Flag\n @37°C Aug", "spt6-YW_Flag\n @37°C Aug","spt6-YW-pob3_Flag\n @37°C Aug"),
    perc=1,genes=genes, img = "K3_Flag_37C_Aug.png",lg=F,row.z.transform = T)
  res[["_Myc_30C_July"]] <- generate.kmeans.change.plots(fllist = fllist,norm = F,prots = c("WT_Myc\n @30°C July", "pob3_Myc\n @30°C July", "spt6-YW_Myc\n @30°C July","spt6-YW-pob3_Myc\n @30°C July"),
    perc=1,genes=genes, img = "K3_Myc_30C_July.png",lg=F,row.z.transform = T)
  res[["Myc_30C_July"]] <- generate.kmeans.change.plots(fllist = fllist,norm = F,prots = c("WT_Myc\n @37°C July", "pob3_Myc\n @37°C July", "spt6-YW_Myc\n @37°C July","spt6-YW-pob3_Myc\n @37°C July"),
    perc=1,genes=genes,img = "K3_Myc_37C_July.png",lg=F,row.z.transform = T)
  save(res,file="K3_Clusters on Total signalAlong the Genes.RData")
  
  for(i in 1:length(res)){
    myname =paste0("K3_",names(res)[i],"_TotSignal.txt")
    cat(names(res)[i],"\n")
    write.table(res[[i]],file=paste0(myname),sep="\t",quote = F,row.names = F)
  }
  
  
  
  ## RNAPII normalized
  generate.kmeans.change.plots(fllist = fllist,norm = F,
     prots = c("WT_Myc\n @30°C July", "pob3_Myc\n @30°C July", "spt6-YW_Myc\n @30°C July","spt6-YW-pob3_Myc\n @30°C July"),
     control = c("WT_8WG16\n @30°C July", "pob3_8WG16\n @30°C July", "spt6-YW_8WG16\n @30°C July","spt6-YW-pob3_8WG16\n @30°C July"),
     perc=1,genes=genes,row_order = order(genes$pol2,decreasing = T), img = "Myc_30C_July_RNAPIInorm.png",lg=F,row.z.transform = T)
  generate.kmeans.change.plots(fllist = fllist,norm = F,
     prots = c("WT_Myc\n @37°C July", "pob3_Myc\n @37°C July", "spt6-YW_Myc\n @37°C July","spt6-YW-pob3_Myc\n @37°C July"),
     control = c("WT_8WG16\n @37°C July", "pob3_8WG16\n @37°C July", "spt6-YW_8WG16\n @37°C July","spt6-YW-pob3_8WG16\n @37°C July"),
     perc=1,genes=genes,row_order = order(genes$pol2,decreasing = T), img = "Myc_37C_July_RNAPIInorm.png",lg=F,row.z.transform = T)
  generate.kmeans.change.plots(fllist = fllist,norm = F,
     prots = c("WT_Flag\n @30°C Aug", "pob3_Flag\n @30°C Aug", "spt6-YW_Flag\n @30°C Aug","spt6-YW-pob3_Flag\n @30°C Aug"),
     control = c("WT_8WG16\n @30°C Aug", "pob3_8WG16\n @30°C Aug", "spt6-YW_8WG16\n @30°C Aug","spt6-YW-pob3_8WG16\n @30°C Aug"),
     perc=1,genes=genes,row_order = order(genes$pol2,decreasing = T), img = "Flag_30C_Aug_RNAPIInorm.png",lg=F,row.z.transform = T)
  generate.kmeans.change.plots(fllist = fllist,norm = F,
     prots = c("WT_Flag\n @37°C Aug", "pob3_Flag\n @37°C Aug", "spt6-YW_Flag\n @37°C Aug","spt6-YW-pob3_Flag\n @37°C Aug"),
     control = c("WT_8WG16\n @37°C Aug", "pob3_8WG16\n @37°C Aug", "spt6-YW_8WG16\n @37°C Aug","spt6-YW-pob3_8WG16\n @37°C Aug"),
     perc=1,genes=genes,row_order = order(genes$pol2,decreasing = T), img = "Flag_37C_Aug_RNAPIInorm.png",lg=F,row.z.transform = T)
  
  ## WT normalized
  generate.kmeans.change.plots(fllist = fllist,norm = F,
                               prots = c("pob3_8WG16\n @30°C Aug", "spt6-YW_8WG16\n @30°C Aug","spt6-YW-pob3_8WG16\n @30°C Aug"),
                               control = c("WT_8WG16\n @30°C Aug", "WT_8WG16\n @30°C Aug","WT_8WG16\n @30°C Aug"),
                               perc=1,genes=genes,row_order = order(genes$pol2,decreasing = T), img = "8WG16_30C_Aug_WTnorm.png",lg=F,row.z.transform = T)
  generate.kmeans.change.plots(fllist = fllist,norm = F,
                               prots = c("pob3_8WG16\n @37°C Aug", "spt6-YW_8WG16\n @37°C Aug","spt6-YW-pob3_8WG16\n @37°C Aug"),
                               control = c("WT_8WG16\n @37°C Aug", "WT_8WG16\n @37°C Aug","WT_8WG16\n @37°C Aug"),
                               perc=1,genes=genes,row_order = order(genes$pol2,decreasing = T), img = "8WG16_37C_Aug_WTnorm.png",lg=F,row.z.transform = T)
  generate.kmeans.change.plots(fllist = fllist,norm = F,
                               prots = c("pob3_8WG16\n @30°C July", "spt6-YW_8WG16\n @30°C July","spt6-YW-pob3_8WG16\n @30°C July"),
                               control = c("WT_8WG16\n @30°C July", "WT_8WG16\n @30°C July","WT_8WG16\n @30°C July"),
                               perc=1,genes=genes,row_order = order(genes$pol2,decreasing = T), img = "8WG16_30C_July_WTnorm.png",lg=F,row.z.transform = T)
  generate.kmeans.change.plots(fllist = fllist,norm = F,
                               prots = c("pob3_8WG16\n @37°C July", "spt6-YW_8WG16\n @37°C July","spt6-YW-pob3_8WG16\n @37°C July"),
                               control = c("WT_8WG16\n @37°C July", "WT_8WG16\n @37°C July","WT_8WG16\n @37°C July"),
                               perc=1,genes=genes,row_order = order(genes$pol2,decreasing = T), img = "8WG16_37C_July_WTnorm.png",lg=F,row.z.transform = T)
  generate.kmeans.change.plots(fllist = fllist,norm = F,
    prots = c("pob3_Flag\n @30°C Aug", "spt6-YW_Flag\n @30°C Aug","spt6-YW-pob3_Flag\n @30°C Aug"),
    control = c("WT_Flag\n @30°C Aug", "WT_Flag\n @30°C Aug","WT_Flag\n @30°C Aug"),
    perc=1,genes=genes,row_order = order(genes$pol2,decreasing = T), img = "Flag_30C_Aug_WTnorm.png",lg=F,row.z.transform = T)
  generate.kmeans.change.plots(fllist = fllist,norm = F,
    prots = c("pob3_Flag\n @37°C Aug", "spt6-YW_Flag\n @37°C Aug","spt6-YW-pob3_Flag\n @37°C Aug"),
    control = c("WT_Flag\n @37°C Aug", "WT_Flag\n @37°C Aug","WT_Flag\n @37°C Aug"),
    perc=1,genes=genes,row_order = order(genes$pol2,decreasing = T), img = "Flag_37C_Aug_WTnorm.png",lg=F,row.z.transform = T)
  generate.kmeans.change.plots(fllist = fllist,norm = F,
    prots = c("pob3_Myc\n @30°C July", "spt6-YW_Myc\n @30°C July","spt6-YW-pob3_Myc\n @30°C July"),
    control = c("WT_Myc\n @30°C July", "WT_Myc\n @30°C July","WT_Myc\n @30°C July"),
    perc=1,genes=genes,row_order = order(genes$pol2,decreasing = T), img = "Myc_30C_July_WTnorm.png",lg=F,row.z.transform = T)
  generate.kmeans.change.plots(fllist = fllist,norm = F,
    prots = c("pob3_Myc\n @37°C July", "spt6-YW_Myc\n @37°C July","spt6-YW-pob3_Myc\n @37°C July"),
    control = c("WT_Myc\n @37°C July", "WT_Myc\n @37°C July","WT_Myc\n @37°C July"),
    perc=1,genes=genes,row_order = order(genes$pol2,decreasing = T), img = "Myc_37C_July_WTnorm.png",lg=F,row.z.transform = T)
 
  ## WT normalized; K-3
  res <- list()
  res[["FlagChIP_30"]] <- generate.kmeans.change.plots(k=3,fllist = fllist,norm = F,
                               prots = c("pob3_Flag\n @30°C Aug", "spt6-YW_Flag\n @30°C Aug","spt6-YW-pob3_Flag\n @30°C Aug"),
                               control = c("WT_Flag\n @30°C Aug", "WT_Flag\n @30°C Aug","WT_Flag\n @30°C Aug"),
                               perc=1,genes=genes, img = "K3_Flag_30C_Aug_WTnorm.png",lg=F,row.z.transform = T)
  res[["FlagChIP_37"]] <- generate.kmeans.change.plots(k=3,fllist = fllist,norm = F,
                               prots = c("pob3_Flag\n @37°C Aug", "spt6-YW_Flag\n @37°C Aug","spt6-YW-pob3_Flag\n @37°C Aug"),
                               control = c("WT_Flag\n @37°C Aug", "WT_Flag\n @37°C Aug","WT_Flag\n @37°C Aug"),
                               perc=1,genes=genes,img = "K3_Flag_37C_Aug_WTnorm.png",lg=F,row.z.transform = T)
  res[["MycChIP_30"]] <- generate.kmeans.change.plots(k=3,fllist = fllist,norm = F,
                               prots = c("pob3_Myc\n @30°C July", "spt6-YW_Myc\n @30°C July","spt6-YW-pob3_Myc\n @30°C July"),
                               control = c("WT_Myc\n @30°C July", "WT_Myc\n @30°C July","WT_Myc\n @30°C July"),
                               perc=1,genes=genes, img = "K3_Myc_30C_July_WTnorm.png",lg=F,row.z.transform = T)
  res[["MycChIP_37"]] <- generate.kmeans.change.plots(k=3,fllist = fllist,norm = F,
                               prots = c("pob3_Myc\n @37°C July", "spt6-YW_Myc\n @37°C July","spt6-YW-pob3_Myc\n @37°C July"),
                               control = c("WT_Myc\n @37°C July", "WT_Myc\n @37°C July","WT_Myc\n @37°C July"),
                               perc=1,genes=genes,img = "K3_Myc_37C_July_WTnorm.png",lg=F,row.z.transform = T)
  res[["Pol2ChIP_30_Aug"]] <- generate.kmeans.change.plots(k=3,fllist = fllist,norm = F,
                               prots = c("pob3_8WG16\n @30°C Aug", "spt6-YW_8WG16\n @30°C Aug","spt6-YW-pob3_8WG16\n @30°C Aug"),
                               control = c("WT_8WG16\n @30°C Aug", "WT_8WG16\n @30°C Aug","WT_8WG16\n @30°C Aug"),
                               perc=1,genes=genes, img = "K3_8WG16_30C_Aug_WTnorm.png",lg=F,row.z.transform = T)
  res[["Pol2ChIP_37_Aug"]] <- generate.kmeans.change.plots(k=3,fllist = fllist,norm = F,
                               prots = c("pob3_8WG16\n @37°C Aug", "spt6-YW_8WG16\n @37°C Aug","spt6-YW-pob3_8WG16\n @37°C Aug"),
                               control = c("WT_8WG16\n @37°C Aug", "WT_8WG16\n @37°C Aug","WT_8WG16\n @37°C Aug"),
                               perc=1,genes=genes, img = "K3_8WG16_37C_Aug_WTnorm.png",lg=F,row.z.transform = T)
  res[["Pol2ChIP_30_July"]] <- generate.kmeans.change.plots(k=3,fllist = fllist,norm = F,
                               prots = c("pob3_8WG16\n @30°C July", "spt6-YW_8WG16\n @30°C July","spt6-YW-pob3_8WG16\n @30°C July"),
                               control = c("WT_8WG16\n @30°C July", "WT_8WG16\n @30°C July","WT_8WG16\n @30°C July"),
                               perc=1,genes=genes, img = "K3_8WG16_30C_July_WTnorm.png",lg=F,row.z.transform = T)
  res[["Pol2ChIP_37_July"]] <- generate.kmeans.change.plots(k=3,fllist = fllist,norm = F,
                               prots = c("pob3_8WG16\n @37°C July", "spt6-YW_8WG16\n @37°C July","spt6-YW-pob3_8WG16\n @37°C July"),
                               control = c("WT_8WG16\n @37°C July", "WT_8WG16\n @37°C July","WT_8WG16\n @37°C July"),
                               perc=1,genes=genes, img = "K3_8WG16_37C_July_WTnorm.png",lg=F,row.z.transform = T)
  save(res,file="K3_Clusters on Total signalAlong the Genes_WTnormalizezd.RData")
  
  for(i in 1:length(res)){
    myname =paste0("K3_",names(res)[i],"_TotSignal_WTNorm.txt")
    cat(names(res)[i],"\n")
    write.table(res[[i]],file=paste0(myname),sep="\t",quote = F,row.names = F)
  }
  
  ### change with respect to temparature
  generate.kmeans.change.plots(fllist = fllist,norm = F,
    prots = c("WT_Myc\n @37°C July", "pob3_Myc\n @37°C July", "spt6-YW_Myc\n @37°C July","spt6-YW-pob3_Myc\n @37°C July"),
    control = c("WT_Myc\n @30°C July", "pob3_Myc\n @30°C July", "spt6-YW_Myc\n @30°C July","spt6-YW-pob3_Myc\n @30°C July"),
    perc=1,genes=genes,row_order = order(genes$pol2,decreasing = T), img = "Myc_ChIP_normalized_to_30C.png",lg=F,row.z.transform = T)
  generate.kmeans.change.plots(fllist = fllist,norm = F,
    prots = c("WT_8WG16\n @37°C July", "pob3_8WG16\n @37°C July", "spt6-YW_8WG16\n @37°C July","spt6-YW-pob3_8WG16\n @37°C July"),
    control = c("WT_8WG16\n @30°C July", "pob3_8WG16\n @30°C July", "spt6-YW_8WG16\n @30°C July","spt6-YW-pob3_8WG16\n @30°C July"),
    perc=1,genes=genes,row_order = order(genes$pol2,decreasing = T), img = "8WG16_ChIP_normalized_to_30C_July.png",lg=F,row.z.transform = T)
  generate.kmeans.change.plots(fllist = fllist,norm = F,
    prots = c("WT_Flag\n @37°C Aug", "pob3_Flag\n @37°C Aug", "spt6-YW_Flag\n @37°C Aug","spt6-YW-pob3_Flag\n @37°C Aug"),
    control = c("WT_Flag\n @30°C Aug", "pob3_Flag\n @30°C Aug", "spt6-YW_Flag\n @30°C Aug","spt6-YW-pob3_Flag\n @30°C Aug"),
    perc=1,genes=genes,row_order = order(genes$pol2,decreasing = T), img = "Flag_ChIP_normalized_to_30C.png",lg=F,row.z.transform = T)
  generate.kmeans.change.plots(fllist = fllist,norm = F,
    prots = c("WT_8WG16\n @37°C Aug", "pob3_8WG16\n @37°C Aug", "spt6-YW_8WG16\n @37°C Aug","spt6-YW-pob3_8WG16\n @37°C Aug"),
    control = c("WT_8WG16\n @30°C Aug", "pob3_8WG16\n @30°C Aug", "spt6-YW_8WG16\n @30°C Aug","spt6-YW-pob3_8WG16\n @30°C Aug"),
    perc=1,genes=genes,row_order = order(genes$pol2,decreasing = T), img = "8WG16_ChIP_normalized_to_30C_Aug.png",lg=F,row.z.transform = T)
  
  
  ### Cluster on signal normalized to 30C
  generate.kmeans.change.plots(k = 3,fllist = fllist,norm = F,
     prots = c("WT_8WG16\n @37°C Aug", "pob3_8WG16\n @37°C Aug", "spt6-YW_8WG16\n @37°C Aug","spt6-YW-pob3_8WG16\n @37°C Aug"),
     control = c("WT_8WG16\n @30°C Aug", "pob3_8WG16\n @30°C Aug", "spt6-YW_8WG16\n @30°C Aug","spt6-YW-pob3_8WG16\n @30°C Aug"),
     perc=1,genes=genes,img = "K3_8WG16_ChIP_normalized_to_30C_Aug.png",lg=F,row.z.transform = T)
  generate.kmeans.change.plots(k = 3,fllist = fllist,norm = F,
     prots = c("WT_8WG16\n @37°C July", "pob3_8WG16\n @37°C July", "spt6-YW_8WG16\n @37°C July","spt6-YW-pob3_8WG16\n @37°C July"),
     control = c("WT_8WG16\n @30°C July", "pob3_8WG16\n @30°C July", "spt6-YW_8WG16\n @30°C July","spt6-YW-pob3_8WG16\n @30°C July"),
     perc=1,genes=genes,img = "K3_8WG16_ChIP_normalized_to_30C_July.png",lg=F,row.z.transform = T)
  
  
}

### 5' enrichment of the signal :::: Using Ginni formula DF calculations
if(F){
  rm(list=ls())
  setwd("A:/work/WinstonLab/Olga/ChIPSeq")
  get.3prime.fragments <- function(v) {
    v1 <- c()
    for (i in 1:length(v)) { 
      v1[i] <- sum(v[i:length(v)])
    }
    return(rev(v1))
  }
  get.3prime.enrichment.scores <- function(wt,exp,genes,Margin=50,bw=10){
    require(flux)
    exp <- exp[match(genes$tracking_id,rownames(exp)),]
    wt <- wt[match(genes$tracking_id,rownames(wt )),]
    
    m <- data.frame(tracking_id=rownames(exp),TPE=NA,aucscore=NA,tot=NA)
    for (i in 1:nrow(genes)){
      j = ceiling(genes[i,]$width/bw)-1
      v1 <- wt[i,(Margin+1):(Margin+j)]
      wt1 <- get.3prime.fragments(v1)
      
      v2 <- exp[i,(Margin+1):(Margin+j)]
      exp1 <- get.3prime.fragments(v2)
      m$TPE[i] <- round(auc(wt1,exp1)/(sum(v1)*sum(v2)/2),4)
      m$aucscore[i] <- auc(wt1,exp1)
      m$tot[i] <- sum(v1)*sum(v2)
      rm(wt1,exp1,v1,v2,j)
    }
    return(m)
  }
  
  load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenesFull.RData")
  chr.flt <- c("chrI","chrII","chrIII","chrIV","chrIX","chrV","chrVI","chrVII","chrVIII","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
  genes <- genes[genes$chr%in%chr.flt & genes$verification=="Verified",]
  genes <- genes[!is.na(genes$chr),]
  
  ### August
  resl <- list()
  protein<- c("WT_8WG16",  "WT_Flag")
   for( prot in protein){
    f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/MultigeneHeatmap_",prot,".RData",sep="")
    if(file.exists(f)){
      cat(prot,"\n")
      load(f)
      cat(names(l)[1:2],"\n",names(l)[3:4],"\n")
      resl[[ names(l)[3] ]] <- get.3prime.enrichment.scores(wt = l[[1]],exp = l[[3]],genes = genes)
      resl[[ names(l)[4] ]] <- get.3prime.enrichment.scores(wt = l[[2]],exp = l[[4]],genes = genes)
    }
  }
  protein <- c("spt6-YW_8WG16","pob3_8WG16","spt6-YW-pob3_8WG16")
   for( prot in protein){
    f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/MultigeneHeatmap_",prot,".RData",sep="")
    fc = "A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/MultigeneHeatmap_WT_8WG16.RData"
    if(file.exists(f)){
      cat(prot,"\n")
      load(fc)
      ctr <- l
      rm(l)
      load(f)
      cat(names(l)[1:2],"\n",names(l)[3:4],"\n")
      resl[[ names(l)[1] ]] <- get.3prime.enrichment.scores(wt = ctr[[1]],exp = l[[1]],genes = genes)
      resl[[ names(l)[2] ]] <- get.3prime.enrichment.scores(wt = ctr[[2]],exp = l[[2]],genes = genes)
      resl[[ names(l)[3] ]] <- get.3prime.enrichment.scores(wt = ctr[[3]],exp = l[[3]],genes = genes)
      resl[[ names(l)[4] ]] <- get.3prime.enrichment.scores(wt = ctr[[4]],exp = l[[4]],genes = genes)
    }
  }
  protein <- c("pob3_Flag",  "spt6-YW_Flag",  "spt6-YW-pob3_Flag")  
   for( prot in protein){
    f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/MultigeneHeatmap_",prot,".RData",sep="")
    fc = "A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/MultigeneHeatmap_WT_Flag.RData"
    if(file.exists(f)){
      cat(prot,"\n")
      load(fc)
      ctr <- l
      rm(l)
      load(f)
      cat(names(l)[1:2],"\n",names(l)[3:4],"\n")
      resl[[ names(l)[1] ]] <- get.3prime.enrichment.scores(wt = ctr[[1]],exp = l[[1]],genes = genes)
      resl[[ names(l)[2] ]] <- get.3prime.enrichment.scores(wt = ctr[[2]],exp = l[[2]],genes = genes)
      resl[[ names(l)[3] ]] <- get.3prime.enrichment.scores(wt = ctr[[3]],exp = l[[3]],genes = genes)
      resl[[ names(l)[4] ]] <- get.3prime.enrichment.scores(wt = ctr[[4]],exp = l[[4]],genes = genes)
    }
  }
  resq <- resl
  rm(resl,l,ctr,f,fc,prot,protein)
      
  ### July  
  resl <- list()
  protein<- c("WT_8WG16",  "WT_Myc")
  for( prot in protein){
    f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/MultigeneHeatmap_",prot,".RData",sep="")
    if(file.exists(f)){
      cat(prot,"\n")
      load(f)
      cat(names(l)[1:2],"\n",names(l)[3:4],"\n")
      resl[[ names(l)[3] ]] <- get.3prime.enrichment.scores(wt = l[[1]],exp = l[[3]],genes = genes)
      resl[[ names(l)[4] ]] <- get.3prime.enrichment.scores(wt = l[[2]],exp = l[[4]],genes = genes)
    }
  }
  protein <- c("spt6-YW_8WG16","pob3_8WG16","spt6-YW-pob3_8WG16")
  for( prot in protein){
    f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/MultigeneHeatmap_",prot,".RData",sep="")
    fc = "A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/MultigeneHeatmap_WT_8WG16.RData"
    if(file.exists(f)){
      cat(prot,"\n")
      load(fc)
      ctr <- l
      rm(l)
      load(f)
      cat(names(l)[1:2],"\n",names(l)[3:4],"\n")
      resl[[ names(l)[1] ]] <- get.3prime.enrichment.scores(wt = ctr[[1]],exp = l[[1]],genes = genes)
      resl[[ names(l)[2] ]] <- get.3prime.enrichment.scores(wt = ctr[[2]],exp = l[[2]],genes = genes)
      resl[[ names(l)[3] ]] <- get.3prime.enrichment.scores(wt = ctr[[3]],exp = l[[3]],genes = genes)
      resl[[ names(l)[4] ]] <- get.3prime.enrichment.scores(wt = ctr[[4]],exp = l[[4]],genes = genes)
    }
  }
  protein <- c("pob3_Myc",  "spt6-YW_Myc",  "spt6-YW-pob3_Myc")  
  for( prot in protein){
    f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/MultigeneHeatmap_",prot,".RData",sep="")
    fc = "A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/MultigeneHeatmap_WT_Myc.RData"
    if(file.exists(f)){
      cat(prot,"\n")
      load(fc)
      ctr <- l
      rm(l)
      load(f)
      cat(names(l)[1:2],"\n",names(l)[3:4],"\n")
      resl[[ names(l)[1] ]] <- get.3prime.enrichment.scores(wt = ctr[[1]],exp = l[[1]],genes = genes)
      resl[[ names(l)[2] ]] <- get.3prime.enrichment.scores(wt = ctr[[2]],exp = l[[2]],genes = genes)
      resl[[ names(l)[3] ]] <- get.3prime.enrichment.scores(wt = ctr[[3]],exp = l[[3]],genes = genes)
      resl[[ names(l)[4] ]] <- get.3prime.enrichment.scores(wt = ctr[[4]],exp = l[[4]],genes = genes)
    }
  }
  rm(l,ctr,f,fc,prot,protein)
  save(resq,resl,file="TPE_scores_raw.RData")
  
  load(file="TPE_scores_raw.RData")
  pl <- resq[[1]][,1:2]
  names(pl)[2] <- names(resq)[1]
  for(i in 2:length(resq)){
    zz <- resq[[i]][,1:2]
    names(zz)[2] <- names(resq)[i]
    pl <- merge(pl,zz,by="tracking_id")
    rm(zz)
  }
  gl <- resl[[1]][,1:2]
  names(gl)[2] <- names(resl)[1]
  for(i in 2:length(resl)){
    zz <- resl[[i]][,1:2]
    names(zz)[2] <- names(resl)[i]
    gl <- merge(gl,zz,by="tracking_id")
    rm(zz)
  }
  save(pl,gl,file="TPE_scores.RData")
}

### RLE boxplots
if(F){
  rm(list=ls())
  load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenesFull.RData")
  chr.flt <- c("chrI","chrII","chrIII","chrIV","chrIX","chrV","chrVI","chrVII","chrVIII","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
  spiked.chr <- c("chrXVII","chrXVIII","chrXIX")
  sp <- genes[genes$species=="Spom",]
  genes <- genes[genes$chr%in%chr.flt & genes$verification=="Verified",]
  genes <- genes[!is.na(genes$chr),]
  
  getdf <- function(OUTDIR){
    norm.file=paste0(OUTDIR,"/NormalizationFactors.txt")
    load(paste0(OUTDIR,"/ChIPSignalAlongGenes.RData")) 
    rownames(cf) <- cf$tracking_id
    m <- cf[!is.na(cf$verification),]
    m <- m[m$seqnames%in%c(chr.flt,spiked.chr),]
    normF <- read.delim(paste0(norm.file),header=T,stringsAsFactors = F)
    #normF$sample <- gsub("-",".",normF$sample)
    normF$sample%in%names(m)
    m <- m[,normF$sample]
    m <- as.data.frame(apply(m,2,function(x) round(x*1e6/sum(x),2) ))
    normF$sample==names(m)
    m1 <- matrix(NA,nrow = nrow(m),ncol = ncol(m))
    for(i in 1:ncol(m)){
      m1[,i] <- m[,i]*normF$normF[i]
    }
    m1 <- as.data.frame(m1)
    h <- names(m)
    h <- gsub("8WG16ChIP","Pol2",h)
    h <- gsub("MycChIP","Spt16",h)
    h <- gsub("FlagChIP","Spt6",h)
    h <- gsub("spt6.YW","spt6-YW",h)
    h <- gsub("_r1_30",",30°C,R1",h)
    h <- gsub("_r2_30",",30°C,R2",h)
    h <- gsub("_r1_37",",37°C,R1",h)
    h <- gsub("_r2_37",",37°C,R2",h)
    h <- gsub("spt6-YW.pob3","double",h)
    names(m) <- names(m1) <- h
    return(list(norm=m,spikein=m1))
  }
  
  l1 <- getdf(OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr")
  l2 <- getdf(OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr")
  names(l1[[1]]) <- gsub("Pol2","Pol2-Aug",names(l1[[1]]))
  names(l1[[2]]) <- gsub("Pol2","Pol2-Aug",names(l1[[2]]))
  names(l2[[1]]) <- gsub("Pol2","Pol2-July",names(l2[[1]]))
  names(l2[[2]]) <- gsub("Pol2","Pol2-July",names(l2[[2]]))
  cf1 <- cbind(l1[[1]],l2[[1]])
  cf2 <- cbind(l1[[2]],l2[[2]])
  rownames(cf2) <- rownames(cf1)
  rm(l1,l2)
  
  #cf1 <- log2(cf1/apply(cf1,1,function(x) median(x,na.rm=T)))
  #cf2 <- log2(cf2/apply(cf2,1,function(x) median(x,na.rm=T)))
  cf1 <- cf1[!rownames(cf1)%in%genes$tracking_id,]
  cf2 <- cf2[!rownames(cf2)%in%genes$tracking_id,]
  
  
  chip="Pol2-July"
  tmp=30
  pl.list <- list()
  for (chip in c("Pol2-July","Pol2-Aug","Spt6","Spt16")){
    z <- cf1[,grep(chip,names(cf1))]
    z1<- cf2[,grep(chip,names(cf2))]
    for(tmp in c(30,37)){
      exp = (1:ncol(z))[grepl(tmp,names(z),fixed = T)]
      myname= paste0( tmp,"_",chip  )
      cat(myname,"\n")
      cat(names(z)[exp],"\n\n")
      cntr=grep("wt|WT",names(z)[exp])
      zz <- z[,exp]
      zz <- zz[rowSums(zz)>10,]
      zz1 <- z1[,exp]
      zz1 <- subset(zz1,rownames(zz1)%in%rownames(zz))
      
      zz <- log2(zz/apply(zz,1,function(x) median(x,na.rm=T)))
      zz$status <- "autosomes"
      #zz[rownames(zz)%in%sp$tracking_id,]$status <- "heterosomes"
      zz <- reshape::melt(zz,measure.var=names(zz)[1:(ncol(zz)-1)] )
      zz$condition <- 'LibSizeNormalized'
      zz1 <- log2(zz1/apply(zz1,1,function(x) median(x,na.rm=T)))
      zz1$status <- "autosomes"
      #zz1[rownames(zz1)%in%sp$tracking_id,]$status <- "heterosomes"
      zz1 <- reshape::melt(zz1,measure.var=names(zz1)[1:(ncol(zz1)-1)] )
      zz1$condition <- 'SpikeinNormalized'
      zz <- rbind(zz,zz1)
      
      q <- ggplot(zz, aes(x=variable,y=value,fill=variable)) + geom_violin(col='black') +
        geom_boxplot(width=0.1,fill='white', color="black",outlier.shape = NA,outlier.size = 0)+
        geom_hline(yintercept = 0,col="black",size=0.7,linetype="dashed")+
        facet_wrap(~condition+status, nrow = 1,scales = "free_x")+
        ylim(-2,2)+ xlab("") + ylab("relative log occupancy")+
        coord_flip()+
        theme_bw()+theme(legend.position = "none",
          axis.title = element_text(size=14,colour = "black"),
          axis.text = element_text(size=14,color="black"))
      pl.list[[myname]] <- print(q)
      rm(zz,zz1,q)
    }
    rm(z,exp,cntr,myname)
  }
  
 pdf("SpikeinNorm_effectAssessment_RLEBoxes.pdf",width = 8,height = 6)
 for(i in 1:length(pl.list)){
   print(pl.list[[i]])
 }
 dev.off()
  
}

### 5' enrichment of the signal :::: length normalized first 500bp/whole gene
if(F){
  rm(list=ls())
  
  OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr"
  load(paste0(OUTDIR,"/ChIPSignalAlongGenes.RData")) 
  cf1 <- cf;rm(cf)
  load(paste0(OUTDIR,"/ChIPSignalAlongGenes_1st500bp.RData")) 
  cf[,11:ncol(cf)] <- round(cf[,11:ncol(cf)]/cf1[,11:ncol(cf1)],2)
  
  cf$width <- cf1$width
  cf$ml <- 500/cf$width
  cf[,11:58] <- cf[,11:58]*cf$ml
  cf$ml <- NULL
  rm(cf1)
  names(cf)[11:58] <- paste0("Aug_",names(cf)[11:58])
  pl <- cf
  rm(cf)
  
  OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr"
  load(paste0(OUTDIR,"/ChIPSignalAlongGenes.RData")) 
  cf1 <- cf;rm(cf)
  load(paste0(OUTDIR,"/ChIPSignalAlongGenes_1st500bp.RData")) 
  cf[,11:ncol(cf)] <- round(cf[,11:ncol(cf)]/cf1[,11:ncol(cf1)],2)
  cf$width <- cf1$width
  cf$ml <- 500/cf$width
  cf[,11:58] <- cf[,11:58]*cf$ml
  cf$ml <- NULL
  rm(cf1)
  names(cf)[11:58] <- paste0("July_",names(cf)[11:58])
  gl <- cf
  rm(cf)
  
  setwd("A:/work/WinstonLab/Olga/ChIPSeq")
  save(pl,gl,file="FPE_scores.RData")
  
}

 # plots and spreadsheets
if(F){
  rm(list=ls())
  setwd("A:/work/WinstonLab/Olga/ChIPSeq")
  test.diff.func <- function(n=pl,set1=c(14,16,30,32),set2=c(6,8,22,24),name="H3K36me2"){
    x <- combn(set1,2)
    y <- expand.grid(set1,set2)
    
    
    null.dist <- c()
    for(i in 1:dim(x)[2]){
      null.dist <- c(null.dist, (n[,x[1,i]]-n[,x[2,i]]))
    }
    null.dist[!is.finite(null.dist)] <- NA
    xbar <- mean(null.dist,na.rm=T)
    mu <- sd(null.dist,na.rm = T)
    
    
    z <- c()
    for(i in 1:dim(y)[1]){
      z <- cbind(z,
                 pnorm((n[,y[i,1]]-n[,y[i,2]]),mean = xbar,sd = mu,lower.tail = T))
    }
    z <- as.data.frame(z)
    z <- as.data.frame(apply(z,1,function(a) max(a,na.rm = T)))  ## max of hte probability
    #z <- as.data.frame(apply(z,1,function(a) {
    #  a <- a[is.finite(a)]
    #  prod(a)^(1/length(a))}   ))  ## geometric mean- joint probability followed by P-value adjustment
    names(z) <- name
    z[,1] <- ifelse(!is.finite(z[,1]),1,z[,1])
    #z[,1] <- p.adjust(z[,1])
    return(z)
  }
  
  if(T){
    load("TPE_scores.RData")
    names(pl)[2:ncol(pl)] <- paste0("Aug_",names(pl)[2:ncol(pl)])
    names(gl)[2:ncol(gl)] <- paste0("July_",names(gl)[2:ncol(gl)])
    
    OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr"
    load(paste0(OUTDIR,"/ChIPSignalAlongGenes_1st500bp.RData")) 
    cf <- cf[cf$species=="Scer",]
    rownames(cf) <- cf$tracking_id
    cf[,11:58] <- apply(cf[,11:58],2,function(x) x*1e6/sum(x))
    cf <- cf[,11:58]
    cf[cf<10] <- NA
    cf <- cf[complete.cases(cf),]
    pl <- subset(pl,pl$tracking_id%in%rownames(cf))
    
    OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr"
    load(paste0(OUTDIR,"/ChIPSignalAlongGenes_1st500bp.RData")) 
    cf <- cf[cf$species=="Scer",]
    rownames(cf) <- cf$tracking_id
    cf[,11:58] <- apply(cf[,11:58],2,function(x) x*1e6/sum(x))
    cf <- cf[,11:58]
    cf[cf<10] <- NA
    cf <- cf[complete.cases(cf),]
    gl <- subset(gl,gl$tracking_id%in%rownames(cf))
    
    pl <- merge(pl,gl,by="tracking_id")
    pl[,2:ncol(pl)] <- round(1/pl[,2:ncol(pl)],3)
    pl <- pl[complete.cases(pl),]
    rm(gl,cf)
  }
  if(F){
    load("FPE_scores.RData")
    pl <- cbind(pl,gl[,11:58])
    chr.flt <- c("chrI","chrII","chrIII","chrIV","chrIX","chrV","chrVI","chrVII","chrVIII","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
    rm(gl)
    pl <- pl[pl$width>=600,]
    pl <- pl[pl$verification=="Verified",]
    pl <- pl[pl$seqnames%in%chr.flt,]
  }
  
  names(pl) <- gsub("_wt_","_WT_",names(pl))
  
  res <- data.frame(tracking_id=pl$tracking_id)
  pl.list<- list()
  pl.list1<- list()
  for (chip in c("8WG16","Myc","Flag")){
    z <- pl[,c(1,grep(chip,names(pl)))]
    for(prof in c("_pob3_","_spt6-YW_","_spt6-YW-pob3_")){
      for(tmp in c(30,37)){
        cntr = grep("WT",names(z))
        #cntr = (1:ncol(z))[grepl("WT",names(z),fixed = T) & grepl(tmp,names(z),fixed = T)]
        exp = (1:ncol(z))[grepl(prof,names(z),fixed = T) & grepl(tmp,names(z),fixed = T)]
        myname= paste0( gsub("_","",prof),"_",tmp,"_",chip  )
        cat(names(z)[cntr],"\n")
        cat(names(z)[exp],"\n\n")
       
        if(length(cntr)>3 & length(exp)>3){
          cat("performing t-test\n")
          z1 <- z[,c(cntr,exp)]
          out = apply(z1,1,function(x) {
             t.test(x[1:4], x[5:8],na.action=na.omit,alternative = "greater")$p.val
          })
          out <- as.data.frame(out)
        }else{
          out = test.diff.func(n = z,set1 = cntr,set2 = exp,name = myname)
        }
        names(out) <- myname
        res <- cbind(res,out)
       
        col <- c(rev(brewer.pal(8,name = "Reds"))[1:length(exp)],
                   rev(brewer.pal(8,name = "Greens"))[1:length(cntr)])
        names(col) <- c( names(z[,exp]),names(z[,cntr]) )
        qq <- ggplot(melt(z[,c(cntr,exp)])) + theme_bw()+
            stat_density(aes(value, fill = variable,y = (..count..)/20), 
                         position = "identity", color = "gray",alpha=0.2)+
            scale_fill_manual(values = col)+geom_vline(xintercept = c(1))+
            ylab("counts")+xlab("5' enrichment score")+ ggtitle(myname)+
            theme(axis.text = element_text(size=14,color="black"),
                  axis.title = element_text(size=16,colour = "black"),
                  plot.title = element_text(size=18,colour = "black",face = "italic",hjust = 0.5,vjust = 0.5),
                  legend.position = "none")+
            xlim(0.75,1.25)
        pl.list[[paste0(myname)]] <- print(qq)
        rm(qq)
        
       qq <- ggplot(melt(z[,c(cntr,exp)]),aes(col= variable,x=value)) + theme_bw()+
          stat_ecdf(geom="step",size=1.5)+
          scale_color_manual(values = col)+geom_vline(xintercept = c(1))+
          ylab("fraction of total genes")+xlab("5' enrichment score, ECDF")+ ggtitle(myname)+
          theme(axis.text = element_text(size=14,color="black"),
                axis.title = element_text(size=16,colour = "black"),
                plot.title = element_text(size=18,colour = "black",face = "italic",hjust = 0.5,vjust = 0.5),
                legend.position = "none")+
          xlim(0.75,1.25)
       pl.list1[[paste0(myname)]] <- print(qq)
       rm(qq)
       
      }
    }
  }
  rm(out,chip,cntr,col,exp,myname,prof,tmp)
  
  
  lay <- rbind(c(1:5),
               c(6:10),
               c(11:15),
               c(16:20))
  pdf(file = "A:/work/WinstonLab/Olga/ChIPSeq/5primeEnrich_DensityPLots.pdf",height = 12,width = 18)
  #pdf(file = "A:/work/WinstonLab/Olga/ChIPSeq/5primeRatioEnrich_DensityPLots.pdf",height = 12,width = 18)
  grid.arrange(grobs = pl.list, layout_matrix = lay)
  grid.arrange(grobs = pl.list1, layout_matrix = lay)
  dev.off()
  
  pdf(file = "A:/work/WinstonLab/Olga/ChIPSeq/5primeEnrich_DensityPLots_GeneCounts.pdf",height = 8,width = 6)
  p <- data.frame(profile=names(res)[2:19], count=apply(res[,2:19],2,function(x) sum(x<0.05)))
  ggplot(p,aes(x=profile,y=count))+geom_bar(stat="identity")+coord_flip()+theme_bw()+
    xlab("")
  dev.off()
  
  load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenesFull.RData")
  chr.flt <- c("chrI","chrII","chrIII","chrIV","chrIX","chrV","chrVI","chrVII","chrVIII","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
  genes <- genes[genes$chr%in%chr.flt & genes$verification=="Verified",]
  genes <- genes[!is.na(genes$chr),]
  res <- merge(res,genes,by="tracking_id",all.x=T)
  library(xlsx)
  write.xlsx(res, file="FPEscore_pvalues.xlsx")
}

## beta-binomial
if(F){
  rm(list=ls())
  load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenesFull.RData")
  chr.flt <- c("chrI","chrII","chrIII","chrIV","chrIX","chrV","chrVI","chrVII","chrVIII","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
  genes <- genes[genes$chr%in%chr.flt & genes$verification=="Verified",]
  genes <- genes[!is.na(genes$chr),]
  
  my_betabibomial_function <- function(z){
    require(gamlss)
    require(dplyr)
    z$H <- z$up/z$tot
    z$H <- ifelse(z$H==1,0.999,z$H)
    m <- MASS::fitdistr(z$H, dbeta,start = list(shape1 = 1, shape2 = 10))
    alpha0 <- m$estimate[1]
    beta0 <- m$estimate[2]
    prior_mu <- alpha0 / (alpha0 + beta0)
    
    career_eb <- z %>%
      mutate(eb_estimate = (up + alpha0) / (tot + alpha0 + beta0)) %>%
      mutate(alpha1 = up + alpha0, beta1 = tot - up + beta0) %>%
      arrange(desc(eb_estimate))
    
    fit <- gamlss(cbind(up, tot - up) ~ log(tot),data = career_eb[1:1000,],
                  family = BB(mu.link = "logit"))
    mu <- fitted(fit, parameter = "mu")
    sigma <- fitted(fit, parameter = "sigma")
    
    career_eb_wAB <- career_eb %>%
      dplyr::select(tracking_id, H,tot, up,original_eb = eb_estimate) %>%
      mutate(mu = mu,
             alpha0 = mu / sigma,
             beta0 = (1 - mu) / sigma,
             alpha1 = alpha0 + up,
             beta1 = beta0 + tot - up,
             new_eb = alpha1 / (alpha1 + beta1))
    return(career_eb_wAB)
    
  }
  
  ### August samples
  OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr"
  load(paste0(OUTDIR,"/ChIPSignalAlongGenes.RData")) 
  rownames(cf) <- cf$tracking_id
  cf <- apply(cf[,11:58],2,function(x) x*1e6/sum(x))
  cf[cf<10] <- NA
  cf <- cf[complete.cases(cf),]
  cf <- cf[rownames(cf)%in%genes$tracking_id,]
  gnames <- rownames(cf)
  OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr"
  load(paste0(OUTDIR,"/ChIPSignalAlongGenes.RData")) 
  rownames(cf) <- cf$tracking_id
  cf <- cf[cf$tracking_id%in%gnames,]
  cf1 <- cf;rm(cf)
  load(paste0(OUTDIR,"/ChIPSignalAlongGenes_1st500bp.RData")) 
  cf <- cf[cf$tracking_id%in%gnames,]
  cf <- cf[,-grep("Inp",names(cf))]
  cf1 <- cf1[,-grep("Inp",names(cf1))]
  
  alpha1 <- c()
  beta1 <- c()
  for(i in 11:42){
    cat(names(cf1)[i],"\n")
    z <- cf1[,c(1:10,i)]
    z$up <- cf[,i]
    names(z)[11] <- "tot"
    z <- z[z$species=="Scer",]
    z <- z[z$width>=600,]
    z <- z[z$verification=="Verified",]
    z <- z[!is.na(z$seqnames),]
    z <- z[,c("tracking_id","tot","up")]
    z <- z[complete.cases(z),]
    z <- z[z$up>=50,]
    z <- my_betabibomial_function(z)
    alpha1 <- rbind(alpha1, cbind(z[,c("tracking_id","alpha1")], prof=names(cf1)[i] ))
    beta1 <-  rbind(beta1,  cbind(z[,c("tracking_id","beta1")], prof=names(cf1)[i] ))
    rm(z)
  } 
  alpha1 <- reshape2::dcast(alpha1,tracking_id~ prof,value.var = "alpha1")
  beta1 <- reshape2::dcast(beta1,tracking_id~ prof,value.var = "beta1")
  save(alpha1,beta1,file="August_betaBinomPriors_5PrimeSignalEnrichment.RData")
  
  #### July samples
  OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr"
  load(paste0(OUTDIR,"/ChIPSignalAlongGenes.RData")) 
  rownames(cf) <- cf$tracking_id
  cf <- apply(cf[,11:58],2,function(x) x*1e6/sum(x))
  cf[cf<10] <- NA
  cf <- cf[complete.cases(cf),]
  cf <- cf[rownames(cf)%in%genes$tracking_id,]
  gnames <- rownames(cf)
  OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr"
  load(paste0(OUTDIR,"/ChIPSignalAlongGenes.RData")) 
  rownames(cf) <- cf$tracking_id
  cf <- cf[cf$tracking_id%in%gnames,]
  cf1 <- cf;rm(cf)
  load(paste0(OUTDIR,"/ChIPSignalAlongGenes_1st500bp.RData")) 
  cf <- cf[cf$tracking_id%in%gnames,]
  cf <- cf[,-grep("Inp",names(cf))]
  cf1 <- cf1[,-grep("Inp",names(cf1))]
  
  alpha2 <- c()
  beta2 <- c()
  for(i in 11:42){
    cat(names(cf1)[i],"\n")
    z <- cf1[,c(1:10,i)]
    z$up <- cf[,i]
    names(z)[11] <- "tot"
    z <- z[z$species=="Scer",]
    z <- z[z$width>=600,]
    z <- z[z$verification=="Verified",]
    z <- z[!is.na(z$seqnames),]
    z <- z[,c("tracking_id","tot","up")]
    z <- z[complete.cases(z),]
    z <- z[z$up>=50 & z$tot>=50,]
    z <- my_betabibomial_function(z)
    alpha2 <- rbind(alpha2, cbind(z[,c("tracking_id","alpha1")], prof=names(cf1)[i] ))
    beta2 <-  rbind(beta2,  cbind(z[,c("tracking_id","beta1")], prof=names(cf1)[i] ))
    rm(z)
  } 
  alpha2 <- reshape2::dcast(alpha2,tracking_id~ prof,value.var = "alpha1")
  beta2 <- reshape2::dcast(beta2,tracking_id~ prof,value.var = "beta1")
  save(alpha2,beta2,file="July_betaBinomPriors_5PrimeSignalEnrichment.RData")
  
  
  rm(list=ls())
  setwd("A:/work/WinstonLab/Olga/ChIPSeq")
  load(file="August_betaBinomPriors_5PrimeSignalEnrichment.RData")
  names(alpha1)[2:33] <- names(beta1)[2:33] <- paste0("Aug_",names(alpha1)[2:33])
  load(file="July_betaBinomPriors_5PrimeSignalEnrichment.RData")
  names(alpha2)[2:33] <- names(beta2)[2:33] <- paste0("July_",names(alpha2)[2:33])
  
  a1 <- merge(alpha1,alpha2,by="tracking_id")
  b1 <- merge(beta1,beta2,by="tracking_id")
  names(a1) <- gsub("_wt_","_WT_",names(a1))
  names(b1) <- gsub("_wt_","_WT_",names(b1))
  rm(alpha1,alpha2,beta1,beta2)
  
  res <- data.frame(tracking_id=a1$tracking_id)
  for (chip in c("8WG16","Flag","Myc")){
    a <- a1[,c(1,grep(chip,names(a1)))]
    b <- b1[,c(1,grep(chip,names(b1)))]
    for(prof in c("_pob3_","_spt6-YW_","_spt6-YW-pob3_")){
      for(tmp in c(30,37)){
        cntr = grep("WT",names(a))
        cntr = (1:ncol(a))[grepl("WT",names(a),fixed = T) & grepl(tmp,names(a),fixed = T)]
        exp = (1:ncol(a))[grepl(prof,names(a),fixed = T) & grepl(tmp,names(a),fixed = T)]
        myname= paste0( gsub("_","",prof),"_",tmp,"_",chip  )
        cat(names(a)[cntr],"\n")
        cat(names(a)[exp],"\n\n")
        
        ac <- c()
        for(i in cntr){
          ac <- cbind(ac, a[,i])
          ac <- cbind(ac, b[,i])
        }
        for(i in exp){
          ac <- cbind(ac, a[,i])
          ac <- cbind(ac, b[,i])
        }
        ac <- as.data.frame(ac)
        
        pvals <- apply(ac,1,function(i){
          a.samples <- c()
          b.samples <- c()
          for(j in seq(1,length(i),2) ){
            if(j < 2*length(cntr)){
              a.samples <- c(a.samples, rbeta(1e4,shape1 = i[j],shape2 = i[j+1]) ) 
            }else{
              b.samples <- c(b.samples, rbeta(1e4,shape1 = i[j],shape2 = i[j+1]) ) 
            }
          }
          pvl <- sum(a.samples < b.samples)/length(a.samples)
          #cat(length(a.samples),"\t",length(b.samples),"\t",pvl,"\n")
          pvl
        })
        pvals <- as.data.frame(pvals)
        names(pvals) <- myname
        res <- cbind(res,pvals)
        rm(ac,pvals)
      }
    }
    rm(a,b,cntr,exp,myname)
  }
  save(res,file="5PrimeEnrichmentRaio_Beysian_Pvalues.RData")
  
  apply(res[,2:19],2,function(x) sum(x<0.05,na.rm = T))
  
}

### high resolution heatmaps (length sorted, with raw signals and WT normalized signals)
if(F){
  rm(list=ls())
  setwd("A:/work/WinstonLab/Olga/ChIPSeq")
  rowscale_zscore <- function(m){
    a3 <- t(apply(m,1,function(x){
      x[x==0] <- NA
      return((x-mean(x,na.rm = T))/sd(x,na.rm = T))
    }))
    return(a3)
  }
  rowscale_proportions <- function(m){
    a3 <- t(apply(m,1,function(x) round(x*100/sum(x),2)))
    a3[!is.finite(a3)] <- 0
    return(a3)
  }
  get_divsion <- function(m,n){
    m <- m*100+1
    n <- n*100+1
    m <- (m-n)/(m+n)
    return(m)
  }
  heatplotdf <- function(l,genes,g,prot="RNAPII",alignment="TSS",rscale="proportions",
                          normF=c(1,1),max.gene.len=4000,bw=10,Margin=500){
    
     ## spikein normalize the signal
     l[[1]] <- l[[1]]*normF[1]
     l[[2]] <- l[[2]]*normF[2]
     
     ## 
     max.gene.len = round(max.gene.len/bw)
     Margin = round(Margin/bw)
     nrow= length(genes$tracking_id)
     mcol <- max.gene.len+2*Margin
     a1 <- matrix(NA,ncol=mcol,nrow = nrow)  ## TSS
     b1 <- matrix(NA,ncol=mcol,nrow = nrow)  ## CPS
     ncol <- dim(l[[1]])[2]
     for(i in 1:length(genes$tracking_id)){
       v = l[[1]][i,]+l[[2]][i,]
       x <- tail(v[(Margin+1):(Margin+ceiling(genes$width[i]/bw)-1) ],max.gene.len)
       x <- c(v[1:Margin],x,v[(ncol-Margin+1):(ncol-1)])
       b1[i,] <- c(rep(0,(mcol-length(x))),x)
       
       v = l[[1]][i,]+l[[2]][i,]
       x <- head(v[(Margin+1):(Margin+ceiling(genes$width[i]/bw)-1)],max.gene.len)
       x <- c(v[1:Margin],x,v[(ncol-Margin+1):(ncol-1)])
       a1[i,] <- c(x,rep(0,(mcol-length(x))))
     }
     rownames(a1) <- rownames(b1) <-  rownames(l[[1]])
     genes <- subset(genes,genes$tracking_id%in%g$tracking_id)
     a1 <- subset(a1,rownames(a1)%in%genes$tracking_id)
     b1 <- subset(b1,rownames(b1)%in%genes$tracking_id)
     
     if(rscale=="proportions"){
       a1 <- rowscale_proportions(a1)
       b1 <- rowscale_proportions(b1)
     }else if(rscale=="zscore"){
       a1 <- rowscale_zscore(a1)
       b1 <- rowscale_zscore(b1)
     }

     if(alignment=="TSS"){
       a1[!is.finite((a1))] <- NA
       a1 <- as.data.frame(a1)
       a1$tracking_id <- rownames(a1)
       return(a1)
     }else{
       b1[!is.finite((b1))] <- NA
       b1 <- as.data.frame(b1)
       b1$tracking_id <- rownames(b1)
       return(b1)
     }
   }
  lengthsortedhetmaps <- function(m,genes,heat.ylab="verified protein coding genes",
                                  xlab=c("","TSS","1kb","2kb","3kb","4kb",""),
                                  downscale=T,colscale=2,
                                  max.gene.len=4000,bw=10,Margin=500,
                                  outfile="test.png"){
    max.gene.len = round(max.gene.len/bw)
    Margin = round(Margin/bw)
    
    m$variable <- gsub("V","",m$variable)
    m$variable <- as.numeric(m$variable)
    m$value <- ifelse(m$value==0,NA,m$value)
    names(m)[3] <- "scaled_signal"
    m$tracking_id <- factor(m$tracking_id,levels=rev(genes$tracking_id))
    
    h <- Margin+ceiling(genes$width/bw)
    h[h>(max.gene.len+Margin)] <- (max.gene.len+Margin)
    h <- data.frame(ed=h,y=genes$tracking_id)
    h$y <- factor(h$y,levels=rev(genes$tracking_id))
    h$b=Margin
    if(length(xlab)>4){
      breaks= c(1,Margin,100+Margin,200+Margin,300+Margin, (max.gene.len+Margin),(max.gene.len+2*Margin))
    }else{
      breaks= c(1,Margin, (max.gene.len+Margin),(max.gene.len+2*Margin))
    }
    if(downscale==T){
      qq<- quantile(m$scaled_signal,probs=0.9,na.rm=T)
      m$scaled_signal <- ifelse(m$scaled_signal>qq,qq,m$scaled_signal)
    }
    #scale_fill_gdistiller(palette = "YlOrBr",direction = 1,space = "Lab",na.value = "white")+
    
    heat_plot <- ggplot(m) + 
      geom_raster(aes(x=variable, y=tracking_id, fill=scaled_signal),interpolate=F)
  
    if(colscale==2){
      heat_plot <- heat_plot +     
        scale_fill_gradient(low="white",high = "steelblue4",na.value = "white",space = "Lab")
    }else{
      heat_plot <- heat_plot +     
        scale_fill_gradient2(low="steelblue4",mid="white",high = "red4",midpoint = 0,na.value = "white",space = "Lab")
    }
      
    heat_plot <- heat_plot + facet_wrap(~profile,nrow=1)+
      scale_x_continuous(breaks = breaks,labels = xlab,expand = c(0,0))+
      geom_point(data = h,inherit.aes = F,aes(x = ed,y = y),col="black",size=0.05,stroke=0,alpha=0.1, shape=20)+
      geom_point(data = h,inherit.aes = F,aes(x = b,y = y),col="black",size=0.05,stroke=0,alpha=0.1,shape=20)+
      theme(panel.background = element_blank())+
      theme(axis.title.x = element_blank(),
            axis.text.x = element_text(size=16,colour = "black"),
            axis.ticks.y = element_blank(),
            axis.text.y=element_blank(),
            axis.title.y = element_text(color="black",size=18),
            strip.background = element_blank(),
            strip.text = element_text(color="black",size=22))+
      theme(legend.position = "bottom",
            legend.title = element_text(size=16, color="black"),
            legend.text = element_text(size=14, face="plain"),
            legend.margin = margin(0,0,0,0),
            legend.key.width  = unit(1,"cm"),
            legend.box.margin = margin(0,0,0,0))+
      theme(axis.line.x = element_line(size = 0.7, colour = "black"),
            axis.line.y = element_blank())+
      ylab(paste0(nrow(genes)," ",heat.ylab))
    
    png(filename = paste0(outfile),width = 1000,height = 800)
    print(heat_plot)
    dev.off()
    
  }
  
  load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenes.RData")
  genes <- subset(genes,genes$verification=="Verified" & genes$species=="Scer" & genes$chr!="chrM")
  genes <- genes[order(genes$width,decreasing = F),]
  
  
  ################################################################################################33
  ################### heatmap
  ################################################################################################33
  ########### row scaled profiles
  ######
  genotypes <- c("WT",  "pob3","spt6-YW", "spt6-YW-pob3") #,   "spt6-YW_Flag","WT_Flag",  "pob3_Flag",     "spt6-YW-pob3_Flag")
  ## July
  chips <- c("8WG16",  "Myc") #,   "spt6-YW_Flag","WT_Flag",  "pob3_Flag",     "spt6-YW-pob3_Flag")
  for(tmp in c(30,37)){
    for(chip in chips){
      pl <- c()
      for(gt in genotypes){
        f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/MultigeneHeatmap_",gt,"_",chip,".RData",sep="")
        if(file.exists(f)){
          cat(f,"\n")
          load(f)
          if(tmp==30){
            myname = paste0(gt,"_",chip,"\n30°C")
            l <- l[1:2]
          }else{
            myname = paste0(gt,"_",chip,"\n37°C")
            l <- l[3:4]
          }
          m <- heatplotdf(l,genes,g = genes,prot = myname,alignment = "TSS")
          m <- melt(m)
          m$profile <- myname
          pl <- rbind(pl,m)
          rm(m,l,f)
        }else{
          cat("no such file: ", f,"\n")
        }
      }
      if(tmp==30){
        lengthsortedhetmaps(m = pl,genes=genes,outfile=paste0("LgSorted_July_30C_",chip,".png"))
      }else{
        lengthsortedhetmaps(m = pl,genes=genes,outfile=paste0("LgSorted_July_37C_",chip,".png"))
      }
      rm(pl)
    }
  }
  ## Aug
  chips <- c("8WG16",  "Flag") #,   "spt6-YW_Flag","WT_Flag",  "pob3_Flag",     "spt6-YW-pob3_Flag")
  for(tmp in c(30,37)){
    for(chip in chips){
      pl <- c()
      for(gt in genotypes){
        f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/MultigeneHeatmap_",gt,"_",chip,".RData",sep="")
        if(file.exists(f)){
          cat(f,"\n")
          load(f)
          if(tmp==30){
            myname = paste0(gt,"_",chip,"\n30°C")
            l <- l[1:2]
          }else{
            myname = paste0(gt,"_",chip,"\n37°C")
            l <- l[3:4]
          }
          m <- heatplotdf(l,genes,g = genes,prot = myname,alignment = "TSS")
          m <- melt(m)
          m$profile <- myname
          pl <- rbind(pl,m)
          rm(m,l,f)
        }else{
          cat("no such file: ", f,"\n")
        }
      }
      if(tmp==30){
        lengthsortedhetmaps(m = pl,genes=genes,outfile=paste0("Aug_30C_",chip,".png"))
      }else{
        lengthsortedhetmaps(m = pl,genes=genes,outfile=paste0("Aug_37C_",chip,".png"))
      }
      rm(pl)
    }
    
  }
  
  ########### Divided by WT profiles
  ######
  genotypes <- c("pob3","spt6-YW", "spt6-YW-pob3") #,   "spt6-YW_Flag","WT_Flag",  "pob3_Flag",     "spt6-YW-pob3_Flag")
  ## July
  chips <- c("8WG16",  "Myc") #,   "spt6-YW_Flag","WT_Flag",  "pob3_Flag",     "spt6-YW-pob3_Flag")
  for(tmp in c(30,37)){
    for(chip in chips){
      pl <- c()
      f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/MultigeneHeatmap_WT_",chip,".RData",sep="")
      load(f)
      if(tmp==30){
        l <- l[1:2]
      }else{
        l <- l[3:4]
      }
      ct <- heatplotdf(l,genes,g = genes,prot = myname,alignment = "TSS",rscale = "a")
      ct[,1:500] <- ct[,1:500]*100+1
      rm(l)
      for(gt in genotypes){
        f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/MultigeneHeatmap_",gt,"_",chip,".RData",sep="")
        if(file.exists(f)){
          load(f)
          cat(f,"\n")
          if(tmp==30){
            l <- l[1:2]
            myname = paste0(gt,"_",chip,"\n30°C")
          }else{
            l <- l[3:4]
            myname = paste0(gt,"_",chip,"\n37°C")
          }
          load(f)
          m <- heatplotdf(l,genes,g = genes,prot = myname,alignment = "TSS",rscale = "a")
          m[,1:500] <- m[,1:500]*100+1
          m[,1:500] <- round( (m[,1:500]-ct[,1:500])/(m[,1:500]+ct[,1:500]),3)
          m <- melt(m)
          m$profile <- myname
          m$value <- ifelse(!is.finite(m$value),NA,m$value)
          pl <- rbind(pl,m)
          rm(m,l,f)
        }else{
          cat("no such file: ", f,"\n")
        }
      }
      if(tmp==30){
        lengthsortedhetmaps(m = pl,genes=genes,outfile=paste0("July_30C_",chip,"_WTnormalized.png"),downscale = F,colscale = 3)
      }else{
        lengthsortedhetmaps(m = pl,genes=genes,outfile=paste0("July_37C_",chip,"_WTnormalized.png"),downscale = F,colscale = 3)
      }
      rm(pl,ct)
    }
  }
  ## Aug
  chips <- c("8WG16",  "Flag") #,   "spt6-YW_Flag","WT_Flag",  "pob3_Flag",     "spt6-YW-pob3_Flag")
  for(tmp in c(30,37)){
    for(chip in chips){
      pl <- c()
      f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/MultigeneHeatmap_WT_",chip,".RData",sep="")
      load(f)
      if(tmp==30){
        l <- l[1:2]
      }else{
        l <- l[3:4]
      }
      ct <- heatplotdf(l,genes,g = genes,prot = myname,alignment = "TSS",rscale = "a")
      ct[,1:500] <- ct[,1:500]*100+1
      rm(l)
      for(gt in genotypes){
        f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/MultigeneHeatmap_",gt,"_",chip,".RData",sep="")
        if(file.exists(f)){
          load(f)
          cat(f,"\n")
          if(tmp==30){
            l <- l[1:2]
            myname = paste0(gt,"_",chip,"\n30°C")
          }else{
            l <- l[3:4]
            myname = paste0(gt,"_",chip,"\n37°C")
          }
          load(f)
          m <- heatplotdf(l,genes,g = genes,prot = myname,alignment = "TSS",rscale = "a")
          m[,1:500] <- m[,1:500]*100+1
          m[,1:500] <- round( (m[,1:500]-ct[,1:500])/(m[,1:500]+ct[,1:500]),3)
          m <- melt(m)
          m$profile <- myname
          m$value <- ifelse(!is.finite(m$value),NA,m$value)
          pl <- rbind(pl,m)
          rm(m,l,f)
        }else{
          cat("no such file: ", f,"\n")
        }
      }
      if(tmp==30){
        lengthsortedhetmaps(m = pl,genes=genes,outfile=paste0("Aug_30C_",chip,"_WTnormalized.png"),downscale = F,colscale = 3)
      }else{
        lengthsortedhetmaps(m = pl,genes=genes,outfile=paste0("Aug_37C_",chip,"_WTnormalized.png"),downscale = F,colscale = 3)
      }
      rm(pl,ct)
    }
  }
  
  ########### Divided by Pol2 profiles
  ######
  genotypes <- c("WT","pob3","spt6-YW", "spt6-YW-pob3") #,   "spt6-YW_Flag","WT_Flag",  "pob3_Flag",     "spt6-YW-pob3_Flag")
  ## July
  chips <- c("Myc") #,   "spt6-YW_Flag","WT_Flag",  "pob3_Flag",     "spt6-YW-pob3_Flag")
  for(tmp in c(30,37)){
    for(chip in chips){
      pl <- c()
      for(gt in genotypes){
        #####
        f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/MultigeneHeatmap_",gt,"_8WG16.RData",sep="")
        load(f)
        if(tmp==30){
          l <- l[1:2]
        }else{
          l <- l[3:4]
        }
        ct <- heatplotdf(l,genes,g = genes,prot = myname,alignment = "TSS",rscale = "a")
        ct[,1:500] <- ct[,1:500]*100+1
        rm(l)
        
        ####
        f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/MultigeneHeatmap_",gt,"_",chip,".RData",sep="")
        if(file.exists(f)){
          load(f)
          cat(f,"\n")
          if(tmp==30){
            l <- l[1:2]
            myname = paste0(gt,"_",chip,"\n30°C")
          }else{
            l <- l[3:4]
            myname = paste0(gt,"_",chip,"\n37°C")
          }
          load(f)
          m <- heatplotdf(l,genes,g = genes,prot = myname,alignment = "TSS",rscale = "a")
          m[,1:500] <- m[,1:500]*100+1
          m[,1:500] <- round( (m[,1:500]-ct[,1:500])/(m[,1:500]+ct[,1:500]),3)
          m <- melt(m)
          m$profile <- myname
          m$value <- ifelse(!is.finite(m$value),NA,m$value)
          pl <- rbind(pl,m)
          rm(m,l,f)
        }else{
          cat("no such file: ", f,"\n")
        }
      }
      if(tmp==30){
        lengthsortedhetmaps(m = pl,genes=genes,outfile=paste0("July_30C_",chip,"_Pol2normalized.png"),downscale = F,colscale = 3)
      }else{
        lengthsortedhetmaps(m = pl,genes=genes,outfile=paste0("July_37C_",chip,"_Pol2normalized.png"),downscale = F,colscale = 3)
      }
      rm(pl,ct)
    }
  }
  ## Aug
  chips <- c("Flag") #,   "spt6-YW_Flag","WT_Flag",  "pob3_Flag",     "spt6-YW-pob3_Flag")
  for(tmp in c(30,37)){
    for(chip in chips){
      pl <- c()
      for(gt in genotypes){
        #####
        f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/MultigeneHeatmap_",gt,"_8WG16.RData",sep="")
        load(f)
        if(tmp==30){
          l <- l[1:2]
        }else{
          l <- l[3:4]
        }
        ct <- heatplotdf(l,genes,g = genes,prot = myname,alignment = "TSS",rscale = "a")
        ct[,1:500] <- ct[,1:500]*100+1
        rm(l)
        
        ####
        f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/MultigeneHeatmap_",gt,"_",chip,".RData",sep="")
        if(file.exists(f)){
          load(f)
          cat(f,"\n")
          if(tmp==30){
            l <- l[1:2]
            myname = paste0(gt,"_",chip,"\n30°C")
          }else{
            l <- l[3:4]
            myname = paste0(gt,"_",chip,"\n37°C")
          }
          load(f)
          m <- heatplotdf(l,genes,g = genes,prot = myname,alignment = "TSS",rscale = "a")
          m[,1:500] <- m[,1:500]*100+1
          m[,1:500] <- round( (m[,1:500]-ct[,1:500])/(m[,1:500]+ct[,1:500]),3)
          m <- melt(m)
          m$profile <- myname
          m$value <- ifelse(!is.finite(m$value),NA,m$value)
          pl <- rbind(pl,m)
          rm(m,l,f)
        }else{
          cat("no such file: ", f,"\n")
        }
      }
      if(tmp==30){
        lengthsortedhetmaps(m = pl,genes=genes,outfile=paste0("Aug_30C_",chip,"_Pol2normalized.png"),downscale = F,colscale = 3)
      }else{
        lengthsortedhetmaps(m = pl,genes=genes,outfile=paste0("Aug_37C_",chip,"_Pol2normalized.png"),downscale = F,colscale = 3)
      }
      rm(pl,ct)
    }
  }
  
  
  
  
  ################################################################################################33
  ################### metagenes 
  ################################################################################################33
  rm(list=ls())
  load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenes.RData")
  pomg<- subset(genes,genes$verification=="Verified" & genes$species=="Spom" & genes$chr!="chrMsp")
  pomg <- pomg[order(pomg$width,decreasing = F),]
  genes <- subset(genes,genes$verification=="Verified" & genes$species=="Scer" & genes$chr!="chrM")
  genes <- genes[order(genes$width,decreasing = F),]
  
  mymetagenefun <- function(l,genes,normF=c(1,1),max.gene.len=4000, bw=10,Margin=500,facet=NULL){
    l[[1]] <- subset(l[[1]],rownames(l[[1]])%in%genes$tracking_id)
    l[[2]] <- subset(l[[2]],rownames(l[[2]])%in%genes$tracking_id)
    
    l[[1]] <- l[[1]]*normF[1]
    l[[2]] <- l[[2]]*normF[2]
    
    l[[1]] <- l[[1]][match(genes$tracking_id,rownames(l[[1]])),]
    l[[2]] <- l[[2]][match(genes$tracking_id,rownames(l[[2]])),]
    
    l[[1]][l[[1]]<0]<- 0
    l[[2]][l[[2]]<0]<- 0
    
    max.gene.len = round(max.gene.len/bw)
    Margin = round(Margin/bw)
    mcol= max.gene.len+2*Margin-1 ## scaled metagene, adjust last bin of the gene when bw is >1
    
    depl1 <- matrix(NA,ncol=mcol,nrow = nrow(genes))
    depl2 <- matrix(NA,ncol=mcol,nrow = nrow(genes))
    ncol=(dim(l[[1]])[2])
    
    if(T){ ## scaled metagene, adjust last bin of the gene when bw is >1 
      for(i in 1:nrow(genes)){
        j = ceiling(genes[i,]$width/bw)-1
        v <- l[[1]][i,(Margin+1):(Margin+j)]
        v <- vector.resizing(v,max.gene.len)
        v <- c(l[[1]][i,1:Margin],v,l[[1]][i,(ncol-Margin+1):(ncol-1)])
        #plot(1:length(v),v)
        depl1[i,] <- v
        v <- l[[2]][i,(Margin+1):(Margin+j)]
        v <- vector.resizing(v,max.gene.len)
        v <- c(l[[2]][i,1:Margin],v,l[[2]][i,(ncol-Margin+1):(ncol-1)])
        depl2[i,] <- v
      }
    }
    rownames(depl1) <- rownames(depl2) <-  rownames(l[[1]])
    #plot(1:499,apply(depl1,2,sum))
    
    avgdepl <- round((depl1+depl2)/2,3)  ## average signal after spikein normalizationnnnnn
    
    if(!is.null(facet)){
      depl1 <-subset(depl1,rownames(depl1)%in%facet$tracking_id)
      depl2 <-subset(depl2,rownames(depl2)%in%facet$tracking_id)
      avgdepl <-subset(avgdepl,rownames(avgdepl)%in%facet$tracking_id)
    }
    
    get.avg <- function(m,nF){
      df_m <- cbind(1:ncol(m),
                    as.data.frame(apply(m,2,mean,na.rm=TRUE)),
                    as.data.frame(apply(m,2,sd,na.rm=TRUE)), 
                    as.data.frame(apply(m,2, function(x) sd(x,na.rm=TRUE)/sqrt( sum(!is.na(x)) ))),
                    as.data.frame(apply(m,2,median,na.rm=TRUE)),
                    as.data.frame(apply(m,2,IQR,na.rm=TRUE)),
                    as.data.frame(apply(m,2,min,na.rm=TRUE)),
                    as.data.frame(apply(m,2,max,na.rm=TRUE))
      )
      names(df_m) <- c("pos","signal", "sd","se","median","iqr","min","max")
      return(df_m)
    }
    pl <- rbind(cbind(get.avg(depl1),condition="Repl-1"),
                cbind(get.avg(depl2),condition="Repl-2"),
                cbind(get.avg(avgdepl),condition="Average")
    )
    return(pl)
  }
  vector.resizing <- function(x,final.len){
    y <- vector()
    len <- length(x)
    y <-spline(1:len,x,n=final.len)$y
    return(y)
  }
  my_relative_signal <- function(l,l2,normF=c(1,1),normF2=c(1,1)){
    ## spikein normalization
    l[[1]] <- l[[1]]*normF[1]
    l[[2]] <- l[[2]]*normF[2]
    
    l2[[1]] <- l2[[1]]*normF2[1]
    l2[[2]] <- l2[[2]]*normF2[2]
    
    for(i in 1:2){
      l[[i]] = (l[[i]]*100)/(l2[[i]]*100)
      l[[i]][!is.finite(l[[i]])] <- 1
      rownames(l[[i]]) <- rownames(l2[[i]])
    }
    return(l)
    
  }
  
  ########## Raw signal
  genotypes <- c("WT",  "pob3","spt6-YW", "spt6-YW-pob3") #,   "spt6-YW_Flag","WT_Flag",  "pob3_Flag",     "spt6-YW-pob3_Flag")
  pl <- c()
  ## July
  res <- read.delim("A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/NormalizationFactors.txt",header=T,stringsAsFactors = F)
  chips <- c("8WG16",  "Myc") #,   "spt6-YW_Flag","WT_Flag",  "pob3_Flag",     "spt6-YW-pob3_Flag")
  for(tmp in c(30,37)){
    for(chip in chips){
      for(gt in genotypes){
        f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/MultigeneHeatmap_",gt,"_",chip,".RData",sep="")
        if(file.exists(f)){
          cat(f,"\n")
          load(f)
          if(tmp==30){
            myname = paste0(gt,"_",chip,"_30°C")
            l <- l[1:2]
            normF = res[res$Factor==paste0(gt,"_",chip),]$normF[1:2]
          }else{
            myname = paste0(gt,"_",chip,"_37°C")
            l <- l[3:4]
            normF = res[res$Factor==paste0(gt,"_",chip),]$normF[3:4]
          }
          cat(normF,"\n")
          lx <- l
          lx[[1]] <- lx[[1]]*1e6/sum(lx[[1]])
          lx[[2]] <- lx[[2]]*1e6/sum(lx[[2]])
          m <- mymetagenefun(lx,genes=genes)
          m$profile <- myname
          m$month="July"
          m$spike="LibSizeNorm"
          pl <- rbind(pl,m)
          m <- mymetagenefun(l,genes=genes,normF=normF)
          m$profile <- myname
          m$month="July"
          m$spike="SpininNorm"
          pl <- rbind(pl,m)
          rm(m,l,f)
        }else{
          cat("no such file: ", f,"\n")
        }
      }
    }
  }
  ## Aug
  res <- read.delim("A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/NormalizationFactors.txt",header=T,stringsAsFactors = F)
  chips <- c("8WG16",  "Flag") #,   "spt6-YW_Flag","WT_Flag",  "pob3_Flag",     "spt6-YW-pob3_Flag")
  for(tmp in c(30,37)){
    for(chip in chips){
      for(gt in genotypes){
        f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/MultigeneHeatmap_",gt,"_",chip,".RData",sep="")
        if(file.exists(f)){
          cat(f,"\n")
          load(f)
          if(tmp==30){
            myname = paste0(gt,"_",chip,"_30°C")
            l <- l[1:2]
            normF = res[res$Factor==paste0(gt,"_",chip),]$normF[1:2]
          }else{
            myname = paste0(gt,"_",chip,"_37°C")
            l <- l[3:4]
            normF = res[res$Factor==paste0(gt,"_",chip),]$normF[3:4]
          }
          cat(normF,"\n")
          lx <- l
          lx[[1]] <- lx[[1]]*1e6/sum(lx[[1]])
          lx[[2]] <- lx[[2]]*1e6/sum(lx[[2]])
          m <- mymetagenefun(lx,genes=genes)
          m$profile <- myname
          m$month="Aug"
          m$spike="LibSizeNorm"
          pl <- rbind(pl,m)
          m <- mymetagenefun(l,genes=genes,normF=normF)
          m$profile <- myname
          m$month="Aug"
          m$spike="SpininNorm"
          pl <- rbind(pl,m)
          rm(m,l,f)
        }else{
          cat("no such file: ", f,"\n")
        }
      }
    }
  }
  
  ###################### WT normalized signal
  gl <- c()
  genotypes <- c("pob3","spt6-YW", "spt6-YW-pob3") #,   "spt6-YW_Flag","WT_Flag",  "pob3_Flag",     "spt6-YW-pob3_Flag")
  ## July
  chips <- c("8WG16",  "Myc") #,   "spt6-YW_Flag","WT_Flag",  "pob3_Flag",     "spt6-YW-pob3_Flag")
  res <- read.delim("A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/NormalizationFactors.txt",header=T,stringsAsFactors = F)
  for(tmp in c(30,37)){
    for(chip in chips){
      f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/MultigeneHeatmap_WT_",chip,".RData",sep="")
      load(f)
      if(tmp==30){
        normF1 = res[res$Factor==paste0("WT_",chip),]$normF[1:2]
        l1 <- l[1:2]
      }else{
        normF1 = res[res$Factor==paste0("WT_",chip),]$normF[3:4]
        l1 <- l[3:4]
      }
      l1x <- l1
      l1x[[1]] <- l1x[[1]]*1e6/sum(l1x[[1]])
      l1x[[2]] <- l1x[[2]]*1e6/sum(l1x[[2]])
      rm(l)
      for(gt in genotypes){
        f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/MultigeneHeatmap_",gt,"_",chip,".RData",sep="")
        if(file.exists(f)){
          cat(f,"\n")
          load(f)
          if(tmp==30){
            l <- l[1:2]
            normF = res[res$Factor==paste0(gt,"_",chip),]$normF[1:2]
            myname = paste0(gt,"_",chip,"_30C")
          }else{
            normF = res[res$Factor==paste0(gt,"_",chip),]$normF[3:4]
            l <- l[3:4]
            myname = paste0(gt,"_",chip,"_37C")
          }
          lx <- l
          lx[[1]] <- lx[[1]]*1e6/sum(lx[[1]])
          lx[[2]] <- lx[[2]]*1e6/sum(lx[[2]])
          m <- mymetagenefun(my_relative_signal(lx,l1x),genes=genes)
          m$profile <- myname
          m$month="July"
          m$spike="LibSizeNorm"
          gl <- rbind(gl,m)
          m <- mymetagenefun(my_relative_signal(l,l1,normF = normF,normF2 = normF1),genes=genes)
          m$profile <- myname
          m$month="July"
          m$spike="SpikeinNorm"
          gl <- rbind(gl,m)
          rm(m,l,f)
        }else{
          cat("no such file: ", f,"\n")
        }
        rm(l,lx,normF,m)
      }
      rm(l1,l1x,normF1)
    }
  }
  ## Aug
  chips <- c("8WG16",  "Flag") #,   "spt6-YW_Flag","WT_Flag",  "pob3_Flag",     "spt6-YW-pob3_Flag")
  res <- read.delim("A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/NormalizationFactors.txt",header=T,stringsAsFactors = F)
  for(tmp in c(30,37)){
    for(chip in chips){
      f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/MultigeneHeatmap_WT_",chip,".RData",sep="")
      load(f)
      if(tmp==30){
        normF1 = res[res$Factor==paste0("WT_",chip),]$normF[1:2]
        l1 <- l[1:2]
      }else{
        normF1 = res[res$Factor==paste0("WT_",chip),]$normF[3:4]
        l1 <- l[3:4]
      }
      l1x <- l1
      l1x[[1]] <- l1x[[1]]*1e6/sum(l1x[[1]])
      l1x[[2]] <- l1x[[2]]*1e6/sum(l1x[[2]])
      rm(l)
      for(gt in genotypes){
        f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/MultigeneHeatmap_",gt,"_",chip,".RData",sep="")
        if(file.exists(f)){
          cat(f,"\n")
          load(f)
          if(tmp==30){
            l <- l[1:2]
            normF = res[res$Factor==paste0(gt,"_",chip),]$normF[1:2]
            myname = paste0(gt,"_",chip,"_30C")
          }else{
            normF = res[res$Factor==paste0(gt,"_",chip),]$normF[3:4]
            l <- l[3:4]
            myname = paste0(gt,"_",chip,"_37C")
          }
          lx <- l
          lx[[1]] <- lx[[1]]*1e6/sum(lx[[1]])
          lx[[2]] <- lx[[2]]*1e6/sum(lx[[2]])
          m <- mymetagenefun(my_relative_signal(lx,l1x),genes=genes)
          m$profile <- myname
          m$month="Aug"
          m$spike="LibSizeNorm"
          gl <- rbind(gl,m)
          m <- mymetagenefun(my_relative_signal(l,l1,normF = normF,normF2 = normF1),genes=genes)
          m$profile <- myname
          m$month="Aug"
          m$spike="SpikeinNorm"
          gl <- rbind(gl,m)
          rm(m,l,f)
        }else{
          cat("no such file: ", f,"\n")
        }
        rm(l,lx,normF,m)
      }
      rm(normF1,l1,l1x)
    }
  }
  
  
  ###################### Pol2 normalized signal
  tl <- c()
  genotypes <- c("WT","pob3","spt6-YW", "spt6-YW-pob3") #,   "spt6-YW_Flag","WT_Flag",  "pob3_Flag",     "spt6-YW-pob3_Flag")
  ## July
  chips <- c("Myc") #,   "spt6-YW_Flag","WT_Flag",  "pob3_Flag",     "spt6-YW-pob3_Flag")
  res <- read.delim("A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/NormalizationFactors.txt",header=T,stringsAsFactors = F)
  for(tmp in c(30,37)){
    for(chip in chips){
      for(gt in genotypes){
        #### 
        f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/MultigeneHeatmap_",gt,"_8WG16.RData",sep="")
        load(f)
        if(tmp==30){
          normF1 = res[res$Factor==paste0(gt,"_8WG16"),]$normF[1:2]
          l1 <- l[1:2]
        }else{
          normF1 = res[res$Factor==paste0(gt,"_8WG16"),]$normF[3:4]
          l1 <- l[3:4]
        }
        l1x <- l1
        l1x[[1]] <- l1x[[1]]*1e6/sum(l1x[[1]])
        l1x[[2]] <- l1x[[2]]*1e6/sum(l1x[[2]])
        rm(l)
        
        
        ####
        f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt16/gr/MultigeneHeatmap_",gt,"_",chip,".RData",sep="")
        if(file.exists(f)){
          cat(f,"\n")
          load(f)
          if(tmp==30){
            l <- l[1:2]
            normF = res[res$Factor==paste0(gt,"_",chip),]$normF[1:2]
            myname = paste0(gt,"_",chip,"_30C")
          }else{
            normF = res[res$Factor==paste0(gt,"_",chip),]$normF[3:4]
            l <- l[3:4]
            myname = paste0(gt,"_",chip,"_37C")
          }
          lx <- l
          lx[[1]] <- lx[[1]]*1e6/sum(lx[[1]])
          lx[[2]] <- lx[[2]]*1e6/sum(lx[[2]])
          m <- mymetagenefun(my_relative_signal(lx,l1x),genes=genes)
          m$profile <- myname
          m$month="July"
          m$spike="LibSizeNorm"
          tl <- rbind(tl,m)
          m <- mymetagenefun(my_relative_signal(l,l1,normF = normF,normF2 = normF1),genes=genes)
          m$profile <- myname
          m$month="July"
          m$spike="SpikeinNorm"
          tl <- rbind(tl,m)
          rm(m,l,f)
        }else{
          cat("no such file: ", f,"\n")
        }
        rm(l,lx,normF,m,l1,l1x,normF1)
      }
    }
  }
  ## Aug
  chips <- c("Flag") #,   "spt6-YW_Flag","WT_Flag",  "pob3_Flag",     "spt6-YW-pob3_Flag")
  res <- read.delim("A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/NormalizationFactors.txt",header=T,stringsAsFactors = F)
  for(tmp in c(30,37)){
    for(chip in chips){
      for(gt in genotypes){
        #### 
        f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/MultigeneHeatmap_",gt,"_8WG16.RData",sep="")
        load(f)
        if(tmp==30){
          normF1 = res[res$Factor==paste0(gt,"_8WG16"),]$normF[1:2]
          l1 <- l[1:2]
        }else{
          normF1 = res[res$Factor==paste0(gt,"_8WG16"),]$normF[3:4]
          l1 <- l[3:4]
        }
        l1x <- l1
        l1x[[1]] <- l1x[[1]]*1e6/sum(l1x[[1]])
        l1x[[2]] <- l1x[[2]]*1e6/sum(l1x[[2]])
        rm(l)
        
        
        ####
        f = paste("A:/work/WinstonLab/Olga/ChIPSeq/Spt6/gr/MultigeneHeatmap_",gt,"_",chip,".RData",sep="")
        if(file.exists(f)){
          cat(f,"\n")
          load(f)
          if(tmp==30){
            l <- l[1:2]
            normF = res[res$Factor==paste0(gt,"_",chip),]$normF[1:2]
            myname = paste0(gt,"_",chip,"_30C")
          }else{
            normF = res[res$Factor==paste0(gt,"_",chip),]$normF[3:4]
            l <- l[3:4]
            myname = paste0(gt,"_",chip,"_37C")
          }
          lx <- l
          lx[[1]] <- lx[[1]]*1e6/sum(lx[[1]])
          lx[[2]] <- lx[[2]]*1e6/sum(lx[[2]])
          m <- mymetagenefun(my_relative_signal(lx,l1x),genes=genes)
          m$profile <- myname
          m$month="Aug"
          m$spike="LibSizeNorm"
          tl <- rbind(tl,m)
          m <- mymetagenefun(my_relative_signal(l,l1,normF = normF,normF2 = normF1),genes=genes)
          m$profile <- myname
          m$month="Aug"
          m$spike="SpikeinNorm"
          tl <- rbind(tl,m)
          rm(m,l,f)
        }else{
          cat("no such file: ", f,"\n")
        }
        rm(l,lx,normF,m,l1,l1x,normF1)
      }
    }
  }
  save(pl,gl,tl,file='ScaledMetagenePLots.RData')
  
  
  
  ################# metagene plots
  rm(list=ls())
  load(file='ScaledMetagenePLots.RData')
  myplot <- function(pl1,titl){
    ggplot(pl1, aes(x=pos,y=signal,ymax=signal+se, ymin=signal-se,col=as.factor(profile)))+
      facet_wrap(~spike,nrow = 1,scales = "free")+
      geom_vline(xintercept = c(50,450),colour="gray50", linetype = "dashed") +
      geom_linerange(col='gray') + geom_line(size=1.2) +
      #scale_color_viridis(option = "viridis",alpha = 0.6) +
      scale_color_manual(values =  col) +
      #scale_color_brewer(palette = "Dark2") +
      scale_y_continuous(expand = c(0,0))+
      scale_x_continuous(breaks =c(1,50,150,250,350,450,500), labels=c("","TSS","","","","CPS",""),expand = c(0,0)) +ylab("normalized signal") + xlab("") +
      theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = NA),axis.line = element_line(size = 0.7, colour = "black"))+
      theme(axis.text.x = element_text(colour = 'black',size=15),
            axis.text.y = element_text(colour = 'black',size=15), 
            axis.title.y = element_text(colour = 'black',size=18),
            strip.text.x = element_text(colour = 'black',size=16),
            strip.background = element_blank())+
      theme(legend.text =element_text(colour = 'black',size=14),
            legend.title = element_blank(),
            legend.position = "bottom",
            plot.title = element_text(size=20,hjust = 0.5,vjust = 0.5))+
      ggtitle(titl)+guides(fill=guide_legend(nrow=2,byrow=TRUE))
    
  }
  
  pl$profile <- gsub("spt6-YW-pob3","double",pl$profile)
  pl$profile <- gsub("8WG16","Pol2",pl$profile)
  pl$profile <- gsub("Myc","Spt16",pl$profile)
  pl$profile <- gsub("Flag","Spt6",pl$profile)
  pl$spike <- gsub("SpininNorm","Spike-in scaled",pl$spike)
  pl$spike <- gsub("LibSizeNorm","Library size scaled",pl$spike)
  pdf("ScaledMetageneplots_RawSignal.pdf",width = 8,height = 5)
  for(chip in c("Pol2","Spt16")){
    for(tmp in c("30°C","37°C")){
      pl1 <- pl[pl$condition=="Average" & pl$month=="July",]
      pl1 <- pl1[grepl(chip,pl1$profile,fixed = T) & grepl(tmp,pl1$profile,fixed = T),]
      pl1$profile <- gsub(paste0("_",chip,"_",tmp),"",pl1$profile)
      print(myplot(pl1,paste0(chip," ChIP (July),",tmp)))
      
      pl1 <- pl[pl$condition!="Average" & pl$month=="July",]
      pl1 <- pl1[grepl(chip,pl1$profile,fixed = T) & grepl(tmp,pl1$profile,fixed = T),]
      pl1$profile <- gsub(paste0("_",chip,"_",tmp),"",pl1$profile)
      pl1$condition <- gsub("Repl-","",pl1$condition)
      pl1$profile <- paste0(pl1$profile,"-",pl1$condition)
      print(myplot(pl1,paste0(chip," ChIP (July),",tmp)))
    }
  }
  for(chip in c("Pol2","Spt6")){
    for(tmp in c("30°C","37°C")){
      pl1 <- pl[pl$condition=="Average" & pl$month=="Aug",]
      pl1 <- pl1[grepl(chip,pl1$profile,fixed = T) & grepl(tmp,pl1$profile,fixed = T),]
      pl1$profile <- gsub(paste0("_",chip,"_",tmp),"",pl1$profile)
      print(myplot(pl1,paste0(chip," ChIP (August),",tmp)))
      
      pl1 <- pl[pl$condition!="Average" & pl$month=="Aug",]
      pl1 <- pl1[grepl(chip,pl1$profile,fixed = T) & grepl(tmp,pl1$profile,fixed = T),]
      pl1$profile <- gsub(paste0("_",chip,"_",tmp),"",pl1$profile)
      pl1$condition <- gsub("Repl-","",pl1$condition)
      pl1$profile <- paste0(pl1$profile,"-",pl1$condition)
      print(myplot(pl1,paste0(chip," ChIP (August),",tmp)))
    }
  }
  dev.off()
  
  gl$profile <- gsub("spt6-YW-pob3","double",gl$profile)
  gl$profile <- gsub("8WG16","Pol2",gl$profile)
  gl$profile <- gsub("Myc","Spt16",gl$profile)
  gl$profile <- gsub("Flag","Spt6",gl$profile)
  gl$profile <- gsub("_30C","_30°C",gl$profile)
  gl$profile <- gsub("_37C","_37°C",gl$profile)
  gl$spike <- gsub("SpikeinNorm","Spike-in scaled",gl$spike)
  gl$spike <- gsub("LibSizeNorm","Library size scaled",gl$spike)
  pdf("ScaledMetageneplots_RawSignal_WTnormalized.pdf",width = 8,height = 5)
  for(chip in c("Pol2","Spt16")){
    for(tmp in c("30°C","37°C")){
      pl1 <- gl[gl$condition=="Average" & gl$month=="July",]
      pl1 <- pl1[grepl(chip,pl1$profile,fixed = T) & grepl(tmp,pl1$profile,fixed = T),]
      pl1$profile <- gsub(paste0("_",chip,"_",tmp),"",pl1$profile)
      print(myplot(pl1,paste0(chip," ChIP (July),",tmp)))
      
      pl1 <- gl[gl$condition!="Average" & gl$month=="July",]
      pl1 <- pl1[grepl(chip,pl1$profile,fixed = T) & grepl(tmp,pl1$profile,fixed = T),]
      pl1$profile <- gsub(paste0("_",chip,"_",tmp),"",pl1$profile)
      pl1$condition <- gsub("Repl-","",pl1$condition)
      pl1$profile <- paste0(pl1$profile,"-",pl1$condition)
      print(myplot(pl1,paste0(chip," ChIP (July),",tmp)))
    }
  }
  for(chip in c("Pol2","Spt6")){
    for(tmp in c("30°C","37°C")){
      pl1 <- gl[gl$condition=="Average" & gl$month=="Aug",]
      pl1 <- pl1[grepl(chip,pl1$profile,fixed = T) & grepl(tmp,pl1$profile,fixed = T),]
      pl1$profile <- gsub(paste0("_",chip,"_",tmp),"",pl1$profile)
      print(myplot(pl1,paste0(chip," ChIP (August),",tmp)))
      
      pl1 <- gl[gl$condition!="Average" & gl$month=="Aug",]
      pl1 <- pl1[grepl(chip,pl1$profile,fixed = T) & grepl(tmp,pl1$profile,fixed = T),]
      pl1$profile <- gsub(paste0("_",chip,"_",tmp),"",pl1$profile)
      pl1$condition <- gsub("Repl-","",pl1$condition)
      pl1$profile <- paste0(pl1$profile,"-",pl1$condition)
      print(myplot(pl1,paste0(chip," ChIP (August),",tmp)))
    }
  }
  dev.off()
  
  
  tl$profile <- gsub("spt6-YW-pob3","double",tl$profile)
  tl$profile <- gsub("8WG16","Pol2",tl$profile)
  tl$profile <- gsub("Myc","Spt16",tl$profile)
  tl$profile <- gsub("Flag","Spt6",tl$profile)
  tl$profile <- gsub("_30C","_30°C",tl$profile)
  tl$profile <- gsub("_37C","_37°C",tl$profile)
  tl$spike <- gsub("SpikeinNorm","Spike-in scaled",tl$spike)
  tl$spike <- gsub("LibSizeNorm","Library size scaled",tl$spike)
  pdf("ScaledMetageneplots_RawSignal_Pol2normalized.pdf",width = 8,height = 5)
  for(chip in c("Spt16")){
    for(tmp in c("30°C","37°C")){
      pl1 <- tl[tl$condition=="Average" & tl$month=="July",]
      pl1 <- pl1[grepl(chip,pl1$profile,fixed = T) & grepl(tmp,pl1$profile,fixed = T),]
      pl1$profile <- gsub(paste0("_",chip,"_",tmp),"",pl1$profile)
      print(myplot(pl1,paste0(chip," ChIP (July),",tmp)))
      
      pl1 <- tl[tl$condition!="Average" & tl$month=="July",]
      pl1 <- pl1[grepl(chip,pl1$profile,fixed = T) & grepl(tmp,pl1$profile,fixed = T),]
      pl1$profile <- gsub(paste0("_",chip,"_",tmp),"",pl1$profile)
      pl1$condition <- gsub("Repl-","",pl1$condition)
      pl1$profile <- paste0(pl1$profile,"-",pl1$condition)
      print(myplot(pl1,paste0(chip," ChIP (July),",tmp)))
    }
  }
  for(chip in c("Spt6")){
    for(tmp in c("30°C","37°C")){
      pl1 <- tl[tl$condition=="Average" & tl$month=="Aug",]
      pl1 <- pl1[grepl(chip,pl1$profile,fixed = T) & grepl(tmp,pl1$profile,fixed = T),]
      pl1$profile <- gsub(paste0("_",chip,"_",tmp),"",pl1$profile)
      print(myplot(pl1,paste0(chip," ChIP (August),",tmp)))
      
      pl1 <- tl[tl$condition!="Average" & tl$month=="Aug",]
      pl1 <- pl1[grepl(chip,pl1$profile,fixed = T) & grepl(tmp,pl1$profile,fixed = T),]
      pl1$profile <- gsub(paste0("_",chip,"_",tmp),"",pl1$profile)
      pl1$condition <- gsub("Repl-","",pl1$condition)
      pl1$profile <- paste0(pl1$profile,"-",pl1$condition)
      print(myplot(pl1,paste0(chip," ChIP (August),",tmp)))
    }
  }
  dev.off()
  
}












