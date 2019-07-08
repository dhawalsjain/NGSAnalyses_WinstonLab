#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (length(args)!=3) {
  stop("Input files not specified. Usage: winstonLab_df.R [outdir] [wd] [scriptdir]", call.=FALSE)
}

require(GenomicRanges)
require(parallel)
require(ShortRead)

if(F){
  OUTDIR="/n/data1/hms/dbmi/park/Dhawal/Project_olga/Combined_MNase_Seq"
  WD="/n/data1/hms/dbmi/park/DATA/Winston/rawdata/Olga_MNase_Mar2019_DJ/bam"  
}
chr.flt <- c("chrI","chrII","chrIII","chrIV","chrIX","chrV","chrVI","chrVII","chrVIII","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
spiked.chr <- c("chrXVII","chrXVIII","chrXIX")

OUTDIR=args[1]
WD=args[2]
SCRIPTDIR=args[3]
source(SCRIPTDIR)
dir.create(paste0(OUTDIR),recursive = T,showWarnings = F)

grdir <- paste0(OUTDIR,"/gr")
dfdir <- paste0(OUTDIR,"/dfs")
visdir <- paste0(OUTDIR,"/vis")

create.GRanges=T
frag_size_distr=T
create_bedgraph=T

## create genomic ranges (cluster)
if(create.GRanges){
  require(GenomicRanges)
  require(parallel)
  require(ShortRead)
  dir.create(paste0(grdir),recursive = T,showWarnings = F)
  setwd(WD)
  files <- Sys.glob("*.bam")
  
  l <- mclapply(files,function(f){
    cat("working on file: ",f,"\n")
    gr <- readGAlignmentPairs(f)
    gr <- grglist(gr)
    gr <- unlist(range(gr))
    f <-  sub(".bam",".gr",f)
    f <- paste0(grdir,"/",f)
    save(gr,file=paste0(f))
    return(f)
  },mc.cores = 6)
  
}

## create fragment size distribution
if(frag_size_distr){
  dir.create(paste0(dfdir),recursive = T,showWarnings = F)
  setwd(grdir)
  files <- Sys.glob("*.gr")
  frag.size.plist <- list()
  res <- c()
  for(f in files){
    cat(f,"\n")
    load(f)
    gr <- subset(gr,seqnames(gr)%in%c(chr.flt,spiked.chr))
    gr <- keepSeqlevels(gr,c(chr.flt,spiked.chr))
    d <- as.data.frame(table(seqnames(gr)))
    res <- rbind(res,
                 data.frame(sample=sub(".gr","",f),
                            autosomes=sum(d[d$Var1%in%chr.flt,]$Freq),
                            heterosomes=sum(d[d$Var1%in%spiked.chr,]$Freq)
                 ))
    f <- sub(".gr","",f)
    f <- sub("_[ATGC]+$","",f)
    z <- density(width(gr))
    z$y <- z$y*z$n
    frag.size.plist[[sub(".gr","",f)]] <- z
  }
  save(frag.size.plist,file=paste0(dfdir,"/fragment.size.distribution.RData"))
  save(res,file=paste0(dfdir,"/proportion.spikein.RData"))
}

## create vis bigWig
if(create_bedgraph){
  library(GenomicRanges)
  library(parallel)
  
  dir.create(paste0(visdir),recursive = T,showWarnings = F)
  setwd(grdir)
  files <- Sys.glob("*.gr")
  
  mclapply(files,function(f){
    load(f)
    cat(f,"\n")
    gr <- subset(gr,seqnames(gr)%in%chr.flt)
    gr <- keepSeqlevels(gr,chr.flt)
    g <- gaussian.smoothing.from.gr.MNase(resize(gr,1,fix="center"),bw = 50)
    #save(g,file=sub(".gr",".grview",f))
    f <- sub(".gr",".bedGraph",f)
    f <- paste0(visdir,"/",f)
    gg <- data.frame(seqnames(g),start(g),end(g),g$score)
    write.table(gg,file=paste0(f),quote = F,sep = "\t",row.names = F,col.names = F)
    return(f)
  },mc.cores = 4)
}





  
  