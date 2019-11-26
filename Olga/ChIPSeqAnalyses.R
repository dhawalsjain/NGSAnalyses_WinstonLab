######################################################################################
### Functions
######################################################################################
library(GenomicRanges)
library(ShortRead)
library(parallel)
library(Rsamtools)

rm(list=ls())

get.unscaled.whole.genebody.matrix <- function (covs, points, Margin =500,bw=10){
  #Multigene heatmaps were produced by binning profiles every 10 bp, sorting genes by length, then aligning all at the 50 TSS.
  #Flanking regions 500 bp upstream and downstream are also shown. Genes longer than 4 kb are allowed to run off the right side of the plot
  
  ## Scaling the signal to 10Mio
  if(!"width"%in%names(points)){
    points$width <- points$end-points$start
  }
  points <- points[order(points$width,decreasing = F),]
  
  ncol=max(points$width)+2*Margin
  max.gene.len = max(points$width)
  
  m <- matrix(NA,nrow=dim(points)[1],ncol=ncol)
  r <- c()
  
  ##### Collecting gene upsream/downstream/gene-body density
  for (i in 1:length(points$chr)) {
    if (points$start[i] > Margin  & (points$end[i] + Margin) < length(covs[[as.character(points$chr[i])]]) & points$chr[i] %in% names(covs)) {
      
      if (points$strand[i]=="+") {
        up<- as.vector( window( covs[[as.character(points$chr[i])]], points$start[i]-Margin, points$start[i]-1 ))
        v <- as.vector( window( covs[[as.character(points$chr[i])]], points$start[i], points$end[i]))
        dn<- as.vector( window( covs[[as.character(points$chr[i])]], points$end[i]+1, points$end[i]+Margin ))
        
      }else {
        up<- rev(as.vector( window( covs[[as.character(points$chr[i])]], points$end[i]+1, points$end[i]+Margin )))
        v <- rev(as.vector( window( covs[[as.character(points$chr[i])]], points$start[i], points$end[i])))
        dn<- rev(as.vector( window( covs[[as.character(points$chr[i])]], points$start[i]-Margin, points$start[i]-1 )))
      }
      
      v <- v[1:max.gene.len]
      l <- c(up,v,dn)
      m[i,1:length(l)] <- l
    }
  }
  
  if(bw>1){
    for(i in seq(1,ncol-1,by=bw)){
      j=(i+bw-1)
      if(j>ncol){
        j=ncol
      }
      r <- cbind(r,rowSums(m[,i:j],na.rm = T)  )
    }
  }else{
    r <- m
    r[is.na(r)] <- 0
  }
  rownames(r) <- points$tracking_id
  return(r)
}

merge_files <- function(...,genes=introns,bw=10){
  ff <- c(...)
  ff <- paste(ff,".gr",sep="")
  chr.flt <- c("chrI","chrII","chrIII","chrIV","chrIX","chrV","chrVI","chrVII","chrVIII","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
  spiked.chr <- c("chrXVII","chrXVIII","chrXIX")
  
  g <- GRanges()
  for(i in 1:length(ff)){
    load(ff[i])
    gr <- subset(gr,seqnames(gr)%in%c(chr.flt,spiked.chr))
    gr <- keepSeqlevels(gr,c(chr.flt,spiked.chr))
    g <- c(g,gr)
  }
  g <- resize(g,200,fix="start")
  g <- trim(g)
  cov <- coverage(g)
  cov <- standardize.cov(cov)
  r <- get.unscaled.whole.genebody.matrix(cov,genes,bw=bw)
  return(r)
}

merge_files_tss <- function(...,genes=genes){
  ff <- c(...)
  ff <- paste(ff,".gr",sep="")
  chr.flt <- c("chrI","chrII","chrIII","chrIV","chrIX","chrV","chrVI","chrVII","chrVIII","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
  spiked.chr <- c("chrXVII","chrXVIII","chrXIX")
  
  g <- GRanges()
  for(i in 1:length(ff)){
    load(ff[i])
    gr <- subset(gr,seqnames(gr)%in%c(chr.flt,spiked.chr))
    gr <- keepSeqlevels(gr,c(chr.flt,spiked.chr))
    g <- c(g,gr)
  }
  g <- resize(g,200,fix="start")
  g <- trim(g)
  cov <- coverage(g)
  cov <- standardize.cov(cov)
  r <- getMatrix_TSS(cov,genes,win=c(-500,4000),cumul.df=F,bin.bw=10)
  return(r)
}

merge_files_tts <- function(...,genes=genes){
  ff <- c(...)
  ff <- paste(ff,".gr",sep="")
  chr.flt <- c("chrI","chrII","chrIII","chrIV","chrIX","chrV","chrVI","chrVII","chrVIII","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
  spiked.chr <- c("chrXVII","chrXVIII","chrXIX")
  
  g <- GRanges()
  for(i in 1:length(ff)){
    load(ff[i])
    gr <- subset(gr,seqnames(gr)%in%c(chr.flt,spiked.chr))
    gr <- keepSeqlevels(gr,c(chr.flt,spiked.chr))
    g <- c(g,gr)
  }
  g <- resize(g,200,fix="start")
  cov <- coverage(g)
  cov <- standardize.cov(cov)
  r <- getMatrix_TTS(cov,genes,win=c(-4000,500),cumul.df=F,bin.bw=10)
  return(r)
}

merge_files_scaleHeat <- function(...,genes=genes){
  ff <- c(...)
  ff <- paste(ff,".gr",sep="")
  chr.flt <- c("chrI","chrII","chrIII","chrIV","chrIX","chrV","chrVI","chrVII","chrVIII","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
  spiked.chr <- c("chrXVII","chrXVIII","chrXIX")
  
  g <- GRanges()
  for(i in 1:length(ff)){
    load(ff[i])
    #gr <- keepSeqlevels(gr,c(chr.flt,spiked.chr))
    gr <- keepSeqlevels(gr,chr.flt)
    g <- c(g,gr)
  }
  g <- resize(g,200,fix="start")
  cov <- coverage(g)
  cov <- standardize.cov(cov)
  r <- get.feature.scaled.matrix(cov,genes,cumul.df=F, Margin =500, bin.bw=10)
  return(r)
}

merge_files_cov <- function(ff,genes=genes){
  chr.flt <- c("chrI","chrII","chrIII","chrIV","chrIX","chrV","chrVI","chrVII","chrVIII","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
  spiked.chr <- c("chrXVII","chrXVIII","chrXIX")
  load(ff)
  r <- get.unscaled.whole.genebody.matrix(cov,genes)
  return(r)
}

standardize.cov <- function(covs){
  total.cov <- sum(as.numeric(unlist(covs)))
  covss <- round(((covs * 1e7) / total.cov), 4)
  return(covss)
}

get.feature.scaled.matrix <- function (covs, points, Margin =200, Gene.scale=2000,cumul.df=T,bin.bw=0){
  require(stats)
  
  m_body <- matrix(nrow=length(points$chr), ncol=Gene.scale)
  m_up <- matrix(nrow=length(points$chr), ncol=Margin)
  m_down <- matrix(nrow=length(points$chr), ncol=Margin)
  Margin = Margin -1
  
  ##### Collecting gene upsream/downstream/gene-body density
  for (i in 1:length(points$chr)) {
    if (points$start[i] > Margin  & (points$end[i] + Margin) < length(covs[[as.character(points$chr[i])]]) & points$chr[i] %in% names(covs)) {
      v <- vector()
      if (points$strand[i]=="+") {
        m_up[i,] <- as.vector( window( covs[[as.character(points$chr[i])]], (points$start[i] - Margin), points$start[i] ) )
        m_down[i,] <- as.vector( window( covs[[as.character(points$chr[i])]], points$end[i], (points$end[i]+Margin) ) )
        v <- as.vector( window( covs[[as.character(points$chr[i])]], points$start[i], points$end[i] ) )
      }else {
        m_up[i,] <- rev(as.vector( window( covs[[as.character(points$chr[i])]], points$end[i], (points$end[i] + Margin) )) )
        m_down[i,] <- rev(as.vector( window( covs[[as.character(points$chr[i])]], (points$start[i]- Margin), points$start[i] )) )
        v <- rev(as.vector( window( covs[[as.character(points$chr[i])]], points$start[i], points$end[i] )))
      }
      m_body[i,] <- vector.resizing(v, Gene.scale) ## interpolating the vector
    }
  }
  
  m <- cbind(m_up,m_body, m_down)
  
  if(cumul.df==T & bin.bw==0){
    df_m <- cbind(1:ncol(m),
                  as.data.frame(apply(m,2,mean,na.rm=TRUE)),
                  as.data.frame(apply(m,2,sd,na.rm=TRUE)), 
                  as.data.frame(apply(m,2, function(x) sd(x,na.rm=TRUE)/sqrt(length(x))))
    )
    names(df_m) <- c("Position","Peak", "SD","SE")
    return(df_m)
  }
  else if(cumul.df==F & bin.bw>0){
    r <- c()
    ncol <- ncol(m)
    for(i in seq(1,ncol-1,by=bin.bw)) r <- cbind(r,rowSums(m[,i:(i+bin.bw-1)],na.rm = T)  )
    rownames(r) <- points$tracking_id
    colnames(r) <- 1:ncol(r)
    return(r)
  }
  else if(bin.bw==0 & cumul.df==F){
    return(m)
  }
  else if(cumul.df==T & bin.bw>0){
    r <- c()
    ncol <- ncol(m)
    for(i in seq(1,ncol-1,by=bin.bw)) r <- cbind(r,rowSums(m[,i:(i+bin.bw-1)],na.rm = T)  )
    rownames(r) <- points$tracking_id
    colnames(r) <- 1:ncol(r)
    df_m <- cbind(1:ncol(r),
                  as.data.frame(apply(r,2,mean,na.rm=TRUE)),
                  as.data.frame(apply(r,2,sd,na.rm=TRUE)), 
                  as.data.frame(apply(r,2, function(x) sd(x,na.rm=TRUE)/sqrt(length(x))))
    )
    names(df_m) <- c("Position","Peak", "SD","SE")
    return(df_m)
  }
  
}

getMatrix_TSS <- function(covs, points, win=c(-200,4000),cumul.df=T,bin.bw=0) {
  m <- matrix(NA,nrow=length(points$chr), ncol=win[2]-win[1]+1)
  
  for (i in 1:length(points$chr)) {
    if (points$TSS[i] > win[2] & (points$TSS[i] + win[2]) < length(covs[[as.character(points$chr[i])]]) & points$chr[i] %in% names(covs)) {
      if (points$strand[i]=="+") {
        m[i,] <- as.vector( window( covs[[as.character(points$chr[i])]], points$TSS[i]+win[1], points$TSS[i]+win[2] ) )
      } else {
        m[i,] <- z <-rev(as.vector( window( covs[[as.character(points$chr[i])]], points$TSS[i] - win[2], points$TSS[i] - win[1] )) )
      }
    }
  }
  rownames(m) <- points$tracking_id
  colnames(m) <- 1:ncol(m)
  
  if(cumul.df==T & bin.bw==0){
    df_m <- cbind(1:ncol(m),
                  as.data.frame(apply(m,2,mean,na.rm=TRUE)),
                  as.data.frame(apply(m,2,sd,na.rm=TRUE)), 
                  as.data.frame(apply(m,2, function(x) sd(x,na.rm=TRUE)/sqrt(length(x))))
    )
    names(df_m) <- c("Position","Peak", "SD","SE")
    return(df_m)
  }
  else if(cumul.df==F & bin.bw>0){
    r <- c()
    ncol <- ncol(m)
    for(i in seq(1,ncol-1,by=bin.bw)) r <- cbind(r,rowSums(m[,i:(i+bin.bw-1)],na.rm = T)  )
    rownames(r) <- points$tracking_id
    colnames(r) <- seq(win[1],win[2]-1,by=bin.bw)+(round(bin.bw/2))
    return(r)
  }
  else if(bin.bw==0 & cumul.df==F){
    return(m)
  }
  else if(cumul.df==T & bin.bw>0){
    r <- c()
    ncol <- ncol(m)
    for(i in seq(1,ncol-1,by=bin.bw)) r <- cbind(r,rowSums(m[,i:(i+bin.bw-1)],na.rm = T)  )
    rownames(r) <- points$tracking_id
    colnames(r) <- seq(win[1],win[2]-1,by=bin.bw)+(round(bin.bw/2))
    df_m <- cbind(1:ncol(r),
                  as.data.frame(apply(r,2,mean,na.rm=TRUE)),
                  as.data.frame(apply(r,2,sd,na.rm=TRUE)), 
                  as.data.frame(apply(r,2, function(x) sd(x,na.rm=TRUE)/sqrt(length(x))))
    )
    names(df_m) <- c("Position","Peak", "SD","SE")
    return(df_m)
  }
}

getMatrix_TTS <- function(covs, points, win=c(-4000,200),cumul.df=T,bin.bw=0) {
  m <- matrix(nrow=length(points$chr), ncol=win[2]-win[1]+1)
  
  for (i in 1:length(points$chr)) {
    if (points$TTS[i] > abs(win[1]) & (points$TTS[i]+abs(win[1])) < length(covs[[as.character(points$chr[i])]]) & points$chr[i] %in% names(covs)) {
      if (points$strand[i]=="+") {
        m[i,] <- as.vector( window( covs[[as.character(points$chr[i])]], points$TTS[i]+win[1], points$TTS[i]+win[2] ) )
      } else {
        m[i,] <- rev(as.vector( window( covs[[as.character(points$chr[i])]], points$TTS[i] - win[2], points$TTS[i] - win[1] )) )
      }
    }
  }
  
  rownames(m) <- points$tracking_id
  colnames(m) <- 1:ncol(m)
  
  if(cumul.df==T & bin.bw==0){
    df_m <- cbind(1:ncol(m),
                  as.data.frame(apply(m,2,mean,na.rm=TRUE)),
                  as.data.frame(apply(m,2,sd,na.rm=TRUE)), 
                  as.data.frame(apply(m,2, function(x) sd(x,na.rm=TRUE)/sqrt(length(x))))
    )
    names(df_m) <- c("Position","Peak", "SD","SE")
    return(df_m)
  }
  else if(cumul.df==F & bin.bw>0){
    r <- c()
    ncol <- ncol(m)
    for(i in seq(1,ncol-1,by=bin.bw)) r <- cbind(r,rowSums(m[,i:(i+bin.bw-1)],na.rm = T)  )
    rownames(r) <- points$tracking_id
    colnames(r) <- seq(win[1],win[2]-1,by=bin.bw)+(round(bin.bw/2))
    return(r)
  }
  else if(bin.bw==0 & cumul.df==F){
    return(m)
  }
  else if(cumul.df==T & bin.bw>0){
    r <- c()
    ncol <- ncol(m)
    for(i in seq(1,ncol-1,by=bin.bw)) r <- cbind(r,rowSums(m[,i:(i+bin.bw-1)],na.rm = T)  )
    rownames(r) <- points$tracking_id
    colnames(r) <- seq(win[1],win[2]-1,by=bin.bw)+(round(bin.bw/2))
    df_m <- cbind(1:ncol(r),
                  as.data.frame(apply(r,2,mean,na.rm=TRUE)),
                  as.data.frame(apply(r,2,sd,na.rm=TRUE)), 
                  as.data.frame(apply(r,2, function(x) sd(x,na.rm=TRUE)/sqrt(length(x))))
    )
    names(df_m) <- c("Position","Peak", "SD","SE")
    return(df_m)
  }
  
}

vector.resizing <- function(x,final.len){
  y <- vector()
  len <- length(x)
  y <-spline(1:len,x,n=final.len)$y
  return(y)
}

reduce.cov.resolution <- function(covs,bw=10){
  gr <-GRanges()
  for(i in names(covs)){
    x <- covs[[i]]
    v <- seq(from=1,to=length(x),by=bw)
    gd <- GRanges(i,IRanges(start=v,end=v),"*",score=x[v])
    gr <- c(gr,gd)
  }
  return(coverage(gr,weight=score(gr)))
}

SPP_peakFinder <- function (chipfile, inputfile, project,peaks=T, tag.density=T,chr.flt=T) {
 
  read.bam.tags <- function (filename, read.tag.names = F, fix.chromosome.names = F){
    #if (!is.element("fastcluster", installed.packages()[, 1])) {
    #  stop("Rsamtools Bioconductor package is now required for BAM file support. Please install")
    #}
    ww <- c("flag", "rname", "pos", "isize", "strand", "mapq",
            "qwidth")
    if (read.tag.names) {
      ww <- c(ww, "qname")
    }
    bam <- Rsamtools::scanBam(filename, param = Rsamtools::ScanBamParam(what = ww,
                                                                        flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE)))[[1]]
    strm <- as.integer(bam$strand == "+")
    if (any(bitwAnd(bam$flag, 1))) {
      rl <- list(tags = tapply(1:length(bam$pos), bam$rname,
                               function(ii) as.numeric(na.omit(strm[ii] * bam$pos[ii] -
                                                                 (1 - strm[ii]) * (bam$pos[ii] + bam$isize[ii])))),
                 flen = tapply(1:length(bam$pos), bam$rname, function(ii) as.numeric(na.omit(abs(bam$isize[ii])))))
    }
    else {
      rl <- list(tags = tapply(1:length(bam$pos), bam$rname,
                               function(ii) bam$pos[ii] * strm[ii] - (1 - strm[ii]) *
                                 (bam$pos[ii] + bam$qwidth[ii])))
    }
    rl <- c(rl, list(quality = tapply(1:length(bam$pos), bam$rname,
                                      function(ii) bam$mapq[ii])))
    if (read.tag.names) {
      rl <- c(rl, list(names = tapply(1:length(bam$pos), bam$rname,
                                      function(ii) bam$qname[ii])))
    }
    if (fix.chromosome.names) {
      names(rl) <- gsub("\\.fa", "", names(rl))
    }
    return(rl)
  }
  environment(read.bam.tags) <- asNamespace('spp')
  
  require(snow)
  library(spp,lib.loc = "/n/data1/hms/dbmi/park/Dhawal/R")
  library(fastcluster,lib.loc = "/n/data1/hms/dbmi/park/Dhawal/R")
  
  cluster <- makeCluster(4)
  chrs <- c("chrI","chrII","chrIII","chrIV","chrIX","chrV","chrVI","chrVII","chrVIII","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
  #chrs <- c("chrI","chrII","chrIII","chrIV","chrIX","chrV","chrVI","chrVII","chrVIII","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI","chrXVII","chrXVIII","chrXIX")
  
  # read in tag alignments
  chip <- read.bam.tags(chipfile);                                       ############## FileName
  input <- read.bam.tags(inputfile);                                     ############## FileName
  
  if(chr.flt==T){
    chip[[1]] <- subset(chip[[1]],names(chip[[1]])%in%chrs)
    chip[[2]] <- subset(chip[[2]],names(chip[[2]])%in%chrs)
    input[[1]] <- subset(input[[1]],names(input[[1]])%in%chrs)
    input[[2]] <- subset(input[[2]],names(input[[2]])%in%chrs)
    
  }
  
  
  # get binding info from cross-correlation profile
  binding.characteristics <- get.binding.characteristics(chip,srange=c(50,500),bin=5,cluster=cluster);
  
  # plot cross-correlation profile
  file1 <- paste(project,"ChIP_crosscorrelation.pdf",sep="")
  pdf(file=file1,width=5,height=5)                                    ############## FileName
  par(mar = c(3.5,3.5,1.0,0.5), mgp = c(2,0.65,0), cex = 0.8);
  plot(binding.characteristics$cross.correlation,type='l',xlab="strand shift",ylab="cross-correlation");
  abline(v=binding.characteristics$peak$x,lty=2,col=2)
  dev.off();
  
  file5 <- paste(project,"_binding.characteristics", sep="")
  save(binding.characteristics, file = file5)         	############## FileName
  
  # select informative tags based on the binding characteristics
  CHIP <- select.informative.tags(chip,binding.characteristics);
  INPUT <- select.informative.tags(input,binding.characteristics);
  
  if(peaks==T){
    #binding detection parameters
    fdr <- 1e-2;
    # the binding.characteristics contains the optimized half-size for binding detection window
    detection.window.halfsize <- binding.characteristics$whs;
    # determine binding positions using wtd method
    bp <- find.binding.positions(signal.data=CHIP,control.data=INPUT,fdr=fdr,whs=detection.window.halfsize, cluster=cluster)
    print(paste("detected",sum(unlist(lapply(bp$npl,function(d) length(d$x)))),"peaks"));
    # output detected binding positions
    file2 <- paste(project,"ChIP-Peaks.txt",sep="")
    output.binding.results(bp,file2)                                    ############## FileName
  }
  
  if(tag.density==T){
    #smothened profile writing step
    smoothed.density <- get.smoothed.tag.density(CHIP,control.tags=INPUT,bandwidth=200,step=25,tag.shift=round(binding.characteristics$peak$x/2),scale.by.dataset.size=T);
    file3 <- paste(project, "_density.wig", sep="")
    file4 <- paste(project, "_smooth-bg", sep="")
    writewig(smoothed.density,file3,file4);        ############## FileName
    
  }
  
  
  #rm(chip, input, smoothed.density, bp, detection.window.halfsize, CHIP, INPUT, binding.characteristics);
  
}


#sbatch -t 0-11:59 -n 6 --mem-per-cpu=4G -p short --wrap "Rscript ChIPSeqAnalyses.R"

######################################################################################
### ChIPSeq paper figures, calculate data frames
######################################################################################

## derive signal regions at 95%
if(F){
  cat("deriving SERs\n")
  sampl.correlation.using.percentile.cutoff <- function(files,cutoff=0.99){
    chr.flt <- c("chrI","chrII","chrIII","chrIV","chrIX","chrV","chrVI","chrVII",
                 "chrVIII","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI","chrXVII","chrXVIII","chrXIX")
    
    ## (1) Identify signal regions in the profile with signal above cutoff
    g <- GRangesList()
    g <- mclapply(files,function(f){
      cat(f,"\n")
      load(f)
      gr <- keepSeqlevels(gr,chr.flt)
      gr <- resize(gr,200,fix="start")
      gr <- coverage(gr)
      gr <- as(gr,"GRanges")
      gr <- gr[gr$score >  quantile(gr$score,cutoff)]
      gr <- reduce(gr)
      return(gr)
    },mc.cores = 8)
    
    g1 <- GRanges()
    for(i in 1:length(g)){
      g1 <- c(g1,g[[i]])
    }
    g <- g1
    rm(g1)
    g <- reduce(g)
    g <- g[width(g)>10]
    
    ## (2) Collect signal along the signal regions
    cf <- c()
    sampl <- c()
    for(f in files){
      cat(f,"\n")
      load(f)
      #gr <- keepSeqlevels(gr,chr.flt)
      gr <- resize(gr,200,fix="start")
      cf <- cbind(cf,countOverlaps(g,gr,ignore.strand=T))
      sampl <- append(sampl, sub(".gr","",f))
    }
    colnames(cf) <- sampl
    cf <- as.data.frame(cf)
    values(g) <- cbind(values(g),data.frame(cf))
    
    return(as.data.frame(g,row.names = NULL))
  }
  
  #files <- Sys.glob("*.gr")
  #cf <- sampl.correlation.using.percentile.cutoff(files,cutoff = 0.9)
  #save(cf,file="CountMatrix_SignalRegions_10p.RData")
  #cf <- sampl.correlation.using.percentile.cutoff(files,cutoff = 0.95)
  #save(cf,file="CountMatrix_SignalRegions_5p.RData")
  #cf <- sampl.correlation.using.percentile.cutoff(files,cutoff = 0.99)
  #save(cf,file="CountMatrix_SignalRegions_1p.RData")
  
  files <- c(Sys.glob("*.gr"), 
             Sys.glob("/n/data1/hms/dbmi/park/DATA/Winston/rawdata/Olga_ChIPSeq_July2019_DJ/analyses/gr/*.gr"))
  cf <- sampl.correlation.using.percentile.cutoff(files,cutoff = 0.95)
  save(cf,file="CountMatrix_SignalRegions_CombinedDatasets_5p.RData")
  
}

## count signal along genes
if(F){
  count.grfilelist.genomicBins <- function(grfiles,genomicBins,filtgr,gwidth=200,fix="start",ignore.strand=F){
    cf <- c()
    for(f in grfiles){
      cat(f,"\n")
      load(f)
      if(length(filtgr)>0){
        cat("filtering\n")
        gr <- gr[!overlapsAny(gr,filtgr,ignore.strand=T)]
      }
      if(!is.null(gwidth)){
        cat("resizing\t",gwidth, "  ",fix,"\n")
        gr <- resize(gr,width=gwidth,fix=fix)
      }
      cf <- cbind(cf,countOverlaps(genomicBins,gr,ignore.strand=ignore.strand))
      rm(gr)
    }
    colnames(cf) <- gsub(".gr","",files)
    cf <- cbind(as.data.frame(genomicBins),cf)
    return(cf)
  }
  
  anno <- read.delim("/n/data1/hms/dbmi/park/Dhawal/Genomes/yeast_Annotations/Scer_fromBurak.gtf",header = F,comment.char = "#") 
  anno <- with(anno,GRanges(anno$V1,IRanges(anno$V4,anno$V5),anno$V7,name=anno$V9,biotype=anno$V2,annotation=anno$V3))
  anno$name <- stringr::str_replace(anno$name,"; transcript.*","")
  anno$name <- stringr::str_replace(anno$name,"gene_id ","")
  anno <- subset(anno,anno$biotype=="rRNA" & seqnames(anno)=="chrXII")
  anno <- range(reduce(anno,ignore.strand=T))
  
  
  files <- Sys.glob("*.gr")
  #load("A:/work/yeast_Annotations/ScSp_RevisedGenes.RData")
  load("/n/data1/hms/dbmi/park/Dhawal/Genomes/yeast_Annotations/ScSp_RevisedGenes.RData")
  genes <- subset(genes,!genes$verification%in%c("transposable_element_gene","Verified|silenced_gene","Dubious",""))
  genes <- as(genes,"GRanges")
  cf <- count.grfilelist.genomicBins(grfiles = files,genomicBins = genes,filtgr = anno,ignore.strand = T)
  save(cf,file="ChIPSignalAlongGenes.RData")
  
  load("/n/data1/hms/dbmi/park/Dhawal/Genomes/yeast_Annotations/ScSp_RevisedGenes.RData")
  genes <-resize(genes,500,"start")
  cf <- count.grfilelist.genomicBins(grfiles = files,genomicBins = genes,filtgr = anno,ignore.strand = T)
  save(cf,file="ChIPSignalAlongGenes_1st500bp.RData")
  
  load("/n/data1/hms/dbmi/park/Dhawal/Genomes/yeast_Annotations/ScSp_RevisedGenes.RData")
  genes <-resize(genes,500,"end")
  cf <- count.grfilelist.genomicBins(grfiles = files,genomicBins = genes,filtgr = anno,ignore.strand = T)
  save(cf,file="ChIPSignalAlongGenes_last500bp.RData")
  
}

## Multigene whole gene body 
if(T){
  files <- Sys.glob("*.gr")
  keys <- read.delim("SampleKey.txt")
  z <- unique(as.character(keys$Input))
  z <- z[z!=""]
  
  load("/n/data1/hms/dbmi/park/Dhawal/Genomes/yeast_Annotations/ScSp_SGD.RevisedGenes.RData")
  #load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenes.RData")
  genes <- subset(genes,genes$verification=="Verified" & genes$species=="Scer" & genes$chr!="chrM")
  genes <- genes[order(genes$width,decreasing = F),]
  keys <- split(keys,keys$Factor)
  
  if(F){
    cat("unscaled multigene heatmap dataframes\n")
    mclapply(1:length(keys),function(i){
      l <- list()
      cat("profile: ",names(keys)[i],"\n")
      depl <- as.character(keys[[i]]$ChIP)
      nondepl <- depl[-grep("_37",depl)] 
      depl <- depl[grep("_37",depl)]
      l[[paste0(nondepl[1])]] <- merge_files(nondepl[1],genes = genes)
      l[[paste0(nondepl[2])]] <- merge_files(nondepl[2],genes = genes)
      l[[paste0(depl[1])]] <- merge_files(depl[1],genes = genes)
      l[[paste0(depl[2])]] <- merge_files(depl[2],genes = genes)
      save(l,file = paste("MultigeneHeatmap_",names(keys)[i],".RData",sep=""))
      return(1)
    },mc.cores = 6)
    
  }
  if(F){
    cat("TSS-Centered multigene heatmap dataframes\n")
    
    mclapply(1:length(keys),function(i){
      l <- list()
      cat("profile: ",names(keys)[i],"\n")
      depl <- as.character(keys[[i]]$ChIP)
      nondepl <- depl[-grep("_37",depl)]
      depl <- depl[grep("_37",depl)]
      l[[paste0(nondepl[1])]] <- merge_files_tss(nondepl[1],genes = genes)
      l[[paste0(nondepl[2])]] <- merge_files_tss(nondepl[2],genes = genes)
      l[[paste0(depl[1])]] <- merge_files_tss(depl[1],genes = genes)
      l[[paste0(depl[2])]] <- merge_files_tss(depl[2],genes = genes)
      save(l,file = paste("MultigeneHeatmap_TSS_",names(keys)[i],".RData",sep=""))
      return(1)
    },mc.cores = 6)
    
  }
  if(F){
    cat("CPS-Centered multigene heatmap dataframes\n")
    
    mclapply(1:length(keys),function(i){
      l <- list()
      cat("profile: ",names(keys)[i],"\n")
      depl <- as.character(keys[[i]]$ChIP)
      nondepl <- depl[-grep("_37",depl)]
      depl <- depl[grep("_37",depl)]
      l[[paste0(nondepl[1])]] <- merge_files_tts(nondepl[1],genes = genes)
      l[[paste0(nondepl[2])]] <- merge_files_tts(nondepl[2],genes = genes)
      l[[paste0(depl[1])]] <- merge_files_tts(depl[1],genes = genes)
      l[[paste0(depl[2])]] <- merge_files_tts(depl[2],genes = genes)
      save(l,file = paste("MultigeneHeatmap_CPS_",names(keys)[i],".RData",sep=""))
      return(1)
    },mc.cores = 6)
    
  }
  
  if(T){
    load("/n/data1/hms/dbmi/park/Dhawal/Genomes/yeast_Annotations/ScSp_SGD.RevisedGenes.RData")
    #load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenes.RData")
    genes <- subset(genes,genes$verification=="Verified" & genes$species=="Spom" & genes$chr!="chrMsp")
    genes <- genes[order(genes$width,decreasing = F),]
    cat("unscaled multigene heatmap dataframes on pombe genes\n")
    mclapply(1:length(keys),function(i){
      l <- list()
      cat("profile: ",names(keys)[i],"\n")
      depl <- as.character(keys[[i]]$ChIP)
      nondepl <- depl[-grep("_37",depl)] 
      depl <- depl[grep("_37",depl)]
      l[[paste0(nondepl[1])]] <- merge_files(nondepl[1],genes = genes)
      l[[paste0(nondepl[2])]] <- merge_files(nondepl[2],genes = genes)
      l[[paste0(depl[1])]] <- merge_files(depl[1],genes = genes)
      l[[paste0(depl[2])]] <- merge_files(depl[2],genes = genes)
      save(l,file = paste("Spombe_MultigeneHeatmap_",names(keys)[i],".RData",sep=""))
      return(1)
    },mc.cores = 6)
  }
  
  
  if(F){
    mclapply(1:length(keys),function(i){
      l <- list()
      cat("profile: ",names(keys)[i],"\n")
      depl <- as.character(keys[[i]]$ChIP)
      nondepl <- depl[-grep("_37",depl)]
      depl <- depl[grep("_37",depl)]
      l[["37"]] <- merge_files_tss(depl[1],depl[2],genes = genes)
      l[["30"]] <- merge_files_tss(nondepl[1],nondepl[2],genes = genes)
      save(l,file = paste("MultigeneHeatmap_TSScomb_",names(keys)[i],".RData",sep=""))
      return(1)
    },mc.cores = 4)
  }
  if(F){
    mclapply(1:length(keys),function(i){
      l <- list()
      cat("profile: ",names(keys)[i],"\n")
      depl <- as.character(keys[[i]]$ChIP)
      nondepl <- depl[-grep("_37",depl)]
      depl <- depl[grep("_37",depl)]
      l[["30"]] <- merge_files_scaleHeat(depl[1],depl[2],genes = genes)
      l[["37"]] <- merge_files_scaleHeat(nondepl[1],nondepl[2],genes = genes)
      save(l,file = paste("MultigeneHeatmap_ScHeat_",names(keys)[i],".RData",sep=""))
      return(1)
    },mc.cores = 4)
    
  }
  if(F){ # introns
    introns <- read.delim(file="A:/work/WinstonLab/Natalia_Reim/October2016/gr/IR_Coveragebased/FeaturesForIRR.bed",header = F)
    names(introns) <- c("chr","start","end","id","tracking_id","strand","anno")
    introns <- introns[introns$anno=="introns",]
    introns$width <- introns$end-introns$start
    introns <- introns[order(introns$width,decreasing = F),]
    introns$tracking_id <- introns$id
    
    
    for(i in 1:length(keys)){
      l <- list()
      cat("profile: ",names(keys)[i],"\n")
      depl <- as.character(keys[[i]]$ChIP)
      nondepl <- depl[-grep("Depl",depl)]
      depl <- depl[grep("Depl",depl)]
      cat(depl,"\n",nondepl,"\n\n")
      l[["Depleted-1"]] <- merge_files(depl[1],genes = introns,bw = 10)
      l[["Depleted-2"]] <- merge_files(depl[2],genes = introns,bw = 10)
      l[["Non-depleted-1"]] <- merge_files(nondepl[1],genes = introns,bw = 10)
      l[["Non-depleted-2"]] <- merge_files(nondepl[2],genes = introns,bw = 10)
      save(l,file = paste("Revised_MultigeneHeatmap_Introns_10bp",names(keys)[i],".RData",sep=""))
    }
    
  }
  
}

### Nucleosome signal change
if(F){
    cat("Counting Signal in intervals for Shift ratios and nucleosomes\n")
    mc.cores=6
    library(GenomicRanges)
    files <- Sys.glob("*.gr")
    
    load("/n/data1/hms/dbmi/park/Dhawal/Genomes/yeast_Annotations/ScSp_RevisedGenes.RData")
    #load("A:/work/yeast_Annotations/ScSp_RevisedGenes.RData")
    genes <- subset(genes,!genes$verification%in%c("transposable_element_gene","Verified|silenced_gene","Dubious",""))
    genes <- genes[genes$width>499,]
    genes <- as(genes,"GRanges")
    tss <- resize(genes,1,"start")
    tts <- resize(genes,1,"end")
    
    g5 <- flank(tss,500,start = F) ## first 500bp
    g3 <- flank(tts,500) ## last 500bp
    n1 <- flank(GenomicRanges::shift(tss,-55),165,start=F) ## first nucleosome
    n2 <- flank(GenomicRanges::shift(tss,110),165,start=F) ## 2nd nucleosome
    n3 <- flank(GenomicRanges::shift(tss,275),165,start=F) ## 2nd nucleosome
    n4 <- flank(GenomicRanges::shift(tss,440),165,start=F) ## 2nd nucleosome
    n5 <- flank(GenomicRanges::shift(tss,605),165,start=F) ## 2nd nucleosome
    n6 <- flank(GenomicRanges::shift(tss,770),165,start=F) ## 2nd nucleosome
    
    m <- c()
    for(f in files){
      cat(f,"\n")
      load(f)
      gr <- resize(gr,200)
      
      d <- data.frame(gene.body=countOverlaps(genes,gr,ignore.strand=T),
                      up.500=countOverlaps(g5,gr,ignore.strand=T),
                      down.500=countOverlaps(g3,gr,ignore.strand=T),
                      Nucl1=countOverlaps(n1,gr,ignore.strand=T),
                      Nucl2=countOverlaps(n2,gr,ignore.strand=T),
                      Nucl3=countOverlaps(n3,gr,ignore.strand=T),
                      Nucl4=countOverlaps(n4,gr,ignore.strand=T),
                      Nucl5=countOverlaps(n5,gr,ignore.strand=T),
                      Nucl6=countOverlaps(n6,gr,ignore.strand=T))
      d <- cbind(as.data.frame(genes),d,profile=sub(".gr","",f))
      names(d)[1] <- "chr"
      m <- rbind(m,d)
      rm(d,gr)
    }
    
    save(m,file="RawCounts_ShiftRatioNucls.RData")   ## Spikein normalization, after input normalization
    cat("\n\n")
}

## CNV
if(F){
  cat("performing CNV analysis\n")
  copy.number.variation.chr <- function(file.exp, file.wt,f.out= "test.png",smooth=1000,chr=c("chrI", "chrII", "chrIII", "chrIV", "chrIX","chrV", "chrVI", "chrVII", "chrVIII","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")){
    require(zoo)
    require(GenomicRanges)
    load(file.exp)
    gr = resize(gr,120);exp= coverage(gr); rm(gr);
    exp= exp[which(names(exp)%in%chr)]
    load(file.wt)
    gr = resize(gr,120);wt=coverage(gr);rm(gr);
    wt= wt[which(names(wt)%in%chr)]
    
    
    exp=reduce.cov.resolution(standardize.cov(exp),bw = 50)*100
    wt=reduce.cov.resolution(standardize.cov(wt),bw = 50)*100
    
    #  options(bitmapType="cairo")
    png(filename = f.out,height = 2500,width = 2500)
    par(mfrow=c(4,4))
    par(mar=c(10,15,6,5))
    for (i in 1:length(wt)) { #
      plot(x=rollmean((1:length(wt[[i]])),smooth), y= rollmean(as.vector((exp[[i]])-(wt[[i]])),smooth),main="",type="l",lwd=2,col=rgb(1,0,0,0.4),ylim=c(-0.3,0.5),axes=FALSE,ann=F);box(lty = 0);
      xlab = round(quantile(1:length(wt[[i]]),probs=c(0,0.33,0.66,1)))
      axis(1, at = xlab,label = paste(round(xlab/1e6,2),"M",sep="") , tck = 0.01,cex.axis = 2.5,col = "black")
      mtext(names(wt)[i],1,5,cex = 2.5)
      axis(2, at = c(-0.3,0,0.3),label = c(-0.3,0,0.3), tck = 0.01,cex.axis = 2.5,col = "black")
      mtext("Scaled signal wrt WT",2,5,cex = 2.5)
      abline(h=0,col="black")
    }
    dev.off()
  }
  
  files <- Sys.glob("*.gr")
  files <- files[grep("Inp",files)]
  file.wt="WT_Inp_r1_30.gr"
  files <- files[!files%in%file.wt]
  
  for(f in files){
    cat(f,"\n")
    copy.number.variation.chr(f,file.wt,sub(".gr",".png",f))
  }
  
}

## SPP derived Wig tracks (input subtracted)
if(F){
  #WD="/n/data1/hms/dbmi/park/DATA/Winston/rawdata/Olga_ChIPSeq_Aug2019_DJ/bam/"
  WD="/n/data1/hms/dbmi/park/DATA/Winston/rawdata/Olga_ChIPSeq_July2019_DJ/bam/"
  keys <- read.delim("SampleKey.txt")
  
  for(i in 1:nrow(keys)){
    cat(paste0(WD,as.character(keys$ChIP[i]),".sort.bam"),"\t")
    cat(paste0(WD,as.character(keys$Input[i]),".sort.bam"),"\t")
    cat(as.character(keys$ChIP[i]),"\n")
    
    if(!file.exists(paste0(WD,as.character(keys$ChIP[i]),".sort.bam"))){
      cat("chip file not found\n")
    }
    if(!file.exists(paste0(WD,as.character(keys$Input[i]),".sort.bam"))){
      cat("input file not found\n")
    }
    
    
    SPP_peakFinder(chipfile =paste0(WD,as.character(keys$ChIP[i]),".sort.bam"),
                   inputfile = paste0(WD,as.character(keys$Input[i]),".sort.bam"),
                   project = as.character(keys$ChIP[i]),peaks = F,tag.density = T,chr.flt = T)
                    
  }
  
}


