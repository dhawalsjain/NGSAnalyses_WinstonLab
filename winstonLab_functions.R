library(psych)
library(ggplot2)
library(viridis)
library(GenomicRanges)
library(GGally)
library(ggpubr)
library(grid)
library(gridExtra)
######################################################
## Gaussian smoothening of the resized bam file
######################################################
gaussian.smoothing.from.gr <- function(gr,bw=50,ignore.strand=F,fix="start",width=1){
  require(GenomicRanges)
  
  gr <- resize(gr,width,fix=fix)
  
  if(ignore.strand==F){
    grp <- subset(gr,strand(gr)=="+")
  }else{
    grp <- gr
  }
  covp <- coverage(grp)
  pos <- lapply(covp,function(x){
    x <- as.vector(x)
    cov <- as.vector(ksmooth(1:length(x), x,kernel = "normal",bandwidth = bw)[[2]])
    cov <- round(cov*10,3)
    return(cov)
  })
  pos <- as(pos, "SimpleRleList")
  grp <- as(pos,"GRanges")
  strand(grp) <- "+"
  
  if(ignore.strand==F){
    grn <- subset(gr,strand(gr)=="-")
    covn <- coverage(grn)
    neg <- lapply(covn,function(x){
      x <- as.vector(x)
      cov <- as.vector(ksmooth(1:length(x), x,kernel = "normal",bandwidth = bw)[[2]])
      cov <- round(cov*10,3)
      return(cov)
    })
    neg <- as(neg, "SimpleRleList")
    grn <- as(neg,"GRanges")
    strand(grn) <- "-"
    gr <- c(grp,grn)
  }else{
    gr <- grp
  }
  return(gr)
}

gaussian.smoothing.from.gr.MNase <- function(gr,bw=50){
  require(GenomicRanges)
  
  g <- GRanges()
  
  gr <- resize(gr,1,fix="center")
  cov <- coverage(gr)
  cov <- standardize.cov(cov)  ## standardize the coverage to 10Mio
  
  g <- lapply(names(cov),function(n){
    x <- as.vector(cov[[n]])
    if(length(x)<1){
      return(GRanges())
    }
    
    x <- as.vector(ksmooth(1:length(x), x,kernel = "normal",bandwidth = bw)[[2]])
    x <- round(x*10,2)
    v <- seq(from=1,to=length(x),by=10)
    gd <- GRanges(n,IRanges(start=v,end=v+9),"+",score=unname(tapply(x, (seq_along(x)-1) %/% 10, max)))
    return(gd)
  })
  
  grs <- GRanges()
  for(i in 1:length(g)){
    grs <- c(grs,g[[i]])
  }
  grs <- subset(grs,grs$score>0)
  return(grs)
}

######################################################
## GenomincRanges to Bed file
######################################################
gr2bed <- function(gr){
  gr <- as.data.frame(gr)
  if(length(grep("score",names(gr)))==0){
    gr$score=255
  }
  if(length(grep("name",names(gr)))==0){
    gr$name="."
  }
  gr <- gr[,c(1,2,3,7,6,5)]
  return(gr)
}


######################################################
## resize vector to a unit length
######################################################
vector.resizing <- function(x,final.len){
  y <- vector()
  len <- length(x)
  y <-spline(1:len,x,n=final.len)$y
  return(y)
}


######################################################
## find peaks from the vector
######################################################
find_peaks <- function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}


######################################################
## multiggplot
######################################################
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  require(ggplot2)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


######################################################
## SPP peak finding from ChIP and Inp bam files
######################################################
SPP_peakFinder <- function (chipfile, inputfile, project,peaks=F, tag.density=T,chrs=chrs) {
  
  require(snow)
  require(Rsamtools)
  library(spp,lib.loc = "/n/data1/hms/dbmi/park/Dhawal/R")
  
  cluster <- makeCluster(1)
  
  # read in tag alignments
  chip <- read.bam.tags(chipfile);                                       ############## FileName
  input <- read.bam.tags(inputfile);                                     ############## FileName
  
  if(length(chrs)>1){
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

######################################################
## standardize coverage to 10M reads
######################################################
standardize.cov <- function(covs,seq.depth=1e7){
  total.cov <- sum(as.numeric(unlist(covs)))
  covss <- round(((covs * seq.depth) / total.cov), 4)
  return(covss)
}

######################################################
## reduce coverage resolution
######################################################
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

######################################################
## copy number variation plots
######################################################
copy.number.variation.chr <- function(file.exp, file.wt,fix="start",gwidth=200,bw=50,f.out= "test.png",
                                      smooth=1000,plot.figure=T,titl="project",
                                      chr=c("chrI", "chrII", "chrIII", "chrIV", "chrIX","chrV", "chrVI", "chrVII", "chrVIII","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")){
  require(zoo)
  require(GenomicRanges)
  load(file.exp)
  #gr <- grs
  gr <- subset(gr,seqnames(gr)%in%chr)
  gr <- keepSeqlevels(gr,chr)
  if(!is.null(gwidth)){
    gr = resize(gr,width = gwidth,fix = fix)
  }
  exp= standardize.cov(coverage(gr))
  rm(gr)
  
  load(file.wt)
  #gr <- grs
  gr <- subset(gr,seqnames(gr)%in%chr)
  gr <- keepSeqlevels(gr,chr)
  if(!is.null(gwidth)){
    gr = resize(gr,width = gwidth,fix = fix)
  }
  wt=standardize.cov(coverage(gr))
  rm(gr)
  
  
  exp=log2(reduce.cov.resolution(standardize.cov(exp),bw = bw)*smooth +1)
  wt=log2(reduce.cov.resolution(standardize.cov(wt),bw = bw)*smooth + 1)
  
  #  options(bitmapType="cairo")
  png(filename = f.out,height = 2000,width = 2000)
  par(mfrow=c(4,4))
  par(mar=c(10,15,8,5))
  for (i in 1:length(wt)) { #
    plot(x=rollmean((1:length(wt[[i]])),smooth), y= rollmean(as.vector((exp[[i]])-(wt[[i]])),smooth),main="",type="l",lwd=2,col=rgb(1,0,0,0.4),ylim=c(-0.2,0.3),axes=FALSE,ann=F);box(lty = 0);
    xlab = round(quantile(1:length(wt[[i]]),probs=c(0,0.33,0.66,1)))
    axis(1, at = xlab,label = paste(round(xlab/1e6,2),"M",sep="") , tck = 0.01,cex.axis = 2.5,col = "black")
    mtext(names(wt)[i],1,5,cex = 2.5)
    axis(2, at = c(-0.2,0,0.2),label = c(-0.2,0,0.2), tck = 0.01,cex.axis = 2.5,col = "black")
    mtext("log(Signal wrt WT)",2,5,cex = 2.5)
    mtext(paste0(titl),side = 3, cex = 2.5)
    abline(h=0,col="black")
  }
  
  
  dev.off()
}



######################################################
## copy number variation plots
######################################################
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

######################################################
## DESeq2 run on count matrix
######################################################
myDeSeq2.funtion <- function(df2,type,pscnt=1,condition,path,image.path){
  require("DESeq2")
  require("RColorBrewer")
  require("gplots")
  require("corrplot")
  
  colData <- data.frame(condition,type)
  row.names(colData) <- colData$type
  dds <- DESeqDataSetFromMatrix(countData = df2,colData = colData,design = ~ condition)
  
  dds <- DESeq(dds)
  rld <- rlog(dds)
  m <- DESeq2::counts(dds,normalized=TRUE)
  m <- m[rowSums(m)>0,]
  m <- m+pscnt
  cormat <- cor(log2(m+pscnt),method="pearson")
  cormat <- round(cormat,2)
  
  
  ## (1)
  ## (2)
  hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  distsRL <- dist(t(assay(rld)))
  mat <- as.matrix(distsRL)
  rownames(mat) <- colnames(mat) <- with(colData(dds),paste(type))
  hc <- hclust(distsRL)
  
  pdf(image.path,width = 12,height = 12)
  
  hmcols <- palette(brewer.pal(n = 11, name = "RdYlBu"))
  heatmap.2(cormat,notecol="black",notecex =1,density.info="none",trace="none",margins =c(14,14),
            cexRow = 1.2,cexCol = 1.2,key.ylab = "",keysize = 0.5,col=hmcols,
            sepwidth = c(0.0,0.0),sepcolor = "black",
            key.title = "",key.xlab = "Correlation Coefficient",main="Pearson correlations",cex=4) 
  hmcols <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
  heatmap.2(cormat,notecol="black",notecex =1,density.info="none",trace="none",margins =c(14,14),cexRow = 1.2,cexCol = 1.2,key.ylab = "",keysize = 0.5,col=hmcols,sepwidth = c(0.0,0.0),sepcolor = "black",
            key.title = "",key.xlab = "Correlation Coefficient",main="Pearson correlations",cex=4) 
  
  par(mfrow=c(1,1))
  heatmap(mat, Rowv=as.dendrogram(hc),symm=TRUE, trace="none",col = rev(hmcol), margin=c(14, 14))
  plotPCA(rld, intgroup=c("condition"))
  ## (3)
  corrplot(cormat, type="lower", order = "hclust",mar = c(2,5,2,2))
  ## (4)
  corrplot(cormat, type="lower", order = "hclust", method="number",mar = c(2,5,2,2))
  dev.off()
  
  save(dds,rld,cormat,mat,hmcol,df2,file = path)
}


######################################################
## RLE transformations
######################################################
rle.transform <- function(m,ps=1){
  m <- t(apply(m, 1,function(x){round(x/median(x),2)}))
  m <- log2(m+ps)
  m[!is.finite(m)] <- 0
  m
}  


######################################################
## asine transforamtion of the coverage list
######################################################
arcsine.transform <- function(coverage) {
  require(IRanges)
  
  atransform <- function(coverage,total.cov) {
    asin(sqrt(coverage/total.cov))
  }
  
  total.cov <- sum(as.numeric(unlist(lapply(coverage,sum))))
  cov.norm <- new("SimpleRleList") 
  for (chr in names(coverage)) {
    x.c <-  as.vector(coverage[[chr]])
    cov.norm[[chr]] <- Rle(round(atransform(x.c, total.cov),10))
  }
  return(cov.norm)
}


######################################################
## zscore transforamtion ofhte coverage
######################################################
zscore.cov <- function(covs) {
  v <- as.vector(unlist(covs))
  flt <- v!=0
  m <- mean(v[flt])
  std <- sd(v[flt])
  covss <- round(((covs-m) / std), 2)
  return(covss)
}


######################################################
## merge GRanges to a coverage list
######################################################
merge_grfiles2cov <- function(...,name="out.cov",gwidth=200,fix="start",
                              chrs=c("chrXVII","chrXVIII","chrXIX","chrI","chrII","chrIII","chrIV","chrIX","chrV","chrVI","chrVII","chrVIII","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")){
  require(GenomicRanges)
  ff <- c(...)
  if(length(grep(".gr",ff))==0){
    ff <- paste(ff,".gr",sep="")
  }

  g <- GRanges()
  for(i in 1:length(ff)){
    load(ff[i])
    gr <- subset(gr,seqnames(gr)%in%chrs)
    gr <- keepSeqlevels(gr,c(chr.flt,spiked.chr))
    g <- c(g,gr)
  }
  if(!is.null(gwidth)){
    g <- resize(g,width = gwidth,fix=fix)
  }
  cov <- coverage(g)
  cov <- standardize.cov(cov)
  save(cov,file=paste(name))
}


######################################################
## TSS centered matrix at bp resolution
######################################################
getMatrix_TSS <- function(covs, points, win=c(-200,4000),cumul.df=T,bin.bw=0) {
  m <- matrix(nrow=length(points$chr), ncol=win[2]-win[1]+1)
  
  for (i in 1:length(points$chr)) {
    if (points$TSS[i] > win[2] & (points$TSS[i] + win[2]) < length(covs[[as.character(points$chr[i])]]) & points$chr[i] %in% names(covs)) {
      if (points$strand[i]=="+") {
        m[i,] <- as.vector( window( covs[[as.character(points$chr[i])]], points$TSS[i]+win[1], points$TSS[i]+win[2] ) )
      } else {
        m[i,] <- rev(as.vector( window( covs[[as.character(points$chr[i])]], points$TSS[i] - win[2], points$TSS[i] - win[1] )) )
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

getMatrix_TSS.stranded <- function(covsp,covsn, points, win=c(-200,4000),cumul.df=T,bin.bw=0) {
  m <- matrix(nrow=length(points$chr), ncol=win[2]-win[1]+1)
  
  for (i in 1:length(points$chr)) {
    if (points$TSS[i] > abs(win[1]) & (points$TSS[i]+win[2]) < length(covsp[[as.character(points$chr[i])]]) & (points$TSS[i]+win[2]) < length(covsn[[as.character(points$chr[i])]]) & points$chr[i] %in% names(covsp)) {
      if (points$strand[i]=="+") {
        m[i,] <- as.vector( window( covsp[[as.character(points$chr[i])]], points$TSS[i]+win[1], points$TSS[i]+win[2] ) )
      }
      else{
        m[i,] <- rev(as.vector( window( covsn[[as.character(points$chr[i])]], points$TSS[i] - win[2], points$TSS[i] - win[1] )) )
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

######################################################
## TTS centered matrix at bp resolution
######################################################
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

getMatrix_TTS.stranded <- function(covsp,covsn, points, win=c(-200,4000),cumul.df=T,bin.bw=0) {
  m <- matrix(nrow=length(points$chr), ncol=win[2]-win[1]+1)
  
  for (i in 1:length(points$chr)) {
    if (points$TTS[i] > abs(win[1]) & (points$TTS[i]+win[2]) < length(covsp[[as.character(points$chr[i])]]) & (points$TTS[i]+win[2]) < length(covsn[[as.character(points$chr[i])]]) & points$chr[i] %in% names(covsp)) {
      if (points$strand[i]=="+") {
        m[i,] <- as.vector( window( covsp[[as.character(points$chr[i])]], points$TTS[i]+win[1], points$TTS[i]+win[2] ) )
      }
      else{
        m[i,] <- rev(as.vector( window( covsn[[as.character(points$chr[i])]], points$TTS[i] - win[2], points$TTS[i] - win[1] )) )
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

######################################################
## feature scaled bp resolution mattrix
######################################################
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


######################################################
## unscaled bp resolution matrices with length.cut
######################################################
get.unscaled.genebody.matrix <- function (covs, points, Margin =500,gene.len.cut=4000,bw=10,align="TSS"){
  #Multigene heatmaps were produced by binning profiles every 10 bp, sorting genes by length, then aligning all at the 50 TSS.
  #Flanking regions 500 bp upstream and downstream are also shown. Genes longer than 4 kb are allowed to run off the right side of the plot
  
  ## Scaling the signal to 10Mio
  points$width <- points$end-points$start
  points <- points[order(points$width,decreasing = F),]
  
  ncol=gene.len.cut+2*Margin
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
      
      if(align=="TSS" & length(v)>gene.len.cut){
        v <- v[1:gene.len.cut]      
      }
      else if(align=="TTS" & length(v)>gene.len.cut){ 
        v <- tail(v,4000)      
      }
      
      l <- c(up,v,dn)
      m[i,1:length(l)] <- l
    }
  }
  for(i in seq(1,ncol-1,by=bw)){
    r <- cbind(r,rowSums(m[,i:(i+bw-1)],na.rm = T)  )
  }
  rownames(r) <- points$tracking_id
  
  return(r)
}

get.unscaled.directional.genebody.matrix <- function (covsp,covsn, points, Margin =500,gene.len.cut=4000,bw=10,align="TSS"){
  #Multigene heatmaps were produced by binning profiles every 10 bp, sorting genes by length, then aligning all at the 50 TSS.
  #Flanking regions 500 bp upstream and downstream are also shown. Genes longer than 4 kb are allowed to run off the right side of the plot
  
  ## Scaling the signal to 10Mio
  points$width <- points$end-points$start
  points <- points[order(points$width,decreasing = F),]
  
  ncol=gene.len.cut+2*Margin
  m <- matrix(NA,nrow=dim(points)[1],ncol=ncol)
  r <- c()
  
  ##### Collecting gene upsream/downstream/gene-body density
  for (i in 1:length(points$chr)) {
    if (points$start[i] > Margin  & (points$end[i] + Margin) < length(covsp[[as.character(points$chr[i])]]) & (points$end[i] + Margin) < length(covsn[[as.character(points$chr[i])]]) & points$chr[i] %in% names(covsp)) {
      
      if (points$strand[i]=="+") {
        up<- as.vector( window( covsp[[as.character(points$chr[i])]], points$start[i]-Margin, points$start[i]-1 ))
        v <- as.vector( window( covsp[[as.character(points$chr[i])]], points$start[i], points$end[i]))
        dn<- as.vector( window( covsp[[as.character(points$chr[i])]], points$end[i]+1, points$end[i]+Margin ))
        
      }else {
        up<- rev(as.vector( window( covsn[[as.character(points$chr[i])]], points$end[i]+1, points$end[i]+Margin )))
        v <- rev(as.vector( window( covsn[[as.character(points$chr[i])]], points$start[i], points$end[i])))
        dn<- rev(as.vector( window( covsn[[as.character(points$chr[i])]], points$start[i]-Margin, points$start[i]-1 )))
      }
      
      if(align=="TSS" & length(v)>gene.len.cut){
        v <- v[1:gene.len.cut]      
      }
      else if(align=="TTS" & length(v)>gene.len.cut){ 
        v <- tail(v,4000)      
      }
      
      l <- c(up,v,dn)
      m[i,1:length(l)] <- l
    }
  }
  for(i in seq(1,ncol-1,by=bw)){
    r <- cbind(r,rowSums(m[,i:(i+bw-1)],na.rm = T)  )
  }
  rownames(r) <- points$tracking_id
  return(r)
}

get.unscaled.whole.genebody.matrix <- function (covs, points, Margin =500,bw=10){
  #Multigene heatmaps were produced by binning profiles every 10 bp, sorting genes by length, then aligning all at the 50 TSS.
  #Flanking regions 500 bp upstream and downstream are also shown. Genes longer than 4 kb are allowed to run off the right side of the plot
  
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


######################################################
## Matrix to metagene plot
######################################################
metagene.signal.plot <- function(m,smooth=T){
  df_m <- data.frame(pos = 1:ncol(m),
                     signal = apply(m,2,mean,na.rm=TRUE),
                     se=apply(m,2, function(x) sd(x,na.rm=TRUE)/sqrt(length(x))),
                     sd=apply(m,2,sd,na.rm=TRUE)
                     )
  if(smooth==T){
    df_m$signal <- ksmooth(1:nrow(df_m),df_m$signal,bandwidth = 50)$y
  }
  df_m$ci <- df_m$se*1.96
  df_m
}

######################################################
## Venn diagrams
######################################################
Venn2 <- function(set1, set2, names,plot=TRUE,return.df=T) {
  require(limma)
  stopifnot( length(names) == 2)
  set1 <- unique(set1)
  set2 <- unique(set2)
  
  # Form universe as union of all three sets
  universe <- sort( unique( c(set1, set2) ) )
  
  Counts <- matrix(0, nrow=length(universe), ncol=2)
  colnames(Counts) <- names
  
  for (i in 1:length(universe))
  {
    Counts[i,1] <- universe[i] %in% set1
    Counts[i,2] <- universe[i] %in% set2
  }
  
  #if(plot==TRUE){ 
  vennDiagram( vennCounts(Counts),circle.col = c(rgb(1,0,0,0.2),rgb(0,0,1,0.2)),cex=c(1.5,1.4,1.4)  ) 
  #}
  
  rownames(Counts) <- universe
  
  df1 <- list(AB=rownames(subset(Counts,Counts[,1]==1 & Counts[,2]==1)),
              A=rownames(subset(Counts,Counts[,1]==1 & Counts[,2]==0)),
              B=rownames(subset(Counts,Counts[,1]==0 & Counts[,2]==1))
  )
  
  names(df1) <- c(paste(names[1],names[2]), names[1],names[2])
  df2 <- c()
  for (name in names(df1)) {
    df2 <- rbind(df2, data.frame(condition=rep(name,length(df1[[name]])),tracking_id=df1[[name]]))
  }
  df2 <- unique(df2)
  if(return.df==T){
    return(df2)
  }
}

Venn3 <- function(set1, set2, set3, names) {
  require(limma)
  stopifnot( length(names) == 3)
  set1 <- unique(set1)
  set2 <- unique(set2)
  set3 <- unique(set3)
  
  # Form universe as union of all three sets
  universe <- sort( unique( c(set1, set2, set3) ) )
  
  Counts <- matrix(0, nrow=length(universe), ncol=3)
  colnames(Counts) <- names
  
  for (i in 1:length(universe))
  {
    Counts[i,1] <- universe[i] %in% set1
    Counts[i,2] <- universe[i] %in% set2
    Counts[i,3] <- universe[i] %in% set3
  }
  
  vennDiagram( vennCounts(Counts),circle.col = c(rgb(1,0,0,0.2),rgb(0,1,0,0.2),rgb(0,0,1,0.2)) )
  rownames(Counts) <- universe
  
  df1 <- list(ABC=rownames(subset(Counts,Counts[,1]==1 & Counts[,2]==1 & Counts[,3]==1)),
              A=rownames(subset(Counts,Counts[,1]==1 & Counts[,2]==0 & Counts[,3]==0)),
              B=rownames(subset(Counts,Counts[,1]==0 & Counts[,2]==1 & Counts[,3]==0)),
              C=rownames(subset(Counts,Counts[,1]==0 & Counts[,2]==0 & Counts[,3]==1)),
              AB=rownames(subset(Counts,Counts[,1]==1 & Counts[,2]==1 & Counts[,3]==0)),
              AC=rownames(subset(Counts,Counts[,1]==1 & Counts[,2]==0 & Counts[,3]==1)),
              BC=rownames(subset(Counts,Counts[,1]==0 & Counts[,2]==1 & Counts[,3]==1))
  )
  
  names(df1) <- c(paste(names[1],names[2],names[3]), names[1],names[2],names[3],
                  paste(names[1],names[2]),paste(names[1],names[3]), paste(names[2],names[3]))
  
  return(df1)
}

Venn4 <- function(set1, set2, set3, set4,names) {
  require(limma)
  stopifnot( length(names) == 4)
  set1 <- unique(set1)
  set2 <- unique(set2)
  set3 <- unique(set3)
  set4 <- unique(set4)
  set1 <- set1[complete.cases(set1)]
  set2 <- set2[complete.cases(set2)]
  set3 <- set3[complete.cases(set3)]
  set4 <- set4[complete.cases(set4)]
  
  # Form universe as union of all three sets
  universe <- sort( unique( c(set1, set2, set3,set4) ) )
  
  Counts <- matrix(0, nrow=length(universe), ncol=4)
  colnames(Counts) <- names
  
  for (i in 1:length(universe))
  {
    Counts[i,1] <- universe[i] %in% set1
    Counts[i,2] <- universe[i] %in% set2
    Counts[i,3] <- universe[i] %in% set3
    Counts[i,4] <- universe[i] %in% set4
    
  }
  
  col = c(rgb(1,0,0,0.2),rgb(0,1,0,0.2),rgb(0,0,1,0.2),rgb(0,0,0,0.2))
  vennDiagram( vennCounts(Counts),circle.col = col,cex=c(1.5,1.4,1.4) )
  rownames(Counts) <- universe
  
  df1 <- list(ABCD=rownames(subset(Counts,Counts[,1]==1 & Counts[,2]==1 & Counts[,3]==1 & Counts[,4]==1)),
              A=rownames(subset(Counts,Counts[,1]==1 & Counts[,2]==0 & Counts[,3]==0 & Counts[,4]==0)),
              B=rownames(subset(Counts,Counts[,1]==0 & Counts[,2]==1 & Counts[,3]==0 & Counts[,4]==0)),
              C=rownames(subset(Counts,Counts[,1]==0 & Counts[,2]==0 & Counts[,3]==1 & Counts[,4]==0)),
              D=rownames(subset(Counts,Counts[,1]==0 & Counts[,2]==0 & Counts[,3]==0 & Counts[,4]==1)),
              
              AB=rownames(subset(Counts,Counts[,1]==1 & Counts[,2]==1 & Counts[,3]==0 & Counts[,4]==0)),
              AC=rownames(subset(Counts,Counts[,1]==1 & Counts[,2]==0 & Counts[,3]==1 & Counts[,4]==0)),
              BC=rownames(subset(Counts,Counts[,1]==0 & Counts[,2]==1 & Counts[,3]==1 & Counts[,4]==0)),
              AD=rownames(subset(Counts,Counts[,1]==1 & Counts[,2]==0 & Counts[,3]==0 & Counts[,4]==1)),
              BD=rownames(subset(Counts,Counts[,1]==0 & Counts[,2]==1 & Counts[,3]==0 & Counts[,4]==1)),
              CD=rownames(subset(Counts,Counts[,1]==0 & Counts[,2]==0 & Counts[,3]==1 & Counts[,4]==1)),
              
              ABC=rownames(subset(Counts,Counts[,1]==1 & Counts[,2]==1 & Counts[,3]==1 & Counts[,4]==0)),
              ABD=rownames(subset(Counts,Counts[,1]==1 & Counts[,2]==1 & Counts[,3]==0 & Counts[,4]==1)),
              ACD=rownames(subset(Counts,Counts[,1]==1 & Counts[,2]==0 & Counts[,3]==1 & Counts[,4]==1)),
              BCD=rownames(subset(Counts,Counts[,1]==0 & Counts[,2]==1 & Counts[,3]==1 & Counts[,4]==1))
  )
  
  names(df1) <- c(paste(names[1],names[2],names[3],names[4]),
                  names[1],
                  names[2],
                  names[3],
                  names[4],
                  paste(names[1],names[2]),
                  paste(names[1],names[3]),
                  paste(names[2],names[3]),
                  paste(names[1],names[4]),
                  paste(names[2],names[4]),
                  paste(names[3],names[4]),
                  paste(names[1],names[2],names[3]),
                  paste(names[1],names[2],names[4]),
                  paste(names[1],names[3],names[4]),
                  paste(names[2],names[3],names[4])
  )
  
  return(df1)
  
}


######################################################
## Gene Ontology enrichments
######################################################
GOStats.enrichment <- function(genes,universe,p.val,species="Scer"){
  require("GOstats")
  require("GSEABase")
  
  if(species=="Scer"){
    #go <- read.delim("A:/work/yeast_Annotations/Scer_GO.txt",header=F)
    #frame <- read.delim("A:/work/yeast_Annotations/YeastGenome_GoTerms.txt",header=T)
    #frame <- subset(frame,frame$GO.Term.Accession%in%go$V6)
    #frame <- frame[complete.cases(frame),]
    #frame <- read.delim("A:/work/yeast_Annotations/YeastGenome_GoTerms.txt",header=T)
    #universe=levels((frame$Ensembl.Gene.ID))
    
    #goframeData = frame[,c(1,2,3)]
    #names(goframeData) = c("go_id", "Evidence", "gene_id")
    #goFrame=GOFrame(goframeData,organism="Saccharomyces cerevisiae")
    #goAllFrame=GOAllFrame(goFrame)
    #gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
    #save(gsc,file="A:/work/yeast_Annotations/SacCer3 gene ontologies_slim.RData")
    
    load("A:/work/yeast_Annotations/SacCer3 gene ontologies_slim.RData")
  }else if(species=="mm9" | species=="mm10"){
    #d <- read.delim("A:/work/mm9_Annotations/mm9_MGI2ucsc.txt")
    #d <- unique(d[,c(2,3)])
    #b <- read.delim("A:/work/mm9_Annotations/mgigene_association.mgi",comment.char = "!",header=F)
    #b <- unique(b[,c(2,5,7)])
    #names(d)
    #names(b)[1] <- "MGI.ID"
    #d <- merge(d,b,by="MGI.ID",all.x=T)
    #d <- unique(d[,c(2:4)])
    #d$UCSC.ID <- ifelse(d$UCSC.ID=="",NA,as.character(d$UCSC.ID))
    #d <- d[complete.cases(d),]
    #write.table(d,file="A:/work/mm9_Annotations/mm9_GOTerms.txt",quote = F,sep = "\t",row.names = F)
    #frame <- read.delim("A:/work/mm9_Annotations/mm9_GOTerms.txt",header=T)
    #frame <- frame[complete.cases(frame),]
    #goframeData = frame[,c(2,3,1)]
    #names(goframeData) = c("go_id", "Evidence", "gene_id")
    #frame <- frame[!frame$V7%in%c("ND","IC","TAS","NAS"),]
    #goFrame=GOFrame(goframeData,organism="Mus musculus")
    #goAllFrame=GOAllFrame(goFrame)
    #gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
    #save(gsc,file="A:/work/mm9_Annotations/mm9 gene ontologies")
    
    load("A:/work/mm9_Annotations/mm9 gene ontologies")
    
  }
  mf <- bp <- cc <- c()
  params <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params",geneSetCollection=gsc,
                               geneIds = genes,universeGeneIds = universe,
                               ontology = "MF",pvalueCutoff = p.val,
                               conditional = FALSE,testDirection = "over")
  mf <- summary(hyperGTest(params))
  if(length(mf)>0){
    mf <- cbind(mf,Category="Molecular Function")
    names(mf)[1] <- c("ID")
  }
  
  params <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params",geneSetCollection=gsc,
                               geneIds = genes,universeGeneIds = universe,
                               ontology = "BP",pvalueCutoff = p.val,
                               conditional = FALSE,testDirection = "over")
  bp <- summary(hyperGTest(params))
  if(length(bp)>0){
    bp <- cbind(bp,Category="Biological Process")
    names(bp)[1] <- c("ID")
  }
  
  params <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params",geneSetCollection=gsc,
                               geneIds = genes,universeGeneIds = universe,
                               ontology = "CC",pvalueCutoff = p.val,
                               conditional = FALSE,testDirection = "over")
  cc <- summary(hyperGTest(params))
  if(length(cc$Pvalue)>0){
    cc <- cbind(cc,Category="Cellular Component")
    names(cc)[1] <- c("ID")
  }
  
  if(length(mf$Pvalue)==0){mf = data.frame(ID=NA,Pvalue=NA,OddsRatio=NA,ExpCount=NA,Count=NA,Size=NA,Term=NA,Category=NA)}
  if(length(bp$Pvalue)==0){bp = data.frame(ID=NA,Pvalue=NA,OddsRatio=NA,ExpCount=NA,Count=NA,Size=NA,Term=NA,Category=NA)}
  if(length(cc$Pvalue)==0){cc = data.frame(ID=NA,Pvalue=NA,OddsRatio=NA,ExpCount=NA,Count=NA,Size=NA,Term=NA,Category=NA)}
  g <- rbind(mf,bp,cc)
  g <- g[complete.cases(g),]
  g <- g[order(g$Pvalue),]
  
  return(g)
}


######################################################
## Gini index for grRanges
######################################################
calculate.Ginni.Index.forGR <- function(gr,genomicBins,filt.gr,gwidth=1,fix="center"){
  if(gwidth>0){
    gr <- resize(gr,width = gwidth,fix=fix)
  }
  if(length(filt.gr)>0){
    gr <- gr[!overlapsAny(gr,filt.gr,ignore.strand=T)]
  }
  d <- countOverlaps(genomicBins,gr,ignore.strand=T)
  d <- sort(d,decreasing = F)
  n <- length(d)
  p <- n:1*d
  p <- (n+1 - 2*(sum(as.numeric(n:1*d))/sum(d)))/n
  return(p)
}


######################################################
## Pairwise scatter plots
######################################################
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
