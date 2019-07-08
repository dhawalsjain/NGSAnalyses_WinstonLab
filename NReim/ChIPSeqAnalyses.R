######################################################################################
### Functions
######################################################################################
library(GenomicRanges)
library(parallel)
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
    #gr <- subset(gr,seqnames(gr)%in%c(chr.flt,spiked.chr))
    #gr <- keepSeqlevels(gr,c(chr.flt,spiked.chr))
    gr <- subset(gr,seqnames(gr)%in%c(chr.flt))
    gr <- keepSeqlevels(gr,chr.flt)
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
    #gr <- keepSeqlevels(gr,c(chr.flt,spiked.chr))
    gr <- keepSeqlevels(gr,chr.flt)
    g <- c(g,gr)
  }
  g <- resize(g,200,fix="start")
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
    #gr <- keepSeqlevels(gr,c(chr.flt,spiked.chr))
    gr <- keepSeqlevels(gr,chr.flt)
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

image.scale <- function(z,zlim=NULL,mid.point=0,cols=c("blue","white","red"),cut.at = c(1,1),enforce.lims=T,length.out=256) {
  ## expected z formats:
  ## z>=0, mid.point=0, one-sided, black&while
  ## z>=0, mid.point=1, blue&white%red
  ## -zlim<z<zlim, mid.point=0,blue&white%red
  
  r=range(z,na.rm=T)
  if(is.null(zlim)) {
    zlim=r
    enforce.lims=F
  }
  if(!enforce.lims) {
    if(zlim[1] < r[1]) zlim[1]=r[1]
    if(zlim[2] > r[2]) zlim[2]=r[2]
  }
  
  z[which(z<zlim[1])] = zlim[1]
  z[which(z>zlim[2])] = zlim[2]
  
  scale=list()
  
  if(length(cut.at)==0) cut.at=c(1,1)
  if(length(cut.at)==1) cut.at=rep(cut.at,2)
  
  plotmin = mid.point+cut.at[1]*(zlim[1]-mid.point)
  if(!enforce.lims & r[1]>plotmin) plotmin=r[1]
  
  plotmax = mid.point+cut.at[2]*(zlim[2]-mid.point)
  if(!enforce.lims & r[2]<plotmax) plotmax=r[2]
  
  scale$x=seq(plotmin,plotmax,length.out=length.out)
  scale$y=0:1
  scale$z=matrix(mid.points(scale$x),ncol=1)
  scale$z[which(scale$z<zlim[1])] = zlim[1]
  scale$z[which(scale$z>zlim[2])] = zlim[2]
  
  bias=1
  if(mid.point>zlim[1]) {
    x.mid = (mid.point-zlim[1])/(zlim[2]-zlim[1])
    if(x.mid>0) {
      bias=(-log2(x.mid))^1.03
    }
  }
  col=colorRampPalette(cols,space="Lab",bias=bias)(length.out)
  col.z = col
  if(enforce.lims) {
    z100 = seq(zlim[1],zlim[2],length.out=length.out)
    ind = which.min(abs(r[1]-z100)):which.min(abs(r[2]-z100))
    col.z=col.z[ind]
  }
  
  return(list(z=z,scale=scale,col.scale=col,col.z=col.z))
}

mid.points <-function(x) {
  l=length(x)
  return((x[2:l-1] + x[2:l] )/2)
}

get.unscaledgene.heatplot.list <- function(r,mid.point=0,zlim=c(0,1),cols=c("blue4","white","red2"),genes,x= seq(-500,4500,10),gene.len.cut=4000,align="TSS"){
  if(!"width"%in%names(genes)){
    genes$width <- genes$end- genes$start
  }
  mid.point = mid.point
  cols=cols#c("white","black")
  z=image.scale(z = apply(r,2,rev),zlim = zlim,mid.point = mid.point,cols=cols)
  border.lines=list()
  l <- genes$width
  l[which(l>gene.len.cut)]=gene.len.cut
  
  if(align=="CPS"){
    border.lines[[1]] = rep(gene.len.cut,length(l))
    border.lines[[2]] = gene.len.cut-l[order(l,decreasing = F)]
  }else{
    border.lines[[1]] = rep(0,length(l))
    border.lines[[2]] = l[order(l,decreasing = F)]
  }
  
  h <- list(x=x,z=z,zlim=zlim,border.lines=border.lines)
  return(h)
}

#sbatch -t 0-11:59 -n 6 --mem-per-cpu=2G -p short --wrap "Rscript Natalia_ChIP_Final_filt.R"

######################################################################################
### Natalia, ChIPSeq paper figures, calculate data frames
######################################################################################

## Multigene whole gene body 
if(T){
  files <- Sys.glob("*.gr")
  keys <- read.delim("SampleKey.txt")
  z <- unique(as.character(keys$Input))
  z <- z[z!=""]
  names(z) <- c("Spn1_round2_Inp", "Spn1_round2_Inp", "Spn1_round2_Inp", "Spn1_round2_Inp", 
                "Spn1_round1_Inp", "Spn1_round1_Inp", "Spn1_round1_Inp", "Spn1_round1_Inp", 
                "Spn1_HA.Set2_Inp", "Spn1_HA.Set2_Inp", "Spn1_HA.Set2_Inp", 
                "Spn1_HA.Set2_Inp", "Spn1_HA.Spt6_Inp", "Spn1_HA.Spt6_Inp", 
                "Spn1_HA.Spt6_Inp", "Spn1_HA.Spt6_Inp")
  
  load("/n/data1/hms/dbmi/park/Dhawal/Genomes/yeast_Annotations/ScSp_SGD.RevisedGenes.RData")
  #load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenes.RData")
  
  genes <- subset(genes,genes$verification=="Verified" & genes$species=="Scer" & genes$chr!="chrM")
  genes <- genes[order(genes$width,decreasing = F),]
  keys <- split(keys,keys$Factor)
  
  if(F){
    mclapply(1:length(keys),function(i){
      l <- list()
      cat("profile: ",names(keys)[i],"\n")
      depl <- as.character(keys[[i]]$ChIP)
      nondepl <- depl[-grep("Depl",depl)]
      depl <- depl[grep("Depl",depl)]
      l[["Depleted-1"]] <- merge_files(depl[1],genes = genes)
      l[["Depleted-2"]] <- merge_files(depl[2],genes = genes)
      l[["Non-depleted-1"]] <- merge_files(nondepl[1],genes = genes)
      l[["Non-depleted-2"]] <- merge_files(nondepl[2],genes = genes)
      save(l,file = paste("Revised_MultigeneHeatmap_",names(keys)[i],".RData",sep=""))
      return(1)
    },mc.cores = 12)
    
  }
  if(F){
    mclapply(1:length(keys),function(i){
      l <- list()
      cat("profile: ",names(keys)[i],"\n")
      depl <- as.character(keys[[i]]$ChIP)
      nondepl <- depl[-grep("Depl",depl)]
      depl <- depl[grep("Depl",depl)]
      l[["Depleted-1"]] <- merge_files_tss(depl[1],genes = genes)
      l[["Depleted-2"]] <- merge_files_tss(depl[2],genes = genes)
      l[["Non-depleted-1"]] <- merge_files_tss(nondepl[1],genes = genes)
      l[["Non-depleted-2"]] <- merge_files_tss(nondepl[2],genes = genes)
      save(l,file = paste("Revised_MultigeneHeatmap_TSS_",names(keys)[i],".RData",sep=""))
      return(1)
    },mc.cores = 6)
    
  }
  if(T){
    mclapply(1:length(keys),function(i){
      l <- list()
      cat("profile: ",names(keys)[i],"\n")
      depl <- as.character(keys[[i]]$ChIP)
      nondepl <- depl[-grep("Depl",depl)]
      depl <- depl[grep("Depl",depl)]
      l[["Depleted-1"]] <- merge_files_tts(depl[1],genes = genes)
      l[["Depleted-2"]] <- merge_files_tts(depl[2],genes = genes)
      l[["Non-depleted-1"]] <- merge_files_tts(nondepl[1],genes = genes)
      l[["Non-depleted-2"]] <- merge_files_tts(nondepl[2],genes = genes)
      save(l,file = paste("Revised_MultigeneHeatmap_CPS_",names(keys)[i],".RData",sep=""))
      return(1)
    },mc.cores = 6)
    
  }
  if(F){
    mclapply(1:length(keys),function(i){
      l <- list()
      cat("profile: ",names(keys)[i],"\n")
      depl <- as.character(keys[[i]]$ChIP)
      nondepl <- depl[-grep("Depl",depl)]
      depl <- depl[grep("Depl",depl)]
      l[["Depleted"]] <- merge_files_tss(depl[1],depl[2],genes = genes)
      l[["Non-depleted"]] <- merge_files_tss(nondepl[1],nondepl[2],genes = genes)
      save(l,file = paste("Revised_MultigeneHeatmap_TSScomb_",names(keys)[i],".RData",sep=""))
      return(1)
    },mc.cores = 4)
    
  }
  if(F){
    files <- c("HA.Spt6_ChIP_R1","S5P_ChIP_R2","S2P_ChIP_R2",
               "H3K4me3_ChIP_R2","H3K36me2_ChIP_R2","H3K36me3_ChIP_R1")
    mclapply(1:6,function(i){
      l <- list()
      cat("profile: ",files[i],"\n")
      depl <- paste("/n/data1/hms/dbmi/park/Dhawal/Project_nreim/MayNov17/cov/Spn1.Depl_",files[i],"_combRelIN.cov",sep="")
      nondepl <- paste("/n/data1/hms/dbmi/park/Dhawal/Project_nreim/MayNov17/cov/Spn1_",files[i],"_combRelIN.cov",sep="")
      l[["Depleted"]] <- merge_files_cov(depl,genes = genes)
      l[["Non-depleted"]] <- merge_files_cov(nondepl,genes = genes)
      save(l,file = paste("Rel_Revised_MultigeneHeatmap_",names(keys)[i],".RData",sep=""))
      return(1)
    },mc.cores = 1)
  }
  if(F){
    mclapply(1:length(keys),function(i){
      l <- list()
      cat("profile: ",names(keys)[i],"\n")
      depl <- as.character(keys[[i]]$ChIP)
      nondepl <- depl[-grep("Depl",depl)]
      depl <- depl[grep("Depl",depl)]
      l[["Depleted"]] <- merge_files_scaleHeat(depl[1],depl[2],genes = genes)
      l[["Non-depleted"]] <- merge_files_scaleHeat(nondepl[1],nondepl[2],genes = genes)
      save(l,file = paste("Revised_MultigeneHeatmap_ScHeat_",names(keys)[i],".RData",sep=""))
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

## Metagene plots for splike in normalizations
if(F){
  cat("Counting metagene signal\n")
  covdir="/n/data1/hms/dbmi/park/Dhawal/Project_nreim/MayNov17/cov"
  mc.cores=12
  
  chr.flt <- c("chrI","chrII","chrIII","chrIV","chrIX","chrV","chrVI","chrVII","chrVIII","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
  spiked.chr <- c("chrXVII","chrXVIII","chrXIX")
  load("/n/data1/hms/dbmi/park/Dhawal/Genomes/yeast_Annotations/ScSp_SGD.RevisedGenes.RData")
  load("/n/data1/hms/dbmi/park/Dhawal/Genomes/yeast_Annotations/Natalia_genes10FPKM.RData")
  genes <- subset(genes,genes$verification=="Verified" & genes$chr!="chrM")
  genes <- genes[genes$width>=1000,]
  g1 <- genes[genes$species=="Spom",]
  genes <- genes[genes$tracking_id%in%as.character(expressed.genes),]
  genes <- rbind(genes,g1)
  rm(g1)
  
  mytss <- function(covs,genes,profile){
    a <- getMatrix_TSS(covs,genes[genes$width>500 & genes$species=="Scer",],win = c(-500,4000),cumul.df = T,bin.bw = 10)  ## TSS
    b <- getMatrix_TSS(covs,genes[genes$width>500 & genes$species=="Spom",],win = c(-500,4000),cumul.df = T,bin.bw = 10)  ## TSS
    a$species <- "S.cer"
    b$species <- "S.pom"
    a <- rbind(a,b)
    a$profile  <- profile
    return(a)
  }
  
  mytts <- function(covs,genes,profile){
    a <- getMatrix_TTS(covs,genes[genes$width>500 & genes$species=="Scer",],win = c(-4000,500),cumul.df = T,bin.bw = 10)  ## TTS
    b <- getMatrix_TTS(covs,genes[genes$width>500 & genes$species=="Spom",],win = c(-4000,500),cumul.df = T,bin.bw = 10)  ## TTS
    a$species <- "S.cer"
    b$species <- "S.pom"
    a <- rbind(a,b)
    a$profile  <- profile
    return(a)
  }
  
  myscaledgene <- function(covs,genes,profile){
    a <- get.feature.scaled.matrix(covs = covs,points = genes[genes$width>500 & genes$species=="Scer",],cumul.df = T,bin.bw = 10)  ## Scaled genes, that are >500bp long
    b <- get.feature.scaled.matrix(covs = covs,points = genes[genes$width>500 & genes$species=="Spom",],cumul.df = T,bin.bw = 10)  ## Scaled genes, that are >500bp long
    a$species <- "S.cer"
    b$species <- "S.pom"
    a <- rbind(a,b)
    a$profile  <- profile
    return(a)
  }
  
  cat("Along Scaled ChIP and Input samples\n")
  setwd(covdir)
  files <- Sys.glob("*_st.cov")
  pl1 <- c()  ## TSS
  pl2 <- c()  ## TTS
  pl3 <- c()  ## Scaled Gene
  l <- mclapply(files,function(f){
    cat("TSS: ",f,"\n")
    load(f)
    return(mytss(cov,genes,profile = sub("_st.cov","",f)))
  },mc.cores = mc.cores)
  for(i in 1:length(l)){ pl1 <- rbind(pl1,l[[i]]) }
  rm(l)
  l <- mclapply(files,function(f){
    cat("TTS: ",f,"\n")
    load(f)
    return(mytts(cov,genes,profile = sub("_st.cov","",f)))
  },mc.cores = mc.cores)
  for(i in 1:length(l)){ pl2 <- rbind(pl2,l[[i]]) }
  rm(l)
  l <- mclapply(files,function(f){
    cat("ScaledGene: ",f,"\n")
    load(f)
    return(myscaledgene(cov,genes,profile = sub("_st.cov","",f)))
  },mc.cores = mc.cores)
  for(i in 1:length(l)){ pl3 <- rbind(pl3,l[[i]]) }
  rm(l)
  cat("\n")
  
  cat("Along Input divided samples\n")
  setwd(covdir)
  files <- Sys.glob("*_IN.cov")
  pl13 <- c()  ## TSS
  pl14 <- c()  ## TTS
  pl15 <- c()  ## Scaled Gene
  l <- mclapply(files,function(f){
    cat("TSS: ",f,"\n")
    load(f)
    return(mytss(cov,genes,profile = sub("_IN.cov","",f)))
  },mc.cores = mc.cores)
  for(i in 1:length(l)){ pl13 <- rbind(pl13,l[[i]]) }
  rm(l)
  l <- mclapply(files,function(f){
    cat("TTS: ",f,"\n")
    load(f)
    return(mytts(cov,genes,profile = sub("_IN.cov","",f)))
  },mc.cores = mc.cores)
  for(i in 1:length(l)){ pl14 <- rbind(pl14,l[[i]]) }
  rm(l)
  l <- mclapply(files,function(f){
    cat("ScaledGene: ",f,"\n")
    load(f)
    return(myscaledgene(cov,genes,profile = sub("_IN.cov","",f)))
  },mc.cores = mc.cores)
  for(i in 1:length(l)){ pl15 <- rbind(pl15,l[[i]]) }
  rm(l)
  cat("\n")
  
  ## Combined replicates
  cat("Along Input divided samples\n")
  setwd(covdir)
  files <- Sys.glob("*_comb.cov")
  pl4 <- c()  ## TSS
  pl5 <- c()  ## TTS
  pl6 <- c()  ## Scaled Gene
  l <- mclapply(files,function(f){
    cat("TSS: ",f,"\n")
    load(f)
    return(mytss(cov,genes,profile = sub("_comb.cov","",f)))
  },mc.cores = mc.cores)
  for(i in 1:length(l)){ pl4 <- rbind(pl4,l[[i]]) }
  rm(l)
  l <- mclapply(files,function(f){
    cat("TTS: ",f,"\n")
    load(f)
    return(mytts(cov,genes,profile = sub("_comb.cov","",f)))
  },mc.cores = mc.cores)
  for(i in 1:length(l)){ pl5 <- rbind(pl5,l[[i]]) }
  rm(l)
  l <- mclapply(files,function(f){
    cat("ScaledGene: ",f,"\n")
    load(f)
    return(myscaledgene(cov,genes,profile = sub("_comb.cov","",f)))
  },mc.cores = mc.cores)
  for(i in 1:length(l)){ pl6 <- rbind(pl6,l[[i]]) }
  rm(l)
  cat("\n")
  
  cat("Along SPP genrated ba kgrorund subtracted signal\n")
  setwd("/n/data1/hms/dbmi/park/Dhawal/Project_nreim/MayNov17/spp")
  files <- Sys.glob("*_density.wig")
  pl7 <- c()  ## TSS
  pl8 <- c()  ## TTS
  pl9 <- c()  ## Scaled Gene
  l <- mclapply(files,function(f){
    cat("TSS: ",f,"\n")
    gr <- read.delim(f,header=F,comment.char = "t",sep = " ")
    gr <- with(gr,GRanges(V1,IRanges(V2,V3),"+",score=V4))
    gr <- keepSeqlevels(gr,c(chr.flt,spiked.chr))
    cov <- standardize.cov(coverage(gr,weight = gr$score))
    return(mytss(cov,genes,profile = sub("_density.wig","",f)))
  },mc.cores = mc.cores)
  for(i in 1:length(l)){ pl7 <- rbind(pl7,l[[i]]) }
  rm(l)
  l <- mclapply(files,function(f){
    cat("TTS: ",f,"\n")
    gr <- read.delim(f,header=F,comment.char = "t",sep = " ")
    gr <- with(gr,GRanges(V1,IRanges(V2,V3),"+",score=V4))
    gr <- keepSeqlevels(gr,c(chr.flt,spiked.chr))
    cov <- standardize.cov(coverage(gr,weight = gr$score))
    return(mytts(cov,genes,profile = sub("_density.wig","",f)))
  },mc.cores = mc.cores)
  for(i in 1:length(l)){ pl8 <- rbind(pl8,l[[i]]) }
  rm(l)
  l <- mclapply(files,function(f){
    cat("ScaledGene: ",f,"\n")
    gr <- read.delim(f,header=F,comment.char = "t",sep = " ")
    gr <- with(gr,GRanges(V1,IRanges(V2,V3),"+",score=V4))
    gr <- keepSeqlevels(gr,c(chr.flt,spiked.chr))
    cov <- standardize.cov(coverage(gr,weight = gr$score))
    return(myscaledgene(cov,genes,profile = sub("_density.wig","",f)))
  },mc.cores = mc.cores)
  for(i in 1:length(l)){ pl9 <- rbind(pl9,l[[i]]) }
  rm(l)
  cat("\n")
  
  cat("Along Input divided samples\n")
  setwd(covdir)
  files <- Sys.glob("*_combIN.cov")
  pl10 <- c()  ## TSS
  pl11 <- c()  ## TTS
  pl12 <- c()  ## Scaled Gene
  l <- mclapply(files,function(f){
    cat("TSS: ",f,"\n")
    load(f)
    return(mytss(cov,genes,profile = sub("_combIN.cov","",f)))
  },mc.cores = mc.cores)
  for(i in 1:length(l)){ pl10 <- rbind(pl10,l[[i]]) }
  rm(l)
  l <- mclapply(files,function(f){
    cat("TTS: ",f,"\n")
    load(f)
    return(mytts(cov,genes,profile = sub("_combIN.cov","",f)))
  },mc.cores = mc.cores)
  for(i in 1:length(l)){ pl11 <- rbind(pl11,l[[i]]) }
  rm(l)
  l <- mclapply(files,function(f){
    cat("ScaledGene: ",f,"\n")
    load(f)
    return(myscaledgene(cov,genes,profile = sub("_combIN.cov","",f)))
  },mc.cores = mc.cores)
  for(i in 1:length(l)){ pl12 <- rbind(pl12,l[[i]]) }
  rm(l)
  cat("\n")
  
  
  ## Combined replicates
  cat("Along Input divided samples\n")
  setwd(covdir)
  files <- Sys.glob("*_combRelIN.cov")
  pl16 <- c()  ## TSS
  pl17 <- c()  ## TTS
  pl18 <- c()  ## Scaled Gene
  l <- mclapply(files,function(f){
    cat("TSS: ",f,"\n")
    load(f)
    return(mytss(cov,genes,profile = sub("_combRelIN.cov","",f)))
  },mc.cores = mc.cores)
  for(i in 1:length(l)){ pl16 <- rbind(pl16,l[[i]]) }
  rm(l)
  l <- mclapply(files,function(f){
    cat("TTS: ",f,"\n")
    load(f)
    return(mytts(cov,genes,profile = sub("_combRelIN.cov","",f)))
  },mc.cores = mc.cores)
  for(i in 1:length(l)){ pl17 <- rbind(pl17,l[[i]]) }
  rm(l)
  l <- mclapply(files,function(f){
    cat("ScaledGene: ",f,"\n")
    load(f)
    return(myscaledgene(cov,genes,profile = sub("_combRelIN.cov","",f)))
  },mc.cores = mc.cores)
  for(i in 1:length(l)){ pl18 <- rbind(pl18,l[[i]]) }
  rm(l)
  cat("\n")
  
  #save.image(file=paste("/n/data1/hms/dbmi/park/Dhawal/Project_nreim/MayNov17/gr/MetagenePlots_revised.RData",sep=""))
  save.image(file=paste("/n/data1/hms/dbmi/park/Dhawal/Project_nreim/MayNov17/gr/MetagenePlots_DF_binned_revised.RData",sep=""))
  
  
}


## Shift ratio calculation
if(F){
  if(F){
    library(GenomicRanges)
    getMatrix_IntervalSum <- function(covs, points) {
      
      m <- matrix(nrow=length(points$chr), ncol=1)
      
      for (i in 1:length(points$chr)) {
        if (points$start[i] > 0 & (points$end[i]) < length(covs[[as.character(points$chr[i])]]) & points$chr[i] %in% names(covs)) {
          m[i,1] <- sum(as.vector( window( covs[[as.character(points$chr[i])]], points$start[i], points$end[i] ) ))
        }
      }
      names(m) <- points$tracking_id
      return(m)
    }
    
    keys <- read.delim("SampleKey_relative.txt")
    load("/n/data1/hms/dbmi/park/Dhawal/Genomes/yeast_Annotations/ScSp_SGD.RevisedGenes.RData")
    #load("A:/work/yeast_Annotations/ScSp_RevisedGenes.RData")
    genes <- subset(genes,!genes$verification%in%c("transposable_element_gene","Verified|silenced_gene","Dubious",""))
    genes <- genes[genes$width>499,]
    genes <- as(genes,"GRanges")
    tss <- resize(genes,1,"start")
    tts <- resize(genes,1,"end")
    
    g5 <- flank(tss,500,start = F) ## first 500bp
    g3 <- flank(tts,500) ## last 500bp
    g5 <- as.data.frame(g5)
    g3 <- as.data.frame(g3)
    names(g3)[1] <- "chr"
    names(g5)[1] <- "chr"
    
    files<- Sys.glob("*_rel_IN.cov")
    l <- mclapply(files,function(f){
      cat(f,"\n")
      load(f)
      a <- getMatrix_IntervalSum(cov,g5)
      b <- getMatrix_IntervalSum(cov,g3)
      return(cbind(g5[,c(1,6,7)],five.prime=a,three.prime=b,profile=sub("_re;_IN.cov","",f)))
    },mc.cores = 4)
    save(l,file="ShiftRatios_RelativeNorm.RData")
    
    
    d <- c()
    for(i in 1:length(l)) d <- rbind(d,l[[i]])
    d$profile <- gsub("_rel_IN.cov","",d$profile)
    save(d,l,file="ShiftRatios_RelativeNorm.RData")
    
  }
  
}


### Nucleosome signal change
if(F){
  if(F){
    cat("Counting Signal in intervals for Shift ratios and nucleosomes\n")
    mc.cores=6
    library(GenomicRanges)
    setwd(grdir)
    files <- Sys.glob("*.gr")
    
    #load("/n/data1/hms/dbmi/park/Dhawal/Genomes/yeast_Annotations/ScSp_RevisedGenes.RData")
    load("A:/work/yeast_Annotations/ScSp_RevisedGenes.RData")
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
  
}


######################################################################################
### Natalia, ChIPSeq paper figures
######################################################################################

## Correlations for the paper
if(F){
  setwd("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr")
  rm(list = ls())
  load("RawCounts_SignalRegions_5p.RData") 
  m <- cf[,6:69]
  normF <- read.delim("Normalization.factors.NataliaMayNov2017_revised.txt",header=T)
  m1 <- matrix(NA,nrow = 31467,ncol = 64)
  for(i in 1:64){
    m1[,i] <- m[,i]*normF$normF[i]
  }
  m1 <- as.data.frame(m1)
  names(m1) <- names(m)
  
  h <- names(m)
  h <- gsub("RNAPII","Rpb1",h)
  h <- gsub("S2P","Ser2-P",h)
  h <- gsub("S5P","Ser5-P",h)
  h <- paste0("(ND) ",h)
  h <- gsub("[(]ND[])] Spn1.Depl_","(D) ",h)
  h <- gsub("HA.","",h)
  h <- gsub("[(]ND[])] Spn1_","(ND) ",h)
  names(m) <- names(m1) <- h
  
  pdf("Sample correlations _5p signal regions_final.pdf",width = 14,height = 14)
  {
    ## All samples (Raw, SpikeNorm)
    z1 <- round(cor(log2(m+1)),2)
    
    hmcols <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
    heatmap.2(z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.0,0.0),sepcolor = "black",
              key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
    hmcols <- palette(brewer.pal(n = 11, name = "RdYlBu"))
    heatmap.2(z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.0,0.0),sepcolor = "black",
              key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
    hmcols <- (colorRampPalette(c("white", "black"))(256))
    heatmap.2(z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.0,0.0),sepcolor = "black",
              key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
      ## post spikein
    z1 <- round(cor(log2(m1+1)),2)
    hmcols <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
    heatmap.2(z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.0,0.0),sepcolor = "black",
              key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
    hmcols <- palette(brewer.pal(n = 11, name = "RdYlBu"))
    heatmap.2(z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.0,0.0),sepcolor = "black",
              key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
    hmcols <- (colorRampPalette(c("white", "black"))(256))
    heatmap.2(z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.0,0.0),sepcolor = "black",
              key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
    
  }
  dev.off()
  
  
  
  ## All ChIP samples (Raw, SpikeNorm)
  n <- names(m) ## total 16 inputs + 48 ChIPs
  n <- n[-grep("Inp",n)]
  z1 <- round(cor(log2(m[ ,names(m)%in%n]+1)),2)
  hmcols <- palette(brewer.pal(n = 11, name = "RdYlBu"))
  heatmap.2(main = "ChIP-Seq correlations (pre-normalization)",z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.01,0.01),sepcolor = "black",colsep=0:ncol(z1),rowsep=0:nrow(z1),
            key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
  z1 <- round(cor(log2(m1[ ,names(m1)%in%n]+1)),2)
  hmcols <- palette(brewer.pal(n = 11, name = "RdYlBu"))
  heatmap.2(main = "ChIP-Seq correlations (post-normalization)",z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.01,0.01),sepcolor = "black",colsep=0:ncol(z1),rowsep=0:nrow(z1),
            key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
  
  ## All Histone samples (Raw, SpikeNorm)
  n <- names(m) ## total 16 inputs + 48 ChIPs
  n <- n[grep("H3|Set2",n)]
  z1 <- round(cor(log2(m[ ,names(m)%in%n]+1)),2)
  hmcols <- palette(brewer.pal(n = 11, name = "RdYlBu"))
  heatmap.2(main = "ChIP-Seq correlations (pre-normalization)",z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.01,0.01),sepcolor = "black",colsep=0:ncol(z1),rowsep=0:nrow(z1),
            key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
  z1 <- round(cor(log2(m1[ ,names(m1)%in%n]+1)),2)
  hmcols <- palette(brewer.pal(n = 11, name = "RdYlBu"))
  heatmap.2(main = "ChIP-Seq correlations (post-normalization)",z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.01,0.01),sepcolor = "black",colsep=0:ncol(z1),rowsep=0:nrow(z1),
            key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
  
  
  ## All Histone samples (Raw, SpikeNorm)
  n <- names(m) ## total 16 inputs + 48 ChIPs
  n <- n[grep("RNAP|S2P|S5P|Spt6",n)]
  z1 <- round(cor(log2(m[ ,names(m)%in%n]+1)),2)
  hmcols <- palette(brewer.pal(n = 11, name = "RdYlBu"))
  heatmap.2(main = "ChIP-Seq correlations (pre-normalization)",z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.01,0.01),sepcolor = "black",colsep=0:ncol(z1),rowsep=0:nrow(z1),
            key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
  z1 <- round(cor(log2(m1[ ,names(m1)%in%n]+1)),2)
  hmcols <- palette(brewer.pal(n = 11, name = "RdYlBu"))
  heatmap.2(main = "ChIP-Seq correlations (post-normalization)",z1,notecol="black",notecex =1,density.info="none",trace="none",margins =c(15,15),col=hmcols,sepwidth = c(0.01,0.01),sepcolor = "black",colsep=0:ncol(z1),rowsep=0:nrow(z1),
            key.title = "",key.xlab = "Correlation Coefficient",key.ylab = "") 
  
  
  
}

## Multigene heatmap and metagenes for gene categories
if(F){
  rm(list=ls())
  load("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/RawCounts_Signal_genes.RData")
  cf <- cf[cf$species=="Scer",]
  cf <- cf[,c(6,65:68)]
  cf$RNAPII <- rowSums(cf[,2:5])
  cf <- cf[,c(1,6)]
  
  load(file="A:/work/WinstonLab/Natalia_Reim/Analysis_Feb18/DE_Genes_Filtered.RData")
  
  load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenes.RData")
  genes <- subset(genes,genes$verification=="Verified" & genes$species=="Scer" & genes$chr!="chrM")
  genes <- merge(genes,cf,by="tracking_id",all.x=T)
  genes <- merge(genes,de.genes[,c("tracking_id","log2FoldChange")],by="tracking_id",all.x=T)
  genes$RNAPII <- ifelse(is.na(genes$RNAPII),0,genes$RNAPII)
  genes$log2FoldChange <- ifelse(is.na(genes$log2FoldChange),0,genes$log2FoldChange)
  
  genes <- genes[order(genes$width,decreasing = F),]
  de.genes <- de.genes[de.genes$log2FoldChange<0 & de.genes$padj<0.05,]
  
  genes <- genes[order(genes$width,decreasing = F),]
  
  rp <- read.delim("A:/work/WinstonLab/Natalia_Reim/Analysis_Feb18/All_RP_proteincodingVErifiedtranscripts.txt",header=T)
  
  load("A:/work/WinstonLab/Natalia_Reim/Analysis_Feb18/intron.stats.5prime.q50.new.RData")
  introns <- unique(introns)
  introns <- subset(introns,introns$reference%in%genes[genes$verification=="Verified" & genes$chr!="chrM",]$tracking_id)
  introns$width <- abs(introns$start-introns$end)
  introns <- introns[introns$width>=50,]
  rm(intron.stats)
  names(introns)[6] <- "tracking_id"
  
  opndpn <- read.delim("A:/work/yeast_Annotations/opn_dpn_genes.txt",header = T)
  opn <- data.frame(tracking_id=opndpn$OPEN_Zuagg)
  dpn <- data.frame(tracking_id=opndpn$CLOSED_Zaugg)
  opn <- as.data.frame(opn[opn$tracking_id!="",])
  dpn <- as.data.frame(dpn[dpn$tracking_id!="",])
  names(opn)[1] <- "tracking_id"
  names(dpn)[1] <- "tracking_id"
  opn <- subset(opn,opn$tracking_id%in%genes$tracking_id)
  dpn <- subset(dpn,dpn$tracking_id%in%genes$tracking_id)
  
  h <- quantile(genes$RNAPII,c(0.33,0.66))
  genes$Pol2OccGp <- 1
  genes$Pol2OccGp <- ifelse(genes$RNAPII>h[2],3,genes$Pol2OccGp)
  genes$Pol2OccGp <- ifelse(genes$RNAPII<h[2] & genes$RNAPII>h[1],2,genes$Pol2OccGp)
  
  genes$lengp <- 1
  genes$lengp <- ifelse(genes$width>2000,4,genes$lengp)
  genes$lengp <- ifelse(genes$width<2000 & genes$width>1000,3,genes$lengp)
  genes$lengp <- ifelse(genes$width<1000 & genes$width>500,2,genes$lengp)
  
  #save.image(file="Gene groups.RData")
  
  rm(list=ls())
  load(file="Gene groups.RData")
  protein<- c("H3 (round1)","H3","H3K4me3","H3K36me2","H3K36me3","RNAPII (round1)","RNAPII","S2P","S5P","Set2","Spn1","Spt6")
  
  ## Heatmaps
  res <- read.delim("SpikeinNormalization_factors_revised.txt",header = T)
  myheatplotfun <- function(l,genes,g,outfile="temp.jpeg",file=c("H3K36me3"),normF=c(1,1)){
    max.gene.len=501
    nrow= length(genes$tracking_id)
    a1 <- matrix(NA,ncol=max.gene.len,nrow = nrow)  ## Depl, TSS
    b1 <- matrix(NA,ncol=max.gene.len,nrow = nrow)  ## Depl, CPS
    a2 <- matrix(NA,ncol=max.gene.len,nrow = nrow)  ## Non-depl, TSS
    b2 <- matrix(NA,ncol=max.gene.len,nrow = nrow)  ## Non-depl, CPS
    for(i in 1:length(genes$tracking_id)){
      v = l[[1]][i,]+l[[2]][i,]
      x <- tail(v[51:(50+round(genes$width[i]/10))],400)
      x <- c(v[1:50],x,v[1524:1574])
      b1[i,] <- c(rep(0,(max.gene.len-length(x))),x)
      
      v = l[[1]][i,]+l[[2]][i,]
      x <- head(v[51:(50+round(genes$width[i]/10))],400)
      x <- c(v[1:50],x,v[1524:1574])
      a1[i,] <- c(x,rep(0,(max.gene.len-length(x))))
      
      
      v = l[[3]][i,]+l[[4]][i,]
      x <- tail(v[51:(50+round(genes$width[i]/10))],400)
      x <- c(v[1:50],x,v[1524:1574])
      b2[i,] <- c(rep(0,(max.gene.len-length(x))),x)
      
      v = l[[3]][i,]+l[[4]][i,]
      x <- head(v[51:(50+round(genes$width[i]/10))],400)
      x <- c(v[1:50],x,v[1524:1574])
      a2[i,] <- c(x,rep(0,(max.gene.len-length(x))))
    }
    rownames(a1) <- rownames(a2) <-  rownames(b1) <- rownames(b2) <-  rownames(l[[1]])
    
    
    genes <- subset(genes,genes$tracking_id%in%g$tracking_id)
    a1 <- subset(a1,rownames(a1)%in%genes$tracking_id)
    a2 <- subset(a2,rownames(a2)%in%genes$tracking_id)
    b1 <- subset(b1,rownames(b1)%in%genes$tracking_id)
    b2 <- subset(b2,rownames(b2)%in%genes$tracking_id)
    
    max.gene.len = round(max(genes$width)/10)
    if(max(genes$width)%%10 < 5){
      max.gene.len=max.gene.len+1
    }
    if(max.gene.len<400){
      gene.len.cut = max.gene.len*10
      x <- seq(-500,(gene.len.cut+500),10)
      a1 <- a1[,1:(max.gene.len+100)]
      a2 <- a2[,1:(max.gene.len+100)]
      b1 <- b1[,(400-max.gene.len+2):501]
      b2 <- b2[,(400-max.gene.len+2):501]
      at = seq(0,max.gene.len*10,length.out = 4)
      label = rep("",length(at))
      label[1] = "TSS"
      label[length(label)] = "CPS"
      at = c(-500,at,at[length(at)]+500)
      label = c("",label,"")
    }else{
      gene.len.cut = 400*10
      x <- seq(-500,(gene.len.cut+500),10)
      at = c(-500,0,1000,2000,3000,4000,4500)
      label=c("","TSS","","","","CPS","")
    }
    
    
    scale_row <- function(m){
      a3 <- t(apply(m,1,function(x){
        x[x==0] <- NA
        return((x-mean(x,na.rm = T))/sd(x,na.rm = T))
      }))
      return(a3)
    }
    get_divsion <- function(m,n){
      m <- m*100+1
      n <- n*100+1
      m <- (m-n)/(m+n)
      return(m)
    }
    
    a1 <- a1*normF[1]
    b1 <- b1*normF[1]
    a2 <- a2*normF[2]
    b2 <- b2*normF[2]
    
    a <- get_divsion(a1,a2) ## TSS centered
    b <- get_divsion(b1,b2) ## TTS centered
    
    a1 <- scale_row(a1)
    a2 <- scale_row(a2)
    b1 <- scale_row(b1)
    b2 <- scale_row(b2)
    
    zlim1=range(a1,a2,b1,b2,na.rm = T)
    zlim2=range(a,b,na.rm = T)
    

    x11()
    layout.show(layout(matrix(c(1:12),ncol=3,byrow = T),heights = c(5,1,5,1,5,1)))
    
    # 1,1 and 1,2
    par(mar=c(3,4,4,1))
    h <- get.unscaledgene.heatplot.list(r=a1,mid.point = 0,genes = genes,zlim = zlim1,cols = c("white","gray90","black"),x = x,gene.len.cut =gene.len.cut )
    image(h$x,1:nrow(h$z$z),t(h$z$z),zlim=h$zlim,col=h$z$col.z,axes=FALSE,xlab=" ",ylab=" ",main="(I) Depleted",cex.main=4);box();
    axis(1,at=at,labels = label,cex.axis=3,padj = 0.5)
    lines(h$border.lines[[1]],1:nrow(h$z$z),lty=3,lwd=1,col="green4",xlab="",ylab="")
    lines(h$border.lines[[2]],rev(1:nrow(h$z$z)),lty=3,lwd=1,col="green4",xlab="",ylab="")
    h <- get.unscaledgene.heatplot.list(r=a2,mid.point = 0,genes = genes,zlim = zlim1,cols = c("white","gray90","black"),x = x,gene.len.cut =gene.len.cut)
    image(h$x,1:nrow(h$z$z),t(h$z$z),zlim=h$zlim,col=h$z$col.z,axes=FALSE,xlab=" ",ylab=" ",main="(II) Non-depleted",cex.main=4);box();
    axis(1,at=at,labels = label,cex.axis=3,padj = 0.5)
    lines(h$border.lines[[1]],1:nrow(h$z$z),lty=3,lwd=1,col="green4",xlab="",ylab="")
    lines(h$border.lines[[2]],rev(1:nrow(h$z$z)),lty=3,lwd=1,col="green4",xlab="",ylab="")
    h <- get.unscaledgene.heatplot.list(r=a,mid.point = 0,genes = genes,zlim = zlim2,cols = c("navyblue","white","red4"),x = x,gene.len.cut =gene.len.cut)
    image(h$x,1:nrow(h$z$z),t(h$z$z),zlim=h$zlim,col=h$z$col.z,axes=FALSE,xlab=" ",ylab=" ",main="Difference, (I-II)/(I+II)",cex.main=4);box();
    axis(1,at=at,labels = label,cex.axis=3,padj = 0.5)
    lines(h$border.lines[[1]],1:nrow(h$z$z),lty=3,lwd=1,col="green4",xlab="",ylab="")
    lines(h$border.lines[[2]],rev(1:nrow(h$z$z)),lty=3,lwd=1,col="green4",xlab="",ylab="")
    
    par(mar=c(3,3,1,3))
    h <- get.unscaledgene.heatplot.list(r=a1,mid.point = 0,genes = genes,zlim = zlim1,cols = c("white","gray90","black"),x = x,gene.len.cut =gene.len.cut)
    image(h$z$scale$x,h$z$scale$y,h$z$scale$z,xlab="",ylab="",col=h$z$col.scale,axes=F)
    axis(1,at=c(min(h$z$scale$x),1,mean(h$z$scale$x),max(h$z$scale$x)),round(c(min(h$z$scale$x),1,mean(h$z$scale$x),max(h$z$scale$x))),cex.axis=3)
    image(h$z$scale$x,h$z$scale$y,h$z$scale$z,xlab="",ylab="",col=h$z$col.scale,axes=F)
    axis(1,at=c(min(h$z$scale$x),1,mean(h$z$scale$x),max(h$z$scale$x)),round(c(min(h$z$scale$x),1,mean(h$z$scale$x),max(h$z$scale$x))),cex.axis=3)
    h <- get.unscaledgene.heatplot.list(r=a,mid.point = 0,genes = genes,zlim = zlim2,cols = c("navyblue","white","red4"))
    image(h$z$scale$x,h$z$scale$y,h$z$scale$z,xlab="",ylab="",col=h$z$col.scale,axes=F)
    axis(1,at=c(min(h$z$scale$x),1,mean(h$z$scale$x),max(h$z$scale$x)),round(c(min(h$z$scale$x),1,mean(h$z$scale$x),max(h$z$scale$x))),cex.axis=3)
    
    
    par(mar=c(3,4,4,1))
    h <- get.unscaledgene.heatplot.list(r=b1,mid.point = 0,genes = genes,zlim = zlim1,cols = c("white","gray90","black"),align = "CPS",x = x,gene.len.cut =gene.len.cut)
    image(h$x,1:nrow(h$z$z),t(h$z$z),zlim=h$zlim,col=h$z$col.z,axes=FALSE,xlab=" ",ylab=" ",main="(I) Depleted",cex.main=4);box();
    axis(1,at=at,labels = label,cex.axis=3,padj = 0.5)
    lines(h$border.lines[[1]],1:nrow(h$z$z),lty=3,lwd=1,col="green4",xlab="",ylab="")
    lines(h$border.lines[[2]],rev(1:nrow(h$z$z)),lty=3,lwd=1,col="green4",xlab="",ylab="")
    h <- get.unscaledgene.heatplot.list(r=b2,mid.point = 0,genes = genes,zlim = zlim1,cols = c("white","gray90","black"),align = "CPS",x = x,gene.len.cut =gene.len.cut)
    image(h$x,1:nrow(h$z$z),t(h$z$z),zlim=h$zlim,col=h$z$col.z,axes=FALSE,xlab=" ",ylab=" ",main="(II) Non-depleted",cex.main=4);box();
    axis(1,at=at,labels = label,cex.axis=3,padj = 0.5)
    lines(h$border.lines[[1]],1:nrow(h$z$z),lty=3,lwd=1,col="green4",xlab="",ylab="")
    lines(h$border.lines[[2]],rev(1:nrow(h$z$z)),lty=3,lwd=1,col="green4",xlab="",ylab="")
    h <- get.unscaledgene.heatplot.list(r=b,mid.point = 0,genes = genes,zlim = zlim2,cols = c("navyblue","white","red4"),align = "CPS",x = x,gene.len.cut =gene.len.cut)
    image(h$x,1:nrow(h$z$z),t(h$z$z),zlim=h$zlim,col=h$z$col.z,axes=FALSE,xlab=" ",ylab=" ",main="Difference, (I-II)/(I+II)",cex.main=4);box();
    axis(1,at=at,labels = label,cex.axis=3,padj = 0.5)
    lines(h$border.lines[[1]],1:nrow(h$z$z),lty=3,lwd=1,col="green4",xlab="",ylab="")
    lines(h$border.lines[[2]],rev(1:nrow(h$z$z)),lty=3,lwd=1,col="green4",xlab="",ylab="")
    
    par(mar=c(3,3,1,3))
    h <- get.unscaledgene.heatplot.list(r=b1,mid.point = 0,genes = genes,zlim = zlim1,cols = c("white","gray90","black"),x = x,gene.len.cut =gene.len.cut)
    image(h$z$scale$x,h$z$scale$y,h$z$scale$z,xlab="",ylab="",col=h$z$col.scale,axes=F)
    axis(1,at=c(min(h$z$scale$x),1,mean(h$z$scale$x),max(h$z$scale$x)),round(c(min(h$z$scale$x),1,mean(h$z$scale$x),max(h$z$scale$x))),cex.axis=3)
    image(h$z$scale$x,h$z$scale$y,h$z$scale$z,xlab="",ylab="",col=h$z$col.scale,axes=F)
    axis(1,at=c(min(h$z$scale$x),1,mean(h$z$scale$x),max(h$z$scale$x)),round(c(min(h$z$scale$x),1,mean(h$z$scale$x),max(h$z$scale$x))),cex.axis=3)
    h <- get.unscaledgene.heatplot.list(r=b,mid.point = 0,genes = genes,zlim = zlim2,cols = c("navyblue","white","red4"))
    image(h$z$scale$x,h$z$scale$y,h$z$scale$z,xlab="",ylab="",col=h$z$col.scale,axes=F)
    axis(1,at=c(min(h$z$scale$x),1,mean(h$z$scale$x),max(h$z$scale$x)),round(c(min(h$z$scale$x),1,mean(h$z$scale$x),max(h$z$scale$x))),cex.axis=2.2)
    
    dev.copy(device=jpeg,file=paste0(outfile),width=1500,height=1800);dev.off();
    
  }
  for(prot in protein){
    f = paste("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_",prot,".RData",sep="")
    if(file.exists(f)){
      cat(prot,"\n")
      cat(res[res$Factor==prot,]$normF,"\n")
      load(f)
      #myheatplotfun(l,genes,g=genes,outfile=paste(prot,"_Allgenes_revised.jpeg",sep=""), normF = res[res$Factor==prot,]$normF)
      #myheatplotfun(l,genes,g=rp,outfile=paste(prot,"_RPgenes_revised.jpeg",sep=""), normF = res[res$Factor==prot,]$normF)
      #myheatplotfun(l,genes,g=introns,outfile=paste(prot,"_Introngenes_revised.jpeg",sep=""), normF = res[res$Factor==prot,]$normF)
      myheatplotfun(l,genes,g=de.genes,outfile=paste(prot,"_DEgenes_revised.jpeg",sep=""), normF = res[res$Factor==prot,]$normF)
      #myheatplotfun(l,genes,g=opn,outfile=paste(prot,"_OPEN_ZuaggGenes_revised.jpeg",sep=""), normF = res[res$Factor==prot,]$normF)
      #myheatplotfun(l,genes,g=dpn,outfile=paste(prot,"_CLOSED_ZuaggGenes_revised.jpeg",sep=""), normF = res[res$Factor==prot,]$normF)
    }
  }
  
  ## Metagne plots
  ## 1
  res <- read.delim("Normalization.factors.NataliaMayNov2017_revised.txt",header = T)
  vector.resizing <- function(x,final.len){
    y <- vector()
    len <- length(x)
    y <-spline(1:len,x,n=final.len)$y
    return(y)
  }
  myplot <- function(pl1){
    require(ggplot2)
    q <- ggplot(pl1, aes(x=pos,y=signal,ymax=signal+se, ymin=signal-se,col=condition))
    q <- q + geom_vline(xintercept = c(50,150),colour="gray50", linetype = "dashed")
    q <- q + geom_linerange(col='gray') + geom_line() #+ ylim(-0.03,0.15) #ylim(-0.01,0.045) #ylim(-0.04,0.1)
    q <- q + scale_color_manual(values = c("Depleted-1" = "pink","Depleted-2" = "red2",
                                           "Non-depleted-1" = "lightblue","Non-depleted-2" = "slateblue2"))
    q <- q + facet_wrap(~profile, ncol=4,scales = "free")
    q <- q + scale_x_continuous(breaks = c(0,50,150,199), labels=c("","TSS","CPS","")) +ylab("average signal") + xlab("") 
    q <- q + theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = NA),axis.line = element_line(size = 0.7, colour = "gray50"))
    q <- q + theme(axis.text.x = element_text(colour = 'black',size=14),
                   axis.text.y = element_text(colour = 'black',size=14),
                   axis.title = element_text(colour = 'black',size=15,face="bold"),
                   strip.text.x = element_text(colour = 'black',size=14,face="bold"))
    q <- q + theme(legend.text =element_text(colour = 'black',size=14),
                   legend.title = element_text(colour = 'black',size=14,face="bold"))
    q <- q + guides(color=guide_legend(title="Spn1 levels"))
    return(q)
    
  }
  mymetagenefun <- function(l,genes,g,normF=c(1,1,1,1)){
    depl1 <- matrix(NA,ncol=199,nrow = 5134)
    depl2 <- matrix(NA,ncol=199,nrow = 5134)
    nondepl1 <- matrix(NA,ncol=199,nrow = 5134)
    nondepl2 <- matrix(NA,ncol=199,nrow = 5134)
    for(i in 1:length(genes$tracking_id)){
      j = round(genes[i,]$width/10)
      v <- l[[1]][i,51:(49+j)]
      v <- vector.resizing(v,100)
      v <- c(l[[1]][i,1:50],v,l[[1]][i,1525:1573])
      #plot(v)
      depl1[i,] <- v
      v <- l[[2]][i,51:(49+j)]
      v <- vector.resizing(v,100)
      v <- c(l[[2]][i,1:50],v,l[[2]][i,1525:1573])
      depl2[i,] <- v
      
      v <- l[[3]][i,51:(49+j)]
      v <- vector.resizing(v,100)
      v <- c(l[[3]][i,1:50],v,l[[3]][i,1525:1573])
      nondepl1[i,] <- v
      v <- l[[4]][i,51:(49+j)]
      v <- vector.resizing(v,100)
      v <- c(l[[4]][i,1:50],v,l[[4]][i,1525:1573])
      nondepl2[i,] <- v
    }
    rownames(depl1) <- rownames(depl2) <- rownames(nondepl1) <- rownames(nondepl2) <-  rownames(l[[1]])
    
    depl1 <- subset(depl1,rownames(depl1)%in%g$tracking_id)
    depl2 <- subset(depl2,rownames(depl2)%in%g$tracking_id)
    nondepl1 <- subset(nondepl1,rownames(nondepl1)%in%g$tracking_id)
    nondepl2 <- subset(nondepl2,rownames(nondepl2)%in%g$tracking_id)
    
    get.avg <- function(m,nF){
      df_m <- cbind(1:ncol(m),
                    as.data.frame(apply(m,2,mean,na.rm=TRUE)),
                    as.data.frame(apply(m,2,sd,na.rm=TRUE)), 
                    as.data.frame(apply(m,2, function(x) sd(x,na.rm=TRUE)/sqrt(length(x))))
      )
      names(df_m) <- c("pos","signal", "sd","se")
      df_m$signal <- df_m$signal*nF
      df_m$sd <- df_m$sd*nF
      return(df_m)
    }
    pl <- rbind(cbind(get.avg(depl1,normF[1]),condition="Depleted-1"),
                cbind(get.avg(depl2,normF[2]),condition="Depleted-2"),
                cbind(get.avg(nondepl1,normF[3]),condition="Non-depleted-1"),
                cbind(get.avg(nondepl2,normF[4]),condition="Non-depleted-2"))
    return(pl)
  }
  pl <- c()
  for(prot in protein){
    f = paste("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_",prot,".RData",sep="")
    if(file.exists(f)){
      cat(prot,"\n")
      cat(res[res$Factor==prot,]$normF,"\n\n")
      load(f)
      pl <- rbind(pl,
                  cbind(profile=prot,genetype = "all genes", mymetagenefun(l,genes,g=genes,normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "RP genes", mymetagenefun(l,genes,g=rp,normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "Intron genes", mymetagenefun(l,genes,g=introns, normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "DE genes", mymetagenefun(l,genes,g=de.genes, normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "OPEN genes", mymetagenefun(l,genes,g=opn, normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "CLOSED genes", mymetagenefun(l,genes,g=dpn, normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "High Pol2", mymetagenefun(l,genes,g=genes[genes$Pol2OccGp==3,],normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "Medium Pol2", mymetagenefun(l,genes,g=genes[genes$Pol2OccGp==2,],normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "Low Pol2", mymetagenefun(l,genes,g=genes[genes$Pol2OccGp==1,], normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "<0.5kb", mymetagenefun(l,genes,g=genes[genes$lengp==1,], normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "0.5kb-1kb", mymetagenefun(l,genes,g=genes[genes$lengp==2,], normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "1kb-2kb", mymetagenefun(l,genes,g=genes[genes$lengp==3,], normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = ">2kb", mymetagenefun(l,genes,g=genes[genes$lengp==4,], normF = res[res$Factor==prot,]$normF))
      )
    }
  }
  
  #pdf("MetagenePlots_SpikeNorm_AllFactors.pdf",width = 15,height = 9)
  pl$signal <- pl$signal/2
  pdf("MetagenePlots_SpikeNorm_AllFactors_avg.pdf",width = 15,height = 9)
  myplot(pl[pl$genetype=="all genes",]) + ggtitle("All verified protein coding genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype=="RP genes",]) + ggtitle("RP[PSL] genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype=="Intron genes",]) + ggtitle("Intron containing genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype=="DE genes",]) + ggtitle("Downregulated genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype=="OPEN genes",]) + ggtitle("OPEN genes (Zuagg et al)")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype=="CLOSED genes",]) + ggtitle("CLOSED genes (Zuagg et al)")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype=="High Pol2",]) + ggtitle("High Pol2 occupancy genes (n=1746)")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype=="Medium Pol2",]) + ggtitle("Medium Pol2 occupancy genes (n=1694)")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype=="Low Pol2",]) + ggtitle("Low Pol2 occupancy genes (n=1694)")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype=="<0.5kb",]) + ggtitle("Gene length <0.5kb (n=193)")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype=="0.5kb-1kb",]) + ggtitle("Gene length 0.5kb-1kb (n=1166)")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype=="1kb-2kb",]) + ggtitle("Gene length 1kb-2kb (n=2244)")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype==">2kb",]) + ggtitle("Gene length >2kb (n=1531)")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  dev.off()
  

  ## 3. Metagene plot for merged replicates
  res <- read.delim("Normalization.factors.merged.txt",header = T)
  myplot <- function(pl1){
    require(ggplot2)
    q <- ggplot(pl1, aes(x=pos,y=signal,ymax=signal+se, ymin=signal-se,col=condition))
    q <- q + geom_vline(xintercept = c(50,150),colour="gray50", linetype = "dashed")
    q <- q + geom_linerange(col='gray') + geom_line() #+ ylim(-0.03,0.15) #ylim(-0.01,0.045) #ylim(-0.04,0.1)
    q <- q + scale_color_manual(values = c("Depleted-1" = "pink","Depleted-2" = "red2",
                                           "Non-depleted-1" = "lightblue","Non-depleted-2" = "slateblue2",
                                           "Depleted" = "red","Non-depleted" = "slateblue"))
    q <- q + facet_wrap(~profile, ncol=4,scales = "free")
    q <- q + scale_x_continuous(breaks = c(0,50,150,199), labels=c("","TSS","CPS","")) +ylab("average signal") + xlab("") 
    q <- q + theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = NA),axis.line = element_line(size = 0.7, colour = "gray50"))
    q <- q + theme(axis.text.x = element_text(colour = 'black',size=14),
                   axis.text.y = element_text(colour = 'black',size=14),
                   axis.title = element_text(colour = 'black',size=15,face="bold"),
                   strip.text.x = element_text(colour = 'black',size=14,face="bold"))
    q <- q + theme(legend.text =element_text(colour = 'black',size=14),
                   legend.title = element_text(colour = 'black',size=14,face="bold"))
    q <- q + guides(color=guide_legend(title="Spn1 levels"))
    return(q)
    
  }
  mymetagenefun <- function(l,genes,g,normF=c(1,1)){
    depl1 <- matrix(NA,ncol=199,nrow = 5134)
    nondepl1 <- matrix(NA,ncol=199,nrow = 5134)
    for(i in 1:length(genes$tracking_id)){
      j = round(genes[i,]$width/10)
      v <- l[[1]][i,51:(49+j)] + l[[2]][i,51:(49+j)]
      v <- vector.resizing(v,100)
      v <- c(l[[1]][i,1:50]+l[[2]][i,1:50],
             v,
             l[[1]][i,1525:1573]+l[[2]][i,1525:1573]
             )
      depl1[i,] <- v

      v <- l[[3]][i,51:(49+j)] + l[[4]][i,51:(49+j)]
      v <- vector.resizing(v,100)
      v <- c(l[[3]][i,1:50]+l[[4]][i,1:50],
             v,
             l[[3]][i,1525:1573]+l[[4]][i,1525:1573]
             )
      nondepl1[i,] <- v
    }
    rownames(depl1) <- rownames(nondepl1) <- rownames(l[[1]])
    
    depl1 <- subset(depl1,rownames(depl1)%in%g$tracking_id)
    nondepl1 <- subset(nondepl1,rownames(nondepl1)%in%g$tracking_id)
    
    get.avg <- function(m,nF){
      df_m <- cbind(1:ncol(m),
                    as.data.frame(apply(m,2,mean,na.rm=TRUE)),
                    as.data.frame(apply(m,2,sd,na.rm=TRUE)), 
                    as.data.frame(apply(m,2, function(x) sd(x,na.rm=TRUE)/sqrt(length(x))))
      )
      names(df_m) <- c("pos","signal", "sd","se")
      df_m$signal <- df_m$signal*nF
      df_m$sd <- df_m$sd*nF
      return(df_m)
    }
    pl <- rbind(cbind(get.avg(depl1,normF[1]),condition="Depleted"),
                cbind(get.avg(nondepl1,normF[2]),condition="Non-depleted")
                )
    return(pl)
  }
  pl <- c()
  for(prot in protein){
    f = paste("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_",prot,".RData",sep="")
    if(file.exists(f)){
      cat(prot,"\n")
      cat(res[res$Factor==prot,]$normF,"\n\n")
      load(f)
      pl <- rbind(pl,
                  cbind(profile=prot,genetype = "all genes", mymetagenefun(l,genes,g=genes,normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "RP genes", mymetagenefun(l,genes,g=rp,normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "Intron genes", mymetagenefun(l,genes,g=introns, normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "DE genes", mymetagenefun(l,genes,g=de.genes, normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "OPEN genes", mymetagenefun(l,genes,g=opn, normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "CLOSED genes", mymetagenefun(l,genes,g=dpn, normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "High Pol2", mymetagenefun(l,genes,g=genes[genes$Pol2OccGp==3,],normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "Medium Pol2", mymetagenefun(l,genes,g=genes[genes$Pol2OccGp==2,],normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "Low Pol2", mymetagenefun(l,genes,g=genes[genes$Pol2OccGp==1,], normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "<0.5kb", mymetagenefun(l,genes,g=genes[genes$lengp==1,], normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "0.5kb-1kb", mymetagenefun(l,genes,g=genes[genes$lengp==2,], normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "1kb-2kb", mymetagenefun(l,genes,g=genes[genes$lengp==3,], normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = ">2kb", mymetagenefun(l,genes,g=genes[genes$lengp==4,], normF = res[res$Factor==prot,]$normF))
      )
    }
  }
  
  #pdf("MetagenePlots_SpikeNorm_AllFactors_Merged.pdf",width = 16,height = 8)
  pl$signal <- pl$signal/2
  pdf("MetagenePlots_SpikeNorm_AllFactors_Merged_avg.pdf",width = 16,height = 8)
  myplot(pl[pl$genetype=="all genes",]) + ggtitle("All verified protein coding genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype=="RP genes",]) + ggtitle("RP[PSL] genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype=="Intron genes",]) + ggtitle("Intron containing genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype=="DE genes",]) + ggtitle("Downregulated genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype=="OPEN genes",]) + ggtitle("OPEN genes (Zuagg et al)")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype=="CLOSED genes",]) + ggtitle("CLOSED genes (Zuagg et al)")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype=="High Pol2",]) + ggtitle("High Pol2 occupancy genes (n=1746)")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype=="Medium Pol2",]) + ggtitle("Medium Pol2 occupancy genes (n=1694)")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype=="Low Pol2",]) + ggtitle("Low Pol2 occupancy genes (n=1694)")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype=="<0.5kb",]) + ggtitle("Gene length <0.5kb (n=193)")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype=="0.5kb-1kb",]) + ggtitle("Gene length 0.5kb-1kb (n=1166)")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype=="1kb-2kb",]) + ggtitle("Gene length 1kb-2kb (n=2244)")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype==">2kb",]) + ggtitle("Gene length >2kb (n=1531)")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  dev.off()
  save(pl,file="MergedMetagenes_ChIPSpike_final.RData")
  
  
}

## relative signal in heatmap and metagene plots
if(F){
  getmatrices <- function(l,genes,g,normF=c(1,1)){
    max.gene.len=501
    nrow= length(genes$tracking_id)
    a1 <- matrix(NA,ncol=max.gene.len,nrow = nrow)  ## Depl, TSS
    b1 <- matrix(NA,ncol=max.gene.len,nrow = nrow)  ## Depl, CPS
    a2 <- matrix(NA,ncol=max.gene.len,nrow = nrow)  ## Non-depl, TSS
    b2 <- matrix(NA,ncol=max.gene.len,nrow = nrow)  ## Non-depl, CPS
    for(i in 1:length(genes$tracking_id)){
      v = l[[1]][i,]+l[[2]][i,]
      x <- tail(v[51:(50+round(genes$width[i]/10))],400)
      x <- c(v[1:50],x,v[1524:1574])
      b1[i,] <- c(rep(0,(max.gene.len-length(x))),x)
      
      v = l[[1]][i,]+l[[2]][i,]
      x <- head(v[51:(50+round(genes$width[i]/10))],400)
      x <- c(v[1:50],x,v[1524:1574])
      a1[i,] <- c(x,rep(0,(max.gene.len-length(x))))
      
      
      v = l[[3]][i,]+l[[4]][i,]
      x <- tail(v[51:(50+round(genes$width[i]/10))],400)
      x <- c(v[1:50],x,v[1524:1574])
      b2[i,] <- c(rep(0,(max.gene.len-length(x))),x)
      
      v = l[[3]][i,]+l[[4]][i,]
      x <- head(v[51:(50+round(genes$width[i]/10))],400)
      x <- c(v[1:50],x,v[1524:1574])
      a2[i,] <- c(x,rep(0,(max.gene.len-length(x))))
    }
    rownames(a1) <- rownames(a2) <-  rownames(b1) <- rownames(b2) <-  rownames(l[[1]])
    genes <- subset(genes,genes$tracking_id%in%g$tracking_id)
    a1 <- subset(a1,rownames(a1)%in%genes$tracking_id)
    a2 <- subset(a2,rownames(a2)%in%genes$tracking_id)
    b1 <- subset(b1,rownames(b1)%in%genes$tracking_id)
    b2 <- subset(b2,rownames(b2)%in%genes$tracking_id)
    a1 <- a1*normF[1]
    b1 <- b1*normF[1]
    a2 <- a2*normF[2]
    b2 <- b2*normF[2]
    return(list(a1=a1,b1=b1,a2=a2,b2=b2))
  }
  myheatplotfun <- function(l1,l2,genes,g,outfile="temp.jpeg",normF1=c(1,1),normF2=c(1,1)){
    max.gene.len=501
    
    p <- getmatrices(l1,genes,g,normF1)
    ax1 <- p[[1]]
    bx1 <- p[[2]]
    ax2 <- p[[3]]
    bx2 <- p[[4]]
    p <- getmatrices(l2,genes,g,normF2)
    ay1 <- p[[1]]
    by1 <- p[[2]]
    ay2 <- p[[3]]
    by2 <- p[[4]]
    
    a1 = (ax1)/(ay1)
    a1[!is.finite(a1)] <-NA
    a2 = (ax2)/(ay2)
    a2[!is.finite(a2)] <-NA
    b1 = (bx1)/(by1)
    b1[!is.finite(b1)] <-NA
    b2 = (bx2)/(by2)
    b2[!is.finite(b2)] <-NA
    rm(p,ax1,ay1,ax2,ay2,bx1,by1,bx2,by2)
    
    genes <- subset(genes,genes$tracking_id%in%g$tracking_id)
    
    max.gene.len = round(max(genes$width)/10)
    if(max(genes$width)%%10 < 5){
      max.gene.len=max.gene.len+1
    }
    if(max.gene.len<400){
      gene.len.cut = max.gene.len*10
      x <- seq(-500,(gene.len.cut+500),10)
      a1 <- a1[,1:(max.gene.len+100)]
      a2 <- a2[,1:(max.gene.len+100)]
      b1 <- b1[,(400-max.gene.len+2):501]
      b2 <- b2[,(400-max.gene.len+2):501]
      at = seq(0,max.gene.len*10,length.out = 4)
      label = rep("",length(at))
      label[1] = "TSS"
      label[length(label)] = "CPS"
      at = c(-500,at,at[length(at)]+500)
      label = c("",label,"")
    }else{
      gene.len.cut = 400*10
      x <- seq(-500,(gene.len.cut+500),10)
      at = c(-500,0,1000,2000,3000,4000,4500)
      label=c("","TSS","","","","CPS","")
    }
    
    scale_row <- function(m){
      a3 <- t(apply(m,1,function(x){
        x[x==0] <- NA
        return((x-mean(x,na.rm = T))/sd(x,na.rm = T))
      }))
      return(a3)
    }
    get_divsion <- function(m,n){
      m <- m*100+1
      n <- n*100+1
      m <- (m-n)/(m+n)
      return(m)
    }
    
    a <- get_divsion(a1,a2) ## TSS centered
    b <- get_divsion(b1,b2) ## TTS centered
    
    a1 <- scale_row(a1)
    a2 <- scale_row(a2)
    b1 <- scale_row(b1)
    b2 <- scale_row(b2)
    
    zlim1=range(a1,a2,b1,b2,na.rm = T)
    zlim2=range(a,b,na.rm = T)
    
    
    x11()
    layout.show(layout(matrix(c(1:12),ncol=3,byrow = T),heights = c(5,1,5,1,5,1)))
    
    # 1,1 and 1,2
    par(mar=c(3,4,4,1))
    h <- get.unscaledgene.heatplot.list(r=a1,mid.point = 0,genes = genes,zlim = zlim1,cols = c("white","gray90","black"),x = x,gene.len.cut =gene.len.cut )
    image(h$x,1:nrow(h$z$z),t(h$z$z),zlim=h$zlim,col=h$z$col.z,axes=FALSE,xlab=" ",ylab=" ",main="(I) Depleted",cex.main=4);box();
    axis(1,at=at,labels = label,cex.axis=3,padj = 0.5)
    lines(h$border.lines[[1]],1:nrow(h$z$z),lty=3,lwd=1,col="green4",xlab="",ylab="")
    lines(h$border.lines[[2]],rev(1:nrow(h$z$z)),lty=3,lwd=1,col="green4",xlab="",ylab="")
    h <- get.unscaledgene.heatplot.list(r=a2,mid.point = 0,genes = genes,zlim = zlim1,cols = c("white","gray90","black"),x = x,gene.len.cut =gene.len.cut)
    image(h$x,1:nrow(h$z$z),t(h$z$z),zlim=h$zlim,col=h$z$col.z,axes=FALSE,xlab=" ",ylab=" ",main="(II) Non-depleted",cex.main=4);box();
    axis(1,at=at,labels = label,cex.axis=3,padj = 0.5)
    lines(h$border.lines[[1]],1:nrow(h$z$z),lty=3,lwd=1,col="green4",xlab="",ylab="")
    lines(h$border.lines[[2]],rev(1:nrow(h$z$z)),lty=3,lwd=1,col="green4",xlab="",ylab="")
    h <- get.unscaledgene.heatplot.list(r=a,mid.point = 0,genes = genes,zlim = zlim2,cols = c("navyblue","white","red4"),x = x,gene.len.cut =gene.len.cut)
    image(h$x,1:nrow(h$z$z),t(h$z$z),zlim=h$zlim,col=h$z$col.z,axes=FALSE,xlab=" ",ylab=" ",main="Difference, (I-II)/(I+II)",cex.main=4);box();
    axis(1,at=at,labels = label,cex.axis=3,padj = 0.5)
    lines(h$border.lines[[1]],1:nrow(h$z$z),lty=3,lwd=1,col="green4",xlab="",ylab="")
    lines(h$border.lines[[2]],rev(1:nrow(h$z$z)),lty=3,lwd=1,col="green4",xlab="",ylab="")
    
    par(mar=c(3,3,1,3))
    h <- get.unscaledgene.heatplot.list(r=a1,mid.point = 0,genes = genes,zlim = zlim1,cols = c("white","gray90","black"),x = x,gene.len.cut =gene.len.cut)
    image(h$z$scale$x,h$z$scale$y,h$z$scale$z,xlab="",ylab="",col=h$z$col.scale,axes=F)
    axis(1,at=c(min(h$z$scale$x),1,mean(h$z$scale$x),max(h$z$scale$x)),round(c(min(h$z$scale$x),1,mean(h$z$scale$x),max(h$z$scale$x))),cex.axis=3)
    image(h$z$scale$x,h$z$scale$y,h$z$scale$z,xlab="",ylab="",col=h$z$col.scale,axes=F)
    axis(1,at=c(min(h$z$scale$x),1,mean(h$z$scale$x),max(h$z$scale$x)),round(c(min(h$z$scale$x),1,mean(h$z$scale$x),max(h$z$scale$x))),cex.axis=3)
    h <- get.unscaledgene.heatplot.list(r=a,mid.point = 0,genes = genes,zlim = zlim2,cols = c("navyblue","white","red4"))
    image(h$z$scale$x,h$z$scale$y,h$z$scale$z,xlab="",ylab="",col=h$z$col.scale,axes=F)
    axis(1,at=c(min(h$z$scale$x),1,mean(h$z$scale$x),max(h$z$scale$x)),round(c(min(h$z$scale$x),1,mean(h$z$scale$x),max(h$z$scale$x))),cex.axis=3)
    
    
    par(mar=c(3,4,4,1))
    h <- get.unscaledgene.heatplot.list(r=b1,mid.point = 0,genes = genes,zlim = zlim1,cols = c("white","gray90","black"),align = "CPS",x = x,gene.len.cut =gene.len.cut)
    image(h$x,1:nrow(h$z$z),t(h$z$z),zlim=h$zlim,col=h$z$col.z,axes=FALSE,xlab=" ",ylab=" ",main="(I) Depleted",cex.main=4);box();
    axis(1,at=at,labels = label,cex.axis=3,padj = 0.5)
    lines(h$border.lines[[1]],1:nrow(h$z$z),lty=3,lwd=1,col="green4",xlab="",ylab="")
    lines(h$border.lines[[2]],rev(1:nrow(h$z$z)),lty=3,lwd=1,col="green4",xlab="",ylab="")
    h <- get.unscaledgene.heatplot.list(r=b2,mid.point = 0,genes = genes,zlim = zlim1,cols = c("white","gray90","black"),align = "CPS",x = x,gene.len.cut =gene.len.cut)
    image(h$x,1:nrow(h$z$z),t(h$z$z),zlim=h$zlim,col=h$z$col.z,axes=FALSE,xlab=" ",ylab=" ",main="(II) Non-depleted",cex.main=4);box();
    axis(1,at=at,labels = label,cex.axis=3,padj = 0.5)
    lines(h$border.lines[[1]],1:nrow(h$z$z),lty=3,lwd=1,col="green4",xlab="",ylab="")
    lines(h$border.lines[[2]],rev(1:nrow(h$z$z)),lty=3,lwd=1,col="green4",xlab="",ylab="")
    h <- get.unscaledgene.heatplot.list(r=b,mid.point = 0,genes = genes,zlim = zlim2,cols = c("navyblue","white","red4"),align = "CPS",x = x,gene.len.cut =gene.len.cut)
    image(h$x,1:nrow(h$z$z),t(h$z$z),zlim=h$zlim,col=h$z$col.z,axes=FALSE,xlab=" ",ylab=" ",main="Difference, (I-II)/(I+II)",cex.main=4);box();
    axis(1,at=at,labels = label,cex.axis=3,padj = 0.5)
    lines(h$border.lines[[1]],1:nrow(h$z$z),lty=3,lwd=1,col="green4",xlab="",ylab="")
    lines(h$border.lines[[2]],rev(1:nrow(h$z$z)),lty=3,lwd=1,col="green4",xlab="",ylab="")
    
    par(mar=c(3,3,1,3))
    h <- get.unscaledgene.heatplot.list(r=b1,mid.point = 0,genes = genes,zlim = zlim1,cols = c("white","gray90","black"),x = x,gene.len.cut =gene.len.cut)
    image(h$z$scale$x,h$z$scale$y,h$z$scale$z,xlab="",ylab="",col=h$z$col.scale,axes=F)
    axis(1,at=c(min(h$z$scale$x),1,mean(h$z$scale$x),max(h$z$scale$x)),round(c(min(h$z$scale$x),1,mean(h$z$scale$x),max(h$z$scale$x))),cex.axis=3)
    image(h$z$scale$x,h$z$scale$y,h$z$scale$z,xlab="",ylab="",col=h$z$col.scale,axes=F)
    axis(1,at=c(min(h$z$scale$x),1,mean(h$z$scale$x),max(h$z$scale$x)),round(c(min(h$z$scale$x),1,mean(h$z$scale$x),max(h$z$scale$x))),cex.axis=3)
    h <- get.unscaledgene.heatplot.list(r=b,mid.point = 0,genes = genes,zlim = zlim2,cols = c("navyblue","white","red4"))
    image(h$z$scale$x,h$z$scale$y,h$z$scale$z,xlab="",ylab="",col=h$z$col.scale,axes=F)
    axis(1,at=c(min(h$z$scale$x),1,mean(h$z$scale$x),max(h$z$scale$x)),round(c(min(h$z$scale$x),1,mean(h$z$scale$x),max(h$z$scale$x))),cex.axis=2.2)
    
    dev.copy(device=jpeg,file=paste0(outfile),width=1500,height=1800);dev.off();
    
  }
  
  
  load("Revised_MultigeneHeatmap_H3K36me3.RData")
  l1 <- l
  load("Revised_MultigeneHeatmap_H3 (round1).RData")
  l2 <- l;rm(l)
  normF1 <- res[res$Factor=="H3K36me3",]$normF
  normF2 <- res[res$Factor=="H3 (round1)",]$normF
  myheatplotfun(l1,l2,genes,genes,outfile="H3K36me3_AllGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,rp,outfile="H3K36me3_RPGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,de.genes,outfile="H3K36me3_DEGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,introns,outfile="H3K36me3_IntronGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,opn,outfile="H3K36me3_OPENGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,dpn,outfile="H3K36me3_CLOSEDGenes_revised_Rel.jpeg",normF1,normF2)
  rm(l1,l,l2,normF1,normF2)
  
  load("Revised_MultigeneHeatmap_H3K4me3.RData")
  l1 <- l
  load("Revised_MultigeneHeatmap_H3.RData")
  l2 <- l
  normF1 <- res[res$Factor=="H3K4me3",]$normF
  normF2 <- res[res$Factor=="H3",]$normF
  myheatplotfun(l1,l2,genes,genes,outfile="H3K4me3_AllGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,rp,outfile="H3K4me3_RPGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,de.genes,outfile="H3K4me3_DEGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,introns,outfile="H3K4me3_IntronGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,opn,outfile="H3K4me3_OPENGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,dpn,outfile="H3K4me3_CLOSEDGenes_revised_Rel.jpeg",normF1,normF2)
  rm(l1,l,l2,normF1,normF2)
  
  load("Revised_MultigeneHeatmap_H3K36me2.RData")
  l1 <- l
  load("Revised_MultigeneHeatmap_H3.RData")
  l2 <- l
  normF1 <- res[res$Factor=="H3K36me2",]$normF
  normF2 <- res[res$Factor=="H3",]$normF
  myheatplotfun(l1,l2,genes,genes,outfile="H3K36me2_AllGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,rp,outfile="H3K36me2_RPGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,de.genes,outfile="H3K36me2_DEGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,introns,outfile="H3K36me2_IntronGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,opn,outfile="H3K36me2_OPENGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,dpn,outfile="H3K36me2_CLOSEDGenes_revised_Rel.jpeg",normF1,normF2)
  rm(l1,l,l2,normF1,normF2)
  
  load("Revised_MultigeneHeatmap_S2P.RData")
  l1 <- l
  load("Revised_MultigeneHeatmap_RNAPII.RData")
  l2 <- l
  normF1 <- res[res$Factor=="S2P",]$normF
  normF2 <- res[res$Factor=="RNAPII",]$normF
  myheatplotfun(l1,l2,genes,genes,outfile="S2P_AllGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,rp,outfile="S2P_RPGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,de.genes,outfile="S2P_DEGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,introns,outfile="S2P_IntronGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,opn,outfile="S2P_OPENGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,dpn,outfile="S2P_CLOSEDGenes_revised_Rel.jpeg",normF1,normF2)
  rm(l1,l,l2,normF1,normF2)
  
  load("Revised_MultigeneHeatmap_S5P.RData")
  l1 <- l
  load("Revised_MultigeneHeatmap_RNAPII.RData")
  l2 <- l
  normF1 <- res[res$Factor=="S5P",]$normF
  normF2 <- res[res$Factor=="RNAPII",]$normF
  myheatplotfun(l1,l2,genes,genes,outfile="S5P_AllGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,rp,outfile="S5P_RPGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,de.genes,outfile="S5P_DEGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,introns,outfile="S5P_IntronGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,opn,outfile="S5P_OPENGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,dpn,outfile="S5P_CLOSEDGenes_revised_Rel.jpeg",normF1,normF2)
  rm(l1,l,l2,normF1,normF2)
  
  load("Revised_MultigeneHeatmap_Spt6.RData")
  l1 <- l
  load("Revised_MultigeneHeatmap_RNAPII (round1).RData")
  l2 <- l
  normF1 <- res[res$Factor=="Spt6",]$normF
  normF2 <- res[res$Factor=="RNAPII (round1)",]$normF
  myheatplotfun(l1,l2,genes,genes,outfile="Spt6_AllGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,rp,outfile="Spt6_RPGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,de.genes,outfile="Spt6_DEGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,introns,outfile="Spt6_IntronGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,opn,outfile="Spt6_OPENGenes_revised_Rel.jpeg",normF1,normF2)
  myheatplotfun(l1,l2,genes,dpn,outfile="Spt6_CLOSEDGenes_revised_Rel.jpeg",normF1,normF2)
  rm(l1,l,l2,normF1,normF2)
  
  
  #### Relative metagene
  rm(list=ls())
  load(file="Gene groups.RData")
  vector.resizing <- function(x,final.len){
    y <- vector()
    len <- length(x)
    y <-spline(1:len,x,n=final.len)$y
    return(y)
  }
  
  res <- read.delim("Normalization.factors.merged.txt",header = T)
  myplot <- function(pl1){
    require(ggplot2)
    q <- ggplot(pl1, aes(x=pos,y=signal,ymax=signal+se, ymin=signal-se,col=condition))
    q <- q + geom_vline(xintercept = c(50,150),colour="gray50", linetype = "dashed")
    q <- q + geom_linerange(col='gray') + geom_line() #+ ylim(-0.03,0.15) #ylim(-0.01,0.045) #ylim(-0.04,0.1)
    q <- q + scale_color_manual(values = c("Depleted-1" = "pink","Depleted-2" = "red2",
                                           "Non-depleted-1" = "lightblue","Non-depleted-2" = "slateblue2",
                                           "Depleted" = "red","Non-depleted" = "slateblue"))
    q <- q + facet_wrap(~profile, ncol=3,scales = "free")
    q <- q + scale_x_continuous(breaks = c(0,50,150,199), labels=c("","TSS","CPS","")) +ylab("average signal") + xlab("") 
    q <- q + theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = NA),axis.line = element_line(size = 0.7, colour = "gray50"))
    q <- q + theme(axis.text.x = element_text(colour = 'black',size=14),
                   axis.text.y = element_text(colour = 'black',size=14),
                   axis.title = element_text(colour = 'black',size=15,face="bold"),
                   strip.text.x = element_text(colour = 'black',size=14,face="bold"))
    q <- q + theme(legend.text =element_text(colour = 'black',size=14),
                   legend.title = element_text(colour = 'black',size=14,face="bold"))
    q <- q + guides(color=guide_legend(title="Spn1 levels"))
    return(q)
    
  }
  mymetagenefun <- function(prot,cnt,genes,g,normF=c(1,1),normF2=c(1,1)){
    f1 = paste("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_",prot,".RData",sep="")
    f2 = paste("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_",cnt,".RData",sep="")
    
    load(f1)
    l1 <- l
    load(f2)
    l2 <- l
    rm(l)
    
    
    depl1 <- matrix(NA,ncol=199,nrow = 5134)
    nondepl1 <- matrix(NA,ncol=199,nrow = 5134)
    for(i in 1:length(genes$tracking_id)){
      j = round(genes[i,]$width/10)
      v <- l1[[1]][i,51:(49+j)] + l1[[2]][i,51:(49+j)]
      v <- vector.resizing(v,100)
      v <- c(l1[[1]][i,1:50]+l1[[2]][i,1:50],
             v,
             l1[[1]][i,1525:1573]+l1[[2]][i,1525:1573]
      )
      depl1[i,] <- v
      
      v <- l1[[3]][i,51:(49+j)] + l1[[4]][i,51:(49+j)]
      v <- vector.resizing(v,100)
      v <- c(l1[[3]][i,1:50]+l1[[4]][i,1:50],
             v,
             l1[[3]][i,1525:1573]+l1[[4]][i,1525:1573]
      )
      nondepl1[i,] <- v
    }
    rownames(depl1) <- rownames(nondepl1) <- rownames(l1[[1]])
    depl1 <- subset(depl1,rownames(depl1)%in%g$tracking_id)
    nondepl1 <- subset(nondepl1,rownames(nondepl1)%in%g$tracking_id)
    
    depl2 <- matrix(NA,ncol=199,nrow = 5134)
    nondepl2 <- matrix(NA,ncol=199,nrow = 5134)
    for(i in 1:length(genes$tracking_id)){
      j = round(genes[i,]$width/10)
      v <- l2[[1]][i,51:(49+j)] + l2[[2]][i,51:(49+j)]
      v <- vector.resizing(v,100)
      v <- c(l2[[1]][i,1:50]+l2[[2]][i,1:50],
             v,
             l2[[1]][i,1525:1573]+l2[[2]][i,1525:1573]
      )
      depl2[i,] <- v
      
      v <- l2[[3]][i,51:(49+j)] + l2[[4]][i,51:(49+j)]
      v <- vector.resizing(v,100)
      v <- c(l2[[3]][i,1:50]+l2[[4]][i,1:50],
             v,
             l2[[3]][i,1525:1573]+l2[[4]][i,1525:1573]
      )
      nondepl2[i,] <- v
    }
    rownames(depl2) <- rownames(nondepl2) <- rownames(l2[[1]])
    depl2 <- subset(depl2,rownames(depl2)%in%g$tracking_id)
    nondepl2 <- subset(nondepl2,rownames(nondepl2)%in%g$tracking_id)
    
    depl1 <- depl1*normF[1]
    nondepl1 <- nondepl1*normF[2]
    depl2 <- depl2*normF2[1]
    nondepl2 <- nondepl2*normF2[2]
    
    #depl1[depl1==0] <- NA
    #depl2[depl2==0] <- NA
    #nondepl1[nondepl1==0] <- NA
    #nondepl2[nondepl2==0] <- NA
    
    depl1 <- (depl1*100+1)/(depl2*100+1)
    nondepl1 <- (nondepl1*100+1)/(nondepl2*100+1)
    depl1[!is.finite(depl1)] <-  NA
    nondepl1[!is.finite(nondepl1)] <- NA
    
    get.avg <- function(m){
      df_m <- cbind(1:ncol(m),
                    as.data.frame(apply(m,2,function(x){
                      x <- x[x>quantile(x,0.05,na.rm=T) & x<quantile(x,0.99,na.rm=T)]
                      return(mean(x,na.rm=T))
                    })),
                    as.data.frame(apply(m,2,function(x){
                      x <- x[x>quantile(x,0.05,na.rm=T) & x<quantile(x,0.99,na.rm=T)]
                      return(sd(x,na.rm=T))
                    })), 
                    as.data.frame(apply(m,2, function(x){
                      x <- x[x>quantile(x,0.05,na.rm=T) & x<quantile(x,0.99,na.rm=T)]
                      return(sd(x,na.rm=TRUE)/sqrt(length(x)))
                    } ))
      )
      names(df_m) <- c("pos","signal", "sd","se")
      #df_m$signal <- df_m$signal*nF
      #df_m$sd <- df_m$sd*nF
      return(df_m)
    }
    
    pl <- rbind(cbind(get.avg(depl1),condition="Depleted"),
                cbind(get.avg(nondepl1),condition="Non-depleted")
    )
    return(pl)
  }
  pl <- c()
  
  profiles <- c("H3K36me2","H3","H3K36me3","H3 (round1)",
                "H3K4me3","H3","S2P","RNAPII","S5P","RNAPII","Spt6","RNAPII (round1)")
  #profiles <- c("Spt6","RNAPII (round1)","Spt6","RNAPII")
  pl <- c()
  for(i in seq(1,11,2)){
    prot = profiles[i]
    cnt = profiles[i+1]
    #pl <- rbind(pl,cbind(profile=prot,genetype = "all genes", mymetagenefun(prot,cnt,genes,g=genes,normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    #pl <- rbind(pl,cbind(profile=prot,genetype = "RP genes", mymetagenefun(prot,cnt,genes,g=rp,normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    #pl <- rbind(pl,cbind(profile=prot,genetype = "Intron genes", mymetagenefun(prot,cnt,genes,g=introns,normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    #pl <- rbind(pl,cbind(profile=prot,genetype = "DE genes", mymetagenefun(prot,cnt,genes,g=de.genes,normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    #pl <- rbind(pl,cbind(profile=prot,genetype = "OPEN genes", mymetagenefun(prot,cnt,genes,g=opn,normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    #pl <- rbind(pl,cbind(profile=prot,genetype = "CLOSED genes", mymetagenefun(prot,cnt,genes,g=dpn,normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    pl <- rbind(pl,cbind(profile=prot,genetype = "<0.5kb", mymetagenefun(prot,cnt,genes,g=genes[genes$lengp==1,],normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    pl <- rbind(pl,cbind(profile=prot,genetype = "0.5kb-1kb", mymetagenefun(prot,cnt,genes,g=genes[genes$lengp==2,],normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    pl <- rbind(pl,cbind(profile=prot,genetype = "1kb-2kb", mymetagenefun(prot,cnt,genes,g=genes[genes$lengp==3,],normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    pl <- rbind(pl,cbind(profile=prot,genetype = ">2kb", mymetagenefun(prot,cnt,genes,g=genes[genes$lengp==4,],normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
  }
  
  #load(file="Metagenplots_ChIPSpike_final_relative.RData")
  save(pl,file="Metagenplots_ChIPSpike_final_relative.RData")
  
  pdf("MetagenePlots_SpikeNorm_AllFactors_Merged_rel.pdf",width = 14,height = 8)
  myplot(pl[pl$genetype=="all genes",]) + ggtitle("All verified protein coding genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype=="RP genes",]) + ggtitle("RP[PSL] genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype=="Intron genes",]) + ggtitle("Intron containing genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype=="DE genes",]) + ggtitle("Downregulated genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype=="OPEN genes",]) + ggtitle("OPEN genes (Zuagg et al)")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype=="CLOSED genes",]) + ggtitle("CLOSED genes (Zuagg et al)")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  
  myplot(pl[pl$genetype=="<0.5kb",]) + ggtitle("Genes with length <0.5kb")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype=="0.5kb-1kb",]) + ggtitle("Genes with length 0.5kb-1kb")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype=="1kb-2kb",]) + ggtitle("Genes with length 1kb-2kb")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype==">2kb",]) + ggtitle("Genes with length >2kb")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  dev.off()
  
}

## Genes grouped in class
if(F){
  rm(list=ls())
  load(file="Gene groups.RData")
  vector.resizing <- function(x,final.len){
    y <- vector()
    len <- length(x)
    y <-spline(1:len,x,n=final.len)$y
    return(y)
  }
  mymetagenefun <- function(l,genes,g,normF=c(1,1)){
    depl1 <- matrix(NA,ncol=199,nrow = 5134)
    nondepl1 <- matrix(NA,ncol=199,nrow = 5134)
    for(i in 1:length(genes$tracking_id)){
      j = round(genes[i,]$width/10)
      v <- l[[1]][i,51:(49+j)] + l[[2]][i,51:(49+j)]
      v <- vector.resizing(v,100)
      v <- c(l[[1]][i,1:50]+l[[2]][i,1:50],
             v,
             l[[1]][i,1525:1573]+l[[2]][i,1525:1573]
      )
      depl1[i,] <- v
      
      v <- l[[3]][i,51:(49+j)] + l[[4]][i,51:(49+j)]
      v <- vector.resizing(v,100)
      v <- c(l[[3]][i,1:50]+l[[4]][i,1:50],
             v,
             l[[3]][i,1525:1573]+l[[4]][i,1525:1573]
      )
      nondepl1[i,] <- v
    }
    rownames(depl1) <- rownames(nondepl1) <- rownames(l[[1]])
    
    depl1 <- subset(depl1,rownames(depl1)%in%g$tracking_id)
    nondepl1 <- subset(nondepl1,rownames(nondepl1)%in%g$tracking_id)
    
    depl1 <- depl1*normF[1]
    nondepl1 <- nondepl1*normF[2]
    
    rat <- (depl1*100+1)/(nondepl1*100+1)
    
    
    get.avg <- function(m){
      df_m <- cbind(1:ncol(m),
                    as.data.frame(apply(m,2,function(x){
                      x <- x[x>quantile(x,0.02,na.rm=T) & x<quantile(x,0.99,na.rm=T)]
                      return(mean(x,na.rm=T))
                    })),
                    as.data.frame(apply(m,2,function(x){
                      x <- x[x>quantile(x,0.02,na.rm=T) & x<quantile(x,0.99,na.rm=T)]
                      return(sd(x,na.rm=T))
                    })), 
                    as.data.frame(apply(m,2, function(x){
                      x <- x[x>quantile(x,0.05,na.rm=T) & x<quantile(x,0.99,na.rm=T)]
                      return(sd(x,na.rm=TRUE)/sqrt(length(x)))
                    } ))
      )
      names(df_m) <- c("pos","signal", "sd","se")
      #df_m$signal <- df_m$signal*nF
      #df_m$sd <- df_m$sd*nF
      return(df_m)
    }
    pl <- cbind(get.avg(rat),condition="Depl./Non-depl.")
    return(pl)
  }
  myplot <- function(pl1){
    require(ggplot2)
    q <- ggplot(pl1, aes(x=pos,y=signal,ymax=signal+se, ymin=signal-se,col=genetype))
    q <- q + geom_hline(yintercept = c(1),col="purple")
    q <- q + geom_vline(xintercept = c(50,150),colour="gray50", linetype = "dashed")
    q <- q + geom_linerange(col='gray') + geom_line() #+ ylim(-0.03,0.15) #ylim(-0.01,0.045) #ylim(-0.04,0.1)
    q <- q + scale_color_manual(values = c("Depleted-1" = "pink","Depleted-2" = "red2",
                                           "Non-depleted-1" = "lightblue","Non-depleted-2" = "slateblue2",
                                           "Depleted" = "red","Non-depleted" = "slateblue",
                                           "OPEN genes"="green","CLOSED genes"="red",
                                           "High Pol2"="red4","Medium Pol2"="red","Low Pol2"="pink",
                                           "<0.5kb"="yellow","0.5kb-1kb"="orange","1kb-2kb"="red",">2kb"="red4",
                                           "rest genes"="pink","DE genes"= "red", "RP genes"="red"))
    q <- q + facet_wrap(~profile, ncol=4,scales = "free")
    q <- q + scale_x_continuous(breaks = c(0,50,150,199), labels=c("","TSS","CPS","")) +ylab("average signal (depleted/non-depleted)") + xlab("") 
    q <- q + theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = NA),axis.line = element_line(size = 0.7, colour = "gray50"))
    q <- q + theme(axis.text.x = element_text(colour = 'black',size=14),
                   axis.text.y = element_text(colour = 'black',size=14),
                   axis.title = element_text(colour = 'black',size=15,face="bold"),
                   strip.text.x = element_text(colour = 'black',size=14,face="bold"))
    q <- q + theme(legend.text =element_text(colour = 'black',size=14),
                   legend.title = element_text(colour = 'black',size=14,face="bold"))
    q <- q + guides(color=guide_legend(title="Gene class"))
    return(q)
    
  }
  
  pl <- c()
  res <- read.delim("Normalization.factors.merged.txt",header = T)
  for(prot in protein){
    f = paste("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_",prot,".RData",sep="")
    if(file.exists(f)){
      cat(prot,"\n")
      cat(res[res$Factor==prot,]$normF,"\n\n")
      load(f)
      pl <- rbind(pl,
                  cbind(profile=prot,genetype = "all genes", mymetagenefun(l,genes,g=genes,normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "RP genes", mymetagenefun(l,genes,g=rp,normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "Intron genes", mymetagenefun(l,genes,g=introns, normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "DE genes", mymetagenefun(l,genes,g=de.genes, normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "OPEN genes", mymetagenefun(l,genes,g=opn, normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "CLOSED genes", mymetagenefun(l,genes,g=dpn, normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "High Pol2", mymetagenefun(l,genes,g=genes[genes$Pol2OccGp==3,],normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "Medium Pol2", mymetagenefun(l,genes,g=genes[genes$Pol2OccGp==2,],normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "Low Pol2", mymetagenefun(l,genes,g=genes[genes$Pol2OccGp==1,], normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "<0.5kb", mymetagenefun(l,genes,g=genes[genes$lengp==1,], normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "0.5kb-1kb", mymetagenefun(l,genes,g=genes[genes$lengp==2,], normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "1kb-2kb", mymetagenefun(l,genes,g=genes[genes$lengp==3,], normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = ">2kb", mymetagenefun(l,genes,g=genes[genes$lengp==4,], normF = res[res$Factor==prot,]$normF)),
                  cbind(profile=prot,genetype = "rest genes", mymetagenefun(l,genes,g=rest,normF = res[res$Factor==prot,]$normF))
      )
    }
  }
  save(pl,file="Metagene_ChIPSpike_merged_deplbynondepl.RData")
  
  pdf("MetagenePlots_SpikeNorm_AllFactors_Merged_geneByGroups.pdf",width = 16,height = 8)
  myplot(pl[pl$genetype%in%c("OPEN genes","CLOSED genes"),]) + ggtitle("OPEN and CLOSED genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype%in%c("High Pol2","Medium Pol2","Low Pol2"),]) + ggtitle("Genes grouped by WT RNAPII occupancy")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype%in%c("<0.5kb","0.5kb-1kb","1kb-2kb",">2kb"),]) + ggtitle("Genes grouped by length")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype%in%c("rest genes","DE genes"),]) + ggtitle("Genes grouped expression status")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  dev.off()
  
  
  
}

### depl/non.depl for H3and RNAPII normalzied profiles
if(F){
  rm(list=ls())
  load(file="Gene groups.RData")
  vector.resizing <- function(x,final.len){
    y <- vector()
    len <- length(x)
    y <-spline(1:len,x,n=final.len)$y
    return(y)
  }
  
  res <- read.delim("Normalization.factors.merged.txt",header = T)
  myplot <- function(pl1){
    require(ggplot2)
    q <- ggplot(pl1, aes(x=pos,y=signal,ymax=signal+se, ymin=signal-se,col=genetype))
    q <- q + geom_hline(yintercept = c(1),col="purple")
    q <- q + geom_vline(xintercept = c(50,150),colour="gray50", linetype = "dashed")
    q <- q + geom_linerange(col='gray') + geom_line() #+ ylim(-0.03,0.15) #ylim(-0.01,0.045) #ylim(-0.04,0.1)
    q <- q + scale_color_manual(values = c("Depleted-1" = "pink","Depleted-2" = "red2",
                                           "Non-depleted-1" = "lightblue","Non-depleted-2" = "slateblue2",
                                           "Depleted" = "red","Non-depleted" = "slateblue",
                                           "OPEN genes"="green","CLOSED genes"="red",
                                           "High Pol2"="red4","Medium Pol2"="red","Low Pol2"="pink",
                                           "<0.5kb"="yellow","0.5kb-1kb"="orange","1kb-2kb"="red",">2kb"="red4",
                                           "non-DE genes"="pink","non-RP genes"="pink","DE genes"= "red", "RP genes"="red",
                                           "Non-Intron genes"="pink","Intron genes"="red",
                                           "RP-Intron genes"="red","nonRP-Intron genes"="pink"
                                           ))
    q <- q + facet_wrap(~profile, ncol=3,scales = "free")
    q <- q + scale_x_continuous(breaks = c(0,50,150,199), labels=c("","TSS","CPS","")) +ylab("average signal (depleted/non-depleted)") + xlab("") 
    q <- q + theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = NA),axis.line = element_line(size = 0.7, colour = "gray50"))
    q <- q + theme(axis.text.x = element_text(colour = 'black',size=14),
                   axis.text.y = element_text(colour = 'black',size=14),
                   axis.title = element_text(colour = 'black',size=15,face="bold"),
                   strip.text.x = element_text(colour = 'black',size=14,face="bold"))
    q <- q + theme(legend.text =element_text(colour = 'black',size=14),
                   legend.title = element_text(colour = 'black',size=14,face="bold"))
    q <- q + guides(color=guide_legend(title="Gene class"))
    return(q)
    
  }
  mymetagenefun <- function(prot,cnt,genes,g,normF=c(1,1),normF2=c(1,1)){
    f1 = paste("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_",prot,".RData",sep="")
    f2 = paste("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_",cnt,".RData",sep="")
    
    load(f1)
    l1 <- l
    load(f2)
    l2 <- l
    rm(l)
    
    
    depl1 <- matrix(NA,ncol=199,nrow = 5134)
    nondepl1 <- matrix(NA,ncol=199,nrow = 5134)
    for(i in 1:length(genes$tracking_id)){
      j = round(genes[i,]$width/10)
      v <- l1[[1]][i,51:(49+j)] + l1[[2]][i,51:(49+j)]
      v <- vector.resizing(v,100)
      v <- c(l1[[1]][i,1:50]+l1[[2]][i,1:50],
             v,
             l1[[1]][i,1525:1573]+l1[[2]][i,1525:1573]
      )
      depl1[i,] <- v
      
      v <- l1[[3]][i,51:(49+j)] + l1[[4]][i,51:(49+j)]
      v <- vector.resizing(v,100)
      v <- c(l1[[3]][i,1:50]+l1[[4]][i,1:50],
             v,
             l1[[3]][i,1525:1573]+l1[[4]][i,1525:1573]
      )
      nondepl1[i,] <- v
    }
    rownames(depl1) <- rownames(nondepl1) <- rownames(l1[[1]])
    depl1 <- subset(depl1,rownames(depl1)%in%g$tracking_id)
    nondepl1 <- subset(nondepl1,rownames(nondepl1)%in%g$tracking_id)
    
    depl2 <- matrix(NA,ncol=199,nrow = 5134)
    nondepl2 <- matrix(NA,ncol=199,nrow = 5134)
    for(i in 1:length(genes$tracking_id)){
      j = round(genes[i,]$width/10)
      v <- l2[[1]][i,51:(49+j)] + l2[[2]][i,51:(49+j)]
      v <- vector.resizing(v,100)
      v <- c(l2[[1]][i,1:50]+l2[[2]][i,1:50],
             v,
             l2[[1]][i,1525:1573]+l2[[2]][i,1525:1573]
      )
      depl2[i,] <- v
      
      v <- l2[[3]][i,51:(49+j)] + l2[[4]][i,51:(49+j)]
      v <- vector.resizing(v,100)
      v <- c(l2[[3]][i,1:50]+l2[[4]][i,1:50],
             v,
             l2[[3]][i,1525:1573]+l2[[4]][i,1525:1573]
      )
      nondepl2[i,] <- v
    }
    rownames(depl2) <- rownames(nondepl2) <- rownames(l2[[1]])
    depl2 <- subset(depl2,rownames(depl2)%in%g$tracking_id)
    nondepl2 <- subset(nondepl2,rownames(nondepl2)%in%g$tracking_id)
    
    depl1 <- depl1*normF[1]
    nondepl1 <- nondepl1*normF[2]
    depl2 <- depl2*normF2[1]
    nondepl2 <- nondepl2*normF2[2]
    
    #depl1[depl1==0] <- NA
    #depl2[depl2==0] <- NA
    #nondepl1[nondepl1==0] <- NA
    #nondepl2[nondepl2==0] <- NA
    
    depl1 <- (depl1*100+1)/(depl2*100+1)
    nondepl1 <- (nondepl1*100+1)/(nondepl2*100+1)
    depl1[!is.finite(depl1)] <-  NA
    nondepl1[!is.finite(nondepl1)] <- NA
    
    rat <- (depl1*100+1)/(nondepl1*100+1)
    rat[!is.finite(rat)] <-  NA
    
    get.avg <- function(m){
      df_m <- cbind(1:ncol(m),
                    as.data.frame(apply(m,2,function(x){
                      x <- x[x>quantile(x,0.05,na.rm=T) & x<quantile(x,0.99,na.rm=T)]
                      return(mean(x,na.rm=T))
                    })),
                    as.data.frame(apply(m,2,function(x){
                      x <- x[x>quantile(x,0.05,na.rm=T) & x<quantile(x,0.99,na.rm=T)]
                      return(sd(x,na.rm=T))
                    })), 
                    as.data.frame(apply(m,2, function(x){
                      x <- x[x>quantile(x,0.05,na.rm=T) & x<quantile(x,0.99,na.rm=T)]
                      return(sd(x,na.rm=TRUE)/sqrt(length(x)))
                    } ))
      )
      names(df_m) <- c("pos","signal", "sd","se")
      #df_m$signal <- df_m$signal*nF
      #df_m$sd <- df_m$sd*nF
      return(df_m)
    }
    
    pl <- cbind(get.avg(rat),condition="Depl./Non-depl.")
    return(pl)
  }
  pl <- c()
  
  profiles <- c("H3K36me2","H3","H3K36me3","H3 (round1)",
                "H3K4me3","H3","S2P","RNAPII","S5P","RNAPII","Spt6","RNAPII (round1)")
  load(file="Metagenplots_ChIPSpike_final_relative_deplbynondepl.RData")
  for(i in seq(1,11,2)){
    prot = profiles[i]
    cnt = profiles[i+1]
    #pl <- rbind(pl,cbind(profile=prot,genetype = "non-DE genes", mymetagenefun(prot,cnt,genes,g=genes[!genes$tracking_id%in%de.genes$tracking_id,],normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    #pl <- rbind(pl,cbind(profile=prot,genetype = "non-RP genes", mymetagenefun(prot,cnt,genes,g=genes[!genes$tracking_id%in%rp$tracking_id,],normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    
    #pl <- rbind(pl,cbind(profile=prot,genetype = "all genes", mymetagenefun(prot,cnt,genes,g=genes,normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    #pl <- rbind(pl,cbind(profile=prot,genetype = "RP genes", mymetagenefun(prot,cnt,genes,g=rp,normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    #pl <- rbind(pl,cbind(profile=prot,genetype = "Intron genes", mymetagenefun(prot,cnt,genes,g=introns,normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    #pl <- rbind(pl,cbind(profile=prot,genetype = "DE genes", mymetagenefun(prot,cnt,genes,g=de.genes,normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    #pl <- rbind(pl,cbind(profile=prot,genetype = "OPEN genes", mymetagenefun(prot,cnt,genes,g=opn,normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    #pl <- rbind(pl,cbind(profile=prot,genetype = "CLOSED genes", mymetagenefun(prot,cnt,genes,g=dpn,normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    
    #pl <- rbind(pl,cbind(profile=prot,genetype = "High Pol2", mymetagenefun(prot,cnt,genes,g=genes[genes$Pol2OccGp==3,],normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    #pl <- rbind(pl,cbind(profile=prot,genetype = "Medium Pol2", mymetagenefun(prot,cnt,genes,g=genes[genes$Pol2OccGp==2,],normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    #pl <- rbind(pl,cbind(profile=prot,genetype = "Low Pol2", mymetagenefun(prot,cnt,genes,g=genes[genes$Pol2OccGp==1,],normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    
    #pl <- rbind(pl,cbind(profile=prot,genetype = "<0.5kb", mymetagenefun(prot,cnt,genes,g=genes[genes$lengp==1,],normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    #pl <- rbind(pl,cbind(profile=prot,genetype = "0.5kb-1kb", mymetagenefun(prot,cnt,genes,g=genes[genes$lengp==2,],normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    #pl <- rbind(pl,cbind(profile=prot,genetype = "1kb-2kb", mymetagenefun(prot,cnt,genes,g=genes[genes$lengp==3,],normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    #pl <- rbind(pl,cbind(profile=prot,genetype = ">2kb", mymetagenefun(prot,cnt,genes,g=genes[genes$lengp==4,],normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    
    pl <- rbind(pl,cbind(profile=prot,genetype = "Non-Intron genes", mymetagenefun(prot,cnt,genes,g=genes[!genes$tracking_id%in%introns$tracking_id,],normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    pl <- rbind(pl,cbind(profile=prot,genetype = "RP-Intron genes", mymetagenefun(prot,cnt,genes,g=introns[introns$tracking_id%in%rp$tracking_id,],normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    pl <- rbind(pl,cbind(profile=prot,genetype = "nonRP-Intron genes", mymetagenefun(prot,cnt,genes,g=introns[!introns$tracking_id%in%rp$tracking_id,],normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    
    
  }
  
  save(pl,file="Metagenplots_ChIPSpike_final_relative_deplbynondepl.RData")
  
  pdf("MetagenePlots_SpikeNorm_AllFactors_Merged_rel_deplbynondepl.pdf",width = 14,height = 8)
  myplot(pl[pl$genetype%in%c("OPEN genes","CLOSED genes"),]) + ggtitle("OPEN and CLOSED genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype%in%c("High Pol2","Medium Pol2","Low Pol2"),]) + ggtitle("Genes grouped by WT RNAPII occupancy")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype%in%c("<0.5kb","0.5kb-1kb","1kb-2kb",">2kb"),]) + ggtitle("Genes grouped by length")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype%in%c("non-DE genes","DE genes"),]) + ggtitle("Genes grouped expression status")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype%in%c("non-RP genes","RP genes"),]) + ggtitle("Genes grouped RP status")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  
  myplot(pl[pl$genetype%in%c("Intron genes","Non-Intron genes"),]) + ggtitle("Intron containing genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype%in%c("RP-Intron genes","nonRP-Intron genes"),]) + ggtitle("RP-Intron containing genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  
  dev.off()
  
}

## Spt6 depl/nondepl RNAPII normalized
### depl/non.depl for H3and RNAPII normalzied profiles
if(F){
  rm(list=ls())
  load(file="Gene groups.RData")
  vector.resizing <- function(x,final.len){
    y <- vector()
    len <- length(x)
    y <-spline(1:len,x,n=final.len)$y
    return(y)
  }
  
  res <- read.delim("Normalization.factors.merged.txt",header = T)
  myplot <- function(pl1){
    require(ggplot2)
    q <- ggplot(pl1, aes(x=pos,y=signal,ymax=signal+se, ymin=signal-se,col=genetype))
    q <- q + geom_hline(yintercept = c(1),col="purple")
    q <- q + geom_vline(xintercept = c(50,150),colour="gray50", linetype = "dashed")
    q <- q + geom_linerange(col='gray') + geom_line() #+ ylim(-0.03,0.15) #ylim(-0.01,0.045) #ylim(-0.04,0.1)
    q <- q + scale_color_manual(values = c("Depleted-1" = "pink","Depleted-2" = "red2",
                                           "Non-depleted-1" = "lightblue","Non-depleted-2" = "slateblue2",
                                           "Depleted" = "red","Non-depleted" = "slateblue",
                                           "OPEN genes"="green","CLOSED genes"="red",
                                           "High Pol2"="red4","Medium Pol2"="red","Low Pol2"="pink",
                                           "<0.5kb"="yellow","0.5kb-1kb"="orange","1kb-2kb"="red",">2kb"="red4",
                                           "non-DE genes"="pink","non-RP genes"="pink","DE genes"= "red", "RP genes"="red"))
    q <- q + facet_wrap(~profile, ncol=3,scales = "free")
    q <- q + scale_x_continuous(breaks = c(0,50,150,199), labels=c("","TSS","CPS","")) +ylab("average signal (depleted/non-depleted)") + xlab("") 
    q <- q + theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = NA),axis.line = element_line(size = 0.7, colour = "gray50"))
    q <- q + theme(axis.text.x = element_text(colour = 'black',size=14),
                   axis.text.y = element_text(colour = 'black',size=14),
                   axis.title = element_text(colour = 'black',size=15,face="bold"),
                   strip.text.x = element_text(colour = 'black',size=14,face="bold"))
    q <- q + theme(legend.text =element_text(colour = 'black',size=14),
                   legend.title = element_text(colour = 'black',size=14,face="bold"))
    q <- q + guides(color=guide_legend(title="Gene class"))
    return(q)
    
  }
  mymetagenefun <- function(prot,cnt,genes,g,normF=c(1,1),normF2=c(1,1)){
    f1 = paste("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_",prot,".RData",sep="")
    f2 = paste("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_",cnt,".RData",sep="")
    
    load(f1)
    l1 <- l
    load(f2)
    l2 <- l
    rm(l)
    
    
    depl1 <- matrix(NA,ncol=199,nrow = 5134)
    nondepl1 <- matrix(NA,ncol=199,nrow = 5134)
    for(i in 1:length(genes$tracking_id)){
      j = round(genes[i,]$width/10)
      v <- l1[[1]][i,51:(49+j)] + l1[[2]][i,51:(49+j)]
      v <- vector.resizing(v,100)
      v <- c(l1[[1]][i,1:50]+l1[[2]][i,1:50],
             v,
             l1[[1]][i,1525:1573]+l1[[2]][i,1525:1573]
      )
      depl1[i,] <- v
      
      v <- l1[[3]][i,51:(49+j)] + l1[[4]][i,51:(49+j)]
      v <- vector.resizing(v,100)
      v <- c(l1[[3]][i,1:50]+l1[[4]][i,1:50],
             v,
             l1[[3]][i,1525:1573]+l1[[4]][i,1525:1573]
      )
      nondepl1[i,] <- v
    }
    rownames(depl1) <- rownames(nondepl1) <- rownames(l1[[1]])
    depl1 <- subset(depl1,rownames(depl1)%in%g$tracking_id)
    nondepl1 <- subset(nondepl1,rownames(nondepl1)%in%g$tracking_id)
    
    depl2 <- matrix(NA,ncol=199,nrow = 5134)
    nondepl2 <- matrix(NA,ncol=199,nrow = 5134)
    for(i in 1:length(genes$tracking_id)){
      j = round(genes[i,]$width/10)
      v <- l2[[1]][i,51:(49+j)] + l2[[2]][i,51:(49+j)]
      v <- vector.resizing(v,100)
      v <- c(l2[[1]][i,1:50]+l2[[2]][i,1:50],
             v,
             l2[[1]][i,1525:1573]+l2[[2]][i,1525:1573]
      )
      depl2[i,] <- v
      
      v <- l2[[3]][i,51:(49+j)] + l2[[4]][i,51:(49+j)]
      v <- vector.resizing(v,100)
      v <- c(l2[[3]][i,1:50]+l2[[4]][i,1:50],
             v,
             l2[[3]][i,1525:1573]+l2[[4]][i,1525:1573]
      )
      nondepl2[i,] <- v
    }
    rownames(depl2) <- rownames(nondepl2) <- rownames(l2[[1]])
    depl2 <- subset(depl2,rownames(depl2)%in%g$tracking_id)
    nondepl2 <- subset(nondepl2,rownames(nondepl2)%in%g$tracking_id)
    
    depl1 <- depl1*normF[1]
    nondepl1 <- nondepl1*normF[2]
    depl2 <- depl2*normF2[1]
    nondepl2 <- nondepl2*normF2[2]
    
    #depl1[depl1==0] <- NA
    #depl2[depl2==0] <- NA
    #nondepl1[nondepl1==0] <- NA
    #nondepl2[nondepl2==0] <- NA
    
    depl1 <- (depl1*100+1)/(depl2*100+1)
    nondepl1 <- (nondepl1*100+1)/(nondepl2*100+1)
    depl1[!is.finite(depl1)] <-  NA
    nondepl1[!is.finite(nondepl1)] <- NA
    
    rat <- (depl1*100+1)/(nondepl1*100+1)
    rat[!is.finite(rat)] <-  NA
    
    get.avg <- function(m){
      df_m <- cbind(1:ncol(m),
                    as.data.frame(apply(m,2,function(x){
                      x <- x[x>quantile(x,0.05,na.rm=T) & x<quantile(x,0.99,na.rm=T)]
                      return(mean(x,na.rm=T))
                    })),
                    as.data.frame(apply(m,2,function(x){
                      x <- x[x>quantile(x,0.05,na.rm=T) & x<quantile(x,0.99,na.rm=T)]
                      return(sd(x,na.rm=T))
                    })), 
                    as.data.frame(apply(m,2, function(x){
                      x <- x[x>quantile(x,0.05,na.rm=T) & x<quantile(x,0.99,na.rm=T)]
                      return(sd(x,na.rm=TRUE)/sqrt(length(x)))
                    } ))
      )
      names(df_m) <- c("pos","signal", "sd","se")
      #df_m$signal <- df_m$signal*nF
      #df_m$sd <- df_m$sd*nF
      return(df_m)
    }
    
    pl <- cbind(get.avg(rat),condition="Depl./Non-depl.")
    return(pl)
  }
  
  #profiles <- c("H3K36me2","H3","H3K36me3","H3 (round1)",
  #              "H3K4me3","H3","S2P","RNAPII","S5P","RNAPII","Spt6","RNAPII (round1)")
  profiles <- c("Spt6","RNAPII (round1)","Spt6","RNAPII")
  pl <- c()
  for(i in seq(1,3,2)){
    prot = profiles[i]
    cnt = profiles[i+1]
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = "non-DE genes", mymetagenefun(prot,cnt,genes,g=genes[!genes$tracking_id%in%de.genes$tracking_id,],normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = "non-RP genes", mymetagenefun(prot,cnt,genes,g=genes[!genes$tracking_id%in%rp$tracking_id,],normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = "all genes", mymetagenefun(prot,cnt,genes,g=genes,normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = "RP genes", mymetagenefun(prot,cnt,genes,g=rp,normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = "Intron genes", mymetagenefun(prot,cnt,genes,g=introns,normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = "DE genes", mymetagenefun(prot,cnt,genes,g=de.genes,normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = "OPEN genes", mymetagenefun(prot,cnt,genes,g=opn,normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = "CLOSED genes", mymetagenefun(prot,cnt,genes,g=dpn,normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = "High Pol2", mymetagenefun(prot,cnt,genes,g=genes[genes$Pol2OccGp==3,],normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = "Medium Pol2", mymetagenefun(prot,cnt,genes,g=genes[genes$Pol2OccGp==2,],normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = "Low Pol2", mymetagenefun(prot,cnt,genes,g=genes[genes$Pol2OccGp==1,],normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = "<0.5kb", mymetagenefun(prot,cnt,genes,g=genes[genes$lengp==1,],normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = "0.5kb-1kb", mymetagenefun(prot,cnt,genes,g=genes[genes$lengp==2,],normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = "1kb-2kb", mymetagenefun(prot,cnt,genes,g=genes[genes$lengp==3,],normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = ">2kb", mymetagenefun(prot,cnt,genes,g=genes[genes$lengp==4,],normF = res[res$Factor==prot,]$normF,normF2 = res[res$Factor==cnt,]$normF)))
  }
  pl1 <- pl
  
  pl <- c()
  for(i in seq(1,3,2)){
    prot = profiles[i]
    cnt = profiles[i+1]
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = "non-DE genes", mymetagenefun(prot,cnt,genes,g=genes[!genes$tracking_id%in%de.genes$tracking_id,])))
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = "non-RP genes", mymetagenefun(prot,cnt,genes,g=genes[!genes$tracking_id%in%rp$tracking_id,])))
    
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = "all genes", mymetagenefun(prot,cnt,genes,g=genes)))
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = "RP genes", mymetagenefun(prot,cnt,genes,g=rp)))
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = "Intron genes", mymetagenefun(prot,cnt,genes,g=introns)))
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = "DE genes", mymetagenefun(prot,cnt,genes,g=de.genes)))
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = "OPEN genes", mymetagenefun(prot,cnt,genes,g=opn)))
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = "CLOSED genes", mymetagenefun(prot,cnt,genes,g=dpn)))
    
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = "High Pol2", mymetagenefun(prot,cnt,genes,g=genes[genes$Pol2OccGp==3,])))
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = "Medium Pol2", mymetagenefun(prot,cnt,genes,g=genes[genes$Pol2OccGp==2,])))
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = "Low Pol2", mymetagenefun(prot,cnt,genes,g=genes[genes$Pol2OccGp==1,])))
    
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = "<0.5kb", mymetagenefun(prot,cnt,genes,g=genes[genes$lengp==1,])))
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = "0.5kb-1kb", mymetagenefun(prot,cnt,genes,g=genes[genes$lengp==2,])))
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = "1kb-2kb", mymetagenefun(prot,cnt,genes,g=genes[genes$lengp==3,])))
    pl <- rbind(pl,cbind(profile=paste0(prot,"/",cnt),genetype = ">2kb", mymetagenefun(prot,cnt,genes,g=genes[genes$lengp==4,])))
  }
  save(pl,pl1,file="Spt6_Metagenplots_ChIPSpike_final_relative_deplbynondepl.RData")
  
  pdf("Spt6_MetagenePlots_SpikeNorm_AllFactors_Merged_rel_deplbynondepl.pdf",width = 10,height = 5)
  myplot(pl1[pl1$genetype%in%c("OPEN genes","CLOSED genes"),]) + ggtitle("OPEN and CLOSED genes (SpikeNormalized)")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl1[pl1$genetype%in%c("High Pol2","Medium Pol2","Low Pol2"),]) + ggtitle("Genes grouped by WT RNAPII occupancy (SpikeNormalized)")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl1[pl1$genetype%in%c("<0.5kb","0.5kb-1kb","1kb-2kb",">2kb"),]) + ggtitle("Genes grouped by length (SpikeNormalized)")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl1[pl1$genetype%in%c("non-DE genes","DE genes"),]) + ggtitle("Genes grouped expression status (SpikeNormalized)")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl1[pl1$genetype%in%c("non-RP genes","RP genes"),]) + ggtitle("Genes grouped RP status (SpikeNormalized)")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  
  myplot(pl[pl$genetype%in%c("OPEN genes","CLOSED genes"),]) + ggtitle("OPEN and CLOSED genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype%in%c("High Pol2","Medium Pol2","Low Pol2"),]) + ggtitle("Genes grouped by WT RNAPII occupancy")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype%in%c("<0.5kb","0.5kb-1kb","1kb-2kb",">2kb"),]) + ggtitle("Genes grouped by length")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype%in%c("non-DE genes","DE genes"),]) + ggtitle("Genes grouped expression status")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  myplot(pl[pl$genetype%in%c("non-RP genes","RP genes"),]) + ggtitle("Genes grouped RP status")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  dev.off()
  
}

## Log2FC comparison between ChIP-seq and RNA-seq
if(F){
  rm(list=ls())
  load("Gene groups.RData")
  load("Log2FC values for the ChIPFactors.RData")
  m0.spike$tracking_id <- rownames(m0.spike)
  m.spike$tracking_id <- rownames(m.spike)
  load(file="A:/work/WinstonLab/Natalia_Reim/Analysis_Feb18/DE_Genes_Filtered.RData")
  
  m0.spike <- merge(m0.spike,de.genes[,c("tracking_id","log2FoldChange","padj")],by="tracking_id",all.x=T)
  m.spike <- merge(m.spike,de.genes[,c("tracking_id","log2FoldChange","padj")],by="tracking_id",all.x=T)
  
  myplotfun <- function(m.spike,i,j){
    plot(m.spike[,i],m.spike[,j],main=names(m.spike)[i],xlab="log2 FC, factor occupancy", ylab="log2 FC,RNA",cex.lab=1.5,cex.axis=1.5,cex=1.2,cex.main=1.5,xaxt='n')
    axis(1,at =seq(-5,5,1),labels = seq(-5,5,1),tick = T ,cex=1.5,cex.axis=1.5)
    abline(h=seq(-5,5,1),v = seq(-5,5,1),col="gray",lty=2)
    mtext(paste("r=",round(cor(m.spike[,i],m.spike[,j],use="complete.obs"),3),sep=""),3,-1.5,col="red",cex=1.25)
    r <- lm(m.spike[,j]~m.spike[,i],na.action = na.omit)
    abline(r,col="red",lty=1)
  }
  
  pdf("FactorOccupancy and expression correlations.pdf",width = 12,height = 15)
  par(mfrow=c(4,3))
  for(i in 2:13){
    myplotfun(m.spike,i,14)
  }
  dev.off()
  
  pdf("FactorOccupancy and expression correlations_RPgenes.pdf",width = 12,height = 15)
  par(mfrow=c(4,3))
  for(i in 2:13){
    myplotfun(m.spike[m.spike$tracking_id%in%rp$tracking_id,],i,14)
  }
  dev.off()
  
  de.genes <- de.genes[de.genes$padj<0.05 & de.genes$log2FoldChange<0,]
  pdf("FactorOccupancy and expression correlations_downgenes.pdf",width = 12,height = 15)
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
  
  pdf("FactorOccupancy and expression sensitivity correlations.pdf",width = 12,height = 15)
  par(mfrow=c(4,3))
  for(i in 2:13){
    myplotfun(m.spike[m.spike$tracking_id%in%de.genes$tracking_id,],i,15)
  }
  dev.off()
  
  boxplot(m0.spike[,2:25])
  abline(h=0)
  ######## Which genes have decreased factor ocupancies?
  m <- m0.spike
  
  
  
}

## Differential binding of the factors
## I use DESeq2 to get pvalues
if(F){
  rm(list=ls())
  files <- Sys.glob("*.gr")
  
  load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenes.RData")
  genes <- subset(genes,genes$verification=="Verified"  & genes$chr!="chrM")
  genes <- as(genes,"GRanges")
  
  count <- c()
  for(f in files){
    cat(f,"\n")
    load(f)
    gr <-resize(gr,200)
    count <- cbind(count,countOverlaps(genes,gr,ignore.strand=T))
    rm(gr,f)
  }
  colnames(count) <- gsub(".gr","",files)
  cf <- cbind(as.data.frame(genes),count)
  save(cf,file="RawCounts_Signal_genesSGD.RData")
  
  write.table(cf,file="RawCountsperGene_ChIPAssay.txt",quote = F,sep = "\t",row.names = F)
  
  
  ################## NOTE
  ## 1) Orlando et al suggest computing spikein normalization factors by incorporating % of mixed spikein chromatin
  ## 2) Usually for RNASeq data is it considered constant
  ## 3) For ChcIPSeq data it cn be obtained from input sequencing
  ## 4) If DeSeq2 is used to calculate spikein size-factors, this contribution is lost
  ## 5) Howver, log2FC values for spikein genome center around 0 (as we want)
  ## 6) If I use Orlando et al spikein factors for DESeq2 computation, the log2FC values for spikein are thrown off in positive direction
  ## 7) It is unclear which is correct approach - as for metagene plots, I use earlier Orlando approach (which is more logical)
  ## 8) For computing log2FC values, i decide to use Deseq2 spikein size factor approach as this makes sure taht spikein log2FC values are centered around 0
  
  ## Approach:
  ## 1) Get raw ChIPseq signal along genes
  ## 2) RPM normalize the counts
  ## 3) Inflate teh counts to interger value by multiplying them by 10
  ## 4) Perform DESeq2 analysis using revised spikein normaization factors from MayNov_2017.txt
  
  rm(list=ls()) 
  load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenes.RData")
  genes <- subset(genes,genes$verification=="Verified" & genes$species=="Scer" & genes$chr!="chrM")
  
  load("A:/work/WinstonLab/Natalia_Reim/DataFrames/RawCounts_ChIPSignal_genesSGD.RData")
  load("A:/work/yeast_Annotations/refs.RData")
  genes <- merge(genes,unique(refs[,c("tracking_id","gene","link")]),by="tracking_id",all.x=T)
  rownames(cf) <- cf$tracking_id
  
  cf[,16:79] <- apply(cf[,16:79],2,function(x) round(x*1e7/sum(x)))  ## RPM normalize the counts
  m0 <- cf[cf$species=="Scer",16:79]
  mr <- cf[cf$species=="Spom",16:79]
  mr <- mr[rowSums(mr)>0,]
  m0 <- m0[rowSums(m0)>0,]
  

  res <- read.delim("Normalization.factors.NataliaMayNov2017_revised.txt",header=T)
  normF <- res$normF
  names(normF) <- res$sample
  normF <- 1/normF
  m0 <- m0[,names(normF)]
  mr <- mr[,names(normF)]
  
  
  
  #type=names(mr)
  #condition <- gsub("_[1234]","",type)
  #colData <- data.frame(condition,type)
  #row.names(colData) <- colData$type
  #dds <- DESeqDataSetFromMatrix(countData = as.matrix(mr),colData = colData,design = ~ condition)
  #sizeF <- estimateSizeFactors(dds)
  #normF = sizeF@colData@listData$sizeFactor
  #m0 <- m0[,names(normF)]
  #mr <- mr[,names(normF)]
  
  mydeseq2fun <- function(m,type,sizeF){
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
    
  }
  save(l,ln,file="Diff_Factor_Occupancy_DESeq2SGD.RData")
  
  
  load(file="Diff_Factor_Occupancy_DESeq2SGD.RData")
  pl <- list()
  pln <- list()
  for(i in 1:length(l)){
    pl[[paste(names(l)[i])]] <- l[[i]][,3]
    pln[[paste(names(ln)[i])]] <- ln[[i]][,3]
  }
  boxplot(pl,ylim=c(-3,3))
  boxplot(pln,ylim=c(-3,2))
  abline(h=0)
  #load("Diff_Factor_Occupancy_DESeq2SGD.RData")
  for(i in 1:12){
    write.table(ln[[i]],file= paste(names(ln)[i],"ChIP_OccupancyLog2FC.txt",sep=""),quote = F,sep = "\t",row.names = F)
  }
  
  m.spike <- ln[[1]]$`log2 (Spn1.Depl_H3_ChIP/Spn1_H3_ChIP)`
  for(i in 2:length(ln)){
    m.spike <- cbind(m.spike,ln[[i]][,3])
  }
  m.spike <- as.data.frame(m.spike)
  rownames(m.spike) <- ln[[1]]$tracking_id
  names(m.spike) <- names(ln)
  save(m.spike,file="Log2FC values for the ChIPFactors_DESeq2.RData")
  
  
  
  ## Comparing the DE expressed genes with factor occupancy changes
  rm(list=ls())
  load("Diff_Factor_Occupancy_DESeq2.RData")
  load(file="A:/work/WinstonLab/Natalia_Reim/Analysis_Feb18/DE_Genes_Filtered.RData")
  de.genes <- de.genes[de.genes$tracking_id!="YNL099C",]
  
  pdf("DiffFactor and RNA-expression_Vennplots.pdf",width = 10,height = 12)
  comgenes <- list()
  par(mfrow=c(4,3))
  for(i in 1:12){
    z <- l[[i]]
    z <- subset(z,!is.na(z$padj) & z$padj<0.05 & z$log2FoldChange<0)
    d <- FP22.Venn2(as.character(unique(de.genes[de.genes$padj<0.05 & de.genes$log2FoldChange<0,]$tracking_id)),
               as.character(unique(z$tracking_id)),names = c("RNA",names(l)[i]))
    comgenes[[paste(names(l)[i])]] <- d[grep(" ",d$condition),]$tracking_id
  }
  dev.off()
  ## Calculate Fusher P values
  load("Gene groups.RData")
  de.genes <- de.genes[de.genes$tracking_id!="YNL099C",]
  h <- quantile(genes$RNAPII,c(0.33,0.66))
  genes$Pol2OccGp <- 1
  genes$Pol2OccGp <- ifelse(genes$RNAPII>h[2],3,genes$Pol2OccGp)
  genes$Pol2OccGp <- ifelse(genes$RNAPII<h[2] & genes$RNAPII>h[1],2,genes$Pol2OccGp)
  genes$lengp <- 1
  genes$lengp <- ifelse(genes$width>2000,4,genes$lengp)
  genes$lengp <- ifelse(genes$width<2000 & genes$width>1000,3,genes$lengp)
  genes$lengp <- ifelse(genes$width<1000 & genes$width>500,2,genes$lengp)
  
  calc.fisher <- function(genes,cmn, gclass){
    d <- length(genes$tracking_id)
    b <- length(cmn)
    
    c <- sum(gclass%in%genes$tracking_id)
    d <- d-c
    a <- sum(gclass%in%cmn)
    b <- b-a
    p <- fisher.test(matrix(c(a,c,b,d),nrow=2,byrow = T),alternative = "greater")$p.value
    if(b>50){
      return(p)
    }else{
      return(1)
    }
  }
  
  myres <- matrix(NA,nrow = 12,ncol = 11)
  for(i in 1:12){
    myres[i,] <- c(calc.fisher(genes,cmn=as.character(comgenes[[i]]),gclass=as.character(introns$tracking_id)),
        calc.fisher(genes,cmn=as.character(comgenes[[i]]),gclass=as.character(rp$tracking_id)),
        calc.fisher(genes,cmn=as.character(comgenes[[i]]),gclass=as.character(opn$tracking_id)),
        calc.fisher(genes,cmn=as.character(comgenes[[i]]),gclass=as.character(dpn$tracking_id)),
    
        calc.fisher(genes,cmn=as.character(comgenes[[i]]),gclass=as.character(genes[genes$Pol2OccGp==1,]$tracking_id)),
        calc.fisher(genes,cmn=as.character(comgenes[[i]]),gclass=as.character(genes[genes$Pol2OccGp==2,]$tracking_id)),
        calc.fisher(genes,cmn=as.character(comgenes[[i]]),gclass=as.character(genes[genes$Pol2OccGp==3,]$tracking_id)),
    
        calc.fisher(genes,cmn=as.character(comgenes[[i]]),gclass=as.character(genes[genes$lengp==1,]$tracking_id)),
        calc.fisher(genes,cmn=as.character(comgenes[[i]]),gclass=as.character(genes[genes$lengp==2,]$tracking_id)),
        calc.fisher(genes,cmn=as.character(comgenes[[i]]),gclass=as.character(genes[genes$lengp==3,]$tracking_id)),
        calc.fisher(genes,cmn=as.character(comgenes[[i]]),gclass=as.character(genes[genes$lengp==4,]$tracking_id)))
  }
  colnames(myres) <- c("Intron","RP","OPEN","CLOSED","High-Pol2",'Medium-Pol2',"Low-Pol2",
                       "<0.5kb","0.5kb-1kb","1kb-2kb",'>2kb')
  rownames(myres) <- names(l)
  
  pdf("DEfactor_DERNA_common genes Pvalues.pdf",width = 15,height = 8)
  myres <- formatC(myres, format = "e", digits = 2)
  grid.table(myres)
  dev.off()
}

## Log2FC comparison between ChIP-seq and RNA-seq with DESeq2
if(F){
  rm(list=ls())
  load("Gene groups.RData")
  load("Log2FC values for the ChIPFactors_DESeq2.RData")
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
  
  pdf("FactorOccupancy and expression correlations_DESeq2.pdf",width = 12,height = 15)
  par(mfrow=c(4,3))
  for(i in 2:13){
    myplotfun(m.spike,i,14)
  }
  dev.off()
  
  pdf("FactorOccupancy and expression correlations_RPgenes_DESeq2.pdf",width = 12,height = 15)
  par(mfrow=c(4,3))
  for(i in 2:13){
    myplotfun(m.spike[m.spike$tracking_id%in%rp$tracking_id,],i,14)
  }
  dev.off()
  
  de.genes <- de.genes[de.genes$padj<0.05 & de.genes$log2FoldChange<0,]
  pdf("FactorOccupancy and expression correlations_downgenes_DESeq2.pdf",width = 12,height = 15)
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
  
  pdf("FactorOccupancy and expression sensitivity correlations_DESeq2.pdf",width = 12,height = 15)
  par(mfrow=c(4,3))
  for(i in 2:13){
    myplotfun(m.spike[m.spike$tracking_id%in%de.genes$tracking_id,],i,15)
  }
  dev.off()
  
}


## Pol2 occupancy vs expression in WT
if(F){
  rm(list=ls())
  load("A:/work/WinstonLab/Natalia_Reim/Analysis_Feb18/Sense_Antisese_ReadCounts_revised.RData")
  load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenes.RData")
  s <- subset(s,rownames(s)%in%genes[genes$verification=="Verified",]$tracking_id)
  as <- subset(as,rownames(as)%in%genes[genes$verification=="Verified",]$tracking_id)
  mytempfun <- function(df1,fdr=0.05,remove.spike=T){
    if(remove.spike==T){
      spiked.chr <- c("chrXVII","chrXVIII","chrXIX","chrMsp","chrM","Sp_chrAB325691")
      df1.e <- subset(df1, !df1$chr%in%spiked.chr)
    }else{
      df1.e <- df1
    }
    df1.e <-as.matrix( round (df1.e[,2:5] ))
    type <- colnames(df1.e)
    condition <- sub("-[012]","",type)
    colData <- data.frame(condition,type)
    row.names(colData) <- colData$type
    dds <- DESeqDataSetFromMatrix(countData = df1.e,colData = colData,design = ~ condition)
    sizeFactors(dds) <- dd[1:4]
    dds <- DESeq(dds)
    res1 <- as.data.frame(results(dds))
    plot(log(res1$baseMean),res1$log2FoldChange,pch=20,main="DESeq2 MA plot",xlab="log mean-expression", ylab="log2 FC",xlim=c(-5,15),cex.axis=1.5,cex.lab=1.5)
    res1 <- res1[res1$padj<fdr,]
    res1 <- res1[!is.na(res1$padj),]
    points(log(res1$baseMean),res1$log2FoldChange,pch=20,col="red")
    return(dds)
  }
  s <- s[complete.cases(s),]
  dds <- mytempfun(df1=s[,c(13,1:4)],fdr=0.05)
  dds <- DESeq2::counts(dds,normalized=T)
  dds <- as.data.frame(dds)
  dds$Expn <- apply(dds[,1:2],1,function(x) round(sqrt((prod(x))),2))
  dds$tracking_id <- rownames(dds)
  dds <- dds[,5:6]
  rm(as,s,genes,ad,dd)
  
  load(file="Gene groups.RData")
  genes <- merge(genes,dds,by="tracking_id",all.x=T)
  genes <- genes[genes$RNAPII>25,]

  pdf("RNAPII occupancy and expression correlation.pdf",5,5)
  plot(log2(genes$RNAPII+1),log2(genes$Expn+1),col=rgb(1,0,0,0.2),xlab="RNAPII occupancy, log2", ylab="Gene expression, log2")
  abline(h = c(5,10,15),v = c(5,10,15),col="gray",lty=2)
  r <- lm(log2(genes$Expn+1)~log2(genes$RNAPII+1),na.action = na.omit)
  abline(r,col="red4",lwd=2)
  mtext(paste("r=",round(cor(log2(genes$RNAPII+1),log2(genes$Expn+1),use="complete.obs"),2),sep=""),3,-1.5,cex=1.5)
  dev.off()
}


### Factor occupancy change quantification using boxplot
if(F){
  rm(list=ls())
  load("Gene groups.RData")
  load("Log2FC values for the ChIPFactors_NormOccup.RData")
  
  res<- read.delim("Normalization.factors.merged.txt",header=T)
  m.spike <- m.spike[,-grep("Inp",names(m.spike))]
  m0.spike <- m0.spike[,-grep("Inp",names(m0.spike))]
  
  names(m.spike) <- c("Depl.H3 (round1)", "Depl.H3", "Depl.H3K36me2", 
                      "Depl.H3K36me3", "Depl.H3K4me3", "Depl.Set2", 
                      "Depl.Spt6", "Depl.RNAPII (round1)", "Depl.RNAPII", 
                      "Depl.S2P", "Depl.S5P", "Depl.Spn1", 
                      "H3 (round1)", "H3", "H3K36me2", "H3K36me3", 
                      "H3K4me3", "Set2", "Spt6", 
                      "RNAPII", "RNAPII", "S2P", 
                      "S5P", "Spn1")
  
  fact <- as.character(unique(res$Factor))
  fact <- fact[-7]
  
  pdf("Factor_occupancy changes boxplots.pdf",width = 8,height = 12)
  par(mfrow=c(4,3))
  for(i in 1:12){
    #boxplot(log2(m.spike[,grep(f,names(m.spike))]+1),outline=F)
    cat(names(m.spike)[i],"\t",names(m.spike)[i+12],"\n")
    z <- log2(m.spike[,c(i,(i+12))]+1)
    names(z) <- c("Depleted","Non-depleted")
    pval <- formatC(t.test(z[,1],z[,2])$p.value,format="e",digits = 2)
    boxplot(z,ylab="Factor occupancy,log2",main=paste0(names(m.spike)[i+12]),col=c("steelblue","lightblue"),ylim=c(min(z[,1],z[,2],na.rm=T),max(z[,1],z[,2],na.rm=T)+1),cex.axis=1.5,cex.lab=1.5,cex.main=2)
    mtext(paste0("T-test, p = ", pval),3,-1,adj = 1,cex=1.3,col="red")
  }
  dev.off()
}

### Genes overlap
if(F){
  rm(list=ls())
  load(file="Gene groups.RData")
  
  pdf("GeneOlaps with DE genes.pdf",width = 8,height = 10)
  par(mfrow=c(3,2))
  FP22.Venn2(as.character(de.genes$tracking_id),as.character(dpn$tracking_id),names = c("DownGenes","CLOSED-prom"))
  FP22.Venn2(as.character(de.genes$tracking_id),as.character(opn$tracking_id),names = c("DownGenes","OPEN-prom"))
  FP22.Venn2(as.character(de.genes$tracking_id),as.character(introns$tracking_id),names = c("DownGenes","Intron-cont."))
  FP22.Venn2(as.character(de.genes$tracking_id),as.character(rp$tracking_id),names = c("DownGenes","RP-genes"))
  
  FP22.Venn2(as.character(de.genes$tracking_id),as.character(introns[introns$tracking_id%in%rp$tracking_id,]$tracking_id),names = c("DownGenes","Intron-RP"))
  dev.off()
  
  
}
## 

######################################################################################
### TSS centered heatmap sorted by change in Spt6
######################################################################################
if(F){
  rm(list = ls())
  library(EnrichedHeatmap)
  library(circlize)
  library(ComplexHeatmap)
  
  load("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/RawCounts_Signal_genes.RData")
  cf <- cf[cf$species=="Scer",]
  cf <- cf[,c(6,65:68)]
  cf$RNAPII <- rowSums(cf[,2:5])
  cf <- cf[,c(1,6)]
  load(file="A:/work/WinstonLab/Natalia_Reim/Analysis_Feb18/DE_Genes_Filtered.RData")
  load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenes.RData")
  genes <- subset(genes,genes$verification=="Verified" & genes$species=="Scer" & genes$chr!="chrM")
  genes <- merge(genes,cf,by="tracking_id",all.x=T)
  genes <- merge(genes,de.genes[,c("tracking_id","log2FoldChange")],by="tracking_id",all.x=T)
  genes$RNAPII <- ifelse(is.na(genes$RNAPII),0,genes$RNAPII)
  genes$log2FoldChange <- ifelse(is.na(genes$log2FoldChange),0,genes$log2FoldChange)
  genes$fpkm <- round(genes$RNAPII*1000/genes$width,2)
 
  g <- as(genes,"GRanges")
  d <- as.data.frame(findOverlaps(g,ignore.strand=T))
  d <- d[d$queryHits!=d$subjectHits,] 
  d <- unique(d$queryHits,d$subjectHits)
  
  myheatplotfun <- function(path,normF=c(1,1)){
    load(path)
    a1 <- l[["Depleted"]]
    a2 <- l[["Non-depleted"]]
    a1 <- a1*normF[1]
    a2 <- a2*normF[2]
    get_divsion <- function(m,n){
      m <- m*100+1
      n <- n*100+1
      m <- (m-n)/(m+n)
      return(m)
    }
    
    #a <- log2((a1+1)/(a2+1))
    #a[!is.finite(a)] <-NA
    #range(a,na.rm=T)
    a <- get_divsion(a1,a2)
    a[!is.finite(a)] <-0
    a <- a[-d,]
    
    dim(a) = dim(a)
    attr(a, "upstream_index") = 1:100
    attr(a, "target_index") = integer(0)
    attr(a, "downstream_index") = 101:200
    attr(a, "extend") = c(1000, 1000)  # it must be a vector of length two
    class(a) = c("normalizedMatrix", "matrix")
    attr(a, "signal_name") = "Spt6"
    attr(a, "target_name") = "TSS"
    attr(a,"target_is_single_point") = FALSE
    return(a)
  }
  
  res <- read.delim("SpikeinNormalization_factors_revised.txt",header = T)

  spt6 <- myheatplotfun("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_TSScomb_Spt6.RData",res[res$Factor=="Spt6",]$normF)
  spn1 <- myheatplotfun("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_TSScomb_Spn1.RData",res[res$Factor=="Spn1",]$normF)
  k36me3 <- myheatplotfun("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_TSScomb_H3K36me3.RData",res[res$Factor=="H3K36me3",]$normF)
  k36me2 <- myheatplotfun("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_TSScomb_H3K36me2.RData",res[res$Factor=="H3K36me2",]$normF)
  k4me3 <- myheatplotfun("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_TSScomb_H3K4me3.RData",res[res$Factor=="H3K4me3",]$normF)
  h3 <- myheatplotfun("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_TSScomb_H3.RData",res[res$Factor=="H3",]$normF)
  pol2 <- myheatplotfun("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_TSScomb_RNAPII.RData",res[res$Factor=="RNAPII",]$normF)
  set2 <- myheatplotfun("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_TSScomb_Set2.RData",res[res$Factor=="Set2",]$normF)
  s2p <- myheatplotfun("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_TSScomb_S2P.RData",res[res$Factor=="S2P",]$normF)
  s5p <- myheatplotfun("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_TSScomb_S5P.RData",res[res$Factor=="S5P",]$normF)
  
  load("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_Spt6.RData")
  z <- data.frame(depl=apply(l[[1]][,51:1523],1,sum)+apply(l[[2]][,51:1523],1,sum),
             nondepl=apply(l[[3]][,51:1523],1,sum)+apply(l[[4]][,51:1523],1,sum))
  z$o <- (z$depl-z$nondepl)/(z$depl+z$nondepl)
  z <- z[-d,]
  o <- order(z$o,decreasing = T)
  
  
  load("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_H3K36me3.RData")
  z <- data.frame(depl=apply(l[[1]][,51:1523],1,sum)+apply(l[[2]][,51:1523],1,sum),
                  nondepl=apply(l[[3]][,51:1523],1,sum)+apply(l[[4]][,51:1523],1,sum))
  z$tracking_id<- rownames(z)
  z <- merge(z,genes[,c("tracking_id","width")],by="tracking_id",all.x=T)
  z$o <- (z$depl-z$nondepl)/(z$depl+z$nondepl)
  z <- z[-d,]
  o1 <- order(z$o,decreasing = T)
  
  load("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_TSScomb_Spt6.RData")
  z <- data.frame(depl = apply(l[[1]][,101:200],1,sum),
                  nondepl = apply(l[[2]][,101:200],1,sum))
  z$o <- (z$depl-z$nondepl)/(z$depl+z$nondepl)
  z <- z[-d,]
  o2 <- order(z$o,decreasing = T)
  
  
  set.seed(123)
  #col_fun = colorRamp2(c(-1,0,1), c("#2166ac","#f7f7f7","#b2182b"))
  col_fun_spt6 = colorRamp2(c(min(spt6),0, max(spt6)), c("blue4","white","red4"))
  col_fun_pol2 = colorRamp2(c(min(pol2),0, max(pol2)), c("blue4","white","red4"))
  col_fun_s2p = colorRamp2(c(min(s2p),0, max(s2p)), c("blue4","white","red4"))
  col_fun_s5p = colorRamp2(c(min(s5p),0, max(s5p)), c("blue4","white","red4"))
  col_fun_spn1 = colorRamp2(c(min(spn1), 0,max(spn1)), c("blue4","white","red4"))
  col_fun_k36me3 = colorRamp2(c(min(k36me3),0, max(k36me3)),  c("blue4","white","red4"))
  col_fun_k36me2 = colorRamp2(c(min(k36me2),0, max(k36me2)),  c("blue4","white","red4"))
  col_fun_k4me3 = colorRamp2(c(min(k4me3),0, max(k4me3)),  c("blue4","white","red4"))
  col_fun_h3 = colorRamp2(c(min(h3),0, max(h3)),  c("blue4","white","red4"))
  col_fun_set2 = colorRamp2(c(min(set2),0, max(set2)),  c("blue4","white","red4"))
  
  
  ## ordering by Spt6 change
    o2 <- order(log2(genes[-d,]$fpkm+1))
    
    ht_list = EnrichedHeatmap(spt6,row_order = o2, col = col_fun_spt6, name = "Spt6",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                    column_title = "Spt6") +
    EnrichedHeatmap(pol2, col = col_fun_pol2,row_order = o2, name = "RNAPII",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                    column_title = "RNAPII") +
    EnrichedHeatmap(s2p, col = col_fun_s2p,row_order = o2, name = "Ser2-P",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                    column_title = "Ser2-P") +
    EnrichedHeatmap(s5p, col = col_fun_s5p,row_order = o2, name = "Ser5-P",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                    column_title = "Ser5-P") +
    EnrichedHeatmap(k36me3, col = col_fun_k36me3,row_order = o2, name = "H3K36me3",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                    column_title = "H3K36me3") +
    EnrichedHeatmap(k36me2, col = col_fun_k36me2,row_order = o2, name = "H3K36me2",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                    column_title = "H3K36me2") +
    EnrichedHeatmap(k4me3, col = col_fun_k4me3,row_order = o2, name = "H3K4me3",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                    column_title = "H3K4me3") +
    EnrichedHeatmap(set2, col = col_fun_set2,row_order = o2, name = "Set2",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                    column_title = "Set2") +
    EnrichedHeatmap(h3, col = col_fun_h3,row_order = o2, name = "H3",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                    column_title = "H3") +
    EnrichedHeatmap(spn1, col = col_fun_spn1,row_order = o2, name = "Spn1",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4))), 
                    column_title = "Spn1") +
    Heatmap(log2(genes[-d,]$fpkm+1),row_order = o2, col = c("white","lightpink","pink", "red4"), name = "log2(fpkm+1)", 
            show_row_names = FALSE, width = unit(4, "mm"))
    png("Spt6 sorted_signal_DeplNonDepl_Expression.png",width = 1500,height = 800)
    draw(ht_list,heatmap_legend_side = "bottom", gap = unit(2, "mm"))
    dev.off()
  

    ## ordering by H3K36me3 change
    load("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_TSScomb_H3K36me3.RData")
    z <- data.frame(depl = apply(l[[1]][,101:200],1,sum),
                    nondepl = apply(l[[2]][,101:200],1,sum))
    z$o <- (z$depl-z$nondepl)/(z$depl+z$nondepl)
    z <- z[-d,]
    o3 <- order(z$o,decreasing = T)

    ht_list = EnrichedHeatmap(k36me3, col = col_fun_k36me3,row_order = o3, name = "H3K36me3",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                      column_title = "H3K36me3") +
      EnrichedHeatmap(k36me2, col = col_fun_k36me2,row_order = o3, name = "H3K36me2",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                      column_title = "H3K36me2") +
      EnrichedHeatmap(k4me3, col = col_fun_k4me3,row_order = o3, name = "H3K4me3",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                      column_title = "H3K4me3") +
      EnrichedHeatmap(set2, col = col_fun_set2,row_order = o3, name = "Set2",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                      column_title = "Set2") +
      EnrichedHeatmap(h3, col = col_fun_h3,row_order = o3, name = "H3",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                      column_title = "H3") +
      EnrichedHeatmap(spt6,row_order = o3, col = col_fun_spt6, name = "Spt6",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                      column_title = "Spt6") +
      EnrichedHeatmap(pol2, col = col_fun_pol2,row_order = o3, name = "RNAPII",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                      column_title = "RNAPII") +
      EnrichedHeatmap(s2p, col = col_fun_s2p,row_order = o3, name = "Ser2-P",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                      column_title = "Ser2-P") +
      EnrichedHeatmap(s5p, col = col_fun_s5p,row_order = o3, name = "Ser5-P",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                      column_title = "Ser5-P") +
      EnrichedHeatmap(spn1, col = col_fun_spn1,row_order = o3, name = "Spn1",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4))), 
                      column_title = "Spn1") +
      Heatmap(log2(genes[-d,]$fpkm+1),row_order = o3, col = c("white","lightpink","pink", "red4"), name = "log2(fpkm+1)", 
              show_row_names = FALSE, width = unit(4, "mm"))
    png("H3K36me3 sorted_signal_DeplNonDepl_1kb.png",width = 1500,height = 800)
    draw(ht_list,heatmap_legend_side = "bottom", gap = unit(2, "mm"))
    dev.off()
    
      
  
  
  
  ## visualize (3 clusters)
  partition = paste0("cluster", kmeans(spt6, centers = 3)$cluster)
  lgd = Legend(at = c("cluster1", "cluster2", "cluster3"), title = "Clusters", 
               type = "lines", legend_gp = gpar(col = 2:4))
  ht_list = Heatmap(partition, col = structure(2:4, names = paste0("cluster", 1:3)), name = "",
                    show_row_names = FALSE, width = unit(1, "mm")) +
            EnrichedHeatmap(spt6, col = col_fun_spt6, name = "Spt6",
                    top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                    column_title = "Spt6") + 
            EnrichedHeatmap(pol2, col = col_fun_pol2, name = "RNAPII",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                    column_title = "RNAPII") +
            EnrichedHeatmap(s2p, col = col_fun_s2p, name = "Ser2-P",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                    column_title = "Ser2-P") +
            EnrichedHeatmap(s5p, col = col_fun_s5p, name = "Ser5-P",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                    column_title = "Ser5-P") +
            EnrichedHeatmap(k36me3, col = col_fun_k36me3, name = "H3K36me3",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                    column_title = "H3K36me3") +
            EnrichedHeatmap(k36me2, col = col_fun_k36me2, name = "H3K36me2",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                    column_title = "H3K36me2") +
            EnrichedHeatmap(k4me3, col = col_fun_k4me3, name = "H3K4me3",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                    column_title = "H3K4me3") +
            EnrichedHeatmap(set2, col = col_fun_set2, name = "Set2",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                    column_title = "Set2") +
            EnrichedHeatmap(h3, col = col_fun_h3, name = "H3",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                    column_title = "H3") +
            EnrichedHeatmap(spn1, col = col_fun_spn1, name = "Spn1",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4))), 
                    column_title = "Spn1") +
            Heatmap(log2(genes$fpkm+1), col = c("white","lightpink","pink", "red4"), name = "log2(fpkm+1)", 
                     show_row_names = FALSE, width = unit(4, "mm"))
  
  png("Spt6 sorted_signal_DeplNonDepl_3clusters.png",width = 1500,height = 800)
  draw(ht_list, split = partition, annotation_legend_list = list(lgd)
         , heatmap_legend_side = "bottom", gap = unit(2, "mm"))
  dev.off()
  
  
  ### 5 clusters
  partition = paste0("cluster", kmeans(spt6, centers = 4)$cluster)
  lgd = Legend(at = c("cluster1", "cluster2", "cluster3","cluster4"), title = "Clusters", 
               type = "lines", legend_gp = gpar(col = 2:5))
  ht_list = Heatmap(partition, col = structure(2:5, names = paste0("cluster", 1:4)), name = "",
                    show_row_names = FALSE, width = unit(1, "mm")) +
    EnrichedHeatmap(spt6, col = col_fun_spt6, name = "Spt6",
                    top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                    column_title = "Spt6") + 
    EnrichedHeatmap(pol2, col = col_fun_pol2, name = "RNAPII",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                    column_title = "RNAPII") +
    EnrichedHeatmap(s2p, col = col_fun_s2p, name = "Ser2-P",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                    column_title = "Ser2-P") +
    EnrichedHeatmap(s5p, col = col_fun_s5p, name = "Ser5-P",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                    column_title = "Ser5-P") +
    EnrichedHeatmap(k36me3, col = col_fun_k36me3, name = "H3K36me3",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                    column_title = "H3K36me3") +
    EnrichedHeatmap(k36me2, col = col_fun_k36me2, name = "H3K36me2",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                    column_title = "H3K36me2") +
    EnrichedHeatmap(k4me3, col = col_fun_k4me3, name = "H3K4me3",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                    column_title = "H3K4me3") +
    EnrichedHeatmap(set2, col = col_fun_set2, name = "Set2",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                    column_title = "Set2") +
    EnrichedHeatmap(h3, col = col_fun_h3, name = "H3",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                    column_title = "H3") +
    EnrichedHeatmap(spn1, col = col_fun_spn1, name = "Spn1",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5))), 
                    column_title = "Spn1") +
    Heatmap(log2(genes$fpkm+1), col = c("white","lightpink","pink", "red4"), name = "log2(fpkm+1)", 
            show_row_names = FALSE, width = unit(4, "mm"))
  
  png("Spt6 sorted_signal_DeplNonDepl_4clusters.png",width = 1500,height = 800)
  draw(ht_list, split = partition, annotation_legend_list = list(lgd)
       , heatmap_legend_side = "bottom", gap = unit(2, "mm"))
  dev.off()
  
  cf1 <- data.frame(cluster=partition,gene=rownames(m))
  
  
  
  ## sorting based on H3K36me3 change
  partition = paste0("cluster", kmeans(k36me3, centers = 4)$cluster)
  lgd = Legend(at = c("cluster1", "cluster2", "cluster3","cluster4"), title = "Clusters", 
               type = "lines", legend_gp = gpar(col = 2:5))
  ht_list = Heatmap(partition, col = structure(2:5, names = paste0("cluster", 1:4)), name = "",
                    show_row_names = FALSE, width = unit(1, "mm")) +
    
    EnrichedHeatmap(k36me3, col = col_fun_k36me3, name = "H3K36me3",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                    column_title = "H3K36me3") +
    EnrichedHeatmap(k36me2, col = col_fun_k36me2, name = "H3K36me2",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                    column_title = "H3K36me2") +
    EnrichedHeatmap(k4me3, col = col_fun_k4me3, name = "H3K4me3",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                    column_title = "H3K4me3") +
    EnrichedHeatmap(set2, col = col_fun_set2, name = "Set2",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                    column_title = "Set2") +
    EnrichedHeatmap(h3, col = col_fun_h3, name = "H3",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                    column_title = "H3") +
    
    EnrichedHeatmap(spt6, col = col_fun_spt6, name = "Spt6",
                    top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                    column_title = "Spt6") + 
    EnrichedHeatmap(pol2, col = col_fun_pol2, name = "RNAPII",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                    column_title = "RNAPII") +
    EnrichedHeatmap(s2p, col = col_fun_s2p, name = "Ser2-P",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                    column_title = "Ser2-P") +
    EnrichedHeatmap(s5p, col = col_fun_s5p, name = "Ser5-P",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                    column_title = "Ser5-P") +
    
    EnrichedHeatmap(spn1, col = col_fun_spn1, name = "Spn1",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5))), 
                    column_title = "Spn1") +
    Heatmap(log2(genes$fpkm+1), col = c("white","lightpink","pink", "red4"), name = "log2(fpkm+1)", 
            show_row_names = FALSE, width = unit(4, "mm"))
  
  png("H3K36me3 sorted_signal_DeplNonDepl_4clusters.png",width = 1500,height = 800)
  draw(ht_list, split = partition, annotation_legend_list = list(lgd)
       , heatmap_legend_side = "bottom", gap = unit(2, "mm"))
  dev.off()
  
  cf2 <- data.frame(cluster=partition,gene=rownames(m))
  
  
  
  cf1 <- merge(cf1,genes,by.x="gene",by.y="tracking_id",all.x=T)  
  cf2 <- merge(cf2,genes,by.x="gene",by.y="tracking_id",all.x=T)  
  
  write.table(cf1,file="Spt6_4Clusters.txt",quote = F,sep = "\t",row.names = F)
  write.table(cf2,file="H3K36me3_4Clusters.txt",quote = F,sep = "\t",row.names = F)
  save(cf1,cf2,file="Spt6_H3K36me3Clusters.RData")
  
  ## 
  rm(list=ls())
  load(file="Gene groups.RData")
  load(file="Spt6_H3K36me3Clusters.RData")
  
  cf1$rp <- ifelse(cf1$gene%in%rp$tracking_id,1,0)
  cf1$introns <- ifelse(cf1$gene%in%introns$tracking_id,1,0)
  cf1$de <- ifelse(cf1$gene%in%de.genes$tracking_id,1,0)
  cf1$opn <- ifelse(cf1$gene%in%opn$tracking_id,1,0)
  cf1$dpn <- ifelse(cf1$gene%in%dpn$tracking_id,1,0)
  
  res <- c()
  
  x <- cbind(table(cf1[cf1$rp==1,]$cluster),
  table(cf1$cluster)-table(cf1[cf1$rp==1,]$cluster),
  round(table(cf1$cluster)*(181/5134)),
  table(cf1$cluster)-round(table(cf1$cluster)*(181/5134)))
  x <- as.data.frame(x)
  res <-  rbind(res,data.frame(category = "RPgenes" ,apply(x, 1, function(x){
       fisher.test(matrix(x,nrow = 2,byrow = T),alternative = "greater")$p.value
     })))
  
  x <- cbind(table(cf1[cf1$introns==1,]$cluster),
             table(cf1$cluster)-table(cf1[cf1$introns==1,]$cluster),
             round(table(cf1$cluster)*(249/5134)),
             table(cf1$cluster)-round(table(cf1$cluster)*(249/5134)))
  x <- as.data.frame(x)
  res <-  rbind(res,data.frame(category = "Introngenes" ,apply(x, 1, function(x){
    fisher.test(matrix(x,nrow = 2,byrow = T),alternative = "greater")$p.value
  })))
  
  x <- cbind(table(cf1[cf1$de==1,]$cluster),
             table(cf1$cluster)-table(cf1[cf1$de==1,]$cluster),
             round(table(cf1$cluster)*(1243/5134)),
             table(cf1$cluster)-round(table(cf1$cluster)*(1243/5134)))
  x <- as.data.frame(x)
  res <-  rbind(res,data.frame(category = "DEgenes" ,apply(x, 1, function(x){
    fisher.test(matrix(x,nrow = 2,byrow = T),alternative = "greater")$p.value
  })))
  
  x <- cbind(table(cf1[cf1$opn==1,]$cluster),
             table(cf1$cluster)-table(cf1[cf1$opn==1,]$cluster),
             round(table(cf1$cluster)*(2597/5134)),
             table(cf1$cluster)-round(table(cf1$cluster)*(2597/5134)))
  x <- as.data.frame(x)
  res <-  rbind(res,data.frame(category = "OPENgenes" ,apply(x, 1, function(x){
    fisher.test(matrix(x,nrow = 2,byrow = T),alternative = "greater")$p.value
  })))
  
  x <- cbind(table(cf1[cf1$dpn==1,]$cluster),
             table(cf1$cluster)-table(cf1[cf1$dpn==1,]$cluster),
             round(table(cf1$cluster)*(1030/5134)),
             table(cf1$cluster)-round(table(cf1$cluster)*(1030/5134)))
  x <- as.data.frame(x)
  res <-  rbind(res,data.frame(category = "CLOSEDgenes" ,apply(x, 1, function(x){
    fisher.test(matrix(x,nrow = 2,byrow = T),alternative = "greater")$p.value
  })))
  

  ## cf2
  cf2$rp <- ifelse(cf2$gene%in%rp$tracking_id,1,0)
  cf2$introns <- ifelse(cf2$gene%in%introns$tracking_id,1,0)
  cf2$de <- ifelse(cf2$gene%in%de.genes$tracking_id,1,0)
  cf2$opn <- ifelse(cf2$gene%in%opn$tracking_id,1,0)
  cf2$dpn <- ifelse(cf2$gene%in%dpn$tracking_id,1,0)
  
  res2 <- c()
  
  x <- cbind(table(cf2[cf2$rp==1,]$cluster),
             table(cf2$cluster)-table(cf2[cf2$rp==1,]$cluster),
             round(table(cf2$cluster)*(181/5134)),
             table(cf2$cluster)-round(table(cf2$cluster)*(181/5134)))
  x <- as.data.frame(x)
  res2 <-  rbind(res2,data.frame(category = "RPgenes" ,apply(x, 1, function(x){
    fisher.test(matrix(x,nrow = 2,byrow = T),alternative = "greater")$p.value
  })))
  
  x <- cbind(table(cf2[cf2$introns==1,]$cluster),
             table(cf2$cluster)-table(cf2[cf2$introns==1,]$cluster),
             round(table(cf2$cluster)*(249/5134)),
             table(cf2$cluster)-round(table(cf2$cluster)*(249/5134)))
  x <- as.data.frame(x)
  res2 <-  rbind(res2,data.frame(category = "Introngenes" ,apply(x, 1, function(x){
    fisher.test(matrix(x,nrow = 2,byrow = T),alternative = "greater")$p.value
  })))
  
  x <- cbind(table(cf2[cf2$de==1,]$cluster),
             table(cf2$cluster)-table(cf2[cf2$de==1,]$cluster),
             round(table(cf2$cluster)*(1243/5134)),
             table(cf2$cluster)-round(table(cf2$cluster)*(1243/5134)))
  x <- as.data.frame(x)
  res2 <-  rbind(res2,data.frame(category = "DEgenes" ,apply(x, 1, function(x){
    fisher.test(matrix(x,nrow = 2,byrow = T),alternative = "greater")$p.value
  })))
  
  x <- cbind(table(cf2[cf2$opn==1,]$cluster),
             table(cf2$cluster)-table(cf2[cf2$opn==1,]$cluster),
             round(table(cf2$cluster)*(2597/5134)),
             table(cf2$cluster)-round(table(cf2$cluster)*(2597/5134)))
  x <- as.data.frame(x)
  res2 <-  rbind(res2,data.frame(category = "OPENgenes" ,apply(x, 1, function(x){
    fisher.test(matrix(x,nrow = 2,byrow = T),alternative = "greater")$p.value
  })))
  
  x <- cbind(table(cf2[cf2$dpn==1,]$cluster),
             table(cf2$cluster)-table(cf2[cf2$dpn==1,]$cluster),
             round(table(cf2$cluster)*(1030/5134)),
             table(cf2$cluster)-round(table(cf2$cluster)*(1030/5134)))
  x <- as.data.frame(x)
  res2 <-  rbind(res2,data.frame(category = "CLOSEDgenes" ,apply(x, 1, function(x){
    fisher.test(matrix(x,nrow = 2,byrow = T),alternative = "greater")$p.value
  })))
  
    
  
  
  names(res)[2] <- names(res2)[2] <- 'pvalue'
  res$cluster <- res2$cluster <- rep(c("Cluster1","Cluster2","Cluster3","Cluster4"),5)
  
  library(ggplot2)
  q <- ggplot(res) + geom_bar(aes(x=cluster,y=-log10(pvalue),fill=cluster),stat="identity",position="dodge") 
  q <- q + xlab("") + ylab("-log10(Pvalue)") +facet_wrap(~category,nrow = 2)  
  q <- q + scale_color_manual(values = c("Cluster1"="red","Cluster2"="green","Cluster3"="blue","Cluster4"="cyan"))
  q 
  
  q <- ggplot(res2) + geom_bar(aes(x=cluster,y=-log10(pvalue),fill=cluster),stat="identity",position="dodge") 
  q <- q + xlab("") + ylab("-log10(Pvalue)") +facet_wrap(~category,nrow = 2)  
  q <- q + scale_color_manual(values = c("Cluster1"="red","Cluster2"="green","Cluster3"="blue","Cluster4"="cyan"))
  q 
  
  
  
  
}

if(F){
  rm(list = ls())
  library(EnrichedHeatmap)
  library(circlize)
  library(ComplexHeatmap)
  
  load("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/RawCounts_Signal_genes.RData")
  cf <- cf[cf$species=="Scer",]
  cf <- cf[,c(6,65:68)]
  cf$RNAPII <- rowSums(cf[,2:5])
  cf <- cf[,c(1,6)]
  load(file="A:/work/WinstonLab/Natalia_Reim/Analysis_Feb18/DE_Genes_Filtered.RData")
  load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenes.RData")
  genes <- subset(genes,genes$verification=="Verified" & genes$species=="Scer" & genes$chr!="chrM")
  genes <- merge(genes,cf,by="tracking_id",all.x=T)
  genes <- merge(genes,de.genes[,c("tracking_id","log2FoldChange")],by="tracking_id",all.x=T)
  genes$RNAPII <- ifelse(is.na(genes$RNAPII),0,genes$RNAPII)
  genes$log2FoldChange <- ifelse(is.na(genes$log2FoldChange),0,genes$log2FoldChange)
  genes$fpkm <- round(genes$RNAPII*1000/genes$width,2)
  load("A:/work/yeast_Annotations/refs.RData")
  genes <- merge(genes,unique(refs[,c("tracking_id","gene","link")]),by="tracking_id",all.x=T)
  
  g <- as(genes,"GRanges")
  d <- as.data.frame(findOverlaps(g,ignore.strand=T))
  d <- d[d$queryHits!=d$subjectHits,] 
  d <- unique(d$queryHits,d$subjectHits)
  
  myheatplotfun <- function(path,normF=c(1,1)){
    load(path)
    a1 <- l[["Depleted"]]
    a2 <- l[["Non-depleted"]]
    a1 <- a1*normF[1]
    a2 <- a2*normF[2]
    get_divsion <- function(m,n){
      m <- m*100+1
      n <- n*100+1
      m <- (m-n)/(m+n)
      return(m)
    }
    
    #a <- log2((a1+1)/(a2+1))
    #a[!is.finite(a)] <-NA
    #range(a,na.rm=T)
    a <- get_divsion(a1,a2)
    a[!is.finite(a)] <-0
    a <- a[-d,]
    
    dim(a) = dim(a)
    attr(a, "upstream_index") = 1:50
    attr(a, "target_index") = 51:250
    attr(a, "downstream_index") = 251:300
    attr(a, "extend") = c(0.5,0.5)  # it must be a vector of length two
    class(a) = c("normalizedMatrix", "matrix")
    attr(a, "signal_name") = "Spt6"
    attr(a, "axis_name") = c(-0.5,"TSS","CPS",0.5)
    attr(a,"axis_name_rot") = 90
    attr(a, "target_name") = "TSS"
    attr(a,"target_is_single_point") = TRUE
    return(a)
  }
  
  res <- read.delim("SpikeinNormalization_factors_revised.txt",header = T)
  
  spt6 <- myheatplotfun(path="A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_ScHeat_Spt6.RData",res[res$Factor=="Spt6",]$normF)
  spn1 <- myheatplotfun("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_ScHeat_Spn1.RData",res[res$Factor=="Spn1",]$normF)
  k36me3 <- myheatplotfun("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_ScHeat_H3K36me3.RData",res[res$Factor=="H3K36me3",]$normF)
  k36me2 <- myheatplotfun("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_ScHeat_H3K36me2.RData",res[res$Factor=="H3K36me2",]$normF)
  k4me3 <- myheatplotfun("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_ScHeat_H3K4me3.RData",res[res$Factor=="H3K4me3",]$normF)
  h3 <- myheatplotfun("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_ScHeat_H3.RData",res[res$Factor=="H3",]$normF)
  pol2 <- myheatplotfun("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_ScHeat_RNAPII.RData",res[res$Factor=="RNAPII",]$normF)
  set2 <- myheatplotfun("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_ScHeat_Set2.RData",res[res$Factor=="Set2",]$normF)
  s2p <- myheatplotfun("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_ScHeat_S2P.RData",res[res$Factor=="S2P",]$normF)
  s5p <- myheatplotfun("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_ScHeat_S5P.RData",res[res$Factor=="S5P",]$normF)
  
  #EnrichedHeatmap(spt6,col = col_fun_spt6,axis_name = c("","TSS","CPS",""), name = "Spt6",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), column_title = "Spt6")
  set.seed(123)
  col_fun_spt6 = colorRamp2(c(min(spt6),0, max(spt6)), c("blue4","white","red4"))
  col_fun_pol2 = colorRamp2(c(min(pol2),0, max(pol2)), c("blue4","white","red4"))
  col_fun_s2p = colorRamp2(c(min(s2p),0, max(s2p)), c("blue4","white","red4"))
  col_fun_s5p = colorRamp2(c(min(s5p),0, max(s5p)), c("blue4","white","red4"))
  col_fun_spn1 = colorRamp2(c(min(spn1), 0,max(spn1)), c("blue4","white","red4"))
  col_fun_k36me3 = colorRamp2(c(min(k36me3),0, max(k36me3)),  c("blue4","white","red4"))
  col_fun_k36me2 = colorRamp2(c(min(k36me2),0, max(k36me2)),  c("blue4","white","red4"))
  col_fun_k4me3 = colorRamp2(c(min(k4me3),0, max(k4me3)),  c("blue4","white","red4"))
  col_fun_h3 = colorRamp2(c(min(h3),0, max(h3)),  c("blue4","white","red4"))
  col_fun_set2 = colorRamp2(c(min(set2),0, max(set2)),  c("blue4","white","red4"))
  
  
  ## ordering by Spt6 change
  load("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_Spt6.RData")
  z <- data.frame(depl=apply(l[[1]][,51:1523],1,sum)+apply(l[[2]][,51:1523],1,sum),
                  nondepl=apply(l[[3]][,51:1523],1,sum)+apply(l[[4]][,51:1523],1,sum))
  z$depl <- z$depl*res[res$Factor=="Spt6",]$normF[1]
  z$nondepl <- z$nondepl*res[res$Factor=="Spt6",]$normF[2]
  z$o <- (z$depl-z$nondepl)/(z$depl+z$nondepl)
  z <- z[-d,]
  
  load("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_H3K36me3.RData")
  z1 <- data.frame(depl=apply(l[[1]][,51:1523],1,sum)+apply(l[[2]][,51:1523],1,sum),
                  nondepl=apply(l[[3]][,51:1523],1,sum)+apply(l[[4]][,51:1523],1,sum))
  z1$depl <- z1$depl*res[res$Factor=="H3K36me3",]$normF[1]
  z1$nondepl <- z1$nondepl*res[res$Factor=="H3K36me3",]$normF[2]
  z1$o <- (z1$depl-z1$nondepl)/(z1$depl+z1$nondepl)
  z1 <- z1[-d,]
  
  load("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_Spn1.RData")
  z2 <- data.frame(depl=apply(l[[1]][,51:1523],1,sum)+apply(l[[2]][,51:1523],1,sum),
                   nondepl=apply(l[[3]][,51:1523],1,sum)+apply(l[[4]][,51:1523],1,sum))
  z2$depl <- z2$depl*res[res$Factor=="Spn1",]$normF[1]
  z2$nondepl <- z2$nondepl*res[res$Factor=="Spn1",]$normF[2]
  z2$o <- (z2$depl-z2$nondepl)/(z2$depl+z2$nondepl)
  z2 <- z2[-d,]
  rm(l)
  #EnrichedHeatmap(spt6,col = col_fun_spt6, name = "Spt6",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), column_title = "Spt6")
  
  myplotfun1 <- function(row_order){
    ht_list =  EnrichedHeatmap(spt6,row_order = row_order,axis_name = c("","TSS","CPS",""), col = col_fun_spt6, name = "Spt6",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                               column_title = "Spt6")+
      EnrichedHeatmap(pol2, col = col_fun_pol2,axis_name = c("","TSS","CPS",""),row_order = row_order, name = "RNAPII",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                      column_title = "RNAPII") +
      EnrichedHeatmap(s2p, col = col_fun_s2p,axis_name = c("","TSS","CPS",""),row_order = row_order, name = "Ser2-P",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                      column_title = "Ser2-P") +
      EnrichedHeatmap(s5p, col = col_fun_s5p,axis_name = c("","TSS","CPS",""),row_order = row_order, name = "Ser5-P",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                      column_title = "Ser5-P") +
      EnrichedHeatmap(k36me3, col = col_fun_k36me3,row_order = row_order,axis_name = c("","TSS","CPS",""), name = "H3K36me3",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                      column_title = "H3K36me3") +
      EnrichedHeatmap(k36me2, col = col_fun_k36me2,axis_name = c("","TSS","CPS",""),row_order = row_order, name = "H3K36me2",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                      column_title = "H3K36me2") +
      EnrichedHeatmap(k4me3, col = col_fun_k4me3,axis_name = c("","TSS","CPS",""),row_order = row_order, name = "H3K4me3",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                      column_title = "H3K4me3") +
      EnrichedHeatmap(set2, col = col_fun_set2,axis_name = c("","TSS","CPS",""),row_order = row_order, name = "Set2",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                      column_title = "Set2") +
      EnrichedHeatmap(h3, col = col_fun_h3,axis_name = c("","TSS","CPS",""),row_order = row_order, name = "H3",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                      column_title = "H3") +
      EnrichedHeatmap(spn1, col = col_fun_spn1,axis_name = c("","TSS","CPS",""),row_order = row_order, name = "Spn1",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4))), 
                      column_title = "Spn1") +
      Heatmap(log2(genes[-d,]$fpkm+1),row_order = row_order, col = c("white","lightpink","pink", "red4"), name = "log2(fpkm+1)", 
              show_row_names = FALSE, width = unit(4, "mm"))
    return(ht_list)
  }
  
  png("Spt6 sorted_signal_DeplNonDepl_ScHeat.png",width = 1500,height = 800)
  ht_list = myplotfun1(order(z$o,decreasing = T))
  draw(ht_list,heatmap_legend_side = "bottom", gap = unit(2, "mm"))
  dev.off()
  
  png("H3K36me3 sorted_signal_DeplNonDepl_ScHeat.png",width = 1500,height = 800)
  ht_list = myplotfun1(order(z1$o,decreasing = T))
  draw(ht_list,heatmap_legend_side = "bottom", gap = unit(2, "mm"))
  dev.off()
  
  png("Spn1 sorted_signal_DeplNonDepl_ScHeat.png",width = 1500,height = 800)
  ht_list = myplotfun1(order(z2$o,decreasing = T))
  draw(ht_list,heatmap_legend_side = "bottom", gap = unit(2, "mm"))
  dev.off()
  
  
  ## clustering
  myplotfun2 <- function(row_order,zz){
    partition = paste0("cl", kmeans(zz[,51:250], centers = 4)$cluster)
    #lgd = Legend(at = c("cluster1", "cluster2", "cluster3","cluster4"), title = "Clusters", type = "lines", legend_gp = gpar(col = 2:5))
    #o <- data.frame(ord=1:length(row_order),cl=partition)
    #o1 <- data.frame(ord=row_order,cl1=paste0("clt", kmeans(zz[row_order,51:250], centers = 4)$cluster))
    #o <- merge(o,o1,by="ord")
    #o1 <- as.matrix(table(o$cl,o$cl1))
    o <- data.frame(ord=1:length(row_order),cl=partition)
    o <- o[row_order,]
    ord <- c()
    v <- o$cl
    for(i in 1:4){
      ord <- c(ord,as.character(v[1]))
      v <- subset(v,v!= v[1])
    }
    partition = gsub(ord[1],"cluster1",partition)
    partition = gsub(ord[2],"cluster2",partition)
    partition = gsub(ord[3],"cluster3",partition)
    partition = gsub(ord[4],"cluster4",partition)
    partition <- factor(partition,levels=c("cluster1", "cluster2", "cluster3","cluster4"))
    
    ht_list = Heatmap(partition, col = structure(2:5, names = paste0("cluster", 1:4)), name = "",
                      show_row_names = FALSE, width = unit(1, "mm")) +
      EnrichedHeatmap(spt6, col = col_fun_spt6,axis_name = c("","TSS","CPS",""), name = "Spt6",
                      top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                      column_title = "Spt6") + 
      EnrichedHeatmap(pol2, col = col_fun_pol2,axis_name = c("","TSS","CPS",""), name = "RNAPII",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                      column_title = "RNAPII") +
      EnrichedHeatmap(s2p, col = col_fun_s2p,axis_name = c("","TSS","CPS",""), name = "Ser2-P",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                      column_title = "Ser2-P") +
      EnrichedHeatmap(s5p, col = col_fun_s5p,axis_name = c("","TSS","CPS",""), name = "Ser5-P",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                      column_title = "Ser5-P") +
      EnrichedHeatmap(k36me3, col = col_fun_k36me3,axis_name = c("","TSS","CPS",""), name = "H3K36me3",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                      column_title = "H3K36me3") +
      EnrichedHeatmap(k36me2, col = col_fun_k36me2,axis_name = c("","TSS","CPS",""), name = "H3K36me2",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                      column_title = "H3K36me2") +
      EnrichedHeatmap(k4me3, col = col_fun_k4me3,axis_name = c("","TSS","CPS",""), name = "H3K4me3",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                      column_title = "H3K4me3") +
      EnrichedHeatmap(set2, col = col_fun_set2,axis_name = c("","TSS","CPS",""),name = "Set2",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                      column_title = "Set2") +
      EnrichedHeatmap(h3, col = col_fun_h3,axis_name = c("","TSS","CPS",""),name = "H3",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                      column_title = "H3") +
      EnrichedHeatmap(spn1, col = col_fun_spn1,axis_name = c("","TSS","CPS",""), name = "Spn1",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5))), 
                      column_title = "Spn1") +
      Heatmap(log2(genes[-d,]$fpkm+1), col = c("white","lightpink","pink", "red4"),name = "log2(fpkm+1)", 
              show_row_names = FALSE, width = unit(4, "mm"))
    
    
    draw(ht_list,  split = partition,row_order=row_order,
         heatmap_legend_side = "bottom", gap = unit(2, "mm"))
    
    load("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_H3K36me3.RData")
    m <- l[[1]]
    
    cf1 <- data.frame(tracking_id=as.character(rownames(m)))
    cf1$tracking_id <- as.character(cf1$tracking_id)
    cf1 <- as.data.frame(cf1[-d,])
    names(cf1) <- "tracking_id"
    cf1$cluster <- partition
    cf1 <- merge(cf1,genes,by="tracking_id",all.x=T)
    return(cf1)
    
  }
  
  png("Spt6 sorted_signal_DeplNonDepl_ScHeat_4clusters.png",width = 1500,height = 800)
  cf1 <- myplotfun2(row_order = order(z$o,decreasing = T),spt6)
  dev.off()
  
  png("H3K36me3 sorted_signal_DeplNonDepl_ScHeat_4clusters.png",width = 1500,height = 800)
  cf2 <- myplotfun2(order(z1$o,decreasing = T),k36me3)
  dev.off()
  
  png("Spn1 sorted_signal_DeplNonDepl_ScHeat_4clusters.png",width = 1500,height = 800)
  cf3 <- myplotfun2(order(z2$o,decreasing = T),spn1)
  dev.off()
  
  
  ##How many clusters?
  library(factoextra)
  library(NbClust)
  # Elbow method
  pdf("Optimal Number of Clusters for Spt6_H3K36me3_Spn1.pdf",5,5)
  mydata <- as.matrix(spt6)
  fviz_nbclust(mydata[,51:250], kmeans, method = "wss") +
    geom_vline(xintercept = 4, linetype = 2)+
    labs(subtitle = "Spt6")
  
  mydata <- as.matrix(k36me3)
  fviz_nbclust(mydata[,51:250], kmeans, method = "wss") +
    geom_vline(xintercept = 4, linetype = 2)+
    labs(subtitle = "H3K36me3")
  
  mydata <- as.matrix(spn1)
  fviz_nbclust(mydata[,51:250], kmeans, method = "wss") +
    geom_vline(xintercept = 4, linetype = 2)+
    labs(subtitle = "Spn1")
  
  dev.off()
  
  write.table(cf1,file = "Spt6_4Clusters.txt",quote = F,sep = "\t",row.names = F)
  write.table(cf2,file = "H3K36me3_4Clusters.txt",quote = F,sep = "\t",row.names = F)
  write.table(cf3,file = "Spn1_4Clusters.txt",quote = F,sep = "\t",row.names = F)
  
  ## cluster4 of both H3K36me3 and Spt6 are the genes with largest change
  
  table(cf1$cluster)
  table(cf2$cluster)
  pdf("Spt6_H3K36me3 strongest loss.pdf",5,5)
  p <- FP22.Venn2(as.character(cf1[cf1$cluster=="cluster4",]$tracking_id),
             as.character(cf1[cf2$cluster=="cluster4",]$tracking_id),names = c("Spt6","H3K36me3"))
  dev.off()
  cf4 <- p[p$condition=="Spt6 H3K36me3",]
  cf4 <- merge(cf4,genes,by="tracking_id",all.x=T)
  
  write.table(cf4,file = "Strongest loss of Spt6 and H3K36me3.txt",quote = F,sep = "\t",row.names = F)
  
}

## relative signals
if(F){
  setwd("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr")
  rm(list = ls())
  library(EnrichedHeatmap)
  library(circlize)
  library(ComplexHeatmap)
  
  
  load("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/RawCounts_Signal_genes.RData")
  cf <- cf[cf$species=="Scer",]
  cf <- cf[,c(6,65:68)]
  cf$RNAPII <- rowSums(cf[,2:5])
  cf <- cf[,c(1,6)]
  load(file="A:/work/WinstonLab/Natalia_Reim/Analysis_Feb18/DE_Genes_Filtered.RData")
  load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenes.RData")
  genes <- subset(genes,genes$verification=="Verified" & genes$species=="Scer" & genes$chr!="chrM")
  genes <- merge(genes,cf,by="tracking_id",all.x=T)
  genes <- merge(genes,de.genes[,c("tracking_id","log2FoldChange")],by="tracking_id",all.x=T)
  genes$RNAPII <- ifelse(is.na(genes$RNAPII),0,genes$RNAPII)
  genes$log2FoldChange <- ifelse(is.na(genes$log2FoldChange),0,genes$log2FoldChange)
  genes$fpkm <- round(genes$RNAPII*1000/genes$width,2)
  load("A:/work/yeast_Annotations/refs.RData")
  genes <- merge(genes,unique(refs[,c("tracking_id","gene","link")]),by="tracking_id",all.x=T)
  
  g <- as(genes,"GRanges")
  d <- as.data.frame(findOverlaps(g,ignore.strand=T))
  d <- d[d$queryHits!=d$subjectHits,] 
  d <- unique(d$queryHits,d$subjectHits)
  
  
  myheatplotfun <- function(path1,path2,normF1=c(1,1),normF2=c(1,1)){
    load(path1)  ## Spt6    
    l1 <- l
    rm(l)
    
    load(path2)  ## RNAPII
    l2 <- l
    rm(l)
    
    l1[["Depleted"]] <- l1[["Depleted"]]*normF1[1]
    l1[["Non-depleted"]] <- l1[["Non-depleted"]]*normF1[2]
    l2[["Depleted"]] <- l2[["Depleted"]]*normF2[1]
    l2[["Non-depleted"]] <- l2[["Non-depleted"]]*normF2[2]
    
    get_divsion <- function(m,n){
      m <- m*100+1
      n <- n*100+1
      m <- (m-n)/(m+n)
      return(m)
    }

    a <- (l1[["Depleted"]]*100+1)/(l2[["Depleted"]]*100+1)
    b <- (l1[["Non-depleted"]]*100+1)/(l2[["Non-depleted"]]*100+1)
    
    a <- get_divsion(a,b)
    a[!is.finite(a)] <-0
    
    dim(a) = dim(a)
    rownames(a) <- rownames(l1[[1]])
    attr(a, "upstream_index") = 1:50
    attr(a, "target_index") = 51:250
    attr(a, "downstream_index") = 251:300
    attr(a, "extend") = c(0.5,0.5)  # it must be a vector of length two
    class(a) = c("normalizedMatrix", "matrix")
    attr(a, "signal_name") = "Spt6"
    attr(a, "axis_name") = c(-0.5,"TSS","CPS",0.5)
    attr(a,"axis_name_rot") = 90
    attr(a, "target_name") = "TSS"
    attr(a,"target_is_single_point") = TRUE
    
    a <- a[-d,]
    return(a)
  }
  
  res <- read.delim("SpikeinNormalization_factors_revised.txt",header = T)
  
  spt6 <- myheatplotfun(path1="A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_ScHeat_Spt6.RData",
                        path2="A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_ScHeat_RNAPII.RData",
                        res[res$Factor=="Spt6",]$normF,res[res$Factor=="RNAPII",]$normF)
  s2p <- myheatplotfun(path1="A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_ScHeat_S2P.RData",
                        path2="A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_ScHeat_RNAPII.RData",
                        res[res$Factor=="S2P",]$normF,res[res$Factor=="RNAPII",]$normF)
  s5p <- myheatplotfun(path1="A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_ScHeat_S5P.RData",
                        path2="A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_ScHeat_RNAPII.RData",
                        res[res$Factor=="S5P",]$normF,res[res$Factor=="RNAPII",]$normF)
  k36me3 <- myheatplotfun(path1="A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_ScHeat_H3K36me3.RData",
                       path2="A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_ScHeat_H3.RData",
                       res[res$Factor=="H3K36me3",]$normF,res[res$Factor=="H3",]$normF)
  k36me2 <- myheatplotfun(path1="A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_ScHeat_H3K36me2.RData",
                       path2="A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_ScHeat_H3.RData",
                       res[res$Factor=="H3K36me2",]$normF,res[res$Factor=="H3",]$normF)
  k4me3 <- myheatplotfun(path1="A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_ScHeat_H3K4me3.RData",
                       path2="A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_ScHeat_H3.RData",
                       res[res$Factor=="H3K4me3",]$normF,res[res$Factor=="H3",]$normF)
  spn1 <- myheatplotfun(path1="A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_ScHeat_Spn1.RData",
                        path2="A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_ScHeat_RNAPII.RData",
                        res[res$Factor=="Spn1",]$normF,res[res$Factor=="RNAPII",]$normF)
  
  col_fun_spt6 = colorRamp2(c(min(spt6),0, max(spt6)), c("blue4","white","red4"))
  col_fun_s2p = colorRamp2(c(min(s2p),0, max(s2p)), c("blue4","white","red4"))
  col_fun_s5p = colorRamp2(c(min(s5p),0, max(s5p)), c("blue4","white","red4"))
  col_fun_k36me3 = colorRamp2(c(min(k36me3),0, max(k36me3)),  c("blue4","white","red4"))
  col_fun_k36me2 = colorRamp2(c(min(k36me2),0, max(k36me2)),  c("blue4","white","red4"))
  col_fun_k4me3 = colorRamp2(c(min(k4me3),0, max(k4me3)),  c("blue4","white","red4"))
  col_fun_spn1 = colorRamp2(c(min(spn1), 0,max(spn1)), c("blue4","white","red4"))
  
  ## ordering by Spt6 relative change
  load("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_Spt6.RData")
  z <- data.frame(depl=apply(l[[1]][,51:1523],1,sum)+apply(l[[2]][,51:1523],1,sum),
                  nondepl=apply(l[[3]][,51:1523],1,sum)+apply(l[[4]][,51:1523],1,sum))
  load("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_RNAPII.RData")
  z <- cbind(z, cdepl=apply(l[[1]][,51:1523],1,sum)+apply(l[[2]][,51:1523],1,sum),
                  cnondepl=apply(l[[3]][,51:1523],1,sum)+apply(l[[4]][,51:1523],1,sum))
  z$depl <- z$depl*res[res$Factor=="Spt6",]$normF[1]
  z$nondepl <- z$nondepl*res[res$Factor=="Spt6",]$normF[2]
  z$cdepl <- z$depl*res[res$Factor=="RNAPII",]$normF[1]
  z$cnondepl <- z$cnondepl*res[res$Factor=="RNAPII",]$normF[2]
  z$depl<- z$depl/z$cdepl
  z$nondepl<-z$nondepl/z$cnondepl
  z$o <- (z$depl-z$nondepl)/(z$depl+z$nondepl)
  z$o <- ifelse(!is.finite(z$o),0,z$o)
  z <- z[-d,]
  
  
  load("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_H3K36me3.RData")
  z1 <- data.frame(depl=apply(l[[1]][,51:1523],1,sum)+apply(l[[2]][,51:1523],1,sum),
                   nondepl=apply(l[[3]][,51:1523],1,sum)+apply(l[[4]][,51:1523],1,sum))
  load("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_H3.RData")
  z1 <- cbind(z1, cdepl=apply(l[[1]][,51:1523],1,sum)+apply(l[[2]][,51:1523],1,sum),
             cnondepl=apply(l[[3]][,51:1523],1,sum)+apply(l[[4]][,51:1523],1,sum))
  z1$depl <- z1$depl*res[res$Factor=="H3K36me3",]$normF[1]
  z1$nondepl <- z1$nondepl*res[res$Factor=="H3K36me3",]$normF[2]
  z1$cdepl <- z1$cdepl*res[res$Factor=="H3",]$normF[1]
  z1$cnondepl <- z1$cnondepl*res[res$Factor=="H3",]$normF[2]
  z1$depl<- z1$depl/z1$cdepl
  z1$nondepl<-z1$nondepl/z1$cnondepl
  z1$o <- (z1$depl-z1$nondepl)/(z1$depl+z1$nondepl)
  z1$o <- ifelse(!is.finite(z1$o),0,z1$o)
  z1 <- z1[-d,]
  
  load("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_Spn1.RData")
  z2 <- data.frame(depl=apply(l[[1]][,51:1523],1,sum)+apply(l[[2]][,51:1523],1,sum),
                   nondepl=apply(l[[3]][,51:1523],1,sum)+apply(l[[4]][,51:1523],1,sum))
  load("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_RNAPII.RData")
  z2 <- cbind(z2, cdepl=apply(l[[1]][,51:1523],1,sum)+apply(l[[2]][,51:1523],1,sum),
             cnondepl=apply(l[[3]][,51:1523],1,sum)+apply(l[[4]][,51:1523],1,sum))
  z2$depl <- z2$depl*res[res$Factor=="Spt6",]$normF[1]
  z2$nondepl <- z2$nondepl*res[res$Factor=="Spt6",]$normF[2]
  z2$cdepl <- z2$cdepl*res[res$Factor=="RNAPII",]$normF[1]
  z2$cnondepl <- z2$cnondepl*res[res$Factor=="RNAPII",]$normF[2]
  z2$depl<- z2$depl/z2$cdepl
  z2$nondepl<-z2$nondepl/z2$cnondepl
  z2$o <- (z2$depl-z2$nondepl)/(z2$depl+z2$nondepl)
  z2$o <- ifelse(!is.finite(z2$o),0,z2$o)
  z2 <- z2[-d,]
  rm(l)

  myplotfun1 <- function(row_order){
    ht_list =  EnrichedHeatmap(spt6,row_order = row_order,axis_name = c("","TSS","CPS",""), col = col_fun_spt6, name = "Spt6",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                               column_title = "Spt6")+
      EnrichedHeatmap(s2p, col = col_fun_s2p,axis_name = c("","TSS","CPS",""),row_order = row_order, name = "Ser2-P",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                      column_title = "Ser2-P") +
      EnrichedHeatmap(s5p, col = col_fun_s5p,axis_name = c("","TSS","CPS",""),row_order = row_order, name = "Ser5-P",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                      column_title = "Ser5-P") +
      EnrichedHeatmap(k36me3, col = col_fun_k36me3,row_order = row_order,axis_name = c("","TSS","CPS",""), name = "H3K36me3",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                      column_title = "H3K36me3") +
      EnrichedHeatmap(k36me2, col = col_fun_k36me2,axis_name = c("","TSS","CPS",""),row_order = row_order, name = "H3K36me2",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                      column_title = "H3K36me2") +
      EnrichedHeatmap(k4me3, col = col_fun_k4me3,axis_name = c("","TSS","CPS",""),row_order = row_order, name = "H3K4me3",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4),yaxis_facing = "left")), 
                      column_title = "H3K4me3") +
      EnrichedHeatmap(spn1, col = col_fun_spn1,axis_name = c("","TSS","CPS",""),row_order = row_order, name = "Spn1",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4))), 
                      column_title = "Spn1") +
      Heatmap(log2(genes[-d,]$fpkm+1),row_order = row_order, col = c("white","lightpink","pink", "red4"), name = "log2(fpkm+1)", 
              show_row_names = FALSE, width = unit(4, "mm"))
    return(ht_list)
  }
  png("Spt6 sorted_signal_DeplNonDepl_ScHeat_relSignal.png",width = 1500,height = 800)
  ht_list=myplotfun1(row_order = order(z$o,decreasing = T))
  draw(ht_list,heatmap_legend_side = "bottom", gap = unit(2, "mm"))
  dev.off()
  png("H3K36me3 sorted_signal_DeplNonDepl_ScHeat_relSignal.png",width = 1500,height = 800)
  ht_list=myplotfun1(order(z1$o,decreasing = T))
  draw(ht_list,heatmap_legend_side = "bottom", gap = unit(2, "mm"))
  dev.off()
  png("Spn1 sorted_signal_DeplNonDepl_ScHeat_relSignal.png",width = 1500,height = 800)
  ht_list=myplotfun1(order(z2$o,decreasing = T))
  draw(ht_list,heatmap_legend_side = "bottom", gap = unit(2, "mm"))
  dev.off()
  
  pdf("Optimal Number of Clusters for Spt6_H3K36me3_Spn1_relsignal.pdf",5,5)
  mydata <- as.matrix(spt6)
  fviz_nbclust(mydata[,51:250], kmeans, method = "wss") +
    geom_vline(xintercept = 4, linetype = 2)+
    labs(subtitle = "Spt6")
  
  mydata <- as.matrix(k36me3)
  fviz_nbclust(mydata[,51:250], kmeans, method = "wss") +
    geom_vline(xintercept = 4, linetype = 2)+
    labs(subtitle = "H3K36me3")
  
  mydata <- as.matrix(spn1)
  fviz_nbclust(mydata[,51:250], kmeans, method = "wss") +
    geom_vline(xintercept = 4, linetype = 2)+
    labs(subtitle = "Spn1")
  
  dev.off()
  
  
  myplotfun2 <- function(row_order,zz){
    if(T){
    partition = paste0("cl", kmeans(zz[,51:250], centers = 4)$cluster)
    #lgd = Legend(at = c("cluster1", "cluster2", "cluster3","cluster4"), title = "Clusters", type = "lines", legend_gp = gpar(col = 2:5))
    #o <- data.frame(ord=1:length(row_order),cl=partition)
    #o1 <- data.frame(ord=row_order,cl1=paste0("clt", kmeans(zz[row_order,51:250], centers = 4)$cluster))
    #o <- merge(o,o1,by="ord")
    #o1 <- as.matrix(table(o$cl,o$cl1))
    o <- data.frame(ord=1:length(row_order),cl=partition)
    o <- o[row_order,]
    ord <- c()
    v <- o$cl
    for(i in 1:4){
      ord <- c(ord,as.character(v[1]))
      v <- subset(v,v!= v[1])
    }
    partition = gsub(ord[1],"cluster1",partition)
    partition = gsub(ord[2],"cluster2",partition)
    partition = gsub(ord[3],"cluster3",partition)
    partition = gsub(ord[4],"cluster4",partition)
    partition <- factor(partition,levels=c("cluster1", "cluster2", "cluster3","cluster4"))
    }
   if(F){    
    partition = paste0("cl", kmeans(zz[,51:250], centers = 4)$cluster)
    #lgd = Legend(at = c("cluster1", "cluster2", "cluster3","cluster4"), title = "Clusters",type = "lines", legend_gp = gpar(col = 2:5))
    o <- data.frame(ord=1:length(row_order),cl=partition)
    o <- o[row_order,]
    ord <- c()
    v <- as.character(o$cl)
    for(i in 1:4){
      vx = names(sort(table(as.character(v[1:10])),decreasing = T))[1]
      ord <- c(ord,vx)
      v <- subset(v,v!= vx)
    }
    partition = gsub(ord[1],"cluster1",partition)
    partition = gsub(ord[2],"cluster2",partition)
    partition = gsub(ord[3],"cluster3",partition)
    partition = gsub(ord[4],"cluster4",partition)
    #partition <- factor(partition,levels=c("cluster1", "cluster2", "cluster3","cluster4"))
   }
    cf1 <- data.frame(num=1:length(rownames(zz)),tracking_id=rownames(zz),cluster = partition)
    cf1 <- cf1[row_order,]
    cf1 <- split(cf1,cf1$cluster)
    
    cf1 <- lapply(cf1,function(x){
      x <- as.data.frame(x)
      ro <- subset(row_order,row_order%in%x$num)
      x <- x[match(ro,x$num),]
      x
    })
    cf <- c()
    for(i in 1:length(cf1)){
        cf <- rbind(cf,cf1[[i]])    
    }
    
    row_order = cf$num    
    cf$tracking_id<- as.character(cf$tracking_id)
    cf <- merge(cf,genes,by="tracking_id",all.x=T)
    cf <- cf[match(row_order,cf$num),]
    
    ht_list = Heatmap(partition, col = structure(2:5, names = paste0("cluster", 1:4)), name = "",
                      show_row_names = FALSE, width = unit(1, "mm")) +
      EnrichedHeatmap(spt6, col = col_fun_spt6,axis_name = c("","TSS","CPS",""), name = "Spt6",
                      top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                      column_title = "Spt6") + 
      EnrichedHeatmap(s2p, col = col_fun_s2p,axis_name = c("","TSS","CPS",""), name = "Ser2-P",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                      column_title = "Ser2-P") +
      EnrichedHeatmap(s5p, col = col_fun_s5p,axis_name = c("","TSS","CPS",""), name = "Ser5-P",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                      column_title = "Ser5-P") +
      EnrichedHeatmap(k36me3, col = col_fun_k36me3,axis_name = c("","TSS","CPS",""), name = "H3K36me3",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                      column_title = "H3K36me3") +
      EnrichedHeatmap(k36me2, col = col_fun_k36me2,axis_name = c("","TSS","CPS",""), name = "H3K36me2",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                      column_title = "H3K36me2") +
      EnrichedHeatmap(k4me3, col = col_fun_k4me3,axis_name = c("","TSS","CPS",""), name = "H3K4me3",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5),yaxis_facing = "left")), 
                      column_title = "H3K4me3") +
      EnrichedHeatmap(spn1, col = col_fun_spn1,axis_name = c("","TSS","CPS",""), name = "Spn1",top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:5))), 
                      column_title = "Spn1") +
      Heatmap(log2(genes[-d,]$fpkm+1), col = c("white","lightpink","pink", "red4"),row_order = row_order, name = "log2(fpkm+1)", 
              show_row_names = FALSE, width = unit(4, "mm"))
    draw(ht_list, split = partition, row_order=row_order
         , heatmap_legend_side = "bottom", gap = unit(2, "mm"))
    
    return(cf)
    
  }
  png("Spt6 sorted_signal_DeplNonDepl_ScHeat_4clusters_relSignal.png",width = 1500,height = 800)
  cf1 <- myplotfun2(row_order = order(z$o,decreasing = T),zz = spt6)
  dev.off()
  
  png("H3K36me3 sorted_signal_DeplNonDepl_ScHeat_4clusters_relSignal.png",width = 1500,height = 800)
  cf2 <- myplotfun2(order(z1$o,decreasing = T),k36me3)
  dev.off()
  png("Spn1 sorted_signal_DeplNonDepl_ScHeat_4clusters_relSignal.png",width = 1500,height = 800)
  cf3 <- myplotfun2(order(z2$o,decreasing = T),spn1)
  dev.off()
  
  
  write.table(cf1,file = "Spt6_4Clusters_RelSignal.txt",quote = F,sep = "\t",row.names = F)
  write.table(cf2,file = "H3K36me3_4Clusters_RelSignal.txt",quote = F,sep = "\t",row.names = F)
  write.table(cf3,file = "Spn1_4Clusters_RelSignal.txt",quote = F,sep = "\t",row.names = F)
  
  
  table(cf1$cluster)
  table(cf2$cluster)
  table(cf3$cluster)
  pdf("Spt6_H3K36me3 strongest loss_relsignal.pdf",5,5)
  p <- FP22.Venn2(as.character(cf1[cf1$cluster=="cluster4",]$tracking_id),
                  as.character(cf1[cf2$cluster=="cluster4",]$tracking_id),names = c("Spt6","H3K36me3"))
  dev.off()
  cf4 <- p[p$condition=="Spt6 H3K36me3",]
  cf4 <- merge(cf4,genes,by="tracking_id",all.x=T)
  
  write.table(cf4,file = "Strongest loss of Spt6 and H3K36me3_relsingal.txt",quote = F,sep = "\t",row.names = F)
  
  

}

#######################################################################################