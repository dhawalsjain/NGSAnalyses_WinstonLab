rm(list=ls())
library(ggplot2)
library(gridExtra)
library(grid)
library(GenomicRanges)
library(gtable)
library(cowplot)
library(reshape)
cols = c("Spn1-depleted" = "#BB5566","non-depleted" = "#4477AA",
         "Spn1-depleted (down)" = "#f2d5da","non-depleted (down)"="#abc9e8",
         "depl." = "#BB5566","non-depl." = "#4477AA",
         "SPN1-AID_DMSO-1"="#4477AA", "SPN1-AID_DMSO-2"="#4477AA", 
         "SPN1-AID_IAA-1"="#BB5566", "SPN1-AID_IAA-2"="#BB5566")
OUTDIR="A:/work/WinstonLab/Natalia_Reim/Paper_Thesis_Figs/"
setwd("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr")

######-------------- plot functions
vector.resizing <- function(x,final.len){
  y <- vector()
  len <- length(x)
  y <-spline(1:len,x,n=final.len)$y
  return(y)
}

mytotalgenesignal <- function(l,genes,max.gene.len=4000, bw=10,Margin=500){
  l[[1]] <- subset(l[[1]],rownames(l[[1]])%in%genes$tracking_id)
  l[[2]] <- subset(l[[2]],rownames(l[[2]])%in%genes$tracking_id)
  l[[3]] <- subset(l[[3]],rownames(l[[3]])%in%genes$tracking_id)
  l[[4]] <- subset(l[[4]],rownames(l[[4]])%in%genes$tracking_id)
  
  max.gene.len = round(max.gene.len/bw)
  Margin = round(Margin/bw)
  mcol= max.gene.len+2*Margin-1
  
  depl1 <- matrix(NA,ncol=1,nrow = nrow(genes))
  depl2 <- matrix(NA,ncol=1,nrow = nrow(genes))
  nondepl1 <- matrix(NA,ncol=1,nrow = nrow(genes))
  nondepl2 <- matrix(NA,ncol=1,nrow = nrow(genes))
  ncol=dim(l[[1]])[2] 
  
  for(i in 1:nrow(genes)){
    j = ceiling(genes[i,]$width/bw)-1
    v <- l[[1]][i,(Margin+1):(Margin+j)]
    depl1[i,1] <- sum(v,na.rm = T)
    v <- l[[2]][i,(Margin+1):(Margin+j)]
    depl2[i,1] <- sum(v,na.rm = T)
    v <- l[[3]][i,(Margin+1):(Margin+j)]
    nondepl1[i,1] <- sum(v,na.rm = T)
    v <- l[[4]][i,(Margin+1):(Margin+j)]
    nondepl2[i,1] <- sum(v,na.rm=T)
  }
  depl1 <- round((depl1+depl2)/2,2)  ## average signal after spikein normalization
  nondepl1 <- round((nondepl1+nondepl2)/2,2)  ## average signal after spikein normalizationnnnnn
  
  pl <- rbind(data.frame(tracking_id=rownames(l[[1]]),condition="Depleted",signal=depl1[,1]),
              data.frame(tracking_id=rownames(l[[1]]),condition="Non-depleted",signal=nondepl1[,1]))
  
  return(pl)  
}

mymetagenefun <- function(l,genes,max.gene.len=4000, bw=10,Margin=500,facet=NULL){
  l[[1]] <- subset(l[[1]],rownames(l[[1]])%in%genes$tracking_id)
  l[[2]] <- subset(l[[2]],rownames(l[[2]])%in%genes$tracking_id)
  l[[3]] <- subset(l[[3]],rownames(l[[3]])%in%genes$tracking_id)
  l[[4]] <- subset(l[[4]],rownames(l[[4]])%in%genes$tracking_id)
  
  l[[1]] <- l[[1]][match(genes$tracking_id,rownames(l[[1]])),]
  l[[2]] <- l[[2]][match(genes$tracking_id,rownames(l[[2]])),]
  l[[3]] <- l[[3]][match(genes$tracking_id,rownames(l[[3]])),]
  l[[4]] <- l[[4]][match(genes$tracking_id,rownames(l[[4]])),]
  
  max.gene.len = round(max.gene.len/bw)
  Margin = round(Margin/bw)
  mcol= max.gene.len+2*Margin 
  
  depl1 <- matrix(NA,ncol=mcol,nrow = nrow(genes))
  depl2 <- matrix(NA,ncol=mcol,nrow = nrow(genes))
  nondepl1 <- matrix(NA,ncol=mcol,nrow = nrow(genes))
  nondepl2 <- matrix(NA,ncol=mcol,nrow = nrow(genes))
  ncol=(dim(l[[1]])[2])
  
  if(T){
    for(i in 1:nrow(genes)){
      j = ceiling(genes[i,]$width/bw)
      v <- l[[1]][i,(Margin+1):(Margin+j)]
      v <- vector.resizing(v,max.gene.len)
      v <- c(l[[1]][i,1:Margin],v,l[[1]][i,(ncol-Margin+1):(ncol)])
      #plot(1:length(v),v)
      depl1[i,] <- v
      v <- l[[2]][i,(Margin+1):(Margin+j)]
      v <- vector.resizing(v,max.gene.len)
      v <- c(l[[2]][i,1:Margin],v,l[[2]][i,(ncol-Margin+1):(ncol)])
      depl2[i,] <- v
      
      v <- l[[3]][i,(Margin+1):(Margin+j)]
      v <- vector.resizing(v,max.gene.len)
      v <- c(l[[3]][i,1:Margin],v,l[[3]][i,(ncol-Margin+1):(ncol)])
      nondepl1[i,] <- v
      v <- l[[4]][i,(Margin+1):(Margin+j)]
      v <- vector.resizing(v,max.gene.len)
      v <- c(l[[4]][i,1:Margin],v,l[[4]][i,(ncol-Margin+1):(ncol)])
      nondepl2[i,] <- v
    }
  }
 
  rownames(depl1) <- rownames(depl2) <- rownames(nondepl1) <- rownames(nondepl2) <-  rownames(l[[1]])
  depl1 <- round((depl1+depl2)/2,2)  ## average signal after spikein normalization
  nondepl1 <- round((nondepl1+nondepl2)/2,2)  ## average signal after spikein normalizationnnnnn
  
  if(!is.null(facet)){
    depl1 <-subset(depl1,rownames(depl1)%in%facet$tracking_id)
    nondepl1 <-subset(nondepl1,rownames(nondepl1)%in%facet$tracking_id)
  }
  
  get.avg <- function(m,nF){
    df_m <- cbind(1:ncol(m),
                  as.data.frame(apply(m,2,mean,na.rm=TRUE)),
                  as.data.frame(apply(m,2,sd,na.rm=TRUE)), 
                  as.data.frame(apply(m,2, function(x) sd(x,na.rm=TRUE)/sqrt(length(x))))
    )
    names(df_m) <- c("pos","signal", "sd","se")
    return(df_m)
  }
  pl <- rbind(cbind(get.avg(depl1),condition="Depleted"),
              cbind(get.avg(nondepl1),condition="Non-depleted")
  )
  return(pl)
}

mymetagenefun_tss <- function(l,genes,facet=NULL){
  l[[1]] <- subset(l[[1]],rownames(l[[1]])%in%genes$tracking_id)
  l[[2]] <- subset(l[[2]],rownames(l[[2]])%in%genes$tracking_id)
  l[[3]] <- subset(l[[3]],rownames(l[[3]])%in%genes$tracking_id)
  l[[4]] <- subset(l[[4]],rownames(l[[4]])%in%genes$tracking_id)
  
  depl1 <- l[[1]]
  depl2 <- l[[2]]
  nondepl1 <- l[[3]]
  nondepl2 <- l[[4]]

  depl1 <- round((depl1+depl2)/2,2)  ## average signal after spikein normalization
  nondepl1 <- round((nondepl1+nondepl2)/2,2)  ## average signal after spikein normalizationnnnnn
  
  if(!is.null(facet)){
    depl1 <-subset(depl1,rownames(depl1)%in%facet$tracking_id)
    nondepl1 <-subset(nondepl1,rownames(nondepl1)%in%facet$tracking_id)
  }
  
  get.avg <- function(m,nF){
    df_m <- cbind(1:ncol(m),
                  as.data.frame(apply(m,2,mean,na.rm=TRUE)),
                  as.data.frame(apply(m,2,sd,na.rm=TRUE)), 
                  as.data.frame(apply(m,2, function(x) sd(x,na.rm=TRUE)/sqrt(length(x))))
    )
    names(df_m) <- c("pos","signal", "sd","se")
    return(df_m)
  }
  pl <- rbind(cbind(get.avg(depl1),condition="Depleted"),
              cbind(get.avg(nondepl1),condition="Non-depleted")
  )
  return(pl)
}

myheatplotfun <- function(l,genes,g,outfile="temp.pdf",prot="RNAPII",
                          normF=c(1,1,1,1),max.gene.len=4000,bw=10,Margin=500,heat.ylab="verified protein coding genes",facet=NULL,xlab=c("","TSS","1kb","2kb","3kb","4kb","")){
  
  ## spikein normalize the signal
  l[[1]] <- l[[1]]*normF[1]
  l[[2]] <- l[[2]]*normF[2]
  l[[3]] <- l[[3]]*normF[3]
  l[[4]] <- l[[4]]*normF[4]
  ## These dataframes have default/hardcoded 500bp unscaled region coverage
  if(Margin<500){
    tyu <- 500-Margin
    ncol.temp <- ncol(l[[1]])
    l[[1]] <- l[[1]][,tyu:(ncol.temp-tyu)]
    l[[2]] <- l[[2]][,tyu:(ncol.temp-tyu)]
    l[[3]] <- l[[3]][,tyu:(ncol.temp-tyu)]
    l[[4]] <- l[[4]][,tyu:(ncol.temp-tyu)]
    rm(tyu, ncol.temp)
  }
  
  ## 
  max.gene.len = round(max.gene.len/bw)
  Margin = round(Margin/bw)
  
  nrow= length(genes$tracking_id)
  mcol <- max.gene.len+2*Margin
  a1 <- matrix(NA,ncol=mcol,nrow = nrow)  ## Depl, TSS
  b1 <- matrix(NA,ncol=mcol,nrow = nrow)  ## Depl, CPS
  a2 <- matrix(NA,ncol=mcol,nrow = nrow)  ## Non-depl, TSS
  b2 <- matrix(NA,ncol=mcol,nrow = nrow)  ## Non-depl, CPS
  ncol <- dim(l[[1]])[2]
  
  for(i in 1:length(genes$tracking_id)){
    v = l[[1]][i,]+l[[2]][i,]
    x <- tail(v[(Margin+1):(Margin+ceiling(genes$width[i]/bw)) ],max.gene.len)
    x <- c(v[1:Margin],x,v[(ncol-Margin+1):ncol])
    b1[i,] <- c(rep(0,(mcol-length(x))),x)
    
    v = l[[1]][i,]+l[[2]][i,]
    x <- head(v[(Margin+1):(Margin+ceiling(genes$width[i]/bw))],max.gene.len)
    x <- c(v[1:Margin],x,v[(ncol-Margin+1):ncol])
    a1[i,] <- c(x,rep(0,(mcol-length(x))))
    
    v = l[[3]][i,]+l[[4]][i,]
    x <- tail(v[(Margin+1):(Margin+ceiling(genes$width[i]/bw))],max.gene.len)
    x <- c(v[1:Margin],x,v[(ncol-Margin+1):ncol])
    b2[i,] <- c(rep(0,(mcol-length(x))),x)
    
    v = l[[3]][i,]+l[[4]][i,]
    x <- head(v[(Margin+1):(Margin+ceiling(genes$width[i]/bw))],max.gene.len)
    x <- c(v[1:Margin],x,v[(ncol-Margin+1):ncol])
    a2[i,] <- c(x,rep(0,(mcol-length(x))))
  }
  rownames(a1) <- rownames(a2) <-  rownames(b1) <- rownames(b2) <-  rownames(l[[1]])
  
  genes <- subset(genes,genes$tracking_id%in%g$tracking_id)
  a1 <- subset(a1,rownames(a1)%in%genes$tracking_id)
  a2 <- subset(a2,rownames(a2)%in%genes$tracking_id)
  b1 <- subset(b1,rownames(b1)%in%genes$tracking_id)
  b2 <- subset(b2,rownames(b2)%in%genes$tracking_id)
  
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
  
  a1 <- scale_row(a1) ##DEpl,TSS
  a2 <- scale_row(a2) ##Non-depl,TSS
  b1 <- scale_row(b1) ##
  b2 <- scale_row(b2)
  
  ## metagene
  pl <- mymetagenefun(l,genes,facet = NULL,max.gene.len = (max.gene.len*bw),bw = bw,Margin = (Margin*bw))
  pl$condition <- gsub("Depleted","Spn1-depleted",pl$condition)
  pl$condition <- gsub("Non-depleted","non-depleted",pl$condition)
  
  if(!is.null(facet)){
    pl1 <- mymetagenefun(l,genes,facet)
    pl1$condition <- gsub("Depleted","Spn1-depleted",pl1$condition)
    pl1$condition <- gsub("Non-depleted","non-depleted",pl1$condition)
    pl1$condition <- paste0(pl1$condition," (down)")
    pl <- rbind(pl,pl1)
  }
  
  ## heatmap
  m <- a
  m <- as.data.frame(m)
  m$tracking_id <- rownames(m)
  m <- melt(m)
  m$variable <- gsub("V","",m$variable)
  m$variable <- as.numeric(m$variable)
  m$value <- ifelse(m$value==0,NA,m$value)
  names(m)[3] <- "change"
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
  
  #scale_fill_gdistiller(palette = "YlOrBr",direction = 1,space = "Lab",na.value = "white")+
  heat_plot <- ggplot(m) + 
    geom_raster(aes(x=variable, y=tracking_id, fill=change),interpolate=F) +
    #scale_fill_gradient2(low="#01665e",mid="#f5f5f5",high = "#8c510a",midpoint = 0,na.value = "white",space = "Lab") +
    scale_fill_gradient2(low="#1b7837",mid="#f7f7f7",high = "#762a83",midpoint = 0,na.value = "white",space = "Lab") +
    scale_x_continuous(breaks = breaks,labels = xlab,expand = c(0,0))+
    geom_point(data = h,inherit.aes = F,aes(x = ed,y = y),col="black",size=0.05,stroke=0,alpha=0.1, shape=20)+
    geom_point(data = h,inherit.aes = F,aes(x = b,y = y),col="black",size=0.05,stroke=0,alpha=0.1,shape=20)+
    theme(panel.background = element_blank())+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(size=16,colour = "black"),
          axis.ticks.y = element_blank(),
          axis.text.y=element_blank(),
          axis.title.y = element_text(color="black",size=18))+
    theme(legend.position = "bottom",
          legend.title = element_text(size=16, color="black"),
          legend.text = element_text(size=14, face="plain"),
          legend.margin = margin(0,0,0,0),
          legend.key.width  = unit(1,"cm"),
          legend.box.margin = margin(0,0,0,0))+
    theme(axis.line.x = element_line(size = 0.7, colour = "black"),
          axis.line.y = element_blank())+
    ylab(paste0(nrow(a)," ",heat.ylab))
  
  meta <- ggplot(pl, aes(x=pos,y=signal,ymax=signal+se, ymin=signal-se,col=condition))+
    geom_vline(xintercept = c(Margin,(Margin+max.gene.len)),colour="gray50", linetype = "dashed") +
    geom_linerange(col='gray') + geom_line(size=1.2) +
    scale_color_manual(values = cols) +
    scale_y_continuous(expand = c(0,0),limits = c(0,NA))+
    scale_x_continuous(breaks =breaks, labels=xlab,expand = c(0,0)) +ylab("normalized signal") + xlab("") +
    theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = NA),axis.line = element_line(size = 0.7, colour = "black"))+
    theme(axis.text.x = element_text(colour = 'black',size=15),
          axis.text.y = element_text(colour = 'black',size=15), 
          axis.title.y = element_text(colour = 'black',size=16),
          strip.text.x = element_text(colour = 'black',size=16))+
    theme(legend.text =element_text(colour = 'black',size=14),
          legend.title = element_text(colour = 'black',size=15))+
    guides(color=guide_legend(title="Spn1 levels"))+
    theme(legend.text = element_text(size = 14,colour = "black"),legend.title = element_blank(),
          legend.key = element_blank(),legend.background = element_blank(),
          legend.position = c(1,1), 
          legend.justification = c(1, 1))
  
  #meta
  prot <- gsub("RNAPII",'Rpb1',prot)
  prot <- gsub("S2P",'Ser2-P',prot)
  prot <- gsub("S5P",'Ser5-P',prot)
  
  meta <- ggplotGrob(meta)
  heat_plot <- ggplotGrob(heat_plot)
  maxWidth = grid::unit.pmax(heat_plot$widths[2:5], meta$widths[2:5])
  meta$widths[2:5] <- as.list(maxWidth)
  heat_plot$widths[2:5] <- as.list(maxWidth)
  
  #ggsave(filename = paste0(outfile),)
  
  pdf(file = paste0(outfile),width = 6,height = 12)
  #grid.arrange(heat_plot,ncol=1,top= textGrob(paste0(prot),gp=gpar(fontsize=20,font=2)))
  grid.arrange(meta,heat_plot,ncol=1,heights=c(3,10),top= textGrob(paste0(prot),gp=gpar(fontsize=20,font=2)))
  #plot_grid(meta,heat_plot,ncol = 1,align = "v",rel_heights = c(3,10),labels = paste0(prot),label_colour = "black",label_size = 20,label_fontface = "bold")
  dev.off()
  
}

my_relative_signal <- function(l,l2,normF=c(1,1,1,1),normF2=c(1,1,1,1)){
  ## spikein normalization
  l[[1]] <- l[[1]]*normF[1]
  l[[2]] <- l[[2]]*normF[2]
  l[[3]] <- l[[3]]*normF[3]
  l[[4]] <- l[[4]]*normF[4]
  
  l2[[1]] <- l2[[1]]*normF2[1]
  l2[[2]] <- l2[[2]]*normF2[2]
  l2[[3]] <- l2[[3]]*normF2[3]
  l2[[4]] <- l2[[4]]*normF2[4]
  
  for(i in 1:4){
    l[[i]] = (l[[i]])/(l2[[i]])
    l[[i]][!is.finite(l[[i]])] <-0
    rownames(l[[i]]) <- rownames(l2[[i]])
  }
  
  return(l)
  
}

load(file="A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Gene groups.RData")
load("A:/work/WinstonLab/Natalia_Reim/ThesisFigures_March2019/NR_March19_10pFDR_downGenes.RData")
res <- read.delim("Normalization.factors.NataliaMayNov2017_CorrOrlando.txt",header = T)
rp <- rp[rp$gene!="RPP1",]
rp <- rp[-grep("MRP",rp$gene),]

#######-------------- Main body
########## (1) Merged replicates
if(F){
  prot="H3K4me3"
  protein<- c("H3 (round1)","H3","H3K4me3","H3K36me2","H3K36me3","RNAPII (round1)","RNAPII","S2P","S5P","Set2","Spn1","Spt6")
  protein<- c("RNAPII (round1)","RNAPII","S2P","S5P")
  for( prot in protein){
    f = paste("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_",prot,".RData",sep="")
    if(file.exists(f)){
      cat(prot,"\n")
      cat(res[res$Factor==prot,]$normF,"\n")
      load(f)
      normF=res[res$Factor==prot,]$normF
      if(prot %in% c("RNAPII (round1)","RNAPII","S2P","S5P")){
        normF=c(1,1,1,1)
        cat ("  using unit normalization factors\n")
      }
      #myheatplotfun(l = l,genes = genes,g = genes,prot = prot,outfile = paste0(OUTDIR,"PF_",prot,"_MetaMultiSpikeIn_ALLGenes_DEfacet.pdf"),normF = normF,facet=de.genes)
      myheatplotfun(l = l,genes = genes,g = genes,prot = prot,outfile = paste0(OUTDIR,"PF_",prot,"_MetaMultiSpikeIn_ALLGenes.pdf"),normF = normF,facet=NULL)
      #myheatplotfun(l = l,genes = genes,g = genes,prot = prot,outfile = paste0(OUTDIR,"PF_",prot,"_MetaMultiSpikeIn_ALLGenes_Scheme1.pdf"),normF = normF,facet=NULL)
      rm(l,f,prot,normF)
    }else{
      cat("file corresponding to ",prot, " does not exists! \n")
    }
  }
  
}

########## (2) Relative changes
if(F){
  protein <- list()
  protein[["S2P"]] <- "RNAPII"
  protein[["S5P"]] <- "RNAPII"
  protein[["H3K4me3"]] <- "H3"
  protein[["H3K36me2"]] <- "H3"
  protein[["H3K36me3"]] <- "H3 (round1)"
  #protein[["Spt6"]] <- "RNAPII (round1)"
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
      
      if(names(protein)[i]%in%c("S2P","S5P")){
        cat("using unit size factors\n")
        normF1=c(1,1,1,1)
        normF2=c(1,1,1,1)
      }else{
        normF1=res[res$Factor==names(protein)[[i]],]$normF
        normF2=res[res$Factor==protein[[i]],]$normF
      }
      
      l <-my_relative_signal(l1,l2,normF1,normF2)
      rm(l1,l2)
      myheatplotfun(l = l,genes = genes,g = genes,prot=names(protein)[i],outfile = paste0(OUTDIR,"PF_",names(protein)[i],"_MetaMultiSpikeIn_ALLGenes_RelSignal.pdf"),normF = c(1,1,1,1))
      #myheatplotfun(l = l,genes = genes,g = genes,prot=names(protein)[i],outfile = paste0(OUTDIR,"PF_",names(protein)[i],"_MetaMultiSpikeIn_ALLGenes_RelSignal_DEfacet.pdf"),normF = c(1,1,1,1),facet = de.genes)
      rm(l,f,prot,normF)
    }else{
      cat("file corresponding to ",names(protein)[i], " does not exists! \n")
    }
  }
  
}

########### (3) Merged replicates _ intron containing genes
if(F){
  myplot <- function(pl1,ylab="",legend.titl="",hline=F){
    require(ggplot2)
    q <- ggplot(pl1, aes(x=pos,y=signal,ymax=signal+se, ymin=signal-se,col=prot))
    if(hline==T){
      q <- q + geom_hline(yintercept = c(1),col="purple",size=1.2)
    }
    q <- q + geom_vline(xintercept = c(500,1000),colour="gray50", linetype = "dashed")
    q <- q + geom_linerange(col='gray') + geom_line(size=1.2) #+ ylim(-0.03,0.15) #ylim(-0.01,0.045) #ylim(-0.04,0.1)
    q <- q + scale_color_manual(values = cols)
    q <- q + facet_wrap(~condition, ncol=2,scales = "free")
    q <- q + scale_x_continuous(breaks = c(0,500,1000,1500), labels=c("","3'SS","5'SS","")) +ylab(ylab) + xlab("") 
    q <- q + theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = NA),axis.line = element_line(size = 0.7, colour = "gray50"))
    q <- q + theme(axis.text.x = element_text(colour = 'black',size=14),
                   axis.text.y = element_text(colour = 'black',size=14),
                   axis.title = element_text(colour = 'black',size=15,face="bold"),
                   strip.text.x = element_text(colour = 'black',size=16))
    q <- q + theme(legend.text =element_text(colour = 'black',size=14),
                   legend.title = element_text(colour = 'black',size=14,face="bold"),
                   legend.position = "bottom",
                   legend.justification = "center")
    q <- q + guides(color=guide_legend(title=legend.titl))
    return(q)
    
  }
  myplot_v2 <- function(pl1,ylab="",legend.titl="",hline=F,double.facet=F){
    require(ggplot2)
    q <- ggplot(pl1, aes(x=pos,y=signal,ymax=signal+se, ymin=signal-se,col=condition))
    if(hline==T){
      q <- q + geom_hline(yintercept = c(1),col="purple",size=1.2)
    }
    q <- q + geom_vline(xintercept = c(500,1000),colour="gray50", linetype = "dashed")
    q <- q + geom_linerange(col='gray') + geom_line(size=1.2) #+ ylim(-0.03,0.15) #ylim(-0.01,0.045) #ylim(-0.04,0.1)
    q <- q + scale_color_manual(values = cols)
    if(double.facet==T){
      q <- q + facet_wrap(~prot+category, ncol=6)
    }else{
      q <- q + facet_wrap(~prot, ncol=3,scales = "free")
    }
    q <- q + scale_x_continuous(breaks = c(0,500,1000,1500), labels=c("","3'SS","5'SS","")) +ylab(ylab) + xlab("") 
    q <- q + theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = NA),axis.line = element_line(size = 0.7, colour = "gray50"))
    q <- q + theme(axis.text.x = element_text(colour = 'black',size=14),
                   axis.text.y = element_text(colour = 'black',size=14),
                   axis.title = element_text(colour = 'black',size=15,face="bold"),
                   strip.text.x = element_text(colour = 'black',size=16))
    q <- q + theme(legend.text =element_text(colour = 'black',size=14),
                   legend.title = element_text(colour = 'black',size=14,face="bold"),
                   legend.position = "bottom",
                   legend.justification = "center")
    q <- q + guides(color=guide_legend(title=legend.titl))
    return(q)
    
  }
  
  introns <- read.delim(file="A:/work/WinstonLab/Natalia_Reim/October2016/gr/IR_Coveragebased/FeaturesForIRR.bed",header = F)
  names(introns) <- c("chr","start","end","id","tracking_id","strand","anno")
  introns <- introns[introns$anno=="introns",]
  introns$width <- introns$end-introns$start
  introns <- introns[order(introns$width,decreasing = F),]
  intr <- introns[!introns$chr%in%c("chrM") & introns$width>50,]
  intr$gene <- intr$tracking_id
  introns$gene <- introns$tracking_id
  intr$tracking_id <- intr$id
  introns$tracking_id <- introns$id
  
  protein<- c("H3 (round1)","H3","H3K4me3","H3K36me2","H3K36me3","RNAPII (round1)","RNAPII","S2P","S5P","Set2","Spn1","Spt6")
  pl <- c()
  gl <- c()
  for( prot in protein){
    f = paste("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_Introns_CCbp",prot,".RData",sep="")
    #f = paste("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_Introns_bp",prot,".RData",sep="")
    if(file.exists(f)){
      cat(prot,"\n")
      cat(res[res$Factor==prot,]$normF,"\n")
      load(f)
      normF=res[res$Factor==prot,]$normF
      if(prot %in% c("RNAPII (round1)","RNAPII","S2P","S5P")){
        normF=c(1,1,1,1)
        cat ("  using unit normalization factors\n")
      }
      
      
      rownames(l[[1]]) <- rownames(l[[2]]) <- rownames(l[[3]]) <- rownames(l[[4]]) <- introns$id
      l[[1]] <- l[[1]]*normF[1]
      l[[2]] <- l[[2]]*normF[2]
      l[[3]] <- l[[3]]*normF[3]
      l[[4]] <- l[[4]]*normF[4]
      myheatplotfun(l = l,genes = introns,g = intr,prot = prot,outfile = paste0(OUTDIR,"PF_",prot,"_MetaMultiSpikeIn_Introns_CC1bp_Signal.pdf"),normF = c(1,1,1,1),facet=NULL,max.gene.len=500,bw=1,Margin=500,heat.ylab="intron containing genes",xlab = c("","3'SS","5'SS",""))
      #myheatplotfun(l = l,genes = introns,g = intr,prot = prot,outfile = paste0(OUTDIR,"PF_",prot,"_MetaMultiSpikeIn_Introns_1bp_Signal.pdf"),normF = c(1,1,1,1),facet=NULL,max.gene.len=500,bw=1,Margin=500,heat.ylab="intron containing genes",xlab = c("","3'SS","5'SS",""))
      
      if(T){
        pl <- rbind(pl,cbind(mymetagenefun(l,genes=intr,facet = NULL,max.gene.len = 500,bw = 1,Margin = 500),prot=prot,category="All-ICG",quant="absolute"))
        pl <- rbind(pl,cbind(mymetagenefun(l,genes=subset(intr,intr$gene%in%rp$tracking_id),facet = NULL,max.gene.len = 500,bw = 1,Margin = 500),prot=prot,category="RP-ICG",quant="absolute"))
        pl <- rbind(pl,cbind(mymetagenefun(l,genes=subset(intr,!intr$gene%in%rp$tracking_id),facet = NULL,max.gene.len = 500,bw = 1,Margin = 500),prot=prot,category="nonRP-ICG",quant="absolute"))
        
        gl <- rbind(gl,cbind(mytotalgenesignal(l,genes=intr,max.gene.len = 500,bw = 1,Margin = 500),prot=prot,category="All-ICG",quant="absolute"))
        gl <- rbind(gl,cbind(mytotalgenesignal(l,genes=subset(intr,intr$gene%in%rp$tracking_id),max.gene.len = 500,bw = 1,Margin = 500),prot=prot,category="RP-ICG",quant="absolute"))
        gl <- rbind(gl,cbind(mytotalgenesignal(l,genes=subset(intr,!intr$gene%in%rp$tracking_id),max.gene.len = 500,bw = 1,Margin = 500),prot=prot,category="nonRP-ICG",quant="absolute"))
        
        l1 <- list()
        m <- l[[1]]/l[[3]]
        m[!is.finite(m)] <- 1
        l1[["Depleted-1"]] <- m
        m <- l[[2]]/l[[4]]
        m[!is.finite(m)] <- 1
        l1[["Depleted-2"]] <- m
        l1[["Non-depleted-1"]] <- (l[[1]]*0+1)
        l1[["Non-depleted-2"]] <- (l[[1]]*0+1)
        pl <- rbind(pl,cbind(mymetagenefun(l1,genes=intr,facet = NULL,max.gene.len = 500,bw = 1,Margin = 500),prot=prot,category="All-ICG",quant="depl/nondepl"))
        pl <- rbind(pl,cbind(mymetagenefun(l1,genes=subset(intr,intr$gene%in%rp$tracking_id),facet = NULL,max.gene.len = 500,bw = 1,Margin = 500),prot=prot,category="RP-ICG",quant="depl/nondepl"))
        pl <- rbind(pl,cbind(mymetagenefun(l1,genes=subset(intr,!intr$gene%in%rp$tracking_id),facet = NULL,max.gene.len = 500,bw = 1,Margin = 500),prot=prot,category="nonRP-ICG",quant="depl/nondepl"))
      }
      rm(l,l1,f,prot,normF)
    }else{
      cat("file corresponding to ",prot, " does not exists! \n")
    }
  }
  
  cols=c("H3"="darkolivegreen1", "H3 (round1)"="darkolivegreen3", 
         "H3K36me2"="royalblue1", "H3K36me3"="royalblue3", 
         "H3K4me3"="royalblue4", "Spt6"="black", "Ser2-P"="firebrick1", 
         "Ser5-P"="firebrick2","Rpb1"="firebrick3","Rpb1 (round1)"="firebrick4",
         "Spn1"="gray","Set2"="brown",
         "Spn1-depleted" = "#BB5566","non-depleted" = "#4477AA",
         "Non-Intron genes"="pink","All-ICG"="brown",
         "RP-ICG"="red","nonRP-ICG"="pink"         )
  
  gl$prot <- gsub("S2P","Ser2-P",gl$prot) 
  gl$prot <- gsub("S5P","Ser5-P",gl$prot) 
  gl$prot <- gsub("RNAPII","Rpb1",gl$prot) 
  gl$condition <- gsub("Depleted","Spn1-depleted",gl$condition) 
  gl$condition <- gsub("Non-depleted","non-depleted",gl$condition) 
  pl$prot <- gsub("S2P","Ser2-P",pl$prot) 
  pl$prot <- gsub("S5P","Ser5-P",pl$prot) 
  pl$prot <- gsub("RNAPII","Rpb1",pl$prot) 
  pl$condition <- gsub("Depleted","Spn1-depleted",pl$condition) 
  pl$condition <- gsub("Non-depleted","non-depleted",pl$condition) 
  
  ## Plot 1, facet by Spn1 depletion
  pdf(paste0(OUTDIR,"PF_ScaledIntronMetagenes_AllProteinSignal.pdf"),width = 9,height = 5)
   pl1 <-pl[pl$prot%in%c("H3 (round1)", "H3", "H3K4me3", "H3K36me2", "H3K36me3","Set2") & pl$category=="All-ICG" & pl$quant == "absolute",]
   myplot(pl1,ylab = "normalized signal",legend.titl = "Spn1 levels")+ggtitle("All intron containing genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
   pl1 <-pl[pl$prot%in%c("Rpb1 (round1)", "Rpb1", "Ser2-P", "Ser5-P", "Spn1", "Spt6") & pl$category=="All-ICG" & pl$quant == "absolute",]
   myplot(pl1,ylab = "normalized signal",legend.titl = "Spn1 levels")+ggtitle("All intron containing genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
   pl1 <-pl[pl$prot%in%c("H3 (round1)", "H3", "H3K4me3", "H3K36me2", "H3K36me3","Set2") & pl$category=="RP-ICG" & pl$quant == "absolute",]
   myplot(pl1,ylab = "normalized signal",legend.titl = "Spn1 levels")+ggtitle("RP[LS] intron containing genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
   pl1 <-pl[pl$prot%in%c("Rpb1 (round1)", "Rpb1", "Ser2-P", "Ser5-P", "Spn1", "Spt6") & pl$category=="RP-ICG" & pl$quant == "absolute",]
   myplot(pl1,ylab = "normalized signal",legend.titl = "Spn1 levels")+ggtitle("RP[LS] intron containing genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
   pl1 <-pl[pl$prot%in%c("H3 (round1)", "H3", "H3K4me3", "H3K36me2", "H3K36me3","Set2") & pl$category=="nonRP-ICG" & pl$quant == "absolute",]
   myplot(pl1,ylab = "normalized signal",legend.titl = "Spn1 levels")+ggtitle("nonRP[LS] intron containing genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
   pl1 <-pl[pl$prot%in%c("Rpb1 (round1)", "Rpb1", "Ser2-P", "Ser5-P", "Spn1", "Spt6") & pl$category=="nonRP-ICG" & pl$quant == "absolute",]
   myplot(pl1,ylab = "normalized signal",legend.titl = "Spn1 levels")+ggtitle("nonRP[LS] intron containing genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  dev.off()
  
  ## Plot2, facet by factors
  #pdf(paste0(OUTDIR,"PF_ScaledIntronMetagenes_AllFactorsCC.pdf"),width = 9,height = 9)
  pdf(paste0(OUTDIR,"PF_ScaledIntronMetagenes_AllFactors.pdf"),width = 9,height = 9)
  ## all genes normalized signal
    pl1 <-pl[pl$quant=="absolute" & pl$category=="All-ICG",]
    myplot_v2(pl1,ylab = "normalized signal",legend.titl = "Spn1 levels")+ggtitle("All intron containing genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
    ## all genes normalized signal
    pl1 <-pl[pl$quant=="absolute" & pl$category=="RP-ICG",]
    myplot_v2(pl1,ylab = "normalized signal",legend.titl = "Spn1 levels")+ggtitle("RP[LS] intron containing genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
    pl1 <-pl[pl$quant=="absolute" & pl$category=="nonRP-ICG",]
    myplot_v2(pl1,ylab = "normalized signal",legend.titl = "Spn1 levels")+ggtitle("nonRP[LS] intron containing genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
    pl1 <-pl[pl$quant=="depl/nondepl" & pl$condition=="Spn1-depleted" & pl$category!="All-ICG",]
    pl1$condition <- pl1$category
    myplot_v2(pl1,ylab = "normalized signal (depleted/non-depleted)",legend.titl = "Gene class",hline = T)+ggtitle("RP[LS] intron containing genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
    pl1 <-pl[pl$quant=="depl/nondepl" & pl$condition=="Spn1-depleted",]
    pl1$condition <- pl1$category
    myplot_v2(pl1,ylab = "normalized signal (depleted/non-depleted)",legend.titl = "Gene class",hline = T)+ggtitle("RP[LS] intron containing genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  dev.off()
  
  ## occupancy
  gl <- reshape2::dcast(gl,tracking_id+prot+category+quant~condition,value.var = "signal")
  gl$l2fc <- log2((gl$`Spn1-depleted`+1)/(gl$`non-depleted`+1))
  w <- ggplot(gl, aes(x=category,y=l2fc,fill=category))
  w <- w + geom_violin(draw_quantiles = c(0.5))
  w <- w + scale_fill_manual(values = cols)
  w <- w + xlab("") + ylab("occupancy fold change, log2")
  w <- w + facet_wrap(~prot,nrow = 2,scales = "free_x")+theme(strip.text = element_text(size=16,colour = "black"))
  w <- w + geom_hline(yintercept = c(0,seq(min(gl$l2fc), max(gl$l2fc),1)),linetype="dashed",col="gray")
  w <- w + theme(panel.background = element_rect(fill = NA),
                 axis.line = element_line(size = 0.7, colour = "gray50"))
  w <- w + theme(axis.text.x = element_blank(),
                 axis.text.y = element_text(colour = 'black',size=15),
                 axis.title = element_text(colour = 'black',size=16))
  w <- w + theme(legend.text =element_text(colour = 'black',size=15),
                 legend.title = element_text(colour = 'black',size=15,face="bold"),
                 legend.position = "bottom",
                 legend.justification = "center")
  pdf(paste0(OUTDIR,"PF_ChIPSeqOccup_Violine_Introns.pdf"),width = 10,height = 7)
  w
  dev.off()
  
  g1 <- gl[gl$category=="All-ICG",]
  g1 <- reshape2::dcast(g1,tracking_id~prot,value.var = "l2fc")
  g1$rp <- ifelse(g1$tracking_id %in% as.character(gl[gl$category=="nonRP-ICG",]$tracking_id),as.character("RP-gene"),as.character("nonRP-gene") )
  g1 <- g1[,c(1,14,2:13)]
  require(xlsx)
  write.xlsx(g1,file=paste0(OUTDIR,"ICG_IntronLevelFactorOccupancy_CorrOrlando.xlsx"),sheetName = "AbsoluteSignal",showNA = T,row.names = F)
  
}


########### (4) Merged replicates _ intron containing genes relative change anlysis
if(F){
  introns <- read.delim(file="A:/work/WinstonLab/Natalia_Reim/October2016/gr/IR_Coveragebased/FeaturesForIRR.bed",header = F)
  names(introns) <- c("chr","start","end","id","tracking_id","strand","anno")
  introns <- introns[introns$anno=="introns",]
  introns$width <- introns$end-introns$start
  introns <- introns[order(introns$width,decreasing = F),]
  intr <- introns[!introns$chr%in%c("chrM") & introns$width>50,]
  intr$gene <- intr$tracking_id
  introns$gene <- introns$tracking_id
  intr$tracking_id <- intr$id
  introns$tracking_id <- introns$id
  
  protein <- list()
  protein[["S2P"]] <- "RNAPII"
  protein[["S5P"]] <- "RNAPII"
  protein[["H3K4me3"]] <- "H3"
  protein[["H3K36me2"]] <- "H3"
  protein[["H3K36me3"]] <- "H3 (round1)"
  #protein[["Spt6"]] <- "RNAPII (round1)"
  
  pl <- c()
  gl <- c()
  for(i in  1:length(protein)){
    #f = paste("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_Introns_CCbp",names(protein)[i],".RData",sep="")
    #fc = paste("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_Introns_CCbp",protein[[i]],".RData",sep="")
    f = paste("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_Introns_bp",names(protein)[i],".RData",sep="")
    fc = paste("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_Introns_bp",protein[[i]],".RData",sep="")
    
    if(file.exists(f) & file.exists(fc)){
      cat(names(protein)[i],"\n")
      cat(res[res$Factor==names(protein)[i],]$normF,"\n")
      cat(protein[[i]],"\n")
      cat(res[res$Factor==protein[[i]],]$normF,"\n")
      
      load(f)
      l1 <- l; rm(l);
      load(fc)
      l2 <- l; rm(l);
      
      if(names(protein)[i]%in%c("S2P","S5P")){
        normF1=c(1,1,1,1)
        normF2=c(1,1,1,1)
      }else{
        normF1=res[res$Factor==names(protein)[[i]],]$normF
        normF2=res[res$Factor==protein[[i]],]$normF
      }
      l <-my_relative_signal(l1,l2,normF1,normF2)
      rm(l1,l2)
      
      rownames(l[[1]]) <- rownames(l[[2]]) <- rownames(l[[3]]) <- rownames(l[[4]]) <- introns$id
      #myheatplotfun(l = l,genes = introns,g = intr,prot=names(protein)[i],outfile = paste0(OUTDIR,"PF_",names(protein)[i],"_MetaMultiSpikeIn_Introns_CC1bp_RelSignal.pdf"),normF = c(1,1,1,1),max.gene.len=500,bw=1,Margin=500,heat.ylab="intron containing genes",xlab = c("","3'SS","5'SS",""))
      myheatplotfun(l = l,genes = introns,g = intr,prot=names(protein)[i],outfile = paste0(OUTDIR,"PF_",names(protein)[i],"_MetaMultiSpikeIn_Introns_1bp_RelSignal.pdf"),normF = c(1,1,1,1),max.gene.len=500,bw=1,Margin=500,heat.ylab="intron containing genes",xlab = c("","3'SS","5'SS",""))
      
      if(T){
        pl <- rbind(pl,cbind(mymetagenefun(l,genes=intr,facet = NULL,max.gene.len = 500,bw = 1,Margin = 500),prot=names(protein)[i],category="All-ICG",quant="absolute"))
        pl <- rbind(pl,cbind(mymetagenefun(l,genes=subset(intr,intr$gene%in%rp$tracking_id),facet = NULL,max.gene.len = 500,bw = 1,Margin = 500),prot=names(protein)[i],category="RP-ICG",quant="absolute"))
        pl <- rbind(pl,cbind(mymetagenefun(l,genes=subset(intr,!intr$gene%in%rp$tracking_id),facet = NULL,max.gene.len = 500,bw = 1,Margin = 500),prot=names(protein)[i],category="nonRP-ICG",quant="absolute"))
        
        gl <- rbind(gl,cbind(mytotalgenesignal(l,genes=intr,max.gene.len = 500,bw = 1,Margin = 500),prot=names(protein)[i],category="All-ICG",quant="absolute"))
        gl <- rbind(gl,cbind(mytotalgenesignal(l,genes=subset(intr,intr$gene%in%rp$tracking_id),max.gene.len = 500,bw = 1,Margin = 500),prot=names(protein)[i],category="RP-ICG",quant="absolute"))
        gl <- rbind(gl,cbind(mytotalgenesignal(l,genes=subset(intr,!intr$gene%in%rp$tracking_id),max.gene.len = 500,bw = 1,Margin = 500),prot=names(protein)[i],category="nonRP-ICG",quant="absolute"))
        
        l1 <- list()
        m <- l[[1]]/l[[3]]
        m[!is.finite(m)] <- 1
        l1[["Depleted-1"]] <- m
        m <- l[[2]]/l[[4]]
        m[!is.finite(m)] <- 1
        l1[["Depleted-2"]] <- m
        l1[["Non-depleted-1"]] <- (l[[1]]*0+1)
        l1[["Non-depleted-2"]] <- (l[[1]]*0+1)
        pl <- rbind(pl,cbind(mymetagenefun(l1,genes=intr,facet = NULL,max.gene.len = 500,bw = 1,Margin = 500),prot=names(protein)[i],category="All-ICG",quant="depl/nondepl"))
        pl <- rbind(pl,cbind(mymetagenefun(l1,genes=subset(intr,intr$gene%in%rp$tracking_id),facet = NULL,max.gene.len = 500,bw = 1,Margin = 500),prot=names(protein)[i],category="RP-ICG",quant="depl/nondepl"))
        pl <- rbind(pl,cbind(mymetagenefun(l1,genes=subset(intr,!intr$gene%in%rp$tracking_id),facet = NULL,max.gene.len = 500,bw = 1,Margin = 500),prot=names(protein)[i],category="nonRP-ICG",quant="depl/nondepl"))
        
      }

      rm(l,l1,f,prot,normF)
    }else{
      cat("file corresponding to ",names(protein)[i], " does not exists! \n")
    }
  }

  gl$prot <- gsub("S2P","Ser2-P",gl$prot) 
  gl$prot <- gsub("S5P","Ser5-P",gl$prot) 
  gl$prot <- gsub("RNAPII","Rpb1",gl$prot) 
  gl$condition <- gsub("Depleted","Spn1-depleted",gl$condition) 
  gl$condition <- gsub("Non-depleted","non-depleted",gl$condition) 
  pl$prot <- gsub("S2P","Ser2-P",pl$prot) 
  pl$prot <- gsub("S5P","Ser5-P",pl$prot) 
  pl$prot <- gsub("RNAPII","Rpb1",pl$prot) 
  pl$condition <- gsub("Depleted","Spn1-depleted",pl$condition) 
  pl$condition <- gsub("Non-depleted","non-depleted",pl$condition) 
  
  cols = c("Depleted-1" = "pink","Depleted-2" = "red2",
    "Non-depleted-1" = "lightblue","Non-depleted-2" = "slateblue2",
    "Depleted" = "red","Non-depleted" = "slateblue",
    "OPEN genes"="green","CLOSED genes"="red",
    "High Pol2"="red4","Medium Pol2"="red","Low Pol2"="pink",
    "<0.5kb"="yellow","0.5kb-1kb"="orange","1kb-2kb"="red",">2kb"="red4",
    "non-DE genes"="pink","non-RP genes"="pink","DE genes"= "red", "RP genes"="red",
    "Non-Intron genes"="pink","All-ICG"="brown",
    "RP-ICG"="red","nonRP-ICG"="pink",
    "Spn1-depleted" = "#BB5566","non-depleted" = "#4477AA")
  myplot <- function(pl1,ylab="",legend.titl="",hline=F){
    require(ggplot2)
    q <- ggplot(pl1, aes(x=pos,y=signal,ymax=signal+se, ymin=signal-se,col=condition))
    if(hline==T){
      q <- q + geom_hline(yintercept = c(1),col="purple",size=1.2)
    }
    q <- q + geom_vline(xintercept = c(500,1000),colour="gray50", linetype = "dashed")
    q <- q + geom_linerange(col='gray') + geom_line(size=1.2) #+ ylim(-0.03,0.15) #ylim(-0.01,0.045) #ylim(-0.04,0.1)
    q <- q + scale_color_manual(values = cols)
    q <- q + facet_wrap(~prot, ncol=3,scales = "free")
    q <- q + scale_x_continuous(breaks = c(0,500,1000,1500), labels=c("","3'SS","5'SS","")) +ylab(ylab) + xlab("") 
    q <- q + theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = NA),axis.line = element_line(size = 0.7, colour = "gray50"))
    q <- q + theme(axis.text.x = element_text(colour = 'black',size=14),
                   axis.text.y = element_text(colour = 'black',size=14),
                   axis.title = element_text(colour = 'black',size=15,face="bold"),
                   strip.text.x = element_text(colour = 'black',size=16))
    q <- q + theme(legend.text =element_text(colour = 'black',size=14),
                   legend.title = element_text(colour = 'black',size=14,face="bold"),
                   legend.position = "bottom",
                   legend.justification = "center")
    q <- q + guides(color=guide_legend(title=legend.titl))
    return(q)
    
  }
  
  #pdf(paste0(OUTDIR,"PF_ScaledIntronMetagenes_CC.pdf"),width = 9,height = 6)
  pdf(paste0(OUTDIR,"PF_ScaledIntronMetagenes.pdf"),width = 9,height = 6)
  ## all genes normalized signal
  pl1 <-pl[pl$quant=="absolute" & pl$category=="All-ICG",]
  myplot(pl1,ylab = "normalized signal",legend.titl = "Spn1 levels")+ggtitle("All intron containing genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  ## all genes normalized signal
  pl1 <-pl[pl$quant=="absolute" & pl$category=="RP-ICG",]
  myplot(pl1,ylab = "normalized signal",legend.titl = "Spn1 levels")+ggtitle("RP[LS] intron containing genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  pl1 <-pl[pl$quant=="absolute" & pl$category=="nonRP-ICG",]
  myplot(pl1,ylab = "normalized signal",legend.titl = "Spn1 levels")+ggtitle("nonRP[LS] intron containing genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  pl1 <-pl[pl$quant=="depl/nondepl" & pl$condition=="Spn1-depleted" & pl$category!="All-ICG",]
  pl1$condition <- pl1$category
  myplot(pl1,ylab = "normalized signal (depleted/non-depleted)",legend.titl = "Gene class",hline = T)+ggtitle("RP[LS] intron containing genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  pl1 <-pl[pl$quant=="depl/nondepl" & pl$condition=="Spn1-depleted",]
  pl1$condition <- pl1$category
  myplot(pl1,ylab = "normalized signal (depleted/non-depleted)",legend.titl = "Gene class",hline = T)+ggtitle("RP[LS] intron containing genes")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  dev.off()
  
  
  ## relative ChIP occupancy
  gl <- reshape2::dcast(gl,tracking_id+prot+category+quant~condition,value.var = "signal")
  gl$l2fc <- log2((gl$`Spn1-depleted`+1)/(gl$`non-depleted`+1))
  
  w <- ggplot(gl, aes(x=category,y=l2fc,fill=category))
  w <- w + geom_violin(draw_quantiles = c(0.5))
  w <- w + scale_fill_manual(values = cols)
  w <- w + xlab("") + ylab("occupancy fold change, log2")
  w <- w + facet_wrap(~prot,nrow = 2,scales = "free_x")+theme(strip.text = element_text(size=16,colour = "black"))
  w <- w + geom_hline(yintercept = c(0,seq(min(gl$l2fc), max(gl$l2fc),1)),linetype="dashed",col="gray")
  w <- w + theme(panel.background = element_rect(fill = NA),
                 axis.line = element_line(size = 0.7, colour = "gray50"))
  w <- w + theme(axis.text.x = element_blank(),
                 axis.text.y = element_text(colour = 'black',size=15),
                 axis.title = element_text(colour = 'black',size=16))
  w <- w + theme(legend.text =element_text(colour = 'black',size=15),
                 legend.title = element_text(colour = 'black',size=15,face="bold"),
                 legend.position = "bottom",
                 legend.justification = "center")
  w
  pdf(paste0(OUTDIR,"PF_ChIPSeqOccup_Violine_Introns_relSignal.pdf"),width = 9,height = 6)
  w
  dev.off()
  
  
  g1 <- gl[gl$category=="All-ICG",]
  g1 <- reshape2::dcast(g1,tracking_id~prot,value.var = "l2fc")
  g1$rp <- ifelse(g1$tracking_id %in% as.character(gl[gl$category=="nonRP-ICG",]$tracking_id),as.character("RP-gene"),as.character("nonRP-gene") )
  g1 <- g1[,c(1,7,2:6)]
  write.xlsx(g1,file=paste0(OUTDIR,"ICG_IntronLevelFactorOccupancy_CorrOrlando.xlsx"),sheetName = "RelativeNormalized",showNA = T,row.names = F,append = T)
  
}


########### (5) Assessing specificities of intron containing genes using relative changes over entire genes
if(F){
  expn <- read.delim("A:/work/WinstonLab/Natalia_Reim/ThesisFigures_March2019/RNASeq_WTComparison.txt",header=T)
  expn <- expn[,c("tracking_id","WT_DMSO")]
  expn$WT_DMSO <- log2(expn$WT_DMSO)
  expn$WT_DMSO <- ifelse(!is.finite(expn$WT_DMSO),0,expn$WT_DMSO)
  
  genes <- merge(genes,expn,by="tracking_id")
  names(genes)[20] <- "wt_expn"
  genes$rp <- ifelse(genes$tracking_id%in%rp$tracking_id,1,0)
  genes$intron <- ifelse(genes$tracking_id%in%introns$tracking_id,1,0)
  genes$rpintrons <- ifelse(genes$rp>0&genes$intron>0,1,0)
  
  genes$length_class <- Hmisc::cut2(log2(genes$width),g=50)
  genes$expn_class <- Hmisc::cut2(log2(genes$wt_expn),g=50)
  genes$comb_class <- paste0(Hmisc::cut2(log2(genes$wt_expn),g=25), "_",Hmisc::cut2(log2(genes$wt_expn),g=25))
  
  get.gene.list <- function(rp1,g1,mult=1){
    len.class.genes <- c()
    for(i in 1:length(rp1)){
      if(nrow(rp1[[i]])>0){
        len = nrow(rp1[[i]])
        len=len*mult
        ids =sample(nrow(g1[[names(rp1)[i]]]),size = len,replace = F)
        len.class.genes <- c(len.class.genes,g1[[names(rp1)[i]]][ids,]$tracking_id)
        rm(len,ids)
      }
    }
    return(len.class.genes)
  }
  
  ### All RP genes
  rp <- genes[genes$rp>0,]
  g <- genes[genes$rp==0,]
   ## length class as RP genes
  rp1 <- split(rp,rp$length_class)
  g1 <- split(g,g$length_class)
  len.class.genes <- get.gene.list(rp1,g1)
    ## expn class as RP genes
  rp1 <- split(rp,rp$expn_class)
  g1 <- split(g,g$expn_class)
  expn.class.genes <- get.gene.list(rp1,g1)
    ## expn class as RP genes
  rp1 <- split(rp,rp$comb_class)
  g1 <- split(g,g$comb_class)
  comb.class.genes <- get.gene.list(rp1,g1)
  
  
  ### Only RP intron genes
  rp <- genes[genes$rpintrons>0,]
  g <- genes[genes$rpintrons==0,]
  ## length class as RP genes
  rp1 <- split(rp,rp$length_class)
  g1 <- split(g,g$length_class)
  len.class.genes.rpI <- get.gene.list(rp1,g1)
  ## expn class as RP genes
  rp1 <- split(rp,rp$expn_class)
  g1 <- split(g,g$expn_class)
  expn.class.genes.rpI <- get.gene.list(rp1,g1)
  ## expn class as RP genes
  rp1 <- split(rp,rp$comb_class)
  g1 <- split(g,g$comb_class)
  comb.class.genes.rpI <- get.gene.list(rp1,g1)
  rm(rp,rp1,g,g1)
  rp <- genes[genes$rp>0,]
  
  
  protein <- list()
  protein[["S2P"]] <- "RNAPII"
  protein[["S5P"]] <- "RNAPII"
  protein[["H3K4me3"]] <- "H3"
  protein[["H3K36me2"]] <- "H3"
  protein[["H3K36me3"]] <- "H3 (round1)"
  #protein[["Spt6"]] <- "RNAPII (round1)"
  rp$introns <- ifelse(rp$tracking_id%in%introns$tracking_id,1,0)
  table(rp$introns)
  
  pl <- c()
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
      
      if(names(protein)[i]%in%c("S2P","S5P")){
        normF1=c(1,1,1,1)
        normF2=c(1,1,1,1)
      }else{
        normF1=res[res$Factor==names(protein)[[i]],]$normF
        normF2=res[res$Factor==protein[[i]],]$normF
      }
      l <-my_relative_signal(l1,l2,normF1,normF2)
      rm(l1,l2)
      
      pl <- rbind(pl,cbind(mymetagenefun(l,genes=rp),prot=names(protein)[i],category="RP genes",quant="abso"))
      pl <- rbind(pl,cbind(mymetagenefun(l,genes=rp[rp$introns>0,]),prot=names(protein)[i],category="RP-ICG",quant="abso"))
      pl <- rbind(pl,cbind(mymetagenefun(l,genes=rp[rp$introns==0,]),prot=names(protein)[i],category="RP-nonICG",quant="abso"))
      pl <- rbind(pl,cbind(mymetagenefun(l,genes=genes[genes$tracking_id%in%len.class.genes,]),prot=names(protein)[i],category="length-matched nonRP genes",quant="abso"))
      pl <- rbind(pl,cbind(mymetagenefun(l,genes=genes[genes$tracking_id%in%expn.class.genes,]),prot=names(protein)[i],category="expn-matched nonRP genes",quant="abso"))
      pl <- rbind(pl,cbind(mymetagenefun(l,genes=genes[genes$tracking_id%in%comb.class.genes,]),prot=names(protein)[i],category="length+expn-matched nonRP genes",quant="abso"))
      
      l1 <- list()
      m <- l[[1]]/l[[3]]
      m[!is.finite(m)] <- 1
      l1[["Depleted-1"]] <- m
      m <- l[[2]]/l[[4]]
      m[!is.finite(m)] <- 1
      l1[["Depleted-2"]] <- m
      l1[["Non-depleted-1"]] <- (l[[1]]*0+1)
      l1[["Non-depleted-2"]] <- (l[[1]]*0+1)
      pl <- rbind(pl,cbind(mymetagenefun(l1,genes=rp),prot=names(protein)[i],category="RP genes",quant="depl/nondepl"))
      pl <- rbind(pl,cbind(mymetagenefun(l1,genes=rp[rp$introns>0,]),prot=names(protein)[i],category="RP-ICG",quant="depl/nondepl"))
      pl <- rbind(pl,cbind(mymetagenefun(l1,genes=rp[rp$introns==0,]),prot=names(protein)[i],category="RP-nonICG",quant="depl/nondepl"))
      pl <- rbind(pl,cbind(mymetagenefun(l1,genes=genes[genes$tracking_id%in%len.class.genes,]),prot=names(protein)[i],category="length-matched nonRP genes",quant="depl/nondepl"))
      pl <- rbind(pl,cbind(mymetagenefun(l1,genes=genes[genes$tracking_id%in%expn.class.genes,]),prot=names(protein)[i],category="expn-matched nonRP genes",quant="depl/nondepl"))
      pl <- rbind(pl,cbind(mymetagenefun(l1,genes=genes[genes$tracking_id%in%comb.class.genes,]),prot=names(protein)[i],category="length+expn-matched nonRP genes",quant="depl/nondepl"))
      rm(l,l1,f,prot,normF)
    }else{
      cat("file corresponding to ",names(protein)[i], " does not exists! \n")
    }
  }
  pl$prot <- gsub("S2P","Ser2-P",pl$prot) 
  pl$prot <- gsub("S5P","Ser5-P",pl$prot) 
  pl$prot <- gsub("RNAPII","Rpb1",pl$prot) 
  pl$condition <- gsub("Depleted","Spn1-depleted",pl$condition) 
  pl$condition <- gsub("Non-depleted","non-depleted",pl$condition) 
  
  myplot <- function(pl1,ylab="",legend.titl="",hline=F){
    require(ggplot2)
    q <- ggplot(pl1, aes(x=pos,y=signal,ymax=signal+1.96*se, ymin=signal-1.96*se,col=category))
    if(hline==T){
      q <- q + geom_hline(yintercept = c(1),col="purple",size=1.2)
    }
    q <- q + facet_wrap(~prot, ncol=3,scales = "free")
    q <- q + geom_vline(xintercept = c(50,450),colour="gray50", linetype = "dashed")
    q <- q + geom_linerange(col='gray') + geom_line(size=1.2) #+ ylim(-0.03,0.15) #ylim(-0.01,0.045) #ylim(-0.04,0.1)
    #q <- q + scale_color_manual(values = cols)
    q <- q + scale_color_brewer(palette = "Dark2")
    q <- q + scale_x_continuous(breaks = c(0,50,450,500), labels=c("","TSS","CPS","")) +ylab(ylab) + xlab("") 
    q <- q + theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = NA),axis.line = element_line(size = 0.7, colour = "gray50"))
    q <- q + theme(axis.text.x = element_text(colour = 'black',size=14),
                   axis.text.y = element_text(colour = 'black',size=14),
                   axis.title = element_text(colour = 'black',size=15,face="bold"),
                   strip.text.x = element_text(colour = 'black',size=16))
    q <- q + theme(legend.text =element_text(colour = 'black',size=14),
                   legend.title = element_text(colour = 'black',size=14,face="bold"),
                   legend.position = "bottom",
                   legend.justification = "center")
    q <- q + guides(color=guide_legend(title=legend.titl),
                    fill=guide_legend(nrow=2,byrow=F))
    return(q)
    
  }
  
  pdf(paste0(OUTDIR,"PF_ScaledIntronMetagenes_RPICGsspecificity.pdf"),width = 9,height = 6)
  ## all genes normalized signal
  pl1 <- pl[-grep("matched",pl$category),]
  pl1 <- pl1[pl1$category!="RP genes",]
  pl1 <-pl1[pl1$quant=="abso" & pl1$condition=="Spn1-depleted",]
  myplot(pl1,ylab = "normalized signal",legend.titl = "Spn1 levels")+ggtitle("Spn1-depleted")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  
  pl1 <- pl[-grep("matched",pl$category),]
  pl1 <- pl1[pl1$category!="RP genes",]
  pl1 <-pl1[pl1$quant=="abso" & pl1$condition=="non-depleted",]
  myplot(pl1,ylab = "normalized signal",legend.titl = "Spn1 levels")+ggtitle("Spn1 non-depleted")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  
  ## all genes normalized signal
  pl1 <- pl[-grep("matched",pl$category),]
  pl1 <- pl1[pl1$category!="RP genes",]
  pl1 <-pl1[pl1$quant=="depl/nondepl" & pl1$condition=="Spn1-depleted",]
  myplot(pl1,ylab = "normalized signal",legend.titl = "Spn1 levels",hline = T)+ggtitle("Spn1-depleted over non-depleted")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  
  #################
  pl1 <-pl[pl$quant=="abso" & pl$condition=="non-depleted" & 
             pl$category%in%c("RP genes","length-matched nonRP genes","expn-matched nonRP genes","length+expn-matched nonRP genes"),]
  myplot(pl1,ylab = "normalized signal",legend.titl = "Spn1 levels")+ggtitle("Spn1 non-depleted")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  
  pl1 <-pl[pl$quant=="abso" & pl$condition=="Spn1-depleted" & 
             pl$category%in%c("RP genes","length-matched nonRP genes","expn-matched nonRP genes","length+expn-matched nonRP genes"),]
  myplot(pl1,ylab = "normalized signal",legend.titl = "Spn1 levels")+ggtitle("Spn1-depleted")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  
  pl1 <-pl[pl$quant=="depl/nondepl" & pl$condition=="Spn1-depleted" & 
             pl$category%in%c("RP genes","length-matched nonRP genes","expn-matched nonRP genes","length+expn-matched nonRP genes"),]
  myplot(pl1,ylab = "normalized signal",legend.titl = "Spn1 levels",hline = T)+ggtitle("Spn1-depleted over non-depleted")+theme(plot.title = element_text(color="black",size=16,vjust=0.5,hjust = 0.5,face = "bold"))
  
  dev.off()
  
}



###### (6) Metagene signals centered at TSS (-500bp,TSS,+4kb)
######   spikein normalized signal and H3 or Rpb1 normalized signals
if(F){
  pl <-c()
  gl <- c()
  
  prot <- "H3 (round1)"
  f = paste("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_TSS_",prot,".RData",sep="")
  load(f)
  normF=res[res$Factor==prot,]$normF
  l[[1]] <- l[[1]]*normF[1]
  l[[2]] <- l[[2]]*normF[2]
  l[[3]] <- l[[3]]*normF[3]
  l[[4]] <- l[[4]]*normF[4]
  for(i in 1:length(l)){
    #pl <- rbind(pl,cbind(metagene.signal.plot(l[[i]]),category=names(l)[i],round="R1"))
    gl <- rbind(gl,cbind(metagene.signal.plot(l[[i]]),category=names(l)[i],round="R1"))
  }
  
  prot <- "H3"
  f = paste("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_TSS_",prot,".RData",sep="")
  load(f)
  normF=res[res$Factor==prot,]$normF
  l[[1]] <- l[[1]]*normF[1]
  l[[2]] <- l[[2]]*normF[2]
  l[[3]] <- l[[3]]*normF[3]
  l[[4]] <- l[[4]]*normF[4]
  for(i in 1:length(l)){
    #pl <- rbind(pl,cbind(metagene.signal.plot(l[[i]]),category=names(l)[i],round="R2"))
    gl <- rbind(gl,cbind(metagene.signal.plot(l[[i]]),category=names(l)[i],round="R2"))
  }
  
  myplot <- function(pl){
    p <- ggplot(pl, aes(x=pos,y=signal,ymax=signal+se, ymin=signal-se,col=category))+
      facet_wrap(~round,nrow = 2,scales = "free")+
      geom_vline(xintercept = c(50),colour="gray50", linetype = "dashed") +
      geom_linerange(col='gray') + geom_line(size=1.2) +
      scale_x_continuous(breaks =c(1,50,150,250,350,450), labels=c("","TSS","1kb","2kb","3kb","4kb"),expand = c(0,0)) +ylab("normalized signal") + xlab("") +
      theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = NA),axis.line = element_line(size = 0.7, colour = "black"))+
      theme(axis.text.x = element_text(colour = 'black',size=15),
            axis.text.y = element_text(colour = 'black',size=15), 
            axis.title.y = element_text(colour = 'black',size=18),
            strip.text.x = element_text(colour = 'black',size=16),
            strip.background = element_blank())+
      theme(legend.text =element_text(colour = 'black',size=14),
            legend.title = element_text(colour = 'black',size=15))+
      guides(color=guide_legend(title="Spn1 levels"))+
      theme(legend.text = element_text(size = 14,colour = "black"),legend.title = element_blank(),
            legend.key = element_blank(),legend.background = element_blank(),
            legend.position = c(1,0.2), 
            legend.justification = c(1, 1))
    return(p)
  }
  p <- myplot(pl) +ggtitle("After Spikein normalization")+theme(plot.title = element_text(size=18,color="black",hjust = 0.5))
  q <- myplot(gl)+ggtitle("Before Spikein normalization")+theme(plot.title = element_text(size=18,color="black",hjust = 0.5))
  
  pdf(paste0(OUTDIR,"H3Metagenes_Tests.pdf"),width = 8,height = 8)
  grid.arrange(p,q,ncol=2)
  dev.off()
}

## TSS centered
if(T){
  pl <- c()
  ## spike normalized
  prot <- "H3 (round1)"
  protein<- c("H3 (round1)","H3","H3K4me3","H3K36me2","H3K36me3","RNAPII (round1)","RNAPII","S2P","S5P","Set2","Spn1","Spt6")
  
  for( prot in protein){
    f = paste("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_TSS_",prot,".RData",sep="")
    if(file.exists(f)){
      cat(prot,"\n")
      cat(res[res$Factor==prot,]$normF,"\n")
      load(f)
      normF=res[res$Factor==prot,]$normF
      if(prot %in% c("RNAPII (round1)","RNAPII","S2P","S5P")){
        normF=c(1,1,1,1)
        cat ("  using unit normalization factors\n")
      }
      
      l[[1]] <- l[[1]]*normF[1]
      l[[2]] <- l[[2]]*normF[2]
      l[[3]] <- l[[3]]*normF[3]
      l[[4]] <- l[[4]]*normF[4]
      pl <- rbind(pl, 
                  cbind(mymetagenefun_tss(l,genes,facet=NULL),factor=paste0(prot),category="spikenorm"))
      rm(l,f,prot,normF)
    }else{
      cat("file corresponding to ",prot, " does not exists! \n")
    }
  }
  ## H3 orRpb1 normalized
  protein <- list()
  protein[["S2P"]] <- "RNAPII"
  protein[["S5P"]] <- "RNAPII"
  protein[["H3K4me3"]] <- "H3"
  protein[["H3K36me2"]] <- "H3"
  protein[["H3K36me3"]] <- "H3 (round1)"
  #protein[["Spt6"]] <- "RNAPII (round1)"
  for(i in  1:length(protein)){
    f = paste("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_TSS_",names(protein)[i],".RData",sep="")
    fc = paste("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_TSS_",protein[[i]],".RData",sep="")
    
    if(file.exists(f) & file.exists(fc)){
      cat(names(protein)[i],"\n")
      cat(res[res$Factor==names(protein)[i],]$normF,"\n")
      cat(protein[[i]],"\n")
      cat(res[res$Factor==protein[[i]],]$normF,"\n")
      
      load(f)
      l1 <- l; rm(l);
      load(fc)
      l2 <- l; rm(l);
      
      if(names(protein)[i]%in%c("S2P","S5P")){
        cat("using unit size factors\n")
        normF1=c(1,1,1,1)
        normF2=c(1,1,1,1)
      }else{
        normF1=res[res$Factor==names(protein)[[i]],]$normF
        normF2=res[res$Factor==protein[[i]],]$normF
      }
      
      l <-my_relative_signal(l1,l2,normF1,normF2)
      pl <- rbind(pl, 
                  cbind(mymetagenefun_tss(l,genes,facet=NULL),factor=paste0(names(protein)[i],"/",protein[[i]]),category="relative"))
      rm(l1,l2)
    }else{
      cat("files do not exists\n")
    }
  }
  
  pl$condition <- gsub("Depleted","Spn1-depleted",pl$condition)
  pl$condition <- gsub("Non-depleted","non-depleted",pl$condition)
  pl$factor <- gsub("RNAPII","Rpb1",pl$factor)
  pl$factor <- paste0(pl$factor," (R2)")
  pl$factor <- gsub("[(]round1[)] [(]R2[)]","(R1)",pl$factor)
  pl$factor <- gsub("Spt6 [(]R2[)]","Spt6 (R1)",pl$factor)
  pl$factor <- gsub("Set2 [(]R2[)]","Set2 (R1)",pl$factor)
  pl$factor <- gsub("H3K36me3 [(]R2[)]","H3K36me3 (R1)",pl$factor)
  pl$factor <- gsub("S2P","Ser2-P",pl$factor)
  pl$factor <- gsub("S5P","Ser5-P",pl$factor)
  
  
  gpl <- ggplot(pl, aes(x=pos,y=signal,ymax=signal+se, ymin=signal-se,col=condition))+
    facet_wrap(~factor,nrow = 5,scales = "free")+
    geom_vline(xintercept = c(50),colour="gray50", linetype = "dashed") +
    geom_linerange(col='gray') + geom_line(size=1.2) +
    scale_color_manual(values = cols) +scale_y_continuous(expand = c(0,0),limits = c(0,NA))+
    scale_x_continuous(breaks =c(1,50,150,250,350,450), labels=c("","TSS","1kb","2kb","3kb","4kb"),expand = c(0,0)) +ylab("normalized signal") + xlab("") +
    theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = NA),axis.line = element_line(size = 0.7, colour = "black"))+
    theme(axis.text.x = element_text(colour = 'black',size=15),
          axis.text.y = element_text(colour = 'black',size=15), 
          axis.title.y = element_text(colour = 'black',size=18),
          strip.text.x = element_text(colour = 'black',size=16),
          strip.background = element_blank())+
    theme(legend.text =element_text(colour = 'black',size=14),
          legend.title = element_text(colour = 'black',size=15))+
    guides(color=guide_legend(title="Spn1 levels"))+
    theme(legend.text = element_text(size = 14,colour = "black"),legend.title = element_blank(),
          legend.key = element_blank(),legend.background = element_blank(),
          legend.position = c(1,0.1), 
          legend.justification = c(1, 1))
  
  pdf(file=paste0(OUTDIR,"MetageneSignals centerd atTSS.pdf"),width = 15,height = 12)
  gpl
  dev.off()
  save(pl,file="A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/MergedMetagenes_ChIPSpike_final_unscaled.RData")
  
}

## CPS centered
if(T){
  pl <- c()
  ## spike normalized
  prot <- "H3 (round1)"
  protein<- c("H3 (round1)","H3","H3K4me3","H3K36me2","H3K36me3","RNAPII (round1)","RNAPII","S2P","S5P","Set2","Spn1","Spt6")
  
  for( prot in protein){
    f = paste("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_CPS_",prot,".RData",sep="")
    if(file.exists(f)){
      cat(prot,"\n")
      cat(res[res$Factor==prot,]$normF,"\n")
      load(f)
      normF=res[res$Factor==prot,]$normF
      if(prot %in% c("RNAPII (round1)","RNAPII","S2P","S5P")){
        normF=c(1,1,1,1)
        cat ("  using unit normalization factors\n")
      }
      l[[1]] <- l[[1]]*normF[1]
      l[[2]] <- l[[2]]*normF[2]
      l[[3]] <- l[[3]]*normF[3]
      l[[4]] <- l[[4]]*normF[4]
      pl <- rbind(pl, 
                  cbind(mymetagenefun_tss(l,genes,facet=NULL),factor=paste0(prot),category="spikenorm"))
      rm(l,f,prot,normF)
    }else{
      cat("file corresponding to ",prot, " does not exists! \n")
    }
  }
  
  ## H3 orRpb1 normalized
  protein <- list()
  protein[["S2P"]] <- "RNAPII"
  protein[["S5P"]] <- "RNAPII"
  protein[["H3K4me3"]] <- "H3"
  protein[["H3K36me2"]] <- "H3"
  protein[["H3K36me3"]] <- "H3 (round1)"
  #protein[["Spt6"]] <- "RNAPII (round1)"
  for(i in  1:length(protein)){
    f = paste("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_CPS_",names(protein)[i],".RData",sep="")
    fc = paste("A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Revised_MultigeneHeatmap_CPS_",protein[[i]],".RData",sep="")
    
    if(file.exists(f) & file.exists(fc)){
      cat(names(protein)[i],"\n")
      cat(res[res$Factor==names(protein)[i],]$normF,"\n")
      cat(protein[[i]],"\n")
      cat(res[res$Factor==protein[[i]],]$normF,"\n")
      
      load(f)
      l1 <- l; rm(l);
      load(fc)
      l2 <- l; rm(l);
      
      if(names(protein)[i]%in%c("S2P","S5P")){
        cat("using unit size factors\n")
        normF1=c(1,1,1,1)
        normF2=c(1,1,1,1)
      }else{
        normF1=res[res$Factor==names(protein)[[i]],]$normF
        normF2=res[res$Factor==protein[[i]],]$normF
      }
      
      l <-my_relative_signal(l1,l2,normF1,normF2)
      pl <- rbind(pl, 
                  cbind(mymetagenefun_tss(l,genes,facet=NULL),factor=paste0(names(protein)[i],"/",protein[[i]]),category="relative"))
      rm(l1,l2)
    }else{
      cat("files do not exists\n")
    }
  }
  
  pl$condition <- gsub("Depleted","Spn1-depleted",pl$condition)
  pl$condition <- gsub("Non-depleted","non-depleted",pl$condition)
  pl$factor <- gsub("RNAPII","Rpb1",pl$factor)
  pl$factor <- paste0(pl$factor," (R2)")
  pl$factor <- gsub("[(]round1[)] [(]R2[)]","(R1)",pl$factor)
  pl$factor <- gsub("Spt6 [(]R2[)]","Spt6 (R1)",pl$factor)
  pl$factor <- gsub("Set2 [(]R2[)]","Set2 (R1)",pl$factor)
  pl$factor <- gsub("H3K36me3 [(]R2[)]","H3K36me3 (R1)",pl$factor)
  pl$factor <- gsub("S2P","Ser2-P",pl$factor)
  pl$factor <- gsub("S5P","Ser5-P",pl$factor)
  
  
  gpl <- ggplot(pl, aes(x=pos,y=signal,ymax=signal+se, ymin=signal-se,col=condition))+
    facet_wrap(~factor,nrow = 5,scales = "free")+
    geom_vline(xintercept = c(400),colour="gray50", linetype = "dashed") +
    geom_linerange(col='gray') + geom_line(size=1.2) +
    scale_color_manual(values = cols) +scale_y_continuous(expand = c(0,0),limits = c(0,NA))+
    scale_x_continuous(breaks =c(1,100,200,300,400,450), labels=c("-4kb","-3kb","-2kb","-1kb","CPS",""),expand = c(0,0)) +ylab("normalized signal") + xlab("") +
    theme(panel.background = element_rect(fill = NA),panel.grid.major = element_line(colour = NA),axis.line = element_line(size = 0.7, colour = "black"))+
    theme(axis.text.x = element_text(colour = 'black',size=15),
          axis.text.y = element_text(colour = 'black',size=15), 
          axis.title.y = element_text(colour = 'black',size=18),
          strip.text.x = element_text(colour = 'black',size=16),
          strip.background = element_blank())+
    theme(legend.text =element_text(colour = 'black',size=14),
          legend.title = element_text(colour = 'black',size=15))+
    guides(color=guide_legend(title="Spn1 levels"))+
    theme(legend.text = element_text(size = 14,colour = "black"),legend.title = element_blank(),
          legend.key = element_blank(),legend.background = element_blank(),
          legend.position = c(1,0.1), 
          legend.justification = c(1, 1))
  
  pdf(file=paste0(OUTDIR,"MetageneSignals centerd atCPS.pdf"),width = 15,height = 12)
  gpl
  dev.off()
  #save(pl,file="A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/MergedMetagenes_ChIPSpike_final_unscaled.RData")
}






