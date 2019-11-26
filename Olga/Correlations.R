
### get pl dataframe
if(F){
  rm(list=ls())
  chr.flt <- c("chrI","chrII","chrIII","chrIV","chrIX","chrV","chrVI","chrVII","chrVIII","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI")
  spiked.chr <- c("chrXVII","chrXVIII","chrXIX")
  
  ### TSS-seq
  load(file="A:/work/WinstonLab/Olga/TSS_Seq/Steinmetz/DE_promoters_PostSpikein_L2FC.RData")
  m.occup <- data.frame(tracking_id=l[[1]]$tracking_id)
  for(i in 1:length(l)){
    cat(names(l)[i],"\n")
    cat(" ",paste0(names(l[[i]][,5:8]),collapse = ","),"\n")
    qwe <- l[[i]][,5:8]
    qwe$x <- apply(qwe[,1:2],1,function(x) sqrt(prod(x)))
    qwe$y <- apply(qwe[,3:4],1,function(x) sqrt(prod(x)))
    names(qwe)[5:6] <- names(qwe)[c(1,3)]
    qwe <- qwe[,5:6]
    names(qwe) <- gsub("_[01234567]$","",names(qwe))
    m.occup <- cbind(m.occup,qwe)
    rm(qwe)
    cat("\n\n")
  }
  m <- data.frame(tracking_id=l[[1]]$tracking_id)
  for(i in 1:length(l)){
    m <- cbind(m,l[[i]][,3])
  }
  m.occup <- m.occup[,unique(names(m.occup))]
  names(m)[2:ncol(m)] <- names(l)
  names(m)[2:ncol(m)] <- paste0("tssl2fc_",names(m)[2:ncol(m)])
  names(m.occup)[2:ncol(m.occup)] <- paste0("tssOccup_",names(m.occup)[2:ncol(m.occup)])
  z <- merge(m,m.occup)
  rm(l,m,m.occup)
  
  ## NETSeq
  load("A:/work/WinstonLab/Olga/Netseq/DE_promoters_PostSpikein_L2FC.RData")
  m.occup <- data.frame(tracking_id=l[[1]]$tracking_id)
  for(i in 1:length(l)){
    cat(names(l)[i],"\n")
    cat(" ",paste0(names(l[[i]][,5:8]),collapse = ","),"\n")
    qwe <- l[[i]][,5:8]
    qwe$x <- apply(qwe[,1:2],1,function(x) sqrt(prod(x)))
    qwe$y <- apply(qwe[,3:4],1,function(x) sqrt(prod(x)))
    names(qwe)[5:6] <- names(qwe)[c(1,3)]
    qwe <- qwe[,5:6]
    names(qwe) <- gsub("_[01234567]$","",names(qwe))
    names(qwe) <- gsub("_I|_II","",names(qwe))
    m.occup <- cbind(m.occup,qwe)
    rm(qwe)
    cat("\n\n")
  }
  m <- data.frame(tracking_id=l[[1]]$tracking_id)
  for(i in 1:length(l)){
    m <- cbind(m,l[[i]][,3])
  }
  names(m.occup) <- gsub("netseq_","",names(m.occup))
  m.occup <- m.occup[,unique(names(m.occup))]
  names(m)[2:ncol(m)] <- names(l)
  names(m)[2:ncol(m)] <- paste0("netl2fc_",names(m)[2:ncol(m)])
  names(m.occup)[2:ncol(m.occup)] <- paste0("netOccup_",names(m.occup)[2:ncol(m.occup)])
  z1 <- merge(m,m.occup)
  rm(l,m,m.occup)
  z <- merge(z,z1,by="tracking_id")
  rm(z1)
  
  ### ChIP-seq
  load(file="A:/work/WinstonLab/Olga/ChIPSeq/DESeq2_ChIPSeq_L2FC_DFs.RData")
  m.occup <- m.occup[,unique(names(m.occup))]
  names(m)[2:25] <- paste0("chipl2fc_",names(m)[2:25])
  names(m.occup)[2:49] <- paste0("chipOccup_",names(m.occup)[2:49])
  pl <- merge(m,z,by="tracking_id")
  pl <- merge(pl,m.occup,by="tracking_id")
  rm(z,m,m.spike,m.occup,m.spike.occup)
  pl <- pl[,-grep("_pob3|_double",names(pl))]
  pl <- pl[,-grep("spt6-1004",names(pl))]
  pl <- pl[,-grep("Spt6_1004",names(pl))]
}

### ribbon plots for l2fc with nd without spikein normalization
if(F){
  rm(list=ls())
  load(file="A:/work/WinstonLab/Olga/ChIPSeq/Diff_Factor_Occupancy_DESeq2SGD.RData")
  m.spike <- data.frame(tracking_id=l[[1]]$tracking_id)
  m <- data.frame(tracking_id=l[[1]]$tracking_id)
  for(i in 1:length(ln)){
    m.spike <- cbind(m.spike,l[[i]][,3])
    m <- cbind(m,w[[i]][,3])
  }
  names(m)[2:25] <- names(m.spike)[2:25] <- names(l)
  
  m.occup <- data.frame(tracking_id=w[[1]]$tracking_id)
  m.spike.occup <- data.frame(tracking_id=l[[1]]$tracking_id)
  for(i in 1:length(ln)){
    cat(names(ln)[i],"\n")
    cat(" ",paste0(names(l[[i]][,5:8]),collapse = ","),"\n")
    qwe <- l[[i]][,5:8]
    qwe$x <- log2(qwe[,3]/qwe[,1])
    qwe$y <- log2(qwe[,4]/qwe[,2])
    names(qwe)[5:6] <- names(qwe)[c(3,4)]
    qwe <- qwe[,5:6]
    m.spike.occup <- cbind(m.spike.occup,qwe)
    rm(qwe)
    
    qwe <- w[[i]][,5:8]
    cat(" ",paste0(names(w[[i]][,5:8]),collapse = ","),"\n")
    qwe$x <- log2(qwe[,3]/qwe[,1])
    qwe$y <- log2(qwe[,4]/qwe[,2])
    names(qwe)[5:6] <- names(qwe)[c(3,4)]
    qwe <- qwe[,5:6]
    m.occup <- cbind(m.occup,qwe)
    rm(qwe)
    cat("\n\n")
  }
  rm(l,ln,w,wn)
  
  
  l <- list()
  for(gt in c("pob3","spt6-YW","double")){
    for(chip in c("Pol2-July","Pol2-Aug","Spt6","Spt16")){
      for(tmp in c(30,37)){
        cat(gt,"\t",chip,"\t",tmp,"\n")
        pl <- m.occup
        z <- pl[,grepl(gt,names(pl),fixed = T) & grepl(chip,names(pl),fixed = T) & grepl(tmp,names(pl),fixed = T)]
        pl <- m
        z$avg <- pl[,grepl(gt,names(pl),fixed = T) & grepl(chip,names(pl),fixed = T) & grepl(tmp,names(pl),fixed = T)]
        names(z) <- c("x","y","avg")
        z$se1 <- apply(z[,1:2],1,function(x) min(x,na.rm=T))
        z$se2 <- apply(z[,1:2],1,function(x) max(x,na.rm=T))
        
        pl <- m.spike.occup
        z1 <- pl[,grepl(gt,names(pl),fixed = T) & grepl(chip,names(pl),fixed = T) & grepl(tmp,names(pl),fixed = T)]
        pl <- m.spike
        z1$avg <- pl[,grepl(gt,names(pl),fixed = T) & grepl(chip,names(pl),fixed = T) & grepl(tmp,names(pl),fixed = T)]
        names(z1) <- c("x","y","avg")
        z1$se1 <- apply(z1[,1:2],1,function(x) min(x,na.rm=T))
        z1$se2 <- apply(z1[,1:2],1,function(x) max(x,na.rm=T))
        
        z$cond="LibSize"
        z1$cond="Spikein"
        z$id = 1:nrow(z)
        z1$id = 1:nrow(z1)
        z <- z[order(z$avg,decreasing = F),]
        z1 <- z1[match(z$id,z1$id),]
        z$id = 1:nrow(z)
        z1$id = 1:nrow(z1)
        z <- rbind(z[,3:7],z1[,3:7])
        head(z)
        
        q <- ggplot(z,aes(x=id,y=avg,col=cond))+
          geom_linerange(alpha=0.05,aes(ymin=se1, ymax=se2,col=cond)) + theme_bw()+ylim(-2,2)+geom_point(size=0.1,alpha=0.5)+
          scale_color_brewer(palette = "Dark2")+
          xlab("5134 genes")+ylab("log2 fold change values")+
          ggtitle(paste(gt,chip,tmp,sep="_"))+
          theme(axis.title = element_text(colour = "black",size = 14),
                plot.title = element_text(hjust = 0.5,colour = "black",size = 18))+
          geom_hline(yintercept = 0,col='black',size=0.5)+
          theme(legend.position = c(0.5, 0.9),
                legend.title = element_blank())
        l[[paste(gt,chip,tmp,sep="_")]] <- print(q)
      }
    }
  }
  rm(q,z,z1,pl)  
  
  
  OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/"
  png(file = paste0(OUTDIR,"Spombe_l2fc-comparisn under Spikein LibSize normalization.png"),width = 800,height = 1200)
  grid.arrange(grobs=l,ncol=4,top= textGrob(paste0("Log2FC comparison under LibSize and Spikein Normalization"),gp=gpar(fontsize=20,font=2)))
  dev.off()  
  
  
  
  }

### plot correlation scatter plots
if(T){
  OUTDIR="A:/work/WinstonLab/Olga/ChIPSeq/"
  gg <- theme_minimal() +theme(axis.text.x = element_text(size=14,colour = "black"),
                               axis.text.y = element_text(size=14,colour = "black"),
                               axis.title =  element_text(size=16,color="black"),
                               axis.line = element_line(color="gray",size = 0.75),
                               axis.ticks = element_line(color="gray",size = 1),
                               legend.position = "none",
                               plot.title = element_text(colour ="black",size = 18,hjust = 0.5))
  
  require(ggpubr)
  myplotfun <- function(l2fc,xlab,ylab,main,x,y){
      q <-ggplot(l2fc, aes(x=l2fc[,x], y=l2fc[,y])) +
        xlab(paste0(xlab))+ylab(paste0(ylab))+ggtitle(paste(main))+
        stat_bin_hex(geom="point", aes(color=log2(..count..)), binwidth=c(.05,.05), size=2, shape=16, stroke=0)+
        scale_color_viridis(option="inferno",direction = 1)+
        geom_hline(yintercept = 0,col="black")+geom_vline(xintercept = 0,col="black")+
        geom_smooth(method='lm',formula=y~x,col="red")+
        stat_cor(method = "pearson",size=5,col="red",label.sep = "\n",
                 label.x.npc = "right", label.y.npc = "bottom")+gg+
        scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))
      return(q)
    }
  
  ### (1).  with respect to tss change,sense, 30C  
  bs <- "tss_spt6-YW_30"
  qs <- names(pl)[grepl("chip_",names(pl),fixed = T) & grepl("30",names(pl),fixed = T)]
  pl.list <- list()
  for(i in 1:length(qs)){
    xlab="TSS-seq, log2FC"
    ylab="ChIP-seq, log2FC"
    x=as.character(bs)
    y=as.character(qs[i])
    main=gsub("chip_","",as.character(qs[i]))
    cat(main,"\n")
    pl.list[[paste0(main)]] <- print(myplotfun(l2fc = pl,xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
    rm(xlab,ylab,main,y,x)
    }
  pdf(file = paste0(OUTDIR,"senseTSS-seq vs ChIPSeq_30C.pdf"),width = 16,height = 12)
  grid.arrange(grobs=pl.list,ncol=4,top= textGrob(paste0("occupancy change v/s TSS-seq log2FC, 30°C"),gp=gpar(fontsize=20,font=2)))
  dev.off()  
  
  ### (2).  with respect to tss change,sense, 37C  
  bs <- "tss_spt6-YW_37"
  qs <- names(pl)[grepl("chip_",names(pl),fixed = T) & grepl("37",names(pl),fixed = T)]
  pl.list <- list()
  for(i in 1:length(qs)){
    xlab="TSS-seq, log2FC"
    ylab="ChIP-seq, log2FC"
    x=as.character(bs)
    y=as.character(qs[i])
    main=gsub("chip_","",as.character(qs[i]))
    cat(main,"\n")
    pl.list[[paste0(main)]] <- print(myplotfun(l2fc = pl,xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
    rm(xlab,ylab,main,y,x)
  }
  pdf(file = paste0(OUTDIR,"senseTSS-seq vs ChIPSeq_37C.pdf"),width = 16,height = 12)
  grid.arrange(grobs=pl.list,ncol=4,top= textGrob(paste0("occupancy change v/s TSS-seq log2FC, 37°C"),gp=gpar(fontsize=20,font=2)))
  dev.off()  
  
  ### (3).  with respect to netseq change,sense, 30C  
  bs <- "net_spt6-YW_30"
  qs <- names(pl)[grepl("chip_",names(pl),fixed = T) & grepl("30",names(pl),fixed = T)]
  pl.list <- list()
  for(i in 1:length(qs)){
    xlab="NET-seq, log2FC"
    ylab="ChIP-seq, log2FC"
    x=as.character(bs)
    y=as.character(qs[i])
    main=gsub("chip_","",as.character(qs[i]))
    cat(main,"\n")
    pl.list[[paste0(main)]] <- print(myplotfun(l2fc = pl,xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
    rm(xlab,ylab,main,y,x)
  }
  pdf(file = paste0(OUTDIR,"senseNET-seq vs ChIPSeq_30C.pdf"),width = 16,height = 12)
  grid.arrange(grobs=pl.list,ncol=4,top= textGrob(paste0("occupancy change v/s NET-seq log2FC, 30°C"),gp=gpar(fontsize=20,font=2)))
  dev.off() 
  
  ### (4).  with respect to netseq change,sense, 30C  
  bs <- "net_spt6-YW_37"
  qs <- names(pl)[grepl("chip_",names(pl),fixed = T) & grepl("37",names(pl),fixed = T)]
  pl.list <- list()
  for(i in 1:length(qs)){
    xlab="NET-seq, log2FC"
    ylab="ChIP-seq, log2FC"
    x=as.character(bs)
    y=as.character(qs[i])
    main=gsub("chip_","",as.character(qs[i]))
    cat(main,"\n")
    pl.list[[paste0(main)]] <- print(myplotfun(l2fc = pl,xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
    rm(xlab,ylab,main,y,x)
  }
  pdf(file = paste0(OUTDIR,"senseNET-seq vs ChIPSeq_37C.pdf"),width = 16,height = 12)
  grid.arrange(grobs=pl.list,ncol=4,top= textGrob(paste0("occupancy change v/s NET-seq log2FC, 37°C"),gp=gpar(fontsize=20,font=2)))
  dev.off() 
  
  ### (5).  with respect to tss change,sense, 30C  
  bs <- "tssAS_spt6-YW_30"
  qs <- names(pl)[grepl("chip_",names(pl),fixed = T) & grepl("30",names(pl),fixed = T)]
  pl.list <- list()
  for(i in 1:length(qs)){
    xlab="AS TSS-seq, log2FC"
    ylab="ChIP-seq, log2FC"
    x=as.character(bs)
    y=as.character(qs[i])
    main=gsub("chip_","",as.character(qs[i]))
    cat(main,"\n")
    pl.list[[paste0(main)]] <- print(myplotfun(l2fc = pl,xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
    rm(xlab,ylab,main,y,x)
  }
  pdf(file = paste0(OUTDIR,"AntisenseTSS-seq vs ChIPSeq_30C.pdf"),width = 16,height = 12)
  grid.arrange(grobs=pl.list,ncol=4,top= textGrob(paste0("occupancy change v/s Antisense TSS-seq log2FC, 30°C"),gp=gpar(fontsize=20,font=2)))
  dev.off()  
  
  ### (6).  with respect to tss change,sense, 37C  
  bs <- "tssAS_spt6-YW_37"
  qs <- names(pl)[grepl("chip_",names(pl),fixed = T) & grepl("37",names(pl),fixed = T)]
  pl.list <- list()
  for(i in 1:length(qs)){
    xlab="AS TSS-seq, log2FC"
    ylab="ChIP-seq, log2FC"
    x=as.character(bs)
    y=as.character(qs[i])
    main=gsub("chip_","",as.character(qs[i]))
    cat(main,"\n")
    pl.list[[paste0(main)]] <- print(myplotfun(l2fc = pl,xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
    rm(xlab,ylab,main,y,x)
  }
  pdf(file = paste0(OUTDIR,"AntisenseTSS-seq vs ChIPSeq_37C.pdf"),width = 16,height = 12)
  grid.arrange(grobs=pl.list,ncol=4,top= textGrob(paste0("occupancy change v/s Antisense TSS-seq log2FC, 37°C"),gp=gpar(fontsize=20,font=2)))
  dev.off()  
  
  ### (7).  with respect to netseq change,sense, 30C  
  bs <- "netAS_spt6-YW_30"
  qs <- names(pl)[grepl("chip_",names(pl),fixed = T) & grepl("30",names(pl),fixed = T)]
  pl.list <- list()
  for(i in 1:length(qs)){
    xlab="AS NET-seq, log2FC"
    ylab="ChIP-seq, log2FC"
    x=as.character(bs)
    y=as.character(qs[i])
    main=gsub("chip_","",as.character(qs[i]))
    cat(main,"\n")
    pl.list[[paste0(main)]] <- print(myplotfun(l2fc = pl,xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
    rm(xlab,ylab,main,y,x)
  }
  pdf(file = paste0(OUTDIR,"AntisenseNET-seq vs ChIPSeq_30C.pdf"),width = 16,height = 12)
  grid.arrange(grobs=pl.list,ncol=4,top= textGrob(paste0("occupancy change v/s Antisense NET-seq log2FC, 30°C"),gp=gpar(fontsize=20,font=2)))
  dev.off() 
  
  ### (8).  with respect to netseq change,sense, 30C  
  bs <- "netAS_spt6-YW_37"
  qs <- names(pl)[grepl("chip_",names(pl),fixed = T) & grepl("37",names(pl),fixed = T)]
  pl.list <- list()
  for(i in 1:length(qs)){
    xlab="AS NET-seq, log2FC"
    ylab="ChIP-seq, log2FC"
    x=as.character(bs)
    y=as.character(qs[i])
    main=gsub("chip_","",as.character(qs[i]))
    cat(main,"\n")
    pl.list[[paste0(main)]] <- print(myplotfun(l2fc = pl,xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
    rm(xlab,ylab,main,y,x)
  }
  pdf(file = paste0(OUTDIR,"AntisenseNET-seq vs ChIPSeq_37C.pdf"),width = 16,height = 12)
  grid.arrange(grobs=pl.list,ncol=4,top= textGrob(paste0("occupancy change v/s Antisense NET-seq log2FC, 37°C"),gp=gpar(fontsize=20,font=2)))
  dev.off()
  
  
  
  
  #### (9) TSS-seq vs NET-seq
  pdf(file = paste0(OUTDIR,"NET-Seq vs TSS-seq log2FC comparisons.pdf"),width = 5,height = 5)
  myplotfun(l2fc = pl,ylab="NET-seq, log2FC",xlab='TSS-seq, log2FC',main="spt6-YW 30°C",
            y="net_spt6-YW_30",x="tss_spt6-YW_30" )
  myplotfun(l2fc = pl,ylab="NET-seq, log2FC",xlab='TSS-seq, log2FC',main="spt6-YW 37°C",
            y="net_spt6-YW_37",x="tss_spt6-YW_37" )
  
  myplotfun(l2fc = pl,ylab="NET-seq, log2FC",xlab='TSS-seq, log2FC',main="spt6-YW 30°C (Antisense)",
            y="netAS_spt6-YW_30",x="tssAS_spt6-YW_30" )
  myplotfun(l2fc = pl,ylab="NET-seq, log2FC",xlab='TSS-seq, log2FC',main="spt6-YW 37°C (Antisense)",
            y="netAS_spt6-YW_37",x="tssAS_spt6-YW_37" )
  dev.off()
  
  #### 
}







