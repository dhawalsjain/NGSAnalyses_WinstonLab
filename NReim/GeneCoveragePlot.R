#rm(list=ls())
library(ggplot2)
library(grid)
library(gridExtra)
library(GenomicRanges)
library(dplyr)
gene.coverage.plot.df <- function(path="A:/work/WinstonLab/Natalia_Reim/October2016/BedGraph/SPN1_DMSO_Mergeq50WR.bedGraph",
                                profile="SPN1-AID_DMSO",genes,Margin=500,glist){
  #load(path)
  #gr$score <- gr$score*(-1) ## gr is derived directly from the bam file, where read strand info is reversed
  #grn <- gr[gr$score<0,]
  #grp <- gr[gr$score>0,]
  #grp <- coverage(grp,weight = grp$score)
  #grn$score <- grn$score*(-1)
  #grn <- coverage(grn,weight = grn$score)
  
  gr <- read.delim(paste0(path),header=F)
  gr <- with(gr,GRanges(V1,IRanges(V2,V3),"+",score=V4))
  grn <- gr[gr$score<0,]
  grp <- gr[gr$score>0,]
  grp <- coverage(grp,weight = grp$score)
  grn$score <- grn$score*(-1)
  grn <- coverage(grn,weight = grn$score)
  
  
  get.plot.df <- function(gene,Margin=200){
    g <- genes[genes$tracking_id==gene,]
    if(length(g$tracking_id)==0){
      cat("Gene:", gene,"Error\n")
      return(NULL)
    }
    if(as.character(g$strand)=="+"){
      s <- as.vector(window(grp[[as.character(g$chr)]],g$start-Margin, g$end+Margin))
      as <- as.vector(window(grn[[as.character(g$chr)]],g$start-Margin, g$end+Margin))
      
    }else{
      s <- (as.vector(window(grn[[as.character(g$chr)]],g$start-Margin, g$end+Margin)))
      as <- (as.vector(window(grp[[as.character(g$chr)]],g$start-Margin, g$end+Margin)))
      
    }
    pl <- rbind(data.frame(pos=(g$start-Margin):(g$end+Margin),signal=s,ori="sense"),
                data.frame(pos=(g$start-Margin):(g$end+Margin),signal=as*(-1),ori="antisense"))
    pl$tracking_id <- g$tracking_id
    pl$gene <- g$gene
    return(pl)
  }
  
  glist <- unlist(strsplit(glist,","))
  pll <- c()
  for(gene in glist){
    pll <- rbind(pll, cbind(get.plot.df(gene,Margin = Margin), profile=paste0(profile)))
  }
  return(pll)
  
}

gene.coverage.plot.ChIPSeq.df <- function(path="A:/work/WinstonLab/Natalia_Reim/MayNov17/vis/Depleted_RNAPII.merged.bedGraph",profile="SPN1-AID_DMSO",genes,Margin=500,glist){
  
  gr <- read.delim(paste0(path),header=F)
  gr <- with(gr,GRanges(V1,IRanges(V2,V3),"+",score=V4))
  cov <- coverage(gr,weight = gr$score)
  
  get.plot.df <- function(gene,Margin=200){
    g <- genes[genes$tracking_id==gene,]
    if(length(g$tracking_id)==0){
      cat("Gene:", gene,"Error\n")
      return(NULL)
    }
    if(as.character(g$strand)=="+"){
      s <- as.vector(window(cov[[as.character(g$chr)]],g$start-Margin, g$end+Margin))
    }else{
      s <- (as.vector(window(cov[[as.character(g$chr)]],g$start-Margin, g$end+Margin)))
    }
    pl <- data.frame(pos=(g$start-Margin):(g$end+Margin),signal=s,ori="sense")
    pl$tracking_id <- g$tracking_id
    pl$gene <- g$gene
    return(pl)
  }
  
  glist <- unlist(strsplit(glist,","))
  pll <- c()
  for(gene in glist){
    pll <- rbind(pll, cbind(get.plot.df(gene,Margin = Margin), profile=paste0(profile)))
  }
  return(pll)
}

alterGFF= function(gff=gff){
  ## this function adds columns to a getGFF() dataframe fro plotting by ggplot2
  ## this function requires specific named columns so is very customisable
  
  ## R does not like variables or column names that begin with numbers e.g. 5_utr_start
  
  #colnames(gff)[4:7]=c('utr5_start', 'utr5_end', 'utr3_start', 'utr3_end')
  
  ## categorise which rows are exons and which utrs and add a column
  type=rep('exon', dim(gff)[1])
  #type[which(!is.na(gff$utr5_start))]='utr5'
  #type[which(!is.na(gff$utr3_start))]='utr3'
  gff=data.frame(gff,type=type)
  
  
  ## to plot a line the length of each gene I need to add this data to each row
  ## so I need to find the max and min exons for each transcript
  ## the merge them back into each row of the original data
  ag.min=aggregate(gff, list(gff$tracking_id), function(x)min(as.numeric(x),na.rm=T))
  ag.max=aggregate(gff, list(gff$tracking_id), function(x)max(as.numeric(x),na.rm=T))
  ag=data.frame(tracking_id= ag.min$Group.1, start=ag.min$gene_start, end=ag.max$gene_end, strand=ag.max$strand)
  ag= ag[order(ag$start),]
  
  ## We are adding a y offset for each transcript element 
  ## So we have to see whether they are overlapping i.e. if so it needs a new row	##
  ## If it is on the pos or neg strand											##
  ## good luck to anyone who wants to try unlooping this							##
  
  yneg=0; ypos=0; mem=vector(); y=vector()
  for(i in 1:length(ag$start))
  {
    genevec=ag$start[i]:ag$end[i]
    ynew=ifelse(length(intersect(genevec, mem))>0, TRUE, FALSE)
    strandpos=ifelse(ag$strand[i]>0, TRUE, FALSE)
    
    if(ynew && strandpos)
      ypos=ypos+1
    if(ynew && !strandpos)
      yneg=yneg-1
    if(!ynew && strandpos)	
      ypos= 1
    if(!ynew && !strandpos)	
      yneg= -1	
    
    if(strandpos)
      y[i]=ypos
    if(!strandpos)
      y[i]=yneg
    
    mem=c(mem,genevec)
  }
  
  ag=data.frame(ag,y=y)
  
  ## here we are merging the start end and the y offset for each transcript back into the original data
  ## so it is associated with every feature-row
  gff=merge(ag,gff) 
  
  ## to plot the line we also need to show arrow direction which requires an option first or last
  ardir=ifelse(gff$strand==1,'last','first')
  gff=data.frame(gff,ardir=ardir)
  
  ## OK there is an annoying glitch- that if you use strand for facet_grid then -1 will be on top
  ## the easiest workround is to create a dummy char variable that is alphabetically the right way!!!!!
  
  strandChar=gff$strand
  strandChar[strandChar==1]<-'forward'
  strandChar[strandChar==(-1)]<-'reverse'
  gff.new=data.frame(gff, strandChar=strandChar)
  
  return(gff.new)
}

gene.coverage.plot <- function(pl,samples=c("WT_DMSO"),mygene="SER3",ori="",
                               cols=cols,anno=anno,genes=genes,
                               Margin=c(200,200),plot.coverage=F,only.gene=F,stranded=F,
                               log.tranform=F,ylim=c(0,Inf)){
  cat(mygene,"\n")
  #mygene=c("SER3","HSP82")
  if(length(Margin)==1){
    Margin = c(Margin,Margin)
  }
  #mygene="RPS27A"
  multfactor <- 1
  if(ori=="sense"){
    pl1 <- pl[pl$gene%in%mygene & pl$ori=="sense" ,]
  }else if(ori=="antisense"){
    pl1 <- pl[pl$gene%in%mygene & pl$ori=="antisense" ,]
    pl1$signal <- pl1$signal*(-1)
  }else{
    pl1 <- pl[pl$gene%in%mygene,]
  }
  pl1 <- pl1[pl1$profile%in%samples,]
  g <- genes[genes$gene==mygene,]
  g <- as(g,"GRanges")
  #g <- resize(g,width(g)+(2*Margin+1),fix="center")
  g <- resize(g,width(g)+(Margin[1]+1),fix="end")
  g <- resize(g,width(g)+Margin[2],fix="start")
  
  
  pl1$pos <- pl1$pos-start(g)
  pl2 <- pl1[pl1$pos>=0 & pl1$pos<=width(g),]
  if(ori =="" & max(pl1[pl1$ori=="sense",]$signal)>10){
    multfactor <- ceiling(max(pl2[pl2$ori=="sense",]$signal)/10)
    pl1[pl1$ori=="sense",]$signal <- pl1[pl1$ori=="sense",]$signal/multfactor
  }
  if(as.character(strand(g))=="-"){
    p <- ggplot(pl1,aes(x=rev(pos),y=signal,fill=colguide)) +geom_area()+ facet_wrap(~profile,ncol = 1)
    
  }else{
    p <- ggplot(pl1,aes(x=(pos),y=signal,fill=colguide)) +geom_area()+ facet_wrap(~profile,ncol = 1)
  }
  p <- p+ ylab("") #normalized counts
  p <- p + xlim(0,width(g))
  p <- p + scale_fill_manual(values=cols)+theme_minimal() 
  #p <- p + scale_x_continuous(expand = c(0, 0))
  p <- p + theme(axis.text.y = element_text(size=14,colour ="black"),
                 axis.title.y = element_text(size=16,colour="black"),
                 axis.title.x = element_blank(),
                 axis.text.x = element_text(size=14,colour ="black"))
  p <- p + theme(axis.ticks = element_line(size=1,colour = "black"))
  p <- p + theme(panel.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(size = 0.7, colour = "black"))
  p <- p + theme(legend.position = "none")+ theme(strip.background = element_blank(), 
                                                  strip.text = element_blank()) #element_text(color="black",size=18,hjust = 0)
  kb <- scales::unit_format(unit = "", scale = 1e-3)
  if(width(g)>=1000){
    breaks = seq(Margin[1],width(g)+Margin[1],by=1000)
    labels= kb((breaks-Margin[1]))
  }else{
    breaks = c(Margin[1],500+Margin[1])
    labels= kb((breaks-Margin[1]))
  }
  p <- p + scale_x_continuous(limits = c(0,width(g)),breaks=breaks,
                              labels = labels) #labels= scales::unit_format(unit = "", sep = "",scale = 1e-3)
  glength = paste0(round((width(g)- (Margin[1]+Margin[2] )  )/1000,1)," kb")
  
  #p <- p + scale_x_continuous(limits = c(0,width(g)),breaks = c(Margin,(width(g)-Margin)),
  #                            labels = c("0 kb",glength))
  p <- p + geom_hline(yintercept = 0,col="gray")
  if(ori ==""){
    y.labels <- round(c(0,seq(min(pl2$signal),max(pl2$signal),length.out = 5)))
    y.labels <- sort(y.labels)
    y.breaks <- ifelse(y.labels>0,round(y.labels/multfactor),y.labels)
    p <- p + scale_y_continuous(breaks=y.breaks,labels = y.labels)
  }
  p
  
  gff <- subsetByOverlaps(anno,g,ignore.strand=T)
  gff <- as.data.frame(gff)
  gff$strand <- ifelse(gff$strand=="+",1,-1)
  gff <- gff[,c(1,5:12)]
  gff.new <- alterGFF(gff)
  
  gff.new$arrow.start <- gff.new$start
  gff.new$arrow.end <- gff.new$end
  gff.new$arrow.start <- ifelse(gff.new$strand==-1, gff.new$end,gff.new$start)
  gff.new$arrow.end <- ifelse(gff.new$strand==-1, gff.new$start,gff.new$end)
  
  mystrand = unique(gff.new[gff.new$gene==mygene,]$strand)
  if(stranded==T){
    gff.new <- gff.new[gff.new$strand==mystrand,]
  }
  if(mygene=="SER3"){
    gff.new <- gff.new[gff.new$gene%in%c(mygene,"SRG1"),]
  }
  if(only.gene==T){
    gff.new <- gff.new[gff.new$gene==mygene,]
  }
  
  gff.new$exon_start <- ifelse(gff.new$exon_start<start(g),start(g),gff.new$exon_start)
  gff.new$exon_end <- ifelse(gff.new$exon_end>end(g),end(g),gff.new$exon_end)
  gff.new$start <- ifelse(gff.new$start<start(g),start(g),gff.new$start)
  gff.new$end <- ifelse(gff.new$end>end(g),end(g),gff.new$end)
  
  q <- ggplot(aes(xmin=exon_start, xmax=exon_end, ymin=y+0.4, ymax= y-0.4,col=strandChar,fill=strandChar, label=gene,col="black"), data=gff.new) 
  q <- q + geom_rect(color=NA)
  q <- q + scale_x_discrete(expand=c(0,0))
  if(as.character(strand(g))=="-"){
    q <- q + xlim(end(g),start(g))#
  }else{
    q <- q + xlim(start(g),end(g))#
  }
  q <- q + geom_segment(aes(x=arrow.start,xend=arrow.end, y=y, yend= y),color="gray60",arrow=arrow(),size=2,lineend = "round", linejoin = "bevel",inherit.aes = T) 
  q <- q + geom_label(aes(label=gene, x=(start + end)/2,y=y,fill=NA),color="black",size=5,fontface="italic")
   
  q <- q + scale_fill_manual(values=cols) + ylab("")
  #q <- q + scale_color_manual(values=c("reverse"="black","forward"="black")) 
  q <- q + theme(panel.background = element_blank(),
                                  panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                                  axis.line.y = element_blank(),
                                  axis.ticks = element_blank())
  q <- q + theme(axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.line.x = element_blank(),
                 axis.title.y = element_text(colour = "white"),
                 axis.text.y = element_text(colour ="white"))
  q <- q + theme(legend.position = "none")+ theme(strip.background = element_blank(), strip.text = element_blank()) 
  #x.labels <- c( start(g),(start(g)+Margin), (end(g)-Margin),end(g)) 
  q
  p <- ggplotGrob(p)
  q <- ggplotGrob(q)
  maxWidth = grid::unit.pmax(p$widths[2:5], q$widths[2:5])
  p$widths[2:5] <- as.list(maxWidth)
  q$widths[2:5] <- as.list(maxWidth)
  
  mylen = length(unique(gff.new$y))
  plot.list <- list()
  plot.list[[paste0(mygene,"_coverage_",mylen)]] <- p
  plot.list[[paste0(mygene,"_annotation_",mylen)]] <- q
  
  if(plot.coverage==T){
    if(length(unique(gff.new$y))==1){
      grid.arrange(p,q,ncol=1,heights=c(8,1.5)) #,top= textGrob(paste0(mygene),gp=gpar(fontsize=22,font=3))
    }else if(length(unique(gff.new$y))==2){
      grid.arrange(p,q,ncol=1,heights=c(7,2)) #,top= textGrob(paste0(mygene),gp=gpar(fontsize=22,font=3))
    }else{
      grid.arrange(p,q,ncol=1,heights=c(6,3)) #,top= textGrob(paste0(mygene),gp=gpar(fontsize=22,font=3))
    }
  }else{
    return(plot.list)
  }
} 







pdf("A:/work/WinstonLab/Natalia_Reim/ThesisFigures_March2019/STF1-RPP1B_ASplots.pdf",width = 8,height = 4)
gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = T,
                   mygene = "STF1",ori = "",only.gene = F,cols = cols,anno = anno,Margin = c(1000))
dev.off()

if(F){
  ### Parameters
  cols = c("non-depleted"="#4477AA", 
           "Spn1-depleted"="#BB5566", 
           "WT_DMSO"="#DDAA33", 
           "WT_IAA"="#dd8033",
           "antisense_SPN1-AID_DMSO"="#4477AA", 
           "antisense_SPN1-AID_IAA"="#BB5566", 
           "antisense_WT_DMSO"="#DDAA33", 
           "sense_SPN1-AID_DMSO"="#4477AA", 
           "sense_SPN1-AID_IAA"="#BB5566", 
           "sense_WT_DMSO"="#DDAA33",
           "reverse"="gray80","forward"="gray80")
  gene.anno.path="A:/work/yeast_Annotations/ScSp_SGD.RevisedGenesFull.RData"
  genome.anno.path="A:/work/yeast_Annotations/Yeast_ref_Annotations_DJ.gff"
  gene.names <- c("STF1")
  OUTDIR="A:/work/WinstonLab/Natalia_Reim/ThesisFigures_March2019"
  OUTPREFIX="NR_March19"
  
  
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
  plotlist[["WT"]] <- gene.coverage.plot.df(path="A:/work/WinstonLab/Natalia_Reim/October2016/BedGraph/whole_read_replicates merged/WT_DMSO_Mergeq50WR.bedGraph",profile = "WT_DMSO",genes = genes,glist = glist,Margin = 2000)
  plotlist[["SPN1_DMSO"]] <- gene.coverage.plot.df(path="A:/work/WinstonLab/Natalia_Reim/October2016/BedGraph/whole_read_replicates merged/SPN1_DMSO_Mergeq50WR.bedGraph",profile = "SPN1-AID_DMSO",genes = genes,glist = glist,Margin = 2000)
  plotlist[["SPN1_IAA"]] <- gene.coverage.plot.df(path="A:/work/WinstonLab/Natalia_Reim/October2016/BedGraph/whole_read_replicates merged/SPN1_IAA_Mergeq50WR.bedGraph",profile = "SPN1-AID_IAA",genes = genes,glist = glist,Margin = 2000)
  pl <- rbind(plotlist[[1]],plotlist[[2]],plotlist[[3]])
  pl$colguide <- paste0(pl$ori,"_",pl$profile)
  pl$profile <- factor(pl$profile,levels=c("WT_DMSO", "SPN1-AID_DMSO","SPN1-AID_IAA"))
  
  p <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),genes=genes,plot.coverage = F,
                          mygene = "STF1",ori = "",only.gene = F,cols = cols,anno = anno,Margin = c(200,1000))
  
  pdf("A:/work/WinstonLab/Natalia_Reim/ThesisFigures_March2019/STF1-RPP1B_ASplots_v1.pdf",width = 4,height = 4)
  grid.arrange(grobs=p,heights=c(8,2)) #,top= textGrob(paste0(mygene),gp=gpar(fontsize=22,font=3))
  dev.off()
  
  
  ## Antisense expression analysis
  myl <- c("MAF1","CUS2","RRP1","KAR1")
  pdf(paste0(OUTDIR,"/",OUTPREFIX,"_AS_RNASeqCoverageExamples.pdf"),width = 8,height = 6)
  for(mygene in myl){
    cat(mygene,"\n")
    print(gene.coverage.plot(pl=pl,samples = c("WT_DMSO","SPN1-AID_DMSO","SPN1-AID_IAA"),ylim = c(0,20),plot.coverage = T,only.gene=T,stranded = F,mygene = mygene,ori = "antisense",cols = cols,anno = anno,Margin = 200))
    print(gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),ylim = c(0,20),plot.coverage = T,only.gene=T,stranded = F,mygene = mygene,ori = "antisense",cols = cols,anno = anno,Margin = 200))
  }  
  dev.off()
  
  ## sense expression analysis
  myl <- c("HSP82","SER3","UBI4","KAP120","GCV3")
  pdf(paste0(OUTDIR,"/",OUTPREFIX,"_Sense_RNASeqCoverageExamples.pdf"),width = 6,height = 5)
  for(mygene in myl){
    cat(mygene,"\n")
    only.gene=T
    if(mygene=="SER3"){
      only.gene=F
    }
    print(gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),plot.coverage = T,stranded = T,only.gene = only.gene,mygene = mygene,ori = "sense",cols = cols,anno = anno,Margin = 200))
  }  
  dev.off()
  
  ## Intron retention analysis
  myl <- c("RPS27A","RPL43B")
  pdf(paste0(OUTDIR,"/",OUTPREFIX,"_IntronRetention_RNASeqCoverageExamples.pdf"),width = 8,height = 6)
  for(mygene in myl){
    cat(mygene,"\n")
    gene.coverage.plot(pl=pl,samples = c("WT_DMSO","SPN1-AID_DMSO","SPN1-AID_IAA"),plot.coverage = T,only.gene = T,mygene = mygene,ori = "sense",cols = cols,anno = anno,Margin = 200)
  }  
  dev.off()
  
  
  p1 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),plot.coverage = F,only.gene = T,
                           mygene = "KAP120",ori = "sense",cols = cols,anno = anno,Margin = 200)
  p2 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),plot.coverage = F,only.gene = T,
                           mygene = "GCV3",ori = "sense",cols = cols,anno = anno,Margin = 200)
  p3 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),plot.coverage = F,only.gene = T,
                           mygene = "RPS27A",ori = "",cols = cols,anno = anno,Margin = 200)
  p4 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),plot.coverage = F,only.gene = T,
                           mygene = "SER3",ori = "",cols = cols,anno = anno,Margin = 200)
  p5 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),plot.coverage = F,only.gene = T,
                           mygene = "UBI4",ori = "",cols = cols,anno = anno,Margin = 200)
  p6 <- gene.coverage.plot(pl=pl,samples = c("SPN1-AID_DMSO","SPN1-AID_IAA"),plot.coverage = F,only.gene = T,
                           mygene = "HSP82",ori = "",cols = cols,anno = anno,Margin = 200)
  pl.list <- list(p1[[1]],p2[[1]],p3[[1]],p1[[2]],p2[[2]],p3[[2]],
                  p4[[1]],p5[[1]],p6[[1]],p4[[2]],p5[[2]],p6[[2]])
  #save(pl.list,file=paste0(OUTDIR,"/",OUTPREFIX,"_SenseCoveragePLotLists.RData"))
  
  pdf(file=paste0(OUTDIR,"/",OUTPREFIX,"_CombinedCovSensePLots.pdf"),width = 15,height = 10)
  grid.arrange(grobs=pl.list,nrow=4,heights=c(8,1.5,8,1.5)) #,top= textGrob(paste0(mygene),gp=gpar(fontsize=22,font=3))
  dev.off()
  
}

