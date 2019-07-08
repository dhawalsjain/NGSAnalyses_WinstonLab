## l2fc dataframe

if(T){
  rm(list=ls())
  library(psych)
  library(ggplot2)
  library(viridis)
  library(GenomicRanges)
  library(GGally)
  library(ggpubr)
  library(grid)
  library(gridExtra)
  ### add new function with name ggscat to GGally namespace!!
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
   ps=1  
  OUTDIR="A:/work/WinstonLab/Natalia_Reim/ThesisFigures_March2019/"
  
  
  ## RNA expression
  DEExpression.filepath="A:/work/WinstonLab/Natalia_Reim/ThesisFigures_March2019/NR_March19_SenseDE.Data.RData"
  load(DEExpression.filepath)
  res$is_down <- ifelse(res$padj<0.1 & !is.na(res$padj) & res$log2FoldChange<0,1,0)
  expn <- res[,c("tracking_id","log2FoldChange","SPN1-AID_DMSO-1","SPN1-AID_DMSO-2",'width',"is_down")]
  expn[,3:4] <- apply(expn[,3:4],2,function(x) x*1e6/sum(x))
  expn[,3:4] <-  round(expn[,3:4]*1000/expn$width,3)
  expn$basal_rna_exp <- apply(expn[,3:4],1,function(x) sqrt(prod(x))) 
  expn <- expn[,c("tracking_id",'log2FoldChange',"basal_rna_exp","is_down","width")]
  names(expn)[2] <- "expn_l2fc"
  names(expn)[5] <- "gene_width"
  expn$basal_rna_exp <- log2(expn$basal_rna_exp+ps)
  rm(res)
  
  load("A:/work/WinstonLab/Natalia_Reim/DataFrames/Log2FC values for the ChIPFactors_DESeq2.RData")
  load("A:/work/WinstonLab/Natalia_Reim/DataFrames/relOccupancyChIPSeq.RData")
  m.spike$tracking_id <- rownames(m.spike)
  m$tracking_id <- rownames(m)
  fdr.spike$tracking_id <- rownames(fdr.spike)
  fdr$tracking_id <- rownames(fdr)
  ## occupancies
  ##  Rpb1 and modification not spikein normalized
  ##  H3 mods are H3 normalized + spikein norm
  ##  Other modifications are spikein normalized
  o1 <- m.rel.occup
  o2 <- m.occup[,c("tracking_id", "Spn1.Depl_RNAPII_ChIP", "Spn1_RNAPII_ChIP", 
                        "Spn1.Depl_RNAPII (round1)_ChIP", "Spn1_RNAPII (round1)_ChIP", 
                        "Spn1.Depl_S2P_ChIP", "Spn1_S2P_ChIP", "Spn1.Depl_S5P_ChIP", 
                        "Spn1_S5P_ChIP")]
  o3 <- m.spike.occup[,c("tracking_id", "Spn1.Depl_H3_ChIP", "Spn1_H3_ChIP", "Spn1.Depl_H3 (round1)_ChIP", 
                                    "Spn1_H3 (round1)_ChIP", "Spn1.Depl_HA.Set2_ChIP", "Spn1_HA.Set2_ChIP", 
                                    "Spn1.Depl_Spn1_ChIP", "Spn1_Spn1_ChIP", "Spn1.Depl_HA.Spt6_ChIP", 
                                    "Spn1_HA.Spt6_ChIP")]
  o1 <- merge(o1,o2,by="tracking_id")
  o1 <- merge(o1,o3,by="tracking_id")
  occup <- o1
  rm(o2,o3,o1)
  occup[,2:25] <- log2(occup[,2:25]+ps)
  
  ## Log2FC ChIP
  ##  Rpb1 and modification not spikein normalized
  ##  H3 mods are H3 normalized + spikein norm
  ##  Other modifications not all spikein normalized
  q0 <- m.spike[,c("tracking_id","H3", "H3 (round1)", "H3K36me2", "H3K36me3", "H3K4me3")]
  q0$H3K36me2 <- q0$H3K36me2-q0$H3
  q0$H3K4me3 <- q0$H3K4me3-q0$H3
  q0$H3K36me3 <- q0$H3K36me3-q0$`H3 (round1)`
  q1 <- m.spike[,c("tracking_id","Set2", "Spn1", "Spt6")]
  q2 <- m[,c("tracking_id","RNAPII", "RNAPII (round1)", "S2P", "S5P")]
  q0 <- merge(q0,q1,by="tracking_id")
  q0 <- merge(q0,q2,by="tracking_id")
  l2fc <- q0
  rm(q0,q1,q2)
  
  ## Shifts, relative
  ##  Rpb1 and modification not spikein normalized
  ##  H3 mods are H3 normalized + spikein norm
  ##  Other modifications not all spikein normalized
  sh <- read.delim("A:/work/WinstonLab/Natalia_Reim/Paper_Thesis_Figs/ShiftRatio_Pvalues_l2fc_1000.txt",header=T)
  sh <- sh[,c("tracking_id","l2fc_H3","l2fc_H3..round1.")]
  names(sh) <- c("tracking_id","shift_H3","shift_H3_r1")
  shifts <- read.delim("A:/work/WinstonLab/Natalia_Reim/Paper_Thesis_Figs/ShiftRatio_PvaluesL2fc_Rel_1000.txt",header=T)
  shifts <- shifts[,c("tracking_id", "l2fc_H3K36me2", "l2fc_H3K36me3","l2fc_H3K4me3","l2fc_S2P","l2fc_S5P")]
  names(shifts) <- gsub("l2fc","shift",names(shifts))
  shifts <- merge(shifts,sh,by="tracking_id")
  rm(sh)
  
  rm(fdr,fdr.spike,m,m.occup,m.spike,m.spike.occup,m.rel.occup,DEExpression.filepath)
  
  
  ##adjust the names
  h <-   names(l2fc)
  h <- c("tracking_id", "l2fc_H3 (R2)", "l2fc_H3 (R1)", "l2fc_H3K36me2 (R2)", "l2fc_H3K36me3 (R1)", 
         "l2fc_H3K4me3 (R2)", "l2fc_Set2 (R1)", "l2fc_Spn1 (R2)", "l2fc_Spt6 (R1)", "l2fc_RNAPII (R2)",
         "l2fc_RNAPII (R1)", "l2fc_Ser2-P (R2)", "l2fc_Ser5-P (R2)")
  names(l2fc) <- h
  
  i <- names(occup)
  i <- c("tracking_id", "occup_Depl.H3K4me3 (R2)", "occup_H3K4me3 (R2)", 
         "occup_Depl.H3K36me2 (R2)", "occup_H3K36me2 (R2)", 
         "occup_Depl.H3K36me3 (R1)", "occup_H3K36me3 (R1)", "occup_Depl.RNAPII (R2)", "occup_RNAPII (R2)", 
         "occup_Depl.RNAPII (R1)", "occup_RNAPII (R1)", 
         "occup_Depl.Ser2-P (R2)", "occup_Ser2-P (R2)", 
         "occup_Depl.Ser5-P (R2)", "occup_Ser5-P (R2)", "occup_Depl.H3 (R2)", "occup_H3 (R2)", 
         "occup_Depl.H3 (R1)", "occup_H3 (R1)", 
         "occup_Depl.Set2 (R1)", "occup_Set2 (R1)", "occup_Depl.Spn1 (R2)", 
         "occup_Spn1 (R2)", "occup_Depl.Spt6 (R1)", "occup_Spt6 (R1)")
  names(occup) <- i 
  
  j <- names(shifts)
  j <- c("tracking_id", "shift_H3K36me2 (R2)", "shift_H3K36me3 (R1)", "shift_H3K4me3 (R2)",
         "shift_Ser2-P (R2)", "shift_Ser5-P (R2)",
         "shift_H3 (R2)", "shift_H3 (R1)")
  names(shifts) <- j
  rm(h,i,j)
  
  pl <- merge(expn,l2fc,by="tracking_id")
  pl <- merge(pl,occup,by="tracking_id")
  pl <- merge(pl,shifts,by="tracking_id",all.x=T)
  rm(expn,l2fc,occup,shifts)
  
  
  
  #### For introns and exons, use the log2FC values from following
  load("A:/work/WinstonLab/Natalia_Reim/DataFrames/Log2FC values for the ChIPFactors_DESeq2_IntronsExons.RData")
  ## Log2FC IntronsExons
  ##  Rpb1 and modification not spikein normalized
  ##  H3 mods are H3 normalized + spikein norm
  ##  Other modifications not all spikein normalized
  m.spike$tracking_id <- rownames(m.spike)
  q0 <- m.spike[,c("tracking_id","H3", "H3 (round1)", "H3K36me2", "H3K36me3", "H3K4me3")]
  q0$H3K36me2 <- q0$H3K36me2-q0$H3
  q0$H3K4me3 <- q0$H3K4me3-q0$H3
  q0$H3K36me3 <- q0$H3K36me3-q0$`H3 (round1)`
  q1 <- m.spike[,c("tracking_id","Set2", "Spn1", "Spt6")]
  m$tracking_id <- rownames(m)
  q2 <- m[,c("tracking_id","RNAPII", "RNAPII (round1)", "S2P", "S5P")]
  
  q0 <- merge(q0,q1,by="tracking_id")
  q0 <- merge(q0,q2,by="tracking_id")
  l2fci <- q0
  rm(q0,q1,q2,m.rel,m,m.spike,m.spike.occup,fdr, fdr.spike,m.occup)
  h <-   names(l2fci)
  h <- c("tracking_id", "l2fc_H3 (R2)", "l2fc_H3 (R1)", "l2fc_H3K36me2 (R2)", "l2fc_H3K36me3 (R1)", 
         "l2fc_H3K4me3 (R2)", "l2fc_Set2 (R1)", "l2fc_Spn1 (R2)", "l2fc_Spt6 (R1)", "l2fc_RNAPII (R2)",
         "l2fc_RNAPII (R1)", "l2fc_Ser2-P (R2)", "l2fc_Ser5-P (R2)")
  names(l2fci) <- h
  genes <- as.data.frame(genes)
  genes$tracking_id <- paste0(genes$tracking_id,"_",genes$feature,genes$number)
  gl <- merge(l2fci,genes,by="tracking_id",all.x=T)
  rm(h,l2fci,genes)
}
#save(pl,gl,file="A:/work/WinstonLab/Natalia_Reim/DataFrames/NataliaChIPOccupShiftsExpression.RData")

gg <- theme_minimal() +theme(axis.text.x = element_text(size=14,colour = "black"),
            axis.text.y = element_text(size=14,colour = "black"),
            axis.title =  element_text(size=16,color="black"),
            axis.line = element_line(color="gray",size = 0.75),
            axis.ticks = element_line(color="gray",size = 1),
            legend.position = "none",
            plot.title = element_text(colour ="black",size = 18,hjust = 0.5))

### writedown the reports
if(F){
  
}


#########################
### Thesis scatter plots
if(F){
  pl <- pl[pl$gene_width>500,]
  ## 1. Scatterplots, basal occupancy vs basal expression
  if(F){
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
    
    bs <- "basal_rna_exp"
    qs <- c("occup_H3 (R1)", "occup_H3 (R2)", "occup_H3K4me3 (R2)", 
            "occup_H3K36me2 (R2)", "occup_H3K36me3 (R1)", 
            "occup_RNAPII (R1)","occup_RNAPII (R2)", 
            "occup_Ser2-P (R2)", "occup_Ser5-P (R2)",
            "occup_Set2 (R1)", "occup_Spt6 (R1)","occup_Spn1 (R2)")
    pl.list <- list()
    pl.listd <- list()
    for(i in 1:length(qs)){
      xlab="basal expression (log2 RPKM)"
      ylab=paste0(gsub("occup_", "",qs[i]),", log2 occupancy")
      x=as.character(bs)
      y=as.character(qs[i])
      main=gsub("occup_","",as.character(qs[i]))
      cat(main,"\n")
      pl.list[[paste0(main)]] <- print(myplotfun(l2fc = pl,xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
      pl.listd[[paste0(main)]] <- print(myplotfun(l2fc = pl[pl$is_down>0,],xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
      rm(xlab,ylab,main,y,x)
    }
    pdf(file = paste0(OUTDIR,"BasalExpnVsBasalOccupancy.pdf"),width = 16,height = 12)
    grid.arrange(grobs=pl.list,ncol=4,top= textGrob(paste0("basal expression v/s basal occupancy"),gp=gpar(fontsize=20,font=2)))
    grid.arrange(grobs=pl.listd,ncol=4,top= textGrob(paste0("basal expression v/s basal occupancy (downreg.)"),gp=gpar(fontsize=20,font=2)))
    dev.off()  
  }
  ## 2. Scatterplots, basal occupancy vs basal expression
  if(F){
    myplotfun <- function(l2fc,xlab,ylab,x,y,main){
      q <-ggplot(l2fc, aes(x=l2fc[,x], y=l2fc[,y])) +
        xlab(paste0(xlab))+ylab(paste0(ylab))+ggtitle(paste(main))+
        stat_bin_hex(geom="point", aes(color=log2(..count..)), binwidth=c(.05,.05), size=2, shape=16, stroke=0)+
        scale_color_viridis(option="inferno",direction = 1)+
        geom_hline(yintercept = 0,col="black")+geom_vline(xintercept = 0,col="black")+
        geom_smooth(method='lm',formula=y~x,col="red")+
        stat_cor(method = "pearson", 
                 size=5,col="red",label.sep = "\n",
                 label.x.npc = "right", label.y.npc = "bottom")+gg+
        scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))
      return(q)
      #sort(unique(l2fc[,x]),decreasing = T)[2]
    }
    
    bs <- "expn_l2fc"
    qs <- c("occup_H3 (R1)", "occup_H3 (R2)", "occup_H3K4me3 (R2)", 
            "occup_H3K36me2 (R2)", "occup_H3K36me3 (R1)", 
            "occup_RNAPII (R1)","occup_RNAPII (R2)", 
            "occup_Ser2-P (R2)", "occup_Ser5-P (R2)",
            "occup_Set2 (R1)", "occup_Spt6 (R1)","occup_Spn1 (R2)")
    pl.list <- list()
    pl.listd <- list()
    for(i in 1:length(qs)){
      xlab="log2FC, expression"
      ylab=paste0(gsub("occup_", "",qs[i]),", log2 occupancy")
      x=bs
      y=qs[i]
      main=gsub("occup_","",qs[i])
      cat(main,"\n")
      pl.list[[paste0(main)]] <- print(myplotfun(pl,xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
      pl.listd[[paste0(main)]] <- print(myplotfun(pl[pl$is_down>0,],xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
      rm(xlab,ylab,main,x,y)
    }
    pl.list[[1]]
    pdf(file = paste0(OUTDIR,"log2FC ExpressionVsBasalOccupancy.pdf"),width = 16,height = 12)
    grid.arrange(grobs=pl.list,ncol=4,top= textGrob(paste0("log2 fold change, expression v/s basal occupancy"),gp=gpar(fontsize=20,font=2)))
    grid.arrange(grobs=pl.listd,ncol=4,top= textGrob(paste0("log2 fold change, expression v/s basal occupancy (downreg.)"),gp=gpar(fontsize=20,font=2)))
    dev.off()  
  }
  ## 3. Scatterplots, basal occupancy vs basal expression
  if(F){
    myplotfun <- function(l2fc,xlab,ylab,x,y,main){
      q <-ggplot(l2fc, aes(x=l2fc[,x], y=l2fc[,y])) +
        xlab(paste0(xlab))+ylab(paste0(ylab))+ggtitle(paste(main))+
        stat_bin_hex(geom="point", aes(color=log2(..count..)), binwidth=c(.05,.05), size=2, shape=16, stroke=0)+
        scale_color_viridis(option="inferno",direction = 1)+
        geom_hline(yintercept = 0,col="black")+geom_vline(xintercept = 0,col="black")+
        geom_smooth(method='lm',formula=y~x,col="red")+
        stat_cor(method = "pearson", size=5,col="red",label.sep = "\n",
                 label.x.npc = "right", label.y.npc = "bottom")+gg+
        scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))
      return(q)
      #sort(unique(l2fc[,x]),decreasing = T)[2]
    }
    
    bs <- "basal_rna_exp"
    qs <- names(pl)[grep("l2fc",names(pl))]
    qs <- qs[2:13]
    
    pl.list <- list()
    pl.listd <- list()
    for(i in 1:length(qs)){
      xlab="basal expression (log2 RPKM)"
      ylab=paste0(gsub("l2fc_", "",qs[i]),", log2 FC")
      x=bs
      y=qs[i]
      main=gsub("l2fc_","",qs[i])
      cat(main,"\n")
      pl.list[[paste0(main)]] <- print(myplotfun(pl,xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
      pl.listd[[paste0(main)]] <- print(myplotfun(pl[pl$is_down>0,],xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
      rm(xlab,ylab,main,x,y)
    }
    pl.list[[1]]
    pdf(file = paste0(OUTDIR,"BasalExpression vs log2FCOccupancy.pdf"),width = 16,height = 12)
    grid.arrange(grobs=pl.list,ncol=4,top= textGrob(paste0("basal expression v/s log2FC in occupancy"),gp=gpar(fontsize=20,font=2)))
    grid.arrange(grobs=pl.listd,ncol=4,top= textGrob(paste0("basal expression v/s log2FC in occupancy (downreg.)"),gp=gpar(fontsize=20,font=2)))
    dev.off()  
  }
  ## 4. Scatterplots, basal occupancy vs basal expression
  if(F){
    myplotfun <- function(l2fc,xlab,ylab,x,y,main){
      q <-ggplot(l2fc, aes(x=l2fc[,x], y=l2fc[,y])) +
        xlab(paste0(xlab))+ylab(paste0(ylab))+ggtitle(paste(main))+
        stat_bin_hex(geom="point", aes(color=log2(..count..)), binwidth=c(.05,.05), size=2, shape=16, stroke=0)+
        scale_color_viridis(option="inferno",direction = 1)+
        geom_hline(yintercept = 0,col="black")+geom_vline(xintercept = 0,col="black")+
        geom_smooth(method='lm',formula=y~x,col="red")+
        stat_cor(method = "pearson", size=5,col="red",label.sep = "\n",
                 label.x.npc = "right", label.y.npc = "bottom")+gg+
        scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))
      return(q)
    }
    
    bs <- "basal_rna_exp"
    qs <- names(pl)[grep("shift",names(pl))]
    
    pl.list <- list()
    pl.listd <- list()
    for(i in 1:length(qs)){
      xlab="basal expression (log2 RPKM)"
      ylab=paste0(gsub("shift_", "",qs[i]),", log2 FC")
      x=bs
      y=qs[i]
      main=gsub("shift_","",qs[i])
      cat(main,"\n")
      pl.list[[paste0(main)]] <- print(myplotfun(pl,xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
      pl.listd[[paste0(main)]] <- print(myplotfun(pl[pl$is_down>0,],xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
      rm(xlab,ylab,main,x,y)
    }
    pl.list[[1]]
    pdf(file = paste0(OUTDIR,"BasalExpression vs Shifts.pdf"),width = 16,height = 8)
    grid.arrange(grobs=pl.list,ncol=4,top= textGrob(paste0("basal expression v/s log2FC in shifts"),gp=gpar(fontsize=20,font=2)))
    grid.arrange(grobs=pl.listd,ncol=4,top= textGrob(paste0("basal expression v/s log2FC in shifts (downreg.)"),gp=gpar(fontsize=20,font=2)))
    dev.off()  
  }
  ## 5. Scatterplots, basal occupancy vs basal expression
  if(F){
    myplotfun <- function(l2fc,xlab,ylab,x,y,main){
      q <-ggplot(l2fc, aes(x=l2fc[,x], y=l2fc[,y])) +
        xlab(paste0(xlab))+ylab(paste0(ylab))+ggtitle(paste(main))+
        stat_bin_hex(geom="point", aes(color=log2(..count..)), binwidth=c(.05,.05), size=2, shape=16, stroke=0)+
        scale_color_viridis(option="inferno",direction = 1)+
        geom_hline(yintercept = 0,col="black")+geom_vline(xintercept = 0,col="black")+
        geom_smooth(method='lm',formula=y~x,col="red")+
        stat_cor(method = "pearson", size=5,col="red",label.sep = "\n",
                 label.x.npc = "right", label.y.npc = "bottom")+gg+
        scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))
      return(q)
      #sort(unique(l2fc[,x]),decreasing = T)[2]
    }
    
    bs <- "expn_l2fc"
    qs <- names(pl)[grep("l2fc",names(pl))]
    qs <- qs[2:13]
    
    pl.list <- list()
    pl.listd <- list()
    for(i in 1:length(qs)){
      xlab="log2FC, expression"
      ylab=paste0(gsub("l2fc_", "",qs[i]),", log2 FC")
      x=bs
      y=qs[i]
      main=gsub("l2fc_","",qs[i])
      cat(main,"\n")
      pl.list[[paste0(main)]] <- print(myplotfun(pl,xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
      pl.listd[[paste0(main)]] <- print(myplotfun(pl[pl$is_down>0,],xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
    }
    pl.list[[1]]
    pdf(file = paste0(OUTDIR,"ExpressionL2FC vs OccupancyL2FC.pdf"),width = 16,height = 12)
    grid.arrange(grobs=pl.list,ncol=4,top= textGrob(paste0("log2FC expression change v/s log2FC occupancy change"),gp=gpar(fontsize=20,font=2)))
    grid.arrange(grobs=pl.listd,ncol=4,top= textGrob(paste0("log2FC expression change v/s log2FC occupancy change (downreg.)"),gp=gpar(fontsize=20,font=2)))
    dev.off()  
  }
  ## 6. Scatterplots, basal occupancy vs basal expression
  if(F){
    myplotfun <- function(l2fc,xlab,ylab,x,y,main){
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
      #sort(unique(l2fc[,x]),decreasing = T)[2]
    }
    
    bs <- "expn_l2fc"
    qs <- names(pl)[grep("shift",names(pl))]
    
    pl.list <- list()
    pl.listd <- list()
    for(i in 1:length(qs)){
      xlab="log2FC, expression"
      ylab=paste0(gsub("shift_", "",qs[i]),", log2 FC")
      x=bs
      y=qs[i]
      main=gsub("shift_","",qs[i])
      cat(main,"\n")
      pl.list[[paste0(main)]] <- print(myplotfun(pl,xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
      pl.listd[[paste0(main)]] <- print(myplotfun(pl[pl$is_down>0,],xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
      #rm(xlab,ylab,main,x,y)
    }
    pl.list[[1]]
    pdf(file = paste0(OUTDIR,"log2FCExpression vs Shifts.pdf"),width = 16,height = 8)
    grid.arrange(grobs=pl.list,ncol=4,top= textGrob(paste0("log2FC expression v/s log2FC in shifts"),gp=gpar(fontsize=20,font=2)))
    grid.arrange(grobs=pl.listd,ncol=4,top= textGrob(paste0("log2FC expression v/s log2FC in shifts (downreg.)"),gp=gpar(fontsize=20,font=2)))
    dev.off()  
  }
  ## 7. Pairwise comparisons
  if(T){
    qs <- names(pl)[grep("occup_",names(pl))]
    qs <- qs[-grep("Depl",qs)]
    
    m <- pl[,qs]
    names(m) <- gsub("occup_","",names(m))
    q <- ggscat(m, columns = 1:ncol(m),color = NULL,size=10)+
      stat_bin_hex(geom="point", aes(color=log2(..count..)), binwidth=c(.05,.05), size=2, shape=16, stroke=0)+
      scale_color_viridis(option="inferno",direction = 1)+
      theme(axis.text.x = element_text(size=14,colour = "black"),
            axis.text.y = element_text(size=14,colour = "black"),
            axis.title =  element_text(size=16,color="black"),
            strip.text = element_text(color="black",size=16))+
      theme(legend.position = "none")+
      theme(panel.background = element_rect(fill = "white", colour = "black"),
            panel.grid = element_blank())+
      geom_smooth(method='lm',formula=y~x,col="red")+
      xlab("factor occupancy,log2")+ylab("factor occupancy, log2")+
      ggtitle("Pairwise basal occupancy comparison")+
      theme(plot.title = element_text(colour ="black",size = 18,hjust = 0.5))
    png(file = paste0(OUTDIR,"Occupancy log2 vs Occupancy log2.png"),width = 2000,height = 2000)
    q
    dev.off()
    
    m <- pl[pl$is_down>0,qs]
    names(m) <- gsub("occup_","",names(m))
    q <- ggscat(m, columns = 1:ncol(m),color = NULL,size=10)+
      stat_bin_hex(geom="point", aes(color=log2(..count..)), binwidth=c(.05,.05), size=2, shape=16, stroke=0)+
      scale_color_viridis(option="inferno",direction = 1)+
      theme(axis.text.x = element_text(size=14,colour = "black"),
            axis.text.y = element_text(size=14,colour = "black"),
            axis.title =  element_text(size=16,color="black"),
            strip.text = element_text(color="black",size=16))+
      theme(legend.position = "none")+
      theme(panel.background = element_rect(fill = "white", colour = "black"),
            panel.grid = element_blank())+
      geom_smooth(method='lm',formula=y~x,col="red")+
      xlab("factor occupancy,log2")+ylab("factor occupancy, log2")+
      ggtitle("Pairwise basal occupancy comparison")+
      theme(plot.title = element_text(colour ="black",size = 18,hjust = 0.5))
    png(file = paste0(OUTDIR,"Occupancy log2 vs Occupancy log2_DownReg.png"),width = 2000,height = 2000)
    q
    dev.off()
    
  }
  ## 8. Basal occupancy vs Occupancy change
  if(T){
    myplotfun <- function(l2fc,xlab,ylab,x,y,main){
      q <-ggplot(l2fc, aes(x=l2fc[,x], y=l2fc[,y])) +
        xlab(paste0(xlab))+ylab(paste0(ylab))+ggtitle(paste(main))+
        stat_bin_hex(geom="point", aes(color=log2(..count..)), binwidth=c(.05,.05), size=2, shape=16, stroke=0)+
        scale_color_viridis(option="inferno",direction = 1)+
        geom_hline(yintercept = 0,col="black")+geom_vline(xintercept = 0,col="black")+
        geom_smooth(method='lm',formula=y~x,col="red")+
        stat_cor(method = "pearson",size=5,col="red",label.sep = "\n",
                 label.x.npc = "left", label.y.npc = "top")+gg+
        scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))
      return(q)
      #sort(unique(l2fc[,x]),decreasing = T)[2]
    }
    
    qs <- names(pl)[grep("occup_",names(pl))]
    qs <- qs[-grep("Depl",qs)]
    bs <- gsub("occup","l2fc",qs)
    
    
    pl.list <- list()
    pl.listd <- list()
    for(i in 1:length(qs)){
      for(j in 1:length(bs)){
        xlab=paste0(gsub("occup_","",qs[i])," occupancy, log2")
        ylab=paste0("log2FC ",gsub("l2fc_","",bs[j]))
        x=qs[i]
        y=bs[j]
        main=""
        cat(qs[i]," ",bs[j],"\n")
        pl.list[[paste0(qs[i]," ",bs[j])]] <- print(myplotfun(pl,xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
        pl.listd[[paste0(qs[i]," ",bs[j])]] <- print(myplotfun(pl[pl$is_down>0,],xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
        rm(xlab,ylab,main,x,y)
      }
    }
    
    
    pdf(file = paste0(OUTDIR,"BasalOccupancy vs OccupancyChange.pdf"),width = 16,height = 12)
    for(i in seq(1,144,12)){
      h <- gsub(" l2fc_.*$","",sub("occup_","",names(pl.list)[i]))
      grid.arrange(grobs=pl.list[i:(i+11)],ncol=4,
                   top= textGrob(paste0(h)," occupancy v/s occupancy change"),gp=gpar(fontsize=20,font=2))
    }
    dev.off()  
    
    pdf(file = paste0(OUTDIR,"BasalOccupancy vs OccupancyChange downreg.pdf"),width = 16,height = 12)
    for(i in seq(1,144,12)){
      h <- gsub(" l2fc_.*$","",sub("occup_","",names(pl.listd)[i]))
      grid.arrange(grobs=pl.listd[i:(i+11)],ncol=4,
                   top= textGrob(paste0(h)," occupancy v/s occupancy change (downreg.)"),gp=gpar(fontsize=20,font=2))
    }
    dev.off()  
    
    
  }
  ## 9. Basal occupancy vs Occupancy change
  if(T){
    myplotfun <- function(l2fc,xlab,ylab,x,y,main){
      q <-ggplot(l2fc, aes(x=l2fc[,x], y=l2fc[,y])) +
        xlab(paste0(xlab))+ylab(paste0(ylab))+ggtitle(paste(main))+
        stat_bin_hex(geom="point", aes(color=log2(..count..)), binwidth=c(.05,.05), size=2, shape=16, stroke=0)+
        scale_color_viridis(option="inferno",direction = 1)+
        geom_hline(yintercept = 0,col="black")+geom_vline(xintercept = 0,col="black")+
        geom_smooth(method='lm',formula=y~x,col="red")+
        stat_cor(method = "pearson",size=5,col="red",label.sep = "\n",
                 label.x.npc = "left", label.y.npc = "top")+gg+
        scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))
      return(q)
      #sort(unique(l2fc[,x]),decreasing = T)[2]
    }
    
    bs <- names(pl)[grep("shift_",names(pl))]
    qs <- names(pl)[grep("occup_",names(pl))]
    qs <- qs[-grep("Depl",qs)]
    pl.list <- list()
    pl.listd <- list()
    
    for(j in 1:length(bs)){
      for(i in 1:length(qs)){
        xlab=paste0(gsub("occup_","",qs[i]), " occupancy, log2")
        ylab=paste0(gsub("shift_","",bs[j]), " shift, log2")
        y=bs[j]
        x=qs[i]
        main=""
        cat(x,"\n")
        pl.list[[paste0(main,"_",i,"_",j)]] <- print(myplotfun(pl,xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
        pl.listd[[paste0(main,"_",i,"_",j)]] <- print(myplotfun(pl[pl$is_down>0,],xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
        rm(xlab,ylab,main,x,y)
      }
      rm(i,j)
    }
    
    pdf(file = paste0(OUTDIR,"BasalOccupancy vs signal shift.pdf"),width = 16,height = 12)
    for(i in seq(1,73,12)){
      j <- i+11
      grid.arrange(grobs=pl.list[i:j],ncol=4,top= textGrob(paste0("basal occupancy v/s shift"),gp=gpar(fontsize=20,font=2)))
    }
    grid.newpage()
    for(i in seq(1,73,12)){
      j <- i+11
      grid.arrange(grobs=pl.listd[i:j],ncol=4,top= textGrob(paste0("basal occupancy v/s shift (downreg.)"),gp=gpar(fontsize=20,font=2)))
    }
    dev.off()  
    
    
    
  }
  ## 10. occup fc vs occup fc
  if(T){
    qs <- names(pl)[grep("l2fc_",names(pl))]
    
    m <- pl[,qs]
    names(m) <- gsub("l2fc_","",names(m))
    q <- ggscat(m, columns = 1:ncol(m),color = NULL,size=10)+
      stat_bin_hex(geom="point", aes(color=log2(..count..)), binwidth=c(.05,.05), size=2, shape=16, stroke=0)+
      scale_color_viridis(option="inferno",direction = 1)+
      theme(axis.text.x = element_text(size=14,colour = "black"),
            axis.text.y = element_text(size=14,colour = "black"),
            axis.title =  element_text(size=16,color="black"),
            strip.text = element_text(color="black",size=16))+
      theme(legend.position = "none")+
      theme(panel.background = element_rect(fill = "white", colour = "black"),
            panel.grid = element_blank())+
      geom_smooth(method='lm',formula=y~x,col="red")+
      xlab("Log2FC, factor occupancy")+ylab("Log2FC, factor occupancy")+
      ggtitle("Pairwise Occupancy change comparison")+
      theme(plot.title = element_text(colour ="black",size = 18,hjust = 0.5))
    png(file = paste0(OUTDIR,"Occupancy log2fc vs Occupancy log2fc.png"),width = 2000,height = 2000)
    q
    dev.off()
    
    m <- pl[pl$is_down>0,qs]
    names(m) <- gsub("l2fc_","",names(m))
    q <- ggscat(m, columns = 1:ncol(m),color = NULL,size=10)+
      stat_bin_hex(geom="point", aes(color=log2(..count..)), binwidth=c(.05,.05), size=2, shape=16, stroke=0)+
      scale_color_viridis(option="inferno",direction = 1)+
      theme(axis.text.x = element_text(size=14,colour = "black"),
            axis.text.y = element_text(size=14,colour = "black"),
            axis.title =  element_text(size=16,color="black"),
            strip.text = element_text(color="black",size=16))+
      theme(legend.position = "none")+
      theme(panel.background = element_rect(fill = "white", colour = "black"),
            panel.grid = element_blank())+
      geom_smooth(method='lm',formula=y~x,col="red")+
      xlab("Log2FC, factor occupancy")+ylab("Log2FC, factor occupancy")+
      ggtitle("Pairwise Occupancy change comparison")+
      theme(plot.title = element_text(colour ="black",size = 18,hjust = 0.5))
    png(file = paste0(OUTDIR,"Occupancy log2fc vs Occupancy log2fc_DownReg.png"),width = 2000,height = 2000)
    q
    dev.off()
    
  }
  ## 11. mods vs mods 
  if(T){
    qs <- names(pl)[grep("shift_",names(pl))]
    
    m <- pl[,qs]
    names(m) <- gsub("shift_","",names(m))
    q <- ggscat(m, columns = 1:ncol(m),color = NULL,size=10)+
      stat_bin_hex(geom="point", aes(color=log2(..count..)), binwidth=c(.05,.05), size=2, shape=16, stroke=0)+
      scale_color_viridis(option="inferno",direction = 1)+
      theme(axis.text.x = element_text(size=14,colour = "black"),
            axis.text.y = element_text(size=14,colour = "black"),
            axis.title =  element_text(size=16,color="black"),
            strip.text = element_text(color="black",size=16))+
      theme(legend.position = "none")+
      theme(panel.background = element_rect(fill = "white", colour = "black"),
            panel.grid = element_blank())+
      geom_smooth(method='lm',formula=y~x,col="red")+
      xlab("signal shift, log2FC")+ylab("signal shift, log2FC")+
      ggtitle("Pairwise shift comparison")+
      theme(plot.title = element_text(colour ="black",size = 18,hjust = 0.5))
    png(file = paste0(OUTDIR,"shifts vs shifts.png"),width = 1000,height = 1000)
    q
    dev.off()
    
    m <- pl[pl$is_down>0,qs]
    names(m) <- gsub("shift_","",names(m))
    q <- ggscat(m, columns = 1:ncol(m),color = NULL,size=10)+
      stat_bin_hex(geom="point", aes(color=log2(..count..)), binwidth=c(.05,.05), size=2, shape=16, stroke=0)+
      scale_color_viridis(option="inferno",direction = 1)+
      theme(axis.text.x = element_text(size=14,colour = "black"),
            axis.text.y = element_text(size=14,colour = "black"),
            axis.title =  element_text(size=16,color="black"),
            strip.text = element_text(color="black",size=16))+
      theme(legend.position = "none")+
      theme(panel.background = element_rect(fill = "white", colour = "black"),
            panel.grid = element_blank())+
      geom_smooth(method='lm',formula=y~x,col="red")+
      xlab("signal shift, log2FC")+ylab("signal shift, log2FC")+
      ggtitle("Pairwise shift comparison")+
      theme(plot.title = element_text(colour ="black",size = 18,hjust = 0.5))
    png(file = paste0(OUTDIR,"shifts vs shifts_DownReg.png"),width = 1000,height = 1000)
    q
    dev.off()
    
  }
  ## 12. occupancy change vs shift change
  if(T){
    myplotfun <- function(l2fc,xlab,ylab,x,y,main){
      q <-ggplot(l2fc, aes(x=l2fc[,x], y=l2fc[,y])) +
        xlab(paste0(xlab))+ylab(paste0(ylab))+ggtitle(paste(main))+
        stat_bin_hex(geom="point", aes(color=log2(..count..)), binwidth=c(.05,.05), size=2, shape=16, stroke=0)+
        scale_color_viridis(option="inferno",direction = 1)+
        geom_hline(yintercept = 0,col="black")+geom_vline(xintercept = 0,col="black")+
        geom_smooth(method='lm',formula=y~x,col="red")+
        stat_cor(method = "pearson",size=5,col="red",label.sep = "\n",
                 label.x.npc = "right", label.y.npc = "top")+gg+
        scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))
      return(q)
      #sort(unique(l2fc[,x]),decreasing = T)[2]
    }
    
    bs <- names(pl)[grep("shift_",names(pl))]
    qs <- names(pl)[grep("l2fc_",names(pl))]
    pl.list <- list()
    pl.listd <- list()
    
    for(j in 1:length(bs)){
      for(i in 1:length(qs)){
        xlab=paste0(gsub("l2fc_","",qs[i]), " log2FC")
        ylab=paste0(gsub("shift_","",bs[j]), " shift, log2")
        y=bs[j]
        x=qs[i]
        main=""
        cat(x,"\n")
        pl.list[[paste0(main,"_",i,"_",j)]] <- print(myplotfun(pl,xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
        pl.listd[[paste0(main,"_",i,"_",j)]] <- print(myplotfun(pl[pl$is_down>0,],xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
        rm(xlab,ylab,main,x,y)
      }
      rm(i,j)
    }
    
    pdf(file = paste0(OUTDIR,"l2fc occupancy change vs signal shift.pdf"),width = 16,height = 12)
    for(i in seq(1,73,12)){
      j <- i+11
      grid.arrange(grobs=pl.list[i:j],ncol=4,top= textGrob(paste0("occupancy change v/s shift"),gp=gpar(fontsize=20,font=2)))
    }
    grid.newpage()
    for(i in seq(1,73,12)){
      j <- i+11
      grid.arrange(grobs=pl.listd[i:j],ncol=4,top= textGrob(paste0("occupancy change v/s shift (downreg.)"),gp=gpar(fontsize=20,font=2)))
    }
    dev.off()  
  }
  
}



### Boxplot quantification of factor Shift and Occupancy for intron containing genes
if(F){
  ## get pl dataframe ffrom above. keep gense <500bp
  load(file="A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Gene groups.RData")
  load("A:/work/WinstonLab/Natalia_Reim/ThesisFigures_March2019/NR_March19_10pFDR_downGenes.RData")
  rp <- rp[rp$gene!="RPP1",]
  rp <- rp[-grep("MRP",rp$gene),]
  
  m <- pl[,c(1,grep("l2fc|shift",names(pl)))]
  m$rp <- ifelse(m$tracking_id%in%rp$tracking_id,1,0)
  sum(m$rp>0)
  m$introns <- ifelse(m$tracking_id%in%introns$tracking_id,1,0)
  sum(m$introns>0)
  cf <- merge(m,genes,by="tracking_id",all.x=T)
  write.table(cf,file=paste0(OUTDIR,"ICG-genes with AlteredOccupancy and Shifts.txt"),quote = F,sep = "\t",row.names = F)
  
  #m <- m[,-grep("occup",names(m))]
  a <- m[m$rp==1 & m$introns==1,] ## RP-ICG
  b <- m[m$rp==0 & m$introns==1,] ## nonRP-ICG
  d <- m[m$introns==1,] ## All-ICG
  a <- reshape::melt(a,measure.var=names(a)[2:21])
  a$category <- "RP-ICG"
  b <- reshape::melt(b,measure.var=names(b)[2:21])
  b$category <- "nonRP-ICG"
  d <- reshape::melt(d,measure.var=names(d)[2:21])
  d$category <- "All-ICG"
  z <- rbind(a,b,d)
  rm(a,b,d,m)
  
  z$quant <- NA
  z[grep("shift_",z$variable),]$quant <- "Shifts"
  z[grep("l2fc_",z$variable),]$quant <- "Occupancy"
  z$rp <- z$introns <- NULL

  names(z) <- c("tracking_id","factor","value","gene_class","quant")
  z$factor <- gsub("shift_","",z$factor)
  z$factor <- gsub("l2fc_","",z$factor)
  z$factor <- gsub("S2P","Ser2-P",z$factor)
  z$factor <- gsub("S5P","Ser5-P",z$factor)
  z$factor <- gsub("RNAPII","Rpb1",z$factor)
  table(z$factor)
  
  cols = c("All-ICG"="brown",
           "RP-ICG"="red","nonRP-ICG"="pink",
           "Spn1-depleted" = "#BB5566","non-depleted" = "#4477AA")
  z$gene_class <- factor(z$gene_class,levels=c("All-ICG","RP-ICG", "nonRP-ICG"))
  z <- z[!is.na(z$quant),]
  
  q <- ggplot(z[z$quant=="Occupancy",], aes(x=gene_class,y=value,fill=gene_class))
  q <- q + geom_violin(draw_quantiles = c(0.5))
  q <- q + scale_fill_manual(values = cols)
  q <- q + xlab("") + ylab("occupancy fold change, log2")
  q <- q + facet_wrap(~factor,nrow = 2,scales = "free_x")+theme(strip.text = element_text(size=16,colour = "black"))
  q <- q + geom_hline(yintercept = c(0,seq(min(z[z$quant=="Occupancy",]$value), max(z[z$quant=="Occupancy",]$value),1)),linetype="dashed",col="gray")
  q <- q + theme(panel.background = element_rect(fill = NA),
                 axis.line = element_line(size = 0.7, colour = "gray50"))
  q <- q + theme(axis.text.x = element_blank(),
                 axis.text.y = element_text(colour = 'black',size=15),
                 axis.title = element_text(colour = 'black',size=16))
  q <- q + theme(legend.text =element_text(colour = 'black',size=15),
                 legend.title = element_text(colour = 'black',size=15,face="bold"),
                 legend.position = "bottom",
                 legend.justification = "center")
  q <- q + ggtitle("Occupancy changes over genes") + theme(title = element_text(color="black",size = 18))
  
  z <- z[!is.na(z$value),]
  w <- ggplot(z[z$quant=="Shifts",], aes(x=gene_class,y=value,fill=gene_class))
  w <- w + geom_violin(draw_quantiles = c(0.5))
  w <- w + scale_fill_manual(values = cols)
  w <- w + xlab("") + ylab("shift ratio fold change, log2")
  w <- w + facet_wrap(~factor,nrow = 2,scales = "free_x")+theme(strip.text = element_text(size=16,colour = "black"))
  w <- w + geom_hline(yintercept = c(0,seq(min(z[z$quant=="Shifts",]$value), max(z[z$quant=="Shifts",]$value),1)),linetype="dashed",col="gray")
  w <- w + theme(panel.background = element_rect(fill = NA),
                 axis.line = element_line(size = 0.7, colour = "gray50"))
  w <- w + theme(axis.text.x = element_blank(),
                 axis.text.y = element_text(colour = 'black',size=15),
                 axis.title = element_text(colour = 'black',size=16))
  w <- w + theme(legend.text =element_text(colour = 'black',size=15),
                 legend.title = element_text(colour = 'black',size=15,face="bold"),
                 legend.position = "bottom",
                 legend.justification = "center")
  w <- w + ggtitle("Shift ratio changes over genes") + theme(title = element_text(color="black",size = 18))
  w
  
  pdf(file = paste0(OUTDIR,"OccupancyNShift changes on IGCS.pdf"),height = 6,width = 8)
  q
  w
  dev.off()
  ## 
  ## 
  
}


### l2fc over introns and exons
if(T){
  load(file="A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Gene groups.RData")
  load("A:/work/WinstonLab/Natalia_Reim/ThesisFigures_March2019/NR_March19_10pFDR_downGenes.RData")
  load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenesFull.RData")
  rp <- rp[rp$gene!="RPP1",]
  rp <- rp[-grep("MRP",rp$gene),]
  
  genes$gene <- ifelse(is.na(genes$gene),as.character(genes$tracking_id),as.character(genes$gene))
  genes$rp <- ifelse(genes$tracking_id%in%rp$tracking_id,1,0)
  genes$de.genes <- ifelse(genes$tracking_id%in%de.genes$tracking_id,1,0)
  z <- genes[,c("tracking_id","width","rp","de.genes")]
  names(z) <- c("id","g.width","rp","de.genes")
  gl$id <- as.character(gsub("_\\S*","",gl$tracking_id))
  gl <- merge(gl,z,by="id",all.x=T)
  rm(z,genes,cf,dpn,introns,opn,rp,de.genes,opndpn)
  gl$feature <- paste0(gl$feature,gl$number)
  
  m <- gl
  
  a <- m[m$rp==1,] ## RP-ICG
  b <- m[m$rp==0,] ## nonRP-ICG
  d <- m ## All-ICG
  
  a <- reshape::melt(a,measure.var=names(a)[3:14])
  a$category <- "RP-ICG"
  b <- reshape::melt(b,measure.var=names(b)[3:14])
  b$category <- "nonRP-ICG"
  d <- reshape::melt(d,measure.var=names(d)[3:14])
  d$category <- "All-ICG"
  
  z <- rbind(a,b,d)
  z <- z[,-c(3:8,10:13)]
    
    
  names(z) <- c("tracking_id","id","feature","factor","value","gene_class")
  z$factor <- gsub("l2fc_","",z$factor)
  z$factor <- gsub("S2P","Ser2-P",z$factor)
  z$factor <- gsub("S5P","Ser5-P",z$factor)
  z$factor <- gsub("RNAPII","Rpb1",z$factor)
  table(z$factor)
  cols = c("All-ICG"="brown",
           "RP-ICG"="red","nonRP-ICG"="pink",
           "Spn1-depleted" = "#BB5566","non-depleted" = "#4477AA")
  
  z$gene_class <- factor(z$gene_class,levels=c("All-ICG","RP-ICG", "nonRP-ICG"))
  z$value <- ifelse(is.na(z$value),0,z$value)
  
  w1 <- ggplot(z[z$feature=="intron1",], aes(x=gene_class,y=value,fill=gene_class))+
        geom_violin(draw_quantiles = c(0.5))+
        scale_fill_manual(values = cols)+
        xlab("") + ylab("occupancy fold change, log2")+
        facet_wrap(~factor,nrow = 2,scales = "free_x")+
        geom_hline(yintercept = c(0,seq(min(z[z$feature=="intron1",]$value), max(z[z$feature=="intron1",]$value),1)),linetype="dashed",col="gray")+
        theme(strip.text = element_text(size=16,colour = "black"),
              panel.background = element_rect(fill = NA),
              axis.line = element_line(size = 0.7, colour = "gray50"),
              axis.text.x = element_blank(),
              axis.text.y = element_text(colour = 'black',size=15),
              axis.title = element_text(colour = 'black',size=16),
              legend.text =element_text(colour = 'black',size=15),
              legend.title = element_text(colour = 'black',size=15,face="bold"),
              legend.position = "bottom",
              legend.justification = "center",
              title = element_text(color="black",size = 18))+
         ggtitle("Occupancy change over intron1")
  w2 <- ggplot(z[z$feature=="intron2",], aes(x=gene_class,y=value,fill=gene_class))+
    geom_violin(draw_quantiles = c(0.5))+
    scale_fill_manual(values = cols)+
    xlab("") + ylab("occupancy fold change, log2")+
    facet_wrap(~factor,nrow = 2,scales = "free_x")+
    geom_hline(yintercept = c(0,seq(min(z[z$feature=="intron2",]$value), max(z[z$feature=="intron2",]$value),1)),linetype="dashed",col="gray")+
    theme(strip.text = element_text(size=16,colour = "black"),
          panel.background = element_rect(fill = NA),
          axis.line = element_line(size = 0.7, colour = "gray50"),
          axis.text.x = element_blank(),
          axis.text.y = element_text(colour = 'black',size=15),
          axis.title = element_text(colour = 'black',size=16),
          legend.text =element_text(colour = 'black',size=15),
          legend.title = element_text(colour = 'black',size=15,face="bold"),
          legend.position = "bottom",
          legend.justification = "center",
          title = element_text(color="black",size = 18))+
    ggtitle("Occupancy change over intron2")
  
  w3 <- ggplot(z[z$feature=="exon1",], aes(x=gene_class,y=value,fill=gene_class))+
           geom_violin(draw_quantiles = c(0.5))+
           scale_fill_manual(values = cols)+
           xlab("") + ylab("occupancy fold change, log2")+
           facet_wrap(~factor,nrow = 2,scales = "free_x")+
           geom_hline(yintercept = c(0,seq(min(z[z$feature=="exon1",]$value), max(z[z$feature=="exon1",]$value),1)),linetype="dashed",col="gray")+
           theme(strip.text = element_text(size=16,colour = "black"),
                 panel.background = element_rect(fill = NA),
                 axis.line = element_line(size = 0.7, colour = "gray50"),
                 axis.text.x = element_blank(),
                 axis.text.y = element_text(colour = 'black',size=15),
                 axis.title = element_text(colour = 'black',size=16),
                 legend.text =element_text(colour = 'black',size=15),
                 legend.title = element_text(colour = 'black',size=15,face="bold"),
                 legend.position = "bottom",
                 legend.justification = "center",
                 title = element_text(color="black",size = 18))+
         ggtitle("Occupancy change over exon1")
         
   w4 <-ggplot(z[z$feature=="exon2",], aes(x=gene_class,y=value,fill=gene_class))+
           geom_violin(draw_quantiles = c(0.5))+
           scale_fill_manual(values = cols)+
           xlab("") + ylab("occupancy fold change, log2")+
           facet_wrap(~factor,nrow = 2,scales = "free_x")+
           geom_hline(yintercept = c(0,seq(min(z[z$feature=="exon2",]$value), max(z[z$feature=="exon2",]$value),1)),linetype="dashed",col="gray")+
           theme(strip.text = element_text(size=16,colour = "black"),
                 panel.background = element_rect(fill = NA),
                 axis.line = element_line(size = 0.7, colour = "gray50"),
                 axis.text.x = element_blank(),
                 axis.text.y = element_text(colour = 'black',size=15),
                 axis.title = element_text(colour = 'black',size=16),
                 legend.text =element_text(colour = 'black',size=15),
                 legend.title = element_text(colour = 'black',size=15,face="bold"),
                 legend.position = "bottom",
                 legend.justification = "center",
                 title = element_text(color="black",size = 18))+
         ggtitle("Occupancy change over exon2")
  

   pdf(file = paste0(OUTDIR,"PF_ShiftsOccupancy_QuantViolinePlots_I2I2E1E2.pdf"),width = 9,height = 6)
   print(w1)
   print(w2)
   print(w3)
   print(w4)
   dev.off() 
                
}




### Scatter plot for ShiftChange vs Intron length/GeneLength
if(F){
  load("A:/work/WinstonLab/Natalia_Reim/ThesisFigures_March2019/NR_March19_10pFDR_downGenes.RData")
  load("A:/work/yeast_Annotations/ScSp_SGD.RevisedGenesFull.RData")
  load(file="A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Gene groups.RData")
  rp <- rp[rp$gene!="RPP1",]
  rp <- rp[-grep("MRP",rp$gene),]
  
  genes$gene <- ifelse(is.na(genes$gene),as.character(genes$tracking_id),as.character(genes$gene))
  genes$rp <- ifelse(genes$tracking_id%in%rp$tracking_id,1,0)
  genes$de.genes <- ifelse(genes$tracking_id%in%de.genes$tracking_id,1,0)
  z <- genes[,c("tracking_id","width","rp","de.genes")]
  names(z) <- c("id","g.width","rp","de.genes")
  gl$id <- as.character(gsub("_\\S*","",gl$tracking_id))
  gl <- merge(gl,z,by="id",all.x=T)
  rm(z,genes,cf,dpn,introns,opn,de.genes,opndpn)
  
  names(gl) <- gsub("l2fc","l2fcIE",names(gl))
  
  z <- pl[,grep("shift_|l2fc_|tracking_id",names(pl))]
  names(z)[1] <- "id"
  gl <- merge(gl,z,by="id",all.x=T)
  
  gl$g.width <- log2(gl$g.width)
  gl$width <- log2(gl$width+1)
  gl$gene <- ifelse(is.na(gl$gene),as.character(gl$id),as.character(gl$gene))
  pl$gene_width <- log2(pl$gene_width)
  names(gl)[c(18,23)] <- c("intron_width","gene_width")
  pl$is_icg <- ifelse(pl$tracking_id%in%unique(gl$id),1,0)
  gl$feature <- paste0(gl$feature,gl$number)
  
  gl$rp <- ifelse(gl$rp==1, as.character("RP-ICG"),as.character("nonRP-ICG"))
  #All ChIP occupancies over entire ICG
  #Shifts over entire ICG
  #ChIP occupancies over:
  #  First exon
  #  Intron
  #  Second exon
  #Is there a correlation between modification defects and lengths of first exon, intron, and/or second exon?
  
  myplotfun <- function(l2fc,xlab,ylab,main,x,y){
    q <-ggplot(l2fc, aes(x=l2fc[,x], y=l2fc[,y])) +
      xlab(paste0(xlab))+ylab(paste0(ylab))+ggtitle(paste(main))+
      stat_bin_hex(geom="point", aes(color=log2(..count..)), binwidth=c(.05,.05), size=2, shape=16, stroke=0)+
      scale_color_viridis(option="inferno",direction = 1)+
      geom_hline(yintercept = 0,col="black")+#geom_vline(xintercept = 0,col="black")+
      geom_smooth(method='lm',formula=y~x,col="red")+
      stat_cor(method = "pearson",size=5,col="red",label.sep = "\n",
               label.x.npc = "right", label.y.npc = "bottom")+gg+
      scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))
    return(q)
  }
  myplotfun1 <- function(l2fc,xlab,ylab,main,x,y){
    q <-ggplot(l2fc, aes(x=l2fc[,x], y=l2fc[,y],col=as.factor(rp))) +
      xlab(paste0(xlab))+ylab(paste0(ylab))+ggtitle(paste(main))+
      geom_point()+
      geom_hline(yintercept = 0,col="black")+#geom_vline(xintercept = 0,col="black")+
      geom_smooth(method='lm',formula=y~x,col="black")+
      stat_cor(method = "pearson",size=5,col="black",label.sep = "\n",
               label.x.npc = "right", label.y.npc = "bottom")+gg+
      scale_x_continuous(expand = c(0,0))+scale_y_continuous(expand = c(0,0))+
      theme(legend.text = element_text(size = 14,colour = "black"),legend.title = element_blank(),
            legend.key = element_blank(),legend.background = element_blank(),
            legend.position = c(0.5,1), 
            legend.justification = c(1, 1))
    return(q)
  }
  
  ### with respect to shift
  bs <- "gene_width"
  qs <- c("shift_H3K36me2 (R2)", 
          "shift_H3K36me3 (R1)", "shift_H3K4me3 (R2)", "shift_Ser2-P (R2)", 
          "shift_Ser5-P (R2)", "shift_H3 (R2)", "shift_H3 (R1)")
  pl.list <- list()
  pl.listd <- list()
  for(i in 1:length(qs)){
    xlab="gene length, log2"
    ylab=paste0(gsub("shift_", "",qs[i]),", log2")
    x=as.character(bs)
    y=as.character(qs[i])
    main=gsub("shift_","",as.character(qs[i]))
    cat(main,"\n")
    pl.list[[paste0(main)]] <- print(myplotfun(l2fc = pl,xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
    pl.listd[[paste0(main)]] <- print(myplotfun(l2fc = pl[pl$is_icg>0,],xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
    rm(xlab,ylab,main,y,x)
  }
  
  bs <- "intron_width"
  pl.listc <- list()
  for(i in 1:length(qs)){
    xlab="intron length, log2"
    ylab=paste0(gsub("shift_", "",qs[i]),", log2")
    x=as.character(bs)
    y=as.character(qs[i])
    main=gsub("shift_","",as.character(qs[i]))
    cat(main,"\n")
    pl.listc[[paste0(main)]] <- print(myplotfun1(l2fc = gl[gl$feature=="intron1",],xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
    rm(xlab,ylab,main,y,x)
  }
  pdf(file = paste0(OUTDIR,"GeneLengths Vs Shift changes.pdf"),width = 16,height = 8)
  grid.arrange(grobs=pl.list,ncol=4,top= textGrob(paste0("gene length v/s shift changes (all genes)"),gp=gpar(fontsize=20,font=2)))
  grid.arrange(grobs=pl.listd,ncol=4,top= textGrob(paste0("gene length v/s shift changes (ICGs)"),gp=gpar(fontsize=20,font=2)))
  grid.arrange(grobs=pl.listc,ncol=4,top= textGrob(paste0("intron length v/s shift changes (ICGs)"),gp=gpar(fontsize=20,font=2)))
  dev.off()  
  
  ### with respect to l2fc
  bs <- "gene_width"
  qs <- c("l2fc_H3 (R2)", "l2fc_H3 (R1)", "l2fc_H3K36me2 (R2)", 
          "l2fc_H3K36me3 (R1)", "l2fc_H3K4me3 (R2)", "l2fc_Set2 (R1)", 
          "l2fc_Spn1 (R2)", "l2fc_Spt6 (R1)", "l2fc_RNAPII (R2)", 
          "l2fc_RNAPII (R1)", "l2fc_Ser2-P (R2)", "l2fc_Ser5-P (R2)")
  pl.list <- list()
  pl.listd <- list()
  for(i in 1:length(qs)){
    xlab="gene length, log2"
    ylab=paste0(gsub("l2fc_", "",qs[i]),", log2")
    x=as.character(bs)
    y=as.character(qs[i])
    main=gsub("l2fc_","",as.character(qs[i]))
    cat(main,"\n")
    pl.list[[paste0(main)]] <- print(myplotfun(l2fc = pl,xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
    pl.listd[[paste0(main)]] <- print(myplotfun(l2fc = pl[pl$is_icg>0,],xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
    rm(xlab,ylab,main,y,x)
  }
  
  
  bs <- "intron_width"
  pl.listc <- list()
  for(i in 1:length(qs)){
    xlab="intron length, log2"
    ylab=paste0(gsub("l2fcIE_", "",qs[i]),", log2")
    x=as.character(bs)
    y=as.character(qs[i])
    main=gsub("l2fcIE_","",as.character(qs[i]))
    cat(main,"\n")
    pl.listc[[paste0(main)]] <- print(myplotfun1(l2fc = gl[gl$feature=="intron1",],xlab=xlab,ylab=ylab,main=main,x=x,y=y) )
    rm(xlab,ylab,main,y,x)
  }
  pdf(file = paste0(OUTDIR,"IntronLengths Vs IntronOccup changes.pdf"),width = 16,height = 12)
  grid.arrange(grobs=pl.list,ncol=4,top= textGrob(paste0("gene length v/s gene occupancy change (all genes)"),gp=gpar(fontsize=20,font=2)))
  grid.arrange(grobs=pl.listd,ncol=4,top= textGrob(paste0("gene length v/s gene occupancy change (ICGs)"),gp=gpar(fontsize=20,font=2)))
  grid.arrange(grobs=pl.listc,ncol=4,top= textGrob(paste0("intron length v/s intron occupancy change (ICGs)"),gp=gpar(fontsize=20,font=2)))
  dev.off()  
  
  
  ########################################
  load(file="A:/work/WinstonLab/Natalia_Reim/MayNov17/gr/Gene groups.RData")
  rp <- rp[rp$gene!="RPP1",]
  rp <- rp[-grep("MRP",rp$gene),]
  pl$rp <- ifelse(pl$tracking_id%in%rp$tracking_id,1,0)
  pl$intron <- ifelse(pl$tracking_id%in%introns$tracking_id,1,0)
  sum(pl$rp)
  sum(pl$intron)
  
  
  qw <- pl[,c("tracking_id","basal_rna_exp","rp","intron")]
  qw$rpicg <- ifelse(qw$rp==1 & qw$intron==1,1,0)
  qw$nonrpicg <- ifelse(qw$rp==0 & qw$intron==1,1,0)
  qw$rpnonicg <- ifelse(qw$rp==1 & qw$intron==0,1,0)
  qw <- qw[rowSums(qw[,5:7])>0,]
  
  sum(qw$rpicg)
  
  qw <- rbind(cbind(qw[qw$rpicg==1,c("tracking_id","basal_rna_exp")], class="RP-ICG (89)"),
        cbind(qw[qw$nonrpicg==1,c("tracking_id","basal_rna_exp")], class="nonRP-ICG (151)"),
        cbind(qw[qw$rpnonicg==1,c("tracking_id","basal_rna_exp")], class="RP-nonICG (48)"))

  p <- ggplot(qw, aes(x=class,y=basal_rna_exp,fill=class))
  p <- p + geom_violin(color="black")
  p <- p + xlab("") + ylab("sense expression, log2")
  p <- p + geom_hline(yintercept = seq(0,14,2),col="gray",lty=2)
  p <- p + scale_color_brewer(palette = "Dark2")
  p <- p + geom_boxplot(width=0.05,fill='white', color="black",outlier.shape = NA)
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
  
  file=paste0(OUTDIR,"Sense expression for ICGs.pdf")
  pdf(file,width = 6,height = 5)
  p
  dev.off()
  
}



########################################### Old correlation script, can be deleted!!
if(T){
  
  
  ### (1) Shifts v/s expression changes
  pdf(file = paste0(OUTDIR,"PF_ShiftsVsExpressionChanges.pdf"),width = 12,height = 8)
  grid.arrange( myplotfun(l2fc,xlab="log2FC, shift",main= "H3K36me2",x="shift_H3K36me2",y="log2FC, expression"),
                myplotfun(l2fc,xlab="log2FC, shift",main="H3K36me3",x="shift_H3K36me3",y="log2FC, expression"),
                myplotfun(l2fc,xlab="log2FC, shift",main="H3K4me3",x="shift_H3K4me3",y="log2FC, expression"),
                myplotfun(l2fc,xlab="log2FC, shift",main="Ser2-P",x="shift_S2P",y="log2FC, expression"),
                myplotfun(l2fc,xlab="log2FC, shift",main="Ser5-P",x="shift_S5P",y="log2FC, expression"),
                myplotfun(l2fc,xlab="log2FC, shift",main="Spt6",x="shift_Spt6",y="log2FC, expression"),
                ncol=3,top= textGrob(paste0("(relative) shift ratio v/s expression changes"),gp=gpar(fontsize=20,font=2)))
  dev.off() 
  
  ## 5% fdr
  pdf(file = paste0(OUTDIR,"PF_ShiftsVsExpressionChanges_DEgenes.pdf"),width = 12,height = 8)
  grid.arrange( myplotfun(l2fc[l2fc$padj<0.05,],xlab="log2FC, shift",main= "H3K36me2",x="shift_H3K36me2",y="log2FC, expression"),
                myplotfun(l2fc[l2fc$padj<0.05,],xlab="log2FC, shift",main="H3K36me3",x="shift_H3K36me3",y="log2FC, expression"),
                myplotfun(l2fc[l2fc$padj<0.05,],xlab="log2FC, shift",main="H3K4me3",x="shift_H3K4me3",y="log2FC, expression"),
                myplotfun(l2fc[l2fc$padj<0.05,],xlab="log2FC, shift",main="Ser2-P",x="shift_S2P",y="log2FC, expression"),
                myplotfun(l2fc[l2fc$padj<0.05,],xlab="log2FC, shift",main="Ser5-P",x="shift_S5P",y="log2FC, expression"),
                myplotfun(l2fc[l2fc$padj<0.05,],xlab="log2FC, shift",main="Spt6",x="shift_Spt6",y="log2FC, expression"),
                ncol=3,top= textGrob(paste0("(relative) shift ratio v/s differential expression changes"),gp=gpar(fontsize=20,font=2)))
  dev.off() 
  ## 10% fdr
  pdf(file = paste0(OUTDIR,"PF_ShiftsVsExpressionChanges_DEgenes_10pFDR.pdf"),width = 12,height = 8)
  grid.arrange( myplotfun(l2fc[l2fc$padj<0.1,],xlab="log2FC, shift",main= "H3K36me2",x="shift_H3K36me2",y="log2FC, expression"),
                myplotfun(l2fc[l2fc$padj<0.1,],xlab="log2FC, shift",main="H3K36me3",x="shift_H3K36me3",y="log2FC, expression"),
                myplotfun(l2fc[l2fc$padj<0.1,],xlab="log2FC, shift",main="H3K4me3",x="shift_H3K4me3",y="log2FC, expression"),
                myplotfun(l2fc[l2fc$padj<0.1,],xlab="log2FC, shift",main="Ser2-P",x="shift_S2P",y="log2FC, expression"),
                myplotfun(l2fc[l2fc$padj<0.1,],xlab="log2FC, shift",main="Ser5-P",x="shift_S5P",y="log2FC, expression"),
                myplotfun(l2fc[l2fc$padj<0.1,],xlab="log2FC, shift",main="Spt6",x="shift_Spt6",y="log2FC, expression"),
                ncol=3,top= textGrob(paste0("(relative) shift ratio v/s differential expression changes"),gp=gpar(fontsize=20,font=2)))
  dev.off() 
  
  ### (2) Factor occupancy v/s expression changes
  plist <- list()
  for(i in 11:22){
    plist[[ paste0(names(l2fc)[i]) ]] <- myplotfun(l2fc,xlab="log2FC, occupancy",main= names(l2fc)[i],x=names(l2fc)[i],y="log2FC, expression")
  }
  
  pdf(file = paste0(OUTDIR,"PF_OccupancyVsExpressionChanges.pdf"),width = 16,height = 12)
  grid.arrange(grobs=plist,ncol=4,top= textGrob(paste0("factor occupancy v/s expression changes"),gp=gpar(fontsize=20,font=2)))
  dev.off() 
  
  
  ### (3) Pairwise factor shift comparisons
  m <- l2fc[,2:7]
  names(m) <- gsub("shift_","",names(m))
  q <- ggscat(m, columns = 1:ncol(m),color = NULL,size=10)+
    stat_bin_hex(geom="point", aes(color=log2(..count..)), binwidth=c(.05,.05), size=2, shape=16, stroke=0)+
    scale_color_viridis(option="inferno",direction = 1)+
    theme(axis.text.x = element_text(size=14,colour = "black"),
          axis.text.y = element_text(size=14,colour = "black"),
          axis.title =  element_text(size=16,color="black"),
          strip.text = element_text(color="black",size=16))+
    theme(legend.position = "none")+
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid = element_blank())+
    geom_smooth(method='lm',formula=y~x,col="red")+
    xlab("Log2FC, absolute shift")+ylab("Log2FC, absolute shift")+
    ggtitle("Pairwise Shift change comparison")+
    theme(plot.title = element_text(colour ="black",size = 18,hjust = 0.5))
  
  pdf(file = paste0(OUTDIR,"PF_PairwiseShiftComparison.pdf"),width = 10,height = 10)
  q
  dev.off()
  
  
  ### (4) Pairwise factor Occupancy comparisons
  m <- l2fc[,11:22]
  q <- ggscat(m, columns = 1:ncol(m),color = NULL,size=10)+
    stat_bin_hex(geom="point", aes(color=log2(..count..)), binwidth=c(.05,.05), size=2, shape=16, stroke=0)+
    scale_color_viridis(option="inferno",direction = 1)+
    theme(axis.text.x = element_text(size=14,colour = "black"),
          axis.text.y = element_text(size=14,colour = "black"),
          axis.title =  element_text(size=16,color="black"),
          strip.text = element_text(color="black",size=16))+
    theme(legend.position = "none")+
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          panel.grid = element_blank())+
    geom_smooth(method='lm',formula=y~x,col="red")+
    xlab("Log2FC, absolute occupancy")+ylab("Log2FC, absolute occupancy")+
    ggtitle("Pairwise Occupancy change comparison")+
    theme(plot.title = element_text(colour ="black",size = 18,hjust = 0.5))
  
  pdf(file = paste0(OUTDIR,"PF_PairwiseOccupancyComparison.pdf"),width = 16,height = 16)
  q
  dev.off()
  
  
  
  ### (5) Shift v/s selected factor occupancy changes
  m <- merge(sh,l2fc,by="tracking_id",all.y=T)
  myplotfun <- function(l2fc,xlab,main,x,y){
    q <-ggplot(l2fc, aes(x=l2fc[,x], y=l2fc[,y])) +
      stat_bin_hex(geom="point", aes(color=log2(..count..)), binwidth=c(.05,.05), size=2, shape=16, stroke=0)+
      scale_color_viridis(option="inferno",direction = 1)+
      theme(axis.text.x = element_text(size=14,colour = "black"),
            axis.text.y = element_text(size=14,colour = "black"),
            axis.title =  element_text(size=16,color="black"))+
      theme(legend.position = "none")+
      theme(panel.background = element_rect(fill = "white", colour = "black"),
            panel.grid = element_blank())+
      geom_hline(yintercept = 0,col="black")+geom_vline(xintercept = 0,col="black")+
      xlab(paste0(xlab))+ylab(paste0(y))+
      geom_smooth(method='lm',formula=y~x,col="red")+
      ggtitle(paste(main))+
      theme(plot.title = element_text(colour ="black",size = 18,hjust = 0.5))+
      stat_cor(method = "pearson", label.sep = "\n",label.x = max(l2fc[,x],na.rm=T)*0.01, label.y = -2,size=5)
    return(q)
  }
  
  plist1 <- list()
  plist2 <- list()
  plist3 <- list()
  plist4 <- list()
  plist5 <- list()
  plist6 <- list()
  plist7 <- list()
  plist8 <- list()
  for(i in 20:24){
    plist1[[ paste0(names(m)[i]) ]] <- myplotfun(m,xlab="log2FC, occupancy",main= names(m)[i],x=names(m)[i],y="shift_H3")
    plist2[[ paste0(names(m)[i]) ]] <- myplotfun(m,xlab="log2FC, occupancy",main= names(m)[i],x=names(m)[i],y="shift_H3_r1")
    plist3[[ paste0(names(m)[i]) ]] <- myplotfun(m,xlab="log2FC, occupancy",main= names(m)[i],x=names(m)[i],y="shift_H3K36me2")
    plist4[[ paste0(names(m)[i]) ]] <- myplotfun(m,xlab="log2FC, occupancy",main= names(m)[i],x=names(m)[i],y="shift_H3K36me3")
    plist5[[ paste0(names(m)[i]) ]] <- myplotfun(m,xlab="log2FC, occupancy",main= names(m)[i],x=names(m)[i],y="shift_H3K4me3")
    plist6[[ paste0(names(m)[i]) ]] <- myplotfun(m,xlab="log2FC, occupancy",main= names(m)[i],x=names(m)[i],y="shift_S2P")
    plist7[[ paste0(names(m)[i]) ]] <- myplotfun(m,xlab="log2FC, occupancy",main= names(m)[i],x=names(m)[i],y="shift_S5P")
    plist8[[ paste0(names(m)[i]) ]] <- myplotfun(m,xlab="log2FC, occupancy",main= names(m)[i],x=names(m)[i],y="shift_Spt6")
  }
  pdf(file = paste0(OUTDIR,"PF_OccupancyVsRelativeShiftChange_SelectedFactors.pdf"),width = 20,height = 5)
  grid.arrange(grobs=plist1,ncol=5,top= textGrob(paste0("factor occupancy v/s H3 shift ratio"),gp=gpar(fontsize=20,font=2)))
  grid.arrange(grobs=plist2,ncol=5,top= textGrob(paste0("factor occupancy v/s H3 (round1) shift ratio"),gp=gpar(fontsize=20,font=2)))
  grid.arrange(grobs=plist3,ncol=5,top= textGrob(paste0("factor occupancy v/s H3K36me2 shift ratio"),gp=gpar(fontsize=20,font=2)))
  grid.arrange(grobs=plist4,ncol=5,top= textGrob(paste0("factor occupancy v/s H3K36me3 shift ratio"),gp=gpar(fontsize=20,font=2)))
  grid.arrange(grobs=plist5,ncol=5,top= textGrob(paste0("factor occupancy v/s H3K4me3 shift ratio"),gp=gpar(fontsize=20,font=2)))
  grid.arrange(grobs=plist6,ncol=5,top= textGrob(paste0("factor occupancy v/s Ser2-P shift ratio"),gp=gpar(fontsize=20,font=2)))
  grid.arrange(grobs=plist7,ncol=5,top= textGrob(paste0("factor occupancy v/s Ser5-P shift ratio"),gp=gpar(fontsize=20,font=2)))
  grid.arrange(grobs=plist8,ncol=5,top= textGrob(paste0("factor occupancy v/s Ser5-P shift ratio"),gp=gpar(fontsize=20,font=2)))
  dev.off() 
  
}



