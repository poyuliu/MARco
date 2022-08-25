#Basic 01
aggregate2df <- function(data,groups,FUN){
  agg.data <- aggregate(data,list(groups),FUN)
  rownames(agg.data) <- agg.data[,1]
  agg.data <- as.data.frame(t(agg.data[,-1]))
  return(agg.data)
}


# > prevalence(data,prevalence)
# filter feature table by prevalence (in 0-1 decimal)
prevalence <- function(data,prevalence=0.1){
  prev <- apply(data,1,function(x) sum(x>0)/ncol(data))
  data <- data[prev>=prevalence,]
  return(data)
}

# Betadiversity 01
ordNtest <- function(data,group,dist=NULL){
  if(is.null(dist)){
    bc.dist <- vegan::vegdist(t(data))
  } else if(!is.null(dist)){
    bc.dist <- dist
  }
  pcoa <- cmdscale(bc.dist,eig = T)
  set.seed(123); permanova <- adonisx(bc.dist~group)
  return(list(pcoa=pcoa,permanova=permanova))
}


# Betadicersity 02
plot.ordi <- function (ordN, group, colset, ordiellipse = FALSE, ordispider = FALSE, 
                       plot = TRUE, pch = 19, legend = NULL, test = FALSE, ordilabel = NULL,
                       right_frame=10.1,inset=c(-0.3,-0.2)) 
{
  ord <- ordN$pcoa
  pct <- round(prop.table(ord$eig[ord$eig > 0]) * 100, 2)
  par(mar = c(5.1, 4.1, 4.1, right_frame))
  if (length(pct) > 0) {
    plot(vegan::scores(ord), type = "n", las = 1, xlab = paste0("PCoA1 (", 
                                                                round(pct[1], 2), "%)"), ylab = paste0("PCoA2 (", 
                                                                                                       round(pct[2], 2), "%)"))
  }
  else if (length(pct) == 0) 
    plot(vegan::scores(ord), type = "n", las = 1, xlab = "Dim1", 
         ylab = "Dim2")
  if (ordiellipse == TRUE) 
    vegan::ordiellipse(ord, group, col = colset, lwd = 0.5)
  if (ordispider == TRUE) 
    vegan::ordispider(ord, group, col = colset, lwd = 0.5)
  if (plot == TRUE & length(pch) > 1) {
    points(vegan::scores(ord), pch = pch[as.numeric(group)], 
           cex = 1, col = paste0(colset, "aa")[as.numeric(group)])
  }
  else if (plot == TRUE & length(pch) == 1) 
    points(vegan::scores(ord), pch = pch, cex = 1, col = paste0(colset, 
                                                                "aa")[as.numeric(group)])
  if (!is.null(legend)) 
    legend("topright", legend = legend, pch = pch, col = paste0(colset, 
                                                                "aa"), bty = "n", cex = 0.8, xpd = T, inset = c(inset[1], 
                                                                                                                0))
  if (test == TRUE) {
    p <- format(ordN$permanova$aov.tab$`Pr(>F)`[1], digits = 2)
    r2 <- format(ordN$permanova$aov.tab$R2[1], digits = 2)
    legend("bottomright", legend = bquote(R^2 == .(r2)), 
           bty = "n", cex = 0.8, xpd = T, inset = c(inset[2], 0.06))
    legend("bottomright", legend = bquote(italic(p == .(p))), 
           bty = "n", cex = 0.8, xpd = T, inset = c(inset[2], 0))
  }
  if (!is.null(ordilabel)) {
    vegan::ordilabel(scores(ordN$pcoa), cex = 0.7, fill = "#cacaca55", 
                     border = 1, col = "#1b262c", select = which(rownames(scores(ordN$pcoa)) %in% 
                                                                   ordilabel))
  }
}

# Betadiversity 03-1
ordEnvfit <- function(data,ENV,...){
  bc.dist <- vegan::vegdist(t(data))
  pcoa <- cmdscale(bc.dist,eig = T)
  fit <-  vegan::envfit(pcoa,ENV,...)
  return(list(pcoa=pcoa,fit=fit))
}

# Betadicersity 03-2
plot.ordi.env <- function(ordN,select,colset,plot=TRUE,pch=19,ordilabel=NULL,fit.p=1,legend=TRUE,legend.labels=NULL,rank=TRUE,select.round=NULL){
  ord <- ordN$pcoa
  pct <- round(prop.table(ord$eig[ord$eig>0])*100,2)
  fit.corrd <- ordN$fit$vectors$arrows*ordN$fit$vectors$r
  select.raw <- select
  if(rank==TRUE){
    if(is.null(select.round)){
      select <- ceiling(rank(select))
    } else select <- ceiling(rank(round(select,select.round)))
  } else select <- as.factor(select)
  par(mar=c(5.1,4.1,4.1,10.1),bty='n',xaxt='n',yaxt='n')
  plot(vegan::scores(ord),type="n",las=1,
       xlab=paste0("PCoA1 (",round(pct[1],2),"%)"),ylab=paste0("PCoA2 (",round(pct[2],2),"%)"))
  abline(h=0,v=0,lty=3)
  if(plot==TRUE & length(pch)>1){
    points(vegan::scores(ord),pch=pch[as.numeric(select)],cex=1,col=paste0(colset,"aa")[as.numeric(select)])
  } else if(plot==TRUE & length(pch)==1) points(vegan::scores(ord),pch=pch,cex=1,col=paste0(colset,"aa")[as.numeric(select)])
  arrows(x0 = 0,y0 = 0,x1 = fit.corrd[which(ordN$fit$vectors$pvals < fit.p),1],y1 = fit.corrd[which(ordN$fit$vectors$pvals < fit.p),2],length = 0)
  vegan::ordilabel(fit.corrd[which(ordN$fit$vectors$pvals < fit.p),],labels = rownames(fit.corrd)[which(ordN$fit$vectors$pvals < fit.p)],cex = 0.7,fill = "#cacaca55",border = 1,col = "#1b262c")
  if(!is.null(ordilabel)){
    vegan::ordilabel(vegan::scores(ordN$pcoa),cex = 0.7,fill = "#cacaca55",border = 1,col = "#1b262c",
                     select = which(rownames(vegan::scores(ordN$pcoa)) %in% ordilabel))
  }
  if(legend==TRUE & rank==TRUE){
    par(mar=c(15.1,40.1,4.1,5.1),new=T)
    image(t(seq_len(length(unique(sort(round(select,1)))))),col=paste0(colset[unique(ceiling(rank(select))[order(round(select,1))])],"aa"))
    text(x=1,y=seq(0,1,length.out = length(unique(sort(round(select.raw,1))))),labels = unique(sort(round(select.raw,1))),cex=0.7,pos=4,xpd=T)
    text(x = 0.5,y = 1.1,labels = legend.labels,font = 2,cex=0.8,adj = 0.8,xpd=T)
  } else if(legend==TRUE & rank==FALSE){
    legend("topright",legend =levels(select),title = legend.labels ,pch=pch,col=paste0(colset,"aa"),bty="n",cex=0.8,xpd=T,inset = c(-0.3,0))
  }
}

# Betadiversity 04-1
dbrdaNtest <- function(data,ENV){
  dbRDA <- vegan::capscale(t(data)~.,data = ENV,distance = "bray")
  set.seed(123); permucapscale <- vegan::anova.cca(dbRDA,by="terms")
  return(list(dbRDA=dbRDA,permucapscale=permucapscale))
}

# Betadicersity 04-2
plot.dbRDA <- function(ordN,select,colset,plot=TRUE,pch=19,ordilabel=NULL,fit.p=1,legend=TRUE,legend.labels=NULL,rank=TRUE,select.round=NULL){
  ord <- ordN$dbRDA
  eig <- c(ordN$dbRDA$CCA$eig,ordN$dbRDA$CA$eig)
  pct <- round(prop.table(eig[eig>0])*100,2)
  fit.corrd <- ordN$dbRDA$CCA$biplot
  select.raw <- select
  if(rank==TRUE){
    if(is.null(select.round)){
      select <- ceiling(rank(select))
    } else select <- ceiling(rank(round(select,select.round)))
  } else select <- as.factor(select)
  par(mar=c(5.1,4.1,4.1,10.1),bty='n',xaxt='n',yaxt='n')
  plot(vegan::scores(ordN$dbRDA,display = "sites"),type="n",las=1,
       xlab=paste0("CAP1 (",round(pct[1],2),"%)"),ylab=paste0("CAP2 (",round(pct[2],2),"%)"))
  abline(h=0,v=0,lty=3)
  if(plot==TRUE & length(pch)>1){
    points(vegan::scores(ordN$dbRDA,display = "sites"),pch=pch[as.numeric(select)],cex=1,
           col=paste0(colset,"aa")[as.numeric(select)])
  } else if(plot==TRUE & length(pch)==1) points(vegan::scores(ordN$dbRDA,display = "sites"),pch=pch,cex=1,col=paste0(colset,"aa")[as.numeric(select)])
  arrows(x0 = 0,y0 = 0,x1 = fit.corrd[which(ordN$permucapscale$`Pr(>F)` < fit.p),1],y1 = fit.corrd[which(ordN$permucapscale$`Pr(>F)` < fit.p),2],length = 0)
  vegan::ordilabel(fit.corrd[which(ordN$permucapscale$`Pr(>F)` < fit.p),],labels = rownames(fit.corrd)[which(ordN$permucapscale$`Pr(>F)` < fit.p)],cex = 0.7,fill = "#cacaca55",border = 1,col = "#1b262c")
  if(!is.null(ordilabel)){
    vegan::ordilabel(vegan::scores(ordN$dbRDA,display = "sites"),cex = 0.7,fill = "#cacaca55",border = 1,col = "#1b262c",
                     select = which(rownames(vegan::scores(ordN$dbRDA,display = "sites")) %in% ordilabel))
  }
  if(legend==TRUE & rank==TRUE){
    par(mar=c(15.1,40.1,4.1,5.1),new=T)
    image(t(seq_len(length(unique(sort(round(select,1)))))),col=paste0(colset[unique(ceiling(rank(select))[order(round(select,1))])],"aa"))
    text(x=1,y=seq(0,1,length.out = length(unique(sort(round(select.raw,1))))),labels = unique(sort(round(select.raw,1))),cex=0.7,pos=4,xpd=T)
    text(x = 0.5,y = 1.1,labels = legend.labels,font = 2,cex=0.8,adj = 0.8,xpd=T)
  } else if(legend==TRUE & rank==FALSE){
    legend("topright",legend =levels(select),title = legend.labels ,pch=pch,col=paste0(colset,"aa"),bty="n",cex=0.8,xpd=T,inset = c(-0.3,0))
  }
}

# Alpha diversity 01
alphdiv <- function(data,group,test=FALSE){
  alpha.0 <- colSums(data>0) # observed OTUs
  alpha.1 <- vegan::diversity(data,index = "shannon",MARGIN = 2)
  alpha.2 <- vegan::diversity(data,index = "simpson",MARGIN = 2)
  alpha.3 <- apply(data,2,fossil::chao1)
  if(test==TRUE & !missing(group)){
    if(length(group)==2){
      test.0 <- wilcox.test(alpha.0~group)
      test.1 <- wilcox.test(alpha.1~group)
      test.2 <- wilcox.test(alpha.2~group)
      test.3 <- wilcox.test(alpha.3~group)
    } else if(length(group)>2){
      test.0 <- kruskal.test(alpha.0~group)
      test.1 <- kruskal.test(alpha.1~group)
      test.2 <- kruskal.test(alpha.2~group)
      test.3 <- kruskal.test(alpha.3~group)
    }
    return(list(obsOTU=alpha.0,shannon=alpha.1,simpson=alpha.2,chao1=alpha.3,
                obsOTU.stat=test.0,shannon.stat=test.1,simpson.stat=test.2,chao1.stat=test.3))
  } else if(test==TRUE & missing(group)){
    warning("'group is missing!!!'")
  } else return(list(obsOTU=alpha.0,shannon=alpha.1,simpson=alpha.2,chao1=alpha.3))
  
}

# Alpha diversity 02
jitter.points <- function (data, groups, colset) 
{
  x <- list()
  for (i in 1:length(groups)) {
    #set.seed(123)
    x[[i]] <- jitter(rep(1, sum(as.numeric(groups) == i)), 8)
  }
  points(x[[1]], data[as.numeric(groups) == 1], pch = 21, bg = colset[1], cex = 1.1)
  for(i in 2:length(groups)){
    points((x[[i]]+(i-1)), data[as.numeric(groups) == i], pch = 21, bg = colset[i], cex = 1.1)
  }
}


# Alpha diversity 03
plot.alpha <- function (alpha, group, colset, names = NULL, Richness = TRUE, 
                        Shannon = TRUE, Simpson = TRUE, Chao1 = TRUE, test = FALSE,...) 
{
  par(bty = "L")
  if (Richness == TRUE) {
    boxplot(alpha$obsOTU ~ group, ..., xlab = "", ylab = "Richness", 
            col = paste0(colset, "99"), names = names, pch = NA)
    jitter.points(alpha$obsOTU, group, colset)
    if (test == TRUE & !is.null(alpha$obsOTU.stat)) {
      legend("topright", legend = paste("p ", ifelse(alpha$obsOTU.stat$p.value >= 
                                                       0.001, paste("=", round(alpha$obsOTU.stat$p.value, 
                                                                               3)), "< 0.001")), bty = "n", xpd = T, inset = -0.05)
    }
    else if (test == TRUE & is.null(alpha$obsOTU.stat)) 
      warning("Alpha lacks testings!")
  }
  if (Shannon == TRUE) {
    boxplot(alpha$shannon ~ group, ..., xlab = "", ylab = "Shannon's index", 
            col = paste0(colset, "99"), names = names, pch = NA)
    jitter.points(alpha$shannon, group, colset)
    if (test == TRUE & !is.null(alpha$shannon.stat)) {
      legend("topright", legend = paste("p ", ifelse(alpha$shannon.stat$p.value >= 
                                                       0.001, paste("=", round(alpha$shannon.stat$p.value, 
                                                                               3)), "< 0.001")), bty = "n", xpd = T, inset = -0.05)
    }
    else if (test == TRUE & is.null(alpha$shannon.stat)) 
      warning("Alpha lacks testings!")
  }
  if (Simpson == TRUE) {
    boxplot(alpha$simpson ~ group, ..., xlab = "", ylab = "Simpson's index", 
            col = paste0(colset, "99"), names = names, pch = NA)
    jitter.points(alpha$simpson, group, colset)
    if (test == TRUE & !is.null(alpha$simpson.stat)) {
      legend("topright", legend = paste("p ", ifelse(alpha$simpson.stat$p.value >= 
                                                       0.001, paste("=", round(alpha$simpson.stat$p.value, 
                                                                               3)), "< 0.001")), bty = "n", xpd = T, inset = -0.05)
    }
    else if (test == TRUE & is.null(alpha$simpson.stat)) 
      warning("Alpha lacks testings!")
  }
  if (Chao1 == TRUE) {
    boxplot(alpha$chao1 ~ group, ..., xlab = "", ylab = "Chao1 index", 
            col = paste0(colset, "99"), names = names, pch = NA)
    jitter.points(alpha$chao1, group, colset)
    if (test == TRUE & !is.null(alpha$chao1.stat)) {
      legend("topright", legend = paste("p ", ifelse(alpha$chao1.stat$p.value >= 
                                                       0.001, paste("=", round(alpha$chao1.stat$p.value, 
                                                                               3)), "< 0.001")), bty = "n", xpd = T, inset = -0.05)
    }
    else if (test == TRUE & is.null(alpha$chao1.stat)) 
      warning("Alpha lacks testings!")
  }
}

# Alpha diversity 04
vioplot.alpha <- function(alpha,group,colset,names=NULL,Richness=TRUE,Shannon=TRUE,Simpson=TRUE,Chao1=TRUE,test=FALSE,...){
  if(Richness==TRUE) {
    vioplot::vioplot(alpha$obsOTU~group,...,xlab="",ylab="Richness",col=paste0(colset,"99"),names=names)
    jitter.points(alpha$obsOTU,group,colset)
    if(test==TRUE & !is.null(alpha$obsOTU.stat)){
      legend("topright",
             legend = paste("p ",ifelse(alpha$obsOTU.stat$p.value>=0.001,paste("=",round(alpha$obsOTU.stat$p.value,3)),"< 0.001")),
             bty="n",xpd=T,inset = -0.05)
    } else if(test==TRUE & is.null(alpha$obsOTU.stat)) warning("Alpha lacks testings!")
  }
  if(Shannon==TRUE) {
    vioplot::vioplot(alpha$shannon~group,...,xlab="",ylab="Shannon's index",col=paste0(colset,"99"),names=names)
    jitter.points(alpha$shannon,group,colset)
    if(test==TRUE & !is.null(alpha$shannon.stat)){
      legend("topright",
             legend = paste("p ",ifelse(alpha$shannon.stat$p.value>=0.001,paste("=",round(alpha$shannon.stat$p.value,3)),"< 0.001")),
             bty="n",xpd=T,inset = -0.05)
    } else if(test==TRUE & is.null(alpha$shannon.stat)) warning("Alpha lacks testings!")
  }
  if(Simpson==TRUE) {
    vioplot::vioplot(alpha$simpson~group,...,xlab="",ylab="Simpson's index",col=paste0(colset,"99"),names=names)
    jitter.points(alpha$simpson,group,colset)
    if(test==TRUE & !is.null(alpha$simpson.stat)){
      legend("topright",
             legend = paste("p ",ifelse(alpha$simpson.stat$p.value>=0.001,paste("=",round(alpha$simpson.stat$p.value,3)),"< 0.001")),
             bty="n",xpd=T,inset = -0.05)
    } else if(test==TRUE & is.null(alpha$simpson.stat)) warning("Alpha lacks testings!")
  }
  if(Chao1==TRUE) {
    vioplot::vioplot(alpha$chao1~group,...,xlab="",ylab="Chao1 index",col=paste0(colset,"99"),names=names)
    jitter.points(alpha$chao1,group,colset)
    if(test==TRUE & !is.null(alpha$chao1.stat)){
      legend("topright",
             legend = paste("p ",ifelse(alpha$chao1.stat$p.value>=0.001,paste("=",round(alpha$chao1.stat$p.value,3)),"< 0.001")),
             bty="n",xpd=T,inset = -0.05)
    } else if(test==TRUE & is.null(alpha$chao1.stat)) warning("Alpha lacks testings!")
  }
  
}

# DESeq2
DESq2 <- function(data,group){
  require(DESeq2)
  cnts <- as.matrix(data)
  cnts <- cnts[rowSums(cnts)!=0,]
  colnames(cnts) <- NULL
  cond <- group
  # object construction
  dds <- DESeq2::DESeqDataSetFromMatrix(cnts, DataFrame(cond), ~ cond)
  dds <- DESeq2::estimateSizeFactors(dds,type="poscounts")
  # standard analysis
  dds <- DESeq2::DESeq(dds)
  res <- as.data.frame(DESeq2::results(dds))
  return(res)
}

# DESeq2 volcano plot
volcano.plot <- function(DESeq.res,colset,p.adj=TRUE,alpha=0.05,...){
  FC2 <- DESeq.res$log2FoldChange
  if(p.adj==TRUE){
    PADJ <- DESeq.res$padj
  
    plot(FC2,-log10(PADJ),
         type="n",xlab=expression(Log[2]~"fold change"),ylab=expression(-Log[10]~"adjusted P"))
    grid()
    abline(v=c(1,-1),h=-log10(alpha),...)
    colX <- colset[as.numeric(FC2>0)+1]
    colX[which(PADJ >= alpha)] <- paste0(colX[which(PADJ >= alpha)],"33")
    colX[which(PADJ < alpha)] <- paste0(colX[which(PADJ < alpha)],"cc")
    points(FC2,-log10(PADJ),pch=21,col="#fcfcfc55",bg=colX,cex=0.9)
  } else if(p.adj==FALSE){
      PADJ <- DESeq.res$pvalue
  
      plot(FC2,-log10(PADJ),
          type="n",xlab=expression(Log[2]~"fold change"),ylab=expression(-Log[10]~"P-values"))
      grid()
      abline(v=c(1,-1),h=-log10(alpha),...)
      colX <- colset[as.numeric(FC2>0)+1]
      colX[which(PADJ >= alpha)] <- paste0(colX[which(PADJ >= alpha)],"33")
      colX[which(PADJ < alpha)] <- paste0(colX[which(PADJ < alpha)],"cc")
      points(FC2,-log10(PADJ),pch=21,col="#fcfcfc55",bg=colX,cex=0.9)
  }
  
  
}

# DESeq2 significant taxa table
DESeq2.sig.table <- function(DESeq.res,taxa,kable=TRUE,caption=NULL,padj=TRUE,alpha=0.05){
  FC2 <- DESeq.res$log2FoldChange
  if(padj==TRUE){
    PVAL <- DESeq.res$padj
    sig.table <- DESeq.res[which(abs(FC2)>1 & PVAL<alpha),]
    } else if(padj==FALSE){
      PVAL <- DESeq.res$pvalue
      sig.table <- DESeq.res[which(abs(FC2)>1 & PVAL<alpha),]
    }
  sig.taxa <- taxa[match(rownames(sig.table),rownames(taxa)),]
  colnames(sig.taxa) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  sig.table <- cbind(sig.table,sig.taxa)
  if(kable==TRUE){
    knitr::kable(sig.table,caption = caption)
    } else return(sig.table)
}

#2 groups heatmap
G2.heatmap <- function(otu.sig, groups,p.val,taxa.sig,otu,metadata,colset=c("#2166AC","#B2182B"),alpha=0.05,...){
  annotation_lab <- data.frame(Groups=groups)
  rownames(annotation_lab) <- metadata$SampleID[match(colnames(otu),metadata$SampleID)]
  ann_colors <- list(
    Groups=colset,
    "P.values"=colorRampPalette(rev(c("#fffdf9","#6bcf00")))(10)
  )
  names(ann_colors$Groups) <- levels(groups)
  
  g2.otu.plot <- otu.sig[order(p.val[p.val< alpha]),]
  g2.tax.plot <- taxa.sig[order(p.val[p.val< alpha]),]
  
  
  annotation_row <- data.frame("P.values"=p.val[p.val< alpha])
  annotation_row$P.values <- annotation_row$P.values[order(annotation_row$P.values)]
  rownames(annotation_row) <- rownames(g2.otu.plot)
  pheatmap::pheatmap(log1p(g2.otu.plot[,order(groups)]),scale = "row",
           annotation_col = annotation_lab,annotation_colors = ann_colors, annotation_legend = T,
           cluster_cols = F,cluster_rows = F,
           labels_row = paste(substr(g2.tax.plot[,6],1,nchar(g2.tax.plot[,6])),
                              substr(g2.tax.plot[,7],1,nchar(g2.tax.plot[,7]))),
           labels_col = colnames(g2.otu.plot[,order(groups)]),na_col = "#deeaee",border_color = "#454140",
           annotation_row = annotation_row,cellwidth = 8,cellheight = 8,fontsize = 8,
           color = colorRampPalette(c("#f7f7f7","#f7fbff","#3c4245"))(10),...)
  
}

# Taxonomy table rearrange
tax.rearrange <- function(taxonomytsv,file=TRUE,sep="; "){
  if(file==TRUE){
    fileinput <- read.delim(taxonomytsv)
    taxa <- as.character(fileinput[,2])
  } else if(file==FALSE) taxa <- as.character(taxonomytsv[,2])
  taxa <- strsplit(taxa,sep)
  for(i in 1:length(taxa)){
    if(length(taxa[[i]]) < 7){
      if(7-length(taxa[[i]]) ==1){
        taxa[[i]][7] <- "Unassigned"
      } else if(7-length(taxa[[i]]) ==2){
        taxa[[i]][7] <- "Unassigned"
        taxa[[i]][6] <- "Unassigned"
      } else if(7-length(taxa[[i]]) ==3){
        taxa[[i]][7] <- "Unassigned"
        taxa[[i]][6] <- "Unassigned"
        taxa[[i]][5] <- "Unassigned"
      } else if(7-length(taxa[[i]]) ==4){
        taxa[[i]][7] <- "Unassigned"
        taxa[[i]][6] <- "Unassigned"
        taxa[[i]][5] <- "Unassigned"
        taxa[[i]][4] <- "Unassigned"
      } else if(7-length(taxa[[i]]) ==5){
        taxa[[i]][7] <- "Unassigned"
        taxa[[i]][6] <- "Unassigned"
        taxa[[i]][5] <- "Unassigned"
        taxa[[i]][4] <- "Unassigned"
        taxa[[i]][3] <- "Unassigned"
      } else if(7-length(taxa[[i]]) ==6){
        taxa[[i]][7] <- "Unassigned"
        taxa[[i]][6] <- "Unassigned"
        taxa[[i]][5] <- "Unassigned"
        taxa[[i]][4] <- "Unassigned"
        taxa[[i]][3] <- "Unassigned"
        taxa[[i]][2] <- "Unassigned"
      } 
    }
  }
  
  taxa <- as.data.frame(do.call(rbind,taxa))
  if(file==TRUE){
    taxa <- cbind(Feature.ID=fileinput[,1],taxa)
  } else if(file==FALSE) taxa <- cbind(Feature.ID=taxonomytsv[,1],taxa)
  taxa <- apply(taxa, 2, as.character)
  rownames(taxa) <- taxa[,1]
  taxa <- taxa[,-1]
  return(taxa)
}


#summarize otu table at each level (in proportion)
taxaLvP <- function(otu,taxafromQ2tsv,level=2){
  taxa <- strsplit(as.character(taxafromQ2tsv[,2]), "; ")
  taxx <- data.frame(Level_1=NA,Level_2=NA,Level_3=NA,Level_4=NA,Level_5=NA,Level_6=NA,Level_7=NA)
  for(i in 1:length(taxa)) taxx[i,1:length(taxa[[i]])] <- taxa[[i]][1:length(taxa[[i]])]
  taxa <- taxx
  taxa.comb <- taxa
  taxa.comb$Level_1 <- paste(taxa$Level_1,sep="; ")
  taxa.comb$Level_2 <- paste(taxa$Level_1,taxa$Level_2,sep="; ")
  taxa.comb$Level_3 <- paste(taxa$Level_1,taxa$Level_2,taxa$Level_3,sep="; ")
  taxa.comb$Level_4 <- paste(taxa$Level_1,taxa$Level_2,taxa$Level_3,taxa$Level_4,sep="; ")
  taxa.comb$Level_5 <- paste(taxa$Level_1,taxa$Level_2,taxa$Level_3,taxa$Level_4,taxa$Level_5,sep="; ")
  taxa.comb$Level_6 <- paste(taxa$Level_1,taxa$Level_2,taxa$Level_3,taxa$Level_4,taxa$Level_5,taxa$Level_6,sep="; ")
  taxa.comb$Level_7 <- paste(taxa$Level_1,taxa$Level_2,taxa$Level_3,taxa$Level_4,taxa$Level_5,taxa$Level_6,taxa$Level_7,sep="; ")
  tx.level <- function(data=otu,level,percentile=T){
    level.n <- taxa.comb[,level]
    sort.sum <- aggregate(data,list(level.n),sum)
    rownames(sort.sum) <- sort.sum[,1]
    sort.sum <- sort.sum[,-1]
    if(percentile){
      sort.sum <- prop.table(as.matrix(sort.sum),2)*100
    } else sort.sum <- sort.sum
    return(sort.sum)
  }
  
  levelsum <- tx.level(data = otu,level)
  return(levelsum)
}


#Env corrrelated OTUs
EnvCorOTU <- function(data,ENV,cor.method="pearson",padj=FALSE,alpha=0.05,toprank=0.95,logdata=FALSE,logENV=FALSE,summary.plot=FALSE){
  ENV.name <- deparse(substitute(ENV))
  if(logdata==TRUE){
    data <- log1p(data)
  } else if(logdata==FALSE) data <- data
  if(logENV==TRUE){
    ENV <- log1p(ENV)
  } else if(logENV==FALSE) ENV <- ENV
  corTest <- apply(data,1,function(x) cor.test(x,ENV,method = cor.method))
  env.corr <- sapply(corTest,"[[",4)
  env.pval <- sapply(corTest,"[[",3)
  env.padj <- p.adjust(env.pval,method = "fdr")
  if(padj==TRUE){
    otu.cor.env <- data[which(env.padj < alpha & abs(env.corr) >= quantile(abs(env.corr),toprank,na.rm = T)),]
    sig.corr <- env.corr[which(env.padj < alpha & abs(env.corr) >= quantile(abs(env.corr),toprank,na.rm = T))]
    sig.p <- env.padj[which(env.padj < alpha & abs(env.corr) >= quantile(abs(env.corr),toprank,na.rm = T))]
  } else if(padj==FALSE){
    otu.cor.env <- data[which(env.pval < alpha & abs(env.corr) >= quantile(abs(env.corr),toprank,na.rm = T)),]
    sig.corr <- env.corr[which(env.pval < alpha & abs(env.corr) >= quantile(abs(env.corr),toprank,na.rm = T))]
    sig.p <- env.pval[which(env.pval < alpha & abs(env.corr) >= quantile(abs(env.corr),toprank,na.rm = T))]
  }
  if(summary.plot==TRUE){
    pseudo.col <- paste0(c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99"
                           ,"#b15928","#88B0C5","#015A96","#94C16C","#15820E","#DD7C7B","#C50000","#DFA151","#E16100","#AC94B8","#4C1F7C"
                           ,"#E1E17B","#933B0A","#6A92A7","#003C78","#76A34E","#006400","#BF5E5D","#A70000","#C18333","#C34300","#8E769A"
                           ,"#2E015E","#C3C35D","#751D00","#4C7489","#001E5A","#588530","#004600","#A1403F","#890000","#A36515","#A52500"
                           ,"#70587C","#100040","#A5A53F","#570000"),"99")
    set.seed(123);pseudo.col <- sample(pseudo.col,nrow(otu.cor.env))
    if(logdata==TRUE & logENV==TRUE){
      plot(y=seq(min(otu.cor.env),max(otu.cor.env),length.out = 10),
           x=seq(min(ENV),max(ENV),length.out = 10),
           las=1,type="n",bty="L",
           ylab="Log OTU counts",xlab=paste("Log", ENV.name))
    } else if(logdata==TRUE & logENV==FALSE){
      plot(y=seq(min(otu.cor.env),max(otu.cor.env),length.out = 10),
           x=seq(min(ENV),max(ENV),length.out = 10),
           las=1,type="n",bty="L",
           ylab="Log OTU counts",xlab=ENV.name)
    } else if(logdata==FALSE & logENV==TRUE){
      plot(y=seq(min(otu.cor.env),max(otu.cor.env),length.out = 10),
           x=seq(min(ENV),max(ENV),length.out = 10),
           las=1,type="n",bty="L",
           ylab="OTU counts",xlab=paste("Log", ENV.name))
    } else if(logdata==FALSE & logENV==FALSE){
      plot(y=seq(min(otu.cor.env),max(otu.cor.env),length.out = 10),
           x=seq(min(ENV),max(ENV),length.out = 10),
           las=1,type="n",bty="L",
           ylab="OTU counts",xlab=ENV.name)
    }
    for(i in 1:nrow(otu.cor.env)){
      abline(lm(t(otu.cor.env[i,])~ENV),col="#CACACA99",lty=3)
      points(y=t(otu.cor.env[i,]),x=ENV,pch=21,col="#CACACA99",bg=pseudo.col[i])
    }
  }
  return(list(otu.cor.env=otu.cor.env,sig.corr=sig.corr,sig.p=sig.p))
}

#PC correlated OTUs
PCENVcorOTU <- function(data,Ord,ENV,cor.method="pearson",alpha.pc.env=0.05,corr.pc.env=0.4,plot.pc.env=FALSE,
                        padj=FALSE,alpha=0.05,toprank=0.95,logdata=FALSE,logENV=FALSE,summary.plot=FALSE){
  ENV.name <- deparse(substitute(ENV))
  Ord <- Ord[[1]]
  PCs <- vegan::scores(Ord,display="sites")
  if(logENV==TRUE){
    ENV <- log1p(ENV)
  } else if(logENV==FALSE) ENV <- ENV
  select.pc.env <- apply(PCs,2,function(x) cor.test(x,ENV,method = cor.method))
  select.parm <- data.frame(corr=sapply(select.pc.env,"[[",4),pval=sapply(select.pc.env,"[[",3))
  i <- which(select.parm$pval < alpha.pc.env & abs(select.parm$corr) > corr.pc.env )
  if(length(i)==0 | length(i)>1){
    i <- which.max(abs(select.parm$corr))
  } else if(length(i)==1) i <- i
  #i <- ifelse(length(i)==0,which.max(abs(select.parm$corr)),i)
  PC <- PCs[,i]
  if(plot.pc.env==TRUE){
    plot(PC,ENV,type="n",las=1,bty="L",xlab=paste("PCoA",i),ylab=ENV.name)
    lm.out <- lm(ENV~PC)
    abline(lm.out,col="gray20")
    newx = seq(min(PC),max(PC),by = 0.05)
    conf_interval <- predict(lm.out, newdata=data.frame("PC"=newx), interval="confidence",
                             level = 0.95)
    lines(newx, conf_interval[,2], col="gray50", lty=2)
    lines(newx, conf_interval[,3], col="gray50", lty=2)
    points(PC,ENV,pch=19)
    legend("right",bty="n",legend = c(paste("Corr =",round(select.parm$corr[i],2)),
                                      paste("pval",ifelse(select.parm$pval[i]>=0.001,paste("=",round(select.parm$pval[i],3)),"< 0.001"))))
  }
  
  if(logdata==TRUE){
    data <- log1p(data)
  } else if(logdata==FALSE) data <- data
  corPC <- apply(data,1,function(x) cor.test(x,PCs[,i],method = cor.method))
  pc.corr <- sapply(corPC,"[[",4)
  pc.pval <- sapply(corPC,"[[",3)
  pc.padj <- p.adjust(pc.pval,method = "fdr")
  if(padj==TRUE){
    otu.cor.pc <- data[which(pc.padj < alpha & abs(pc.corr) >= quantile(abs(pc.corr),toprank,na.rm = T)),]
    sig.corr <- pc.corr[which(pc.padj < alpha & abs(pc.corr) >= quantile(abs(pc.corr),toprank,na.rm = T))]
    sig.p <- pc.padj[which(pc.padj < alpha & abs(pc.corr) >= quantile(abs(pc.corr),toprank,na.rm = T))]
  } else if(padj==FALSE){
    otu.cor.pc <- data[which(pc.pval < alpha & abs(pc.corr) >= quantile(abs(pc.corr),toprank,na.rm = T)),]
    sig.corr <- pc.corr[which(pc.pval < alpha & abs(pc.corr) >= quantile(abs(pc.corr),toprank,na.rm = T))]
    sig.p <- pc.pval[which(pc.pval < alpha & abs(pc.corr) >= quantile(abs(pc.corr),toprank,na.rm = T))]
  }
  if(summary.plot==TRUE){
    pseudo.col <- paste0(c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a","#ffff99"
                           ,"#b15928","#88B0C5","#015A96","#94C16C","#15820E","#DD7C7B","#C50000","#DFA151","#E16100","#AC94B8","#4C1F7C"
                           ,"#E1E17B","#933B0A","#6A92A7","#003C78","#76A34E","#006400","#BF5E5D","#A70000","#C18333","#C34300","#8E769A"
                           ,"#2E015E","#C3C35D","#751D00","#4C7489","#001E5A","#588530","#004600","#A1403F","#890000","#A36515","#A52500"
                           ,"#70587C","#100040","#A5A53F","#570000"),"99")
    set.seed(123);pseudo.col <- sample(pseudo.col,nrow(otu.cor.pc))
    if(logdata==TRUE){
      plot(x=seq(min(otu.cor.pc),max(otu.cor.pc),length.out = 10),
           y=seq(min(PC),max(PC),length.out = 10),
           las=1,type="n",bty="L",
           xlab="Log OTU counts",ylab=paste("Log", "PCoA",i))
    } else if(logdata==FALSE){
      plot(x=seq(min(otu.cor.pc),max(otu.cor.pc),length.out = 10),
           y=seq(min(PC),max(PC),length.out = 10),
           las=1,type="n",bty="L",
           xlab="OTU counts",ylab=paste("PCoA",i))
    } 
    for(i in 1:nrow(otu.cor.pc)){
      abline(lm(PC~t(otu.cor.pc[i,])),col="#CACACA99",lty=3)
      points(x=t(otu.cor.pc[i,]),y=PC,pch=21,col="#CACACA99",bg=pseudo.col[i])
    }
  }
  return(list(otu.cor.pc=otu.cor.pc,sig.corr=sig.corr,sig.p=sig.p))
}


# Barplot for Phylum to family levels
TaxaLvABar <- function(otu,taxafromQ2tsv,level=2,groups=NULL,sample.order=NULL,border = "#E0E0E044",colsets,pctcutoff=1,legend=FALSE,inset=-0.05){
  taxafromQ2tsv <- taxafromQ2tsv[match(rownames(otu),taxafromQ2tsv[,1]),]
  taxa <- strsplit(as.character(taxafromQ2tsv[,2]), "; ")
  taxx <- data.frame(Level_1=NA,Level_2=NA,Level_3=NA,Level_4=NA,Level_5=NA,Level_6=NA,Level_7=NA)
  for(i in 1:length(taxa)) taxx[i,1:length(taxa[[i]])] <- taxa[[i]][1:length(taxa[[i]])]
  taxa <- taxx
  taxa.comb <- taxa
  taxa.comb$Level_3 <- paste(taxa$Level_2,taxa$Level_3,sep="; ")
  taxa.comb$Level_4 <- paste(taxa$Level_2,taxa$Level_3,taxa$Level_4,sep="; ")
  taxa.comb$Level_5 <- paste(taxa$Level_2,taxa$Level_3,taxa$Level_4,taxa$Level_5,sep="; ")
  taxa.comb$Level_6 <- paste(taxa$Level_2,taxa$Level_3,taxa$Level_4,taxa$Level_5,taxa$Level_6,sep="; ")
  taxa.comb$Level_7 <- paste(taxa$Level_2,taxa$Level_3,taxa$Level_4,taxa$Level_5,taxa$Level_6,taxa$Level_7,sep="; ")
  tx.level <- function(data=otu,level,percentile=T){
    level.n <- taxa.comb[,level]
    sort.sum <- aggregate(data,list(level.n),sum)
    rownames(sort.sum) <- sort.sum[,1]
    sort.sum <- sort.sum[,-1]
    if(percentile){
      sort.sum <- prop.table(as.matrix(sort.sum),2)*100
    } else sort.sum <- sort.sum
    return(sort.sum)
  }
  
  phyla <- tx.level(data = otu,2)
  levelsum <- tx.level(data = otu,level)
  
  #par(mar=c(5.1,4.1,4.1,15.1))
  if(!is.null(groups)){
    levelsum <- as.matrix(aggregate2df(t(levelsum),groups,mean))
    colx <- rowMeans(phyla)
    colx[which(rowMeans(phyla) > pctcutoff)] <- colsets[seq_len(sum(rowMeans(phyla) > pctcutoff))]
    colx[which(rowMeans(phyla) < pctcutoff)] <- "#F0F0F0"
    {
      lowlv <- data.frame(taxa=sort(c(unique(taxa.comb$Level_3),unique(taxa.comb$Level_4),unique(taxa.comb$Level_5))))
      
      taxa.ls <- c()
      for(i in 1:nrow(phyla)) taxa.ls[i] <- length(grep(paste0(rownames(phyla)[i],"; "),lowlv$taxa,fixed = T))
      hexs <- list()
      for(j in 1:length(taxa.ls)){
        trans <- col2rgb(colx[j])
        rgbs <- rbind(
          seq(trans[1],0,length.out = taxa.ls[j]+5),
          seq(trans[2],0,length.out = taxa.ls[j]+5),
          seq(trans[3],0,length.out = taxa.ls[j]+5)
        )
        hext <- c()
        for(i in 1:(taxa.ls[j]+5)) {
          hext[i] <-rgb(red = rgbs[1,i],green = rgbs[2,i],blue = rgbs[3,i],maxColorValue = 255)
        }
        hexs[[j]] <- hext[2:(length(hext)-4)]
      }
      
      lowlv$col <- c(unlist(hexs),rep("#FFFFFF",length(grep("^NA;",lowlv$taxa))))
      
    }
    if(level<=2){
      t1 <- barplot(levelsum[nrow(levelsum):1,],col=colx,las=2,ylab="Relative abundance (%)",cex.names = 0.7,border = border,lwd=0.15,yaxt="n");axis(2,las=2)
    } else if(level>2) t1 <- barplot(levelsum[nrow(levelsum):1,],col=lowlv[match(rownames(levelsum), lowlv$taxa),2],las=2,ylab="Relative abundance (%)",cex.names = 0.7,border = border,lwd=0.15,yaxt="n");axis(2,las=2)
  } else if(is.null(groups)){
    colx <- rowMeans(phyla)
    colx[which(rowMeans(phyla) > pctcutoff)] <- colsets[seq_len(sum(rowMeans(phyla) > pctcutoff))]
    colx[which(rowMeans(phyla) < pctcutoff)] <- "#F0F0F0"
    {
      lowlv <- data.frame(taxa=sort(c(unique(taxa.comb$Level_3),unique(taxa.comb$Level_4),unique(taxa.comb$Level_5))))
      
      taxa.ls <- c()
      for(i in 1:nrow(phyla)) taxa.ls[i] <- length(grep(paste0(rownames(phyla)[i],"; "),lowlv$taxa,fixed = T))
      hexs <- list()
      for(j in 1:length(taxa.ls)){
        trans <- col2rgb(colx[j])
        rgbs <- rbind(
          seq(trans[1],0,length.out = taxa.ls[j]+5),
          seq(trans[2],0,length.out = taxa.ls[j]+5),
          seq(trans[3],0,length.out = taxa.ls[j]+5)
        )
        hext <- c()
        for(i in 1:(taxa.ls[j]+5)) {
          hext[i] <-rgb(red = rgbs[1,i],green = rgbs[2,i],blue = rgbs[3,i],maxColorValue = 255)
        }
        hexs[[j]] <- hext[2:(length(hext)-4)]
      }
      
      lowlv$col <- c(unlist(hexs),rep("#FFFFFF",length(grep("^NA;",lowlv$taxa))))
      
    }
    if(!is.null(sample.order)){
      if(level<=2){
        t1 <- barplot(levelsum[nrow(levelsum):1,sample.order],col=colx,las=2,ylab="Relative abundance (%)",cex.names = 0.7,border = border,lwd=0.15,yaxt="n");axis(2,las=2)
      } else if(level>2) t1 <- barplot(levelsum[nrow(levelsum):1,sample.order],col=lowlv[match(rownames(levelsum), lowlv$taxa),2],las=2,ylab="Relative abundance (%)",cex.names = 0.7,border = border,lwd=0.15,yaxt="n");axis(2,las=2)
    } else if(is.null(sample.order) & level<=2){
      t1 <- barplot(levelsum[nrow(levelsum):1,],col=colx,las=2,ylab="Relative abundance (%)",cex.names = 0.7,border = border,lwd=0.15,yaxt="n");axis(2,las=2)
    } else if(is.null(sample.order) & level>2) t1 <- barplot(levelsum[nrow(levelsum):1,],col=lowlv[match(rownames(levelsum), lowlv$taxa),2],las=2,ylab="Relative abundance (%)",cex.names = 0.7,border = border,lwd=0.15,yaxt="n");axis(2,las=2)
  }
  
  legend.names <- c(rownames(levelsum)[which(rowMeans(levelsum) > pctcutoff)],"others")
  if(level<=2){
    legend.cols <- c(colx[which(rowMeans(levelsum) > pctcutoff)],"#F0F0F0")
  } else if(level>2) legend.cols <- c(lowlv[match(rownames(levelsum), lowlv$taxa),2][which(rowMeans(levelsum) > pctcutoff)],"#F0F0F0")
  if(legend==TRUE) legend("topright",legend = legend.names,fill = rev(legend.cols),bty="n",cex=0.8,inset = c(inset,0),xpd=T)
  
}



# oblique gradient
oblique.gd <- function(ref.point.xy,ordi.xy){
  slope <- ref.point.xy[2]/ref.point.xy[1]
  #abline(a=0,b=slope)
  
  ys <- ordi.xy[,1]*slope
  #points(ordi.xy[,1],ys,pch=19,col=2,cex=0.5)
  #segments(x0 = ordi.xy[,1],y0 = ordi.xy[,2], x1 = ordi.xy[,1],y1 = ys)
  
  if(slope!=0){
    minus.a <- -((-1/slope)*ordi.xy[,1]-ordi.xy[,2])
    cross.xy <- matrix(NA,nrow=2,ncol = length(minus.a))
    for(i in 1:length(minus.a)){
      #abline(a=minus.a[i],b=-1/slope,lty=3,col="gray50")
      cross.xy[,i] <- solve(matrix(c(slope,1,-1/slope,1),nrow = 2,byrow = T),matrix(c(0,-minus.a[i]),nrow=2))
    }
    ys2 <- cross.xy[1,]*slope
    #points(cross.xy[1,],ys2,pch=19,col=3,cex=0.5)
    ys3 <- cross.xy[1,]
  } else if(slope==0) ys3 <- ordi.xy[,1]
  
  return(obgradient=ys3)
}

Xaxis.transform <- function(ordi.xy,slope){
  minus.a <- -((-1/slope)*ordi.xy[,1]-ordi.xy[,2])
  cross.xy <- matrix(NA,nrow=2,ncol = length(minus.a))
  for(i in 1:length(minus.a)){
    #abline(a=0,b=-1/slope,lty=3,col="gray50")
    cross.xy[,i] <- solve(matrix(c(slope,1,-1/slope,1),nrow = 2,byrow = T),matrix(c(0,-minus.a[i]),nrow=2))
  }
  ys2 <- cross.xy[1,]*slope
  #points(cross.xy[1,],ys2,pch=19,col=3,cex=0.5)
  ys3 <- cross.xy[1,]
  return(data.frame(x=ys3,y=ys2))
}


###### LEfSe ######
# combine taxa levels
tx.level <- function(data,taxa,level,percentile=T){
  taxa.comb <- taxa
  taxa.comb[,3] <- paste(taxa[,2],taxa[,3],sep="; ")
  taxa.comb[,4] <- paste(taxa[,2],taxa[,3],taxa[,4],sep="; ")
  taxa.comb[,5] <- paste(taxa[,2],taxa[,3],taxa[,4],taxa[,5],sep="; ")
  taxa.comb[,6] <- paste(taxa[,2],taxa[,3],taxa[,4],taxa[,5],taxa[,6],sep="; ")
  taxa.comb[,7] <- paste(taxa[,2],taxa[,3],taxa[,4],taxa[,5],taxa[,6],taxa[,7],sep="; ")
  level.n <- taxa.comb[,level]
  sort.sum <- aggregate(data,list(level.n),sum)
  rownames(sort.sum) <- sort.sum[,1]
  sort.sum <- sort.sum[,-1]
  if(percentile){
    sort.sum <- prop.table(as.matrix(sort.sum),2)*100
  } else sort.sum <- sort.sum
  return(sort.sum)
}

#LDA v1 only for 2-group comparison
lefse <- function(data,taxa,Groups,Lefse.factor=1,FDR=FALSE,alpha=0.05,LDAscore=2,plot=TRUE,col){
  options(warn = -1)
  phyla <- tx.level(data = data,taxa = taxa,level = 2)
  classes <- tx.level(data = data,taxa = taxa,level = 3)
  orders <- tx.level(data = data,taxa = taxa,level = 4)
  families <- tx.level(data = data,taxa = taxa,level = 5)
  genera <- tx.level(data = data,taxa = taxa,level = 6)
  species <- tx.level(data = data,taxa = taxa,level = 7)
  
  lefse.tab <- rbind(phyla,classes,orders,families,genera,species)
  lefse.tab <- (lefse.tab/100)*Lefse.factor
  
  kw.test <- apply(lefse.tab,MARGIN = 1,function(x) kruskal.test(x~Groups))
  p.val <- sapply(kw.test, "[[",3)
  p.adj <- p.adjust(p.val,"fdr")
  
  if(FDR==TRUE){
    lefse.tab.kw <- lefse.tab[which(p.adj < alpha),]
  } else if(FDR==FALSE) lefse.tab.kw <- lefse.tab[which(p.val < alpha),]
  wc.test <- apply(lefse.tab.kw,MARGIN = 1,function(x) wilcox.test(x~Groups))
  p.val.wc <- sapply(wc.test, "[[",3)
  p.adj.wc <- p.adjust(p.val.wc,"fdr")
  if(FDR==TRUE){
    lefse.tab.wc <- lefse.tab.kw[which(p.adj.wc < alpha),]
  } else if(FDR==FALSE) lefse.tab.wc <- lefse.tab.kw[which(p.val.wc < alpha),]
  
  LDA <- MASS::lda(Groups~t(lefse.tab.wc))
  lda.score <- data.frame(taxa=rownames(lefse.tab.wc),score=LDA$scaling,row.names = NULL)
  ls.1 <- lda.score[which(lda.score$LD1<0),]
  ls.1 <- ls.1[order(ls.1$LD1,decreasing = T),]
  ls.2 <- lda.score[which(lda.score$LD1>0),]
  ls.2 <- ls.2[order(ls.2$LD1),]
  lda.score <- rbind(ls.1,ls.2)
  lda.score$LD1 <- lda.score$LD1*Lefse.factor
  rm(ls.1,ls.2)
  lda.score$EffectSize <- log10(abs(lda.score$LD1))
  lda.score$EffectSize[lda.score$LD1<0] <- -1*lda.score$EffectSize[lda.score$LD1<0]
  lda.score <- lda.score[abs(lda.score$EffectSize)>LDAscore,]
  
  if(plot==TRUE){
    par(mar=c(5.1,32.1,4.1,2.1))
    by <- barplot(lda.score$EffectSize,horiz = T,col=col[(lda.score$EffectSize>0)+1],
                  border = col[(lda.score$EffectSize>0)+1],
                  xlab=expression(paste(Log[10], " LDA score")),xlim=c(-5,5))
    axis(2,by,labels = lda.score$taxa,las=1,cex.axis=0.7)
    legend("bottomright",legend = levels(Groups),fill = col,border = col,xpd=T,bty="n",cex = 0.9)
  }
  options(warn = 0)
  return(list(LEfSe=lda.score,LDA.data.pXfactor=lefse.tab.wc))
}

# correlation heatmap
corheatmap <- function(otutable,taxonomy,corr,percent.cf=0.975,percent.cf.pos=TRUE,parm,parm.cf,parm.gp=c("Group1","Group2"),log=TRUE,
                       GroupName="Phenotype",ParmName="Levels",
                       groupcol=c("#d1221d","#003366"),phenocol=c("#ffda53","#cc092f"),corrcol=c("#c0d4eb","#3b96ff"),
                       heatcol=c("#0f054f","#345c66","#fef4ad","#e09e3e","#9d2a2d","#56090f"),...){
  if(percent.cf.pos==TRUE){
    cor.chmp <- corr[which(corr >= quantile(abs(corr),percent.cf))] # correlation coefficients subset with % cutoff
    data.chmp <- otutable[which(corr >= quantile(abs(corr),percent.cf)),] # otu table subset according to correlation coefficients subset
  } else if(percent.cf.pos==FALSE){
    cor.chmp <- corr[which(abs(corr) >= quantile(abs(corr),percent.cf))]
    data.chmp <- otutable[which(abs(corr) >= quantile(abs(corr),percent.cf)),]
  }
  
  
  taxonomy <- taxonomy[order(cor.chmp,decreasing = T),] # sort taxonomy (8 columns: ID + 7 lvs taxonomy)
  data.chmp <- data.chmp[order(cor.chmp,decreasing = T),] # sort otu table row
  cor.chmp <- cor.chmp[order(cor.chmp,decreasing = T)] # sort correlation coefficients
  
  if(log==TRUE){
    df <- log1p(data.chmp[,order(parm)]) # sort otu table column by metadata, w/ log1p transform
  } else if(log==FALSE){
    df <- data.chmp[,order(parm)]
  }
  
  
  annotation_cols <- data.frame(ifelse(parm>parm.cf,parm.gp[1],parm.gp[2]),
                                parm,
                                row.names = colnames(data.chmp))
  colnames(annotation_cols) <- c(GroupName,ParmName)
  annotation_cols <- annotation_cols[order(parm),]
  annotation_rows <- data.frame("r"=cor.chmp,row.names = rownames(data.chmp))
  annotation_colors <- list(c(groupcol[1],groupcol[2]),
                            colorRampPalette(phenocol)(length(unique(parm))),
                            colorRampPalette(corrcol)(length(unique(round(cor.chmp,2)))))
  names(annotation_colors) <- c(GroupName,ParmName,"r")
  names(annotation_colors[[1]]) <- c(parm.gp[1],parm.gp[2])
  
  pheatmap::pheatmap(df,cluster_cols = F,cluster_rows = F,
                     scale = "row",show_rownames = T,show_colnames = T,
                     fontsize_row = 7,fontsize_col = 7,
                     annotation_col = annotation_cols,
                     annotation_row = annotation_rows,
                     annotation_colors = annotation_colors,
                     color = colorRampPalette(heatcol)(100),
                     labels_row = taxonomy[,8],...)
}

##### taxonomy refine for heatmap#####
taxa4hmp <- function(taxa_w_8col){
  taxa_w_8col <- as.data.frame(apply(taxa_w_8col, 2, function(x) gsub("D_.__","",x)))
  for(i in 2:8){
    taxa_w_8col[,i] <- ifelse(taxa_w_8col[,i]=="Unassigned",
                              as.character(taxa_w_8col[,i-1]),
                              as.character(taxa_w_8col[,i]))
    taxa_w_8col[,i] <- ifelse(taxa_w_8col[,i]=="uncultured" | taxa_w_8col[,i]=="uncultured bacterium",
                              as.character(taxa_w_8col[,i-1]),
                              as.character(taxa_w_8col[,i]))
  }
  taxa_w_8col[,8] <- ifelse(taxa_w_8col[,8]==taxa_w_8col[,7],paste(taxa_w_8col[,8],"sp."),taxa_w_8col[,8])
  taxa_w_8col[grep("uncultured",taxa_w_8col[,8]),8] <- 
    paste(taxa_w_8col[grep("uncultured",taxa_w_8col[,8]),7],
          taxa_w_8col[grep("uncultured",taxa_w_8col[,8]),8])
  taxa_w_8col[grep("metagenome",taxa_w_8col[,8]),8] <- 
    paste(taxa_w_8col[grep("metagenome",taxa_w_8col[,8]),7],
          taxa_w_8col[grep("metagenome",taxa_w_8col[,8]),8])
  return(taxa_w_8col)
}

##### Internal functions #####
#old version vegan::adonis
adonisx <- function (formula, data = NULL, permutations = 999, method = "bray", 
                     strata = NULL, contr.unordered = "contr.sum", contr.ordered = "contr.poly", 
                     parallel = getOption("mc.cores"), ...) 
{
  EPS <- sqrt(.Machine$double.eps)
  TOL <- 1e-07
  Terms <- terms(formula, data = data)
  lhs <- formula[[2]]
  lhs <- eval(lhs, data, parent.frame())
  formula[[2]] <- NULL
  rhs.frame <- model.frame(formula, data, drop.unused.levels = TRUE)
  op.c <- options()$contrasts
  options(contrasts = c(contr.unordered, contr.ordered))
  rhs <- model.matrix(formula, rhs.frame)
  options(contrasts = op.c)
  grps <- attr(rhs, "assign")
  qrhs <- qr(rhs)
  rhs <- rhs[, qrhs$pivot, drop = FALSE]
  rhs <- rhs[, 1:qrhs$rank, drop = FALSE]
  grps <- grps[qrhs$pivot][1:qrhs$rank]
  u.grps <- unique(grps)
  nterms <- length(u.grps) - 1
  if (nterms < 1) 
    stop("right-hand-side of formula has no usable terms")
  H.s <- lapply(2:length(u.grps), function(j) {
    Xj <- rhs[, grps %in% u.grps[1:j]]
    qrX <- qr(Xj, tol = TOL)
    Q <- qr.Q(qrX)
    tcrossprod(Q[, 1:qrX$rank])
  })
  if (inherits(lhs, "dist")) {
    if (any(lhs < -TOL)) 
      stop("dissimilarities must be non-negative")
    dmat <- as.matrix(lhs^2)
  }
  else if ((is.matrix(lhs) || is.data.frame(lhs)) && isSymmetric(unname(as.matrix(lhs)))) {
    dmat <- as.matrix(lhs^2)
    lhs <- as.dist(lhs)
  }
  else {
    dist.lhs <- as.matrix(vegdist(lhs, method = method, 
                                  ...))
    dmat <- dist.lhs^2
  }
  n <- nrow(dmat)
  G <- -sweep(dmat, 1, rowMeans(dmat))/2
  SS.Exp.comb <- sapply(H.s, function(hat) sum(G * t(hat)))
  SS.Exp.each <- c(SS.Exp.comb - c(0, SS.Exp.comb[-nterms]))
  H.snterm <- H.s[[nterms]]
  tIH.snterm <- t(diag(n) - H.snterm)
  if (length(H.s) > 1) 
    for (i in length(H.s):2) H.s[[i]] <- H.s[[i]] - H.s[[i - 
                                                           1]]
  SS.Res <- sum(G * tIH.snterm)
  df.Exp <- sapply(u.grps[-1], function(i) sum(grps == i))
  df.Res <- n - qrhs$rank
  if (inherits(lhs, "dist")) {
    beta.sites <- qr.coef(qrhs, as.matrix(lhs))
    beta.spp <- NULL
  }
  else {
    beta.sites <- qr.coef(qrhs, dist.lhs)
    beta.spp <- qr.coef(qrhs, as.matrix(lhs))
  }
  colnames(beta.spp) <- NULL
  colnames(beta.sites) <- NULL
  F.Mod <- (SS.Exp.each/df.Exp)/(SS.Res/df.Res)
  f.test <- function(tH, G, df.Exp, df.Res, tIH.snterm) {
    (sum(G * tH)/df.Exp)/(sum(G * tIH.snterm)/df.Res)
  }
  SS.perms <- function(H, G, I) {
    c(SS.Exp.p = sum(G * t(H)), S.Res.p = sum(G * t(I - 
                                                      H)))
  }
  p <- vegan:::getPermuteMatrix(permutations, n, strata = strata)
  permutations <- nrow(p)
  if (permutations) {
    tH.s <- lapply(H.s, t)
    if (is.null(parallel)) 
      parallel <- 1
    hasClus <- inherits(parallel, "cluster")
    isParal <- hasClus || parallel > 1
    isMulticore <- .Platform$OS.type == "unix" && !hasClus
    if (isParal && !isMulticore && !hasClus) {
      parallel <- makeCluster(parallel)
    }
    if (isParal) {
      if (isMulticore) {
        f.perms <- sapply(1:nterms, function(i) unlist(mclapply(1:permutations, 
                                                                function(j) f.test(tH.s[[i]], G[p[j, ], p[j, 
                                                                                                          ]], df.Exp[i], df.Res, tIH.snterm), mc.cores = parallel)))
      }
      else {
        f.perms <- sapply(1:nterms, function(i) parSapply(parallel, 
                                                          1:permutations, function(j) f.test(tH.s[[i]], 
                                                                                             G[p[j, ], p[j, ]], df.Exp[i], df.Res, tIH.snterm)))
      }
    }
    else {
      f.perms <- sapply(1:nterms, function(i) sapply(1:permutations, 
                                                     function(j) f.test(tH.s[[i]], G[p[j, ], p[j, 
                                                                                               ]], df.Exp[i], df.Res, tIH.snterm)))
    }
    if (isParal && !isMulticore && !hasClus) 
      stopCluster(parallel)
    P <- (rowSums(t(f.perms) >= F.Mod - EPS) + 1)/(permutations + 
                                                     1)
  }
  else {
    f.perms <- P <- rep(NA, nterms)
  }
  SumsOfSqs = c(SS.Exp.each, SS.Res, sum(SS.Exp.each) + SS.Res)
  tab <- data.frame(Df = c(df.Exp, df.Res, n - 1), SumsOfSqs = SumsOfSqs, 
                    MeanSqs = c(SS.Exp.each/df.Exp, SS.Res/df.Res, NA), 
                    F.Model = c(F.Mod, NA, NA), R2 = SumsOfSqs/SumsOfSqs[length(SumsOfSqs)], 
                    P = c(P, NA, NA))
  rownames(tab) <- c(attr(attr(rhs.frame, "terms"), "term.labels")[u.grps], 
                     "Residuals", "Total")
  colnames(tab)[ncol(tab)] <- "Pr(>F)"
  attr(tab, "heading") <- c(vegan:::howHead(attr(p, "control")), "Terms added sequentially (first to last)\n")
  class(tab) <- c("anova", class(tab))
  out <- list(aov.tab = tab, call = match.call(), coefficients = beta.spp, 
              coef.sites = beta.sites, f.perms = f.perms, model.matrix = rhs, 
              terms = Terms)
  class(out) <- "adonis"
  out
}

#MDI: microbial dysbiosis index
MDI <- function(otu,groups){
  message("The healthy control group must be 1 and the disease group must be 2 in grouping factor")
  prop.otu <- as.data.frame(prop.table(t(otu),2))
  group.mean <- aggregate(prop.otu,list(groups),mean)
  rownames(group.mean) <- group.mean$Group.1
  group.mean <- group.mean[,-1]
  
  idx.ics <- which(group.mean[1,] < group.mean[2,])
  idx.dcs <- which(group.mean[1,] > group.mean[2,])
  
  MDI <- c()
  for(i in 1:nrow(prop.otu)){
    MDI[i] <- log(sum(prop.otu[i,idx.ics])/sum(prop.otu[i,idx.dcs]))
  }
  return(MDI)
}

#single microbe logistic regression, output odds ratio and 95% CI
smilr <- function(otu,groups,control=NULL,case=NULL){
  logit.p <- list()
  for(i in 1:nrow(otu)){
    if(is.null(control) & is.null(case)){
      log_data <- data.frame(species=t(log1p(otu))[,i],
                             gx=as.numeric(factor(groups))-1)
    } else if(!is.null(control) & !is.null(case)){
      log_data <- data.frame(species=t(log1p(otu))[,i],
                             gx=as.numeric(factor(groups,levels = c(control,case)))-1)
    }
    logit <- glm(gx~species,data=log_data,family = "binomial")
    coef <- summary(logit)$coefficients[2,c(1,4)]
    OR <- exp(cbind(OR=coef(logit),confint(logit)))[2,]
    logit.p[[i]] <- c(coef,OR)
  }
  logit.p <- as.data.frame(do.call(rbind,logit.p))
  rownames(logit.p) <- rownames(DESeq_output)
  return(logit.p)
}
