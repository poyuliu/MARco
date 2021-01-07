#list.of.packages <- c("gage")
#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) {
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#  BiocManager::install("gage")
#}
#rm(list = c("list.of.packages","new.packages"))
#
#rm(list=ls())
#picrust2 <- read.delim("KEGG/pred_metagenome_unstrat_descrip.tsv",row.names = 1)
#picrust2 <- picrust2[,-1]
#colnames(picrust2) <- gsub("Ke","Exp.",colnames(picrust2))
#saveRDS(picrust2,file = "example_PICRUSt2.RDS")
#gs.table <- KO2path(predtable = picrust2,"module")
#gs.table <- KO2path(predtable = picrust2,"pathway")
#lv2.table <- path.lvs(gs.table,2)

# Example file: K. eiffingeri tadpole gut microbiome, PICRUSt2 metagenome prediction using KEGG database
#picrust2 <- readRDS("example_PICRUSt2.RDS")
#head(picrust2)

# categorized KEGG orthology (KO) IDs [abundance table] to pathways or modules/ Genesets: all pathway, metabolic pathway, modules/ options: convert to proportional table 
KO2path <- function(predtable,GeneSet="pathway",prop=FALSE){
  if(GeneSet=="pathway"){
    genesets <- readRDS(url("http://www.lifescipy.net/RcodeDB/KEGG/AA_KEGG_GENESETS.RDS"))
  } else if(GeneSet=="metabolic"){
    genesets <- readRDS(url("http://www.lifescipy.net/RcodeDB/KEGG/AA_KEGG_METABOLIC.RDS")) 
  } else if(GeneSet=="module"){
    genesets <- readRDS(url("http://www.lifescipy.net/RcodeDB/KEGG/AA_KEGG_MODULES.RDS"))
    genesets <- sapply(genesets,"[[",1)
    genesets <- sapply(genesets,as.character)
  }

  gs.table <- data.frame(matrix(NA,nrow = length(genesets),ncol = ncol(predtable)))
  colnames(gs.table) <- colnames(predtable)
  rownames(gs.table) <- names(genesets)
  for(i in 1:length(genesets)){
    m.idx <- match(genesets[[i]], rownames(predtable))
    m.idx <- m.idx[!is.na(m.idx)]
    if(length(m.idx)>1){
      gs.table[i,] <- colSums(predtable[m.idx,])
    } else if(length(m.idx)==1){
      gs.table[i,] <- predtable[m.idx,]
    } else gs.table[i,] <- 0
  }
  
  gs.table <- gs.table[rowSums(gs.table) != 0,]
  gs.table.p <- as.data.frame(prop.table(as.matrix(gs.table),margin = 2))
  
  if(prop==FALSE){
    return(gs.table)
  } else if(prop==TRUE){
    return(gs.table.p)
  }
  
}

# sort pathway/module categoized table with functional levels
path.lvs <- function(gs.table,aggregate.lv=NULL){
  aggregate2dff <- function(data,groups,FUN){
    agg.data <- aggregate(data,list(groups),FUN)
    rownames(agg.data) <- agg.data[,1]
    agg.data <- as.data.frame(agg.data[,-1])
    return(agg.data)
  }
  if(substr(rownames(gs.table)[1],start = 1,stop = 1)=="M"){
    lvs <- readRDS(url("http://www.lifescipy.net/RcodeDB/KEGG/AA_KEGG_M_LEVELS.RDS"))
    lvs <- lvs[which(lvs$Level_4 %in% substr(rownames(gs.table),1,6)),]
    gs.table <- gs.table[match(lvs$Level_4,substr(rownames(gs.table),1,6)),]
    gs.table <- cbind(lvs,gs.table)
    rownames(gs.table) <- lvs$Level_4
    gs.table <- gs.table[,-4]
    colnames(gs.table)[4] <- "Level_4"
    if(is.null(aggregate.lv)){
      gs.table <- gs.table
    } else if(!is.null(aggregate.lv) & aggregate.lv <4){
      gs.table <- gs.table.agg <- aggregate2df(gs.table[,-c(1:4)],gs.table[,aggregate.lv],sum)
    }
  } else if(substr(rownames(gs.table)[1],start = 1,stop = 1)=="k"){
    lvs <- readRDS(url("http://www.lifescipy.net/RcodeDB/KEGG/AA_KEGG_LEVELS.RDS"))
    lvs <- lvs[which(lvs$Level_3 %in% substr(rownames(gs.table),1,7)),]
    gs.table <- gs.table[match(lvs$Level_3,substr(rownames(gs.table),1,7)),]
    gs.table <- cbind(lvs,gs.table)
    rownames(gs.table) <- lvs$Level_3
    gs.table <- gs.table[,-3]
    colnames(gs.table)[3] <- "Level_3"
    if(is.null(aggregate.lv)){
      gs.table <- gs.table
    } else  if(!is.null(aggregate.lv) & aggregate.lv < 3){
      gs.table <- gs.table.agg <- aggregate2dff(gs.table[,-c(1:3)],gs.table[,aggregate.lv],sum)
    }
  }
    return(gs.table)
}

# FSEA: Functional Sets Enrichment Analysis
# within test: compare to "expected value" within group/individual
# between test: compare to reference group/individual (NOT COMPLETE)
FSEA <- function(Fdata,Groups,FunctionalSets="pathway",SortLv=2,Test="within",alpha=0.05,plot=TRUE,barcol=NULL,legend=TRUE){
  aggregate2dfx <- function(data,groups,FUN){
    agg.data <- aggregate(data,list(groups),FUN)
    rownames(agg.data) <- agg.data[,1]
    agg.data <- as.data.frame(t(agg.data[,-1]))
    return(agg.data)
  }
  
  if(FunctionalSets=="pathway"){
    kegg.gs <- readRDS(url("http://www.lifescipy.net/RcodeDB/KEGG/AA_KEGG_GENESETS.RDS"))
    Levels.KEGG <- readRDS(url("http://www.lifescipy.net/RcodeDB/KEGG/AA_KEGG_LEVELS.RDS"))
  } else if(FunctionalSets=="metabolic"){
    kegg.gs <- readRDS(url("http://www.lifescipy.net/RcodeDB/KEGG/AA_KEGG_METABOLIC.RDS"))
    Levels.KEGG <- readRDS(url("http://www.lifescipy.net/RcodeDB/KEGG/AA_KEGG_LEVELS.RDS"))
  } else if(FunctionalSets=="module"){
    kegg.gs <- readRDS(url("http://www.lifescipy.net/RcodeDB/KEGG/AA_KEGG_MODULES.RDS"))
    kegg.gs <- sapply(kegg.gs,"[[",1)
    kegg.gs <- sapply(kegg.gs,as.character)
    Levels.KEGG <- readRDS(url("http://www.lifescipy.net/RcodeDB/KEGG/AA_KEGG_M_LEVELS.RDS"))
  }
  if(Test=="within"){
    qval <- data.frame(matrix(nrow=length(kegg.gs),ncol=ncol(Fdata)),row.names = names(kegg.gs))
    colnames(qval) <- colnames(Fdata)
    for(i in 1:ncol(Fdata)){
      PropData <- prop.table(as.matrix(Fdata)[,i])
      gagetest <- gage::gage(PropData,gset=kegg.gs,rank.test = T,saaTest = gage::gs.zTest)
      gagetest <- as.data.frame(gagetest[[1]])
      gagetest <- gagetest[order(rownames(gagetest)),]
      qval[,i] <- gagetest$q.val
    }
    EScore <- -log10(qval)
    EScore <- EScore[!is.na(rowSums(EScore)),]
    EScore.mean <- aggregate2dfx(t(EScore),Groups,mean)
    EScore.sem <- aggregate2dfx(t(EScore),Groups,sd)/sqrt(as.vector(table(Groups)))
    EScore.idx <- EScore.mean
    EScore.idx[EScore.mean < (-log10(alpha))] <- 0
    EScore.idx[EScore.mean > (-log10(alpha))] <- 1
    EScore.mean <- EScore.mean[rowSums(EScore.idx)!=0,]
    EScore.sem <- EScore.sem[rowSums(EScore.idx)!=0,]
    EScore.mean.order <- EScore.mean[order(rowSums(EScore.mean)),]
    EScore.sem.order <- EScore.sem[order(rowSums(EScore.mean)),]
    if(FunctionalSets=="module"){
      gagetest.Lv <- as.character(Levels.KEGG[,SortLv])[match(substr(rownames(EScore.mean.order),1,6),Levels.KEGG$Level_4)]
    } else 
      gagetest.Lv <- as.character(Levels.KEGG[,SortLv])[match(substr(rownames(EScore.mean.order),1,6),Levels.KEGG$Level_3)]
    
    EScore.mean.order <- EScore.mean.order[order(gagetest.Lv,decreasing = T),]
    EScore.sem.order <- EScore.sem.order[order(gagetest.Lv,decreasing = T),]
    EScore.order <- list(EScore.mean=EScore.mean.order,EScore.sem=EScore.sem.order)
    
    if(plot==TRUE){
      par(mar=c(5.1,20.1,1.1,2.1),cex.lab=0.7,cex.axis=0.7)
      bars <- barplot(t(EScore.mean.order),beside = T,horiz = T,
                      las=1,xlab="Enrichment score",border = 0,cex.names = 0.7,
                      col = barcol)
      errobar <- as.data.frame(cbind(xaxis=as.vector(bars),
                                     upper=as.vector(t(EScore.mean.order))+as.vector(t(EScore.sem.order)),
                                     lower=as.vector(t(EScore.mean.order))-as.vector(t(EScore.sem.order))))
      segments(x0 = errobar$upper,y0 = errobar$xaxis,x1 = errobar$lower,y1 = errobar$xaxis)
      abline(v=-log10(alpha),lty=3)
      if(legend==TRUE){
        if(!is.null(barcol)){
          legend("bottomright",legend = levels(as.factor(Groups)),bty="n",fill=barcol,cex = 0.7)
        } else legend("bottomright",legend = levels(as.factor(Groups)),bty="n",
                      fill=gray.colors(length(unique(as.factor(Groups)))),cex = 0.7)
      }
    }
    
    return(EScore.order)
  } else if(Test=="between"){
    PropData <- prop.table(as.matrix(Fdata),margin = 2)
    gageX <- list()
    cbm <- cbind(combn(levels(Groups),2),combn(levels(Groups),2)[2:1,])
    for(k in 1:ncol(cbm)){
      gageX[[k]] <- gage::gage(PropData,gset=kegg.gs,ref=which(Groups==cbm[2,k]),samp=which(Groups==cbm[1,k]),rank.test = F,
                         saaTest = gage::gs.tTest,compare = "unpaired")
      gageX[[k]] <- as.data.frame(gageX[[k]][[1]])
      gageX[[k]] <- gageX[[k]][order(rownames(gageX[[k]])),]
    }
    
    qvals <- data.frame(do.call(cbind,lapply(gageX,"[","q.val")))
    colnames(qvals) <- paste0("r:",t(cbm)[,1],"-","s:",t(cbm)[,2])
    qvals[is.na(qvals)] <- 1
    qvals.sig <- qvals[which(rowSums(qvals<alpha)>=1),]
    Escore.sig <- -log10(qvals.sig)
    
    ESx <- list()
    if(length(levels(groups)==2)){
      m <- 1
      ESx[[m]] <- Escore.sig
      
      if(plot==TRUE){
        par(mar=c(4.1,18.1,1.1,2.1),cex.lab=0.7,cex.axis=0.7)
          barplot(t(ESx[[m]]),beside = T,horiz = T,las=1,cex.name=0.7,xlim=c(0,12),col = barcol,border = F)
          abline(v=-log10(alpha),lty=3)
          legend("bottomright",legend = colnames(ESx[[m]]),
                 bty="n",fill=barcol,cex = 0.7)
      }
      
    } else if(length(levels(groups)>2)){
      for (m in 1:length(levels(groups))) {
        ESx[[m]] <- Escore.sig[, grep(paste0("r:", levels(groups)[m], 
                                             "-"), names(Escore.sig))]
        ESx[[m]] <- ESx[[m]][order(rowSums(ESx[[m]])), ]
      }
      if(plot==TRUE){
        par(mar=c(4.1,18.1,1.1,2.1),cex.lab=0.7,cex.axis=0.7)
        for(m in 1:length(levels(Groups))){
          barplot(t(ESx[[m]]),beside = T,horiz = T,las=1,cex.name=0.7,xlim=c(0,12),col = barcol,border = F)
          abline(v=-log10(alpha),lty=3)
          legend("bottomright",legend = colnames(ESx[[m]]),
                 bty="n",fill=barcol,cex = 0.7)
        }
      }
    }
    
    
    return(ESx)
  }
}


# KO to metabolome
KO2metabolite <- function(data.p){
  kegg.rx <- readRDS(url("http://www.lifescipy.net/RcodeDB/KEGG/KO-RN-CPD.db.RDS"))
  kegg.rx <- do.call(rbind,kegg.rx)
  
  cpd.w <- list()
  for(i in 1:nrow(data.p)){
    rx.i <- kegg.rx[which(kegg.rx$KO %in% as.character(data.p[i,1])),]
    if(nrow(rx.i)>1){
      rx.i <- rx.i[!duplicated(paste(rx.i$rn,rx.i$KO)),]
      rx.i[which(rx.i$type=="irreversible"),c("substrate.1","substrate.2","substrate.3","substrate.4","substrate.5")] <- NA
      
      cpd.i <- data.frame(CID=matrix(t(rx.i[,5:14]),ncol = 1),weight=NA)
      cpd.i <- cpd.i[!is.na(cpd.i$CID),]
      if(nrow(cpd.i)>1){
        subject.i <- matrix(apply(!is.na(rx.i[,5:9]),2,as.numeric),ncol = 5)
        product.i <- matrix(apply(!is.na(rx.i[,10:14]),2,as.numeric),ncol = 5)
        
        subject.w <- subject.i/rowSums(subject.i)
        subject.w[is.nan(subject.w)] <- 0
        subject.w <- subject.w / (as.numeric(rowSums(product.i)>0) + as.numeric(rowSums(subject.i)>0))
        product.w <- product.i/rowSums(product.i)
        product.w[is.nan(product.w)] <- 0
        product.w <- product.w  / (as.numeric(rowSums(product.i)>0) + as.numeric(rowSums(subject.i)>0))
        subject.w[subject.w==0] <- NA
        product.w[product.w==0] <- NA
        weight <- cbind(subject.w, product.w)
        weight <- matrix(t(weight),ncol = 1)
        cpd.i$weight <- weight[!is.na(weight),]
        cpd.w[[i]] <- aggregate(cpd.i$weight,list(cpd.i$CID),sum)
      }
    } else cpd.w[[i]] <- data.frame("Group.1"=NA,"x"=NA)
  }
  names(cpd.w) <- data.p$KO
  
  cpd.x <- cpd.w
  sample.names <- names(data.p)
  for(i in 1:length(cpd.x)){
    for(j in 2:length(sample.names)){
      if(!is.null(cpd.x[[i]])){
        cpd.x[[i]]$newcol <- cpd.x[[i]]$x * data.p[i,j]
        names(cpd.x[[i]])[j+1] <- sample.names[j]
      }
    }
  }
  
  cpd.x <- do.call(rbind,cpd.x)
  cpd.x <- aggregate(cpd.x[,-1],list(cpd.x$Group.1),sum)
  return(cpd.x)
}
