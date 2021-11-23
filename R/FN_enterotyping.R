#list.of.packages <- c("clusterSim","ade4")
#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages)
#
##required packages
#require(clusterSim)
#require(ade4)
#require(vegan)


dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}

pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  cluster = as.vector(cluster::pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

noise.removal <- function(dataframe, percent=0.01, top=NULL){
  dataframe->Matrix
  bigones <- rowSums(Matrix)*100/(sum(rowSums(Matrix))) > percent 
  Matrix_1 <- Matrix[bigones,]
  print(percent)
  return(Matrix_1)
}

enterotyping <- function(otu,remove.noise=TRUE,rm.pct=0.01,kx=NULL){
  if(remove.noise==TRUE){
    otu <- noise.removal(dataframe = otu,percent = rm.pct)
  } 
  data.dist=dist.JSD(otu)
  
  nclusters=NULL
  require(clusterSim)
  require(cluster)
  for (k in 1:20) { 
    if (k==1) {
      nclusters[k]=NA 
    } else {
      data.cluster_temp=pam.clustering(data.dist, k)
      nclusters[k]=clusterSim::index.G1(t(otu),data.cluster_temp,  d = data.dist,
                                        centrotypes = "medoids")
    }
  }
  
  best.k <- which.max(nclusters)
  print(best.k)
  plot(nclusters, type="h", xlab="k clusters", ylab="CH index")
  if(is.null(kx)){
    best.k <- best.k
  } else best.k <- kx
   
  data.cluster=pam.clustering(data.dist, k=best.k)
  
  return(list(otu=otu,data.dist=data.dist,data.cluster=data.cluster))
}

plot.ET <- function(enterotyping,colset,...){
  #ordNtest <- function(data,group,dist=NULL){
  #  if(is.null(dist)){
  #    bc.dist <- vegan::vegdist(t(data))
  #  } else if(!is.null(dist)){
  #    bc.dist <- dist
  #  }
  #  pcoa <- cmdscale(bc.dist,eig = T)
  #  set.seed(123); permanova <- adonisx(bc.dist~group)
  #  return(list(pcoa=pcoa,permanova=permanova))
  #}
  #plot.ordi <- function(ordN,group,colset,ordiellipse=FALSE,ordispider=FALSE,plot=TRUE,pch=19,legend=NULL,test=FALSE,ordilabel=NULL){
  #  ord <- ordN$pcoa
  #  pct <- round(prop.table(ord$eig[ord$eig>0])*100,2)
  #  par(mar=c(5.1,4.1,4.1,10.1))
  #  if(length(pct)>0){
  #    plot(vegan::scores(ord),type="n",las=1,
  #         xlab=paste0("PCoA1 (",round(pct[1],2),"%)"),ylab=paste0("PCoA2 (",round(pct[2],2),"%)"))
  #  } else if(length(pct)==0) plot(vegan::scores(ord),type="n",las=1,xlab="Dim1",ylab="Dim2")
  #  if(ordiellipse==TRUE) vegan::ordiellipse(ord,group,col = colset,lwd = 0.5)
  #  if(ordispider==TRUE) vegan::ordispider(ord,group,col = colset,lwd = 0.5)
  #  if(plot==TRUE & length(pch)>1){
  #    points(vegan::scores(ord),pch=pch[as.numeric(group)],cex=1,col=paste0(colset,"aa")[as.numeric(group)])
  #  } else if(plot==TRUE & length(pch)==1) points(vegan::scores(ord),pch=pch,cex=1,col=paste0(colset,"aa")[as.numeric(group)])
  #  if(!is.null(legend)) legend("topright",legend =legend ,pch=pch,col=paste0(colset,"aa"),bty="n",cex=0.8,xpd=T,inset = c(-0.3,0))
  #  if(test==TRUE) {
  #    p <- format(ordN$permanova$aov.tab$`Pr(>F)`[1],digits = 2)
  #    r2 <- format(ordN$permanova$aov.tab$R2[1],digits = 2)
  #    legend("bottomright",legend =bquote(R^2 == .(r2)),bty="n",cex=0.8,xpd=T,inset = c(-0.2,0.06))
  #    legend("bottomright",legend =bquote(italic(p == .(p))),bty="n",cex=0.8,xpd=T,inset = c(-0.2,0))
  #  }
  #  if(!is.null(ordilabel)){
  #    vegan::ordilabel(scores(ordN$pcoa),cex = 0.7,fill = "#cacaca55",border = 1,col = "#1b262c",
  #                     select = which(rownames(scores(ordN$pcoa)) %in% ordilabel))
  #  }
  #}
  
  pcoa <- ordNtest(data = enterotyping$otu,
                   group = as.factor(enterotyping$data.cluster),
                   dist = enterotyping$data.dist)
  
  plot.ordi(ordN = pcoa,group = as.factor(enterotyping$data.cluster),
            colset = colset,legend = paste("cluster",levels(as.factor(enterotyping$data.cluster))),...)
  
}

box.ET <- function(enterotyping,clusterN,colset=NULL,out=FALSE,...){
  jitter.points <- function(data,groups,colset){
    for(i in 1:length(groups)){
      set.seed(123)
      points(jitter(rep(i,sum(as.numeric(groups)==i)),8/(2^(i-1))),data[as.numeric(groups)==i],pch=21,bg=colset[i],cex=1.1)
    }
  }
  ETn <- enterotyping$otu[,enterotyping$data.cluster==clusterN]
  MAX <- which.max(prop.table(rowSums(ETn)))
  MAX.OTU <- (sort(prop.table(rowSums(ETn)),decreasing = T)[1])*100
  print(MAX.OTU)
  ETotu <- as.data.frame(prop.table(t(enterotyping$otu),1))
  boxplot(ETotu[,MAX]~as.factor(enterotyping$data.cluster),
          las=1,col=paste0(colset,"99"),
          xlab="Enterotypes",
          ylab="Relative abundance",...)
  jitter.points(ETotu[,MAX],as.factor(enterotyping$data.cluster),colset)
  
  if(out==TRUE){
    return(ETn)
  }
}