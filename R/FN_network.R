require(SpiecEasi)

###### FUNCTION 01 #####
# > prevalence(data,prevalence)
# filter feature table by prevalence (in 0-1 decimal)
prevalence <- function(data,prevalence=0.1){
  prev <- apply(data,1,function(x) sum(x>0)/ncol(data))
  data <- data[prev>=prevalence,]
  return(data)
}

##### FUNCTION 02 #####
# > corr(data,cutoff,method)
# establish a feature correlation matrix from a feature table
# SparCC is default method; "SpiecEasi" package is required
corr <- function(data,cutoff=0,method="sparcc"){
  data <- data.frame(t(data))
  if(method=="sparcc"){
    w <- SpiecEasi::sparcc(data)[[2]]
    colnames(w) <- rownames(w) <- colnames(data)
  } else  w <- cor(data,method = method)
  diag(w) <- 0
  w[w < cutoff] <- 0
  return(w)
}

##### FUNCTION 02.1 #####
# use after sparccboot function
boot2matrix <- function(bootcors,otutable){
  mempty <- matrix(0, nrow = nrow(otutable), ncol = nrow(otutable))
  mindex <- matrix(1:c(nrow(otutable)^2), nrow =  nrow(otutable), ncol =  nrow(otutable))
  mempty[mindex[upper.tri(mindex)]] <- bootcors
  mempty[lower.tri(mempty)] <- t(mempty)[lower.tri(t(mempty))]
  diag(mempty) <- 1
  return(mempty)
}

##### FUNCTION 02.2 #####
#network matrix
network.pipeline <- function(otu,prevalence=0.1,alpha=0.05,bootstrap=99,cpu=1){
  net.o <- prevalence(data = otu,prevalence = prevalence)
  ww <- SpiecEasi::sparccboot(t(net.o),R=bootstrap,ncpus = cpu) # R: bootsrtapping numbers
  w.boot <- SpiecEasi::pval.sparccboot(ww)
  w.mb <- w.boot$cors
  p.mb <- w.boot$pvals
  w.mb <- boot2matrix(w.mb,net.o)
  p.mb <- boot2matrix(p.mb,net.o)
  diag(w.mb) <- 0
  w.mb[p.mb > alpha] <- 0 # significance cutoff
  rownames(w.mb) <- colnames(w.mb) <- rownames(net.o)
  return(w.mb)
}

##### FUNCTION 03 #####
# > make_network(w,data,cutoff,size,degree,clustering)
# make network using igraph package
make_network <- function(w,data,cutoff=0.5,both=TRUE,size=1,degree=1,clustering=TRUE){
  require(igraph)
  #g.x <- graph.adjacency(w, weighted=TRUE,mode="lower")
  if(both==TRUE){
    g.x <- graph.adjacency(abs(w), weighted=TRUE,mode="lower")
    g.y <- graph.adjacency(w, weighted=TRUE,mode="lower")
    #g.x <- delete.edges(g.x, which(abs(E(g.x)$weight) < cutoff))
    g.x <- delete.edges(g.x, which(E(g.x)$weight < cutoff))
    g.x <- delete.vertices(g.x,which(degree(g.x) <= degree)) 
    g.y <- delete.edges(g.y, which(!E(g.y) %in% E(g.x)))
    w.sign <- E(g.y)$weight #
    w.sign <- w.sign/abs(w.sign) #
  } else if(both==FALSE){
    g.x <- graph.adjacency(w, weighted=TRUE,mode="lower")
    if(cutoff > 0){
      g.x <- delete.edges(g.x, which(E(g.x)$weight < cutoff))
      g.x <- delete.vertices(g.x,which(degree(g.x) <= degree)) #
      w.sign <- E(g.x)$weight #
      w.sign <- w.sign/abs(w.sign) #
    } else if(cutoff <0){
      g.x <- delete.edges(g.x, which(E(g.x)$weight > cutoff))
      E(g.x)$weight <- abs(E(g.x)$weight) #
      g.x <- delete.vertices(g.x,which(degree(g.x) <= degree)) #
      w.sign <- E(g.x)$weight #
      w.sign <- w.sign/abs(w.sign) #
    }
  }
  V(g.x)$size <- log10(rowSums(data)[match(V(g.x)$name,rownames(data))])*size
  #g.x <- delete.vertices(g.x,which(degree(g.x) <= degree))
  if(clustering==TRUE){
    g.xx <- g.x
    E(g.xx)$weight <- abs(E(g.xx)$weight)
    set.seed(100)
    cluster <-cluster_fast_greedy(g.xx)
    clusters <- membership(cluster)
    #sizes(cluster)
    #V(g.x)$color <- c(cols,rep("grey50",length(sizes(cluster))-3))[clusters]
  }
  set.seed(123)
  layout.x <- layout_with_fr(g.x)
  layout.x <- layout.norm(layout.x, ymin=-1, ymax=1, xmin=-1, xmax=1)
  E(g.x)$weight <- E(g.x)$weight * w.sign
  if(clustering==TRUE){
    output <- list(graph=g.x,layout=layout.x,cluster.out=cluster,clusters=clusters)
  } else output <- list(graph=g.x,layout=layout.x)
  return(output)
}

##### FUNCTION 03.1 #####
# > plot_network(graph_from_make_network,positive=TRUE,negative=TRUE,pos.col="#0652DD55",neg.col="#EA202755")
plot_network <- function(graph_from_make_network,positive=TRUE,negative=TRUE,pos.col="#0652DD55",neg.col="#EA202755",curve=0,frame.color="#FAFAFA"){
  require(igraph)
  if(positive==TRUE && negative==TRUE){
    E(graph_from_make_network$graph)$color <- ifelse(E(graph_from_make_network$graph)$weight>0,pos.col,neg.col)
  } else if(positive==TRUE && negative==FALSE){
    E(graph_from_make_network$graph)$color <- pos.col
  } else if(positive==FALSE && negative==TRUE){
    E(graph_from_make_network$graph)$color <- neg.col
  }
  plot(graph_from_make_network$graph,
       vertex.label=NA,
       rescale=F,
       layout = graph_from_make_network$layout, 
       vertex.frame.color= frame.color,edge.curved=curve
  )
}

##### FUNCTION 05 #####
# > scatter_hull(x,y,)
scatter_hull <- function(xcoord, ycoord, lcolor){
  hpts <- chull(x = xcoord, y = ycoord)
  hpts <- c(hpts, hpts[1])
  #lines(xcoord[hpts], ycoord[hpts], col = lcolor)
  polygon(xcoord[hpts], ycoord[hpts],border = lcolor, col=paste0(lcolor,"22"),lty=3)
}  
