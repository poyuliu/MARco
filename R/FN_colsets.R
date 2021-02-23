colset.d.2 <- RColorBrewer::brewer.pal(11,"RdBu")[c(10,2)]
colset.d.4 <- c("#227093","#ffa801","#2c2c54","#b33939")
colset.d.5 <- c("#4285F4","#FBBC05","#EA4335","#0b5e12","#2e2128")
colset.d.5 <- c("#4285F4","#FBBC05","#EA4335","#0b5e12","#2e2128")
colset.d.6 <- c("#003C78","#FFBE00","#46961E","#645AA0","#D25500","#474747")
colset.d.12 <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')


darker.col <- function(colset,degree=30){
    colset.dark <- col2rgb(colset)-degree
    colset.dark <- ifelse(colset.dark<0,0,colset.dark)
    colset.dark <- ifelse(colset.dark>255,255,colset.dark)
    colset.dark <- apply(colset.dark,2,function(x) rgb(x[1],x[2],x[3],maxColorValue = 255))
    return(colset.dark)
}


WhBuRd.g <- colorRampPalette(c("#FFFFFF","#2e96ff", "#ff243f"))
BuWhRd.g <- colorRampPalette(c("#2e96ff","#FFFFFF", "#ff243f"))

#save.image("FN_colsets.RData")
