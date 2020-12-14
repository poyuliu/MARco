---
title: 'MARco: Microbiome Analysis RcodeDB'
author: "Liu, Po-Yu, Ph.D., National Taiwan University College of Medicine"
output:
  html_document:
    toc: yes
    toc_float:
      collapsed: no
      smooth_scroll: yes
  word_document:
    toc: yes
---

<style type="text/css">
.main-container {
  max-width: 1800px;
  margin-left: auto;
  margin-right: auto;
}
</style>

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# MARco

**MARco** is an R code and function set for microbiome analysis with less pain to generate basic plots of microbial ecology, including alpha diversity boxplots, violin plots, beta diversity ordination (PCoA based on Bray-Curtis dissimilarity) with ADONIS test, and several statistical and visualization tools (see below). **MARco** also contains several discrete color sets for applying to your grouped plots. For advanced analysis, co-occurrence network analysis is also able to conduct by **MARco**.

All updated codes will be released on <http://www.lifescipy.net/RcodeDB/MARco.html>.

**Please cites:**  
Liu, P.-Y. _MARco: Microbiome Analysis RcodeDB_ (2020). Available at: <http://www.lifescipy.net/RcodeDB/MARco.html>. (Accessed: Day Month Year)

**Latest updates:** `r format(Sys.Date(), "%B %d, %Y")`

### Required packages:
`vegan` with functions for most community ecology analyses  
`fossil` for calculating Chao1 index  
`vioplot` for ploting violin plots  
`igraph` for constructing network  
`SpiecEasi` for calculating SparCC correlation, how to install: <https://github.com/zdk123/SpiecEasi>  
`DESeq2` RNA-Seq based differential expression test, see details on Bioconductor: <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>  
`pheatmap` for plot heatmap  
`MASS` for LDA (linear discriminant analysis) effect size (LEfSe)  
`gage` for GAGE (Generally Applicable Gene-set Enrichment) analysis  
`cluster` and `clusterSim` for Enterotyping

### Available codes & how to use:
```{r source, message=FALSE,warning=FALSE}
source("http://www.lifescipy.net/RcodeDB/FN_microbiome.R") # Microbiome analysis function set
source("http://www.lifescipy.net/RcodeDB/FN_colsets.R") # color sets
source("http://www.lifescipy.net/RcodeDB/FN_network.R") # network analysis function modules
source("http://www.lifescipy.net/RcodeDB/FN_RandomForest.R") # Random Forest function modules
source("http://www.lifescipy.net/RcodeDB/FN_enterotyping.R") # Enterotyping 
source("http://www.lifescipy.net/RcodeDB/FN_KO2PathwayModule.R") # Reformat KEGG Orthology (KO) to pathways and modules from metagenome, metatranscriptome, or functional prediction data
```

### Colorsets
```{r colorsets, echo=FALSE}
par(xaxt='n',yaxt="n",bty="n",mar=c(0.1,15,0.1,0),mfrow=c(6,1))
image(matrix(1:2),col=colset.d.2)
text(seq(0,length(colset.d.2),length.out = 2)/length(colset.d.2),0,labels = colset.d.2,srt=90,col=0)
text(-0.75,0,labels = 'colset.d.2',font = 2,cex=1.5,xpd=T)
image(matrix(1:4),col=colset.d.4)
text(seq(0,length(colset.d.4),length.out = 4)/length(colset.d.4),0,labels = colset.d.4,srt=90,col=0)
text(-0.33,0,labels = 'colset.d.4',font = 2,cex=1.5,xpd=T)
image(matrix(1:5),col=colset.d.5)
text(seq(0,length(colset.d.5),length.out = 5)/length(colset.d.5),0,labels = colset.d.5,srt=90,col=0)
text(-0.275,0,labels = 'colset.d.5',font = 2,cex=1.5,xpd=T)
image(matrix(1:12),col=colset.d.12)
text(seq(0,length(colset.d.12),length.out = 12)/length(colset.d.12),0,labels = colset.d.12,srt=90)
text(-0.175,0,labels = 'colset.d.12',font = 2,cex=1.5,xpd=T)
image(matrix(1:12),col=darker.col(colset.d.12,degree = 30))
text(0.5,0,labels = 'darker.col(colset,degree = 30)')
text(-0.175,0,labels = 'darker.col()',font = 2,cex=1.5,xpd=T)
image(matrix(1:30),col=WhBuRd.g(30))
text(0.5,0,labels = 'WhBuRd.g(n)')
text(-0.14,0,labels = 'WhBuRd.g()',font = 2,cex=1.5,xpd=T)
```

# Examples
## Regular analyses
Example dataset:  
_Suzuki T., et al._ **Host genetic determinants of the gut microbiota of wild mice.** Molecular ecology 28 (13), 3197-3207 (2019)  

### Alpha diversity
```{r suzuki, echo=FALSE}
load("example_data.RData")
```
```{r suzuki.alpha, warning=F,fig.width=9}
alpha.diversity <- alphdiv(data = otu,group = groups,test = T) # calculate alpha diversity indices (Richness, Shannon, Simpson, Chao1)
par(mfrow=c(1,4))
plot.alpha(alpha = alpha.diversity,group = groups,colset = colset.d.4,names = levels(groups),test = T) # output boxplots
vioplot.alpha(alpha = alpha.diversity,group = groups,colset = colset.d.4,names = levels(groups),test = T) # output violin plots
```

### Relative abundance barplots
```{r barplot.1, echo=TRUE,message=FALSE,fig.height = 5, fig.width=10}
taxaQ2output <- read.delim("suzuki_taxonomy/taxonomy.tsv")
TaxaLvABar(otu = otu,taxafromQ2tsv = taxaQ2output,level = 2,colsets = colset.d.12)
TaxaLvABar(otu = otu,taxafromQ2tsv = taxaQ2output,level = 3,colsets = colset.d.12)
TaxaLvABar(otu = otu,taxafromQ2tsv = taxaQ2output,level = 5,colsets = colset.d.12)
```
```{r barplot.2, echo=TRUE,message=FALSE,fig.height = 5, fig.width=9}
par(mfrow=c(1,3))
TaxaLvABar(otu = otu,taxafromQ2tsv = taxaQ2output,level = 2,colsets = colset.d.12,grouped = TRUE,groups=groups)
TaxaLvABar(otu = otu,taxafromQ2tsv = taxaQ2output,level = 3,colsets = colset.d.12,grouped = TRUE,groups=groups)
TaxaLvABar(otu = otu,taxafromQ2tsv = taxaQ2output,level = 5,colsets = colset.d.12,grouped = TRUE,groups=groups)
```

### Beta diversity
```{r suzuki.beta.1, warning=FALSE,echo=FALSE,message=FALSE}
filt <- ceiling(mean(colSums(otu))*0.001)
otu <- otu[rowSums(otu) >= filt,]
```
```{r suzuki.beta.2, warning=FALSE,message=F,fig.height=6}
library(vegan) # "vegan" is required!
pcoa <- ordNtest(data = otu,group = groups)
plot.ordi(ordN = pcoa,group = groups,colset = colset.d.4,
          ordiellipse = T,ordispider = T,
          legend = levels(groups),test = T)
```

### DESeq2 of 2-group comparison

#### (Statistic/Taxonomy tables &) Volcano plot (Fold change >2, adjust P < 0.05)  
```{r DESeq2.1, echo=FALSE,message=FALSE}
otu.2 <- otu[,groups!="NY"]
otu.2 <- otu.2[rowSums(otu.2)!=0,]
group.2 <- as.factor(as.character(groups)[groups!="NY"])
group.2 <- group.2[colSums(otu.2)!=0]
otu.2 <- otu.2[,colSums(otu.2)!=0]
taxa <- tax.rearrange("suzuki_taxonomy/taxonomy.tsv")
```
```{r DESeq2.2, echo=TRUE,message=FALSE,fig.height = 5, fig.width=5}
library(DESeq2) # "DESeq2" is required!
DESeq2.FL.MP <- DESq2(data = otu.2,group = group.2)
#DESeq2.sig.table(DESeq.res = head(DESeq2.FL.MP),taxa = taxa ,caption = "MP vs. FL")
volcano.plot(DESeq.res = DESeq2.FL.MP,colset = colset.d.4)
```

#### Heatmap with _p_ values
```{r DESeq2Heatmaps, echo=TRUE,message=FALSE,fig.height = 10, fig.width = 18,out.width=1152}
G2.heatmap(otu.sig = otu.2[which(DESeq2.FL.MP$padj<0.0001),],groups = group.2,
           p.val = DESeq2.FL.MP$padj[which(DESeq2.FL.MP$padj<0.0001)],
           taxa.sig = taxa[match(rownames(otu.2),rownames(taxa)),][which(DESeq2.FL.MP$padj<0.0001),],
           otu = otu.2,metadata = as.data.frame(cbind(SampleID=colnames(otu.2))),colset = colset.d.4[1:2])
```

### Co-occurrence network
```{r network.1, echo=TRUE,message=FALSE,fig.height = 7, fig.width=7}
library(igraph)
library(SpiecEasi)
source("http://www.lifescipy.net/RcodeDB/FN_network.R") # network analysis function modules
```
```{r network.2, echo=TRUE,eval=FALSE}
# prepare data and construct a correlation matrix with SparCC
net.o <- prevalance(data = otu,prevalance = 0.4) # prevalenve ≥ 0.4
ww <- sparccboot(t(net.o),R=9,ncpus = 6) # R: bootsrtapping numbers
w.boot <- pval.sparccboot(ww)
w.mb <- w.boot$cors
p.mb <- w.boot$pvals
w.mb <- boot2matrix(w.mb,net.o)
p.mb <- boot2matrix(p.mb,net.o)
diag(w.mb) <- 0
w.mb[p.mb > 0.05] <- 0 # significance cutoff
rownames(w.mb) <- colnames(w.mb) <- rownames(net.o)
```  
or just run a one-step function:  
```{r network.2-2, echo=TRUE,eval=FALSE}
w.mb <- network.pipeline(otu=net.o,prevalence=0.1,alpha=0.05,bootstrap=99,cpu=12)
```
```{r network.3, echo=FALSE,message=FALSE}
load("example_net.RData")
```
```{r network.4, echo=TRUE,message=FALSE,fig.height = 6, fig.width=6}
# construct and plot network
g.mb <- make_network(w = w.mb,
                     data = net.o,
                     both = F, # only positive correlation
                     cutoff = quantile(w.mb[w.mb>0],0.95), # correlation coefficients / weights cutoff for removing nodes
                     degree = 0, # connection degree cutoff for revoming nodes
                     size = 1.5,
                     clustering = T)
V(g.mb$graph)$color <- c(colset.d.12,"#cacaca")[g.mb$clusters] # label clusters with color (13 clusters)
plot_network(g.mb,negative = F,pos.col="#5d5d5d55")
```

## Advanced analyses  
Example dataset:  
_Liu P.-Y., et al._ **Variations in Gut Microbiota of Siberian Flying Squirrels Correspond to Seasonal Phenological Changes in Their Hokkaido Subarctic Forest Ecosystem.** Microbial ecology 78, 223–231 (2019)  

### LEfSe
Linear discriminant analysis (LDA) Effect Size (version 1.0; only for 2-group comparison)
```{r LDA load data, echo=FALSE, message=FALSE}
library(MASS)
load("example_PVo.RData",envir = loaddata <- new.env())
pvo.otu <- loaddata$otu.c97[,7:25]
pvo.taxa <- loaddata$taxa.c97[which(rowSums(pvo.otu)>0),]
pvo.taxa <- apply(pvo.taxa,2,function(x) gsub("D_.__","",as.character(x)))
pvo.otu <- pvo.otu[which(rowSums(pvo.otu)>0),]
meteor <- loaddata$meteor
meteor$hot <- as.factor(ifelse(meteor$Temperature>20,"hot","cool"))
rm(loaddata)
```
```{r LEfSe, message=FALSE,fig.height=5,fig.width=10, out.width=1152}
Lefse.temp <- lefse(data = pvo.otu,taxa = pvo.taxa,Groups = meteor$hot,Lefse.factor = 1,FDR = FALSE,alpha = 0.1,LDAscore = 2,plot = TRUE,col = colset.d.2)

print(Lefse.temp$LEfSe)
```


### Random Forest
```{r RF, echo=FALSE, message=FALSE}
load("example_PVo.RData")
RF <- readRDS("RFout.rds")
RF2 <- readRDS("RF2out.rds")
rownames(pvo.otu) <- paste0("ASV_",rownames(pvo.otu)) # A number of the first character of feature name is not allowed 
rownames(otu) <- paste0("ASV_",rownames(otu))
meteor$hot <- as.factor(ifelse(meteor$Temperature>20,"hot","cool"))
```
Random Forest regression

Random Forest classification for model prediction and feature selection
```{r RF classification, echo=TRUE,message=FALSE,eval=FALSE}
rownames(pvo.otu) <- paste0("ASV_",rownames(pvo.otu)) # A number of the first character of feature name is not allowed 
rownames(otu) <- paste0("ASV_",rownames(otu))
meteor$hot <- as.factor(ifelse(meteor$Temperature>20,"hot","cool"))

RF <- runRF(otu_table_scaled = log1p(pvo.otu),Y = hot,metadata = meteor,model = "classification")
```
```{r RF AUROC, echo=TRUE,message=FALSE}
AUROC <- RF2ROC(RF,plot=T)
```

```{r RF2,eval=FALSE}
pvo.otu.2 <- pvo.otu[match(RF$SortFeature$features[1:10],rownames(pvo.otu)),] # select top 10 features
RF2 <- runRF(otu_table_scaled = log1p(pvo.otu.2),Y = hot,metadata = meteor,model = "classification")
```
```{r RF AUROC2, echo=TRUE,message=FALSE}
AUROC2 <- RF2ROC(RF2,plot=T)
```

Validation cohort cross validation
```{r RF validation, echo=TRUE,message=FALSE}
newd <- log1p(otu[,c(1:6,26:29)])
newdgroup <- as.factor(c(rep("hot",6),rep("cool",4)))
validation.out <- VCCV(RFmodel = RF2,newdata = newd,newdataGroup = newdgroup,RFROC = AUROC2,plot = T)
```


### Longitudinal analysis

### Enterotyping (PAM method)
Reference: Gomez, A. et al. Plasticity in the Human Gut Microbiome Defies Evolutionary Constraints. mSphere 4, doi:10.1128/mSphere.00271-19 (2019).  
Tutorial: <https://enterotype.embl.de/enterotypes.html>
```{r enterotyping, echo=TRUE,eval=FALSE}
source("http://www.lifescipy.net/RcodeDB/FN_enterotyping.R") # Enterotyping
ET <- enterotyping(otu = genus.otu,remove.noise = T,rm.pct = 0.05)
plot.ET(enterotyping = ET,colset = colset.d.4,ordiellipse = T,ordispider = T,test = T)
ET1 <- box.ET(enterotyping = ET,clusterN = 1,out = T,colset = colset.d.4,outline=F)
ET2 <- box.ET(enterotyping = ET,clusterN = 2,out = T,colset = colset.d.4,outline=F)
```  
Visualize enterotypes using network



### Fitting variable gradients   
Distance-based redundancy analysis (db-RDA) principle and examples: <https://archetypalecology.wordpress.com/2018/02/21/distance-based-redundancy-analysis-db-rda-in-r/>  
```{r PVO, echo=FALSE}
load("example_PVo.RData")
```
```{r PVO.ENV}
head(pvo.otu)[,1:7] # view gut microbiota OTU table of Siberian flying squirrels
head(meteor) # view environmental factor metadata
```  

**dbRDA with continuous variable gradients fitting**
```{r dbRDA.PCoA+ENV, echo=TRUE,message=FALSE,fig.height = 6.25, fig.width=9.2}
dbRDA <- dbrdaNtest(data = pvo.otu,ENV = meteor)
print(dbRDA$permucapscale)
plot.dbRDA(ordN = dbRDA, select = meteor$Temperature, colset = WhBuRd.g(20), fit.p = 0.1, legend = TRUE,legend.labels = "Temperature (C)")
terrain20 <- substr(terrain.colors(20,rev = T),1,7) # color for NDVI
plot.dbRDA(ordN = dbRDA, select = meteor$NDVI, colset = terrain20, fit.p = 1, legend = TRUE, legend.labels = "Mean NDVI")
```

**PCoA with continuous variable gradients fitting**
```{r PVO.PCoA+ENV, echo=TRUE,message=FALSE,fig.height = 6.25, fig.width=9.2}
pvo <- ordEnvfit(data = pvo.otu,ENV = meteor) # compute PCoA and environmental factor fitting
print(pvo$fit)
plot.ordi.env(ordN = pvo,select = meteor$Temperature,colset = WhBuRd.g(20),fit.p = 0.1,legend = TRUE,legend.labels = "Temperature (C)")

terrain20 <- substr(terrain.colors(20,rev = T),1,7) # color for NDVI
plot.ordi.env(ordN = pvo,select = meteor$NDVI,colset = terrain20,fit.p = 1,legend = TRUE,legend.labels = "Mean NDVI")
```

### OTU correlations
**Simple Env correlation**
```{r simpleEnvCor, echo=TRUE, message=FALSE, fig.width=9.2}
Temperature <-  meteor$Temperature
par(mfrow=c(1,2))
SimpleEnvCor <- EnvCorOTU(data = pvo.otu,ENV = Temperature,summary.plot = T)
SimpleEnvCor.loglog <- EnvCorOTU(data = pvo.otu,ENV = Temperature,logdata = TRUE,logENV = TRUE,summary.plot = T)
head(SimpleEnvCor$otu.cor.env)[,1:7]
taxa[match(rownames(SimpleEnvCor$otu.cor.env),rownames(taxa)),][1:3,]
head(SimpleEnvCor$sig.corr)
```
**Coordinates x Env correlation**  
```{r coordinateEnvCor.1, echo=TRUE, message=FALSE, fig.width=6}
CoordinateEnvCor <- PCENVcorOTU(data=pvo.otu,Ord=pvo,ENV=Temperature,plot.pc.env = TRUE)
head(CoordinateEnvCor$otu.cor.pc,n = 3)
head(CoordinateEnvCor$sig.corr,n = 5)
head(CoordinateEnvCor$sig.p,n = 5)
```

```{r coordinateEnvCor.2, echo=TRUE, message=FALSE, fig.width=6}
CoordinateEnvCor <- PCENVcorOTU(data=pvo.otu,Ord=pvo,ENV=Temperature,summary.plot = T)
```

# Functional analysis
## KO2pathway
Reformat KEGG Orthology (KO) annoated data to pathways or modules

**Reformat**

(example) original data:
```{r picrust.load, echo=FALSE}
picrust2 <- readRDS("example_PICRUSt2.RDS")
```
```{r picrust.print}
head(picrust2) # functional data input format
```
Convert to metabolic pathways
```{r picrust.met}
gs.table <- KO2path(predtable = picrust2,GeneSet = "metabolic") # use "metabolic pathway" gene sets
# other gene set options: "pathway" [all KEGG pathways], "module" [KEGG modules]
head(gs.table)
gs.table.p <- KO2path(predtable = picrust2,GeneSet = "metabolic",prop = TRUE) # present with proportion
head(gs.table.p)
```
Convert to KEGG modules
```{r picrust.mod}
md.table <- KO2path(predtable = picrust2,GeneSet = "module")
head(md.table)
```
**Annotation**

Annotate with upper levels
```{r picrust.lvs.1}
lv.table <- path.lvs(gs.table)
head(lv.table)
```
Summarize to level 2
```{r picrust.lvs.2}
lv2.table <- path.lvs(gs.table,aggregate.lv = 2)
head(lv2.table)
```

## FSEA  
**Functional Sets Enrichment Analysis**  
Methodology reference:  
_Liu, P-Y, et al._ **Body-size Scaling is Related to Gut Microbial Diversity, Metabolism and Dietary niche of Arboreal folivorous flying Squirrels.** Scientific reports 10, 1-12 (2020)

```{r gage.load, echo=FALSE}
library(gage)
picrust2 <- readRDS("example_PICRUSt2.RDS")
tax4fun <- readRDS("example_4FS_Tax4Fun.RDS")
```

**__Within__ group test**
```{r gage.within}
Clutch <- factor(c("clutch.1","clutch.1","clutch.2","clutch.2"))
FSEA.metabolic <- FSEA(Fdata = picrust2,Groups = Clutch,FunctionalSets = "metabolic",SortLv = 2,Test = "within",barcol = colset.d.4[1:2],plot = TRUE)
FSEA.module <- FSEA(Fdata = picrust2,Groups = Clutch,FunctionalSets = "module",SortLv = 3,Test = "within",barcol = colset.d.4[1:2],plot = TRUE)
```

**__Between__ group test**
```{r gage.between, fig.height=10}
FS4gp <- substr(colnames(tax4fun),1,3)
FS4gp <- gsub("PVO","Small FS",FS4gp)
FS4gp <- gsub("TX0","Medium FS",FS4gp)
FS4gp <- gsub("PAL","Large FS",FS4gp)
FS4gp <- gsub("PPG","Large FS",FS4gp)
FS4gp <- as.factor(FS4gp)

par(mfcol=c(3,1))
FS.FSEA.between.metabolic <- FSEA(Fdata = tax4fun,Groups = FS4gp,FunctionalSets = "metabolic",Test = "between",plot = T,barcol = colset.d.2)
```

## Predict metabolome profile
```{r predict_metabolome_code, eval=FALSE}
KO2metabolite(petaurista)
```
```{r predict_metabolome_example, echo=FALSE,fig.width=7,fig.height=7.5}
source("example_predict_metabolome.R")
```

## Reconstruct metabolic networks (reactions)