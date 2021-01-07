# MARco

**MARco** is an R code and function set for microbiome analysis with less pain to generate basic plots of microbial ecology, including alpha diversity boxplots, violin plots, beta diversity ordination (PCoA based on Bray-Curtis dissimilarity) with ADONIS test, and several statistical and visualization tools (see below). **MARco** also contains several discrete color sets for applying to your grouped plots. For advanced analysis, co-occurrence network analysis is also able to conduct by **MARco**.

All updated codes will be released on <http://www.lifescipy.net/RcodeDB/MARco.html>.

**Please cites:**  
Liu, P.-Y. _MARco: Microbiome Analysis RcodeDB_ (2020). Available at: <https://github.com/poyuliu/MARco/>. (Accessed: Day Month Year)


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
`randomForest`, `plyr`, `rfUtilities`, `caret`, `ROCR` and `pROC` for Random Forest analysis

### Available codes & how to use:
~~source("http://www.lifescipy.net/RcodeDB/FN_microbiome.R") # Microbiome analysis function set~~
~~source("http://www.lifescipy.net/RcodeDB/FN_colsets.R") # color sets~~
~~source("http://www.lifescipy.net/RcodeDB/FN_network.R") # network analysis function modules~~
~~source("http://www.lifescipy.net/RcodeDB/FN_RandomForest.R") # Random Forest function modules~~
~~source("http://www.lifescipy.net/RcodeDB/FN_enterotyping.R") # Enterotyping~~ 
~~source("http://www.lifescipy.net/RcodeDB/FN_KO2PathwayModule.R") # Reformat KEGG Orthology (KO) to pathways and modules from metagenome, metatranscriptome, or functional prediction data~~

install package from github:
```
library(devtools)
install_github("poyuliu/MARco")
``` 
if fail to install rgdal package automatically
<https://r-forge.r-project.org/R/?group_id=884>
```
 install.packages("rgdal", repos="http://R-Forge.R-project.org")
```
sourcing MARco R library 
```
library(MARco)
```
