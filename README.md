# MARco

**MARco** is an R code and function set for microbiome analysis with less pain to generate basic plots of microbial ecology, including alpha diversity boxplots, violin plots, beta diversity ordination (PCoA based on Bray-Curtis dissimilarity) with ADONIS test, and several statistical and visualization tools (see below). **MARco** also contains several discrete color sets for applying to your grouped plots. For advanced analysis, co-occurrence network analysis is also able to conduct by **MARco**.

Tutorial and all updated codes will be released on <http://www.lifescipy.net/MARco/MARco.html>.

**Please cites:**  
Liu, P.-Y. _poyuliu/MARco: MARco: Microbiome Analysis RcodeDB (Version v1.0)_. Zenodo. Available online: http://doi.org/10.5281/zenodo.4589897 (accessed on Day Month Year).

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
`doParallel` for Windows version Random Forest analysis (`runRFwin()`)  

### Available codes & how to use:
1. Microbiome analysis function set
2. color sets
3. network analysis function modules
4. Random Forest function modules
5. Enterotyping 
6. Reformat KEGG Orthology (KO) to pathways and modules from metagenome, metatranscriptome, or functional prediction data  


install dependencies from Bioconductor  
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("gage","DESeq2"))
```  
install dependencies from GitHub  
```
library(devtools)
install_github("zdk123/SpiecEasi")
```  
if fail to install rgl package automatically  
please refer to: <https://r-forge.r-project.org/R/?group_id=234>
```
 install.packages("rgl", repos="http://R-Forge.R-project.org")
```

install package from github:
```
library(devtools)
install_github("poyuliu/MARco")
``` 
sourcing MARco R library 
```
library(MARco)
```
