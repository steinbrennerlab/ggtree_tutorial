---
title: "custom ggtree"
output:
  github_document:
    toc: true
    toc_depth: 2
---
## About the scripts
Adam Steinbrenner <br>
astein10@uw.edu <br>
http://steinbrennerlab.org <br>
Updated 1/29/2019 <br>
<br>
The following is the Steinbrenner lab ggtree pipeline for phylogenetics tree visualization incorporating other attributes and heatmaps

---

## Download ggtree.R and data
Find ggtree.R in main directory, and all input files in /data/ggtree

---

## Install packages
The function tree_subset may require the development version of treeio; make sure it is current from github
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.10")

BiocManager::install("EBImage")
BiocManager::install("treeio")
BiocManager::install("ggtree")

install.packages("phytools")
install.packages("optparse")
install.packages("tidyselect")
install.packages("tidyselect")
install.packages("labeling")
install.packages("devtools")
devtools::install_github("GuangchuangYu/treeio")
```


---

## Load libraries
```{r echo=FALSE, message=FALSE, warning=FALSE}
library("ggplot2")
library("treeio")
library("phytools") # for sims and ASRs
library("EBImage") # for images
library("ggtree")
library("optparse")
```
sessionInfo()


---

## Usage and Options
The script ggtree.R can be run separately with options
e.g.
`Rscript ggtree.R --entry fig1_rlck --write test28 --width 2 --height 3.2 -size 2 -line 0.2 --push 5 --label_offset  1.5 --symbol_size 0.5`

```
option_list <- list( 
    make_option(c("-e", "--entry"), action="store", type="character",
        help="entry"),
	make_option(c("-b", "--write"), action="store", type="character", 
        help="parameters term for file writing"),
	make_option(c("-n", "--node"), action="store", default=0, type="integer", 
        help="feed script a node and it will subset one node higher"),
	make_option(c("-f", "--width"), action="store", default=0, type="numeric", 
        help="width for pdf, optional"),
	make_option(c("-g", "--height"), action="store", default=0, type="numeric", 
        help="height for pdf, optional"),
	make_option(c("-i", "--line"), action="store", default=0.3, type="numeric", 
        help="line_width"),
	make_option(c("-x", "--push"), action="store", default=0.3, type="numeric", 
        help="push right side to make room"),
	make_option(c("-r", "--label_offset"), action="store", default=1, type="numeric", 
        help="offset for gene symbols"),
	make_option(c("-m", "--symbol_size"), action="store", default=1, type="numeric", 
        help="offset for gene symbols"),
	make_option(c("-z", "--size"), action="store", default=2, type="numeric", 
        help="size of font")
		)
message(option_list)
opt <- parse_args(OptionParser(option_list=option_list))
```

For the tutorial, set these values in Rstudio manually by running the following code chunk
```{r}
option_list <- list( 
    make_option(c("-e", "--entry"), action="store", type="character",
        help="entry"),
	make_option(c("-b", "--write"), action="store", type="character", default="output", 
        help="parameters term for file writing"),
	make_option(c("-n", "--node"), action="store", default=0, type="integer", 
        help="feed script a node and it will subset one node higher"),
	make_option(c("-f", "--width"), action="store", default=0, type="numeric", 
        help="width for pdf, optional"),
	make_option(c("-g", "--height"), action="store", default=0, type="numeric", 
        help="height for pdf, optional"),
	make_option(c("-i", "--line"), action="store", default=0.3, type="numeric", 
        help="line_width"),
	make_option(c("-x", "--push"), action="store", default=0.3, type="numeric", 
        help="push right side to make room"),
	make_option(c("-r", "--label_offset"), action="store", default=1, type="numeric", 
        help="offset for gene symbols"),
	make_option(c("-m", "--symbol_size"), action="store", default=1, type="numeric", 
        help="offset for gene symbols"),
	make_option(c("-z", "--size"), action="store", default=1, type="numeric", 
        help="size of font")
		)
message(option_list)
opt <- parse_args(OptionParser(option_list=option_list))
```

---

## Inputs: newick tree, attribute files, heatmap data.  Replace the directory with your own
Newick tree
```{r}
dir<-paste("C:/Users/Adam/Dropbox/github/ggtree_tutorial/ggtree_tutorial")
message(dir)
tree <- read.tree(paste(dir,"/ggtree_files/tree_input.nwk",sep=""))
```
Attribute files, tab delimited with columns taxa and atttribute
Species, specifying from which genome each gene entry came from
```{r}
dd <- read.table(paste(dir,"/ggtree_files/attribute_species.txt",sep=""), header = TRUE, stringsAsFactor=F)
```
HMM, an output of significant hits from an HMMer pipeline
```{r}
dd2 <- read.table(paste(dir,"/ggtree_files/attribute_hmm.txt", sep=""), header = TRUE, stringsAsFactor=F)
```
Heatmap, e.g. a list of log2 Fold Changes associated with the genes
```{r}
counts_file <- read.table(paste(dir,"/ggtree_files/log2FoldChanges.txt", sep=""), row.names = 1, header = TRUE, sep="\t", stringsAsFactor=F)
```
Set upper and lower limits for the heatmap
```{r}
upper <- 5
lower <- -5
```

---

## ggTree object
ggTree creates a tree visualization using ggplot and feeds in the various dataframes (dd, dd2).  You can call columns from these inputs to specify different aspects of the visualization
```{r}
p <- ggtree(tree, size=opt$line) #size specifies line size
###OPTIONAL: adds hmm coding to ggtree object q
###OPTIONAL: takes hmm_coding and adds to a separate dataframe dd2

p <- p %<+% dd2

#Takes the messy tree object and converts it to a dataframe using fortify.  Simplifies it down using data.frame and then reorders it according to the graphical position!
#Apparently fortify might deprecate and switch to the "broom" package for tidying data.  Try to figure this out on the ggtree object "q", not "tree", so that flip functions will be reflected in the output
tips <- fortify(tree)
tips <- data.frame(tips$label,tips$y,tips$isTip)
tips <- tips[order(tips[,3],tips[,2],decreasing=T),]
#Writes the tips to a csv file.  Name is based on the option -b specified when the script is called


dir.create(paste(getwd(),"/",opt$entry,"output", sep=''))
file_csv <- paste(getwd(),"/",opt$entry,"output/",opt$write,".csv", sep='')
message(file_csv)


node_count <- length(tree$tip.label)
for (i in 1:node_count) {
write(as.matrix(tips)[,1][i],file=file_csv,sep=",",append=T)
}

#Specifies tip shape using hmm, a column that comes from the attributes file
p <- p %<+% dd + geom_tiplab(size=opt$size,offset=0.05,aes(color=species)) + geom_tippoint(aes(size=opt$size,shape=hmm)) + scale_size_identity() #you need scale_size_identity! https://groups.google.com/forum/#!topic/bioc-ggtree/XHaq9Sk3b00

#Generates figure using ggtree, counts_file is a formatted list of numbers within the limits (e.g. RNAseq fold changes)
figure <- gheatmap(p,counts_file, offset = 1.5, width=0.6, font.size=1.5, colnames_angle=-45, hjust=0) + 
  #Color by species
	geom_tiplab(size=opt$size,offset=0.05,aes(color=species)) +
	scale_fill_gradient2(low="#000099",high="#FF0000",mid="white",limits=c(lower,upper)) +
	#node labels
	geom_text2(aes(subset=!isTip, label=node), hjust=-.3, size=opt$size) + 
	scale_colour_manual(values=c("black","red","blue","orange","purple","darkgreen","cadetblue","deeppink","darkgoldenrod","brown4","olivedrab2"))
figure
```

---

## pdf output
```{r eval = FALSE}
file_pdf <- paste(getwd(),"/",opt$entry,"output/",opt$write,".pdf", sep='')
message(file_pdf)

pdf(file_pdf)
figure
dev.off()
```

