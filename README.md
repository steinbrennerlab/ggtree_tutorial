---
title: "ggtree tutorial"
output:
  github_document:
    toc: true
    toc_depth: 2
---
## About the scripts
Adam Steinbrenner <br>
astein10@uw.edu <br>
http://steinbrennerlab.org <br>
Updated 1/28/2020 <br>
<br>
The following describes ggtree.R, a custom script for using ggtree for phylogenetics tree visualization incorporating other attributes and heatmaps

---

## Download ggtree.R and data
Find ggtree.R in main directory, and all input files in /data/ggtree

---

## Install packages
The function tree_subset may require the development version of treeio; make sure it is current from github
```
source("https://bioconductor.org/biocLite.R")
biocLite("EBImage")
biocLite("treeio")
biocLite("ggtree")
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

##Usage and Options
Use options -e, the query ENTRY, and -o, the filename specifier OUTPUT. -n is optional to specify a specific node
e.g.
`Rscript tree.R -e Vigun07g219600.1 -o output_testing`

first two options specify the folder, blasted specifies the filename
```
option_list <- list( 
  make_option(c("-e", "--entry"), action="store", type="character",
        help="entry"),
	make_option(c("-o", "--output"), action="store", type="character", 
        help="name for file writing"),
	make_option(c("-n", "--node"), action="store", default=0, type="integer", 
        help="number of total sequences!  used to compute font size")
	)
opt <- parse_args(OptionParser(option_list=option_list))
```

---

##Inputs: newick tree, attribute files, heatmap data
Newick tree
```{r}
dir<-paste("C:/Users/Adam/Dropbox/github/alluvial_diagrams/alluvial_diagrams/data/ggtree/")
message(dir)
#tree <- paste(dir,opt$entry,sep="")
tree <- read.tree(paste(dir,"tree_input.nwk",sep=""))
```
Attribute files, tab delimited with columns taxa and atttribute
Species, specifying from which genome each gene entry came from
```{r}
dd <- read.table("C:/Users/Adam/Dropbox/github/alluvial_diagrams/alluvial_diagrams/data/ggtree/attribute_species.txt", sep="\t", header = TRUE, stringsAsFactor=F)
```
HMM, an output of significant hits from an HMMer pipeline
```{r}
dd2 <- read.table("C:/Users/Adam/Dropbox/github/alluvial_diagrams/alluvial_diagrams/data/ggtree/attribute_hmm.txt", sep="\t", header = TRUE, stringsAsFactor=F)
```
Heatmap, e.g. a list of log2 Fold Changes associated with the genes
```{r}
counts_file <- read.table("C:/Users/Adam/Dropbox/github/alluvial_diagrams/alluvial_diagrams/data/ggtree/log2FoldChanges.txt", sep="\t", row.names = 1, header = TRUE, stringsAsFactor=F)
```
Set upper and lower limits for the heatmap
```{r}
upper <- 5
lower <- -5
```

---

##ggTree object
ggTree creates a tree visualization using ggplot and feeds in the various dataframes (dd, dd2).  You can call columns from these inputs to specify different aspects of the visualization
```{r}
q <- ggtree(tree, size=0.1) #size specifies line size
###OPTIONAL: adds hmm coding to ggtree object q
###OPTIONAL: takes hmm_coding and adds to a separate dataframe dd2

q <- q %<+% dd2

#number of nodes
node_count <- length(tree$tip.label)

#sets tip and node label size
size <-  3.63 - (0.484*log(node_count)) #computes appropriate font size for tree based on good sizes for 10, 30, and 100
#size <- 2.5 #use a default font size instead
size2 <- (size/2)

#Specifies tip shape using hmm, a column that comes from the attributes file
q <- q %<+% dd + geom_tiplab(size=size,offset=0.05,aes(color=species)) + geom_tippoint(aes(size=size2,shape=hmm)) + scale_size_identity() #you need scale_size_identity! https://groups.google.com/forum/#!topic/bioc-ggtree/XHaq9Sk3b00

#Generates figure using ggtree, counts_file is a formatted list of numbers within the limits (e.g. RNAseq fold changes)
figure <- gheatmap(q,counts_file, offset = 1.5, width=0.6, font.size=1.5, colnames_angle=-45, hjust=0) + 
  #Color by species
	geom_tiplab(size=size,offset=0.05,aes(color=species)) +
	scale_fill_gradient2(low="#000099",high="#FF0000",mid="white",limits=c(lower,upper)) +
	#node labels
	geom_text2(aes(subset=!isTip, label=node), hjust=-.3, size=size2) + 
	scale_colour_manual(values=c("black","red","blue","orange","purple","darkgreen","cadetblue","deeppink","darkgoldenrod","brown4","olivedrab2"))
figure
```

---

##pdf output
```{r eval = FALSE}
file <- "C:/Users/Adam/Dropbox/github/alluvial_diagrams/alluvial_diagrams/data/ggtree/ggtree_output.pdf"
message(file)

pdf(file)
figure
dev.off()
```

---

##Fortify to write list of genes
Takes the tree object and converts it to a tidy dataframe using fortify, then reorders it according to the graphical position
Apparently fortify might deprecate and switch to the "broom" package for tidying data
```{r eval = FALSE}
tips <- fortify(tree)
tips <- data.frame(tips$label,tips$y,tips$isTip)
tips <- tips[order(tips[,3],tips[,2],decreasing=T),]
#Writes the tips to a csv file.  Name is based on the option -b specified when the script is called
file_csv <- "C:/Users/Adam/Dropbox/github/alluvial_diagrams/alluvial_diagrams/data/ggtree/ggtree_output.csv"
message(file_csv)

for(i in 1:node_count) {
write(as.matrix(tips)[,1][i],file=file_csv,sep=",",append=T)
}
```
Using the output csv, you can call a python script to get the original nucleotide sequences from the parsed, merged, fasta file!  This output fa is in the same order as the tree
```
system(paste("python extract_seq.py ",opt$entry," ",opt$output,".csv",sep=""))
```

##Subtrees
OPTIONAL: can take a subset of the original tree; use a node 1 deeper than you want!  For example, try the NIK3 clade using `node <- 34`
```{r}
#node<-opt$node
node<-34
###
if (node>0) {
#nodenum <- opt$node
nodenum <- node
tree <- tree_subset(tree, nodenum, levels_back = 1) #Right now a bug with levels_back=0 is preventing me from specifying the node ITSELF
}

q <- ggtree(tree, size=0.1) #size specifies line size
###OPTIONAL: adds hmm coding to ggtree object q
###OPTIONAL: takes hmm_coding and adds to a separate dataframe dd2
dd2 <- read.table("C:/Users/Adam/Dropbox/github/alluvial_diagrams/alluvial_diagrams/data/ggtree/attribute_hmm.txt", sep="\t", header = TRUE, stringsAsFactor=F)
q <- q %<+% dd2

#number of nodes
node_count <- length(tree$tip.label)

#sets tip and node label size
size <-  3.63 - (0.484*log(node_count)) #computes appropriate font size for tree based on good sizes for 10, 30, and 100
#size <- 2.5 #use a default font size instead
size2 <- (size/2)

q <- q %<+% dd + geom_tiplab(size=size,offset=0.05,aes(color=species)) + geom_tippoint(aes(size=size2,shape=hmm)) + scale_size_identity() #you need scale_size_identity! https://groups.google.com/forum/#!topic/bioc-ggtree/XHaq9Sk3b00

figure <- gheatmap(q,counts_file, offset = 1.5, width=0.6, font.size=1.5, colnames_angle=-45, hjust=0) + 
	geom_tiplab(size=size,offset=0.05,aes(color=species)) +
	scale_fill_gradient2(low="#000099",high="#FF0000",mid="white",limits=c(lower,upper)) +
	#node labels
	geom_text2(aes(subset=!isTip, label=node), hjust=-.3, size=size2) + 
	scale_colour_manual(values=c("black","red","blue","orange","purple","darkgreen","cadetblue","deeppink","darkgoldenrod","brown4","olivedrab2"))
figure

```

Old code is included in the R markdown file below but hidden in the html rendering.  It might contain some useful snippets
```{r echo = FALSE, eval = FALSE}
#Below is idea for iterating it over every ndode
#for (i in 1:node_count) {
#	subtree <- tree_subset(tree,i)
#}
#b <- as.data.frame(table(p$data$hit))
#b
#8/22: get this to only append the last two columns, include the node number, and do it horizontally!
#c <- as.data.frame(table(p$data$hit),responseName = "node")
#c
#d <- merge(b,c,by="Var1",all=TRUE) # all=TRUE fills missing rows if the dataframes are unequal
#d[is.na(d)] <- 0 # replaces n/a with 0
#d
#write.table(d,file=file_csv,sep=",",append=T)


#Idea, for each node output the count of each type of column "hit".  
##The ggtree object has a dataframe, data, with hit added to it.  A subtree should have all these features
```