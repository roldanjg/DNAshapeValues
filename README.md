# DNAshapeValues
Pipeline to get the values of Roll, MGW, ProT and HelT from a genome.fa file
#### **Previous steps in R**

'#if (!requireNamespace("BiocManager", quietly = TRUE))
'#  install.packages("BiocManager")

'# The following initializes usage of Bioc devel
'# BiocManager::install(version='devel')

'# BiocManager::install("DNAshapeR")
setwd('C:/Users/JOAQUINGR/PycharmProjects/untitled/tomato_genome/')
library(DNAshapeR)
library(Biostrings)

pred3 <- getShape("tomato_SL4.0ch12")


fn3 <- system.file("extdata", "cortito.fa", package = "DNAshapeR")
fn3
write.table(pred3[["HelT"]], file = "Sol.Helt.txt",sep = "\t",row.names=TRUE, col.names=TRUE)
write.table(pred3[["Roll"]], file = "Sol.Roll.txt",sep = "\t",row.names=TRUE, col.names=TRUE)
write.table(pred3[["ProT"]], file = "Sol.ProT.txt",sep = "\t",row.names=TRUE, col.names=TRUE)
write.table(pred3[["MGW"]], file = "Sol.MGW.txt",sep = "\t",row.names=TRUE, col.names=TRUE)

head(pred3)
featureType <- c("1-mer", "1-shape")
featureVector <- encodeSeqShape(fn3, pred3, featureType)
head(featureVector)

head(fn3)
<br>
You might have to do some cleaning to get the shape values, if genome is too large, split the genome firs,
then use the R script to get the values and the getting_shapevals.py to get the wig files that you will use to create BigWig 
using the joinfiles.sh command 

