source("http://bioconductor.org/biocLite.R")
biocLite("GenomicScores")
library(GenomicScores)
library(phastCons100way.UCSC.hg19)


data = read.csv('Crispr_Data.csv', header = TRUE)
data[,5]

numb = scores(phastCons100way.UCSC.hg19, GRanges(seqnames=data[,5], IRanges(start=data[,6], width=23)))

write.csv(numb,"phastCons.csv")


gsco <- getGScores("phyloP100way.UCSC.hg19")
numb2 = scores(gsco, GRanges(seqnames=data[,5], IRanges(start=data[,6], width=23)))

write.csv(numb2,"phyolop.csv")

