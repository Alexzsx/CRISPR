source("https://bioconductor.org/biocLite.R")
biocLite("devtools")
biocLite("BiocUpgrade")
biocLite("GenomeRanges")
library(devtools)
install_github("BIMSBbioinfo/genomation",build_vignettes=FALSE)

p300.url="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/wgEncodeAwgTfbsHaibH1hescP300V0416102UniPk.narrowPeak.gz"

Nanog.url="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/wgEncodeAwgTfbsHaibH1hescNanogsc33759V0416102UniPk.narrowPeak.gz"

SP1.url="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgTfbsUniform/wgEncodeAwgTfbsHaibH1hescSp1Pcr1xUniPk.narrowPeak.gz"

chrHMM.url="http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmH1hescHMM.bed.gz"

Segway.url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgSegmentation/wgEncodeAwgSegmentationSegwayH1hesc.bed.gz"
#http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgSegmentation/wgEncodeAwgSegmentationChromhmmH1hesc.bed.gz"

# get chromHMM annotation and make a list out of it
library("genomation")
chrHMM = readBed(chrHMM.url)
chrHMM.list = GenomicRanges::split(chrHMM, chrHMM$name,drop=TRUE)

Segway =readBed(Segway.url)
Segway.list = GenomicRanges::split(Segway, Segway$name,drop=TRUE)

# get peaks
p300=readBed(p300.url)
SP1=readBed(SP1.url)
NANOG=readBed(Nanog.url)

peak2ann=annotateWithFeatures(p300,chrHMM.list)
peak2ann

library("genomation")
library("GenomicRanges")
#chrList = paste('chr', FilteredData$CHROM, sep='')
#startList = FilteredData$POS
#endList = FilteredData$POS+1
#rangeList = IRanges(startList,endList)
#gr <- GRanges(chrList,rangeList)
data = read.csv('merge_Cho_c4bpb.csv', header = TRUE)
#data[,5]

gr <- GRanges(seqnames=data[,5], IRanges(start=data[,6], width=23))

snv2ann=annotateWithFeatures(gr,chrHMM.list)

snv2ann.segway = annotateWithFeatures(gr,Segway.list)

promoterMask = (snv2ann@members[,'1_Active_Promoter']==1)
enhancerMask = (snv2ann@members[,'4_Strong_Enhancer']==1)
enhancerMask2 = (snv2ann@members[,'5_Strong_Enhancer']==1)

data1 = data[promoterMask & enhancerMask & enhancerMask2,]

snv2ann.segway@members[,'EnhWf2']
Promoter = snv2ann@members[,'1_Active_Promoter']
Enhancer = snv2ann@members[,'4_Strong_Enhancer'] + snv2ann@members[,'5_Strong_Enhancer']
chrom = cbind(Promoter, Enhancer)
chromnew = cbind(data, chrom)

write.csv(chromnew,"merge_Cho_c4bpb_HMM.csv")




segwayPromoter = snv2ann.segway@members[,'Tss'] + snv2ann.segway@members[,'DnaseD'] 
segwayPromoter[which(segwayPromoter == 2)] = 1
segwayEnhancer =snv2ann.segway@members[,'EnhPr']

Seg = cbind(segwayPromoter, segwayEnhancer)
Segnew = cbind(data, Seg)

write.csv(Segnew,"merge_Cho_c4bpb_Seg.csv")

