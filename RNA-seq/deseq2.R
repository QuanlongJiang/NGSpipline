library("DESeq2")
library("BiocParallel")
register(MulticoreParam(8))

path_in <- "../results/htseq/"
path_o <- '../results/deseq2/'
sampleFiles = grep('F' ,list.files(path_in),value=TRUE)
sampleCondition = sub('(.*day).*','\\1', sampleFiles )
sampleTable = data.frame(sampleName = sampleFiles, fileName = sampleFiles, condition = sampleCondition )

dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,directory = directory,design= ~ condition)
dds <- dds[ rowSums(counts(dds)) > 1, ] #

dds <- DESeq(dds)
#res <- results(dds,addMLE=TRUE)
day3 <- results(dds, contrast=c("condition","spd-_3day","spd+_3day"))
day10 <- results(dds, contrast=c("condition","spd_-_10day","spd_+_10day"))
 
#dds <- estimateSizeFactors(dds)
#counts(dds, normalized=TRUE)

write.csv(day3, paste(path_o,'degs_day3.csv',sep=''))
write.csv(day10,paste(path_o,'degs_day10.csv',sep=''))


