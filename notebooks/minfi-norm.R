options(stringsAsFactors=FALSE)

targets = read.csv('covariates.csv', header=T, row.names=1)
targets$Basename = paste(rownames(targets), targets$Basename, sep="_")
print(head(targets))
targets$id = rownames(targets)

library(minfi)


rgSet = read.450k.exp(base=".", targets, extended=TRUE)
#qcReport(rgSet)

failed = detectionP(rgSet) > 0.01
bad_probes = rowMeans(failed) > 0.10
m = preprocessQuantile(rgSet)[!bad_probes,]
rm(rgSet, failed, bad_probes); gc()
message('B')


b.norm = getBeta(m)
locs = getLocations(m)
rm(m); gc()
message('C')
mdsPlot(b.norm)

message(sprintf("seqnames length: %d, dim: %d, %d", length(locs@seqnames),
                        nrow(locs@ranges@start), ncol(locs@ranges@start)))

xprobes = read.delim('notebooks/non-specific.txt', header=T)[,1]

b.norm = b.norm[!rownames(b.norm) %in% xprobes,]
    
rownames(b.norm) = paste(locs@seqnames, locs@ranges@start, sep=":")
colnames(b.norm) = targets$id
snames = as.character(locs@seqnames)
starts = as.integer(locs@ranges@start)
b.norm = b.norm[order(snames, starts),]
write.table(data.frame(probe=rownames(b.norm), format(b.norm, digits=4)), file="norm.beta.txt", sep="\t", quote=F, row.names=F)


