options(stringsAsFactors=FALSE)


targets = read.csv('notebooks/d0-450K-12.covariates.csv', header=T, row.names=1)
for(f in c('notebooks/d1-450K-12.covariates.csv',
           'notebooks/v0-450K-12.covariates.csv',
           'notebooks/v1-450K-12.covariates.csv')){
    targets = rbind(targets, read.csv(f, header=T, row.names=1))
}

targets$id = rownames(targets)
targets$Basename = paste(targets$id, targets$Basename, sep="_")
#write(paste("/drive/450k/", paste0(targets$Basename, '_Grn.idat.gz', sep=""), sep=""), stdout())
#write(paste("/drive/450k/", paste0(targets$Basename, '_Red.idat.gz', sep=""), sep=""), stdout())

library(minfi)

files = paste0("/drive/450k/", targets$Basename, "_Red.idat")

targets$Basename = paste0("/drive/450k/", targets$Basename)

rgSet = read.450k.exp(base="/drive/450k/", targets) #, extended=TRUE)
failed = detectionP(rgSet) > 0.01
bad_probes = rowMeans(failed) > 0.10
#rgSet = rgSet[rowMeans(failed) < 0.1,]
#rm(bad_probes, failed); gc()
xprobes = read.delim('notebooks/non-specific.txt', header=T)[,1]
#rgSet = rgSet[!rownames(rgSet) %in% xprobes,]

methods = c('raw', 'quantile', 'swan', 'funnorm')
i = 1
for(method in c(preprocessRaw, preprocessQuantile, preprocessSWAN, preprocessFunnorm)){
    print(methods[i])
    m = method(rgSet)
    m = m[!bad_probes,]
    m = m[!rownames(m) %in% xprobes,]

    b.norm = getBeta(m)
    # squeeze to prevent 0, 1 according to Smithson (lemon squeezer)
    b.norm = ((b.norm - 0.5) * 0.99) + 0.5
    print("OK")
    locs = getLocations(m)
    rownames(b.norm) = paste(locs@seqnames, locs@ranges@start, sep=":")
    colnames(b.norm) = targets$id
    snames = as.character(locs@seqnames)
    starts = as.integer(locs@ranges@start)
    b.norm = b.norm[order(snames, starts),]
    write.table(data.frame(probe=rownames(b.norm), format(b.norm, digits=4)),
                file=sprintf("norm.%s.txt", methods[i]), sep="\t", quote=F, row.names=F)

    i = i + 1
}
