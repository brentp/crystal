suppressPackageStartupMessages(library(BiSeq))
suppressPackageStartupMessages(library(data.table))
library(data.table)

args = commandArgs(TRUE)

cov = read.delim(args[1])
print(head(cov))
rownames(cov) = cov$Sample

dcounts = read.delim(args[2], row.names=1, quote='')
dM = read.delim(args[3], row.names=1, quote='')


chroms = unlist(lapply(strsplit(rownames(dM), ":", fixed=TRUE), function(r) r[[1]]))
posns = unlist(lapply(strsplit(rownames(dM), ":", fixed=TRUE), function(r) as.integer(r[[2]])))

M = as.matrix(dM)
counts = as.matrix(dcounts)
print(M[1:5, 1:5])
print(counts[1:5, 1:5])

raw = BSraw(colData=DataFrame(cov),
            rowData=GRanges(chroms, IRanges(start=posns, width=1)),
            totalReads=counts,
            methReads=M)

print(colData(raw)$Group)
raw.clust <- clusterSites(object = raw,
                                groups = colData(raw)$Group,
                                perc.samples = 4/5,
                                min.sites = 6,
                                max.dist = 200)

ind.cov <- totalReads(raw.clust) > 0
quant <- quantile(totalReads(raw.clust)[ind.cov], 0.9)

raw.lim = limitCov(raw.clust, maxCov=quant)
predictedMeth = predictMeth(object=raw.lim)
br = betaRegression(formula = ~ Group, link="probit", object=predictedMeth, type="BR", mc.cores=4)


s = sample(cov$Group)
print(s)
colData(predictedMeth)$group.null = s

br.null = betaRegression(formula = ~group.null,
                link = "probit",
                object = predictedMeth,
                type="BR", mc.cores=4)

vario <- makeVariogram(br.null)
vario.sm <- smoothVariogram(vario, sill = 0.9)
print("done vario null")

vario.aux <- makeVariogram(br, make.variogram=FALSE)
print("done makeVariogram(br)")
vario.sm$pValsList = vario.aux$pValsList
print("done vario")

locCor <- estLocCor(vario.sm)
print("done estLocCor")
clusters.rej <- testClusters(locCor, FDR.cluster = 0.9999)
print("done testClusters")
print(dim(clusters.rej))
clusters.trimmed <- trimClusters(clusters.rej, FDR.loc = 0.9999)
print("done trim Clusters")
DMRs = findDMRs(clusters.trimmed, max.dist=100, diff.dir=TRUE)

write.table(as.data.frame(DMRs), file='biseq.res.txt', row.names=F, sep="\t", quote=F)

