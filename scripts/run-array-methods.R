library(data.table)


run.dmrcate = function(M, design, coef="genderM", outf){
    library(DMRcate)
    ilmn = readRDS('manifest.450.rds')
    rownames(M) = ilmn[rownames(M), "IlmnID"]
    myanno = cpg.annotate(M, analysis.type="differential", design=design, coef=coef, pcutoff=1)
    dmrs = dmrcate(myanno, lambda=300, pcutoff=1)
    df = dmrs$results
    pos = df$hg19coord
    chroms = unlist(lapply(strsplit(pos, ":", fixed=TRUE), function(r){ r[1] }))
    starts = unlist(lapply(strsplit(pos, ":", fixed=TRUE), function(r){
                     as.integer(unlist(strsplit(r[2], "-"))[[1]]) }))
    ends =   unlist(lapply(strsplit(pos, ":", fixed=TRUE), function(r){
                     as.integer(unlist(strsplit(r[2], "-"))[[2]]) }))

    write.table(data.frame(chrom=chroms, start=starts-1, end=ends, p=df$meanpval, df),
                quote=F, sep="\t", row.names=F, file=outf)
}


run.bumphunter = function(M, design, coef="genderM", outf){
    library(bumphunter)
    library(doParallel)

    chroms = unlist(lapply(strsplit(rownames(M), ":", fixed=TRUE), function(v) v[[1]]))
    posns = unlist(lapply(strsplit(rownames(M), ":", fixed=TRUE), function(v) as.integer(v[[2]])))

    coef = grep(coef, colnames(design), fixed=TRUE)
    if (!any(coef)){ stop(sprintf("%s not found in design", coef)) }
    stopifnot(length(coef) == 1)

    cl = clusterMaker(chroms, posns, maxGap=300)
    registerDoParallel(cores=3)
    tab = bumphunter(M, design, chroms, posns, cl, smooth=TRUE, B=200, cutoff=0.01, coef=coef)#pickCutoff=T)
    colnames(tab$table)[1] = "chrom"
    tab$table$p = tab$table$p.value # p.valueArea
    write.table(tab$table, file=outf, quote=F, sep="\t", row.names=F)
}

run.limma = function(M, design, coef="genderM", outf){
    library(limma)

    fit = eBayes(lmFit(M, design), trend=FALSE)
    top = topTable(fit, coef=coef, n=Inf, sort.by='none')

    chroms = unlist(lapply(strsplit(as.character(rownames(top)), ":", fixed=TRUE), function(r){ r[1] }))
    posns = unlist(lapply(strsplit(as.character(rownames(top)), ":", fixed=TRUE), function(r){ as.integer(r[2]) }))

    write.table(data.frame(chrom=chroms, start=posns - 1, end=posns, p=top$P.Value, top), quote=F, sep="\t", row.names=F,
                file=outf)
}

run.champ = function(M, design, coef, outf){
    library(ChAMP)
    M = 1 / (1 + exp(-M))
    design = data.frame(design)
    design$Sample_Group = design[,coef]
    design$Sample_Name = rownames(design)
    r = champ.lasso(beta.norm=M, pd=data.frame(design), filterXY=FALSE,
                    image=FALSE, mafPol.upper=0, minDmrSep=300,
                    minSigProbesLasso=1,  bedFile=TRUE,
                    resultsDir=dirname(outf),
                    adjPVal=0.2,
                    DMRpval=1)
    write.table(r, file=outf, quote=F, sep="\t", row.names=F)
}


run.betabinomial = function(M, design, coef, outf, mc.cores=9){
    library(betareg)
    library(parallel)
    model = paste(grep("Intercept", colnames(design), invert=T, value=T), collapse=" + ")
    model = sprintf("y ~ %s", model)
    df = data.frame(design)

    chroms = unlist(lapply(strsplit(as.character(rownames(M)), ":", fixed=TRUE), function(r){ r[1] }))
    posns = unlist(lapply(strsplit(as.character(rownames(M)), ":", fixed=TRUE), function(r){ as.integer(r[2]) }))

    M = 1 / (1 + exp(-M))

    res = mclapply(1:nrow(M), function(i){
        if(i %% 10000 == 0){ message(paste("at record", i)) }
        df$y = M[i,]
        r = summary(betareg(as.formula(model), df))$coefficients$mean[coef,]
        data.frame(chrom=chroms[i], start=posns[i] - 1,
                   end=posns[i], p=r[['Pr(>|z|)']])
  }, mc.cores=mc.cores)
  res = rbindlist(res)
  res$qvalue = p.adjust(res$p, "fdr")
  write.table(res, quote=F, row.names=F, sep="\t", file=outf)
}

args = commandArgs(TRUE)

covs = read.csv(args[1])
rownames(covs) = covs$id
M = as.matrix(read.delim(args[2], row.names=1))#, nrow=2000))
prefix = args[3]

design = model.matrix(~ gender + age, covs)
coef = ifelse(is.na(args[4]), "genderM", args[4])

run.dmrcate(M, design, coef, outf=sprintf("%s%s", prefix, "dmrcate.output.bed"))
run.limma(M, design, coef, outf=sprintf("%s%s", prefix, 'limma.output.bed'))
run.bumphunter(M, design, coef, outf=sprintf("%s%s", prefix, 'bumphunter.output.bed'))
run.betabinomial(M, design, coef, outf=sprintf("%s%s", prefix, 'betabinomial.output.bed'))
run.champ(M, design, coef, outf=sprintf("%s%s", prefix, 'champ.output.bed'))
