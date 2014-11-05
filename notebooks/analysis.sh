### Simulation
Rscript scripts/run-array-methods.R notebooks/v0-450K-12.covariates.csv notebooks/ex-sim-meth.txt work/
comb-p pipeline -c p --dist 300 --step 100 -p work/comb-p/cpv-limma work/limma.output.txt

### Replication
Rscript scripts/run-array-methods.R notebooks/d0-450K-12.covariates.csv notebooks/d0-450K-12.meth.txt.gz work/discovery age
Rscript scripts/run-array-methods.R notebooks/v0-450K-12.covariates.csv notebooks/v0-450K-12.meth.txt.gz work/validation age
for cohort in discovery validation; do
    comb-p pipeline -c p --dist 300 --step 100 -p work/comb-p-replication/$cohort work/${cohort}limma.output.bed
done


