from __future__ import print_function
import sys
import toolshed as ts
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from aclust import mclust
import crystal

covariates_file = '../../../crystal/tests/covs.csv'
methylation_file = '../../../crystal/tests/meth.txt.gz'
formula = 'methylation ~ age + gender'
coef = 'gender'

covs = pd.read_csv(covariates_file)

def feature_gen(fname):
    for i, d in enumerate(ts.reader(fname, header=False)):
        if i == 0: continue
        chrom, pos = d[0].split(":")
        yield crystal.Feature(chrom, int(pos), crystal.logit(np.array(map(float, d[1:]))))


cluster_iter = mclust(feature_gen(methylation_file), max_dist=100, max_skip=0)


fmt = "{chrom}\t{start}\t{end}\t{p:.4g}\t{coef:.3f}\t{n_sites:d}"
print(ts.fmt2header(fmt))
for i, c in enumerate(crystal.model_clusters(cluster_iter,
                                             covs, formula, coef,
                                            model_fn=crystal.zscore_cluster,
                                            n_cpu=1)):
    print(fmt.format(**c))
    if c['p'] < 1e-3 and abs(c['coef']) > 0.2 and c['n_sites'] > 3:
        crystal.plot.spaghetti_plot(c, covs)
        plt.savefig('/tmp/figure-1.eps')
        break
