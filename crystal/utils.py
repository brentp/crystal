from collections import defaultdict

import toolshed as ts
import numpy as np
import pandas as pd
from .crystal import Feature, CountFeature
import os

def example_random_cluster(n_samples, n_sites, seed=42):
    np.random.seed(seed)
    if n_samples % 2 != 0: n_samples += 1

    covs = pd.DataFrame({'gender': ['F'] * (n_samples / 2) + ['M'] * (n_samples / 2),
                         'age': np.random.uniform(10, 25, size=n_samples) })

    methylation = np.random.normal(0.15, 0.75, size=(n_sites, n_samples))

    cluster = [Feature('chr1', (i + 1) * 10, m) for i, m in enumerate(methylation)]
    covs['id'] = ['id_%i' %i for i in range(len(covs))]

    return covs, cluster

def real_cluster():

    path = os.path.join(os.path.dirname(__file__), "tests")

    meth = pd.read_csv('%s/real_cluster.csv' % path, index_col=0)
    chroms = [x.split(":")[0] for x in meth.index]
    starts = [int(x.split(":")[1]) for x in meth.index]

    cluster = [Feature(chroms[i], starts[i], np.array(meth.ix[i, :])) for i in
            range(len(meth))]

    covs = pd.read_csv('%s/covs.csv' % path)
    return covs, cluster

def real_count_cluster():

    path = os.path.join(os.path.dirname(__file__), "tests")

    c = pd.read_csv('%s/m.counts.csv' % path, index_col=0)
    m = pd.read_csv('%s/m.methylated.csv' % path, index_col=0)

    chroms = [x.split(":")[0] for x in m.index]
    starts = [int(x.split(":")[1]) for x in c.index]

    cluster = [CountFeature(chroms[i], starts[i],
                            np.array(m.ix[i, :]),
                            np.array(c.ix[i, :]))
                  for i in range(len(m))]

    covs = pd.read_table('%s/m.covs.txt' % path)
    return covs, cluster


def write_cluster(cluster, fh, float_format="%.4f", count_fh=None):
    """
    Write a cluster to file.

    Parameters
    ----------

    cluster : cluster
              a cluster from aclust (or just a list of features)

    fh : filehandle

    count_fh : filehandle
               if cluster is of `CountFeature` then a count_fh must be
               specified so that the counts can be written to file as
               well.

    """

    fmt = "{chrom}:{position}\t{values}\n"
    if isinstance(cluster[0], Feature):
        for f in cluster:
            fh.write(fmt.format(chrom=f.chrom, position=f.position,
                    values="\t".join((float_format % v for v in f.values))))
        return


    elif isinstance(cluster[0], CountFeature):
        assert count_fh is not None
        for f in cluster:
            fh.write(fmt.format(chrom=f.chrom, position=f.position,
                    values="\t".join(f.methylated)))
            count_fh.write(fmt.format(chrom=f.chrom, position=f.position,
                    values="\t".join(f.counts)))

def roc_out(p_bed, p_col, truth_region_bed, exclude=('-1', 'NA', 'nan')):
    """Create ROC for a bed file of p-values given known truth regions.

    Parameters
    ----------

    p_bed : file

    p_col : int
            column containing the p-value from `p_bed`

    truth_region_bed : file
                       contains the true regions
    """
    p_col -= 1 # 0-based

    regions = defaultdict(list)
    for toks in ts.reader(truth_region_bed, header=False):
        if not (toks[1] + toks[2]).isdigit(): continue
        regions[toks[0]].append((int(toks[1]), int(toks[2])))

    truths = []
    vals = []
    for toks in ts.reader(p_bed, header=False):
        if not (toks[1] + toks[2]).isdigit(): continue
        reg = regions[toks[0]]

        s, e = int(toks[1]), int(toks[2])

        p = toks[p_col]
        if p in exclude: continue
        vals.append(1.0 - float(p))

        truth = any(rs <= s <= re or rs <= e <= re for rs, re in reg)
        truths.append(truth)

    return np.array(truths).astype(int), np.array(vals)
