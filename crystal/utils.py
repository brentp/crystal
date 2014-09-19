import numpy as np
import pandas as pd
from .crystal import Feature
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
