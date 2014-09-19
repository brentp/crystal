import numpy as np
import pandas as pd
from .crystal import Feature

def example_random_cluster(n_samples, n_sites, seed=42):
    np.random.seed(seed)
    if n_samples % 2 != 0: n_samples += 1
     
    covs = pd.DataFrame({'gender': ['F'] * (n_samples / 2) + ['M'] * (n_samples / 2),
                         'age': np.random.uniform(10, 25, size=n_samples) })

    methylation = np.random.normal(size=(n_sites, n_samples))

    cluster = [Feature('chr1', i* 10, m) for i, m in enumerate(methylation)]

    return covs, cluster
