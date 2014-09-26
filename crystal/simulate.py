import sys
import toolshed as ts
import numpy as np
import scipy.stats as ss
from copy import deepcopy
from collections import defaultdict
from . import CountFeature
from statsmodels.genmod.families import Poisson

SIZES = dict.fromkeys([1, 2] + range(4, 12, 2), 100)
SIZES[2] = SIZES[1] = 200

def rr_cluster(cluster, covs, formula="methylation ~ age"):
    """Set cluster values to reduced-residuals."""
    cluster = deepcopy(cluster)
    from statsmodels.formula.api import ols, glm

    if isinstance(cluster[0], CountFeature):
        for f in cluster:
            covs['methylation'] = f.methylated
            f.methylated[:] = np.round(glm(formula,
                                           covs,
                                           exposure=f.counts,
                                           family=Poisson()
                                          ).fit().resid
                                       ).astype(int)
            f.values[:] = f.methylated.astype(float) / f.counts
    else:
        for f in cluster:
            covs['methylation'] = f.values
            f.values[:] = ols(formula, covs).fit().resid
    return cluster


choice = np.random.choice

def simulate_cluster(cluster, w=0, class_order=None, attrs=(),
                     get_reduced_residuals=None, grr_args=()):
    """Modify the data in existing clusters to create or remove and effect.

    Parameters
    ----------

    cluster : list of clusters
              *should* include clusters of length 1.

    w : float
        w = 0 generates random data. Higher values are
        more likely to separate the groups.

    class_order : np.array of 0/1
        list of same length as cluster[*].values indicating
        which group (0 or 1) each sample belongs to.

    get_reduced_residuals : function
        optional, see :func:`~simulate_regions`

    """
    if isinstance(cluster[0], CountFeature) and attrs==():
        attrs = ('methylated', 'counts')

    # copy since we modify the data in-place.
    if get_reduced_residuals is None:
        cluster = [deepcopy(f) for f in cluster]
    else:
        cluster = get_reduced_residuals(cluster, *grr_args)

    # make this one faster since we don't need to keep track of high/low
    if w == 0:
        idxs = np.arange(len(cluster[0].values))
        np.random.shuffle(idxs)
        for feature in cluster:
            feature.values = feature.values[idxs]
            for attr in attrs:
                setattr(feature, attr, getattr(feature, attr)[idxs])
        return cluster


    N = len(cluster[0].values)
    assert N % 2 == 0
    n = N / 2

    if class_order is None:
        class_order = np.zeros_like(cluster[0].values).astype(int)
        class_order[n:] = 1

    n_probes = len(cluster)
    new_data = np.zeros((n_probes, 2 * n))

    # choose a random probe from the set.
    # the values across the cluster will be determined
    # by the values in this randomly chose probe.
    i = choice(range(n_probes))#[0] if n_probes > 1 else 0
    c = cluster[i]

    idxs = np.arange(N)

    # just pull based on the index. so we need to sort the values
    # as well.
    idx_order = np.argsort(c.values)

    ords = np.arange(1, N + 1) / (N + 1.0)
    ords = (1.0 - ords)**w
    h_idxs = choice(idxs, replace=False, p=ords/ords.sum(), size=n)

    idxs = np.setdiff1d(idxs, h_idxs, assume_unique=True)
    idxs.sort()

    ords = np.arange(1, N + 1 - n) / (N + 1.0 - n)
    assert ords.shape == idxs.shape
    ords = (ords)**w
    l_idxs = choice(idxs, replace=False, p=ords/ords.sum(), size=n)

    assert len(np.intersect1d(h_idxs, l_idxs)) == 0
    for j in range(n_probes):
        tmph = np.array(cluster[j].values[idx_order][h_idxs])
        tmpl = np.array(cluster[j].values[idx_order][l_idxs])
        cluster[j].values[class_order == 0] = tmph
        cluster[j].values[class_order == 1] = tmpl
        for attr in attrs:
            vals = getattr(cluster[j], attr)[idx_order]
            hi_vals = np.array(vals[idx_order][h_idxs])
            lo_vals = np.array(vals[idx_order][l_idxs])
            vals[class_order == 0] = hi_vals
            vals[class_order == 1] = low_vals

    return cluster


def simulate_regions(clust_list, region_fh, sizes=SIZES, class_order=None,
        seed=42, get_reduced_residuals=None, get_reduced_residuals_args=()):
    """Simulate regions and randomize others.

    Parameters
    ----------

    clust_list : list of clusters
        should include clusters of length 1.

    region_fh : filehandle
        a BED file of all position will be written to this
        file. The 4th column will indicate true/false indicating
        if it was simulated to have a difference. The fifth
        column will indicate the size of the cluster it was in.

    size : dict
        keys of the clust_size and values of how many clusters
        to create of that size. Default is to create 100 of each
        size from 3 to 8 and 200 clusters of size one and 2. All
        others are randomized.

    classes : np.array
        same length as cluster[*].values indicating
        which group each sample belongs to.

    seed: int

    get_reduced_residuals : function
        If this parameter is None, then they values are shuffled as they
        are received.
        A function that accepts a cluster and returns residuals of the reduced
        model. e.g. if the full model of interest is:
            methylation ~ disease + age + gender
        the reduced model would be:
            methylation ~ age + gender
        so that only the residuals of the reduced model are shuffled and the
        other effects should remain. This will implement the bootsrap for
        linear models from Efron and Tibshirani.
        An Example function would be: :func:`~rr_cluster`


    Returns
    -------

    generator of clusters in the same order as clust_list.

    """
    np.random.seed(seed)
    assert isinstance(clust_list, list), ("need a list due to multiple \
            iterations")

    if class_order is not None:
        class_order = np.array(class_order)
        classes = np.unique(class_order)
        assert len(classes) == 2, (classes, "should have 2 unique")
        classes = {classes[0]: 0, classes[1]: 1}
        class_order = np.array([classes[c] for c in class_order])

    clusts = defaultdict(list)
    for clust in clust_list:
        clusts[len(clust)].append(clust)
    clusts = dict(clusts)

    sim_idxs = {}
    # for each size of clust, choose n random indices based on how
    # many of that cluster we saw.
    for size, n in sizes.items():
        idxs = np.arange(len(clusts[size]))
        # get the indexes of the clusters we want
        sim_idxs[size] = frozenset(np.random.choice(idxs, size=min(n,
                                                    len(idxs)), replace=False))

    fmt = "{chrom}\t{start}\t{end}\t{truth}\t{size}\n"
    region_fh.write(ts.fmt2header(fmt))

    seen = defaultdict(int)
    for c in clust_list:
        l = len(c)
        w = 0
        # need this if block in case we get a cluster longer
        # than we have in sizes
        if l in sim_idxs:
            s = seen[l]
            w = int(s in sim_idxs[l])
            seen[l] += 1

        truth = l in sim_idxs and s in sim_idxs[l]
        for f in c:
            region_fh.write(fmt.format(chrom=f.chrom,
                                       start=f.position - 1,
                                       end=f.position,
                                       truth="true" if truth else "false",
                                       size=len(c)))
        yield simulate_cluster(c, w, class_order, get_reduced_residuals)

    region_fh.flush()

