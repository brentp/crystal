import sys
import toolshed as ts
import numpy as np
import scipy.stats as ss
from copy import deepcopy
from collections import defaultdict
from . import CountFeature
from statsmodels.genmod.families import Poisson

SIZES = {}
SIZES[3] = 600
SIZES[5] = 200
SIZES[8] = 150
SIZES[9] = 100
SIZES[10] = 50

def scale(values):
    values = np.asarray(values).astype(float)
    values += values.min()
    values /= values.max()
    # squeeze
    values = 0.5 + 0.90 * (values - 0.5)
    return values / values.sum()

def rr_cluster(cluster, covs, formula):
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
            fit = ols(formula, covs).fit()
            f.values[:] = fit.resid
            f.ovalues = fit.fittedvalues
    return cluster


choice = np.random.choice

def simulate_cluster(cluster, w=0, class_order=None,
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
    if isinstance(cluster[0], CountFeature):
        attrs = ('methylated', 'counts')
    else:
        attrs = []

    # copy since we modify the data in-place.
    has_fitted = False
    if get_reduced_residuals is None:
        cluster = [deepcopy(f) for f in cluster]
    else:
        has_fitted = True
        cluster = get_reduced_residuals(cluster, *grr_args)

    # make this one faster since we don't need to keep track of high/low
    if w == 0:
        idxs = np.arange(len(cluster[0].values))
        np.random.shuffle(idxs)
        for feature in cluster:
            feature.values = feature.values[idxs]
            for attr in attrs:
                setattr(feature, attr, getattr(feature, attr)[idxs])
            # add the shuffled residuals to the fitted values
            if has_fitted:
                feature.values += feature.ovalues
        return cluster


    N = len(cluster[0].values)
    nH, r = divmod(N, 2)
    nL = nH + 1
    if nL + nH > N:
        nL -= 1

    if class_order is None:
        class_order = np.zeros_like(cluster[0].values).astype(int)
        class_order[nH:] = 1

    else:
        nL = (class_order == 0).sum()
        nH = (class_order == 1).sum()

    assert nL + nH == N, (nL, nH, N)
    assert nL * nH > 0

    n_probes = len(cluster)

    # choose a random probe from the set.
    # the values across the cluster will be determined
    # by the values in this randomly chose probe.
    #i = choice(range(n_probes))
    #c = cluster[i]

    # choose the probe with the highest variance from the set.
    i = np.argmax([np.var(f.values) for f in cluster])
    c = cluster[i]

    idxs = np.arange(N)

    # just pull based on the index. so we need to sort the values
    # as well.
    idx_order = np.argsort(c.values)
    p_values = scale(c.values)

    # HI : TODO: instead of using ords, use weights of the actual
    # methylation values
    #ords = np.arange(1, N + 1) / (N + 1.0)
    #ords = (1.0 - ords)**w
    # sample with replacement to allow for stronger differences.
    #h_idxs = choice(idxs, replace=True, p=ords/ords.sum(), size=nH)
    #h_idxs.sort()
    h_idxs = choice(idxs, replace=True, p=p_values, size=nH)
    l_idxs = choice(idxs, replace=True, p=1 - p_values, size=nL)


    """
    # LO -without replacement
    l_idxs = np.setdiff1d(idxs, h_idxs, assume_unique=True)
    l_idxs.sort()

    if len(l_idxs) + len(h_idxs) != N:
        # only need to do this if there's a choice.
        l_ords = np.arange(1, nL + 1.) / (nL + 1.0)
        assert l_ords.shape == l_idxs.shape
        l_ords = (l_ords)**w
        l_idxs = choice(l_idxs, replace=False, p=l_ords/l_ords.sum(), size=nL)
    assert len(np.intersect1d(h_idxs, l_idxs)) == 0
    """

    #l_ords = np.arange(1, N + 1) / (N + 1.0)
    #l_ords = l_ords ** w
    #l_idxs = choice(idxs, replace=True, p=l_ords/l_ords.sum(), size=nL)

    for j in range(n_probes):
        #tmpl = np.array(cluster[j].values[idx_order][l_idxs])
        #tmph = np.array(cluster[j].values[idx_order][h_idxs])
        tmpl = np.array(cluster[j].values[l_idxs])
        tmph = np.array(cluster[j].values[h_idxs])
        cluster[j].values[class_order == 0] = tmpl
        cluster[j].values[class_order == 1] = tmph
        # add the shuffled residuals to the fitted values
        # actually add the fitted values back on the the shuffled residuals
        if has_fitted:
            cluster[j].values += cluster[j].ovalues

        for attr in attrs:
            vals = getattr(cluster[j], attr)#[idx_order]
            lo_vals = np.array(vals[l_idxs])
            hi_vals = np.array(vals[h_idxs])
            vals[class_order == 0] = lo_vals
            vals[class_order == 1] = hi_vals
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
        same length as cluster[i].values indicating
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
    from math import log
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
    del clusts
    fmt = "{chrom}\t{start}\t{end}\t{truth}\t{size}\n"
    region_fh.write(ts.fmt2header(fmt))

    seen = defaultdict(int)
    changed = defaultdict(int)
    for c in clust_list:
        l = len(c)
        w = 0
        # need this if block in case we get a cluster longer
        # than we have in sizes
        if l in sim_idxs:
            s = seen[l]
            # denominator sets larger DMRs to have a smaller per-probe effect.
            #w = 2 * int(s in sim_idxs[l])
            #w = int(s in sim_idxs[l]) / (log(l + 1))
            #w = 2.0 * int(s in sim_idxs[l]) / (1 + log(l, 2))
            #w = int(s in sim_idxs[l]) * 2 ** (1.0 / l)
            #w = (3. * int(s in sim_idxs[l])) ** (1.0 / l)
            w = 2 * int(s in sim_idxs[l]) / log(l * 2)

            seen[l] += 1
            if w > 0: changed[l] += 1

        truth = l in sim_idxs and s in sim_idxs[l]
        for f in c:
            region_fh.write(fmt.format(chrom=f.chrom,
                                       start=f.position - 1,
                                       end=f.position,
                                       truth="true" if truth else "false",
                                       size=len(c)))
        yield simulate_cluster(c, w, class_order, get_reduced_residuals,
                get_reduced_residuals_args)
    sys.stderr.write("changed:" + str(dict(changed)) + "\n")
    sys.stderr.write("total:" + str(dict(seen)) + "\n")
    region_fh.flush()

