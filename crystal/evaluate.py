from collections import defaultdict
from .crystal import model_clusters
import toolshed as ts
import itertools as it
import seaborn as sns
import pandas as pd
import numpy as np
import sys



def evaluate_replication(discovery_clusters, replication_feats,
                         discovery_covs, replication_covs, formula,
                         coef, model_fn, pool=None, kwargs=None):
    """
    Evaluate a method (`model_fn`) by reporting p-values for discovery and
    replication cohorts.

    Parameters
    ----------

    discovery_clusters : list
        list of clusters from `aclust` from the discovery cohort.

    replication_feats : list
        list of features from the replication cohort. Clusters are created by
        mirroring those found in the discovery cohort.

    discovery_covs : DataFrame
        pandas.DataFrame of covariates for the discovery cohort. Rows must be
        in same order as the feature values in the discovery_cluster.

    replication_covs : DataFrame
        pandas.DataFrame of covariates for the repliation cohort. Rows must be
        in same order as the feature values in the replication_cluster.

    formula : str
        R (patsy) style formula. Must contain 'methylation': e.g.:
        methylation ~ age + gender + race

    coef : str
        The coefficient of interest in the model, e.g. 'age'

    model_fn : fn
        A function with signature
        fn(formula, methylation, covs, coef, kwargs)
        that returns a dictionary with at least p-value and coef

    pool : multiprocessing.Pool
        a pool to re-use.

    kwargs : dict
        additional args sent to `model_fn`
    """
    print model_fn.func_name
    if kwargs is None: kwargs = {}

    assert isinstance(discovery_clusters, list)
    if not isinstance(replication_feats, list):
        sys.stderr.write("warning: feats is not a list. we will exhaust a \
                generator on first run\n")

    # make replication clusters that are the same as the discovery_clusters
    from .regions import region_cluster_gen
    def cluster_bedder(clusts):
        for c in clusts:
            yield (c[0].chrom, str(c[0].position - 1), str(c[-1].position))

    replication_clusters = list(region_cluster_gen(
            replication_feats,
            cluster_bedder(discovery_clusters)))

    assert len(replication_clusters) == len(discovery_clusters)

    dclusters = list(model_clusters(discovery_clusters, discovery_covs,
                                    formula, coef, model_fn=model_fn, pool=pool,
                                    **kwargs))
    tot_time = sum(c['time'] for c in dclusters)

    rclusters = list(model_clusters(replication_clusters, replication_covs,
                                    formula, coef, model_fn=model_fn, pool=pool,
                                    **kwargs))
    tot_time += sum(c['time'] for c in rclusters)
    r = dict(method=model_fn.func_name, formula=formula,
        time=tot_time)
    r['replication_p'] = np.array([c['p'] for c in rclusters])
    r['discovery_p'] = np.array([c['p'] for c in dclusters])
    assert len(r['replication_p']) == len(r['discovery_p'])
    return r


def evaluate_regions(clust_list, true_regions, df, formula, coef, model_fn,
        pool=None, kwargs=None):
    """
    Evaluate the accuracy of a method (`model_fn`) by comparing how many
    fall into a set of "truth" regions vs. outside.

    Parameters
    ----------

    clust_list : list
        list of clusters

    true_regions : file
        BED file of regions we expect to find DMRs.

    df: pandas.DataFrame
        rontains columns of the covariates listed in formula. Rows
        must be in the same order as they are in clust_list

    formula : str
        R (patsy) style formula. Must contain 'methylation': e.g.:
        methylation ~ age + gender + race

    coef : str
        The coefficient of interest in the model, e.g. 'age'

    model_fn : fn
        A function with signature
        fn(formula, methylation, covs, coef, kwargs)
        that returns a dictionary with at least p-value and coef

    kwargs: dict
        extra arguments sent to model_fn and a sub-dict of plot_kwargs
        to send to the plotting if ax is not None
    """
    assert isinstance(clust_list, list)

    regions = defaultdict(list)
    for i, toks in enumerate(ts.reader(true_regions, header=False)):
        # see if it's a header.
        if i == 0 and not (toks[1] + toks[2]).isdigit(): continue
        # seen used keep track of the regions we've found
        chrom, start, end = toks[0], int(toks[1]), int(toks[2])
        regions[chrom].append((start, end, (chrom, start, end)))

    def is_in(c, regions=regions):
        r = regions[c[0].chrom]
        return any(s <= c[-1].position and e >= c[0].position for s, e, _ in r)

    out_clusters = [c for c in clust_list if not is_in(c)]
    in_clusters = [c for c in clust_list if is_in(c)]
    print "%i clusters in regions of interest" % len(in_clusters)
    print "%i clusters outside regions of interest" % len(out_clusters)
    return evaluate_method(in_clusters + out_clusters,
                           len(in_clusters),
                           df, formula, coef, model_fn, pool=pool,
                           kwargs=kwargs)

def evaluate_method(clust_list, n_true, df, formula, coef, model_fn, pool=None,
        kwargs=None):
    """
    Evaluate the accuracy of a method (`model_fn`) by see checking
    the number of DMRs found at various cutoffs for the real data as sent
    in and data shuffled as in A-clustering

    Parameters
    ----------

    clust_list : list
        list of clusters

    n_true : int
        number of clusters at start of `clust_list` that are true positives

    df: pandas.DataFrame
        rontains columns of the covariates listed in formula. Rows
        must be in the same order as they are in clust_list

    formula : str
        R (patsy) style formula. Must contain 'methylation': e.g.:
        methylation ~ age + gender + race

    coef : str
        The coefficient of interest in the model, e.g. 'age'

    model_fn : fn
        A function with signature
        fn(formula, methylation, covs, coef, kwargs)
        that returns a dictionary with at least p-value and coef

    kwargs: dict
        extra arguments sent to model_fn and a sub-dict of plot_kwargs
        to send to the plotting if ax is not None
    """
    if kwargs is None: kwargs = {}
    plot_kwargs = kwargs.pop('plot_kwargs', {})
    ax = kwargs.pop('ax', None)
    print model_fn.func_name, kwargs

    clusters = model_clusters(clust_list, df, formula, coef,
                              model_fn=model_fn, pool=pool, **kwargs)

    trues, falses = [], []

    tot_time = 0
    for i, c in enumerate(clusters):
        tot_time += c['time']
        (trues if i < n_true else falses).append(c['p'])

    r = dict(method=model_fn.func_name, formula=formula, time=tot_time)

    # find number less than each alpha
    for e in range(8):
        v = 10**-e
        r['true_%i' % e] = sum(t <= v for t in trues)
        r['false_%i' % e] = sum(f <= v for f in falses)

    r['null-ps'] = np.array(falses) # to get sense of distributions
    r['true-ps'] = np.array(trues)
    r['ps'] = np.array(trues + falses)
    r['label'] = plot_kwargs.pop('label',
             r['method'].replace('_cluster', '').replace('_', ' '))
    if ax is not None:
        plot_roc(ax, r, plot_kwargs)
    return r

def plot_alphas(axs, results, colors=None, cutoff=1e-2, labels=None):
    """
    Plot true and false positives from :func:`~evaluate_method`.
    This gives a different view from the ROC plot because it
    may be important to know the performance at a particular \$alpha$
    cutoff.

    Parameters
    ----------

    axs: matplotlib axes
        axs must be of length 2 from plt.subplots(2, 1)

    results: list
        list of results from evaluate_method

    colors: list
         list of colors the same length as results

    cutoff : float
         only plot powers of negative 10 less than this

    """
    if colors is None:
        colors = sns.color_palette("Set1", len(results))
    if labels is None:
        labels = [r.get('label', r['method']) for r in results]

    assert len(labels) == len(results)

    dr = pd.DataFrame(results)
    dr = pd.melt(dr,
        id_vars=[c for c in dr.columns if not ('false_' in c or 'true_' in c)],
        value_vars=[c for c in dr.columns if 'false_' in c or 'true_' in
            c], value_name='n_lt_alpha')

    dr['alpha'] = [10**-int(x.split("_")[1]) for x in dr['variable']]
    dr['truth'] = [x.split('_')[0] for x in dr['variable']]

    # only care about stuff with at least p < cutoff
    r = dr[dr.alpha < cutoff]
    xmin = 1 + int(-np.log10(cutoff))

    for j, truth in enumerate(("true", "false")):
        shapes = []
        for i, result in enumerate(results):
            xx = -np.log10(r.alpha[(r.truth == truth) & (r.method == result['method'])])
            y = r.n_lt_alpha[(r.truth == truth) & (r.method == result['method'])]
            f = axs[j].bar(left=xx+i/9. - 0.28, height=y, width=0.1, fc=colors[i], ec=colors[i])
            shapes.append(f[0])
            axs[j].set_xlim(xmin - 0.5, 7.5)

    axs[1].legend(shapes, labels)
    axs[0].set_xticks([])
    axs[1].set_xticks(xx.unique())
    axs[1].set_xticklabels(['1e-%s' % x for x in range(xmin, 8)])
    axs[1].set_xlabel('$a$')

    axs[0].set_ylabel('true positivies')
    axs[1].set_ylabel('false positives')


def plot_roc(ax, r, plot_kwargs={}):
    """
    Plot ROC for a given result.

    Parameters
    ----------

    ax: matplotlib axis

    r: dict
        return value from :func:`~evaluate_method`
    """
    from sklearn.metrics import roc_curve, auc
    truth = np.array([1] * len(r['true-ps']) + [0] * len(r['null-ps']))
    vals = 1 - r['ps']
    vals = vals[~np.isnan(vals)]
    truth = truth[~np.isnan(r['ps'])]
    fpr, tpr, _ = roc_curve(truth, vals)
    label = ("AUC: %.4f | " % auc(fpr, tpr)) + r['label']
    ax.plot(fpr[1:], tpr[1:], label=label, **plot_kwargs)
    ax.set_xlabel('1 - specificity')
    ax.set_ylabel('sensitivity')

def plot_pvalue_grid(results):
    from matplotlib import pyplot as plt
    # http://web.stanford.edu/~mwaskom/software/seaborn/tutorial/axis_grids.html
    grid = np.array([r['null-ps'] for r in results]).T
    grid = pd.DataFrame(grid, columns=[r.get('label', r['method']) for r in results])

    g = sns.PairGrid(grid, diag_sharey=True)#, hue='truth')

    #g.map_diag(plt.hist)
    g.map_offdiag(lambda x, y, *args, **kwargs: plt.scatter(-np.log10(x),
                                                            -np.log10(y), *args, **kwargs))
    return g


def plot_times(ax, results, colors=None):
    """
    Plot the times taken for each result in results
    from :func:`evalute_method`

    Parameters
    ----------

    ax: matplotlib axis

    results: list of dict
        list of return values from :func:`~evaluate_method`

    colors: list
        list of colors for matplotlib

    """

    if colors is None:
        colors = sns.color_palette("Set1", len(results))

    tmax = int(0.5 + np.log10(1 + max(float(m['time']) for m in results) / 60.))
    for i, m in enumerate(results):
        t = float(m['time']) / 60.
        s = ax.bar(i + 0.08, height=np.log10(1 + t), width=0.85, fc=colors[i], label=m['label'] + ' - %.1f min.' % t)
    ax.legend(loc="best")

    ax.set_yticklabels(["%i" % (10**i) for i in range(tmax)])
    ax.set_title('CPU time (minutes)')
    ax.set_xticks([])


def plot_replication_roc(ax, reps, colors=None, labels=None, alpha_discovery=0.05):
    """
    Plot the ROC curve of results from from :func:`evalute_replication`

    Parameters
    ----------

    ax: matplotlib axis

    reps: list of dict
        list of return values from :func:`~evaluate_replication`

    colors: list
        list of colors for matplotlib

    labels: list
        list of labels to use. same length as `reps`.

    alpha_discovery : float
        alpha cutoff to use a "true" positive in discovery cohort.
    """
    if colors is None:
        colors = sns.color_palette("Set1", len(reps))

    from sklearn.metrics import roc_curve, auc
    for i, rep in enumerate(reps):
        d, r = rep['discovery_p'], rep['replication_p']
        fpr, tpr, _ = roc_curve(d < alpha_discovery, 1 - r)
        label = "AUC: %.4f - %s" % (auc(fpr, tpr), labels[i] if labels else reps['method'].replace("_", " "))
        ax.plot(fpr, tpr, label=label, color=colors[i])
    ax.legend(loc="best")
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.5)
