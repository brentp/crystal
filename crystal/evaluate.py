from .crystal import model_clusters
import numpy as np
import sys

def evaluate_replication(discovery_clusters, replication_feats,
                         discovery_covs, replication_covs, formula,
                         coef, model_fn, pool=None, kwargs=None):
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

def plot_alphas(axs, results):
    """
    axs must be of length 2 from plt.subplots(2, 1)
    """
    dr = pd.DataFrame(results)
    dr = pd.melt(dr,
        id_vars=[c for c in dr.columns if not ('false_' in c or 'true_' in c)],
        value_vars=[c for c in dr.columns if 'false_' in c or 'true_' in
            c], value_name='n_lt_alpha')

    dr['alpha'] = [10**-int(x.split("_")[1]) for x in dr['variable']]
    dr['truth'] = [x.split('_')[0] for x in dr['variable']]

    # only care about stuff with at least p < 0.01
    r = dr[dr.alpha < 0.01]

    for j, truth in enumerate(("true", "false")):
        shapes = []
        for i, m in enumerate(methods):
            xx = -np.log10(r.alpha[(r.truth == truth) & (r.method == m.func_name)])
            y = r.n_lt_alpha[(r.truth == truth) & (r.method == m.func_name)]
            f = axs[j].bar(left=xx+i/9. - 0.28, height=y, width=0.1, fc=colors[i], ec=colors[i])
            shapes.append(f[0])
            axs[j].set_xlim(2.5, 7.5)

    axs[1].legend(shapes, clean_methods)
    axs[1].set_xticks(xx.unique())
    axs[1].set_xticklabels(['1e-%s' % x for x in range(3, 8)])
    axs[1].set_xlabel('$a$')

    axs[1].set_ylabel('false positives')
    axs[0].set_ylabel('true positivies')
    return fig, axs




def plot_roc(ax, r, plot_kwargs):
    from sklearn.metrics import roc_curve, auc
    truth = np.array([1] * len(r['true-ps']) + [0] * len(r['null-ps']))
    vals = 1 - r['ps']
    vals = vals[~np.isnan(vals)]
    truth = truth[~np.isnan(truth)]
    fpr, tpr, _ = roc_curve(truth, vals)
    label = ("AUC: %.4f | " % auc(fpr, tpr)) + r['label']
    ax.plot(fpr[1:], tpr[1:], label=label, **plot_kwargs)
    ax.set_xlabel('1 - specificity')
    ax.set_ylabel('sensitivity')

def plot_times(ax, results, colors):
    tmax = int(0.5 + np.log10(1 + max(float(m['time']) for m in results) / 60.))
    for i, m in enumerate(results):
        t = float(m['time']) / 60.
        s = ax.bar(i + 0.08, height=np.log10(1 + t), width=0.85, fc=colors[i], label=m['label'] + ' - %.1f min.' % t)
    ax.legend(loc="best")

    ax.set_yticklabels(["%i" % (10**i) for i in range(tmax)])
    ax.set_title('CPU time (minutes)')
    ax.set_xticks([])


def plot_replication_roc(ax, reps, colors, labels=None, alpha_discovery=0.05):
    from sklearn.metrics import roc_curve, auc
    for i, rep in enumerate(reps):
        d, r = rep['discovery_p'], rep['replication_p']
        fpr, tpr, _ = roc_curve(d < alpha_discovery, 1 - r)
        label = "AUC: %.4f - %s" % (auc(fpr, tpr), labels[i] if labels else reps['method'].replace("_", " "))
        ax.plot(fpr, tpr, label=label, color=colors[i])
    ax.legend(loc="best")
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.5)
