"""Plotting functions for Clusters.

Note that some visualizations work best with datasets that
have fewer samples while some are more informative with more
samples.

.. testsetup:: *

    import crystal
    import crystal.utils as cu

"""


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

def logit(a): return np.log(a / (1 - a))
def ilogit(m): return 1.0 / (1 + np.exp(-m))

sns.set(style="white", context="talk")
colors = sns.color_palette("Set1", 8)

def _plot_continuous(feature, var, ax):
    ax.scatter(ilogit(feature.values), var, s=(2 if len(var) > 40 else 4),
            c=colors[0])

def _plot_dichotomous(feature, var, ax, normed=False):
    var = np.asarray(var)
    cats = np.unique(var)
    cats.sort()

    xvals = []
    for i, cat in enumerate(cats):
        vals = ilogit(feature.values[var == cat])
        xvals.append(vals)

    ax.hist(xvals, alpha=0.8, color=colors[:len(cats)], label=list(cats), normed=normed)

def is_dichotomous(col):
    dichotomous = not np.issubdtype(col, float) or len(col.unique()) == 2
    return dichotomous

def plot_cluster(cluster, covs, normed=False):
    """
    Plot a cluster (output from `crystal.model_cluster`)

    Plot will vary depending on if cluster['var'] is dichotomous
    or continuous.

    .. plot::
        :include-source: true

        >>> import crystal
        >>> import crystal.utils as cu
        >>> covs, cluster = cu.real_cluster()
        >>> formula = "methylation ~ age + gender"
        >>> c = crystal.wrapper(crystal.zscore_cluster, formula, cluster, covs, "gender")
        >>> crystal.plot.plot_cluster(c, covs)

    We will get a different plot if the cluster is a continuous variable.
    """

    dichotomous = is_dichotomous(covs[cluster['var']])
    fig, axs = plt.subplots(cluster['n_sites'], sharey=not dichotomous)

    for i, f in enumerate(cluster['cluster']):
        ax = axs[i]
        ax.set_xlim(0, 1)
        ax.set_xticklabels([])

        ax.locator_params(tight=True, nbins=4)
        #ax.set_yticks(np.arange(0, 1.001, 0.25))
        ax.set_title(f.spos)

        if dichotomous:
            _plot_dichotomous(f, covs[cluster['var']], ax, normed=normed)
            if normed:
                ax.set_yticklabels([])
        else:
            _plot_continuous(f, covs[cluster['var']], ax)
            ax.set_ylabel(cluster['var'])

    ax.set_xticks(np.arange(0, 1.001, 0.25))
    ax.set_xticklabels(np.arange(0, 1.001, 0.25))
    if dichotomous:
        axs[0].legend(loc="best")
        l = axs[0].get_legend()
        l.set_frame_on(True)
        l.get_frame().set_facecolor('white')
        l.get_frame().set_alpha(0.5)
    plt.tight_layout()
    return fig, axs

def factorplot_cluster(cluster, cov, palette='Set1', ilogit=False):
    # http://web.stanford.edu/~mwaskom/software/seaborn/tutorial/categorical_linear_models.html
    """
    Create a factor plot.

    .. plot::
        :include-source: true

        >>> import crystal
        >>> import crystal.utils as cu
        >>> covs, cluster = cu.real_cluster()
        >>> formula = "methylation ~ age + gender"
        >>> c = crystal.wrapper(crystal.zscore_cluster, formula, cluster, covs, "gender")
        >>> crystal.plot.factorplot_cluster(c, covs)

    """
    features = cluster['cluster']
    methylation = np.array([f.values for f in features]).T
    if ilogit:
        methylation = 1 / (1 + np.exp(-methylation))
    df = pd.DataFrame(methylation, columns=[f.spos for f in features])
    df['id'] = [str(x) for x in cov['id']]
    var = cluster['var']
    df[var] = [str(x) for x in cov[var]]
    df = pd.melt(df, id_vars=(var, 'id'), var_name='position')

    fp = sns.factorplot('position', 'value', var, df, kind='box',
            palette=palette, legend=False)
    plt.draw()
    ax = fp.fig.axes[0]
    ax.legend(loc='best')
    return fp

def barplot_cluster(cluster, covs, normed=False, n_bins=50):
    """
    Make a barplot of a cluster. Only works for 2-class data,
    e.g. gender or case-control

    .. plot::
        :include-source: true

        >>> import crystal
        >>> import crystal.utils as cu
        >>> covs, cluster = cu.real_cluster()
        >>> formula = "methylation ~ age + gender"
        >>> c = crystal.wrapper(crystal.zscore_cluster, formula, cluster, covs, "gender")
        >>> crystal.plot.barplot_cluster(c, covs)
    """
    group = np.array(covs[cluster['var']])
    grps = sorted(np.unique(group))
    assert len(grps) == 2

    fig, ax = plt.subplots(1)

    # get min and max for all features so we can use same scale.
    dmin = ilogit(min(f.values.min() for f in cluster['cluster']))
    dmax = ilogit(max(f.values.max() for f in cluster['cluster']))

    for i, feature in enumerate(cluster['cluster']):
        g0 = ilogit(feature.values[group == grps[0]])
        g1 = ilogit(feature.values[group == grps[1]])

        shape0, max0 = half_horizontal_bar(g0, i , ax, True,
                                     dmin=dmin, dmax=dmax, edgecolor='0.4',
                                     facecolor=colors[0], n_bins=n_bins)
        shape1, max1 = half_horizontal_bar(g1, i, ax, False,
                                     dmin=dmin, dmax=dmax, edgecolor='0.4',
                                     facecolor=colors[1], n_bins=n_bins,
                                     norm=max0)

    ax.set_xticks(range(len(cluster['cluster'])))
    ax.set_xticklabels([f.position for f in cluster['cluster']])
    l = ax.legend((shape0, shape1),
            ("%s - %s" % (cluster['var'], grps[0]),
            ("%s - %s" % (cluster['var'], grps[1]))), fancybox=True, loc='best')
    l.set_frame_on(True)
    l.get_frame().set_facecolor('white')
    l.get_frame().set_alpha(0.5)
    ax.xaxis.grid(linewidth=0.25, color="0.02")

    return fig, ax

def half_horizontal_bar(data, pos, ax, left=False, dmin=0, dmax=1, n_bins=70,
        norm=None, **kwargs):

    bins = np.linspace(dmin, dmax, n_bins + 1)
    counts, edges = np.histogram(data, bins=bins, density=True)
    counts = (0 + counts)
    cnorm = counts.sum()
    edges = edges[:n_bins]

    bsize = edges[1] - edges[0]
    counts /= (2.5 * float(counts.max()))

    keep = counts > 0 # dont draw 0-height bars.
    counts = counts[keep]
    edges = edges[keep]

    if norm is not None:

        counts *= norm / cnorm

    if left:
        counts *= -1

    pos += (-0.0002 if left else 0.0002)
    return ax.barh(edges, counts, bsize, left=pos, **kwargs)[0], cnorm


def spaghetti_plot(cluster, cov, ax=None, ilogit=False, palette='Set1'):
    """Create a spaghetti plot of a modeled cluster. This is best when the
       number of samples is less than about 20. Otherwise, use
       :func:`~plot_cluster`

    .. plot::
        :include-source: true

        >>> import crystal
        >>> import crystal.utils as cu
        >>> covs, cluster = cu.real_cluster()
        >>> formula = "methylation ~ age + gender"
        >>> c = crystal.wrapper(crystal.zscore_cluster, formula, cluster, covs, "gender")
        >>> crystal.plot.spaghetti_plot(c, covs)
    """
    from pandas.tools.plotting import parallel_coordinates
    from crystal import CountFeature
    features = cluster['cluster']
    methylation = np.array([f.values for f in features]).T
    if ilogit:
        methylation = 1 / (1 + np.exp(-methylation))
    df = pd.DataFrame(methylation, columns=[f.spos for f in features])
    var = cluster['var']
    df[var] = [str(x) for x in cov[var]]

    mmax = methylation.max().max()
    mmin = methylation.min().min()
    if ax is None:
        fig, ax = plt.subplots(1)

    colors=sns.color_palette(palette)
    ax = parallel_coordinates(df, var, color=colors, ax=ax, use_columns=False)
    # pandas adds dark axvline, this is to remove that.
    lines = ax.get_lines()
    for i in range(len(features)):
        lines.pop().remove()
    lbls = ax.get_legend().get_texts()
    l = ax.get_legend()
    l.set_frame_on(True)
    l.get_frame().set_facecolor('white')
    l.get_frame().set_alpha(0.5)

    if isinstance(features[0], CountFeature):
        counts = np.array([f.counts for f in features]).T
        for icol, f in enumerate(features):
            for group in sorted(df[var].unique()):
                ax.scatter([icol] * sum(df[var] == group),
                       methylation[df[var] == group, icol],
                       edgecolors=colors[icol],
                       facecolors=colors[j],
                       alpha=0.5,
                       s=counts[df[var] == group, icol])
    plt.draw()
    xmin, xmax = ax.get_xlim()
    ax.set_xlim(int(xmin) - 0.05, int(xmax) + 0.05)
    ax.get_legend().set_title(var)
    ax.set_ylabel('methylation')
    sns.axes_style({'axes.linewidth': 0, 'axes.grid': False})
    sns.despine()
    sns.set_style("ticks")
    return ax
