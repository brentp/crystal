"""Plotting functions for Clusters."""


import matplotlib.pyplot as plt
import numpy as np
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
        else:
            _plot_continuous(f, covs[cluster['var']], ax)
            ax.set_ylabel(cluster['var'])

    ax.set_xticks(np.arange(0, 1.001, 0.25))
    ax.set_xticklabels(np.arange(0, 1.001, 0.25))
    if dichotomous:
        axs[0].legend()
    return fig, axs


def barplot_cluster(cluster, covs, normed=False, n_bins=50):
    # this only works for 2-class data.
    group = np.array(covs[cluster['var']])
    grps = sorted(np.unique(group))
    assert len(grps) == 2

    fig, ax = plt.subplots(1)

    # get min and max for all features so we can use same scale.
    dmin = min(f.values.min() for f in cluster['cluster'])
    dmax = max(f.values.max() for f in cluster['cluster'])

    for i, feature in enumerate(cluster['cluster']):
        g0 = ilogit(feature.values[group == grps[0]])
        g1 = ilogit(feature.values[group == grps[1]])

        shape0 = half_horizontal_bar(g0, i , ax, True,
                                     dmin=dmin, dmax=dmax, edgecolor='0.4',
                                     facecolor=colors[0], n_bins=n_bins)
        shape1 = half_horizontal_bar(g1, i, ax, False,
                                     dmin=dmin, dmax=dmax, edgecolor='0.4',
                                     facecolor=colors[1], n_bins=n_bins)

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

def half_horizontal_bar(data, pos, ax, left=False, dmin=0, dmax=1, n_bins=70, **kwargs):

    bins = np.linspace(dmin, dmax, n_bins + 1)
    counts, edges = np.histogram(data, bins=bins, density=True)
    counts = (0 + counts)
    edges = edges[:n_bins]

    bsize = edges[1] - edges[0]
    counts /= (2.5 * float(counts.max()))

    keep = counts > 0 # dont draw 0-height bars.
    counts = counts[keep]
    edges = edges[keep]


    if left:
        counts *= -1

    pos += (-0.0002 if left else 0.0002)
    return ax.barh(edges, counts, bsize, left=pos, **kwargs)[0]
