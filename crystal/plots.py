import matplotlib.pyplot as plt
import numpy as np

def logit(a): return np.log(a / (1 - a))
def ilogit(m): return 1.0 / (1 + np.exp(-m))

import seaborn as sns
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

def plot_cluster(cluster, covs, normed=False):

    dichotomous = not np.issubdtype(covs[cluster['var']], float) or len(covs[cluster['var']].unique()) == 2
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
