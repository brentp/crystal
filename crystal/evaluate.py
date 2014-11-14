from collections import defaultdict
from interlap import InterLap
import toolshed as ts
import itertools as it
import seaborn as sns
sns.set(style="white", context="talk")
colors = sns.color_palette("Set1", 8)

import pandas as pd
import numpy as np
import sys
from sklearn.metrics import roc_curve, auc


def cluster_bediter(modeled_clusters):
    for m in modeled_clusters:
        c = m['cluster']
        yield (c[0].chrom, c[0].position - 1, c[-1].position, m['p'], len(c))

def evaluate(bed_iter, regions, sizes=None, label='', ax=None, **plot_kwargs):
    """
    Evaluate the accuracy of a method (`model_fn`) by comparing how many
    fall into a set of "truth" regions vs. outside.

    Parameters
    ----------

    bed_iter : iterable
        iterable with each element of (chrom, start, end, p-value, n_sites)
        can be regions or features.

    regions : file
        BED file of regions.

    sizes : tuple(int,)
        which size regions to test. default (None) is all sizes.

    label : str
        label for plot legend

    ax : axis

    """

    true_regions = defaultdict(InterLap)
    false_regions = defaultdict(InterLap)
    for i, toks in enumerate(ts.reader(regions, header=False)):
        # see if it's a header.
        if i == 0 and not (toks[1] + toks[2]).isdigit(): continue
        # if they request a specific size, only regions of exactly
        # that length are kept.
        if sizes is not None and int(toks[4]) not in sizes: continue
        # seen used keep track of the regions we've found
        chrom, start, end = toks[0], int(toks[1]), int(toks[2])
        if len(toks) <= 3 or toks[3][:3] == "tru":
            true_regions[chrom].add((start, end))
        else:
            false_regions[chrom].add((start, end))

    seen_true = defaultdict(set)
    seen_false = defaultdict(set)

    def is_in(b, regions, check_regions):
        found = list(regions[b[0]].find((b[1], b[2])))
        for s, e in found:
            check_regions[b[0]].add((s, e))
        return len(found)

    no_false = len(false_regions) == 0
    if no_false: assert sizes is None, ('cant have sizes without false regions')

    ps, truths = [], []
    for b in bed_iter:
        chrom, start, end, p, n_sites = b[:5]
        assert isinstance(start, int)
        n_true = is_in(b, true_regions, seen_true)
        # also keep track of which false regions have been seen
        if no_false:
            n_false = n_sites - n_true
        else:
            n_false = is_in(b, false_regions, seen_false)
            #assert n_true + n_false == n_sites, (n_true, n_false, b)

        # here, we multiply because each region can overlap multiple sites
        truths.extend([1] * n_true)
        truths.extend([0] * n_false)

        ps.extend([p] * (n_true + n_false))

    # add in the missed trues and missed falses.
    #"""
    one = 1.0 - 1e-16
    for chrom in true_regions:
        regs = set(true_regions[chrom])
        seen = seen_true[chrom]
        assert not seen - regs
        missed = regs - seen
        # add a p-value of 1 for missed regions
        ps.extend([one] * len(missed))
        truths.extend([1] * len(missed))

    # blah code duplication.
    for chrom in false_regions:
        regs = set(false_regions[chrom])
        seen = seen_false[chrom]
        assert not seen - regs
        missed = regs - seen

        ps.extend([one] * len(missed))
        truths.extend([0] * len(missed))
    #"""
    truths, ps = np.array(truths), np.array(ps)

    if ax is None:
        return truths, ps

    from sklearn.metrics import classification_report, precision_recall_curve, roc_curve
    subset = ~np.isnan(ps)
    precision, recall, _ = precision_recall_curve(truths[subset], 1.0 - ps[subset])
    fpr, tpr, _ = roc_curve(truths[subset], 1.0 - ps[subset])

    print label, len(truths)
    ax.plot(recall, precision, label=label, **plot_kwargs)
    #ax.plot(fpr, tpr, label=label, **plot_kwargs)
    #ax.set_xlabel('recall')
    #ax.set_ylabel('precision')
    return truths, ps

def write_region_bed(feature_iter, true_regions, out_fh):
    """
    Write a region bed file suitable for use in :func:`~evaluate`.
    given true regions (likely from an external program, otherwise use
    :func:`~write_modeled_regions`).

    Parameters
    ----------

    feature_iter : iterable of Features

    true_regions : file
        BED file containing true regions

    out_fh : filehandle
        where to write the data
    """
    fmt = "{chrom}\t{start}\t{end}\t{truth}\t{size}\n"
    out_fh.write(ts.fmt2header(fmt))

    regions = defaultdict(InterLap)

    for i, toks in enumerate(ts.reader(true_regions, header=False)):
        # see if it's a header.
        if i == 0 and not (toks[1] + toks[2]).isdigit(): continue
        chrom, start, end = toks[0], int(toks[1]), int(toks[2])
        regions[chrom].add((start, end))

    for f in feature_iter:
        truth = 'true' if (f.position, f.position) in regions[f.chrom] else 'false'
        out_fh.write(fmt.format(chrom=f.chrom, start=f.position - 1,
                    end=f.position, truth=truth, size=1))
    out_fh.flush()

def write_modeled_regions(modeled_clusters, p_cutoff, out_fh):
    """
    Write a region bed file suitable for use in :func:`~evaluate_modeled_regions`.

    Parameters
    ----------

    modeled_clusters : list
        output from :func:`~model_clusters`

    p_cutoff : float
        values < this are set as true

    out_fh : filehandle
        where to write the data
    """
    fmt = "{chrom}\t{start}\t{end}\t{truth}\t{size}\n"
    out_fh.write(ts.fmt2header(fmt))
    for mc in modeled_clusters:
        c = mc['cluster']
        truth = ['false', 'true'][int(mc['p'] < p_cutoff)]
        for f in c:
            out_fh.write(fmt.format(**dict(chrom=f.chrom, start=f.position - 1,
                         end=f.position,
                         truth=truth,
                         size=len(c))))
