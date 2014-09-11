from __future__ import print_function
from collections import defaultdict
import toolshed as ts
import sys

def region_cluster_gen(features, regions_bed):
    r"""Generate clusters from regions defined in a BED file.

    Useful for re-testing existing regions--no need to use a
    clustering method, just yield 'clusters' of sites that fall
    within the regions in regions_bed.

    features: a generator of features
    regions_bed: a bed file

    regions are read into memory, but feature_gen is sorted.

    >>> from crystal import Feature
    >>> feats = [Feature('chr1', i, []) for i in [1, 10, 20, 25, 200, 1000, 1001]]
    >>> feats.append(Feature('chr2', 999, []))
    >>> with open('/tmp/zxzz.bed', 'w') as fh:
    ...  fh.write('chr1\t5\t25\nchr1\t25\t26\nchr1\t999\t1002\nchr2\t999\t1002\n')

    >>> for cluster in region_cluster_gen(feats, fh.name):
    ...    print(cluster)
    [Feature(chr1:25)]
    [Feature(chr1:10), Feature(chr1:20), Feature(chr1:25)]
    [Feature(chr1:1000), Feature(chr1:1001)]
    [Feature(chr2:999)]
    """

    regions = defaultdict(list)
    seen = set()

    for i, toks in enumerate(ts.reader(regions_bed, header=False)):
        # see if it's a header.
        if i == 0 and not (toks[1] + toks[2]).isdigit(): continue
        # seen used keep track of the regions we've found
        chrom, start, end = toks[0], int(toks[1]), int(toks[2])
        seen.add((chrom, start, end))
        regions[chrom].append((start, end, (chrom, start, end)))
    regions = dict(regions)
    # TODO: handle overlapping regions

    dmrs = {}
    for f in features:
        cr = regions.get(f.chrom)
        if cr is None:
            for reg in dmrs.keys():
                seen.remove(reg)
                yield dmrs[reg]
                dmrs.pop(reg)
            continue

        # loop over all regions of overlap and add the feature
        for reg in (rid for (s, e, rid) in cr if s <= f.position <= e):
            if not reg in dmrs: dmrs[reg] = []
            dmrs[reg].append(f)

        # loop all existing regions and see if the current feature start >
        # region end. if so, we can pop and yield that region.
        for rchrom, rstart, rend in dmrs.keys():
            if rchrom != f.chrom or rend < f.position:
                seen.remove((rchrom, rstart, rend))
                yield dmrs.pop((rchrom, rstart, rend))

    for reg in dmrs.keys():
        seen.remove(reg)
        yield dmrs.pop(reg)

    if seen:
        print(seen, file=sys.stderr)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
