from collections import defaultdict
import toolshed as ts

def region_cluster_gen(features, regions_bed):
    r"""Generate clusters from regions defined in a BED file.

    Useful for re-testing existing regions--no need to use a
    clustering method, just yield 'clusters' of sites that fall
    within the regions in regions_bed.

    features: a generator of features
    regions_bed: a bed file

    regions are read into memory, but feature_gen is sorted.

    >>> from crystal import Feature
    >>> feats = [Feature('chr1', i, []) for i in [1, 10, 20, 200, 1000, 1001]]
    >>> with open('/tmp/zxzz.bed', 'w') as fh:
    ...     fh.write('chr1\t5\t25\nchr1\t999\t1002\n')

    >>> for cluster in region_cluster_gen(feats, fh.name):
    ...    print cluster
    [Feature(chr1:10), Feature(chr1:20)]
    [Feature(chr1:1000), Feature(chr1:1001)]
    """

    regions = defaultdict(list)

    for i, toks in enumerate(ts.reader(regions_bed, header=False)):
        # see if it's a header.
        if i == 0 and not (toks[1] + toks[2]).isdigit(): continue
        regions[toks[0]].append((int(toks[1]), int(toks[2])))
    regions = dict(regions)

    dmr = []
    for f in features:
        cr = regions.get(f.chrom)
        if cr is None:
            if dmr != []:
                yield dmr
                dmr = []
            continue

        if any(s <= f.position <= e for (s, e) in cr):
            dmr.append(f)
        elif dmr != []:
                yield dmr
                dmr = []

    if dmr != []: yield dmr

if __name__ == "__main__":
    import doctest
    doctest.testmod()
