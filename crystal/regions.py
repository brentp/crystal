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
    >>> feats = [Feature('chr1', i, []) for i in [1, 10, 20, 25, 200, 1000, 1001]]
    >>> with open('/tmp/zxzz.bed', 'w') as fh:
    ...     fh.write('chr1\t5\t24\nchr1\t25\t26\nchr1\t999\t1002\n')

    >>> for cluster in region_cluster_gen(feats, fh.name):
    ...    print cluster
    [Feature(chr1:10), Feature(chr1:20)]
    [Feature(chr1:25)]
    [Feature(chr1:1000), Feature(chr1:1001)]
    """

    regions = defaultdict(list)

    for i, toks in enumerate(ts.reader(regions_bed, header=False)):
        # see if it's a header.
        if i == 0 and not (toks[1] + toks[2]).isdigit(): continue
        regions[toks[0]].append((int(toks[1]), int(toks[2]), "-".join(toks)))
    regions = dict(regions)
    # TODO: handle overlapping regions

    dmr = []
    for f in features:
        cr = regions.get(f.chrom)
        if cr is None:
            if dmr != []:
                yield [x[0] for x in dmr]
                dmr = []
            continue

        try:
            reg = next(rid for (s, e, rid) in cr if s <= f.position <= e)
        except StopIteration:
            continue

        if dmr != [] and dmr[-1][-1] != reg:
            yield [x[0] for x in dmr]
            dmr = []

        dmr.append((f, reg))

    if dmr != []: yield [x[0] for x in dmr]

if __name__ == "__main__":
    import doctest
    doctest.testmod()
