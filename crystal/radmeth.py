"""
Some functions to make radmeth callable from python
"""
from sklearn.metrics import roc_curve, auc
import itertools as it
import numpy as np
import os
import pandas as pd
import toolshed as ts
from collections import defaultdict
import tempfile
import atexit
import shutil

def rad_format(fmethylated, fcounts, fout):
    if isinstance(fout, basestring):
        fout = ts.nopen(fout, "w")
    for i, (m, c) in enumerate(it.izip(ts.reader(fmethylated, header=False),
                                       ts.reader(fcounts, header=False))):
        if i == 0:
            fout.write("\t" + "\t".join(m[1:]) + "\n")
        else:
            assert m[0] == c[0]
            methyls = m[1:]
            counts = c[1:]
            pairs = "\t".join("%s %s" % (ci, mi) for mi, ci in zip(methyls, counts))
            chrom, pos = c[0].split(":")
            pos = int(pos)
            site = "%s:%i:%i" % (chrom, pos, pos + 1)
            fout.write("%s\t%s\n" % (site, pairs))
    return fout.name

def mktemp():
    _, f = tempfile.mkstemp()
    atexit.register(os.unlink, f)
    return f

def run_radmeth(design, fmethylated, fcounts, dmrs_out, bins="1:200:1"):
    for c in design.columns:
        design[c] = design[c].astype(int)


    designf = mktemp()
    design.to_csv(designf, sep="\t", index=True, index_label=False)

    radf = mktemp()
    rad_format(fmethylated, fcounts, radf)

    cpg = mktemp()

    factor = design.columns[1]
    print factor
    list(ts.nopen("|wand -factor %s %s %s > %s" % (factor, designf, radf, cpg)))
    list(ts.nopen("|adjust -bins %s %s > %s" % (bins, cpg, dmrs_out)))

    return dmrs_out

if __name__ == "__main__":
    import patsy
    design = pd.read_table("/drive/450k/mthfr/covs.Control_Saline-vs-MTHFR-Saline.txt")
    design.index = design.Sample
    X = patsy.dmatrix("~ Group", design, return_type="dataframe")
    dmrs_out = "cpgs.adjust.bed"

    fmethylated = "/drive/450k/mthfr/methylated4.txt.gz"
    fcounts = "/drive/450k/mthfr/counts4.txt.gz"

    bed = run_radmeth(X, fmethylated, fcounts, dmrs_out, bins="1:200:1")

    t, v = roc_out(bed, 5, "../crystal/notebooks/mthfr-region.bed")
    fpr, tpr, _ = roc_curve(t, v)
    print auc(fpr, tpr)
