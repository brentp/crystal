import crystal
import crystal.utils
import tempfile
import sys
import StringIO as io


def test_roc_bed():
    with open('/tmp/p.bed', 'w') as p_bed:
        p_bed.write("""\
chr1\t1\t1\t1e-3
chr1\t3\t3\t1e-3""")
    with open('/tmp/r.bed', 'w') as r_bed:
        r_bed.write("""\
chr1\t1\t1""")

    res = crystal.utils.roc_out(p_bed.name, 4, r_bed.name)
    assert len(res) == 2
    assert len(res[0]) == 2
    assert len(res[1]) == 2

def test_write_cluster():
    import locale
    covs, cluster = crystal.utils.example_random_cluster()
    fh = io.StringIO()
    crystal.utils.write_cluster(cluster, fh)
    s = fh.getvalue()
    assert cluster[0].chrom in s, s
    assert s.startswith(cluster[0].chrom + ":" + str(cluster[0].position) + "\t"), s

       

