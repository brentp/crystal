import crystal
from crystal.simulate import rr_cluster, simulate_cluster, simulate_regions
import crystal.utils


def test_rr_cluster():
    covs, cluster = crystal.utils.real_cluster()

    shuffled = rr_cluster(cluster, covs, formula="methylation ~ age")
    assert len(shuffled) == len(cluster)
    assert len(shuffled[0].values) == len(cluster[0].values)
    assert shuffled[0].values[0] != cluster[0].values[0]
    assert shuffled[0].position == cluster[0].position

def test_rr_count_cluster():
    covs, cluster = crystal.utils.real_count_cluster()

    shuffled = rr_cluster(cluster, covs, formula="methylation ~ Eos")
    assert len(shuffled) == len(cluster)
    assert len(shuffled[0].values) == len(cluster[0].values)
    assert len(shuffled[0].counts) == len(cluster[0].counts)
    assert shuffled[0].values[0] != cluster[0].values[0]
    assert shuffled[0].counts[0] == cluster[0].counts[0]
    assert shuffled[0].methylated[0] == cluster[0].methylated[0]

    assert shuffled[0].position == cluster[0].position

def test_rr_cluster():
    covs, cluster = crystal.utils.real_cluster()

    c2 = simulate_cluster(cluster)
    assert c2[0].values[0] != cluster[0].values[0]

def test_rr_count_cluster():
    covs, cluster = crystal.utils.real_count_cluster()

    c2 = simulate_cluster(cluster)
    assert c2[0].values[0] != cluster[0].values[0]
    assert c2[0].counts[0] != cluster[0].counts[0]
    assert c2[0].methylated[0] != cluster[0].methylated[0]

def test_simulate_region():
    covs, cluster = crystal.utils.real_cluster()
    
    r = open('/tmp/_rr.bed', 'w')
    clusts = list(simulate_regions([cluster], r, sizes={}))
    assert len(clusts) == 1
