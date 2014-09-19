import pandas as pd
import numpy as np
import crystal

np.random.seed(10)
covs = pd.DataFrame({'gender': ['F'] * 10 + ['M'] * 10,
                     'age': np.random.uniform(10, 25) })

methylation = np.random.normal(size=(5, covs.shape[0]))

cluster = [crystal.Feature('chr1', i* 10, m) for i, m in enumerate(methylation)]

formula = "methylation ~ age + gender"
formula_cpg = "methylation ~ age + gender + CpG"


def check_model(m, formula, coef):
    """testing models on cluster"""
    r = m(formula, cluster, covs, coef)
    assert 'p' in r
    assert 't' in r
    assert coef in r['covar']
    print m.func_name, r
    assert 'coef' in r

def check_wrapper(m, formula, coef):

    r = crystal.wrapper(m, formula, cluster, covs, coef)
    for v in ('var', 'n_sites', 'sites', 'time', 'chrom', 'start', 'end'):
        assert v in r

def test_one():

    r = crystal.one_cluster(formula, cluster, covs, "age")
    for k in ('p', 't', 'covar', 'coef'):
        assert k in r, r

def test_models():
    """testing models on cluster"""
    for coef in ('age', 'gender'):
        for m, form in ((crystal.zscore_cluster, formula),
                  (crystal.liptak_cluster, formula),
    #              (crystal.bump_cluster, formula),
                  (crystal.gee_cluster, formula_cpg),
                  (crystal.mixed_model_cluster, formula_cpg),
                  (crystal.glsar_cluster, formula_cpg),
                  (crystal.ols_cluster_robust, formula_cpg)):
            yield check_model, m, form, coef
            yield check_wrapper, m, form, coef

def test_model_clusters():

    r = list(crystal.model_clusters([cluster], covs, formula, "age",
        crystal.zscore_cluster, n_cpu=1))
    assert len(r) == 1
    r = r[0]
    for v in ('var', 'n_sites', 'sites', 'time', 'chrom', 'start', 'end'):
        assert v in r


def test_long_covs():

    dflong = crystal.long_covs(covs, methylation)
    assert dflong.shape[0] == covs.shape[0] * len(methylation)

    assert len(dflong['id'].unique()) == len(methylation[0])
    assert len(dflong['CpG'].unique()) == len(methylation)