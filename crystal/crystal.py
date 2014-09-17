"""
Model clustered, correlated data.
"""

import sys
from aclust import mclust
import toolshed as ts
from itertools import starmap, izip
import time

import pandas as pd
import numpy as np
import scipy.stats as ss
import patsy

from scipy.stats import norm
from numpy.linalg import cholesky as chol

from statsmodels.api import GEE, GLM, MixedLM, RLM, GLS, OLS, GLSAR
from statsmodels.genmod.dependence_structures import Exchangeable
from statsmodels.genmod.families import Gaussian
from statsmodels.discrete.discrete_model import NegativeBinomial

def one_cluster(formula, methylation, covs, coef, robust=False):
    """used when we have a "cluster" with 1 probe."""
    c = covs.copy()
    c['methylation'] = methylation
    res = (RLM if robust else GLS).from_formula(formula, data=c).fit()
    return get_ptc(res, coef)

def long_covs(covs, methylation):
    covs['id'] = ['id_%i' % i for i in range(len(covs))]
    cov_rep = pd.concat((covs for i in range(len(methylation))))
    #nr, nc = methylation.shape
    cov_rep['CpG'] = np.repeat(['CpG_%i' % i for i in
        range(methylation.shape[0])], methylation.shape[1])

    #cov_rep['CpG'] = np.repeat(['CpG_%i' % i for i in range(nr)], nc)
    cov_rep['methylation'] = np.concatenate(methylation)
    return cov_rep

def corr(methylations):
    if len(methylations) == 0: return np.nan
    if len(methylations) == 1: return [[1]]
    c = np.abs(ss.spearmanr(methylations.T))
    if len(c.shape) == 1:
        assert len(methylations) == 2
        return np.array([[1, c[0]], [c[0], 1]])
    return c[0]


def nb_cluster(formula, cluster, covs, coef):
    """Model a cluster of correlated features with the negative binomial"""
    methylated = np.array([f.methylated for f in cluster])
    counts = np.array([f.counts for f in cluster])
    try:
        # intentionally flipped so that we model count ~ disease with
        # methylated offset.
        res = [NegativeBinomial.from_formula(formula, covs, offset=count)\
               .fit(disp=0) for methylation, count in izip(methylated, counts)]
    except:
        return dict(t=np.nan, coef=np.nan, covar="NA", p=np.nan, corr=np.nan)

    methylations = methylated / counts # for correlation below.

    r = get_ptc(res, coef)
    nan = np.isnan(r['p'])
    r['p'] = zscore_combine(r['p'][~nan], corr(methylations[~nan]))
    r['t'], r['coef'] = np.mean(r['t']), np.mean(r['coef'])
    return r

def gee_cluster(formula, methylation, covs, coef, cov_struct=Exchangeable(),
        family=Gaussian()):
    """An example of a `model_fn`; any function with a similar signature
    can be used.

    Parameters
    ----------

    formula : str
        R (patsy) style formula. Must contain 'methylation': e.g.:
        methylation ~ age + gender + race

    methylation : numpy.ndarray
        n_sites * n_samples array of (methylation) measurements.

    covs : pandas.DataFrame
        Contains covariates from `formula`

    coef: str
        coefficient of interest, e.g. 'age'

    Returns
    -------

    result : dict
        dict with values (keys) of at least p-value ('p'), coefficient
        estimate ('coef') and any other information desired.
    """

    cov_rep = long_covs(covs, methylation)
    res = GEE.from_formula(formula, groups=cov_rep['id'], data=cov_rep, cov_struct=cov_struct).fit()
    return get_ptc(res, coef)

# see: https://gist.github.com/josef-pkt/89585d0b084739a4ed1c
def ols_cluster_robust(formula, methylation, covs, coef):
    """Model clusters with cluster-robust OLS, same signature as
    :func:`~gee_cluster`"""
    cov_rep = long_covs(covs, methylation)
    res = OLS.from_formula(formula, data=cov_rep).fit(cov_type='cluster',
            cov_kwds=dict(groups=cov_rep['id']))
    return get_ptc(res, coef)

def glsar_cluster(formula, methylation, covs, coef, rho=6):
    cov_rep = long_covs(covs, methylation)
    # group by id, then sort by CpG so that AR is to adjacent CpGs
    cov_rep.sort(['id', 'CpG'], inplace=True)
    res = GLSAR.from_formula(formula, data=cov_rep, rho=rho).iterative_fit(maxiter=5)
    return get_ptc(res, coef)

def gls_cluster(formula, methylation, covs, coef):
    # TODO: currently this is very, very slow
    cov_rep = long_covs(covs, methylation)
    z = np.cov(methylation.T)
    sigma = np.repeat(np.repeat(z, len(methylation), axis=0), len(methylation), axis=1)
    res = GLS.from_formula(formula, data=cov_rep, sigma=sigma).fit()
    return get_ptc(res, coef)

def mixed_model_cluster(formula, methylation, covs, coef):
    """Model clusters with a mixed-model, same signature as
    :func:`~gee_cluster`"""
    cov_rep = long_covs(covs, methylation)
    # TODO: remove this once newer version of statsmodels is out.
    # speeds convergence by using fixed estimates from OLS
    params = OLS.from_formula(formula, data=cov_rep).fit().params

    res = MixedLM.from_formula(formula, groups='id',
            data=cov_rep).fit(start_params=dict(fe=params), reml=False,
                    method='bfgs')

    return get_ptc(res, coef)

def get_ptc(fit, coef):
    if isinstance(fit, list):
        result, res = [get_ptc(f, coef) for f in fit], {}
        for c in ('p', 't', 'coef'):
            res[c] = np.array([r[c] for r in result])
        res['covar'] = result[0]['covar']
        return res

    idx = [i for i, par in enumerate(fit.model.exog_names)
                        if par.startswith(coef)]
    assert len(idx) == 1, ("too many params like", coef)
    try:
        return {'p': fit.pvalues[idx[0]],
            't': fit.tvalues[idx[0]],
            'coef': fit.params[idx[0]],
            'covar': fit.model.exog_names[idx[0]]}
    except ValueError:
        return dict(p=np.nan, t=np.nan, coef=np.nan,
                covar=fit.model.exog_names[idx[0]])


def stouffer_liptak_combine(pvals, sigma):
    """Combine p-values accounting for correlation."""
    qvals = norm.isf(pvals).reshape(len(pvals), 1)
    try:
        C = np.asmatrix(chol(sigma)).I
    except np.linalg.linalg.LinAlgError:
          # for non positive definite matrix default to z-score correction.
          return zscore_combine(pvals, sigma)

    qvals = C * qvals
    Cp = qvals.sum() / np.sqrt(len(qvals))
    return norm.sf(Cp)

def zscore_combine(pvals, sigma):
    if np.all(np.isnan(sigma)): return np.nan
    z, L = np.mean(norm.isf(pvals)), len(pvals)
    sz = 1.0 / L * np.sqrt(L + 2 * np.tril(sigma, k=-1).sum())
    return norm.sf(z / sz)

def _combine_cluster(formula, methylations, covs, coef, robust=False):
    """function called by z-score and liptak to get pvalues"""
    res = [(RLM if robust else GLS).from_formula(formula, covs).fit()
            for methylation in methylations]
    idx = [i for i, par in enumerate(res[0].model.exog_names)
                   if par.startswith(coef)][0]
    pvals = np.array([r.pvalues[idx] for r in res], dtype=np.float64)
    pvals[pvals == 1] = 1.0 - 9e-16
    return dict(t=np.array([r.tvalues[idx] for r in res]),
                coef=np.array([r.params[idx] for r in res]),
                covar=res[0].model.exog_names[idx],
                p=pvals[~np.isnan(pvals)],
                corr=corr(methylations[~np.isnan(pvals)]))

def bayes_cluster():
    pass

def liptak_cluster(formula, methylations, covs, coef, robust=False):
    """Model clusters by fitting model at each site and then
    combining using :func:`~stouffer_liptak`. same signature as
    :func:`~gee_cluster`"""
    r = _combine_cluster(formula, methylations, covs, coef, robust=robust)
    r['p'] = stouffer_liptak_combine(r['p'], r['corr'])
    r['t'], r['coef'] = r['t'].mean(), r['coef'].mean()
    r.pop('corr')
    return r

def zscore_cluster(formula, methylations, covs, coef, robust=False):
    """Model clusters by fitting model at each site and then
    combining using the z-score method. Same signature as
    :func:`~gee_cluster`"""
    r = _combine_cluster(formula, methylations, covs, coef, robust=robust)
    r['p'] = zscore_combine(r['p'], r.pop('corr'))
    r['t'], r['coef'] = r['t'].mean(), r['coef'].mean()
    return r

def wrapper(model_fn, formula, cluster, clin_df, coef, kwargs,
            as_features=False):
    """wrap the user-defined functions to return everything we expect and
    to call just GLS when there is a single probe."""
    t = time.time()
    if len(cluster) > 1:
        r = model_fn(formula,
                cluster if as_features else np.array([c.values for c in cluster]),
                clin_df,
                coef, **kwargs)
    else:
        r = one_cluster(formula, cluster[0] if as_features else cluster[0].values, clin_df, coef)
    r['time'] = time.time() - t
    r['chrom'] = cluster[0].chrom
    r['start'] = cluster[0].position - 1
    r['end'] = cluster[-1].position
    r['n_sites'] = len(cluster)
    r['sites'] = [c.spos for c in cluster]
    r['var'] = coef
    r['cluster'] = cluster
    return r


# function for comparing with bump_cluster
# takes the return value of _combine_cluster and returns a single numeric value
def coef_sum(c, cutoff=0.005):
    coefs = c['coef']
    return sum(min(0, c + cutoff) if c < 0 else max(0, c - cutoff) for c in coefs)

# function for comparing with bump_cluster
def t_sum(c, cutoff=2):
    coefs = c['t']
    return sum(min(0, c + cutoff) if c < 0 else max(0, c - cutoff) for c in coefs)

# function for comparing with bump_cluster
def coef_t_prod(coefs):
    return np.median([coefs['t'][i] * coefs['coef'][i]
                        for i in range(len(coefs['coef']))])

def bump_cluster(formula, methylations, covs, coef, nsims=100000,
        value_fn=coef_sum, robust=False):
    """Model clusters by fitting model at each site and then comparing some
    metric to the same metric from models fit to simulated data.
    Uses sequential Monte-carlo to stop once we know the simulated p-value is
    high (since we are always interested in low p-values).

    Same signature as :func:`~gee_cluster`
    """
    orig = _combine_cluster(formula, methylations, covs, coef, robust=robust)
    obs_coef = value_fn(orig)

    reduced_residuals, reduced_fitted = [], []

    # get the reduced residuals and models so we can shuffle
    for i, methylation in enumerate(methylations):
        y, X = patsy.dmatrices(formula, covs, return_type='dataframe')
        idxs = [par for par in X.columns if par.startswith(coef)]
        assert len(idxs) == 1, ('too many coefficents like', coef)
        X.pop(idxs[0])
        fitr = (RLM if robust else GLS)(y, X).fit()

        reduced_residuals.append(np.array(fitr.resid))
        reduced_fitted.append(np.array(fitr.fittedvalues))

    ngt, idxs = 0, np.arange(len(methylations[0]))

    for isim in range(nsims):
        np.random.shuffle(idxs)

        fakem = np.array([rf + rr[idxs] for rf, rr in izip(reduced_fitted,
            reduced_residuals)])
        assert fakem.shape == methylations.shape

        sim = _combine_cluster(formula, fakem, covs, coef, robust=robust)
        ccut = value_fn(sim)
        ngt += abs(ccut) > abs(obs_coef)
        # sequential monte-carlo.
        if ngt > 5: break

    p = (1.0 + ngt) / (2.0 + isim) # extra 1 in denom for 0-index
    orig.pop('corr')
    orig['p'] = p
    orig['coef'], orig['t'] = orig['coef'].mean(), orig['t'].mean()
    return orig


def model_clusters(clust_iter, clin_df, formula, coef, model_fn=gee_cluster,
        pool=None,
        n_cpu=None,
        as_feature=False,
        **kwargs):
    """For each cluster in an iterable, evaluate the chosen model and
    yield a dictionary of information

    Parameters
    ----------

    clust_iter : iterable
        iterable of clusters

    clin_df : pandas.DataFrame
        Contains covariates from `formula`

    formula : str
        R (patsy) style formula. Must contain 'methylation': e.g.:
        methylation ~ age + gender + race

    coef : str
        The coefficient of interest in the model, e.g. 'age'

    model_fn : fn
        A function with signature
        fn(formula, methylation, covs, coef, kwargs)
        that returns a dictionary with at least p-value and coef

    n_cpu : int

    kwargs: dict
        arguments sent to `model_fn`
    """
    as_feature = as_feature or kwargs.pop('as_feature', False)

    for r in ts.pmap(wrapper, ((model_fn, formula, cluster, clin_df, coef,
                                kwargs, as_feature) for cluster in clust_iter), n_cpu,
                                p=pool):
        yield r


class Feature(object):
    """
    A feature object that can and likely should be used by all programs that
    call `crystal`. Takes a chromosome, a position and a list of float
    values that are the methylation measurements (should be logit transformed).

    Attributes
    ----------

    chrom: str

    position: int

    values: list

    spos: str
        string position (chr1:12354)

    rho_min : float
        minimum spearman's R to be considered correlated

    ovalues : list
        other values potentially used by modeling functions
    """

    __slots__ = "chrom position values spos rho_min ovalues".split()

    def __init__(self, chrom, pos, values, ovalues=None, rho_min=0.5):
        self.spos = "%s:%i" % (chrom, pos)
        self.chrom, self.position, self.values = chrom, pos, np.asarray(values)
        self.rho_min = rho_min
        self.ovalues = ovalues

    def distance(self, other):
        """Distance between this feature and another."""
        if self.chrom != other.chrom: return sys.maxint
        return self.position - other.position

    def is_correlated(self, other):
        """Return boolean indicating  correlation with other."""
        rho, p = ss.spearmanr(self.values, other.values)
        return rho > self.rho_min

    def __repr__(self):
        return "{cls}({spos})".format(spos=self.spos,
                                        cls=self.__class__.__name__)

    def __cmp__(self, other):
        return cmp(self.chrom, other.chrom) or cmp(self.position,
                                                   other.position)

class CountFeature(Feature):
    def __init__(self, chrom, pos, methylated, counts, rho_min=0.5):
        Feature.__init__(self, chrom, pos, methylated, counts)
        self.counts = counts
        # use ratios for correlation.
        self.values = methylated / counts.astype(float)
        self.methylated = methylated

if __name__ == "__main__":

    def feature_gen(fname):
        for i, toks in enumerate(ts.reader(fname, header=False)):
            if i == 0: continue
            chrom, pos = toks[0].split(":")
            yield Feature(chrom, int(pos), map(float, toks[1:]))

    #fmt = "{chrom}\t{start}\t{end}\t{n_probes}\t{p:5g}\t{t:.4f}\t{coef:.4f}\t{var}"
    #print ts.fmt2header(fmt)

    clust_iter = (c for c in mclust(feature_gen(sys.argv[1]),
                                    max_dist=400, max_skip=1) if len(c) > 2)


    df = pd.read_table('meth.clin.txt')
    df['id'] = df['StudyID']
    formula = "methylation ~ asthma + age + gender"

    np.random.seed(10)
    ntrue, nfalse = 10, 10

    results = []
    for fn in (bump_cluster, liptak_cluster, zscore_cluster):
        results.append(evaluate_method(clust_iter, df, formula, 'asthma', fn,
            ntrue, nfalse))

    formula = "methylation ~ asthma + age + gender"
    for fn in (gee_cluster, mixed_model_cluster):
        results.append(evaluate_method(clust_iter, df, formula, 'asthma', fn,
            ntrue, nfalse))

    results = pd.DataFrame(results)
    print pd.melt(results,
            id_vars=[c for c in results.columns if not ('false' in c or 'true' in c)],
            value_vars=[c for c in results.columns if 'false' in c or 'true' in
                c], value_name='n_lt_alpha')


