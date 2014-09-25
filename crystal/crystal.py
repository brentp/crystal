"""
Model clustered, correlated data.
"""

import sys
import toolshed as ts
from itertools import izip
import time
import re

import pandas as pd
import numpy as np
import scipy.stats as ss
import patsy

from scipy.stats import norm
from numpy.linalg import cholesky as chol, lstsq

from statsmodels.api import GEE, GLM, MixedLM, RLM, GLS, OLS, GLSAR
from statsmodels.genmod.dependence_structures import Exchangeable, Independence
from statsmodels.genmod.families import (Gaussian, Poisson,
        NegativeBinomial as NB)
from statsmodels.discrete.discrete_model import NegativeBinomial

def long_covs(covs, methylation, **kwargs):
    covs['id'] = ['id_%i' % i for i in range(len(covs))]
    cov_rep = pd.concat((covs for i in range(len(methylation))))
    #nr, nc = methylation.shape
    cov_rep['CpG'] = np.repeat(['CpG_%i' % i for i in
        range(methylation.shape[0])], methylation.shape[1])

    #cov_rep['CpG'] = np.repeat(['CpG_%i' % i for i in range(nr)], nc)
    cov_rep['methylation'] = np.concatenate(methylation)
    for k, arr in kwargs.items():
        assert arr.shape == methylation.shape
        cov_rep[k] = np.concatenate(arr)
    cov_rep.index = np.arange(len(cov_rep), dtype=int)
    return cov_rep

def corr(methylations):
    if len(methylations) == 0: return np.nan
    if len(methylations) == 1: return [[1]]
    c = np.abs(ss.spearmanr(methylations.T))
    if len(c.shape) == 1:
        assert len(methylations) == 2
        return np.array([[1, c[0]], [c[0], 1]])
    return c[0]

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
    except (ValueError, np.linalg.linalg.LinAlgError) as e:
        return dict(p=np.nan, t=np.nan, coef=np.nan,
                covar=fit.model.exog_names[idx[0]])

def one_cluster(formula, feature, covs, coef, robust=False,
        _pat=re.compile("\+\s*CpG")):
    """used when we have a "cluster" with 1 probe."""
    c = covs.copy()
    # remove the CpG in the formula
    formula = _pat.sub("", formula)
    if isinstance(feature, CountFeature):
        c['methylation'] = feature.methylated
        c['counts'] = feature.counts
        c = c[c['counts'] > 0]
        return get_ptc(GLM.from_formula(formula, data=c,
                                        exposure=c['counts'],
                                        family=Poisson()).fit(), coef)
    else:
        c['methylation'] = feature.values
        res = (RLM if robust else OLS).from_formula(formula, data=c).fit()
        return get_ptc(res, coef)


def nb_cluster(formula, cluster, covs, coef):
    """Model a cluster of correlated features with the negative binomial"""
    methylated = np.array([f.methylated for f in cluster])
    counts = np.array([f.counts for f in cluster])
    try:
        res = [NegativeBinomial.from_formula(formula, covs, offset=np.log(count))\
               .fit(disp=0) for methylation, count in izip(methylated, counts)]
    except:
        return dict(t=np.nan, coef=np.nan, covar="NA", p=np.nan, corr=np.nan)

    methylations = methylated / counts # for correlation below.

    r = get_ptc(res, coef)
    nan = np.isnan(r['p'])
    r['p'] = zscore_combine(r['p'][~nan], corr(methylations[~nan]))
    r['t'], r['coef'] = np.mean(r['t']), np.mean(r['coef'])
    return r

def gee_cluster(formula, cluster, covs, coef, cov_struct=Exchangeable(),
        family=Gaussian()):
    """An example of a `model_fn`; any function with a similar signature
    can be used.

    Parameters
    ----------

    formula : str
        R (patsy) style formula. Must contain 'methylation': e.g.:
        methylation ~ age + gender + race

    cluster : list of Features
        cluster of features from clustering or a region.
        most functions will create a methylation matrix with:
        >> meth = np.array([f.values for f in features])

    covs : pandas.DataFrame
        Contains covariates from `formula`

    coef: str
        coefficient of interest, e.g. 'age'

    cov_struct: object
        one of the covariance structures provided by statsmodels.
        Likely either Exchangeable() or Independence()

    family: object
        one of the familyies provided by statsmodels. If Guassian(),
        then methylation is assumed to be count-based (clusters of
        CountFeatures.

    Returns
    -------

    result : dict
        dict with values (keys) of at least p-value ('p'), coefficient
        estimate ('coef') and any other information desired.
    """
    if isinstance(family, Gaussian):
        cov_rep = long_covs(covs, np.array([f.values for f in cluster]))
        res = GEE.from_formula(formula, groups=cov_rep['id'], data=cov_rep,
            cov_struct=cov_struct, family=family).fit()
    elif isinstance(family, (NB, Poisson)):
        cov_rep = long_covs(covs, np.array([f.methylated for f in cluster]),
                counts = np.array([f.counts for f in cluster]))
        res = GEE.from_formula(formula, groups=cov_rep['id'], data=cov_rep,
                               cov_struct=cov_struct, family=family,
                               offset=np.log(cov_rep['counts'])).fit()
    else:
        raise Exception("Only guassian and poisson are supported")

    return get_ptc(res, coef)

# see: https://gist.github.com/josef-pkt/89585d0b084739a4ed1c
def ols_cluster_robust(formula, cluster, covs, coef):
    """Model clusters with cluster-robust OLS, same signature as
    :func:`~gee_cluster`"""
    cov_rep = long_covs(covs, np.array([f.values for f in cluster]))
    res = OLS.from_formula(formula, data=cov_rep).fit(cov_type='cluster',
            cov_kwds=dict(groups=cov_rep['id']))
    return get_ptc(res, coef)

def glsar_cluster(formula, cluster, covs, coef, rho=6):
    cov_rep = long_covs(covs, np.array([f.values for f in cluster]))
    # group by id, then sort by CpG so that AR is to adjacent CpGs
    cov_rep.sort(['id', 'CpG'], inplace=True)
    res = GLSAR.from_formula(formula, data=cov_rep, rho=rho).iterative_fit(maxiter=5)
    return get_ptc(res, coef)

def gls_cluster(formula, cluster, covs, coef):
    methylation = np.array([f.values for f in cluster])
    # TODO: currently this is very, very slow
    cov_rep = long_covs(covs, methylation)
    c = np.array(pd.DataFrame(methylation).corr())
    n = methylation.shape[0]
    sigma = np.repeat(np.repeat(c, n, axis=0), n, axis=1)

    sigma[np.tril_indices_from(sigma)] = 0
    sigma[np.diag_indices_from(sigma)] = 1

    res = GLS.from_formula(formula, data=cov_rep, sigma=sigma).fit()
    return get_ptc(res, coef)

def mixed_model_cluster(formula, cluster, covs, coef):
    """Model clusters with a mixed-model, same signature as
    :func:`~gee_cluster`"""
    cov_rep = long_covs(covs, np.array([f.values for f in cluster]))
    # TODO: remove this once newer version of statsmodels is out.
    # speeds convergence by using fixed estimates from OLS
    params = OLS.from_formula(formula, data=cov_rep).fit().params

    res = MixedLM.from_formula(formula, groups='id',
            data=cov_rep).fit(start_params=dict(fe=params), reml=False,
                    method='bfgs')

    return get_ptc(res, coef)

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

def _combine_cluster(formula, methylations, covs, coef, robust=False,
        _corr=True):
    """function called by z-score and liptak to get pvalues"""
    res = [(RLM if robust else OLS).from_formula(formula, covs).fit()
            for methylation in methylations]
    idx = [i for i, par in enumerate(res[0].model.exog_names)
                   if par.startswith(coef)][0]
    pvals = np.array([r.pvalues[idx] for r in res], dtype=np.float64)
    pvals[pvals == 1] = 1.0 - 9e-16
    res = dict(t=np.array([r.tvalues[idx] for r in res]),
                coef=np.array([r.params[idx] for r in res]),
                covar=res[0].model.exog_names[idx],
                p=pvals[~np.isnan(pvals)])
    # save some time for bumphunting where we don't need
    # the correlation.
    if _corr:
        res['corr'] = corr(methylations[~np.isnan(pvals)])
    return res

def bayes_cluster():
    pass

def liptak_cluster(formula, cluster, covs, coef, robust=False):
    """Model clusters by fitting model at each site and then
    combining using :func:`~stouffer_liptak`. same signature as
    :func:`~gee_cluster`"""
    methylations = np.array([f.values for f in cluster])
    r = _combine_cluster(formula, methylations, covs, coef, robust=robust)
    r['p'] = stouffer_liptak_combine(r['p'], r['corr'])
    r['t'], r['coef'] = r['t'].mean(), r['coef'].mean()
    r.pop('corr')
    return r

def zscore_cluster(formula, cluster, covs, coef, robust=False):
    """Model clusters by fitting model at each site and then
    combining using the z-score method. Same signature as
    :func:`~gee_cluster`"""
    methylations = np.array([f.values for f in cluster])
    r = _combine_cluster(formula, methylations, covs, coef, robust=robust)
    r['p'] = zscore_combine(r['p'], r.pop('corr'))
    r['t'], r['coef'] = r['t'].mean(), r['coef'].mean()
    return r

# function for comparing with bump_cluster
# takes the return value of _combine_cluster and returns a single numeric value
def coef_sum(c, cutoff=0.01):
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


def _cluster_coefs(formula, y, covs, coef):
    """Fast coefficient calculation for bumping."""
    X = patsy.dmatrix(formula.split("~")[-1], covs, return_type="dataframe")
    idx = [i for i, c in enumerate(X.columns) if c.startswith(coef)]
    assert len(idx) == 1
    x, resids, rank, s = lstsq(X, y.T)
    coefs = x[idx[0]]
    return coefs

def bump_cluster(formula, cluster, covs, coef, nsims=20000,
        value_fn=coef_sum, robust=False):
    """Model clusters by fitting model at each site and then comparing some
    metric to the same metric from models fit to simulated data.
    Uses sequential Monte-carlo to stop once we know the simulated p-value is
    high (since we are always interested in low p-values).

    Same signature as :func:`~gee_cluster`
    """
    methylations = np.array([f.values for f in cluster])
    orig = _combine_cluster(formula, methylations, covs, coef, robust=robust)
    obs_coef = value_fn(orig)

    reduced_residuals, reduced_fitted = [], []

    # get the reduced residuals and models so we can shuffle
    for i, methylation in enumerate(methylations):
        y, X = patsy.dmatrices(formula, covs, return_type='dataframe')
        idxs = [par for par in X.columns if par.startswith(coef)]
        assert len(idxs) == 1, ('too many coefficents like', coef)
        X.pop(idxs[0])
        fitr = (RLM if robust else OLS)(y, X).fit()

        reduced_residuals.append(np.array(fitr.resid))
        reduced_fitted.append(np.array(fitr.fittedvalues))

    ngt, idxs = 0, np.arange(len(methylations[0]))

    for isim in range(nsims):
        np.random.shuffle(idxs)

        fakem = np.array([rf + rr[idxs] for rf, rr in izip(reduced_fitted,
            reduced_residuals)])
        assert fakem.shape == methylations.shape

        coefs = _cluster_coefs(formula, fakem, covs, coef)
        ccut = value_fn(dict(coef=coefs))

        #sim = _combine_cluster(formula, fakem, covs, coef, robust=robust,
        #        _corr=False)
        #ccut = value_fn(sim)

        ngt += abs(ccut) > abs(obs_coef)
        # sequential monte-carlo.
        if ngt > 6: break

    orig.pop('corr')
    # extra 1 in denom for 0-index
    orig['n_sim'] = isim + 1
    # extra 1 so we can't get 0 p-value
    orig['p'] = (1.0 + ngt) / (1.0 + orig['n_sim'])
    orig['coef'], orig['t'] = orig['coef'].mean(), orig['t'].mean()
    return orig

def wrapper(model_fn, formula, cluster, clin_df, coef, kwargs=None):
    """wrap the user-defined functions to return everything we expect and
    to call just OLS when there is a single probe."""
    if kwargs is None: kwargs = {}
    t = time.time()
    if len(cluster) > 1:
        r = model_fn(formula,
                cluster,
                clin_df,
                coef, **kwargs)
    else:
        r = one_cluster(formula, cluster[0], clin_df, coef)
    r['time'] = time.time() - t
    r['chrom'] = cluster[0].chrom
    r['start'] = cluster[0].position - 1
    r['end'] = cluster[-1].position
    r['n_sites'] = len(cluster)
    r['sites'] = [c.spos for c in cluster]
    r['var'] = coef
    r['cluster'] = cluster
    return r

def model_clusters(clust_iter, clin_df, formula, coef, model_fn=gee_cluster,
        pool=None,
        n_cpu=None,
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
    for r in ts.pmap(wrapper, ((model_fn, formula, cluster, clin_df, coef,
                                kwargs) for cluster in clust_iter), n_cpu,
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

    """Feature Class that supports count data."""

    def __init__(self, chrom, pos, methylated, counts, rho_min=0.5):
        Feature.__init__(self, chrom, pos, methylated, counts)
        self.counts = counts
        # use ratios for correlation.
        self.values = methylated / counts.astype(float)
        self.methylated = methylated
