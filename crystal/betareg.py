import numpy as np
from scipy.special import gammaln as lgamma
from statsmodels.base.model import GenericLikelihoodModel
from statsmodels.api import GLM
from statsmodels.genmod.families import Binomial
import statsmodels.api as sm
import pandas as pd

#see http://cran.r-project.org/web/packages/betareg/vignettes/betareg-ext.pdf
# nice reference:
# http://psychology3.anu.edu.au/people/smithson/details/Pubs/Smithson_Verkuilen06.pdf

class Logit(sm.families.links.Logit):
    def inverse(self, z):
        return 1 / (1. + np.exp(-z))

class BetaReg(GenericLikelihoodModel):
    def __init__(self, endog, exog, Z=None, link=Logit(),
            link_phi=sm.families.links.Log(), **kwds):
        assert np.all((0 < endog) & (endog < 1))

        super(BetaReg, self).__init__(endog, exog, **kwds)
        self.link = link
        self.link_phi = link_phi
        # how to set default Z?
        if Z is None:
            self.Z = np.ones((self.endog.shape[0], 1), dtype='f')
        else:
            self.Z = np.asarray(Z)
            assert len(self.Z) == len(self.endog)
    def nloglikeobs(self, params):
        return -self._ll_br(self.endog, self.exog, self.Z, params)

    def fit(self, start_params=None, maxiter=1000000, maxfun=50000,
            disp=False,
            method='bfgs', **kwds):
        if start_params is None:
            start_params = GLM(self.endog, self.exog, family=Binomial()
                              ).fit(disp=False).params
            start_params = np.append(start_params, [0.5] * self.Z.shape[1])
            #start_params = np.append(np.zeros(self.exog.shape[1]), 0.5)
        #self.exog[0] = np.mean(self.endog)

        return super(BetaReg, self).fit(start_params=start_params,
                                             maxiter=maxiter,
                                             maxfun=maxfun,
                                             method=method,
                                             disp=disp,
                                             **kwds)
    def _ll_br(self, y, X, Z, params):
        nz = self.Z.shape[1]

        Xparams = params[:-nz]
        Zparams = params[-nz:]

        mu = self.link.inverse(np.dot(X, Xparams))
        phi = self.link_phi.inverse(np.dot(Z, Zparams))

        if np.any(phi <= np.finfo(float).eps): return np.array(-inf)

        ll = lgamma(phi) - lgamma(mu * phi) - lgamma((1 - mu) * phi) \
                + (mu * phi - 1) * np.log(y) + (((1 - mu) * phi) - 1) * np.log(1 - y)
        return ll

if __name__ == "__main__":

    import pandas as pd
    import patsy
    dat = pd.read_table('gasoline.txt')
    Z = patsy.dmatrix('~ temp', dat, return_type='dataframe')
    # using other precison params with
    m = BetaReg.from_formula('iyield ~ C(batch, Treatment(10)) + temp', dat,
            Z=Z, link_phi=sm.families.links.identity())
    print m.fit().summary()

    fex = pd.read_csv('foodexpenditure.csv')
    m = BetaReg.from_formula(' I(food/income) ~ income + persons', fex)
    print m.fit().summary()
    #print GLM.from_formula('iyield ~ C(batch) + temp', dat, family=Binomial()).fit().summary()

    dev = pd.read_csv('methylation-test.csv')
    dev['methylation'] = sm.families.links.Logit().inverse(dev['methylation'])
    m = BetaReg.from_formula('methylation ~ age + gender', dev,
            link_phi=sm.families.links.identity())
    print m.fit().summary()
