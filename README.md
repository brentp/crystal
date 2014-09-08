crystal
-------

A python module to evaluate methylation modeling strategies and find differentially methylated regions.

quick-start
-----------

All functions to model clusters take:

1. formula: a string model formula in R (or [patsy](http://patsy.readthedocs.org/en/latest/) syntax).
2. a methylation cluster from [aclust](https://github.com/brentp/aclust)
3. a pandas datafame containing the covariates of interest. The samples in this
dataframe should be in the **same order** as they are in the methylation cluster.
4. a coefficient of interest from the model.

So a single call would look like:

```Python
import pandas as pd
import aclust
import crystal

covs_df = pd.read_csv('covariates.csv')
cluster_iter = aclust.mclust(feature_gen(), max_dist=100)
crystal.zscore_cluster("methylation ~ disease + age + gender",
                        next(cluster_iter),
                        covs_df,
                        "disease")
```
Which will return a dictionary containing the `p` value, the `coef`
and the `t` statistic for 'disease' status. However it is more likely
that users will want a simple way to model every cluster:

```Python
for cluster in crystal.model_clusters(cluster_iter, covs_df,
                                      formula, "disease",
                                      crystal.zscore_cluster, n_cpu=10):
    print cluster
```

See the [introduction notebook](http://nbviewer.ipython.org/github/brentp/crystal/blob/master/notebooks/Introduction.ipynb) for more detail

Installation
------------

users unaccustomed to installing their own python packages should
download [anaconda](https://store.continuum.io/cshop/anaconda/) and
then install additional modules with `pip`.

Running `python setup.py install` should work on systems with the
scientific python stack.

How
---

*Crystal* works by generating clusters using [aclust](https://github.com/brentp/aclust) and then testing them against a chosen method. New methods can be implemented easily, but currently, there is:

1. GEE: generalized estimating equation using sample as the grouping variable
2. Mixed Effect Model: with a random intercept by sample
3. Combining P-values: test each site in a region and combine the resulting p-values using the Stouffer-Liptak or Z-score method
4. Bump(hunt)ing: A modification of bumphunting that works on local clusters and allows arbitrary metrics (not just the sum of the smoothed coefficients).

Methods using Robust regression are also available for the above.

Note that these are cleanly implemented in python thanks to the **excellent** [statsmodels package](https://github.com/statsmodels/statsmodels)


Evaluation
----------

I have evaluated a number of methods for modeling correlated data (clusters).
See the [evaluation notebook](http://nbviewer.ipython.org/github/brentp/crystal/blob/master/notebooks/crystal-methods-evaluation.ipynb) for more detail. Here
is the summary image--higher is better for the top plot and lower is better for the bottom plot:

![evaluation](https://gist.githubusercontent.com/brentp/bf7d3c3d3f23cc319ed8/raw/4fa70b33c0ad65c467dec42ef8bff77856dcc114/eval.png "Evaluation Plot")
