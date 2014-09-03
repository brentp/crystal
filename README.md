crystal
-------

A framework to evaulate methylation modeling strategies

quickstart
----------

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
cluster_iter = aclust.clust(feature_gen(), max_dist=100)
crystal.zscore_cluster("methylation ~ disease + age + gender",
                        next(cluster_iter),
                        covs_df,
                        "disease)
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

See the [introduction notebook TODO](nbviewer.org) for more detail

How
---

*Crystal* works by generating clusters using [aclust](https://github.com/brentp/aclust) and then testing them against a chosen method. New methods can be implemented easily, but currently, there is:

1. GEE: generalized estimating equation using sample as the grouping variable
2. Mixed Effect Model: with a random intercept by sample
3. Combining P-values: test each site in a region and combine the resulting p-values using the Stouffer-Liptak or Z-score method
4. Bump(hunt)ing: A modification of bumphunting that works on local clusters and allows arbitrary metrics (not just the sum of the smoothed coefficients).

Methods using Robust regression are also available for the above.

Note that these are cleanly implemented in python thanks to the **excellent** [statsmodels package](https://github.com/statsmodels/statsmodels)


