.. crystal documentation master file, created by
   sphinx-quickstart on Mon Sep  8 11:12:44 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

crystal documentation
=====================

`crystal` is a python module to assign significance to clusters of
correlated data. The most likely use-case is for DNA methylation data, but
it can be use for any correlated data.

New users should check-out the introductory IPython notebook `here <http://nbviewer.ipython.org/github/brentp/crystal/blob/master/notebooks/Introduction.ipynb>`_.

API Docs
========

.. currentmodule:: sphinx

.. toctree::
   :maxdepth: 3

   crystal
   plot
   regions

Introduction
============

crystal models clusters of correlated data.

quick-start
-----------

All functions to model clusters take:

1. formula: a string model formula in R (or `patsy <http://patsy.readthedocs.org/en/latest/>`_ syntax). 
2. a methylation cluster from `aclust <https://github.com/brentp/aclust>`_
3. a pandas datafame containing the covariates of interest. The samples in this
   dataframe should be in the **same order** as they are in the methylation cluster.
4. a coefficient of interest from the model.

So a single call would look like:

.. code-block:: python

    import pandas as pd
    import aclust
    import crystal

    covs_df = pd.read_csv('covariates.csv')
    cluster_iter = aclust.mclust(feature_gen(), max_dist=100)
    crystal.zscore_cluster("methylation ~ disease + age + gender",
                            next(cluster_iter),
                            covs_df,
                            "disease")

Which will return a dictionary containing the `p` value, the `coef`
and the `t` statistic for 'disease' status. However it is more likely
that users will want a simple way to model every cluster:

.. code-block:: python 

    for cluster in crystal.model_clusters(cluster_iter,
                                          covs_df,
                                          formula,
                                          "disease",
                                          crystal.zscore_cluster,
                                          n_cpu=10):
        print cluster

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

